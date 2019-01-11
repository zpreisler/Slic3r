#include "SLABasePool.hpp"
#include "SLABoilerPlate.hpp"

#include "boost/log/trivial.hpp"
#include "SLABoostAdapter.hpp"
#include "ClipperUtils.hpp"

// We need the EdgeCache class from libnest2d's nfpplacer
#include <libnest2d/backends/clipper/geometries.hpp>
#include <libnest2d/placers/nfpplacer.hpp>

//#include "SVG.hpp"
//#include "benchmark.h"

namespace Slic3r { namespace sla {

/// Convert the triangulation output to an intermediate mesh.
Contour3D convert(const Polygons& triangles, coord_t z, bool dir) {

    Pointf3s points;
    points.reserve(3*triangles.size());
    Indices indices;
    indices.reserve(points.size());

    for(auto& tr : triangles) {
        auto c = coord_t(points.size()), b = c++, a = c++;
        if(dir) indices.emplace_back(a, b, c);
        else indices.emplace_back(c, b, a);
        for(auto& p : tr.points) {
            points.emplace_back(unscale(x(p), y(p), z));
        }
    }

    return {points, indices};
}

Contour3D walls(const ExPolygon& floor_plate, const ExPolygon& ceiling,
                double floor_z_mm, double ceiling_z_mm,
                ThrowOnCancel thr)
{
    using std::transform; using std::back_inserter;

    ExPolygon poly;
    poly.contour.points = floor_plate.contour.points;
    poly.holes.emplace_back(ceiling.contour);
    auto& h = poly.holes.front();
    std::reverse(h.points.begin(), h.points.end());
    Polygons tri = triangulate(poly);

    Contour3D ret;
    ret.points.reserve(tri.size() * 3);

    double fz = floor_z_mm;
    double cz = ceiling_z_mm;
    auto& rp = ret.points;
    auto& rpi = ret.indices;
    ret.indices.reserve(tri.size() * 3);

    coord_t idx = 0;

    auto hlines = h.lines();
    auto is_upper = [&hlines](const Point& p) {
        return std::any_of(hlines.begin(), hlines.end(),
                               [&p](const Line& l) {
            return l.distance_to(p) < mm(1e-6);
        });
    };

    std::for_each(tri.begin(), tri.end(),
                  [&rp, &rpi, thr, &idx, is_upper, fz, cz](const Polygon& pp)
    {
        thr(); // may throw if cancellation was requested

        for(auto& p : pp.points)
            if(is_upper(p))
                rp.emplace_back(unscale(x(p), y(p), mm(cz)));
            else rp.emplace_back(unscale(x(p), y(p), mm(fz)));

        coord_t a = idx++, b = idx++, c = idx++;
        if(fz > cz) rpi.emplace_back(c, b, a);
        else rpi.emplace_back(a, b, c);
    });

    return ret;
}

/// Offsetting with clipper and smoothing the edges into a curvature.
void offset(ExPolygon& sh, coord_t distance, bool edgerounding = true) {
    using ClipperLib::ClipperOffset;
    using ClipperLib::jtRound;
    using ClipperLib::jtMiter;
    using ClipperLib::etClosedPolygon;
    using ClipperLib::Paths;
    using ClipperLib::Path;

    auto&& ctour = Slic3rMultiPoint_to_ClipperPath(sh.contour);
    auto&& holes = Slic3rMultiPoints_to_ClipperPaths(sh.holes);

    // If the input is not at least a triangle, we can not do this algorithm
    if(ctour.size() < 3 ||
       std::any_of(holes.begin(), holes.end(),
                   [](const Path& p) { return p.size() < 3; })
            ) {
        BOOST_LOG_TRIVIAL(error) << "Invalid geometry for offsetting!";
        return;
    }

    auto jointype = edgerounding? jtRound : jtMiter;

    ClipperOffset offs;
    offs.ArcTolerance = 0.01*mm(1);
    Paths result;
    offs.AddPath(ctour, jointype, etClosedPolygon);
    offs.AddPaths(holes, jointype, etClosedPolygon);
    offs.Execute(result, static_cast<double>(distance));

    // Offsetting reverts the orientation and also removes the last vertex
    // so boost will not have a closed polygon.

    bool found_the_contour = false;
    sh.holes.clear();
    for(auto& r : result) {
        if(ClipperLib::Orientation(r)) {
            // We don't like if the offsetting generates more than one contour
            // but throwing would be an overkill. Instead, we should warn the
            // caller about the inability to create correct geometries
            if(!found_the_contour) {
                auto rr = ClipperPath_to_Slic3rPolygon(r);
                sh.contour.points.swap(rr.points);
                found_the_contour = true;
            } else {
                BOOST_LOG_TRIVIAL(warning)
                        << "Warning: offsetting result is invalid!";
            }
        } else {
            // TODO If there are multiple contours we can't be sure which hole
            // belongs to the first contour. (But in this case the situation is
            // bad enough to let it go...)
            sh.holes.emplace_back(ClipperPath_to_Slic3rPolygon(r));
        }
    }
}

/// Unification of polygons (with clipper) preserving holes as well.
ExPolygons unify(const ExPolygons& shapes) {
    using ClipperLib::ptSubject;

    ExPolygons retv;

    bool closed = true;
    bool valid = true;

    ClipperLib::Clipper clipper;

    for(auto& path : shapes) {
        auto clipperpath = Slic3rMultiPoint_to_ClipperPath(path.contour);

        if(!clipperpath.empty())
            valid &= clipper.AddPath(clipperpath, ptSubject, closed);

        auto clipperholes = Slic3rMultiPoints_to_ClipperPaths(path.holes);

        for(auto& hole : clipperholes) {
            if(!hole.empty())
                valid &= clipper.AddPath(hole, ptSubject, closed);
        }
    }

    if(!valid) BOOST_LOG_TRIVIAL(warning) << "Unification of invalid shapes!";

    ClipperLib::PolyTree result;
    clipper.Execute(ClipperLib::ctUnion, result, ClipperLib::pftNonZero);

    retv.reserve(static_cast<size_t>(result.Total()));

    // Now we will recursively traverse the polygon tree and serialize it
    // into an ExPolygon with holes. The polygon tree has the clipper-ish
    // PolyTree structure which alternates its nodes as contours and holes

    // A "declaration" of function for traversing leafs which are holes
    std::function<void(ClipperLib::PolyNode*, ExPolygon&)> processHole;

    // Process polygon which calls processHoles which than calls processPoly
    // again until no leafs are left.
    auto processPoly = [&retv, &processHole](ClipperLib::PolyNode *pptr) {
        ExPolygon poly;
        poly.contour.points = ClipperPath_to_Slic3rPolygon(pptr->Contour);
        for(auto h : pptr->Childs) { processHole(h, poly); }
        retv.push_back(poly);
    };

    // Body of the processHole function
    processHole = [&processPoly](ClipperLib::PolyNode *pptr, ExPolygon& poly)
    {
        poly.holes.emplace_back();
        poly.holes.back().points = ClipperPath_to_Slic3rPolygon(pptr->Contour);
        for(auto c : pptr->Childs) processPoly(c);
    };

    // Wrapper for traversing.
    auto traverse = [&processPoly] (ClipperLib::PolyNode *node)
    {
        for(auto ch : node->Childs) {
            processPoly(ch);
        }
    };

    // Here is the actual traverse
    traverse(&result);

    return retv;
}

/// Only a debug function to generate top and bottom plates from a 2D shape.
/// It is not used in the algorithm directly.
inline Contour3D roofs(const ExPolygon& poly, coord_t z_distance) {
    Polygons triangles = triangulate(poly);

    auto lower = convert(triangles, 0, false);
    auto upper = convert(triangles, z_distance, true);
    lower.merge(upper);
    return lower;
}

Contour3D round_edges(const ExPolygon& base_plate,
                      double radius_mm,
                      double degrees,
                      double ceilheight_mm,
                      bool dir,
                      ThrowOnCancel throw_on_cancel,
                      ExPolygon& last_offset, double& last_height)
{
    auto ob = base_plate;
    auto ob_prev = ob;
    double wh = ceilheight_mm, wh_prev = wh;
    Contour3D curvedwalls;

    int steps = 30;
    double stepx = radius_mm / steps;
    coord_t s = dir? 1 : -1;
    degrees = std::fmod(degrees, 180);

    // we use sin for x distance because we interpret the angle starting from
    // PI/2
    int tos = degrees < 90?
               int(radius_mm*std::cos(degrees * PI / 180 - PI/2) / stepx) : steps;

    for(int i = 1; i <= tos; ++i) {
        throw_on_cancel();

        ob = base_plate;

        double r2 = radius_mm * radius_mm;
        double xx = i*stepx;
        double x2 = xx*xx;
        double stepy = std::sqrt(r2 - x2);

        offset(ob, s*mm(xx));
        wh = ceilheight_mm - radius_mm + stepy;

        Contour3D pwalls;
        pwalls = walls(ob, ob_prev, wh, wh_prev, throw_on_cancel);

        curvedwalls.merge(pwalls);
        ob_prev = ob;
        wh_prev = wh;
    }

    if(degrees > 90) {
        double tox = radius_mm - radius_mm*std::cos(degrees * PI / 180 - PI/2);
        int tos = int(tox / stepx);

        for(int i = 1; i <= tos; ++i) {
            throw_on_cancel();
            ob = base_plate;

            double r2 = radius_mm * radius_mm;
            double xx = radius_mm - i*stepx;
            double x2 = xx*xx;
            double stepy = std::sqrt(r2 - x2);
            offset(ob, s*mm(xx));
            wh = ceilheight_mm - radius_mm - stepy;

            Contour3D pwalls;
            pwalls = walls(ob_prev, ob, wh_prev, wh, throw_on_cancel);

            curvedwalls.merge(pwalls);
            ob_prev = ob;
            wh_prev = wh;
        }
    }

    last_offset = std::move(ob);
    last_height = wh;

    return curvedwalls;
}

/// Generating the concave part of the 3D pool with the bottom plate and the
/// side walls.
Contour3D inner_bed(const ExPolygon& poly, double depth_mm,
                           double begin_h_mm = 0) {

    Polygons triangles = triangulate(poly);

    coord_t depth = mm(depth_mm);
    coord_t begin_h = mm(begin_h_mm);

    auto bottom = convert(triangles, -depth + begin_h, false);
    auto lines = poly.lines();

    // Generate outer walls
    auto fp = [](const Point& p, Point::coord_type z) {
        return unscale(x(p), y(p), z);
    };

    for(auto& l : lines) {
        auto s = coord_t(bottom.points.size());

        bottom.points.emplace_back(fp(l.a, -depth + begin_h));
        bottom.points.emplace_back(fp(l.b, -depth + begin_h));
        bottom.points.emplace_back(fp(l.a, begin_h));
        bottom.points.emplace_back(fp(l.b, begin_h));

        bottom.indices.emplace_back(s + 3, s + 1, s);
        bottom.indices.emplace_back(s + 2, s + 3, s);
    }

    return bottom;
}

inline Point centroid(Points& pp) {
    Point c;
    switch(pp.size()) {
    case 0: break;
    case 1: c = pp.front(); break;
    case 2: c = (pp[0] + pp[1]) / 2; break;
    default: {
        auto MAX = std::numeric_limits<Point::coord_type>::max();
        auto MIN = std::numeric_limits<Point::coord_type>::min();
        Point min = {MAX, MAX}, max = {MIN, MIN};

        for(auto& p : pp) {
            if(p(0) < min(0)) min(0) = p(0);
            if(p(1) < min(1)) min(1) = p(1);
            if(p(0) > max(0)) max(0) = p(0);
            if(p(1) > max(1)) max(1) = p(1);
        }
        c(0) = min(0) + (max(0) - min(0)) / 2;
        c(1) = min(1) + (max(1) - min(1)) / 2;

        // TODO: fails for non convex cluster
//        c = std::accumulate(pp.begin(), pp.end(), Point{0, 0});
//        x(c) /= coord_t(pp.size()); y(c) /= coord_t(pp.size());
        break;
    }
    }

    return c;
}

inline Point centroid(const ExPolygon& poly) {
    return poly.contour.centroid();
}

/// A fake concave hull that is constructed by connecting separate shapes
/// with explicit bridges. Bridges are generated from each shape's centroid
/// to the center of the "scene" which is the centroid calculated from the shape
/// centroids (a star is created...)
ExPolygons concave_hull(const ExPolygons& polys, double max_dist_mm = 50,
                        ThrowOnCancel throw_on_cancel = [](){})
{
    namespace bgi = boost::geometry::index;
    using SpatElement = std::pair<BoundingBox, unsigned>;
    using SpatIndex = bgi::rtree< SpatElement, bgi::rstar<16, 4> >;

    if(polys.empty()) return ExPolygons();

    ExPolygons punion = unify(polys);   // could be redundant

    if(punion.size() == 1) return punion;

    // We get the centroids of all the islands in the 2D slice
    Points centroids; centroids.reserve(punion.size());
    std::transform(punion.begin(), punion.end(), std::back_inserter(centroids),
                   [](const ExPolygon& poly) { return centroid(poly); });


    SpatIndex boxindex; unsigned idx = 0;
    std::for_each(punion.begin(), punion.end(),
                  [&boxindex, &idx](const ExPolygon& expo) {
        BoundingBox bb(expo);
        boxindex.insert(std::make_pair(bb, idx++));
    });


    // Centroid of the centroids of islands. This is where the additional
    // connector sticks are routed.
    Point cc = centroid(centroids);

    punion.reserve(punion.size() + centroids.size());

    idx = 0;
    std::transform(centroids.begin(), centroids.end(),
                   std::back_inserter(punion),
                   [&punion, &boxindex, cc, max_dist_mm, &idx, throw_on_cancel]
                   (const Point& c)
    {
        throw_on_cancel();
        double dx = x(c) - x(cc), dy = y(c) - y(cc);
        double l = std::sqrt(dx * dx + dy * dy);
        double nx = dx / l, ny = dy / l;
        double max_dist = mm(max_dist_mm);

        ExPolygon& expo = punion[idx++];
        BoundingBox querybb(expo);

        querybb.offset(max_dist);
        std::vector<SpatElement> result;
        boxindex.query(bgi::intersects(querybb), std::back_inserter(result));
        if(result.size() <= 1) return ExPolygon();

        ExPolygon r;
        auto& ctour = r.contour.points;

        ctour.reserve(3);
        ctour.emplace_back(cc);

        Point d(coord_t(mm(1)*nx), coord_t(mm(1)*ny));
        ctour.emplace_back(c + Point( -y(d),  x(d) ));
        ctour.emplace_back(c + Point(  y(d), -x(d) ));
        offset(r, mm(1));

        return r;
    });

    punion = unify(punion);

    return punion;
}

void base_plate(const TriangleMesh &mesh, ExPolygons &output, float h,
                float layerh, ThrowOnCancel thrfn)
{
    TriangleMesh m = mesh;
    TriangleMeshSlicer slicer(&m);

    auto bb = mesh.bounding_box();
    float gnd = float(bb.min(Z));
    std::vector<float> heights = {float(bb.min(Z))};
    for(float hi = gnd + layerh; hi <= gnd + h; hi += layerh)
        heights.emplace_back(hi);

    std::vector<ExPolygons> out; out.reserve(size_t(std::ceil(h/layerh)));
    slicer.slice(heights, &out, thrfn);

    size_t count = 0; for(auto& o : out) count += o.size();
    ExPolygons tmp; tmp.reserve(count);
    for(auto& o : out) for(auto& e : o) tmp.emplace_back(std::move(e));

    ExPolygons utmp = unify(tmp);
    for(auto& o : utmp) {
        auto&& smp = o.simplify(0.1/SCALING_FACTOR);
        output.insert(output.end(), smp.begin(), smp.end());
    }
}

void offset_with_breakstick_holes(ExPolygon& expoly,
                                  double padding,
                                  double stride,
                                  double stick_width) {

    // We included libnest2d clipper backend so PolygonImpl is compatible with
    // clipper Path and PointImpl is basically ClipperLib::IntPoint
    using libnest2d::PolygonImpl;
    using libnest2d::PointImpl;
    using ClipperLib::ctDifference;
    using ClipperLib::ptSubject;
    using ClipperLib::ptClip;

    // Summon our main tool from libnest2d
    using EdgeCache = libnest2d::placers::EdgeCache<PolygonImpl>;

    // We do the basic offsetting first
    const bool dont_round_edges = false;
    offset(expoly, coord_t(padding / SCALING_FACTOR), dont_round_edges);

    // Ok, we need the edge-cache...
    // go around the polygon and add additional vertices. Four points for
    // each breakstick. Care must be taken for the right orientation of the
    // added points.

    // We included libnest2d clipper backend so PolygonImpl is compatible with
    // clipper Path
    PolygonImpl poly;
    poly.Contour = Slic3rMultiPoint_to_ClipperPath(expoly.contour);
    // Holes will be ignored
    // std::reverse(poly.Contour.begin(), poly.Contour.end());
    // poly.Holes = Slic3rMultiPoints_to_ClipperPaths(expoly.holes);
    // for(auto& h : poly.Holes) std::reverse(h.begin(), h.end());

    EdgeCache ecache(poly);

    // still in clipper coordinates
    double circ = ecache.circumference() * SCALING_FACTOR;
    auto count = unsigned(circ / stride);
    double q = 1.0 / circ;
    double dwidth = stick_width * q ;
    auto swidth   = coord_t(stick_width / SCALING_FACTOR);
    auto spadding = coord_t(padding / SCALING_FACTOR);
    bool polygon_is_closed = true;

    ClipperLib::Clipper clipper;
    clipper.AddPath(poly.Contour, ptSubject, polygon_is_closed);

    for(unsigned i = 0; i < count; ++i) {
        double loc = i * stride * q;

        PointImpl p1 = ecache.coords(loc - dwidth);
        PointImpl pq = ecache.coords(loc + dwidth); // just for the normal

        Vec2d p1d(p1.X, p1.Y); p1d *= SCALING_FACTOR;
        Vec2d pqd(pq.X, pq.Y); pqd *= SCALING_FACTOR;
        auto d = (pqd - p1d).normalized();          // direction vector
        Vec2d n(-d(Y), d(X));                       // normal

        // Now we have the two points

        PolygonImpl stick;
        stick.Contour.emplace_back(p1); // emplace the starting point
        PointImpl ds(coord_t(d(X)*swidth), coord_t(d(Y)*swidth));
        PointImpl ns(coord_t(n(X)*spadding), coord_t(n(Y)*spadding));

        auto p2 = p1 + ds;
        auto p3 = p2 + ns;
        auto p4 = p1 + ns;

        clipper.AddPath({p1, p2, p3, p4}, ptClip, polygon_is_closed);
    }

    ClipperLib::Paths sol;
    clipper.Execute(ctDifference, sol);

//    SVG svg("bridgestick_plate.svg");
//    svg.draw(sol, 1);
//    svg.Close();
    if(!sol.empty()) expoly.contour = ClipperPath_to_Slic3rPolygon(sol.front());
}

void create_base_pool(const ExPolygons &ground_layer,
                      const Polygon &object_self_pad,
                      TriangleMesh& out,
                      const PoolConfig& cfg)
{

    double mergedist = 2*(1.8*cfg.min_wall_thickness_mm + 4*cfg.edge_radius_mm)+
                       cfg.max_merge_distance_mm;

    // Here we get the base polygon from which the pad has to be generated.
    // We create an artificial concave hull from this polygon and that will
    // serve as the bottom plate of the pad. We will offset this concave hull
    // and then offset back the result with clipper with rounding edges ON. This
    // trick will create a nice rounded pad shape.
    auto concavehs = concave_hull(ground_layer, mergedist, cfg.throw_on_cancel);

    const double thickness      = cfg.min_wall_thickness_mm;
    const double wingheight     = cfg.min_wall_height_mm;
    const double fullheight     = wingheight + thickness;
    const double tilt = PI/4;
    const double wingdist       = wingheight / std::tan(tilt);

    // scaled values
    const coord_t s_thickness   = mm(thickness);
    const coord_t s_eradius     = mm(cfg.edge_radius_mm);
    const coord_t s_safety_dist = 2*s_eradius + coord_t(0.8*s_thickness);
    // const coord_t wheight    = mm(cfg.min_wall_height_mm);
    coord_t s_wingdist          = mm(wingdist);
    coord_t s_wingheight        = mm(wingheight);

    auto& thrcl = cfg.throw_on_cancel;

    for(ExPolygon& concaveh : concavehs) {
        if(concaveh.contour.points.empty()) return;

        // Get rif of any holes in the concave hull output.
        concaveh.holes.clear();

        // Here lies the trick that does the smooting only with clipper offset
        // calls. The offset is configured to round edges. Inner edges will
        // be rounded because we offset twice: ones to get the outer (top) plate
        // and again to get the inner (bottom) plate
        auto outer_base = concaveh;
        outer_base.holes.clear();
        offset(outer_base, s_safety_dist + s_wingdist + s_thickness);
        auto inner_base = outer_base;
        offset(inner_base, -(s_thickness + s_wingdist));

        // Punching a hole in the top plate for the cavity
        ExPolygon top_poly;
        ExPolygon middle_base;
        top_poly.contour = outer_base.contour;

        if(wingheight > 0) {
            middle_base = outer_base;
            offset(middle_base, -s_thickness);
            top_poly.holes.emplace_back(middle_base.contour);
            auto& tph = top_poly.holes.back().points;
            std::reverse(tph.begin(), tph.end());
        }

        Contour3D pool;

        ExPolygon ob = outer_base; double wh = 0;

        // now we will calculate the angle or portion of the circle from
        // pi/2 that will connect perfectly with the bottom plate.
        // this is a tangent point calculation problem and the equation can
        // be found for example here:
        // http://www.ambrsoft.com/TrigoCalc/Circles2/CirclePoint/CirclePointDistance.htm
        // the y coordinate would be:
        // y = cy + (r^2*py - r*px*sqrt(px^2 + py^2 - r^2) / (px^2 + py^2)
        // where px and py are the coordinates of the point outside the circle
        // cx and cy are the circle center, r is the radius
        // We place the circle center to (0, 0) in the calculation the make
        // things easier.
        // to get the angle we use arcsin function and subtract 90 degrees then
        // flip the sign to get the right input to the round_edge function.
        double r = cfg.edge_radius_mm;
        double cy = 0;
        double cx = 0;
        double px = thickness + wingdist;
        double py = r - fullheight;

        double pxcx = px - cx;
        double pycy = py - cy;
        double b_2 = pxcx*pxcx + pycy*pycy;
        double r_2 = r*r;
        double D = std::sqrt(b_2 - r_2);
        double vy = (r_2*pycy - r*pxcx*D) / b_2;
        double phi = -(std::asin(vy/r) * 180 / PI - 90);


        // Generate the smoothed edge geometry
        auto walledges = round_edges(ob,
                                     r,
                                     phi,
                                     0,    // z position of the input plane
                                     true,
                                     thrcl,
                                     ob, wh);
        pool.merge(walledges);

        // Now that we have the rounded edge connencting the top plate with
        // the outer side walls, we can generate and merge the sidewall geometry
        auto pwalls = walls(ob, inner_base, wh, -fullheight, thrcl);
        pool.merge(pwalls);

        if(wingheight > 0) {
            // Generate the smoothed edge geometry
            auto cavityedges = round_edges(middle_base,
                                           r,
                                           phi - 90, // from tangent lines
                                           0,
                                           false,
                                           thrcl,
                                           ob, wh);
            pool.merge(cavityedges);

            // Next is the cavity walls connecting to the top plate's
            // artificially created hole.
            auto cavitywalls = walls(inner_base, ob, -wingheight, wh, thrcl);
            pool.merge(cavitywalls);
        }

        // Now we need to triangulate the top and bottom plates as well as the
        // cavity bottom plate which is the same as the bottom plate but it is
        // eleveted by the thickness.
        Polygons top_triangles, bottom_triangles;

        if(!object_self_pad.empty()) {
            // cutting the object shape into the pad with small breakable sticks

            ExPolygon object_base;
            object_base.contour = object_self_pad;
            offset_with_breakstick_holes(object_base, 0.5, 10, 0.3);

            // We can cut a hole in the pad corresponding to the object shape:
            inner_base.holes.emplace_back(object_base);
            inner_base.holes.back().reverse();

            if(wingheight <= 0) {
                // if no wings then cut the hole in the upper plate as well
                top_poly.holes.emplace_back(object_base);
                top_poly.holes.back().reverse();
            }

            triangulate(top_poly, top_triangles);
            triangulate(inner_base, bottom_triangles);

            auto lines = object_base.lines();

            // Generate outer walls
            auto fp = [](const Point& p, Point::coord_type z) {
                return unscale(x(p), y(p), z);
            };

            for(auto& l : lines) {
                auto s = coord_t(pool.points.size());

                pool.points.emplace_back(fp(l.a, -s_thickness - s_wingheight));
                pool.points.emplace_back(fp(l.b, -s_thickness - s_wingheight));
                pool.points.emplace_back(fp(l.a, -s_wingheight));
                pool.points.emplace_back(fp(l.b, -s_wingheight));

                pool.indices.emplace_back(s + 3, s + 1, s);
                pool.indices.emplace_back(s + 2, s + 3, s);
            }
        }

        auto top_plate = convert(top_triangles, 0, false);
        auto bottom_plate = convert(bottom_triangles, -mm(fullheight), true);

        pool.merge(top_plate);
        pool.merge(bottom_plate);

        if(wingheight > 0) {
            Polygons middle_triangles;
            triangulate(inner_base, middle_triangles);
            auto middle_plate = convert(middle_triangles, -mm(wingheight), false);
            pool.merge(middle_plate);
        }

        out.merge(mesh(pool));
    }
}

void create_base_pool(const ExPolygons &base_plate,
                      TriangleMesh &output_mesh,
                      const PoolConfig & cfg)
{
    create_base_pool(base_plate, {}, output_mesh, cfg);
}

}
}
