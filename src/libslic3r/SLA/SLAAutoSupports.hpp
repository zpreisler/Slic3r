#ifndef SLAAUTOSUPPORTS_HPP_
#define SLAAUTOSUPPORTS_HPP_

#include <chrono>

#include <libslic3r/Point.hpp>
#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/SLA/SLASupportTree.hpp>

 #define SLA_AUTOSUPPORTS_DEBUG

namespace Slic3r {

class SLAAutoSupports {
public:
    struct Config {
            float density_at_horizontal;
            float density_at_45;
            float minimal_z;
        };

    SLAAutoSupports(const TriangleMesh& mesh, const sla::EigenMesh3D& emesh, const std::vector<ExPolygons>& slices, 
        const std::vector<float>& heights, const Config& config, std::function<void(void)> throw_on_cancel);
    const std::vector<Vec3d>& output() { return m_output; }

private:
    std::vector<Vec3d> m_output;
    std::vector<Vec3d> m_normals;
    static float angle_from_normal(const stl_normal& normal) { return acos((-normal.normalized())(2)); }
    float get_required_density(float angle) const;
    float distance_limit(float angle) const;
    static float approximate_geodesic_distance(const Vec3d& p1, const Vec3d& p2, Vec3d& n1, Vec3d& n2);
    std::vector<std::pair<ExPolygon, coord_t>> find_islands(const std::vector<ExPolygons>& slices, const std::vector<float>& heights) const;
    void sprinkle_mesh(const TriangleMesh& mesh);
    std::vector<Vec3d> uniformly_cover(const std::pair<ExPolygon, coord_t>& island);
    void project_upward_onto_mesh(std::vector<Vec3d>& points) const;

#ifdef SLA_AUTOSUPPORTS_DEBUG
    void output_expolygons(const ExPolygons& expolys, std::string filename) const;
#endif /* SLA_AUTOSUPPORTS_DEBUG */

    SLAAutoSupports::Config m_config;
    std::function<void(void)> m_throw_on_cancel;
    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


class SLAAutoSupports2 {
public:
    struct Config {
            float density_at_horizontal;
            float density_at_45;
            float minimal_z;
            ///////////////
            float support_force = 1; // a force unit
            float tear_pressure = 1; // the force unit per mm2
        };

    SLAAutoSupports2(const TriangleMesh& mesh, const sla::EigenMesh3D& emesh, const std::vector<ExPolygons>& slices,
                     const std::vector<float>& heights, const Config& config, std::function<void(void)> throw_on_cancel);
    const std::vector<Vec3d>& output() { return m_output; }

private:
    std::vector<Vec3d> m_output;

#ifdef SLA_AUTOSUPPORTS_DEBUG
    void output_expolygons(const ExPolygons& expolys, std::string filename) const;
    void output_structures() const;
#endif /* SLA_AUTOSUPPORTS_DEBUG */

    SLAAutoSupports2::Config m_config;

    struct Structure {
        Structure(const ExPolygon& polygon, coord_t h) : height(h), supports_force(0.f), unique_id(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch())) {
            top_slice.push_back(polygon);
            new_slice.push_back(polygon);
        }
        ExPolygons top_slice;
        ExPolygons new_slice;
        coord_t height;
        Point supports_force_point;
        float supports_force;
        std::chrono::milliseconds unique_id;
    };
    std::vector<Structure> m_structures;

    std::vector<std::pair<ExPolygon, coord_t>> find_islands(const std::vector<ExPolygons>& slices, const std::vector<float>& heights);
    void polygon_to_structures(const ExPolygon& polygon, coord_t height, std::vector<unsigned int>& indices);
    void delete_finished_structures(coord_t height);


    std::function<void(void)> m_throw_on_cancel;
    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;
};


} // namespace Slic3r


#endif // SLAAUTOSUPPORTS_HPP_