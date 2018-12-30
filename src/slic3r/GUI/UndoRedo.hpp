#ifndef UNDOREDO_HPP
#define UNDOREDO_HPP

#include <libslic3r/Model.hpp>

#include <memory>

namespace Slic3r {


class UndoRedo {
public:
    UndoRedo(Model* model) : m_model(model) {}
    
    class Command {
    public:
        virtual void redo() = 0;
        virtual void undo() = 0;
        virtual ~Command() {}
        bool bound_to_previous = false;
        std::string description;
        Model* m_model;
    };
    
    
    class Change : public Command {
    public:
        Change(ModelInstance* inst, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans);
        Change(ModelVolume* vol, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans,
                             const std::string old_name, const std::string new_name, ModelVolume::Type old_type, ModelVolume::Type new_type);
        void redo() override;
        void undo() override;
    private:
        int m_mo_idx = -1;
        int m_mi_idx = -1;
        int m_mv_idx = -1;
        Geometry::Transformation m_old_trans;
        Geometry::Transformation m_new_trans;
        std::string m_old_name;
        std::string m_new_name;
        ModelVolume::Type m_old_type;
        ModelVolume::Type m_new_type;
    };

    class AddInstance : public Command {
    public:
        AddInstance(ModelInstance* mi, unsigned int mo_idx);
        void redo() override;
        void undo() override;
    private:
        unsigned int m_mo_idx;
        Geometry::Transformation m_trans;
    };

    class RemoveInstance : public Command {
    public:
        RemoveInstance(Model* model, Geometry::Transformation trans, unsigned int mo_idx, unsigned int mi_idx);
        void redo() override;
        void undo() override;
    private:
        unsigned int m_mo_idx;
        unsigned int m_mi_idx;
        Geometry::Transformation m_trans;
    };

    enum class CommandType {
        None,
        NewModelObject,
        ModelObjectManipulation,
        ModelInstanceManipulation,
        ModelVolumeManipulation
    };



    void begin_batch(const std::string& desc);
    void end_batch();

    void begin();               // new ModelObject is about to be created
    void begin(ModelObject*);   // ModelObject - deletion / name change / new instance / new ModelVolume / config change / layer_editing
    void begin(ModelInstance*); // ModelInstance - deletion / transformation matrix change
    void begin(ModelVolume*);   // ModelVolume - deletion / transformation change / ModelType change / name change / config change
    void end();                 // this step is finished
     
    bool anything_to_redo() const {
        return (!m_stack.empty() && m_index != m_stack.size());
    }
     
    bool anything_to_undo() const {
        return (!m_stack.empty() && m_index != 0);
    }
    
    void push(Command* command);
    
    void undo();
    void redo();
    
    bool working() const;

    std::string get_undo_description() const {
        return m_stack[m_index-1]->description;
    }
    std::string get_redo_description() const {
        return m_stack[m_index]->description;
    }

private:
	std::vector<std::unique_ptr<Command>> m_stack;
	unsigned int m_index = 0;
    std::string m_batch_desc = "";
    bool m_batch_start = false;
    bool m_batch_running = false;
    CommandType m_current_command_type = CommandType::None;
    bool m_lock = false;
    Model* m_model;
    
    struct ModelInstanceData {
        unsigned int mo_idx;
        unsigned int mi_idx;
        unsigned int inst_num;
        Geometry::Transformation transformation;
    }m_model_instance_data;

    struct ModelVolumeData {
        unsigned int mo_idx;
        unsigned int mv_idx;
        unsigned int vols_num;
        Geometry::Transformation transformation;
        ModelVolume::Type type;
        std::string name;
    }m_model_volume_data;

    struct ModelObjectData {
        unsigned int mo_idx;
        unsigned int inst_num;
    }m_model_object_data;
};


} // namespace Slic3r

#endif // UNDOREDO_HPP