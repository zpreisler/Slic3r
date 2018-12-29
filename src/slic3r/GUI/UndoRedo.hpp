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
    };
    
    
    class ChangeTransformation : public Command {
    public:
        ChangeTransformation(const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans);
        ChangeTransformation(ModelInstance* inst, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans);
        ChangeTransformation(ModelVolume* vol, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans);
        
        void redo() override;
        void undo() override;

    private:
        ModelInstance* m_instance = nullptr;
        ModelVolume* m_volume = nullptr;
        Geometry::Transformation m_old_trans;
        Geometry::Transformation m_new_trans;
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
    
    void begin();                // new ModelObject is about to be created
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
    
    std::string get_undo_description() const {
        return m_stack[m_index-1]->description;
    }
    std::string get_redo_description() const {
        return m_stack[m_index]->description;
    }

private:
	std::vector<std::unique_ptr<Command>> m_stack;
	unsigned int m_index = 0; // points to a command to redo
    std::string m_batch_desc = "";
    bool m_batch_start = false;
    bool m_batch_running = false;
    CommandType m_current_command_type = CommandType::None;
    Model* m_model;
    
    struct ModelInstanceData {
        unsigned int mo_idx;
        unsigned int mi_idx;
        unsigned int inst_num;
        Geometry::Transformation transformation;
    }m_model_instance_data;
};


} // namespace Slic3r

#endif // UNDOREDO_HPP