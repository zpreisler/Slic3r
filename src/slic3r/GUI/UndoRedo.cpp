#include "UndoRedo.hpp"

namespace Slic3r {
    
UndoRedo::ChangeTransformation::ChangeTransformation(const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans)
    : m_old_trans(old_trans), m_new_trans(new_trans) {}

UndoRedo::ChangeTransformation::ChangeTransformation(ModelInstance* inst, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans)
    : ChangeTransformation(old_trans, new_trans)
{
    m_instance = inst;
}

UndoRedo::ChangeTransformation::ChangeTransformation(ModelVolume* vol, const Geometry::Transformation& old_trans, const Geometry::Transformation& new_trans)
    : ChangeTransformation(old_trans, new_trans)
{
    m_volume = vol;
}

void UndoRedo::ChangeTransformation::redo()
{
    if (m_instance)
        m_instance->set_transformation(m_new_trans);
    else
        m_volume->set_transformation(m_new_trans);
}

void UndoRedo::ChangeTransformation::undo()
{
    if (m_instance)
        m_instance->set_transformation(m_old_trans);
    else
        m_volume->set_transformation(m_old_trans);
}


void UndoRedo::begin_batch(const std::string& desc)
{
    if (m_batch_start || m_batch_running)
        throw std::runtime_error("UndoRedo does not allow nested batches.");
    m_batch_desc = desc;
    m_batch_start = true;
}


void UndoRedo::end_batch() {
    if (!m_batch_start && !m_batch_running)
        throw std::runtime_error("UndoRedo batch closed when none was started.");

    m_batch_start = false;
    m_batch_running = false;
}

void UndoRedo::begin(ModelInstance* inst)
{
    if (m_current_command_type != CommandType::None)
        throw std::runtime_error("Undo/Redo stack - attemted to nest commands.");
    
    const ModelObject* mo = inst->get_object();
    unsigned int i=0;
    for (; i<mo->instances.size(); ++i)
        if (mo->instances[i] == inst)
            break;
    m_model_instance_data.mi_idx = i;
    m_model_instance_data.inst_num = mo->instances.size();
    
    i=0;
    for (; i<m_model->objects.size(); ++i)
        if (m_model->objects[i] == mo)
            break;
    m_model_instance_data.mo_idx = i;
    m_model_instance_data.transformation = inst->get_transformation();
    m_current_command_type = CommandType::ModelInstanceManipulation;
    std::cout << "UndoRedo::begin(ModelInstance*) konci..." << std::endl;
    std::cout << m_model_instance_data.mi_idx << " " << m_model_instance_data.mo_idx << std::endl;
}


void UndoRedo::end() 
{
    std::cout << "undoRedo::end()" << std::endl;
    if (m_current_command_type == CommandType::None)
        throw std::runtime_error("Undo/Redo stack - unexpected command end.");
        
    if (m_current_command_type == CommandType::ModelInstanceManipulation) {
        if (m_model->objects[m_model_instance_data.mo_idx]->instances.size() != m_model_instance_data.inst_num) {
            // the instance was deleted

        }
        else {
            // the transformation matrix was changed
            push(new ChangeTransformation(m_model->objects[m_model_instance_data.mo_idx]->instances[m_model_instance_data.mi_idx],
                                          m_model_instance_data.transformation,
                                          m_model->objects[m_model_instance_data.mo_idx]->instances[m_model_instance_data.mi_idx]->get_transformation()));
        }
    }
    else
        std::cout << "blby CT" << std::endl;

    m_current_command_type = CommandType::None;
}

void UndoRedo::push(Command* command) {
    std::cout << "UndoRedo::push" << std::endl;
    if (m_batch_running)
        command->bound_to_previous = true;
    if (m_batch_start)
        m_batch_running = true;
    if (m_batch_running)
        command->description = m_batch_desc;

     m_stack.resize(m_index); // clears the redo part of the stack
     m_stack.emplace_back(command);
     m_index = m_stack.size();
}


void UndoRedo::undo()
{
    do {
        if (!anything_to_undo())
            return;

        --m_index;
        m_stack[m_index]->undo();
    } while (m_stack[m_index]->bound_to_previous);
}

void UndoRedo::redo()
{
    do {
        if (!anything_to_redo())
            return;

        m_stack[m_index]->redo();
        ++m_index;
    } while (m_index < m_stack.size() && m_stack[m_index]->bound_to_previous);
}

} // namespace Slic3r