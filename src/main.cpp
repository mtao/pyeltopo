#include <memory>
#include <algorithm>
#include <iterator>
#include <fstream>
#include "eltopo.h"
#include <iostream>

#include <eltopo3d/eltopo.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;


PYBIND11_MODULE(pyeltopo, m) {
    py::class_<ElTopoTracker>(m, "ElTopoTracker")
        //.def(py::init<const ElTopoTracker::CRefCV3d&, const ElTopoTracker::CRefCV3i&>())
        .def(py::init<const ElTopoTracker::CRefCV3d&, const ElTopoTracker::CRefCV3i&,bool, bool, bool,bool>(), py::arg("vertices"), py::arg("faces") ,py::arg("collision_safety")=true, py::arg("defrag_mesh")=true,py::arg("verbose") = false, py::arg("initialize")=true)
        //.def(py::init<const ElTopoTracker::CRefCV3d&, const ElTopoTracker::CRefCV3i&,bool,bool>())
        .def_readwrite("init_params", &ElTopoTracker::m_init_params)
        .def_readwrite("enable_mesh_defrag", &ElTopoTracker::m_defrag_mesh)
        .def("get_mesh",&ElTopoTracker::get_mesh)
        .def("get_triangles",&ElTopoTracker::get_triangles)
        .def("get_vertices",&ElTopoTracker::get_vertices)
        .def("integrate",&ElTopoTracker::integrate_py)
        .def("improve",&ElTopoTracker::improve)
        .def("initialize",&ElTopoTracker::initialize)
        .def("step",&ElTopoTracker::step_py)
        .def("defrag_mesh",&ElTopoTracker::defrag_mesh)
        .def("split_edge",&ElTopoTracker::split_edge)
        .def("split_triangle",&ElTopoTracker::split_triangle);
    py::class_<SurfTrackInitializationParameters>(m, "SurfTrackInitializationParameters")
        .def(py::init<>())
        .def_readwrite("proximity_epsilon",&SurfTrackInitializationParameters::m_proximity_epsilon)
        .def_readwrite("friction_coefficient",&SurfTrackInitializationParameters::m_friction_coefficient)
        .def_readwrite("min_triangle_area",&SurfTrackInitializationParameters::m_min_triangle_area)
        .def_readwrite("improve_collision_epsilon",&SurfTrackInitializationParameters::m_improve_collision_epsilon)
        .def_readwrite("use_fraction",&SurfTrackInitializationParameters::m_use_fraction)
        .def_readwrite("min_edge_length",&SurfTrackInitializationParameters::m_min_edge_length)
        .def_readwrite("max_edge_length",&SurfTrackInitializationParameters::m_max_edge_length) 
        .def_readwrite("max_volume_change",&SurfTrackInitializationParameters::m_max_volume_change)
        .def_readwrite("min_triangle_angle",&SurfTrackInitializationParameters::m_min_triangle_angle)
        .def_readwrite("max_triangle_angle",&SurfTrackInitializationParameters::m_max_triangle_angle)   
        .def_readwrite("use_curvature_when_splitting",&SurfTrackInitializationParameters::m_use_curvature_when_splitting)
        .def_readwrite("use_curvature_when_collapsing",&SurfTrackInitializationParameters::m_use_curvature_when_collapsing)
        .def_readwrite("min_curvature_multiplier",&SurfTrackInitializationParameters::m_min_curvature_multiplier)
        .def_readwrite("max_curvature_multiplier",&SurfTrackInitializationParameters::m_max_curvature_multiplier)
        .def_readwrite("allow_vertex_movement",&SurfTrackInitializationParameters::m_allow_vertex_movement)
        .def_readwrite("edge_flip_min_length_change",&SurfTrackInitializationParameters::m_edge_flip_min_length_change)
        .def_readwrite("merge_proximity_epsilon",&SurfTrackInitializationParameters::m_merge_proximity_epsilon)
        .def_readwrite("collision_safety",&SurfTrackInitializationParameters::m_collision_safety)
        .def_readwrite("allow_topology_changes",&SurfTrackInitializationParameters::m_allow_topology_changes)
        .def_readwrite("allow_non_manifold",&SurfTrackInitializationParameters::m_allow_non_manifold)
        .def_readwrite("perform_improvement",&SurfTrackInitializationParameters::m_perform_improvement);
}


