#pragma once


//based off of ihttps://github.com/alecjacobson/gptoolbox/blob/master/mex/eltopo.cpp




#include <eltopo3d/surftrack.h>
#include <Eigen/Dense>
#include <memory>

class SurfTrack;

class ElTopoTracker {
    public:
        using EigenVec3d = Eigen::Vector3d;
        using EigenVec3i = Eigen::Vector3i;
        using ColVectors3d = Eigen::Matrix<double,3,Eigen::Dynamic>;
        using ColVectors3i = Eigen::Matrix<int,3,Eigen::Dynamic>;
        using RefCV3d= Eigen::Ref<ColVectors3d>;
        using RefCV3i = Eigen::Ref<ColVectors3i>;
        using CRefCV3d= Eigen::Ref<const ColVectors3d>;
        using CRefCV3i = Eigen::Ref<const ColVectors3i>;
        ElTopoTracker& operator=(ElTopoTracker&&) = default;
        ~ElTopoTracker();
        ElTopoTracker(const CRefCV3d& V, const CRefCV3i& F, bool collision_safety=true, bool defrag_mesh = true, bool verbose = false, bool do_initialize=true);
        static Eigen::VectorXd dual_volumes(const CRefCV3d& V, const CRefCV3i& F);

        void split_edge(size_t edge_index);
        //split_triangle splits by the longest edge
        void split_triangle(size_t triangle_index);

        void initialize();
        void improve();
        double integrate(const CRefCV3d& Vnew, double dt);
        double step(const CRefCV3d& Vnew, double dt);
        double integrate_py(const Eigen::Ref<const Eigen::MatrixXd>& Vnew, double dt) {
            return integrate(Vnew,dt);
        }
        double step_py(const Eigen::Ref<const Eigen::MatrixXd>& Vnew, double dt) {
            return step(Vnew,dt);
        }

        void defrag_mesh();

        std::tuple<ColVectors3d,ColVectors3i> get_mesh() const;


        SurfTrackInitializationParameters& init_params() { return m_init_params; }
        ColVectors3d get_vertices() const;
        ColVectors3i get_triangles() const;
        SurfTrackInitializationParameters m_init_params;
        bool m_verbose = false;
        bool m_defrag_mesh = false;
    private:
        std::vector<Vec3st> tris;
        std::vector<Vec3d> verts;
        std::vector<double> masses;
        std::unique_ptr<SurfTrack> m_surf = nullptr;
        std::unique_ptr<SubdivisionScheme> m_subdivision_scheme;
        bool m_defrag_dirty = false;

};
