// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

namespace util
{
  /**
   * Extract communicator of @p mesh.
   */
  template <typename MeshType>
  MPI_Comm
  get_mpi_comm(const MeshType &mesh)
  {
    const auto *tria_parallel = dynamic_cast<
      const parallel::TriangulationBase<MeshType::dimension,
                                        MeshType::space_dimension> *>(
      &(mesh.get_triangulation()));

    return tria_parallel != nullptr ? tria_parallel->get_communicator() :
                                      MPI_COMM_SELF;
  }

  /**
   * Create a deal.II partitioner.
   */
  template <int dim, int spacedim>
  std::shared_ptr<const Utilities::MPI::Partitioner>
  create_dealii_partitioner(
    const DoFHandler<dim, spacedim> &dof_handler,
    unsigned int                     mg_level = numbers::invalid_unsigned_int)
  {
    IndexSet locally_relevant_dofs;

    if (mg_level == numbers::invalid_unsigned_int)
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
    else
      DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                    mg_level,
                                                    locally_relevant_dofs);

    return std::make_shared<const Utilities::MPI::Partitioner>(
      mg_level == numbers::invalid_unsigned_int ?
        dof_handler.locally_owned_dofs() :
        dof_handler.locally_owned_mg_dofs(mg_level),
      locally_relevant_dofs,
      get_mpi_comm(dof_handler));
  }

  /**
   * Initialize dof vector @p vec.
   */
  template <int dim, typename Number>
  void
  initialize_dof_vector(const DoFHandler<dim> &                     dof_handler,
                        LinearAlgebra::distributed::Vector<Number> &vec)
  {
    const auto partitioner_dealii = create_dealii_partitioner(dof_handler);
    vec.reinit(partitioner_dealii);
  }

  /**
   * Initialize dof vector @p vec.
   */
  template <int dim>
  void
  initialize_system_matrix(const DoFHandler<dim> &          dof_handler,
                           const AffineConstraints<double> &constraints,
                           TrilinosWrappers::SparseMatrix & system_matrix)
  {
    TrilinosWrappers::SparsityPattern dsp(dof_handler.locally_owned_dofs(),
                                          get_mpi_comm(dof_handler));
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    dsp.compress();

    system_matrix.reinit(dsp);
  }


  /**
   * Return stress-strain tensor.
   */
  template <int dim>
  SymmetricTensor<4, dim>
  get_stress_strain_tensor(const double lambda, const double mu)
  {
    SymmetricTensor<4, dim> tmp;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
                               ((i == l) && (j == k) ? mu : 0.0) +
                               ((i == j) && (k == l) ? lambda : 0.0));
    return tmp;
  }


  /**
   * Compute strain.
   */
  template <int dim>
  inline SymmetricTensor<2, dim>
  get_strain(const FEValues<dim> &fe_values,
             const unsigned int   shape_func,
             const unsigned int   q_point)
  {
    SymmetricTensor<2, dim> tmp;

    for (unsigned int i = 0; i < dim; ++i)
      tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i + 1; j < dim; ++j)
        tmp[i][j] =
          (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
           fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
          2;

    return tmp;
  }

} // namespace util


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  const unsigned int dim           = 2; // only 2D is working for now
  const unsigned int degree        = 1;
  const unsigned int n_refinements = 0;
  const unsigned int n_patches     = 3;

  // create mesh, select relevant FEM ingredients, and set up DoFHandler
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_rectangle(
    tria, {10, 2}, Point<dim>(0, 0), Point<dim>(1, 0.2), true);
  tria.refine_global(n_refinements);

  FESystem<dim>        fe(FE_Q<dim>(degree), dim);
  QGauss<dim>          quad(degree + 1);
  QGauss<dim - 1>      face_quad(degree + 1);
  MappingQGeneric<dim> mapping(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Create constraint matrix
  AffineConstraints<double> constraints;

  // ... fill constraint matrix
  if (false)
    {
      VectorTools::interpolate_boundary_values(
        dof_handler,
        0, // left face
        Functions::ConstantFunction<dim>(std::vector<double>{0.0, 0.0}),
        constraints,
        ComponentMask(std::vector<bool>{true, false}));

      VectorTools::interpolate_boundary_values(
        dof_handler,
        2, // bottom face
        Functions::ConstantFunction<dim>(std::vector<double>{0.0, 0.0}),
        constraints,
        ComponentMask(std::vector<bool>{false, true}));
    }
  else
    {
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0, // left face
                                               Functions::ConstantFunction<dim>(
                                                 std::vector<double>{0.0, 0.0}),
                                               constraints);
    }

  // ... finalize constraint matrix
  constraints.close();

  // initialize vectors and system matrix
  LinearAlgebra::distributed::Vector<double> x, b;
  TrilinosWrappers::SparseMatrix             A;
  {
    util::initialize_dof_vector(dof_handler, x);
    util::initialize_dof_vector(dof_handler, b);
    util::initialize_system_matrix(dof_handler, constraints, A);
  }

  // assemble right-hand side and system matrix
  {
    FEValues<dim> fe_values(mapping,
                            fe,
                            quad,
                            update_gradients | update_JxW_values);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     face_quad,
                                     update_values | update_JxW_values);


    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    const auto stress_strain_tensor =
      util::get_stress_strain_tensor<dim>(9.695e10, 7.617e10);

    // loop over all cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        fe_values.reinit(cell);

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);

        // loop over cell dofs
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              {
                const SymmetricTensor<2, dim>
                  eps_phi_i = util::get_strain(fe_values, i, q),
                  eps_phi_j = util::get_strain(fe_values, j, q);

                cell_matrix(i, j) += (eps_phi_i *            //
                                      stress_strain_tensor * //
                                      eps_phi_j              //
                                      ) *                    //
                                     fe_values.JxW(q);       //
              }

        // loop over all cell faces and their dofs
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            if (!cell->at_boundary(face) ||
                cell->face(face)->boundary_id() != 1) // we only want to apply
              continue;                               // NBC on the right face

            fe_face_values.reinit(cell, face);

            Tensor<1, dim> force;
            force[0] = 1e9;
            force[1] = 0e9;

            for (unsigned int q = 0; q < fe_face_values.n_quadrature_points;
                 ++q)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += fe_face_values.shape_value(i, q) *
                               force[fe.system_to_component_index(i).first] *
                               fe_face_values.JxW(q);
          }


        local_dof_indices.resize(cell->get_fe().dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, A, b);
      }

    b.compress(VectorOperation::values::add);
    A.compress(VectorOperation::values::add);
  }

  // solve linear equation system
  {
    ReductionControl                                     reduction_control;
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(
      reduction_control);
    solver.solve(A, x, b, PreconditionIdentity());

    if (Utilities::MPI::this_mpi_process(util::get_mpi_comm(tria)) == 0)
      printf("Solved in %d iterations.\n", reduction_control.last_step());

    constraints.distribute(x);
  }

  // output results
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    x.update_ghost_values();
    data_out.add_data_vector(
      dof_handler,
      x,
      "solution",
      std::vector<DataComponentInterpretation::DataComponentInterpretation>(
        dim, DataComponentInterpretation::component_is_part_of_vector));
    data_out.build_patches(mapping, n_patches);
    data_out.write_vtu_with_pvtu_record("./",
                                        "result",
                                        0,
                                        util::get_mpi_comm(tria));
  }
}