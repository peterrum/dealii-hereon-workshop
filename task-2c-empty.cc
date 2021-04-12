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

  template <int dim, int spacedim>
  void
  create_reentrant_corner(Triangulation<dim, spacedim> &tria)
  {
    const unsigned int n_refinements = 1;

    std::vector<unsigned int> repetitions(dim, 2);
    Point<dim>                bottom_left, top_right;
    for (unsigned int d = 0; d < dim; ++d)
      {
        bottom_left[d] = -1.;
        top_right[d]   = 1.;
      }
    std::vector<int> cells_to_remove(dim, 1);
    cells_to_remove[0] = -1;

    GridGenerator::subdivided_hyper_L(
      tria, repetitions, bottom_left, top_right, cells_to_remove);

    tria.refine_global(n_refinements);
  }

} // namespace util


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  const unsigned int dim                  = 2; // only 2D is working for now
  const unsigned int degree               = 3;
  const unsigned int n_global_refinements = 3;

  // create mesh, select relevant FEM ingredients, and set up DoFHandler
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(n_global_refinements);

  FE_Q<dim>            fe(degree);
  QGauss<dim>          quad(degree + 1);
  QGauss<dim - 1>      quad_face(degree + 1);
  MappingQGeneric<dim> mapping(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Create constraint matrix
  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  // initialize vectors and system matrix
  LinearAlgebra::distributed::Vector<double> x, b;
  TrilinosWrappers::SparseMatrix             A;
  util::initialize_dof_vector(dof_handler, x);
  util::initialize_dof_vector(dof_handler, b);
  util::initialize_system_matrix(dof_handler, constraints, A);

  // assemble right-hand side and system matrix
  FEValues<dim>     fe_values(mapping,
                          fe,
                          quad,
                          update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   quad_face,
                                   update_values | update_JxW_values);

  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

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
      for (const auto q : fe_values.quadrature_point_indices())
        {
          for (const auto i : fe_values.dof_indices())
            for (const auto j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q) * // grad phi_j(x_q)
                 fe_values.JxW(q));           // dx

          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q) * // phi_i(x_q)
                            1. *                          // f(x_q)
                            fe_values.JxW(q));            // dx
        }

      for (const auto &face : cell->face_iterators())
        {
          if (face->at_boundary() == false || face->boundary_id() != 1)
            continue;

          fe_face_values.reinit(cell, face);

          for (const auto q : fe_face_values.quadrature_point_indices())
            for (const unsigned int i : fe_face_values.dof_indices())
              cell_rhs(i) += (fe_face_values.shape_value(i, q) * // phi_i(x_q)
                              1. *                               // f(x_q)
                              fe_face_values.JxW(q));            // dx
        }

      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, A, b);
    }

  b.compress(VectorOperation::values::add);
  A.compress(VectorOperation::values::add);

  // solve linear equation system
  ReductionControl reduction_control(100, 1e-10, 1e-4);
  SolverCG<LinearAlgebra::distributed::Vector<double>> solver(
    reduction_control);
  solver.solve(A, x, b, PreconditionIdentity());

  if (Utilities::MPI::this_mpi_process(util::get_mpi_comm(tria)) == 0)
    printf("Solved in %d iterations.\n", reduction_control.last_step());

  constraints.distribute(x);

  // output results
  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);
  x.update_ghost_values();
  data_out.add_data_vector(dof_handler, x, "solution");
  data_out.build_patches(mapping, degree + 1);
  data_out.write_vtu_with_pvtu_record("./",
                                      "result",
                                      0,
                                      util::get_mpi_comm(tria));
}