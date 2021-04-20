
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

using namespace dealii;


int
main()
{
  const unsigned int dim    = 2;
  const unsigned int degree = 1;

  FE_Q<dim>      fe(degree);
  QGauss<dim>    quad(degree + 1);
  MappingQ1<dim> mapping;

  FEValues<dim> fe_values(mapping, fe, quad, update_values | update_JxW_values);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  FullMatrix<double> matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());

  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values.reinit(cell);

      for (const auto i : fe_values.dof_indices())
        for (const auto j : fe_values.dof_indices())
          for (const auto q : fe_values.quadrature_point_indices())
            matrix(i, j) += fe_values.shape_value(i, q) *
                            fe_values.shape_value(j, q) * fe_values.JxW(q);
    }

  matrix.print_formatted(std::cout);
}
