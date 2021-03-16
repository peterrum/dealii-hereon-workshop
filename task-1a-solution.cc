#include <deal.II/grid/tria.h>

// necessary includes
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

const int dim = 2;

int
main()
{
  Triangulation<dim> tria;

  {
    // read mesh with GridIn
    GridIn<dim> grid_in(tria);
    grid_in.read("beam.msh");
  }

  {
    std::ofstream output_1("task-1a-grid.vtk");

    // write mesh with GridOut in VTK format
    GridOut grid_out;
    grid_out.write_vtk(tria, output_1);
  }

  {
    AssertDimension(tria.get_reference_cells().size(), 1);
    const auto &mapping = tria.get_reference_cells()[0]
                            .template get_default_linear_mapping<dim, dim>();

    std::ofstream output_2("task-1a-data.vtk");

    // write mesh with DataOut in VTK format
    DataOut<dim> data_out;
    data_out.attach_triangulation(tria);
    data_out.build_patches(mapping);
    data_out.write_vtk(output_2);
  }
}
