#include <deal.II/grid/tria.h>

// necessary includes
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

const int         dim                = 2;
const std::string grid_in_file_name  = "";
const std::string grid_out_file_name = "";
const std::string data_out_file_name = "";

int
main()
{
  Triangulation<dim> tria;

  {
    // read mesh with GridIn
    GridIn<dim> grid_in(tria);
    grid_in.read(grid_in_file_name);
  }

  {
    std::ofstream output_1(grid_out_file_name);

    // write mesh with GridOut in VTK format
    GridOut grid_out;
    grid_out.write_vtk(tria, output_1);
  }

  {
    AssertDimension(tria.get_reference_cells().size(), 1);
    const auto &mapping = tria.get_reference_cells()[0]
                            .template get_default_linear_mapping<dim, dim>();

    std::ofstream output_2(grid_out_file_name);

    // write mesh with DataOut in VTK format
    DataOut<dim> data_out;
    data_out.attach_triangulation(tria);
    data_out.build_patches(mapping);
    data_out.write_vtk(output_2);
  }
}