#include <deal.II/grid/tria.h>

// TODO: include the right files

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
    // TODO: read mesh with GridIn
    (void)tria;
  }

  {
    std::ofstream output_1(grid_out_file_name);

    // TODO: write mesh with GridOut in VTK format
    (void)tria;
    (void)output_1;
  }

  {
    AssertDimension(tria.get_reference_cells().size(), 1);
    const auto &mapping = tria.get_reference_cells()[0]
                            .template get_default_linear_mapping<dim, dim>();

    std::ofstream output_2(grid_out_file_name);

    // TODO: write mesh with DataOut in VTK format
    (void)tria;
    (void)mapping;
    (void)output_2;
  }
}