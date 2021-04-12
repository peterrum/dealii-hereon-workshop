#include <deal.II/grid/tria.h>

// TODO: include the right files

#include <fstream>

using namespace dealii;

const int dim = 2;

int
main()
{
  Triangulation<dim> tria;

  {
    // TODO: read mesh with GridIn
    (void)tria;
  }

  {
    std::ofstream output_1("task-1a-grid.vtk");

    // TODO: write mesh with GridOut in VTK format
    (void)tria;
    (void)output_1;
  }

  {
    std::ofstream output_2("task-1a-data.vtk");

    // TODO: write mesh with DataOut in VTK format
    (void)tria;
    (void)output_2;
  }
}
