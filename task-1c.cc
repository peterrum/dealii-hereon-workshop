#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

const int dim = 2;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD, Triangulation<dim>::MeshSmoothing::none, true);

  {
    GridIn<dim> grid_in(tria);
    grid_in.read("beam.msh");
  }

  {
    GridOut grid_out;
    grid_out.write_mesh_per_processor_as_vtu(tria, "task-1c");
  }
}
