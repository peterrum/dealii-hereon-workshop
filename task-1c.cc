#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_in.h>

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
    AssertDimension(tria.get_reference_cells().size(), 1);
    const auto &mapping = tria.get_reference_cells()[0]
                            .template get_default_linear_mapping<dim, dim>();

    std::ofstream output(
      "task-1c." +
      std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) +
      ".vtk");

    Vector<double> ranks(tria.n_active_cells());

    for (const auto &cell : tria.active_cell_iterators())
      ranks[cell->active_cell_index()] = cell->subdomain_id();

    DataOut<dim> data_out;
    data_out.attach_triangulation(tria);
    data_out.set_cell_selection(
      [](const auto &cell) { return cell->is_artificial() == false; });
    data_out.add_data_vector(ranks, "ranks");
    data_out.build_patches(mapping, 1);
    data_out.write_vtk(output);
  }
}
