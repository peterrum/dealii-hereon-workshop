// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#include <deal.II/grid/tria.h>

// necessary includes
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

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
    std::ofstream output_2("task-1a-data.vtk");

    // write mesh with DataOut in VTK format
    DataOut<dim> data_out;
    data_out.attach_triangulation(tria);
    data_out.build_patches();
    data_out.write_vtk(output_2);
  }
}
