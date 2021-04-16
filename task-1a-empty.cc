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
