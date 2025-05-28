program test_vtk
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  type(iu_grid_t)    :: ug

  call iu_read_grid('test_data/triangle.binda', ug)
  call iu_write_vtk(ug, 'test_data/test_vtk_triangle.vtu')

  call iu_read_grid('test_data/quad.binda', ug)
  call iu_write_vtk(ug, 'test_data/test_vtk_quad.vtu')

  call iu_read_grid('test_data/tetra.binda', ug)
  call iu_write_vtk(ug, 'test_data/test_vtk_tetra.vtu')

end program test_vtk
