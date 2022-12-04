#ifort f90mtranu.f90 mod_param.f90 LJ.f90 -o soft_bending.exe
gfortran -g -fbounds-check f90mtranu.f90 mod_param.f90 mod_force.f90 LJ.f90 -o invasion3D.exe
