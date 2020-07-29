PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  USE paths

  IMPLICIT NONE

  REAL(KIND=wp) :: nel
  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: drz, dqz, ciccio
  TYPE (cube) :: acube, bcube, abcube, ref, drho
  CHARACTER (LEN=72) :: outfile1, outfile2, outfile3,outfile4, outfile5
  TYPE(path) :: geom


! file where to write cube_print subroutine 
  outfile1 = "/home/urekmazino/Downloads/F95_tries/cic1819-lab/cda/aplusb.cube"
  outfile2 = "/home/urekmazino/Downloads/F95_tries/cic1819-lab/cda/drho.cube"
  outfile3 = "/home/urekmazino/Downloads/F95_tries/cic1819-lab/cda/drzplot.dat"
  outfile4 = "/home/urekmazino/Downloads/F95_tries/cic1819-lab/cda/dqzplot.dat"
  outfile5 = "/home/urekmazino/Downloads/F95_tries//cic1819-lab/cda/gen_path.dat"

! read cube files
  CALL cube_get(acube,"/home/urekmazino/Downloads/F95_tries/cic1819-lab/test/CuCO+/a.cube")
  CALL cube_get(bcube,"/home/urekmazino/Downloads/F95_tries/cic1819-lab/test/CuCO+/b.cube")
  CALL cube_get(abcube,"/home/urekmazino/Downloads/F95_tries/cic1819-lab/test/CuCO+/ab.cube")

! call function to evaluate the integral of electronic density
!  PRINT *, cube_int(abcube)
!  PRINT *, cube_int(acube)
!  PRINT *, cube_int(bcube)

  ref  = acube + bcube

!  PRINT *, cube_int(ref)

! subrout to print in .cube form
!  CALL cube_print(outfile1,ref)

  drho = abcube - ref

!  PRINT *, cube_int(drho)

! subrout to print in .cube form
!  CALL cube_print(outfile2,drho)

!  compute Delta rho(z)
  drz = cube_drz(drho)

!  CALL z_print (outfile3,drz,drho) 
!  compute dqz from drz
  dqz = cube_dqz(drz,drho)

!  CALL z_print (outfile4,dqz,drho) 
  
!  ciccio = cube_path(drho)

!  PRINT *, cube_path(outfile5,drho)

  CALL path_get(geom,"/home/urekmazino/Downloads/F95_tries/cic1819-lab/ab.xyz")



END PROGRAM cda
