MODULE paths
  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  TYPE, PUBLIC :: path

  PRIVATE
  
!  CHARACTER (LEN=72) :: str1

  INTEGER :: natom 

!  CHARCTER (LEN=*), POINTER :: zahl

  REAL(KIND=wp), DIMENSION(:), POINTER :: xa, ya, za

  
  END TYPE path

  PUBLIC :: path_get, &
            path_print
  CONTAINS


  SUBROUTINE path_get (mypath, infile)

    INTEGER :: i
    CHARACTER(LEN=1) :: junk
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (path), INTENT(OUT) :: mypath
    
    OPEN (UNIT=11, FILE=infile, STATUS="old", ACTION="read")

    READ (UNIT=11, FMT=*) mypath%natom
    PRINT *, mypath%natom

    ALLOCATE ( mypath%xa(mypath%natom) ) 

    ALLOCATE ( mypath%ya(mypath%natom) ) 

    ALLOCATE ( mypath%za(mypath%natom) ) 

!    ALLOCATE (mypath%zahl(mypath%natom))

    do i=1,mypath%natom

      READ (UNIT=11, FMT=*) junk, mypath%xa(i), mypath%ya(i), &
                            mypath%za(i)

    end do

    CLOSE (11)

!    PRINT *, mypath%xa(2)
    
  END SUBROUTINE path_get

  SUBROUTINE path_print(mypath,npoints,centers)

    INTEGER :: i,j,k                             
    INTEGER, INTENT(IN) :: npoints
    INTEGER, DIMENSION(:), INTENT(IN) :: centers
    TYPE (path), INTENT(IN) :: mypath
    REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: outpath
     
    ALLOCATE(outpath(npoints,3))

    






  END SUBROUTINE path_print

END MODULE paths
