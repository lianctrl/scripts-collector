MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  TYPE, PUBLIC :: cube
    PRIVATE

    CHARACTER (LEN=72) :: str1, str2

    INTEGER :: natom, nx, ny, nz 

    INTEGER, DIMENSION(:), POINTER :: zahl

    REAL(KIND=wp) :: xmin, ymin, zmin, dx, dy, dz

    REAL(KIND=wp), DIMENSION(:), POINTER :: zch, xa, ya, za, array

    ! fill
  END TYPE cube

  PUBLIC :: cube_get, &
            cube_print, &
            cube_add, &
            cube_sub, &
            cube_int, &
            z_print, &
            cube_dqz, &
            cube_path, &
            cube_drz

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE


  CONTAINS


  SUBROUTINE cube_get (mycube, infile)

    INTEGER :: i, j, k
    REAL(KIND=wp) :: junk
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (cube), INTENT(OUT) :: mycube
    
    OPEN (UNIT=11, FILE=infile, STATUS="old", ACTION="read")

    READ (UNIT=11, FMT=*) mycube%str1

    READ (UNIT=11, FMT=*) mycube%str2

    READ (UNIT=11, FMT=*) mycube%natom, mycube%xmin, mycube%ymin, &
                          mycube%zmin 

    READ (UNIT=11, FMT=*) mycube%nx, mycube%dx

    READ (UNIT=11, FMT=*) mycube%ny, junk, mycube%dy

    READ (UNIT=11, FMT=*) mycube%nz, junk, junk, mycube%dz

    ALLOCATE ( mycube%zahl(mycube%natom) ) 

    ALLOCATE ( mycube%zch(mycube%natom) ) 

    ALLOCATE ( mycube%xa(mycube%natom) ) 

    ALLOCATE ( mycube%ya(mycube%natom) ) 

    ALLOCATE ( mycube%za(mycube%natom) ) 

    ALLOCATE ( mycube%array(mycube%nx*mycube%ny*mycube%nz) )

    do i=1,mycube%natom

      READ (UNIT=11, FMT=*) mycube%zahl(i), mycube%zch(i), &
                            mycube%xa(i), mycube%ya(i), &
                            mycube%za(i)

    end do

!example how to read till the end of file
!100 CONTINUE
!
!    READ  (UNIT=11, FMT=*, END=200) mycube%array
!
!    GO TO 100
!
!200 CONTINUE

!below array allocated takes care of all
    READ  (UNIT=11, FMT=*) mycube%array

    CLOSE (11)

!    print *, mycube%nx
!    print *, mycube%nx*mycube%ny*mycube%nz
!    print *, mycube%zmin
!    print *, mycube%dx
!    print *, mycube%zahl(1)

!    do i=1,40
!
!      print *, mycube%array(i)
!
!    end do
 
!    access attributes as follows:
!    mycube%str1 = ...

  END SUBROUTINE cube_get


  SUBROUTINE cube_print(outfile,pcube)

    INTEGER :: i, j, k
    REAL(KIND=wp) :: junk
    CHARACTER(LEN=*), INTENT(OUT) :: outfile
    TYPE (cube), INTENT(IN) :: pcube

    junk = 0.0_wp

    OPEN (UNIT=11, FILE=outfile, STATUS="new", ACTION="write")

    WRITE (UNIT=11, FMT=*) pcube%str1

    WRITE (UNIT=11, FMT=*) pcube%str2

    WRITE (UNIT=11, FMT=*) pcube%natom, pcube%xmin, pcube%ymin, &
                           pcube%zmin 

    WRITE (UNIT=11, FMT=*) pcube%nx, pcube%dx, junk, junk 

    WRITE (UNIT=11, FMT=*) pcube%ny, junk, pcube%dy, junk

    WRITE (UNIT=11, FMT=*) pcube%nz, junk, junk, pcube%dz

    do i=1,pcube%natom

      WRITE (UNIT=11, FMT=*) pcube%zahl(i), pcube%zch(i), &
                            pcube%xa(i), pcube%ya(i), &
                            pcube%za(i)

    end do

    WRITE (UNIT=11, FMT=*) pcube%array

    CLOSE (11)

  END SUBROUTINE cube_print


  FUNCTION cube_add (mycube1, mycube2)

    INTEGER :: i,j,k
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    
    cube_add%str1 = "A + B"
    cube_add%str2 = "Perform sum of the first cube and second cube"

    cube_add%natom = mycube1%natom+mycube2%natom

! this only if cubes topologies are comparable
    cube_add%xmin = mycube1%xmin

    cube_add%ymin = mycube1%ymin

    cube_add%zmin = mycube1%zmin

    cube_add%nx = mycube1%nx
    cube_add%dx = mycube1%dx

    cube_add%ny = mycube1%ny
    cube_add%dy = mycube1%dy

    cube_add%nz = mycube1%nz
    cube_add%dz = mycube1%dz
! end inconsistency if cubes are different

    ALLOCATE ( cube_add%zahl(mycube1%natom+mycube2%natom) )

    ALLOCATE ( cube_add%zch(mycube1%natom+mycube2%natom) )

    ALLOCATE ( cube_add%xa(mycube1%natom+mycube2%natom) )

    ALLOCATE ( cube_add%ya(mycube1%natom+mycube2%natom) )

    ALLOCATE ( cube_add%za(mycube1%natom+mycube2%natom) )

    ALLOCATE ( cube_add%array(cube_add%nx*cube_add%ny*cube_add%nz) )


    do i=1,mycube1%natom

      cube_add%zahl(i) = mycube1%zahl(i)
      cube_add%zch(i) = mycube1%zch(i)
      cube_add%xa(i) = mycube1%xa(i) 
      cube_add%ya(i) = mycube1%ya(i)
      cube_add%za(i) = mycube1%za(i)

    end do

    do i=1+mycube1%natom,cube_add%natom

      cube_add%zahl(i) = mycube2%zahl(i-mycube1%natom)
      cube_add%zch(i) = mycube2%zch(i-mycube1%natom)
      cube_add%xa(i) = mycube2%xa(i-mycube1%natom) 
      cube_add%ya(i) = mycube2%ya(i-mycube1%natom)
      cube_add%za(i) = mycube2%za(i-mycube1%natom)

    end do

    cube_add%array = mycube1%array + mycube2%array

!    do i=1,40
!
!      print *, cube_add%array(i)
!    end do

    ! ...
  END FUNCTION cube_add


  FUNCTION cube_sub (mycube1, mycube2)

    INTEGER :: i,j,k
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2

    cube_sub%str1 = "A - B"
    cube_sub%str2 = "Perform the subtraction of the first cube and second cube"

! control if num of atoms are equal

    cube_sub%natom = mycube1%natom

! this only if cubes topologies are comparable
    cube_sub%xmin = mycube1%xmin

    cube_sub%ymin = mycube1%ymin

    cube_sub%zmin = mycube1%zmin

    cube_sub%nx = mycube1%nx
    cube_sub%dx = mycube1%dx

    cube_sub%ny = mycube1%ny
    cube_sub%dy = mycube1%dy

    cube_sub%nz = mycube1%nz
    cube_sub%dz = mycube1%dz
! end inconsistency if cubes are different

    ALLOCATE ( cube_sub%zahl(mycube1%natom) )

    ALLOCATE ( cube_sub%zch(mycube1%natom) )

    ALLOCATE ( cube_sub%xa(mycube1%natom) )

    ALLOCATE ( cube_sub%ya(mycube1%natom) )

    ALLOCATE ( cube_sub%za(mycube1%natom) )

    ALLOCATE ( cube_sub%array(cube_sub%nx*cube_sub%ny*cube_sub%nz) )


    do i=1,mycube1%natom

      cube_sub%zahl(i) = mycube1%zahl(i)
      cube_sub%zch(i) = mycube1%zch(i)
      cube_sub%xa(i) = mycube1%xa(i) 
      cube_sub%ya(i) = mycube1%ya(i)
      cube_sub%za(i) = mycube1%za(i)

    end do

    cube_sub%array = mycube1%array - mycube2%array

!    do i=1,40
!
!      print *, cube_add%array(i)
!    end do

    ! ...
  END FUNCTION cube_sub


  FUNCTION cube_int (mycube)

    INTEGER :: i,j,k
    REAL (KIND=wp) :: cube_int, hx, hy, hz, igrl
    TYPE (cube), INTENT(IN) :: mycube

    hx = mycube%dx
    hy = mycube%dy
    hz = mycube%dz
    igrl = 0.0_wp

    do i = 0,(mycube%nx-1) 
      do j = 0,(mycube%ny-1)
        do k = 1,(mycube%nz)

!                xi = mycube%xmin + hx/2 + i*hx;
!
!                yj = mycube%ymin + hy/2 + j*hy;

!                zk = mycube%zmin + hz/2 + k*hz;

                igrl = igrl + hx*hy*hz*mycube%array(i*(mycube%ny*mycube%nz)+j*(mycube%nz)+k);
 
       end do
      end do
    end do

    cube_int = igrl;

    ! ...
  END FUNCTION cube_int


  FUNCTION cube_drz (mycube)

    INTEGER :: i,j,k                               
    REAL (KIND=wp) :: hx, hy, hz, igrl
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cube_drz
    TYPE (cube), INTENT(IN) :: mycube

    ALLOCATE(cube_drz(mycube%nz))

    hx = mycube%dx
    hy = mycube%dy
!    hz = mycube%dz
    cube_drz(:) = 0.0_wp

    do i = 1,(mycube%nz) 
      do j = 0,(mycube%nx*mycube%ny-1)
!        do k = 1,(mycube%nz)

        cube_drz(i)= cube_drz(i) + hx*hy*mycube%array((j*(mycube%nz)+i));

!        print *, cube_drz(i)
!       end do
      end do
    end do

    ! ...
  END FUNCTION cube_drz

  SUBROUTINE z_print(outfile,drz,dcube)

    INTEGER :: i, j, k
    REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: drz
    CHARACTER(LEN=*), INTENT(OUT) :: outfile
    TYPE (cube), INTENT(IN) :: dcube

    OPEN (UNIT=11, FILE=outfile, STATUS="new", ACTION="write")

    do i=1,dcube%nz

      WRITE (UNIT=11, FMT=*) (dcube%zmin + dcube%dz*i), drz(i) 

    end do

    CLOSE (11)

  END SUBROUTINE z_print

  FUNCTION cube_dqz(drz,dcube)

    REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: drz
    INTEGER :: i,j,k                               
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cube_dqz
    TYPE (cube), INTENT(IN) :: dcube

    ALLOCATE(cube_dqz(dcube%nz))
    cube_dqz(:) = 0.0_wp

    do i=1,dcube%nz
      do j=1,i

        cube_dqz(i) = cube_dqz(i) + dcube%dz*drz(j)

      end do
    end do

  END FUNCTION cube_dqz

  FUNCTION cube_path(outfile,mycube)

  CHARACTER(LEN=*), INTENT(OUT) :: outfile
  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: cube_path
  TYPE (cube), INTENT(IN) :: mycube
  INTEGER :: i,j,k                       
  REAL (KIND=wp) :: r, theta, phi

  ALLOCATE(cube_path(mycube%nz,3))


  r = SQRT(mycube%dx*mycube%dx+mycube%dy*mycube%dy+mycube%dz*mycube%dz)

  theta = ATAN(SQRT(mycube%dx*mycube%dx+mycube%dy*mycube%dy)/mycube%dz)

  phi = ATAN(mycube%dy/mycube%dx)

  do i = 1,mycube%nz

    cube_path(i,1) = mycube%xmin + i*r + i*SIN(theta) + i*COS(phi)
    cube_path(i,2) = mycube%ymin + i*r + i*COS(theta) + i*SIN(phi)
    cube_path(i,3) = mycube%zmin + i*r + i*COS(phi)

  end do

  OPEN (UNIT=11, FILE=outfile, STATUS="new", ACTION="write")
  WRITE (UNIT=11, FMT=*) cube_path 

!  do i=1,mycube%nz
!    PRINT *, cube_path(i,1)
!    PRINT *, cube_path(i,2)
!     *, cube_path(i,3)
!  end do

  END FUNCTION cube_path

!  FUNCTION cube_voro(path,mycube)
!
!  REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: path
!  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: drhos
!
!  TYPE (cube), INTENT(IN) :: mycube
!  INTEGER :: i,j,k,m, indexarr, indexs                       
!  REAL (KIND=wp) :: ds, dV, dist, r
!
!  ds = SQRT((path(1,1)-path(2,1))**2+(path(1,2)-path(2,2))**2+(path(1,3)-path(2,3))**2)
!  
!  dV = mycube%dx*mycube%dy*mycube%dz
!
!  ALLOCATE(drhos(mycube%nz))
!
!  drhos=0.0_wp
!
!  dist=1.0E10_wp
!
!  indexarr=-1
!
!  indexs=-1
!
!  do i=1,mycube%nx
!    do j=1,mycube%ny
!      do k=1,mycube%nz
!        do m=1,mycube%nz
!        
!          r = SQRT(((mycube%xmin+i*mycube%nx)-path(m,1))**2 &
!              +((mycube%ymin+j*mycube%ny)-path(m,2))**2 &
!              +((mycube%zmin+k*mycube%nx)-path(m,3))**2)
!
!          if r < dist 
!
!            dist = r
!            indexarr = i*j*k    
!            indexs = m
!
!          end if         
!         
!
!        end do
!
!        drhos(indexs)=drhos(indexs)+(mycube%array(indexarr)*(dV/ds))
!
!      end do
!    end do
!  end do
!
!  END FUNCTION cube_voro




END MODULE cubes
