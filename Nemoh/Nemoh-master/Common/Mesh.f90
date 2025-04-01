!--------------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!   Contributors list:
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
MODULE MMesh

  IMPLICIT NONE

  PUBLIC

  ! A whole mesh
  TYPE TMesh
    INTEGER              :: Isym        ! Symmetry about the xOz plane (1 for yes)
    INTEGER              :: Npoints     ! Total number of points in mesh
    INTEGER              :: Npanels     ! Total number of panels in mesh
    INTEGER              :: Nbodies     ! Total number of bodies in mesh
    REAL, ALLOCATABLE    :: X(:, :)     ! Nodes coordinates
    REAL, ALLOCATABLE    :: N(:, :)     ! Normal vectors
    REAL, ALLOCATABLE    :: XM(:, :)    ! Centre of panels
    INTEGER, ALLOCATABLE :: P(:, :)     ! Connectivities
    INTEGER, ALLOCATABLE :: cPanel(:)   ! To which body belongs the panel
    REAL, ALLOCATABLE    :: A(:)        ! Area of panel
    REAL                 :: xy_diameter ! Maximal distance between two points in the mesh projected on (x, y) plane
  END TYPE TMesh

CONTAINS

  ! Operators for creation, copy, initialisation and destruction

  SUBROUTINE CreateTMesh(Mesh,Npoints,Npanels,Nbodies)
    TYPE(TMesh) :: Mesh
    INTEGER :: Npoints,Npanels,Nbodies
    Mesh%Npoints = Npoints
    Mesh%Npanels = Npanels
    Mesh%Nbodies = Nbodies
    ALLOCATE(Mesh%X(3,Mesh%Npoints),Mesh%N(3,Mesh%Npanels),Mesh%XM(3,Mesh%Npanels))
    ALLOCATE(Mesh%P(4,Mesh%Npanels),Mesh%cPanel(Mesh%Npanels),Mesh%A(Mesh%Npanels))
    Mesh%X(:,:)      = 0.0
    Mesh%N(:,:)      = 0.0
    Mesh%XM(:,:)     = 0.0
    Mesh%P(:,:)      = 0
    Mesh%cPanel(:)   = 0
    Mesh%A(:)        = 0.0
    Mesh%xy_diameter = 0.0
  END SUBROUTINE CreateTMesh


  SUBROUTINE CopyTMesh(MeshTarget,MeshSource)
    INTEGER :: i,j
    TYPE(TMesh) :: MeshTarget,MeshSource
    CALL CreateTMesh(MeshTarget,MeshSource%Npoints,MeshSource%Npanels,MeshSource%Nbodies)
    MeshTarget%Isym=MeshSource%Isym
    DO j=1,MeshTarget%Npoints
      DO i=1,3
        MeshTarget%X(i,j)=MeshSource%X(i,j)
      END DO
    END DO
    DO j=1,MeshTarget%Npanels
      DO i=1,3
        MeshTarget%N(i,j)=MeshSource%N(i,j)
        MeshTarget%XM(i,j)=MeshSource%XM(i,j)
      END DO
      DO i=1,4
        MeshTarget%P(i,j)=MeshSource%P(i,j)
      END DO
      MeshTarget%cPanel(j)=MeshSource%cPanel(j)
      MeshTarget%A(j)=MeshSource%A(j)
    END DO
    MeshTarget%xy_diameter=MeshSource%xy_diameter
  END SUBROUTINE CopyTMesh


  SUBROUTINE ReadTMesh(Mesh, mesh_directory)

    CHARACTER(LEN=*), INTENT(IN)  :: mesh_directory
    TYPE(TMesh),      INTENT(OUT) :: Mesh

    ! Local variables
    INTEGER :: i,j,k
    INTEGER :: Npoints,Npanels,Nbodies

    ! Initialize the mesh structure
    Npoints = 0
    Npanels = 0
    Nbodies = 0
    OPEN(10, FILE=TRIM(mesh_directory)//'L10.dat')
    READ(10,*)
    READ(10,*) i,Npoints,Npanels,Nbodies
    CLOSE(10)
    CALL CreateTMesh(Mesh,Npoints,Npanels,Nbodies)

    ! Read mesh
    OPEN(10, FILE=TRIM(mesh_directory)//'L12.dat')
    READ(10,*) i,Mesh%Isym
    IF (i.NE.2) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' The mesh file format is not correct. '
      WRITE(*,*) ' Stopping'
      STOP
    END IF
    DO i=1,Npoints
      READ(10,*) j,(Mesh%X(k,i),k=1,3)
    END DO
    READ(10,*)
    DO i=1,Npanels
      READ(10,*) (Mesh%P(k,i),k=1,4)
    END DO
    READ(10,*)
    CLOSE(10)

    OPEN(10, FILE=TRIM(mesh_directory)//'L10.dat')
    READ(10,*)
    READ(10,*) j
    IF (j.NE.Mesh%Isym) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' The mesh file format is not correct. '
      WRITE(*,*) ' Stopping'
      STOP
    END IF
    DO i=1,Npanels
      READ(10,*) Mesh%cPanel(i),(Mesh%XM(k,i),k=1,3),(Mesh%N(k,i),k=1,3),Mesh%A(i)
    END DO
    CLOSE(10)

    ! Search maximal distance between two points of the mesh
    Mesh%xy_diameter = 0.0
    DO i = 1, Mesh%NPoints-1
      DO j = i+1, Mesh%NPoints
        Mesh%xy_diameter = MAX(Mesh%xy_diameter, NORM2(Mesh%X(1:2, i) - Mesh%X(1:2, j)))
      END DO
      IF (Mesh%ISym == 1) THEN
        DO j = i+1, Mesh%NPoints
          Mesh%xy_diameter = MAX(Mesh%xy_diameter, NORM2(Mesh%X(1:2, i) - y_symmetry(Mesh%X(1:2, j))))
        END DO
      END IF
    END DO

  CONTAINS
    FUNCTION y_symmetry (x) result (y)
      REAL, DIMENSION(1:2) :: x, y
      y(:) = x(:)
      y(2) = -x(2)
    END FUNCTION y_symmetry
  END SUBROUTINE ReadTMesh


  SUBROUTINE DeleteTMesh(Mesh)
    TYPE(TMesh) :: Mesh
    DEALLOCATE(Mesh%X,Mesh%N,Mesh%P,Mesh%XM,Mesh%A,Mesh%cPanel)
  END SUBROUTINE DeleteTMesh



END MODULE MMesh
