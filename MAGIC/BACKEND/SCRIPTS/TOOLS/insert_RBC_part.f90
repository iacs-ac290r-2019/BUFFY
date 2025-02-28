MODULE all_mod

USE ISO_C_BINDING

IMPLICIT NONE

INTEGER :: nx,ny,nz
INTEGER :: nx2 = 0, nxy2 = 0

CONTAINS

INTEGER FUNCTION ijk(idir, i4)
INTEGER :: idir,iy,iz
INTEGER(C_LONG) :: i4

      iz = i4/INT(nxy2)
      IF(idir==3) THEN
         ijk = INT(iz)
         RETURN
      ENDIF

      iy = (i4 - iz*INT(nxy2))/INT(nx2)
      IF(idir==2) THEN
         ijk = INT(iy)
         RETURN
      ENDIF

      ijk = INT(i4 - iz*INT(nxy2) - iy*INT(nx2))


ENDFUNCTION

FUNCTION random_ijk(itp, n) RESULT(out)
INTEGER(C_LONG) :: itp(0:)
INTEGER :: n
INTEGER :: nrnd
REAL :: out(3)
REAL :: rnd
INTEGER(C_LONG) :: i4
INTEGER :: i,j,k

    CALL RANDOM_NUMBER(rnd)
    nrnd = INT(rnd * SIZE(itp))
    i4 = itp(nrnd)
    i = ijk(1,i4)
    j = ijk(2,i4)
    k = ijk(3,i4)
    print *,'nrnd:',nrnd,i,j,k
    out = [i,j,k]

ENDFUNCTION

REAL FUNCTION dist2(r1, r2, ialign)
    REAL,INTENT(in) :: r1(:), r2(:)
    INTEGER,INTENT(in) :: ialign
    REAL :: dx,dy,dz

    dx = r1(1) - r2(1)
    dy = r1(2) - r2(2)
    dz = r1(3) - r2(3)

    IF      (ialign==1) THEN
        dx = dx - NINT(dx/NX)*NX

    ELSE IF (ialign==2) THEN
        dy = dy - NINT(dy/NY)*NY

    ELSE IF (ialign==3) THEN
        dz = dz - NINT(dz/NZ)*NZ

    ENDIF

    dist2 = dx**2 + dy**2 + dz**2

ENDFUNCTION

ENDMODULE

SUBROUTINE preproc_part(ialign, &
                        nx_l, ny_l, nz_l, &
                        nfl, nwl, &
                        nrbc, nsph, &
                        cutoff_insertion, cutoff_wall, &
                        ii_f,jj_f,kk_f, &
                        ii_w,jj_w,kk_w, &
                        xo, yo, zo) BIND(C,NAME='preproc_part')

USE all_mod

IMPLICIT NONE

INTEGER(C_INT),INTENT(in),VALUE :: ialign
INTEGER(C_INT),INTENT(in),VALUE :: nx_l,ny_l,nz_l,nfl,nwl,nrbc,nsph
REAL(C_FLOAT),INTENT(in),VALUE :: cutoff_insertion, cutoff_wall
INTEGER(C_LONG),DIMENSION(0:nfl-1) :: ii_f, jj_f, kk_f
INTEGER(C_LONG),DIMENSION(0:nfl-1) :: ii_w, jj_w, kk_w
! REAL(C_FLOAT),DIMENSION(1:nfl) :: xo, yo, zo
REAL(C_DOUBLE),DIMENSION(1:nfl) :: xo, yo, zo

! REAL :: cy,cz
REAL :: d2
LOGICAL :: tooclose
REAL :: rr(3), ro(3), rw(3)
INTEGER :: nold,i,j,k,nrnd
REAL :: cut2,cutw2
REAL :: rnd
INTEGER :: ib
INTEGER,ALLOCATABLE,SAVE :: ncx(:),indx(:,:)

    CALL RANDOM_SEED()

    nx = nx_l
    ny = ny_l
    nz = nz_l
    nx2 = nx + 2
    nxy2 = (nx + 2) * (ny + 2)

#if defined(LINKCELL)
    IF     (ialign==1) THEN ! x-axis
        print *,'x-alignment'
        ALLOCATE(ncx(1:nx), indx(1:nx, 1:nfl))

    ELSE IF(ialign==2) THEN ! y-axis
        print *,'y-alignment'
        ALLOCATE(ncx(1:ny), indx(1:ny, 1:nfl))

    ELSE IF(ialign==3) THEN ! z-axis
        print *,'z-alignment'
        ALLOCATE(ncx(1:nz), indx(1:nz, 1:nfl))

    ENDIF
    ncx(:) = 0
#endif

    print *,'...Fast particle placement...'
    !print *,'VALS:',nx_l,ny_l,nz_l,nfl,nwl,nrbc,nsph
    print *,'nx,ny,nz:',nx,ny,nz
    print *,'part-part cutoff:',cutoff_insertion,'part-wall cutoff:',cutoff_wall
    !print *,ii_f(0:10)
    !print *,jj_f(0:10)
    !print *,kk_f(0:10)

    cut2 = cutoff_insertion**2
    cutw2 = cutoff_wall**2

    nold = 0
    ! cy = NY/2.; cz = NZ/2.
    DO WHILE (.true.)

        CALL RANDOM_NUMBER(rnd)

        ! nrnd = INT(rnd * SIZE(ii_f))
        nrnd = INT(rnd * nfl)

        rr = [ii_f(nrnd), jj_f(nrnd), kk_f(nrnd)]
        ! print *,'rr:',rr

        tooclose = .false.

#if 0
        rw2 = radius - SQRT(cutw)
        d2 = (rr(2) - cy)**2 + (rr(3) - cz)**2
        tooclose = (d2 > rw2) 
#else
        DO i = 1, SIZE(ii_w)
            !print *, i, ii_w(i), jj_w(i), kk_w(i)
            rw(1) = ii_w(i)
            rw(2) = jj_w(i)
            rw(3) = kk_w(i)
            IF ( dist2(rw, rr, ialign) < cutw2 ) THEN
                tooclose = .true.
                ! print *,'too close to wall', dist2(ro, rr, ialign)
                EXIT
            ENDIF
        ENDDO
#endif

        IF (.NOT. tooclose) THEN

#if !defined(LINKCELL)
            DO i = 1, nold
                ro = [xo(i), yo(i), zo(i)]
                IF ( dist2(ro, rr, ialign) < cut2 ) THEN
                    tooclose = .true.
                    EXIT
                ENDIF
            ENDDO
#else
            ! little link cell
            DO i = ii_f(nrnd) - 1, ii_f(nrnd) + 1 
              ib = i

              IF      (ialign==1) THEN
                IF ( ib < 1 )  ib = ib + nx
                IF ( ib > nx ) ib = ib - nx

              ELSE IF (ialign==2) THEN
                IF ( ib < 1 )  ib = ib + ny
                IF ( ib > ny ) ib = ib - ny

              ELSE IF (ialign==3) THEN
                IF ( ib < 1 )  ib = ib + nz
                IF ( ib > nz ) ib = ib - nz

              ENDIF

              DO j = 1, ncx(ib)

                k = indx(ib,j)
                ro = [xo(k), yo(k), zo(k)]

                IF ( dist2(ro, rr, ialign) < cut2 ) THEN
                    tooclose = .true.
                    EXIT
                ENDIF
              ENDDO
            ENDDO
#endif
        ENDIF
    
        IF (.NOT. tooclose) THEN
            IF(MOD(nold,1000)==0) print *,'accept nold:',nold,'/',nrbc + nsph, &
                                           INT(nold*100./(nrbc + nsph)),'%'

            nold = nold + 1
            ! rro(:,nold) = rr(:)
            xo(nold) = rr(1)
            yo(nold) = rr(2)
            zo(nold) = rr(3)

#if defined(LINKCELL)
            ib = rr(ialign)
            ncx(ib) = ncx(ib) + 1
            indx(ib,ncx(ib)) = nold
#endif

            IF (nold >= nrbc + nsph) EXIT
        ENDIF

    ENDDO

    ! xo(1:nold) = rro(1,1:nold)
    ! yo(1:nold) = rro(2,1:nold)
    ! zo(1:nold) = rro(3,1:nold)

    !DO i=1, 10 ! SIZE(xo)
    !    print *,'xo_F:',i,xo(i)
    !ENDDO

#if 0
    ! ctypes does not like to open files...crash
    OPEN(UNIT=1000, STATUS='unknown', file='FettaSalame.xyz')
    WRITE(1000,*) nold
    WRITE(1000,*) 
    DO i=1,nold
        print *,i,'...write'
        WRITE(1000,*) 'O', xo(i), yo(i), zo(i)
    ENDDO
    CLOSE(1000)
    STOP
#endif

ENDSUBROUTINE
