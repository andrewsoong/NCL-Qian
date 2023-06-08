C NCLFORTSTART
      SUBROUTINE calc_clfr(west_east_dim,south_north_dim,
     1   bottom_top_dim,
     2   PRES, RH_in, SCRa, SCRb, SCRc, SCRt)
      IMPLICIT NONE
      !Arguments
      integer west_east_dim,south_north_dim,bottom_top_dim
      real SCRa(west_east_dim,south_north_dim)
      real SCRb(west_east_dim,south_north_dim)
      real SCRc(west_east_dim,south_north_dim)
      real SCRt(west_east_dim,south_north_dim)
      real PRES(west_east_dim,south_north_dim,bottom_top_dim) 
CC  Unit: Pa
      real RH_in(west_east_dim,south_north_dim,bottom_top_dim) 
CC  Unit: %
C NCLEND
CC Diagnostics: low / mid / high cloud fraction
CC Code originally from MM5toGrADS
CC    ! Local variables
      integer                  :: i, j, k, kclo, kcmi, kchi

      SCRa = -9999.
      SCRb = -9999.
      SCRc = -9999.
      SCRt = -9999.
      DO j = 1,south_north_dim   
        DO i = 1,west_east_dim
          DO k = 1,bottom_top_dim
            IF ( PRES(i,j,k) .gt. 97000. ) kclo=k
            IF ( PRES(i,j,k) .gt. 80000. ) kcmi=k
            IF ( PRES(i,j,k) .gt. 45000. ) kchi=k
          END DO
          DO k = 1,bottom_top_dim
CC  low cloud
            IF ( k .ge. kclo .AND. k .lt. kcmi )
     &          SCRa(i,j) = AMAX1(RH_in(i,j,k),SCRa(i,j))
CC  mid cloud
            IF ( k .ge. kcmi .AND. k .lt. kchi )
     &         SCRb(i,j) = AMAX1(RH_in(i,j,k),SCRb(i,j))
CC  high cloud
            IF ( k .ge. kchi )
     &         SCRc(i,j) = AMAX1(RH_in(i,j,k),SCRc(i,j))
          END DO

          SCRa(i,j)=4.0*SCRa(i,j)/100.-3.0
          SCRb(i,j)=4.0*SCRb(i,j)/100.-3.0
          SCRc(i,j)=2.5*SCRc(i,j)/100.-1.5

          SCRa(i,j)=amin1(SCRa(i,j),1.0)
          SCRa(i,j)=amax1(SCRa(i,j),0.0)
          SCRb(i,j)=amin1(SCRb(i,j),1.0)
          SCRb(i,j)=amax1(SCRb(i,j),0.0)
          SCRc(i,j)=amin1(SCRc(i,j),1.0)
          SCRc(i,j)=amax1(SCRc(i,j),0.0)

          SCRt(i,j)=1.-(1.-SCRa(i,j))*(1.-SCRb(i,j))*(1.-SCRc(i,j))

        END DO 
      END DO 


      END SUBROUTINE calc_clfr
