      program matinv

      implicit none

      include 'add.inc'
      include 'mode.inc'
      include 'grid.inc'

      integer :: I,J,K,N,MUP,info,L
      double precision :: R(NR),RR(NR),POTU(NR),R0,NH,work
      parameter(mup=nr*nlu)
      double complex one,zero,ione,hamup(mup,mup),haminv(mup,mup)
      double complex Hfinal(mup,mup),identity(mup,mup),MAINDIAGU(MUP),
     &                  offdiagu(mup)
      parameter(IONE=(0.D0,1.D0),ONE=(1.D0,0.D0),ZERO=(0.D0,0.D0))


      DO I=1,NR
         R(I) = (DBLE(I)-N0)*DR
      ENDDO

      DO I=1,NR
         RR(I) = DEXP(R(I))
      ENDDO
             
      DO I=1,NR
         POTU(I) = 0.5D0*RR(I)**2
      ENDDO

      DO J=1,NLU
         DO I=1,NR
            MAINDIAGU((J-1)*NR+I)=ZERO
         enddo       
      ENDDO


      NH = (NLU+1)/2
        DO J=1,NLU
           L=J-NH

           DO I=1,NR
               MAINDIAGU((J-1)*NR+I) =
     &              (1.D0/DR**2+L**2*0.5D0)*DEXP(-2.D0*R(I)) + POTU(I)
           enddo
 

                   MAINDIAGU((J-1)*NR+1) = MAINDIAGU((J-1)*NR+1)
     &                           - (0.5D0/DR**2)*DEXP(-2.D0*R(1))
     &                           * DEXP(-ABS(L)*(R(1)-R0))

        ENDDO

        DO J=1,NLU

           OFFDIAGU((J-1)*NR+NR)= ZERO

           DO I=1,NR-1
              OFFDIAGU((J-1)*NR+I) = -(0.5D0/DR**2)*DEXP(-R(I)-R(I+1))
           enddo
        ENDDO

      DO I=1,MUP
         DO K=1,MUP
            HAMUP(I,K) = ZERO
         enddo
      ENDDO

      DO I=1,MUP
         HAMUP(I,I) = MAINDIAGU(I)*IONE*0.5D0*DT
      ENDDO

      DO I=1,MUP-1
         HAMUP(I,I+1) = OFFDIAGU(I)*IONE*0.5D0*DT
      ENDDO
       
      do i=1,mup
         do k=1,mup
            Hfinal(i,k) = identity(i,k)-ione*hamup(i,k)/2.d0
         enddo
      enddo
  
      do i=1,mup
         do k=1,mup
            haminv(i,k) = identity(i,k)+ione*hamup(i,k)/2.d0
         enddo
      enddo  
 
      hamup = hfinal/haminv

      write(99,*) dble(hamup)

      END
