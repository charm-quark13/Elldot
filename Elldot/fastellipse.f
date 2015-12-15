       PROGRAM FASTELLIPSE

       IMPLICIT NONE

       DOUBLE PRECISION Wy(3,3), Wx(3,3),J(3,3)
       INTEGER I,K

       DO 10 I=1,3
       DO 10 K=1,3
          Wy(I,K) = 0.D0
          Wx(I,K) = 0.D0
10     CONTINUE
         
          Wy(1,3) = 2.D0 
         
          Wx(1,1) = 2.D0
          Wx(1,2) = 1.D0
          Wx(2,2) = 1.D0
          Wx(3,3) = 1.D0
 
       DO 20 I=1,3
          Wy(I,I) = 1.D0
20     CONTINUE

       DO 40 I=1,3
       DO 40 K=1,3
          J(I,K) = 0.D0
40     CONTINUE
       
       J=MATMUL(Wx,Wy)  

       WRITE(1111,*)J 
       END
