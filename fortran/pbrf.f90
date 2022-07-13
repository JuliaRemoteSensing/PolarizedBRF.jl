! This code is a modified version of M. Mishchenko et al.'s old Fortran 77 code.
! The original version can be downloadable from
! https://www.giss.nasa.gov/staff/mmishchenko/brf/)
!
! Read the following paper for more details:
!
! Mishchenko, M.I., Dlugach, J.M., Chowdhary, J., Zakharova, N.T., 2015.
! Polarized bidirectional reflectance of optically thick sparse particulate
! layers: An efficient numerically exact radiative-transfer solution. Journal
! of Quantitative Spectroscopy and Radiative Transfer 156, 97–108.
! https://doi.org/10.1016/j.jqsrt.2015.02.003

module pbrf
   implicit none
   integer, private :: KONTR, NG, NPH
   integer, allocatable :: II(:)
   real(8), private :: PI, EP
   real(8), private, parameter :: DDD = 1.0d0/dsqrt(2.0d0)
   real(8), allocatable :: PH1(:,:), PH2(:,:), PH3(:,:), PH4(:,:), &
      RR(:,:),F3(:,:),XX(:),XXX(:,:)
   data KONTR/4/
   data PI/3.1415926535897932385d0/

contains
   subroutine run_pbrf(NGAUSS, LMAX1, NQUADR, EPS, ALB, AL1, AL2, &
      AL3, AL4, BET1, BET2, X, R)
      integer, intent(in) :: NGAUSS, LMAX1, NQUADR
      integer :: I, M, M1
      real(8), intent(in) :: EPS, ALB, AL1(LMAX1), AL2(LMAX1), AL3(LMAX1), &
         AL4(LMAX1), BET1(LMAX1), BET2(LMAX1)
      real(8), intent(out) :: X(NGAUSS), R(4,4,NG,NG,LMAX1)
      real(8) :: W(NGAUSS), WW(NGAUSS), X1(NGAUSS), W1(NGAUSS)

      NG=NGAUSS
      NPH=4*NG
      EP=EPS

      ! In the original version, NQUADR=1 is the same as NQUADR=2
      ! due to a bug. This has been fixed in this modified version.
      if (NQUADR==1) then
         call gauss(NG, 0, X, W)
         do I=1,NG
            X1(I)=PI*X(I)/4.0d0+PI/4.0d0
            W1(I)=PI*W(I)/4.0d0
            W1(I)=W1(I)*dsin(X1(I))
            X1(I)=dcos(X1(I))
         end do

         do I=1,NG
            X(I)=X1(NG-I+1)
            W(I)=W1(NG-I+1)
         end do
      end if

      if (NQUADR==2) then
         call gauss(NG, 1, X, W)
      end if

      if (NQUADR==3) then
         call markov(NG, 1, X, W)
      end if

      do I=1,NG
         WW(I)=dsqrt(W(I))
      end do

      allocate(PH1(NPH,NPH), PH2(NPH,NPH), PH3(NPH,NPH), PH4(NPH,NPH))
      allocate(RR(NPH,NPH), F3(NPH,NPH), XX(NG), XXX(NG,NG), II(NPH))

      do M1=1,LMAX1
         M=M1-1
         call phase(M, LMAX1, AL1, AL2, AL3, AL4, BET1, BET2, X, WW)
         if (M==0) call renorm1(M, LMAX1, AL1, X, W, WW)
         call rminf(M, LMAX1, ALB, AL1, X, WW, R)
      end do
   end subroutine

   subroutine phase(M, LMAX1, AL1, AL2, AL3, AL4, BET1, BET2, X, WW)
      integer, intent(in) :: M, LMAX1
      integer :: I, I1, J, J1, L, LMIN, LMAX, N1, N2, NN1, NN2
      real(8), intent(in) :: AL1(LMAX1), AL2(LMAX1), AL3(LMAX1), &
         AL4(LMAX1), BET1(LMAX1), BET2(LMAX1), X(NG), WW(NG)
      real(8) :: D1(LMAX1), D2(LMAX1), D3(LMAX1), S(4,4), SL(4,4,LMAX1), &
         AA(4,4), BB(4,4), P(4,LMAX1,NG),PP(4,LMAX1,NG), A1(4,4), A2(4,4), &
         A3(4,4), A4(4,4), XI, XXI, WI, DL, PN1I, PPN1I, PN2J, PPN2J, SLL

      LMIN=M+1
      LMAX=LMAX1-1
      do I=1,NG
         XI=X(I)
         XXI=-XI
         WI=WW(I)
         call dd(XI,LMAX1,M,D1,D2,D3)
         do L=LMIN,LMAX1
            DL=D1(L)*WI
            P(1,L,I)=DL
            P(2,L,I)=-D3(L)*WI
            P(3,L,I)=-D2(L)*WI
            P(4,L,I)=DL
         end do
         call dd(XXI,LMAX1,M,D1,D2,D3)
         do L=LMIN,LMAX1
            DL=D1(L)*WI
            PP(1,L,I)=DL
            PP(2,L,I)=-D3(L)*WI
            PP(3,L,I)=-D2(L)*WI
            PP(4,L,I)=DL
         end do
      end do
      S(:,:)=0.0d0
      S(1,1)=1D0
      S(2,2)=DDD
      S(2,3)=DDD
      S(3,2)=DDD
      S(3,3)=-DDD
      S(4,4)=1D0
      do L=LMIN,LMAX1
         BB(:,:)=0.0d0
         BB(1,1)=AL1(L)
         BB(1,2)=BET1(L)
         BB(2,1)=BET1(L)
         BB(2,2)=AL2(L)
         BB(3,3)=AL3(L)
         BB(3,4)=BET2(L)
         BB(4,3)=-BET2(L)
         BB(4,4)=AL4(L)
         call prod(S,BB,AA)
         call prod(AA,S,BB)
         SL(:,:,L)=BB(:,:)
      end do

      do J=1,NG
         J1=(J-1)*KONTR
         do I=1,NG
            I1=(I-1)*KONTR
            A1(:,:)=0.0d0
            A2(:,:)=0.0d0
            A3(:,:)=0.0d0
            A4(:,:)=0.0d0
            do L=LMIN,LMAX1
               do N2=1,4
                  PN2J=P(N2,L,J)
                  PPN2J=PP(N2,L,J)
                  do N1=1,4
                     SLL=SL(N1,N2,L)
                     PN1I=P(N1,L,I)*SLL
                     PPN1I=PP(N1,L,I)*SLL
                     A1(N1,N2)=A1(N1,N2)+PPN1I*PN2J
                     A2(N1,N2)=A2(N1,N2)+PN1I*PN2J
                     A3(N1,N2)=A3(N1,N2)+PPN1I*PPN2J
                     A4(N1,N2)=A4(N1,N2)+PN1I*PPN2J
                  end do
               end do
            end do
            call prod(A1,S,BB)
            call prod(S,BB,A1)
            call prod(A2,S,BB)
            call prod(S,BB,A2)
            call prod(A3,S,BB)
            call prod(S,BB,A3)
            call prod(A4,S,BB)
            call prod(S,BB,A4)
            do N2=1,KONTR
               NN2=J1+N2
               do N1=1,KONTR
                  NN1=I1+N1
                  PH1(NN1,NN2)=A1(N1,N2)
                  PH2(NN1,NN2)=A2(N1,N2)
                  PH3(NN1,NN2)=A3(N1,N2)
                  PH4(NN1,NN2)=A4(N1,N2)
               end do
            end do
         end do
      end do
   end subroutine

   subroutine rminf(M, LMAX1, ALB, AL1, X, WW, R)
      integer, intent(in) :: M, LMAX1
      integer :: I, I1, J, J1, K, K1, L, L1, M1, KPAR, ND, ND1, ND2, NITER
      real(8), intent(in) :: ALB, AL1(LMAX1), X(NG), WW(NG)
      real(8), intent(out) :: R(4,4,NG,NG,LMAX1)
      real(8) :: AA, AAM1, AAM2, SEMALB, XI, XJ, XIJ, A(NG), CC(NG,NG), &
         Q(NPH), DELT, DCOEFF, WI, WJ, WIJ, RRIJ, DRIJ

      M1=M+1
      KPAR=0
      AAM1=0.0d0
      AAM2=0.0d0
      do I=M1,LMAX1
         AAM1=dmax1(AAM1,AL1(I))
         AAM2=dmax1(AAM2,AL1(I)/dfloat(2*I-1))
      end do
      if (AAM1*AAM2>EP) KPAR=1
      ND=NG*KONTR

      SEMALB=0.0d0
      if (M<=2.and.ALB>0.8d0) SEMALB=ALB*0.8d0
      RR(:,:)=0.0d0

      do I=1,NG
         XI=1.0d0/X(I)
         A(I)=XI
         XX(I)=XI*0.5d0
         do J=1,I
            XJ=A(J)
            XIJ=XI*XJ*0.25d0
            XXX(I,J)=XIJ
            XXX(J,I)=XIJ
            XIJ=1.0d0/(XI+XJ)
            CC(I,J)=XIJ
            CC(J,I)=XIJ
         end do
      end do
      ND1=ND-KONTR+1
      ND2=ND1+1
      do I=1,KONTR
         AA=1.0d0
         IF(I==3) AA=-1.0d0
         do J=1,NG
            Q(I+(J-1)*KONTR)=AA
         end do
      end do
      do I=1,NG
         I1=KONTR*(I-1)
         DO J=1,KONTR
            II(I1+J)=I
         end do
      end do

      do I=1,ND
         I1=II(I)
         DO J=1,I
            J1=II(J)
            RR(I,J)=ALB*XXX(I1,J1)*PH1(I,J)*CC(I1,J1)
         end do
      end do

      call symm(ND, Q)

      NITER=0
      DELT=0.0d0
      if (KPAR/=0) then
         do
            call ffr(ND, ALB)
            NITER=NITER+1
            if (NITER>1000) then
               write(6, *) 'WARNING: NUMBER OF ITERATIONS EXCEEDS 1000'
               stop
            end if

            DELT=0.0d0
            DCOEFF=0.0d0
            if (NITER>=4) DCOEFF=SEMALB
            do I=1,ND
               I1=II(I)
               WI=1.0d0/WW(I1)
               do J=1,I
                  J1=II(J)
                  RRIJ=F3(I,J)*CC(I1,J1)
                  DRIJ=RRIJ-RR(I,J)
                  DELT=dmax1(DELT,dabs(DRIJ)*WI/WW(J1))
                  RR(I,J)=RRIJ+DRIJ*DCOEFF
               end do
            end do
            call symm(ND,Q)
            if (DELT<=EP) exit
            if (M==0.AND.dabs(ALB-1.0d0)<=1.0d-12) call renorm2(X, WW)
         end do
      end if

      do J=1,NG
         J1=(J-1)*KONTR
         WJ=WW(J)
         do I=1,NG
            I1=(I-1)*KONTR
            WIJ=1.0d0/(WJ*WW(I))
            do L=1,KONTR
               L1=L+J1
               do K=1,KONTR
                  K1=K+I1
                  R(K,L,I,J,M1)=RR(K1,L1)*WIJ
               end do
            end do
         end do
      end do
   end subroutine

   subroutine renorm1(M, LMAX1, AL1, X, W, WW)
      integer, intent(in) :: M, LMAX1
      integer :: I, I1, J, J1, K, K1
      real(8), intent(in) :: AL1(LMAX1), X(NG), W(NG), WW(NG)
      real(8) :: AM, DS, U, PR(NG,NG), PT(NG,NG), AL(NG), BET(NG), &
         WI, WJ, S, EPS

      AM=AL1(M+1)*2.0d0/dfloat(2*M+1)
      do I=1,NG
         DS=1.0d0
         if (M/=0) then
            U=X(I)
            U=dsqrt(1.0d0-U*U)
            do J=1,M
               DS=DS*U
            end do
         end if
         AL(I)=W(I)*DS
         BET(I)=AM*DS
         WI=1.0d0/WW(I)
         I1=(I-1)*KONTR+1
         do J=1,NG
            J1=(J-1)*KONTR+1
            WJ=WI/WW(J)
            PR(I,J)=PH1(I1,J1)*WJ
            PT(I,J)=PH2(I1,J1)*WJ
         end do
      end do

      do K=1,NG
         S=0.0d0
         do I=1,NG
            S=S+AL(I)*(PR(I,K)+PT(I,K))
         end do
         EPS=(BET(K)-S)/(PT(K,K)*AL(K)) + 1.0d0
         K1=(K-1)*KONTR+1
         PH2(K1,K1)=PH2(K1,K1)*EPS
         PH3(K1,K1)=PH3(K1,K1)*EPS
      end do
   end subroutine

   subroutine renorm2(X, WW)
      integer :: I, I1, I2, J, J1, J2, MRI
      real(8), intent(in) :: X(NG), WW(NG)
      real(8) :: R1(NG,NG), R2(NG,NG), T1(NG), T2(NG), W(NG), &
         DELT, WWJ, WWIJ, U, T1I, T2I, WX

      DELT=EP*0.1d0
      do J=1,NG
         J1=(J-1)*KONTR+1
         WWJ=WW(J)
         W(J)=WWJ*WWJ
         WWJ=1.0d0/WWJ
         do I=1,NG
            WWIJ=WWJ/WW(I)
            I1=(I-1)*KONTR+1
            I2=I1+1
            R1(I,J)=RR(I1,J1)*WWIJ
            R2(I,J)=RR(I2,J1)*WWIJ
         end do
      end do

      MRI=0

      do
         MRI=MRI+1
         if (MRI>200) then
            write(6,*) 'WARNING: RENORMALIZATION DIVERGED'
            stop
         end if

         U=0.0d0
         do I=1,NG
            T1I=1.0d0
            T2I=0.0d0
            do J=1,NG
               WX=W(J)*X(J)*2.0d0
               T1I=T1I-WX*R1(I,J)
               T2I=T2I-WX*R2(I,J)
            end do
            T1(I)=T1I
            T2(I)=T2I
            U=dmax1(U,dabs(T1I),dabs(T2I))
         end do

         if (U<=DELT) exit

         do I=1,NG
            T1I=T1(I)
            T2I=T2(I)*0.5d0
            do J=1,NG
               R1(I,J)=R1(I,J)+0.5d0*(T1I+T1(J))
               R2(I,J)=R2(I,J)+T2I
            end do
         end do
      end do

      do J=1,NG
         J1=(J-1)*KONTR+1
         J2=J1+1
         WWJ=WW(J)
         do I=1,NG
            I1=(I-1)*KONTR+1
            I2=I1+1
            WWIJ=WW(I)*WWJ
            RR(I1,J1)=R1(I,J)*WWIJ
            RR(I2,J1)=R2(I,J)*WWIJ
            RR(I1,J2)=R2(J,I)*WWIJ
         end do
      end do
   end subroutine

   subroutine dd(X,LMAX1,M,D1,D2,D3)
      integer, intent(in) :: M, LMAX1
      integer :: L, L1, L2, LMAX, M1
      real(8), intent(in) :: X
      real(8), intent(out) :: D1(LMAX1), D2(LMAX1), D3(LMAX1)
      real(8) :: DX, DDX, DD1, A, B, C, D, E, F, X1, DM, DM2, &
         DL, DL1, DL2, DL3, DDL

      LMAX=LMAX1-1
      DX=1.0d0-X*X
      DDX=dsqrt(DX)

      if (M==0) then
         D1(1)=1.0d0
         D1(2)=X
         D1(3)=(3.0d0*X*X-1.0d0)*0.5d0
         D2(1)=0.0d0
         D2(2)=0.0d0
         DD1=0.25d0*DSQRT(6.0d0)*DX
         D2(3)=DD1
         D3(1)=0.0d0
         D3(2)=0.0d0
         D3(3)=DD1
         IF (LMAX<=2) return
      end if

      if (M==1) then
         D1(2)=DSQRT(0.5d0*DX)
         D1(3)=DSQRT(1.5d0*DX)*X
         D2(2)=0.0d0
         DD1=0.5d0*DDX
         D2(3)=-DD1*(1.0d0+X)
         D3(2)=0.0d0
         D3(3)=DD1*(1.0d0-X)
         if (LMAX<=2) return
      end if

      if (M>=2) then
         A=DX*0.25d0*dsqrt(6.0d0)
         X1=1.0d0+X
         B=X1*X1*0.25d0
         X1=1.0d0-X
         C=X1*X1*0.25d0
         IF (M/=2) then
            do L=3,M
               A=A*dsqrt(dfloat(2*L-1)/dfloat(2*L))*DDX
               D=dsqrt(dfloat((2*L-1)*L)/dfloat(2*(L*L-4)))*DDX
               B=B*D
               C=C*D
            end do
         end if
         M1=M+1
         D1(M)=0.0d0
         D1(M1)=A
         D2(M)=0.0d0
         D2(M1)=B
         D3(M)=0.0d0
         D3(M1)=C
         if (LMAX<=2) return
      end if

      M1=max(M+1,3)
      DM=dfloat(M*M)
      DM2=dfloat(2*M)
      do L=M1,LMAX
         DL=dfloat(L)
         L1=L+1
         L2=L-1
         DL2=dfloat(L2)
         DL1=dfloat(2*L-1)
         DL3=DL2*DL2
         DDL=DL*DL
         A=dsqrt(DL3-DM)
         B=1.0d0/dsqrt(DDL-DM)
         C=DL1*X
         D=DL*dsqrt(DL3-4.0d0)*A
         E=B/(DL2*dsqrt(DDL-4.0d0))
         F=DL*DL2*X
         D1(L1)=(C*D1(L)-A*D1(L2))*B
         D2(L1)=(DL1*(F-DM2)*D2(L)-D*D2(L2))*E
         D3(L1)=(DL1*(F+DM2)*D3(L)-D*D3(L2))*E
      end do
   end subroutine

   subroutine ffr(ND, ALB)
      integer,intent(in) :: ND
      integer :: I, I1, J, J1, K
      real(8),intent(in) :: ALB
      real(8) :: A, A1, A2, A3, B(NPH), XI, RKJ

      do I=1,ND
         I1=II(I)
         XI=XX(I1)
         do J=1,ND
            A=0.0d0
            do K=1,ND
               A=A+RR(I,K)*PH4(K,J)
            end do
            B(J)=A
         end do
         do J=1,I
            A1=0D0
            A2=0D0
            A3=0D0
            J1=II(J)
            do K=1,ND
               RKJ=RR(K,J)
               A1=A1+RR(I,K)*PH2(K,J)
               A2=A2+PH3(I,K)*RKJ
               A3=A3+B(K)*RKJ
            end do
            F3(I,J)=ALB*(XXX(I1,J1)*PH1(I,J)+XX(J1)*A1+XI*A2+A3)
         end do
      end do
   end subroutine

   subroutine prod(A, B, C)
      integer :: I, J, K
      real(8), intent(in) :: A(4,4), B(4, 4)
      real(8), intent(out) :: C(4,4)
      real(8) :: D
      do I=1,4
         do J=1,4
            D=0.0d0
            do K=1,4
               D=D+A(I,K)*B(K,J)
            end do
            C(I,J)=D
         end do
      end do
   end subroutine

   subroutine symm(ND, Q)
      integer, intent(in) :: ND
      integer :: I, I1, J
      real(8), intent(in) :: Q(NPH)
      real(8) :: QI

      do I=2,ND
         I1=I-1
         QI=Q(I)
         do J=1,I1
            RR(J,I)=RR(I,J)*QI*Q(J)
         end do
      end do
   end subroutine

   subroutine gauss(N, IND1, Z, W)
      integer, intent(in) :: N, IND1
      integer :: IND, I, J, K, M, NITER
      real(8), intent(out) :: Z(N), W(N)
      real(8) :: A, B, C, F, X, CHECK, PA, PB, PC, DJ

      A = 1.0d0
      B = 2.0d0
      C = 3.0d0
      IND = mod(N, 2)
      K = N/2 + IND
      F = dfloat(N)

      do I=1,K
         M=N+1-I
         if (I==1) X=A-B/((F+A)*F)
         if (I==2) X=(Z(N)-A)*4.0d0+Z(N)
         if (I==3) X=(Z(N-1)-Z(N))*1.6d0+Z(N-1)
         if (I>3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
         if (I==K.and.IND==1) X=0.0d0
         NITER=0
         CHECK=1.0d-16
         do
            PB=1.0d0
            NITER=NITER+1
            if (NITER>=100) CHECK=CHECK*10.0d0
            PC=X
            DJ=A
            do J=2,N
               DJ=DJ+A
               PA=PB
               PB=PC
               PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
            end do
            PA=A/((PB-X*PC)*F)
            PB=PA*PC*(A-X*X)
            X=X-PB
            if (dabs(PB)<=CHECK*dabs(X)) exit
         end do
         Z(M)=X
         W(M)=PA*PA*(A-X*X)
         if (IND1==0) W(M)=B*W(M)
         if (.not.(I==K.and.IND==1)) then
            Z(I)=-Z(M)
            W(I)=W(M)
         end if
      end do

      if (IND1==1) then
         do I=1,N
            Z(I)=(A+Z(I))/B
         end do
      end if
   end subroutine

   function f(X, N) result(Y)
      integer, intent(in) :: N
      integer :: I
      real(8), intent(in) :: X
      real(8) :: Y, A, B, C, DI

      A=1.0d0
      B=X
      do I=1,N
         C=B
         DI=dfloat(I)
         B=((2.0d0*DI+1.0d0)*X*B-DI*A)/(DI+1.0d0)
         A=C
      end do
      Y=(A-B)/(1.0d0-X)
   end function

   function zeroin(AX,BX,TOL,N) result(B)
      integer, intent(in) :: N
      real(8), intent(in) :: AX, BX, TOL
      real(8) :: A, B, C, D, E, EPS, TOL1, FA, FB, FC, XM, &
         P, Q, R, S

      EPS=1.0d0
      do
         EPS=0.5d0*EPS
         TOL1=1.0d0+EPS
         if (TOL1<=1.0d0) exit
      end do

      A=AX
      B=BX
      FA=f(A,N)
      FB=f(B,N)
      C=A
      FC=FA
      D=B-A
      E=D

      do
         if (dabs(FC)<dabs(FB)) then
            A=B
            B=C
            C=A
            FA=FB
            FB=FC
            FC=FA
         end if
         TOL1=2.0d0*EPS*dabs(B)+0.5d0*TOL
         XM=0.5d0*(C-B)
         if (dabs(XM)<=TOL1) return
         if (FB==0.0d0) return

         if ((dabs(E)<TOL1).or.(dabs(FA)<=dabs(FB))) then
            D=XM
         else
            if (A/=C) then
               Q=FA/FC
               R=FB/FC
               S=FB/FA
               P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
               Q=(Q-1D0)*(R-1D0)*(S-1D0)
            else
               S=FB/FA
               P=2D0*XM*S
               Q=1D0-S
            end if
            IF (P>0.0d0) Q=-Q
            P=dabs(P)
            if ((2.0d0*P)>=(3.0d0*XM*Q-dabs(TOL1*Q)).or. &
               (P>=dabs(0.5d0*E*Q))) then
               D=XM
            else
               D=P/Q
            end if
         end if
         E=D
         A=B
         FA=FB
         if (dabs(D)>TOL1) B=B+D
         if (dabs(D)<=TOL1) B=B+dsign(TOL1,XM)
         FB=F(B,N)
         if ((FB*(FC/dabs(FC)))>0.0d0) then
            C=A
            FC=FA
            D=B-A
            E=D
         end if
      end do
   end function

   subroutine markov(NG, IND1, X, W)
      integer, intent(in) :: NG, IND1
      integer :: I, K, M, N, N1
      real(8), intent(out) :: X(NG), W(NG)
      real(8) :: BBB, TOL, E, DN, DD, A, B, C, D, &
         YB, Y, X1, DM, WW

      N=NG-1
      BBB=1.0d-2
      IF (N>=30) BBB=1.0d-3
      IF (N>=100) BBB=1.5d-4
      IF (N>=200) BBB=3.0d-5
      TOL=1.0d-15
      E=2.0d0
      DN=dfloat(N)
      DD=2.0d0*DN+1.0d0
      A=-1.0d0
      K=0
      DD=DD*DD
      YB=F(A,N)
      outer: do
         Y=YB
         B=A
         do
            B=B+BBB
            if (1.0d0-B<=1.0d-15) exit outer
            YB=F(B,N)
            if (Y/YB<0.0d0) exit
         end do
         K=K+1
         X1=zeroin(A,B,TOL,N)
         if (K>=2) BBB=(X1-X(K-1))*0.2d0
         X(K)=X1
         N1=N-1
         D=f(X1,N1)
         C=D*2.0d0*(DN+1.0d0)*DN
         WW=4.0d0*(1.0d0+X1)*DD/(C*C)
         W(K)=WW
         E=E-WW
         A=B
      end do outer

      K=K+1
      W(K)=E
      X(K)=1.0d0
      M=2*(K-1)
      A=0.0d0
      DM=dfloat(M)
      do I=1,K
         A=A+W(I)*dabs(X(I))**DM
      end do
      B=0.5d0*dfloat(M+1)
      A=A*B
      if (IND1==1) then
         DO I=1,K
            W(I)=0.5d0*W(I)
            X(I)=0.5d0*(X(I)+1.0d0)
         end do
      end if
      N=K
   end subroutine
end module
