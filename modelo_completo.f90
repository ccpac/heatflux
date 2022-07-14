!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Solução do Modelo Completo!
!Método dos Volumes Finitos!
!	     Implicito         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modelocompleto(hTold,qest)
    use global
    implicit none
    integer::i,j,k
    integer,parameter::tci=nx,tcj=ny,tck=nz
    real(8)::tol,eps,kt,rhoc
    real(8),dimension(tci*tcj)::qest
    real(8),dimension(tci,tcj,tck)::hT,hTold,hTnew,hkt,hrhoc,ap0,ap
    real(8),dimension(tci-1,tcj,tck)::ae
    real(8),dimension(tci,tcj-1,tck)::an
    real(8),dimension(tci,tcj,tck-1)::at

!!!!!!!!!!!!
!Tolerância!
!!!!!!!!!!!!
    tol=1.d-8
!!!!!!!!!!!!!!
! Método de  !
!Gauss-Seidel!
!!!!!!!!!!!!!!
    hTnew=hTold
    eps=1.d0
    do while(eps.gt.tol)
        hT=hTnew
    !!!!!!!!!!!!!!
    !Propriedades!
    !Termofísicas!
    !!!!!!!!!!!!!!
        do i=1,tci
            do j=1,tcj
                do k=1,tck
                    hkt(i,j,k)=kt(hTnew(i,j,k))
                    hrhoc(i,j,k)=rhoc(hTnew(i,j,k))
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Coeficientes!
    !!!!!!!!!!!!!!
        do i=1,tci-1
            do j=1,tcj
                do k=1,tck
                    ae(i,j,k)=2.d0*hkt(i,j,k)*hkt(i+1,j,k)/(hkt(i,j,k)+hkt(i+1,j,k))*dcy*dcz/dcx
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj-1
                do k=1,tck
                    an(i,j,k)=2.d0*hkt(i,j,k)*hkt(i,j+1,k)/(hkt(i,j,k)+hkt(i,j+1,k))*dcx*dcz/dcy
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj
                do k=1,tck-1
                    at(i,j,k)=2.d0*hkt(i,j,k)*hkt(i,j,k+1)/(hkt(i,j,k)+hkt(i,j,k+1))*dcx*dcy/dcz
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj
                do k=1,tck
                  ap0(i,j,k)=hrhoc(i,j,k)*dcx*dcy*dcz/dt
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!!!!!!
    !Pontos Interiores!
    !!!!!!!!!!!!!!!!!!!
        do i=2,tci-1
            do j=2,tcj-1
                do k=2,tck-1
                    ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                    hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !!!!!!!!!!!!!!
        i=1
        do j=2,tcj-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                            an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                            at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !!!!!!!!!!!!!!
        i=tci
        do j=2,tcj-1
            do k=2,tck-1
                ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                            an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                            at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !!!!!!!!!!!!!!
        j=1
        do i=2,tci-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !!!!!!!!!!!!!!
        j=tcj
        do i=2,tci-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno z=0!
    !!!!!!!!!!!!!!
        k=1
        do i=2,tci-1
            do j=2,tcj-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno z=c!
    !!!!!!!!!!!!!!
        k=tck
        do i=2,tci-1
            do j=2,tcj-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !	      y=0!
    !!!!!!!!!!!!!!
        i=1
        j=1
        do k=2,tck-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !!!!!!!!!!!!!!	
        i=1
        j=tcj
        do k=2,tck-1
            ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        k=1
        do j=2,tcj-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        k=tck
        do j=2,tcj-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        do k=2,tck-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        do k=2,tck-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        k=1
        do j=2,tcj-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        k=tck
        do j=2,tcj-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        j=1
        k=1
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        j=1
        k=tck
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        j=tcj
        k=1
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        j=tcj
        k=tck
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        j=1
        k=1
        ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        j=1
        k=tck
        ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        j=tcj
        k=1
        ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        j=tcj
        k=tck
        ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        k=1
        ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        k=tck
        ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        k=1
        ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        k=tck
        ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!!!!!!!!!!
    !Teste da Converg�ncia!
    !!!!!!!!!!!!!!!!!!!!!!!
        eps=sum((hTnew-hT)**2.d0)
    enddo
    hTold=hTnew
end subroutine

!!!!!!!!!!!!!!
!Propriedades!
!Termof�sicas!
!!!!!!!!!!!!!!
real(8) function kt(T)
    real(8) T
    kt=12.45d0+.014d0*T+2.517d-6*(T**2.d0)
end function

real(8) function rhoc(T)
    real(8) T
    rhoc=1324.75d0*T+3557900.d0
end function