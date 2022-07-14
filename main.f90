module global
    integer,parameter::nx=24,ny=24,nz=6,nt=100,nxy=nx*ny
    real(8)::drx,dry,dt
    real(8)::dcx,dcy,dcz
end module global

program KF_welding
    use global
    implicit none
    integer::i,j,k,t
    integer::info
    integer,parameter::lwork=100*nxy
    integer,dimension(nxy)::ipiv
    integer,dimension(4)::iseed
    real(8),dimension(lwork)::work
    real(8)::tt,x,y,r
    real(8)::r11,r12,r21,r22,rx,ry,rt
    real(8)::ktr,rhocr,a,b,c,kt,rhoc,Tref,time,cmh,T0,sT,sq,sy
    real(8),dimension(nxy)::vy,vq,vyp
    real(8),dimension(2*nxy)::vxp,vxm
    real(8),dimension(nxy,nxy)::mR,mInv
    real(8),dimension(nxy,2*nxy)::mH
    real(8),dimension(2*nxy,nxy)::mK,mPHT
    real(8),dimension(2*nxy,2*nxy)::mF1,mQ,mPp,mPm,maux
    real(8),allocatable,dimension(:,:,:)::hT

!!!!!!!!!!!!!!!!
! Problem Data !
!!!!!!!!!!!!!!!!

    a=0.12d0
    b=0.12d0
    c=0.003d0
    time=2.d0
    T0=300.d0

!!!!!!!!!!!!!!!!!!!!!
! Measurement Noise !
!!!!!!!!!!!!!!!!!!!!!

    sy=1.d-2

!!!!!!!!!!!!!!!!!
! Random Number !
!   Generator   !
!!!!!!!!!!!!!!!!!

    iseed(1)=0
    iseed(2)=1500
    iseed(3)=3000
    iseed(4)=4095

!!!!!!!!!!!!!
! Grid size !
!!!!!!!!!!!!!

    drx=a/real(nx,8)
    dry=b/real(ny,8)
    dt=time/real(nt,8)

    dcx=drx
    dcy=dry
    dcz=c/real(nz,8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constant Thermal Properties !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Tref=600.d0
    ktr=kt(Tref)
    rhocr=rhoc(Tref)

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Synthetic Measurements !
! & Reference Heat Flux  !
!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(hT(nx,ny,nz))
    hT=T0
    open(unit=10,file="vq.dat",status="replace")
    open(unit=11,file="vy_raw.dat",status="replace")
    open(unit=12,file="vy.dat",status="replace")
    do t=1,nt
        tt=real(t,8)*dt
        vq=0.d0
        call vqf(vq,t)
        call modelocompleto(hT,vq)
        vy=reshape(hT(:,:,1),(/nxy/))
        write(unit=10,fmt=99)'vq','"Heat Flux [W/m2]"',nx,ny,t,tt
        write(unit=12,fmt=99)'vy','"Temperature [K]"',nx,ny,t,tt
        do i=1,nx
            x=(real(i,8)-0.5d0)*drx
            do j=1,ny
                k=i+(j-1)*nx
                y=(real(j,8)-0.5d0)*dry
                call dlarnv(3,iseed,1,r)
                vy(k)=vy(k)+r*sy
                write(unit=10,fmt=*)x,y,vq(k)
                write(unit=12,fmt=*)x,y,vy(k)
            enddo
        enddo
        write(unit=11,fmt=*)vy
    enddo
    close(unit=10)
    close(unit=11)
    close(unit=12)
    deallocate(hT)

!!!!!!!!!!!!!!!!!!!!
! Evolution Matrix !
!!!!!!!!!!!!!!!!!!!!

    r11=1.d0-ktr*dt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    r12=1.d0-ktr*dt*(drx**(-2.d0)+2.d0*dry**(-2.d0))/rhocr
    r21=1.d0-ktr*dt*(2.d0*drx**(-2.d0)+dry**(-2.d0))/rhocr
    r22=1.d0-2.d0*ktr*dt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    rx=ktr*dt/rhocr/drx**2.d0
    ry=ktr*dt/rhocr/dry**2.d0
    rt=dt/c/rhocr
    mF1=0.d0
    mF1(1,1)=r11
    do i=2,nx-1
        mF1(i,i)=r21
    enddo
    mF1(nx,nx)=r11
    do i=nx+1,nx*(ny-1),nx
        mF1(i,i)=r12
        do j=2,nx-1
            mF1(i+j-1,i+j-1)=r22
        enddo
        mF1(nx+i-1,nx+i-1)=r12
    enddo
    mF1(nx*(ny-1)+1,nx*(ny-1)+1)=r11
    do i=nx*(ny-1)+2,nx*ny-1
        mF1(i,i)=r21
    enddo
    mF1(nx*ny,nx*ny)=r11
    do i=1,nx
        mF1((i-1)*ny+1,(i-1)*ny+2)=rx
        do j=2,ny-1
            mF1((i-1)*ny+j,(i-1)*ny+j-1)=rx
            mF1((i-1)*ny+j,(i-1)*ny+j+1)=rx
        enddo
        mF1(i*ny,i*ny-1)=rx
    enddo
    do i=1,nx
        mF1(i,i+ny)=ry
    enddo
    do i=nx+1,(nx-1)*ny
        mF1(i,i-ny)=ry
        mF1(i,i+ny)=ry
    enddo
    do i=(nx-1)*ny+1,nx*ny
        mF1(i,i-ny)=ry
    enddo
    do i=1,nx*ny
        mF1(i,nx*ny+i)=rt
        mF1(nx*ny+i,nx*ny+i)=1.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!
! Observation Model !
!!!!!!!!!!!!!!!!!!!!!
    
    mH=0.d0
    cmh=-c/(6.d0*ktr)
    do i=1,nxy
        mH(i,i)=1.d0
        mH(i,nxy+i)=cmh
    enddo

!!!!!!!!!!!!!!!!!!!!
! Noise Covariance !
!!!!!!!!!!!!!!!!!!!!

    sT=1.d-1
    sq=1.d2
    mQ=0.d0
    mR=0.d0
    do i=1,nxy
        mQ(i,i)=sT**2.d0
        mQ(nxy+i,nxy+i)=sq**2.d0
        mR(i,i)=sy**2.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!
! Initializing Remaining!
! Vectors and Matrices  !
!!!!!!!!!!!!!!!!!!!!!!!!!

    mK=0.d0
    mPp=mQ
    mPm=mQ
    mPHT=0.d0
    mInv=0.d0
    vxp=0.d0
    vxm=0.d0
    do i=1,nxy
        vxp(i)=T0
        vxp(nxy+i)=0.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! State Estimation Problem!
! Classical Kalman Filter !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    open(unit=10,file="vTest.dat",status="replace")
    open(unit=11,file="vqest.dat",status="replace")
    open(unit=12,file="vy_raw.dat",status="old")
    do t=1,nt
        tt=real(t,8)*dt
    !!!!!!!!!!!!!
    ! Read Data !
    !!!!!!!!!!!!!

        read(unit=12,fmt=*)vy

    !!!!!!!!!!
    ! x = Fx !
    !!!!!!!!!!

        call dgemv('N',2*nxy,2*nxy,1.d0,mF1,2*nxy,vxp,1,0.d0,vxm,1)
        vxp=vxm

    !!!!!!!!!!!!!!!!!
    ! P = FPF^T + Q !
    !!!!!!!!!!!!!!!!!

        mPm=mPp
        call dgemm('N','N',2*nxy,2*nxy,2*nxy,1.d0,mF1,2*nxy,mPm,2*nxy,0.d0,mPp,2*nxy)
        mPm=mQ
        call dgemm('N','T',2*nxy,2*nxy,2*nxy,1.d0,mPp,2*nxy,mF1,2*nxy,1.d0,mPm,2*nxy)

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! K = PH^T(HPH^T + R)^-1 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!

        mPHT=0.d0
        call dgemm('N','T',2*nxy,nxy,2*nxy,1.d0,mPm,2*nxy,mH,nxy,0.d0,mPHT,2*nxy)
        mInv=mR
        call dgemm('N','N',nxy,2*nxy,nxy,1.d0,mH,nxy,mPHT,2*nxy,1.d0,mInv,nxy)
        call dgetrf(nxy,nxy,mInv,nxy,ipiv,info)
        call dgetri(nxy,mInv,nxy,ipiv,work,lwork,info)
        mK=0.d0
        call dgemm('N','N',2*nxy,2*nxy,nxy,1.d0,mPHT,2*nxy,mInv,nxy,0.d0,mK,2*nxy)

    !!!!!!!!!!!!!!!!!!!
    ! x = x + K(y-Hx) !
    !!!!!!!!!!!!!!!!!!!

        call dgemv('N',nxy,2*nxy,1.d0,mH,nxy,vxm,1,0.d0,vyp,1)
        vxp=vxm
        call dgemv('N',2*nxy,nxy,1.d0,mK,2*nxy,vy-vyp,1,1.d0,vxp,1)

    !!!!!!!!!!!!!!!
    ! P = (I-KH)P !
    !!!!!!!!!!!!!!!
        
        call dgemm('N','N',2*nxy,2*nxy,nxy,1.d0,mK,2*nxy,mH,nxy,0.d0,maux,2*nxy)
        do i=1,2*nxy
            maux(i,i)=maux(i,i)-1.d0
        enddo
        call dgemm('N','N',2*nxy,2*nxy,2*nxy,1.d0,-maux,2*nxy,mPm,2*nxy,0.d0,mPp,2*nxy)

    !!!!!!!!!!
    ! Output !
    !!!!!!!!!!
        
        write(unit=10,fmt=99)'vT','"Temperature [K]"',nx,ny,t,tt
        write(unit=11,fmt=99)'vq','"Heat Flux [W/m2]"',nx,ny,t,tt
        do i=1,nx
            x=(real(i,8)-0.5d0)*drx
            do j=1,ny
                y=(real(j,8)-0.5d0)*dry
                k=i+(j-1)*nx
                write(unit=10,fmt=*)x,y,vxp(k)
                write(unit=11,fmt=*)x,y,vxp(nxy+k)
            enddo
        enddo
        write(*,*)maxval(vxp(1:nxy)),maxval(vxp(nxy+1:2*nxy))
        write(*,*)maxval(vy),maxval(vyp)
        write(*,*)
    enddo
    close(unit=10)
    close(unit=11)
    close(unit=12)

!!!!!!!!!!!!!!!!!!
!Formatting Rules!
!!!!!!!!!!!!!!!!!!

    99 format ('TITLE = ',A,/,&
        'Variables = "x [m]", "y [m]",',A,/,&
        'ZONE, i=',i3,', j=',i3,', f=point, STRANDID=',i3,',SOLUTIONTIME=',es14.6)
end program

subroutine vqf(vq,it)
    use global
    implicit none
    integer::i,j,k,it,i0,i1,j0,j1
    real(8)::tt,x,y
    real(8),dimension(nx*ny)::vq
    tt=real(it,8)*dt
    i0=8
    i1=10
    j0=8
    j1=10
    do i=1,nx
        x=(real(i,8)-0.5d0)*drx
        do j=1,ny
            k=i+(j-1)*nx
            y=(real(j,8)-0.5d0)*dry
            if((i.ge.i0).and.(i.le.i1).and.(j.ge.j0).and.(j.le.j1))then
                vq(k)=1.d5
            else
                vq(k)=0.d0
            endif
        enddo
    enddo
end subroutine