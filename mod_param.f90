!! Module to define tha parameters
!! ---------------------------------------------------------
!!
module mod_param
  implicit none
  save
  !
  integer*4::np,Npt,Npy,Npx,it,Nt,nlw,Npx0
  integer*4::Npr,Npr0,Npr1,is,Nrs
  real*8::theta,phi,dth,dphi
  integer*4::fno,ndata,nrun
  real*8::Df,a,kss,kas,Es,kt,mu,gm,gmm,mot0,tau,tauG,delta,rho0,rhom,rhos
  real*8::rlj,alphar                             ! New LJ potential parameter from Wang & Frenkel
  real*8::eps0,epsm,Ar,Aa0,Aamin,rc,area_c,rhoc,l_rho,Aaw,Aaw0   !LJ Potential
  real*8::rn,an,r12,a6,epst,rnw,anw               !LJ Potential
  real*8::pi,R0,R1,ls_cut,rmin
  integer::iseed
  real*8,allocatable,dimension(:)::xx,yy,zz,vv,Pr,xix,xiy,xiz,mot,dtime
  real*8,allocatable,dimension(:)::xx0,yy0,zz0,fx,fy,fz,Aa,rho
  real*8,allocatable,dimension(:)::xrs,yrs,zrs
  real*8,allocatable,dimension(:)::xrs0,yrs0,zrs0
  real*8,allocatable,dimension(:)::fxrs,fyrs,fzrs
  real*8,allocatable,dimension(:) :: ran
  real*8,allocatable,dimension(:)::ttp
  character*500 :: cnpar,formp
  logical :: FGF_grad
  logical :: density_Aa
  logical :: stress_dependent
  logical :: uniform_prol
  logical :: uniform_prol_center
  !
!  private
!  public:: open_track,write_track,close_track
contains
!! ---------------SUBROUINES ---------------
!!

!! Read input parameter file 
subroutine read_input
  implicit none
  pi = 4.d0*datan(1.d0)
  open(unit=11,file='para.in',status='old'); 
  read(11,*) nrun
  read(11,*) Npt
  read(11,*) Nt
  read(11,*) ndata
  read(11,*) Npy
  read(11,*) Npx0
  read(11,*) Npr0
  read(11,*) Npr1
  read(11,*) Npx
  read(11,*) Df
  read(11,*) iseed
  read(11,*) kss
  read(11,*) kas
  read(11,*) Es
  read(11,*) a
  read(11,*) rlj
  read(11,*) gm
  read(11,*) mu
  read(11,*) tau
  read(11,*) tauG
  read(11,*) ls_cut
  read(11,*) delta
  read(11,*) rhos
  read(11,*) rho0
  read(11,*) rhom
  read(11,*) eps0
  read(11,*) epsm
  read(11,*) epst
  read(11,*) Ar
  read(11,*) Aa0
  read(11,*) Aaw0 
  read(11,*) Aamin
  read(11,*) rhoc
  read(11,*) l_rho
  read(11,*) rn
  read(11,*) an 
  read(11,*) rnw
  read(11,*) anw
  read(11,*) r12
  read(11,*) a6 
  read(11,*) FGF_grad
  read(11,*) density_Aa
  read(11,*) stress_dependent
  read(11,*) uniform_prol
  read(11,*) uniform_prol_center
  close(11)
  rmin = 0.7*a;
  alphar = 2.0d0*(rlj/a)**2*(3.0d0/(2*((rlj/a)**2-1.0d0)))**3
end subroutine read_input
!!

!! Allocate arrays 
subroutine allocate_arrays
  implicit none
  integer::ip,fno
  real::rrr
! Allocate the arrays
  allocate(xx(Npt),yy(Npt),zz(Npt),vv(Npt),Pr(Npt),xix(Npt),xiy(Npt),xiz(Npt))
  allocate(dtime(Npt)) !Time of last division
  allocate(xx0(Npt),yy0(Npt),zz0(Npt),mot(Npt),ttp(Npt))
  allocate(xrs(Nrs),yrs(Nrs),zrs(Nrs))
  allocate(xrs0(Nrs),yrs0(Nrs),zrs0(Nrs))
  allocate(fxrs(Nrs),fyrs(Nrs),fzrs(Nrs))
  allocate(fx(Npt),fy(Npt),fz(Npt))
  allocate(Aa(Npt),rho(Npt))
  allocate(ran(3))

end subroutine allocate_arrays
!!


!! Read trajectory of particles from old file
subroutine read_traj
  implicit none
  integer::ip,fno
  real::rrr
  open(unit=13,file='xpt.in',status='old');
  open(unit=14,file='ypt.in',status='old'); 
  read(13,formp)rrr,(xx(ip),ip=1,ndata*Npy);
  read(14,formp)rrr,(yy(ip),ip=1,ndata*Npy);
  read(14,formp)rrr,(zz(ip),ip=1,ndata*Npy);
  close(13)
  close(14)
  np = nint(rrr)
  print*,np
  open(unit=15,file='mot.txt',status='old'); 
  read(15,formp)(mot(ip),ip=1,np);
  read(15,formp)(ttp(ip),ip=1,np);
  close(15)

end subroutine read_traj
!!


!! Open particle trajectory files
subroutine open_track
  implicit none
  integer::ip,fno
  fno = 101;
  open(unit=fno,file='xptrack.out',status='unknown'); fno=fno+1;
  open(unit=fno,file='yptrack.out',status='unknown'); fno=fno+1; 
  open(unit=fno,file='zptrack.out',status='unknown'); fno=fno+1; 
  open(unit=fno,file='mottrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='fxtrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='fytrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='fztrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='rxtrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='rytrack.out',status='unknown'); 
  fno=fno+1; open(unit=fno,file='rztrack.out',status='unknown'); 
end subroutine open_track
!!

!! Write particle trajectory
subroutine write_track
  implicit none
  integer::ip,fno,np1,np2
  fno = 101;
  np1 = np+Nrs;
  xx(np+1:np1) = xrs(1:Nrs);
  yy(np+1:np1) = yrs(1:Nrs);
  zz(np+1:np1) = zrs(1:Nrs);
  write(fno,formp)real(np1),real(Nrs),(xx(ip),ip=1,ndata*Npy); fno = fno+1;
  write(fno,formp)real(np1),real(Nrs),(yy(ip),ip=1,ndata*Npy); fno = fno+1; !Npt/2)
  write(fno,formp)real(np1),real(Nrs),(zz(ip),ip=1,ndata*Npy); fno = fno+1; !Npt/2)
  write(fno,formp)real(np1),real(Nrs),(mot(ip),ip=1,ndata*Npy); !Npt/2)
  fno=fno+1; write(fno,formp)real(np1),real(Nrs),(fx(ip),ip=1,ndata*Npy);
  fno=fno+1; write(fno,formp)real(np1),real(Nrs),(fy(ip),ip=1,ndata*Npy);
  fno=fno+1; write(fno,formp)real(np1),real(Nrs),(fz(ip),ip=1,ndata*Npy);
  fno=fno+1; write(fno,formp)real(np),real(Nrs),(xix(ip),ip=1,ndata*Npy);
  fno=fno+1; write(fno,formp)real(np),real(Nrs),(xiy(ip),ip=1,ndata*Npy);
  fno=fno+1; write(fno,formp)real(np),real(Nrs),(xiz(ip),ip=1,ndata*Npy);

end subroutine write_track
!!

!! Close particle trajectory files
subroutine close_track
  implicit none
  integer::ip,fno
  fno = 101;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)
  fno=fno+1;  close(fno)

end subroutine close_track
!!

!! ---------------FUNCTIONS-----------------
!! Function to calculate gaussian deviate from a uniform 
!! deviate created by Mersenne-Twister
        function gasdev()
        use mtmod
        implicit none
        double precision :: gasdev
        integer::iset
        double precision::fac,gset,rsq,v1,v2,ran1
        save iset,gset
        data iset/0/
        if(iset.eq.0)then
     1  v1=2.0d0*grnd()-1
        v2=2.0d0*grnd()-1
        rsq=v1*v1+v2*v2
        if(rsq.ge.1.or.rsq.eq.0)goto 1
        fac=sqrt(-2.0d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
        else
        gasdev=gset
        iset=0
        endif
        return
        end function gasdev 
!!
!! -----------------------------------------
!!
end module mod_param
