!! Module to define tha parameters
!! ---------------------------------------------------------
!!
module mod_force
  use mod_param
  use mtmod
  implicit none
  real*8::rij,xij,yij,zij,dr0,tt,avg_dx
  real*8::fLJ
  save
  !
contains
!! ---------------SUBROUINES ---------------
!!


subroutine initial_condition
  implicit none
  integer*4::ii,Npsm,np0,jj,kk,nrs0,nstheta,nsphi
  integer*4::ip
  real*8::rrnd,trnd,phrnd,rtheta
  if(nrun==1) then
! Initialize the position of the particle in the chain (boundary of the tissue)
  dth = (2.0*pi)/dble(Nrs)
  is = 0
  do ii = 1,nint(R1-R0);  ! Matrix 
    nrs0 = nint(4.0*pi*(R0+ii-1)**2*rhos)
    nstheta = ceiling(pi*(R0+ii-1)*rhos)
    dth = (pi)/dble(nstheta)
    do jj = 1,nstheta+1;  ! Matrix 
      theta = dble(jj-1)*dth+(ii-1)*dth*0.5
      rtheta = (R0+ii-1)*sin(theta)
      nsphi = ceiling(2*pi*(rtheta))
      dphi = (2*pi)/dble(nsphi)
      do kk = 1,nsphi;  ! Matrix 
        phi = dble(kk)*dphi+(jj-1)*dphi*0.5
        is = is + 1
        xrs(is) = dble(R0+ii-1)*dsin(theta)*dcos(phi)
        yrs(is) = dble(R0+ii-1)*dsin(theta)*dsin(phi)
        zrs(is) = dble(R0+ii-1)*dcos(theta)
      enddo
    enddo
  enddo
  Nrs = is;
!!- Random packing
!  do is = 1,Nrs;  ! Matrix 
!62  rrnd = (R1-R0)*grnd()+R0;  trnd = (2.0*pi)*grnd()     ! For compact packing
!    xrs(is) = rrnd*cos(trnd);
!    yrs(is) = rrnd*sin(trnd);
!    dr0 = 2.0*a;
!    do ii = 1,is-1;
!      yij = yrs(ii)-yrs(is); 
!      xij = xrs(ii)-xrs(is);
!      rij = dsqrt(xij**2+yij**2);
!      if(rij<dr0) dr0 = rij
!    enddo
!    if(dr0<a) goto 62
!  enddo; 
  xrs0 = xrs;
  yrs0 = yrs;
  zrs0 = zrs;
  call write_track
! Populate the rest of PSM 
  ! When whole spheroid is populated
   Npsm = int((4.0/3.0)*pi*(R0-a)**3/a*rho0) 
  do np = 1,Npsm;
63  rrnd = (R0-a)*grnd();  trnd = pi*grnd();  phrnd = 2.0*pi*grnd();     ! For compact packing
!63  rrnd = (R0-4.0*a)+2.0d0*grnd();  trnd = (2.0*pi)*grnd()
    xx(np) = rrnd*dsin(trnd)*dcos(phrnd);
    yy(np) = rrnd*dsin(trnd)*dsin(phrnd);
    zz(np) = rrnd*dcos(trnd);
    dr0 = 2.0*a;
    do ii = 1,np-1;
      zij = zz(np)-zz(ii);
      yij = yy(np)-yy(ii); 
      xij = xx(np)-xx(ii);
      rij = dsqrt(xij**2+yij**2+zij**2);
      if(rij<dr0) dr0 = rij
    enddo
    if(dr0<a) goto 63
  print*,np
  enddo 
  np = Npsm;
  dtime(1:np) = 0.0d0;!-tauG;
  call write_track
! The motility of cells decreases as the cells move away from tailbud
   mot = mot0; ttp = 0.0d0;
   if (FGF_grad) then
   do ip = 1,np;
     ttp(ip) = (R0-xx(ip))/Df
     mot(ip) = mot0*dexp(-ttp(ip)/tau)
   enddo
   endif
   else
! Read from old trajectory file.
     call read_traj
   endif
end subroutine initial_condition


!!-
! Force acting on the cell aggregate 
subroutine cal_force_matrix
  implicit none
  integer::ip,jp
    do is = 1,Nrs;
! Soft repulsion
      do jp = 1,Nrs;
  ! To avoid self counting
        if(jp==is) then
        else
          zij = zrs(jp)-zrs(is); 
          yij = yrs(jp)-yrs(is); 
          xij = xrs(jp)-xrs(is);
          rij = dsqrt(xij**2+yij**2+zij**2);
          rij = max(rij,rmin);
       !   if(rij<a) then
       !     fxrs(is) = fxrs(is) + kss*(rij-a)*xij/rij;
       !     fyrs(is) = fyrs(is) + kss*(rij-a)*yij/rij;
       !   endif
          !fLJ = -4.d0*epsm*((rn*Ar/(rij)**(rn+2)-an*Aa0/(rij)**(an+2)))     !! LJ 12-6
          fLJ = 0.0d0
          if(rij<rlj) then
          fLJ = -(2.0d0*alphar*epsm/rij**4)*((rlj/rij)**2-1)*(3.0d0*(a*rlj/rij)**2-2.0d0*rlj*rlj-a*a)
          endif
! Matrix-matrix LJ interaction
          fxrs(is) = fxrs(is) + fLJ*xij;
          fyrs(is) = fyrs(is) + fLJ*yij;
          fzrs(is) = fzrs(is) + fLJ*zij;
        endif  
      enddo
! Link rearrangement
       zij = zrs0(is)-zrs(is); 
       yij = yrs0(is)-yrs(is); 
       xij = xrs0(is)-xrs(is);
       rij = dsqrt(xij**2+yij**2+zij**2);
       if(rij>ls_cut) then
         xrs0(is) = xrs(is);
         yrs0(is) = yrs(is);
         zrs0(is) = zrs(is);
       endif
! Elastic Contribution
       fxrs(is) = fxrs(is) + Es*xij;
       fyrs(is) = fyrs(is) + Es*yij;
       fzrs(is) = fzrs(is) + Es*zij;
! Inter species repulsion
      do jp = 1,np;
        zij = zz(jp)-zrs(is); 
        yij = yy(jp)-yrs(is); 
        xij = xx(jp)-xrs(is);
        rij = dsqrt(xij**2+yij**2+zij**2);
        if(rij<a) then
          fxrs(is) = fxrs(is) + kas*(rij-a)*xij/rij;
          fyrs(is) = fyrs(is) + kas*(rij-a)*yij/rij;
          fzrs(is) = fzrs(is) + kas*(rij-a)*zij/rij;
        endif
      enddo
    enddo
end subroutine cal_force_matrix
!!

! Force acting on the cell aggregate 
subroutine cal_force_cell
  implicit none
  integer::ip,jp
      rho = 0.0d0
      do ip = 1,np;
        Aa(ip) = Aa0! + 2.0d0*Aa0*(1-dexp(dtime(ip)))
      enddo
    do ip = 1,np;
      ! Only repulsion from matrix
      do jp = 1,Nrs;
          zij = zrs(jp) - zz(ip);
          yij = yrs(jp) - yy(ip);
          xij = xrs(jp) - xx(ip);
          rij = dsqrt(xij**2+yij**2+zij**2);
        if(rij<a) then
          fx(ip) = fx(ip) + kas*(rij-a)*xij/rij;
          fy(ip) = fy(ip) + kas*(rij-a)*yij/rij;
          fz(ip) = fz(ip) + kas*(rij-a)*zij/rij;
        endif
    !      fLJ = 4.d0*epst*((r12*Ar/(rij)**(r12+2)-a6*Aamin/(rij)**(a6+2))) 
    !      fx(ip) = fx(ip) + fLJ*xij/rij;
    !      fy(ip) = fy(ip) + fLJ*yij/rij;
      enddo

! Calculcate density for PSM, which gives us the attraction potential
   if(density_Aa) then
      do jp = 1,np;
! To avoid self counting
        if(jp==ip) then
        else
          zij = zz(ip)-zz(jp); 
          yij = yy(ip)-yy(jp); 
          xij = xx(ip)-xx(jp);
          rij = dsqrt(xij**2+yij**2+zij**2);
          if(rij<=rc) rho(ip) = rho(ip) + 1.d0
        endif
      enddo
      rho(ip) = rho(ip)/area_c
      Aa(ip) = Aamin + (Aa0-Aamin)/(1.d0+dexp((rho(ip)-rhoc)/l_rho))
   endif
!
      do jp = 1,np;
! To avoid self counting
        if(jp==ip) then
        else
          zij = zz(jp)-zz(ip); 
          yij = yy(jp)-yy(ip); 
          xij = xx(jp)-xx(ip);
          rij = dsqrt(xij**2+yij**2+zij**2);
          rij = max(rij,rmin);
          fLJ = 0.0d0
          if(rij<rlj) then
          fLJ = -(2.0d0*alphar*eps0/rij**4)*((rlj/rij)**2-1)*(3.0d0*(a*rlj/rij)**2-2.0d0*rlj*rlj-a*a)
          endif
          !fLJ = -4.d0*eps0*((rn*Ar/(rij)**(rn+2)-an*0.5d0*(Aa(ip)+Aa(jp))/(rij)**(an+2))) !! LJ 12-6
          !fLJ = 4.d0*eps0*((12.d0*Ar/(rij)**14-6.d0*Aa(ip)/(rij)**8)) 
          !fLJ = 4.d0*eps0*((4.d0*Ar/(rij)**6-2.d0*Aa(ip)/(rij)**4)) 
! Cell-cell interaction
          fx(ip) = fx(ip) + fLJ*xij;
          fy(ip) = fy(ip) + fLJ*yij;
          fz(ip) = fz(ip) + fLJ*zij;
        endif  
      enddo
    enddo
end subroutine cal_force_cell
!!

!! Evolve the matrix 
subroutine evolve_matrix
  implicit none
  integer::ip
    do ip = 1,Nrs;
      xrs(ip) = xrs(ip) + delta*(fxrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      yrs(ip) = yrs(ip) + delta*(fyrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      zrs(ip) = zrs(ip) + delta*(fzrs(ip))/mu ! This term can be evolved using second or third order schemes as well
    enddo
end subroutine evolve_matrix
!!

!! Evolve the cell: Deterministic process
subroutine evolve_cell
  implicit none
  integer::ip
   ! Deterministic process
    do ip = 1,np;
      xx(ip) = xx(ip) + delta*(fx(ip))/gm ! This term can be evolved using second or third order schemes as well
      yy(ip) = yy(ip) + delta*(fy(ip))/gm ! This term can be evolved using second or third order schemes as well
      zz(ip) = zz(ip) + delta*(fz(ip))/gm ! This term can be evolved using second or third order schemes as well
    enddo
!!
end subroutine evolve_cell
!!

!! Evolve the cell: Wiener process
subroutine evolve_cell_Weiner
  implicit none
  integer::ip
  real*8::rnumx,rnumy,rnumz
   ! Wiener process
    do ip = 1,np;
if (FGF_grad) mot(ip) = mot0*dexp(-(ttp(ip))/tau)
      ttp(ip) = ttp(ip) + delta
64    rnumx = gasdev(); rnumy = gasdev(); rnumz = gasdev();
      !if(max(abs(rnumx),abs(rnumy))>4.0d0) 
      if(max(abs(rnumx),abs(rnumy),abs(rnumz))>8.0d0) goto 64
      xix(ip) = dsqrt(2.d0*mot(ip)*gm)*rnumx
      xiy(ip) = dsqrt(2.d0*mot(ip)*gm)*rnumy
      xiz(ip) = dsqrt(2.d0*mot(ip)*gm)*rnumz
      xx(ip) = xx(ip) + dsqrt(delta)*(xix(ip))/gm ! Noise term; evolved at the end with a sqrt(delta) term
      yy(ip) = yy(ip) + dsqrt(delta)*(xiy(ip))/gm ! Noise term; evolved at the end with a sqrt(delta) term
      zz(ip) = zz(ip) + dsqrt(delta)*(xiz(ip))/gm ! Noise term; evolved at the end with a sqrt(delta) term
    enddo
end subroutine evolve_cell_Weiner
!!

!!-- Cell proliferation
subroutine cell_proliferation
  implicit none
  integer::ip,ii,np0
  real*8::rhotemp,Unp,dPr
  real*8::dstep
    !rnum = grnd() !White noise for position of particle in y direction
! Check the density if cell generation is needed.
! Cell generation; add one cell, if there is enough space. 
        np0 = np
    do ip = 1,np;
      if((tt-dtime(ip))>=tauG) then
        ii = floor((np0-1)*grnd())+1 !White noise for position of particle in y direction
        theta = pi*grnd()
        phi = (2.0*pi)*grnd()
        xx(np+1) = (xx(ii))+a*dsin(theta)*dcos(phi);
        yy(np+1) = (yy(ii))+a*dsin(theta)*dsin(phi);
        zz(np+1) = (zz(ii))+a*dcos(theta);
        
        np = np+1;
        ttp(np) = 0.d0;
        dr0 = 2.0*a;
        Unp = 0.0d0;
        do ii = 1,np-1;
          zij = zz(np)-zz(ii); 
          yij = yy(np)-yy(ii); 
          xij = xx(np)-xx(ii);
          rij = dsqrt(xij**2+yij**2+zij**2);
          if(rij<dr0) dr0 = rij
          if(rij<rlj) then
            !Unp = Unp + 4.d0*eps0*((Ar/(rij)**(rn)-0.5d0*(Aa(ip)+Aa(ii))/(rij)**(an)))
            Unp = Unp + alphar*eps0*((a/rij)**2-1)*((rlj/rij)**2-1)**2
          endif
        enddo
!  Loop for the wall
        do ii = 1,nrs;
          zij = zz(np)-zrs(ii); 
          yij = yy(np)-yrs(ii); 
          xij = xx(np)-xrs(ii);
          rij = dsqrt(xij**2+yij**2+zij**2);
          if(rij<dr0) dr0 = rij
          if(rij<a) then
            Unp = Unp + 0.5d0*kas*(rij-a)**2;
          endif
        enddo
        dPr = dexp(-Unp/mot0)
        dstep = grnd()
        
        !if(dr0<0.9*a) then
        if(dstep>dPr) then
          np = np-1;  
        else
          print*,dstep,dPr,dr0,np,it
          dtime(ip) = tt
          dtime(np) = tt
        endif
        rhotemp = dble(np)/(pi*R0**2)
      endif
    enddo !Cell proliferation
end subroutine cell_proliferation

!!-- Cell proliferation uniform
subroutine cell_proliferation_uniform
  implicit none
  integer::ip,ii,np0
! Cell generation; add one cell, if it is old enough. 
        np0 = np
    do ip = 1,np;
       ii = floor((np0-1)*grnd())+1 !White noise for picking the random cell
       print*,ii,dtime(ii)
       if((tt-dtime(ip))>=tauG) then
          theta = pi*grnd()
          phi = (2.0*pi)*grnd()
          xx(np+1) = (xx(ii))+a*dsin(theta)*dcos(phi);
          yy(np+1) = (yy(ii))+a*dsin(theta)*dsin(phi);
          zz(np+1) = (zz(ii))+a*dcos(theta);
          np = np+1;
          ttp(np) = 0.d0;
          dtime(ii) = tt
          dtime(np) = tt
       endif
       if(np.gt.np0) goto34
    enddo !Cell proliferation fix
34  print*,np0,np
end subroutine cell_proliferation_uniform

!!-- Cell proliferation uniform at the center
subroutine cell_proliferation_uniform_center
  implicit none
! Cell generation at the center; add one cell, if it is old enough. 
          theta = pi*grnd()
          phi = (2.0*pi)*grnd()
          xx(np+1) = a*dsin(theta)*dcos(phi);
          yy(np+1) = a*dsin(theta)*dsin(phi);
          zz(np+1) = a*dcos(theta);
          np = np+1;
!!
end subroutine cell_proliferation_uniform_center


!!
end module mod_force
