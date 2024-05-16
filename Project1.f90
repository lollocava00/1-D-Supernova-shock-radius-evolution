!*********************************************************************
!  Program for project 1- Calcolo 20121. Lorenzo Cavazzini
!  Solves the hydrostatic equilibrium equation in a NFW+BCG potential
!  Then solves the Fe diffusion eq. (for Perseus)
!*********************************************************************
PROGRAM project1
parameter(jmax=5000)
IMPLICIT REAL*8 (a-h,o-z)
REAL*8, DIMENSION(jmax):: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),rho(jmax),& !PART 1 VARIABLES
            rho2(jmax),rho3(jmax),mhern(jmax),rhonfw(jmax),rhost(jmax),mdark(jmax), &
            grvnfw(jmax),grvnfw2(jmax),lnd(jmax),lnd2(jmax),lnd3(jmax),rho_an(jmax),&
            mgas(jmax),mgas2(jmax),mgas3(jmax),T(jmax)

REAL*8:: msol,mu,mp,rmin,rmax,mvir,rvir,mbcg,ahern,b,Tmed, &
         alphast,alphasn,zfesn,kappa,lturb,vturb,dt,tend,tnow,t1,t2, &
         n,m,fb,fb2,fb3, rho0,rho02,rho03

REAL*8, DIMENSION(jmax):: zfe(jmax),rhofedot(jmax),rhofe(jmax),zfest(jmax),& !PART 2 VARIABLES
       amfeiniz(jmax),amfe(jmax),gradzfe(jmax),zfeobs(jmax), amfeobs(jmax),&
       rhofeobs(jmax), rhofeiniz(jmax)

1100 format(1(1pe12.4))
1200 format(2(1pe12.4))
1300 format(3(1pe12.4))
1400 format(4(1pe12.4))
1500 format(5(1pe12.4))
1600 format(6(1pe12.4))

!  constants

pi=3.14159
msol = 1.989e33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
years=3.156d7
zfesol=1.8e-3
zfesn=0.744/1.4

!! set the grid(uniform)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rmin = 0.*cmkpc
rmax = 3000.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
end do
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
end do
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
! open(10,file='grid.dat',status='unknown')
! do j=1,jmax
!    write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
! enddo
! close(10)

vol(1)=4.1888*r(1)**3  !=0
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
end do

!! parametri del problema

rho0nfw=7.35d-26
rs=435.7*cmkpc
!********!!! I CHOSE THIS TO GET A BARION FRACTION OF 0.16 !!!********
rho0=4.46d-26
rho02=8.76d-26
rho03=1.68d-25    
!*********************************************************************
ticm=8.9e7
rvir=2797.*cmkpc
fc=1.138799      
mvir=1.3e15*msol
mbcg=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))   

!analytical density profile of DM and star density profile
do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
   rhost(j)=mbcg/(2.*pi)*(ahern/rr(j))/(ahern+rr(j))**3
end do

!mass profile
open(20,file='masse.dat')
mnfw(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)   !numerical mass profile of DM
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc  !analytical mass profile
   mhern(j)=mbcg*r(j)**2/(r(j)+ahern)**2  !mass profile of stars (Hernquist profile)
   write(20,1400)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
end do

close(20)
print*, "mass profile data has been saved in: masse.dat"

!deriving temperature profile and saving it to file
open(20,file='temp.dat',status='unknown')
do j=1,jmax
   x=rr(j)/(rvir*0.5)
   T(j)=ticm*(1.35*(((x/0.045)**1.9+0.45)/(((x/0.045)**1.9)+1.))*(1./(1.+(x/0.6)**2)**0.45))
   write(20,1200) rr(j)/cmkpc,T(j)
end do
close(20)
print*, "non constant-temperature profile has been saved in: temp.dat"

!!!! numerical integration of ODE solved for density profile !!!!!!!!!!!!!!!

!grvnfw(j) è la forza gravitazionale al punto j
grvnfw(1)=0.         !! ok per alone NFW, isotermo o beta-model
grvnfw2(1)=0.
do j=2,jmax
   grvnfw(j)=guniv*mnfw(j)/r(j)**2
   grvnfw2(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
end do

lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
lnd2(1)=log(rho02)
lnd3(1)=log(rho03)
DO j=2,jmax
   gg=grvnfw(j)
   gg2=grvnfw2(j)
   Tmed=0.5*(T(j)+T(j-1))
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) 
   lnd2(j)=lnd2(j-1)-gg2*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm)
   lnd3(j)=lnd3(j-1)-gg2*(mu*mp)*(rr(j)-rr(j-1))/(boltz*Tmed)-(log(T(j))-log(T(j-1)))
END DO

do j=1,jmax
   rho(j)=exp(lnd(j))       !isothermal
   rho2(j)=exp(lnd2(j))     !with BCG
   rho3(j)=exp(lnd3(j))     !with BCG +dT/dr
end do
print*, "ODE for hidrostatic equilibrium integrated"

!computes the analytical solution for rho
b=(8.*pi*guniv*mu*mp*rho0nfw*rs**2)/(27.*boltz*ticm)
rho_an(1)=rho0
do j=2,jmax
   !b=(8.*pi*guniv*mu*mp*rho0nfw*rs**2)/(27.*boltz*T(j))
   rho_an(j)=rho_an(1)*exp((-27.*b/2.)*(1.-log(1.+rr(j)/rs)/(rr(j)/rs)))
end do

!saving density profile for DM and gas
open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1600) rr(j)/cmkpc,rho(j),rho2(j),rho3(j),rhonfw(j),rho_an(j)
end do
print*, "density data has been saved in: density.dat"

!in order to obtain the barion fraction we first calculate the gas mass
mgas(1)=0  
mgas2(1)=0 
mgas3(1)=0
do j=2,jmax
mgas(j)=mgas(j-1)+rho(j-1)*vol(j)
mgas2(j)=mgas2(j-1)+rho2(j-1)*vol(j)
mgas3(j)=mgas3(j-1)+rho3(j-1)*vol(j)
end do
!barion fraction calculatated at the virial radius
fb=mgas(jmax)/(mnfw(jmax)+mgas(jmax))
fb2=(mhern(jmax)+mgas2(jmax))/(mhern(jmax)+mnfw(jmax)+mgas2(jmax))
fb3=(mhern(jmax)+mgas3(jmax))/(mhern(jmax)+mnfw(jmax)+mgas3(jmax))

print*, ''
print*, ' baryion fraction for case1' ,real(fb)
print*, ' baryion fraction for case2(+BCG)' ,real(fb2)
print*, ' baryion fraction for case3(+BCG and dT/dr)' ,real(fb3)
print*, ''

!***********************************************************************
!! At this point we have the gas density profile and we can proceed
!! with the integration of the diffusion equation for rhofe
!***********************************************************************

!! Set the initial Fe abundance profile

zfeout=0.4*zfesol   !! this is the background abundance !!

do j=1,jmax
   x=rr(j)/(80.*cmkpc)

   zfeobs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  !Perseus!
   zfeobs(j)=zfeobs(j) - zfeout   !! subtract z_Fe,out which is the background metallicity
   zfeobs(j)=max(zfeobs(j),0.)    !!eliminates negative values

   zfest(j)=1.*zfesol             !!set the stellar abundance
   zfe(j)=0. 
   rhofe(j)=rho3(j)*zfe(j)/1.4
   rhofeobs(j)=rho3(j)*zfeobs(j)/1.4
end do

!this is the (cumulative)MASS abundance of Fe
amfeiniz(1)=rhofeobs(1)*vol(1)
rhofeiniz(1)=rhofeobs(1)
amfeobs(1)=rhofeobs(1)*vol(1)
amfe(1)=rhofe(1)*vol(1)
do j=2,jmax
   rhofeiniz(j)=rhofeobs(j)
   amfeiniz(j)=amfeiniz(j-1)+rhofeobs(j-1)*vol(j)
   amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)    !for diffusion only
   amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)             !for source only   
end do

open(20,file='zfe_initial.dat')
do j=1,jmax
   write(20,1400)rr(j)/cmkpc,zfeobs(j)/zfesol,rhofeobs(j),amfeiniz(j)/msol
end do
close(20)

!! boundary conditions
zfe(1)=zfe(2)
zfe(jmax)=zfe(jmax-1)
rhofe(1)=rho3(1)*zfe(1)/1.4
rhofe(jmax)=rho3(jmax)*zfe(jmax)/1.4

!!  set the diffusion coefficient kappa = C*v*l (for now constant)
 vturb=260.e5   
 lturb=15.*cmkpc  
 rscala=30.*cmkpc
!**************************************************
 kappa=0.11*vturb*lturb
 !kappa=2.4e29
!**************************************************
 print*, 'diffusion coefficient:', real(kappa)
 print*, ''
 !write(*,1200) amfeobs(180)/msol, amfeobs(jmax)/msol
 !write(*,1200) amfe(180)/msol, amfe(jmax)/msol

tnow=13.7*1.e9*years
!*********CHANGE THE TIME OF INTEGRATION**************
time0=tnow-5.*1.e9*years
t1=tnow-4.*1.e9*years   !1Gyr of integrated time
t2=tnow-3.*1.e9*years   !2Gyr of integrated time
!*****************************************************
time=time0
tend=tnow
n=0                     !counter for the number of cycles
dt=0.4*((r(5)-r(4))**2/(2.*kappa)) !!HAS TO SATISFY STABILITY CONDITION!!
m=0                     !counter for diffusion timescale calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Here start the time integration (use FTCS method)
DO WHILE(time<tend)     
   n=n+1
   time=time+dt

!********CHANGE THE SNu TO GET DIFFERENT SOURCE TERM*************
   !snu=0.7*(time/tnow)**(-1.1)
   snu=0.15
   !snu=0.5
!****************************************************************
!! write the source terms (SNIa and stellar winds)
   aml=7.5
   !alphast=4.7e-20*(time/tnow)**(-1.26)
   !alphasn=4.436e-20*(snu/aml)*(time/tnow)**(-slope)
   alphast= 4.7e-20
   alphasn= 5.91e-21*snu
!THIS IS THE SOURCE TERM
   do j=2,jmax-1
      rhofedot(j)=(alphast*zfest(j)/1.4+alphasn*zfesn)*rhost(j)
   end do

!! UPDATE DENSITY FOR SOURCE TERM
   do j=2,jmax-1
      rhofe(j)=rhofe(j) + dt*rhofedot(j)
      zfe(j)=rhofe(j)/rho3(j)*1.4
   end do

 !! set the boundary conditions 
   zfe(1)=zfe(2)
   zfe(jmax)=zfe(jmax-1)
   rhofe(1)=rhofe(2)
   rhofe(jmax)=rhofe(jmax-1)

!! set the boundary conditions (diffusion only)
   ! zfeobs(1)=zfeobs(2)
   ! zfeobs(jmax)=zfeobs(jmax-1)
   ! rhofeobs(1)=rhofeobs(2)
   ! rhofeobs(jmax)=rhofeobs(jmax-1)

!  diffusive step   !  check the Fe conservation !
   do j=2,jmax-1
      gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
   enddo
   gradzfe(1)=0.        !! B.C. !!
   gradzfe(jmax)=0.

   do j=2,jmax-1
      rhojp1=0.5*(rho3(j+1)+rho3(j))  !! rho centered at "j+1" !!
      rhoj=0.5*(rho3(j-1)+rho3(j))    !! rho centered at "j" !!
      rhofe(j)=rhofe(j) &
               + (dt/1.4)*(r(j+1)**2*kappa*rhojp1*gradzfe(j+1) &
               -r(j)**2*kappa*rhoj*gradzfe(j))   &
               / (0.33333333*(r(j+1)**3-r(j)**3))
      zfe(j)=1.4*rhofe(j)/rho3(j)  !! update Z_Fe with the new rho_Fe !!
   enddo

!! ONLY DIFFUSION STEP
   ! do j=2,jmax-1
   !    gradzfe(j)=(zfeobs(j)-zfeobs(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
   ! end do
   ! gradzfe(1)=0.        !! B.C. !!
   ! gradzfe(jmax)=0.

   ! do j=2,jmax-1
   !    rhojp1=0.5*(rho3(j+1)+rho3(j))  !! rho centered at "j+1" !!
   !    rhoj=0.5*(rho3(j-1)+rho3(j))    !! rho centered at "j" !!
   !    rhofeobs(j)=rhofeobs(j) &
   !             + (dt/1.4)*(r(j+1)**2*kappa*rhojp1*gradzfe(j+1) &
   !             -r(j)**2*kappa*rhoj*gradzfe(j))   &
   !             / (0.33333333*(r(j+1)**3-r(j)**3))
   !    zfeobs(j)=1.4*rhofeobs(j)/rho3(j)  !! update Z_Fe with the new rho_Fe !!
   ! end do

   !cumulative mass of iron (DIFFUSION TERM ONLY)
   ! amfeobs(1)=rhofeobs(1)*vol(1)
   ! do j=2,jmax
   !    amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)
   ! end do

   !cumulative mass of iron (SOURCE TERM ONLY)
   ! amfe(1)=rhofe(1)*vol(1)
   ! do j=2,jmax
   !    amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
   ! end do

   !theoretical mass profile of iron (source only)
   !amfe_theo=5.e-22*mbcg*(time-8.7e9*years)


!! data from 1Gyr of integration!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (t1-dt<time .and. time<t1) THEN  !Fe abundance, density and mass after 1 Gyr
      open(26,file='diff1.dat',status='unknown')
      do j=1,jmax
         write(26,1400) zfeobs(j)/zfesol,rhofeobs(j),zfe(j)/zfesol,rhofe(j)
      end do
      close(26)
      !write(*,1200) amfeobs(180)/msol, amfeobs(jmax)/msol              !for diffusion
      !write(*,1300) amfe(180)/msol, amfe(jmax)/msol, amfe_theo/msol      !for source and source+diffusion
      print*, '1 GYR of integrated time has passed'
      print*, 'iron abundance and density saved in: diff1.dat'
      print*, ''
   END IF

!! data from 2Gyr of integration!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (t2-dt<time .and. time<t2) THEN  !Fe abundance, density and mass after 2 Gyr
      open(27,file='diff2.dat',status='unknown')
      do j=1,jmax
         write(27,1400) zfeobs(j)/zfesol,rhofeobs(j), zfe(j)/zfesol,rhofe(j)
      end do
      close(27)
      !write(*,1200) amfeobs(180)/msol, amfeobs(jmax)/msol              !for diffusion
      !write(*,1300) amfe(180)/msol, amfe(jmax)/msol, amfe_theo/msol     !for source and source+diffusion
      print*, '2GYRs of integrated time have passed'
      print*, 'iron abundance and density saved in: diff2.dat'
      print*, ''
   END IF

! cycle for the diffusion time scale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! IF (amfeobs(135)<(0.5*amfeiniz(135)) .and. m==0) THEN  !checks the mass at 80kpc
   !    m=1
   !    print*, 'diffusion time-scale(GYR):', (time-time0)/3.156e16 
   !    print*, 'calculated at Fe peak 80kpc'
   !    print*, ''
   ! END IF


END DO 
!*******************************************************************************

!write(*,1200) amfeobs(180)/msol, amfeobs(jmax)/msol
!write(*,1300) amfe(180)/msol, amfe(jmax)/msol, amfe_theo/msol
print*, 'time integration completed'
print*,'TIME (Gyr) = ',real(time/3.156e16), 'Number of cycles: ', n

!writes metallicity data to file
open(21,file='diff.dat',status='unknown')
do j=1,jmax
   write(21,1500)rr(j)/cmkpc,zfeobs(j)/zfesol,rhofeobs(j), zfe(j)/zfesol,rhofe(j)
enddo
close(21)
print*, 'The final results have been saved in: diff.dat'

STOP
END PROGRAM project1
