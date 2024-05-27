!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    program written by Lorenzo Cavazzini.											!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DATA
integer :: N, i, j, k, cool_onoff
real*8 :: pi,cmpc,cmkpc,yr,kbol,mu,mp
parameter (N=500)
parameter (k=500) 
parameter(pi=3.141592)
parameter(cmpc=3.085d18)
parameter(cmkpc=1000.*cmpc)
parameter(yr=3.156d7)
parameter(kbol=1.38d-16)
parameter(mu=0.61)
parameter(mp=1.67d-24)
END MODULE DATA

PROGRAM ZEUS
USE DATA
IMPLICIT NONE
real*8 :: xa(N), xb(N), xmax, xmin, deltax, dxa(N), dxb(N)
real*8 :: d(N), e(N), v(N), P(N), s(N), Temp(N) !DENSITA', ENERGIAINTERNA, VELOCITA', PRESSIONE, MOMENTO
real*8 :: q(N) !VISCOSITA' ARTIFICIALE
real*8 :: g2a(N), g2b(N), g31a(N), g31b(N), dvl1a(N), dvl1b(N) 
real*8 :: F1(N), F2(N), F3(N), M(N),  dstar(N),  e_dstar(N), vstar(N)
real*8 :: divV(N)
real*8 :: rshock(k), tsh(k), rsedlaw(k), Lx(k), Etot(k), Ek(k), Et(k)
real*8 :: Ekin,Eth,rho0,T0,E0
real*8 :: dtmin, tmax, t, C2, gam, Cv, t1, t2, t3, tt, LumX, cfl, cont
real*8 :: t4, t5, t6, t7, t8
integer :: sdr, ncicli, g
real*8, EXTERNAL :: Cool


!CREAZIONE DOPPIA GRIGLIA (xa e xb)

xmin=0.
xmax=80.*cmpc

!GRIGLIA "a"
xa(1)=xmin+(xmax-xmin)*(i-2.)/(N-1.)
do i=1,N
	xa(i)= xmin+(xmax-xmin)*(i-2.)/(N-1.)
end do

deltax=xa(3)-xa(2)

!GRIGLIA "b"
do i=1, N-1
	xb(i)=0.5*(xa(i)+xa(i+1))
end do
xb(N)=xb(N-1)+(xb(N-1)-xb(N-2))   !! add the last calculated Delta_xb to xb(N-1)

do i=2, N-1
	dxa(i)=xa(i+1)-xa(i)
	dxb(i)=xb(i)-xb(i-1)
end do

dxa(1)=xa(2)-xa(1)
dxa(N)=dxa(N-1)
dxb(1)=dxb(2)
dxb(N)=xb(N)-xb(N-1)


!DEFINIZIONE FATTORI DI SCALA METRICI 

sdr=1    !! this parameter selects the type of coordinates: 0 = Cartesian

if (sdr==0) then  !! Cartesian !!

	do i=1, N
	g2a(i)=1.
	g2b(i)=1.
	g31a(i)=1.
	g31b(i)=1.
 
	end do
	do i=1, N-1
	dvl1a(i)=xa(i+1)-xa(i)   !! Note that is centered in xb(i)
	end do
	dvl1a(N)=dvl1a(N-1)
	do i=2, N
	dvl1b(i)=xb(i)-xb(i-1)  !! Note that it is centered in xa(i)
	end do
	dvl1b(1)=dvl1b(2)


	
else if (sdr==1) then   !! spherical !!
	do i=1, N
	g2a(i)=xa(i)
	g31a(i)=xa(i)
	g2b(i)=xb(i)
	g31b(i)=xb(i)
	end do

	do i=1, N-1
	dvl1a(i)=(xa(i+1)**3-xa(i)**3)/3.

	end do
	dvl1a(N)=dvl1a(N-1)
	do i=2, N
	dvl1b(i)=(xb(i)**3-xb(i-1)**3)/3.
	end do
	dvl1b(1)=dvl1b(2)

end if

open(20,file='grid.dat')
do i=1,N
   write(20,999)xa(i),xb(i),dxa(i),dxb(i),dvl1a(i),dvl1b(i)
enddo
close(20)
999 format(6(1pe12.4))
print*, 'staggered grid created. beggining time integration...'

!IMPLEMENTAZIONE CONDIZIONI INIZIALI

!***********************************
rho0=2.d-24
T0=1.d4
E0=1.d51
!***********************************

gam=5./3.
cv=1.99d8    
t=0.
!***** INTEGRATION TIME ************
t1=2.d4*yr
t2=4.d4*yr
t3=6.d4*yr
t4=8.d4*yr
t5=1.d5*yr
t6=2.d5*yr
t7=3.d5*yr
t8=4.d5*yr
tmax=5.d5*yr
!***********************************
tt=1.d3*yr
!***********************************
c2=3.
cfl=0.01

do i=1, N
    d(i)=rho0
	Temp(i)=T0    
	v(i)=0.
	e(i)=cv*d(i)*Temp(i)
	p(i)=(gam-1.)*e(i)
end do	

!SN ENERGY INJECTION
e(2)=E0/((4./3.)*pi*xa(4)**3)
e(3)=e(2)


p(2)=(gam-1.)*e(2)
p(3)=p(2)
Temp(2)=e(2)/(cv*d(2))
Temp(3)=Temp(2)

CALL BCb(e)
CALL BCb(p)
CALL BCb(Temp)


        ncicli=0
		cont=0
		g=0
!***************************************************************************
do while (t<tmax)      !!!! HERE STARTS THE TIME INTEGRATION !!!!!
        ncicli=ncicli+1

!CALCOLO DTMIN

        dtmin=1.d30   !! any very large value !!
        p=(gam-1.)*e
	do i=2, N-1
		 dtmin=min(dtmin,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gam*p(i)/d(i))))
	end do
        dtmin=cfl*dtmin
        t=t+dtmin
		cfl=min(0.5,cfl*1.1)
        !print*,'ncicli, dtmin = ',ncicli, real(dtmin),real(t)


!SOURCE STEP
!SUBSTEP I: AGGIORNAMENTO DELLA VELOCITÀ PER GRADIENTE DI P

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i))	
	end do
	CALL BCa(v)


!CALCOLO Q
	do i=2, N-1
		if ((v(i+1)-v(i))<0.) then
			q(i)=C2*d(i)*(v(i+1)-v(i))**2
		else 
			q(i)=0.
		end if
	end do
	CALL BCb(q)

!SUBSTEP II: AGGIORNAMENTO PER VISCOSITÀ ARTIFICIALE

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(q(i)-q(i-1))/((d(i)+d(i-1))*dxb(i))
	end do
	CALL BCa(v)

	do i=2, N-1
		e(i)=e(i)-dtmin*q(i)*(v(i+1)-v(i))/dxa(i)
	end do
	CALL BCb(e)

!SUBSTEP III: AGGIORNAMENTO PER RISCALDAMENTO DA COMPRESSIONE
	do i=2,N-1
		divV(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
	end do
	CALL BCa(divV)

	do i=2, N-1
		e(i)=e(i)*(1.-0.5*dtmin*(gam-1.)*divV(i))/(1.+0.5*dtmin*(gam-1.)*divV(i))
	end do
	CALL BCb(e)

	do i=2,N-1
		Temp(i)=e(i)/(cv*d(i))
	end do
	CALL BCb(Temp)

	! AGGIORNAMENTO PER IL RADIATIVE COOLING
!********************************************************************
	cool_onoff=1	!this parameter switches on/off radiative cooling
!********************************************************************
	IF(cool_onoff==1) THEN
		do i=2, N-1
			e(i)=e(i)-Cool(Temp(i),d(i))*dtmin
			Temp(i)=e(i)/(cv*d(i))
			if(Temp(i)<1.d4) then
				Temp(i)=1.d4
				e(i)=cv*Temp(i)*d(i)
			end if
		end do
		!print*, e(100)/d(100), Temp(100), Cool(Temp(100),d(100))
		CALL BCb(e)
		CALL BCb(Temp)
	END IF


!!!!!!TRANSPORT STEP (use Upwind first order only)

	do i=2, N-1       !! here define the momentum density
		s(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
	end do	

	CALL BCa(s)

!AGGIORNAMENTO DENSITÀ

	do i=2, N-1       !! here select the value of the density at the interface "i"
		if (v(i)>0.) then
			dstar(i)=d(i-1)     !! at i !!
		else
			dstar(i)=d(i)
		end if
	end do
	dstar(N)=dstar(N-1)
	dstar(1)=dstar(3)

	do i=2, N
		F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)    !! at i !!	
	end do

!AGGIORNAMENTO ENERGIA

	do i=2, N-1
		M(i)=dstar(i)*v(i)
	end do
	CALL BCa(M)
	
	
	do i=2, N-1
		if (v(i)>0.) then
			e_dstar(i)=e(i-1)/d(i-1)   !! at i !!
		else
			e_dstar(i)=e(i)/d(i)
		end if
	end do
	e_dstar(N)=e_dstar(N-1)
	e_dstar(1)=e_dstar(3)


	!ORA AGGIORNO LA DENSITÀ
	do i=2, N-1
		d(i)=d(i)-dtmin*(F1(i+1)-F1(i))/dvl1a(i)
	end do 
	CALL BCb(d)
	

	do i=2, N
		F2(i)=e_dstar(i)*M(i)*g2a(i)*g31a(i)				
	end do
	CALL BCa(F2)

	do i=2, N-1
		e(i)=e(i)-dtmin*(F2(i+1)-F2(i))/dvl1a(i)
	end do

	CALL BCb(e)


!AGGIORNAMENTO MOMENTO 

	do i=2, N-1
		if ((v(i-1)+v(i))*0.5>0) then
			vstar(i)=v(i-1)       !! at i-1/2  !!
		else
			vstar(i)=v(i)
		end if
	end do

	CALL BCb (vstar)

	do i=1, N-1
		F3(i)=vstar(i+1)*0.5*(M(i)+M(i+1))*g2b(i)*g31b(i)   !! questo e' a i+1/2, occhio !!  
	end do
	
	do i=2, N-1
		s(i)=s(i)-dtmin/dvl1b(i)*(F3(i)-F3(i-1))
	end do

	CALL BCa(s)

	do i=2, N-1
		v(i)=2.*s(i)/(d(i)+d(i-1))
	end do

	CALL BCa(v)


! WRITE FILES AT INTERMEDIATE INTEGRATION TIME
	j=0
	if(t>=t1-dtmin .and. t<t1) then
		j=1
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 2*10^4 yr'
	else if (t>=t2-dtmin .and. t<t2) then
		j=2
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 4*10^4 yr'
	else if (t>=t3-dtmin .and. t<t3) then
		j=3
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 6*10^4 yr'
	else if (t>=t4-dtmin .and. t<t4) then
		j=4
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 8*10^4 yr'
	else if (t>=t5-dtmin .and. t<t5) then
		j=5
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 1*10^5 yr'
	else if (t>=t6-dtmin .and. t<t6) then
		j=6
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 2*10^5 yr'
	else if (t>=t7-dtmin .and. t<t7) then
		j=7
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 3*10^5 yr'
	else if (t>=t8-dtmin .and. t<t8) then
		j=8
		CALL data_strip(d,v,e,p,Temp)
		print*, 'integration time: 4*10^5 yr'
	end if

    
! WRITE R_SHOCK, LUMINOSITY IN X-RAY AND ENERGY EVERY T=10^3 YEARS 
	LumX=0.
	Ekin=0.
	Eth=0.
	if(t>=cont+dtmin .and. t<cont+2*dtmin) then		!ROOT OF ALL EVILS
		g=g+1
		do i=2,N
			if(d(i)==maxval(d)) then
				rshock(g)=xb(i)
			end if
			if(Temp(i)>1.d6) then
				LumX=LumX+4.*pi*dvl1b(i)*Cool(Temp(i),d(i))
			end if
			Ekin=Ekin+0.5*4*pi*dvl1a(i)*d(i)*(v(i))**2
			Eth=Eth+e(i)*4*pi*dvl1a(i)-cv*rho0*T0*4*pi*dvl1b(i)
			!Eth=Eth+(4./3.)*pi*(xa(i)**3-xa(i-1)**3)*e(i-1)
		end do

		tsh(g)=cont
		rsedlaw(g)=(2.*E0/rho0)**(1./5.)*t**(2./5.)
		Lx(g)=LumX
	
		Ek(g)=Ekin
		Et(g)=Eth
		Etot(g)=Ekin+Eth

	 	cont=cont+tt
	end if

enddo       !! here the time integration ends !!
!***************************************************************************************
print*, 'time integration completed'

IF(cool_onoff==1) THEN
	open(20,file='results_cooling.dat')
	do i=1,N  !! write the results in the file "results.dat"
		if (v(i)< 1.d-30) then
			v(i)=0.
		end if
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i)/(mp*mu),v(i)/1.d5,e(i)/d(i),p(i),Temp(i)
	end do
	1000 format(7(1pe12.4))
	close(20)

	open(30,file='r_shock_cooling.dat')
	do g=1,k
		write(30,1001)  tsh(g)/yr, rshock(g)/cmpc, rsedlaw(g)/cmpc,Lx(g), Ek(g)/E0,Et(g)/E0,Etot(g)/E0
	end do
	1001 format(7(1pe12.4))
	close(30)
	print*, '-- radiative cooling is ON --'
ELSE
	open(21,file='results.dat')
	do i=1,N  !! write the results in the file "results.dat"
		if (v(i)< 1.d-30) then
			v(i)=0.
		end if
		write (21,1002) xa(i)/cmpc,xb(i)/cmpc,d(i)/(mp*mu),v(i)/1.d5,e(i)/d(i),p(i),Temp(i)
	end do
	1002 format(7(1pe12.4))
	close(21)

	open(31,file='r_shock.dat')
	do g=1,k
		write(31,1003)  tsh(g)/yr, rshock(g)/cmpc, rsedlaw(g)/cmpc,Lx(g), Ek(g)/E0,Et(g)/E0,Etot(g)/E0
	end do
	1003 format(7(1pe12.4))
	close(31)
	print*, '-- radiative cooling is OFF --'
END IF

print*, 'final results saved to files.'

END PROGRAM ZEUS


SUBROUTINE BCa(z1) !corrette BC per velocità e momento (riflessione)
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z1

z1(2)=0.
z1(1)=-z1(3)
z1(N)=z1(N-1)
! z1(1)=z1(2)       !! ouflow !!
! z1(N)=z1(N-1)

END SUBROUTINE BCa

SUBROUTINE BCb(z2) ! BC di outflow tradizionali
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z2
z2(1)=z2(2)
z2(N)=z2(N-1)
END SUBROUTINE BCb

Real*8 FUNCTION Cool(Temp1, d1)
USE DATA
IMPLICIT NONE
Real*8, intent(IN):: Temp1, d1
real*8:: Tkev
	Tkev=Temp1/1.16d7
		if (Temp1<1.99d4) then
			Cool=((d1/2.17d-24)**2)*1.544d-22*(Tkev/0.0017235)**6
		else if (Temp1<2.32d5 .AND. Temp1>1.99d4) then
			Cool=((d1/2.17d-24)**2)*6.72d-22*(Tkev/0.02)**0.6	
		else if(Temp1>2.32d5) then 
			Cool=((d1/2.17d-24)**2)*1.d-22*(8.6d-3*Tkev**(-1.7)+0.058*Tkev**(0.5)+0.063)
		end if

END FUNCTION Cool

SUBROUTINE data_strip(d,v,e,p,Temp)
USE DATA
IMPLICIT NONE
real*8,dimension(N) :: d,v,e,p,Temp
character :: str

	if (j > 0 .and. j < 9) write (str, "(i0)") j
	! if (j==1) then 
	! 	str='1'
	! else if (j==2) then
	! 	str='2'
	! else if (j==3) then
	! 	str='3'
	! else if (j==4) then
	! 	str='4'
	! else if (j==5) then
	! 	str='5'
	! else if (j==6) then
	! 	str='6'
	! else if (j==7) then
	! 	str='7'
	! else if (j==8) then
	! 	str='8'
	! end if

	IF(cool_onoff==1) THEN
		open(j*100,file='results'//str//'_cooling.dat')
		do i=1,N
			if (v(i)< 1.d-30) then
				v(i)=0.
			end if
			write(j*100,1004) d(i)/(mp*mu),v(i)/1.d5,e(i)/d(i),p(i),Temp(i)
		end do
		1004 format(5(1pe12.4))
		close(j*100)
	ELSE 
		open(j*100,file='results'//str//'.dat')
		do i=1,N
			if (v(i)< 1.d-30) then
				v(i)=0.
			end if
			write(j*100,1005) d(i)/(mp*mu),v(i)/1.d5,e(i)/d(i),p(i),Temp(i)
		end do
		1005 format(5(1pe12.4))
		close(j*100)
	END IF
END SUBROUTINE data_strip