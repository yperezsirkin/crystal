!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calcula el volumen de un polimero centrado en  !!!
!!!                 0, 0, 0                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crea_prot

use system
use protein
use MPI
use molecules, only : vsol
implicit none
!!!variables internas
integer ix, iy, iz
integer lx, ly, lz
integer sub
real *8 r2
real *8 subdim, subvol, subnumb,cellvol
real *8 r_prot, r_prot2, r_protdelta2
real *8 x, y, z
real *8 sumprotpol
real *8 d_prot, d_prot2
!!!
call MPI_Barrier(  MPI_COMM_WORLD, ierr)
r_prot = dble(protR * delta)
r_prot2 = r_prot * r_prot
r_protdelta2 = (r_prot+delta)*(r_prot+delta)
sub = 100
d_prot = 2*r_prot
d_prot2 = d_prot*d_prot 
subdim = dble(delta)/dble(sub)
subvol = subdim*subdim*subdim
!subnumb = sub*sub*sub
subnumb =  dble(sub*sub*sub)*2**3

cellvol = dble(delta)*dble(delta)*dble(delta)
protvol = 0.d0
sumprotvol = 0.d0

!!!!

do ix = -protR, protR
do iy = -protR, protR
do iz = -protR, protR
  

 do lx = -sub, sub - 1
 do ly = -sub, sub - 1
 do lz = -sub, sub - 1

  x = dble(ix)*delta+ dble(lx)/sub*delta/2.0
  y = dble(iy)*delta+ dble(ly)/sub*delta/2.0
  z = dble(iz)*delta+ dble(lz)/sub*delta/2.0

  r2 = (x*x + y*y + z*z) 
  
  if(r2.lt.r_prot2) protvol(ix,iy,iz) = protvol(ix,iy,iz) + 1.d0/dble(subnumb)*cellvol
 enddo !lx
 enddo !ly
 enddo !lz
 sumprotvol = sumprotvol + protvol(ix,iy,iz)
enddo !ix
enddo !iy
enddo !iz

if(rank.eq.0)then
write(6,*) 'volumen protein compute vs geom', sumprotvol,vprot*vsol
endif !rank
!!!!!!!!!!!!!!!!!!!
sumprotpol = 0.d0

do ix = -protR - 1, protR + 1
do iy = -protR - 1, protR + 1
do iz = -protR - 1, protR + 1


 do lx = -sub, sub - 1
 do ly = -sub, sub - 1
 do lz = -sub, sub - 1

  x = dble(ix)*delta+ dble(lx)/sub*delta/2.0
  y = dble(iy)*delta+ dble(ly)/sub*delta/2.0
  z = dble(iz)*delta+ dble(lz)/sub*delta/2.0

  r2 = (x*x + y*y + z*z)

 if(r2.gt.r_prot2.AND.r2.lt.r_protdelta2) protpol(ix,iy,iz) = protpol(ix,iy,iz) + 1.d0/dble(subnumb)


 enddo !lx
 enddo !ly
 enddo !lz
 
  sumprotpol = sumprotpol + protpol(ix,iy,iz)


enddo !ix
enddo !iy
enddo !iz

protpol = protpol * cellvol

if(rank.eq.0) then
write(6,*)'vol de int', sumprotpol*cellvol, 4.d0/3.d0*3.1416*(((protR+1)*delta)**3-(protR*delta)**3)
endif !rank

!test
!call stopundef('rock')

!! matriz interaccion particula-particula
wpartpart=0

do ix = -rcutpart, rcutpart
do iy = -rcutpart, rcutpart
do iz = -rcutpart, rcutpart 

r2 = (ix*ix+iy*iy+iz*iz)*delta*delta

if(r2.gt.d_prot2) then

wpartpart(ix,iy,iz) = -(d_prot2/(r2-d_prot2)+ d_prot2/r2 + 2.d0*log((r2-d_prot2)/r2))/12.d0

endif

enddo !ix
enddo !iy
enddo !iz
call MPI_Barrier(  MPI_COMM_WORLD, ierr)
end subroutine
