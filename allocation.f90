subroutine allocation

use system
use fields_fkfun
use conformations
use chainsdat
use kinsol
use results
use ematrix
use mkinsol
use ellipsoid
use MPI
use kai
use mparameters_monomer
use protein !yamila
implicit none

! fields_fkfun
!ALLOCATE(xtotal(1-Xulimit:dimx+Xulimit, 1-Xulimit:dimy+Xulimit, 1-Xulimit:dimz+Xulimit)) ! xtotal para poor solvent
ALLOCATE(xtotal(dimx, dimy, dimz, N_poorsol))
ALLOCATE(psi(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE(xh(dimx, dimy, dimz))

! kinsol
ALLOCATE (xflag(eqs*dimx*dimy*dimz))
ALLOCATE (xpar(dimx*dimy*dimz))

! results
ALLOCATE (avpol(dimx, dimy, dimz, N_monomer))
ALLOCATE (xpos(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin(dimx, dimy, dimz)) ! OH-
ALLOCATE (fdis(dimx, dimy, dimz, N_monomer))
ALLOCATE (epsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE (Depsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE (xprot(dimx, dimy, dimz)) !yamila
ALLOCATE (fprot(dimx, dimy, dimz)) !yamila

! ematrix
ALLOCATE (volprot(dimx,dimy,dimz))
ALLOCATE (volprot1(dimx,dimy,dimz))
ALLOCATE (voleps(dimx,dimy,dimz))
ALLOCATE (voleps1(dimx,dimy,dimz))
ALLOCATE (volq(dimx,dimy,dimz))
ALLOCATE (volq1(dimx,dimy,dimz))
ALLOCATE (denpart(dimx,dimy,dimz)) !yamila
ALLOCATE (denpart1(dimx,dimy,dimz)) !yamila

! mkinsol
ALLOCATE (pp(eqs*dimx*dimy*dimz))

! chainsdat
allocate(in1(long,3))
allocate(cpp(size))
allocate(cppini(size))

!protein yamila
allocate (protvol(-protR:protR, -protR:protR, -protR:protR))
allocate (protpol(-protR-1:protR+1, -protR-1:protR+1, -protR-1:protR+1))
end subroutine
