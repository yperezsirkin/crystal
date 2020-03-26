subroutine allocateell

use ellipsoid

! ellipsoid
ALLOCATE (rotmatrix(3,3, NNN))
ALLOCATE (Aell(3, NNN))
ALLOCATE (AellS(3, NNN))
ALLOCATE (AellL(3, NNN))
ALLOCATE (AellX(3, NNN))
!ALLOCATE (AellB(3, NNN)) !yamila
ALLOCATE (AAA(3,3, NNN))
ALLOCATE (AAAS(3,3, NNN))
ALLOCATE (AAAL(3,3, NNN))
ALLOCATE (AAAX(3,3, NNN))
!ALLOCATE (AAAB(3,3, NNN)) !yamila
ALLOCATE (Rell(3, NNN))
ALLOCATE (Rellf(3, NNN))
ALLOCATE (orient(3, NNN))
ALLOCATE (echarge(NNN))
ALLOCATE (sigma(NNN))
ALLOCATE (eeps(NNN))

end subroutine
