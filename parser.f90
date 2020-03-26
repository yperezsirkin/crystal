
subroutine readinput
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use transform
use system
use kai
use kaist
use s2d
use channel
use branches
use protein
implicit none

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer ios
integer line, linemax
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi
real*8 ndr

! not defined variables, change if any variable can take the value

seed = 938121
seed2 = 938121
PBC = 1

ndi = -1e5
ndr = -1.0d10

verbose = 5
stdout = 6

electroflag = 1 ! system with electrostatics?

branched = 0 ! branched chains?

sigmar = 0.0 ! random sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
!

vscan = ndi
scx = ndi
scy = ndi
scz = ndi
vtkflag = ndi
systemtype = ndi
dimx = ndi
dimy = ndi
dimz = ndi
long = ndi
cuantas = ndi
readchains = ndi
infile = ndi
randominput = 0
epstype = 0
cutoff = ndr
lseg = ndr
nst = ndi
dielS = ndr
pHbulk = ndr
dielP = ndr
delta = ndr
dx = ndr
dy = ndr
dz = ndr
cdiva = ndr
csalt = ndr
vpol = ndr
fz=ndr 
vsol0 = ndr
gama0 = ndr
benergy = ndr
NNN= ndi 
nsc = 1
scs(1) = 1.0
cprot = ndr !yamila
protR = ndi !yamila
eprotpol = ndr !yamila
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)write(stdout,*) 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


 select case (label)

 case ('PBC')
   read(buffer, *, iostat=ios) PBC(1),PBC(2),PBC(3),PBC(4),PBC(5),PBC(6)
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   do j = 1,5,2
    if((PBC(j).eq.1).and.(PBC(j+1).ne.1)) then 
      write(stdout,*) 'parser:', 'Error in PBC'
      stop
    endif
    if((PBC(j+1).eq.1).and.(PBC(j).ne.1)) then
      write(stdout,*) 'parser:', 'Error in PBC'
      stop
    endif
   enddo

 case ('verbose')
   read(buffer, *, iostat=ios) verbose
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('seed')
   read(buffer, *, iostat=ios) seed2
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('stdout')
   read(buffer, *, iostat=ios) stdout
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vtkflag')
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('electroflag')
   read(buffer, *, iostat=ios) electroflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('branched')
   read(buffer, *, iostat=ios) branched
    
 if(branched.eq.1) then
   read(fh, *) basura
   read(fh, *)longb(1), longb(2), longb(3)
 endif

 if(branched.eq.2) then
   read(fh, *) basura
   read(fh, *)longb(1), longb(2)
 endif


 case ('randominput')
   read(buffer, *, iostat=ios) randominput
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('epstype')
   read(buffer, *, iostat=ios) epstype
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('readchains')
   read(buffer, *, iostat=ios) readchains
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimx')
   read(buffer, *, iostat=ios) dimx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scx')
   read(buffer, *, iostat=ios) scx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scy')
   read(buffer, *, iostat=ios) scy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scz')
   read(buffer, *, iostat=ios) scz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('delta')
   read(buffer, *, iostat=ios) delta
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dx')
   read(buffer, *, iostat=ios) dx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dy')
   read(buffer, *, iostat=ios) dy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dz')
   read(buffer, *, iostat=ios) dz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cdiva')
   read(buffer, *, iostat=ios) cdiva
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('dimy')
   read(buffer, *, iostat=ios) dimy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimz')
   read(buffer, *, iostat=ios) dimz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('long')
   read(buffer, *, iostat=ios) long
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cuantas')
   read(buffer, *, iostat=ios) cuantas
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('lseg')
   read(buffer, *, iostat=ios) lseg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dielS')
   read(buffer, *, iostat=ios) dielS
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('csalt')
   read(buffer, *, iostat=ios) csalt
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vsol')
   read(buffer, *, iostat=ios) vsol0
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('benergy')
   read(buffer, *, iostat=ios) benergy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vpol')
   read(buffer, *, iostat=ios) vpol
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vscan')
   read(buffer, *, iostat=ios) vscan
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('sigmar')
   read(buffer, *, iostat=ios) sigmar
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('gama')
   read(buffer, *, iostat=ios) gama0
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('infile')
   read(buffer, *, iostat=ios) infile
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('hguess')
   read(buffer, *, iostat=ios) hguess
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('hring')
   read(buffer, *, iostat=ios) hring
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('oval')
   read(buffer, *, iostat=ios) oval
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('nkp')
   read(buffer, *, iostat=ios) nkp
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   do i = 1, nkp
   read(fh,*)kps(i)
   enddo

 case ('nst')
   read(buffer, *, iostat=ios) nst
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  
   do i = 1, nst
   read(fh,*)sts(i)
   enddo 


 case ('nsc')
   read(buffer, *, iostat=ios) nsc
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  
   do i = 1, nsc
   read(fh,*)scs(i)
   enddo 

 case ('Xucutoff')
   read(buffer, *, iostat=ios) cutoff
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case('fz') 
   read(buffer, *, iostat=ios) fz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 
 case ('cprot') !yamila
   read(buffer, *, iostat=ios) cprot
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('protR') !yamila
   read(buffer, *, iostat=ios) protR
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('eprotpol') !yamila
   read(buffer, *, iostat=ios) eprotpol
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('systemtype')
   read(buffer, *, iostat=ios) systemtype
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   select case (systemtype) ! TYPE OF SYSTEM
                            ! TYPE = 1 is nanoparticle crystal
                            ! TYPE = 2 is 3D channel
                            ! TYPE = 3 is a 3D channel with chains at specific conditions

    case(2)
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) sigmac
     read(fh, *) basura
     read(fh, *) echargec
     read(fh, *) basura
     read(fh, *) eepsc
  
    case(3, 4, 41)
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) echargec
     read(fh, *) basura
     read(fh, *) eepsc

    case(6) ! planar surface
     read(fh, *) basura
     read(fh, *) Npolx, Npoly ! number of polymers in x and y
     read(fh, *) basura
     read(fh, *) eepsc
     NNN = 0 ! no particles

 
    case(42, 52) ! 42: channel, 52: rod
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) Nrings

     allocate (ringpos(Nrings))

     read(fh, *) basura
     do i = 1, Nrings
       read(fh, *) ringpos(i)
     enddo
      ringpos = ringpos - 0.5
    
     read(fh, *) basura
     read(fh, *) echargec
     read(fh, *) basura
     read(fh, *) eepsc

    case(60) ! channel + particle
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) Nrings

     allocate (ringpos(Nrings))

     read(fh, *) basura
     do i = 1, Nrings
       read(fh, *) ringpos(i)
     enddo
      ringpos = ringpos - 0.5
     read(fh,*) basura
     read(fh,*) NNN ! number of particles  		
  !    NNN = 1 ! only one particle
     call allocateell

     read(fh, *) basura  
     do j = 1, NNN
     read(fh, *) Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo

    

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *) basura

     do j = 1, NNN
     read(fh, *) rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *) rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *) rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         write(stdout,*) 'parser:','Set particle',j,'rotation to:'
         write(stdout,*) 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         write(stdout,*) 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         write(stdout,*) 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) echarge(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'charge to', echarge(j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) eeps(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

    case(1) 
     read(fh, *) basura
     read(fh, *)NNN

     if(NNN.ne.0) then

     call allocateell
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *) rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *) rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         write(stdout,*) 'parser:','Set particle',j,'rotation to:'
         write(stdout,*) 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         write(stdout,*) 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         write(stdout,*) 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) sigma(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'surface coverage to', sigma(j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) echarge(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'charge to', echarge(j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) eeps(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     endif ! NNN
endselect
endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
! 

if(systemtype.eq.2) then
 if((cdiva.ne.1.0).or.(gama0.ne.90.0)) then
  write(stdout,*) 'Channel works only for cdiva = 1 and gama0 = 90.0... ending'
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
 endif
endif

if (branched.eq.1) then
 longbb = long
 long = longbb + longb(1) + longb(2) + longb(3)
endif 

if (branched.eq.2) then
 longbb = long
 long = long + longb(2) - longb(1)
endif 


if(vtkflag.eq.ndi)call stopundef('vtkflag')
if(dimx.eq.ndi)call stopundef('dimx')
if(scx.eq.ndi)call stopundef('scx')
if(scy.eq.ndi)call stopundef('scy')
if(scz.eq.ndi)call stopundef('scz')
if(dimy.eq.ndi)call stopundef('dimy')
if(dimz.eq.ndi)call stopundef('dimz')
if(ncha.eq.ndi)call stopundef('ncha')
if(long.eq.ndi)call stopundef('long')
if(cuantas.eq.ndi)call stopundef('cuantas')
if(infile.eq.ndi)call stopundef('infile')
if(cutoff.eq.ndr)call stopundef('Xucutoff')
if(readchains.eq.ndi)call stopundef('readchains')
if(systemtype.eq.ndi)call stopundef('systemtype')
if(nst.eq.ndi)call stopundef('nst')

if(delta.eq.ndr)call stopundef('delta')
if(dx.eq.ndr)call stopundef('dx')
if(dy.eq.ndr)call stopundef('dy')
if(dz.eq.ndr)call stopundef('dz')
if(cdiva.eq.ndr)call stopundef('cdiva')
if(dielS.eq.ndr)call stopundef('dielS')
if(dielP.eq.ndr)call stopundef('dielP')
if(lseg.eq.ndr)call stopundef('lseg')
if(csalt.eq.ndr)call stopundef('csalt')
if(pHbulk.eq.ndr)call stopundef('pHbulk')
if(vpol0.eq.ndr)call stopundef('vpol')
if(vsol0.eq.ndr)call stopundef('vsol')
if(benergy.eq.ndr)call stopundef('benergy')
if(gama0.eq.ndr)call stopundef('gama')
if(fz.eq.ndr)call stopundef('fz') 
if(NNN.eq.ndi) call stopundef('NNN') 
if(cprot.eq.ndr)call stopundef('cprot') ! yamila
if(protR.eq.ndi)call stopundef('protR') ! yamila
if(eprotpol.eq.ndr)call stopundef('cprot') ! yamila
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine stopundef(namevar)
use const
character(len=*) :: namevar
write(stdout,*) 'parser:', 'Variable ', namevar, ' is undefined '
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

