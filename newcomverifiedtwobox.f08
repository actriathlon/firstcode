program move

  implicit none
integer::topprotcoords,totpoints,timeofinterest
  double precision :: totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength,choose2
  double precision :: actualchainlength,random,totiso,sumsep,bulksumsep
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug,testdummy,polymerrog,polymerl
  integer :: rept,datayes,maxclussize,minl,outputrate,backtime,magiccontent,AB,AC,BC,ABC,bAB,bAC,bBC,bABC
  integer ::lastruntime,timebase,steplength,macc,mrej,roatt,roacc,clientlength,telacc,telrej
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count,avecluspop,aveclusnumber,vtype,gridsizex
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib,cltacc,cltatt,binsize,normaliser
  integer :: successful,reject,maxlength1,maxlength2,q,d,u,sm,qt,ft,gt,scans,nclients,phasecount
  integer:: rerej,reacc,nprotein1,nprotein2,nspecies,mtype,freezetime,nclientspecies,nscaffold
  integer,dimension(:,:),allocatable ::bonddd,study,sumclient,visit,visit2,unbound
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen,clientlist,specieslinker
  integer,dimension(:),allocatable:: clnos,totalbindevents
  !integer,dimension(:),allocatable::
  double precision :: trackx,tracky,trackz,intraen,intdummy,totdire,r
  double precision,dimension(:),allocatable::chen
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy
  real, external :: ran2
!character(len = 10) :: commandread  
integer::preinterest,ttt
  character(len = 10) ::commandread1
integer,dimension(4):: clusneigh,bulkneigh

integer,dimension(:),allocatable::locb,loc
real,dimension(:),allocatable :: variance,disp
  real,dimension(:,:),allocatable:: summing,vdist,freesites,clusdist
  logical :: exist,fail,finalfail,debugyes,film,isbond,clt,clr,noinfo,scalinginfo,nobonds,restart
  real::start,finish,midpoint,choose1
  !integer::removethis,pivde
  logical :: suffclust,settingup,twosys
        real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
integer::pointless1,pointless2,dright,dleft,noROGS
double precision:: averageROG
integer,dimension(:),allocatable::totalfree

  !cluster COM may fall foul of pbcs

  !Sort MSD outputs

  
  type protein
     integer :: x,y,z,species,type,linker,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl
  end type protein

  type rprot
     real :: x,y,z
  end type rprot

  type basicp
     integer :: x,y,z
  end type basicp

  type prottemp
     integer :: x,y,z,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl
  end type prottemp

  type centremass
     double precision :: x,y,z
  end type centremass

  type lbond
integer:: bm,bl
end type lbond

  type(protein),dimension(:,:),allocatable :: protcoords,client
  type(centremass),dimension(:,:),allocatable :: com
   type(centremass) :: dpcomclusterinit


  call get_command_argument(1,commandread1)
  read(commandread1, '(i4)') timeofinterest




sumsep = 0.0d0
bulksumsep = 0.0d0
AB=0
AC=0
BC = 0
ABC = 0
bAB=0
bAC=0
bBC = 0
bABC = 0
telacc=0
telrej = 0

pointless1= 0
pointless2=0
dright = 0
dleft = 0
noROGS = 0
averageROG = 0.0d0

allocate(locb(2))
locb(1) = 0
locb(2) = 0
allocate(totalbindevents(3))
totalbindevents(:) = 0
  allocate(client(1,5))
  mtype = 0
vtype=0        
!noinfo = .false.
  call read_setup

gridsizex = 2*gridsize
 !call get_command_argument(1,commandread)
  !read(commandread, '(i3)') freezetime

clusneigh(:) = 0
bulkneigh(:) = 0

  write(6,*) 'CLIENTINFO',clientlength,nclients!,freezetime 
  allocate(sumclient(nclients,2*(clientlength+1)))
allocate(unbound(nclients,clientlength))
unbound(:,:) = 0  
sumclient(:,:) = 0
suffclust = .false.
  nprotein = sum(speciespop)
  write(6,*) 'nprotein = ', nprotein
  maxlength = maxval(specieslen)
  write(6,*) 'maxlength =',maxlength
  open(21, file = 'newmovetagged.vtf', action = 'read')



  
  finalfail = .false.
  maxclussize = 1
  

if(int(totiso) == 0) isbond = .false.
  if(int(totdire) == 0) nobonds =  .true.

  !interenergy = -1.0

  !write(6,*) maxlength1,maxlength2
  write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT
  write(6,*) 'clr',clr,'clt',clt



  allocate(protcoords(nprotein,maxlength))
  !allocate(protcoords(nprotein,maxlength))
  allocate(com(2,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(bonddd(nprotein,maxlength))
  allocate(chlen(nprotein))
  allocate(chen(nprotein))
  allocate(clientlist(nclients))
  allocate(visit(nprotein,maxlength))
  allocate(clnos(nprotein))
  allocate(vdist(4*clientlength,gridsizex))
allocate(clusdist(2*mtype,gridsizex))
allocate(loc(2))
allocate(totalfree(2*mtype))

totalfree(:) = 0


loc(:) = 0
  phasecount = 0
vdist(:,:) = 0
  visit(:,:) = 0
!summing(:,:) = 0
clusdist(:,:) = 0
  
binsize = 1

  N = maxlength


  do qt = 1,nclients
     clientlist(qt) = (nprotein -nclients) +qt
end do
  count = 0
  time = 0
  totdisp = 0.0d0


  WRITE(6,*) 'ISBOND =',isbond
  write(6,*) 'interen',interen
  write(6,*) 'intraen',intraenergy


  macc = 0
  mrej = 0
  roatt = 0
  roacc = 0
nscaffold = nprotein-nclients

write(6,*) 'NSCAFFOLD',nscaffold

allocate(study(nclients,maxlength))
  !sets the limit of the pivot move


  call CPU_TIME(start)

  write(6,*) 'foundation start'
     
 call foundation

  write(6,*) 'foundation complete'

  minl = minval(specieslen)


  allocate(summing(2,(protcoords(nprotein,1)%linker*binsize)+1))
  summing(:,:) = 0
  allocate(freesites(mtype*6,gridsizex*binsize))
  freesites(:,:) = 0

  do qt = 1,nprotein
     call comfind(qt,.false.)
     com(2,qt)%x = com(1,qt)%x
     com(2,qt)%y = com(1,qt)%y
     com(2,qt)%z = com(1,qt)%z
  end do
  write(6,*) 'b'


  settingup = .true.
  !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',content
settingup = .false.
normaliser = 0
call freesitescalc
!call neighbour

write(6,*) 'centre of mass details:',dpcomclusterinit

  
  actualchainlength = 0.0

do ttt = 1,timeofinterest-1

call nextstep
 call freesitescalc
!call neighbour
end do

do ft = 1,nscaffold
do qt = 1,chlen(ft)
if(protcoords(ft,qt)%x > gridsize) dright = dright+1
if(protcoords(ft,qt)%x <= gridsize) dleft = dleft+1
end do
end do


write(1771,*) clusneigh(2),bulkneigh(2),totalfree(2),totalfree(4)


contains

subroutine neighbour(bound)
integer::m,l,maxl,xp,yp,zp,dummymax
    integer,dimension(:,:),intent(in)::bound

dummymax = speciespop(1) + speciespop(2)

do m = 1,dummymax,1
maxl = chlen(m)
do l = 1,maxl
if (bound(m,l) == 0) then
if(protcoords(m,l)%am == 0)then
 bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%bm == 0) then
bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%cm == 0) then
bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%dm == 0) then
bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%em == 0) then
bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%fm == 0) then
bulkneigh(protcoords(m,l)%type) = bulkneigh(protcoords(m,l)%type)+1
end if
end if
end do
end do



do m = dummymax+1,nprotein-1
maxl = chlen(m)
do l = 1,maxl
if(bound(m,l) == 0) then
if(protcoords(m,l)%am == 0)then
 clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%bm == 0) then
clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%cm == 0) then
clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%dm == 0) then
clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%em == 0) then
clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
if(protcoords(m,l)%fm == 0) then
clusneigh(protcoords(m,l)%type) = clusneigh(protcoords(m,l)%type)+1
end if
end if
end do
end do


end subroutine neighbour




  subroutine restdata
    integer::m,l,maxl,binsize,lx
    character(len=10)::BIN
    open(119,file='visit.dat',action  = 'read')
    open(108,file='bonddate.dat',action  = 'read')

    do m= 1,nprotein,1
       maxl = chlen(m)
       do l = 1,maxl
          read(119,*) BIN,BIN,visit(m,l)
       end do
    end do
    read(108,*)sumclient(:,:),BIN,BIN


        close(108)
        close(119)
  end subroutine restdata



 subroutine rdftot(m)
   integer,intent(in)::m
    integer :: l,maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


!phasecountdoub = phasecountdoub+1

       maxl = chlen(m)
       do l = 1,maxl
          !phasecount = phasecount+1
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsizex - abs(protcoords(m,l)%x - dpcomclusterinit%x))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - dpcomclusterinit%y))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - dpcomclusterinit%z))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


   vdist(clientlength+l,lx+1) = vdist(clientlength+l,lx+1)+1.0
   !vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0
       end do

    
  end subroutine rdftot


  

 subroutine rdf(m,l,left)
   integer,intent(in)::m,l
logical,intent(in)::left
    integer :: maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


phasecount = phasecount+1

if(left .eqv. .false.) then

       maxl = chlen(m)
          dx = min(abs(protcoords(m,l)%x - ((3.0/4.0)*gridsize)), gridsizex - abs(protcoords(m,l)%x - ((3.0/4.0)*gridsize)))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - (gridsize/2.0)))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - (gridsize/2.0)))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            vdist(l,lx+1) = vdist(l,lx+1)+1.0
   vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0

else if(left .eqv. .true.) then

       maxl = chlen(m)
          dx = min(abs(protcoords(m,l)%x - (gridsize/2.0)), gridsizex - abs(protcoords(m,l)%x - ((3.0/4.0)*gridsize)))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - (gridsize/2.0)))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - (gridsize/2.0)))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            vdist(l,lx+1) = vdist(l,lx+1)+1.0
   vdist((3*clientlength)+1,lx+1) = vdist((3*clientlength)+1,lx+1) + 1.0


end if

    
 end subroutine rdf


 
 subroutine rdfcomplex(m,control,attached)
   integer,intent(in)::m,control
   integer,dimension(:),intent(in)::attached
    integer :: l,maxl,lx,j
    double precision :: dx,dy,dz
    double precision :: delta


    do j = 1,control,1
       l = attached(j)
!phasecount = phasecount+1

       maxl = chlen(m)
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsizex - abs(protcoords(m,l)%x - dpcomclusterinit%x))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - dpcomclusterinit%y))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - dpcomclusterinit%z))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            vdist((3*clientlength)+l,lx+1) = vdist((3*clientlength)+l,lx+1)+1.0
   !vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0
end do

    
  end subroutine rdfcomplex


 subroutine rdfcluster(m,control,attached)
   integer,intent(in)::m,control
   integer,dimension(:),intent(in)::attached
    integer :: l,maxl,lx,j
    double precision :: dx,dy,dz
    double precision :: delta


    do j = 1,control,1
       l = attached(j)
!phasecount = phasecount+1

       maxl = chlen(m)
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsizex - abs(protcoords(m,l)%x - dpcomclusterinit%x))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - dpcomclusterinit%y))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - dpcomclusterinit%z))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            clusdist(protcoords(m,l)%type,lx+1) = clusdist(protcoords(m,l)%type,lx+1)+1.0
   !vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0
end do

    
  end subroutine rdfcluster

  

subroutine foundation
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z,tet
    character(len = 10) :: BIN,BIN8,nonsensea,nonsenseb,nonsensec
integer::dummymax,g,check,bondsneeded
doubleprecision :: rx,ry,rz,bondx,bondy,bondz
tet = 1
check = 0
topprotcoords = 0
  totpoints = 0
dummymax=0
  write(6,*) 'nspecies',nspecies,speciespop(5),specieslen(5)
  do n=1,5,1
	topprotcoords = topprotcoords + speciespop(n)
     totpoints = totpoints + (speciespop(n)*specieslen(n))
  end do

dummymax= speciespop(1) + speciespop(2)

  write(6,*) 'totpoints',totpoints,nprotein,nclients
  
  !read(21,*) BIN

totpoints = totpoints*2


write(6,*) 'totpoints post double',totpoints,topprotcoords

 read(21,*) BIN,BIN,BIN,BIN,BIN,BIN
 do m = 1,topprotcoords,1 !nprotein-nclients

!write(6,*) 'm',tet,m
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,1)%type
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        !write(6,*) 'BIN8 ', BIN8
        !tet = int(BIN8)
        !if((BIN8 /= 1) .and. (BIN8 /=2)) write(6,*) 'found it!!!!',BIN8

if(m<=(dummymax)) then

protcoords(m,1)%species = protcoords(m,1)%type
     maxl = specieslen(protcoords(m,1)%type)
     chlen(m) = specieslen(protcoords(m,1)%type)
check = check+maxl

else if((m>dummymax) .and. (m<topprotcoords)) then
        protcoords(m,1)%species = protcoords(m,1)%type
        maxl = specieslen(protcoords(m,1)%type+2)
        chlen(m) = specieslen(protcoords(m,1)%type+2)
check = check+maxl

else if (m==topprotcoords) then

        protcoords(m,1)%species = protcoords(m,1)%type
        maxl = specieslen(protcoords(m,1)%type+4)
        chlen(m) = specieslen(protcoords(m,1)%type+4)
!write(6,*) 'here is visited!!'
write(6,*) 'check',check
check = check+maxl
write(6,*) 'check',check

end if


     count = count+1
tet = tet+maxl
!write(6,*) 'count',count
     do l = 2,maxl
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,l)%type
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        protcoords(m,l)%species = protcoords(m,l)%type
        count = count+1
     end do

  end do
write(6,*) 'tet',tet,totpoints,check
 ! dummymax = totpoints - topprotcoords !(nprotein-nclients)
!read(21,*) BIN,BIN8
!read(21,*) BIN
!write(6,*) 'Transition ',BIN , BIN8
!STOP
 

bondsneeded=0
do n = 1,5,1

bondsneeded = bondsneeded+(((2*specieslen(n))-1)*speciespop(n))
end do

 do m = 1,bondsneeded
    read(21,*) BIN8,BIN,BIN
!        write(6,*) 'BIN8 ',BIN8,m
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !accounts for 2 extra lines in energy than in move
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
     read(21,*) BIN,BIN8
        write(6,*) 'BIN',BIN,BIN8
     read(21,*) BIN
  
     do m = 1,topprotcoords,1 !nprotein-nclients,1
        read(21,*)  BIN,rx,ry,rz
        !write(6,*) BIN
        read(21,*) BIN,bondx,bondy,bondz
protcoords(m,1)%linker = specieslinker(protcoords(m,1)%type)
        protcoords(m,1)%x = rx/4
        protcoords(m,1)%y = ry/4
        protcoords(m,1)%z = rz/4


        if((protcoords(m,1)%x-(bondx/4.0) == -0.5) .or. (protcoords(m,1)%x-(bondx/4.0) == (gridsizex-0.5))) bonddd(m,1) = 1
        if((protcoords(m,1)%x-(bondx/4.0) == 0.5) .or. (protcoords(m,1)%x-(bondx/4.0) == -(gridsizex-0.5))) bonddd(m,1) = -1
        if((protcoords(m,1)%y-(bondy/4.0) == -0.5) .or. (protcoords(m,1)%y-(bondy/4.0) == (gridsize-0.5))) bonddd(m,1) = 2
        if((protcoords(m,1)%y-(bondy/4.0) == 0.5) .or. (protcoords(m,1)%y-(bondy/4.0) == -(gridsize-0.5))) bonddd(m,1) = -2
        if((protcoords(m,1)%z-(bondz/4.0) == -0.5) .or. (protcoords(m,1)%z-(bondz/4.0) == (gridsize-0.5))) bonddd(m,1) = 3
        if((protcoords(m,1)%z-(bondz/4.0) == 0.5) .or. (protcoords(m,1)%z-(bondz/4.0) == -(gridsize-0.5))) bonddd(m,1) = -3
        !chlen(m) = specieslen(protcoords(m,1)%type)
        maxl = chlen(m)

        if(bonddd(m,1) ==0) write(6,*) 'no bond allocation',m,1
        do l = 2,maxl,1
        protcoords(m,l)%linker = specieslinker(protcoords(m,l)%type)
           read(21,*)  BIN,rx,ry,rz
           read(21,*) BIN,bondx,bondy,bondz

           protcoords(m,l)%x = rx/4
           protcoords(m,l)%y = ry/4
           protcoords(m,l)%z = rz/4

           if((protcoords(m,l)%x-(bondx/4.0) == -0.5) .or. (protcoords(m,l)%x-(bondx/4.0) == (gridsizex-0.5))) bonddd(m,l) = 1
           if((protcoords(m,l)%x-(bondx/4.0) == 0.5) .or. (protcoords(m,l)%x-(bondx/4.0) == -(gridsizex-0.5))) bonddd(m,l) = -1
           if((protcoords(m,l)%y-(bondy/4.0) == -0.5) .or. (protcoords(m,l)%y-(bondy/4.0) == (gridsize-0.5))) bonddd(m,l) = 2
           if((protcoords(m,l)%y-(bondy/4.0) == 0.5) .or. (protcoords(m,l)%y-(bondy/4.0) == -(gridsize-0.5))) bonddd(m,l) = -2
           if((protcoords(m,l)%z-(bondz/4.0) == -0.5) .or. (protcoords(m,l)%z-(bondz/4.0) == (gridsize-0.5))) bonddd(m,l) = 3
           if((protcoords(m,l)%z-(bondz/4.0) == 0.5) .or. (protcoords(m,l)%z-(bondz/4.0) == -(gridsize-0.5))) bonddd(m,l) = -3
           !write(6,*) 'bond',protcoords(m,l)%bond
           if(bonddd(m,l) ==0) write(6,*) 'no bond allocation',m,l
        end do
     end do

    call debug



  end subroutine foundation

subroutine nextstep
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z
    character(len = 10) :: BIN,nonsensea,nonsenseb,nonsensec
integer::dummymax,g
doubleprecision :: rx,ry,rz,bondx,bondy,bondz




   read(21,*) BIN
        write(6,*) 'BIN',BIN,topprotcoords
     read(21,*) BIN

     do m = 1,topprotcoords,1 !nprotein-nclients,1
         
       read(21,*)  BIN,rx,ry,rz
        read(21,*) BIN,bondx,bondy,bondz
        protcoords(m,1)%x = rx/4
        protcoords(m,1)%y = ry/4
        protcoords(m,1)%z = rz/4


        if((protcoords(m,1)%x-(bondx/4.0) == -0.5) .or.(protcoords(m,1)%x-(bondx/4.0) == (gridsizex-0.5))) bonddd(m,1) = 1
        if((protcoords(m,1)%x-(bondx/4.0) == 0.5) .or.(protcoords(m,1)%x-(bondx/4.0) == -(gridsizex-0.5))) bonddd(m,1) = -1
        if((protcoords(m,1)%y-(bondy/4.0) == -0.5) .or.(protcoords(m,1)%y-(bondy/4.0) == (gridsize-0.5))) bonddd(m,1) = 2
        if((protcoords(m,1)%y-(bondy/4.0) == 0.5) .or.(protcoords(m,1)%y-(bondy/4.0) == -(gridsize-0.5))) bonddd(m,1) = -2
        if((protcoords(m,1)%z-(bondz/4.0) == -0.5) .or.(protcoords(m,1)%z-(bondz/4.0) == (gridsize-0.5))) bonddd(m,1) = 3
        if((protcoords(m,1)%z-(bondz/4.0) == 0.5) .or.(protcoords(m,1)%z-(bondz/4.0) == -(gridsize-0.5))) bonddd(m,1) = -3
        chlen(m) = specieslen(protcoords(m,1)%type)
        maxl = chlen(m)

        if(bonddd(m,1) ==0) write(6,*) 'no bond allocation',m,1
        do l = 2,maxl,1
           read(21,*)  BIN,rx,ry,rz
           read(21,*) BIN,bondx,bondy,bondz

           protcoords(m,l)%x = rx/4
           protcoords(m,l)%y = ry/4
           protcoords(m,l)%z = rz/4

           if((protcoords(m,l)%x-(bondx/4.0) == -0.5) .or.(protcoords(m,l)%x-(bondx/4.0) == (gridsizex-0.5))) bonddd(m,l) = 1
           if((protcoords(m,l)%x-(bondx/4.0) == 0.5) .or.(protcoords(m,l)%x-(bondx/4.0) == -(gridsizex-0.5))) bonddd(m,l) = -1
           if((protcoords(m,l)%y-(bondy/4.0) == -0.5) .or.(protcoords(m,l)%y-(bondy/4.0) == (gridsize-0.5))) bonddd(m,l) = 2
           if((protcoords(m,l)%y-(bondy/4.0) == 0.5) .or.(protcoords(m,l)%y-(bondy/4.0) == -(gridsize-0.5))) bonddd(m,l) = -2
           if((protcoords(m,l)%z-(bondz/4.0) == -0.5) .or.(protcoords(m,l)%z-(bondz/4.0) == (gridsize-0.5))) bonddd(m,l) = 3
           if((protcoords(m,l)%z-(bondz/4.0) == 0.5) .or.(protcoords(m,l)%z-(bondz/4.0) == -(gridsize-0.5))) bonddd(m,l) = -3
           !write(6,*) 'bond',protcoords(m,l)%bond
           if(bonddd(m,l) ==0) write(6,*) 'no bond allocation',m,l
        end do
     end do

    call debug

end subroutine nextstep

  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,tbo,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    integer,dimension(:),intent(in)::tbo
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding


    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
       bonddd(chainnum,beadnum) = tbo(beadnum)
     end do

    call comfind(chainnum,.TRUE.)


    !call rms(chainnum)
  end subroutine updatepos

  
  subroutine removeenergyintra(m,l,deltaenergy,tempcoord)
    integer,intent(in)::m,l
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout)::tempcoord
    integer :: g,f

    g = m
    if(protcoords(m,l)%am == m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%al
       tempcoord(l)%am = -1
       if(bonddd(m,l) == 1 .and. bonddd(g,f) == -1)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if(protcoords(m,l)%bm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%bl
       tempcoord(l)%bm = -1
       if(bonddd(m,l) == -1 .and. bonddd(g,f) == 1)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if(protcoords(m,l)%cm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
       if(bonddd(m,l) == 2 .and. bonddd(g,f) == -2)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if(protcoords(m,l)%dm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%dl
       tempcoord(l)%dm = -1
       if(bonddd(m,l) == -2 .and. bonddd(g,f) == 2)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if(protcoords(m,l)%em ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%el
       tempcoord(l)%em = -1
       if(bonddd(m,l) == 3 .and. bonddd(g,f) == -3)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if(protcoords(m,l)%fm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%fl
       tempcoord(l)%fm = -1
       if(bonddd(m,l) == -3 .and. bonddd(g,f) == 3)   deltaenergy = &
            deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

  end subroutine removeenergyintra


  subroutine removeenergyinter(m,l,deltaenergy,tempcoord)
    integer,intent(in)::m,l
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout)::tempcoord
    integer :: g,f

    !write(6,*) 'a'
    if((protcoords(m,l)%am /= 0) .and. (protcoords(m,l)%am /= m)) then
       g = protcoords(m,l)%am
       f = protcoords(m,l)%al
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%am = -1
       if((bonddd(m,l) == 1) .and. (bonddd(g,f) == -1))then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

    if((protcoords(m,l)%bm /= 0)  .and. (protcoords(m,l)%bm /= m)) then
       g = protcoords(m,l)%bm
       f = protcoords(m,l)%bl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%bm = -1
       if((bonddd(m,l) == -1) .and. (bonddd(g,f) == 1)) then
          deltaenergy = deltaenergy  - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

    if((protcoords(m,l)%cm /= 0) .and. (protcoords(m,l)%cm /= m)) then
       g = protcoords(m,l)%cm
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       if((bonddd(m,l) == 2) .and. (bonddd(g,f) == -2)) then
          deltaenergy = deltaenergy  - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

    if((protcoords(m,l)%dm /= 0) .and. (protcoords(m,l)%dm /= m)) then
       g = protcoords(m,l)%dm
       f = protcoords(m,l)%dl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%dm = -1
       if((bonddd(m,l) == -2) .and. (bonddd(g,f) == 2)) then
          deltaenergy = deltaenergy  - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

    if((protcoords(m,l)%em /= 0) .and. (protcoords(m,l)%em /= m)) then
       g = protcoords(m,l)%em
       f = protcoords(m,l)%el
       !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1)) write(6,*) 'successful bond rupture',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%em = -1
       if((bonddd(m,l) == 3) .and. (bonddd(g,f) == -3)) then
          deltaenergy = deltaenergy  - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

    if((protcoords(m,l)%fm /= 0) .and. (protcoords(m,l)%fm /= m)) then
       g = protcoords(m,l)%fm
       f = protcoords(m,l)%fl
       !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'successful rupture of bond',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%fm = -1
       if((bonddd(m,l) == -3) .and. (bonddd(g,f) == 3)) then
          deltaenergy = deltaenergy  - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          !write(6,*) 'bond break',m,l,g,f
       end if
    end if

  end subroutine removeenergyinter

  subroutine bondformintra(m,l,st,tempcoord,tbo,bdir,deltaenergy)
    integer,intent(in)::m,l,st,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    integer,dimension(:),allocatable::tbo
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'start fail'
    if(bdir == -1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%am = m
       tempcoord(l)%al = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if

    if(bdir == 1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%bm = m
       tempcoord(l)%bl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if

    if(bdir == -2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%cm = m
       tempcoord(l)%cl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if

    if(bdir == 2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%dm = m
       tempcoord(l)%dl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if

    if(bdir == -3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%em = m
       tempcoord(l)%el = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if

    if(bdir == 3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%fm = m
       tempcoord(l)%fl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(m,st))) then
          deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)
       end if
    end if
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'fail'
  end subroutine bondformintra

  subroutine bondforminter(m,g,l,st,tempcoord,tbo,bdir,deltaenergy)
    integer,intent(in)::m,l,st,g,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    integer,dimension(:),allocatable::tbo
    !check this and cluster bonds update


    !write(6,*) 'm',m,'l',l,'g',g,'st',st
    deltaenergy = deltaenergy + interen(protcoords(m,l)%type,protcoords(g,st)%type)
    if(bdir == -1) then
       tempcoord(l)%am = g
       tempcoord(l)%al = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then
          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    if(bdir == 1) then
       tempcoord(l)%bm = g
       tempcoord(l)%bl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then
          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    if(bdir == -2) then
       tempcoord(l)%cm = g
       tempcoord(l)%cl = st
       !if(modulo(time,1000)==0) write(6,*) 'interform',bdir,tbo(l),bonddd(g,st)
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then
          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    if(bdir == 2) then
       tempcoord(l)%dm = g
       tempcoord(l)%dl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then
          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    if(bdir == -3) then
       tempcoord(l)%em = g
       tempcoord(l)%el = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then

          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    if(bdir == 3) then
       tempcoord(l)%fm = g
       tempcoord(l)%fl = st
       if((bdir== -1*tbo(l)) .and. (bdir== 1*bonddd(g,st))) then
          deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,st)%type)
       end if
    end if

    !if(modulo(time,1000)==0) write(6,*) 'interbond',bdir,deltaenergy,tbo(l)
  end subroutine bondforminter










  subroutine bonddirection(tempbonddd)
    integer,intent(inout) :: tempbonddd
    integer :: xp,xn,yp,yn,zp,zn,chdir 

    chdir = int(ran2(seed)*6)+1


    if(chdir == 1) tempbonddd = 1
     if(chdir == 2) tempbonddd = -1
    if(chdir == 3) tempbonddd = 2
    if(chdir == 4) tempbonddd = -2
    if(chdir == 5) tempbonddd = 3
    if(chdir == 6) tempbonddd = -3

  end subroutine bonddirection

  

  logical Function overlaps (pr,ch1,ch2,tempcoords)
    integer,intent(in) :: pr,ch1,ch2
    Type(prottemp),dimension(:),intent(in) :: tempcoords
    overlaps = .true.
    if((protcoords(pr,ch2)%x == tempcoords(ch1)%x) .and. &
         (protcoords(pr,ch2)%y == tempcoords(ch1)%y) .and. &
         (protcoords(pr,ch2)%z == tempcoords(ch1)%z)) then
       overlaps = .false.
    end if
  end Function Overlaps

  logical Function overlapsclus (m,pr,ch1,ch2,tempcluscoord)
    integer,intent(in) :: pr,ch1,ch2,m
    Type(prottemp),dimension(:,:),intent(in) :: tempcluscoord
    overlapsclus = .true.
    if((protcoords(pr,ch2)%x == tempcluscoord(m,ch1)%x) .and. &
         (protcoords(pr,ch2)%y == tempcluscoord(m,ch1)%y) .and. &
         (protcoords(pr,ch2)%z == tempcluscoord(m,ch1)%z)) then
       overlapsclus = .false.
    end if
  end Function Overlapsclus

 subroutine updatebondold(m,beadmin,beadmax,tempcoord)
    integer,intent(in):: m,beadmin,beadmax
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f,bg,bf



    do l = beadmin,beadmax,1


          protcoords(m,l)%am = tempcoord(l)%am
          protcoords(m,l)%al = tempcoord(l)%al
          protcoords(m,l)%bm = tempcoord(l)%bm
          protcoords(m,l)%bl = tempcoord(l)%bl
          protcoords(m,l)%cm = tempcoord(l)%cm
          protcoords(m,l)%cl = tempcoord(l)%cl
          protcoords(m,l)%dm = tempcoord(l)%dm
          protcoords(m,l)%dl = tempcoord(l)%dl
          protcoords(m,l)%em = tempcoord(l)%em
          protcoords(m,l)%el = tempcoord(l)%el
          protcoords(m,l)%fm = tempcoord(l)%fm
          protcoords(m,l)%fl = tempcoord(l)%fl

    end do

  end subroutine updatebondold






  subroutine energy(init)
    integer :: dx,dy,dz,m,l,f,g,bdir,delx,dely,delz,maxlengthss,maxl
    double precision :: initialenergy,olderenergy,dumoldenergy
    Type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: tbo
  !double precision,dimension(:),intent(inout)::chen
    logical :: adjver
    logical,intent(in)::init
    initialenergy = 0.0d0
    !totalenergy = 0.0d0
    allocate(tempcoord(maxlength))
    allocate(tbo(maxlength))


    if(init .eqv. .true.) then
       do m = 1,nprotein
          maxl = chlen(m)
          do l = 1,maxl
             protcoords(m,l)%am = 0
             protcoords(m,l)%al = 0
             protcoords(m,l)%bm = 0
             protcoords(m,l)%bl = 0
             protcoords(m,l)%cm = 0
             protcoords(m,l)%cl = 0
             protcoords(m,l)%dm = 0
             protcoords(m,l)%dl = 0
             protcoords(m,l)%em = 0
             protcoords(m,l)%el = 0
             protcoords(m,l)%fm = 0
             protcoords(m,l)%fl = 0

             tempcoord(l)%am = 0
             tempcoord(l)%bm = 0
             tempcoord(l)%cm = 0
             tempcoord(l)%dm = 0
             tempcoord(l)%em = 0
             tempcoord(l)%fm = 0
             tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y
             tempcoord(l)%z = protcoords(m,l)%z
             tbo(l) = bonddd(m,l)
          end do
          do l = 1,maxl,1
             do f = 1,maxl
                if(abs(f-l)>2) then
                   call adjacent(m,l,f,tempcoord,adjver,bdir)
                   if((adjver .eqv. .true.)) then
                      call bondformintra(m,l,f,tempcoord,tbo,bdir,initialenergy)
                   end if
                end if
             end do
          end do

          !write(6,*) 'm',m
          do g = 1,nprotein,1
             if(g/=m) then
                maxlengthss = chlen(g)
                do l = 1,maxl,1
                   do f = 1,maxlengthss,1
                      call adjacent(g,l,f,tempcoord,adjver,bdir)
                      if((adjver .eqv. .true.)) then
                         olderenergy = initialenergy
                         call bondforminter(m,g,l,f,tempcoord,tbo,bdir,initialenergy)
                         if(initialenergy < olderenergy) then
                            visit(m,l) = 1
                            visit(g,f) = 1
                            end if
                      end if
                   end do
                end do
             end if
          end do
          call updatebondold(m,1,maxl,tempcoord)
       end do
    end if

    if(init .eqv. .false.) then
       do m=1,nprotein,1
chen(m)= 0.0d0
       end do
       do m =1,nprotein,1
          maxl = chlen(m)
          do l =1,maxl,1
             tempcoord(l)%am = 0
             tempcoord(l)%bm = 0
             tempcoord(l)%cm = 0
             tempcoord(l)%dm = 0
             tempcoord(l)%em = 0
             tempcoord(l)%fm = 0
             tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y
             tempcoord(l)%z = protcoords(m,l)%z
             tbo(l) = bonddd(m,l)
          end do
          do l = 1,maxl-3,1
             do f = l+3,maxl
                call adjacent(m,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.)) then
                   dumoldenergy = initialenergy
                   call bondformintra(m,l,f,tempcoord,tbo,bdir,initialenergy)
                   chen(m) = chen(m) + (initialenergy-dumoldenergy)
                end if
             end do
          end do

          if(m<nprotein) then
             do g = m+1,nprotein,1
                maxlengthss = chlen(g)
                do l = 1,maxl,1
                   do f = 1,maxlengthss,1
                      call adjacent(g,l,f,tempcoord,adjver,bdir)
                      if((adjver .eqv. .true.)) then
                         dumoldenergy = initialenergy
                         call bondforminter(m,g,l,f,tempcoord,tbo,bdir,initialenergy)
                         chen(m) = chen(m) + (initialenergy-dumoldenergy)
                      end if
                   end do
                end do
             end do
          end if
       end do
    end if


    if(init .eqv. .true.) initialenergy = initialenergy/2
    if((time > 0) .and. (nint(initialenergy) /= nint(totalenergy))) then
       write(6,*) 'energyfail',totalenergy,initialenergy,time
       finalfail = .true.
    end if

    if(init .eqv. .true.) totalenergy = initialenergy

!write(6,*) 'ENERGYSUBROUTINE',totalenergy
  end subroutine energy








  subroutine adjacent(pr,ch1,ch2,tempory,adjver,bdir)
    integer,intent(in):: pr,ch1,ch2
    integer,intent(inout)::bdir
    integer::dx,dy,dz,delx,dely,delz,signs
    logical,intent(inout)::adjver
    Type(prottemp),dimension(:),intent(in) :: tempory
    dx = 4
    dy = 4
    dz = 4

    delx = tempory(ch1)%x - protcoords(pr,ch2)%x
    if((abs(delx) == 1) .or. (abs(delx) == (gridsizex - 1))) then
       dx = 1
       bdir = delx
if(abs(delx) == (gridsizex - 1)) bdir = (delx/(gridsizex-1))*(gridsize-1)
    else if(delx == 0) then
       dx = 0
    end if
    if(dx == 4) then
       adjver = .false.
       goto 31
    end if

    dely = tempory(ch1)%y - protcoords(pr,ch2)%y
    if((abs(dely) == 1) .or.(abs(dely) == (gridsize - 1))) then
       dy = 1
       bdir = 2*dely
    else if(dely == 0) then
       dy = 0
    end if
    if(dy == 4)then
       adjver = .false.
       goto 31
    end if

    delz = tempory(ch1)%z - protcoords(pr,ch2)%z
    if((abs(delz) == 1) .or. (abs(delz) == (gridsize -1))) then
       dz = 1
       bdir = 3*delz
    else if(delz == 0) then
       dz = 0
    end if
    if(dz == 4) then
       adjver = .false.
       goto 31
    end if

31  if(dx + dy + dz ==1) then
       adjver = .true.
       if(abs(bdir) > 3) then
          signs = int(bdir/(gridsize-1))
          bdir = -1*signs
       end if
       !write(6,*) bdir
    else if(dx + dy + dz /= 1) then
       adjver = .false.
    end if

  end subroutine adjacent


  subroutine oldadjacent(pr,ch1,ch2,tempory,adjver,bdir)
    integer,intent(in):: pr,ch1,ch2
    integer,intent(inout)::bdir
    integer::dx,dy,dz,delx,dely,delz,signs
    logical,intent(inout)::adjver
    Type(prottemp),dimension(:),intent(in) :: tempory
    dx = 4
    dy = 4
    dz = 4

    delx = tempory(ch1)%x - protcoords(pr,ch2)%x
    if((abs(delx) == 1) .or. (abs(delx) == (gridsize - 1))) then
       dx = 1
       bdir = delx
if(abs(delx) == (gridsize - 1)) bdir = (delx/(gridsize-1))*(gridsizex-1)
    else if(delx == 0) then
       dx = 0
    end if
    if(dx == 4) then
       adjver = .false.
       goto 31
    end if

    dely = tempory(ch1)%y - protcoords(pr,ch2)%y
    if((abs(dely) == 1) .or.(abs(dely) == (gridsize - 1))) then
       dy = 1
       bdir = 2*dely
    else if(dely == 0) then
       dy = 0
    end if
    if(dy == 4)then
       adjver = .false.
       goto 31
    end if

    delz = tempory(ch1)%z - protcoords(pr,ch2)%z
    if((abs(delz) == 1) .or. (abs(delz) == (gridsize -1))) then
       dz = 1
       bdir = 3*delz
    else if(delz == 0) then
       dz = 0
    end if
    if(dz == 4) then
       adjver = .false.
       goto 31
    end if

31  if(dx + dy + dz ==1) then
       adjver = .true.
       if(abs(bdir) > 3) then
          signs = int(bdir/(gridsize-1))
          bdir = -1*signs
       end if
       !write(6,*) bdir
    else if(dx + dy + dz /= 1) then
       adjver = .false.
    end if

  end subroutine oldadjacent


  subroutine comfind(m,track)
    integer :: dcomx,dcomy,dcomz,maxl,comx,comy,comz
    integer :: l,a
    logical,intent(in) :: track
    integer,intent(in)::m


    if(time >0) then
       com(2,m)%x = com(1,m)%x
       com(2,m)%y = com(1,m)%y
       com(2,m)%z = com(1,m)%z
    end if

       comx = 0
       comy = 0
       comz = 0
        maxl = chlen(m)
!linkerlength = 7

        !do a =1,maxl
           !write(6,*) 'bead',a,protcoords(m,a)%x,protcoords(m,a)%y,protcoords(m,a)%z
        !end do
        do l = 2,maxl,1
           dcomx = protcoords(m,l)%x-protcoords(m,l-1)%x
           if (dcomx >protcoords(m,1)%linker) dcomx = dcomx-gridsizex
           if (dcomx <(-protcoords(m,1)%linker)) dcomx = gridsizex+dcomx
           comx = comx + ((maxl+1-l)*dcomx)
           dcomy = protcoords(m,l)%y-protcoords(m,l-1)%y
           !if(abs(dcomy) > 7) write(6,*) 'over pbc',gridsize,linkerlength,dcomy
           if (dcomy >protcoords(m,1)%linker) dcomy = dcomy-gridsize
           if (dcomy <(-protcoords(m,1)%linker)) dcomy = gridsize+dcomy
           !if(abs(dcomy) > 7) write(6,*) 'faillllll',gridsize,linkerlength,dcomy
           comy = comy + ((maxl+1-l)*dcomy)
           dcomz = protcoords(m,l)%z-protcoords(m,l-1)%z
           if (dcomz >protcoords(m,1)%linker) dcomz = dcomz-gridsize
           if (dcomz <(-protcoords(m,1)%linker)) dcomz = gridsize+dcomz
           comz = comz + ((maxl+1-l)*dcomz)
        end do
        !write(6,*) 'deltas',comx,comy,comz
        com(1,m)%x = modulo(protcoords(m,1)%x +(real(comx)/(maxl)) -1,real(gridsizex))+1
        com(1,m)%y = modulo(protcoords(m,1)%y +(real(comy)/(maxl)) -1,real(gridsize))+1
        com(1,m)%z = modulo(protcoords(m,1)%z +(real(comz)/(maxl)) -1,real(gridsize))+1
        !end do
        !if((com(1,1)%y>1500) .and.(com(1,1)%y<2500))write(6,*) 'com',time,com(1,1)%y,protcoords(1,1)%y,comy
        !write(6,*) 'commmm 2'
        if(time > 0 .and. (track .eqv. .true.)) then
           call tracking(m)
        end if
        
        if(track .eqv. .true.) call rms(m)
    
  end subroutine comfind

  subroutine tracking(m)
    integer,intent(in)::m
    double precision :: dcx,dcy,dcz

          dcx = 0.0d0
          dcy = 0.0d0
          dcz = 0.0d0
         dcx = com(1,m)%x-com(2,m)%x
          if (dcx > gridsizex/2) dcx = dcx-gridsizex
          if (dcx < -gridsizex/2) dcx = gridsizex + dcx
          trackx = trackx + dcx
          dcy = com(1,m)%y-com(2,m)%y
          if (dcy > gridsize/2) dcy = dcy - gridsize
          if (dcy < -gridsize/2) dcy = dcy +gridsize
          tracky = tracky + dcy
          dcz = com(1,m)%z-com(2,m)%z
          if (dcz > gridsize/2) dcz = dcz - gridsize
          if (dcz < -gridsize/2) dcz = dcz + gridsize
          trackz = trackz + dcz
    
  end subroutine tracking


  
  subroutine rms(m)
    double precision :: msdsum,msd,msdx,msdy,msdz
    integer,intent(in)::m

    if(time<= equilib) return 

    msdx = (min(abs(com(1,m)%x - com(2,m)%x), &
         (gridsizex - abs(com(1,m)%x -com(2,m)%x))))**2
    msdy = (min(abs(com(1,m)%y -com(2,m)%y), &
         (gridsize - abs(com(1,m)%y -com(2,m)%y))))**2
    msdz = (min(abs(com(1,m)%z - com(2,m)%z), &
         (gridsize - abs(com(1,m)%z -com(2,m)%z))))**2

   
    totrmsbrute = totrmsbrute + msdx + msdy + msdz

  end subroutine rms

subroutine radiusofgyration
    integer::m,l,maxl,nobeads
    double precision :: rog,totrog!,avex,avey,avez
totrog = 0.0d0
nobeads = 0
!avex = 0.0d0
!avey = 0.0d0
!avez = 0.0d0
    do m = nprotein-nclients,nprotein,1
        maxl = chlen(m)
!do l = 1,maxl,1
!avex = avex + protcoords(m,l)%x
!avey = avey + protcoords(m,l)%y
!avez = avez + protcoords(m,l)%z
!end do
!avex = avex/maxl
!avey = avey/maxl
!avez = avez/maxl

       do l = 1, maxl,1
          rog = (min(modulo(protcoords(m,l)%x-com(1,m)%x,gridsizex*1.0),(gridsizex - &
               modulo(protcoords(m,l)%x-com(1,m)%x,gridsizex*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0))))**2
          totrog = totrog + rog
       end do
       nobeads = nobeads+maxl
    end do

    runningaveROG = runningaveROG + totrog
averageROG = averageROG + totrog       
 noROGS = noROGS+1
    polymerrog = polymerrog + (totrog/nprotein)
    !if(noinfo .eqv. .false.) write(97,*) time-equilib, &
      !   totrog/(nobeads),SQRT(totrog/(nobeads))

  end subroutine radiusofgyration
  
  

  subroutine debug
    integer::m,l,f,g,dx,dy,dz,maxlengthss,maxl,zz
 

    !write(6,*) 'debug start'
    do m = 1,nprotein,1
       maxl = chlen(m)
       !write(6,*) 'max', maxl
       do l = 1,maxl,1
          if (protcoords(m,l)%x == 0) write(6,*) 'x fail',l
          if (protcoords(m,l)%y == 0) write(6,*) 'y fail',l
          if (protcoords(m,l)%z == 0) write(6,*) 'z fail',l


          do f = 1,nprotein,1
             maxlengthss = chlen(f)
             do g = 1,maxlengthss,1
                if(m ==f .and. l == g) then
                   continue
                else if ((m/= f) .or. (l/= g)) then
                   !write(6,*) protcoords(m,l)%x, protcoords(f,g)%x
                   if((protcoords(m,l)%x == protcoords(f,g)%x)  .and. (protcoords(m,l)%y == protcoords(f,g)%y ) &
                        .and. (protcoords(m,l)%z == protcoords(f,g)%z)) then
                      write(6,*) protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z
                      write(6,*) protcoords(f,g)%x,protcoords(f,g)%y,protcoords(f,g)%z
                      write(6,*) 'fail is due to', m,l,f,g
                      finalfail = .true.
                      fail = .true.
                      write(6,*) 'FAILLLLLLLLLLLLLLLLLL',time
                      exit
                   end if
                end if
             end do
          end do

        if(l>1) then

if(m<=nscaffold) then
 dx = min(abs(protcoords(m,l)%x - protcoords(m,l-1)%x),gridsize-abs(protcoords(m,l)%x &
                  - protcoords(m,l-1)%x))
             dy = min(abs(protcoords(m,l)%y - protcoords(m,l-1)%y),gridsize-abs(protcoords(m,l)%y &
                  - protcoords(m,l-1)%y)) 
             dz = min(abs(protcoords(m,l)%z - protcoords(m,l-1)%z),gridsize-abs(protcoords(m,l)%z &
                  - protcoords(m,l-1)%z)) 

              sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
             if(sumdebug > protcoords(m,1)%linker) then
                write(6,*) 'SPLIT',m,l,l-1,sumdebug
                finalfail = .true.
             end if
          end if

if(m>nscaffold) then
 dx = min(abs(protcoords(m,l)%x - protcoords(m,l-1)%x),gridsizex-abs(protcoords(m,l)%x &
                  - protcoords(m,l-1)%x))
             dy = min(abs(protcoords(m,l)%y - protcoords(m,l-1)%y),gridsize-abs(protcoords(m,l)%y &
                  - protcoords(m,l-1)%y))
             dz = min(abs(protcoords(m,l)%z - protcoords(m,l-1)%z),gridsize-abs(protcoords(m,l)%z &
                  - protcoords(m,l-1)%z))

              sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
             if(sumdebug > protcoords(m,1)%linker) then
                write(6,*) 'SPLIT',m,l,l-1,sumdebug
                finalfail = .true.
             end if
          end if
end if


          if(bonddd(m,l) == 0) write(6,*) 'bond vector fail'
       end do
    end do

   

    
         
    do m =1,nprotein
       maxl = chlen(m)
       do l=1,maxl


          if(protcoords(m,l)%am /= 0) then
             g = protcoords(m,l)%am
             f = protcoords(m,l)%al
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'a'
                finalfail = .true.
             end if
             if((protcoords(g,f)%bm /=m) .or. (protcoords(g,f)%bl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'a',protcoords(g,f)%bm,&
                     protcoords(g,f)%bl
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
          end if

          if(protcoords(m,l)%bm /= 0) then
             g = protcoords(m,l)%bm
             f = protcoords(m,l)%bl
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'b'
                finalfail = .true.
             end if
             if((protcoords(g,f)%am /=m) .or. (protcoords(g,f)%al /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'b',protcoords(g,f)%am,&
                     protcoords(g,f)%al
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
          end if


          if(protcoords(m,l)%cm /= 0) then
             g = protcoords(m,l)%cm
             f = protcoords(m,l)%cl
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'c'
                finalfail = .true.
             end if
             if((protcoords(g,f)%dm /=m) .or. (protcoords(g,f)%dl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'c',protcoords(g,f)%dm,&
                     protcoords(g,f)%dl
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
          end if

          if(protcoords(m,l)%dm /= 0) then
             g = protcoords(m,l)%dm
             f = protcoords(m,l)%dl
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'d'
                finalfail = .true.
             end if
             if((protcoords(g,f)%cm /=m) .or. (protcoords(g,f)%cl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'d',protcoords(g,f)%cm,&
                     protcoords(g,f)%cl
                finalfail = .true.
             end if
          end if

          if(protcoords(m,l)%em /= 0) then
             g = protcoords(m,l)%em
             f = protcoords(m,l)%el
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'e'
                finalfail = .true.
             end if
             if((protcoords(g,f)%fm /=m) .or. (protcoords(g,f)%fl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'e',protcoords(g,f)%fm,&
                     protcoords(g,f)%fl
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
          end if

          if(protcoords(m,l)%fm /= 0) then
             g = protcoords(m,l)%fm
             f = protcoords(m,l)%fl
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'f'
                finalfail = .true.
             end if
             if((protcoords(g,f)%em /=m) .or. (protcoords(g,f)%el /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'f',protcoords(g,f)%em,&
                     protcoords(g,f)%el
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
          end if

       end do
    end do


    if(finalfail .eqv. .true.) then
       write(6,*) time
       
       stop 9
    end if

  end subroutine debug







  subroutine freesitescalc
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn
    integer,dimension(:),allocatable:: countoffree
    logical :: clusyes,adjver,newcl
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize,maxl
    integer,dimension(:,:),allocatable::cb,bound
    integer,dimension(:),allocatable::maxcluslist,concount
    double precision::maxclusenergy,tryen
    allocate(cb(nprotein,nprotein))
    !allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(bound(nprotein,maxlength))
    allocate(concount(nprotein))
allocate(countoffree(2*mtype))
    !allocate(mashist(10))
    maxclusenergy = 0.0d0
tryen = 0.0d0
    normalisedsize = 0
    clcount = 0
    clusterpop = 0
    clustcount = 0
    do m = 1,nprotein
       clnos(m) = m
       cllist(m) = 0
       histcl(m) = 0
       concount(m) =0
       do a =1,maxlength
          bound(m,a)=0
          end do
           end do

           a = 1

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do




    do m = 1,nprotein-1,1
       maxback = chlen(m)
       do g = m+1,nprotein,1
          !if((g/=m) .and. (clnos(m) /= clnos(g))) then
             maxlengthss = chlen(g)
             do l = 1,maxback
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z
                do f = 1,maxlengthss
                   call oldadjacent(g,l,f,tempcoord,adjver,bdir)              
                   if(adjver.eqv. .true.) then
                      if(isbond .eqv. .true.) then
                         cb(m,g) = 1
                         cb(g,m) = 1
bound(m,l) = 1
bound(g,f) = 1                         
if(interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0) then

                                 
                            if(clnos(g) < clnos(m)) then
                               oldcl = clnos(m)
                               do zzz = 1,nprotein
                                  if(clnos(zzz) == oldcl) clnos(zzz) = clnos(g)
                               end do
                            else if(clnos(g)> clnos(m)) then
                               oldcl =clnos(g)
                               do zzz = 1,nprotein
                                  if(clnos(zzz) == oldcl) clnos(zzz) = clnos(m)
                               end do
                            end if
                            clusyes = .true.
                         end if
                      end if
                      if((bdir == bonddd(g,f)) .and. (-1.0*bdir == (bonddd(m,l))) .and. &
                           (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0)) then
tryen = tryen + (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)) 
                      bound(m,l) = 1           
                        bound(g,f) = 1
                             
                         cb(m,g) = 1
                         cb(g,m) = 1
                         if(clnos(g) < clnos(m)) then
                            oldcl = clnos(m)
                            do zzz = 1,nprotein
                               if(clnos(zzz) == oldcl) clnos(zzz) = clnos(g)
                            end do
                         else if(clnos(g)> clnos(m)) then
                            oldcl =clnos(g)
                            do zzz = 1,nprotein
                               if(clnos(zzz) == oldcl) clnos(zzz) = clnos(m)
                            end do
                         end if
                         clusyes = .true.                 
                      end if
                   end if
                end do
             end do
          !end if
       end do
    end do






    write(6,*) 'initial energy',tryen


call energy(.false.)

    do g = 1,nprotein-nclients
       do m = g,nprotein-nclients
          if(clnos(m) ==g) then
             histcl(g) = histcl(g) +1
          end if
       end do
       !if(time>10000) write(6,*) 'hist', histcl(g)
       normalisedsize =  normalisedsize + ((histcl(g))**2)
    end do


    !maxclus = maxval(histcl)


    maxclus = 0
    do g = 1,nprotein-nclients
       if(histcl(g) > maxclus) then
          maxclus = histcl(g)
          content = g
       end if
    end do

countoffree(:) = 0

write(6,*) 'freesites COM',maxclus,content,dpcomclusterinit
    do g = 1,nprotein-nclients,1
       if(protcoords(g,1)%x<=gridsize) then
maxl = chlen(g)
do l = 1,maxl,1
if(bound(g,l) == 0) then
call leftdistribution(g,l)
countoffree(protcoords(g,l)%type) = countoffree(protcoords(g,l)%type) + 1
end if
end do       
else if(protcoords(g,1)%x > gridsize) then
maxl = chlen(g)
do l = 1,maxl,1      
if(bound(g,l) == 0) then
  call rightdistribution(g,l)
countoffree(protcoords(g,l)%type+2) = countoffree(protcoords(g,l)%type+2) + 1
end if
end do
       end if
    end do
    
    write(6,*) 'unbound species in left and right',countoffree(:)

do l= 1,2*mtype
totalfree(l) = totalfree(l) +countoffree(l)
end do
  
   
    magiccontent = content


call neighbour(bound)

  
    deallocate(cllist)
    !deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
    !write(6,*) 'CONTENTTTTTTTTT',content
  end subroutine freesitescalc


subroutine leftdistribution(m,l)

   integer,intent(in)::m,l
    integer :: maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


          normaliser = normaliser + 1
          !phasecount = phasecount+1
          dx = min(abs(protcoords(m,l)%x - (gridsize/2.0)), gridsizex - abs(protcoords(m,l)%x - (gridsize/2.0)))
          dy = min(abs(protcoords(m,l)%y - (gridsize/2.0)), gridsize - abs(protcoords(m,l)%y - (gridsize/2.0)))
          dz = min(abs(protcoords(m,l)%z - (gridsize/2.0)), gridsize - abs(protcoords(m,l)%z - (gridsize/2.0)))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))

!write(6,*) 'type', protcoords(m,l)%type
   freesites(protcoords(m,l)%type,lx+1) = freesites(protcoords(m,l)%type,lx+1)+1.0
   freesites(protcoords(m,l)%type+mtype,lx+1) = freesites(protcoords(m,l)%type+mtype,lx+1)+1.0

    
     end subroutine leftdistribution

     subroutine rightdistribution(m,l)

   integer,intent(in)::m,l
    integer :: maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


!phasecountdoub = phasecountdoub+1

          normaliser = normaliser + 1
          !phasecount = phasecount+1
          dx = min(abs(protcoords(m,l)%x - (gridsize*(3.0/4.0))), gridsize - abs(protcoords(m,l)%x - (gridsize*(3.0/4.0))))
          dy = min(abs(protcoords(m,l)%y - (gridsize/2.0)), gridsize - abs(protcoords(m,l)%y - (gridsize/2.0)))
          dz = min(abs(protcoords(m,l)%z - (gridsize/2.0)), gridsize - abs(protcoords(m,l)%z - (gridsize/2.0)))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


   freesites(protcoords(m,l)%type,lx+1) = freesites(protcoords(m,l)%type,lx+1)+1.0
   freesites(protcoords(m,l)%type+(2*mtype),lx+1) = freesites(protcoords(m,l)%type+(2*mtype),lx+1)+1.0

    
     end subroutine rightdistribution

  
  SUBROUTINE read_setup

    USE keywords
    !USE setup_vars
    !USE stats_vars
    IMPLICIT NONE

    CHARACTER(LEN=20) :: keyword, option, argument
    CHARACTER(LEN=18), PARAMETER :: param_fmt1='(A16, 1PE20.10, A)'
    CHARACTER(LEN=16), PARAMETER :: param_fmt2='(A16, F16.10)'
    INTEGER :: err, i, j,runtype,dummy2,dummytype,xl,f,dum,maxlengthstart
    LOGICAL :: success
    double precision:: intra


    !  Defaults and descriptions.
    !  S.T.A.T. stands for short-time averaged temperature.




    seed = 6                      !seed for random number generation
    maxlength = 20              !number of beads in chain
    nprotein = 1                 !number of protcoords
    gridsize= 100        !dimensions of the lattice
    maxtime= 10000                !length of simulation
    equilib = 5000              !equilibration time
    rept = 1                  !1 = reptation 0 = no reptations
    kT = 1.0
    film = .true.
    debugyes = .true.
nclientspecies = 1
nclients = 0

    OPEN (UNIT=20, FILE='2systemsetupclient.txt', STATUS='OLD', IOSTAT=err)

    IF (err /= 0) THEN
       WRITE (6, '(/, A, /)') 'ERROR: Could not find setup file.'
       STOP
    ENDIF

    DO
       CALL read_line(20, success)
       IF (.NOT.success) EXIT
       CALL get_string(keyword)
       CALL upper_case(keyword)

       IF (keyword(1:1) == '#') keyword='#'

       SELECT CASE (TRIM(keyword))

       CASE ('#')
          !        Comment line; do nothing.

       CASE ('SEED')
          CALL get_integer(seed)
          seed = -ABS(seed)

       CASE ('TWOSYS')
          call get_logical(twosys)

          
          !CASE ('CHAIN_LENGTH')
          !CALL get_integer(maxlength1)

          !CASE('CHAIN_2_LENGTH')
          !call get_integer(maxlength2)

          !CASE ('NUM_CHAINS')
          !CALL get_integer(nprotein1)

          !CASE ('NUM_CHAINS_2')
          !CALL get_integer(nprotein2)
       CASE ('SCALEINFO')
          call get_logical(scalinginfo)

       CASE ('TOTTYPES')
          CALL get_integer(nspecies)
          call get_integer(maxlengthstart)
          allocate (speciespop(nspecies))
          allocate(specieslen(nspecies))
          allocate(specieslinker(nspecies))

       CASE ('TOTCHAINS')
          !CALL get_integer(nprotein)

          runtype = 1
       CASE ('NEWCHAIN')
          CALL get_integer(specieslen(runtype))
          CALL get_integer(speciespop(runtype))
          call get_integer(specieslinker(runtype))   
          do f = 1,specieslen(runtype)
             call get_integer(dummytype)
             if((dummytype>mtype)) mtype = dummytype
          end do
          write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
          runtype = runtype + 1
                vtype = vtype+1
        nclients = 0
          nclientspecies = 1 !nspecies- mtype
       
deallocate(client)
allocate(client(nspecies,maxlengthstart))

       CASE ('CLIENT')
          mtype = mtype +1
        vtype = vtype+1
          CALL get_integer(specieslen(vtype))
          CALL get_integer(speciespop(vtype))
          clientlength = specieslen(vtype)
          client(vtype,1)%species = vtype
          call get_integer(client(vtype,1)%linker)   
          do f = 1,specieslen(vtype)
             call get_integer(client(vtype,f)%type)
             client(vtype,f)%linker = client(vtype,1)%linker
             write(6,*) 'client details',client(vtype,f)%type
                          !if(dummytype>mtype) mtype = dummytype
          end do
          nclients = nclients+speciespop(vtype)

          write(6,*) 'runtype',vtype,specieslen(vtype),speciespop(vtype)
          
          

       CASE ('LATTICE_DIMENSIONS')
          CALL get_integer(gridsize)


          allocate(interen(mtype-nclientspecies,mtype-nclientspecies))
          allocate(interenergy(mtype-nclientspecies,mtype-nclientspecies))
          allocate(intraenergy(mtype-nclientspecies,mtype-nclientspecies))
          !allocate(tttt())

       CASE ('SIM_TIME')
          CALL get_integer(maxtime)
          maxtime = maxtime*1000
       CASE ('EQUIB_TIME')
          CALL get_integer(equilib)
          equilib = equilib*1000
       CASE ('RIGHTANGLE')
          !CALL get_integer(right)

       CASE ('CRANK')
          !CALL get_integer(crank)

       CASE ('REPTATION')
          CALL get_integer(rept)

       CASE ('PIVOT')
          !CALL get_integer(piv)

       CASE ('KT')
          CALL get_dp(kT)

       CASE ('FILM')
          CALL get_logical(film)

       CASE ('DEBUGYES')
          CALL get_logical(debugyes)


          
       CASE ('FREEZE')
          CALL get_integer(freezetime)
          
          !CASE ('ISBOND')
          !CALL get_logical(isbond)

       CASE ('MAXPIV')
          !CALL get_integer(maxpiv)

       CASE ('CRANKLIMIT')
          !CALL get_integer(cranklimit)



       CASE ('INTRAEN')
          CALL get_dp(intraen)
          intraen = -ABS(intraen)


       CASE ('INTERACTION')
          call get_integer(xl)
          call get_dp(intra)
          write(6,*) 'mtype',mtype
          do f = 1,mtype-nclientspecies
             call get_dp(interenergy(xl,f))
             intraenergy(xl,f) = -abs(intra)
             interenergy(xl,f) = -abs(interenergy(xl,f))
             call get_dp(interen(xl,f))
             interen(xl,f) = -abs(interen(xl,f))
          end do

       CASE ('CLROTATION')
          CALL get_logical(clr)

       CASE ('CLTRANSLATION')
          CALL get_logical(clt)

       CASE ('NOINFO')
          CALL get_logical(noinfo)

          !CASE ('SCALEINFO')
          !CALL get_logical(scalinginfo)

       CASE ('RESTART')
          call get_logical(restart)

       CASE ('LINK')
          !call get_integer(linkerlength)

CASE ('CLUSSTEP')
CALL get_integer(steplength)

       CASE DEFAULT
          WRITE (6, '(/, 3A, /)') 'ERROR: Keyword not recognised in setup file: ', TRIM(keyword), '.'
          STOP

       END SELECT
    END DO
    CLOSE(20)


  END SUBROUTINE read_setup

  !................................................................................!

!!!
!!! Little subroutine to print an error if a keyword is followed by an
!!! unrecognised option.
!!!

  SUBROUTINE bad_option(keyword, option)

    IMPLICIT NONE
    CHARACTER(LEN=*) :: keyword, option

    WRITE (6, '(/, 5A, /)') 'ERROR: ', TRIM(option), &
         & ' is not a valid option after the keyword ', TRIM(keyword), &
         &' in the setup file.'
    STOP

  END SUBROUTINE bad_option



end program move


FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
       NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do 11 j=NTAB+8,1,-1

        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
     endif
     k=idum/IQ1
     idum=IA1*(idum-k*IQ1)-k*IR1
     if (idum.lt.0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
     if (idum2.lt.0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
     if(iy.lt.1)iy=iy+IMM1
     ran2=min(AM*iy,RNMX)
     return
   END Function ran2



   MODULE keywords
     IMPLICIT NONE
     INTEGER, PARAMETER :: max_lines=2, max_length=100
     INTEGER, PARAMETER :: tot_length=max_lines*max_length
     INTEGER :: position
     CHARACTER(LEN=tot_length) :: input
     SAVE
   CONTAINS
     !!
     SUBROUTINE read_line(u, success)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: u
       LOGICAL, INTENT(INOUT), OPTIONAL :: success
       INTEGER :: i, lines, next, err
       LOGICAL :: continue
       CHARACTER(LEN=7) :: fmt_string
       !     Generate format string of max_length characters.
       WRITE (fmt_string, '(I4)') max_length
       fmt_string = '(A' // TRIM(ADJUSTL(fmt_string)) // ')'
       DO
          input=' '
          next = 1
          !        Read in a logical line consisting of up to max_lines of input file.
          DO lines=1, max_lines
             continue = .FALSE.
             READ (u, fmt_string, IOSTAT=err) input(next:next+max_length-1)
             IF (err == 0) THEN
                IF (PRESENT(success)) success = .TRUE.
             ELSE
                IF (PRESENT(success)) success = .FALSE.
                EXIT
             ENDIF
             !           Check for continuation symbol (&).
             DO i=next, next+max_length-1
                IF (input(i:i)=='&') THEN
                   continue = .TRUE.
                   next = i
                   EXIT
                ENDIF
             END DO
             IF (.NOT.continue) EXIT
          END DO
          IF (err /= 0) EXIT
          IF (TRIM(input) /= '') EXIT   ! Only read in next line if this one's empty.
       END DO
       position = 1
     END SUBROUTINE read_line
     !!
     SUBROUTINE upper_case(string)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(INOUT) :: string
       INTEGER, PARAMETER :: lower_to_upper = ICHAR("A")-ICHAR("a")
       INTEGER :: i
       DO i=1, LEN_TRIM(string)
          IF (LGE(string(i:i), 'a').AND.LLE(string(i:i), 'z')) THEN
             string(i:i) = ACHAR(IACHAR(string(i:i))+lower_to_upper)
          ENDIF
       END DO
     END SUBROUTINE upper_case
     !!
     SUBROUTINE get_string(string, success)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(INOUT) :: string
       LOGICAL, INTENT(INOUT), OPTIONAL :: success
       CHARACTER(LEN=tot_length) :: temp
       INTEGER :: outcome
       CALL next_item(temp, outcome)
       IF (outcome == 3) THEN
          string = temp
          IF (PRESENT(success)) success = .TRUE.
       ELSE
          IF (PRESENT(success)) success = .FALSE.
       ENDIF
     END SUBROUTINE get_string
     !!
     SUBROUTINE get_integer(value, success)
       IMPLICIT NONE
       INTEGER, INTENT(INOUT) :: value
       LOGICAL, INTENT(INOUT), OPTIONAL :: success
       INTEGER :: temp, outcome, err
       CHARACTER(LEN=tot_length) :: item
       CALL next_item(item, outcome)
       READ (UNIT=item, FMT=*, IOSTAT=err) temp
       IF ((err == 0).AND.(outcome==3)) THEN
          value = temp
          IF (PRESENT(success)) success = .TRUE.
       ELSE
          IF (PRESENT(success)) success = .FALSE.
       ENDIF
     END SUBROUTINE get_integer
     !!
     SUBROUTINE get_dp(value, success)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(INOUT) :: value
       LOGICAL, INTENT(INOUT), OPTIONAL :: success
       DOUBLE PRECISION :: temp
       INTEGER :: outcome, err
       CHARACTER(LEN=tot_length) :: item
       CALL next_item(item, outcome)
       READ (UNIT=item, FMT=*, IOSTAT=err) temp
       IF ((err == 0).AND.(outcome==3)) THEN
          value = temp
          IF (PRESENT(success)) success = .TRUE.
       ELSE
          IF (PRESENT(success)) success = .FALSE.
       ENDIF
     END SUBROUTINE get_dp
     !!
     SUBROUTINE get_logical(lgcl, success)
       IMPLICIT NONE
       LOGICAL, INTENT(INOUT) :: lgcl
       LOGICAL, INTENT(INOUT), OPTIONAL :: success
       INTEGER :: outcome
       CHARACTER(LEN=tot_length) :: item
       CALL next_item(item, outcome)
       CALL upper_case(item)
       IF ((TRIM(item)=='TRUE').OR.(TRIM(item)=='T').OR.(TRIM(item)=='.TRUE.') &
            & .OR.(TRIM(item)=='.T.').OR.(TRIM(item)=='ON')) THEN
          lgcl = .TRUE.
          IF (PRESENT(success)) success = .TRUE.
       ELSE IF ((TRIM(item)=='FALSE').OR.(TRIM(item)=='F').OR.(TRIM(item)=='.FALSE.') &
            & .OR.(TRIM(item)=='.F.').OR.(TRIM(item)=='OFF')) THEN
          lgcl = .FALSE.
          IF (PRESENT(success)) success = .TRUE.
       ELSE
          IF (PRESENT(success)) success = .FALSE.
       ENDIF
     END SUBROUTINE get_logical
     !!
     SUBROUTINE next_item(item, outcome)
       IMPLICIT NONE
       CHARACTER(LEN=tot_length), INTENT(OUT) :: item
       INTEGER, INTENT(OUT) :: outcome
       !     Values of outcome:
       !      1: null string read
       !      2: end of line reached with no string
       !      3: correctly read a string of at least one character
       INTEGER :: i, j
       item = ''
       outcome=1
       !     Check we've not already reached the end of the input string.
       IF (position > tot_length) THEN
          outcome = 2
       ELSE
          !        Read past leading blanks.
          DO i=position, tot_length
             IF (input(i:i) /= ' ') EXIT
          END DO
          !        Check that this hasn't brought us to the end of the input string.
          IF (i==tot_length+1) THEN
             outcome=2
          ELSE
             position = i
             j = 1
             !           Read until the next space or comma.
             DO i=position, tot_length
                SELECT CASE(input(i:i))
                   !              If we've reached a comma, record the position for the next
                   !              item and exit loop.
                CASE(',')
                   position = i+1
                   EXIT
                   !              If we've reached a space, check for a comma preceded by some
                   !              blanks, and record the position for the next item as after the
                   !              comma if one is found.
                CASE(' ')
                   DO j=position+1, tot_length
                      SELECT CASE(input(j:j))
                      CASE(',')
                         position = j+1
                         EXIT
                      CASE (' ')   ! Do nothing.
                      CASE DEFAULT
                         position = j
                         EXIT
                      END SELECT
                   END DO
                   EXIT
                   !              Any other character is the next character of the item being read.
                CASE DEFAULT
                   item(j:j) = input(i:i)
                   j = j + 1
                   outcome=3
                   position = i+1
                END SELECT
             END DO
          ENDIF
       ENDIF
     END SUBROUTINE next_item
     !!
   END MODULE keywords



