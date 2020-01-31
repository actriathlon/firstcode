program move

  implicit none

  double precision :: totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength,choose2
  double precision :: actualchainlength,random,totiso,sumsep,bulksumsep
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug,testdummy,polymerrog,polymerl
  integer :: rept,datayes,maxclussize,minl,outputrate,backtime,magiccontent,AB,AC,BC,ABC,bAB,bAC,bBC,bABC
  integer ::lastruntime,timebase,steplength,macc,mrej,roatt,roacc,clientlength,telacc,telatt
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count,avecluspop,aveclusnumber,vtype
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib,cltacc,cltatt,binsize,normaliser
  integer :: successful,reject,maxlength1,maxlength2,q,d,u,sm,qt,ft,gt,scans,nclients,phasecount
  integer:: rerej,reacc,nprotein1,nprotein2,nspecies,mtype,freezetime,nclientspecies,nscaffold
  integer,dimension(:,:),allocatable ::study,sumclient,visit,visit2,unbound
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen,clientlist,specieslinker
  integer,dimension(:),allocatable:: clnos,totalbindevents
  integer::count_rate,count_max,start,finish,midpoint
  !integer,dimension(:),allocatable::
  double precision :: trackx,tracky,trackz,intraen,intdummy,totdire,r
  double precision,dimension(:),allocatable::chen
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy
  real, external :: ran2
!character(len = 10) :: commandread  
integer,dimension(:),allocatable::locb,loc
real,dimension(:),allocatable :: variance,disp
  real,dimension(:,:),allocatable:: summing,vdist,freesites,clusdist
  logical :: exist,fail,finalfail,debugyes,film,isbond,clt,clr,noinfo,scalinginfo,nobonds,restart
  real::choose1
  !integer::removethis,pivde
  logical :: suffclust,settingup,twosys
        real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
integer::pointless1,pointless2,dright,dleft,noROGS,a111,a11n
double precision:: averageROG
integer,dimension(:),allocatable::zloc,bindinghist
integer::boundcount,unboundcount
integer,dimension(:),allocatable::unboundarray,boundarray
integer::bondformacccount,bonddeleteacccount,bondformrejcount,bonddeleterejcount
integer,dimension(:,:),allocatable::scaffoldarray
character(len=50)::pathtofile
integer::avesepcount
double precision::avesep,aveete

  !cluster COM may fall foul of pbcs

  !Sort MSD outputs
  
  type protein
     integer ::x,y,z,species,type,linker,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl,a111,a11n
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
  type(lbond),dimension(:,:,:),allocatable::cell
     call system_clock(start,count_rate,count_max)


write(6,*) 'NEW CODE'
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
telatt = 0

boundcount=0
unboundcount=0
pointless1= 0
pointless2=0
dright = 0
dleft = 0
noROGS = 0
averageROG = 0.0d0

allocate(locb(2))
locb(1) = 0
locb(2) = 0
allocate(totalbindevents(2))
totalbindevents(:) = 0
  allocate(client(1,5))
  mtype = 0
vtype=0        
!noinfo = .false.
  call read_setup

 !call get_command_argument(1,commandread)
  !read(commandread, '(i3)') freezetime



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
  open(21, file = pathtofile, action = 'read')

  if(restart .eqv. .false.) then
     !open(23, file = 'initialtake2.xyz', action = 'read')
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',action = 'write')
        open(82,file = 'newclusterdata.dat', action = 'write')
     else if(noinfo .eqv. .false.) then
        if(scalinginfo .eqv. .true.) then
           open(11,file='scalingdata.dat',action = 'write')
           open(29, file = 'rms.dat', action = 'write')
           open(91, file = 'avechainlength.dat', action = 'write')
           open(97, file = 'radiusofgyration.dat', action = 'write')
           open(79, file = 'runningave.dat', action = 'write')
           open(93, file = 'energy.dat', action = 'write')
           open(13,file='timedata.dat',action = 'write')
           open(19,file='acceptance.dat',action  = 'write')
        else if(scalinginfo .eqv. .false.) then
           open(29, file = 'rms.dat', action = 'write')
           open(93, file = 'energy.dat', action = 'write')
           open(77,file='phasetrans.dat',action = 'write')
           open(82,file = 'newclusterdata.dat', action = 'write')
           open(13,file='timedata.dat',action = 'write')
           open(19,file='acceptance.dat',action  = 'write')
           !open(38,file = 'histclust.dat',action = 'write')
           open(3,file='coms.dat',action = 'write')
           !open(109,file='bondcount.dat',action='write')
           !open(67,file='bondddd.dat',action  = 'write')
        end if
     end if
  else if(restart .eqv. .true.) then
     open(23, file = 'checkpoint.xyz', action = 'read')
     write(6,*) 'logicals',scalinginfo,noinfo
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',access = 'append')
        open(82,file = 'newclusterdata.dat', access = 'append')
     else if(noinfo .eqv. .false.) then
        if(scalinginfo .eqv. .true.) then
           open(11,file='scalingdata.dat',access = 'append')
           open(29, file = 'rms.dat', access = 'append')
           open(91, file = 'avechainlength.dat', access = 'append')
           open(97, file = 'radiusofgyration.dat', access = 'append')
           open(79, file = 'runningave.dat', access = 'append')
           open(93, file = 'energy.dat', access = 'append')
           open(13,file='timedata.dat',access = 'append')
           open(19,file='acceptance.dat',access  = 'append')
        else if(scalinginfo .eqv. .false.) then
           open(29, file = 'rms.dat', access = 'append')
           open(93, file = 'energy.dat', access = 'append')
           open(82,file = 'newclusterdata.dat', access = 'append')
           open(13,file='timedata.dat',access = 'append')
           open(77,file='phasetrans.dat',access = 'append')
           open(19,file='acceptance.dat',access  = 'append')
           !open(38,file = 'histclust.dat',access = 'append')
           open(3,file='coms.dat',action = 'write')
           !open(109,file='bondcount.dat',action='write')
           !open(67,file='bondddd.dat',action  = 'write')
        end if
     end if
  end if


  
  reptforward = 0
  reptbackward = 0
  totrmsbrute = 0.0d0
  totalenergy = 0.0d0
  runningaveROG = 0.0d0
  runningaveEtE = 0.0d0
  avecluspop = 0
  aveclusnumber = 0
  finalfail = .false.
  maxclussize = 1
  isbond = .true.
  trackx = 0.0d0
  tracky = 0.0d0
  trackz = 0.0d0
  polymerl = 0.0d0
  polymerrog = 0.0d0
avesepcount=0
avesep=0.0d0
aveete=0.0d0

  if(restart .eqv. .true.) equilib = 0
  maxtime = maxtime + equilib
  outputrate = (maxtime-equilib)/10
  
  cltacc = 0
  cltatt = 0
  !kT = 50.0
  successful = 0
  reject = 0
  backtime = maxtime


  totdire = 0.0
  totiso = 0.0d0
  do ft = 1,mtype!-nclientspecies
     do gt = 1,mtype!-nclientspecies
        totdire = totdire + interenergy(ft,gt)
        totiso = totiso + interen(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
  end do
  if(int(totiso) == 0) isbond = .false.
  if(int(totdire) == 0) nobonds =  .true.

  !interenergy = -1.0

  !write(6,*) maxlength1,maxlength2
  write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT
  write(6,*) 'clr',clr,'clt',clt

  if(restart .eqv. .false.) then
     open(67, file = 'move.vtf', action = 'write')
open(88,file = 'newmovetagged.vtf',action='write')
  else if(restart .eqv. .true.) then
     open(67, file = 'move.vtf', access = 'append')
     open(88, file = 'newmovetagged.vtf', access = 'append')
  end if


open(141,file='linkerhistogram.dat',action='write')
open(160,file='origdimens.dat',action='write')
open(161,file='dimens.dat',action='write')


  allocate(protcoords(nprotein,maxlength))
  !allocate(protcoords(nprotein,maxlength))
  allocate(com(2,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(chlen(nprotein))
  allocate(chen(nprotein))
  allocate(clientlist(nclients))
  allocate(visit(nprotein,maxlength))
  allocate(clnos(nprotein))
  allocate(vdist(4*clientlength,gridsize))
allocate(clusdist(2*mtype,gridsize))
allocate(loc(2))
allocate(cell(gridsize,gridsize,gridsize))
allocate(bindinghist(clientlength+1))
allocate(zloc(gridsize))
allocate(unboundarray(clientlength*nclients*6))
allocate(boundarray(clientlength*nclients*6))
!allocate(supercala(4,clientlength))
!supercala(:,:) = 0

zloc(:) = 0
bindinghist(:) = 0

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

reacc = 0
rerej = 0
  macc = 0
  mrej = 0
  roatt = 0
  roacc = 0
  bondformacccount=0
  bondformrejcount=0
  bonddeleterejcount=0
  bonddeleteacccount=0
nscaffold = nprotein-nclients

write(6,*) 'NSCAFFOLD',nscaffold

allocate(study(nclients,maxlength))
allocate(scaffoldarray(nscaffold,maxlength))
scaffoldarray(:,:)=0

  !sets the limit of the pivot move


  !call CPU_TIME(start)




  write(6,*) 'foundation start'
  if(restart .eqv. .false.) then
     
	if(vtype>1) call foundation
     if(twosys .eqv. .true.) call foundationtwo
     call clientfoundation
call debug
  else if (restart .eqv. .true.) then
     call foundationrestart
     call restdata
  end if
  write(6,*) 'foundation complete'

  if((scalinginfo .eqv. .false.) .and. (restart .eqv. .false.)) call pdbsetup
  minl = minval(specieslen)

call initialisearray

  allocate(summing(2,(protcoords(nprotein,1)%linker*binsize)+1))
  summing(:,:) = 0
  allocate(freesites(mtype*6,gridsize*binsize))
  freesites(:,:) = 0

  write(6,*) 'a'
  !call dataout
  do qt = 1,nprotein
     call comfind(qt,.false.)
     com(2,qt)%x = com(1,qt)%x
     com(2,qt)%y = com(1,qt)%y
     com(2,qt)%z = com(1,qt)%z
  end do
  write(6,*) 'b'


  settingup = .true.

             dpcomclusterinit%x = real(gridsize)/2
             dpcomclusterinit%y = real(gridsize)/2
             dpcomclusterinit%z = real(gridsize)/2

  !call clustercount
  !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',content
settingup = .false.
normaliser = 0
!call freesitescalc


write(6,*) 'centre of mass details:',dpcomclusterinit

  
  actualchainlength = 0.0
  do timebase = 1,maxtime,1


!if(timebase == 3) call freesitescalc
     if(restart .eqv. .false.) then
        time = timebase
     else
        time = timebase+lastruntime
     end if
     if(time < 0) time = abs(time)
     fail = .false.

     if ((mod(time,outputrate) == 0) .and. (film .eqv. .true.) ) then
        call dataout
     end if
!write(6,*) 'c'

      call pickmoves

if(modulo(time,100)==0) then
      call analyse
        call sepcalc
end if
     
     if (mod(time,outputrate) == 0) then
        if(noinfo .eqv. .false.)  write(93,*) time,totalenergy
     end if
     if(modulo(time,outputrate) == 0 .and. (scalinginfo .eqv. .false.)) then 
        call energy(.false.)
        !call clustercount
     end if

     if((time > equilib) .and. (modulo(time,outputrate) == 0) .and. (scalinginfo .eqv. .true.)) then
        call length
        !call radiusofgyration
        if((noinfo .eqv. .false.) .and. (scalinginfo .eqv. .true.)) write(79,*) time-equilib,&
             runningaveEtE/((time-equilib)/outputrate),&
             runningaveROG/((time-equilib)/outputrate)
     end if
     

     if (modulo(time,outputrate) == 0 .and. debugyes .eqv. .true.) then
        call debug
        !do ft = 1,nclients
           !write(109,*) time,sumclient(ft,1),sumclient(ft,2),sumclient(ft,3),&
                !sumclient(ft,4),sumclient(ft,5),sumclient(ft,6)
        !end do
        !call energy(.false.)
     end if
     
     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else if(modulo(time,outputrate/10) == 0) then
        write(6,*) 'step =',time
    write(9999,*) time,loc(:),locb(:),avesep/avesepcount,aveete/avesepcount
	 end if
     
     
     if((modulo(time,maxtime/10) == 0) .and.(noinfo .eqv. .false.)) then
        write(19,*) 'moves',macc,mrej, real(macc)/(macc+mrej)
        write(19,*) 'reptation',reacc,rerej,real(reacc)/(rerej+reacc)
        write(19,*) 'cluster translate',cltacc,cltatt,real(cltacc)/cltatt
        write(19,*) 'vector moves',roacc,roatt,real(roacc)/roatt   
        write(19,*) 'bond formations',bondformacccount,bondformrejcount,&
real(bondformacccount)/(bondformacccount+bondformrejcount)
        write(19,*) 'bond breaks',bonddeleteacccount,bonddeleterejcount,&
real(bonddeleteacccount)/(bonddeleterejcount+bonddeleteacccount)
     end if
     
     call system_clock(midpoint,count_rate,count_max)

     if((real(modulo(midpoint-start,count_max))/count_rate)> 216000) then
        backtime = time
        goto 93
     end if
     !write(6,*) 'sweep complete'
  end do


 if(restart .eqv. .false.) then
        backtime = timebase
     else
        backtime = timebase+lastruntime
     end if


93 if((film .eqv. .false.) .and. (scalinginfo .eqv. .false.)) call dataout
  write(6,*) 'deltax =',trackx, 'deltay =',tracky, 'deltaz =',trackz
     call system_clock(finish,count_rate,count_max)
  
call dataoutrestart
  if(noinfo .eqv. .false.) write(13,*) real(modulo(finish-start,count_max))/count_rate
  !call error
  call energy(.false.)
  if((noinfo .eqv. .true.) .or. (scalinginfo .eqv. .true.)) write(11,*) maxlength1, &
       (polymerl/((backtime-equilib)/outputrate)), &
       (polymerrog/((backtime-equilib)/outputrate))
  write(6,*) 2*log(actualchainlength/backtime)/log(real(maxlength))
  write(6,*) 'average clusper pop', (real(avecluspop)/aveclusnumber)/((backtime-equilib)/outputrate), &
       real(avecluspop)/((backtime-equilib)/outputrate), real(aveclusnumber)/((backtime-equilib)/outputrate)
  !inquire(file = "lengthdataaverages.dat", exist = exist)
  !if (exist) then
  ! open(30, file = "lengthdataaverages.dat", status = "old", position = "append",action ="write")
  !else
  !open(30, file = "lengthdataaverages.dat", status = "new",action ="write")
  !end if
  !write(30,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))!,2*(SQRT(totvar)/maxtime)/(actualchainlength/maxtime)

call cleanup

open(1717,file='radiusofgclient.dat',action='write')
write(1717,*) averageROG/noROGS
close(1717)
  
  write(6,*) 'reject = ', reject
  write(6,*) 'successful =', successful
  !write(6,*) 1.0*successful/(reject + successful)
  write(6,*) 'acceptance ratio =', 1.0*successful/(reject + successful)
  write(6,*) 'teleport acceptance=', telacc,telatt,real(telacc)/telatt
  write(6,*) 'count =', count
  write(6,*) 'maxlength =', maxlength
  write(6,*) 'forward =', reptforward
  write(6,*) 'backward =' , reptbackward
!  write(6,*) 'direct energy',interenergy(1,1),interenergy(1,2),interenergy(2,1),interenergy(2,2)
  if (finalfail .eqv. .true.) then
     write(6,*) '***********^^^^^^^^^^^^^^^^^^^This simulation FAILED!^^^^^^^^^^***********'
  end if
  !write(6,*) 'end'


do ft = 1,nscaffold
do qt = 1,chlen(ft)
if(protcoords(ft,qt)%x > gridsize) dright = dright+1
if(protcoords(ft,qt)%x <= gridsize) dleft = dleft+1
end do
end do

 open(1118,file = 'location.dat',action = 'write')
write(1118,*) loc(:)
close(1118)


!do qt = 1,clientlength,1
!write(6,*) 'Bindatt',qt,supercala(1,qt)
!end do

!do qt = 1,clientlength,1
!write(6,*) 'Positioning',qt,supercala(2,qt)
!end do


!do qt = 1,clientlength,1
!write(6,*) 'After Stretch ositioning',qt,supercala(3,qt)
!end do



write(6,*) 'LOC 1',loc(:),gridsize,gridsize
write(6,*) 'teleport attempts',pointless1,pointless2
write(6,*) 'Density Right then Left',dright,dleft
contains


subroutine unboundcounter
integer::m,l,maxl

unboundcount=0
boundcount=0

scaffoldarray(:,:)=0

    do m= nscaffold+1,nprotein,1
   maxl = chlen(m)
        do l=1,maxl,1

                if((protcoords(m,l)%am /=0) .and. (protcoords(m,l)%am <=(nscaffold))) then 
                        if(protcoords(m,l)%a111==0) then 
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+1
                        else if(protcoords(m,l)%a111==1)then
                                boundcount = boundcount+1        
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+1
                                scaffoldarray(protcoords(m,l)%am,protcoords(m,l)%al)=&
                                 scaffoldarray(protcoords(m,l)%am,protcoords(m,l)%al) +1    
                        end if
                end if

                if((protcoords(m,l)%bm /=0) .and. (protcoords(m,l)%bm<=(nscaffold)))then
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+2
                        else if(protcoords(m,l)%a111==2)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+2
                                scaffoldarray(protcoords(m,l)%bm,protcoords(m,l)%bl)=&
                                 scaffoldarray(protcoords(m,l)%bm,protcoords(m,l)%bl) +1    
                        end if
                end if

                if((protcoords(m,l)%cm /=0) .and.(protcoords(m,l)%cm<=(nscaffold)))then
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+3
                        else if(protcoords(m,l)%a111==3)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+3
                                scaffoldarray(protcoords(m,l)%cm,protcoords(m,l)%cl)=&
                                 scaffoldarray(protcoords(m,l)%cm,protcoords(m,l)%cl) +1    
                        end if
                end if


                if((protcoords(m,l)%dm /=0) .and.(protcoords(m,l)%dm<=(nscaffold)))then
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+4
                        else if(protcoords(m,l)%a111==4)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+4
                                scaffoldarray(protcoords(m,l)%dm,protcoords(m,l)%dl)=&
                                 scaffoldarray(protcoords(m,l)%dm,protcoords(m,l)%dl) +1    
                        end if
                end if

                if((protcoords(m,l)%em /=0) .and.(protcoords(m,l)%em<=(nscaffold)))then
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+5
                        else if(protcoords(m,l)%a111==5)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+5
                                scaffoldarray(protcoords(m,l)%em,protcoords(m,l)%el)=&
                                 scaffoldarray(protcoords(m,l)%em,protcoords(m,l)%el) +1    
                        end if
                end if

                if((protcoords(m,l)%fm /=0) .and.(protcoords(m,l)%fm<=(nscaffold)))then
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+6
                        else if(protcoords(m,l)%a111==6)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+6
                                scaffoldarray(protcoords(m,l)%fm,protcoords(m,l)%fl)=&
                                 scaffoldarray(protcoords(m,l)%fm,protcoords(m,l)%fl) +1    
                        end if
                end if



        end do
    end do


end subroutine unboundcounter


subroutine bondformatt
integer::pick,locator,m,l,dir,g,f,delbond,dir2,newlocator,maxl
double precision:: deltaenergy,denergy
logical::formlogic
integer::sumscaff
deltaenergy = 0.0d0



if(unboundcount==0) return


formlogic=.false.
pick = int(ran2(seed)*unboundcount)+1
locator=unboundarray(pick)
m=floor(real(locator-1)/(6*clientlength))+1
newlocator=locator-(6*clientlength*(m-1))
l=floor(real(newlocator-1)/6)+1

!supercala(2,l)=supercala(2,l)+1

if(protcoords(m,l)%a111>0) then
bondformrejcount= bondformrejcount+1
return
end if

dir=newlocator-(6*(l-1))

SELECT CASE (dir)

CASE (1)
g=protcoords(m,l)%am
f=protcoords(m,l)%al
dir2=2
CASE (2)
g=protcoords(m,l)%bm
f=protcoords(m,l)%bl
dir2=1
CASE (3)
g=protcoords(m,l)%cm
f=protcoords(m,l)%cl
dir2=4
CASE (4)
g=protcoords(m,l)%dm
f=protcoords(m,l)%dl
dir2=3
CASE(5)
g=protcoords(m,l)%em
f=protcoords(m,l)%el
dir2=6
CASE(6)
g=protcoords(m,l)%fm
f=protcoords(m,l)%fl
dir2=5
END SELECT

if(scaffoldarray(g,f)>0) then
bondformrejcount= bondformrejcount+1
!write(6,*) 'FAILLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Scaffold'
!stop 99
return
end if

if(g==m) then
bondformrejcount= bondformrejcount+1
return
end if



if(protcoords(g,f)%a111>0) then
bondformrejcount= bondformrejcount+1
return
end if

!if(interenergy(protcoords(m,l)%type,protcoords(g,f)%type) >=0) then
if(protcoords(m,l)%type==protcoords(g,f)%type) then
bondformrejcount= bondformrejcount+1
return
end if



    deltaenergy = interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    denergy=(real(unboundcount)/real(boundcount+1))

    !if(deltaenergy <= 0.0) then
    !    formlogic=.true.
    !end if
    if((denergy*exp(-deltaenergy/kT)) > ran2(seed)) then
        formlogic=.true.
    else
        bondformrejcount=bondformrejcount+1
    return
    end if



if (formlogic .eqv. .true.) then
    boundcount=boundcount+1
    unboundcount=unboundcount-1
    protcoords(m,l)%a111=dir
    protcoords(m,l)%a11n=g
    protcoords(g,f)%a111=dir2
    protcoords(g,f)%a11n=m
        maxl=chlen(m)
    boundarray(boundcount)=((m-1)*maxl*6)+((l-1)*6)+dir
totalenergy=totalenergy+deltaenergy
scaffoldarray(g,f)=scaffoldarray(g,f)+1

    do delbond=pick,unboundcount,1
        unboundarray(delbond)=unboundarray(delbond+1)
    end do
        unboundarray(unboundcount+1)=0
        bondformacccount=bondformacccount+1
end if

sumscaff=0
do g=1,nprotein-1,1
        maxl = chlen(g)
        do f= 1,maxl,1
          sumscaff=sumscaff+scaffoldarray(g,f)      
        end do
end do

if(sumscaff/=int(totalenergy/interenergy(1,2))) then
write(6,*) 'SCAFFOLD FAIL',totalenergy,sumscaff 
end if

end subroutine bondformatt

    subroutine sepcalc
    integer::m,l
    double precision::sep1,sep2,dx,dy,dz
    integer::maxl,etex,etey,etez,dxtest,dytest,dztest


    m = ((nprotein-nclients) +int(ran2(seed)*nclients))+1
 
etex=0
etey=0
etez=0

    maxl = chlen(m)

    do l=1,maxl-1,1
    if(l ==maxl) return


dxtest=0
dytest=0
dztest=0

dxtest=protcoords(m,l)%x-protcoords(m,l+1)%x
if(dxtest>protcoords(m,l)%linker) dxtest=dxtest-gridsize
if(dxtest<(-protcoords(m,l)%linker)) dxtest=gridsize-dxtest

dytest=protcoords(m,l)%y-protcoords(m,l+1)%y
if(dytest>protcoords(m,l)%linker) dytest=dytest-gridsize
if(dytest<(-protcoords(m,l)%linker)) dytest=gridsize-dytest

dztest=protcoords(m,l)%z-protcoords(m,l+1)%z
if(dztest>protcoords(m,l)%linker) dztest=dztest-gridsize
if (dztest<(-protcoords(m,l)%linker))dztest=gridsize-dztest


          dx = min(abs(protcoords(m,l)%x - protcoords(m,l+1)%x), gridsize-abs(protcoords(m,l)%x - protcoords(m,l+1)%x))
          dy = min(abs(protcoords(m,l)%y - protcoords(m,l+1)%y), gridsize-abs(protcoords(m,l)%y - protcoords(m,l+1)%y))
          dz = min(abs(protcoords(m,l)%z - protcoords(m,l+1)%z), gridsize-abs(protcoords(m,l)%z - protcoords(m,l+1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))

etex=etex+dxtest
etey=etey+dytest
etez=etez+dztest

                avesep=avesep+sep1
                avesepcount=avesepcount+1
     end do


sep2=sqrt(real(etex**2)+real(etey**2)+real(etez**2))
                aveete=aveete+sep2

  end subroutine sepcalc



subroutine bonddeleteatt
integer::pick,locator,m,l,dir,g,f,delbond,dir2,newlocator,maxl
logical:: formlogic
double precision:: deltaenergy,denergy

if(boundcount==0) return


deltaenergy = 0.0d0
formlogic = .false.
pick = int(ran2(seed)*boundcount)+1
locator=boundarray(pick)
m=floor(real(locator-1)/(6*clientlength))+1
newlocator=locator-(6*clientlength*(m-1))
l=floor(real(newlocator-1)/6)+1
if(protcoords(m,l)%a111==0) then
bonddeleterejcount= bonddeleterejcount+1
return
end if


!supercala(1,l) = supercala(1,l) +1

dir=newlocator-(6*(l-1))



SELECT CASE (dir)

CASE (1)
g=protcoords(m,l)%am
f=protcoords(m,l)%al
dir2=2
CASE (2)
g=protcoords(m,l)%bm
f=protcoords(m,l)%bl
dir2=1
CASE (3)
g=protcoords(m,l)%cm
f=protcoords(m,l)%cl
dir2=4
CASE (4)
g=protcoords(m,l)%dm
f=protcoords(m,l)%dl
dir2=3
CASE(5)
g=protcoords(m,l)%em
f=protcoords(m,l)%el
dir2=6
CASE(6)
g=protcoords(m,l)%fm
f=protcoords(m,l)%fl
dir2=5
END SELECT


if(g==m) then
bonddeleterejcount= bonddeleterejcount+1
return
end if



if(protcoords(g,f)%a111==0) then
bonddeleterejcount= bonddeleterejcount+1
return
end if

!if(interenergy(protcoords(m,l)%type,protcoords(g,f)%type) >=0) then
if(protcoords(m,l)%type==protcoords(g,f)%type)  then
bonddeleterejcount= bonddeleterejcount+1
return
end if



    deltaenergy = -1*interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    denergy=(real(boundcount)/real(unboundcount+1))

    !if(deltaenergy <= 0.0) then
    !    formlogic=.true.
    !end if
    if((denergy*exp(-deltaenergy/kT)) > ran2(seed)) then
        formlogic=.true.
    else
        bonddeleterejcount=bonddeleterejcount+1
    return
    end if



        if(formlogic .eqv. .true.) then
    unboundcount=unboundcount+1
    boundcount=boundcount-1
    protcoords(m,l)%a111=0
    protcoords(m,l)%a11n=g
    protcoords(g,f)%a111=0
    protcoords(g,f)%a11n=m
        maxl=chlen(m)
    unboundarray(unboundcount)=((m-1)*maxl*6)+((l-1)*6)+dir
totalenergy=totalenergy+deltaenergy
scaffoldarray(g,f)=scaffoldarray(g,f)-1
        
    do delbond=pick,boundcount,1
        boundarray(delbond)=boundarray(delbond+1)
    end do
        boundarray(boundcount+1)=0
        bonddeleteacccount= bonddeleteacccount+1
    end if


end subroutine bonddeleteatt

  subroutine cleanup
    integer::m,l,maxl,a

    if(restart .eqv. .false.) then
       open(119,file='visit.dat',action  = 'write')
       open(108,file = 'bonddate.dat',action = 'write')
       open(212,file = 'clientrdf.dat',action = 'write')
       !open(213,file = 'clusterclientrdf.dat',action = 'write')       
 else if (restart .eqv. .true.) then
    open(119,file='visit.dat',access = 'append')
          open(212,file = 'clientrdf.dat',access = 'append')
     !open(213,file = 'clusterclientrdf.dat',access = 'append')
          !open(108,file = 'bonddate.dat',access = 'append')
    !open(119,file='visit.dat',action  = 'write')
    
    open(108,file = 'bonddate.dat',action = 'write')
    end if
    do m= 1,nprotein,1
       maxl = chlen(m)
       do l = 1,maxl
          write(119,*) m,l,visit(m,l)
       end do
    end do
    do a = 1,nclients
       write(108,*) sumclient(a,:),sumsep/sumclient(a,clientlength+1),bulksumsep/sumclient(a,2*(clientlength+1)),unbound(a,:),&
        totalbindevents(1),totalbindevents(2),locb(1),locb(2),bindinghist(:)
    end do


    if(clientlength ==2) then
       write(178,*) sumclient(1,1)-sumclient(1,3),sumclient(1,2)-sumclient(1,3),&
            sumclient(1,4)-sumclient(1,6),sumclient(1,5)-sumclient(1,6)            
    end if
    
    if(clientlength == 3) then
    write(178,*) AB,AC,BC,ABC,bAB,bAC,bBC,bABC
    end if
    do l = 1,binsize*protcoords(nprotein,1)%linker
       write(141,*) real(l)/binsize,summing(1,l),summing(2,l)
    end do

    do l = 1,gridsize*binsize
       r = real(l)/binsize
       write(212,*) r,vdist(:,l)/(phasecount*(4.0/3.0)*PI_8*((r**3)-((real(l-1)/binsize)**3)))
       !write(213,*) r,clusdist(:,l)/(phasecount*(4.0/3.0)*PI_8*((r**3)-((real(l-1)/binsize)**3)))     
    end do

!open(251,file = 'freedist.dat',action = 'write')
!    do l = 1,(gridsize/2)*binsize
!       r = real(l)/binsize
!       write(251,*) r,freesites(:,l)/(normaliser*(4.0/3.0)*PI_8*((r**3)-((real(l-1)/binsize)**3)))    
!    end do

open(1003,file='z.dat',action='write')
write(1003,*) zloc(:)
close(1003)
    
  end subroutine cleanup

subroutine initialisearray
integer::m,l,maxl,axx,bxx,cxx



do axx=1,gridsize,1
do bxx=1,gridsize,1
do cxx=1,gridsize,1

cell(axx,bxx,cxx)%bm=0
end do
end do
end do

do m = 1,nscaffold,1
maxl = chlen(m)
do l=1,maxl

        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bm = m
        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl = l


end do
end do



end subroutine initialisearray


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


 subroutine singlemove(m)
    integer,intent(in)::m
    integer :: l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz,bdir,str
    double precision :: deltaenergy,tele
    logical :: rac,racc,adjver
    integer:: st,pr,maxl,maxback
    type(prottemp),dimension(:),allocatable :: tempcoord

    allocate(tempcoord(1))
    !return    
    rac = .true.


    maxl = 1
    l=1


if(protcoords(m,l)%a111>0) then
rac = .false.
goto 37
end if


tele = -0.5 !ran2(seed)-0.5

if(tele> 0.0) then


    
    l =1
    dx = 0
    dy = 0
    dz = 0



15  continue

    dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
    dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
    dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))



    if((dx**2 + dy**2 + dz**2) > protcoords(m,1)%linker) goto 15



    tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
    tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
    tempcoord(l)%z = protcoords(m,l)%z+dz

if(tempcoord(l)%z> gridsize) then
rac = .false.
goto 37
end if

if(tempcoord(l)%z<2) then
rac = .false.
goto 37
end if


else if (tele < 0.0) then

         tempcoord(l)%x = int(ran2(seed)*gridsize)+1
         tempcoord(l)%y = int(ran2(seed)*gridsize)+1
         tempcoord(l)%z = int(ran2(seed)*(gridsize-1))+2

end if


    deltaenergy = 0.0d0        
  
         
  

if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bm /=0) then
rac = .false.
goto 37
end if


    tempcoord(l)%am = 0
    tempcoord(l)%bm = 0
    tempcoord(l)%cm = 0
    tempcoord(l)%dm = 0
    tempcoord(l)%em = 0
    tempcoord(l)%fm = 0



    call removeenergyintra(m,l,deltaenergy,tempcoord)
    call removeenergyinter(m,l,deltaenergy,tempcoord)

  !  do pr = 1,nprotein
  !     if(pr /= m) then
  !        maxback = chlen(pr)
  !        do str = 1,maxback
  !           call adjacent(pr,l,str,tempcoord,adjver,bdir)
  !           if((adjver.eqv. .true.)) then
  !              call bondforminter(m,pr,l,str,tempcoord,tbo,bdir,deltaenergy)
  !           end if
  !        end do
  !     end if
  !  end do



                if((tempcoord(l)%x < gridsize).or.(tempcoord(l)%x>(gridsize-1))) then
if(cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0)then
                pr=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

                if((tempcoord(l)%x > 1).and.(tempcoord(l)%x < gridsize+2)) then
if(cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0)then
                pr=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if



        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm/=0) then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm/=0)then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

if(tempcoord(l)%z<gridsize) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bl
        bdir = -3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if


if(tempcoord(l)%z >1) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bl
        bdir = 3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if




    if(Energydecision(deltaenergy) .eqv. .true.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l,l,tempcoord,maxl)
       call updatebond(m,l,l,tempcoord)
       macc = macc + 1
        call unboundcounter
    else
       rac = .false.
       goto 37
    end if


37  if (rac .eqv. .false.) then
       mrej = mrej + 1
    end if
  end subroutine singlemove



   subroutine teleport

     integer :: l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz,bdir,str,m
     integer::xbk,ybk,zbk,i,j,k,g,f,hold
    double precision :: deltaenergy
    logical :: rac,racc,adjver
    integer:: st,pr,maxl,maxback
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer::gate,control
telatt=telatt+1
       m = ((nprotein-nclients) +int(ran2(seed)*nclients))+1
    maxl = chlen(m)

!return

do l=1,maxl
        if(protcoords(m,l)%a111>0) then
                return
        end if
end do

    allocate(tempcoord(maxl))
    !return    
    rac = .true.
deltaenergy = 0.0d0


    dx = 0
    dy = 0
    dz = 0







      hold = m
      !write(6,*) 'oooooh',m,protcoords(m,1)%linker
      
      control = 1
      gate = 1


      
      
35    continue
         l = 1
         i = int(ran2(seed)*gridsize)+1
         j = int(ran2(seed)*gridsize)+1

!CHANGE THIS IN OTHER CODES !!! %%%%%
         k = int(ran2(seed)*(gridsize-1))+2
!CHANGE THIS IN OTHER CODES !!! %%%%%

         !write(6,*) 'e'
        
!if(i>gridsize) pointless1 = pointless1 + 1
!if(i<=gridsize) pointless2 = pointless2 + 1

if(cell(i,j,k)%bm /=0) then
return
end if


 
!if(gridsizex /= 30) write(6,*) 'GRIDSIZE HAS CHANGED'

         tempcoord(l)%x = i
         tempcoord(l)%y = j
         tempcoord(l)%z = k

         !write(6,*) 'here'

         do l = 2,maxl

            xbk = tempcoord(l-1)%x
            ybk = tempcoord(l-1)%y
            zbk = tempcoord(l-1)%z


            control = 1

            
15          continue

!write(6,*) 'linkerlen',protcoords(m,l)%linker
               dx = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))
               dy = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))
               dz = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))



               if(sqrt(real((dx**2) + (dy**2) + (dz**2))) > protcoords(m,l-1)%linker) goto 15
               

                     i = modulo(xbk+dx-1,gridsize)+1
                     j = modulo(ybk +dy-1,gridsize)+1
                     k  = zbk+dz

if(k<2) return
if(k>gridsize) return


if(cell(i,j,k)%bm /=0) then
return
end if

               !do f = 1,l-1,1
                  !if (i == tempcoord(f)%x .and. j == tempcoord(f)%y &
                       !.and. k == tempcoord(f)%z) then
                     !return
                     !goto 15
                  !end if
               !end do

               tempcoord(l)%x = i
               tempcoord(l)%y = j
               tempcoord(l)%z = k

       
  
             end do
  
    
do l =1,maxl


    tempcoord(l)%am = 0
    tempcoord(l)%bm = 0
    tempcoord(l)%cm = 0
    tempcoord(l)%dm = 0
    tempcoord(l)%em = 0
    tempcoord(l)%fm = 0

!if(protcoords(m,l)%z==1) write(6,*) 'FAILLLLL Client on the bottom plane'


    call removeenergyintra(m,l,deltaenergy,tempcoord)
    call removeenergyinter(m,l,deltaenergy,tempcoord)



                if((tempcoord(l)%x < gridsize) .or.(tempcoord(l)%x>(gridsize-1))) then
if(cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0) then
                pr=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

                if((tempcoord(l)%x > 1).and.(tempcoord(l)%x < gridsize+2)) then
if(cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0) then
                pr=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if


        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm/=0) then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm/=0)then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

if(tempcoord(l)%z<gridsize) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bl
        bdir = -3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if


if(tempcoord(l)%z >1) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bl
        bdir = 3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if






 end do
 
    if(Energydecision(deltaenergy) .eqv. .true.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,1,maxl,tempcoord,maxl)
       call updatebond(m,1,maxl,tempcoord)
       telacc = telacc + 1
        call unboundcounter
    else
        return
        end if

!if(deltaenergy>0) then
!write(6,*) 'deltaenergy',deltaenergy
!call energy(.false.)
!end if

!37  if (rac .eqv. .false.) then
!       telrej = telrej + 1
!    end if
  end subroutine teleport

 subroutine rdftot(m)
   integer,intent(in)::m
    integer :: l,maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


!phasecountdoub = phasecountdoub+1

       maxl = chlen(m)
       do l = 1,maxl
          !phasecount = phasecount+1
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsize - abs(protcoords(m,l)%x - dpcomclusterinit%x))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - dpcomclusterinit%y))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - dpcomclusterinit%z))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


   vdist(clientlength+l,lx+1) = vdist(clientlength+l,lx+1)+1.0
   !vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0
       end do

    
  end subroutine rdftot


  

 subroutine rdf(m,l)
   integer,intent(in)::m,l
    integer :: maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


phasecount = phasecount+1

       maxl = chlen(m)
          dx = min(abs(protcoords(m,l)%x - ((3.0/4.0)*gridsize)), gridsize - abs(protcoords(m,l)%x - ((3.0/4.0)*gridsize)))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - (gridsize/2.0)))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - (gridsize/2.0)))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            vdist(l,lx+1) = vdist(l,lx+1)+1.0
   vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0

    
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
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsize - abs(protcoords(m,l)%x - dpcomclusterinit%x))
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
          dx = min(abs(protcoords(m,l)%x - dpcomclusterinit%x), gridsize - abs(protcoords(m,l)%x - dpcomclusterinit%x))
          dy = min(abs(protcoords(m,l)%y - dpcomclusterinit%y), gridsize - abs(protcoords(m,l)%y - dpcomclusterinit%y))
          dz = min(abs(protcoords(m,l)%z - dpcomclusterinit%z), gridsize - abs(protcoords(m,l)%z - dpcomclusterinit%z))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            clusdist(protcoords(m,l)%type,lx+1) = clusdist(protcoords(m,l)%type,lx+1)+1.0
   !vdist((2*clientlength)+1,lx+1) = vdist((2*clientlength)+1,lx+1) + 1.0
end do

    
  end subroutine rdfcluster

  
  subroutine analyse
    integer::m,l,g,f,control,dx,dy,dz,maxl,a,lx,dj
    double precision:: separation
    integer,dimension(clientlength) :: attached


    do a =1,clientlength
       attached = 0
    end do

!write(6,*) 'content',content,clnos(10)
    
    control = 0
    do a = 1,nclients
       m = a+nscaffold
    control = 0

       do g = 1,clientlength
           study(a,g) = 0
      end do
       maxl = chlen(m)
call comfind(m,.false.)

!  if(protcoords(m,1)%x >gridsize) then
!     loc(2) = loc(2)+1
!        lc = .false.
!  else if(protcoords(m,1)%x <=gridsize) then
!     loc(1) = loc(1)+1
!        lc = .true.
!     end if



do l = 1,maxl

!if(protcoords(m,l)%z==1) write(6,*) 'FAILLLLL Client on the bottom plane'

zloc(protcoords(m,l)%x) = zloc(protcoords(m,l)%x) + 1
                      if(protcoords(m,l)%a111>0) then

         dj = l
!call rdf(m,l)


!%%%%CHANGES - THE TWO LINES BELOW ARE NEW
g=protcoords(m,l)%a11n
visit(g,1)=visit(g,1)+1
         attached(control +1) = l
         study(a,l) = 1
         control = control + 1
         visit(m,l) = visit(m,l)+1
locb(1) = locb(1) +1


   end if

   
         
      end do

     bindinghist(control+1) = bindinghist(control+1) +1


if(control >0) then
!left = .true.


call radiusofgyration


do g = 1,clientlength,1
if(study(a,g) == 1) call rdf(m,g)
end do

end if
         


if(control == 0) then

totalbindevents(2) = totalbindevents(2) + 1
else if (control>0) then
totalbindevents(1) = totalbindevents(1) +1
end if


do g = 1,clientlength
if(study(a,g) ==0) unbound(a,g) =unbound(a,g)+1
end do



 
   end do
   

  end subroutine analyse

  subroutine phase(mz,clsize,clnos,cb)
    integer :: m,l,maxl,h,a,totcluspop,rogpop,maxrogpop,g
    integer::maxx,minx,maxy,miny,maxz,minz
    double precision :: radofgy,maxradofgy
    type(centremass)::avepos,sumsq
    integer,intent(inout)::clsize,mz
    integer,dimension(:),allocatable::cllist
    integer,dimension(:,:),allocatable:: route
    integer,dimension(:),allocatable::rcounter,checklist
    integer,dimension(:,:),intent(inout)::cb
    integer,dimension(:),intent(inout)::clnos
    type(basicp) :: comcluster,dummycom
    logical :: toolarge

    maxx=0
    minx=0
    maxy=0
    miny = 0
    maxz =0
    minz =0

  
    
    
    toolarge = .FALSE.
    rogpop = 0
    allocate(cllist(clsize))
    allocate(route(nprotein,clsize))
    allocate(rcounter(nprotein))
    allocate(checklist(clsize))

    cllist(1) =mz
    h = 2

    do a = 1,nprotein
       if(clnos(a) == clnos(mz)) then
          if(a/=mz) then
             cllist(h) = a
             h = h + 1
          end if
       end if
    end do
    call clustercom(mz,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge)
    !   write(6,*) 'comcluster',comcluster
    !  write(6,*) toolarge
    !if(toolarge .eqv. .false.) call clustercomcheck(mz,clsize,comcluster,cb,cllist)
    !write(6,*) 'check comcluster',comcluster
 !call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)

    
    
    radofgy = 0.0d0
    do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          sumsq%x = min(abs(protcoords(m,l)%x-comcluster%x),gridsize-abs(protcoords(m,l)%x-comcluster%x))
          sumsq%y = min(abs(protcoords(m,l)%y-comcluster%y),gridsize-abs(protcoords(m,l)%y-comcluster%y))
          sumsq%z = min(abs(protcoords(m,l)%z-comcluster%z),gridsize-abs(protcoords(m,l)%z-comcluster%z))
          radofgy = radofgy + (sumsq%x**2) + (sumsq%y**2) + (sumsq%z**2)
       end do
       rogpop = rogpop + maxl
    end do
    radofgy = radofgy/rogpop


    maxradofgy = 0.0d0
    do g = 1,clsize
       m = cllist(g)
       maxl = chlen(m)
       do l = 1,maxl
          sumsq%x = min(abs(protcoords(m,l)%x-comcluster%x),gridsize-abs(protcoords(m,l)%x-comcluster%x))
          sumsq%y = min(abs(protcoords(m,l)%y-comcluster%y),gridsize-abs(protcoords(m,l)%y-comcluster%y))
          sumsq%z = min(abs(protcoords(m,l)%z-comcluster%z),gridsize-abs(protcoords(m,l)%z-comcluster%z))
          maxradofgy = maxradofgy + (sumsq%x**2) + (sumsq%y**2) + (sumsq%z**2)
       end do
       maxrogpop = maxrogpop + maxl
    end do
    maxradofgy = maxradofgy/maxrogpop

    do g = 1,clsize
       m = cllist(g)
       maxl = chlen(m)
       do l = 1,maxl

          if(((protcoords(m,l)%x-comcluster%x)>maxx) .and. &
               ((protcoords(m,l)%x-comcluster%x)<(gridsize/2))) maxx =(protcoords(m,l)%x-comcluster%x)
          if((((gridsize+protcoords(m,l)%x)-comcluster%x)>maxx).and. &
               (((gridsize+protcoords(m,l)%x)-comcluster%x)<(gridsize/2))) then
          maxx =((gridsize+protcoords(m,l)%x)-comcluster%x)
       end if
       if(((comcluster%x-protcoords(m,l)%x)>minx) .and. &
           ((comcluster%x-protcoords(m,l)%x)<(gridsize/2))) minx =(comcluster%x-protcoords(m,l)%x)
          if((((gridsize+comcluster%x)-protcoords(m,l)%x)>minx).and. &
               (((gridsize+comcluster%x)-protcoords(m,l)%x)<(gridsize/2))) then
          minx =((gridsize+comcluster%x)-protcoords(m,l)%x)
       end if       

       if(((protcoords(m,l)%y-comcluster%y)>maxy) .and. &
           ((protcoords(m,l)%y-comcluster%y)<(gridsize/2))) maxy =(protcoords(m,l)%y-comcluster%y)
       if((((gridsize+protcoords(m,l)%y)-comcluster%y)>maxy).and. &
            (((gridsize+protcoords(m,l)%y)-comcluster%y)<(gridsize/2))) then
          maxy =((gridsize+protcoords(m,l)%y)-comcluster%y)
       end if
       if(((comcluster%y-protcoords(m,l)%y)>miny) .and. &
            ((comcluster%y-protcoords(m,l)%y)<(gridsize/2))) miny =(comcluster%y-protcoords(m,l)%y)
       if((((gridsize+comcluster%y)-protcoords(m,l)%y)>miny).and. &
            (((gridsize+comcluster%y)-protcoords(m,l)%y)<(gridsize/2))) then
          miny =((gridsize+comcluster%y)-protcoords(m,l)%y)
       end if

       if(((protcoords(m,l)%z-comcluster%z)>maxz) &
            .and. ((protcoords(m,l)%z-comcluster%z)<(gridsize/2))) maxz =(protcoords(m,l)%z-comcluster%z)
       if((((gridsize+protcoords(m,l)%z)-comcluster%z)>maxz).and. &
            (((gridsize+protcoords(m,l)%z)-comcluster%z)<(gridsize/2))) then
          maxz =((gridsize+protcoords(m,l)%z)-comcluster%z)
       end if
       if(((comcluster%z-protcoords(m,l)%z)>minz) &
            .and. ((comcluster%z-protcoords(m,l)%z)<(gridsize/2))) minz =(comcluster%z-protcoords(m,l)%z)
       if((((gridsize+comcluster%z)-protcoords(m,l)%y)>minz).and. &
            (((gridsize+comcluster%z)-protcoords(m,l)%z)<(gridsize/2))) then
          minz =((gridsize+comcluster%z)-protcoords(m,l)%z)
       end if
       
    end do
 end do

    write(160,*) maxx+minx,maxy+miny,maxz+minz
    
    
    write(77,*) time,radofgy,sqrt(radofgy),maxradofgy,sqrt(maxradofgy)

    
  end subroutine phase


     
  
  
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


subroutine foundationtwo
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z,topchains
    character(len = 10) :: BIN,nonsensea,nonsenseb,nonsensec
integer::totpoints,dummymax,g
doubleprecision :: rx,ry,rz,bondx,bondy,bondz


  open(121, file = pathtofile, action = 'read')  

!new way of finding totpoints
!make gas new types

topchains = 0
  totpoints = 0
  !write(6,*) 'nspecies',nspecies,speciespop(1),specieslen(1)
	topchains = topchains+speciespop(2)
     totpoints = totpoints + (speciespop(2)*specieslen(2))

  write(6,*) 'totpoints',totpoints,nprotein,nclients
    !read(21,*) BIN

totpoints = totpoints*2
 read(121,*) BIN,BIN,BIN,BIN,BIN,BIN
 do m = nprotein+1-(topchains+nclients),nprotein-nclients



    
     read(121,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,1)%type
!     read(121,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
protcoords(m,1)%species = protcoords(m,1)%type+2
     maxl = specieslen(protcoords(m,1)%type+2)
     chlen(m) = specieslen(protcoords(m,1)%type+2)
     count = count+1


     do l = 2,maxl

        read(121,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,l)%type
 !       read(121,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        protcoords(m,l)%species = protcoords(m,l)%type+2
        count = count+1
     end do

  end do


  dummymax = totpoints - topchains !(nprotein-nclients)
  do m = 1,dummymax
     read(121,*)BIN,BIN,BIN
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !accounts for 2 extra lines in energy than in move
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
  do t = 1,freezetime-1
     read(121,*) BIN
     read(121,*) nonsensea !BIN
     !write(6,*) 'a',nonsensea
     do m = 1,totpoints
        read(121,*)nonsenseb
        !write(6,*) 'b',nonsenseb
     end do
  end do

     read(121,*) BIN
        write(6,*) 'BIN',BIN
     read(121,*) BIN
  
     do m = nprotein+1-(topchains+nclients),nprotein-nclients,1
        read(121,*)  BIN,rx,ry,rz
        !write(6,*) BIN
protcoords(m,1)%linker = specieslinker(protcoords(m,1)%type+2)
        protcoords(m,1)%x = rx/4
        protcoords(m,1)%y = ry/4
        protcoords(m,1)%z = rz/4
protcoords(m,1)%a111=0
protcoords(m,1)%a11n=0



        protcoords(m,1)%x = (rx/4)+gridsize
        protcoords(m,1)%y = (ry/4)
        protcoords(m,1)%z = (rz/4)
        chlen(m) = specieslen(protcoords(m,1)%type+2)
        maxl = chlen(m)

        do l = 2,maxl,1
        protcoords(m,l)%linker = specieslinker(protcoords(m,l)%type+2)
           read(121,*)  BIN,rx,ry,rz

           protcoords(m,l)%x = rx/4
           protcoords(m,l)%y = ry/4
           protcoords(m,l)%z = rz/4
protcoords(m,1)%a111=0
protcoords(m,1)%a11n=0


        protcoords(m,l)%x = (rx/4)+gridsize
        protcoords(m,l)%y = (ry/4)
        protcoords(m,l)%z = (rz/4)
           !write(6,*) 'bond',chains(m,l)%bond
        end do
     end do

    call debug


   



  end subroutine foundationtwo

  subroutine clientfoundation
integer::g,m,l

    
        if(scalinginfo .eqv. .true.) call pdbsetup

    do g = (vtype-nclientspecies)+1,vtype,1
	!write(6,*) 'g',g,vtype,nclientspecies
       write(6,*) 'g !!!!!!!!!!!!!!!!!!!',g,client(g,1)%type,client(g,1)%species

       do m = (nprotein-nclients)+1,(nprotein-nclients)+speciespop(g),1
                 chlen(m) = specieslen(client(g,1)%species)
       protcoords(m,1)%species = client(g,1)%species
       protcoords(m,1)%type = client(g,1)%type
       protcoords(m,1)%linker = client(g,1)%linker
protcoords(m,1)%a111=0
protcoords(m,1)%a11n=0      
 
       do l =2,chlen(m)
          protcoords(m,l)%species = client(g,1)%species
          protcoords(m,l)%type = client(g,l)%type
          protcoords(m,l)%linker = client(g,l)%linker
protcoords(m,l)%a111=0
protcoords(m,l)%a11n=0

       end do
       
    end do
!write(6,*) 'clienttype',protcoords(441,2)%type

!    write(6,*) 'data',protcoords(441,3)%linker
    do m = (nprotein-nclients)+1,(nprotein-nclients)+speciespop(g),1
       write(6,*) 'm',m,'vtype',vtype
       call setupclient(m)
       write(6,*) 'finished client 1'
end do


   

end do


call energy(.true.)
    
  end subroutine clientfoundation
  
subroutine foundation
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z,topchains
    character(len = 10) :: BIN,nonsensea,nonsenseb,nonsensec
integer::totpoints,dummymax,g
doubleprecision :: rx,ry,rz,bondx,bondy,bondz

topchains = 0
  totpoints = 0
  write(6,*) 'nspecies',nspecies,speciespop(1),specieslen(1)
	topchains = topchains + speciespop(1)
     totpoints = totpoints + (speciespop(1)*specieslen(1))

  write(6,*) 'totpoints',totpoints,nprotein,nclients
    !read(21,*) BIN

totpoints = totpoints*2
 read(21,*) BIN,BIN,BIN,BIN,BIN,BIN
 do m = 1,topchains,1 !nprotein-nclients


    
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,1)%type
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
protcoords(m,1)%species = protcoords(m,1)%type
     maxl = specieslen(protcoords(m,1)%type)
     chlen(m) = specieslen(protcoords(m,1)%type)
     count = count+1
protcoords(m,1)%a111=0
protcoords(m,1)%a11n=0

     do l = 2,maxl
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,protcoords(m,l)%type
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        protcoords(m,l)%species = protcoords(m,l)%type
        count = count+1
protcoords(m,l)%a111=0
protcoords(m,l)%a11n=0

     end do

  end do


  dummymax = totpoints - topchains !(nprotein-nclients)

  do m = 1,dummymax
     read(21,*)BIN,BIN,BIN
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !accounts for 2 extra lines in energy than in move
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
  do t = 1,freezetime-1
     read(21,*) BIN
     read(21,*) nonsensea !BIN
     !write(6,*) 'a',nonsensea
     do m = 1,totpoints
        read(21,*)nonsenseb
        !write(6,*) 'b',nonsenseb
     end do
  end do

     read(21,*) BIN
        write(6,*) 'BIN',BIN
     read(21,*) BIN
  
     do m = 1,topchains,1 !nprotein-nclients,1
        read(21,*)  BIN,rx,ry,rz
        read(21,*) BIN,BIN,BIN,BIN
              !write(6,*) BIN
protcoords(m,1)%linker = specieslinker(protcoords(m,1)%type)
        protcoords(m,1)%x = rx/4
        protcoords(m,1)%y = ry/4
        protcoords(m,1)%z = rz/4


        chlen(m) = specieslen(protcoords(m,1)%type)
        maxl = chlen(m)

        do l = 2,maxl,1
        protcoords(m,l)%linker = specieslinker(protcoords(m,l)%type)
           read(21,*)  BIN,rx,ry,rz
                read(21,*) BIN,BIN,BIN,BIN

           protcoords(m,l)%x = rx/4
           protcoords(m,l)%y = ry/4
           protcoords(m,l)%z = rz/4

        end do
     end do

    !call debug
    if(scalinginfo .eqv. .true.) call pdbsetup

  end subroutine foundation



  subroutine setupclient(m)
    integer,intent(in)::m
    integer::maxl,hold,protpass,control,gate,i,j,k,dx,dy,dz,g,l,f
    integer::xbk,ybk,zbk,maxback
    
maxl = chlen(m)



      hold = m
      write(6,*) 'oooooh',m,protcoords(m,1)%linker,hold
      
      control = 1
      gate = 1
      protpass = 1
47    if (protpass ==1) then
         continue
      end if

      
      
35    continue
         l = 1
         i = int(ran2(seed)*gridsize)+1
         j = int(ran2(seed)*gridsize)+1
         k = int(ran2(seed)*(gridsize-1))+2
         !write(6,*) 'e'
         if (hold /= 1) then
            do g = 1, hold-1
        !write(6,*) 'g',g,protcoords(g,1)%species     
	  maxback = specieslen(protcoords(g,1)%species)
               do f = 1,maxback
                  if (i == protcoords(g,f)%x .and. j == protcoords(g,f)%y &
                       .and. k == protcoords(g,f)%z ) then
                     goto 35
                  else
                     continue
                  end if
               end do
            end do
         else
            continue
         end if

         protcoords(m,l)%x = i
         protcoords(m,l)%y = j
         protcoords(m,l)%z = k

         write(6,*) 'here'

         do l = 2,maxl

            xbk = protcoords(m,l-1)%x
            ybk = protcoords(m,l-1)%y
            zbk = protcoords(m,l-1)%z


            control = 1

            
15          continue

!write(6,*) 'linkerlen',protcoords(m,l)%linker
               dx = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))
               dy = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))
               dz = nint((ran2(seed)-0.5)*((2*protcoords(m,l-1)%linker)+1))



               if(sqrt(real((dx**2) + (dy**2) + (dz**2))) > protcoords(m,l-1)%linker) goto 15
               

                     i = modulo(xbk+dx-1,gridsize)+1
                     j = modulo(ybk +dy-1,gridsize)+1
                     k  = zbk+dz

                if(k<2) goto 15
                if(k>gridsize) goto 15

               do g = 1, hold-1
               maxback = specieslen(protcoords(g,1)%species)
                  do f = 1,maxback
                     if (i == protcoords(g,f)%x .and. j == protcoords(g,f)%y &
                          .and. k == protcoords(g,f)%z) then
                           goto 15
                        end if
                  end do
               end do

               do f = 1,maxl
                  if (i == protcoords(hold,f)%x .and. j == protcoords(hold,f)%y &
                       .and. k == protcoords(hold,f)%z) then
                        goto 15
                  end if
               end do

               protcoords(m,l)%x = i
               protcoords(m,l)%y = j
               protcoords(m,l)%z = k

       
  
             end do


do l=1,maxl,1
protcoords(m,l)%a111=0
protcoords(m,l)%a11n=0
end do

    
     end subroutine setupclient

    subroutine foundationrestart
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z,a
        double precision :: vb,vc
    character(len = 10) :: BIN
    read(23,*) lastruntime,totrmsbrute,suffclust
lastruntime = abs(lastruntime)
    do m = 1,nprotein,1
       read(23,*) protcoords(m,1)%species, protcoords(m,1)%x, protcoords(m,1)%y, protcoords(m,1)%z,&
            protcoords(m,1)%type,chlen(m),protcoords(m,1)%linker,protcoords(m,1)%a111,protcoords(m,1)%a11n

       maxl = chlen(m)

       do l = 2,maxl,1
          !write(6,*) l
          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z,&
               protcoords(m,l)%type,protcoords(m,l)%a111,protcoords(m,l)%a11n
       end do
    end do
open(108,file= 'bonddate.dat',action= 'read')
    do a = 1,nclients
       read(108,*)sumclient(a,:),vc,vb,unbound(a,:),totalbindevents(1),totalbindevents(2),locb(1),locb(2),&
bindinghist(:)
    end do
close(108)

open(108,file= 'location.dat',action= 'read')
       read(108,*) loc(:)
close(108)

    
do a = 1,nclients
sumsep = vc*sumclient(a,clientlength+1)
bulksumsep = vb*sumclient(a,2*(clientlength+1))
end do


open(108,file= 'z.dat',action= 'read')
    do a = 1,nclients
        read(108,*) zloc(:)
    end do
close(108)



    call energy(.TRUE.)
    call debug

    close(23)
  end subroutine foundationrestart

  

  subroutine rotate(m,l,b,tempcoords,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz)
    integer,intent(in):: m,l,b,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz
    type(prottemp),dimension(:),intent(inout) :: tempcoords
    !write(6,*) dx,dy,dz
    tempcoords(b)%x = modulo(protcoords(m,l)%x + (xx*dx) + (xy*dy) + (xz*dz) -1,gridsize)+1
    tempcoords(b)%y = modulo(protcoords(m,l)%y + (yx*dx) + (yy*dy) + (yz*dz)-1,gridsize)+1
    tempcoords(b)%z = modulo(protcoords(m,l)%z + (zx*dx) + (zy*dy) + (zz*dz)-1,gridsize)+1

  end subroutine rotate

 subroutine pickmoves
    integer:: scans,m,maxl,decide,l,a
real::bondformdelete

scans = 1

    decide = int(ran2(seed)*12)+1

    if(decide <= 6) then
       do a=1,clientlength
       bondformdelete = ran2(seed)-0.5
        if(bondformdelete>0.0) then
        call bondformatt
        else if(bondformdelete<0.0) then
        call bonddeleteatt
        end if

       end do
    end if
    if((decide > 6) .and. (decide <= 8)) call reptation
    if((decide > 8) .and. (decide <= 10)) call positioning
   if((decide > 10) .and. (decide <= 12)) call teleport


    
        if((modulo(time,maxtime/100) == 0)) then
        write(29,*) (time-equilib)*scans, real(totrmsbrute)/nprotein
     end if

   end subroutine pickmoves

   
  logical Function Energydecision(denergy)
    double precision, intent(in) :: denergy
    if(denergy <= 0.0) then
       Energydecision = .True.
    else if(exp(-denergy/kT) > ran2(seed)) then
       Energydecision = .True.
    else
       Energydecision = .False.
    end if
  end Function Energydecision


 

  subroutine positioning
    !performs a diagonal flip on a bead that has its two connecting beads perpendicular to each other
    integer :: m,l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz,bdir,str
    double precision :: deltaenergy
    logical :: rac,racc,adjver
    integer:: st,pr,maxl,maxback,decide
    !double precision::sumdebug
    type(prottemp),dimension(:),allocatable :: tempcoord
    !integer::fakex,fakey,fakez
    allocate(tempcoord(maxlength))
!return    
    rac = .true.

    
    m = ((nprotein-nclients) +int(ran2(seed)*nclients))+1
 
   maxl = chlen(m)

if(maxl ==1) then
       call singlemove(m)
       return
    end if

   l = int(ran2(seed)*maxl)+1
!supercala(2,l)=supercala(2,l)+1
!supercala(3,l)=supercala(3,l)+1


 
if(protcoords(m,l)%a111>0) then
       rac = .false.
        goto 37
    end if


!supercala(2,l)=supercala(2,l)+1

       dx = 0
       dy = 0
       dz = 0
       !write(6,*) 'positioning'


!15     continue
       
!               dx = nint((ran2(seed)-0.5)*(2*protcoords(m,1)%linker))
 !              dy = nint((ran2(seed)-0.5)*(2*protcoords(m,1)%linker))
 !              dz = nint((ran2(seed)-0.5)*(2*protcoords(m,1)%linker))

decide = int(ran2(seed)*3)+1

SELECT CASE (decide)

CASE (1)
dx = 1
CASE (2)
dx = 0
CASE (3)
dx = -1

END SELECT

decide = int(ran2(seed)*3)+1

SELECT CASE (decide)

CASE (1)
dy = 1
CASE (2)
dy = 0
CASE (3)
dy = -1

END SELECT

decide = int(ran2(seed)*3)+1

SELECT CASE(decide)

CASE (1)
dz = 1
CASE (2)
dz = 0
CASE (3)
dz = -1

END SELECT


               !if((dx**2 + dy**2 + dz**2) > protcoords(m,1)%linker) goto 15
               

!supercala(3,l)=supercala(3,l)+1


    tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
    tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
    tempcoord(l)%z = protcoords(m,l)%z+dz

        if(tempcoord(l)%z>gridsize) then
                rac = .false.
                goto 37
        end if

        if(tempcoord(l)%z<2) then
                rac = .false.
                goto 37
        end if
!supercala(3,l)=supercala(3,l)+1


if(l>1) then
    dx = min(abs(tempcoord(l)%x - protcoords(m,l-1)%x),gridsize-abs(tempcoord(l)%x &
         - protcoords(m,l-1)%x)) 
    dy = min(abs(tempcoord(l)%y - protcoords(m,l-1)%y),gridsize- abs(tempcoord(l)%y &
         - protcoords(m,l-1)%y)) 
    dz = min(abs(tempcoord(l)%z - protcoords(m,l-1)%z),gridsize - abs(tempcoord(l)%z &
         - protcoords(m,l-1)%z)) 

 sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
    if(sumdebug > protcoords(m,1)%linker) then
       rac = .false.
        goto 37
    end if
    end if

    if(l<maxl) then

     dx = min(abs(tempcoord(l)%x - protcoords(m,l+1)%x),gridsize-abs(tempcoord(l)%x &
         - protcoords(m,l+1)%x)) 
    dy = min(abs(tempcoord(l)%y - protcoords(m,l+1)%y),gridsize- abs(tempcoord(l)%y &
         - protcoords(m,l+1)%y)) 
    dz = min(abs(tempcoord(l)%z - protcoords(m,l+1)%z),gridsize - abs(tempcoord(l)%z &
         - protcoords(m,l+1)%z)) 

    sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
    if(sumdebug > protcoords(m,1)%linker) then
        rac = .false.
        goto 37
    end if
end if

!supercala(3,l)=supercala(3,l)+1


deltaenergy = 0.0d0        
       !do st = 1,maxl,1
          !if(st /=l) then
             !if(overlaps(m,l,st,tempcoord) .eqv. .false.) then
                !rac = .false.
                !goto 37
             !end if
          !end if
       !end do


if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bm /=0) then
rac = .false.
goto 37
end if



          tempcoord(l)%am = 0
          tempcoord(l)%bm = 0
          tempcoord(l)%cm = 0
          tempcoord(l)%dm = 0
          tempcoord(l)%em = 0
          tempcoord(l)%fm = 0

       

          call removeenergyintra(m,l,deltaenergy,tempcoord)
          call removeenergyinter(m,l,deltaenergy,tempcoord)

          do str = 1,l-3,1
             call adjacent(m,l,str,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                call bondformintra(m,l,str,tempcoord,bdir,deltaenergy)
             end if
          end do
          do str = l+3,maxl,1
             call adjacent(m,l,str,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                call bondformintra(m,l,str,tempcoord,bdir,deltaenergy)
             end if
          end do



                if((tempcoord(l)%x < gridsize) .or. (tempcoord(l)%x>(gridsize-1))) then
if(cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0) then
                pr=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

                if((tempcoord(l)%x > 1).and.(tempcoord(l)%x < gridsize+2)) then
if(cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0) then
                pr=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if



        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm/=0) then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm/=0)then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if

if(tempcoord(l)%z<gridsize) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z+1)%bl
        bdir = -3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if


if(tempcoord(l)%z >1) then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm/=0)then
                pr=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z-1)%bl
        bdir = 3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if



   !if(fakex /= tempcoord(l)%x) write(6,*) 'x change'
   !if(fakey /= tempcoord(l)%y) write(6,*) 'y change'
    !if(fakez /= tempcoord(l)%z) write(6,*) 'z change'

   
       
       !if(NINT(deltaenergy) /= 0) write(6,*) 'pivot',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l,l,tempcoord,maxl)
          call updatebond(m,l,l,tempcoord)
macc = macc + 1

call unboundcounter
!supercala(2,l)=supercala(2,l)+1

       else
          rac = .false.
          goto 37
       end if


37  if (rac .eqv. .false.) then
  mrej = mrej + 1
    end if
  end subroutine positioning



  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding


    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
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
    end if

    if(protcoords(m,l)%bm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%bl
       tempcoord(l)%bm = -1
    end if
    if(protcoords(m,l)%cm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
    end if
    if(protcoords(m,l)%dm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%dl
       tempcoord(l)%dm = -1
    end if
    if(protcoords(m,l)%em ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%el
       tempcoord(l)%em = -1
    end if
    if(protcoords(m,l)%fm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%fl
       tempcoord(l)%fm = -1
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
    end if

    if((protcoords(m,l)%bm /= 0)  .and. (protcoords(m,l)%bm /= m)) then
       g = protcoords(m,l)%bm
       f = protcoords(m,l)%bl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%bm = -1
    end if

    if((protcoords(m,l)%cm /= 0) .and. (protcoords(m,l)%cm /= m)) then
       g = protcoords(m,l)%cm
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if((protcoords(m,l)%dm /= 0) .and. (protcoords(m,l)%dm /= m)) then
       g = protcoords(m,l)%dm
       f = protcoords(m,l)%dl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%dm = -1
    end if

    if((protcoords(m,l)%em /= 0) .and. (protcoords(m,l)%em /= m)) then
       g = protcoords(m,l)%em
       f = protcoords(m,l)%el
       !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1)) write(6,*) 'successful bond rupture',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%em = -1
    end if

    if((protcoords(m,l)%fm /= 0) .and. (protcoords(m,l)%fm /= m)) then
       g = protcoords(m,l)%fm
       f = protcoords(m,l)%fl
       !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'successful rupture of bond',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%fm = -1
    end if

  end subroutine removeenergyinter

  subroutine bondformintra(m,l,st,tempcoord,bdir,deltaenergy)
    integer,intent(in)::m,l,st,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'start fail'
    if(bdir == -1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%am = m
       tempcoord(l)%al = st
    end if

    if(bdir == 1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%bm = m
       tempcoord(l)%bl = st
    end if

    if(bdir == -2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%cm = m
       tempcoord(l)%cl = st
    end if

    if(bdir == 2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%dm = m
       tempcoord(l)%dl = st
    end if

    if(bdir == -3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%em = m
       tempcoord(l)%el = st
    end if

    if(bdir == 3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%fm = m
       tempcoord(l)%fl = st
    end if
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'fail'
  end subroutine bondformintra

  subroutine bondforminter(m,g,l,st,tempcoord,bdir,deltaenergy)
    integer,intent(in)::m,l,st,g,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    !check this and cluster bonds update


    !write(6,*) 'm',m,'l',l,'g',g,'st',st
    if(bdir == -1) then
       tempcoord(l)%am = g
       tempcoord(l)%al = st
    end if

    if(bdir == 1) then
       tempcoord(l)%bm = g
       tempcoord(l)%bl = st
    end if

    if(bdir == -2) then
       tempcoord(l)%cm = g
       tempcoord(l)%cl = st
       !if(modulo(time,1000)==0) write(6,*) 'interform',bdir,tbo(l),bonddd(g,st)
    end if

    if(bdir == 2) then
       tempcoord(l)%dm = g
       tempcoord(l)%dl = st
    end if

    if(bdir == -3) then
       tempcoord(l)%em = g
       tempcoord(l)%el = st
    end if

    if(bdir == 3) then
       tempcoord(l)%fm = g
       tempcoord(l)%fl = st
    end if

  end subroutine bondforminter

  subroutine updatebondother(m,bead,tempcoord)
    integer,intent(in):: m,bead
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f


    if(tempcoord(bead)%am == -1) then
       g = protcoords(m,bead)%am
       f = protcoords(m,bead)%al
       protcoords(g,f)%bl = 0
       protcoords(g,f)%bm = 0
    end if

    if(tempcoord(bead)%bm == -1) then
       g = protcoords(m,bead)%bm
       f = protcoords(m,bead)%bl
       protcoords(g,f)%al = 0
       protcoords(g,f)%am = 0
    end if


    if(tempcoord(bead)%cm == -1) then
       g = protcoords(m,bead)%cm
       f = protcoords(m,bead)%cl
       protcoords(g,f)%dl = 0
       protcoords(g,f)%dm = 0
    end if

    if(tempcoord(bead)%dm == -1) then
       g = protcoords(m,bead)%dm
       f = protcoords(m,bead)%dl
       protcoords(g,f)%cl = 0
       protcoords(g,f)%cm = 0
    end if

    if(tempcoord(bead)%em == -1) then
       g = protcoords(m,bead)%em
       f = protcoords(m,bead)%el
       protcoords(g,f)%fl = 0
       protcoords(g,f)%fm = 0
    end if
    if(tempcoord(bead)%fm == -1) then
       g = protcoords(m,bead)%fm
       f = protcoords(m,bead)%fl
       protcoords(g,f)%el = 0
       protcoords(g,f)%em = 0
    end if


  end subroutine updatebondother




  subroutine updatebond(m,beadmin,beadmax,tempcoord)
    integer,intent(in):: m,beadmin,beadmax
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f,bg,bf

    !write(6,*) 'updateeeee'

    do l = beadmin,beadmax,1

       if(tempcoord(l)%am == -1) then
          g = protcoords(m,l)%am
          f = protcoords(m,l)%al
          !write(6,*) 'bond kill',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%bm == m) .and. (protcoords(g,f)%bl == l)) then
                protcoords(g,f)%bl = 0
                protcoords(g,f)%bm = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%bm == m) .and. (protcoords(g,f)%bl == l)) then
                   protcoords(g,f)%bl = 0
                   protcoords(g,f)%bm = 0
                end if
             end if
          end if
          protcoords(m,l)%al = 0
          protcoords(m,l)%am = 0
       else if(tempcoord(l)%am == 0) then
          protcoords(m,l)%am = tempcoord(l)%am
          protcoords(m,l)%al = tempcoord(l)%al
       else if(tempcoord(l)%am > 0) then
          g = tempcoord(l)%am
          f = tempcoord(l)%al
          bg = protcoords(m,l)%am
          bf = protcoords(m,l)%al
          if(bg /=0)then
             if((bg /= g) .or. (bf /=f)) then
                !does this work? may cause cancellation of newly formed bonds
                if((bg /= m)) then
                   if((protcoords(bg,bf)%bm == m) .and. (protcoords(bg,bf)%bl == l)) then
                      protcoords(bg,bf)%bm = 0
                      protcoords(bg,bf)%bl = 0
                   end if
                else if(bg ==m) then
                   if((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%bm == m) .and. (protcoords(bg,bf)%bl == l)) then
                         protcoords(bg,bf)%bm = 0
                         protcoords(bg,bf)%bl = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%am = tempcoord(l)%am
          protcoords(m,l)%al = tempcoord(l)%al

          !write(6,*) 'bond form',m,l,g,f
          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%bl = l
                protcoords(g,f)%bm = m
             end if
          else if(g/=m) then
             protcoords(g,f)%bl = l
             protcoords(g,f)%bm = m 
          end if
       end if


       if(tempcoord(l)%bm == -1) then
          g = protcoords(m,l)%bm
          f = protcoords(m,l)%bl
          !write(6,*) 'help',m,l,g,f,tempcoord(l)%bm
          !write(6,*) 'bond kill',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%am == m) .and. (protcoords(g,f)%al == l)) then
                protcoords(g,f)%al = 0
                protcoords(g,f)%am = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%am == m) .and. (protcoords(g,f)%al == l)) then             
                   protcoords(g,f)%al = 0
                   protcoords(g,f)%am = 0
                end if
             end if
          end if
          protcoords(m,l)%bl = 0
          protcoords(m,l)%bm = 0
       else if(tempcoord(l)%bm == 0) then
          protcoords(m,l)%bm = tempcoord(l)%bm
          protcoords(m,l)%bl = tempcoord(l)%bl
       else if(tempcoord(l)%bm > 0) then
          g = tempcoord(l)%bm
          f = tempcoord(l)%bl
          bg = protcoords(m,l)%bm
          bf = protcoords(m,l)%bl
          if(bg /=0) then
             if((bg /= g) .or. (bf /=f)) then
                if((bg /= m)) then
                   if((protcoords(bg,bf)%am == m) .and. (protcoords(bg,bf)%al == l)) then
                      protcoords(bg,bf)%am = 0
                      protcoords(bg,bf)%al = 0
                   end if
                else if(bg ==m) then
                   if((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%am == m) .and. (protcoords(bg,bf)%al == l)) then
                         protcoords(bg,bf)%am = 0
                         protcoords(bg,bf)%al = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%bm = tempcoord(l)%bm
          protcoords(m,l)%bl = tempcoord(l)%bl

          !write(6,*) 'bond form',m,l,g,f
          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%al = l
                protcoords(g,f)%am = m
             end if
          else if(g/=m) then
             protcoords(g,f)%al = l
             protcoords(g,f)%am = m 
          end if
       end if


       if(tempcoord(l)%cm == -1) then
          g = protcoords(m,l)%cm
          f = protcoords(m,l)%cl
          !write(6,*) 'bond kill',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%dm == m) .and. (protcoords(g,f)%dl == l)) then
                protcoords(g,f)%dl = 0
                protcoords(g,f)%dm = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%dm == m) .and. (protcoords(g,f)%dl == l)) then
                   protcoords(g,f)%dl = 0
                   protcoords(g,f)%dm = 0
                end if
             end if
          end if
          protcoords(m,l)%cl = 0
          protcoords(m,l)%cm = 0
       else if(tempcoord(l)%cm == 0) then
          protcoords(m,l)%cm = tempcoord(l)%cm
          protcoords(m,l)%cl = tempcoord(l)%cl
       else if(tempcoord(l)%cm > 0) then
          g = tempcoord(l)%cm
          f = tempcoord(l)%cl
          bg = protcoords(m,l)%cm
          bf = protcoords(m,l)%cl
          if(bg /=0) then
             if((bg /= g) .or. (bf /=f)) then
                if((bg /= m)) then
                   if((protcoords(bg,bf)%dm == m) .and. (protcoords(bg,bf)%dl == l)) then
                      protcoords(bg,bf)%dm = 0
                      protcoords(bg,bf)%dl = 0
                   end if
                else if(bg ==m) then
                   if ((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%dm == m) .and. (protcoords(bg,bf)%dl == l)) then
                         protcoords(bg,bf)%dm = 0
                         protcoords(bg,bf)%dl = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%cm = tempcoord(l)%cm
          protcoords(m,l)%cl = tempcoord(l)%cl
          !write(6,*) 'bond form',m,l,g,f

          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%dl = l
                protcoords(g,f)%dm = m
             end if
          else if(g/=m) then
             protcoords(g,f)%dl = l
             protcoords(g,f)%dm = m 
          end if
          end if


       if(tempcoord(l)%dm == -1) then
          g = protcoords(m,l)%dm
          f = protcoords(m,l)%dl
          !write(6,*) 'bond kill',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%cm == m) .and. (protcoords(g,f)%cl == l)) then
                protcoords(g,f)%cl = 0
                protcoords(g,f)%cm = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%cm == m) .and. (protcoords(g,f)%cl == l)) then
                   protcoords(g,f)%cl = 0
                   protcoords(g,f)%cm = 0
                end if
             end if
          end if
          protcoords(m,l)%dl = 0
          protcoords(m,l)%dm = 0
       else if(tempcoord(l)%dm == 0) then
          protcoords(m,l)%dm = tempcoord(l)%dm
          protcoords(m,l)%dl = tempcoord(l)%dl
       else if(tempcoord(l)%dm > 0) then
          g = tempcoord(l)%dm
          f = tempcoord(l)%dl
          bg = protcoords(m,l)%dm
          bf = protcoords(m,l)%dl
          if(bg /=0) then
             if((bg /= g) .or. (bf /=f)) then
                if((bg /= m)) then
                   if((protcoords(bg,bf)%cm == m) .and. (protcoords(bg,bf)%cl == l)) then
                      protcoords(bg,bf)%cm = 0
                      protcoords(bg,bf)%cl = 0
                   end if
                else if(bg==m) then
                   if ((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%cm == m) .and. (protcoords(bg,bf)%cl == l)) then
                         protcoords(bg,bf)%cm = 0
                         protcoords(bg,bf)%cl = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%dm = tempcoord(l)%dm
          protcoords(m,l)%dl = tempcoord(l)%dl
          !write(6,*) 'bond form',m,l,g,f
          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%cl = l
                protcoords(g,f)%cm = m
             end if
          else if(g/=m) then
             protcoords(g,f)%cl = l
             protcoords(g,f)%cm = m
             end if
       end if


       if(tempcoord(l)%em == -1) then
          g = protcoords(m,l)%em
          f = protcoords(m,l)%el
          !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1)) write(6,*) 'bond kill',m,l,g,f
          !write(6,*) 'help',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%fm == m) .and. (protcoords(g,f)%fl == l)) then
                protcoords(g,f)%fl = 0
                protcoords(g,f)%fm = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%fm == m) .and. (protcoords(g,f)%fl == l)) then
                   protcoords(g,f)%fl = 0
                   protcoords(g,f)%fm = 0
                end if
             end if
          end if
          protcoords(m,l)%el = 0
          protcoords(m,l)%em = 0

       else if(tempcoord(l)%em == 0) then
          protcoords(m,l)%em = tempcoord(l)%em
          protcoords(m,l)%el = tempcoord(l)%el
       else if(tempcoord(l)%em > 0) then
          g = tempcoord(l)%em
          f = tempcoord(l)%el
          bg = protcoords(m,l)%em
          bf = protcoords(m,l)%el
          if(bg /= 0) then
             if((bg /= g) .or. (bf /=f)) then
                if(bg /= m) then
                   if((protcoords(bg,bf)%fm == m) .and. (protcoords(bg,bf)%fl == l)) then
                      protcoords(bg,bf)%fm = 0
                      protcoords(bg,bf)%fl = 0
                   end if
                else if (bg ==m) then
                   if((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%fm == m) .and. (protcoords(bg,bf)%fl == l)) then
                         protcoords(bg,bf)%fm = 0
                         protcoords(bg,bf)%fl = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%em = tempcoord(l)%em
          protcoords(m,l)%el = tempcoord(l)%el
          !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1))write(6,*) 'bond form',m,l,g,f
          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%fl = l
                protcoords(g,f)%fm = m
             end if
         else if(g/=m) then
             protcoords(g,f)%fl = l
             protcoords(g,f)%fm = m
             end if
       end if

       if(tempcoord(l)%fm == -1) then
          g = protcoords(m,l)%fm
          f = protcoords(m,l)%fl
          !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'bond kill',m,l,g,f
          if(g/=m) then
             if((protcoords(g,f)%em == m) .and. (protcoords(g,f)%el == l)) then
                protcoords(g,f)%el = 0
                protcoords(g,f)%em = 0
             end if
          else if (g == m) then
             if((f<beadmin) .or. (f>beadmax)) then
                if((protcoords(g,f)%em == m) .and. (protcoords(g,f)%el == l)) then
                   protcoords(g,f)%el = 0
                   protcoords(g,f)%em = 0
                end if
             end if
          end if
          protcoords(m,l)%fl = 0
          protcoords(m,l)%fm = 0
       else if(tempcoord(l)%fm == 0) then
          protcoords(m,l)%fm = tempcoord(l)%fm
          protcoords(m,l)%fl = tempcoord(l)%fl
       else if(tempcoord(l)%fm > 0) then
          g = tempcoord(l)%fm
          f = tempcoord(l)%fl
          bg = protcoords(m,l)%fm
          bf = protcoords(m,l)%fl
          !if(l == 309) write(6,*) 'successssss',bg,bf
          if(bg /=0) then
             if((bg /= g) .or. (bf /=f)) then
                !if(bf == 298) write(6,*) 'successssss'
                if((bg /= m)) then
                   if((protcoords(bg,bf)%em == m) .and. (protcoords(bg,bf)%el == l)) then
                      protcoords(bg,bf)%em = 0
                      protcoords(bg,bf)%el = 0
                   end if
                else if(bg ==m) then
                   if((bf<beadmin) .or. (bf>beadmax)) then
                      if((protcoords(bg,bf)%em == m) .and. (protcoords(bg,bf)%el == l)) then
                         protcoords(bg,bf)%em = 0
                         protcoords(bg,bf)%el = 0
                      end if
                   end if
                end if
             end if
          end if
          protcoords(m,l)%fm = tempcoord(l)%fm
          protcoords(m,l)%fl = tempcoord(l)%fl
          !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'bond form',m,l,g,f
          if(g==m) then
             if((f<beadmin) .or. (f>beadmax)) then
                protcoords(g,f)%el = l
                protcoords(g,f)%em = m
                !                if(f == 310) write(6,*) 'results',protcoords(1,310)%el,protcoords(1,310)%em
             end if
         else if(g/=m) then
             protcoords(g,f)%el = l
             protcoords(g,f)%em = m
             end if
       end if
       !if(l==38) write(6,*) protcoords(1,43)%am,protcoords(1,43)%al
    end do

  end subroutine updatebond


  subroutine updatebondrep(m,beadmin,beadmax,tempcoord,del)
    integer,intent(in):: m,beadmin,beadmax,del
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f


    do l = beadmin,beadmax,1

       protcoords(m,l)%am = tempcoord(l)%am
       if(tempcoord(l)%am == m) then
          protcoords(m,l)%al = tempcoord(l)%al+del
       else
          protcoords(m,l)%al = tempcoord(l)%al
       end if
       if(tempcoord(l)%am > 0) then
          g = tempcoord(l)%am
          f = tempcoord(l)%al
          if(g ==m) f = f+del
          protcoords(g,f)%bl = l
          protcoords(g,f)%bm = m
       end if


       protcoords(m,l)%bm = tempcoord(l)%bm
       if(tempcoord(l)%bm == m) then
          protcoords(m,l)%bl = tempcoord(l)%bl+del
       else
          protcoords(m,l)%bl = tempcoord(l)%bl
       end if
       if(tempcoord(l)%bm > 0) then
          g = tempcoord(l)%bm
          f = tempcoord(l)%bl
          if(g ==m) f = f+del
          protcoords(g,f)%al = l
          protcoords(g,f)%am = m
       end if



       protcoords(m,l)%cm = tempcoord(l)%cm
       if(tempcoord(l)%cm == m) then
          protcoords(m,l)%cl = tempcoord(l)%cl+del
       else
          protcoords(m,l)%cl = tempcoord(l)%cl
       end if
       if(tempcoord(l)%cm > 0) then
          g = tempcoord(l)%cm
          f = tempcoord(l)%cl
          if(g ==m) f = f+del          
          protcoords(g,f)%dl = l
          protcoords(g,f)%dm = m
       end if



       protcoords(m,l)%dm = tempcoord(l)%dm
       if(tempcoord(l)%dm == m) then
          protcoords(m,l)%dl = tempcoord(l)%dl+del
       else
          protcoords(m,l)%dl = tempcoord(l)%dl
       end if
       if(tempcoord(l)%dm > 0) then
          g = tempcoord(l)%dm
          f = tempcoord(l)%dl
          if(g ==m) f = f+del
          protcoords(g,f)%cl = l
          protcoords(g,f)%cm = m
       end if


       protcoords(m,l)%em = tempcoord(l)%em
       if(tempcoord(l)%em == m) then
          protcoords(m,l)%el = tempcoord(l)%el+del
       else
          protcoords(m,l)%el = tempcoord(l)%el
       end if
       if(tempcoord(l)%em > 0) then
          g = tempcoord(l)%em
          f = tempcoord(l)%el
          if(g ==m) f = f+del
          protcoords(g,f)%fl = l
          protcoords(g,f)%fm = m
       end if



       protcoords(m,l)%fm = tempcoord(l)%fm
       if(tempcoord(l)%fm == m) then
          protcoords(m,l)%fl = tempcoord(l)%fl+del
       else
          protcoords(m,l)%fl = tempcoord(l)%fl
       end if
       if(tempcoord(l)%fm > 0) then
          g = tempcoord(l)%fm
          f = tempcoord(l)%fl
          if(g ==m) f = f+del
          protcoords(g,f)%el = l
          protcoords(g,f)%em = m
       end if

    end do

  end subroutine updatebondrep



  subroutine reptation
    !perform reptation
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir,maxl,maxback,str
    logical :: reptcont,adjver
    type(prottemp),dimension(:),allocatable :: tempcoord
!return

    !write(6,*) 'for information',protcoords(20,6)%dm,protcoords(20,6)%dl
    t = time + 1

    !choose = -0.5
    choose = ran2(seed) - 0.5
    !write(6,*) 'reptation start',choose

    reptcont = .true.
     m = ((nprotein-nclients) +int(ran2(seed)*nclients))+1



     maxl = chlen(m)

        do l=1,maxl
                if(protcoords(m,l)%a111>0) then
                rerej=rerej+1  
                reject = reject+1      
                return

                end if
        end do

     if(maxl ==1) then
       call singlemove(m)
       return
    end if
    !write(6,*) 'maxl',maxl,m
    allocate(tempcoord(maxl))


   
    direc = int(ran2(seed)*5) +1
    !if(time == 3387) write(6,*) 'choooooose =', choose,direc,m
    if (choose > 0.0) then
	!if(choose == 1) then
19     continue

       
               dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))



               if(sqrt(real((dx**2 + dy**2 + dz**2))) > protcoords(m,1)%linker) goto 19



      
       tempcoord(1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(1)%z = protcoords(m,1)%z + dz

if(tempcoord(1)%z>gridsize) then
       reptcont = .false.
       goto 83
    end if

if(tempcoord(1)%z<2) then
       reptcont = .false.
       goto 83
    end if


     dx = min(abs(tempcoord(1)%x - protcoords(m,2)%x),gridsize-abs(tempcoord(1)%x &
         - protcoords(m,2)%x)) 
    dy = min(abs(tempcoord(1)%y - protcoords(m,2)%y),gridsize- abs(tempcoord(1)%y &
         - protcoords(m,2)%y)) 
    dz = min(abs(tempcoord(1)%z - protcoords(m,2)%z),gridsize - abs(tempcoord(1)%z &
         - protcoords(m,2)%z)) 

    sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
    if(sumdebug > real(protcoords(m,1)%linker)) then
       reptcont = .false.
       goto 83
    end if
 
       ! write(6,*) tbdi(1),dx,dy,dz
       ! write(6,*) tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z
       ! write(6,*) protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z

       deltaenergy = 0.0d0

       !do st = 1,maxl-1,1
          !if(overlaps(m,1,st,tempcoord) .eqv. .false.) then
             !reptcont = .false.
             !goto 83
          !end if
       !end do

  

if(cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z)%bm /=0) then
reptcont = .false.
goto 83
end if


       
       tempcoord(maxl)%am = 0
       tempcoord(maxl)%bm = 0
       tempcoord(maxl)%cm = 0
       tempcoord(maxl)%dm = 0
       tempcoord(maxl)%em = 0
       tempcoord(maxl)%fm = 0
       call removeenergyintra(m,maxl,deltaenergy,tempcoord)
       call removeenergyinter(m,maxl,deltaenergy,tempcoord)

       tempcoord(1)%am = 0
       tempcoord(1)%bm = 0
       tempcoord(1)%cm = 0
       tempcoord(1)%dm = 0
       tempcoord(1)%em = 0
       tempcoord(1)%fm = 0

       do st = 3,maxl-1
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,1,st,tempcoord,bdir,deltaenergy)
          end if
       end do

       !if(deltaenergy<0.0)write(6,*) 'energy post intra=',deltaenergy



        if((tempcoord(1)%x < gridsize) .or. (tempcoord(1)%x>(gridsize-1))) then

if(cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=0)then
                pr=cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm
                str=cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bl
        bdir = -1

                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if


        if((tempcoord(1)%x > 1).and.(tempcoord(1)%x < gridsize+2)) then
if(cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=0) then
                pr=cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm
                str=cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)

end if
end if



        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bm/=0) then
                pr=cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bm
                str=cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bl
        bdir = -2

                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if

        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bm/=0)then
                pr=cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bm
                str=cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if


if(tempcoord(1)%z<gridsize) then
        if(cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z+1)%bm/=0) then
                pr=cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z+1)%bm
                str=cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z+1)%bl
        bdir = -3
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if

if(tempcoord(1)%z>1) then
        if(cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z-1)%bm/=0) then
                pr=cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z-1)%bm
                str=cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z-1)%bl
        bdir = 3
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if









       if(Energydecision(deltaenergy) .eqv. .true.) then

          totalenergy = totalenergy + deltaenergy
          do l =maxl,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z
          end do
          call updatebondother(m,maxl,tempcoord)
          call updatepos(m,1,1,tempcoord,maxl)
          call reptationadjust(m,maxl,2,tempcoord,-1)
          call updatebondrep(m,1,1,tempcoord,1)
          successful = successful + 1
          reptforward = reptforward + 1
reacc = reacc +1 
call unboundcounter     
 else
          reptcont = .false.
       end if
    else if (choose < 0.0) then
	!else if (choose ==2) then
      
15     continue
               
               dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))



               if(sqrt(real((dx**2 + dy**2 + dz**2))) > protcoords(m,1)%linker) goto 15


       tempcoord(maxl)%x = modulo(protcoords(m,maxl)%x + dx-1,gridsize)+1
       tempcoord(maxl)%y = modulo(protcoords(m,maxl)%y + dy-1,gridsize)+1
       tempcoord(maxl)%z = protcoords(m,maxl)%z + dz

if(tempcoord(maxl)%z > gridsize) then
reptcont=.false.
goto 83
end if

if(tempcoord(maxl)%z <2) then
reptcont=.false.
goto 83
end if


     dx = min(abs(tempcoord(maxl)%x - protcoords(m,maxl-1)%x),gridsize-abs(tempcoord(maxl)%x &
         - protcoords(m,maxl-1)%x)) 
    dy = min(abs(tempcoord(maxl)%y - protcoords(m,maxl-1)%y),gridsize- abs(tempcoord(maxl)%y &
         - protcoords(m,maxl-1)%y)) 
    dz = min(abs(tempcoord(maxl)%z - protcoords(m,maxl-1)%z),gridsize - abs(tempcoord(maxl)%z &
         - protcoords(m,maxl-1)%z)) 

    sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
    if(sumdebug > real(protcoords(m,1)%linker)) then
       reptcont = .false.
       goto 83
    end if
       
       deltaenergy = 0.0d0


       !do st = 2,maxl,1
          !if(overlaps(m,maxl,st,tempcoord) .eqv. .false.) then
             !reptcont = .false.
             !goto 83
          !end if
       !end do



if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm /=0) then
reptcont = .false.
goto 83
end if



       
       tempcoord(1)%am = 0
       tempcoord(1)%bm = 0
       tempcoord(1)%cm = 0
       tempcoord(1)%dm = 0
       tempcoord(1)%em = 0
       tempcoord(1)%fm = 0

       call removeenergyintra(m,1,deltaenergy,tempcoord)
       call removeenergyinter(m,1,deltaenergy,tempcoord)

       tempcoord(maxl)%am = 0
       tempcoord(maxl)%bm = 0
       tempcoord(maxl)%cm = 0
       tempcoord(maxl)%dm = 0
       tempcoord(maxl)%em = 0
       tempcoord(maxl)%fm = 0

       do st = 2,maxl-2
          call adjacent(m,maxl,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,maxl,st,tempcoord,bdir,deltaenergy)
          end if
       end do


     if((tempcoord(maxl)%x < gridsize) .or. (tempcoord(maxl)%x>(gridsize-1)))then
        if(cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=0)then
                pr=cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm
                str=cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bl
        bdir = -1
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if

                if((tempcoord(maxl)%x > 1).and.(tempcoord(maxl)%x <(gridsize+2))) then
if(cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=0)then
                pr=cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm
                str=cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bl
        bdir = 1
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if


        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bm/=0)then
                pr=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bm
                str=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bl
        bdir = -2
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)

end if

        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bm/=0)then
                pr=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bm
                str=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bl
        bdir = 2
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)

end if


if(tempcoord(maxl)%z<gridsize) then
        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z+1)%bm/=0)then
                pr=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z+1)%bm
                str=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z+1)%bl
        bdir = -3
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)

end if
end if

if(tempcoord(maxl)%z>1) then
        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z-1)%bm/=0)then
                pr=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z-1)%bm
                str=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z-1)%bl
        bdir = 3
                call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)

end if
end if







       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy

          do l =1,maxl-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z
          end do
          
          call updatebondother(m,1,tempcoord)
          tempcoord(maxl-1)%x = protcoords(m,maxl-1)%x
          tempcoord(maxl-1)%y = protcoords(m,maxl-1)%y
          tempcoord(maxl-1)%z = protcoords(m,maxl-1)%z
        
          call updatepos(m,maxl-1,maxl,tempcoord,maxl)
          call reptationadjust(m,1,maxl-1,tempcoord,1)
          call updatebondrep(m,maxl,maxl,tempcoord,-1)
          successful = successful + 1
          reptbackward = reptbackward + 1
reacc = reacc + 1  
call unboundcounter     
else
          reptcont = .false.
       end if
    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
rerej = rerej +1
    end if


  end subroutine reptation



  subroutine reptationadjust(m,beadmin,beadmax,tempcoord,del)
    integer,intent(in)::m,beadmin,beadmax,del
    type(prottemp),dimension(:),intent(inout)::tempcoord
    integer :: g,f,l,buzz

    if(del == -1) buzz = chlen(m)
    if(del ==1) buzz = 1

    do l =beadmin,beadmax,del
       protcoords(m,l)%am = protcoords(m,l+del)%am
       g = protcoords(m,l+del)%am
       if(g /=0) then
          if(g == m) then
             if(protcoords(m,l+del)%al /= buzz) then
                protcoords(m,l)%al = protcoords(m,l+del)%al - del
             else
                protcoords(m,l)%am = 0
             end if
          else if(g/=m) then
             protcoords(m,l)%al = protcoords(m,l+del)%al
             f = protcoords(m,l)%al
             protcoords(g,f)%bl = l
          end if
       end if

       protcoords(m,l)%bm = protcoords(m,l+del)%bm
       g = protcoords(m,l+del)%bm
       if(g /=0) then
          if(g == m) then

             if(protcoords(m,l+del)%bl /= buzz) then
                protcoords(m,l)%bl = protcoords(m,l+del)%bl - del
             else
                protcoords(m,l)%bm = 0
             end if
          else if(g /=m) then
             protcoords(m,l)%bl = protcoords(m,l+del)%bl
             f = protcoords(m,l)%bl
             protcoords(g,f)%al = l
          end if
       end if

       protcoords(m,l)%cm = protcoords(m,l+del)%cm
       g = protcoords(m,l+del)%cm
       if(g /=0) then
          if(g == m) then
             if(protcoords(m,l+del)%cl /= buzz) then
                protcoords(m,l)%cl = protcoords(m,l+del)%cl - del
             else
                protcoords(m,l)%cm = 0
             end if

          else if (g/=m) then
             protcoords(m,l)%cl = protcoords(m,l+del)%cl
             f = protcoords(m,l)%cl
             protcoords(g,f)%dl = l
          end if
       end if

       protcoords(m,l)%dm = protcoords(m,l+del)%dm
       g = protcoords(m,l+del)%dm
       if(g /=0) then
          if(g == m) then
           
             if(protcoords(m,l+del)%dl /= buzz) then
                protcoords(m,l)%dl = protcoords(m,l+del)%dl - del
             else
                protcoords(m,l)%dm = 0
             end if

          else if(g/=m) then
             protcoords(m,l)%dl = protcoords(m,l+del)%dl
             f = protcoords(m,l)%dl
             protcoords(g,f)%cl = l
          end if
       end if

       protcoords(m,l)%em = protcoords(m,l+del)%em
       g = protcoords(m,l+del)%em
       if(g /=0) then
          if(g == m) then
            
             if(protcoords(m,l+del)%el /= buzz) then
                protcoords(m,l)%el = protcoords(m,l+del)%el - del
             else
                protcoords(m,l)%em = 0
             end if


          else if(g/=m) then
             protcoords(m,l)%el = protcoords(m,l+del)%el
             f = protcoords(m,l)%el
             protcoords(g,f)%fl = l
          end if
       end if

       protcoords(m,l)%fm = protcoords(m,l+del)%fm
       g = protcoords(m,l+del)%fm
       if(g /=0) then
          if(g == m) then
            

             if(protcoords(m,l+del)%fl /= buzz) then
                protcoords(m,l)%fl = protcoords(m,l+del)%fl - del
             else
                protcoords(m,l)%fm = 0
             end if

          else if(g/=m) then
             protcoords(m,l)%fl = protcoords(m,l+del)%fl
             f = protcoords(m,l)%fl
             protcoords(g,f)%el = l
          end if
       end if


    end do
  end subroutine reptationadjust

  

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







  subroutine energy(init)
    integer :: dx,dy,dz,m,l,f,g,bdir,delx,dely,delz,maxlengthss,maxl
    double precision :: initialenergy,olderenergy,dumoldenergy
    Type(prottemp),dimension(:),allocatable :: tempcoord
  !double precision,dimension(:),intent(inout)::chen
    logical :: adjver
    logical,intent(in)::init
integer::bdir2,bdir3
    initialenergy = 0.0d0
    !totalenergy = 0.0d0
    allocate(tempcoord(maxlength))


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
          end do
          do l = 1,maxl,1
             do f = 1,maxl
                if(abs(f-l)>2) then
                   call adjacent(m,l,f,tempcoord,adjver,bdir)
                   if((adjver .eqv. .true.)) then
                      call bondformintra(m,l,f,tempcoord,bdir,initialenergy)
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
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy)
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
          end do

          if(m<nprotein) then
             do g = m+1,nprotein,1
                maxlengthss = chlen(g)
                do l = 1,maxl,1
                   do f = 1,maxlengthss,1
                      call adjacent(g,l,f,tempcoord,adjver,bdir)
                      if((adjver .eqv. .true.)) then


SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
CASE (1)
bdir2=2
bdir3=1
CASE (-2)
bdir2=3
bdir3=4
CASE (2)
bdir2=4
bdir3=3
CASE(-3)
bdir2=5
bdir3=6
CASE(3)
bdir2=6
bdir3=5
END SELECT



                         dumoldenergy = initialenergy
                         if((protcoords(m,l)%a111==bdir2).and.(protcoords(g,f)%a111==bdir3) .and. &
                        (protcoords(m,l)%a11n==g) .and.(protcoords(g,f)%a11n==m)) then
                                initialenergy =initialenergy+interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                           end if
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy)
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

if((tempory(ch1)%z == gridsize) .and. (protcoords(pr,ch2)%z==1)) then
        adjver = .false.
goto 31
else if((tempory(ch1)%z==1) .and. (protcoords(pr,ch2)%z==gridsize)) then
        adjver=.false.
        goto 31
end if



    delx = tempory(ch1)%x - protcoords(pr,ch2)%x
    if((abs(delx) == 1) .or. (abs(delx) == (gridsize - 1))) then
       dx = 1
       bdir = delx
if(abs(delx) == (gridsize - 1)) bdir = (delx/(gridsize-1))*(gridsize-1)
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

if((tempory(ch1)%z == gridsize) .and. (protcoords(pr,ch2)%z==1)) then
        adjver = .false.
goto 31
else if((tempory(ch1)%z==1) .and. (protcoords(pr,ch2)%z==gridsize)) then
        adjver=.false.
        goto 31
end if



    delx = tempory(ch1)%x - protcoords(pr,ch2)%x
    if((abs(delx) == 1) .or. (abs(delx) == (gridsize - 1))) then
       dx = 1
       bdir = delx
!%%%%CHANGES - CHANGED GRIDSIZE TO GRIDSIZEX
if(abs(delx) == (gridsize - 1)) bdir = (delx/(gridsize-1))*(gridsize-1)
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
           if (dcomx >protcoords(m,1)%linker) dcomx = dcomx-gridsize
           if (dcomx <(-protcoords(m,1)%linker)) dcomx = gridsize+dcomx
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
        com(1,m)%x = modulo(protcoords(m,1)%x +(real(comx)/(maxl)) -1,real(gridsize))+1
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
          if (dcx > gridsize/2) dcx = dcx-gridsize
          if (dcx < -gridsize/2) dcx = gridsize + dcx
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
         (gridsize - abs(com(1,m)%x -com(2,m)%x))))**2
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
          rog = (min(modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0))))**2+ &
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
  
  
  subroutine length
    integer :: ddx,ddy,ddz,dxsum,dysum,dzsum,m,l,maxl
    double precision :: normchainl
    t = time
    totchainlength = 0.0
    normchainl = 0.0d0

    do m = 1,nprotein,1
       
        maxl = chlen(m)
        dxsum = min(abs(protcoords(m,maxl)%x - protcoords(m,1)%x),gridsize-abs(protcoords(m,maxl)%x - protcoords(m,1)%x))
        dysum = min(abs(protcoords(m,maxl)%y - protcoords(m,1)%y),gridsize-abs(protcoords(m,maxl)%y - protcoords(m,1)%y))
        dzsum = min(abs(protcoords(m,maxl)%z - protcoords(m,1)%z),gridsize-abs(protcoords(m,maxl)%z - protcoords(m,1)%z))
        
       chainlength = ((dxsum**2) + (dysum**2) + (dzsum**2))
       totchainlength = totchainlength + chainlength
    end do
    !totchainlength = sum(chainlength)
    runningaveEtE = runningaveEtE + totchainlength  !sort this out
    normchainl = totchainlength/nprotein
polymerl = polymerl + normchainl
    !actualchainlength = actualchainlength + avechainlength
    write(91,*) t,normchainl
end subroutine length

  subroutine dataout
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit
    integer :: choice
    allocate(bondposit(nprotein,maxlength))
    !write(6,*) 'data'
    t = time
    if(suffclust .eqv. .false.) choice = 67
    if(suffclust .eqv. .true.) choice = 88
    write(choice,*) ' '
    write(choice,*) 'timestep ', 'indexed'
    write(choice,*) 'pbc', real(4*gridsize),real(4*gridsize),real(4*gridsize)
    !write(67,*) ((nprotein1*maxlength1) + (nprotein2*maxlength2))
    acount = 0

    do m = 1,nprotein,1
       maxl = chlen(m)
       do l = 1,maxl,1
          write(choice,*) acount,real(4*protcoords(m,l)%x), real(4*protcoords(m,l)%y),real(4*protcoords(m,l)%z) 
          acount = acount +1

71        if(m == m) continue
       end do
    end do



  end subroutine dataout

  subroutine dataoutrestart
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit

    open(59, file = 'checkpoint.xyz', action = 'write')
    write(59,*) backtime,totrmsbrute,suffclust
    do m = 1,nprotein,1
       maxl = chlen(m)
       write(59,*) protcoords(m,1)%species,protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z,&
            protcoords(m,1)%type,maxl,protcoords(m,1)%linker,protcoords(m,1)%a111 &
,protcoords(m,1)%a11n
       do l = 2,maxl,1
          write(59,*) protcoords(m,l)%species,protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z,protcoords(m,l)%type, &
protcoords(m,l)%a111,protcoords(m,l)%a11n 
       end do
    end do


  end subroutine dataoutrestart





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


          !do f = 1,nprotein,1
             !maxlengthss = chlen(f)
             !do g = 1,maxlengthss,1
                !if(m ==f .and. l == g) then
                   !continue
                !else if ((m/= f) .or. (l/= g)) then
                   !write(6,*) protcoords(m,l)%x, protcoords(f,g)%x
                   !if((protcoords(m,l)%x == protcoords(f,g)%x)  .and. (protcoords(m,l)%y == protcoords(f,g)%y ) &
                   !     .and. (protcoords(m,l)%z == protcoords(f,g)%z)) then
                   !   write(6,*) protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z
                   !   write(6,*) protcoords(f,g)%x,protcoords(f,g)%y,protcoords(f,g)%z
                   !   write(6,*) 'fail is due to', m,l,f,g
                   !   finalfail = .true.
                   !   fail = .true.
                   !   write(6,*) 'FAILLLLLLLLLLLLLLLLLL',time
                   !   exit
                   !end if
                !end if
             !end do
          !end do

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
end if


       end do
    end do

   

    
         

    if(finalfail .eqv. .true.) then
       write(6,*) time
       
       stop 9
    end if

  end subroutine debug





  
  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn
    !integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize,bdir2,bdir3
    integer,dimension(:,:),allocatable::cb,conn,bound
    integer,dimension(:),allocatable::maxcluslist,concount
    double precision::maxclusenergy
    allocate(cb(nprotein,nprotein))
    !allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(conn(nprotein,maxlength))
    allocate(concount(nprotein))
allocate(bound(nprotein,maxlength))    
!allocate(mashist(10))
    maxclusenergy = 0.0d0
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
          conn(m,a)=0
bound(m,a)= 0
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
  if (g<=nscaffold) call adjacent(g,l,f,tempcoord,adjver,bdir)            
!   if (g<=nscaffold) call oldadjacent(g,l,f,tempcoord,adjver,bdir)              
 !               if(g==nscaffold)  call adjacent(g,l,f,tempcoord,adjver,bdir)   
		 if(adjver.eqv. .true.) then
SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
CASE (1)
bdir2=2
bdir3=1
CASE (-2)
bdir2=3
bdir3=4
CASE (2)
bdir2=4
bdir3=3
CASE(-3)
bdir2=5
bdir3=6
CASE(3)
bdir2=6
bdir3=5
END SELECT

                      if(isbond .eqv. .true.) then
                         cb(m,g) = 1
                         cb(g,m) = 1
                         if(interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0) then

                            if(any(conn(m,:) .eq. g)) then
                               goto 63
                            else
                               concount(m) = concount(m) + 1
                               conn(m,concount(m)) = g
                            end if
63                          continue

                            if(any(conn(g,:) .eq. m)) then
                               goto 83
                            else
                               concount(g) = concount(g) + 1
                               conn(g,concount(g)) = m
                            end if
83                          continue
                            
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
                      if((protcoords(m,l)%a111==bdir2) .and. (protcoords(m,l)%a11n==g) .and. &
                        (protcoords(g,f)%a111>0) .and. (protcoords(g,f)%a11n==m) .and. &
                           (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0)) then
if((bound(m,l) ==1) .or. (bound(g,f) ==1))  goto 1111
bound(m,l) = 1
bound(g,f) = 1
                     if(any(conn(m,:) .eq. g)) then
                            goto 65
                         else
                            concount(m) = concount(m) + 1
                            conn(m,concount(m)) = g
                         end if
65                       continue

                         if(any(conn(g,:) .eq. m)) then
                            goto 85
                         else
                            concount(g) = concount(g) + 1
                            conn(g,concount(g)) = m
                         end if
85                       continue
                         
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
                      1111 continue
end if
                   end if
                end do
             end do
          !end if
       end do
    end do


    do m = 1,nprotein
       if(clnos(m) /=m) then
          clusterpop = clusterpop + 1
          if(ANY(cllist(:) .eq. clnos(m))) then
             continue
          else
             clustcount = clustcount + 1
             cllist(a) = clnos(m)
             a = a+1
          end if
       end if
    end do


    if(clustcount /= 0 .and.  (time > equilib)) then
       avecluspop = avecluspop + clusterpop
       aveclusnumber = aveclusnumber + clustcount
    end if

    clusterpop = clusterpop + clustcount

    do g = 1,nprotein
       do m = g,nprotein
          if(clnos(m) ==g) then
             histcl(g) = histcl(g) +1
          end if
       end do
       !if(time>10000) write(6,*) 'hist', histcl(g)
       normalisedsize =  normalisedsize + ((histcl(g))**2)
    end do


    !maxclus = maxval(histcl)


    maxclus = 0
    do g = 1,nprotein
       if(histcl(g) > maxclus) then
          maxclus = histcl(g)
          content = g
       end if
    end do
allocate(maxcluslist(maxclus))

z = 1
atconn= 0
do m = 1,nprotein
   if(protcoords(m,1)%species==1) atconn = atconn+concount(m)
       if(clnos(m) == content)then
          maxcluslist(z) = m
          z = z+1
       end if
    end do

    cconn =0
    acconn =0
    
    do m = 1,maxclus   
       g = maxcluslist(m)
       cconn = cconn + concount(g)
       if(protcoords(g,1)%species ==1) acconn = acconn+concount(g)
       maxclusenergy = maxclusenergy + chen(g)
    end do

    tconn=sum(concount(:))
    
    if((maxclus >10) .and. (suffclust .eqv. .false.)) then
       suffclust = .true.
       call pdbsetupalt(content,clnos)
    end if

    magiccontent = content

    !write(6,*) 'CONTENTTTTTTTTT',content
    !write(6,*) 'maxcl',maxclus
    !do g = 1,10,1
    !do m = 1,nprotein
    !if((histcl(m) > (((g-1)/10.0)*(maxclus))) .and. (histcl(m) <= ((g/10.0)*maxclus))) then
    !mashist(g) = mashist(g) + histcl(m)**2
    !end if
    !end do
    !write(38,*) time, (g/10.0)*(maxclus),mashist(g)
    !end do

   ! write(6,*) 'tconn',tconn,'atconn',atconn,'acconn',acconn,'cconn',cconn,maxclus
onepop = 0
do m = 1,maxclus
   g = maxcluslist(m)
if(protcoords(g,1)%species ==1) onepop = onepop+1
end do

    !do m = 1,nprotein
    !   write(38,*) time,histcl(m),histcl(m)**2
    !end do


    !if(noinfo .eqv. .false.) then
    !if(clustcount /= 0) write(6,*) time,clustcount, 'cluster!'
    if(clustcount /= 0) write(82,*) time,clustcount, clusterpop, real(clusterpop)/clustcount &
         , real((normalisedsize) +(nprotein-clusterpop))/nprotein,maxclus,onepop,maxclus-onepop,maxclusenergy&
,sum(chen)-maxclusenergy,acconn,cconn-acconn,atconn-acconn,tconn-((atconn-acconn)+cconn) 

    !,clcount
    !if(time>10000) write(6,*) 'output',time, real((normalisedsize)+(nprotein-clusterpop))/nprotein,real((normalisedsize)+ &
    !(nprotein-clusterpop))
    if(clustcount == 0)  write(82,*) time,clustcount, clusterpop,0.0,0.0,0.0,0,0,0,0,0,0,0,0
    !end if
    !write(6,*) 'CONTENTTTTTTTTT',content

    call phase(content,maxclus,clnos,cb)
    deallocate(cllist)
    !deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
    !write(6,*) 'CONTENTTTTTTTTT',content
  end subroutine clustercount


  subroutine clustercom(m,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge)
          integer,intent(in)::m
          integer,intent(inout)::clsize,totcluspop
          integer,dimension(:),intent(inout) ::cllist,rcounter,checklist
          integer,dimension(:,:),intent(inout) ::cb,route
          integer ::  a,b,base,maxl,maxb,ra,dint,centralchain,f,g,cench,rc,zx,maxltest
          double precision :: delx,dely,delz
          double precision,dimension(:,:),allocatable::comx,comy,comz
          logical,intent(inout) :: toolarge

          logical :: compass
          type(basicp),intent(inout) :: comcluster
          type(rprot) :: dummycom,dpcomcluster

          allocate(comx(nprotein,nprotein))
          allocate(comy(nprotein,nprotein))
          allocate(comz(nprotein,nprotein))

          t = time+1

          !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          compass = .false.
          comcluster%x = 0
          comcluster%y = 0
          comcluster%z = 0


          maxl = chlen(m)
          totcluspop = maxl 

          comx = 0.0d0
          comy = 0.0d0
          comz = 0.0d0
          dummycom%x = 0
          dummycom%y = 0
          dummycom%z = 0
          if(clsize ==1) goto 76


          do dint =1,nprotein
             rcounter(dint) = 0
          end do

          do b = 1,clsize
             do dint = 1,nprotein
                route(dint,b) = 0
             end do
             checklist(b) = 0
          end do

          totcluspop = 0

          !write(6,*) 'backup 1',rcounter(4)
          do a = 1,clsize-1,1
             b = cllist(a)
             !write(6,*) 'a'
             !if(clsize>1) write(6,*) 'component',b,com(1,b)%x,com(1,b)%y,com(1,b)%z
             do g = a+1,clsize
                f = cllist(g)
                comx(b,f) = 0.0d0
                comy(b,f) = 0.0d0
                comz(b,f) = 0.0d0
                comx(f,b) = 0.0d0
                comy(f,b) = 0.0d0
                comz(f,b) = 0.0d0
                if(cb(b,f)==1) then
                   maxb = chlen(b)
                   maxltest = chlen(f)
                   delx = com(1,f)%x - com(1,b)%x
                   if (delx > (gridsize)/2) delx = delx-gridsize
                   if (delx < -(gridsize)/2) delx = delx + gridsize
                   comx(b,f) =  delx
                   comx(f,b) = -delx
                   dely = com(1,f)%y - com(1,b)%y
                   if (dely > real(gridsize)/2) dely = dely-gridsize
                   if (dely < -real(gridsize)/2) dely = dely + gridsize
                   comy(b,f) =  dely
                   comy(f,b) = -dely
                   delz = com(1,f)%z - com(1,b)%z
                   if (delz > real(gridsize)/2) delz = delz-gridsize
                   if (delz < -real(gridsize)/2) delz = delz + gridsize
                   comz(b,f) =  delz
                   comz(f,b) = -delz
                end if
             end do
          end do

          centralchain = m
          cench = m
          checklist(1) = m
          zx = 2


          do dint = 1,clsize
             cench = checklist(dint)
             !write(6,*) 'checklist',checklist(dint)

             maxl = chlen(dint)
             totcluspop = totcluspop + maxl

             do b = 1,clsize
                a = cllist(b)
                if(a==cench) goto 16
                if(cb(cench,a) == 1) then
                   if(ANY(checklist .eq. a)) goto 16
                   rcounter(a) = rcounter(a) + 1
                   route(a,rcounter(a)) = a

                   do rc = rcounter(a),1,-1
                      route(a,rc+rcounter(cench)) = route(a,rc)
                   end do

                   do rc = 1,rcounter(cench),1
                      route(a,rc) = route(cench,rc)
                   end do
                   rcounter(a) = rcounter(a) + rcounter(cench)
                   !write(6,*) 'rcounts',rcounter(b),rcounter(dummyint),a,cench
                   checklist(zx) = a
                   zx = zx +1
                end if
16              if(a == a) continue
             end do

          end do


          !write(6,*) 'backup 2',rcounter(4)
          !write(6,*) 'orig central chain',centralchain,toolarge
          call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
          !write(6,*) 'key'


          !write(6,*) 'backup 3',rcounter(4)

76        if(clsize == 1) then
             maxl = chlen(m)
             centralchain = m
             totcluspop = maxl
             ! write(6,*) 'individual chain',com(1,m)%x,com(1,m)%y,com(1,m)%z
          end if

          !if(clsize>1) write(6,*) 'totcluspop',totcluspop
          !write(6,*) centralchain
          !write(6,*) 'runtime', com(1,centralchain)%x
          !write(6,*) dummycom%x
          comcluster%x = INT(modulo((com(1,centralchain)%x) + (dummycom%x/totcluspop)-1,real(gridsize))+1)
          comcluster%y = INT(modulo((com(1,centralchain)%y) + (dummycom%y/totcluspop)-1,real(gridsize))+1)
          comcluster%z = INT(modulo((com(1,centralchain)%z) + (dummycom%z/totcluspop)-1,real(gridsize))+1)


          dpcomcluster%x = modulo((com(1,centralchain)%x) + (dummycom%x/totcluspop)-1,real(gridsize))+1
          dpcomcluster%y = modulo((com(1,centralchain)%y) + (dummycom%y/totcluspop)-1,real(gridsize))+1
          dpcomcluster%z = modulo((com(1,centralchain)%z) + (dummycom%z/totcluspop)-1,real(gridsize))+1
          write(3,*) time,dpcomcluster%x,dpcomcluster%y,dpcomcluster%z

          if(settingup .eqv. .true.) then
             dpcomclusterinit%x = dpcomcluster%x
             dpcomclusterinit%y = dpcomcluster%y
             dpcomclusterinit%z = dpcomcluster%z
             end if
53        if(compass .eqv. .true.) continue


comcluster%x = int(real(gridsize)/2)
comcluster%y = int(real(gridsize)/2)
comcluster%z = int(real(gridsize)/2)

        end subroutine clustercom



        subroutine clustercomcheck(m,clsize,comcluster,cb,cllist)
          integer,intent(in)::m,clsize
          integer,dimension(:),intent(in) ::cllist
          integer,dimension(:,:),intent(in) ::cb
          integer ::  a,b,base,maxl,totcluspop,maxb,ra
          double precision :: delx,dely,delz,comx,comy,comz
          logical :: compass
          type(basicp),intent(inout) :: comcluster
          type(rprot) :: dummycom
          t = time+1
          !write(6,*) 'com start'
          compass = .false.
          comcluster%x = 0
          comcluster%y = 0
          comcluster%z = 0
          !clsize = 1

          maxl = chlen(m)
          totcluspop = maxl 

          comx = 0.0d0
          comy = 0.0d0
          comz = 0.0d0

          do ra = 1,clsize
             a = cllist(ra)

             if(a/=m) then
                maxb = chlen(a)
                delx = com(1,a)%x - com(1,m)%x
                if (delx > real(gridsize)/2) delx = delx-gridsize
                if (delx < -real(gridsize)/2) delx = delx + gridsize
                comx = comx + (delx*maxb)
                dely = com(1,a)%y -com(1,m)%y
                if (dely > real(gridsize)/2) dely = dely-gridsize
                if (dely < -real(gridsize)/2) dely = dely + gridsize
                comy = comy + (dely*maxb)
                delz = com(1,a)%z-com(1,m)%z
                if (delz > real(gridsize)/2) delz = delz-gridsize
                if (delz < -real(gridsize)/2) delz = delz + gridsize
                comz = comz + (delz*maxb)
                totcluspop = totcluspop + maxb
             end if
          end do


          dummycom%x = modulo((com(1,m)%x) + (comx/totcluspop)-1,real(gridsize))+1
          comcluster%x = INT(dummycom%x)
          dummycom%y = modulo((com(1,m)%y) + (comy/totcluspop)-1,real(gridsize))+1
          comcluster%y = INT(dummycom%y)
          dummycom%z = modulo((com(1,m)%z) + (comz/totcluspop)-1,real(gridsize))+1
          comcluster%z = INT(dummycom%z)


        end subroutine clustercomcheck

  
  subroutine pdbsetup
    integer::m,l,maxl,acount,dumres
    character(len=8) :: atnum
    acount = 0
    write(67,*) 'atom ', 'default ', 'radius ', 1.00000, 'name ','C'
    do m= 1,nprotein
       maxl = chlen(m)
       !if(protcoords(m,1)%species == 1) dumres = 1
       !if(protcoords(m,1)%species == 2) dumres = 2
       !if(dumres == 1) then
       do l = 1,maxl
          write(67,*) 'atom', acount, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species), &
               'resid ',protcoords(m,l)%type
          acount = acount+1
       end do
       !else if(dumres == 2) then
       !do l = 1,maxl
       !write(67,*) 'atom', acount, 'radius', 1.00000, 'name ','O ', 'resid ',dumres
       !write(67,*) 'atom', acount+1, 'radius', 1.00000, 'name ','N ','resid ',dumres
       !acount = acount+2
       !end do
       !end if

    end do
    write(67,*) ' '
    acount = 0
    do m= 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl-1
          write(atnum,'(i7)') acount
          write(67,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+1
       end do
       acount = acount +1
    end do

  end subroutine pdbsetup


  
        subroutine comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
          integer,intent(in)::clsize
          integer,intent(inout),dimension(:,:)::route
          integer,intent(inout),dimension(:)::checklist,rcounter
          integer :: zna,a,b,rc,maxb
          double precision,intent(in),dimension(:,:) ::comx,comy,comz
          type(rprot),intent(inout) :: dummycom
          logical,intent(inout)::toolarge
          double precision :: bx,by,bz,maxx,maxy,maxz,minx,miny,minz
          !write(6,*) 'centralchain',checklist(1)


          minx = 0.0d0
          miny = 0.0d0
          minz = 0.0d0
          maxx = 0.0d0
          maxy = 0.0d0
          maxz = 0.0d0

          dummycom%x = 0.0d0
          dummycom%y = 0.0d0
          dummycom%z = 0.0d0

          do zna = 2,clsize

             bx = 0.0d0
             by = 0.0d0
             bz = 0.0d0
             a= checklist(zna)
             !write(6,*) 'checkking',a,zna,rcounter(zna)
             maxb = chlen(a)
             !write(6,*) 'route', route(zna,1),a
             !write(6,*) 'finalcoms', comx(centralchain,route(zna,1)),rcounter(zna)

             bx = comx(checklist(1),route(a,1))
             by =comy(checklist(1),route(a,1))
             bz = comz(checklist(1),route(a,1))

             ! write(6,*) 'dummycom', dummycom%x,dummycom%y,dummycom%z
             !write(6,*) 'x',comx(checklist(1),route(a,1)),checklist(1),route(a,1),&
             !  com(1,checklist(1))%x,com(1,route(a,1))%x
             ! write(6,*) 'y',comy(checklist(1),route(a,1)),checklist(1),route(a,1),&
             !      com(1,checklist(1))%y,com(1,route(a,1))%y
             !write(6,*) 'z',comz(checklist(1),route(a,1)),checklist(1),route(a,1),&
             !    com(1,checklist(1))%z,com(1,route(a,1))%z


             if(rcounter(a)>1) then
                do rc = 1,rcounter(a)-1,1
                   !write(6,*) 'finalcoms', comz(centralchain,route(a,rc))
                   !write(6,*) 'route', route(zna,rc)
                   !write(6,*) 'total count =',rcounter(zna)
                   bx  = bx + comx(route(a,rc),route(a,rc+1))
                   by = by + comy(route(a,rc),route(a,rc+1))
                   bz = bz + comz(route(a,rc),route(a,rc+1))
                   !write(6,*) 'x',comx(route(a,rc),route(a,rc+1))
                   !write(6,*) 'y',comy(route(a,rc),route(a,rc+1))
                   !write(6,*) 'z',comz(route(a,rc),route(a,rc+1))
                   !write(6,*) 'dummycom', dummycom%x,dummycom%y,dummycom%z
                end do
             end if
             !write(6,*) 'finito'

             if(bx > maxx) maxx = bx
             if(by > maxy) maxy = by
             if(bz > maxz) maxz = bz
             if(bx < minx) minx = bx
             if(by < miny) miny = by
             if(bz < minz) minz = bz

             dummycom%x = dummycom%x + (maxb*bx)
             dummycom%y = dummycom%y + (maxb*by)
             dummycom%z = dummycom%z + (maxb*bz)

          end do


          if((maxx - minx) > real(gridsize)/2) toolarge = .true.
          if((maxy - miny) > real(gridsize)/2) toolarge = .true.
          if((maxz - minz) > real(gridsize)/2) toolarge = .true.


          write(161,*) maxx-minx,maxy-miny,maxz-minz
          !do a = 1,clsize
          !b= checklist(a)
          !write(6,*) 'individual',b,com(1,b)%x,com(1,b)%y,com(1,b)%z
          ! end do

        end subroutine comlargecomp

  subroutine pdbsetupalt(g,clnos)
    integer,intent(in)::g
    integer,dimension(:),intent(in):: clnos
    integer::m,l,maxl,acount,dumres,infoout
    character(len=8) :: atnum
    acount = 0
    write(88,*) 'atom ', 'default ', 'radius ', 1.00000, 'name ','C'
    do m= 1,nprotein
       maxl = chlen(m)
       !if(protcoords(m,1)%species == 1) dumres = 1
       !if(protcoords(m,1)%species == 2) dumres = 2
       !if(dumres == 1) then
       if(clnos(m) == clnos(g)) then
          infoout = 100
       else
          infoout = protcoords(m,1)%species
       end if

       
       do l = 1,maxl
          write(88,*) 'atom', acount, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species), &
               'resid ',protcoords(m,l)%type
          acount = acount+1
       end do
       !else if(dumres == 2) then
       !do l = 1,maxl
       !write(67,*) 'atom', acount, 'radius', 1.00000, 'name ','O ', 'resid ',dumres
       !write(67,*) 'atom', acount+1, 'radius', 1.00000, 'name ','N ','resid ',dumres
       !acount = acount+2
       !end do
       !end if

    end do
    write(67,*) ' '
    acount = 0
    do m= 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl-1
          write(atnum,'(i7)') acount
          write(88,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+1
       end do
       acount = acount +1
    end do

  end subroutine pdbsetupalt








  subroutine freesitescalc
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn
    integer,dimension(:),allocatable:: countoffree
    logical :: clusyes,adjver,newcl
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize,maxl,bdir2,bdir3
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
SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
CASE (1)
bdir2=2
bdir3=1
CASE (-2)
bdir2=3
bdir3=4
CASE (2)
bdir2=4
bdir3=3
CASE(-3)
bdir2=5
bdir3=6
CASE(3)
bdir2=6
bdir3=5
END SELECT
           
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
                      if((protcoords(m,l)%a111==bdir2) .and. (protcoords(m,l)%a11n==g) .and. &
                        (protcoords(g,f)%a111>0) .and. (protcoords(g,f)%a11n==m) .and. &
                           (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)< 0.0)) then
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
countoffree(protcoords(g,l)%type+mtype) = countoffree(protcoords(g,l)%type+mtype) + 1
end if
end do
       end if
    end do
    
    write(6,*) 'unbound species in left and right',countoffree(:)


  
   
    magiccontent = content

  
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
          dx = min(abs(protcoords(m,l)%x - (gridsize/2.0)), gridsize - abs(protcoords(m,l)%x - (gridsize/2.0)))
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
    nprotein = 1                 !number of chains
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


          allocate(interen(mtype,mtype)) !-nclientspecies,mtype-nclientspecies))
          allocate(interenergy(mtype,mtype))!-nclientspecies,mtype-nclientspecies))
          allocate(intraenergy(mtype,mtype))!-nclientspecies,mtype-nclientspecies))
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
          do f = 1,mtype!-nclientspecies
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

CASE ('PATH')
CALL get_string(pathtofile)

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



