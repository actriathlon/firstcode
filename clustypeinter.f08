program move

  implicit none

  double precision :: totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength
  double precision :: actualchainlength,random
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug,testdummy,polymerrog,polymerl
  integer :: piv,right,endm,crank,rept,datayes,maxclussize,cranklimit,minl,outputrate
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count,avecluspop,aveclusnumber
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib,cltacc,cltatt,clracc,clratt
  integer :: successful,reject,pivotlimit,maxpiv,maxlength1,maxlength2,q,d,u,sm,qt,ft,gt,scans
  integer:: cracc,crrej,raacc,rarej,rerej,reacc,pvacc,pvrej,nprotein1,nprotein2,nspecies,mtype
  integer,dimension(:,:),allocatable :: bonddd,bdi
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen
  double precision :: crrat,emrat,rarat,pvrat,rerat,trackx,tracky,trackz,intraen,intdummy,totdire
  double precision,dimension(:,:),allocatable :: interenergy,interen,tttt  
  double precision,dimension(:),allocatable :: intraenergy
  integer,dimension(:,:,:), allocatable :: isobond,ba
  integer,dimension(:,:,:,:),allocatable :: intbond,br
  real, external :: ran2
  real,dimension(:),allocatable :: variance,disp
  logical :: exist,fail,finalfail,debugyes,film,isbond,clt,clr,noinfo,scalinginfo,nobonds
  real::start,finish


  !change COM code for large protein in small box to stop pbc effect - similar to old version - but still
  !use dispacement from central bead

  !cluster COM may fall foul of pbcs

  !Sort MSD outputs
  
  type protein
     integer :: x,y,z,species,type
  end type protein

  type rprot
     real :: x,y,z
  end type rprot

  type basicp
     integer :: x,y,z
  end type basicp

  type centremass
     double precision :: x,y,z
  end type centremass


  type(protein),dimension(:,:),allocatable :: protcoords
  type(centremass),dimension(:,:),allocatable :: com
 
mtype = 0
  !noinfo = .false.
  call read_setup

  
  nprotein = sum(speciespop)
  write(6,*) 'nprotein = ', nprotein
  maxlength = maxval(specieslen)
  write(6,*) 'maxlength =',maxlength
  !open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')
write(6,*) 'logicals',scalinginfo,noinfo
  if(noinfo .eqv. .true.) then
     open(11,file='scalingdata.dat',action = 'write')
       open(82,file = 'clusterdata.dat', action = 'write')
    else if(noinfo .eqv. .false.) then
       if(scalinginfo .eqv. .true.) then
          open(11,file='scalingdata.dat',action = 'write')
          open(49,file = 'structure.dat',action = 'write')
  open(29, file = 'rms.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  open(79, file = 'runningave.dat', action = 'write')
  open(93, file = 'energy.dat', action = 'write')
  open(13,file='timedata.dat',action = 'write')
  open(19,file='acceptance.dat',action  = 'write')
else if(scalinginfo .eqv. .false.) then
   open(49,file = 'structure.dat',action = 'write')
     open(29, file = 'rms.dat', action = 'write')
  open(93, file = 'energy.dat', action = 'write')
  open(82,file = 'clusterdata.dat', action = 'write')
  open(13,file='timedata.dat',action = 'write')
  open(19,file='acceptance.dat',action  = 'write')
    open(38,file = 'histclust.dat',action = 'write')
  !open(67,file='bondddd.dat',action  = 'write')
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
  trackx = 0.0d0
  tracky = 0.0d0
  trackz = 0.0d0
  polymerl = 0.0d0
  polymerrog = 0.0d0
  outputrate = maxtime/100
  clracc = 0
  clratt = 0
  cltacc = 0
  cltatt = 0
  !kT = 50.0
  successful = 0
  reject = 0

  
totdire = 0.0d0
  do ft = 1,mtype
     do gt = 1,mtype
        totdire = totdire + interenergy(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
     end do
 
if(totdire == 0.0d0) nobonds =  .true.



  write(6,*) maxlength1,maxlength2
  write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT
  write(6,*) 'clr',clr,'clt',clt

  open(67, file = 'move.vtf', action = 'write')
  allocate(protcoords(nprotein,maxlength))
  allocate(com(2,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(isobond(nprotein,maxlength,maxlength))
  allocate(intbond(nprotein,nprotein,maxlength,maxlength))
  allocate(ba(nprotein,maxlength,maxlength))
  allocate(br(nprotein,nprotein,maxlength,maxlength))
  allocate(bonddd(nprotein,maxlength))
  allocate(bdi(nprotein,maxlength-1))
  allocate(chlen(nprotein))

!write(6,*) 'modulotest', modulo(11,10)
  N = maxlength
  !interenergy(1,1) =  tttt(1,1) 
  !interenergy(1,2) =  tttt(1,2)
  !interenergy(2,1) =  tttt(2,1) 
  !interenergy(2,2) =  tttt(2,2) 
  intraenergy(1) = -0.0  
  intraenergy(2) = -0.0
  intraenergy(3) = -0.0  
  !write(6,*) 'direct energy',interenergy(1,1),interenergy(1,2),interenergy(2,1),interenergy(2,2),&
       !interenergy(1,3),interenergy(2,3),interenergy(3,3)
  !interen = -2.0
  !intraen = -1.0
  count = 0
  time = 0
  totdisp = 0.0d0

  if(abs(interen(1,1)) + abs(interen(2,2))+ abs(interen(1,2)) /= 0.0) then
     isbond = .true.
  else
     isbond = .false.
  end if
  WRITE(6,*) 'ISBOND =',isbond
  write(6,*) 'interen',interen

  cracc = 0
  crrej = 0
  rerej = 0
  reacc = 0
  raacc = 0
  rarej = 0
  pvacc = 0
  pvrej = 0


  !sets the limit of the pivot move

  do sm = 1,nprotein
     do  u = 1,maxlength
        do d = 1,maxlength
           ba(sm,u,d) = 0
        end do
        do q = 1,nprotein
           do d = 1,maxlength
              br(sm,q,u,d) = 0
           end do
        end do
     end do
  end do
     call CPU_TIME(start)

  call foundation

  if(scalinginfo .eqv. .false.) call pdbsetup
  minl = minval(specieslen)
  pivotlimit = min((minl/2),maxpiv)

  write(6,*) 'a'
  !call dataout
  do qt = 1,nprotein
     call comfind(qt,.false.)
  end do
  write(6,*) 'b'

  actualchainlength = 0.0
  do time = 1,maxtime,1
!write(6,*) 'c'
     fail = .false.
         !if(time > 288000 .and. (time < 290000)) then
     if ((mod(time,outputrate) == 0) .and. (film .eqv. .true.) ) then

     !if(time>5300)
     call dataout
     end if
     !end if
     if((nobonds .eqv. .true.) .and. (clt .eqv. .false.) .and. (clr .eqv. .false.)) then
!do scans = 1,1000
   call positioning
 !  end do
else
  ! do scans = 1,1000
      call pickmoves
   !   end do
end if
     if (mod(time,outputrate) == 0) then
        if(noinfo .eqv. .false.)  write(93,*) time,totalenergy
        !call energy
        !call bondcheck
     end if
     if(modulo(time,outputrate) == 0 .and. (scalinginfo .eqv. .false.)) call clustercount
!call clustercount
     if((time > equilib) .and. (modulo(time,outputrate) == 0) .and. (scalinginfo .eqv. .true.)) then
        call length
        call radiusofgyration
        if((noinfo .eqv. .false.) .and. (scalinginfo .eqv. .true.)) write(79,*) time-equilib,&
             runningaveEtE/(time-equilib), runningaveROG/(time-equilib)
     end if
     
     !write(6,*) 'gfd'
     if (modulo(time,outputrate) == 0 .and. debugyes .eqv. .true.) then
        !write(6,*) 'debug start'
        call debug
        !write(6,*) 'debug fin'
        end if
        if(modulo(time,100000) ==0) call energy
!if((time> 13600) .and. (modulo(time,20) == 0)) call energy
        !write(6,*) 'debug success'
     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else if(modulo(time,outputrate) == 0) then
        write(6,*) 'step =',time
     end if
    
     if((modulo(time,maxtime/10) == 0) .and.(noinfo .eqv. .false.)) then
        write(19,*) 'pivot',pvacc,pvrej, real(pvacc)/(pvacc+pvrej)
        write(19,*) 'reptation',reacc,rerej, real(reacc)/(reacc+rerej)
        write(19,*) 'right angle',raacc,rarej, real(raacc)/(raacc+rarej)
        write(19,*) 'crankshaft',cracc,crrej, real(cracc)/(cracc+crrej)
        write(19,*) 'cluster rotate',clracc,clratt,real(clracc)/clratt
        write(19,*) 'cluster translation',cltacc,cltatt,real(cltacc)/cltatt
     end if

     !write(6,*) 'sweep complete'
  end do

  if((film .eqv. .false.) .and. (scalinginfo .eqv. .false.)) call dataout
  write(6,*) 'deltax =',trackx, 'deltay =',tracky, 'deltaz =',trackz
   call CPU_TIME(finish)
      if(noinfo .eqv. .false.) write(13,*) (finish-start)
  !call error
  call energy
  if((noinfo .eqv. .true.) .or. (scalinginfo .eqv. .true.)) write(11,*) maxlength1, &
       (polymerl/((maxtime-equilib)/outputrate)), &
      (polymerrog/((maxtime-equilib)/outputrate))
  write(6,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))
  write(6,*) 'average clusper pop', (real(avecluspop)/aveclusnumber)/((time-equilib)/outputrate), &
     real(avecluspop)/((time-equilib)/outputrate), real(aveclusnumber)/((time-equilib)/outputrate)
  !inquire(file = "lengthdataaverages.dat", exist = exist)
  !if (exist) then
  ! open(30, file = "lengthdataaverages.dat", status = "old", position = "append",action ="write")
  !else
  !open(30, file = "lengthdataaverages.dat", status = "new",action ="write")
  !end if
  !write(30,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))!,2*(SQRT(totvar)/maxtime)/(actualchainlength/maxtime)
  write(6,*) 'reject = ', reject
  write(6,*) 'successful =', successful
  !write(6,*) 1.0*successful/(reject + successful)
  write(6,*) 'acceptance ratio =', 1.0*successful/(reject + successful)
  write(6,*) 'count =', count
  write(6,*) 'maxlength =', maxlength
  write(6,*) 'forward =', reptforward
  write(6,*) 'backward =' , reptbackward
write(6,*) 'direct energy',interenergy(1,1),interenergy(1,2),interenergy(2,1),interenergy(2,2)
  if (finalfail .eqv. .true.) then
     write(6,*) '***********^^^^^^^^^^^^^^^^^^^This simulation FAILED!^^^^^^^^^^***********'
  end if
contains


 
subroutine vectorspin(m,l,maxl)
    integer,intent(in)::maxl,m,l
    integer::p,f,dum1,dum2,no1,no2,maxt,g
    integer,dimension(:),allocatable:: tbo,tbdi
    type(basicp),dimension(:),allocatable::tempcoord
    double precision :: deltaenergy
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:), allocatable :: tempintbond,tbr
    allocate(tempisobond(maxlength,maxl))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tbr(nprotein,maxlength,maxlength))
    allocate(tbo(maxl))
    allocate(tbdi(maxl-1))
    allocate(tempcoord(maxl))

    !write(6,*) maxl,maxlength1,maxlength2,protcoords(m,1)%species
    !do a = 1,nprotein*maxlength

    !do l = p,f
    !write(6,*) 'vector' ! totalenergy
    if(l==maxl) then
       call bonddirection(5,7,-1*bdi(m,l-1),tbo(l))
    tempcoord(l)%x = protcoords(m,l)%x
    tempcoord(l)%y = protcoords(m,l)%y
    tempcoord(l)%z = protcoords(m,l)%z
       else
    if(l == 1) call bonddirection(5,bdi(m,l),7,tbo(l))
    if((l >1) .and. (l <maxl)) call bonddirection(4,bdi(m,l),-1*bdi(m,l-1),tbo(l))
    tempcoord(l)%x = protcoords(m,l)%x
    tempcoord(l)%y = protcoords(m,l)%y
    tempcoord(l)%z = protcoords(m,l)%z
    tbdi(l) = bdi(m,l)
end if

    deltaenergy = 0.0d0

    call intrabondvector(m,l,tempcoord,tempisobond,deltaenergy,tbo,maxl)
    !write(6,*) 'up to here'
    call intervector(m,l,l,tempcoord,tempintbond,tbr,deltaenergy,tbo)
    !error here
    !write(6,*) 'survived to here'

    if(Energydecision(deltaenergy) .eqv. .True.) then
       totalenergy = totalenergy + deltaenergy

       call updatepos(m,l,l,tempcoord,tbo,tbdi,maxl)
       call updateintrabondvector(m,l,tempisobond,maxl)
       call updateinterbond(m,l,l,tempintbond,tbr)
    else 
       continue
    end if
    !write(6,*) 'finished vector'
  end subroutine vectorspin

  subroutine updateintrabondvector(chainnum,beadnum,tempisobond,maxl)
    integer,intent(in):: chainnum,maxl,beadnum
    integer ::g
    integer,dimension(:,:),intent(in) :: tempisobond

       do g = 1,maxl
          if(g<beadnum-2) isobond(chainnum,g,beadnum) = tempisobond(g,beadnum)
          if(g> beadnum+2) isobond(chainnum,beadnum,g) = tempisobond(beadnum,g)
       end do

  end subroutine updateintrabondvector

  subroutine intervector(chain1,beadmin,beadmax,tempc,tempintbond,tbr,deltaenergy,tbo)

    integer,intent(in)::chain1,beadmin,beadmax
    integer,dimension(:),intent(in)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:,:),intent(inout) :: tempintbond,tbr
    integer::bead1,bead2,chain2,bdir,maxback,f
    Type(basicp),dimension(:),intent(in) :: tempc
    logical :: adjver


    if(nprotein == 1) return

    do chain2 = 1,chain1-1,1
       do bead1 = beadmin,beadmax,1

        maxback = chlen(chain2)

          do bead2 = 1,maxback,1
             if(intbond(chain2,chain1,bead2,bead1) == 1) then
                deltaenergy = deltaenergy -interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             end if
             if(br(chain2,chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - &
                  interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             call adjacent(chain2,bead1,bead2,tempc,adjver,bdir)
             if((adjver .eqv. .true.)) then
                deltaenergy = deltaenergy + &
                     interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                tbr(chain2,bead1,bead2) = 1
                if((bdir == -1*tbo(bead1)).and. (bdir ==bonddd(chain2,bead2))) then
                   deltaenergy = deltaenergy + &
                        interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                   tempintbond(chain2,bead1,bead2) = 1
                else
                   tempintbond(chain2,bead1,bead2) = 0
                end if
             else if((adjver .eqv. .false.)) then
                tbr(chain2,bead1,bead2) = 0
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do


    do chain2 = chain1+1,nprotein,1
       do bead1 = beadmin,beadmax,1
        maxback = chlen(chain2)
          do bead2 = 1,maxback,1
             if(intbond(chain1,chain2,bead1,bead2) == 1) deltaenergy = deltaenergy - &
                  interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             if(br(chain1,chain2,bead1,bead2) == 1) deltaenergy = deltaenergy - &
                  interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             call adjacent(chain2,bead1,bead2,tempc,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + &
                     interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                tbr(chain2,bead1,bead2) = 1
                if((bdir == -1*tbo(bead1)) .and. (bdir == bonddd(chain2,bead2))) then
                   deltaenergy = deltaenergy + &
                        interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                   tempintbond(chain2,bead1,bead2) = 1
                else
                   tempintbond(chain2,bead1,bead2) = 0
                end if
             else if(adjver .eqv. .false.) then
                tbr(chain2,bead1,bead2) = 0
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do

  end subroutine intervector

  subroutine intrabondvector(chain1,bead1,tempcoord,tempisobond,deltaenergy,tbo,maxl)

    integer,intent(in)::chain1,maxl,bead1
    integer,dimension(:),intent(inout)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer::bead2,bdir,f
    Type(basicp),dimension(:),intent(inout) :: tempcoord
    logical :: adjver


       do bead2 = bead1+3,maxl
          if(isobond(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(chain1,bead1)%type)
          call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
          if ((adjver.eqv. .true.) .and. &
               (bdir == -1*tbo(bead1)) .and. (bdir == tbo(bead2))) then
             deltaenergy = deltaenergy + intraenergy(protcoords(chain1,bead1)%type)
             tempisobond(bead1,bead2) = 1
          else 
             tempisobond(bead1,bead2) = 0
          end if
       end do


       do bead2 = 1,bead1-3
          if(isobond(chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(chain1,bead2)%type)
          call adjacent(chain1,bead2,bead1,tempcoord,adjver,bdir)
          if ((adjver.eqv. .true.) .and. &
               (bdir == -1*tbo(bead2)) .and. (bdir == tbo(bead1))) then
             deltaenergy = deltaenergy + intraenergy(protcoords(chain1,bead2)%type)
             tempisobond(bead2,bead1) = 1
          else 
             tempisobond(bead2,bead1) = 0
          end if
       end do
       
  end subroutine intrabondvector



  subroutine foundation
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z
    character(len = 10) :: BIN
    read(23,*) BIN
    !write(6,*) maxlength

    !write(67,*) ' '
    do m = 1,nprotein,1
       read(23,*) protcoords(m,1)%species, protcoords(m,1)%x, protcoords(m,1)%y, protcoords(m,1)%z,protcoords(m,1)%type,chlen(m)

       maxl = chlen(m)

       do l = 2,maxl,1
          !write(6,*) l
          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z,protcoords(m,l)%type
       end do
    end do


    do m = 1,nprotein
       l = 1
       !bonddd(m,l) = 0

       if((protcoords(m,l+1)%x - protcoords(m,l)%x == 1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
            -(gridsize-1))) bdi(m,l) = 1 
       if((protcoords(m,l+1)%x - protcoords(m,l)%x == -1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
            (gridsize-1))) bdi(m,l) = -1
       if((protcoords(m,l+1)%y - protcoords(m,l)%y == 1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
            -(gridsize-1))) bdi(m,l) = 2
       if((protcoords(m,l+1)%y - protcoords(m,l)%y == -1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
            (gridsize-1))) bdi(m,l) = -2
       if((protcoords(m,l+1)%z - protcoords(m,l)%z == 1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
            -(gridsize-1))) bdi(m,l) = 3
       if((protcoords(m,l+1)%z - protcoords(m,l)%z == -1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
            (gridsize-1))) bdi(m,l) = -3
       !if(modulo(protcoords(m,l+1)%x - protcoords(m,l)%x,gridsize)+1 == 1) bdi(m,l) = 1 
       !if(modulo(protcoords(m,l+1)%x - protcoords(m,l)%x,gridsize)+1 == -1) bdi(m,l) = -1
       !if(modulo(protcoords(m,l+1)%y - protcoords(m,l)%y,gridsize)+1 == 1) bdi(m,l) = 2
       !if(modulo(protcoords(m,l+1)%y - protcoords(m,l)%y,gridsize)+1 == -1) bdi(m,l) = -2
       !if(modulo(protcoords(m,l+1)%z - protcoords(m,l)%z,gridsize)+1 == 1) bdi(m,l) = 3
       !if(modulo(protcoords(m,l+1)%z - protcoords(m,l)%z,gridsize)+1 == -1) bdi(m,l) = -3
       no2 = 7
       call bonddirection(5,bdi(m,l),no2,bonddd(m,l))
       !write(6,*) bonddd(m,l)
        maxl = chlen(m)
       do l =2,maxl-1
          !bonddd(m,l) = 0
          if((protcoords(m,l+1)%x - protcoords(m,l)%x == 1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
               -(gridsize-1))) bdi(m,l) = 1 
          if((protcoords(m,l+1)%x - protcoords(m,l)%x == -1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
               (gridsize-1))) bdi(m,l) = -1
          if((protcoords(m,l+1)%y - protcoords(m,l)%y == 1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
               -(gridsize-1))) bdi(m,l) = 2
          if((protcoords(m,l+1)%y - protcoords(m,l)%y == -1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
               (gridsize-1))) bdi(m,l) = -2
          if((protcoords(m,l+1)%z - protcoords(m,l)%z == 1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
               -(gridsize-1))) bdi(m,l) = 3
          if((protcoords(m,l+1)%z - protcoords(m,l)%z == -1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
               (gridsize-1))) bdi(m,l) = -3
          !if(modulo(protcoords(m,l+1)%x - protcoords(m,l)%x,gridsize)+1 == 1) bdi(m,l) = 1 
          !if(modulo(protcoords(m,l+1)%x - protcoords(m,l)%x,gridsize)+1 == -1) bdi(m,l) = -1
          ! if(modulo(protcoords(m,l+1)%y - protcoords(m,l)%y,gridsize)+1 == 1) bdi(m,l) = 2
          !if(modulo(protcoords(m,l+1)%y - protcoords(m,l)%y,gridsize)+1 == -1) bdi(m,l) = -2
          !if(modulo(protcoords(m,l+1)%z - protcoords(m,l)%z,gridsize)+1 == 1) bdi(m,l) = 3
          !if(modulo(protcoords(m,l+1)%z - protcoords(m,l)%z,gridsize)+1 == -1) bdi(m,l) = -3
          call bonddirection(4,bdi(m,l),-1*bdi(m,l-1),bonddd(m,l))
       end do
       l = maxl
       no1 = 7
       call bonddirection(5,no1,-1*bdi(m,l-1),bonddd(m,l))

       !if(bonddd(m,l) == 0) write(6,*) 'fail'
    end do

    do m =1,nprotein
        maxl = chlen(m)
       do l =1,maxl-1
          if((bdi(m,l) == 0) .or. (abs(bdi(m,l)) >3)) write(6,*) 'fail',m,l
       end do
    end do


    call energy
    call debug

  end subroutine foundation

  subroutine clusterposition
    integer::scans,m,clsize,clback,h,a
    real :: decider
    integer,dimension(:),allocatable:: clno,cllist
    integer,dimension(:,:),allocatable::cb
    allocate(clno(nprotein))
    allocate(cb(nprotein,nprotein))
    m =int(ran2(seed)*(nprotein))+1
    !call clustercount(clno)
    !write(6,*) 'cluster'


    if((clr .eqv. .false.) .and. (clt .eqv. .false.)) then
       call positioning
       return
       end if

    call clusterassign(m,clno,clsize,cb)
    allocate(cllist(clsize))
    !clback = clsize
!write(6,*) 'clsize',clsize
  if((1.0/clsize) > (0.1*ran2(seed))) then
    cllist(1) =m
    h = 2
    do a = 1,nprotein
       if(clno(a) == clno(m)) then
          if(a/=m) then
             cllist(h) = a
             h = h + 1
          end if
       end if
    end do



    if(clr .eqv. .false.) then
       cltatt = cltatt+1
       call clustertranslation(m,clno,clsize,cb,cllist)
    else if (clt .eqv. .false.) then
       clratt = clratt+1
       call clusterrotation(m,clno,clsize,cb,cllist)
    else if((clr .eqv. .true.) .and. (clt .eqv. .true.)) then
  
       
       !write(6,*) 'vects',time
       decider = ran2(seed) - 0.5
       !do scans = 1,nprotein
       !if(time > 42000 .and. (time < 47000))write(6,*) 'clustermove start'
       if(decider > 0.0) then
          !return
          cltatt = cltatt+1
          call clustertranslation(m,clno,clsize,cb,cllist)
          !if(time> 3380 .and. (time < 3390)) write(6,*) 'translate',time
       else if(decider <= 0.0) then
          clratt = clratt+1
          call clusterrotation(m,clno,clsize,cb,cllist)
          !if(time> 3380 .and. (time < 3390)) write(6,*) 'rotate',time
       end if
       !end do
        !if(time > 42000 .and. (time < 47000))write(6,*) 'clustermove fin'
       !call clusterassign(m,clno,clback,cb)

       !if(clback /= clsize) then
       !write(6,*) 'cluster fail -------------',clback,clsize,time
       !finalfail = .true.
       !end if
    end if
    end if
!write(6,*) 'cluster complete'
  end subroutine clusterposition

  subroutine pickmoves
    integer:: scans,m,maxl,decide,l,a
    !double precision :: decide
scans = 1
    !do scans = 1,(nprotein1*maxlength1 + nprotein2*maxlength2)
    decide = int(ran2(seed)*10)+1
        !if(time > 570000 .and. (time < 573000))write(6,*) 'decide',decide,time

    !write(6,*) 'start loop',decide

    if(decide <= 6) then
       do a=1,maxlength
       m =int(ran2(seed)*(nprotein))+1
        maxl = chlen(m)
       l = int(ran2(seed)*(maxl))+1
       call vectorspin(m,l,maxl)
       end do
    end if
    if(decide > 6 .and. decide <= 8) call positioning
    if(decide > 8) call clusterposition
    !if(time > 998 .and. time < 1002) write(6,*) 'information',decide,time
    !if(time == 1166) write(6,*) 'hello'
    !end do

    
        if((modulo(time,maxtime/100) == 0) .and. (time>equilib)) then
        write(29,*) (time-equilib)*scans, real(totrmsbrute/nprotein)
     end if

  end subroutine pickmoves


  subroutine rotate(m,l,b,tempcoords,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz)
    integer,intent(in):: m,l,b,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz
    type(basicp),dimension(:),intent(inout) :: tempcoords
    !write(6,*) dx,dy,dz
    tempcoords(b)%x = modulo(protcoords(m,l)%x + (xx*dx) + (xy*dy) + (xz*dz) -1,gridsize)+1
    tempcoords(b)%y = modulo(protcoords(m,l)%y + (yx*dx) + (yy*dy) + (yz*dz)-1,gridsize)+1
    tempcoords(b)%z = modulo(protcoords(m,l)%z + (zx*dx) + (zy*dy) + (zz*dz)-1,gridsize)+1

  end subroutine rotate



  subroutine positioning
    !Selects move to be carried out on the proteins
    integer :: nmoves,scan,m,l,maxl
    logical :: run,run2
    double precision:: choose
    integer ::randomz,pivlen
    !if(time > 42000 .and. (time < 47000)) write(6,*) 'start positioning'
    t = time + 1

    count = count + 1
    run = .true.
    run2 = .true.
    !return

    nmoves = piv + crank + right + rept
    randomz = int(ran2(seed)*nmoves)+1

    if (randomz == crank .and. crank == 1) then
       !if(time == 49022) write(6,*) 'crank'
       call crankshaftmove
    else if (randomz == (crank + right) .and. right ==1) then
       !if(time== 49022) write(6,*) 'right' 
       call rightanglemove
    else if (randomz == (crank + right + rept) .and. rept ==1) then
       !if(time== 49022) write(6,*) 'reptation' 
       call reptation
    else if (randomz == (crank + right + rept + piv) .and. piv == 1) then 
       m =int(ran2(seed)*(nprotein))+1

17             maxl = chlen(m)
       !write(6,*) maxl
       
       pivlen = int(ran2(seed)*(maxl/2))+1
       if(pivlen >= gridsize/2) goto 17
       
       if(ran2(seed) -0.5 > 0) then
          l = pivlen
       else
          l = maxl -  pivlen
       end if
       !if(time == 49022) write(6,*) 'pivot',m,l
       call pivot(m,l,maxl)
       call comfind(m,.TRUE.)
    end if

    !write(6,*) 'finish positioning'
  end subroutine positioning

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

  subroutine crankshaftmove
    !performs crankshaft move on section of chain around axis which the section lies upon
    logical :: crankcont,cranksep,pivotx,pivoty,pivotz
    double precision :: deltaenergy
    integer :: dummy,dummy2,p,s,m,l,str,st,pr,sign,probs,maxl,maxback
    integer,dimension(:),allocatable::dx1,dy1,dz1
    type(basicp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbo,tbdi
    integer,dimension(:,:), allocatable :: tempisobond,tba
    integer,dimension(:,:,:),allocatable :: tempintbond,tbr
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tba(maxlength,maxlength))
    allocate(tbr(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(dx1(maxlength))
    allocate(dy1(maxlength))
    allocate(dz1(maxlength))
    allocate(tbo(maxlength))
    allocate(tbdi(maxlength))
    !write(6,*) 'start crank'
    !this means that m/=0 and is up to nprotein
    m =int(ran2(seed)*(nprotein))+1

    !if(time == 236) write(6,*) 'crank' ,m
        maxl = chlen(m)
    if(maxl < 4) return
    !this means that l cannot be maxlength,maxlength - 1 or l
    cranksep = .true.

    l = int(ran2(seed)*(maxl-4))+1
    !dummy2 = int(ran2(seed)*(maxl-1))+1
    !l = min(dummy,dummy2)
    !p = max(dummy,dummy2)
    !if ((p-l < 3) .or.( p-l > cranklimit)) then
       !goto 71
    !end if
 p = l + 3
    probs = int(ran2(seed)*3) +1
    pivotx = .false.
    pivoty = .false.
    pivotz = .false.

    !check to see if 4 atoms are in the same plane
    if ((protcoords(m,l)%x /= protcoords(m,p)%x) &
         .and.  (protcoords(m,l)%y == protcoords(m,p)%y) &
         .and. (protcoords(m,l)%z == protcoords(m,p)%z)) then
       pivotx = .true.                  
    else if ((protcoords(m,l)%x == protcoords(m,p)%x) .and. &
         (protcoords(m,l)%y /= protcoords(m,p)%y) &
         .and. (protcoords(m,l)%z == protcoords(m,p)%z) ) then
       pivoty = .true.                               
    else if ((protcoords(m,l)%x == protcoords(m,p)%x) &
         .and. (protcoords(m,l)%y == protcoords(m,p)%y) &
         .and. (protcoords(m,l)%z /= protcoords(m,p)%z)) then
       pivotz = .true.
    end if


    if ((pivotx .eqv. .false.) .and. (pivoty .eqv. .false.) .and. (pivotz .eqv. .false.)) then
       crankcont = .false.
       goto 43
    end if

    crankcont = .true.

    if (pivotx .eqv. .true.) then
       !write(6,*) 'x start'
       do s= l+1,p-1,1

          dx1(s) = protcoords(m,s)%x - protcoords(m,l)%x
          if (dx1(s) > gridsize/2) dx1(s) = dx1(s) - gridsize
          if (dx1(s) < -gridsize/2) dx1(s) = dx1(s) + gridsize

          dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if(dy1(s) > gridsize/2) dy1(s) = dy1(s)-gridsize
          if(dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize


          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if(dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if(dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize

       end do

       if (probs == 1) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,1,0,0,0,0,-1,0,1,0,dx1(s),dy1(s),dz1(s))
          end do
          do s =l,p,1
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,-1)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,1,-1)
          end do

       else if (probs == 2) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,1,0,0,0,0,1,0,-1,0,dx1(s),dy1(s),dz1(s))
          end do
          do s =l,p,1
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,1)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,1,1)
          end do


       else if (probs == 3) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,1,0,0,0,-1,0,0,0,-1,dx1(s),dy1(s),dz1(s))
          end do
          do s = l,p,1
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,0)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,1,0)
          end do

       end if
       !write(6,*) 'x done'
    else if (pivoty .eqv. .true.) then
       !write(6,*) 'y start'
       do s= l+1,p-1,1

          dx1(s) = protcoords(m,s)%x - protcoords(m,l)%x
          if (dx1(s) > gridsize/2) dx1(s) = dx1(s) - gridsize
          if (dx1(s) < -gridsize/2) dx1(s) = dx1(s) + gridsize
          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if (dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if (dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize
          dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if (dy1(s) > gridsize/2) dy1(s) = dy1(s) - gridsize
          if (dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize

       end do
       !write(6,*) probs
       if (probs == 1) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,0,-1,0,1,0,1,0,0,dx1(s),dy1(s),dz1(s))
          end do
          !write(6,*) 'done rotate 1'
          do s = l,p
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,2,1)
             end if
          end do
          !write(6,*) 'done rotate 2'
          do s = l,p-1,1
             call ribond(m,s,tbdi,2,1)
          end do

          !write(6,*) 'done rotate 3'
       else if (probs ==2) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,0,1,0,1,0,-1,0,0,dx1(s),dy1(s),dz1(s))
          end do
          do s = l,p,1
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,2,-1)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,2,-1)
          end do

       else if (probs == 3) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,-1,0,0,0,1,0,0,0,-1,dx1(s),dy1(s),dz1(s))
          end do

          do s = l,p
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,2,0)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,2,0)
          end do

       end if
       !!pivot y!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
       !write(6,*) 'y done'
    else if (pivotz .eqv. .true.) then
       !write(6,*) 'z start'
       do s= l+1,p-1,1

          dx1(s) = protcoords(m,s)%x - protcoords(m,l)%x
          if (dx1(s) > gridsize/2) dx1(s) = dx1(s) - gridsize
          if (dx1(s) < -gridsize/2) dx1(s) = dx1(s) + gridsize
          dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if (dy1(s) > gridsize/2) dy1(s) = dy1(s) - gridsize
          if (dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize
          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if (dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if (dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize
       end do
       !write(6,*) probs
       if (probs == 1) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,1,0,-1,0,0,0,0,1,dx1(s),dy1(s),dz1(s))
          end do

          do s=l,p,1
             if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,3,1)
             end if
             !write(6,*) m,s,tbo(s)
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,3,1)
          end do



       else if (probs ==2) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,-1,0,1,0,0,0,0,1,dx1(s),dy1(s),dz1(s))
          end do
          do s = l,p,1
             if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,3,-1)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,3,-1)
          end do

       else if (probs ==3 ) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,-1,0,0,0,-1,0,0,0,1,dx1(s),dy1(s),dz1(s))
          end do
          do s = l,p,1
             if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,3,0)
             end if
          end do
          do s = l,p-1,1
             call ribond(m,s,tbdi,3,0)
          end do

       end if
       !write(6,*) 'z done'
    end if


if(l>1)then
    if((bdi(m,l-1) == -1*tbo(l))) then 
       crankcont = .false.
       goto 43
    end if
 end if

 if(p<maxl) then
    if((bdi(m,p) == tbo(p))) then
       crankcont = .false.
       goto 43
    end if
 end if
       
    tempcoord(l)%x = protcoords(m,l)%x
    tempcoord(p)%x = protcoords(m,p)%x
    tempcoord(l)%y = protcoords(m,l)%y
    tempcoord(p)%y = protcoords(m,p)%y
    tempcoord(l)%z = protcoords(m,l)%z
    tempcoord(p)%z = protcoords(m,p)%z

    if(p<maxl) tbdi(p) = bdi(m,p)


    !if (pivotx .eqv. .false. .and. pivoty .eqv. .false. .and. pivotz .eqv. .false.) then
    !crankcont = .false.
    !end if

    deltaenergy = 0.0d0    
    !pr = m

    do str = l+1,p-1,1
       do st = 1,maxl,1
          if((st < l+1) .or. (st > p-1)) then
             if (overlaps(m,str,st,tempcoord) .eqv. .false.) then
                crankcont = .false.
                goto 43
             end if
          end if
       end do
    end do

    do pr = 1,nprotein,1
       if(pr /= m) then
          do str = l+1,p-1,1
         maxback = chlen(pr)
             do st = 1,maxback,1
                if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                   crankcont = .false.
                   goto 43
                end if
             end do
          end do
       end if
    end do

    !write(6,*) 'a'
    call intrabondallocateupper(m,l,p,1,l-1,tempcoord,tempisobond,tba,deltaenergy,tbo)
    !write(6,*) 'b'
    call intrabondallocatelower(m,l,p,p+1,maxl,tempcoord,tempisobond,tba,deltaenergy,tbo)
    !write(6,*) 'c'
    call interbondallocate(m,l,p,tempcoord,tempintbond,tbr,deltaenergy,tbo)
    !write(6,*) 'd'
    if(Energydecision(deltaenergy) .eqv. .True.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l,p,tempcoord,tbo,tbdi,maxl)
       call updateintrabond(m,l,p,1,l-1,tempisobond,tba)
       call updateintrabond(m,l,p,p+1,maxl,tempisobond,tba)
       call updateinterbond(m,l,p,tempintbond,tbr)
    else 
       crankcont = .false.
       goto 43
    end if
!call comfind(m,.TRUE.)
    !       write(6,*) time, 'crank'
    successful = successful + 1
    call counts(1,0,0,0)
43  if (crankcont .eqv. .false.) then
       reject = reject + 1
       call counts(2,0,0,0)
    end if

  end subroutine crankshaftmove


  !call rbond(m,s,bonddd(m,s),tbo,3,1)
  subroutine rbond(m,s,bnd,tbo,axis,rotor)
    integer,intent(in):: m,s,bnd,axis,rotor
    integer,dimension(:),intent(inout)::tbo
    logical :: bm
    integer :: sign

    !write(6,*) 'fd'
    bm = .true.
    if(abs(bnd) == axis) then
       bm = .true.
       tbo(s) = bnd
       return
    end if
    !write(6,*) 'ab'

    if(rotor == 0) then
       bm = .true.
       tbo(s) = -1*bnd
       return
    end if
    !write(6,*) 'cd'
    !write(6,*) 'abc'
    !write(6,*) bonddd(m,s),m,s
    sign = bonddd(m,s)/abs(bonddd(m,s))
    !write(6,*) 'def'
    if(axis == 1) then
       if(rotor == -1) then
          if(abs(bonddd(m,s)) == 2) tbo(s) = sign*(abs(bonddd(m,s)) + 1)
          if(abs(bonddd(m,s)) == 3) tbo(s) = -1*sign*(abs(bonddd(m,s)) - 1)
          return
       else if(rotor == 1) then
          if(abs(bonddd(m,s)) == 2) tbo(s) = -1*sign*(abs(bonddd(m,s)) + 1)
          if(abs(bonddd(m,s)) == 3) tbo(s) = sign*(abs(bonddd(m,s)) - 1)
          return
       end if
    else if(axis ==2 ) then
       if(rotor == -1) then
          if(abs(bonddd(m,s)) == 1) tbo(s) = -1*sign*(abs(bonddd(m,s)) + 2)
          if(abs(bonddd(m,s)) == 3) tbo(s) = sign*(abs(bonddd(m,s)) - 2)        
          return
       else if(rotor == 1) then
          if(abs(bonddd(m,s)) == 1) tbo(s) = sign*(abs(bonddd(m,s)) + 2)
          if(abs(bonddd(m,s)) == 3) tbo(s) = -1*sign*(abs(bonddd(m,s)) - 2)
          return
       end if

    else if(axis ==3) then
       if(rotor == -1) then
          if(abs(bonddd(m,s)) == 1) tbo(s) = sign*(abs(bonddd(m,s)) + 1)
          if(abs(bonddd(m,s)) == 2) tbo(s) = -1*sign*(abs(bonddd(m,s)) - 1)        
          return
       else if(rotor == 1) then
          if(abs(bonddd(m,s)) == 1) tbo(s) = -1*sign*(abs(bonddd(m,s)) + 1)
          if(abs(bonddd(m,s)) == 2) tbo(s) = sign*(abs(bonddd(m,s)) - 1)         
       end if

    end if


    !87   if(bm .eqv. .true.) continue


  end subroutine rbond


  subroutine ribond(m,s,tbdi,axis,rotor)
    integer,intent(in):: m,s,axis,rotor
    integer,dimension(:),intent(inout)::tbdi
    logical :: bm
    integer :: sign

    !write(6,*) 'fd'
    bm = .true.
    if(abs(bdi(m,s)) == axis) then
       bm = .true.
       tbdi(s) = bdi(m,s)
       return
    end if
    !write(6,*) 'ab'

    if(rotor == 0) then
       bm = .true.
       tbdi(s) = -1*bdi(m,s)
       return
    end if
    !write(6,*) 'beadno',m,s, bdi(m,s)
    sign = bdi(m,s)/abs(bdi(m,s))

    if(axis == 1) then
       if(rotor == -1) then
          if(abs(bdi(m,s)) == 2) tbdi(s) = sign*(abs(bdi(m,s)) + 1)
          if(abs(bdi(m,s)) == 3) tbdi(s) = -1*sign*(abs(bdi(m,s)) - 1)
          return
       else if(rotor == 1) then
          if(abs(bdi(m,s)) == 2) tbdi(s) = -1*sign*(abs(bdi(m,s)) + 1)
          if(abs(bdi(m,s)) == 3) tbdi(s) = sign*(abs(bdi(m,s)) - 1)
          return
       end if
    else if(axis ==2 ) then
       if(rotor == -1) then
          if(abs(bdi(m,s)) == 1) tbdi(s) = -1*sign*(abs(bdi(m,s)) + 2)
          if(abs(bdi(m,s)) == 3) tbdi(s) = sign*(abs(bdi(m,s)) - 2)        
          return
       else if(rotor == 1) then
          if(abs(bdi(m,s)) == 1) tbdi(s) = sign*(abs(bdi(m,s)) + 2)
          if(abs(bdi(m,s)) == 3) tbdi(s) = -1*sign*(abs(bdi(m,s)) - 2)
          return
       end if

    else if(axis ==3) then
       if(rotor == -1) then
          if(abs(bdi(m,s)) == 1) tbdi(s) = sign*(abs(bdi(m,s)) + 1)
          if(abs(bdi(m,s)) == 2) tbdi(s) = -1*sign*(abs(bdi(m,s)) - 1)        
          return
       else if(rotor == 1) then
          if(abs(bdi(m,s)) == 1) tbdi(s) = -1*sign*(abs(bdi(m,s)) + 1)
          if(abs(bdi(m,s)) == 2) tbdi(s) = sign*(abs(bdi(m,s)) - 1)         
       end if

    end if
    !write(6,*) 'finish ribond'

    !87   if(bm .eqv. .true.) continue


  end subroutine ribond



  subroutine counts(cs,ra,pv,re)
    integer,intent(in):: cs,ra,pv,re
    if(cs /= 0) then
       if (cs == 1) cracc = cracc + 1
       if (cs == 2) crrej = crrej + 1
       !write(6,*) crrej,cracc
    else if(ra /= 0) then
       !write(6,*) 'counted'
       if (ra == 1) raacc = raacc + 1
       if (ra == 2) rarej = rarej + 1
    else if(pv /= 0) then
       if (pv == 1) pvacc = pvacc + 1
       if (pv == 2) pvrej = pvrej + 1
    else if(re /= 0) then
       if(re ==1 ) reacc = reacc + 1
       if(re == 2) rerej = rerej + 1
    end if


  end subroutine counts

  subroutine rightanglemove
    !performs a diagonal flip on a bead that has its two connecting beads perpendicular to each other
    integer :: m,l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz
    double precision :: deltaenergy
    logical :: rac,racc
    integer:: st,pr,maxl,maxback
    type(basicp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbo,tbdi
    integer,dimension(:,:), allocatable :: tempisobond,tba
    integer,dimension(:,:,:),allocatable :: tempintbond,tbr
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tba(maxlength,maxlength))
    allocate(tbr(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(tbo(maxlength))
    allocate(tbdi(maxlength))


    m =int(ran2(seed)*(nprotein))+1
        maxl = chlen(m)
    ! l cannot be equal to 1 or maxlength
    l = int(ran2(seed)*(maxl-2))+2
    !if(time> 42900) write(6,*) 'right', m

    rac = .true.

    if ((abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= 2.0) &
         .and. (abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= 2.0) &
         .and. (abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= 2.0) &
         .and. (abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= gridsize-2.0) &
         .and. (abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= gridsize-2.0) &
         .and. (abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= gridsize-2.0)) then
       rac = .true.
    else 
       rac = .false.
    end if

    if (rac .eqv. .true.) then    

       deltax1 = protcoords(m,l+1)%x - protcoords(m,l)%x
       if(deltax1 == gridsize-1) deltax1 = -1
       if(deltax1 == 1-gridsize) deltax1 = 1
       deltax2 =protcoords(m,l-1)%x - protcoords(m,l)%x
       if(deltax2 == gridsize-1) deltax2 = -1
       if(deltax2 == 1-gridsize) deltax2 = 1
       deltay1=protcoords(m,l+1)%y - protcoords(m,l)%y
       if(deltay1 == gridsize-1) deltay1 = -1
       if(deltay1 == 1-gridsize) deltay1 = 1
       deltay2 =protcoords(m,l-1)%y - protcoords(m,l)%y
       if(deltay2 == gridsize-1) deltay2 = -1
       if(deltay2 == 1-gridsize) deltay2 = 1
       deltaz1 = protcoords(m,l+1)%z - protcoords(m,l)%z
       if(deltaz1 == gridsize-1.0) deltaz1 = -1
       if(deltaz1 == 1-gridsize) deltaz1 = 1
       deltaz2= protcoords(m,l-1)%z - protcoords(m,l)%z
       if(deltaz2 == gridsize-1) deltaz2 = -1
       if(deltaz2 == 1-gridsize) deltaz2 = 1

       dx = 0
       dy = 0
       dz = 0

       dx = (deltax1 + deltax2)
       dy = (deltay1 + deltay2)
       dz = (deltaz1 + deltaz2)


       call rightvector(m,l,dx,dy,dz,tbo,racc,tbdi,maxl)

       if(racc .eqv. .false.) then
          rac = .false.
          goto 37
       end if
       !write(6,*) tbo(l-1),tbo(l),tbo(l+1)
       tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
       tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
       tempcoord(l)%z = modulo(protcoords(m,l)%z+dz-1,gridsize)+1

       tempcoord(l-1)%x = protcoords(m,l-1)%x
       tempcoord(l+1)%x = protcoords(m,l+1)%x
       tempcoord(l-1)%y = protcoords(m,l-1)%y
       tempcoord(l+1)%y = protcoords(m,l+1)%y
       tempcoord(l-1)%z = protcoords(m,l-1)%z
       tempcoord(l+1)%z = protcoords(m,l+1)%z
if(l<maxl-1) then
       tbdi(l+1) = bdi(m,l+1)
end if
       deltaenergy = 0.0d0        
       do st = 1,maxl,1
          if(st /=l) then
             if(overlaps(m,l,st,tempcoord) .eqv. .false.) then
                rac = .false.
                goto 37
             end if
          end if
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
        maxback = chlen(pr)
             do st = 1,maxback,1
                if(overlaps(pr,l,st,tempcoord) .eqv. .false.) then
                   rac = .false.
                   goto 37
                end if
             end do
          end if
       end do

       call intrabondallocateupper(m,l-1,l+1,1,l-3,tempcoord,tempisobond,tba,deltaenergy,tbo)
       call intrabondallocatelower(m,l-1,l+1,l+3,maxl,tempcoord,tempisobond,tba,deltaenergy,tbo)

       call interbondallocate(m,l-1,l+1,tempcoord,tempintbond,tbr,deltaenergy,tbo)


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l-1,l+1,tempcoord,tbo,tbdi,maxl)
          call updateintrabond(m,l-1,l+1,1,l-3,tempisobond,tba)
          call updateintrabond(m,l-1,l+1,l+3,maxl,tempisobond,tba)
          call updateinterbond(m,l-1,l+1,tempintbond,tbr)
          !if(m==5 .and. l==8)write(6,*) m,l,'right'
          !call comfind(m,.TRUE.)
          successful = successful + 1
          call counts(0,1,0,0)
       else
          rac = .false.
          goto 37
       end if
    end if

37  if (rac .eqv. .false.) then
       reject = reject + 1

       call counts(0,2,0,0)
    end if

  end subroutine rightanglemove

  subroutine rightvector(m,l,dx,dy,dz,tbo,racc,tbdi,maxl)
    integer,intent(in):: m,l,dx,dy,dz,maxl
    integer,dimension(:),intent(inout)::tbo,tbdi
    logical,intent(inout)::racc
    integer :: runn,k
    if(dx == 0) runn =1
    if(dy == 0) runn = 2
    if(dz == 0) runn =3
    racc = .true.
    if(runn ==1) then
       if(dy + dz ==0) then
          do k=l-1,l+1,1
             if(abs(bonddd(m,k)) == 1) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 2) tbo(k) = 3*(bonddd(m,k)/2)
             if(abs(bonddd(m,k)) ==3) tbo(k) = 2*(bonddd(m,k)/3)
          end do
       else if(dy+dz /= 0) then
          do k=l-1,l+1,1
             if(abs(bonddd(m,k)) == 1) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 2) tbo(k) = -3*(bonddd(m,k)/2)
             if(abs(bonddd(m,k)) ==3) tbo(k) = -2*(bonddd(m,k)/3)
          end do
       end if

       goto 52
    else if(runn == 2) then
       if(dx + dz == 0) then
          do k = l-1,l+1,1
             if(abs(bonddd(m,k)) == 2) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 1) tbo(k) = 3*(bonddd(m,k))
             if(abs(bonddd(m,k)) ==3) tbo(k) = 1*(bonddd(m,k)/3)
          end do
       else if(dx + dz /= 0) then
          do k = l-1,l+1,1
             if(abs(bonddd(m,k)) == 2) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 1) tbo(k) = -3*(bonddd(m,k))
             if(abs(bonddd(m,k)) ==3) tbo(k) = -(bonddd(m,k)/3)
          end do
       end if
       goto 52
    else if(runn ==3) then
       if(dx+dy == 0) then
          do k = l-1,l+1,1
             if(abs(bonddd(m,k)) == 3) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 2) tbo(k) = (bonddd(m,k)/2)
             if(abs(bonddd(m,k)) ==1) tbo(k) = 2*(bonddd(m,k))
          end do
       else if(dx+dy /= 0) then
          do k = l-1,l+1,1
             if(abs(bonddd(m,k)) == 3) tbo(k) = -1*bonddd(m,k)
             if(abs(bonddd(m,k)) == 2) tbo(k) = -(bonddd(m,k)/2)
             if(abs(bonddd(m,k)) ==1) tbo(k) = -2*(bonddd(m,k))
          end do
       end if
       goto 52
    end if



52  if(racc .eqv. .true.) then
       !write(6,*) 'l',l,maxlength,protcoords(m,l)%species
       if(l<maxl-1) then
          if(bdi(m,l+1)== tbo(l+1)) racc = .false.
       end if
       if(l>2) then
          if(-(bdi(m,l-2))== tbo(l-1)) racc = .false.
       end if
    end if

    if(racc .eqv. .false.) return

    if(runn ==1) then
       if(dy + dz ==0) then
          do k=l-1,l,1
             if(abs(bdi(m,k)) == 2) tbdi(k) = 3*(bdi(m,k)/2)
             if(abs(bdi(m,k)) ==3) tbdi(k) = 2*(bdi(m,k)/3)
          end do
       else if(dy+dz /= 0) then
          do k=l-1,l,1
             if(abs(bdi(m,k)) == 2) tbdi(k) = -3*(bdi(m,k)/2)
             if(abs(bdi(m,k)) ==3) tbdi(k) = -2*(bdi(m,k)/3)
          end do
       end if

       goto 51
    else if(runn == 2) then
       if(dx + dz == 0) then
          do k = l-1,l,1
             if(abs(bdi(m,k)) == 1) tbdi(k) = 3*(bdi(m,k))
             if(abs(bdi(m,k)) ==3) tbdi(k) = 1*(bdi(m,k)/3)
          end do
       else if(dx + dz /= 0) then
          do k = l-1,l,1
             if(abs(bdi(m,k)) == 1) tbdi(k) = -3*(bdi(m,k))
             if(abs(bdi(m,k)) ==3) tbdi(k) = -(bdi(m,k)/3)
          end do
       end if
       goto 51
    else if(runn ==3) then
       if(dx+dy == 0) then
          do k = l-1,l,1
             if(abs(bdi(m,k)) == 2) tbdi(k) = (bdi(m,k)/2)
             if(abs(bdi(m,k)) ==1) tbdi(k) = 2*(bdi(m,k))
          end do
       else if(dx+dy /= 0) then
          do k = l-1,l,1
             if(abs(bdi(m,k)) == 2) tbdi(k) = -(bdi(m,k)/2)
             if(abs(bdi(m,k)) ==1) tbdi(k) = -2*(bdi(m,k))
          end do
       end if
       goto 51
    end if

51  if(racc .eqv. .true.) continue
  end subroutine rightvector

  subroutine pivot(m,l,maxl)
    !rotates a section of chain around the selected bead
    integer,intent(in)::m,l,maxl
    integer :: b,g,xhold,yhold,zhold,str,st,pr,choose3,maxback
    logical :: pivcont
    double precision :: deltaenergy
    integer,dimension(:),allocatable :: delx,dely,delz
    type(basicp),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond,tba
    integer,dimension(:,:,:),allocatable :: tempintbond,tbr
    integer,dimension(:),allocatable::tbo,tbdi
    allocate(tempisobond(maxl,maxl))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tba(maxl,maxl))
    allocate(tbr(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxl))
    allocate(delx(maxl))
    allocate(dely(maxl))   
    allocate(delz(maxl))
    allocate(tbo(maxl))
    allocate(tbdi(maxl-1))
    choose3 = int(ran2(seed)*5) + 1
    pivcont = .true. 

    t = time + 1

    if(l > maxl/2) then
       do g = l+1,maxl,1
          delx(g) = protcoords(m,g)%x - protcoords(m,l)%x
          xhold = delx(g)
          if(delx(g) > gridsize/2) delx(g) =  xhold-gridsize
          if(delx(g) < -gridsize/2) delx(g) = gridsize + xhold
          dely(g) = protcoords(m,g)%y - protcoords(m,l)%y
          yhold = dely(g)
          if(dely(g) > gridsize/2) dely(g) = yhold - gridsize
          if(dely(g) < -gridsize/2) dely(g) = gridsize + yhold
          delz(g) = protcoords(m,g)%z - protcoords(m,l)%z
          zhold = delz(g)
          if(delz(g) > gridsize/2) delz(g) =  zhold-gridsize
          if(delz(g) < -gridsize/2) delz(g) = gridsize + zhold

       end do

       !write(6,*) choose3
       if (choose3 ==1) then

          do b = l+1,maxl,1
             call rotate(m,l,b,tempcoord,0,-1,0,1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,-1)
             end if
          end do
          do b = l,maxl-1,1
             call ribond(m,b,tbdi,3,-1)
          end do

       else if (choose3 == 2) then
          do b =l+1,maxl,1
             call rotate(m,l,b,tempcoord,0,1,0,-1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,1)
             end if
          end do
          do b = l,maxl-1,1
             call ribond(m,b,tbdi,3,1)
          end do

       else if (choose3 ==3) then

          do b =l+1,maxl,1
             call rotate(m,l,b,tempcoord,-1,0,0,0,-1,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,0)
             end if
          end do
          do b = l,maxl-1,1

             call ribond(m,b,tbdi,3,0)

          end do
       else if (choose3 ==4) then

          do b =l+1,maxl,1
             call rotate(m,l,b,tempcoord,0,0,-1,0,1,0,1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,2,1)
             end if
          end do
          do b = l,maxl-1,1
             call ribond(m,b,tbdi,2,1)
          end do
       else if (choose3 == 5) then

          do b = l+1,maxl,1
             call rotate(m,l,b,tempcoord,0,0,1,0,1,0,-1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,2,-1)
             end if
          end do
          do b = l,maxl-1,1
             call ribond(m,b,tbdi,2,-1)
          end do

       end if
if(l<maxl) then
       if(tbdi(l) == bonddd(m,l))  then
          pivcont = .false.
          goto 75
       end if
end if
       tempcoord(l)%x = protcoords(m,l)%x
       tempcoord(l)%y = protcoords(m,l)%y
       tempcoord(l)%z = protcoords(m,l)%z
       tbo(l) = bonddd(m,l)



       do str = l+1,maxl,1
          do st = 1,l,1
             if(overlaps(m,str,st,tempcoord) .eqv. .false.) then
                pivcont = .false.
                goto 75
             end if
          end do
       end do

       deltaenergy = 0.0d0
       do pr = 1,nprotein,1
          if(pr /= m) then
             maxback = chlen(pr)
             do str = l+1,maxl,1
                do st = 1,maxback,1
                   if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end if
       end do


       call intrabondallocateupper(m,l+1,maxl,1,l,tempcoord,tempisobond,tba,deltaenergy,tbo)

       call interbondallocate(m,l+1,maxl,tempcoord,tempintbond,tbr,deltaenergy,tbo)

!if(NINT(deltaenergy) /= 0) write(6,*) 'pivot',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then

          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l,maxl,tempcoord,tbo,tbdi,maxl)
          call updateintrabond(m,l+1,maxl,1,l,tempisobond,tba)
          call updateinterbond(m,l+1,maxl,tempintbond,tbr)
          !call comfind(m,.TRUE.)
          successful = successful + 1
          call counts(0,0,1,0)

       else
          pivcont = .false.
          goto 75
       end if

    else if (l <= maxl/2) then

       do g = 1,l-1,1
          delx(g) = protcoords(m,g)%x - protcoords(m,l)%x
          xhold = delx(g)
          if(delx(g) > gridsize/2) delx(g) =  xhold-gridsize
          if(delx(g) < -1.0*gridsize/2) delx(g) = gridsize + xhold
          dely(g) = protcoords(m,g)%y - protcoords(m,l)%y
          yhold = dely(g)
          if(dely(g) > gridsize/2) dely(g) = yhold - gridsize
          if(dely(g) < -1.0*gridsize/2) dely(g) = gridsize + yhold
          delz(g) = protcoords(m,g)%z - protcoords(m,l)%z
          zhold = delz(g)
          if(delz(g) > gridsize/2) delz(g) =  zhold-gridsize
          if(delz(g) < -1.0*gridsize/2) delz(g) = gridsize + zhold
       end do

       if (choose3 ==1) then

          do b = 1,l-1,1
             call rotate(m,l,b,tempcoord,0,-1,0,1,0,0,0,0,1,delx(b),dely(b),delz(b))

             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,-1)
             end if
          end do
          do b = 1,l-1,1
             call ribond(m,b,tbdi,3,-1)
          end do

       else if (choose3 == 2) then 
          do b =1,l-1,1
             call rotate(m,l,b,tempcoord,0,1,0,-1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,1)
             end if
          end do
          do b = 1,l-1,1
             call ribond(m,b,tbdi,3,1)
          end do

       else if (choose3 == 3) then
          do b =1,l-1,1
             call rotate(m,l,b,tempcoord,-1,0,0,0,-1,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,0)
             end if
          end do
          do b = 1,l-1,1
             call ribond(m,b,tbdi,3,0)
          end do

       else if (choose3==4) then
          do b =1,l-1,1
             call rotate(m,l,b,tempcoord,0,0,-1,0,1,0,1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,2,1)
             end if
          end do
          do b = 1,l-1,1
             call ribond(m,b,tbdi,2,1)
          end do

       else if (choose3==5) then
          do b = 1,l-1,1
             call rotate(m,l,b,tempcoord,0,0,1,0,1,0,-1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,2,-1)
             end if
          end do
          do b = 1,l-1,1
             call ribond(m,b,tbdi,2,-1)
          end do
       end if
if(l>1) then
       if(-(tbdi(l-1)) == bonddd(m,l)) then
          pivcont = .false.
          goto 75
       end if
end if


       deltaenergy = 0.0d0
       do str = 1,l-1,1
          do st = l,maxl,1
             if(overlaps(m,str,st,tempcoord) .eqv. .false.) then
                pivcont = .false.
                goto 75
             end if
          end do
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
             do str = 1,l-1,1
                maxback = chlen(pr)
                do st = 1,maxback,1
                   if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end if
       end do

       call intrabondallocatelower(m,1,l-1,l,maxl,tempcoord,tempisobond,tba,deltaenergy,tbo)
       call interbondallocate(m,1,l-1,tempcoord,tempintbond,tbr,deltaenergy,tbo)
!if(NINT(deltaenergy) /= 0) write(6,*) 'pivot',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then

          totalenergy = totalenergy + deltaenergy
          call updatepos(m,1,l-1,tempcoord,tbo,tbdi,maxl)
          call updateintrabond(m,1,l-1,l,maxl,tempisobond,tba)
          call updateinterbond(m,1,l-1,tempintbond,tbr)
          !call comfind(m,.TRUE.)
          successful = successful + 1
          call counts(0,0,1,0)
       else
          pivcont = .false.
          goto 75
       end if
    end if

75  if (pivcont .eqv. .false.) then
       reject = reject + 1
       call counts(0,0,2,0)
    end if

  end subroutine pivot


  subroutine intrabondallocatelower(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,tba,deltaenergy,tbo)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    integer,dimension(:),intent(in)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond,tba
    integer::bead1,bead2,bdir
    Type(basicp),dimension(:),intent(inout) :: tempcoord
    logical :: adjver
    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead2 > bead1 + 2) then
             if(isobond(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(chain1,bead1)%type)
             if(ba(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraen
             call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
             if ((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + intraen
                tba(bead1,bead2) = 1
                if((bdir == -1*tbo(bead1)) .and. (bdir == bonddd(chain1,bead2))) then
                   deltaenergy = deltaenergy + intraenergy(protcoords(chain1,bead1)%type)
                   tempisobond(bead1,bead2) = 1
                else 
                   tempisobond(bead1,bead2) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tba(bead1,bead2) = 0
                tempisobond(bead1,bead2) = 0
             end if

          end if
       end do
    end do
    !write(6,*) 'intra lower energy', deltaenergy 

  end subroutine intrabondallocatelower

  subroutine intrabondallocateupper(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,tba,deltaenergy,tbo)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond,tba
    integer,dimension(:),intent(in)::tbo
    integer::bead1,bead2,bdir
    Type(basicp),dimension(:),intent(inout) :: tempcoord
    logical :: adjver
    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead1 > bead2 + 2) then
             if(isobond(chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(chain1,bead1)%type)
             if(ba(chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - intraen
             call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + intraen
                tba(bead2,bead1) = 1
                if((bdir == -1*tbo(bead1)) .and. bdir== bonddd(chain1,bead2)) then
                   !write(6,*) bdir,tbo(bead1),bonddd(chain1,bead2),adjver
                   deltaenergy = deltaenergy + intraenergy(protcoords(chain1,bead1)%type)
                   tempisobond(bead2,bead1) = 1
                else 
                   tempisobond(bead2,bead1) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tba(bead2,bead1) = 0
                tempisobond(bead2,bead1) = 0
             end if
          end if
       end do
    end do
    !write(6,*) 'intra upper energy', deltaenergy 

  end subroutine intrabondallocateupper

  subroutine interbondallocate(chain1,beadmin,beadmax,tempc,tempintbond,tbr,deltaenergy,tbo)

    integer,intent(in)::chain1,beadmin,beadmax
    integer,dimension(:),intent(in)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:,:),intent(inout) :: tempintbond,tbr
    integer::bead1,bead2,chain2,bdir,maxback,f
    Type(basicp),dimension(:),intent(in) :: tempc
    logical :: adjver

    if(nprotein == 1) return
    !write(6,*) 'abcde',beadmin,beadmax
    do chain2 = 1,chain1-1,1

             maxback = chlen(chain2)
       do bead1 = beadmin,beadmax,1
          do bead2 = 1,maxback,1
             !tempintbond(chain2,bead1,bead2) = 0
             if(intbond(chain2,chain1,bead2,bead1) == 1)then
                deltaenergy = deltaenergy - &
                     interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                !write(6,*) 'bondbreak',intbond(chain2,chain1,bead2,bead1),chain1,chain2,bead1,bead2
             end if
             if(br(chain2,chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - &
                  interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             call adjacent(chain2,bead1,bead2,tempc,adjver,bdir)
             if((adjver .eqv. .true.)) then
                deltaenergy = deltaenergy + &
                     interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                tbr(chain2,bead1,bead2) = 1
                if((bdir == -1*tbo(bead1)) .and. (bdir ==bonddd(chain2,bead2))) then
                   deltaenergy = deltaenergy + &
                        interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                   tempintbond(chain2,bead1,bead2) = 1
                   !write(6,*) 'bondform',intbond(chain2,chain1,bead2,bead1),chain1,chain2,bead1,bead2
                else
                   tempintbond(chain2,bead1,bead2) = 0
                end if
             else if(adjver .eqv. .false.) then
                tbr(chain2,bead1,bead2) = 0
                tempintbond(chain2,bead1,bead2) = 0
             end if

          end do
       end do
    end do


    do chain2 = chain1+1,nprotein,1

        maxback = chlen(chain2)
       do bead1 = beadmin,beadmax,1
          do bead2 = 1,maxback,1
             !tempintbond(chain2,bead1,bead2) = 0
             if(intbond(chain1,chain2,bead1,bead2) == 1) then
                deltaenergy = deltaenergy - &
                     interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                !write(6,*) 'bondbreak',intbond(chain1,chain2,bead1,bead2),chain1,chain2,bead1,bead2
             end if
             if(br(chain1,chain2,bead1,bead2) == 1) deltaenergy = deltaenergy - &
                  interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
             call adjacent(chain2,bead1,bead2,tempc,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + &
                     interen(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                tbr(chain2,bead1,bead2) = 1
                if((bdir == -1*tbo(bead1)) .and. (bdir == bonddd(chain2,bead2))) then
                   deltaenergy = deltaenergy + &
                        interenergy(protcoords(chain1,bead1)%type,protcoords(chain2,bead2)%type)
                   tempintbond(chain2,bead1,bead2) = 1
                   !write(6,*) 'bondform',intbond(chain1,chain2,bead1,bead2),chain1,chain2,bead1,bead2
                else
                   tempintbond(chain2,bead1,bead2) = 0
                end if
             else if(adjver .eqv. .false.) then
                tbr(chain2,bead1,bead2) =0
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do
    !write(6,*) 'a',beadmin,beadmax
    !do chain2 = 1,nprotein,1
    !if(chain2 /= chain1) then
    !if(protcoords(chain2,1)%species == 1) maxback = maxlength1
    !if(protcoords(chain2,1)%species == 2) maxback = maxlength2
    !do bead1 = beadmin,beadmax,1
    !  do bead2 = 1,maxback
    !    write(6,*) 'bond set',tempintbond(chain2,bead1,bead2),chain2,bead1,bead2
    ! end do
    !end do
    !end if
    !end do

    !write(6,*) 'inter energy', deltaenergy 
  end subroutine interbondallocate

  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,tbo,tbdi,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    integer,dimension(:),intent(inout)::tbo,tbdi
    Type(basicp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding

    
    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
       bonddd(chainnum,beadnum) = tbo(beadnum)
       !if(tbdi(beadnum) == 0) write(6,*) 'fail',chainnum,beadnum
       if(beadnum<maxl) bdi(chainnum,beadnum) = tbdi(beadnum)
    end do

    call comfind(chainnum,.TRUE.)
!call rms(chainnum)
  end subroutine updatepos

  subroutine updatecluspos(tempcluscoord,ctbo,ctbdi,tcom,cllist,clsize)
    !updates bead positions
    integer,dimension(:),intent(in) :: cllist
    integer,dimension(:,:),intent(in)::ctbo,ctbdi
    integer,intent(in) :: clsize
    Type(basicp),dimension(:,:),intent(in) :: tempcluscoord
    type(centremass),dimension(:),intent(in)::tcom
    integer ::b,l,a,maxback,ra
    !moves beads to new positions and reassigns isobonding

    do ra = 1,clsize
a = cllist(ra)
        maxback = chlen(a)
          call rms(a)
          do l = 1,maxback,1
             protcoords(a,l)%x = tempcluscoord(a,l)%x
             protcoords(a,l)%y = tempcluscoord(a,l)%y
             protcoords(a,l)%z = tempcluscoord(a,l)%z
             bonddd(a,l) = ctbo(a,l)
             if(l<maxback) then
                bdi(a,l) = ctbdi(a,l)
             if(ctbdi(a,l) == 0) write(6,*) 'update fail',a,l
          end if
          end do
          com(2,a)%x = com(1,a)%x
          com(2,a)%y = com(1,a)%y
          com(2,a)%z = com(1,a)%z
          com(1,a)%x = tcom(a)%x
          com(1,a)%y = tcom(a)%y
          com(1,a)%z = tcom(a)%z
          !call comfind(a,.TRUE.)
          call rms(a)
          call tracking(a)
       end do

    !write(6,*) 'successful cluster move'
  end subroutine updatecluspos

  subroutine bondcheck

    integer:: m,l,f
    !open(12, file = 'bondintra.dat', action = 'write')
    do m = 1,nprotein
       do l = 1,maxlength,1
          do f = 1,maxlength,1
             if ((f < l -2) .or. (f > l+ 2)) write(12,*) l,f,isobond(m,l,f)
          end do
       end do
    end do

  end subroutine bondcheck

  subroutine updateintrabond(chainnum,beadmin,beadmax,statmin,statmax,tempisobond,tba)
    integer,intent(in):: chainnum,beadmin,beadmax,statmin,statmax
    integer ::g,beadnum
    integer,dimension(:,:),intent(in) :: tempisobond,tba
    do beadnum = beadmin,beadmax,1
       do g = statmin,statmax,1
          if(g<beadnum-2) then
             isobond(chainnum,g,beadnum) = tempisobond(g,beadnum)
             ba(chainnum,g,beadnum) = tba(g,beadnum) 
          else if(g> beadnum+2) then
             isobond(chainnum,beadnum,g) = tempisobond(beadnum,g)
             ba(chainnum,beadnum,g) = tba(beadnum,g)
          end if
       end do
    end do

  end subroutine updateintrabond


  subroutine clusterassign(chainnum,clno,clsize,cb)
    integer,intent(in) :: chainnum
    integer,dimension(:),intent(inout):: clno
        integer,dimension(:,:),intent(inout)::cb
    integer,intent(inout):: clsize
    integer:: m,l,g,f,clcount,bdir,maxlengthss,maxback,zzz,dum2,dum1,oldcl
    logical :: clusyes,adjver
    type(basicp),dimension(:),allocatable :: tempcoord

    allocate(tempcoord(maxlength))
    clusyes = .true.
    do m = 1,nprotein
       clno(m) = m
       do g = 1,nprotein
          cb(m,g) = 0
          end do
    end do

    do m = 1,nprotein-1,1
         maxback = chlen(m)
       do g = m+1,nprotein,1

           maxlengthss = chlen(g)

          do l = 1,maxback
             tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y
             tempcoord(l)%z = protcoords(m,l)%z
             do f = 1,maxlengthss
                call adjacent(g,l,f,tempcoord,adjver,bdir)              
                if(adjver.eqv. .true.) then
                   if(isbond .eqv. .true.) then
                      if(interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0) then
                      cb(m,g) = 1
                      cb(g,m) = 1
                      if(clno(g) < clno(m)) then
                         oldcl = clno(m)
                         do zzz = 1,nprotein
                            if(clno(zzz) == oldcl) clno(zzz) = clno(g)
                         end do
                      else if(clno(g)> clno(m)) then
                         oldcl =clno(g)
                         do zzz = 1,nprotein
                            if(clno(zzz) == oldcl) clno(zzz) = clno(m)
                         end do
                      end if
                      clusyes = .true.
end if
                   end if
                   if((bdir == bonddd(g,f)) .and. (bdir == (-1*bonddd(m,l))) .and. &
                        (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  /= 0.0)) then
                      if((m == 2)) write(49,*) time,g
                      if((g == 2)) write(49,*) time,m
                      cb(m,g) = 1
                      cb(g,m) = 1
                      if(clno(g) < clno(m)) then
                         oldcl = clno(m)
                         do zzz = 1,nprotein
                            if(clno(zzz) == oldcl) clno(zzz) = clno(g)
                         end do
                      else if(clno(g)> clno(m)) then
                         oldcl =clno(g)
                         do zzz = 1,nprotein
                            if(clno(zzz) == oldcl) clno(zzz) = clno(m)
                         end do
                      end if
                      clusyes = .true.                 
                   end if
                end if
             end do
          end do
       end do
    end do

    clsize = 0

    do m = 1,nprotein
       if(clno(m) == clno(chainnum)) clsize = clsize + 1
    end do
  end subroutine clusterassign


  
  subroutine updateinterbond(chainnum,beadmin,beadmax,tempintbond,tbr)
    integer,intent(in):: chainnum,beadmin,beadmax
    integer,dimension(:,:,:),intent(in) :: tempintbond,tbr
    integer ::g,beadnum,chain2,maxback
    !write(6,*) 'start'
    !do chain2 = 1,nprotein,1
    !if(chain2 /= chainnum) then
    !if(protcoords(chain2,1)%species == 1) maxback = maxlength1
    !if(protcoords(chain2,1)%species == 2) maxback = maxlength2
    !do beadnum = beadmin,beadmax,1
    !do g = 1,maxback
    !write(6,*) 'bond detail',tempintbond(chain2,beadnum,g),chain2,beadnum,g
    !end do
    !end do
    !end if
    !end do


    do chain2 = 1,chainnum-1,1
        maxback = chlen(chain2)
       do beadnum = beadmin,beadmax,1
          do g = 1,maxback,1
             intbond(chain2,chainnum,g,beadnum) = tempintbond(chain2,beadnum,g)
             br(chain2,chainnum,g,beadnum) = tbr(chain2,beadnum,g)
             !if(intbond(chain2,chainnum,g,beadnum) /= 0) write(6,*)'bondf1',intbond(chain2,chainnum,g,beadnum),&
             !tempintbond(chain2,beadnum,g),chain2,beadnum,g
          end do
       end do
    end do
    do chain2 = chainnum+1,nprotein,1
        maxback = chlen(chain2)
       do beadnum = beadmin,beadmax,1
          do g = 1,maxback,1
             intbond(chainnum,chain2,beadnum,g) = tempintbond(chain2,beadnum,g)
             br(chainnum,chain2,beadnum,g) = tbr(chain2,beadnum,g)
             !if(intbond(chain2,chainnum,g,beadnum) /= 0) write(6,*)'bondf2',intbond(chainnum,chain2,beadnum,g),&
             !tempintbond(chain2,beadnum,g),chain2,beadnum,g
          end do
       end do
    end do

  end subroutine updateinterbond

  subroutine reptation
    !perform reptation
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir,maxl,maxback,str
    logical :: reptcont,adjver
    type(basicp),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond,tba
    integer,dimension(:),allocatable::tbo,tbdi
    integer,dimension(:,:,:),allocatable :: tempintbond,tbr
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tba(maxlength,maxlength))
    allocate(tbr(nprotein,maxlength,maxlength))


    t = time + 1

!choose = -0.5
    choose =  ran2(seed) - 0.5
    !write(6,*) 'reptation start',choose
    
    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1
        maxl = chlen(m)

    allocate(tempcoord(maxl))
    allocate(tbo(maxl))
    allocate(tbdi(maxl-1))

    g1 = 1
    g2 = 1
    g3 = 1
    g4 = 1
    g5 = 1
    g6 = 1
    dx = 0
    dy = 0
    dz = 0
    direc = int(ran2(seed)*5) +1
        !if(time == 3387) write(6,*) 'choooooose =', choose,direc,m
    if (choose > 0.0) then

       !if(bonddd(m,maxl) == -(bdi(m,maxl-1))) then
       !reptcont = .false.
       !goto 83
       !end if
       dxx = protcoords(m,2)%x - protcoords(m,1)%x
       if (dxx > 1) dxx = -1
       if(dxx <-1) dxx = 1
       dyx = protcoords(m,2)%y - protcoords(m,1)%y
       if (dyx > 1) dyx = -1
       if (dyx < -1) dyx = 1
       dzx = protcoords(m,2)%z - protcoords(m,1)%z
       if (dzx > 1) dzx = -1
       if (dzx < -1) dzx = 1

       if (dxx == 1) g1 = 0
       if (dxx == -1) g2 = 0
       if (dyx == 1) g3 = 0
       if (dyx == -1) g4 = 0
       if (dzx == 1) g5 =0
       if (dzx == -1) g6 = 0


       if((direc == g1) .and. (g1 == 1)) then
          dx = 1
          tbdi(1) = -1
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == 1) reptcont = .false.
       else if((direc == (g1 +g2))  .and.  (g2 == 1)) then
          dx = -1
          tbdi(1) = 1
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == -1) reptcont = .false.
       else if((direc == (g1+g2+g3)) .and.  (g3 == 1)) then
          dy = 1
          tbdi(1) = -2
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == 2) reptcont = .false.
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          tbdi(1) = 2
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == -2) reptcont = .false.
       else if((direc == (g1+g2+g3+g4+g5)) .and. (g5 == 1)) then
          dz = 1
          tbdi(1) = -3
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == 3) reptcont = .false.
       else if((direc == (g1+g2+g3+g4+g5+g6)) .and.  (g6 == 1)) then
          dz = -1
          tbdi(1) = 3
          call bonddirection(5,tbdi(1),8,tbo(1))
          if(bonddd(m,1) == -3) reptcont = .false.
       end if

       if(reptcont .eqv. .false.) goto 83



       tempcoord(1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(1)%z = modulo(protcoords(m,1)%z + dz-1,gridsize)+1
       ! write(6,*) tbdi(1),dx,dy,dz
       ! write(6,*) tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z
       ! write(6,*) protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z

       deltaenergy = 0.0d0

       do st = 1,maxl-1,1
          if(overlaps(m,1,st,tempcoord) .eqv. .false.) then
             reptcont = .false.
             goto 83
          end if
       end do

       do pr =1,nprotein    
        maxback = chlen(pr)
          do st = 1,maxback,1
             if((pr /=m))then
                if(overlaps(pr,1,st,tempcoord) .eqv. .false.) then
                   reptcont = .false.
                   goto 83
                end if
             end if
          end do
       end do

       do st = 1,maxl-3,1
          if(isobond(m,st,maxl) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(m,st)%type)
          if(ba(m,st,maxl) == 1) deltaenergy = deltaenergy - intraen
       end do



       do st = 3,maxl-1
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             deltaenergy = deltaenergy + intraen
             tba(1,st+1) = 1
             if((bdir== -1*tbo(1)) .and. (bdir== bonddd(m,st))) then
                deltaenergy = deltaenergy + intraenergy(protcoords(m,st)%type)
                tempisobond(1,st+1) = 1
             else
                tempisobond(1,st+1) = 0
             end if
          else if((adjver.eqv. .false.)) then
             tba(1,st+1) = 0
             tempisobond(1,st+1) = 0
          end if
       end do

       do pr = 1,m-1,1
        maxback = chlen(pr)
          do st = 1,maxback
             if(intbond(pr,m,st,maxl) == 1) deltaenergy = deltaenergy - &
                  interenergy(protcoords(m,maxl)%type,protcoords(pr,st)%type)
             if(br(pr,m,st,maxl) == 1) deltaenergy = deltaenergy - interen(protcoords(m,maxl)%type,protcoords(pr,st)%type)
             call adjacent(pr,1,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + interen(protcoords(m,1)%type,protcoords(pr,st)%type)
                tbr(pr,1,st) = 1
                if((bdir== -1*tbo(1)) .and. (bdir== bonddd(pr,st))) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,1)%type,protcoords(pr,st)%type)
                   tempintbond(pr,1,st) = 1
                else
                   tempintbond(pr,1,st) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tbr(pr,1,st) = 0
                tempintbond(pr,1,st) = 0
             end if
          end do
       end do

       do pr = m+1,nprotein,1
        maxback = chlen(pr)
          do st = 1,maxback
             if(intbond(m,pr,maxl,st) == 1) deltaenergy = deltaenergy - &
                  interenergy(protcoords(m,maxl)%type,protcoords(pr,st)%type)
             if(br(m,pr,maxl,st) == 1) deltaenergy = deltaenergy -&
                  interen(protcoords(m,maxl)%type,protcoords(pr,st)%type)
             call adjacent(pr,1,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + interen(protcoords(m,1)%type,protcoords(pr,st)%type)
                tbr(pr,1,st) = 1
                if((bdir== -1*tbo(1)) .and. (bdir== bonddd(pr,st))) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,1)%type,protcoords(pr,st)%type)
                   tempintbond(pr,1,st) = 1
                else
                   tempintbond(pr,1,st) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tbr(pr,1,st) = 0
                tempintbond(pr,1,st) = 0
             end if
          end do
       end do



       !section to adjust for change in types


       do str = 1,maxl-1,1
          do pr = 1,m-1,1
             maxback = chlen(pr)
             do st = 1,maxback
                if (intbond(pr,m,st,str) ==1) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,str+1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interenergy(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if

                if (br(pr,m,st,str) ==1) then
                   deltaenergy = deltaenergy + interen(protcoords(m,str+1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interen(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if                
             end do
          end do
       end do

       do str = 1,maxl-1,1
          do pr = m+1,nprotein,1
             maxback = chlen(pr)
             do st = 1,maxback
                if (intbond(m,pr,str,st) ==1) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,str+1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interenergy(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if

                if (br(m,pr,str,st) ==1) then
                   deltaenergy = deltaenergy + interen(protcoords(m,str+1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interen(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if                
             end do
          end do
       end do


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          do l =maxl,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z
             bonddd(m,l) = bonddd(m,l-1)
             if(l<maxl) bdi(m,l) = bdi(m,l-1)

             do g = l-3,2,-1
                isobond(m,g,l) = isobond(m,g-1,l-1)
                ba(m,g,l) = ba(m,g-1,l-1)
             end do

             do chain2 = 1,m-1,1
           maxback = chlen(chain2)             
                do beads2 = 1,maxback
                   intbond(chain2,m,beads2,l) = intbond(chain2,m,beads2,l-1)
                   br(chain2,m,beads2,l) = br(chain2,m,beads2,l-1)
                end do
             end do

             do chain2 = m+1,nprotein,1
                maxback = chlen(chain2) 
                do beads2 = 1,maxback,1
                   intbond(m,chain2,l,beads2) = intbond(m,chain2,l-1,beads2)
                   br(m,chain2,l,beads2) = br(m,chain2,l-1,beads2)
                end do
             end do
          end do

          call updatepos(m,1,1,tempcoord,tbo,tbdi,maxl)
          call updateintrabond(m,1,1,4,maxl,tempisobond,tba)
          call updateinterbond(m,1,1,tempintbond,tbr)
          
          successful = successful + 1
          reptforward = reptforward + 1
          call counts(0,0,0,1)
       else
          reptcont = .false.
       end if
    else if (choose < 0.0) then

       !if(bonddd(m,1) == bdi(m,2)) then
       !reptcont = .false.
       !goto 83
       !end if


       dxx = protcoords(m,maxl-1)%x - protcoords(m,maxl)%x
       if (dxx > 1) dxx = -1
       if(dxx <-1) dxx = 1
       dyx = protcoords(m,maxl-1)%y - protcoords(m,maxl)%y
       if (dyx > 1) dyx = -1
       if (dyx < -1) dyx = 1
       dzx = protcoords(m,maxl-1)%z - protcoords(m,maxl)%z
       if (dzx > 1) dzx = -1
       if (dzx < -1) dzx = 1

       if (dxx == -1) g1 = 0
       if (dxx == 1) g2 = 0
       if (dyx == -1) g3 = 0
       if (dyx == 1) g4 = 0
       if (dzx == -1) g5 =0
       if (dzx == 1) g6 = 0


       if((direc == g1) .and. (g1 == 1)) then
          dx = 1
          tbdi(maxl-1) = 1
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == 1) reptcont = .false.
       else if((direc ==(g1 +g2)) .and.  (g2 == 1)) then
          dx = -1
          tbdi(maxl-1) = -1
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == -1) reptcont = .false.
       else if((direc == (g1+g2+g3)).and.  (g3 == 1)) then
          dy = 1
          tbdi(maxl-1) = 2
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == 2) reptcont = .false.
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          tbdi(maxl-1) = -2
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == -2) reptcont = .false.
       else if((direc == (g1+g2+g3+g4+g5)) .and.( g5 == 1)) then
          dz = 1
          tbdi(maxl-1) = 3
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == 3) reptcont = .false.
       else if((direc == (g1+g2+g3+g4+g5+g6)).and.  (g6 == 1)) then
          dz = -1
          tbdi(maxl-1) = -3
          call bonddirection(5,-tbdi(maxl-1),8,tbo(maxl))
          if(bonddd(m,maxl) == -3) reptcont = .false.
       end if
       if(reptcont .eqv. .false.) goto 83


       tempcoord(maxl)%x = modulo(protcoords(m,maxl)%x + dx-1,gridsize)+1
       tempcoord(maxl)%y = modulo(protcoords(m,maxl)%y + dy-1,gridsize)+1
       tempcoord(maxl)%z = modulo(protcoords(m,maxl)%z + dz-1,gridsize)+1

       !write(6,*) tbdi(maxl-1),dx,dy,dz
       !write(6,*) tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z
       !write(6,*) protcoords(m,maxl)%x,protcoords(m,maxl)%y,protcoords(m,maxl)%z




       deltaenergy = 0.0d0


       do st = 2,maxl,1
          if(overlaps(m,maxl,st,tempcoord) .eqv. .false.) then
             reptcont = .false.
             goto 83
          end if
       end do


       do pr =1,nprotein
          if(pr /=m)then
        maxback = chlen(pr)
             !write(6,*) pr,maxback,maxl
             do st = 1,maxback,1
                if(overlaps(pr,maxl,st,tempcoord) .eqv. .false.) then
                   reptcont = .false.
                   goto 83
                end if
             end do
          end if
       end do

       do st = 4,maxl,1
          if(isobond(m,1,st) == 1) deltaenergy = deltaenergy - intraenergy(protcoords(m,st)%type)
          if(ba(m,1,st) == 1) deltaenergy = deltaenergy - intraen
       end do

       do st = 2,maxl-2
          call adjacent(m,maxl,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             deltaenergy = deltaenergy + intraen
             tba(st-1,maxl) = 1
             if((bdir == -1*tbo(maxl)) .and. (bdir ==bonddd(m,st))) then
                deltaenergy = deltaenergy + intraenergy(protcoords(m,st)%type)
                tempisobond(st-1,maxl) = 1
             else
                tempisobond(st-1,maxl) = 0
             end if
          else if((adjver.eqv. .false.)) then
             tba(st-1,maxl) = 0
             tempisobond(st-1,maxl) = 0
          end if
       end do

       do pr = 1,m-1,1
         maxback = chlen(pr)
          do st = 1,maxback,1
             if(intbond(pr,m,st,1) == 1) deltaenergy = deltaenergy - &
                  interenergy(protcoords(m,1)%type,protcoords(pr,st)%type)
             if(br(pr,m,st,1) == 1) deltaenergy = deltaenergy - interen(protcoords(m,1)%type,protcoords(pr,st)%type)
             call adjacent(pr,maxl,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + interen(protcoords(m,maxl)%type,protcoords(pr,st)%type)
                tbr(pr,maxl,st) = 1
                if((bdir == (-1*tbo(maxl))) .and. (bdir ==bonddd(pr,st))) then
                   deltaenergy = deltaenergy + &
                        interenergy(protcoords(m,maxl)%type,protcoords(pr,st)%type)
                   !if(time == 3387) write(6,*) 'bond',pr,st,m
                   tempintbond(pr,maxl,st) = 1
                else
                   tempintbond(pr,maxl,st) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tbr(pr,maxl,st) = 0
                tempintbond(pr,maxl,st) = 0
             end if
          end do
       end do

       do pr = m+1,nprotein,1
        maxback = chlen(pr)
          do st = 1,maxback,1
             if(intbond(m,pr,1,st) == 1) deltaenergy = deltaenergy - &
                  interenergy(protcoords(m,1)%type,protcoords(pr,st)%type)
             if(br(m,pr,1,st) == 1) deltaenergy = deltaenergy -interen(protcoords(m,1)%type,protcoords(pr,st)%type)
             call adjacent(pr,maxl,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                deltaenergy = deltaenergy + interen(protcoords(m,maxl)%type,protcoords(pr,st)%type)
                tbr(pr,maxl,st) = 1
                if((bdir == (-1*tbo(maxl))) .and. (bdir ==bonddd(pr,st))) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,maxl)%type,protcoords(pr,st)%type)
                   tempintbond(pr,maxl,st) = 1
                   !if(time == 3387) write(6,*) 'bond',pr,st,m
                else
                   tempintbond(pr,maxl,st) = 0
                end if
             else if((adjver.eqv. .false.)) then
                tbr(pr,maxl,st) = 0
                tempintbond(pr,maxl,st) = 0
             end if
          end do
       end do


       do str = 2,maxl,1
          do pr = 1,m-1,1
             maxback = chlen(pr)
             do st = 1,maxback
                if (intbond(pr,m,st,str) ==1) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,str-1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interenergy(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if
                
                if (br(pr,m,st,str) ==1) then
                   deltaenergy = deltaenergy + interen(protcoords(m,str-1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interen(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if
             end do
          end do
       end do

       do str = 2,maxl,1
          do pr = m+1,nprotein,1
             maxback = chlen(pr)
             do st = 1,maxback
                if (intbond(m,pr,str,st) ==1) then
                   deltaenergy = deltaenergy + interenergy(protcoords(m,str-1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interenergy(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if
                
                if (br(m,pr,str,st) ==1) then
                   deltaenergy = deltaenergy + interen(protcoords(m,str-1)%type,protcoords(pr,st)%type)
                   deltaenergy = deltaenergy - interen(protcoords(m,str)%type,protcoords(pr,st)%type)
                end if
             end do
          end do
       end do

       
    !   write(6,*) 'delta energy', deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          do l =1,maxl-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z
             bonddd(m,l) = bonddd(m,l+1)
             if(l<maxl-1) bdi(m,l) = bdi(m,l+1)
             do g= l+3,maxl-1,1
                isobond(m,l,g) = isobond(m,l+1,g+1)
                ba(m,l,g) = ba(m,l+1,g+1)
             end do

             do chain2 = 1,m-1,1
        maxback = chlen(chain2)
                do beads2 = 1,maxback
                   intbond(chain2,m,beads2,l) = intbond(chain2,m,beads2,l+1)
                   br(chain2,m,beads2,l) = br(chain2,m,beads2,l+1)
                end do
             end do

             do chain2 = m+1,nprotein,1
        maxback = chlen(chain2)
                do beads2 = 1,maxback
                   intbond(m,chain2,l,beads2) = intbond(m,chain2,l+1,beads2)
                   br(m,chain2,l,beads2) = br(m,chain2,l+1,beads2)
                end do
             end do

          end do


          tempcoord(maxl-1)%x = protcoords(m,maxl-1)%x
          tempcoord(maxl-1)%y = protcoords(m,maxl-1)%y
          tempcoord(maxl-1)%z = protcoords(m,maxl-1)%z
          tbo(maxl-1) = bonddd(m,maxl-1)
          !tbdi(maxl-1) = bdi(m,maxl-1)

          call updatepos(m,maxl-1,maxl,tempcoord,tbo,tbdi,maxl)
          call updateintrabond(m,maxl,maxl,1,maxl-3,tempisobond,tba)
          call updateinterbond(m,maxl,maxl,tempintbond,tbr)
          !write(6,*) 'backwards rep!!!!!!!!!!!!'
          successful = successful + 1
          reptbackward = reptbackward + 1
          call counts(0,0,0,1)
       else
          reptcont = .false.
       end if
    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
       call counts(0,0,0,2)
    end if

  end subroutine reptation


  subroutine bonddirection(ndir,no1,no2,tempbonddd)
    integer,intent(in) :: ndir,no1,no2
    integer,intent(inout) :: tempbonddd
    integer :: xp,xn,yp,yn,zp,zn,chdir 

    chdir = int(ran2(seed)*ndir)+1

    xp = 1
    xn = 1
    yp = 1
    yn = 1
    zp = 1
    zn = 1

    if(no1 == 1 .or. no2 ==1) xp = 0 
    if(no1 == 2 .or. no2 ==2) yp = 0
    if(no1 == 3 .or. no2 ==3) zp = 0
    if(no1 == -3 .or. no2 ==-3) zn = 0
    if(no1 == -2 .or. no2 ==-2) yn = 0
    if(no1 == -1 .or. no2 ==-1) xn = 0

    if((chdir == xp) .and. (xp ==1)) tempbonddd = 1
    if( (chdir == xp + yp) .and. (yp ==1)) tempbonddd = 2
    if( (chdir == xp + yp+ zp) .and. (zp ==1)) tempbonddd = 3
    if((chdir == xp + yp + zp+zn) .and. (zn ==1)) tempbonddd = -3
    if((chdir == xp + yp + zp+zn+yn) .and. (yn ==1)) tempbonddd = -2
    if((chdir == ndir) .and. (xn ==1)) tempbonddd = -1

  end subroutine bonddirection

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


            do a =1,maxl
       !write(6,*) 'bead',a,protcoords(m,a)%x,protcoords(m,a)%y,protcoords(m,a)%z
       end do
       do l = 2,maxl,1
          dcomx = protcoords(m,l)%x-protcoords(m,l-1)%x
          if (dcomx == gridsize-1) dcomx = -1
          if (dcomx == 1-gridsize) dcomx = 1
          comx = comx + ((maxl+1-l)*dcomx)
          dcomy = protcoords(m,l)%y-protcoords(m,l-1)%y
          if (dcomy == gridsize-1) dcomy = -1
          if (dcomy == 1-gridsize) dcomy = 1
          comy = comy + ((maxl+1-l)*dcomy)
          dcomz = protcoords(m,l)%z-protcoords(m,l-1)%z
          if (dcomz == gridsize-1) dcomz = -1
          if (dcomz == 1-gridsize) dcomz = 1
          comz = comz + ((maxl+1-l)*dcomz)
       end do
       !write(6,*) 'deltas',comx,comy,comz
       com(1,m)%x = modulo(protcoords(m,1)%x +(real(comx)/(maxl)) -1,real(gridsize))+1
       com(1,m)%y = modulo(protcoords(m,1)%y +(real(comy)/(maxl)) -1,real(gridsize))+1
       com(1,m)%z = modulo(protcoords(m,1)%z +(real(comz)/(maxl)) -1,real(gridsize))+1
    !end do
           !write(6,*) 'com',com(1,m)%x,com(1,m)%y,com(1,m)%z
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
    type(rprot) :: dummycom

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
             if (delx > (maxb+maxltest)/2) delx = delx-gridsize
             if (delx < -(maxb+maxltest)/2) delx = delx + gridsize
             comx(b,f) =  delx
             comx(f,b) = -delx
             dely = com(1,f)%y - com(1,b)%y
             if (dely > (maxb+maxltest)/2) dely = dely-gridsize
             if (dely < -(maxb+maxltest)/2) dely = dely + gridsize
             comy(b,f) =  dely
             comy(f,b) = -dely
             delz = com(1,f)%z - com(1,b)%z
             if (delz > (maxb+maxltest)/2) delz = delz-gridsize
             if (delz < -(maxb+maxltest)/2) delz = delz + gridsize
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
16           if(a == a) continue
          end do

       end do

           
    !write(6,*) 'backup 2',rcounter(4)
!write(6,*) 'orig central chain',centralchain,toolarge
   call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
!write(6,*) 'key'

    
    !write(6,*) 'backup 3',rcounter(4)

76  if(clsize == 1) then
      maxl = chlen(m)
    centralchain = m
    totcluspop = maxl
   ! write(6,*) 'individual chain',com(1,m)%x,com(1,m)%y,com(1,m)%z
    end if

    !if(clsize>1) write(6,*) 'totcluspop',totcluspop
    !write(6,*) centralchain
    !write(6,*) 'runtime', com(1,centralchain)%x
    !write(6,*) dummycom%x
comcluster%x = NINT(modulo(com(1,centralchain)%x + (dummycom%x/totcluspop)-1,real(gridsize))+1)

comcluster%y = NINT(modulo(com(1,centralchain)%y + (dummycom%y/totcluspop)-1,real(gridsize))+1)

comcluster%z = NINT(modulo(com(1,centralchain)%z + (dummycom%z/totcluspop)-1,real(gridsize))+1)



53  if(compass .eqv. .true.) continue

  end subroutine clustercom


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


    if((maxx - minx) > real(gridsize/2)) toolarge = .true.
    if((maxy - miny) > real(gridsize/2)) toolarge = .true.
    if((maxz - minz) > real(gridsize/2)) toolarge = .true.
    
    !do a = 1,clsize
       !b= checklist(a)
       !write(6,*) 'individual',b,com(1,b)%x,com(1,b)%y,com(1,b)%z
   ! end do
    
  end subroutine comlargecomp


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
          if (delx > gridsize/2) delx = delx-gridsize
          if (delx < -gridsize/2) delx = delx + gridsize
          comx = comx + (delx*maxb)
          dely = com(1,a)%y -com(1,m)%y
          if (dely > gridsize/2) dely = dely-gridsize
          if (dely < -gridsize/2) dely = dely + gridsize
          comy = comy + (dely*maxb)
          delz = com(1,a)%z-com(1,m)%z
          if (delz > gridsize/2) delz = delz-gridsize
          if (delz < -gridsize/2) delz = delz + gridsize
          comz = comz + (delz*maxb)
          totcluspop = totcluspop + maxb
       end if
    end do
       
    
    dummycom%x = modulo(com(1,m)%x + (comx/totcluspop)-1,real(gridsize))+1
    comcluster%x = NINT(dummycom%x)
    dummycom%y = modulo(com(1,m)%y + (comy/totcluspop)-1,real(gridsize))+1
    comcluster%y = NINT(dummycom%y)
    dummycom%z = modulo(com(1,m)%z + (comz/totcluspop)-1,real(gridsize))+1
    comcluster%z = NINT(dummycom%z)
    
    
  end subroutine clustercomcheck

            
  logical Function overlaps (pr,ch1,ch2,tempcoords)
    integer,intent(in) :: pr,ch1,ch2
    Type(basicp),dimension(:),intent(in) :: tempcoords
    overlaps = .true.
    if((protcoords(pr,ch2)%x == tempcoords(ch1)%x) .and. &
         (protcoords(pr,ch2)%y == tempcoords(ch1)%y) .and. &
         (protcoords(pr,ch2)%z == tempcoords(ch1)%z)) then
       overlaps = .false.
    end if
  end Function Overlaps

  logical Function overlapsclus (m,pr,ch1,ch2,tempcluscoord)
    integer,intent(in) :: pr,ch1,ch2,m
    Type(basicp),dimension(:,:),intent(in) :: tempcluscoord
    overlapsclus = .true.
    if((protcoords(pr,ch2)%x == tempcluscoord(m,ch1)%x) .and. &
         (protcoords(pr,ch2)%y == tempcluscoord(m,ch1)%y) .and. &
         (protcoords(pr,ch2)%z == tempcluscoord(m,ch1)%z)) then
       overlapsclus = .false.
    end if
  end Function Overlapsclus

  subroutine energy
    integer :: dx,dy,dz,m,l,f,g,bdir,delx,dely,delz,maxlengthss,maxl
    double precision :: initialenergy
    Type(basicp),dimension(:),allocatable :: tempcoord
    logical :: adjver
    initialenergy = 0.0d0
    !totalenergy = 0.0d0
   allocate(tempcoord(maxlength))
    do m= 1,nprotein,1
        maxl = chlen(m)
       do l = 1,maxl-3,1
          do f = l+3,maxl,1
             if((f < l-2) .or. (f > l + 2)) then
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z
                call adjacent(m,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.)) then
                   if(time == 0) ba(m,l,f) = 1
                   initialenergy = initialenergy + intraen
                   if((bdir == bonddd(m,f)) .and. (bdir == -1*bonddd(m,l))) then
                      if(time == 0) isobond(m,l,f) = 1
                      initialenergy = initialenergy + intraenergy(protcoords(m,l)%type)
                   else
                      if(time == 0) isobond(m,l,f) = 0
                   end if
                else if((adjver .eqv. .false.)) then
                   if(time == 0) isobond(m,l,f) = 0
                   if(time == 0) ba(m,l,f) = 0
                end if
             end if

          end do
       end do
    end do


!!!!!!! interatomic
    do m= 1,nprotein-1,1
        maxl = chlen(m)
       do g = m+1,nprotein,1
                  maxlengthss = chlen(g)
          do l = 1,maxl,1
             do f = 1,maxlengthss,1
                !if(g /=m) then

                   tempcoord(l)%x = protcoords(m,l)%x
                   tempcoord(l)%y = protcoords(m,l)%y
                   tempcoord(l)%z = protcoords(m,l)%z



                   call adjacent(g,l,f,tempcoord,adjver,bdir)
                   if((adjver .eqv. .true.)) then
                      if(time == 0) br(m,g,l,f) = 1
                      !if(isbond .eqv. .true.) then
                         initialenergy = initialenergy + interen(protcoords(m,l)%type,protcoords(g,f)%type)
                      !end if
                      if((bdir == bonddd(g,f)) .and. (bdir == (-1*bonddd(m,l)))) then
                         if(time == 0) intbond(m,g,l,f) = 1
                         initialenergy = initialenergy + interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                         !write(6,*) 'energy check',interen(protcoords(m,l)%type,protcoords(g,f)%type),&
                              !interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                      else
                         if(time == 0) intbond(m,g,l,f) = 0
                      end if
                   else if((adjver .eqv. .false.)) then
                      br(m,g,l,f) = 0
                       if(time == 0) intbond(m,g,l,f) = 0
                   end if

                !end if
             end do
          end do
       end do
    end do

    if((time > 0) .and. (initialenergy /= totalenergy)) then
       write(6,*) 'energyfail',totalenergy,initialenergy,time
       finalfail = .true.
    end if
    !if(modulo(time,1000) == 0 .and. (noinfo .eqv. .false.)) write(93,*) 'manual check',time,totalenergy
    !if(time == 49021) write(6,*) 'energy =',totalenergy,time
    if(time == 0) totalenergy = initialenergy


  end subroutine energy

  subroutine radiusofgyration
    integer::m,l,maxl
    double precision :: rog,totrog
totrog = 0.0d0
    do m = 1,nprotein,1
        maxl = chlen(m)
       do l = 1, maxl,1
          rog = (min(modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0))))**2
          totrog = totrog + rog/maxl
       end do
    end do

    runningaveROG = runningaveROG + totrog
    polymerrog = polymerrog + (totrog/nprotein)
    if(noinfo .eqv. .false.) write(97,*) time-equilib, &
         totrog/(nprotein*maxlength),SQRT(totrog/(nprotein*maxlength))

  end subroutine radiusofgyration


  
  subroutine adjacent(pr,ch1,ch2,tempory,adjver,bdir)
    integer,intent(in):: pr,ch1,ch2
    integer,intent(inout)::bdir
    integer::dx,dy,dz,delx,dely,delz,signs
    logical,intent(inout)::adjver
    Type(basicp),dimension(:),intent(in) :: tempory
    dx = 4
    dy = 4
    dz = 4

    delx = tempory(ch1)%x - protcoords(pr,ch2)%x
    if((abs(delx) == 1) .or. (abs(delx) == (gridsize - 1))) then
       dx = 1
       bdir = delx
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

       !msdx = (min(modulo(com(1,m)%x - com(2,m)%x,gridsize*1.0), &
        !    (gridsize - modulo(com(1,m)%x -com(2,m)%x,gridsize*1.0))))**2
       !msdy = (min(modulo(com(1,m)%y -com(2,m)%y,gridsize*1.0), &
       !     (gridsize - modulo(com(1,m)%y -com(2,m)%y,gridsize*1.0))))**2
       !msdz = (min(modulo(com(1,m)%z - com(2,m)%z,gridsize*1.0), &
       !     (gridsize - modulo(com(1,m)%z -com(2,m)%z,gridsize*1.0))))**2
       
       !total = msdx + msdy + msdz
       
       !if(time > 42000 .and. (time < 47000) .and. total > 1.0) then
          !write(6,*) 'MSD', total,time
       !write(6,*) 'x com',com(1,m)%x,com(2,m)%x,msdx
          !write(6,*) 'y com',com(1,m)%y,com(2,m)%y,msdy
          !write(6,*) 'z com',com(1,m)%z,com(2,m)%z,msdz
          !end if
       
       
       totrmsbrute = totrmsbrute + msdx + msdy + msdz

  end subroutine rms


  subroutine dataout
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit

    allocate(bondposit(nprotein,maxlength))
    !write(6,*) 'data'
    t = time
    write(67,*) ' '
    write(67,*) 'timestep ', 'indexed'
    write(67,*) 'pbc', real(4*gridsize),real(4*gridsize),real(4*gridsize)
    !write(67,*) ((nprotein1*maxlength1) + (nprotein2*maxlength2))
    acount = 0

    do m = 1,nprotein,1
        maxl = chlen(m)
       do l = 1,maxl,1
          write(67,*) acount,real(4*protcoords(m,l)%x), real(4*protcoords(m,l)%y),real(4*protcoords(m,l)%z) 
          acount = acount +1

          bondposit(m,l)%x = real(protcoords(m,l)%x)
          bondposit(m,l)%y = real(protcoords(m,l)%y)
          bondposit(m,l)%z = real(protcoords(m,l)%z)

          if(bonddd(m,l) == 1) bondposit(m,l)%x = modulo(bondposit(m,l)%x + 0.5,real(gridsize))
          if(bonddd(m,l) == -1) bondposit(m,l)%x = modulo(bondposit(m,l)%x - 0.5,real(gridsize))
          if(bonddd(m,l) == 2) bondposit(m,l)%y = modulo(bondposit(m,l)%y + 0.5,real(gridsize))
          if(bonddd(m,l) == -2) bondposit(m,l)%y = modulo(bondposit(m,l)%y - 0.5,real(gridsize))
          if(bonddd(m,l) == 3) bondposit(m,l)%z = modulo(bondposit(m,l)%z + 0.5,real(gridsize))
          if(bonddd(m,l) == -3) bondposit(m,l)%z = modulo(bondposit(m,l)%z - 0.5,real(gridsize))


          do f= 1,l-1
             if((isobond(m,f,l) ==1)) then
                write(67,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
                4*bondposit(m,l)%z
                goto 71           
             end if
          end do

          do f= l+1,maxl
             if((isobond(m,f,l) ==1)) then
                write(67,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
                4*bondposit(m,l)%z
                goto 71           
             end if
          end do
          do g = m+1,nprotein
        maxback = chlen(g)
             do f= 1,maxback
                if((intbond(m,g,l,f) ==1)) then
                   write(67,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
                   4*bondposit(m,l)%z
                   goto 71
                end if
             end do
          end do

    do g = 1,m-1
        maxback = chlen(g)
             do f= 1,maxback
                if((intbond(g,m,f,l) ==1)) then
                   write(67,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
                   4*bondposit(m,l)%z
                   goto 71
                end if
             end do
          end do

          
          write(67,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
          4*bondposit(m,l)%z

71        if(m == m) continue
          acount = acount + 1
       end do
    end do



  end subroutine dataout


  subroutine clusterenergy(tempcluscoord,clno,clusen,m,ctbo,clsize,cllist)
    integer,dimension(:),intent(in) ::cllist,clno
    integer,dimension(:,:),intent(inout)::ctbo
    integer,intent(in) :: clsize,m
    double precision,intent(inout)::clusen
    type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    type(basicp),dimension(:),allocatable :: tempcoord
    integer:: a,l,g,f,bdir,maxlengthss,maxl,ra
    logical :: clusmove,adjver
    allocate(tempcoord(maxlength))

    clusen = 0.0
!write(6,*) 'clsize',clsize
    do ra = 1,clsize
       a = cllist(ra)
      ! write(6,*) 'a',a
        maxl = chlen(a)
       do g = 1,nprotein
          if(clno(a) /= clno(g)) then
             maxlengthss = chlen(g)

             clusmove = .true.
             do l = 1,maxl,1
                tempcoord(l)%x = tempcluscoord(a,l)%x
                tempcoord(l)%y = tempcluscoord(a,l)%y
                tempcoord(l)%z = tempcluscoord(a,l)%z

                do f = 1,maxlengthss,1
                   call adjacent(g,l,f,tempcoord,adjver,bdir)
                   if(adjver.eqv. .true.) then
                      if((isbond .eqv. .true.) .and. (interen(protcoords(a,l)%type,protcoords(g,f)%type) < 0.0)) then
                         !energyofcluster = energyofcluster + interen
                         clusen = 100
                         return
                      end if
                      if((bdir == (-1*ctbo(a,l))) .and. (bdir == bonddd(g,f)) .and. &
                           (interenergy(protcoords(a,l)%type,protcoords(g,f)%type)  < 0.0)) then
                         !energyofcluster = energyofcluster + interenergy(protcoords(a,l)%species,protcoords(g,f)%species)
                         clusen = 100
                         return
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do
    !if(time> 1490 .and. (time < 1510)) write(6,*) 'clusterenergy =',energyofcluster
    !write(6,*) 'aaaaaa'
  end subroutine clusterenergy


  subroutine length
    integer :: ddx,ddy,ddz,dxsum,dysum,dzsum,m,l,maxl
    double precision :: normchainl
    t = time
    totchainlength = 0.0
    normchainl = 0.0d0

    do m = 1,nprotein,1
       !chainlength = 0.0
       ddx = 0
       ddy = 0
       ddz = 0
       dxsum = 0
       dysum = 0
       dzsum = 0
        maxl = chlen(m)
       do l = 2,maxl,1

          ddx = protcoords(m,l)%x-protcoords(m,l-1)%x
          if (ddx > 1) ddx =-1
          if (ddx < -1) ddx = 1


          dxsum = dxsum + ddx
          ddy = protcoords(m,l)%y-protcoords(m,l-1)%y
          if (ddy > 1) ddy =-1
          if (ddy < -1) ddy = 1

          dysum = dysum + ddy
          ddz = protcoords(m,l)%z-protcoords(m,l-1)%z
          if (ddz > 1) ddz =-1
          if (ddz < -1) ddz = 1

          dzsum = dzsum + ddz
       end do
       chainlength = ((dxsum**2) + (dysum**2) + (dzsum**2))
       totchainlength = totchainlength + chainlength
    end do
    !totchainlength = sum(chainlength)
    runningaveEtE = runningaveEtE + totchainlength  !sort this out
    normchainl = totchainlength/nprotein
polymerl = polymerl + normchainl
    actualchainlength = actualchainlength + avechainlength
    !write(91,*) t,avechainlength,actualchainlength/((time-equilib)/10000)
  end subroutine length

  subroutine clusrbond(m,ctbo,axis,rotor,maxl)
    integer,intent(in):: m,axis,rotor,maxl
    integer,dimension(:,:),intent(inout)::ctbo
    logical :: bm
    integer :: sign,s


    do s = 1,maxl
       bm = .true.
       if(abs(bonddd(m,s)) == axis) then
          bm = .true.
          ctbo(m,s) = bonddd(m,s)
          goto 87
       end if


       if(rotor == 0) then
          bm = .true.
          ctbo(m,s) = -1*bonddd(m,s)
          !if(abs(bonddd(m,s)) /= (modulo(axis,3) + 1)) ctbo(m,s) = 1*bonddd(m,s)
          goto 87
       end if
       !write(6,*) bonddd(m,s),m,s
       sign = bonddd(m,s)/abs(bonddd(m,s))
       if(axis == 1) then
          if(rotor == 1) then
             if(abs(bonddd(m,s)) == 2) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 1)
             if(abs(bonddd(m,s)) == 3) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 1)
             goto 87
          else if(rotor == -1) then
             if(abs(bonddd(m,s)) == 2) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 1)
             if(abs(bonddd(m,s)) == 3) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 1)
             goto 87
          end if
       else if(axis ==2 ) then
          if(rotor == 1) then
             if(abs(bonddd(m,s)) == 1) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 2)
             if(abs(bonddd(m,s)) == 3) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 2)        
             goto 87
          else if(rotor == -1) then
             if(abs(bonddd(m,s)) == 1) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 2)
             if(abs(bonddd(m,s)) == 3) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 2)         
             goto 87
          end if

       else if(axis ==3) then
          if(rotor == 1) then
             if(abs(bonddd(m,s)) == 1) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 1)
             if(abs(bonddd(m,s)) == 2) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 1)        
             goto 87
          else if(rotor == -1) then
             if(abs(bonddd(m,s)) == 1) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 1)
             if(abs(bonddd(m,s)) == 2) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 1)         
          end if

       end if


87     if(bm .eqv. .true.) continue
       if(ctbo(m,s) == 0) write(6,*) 'fail cluster bond',m,s
    end do
  end subroutine clusrbond

  subroutine clustertranslation(m,clno,clsize,cb,cllist)
    integer,intent(inout)::m
    integer,intent(inout):: clsize
    type(centremass),dimension(:),allocatable:: tcom
    integer,dimension(:),intent(inout) ::clno,cllist
    integer,dimension(:,:),intent(inout) ::cb
    integer :: steplength,dx,dy,dz,vv,rvv,a
    type(basicp),dimension(:,:),allocatable :: tempcluscoord
    integer,dimension(:,:),allocatable::ctbo,ctbdi
    double precision :: clusen,direction
    logical :: moveallow
    allocate(tempcluscoord(nprotein,maxlength))
    allocate(ctbo(nprotein,maxlength))
    allocate(ctbdi(nprotein,maxlength-1))
    allocate(tcom(nprotein))
    !return
   !write(6,*) 'translate'
    direction = ran2(seed)
    steplength = 1
    moveallow = .true.
    if(direction <= 1.0/6) then
       dx = steplength
       dy = 0
       dz = 0
    else if((direction > 1.0/6) .and. (direction <= 1.0/3)) then
       dx = -(steplength)
       dy = 0
       dz = 0
    else if((direction > 1.0/3) .and. (direction <= 1.0/2)) then
       dx = 0
       dy = steplength
       dz = 0
    else if((direction > 1.0/2) .and. (direction <= 2.0/3)) then
       dx = 0
       dy = -(steplength)
       dz = 0
    else  if((direction > 2.0/3) .and. (direction <= 5.0/6)) then
       dx = 0
       dy = 0
       dz = (steplength)
    else if((direction > 5.0/6) .and. (direction <= 1.0/1)) then
       dx = 0
       dy = 0
       dz = -steplength
    end if
! write(6,*) 'hereheee'
    do rvv = 1,clsize
       vv = cllist(rvv)
       !write(6,*) rvv,vv,clsize
          !if(time > 1750 .and. time<1760) write(6,*) time,'b'
          call clustermove(clno,vv,dx,dy,dz,tempcluscoord,moveallow,ctbo,ctbdi)
          if(moveallow .eqv. .false.) return
          tcom(vv)%x = modulo(com(1,vv)%x + dx -1,real(gridsize))+1
          tcom(vv)%y = modulo(com(1,vv)%y + dy -1,real(gridsize))+1
          tcom(vv)%z = modulo(com(1,vv)%z + dz -1,real(gridsize))+1
!write(6,*) tcom(vv)%z
    end do
 !write(6,*) 'pass'
    !if(moveallow .eqv. .true.) then
    clusen = 0.0
 

    
    !write(6,*) 'a'
    call clusterenergy(tempcluscoord,clno,clusen,m,ctbo,clsize,cllist)
    !write(6,*) 'translation cluster energy', clusen
    if(clusen == 0.0) then
       call updatecluspos(tempcluscoord,ctbo,ctbdi,tcom,cllist,clsize)
       call updateclusbonds(clsize,cllist,clno)
       cltacc = cltacc + 1
       end if
       !write(6,*) 'b'


  end subroutine clustertranslation


  
  subroutine clusterrotation(m,clno,clsize,cb,cllist)
    integer,intent(inout) :: m
    integer,intent(inout)::clsize
    integer,dimension(:),intent(inout) ::clno,cllist
    integer,dimension(:,:),intent(inout) ::cb
    integer,dimension(:,:),allocatable:: route
    integer,dimension(:),allocatable::rcounter,checklist
    integer :: a,b,l,pr,st,maxl,maxback,dummy,rotatechoose,xr,yr,zr,rb,chnum,totcluspop
    double precision :: clusen
    integer,dimension(:,:),allocatable :: deltax,deltay,deltaz
    Type(basicp),dimension(:,:),allocatable :: tempcluscoord
    type(centremass),dimension(:),allocatable:: tcom
    integer,dimension(:,:),allocatable::ctbo,ctbdi
    logical :: clusmove,clusrotpass,toolarge
    type(basicp) :: comcluster,comclustemp,mockcomcluster
    clusmove = .true.

    allocate(ctbo(nprotein,maxlength))
    allocate(ctbdi(nprotein,maxlength-1))
    allocate(tempcluscoord(nprotein,maxlength))
    allocate(tcom(nprotein))
    allocate(deltax(nprotein,maxlength))
    allocate(deltay(nprotein,maxlength))
    allocate(deltaz(nprotein,maxlength))
    allocate(route(nprotein,clsize))
    allocate(rcounter(nprotein))
    allocate(checklist(clsize))
  
    t = time + 1
    !write(6,*) 'rotate'

     do b = 1,clsize
rb = cllist(b)
           call comfind(rb,.false.)
     end do
toolarge = .false.
!write(6,*) 'clsize',clsize
if((clsize*maxlength)>gridsize) then

   call clustercom(m,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge)
      !write(6,*) 'b',toolarge
        if(toolarge .eqv. .true.) return 
else
   call clustercomcheck(m,clsize,comcluster,cb,cllist)
end if



!if(comcluster%x > gridsize) comcluster%x = modulo(comcluster%x,gridsize)
!if(comcluster%y > gridsize) comcluster%y = modulo(comcluster%y,gridsize)
!if(comcluster%z > gridsize) comcluster%z = modulo(comcluster%x,gridsize)

    !if((clsize*maxlength)>gridsize) then
                   ! call clustercomcheck(m,clsize,comclustemp,cb,cllist)
                 !if(((comclustemp%y /= comcluster%y)))  write(6,*) 'fail\\\\\\\\\\\\\\\\\\\\\\\\\\'
 !end if
     !if(clsize>1) then
        !write(6,*) 'new com 1', comcluster%x,comcluster%y,comcluster%z
        !write(6,*) 'old com 2', mockcomcluster%x,mockcomcluster%y,mockcomcluster%z    
!end if
     rotatechoose =  int(ran2(seed)*9)+1

     do rb = 1,clsize
        b = cllist(rb)
     
        maxl = chlen(b)
          do l = 1, maxl
             deltax(b,l) =    protcoords(b,l)%x -comcluster%x 
             if(deltax(b,l) > gridsize/2) deltax(b,l) =  deltax(b,l)-gridsize 
             if(deltax(b,l) < -gridsize/2) deltax(b,l) = gridsize + deltax(b,l)
             deltay(b,l) =    protcoords(b,l)%y - comcluster%y
             if(deltay(b,l) > gridsize/2) deltay(b,l) =  deltay(b,l)-gridsize 
             if(deltay(b,l) < -gridsize/2) deltay(b,l) = gridsize + deltay(b,l) 
             deltaz(b,l) =    protcoords(b,l)%z - comcluster%z 
             if(deltaz(b,l) > gridsize/2) deltaz(b,l) =  deltaz(b,l)-gridsize
             if(deltaz(b,l) < -gridsize/2) deltaz(b,l) = gridsize + deltaz(b,l) 
          end do

          if(rotatechoose == 1) then
             call rotatex(b,1,-1,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 2) then
             call rotatex(b,-1,1,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 3 ) then
             call clusterflip(b,1,-1,-1,deltax,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)    
          else if(rotatechoose == 4) then
             call rotatey(b,1,-1,deltax,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 5) then
             call rotatey(b,-1,1,deltax,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 6) then
             call clusterflip(b,-1,1,-1,deltax,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 7) then
             call rotatez(b,1,-1,deltax,deltay,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 8) then
             call rotatez(b,-1,1,deltax,deltay,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          else if(rotatechoose == 9) then
             call clusterflip(b,-1,-1,1,deltax,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
          end if
          
          do pr = 1,nprotein
             if(clno(pr) /= clno(m)) then
                do l = 1,maxl
                   maxback = chlen(pr)
                   do st = 1,maxback,1
                      if(overlapsclus(b,pr,l,st,tempcluscoord) .eqv. .false.) then
                         clusmove = .false.
                         goto 31
                      end if
                   end do
                end do
                end if
          end do
    end do
      

    do rb = 1,clsize
       b = cllist(rb)
          call tempcomfind(b,tcom,tempcluscoord)
       end do

  
if((clsize*maxlength)>gridsize) then
       call tempclustercom(m,clsize,tcom,comclustemp,cb,cllist,route,rcounter,totcluspop,checklist)
    else
       call shorttempclustercom(m,clsize,tcom,comclustemp,cb,cllist)
    end if

    !if(abs(comclustemp%x - comcluster%x)> 2 .and. (abs(comclustemp%x - comcluster%x)<gridsize-1)) then
       !write(6,*) 'rotatechoose',rotatechoose
       !write(6,*) 'delta', protcoords(m,5)%x -comcluster%x, deltax(m,5),protcoords(m,5)%x,tempcluscoord(m,5)%x
       !write(6,*) 'com x',(comclustemp%x - comcluster%x),comclustemp%x,comcluster%x,com(1,m)%x,tcom(m)%x
       !write(6,*) 'com y',(comclustemp%y - comcluster%y),comclustemp%y,comcluster%y,com(1,m)%y
       !write(6,*) 'com z',(comclustemp%z - comcluster%z),comclustemp%z,comcluster%z,com(1,m)%z
    !end if
 

       xr = 0
       yr = 0
       zr = 0
!if(comclustemp%x > gridsize) comclustemp%x = modulo(comclustemp%x,gridsize)
!if(comclustemp%y > gridsize) comclustemp%y = modulo(comclustemp%y,gridsize)
!if(comclustemp%z > gridsize) comclustemp%z = modulo(comclustemp%x,gridsize)

      !write(6,*) 'change in com',(comclustemp%x - comcluster%x)
            !(comclustemp%y - comcluster%y),(comclustemp%z - comcluster%z)
       if(((comclustemp%x - comcluster%x) == 1) .or.((comclustemp%x - comcluster%x)== 1-gridsize)) xr = -1
       if(((comclustemp%y - comcluster%y) == 1) .or. ((comclustemp%y - comcluster%y)== 1-gridsize)) yr = -1
       if(((comclustemp%z - comcluster%z) == 1) .or. ((comclustemp%z - comcluster%z)== 1-gridsize)) zr = -1



       
clusrotpass = .true.
       if((xr+yr+zr) /= 0) &
            call rotcorrect(m,clno,clsize,xr,yr,zr,tempcluscoord,tcom,cllist,clusrotpass)
       
if(clusrotpass .eqv. .false.) return


       clusen = 0.0

    call clusterenergy(tempcluscoord,clno,clusen,m,ctbo,clsize,cllist)
    
    if(clusen == 0.0) then
       call updatecluspos(tempcluscoord,ctbo,ctbdi,tcom,cllist,clsize)
       call updateclusbonds(clsize,cllist,clno)
       clracc = clracc + 1
    end if


    !totcluspop = 0
         !call clustercom(m,clsize,comclustemp,cb,cllist,route,rcounter,totcluspop,checklist,toolarge)
         !if((abs(comclustemp%y - comcluster%y)>1) .and. (abs(comclustemp%y - comcluster%y)<gridsize-1)) then
            !write(6,*) 'yfail actual', (comclustemp%y - comcluster%y),comcluster%y,comclustemp%y
!end if
31  if(clusmove .eqv. .false.) continue
    !calculate the energy penalty and insure no overlaps
!write(6,*) 'rorate done'
  end subroutine clusterrotation


  subroutine updateclusbonds(clsize,cllist,clno)
integer :: m,rb,pr,l,st,maxl,maxback,bdir
  integer,intent(inout)::clsize
    integer,dimension(:),intent(inout) ::cllist,clno
type(basicp),dimension(:),allocatable:: tempc
    logical :: adjver

allocate(tempc(maxlength))

 do rb = 1,clsize
        m = cllist(rb)
        maxl = chlen(m)

        do l = 1,maxl
           tempc(l)%x = protcoords(m,l)%x
           tempc(l)%y = protcoords(m,l)%y
           tempc(l)%z = protcoords(m,l)%z
           end do
        
        do pr = 1,m-1,1
           if(clno(pr) /= clno(m)) then
              do l = 1,maxl
                 maxback = chlen(pr)
                 do st = 1,maxback,1
             call adjacent(pr,l,st,tempc,adjver,bdir)
             if((adjver .eqv. .true.)) then
                br(pr,m,st,l) = 1
                if((bdir == -1*bonddd(m,l)) .and. (bdir ==bonddd(pr,st))) then
                   intbond(pr,m,st,l) = 1
                else
                   intbond(pr,m,st,l) = 0
                end if
             else if(adjver .eqv. .false.) then
                intbond(pr,m,st,l) = 0
                br(pr,m,st,l) = 0
             end if
          end do
       end do
       end if
    end do


    do pr = m+1,nprotein,1
           if(clno(pr) /= clno(m)) then
              do l = 1,maxl
                 maxback = chlen(pr)
                 do st = 1,maxback,1
             call adjacent(pr,l,st,tempc,adjver,bdir)
             if((adjver .eqv. .true.)) then
                br(m,pr,l,st) = 1
                if((bdir == -1*bonddd(m,l)) .and. (bdir ==bonddd(pr,st))) then
                   intbond(m,pr,l,st) = 1
                else
                   intbond(m,pr,l,st) = 0
                end if
             else if(adjver .eqv. .false.) then
                intbond(m,pr,l,st) = 0
                br(m,pr,l,st) = 0
             end if

          end do
       end do
       end if
    end do 


        end do

  end subroutine updateclusbonds
  
    subroutine tempcomfind(m,tcom,tempcluscoord)
    integer :: dcomx,dcomy,dcomz,maxl,comx,comy,comz
    double precision :: dcx,dcy,dcz
    integer :: l
    integer,intent(inout)::m
    type(basicp),dimension(:,:),intent(inout):: tempcluscoord
    type(centremass),dimension(:),intent(inout)::tcom

    comx = 0
    comy = 0
    comz = 0
        maxl = chlen(m)



    do l = 2,maxl,1
       dcomx = tempcluscoord(m,l)%x-tempcluscoord(m,l-1)%x
       if (dcomx == gridsize-1) dcomx = -1
       if (dcomx == 1-gridsize) dcomx = 1
       comx = comx + ((maxl+1-l)*dcomx)
       dcomy = tempcluscoord(m,l)%y-tempcluscoord(m,l-1)%y
       if (dcomy == gridsize-1) dcomy = -1
       if (dcomy == 1-gridsize) dcomy = 1
       comy = comy + ((maxl+1-l)*dcomy)
       dcomz = tempcluscoord(m,l)%z-tempcluscoord(m,l-1)%z
       if (dcomz == gridsize-1) dcomz = -1
       if (dcomz == 1-gridsize) dcomz = 1
       comz = comz + ((maxl+1-l)*dcomz)
    end do
    tcom(m)%x = modulo(tempcluscoord(m,1)%x +(real(comx)/(maxl)) -1,real(gridsize))+1
    tcom(m)%y = modulo(tempcluscoord(m,1)%y +(real(comy)/(maxl)) -1,real(gridsize))+1
    tcom(m)%z = modulo(tempcluscoord(m,1)%z +(real(comz)/(maxl)) -1,real(gridsize))+1


  end subroutine tempcomfind

 subroutine shorttempclustercom(m,clsize,tcom,comclustemp,cb,cllist)
    integer,intent(in)::m,clsize
    type(centremass),dimension(:),intent(inout) :: tcom
    integer,dimension(:),intent(in) :: cllist
        integer,dimension(:,:),intent(inout) ::cb
    integer ::  a,b,base,maxl,totcluspop,maxb,ra
    double precision :: delx,dely,delz,comx,comy,comz
    logical :: compass
    type(basicp),intent(inout) :: comclustemp
    type(rprot) :: dummycom
    t = time+1
    !write(6,*) 'com start'
    compass = .false.
    comclustemp%x = 0
    comclustemp%y = 0
    comclustemp%z = 0

        maxl = chlen(m)
    totcluspop = maxl 

    comx = 0.0d0
    comy = 0.0d0
    comz = 0.0d0

    do ra = 1,clsize
    a = cllist(ra)
          if(a/=m) then
        maxb = chlen(a)
             delx = tcom(a)%x - tcom(m)%x
             if (delx > gridsize/2) delx = delx-gridsize
             if (delx < -gridsize/2) delx = delx + gridsize
             comx = comx + (delx*maxb)
             dely = tcom(a)%y -tcom(m)%y
             if (dely > gridsize/2) dely = dely-gridsize
             if (dely < -gridsize/2) dely = dely + gridsize
             comy = comy + (dely*maxb)
             delz = tcom(a)%z-tcom(m)%z
             if (delz > gridsize/2) delz = delz-gridsize
             if (delz < -gridsize/2) delz = delz + gridsize
             comz = comz + (delz*maxb)
             totcluspop = totcluspop + maxb

          end if
    end do
    
    
    !write(6,*) 'stat'
    dummycom%x = modulo(tcom(m)%x + (comx/totcluspop)-1,real(gridsize))+1
    comclustemp%x = NINT(dummycom%x)
    !write(6,*) 'x coords',(com(1,m)%x + (comx/totcluspop)),comcluster%x
    dummycom%y = modulo(tcom(m)%y + (comy/totcluspop)-1,real(gridsize))+1
    comclustemp%y = NINT(dummycom%y)
    !write(6,*) 'y coords',(com(1,m)%y + (comx/totcluspop)),comcluster%y
    dummycom%z = modulo(tcom(m)%z + (comz/totcluspop)-1,real(gridsize))+1
    comclustemp%z = NINT(dummycom%z)

    
53  if(compass .eqv. .true.) continue

  end subroutine shorttempclustercom
  
  subroutine tempclustercom(m,clsize,tcom,comclustemp,cb,cllist,route,rcounter,totcluspop,checklist)
    integer,intent(in)::m
    integer,intent(inout)::clsize,totcluspop
    type(centremass),dimension(:),intent(inout) :: tcom
    integer,dimension(:),intent(in) :: cllist
    integer,dimension(:),intent(inout) :: checklist
    integer,dimension(:),intent(inout) ::rcounter
    integer,dimension(:,:),intent(inout) ::cb,route
    integer ::  a,b,base,maxl,maxb,ra,centralchain,f,g,cench
    double precision :: delx,dely,delz
    double precision,dimension(:,:),allocatable:: comx,comy,comz
    logical :: compass,dummyl
    type(basicp),intent(inout) :: comclustemp
    type(rprot) :: dummycom

    allocate(comx(nprotein,nprotein))
    allocate(comy(nprotein,nprotein))
    allocate(comz(nprotein,nprotein))
    !allocate(checklist(clsize))
  
    t = time+1
    !write(6,*) 'com start'
    compass = .false.
    comclustemp%x = 0
    comclustemp%y = 0
    comclustemp%z = 0

    centralchain = m
    dummycom%x = 0
    dummycom%y = 0
    dummycom%z = 0
dummyl = .false.
    
    if(clsize ==1) goto 76

   
    do a = 1,clsize-1,1
       b = cllist(a)
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
                   maxl = chlen(f)
                   delx = tcom(f)%x - tcom(b)%x
                   if (delx > (maxb+maxl)/2) delx = delx-gridsize
                   if (delx < -(maxb+maxl)/2) delx = delx + gridsize
                   comx(b,f) =  delx
                   comx(f,b) = -delx
                   dely = tcom(f)%y - tcom(b)%y
                   if (dely > (maxb+maxl)/2) dely = dely-gridsize
                   if (dely < -(maxb+maxl)/2) dely = dely + gridsize
                   comy(b,f) =  dely
                   comy(f,b) = -dely
                   delz = tcom(f)%z - tcom(b)%z
                   if (delz > (maxb+maxl)/2) delz = delz-gridsize
                   if (delz < -(maxb+maxl)/2) delz = delz + gridsize
                   comz(b,f) =  delz
                   comz(f,b) = -delz
                end if
             end do
             end do


   call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,dummyl)

    
76  if(clsize == 1) then


    end if


    comclustemp%x = NINT(modulo(tcom(centralchain)%x + (dummycom%x/totcluspop)-1,real(gridsize))+1)
    comclustemp%y = NINT(modulo(tcom(centralchain)%y + (dummycom%y/totcluspop)-1,real(gridsize))+1)
    comclustemp%z = NINT(modulo(tcom(centralchain)%z + (dummycom%z/totcluspop)-1,real(gridsize))+1)

    
53  if(compass .eqv. .true.) continue

  end subroutine tempclustercom

    subroutine rotcorrect(m,clno,clsize,dx,dy,dz,tempcluscoord,tcom,cllist,clusrotpass)
    integer,intent(in)::m,dx,dy,dz
    integer,intent(in):: clsize
    type(centremass),dimension(:),intent(inout)::tcom
    logical,intent(inout)::clusrotpass
    integer,dimension(:),intent(inout) ::clno,cllist
    integer :: vv,rvv
    type(basicp),dimension(:,:),intent(inout):: tempcluscoord
    logical :: moveallow
  
    moveallow = .true.
    do rvv = 1,clsize
       vv = cllist(rvv)
          call tempclustermove(vv,clno,dx,dy,dz,tempcluscoord,moveallow)
          if(moveallow .eqv. .false.) then
             clusrotpass = .false.
             return
          end if
          tcom(vv)%x = modulo(tcom(vv)%x + dx -1,real(gridsize))+1
          tcom(vv)%y = modulo(tcom(vv)%y + dy -1,real(gridsize))+1
          tcom(vv)%z = modulo(tcom(vv)%z + dz -1,real(gridsize))+1
    end do
  end subroutine rotcorrect

  subroutine cribond(m,ctbdi,axis,rotor,maxl)
    integer,intent(in):: m,axis,rotor,maxl
    integer,dimension(:,:),intent(inout)::ctbdi
    logical :: bm
    integer :: sign,s


    do s = 1,maxl-1

       !write(6,*) 'fd'
       bm = .true.
       if(abs(bdi(m,s)) == axis) then
          bm = .true.
          ctbdi(m,s) = bdi(m,s)
          goto 61
       end if
       !write(6,*) 'ab'

       if(rotor == 0) then
          bm = .true.
          ctbdi(m,s) = -1*bdi(m,s)
          !if(abs(bdi(m,s)) /= (modulo(axis,3) + 1)) ctbdi(m,s) = 1*bdi(m,s)
          goto 61
       end if
       !write(6,*) bdi(m,s),m,s 
       sign = bdi(m,s)/abs(bdi(m,s))

       if(axis == 1) then
          if(rotor == 1) then
             if(abs(bdi(m,s)) == 2) ctbdi(m,s) = sign*(abs(bdi(m,s)) + 1)
             if(abs(bdi(m,s)) == 3) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) - 1)
             goto 61
          else if(rotor == -1) then
             if(abs(bdi(m,s)) == 2) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) + 1)
             if(abs(bdi(m,s)) == 3) ctbdi(m,s) = sign*(abs(bdi(m,s)) - 1)
             goto 61
          end if
       else if(axis ==2 ) then
          if(rotor == 1) then
             if(abs(bdi(m,s)) == 1) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) + 2)
             if(abs(bdi(m,s)) == 3) ctbdi(m,s) = sign*(abs(bdi(m,s)) - 2)        
             goto 61
          else if(rotor == -1) then
             if(abs(bdi(m,s)) == 1) ctbdi(m,s) = sign*(abs(bdi(m,s)) + 2)
             if(abs(bdi(m,s)) == 3) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) - 2)
             goto 61
          end if

       else if(axis ==3) then
          if(rotor == 1) then
             if(abs(bdi(m,s)) == 1) ctbdi(m,s) = sign*(abs(bdi(m,s)) + 1)
             if(abs(bdi(m,s)) == 2) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) - 1)        
             goto 61
          else if(rotor == -1) then
             if(abs(bdi(m,s)) == 1) ctbdi(m,s) = -1*sign*(abs(bdi(m,s)) + 1)
             if(abs(bdi(m,s)) == 2) ctbdi(m,s) = sign*(abs(bdi(m,s)) - 1)         
          end if
          goto 61
       end if

61     if(bm .eqv. .true.) continue
    end do

    !if(ctbdi(m,s) == 0) write(6,*) 'bdi faillll'
  end subroutine cribond

  subroutine clusterflip(protno,cx,cy,cz,deltax,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
    integer,intent(in)::protno,cx,cy,cz,maxl
    integer,dimension(:,:),intent(inout)::ctbo,ctbdi
    !integer,dimension(:),intent(in)::clno
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay,deltaz
    integer :: l
    type(basicp),intent(inout) :: comcluster

    do l = 1,maxl
       tempcluscoord(protno,l)%x = modulo(comcluster%x + (cx*deltax(protno,l)) -1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(comcluster%y + (cy*deltay(protno,l)) -1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z + (cz*deltaz(protno,l)) -1,gridsize)+1
    end do

    if(cx == 1) then
       call clusrbond(protno,ctbo,1,0,maxl)
       call cribond(protno,ctbdi,1,0,maxl)
    else if(cy == 1) then
       call clusrbond(protno,ctbo,2,0,maxl)
       call cribond(protno,ctbdi,2,0,maxl)
    else if(cz == 1) then
       call clusrbond(protno,ctbo,3,0,maxl)
       call cribond(protno,ctbdi,3,0,maxl)
    end if

  end subroutine clusterflip

  subroutine rotatex(protno,cy,cz,deltay,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
    integer,intent(in)::protno,cy,cz,maxl
    integer,dimension(:,:),intent(inout)::ctbo,ctbdi
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltay,deltaz
    integer :: l
 type(basicp),intent(inout) :: comcluster
    do l = 1,maxl
       tempcluscoord(protno,l)%x = protcoords(protno,l)%x
       tempcluscoord(protno,l)%y = modulo(comcluster%y +((cy*deltaz(protno,l)) -1),gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z +((cz*deltay(protno,l)) -1),gridsize)+1
    end do
    !flipped rotation
    if(cy ==-1) then
       call clusrbond(protno,ctbo,1,1,maxl)
       call cribond(protno,ctbdi,1,1,maxl)
    else if(cz == -1) then
       call clusrbond(protno,ctbo,1,-1,maxl)
       call cribond(protno,ctbdi,1,-1,maxl)
    end if
  end subroutine rotatex

  subroutine rotatey(protno,cx,cz,deltax,deltaz,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
    integer,intent(in)::protno,cx,cz,maxl
    integer,dimension(:,:),intent(inout)::ctbo,ctbdi
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltaz
    integer :: l
     type(basicp),intent(inout) :: comcluster
    !write(6,*) 'aaaaaa'
    do l = 1,maxl
       tempcluscoord(protno,l)%y = protcoords(protno,l)%y
       tempcluscoord(protno,l)%x = modulo(comcluster%x +(cx*deltaz(protno,l))-1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z +(cz*deltax(protno,l))-1,gridsize)+1
    end do

    if(cx ==-1) then
       call clusrbond(protno,ctbo,2,-1,maxl)
       call cribond(protno,ctbdi,2,-1,maxl)
    else if(cz == -1) then
       call clusrbond(protno,ctbo,2,1,maxl)
       call cribond(protno,ctbdi,2,1,maxl)
    end if


  end subroutine rotatey

  subroutine rotatez(protno,cx,cy,deltax,deltay,tempcluscoord,ctbo,maxl,ctbdi,comcluster)
    integer,intent(in)::protno,cx,cy,maxl
    integer,dimension(:,:),intent(inout)::ctbo,ctbdi
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay
    integer :: l
 type(basicp),intent(inout) :: comcluster
 do l = 1,maxl
       tempcluscoord(protno,l)%z = protcoords(protno,l)%z
       tempcluscoord(protno,l)%x = modulo(comcluster%x +(cx*deltay(protno,l))-1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(comcluster%y +(cy*deltax(protno,l))-1,gridsize)+1
      ! if(time == 1163 .and. l ==1) write(6,*) cx,deltay(protno,l),comcluster%x,tempcluscoord(protno,l)%x
       !if(time == 1163 .and. l ==1) write(6,*) cy,deltax(protno,l),comcluster%y,tempcluscoord(protno,l)%y
    end do

    if(cx ==-1) then
       call clusrbond(protno,ctbo,3,1,maxl)
       call cribond(protno,ctbdi,3,1,maxl)
    else if(cy == -1) then
       call clusrbond(protno,ctbo,3,-1,maxl)
       call cribond(protno,ctbdi,3,-1,maxl)
    end if
  end subroutine rotatez
  
 subroutine tempclustermove(chainno,clno,delx,dely,delz,tempcluscoord,moveallow)
    integer,intent(in) :: chainno,delx,dely,delz
    logical,intent(inout)::moveallow
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::clno
    integer :: b,pr,st,maxback,maxl
    logical :: clusmove

    clusmove = .true.
    moveallow = .true.
    maxl = chlen(chainno)

    do b =1,maxl
       tempcluscoord(chainno,b)%x = modulo(tempcluscoord(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(tempcluscoord(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(tempcluscoord(chainno,b)%z + delz-1,gridsize)+1

       do pr = 1,nprotein
          if(clno(pr) == clno(chainno)) goto 39
maxback = chlen(pr)
          do st = 1,maxback,1
             if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
                clusmove = .false.
                !write(6,*) 'translation fail'
                goto 31
             end if
          end do
39        if(clusmove .eqv. .true.) continue
       end do

    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine tempclustermove

  subroutine clustermove(clno,chainno,delx,dely,delz,tempcluscoord,moveallow,ctbo,ctbdi)
    integer,intent(in) :: chainno,delx,dely,delz
    integer,dimension(:,:),intent(inout)::ctbo,ctbdi
    logical,intent(inout)::moveallow
    Type(basicp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::clno
    integer :: b,pr,st,maxback,maxl
    logical :: clusmove

    clusmove = .true.
    moveallow = .true.
        maxl = chlen(chainno)
    do b =1,maxl
       tempcluscoord(chainno,b)%x = modulo(protcoords(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(protcoords(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(protcoords(chainno,b)%z + delz-1,gridsize)+1
       ctbo(chainno,b) = bonddd(chainno,b)
       if(b<maxl) ctbdi(chainno,b) = bdi(chainno,b)

       do pr = 1,nprotein
          if(clno(pr) == clno(chainno)) goto 39
        maxback = chlen(pr)
          do st = 1,maxback,1
             if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
                clusmove = .false.
                !write(6,*) 'translation fail'
                goto 31
             end if
          end do
39        if(clusmove .eqv. .true.) continue
       end do

    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove

  subroutine debug
    integer::m,l,f,g,dx,dy,dz,maxlengthss,maxl,zz
    integer,dimension(:,:),allocatable:: dbbdi
    allocate(dbbdi(nprotein,maxlength))

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

             !write(6,*) 'l',l,'f',f

             if((abs(protcoords(m,l)%x - protcoords(m,l-1)%x) == 1) .or. (abs(protcoords(m,l)%x - &
                  protcoords(m,l-1)%x) == (gridsize - 1))) then
                dx = 1
             else if(abs(protcoords(m,l)%x - protcoords(m,l-1)%x) == 0) then
                dx = 0
             end if
             !energypass = .false.
             if((abs(protcoords(m,l)%y - protcoords(m,l-1)%y) == 1) .or. (abs(protcoords(m,l)%y - &
                  protcoords(m,l-1)%y) == (gridsize - 1))) then
                dy = 1
             else if(abs(protcoords(m,l)%y - protcoords(m,l-1)%y) == 0) then
                dy = 0
             end if
             !energypass = .false.
             if((abs(protcoords(m,l)%z - protcoords(m,l-1)%z) == 1) .or. (abs(protcoords(m,l)%z &
                  - protcoords(m,l-1)%z) == (gridsize -1))) then
                dz = 1
             else if(abs(protcoords(m,l)%z - protcoords(m,l-1)%z) == 0) then
                dz = 0
             end if
             sumdebug = dx + dy + dz
             if(sumdebug /= 1.0) then
                write(6,*) 'SPLIT'
                finalfail = .true.
             end if
          end if

          if(bonddd(m,l) == 0) write(6,*) 'bond vector fail'
       end do
    end do

    do m=1,nprotein
        maxl = chlen(m)
       l = 1
       if((bonddd(m,l) == bdi(m,l))) write(6,*) 'bond on bond 111'
       do l=2,maxl-1,1
          if((bonddd(m,l) == bdi(m,l)) .or. (bonddd(m,l) == -1*bdi(m,l-1))) &
               write(6,*) 'bond on bond',m,l,bonddd(m,l),bdi(m,l),-bdi(m,l-1),time
       end do
       l = maxl
       if((bonddd(m,l) == -1*bdi(m,l-1))) write(6,*) 'bond on bond 222',m,l
    end do

    do m = 1,nprotein
        maxl = chlen(m)
       do l= 1,maxl-1

          if((protcoords(m,l+1)%x - protcoords(m,l)%x == 1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
               -1*(gridsize-1))) then
             dbbdi(m,l) = 1 
          else if((protcoords(m,l+1)%x - protcoords(m,l)%x == -1) .or. (protcoords(m,l+1)%x - protcoords(m,l)%x ==&
               (gridsize-1))) then
             dbbdi(m,l) = -1
          else if((protcoords(m,l+1)%y - protcoords(m,l)%y == 1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
               -1*(gridsize-1)))then
             dbbdi(m,l) = 2
          else if((protcoords(m,l+1)%y - protcoords(m,l)%y == -1) .or. (protcoords(m,l+1)%y - protcoords(m,l)%y ==&
               (gridsize-1))) then
             dbbdi(m,l) = -2
          else if((protcoords(m,l+1)%z - protcoords(m,l)%z == 1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
               -1*(gridsize-1))) then
             dbbdi(m,l) = 3
          else if((protcoords(m,l+1)%z - protcoords(m,l)%z == -1) .or. (protcoords(m,l+1)%z - protcoords(m,l)%z ==&
               (gridsize-1))) then
             dbbdi(m,l) = -3
          else
             write(6,*) 'no assignment'
          end if
          if(dbbdi(m,l) /= bdi(m,l)) write(6,*) 'backbone chain fail',m,l,bdi(m,l),dbbdi(m,l)

       end do
       !do zz = 3,5
       !write(6,*) protcoords(m,zz)%x,protcoords(m,zz)%y,protcoords(m,zz)%z
       !end do
    end do
    !write(6,*) 'dubug success'

  end subroutine debug

  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(basicp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    !allocate(mashist(10))
    normalisedsize = 0
    clcount = 0
    clusterpop = 0
    clustcount = 0
    do m = 1,nprotein
       clnos(m) = m
       cllist(m) = 0
       histcl(m) = 0
    end do
    a = 1



    do m = 1,nprotein-1,1
        maxback = chlen(m)
       do g = m+1,nprotein,1
          if((g/=m) .and. (clnos(m) /= clnos(g))) then
                     maxlengthss = chlen(g)
             do l = 1,maxback
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z
                do f = 1,maxlengthss
                   call adjacent(g,l,f,tempcoord,adjver,bdir)              
                   if(adjver.eqv. .true.) then
                      if(isbond .eqv. .true.) then
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
                      if((bdir == bonddd(g,f)) .and. (bdir == (-1*bonddd(m,l))) .and. &
                           (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  /= 0.0)) then

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
          end if
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


    maxclus = maxval(histcl)
!write(6,*) 'maxcl',maxclus
    !do g = 1,10,1
      !do m = 1,nprotein
          !if((histcl(m) > (((g-1)/10.0)*(maxclus))) .and. (histcl(m) <= ((g/10.0)*maxclus))) then
             !mashist(g) = mashist(g) + histcl(m)**2
             !end if
          !end do
          !write(38,*) time, (g/10.0)*(maxclus),mashist(g)
    !end do

  
    do m = 1,nprotein
write(38,*) time,histcl(m),histcl(m)**2
       end do
    
!if(noinfo .eqv. .false.) then
    !if(clustcount /= 0) write(6,*) time,clustcount, 'cluster!'
    if(clustcount /= 0) write(82,*) time,clustcount, clusterpop, real(clusterpop)/clustcount, real((normalisedsize) &
         +(nprotein-clusterpop))/nprotein,maxclus !,clcount
    !if(time>10000) write(6,*) 'output',time, real((normalisedsize)+(nprotein-clusterpop))/nprotein,real((normalisedsize)+ &
         !(nprotein-clusterpop))
    if(clustcount == 0)  write(82,*) time,clustcount, clusterpop, 0.0
!end if
        deallocate(cllist)
  end subroutine clustercount


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
          write(67,*) 'atom', acount, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species), 'resid ',protcoords(m,1)%species
          write(67,*) 'atom', acount+1, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species)-1,'resid ',protcoords(m,1)%species
          acount = acount+2
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
          write(67,*) 'bond', adjustr(atnum)//':',acount+2
          write(67,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+2
       end do
       write(atnum,'(i7)') acount
       write(67,*) 'bond', adjustr(atnum)//':',acount+1
       acount = acount +2
    end do
    
  end subroutine pdbsetup


  

  SUBROUTINE read_setup

    USE keywords
    !USE setup_vars
    !USE stats_vars
    IMPLICIT NONE

    CHARACTER(LEN=20) :: keyword, option, argument
    CHARACTER(LEN=18), PARAMETER :: param_fmt1='(A16, 1PE20.10, A)'
    CHARACTER(LEN=16), PARAMETER :: param_fmt2='(A16, F16.10)'
    INTEGER :: err, i, j,runtype,dummy2,dummytype,xl,f
    LOGICAL :: success

 

    !  Defaults and descriptions.
    !  S.T.A.T. stands for short-time averaged temperature.




    seed = 6                      !seed for random number generation
    maxlength = 20              !number of beads in chain
    nprotein = 1                 !number of chains
    gridsize= 100        !dimensions of the lattice
    maxtime= 10000                !length of simulation
    equilib = 5000              !equilibration time
    right = 1                 !1 = rightangle move 0 = no rightangle move
    crank = 1                      !1 = crank move 0 = no crank move
    rept = 1                  !1 = reptation 0 = no reptations
    piv = 1                      !1 = pivot 0 = no pivots
    kT = 1.0
    film = .true.
    debugyes = .true.
    maxpiv = 10
    cranklimit = 5
    !interenergy(1,2) =  -5.0 

    !traj_file = 'trajectory'       ! Stem of name for trajectory files.
    !vgen = 1                       ! 1=random initial velocities (no momentum), 2=Maxwell-Boltzmann
    !v_init_file = ' '              ! Name of starting velocities file.
    !v_rescale = .FALSE.            ! If reading velocities from a file, .T.=rescale to get specified K.E.
    !wellstats_file = 'wellstats'   ! Stem of name for individual well statistics files.
    !well_tol = 1.0D-6              ! Energy tolerance for quenches to be considered the same.
    !xmol_type(1) = ''              ! Atom type for A (or charged) atoms in xmol dumps.
    !xmol_type(2) = ''              ! Atom type for B (or uncharged) atoms in xmol dumps.


    OPEN (UNIT=20, FILE='interactionsetup.txt', STATUS='OLD', IOSTAT=err)

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

       !CASE ('CHAIN_LENGTH')
          !CALL get_integer(maxlength1)

       !CASE('CHAIN_2_LENGTH')
          !call get_integer(maxlength2)

       !CASE ('NUM_CHAINS')
          !CALL get_integer(nprotein1)

       !CASE ('NUM_CHAINS_2')
          !CALL get_integer(nprotein2)

       CASE ('TOTTYPES')
         CALL get_integer(nspecies)
         allocate (speciespop(nspecies))
         allocate(specieslen(nspecies))
         
       CASE ('TOTCHAINS')
            !CALL get_integer(nprotein)
runtype = 1
       CASE ('NEWCHAIN')
          CALL get_integer(specieslen(runtype))
          CALL get_integer(speciespop(runtype))
          do f = 1,specieslen(runtype)
             call get_integer(dummytype)
             if(dummytype>mtype) mtype = dummytype
          end do
          write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
          runtype = runtype + 1
          
       CASE ('LATTICE_DIMENSIONS')
          CALL get_integer(gridsize)

         
 allocate(interen(mtype,mtype))
  allocate(interenergy(mtype,mtype))
  allocate(intraenergy(mtype))
  !allocate(tttt())
          
       CASE ('SIM_TIME')
          CALL get_integer(maxtime)
maxtime = maxtime*1000
       CASE ('EQUIB_TIME')
          CALL get_integer(equilib)
equilib =equilib*1000
       CASE ('RIGHTANGLE')
          CALL get_integer(right)

       CASE ('CRANK')
          CALL get_integer(crank)

       CASE ('REPTATION')
          CALL get_integer(rept)

       CASE ('PIVOT')
          CALL get_integer(piv)

       CASE ('KT')
          CALL get_dp(kT)

       CASE ('FILM')
          CALL get_logical(film)

       CASE ('DEBUGYES')
          CALL get_logical(debugyes)

          !CASE ('ISBOND')
          !CALL get_logical(isbond)

       CASE ('MAXPIV')
          CALL get_integer(maxpiv)

       CASE ('CRANKLIMIT')
          CALL get_integer(cranklimit)



       CASE ('INTRAEN')
          CALL get_dp(intraen)
          intraen = -ABS(intraen)
          

       CASE ('INTERACTION')
          call get_integer(xl)
          call get_dp(intraenergy(xl))
intraenergy(xl) = -abs(intraenergy(xl))
          do f = 1,mtype
             call get_dp(interenergy(xl,f))
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
          
       CASE ('SCALEINFO')
          CALL get_logical(scalinginfo)
         
CASE ('RESTART')
 
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



!          if((abs(comclustemp%x - comcluster%x)>1) .and. (abs(comclustemp%x - comcluster%x)<gridsize-1)) then
!          do rb = 1,clsize
!             b = cllist(rb)
!             write(6,*) 'com',com(1,b)%x,tcom(b)%x
!          end do
!          write(6,*) 'xfail', (comclustemp%x - comcluster%x),comcluster%x,comclustemp%x
!          end if
!          if((abs(comclustemp%y - comcluster%y)>1) .and. (abs(comclustemp%y - comcluster%y)<gridsize-1)) then
!          do rb = 1,clsize
!             b = cllist(rb)
!             write(6,*) 'com',com(1,b)%y,tcom(b)%y
!          end do
!             write(6,*) 'yfail',(comclustemp%y - comcluster%y),comcluster%y,comclustemp%y
!             end if
!             if((abs(comclustemp%z - comcluster%z)>1) .and. (abs(comclustemp%z - comcluster%z)<gridsize-1)) then
!                          do rb = 1,clsize
!             b = cllist(rb)
!             write(6,*) 'com',com(1,b)%z,tcom(b)%z
!          end do
!                write(6,*) 'zfail', (comclustemp%z - comcluster%z),comcluster%z,comclustemp%z

                !end if
