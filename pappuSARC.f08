program move

  implicit none

  double precision :: totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength,choose2
  double precision :: actualchainlength,random,totiso
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug,testdummy,polymerrog,polymerl
  integer :: piv,right,endm,crank,rept,datayes,maxclussize,cranklimit,minl,outputrate,backtime
  integer ::lastruntime,timebase,steplength
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count,avecluspop,aveclusnumber
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib,cltacc,cltatt,clracc,clratt
  integer :: successful,reject,pivotlimit,maxpiv,maxlength1,maxlength2,q,d,u,sm,qt,ft,gt,scans
  integer:: cracc,crrej,raacc,rarej,rerej,reacc,pvacc,pvrej,nprotein1,nprotein2,nspecies,mtype,endacc,endatt
  integer,dimension(:,:),allocatable :: bdi
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen
  double precision :: crrat,emrat,rarat,pvrat,rerat,trackx,tracky,trackz,intraen,intdummy,totdire
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy
  real, external :: ran2
  real,dimension(:),allocatable :: variance,disp
  logical :: exist,fail,finalfail,debugyes,film,isbond,clt,clr,noinfo,scalinginfo,nobonds,restart
  real::start,finish,midpoint
  integer::removethis,pivde
  logical :: suffclust
  !change COM code for large protein in small box to stop pbc effect - similar to old version - but still
  !use dispacement from central bead

  !cluster COM may fall foul of pbcs

  !Sort MSD outputs

  
 type protein
     integer :: x,y,z,species,type,bstate,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl,bondg,bondf
  end type protein

  type rprot
     real :: x,y,z
  end type rprot

  type basicp
     integer :: x,y,z
  end type basicp


  type prottemp
     integer :: x,y,z,bstate,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl,bondf,bondg
  end type prottemp

  type centremass
     double precision :: x,y,z
  end type centremass


  type(protein),dimension(:,:),allocatable :: protcoords
  type(centremass),dimension(:,:),allocatable :: com

  mtype = 0
  !noinfo = .false.
  call read_setup

suffclust = .false.
  nprotein = sum(speciespop)
  write(6,*) 'nprotein = ', nprotein
  maxlength = maxval(specieslen)
  write(6,*) 'maxlength =',maxlength
  !open(17, file = 'setup2.txt', action = 'read')
  if(restart .eqv. .false.) then
     open(23, file = 'initialtake2.xyz', action = 'read')
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',action = 'write')
        open(82,file = 'clusterdata.dat', action = 'write')
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
           open(82,file = 'clusterdata.dat', action = 'write')
           open(13,file='timedata.dat',action = 'write')
           open(19,file='acceptance.dat',action  = 'write')
           open(38,file = 'histclust.dat',action = 'write')
           !open(67,file='bondddd.dat',action  = 'write')
        end if
     end if
  else if(restart .eqv. .true.) then
     open(23, file = 'checkpoint.xyz', action = 'read')
     write(6,*) 'logicals',scalinginfo,noinfo
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',access = 'append')
        open(82,file = 'clusterdata.dat', access = 'append')
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
           open(82,file = 'clusterdata.dat', access = 'append')
           open(13,file='timedata.dat',access = 'append')
           open(77,file='phasetrans.dat',action = 'write')
           open(19,file='acceptance.dat',access  = 'append')
           open(38,file = 'histclust.dat',access = 'append')
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
  if(restart .eqv. .true.) equilib = 0
  maxtime = maxtime + equilib
  outputrate = (maxtime-equilib)/100
  clracc = 0
  clratt = 0
  cltacc = 0
  cltatt = 0
  !kT = 50.0
  successful = 0
  reject = 0
  backtime = maxtime


  totdire = 0.0
  totiso = 0.0d0
  do ft = 1,mtype
     do gt = 1,mtype
        totdire = totdire + interenergy(ft,gt)
        totiso = totiso + interen(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
  end do
  if(totiso == 0.0d0) isbond = .false.
  if(totdire == 0.0d0) nobonds =  .true.

  !interenergy = -1.0

  !write(6,*) maxlength1,maxlength2
  write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT
  write(6,*) 'clr',clr,'clt',clt

  if(restart .eqv. .false.) then
     open(67, file = 'move.vtf', action = 'write')
open(88,file = 'movetagged.vtf',action='write')
  else if(restart .eqv. .true.) then
     open(67, file = 'move.vtf', access = 'append')
     open(88, file = 'movetagged.vtf', access = 'append')
  end if


  allocate(protcoords(nprotein,maxlength))
  allocate(com(2,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(bdi(nprotein,maxlength-1))
  allocate(chlen(nprotein))

  N = maxlength


  count = 0
  time = 0
  totdisp = 0.0d0


  WRITE(6,*) 'ISBOND =',isbond
  write(6,*) 'interen',interen
  write(6,*) 'intraen',intraenergy
  cracc = 0
  crrej = 0
  rerej = 0
  reacc = 0
  raacc = 0
  rarej = 0
  pvacc = 0
  pvrej = 0


  !sets the limit of the pivot move


  call CPU_TIME(start)

  write(6,*) 'foundation start'
  if(restart .eqv. .false.) then
     call foundation
  else if (restart .eqv. .true.) then
     call foundationrestart
  end if
  write(6,*) 'foundation complete'

  if((scalinginfo .eqv. .false.) .and. (restart .eqv. .false.)) call pdbsetup
  minl = minval(specieslen)
  pivotlimit = min((minl/2),maxpiv)

  write(6,*) 'a'
  !call dataout
  do qt = 1,nprotein
     call comfind(qt,.false.)
  end do
  write(6,*) 'b'
  call clustercount

  actualchainlength = 0.0
  do timebase = 1,maxtime,1
     if(restart .eqv. .false.) then
        time = timebase
     else
        time = timebase+lastruntime
     end if
     fail = .false.

    ! if(time == 18440) STOP
     if ((mod(time,outputrate) == 0) .and. (film .eqv. .true.) ) then
        call dataout
     end if

     if((nobonds .eqv. .true.) .and. (clt .eqv. .false.) .and. (clr .eqv. .false.)) then
                 call positioning
        if((modulo(time,maxtime/100) == 0) .and. (time>equilib)) then
       write(29,*) (time-equilib), real(totrmsbrute)/nprotein
    end if
     else
        call pickmoves
     end if
     if (mod(time,outputrate) == 0) then
        if(noinfo .eqv. .false.)  write(93,*) time,totalenergy
     end if
     if(modulo(time,outputrate*10) == 0 .and. (scalinginfo .eqv. .false.)) call clustercount

     if((time > equilib) .and. (modulo(time,outputrate) == 0) .and. (scalinginfo .eqv. .true.)) then
        call length
        call radiusofgyration
        if((noinfo .eqv. .false.) .and. (scalinginfo .eqv. .true.)) write(79,*) time-equilib,&
             runningaveEtE/((time-equilib)/outputrate),&
             runningaveROG/((time-equilib)/outputrate)
     end if


     if (modulo(time,outputrate) == 0 .and. debugyes .eqv. .true.) then
        call debug
        call energy(.false.)
     end if

     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else if(modulo(time,100) == 0) then
        write(6,*) 'step =',time
     end if


     if((modulo(time,maxtime/1000) == 0) .and.(noinfo .eqv. .false.)) then
        write(19,*) 'endmove',endacc,endatt, real(endacc)/endatt
        write(19,*) 'reptation',reacc,rerej, real(reacc)/(reacc+rerej)
        write(19,*) 'right angle',raacc,rarej, real(raacc)/(raacc+rarej)
        write(19,*) 'cluster translation',cltacc,cltatt,real(cltacc)/cltatt
        
     end if
     call CPU_TIME(midpoint)
     if((midpoint-start)> 252000) then
        backtime = time
        goto 93
     end if
     !write(6,*) 'sweep complete'
  end do

93 if((film .eqv. .false.) .and. (scalinginfo .eqv. .false.)) call dataout
  write(6,*) 'deltax =',trackx, 'deltay =',tracky, 'deltaz =',trackz
  call CPU_TIME(finish)
  call dataoutrestart
  if(noinfo .eqv. .false.) write(13,*) (finish-start)
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
  !write(6,*) 'end'
contains


  subroutine phase(mz,clsize,clnos,cb)
    integer :: m,l,maxl,h,a,totcluspop,rogpop
    double precision :: radofgy
    type(centremass)::avepos,sumsq
    integer,intent(inout)::clsize,mz
    integer,dimension(:),allocatable::cllist
    integer,dimension(:,:),allocatable:: route
    integer,dimension(:),allocatable::rcounter,checklist
    integer,dimension(:,:),intent(inout)::cb
    integer,dimension(:),intent(inout)::clnos
    type(basicp) :: comcluster
    logical :: toolarge

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
    write(77,*) time,radofgy,sqrt(radofgy)


  end subroutine phase
  
subroutine updatebondold(m,beadmin,beadmax,tempcoord,bond)
    integer,intent(in):: m,beadmin,beadmax
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer,dimension(:,:),intent(inout):: bond
    integer ::l,g,f,bg,bf,pick,reshuffle,rez

    !write(6,*) 'updateeeee'

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

       if((protcoords(m,l)%bstate == 2) .or. (protcoords(m,l)%bstate == 3)) goto 18

       29 continue
       if(bond(l,1) ==0) goto 18
       
       pick = int(ran2(seed)*bond(l,1))+1
       if(bond(l,pick)==1) then
          g = protcoords(m,l)%am
          f = protcoords(m,l)%al
       else if(bond(l,pick)==2) then
          g = protcoords(m,l)%bm
          f = protcoords(m,l)%bl
       else if(bond(l,pick)==3) then
          g = protcoords(m,l)%cm
          f = protcoords(m,l)%cl
       else if(bond(l,pick)==4) then
          g = protcoords(m,l)%dm
          f = protcoords(m,l)%dl
       else if(bond(l,pick)==5) then
          g = protcoords(m,l)%em
          f = protcoords(m,l)%el
       else if(bond(l,pick)==6) then
          g = protcoords(m,l)%fm
          f = protcoords(m,l)%fl
       end if


          reshuffle = pick
          if(protcoords(g,f)%bstate == 1) then
             protcoords(m,l)%bstate = 3
             protcoords(g,f)%bstate = 3
             protcoords(m,l)%bondg = g
             protcoords(m,l)%bondf = f
             protcoords(g,f)%bondg = m
             protcoords(g,f)%bondf = l
          else
             do rez = reshuffle,bond(l,1)
                bond(l,rez) = bond(l,rez+1)
             end do
             bond(l,1) = bond(l,1) -1
             goto 29
          end if
    
          
       18 continue

    end do

  end subroutine updatebondold


  subroutine vectorspin(m,l,maxl)
    integer,intent(in)::maxl,m,l
    integer::p,f,dum1,dum2,no1,no2,maxt,g,bdir,maxback,pr,str
    integer,dimension(:),allocatable:: tbdi
    integer,dimension(7) :: bond
    type(prottemp),dimension(:),allocatable::tempcoord
    double precision :: deltaenergy
    logical :: adjver
    allocate(tbdi(maxl))
    allocate(tempcoord(maxl))

    !call energy(.false.)
    !write(6,*) 'pre spin'
    !return
    !roatt = roatt + 1   
    !write(6,*) maxl,maxlength1,maxlength2,protcoords(m,1)%species
    !do a = 1,nprotein*maxlength

    !do l = p,f
    !write(6,*) 'vector' ! totalenergy
    call bonddirection(tempcoord(l)%bstate)
    tempcoord(l)%x = protcoords(m,l)%x
    tempcoord(l)%y = protcoords(m,l)%y
    tempcoord(l)%z = protcoords(m,l)%z

    !if(tempcoord(l)%z == -3) write(6,*) 'failll',time
    deltaenergy = 0.0d0

    tempcoord(l)%am = 0
    tempcoord(l)%bm = 0
    tempcoord(l)%cm = 0
    tempcoord(l)%dm = 0
    tempcoord(l)%em = 0
    tempcoord(l)%fm = 0
    tempcoord(l)%bondg = 0
    tempcoord(l)%bondf = 0
    bond(1) = 0
if(l<maxl) tbdi(l) = bdi(m,l)
    call removeenergyintra(m,l,deltaenergy,tempcoord)
    call removeenergyinter(m,l,deltaenergy,tempcoord)


    do str = 1,l-3,1
       call adjacent(m,l,str,tempcoord,adjver,bdir)
       if((adjver.eqv. .true.)) then
          call bondformintra(m,l,str,tempcoord,bdir,deltaenergy,bond)
       end if
    end do

    do str = l+3,maxl,1
       call adjacent(m,l,str,tempcoord,adjver,bdir)
       if((adjver.eqv. .true.)) then
          call bondformintra(m,l,str,tempcoord,bdir,deltaenergy,bond)
       end if
    end do

    do pr = 1,nprotein
       if(pr /= m) then
          maxback = chlen(pr)
          do str = 1,maxback
             call adjacent(pr,l,str,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy,bond)
             end if
          end do
       end if
    end do

    if(bond(1)>0) then
       call bondchoose(bond,tempcoord,deltaenergy,m,l)
    end if
    if(time == 18432) then
       write(6,*) 'data',m,l
       write(6,*) 'tempcoord',tempcoord(l)
       write(6,*) 'prtocoords',protcoords(m,l)
       write(6,*) 'prtocoords of 74',protcoords(74,7)
       end if
 !write(6,*) 'delta',deltaenergy
    
    if(Energydecision(deltaenergy) .eqv. .true.) then
       !write(6,*) 'update',bond(1),bond(2)
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l,l,tempcoord,tbdi,maxl)
       call updatebond(m,l,l,tempcoord)
       !roacc = roacc + 1
    end if
    deallocate(tempcoord)

    !write(6,*) 'vector end'
        !call energy(.false.)
    write(6,*) 'vector spin'
 call debug
  end subroutine vectorspin

  

  subroutine bondchoose(bond,tempcoord,deltaenergy,m,l)
    integer,dimension(:),intent(inout)::bond
    type(prottemp),dimension(:),intent(inout)::tempcoord
    double precision,intent(inout)::deltaenergy
    integer,intent(in):: m,l
    integer:: pick,g,f,reshuffle,rez

     if((tempcoord(l)%bstate == 2)) goto 18

       29 continue
       if(bond(1) ==0) goto 18
       
       pick = int(ran2(seed)*bond(1))+2
       if(bond(pick)==1) then
          g = tempcoord(l)%am
          f = tempcoord(l)%al
       else if(bond(pick)==2) then
          g = tempcoord(l)%bm
          f = tempcoord(l)%bl
       else if(bond(pick)==3) then
          g = tempcoord(l)%cm
          f = tempcoord(l)%cl
       else if(bond(pick)==4) then
          g = tempcoord(l)%dm
          f = tempcoord(l)%dl
       else if(bond(pick)==5) then
          g = tempcoord(l)%em
          f = tempcoord(l)%el
       else if(bond(pick)==6) then
          g = tempcoord(l)%fm
          f = tempcoord(l)%fl
       end if

!write(6,*) 'info',g,f,bond(1),bond(2),protcoords(m,l)%cm
          reshuffle = pick
          if(protcoords(g,f)%bstate == 1) then
             tempcoord(l)%bstate = 3
             tempcoord(l)%bondg = g
             tempcoord(l)%bondf = f
             if (g==m) then
                deltaenergy = deltaenergy + intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
             else if(g/=m) then
                deltaenergy = deltaenergy + interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
             end if
          else
             do rez = reshuffle,bond(1)
                bond(rez) = bond(rez+1)
             end do
             bond(1) = bond(1) -1
             goto 29
          end if
    

       18 continue


  end subroutine bondchoose
  

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


       call bonddirection(protcoords(m,l)%bstate)
 
       maxl = chlen(m)
       do l =2,maxl-1

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

  call bonddirection(protcoords(m,l)%bstate)
       end do
       l = maxl
  call bonddirection(protcoords(m,l)%bstate)

       !if(bonddd(m,l) == 0) write(6,*) 'fail $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
    end do

    do m =1,nprotein
       maxl = chlen(m)
       do l =1,maxl-1
          if((bdi(m,l) == 0) .or. (abs(bdi(m,l)) >3)) write(6,*) 'fail',m,l
       end do
    end do


    call energy(.true.)
    call debug

  end subroutine foundation

  subroutine foundationrestart
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z
    character(len = 10) :: BIN
    read(23,*) lastruntime,totrmsbrute,suffclust
    do m = 1,nprotein,1
       read(23,*) protcoords(m,1)%species, protcoords(m,1)%x, protcoords(m,1)%y, protcoords(m,1)%z,&
            protcoords(m,1)%type,protcoords(m,1)%bstate,protcoords(m,1)%bondg,protcoords(m,1)%bondf,chlen(m)

       maxl = chlen(m)

       do l = 2,maxl,1

          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z,&
               protcoords(m,l)%type,protcoords(m,l)%bstate,protcoords(m,l)%bondg,protcoords(m,l)%bondf
       end do
    end do


    do m = 1,nprotein
       l = 1

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


       maxl = chlen(m)
       do l =2,maxl-1

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

       end do

    end do

    do m =1,nprotein
       maxl = chlen(m)
       do l =1,maxl-1
          if((bdi(m,l) == 0) .or. (abs(bdi(m,l)) >3)) write(6,*) 'fail',m,l
       end do
    end do


    call energy(.true.)
    call debug

    close(23,status = 'delete')
  end subroutine foundationrestart


  subroutine clusterposition
    integer::scans,m,clsize,clback,h,a !,clsizetest
    real :: decider
    integer,dimension(:),allocatable:: cllist,cllistfind
    integer,dimension(:,:),allocatable::cb
    logical,dimension(:),allocatable:: cl
!integer,dimension(:),allocatable:: clno
    allocate(cllistfind(nprotein))
    allocate(cb(nprotein,nprotein))
    allocate(cl(nprotein))
    !allocate(clno(nprotein))

    cl = .false.
    
    m =int(ran2(seed)*(nprotein))+1
  

    if((clr .eqv. .false.) .and. (clt .eqv. .false.)) then
       call positioning
       return
    end if
!write(6,*) 'preassign'
    call clusterassignmove(m,cllistfind,clsize,cl)


    !call clusterassign(m,clno,clsizetest,cb)
    !if(clsize /= clsizetest) write(6,*) 'allocate fail',clsize,clsizetest
    allocate(cllist(clsize))
!write(6,*) 'assigned'
    if((1.0/clsize) > (0.1*ran2(seed))) then
    do a = 1,clsize
       cllist(a) = cllistfind(a)
    end do



   
          call clustertranslation(m,clsize,cb,cllist,cl)
 
    end if
    !write(6,*) 'cluster complete'
  end subroutine clusterposition

  subroutine pickmoves
    integer:: scans,m,maxl,decide,l,a

    scans = 1


    decide = int(ran2(seed)*10)+1

    if(decide <= 4) then
       do a=1,maxlength
          m =int(ran2(seed)*(nprotein))+1
          maxl = chlen(m)
          l = int(ran2(seed)*(maxl))+1
          call vectorspin(m,l,maxl)
       end do

    end if
    if(decide > 4 .and. decide <= 8) call positioning
    if(decide > 8) call clusterposition
    !if(time > 998 .and. time < 1002) write(6,*) 'information',decide,time
    !if(time == 1166) write(6,*) 'hello'
    !end do

    if((modulo(time,maxtime/100) == 0) .and. (time>equilib)) then
       write(29,*) (time-equilib-1)*scans, real(totrmsbrute)/nprotein
    end if

  end subroutine pickmoves


  subroutine rotate(m,l,b,tempcoords,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz)
    integer,intent(in):: m,l,b,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz
    type(prottemp),dimension(:),intent(inout) :: tempcoords
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
    integer ::randomz,pivlen,end
    !if(time > 42000 .and. (time < 47000)) write(6,*) 'start positioning'
    t = time + 1

    count = count + 1
    run = .true.
    run2 = .true.
    !return
    crank = 0
    piv = 0
end = 1
    nmoves = piv + end + right + rept
    randomz = int(ran2(seed)*nmoves)+1

 
    if (randomz == (right) .and. right ==1) then

       call rightanglemove
    else if (randomz == (right + rept) .and. rept ==1) then

       call reptation
       else if (randomz == (right + rept + end) .and. end ==1) then

       call endmove
           end if


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




  subroutine ribond(m,s,tbdi,axis,rotor)
    integer,intent(in):: m,s,axis,rotor
    integer,dimension(:),intent(inout)::tbdi
    logical :: bm
    integer :: sign


    bm = .true.
    if(abs(bdi(m,s)) == axis) then
       bm = .true.
       tbdi(s) = bdi(m,s)
       return
    end if


    if(rotor == 0) then
       bm = .true.
       tbdi(s) = -1*bdi(m,s)
       return
    end if

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
    integer :: m,l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz,bdir,str
    double precision :: deltaenergy
    logical :: rac,racc,adjver
    integer:: pr,maxl,maxback,st
    integer,dimension(7)::bond
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbdi
    allocate(tempcoord(maxlength))
    allocate(tbdi(maxlength))


    m =int(ran2(seed)*(nprotein))+1
    maxl = chlen(m)
    ! l cannot be equal to 1 or maxlength
    l = int(ran2(seed)*(maxl-2))+2
    !if(time == 2370) write(6,*) 'right', m,l

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


       call rightvector(m,l,dx,dy,dz,racc,tbdi,maxl)

       if(racc .eqv. .false.) then
          rac = .false.
          goto 37
       end if
       !write(6,*) tbo(l-1),tbo(l),tbo(l+1)
       tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
       tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
       tempcoord(l)%z = modulo(protcoords(m,l)%z+dz-1,gridsize)+1

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



          tempcoord(l)%am = 0
          tempcoord(l)%bm = 0
          tempcoord(l)%cm = 0
          tempcoord(l)%dm = 0
          tempcoord(l)%em = 0
          tempcoord(l)%fm = 0
              tempcoord(l)%bondg = 0
    tempcoord(l)%bondf = 0
    bond(1) = 0

       call bonddirection(tempcoord(l)%bstate)
    
          call removeenergyintra(m,l,deltaenergy,tempcoord)
          call removeenergyinter(m,l,deltaenergy,tempcoord)

          do str = 1,l-3,1
             call adjacent(m,l,str,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                call bondformintra(m,l,str,tempcoord,bdir,deltaenergy,bond)
             end if
          end do
          do str = l+3,maxl,1
             call adjacent(m,l,str,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.)) then
                call bondformintra(m,l,str,tempcoord,bdir,deltaenergy,bond)
             end if
          end do

          do pr = 1,nprotein
             if(pr /= m) then
                maxback = chlen(pr)
                do str = 1,maxback
                   call adjacent(pr,l,str,tempcoord,adjver,bdir)
                   if((adjver.eqv. .true.)) then
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy,bond)

                   end if
                end do
             end if
          end do
 
          if(bond(1)>0) then
                    !write(6,*) 'vector right angle'
             call bondchoose(bond,tempcoord,deltaenergy,m,l)
              !         write(6,*) 'vector right angle end'
    end if
       
       !if(NINT(deltaenergy) /= 0) write(6,*) 'pivot',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then

          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l,l,tempcoord,tbdi,maxl)
          call updatebond(m,l,l,tempcoord)
          bdi(m,l-1) = tbdi(l-1)

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
        !call energy(.false.)
    write(6,*) 'rightangle'
call debug
  end subroutine rightanglemove

  subroutine rightvector(m,l,dx,dy,dz,racc,tbdi,maxl)
    integer,intent(in):: m,l,dx,dy,dz,maxl
    integer,dimension(:),intent(inout)::tbdi
    logical,intent(inout)::racc
    integer :: runn,k
    if(dx == 0) runn =1
    if(dy == 0) runn = 2
    if(dz == 0) runn =3
    racc = .true.
    


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

 
 subroutine endmove
    !perform reptation
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir,maxl,maxback,str
    integer,dimension(7)::bond
    logical :: reptcont,adjver
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbdi
!write(6,*) 'end pre move' 
!call energy(.false.)
    !write(6,*) 'for information',protcoords(20,6)%dm,protcoords(20,6)%dl
    t = time + 1

    !choose = -0.5
    choose =  ran2(seed) - 0.5
    !write(6,*) 'reptation start',choose
endatt = endatt+1
    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1
    maxl = chlen(m)
    !write(6,*) 'maxl',maxl
    allocate(tempcoord(maxl))
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
       else if((direc == (g1 +g2))  .and.  (g2 == 1)) then
          dx = -1
          tbdi(1) =1
       else if((direc == (g1+g2+g3)) .and.  (g3 == 1)) then
          dy = 1
          tbdi(1) = -2
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          tbdi(1) = 2
       else if((direc == (g1+g2+g3+g4+g5)) .and. (g5 == 1)) then
          dz = 1
          tbdi(1) = -3
       else if((direc == (g1+g2+g3+g4+g5+g6)) .and.  (g6 == 1)) then
          dz = -1
          tbdi(1) = 3
       end if
call bonddirection(tempcoord(1)%bstate)



       tempcoord(1)%x = modulo(protcoords(m,2)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,2)%y + dy-1,gridsize)+1
       tempcoord(1)%z = modulo(protcoords(m,2)%z + dz-1,gridsize)+1
       ! write(6,*) tbdi(1),dx,dy,dz
       ! write(6,*) tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z
       ! write(6,*) protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z

       deltaenergy = 0.0d0

       do st = 2,maxl,1
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

     tempcoord(1)%am = 0
       tempcoord(1)%bm = 0
       tempcoord(1)%cm = 0
       tempcoord(1)%dm = 0
       tempcoord(1)%em = 0
       tempcoord(1)%fm = 0
       tempcoord(1)%bondg = 0
       tempcoord(1)%bondf = 0
       bond(1) = 0
       call removeenergyintra(m,1,deltaenergy,tempcoord)
       call removeenergyinter(m,1,deltaenergy,tempcoord)

  

       do st = 4,maxl
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,1,st,tempcoord,bdir,deltaenergy,bond)
          end if
       end do

       !if(deltaenergy<0.0)write(6,*) 'energy post intra=',deltaenergy

       do pr = 1,nprotein
          if(pr /= m) then
             maxback = chlen(pr)
             do st = 1,maxback
                call adjacent(pr,1,st,tempcoord,adjver,bdir)
                if((adjver.eqv. .true.)) then
                   call bondforminter(m,pr,1,st,tempcoord,bdir,deltaenergy,bond)
                end if
             end do
          end if
       end do



       if(bond(1)>0) then
          call bondchoose(bond,tempcoord,deltaenergy,m,1)
       end if

       !if(deltaenergy<0.0)write(6,*) 'energy final=',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,1,1,tempcoord,tbdi,maxl)
          call updatebond(m,1,1,tempcoord)
          successful = successful + 1
          endacc = endacc + 1
          call counts(0,0,0,1)
       else
          reptcont = .false.
       end if
    else if (choose < 0.0) then

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
       else if((direc ==(g1 +g2)) .and.  (g2 == 1)) then
          dx = -1
          tbdi(maxl-1) = -1
       else if((direc == (g1+g2+g3)).and.  (g3 == 1)) then
          dy = 1
          tbdi(maxl-1) = 2
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
      tbdi(maxl-1) = -2
       else if((direc == (g1+g2+g3+g4+g5)) .and.( g5 == 1)) then
          dz = 1
         tbdi(maxl-1) = 3
       else if((direc == (g1+g2+g3+g4+g5+g6)).and.  (g6 == 1)) then
          dz = -1
          tbdi(maxl-1) = -3
       end if

       call bonddirection(tempcoord(maxl)%bstate)


       tempcoord(maxl)%x = modulo(protcoords(m,maxl-1)%x + dx-1,gridsize)+1
       tempcoord(maxl)%y = modulo(protcoords(m,maxl-1)%y + dy-1,gridsize)+1
       tempcoord(maxl)%z = modulo(protcoords(m,maxl-1)%z + dz-1,gridsize)+1

       !write(6,*) tbdi(maxl-1),dx,dy,dz
       !write(6,*) tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z
       !write(6,*) protcoords(m,maxl)%x,protcoords(m,maxl)%y,protcoords(m,maxl)%z

       deltaenergy = 0.0d0


       do st = 1,maxl-1,1
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

       tempcoord(maxl)%am = 0
       tempcoord(maxl)%bm = 0
       tempcoord(maxl)%cm = 0
       tempcoord(maxl)%dm = 0
       tempcoord(maxl)%em = 0
       tempcoord(maxl)%fm = 0
       tempcoord(maxl)%bondg = 0
       tempcoord(maxl)%bondf = 0
       bond(1) = 0
       call removeenergyintra(m,maxl,deltaenergy,tempcoord)
       call removeenergyinter(m,maxl,deltaenergy,tempcoord)



       do st = 1,maxl-3
          call adjacent(m,maxl,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,maxl,st,tempcoord,bdir,deltaenergy,bond)
          end if
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
             maxback = chlen(pr)
             do st = 1,maxback,1
                call adjacent(pr,maxl,st,tempcoord,adjver,bdir)
                if((adjver.eqv. .true.)) then
                   call bondforminter(m,pr,maxl,st,tempcoord,bdir,deltaenergy,bond)
                end if
             end do
          end if
       end do

 if(bond(1)>0) then
    call bondchoose(bond,tempcoord,deltaenergy,m,maxl)
    end if


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,maxl,maxl,tempcoord,tbdi,maxl)
          call updatebond(m,maxl,maxl,tempcoord)
            bdi(m,maxl-1)= tbdi(maxl-1)
            endacc = endacc + 1
          call counts(0,0,0,1)
       else
          reptcont = .false.
       end if
    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
       call counts(0,0,0,2)
    end if
!write(6,*) 'm',m,choose
    !call debug
   ! write(6,*) 'end start'
    !call energy(.false.)
    !write(6,*) 'end end'

        !call energy(.false.)
    write(6,*) 'end move'
    call debug
  end subroutine endmove

  

  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,tbdi,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    integer,dimension(:),intent(in)::tbdi
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum,g,f
    !moves beads to new positions and reassigns isobonding


    do beadnum = beadmin,beadmax,1

       if(protcoords(chainnum,beadnum)%bondg /=0) then
          g = protcoords(chainnum,beadnum)%bondg
          f = protcoords(chainnum,beadnum)%bondf
          protcoords(g,f)%bondg = 0
          protcoords(g,f)%bondf = 0
          protcoords(g,f)%bstate = 1
          end if

          
          protcoords(chainnum,beadnum)%bstate = tempcoord(beadnum)%bstate
       protcoords(chainnum,beadnum)%bondg = tempcoord(beadnum)%bondg
       protcoords(chainnum,beadnum)%bondf = tempcoord(beadnum)%bondf
       
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z

   if(protcoords(chainnum,beadnum)%bondg /=0) then
          g = protcoords(chainnum,beadnum)%bondg
          f = protcoords(chainnum,beadnum)%bondf
          protcoords(g,f)%bstate = 3
          protcoords(g,f)%bondg = chainnum
          protcoords(g,f)%bondf = beadnum
       end if
       
       if(beadnum<maxl) bdi(chainnum,beadnum) = tbdi(beadnum)
    end do

    call comfind(chainnum,.TRUE.)


    !call rms(chainnum)
  end subroutine updatepos


   subroutine updateposrep(chainnum,beadmin,beadmax,tempcoord,tbdi,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    integer,dimension(:),intent(in)::tbdi
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum,g,f
    !moves beads to new positions and reassigns isobonding


    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
       protcoords(chainnum,beadnum)%bstate = tempcoord(beadnum)%bstate
       protcoords(chainnum,beadnum)%bondg = tempcoord(beadnum)%bondg
       protcoords(chainnum,beadnum)%bondf = tempcoord(beadnum)%bondf

   if(protcoords(chainnum,beadnum)%bondg /=0) then
          g = protcoords(chainnum,beadnum)%bondg
          f = protcoords(chainnum,beadnum)%bondf
          protcoords(g,f)%bstate = 3
          protcoords(g,f)%bondg = chainnum
          protcoords(g,f)%bondf = beadnum
       end if
       
       if(beadnum<maxl) bdi(chainnum,beadnum) = tbdi(beadnum)
    end do

    call comfind(chainnum,.TRUE.)


    !call rms(chainnum)
  end subroutine updateposrep

  subroutine updatecluspos(tempcluscoord,ctbdi,tcom,cllist,clsize)
    !updates bead positions
    integer,dimension(:),intent(in) :: cllist
    integer,dimension(:,:),intent(in)::ctbdi
    integer,intent(in) :: clsize
    Type(prottemp),dimension(:,:),intent(in) :: tempcluscoord
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
          if(l<maxback) then
             bdi(a,l) = ctbdi(a,l)
             !if(ctbdi(a,l) == 0) write(6,*) 'update fail',a,l
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


  subroutine clusterassignmove(chainnum,cllist,clsize,cl)
    integer,intent(in) :: chainnum
    integer,dimension(:),intent(inout)::cllist
    logical,dimension(:),intent(inout)::cl
    integer,intent(inout):: clsize
    integer:: m,l,g,f,clcount,bdir,maxl,a,z
    logical :: clusyes

    
    clusyes = .true.

    cllist(1) = chainnum
    cl(chainnum) = .true.
    clsize = 1
    a = 1

    
    !do z = 1,clsize
z = 1
    do while(z<=a)
    m = cllist(z)
    !write(6,*) 'cluster assign output',m,a,cllist(a),clsize
    maxl = chlen(m)
    do l = 1,maxl
      if(protcoords(m,l)%am > 0) then
         call clusteraid(cllist,clsize,-1,a,m,l,cl)
      end if      
      if(protcoords(m,l)%bm > 0) then
         call clusteraid(cllist,clsize,1,a,m,l,cl)
      end if
      if(protcoords(m,l)%cm > 0) then
         call clusteraid(cllist,clsize,-2,a,m,l,cl)
      end if
      if(protcoords(m,l)%dm > 0) then
         call clusteraid(cllist,clsize,2,a,m,l,cl)
      end if
      if(protcoords(m,l)%em > 0) then
         call clusteraid(cllist,clsize,-3,a,m,l,cl)
      end if
      if(protcoords(m,l)%fm > 0) then
         call clusteraid(cllist,clsize,3,a,m,l,cl)
      end if
   end do
   z = z+1
end do



  end subroutine clusterassignmove

  subroutine clusteraid(cllist,cluspop,bdir,a,m,l,cl)
    integer,intent(in) :: bdir,m,l
    integer,intent(inout)::a,cluspop
    logical,dimension(:),intent(inout)::cl
    integer,dimension(:),intent(inout):: cllist
    integer :: f,g


    if(bdir == -1) then
    g = protcoords(m,l)%am
    f = protcoords(m,l)%al
 else if(bdir == 1) then
    g = protcoords(m,l)%bm
    f = protcoords(m,l)%bl
    else if(bdir == -2) then
    g = protcoords(m,l)%cm
    f = protcoords(m,l)%cl
    else if(bdir == 2) then
    g = protcoords(m,l)%dm
    f = protcoords(m,l)%dl
    else if(bdir == -3) then
    g = protcoords(m,l)%em
    f = protcoords(m,l)%el
    else if(bdir == 3) then
    g = protcoords(m,l)%fm
    f = protcoords(m,l)%fl
    end if
    

    if(cl(g) .eqv. .true.) goto 31
    
    if(isbond .eqv. .true.) then
       if(interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0) then
          a = a + 1
          cllist(a) = g
          cl(g) = .true.
          cluspop = cluspop + 1
          goto 31
       end if
    end if
    if((protcoords(m,l)%bondg == g) .and. (protcoords(m,l)%bondf ==f) .and. &
         (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0)) then
       if(cl(g) .eqv. .true.) goto 31
       a = a + 1
       cllist(a) = g
       cl(g) = .true.
       cluspop = cluspop + 1
    end if

  

31  continue

  end subroutine clusteraid




  subroutine clusterassign(chainnum,clno,clsize,cb)
    integer,intent(in) :: chainnum
    integer,dimension(:),intent(inout):: clno
    integer,dimension(:,:),intent(inout)::cb
    integer,intent(inout):: clsize
    integer:: m,l,g,f,clcount,bdir,maxlengthss,maxback,zzz,dum2,dum1,oldcl
    logical :: clusyes,adjver
    type(prottemp),dimension(:),allocatable :: tempcoord

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
                   if((protcoords(g,f)%bondg == m) .and. (protcoords(g,f)%bondf==l) .and. &
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
       !if(clno(m) == clno(chainnum))write(6,*) 'cluster info',time,m
    end do
  end subroutine clusterassign

  
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
       if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
    end if

    if(protcoords(m,l)%bm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%bl
       tempcoord(l)%bm = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
    end if

    if(protcoords(m,l)%cm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
    end if

    if(protcoords(m,l)%dm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%dl
       tempcoord(l)%dm = -1
        if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
    end if

    if(protcoords(m,l)%em ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%el
       tempcoord(l)%em = -1
                if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
    end if

    if(protcoords(m,l)%fm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%fl
       tempcoord(l)%fm = -1
    if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
          end if
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
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

    if((protcoords(m,l)%bm /= 0)  .and. (protcoords(m,l)%bm /= m)) then
       g = protcoords(m,l)%bm
       f = protcoords(m,l)%bl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%bm = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

    if((protcoords(m,l)%cm /= 0) .and. (protcoords(m,l)%cm /= m)) then
       g = protcoords(m,l)%cm
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

    if((protcoords(m,l)%dm /= 0) .and. (protcoords(m,l)%dm /= m)) then
       g = protcoords(m,l)%dm
       f = protcoords(m,l)%dl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%dm = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

    if((protcoords(m,l)%em /= 0) .and. (protcoords(m,l)%em /= m)) then
       g = protcoords(m,l)%em
       f = protcoords(m,l)%el
       !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1)) write(6,*) 'successful bond rupture',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%em = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

    if((protcoords(m,l)%fm /= 0) .and. (protcoords(m,l)%fm /= m)) then
       g = protcoords(m,l)%fm
       f = protcoords(m,l)%fl
       !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'successful rupture of bond',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%fm = -1
   if((protcoords(m,l)%bondg ==g)  .and. (protcoords(m,l)%bondf == f))  then
          deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
       end if
    end if

  end subroutine removeenergyinter

  subroutine bondformintra(m,l,st,tempcoord,bdir,deltaenergy,bond)
    integer,intent(in)::m,l,st,bdir
    double precision,intent(inout)::deltaenergy
    integer,dimension(:),intent(inout)::bond
    type(prottemp),dimension(:),intent(inout):: tempcoord
    integer::dir



    deltaenergy = deltaenergy + intraen
    
    if(bdir == -1) then
       tempcoord(l)%am = m
       tempcoord(l)%al = st
       dir = 1
    else if(bdir == 1) then
       tempcoord(l)%bm = m
       tempcoord(l)%bl = st
       dir = 2
    else if(bdir == -2) then
       tempcoord(l)%cm = m
       tempcoord(l)%cl = st
       dir = 3
    else if(bdir == 2) then
       tempcoord(l)%dm = m
       tempcoord(l)%dl = st
       dir = 4
    else if(bdir == -3) then
       tempcoord(l)%em = m
       tempcoord(l)%el = st
       dir = 5
    else if(bdir == 3) then
       tempcoord(l)%fm = m
       tempcoord(l)%fl = st
       dir = 6
    end if

    if(tempcoord(l)%bstate == 1) then
       if (intraenergy(protcoords(m,l)%type,protcoords(m,st)%type)< 0.0) then
          if(protcoords(m,st)%bstate == 1) then
             bond(1) = bond(1)+1
             bond(bond(1)+1) = dir
          else if(protcoords(m,st)%bstate ==3) then
             if((protcoords(m,st)%bondg == m) .and. (protcoords(m,st)%bondf == l)) then
                bond(1) = bond(1)+1
                bond(bond(1)+1) = dir
             end if
          end if
       end if
    end if
    
  end subroutine bondformintra

  subroutine bondforminter(m,g,l,st,tempcoord,bdir,deltaenergy,bond)
    integer,intent(in)::m,l,st,g,bdir
    double precision,intent(inout)::deltaenergy
    integer,dimension(:),intent(inout)::bond
    type(prottemp),dimension(:),intent(inout):: tempcoord
    integer::dir
    !check this and cluster bonds update


    deltaenergy = deltaenergy + interen(protcoords(m,l)%type,protcoords(g,st)%type)
    if(bdir == -1) then
       tempcoord(l)%am = g
       tempcoord(l)%al = st
       dir = 1
    else if(bdir == 1) then
       tempcoord(l)%bm = g
       tempcoord(l)%bl = st
       dir = 2
    else if(bdir == -2) then
       tempcoord(l)%cm = g
       tempcoord(l)%cl = st
       dir = 3
    else if(bdir == 2) then
       tempcoord(l)%dm = g
       tempcoord(l)%dl = st
       dir = 4
    else if(bdir == -3) then
       tempcoord(l)%em = g
       tempcoord(l)%el = st
       dir = 5
    else if(bdir == 3) then
       tempcoord(l)%fm = g
       tempcoord(l)%fl = st
       dir = 6
end if

if(tempcoord(l)%bstate == 1) then
   if (interenergy(protcoords(m,l)%type,protcoords(g,st)%type)< 0.0) then
      if(protcoords(g,st)%bstate == 1) then
         bond(1) = bond(1)+1
         bond(bond(1)+1) = dir
         !write(6,*) 'string',m,l,g,st,dir
      else if(protcoords(g,st)%bstate ==3) then
         if((protcoords(g,st)%bondg == m) .and. (protcoords(g,st)%bondf == l)) then
            bond(1) = bond(1)+1
            bond(bond(1)+1) = dir
         end if
      end if
   end if
   
end if
    !if(modulo(time,1000)==0) write(6,*) 'interbond',bdir,deltaenergy,tbo(l)
  end subroutine bondforminter

  subroutine updatebondother(m,bead,tempcoord)
    integer,intent(in):: m,bead
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f


    if(tempcoord(bead)%am == -1) then
       g = protcoords(m,bead)%am
       f = protcoords(m,bead)%al
       !write(6,*) 'info again 1/2',m,bead,g,f
       protcoords(g,f)%bl = 0
       protcoords(g,f)%bm = 0
    end if

    if(tempcoord(bead)%bm == -1) then
       g = protcoords(m,bead)%bm
       f = protcoords(m,bead)%bl
       !write(6,*) 'info again 2/2',m,bead
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

    if(protcoords(m,bead)%bstate == 3) then
       g = protcoords(m,bead)%bondg
       f = protcoords(m,bead)%bondf
       protcoords(g,f)%bstate = 1
       protcoords(g,f)%bondg = 0
       protcoords(g,f)%bondf = 0
       
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
          !write(6,*) 'update',g,f
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


       protcoords(m,l)%bstate = tempcoord(l)%bstate
       if(protcoords(m,l)%bstate ==3) then
       if(tempcoord(l)%bondg == m) then
          protcoords(m,l)%bondf = tempcoord(l)%bondf+del
       else
          protcoords(m,l)%bondf = tempcoord(l)%bondf
       end if
       if(tempcoord(l)%bondg > 0) then
          g = tempcoord(l)%bondg
          f = tempcoord(l)%bondf
          if(g ==m) f = f+del
          protcoords(g,f)%bondf = l
          protcoords(g,f)%bondg = m
          protcoords(g,f)%bstate = 3
       end if
       end if

       
    end do

  end subroutine updatebondrep


  subroutine reptation
    !perform reptation
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir,maxl,maxback,str
    integer,dimension(7)::bond
    logical :: reptcont,adjver
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbdi


    !write(6,*) 'for information',protcoords(20,6)%dm,protcoords(20,6)%dl
    t = time + 1

    !choose = -0.5
    choose =  ran2(seed) - 0.5
    !write(6,*) 'reptation start',choose

    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1
    maxl = chlen(m)
    !write(6,*) 'maxl',maxl
    allocate(tempcoord(maxl))
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
       else if((direc == (g1 +g2))  .and.  (g2 == 1)) then
          dx = -1
          tbdi(1) =1
       else if((direc == (g1+g2+g3)) .and.  (g3 == 1)) then
          dy = 1
          tbdi(1) = -2
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          tbdi(1) = 2
       else if((direc == (g1+g2+g3+g4+g5)) .and. (g5 == 1)) then
          dz = 1
          tbdi(1) = -3
       else if((direc == (g1+g2+g3+g4+g5+g6)) .and.  (g6 == 1)) then
          dz = -1
          tbdi(1) = 3
       end if
call bonddirection(tempcoord(1)%bstate)



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
       tempcoord(1)%bondg = 0
       tempcoord(1)%bondf = 0
       bond(1) = 0

       do st = 3,maxl-1
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,1,st,tempcoord,bdir,deltaenergy,bond)
          end if
       end do

       !if(deltaenergy<0.0)write(6,*) 'energy post intra=',deltaenergy

       do pr = 1,nprotein
          if(pr /= m) then
             maxback = chlen(pr)
             do st = 1,maxback
                call adjacent(pr,1,st,tempcoord,adjver,bdir)
                if((adjver.eqv. .true.)) then
                   call bondforminter(m,pr,1,st,tempcoord,bdir,deltaenergy,bond)
                end if
             end do
          end if
       end do

       !if(deltaenergy<0.0) write(6,*) 'energy post inter=',deltaenergy
       !section to adjust for change in types


       call reptationtype(m,1,deltaenergy,tempcoord,1)


       if(bond(1)>0) then
          call bondchoose(bond,tempcoord,deltaenergy,m,1)
       end if

       !if(deltaenergy<0.0)write(6,*) 'energy final=',deltaenergy
       if(Energydecision(deltaenergy) .eqv. .true.) then
!write(6,*) 'tempcoordinfo',tempcoord(1)%bstate,tempcoord(1)%bondg,tempcoord(1)%bondf
          totalenergy = totalenergy + deltaenergy
          do l =maxl,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z
             if(l<maxl) bdi(m,l) = bdi(m,l-1)


          end do
          !write(6,*) 'tbdi',tbdi(1)
   
          call updatebondother(m,maxl,tempcoord)
      call reptationadjust(m,maxl,2,tempcoord,-1)
          call updateposrep(m,1,1,tempcoord,tbdi,maxl)

          call updatebondrep(m,1,1,tempcoord,1)
          successful = successful + 1
          reptforward = reptforward + 1
          call counts(0,0,0,1)
       else
          reptcont = .false.
       end if
    else if (choose < 0.0) then

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
       else if((direc ==(g1 +g2)) .and.  (g2 == 1)) then
          dx = -1
          tbdi(maxl-1) = -1
       else if((direc == (g1+g2+g3)).and.  (g3 == 1)) then
          dy = 1
          tbdi(maxl-1) = 2
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
      tbdi(maxl-1) = -2
       else if((direc == (g1+g2+g3+g4+g5)) .and.( g5 == 1)) then
          dz = 1
         tbdi(maxl-1) = 3
       else if((direc == (g1+g2+g3+g4+g5+g6)).and.  (g6 == 1)) then
          dz = -1
          tbdi(maxl-1) = -3
       end if

       call bonddirection(tempcoord(maxl)%bstate)


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
       tempcoord(maxl)%bondg = 0
       tempcoord(maxl)%bondf = 0
       bond(1) = 0

       do st = 2,maxl-2
          call adjacent(m,maxl,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,maxl,st,tempcoord,bdir,deltaenergy,bond)
          end if
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
             maxback = chlen(pr)
             do st = 1,maxback,1
                call adjacent(pr,maxl,st,tempcoord,adjver,bdir)
                if((adjver.eqv. .true.)) then
                   call bondforminter(m,pr,maxl,st,tempcoord,bdir,deltaenergy,bond)
                end if
             end do
          end if
       end do

       call reptationtype(m,maxl,deltaenergy,tempcoord,-1)

 if(bond(1)>0) then
    call bondchoose(bond,tempcoord,deltaenergy,m,maxl)
    end if


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy

          do l =1,maxl-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z
             if(l<maxl-1) bdi(m,l) = bdi(m,l+1)


          end do
          !call checking(1)
          call updatebondother(m,1,tempcoord)
          call reptationadjust(m,1,maxl-1,tempcoord,1)
          call updateposrep(m,maxl,maxl,tempcoord,tbdi,maxl)
          call updatebondrep(m,maxl,maxl,tempcoord,-1)
            bdi(m,maxl-1)= tbdi(maxl-1)
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
!write(6,*) 'm',m,choose
    !call debug
    call energy(.false.)
    !write(6,*) 'reptation'
call debug
    
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
             !  end if
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
             !      protcoords(m,l)%bl = protcoords(m,l+del)%bl !-del
             !end if
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
       
       protcoords(m,l)%bondg = protcoords(m,l+del)%bondg
       protcoords(m,l)%bstate = protcoords(m,l+del)%bstate
       g= protcoords(m,l+del)%bondg
       if(g /=0) then
          if(g == m) then
             if(protcoords(m,l+del)%bondf /= buzz) then
                protcoords(m,l)%bondf = protcoords(m,l+del)%bondf - del
                
             else
                protcoords(m,l)%bondg = 0
                protcoords(m,l)%bondf= 0
                !protcoords(m,l)%bstate = 0
             end if
          else if(g/=m) then
             protcoords(m,l)%bondf = protcoords(m,l+del)%bondf
             f = protcoords(m,l)%bondf
             protcoords(g,f)%bondf = l
          end if
       else
          protcoords(m,l)%bondg =0
          protcoords(m,l)%bondf = 0
       end if
       
    end do
  end subroutine reptationadjust

  subroutine reptationtype(m,st,deltaenergy,tempcoord,del)
    integer,intent(in)::m,st,del
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout)::tempcoord
    integer :: g,f,l,rangel,rangeh,maxl,del2,buzz

    !         del2 = -abs(del)

    maxl = chlen(m)
    if(tempcoord(st)%bstate ==3) then
       !write(6,*) 'ddddebug',tempcoord(st)%am,st,tempcoord(st)%al
       g = tempcoord(st)%bondg
       if(g == m) then
          f = tempcoord(st)%bondf
          !if((protcoords(g,f)%bondg==m) .and. protcoords(g,f)%bondf==st) then
             deltaenergy = deltaenergy + intraenergy(protcoords(m,st)%type,protcoords(g,f+del)%type)
             deltaenergy = deltaenergy - intraenergy(protcoords(m,st)%type,protcoords(g,f)%type)
          !end if
       end if
    end if

   

    if (del == 1) then
       rangel = 1
       buzz = maxl
       rangeh = maxl-1
       del2 = 1
    else if (del == -1) then
       rangel = 2
       buzz = 1
       rangeh = maxl
       del2 = -1
    end if

    do l = rangel,rangeh,1

       if(protcoords(m,l)%bstate == 3) then             
          g = protcoords(m,l)%bondg
          f = protcoords(m,l)%bondf
             if(g == m) then
                if((f /=buzz) .and. (l<f)) then
                   deltaenergy = deltaenergy + intraenergy(protcoords(m,l+del2)%type,protcoords(g,f+del2)%type)
                   deltaenergy = deltaenergy - intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                end if
             else if(g/=m) then
                deltaenergy = deltaenergy + interenergy(protcoords(m,l+del2)%type,protcoords(g,f)%type)
                deltaenergy = deltaenergy - interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
             end if
       end if

    end do
  end subroutine reptationtype
  
   subroutine bonddirection(a)
    integer,intent(inout)::a
    
    a = int(ran2(seed)*2)+1

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

             checklist(zx) = a
             zx = zx +1
          end if
16        if(a == a) continue
       end do

    end do



    call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
  

76  if(clsize == 1) then
       maxl = chlen(m)
       centralchain = m
       totcluspop = maxl
       ! write(6,*) 'individual chain',com(1,m)%x,com(1,m)%y,com(1,m)%z
    end if

    comcluster%x = INT(modulo((com(1,centralchain)%x) + (dummycom%x/totcluspop)-1,real(gridsize))+1)

    comcluster%y = INT(modulo((com(1,centralchain)%y) + (dummycom%y/totcluspop)-1,real(gridsize))+1)

    comcluster%z = INT(modulo((com(1,centralchain)%z) + (dummycom%z/totcluspop)-1,real(gridsize))+1)

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


       if(rcounter(a)>1) then
          do rc = 1,rcounter(a)-1,1
             bx  = bx + comx(route(a,rc),route(a,rc+1))
             by = by + comy(route(a,rc),route(a,rc+1))
             bz = bz + comz(route(a,rc),route(a,rc+1))
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



    dummycom%x = modulo((com(1,m)%x) + (comx/totcluspop)-1,real(gridsize))+1
    comcluster%x = INT(dummycom%x)
    dummycom%y = modulo((com(1,m)%y) + (comy/totcluspop)-1,real(gridsize))+1
    comcluster%y = INT(dummycom%y)
    dummycom%z = modulo((com(1,m)%z) + (comz/totcluspop)-1,real(gridsize))+1
    comcluster%z = INT(dummycom%z)


  end subroutine clustercomcheck


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
    double precision :: initialenergy,olderenergy
    Type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:),allocatable :: bond
    logical :: adjver
    logical,intent(in)::init
    initialenergy = 0.0d0
    !totalenergy = 0.0d0
    allocate(tempcoord(maxlength))



    if(init .eqv. .true.) then
       do m = 1,nprotein
          maxl = chlen(m)
          allocate(bond(maxl,7))
          do l = 1,maxl

             bond(l,1) = 0
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
                      call bondformintra(m,l,f,tempcoord,bdir,initialenergy,bond(l,:))
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
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy,bond(l,:))

                      end if
                   end do
                end do
             end if
          end do
do l = 1,maxl
   if(bond(l,1)>0) then
             call bondchoose(bond(l,:),tempcoord,initialenergy,m,l)
          end if
          end do
          call updatebondold(m,1,maxl,tempcoord,bond)
          deallocate(bond)
       end do
    end if

    if(init .eqv. .false.) then
       do m =1,nprotein,1
          maxl = chlen(m)
                    do l =1,maxl,1


                       if(protcoords(m,l)%bstate ==3) then
                          g = protcoords(m,l)%bondg
                          f = protcoords(m,l)%bondf

                          if(g == m) then
                          initialenergy = initialenergy+intraenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                       else if(g/=m) then
                          initialenergy = initialenergy+interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                          end if
            
                       end if

                    end do
                 end do
                 initialenergy = initialenergy/2


  if(isbond .eqv. .true.) then           
                 do m =1,nprotein,1
                    maxl = chlen(m)
                    
                    allocate(bond(maxl,7))
                    do l =1,maxl,1
          bond(l,1) = 0
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
          do l = 1,maxl-3,1
             do f = l+3,maxl
                call adjacent(m,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.)) then
                      call bondformintra(m,l,f,tempcoord,bdir,initialenergy,bond(l,:))
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
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy,bond(l,:))
                      end if
                   end do
                end do
             end do
          end if
          deallocate(bond)
       end do
       end if
    end if


write(6,*) totalenergy

    

    if(init .eqv. .true.) initialenergy = initialenergy/2
    if((time > 0) .and. (initialenergy /= totalenergy)) then
       write(6,*) 'energyfail',totalenergy,initialenergy,time
       finalfail = .true.
    end if

    if(init .eqv. .true.) totalenergy = initialenergy


  end subroutine energy

  
  subroutine radiusofgyration
    integer::m,l,maxl,nobeads
    double precision :: rog,totrog
    totrog = 0.0d0
nobeads = 0
    do m = 1,nprotein,1
       maxl = chlen(m)
       do l = 1, maxl,1
          rog = (min(modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0))))**2
          totrog = totrog + rog
       end do
       nobeads = nobeads + maxl
    end do

    runningaveROG = runningaveROG + totrog
    polymerrog = polymerrog + (totrog/nprotein)
    if(noinfo .eqv. .false.) write(97,*) time-equilib, &
         totrog/nobeads,SQRT(totrog/nobeads)

  end subroutine radiusofgyration



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

   
    totrmsbrute = totrmsbrute + msdx + msdy + msdz

  end subroutine rms

 subroutine dataout
    integer::m,l,maxl,f,g,maxback,acount
    integer :: choice,dx,dy,dz
    real::xn,yn,zn
    
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


          if(protcoords(m,l)%bstate ==3) then
             g = protcoords(m,l)%bondg
             f = protcoords(m,l)%bondf
             dx = (protcoords(m,l)%x -protcoords(g,f)%x)
             if(dx>1) dx = dx-gridsize
             if(dx<-1) dx = gridsize +dx
             dy = (protcoords(m,l)%y -protcoords(g,f)%y)
             if(dy>1) dy = dy-gridsize
             if(dy<-1) dy = gridsize +dy
             dz = (protcoords(m,l)%z -protcoords(g,f)%z)
             if(dz>1) dz = dz-gridsize
             if(dz<-1) dz = gridsize +dz
             
             xn = modulo(protcoords(m,l)%x -(real(dx)/2.0)-1,real(gridsize))+1
             yn = modulo(protcoords(m,l)%y -(real(dy)/2.0)-1,real(gridsize))+1
             zn = modulo(protcoords(m,l)%z -(real(dz)/2.0)-1,real(gridsize))+1
             
             write(choice,*) acount, 4*xn, 4*yn,4*zn
             goto 71           
          end if

          write(choice,*) acount,real(4*protcoords(m,l)%x), real(4*protcoords(m,l)%y),real(4*protcoords(m,l)%z) 

71        if(m == m) continue
          acount = acount + 1
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
            protcoords(m,1)%type,protcoords(m,1)%bstate,protcoords(m,1)%bondg,protcoords(m,1)%bondf,maxl
       do l = 2,maxl,1
          write(59,*) protcoords(m,l)%species,protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z,protcoords(m,l)%type,&
               protcoords(m,l)%bstate,protcoords(m,l)%bondg,protcoords(m,l)%bondf
       end do
    end do


  end subroutine dataoutrestart

  
  subroutine clusbondform(m,g,l,f,tempcoord,bdir)
    integer,intent(in)::m,g,l,f,bdir
    type(prottemp),dimension(:,:),intent(inout) :: tempcoord
    if(bdir == -1) then
    tempcoord(m,l)%am = g
    tempcoord(m,l)%al = f
    else if(bdir == 1) then
    tempcoord(m,l)%bm = g
    tempcoord(m,l)%bl = f    
 else if(bdir == -2) then
    tempcoord(m,l)%cm = g
    tempcoord(m,l)%cl = f
     else if(bdir == 2) then
    tempcoord(m,l)%dm = g
    tempcoord(m,l)%dl = f
     else if(bdir == -3) then
    tempcoord(m,l)%em = g
    tempcoord(m,l)%el = f
     else if(bdir == 3) then
    tempcoord(m,l)%fm = g
    tempcoord(m,l)%fl = f
    end if

  end subroutine clusbondform

    subroutine clusterenergy(tempcluscoord,clusen,m,clsize,cllist,cl)
    integer,dimension(:),intent(in) ::cllist
    logical,dimension(:),intent(in)::cl
    integer,intent(in) :: clsize,m
    double precision,intent(inout)::clusen
    type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer:: a,l,g,f,bdir,maxlengthss,maxl,ra
    logical :: clusmove,adjver
    allocate(tempcoord(maxlength))
    !error in this!!!

    !write(6,*) 'realsanity',m,l,tempcluscoord(472,4)%x,tempcluscoord(472,4)%y,tempcluscoord(472,4)%z
    clusen = 0.0
    !write(6,*) 'clsize',clsize
    do ra = 1,clsize
       a = cllist(ra)
       !write(6,*) 'a',a
       maxl = chlen(a)
       do g = 1,nprotein
          if(cl(g) .eqv. .true.) goto 23
          maxlengthss = chlen(g)
          clusmove = .true.
          do l = 1,maxl,1
             !tempcoord(l)%x = tempcluscoord(a,l)%x
             !tempcoord(l)%y = tempcluscoord(a,l)%y
             !tempcoord(l)%z = tempcluscoord(a,l)%z
             !write(6,*) 'real life',m,l,tempcluscoord(472,4)%x,tempcluscoord(472,4)%y,tempcluscoord(472,4)%z
             do f = 1,maxlengthss,1
                call adjacent(g,l,f,tempcluscoord(a,:),adjver,bdir)
                if(adjver.eqv. .true.) then
                   call clusbondform(a,g,l,f,tempcluscoord,bdir)
                   if((isbond .eqv. .true.) .and. (interen(protcoords(a,l)%type,protcoords(g,f)%type) < 0.0)) then
                      !energyofcluster = energyofcluster + interen
                      clusen = 100
                      return
                   end if
                   if((protcoords(a,l)%bondg ==g) .and. (protcoords(a,l)%bondf==f) .and. &
                        (interenergy(protcoords(a,l)%type,protcoords(g,f)%type)  < 0.0)) then
                      !energyofcluster = energyofcluster + interenergy(protcoords(a,l)%species,protcoords(g,f)%species)
                      clusen = 100
                      return
                   end if
                end if
             end do
          end do
23        continue
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
    write(91,*) t,normchainl,sqrt(normchainl)
    !write(91,*) t,avechainlength,actualchainlength/((time-equilib)/10000)
  end subroutine length

 

  subroutine clustertranslation(m,clsize,cb,cllist,cl)
    integer,intent(inout)::m
    integer,intent(inout):: clsize
    logical,dimension(:),intent(in)::cl
    type(centremass),dimension(:),allocatable:: tcom
    integer,dimension(:),intent(inout) ::cllist
    integer,dimension(:,:),intent(inout) ::cb
    integer :: dx,dy,dz,vv,rvv,a
    type(prottemp),dimension(:,:),allocatable :: tempcluscoord
    integer,dimension(:,:),allocatable::ctbdi
    double precision :: clusen,direction
    logical :: moveallow
    allocate(tempcluscoord(nprotein,maxlength))
    allocate(ctbdi(nprotein,maxlength-1))
    allocate(tcom(nprotein))
    !return
    !write(6,*) 'translate', totalenergy
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

!write(6,*) 'info',dx,dy,dz
    call clusremove(cllist,clsize,tempcluscoord,cl)
!write(6,*) 'corresponding info',tempcluscoord(464,6)%x,tempcluscoord(464,6)%y,tempcluscoord(464,6)%z
    
    do rvv = 1,clsize
       vv = cllist(rvv)
       call clustermove(vv,dx,dy,dz,tempcluscoord,moveallow,ctbdi,cllist,cl)
       if(moveallow .eqv. .false.) return
       tcom(vv)%x = modulo(com(1,vv)%x + dx -1,real(gridsize))+1
       tcom(vv)%y = modulo(com(1,vv)%y + dy -1,real(gridsize))+1
       tcom(vv)%z = modulo(com(1,vv)%z + dz -1,real(gridsize))+1
    end do
 
    clusen = 0.0

!write(6,*) 'corresponding infopre',tempcluscoord(464,6)%x

!write(6,*) 'intrigue',protcoords(472,4)%dm
!write(6,*) 'protein 6',6,tempcluscoord(472,4)
call clusterenergy(tempcluscoord,clusen,m,clsize,cllist,cl)

!write(6,*) 'all the temporarys',1,tempcluscoord(472,4)
!write(6,*) 'all the temporarys 331111',tempcluscoord(331,4)
!write(6,*) 'all the temporarys',3,tempcluscoord(464,3)
!write(6,*) 'all the temporarys',4,tempcluscoord(464,4)
!write(6,*) 'all the temporarys',5,tempcluscoord(464,5)
!write(6,*) 'all the temporarys',6,tempcluscoord(464,6)


!write(6,*) 'corresponding info after en',tempcluscoord(14,3)%am
    if(clusen == 0.0) then
       call updatecluspos(tempcluscoord,ctbdi,tcom,cllist,clsize)
       call updateclusbonds(clsize,cllist,cl,tempcluscoord)
       !call energy(.true.)
       cltacc = cltacc + 1
    end if


        call energy(.false.)
    write(6,*) 'cluster'
!write(6,*) 'info poster', protcoords(304,4)%am,protcoords(304,4)%bm
    !write(6,*) 'b',clsize,cllist(1:clsize)
call debug


  end subroutine clustertranslation

  subroutine clusremove(cllist,clsize,tempcluscoord,cl)
  Type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    logical,dimension(:),intent(in)::cl
    integer,dimension(:),intent(inout) ::cllist
    integer,intent(in)::clsize
    integer::rvv,m,l,maxl
    integer :: g,f

    do rvv = 1,clsize
       m = cllist(rvv)
       maxl = chlen(m)
       do l = 1,maxl
          if((protcoords(m,l)%am ==0)) then
             tempcluscoord(m,l)%am = 0
             tempcluscoord(m,l)%al = 0
          else
             if(cl(protcoords(m,l)%am) .eqv. .false.) then
                tempcluscoord(m,l)%am = -1
                !g = protcoords(m,l)%am
                !f = protcoords(m,l)%al
                !if((bonddd(m,l) ==1) .and. (bonddd(protcoords(m,l)%am,protcoords(m,l)%al) == -1)&
                 !    .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                 !  write(6,*) 'bond break!!!',m,l,protcoords(m,l)%am,cllist(1:clsize),'a'
                 !  end if
             else if(cl(protcoords(m,l)%am) .eqv. .true.) then
                tempcluscoord(m,l)%am = protcoords(m,l)%am
                tempcluscoord(m,l)%al = protcoords(m,l)%al
             end if
          end if

          if((protcoords(m,l)%bm ==0)) then
             tempcluscoord(m,l)%bm = 0
             tempcluscoord(m,l)%bl = 0
          else
             if(cl(protcoords(m,l)%bm) .eqv. .false.) then
                tempcluscoord(m,l)%bm = -1
                !g = protcoords(m,l)%bm
                !f = protcoords(m,l)%bl
                 !if((bonddd(m,l) ==-1) .and. (bonddd(protcoords(m,l)%bm,protcoords(m,l)%bl) == 1)&
                 !    .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                 !  write(6,*) 'bond break!!!',m,l,protcoords(m,l)%bm,cllist(1:clsize),'b'
                 !  end if
             else if(cl(protcoords(m,l)%bm) .eqv. .true.) then
                tempcluscoord(m,l)%bm = protcoords(m,l)%bm
                tempcluscoord(m,l)%bl = protcoords(m,l)%bl
             end if
          end if

          if((protcoords(m,l)%cm ==0)) then
             tempcluscoord(m,l)%cm = 0
             tempcluscoord(m,l)%cl = 0
          else if((protcoords(m,l)%cm /=0)) then
             if(cl(protcoords(m,l)%cm) .eqv. .false.) then
                tempcluscoord(m,l)%cm = -1
                !g = protcoords(m,l)%cm
                !f = protcoords(m,l)%cl
                 !if((bonddd(m,l) ==2) .and. (bonddd(protcoords(m,l)%cm,protcoords(m,l)%cl) == -2)&
                 !    .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                 !  write(6,*) 'bond break!!!',m,l,protcoords(m,l)%cm,cllist(1:clsize),'c'
                 !  end if
             else if(cl(protcoords(m,l)%cm) .eqv. .true.) then
                tempcluscoord(m,l)%cm = protcoords(m,l)%cm
                tempcluscoord(m,l)%cl = protcoords(m,l)%cl
             end if
          end if

          if((protcoords(m,l)%dm ==0)) then
             tempcluscoord(m,l)%dm = 0
             tempcluscoord(m,l)%dl = 0
          else if((protcoords(m,l)%dm /=0)) then
             if(cl(protcoords(m,l)%dm) .eqv. .false.) then
                tempcluscoord(m,l)%dm = -1
                !g = protcoords(m,l)%dm
                !f = protcoords(m,l)%dl
                 !if((bonddd(m,l) ==-2) .and. (bonddd(protcoords(m,l)%dm,protcoords(m,l)%dl) == 2)&
                 !    .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                 !   write(6,*) 'bond break!!!',m,l,protcoords(m,l)%dm,cllist(1:clsize),'d'
                 !end if
             else if(cl(protcoords(m,l)%dm) .eqv. .true.) then
                tempcluscoord(m,l)%dm = protcoords(m,l)%dm
                tempcluscoord(m,l)%dl = protcoords(m,l)%dl
             end if
          end if

          if((protcoords(m,l)%em ==0)) then
             tempcluscoord(m,l)%em = 0
             tempcluscoord(m,l)%el = 0
          else
             if(cl(protcoords(m,l)%em) .eqv. .false.) then
                tempcluscoord(m,l)%em = -1
                !g = protcoords(m,l)%em
                !f = protcoords(m,l)%el
                !if((bonddd(m,l) ==3) .and. (bonddd(protcoords(m,l)%em,protcoords(m,l)%el) == -3)&
                 !    .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                  ! write(6,*) 'bond break!!!',m,l,protcoords(m,l)%em,cllist(1:clsize),'e'
                  ! end if
             else if(cl(protcoords(m,l)%em) .eqv. .true.) then
                tempcluscoord(m,l)%em = protcoords(m,l)%em
                tempcluscoord(m,l)%el = protcoords(m,l)%el
             end if
          end if

          if((protcoords(m,l)%fm ==0)) then
             tempcluscoord(m,l)%fm = 0
             tempcluscoord(m,l)%fl = 0
          else
             if(cl(protcoords(m,l)%fm) .eqv. .false.) then
                tempcluscoord(m,l)%fm = -1
                !g = protcoords(m,l)%fm
                !f = protcoords(m,l)%fl
                 !if((bonddd(m,l) ==-3) .and. (bonddd(protcoords(m,l)%fm,protcoords(m,l)%fl) == 3)&
                  !   .and. (interen(protcoords(m,l)%type,protcoords(g,f)%type) < 0.0)) then
                  ! write(6,*) 'bond break!!!',m,l,protcoords(m,l)%fm,cllist(1:clsize),'f'
                   !end if
             else if(cl(protcoords(m,l)%fm) .eqv. .true.) then
                tempcluscoord(m,l)%fm = protcoords(m,l)%fm
                tempcluscoord(m,l)%fl = protcoords(m,l)%fl
             end if
          end if

   
         end do
         end do
    

  
  end subroutine clusremove

  
  

  subroutine updateclusbonds(clsize,cllist,cl,tempcluscoord)
    integer :: m,mt,maxl,l,g,f,gt,ft
    integer,intent(inout)::clsize
    integer,dimension(:),intent(inout) ::cllist
    logical,dimension(:),intent(in)::cl
    type(prottemp),dimension(:,:),intent(in):: tempcluscoord


    do mt = 1,clsize
       m = cllist(mt)
       maxl = chlen(m)
       do l = 1,maxl

          if(tempcluscoord(m,l)%am /= -1) then
             if(tempcluscoord(m,l)%am /= 0) then
                if(cl(tempcluscoord(m,l)%am) .eqv. .false.) then
                   g = tempcluscoord(m,l)%am
                   f = tempcluscoord(m,l)%al
                   !write(6,*) 'g and f',g,f
                   if(protcoords(m,l)%am /= 0) then
                      if((g/=protcoords(m,l)%am) .or. (f /=protcoords(m,l)%al)) then
                         gt = protcoords(m,l)%am
                         ft = protcoords(m,l)%al
                         if((protcoords(gt,ft)%bm == m) .and. (protcoords(gt,ft)%bl == l)) then
                            protcoords(gt,ft)%bm = 0
                            protcoords(gt,ft)%bl = 0
                         end if
                      end if
                   end if
                   protcoords(g,f)%bm = m
                   protcoords(g,f)%bl = l
                   !write(6,*) 'confusion',g,f,protcoords(g,f)%bm
                end if
             end if
             protcoords(m,l)%am = tempcluscoord(m,l)%am
             protcoords(m,l)%al = tempcluscoord(m,l)%al
          else if(tempcluscoord(m,l)%am == -1) then
             g = protcoords(m,l)%am
             f = protcoords(m,l)%al
             if((protcoords(g,f)%bm == m) .and. (protcoords(g,f)%bl == l)) then
                protcoords(g,f)%bm = 0
                protcoords(g,f)%bl = 0
             end if
             protcoords(m,l)%am = 0
             protcoords(m,l)%al = 0
          end if
!write(6,*) 'why??', protcoords(14,5)%am,tempcluscoord(14,5)%am
!write(6,*) 'why?? two', protcoords(50,3)%bm
!write(6,*) 'help', protcoords(182,4)



          if(tempcluscoord(m,l)%bm /= -1) then
             if(tempcluscoord(m,l)%bm /= 0) then
                if(cl(tempcluscoord(m,l)%bm) .eqv. .false.) then
                   g = tempcluscoord(m,l)%bm
                   f = tempcluscoord(m,l)%bl
                   if(protcoords(m,l)%bm /= 0) then
                      if((g/=protcoords(m,l)%bm) .or. (f /=protcoords(m,l)%bl)) then
                         gt = protcoords(m,l)%bm
                         ft = protcoords(m,l)%bl
                         if((protcoords(gt,ft)%am == m) .and. (protcoords(gt,ft)%al == l)) then
                            protcoords(gt,ft)%am = 0
                            protcoords(gt,ft)%al = 0
                         end if
                      end if
                   end if
                   protcoords(g,f)%am = m
                   protcoords(g,f)%al = l
                end if
             end if
             protcoords(m,l)%bm = tempcluscoord(m,l)%bm
             protcoords(m,l)%bl = tempcluscoord(m,l)%bl
          else if(tempcluscoord(m,l)%bm == -1) then
             g = protcoords(m,l)%bm
             f = protcoords(m,l)%bl
             if((protcoords(g,f)%am == m) .and. (protcoords(g,f)%al == l)) then
             protcoords(g,f)%am = 0
             protcoords(g,f)%al = 0
             end if
             protcoords(m,l)%bm = 0
             protcoords(m,l)%bl = 0
          end if


           if(tempcluscoord(m,l)%cm /= -1) then
             if(tempcluscoord(m,l)%cm /= 0) then
                if(cl(tempcluscoord(m,l)%cm) .eqv. .false.) then
                   g = tempcluscoord(m,l)%cm
                   f = tempcluscoord(m,l)%cl
                   if(protcoords(m,l)%cm /= 0) then
                   if((g/=protcoords(m,l)%cm) .or. (f /=protcoords(m,l)%cl)) then
                      gt = protcoords(m,l)%cm
                      ft = protcoords(m,l)%cl
                      if((protcoords(gt,ft)%dm == m) .and. (protcoords(gt,ft)%dl == l)) then
                      protcoords(gt,ft)%dm = 0
                      protcoords(gt,ft)%dl = 0
                   end if
                end if
                end if
                   protcoords(g,f)%dm = m
                   protcoords(g,f)%dl = l
                end if
             end if
             protcoords(m,l)%cm = tempcluscoord(m,l)%cm
             protcoords(m,l)%cl = tempcluscoord(m,l)%cl
          else if(tempcluscoord(m,l)%cm == -1) then
             g = protcoords(m,l)%cm
             f = protcoords(m,l)%cl
             if((protcoords(g,f)%dm == m) .and. (protcoords(g,f)%dl == l)) then
             protcoords(g,f)%dm = 0
             protcoords(g,f)%dl = 0
             end if
             protcoords(m,l)%cm = 0
             protcoords(m,l)%cl = 0
          end if


           if(tempcluscoord(m,l)%dm /= -1) then
             if(tempcluscoord(m,l)%dm /= 0) then
                if(cl(tempcluscoord(m,l)%dm) .eqv. .false.) then
                   g = tempcluscoord(m,l)%dm
                   f = tempcluscoord(m,l)%dl
                   if(protcoords(m,l)%dm /= 0) then
                   if((g/=protcoords(m,l)%dm) .or. (f /=protcoords(m,l)%dl)) then
                      gt = protcoords(m,l)%dm
                      ft = protcoords(m,l)%dl
                      if((protcoords(gt,ft)%cm == m) .and. (protcoords(gt,ft)%cl == l)) then
                      protcoords(gt,ft)%cm = 0
                      protcoords(gt,ft)%cl = 0
                   end if
                   end if
                   end if
                   protcoords(g,f)%cm = m
                   protcoords(g,f)%cl = l
                end if
             end if
             protcoords(m,l)%dm = tempcluscoord(m,l)%dm
             protcoords(m,l)%dl = tempcluscoord(m,l)%dl
          else if(tempcluscoord(m,l)%dm == -1) then
             g = protcoords(m,l)%dm
             f = protcoords(m,l)%dl
             if((protcoords(g,f)%cm == m) .and. (protcoords(g,f)%cl == l)) then
             protcoords(g,f)%cm = 0
             protcoords(g,f)%cl = 0
             end if
             protcoords(m,l)%dm = 0
             protcoords(m,l)%dl = 0
          end if

          if(tempcluscoord(m,l)%em /= -1) then
             if(tempcluscoord(m,l)%em /= 0) then
                if(cl(tempcluscoord(m,l)%em) .eqv. .false.) then
                   g = tempcluscoord(m,l)%em
                   f = tempcluscoord(m,l)%el
                   if(protcoords(m,l)%em /= 0) then
                   if((g/=protcoords(m,l)%em) .or. (f /=protcoords(m,l)%el)) then
                      gt = protcoords(m,l)%em
                      ft = protcoords(m,l)%el
                      if((protcoords(gt,ft)%fm == m) .and. (protcoords(gt,ft)%fl == l)) then
                      protcoords(gt,ft)%fm = 0
                      protcoords(gt,ft)%fl = 0
                   end if
                   end if
                   end if
                   protcoords(g,f)%fm = m
                   protcoords(g,f)%fl = l
                end if
             end if
             protcoords(m,l)%em = tempcluscoord(m,l)%em
             protcoords(m,l)%el = tempcluscoord(m,l)%el
          else if(tempcluscoord(m,l)%em == -1) then
             g = protcoords(m,l)%em
             f = protcoords(m,l)%el
             if((protcoords(g,f)%fm == m) .and. (protcoords(g,f)%fl == l)) then
             protcoords(g,f)%fm = 0
             protcoords(g,f)%fl = 0
             end if
             protcoords(m,l)%em = 0
             protcoords(m,l)%el = 0
          end if


           if(tempcluscoord(m,l)%fm /= -1) then
             if(tempcluscoord(m,l)%fm /= 0) then
                if(cl(tempcluscoord(m,l)%fm) .eqv. .false.) then
                   g = tempcluscoord(m,l)%fm
                   f = tempcluscoord(m,l)%fl
                   if(protcoords(m,l)%fm /= 0) then
                   if((g/=protcoords(m,l)%fm) .or. (f /=protcoords(m,l)%fl)) then
                      gt = protcoords(m,l)%fm
                      ft = protcoords(m,l)%fl
                      if((protcoords(gt,ft)%em == m) .and. (protcoords(gt,ft)%el == l)) then
                      protcoords(gt,ft)%em = 0
                      protcoords(gt,ft)%el = 0
                   end if
                   end if
                   end if
                   protcoords(g,f)%em = m
                   protcoords(g,f)%el = l
                end if
             end if
             protcoords(m,l)%fm = tempcluscoord(m,l)%fm
             protcoords(m,l)%fl = tempcluscoord(m,l)%fl
          else if(tempcluscoord(m,l)%fm == -1) then
             g = protcoords(m,l)%fm
             f = protcoords(m,l)%fl
             if((protcoords(g,f)%em == m) .and. (protcoords(g,f)%el == l)) then
             protcoords(g,f)%em = 0
             protcoords(g,f)%el = 0
             end if
             protcoords(m,l)%fm = 0
             protcoords(m,l)%fl = 0
          end if
          
          
       end do
       end do
                
  end subroutine updateclusbonds

  subroutine tempcomfind(m,tcom,tempcluscoord)
    integer :: dcomx,dcomy,dcomz,maxl,comx,comy,comz
    double precision :: dcx,dcy,dcz
    integer :: l
    integer,intent(inout)::m
    type(prottemp),dimension(:,:),intent(inout):: tempcluscoord
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
    dummycom%x = modulo((chlen(m)*tcom(m)%x) + (comx/totcluspop)-1,real(gridsize))+1
    comclustemp%x = NINT(dummycom%x)
    !write(6,*) 'x coords',(com(1,m)%x + (comx/totcluspop)),comcluster%x
    dummycom%y = modulo((chlen(m)*tcom(m)%y) + (comy/totcluspop)-1,real(gridsize))+1
    comclustemp%y = NINT(dummycom%y)
    !write(6,*) 'y coords',(com(1,m)%y + (comx/totcluspop)),comcluster%y
    dummycom%z = modulo((chlen(m)*tcom(m)%z) + (comz/totcluspop)-1,real(gridsize))+1
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


    comclustemp%x = NINT(modulo((chlen(centralchain)*tcom(centralchain)%x) + &
         (dummycom%x/totcluspop)-1,real(gridsize))+1)
    comclustemp%y = NINT(modulo((chlen(centralchain)*tcom(centralchain)%y) + &
         (dummycom%y/totcluspop)-1,real(gridsize))+1)
    comclustemp%z = NINT(modulo((chlen(centralchain)*tcom(centralchain)%z) + &
         (dummycom%z/totcluspop)-1,real(gridsize))+1)


53  if(compass .eqv. .true.) continue

  end subroutine tempclustercom

  

 
  subroutine tempclustermove(chainno,delx,dely,delz,tempcluscoord,moveallow,cllist)
    integer,intent(in) :: chainno,delx,dely,delz
    logical,intent(inout)::moveallow
    Type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::cllist
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
          if(ANY(cllist .eq. pr)) goto 39
          maxback = chlen(pr)
          do st = 1,maxback,1
             if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
                clusmove = .false.
                !write(6,*) 'translation fail'
                goto 31
             end if
          end do
39        continue
       end do

    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine tempclustermove

  subroutine clustermove(chainno,delx,dely,delz,tempcluscoord,moveallow,ctbdi,cllist,cl)
    integer,intent(in) :: chainno,delx,dely,delz
    integer,dimension(:,:),intent(inout)::ctbdi
    logical,dimension(:),intent(in)::cl
    logical,intent(inout)::moveallow
    Type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::cllist
    integer :: b,pr,st,maxback,maxl
    logical :: clusmove

    clusmove = .true.
    moveallow = .true.
    maxl = chlen(chainno)
    do b =1,maxl
       tempcluscoord(chainno,b)%x = modulo(protcoords(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(protcoords(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(protcoords(chainno,b)%z + delz-1,gridsize)+1
       if(b<maxl) ctbdi(chainno,b) = bdi(chainno,b)

       do pr = 1,nprotein
          if(cl(pr) .eqv. .true.) goto 39
          maxback = chlen(pr)
          do st = 1,maxback,1
             if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
                clusmove = .false.
                !write(6,*) 'translation fail'
                goto 31
             end if
          end do
39         continue
       end do

    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove

  subroutine debug
    integer::m,l,f,g,dx,dy,dz,maxlengthss,maxl,zz
    integer,dimension(:,:),allocatable:: dbbdi
    logical:: bondnoneigh
    allocate(dbbdi(nprotein,maxlength))

    bondnoneigh = .false.
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

       end do
    end do

   
    do m = 1,nprotein
       maxl = chlen(m)
       do l= 1,maxl-1

          if((protcoords(m,l+1)%x - protcoords(m,l)%x == 1) .or. (protcoords(m,l+1)%x -&
               protcoords(m,l)%x == -1*(gridsize-1))) then
             dbbdi(m,l) = 1 
          else if((protcoords(m,l+1)%x - protcoords(m,l)%x == -1) .or. (protcoords(m,l+1)%x -&
               protcoords(m,l)%x == (gridsize-1))) then
             dbbdi(m,l) = -1
          else if((protcoords(m,l+1)%y - protcoords(m,l)%y == 1) .or. (protcoords(m,l+1)%y - &
               protcoords(m,l)%y == -1*(gridsize-1)))then
             dbbdi(m,l) = 2
          else if((protcoords(m,l+1)%y - protcoords(m,l)%y == -1) .or. (protcoords(m,l+1)%y - &
               protcoords(m,l)%y == (gridsize-1))) then
             dbbdi(m,l) = -2
          else if((protcoords(m,l+1)%z - protcoords(m,l)%z == 1) .or. (protcoords(m,l+1)%z -&
               protcoords(m,l)%z == -1*(gridsize-1))) then
             dbbdi(m,l) = 3
          else if((protcoords(m,l+1)%z - protcoords(m,l)%z == -1) .or. (protcoords(m,l+1)%z - &
               protcoords(m,l)%z == (gridsize-1))) then
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
    !bond checker
    do m =1,nprotein
       maxl = chlen(m)
       do l=1,maxl
if(protcoords(m,l)%bstate ==3) bondnoneigh = .true.

          if(protcoords(m,l)%am /= 0) then
             g = protcoords(m,l)%am
             f = protcoords(m,l)%al
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'a'
                finalfail = .true.
             end if
             if((protcoords(g,f)%bm /=m) .or. (protcoords(g,f)%bl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'a',protcoords(g,f)%bm,&
                     protcoords(g,f)%bl,removethis,pivde
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
             if((protcoords(m,l)%am == protcoords(m,l)%bondg) .and. (protcoords(m,l)%al &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
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
                     protcoords(g,f)%al,removethis,pivde
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
             if((protcoords(m,l)%bm == protcoords(m,l)%bondg) .and. (protcoords(m,l)%bl &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
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
                     protcoords(g,f)%dl,removethis,pivde
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
             if((protcoords(m,l)%cm == protcoords(m,l)%bondg) .and. (protcoords(m,l)%cl &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
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
                     protcoords(g,f)%cl,removethis,pivde
                finalfail = .true.
             end if
             if((protcoords(m,l)%dm == protcoords(m,l)%bondg) .and. (protcoords(m,l)%dl &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
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
                     protcoords(g,f)%fl,removethis,pivde
                                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
             if((protcoords(m,l)%em == protcoords(m,l)%bondg) .and. (protcoords(m,l)%el &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
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
                     protcoords(g,f)%el,removethis,pivde
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
             if((protcoords(m,l)%fm == protcoords(m,l)%bondg) .and. (protcoords(m,l)%fl &
                  ==protcoords(m,l)%bondf) .and. (bondnoneigh .eqv. .true.)) bondnoneigh = .false.
          end if

          if(bondnoneigh .eqv. .true.) then
             write(6,*) 'FAIL: bonded to bead that is not nearest neighbour'
             finalfail =.true.
             end if

if(protcoords(m,l)%bondg /= 0) then
             g = protcoords(m,l)%bondg
             f = protcoords(m,l)%bondf
             if(protcoords(m,l)%bstate /=3) write(6,*) 'fail: non-zerobondg value for nonbound site'
             if(protcoords(g,f)%bstate /=3) write(6,*) 'fail: bonded bead is unaware'
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -incorrect bonding configuration-bonded to itself',m,l,g,f,'f'
                finalfail = .true.
             end if
             if((protcoords(g,f)%bondg /=m) .or. (protcoords(g,f)%bondf /=l)) then
                write(6,*) 'FAIL-incorrect bonding configuration -bond mismatch',m,l,g,f,'f',protcoords(g,f)%bondg,&
                     protcoords(g,f)%bondf
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

  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize
    integer,dimension(:,:),allocatable::cb
    allocate(cb(nprotein,nprotein))
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

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do


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
                         cb(m,g) = 1
                         cb(g,m) = 1
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
                      if((protcoords(g,f)%bondg == m) .and. (protcoords(g,f)%bondf==l) .and. &
                           (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0)) then
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


    !maxclus = maxval(histcl)


    maxclus = 0
    do g = 1,nprotein
       if(histcl(g) > maxclus) then
          maxclus = histcl(g)
          content = g
       end if
    end do

    if((maxclus >10) .and. (suffclust .eqv. .false.)) then
       suffclust = .true.
       call pdbsetupalt(content,clnos)
    end if
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
    if(clustcount /= 0) write(82,*) time,clustcount, clusterpop, real(clusterpop)/clustcount &
         , real((normalisedsize) +(nprotein-clusterpop))/nprotein,maxclus !,clcount
    !if(time>10000) write(6,*) 'output',time, real((normalisedsize)+(nprotein-clusterpop))/nprotein,real((normalisedsize)+ &
    !(nprotein-clusterpop))
    if(clustcount == 0)  write(82,*) time,clustcount, clusterpop, 0.0
    !end if

    call phase(content,maxclus,clnos,cb)
    deallocate(cllist)
    deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
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
          write(67,*) 'atom', acount, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species), &
               'resid ',protcoords(m,1)%species
          write(67,*) 'atom', acount+1, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species)-1,'resid ', &
               protcoords(m,1)%species
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
               'resid ',infoout !protcoords(m,1)%species
          write(88,*) 'atom', acount+1, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species)-1,'resid ', &
               infoout !protcoords(m,1)%species
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
          write(88,*) 'bond', adjustr(atnum)//':',acount+2
          write(88,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+2
       end do
       write(atnum,'(i7)') acount
       write(88,*) 'bond', adjustr(atnum)//':',acount+1
       acount = acount +2
    end do

  end subroutine pdbsetupalt

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
    double precision:: intra


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
       CASE ('SCALEINFO')
          call get_logical(scalinginfo)

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
          allocate(intraenergy(mtype,mtype))
          !allocate(tttt())

       CASE ('SIM_TIME')
          CALL get_integer(maxtime)
          maxtime = maxtime*1000
       CASE ('EQUIB_TIME')
          CALL get_integer(equilib)
          equilib = equilib*1000
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
          call get_dp(intra)
          write(6,*) 'mtype',mtype
          do f = 1,mtype
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
