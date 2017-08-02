program move
  implicit none

  double precision :: intraenergy,interenergy,totalenergy,kT,totrmsbrute,randomz
  double precision :: totdisp,totchainlength,avechainlength,probs
  double precision :: actualchainlength,random
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug
  integer :: piv,right,end,crank,rept,datayes,maxclussize
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib
  integer :: successful,reject
  character(len = 10) :: commandread,commandread2,comread3,comread4,comread5
  character(len = 10) ::  comread6,comread7,comread8,comread9,comread10,comread11,comread12
  integer,dimension(:,:,:), allocatable :: isobond
  integer,dimension(:,:,:,:),allocatable :: intbond
  real, external :: ran2
  real,dimension(:),allocatable :: variance,disp
  logical :: exist,fail,finalfail
  real::start,finish


  !To do list
  !check energy counting routine
  !check energy running count


  type protein
     integer :: x,y,z
  end type protein

  type centremass
     double precision :: x,y,z
  end type centremass

  type grouping
     integer,dimension(1000) :: clusterlist
     integer :: clustersize
  end type grouping

  !  type(grouping),dimension(:),allocatable :: cluster

  type(protein),dimension(:,:),allocatable :: protcoords
  type(centremass),dimension(:,:),allocatable :: com
  type(protein) :: comcluster


  open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')
  open(29, file = 'rms.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  open(79, file = 'runningave.dat', action = 'write')
  open(93, file = 'energy.dat', action = 'write')
  open(82,file = 'clusterdata.dat', action = 'write')
  open(13,file='timedata.dat',action = 'write')
  reptforward = 0
  reptbackward = 0
  totrmsbrute = 0.0d0
  totalenergy = 0.0d0
  runningaveROG = 0.0d0
  runningaveEtE = 0.0d0
  finalfail = .false.
  maxclussize = 1

  !kT = 50.0
  successful = 0
  reject = 0
  call initial
  call get_command_argument(1,commandread)
  read(commandread, '(i3)') seed
  call get_command_argument(2,commandread2)
  read(commandread2, '(i4)') maxlength
  call get_command_argument(3,comread3)
  read(comread3, '(i1)') right
  call get_command_argument(4,comread4)
  read(comread4, '(i1)') crank
  call get_command_argument(5,comread5)
  read(comread5, '(i1)') end
  call get_command_argument(6,comread6)
  read(comread6, '(i1)') rept
  call get_command_argument(7,comread7)
  read(comread7, '(i1)') piv
  call get_command_argument(8,comread8)
  read(comread8, '(i7)') maxtime
  call get_command_argument(9,comread9)
  read(comread9, '(i5)') equilib 
  call get_command_argument(10,comread10)
  read(comread10, '(i1)') debugging
  call get_command_argument(11,comread11)
  read(comread11, '(i1)') datayes
  call get_command_argument(12,comread12)
  read(comread12,*) kT

  write(6,*) 'kT =', kT

  if (datayes == 1) open(67, file = 'move.xyz', action = 'write')
  allocate(protcoords(nprotein,maxlength))
  allocate(com(maxtime+1,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(isobond(nprotein,maxlength,maxlength))
  allocate(intbond(nprotein,nprotein,maxlength,maxlength))
  !allocate(cluster(nprotein))
  !allocate(comcluster(nprotein))
  N = maxlength
  interenergy =  -1.0 !-1.0 !-1.0d0
  intraenergy = -1.0  !-1.0 !-1.0d0

  count = 0
  time = 0
  totdisp = 0.0d0

  call foundation



  !call dataout
  call comfind

  !call clustercom
  actualchainlength = 0.0
  do time = 1,maxtime,1
     call CPU_TIME(start)
     fail = .false.
     if (mod(time,100) == 0 .and. datayes == 1) then
        call dataout
        call bondcheck
     end if
     call positioning

     call comfind
     !call clustercom


     call clusterposition

     if(modulo(time,100) == 0) call clustercount
    !write(6,*) 'reintegration'
     if(modulo(time,10) == 0) write(93,*) time,totalenergy
     if (time > equilib) then
        call rms
     end if
  !write(6,*) 'maxclus', maxclussize
     if(time >equilib .and. modulo(time,100) == 0) then
        call length
        call radiusofgyration
        write(79,*) time-equilib, runningaveEtE/(time-equilib), runningaveROG/(time-equilib)
     end if
     if (debugging == 1 .and. modulo(time,1000) == 0) then
        call debug
     else
        continue
     end if

     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else if(modulo(time,100) == 0) then
        write(6,*) 'step =',time
     end if
     call CPU_TIME(finish)
     if(modulo(time,100) == 0) write(13,*) (finish-start)*100
     
  end do
  !call error
  call energy
  write(6,*) 'average chain length squared is;' ,(actualchainlength/(maxtime/10)), &
       sqrt(actualchainlength/(maxtime/10))
  write(6,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))
  inquire(file = "lengthdataaverages.dat", exist = exist)
  if (exist) then
     open(19, file = "lengthdataaverages.dat", status = "old", position = "append",action ="write")
  else
     open(19, file = "lengthdataaverages.dat", status = "new",action ="write")
  end if
  write(19,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))!,2*(SQRT(totvar)/maxtime)/(actualchainlength/maxtime)
  write(6,*) 'reject = ', reject
  write(6,*) 'successful =', successful
  write(6,*) 1.0*successful/(reject + successful)
  write(6,*) 'acceptance ratio =', 1.0*successful/(reject + successful)
  write(6,*) 'count =', count
  write(6,*) 'maxlength =', maxlength
  write(6,*) 'forward =', reptforward
  write(6,*) 'backward =' , reptbackward

  if (finalfail .eqv. .true.) then
     write(6,*) 'This simulation FAILED - OVERLAPS were observed'
  end if
contains

  subroutine initial
    !reads in initial parameters from setup test file

    character(len = 10) :: BIN
    read(17,*) BIN
    read(17,*) BIN !seed
    read(17,*) BIN
    read(17,*) BIN !maxlength
    read(17,*) BIN
    read(17,*) nprotein
    read(17,*) BIN
    read(17,*) gridsize
    natoms = gridsize**3
    read(17,*) BIN
    read(17,*) BIN !maxtime
  end subroutine initial

  subroutine foundation
    !read in initial positions of the protein
    integer::m,l
    character(len = 10) :: BIN
    read(23,*) BIN
    write(6,*) maxlength
    !write(67,*) nprotein*maxlength
    !write(67,*) ' '
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          read(23,*) BIN, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
       end do
    end do

    call energy
    call debug
  end subroutine foundation

  subroutine clusterposition
    integer::scans,m
    real :: decider
    type(grouping) :: cluster
    m =int(ran2(seed)*(nprotein))+1
    call clusterassign(m,cluster)
if(1/cluster%clustersize > ran2(seed))then
decider = ran2(seed) - 0.5
    do scans = 1,nprotein
       if(decider > 0.0) then
          !write(6,*) 'move'
          call clustertranslation(m,cluster)
       else if(decider <= 0.0) then
          !write(6,*) 'turn'
          call clusterrotation(m,cluster)
       end if
    end do
    end if
  end subroutine clusterposition

  subroutine positioning
    !Selects move to be carried out on the proteins
    integer :: nmoves,scan
    logical :: run,run2

    !write(6,*) 'a'
    t = time + 1
    do scan = 1,nprotein*maxlength,1
       count = count + 1
       run = .true.
       run2 = .true. 

       randomz = ran2(seed)
       nmoves = piv + crank + end + right + rept

       if (randomz <= (1.0*end)/ nmoves .and. end == 1) then
          !write(6,*) 'end' 
          call endmove
       else if (randomz > (1.0*end)/nmoves .and. randomz <= (1.0*end+crank)/nmoves .and. crank == 1) then
          !write(6,*) 'crank' 
          call crankshaftmove
       else if (randomz > (1.0*end + crank)/nmoves .and. randomz <= (1.0*end + crank + right)/nmoves .and. right ==1) then
          !write(6,*) 'right' 
          call rightanglemove
       else if (randomz > (1.0*end + crank + right)/nmoves .and. randomz <= (1.0*end + crank + right + rept)/nmoves &
            .and. rept ==1) then
          !write(6,*) 'reptation' 
          call reptation
       else if (randomz > (1.0*end + crank + right+rept)/nmoves .and. randomz <= (1.0*end + crank + right + rept + piv)/nmoves &
            .and. piv == 1) then 
          !write(6,*) 'pivot' 
          call pivot
       end if
    end do
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
    integer :: dummy,dummy2,p,s,m,l,str,st,pr
    integer,dimension(:),allocatable::dx1,dy1,dz1
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(dx1(maxlength))
    allocate(dy1(maxlength))
    allocate(dz1(maxlength))

    !this means that m/=0 and is up to nprotein
    m =int(ran2(seed)*(nprotein))+1
    !this means that l cannot be maxlength,maxlength - 1 or l
    cranksep = .true.
71  if(cranksep .eqv. .true.) then
       continue
    end if
    !following section ensures p is at least 3 beads away from l and that p is greater than l
    dummy = int(ran2(seed)*(maxlength-1))+1
    dummy2 = int(ran2(seed)*(maxlength-1))+1
    l = min(dummy,dummy2)
    p = max(dummy,dummy2)
    if (p-l < 3) then
       goto 71
    end if

    probs = ran2(seed)
    pivotx = .false.
    pivoty = .false.
    pivotz = .false.

    !check to see if 4 atoms are in the same plane
    if (protcoords(m,l)%x /= protcoords(m,p)%x &
         .and.  protcoords(m,l)%y == protcoords(m,p)%y &
         .and.protcoords(m,l)%z == protcoords(m,p)%z) then
       pivotx = .true.                  
    else if (protcoords(m,l)%x == protcoords(m,p)%x .and. &
         protcoords(m,l)%y /= protcoords(m,p)%y &
         .and.protcoords(m,l)%z == protcoords(m,p)%z ) then
       pivoty = .true.                               
    else if (protcoords(m,l)%x == protcoords(m,p)%x &
         .and. protcoords(m,l)%y == protcoords(m,p)%y &
         .and.protcoords(m,l)%z /= protcoords(m,p)%z ) then
       pivotz = .true.
    end if


    if (pivotx .eqv. .false. .and. pivoty .eqv. .false. .and. pivotz .eqv. .false.) then
       crankcont = .false.
       goto 43
    end if

    crankcont = .true.

    if (pivotx .eqv. .true.) then

       do s= l+1,p-1,1

          dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if(dy1(s) > gridsize/2) dy1(s) = dy1(s)-gridsize
          if(dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize


          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if(dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if(dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize

       end do

       if (probs <= (1.0/3)) then
          do s = l+1,p-1,1
             tempcoord(s)%x = protcoords(m,s)%x
             tempcoord(s)%y = modulo(protcoords(m,l)%y-dz1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z+dy1(s)-1,gridsize)+1
          end do

       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(s)%x = protcoords(m,s)%x
             tempcoord(s)%y = modulo(protcoords(m,l)%y+dz1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z-dy1(s)-1,gridsize)+1
          end do

       else if (probs > (2.0/3) .and. probs <= (1.0)) then

          do s = l+1,p-1,1
             tempcoord(s)%x = protcoords(m,s)%x
             tempcoord(s)%y = modulo(protcoords(m,l)%y-dy1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z-dz1(s)-1,gridsize)+1
          end do

       end if

    else if (pivoty .eqv. .true.) then

       do s= l+1,p-1,1

          dx1(s) = protcoords(m,s)%x - protcoords(m,l)%x
          if (dx1(s) > gridsize/2) dx1(s) = dx1(s) - gridsize
          if (dx1(s) < -gridsize/2) dx1(s) = dx1(s) + gridsize
          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if (dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if (dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize

       end do

       if (probs <= (1.0/3)) then
          do s = l+1,p-1,1
             tempcoord(s)%y = protcoords(m,s)%y
             tempcoord(s)%x = modulo(protcoords(m,l)%x-dz1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z+dx1(s)-1,gridsize)+1
          end do

       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(s)%y = protcoords(m,s)%y
             tempcoord(s)%x = modulo(protcoords(m,l)%x+dz1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z-dx1(s)-1,gridsize)+1
          end do

       else if (probs > (2.0/3) .and. probs <= (1.0)) then

          do s = l+1,p-1,1
             tempcoord(s)%y = protcoords(m,s)%y
             tempcoord(s)%x = modulo(protcoords(m,l)%x-dx1(s)-1,gridsize)+1
             tempcoord(s)%z = modulo(protcoords(m,l)%z-dz1(s)-1,gridsize)+1
          end do


       end if
       !!pivot y!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    else if (pivotz .eqv. .true.) then

       do s= l+1,p-1,1

          dx1(s) = protcoords(m,s)%x - protcoords(m,l)%x
          if (dx1(s) > gridsize/2) dx1(s) = dx1(s) - gridsize
          if (dx1(s) < -gridsize/2) dx1(s) = dx1(s) + gridsize
          dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if (dy1(s) > gridsize/2) dy1(s) = dy1(s) - gridsize
          if (dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize
       end do

       if (probs <= (1.0/3)) then
          do s = l+1,p-1,1
             tempcoord(s)%z = protcoords(m,s)%z
             tempcoord(s)%y = modulo(protcoords(m,l)%y-dx1(s)-1,gridsize)+1
             tempcoord(s)%x = modulo(protcoords(m,l)%x+dy1(s)-1,gridsize)+1
          end do

       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(s)%z = protcoords(m,s)%z
             tempcoord(s)%y = modulo(protcoords(m,l)%y+dx1(s)-1,gridsize)+1
             tempcoord(s)%x = modulo(protcoords(m,l)%x-dy1(s)-1,gridsize)+1
          end do

       else if (probs > (2.0/3) .and. probs <= (1.0)) then
          do s = l+1,p-1,1
             tempcoord(s)%y = modulo(protcoords(m,l)%y-dy1(s)-1,gridsize)+1
             tempcoord(s)%x = modulo(protcoords(m,l)%x-dx1(s)-1,gridsize)+1
             tempcoord(s)%z = protcoords(m,s)%z
          end do
       end if
    end if


    !if (pivotx .eqv. .false. .and. pivoty .eqv. .false. .and. pivotz .eqv. .false.) then
    !crankcont = .false.
    !end if

    deltaenergy = 0.0d0    
    !pr = m

    do str = l+1,p-1,1
       do st = 1,maxlength,1
          if(st < l+1 .or. st > p-1) then
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
             do st = 1,maxlength,1
                if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                   crankcont = .false.
                   goto 43
                end if
             end do
          end do
       end if
    end do


    call intrabondallocateupper(m,l+1,p-1,1,l,tempcoord,tempisobond,deltaenergy)
    call intrabondallocatelower(m,l+1,p-1,p,maxlength,tempcoord,tempisobond,deltaenergy)

    call interbondallocate(m,l+1,p-1,tempcoord,tempintbond,deltaenergy)

    if(Energydecision(deltaenergy) .eqv. .True.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l+1,p-1,tempcoord)
       call updateintrabond(m,l+1,p-1,1,l,tempisobond)
       call updateintrabond(m,l+1,p-1,p,maxlength,tempisobond)
       call updateinterbond(m,l+1,p-1,tempintbond)
    else 
       crankcont = .false.
       goto 43
    end if

    !       write(6,*) time, 'crank'
    successful = successful + 1

43  if (crankcont .eqv. .false.) then
       reject = reject + 1
    end if

  end subroutine crankshaftmove

  subroutine rightanglemove
    !performs a diagonal flip on a bead that has its two connecting beads perpendicular to each other
    integer :: m,l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz
    double precision :: deltaenergy
    logical :: rac
    integer:: st,pr
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))


    m =int(ran2(seed)*(nprotein))+1
    ! l cannot be equal to 1 or maxlength
    l = int(ran2(seed)*(maxlength-2))+2
    !write(6,*) 'm =', m, 'l = ', l

    rac = .true.

    if (abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= gridsize-2.0 &
         .and. abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= gridsize-2.0 &
         .and. abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= gridsize-2.0) then
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

       tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
       tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
       tempcoord(l)%z = modulo(protcoords(m,l)%z+dz-1,gridsize)+1

       deltaenergy = 0.0d0        
       do st = 1,maxlength,1
          if(st /=l) then
             if(overlaps(m,l,st,tempcoord) .eqv. .false.) then
                rac = .false.
                goto 37
             end if
          end if
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
             do st = 1,maxlength,1
                if(overlaps(pr,l,st,tempcoord) .eqv. .false.) then
                   rac = .false.
                   goto 37
                end if
             end do
          end if
       end do

       call intrabondallocateupper(m,l,l,1,l-3,tempcoord,tempisobond,deltaenergy)
       call intrabondallocatelower(m,l,l,l+3,maxlength,tempcoord,tempisobond,deltaenergy)

       call interbondallocate(m,l,l,tempcoord,tempintbond,deltaenergy)


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l,l,tempcoord)
          call updateintrabond(m,l,l,1,l-3,tempisobond)
          call updateintrabond(m,l,l,l+3,maxlength,tempisobond)
          call updateinterbond(m,l,l,tempintbond)
          successful = successful + 1
       else
          rac = .false.
          goto 37
       end if
    end if

37  if (rac .eqv. .false.) then
       reject = reject + 1
    end if

  end subroutine rightanglemove

  subroutine pivot
    !rotates a section of chain around the selected bead
    integer :: m,l,b,g,xhold,yhold,zhold,str,st,pr
    logical :: pivcont
    double precision :: choose3,deltaenergy
    integer,dimension(:),allocatable :: delx,dely,delz
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(delx(maxlength))
    allocate(dely(maxlength))   
    allocate(delz(maxlength))

    choose3 = ran2(seed)
    pivcont = .true. 
    m =int(ran2(seed)*(nprotein))+1
    l = int(ran2(seed)*(maxlength-2))+2
    t = time + 1

    if(l > maxlength/2) then
       do g = l+1,maxlength,1
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
       if (choose3 <= 1.0/5) then

          do b = l+1,maxlength,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x - dely(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y + delx(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z
          end do

       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then
          do b =l+1,maxlength,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x + dely(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y - delx(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z                    
          end do

       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then

          do b =l+1,maxlength,1                 
             tempcoord(b)%x = modulo(protcoords(m,l)%x - delx(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y - dely(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z
          end do

       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then

          do b =l+1,maxlength,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x - delz(b)-1,gridsize)+1
             tempcoord(b)%z = modulo(protcoords(m,l)%z + delx(b)-1,gridsize)+1
             tempcoord(b)%y = protcoords(m,b)%y    
          end do

       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then

          do b = l+1,maxlength,1                
             tempcoord(b)%x = modulo(protcoords(m,l)%x + delz(b)-1,gridsize)+1
             tempcoord(b)%z = modulo(protcoords(m,l)%z - delx(b)-1,gridsize)+1
             tempcoord(b)%y = protcoords(m,b)%y
          end do

       end if

       do str = l+1,maxlength,1
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
             do str = l+1,maxlength,1
                do st = 1,maxlength,1
                   if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end if
       end do


       call intrabondallocateupper(m,l+1,maxlength,1,l,tempcoord,tempisobond,deltaenergy)

       call interbondallocate(m,l+1,maxlength,tempcoord,tempintbond,deltaenergy)

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l+1,maxlength,tempcoord)
          call updateintrabond(m,l+1,maxlength,1,l,tempisobond)
          call updateinterbond(m,l+1,maxlength,tempintbond)
          successful = successful + 1

       else
          pivcont = .false.
          goto 75
       end if

    else if (l <= maxlength/2) then
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

       if (choose3 <= 1.0/5) then

          do b = 1,l-1,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x -dely(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y + delx(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z
          end do

       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then 
          do b =1,l-1,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x + dely(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y - delx(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z
          end do

       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then
          do b =1,l-1,1                 
             tempcoord(b)%x = modulo(protcoords(m,l)%x -delx(b)-1,gridsize)+1
             tempcoord(b)%y = modulo(protcoords(m,l)%y - dely(b)-1,gridsize)+1
             tempcoord(b)%z = protcoords(m,b)%z
          end do

       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then
          do b =1,l-1,1 
             tempcoord(b)%x = modulo(protcoords(m,l)%x - delz(b)-1,gridsize)+1
             tempcoord(b)%z = modulo(protcoords(m,l)%z + delx(b)-1,gridsize)+1
             tempcoord(b)%y = protcoords(m,b)%y
          end do

       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then
          do b = 1,l-1,1               
             tempcoord(b)%x = modulo(protcoords(m,l)%x + delz(b)-1,gridsize)+1
             tempcoord(b)%z = modulo(protcoords(m,l)%z - delx(b)-1,gridsize)+1
             tempcoord(b)%y = protcoords(m,b)%y
          end do
       end if
       deltaenergy = 0.0d0
       do str = 1,l-1,1
          do st = l,maxlength,1
             if(overlaps(m,str,st,tempcoord) .eqv. .false.) then
                pivcont = .false.
                goto 75
             end if
          end do
       end do

       do pr = 1,nprotein,1
          if(pr /= m) then
             do str = 1,l-1,1
                do st = 1,maxlength,1
                   if(overlaps(pr,str,st,tempcoord) .eqv. .false.) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end if
       end do

       call intrabondallocatelower(m,1,l-1,l,maxlength,tempcoord,tempisobond,deltaenergy)


       call interbondallocate(m,1,l-1,tempcoord,tempintbond,deltaenergy)

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,1,l-1,tempcoord)
          call updateintrabond(m,1,l-1,l,maxlength,tempisobond)
          call updateinterbond(m,1,l-1,tempintbond)
          successful = successful + 1
       else
          pivcont = .false.
          goto 75
       end if
    end if

75  if (pivcont .eqv. .false.) then
       reject = reject + 1
    end if

  end subroutine pivot


  subroutine intrabondallocatelower(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,deltaenergy)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer::bead1,bead2
    Type(protein),dimension(:),intent(in) :: tempcoord

    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead2 > bead1 + 2) then
             if(isobond(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraenergy
             if(adjacent(chain1,bead1,bead2,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + intraenergy
                tempisobond(bead1,bead2) = 1
             else 
                tempisobond(bead1,bead2) = 0
             end if
          end if
       end do
    end do

  end subroutine intrabondallocatelower

  subroutine intrabondallocateupper(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,deltaenergy)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer::bead1,bead2
    Type(protein),dimension(:),intent(in) :: tempcoord

    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead1 > bead2 + 2) then
             if(isobond(chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - intraenergy
             if(adjacent(chain1,bead1,bead2,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + intraenergy
                tempisobond(bead2,bead1) = 1
             else 
                tempisobond(bead2,bead1) = 0
             end if
          end if
       end do
    end do

  end subroutine intrabondallocateupper

  subroutine interbondallocate(chain1,beadmin,beadmax,tempcoord,tempintbond,deltaenergy)

    integer,intent(in)::chain1,beadmin,beadmax
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:,:),intent(inout) :: tempintbond
    integer::bead1,bead2,chain2
    Type(protein),dimension(:),intent(in) :: tempcoord

    do chain2 = 1,chain1-1,1
       do bead1 = beadmin,beadmax,1
          do bead2 = 1,maxlength,1
             if(intbond(chain2,chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(chain2,bead1,bead2,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(chain2,bead1,bead2) = 1
             else
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do

    do chain2 = chain1+1,nprotein,1
       do bead1 = beadmin,beadmax,1
          do bead2 = 1,maxlength,1
             if(intbond(chain1,chain2,bead1,bead2) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(chain2,bead1,bead2,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(chain2,bead1,bead2) = 1
             else
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do
  end subroutine interbondallocate

  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax
    Type(protein),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding

    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
    end do

  end subroutine updatepos

  subroutine updatecluspos(cluster,tempcluscoord)
    !updates bead positions
    type(grouping),intent(inout) ::cluster
    Type(protein),dimension(:,:),intent(in) :: tempcluscoord
    integer ::b,l,a
    !moves beads to new positions and reassigns isobonding

    do a = 1,cluster%clustersize
       b = cluster%clusterlist(a)
       do l = 1,maxlength,1
          protcoords(b,l)%x = tempcluscoord(b,l)%x
          protcoords(b,l)%y = tempcluscoord(b,l)%y
          protcoords(b,l)%z = tempcluscoord(b,l)%z
       end do
    end do
  end subroutine updatecluspos

  subroutine bondcheck

    integer:: m,l,f
    open(12, file = 'bondintra.dat', action = 'write')
    do m = 1,nprotein
       do l = 1,maxlength,1
          do f = 1,maxlength,1
             if (f < l -2 .or. f > l+ 2) write(12,*) l,f,isobond(m,l,f)
          end do
       end do
    end do

  end subroutine bondcheck

  subroutine updateintrabond(chainnum,beadmin,beadmax,statmin,statmax,tempisobond)
    integer,intent(in):: chainnum,beadmin,beadmax,statmin,statmax
    integer ::g,beadnum
    integer,dimension(:,:),intent(in) :: tempisobond
    do beadnum = beadmin,beadmax,1
       do g = statmin,statmax,1
          if(g<beadnum-2) isobond(chainnum,g,beadnum) = tempisobond(g,beadnum) 
          if(g> beadnum+2) isobond(chainnum,beadnum,g) = tempisobond(beadnum,g)
       end do
    end do
  end subroutine updateintrabond


  subroutine clusterassign(m,cluster)
    integer,intent(in):: m
    type(grouping),intent(inout) ::cluster
    type(protein),dimension(:),allocatable :: tempcoord
    integer :: a,n,controlcount,l,g,counter,chainnum
    logical :: newcluster,clusterpass 
    allocate(tempcoord(maxlength))
    controlcount = 1
    cluster%clustersize = 1
    cluster%clusterlist(1) = m
    do a = 2,nprotein
cluster%clusterlist(a) = 0
       end do

    clusterpass = .true.
    counter = 0
    do while (clusterpass .eqv. .true.)
       counter = counter + 1
        chainnum = cluster%clusterlist(counter)
        clusterpass = .false.
       do l = 1,maxlength
          tempcoord(l)%x = protcoords(chainnum,l)%x
          tempcoord(l)%y = protcoords(chainnum,l)%y
          tempcoord(l)%z = protcoords(chainnum,l)%z
          do n = 1,nprotein
             newcluster = .false.
             if (ANY(cluster%clusterlist(:) .eq. n)) then
                goto 25
             else
                do g= 1,maxlength
                   if(adjacent(n,l,g,tempcoord) .eqv. .true.) then
                      newcluster = .true.
                      clusterpass = .true.
                      goto 56
                   end if
                end do
             end if
25           if(newcluster .eqv. .false. .and. cluster%clusterlist(counter) /= n) then
                if(counter < cluster%clustersize) clusterpass = .true.
             end if
          end do
56        if(newcluster .eqv. .true.) then
             controlcount = controlcount + 1
             cluster%clustersize = controlcount
             cluster%clusterlist(controlcount) = n
          end if
       end do
    end do
!write(6,*) 'a'
    if(controlcount > maxclussize) then
       maxclussize = controlcount
    end if
   !write(6,*) 'max',maxclussize
  end subroutine clusterassign


  subroutine updateinterbond(chainnum,beadmin,beadmax,tempintbond)
    integer,intent(in):: chainnum,beadmin,beadmax
    integer,dimension(:,:,:),intent(in) :: tempintbond
    integer ::g,beadnum,chain2
    do beadnum = beadmin,beadmax,1
       do chain2 = 1,chainnum-1,1
          do g = 1,maxlength,1
             intbond(chain2,chainnum,g,beadnum) = tempintbond(chain2,beadnum,g)
          end do
       end do
       do chain2 = chainnum+1,nprotein,1
          do g = 1,maxlength,1
             intbond(chainnum,chain2,beadnum,g) = tempintbond(chain2,beadnum,g)
          end do
       end do
    end do

  end subroutine updateinterbond

  subroutine reptation
    !perform reptation
    double precision :: choose,direc,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx
    integer:: st,pr,dx,dy,dz,chain2,beads2
    logical :: reptcont
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))

    t = time + 1
    choose =  ran2(seed) - 0.5
    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1
    g1 = 1
    g2 = 1
    g3 = 1
    g4 = 1
    g5 = 1
    g6 = 1
    dx = 0
    dy = 0
    dz = 0

    if (choose > 0.0) then
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

       direc = ran2(seed)
       if(direc <= (1.0*g1)/5 .and. g1 == 1) then
          dx = 1
       else if(direc <= (1.0*g1 +g2)/5 .and.direc > (1.0*g1)/5 .and.  g2 == 1) then
          dx = -1
       else if(direc <= (1.0*g1+g2+g3)/5 .and.direc > (1.0*g1+g2)/5 .and.  g3 == 1) then
          dy = 1
       else if(direc <= (1.0*g1+g2+g3+g4)/5 .and.direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
          dy = -1
       else if(direc <= (1.0*g1+g2+g3+g4+g5)/5 .and.direc > (1.0*g1+g2+g3+g4)/5 .and. &
            g5 == 1) then
          dz = 1
       else if(direc <= (1.0*g1+g2+g3+g4+g5+g6)/5 .and. &
            direc > (1.0*g1+g2+g3+g4+g5)/5 .and.  g6 == 1) then
          dz = -1
       end if

       tempcoord(1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(1)%z = modulo(protcoords(m,1)%z + dz-1,gridsize)+1

       deltaenergy = 0.0d0
       do pr =1,nprotein    
          do st = 1,maxlength,1
             if(pr /=m .or. st > 1 .or. st < maxlength)then
                if(overlaps(pr,1,st,tempcoord) .eqv. .false.) then
                   reptcont = .false.
                   goto 83
                end if
             end if
          end do
       end do

       do st = 1,maxlength-3,1
          if(isobond(m,st,maxlength) == 1) deltaenergy = deltaenergy - intraenergy
       end do


       do st = 3,maxlength-1
          if(adjacent(m,1,st,tempcoord) .eqv. .true.) then
             deltaenergy = deltaenergy + intraenergy
             tempisobond(1,st+1) = 1
          else
             tempisobond(1,st+1) = 0
          end if
       end do

       do pr = 1,m-1,1
          do st = 1,maxlength
             if(intbond(pr,m,st,maxlength) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(pr,1,st,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(pr,1,st) = 1
             else
                tempintbond(pr,1,st) = 0
             end if
          end do
       end do

       do pr = m+1,nprotein,1
          do st = 1,maxlength
             if(intbond(m,pr,maxlength,st) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(pr,1,st,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(pr,1,st) = 1
             else
                tempintbond(pr,1,st) = 0
             end if
          end do
       end do

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          do l =maxlength,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z


             do g = l-3,2,-1
                isobond(m,g,l) = isobond(m,g-1,l-1)
             end do

             do chain2 = 1,m-1,1                   
                do beads2 = 1,maxlength
                   intbond(chain2,m,beads2,l) = intbond(chain2,m,beads2,l-1)
                end do
             end do

             do chain2 = m+1,nprotein,1
                do beads2 = 1,maxlength,1
                   intbond(m,chain2,l,beads2) = intbond(m,chain2,l-1,beads2)
                end do
             end do
          end do

          call updatepos(m,1,1,tempcoord)
          call updateintrabond(m,1,1,4,maxlength,tempisobond)
          call updateinterbond(m,1,1,tempintbond)

          successful = successful + 1
          reptforward = reptforward + 1

       else
          reptcont = .false.
       end if
    else if (choose < 0.0) then

       dxx = protcoords(m,maxlength-1)%x - protcoords(m,maxlength)%x
       if (dxx > 1) dxx = -1
       if(dxx <-1) dxx = 1
       dyx = protcoords(m,maxlength-1)%y - protcoords(m,maxlength)%y
       if (dyx > 1) dyx = -1
       if (dyx < -1) dyx = 1
       dzx = protcoords(m,maxlength-1)%z - protcoords(m,maxlength)%z
       if (dzx > 1) dzx = -1
       if (dzx < -1) dzx = 1

       if (dxx == 1) g1 = 0
       if (dxx == -1) g2 = 0
       if (dyx == 1) g3 = 0
       if (dyx == -1) g4 = 0
       if (dzx == 1) g5 =0
       if (dzx == -1) g6 = 0

       direc = ran2(seed)
       if(direc <= (1.0*g1)/5 .and. g1 == 1) then
          dx = 1
       else if(direc <= (1.0*g1 +g2)/5 .and.direc > (1.0*g1)/5 .and.  g2 == 1) then
          dx = -1
       else if(direc <= (1.0*g1+g2+g3)/5 .and.direc > (1.0*g1+g2)/5 .and.  g3 == 1) then
          dy = 1
       else if(direc <= (1.0*g1+g2+g3+g4)/5 .and. direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
          dy = -1
       else if(direc <= (1.0*g1+g2+g3+g4+g5)/5 .and. &
            direc > (1.0*g1+g2+g3+g4)/5 .and.  g5 == 1) then
          dz = 1
       else if(direc <= (1.0*g1+g2+g3+g4+g5+g6)/5 .and. &
            direc > (1.0*g1+g2+g3+g4+g5)/5 .and.  g6 == 1) then
          dz = -1
       end if

       tempcoord(maxlength)%x = modulo(protcoords(m,maxlength)%x + dx-1,gridsize)+1
       tempcoord(maxlength)%y = modulo(protcoords(m,maxlength)%y + dy-1,gridsize)+1
       tempcoord(maxlength)%z = modulo(protcoords(m,maxlength)%z + dz-1,gridsize)+1

       deltaenergy = 0.0d0
       do pr =1,nprotein    
          do st = 1,maxlength,1
             if(pr /=m .or. st >1 .or. st < maxlength)then
                if(overlaps(pr,maxlength,st,tempcoord) .eqv. .false.) then
                   reptcont = .false.
                   goto 83
                end if
             end if
          end do
       end do

       do st = 4,maxlength,1
          if(isobond(m,1,st) == 1) deltaenergy = deltaenergy - intraenergy
       end do

       do st = 2,maxlength-2
          if(adjacent(m,maxlength,st,tempcoord) .eqv. .true.) then
             deltaenergy = deltaenergy + intraenergy
             tempisobond(st-1,maxlength) = 1
          else
             tempisobond(st-1,maxlength) = 0
          end if
       end do

       do pr = 1,m-1,1
          do st = 1,maxlength,1
             if(intbond(pr,m,st,1) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(pr,maxlength,st,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(pr,maxlength,st) = 1
             else
                tempintbond(pr,maxlength,st) = 0
             end if
          end do
       end do

       do pr = m+1,nprotein,1
          do st = 1,maxlength,1
             if(intbond(m,pr,1,st) == 1) deltaenergy = deltaenergy - interenergy
             if(adjacent(pr,maxlength,st,tempcoord) .eqv. .true.) then
                deltaenergy = deltaenergy + interenergy
                tempintbond(pr,maxlength,st) = 1
             else
                tempintbond(pr,maxlength,st) = 0
             end if
          end do
       end do

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          do l =1,maxlength-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z

             do g= l+3,maxlength-1,1
                isobond(m,l,g) = isobond(m,l+1,g+1)
             end do

             do chain2 = 1,m-1,1
                do beads2 = 1,maxlength
                   intbond(chain2,m,beads2,l) = intbond(chain2,m,beads2,l+1)
                end do
             end do

             do chain2 = m+1,nprotein,1
                do beads2 = 1,maxlength
                   intbond(m,chain2,l,beads2) = intbond(m,chain2,l+1,beads2)
                end do
             end do

          end do
          call updatepos(m,maxlength,maxlength,tempcoord)
          call updateintrabond(m,maxlength,maxlength,1,maxlength-3,tempisobond)
          call updateinterbond(m,maxlength,maxlength,tempintbond)
          successful = successful + 1
          reptbackward = reptbackward + 1
       else
          reptcont = .false.
       end if
    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
    end if

  end subroutine reptation

  subroutine endmove
    !performs end chain rotation
    double precision :: choose2,deltaenergy
    logical :: endcont
    integer :: l,m,st,pr
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))

    m = int(ran2(seed)*(nprotein))+1
    random = ran2(seed)
    choose2 = ran2(seed)-0.5
    endcont = .true.

    if (choose2 >= 0.0) then
       l = 1
       tempcoord(l)%x = protcoords(m,l+1)%x
       tempcoord(l)%y = protcoords(m,l+1)%y
       tempcoord(l)%z = protcoords(m,l+1)%z
    else if (choose2 < 0.0) then
       l = maxlength
       tempcoord(l)%x = protcoords(m,l-1)%x
       tempcoord(l)%y = protcoords(m,l-1)%y
       tempcoord(l)%z = protcoords(m,l-1)%z
    end if

    if(random<= 1.0/6) then
       tempcoord(l)%x = modulo(tempcoord(l)%x,gridsize)+1
    else if (random> 1.0/6 .AND. random<=2.0/6) then
       tempcoord(l)%x = modulo(tempcoord(l)%x - 2,gridsize)+1
    else if (random> 2.0/6 .AND. random<=1.0/2) then
       tempcoord(l)%y = modulo(tempcoord(l)%y,gridsize)+1
    else if (random> 1.0/2 .AND. random<=2.0/3) then                
       tempcoord(l)%y = modulo(tempcoord(l)%y -2,gridsize)+1
    else if (random> 2.0/3 .AND. random<=5.0/6) then
       tempcoord(l)%z  = modulo(tempcoord(l)%z,gridsize)+1
    else if (random> 5.0/6 .AND. random<=1.0) then                
       tempcoord(l)%z = modulo(tempcoord(l)%z-2,gridsize)+1
    end if

    deltaenergy = 0.0d0
    do pr =1,nprotein    
       do st = 1,maxlength,1
          if(pr /=m .or. st >1 .or. st < maxlength)then
             if(overlaps(pr,l,st,tempcoord) .eqv. .false.) then
                endcont = .false.
                goto 71
             end if
          end if
       end do
    end do

    if(l == 1)  call intrabondallocatelower(m,1,1,4,maxlength,tempcoord,tempisobond,deltaenergy)
    if(l==maxlength) call intrabondallocateupper(m,maxlength,maxlength,1,maxlength-3,tempcoord,tempisobond,deltaenergy)

    call interbondallocate(m,l,l,tempcoord,tempintbond,deltaenergy)

    if(Energydecision(deltaenergy) .eqv. .true.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l,l,tempcoord)
       if(l==1) call updateintrabond(m,1,1,4,maxlength,tempisobond)
       if(l==maxlength) call updateintrabond(m,maxlength,maxlength,1,maxlength-3,tempisobond) 
       call updateinterbond(m,l,l,tempintbond)
       successful = successful + 1
    else
       endcont = .false.
    end if

71  if (endcont .eqv. .false.) then
       reject = reject + 1
    end if
  end subroutine endmove

  subroutine comfind
    integer :: dcomx,dcomy,dcomz
    double precision :: comx,comy,comz
    integer :: m,l
    t = time +1
    do m = 1,nprotein,1
       comx = 0.0d0
       comy = 0.0d0
       comz = 0.0d0
       do l = 2,maxlength,1
          dcomx = protcoords(m,l)%x-protcoords(m,l-1)%x
          if (dcomx > 1.5) dcomx = -1
          if (dcomx < -1.5) dcomx = 1
          comx = comx + dcomx
          dcomy = protcoords(m,l)%y-protcoords(m,l-1)%y
          if (dcomy > 1.5) dcomy = -1
          if (dcomy < -1.5) dcomy = 1
          comy = comy + dcomy
          dcomz = protcoords(m,l)%z-protcoords(m,l-1)%z
          if (dcomz > 1.5) dcomz = -1
          if (dcomz < -1.5) dcomz = 1
          comz = comz + dcomz
       end do
       com(t,m)%x = modulo(protcoords(m,1)%x +(comx/2),real(gridsize))
       com(t,m)%y = modulo(protcoords(m,1)%y +(comy/2),real(gridsize))
       com(t,m)%z = modulo(protcoords(m,1)%z +(comz/2),real(gridsize))
    end do
  end subroutine comfind

  subroutine clustercom(m,cluster)
    integer,intent(in)::m
    type(grouping),intent(inout) ::cluster
    integer ::  a
    double precision :: delx,dely,delz,comx,comy,comz
    logical :: compass
    t = time+1
!write(6,*) 'com start'
 compass = .false.
       comcluster%x = 0
       comcluster%y = 0
       comcluster%z = 0
       if(cluster%clustersize ==1) then

       comcluster%x = int(com(t,m)%x)
       comcluster%y = int(com(t,m)%y)
       comcluster%z = int(com(t,m)%z)
       compass = .true.
       goto 53
end if
          
       do a = 1,cluster%clustersize-1
          delx = com(t,cluster%clusterlist(a+1))%x - com(t,cluster%clusterlist(a))%x
          if (delx > gridsize/2) delx = delx-gridsize
          if (delx < -1.5) delx = delx + gridsize
          comx = comx + delx
          dely = com(t,cluster%clusterlist(a+1))%y -com(t,cluster%clusterlist(a))%y
          if (dely > gridsize/2) dely = dely-gridsize
          if (dely < -1.5) dely = dely + gridsize
          comy = comy + dely
          delz = com(t,cluster%clusterlist(a+1))%z-com(t,cluster%clusterlist(a))%z
          if (delz > gridsize/2) delz = delz-gridsize
          if (delz < -1.5) delz = delz + gridsize
          comz = comz + delz
       end do
       comcluster%x = int(modulo(com(t,cluster%clusterlist(1))%x + delx/2 -1,real(gridsize)))+1
       comcluster%y = int(modulo(com(t,cluster%clusterlist(1))%y + dely/2 -1,real(gridsize)))+1
       comcluster%z = int(modulo(com(t,cluster%clusterlist(1))%z + delz/2 -1,real(gridsize)))+1
    
    53 if(compass .eqv. .true.) continue
    !write(6,*) 'com finish'
  end subroutine clustercom


  logical Function overlaps (pr,ch1,ch2,tempcoords)
    integer,intent(in) :: pr,ch1,ch2
    Type(protein),dimension(:),intent(in) :: tempcoords
    overlaps = .true.
    if(protcoords(pr,ch2)%x == tempcoords(ch1)%x.and. &
         protcoords(pr,ch2)%y == tempcoords(ch1)%y .and. &
         protcoords(pr,ch2)%z == tempcoords(ch1)%z) then
       overlaps = .false.
    end if
  end Function Overlaps

  logical Function overlapsclus (m,pr,ch1,ch2,tempcluscoord)
    integer,intent(in) :: pr,ch1,ch2,m
    Type(protein),dimension(:,:),intent(in) :: tempcluscoord
    overlapsclus = .true.
    if(protcoords(pr,ch2)%x == tempcluscoord(m,ch1)%x.and. &
         protcoords(pr,ch2)%y == tempcluscoord(m,ch1)%y .and. &
         protcoords(pr,ch2)%z == tempcluscoord(m,ch1)%z) then
       overlapsclus = .false.
    end if
  end Function Overlapsclus

  subroutine energy
    integer :: dx,dy,dz,m,l,f,g
    double precision :: initialenergy
    initialenergy = 0.0d0
    totalenergy = 0.0d0

    do m= 1,nprotein,1
       do l = 1,maxlength-3,1
          do f = l+3,maxlength,1

             dx = 4
             dy = 4
             dz = 4

             if(f < l-2 .or. f > l + 2) then
                if(abs(protcoords(m,l)%x - protcoords(m,f)%x) == 1 .or. abs(protcoords(m,l)%x - &
                     protcoords(m,f)%x) == (gridsize - 1)) then
                   dx = 1
                else if(abs(protcoords(m,l)%x - protcoords(m,f)%x) == 0) then
                   dx = 0
                end if
                if(dx == 4) goto 69 !check this

                if(abs(protcoords(m,l)%y - protcoords(m,f)%y) == 1 .or.abs(protcoords(m,l)%y - &
                     protcoords(m,f)%y) == (gridsize - 1)) then
                   dy = 1

                else if(abs(protcoords(m,l)%y - protcoords(m,f)%y) == 0) then
                   dy = 0
                end if
                if(dy ==4) goto 69
                if(abs(protcoords(m,l)%z - protcoords(m,f)%z) == 1 .or.abs(protcoords(m,l)%z &
                     - protcoords(m,f)%z) == (gridsize -1)) then
                   dz = 1
                else if(abs(protcoords(m,l)%z - protcoords(m,f)%z) == 0) then
                   dz = 0
                end if
                if(dz == 4) goto 69

                if(dx + dy + dz == 1) then
                   isobond(m,l,f) = 1
                   initialenergy = initialenergy + intraenergy
                else 
                   isobond(m,l,f) = 0
                end if
             end if
69           if((dx ==4) .or. (dy ==4) .or. (dz==4)) then
                isobond(m,l,f) = 0
             end if
          end do
       end do
    end do

!!!!!!! interatomic
    do m= 1,nprotein,1
       do g = m,nprotein,1
          do l = 1,maxlength,1
             do f = 1,maxlength,1
                if(g /=m) then
                   dx = 4
                   dy = 4
                   dz = 4
                   if(abs(protcoords(m,l)%x - protcoords(g,f)%x) == 1 .or. abs(protcoords(m,l)%x - &
                        protcoords(g,f)%x) == (gridsize - 1)) then
                      dx = 1
                   else if(abs(protcoords(m,l)%x - protcoords(g,f)%x) == 0) then
                      dx = 0
                   end if
                   if(dx == 4) goto 71
                   if(abs(protcoords(m,l)%y - protcoords(g,f)%y) == 1 .or.abs(protcoords(m,l)%y - &
                        protcoords(g,f)%y) == (gridsize - 1)) then
                      dy = 1
                   else if(abs(protcoords(m,l)%y - protcoords(g,f)%y) == 0) then
                      dy = 0
                   end if
                   if(dy == 4) goto 71
                   if(abs(protcoords(m,l)%z - protcoords(g,f)%z) == 1 .or.abs(protcoords(m,l)%z &
                        - protcoords(g,f)%z) == (gridsize -1)) then
                      dz = 1
                   else if(abs(protcoords(m,l)%z - protcoords(g,f)%z) == 0) then
                      dz = 0
                   end if
                   if(dz == 4) goto 71
                   if(dx + dy + dz == 1) then
                      intbond(m,g,l,f) = 1
                      initialenergy = initialenergy + interenergy
                   else 
                      intbond(m,g,l,f) = 0
                   end if
                end if
71              if(dx ==4 .or. dy ==4 .or. dz==4) then
                   intbond(m,g,l,f) = 0
                end if

             end do
          end do
       end do
    end do
    write(93,*) 'manual chaeck' ,initialenergy
    totalenergy = initialenergy

  end subroutine energy

  subroutine radiusofgyration
    integer::m,l
    double precision :: rog,totrog

    do m = 1,nprotein,1
       do l = 1, maxlength,1
          rog = (min(modulo(protcoords(m,l)%x-com(time,m)%x,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%x-com(time,m)%x,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(time,m)%y,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%y-com(time,m)%y,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(time,m)%z,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%z-com(time,m)%z,gridsize*1.0))))**2
          totrog = totrog + rog
       end do
    end do

    runningaveROG = runningaveROG + totrog
    write(97,*) time,totrog/(nprotein*maxlength),SQRT(totrog/(nprotein*maxlength))

  end subroutine radiusofgyration

  logical Function adjacent(pr,ch1,ch2,tempory)
    integer,intent(in):: pr,ch1,ch2
    integer::dx,dy,dz,delx,dely,delz
    Type(protein),dimension(:),intent(in) :: tempory
    dx = 4
    dy = 4
    dz = 4

    delx = abs(tempory(ch1)%x - protcoords(pr,ch2)%x)
    if(delx == 1 .or. delx == (gridsize - 1)) then
       dx = 1
    else if(delx == 0) then
       dx = 0
    end if
    if(dx == 4) then
       adjacent = .false.
       return
    end if

    dely = abs(tempory(ch1)%y - protcoords(pr,ch2)%y)
    if(dely == 1 .or.dely == (gridsize - 1)) then
       dy = 1
    else if(dely == 0) then
       dy = 0
    end if
    if(dy == 4)then
       adjacent = .false.
       return
    end if

    delz = abs(tempory(ch1)%z - protcoords(pr,ch2)%z)
    if(delz == 1 .or. delz == (gridsize -1)) then
       dz = 1
    else if(delz == 0) then
       dz = 0
    end if
    if(dz == 4) then
       adjacent = .false.
       return
    end if

    if(dx + dy + dz ==1) then
       adjacent = .true.
    else if(dx + dy + dz /= 1) then
       adjacent = .false.
    end if

  end Function adjacent

  subroutine rms
    double precision :: msdsum,msd
    integer::m
    double precision,dimension(:),allocatable :: msdx,msderrorx,msdy,msderrory,msdz,msderrorz
    allocate(msdx(nprotein))
    allocate(msderrorx(nprotein))
    allocate(msdy(nprotein))
    allocate(msderrory(nprotein))
    allocate(msdz(nprotein))
    allocate(msderrorz(nprotein))

    t = time + 1
    msd = 0.0d0
    msdsum = 0.0d0
    do m = 1, nprotein,1

       msdx(m) = (min(modulo(com(t,m)%x - com(t-1,m)%x,gridsize*1.0), &
            (gridsize - modulo(com(t,m)%x -com(t-1,m)%x,gridsize*1.0))))**2
       msdy(m) = (min(modulo(com(t,m)%y -com(t-1,m)%y,gridsize*1.0), &
            (gridsize - modulo(com(t,m)%y -com(t-1,m)%y,gridsize*1.0))))**2
       msdz(m) = (min(modulo(com(t,m)%z - com(t-1,m)%z,gridsize*1.0), &
            (gridsize - modulo(com(t,m)%z -com(t-1,m)%z,gridsize*1.0))))**2
     
       totrmsbrute = totrmsbrute + msdx(m) + msdy(m) + msdz(m)
    end do
    msd = sum(msdx) + sum(msdy) + sum(msdz)
    msdsum = msdsum + msd

    if(modulo(time,100) == 0) then
       write(29,*) nprotein*maxlength*(time-equilib),totrmsbrute
       !write(77,*) log(real(nprotein*maxlength*time)), log(totrmsbrute)
    end if
  end subroutine rms



  subroutine dataout
    integer::m,l
    !write(6,*) 'data'
    write(67,*) nprotein*maxlength
    write(67,*) ' '
    t = time 
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          write(67,*) m, 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
               2*protcoords(m,l)%z
       end do
    end do
  end subroutine dataout


  subroutine clusterenergy(cluster,tempcluscoord,energyofcluster)
    type(grouping),intent(inout) ::cluster
    type(protein),dimension(:,:),intent(in) :: tempcluscoord
    type(protein),dimension(:),allocatable :: tempcoord
    integer:: a,b,l,g,f
    doubleprecision,intent(inout) :: energyofcluster
    logical :: clusmove
    allocate(tempcoord(maxlength))

    energyofcluster = 0.0

    do a = 1,cluster%clustersize
       b = cluster%clusterlist(a)
       do g = 1,nprotein
          clusmove = .true.
          if(ANY(cluster%clusterlist(:) .eq. g)) goto 62
          do l = 1,maxlength,1
             tempcoord(l) = tempcluscoord(b,l)
             do f = 1,maxlength,1
                if(adjacent(g,l,f,tempcoord) .eqv. .true.) energyofcluster = energyofcluster + interenergy
             end do
          end do
62        if(clusmove .eqv. .true.) continue
       end do
    end do

  end subroutine clusterenergy


  subroutine length
    integer :: ddx,ddy,ddz,dxsum,dysum,dzsum,m,l
    t = time
    totchainlength = 0.0

    do m = 1,nprotein,1
       !chainlength = 0.0
       ddx = 0
       ddy = 0
       ddz = 0
       dxsum = 0
       dysum = 0
       dzsum = 0
       do l = 2,maxlength,1

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
    avechainlength = totchainlength/nprotein

    actualchainlength = actualchainlength + avechainlength
    write(91,*) t,avechainlength,actualchainlength/(time/10)
  end subroutine length

  subroutine clustertranslation(m,cluster)
    integer,intent(in)::m
    type(grouping),intent(inout) ::cluster
    integer :: steplength,dx,dy,dz
    type(protein),dimension(:,:),allocatable :: tempcluscoord
    double precision :: clusen,direction
    allocate(tempcluscoord(nprotein,maxlength))

    !m = int(ran2(seed)*(clusterno))+1
    !write(6,*) 'green'
      !write(6,*) 'maxclus', maxclussize
    !write(6,*) cluster%clustersize
    steplength = int(maxclussize/cluster%clustersize)
    direction = ran2(seed)

    if(direction <= 1.0/6) then
       dx = steplength
       dy = 0
       dz = 0
    else if(direction > 1.0/6 .and. direction <= 1.0/3) then
       dx = -steplength
       dy = 0
       dz = 0
    else if(direction > 1.0/3 .and. direction <= 1.0/2) then
       dx = 0
       dy = steplength
       dz = 0
    else if(direction > 1.0/2 .and. direction <= 2.0/3) then
       dx = 0
       dy = -steplength
       dz = 0
    else  if(direction > 2.0/3 .and. direction <= 5.0/6) then
       dx = 0
       dy = 0
       dz = steplength
    else if(direction > 5.0/6 .and. direction <= 1.0/1) then
       dx = 0
       dy = 0
       dz = -steplength
    end if

    call clustermove(cluster,m,dx,dy,dz,tempcluscoord)
    clusen = 0.0
    call clusterenergy(cluster,tempcluscoord,clusen)
    if(clusen == 1)  call updatecluspos(cluster,tempcluscoord)
  end subroutine clustertranslation


  subroutine clusterrotation(m,cluster)
    integer,intent(in) :: m
    type(grouping),intent(inout) ::cluster
    integer :: a,b,l,pr,st
    double precision :: rotatechoose,clusen
    integer,dimension(:,:),allocatable :: deltax,deltay,deltaz
    Type(protein),dimension(:,:),allocatable :: tempcluscoord
    logical :: clusmove
    clusmove = .true.
    allocate(tempcluscoord(nprotein,maxlength))
    t = time + 1
    !m = int(ran2(seed)*(clusterno))+1
    !write(6,*) 'pre com'
    call clustercom(m,cluster)
    !write(6,*) 'compass'
    allocate(deltax(nprotein,maxlength))
    allocate(deltay(nprotein,maxlength))
    allocate(deltaz(nprotein,maxlength))

    do a = 1,cluster%clustersize
       b = cluster%clusterlist(a)
       do l = 1, maxlength
          deltax(b,l) =   comcluster%x - protcoords(b,l)%x 
          if(deltax(b,l) > gridsize/2) deltax(b,l) =  deltax(b,l)-gridsize
          if(deltax(b,l) < -gridsize/2) deltax(b,l) = gridsize + deltax(b,l)
          deltay(b,l) =   comcluster%y - protcoords(b,l)%y
          if(deltay(b,l) > gridsize/2) deltay(b,l) =  deltay(b,l)-gridsize
          if(deltay(b,l) < -gridsize/2) deltay(b,l) = gridsize + deltay(b,l)
          deltaz(b,l) =   comcluster%z - protcoords(b,l)%z 
          if(deltaz(b,l) > gridsize/2) deltaz(b,l) =  deltaz(b,l)-gridsize
          if(deltaz(b,l) < -gridsize/2) deltaz(b,l) = gridsize + deltaz(b,l)
       end do
    end do

    rotatechoose = ran2(seed)

    do a = 1,cluster%clustersize
       b = cluster%clusterlist(a)    
       if(rotatechoose <= 1.0/9) call rotatex(b,1,-1,deltay,deltaz,tempcluscoord)
       if(rotatechoose > 1.0/9 .and. rotatechoose <= 2.0/9) call rotatex(b,-1,1,deltay,deltaz,tempcluscoord)
       if(rotatechoose > 2.0/9 .and. rotatechoose <= 1.0/3) call clusterflip(b,deltax,deltay,deltaz,tempcluscoord)    
       if(rotatechoose > 1.0/3 .and. rotatechoose <= 4.0/9) call rotatey(b,1,-1,deltax,deltaz,tempcluscoord)
       if(rotatechoose > 4.0/9 .and. rotatechoose <= 5.0/9) call rotatey(b,-1,1,deltax,deltaz,tempcluscoord)
       if(rotatechoose > 5.0/9 .and. rotatechoose <= 6.0/3) call clusterflip(b,deltax,deltay,deltaz,tempcluscoord)
       if(rotatechoose>6.0/9 .and. rotatechoose <= 7.0/9) call rotatez(b,1,-1,deltax,deltay,tempcluscoord)
       if(rotatechoose > 7.0/9 .and. rotatechoose <= 8.0/9) call rotatez(b,-1,1,deltax,deltay,tempcluscoord)
       if(rotatechoose > 8.0/9 .and. rotatechoose <= 1.0) call clusterflip(b,deltax,deltay,deltaz,tempcluscoord)


       do pr = 1,b
          if(ANY(cluster%clusterlist(:) .eq. pr)) goto 39
          do l = 1,maxlength
             do st = 1,maxlength,1
                if(overlapsclus(b,pr,l,st,tempcluscoord) .eqv. .false.) then
                   clusmove = .false.
                   goto 31
                end if
             end do
          end do
39        if(clusmove .eqv. .true.) continue
       end do

       do pr = b,nprotein
          if(ANY(cluster%clusterlist(:) .eq. pr)) goto 49
          do l = 1,maxlength
             do st = 1,maxlength,1
                if(overlapsclus(b,pr,l,st,tempcluscoord) .eqv. .false.) then
                   clusmove = .false.
                   goto 31
                end if
             end do
          end do
49        if(clusmove .eqv. .true.) continue
       end do

    end do

    clusen = 0.0
    call clusterenergy(cluster,tempcluscoord,clusen)
    if(clusen == 1)  call updatecluspos(cluster,tempcluscoord)


31  if(clusmove .eqv. .false.) continue
    !calculate the energy penalty and insure no overlaps

  end subroutine clusterrotation

  subroutine clusterflip(protno,deltax,deltay,deltaz,tempcluscoord)
    integer,intent(in)::protno
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%x = modulo(protcoords(protno,l)%x -2*deltax(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(protcoords(protno,l)%y-2*deltay(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(protcoords(protno,l)%z-2*deltaz(protno,l) -1,gridsize)+1
    end do

  end subroutine clusterflip

  subroutine rotatex(protno,cy,cz,deltay,deltaz,tempcluscoord)
    integer,intent(in)::protno,cy,cz
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltay,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%x = protcoords(protno,l)%x
       tempcluscoord(protno,l)%y = modulo(protcoords(protno,l)%y+cy*deltaz(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(protcoords(protno,l)%z+cz*deltay(protno,l) -1,gridsize)+1
    end do

  end subroutine rotatex

  subroutine rotatey(protno,cx,cz,deltax,deltaz,tempcluscoord)
    integer,intent(in)::protno,cx,cz
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%y = protcoords(protno,l)%y
       tempcluscoord(protno,l)%x = modulo(protcoords(protno,l)%x+cx*deltaz(protno,l)-1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(protcoords(protno,l)%z+cz*deltax(protno,l)-1,gridsize)+1
    end do

  end subroutine rotatey

  subroutine rotatez(protno,cx,cy,deltax,deltay,tempcluscoord)
    integer,intent(in)::protno,cx,cy
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay
    integer :: l

       do l = 1,maxlength
       tempcluscoord(protno,l)%z = protcoords(protno,l)%z
       tempcluscoord(protno,l)%x = modulo(protcoords(protno,l)%x+cx*deltay(protno,l)-1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(protcoords(protno,l)%y+cy*deltax(protno,l)-1,gridsize)+1
    end do

  end subroutine rotatez


  subroutine clustermove(cluster,chainno,delx,dely,delz,tempcluscoord)
    integer,intent(in) :: chainno,delx,dely,delz
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    type(grouping),intent(inout) ::cluster
    integer :: b,pr,st
    logical :: clusmove

    clusmove = .true.

    do b =1,maxlength
       tempcluscoord(chainno,b)%x = modulo(protcoords(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(protcoords(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(protcoords(chainno,b)%z + delz-1,gridsize)+1
    end do

    do pr = 1,chainno
       if(ANY(cluster%clusterlist(:) .eq. pr)) goto 39
       do st = 1,maxlength,1
          if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
             clusmove = .false.
             goto 31
          end if
       end do
39     if(clusmove .eqv. .true.) continue
    end do

    do pr = chainno,nprotein
       if(ANY(cluster%clusterlist(:) .eq. pr)) goto 49
       do st = 1,maxlength,1
          if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
             clusmove = .false.
             goto 31
          end if
       end do
49     if(clusmove .eqv. .true.) continue
    end do

   31  if(clusmove .eqv. .false.) continue
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove

  subroutine debug
    integer::m,l,f,g,dx,dy,dz

    do m = 1,nprotein,1
       do l = 1,maxlength,1
          do f = 1,nprotein,1
             do g = 1,maxlength,1
                if(m ==f .and. l == g) then
                   continue
                else if (m/= f .or. l/= g) then
                   !write(6,*) protcoords(m,l)%x, protcoords(f,g)%x
                   if(protcoords(m,l)%x == protcoords(f,g)%x  .and. protcoords(m,l)%y == protcoords(f,g)%y  &
                        .and. protcoords(m,l)%z == protcoords(f,g)%z ) then
                      !write(6,*) 'fail is due to', l,g
                      finalfail = .true.
                      fail = .true.
                      if (protcoords(m,l)%x == 0.0) write(6,*) 'x fail',l
                      if (protcoords(m,l)%y == 0.0) write(6,*) 'y fail',l
                      if (protcoords(m,l)%z == 0.0) write(6,*) 'z fail',l
                      write(6,*) 'FAILLLLLLLLLLLLLLLLLL'
                   end if
                end if
             end do
          end do

          if(l>1) then

             !write(6,*) 'l',l,'f',f

             if(abs(protcoords(m,l)%x - protcoords(m,l-1)%x) == 1 .or. abs(protcoords(m,l)%x - &
                  protcoords(m,l-1)%x) == (gridsize - 1)) then
                dx = 1
             else if(abs(protcoords(m,l)%x - protcoords(m,l-1)%x) == 0) then
                dx = 0
             end if
             !energypass = .false.
             if(abs(protcoords(m,l)%y - protcoords(m,l-1)%y) == 1 .or.abs(protcoords(m,l)%y - &
                  protcoords(m,l-1)%y) == (gridsize - 1)) then
                dy = 1
             else if(abs(protcoords(m,l)%y - protcoords(m,l-1)%y) == 0) then
                dy = 0
             end if
             !energypass = .false.
             if(abs(protcoords(m,l)%z - protcoords(m,l-1)%z) == 1 .or.abs(protcoords(m,l)%z &
                  - protcoords(m,l-1)%z) == (gridsize -1)) then
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
  end subroutine debug

  subroutine clustercount
    integer:: a,c,m,l,g,count,controlcount,counter,totalclusters,clusteredchains,chainnum,overallcount
    type(protein),dimension(:),allocatable :: tempcoord
    type(grouping),dimension(:),allocatable :: clus
    integer,dimension(:),allocatable :: clusterlist
    logical :: clusterpass,newcluster,newclusinit
    allocate(clus(nprotein))
    allocate(tempcoord(maxlength))
    allocate(clusterlist(nprotein))
    !allocate(clus%clusterlist(nprotein))
    totalclusters = 0
    clusteredchains = 0
    count = 1
    overallcount = 0
   do a = 1,nprotein
      clusterlist(a) = 0
      do c = 1,nprotein
         clus(a)%clusterlist(c) = 0
         end do
   end do
    do m = 1,nprotein,1
       if(ANY(clusterlist(:) .eq. m))then
          newclusinit = .false.
          goto 72
       end if
       clus(count)%clusterlist(1) = m
       clus(count)%clustersize = 1
       controlcount = 1
       clusterpass = .true.
       !write(6,*) 'a', count,m
       counter = 1
       do while (clusterpass .eqv. .true.)
          !counter = counter 
          chainnum = clus(count)%clusterlist(counter)
          !write(6,*) chainnum,'chainnum'
          clusterpass = .false.
          do l = 1,maxlength
             tempcoord(l)%x = protcoords(chainnum,l)%x
             tempcoord(l)%y = protcoords(chainnum,l)%y
             tempcoord(l)%z = protcoords(chainnum,l)%z
             do n = 1,nprotein
                
                newcluster = .false.
                if (ANY(clus(count)%clusterlist(:) .eq. n )) then
                   goto 25
                else
                   do g= 1,maxlength
                      if(adjacent(n,l,g,tempcoord) .eqv. .true.) then
                         !write(6,*) 'b1'
                         newcluster = .true.
                         clusterpass = .true.
                         goto 56
                      end if
                   end do
                end if
25              if(newcluster .eqv. .false. .and. clus(count)%clusterlist(counter) /= n) then
                   if(counter < clus(count)%clustersize) clusterpass = .true.
                   !write(6,*) 'b2',n
                end if
             end do
56           if(newcluster .eqv. .true.) then
                !write(6,*) 'c'
                counter = counter + 1
                overallcount = overallcount + 1
                !write(6,*) 'e',overallcount
                clus(count)%clustersize = counter
                !write(6,*) 'f'
                clus(count)%clusterlist(counter) = n
                !write(6,*) 'g'
                clusterlist(overallcount) = n
             !write(6,*) 'd'
             end if

          end do
       end do
       
       if(clus(count)%clustersize > 1) then  
          totalclusters = totalclusters + 1
          clusteredchains = clusteredchains + clus(count)%clustersize
          count = count + 1
       end if

72     if(newclusinit .eqv. .false.) continue
    end do

    if(totalclusters > 0)  write(82,*) time,clusteredchains/(1.0*totalclusters),clusteredchains,totalclusters


!if this crashes change the size of the array of type(grouping)
    
          !write(6,*) 'finish clustercount'
  end subroutine clustercount

  
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
