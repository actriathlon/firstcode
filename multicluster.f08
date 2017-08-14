program move

  implicit none
  
  double precision :: intraenergy,interenergy,totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength
  double precision :: actualchainlength,random
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug
  integer :: piv,right,endm,crank,rept,datayes,maxclussize,cranklimit
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib
  integer :: successful,reject,pivotlimit,maxpiv,maxlength1,maxlength2
  integer:: cracc,crrej,raacc,rarej,rerej,reacc,pvacc,pvrej,nprotein1,nprotein2
  integer,dimension(:,:),allocatable :: bonddd
  double precision :: crrat,emrat,rarat,pvrat,rerat
  integer,dimension(:,:,:), allocatable :: isobond
  integer,dimension(:,:,:,:),allocatable :: intbond
  real, external :: ran2
  real,dimension(:),allocatable :: variance,disp
  logical :: exist,fail,finalfail,debugyes,film
  real::start,finish

  type protein
     integer :: x,y,z,species
  end type protein

  type centremass
     double precision :: x,y,z
  end type centremass


  type(protein),dimension(:,:),allocatable :: protcoords
  type(centremass),dimension(:,:),allocatable :: com
  type(protein) :: comcluster


  !open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')
  open(29, file = 'rms.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  open(79, file = 'runningave.dat', action = 'write')
  open(93, file = 'energy.dat', action = 'write')
  open(82,file = 'clusterdata.dat', action = 'write')
  open(13,file='timedata.dat',action = 'write')
open(19,file='acceptance.dat',action  = 'write')
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
  call read_setup
  nprotein = nprotein1 + nprotein2
    !write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT

   open(67, file = 'move.xyz', action = 'write')
  allocate(protcoords(nprotein,maxlength))
  allocate(com(maxtime+1,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(isobond(nprotein,maxlength,maxlength))
  allocate(intbond(nprotein,nprotein,maxlength,maxlength))
allocate(bonddd(nprotein,maxlength))
  !allocate(cluster(nprotein))
  !allocate(comcluster(nprotein))
  N = maxlength
  interenergy =  -1.0 !-1.0 !-1.0d0
  intraenergy = -1.0  !-1.0 !-1.0d0

  count = 0
  time = 0
  totdisp = 0.0d0


  cracc = 0
  crrej = 0
  rerej = 0
  reacc = 0
  raacc = 0
  rarej = 0
  pvacc = 0
  pvrej = 0

  call foundation
!sets the limit of the pivot move
  pivotlimit = min((maxlength/2),maxpiv)
  !cranklimit = 5
write(6,*) 'a'
  !call dataout
  call comfind
write(6,*) 'b'
  !call clustercom
  actualchainlength = 0.0
  do time = 1,maxtime,1
     call CPU_TIME(start)
     fail = .false.
     if (mod(time,100) == 0) then
        if (film .eqv. .true.) call dataout
        !call bondcheck
     end if
!write(6,*) 'c'
!call positioning
!write(6,*) 'd'
call pickmoves
!write(6,*) 'e'
     call comfind
     !call vectorspin
     !call clustercom

     !call clusterposition

     if(modulo(time,1000) == 0) call clustercount
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
     if (modulo(time,100) == 0 .and. debugyes .eqv. .true.) then
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
     if(modulo(time,1000) == 0) then
        write(19,*) 'pivot',pvacc,pvrej, real(pvacc)/(pvacc+pvrej)
        write(19,*) 'reptation',reacc,rerej, real(reacc)/(reacc+rerej)
        write(19,*) 'right angle',raacc,rarej, real(raacc)/(raacc+rarej)
        write(19,*) 'crankshaft',cracc,crrej, real(cracc)/(cracc+crrej)
        end if
  end do
  !call error
  call energy
  write(6,*) 'average chain length squared is;' ,(actualchainlength/(maxtime/10)), &
       sqrt(actualchainlength/(maxtime/10))
  write(6,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))
  inquire(file = "lengthdataaverages.dat", exist = exist)
  if (exist) then
     open(30, file = "lengthdataaverages.dat", status = "old", position = "append",action ="write")
  else
     open(30, file = "lengthdataaverages.dat", status = "new",action ="write")
  end if
  write(30,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))!,2*(SQRT(totvar)/maxtime)/(actualchainlength/maxtime)
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


  subroutine vectorspin
    integer::m,p,f,dum1,dum2,l,no1,no2
    integer,dimension(:),allocatable:: tbo
    type(protein),dimension(:),allocatable::tempcoord
    double precision :: deltaenergy
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:), allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tbo(maxlength))
    allocate(tempcoord(maxlength))
    m =int(ran2(seed)*(nprotein))+1
    !dum1 =int(ran2(seed)*(maxlength))+1
    !dum2 =int(ran2(seed)*(maxlength))+1
    !p = min(dum1,dum2)
    !f = max(dum1,dum2)
      if(protcoords(m,1)%species == 1) maxlength = maxlength1
      if(protcoords(m,1)%species == 1) maxlength = maxlength2

    !do a = 1,nprotein*maxlength
       
       !do l = p,f
          !write(6,*) totalenergy
    l = 1
    if(protcoords(m,l+1)%x - protcoords(m,l)%x == 1) no1 = 1 
             if(protcoords(m,l+1)%x - protcoords(m,l)%x == -1) no1 = -1
             if(protcoords(m,l+1)%y - protcoords(m,l)%y == 1) no1 = 2
             if(protcoords(m,l+1)%y - protcoords(m,l)%y == -1) no1 = -2
             if(protcoords(m,l+1)%z - protcoords(m,l)%z == 1) no1 = 3
             if(protcoords(m,l+1)%z - protcoords(m,l)%z == -1) no1 = -3
             no2 = 7
             call bonddirection(5,no1,no2,tbo(l))
          tempcoord(l)%x = protcoords(m,l)%x
          tempcoord(l)%y = protcoords(m,l)%y
          tempcoord(l)%z = protcoords(m,l)%z
             do l = 2,maxlength-1
             if(protcoords(m,l+1)%x - protcoords(m,l)%x == 1) no1 = 1 
             if(protcoords(m,l+1)%x - protcoords(m,l)%x == -1) no1 = -1
             if(protcoords(m,l+1)%y - protcoords(m,l)%y == 1) no1 = 2
             if(protcoords(m,l+1)%y - protcoords(m,l)%y == -1) no1 = -2
             if(protcoords(m,l+1)%z - protcoords(m,l)%z == 1) no1 = 3
             if(protcoords(m,l+1)%z - protcoords(m,l)%z == -1) no1 = -3
             if(protcoords(m,l-1)%x - protcoords(m,l)%x == 1) no2 = 1 
             if(protcoords(m,l-1)%x - protcoords(m,l)%x == -1) no2 = -1
             if(protcoords(m,l-1)%y - protcoords(m,l)%y == 1) no2 = 2
             if(protcoords(m,l-1)%y - protcoords(m,l)%y == -1) no2 = -2
             if(protcoords(m,l-1)%z - protcoords(m,l)%z == 1) no2 = 3
             if(protcoords(m,l-1)%z - protcoords(m,l)%z == -1) no2 = -3
             call bonddirection(4,no1,no2,tbo(l))
             !write(6,*) tbo(l)
             tempcoord(l)%x = protcoords(m,l)%x
          tempcoord(l)%y = protcoords(m,l)%y
          tempcoord(l)%z = protcoords(m,l)%z
             end do
                l = maxlength
             if(protcoords(m,l-1)%x - protcoords(m,l)%x == 1) no2 = 1 
             if(protcoords(m,l-1)%x - protcoords(m,l)%x == -1) no2 = -1
             if(protcoords(m,l-1)%y - protcoords(m,l)%y == 1) no2 = 2
             if(protcoords(m,l-1)%y - protcoords(m,l)%y == -1) no2 = -2
             if(protcoords(m,l-1)%z - protcoords(m,l)%z == 1) no2 = 3
             if(protcoords(m,l-1)%z - protcoords(m,l)%z == -1) no2 = -3
             no1 = 7
             call bonddirection(5,no1,no2,tbo(l))
  
          tempcoord(l)%x = protcoords(m,l)%x
          tempcoord(l)%y = protcoords(m,l)%y
          tempcoord(l)%z = protcoords(m,l)%z

          deltaenergy = 0.0d0
         !write(6,*) deltaenergy
          call intrabondvector(m,tempcoord,tempisobond,deltaenergy,tbo)
          !write(6,*) deltaenergy
    call interbondallocate(m,1,maxlength,tempcoord,tempintbond,deltaenergy,tbo) !error here
     !write(6,*) deltaenergy
    if(Energydecision(deltaenergy) .eqv. .True.) then
              totalenergy = totalenergy + deltaenergy
       !write(6,*) totalenergy
       call updatepos(m,1,maxlength,tempcoord,tbo)
       call updateintrabondvector(m,tempisobond)
       call updateinterbond(m,1,maxlength,tempintbond)
    else 
       continue
    end if

  end subroutine vectorspin

    subroutine updateintrabondvector(chainnum,tempisobond)
    integer,intent(in):: chainnum
    integer ::g,beadnum
    integer,dimension(:,:),intent(in) :: tempisobond
    do beadnum = 1,maxlength-1
       do g = beadnum+3,maxlength
          if(g<beadnum-2) isobond(chainnum,g,beadnum) = tempisobond(g,beadnum) 
          if(g> beadnum+2) isobond(chainnum,beadnum,g) = tempisobond(beadnum,g)
       end do
    end do
  end subroutine updateintrabondvector

  subroutine intrabondvector(chain1,tempcoord,tempisobond,deltaenergy,tbo)

    integer,intent(in)::chain1
        integer,dimension(:),intent(inout)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer::bead1,bead2,bdir
    Type(protein),dimension(:),intent(in) :: tempcoord
    logical :: adjver
    
    do bead1 = 1,maxlength-3
       do bead2 = bead1+3,maxlength
             if(isobond(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraenergy
             call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
             if ((adjver.eqv. .true.) .and. &
                  (bdir == -1*tbo(bead1)) .and. (bdir == tbo(bead2))) then
                deltaenergy = deltaenergy + intraenergy
                tempisobond(bead1,bead2) = 1
             else 
                tempisobond(bead1,bead2) = 0
             end if
       end do
    end do

  end subroutine intrabondvector

  
  subroutine foundation
    !read in initial positions of the protein
    integer::m,l,no1,no2
    character(len = 10) :: BIN
    read(23,*) BIN
    write(6,*) maxlength

    !write(67,*) ' '
    do m = 1,nprotein,1
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
       if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do l = 1,maxlength,1
          !write(6,*) l
          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
       end do
    end do
    call energy
   ! call debug

    do m = 1,nprotein
       l = 1
       !bonddd(m,l) = 0
       if(protcoords(m,l+1)%x - protcoords(m,l)%x == 1) no1 = 1 
       if(protcoords(m,l+1)%x - protcoords(m,l)%x == -1) no1 = -1
       if(protcoords(m,l+1)%y - protcoords(m,l)%y == 1) no1 = 2
       if(protcoords(m,l+1)%y - protcoords(m,l)%y == -1) no1 = -2
       if(protcoords(m,l+1)%z - protcoords(m,l)%z == 1) no1 = 3
       if(protcoords(m,l+1)%z - protcoords(m,l)%z == -1) no1 = -3
       no2 = 7
       call bonddirection(5,no1,no2,bonddd(m,l))
       !write(6,*) bonddd(m,l)
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
       if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do l =2,maxlength-1
          !bonddd(m,l) = 0
       if(protcoords(m,l+1)%x - protcoords(m,l)%x == 1) no1 = 1 
       if(protcoords(m,l+1)%x - protcoords(m,l)%x == -1) no1 = -1
       if(protcoords(m,l+1)%y - protcoords(m,l)%y == 1) no1 = 2
       if(protcoords(m,l+1)%y - protcoords(m,l)%y == -1) no1 = -2
       if(protcoords(m,l+1)%z - protcoords(m,l)%z == 1) no1 = 3
       if(protcoords(m,l+1)%z - protcoords(m,l)%z == -1) no1 = -3
       if(protcoords(m,l-1)%x - protcoords(m,l)%x == 1) no2 = 1 
       if(protcoords(m,l-1)%x - protcoords(m,l)%x == -1) no2 = -1
       if(protcoords(m,l-1)%y - protcoords(m,l)%y == 1) no2 = 2
       if(protcoords(m,l-1)%y - protcoords(m,l)%y == -1) no2 = -2
       if(protcoords(m,l-1)%z - protcoords(m,l)%z == 1) no2 = 3
       if(protcoords(m,l-1)%z - protcoords(m,l)%z == -1) no2 = -3
      call bonddirection(4,no1,no2,bonddd(m,l))
    end do
    l = maxlength
           !bonddd(m,l) = 0
       if(protcoords(m,l-1)%x - protcoords(m,l)%x == 1) no2 = 1 
       if(protcoords(m,l-1)%x - protcoords(m,l)%x == -1) no2 = -1
       if(protcoords(m,l-1)%y - protcoords(m,l)%y == 1) no2 = 2
       if(protcoords(m,l-1)%y - protcoords(m,l)%y == -1) no2 = -2
       if(protcoords(m,l-1)%z - protcoords(m,l)%z == 1) no2 = 3
       if(protcoords(m,l-1)%z - protcoords(m,l)%z == -1) no2 = -3
       no1 = 7
       call bonddirection(5,no1,no2,bonddd(m,l))

         !if(bonddd(m,l) == 0) write(6,*) 'fail'
       end do
                  
  end subroutine foundation

  subroutine clusterposition
    integer::scans,m,clsize
    real :: decider
    integer,dimension(:),allocatable:: clno
    allocate(clno(nprotein))
    m =int(ran2(seed)*(nprotein))+1
    !call clustercount(clno)
    call clusterassign(m,clno,clsize)
    if(1/clsize > ran2(seed))then
       decider = ran2(seed) - 0.5
       !do scans = 1,nprotein
          
          if(decider > 0.0) then
             call clustertranslation(m,clno,clsize)
          else if(decider <= 0.0) then
             !write(6,*) 'turn'
             call clusterrotation(m,clno,clsize)
          end if
       !end do
    end if
  end subroutine clusterposition

subroutine pickmoves
  integer:: scans
double precision :: decide
  
  do scans = 1,nprotein*maxlength
     decide = ran2(seed)
     !write(6,*) 'f'
     if(decide <= 8.0/11) call vectorspin
     if(decide > 8.0/11 .and. decide <= 10.0/11) call positioning
     if(decide > 10.0/11 .and. decide <= 1.0) call clusterposition
     end do
end subroutine pickmoves

subroutine rotate(m,l,b,tempcoords,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz)
  integer,intent(in):: m,l,b,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz
  type(protein),dimension(:),intent(inout) :: tempcoords
  !write(6,*) dx,dy,dz
  tempcoords(b)%x = modulo(protcoords(m,l)%x + (xx*dx) + (xy*dy) + (xz*dz) -1,gridsize)+1
  tempcoords(b)%y = modulo(protcoords(m,l)%y + (yx*dx) + (yy*dy) + (yz*dz)-1,gridsize)+1
    tempcoords(b)%z = modulo(protcoords(m,l)%z + (zx*dx) + (zy*dy) + (zz*dz)-1,gridsize)+1

  
end subroutine rotate

  subroutine positioning
    !Selects move to be carried out on the proteins
    integer :: nmoves,scan,m,l
    logical :: run,run2
    double precision:: choose
    integer ::randomz

    t = time + 1
    !do scan = 1,nprotein*maxlength,1
       count = count + 1
       run = .true.
       run2 = .true. 


       nmoves = piv + crank + right + rept
       randomz = int(ran2(seed)*nmoves)+1
       
       if (randomz == crank .and. crank == 1) then
          !write(6,*) 'crank' 
          call crankshaftmove
       else if (randomz == (crank + right) .and. right ==1) then
          !write(6,*) 'right' 
          call rightanglemove
       else if (randomz == (crank + right + rept) .and. rept ==1) then
          !write(6,*) 'reptation' 
          call reptation
       else if (randomz == (crank + right + rept + piv) .and. piv == 1) then 
          m =int(ran2(seed)*(nprotein))+1
          if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
          if(ran2(seed) -0.5 > 0) then
             l = int(ran2(seed)*(pivotlimit))+1
          else
             l = maxlength -  int(ran2(seed)*(pivotlimit))
          end if
              call pivot(m,l)
       end if
    !end do
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
    integer :: dummy,dummy2,p,s,m,l,str,st,pr,sign,probs
    integer,dimension(:),allocatable::dx1,dy1,dz1
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbo
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(dx1(maxlength))
    allocate(dy1(maxlength))
    allocate(dz1(maxlength))
    allocate(tbo(maxlength))

    !this means that m/=0 and is up to nprotein
    m =int(ran2(seed)*(nprotein))+1
    if(protcoords(m,1)%species == 1) maxlength = maxlength1
    if(protcoords(m,1)%species == 1) maxlength = maxlength2
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
    if ((p-l < 3) .or.( p-l > cranklimit)) then
       goto 71
    end if

    probs = int(ran2(seed)*3) +1
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


    if ((pivotx .eqv. .false.) .and. (pivoty .eqv. .false.) .and. (pivotz .eqv. .false.)) then
       crankcont = .false.
       goto 43
    end if

    crankcont = .true.

    if (pivotx .eqv. .true.) then

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
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,-1)
             end if
          end do
          
       else if (probs == 2) then
          
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,1,0,0,0,0,1,0,-1,0,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,1)
             end if
          end do
          
       else if (probs == 3) then
          
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,1,0,0,0,-1,0,0,0,-1,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 1) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,1,0)
             end if
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
                    dy1(s) = protcoords(m,s)%y - protcoords(m,l)%y
          if (dy1(s) > gridsize/2) dy1(s) = dy1(s) - gridsize
          if (dy1(s) < -gridsize/2) dy1(s) = dy1(s) + gridsize

       end do

       if (probs == 1) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,0,-1,0,1,0,1,0,0,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,2,-1)
                end if
          end do

       else if (probs ==2) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,0,1,0,1,0,-1,0,0,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,2,1)
                end if
          end do

       else if (probs == 3) then

          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,-1,0,0,0,1,0,0,0,-1,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 2) then
                tbo(s) = bonddd(m,s)
             else
             call rbond(m,s,bonddd(m,s),tbo,2,0)
          end if
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
          dz1(s) = protcoords(m,s)%z - protcoords(m,l)%z
          if (dz1(s) > gridsize/2) dz1(s) = dz1(s) - gridsize
          if (dz1(s) < -gridsize/2) dz1(s) = dz1(s) + gridsize
       end do

       if (probs == 1) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,1,0,-1,0,0,0,0,1,dx1(s),dy1(s),dz1(s))
             if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,3,1)
             end if
          end do
          
       else if (probs ==2) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,0,-1,0,1,0,0,0,0,1,dx1(s),dy1(s),dz1(s))
             
             if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
                call rbond(m,s,bonddd(m,s),tbo,3,-1)
end if
             end do

       else if (probs ==3 ) then
          do s = l+1,p-1,1
             call rotate(m,l,s,tempcoord,-1,0,0,0,-1,0,0,0,1,dx1(s),dy1(s),dz1(s))
                          if(abs(bonddd(m,s)) == 3) then
                tbo(s) = bonddd(m,s)
             else
             call rbond(m,s,bonddd(m,s),tbo,3,0)
end if
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


    call intrabondallocateupper(m,l+1,p-1,1,l,tempcoord,tempisobond,deltaenergy,tbo)
    call intrabondallocatelower(m,l+1,p-1,p,maxlength,tempcoord,tempisobond,deltaenergy,tbo)

    call interbondallocate(m,l+1,p-1,tempcoord,tempintbond,deltaenergy,tbo)

    if(Energydecision(deltaenergy) .eqv. .True.) then
       totalenergy = totalenergy + deltaenergy
       call updatepos(m,l+1,p-1,tempcoord,tbo)
       call updateintrabond(m,l+1,p-1,1,l,tempisobond)
       call updateintrabond(m,l+1,p-1,p,maxlength,tempisobond)
       call updateinterbond(m,l+1,p-1,tempintbond)
    else 
       crankcont = .false.
       goto 43
    end if

    !       write(6,*) time, 'crank'
    successful = successful + 1
    call counts(1,0,0,0)
43  if (crankcont .eqv. .false.) then
       reject = reject + 1
       call counts(2,0,0,0)
    end if

  end subroutine crankshaftmove

subroutine rbond(m,s,bnd,tbo,axis,rotor)
  integer,intent(in):: m,s,bnd,axis,rotor
  integer,dimension(:),intent(inout)::tbo
  logical :: bm
  integer :: sign

    !write(6,*) 'fd'
bm = .true.
  if(bnd == axis) then
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
    logical :: rac
    integer:: st,pr
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable::tbo
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(tbo(maxlength))

    m =int(ran2(seed)*(nprotein))+1
    if(protcoords(m,1)%species == 1) maxlength = maxlength1
    if(protcoords(m,1)%species == 1) maxlength = maxlength2
    ! l cannot be equal to 1 or maxlength
    l = int(ran2(seed)*(maxlength-2))+2
    !write(6,*) 'm =', m, 'l = ', l

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

       tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
       tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
       tempcoord(l)%z = modulo(protcoords(m,l)%z+dz-1,gridsize)+1
       tbo(l) = -1*bonddd(m,l)
       
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

       call intrabondallocateupper(m,l,l,1,l-3,tempcoord,tempisobond,deltaenergy,tbo)
       call intrabondallocatelower(m,l,l,l+3,maxlength,tempcoord,tempisobond,deltaenergy,tbo)

       call interbondallocate(m,l,l,tempcoord,tempintbond,deltaenergy,tbo)


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l,l,tempcoord,tbo)
          call updateintrabond(m,l,l,1,l-3,tempisobond)
          call updateintrabond(m,l,l,l+3,maxlength,tempisobond)
          call updateinterbond(m,l,l,tempintbond)
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

  subroutine pivot(m,l)
    !rotates a section of chain around the selected bead
    integer,intent(in)::m,l
    integer :: b,g,xhold,yhold,zhold,str,st,pr,choose3
    logical :: pivcont
    double precision :: deltaenergy
    integer,dimension(:),allocatable :: delx,dely,delz
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:,:,:),allocatable :: tempintbond
    integer,dimension(:),allocatable::tbo
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(delx(maxlength))
    allocate(dely(maxlength))   
    allocate(delz(maxlength))
    allocate(tbo(maxlength))
    choose3 = int(ran2(seed)*5) + 1
    pivcont = .true. 

    t = time + 1

    if(l > maxlength/2) then
       !write(6,*) 'one'
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
       if (choose3 ==1) then

          do b = l+1,maxlength,1
             call rotate(m,l,b,tempcoord,0,-1,0,1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
                call rbond(m,b,bonddd(m,b),tbo,3,-1)
             end if
          end do

       else if (choose3 == 2) then
          do b =l+1,maxlength,1
             call rotate(m,l,b,tempcoord,0,1,0,-1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,3,1)
          end if
       end do

       else if (choose3 ==3) then

          do b =l+1,maxlength,1
             call rotate(m,l,b,tempcoord,-1,0,0,0,-1,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,3,0)
end if
          end do

       else if (choose3 ==4) then

          do b =l+1,maxlength,1
             call rotate(m,l,b,tempcoord,0,0,-1,0,1,0,1,0,0,delx(b),dely(b),delz(b))
                  if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,2,1)
          end if
       end do

       else if (choose3 == 5) then

          do b = l+1,maxlength,1
             call rotate(m,l,b,tempcoord,0,0,1,0,1,0,-1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,2,-1)
end if
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


       call intrabondallocateupper(m,l+1,maxlength,1,l,tempcoord,tempisobond,deltaenergy,tbo)

       call interbondallocate(m,l+1,maxlength,tempcoord,tempintbond,deltaenergy,tbo)

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,l+1,maxlength,tempcoord,tbo)
          call updateintrabond(m,l+1,maxlength,1,l,tempisobond)
          call updateinterbond(m,l+1,maxlength,tempintbond)
          successful = successful + 1
          call counts(0,0,1,0)

       else
          pivcont = .false.
          goto 75
       end if

    else if (l <= maxlength/2) then
       !write(6,*) 'two'
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
             !write(6,*) 'a'
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,3,1)
          end if
       end do

       else if (choose3 == 2) then 
          do b =1,l-1,1
             call rotate(m,l,b,tempcoord,0,1,0,-1,0,0,0,0,1,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 3) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,3,-1)
          end if
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

       else if (choose3==4) then
          do b =1,l-1,1
             call rotate(m,l,b,tempcoord,0,0,-1,0,1,0,1,0,0,delx(b),dely(b),delz(b))
             if(abs(bonddd(m,b)) == 2) then
                tbo(b) = bonddd(m,b)
             else
             call rbond(m,b,bonddd(m,b),tbo,2,1)
          end if
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
!write(6,*) '75'
       call intrabondallocatelower(m,1,l-1,l,maxlength,tempcoord,tempisobond,deltaenergy,tbo)


       call interbondallocate(m,1,l-1,tempcoord,tempintbond,deltaenergy,tbo)

       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          call updatepos(m,1,l-1,tempcoord,tbo)
          call updateintrabond(m,1,l-1,l,maxlength,tempisobond)
          call updateinterbond(m,1,l-1,tempintbond)
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


  subroutine intrabondallocatelower(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,deltaenergy,tbo)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    integer,dimension(:),intent(in)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer::bead1,bead2,bdir
    Type(protein),dimension(:),intent(in) :: tempcoord
    logical :: adjver
    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead2 > bead1 + 2) then
             if(isobond(chain1,bead1,bead2) == 1) deltaenergy = deltaenergy - intraenergy
             call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
             if ((adjver.eqv. .true.) .and. &
                  (bdir == -1*tbo(bead1)) .and. (bdir == bonddd(chain1,bead2))) then
                deltaenergy = deltaenergy + intraenergy
                tempisobond(bead1,bead2) = 1
             else 
                tempisobond(bead1,bead2) = 0
             end if
          end if
       end do
    end do
!write(6,*) 'intra lower energy', deltaenergy 
    
  end subroutine intrabondallocatelower

  subroutine intrabondallocateupper(chain1,bead1min,bead1max,bead2min,bead2max,tempcoord,tempisobond,deltaenergy,tbo)

    integer,intent(in)::chain1,bead1min,bead1max,bead2min,bead2max
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:),intent(inout) :: tempisobond
    integer,dimension(:),intent(in)::tbo
    integer::bead1,bead2,bdir
    Type(protein),dimension(:),intent(in) :: tempcoord
logical :: adjver
    do bead1 = bead1min,bead1max,1
       do bead2 = bead2min,bead2max,1
          if(bead1 > bead2 + 2) then
             if(isobond(chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - intraenergy
             call adjacent(chain1,bead1,bead2,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and. &
                  (bdir == -1*tbo(bead1)) .and. bdir== bonddd(chain1,bead2)) then
                !write(6,*) bdir,tbo(bead1),bonddd(chain1,bead2),adjver
                deltaenergy = deltaenergy + intraenergy
                tempisobond(bead2,bead1) = 1
             else 
                tempisobond(bead2,bead1) = 0
             end if
          end if
       end do
    end do
    !write(6,*) 'intra upper energy', deltaenergy 

  end subroutine intrabondallocateupper

  subroutine interbondallocate(chain1,beadmin,beadmax,tempcoord,tempintbond,deltaenergy,tbo)

    integer,intent(in)::chain1,beadmin,beadmax
    integer,dimension(:),intent(in)::tbo
    double precision,intent(inout)::deltaenergy
    integer,dimension(:,:,:),intent(inout) :: tempintbond
    integer::bead1,bead2,chain2,bdir
    Type(protein),dimension(:),intent(in) :: tempcoord
    logical :: adjver

    if(nprotein == 1) return
    
    do chain2 = 1,chain1-1,1
       do bead1 = beadmin,beadmax,1
          if(protcoords(chain2,1)%species == 1) maxlength = maxlength1
          if(protcoords(chain2,1)%species == 1) maxlength = maxlength2
          do bead2 = 1,maxlength,1
             if(intbond(chain2,chain1,bead2,bead1) == 1) deltaenergy = deltaenergy - interenergy
             !write(6,*) tbo(bead1)
             call adjacent(chain2,bead1,bead2,tempcoord,adjver,bdir)
             if((adjver .eqv. .true.) .and. (bdir == -1*tbo(bead1)) &
                  .and. (bdir ==bonddd(chain2,bead2))) then
                deltaenergy = deltaenergy + interenergy
                !write(6,*) bdir ,tbo(bead1),bonddd(chain2,bead2)
                tempintbond(chain2,bead1,bead2) = 1
             else
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do

    do chain2 = chain1+1,nprotein,1
       do bead1 = beadmin,beadmax,1
             if(protcoords(chain2,1)%species == 1) maxlength = maxlength1
          if(protcoords(chain2,1)%species == 1) maxlength = maxlength2
          do bead2 = 1,maxlength,1
             if(intbond(chain1,chain2,bead1,bead2) == 1) deltaenergy = deltaenergy - interenergy
             !write(6,*) tbo(bead1)
             call adjacent(chain2,bead1,bead2,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and. (bdir == -1*tbo(bead1)) .and.&
                  (bdir == bonddd(chain2,bead2))) then
                !write(6,*) bdir ,tbo(bead1),bonddd(chain2,bead2),chain1,chain2
                deltaenergy = deltaenergy + interenergy
                tempintbond(chain2,bead1,bead2) = 1
             else
                tempintbond(chain2,bead1,bead2) = 0
             end if
          end do
       end do
    end do

    !write(6,*) 'inter energy', deltaenergy 
  end subroutine interbondallocate

  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,tbo)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax
    integer,dimension(:),intent(inout)::tbo
    Type(protein),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding

    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
       bonddd(chainnum,beadnum) = tbo(beadnum)
    end do

  end subroutine updatepos

  subroutine updatecluspos(clno,tempcluscoord,chnum,ctbo)
    !updates bead positions
    integer,dimension(:),intent(in) :: clno
    integer,dimension(:,:),intent(in)::ctbo
    integer,intent(in) :: chnum
    Type(protein),dimension(:,:),intent(in) :: tempcluscoord
    integer ::b,l,a
    !moves beads to new positions and reassigns isobonding

    do a = 1,nprotein
       if(clno(a) == clno(chnum)) then
          if(protcoords(a,1)%species == 1) maxlength = maxlength1
          if(protcoords(a,1)%species == 1) maxlength = maxlength2
       do l = 1,maxlength,1
          protcoords(a,l)%x = tempcluscoord(a,l)%x
          protcoords(a,l)%y = tempcluscoord(a,l)%y
          protcoords(a,l)%z = tempcluscoord(a,l)%z
          bonddd(a,l) = ctbo(a,l)
       end do
       end if
    end do
  end subroutine updatecluspos

  subroutine bondcheck

    integer:: m,l,f
    !open(12, file = 'bondintra.dat', action = 'write')
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


  subroutine clusterassign(chainnum,clno,clsize)
    integer,intent(in) :: chainnum
    integer,dimension(:),intent(inout):: clno
    integer,intent(inout):: clsize
    integer:: m,l,g,f,clcount,bdir,maxlengthss
    logical :: clusyes,adjver
    type(protein),dimension(:),allocatable :: tempcoord
    !allocate(clno(nprotein))
    allocate(tempcoord(maxlength))

    do m = 1,nprotein
       clno(m) = m
    end do


    do m = 1,nprotein-1
          if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do g = m,nprotein
          if((g/=m) .and. (clno(m) /= clno(g))) then
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength1
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength2
             do l = 1,maxlength
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z
                do f = 1,maxlengthss
                   call adjacent(g,l,f,tempcoord,adjver,bdir)
                   if(adjver.eqv. .true. .and. bdir == bonddd(g,f) .and. &
                        bdir == -1*bonddd(m,l)) then
                      clno(g) = clno(m)
                      clusyes = .true.
                      goto 13
                   end if
                end do
             end do
             clusyes = .true.
13           if(clusyes .eqv. .true.) continue
          end if
       end do
    end do

    clsize = 1
    do m = 1,nprotein
       if(clno(m) == clno(chainnum)) clsize = clsize + 1
    end do
    
  end subroutine clusterassign


  subroutine updateinterbond(chainnum,beadmin,beadmax,tempintbond)
    integer,intent(in):: chainnum,beadmin,beadmax
    integer,dimension(:,:,:),intent(in) :: tempintbond
    integer ::g,beadnum,chain2
    do beadnum = beadmin,beadmax,1
       do chain2 = 1,chainnum-1,1
             if(protcoords(chain2,1)%species == 1) maxlength = maxlength1
          if(protcoords(chain2,1)%species == 1) maxlength = maxlength2
          do g = 1,maxlength,1
             intbond(chain2,chainnum,g,beadnum) = tempintbond(chain2,beadnum,g)
          end do
       end do
       do chain2 = chainnum+1,nprotein,1
             if(protcoords(chain2,1)%species == 1) maxlength = maxlength1
          if(protcoords(chain2,1)%species == 1) maxlength = maxlength2
          do g = 1,maxlength,1
             intbond(chainnum,chain2,beadnum,g) = tempintbond(chain2,beadnum,g)
          end do
       end do
    end do

  end subroutine updateinterbond

  subroutine reptation
    !perform reptation
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir
    logical :: reptcont,adjver
    type(protein),dimension(:),allocatable :: tempcoord
    integer,dimension(:,:), allocatable :: tempisobond
    integer,dimension(:),allocatable::tbo
    integer,dimension(:,:,:),allocatable :: tempintbond
    allocate(tempisobond(maxlength,maxlength))
    allocate(tempintbond(nprotein,maxlength,maxlength))
    allocate(tempcoord(maxlength))
    allocate(tbo(maxlength))

    t = time + 1
    choose =  ran2(seed) - 0.5
    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
       if(protcoords(m,1)%species == 1) maxlength = maxlength2
    
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

       direc = int(ran2(seed)*5) +1
       if((direc == g1) .and. (g1 == 1)) then
          dx = 1
          call bonddirection(5,1,8,tbo(1))
       else if((direc == (g1 +g2))  .and.  (g2 == 1)) then
          dx = -1
          call bonddirection(5,-1,8,tbo(1))
       else if((direc == (g1+g2+g3)) .and.  (g3 == 1)) then
          dy = 1
          call bonddirection(5,2,8,tbo(1))
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          call bonddirection(5,-2,8,tbo(1))
       else if((direc == (g1+g2+g3+g4+g5)) .and. (g5 == 1)) then
          dz = 1
          call bonddirection(5,3,8,tbo(1))
       else if((direc == (g1+g2+g3+g4+g5+g6)) .and.  (g6 == 1)) then
          dz = -1
          call bonddirection(5,-3,8,tbo(1))
       end if
       
       tempcoord(1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(1)%z = modulo(protcoords(m,1)%z + dz-1,gridsize)+1

       deltaenergy = 0.0d0
       do pr =1,nprotein    
          do st = 1,maxlength,1
             if((pr /=m) .or. (st > 1) .or. (st < maxlength))then
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
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.) .and. (bdir== -1*tbo(1)) .and. (bdir== bonddd(m,st))) then
             deltaenergy = deltaenergy + intraenergy
             tempisobond(1,st+1) = 1
          else
             tempisobond(1,st+1) = 0
          end if
       end do

       do pr = 1,m-1,1
          do st = 1,maxlength
             if(intbond(pr,m,st,maxlength) == 1) deltaenergy = deltaenergy - interenergy
             call adjacent(pr,1,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and.(bdir== -1*tbo(1)) .and. (bdir== bonddd(pr,st))) then
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
             call adjacent(pr,1,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and.(bdir== -1*tbo(1)) .and. (bdir== bonddd(pr,st))) then
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
             bonddd(m,l) = bonddd(m,l-1)

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

          call updatepos(m,1,1,tempcoord,tbo)
          call updateintrabond(m,1,1,4,maxlength,tempisobond)
          call updateinterbond(m,1,1,tempintbond)

          successful = successful + 1
          reptforward = reptforward + 1
          call counts(0,0,0,1)
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
       if((direc == g1) .and. (g1 == 1)) then
          dx = 1
          call bonddirection(5,1,8,tbo(maxlength))
       else if((direc ==(g1 +g2)) .and.  (g2 == 1)) then
          dx = -1
          call bonddirection(5,-1,8,tbo(maxlength))
       else if((direc == (g1+g2+g3)).and.  (g3 == 1)) then
          dy = 1
          call bonddirection(5,2,8,tbo(maxlength))
       else if((direc == (g1+g2+g3+g4)) .and.  (g4 == 1)) then
          dy = -1
          call bonddirection(5,2,8,tbo(maxlength))
       else if((direc == (g1+g2+g3+g4+g5)) .and.( g5 == 1)) then
          dz = 1
          call bonddirection(5,3,8,tbo(maxlength))
       else if((direc == (g1+g2+g3+g4+g5+g6)).and.  (g6 == 1)) then
          dz = -1
          call bonddirection(5,-3,8,tbo(maxlength))
       end if

       tempcoord(maxlength)%x = modulo(protcoords(m,maxlength)%x + dx-1,gridsize)+1
       tempcoord(maxlength)%y = modulo(protcoords(m,maxlength)%y + dy-1,gridsize)+1
       tempcoord(maxlength)%z = modulo(protcoords(m,maxlength)%z + dz-1,gridsize)+1

       deltaenergy = 0.0d0
       do pr =1,nprotein    
          do st = 1,maxlength,1
             if((pr /=m) .or. (st >1) .or. (st < maxlength))then
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
          call adjacent(m,maxlength,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.) .and. (bdir == -1*tbo(maxlength)) &
               .and. (bdir ==bonddd(m,st))) then
             deltaenergy = deltaenergy + intraenergy
             tempisobond(st-1,maxlength) = 1
          else
             tempisobond(st-1,maxlength) = 0
          end if
       end do

       do pr = 1,m-1,1
          do st = 1,maxlength,1
             if(intbond(pr,m,st,1) == 1) deltaenergy = deltaenergy - interenergy
             call adjacent(pr,maxlength,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and. (bdir == -1*tbo(maxlength)) &
               .and. (bdir ==bonddd(pr,st))) then
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
             call adjacent(pr,maxlength,st,tempcoord,adjver,bdir)
             if((adjver.eqv. .true.) .and. (bdir == -1*tbo(maxlength)) &
               .and. (bdir ==bonddd(pr,st))) then
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
             bonddd(m,l) = bonddd(m,l+1)
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
          call updatepos(m,maxlength,maxlength,tempcoord,tbo)
          call updateintrabond(m,maxlength,maxlength,1,maxlength-3,tempisobond)
          call updateinterbond(m,maxlength,maxlength,tempintbond)
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
     if((chdir == ndir .and. xn ==1)) tempbonddd = -1
    
  end subroutine bonddirection

  subroutine comfind
    integer :: dcomx,dcomy,dcomz
    double precision :: comx,comy,comz
    integer :: m,l
    t = time +1
    do m = 1,nprotein,1
       comx = 0.0d0
       comy = 0.0d0
       comz = 0.0d0
          if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
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

  subroutine clustercom(m,clno,clsize)
    integer,intent(in)::m
    integer,dimension(:),intent(inout) ::clno
    integer ::  a,clsize,b,base
    double precision :: delx,dely,delz,comx,comy,comz
    logical :: compass
    t = time+1
    !write(6,*) 'com start'
    compass = .false.
       comcluster%x = 0
       comcluster%y = 0
       comcluster%z = 0
       clsize = 1

    comx = 0
    comy = 0
    comz = 0
          
    do a = 1,nprotein
       if(clno(a) == clno(m)) then
          if(clsize > 1) then
          delx = com(t,a)%x - com(t,b)%x
          if (delx > gridsize/2) delx = delx-gridsize
          if (delx < -1.5) delx = delx + gridsize
          comx = comx + delx
          dely = com(t,a)%y -com(t,b)%y
          if (dely > gridsize/2) dely = dely-gridsize
          if (dely < -1.5) dely = dely + gridsize
          comy = comy + dely
          delz = com(t,a)%z-com(t,b)%z
          if (delz > gridsize/2) delz = delz-gridsize
          if (delz < -1.5) delz = delz + gridsize
          comz = comz + delz
       end if
       clsize = clsize + 1
       b = a
    end if
    if(clsize == 1) base = a
       end do
       comcluster%x = int(modulo(com(t,base)%x + comx/2 -1,real(gridsize)))+1
       comcluster%y = int(modulo(com(t,base)%y + comy/2 -1,real(gridsize)))+1
       comcluster%z = int(modulo(com(t,base)%z + comz/2 -1,real(gridsize)))+1

       !if(clsize ==1) then

       !comcluster%x = int(com(t,m)%x)
       !comcluster%y = int(com(t,m)%y)
       !comcluster%z = int(com(t,m)%z)
       !compass = .true.
       !goto 53
    !end if
       
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
    integer :: dx,dy,dz,m,l,f,g,bdir,delx,dely,delz,maxlengthss
    double precision :: initialenergy
    Type(protein),dimension(:),allocatable :: tempcoord
    logical :: adjver
    initialenergy = 0.0d0
    totalenergy = 0.0d0
    allocate(tempcoord(maxlength))
    do m= 1,nprotein,1
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do l = 1,maxlength-3,1
          do f = l+3,maxlength,1

        
             if(f < l-2 .or. f > l + 2) then
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z

                
                call adjacent(m,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.) .and. (bdir == bonddd(m,f)) .and. &
                     (bdir == -1*bonddd(m,l))) then
                   isobond(m,l,f) = 1
                   initialenergy = initialenergy + intraenergy
                else
                   isobond(m,l,f) = 0
                end if
             end if
             !write(6,*) isobond(m,l,f)
          end do
       end do
    end do

!!!!!!! interatomic
    do m= 1,nprotein,1
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
          do g = m,nprotein,1
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength1
          if(protcoords(g,1)%species == 1) maxlengthss = maxlength2
          do l = 1,maxlength,1
             do f = 1,maxlengthss,1
                if(g /=m) then

                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z



                call adjacent(g,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.) .and. bdir == bonddd(m,f) .and. &
                     bdir == -1*bonddd(m,l)) then
                   intbond(m,g,l,f) = 1
                   initialenergy = initialenergy + interenergy
                else
                   intbond(m,g,l,f) = 0
                end if
                
end if
             end do
          end do
       end do
    end do
    write(93,*) 'manual chaeck' ,initialenergy
    totalenergy = initialenergy
    write(6,*) totalenergy

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

  subroutine adjacent(pr,ch1,ch2,tempory,adjver,bdir)
    integer,intent(in):: pr,ch1,ch2
    integer,intent(inout)::bdir
    integer::dx,dy,dz,delx,dely,delz,signs
    logical,intent(inout)::adjver
    Type(protein),dimension(:),intent(in) :: tempory
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

 31   if(dx + dy + dz ==1) then
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
    write(67,*) (nprotein1*maxlength1) + (nprotein2*maxlength2)
    write(67,*) ' '
    t = time 
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          write(67,*) m, 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
               2*protcoords(m,l)%z !,bonddd(m,l)
       end do
    end do
  end subroutine dataout


  subroutine clusterenergy(clno,tempcluscoord,energyofcluster,chnum,ctbo)
    integer,dimension(:),intent(inout) ::clno
    integer,dimension(:,:),intent(inout)::ctbo
    integer,intent(in) :: chnum
    type(protein),dimension(:,:),intent(in) :: tempcluscoord
    type(protein),dimension(:),allocatable :: tempcoord
    integer:: a,l,g,f,bdir,maxlengthss
    doubleprecision,intent(inout) :: energyofcluster
    logical :: clusmove,adjver
    allocate(tempcoord(maxlength))

    energyofcluster = 0.0

    do a = 1,nprotein
          if(protcoords(a,1)%species == 1) maxlength = maxlength1
          if(protcoords(a,1)%species == 1) maxlength = maxlength2
       if(clno(a) == clno(chnum)) then
          do g = 1,nprotein
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength1
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength2
             clusmove = .true.
             if(clno(a) == clno(g)) goto 62
             do l = 1,maxlength,1
                tempcoord(l) = tempcluscoord(a,l)
                do f = 1,maxlengthss,1
                   call adjacent(g,l,f,tempcoord,adjver,bdir)
                   if(adjver.eqv. .true. .and. bdir == -1*ctbo(a,l) .and. &
                        bdir == ctbo(g,f)) energyofcluster = energyofcluster + interenergy
                end do
             end do
62           if(clusmove .eqv. .true.) continue
          end do
       end if
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
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
       if(protcoords(m,1)%species == 1) maxlength = maxlength2
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

subroutine clusrbond(m,ctbo,axis,rotor)
  integer,intent(in):: m,axis,rotor
  integer,dimension(:,:),intent(inout)::ctbo
  logical :: bm
  integer :: sign,s


  do s = 1,maxlength
bm = .true.
  if(bonddd(m,s) == axis) then
     bm = .true.
        ctbo(m,s) = bonddd(m,s)
     goto 87
  end if


    if(rotor == 0) then
       bm = .true.
           ctbo(m,s) = -1*bonddd(m,s)
      goto 87
  end if

  sign = bonddd(m,s)/abs(bonddd(m,s))
  if(axis == 1) then
     if(rotor == -1) then
        if(abs(bonddd(m,s)) == 2) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 1)
        if(abs(bonddd(m,s)) == 3) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 1)
        goto 87
     else if(rotor == 1) then
        if(abs(bonddd(m,s)) == 2) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 1)
        if(abs(bonddd(m,s)) == 3) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 1)
        goto 87
     end if
  else if(axis ==2 ) then
     if(rotor == -1) then
        if(abs(bonddd(m,s)) == 1) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 2)
        if(abs(bonddd(m,s)) == 3) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 2)        
        goto 87
     else if(rotor == 1) then
        if(abs(bonddd(m,s)) == 1) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 2)
        if(abs(bonddd(m,s)) == 3) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 2)         
goto 87
     end if

  else if(axis ==3) then
             if(rotor == -1) then
        if(abs(bonddd(m,s)) == 1) ctbo(m,s) = sign*(abs(bonddd(m,s)) + 1)
        if(abs(bonddd(m,s)) == 2) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) - 1)        
        goto 87
     else if(rotor == 1) then
        if(abs(bonddd(m,s)) == 1) ctbo(m,s) = -1*sign*(abs(bonddd(m,s)) + 1)
        if(abs(bonddd(m,s)) == 2) ctbo(m,s) = sign*(abs(bonddd(m,s)) - 1)         
     end if

     end if

87   if(bm .eqv. .true.) continue
  end do
end subroutine clusrbond
  
  subroutine clustertranslation(m,clno,clsize)
    integer,intent(in)::m
    integer,intent(inout):: clsize
    integer,dimension(:),intent(inout) ::clno
    integer :: steplength,dx,dy,dz
    type(protein),dimension(:,:),allocatable :: tempcluscoord
    integer,dimension(:,:),allocatable::ctbo
    double precision :: clusen,direction
        logical :: moveallow
    allocate(tempcluscoord(nprotein,maxlength))
    allocate(ctbo(nprotein,maxlength))

 
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

    call clustermove(clno,m,dx,dy,dz,tempcluscoord,moveallow,ctbo)
    if(moveallow .eqv. .true.) then
    clusen = 0.0
    call clusterenergy(clno,tempcluscoord,clusen,m,ctbo)
    if(clusen == 0.0)  call updatecluspos(clno,tempcluscoord,m,ctbo)
    end if
  end subroutine clustertranslation


  subroutine clusterrotation(m,clno,clsize)
    integer,intent(in) :: m
    integer,intent(inout)::clsize
    integer,dimension(:),intent(inout) ::clno
    integer :: a,b,l,pr,st
    double precision :: rotatechoose,clusen
    integer,dimension(:,:),allocatable :: deltax,deltay,deltaz
    Type(protein),dimension(:,:),allocatable :: tempcluscoord
    integer,dimension(:,:),allocatable::ctbo
    logical :: clusmove
    clusmove = .true.
    allocate(ctbo(nprotein,maxlength))
    allocate(tempcluscoord(nprotein,maxlength))
    t = time + 1
    !m = int(ran2(seed)*(clusterno))+1
    !write(6,*) 'pre com'
    call clustercom(m,clno,clsize)
    !write(6,*) 'compass'
    allocate(deltax(nprotein,maxlength))
    allocate(deltay(nprotein,maxlength))
    allocate(deltaz(nprotein,maxlength))

    do a = 1,nprotein
       if(clno(a) == clno(m)) then
          b = clno(m)
             if(protcoords(b,1)%species == 1) maxlength = maxlength1
             if(protcoords(b,1)%species == 1) maxlength = maxlength2
       do l = 1, maxlength
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


    rotatechoose = ran2(seed)

    if(rotatechoose <= 1.0/9) then
       call rotatex(b,1,-1,deltay,deltaz,tempcluscoord,ctbo)
       else if(rotatechoose > 1.0/9 .and. rotatechoose <= 2.0/9) then
          call rotatex(b,-1,1,deltay,deltaz,tempcluscoord,ctbo)
       else if(rotatechoose > 2.0/9 .and. rotatechoose <= 1.0/3) then
          call clusterflip(b,1,-1,1,deltax,deltay,deltaz,tempcluscoord,ctbo)    
       else if(rotatechoose > 1.0/3 .and. rotatechoose <= 4.0/9) then
          call rotatey(b,1,-1,deltax,deltaz,tempcluscoord,ctbo)
       else if(rotatechoose > 4.0/9 .and. rotatechoose <= 5.0/9) then
          call rotatey(b,-1,1,deltax,deltaz,tempcluscoord,ctbo)
       else if(rotatechoose > 5.0/9 .and. rotatechoose <= 6.0/3) then
          call clusterflip(b,1,1,-1,deltax,deltay,deltaz,tempcluscoord,ctbo)
       else if(rotatechoose>6.0/9 .and. rotatechoose <= 7.0/9) then
          call rotatez(b,1,-1,deltax,deltay,tempcluscoord,ctbo)
       else if(rotatechoose > 7.0/9 .and. rotatechoose <= 8.0/9) then
          call rotatez(b,-1,1,deltax,deltay,tempcluscoord,ctbo)
       else if(rotatechoose > 8.0/9 .and. rotatechoose <= 1.0) then
          call clusterflip(b,-1,1,1,deltax,deltay,deltaz,tempcluscoord,ctbo)
end if
       end if
    end do

       do pr = 1,nprotein
          if(clno(pr) == clno(m)) goto 39
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


    clusen = 0.0
    call clusterenergy(clno,tempcluscoord,clusen,m,ctbo)
    if(clusen == 0.0)  call updatecluspos(clno,tempcluscoord,m,ctbo)


31  if(clusmove .eqv. .false.) continue
    !calculate the energy penalty and insure no overlaps

  end subroutine clusterrotation

  subroutine clusterflip(protno,cx,cy,cz,deltax,deltay,deltaz,tempcluscoord,ctbo)
    integer,intent(in)::protno,cx,cy,cz
        integer,dimension(:,:),intent(inout)::ctbo
    !integer,dimension(:),intent(in)::clno
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%x = modulo(comcluster%x + cx*deltax(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(comcluster%y + cy*deltay(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z + cz*deltaz(protno,l) -1,gridsize)+1
    end do

if(cx == -1) call clusrbond(protno,ctbo,1,0)    
if(cy == -1) call clusrbond(protno,ctbo,2,0)
if(cz == -1) call clusrbond(protno,ctbo,3,0)   


end subroutine clusterflip

  subroutine rotatex(protno,cy,cz,deltay,deltaz,tempcluscoord,ctbo)
    integer,intent(in)::protno,cy,cz
        integer,dimension(:,:),intent(inout)::ctbo
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltay,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%x = protcoords(protno,l)%x
       tempcluscoord(protno,l)%y = modulo(comcluster%y +cy*deltaz(protno,l) -1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z +cz*deltay(protno,l) -1,gridsize)+1
    end do

    if(cy ==-1) call clusrbond(protno,ctbo,1,1)
    if(cz == -1) call clusrbond(protno,ctbo,1,-1)

  end subroutine rotatex

  subroutine rotatey(protno,cx,cz,deltax,deltaz,tempcluscoord,ctbo)
    integer,intent(in)::protno,cx,cz
    integer,dimension(:,:),intent(inout)::ctbo
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltaz
    integer :: l

    do l = 1,maxlength
       tempcluscoord(protno,l)%y = protcoords(protno,l)%y
       tempcluscoord(protno,l)%x = modulo(comcluster%x +cx*deltaz(protno,l)-1,gridsize)+1
       tempcluscoord(protno,l)%z = modulo(comcluster%z +cz*deltax(protno,l)-1,gridsize)+1
    end do

        if(cx ==-1) call clusrbond(protno,ctbo,2,1)
    if(cz == -1) call clusrbond(protno,ctbo,2,-1)
    

  end subroutine rotatey

  subroutine rotatez(protno,cx,cy,deltax,deltay,tempcluscoord,ctbo)
    integer,intent(in)::protno,cx,cy
        integer,dimension(:,:),intent(inout)::ctbo
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:,:),intent(in) :: deltax,deltay
    integer :: l

       do l = 1,maxlength
       tempcluscoord(protno,l)%z = protcoords(protno,l)%z
       tempcluscoord(protno,l)%x = modulo(comcluster%x +cx*deltay(protno,l)-1,gridsize)+1
       tempcluscoord(protno,l)%y = modulo(comcluster%y +cy*deltax(protno,l)-1,gridsize)+1
    end do

        if(cx ==-1) call clusrbond(protno,ctbo,3,-1)
    if(cy == -1) call clusrbond(protno,ctbo,3,1)

  end subroutine rotatez


  subroutine clustermove(clno,chainno,delx,dely,delz,tempcluscoord,moveallow,ctbo)
    integer,intent(in) :: chainno,delx,dely,delz
    integer,dimension(:,:),intent(inout)::ctbo
    logical,intent(inout)::moveallow
    Type(protein),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::clno
    integer :: b,pr,st
    logical :: clusmove

    clusmove = .true.
    moveallow = .true.
    if(protcoords(chainno,1)%species == 1) maxlength = maxlength1
    if(protcoords(chainno,1)%species == 1) maxlength = maxlength2
    do b =1,maxlength
       tempcluscoord(chainno,b)%x = modulo(protcoords(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(protcoords(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(protcoords(chainno,b)%z + delz-1,gridsize)+1
       ctbo(chainno,b) = bonddd(chainno,b)

    do pr = 1,chainno
       if(clno(pr) == clno(chainno)) goto 39
       do st = 1,maxlength,1
          if(overlapsclus(chainno,pr,b,st,tempcluscoord) .eqv. .false.) then
             clusmove = .false.
             goto 31
          end if
       end do
39     if(clusmove .eqv. .true.) continue
    end do

end do    
   31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove

  subroutine debug
    integer::m,l,f,g,dx,dy,dz,maxlengthss

    do m = 1,nprotein,1
          if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do l = 1,maxlength,1
          do f = 1,nprotein,1
             if(protcoords(f,1)%species == 1) maxlengthss = maxlength1
             if(protcoords(f,1)%species == 1) maxlengthss = maxlength2
             do g = 1,maxlengthss,1
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

               if(bonddd(m,l) == 0) write(6,*) 'bond vector fail'
       end do
    end do


    
  end subroutine debug

  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver
    type(protein),dimension(:),allocatable :: tempcoord
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    clcount = 0
    clusterpop = 0
    do m = 1,nprotein
       clnos(m) = m
    end do


    do m = 1,nprotein-1
       if(protcoords(m,1)%species == 1) maxlength = maxlength1
          if(protcoords(m,1)%species == 1) maxlength = maxlength2
       do g = m,nprotein
          if((g/=m) .and. (clnos(m) /= clnos(g))) then
             if(protcoords(g,1)%species == 1) maxlengthss = maxlength1
          if(protcoords(g,1)%species == 1) maxlengthss = maxlength2
             do l = 1,maxlength
                tempcoord(l)%x = protcoords(m,l)%x
                tempcoord(l)%y = protcoords(m,l)%y
                tempcoord(l)%z = protcoords(m,l)%z
                do f = 1,maxlengthss
                   call adjacent(g,l,f,tempcoord,adjver,bdir)
                   if(adjver.eqv. .true. .and. bdir == bonddd(g,f) .and. &
                        bdir == -1*bonddd(m,l)) then
                      clnos(g) = clnos(m)
                      !clusterpop = clusterpop + 1
                      clusyes = .true.
                      goto 13
                   end if
                end do
             end do
             clusyes = .true.
13           if(clusyes .eqv. .true.) continue
          end if
       end do
    end do
    
    do m = 1,nprotein
       if(clnos(m)==m) clcount = clcount + 1
       if(clnos(m) /=m) clusterpop = clusterpop + 1
    end do
    
    clusterpop = clusterpop + clcount
    write(82,*) time,clcount, clusterpop, (clusterpop/clcount)
  end subroutine clustercount



SUBROUTINE read_setup

   USE keywords
   !USE setup_vars
   !USE stats_vars
   IMPLICIT NONE

   CHARACTER(LEN=20) :: keyword, option, argument
   CHARACTER(LEN=18), PARAMETER :: param_fmt1='(A16, 1PE20.10, A)'
   CHARACTER(LEN=16), PARAMETER :: param_fmt2='(A16, F16.10)'
   INTEGER :: err, i, j
   LOGICAL :: success

   double precision :: intraenergy,interenergy

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
   
   !traj_file = 'trajectory'       ! Stem of name for trajectory files.
   !vgen = 1                       ! 1=random initial velocities (no momentum), 2=Maxwell-Boltzmann
   !v_init_file = ' '              ! Name of starting velocities file.
   !v_rescale = .FALSE.            ! If reading velocities from a file, .T.=rescale to get specified K.E.
   !wellstats_file = 'wellstats'   ! Stem of name for individual well statistics files.
   !well_tol = 1.0D-6              ! Energy tolerance for quenches to be considered the same.
   !xmol_type(1) = ''              ! Atom type for A (or charged) atoms in xmol dumps.
   !xmol_type(2) = ''              ! Atom type for B (or uncharged) atoms in xmol dumps.

   OPEN (UNIT=20, FILE='setmultup.txt', STATUS='OLD', IOSTAT=err)

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
         
      CASE ('CHAIN_LENGTH')
         CALL get_integer(maxlength1)

               CASE('CHAIN_2_LENGTH')
         call get_integer(maxlength2)

      CASE ('NUM_CHAINS')
         CALL get_integer(nprotein1)

      CASE ('NUM_CHAINS_2')
         CALL get_integer(nprotein2)

      CASE ('LATTICE_DIMENSIONS')
         CALL get_integer(gridsize)
                  
      CASE ('SIM_TIME')
         CALL get_integer(maxtime)
         
      CASE ('EQUIB_TIME')
         CALL get_integer(equilib)

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

      CASE ('MAXPIV')
         CALL get_integer(maxpiv)

      CASE ('CRANKLIMIT')
         CALL get_integer(cranklimit)
         
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
