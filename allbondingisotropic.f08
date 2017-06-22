program move
  implicit none


  integer :: maxno,protno,seed,natoms,maxsize,maxtime,reptbackward,reptforward,equilib,st,fi
  real :: i,j,k,gridsize,randomz,deltax1,deltax2,dx,dy,dz,msd
  integer :: l,m,hold,z1,z2,z3,z4,z5,z6,control,debugging
  real :: deltay1,deltay2,deltaz1,deltaz2,totrmsbrute
  integer :: nprotein,split,section,time,count
  integer :: maxlength,g,xl,yl,zl,t,r,N
  real, external :: ran2
  logical :: rightangle,endchain,run,pivotx,pivoty,pivotz,run2,rac,exist,fail,finalfail
  real :: totdisp,totchainlength,avechainlength,probs,comx,comy,comz,totrmserror
  real :: actualchainlength,totvar,deviate,totdeviate,varian,separation,lengthvar,random
  real,dimension(:),allocatable :: variance,disp,rmserror,dxtot,dytot,dztot,dx1,dy1,dz1
  integer :: gate,xbk,ybk,zbk,protpass,f,scan,c
  Real:: t1,t2,t3,t4,t5,t6,di,dj,dk
  integer :: cont1,cont2,cont3,cont4,totcont,successful,reject
  character(len = 10) :: commandread,commandread2,comread3,comread4,comread5,comread6,comread7,comread8,comread9,comread10,comread11
  real :: runningaveROG,runningaveEtE,chainlength,sumdebug
  !real :: testlength,testtotlength
  integer :: piv,right,end,crank,rept,isoenergy,deltaiso,isocount,datayes,deltaint
  real,dimension(:,:,:), allocatable :: isobond,tempisobond
  real,dimension(:,:,:,:),allocatable :: tempintbond,intbond
  real :: kT
  character(len = 10) :: comread12
  logical :: energypassx,energypassy,energypassz,overlapvar,reptcont
  integer :: pr,ch1,ch2,str,fir,proti,othere


  type protein
     real :: x,y,z
     integer :: location
  end type protein
  type(protein),dimension(:,:),allocatable :: protcoords
  type(protein),dimension(:,:),allocatable :: com
  type(protein),dimension(:,:),allocatable :: tempcoord

  open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')

  open(29, file = 'rms.dat', action = 'write')
  !open(31, file = 'chainlength.dat', action = 'write')
  !open(51, file = 'bonds.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  !open(77, file = 'logmsd.dat', action = 'write')
  open(79, file = 'runningave.dat', action = 'write')
  !open(99, file = 'deltas.dat', action = 'write')
  !open(93, file = 'temporarycoord.dat', action = 'write')
  open(93, file = 'energyiso.dat', action = 'write')
  !open(72, file = 'bonddata.dat', action = 'write')


  !sort out end moves, pivot, right and crank
  reptforward = 0
  reptbackward = 0
  totrmsbrute = 0.0
  runningaveROG = 0.0
  runningaveEtE = 0.0
  finalfail = .false.


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

  !reads in seed from commandline


  if (datayes == 1) open(67, file = 'move.xyz', action = 'write')
  N = maxlength
  allocate(protcoords(nprotein,maxlength))
  allocate(com(maxtime+1,nprotein))
  allocate(tempcoord(nprotein,maxlength))
  !allocate(chainlength(nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(dxtot(nprotein))
  allocate(dytot(nprotein))
  allocate(dztot(nprotein))
  allocate(dx1(maxlength))
  allocate(dy1(maxlength))
  allocate(dz1(maxlength))
  allocate(isobond(nprotein,maxlength,maxlength))
  allocate(tempisobond(nprotein,maxlength,maxlength))
  allocate(intbond(nprotein,nprotein,maxlength,maxlength))
  allocate(tempintbond(nprotein,nprotein,maxlength,maxlength))

  count = 0
  time = 0
  totdisp = 0.0

  call foundation

  !call dataout
  call comfind
  actualchainlength = 0.0
  do time = 1,maxtime,1
     fail = .false.
     if (mod(time,1) == 0 .and. datayes == 1) then
        call dataout
     end if
     call positioning
     call comfind

     if (time > equilib) then
        call rms
     end if

     if(time >equilib .and. modulo(time,100) == 0) then
        call length
        call radiusofgyration
        write(79,*) time-equilib, runningaveEtE/(time-equilib), runningaveROG/(time-equilib)
     end if
     if (debugging == 1) then
        call debug
     else
        continue
     end if

     !write(6,*) 'g'
     !if (fail .eqv. .true.) then
     !write(6,*) 'step= ', time, 'FAIL'
     !else
     !write(6,*) 'step =',time
     !end if

  end do
  !call error
  write(6,*) 'average chain length squared is;' ,(actualchainlength/maxtime)**2, actualchainlength/maxtime
  write(6,*) 2*log(actualchainlength/maxtime)/log(real(maxlength))
  inquire(file = "lengthdataaverages.dat", exist = exist)
  if (exist) then
     open(19, file = "lengthdataaverages.dat", status = "old", position = "append",action ="write")
  else
     open(19, file = "lengthdataaverages.dat", status = "new",action ="write")
  end if
  write(19,*) 2*log(actualchainlength/maxtime)/log(real(maxlength)),2*(SQRT(totvar)/maxtime)/(actualchainlength/maxtime)
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
    read(17,*) BIN
    read(17,*) split
  end subroutine initial

  subroutine foundation
    !read in initial positions of the protein

    character(len = 10) :: BIN
    read(23,*) BIN
    write(6,*) maxlength
    !write(67,*) nprotein*maxlength
    !write(67,*) ' '
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          read(23,*) BIN, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
          !write(67,*) 'C', 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
          !2*protcoords(m,l)%z
          !write(6,*) protcoords(m,l)%x
       end do
    end do
    isocount = 0
    write(6,*) 'a'
    call energy
    write(6,*) 'b'
    isoenergy = isocount
    write(93,*) 0, isoenergy

  end subroutine foundation




  subroutine positioning
    !randomly selects one of the 4 moves

    integer :: contchoose,a1,a2,a3,a4,nmoves


    t = time + 1
    do scan = 1,nprotein*maxlength,1
       count = count + 1
       rac = .true.
       run = .true.
       run2 = .true. 

       overlapvar = .false.
       reptcont = .false.
       randomz = ran2(seed)
       nmoves = piv + crank + end + right + rept

       if (randomz <= (1.0*end)/ nmoves .and. end == 1) then
          !write(6,*) 'end' 
          call endmove

       else if (randomz > (1.0*end)/nmoves .and. randomz <= (1.0*end+crank)/nmoves .and. crank == 1) then
          ! write(6,*) 'crank' 
          call crankshaftmove

       else if (randomz > (1.0*end + crank)/nmoves .and. randomz <= (1.0*end + crank + right)/nmoves .and. right ==1) then
          ! write(6,*) 'right' 
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


  subroutine crankshaftmove
    !performs crankshaft move
    logical :: crankcont,cranksep
    real :: dx11,dx12,dy11,dy12,dz11,dz12,dummy,dummy2
    integer :: p,s
    !this means that m/=0 and is up to nprotein
    m =int(ran2(seed)*(nprotein-1))+1
    !this means that l cannot be maxlength,maxlength - 1 or l
    cranksep = .true.
71  if(cranksep .eqv. .true.) then
       continue
    end if

    proti = m

    dummy = int(ran2(seed)*(maxlength-1))+1
    dummy2 = int(ran2(seed)*(maxlength-1))+1
    l = min(dummy,dummy2)
    p = max(dummy,dummy2)
    if (p-l <= 3) then
       goto 71
    end if

    probs = ran2(seed)
    pivotx = .false.
    pivoty = .false.
    pivotz = .false.


    str = l+1
    fir = p-1
    st = 1
    fi = maxlength

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

    !pivotz = .false.
    !pivoty = .false.

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
             tempcoord(m,s)%x = protcoords(m,s)%x
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y-dz1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z+dy1(s)-1,gridsize)+1
          end do


          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if







          !xl = int((protcoords(t,m,l)%x*10)/split)-1
          !yl = 1+int(((protcoords(t,m,l)%y*10)/split) -1)*N
          !zl = int(((protcoords(t,m,l)%z*10)/split) -1)*(N**2)
          !section = xl + yl + zl
          !protcoords(t,m,l)%location = section

       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%x = protcoords(m,s)%x
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y+dz1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z-dy1(s)-1,gridsize)+1
          end do


          call overlap

          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if


       else if (probs > (2.0/3) .and. probs <= (1.0)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%x = protcoords(m,s)%x
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y-dy1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z-dz1(s)-1,gridsize)+1
          end do


          call overlap

          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if
       end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
             tempcoord(m,s)%y = protcoords(m,s)%y
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x-dz1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z+dx1(s)-1,gridsize)+1
          end do

          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if




       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%y = protcoords(m,s)%y
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x+dz1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z-dx1(s)-1,gridsize)+1
          end do


          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if


       else if (probs > (2.0/3) .and. probs <= (1.0)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%y = protcoords(m,s)%y
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x-dx1(s)-1,gridsize)+1
             tempcoord(m,s)%z = modulo(protcoords(m,l)%z-dz1(s)-1,gridsize)+1
          end do

          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if

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
             tempcoord(m,s)%z = protcoords(m,s)%z
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y-dx1(s)-1,gridsize)+1
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x+dy1(s)-1,gridsize)+1
          end do
          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if


       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%z = protcoords(m,s)%z
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y+dx1(s)-1,gridsize)+1
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x-dy1(s)-1,gridsize)+1
          end do
          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if



       else if (probs > (2.0/3) .and. probs <= (1.0)) then

          do s = l+1,p-1,1
             tempcoord(m,s)%y = modulo(protcoords(m,l)%y-dy1(s)-1,gridsize)+1
             tempcoord(m,s)%x = modulo(protcoords(m,l)%x-dx1(s)-1,gridsize)+1
             tempcoord(m,s)%z = protcoords(m,s)%z
          end do

          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             crankcont = .false.
             goto 43
          end if

       end if
    end if

    !if (pivotx .eqv. .false. .and. pivoty .eqv. .false. .and. pivotz .eqv. .false.) then
    !crankcont = .false.
    !end if

    deltaiso = 0

    str = l+1
    fir = p-1
    st = 1
    fi = maxlength
    !write(6,*) 'start crank'
    !call deltae

    write(6,*) 'deltaiso',deltaiso

98  if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then
 
       do s = l+1,p-1,1
          protcoords(m,s)%x = tempcoord(m,s)%x
          protcoords(m,s)%y = tempcoord(m,s)%y
          protcoords(m,s)%z = tempcoord(m,s)%z
          do g = 1,maxlength,1
             if(g <= l .or. g >= p) then
                isobond(m,s,g) = tempisobond(m,s,g)
             end if
          end do

       end do
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

    !performs a 180 flip on an atoms

    m =int(ran2(seed)*(nprotein-1))+1
    ! l cannot be equal to 1 or maxlength
    l = int(ran2(seed)*(maxlength-3))+2
    !write(6,*) 'm =', m, 'l = ', l

    rac = .true.

    if (abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= 2.0 &
         .and. abs(mod(protcoords(m,l+1)%x - protcoords(m,l-1)%x,gridsize)) /= gridsize-2.0 &
         .and. abs(mod(protcoords(m,l+1)%y - protcoords(m,l-1)%y,gridsize)) /= gridsize-2.0 &
         .and. abs(mod(protcoords(m,l+1)%z - protcoords(m,l-1)%z,gridsize)) /= gridsize-2.0) then
       rightangle = .true.
    else 
       rightangle = .false.
    end if

    if (rightangle .eqv. .true.) then    

       deltax1 = protcoords(m,l+1)%x - protcoords(m,l)%x
       deltax2 =protcoords(m,l-1)%x - protcoords(m,l)%x
       deltay1=protcoords(m,l+1)%y - protcoords(m,l)%y
       deltay2 =protcoords(m,l-1)%y - protcoords(m,l)%y
       deltaz1 = protcoords(m,l+1)%z - protcoords(m,l)%z
       deltaz2= protcoords(m,l-1)%z - protcoords(m,l)%z

       di = 0
       dj = 0
       dk = 0

       di = (deltax1 + deltax2)
       dj = (deltay1 + deltay2)
       dk = (deltaz1 + deltaz2)



       if (abs(di) <= 1 .or. abs(dj) <= 1 .or. abs(dk) <= 1) then

          tempcoord(m,l)%x = modulo(protcoords(m,l)%x+di-1,gridsize)+1
          tempcoord(m,l)%y = modulo(protcoords(m,l)%y+dj-1,gridsize)+1
          tempcoord(m,l)%z = modulo(protcoords(m,l)%z+dk-1,gridsize)+1

          proti = m





          deltaiso = 0
          str = l
          fir = l
          st = 1
          fi = maxlength
          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             rac = .false.
             goto 37
          end if



98       if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then
             protcoords(m,l)%x = tempcoord(m,l)%x
             protcoords(m,l)%y = tempcoord(m,l)%y
             protcoords(m,l)%z = tempcoord(m,l)%z

             do g = 1,maxlength,1
                if(g /= l) then
                   isobond(m,l,g) = tempisobond(m,l,g)
                end if
             end do

             !    write(6,*) time, 'rigth'
             successful = successful + 1
          else
             rac = .false.
             goto 37
          end if
       end if
       end if

37  if (rac .eqv. .false. .or. rightangle .eqv. .false.) then

       reject = reject + 1
    end if


  end subroutine rightanglemove

  subroutine pivot
    integer :: b
    logical :: pivcont
    real :: choose3,xhold,yhold,zhold
    real,dimension(:),allocatable :: delx,dely,delz
    allocate(delx(maxlength))
    allocate(dely(maxlength))   
    allocate(delz(maxlength))

    !write(6,*) 'here comes the pivot!!'

    choose3 = ran2(seed)
    pivcont = .true. 
    m =int(ran2(seed)*(nprotein-1))+1

    l = int(ran2(seed)*(maxlength-3))+2
    !write(6,*) 'm = ', m
    !write(6,*) 'l =', l
    t = time + 1

    proti = m


    if(l > maxlength/2) then
       !continue
       !else if (l == 1000) then

       str = l+1
       fir = maxlength
       st = 1
       fi = l


       do g = l+1,maxlength,1

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
          !write(99,*) delx(g),dely(g),delz(g) !,xhold,yhold,zhold
       end do
       if (choose3 <= 1.0/5) then

          do b = l+1,maxlength,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x - dely(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y + delx(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z
          end do
          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if


       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then
          do b =l+1,maxlength,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x + dely(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y - delx(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z                    
          end do

          call overlap

          if(overlapvar .eqv. .true.) then

             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if


       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then

          do b =l+1,maxlength,1                 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x - delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y - dely(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z
          end do
          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if


       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then

          do b =l+1,maxlength,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x - delz(b)-1,gridsize)+1
             tempcoord(m,b)%z = modulo(protcoords(m,l)%z + delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = protcoords(m,b)%y    
          end do

          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if

       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then

          do b = l+1,maxlength,1                
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x + delz(b)-1,gridsize)+1
             tempcoord(m,b)%z = modulo(protcoords(m,l)%z - delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = protcoords(m,b)%y
          end do

          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 98
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if

       end if




       write(6,*) 'deltaiso', deltaiso


98     if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then

          do b = l+1,maxlength,1 
             protcoords(m,b)%x = tempcoord(m,b)%x
             protcoords(m,b)%y = tempcoord(m,b)%y
             protcoords(m,b)%z = tempcoord(m,b)%z

             do g = 1,l,1

                isobond(m,b,g) = tempisobond(m,b,g)
             end do

          end do
          !          write(6,*) time, 'pivot1'
          successful = successful + 1

       else
          pivcont = .false.
          goto 75
       end if

    else if (l <= maxlength/2) then
       !continue
       !else if (l == 100000) then
       str = 1
       fir = l-1
       st = l
       fi =maxlength
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

          !write(99,*) delx(g),dely(g),delz(g),xhold,yhold,zhold

       end do
       if (choose3 <= 1.0/5) then


          do b = 1,l-1,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x -dely(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y + delx(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z
          end do

          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 92
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if


       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then 
          do b =1,l-1,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x + dely(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y - delx(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z
          end do


          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 92
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if



       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then
          do b =1,l-1,1                 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x -delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = modulo(protcoords(m,l)%y - dely(b)-1,gridsize)+1
             tempcoord(m,b)%z = protcoords(m,b)%z
          end do

          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 92
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if


       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then
          do b =1,l-1,1 
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x - delz(b)-1,gridsize)+1
             tempcoord(m,b)%z = modulo(protcoords(m,l)%z + delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = protcoords(m,b)%y
          end do

          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 92
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if

       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then
          do b = 1,l-1,1               
             tempcoord(m,b)%x = modulo(protcoords(m,l)%x + delz(b)-1,gridsize)+1
             tempcoord(m,b)%z = modulo(protcoords(m,l)%z - delx(b)-1,gridsize)+1
             tempcoord(m,b)%y = protcoords(m,b)%y
          end do
          call overlap
          if(overlapvar .eqv. .true.) then
             call deltae
             goto 92
          else if(overlapvar .eqv. .false.) then
             pivcont = .false.
             goto 75
          end if

       end if

92     if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then
    write(6,*) 'deltaiso', deltaiso

          do b = 1,l-1,1
             protcoords(m,b)%x = tempcoord(m,b)%x
             protcoords(m,b)%y = tempcoord(m,b)%y
             protcoords(m,b)%z = tempcoord(m,b)%z

             do g = l,maxlength,1
                isobond(m,b,g) = tempisobond(m,b,g)
             end do

          end do
          successful = successful + 1
          !  write(6,*) time, 'pivot 2'
       else
          pivcont = .false.
          goto 75
       end if
    end if

75  if (pivcont .eqv. .false.) then
       reject = reject + 1
    end if


  end subroutine pivot

  subroutine reptation
    !perform reptation
    real :: choose,dxx,dyx,dzx,direc
    integer :: endchaincont,g1,g2,g3,g4,g5,g6
    logical :: energypass

    !do m = 1,nprotein
    !do l = 1,maxlength
    !do f = 1,maxlength
    !write(6,*) l,f,isobond(m,l,f)
    !end do
    !end do
    !end do

    t = time + 1
    choose = ran2(seed) - 0.5
    !write(6,*) 'choose', choose
    reptcont = .true.
    m =int(ran2(seed)*(nprotein-1))+1

    proti = m

    g1 = 1
    g2 = 1
    g3 = 1
    g4 = 1
    g5 = 1
    g6 = 1
    dx = 0.0
    dy = 0.0
    dz = 0.0

    if (choose > 0.0) then
       dxx = protcoords(m,2)%x - protcoords(m,1)%x
       if (dxx > 1.5) dxx = -1.0
       if(dxx <-1.5) dxx = 1.0
       dyx = protcoords(m,2)%y - protcoords(m,1)%y
       if (dyx > 1.5) dyx = -1.0
       if (dyx < -1.5) dyx = 1.0
       dzx = protcoords(m,2)%z - protcoords(m,1)%z
       if (dzx > 1.5) dzx = -1.0
       if (dzx < -1.5) dzx = 1.0

       if (dxx == 1.0) g1 = 0
       if (dxx == -1.0) g2 = 0
       if (dyx == 1.0) g3 = 0
       if (dyx == -1.0) g4 = 0
       if (dzx == 1.0) g5 =0
       if (dzx == -1.0) g6 = 0

       direc = ran2(seed)
       if(direc <= (1.0*g1)/5 .and. g1 == 1) then
          dx = 1.0
       else if(direc <= (1.0*g1 +g2)/5 .and.direc > (1.0*g1)/5 .and.  g2 == 1) then
          dx = -1.0
       else if(direc <= (1.0*g1+g2+g3)/5 .and.direc > (1.0*g1+g2)/5 .and.  g3 == 1) then
          dy = 1.0
       else if(direc <= (1.0*g1+g2+g3+g4)/5 .and.direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
          dy = -1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5)/5 .and.direc > (1.0*g1+g2+g3+g4)/5 .and. &
            g5 == 1) then
          dz = 1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5+g6)/5 .and. &
            direc > (1.0*g1+g2+g3+g4+g5)/5 .and.  g6 == 1) then
          dz = -1.0
       end if
       !write(6,*) 'forwrd',dx,dy,dz
       tempcoord(m,1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(m,1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(m,1)%z = modulo(protcoords(m,1)%z + dz-1,gridsize)+1

       str = 1
       fir = 1
       st = 4
       fi = maxlength -1
       othere = maxlength
       call overlap
       if(overlapvar .eqv. .true.) then
          call deltae
         ! write(6,*) 'deltaiso', deltaiso
          goto 92
       else if(overlapvar .eqv. .false.) then
          reptcont = .false.
          goto 83
       end if

       write(6,*) 'pre'
       !       call deltae
       write(6,*) 'post'

92     if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then

              write(6,*) 'deltaiso', deltaiso
          do l =maxlength,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z

             do g = maxlength,2,-1
                isobond(m,l,g) = isobond(m,l-1,g)
             end do
             !successful = successful+1
          end do

          !      write(6,*) 'dx =', dx, 'dy =', dy, 'dz = ', dz
          protcoords(m,1)%x = tempcoord(m,1)%x
          protcoords(m,1)%y = tempcoord(m,1)%y
          protcoords(m,1)%z = tempcoord(m,1)%z

          do g = 4,maxlength-1,1
             isobond(m,1,g) = tempisobond(m,1,g)
          end do
          !write(6,*) time, 'rept1'   
          successful = successful + 1
          reptforward = reptforward + 1
       else
          reptcont = .false.
       end if

       !end if
    else if (choose < 0.0) then


       dxx = protcoords(m,maxlength-1)%x - protcoords(m,maxlength)%x
       if (dxx > 1.5) dxx = -1.0
       if(dxx <-1.5) dxx = 1.0
       dyx = protcoords(m,maxlength-1)%y - protcoords(m,maxlength)%y
       if (dyx > 1.5) dyx = -1.0
       if (dyx < -1.5) dyx = 1.0
       dzx = protcoords(m,maxlength-1)%z - protcoords(m,maxlength)%z
       if (dzx > 1.5) dzx = -1.0
       if (dzx < -1.5) dzx = 1.0

       if (dxx == 1.0) g1 = 0
       if (dxx == -1.0) g2 = 0
       if (dyx == 1.0) g3 = 0
       if (dyx == -1.0) g4 = 0
       if (dzx == 1.0) g5 =0
       if (dzx == -1.0) g6 = 0
       !write(6,*) 'success?', g1+g2+g3+g4+g5+g6
       direc = ran2(seed)
       if(direc <= (1.0*g1)/5 .and. g1 == 1) then
          dx = 1.0
       else if(direc <= (1.0*g1 +g2)/5 .and.direc > (1.0*g1)/5 .and.  g2 == 1) then
          dx = -1.0
       else if(direc <= (1.0*g1+g2+g3)/5 .and.direc > (1.0*g1+g2)/5 .and.  g3 == 1) then
          dy = 1.0

       else if(direc <= (1.0*g1+g2+g3+g4)/5 .and. direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
          dy = -1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5)/5 .and. &
            direc > (1.0*g1+g2+g3+g4)/5 .and.  g5 == 1) then
          dz = 1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5+g6)/5 .and. &
            direc > (1.0*g1+g2+g3+g4+g5)/5 .and.  g6 == 1) then
          dz = -1.0
       end if
       !write(6,*) 'backward',dx,dy,dz
       tempcoord(m,maxlength)%x = modulo(protcoords(m,maxlength)%x + dx-1,gridsize)+1
       tempcoord(m,maxlength)%y = modulo(protcoords(m,maxlength)%y + dy-1,gridsize)+1
       tempcoord(m,maxlength)%z = modulo(protcoords(m,maxlength)%z + dz-1,gridsize)+1


       str = maxlength
       fir = maxlength
       st = 2
       fi = maxlength -3
       othere = 1
       call overlap
       !       write(6,*) 'backward'
       if(overlapvar .eqv. .true.) then
          !write(6,*) 'survived overlap'
          call deltae
!          write(6,*) 'deltaiso', deltaiso
          goto 98
       else if(overlapvar .eqv. .false.) then
          reptcont = .false.
          goto 83
       end if



       !if(reptcont .eqv. .true.) then






98    if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then
          do l =1,maxlength-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z

             do g = 1,maxlength-1,1
                isobond(m,l,g) = isobond(m,l+1,g)
             end do

          end do

          !write(6,*) 'dx =', dx, 'dy =', dy, 'dz =' , dz

          protcoords(m,maxlength)%x = tempcoord(m,maxlength)%x
          protcoords(m,maxlength)%y = tempcoord(m,maxlength)%y
          protcoords(m,maxlength)%z = tempcoord(m,maxlength)%z

          do g = 2,maxlength-3,1
             isobond(m,maxlength,g) = tempisobond(m,maxlength,g)
          end do

          !write(6,*) time, 'rept 2'      

          successful = successful + 1
          reptbackward = reptbackward + 1
       else
          reptcont = .false.
       end if
       !end if


    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
    end if

  end subroutine reptation

  subroutine endmove
    !performs end chain rotation
    real :: choose2
    logical :: endcont

    m = int(ran2(seed)*(nprotein-1))+1
    !if(m == 0) m = 1
    random = ran2(seed)
    choose2 = ran2(seed)-0.5

    proti = m
    endcont = .true.

    if (choose2 >= 0.0) then
       l = 1
       i = protcoords(m,l+1)%x
       j = protcoords(m,l+1)%y
       k = protcoords(m,l+1)%z
       !write(6,*) 'R'
    else if (choose2 < 0.0) then
       l = maxlength
       i = protcoords(m,l-1)%x
       j = protcoords(m,l-1)%y
       k = protcoords(m,l-1)%z
       !write(6,*) 'L'
    end if


    if(random<= 1.0/6) then
       i = modulo(i,gridsize)+1
    else if (random> 1.0/6 .AND. random<=2.0/6) then
       i = modulo(i - 2,gridsize)+1
    else if (random> 2.0/6 .AND. random<=1.0/2) then
       j = modulo(j,gridsize)+1
    else if (random> 1.0/2 .AND. random<=2.0/3) then                
       j = modulo(j -2,gridsize)+1
    else if (random> 2.0/3 .AND. random<=5.0/6) then
       k  = modulo(k,gridsize)+1
    else if (random> 5.0/6 .AND. random<=1.0) then                
       k = modulo(k-2,gridsize)+1
    end if

    tempcoord(m,l)%x = i
    tempcoord(m,l)%y = j
    tempcoord(m,l)%z = k

    str = l
    fir = l
    st = 4
    fi = maxlength -1
    call overlap
    if(overlapvar .eqv. .true.) then
       call deltae
       goto 98
    else if(overlapvar .eqv. .false.) then
       endcont = .false.
       goto 71
    end if

98  if(deltaiso + deltaint > 0.0 .or. ran2(seed) > (0.5*exp(-1.0*(deltaiso+deltaint)/kT))) then

       protcoords(m,l)%x = i
       protcoords(m,l)%y = j
       protcoords(m,l)%z = k
       protcoords(m,l)%location = section
       do g = 1, maxlength,1
          if(g /= l) isobond(m,l,g) = tempisobond(m,l,g)
       end do
       successful = successful + 1
       ! write(6,*) time, 'end'
    else
       endcont = .false.
    end if


71  if (endcont .eqv. .false.) then
       reject = reject + 1
    end if
  end subroutine endmove


  subroutine comfind
    real :: dcomx,dcomy,dcomz
    t = time +1
    do m = 1,nprotein,1
       comx = 0.0
       comy = 0.0
       comz = 0.0
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
       com(t,m)%x = modulo(protcoords(m,1)%x +(comx/2),gridsize)
       com(t,m)%y = modulo(protcoords(m,1)%y +(comy/2),gridsize)
       com(t,m)%z = modulo(protcoords(m,1)%z +(comz/2),gridsize)

       !write(6,*) 'x =', protcoords(t,m,1)%x, 'y = ', protcoords(t,m,1)%y, 'z = ',protcoords(t,m,1)%z
       !write(6,*) 'dx =', comx, 'dy =', comy, 'dz =', comz

       !write(6,*) 'comx =', com(t,m)%x, 'comy =', com(t,m)%y, 'comz =', com(t,m)%z
    end do
  end subroutine comfind

  subroutine overlap

    overlapvar = .true.
    do pr = 1,nprotein
       do ch1 = str,fir,1
          do ch2 = 1,maxlength,1
             if(ch2 < ch1 .or. ch2 > ch1 .or. pr /= proti) then
                if(protcoords(pr,ch2)%x == tempcoord(proti,ch1)%x.and. &
                     protcoords(pr,ch2)%y == tempcoord(proti,ch1)%y .and. &
                     protcoords(pr,ch2)%z == tempcoord(proti,ch1)%z) then
                   overlapvar = .false.
                   goto 89
                else                   
                   continue
                end if
             end if
          end do
       end do
    end do

  
    !write(6,*) 'no overlap'
89  if(overlapvar .eqv. .false.) then
       continue
       !write(6,*) time,'overlap reject',ch1,ch2,protcoords(pr,ch2)%x,tempcoord(proti,ch1)%x
       !reject = reject + 1
    end if

  end subroutine overlap

  subroutine energy
    !logical :: energypass

    isoenergy = 0.0
    do m= 1,nprotein,1
       do l = 1,maxlength,1
          do f = 1,maxlength,1

             energypassx = .false.
             energypassy = .false.
             energypassz = .false.
             dx = 4.0
             dy = 4.0
             dz = 4.0
             isobond(m,l,f) = 0.0
             if(f < l-2 .or. f > l + 2) then
                !write(6,*) 'l',l,'f',f

                if(abs(protcoords(m,l)%x - protcoords(m,f)%x) == 1.0 .or. abs(protcoords(m,l)%x - &
                     protcoords(m,f)%x) == (gridsize - 1.0)) then
                   dx = 1.0
                   energypassx = .true.

                else if(abs(protcoords(m,l)%x - protcoords(m,f)%x) == 0.0) then
                   dx = 0.0
                   energypassx = .true.
                end if
                !energypass = .false.
                if(abs(protcoords(m,l)%y - protcoords(m,f)%y) == 1.0 .or.abs(protcoords(m,l)%y - &
                     protcoords(m,f)%y) == (gridsize - 1.0)) then
                   dy = 1.0
                   energypassy = .true.
                else if(abs(protcoords(m,l)%y - protcoords(m,f)%y) == 0.0) then
                   dy = 0.0
                   energypassy = .true.
                end if
                !energypass = .false.
                if(abs(protcoords(m,l)%z - protcoords(m,f)%z) == 1.0 .or.abs(protcoords(m,l)%z &
                     - protcoords(m,f)%z) == (gridsize -1.0)) then
                   dz = 1.0
                   energypassz = .true.
                else if(abs(protcoords(m,l)%z - protcoords(m,f)%z) == 0.0) then
                   dz = 0.0
                   energypassz = .true.

                end if
                if(energypassx .eqv. .true. .and. energypassy .eqv. .true. .and. energypassz .eqv. .true.) then
                   sumdebug = dx + dy + dz
                   if(sumdebug == 1.0) then
                      !      write(6,*) f,dx,dy,dz
                      write(6,*) l, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
                      write(6,*) f, protcoords(m,f)%x, protcoords(m,f)%y, protcoords(m,f)%z
                      !      write(6,*) 'deltas', modulo(protcoords(m,l)%x - protcoords(m,f)%x,gridsize), &
                      !           modulo(protcoords(m,l)%y - protcoords(m,f)%y,gridsize), modulo(protcoords(m,l)%z - protcoords(m,f)%z,gridsize)
                      isobond(m,l,f) = 1.0
                   else if(sumdebug /= 1.0) then
                      isobond(m,l,f) = 0.0
                   end if
                else if(energypassx .eqv. .false. .or. energypassy .eqv. .false. .or. energypassz .eqv. .false. ) then
                   isobond(m,l,f) = 0.0
                end if


             end if

 !            write(6,*) m,l,f,isobond(m,l,f)
          end do
       end do
    end do

!!!!!!! interatomic
      do m= 1,nprotein,1
    do g = 1,nprotein
       do l = 1,maxlength,1
          do f = 1,maxlength,1
if(g /=m) then
             energypassx = .false.
             energypassy = .false.
             energypassz = .false.
             dx = 4.0
             dy = 4.0
             dz = 4.0
             intbond(g,m,l,f) = 0.0

                !write(6,*) 'l',l,'f',f

                if(abs(protcoords(m,l)%x - protcoords(g,f)%x) == 1.0 .or. abs(protcoords(m,l)%x - &
                     protcoords(g,f)%x) == (gridsize - 1.0)) then
                   dx = 1.0
                   energypassx = .true.

                else if(abs(protcoords(m,l)%x - protcoords(g,f)%x) == 0.0) then
                   dx = 0.0
                   energypassx = .true.
                end if
                !energypass = .false.
                if(abs(protcoords(m,l)%y - protcoords(g,f)%y) == 1.0 .or.abs(protcoords(m,l)%y - &
                     protcoords(g,f)%y) == (gridsize - 1.0)) then
                   dy = 1.0
                   energypassy = .true.
                else if(abs(protcoords(m,l)%y - protcoords(g,f)%y) == 0.0) then
                   dy = 0.0
                   energypassy = .true.
                end if
                !energypass = .false.
                if(abs(protcoords(m,l)%z - protcoords(g,f)%z) == 1.0 .or.abs(protcoords(m,l)%z &
                     - protcoords(g,f)%z) == (gridsize -1.0)) then
                   dz = 1.0
                   energypassz = .true.
                else if(abs(protcoords(m,l)%z - protcoords(g,f)%z) == 0.0) then
                   dz = 0.0
                   energypassz = .true.

                end if
                if(energypassx .eqv. .true. .and. energypassy .eqv. .true. .and. energypassz .eqv. .true.) then
                   sumdebug = dx + dy + dz
                   if(sumdebug == 1.0) then
                      !      write(6,*) f,dx,dy,dz
                      write(6,*) l, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
                      write(6,*) f, protcoords(g,f)%x, protcoords(g,f)%y, protcoords(g,f)%z
                      !      write(6,*) 'deltas', modulo(protcoords(m,l)%x - protcoords(m,f)%x,gridsize), &
                      !           modulo(protcoords(m,l)%y - protcoords(m,f)%y,gridsize), modulo(protcoords(m,l)%z - protcoords(m,f)%z,gridsize)
                      intbond(m,g,l,f) = 1.0
                   else if(sumdebug /= 1.0) then
                      intbond(m,g,l,f) = 0.0
                   end if
                else if(energypassx .eqv. .false. .or. energypassy .eqv. .false. .or. energypassz .eqv. .false. ) then
                   intbond(m,g,l,f) = 0.0
                end if



!             write(6,*) m,g,l,f,intbond(m,g,l,f)
             end if
          end do
       end do
    end do
    end do
    !   write(93,*) time,isoenergy
  end subroutine energy

  subroutine radiusofgyration
    real :: totrog
    real,dimension(:,:),allocatable :: rog
    allocate(rog(nprotein,maxlength))
    !check this
    do m = 1,nprotein,1
       do l = 1, maxlength,1
          rog(m,l) = (min(modulo(protcoords(m,l)%x-com(time,m)%x,gridsize),(gridsize - &
               modulo(protcoords(m,l)%x-com(time,m)%x,gridsize))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(time,m)%y,gridsize),(gridsize - &
               modulo(protcoords(m,l)%y-com(time,m)%y,gridsize))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(time,m)%z,gridsize),(gridsize - &
               modulo(protcoords(m,l)%z-com(time,m)%z,gridsize))))**2
          totrog = totrog + rog(m,l)
       end do
    end do

    runningaveROG = runningaveROG + totrog


    write(97,*) time,totrog/(nprotein*maxlength),SQRT(totrog/(nprotein*maxlength))

  end subroutine radiusofgyration


  subroutine deltae

    deltaiso = 0
!write(6,*) 'proti',proti
    !do pr = 1,nprotein
       do ch1 = str,fir,1
          do ch2 = st,fi,1

             energypassx = .false.
             energypassy = .false.
             energypassz = .false.
             dx = 4.0
             dy = 4.0
             dz = 4.0

             if(ch2 < ch1-2 .or. ch2 > ch1 + 2) then


                if(abs(tempcoord(proti,ch1)%x - protcoords(pr,ch2)%x) == 1.0 .or. abs(tempcoord(proti,ch1)%x - &
                     protcoords(pr,ch2)%x) == (gridsize - 1.0)) then
                   dx = 1.0
                   energypassx = .true.

                else if(abs(tempcoord(proti,ch1)%x - protcoords(pr,ch2)%x) == 0.0) then
                   dx = 0.0
                   energypassx = .true.
                end if
                !energypass = .false.
                if(abs(tempcoord(proti,ch1)%y - protcoords(pr,ch2)%y) == 1.0 .or.abs(tempcoord(proti,ch1)%y - &
                     protcoords(pr,ch2)%y) == (gridsize - 1.0)) then
                   dy = 1.0
                   energypassy = .true.
                else if(abs(tempcoord(pr,ch1)%y - protcoords(pr,ch2)%y) == 0.0) then
                   dy = 0.0
                   energypassy = .true.
                end if
                !energypass = .false.
                if(abs(tempcoord(proti,ch1)%z - protcoords(pr,ch2)%z) == 1.0 .or.abs(tempcoord(proti,ch1)%z &
                     - protcoords(pr,ch2)%z) == (gridsize -1.0)) then
                   dz = 1.0
                   energypassz = .true.
                else if(abs(tempcoord(proti,ch1)%z - protcoords(pr,ch2)%z) == 0.0) then
                   dz = 0.0
                   energypassz = .true.

                end if
                if(energypassx .eqv. .true. .and. energypassy .eqv. .true. .and. energypassz .eqv. .true.) then
                   sumdebug = dx + dy + dz
                   if(sumdebug == 1.0) then
                      tempisobond(pr,ch1,ch2) = 1.0
                      if(isobond(pr,ch1,ch2) == 0) deltaiso = deltaiso + 1
                   end if

                else if(sumdebug /= 1.0) then
                   tempisobond(pr,ch1,ch2) = 0.0
                end if
             else if(energypassx .eqv. .false. .or. energypassy .eqv. .false. .or. energypassz .eqv. .false. ) then
                tempisobond(pr,ch1,ch2) = 0.0
             end if

             if(tempisobond(pr,ch1,ch2) == 0.0 .and. isobond(pr,ch1,ch2) == 1.0) deltaiso = deltaiso - 1

             !write(6,*) 'isobond',isobond(pr,ch1,ch2)
             !write(6,*) 'tempisobond',tempisobond(pr,ch1,ch2)

          end do
       end do
    !end do

    if(reptcont .eqv. .true.) then
       deltaiso = 0
          do ch1 = str,fir,1
             do ch2 = st,fi,1
                if(tempisobond(proti,ch1,ch2) == 1.0) deltaiso = deltaiso +1
!                write(6,*) 'tempisobond', tempisobond(proti,ch1,ch2)
          end do
       end do

       if(othere == 1) then
          do ch2 = 4,maxlength-1,1
 !            write(6,*) 'isobond one', isobond(proti,othere,ch2)
             if(isobond(proti,othere,ch2) == 1.0) deltaiso = deltaiso - 1
          end do
       else if(othere == maxlength) then
          do ch2 = 2,maxlength-3,1
  !           write(6,*) 'isobond maxlength', isobond(proti,othere,ch2)
             if(isobond(proti,othere,ch2) == 1.0) deltaiso = deltaiso - 1
          end do
       end if
    end if
!    write(6,*) 'deltaiso deltae', deltaiso

    call deltainter
  end subroutine deltae


 subroutine deltainter

    deltaint = 0
!write(6,*) 'proti',proti
    do pr = 1,nprotein
       if(pr /= proti) then
       do ch1 = str,fir,1
          do ch2 = 1,maxlength

             energypassx = .false.
             energypassy = .false.
             energypassz = .false.
             dx = 4.0
             dy = 4.0
             dz = 4.0

           

                if(abs(tempcoord(proti,ch1)%x - protcoords(pr,ch2)%x) == 1.0 .or. abs(tempcoord(proti,ch1)%x - &
                     protcoords(pr,ch2)%x) == (gridsize - 1.0)) then
                   dx = 1.0
                   energypassx = .true.

                else if(abs(tempcoord(proti,ch1)%x - protcoords(pr,ch2)%x) == 0.0) then
                   dx = 0.0
                   energypassx = .true.
                end if
                !energypass = .false.
                if(abs(tempcoord(proti,ch1)%y - protcoords(pr,ch2)%y) == 1.0 .or.abs(tempcoord(proti,ch1)%y - &
                     protcoords(pr,ch2)%y) == (gridsize - 1.0)) then
                   dy = 1.0
                   energypassy = .true.
                else if(abs(tempcoord(pr,ch1)%y - protcoords(pr,ch2)%y) == 0.0) then
                   dy = 0.0
                   energypassy = .true.
                end if
                !energypass = .false.
                if(abs(tempcoord(proti,ch1)%z - protcoords(pr,ch2)%z) == 1.0 .or.abs(tempcoord(proti,ch1)%z &
                     - protcoords(pr,ch2)%z) == (gridsize -1.0)) then
                   dz = 1.0
                   energypassz = .true.
                else if(abs(tempcoord(proti,ch1)%z - protcoords(pr,ch2)%z) == 0.0) then
                   dz = 0.0
                   energypassz = .true.

                end if
                if(energypassx .eqv. .true. .and. energypassy .eqv. .true. .and. energypassz .eqv. .true.) then
                   sumdebug = dx + dy + dz
                   if(sumdebug == 1.0) then
                      tempintbond(proti,pr,ch1,ch2) = 1.0
                      if(intbond(proti,pr,ch1,ch2) == 0) deltaint = deltaint + 1
                 
                else if(sumdebug /= 1.0) then
                   tempintbond(proti,pr,ch1,ch2) = 0.0
                end if
             else if(energypassx .eqv. .false. .or. energypassy .eqv. .false. .or. energypassz .eqv. .false. ) then
                tempintbond(proti,pr,ch1,ch2) = 0.0
             end if

             if(tempintbond(proti,pr,ch1,ch2) == 0.0 .and. intbond(proti,pr,ch1,ch2) == 1.0) deltaint = deltaint - 1
          end do
       end do
       end if
    end do

    if(reptcont .eqv. .true.) then
       deltaint = 0
       do pr = 1,nprotein,1
          if(pr /=proti) then
          do ch1 = str,fir,1
             do ch2 = st,fi,1
                if(tempintbond(proti,pr,ch1,ch2) == 1.0) deltaint = deltaint +1
!                write(6,*) 'tempisobond', tempisobond(proti,ch1,ch2)
          end do
       end do
       end if
end do
if(othere == 1) then
   do pr = 1,nprotein,1
      if(pr /= proti) then
          do ch2 = 4,maxlength-1,1
 !            write(6,*) 'isobond one', isobond(proti,othere,ch2)
             if(intbond(proti,pr,othere,ch2) == 1.0) deltaint = deltaint - 1
          end do
       end if
       end do
       else if(othere == maxlength) then
             do pr = 1,nprotein,1
      if(pr /= proti) then
          do ch2 = 2,maxlength-3,1
  !           write(6,*) 'isobond maxlength', isobond(proti,othere,ch2)
             if(intbond(proti,pr,othere,ch2) == 1.0) deltaint = deltaint - 1
          end do
       end if
    end do
    end if
    end if
 !   write(6,*) 'deltaint delint', deltaint

  end subroutine deltainter


  

  subroutine rms
    real :: msdsum,totmsderror,msd
    real,dimension(:),allocatable :: msdx,msderrorx,msdy,msderrory,msdz,msderrorz,msdt
    allocate(msdx(nprotein))
    allocate(msderrorx(nprotein))
    allocate(msdy(nprotein))
    allocate(msderrory(nprotein))
    allocate(msdz(nprotein))
    allocate(msderrorz(nprotein))

    t = time + 1
    dxtot = 0.0
    dytot = 0.0
    dztot = 0.0
    msd = 0.0
    do m = 1, nprotein,1

       msdx(m) = (min(modulo(com(t,m)%x - com(t-1,m)%x,gridsize), &
            (gridsize - modulo(com(t,m)%x -com(t-1,m)%x,gridsize))))**2
       !write(6,*) 'x =',msdx(m)

       !write(6,*) 'y sortr = ',com(t,m)%y,com(t-1,m)%y
       msdy(m) = (min(modulo(com(t,m)%y -com(t-1,m)%y,gridsize), &
            (gridsize - modulo(com(t,m)%y -com(t-1,m)%y,gridsize))))**2
       !write(6,*) 'y t = ', com(t,m)%y,'y t-1 = ', com(t-1,m)%y
       !write(6,*) 'y =',msdy(m)
       msdz(m) = (min(modulo(com(t,m)%z - com(t-1,m)%z,gridsize), &
            (gridsize - modulo(com(t,m)%z -com(t-1,m)%z,gridsize))))**2

       !msdz(m) = mod(com(t,m)%z - com(t-1,m)%z,gridsize)
       !write(6,*) 'z =',msdz(m)
       totrmsbrute = totrmsbrute + msdx(m) + msdy(m) + msdz(m)
    end do
    msd = sum(msdx) + sum(msdy) + sum(msdz)
    !write(6,*) 'msd =', msd
    !write(6,*) 'msdsum = ', msdsum !this is wrong
    !write(6,*) 'msd brute' , totrmsbrute
    msdsum = msdsum + msd

    if(modulo(time,100) == 0) then
       write(29,*) nprotein*maxlength*(time-equilib),totrmsbrute
       !write(77,*) log(real(nprotein*maxlength*time)), log(totrmsbrute)
    end if
  end subroutine rms



  subroutine dataout

    !write(6,*) 'data'
    write(67,*) nprotein*maxlength
    write(67,*) ' '
    t = time 
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          !maybe put in m
          write(67,*) 'C', 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
               2*protcoords(m,l)%z
          !write(19,*) l, protcoords(t,m,l)%x, protcoords(t,m,l)%y, &
          !    protcoords(t,m,l)%z
          !do g = 1,maxlength,1
          !write(72,*) m,l,g,isobond(m,l,g)
          !end do
       end do
    end do

  end subroutine dataout

  subroutine length
    real :: ddx,ddy,ddz,dxsum,dysum,dzsum
    t = time
    totchainlength = 0.0
    !actualchainlength = 0.0

    do m = 1,nprotein,1
       !chainlength = 0.0
       ddx = 0.0
       ddy = 0
       ddz = 0
       dxsum = 0.0
       dysum = 0.0
       dzsum = 0.0
       do l = 2,maxlength,1

          ddx = protcoords(m,l)%x-protcoords(m,l-1)%x
          if (ddx > 1.5) ddx =-1
          if (ddx < -1.5) ddx = 1


          dxsum = dxsum + ddx
          ddy = protcoords(m,l)%y-protcoords(m,l-1)%y
          if (ddy > 1.5) ddy =-1
          if (ddy < -1.5) ddy = 1

          dysum = dysum + ddy
          ddz = protcoords(m,l)%z-protcoords(m,l-1)%z
          if (ddz > 1.5) ddz =-1
          if (ddz < -1.5) ddz = 1

          dzsum = dzsum + ddz

          !write(6,*) 'ddx =', ddx, 'ddy = ', ddy, 'ddz = ', ddz
          !chainlength(m) = SQRT((min, abs(gridsize- &
          ! modulo(protcoords(t,m,1)%x - protcoords(t,m,maxlength)%x,gridsize))))**2 + &
          !(min(modulo(protcoords(t,m,1)%y - protcoords(t,m,maxlength)%y,gridsize), abs(gridsize - &
          !modulo(protcoords(t,m,1)%y - protcoords(t,m,maxlength)%y,gridsize))))**2 + &
          !(min(modulo(protcoords(t,m,1)%z - protcoords(t,m,maxlength)%z,gridsize), abs(gridsize - &
          !modulo(protcoords(t,m,1)%z - protcoords(t,m,maxlength)%z,gridsize))))**2)
       end do
       chainlength = ((dxsum**2) + (dysum**2) + (dzsum**2))
       totchainlength = totchainlength + chainlength
    end do
    !totchainlength = sum(chainlength)
    runningaveEtE = runningaveEtE + totchainlength  !sort this out
    avechainlength = totchainlength/nprotein

    actualchainlength = actualchainlength + avechainlength
    write(91,*) t,avechainlength !chainlength(1)
  end subroutine length

  subroutine error
    lengthvar = 0.0
    do t = 1,maxtime,1
       !lengthvar =
    end do

  end subroutine error


  subroutine bonding
    integer :: bx,by,bz,f
    bx = 0
    by = 0
    bz = 0
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          do f = 1,nprotein,1
             if (f /=m) then
                do g = 1,maxlength
                   !if (protcoords(t,m,l)%location == protcoords(t,f,g)%location) then
                   if (modulo(protcoords(m,l)%x - protcoords(f,g)%x,gridsize) < separation &
                        .and. modulo(protcoords(m,l)%y - protcoords(f,g)%y,gridsize) < separation &
                        .and. modulo(protcoords(m,l)%z - protcoords(f,g)%z,gridsize) < separation) then
                      bz = bz+1
                   end if
                   !end if
                end do
             else
                continue
             end if
          end do
       end do
    end do
    !write(51,*) time,bz

    isocount = 0
    do m = 1,nprotein,1
       do l = 1,maxlength,1
          do f = 1,nprotein,1
             do g = 1,maxlength,1
                if(m /= f .or. l /= g) then
                   sumdebug = abs(protcoords(m,l)%x - protcoords(f,g)%x) + &
                        abs(protcoords(m,l)%y - protcoords(f,g)%y) + abs(protcoords(m,l)%z - protcoords(f,g)%z)
                   if(sumdebug == 1.0 .or. sumdebug == gridsize-1.0) then
                      isobond(m,l,g) = 1.0
                      isocount = isocount + 1
                   end if
                end if
             end do
          end do
       end do
    end do
    isoenergy = isocount/2
    write(6,*) 'Energy =', isoenergy


  end subroutine bonding


  subroutine debug


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

                      !write(6,*) 'FAILLLLLLLLLLLLLLLLLL'
                   end if
                end if
                if(l>1) then
                   sumdebug = abs(protcoords(m,l)%x - protcoords(m,l-1)%x) + &
                        abs(protcoords(m,l)%y - protcoords(m,l-1)%y) + abs(protcoords(m,l)%z - protcoords(m,l-1)%z)
                   if(sumdebug /= 1.0 .and. sumdebug /= gridsize-1.0) then
                      !write(6,*) 'SPLIT'
                      finalfail = .true.
                   end if
                end if
             end do
          end do
       end do
    end do


  end subroutine debug

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
