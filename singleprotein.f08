program move
  implicit none
  
  
  integer :: maxno,protno,seed,natoms,maxsize,maxtime,reptbackward,reptforward,equilib
  real :: i,j,k,gridsize,randomz,deltax1,deltax2,dx,dy,dz,msd
  integer :: l,m,hold,z1,z2,z3,z4,z5,z6,control
  real :: deltay1,deltay2,deltaz1,deltaz2,totrmsbrute
  integer :: nprotein,split,section,time,count
  integer :: maxlength,g,xl,yl,zl,t,r,N
  real, external :: ran2
  logical :: rightangle,endchain,run,pivotx,pivoty,pivotz,run2,rac,exist,fail
  real :: totdisp,totchainlength,avechainlength,probs,comx,comy,comz,totrmserror
  real :: actualchainlength,totvar,deviate,totdeviate,varian,separation,lengthvar,random
  real,dimension(:),allocatable :: variance,disp,rmserror,dxtot,dytot,dztot
  integer :: gate,xbk,ybk,zbk,protpass,f,scan,c
  Real:: t1,t2,t3,t4,t5,t6,di,dj,dk
  integer :: cont1,cont2,cont3,cont4,totcont,successful,reject
  character(len = 10) :: commandread,commandread2,comread3,comread4,comread5,comread6,comread7,comread8,comread9
  real :: runningaveROG,runningaveEtE,chainlength
  !real :: testlength,testtotlength
  integer :: piv,right,end,crank,rept
  
  
  type protein
     real :: x,y,z
     integer :: location
  end type protein
  type(protein),dimension(:,:),allocatable :: protcoords
  type(protein),dimension(:,:),allocatable :: com
  
  open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')
  open(67, file = 'move.xyz', action = 'write')
  open(29, file = 'rms.dat', action = 'write')
  !open(31, file = 'chainlength.dat', action = 'write')
  !open(51, file = 'bonds.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  !open(77, file = 'logmsd.dat', action = 'write')
  open(79, file = 'runningave.dat', action = 'write')
  open(99, file = 'deltas.dat', action = 'write')
  
  totrmsbrute = 0.0
  runningaveROG = 0.0
  runningaveEtE = 0.0
  
  
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
  read(comread8, '(i1)') maxtime
  call get_command_argument(9,comread9)
  read(comread9, '(i1)') equilib 
  
  !reads in seed from commandline
  
  N = maxlength
  allocate(protcoords(nprotein,maxlength))
  allocate(com(maxtime+1,nprotein))
  !allocate(chainlength(nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(dxtot(nprotein))
  allocate(dytot(nprotein))
  allocate(dztot(nprotein))
  count = 0
  time = 0
  totdisp = 0.0
  
  call foundation
  
  !call dataout
  call comfind
  actualchainlength = 0.0
  do time = 1,maxtime
     fail = .false.
     !call comfind
     !write(6,*) 'a'
     !if (mod(time,10) == 0) then
     call dataout
     !end if
     !write(6,*) 'b'
     !call bonding
     !write(6,*) 'c'
     if(time >equilib) then
     call length
     end if
     !write(6,*) 'd'
     
     !write(6,*) 'e'
     call positioning
     call comfind
     
     if (time > equilib) then
        call rms
        call radiusofgyration
        end if
     call debug
     
     write(79,*) time-equilib, runningaveEtE/(time-equilib), runningaveROG/(time-equilib)
     
     !write(6,*) 'g'
     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else
        write(6,*) 'step =',time
     end if
     
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
    do m = 1,nprotein
       do l = 1,maxlength
          read(23,*) BIN, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z
          !write(67,*) 'C', 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
          !2*protcoords(m,l)%z
          !write(6,*) protcoords(m,l)%x
       end do
    end do
    
  end subroutine foundation
  
  
  
 
  subroutine positioning
    !randomly selects one of the 4 moves
    
    integer :: contchoose,a1,a2,a3,a4,nmoves
    
    
    t = time + 1
    do scan = 1, nprotein*maxlength
       count = count + 1
       rac = .true.
       run = .true.
       run2 = .true. 
       
       
       
       randomz = ran2(seed)
       nmoves = piv + crank + end + right + rept
       
       if (randomz <= (1.0*end)/ nmoves .and. end == 1) then
          
          call endmove
          
       else if (randomz > (1.0*end)/nmoves .and. randomz <= (1.0*end+crank)/nmoves .and. crank == 1) then
          call crankshaftmove
          
       else if (randomz > (1.0*end + crank)/nmoves .and. randomz <= (1.0*end + crank + right)/nmoves .and. right ==1) then
          call rightanglemove
          
       else if (randomz > (1.0*end + crank + right)/nmoves .and. randomz <= (1.0*end + crank + right + rept)/nmoves &
            .and. rept ==1) then
          
          call reptation
       else if (randomz > (1.0*end + crank + right+rept)/nmoves .and. randomz <= (1.0*end + crank + right + rept + piv)/nmoves &
            .and. piv == 1) then 
          
          call pivot
       end if
    end do
  end subroutine positioning
  
  
  subroutine crankshaftmove
    !performs crankshaft move
    logical :: crankcont
    real :: dx11,dx12,dy11,dy12,dz11,dz12
    !this means that m/=0 and is up to nprotein
    m =int(ran2(seed)*(nprotein-1))+1
    !this means that l cannot be maxlength,maxlength - 1 or l
    l = int(ran2(seed)*(maxlength-4))+2

    !write(6,*) 'm =',m,'l =', l
    probs = ran2(seed)
    pivotx = .false.
    pivoty = .false.
    pivotz = .false.
    
    
    
    !check to see if 4 atoms are in the same plane
    if (protcoords(m,l-1)%x == protcoords(m,l)%x &
         .and.  protcoords(m,l-1)%x == protcoords(m,l+1)%x &
         .and.protcoords(m,l-1)%x == protcoords(m,l+2)%x) then
       pivotx = .true.                  
    else if (protcoords(m,l-1)%y == protcoords(m,l)%y .and. &
         protcoords(m,l-1)%y == protcoords(m,l+1)%y &
         .and.protcoords(m,l-1)%y == protcoords(m,l+2)%y ) then
       pivoty = .true.                               
    else if (protcoords(m,l-1)%z == protcoords(m,l)%z &
         .and. protcoords(m,l-1)%z == protcoords(m,l+1)%z &
         .and.protcoords(m,l-1)%z == protcoords(m,l+2)%z ) then
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
       
       dy11 = min(protcoords(m,l)%y - protcoords(m,l-1)%y, gridsize - &
            protcoords(m,l)%y - protcoords(m,l-1)%y)
       dy12 = min(protcoords(m,l+1)%y - protcoords(m,l+2)%y, gridsize - &
            protcoords(m,l+1)%y - protcoords(m,l+2)%y)
       
       dz11 = min(protcoords(m,l)%z - protcoords(m,l-1)%z, gridsize - &
            protcoords(m,l)%z - protcoords(m,l-1)%z)
       dz12 = min(protcoords(m,l+1)%z - protcoords(m,l+2)%z, gridsize - &
            protcoords(m,l+1)%z- protcoords(m,l+2)%z)

       dy = min(protcoords(m,l)%y - protcoords(m,l+1)%y, gridsize - &
            protcoords(m,l)%y - protcoords(m,l+1)%y)
       dz = min(protcoords(m,l)%z - protcoords(m,l+1)%z, gridsize - &
            protcoords(m,l)%z - protcoords(m,l+1)%z)

       
       if ( dy11 == dy12 .and. abs(dy12) ==1.0 .and. dz == 0.0 .and. &
            modulo(protcoords(m,l)%y - protcoords(m,l+1)%y,gridsize) ==0.0) then
          dy = protcoords(m,l)%y - protcoords(m,l-1)%y
          if (dy > 1.5) dy = -1.0
          if (dy < 1.5) dy = 1.0
          dz = 0.0 
       else  if (dz11 == dz12 .and. abs(dz12) ==1.0 .and. dy == 0.0 .and. &
            modulo(protcoords(m,l)%z - protcoords(m,l+1)%z,gridsize) ==0.0) then
          dy = 0.0
          dz = protcoords(m,l)%z - protcoords(m,l-1)%z
          if (dz > 1.5) dz = -1.0
          if (dz < 1.5) dz = 1.0
       else
          crankcont = .false.
          goto 43
       end if
       
       if (probs <= (1.0/3)) then
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == modulo(protcoords(m,l)%y-(2*dy)-1,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l)%x .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l)%z-(2*dz)-1,gridsize)+1 .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l+1)%y-(2*dy)-1,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l+1)%x .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l+1)%z-(2*dz)-1,gridsize)+1) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do

          
         protcoords(m,l)%y = modulo(protcoords(m,l)%y-(2*dy)-1,gridsize)+1
          protcoords(m,l)%z = modulo(protcoords(m,l)%z-(2*dz)-1,gridsize)+1
          protcoords(m,l+1)%y = modulo(protcoords(m,l+1)%y-(2*dy)-1,gridsize)+1
          protcoords(m,l+1)%z = modulo(protcoords(m,l+1)%z-(2*dz)-1,gridsize)+1

          
          !xl = int((protcoords(t,m,l)%x*10)/split)-1
          !yl = 1+int(((protcoords(t,m,l)%y*10)/split) -1)*N
          !zl = int(((protcoords(t,m,l)%z*10)/split) -1)*(N**2)
          !section = xl + yl + zl
          !protcoords(t,m,l)%location = section
          
          successful = successful + 1
          
       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == protcoords(m,l-1)%y .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l)%x-2,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l-1)%z .and. &
                     protcoords(g,f)%y == protcoords(m,l+2)%y .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l+1)%x-2,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+2)%z) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do

          
          protcoords(m,l)%y = protcoords(m,l-1)%y
          protcoords(m,l)%z = protcoords(m,l-1)%z
          protcoords(m,l+1)%y = protcoords(m,l+2)%y
          protcoords(m,l+1)%z = protcoords(m,l+2)%z
          protcoords(m,l)%x = modulo(protcoords(m,l)%x-2,gridsize)+1
          protcoords(m,l+1)%x = modulo(protcoords(m,l+1)%x-2,gridsize)+1
          
          
          successful = successful +1 
       else if (probs > (2.0/3) .and. probs <= (1.0)) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == protcoords(m,l-1)%y .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l)%x,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l-1)%z .and. &
                     protcoords(g,f)%y == protcoords(m,l+2)%y .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l+1)%x,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+2)%z) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do

          
          
          protcoords(m,l)%y = protcoords(m,l-1)%y
          protcoords(m,l)%z = protcoords(m,l-1)%z
          protcoords(m,l+1)%y = protcoords(m,l+2)%y
          protcoords(m,l+1)%z = protcoords(m,l+2)%z
          protcoords(m,l)%x = modulo(protcoords(m,l)%x,gridsize)+1
          protcoords(m,l+1)%x = modulo(protcoords(m,l+1)%x,gridsize)+1
          
          successful = successful + 1
       
       end if
    end if
    
    
    !!pivot y
    
    if (pivoty .eqv. .true.) then
       
       dx11 = min(protcoords(m,l)%x - protcoords(m,l-1)%x, gridsize - &
            protcoords(m,l)%x - protcoords(m,l-1)%x)
       dx12 = min(protcoords(m,l+1)%x - protcoords(m,l+2)%x, gridsize - &
            protcoords(m,l+1)%x - protcoords(m,l+2)%x)
       dz11 = min(protcoords(m,l)%z - protcoords(m,l-1)%z, gridsize - &
            protcoords(m,l)%z - protcoords(m,l-1)%z)
       dz12 = min(protcoords(m,l+1)%z - protcoords(m,l+2)%z, gridsize - &
            protcoords(m,l+1)%z - protcoords(m,l+2)%z)
       dx =  min(protcoords(m,l)%x - protcoords(m,l+1)%x, gridsize - &
            protcoords(m,l)%x - protcoords(m,l+1)%x)
       dz = min(protcoords(m,l)%z - protcoords(m,l+1)%z, gridsize - &
            protcoords(m,l)%z - protcoords(m,l+1)%z)
       if (dx11 == dx12 .and. abs(dx12) ==1.0 .and. dz == 0.0 .and.  &
            modulo(protcoords(m,l)%x - protcoords(m,l+1)%x,gridsize) ==0.0) then
          dx = protcoords(m,l)%x - protcoords(m,l-1)%x
          if (dx > 1.5) dx = -1.0
          if (dx < 1.5) dx = 1.0
          dz = 0.0 
       else  if (dz11 == dz12 .and. abs(dz12) ==1.0 .and. dx == 0.0 .and. &
            modulo(protcoords(m,l)%z - protcoords(m,l+1)%z,gridsize) ==0.0) then
          dx = 0.0
          dz = protcoords(m,l)%z - protcoords(m,l-1)%z
          if (dz > 1.5) dz = -1.0
          if (dz < 1.5) dz = 1.0
       else
          crankcont = .false.
          goto 43
          
       end if
       
       if (probs <= (1.0/3)) then
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%x == modulo(protcoords(m,l)%x-2*dx-1,gridsize)+1 .and. &
                     protcoords(g,f)%y == protcoords(m,l)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l)%z-2*dz-1,gridsize)+1 .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l+1)%x-2*dx-1,gridsize)+1 .and. &
                     protcoords(g,f)%y == protcoords(m,l+1)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l+1)%z-2*dz-1,gridsize)+1) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%x = modulo(protcoords(m,l)%x-2*dx-1,gridsize)+1
          protcoords(m,l)%z = modulo(protcoords(m,l)%z-2*dz-1,gridsize)+1
          protcoords(m,l+1)%x = modulo(protcoords(m,l+1)%x-2*dx-1,gridsize)+1
          protcoords(m,l+1)%z = modulo(protcoords(m,l+1)%z-2*dz-1,gridsize)+1
          
          successful = successful +1
       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%x == protcoords(m,l-1)%x .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l)%y-2,gridsize)+1  .and. &
                     protcoords(g,f)%z == protcoords(m,l-1)%z .and. &
                     protcoords(g,f)%x == protcoords(m,l+2)%x .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l+1)%y-2,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+2)%z) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%x = protcoords(m,l-1)%x
          protcoords(m,l)%z = protcoords(m,l-1)%z
          protcoords(m,l+1)%x = protcoords(m,l+2)%x
          protcoords(m,l+1)%z = protcoords(m,l+2)%z
          protcoords(m,l)%y = modulo(protcoords(m,l)%y-2,gridsize)+1
          protcoords(m,l+1)%y = modulo(protcoords(m,l+1)%y-2,gridsize)+1
          
          
          successful = successful + 1
       else if (probs > (2.0/3) .and. probs <= 1.0) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%x == protcoords(m,l-1)%x .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l)%y,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+2)%z .and. &
                     protcoords(g,f)%x == protcoords(m,l+2)%x .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l+1)%y,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+2)%z) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%x = modulo(protcoords(m,l)%x-dx-1,gridsize)+1
          protcoords(m,l)%z = modulo(protcoords(m,l)%z-dz-1,gridsize)+1
          protcoords(m,l+1)%x = modulo(protcoords(m,l+1)%x-dx-1,gridsize)+1
          protcoords(m,l+1)%z = modulo(protcoords(m,l+1)%z-dz-1,gridsize)+1
          protcoords(m,l)%y = modulo(protcoords(m,l)%y,gridsize)+1
          protcoords(m,l+1)%y = modulo(protcoords(m,l+1)%y,gridsize)+1
          
          
          successful = successful +1
       end if
    end if
    
    !pivot z
    
    
    if (pivotz .eqv. .true.) then

       dy11 = min(protcoords(m,l)%y - protcoords(m,l-1)%y, gridsize - &
            protcoords(m,l)%y - protcoords(m,l-1)%y)
       dy12 = min(protcoords(m,l+1)%y - protcoords(m,l+2)%y, gridsize - &
            protcoords(m,l+1)%y - protcoords(m,l+2)%y)
       
       dx11 = min(protcoords(m,l)%x - protcoords(m,l-1)%x, gridsize - &
            protcoords(m,l)%x - protcoords(m,l-1)%x)
       dx12 = min(protcoords(m,l+1)%x - protcoords(m,l+2)%x, gridsize - &
            protcoords(m,l+1)%x - protcoords(m,l+2)%x)
       dx = min(protcoords(m,l)%x - protcoords(m,l+1)%x, gridsize - &
            protcoords(m,l)%x - protcoords(m,l+1)%x)
       dy = min(protcoords(m,l)%y - protcoords(m,l+1)%y, gridsize - &
            protcoords(m,l)%y - protcoords(m,l+1)%y)
       
       
       if (dy11 == dy12 .and. abs(dy12) ==1.0 .and. dx == 0.0 .and.  &
            modulo(protcoords(m,l)%y - protcoords(m,l+1)%y,gridsize) ==0.0) then
          dy = protcoords(m,l)%y - protcoords(m,l-1)%y
          if (dy > 1.5) dy = -1.0
          if (dy < 1.5) dy = 1.0
          dx = 0.0 
       else  if (dx11 ==dx12 .and. abs(dx12) ==1.0 .and. dy == 0.0 .and. &
            mod(protcoords(m,l)%x - protcoords(m,l+1)%x,gridsize) ==0.0) then
          dy = 0.0
          dx = protcoords(m,l)%x - protcoords(m,l-1)%x
          if (dx > 1.5) dx = -1.0
          if (dx < 1.5) dx = 1.0
       else
          crankcont = .false.
          goto 43
       end if
       
       if (probs <= (1.0/3)) then
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == modulo(protcoords(m,l)%y-2*dy-1,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l)%z .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l)%x-2*dx-1,gridsize)+1 .and. &
                     protcoords(g,f)%y == modulo(protcoords(m,l+1)%y-2*dy-1,gridsize)+1 .and. &
                     protcoords(g,f)%z == protcoords(m,l+1)%z .and. &
                     protcoords(g,f)%x == modulo(protcoords(m,l+1)%x-2*dx-1,gridsize)+1) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%y = modulo(protcoords(m,l)%y-2*dy-1,gridsize)+1
          protcoords(m,l)%x = modulo(protcoords(m,l)%x-2*dx-1,gridsize)+1
          protcoords(m,l+1)%y = modulo(protcoords(m,l+1)%y-2*dy-1,gridsize)+1
          protcoords(m,l+1)%x = modulo(protcoords(m,l+1)%x-2*dx-1,gridsize)+1
          
          successful = successful + 1
       else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == protcoords(m,l-1)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l)%z-2,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l-1)%x .and. &
                     protcoords(g,f)%y == protcoords(m,l+2)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l+1)%z-2,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l+2)%x) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%y = protcoords(m,l-1)%y
          protcoords(m,l)%x = protcoords(m,l-1)%x
          protcoords(m,l+1)%y = protcoords(m,l+2)%y
          protcoords(m,l+1)%x = protcoords(m,l+2)%x
          protcoords(m,l)%z = modulo(protcoords(m,l)%z-2,gridsize)+1
          protcoords(m,l+1)%z = modulo(protcoords(m,l+1)%z-2,gridsize)+1
          
          successful = successful + 1
       else if (probs > (2.0/3) .and. probs <= 1.0) then
          
          do g = 1,nprotein
             do f = 1,maxlength
                if(protcoords(g,f)%y == protcoords(m,l-1)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l)%z,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l-1)%x .and. &
                     protcoords(g,f)%y == protcoords(m,l+2)%y .and. &
                     protcoords(g,f)%z == modulo(protcoords(m,l+1)%z,gridsize)+1 .and. &
                     protcoords(g,f)%x == protcoords(m,l+2)%x) then
                   crankcont = .false.
                   goto 43
                else
                   continue
                end if
             end do
          end do
          
          protcoords(m,l)%y = protcoords(m,l-1)%y
          protcoords(m,l)%x = protcoords(m,l-1)%x
          protcoords(m,l+1)%y = protcoords(m,l+2)%y
          protcoords(m,l+1)%x = protcoords(m,l+2)%x
          protcoords(m,l)%z = modulo(protcoords(m,l)%z,gridsize)+1
          protcoords(m,l+1)%z = modulo(protcoords(m,l+1)%z,gridsize)+1
          
          successful = successful + 1
       end if
    end if
    
    
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
          
          i = modulo(protcoords(m,l)%x+di-1,gridsize)+1
          j = modulo(protcoords(m,l)%y+dj-1,gridsize)+1
          k = modulo(protcoords(m,l)%z+dk-1,gridsize)+1
          
          
          do g = 1, nprotein
             do f = 1,maxlength
                if (i == protcoords(g,f)%x .and. j == protcoords(g,f)%y &
                     .and. k == protcoords(g,f)%z) then
                   
                   rac = .false.
                   goto 37              
                else
                   continue                
                end if
             end do
          end do
           
          protcoords(m,l)%x = i
          protcoords(m,l)%y = j
          protcoords(m,l)%z = k
          
           
 

          successful = successful + 1
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
    real,dimension(:,:),allocatable :: delx,dely,delz
    allocate(delx(nprotein,maxlength))
    allocate(dely(nprotein,maxlength))   
    allocate(delz(nprotein,maxlength))
    
    !write(6,*) 'here comes the pivot!!'
    
    choose3 = ran2(seed)
    pivcont = .true. 
    m =int(ran2(seed)*(nprotein-1))+1
    
    l = int(ran2(seed)*(maxlength-3))+2
    !write(6,*) 'm = ', m
    !write(6,*) 'l =', l
    t = time + 1
    
    
    
    if(l > maxlength/2) then
       
       
       do g = l+1,maxlength
          
          delx(m,g) = protcoords(m,g)%x - protcoords(m,l)%x
          xhold = delx(m,g)
          if(delx(m,g) > gridsize/2) delx(m,g) =  xhold-gridsize
          if(delx(m,g) < -1.0*gridsize/2) delx(m,g) = gridsize + xhold
          dely(m,g) = protcoords(m,g)%y - protcoords(m,l)%y
          yhold = dely(m,g)
          if(dely(m,g) > gridsize/2) dely(m,g) = yhold - gridsize
          if(dely(m,g) < -1.0*gridsize/2) dely(m,g) = gridsize + yhold
          delz(m,g) = protcoords(m,g)%z - protcoords(m,l)%z
          zhold = delz(m,g)
          if(delz(m,g) > gridsize/2) delz(m,g) =  zhold-gridsize
          if(delz(m,g) < -1.0*gridsize/2) delz(m,g) = gridsize + zhold
          write(99,*) delx(m,g),dely(m,g),delz(m,g) !,xhold,yhold,zhold
       end do
       if (choose3 <= 1.0/5) then
          do b = l+1,maxlength 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - dely(m,b)-1,gridsize) +1  == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y + delx(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          
          do b = l+1,maxlength 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x - dely(m,b)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y + delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z 
             successful = successful + 1
             
          end do
          !write(6,*) '1a'       
       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then 
          
          do b = l+1,maxlength 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x + dely(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l+1,maxlength 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x + dely(m,b)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y - delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          !write(6,*) '2a'
          
       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then
          
          do b = l+1,maxlength 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y - dely(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l+1,maxlength                 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x - delx(m,b)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y - dely(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          ! write(6,*) '3a'
          
       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then
          
          do b = l+1,maxlength 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - delz(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%z + delx(m,b)-1,gridsize)+1 == protcoords(f,g)%z &
                        .and. protcoords(m,b)%y == protcoords(f,g)%y) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l+1,maxlength 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x - delz(m,b)-1,gridsize)+1
             protcoords(m,b)%z = modulo(protcoords(m,l)%z + delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          !write(6,*) '4a'
          
       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then
          
          do b = l+1,maxlength 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x + delz(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%z - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%z &
                        .and. protcoords(m,b)%y == protcoords(f,g)%y) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b = l+1,maxlength                
             protcoords(m,b)%x = modulo(protcoords(m,l)%x + delz(m,b)-1,gridsize)+1
             protcoords(m,b)%z = modulo(protcoords(m,l)%z - delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          !write(6,*) '5a'
       end if
       
    else if (l <= maxlength/2) then
       !continue
       !else if (l == 100000) then
       do g = l-1,1,-1
          
          delx(m,g) = protcoords(m,g)%x - protcoords(m,l)%x
          xhold = delx(m,g)
          if(delx(m,g) > gridsize/2) delx(m,g) =  xhold-gridsize
          if(delx(m,g) < -1.0*gridsize/2) delx(m,g) = gridsize + xhold
          dely(m,g) = protcoords(m,g)%y - protcoords(m,l)%y
          yhold = dely(m,g)
          if(dely(m,g) > gridsize/2) dely(m,g) = yhold - gridsize
          if(dely(m,g) < -1.0*gridsize/2) dely(m,g) = gridsize + yhold
          delz(m,g) = protcoords(m,g)%z - protcoords(m,l)%z
          zhold = delz(m,g)
          if(delz(m,g) > gridsize/2) delz(m,g) =  zhold-gridsize
          if(delz(m,g) < -1.0*gridsize/2) delz(m,g) = gridsize + zhold
          
          write(99,*) delx(m,g),dely(m,g),delz(m,g),xhold,yhold,zhold
          
       end do
       if (choose3 <= 1.0/5) then
          do b = l-1,1,-1 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - dely(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y + delx(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          
          do b = l-1,1,-1 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x -dely(m,b)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y + delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z 
             successful = successful + 1
             
          end do
          write(99,*) '1a'
       else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then 
          
          do b = l-1,1,-1 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x + dely(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l-1,1,-1 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x + dely(m,b)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y - delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          write(99,*) '2a'
          
          
       else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then
          
          do b = l-1,1,-1 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%y - dely(m,b)-1,gridsize)+1 == protcoords(f,g)%y &
                        .and. protcoords(m,b)%z == protcoords(f,g)%z) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l-1,1,-1                 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x -delx(m,l)-1,gridsize)+1
             protcoords(m,b)%y = modulo(protcoords(m,l)%y - dely(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          write(99,*) '3a'
          
       else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then
          
          do b = l-1,1,-1 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x - delz(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%z + delx(m,b)-1,gridsize)+1 == protcoords(f,g)%z &
                        .and. protcoords(m,b)%y == protcoords(f,g)%y) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b =l-1,1,-1 
             protcoords(m,b)%x = modulo(protcoords(m,l)%x - delz(m,b)-1,gridsize)+1
             protcoords(m,b)%z = modulo(protcoords(m,l)%z + delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          write(99,*) '4a'
       else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then
          
          do b = l-1,1,-1 
             do f = 1,nprotein
                do g = 1,maxlength
                   if(modulo(protcoords(m,l)%x + delz(m,b)-1,gridsize)+1 == protcoords(f,g)%x .and. &
                        modulo(protcoords(m,l)%z - delx(m,b)-1,gridsize)+1 == protcoords(f,g)%z &
                        .and. protcoords(m,b)%y == protcoords(f,g)%y) then
                      pivcont = .false.
                      goto 75
                   end if
                end do
             end do
          end do
          do b = l-1,1,-1                
             protcoords(m,b)%x = modulo(protcoords(m,l)%x + delz(m,b)-1,gridsize)+1
             protcoords(m,b)%z = modulo(protcoords(m,l)%z - delx(m,b)-1,gridsize)+1
             !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
             successful = successful + 1
             
          end do
          write(99,*) '5a'
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
    logical :: reptcont
    
    t = time + 1
    choose = ran2(seed) - 0.5
    reptcont = .true.
    m =int(ran2(seed)*(nprotein-1))+1
    c = m
    
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
       dxx = protcoords(c,2)%x - protcoords(c,1)%x
       if (dxx > 1.5) dxx = -1.0
       if(dxx <-1.5) dxx = 1.0
       dyx = protcoords(c,2)%y - protcoords(c,1)%y
       if (dyx > 1.5) dyx = -1.0
       if (dyx < -1.5) dyx = 1.0
       dzx = protcoords(c,2)%z - protcoords(c,1)%z
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
       
       
       do g = 1,nprotein
          do f = 1,maxlength
             if(protcoords(g,f)%x == modulo(protcoords(c,1)%x + dx -1,gridsize)+1 .and. &
                  protcoords(g,f)%y == modulo(protcoords(c,1)%y + dy -1,gridsize)+1 .and. &
                  protcoords(g,f)%z == modulo(protcoords(c,1)%z + dz -1,gridsize)+1) then
                
                reptcont = .false.
                goto 83
             else
                continue
             end if
          end do
       end do

       if(reptcont .eqv. .true.) then
       do l =maxlength,2,-1
          protcoords(c,l)%x = protcoords(c,l-1)%x
          protcoords(c,l)%y = protcoords(c,l-1)%y
          protcoords(c,l)%z = protcoords(c,l-1)%z
          
          !successful = successful+1
       end do

       !write(6,*) 'dx =', dx, 'dy =', dy, 'dz = ', dz
       protcoords(c,1)%x = modulo(protcoords(c,1)%x + dx-1,gridsize)+1
       protcoords(c,1)%y = modulo(protcoords(c,1)%y + dy-1,gridsize)+1
       protcoords(c,1)%z = modulo(protcoords(c,1)%z + dz-1,gridsize)+1
       
       
       successful = successful + 1
       reptforward = reptforward + 1
       end if
    else if (choose < 0.0) then
       
       
       dxx = protcoords(c,maxlength-1)%x - protcoords(c,maxlength)%x
       if (dxx > 1.5) dxx = -1.0
       if(dxx <-1.5) dxx = 1.0
       dyx = protcoords(c,maxlength-1)%y - protcoords(c,maxlength)%y
       if (dyx > 1.5) dyx = -1.0
       if (dyx < -1.5) dyx = 1.0
       dzx = protcoords(c,maxlength-1)%z - protcoords(c,maxlength)%z
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
          
       else if(direc <= (1.0*g1+g2+g3+g4)/5 .and. direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
          dy = -1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5)/5 .and. &
            direc > (1.0*g1+g2+g3+g4)/5 .and.  g5 == 1) then
          dz = 1.0
       else if(direc <= (1.0*g1+g2+g3+g4+g5+g6)/5 .and. &
            direc > (1.0*g1+g2+g3+g4+g5)/5 .and.  g6 == 1) then
          dz = -1.0
       end if
                
       
       
       do g = 1,nprotein
          do f = 1,maxlength
             if(protcoords(g,f)%x == modulo(protcoords(c,maxlength)%x + dx -1,gridsize)+1 .and. &
                  protcoords(g,f)%y == modulo(protcoords(c,maxlength)%y + dy -1,gridsize)+1 .and. &
                  protcoords(g,f)%z == modulo(protcoords(c,maxlength)%z + dz -1,gridsize)+1) then
                reptcont = .false.
                
                goto 83
             else
                continue
             end if
          end do
       end do

             if(reptcont .eqv. .true.) then
       
       do l =1,maxlength-1
          protcoords(c,l)%x = protcoords(c,l+1)%x
          protcoords(c,l)%y = protcoords(c,l+1)%y
          protcoords(c,l)%z = protcoords(c,l+1)%z
          
       end do

       !write(6,*) 'dx =', dx, 'dy =', dy, 'dz =' , dz
       
       protcoords(c,maxlength)%x = modulo(protcoords(c,maxlength)%x + dx-1,gridsize)+1
       protcoords(c,maxlength)%y = modulo(protcoords(c,maxlength)%y + dy-1,gridsize)+1
       protcoords(c,maxlength)%z = modulo(protcoords(c,maxlength)%z + dz-1,gridsize)+1
       
       
       
       successful = successful + 1
       reptbackward = reptbackward + 1
       end if
       
       
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
    
    
    endcont = .true.
    
    if (choose2 >= 0.0) then
       l = 1
       i = protcoords(m,l+1)%x
       j = protcoords(m,l+1)%y
       k = protcoords(m,l+1)%z
    else if (choose2 <= 0.0) then
       l = maxlength
       i = protcoords(m,l-1)%x
       j = protcoords(m,l-1)%y
       k = protcoords(m,l-1)%z
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
    do g = 1, nprotein
       do f = 1,maxlength
          if (i == protcoords(g,f)%x .and. j == protcoords(g,f)%y &
               .and. k == protcoords(g,f)%z ) then
             endcont = .false.
             goto 71 
          else
             continue
          end if
       end do
    end do
           
    
              
    protcoords(m,l)%x = i
    protcoords(m,l)%y = j
    protcoords(m,l)%z = k
    protcoords(m,l)%location = section

    successful = successful + 1
                
     

71  if (endcont .eqv. .false.) then
       reject = reject + 1
    end if
  end subroutine endmove
  
  
  subroutine comfind
    real :: dcomx,dcomy,dcomz
    t = time +1
    do m = 1,nprotein
       comx = 0.0
       comy = 0.0
       comz = 0.0
       do l = 2,maxlength
          
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
 
 
  subroutine radiusofgyration
    real :: totrog
real,dimension(:,:),allocatable :: rog
allocate(rog(nprotein,maxlength))
!check this
do m = 1,nprotein
   do l = 1, maxlength
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
    do m = 1, nprotein
       
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
    
    if(modulo(time,10) == 0) then
       write(29,*) nprotein*maxlength*(time-equilib),totrmsbrute
       !write(77,*) log(real(nprotein*maxlength*time)), log(totrmsbrute)
    end if
   end subroutine rms



   subroutine dataout

     !write(6,*) 'data'
     write(67,*) nprotein*maxlength
     write(67,*) ' '
     t = time 
     do m = 1,nprotein
        do l = 1,maxlength
           !maybe put in m
           write(67,*) 'C', 2*protcoords(m,l)%x, 2*protcoords(m,l)%y, &
           2*protcoords(m,l)%z
           !write(19,*) l, protcoords(t,m,l)%x, protcoords(t,m,l)%y, &
            !    protcoords(t,m,l)%z
        end do
     end do
   end subroutine dataout

   subroutine length
     real :: ddx,ddy,ddz,dxsum,dysum,dzsum
     t = time
     totchainlength = 0.0
     !actualchainlength = 0.0

     do m = 1,nprotein
        !chainlength = 0.0
        ddx = 0.0
        ddy = 0
        ddz = 0
        dxsum = 0.0
        dysum = 0.0
        dzsum = 0.0
        do l = 2,maxlength
           
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
     do t = 1,maxtime
!lengthvar =
     end do
   
   end subroutine error

 
  subroutine bonding
    integer :: bx,by,bz,f
    bx = 0
    by = 0
    bz = 0
    do m = 1,nprotein
       do l = 1,maxlength
          do f = 1,nprotein
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

  
  end subroutine bonding
 

  subroutine debug

    do m = 1,nprotein
       do l = 1,maxlength
          do f = 1,nprotein
             do g = 1,maxlength
                if(m ==f .and. l == g) then
                   continue
                else if (m/= f .or. l/= g) then
                   !write(6,*) protcoords(m,l)%x, protcoords(f,g)%x
                      if(protcoords(m,l)%x == protcoords(f,g)%x  .and. protcoords(m,l)%y == protcoords(f,g)%y  &
                           .and. protcoords(m,l)%z == protcoords(f,g)%z ) then
                         !write(6,*) 'fail is due to', l,g
                         fail = .true.
                         if (protcoords(m,l)%x == 0.0) write(6,*) 'x fail',l
                         if (protcoords(m,l)%y == 0.0) write(6,*) 'y fail',l
                          if (protcoords(m,l)%z == 0.0) write(6,*) 'z fail',l
                         !write(6,*) 'FAILLLLLLLLLLLLLLLLLL'
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
11       continue
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
