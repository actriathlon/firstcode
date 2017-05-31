program move
  implicit none


  integer :: maxno,protno,seed,natoms,maxsize,maxtime,reptbackward,reptforward
  real :: i,j,k,gridsize,randomz,deltax1,deltax2,dx,dy,dz,msd
  integer :: l,m,hold,z1,z2,z3,z4,z5,z6,control
  real :: deltay1,deltay2,deltaz1,deltaz2,totrmsbrute
  integer :: nprotein,split,section,time,count
  integer :: maxlength,g,xl,yl,zl,t,r,N
  real, external :: ran2
  logical :: rightangle,endchain,run,pivotx,pivoty,pivotz,run2,rac,exist
  real :: totdisp,totchainlength,avechainlength,probs,comx,comy,comz,totrmserror
  real :: actualchainlength,totvar,deviate,totdeviate,varian,separation,lengthvar,random
  real,dimension(:),allocatable :: chainlength,variance,disp,rmserror,dxtot,dytot,dztot
  integer :: gate,xbk,ybk,zbk,protpass,f,scan,c
  Real:: t1,t2,t3,t4,t5,t6,di,dj,dk
  integer :: cont1,cont2,cont3,cont4,totcont,successful,reject
  character(len = 10) :: commandread,commandread2
!real :: testlength,testtotlength
 
  type protein
     real :: x,y,z
     integer :: location
  end type protein
  type(protein),dimension(:,:,:),allocatable :: protcoords
  type(protein),dimension(:,:),allocatable :: com
 
  open(17, file = 'setup2.txt', action = 'read')
  open(23, file = 'initialtake2.xyz', action = 'read')
  open(67, file = 'move.xyz', action = 'write')
  open(29, file = 'rms.dat', action = 'write')
  open(31, file = 'chainlength.dat', action = 'write')
  open(51, file = 'bonds.dat', action = 'write')
  open(91, file = 'avechainlength.dat', action = 'write')
  open(97, file = 'radiusofgyration.dat', action = 'write')
  open(77, file = 'logmsd.dat', action = 'write')
  
 totrmsbrute = 0.0
 
     successful = 0
     reject = 0
     call initial
     call get_command_argument(1,commandread)
     read(commandread, '(i3)') seed
          call get_command_argument(2,commandread2)
          read(commandread2, '(i4)') maxlength
          
     !reads in seed from commandline

  N = maxlength
  allocate(protcoords(maxtime+1,nprotein,maxlength))
  allocate(com(maxtime+1,nprotein))
  allocate(chainlength(nprotein))
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
     !call comfind
     !write(6,*) 'a'
     if (mod(time,10) == 0) then
        call dataout
     end if
     !write(6,*) 'b'
     call bonding
     !write(6,*) 'c'
     call length
     !write(6,*) 'd'
     
     !write(6,*) 'e'
     call positioning
     call comfind
     call rms
     call radiusofgyration
     
!write(6,*) 'g'
     write(6,*) 'step =',time
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
    read(17,*) maxtime
    read(17,*) BIN
    read(17,*) split
  end subroutine initial
 
  subroutine foundation
    !read in initial positions of the protein
    
    character(len = 10) :: BIN
    read(23,*) BIN
    write(6,*) maxlength
    do m = 1,nprotein
       do l = 1,maxlength
          read(23,*) BIN, protcoords(1,m,l)%x, protcoords(1,m,l)%y, protcoords(1,m,l)%z
          write(6,*) protcoords(1,m,l)%x
       end do
    end do
  
  end subroutine foundation
 
 
 
 
  subroutine positioning
    !randomly selects one of the 4 moves

    integer :: contchoose,a1,a2,a3,a4

           do m = 1,nprotein
              do l = 1,maxlength
             protcoords(time+1,m,l)%x = protcoords(time,m,l)%x
             protcoords(time+1,m,l)%y = protcoords(time,m,l)%y
             protcoords(time+1,m,l)%z = protcoords(time,m,l)%z
          end do
       end do
       t = time + 1
    do scan = 1, 4 !nprotein*maxlength
count = count + 1
     rac = .true.
       run = .true.
           run2 = .true. 
        


           randomz = ran2(seed)
           

!if (randomz <= 0.25) then
              
             !call endmove
          !else if (randomz > 0.25 .and. randomz <= 0.5) then
                                !call crankshaftmove
           
     !else if (randomz > 0.5 .and. randomz <= 0.75) then
          !call rightanglemove
        !else if (randomz > 0.75 .and. randomz <= 1.0) then
       
   !call reptation
   
             ! end if
               call pivot   
    end do
  end subroutine positioning


  subroutine crankshaftmove
    !performs crankshaft move
logical :: crankcont
    !this means that m/=0 and is up to nprotein
m =modulo(int(ran2(seed)*nprotein),nprotein)+1
!this means that l cannot be maxlength,maxlength - 1 or l
l = modulo(int(ran2(seed)*maxlength-3),maxlength-3)+2

probs = ran2(seed)
              pivotx = .false.
          pivoty = .false.
          pivotz = .false.

          
          
               !check to see if 4 atoms are in the same plane
             if (protcoords(t,m,l-1)%x == protcoords(t,m,l)%x &
                  .and.  protcoords(t,m,l-1)%x == protcoords(t,m,l+1)%x &
                  .and.protcoords(t,m,l-1)%x == protcoords(t,m,l+2)%x) then
                pivotx = .true.                  
             else if (protcoords(t,m,l-1)%y == protcoords(t,m,l)%y .and. &
                  protcoords(t,m,l-1)%y == protcoords(t,m,l+1)%y &
                  .and.protcoords(t,m,l-1)%y == protcoords(t,m,l+2)%y ) then
                pivoty = .true.                               
             else if (protcoords(t,m,l-1)%z == protcoords(t,m,l)%z &
                  .and. protcoords(t,m,l-1)%z == protcoords(t,m,l+1)%z &
                  .and.protcoords(t,m,l-1)%z == protcoords(t,m,l+2)%z ) then
                pivotz = .true.
             end if


             if (pivotx .eqv. .false. .and. pivoty .eqv. .false. .and. pivotz .eqv. .false.) then
                crankcont = .false.
                goto 43
             end if
             
crankcont = .true.

if (pivotx .eqv. .true.) then

   
             if (mod(protcoords(t,m,l)%y - protcoords(t,m,l-1)%y,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%y - protcoords(t,m,l+2)%y,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%y - protcoords(t,m,l+1)%y,gridsize) ==0.0) then
                dy = 1.0
                dz = 0.0 
             else  if (mod(protcoords(t,m,l)%z - protcoords(t,m,l-1)%z,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%z - protcoords(t,m,l+2)%z,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%z - protcoords(t,m,l+1)%z,gridsize) ==0.0) then
                dy = 0.0
                dz = 1.0
             else
                crankcont = .false.
                goto 43
             end if
           
             if (probs <= (1.0/3)) then
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-2*dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == protcoords(t,m,l)%x .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-2*dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-2*dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == protcoords(t,m,l+1)%x .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-2*dz-1,gridsize)+1) then
                        crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-2*dy-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-2*dz-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-2*dy-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-2*dz-1,gridsize)+1
              
                !xl = int((protcoords(t,m,l)%x*10)/split)-1
                !yl = 1+int(((protcoords(t,m,l)%y*10)/split) -1)*N
                !zl = int(((protcoords(t,m,l)%z*10)/split) -1)*(N**2)
                !section = xl + yl + zl
                !protcoords(t,m,l)%location = section

                successful = successful + 1
              
             else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-2,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-2,gridsize)+1
              

              successful = successful +1 
             else if (probs > (2.0/3) .and. probs <= (1.0)) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x,gridsize)+1
  
              successful = successful + 1
             end if
          end if

          
!!!!!!!!!!pivot y
        
       if (pivoty .eqv. .true.) then
             
             if (mod(protcoords(t,m,l)%x - protcoords(t,m,l-1)%x,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%x - protcoords(t,m,l+2)%x,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%x - protcoords(t,m,l+1)%x,gridsize) ==0.0) then
                dx = 1.0
                dz = 0.0 
             else  if (mod(protcoords(t,m,l)%z - protcoords(t,m,l-1)%z,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%z - protcoords(t,m,l+2)%z,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%z - protcoords(t,m,l+1)%z,gridsize) ==0.0) then
                dx = 0.0
                dz = 1.0
             else
                crankcont = .false.
                goto 43

             end if
           
             if (probs <= (1.0/3)) then
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-2*dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == protcoords(t,m,l)%y .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-2*dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-2*dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == protcoords(t,m,l+1)%y .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-2*dz-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-2*dx-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-2*dz-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-2*dx-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-2*dz-1,gridsize)+1
              
              successful = successful +1
             else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-2,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-2,gridsize)+1
              
  
              successful = successful + 1
             else if (probs > (2.0/3) .and. probs <= 1.0) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-dz-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-dz-1,gridsize)+1
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y,gridsize)+1
              

                 successful = successful +1
             end if
          end if
        
!!!!pivot z
        
        
    if (pivotz .eqv. .true.) then




             
             if (mod(protcoords(t,m,l)%y - protcoords(t,m,l-1)%y,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%y - protcoords(t,m,l+2)%y,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%y - protcoords(t,m,l+1)%y,gridsize) ==0.0) then
                dy = 1.0
                dx = 0.0 
             else  if (mod(protcoords(t,m,l)%x - protcoords(t,m,l-1)%x,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l+1)%x - protcoords(t,m,l+2)%x,gridsize) ==1.0 .and. &
                  mod(protcoords(t,m,l)%x - protcoords(t,m,l+1)%x,gridsize) ==0.0) then
                dy = 0.0
                dx = 1.0
             else
                crankcont = .false.
                goto 43
             end if
           
             if (probs <= (1.0/3)) then
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-2*dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == protcoords(t,m,l)%z .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-2*dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-2*dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == protcoords(t,m,l+1)%z .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-2*dx-1,gridsize)+1) then
                        crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-2*dy-1,gridsize)+1
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-2*dx-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-2*dy-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-2*dx-1,gridsize)+1

              successful = successful + 1
             else if (probs > (1.0/3) .and. probs <= (2.0/3)) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z-2,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1) then
                         crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z-2,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z-2,gridsize)+1
              
                successful = successful + 1
             else if (probs > (2.0/3) .and. probs <= 1.0) then
              
                do g = 1,nprotein
                   do f = 1,maxlength
                      if(protcoords(t,g,f)%y == modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l)%z,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,m,l+1)%z,gridsize)+1 .and. &
                           protcoords(t,g,f)%x == modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1) then
crankcont = .false.
                         goto 43
                      else
                         continue
                      end if
                   end do
                end do
              
                protcoords(t,m,l)%y = modulo(protcoords(t,m,l)%y-dy-1,gridsize)+1
                protcoords(t,m,l)%x = modulo(protcoords(t,m,l)%x-dx-1,gridsize)+1
                protcoords(t,m,l+1)%y = modulo(protcoords(t,m,l+1)%y-dy-1,gridsize)+1
                protcoords(t,m,l+1)%x = modulo(protcoords(t,m,l+1)%x-dx-1,gridsize)+1
                protcoords(t,m,l)%z = modulo(protcoords(t,m,l)%z,gridsize)+1
                protcoords(t,m,l+1)%z = modulo(protcoords(t,m,l+1)%z,gridsize)+1
              
                successful = successful + 1
             end if
          end if
       

43     if (crankcont .eqv. .false.) then
          reject = reject + 1
          end if


        
  end subroutine crankshaftmove

  subroutine rightanglemove

    !performs a 180 flip on an atoms

    m =modulo(int(ran2(seed)*nprotein),nprotein)+1
    ! l cannot be equal to 1 or maxlength
    l = modulo(int(ran2(seed)*maxlength-2),maxlength-2)+2

    rac = .true.

       if (abs(mod(protcoords(t,m,l+1)%x - protcoords(t,m,l-1)%x,gridsize)) /= 2.0 &
            .and. abs(mod(protcoords(t,m,l+1)%y - protcoords(t,m,l-1)%y,gridsize)) /= 2.0 &
            .and. abs(mod(protcoords(t,m,l+1)%z - protcoords(t,m,l-1)%z,gridsize)) /= 2.0 ) then
          rightangle = .true.
       else 
          rightangle = .false.
       end if

       if (rightangle .eqv. .true.) then    
              
 deltax1 = mod(protcoords(t,m,l+1)%x - protcoords(t,m,l)%x,gridsize)
             deltax2 =mod(protcoords(t,m,l-1)%x - protcoords(t,m,l)%x,gridsize)
             deltay1=mod(protcoords(t,m,l+1)%y - protcoords(t,m,l)%y,gridsize)
             deltay2 =mod(protcoords(t,m,l-1)%y - protcoords(t,m,l)%y,gridsize)
             deltaz1 =mod(protcoords(t,m,l+1)%z - protcoords(t,m,l)%z,gridsize)
             deltaz2= mod(protcoords(t,m,l-1)%z - protcoords(t,m,l)%z,gridsize)
             
             di = 0
             dj = 0
             dk = 0
           
             di = (deltax1 + deltax2)
             dj = (deltay1 + deltay2)
             dk = (deltaz1 + deltaz2)

             
           
             if (abs(di) <= 1 .or. abs(dj) <= 1 .or. abs(dk) <= 1) then
                
             i = modulo(protcoords(t,m,l)%x+di-1,gridsize)+1
             j = modulo(protcoords(t,m,l)%y+dj-1,gridsize)+1
             k = modulo(protcoords(t,m,l)%z+dk-1,gridsize)+1
           
           
             do g = 1, nprotein
                do f = 1,maxlength
                   if (i == protcoords(t,g,f)%x .and. j == protcoords(t,g,f)%y &
                        .and. k == protcoords(t,g,f)%z) then

                      rac = .false.
                      goto 37              
                   else
                      continue                
                   end if
                end do
             end do
           
             protcoords(t,m,l)%x = i
             protcoords(t,m,l)%y = j
             protcoords(t,m,l)%z = k

           
 

             successful = successful + 1
          end if
          end if
37 if (rac .eqv. .false. .or. rightangle .eqv. .false.) then

   reject = reject + 1
   end if

        
 end subroutine rightanglemove

  subroutine pivot
                integer :: b
                logical :: pivcont
                real :: choose3

                write(6,*) 'here comes the pivot!!'
                
                choose3 = ran2(seed)
                pivcont = .true. 
                m =modulo(int(ran2(seed)*nprotein),nprotein)+1

               l = modulo(int(ran2(seed)*maxlength-3),maxlength-3)+2
               write(6,*) 'm = ', m
                              write(6,*) 'l =', l
                t = time + 1


                
                if (l> maxlength/2) then
                   if (choose3 <= 1.0/5) then
                   do b = l+1,maxlength 
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%y + &
                                    mod(protcoords(t,m,l)%x-protcoords(t,m,b)%x,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do

                      do b = l+1,maxlength 
                         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,l)%y &
                              -protcoords(t,m,b)%y,gridsize),gridsize)
                         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y + mod(protcoords(t,m,l)%x &
                              -protcoords(t,m,b)%x,gridsize),gridsize)
                         !protcoords(t,m,b)%z = protcoords(t,m,b)%z 
                         successful = successful + 1
                      end do
                      
                      else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then 
                      
                      do b = l+1,maxlength 
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%y-protcoords(t,m,l)%y,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%y - &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b =l+1,maxlength 
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%y- &
              protcoords(t,m,l)%y,gridsize),gridsize)
         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do

                      else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then

                      do b = l+1,maxlength 
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%y - &
                                    mod(protcoords(t,m,b)%y-protcoords(t,m,l)%y,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b =l+1,maxlength                 
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y - mod(protcoords(t,m,b)%y- &
              protcoords(t,m,l)%y,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do
                      else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then

                      do b = l+1,maxlength 
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%z-protcoords(t,m,l)%z,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%z + &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%z &
                                    .and. protcoords(t,m,b)%y == protcoords(t,f,g)%y) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b =l+1,maxlength 
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%z- &
              protcoords(t,m,l)%z,gridsize),gridsize)
         protcoords(t,m,b)%z = modulo(protcoords(t,m,b)%z + mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do
                      else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then

                      do b = l+1,maxlength 
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%z-protcoords(t,m,l)%z,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%z - &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%z &
                                    .and. protcoords(t,m,b)%y == protcoords(t,f,g)%y) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b = l+1,maxlength                
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%z- &
              protcoords(t,m,l)%z,gridsize),gridsize)
         protcoords(t,m,b)%z = modulo(protcoords(t,m,b)%z - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do
           end if
                      

else if (l< maxlength/2) then


                         if (choose3 <= 1.0/5) then
                      do b = l-1,1,-1
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%y-protcoords(t,m,l)%y,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%y + &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do

                      do b = l-1,1,-1
                         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%y- &
                              protcoords(t,m,l)%y,gridsize),gridsize)
                         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y + mod(protcoords(t,m,b)%x- &
                              protcoords(t,m,l)%x,gridsize),gridsize)
                         !protcoords(t,m,b)%z = protcoords(t,m,b)%z 
                         successful = successful + 1
                      end do
                      
                                            else if (choose3 > 1.0/5 .and. choose3 <= 2.0/5) then
                      
                      do b = l-1,1,-1
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%y-protcoords(t,m,l)%y,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%y - &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b = l-1,1,-1                
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%y- &
              protcoords(t,m,l)%y,gridsize),gridsize)
         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
         !protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do

                      else if (choose3 > 2.0/5 .and. choose3 <= 3.0/5) then

                      do b = l-1,1,-1
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and.  modulo(protcoords(t,m,b)%y - &
                                    mod(protcoords(t,m,b)%y-protcoords(t,m,l)%y,gridsize),gridsize) == protcoords(t,f,g)%y &
                                    .and. protcoords(t,m,b)%z == protcoords(t,f,g)%z) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b = l-1,1,-1                
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
         protcoords(t,m,b)%y = modulo(protcoords(t,m,b)%y - mod(protcoords(t,m,b)%y- &
              protcoords(t,m,l)%y,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do

                      else if (choose3 > 3.0/5 .and. choose3 <= 4.0/5) then
                      do b = l-1,1,-1
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%z-protcoords(t,m,l)%z,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%z + &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%z &
                                    .and. protcoords(t,m,b)%y == protcoords(t,f,g)%y) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b = l-1,1,-1                
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x - mod(protcoords(t,m,b)%z- &
              protcoords(t,m,l)%z,gridsize),gridsize)
         protcoords(t,m,b)%z = modulo(protcoords(t,m,b)%z + mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do
                      else if (choose3 > 4.0/5 .and. choose3 <= 1.0) then

                      do b = l-1,1,-1
                         do f = 1,nprotein
                            do g = 1,maxlength
                               if(modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%z-protcoords(t,m,l)%z,gridsize),gridsize) &
                                    == protcoords(t,f,g)%x .and. modulo(protcoords(t,m,b)%z - &
                                    mod(protcoords(t,m,b)%x-protcoords(t,m,l)%x,gridsize),gridsize) == protcoords(t,f,g)%z &
                                    .and. protcoords(t,m,b)%y == protcoords(t,f,g)%y) then
                                  pivcont = .false.
                                  goto 75
                               end if
                            end do
                         end do
                      end do
      do b = l-1,1,-1                
         protcoords(t,m,b)%x = modulo(protcoords(t,m,b)%x + mod(protcoords(t,m,b)%z- &
              protcoords(t,m,l)%z,gridsize),gridsize)
         protcoords(t,m,b)%z = modulo(protcoords(t,m,b)%z - mod(protcoords(t,m,b)%x- &
              protcoords(t,m,l)%x,gridsize),gridsize)
!protcoords(t,m,b)%z = protcoords(t,m,b)%z + mod(protcoords(t,m,l)%y-protcoords(t,m,b)%y,gridsize
successful = successful + 1
end do
end if

                   end if

75                 if (pivcont .eqv. .false.) then
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
m =modulo(int(ran2(seed)*nprotein),nprotein)+1
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
                if (choose > 0) then
                dxx = mod(protcoords(t,c,2)%x - protcoords(t,c,1)%x,gridsize)
                dyx = mod(protcoords(t,c,2)%y - protcoords(t,c,1)%y,gridsize)
                dzx = mod(protcoords(t,c,2)%z - protcoords(t,c,1)%z,gridsize)

                if (dxx == 1.0) g1 = 0
                if (dxx == -1.0) g2 = 0
                if (dyx == 1.0) g3 = 0
                if (dyx == -1.0) g4 = 0
                if (dzx == 1.0) g5 =0
                if (dzx == 1.0) g6 = 0

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
                      if(protcoords(t,g,f)%x == modulo(protcoords(t,c,1)%x + dx -1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,c,1)%y + dy -1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,c,1)%z + dz -1,gridsize)+1) then
                         
                         reptcont = .false.
                         goto 83
                      else
                         continue
                      end if
                   end do
                end do
              
                do l =1,maxlength-1
                   protcoords(t,c,maxlength-l+1)%x = protcoords(t,c,maxlength-l)%x
                   protcoords(t,c,maxlength-l+1)%y = protcoords(t,c,maxlength-l)%y
                   protcoords(t,c,maxlength-l+1)%z = protcoords(t,c,maxlength-l)%z
                 
!successful = successful+1
                end do
              
                protcoords(t,c,1)%x = modulo(protcoords(t,c,1)%x + dx,gridsize)
                protcoords(t,c,1)%y = modulo(protcoords(t,c,1)%y + dy,gridsize)
                protcoords(t,c,1)%z = modulo(protcoords(t,c,1)%z + dz,gridsize)
              
  
                successful = successful + 1
reptforward = reptforward + 1                
             else if (choose < 0) then



                                dxx = mod(protcoords(t,c,maxlength-1)%x - protcoords(t,c,maxlength)%x,gridsize)
                dyx = mod(protcoords(t,c,maxlength-1)%y - protcoords(t,c,maxlength)%y,gridsize)
                dzx = mod(protcoords(t,c,maxlength-1)%z - protcoords(t,c,maxlength)%z,gridsize)

                if (dxx == 1.0) g1 = 0
                if (dxx == -1.0) g2 = 0
                if (dyx == 1.0) g3 = 0
                if (dyx == -1.0) g4 = 0
                if (dzx == 1.0) g5 =0
                if (dzx == 1.0) g6 = 0

                direc = ran2(seed)
                if(direc <= (1.0*g1)/5 .and. g1 == 1) then
                   dx = 1.0
                else if(direc <= (1.0*g1 +g2)/5 .and.direc > (1.0*g1)/5 .and.  g2 == 1) then
                   dx = -1.0
                else if(direc <= (1.0*g1+g2+g3)/5 .and.direc > (1.0*g1+g2)/5 .and.  g3 == 1) then
                   dy = 1.0
                                      
                else if(direc <= (1.0*g1+g2+g3+g4)/5 .and.direc > (1.0*g1+g2+g3)/5 .and.  g4 == 1) then
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
                      if(protcoords(t,g,f)%x == modulo(protcoords(t,c,maxlength)%x + dx -1,gridsize)+1 .and. &
                           protcoords(t,g,f)%y == modulo(protcoords(t,c,maxlength)%y + dy -1,gridsize)+1 .and. &
                           protcoords(t,g,f)%z == modulo(protcoords(t,c,maxlength)%z + dz -1,gridsize)+1) then
                         reptcont = .false.
                                                 
                         goto 83
                      else
                         continue
                      end if
                   end do
                end do
                                        
              
                                         do l =1,maxlength-1
                                            protcoords(t,c,l)%x = protcoords(t,c,l+1)%x
                                            protcoords(t,c,l)%y = protcoords(t,c,l+1)%y
                                            protcoords(t,c,l)%z = protcoords(t,c,l+1)%z
 
                                         end do

              
                protcoords(t,c,maxlength)%x = modulo(protcoords(t,c,maxlength)%x + dx,gridsize)
                protcoords(t,c,maxlength)%y = modulo(protcoords(t,c,maxlength)%y + dy,gridsize)
                protcoords(t,c,maxlength)%z = modulo(protcoords(t,c,maxlength)%z + dz,gridsize)
              


                successful = successful + 1
reptbackward = reptbackward + 1


                
             end if
83           if (reptcont .eqv. .false.) then
                reject = reject + 1
                end if

  end subroutine reptation
 
  subroutine endmove
    !performs end chain rotation
    real :: choose2
    logical :: endcont
    
        m =modulo(int(ran2(seed)*nprotein),nprotein)+1
    !if(m == 0) m = 1
        random = ran2(seed)
        choose2 = ran2(seed)-0.5
           

           endcont = .true.

           if (choose2 >= 0.0) then
              l = 1
                   i = protcoords(t,m,l+1)%x
                   j = protcoords(t,m,l+1)%y
                   k = protcoords(t,m,l+1)%z
                else if (choose2 <= 0.0) then
                   l = maxlength
                   i = protcoords(t,m,l-1)%x
                   j = protcoords(t,m,l-1)%y
                   k = protcoords(t,m,l-1)%z
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
                      if (i == protcoords(t,g,f)%x .and. j == protcoords(t,g,f)%y &
                           .and. k == protcoords(t,g,f)%z ) then
endcont = .false.
                            goto 71 
                      else
                         continue
                      end if
                   end do
                end do
           
    
              
                protcoords(t,m,l)%x = i
                protcoords(t,m,l)%y = j
                protcoords(t,m,l)%z = k
                protcoords(t,m,l)%location = section

                successful = successful + 1
            
     

71        if (endcont .eqv. .false.) then
                   reject = reject + 1
                   end if
  end subroutine endmove
 
 
  subroutine comfind
    t = time +1
    do m = 1,nprotein
       comx = 0.0
       comy = 0.0
       comz = 0.0
       do l = 2,maxlength
          comx = comx + mod(protcoords(t,m,l)%x-protcoords(t,m,l-1)%x,gridsize)
          comy = comy + mod(protcoords(t,m,l)%y-protcoords(t,m,l-1)%y,gridsize)
          comz = comz + mod(protcoords(t,m,l)%z-protcoords(t,m,l-1)%z,gridsize)
 
       end do
       com(t,m)%x = modulo(protcoords(t,m,1)%x +comx,gridsize)
       com(t,m)%y = modulo(protcoords(t,m,1)%y +comy,gridsize)
       com(t,m)%z = modulo(protcoords(t,m,1)%z +comz,gridsize)

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
      rog(m,l) = min((mod(protcoords(time,m,l)%x-com(time,m)%x,gridsize))**2,(gridsize - &
           mod(protcoords(time,m,l)%x-com(time,m)%x,gridsize))**2)+ &
           min((mod(protcoords(time,m,l)%y-com(time,m)%y,gridsize))**2,(gridsize - &
           mod(protcoords(time,m,l)%y-com(time,m)%y,gridsize))**2)+ &
           min((mod(protcoords(time,m,l)%z-com(time,m)%z,gridsize))**2,(gridsize - &
           mod(protcoords(time,m,l)%z-com(time,m)%z,gridsize))**2)
      totrog = totrog + rog(m,l)
   end do
end do



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

       msdx(m) = min((mod(com(t,m)%x - com(t-1,m)%x,gridsize))**2, &
            (gridsize - mod(com(t,m)%x -com(t-1,m)%x,gridsize))**2)
       !write(6,*) 'x =',msdx(m)

       !write(6,*) 'y sortr = ',com(t,m)%y,com(t-1,m)%y
       msdy(m) = min((mod(com(t,m)%y -com(t-1,m)%y,gridsize))**2, &
            (gridsize - mod(com(t,m)%y -com(t-1,m)%y,gridsize))**2)
          !write(6,*) 'y t = ', com(t,m)%y,'y t-1 = ', com(t-1,m)%y
       !write(6,*) 'y =',msdy(m)
       msdz(m) = min((mod(com(t,m)%z - com(t-1,m)%z,gridsize))**2, &
            (gridsize - mod(com(t,m)%z -com(t-1,m)%z,gridsize))**2)
       
       !msdz(m) = mod(com(t,m)%z - com(t-1,m)%z,gridsize)
              !write(6,*) 'z =',msdz(m)
       totrmsbrute = totrmsbrute + msdx(m) + msdy(m) + msdz(m)
   end do
    msd = sum(msdx) + sum(msdy) + sum(msdz)
    write(6,*) 'msd =', msd
    write(6,*) 'msdsum = ', msdsum !this is wrong
    write(6,*) 'msd brute' , totrmsbrute
       msdsum = msdsum + msd
       write(29,*) time,totrmsbrute
       write(77,*) log(real(nprotein*maxlength*time)), log(totrmsbrute)
 
   end subroutine rms



   subroutine dataout

     !write(6,*) 'data'
     write(67,*) nprotein*maxlength
     write(67,*) ' '
     t = time 
     do m = 1,nprotein
        do l = 1,maxlength
           !maybe put in m
           write(67,*) 'C', 2*protcoords(t,m,l)%x, 2*protcoords(t,m,l)%y, &
           2*protcoords(t,m,l)%z
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
           
           ddx = mod(protcoords(t,m,l)%x-protcoords(t,m,l-1)%x,gridsize)
           if (ddx > 1.5) ddx =1
           if (ddx < -1.5) ddx = -1
               

           dxsum = dxsum + ddx
           ddy = mod(protcoords(t,m,l)%y-protcoords(t,m,l-1)%y,gridsize)
           if (ddy > 1.5) ddy =1
           if (ddy < -1.5) ddy = -1
           
           dysum = dysum + ddy
           ddz = mod(protcoords(t,m,l)%z-protcoords(t,m,l-1)%z,gridsize)
           if (ddz > 1.5) ddz =1
           if (ddz < -1.5) ddz = -1
           
           dzsum = dzsum + ddz
           
           !write(6,*) 'ddx =', ddx, 'ddy = ', ddy, 'ddz = ', ddz
        !chainlength(m) = SQRT((min, abs(gridsize- &
            ! modulo(protcoords(t,m,1)%x - protcoords(t,m,maxlength)%x,gridsize))))**2 + &
             !(min(modulo(protcoords(t,m,1)%y - protcoords(t,m,maxlength)%y,gridsize), abs(gridsize - &
             !modulo(protcoords(t,m,1)%y - protcoords(t,m,maxlength)%y,gridsize))))**2 + &
             !(min(modulo(protcoords(t,m,1)%z - protcoords(t,m,maxlength)%z,gridsize), abs(gridsize - &
            !modulo(protcoords(t,m,1)%z - protcoords(t,m,maxlength)%z,gridsize))))**2)
end do
chainlength(m) = ((dxsum**2) + (dysum**2) + (dzsum**2))
        totchainlength = totchainlength + chainlength(m)
     end do
             !totchainlength = sum(chainlength)
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
                if (modulo(protcoords(t,m,l)%x - protcoords(t,f,g)%x,gridsize) < separation &
                     .and. modulo(protcoords(t,m,l)%y - protcoords(t,f,g)%y,gridsize) < separation &
                     .and. modulo(protcoords(t,m,l)%z - protcoords(t,f,g)%z,gridsize) < separation) then
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
    write(51,*) time,bz

  
  end subroutine bonding
 
 
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
