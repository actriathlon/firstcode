program setup
  implicit none


  integer :: maxno,protno,seed,natoms,maxsize,maxtime,maxlength1,maxlength2
  integer :: i,j,k,gridsize,nprotein1,nprotein2,nspecies
  integer :: l,m,hold,control,time,t
  integer :: nprotein,split,section,species,scans
  integer :: maxlength,g,xl,yl,zl,h,p
  real, external :: ran2
  real::x
  character(len = 10) :: commandread,commandread2
  integer,dimension(:,:),allocatable :: master,bonddd
  integer,dimension(:),allocatable :: specieslen,speciespop,link,chlen
  logical :: scalinginfo
   type protein
      integer :: x,y,z,species,type,linker
   end type protein

     type rprot
     real :: x,y,z
  end type rprot

  type(protein),dimension(:,:),allocatable :: protcoords
!linkerlength = 7
write(6,*) 'a'         

  !open(17, file = 'setup.txt', action = 'read')
  open(19, file = 'movetagged.vtf', action = 'write')
! write(6,*) seed,maxlength

write(6,*) 'a'
call read_setup
maxlength = maxval(specieslen)
  write(6,*) 'b'
  write(6,*) 'speciespop',speciespop(1)
  nprotein = sum(speciespop)
  write(6,*) 'proteins',nprotein
  allocate(protcoords(nprotein,maxlength))
  allocate(chlen(nprotein))
  allocate(bonddd(nprotein,maxlength))
  do h = 1,nspecies
     write(6,*) 'species length',specieslen(h)
     do p = 1,specieslen(h)
        write(6,*) 'type',h,p,master(h,p)
     end do
     end do
  
  call foundation
  call pdbsetup

do time = 1,maxtime
  call positioning
end do
     contains






  

 subroutine foundation
   integer::maxl,speciesselect
        integer,dimension(:),allocatable::spcount
   allocate(spcount(nspecies))
   
   write(6,*) maxlength,nprotein
   do m = 1,nprotein
      do l = 1,maxlength
      protcoords(m,l)%x = 0
   protcoords(m,l)%y = 0
    protcoords(m,l)%z = 0
      end do
   end do
  

   do m = 1,nspecies
spcount(m) = 0
      end do
 

   do m = 1,nprotein


32      speciesselect = int(ran2(seed)*nspecies)+1


      maxl = specieslen(speciesselect)
        chlen(m) = maxl

      if(spcount(speciesselect) == speciespop(speciesselect)) goto 32

        do l=1,maxl
                protcoords(m,l)%species = speciesselect
                protcoords(m,l)%type = master(speciesselect,l)
                protcoords(m,1)%linker = link(speciesselect)
        end do

        spcount(speciesselect)=spcount(speciesselect)+1

end do
   
 end subroutine foundation

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


  subroutine dataout
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit
    integer :: choice
    allocate(bondposit(nprotein,maxlength))
    !write(6,*) 'data'
    t = time
    choice = 19
    write(choice,*) ' '
    write(choice,*) 'timestep ', 'indexed'
    write(choice,*) 'pbc', real(4*gridsize),real(4*gridsize),real(4*gridsize)
    !write(67,*) ((nprotein1*maxlength1) + (nprotein2*maxlength2))
    acount = 0

    do m = 1,nprotein,1
       maxl = chlen(m)
       do l = 1,maxl,1
          write(choice,*) acount,real(4*protcoords(m,l)%x),real(4*protcoords(m,l)%y),real(4*protcoords(m,l)%z)
          acount = acount +1

          bondposit(m,l)%x = real(protcoords(m,l)%x)
          bondposit(m,l)%y = real(protcoords(m,l)%y)
          bondposit(m,l)%z = real(protcoords(m,l)%z)

          if(bonddd(m,l) == 1) bondposit(m,l)%x = modulo(bondposit(m,l)%x +0.5,real(gridsize))
          if(bonddd(m,l) == -1) bondposit(m,l)%x = modulo(bondposit(m,l)%x -0.5,real(gridsize))
          if(bonddd(m,l) == 2) bondposit(m,l)%y = modulo(bondposit(m,l)%y +0.5,real(gridsize))
          if(bonddd(m,l) == -2) bondposit(m,l)%y = modulo(bondposit(m,l)%y -0.5,real(gridsize))
          if(bonddd(m,l) == 3) bondposit(m,l)%z = modulo(bondposit(m,l)%z +0.5,real(gridsize))
          if(bonddd(m,l) == -3) bondposit(m,l)%z = modulo(bondposit(m,l)%z -0.5,real(gridsize))




          write(choice,*) acount, 4*bondposit(m,l)%x, 4*bondposit(m,l)%y, &
               4*bondposit(m,l)%z

71        continue
          acount = acount + 1
       end do
    end do



  end subroutine dataout



  subroutine pdbsetup
    integer::m,l,maxl,acount,dumres
    character(len=8) :: atnum
    acount = 0
    write(19,*) 'atom ', 'default ', 'radius ', 1.00000, 'name ','C'
    do m= 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          write(19,*) 'atom', acount, 'radius', 1.00000, 'name',(2*protcoords(m,1)%species), &
               'resid ',protcoords(m,l)%type
          write(19,*) 'atom', acount+1, 'radius', 1.00000, 'name',(2*protcoords(m,1)%species)-1,'resid ', &
               protcoords(m,l)%type
          acount = acount+2
       end do

    end do
    write(19,*) ' '
    acount = 0
    do m= 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl-1
          write(atnum,'(i7)') acount
          write(19,*) 'bond', adjustr(atnum)//':',acount+2
          write(19,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+2
       end do
       write(atnum,'(i7)') acount
       write(19,*) 'bond', adjustr(atnum)//':',acount+1
       acount = acount +2
    end do

  end subroutine pdbsetup





 subroutine positioning
   integer :: gate,xbk,ybk,zbk,protpass,f,speciesselect,spcount1,spcount2,maxl
   integer::maxback,dx,dy,dz
   Real:: t1,t2,t3,t4,t5,t6
   integer,dimension(:),allocatable::spcount
   allocate(spcount(nspecies))

   
   do m = 1,nprotein
      hold = m
      control = 1
      gate = 1
      protpass = 1
47    if (protpass ==1) then
         continue
      end if
      
      maxl = specieslen(protcoords(m,1)%species)

      
35    if (gate == 1) then
         l = 1
         i = int(ran2(seed)*gridsize)+1
         j = int(ran2(seed)*gridsize)+1
         k = int(ran2(seed)*gridsize)+1
         if (hold /= 1) then
            do g = 1, hold-1
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

      call bonddirection(bonddd(m,l))

         
         do l = 2,maxl

            xbk = protcoords(m,l-1)%x
            ybk = protcoords(m,l-1)%y
            zbk = protcoords(m,l-1)%z


            control = 1

            
15          if (control ==1) then


               dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
               dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))

               if(sqrt(real((dx**2) + (dy**2) + (dz**2))) > protcoords(m,1)%linker) goto 15
               

                     i = modulo(xbk+dx-1,gridsize)+1
                     j = modulo(ybk +dy-1,gridsize)+1
                     k  = modulo(zbk+dz-1,gridsize)+1

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
       
          
            end if
             end do
          end if

do l = 1,maxl,1
              call bonddirection(bonddd(m,l))
end do

       end do

        call dataout

     end subroutine positioning
     

SUBROUTINE read_setup

   USE keywords
   !USE setup_vars
   !USE stats_vars
   IMPLICIT NONE

   CHARACTER(LEN=20) :: keyword, option, argument
   CHARACTER(LEN=18), PARAMETER :: param_fmt1='(A16, 1PE20.10, A)'
   CHARACTER(LEN=16), PARAMETER :: param_fmt2='(A16, F16.10)'
   INTEGER :: err, i, j,runtype,typevar
   LOGICAL :: success
   !integer,dimension(:,:),intent(out),allocatable:: master
   !integer,dimension(:),intent(out),allocatable:: specieslen

   double precision :: intraenergy,interenergy

!  Defaults and descriptions.
   !  S.T.A.T. stands for short-time averaged temperature.


  !scalinginfo = .true.

   !seed = 6                      !seed for random number generation
   !maxlength = 20              !number of beads in chain
   !nprotein = 40                 !number of chains
   !gridsize= 100        !dimensions of the lattice
   !maxtime= 10000                !length of simulation
  
   
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
!maxlength = max(maxlength1,maxlength2)

         
      !CASE ('NUM_CHAINS')
         !CALL get_integer(nprotein1)

      !CASE ('NUM_CHAINS_2')
         !CALL get_integer(nprotein2)

      CASE ('TOTTYPES')
         CALL get_integer(nspecies)
         call get_integer(maxlength)
   write(6,*) 'bizarre',maxlength,nspecies
      !maxlength = 15

      CASE ('SCALEINFO')
         call get_logical(scalinginfo)

         !write(6,*) 'nspeciess',nspecies,maxlength
       CASE ('TOTCHAINS')
write(6,*) 'nspecies',nspecies,maxlength          
!CALL get_integer(nprotein)
          allocate(master(nspecies,maxlength))
          allocate(specieslen(nspecies))
          allocate(speciespop(nspecies))
      allocate(link(nspecies))    
       

   runtype = 1
       CASE ('NEWCHAIN')
          CALL get_integer(specieslen(runtype))
          if(specieslen(runtype)>maxlength) then
             write(6,*) 'FAIL: maxlength is too small--------------------'
             STOP 9
          end if
          CALL get_integer(speciespop(runtype))
call get_integer(link(runtype)) 
                    write(6,*) 'scalinginfo = ',scalinginfo
             if(scalinginfo .eqv. .false.) then
          write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
          do typevar = 1,specieslen(runtype)
             call get_integer(master(runtype,typevar))
             end do
             runtype = runtype+1
          else if(scalinginfo .eqv. .true.) then
             write(6,*) 'herererererere'
             do typevar = 1,specieslen(runtype)
                master(runtype,typevar) = 1
                write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
             end do
             runtype = runtype+1
end if
write(6,*) 'runtype =',runtype
if(runtype > nspecies+1) then
write(6,*) 'tottypes is too small'
STOP 9
end if

          


          
       CASE ('LATTICE_DIMENSIONS')
         CALL get_integer(gridsize)
                  
      CASE ('SIM_TIME')
         CALL get_integer(maxtime)
         
           CASE ('EQUIB_TIME')
         !CALL get_integer(equilib)

      CASE ('ENDM')
         !CALL get_integer(endm)

      CASE ('RIGHTANGLE')
         !CALL get_integer(right)

      CASE ('CRANK')
         !CALL get_integer(crank)

      CASE ('REPTATION')
         !CALL get_integer(rept)

      CASE ('PIVOT')
         !CALL get_integer(piv)

      CASE ('KT')
         !CALL get_dp(kT)

      !CASE ('ISBOND')
         
      CASE ('FILM')

      CASE ('DEBUGYES')

      CASE ('MAXPIV')

         CASE ('CRANKLIMIT')

     !CASE ('INTEREN')
     
     CASE ('INTRAEN')

CASE ('INTERACTION')
       
          
       CASE ('CLROTATION')
          !CALL get_logical(clr)
          
       CASE ('CLTRANSLATION')
          !CALL get_logical(clt)

          CASE ('NOINFO')


CASE ('RESTART')

CASE ('LINK')
!call get_integer(linkerlength)

CASE ('CLUSSTEP')

           CASE DEFAULT
         WRITE (6, '(/, 3A, /)') 'ERROR: Keyword not recognised in setup file: ', TRIM(keyword), '.'
         STOP

      END SELECT
   END DO
   CLOSE(20)


END SUBROUTINE read_setup
     
   end program setup





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
11        continue
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
