program setup
  implicit none


  integer :: maxno,protno,seed,natoms,maxsize,maxtime
  integer :: i,j,k,gridsize
  integer :: l,m,hold,z1,z2,z3,z4,z5,z6,control, n
  integer :: nprotein,split,section
  integer :: maxlength,g,xl,yl,zl
  real, external :: ran2
  real::x
  character(len = 10) :: commandread,commandread2
  
   type protein
     integer :: x,y,z
  end type protein
  type(protein),dimension(:,:),allocatable :: protcoords

     

  !open(17, file = 'setup.txt', action = 'read')
  open(19, file = 'initialtake2.xyz', action = 'write')
! write(6,*) seed,maxlength

!write(19,*) 'a'
  call read_setup
  !write(6,*) 'b'
  allocate(protcoords(nprotein,maxlength))
  call foundation
  !write(6,*) 'c'
  write(19,*) nprotein*maxlength
  write(19,*) ' '

  call positioning


contains





  
 subroutine initial
   character(len = 10) :: BIN
   read(17,*) BIN
   read(17,*) BIN !seed
   read(17,*) BIN
   read(17,*) BIN !maxlength
   read(17,*) BIN
   read(17,*) nprotein
   write(19,*) nprotein*maxlength
  ! write(23,*) ' ' 
   read(17,*) BIN
   read(17,*) gridsize
   natoms = gridsize**3
   !write(23,*) natoms
  ! write(23,*) ' '
   read(17,*) BIN
   read(17,*) maxtime
      read(17,*) BIN
   read(17,*) split
   write(19,*) ' '
 end subroutine initial

 subroutine foundation
    !type(protein),dimension(nprotein,maxlength) :: protcoords
   do m = 1,nprotein
      do l = 1,maxlength
      protcoords(m,l)%x = 0
   protcoords(m,l)%y = 0
    protcoords(m,l)%z = 0
      end do
   end do
   
   
 end subroutine foundation
 
 subroutine positioning
   integer :: gate,xbk,ybk,zbk,protpass,f
   Real:: t1,t2,t3,t4,t5,t6
   do m = 1,nprotein
      hold = m
      
      
      
      control = 1
      gate = 1
      protpass = 1
47    if (protpass ==1) then
         continue
      end if
     ! write(6,*) 'd'
35    if (gate == 1) then
         l = 1
         i = int(ran2(seed)*gridsize)+1
         j = int(ran2(seed)*gridsize)+1
         k = int(ran2(seed)*gridsize)+1
         !write(6,*) 'e'
         if (hold /= 1) then
            do g = 1, hold-1
               do f = 1,maxlength
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
         !xl = int((i*10)/split)-1
         !yl = int(((j*10)/split) -1)*N +1
         !zl = int(((k*10)/split) -1)*(N**2)
         !section = xl + yl + zl
         !write(6,*) 'f'
         protcoords(m,l)%x = i
         protcoords(m,l)%y = j
         protcoords(m,l)%z = k
         !protcoords(m,l)%location = section
         !write(19,*) 'O' ,    protcoords(m,l)%x,    protcoords(m,l)%y,  &
              !protcoords(m,l)%z
         
         do l = 2,maxlength

            xbk = protcoords(m,l-1)%x
            ybk = protcoords(m,l-1)%y
            zbk = protcoords(m,l-1)%z

            z1=1
            z2=1
            z3=1
            z4=1
            z5=1
            z6 = 1
            n=6
            
            control = 1
            x = ran2(seed)
            
15          if (control ==1) then

               i = xbk
               j = ybk
               k = zbk
               
               t1 = real(z1)/n
               t2 = real(z1+z2)/n
               t3 = real(z1 + z2 +z3)/n
               t4 = real(z1 + z2 + z3 + z4)/n
               t5 = real(z1 + z2 + z3 + z4 + z5)/n
               t6 = real(z1 + z2 + z3 + z4 + z5 + z6)/n
               if(x<= t1 .AND. z1 == 1) then
                  i = modulo(i,gridsize)+1
                  z1=0
               else if (x> t1 .AND. x<=t2 .AND. z2 == 1) then
                  i = modulo(i - 2,gridsize)+1
                  z2=0
                  
               else if (x> t2 .AND. x<=t3 .AND. z3 == 1) then
                  j = modulo(j,gridsize) +1
                  z3=0
               else if (x> t3 .AND. x<=t4 .AND. z4 == 1) then
                  
                  j = modulo(j -2,gridsize)+1
                  z4=0
               else if (x> t4 .AND. x<=t5 .AND. z5 == 1) then
                  k  = modulo(k,gridsize)+1
                  z5=0
               else if (x> t5 .AND. x<=1.0 .AND. z6 == 1) then
                   
                  k = modulo(k-2,gridsize)+1
                  z6=0
               end if
               do g = 1, hold
                  do f = 1,maxlength
                     if (i == protcoords(g,f)%x .and. j == protcoords(g,f)%y &
                          .and. k == protcoords(g,f)%z) then
                        n=n-1
                        if (n == 0) then
                           goto 47
                        else
                           goto 15
                        end if
                     end if
                  end do
               end do
               !end if
               !xl = int((i*10)/split)-1
               !yl = 1+int(((j*10)/split) -1)*N 
               !zl = int(((k*10)/split) -1)*(N**2)
               !section = xl + yl + zl
               !write(6,*) 'f'
               protcoords(m,l)%x = i
               protcoords(m,l)%y = j
               protcoords(m,l)%z = k
               !protcoords(m,l)%location = section               
               
               
               
               

                end if
             end do
          end if
          !end if
       

          do l = 1,maxlength
             write(19,*) m ,    protcoords(m,l)%x,    protcoords(m,l)%y,    protcoords(m,l)%z
          end do
       end do
     end subroutine positioning
     

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
   nprotein = 40                 !number of chains
   gridsize= 100        !dimensions of the lattice
   maxtime= 10000                !length of simulation
  
   
   !traj_file = 'trajectory'       ! Stem of name for trajectory files.
   !vgen = 1                       ! 1=random initial velocities (no momentum), 2=Maxwell-Boltzmann
   !v_init_file = ' '              ! Name of starting velocities file.
   !v_rescale = .FALSE.            ! If reading velocities from a file, .T.=rescale to get specified K.E.
   !wellstats_file = 'wellstats'   ! Stem of name for individual well statistics files.
   !well_tol = 1.0D-6              ! Energy tolerance for quenches to be considered the same.
   !xmol_type(1) = ''              ! Atom type for A (or charged) atoms in xmol dumps.
   !xmol_type(2) = ''              ! Atom type for B (or uncharged) atoms in xmol dumps.

   OPEN (UNIT=20, FILE='setup2.txt', STATUS='OLD', IOSTAT=err)

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
         CALL get_integer(maxlength)

      CASE ('NUM_CHAINS')
         CALL get_integer(nprotein)

      CASE ('LATTICE_DIMENSIONS')
         CALL get_integer(gridsize)
                  
      CASE ('SIM_TIME')
         CALL get_integer(maxtime)
         
           CASE ('EQUIB_TIME')
         !CALL get_integer(equilib)

      !CASE ('ENDM')
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

      CASE ('FILM')

      CASE ('DEBUGYES')

      CASE ('MAXPIV')

         CASE ('CRANKLIMIT')
         
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
