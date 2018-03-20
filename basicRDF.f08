program move

  implicit none

  character(len=10)::BIN
  integer,dimension(:),allocatable::type
  double precision,dimension(:),allocatable::coordx,coordy,coordz
  integer::totpoints,n,m,l,t,j,lx,timeofinterest,dummymax,count
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen
  double precision::dx,dy,dz,delta,r,dr,rho
  integer,dimension(:,:),allocatable::summing
  integer::binsize,maxl,gridsize,unusedtime,maxtime,nprotein,mtype,nspecies
  character(len = 10) :: commandread,commandread2,commandread3
character(len=10) :: a,b,c
real(8),  parameter :: PI_8  = 4 * atan (1.0_8)


  type rprot
     real :: x,y,z
  end type rprot

  type(rprot),dimension(:),allocatable :: com



  open(21, file = 'movetagged.vtf', action = 'read')
  open(23, file = 'coms.dat', action = 'read')
  open(29, file = 'RDF.dat', action = 'write')
  open(31, file = 'unnormRDF.dat', action = 'write')

  call read_setup
write(6,*) 'mytype',mtype
  nprotein = sum(speciespop)
  !mtype = 3
  totpoints = 0
write(6,*) 'nspecies',nspecies,speciespop(1),specieslen(1)
  do n=1,nspecies
     totpoints = totpoints + speciespop(n)*specieslen(n)
  end do
  totpoints = totpoints*2
  write(6,*) 'totpoints',totpoints,nprotein
  binsize = 10
  unusedtime = 100
  
call get_command_argument(1,commandread)
read(commandread, '(i3)') binsize
  call get_command_argument(2,commandread2)
read(commandread2, '(i4)') unusedtime
  call get_command_argument(3,commandread3)
read(commandread3, '(i4)') timeofinterest



  allocate(summing(mtype,4*gridsize*binsize))
  allocate(com(totpoints/2))
  allocate(coordx(totpoints/2))
  allocate(coordy(totpoints/2))
  allocate(coordz(totpoints/2))

  summing(:,:) = 0

    read(23,*) BIN,BIN,BIN,BIN
  do t = 1,timeofinterest
     read(23,*) BIN,com(t)%x,com(t)%y,com(t)%z
  end do

  allocate(type(totpoints/2))

count = 1

  read(21,*) BIN,BIN,BIN,BIN,BIN,BIN
  do m = 1,nprotein
     read(21,*)BIN,BIN,BIN,BIN,BIN,type(count),BIN,BIN
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
     !write(6,*) 'type',type(count),count
     maxl = specieslen(type(count)/2)
     !write(6,*) 'maxl',maxl
     count = count+1
     
     do l = 2,maxl
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,type(count)
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
     count = count+1
  end do

     end do


  dummymax = totpoints - nprotein
  do m = 1,dummymax
             read(21,*)BIN,BIN,BIN
     end do


     do t = 1,unusedtime
        read(21,*) BIN
        read(21,*) BIN
        !read(21,*) BIN,BIN,BIN
     do m = 1,totpoints
           read(21,*)BIN
        end do
  end do


  do t = 1,timeofinterest

   read(21,*) BIN
        read(21,*) BIN


        do m = 1,totpoints,1


        if(modulo(m,2) == 0) then
read(21,*) BIN,BIN,BIN,BIN
           go to 11 
end if
           read(21,*) j,coordx(j/2+1),coordy(j/2+1),coordz(j/2+1)


        dx = min(abs(com(t)%x-coordx((m+1)/2)),gridsize-abs(com(t)%x-coordx((m+1)/2)))
        dy = min(abs(com(t)%y-coordy((m+1)/2)),gridsize-abs(com(t)%y-coordy((m+1)/2)))
        dz = min(abs(com(t)%z-coordz((m+1)/2)),gridsize-abs(com(t)%z-coordz((m+1)/2)))

        delta = sqrt((dx**2)+(dy**2)+(dz**2))
        
        lx = INT(delta/(1.0/binsize))

        if(abs(lx)>(4*gridsize*binsize)) goto 11

        summing(type((m+1)/2)/2,lx) = summing(type((m+1)/2)/2,lx)+1

11      continue

     end do
  write(6,*) 'finished'
  end do

dr = 1.0/binsize
rho = (totpoints/2)/real(gridsize**3)

  do l = 1,4*gridsize*binsize
      !write(29,'(f8.3)',advance= 'no') real(l)/binsize
!do n = 1,mtype
   !if(n<mtype) write(29,'(f8.3)',advance= 'no') summing(n,l)
   !if(n==mtype) write(6,*) summing(n,l)
     !end do

   r = real(l)/binsize  
   write(29,*) r,summing(1,l)/(4*PI_8*(r**2)*rho),summing(2,l)/(4*PI_8*(r**2)*rho),summing(3,l)/(4*PI_8*(r**2)*rho)
         write(31,*) r,summing(1,l),summing(2,l),summing(3,l)
end do

contains

  SUBROUTINE read_setup

    USE keywords
    !USE setup_vars
    !USE stats_vars
    IMPLICIT NONE

    CHARACTER(LEN=20) :: keyword, option, argument
    CHARACTER(LEN=18), PARAMETER :: param_fmt1='(A16, 1PE20.10, A)'
    CHARACTER(LEN=16), PARAMETER :: param_fmt2='(A16, F16.10)'
    INTEGER :: err, i, j,runtype,dummy2,dummytype,xl,f,dum
    LOGICAL :: success
    double precision:: intra


    !  Defaults and descriptions.
    !  S.T.A.T. stands for short-time averaged temperature.



mtype = 0
 



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
          !CALL get_integer(seed)
          !seed = -ABS(seed)

          !CASE ('CHAIN_LENGTH')
          !CALL get_integer(maxlength1)

          !CASE('CHAIN_2_LENGTH')
          !call get_integer(maxlength2)

          !CASE ('NUM_CHAINS')
          !CALL get_integer(nprotein1)

          !CASE ('NUM_CHAINS_2')
          !CALL get_integer(nprotein2)
       CASE ('SCALEINFO')


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
          call get_integer(dum)   
          do f = 1,specieslen(runtype)
             call get_integer(dummytype)
             if(dummytype>mtype) mtype = dummytype
          end do
          write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
          runtype = runtype + 1

       CASE ('LATTICE_DIMENSIONS')
          CALL get_integer(gridsize)



          !allocate(tttt())

       CASE ('SIM_TIME')
          CALL get_integer(maxtime)
          maxtime = maxtime*1000


       CASE ('EQUIB_TIME')

       CASE ('RIGHTANGLE')
          !CALL get_integer(right)

       CASE ('CRANK')
          !CALL get_integer(crank)

       CASE ('REPTATION')


       CASE ('PIVOT')


       CASE ('KT')


       CASE ('FILM')

       CASE ('DEBUGYES')


          !CASE ('ISBOND')
          !CALL get_logical(isbond)

       CASE ('MAXPIV')
          !CALL get_integer(maxpiv)

       CASE ('CRANKLIMIT')
          !CALL get_integer(cranklimit)



       CASE ('INTRAEN')


       CASE ('INTERACTION')


       CASE ('CLROTATION')


       CASE ('CLTRANSLATION')


       CASE ('NOINFO')


       CASE ('RESTART')


       CASE ('LINK')

       CASE ('CLUSSTEP')


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



