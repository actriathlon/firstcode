program corr
  implicit none
  integer ::a,b,c,d,num_calc,interval,maxint
  integer::nprotein,time,choose
  double precision :: del,sum_corr,su
  double precision :: sig
  integer,dimension(:,:),allocatable::correlation
  double precision,dimension(:),allocatable::averagevals,totcorr,abar,sigma,stdev
  character(len = 10) :: commandread


  call get_command_argument(1,commandread)
  read(commandread, '(i3)') nprotein
  call get_command_argument(2,commandread)
  read(commandread, '(i4)') maxint
  call get_command_argument(3,commandread)
  read(commandread, '(i5)') time
  call get_command_argument(4,commandread)
  read(commandread, '(i1)') choose


  nprotein = 2*nprotein


  allocate(averagevals(nprotein))
  allocate(correlation(time,nprotein))
  allocate(totcorr(maxint+1))
  allocate(abar(maxint+1))
  allocate(sigma(maxint+1))
  allocate(stdev(maxint+1))

  if(choose == 0) then
     open(12,file ='onecorrelation.dat',action='read')
     open(21,file='onecorrelationplot.dat',action='write')
  else if (choose == 1) then
     open(12,file='clustercorrelation.dat',action='read')
     open(21,file='clustercorrelationplot.dat',action = 'write')
  end if

  write(6,*) 'time',time,'maxint',maxint
  write(6,*) 'a'

  do a = 1,time
     read(12,*) correlation(a,:)
  end do

  !write(6,*) correlation(1,:)
  !write(6,*) correlation(2,:)
  abar(:) = 0.0d0
  write(6,*) 'b'

  averagevals(:) = 0.0d0

  do b =1,nprotein,1
     su = 0.0d0
     do a =1,time,1
        su = su +correlation(a,b)
        !averagevals(b) = sum(correlation(:,b))/1000
     end do
     !averagevals(b) = su/time 
  end do

  write(6,*) 'c' !,averagevals(:)

  totcorr(:) = 0.0d0
  sigma(:) = 0.0d0
  stdev(:) = 0.0d0
  do b = 1,nprotein,1
     do interval = 1,maxint+1,1
        num_calc = 0
        sum_corr = 0.0d0
        do a = 1,time-(interval-1),1
           del = ((correlation(a,b))-averagevals(b))*((correlation(a+(interval-1),b))-averagevals(b))
           num_calc = num_calc+1
           sum_corr = sum_corr + del
        end do
        abar(interval) = abar(interval) + sum_corr
        totcorr(interval) = totcorr(interval) + (sum_corr)
     end do
  end do




  do interval = 1,maxint+1,1
     abar(interval) = abar(interval)/(nprotein*(time-(interval-1)))
  end do

  !write(6,*) 'abars',abar(:)

  do b = 1,nprotein,1
     do interval = 1,maxint+1,1
        sig = 0.0d0
        do a = 1,time-(interval-1)
           del =(correlation(a,b)-averagevals(b))*(correlation(a+(interval-1),b)-averagevals(b))
           sig = sig + del
        end do
        sig = sig/(time-(interval-1))
        sigma(interval) = sigma(interval) + ((sig-abar(interval))**2) 
     end do
  end do

  do interval = 1,maxint+1,1
     stdev(interval) = sqrt(sigma(interval)/(nprotein-1))
  end do


if(choose ==1) then
  do interval = 1,maxint+1,1
     write(21,*)interval-1,totcorr(interval)/(nprotein*(time-(interval-1))),stdev(interval)
  end do
else if(choose == 0) then
do interval = 1,maxint+1,1
     write(21,*)interval-1,totcorr(interval)/(0.5*nprotein*(time-(interval-1))),stdev(interval)
  end do
   end if

end program corr
