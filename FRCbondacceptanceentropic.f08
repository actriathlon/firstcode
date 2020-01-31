program move

  implicit none

  double precision :: totalenergy,kT,totrmsbrute!,randomz
  double precision :: totdisp,totchainlength,avechainlength,choose2
  double precision :: actualchainlength,random,totiso
  double precision :: runningaveROG,runningaveEtE,chainlength,sumdebug,testdummy,polymerrog,polymerl
  integer :: rept,datayes,maxclussize,minl,outputrate,backtime
  integer ::lastruntime,timebase,steplength,macc,mrej,roatt,roacc
  integer :: gr,timeofinterest,totpoints,binsize,phasecount
  integer :: gridsize,maxlength,t,N,nprotein,time,debugging,count,avecluspop,aveclusnumber
  integer :: seed,natoms,maxtime,reptbackward,reptforward,equilib,cltacc,cltatt
  integer :: successful,reject,maxlength1,maxlength2,q,d,u,sm,qt,ft,gt,scans
  integer:: rerej,reacc,nprotein1,nprotein2,nspecies,mtype,sumr1,sumr2,yl
  integer::count_rate,count_max,start,finish,midpoint
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen
  double precision :: trackx,tracky,trackz,intraen,intdummy,totdire,rho,rr
  double precision,dimension(:),allocatable::chen
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy,summing
  real, external :: ran2
  real,dimension(:),allocatable :: variance,disp
  logical :: exist,fail,finalfail,debugyes,film,isbond,clt,clr,noinfo,scalinginfo,nobonds,restart
  !real::start,finish,midpoint,choose1,interesttest,realfinal
  !integer::removethis,pivde
  logical :: suffclust
        real::choose1
  real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
integer,dimension(:),allocatable::unboundarray,boundarray
integer::bondformacccount,bonddeleteacccount,bondformrejcount,bonddeleterejcount
integer::boundcount,unboundcount
!change COM code for large protein in small box to stop pbc effect - similar to old version - but still
  !use dispacement from central bead

  !cluster COM may fall foul of pbcs

  !Sort MSD outputs


  type protein
     integer ::x,y,z,species,type,linker,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl,a111,a11n
  end type protein

  type rprot
     real :: x,y,z
  end type rprot

  type basicp
     integer :: x,y,z
  end type basicp

  type prottemp
     integer :: x,y,z,am,al,bm,bl,cm,cl,dm,dl,em,el,fm,fl
  end type prottemp

  type centremass
     double precision :: x,y,z
  end type centremass

  type lbond
integer:: bm,bl
end type lbond


  type(protein),dimension(:,:),allocatable :: protcoords
  type(centremass),dimension(:,:),allocatable :: com
  type(lbond),dimension(:,:,:),allocatable::cell

     call system_clock(start,count_rate,count_max)

write(6,*) 'NEWCODE'

  mtype = 0
  !noinfo = .false.
  call read_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  binsize = 1
  phasecount = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(summing(mtype,(gridsize*binsize)+1))
  summing(:,:) = 0.0
  totpoints = 0
  do gr = 1,mtype
     totpoints = totpoints + (speciespop(gr)*specieslen(gr))
  end do
  write(6,*) 'total points =',totpoints


  suffclust = .true.
  nprotein = sum(speciespop)
  write(6,*) 'nprotein = ', nprotein
  maxlength = maxval(specieslen)
  write(6,*) 'maxlength =',maxlength
  !open(17, file = 'setup2.txt', action = 'read')
  if(restart .eqv. .false.) then
     open(23, file = 'initialtake2.xyz', action = 'read')
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',action = 'write')
        open(82,file = 'newclusterdata.dat', action = 'write')
     else if(noinfo .eqv. .false.) then
        if(scalinginfo .eqv. .true.) then
           open(11,file='scalingdata.dat',action = 'write')
           !open(29, file = 'rms.dat', action = 'write')
           open(91, file = 'avechainlength.dat', action = 'write')
           open(97, file = 'radiusofgyration.dat', action = 'write')
           open(79, file = 'runningave.dat', action = 'write')
           open(93, file = 'energy.dat', action = 'write')
           open(13,file='timedata.dat',action = 'write')
           open(19,file='acceptance.dat',action  = 'write')
        else if(scalinginfo .eqv. .false.) then
           !open(29, file = 'rms.dat', action = 'write')
           open(93, file = 'energy.dat', action = 'write')
           open(77,file='phasetrans.dat',action = 'write')
           open(82,file = 'newclusterdata.dat', action = 'write')
           open(13,file='timedata.dat',action = 'write')
           open(19,file='acceptance.dat',action  = 'write')
           open(38,file = 'histclust.dat',action = 'write')
           open(3,file='coms.dat',action = 'write')
           !open(67,file='bondddd.dat',action  = 'write')
        end if
     end if
  else if(restart .eqv. .true.) then
     open(23, file = 'checkpoint.xyz', action = 'read')
     write(6,*) 'logicals',scalinginfo,noinfo
     if(noinfo .eqv. .true.) then
        open(11,file='scalingdata.dat',access = 'append')
        open(82,file = 'newclusterdata.dat', access = 'append')
     else if(noinfo .eqv. .false.) then
        if(scalinginfo .eqv. .true.) then
           open(11,file='scalingdata.dat',access = 'append')
           !open(29, file = 'rms.dat', access = 'append')
           open(91, file = 'avechainlength.dat', access = 'append')
           open(97, file = 'radiusofgyration.dat', access = 'append')
           open(79, file = 'runningave.dat', access = 'append')
           open(93, file = 'energy.dat', access = 'append')
           open(13,file='timedata.dat',access = 'append')
           open(19,file='acceptance.dat',access  = 'append')
        else if(scalinginfo .eqv. .false.) then
           !open(29, file = 'rms.dat', access = 'append')
           open(93, file = 'energy.dat', access = 'append')
           open(82,file = 'newclusterdata.dat', access = 'append')
           open(13,file='timedata.dat',access = 'append')
           open(77,file='phasetrans.dat',access = 'append')
           open(19,file='acceptance.dat',access  = 'append')
           open(38,file = 'histclust.dat',access = 'append')
           open(3,file='coms.dat',access = 'append')
        open(886,file='onecorrelation.dat',access='append')
        open(887,file= 'clustercorrelation.dat',access='append')
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

  cltacc = 0
  cltatt = 0
  !kT = 50.0
  successful = 0
  reject = 0
  backtime = maxtime
  reacc = 0
  rerej = 0

  totdire = 0.0
  totiso = 0.0d0
  do ft = 1,mtype
     do gt = 1,mtype
        totdire = totdire + interenergy(ft,gt)
        totiso = totiso + interen(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
  end do
  if(int(totiso) == 0) isbond = .false.
  if(int(totdire) == 0) nobonds =  .true.

  !interenergy = -1.0

  !write(6,*) maxlength1,maxlength2
  write(6,*) 'help',nprotein,maxlength
  write(6,*) 'kT =', kT
  write(6,*) 'clr',clr,'clt',clt

  if(restart .eqv. .false.) then
     !open(67, file = 'move.vtf', action = 'write')
     open(88,file = 'movetagged.vtf',action='write')
  else if(restart .eqv. .true.) then
     !open(67, file = 'move.vtf', access = 'append')
     open(88, file = 'movetagged.vtf', access = 'append')
  end if


  allocate(protcoords(nprotein,maxlength))
  allocate(com(2,nprotein))
  allocate(variance(nprotein))
  allocate(disp(maxtime))
  allocate(chlen(nprotein))
  allocate(chen(nprotein))
  allocate(cell(gridsize,gridsize,gridsize))
allocate(unboundarray(maxlength*nprotein*6))
allocate(boundarray(maxlength*nprotein*6))

  N = maxlength


  count = 0
  time = 0
  totdisp = 0.0d0


  WRITE(6,*) 'ISBOND =',isbond
  write(6,*) 'interen',interen
  write(6,*) 'intraen',intraenergy


  macc = 0
  mrej = 0
  roatt = 0
  roacc = 0


  !sets the limit of the pivot move



  write(6,*) 'foundation start'
  if(restart .eqv. .false.) then
     call foundation
  else if (restart .eqv. .true.) then
     call foundationrestart
  end if
  write(6,*) 'foundation complete'

  if((scalinginfo .eqv. .false.) .and. (restart .eqv. .false.)) call pdbsetupalt
  minl = minval(specieslen)

call initialisearray

  write(6,*) 'a'
  !call dataout
  do qt = 1,nprotein
     call comfind(qt,.false.)
     com(2,qt)%x = com(1,qt)%x
     com(2,qt)%y = com(1,qt)%y
     com(2,qt)%z = com(1,qt)%z
  end do
  write(6,*) 'b'



  if(restart .eqv. .false.) call clustercount !call clustercount2

  actualchainlength = 0.0
  do timebase = 1,maxtime,1
     if(restart .eqv. .false.) then
        time = timebase
     else
        time = timebase+lastruntime
     end if
     if(time < 0) time = abs(time)
     fail = .false.

     !if ((mod(time,outputrate) == 0) .and. (film .eqv. .true.) ) then
     !   call dataout
     !end if
     !write(6,*) 'c'
     do qt = 1,1000
        if((nobonds .eqv. .true.) .and. (clt .eqv. .false.)) then
           choose1 = ran2(seed)-0.5
           if(choose1 <=0.0) call positioning
           if((choose1>0.0) .and. (rept ==1)) call reptation
           if((choose1>0.0) .and. (rept ==0)) call positioning
           !if((modulo(time,maxtime/100) == 0)) then
           !write(29,*) (time-equilib), real(totrmsbrute/nprotein)
           !end if

        else
           call pickmoves
        end if
   
        call system_clock(midpoint,count_rate,count_max)

     if((real(modulo(midpoint-start,count_max))/count_rate)> 216000) then
        backtime = time
        goto 93
     end if

     end do


        if((modulo(time,outputrate/100) == 0) .and. (modulo(time,outputrate) /= 0)) then
        call strippedclustercount
        end if
     if (mod(time,outputrate) == 0) then
        if(noinfo .eqv. .false.)  write(93,*) time,totalenergy
     end if
     if(modulo(time,outputrate) == 0 .and. (scalinginfo .eqv. .false.)) then 
        !call strippedclustercount
        call energy(.false.)
        !call clustercount
        call clustercount2
     end if

     if ((mod(time,outputrate) == 0) .and. (film .eqv. .true.) ) then
        bonddeleterejcount=bonddeleterejcount/1000
        bondformrejcount=bondformrejcount/1000
        bonddeleteacccount=bonddeleteacccount/1000
        bondformacccount = bondformacccount/1000
                call dataout
     end if

!write(6,*) 'e'

     if((time > equilib) .and. (modulo(time,outputrate) == 0) .and. (scalinginfo .eqv. .true.)) then
        call length
        call radiusofgyration
        if((noinfo .eqv. .false.) .and. (scalinginfo .eqv. .true.)) write(79,*) time-equilib,&
             runningaveEtE/((time-equilib)/outputrate),&
             runningaveROG/((time-equilib)/outputrate)
     end if


     if (modulo(time,outputrate/100) == 0 .and. debugyes .eqv. .true.) then
        call debug
        !call energy(.false.)
     end if

     if (fail .eqv. .true.) then
        write(6,*) 'step= ', time, 'FAIL'
     else if(modulo(time,100000) == 0) then
        write(6,*) 'step =',time
     end if


     if((modulo(time,maxtime/10) == 0) .and.(noinfo .eqv. .false.)) then
        write(19,*) 'moves',macc,mrej, real(macc)/(macc+mrej)
        write(19,*) 'reptation',reacc,rerej,real(reacc)/(rerej+reacc)
        write(19,*) 'cluster translate',cltacc,cltatt,real(cltacc)/cltatt
        write(19,*) 'vector moves',roacc,roatt,real(roacc)/roatt   
        write(19,*) 'bond forms',bondformacccount,bondformrejcount,&
                real(bondformacccount)/(bondformacccount+bondformrejcount)
        write(19,*) 'bond deletes',bonddeleteacccount,bonddeleterejcount,&
                real(bonddeleteacccount)/(bonddeleteacccount+bonddeleterejcount) 
    end if
          call system_clock(midpoint,count_rate,count_max)

     if((real(modulo(midpoint-start,count_max))/count_rate)> 216000) then
        backtime = time
        goto 93
     end if
!write(6,*) 'sweep complete'
  end do


  if(restart .eqv. .false.) then
     backtime = timebase
  else
     backtime = timebase+lastruntime
  end if


93 if((film .eqv. .false.) .and. (scalinginfo .eqv. .false.)) call dataout
  write(6,*) 'deltax =',trackx, 'deltay =',tracky, 'deltaz =',trackz


     call system_clock(finish,count_rate,count_max)

  call dataoutrestart
  if(noinfo .eqv. .false.) write(13,*)real(modulo(finish-start,count_max))/count_rate
 
 !if(noinfo .eqv. .false.) write(13,*) (finish-start)
  !write(6,*) 'TIME to FINISH',finish-start  

  !call error
  call energy(.false.)
  if((noinfo .eqv. .true.) .or. (scalinginfo .eqv. .true.)) write(11,*) maxlength1, &
       (polymerl/((backtime-equilib)/outputrate)), &
       (polymerrog/((backtime-equilib)/outputrate))
  write(6,*) 2*log(actualchainlength/backtime)/log(real(maxlength))
  write(6,*) 'average clusper pop', (real(avecluspop)/aveclusnumber)/((backtime-equilib)/outputrate), &
       real(avecluspop)/((backtime-equilib)/outputrate), real(aveclusnumber)/((backtime-equilib)/outputrate)



  if(restart .eqv. .false.) then
     continue
  else if (restart .eqv. .true.) then

     timeofinterest = phasecount
     !(backtime - lastruntime)/outputrate
     write(6,*) 'time of interest',timeofinterest,backtime,lastruntime
     open(113,file='comrdf.dat',action = 'write')
     open(114,file='unnormcomrdf.dat',action  = 'write')

     rho = (totpoints)/real(gridsize**3)
     write(113,*) 0,0,0
     do gr = 1,(binsize*(gridsize/2)) + 1
        rr = (real(gr)/binsize)
        !write(6,*) 'rho',rho
        sumr1= 0
        sumr2 = 0
        do yl = 1,gr,1
           sumr1 = sumr1 + summing(1,yl)
           sumr2 = sumr2 + summing(2,yl)
        end do

        write(113,*) rr,sumr1/(timeofinterest*(4.0/3.0)*PI_8*(rr**3)),sumr2/(timeofinterest*(4.0/3.0)*PI_8*(rr**3)),&
             (sumr2+sumr1)/(timeofinterest*(4.0/3.0)*PI_8*(rr**3))
        write(114,*) rr,sumr1/timeofinterest,sumr2/timeofinterest,summing(1,gr)/timeofinterest,summing(2,gr)/timeofinterest
     end do
     close(113)
     close(114)
  end if



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

  !call CPU_TIME(realfinal)
  !write(6,*) 'Very final time:',realfinal-finish
  !if(noinfo .eqv. .false.) write(13,*) (realfinal-start)

contains



subroutine unboundcounter
integer::m,l,maxl,g,f

unboundcount=0
boundcount=0



    do m= 1,nprotein,1
   maxl = chlen(m)
        do l=1,maxl,1
                if((protcoords(m,l)%am /=0) .and. (protcoords(m,l)%am /=m))then
        if(protcoords(m,l)%am>m) then
                        if(int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%am,1)%type))/=0) then
      g=protcoords(m,l)%am
      f=protcoords(m,l)%al
      if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+1
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+2
                        else if(protcoords(m,l)%a111==1)then
                                boundcount = boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+1
                                boundcount = boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+2
                        end if
                        end if
                end if
        end if

                if((protcoords(m,l)%bm /=0).and. (protcoords(m,l)%bm /=m)) then

        if(protcoords(m,l)%bm>m) then
                        if(int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%bm,1)%type))/=0) then
        g=protcoords(m,l)%bm
        f=protcoords(m,l)%bl
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+2
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+1
                        else if(protcoords(m,l)%a111==2)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+2
                                boundcount =boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+1
                        end if
                end if
        end if
end if

                if((protcoords(m,l)%cm /=0) .and. (protcoords(m,l)%cm /=m))  then
        if(protcoords(m,l)%cm>m) then

                        if (int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%cm,1)%type))/=0) then
        g=protcoords(m,l)%cm
        f=protcoords(m,l)%cl
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+3
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+4
                        else if(protcoords(m,l)%a111==3)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+3
                                boundcount =boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+4
                        end if
                end if
        end if
end if


                if((protcoords(m,l)%dm /=0) .and. (protcoords(m,l)%dm /=m))then 
        if(protcoords(m,l)%dm>m) then

                        if(int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%dm,1)%type))/=0) then
        g=protcoords(m,l)%dm
        f=protcoords(m,l)%dl
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+4
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+3
                        else if(protcoords(m,l)%a111==4)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+4
                                boundcount =boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+3
                        end if
                end if
        end if
end if

                if((protcoords(m,l)%em /=0).and. (protcoords(m,l)%em /=m)) then 
        if(protcoords(m,l)%em>m) then

                        if(int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%em,1)%type))/=0) then
        g=protcoords(m,l)%em
        f=protcoords(m,l)%el
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+5
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+6
                        else if(protcoords(m,l)%a111==5)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+5
                                boundcount =boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+6
                        end if
                end if
        end if
end if


                if((protcoords(m,l)%fm /=0).and. (protcoords(m,l)%fm /=m)) then
        if(protcoords(m,l)%fm>m) then

                        if(int(interenergy(protcoords(m,l)%type,protcoords(protcoords(m,l)%fm,1)%type))/=0) then
        g=protcoords(m,l)%fm
        f=protcoords(m,l)%fl
                        if(protcoords(m,l)%a111==0) then
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+6
                                unboundcount = unboundcount+1
                                unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+5
                        else if(protcoords(m,l)%a111==6)then
                                boundcount =boundcount+1
                                boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+6
                                boundcount =boundcount+1
                                boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+5
                        end if
                end if
end if
end if

        end do
    end do



end subroutine unboundcounter



subroutine bondformatt
integer::pick,locator,m,l,dir,g,f,delbond,dir2,newlocator,maxl
double precision:: deltaenergy,denergy
logical::formlogic
real::realpick
deltaenergy = 0.0d0



if(unboundcount==0) return


formlogic=.false.
pick = int(ran2(seed)*unboundcount)+1
realpick=real(pick)
if(int(2*((realpick/2)-floor(realpick/2)))==0) pick = pick-1

locator=unboundarray(pick)
m=floor(real(locator-1)/(6*maxlength))+1
maxl=chlen(m)
newlocator=locator-(6*maxl*(m-1))
l=floor(real(newlocator-1)/6)+1


if(protcoords(m,l)%a111>0) then
bondformrejcount= bondformrejcount+1
return
end if

dir=newlocator-(6*(l-1))

SELECT CASE (dir)

CASE (1)
g=protcoords(m,l)%am
f=protcoords(m,l)%al
dir2=2
CASE (2)
g=protcoords(m,l)%bm
f=protcoords(m,l)%bl
dir2=1
CASE (3)
g=protcoords(m,l)%cm
f=protcoords(m,l)%cl
dir2=4
CASE (4)
g=protcoords(m,l)%dm
f=protcoords(m,l)%dl
dir2=3
CASE(5)
g=protcoords(m,l)%em
f=protcoords(m,l)%el
dir2=6
CASE(6)
g=protcoords(m,l)%fm
f=protcoords(m,l)%fl
dir2=5
END SELECT

if(g==m) then
bondformrejcount= bondformrejcount+1
return
end if

if(protcoords(g,f)%a111>0) then
bondformrejcount= bondformrejcount+1
return
end if

!if(interenergy(protcoords(m,l)%type,protcoords(g,f)%type) >=0) then
if(protcoords(m,l)%type==protcoords(g,f)%type) then
bondformrejcount= bondformrejcount+1
return
end if



    deltaenergy = interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    denergy=(real(unboundcount)/real(boundcount+1))

    !if(deltaenergy <= 0.0) then
    !    formlogic=.true.
    !end if
    if((denergy*exp(-deltaenergy/kT)) > ran2(seed)) then
        formlogic=.true.
    else
        bondformrejcount=bondformrejcount+1
    return
    end if



if (formlogic .eqv. .true.) then




    boundcount=boundcount+1
    unboundcount=unboundcount-1
    protcoords(m,l)%a111=dir
    protcoords(m,l)%a11n=g
    protcoords(g,f)%a111=dir2
    protcoords(g,f)%a11n=m
        maxl=chlen(m)
    boundarray(boundcount)=((m-1)*maxlength*6)+((l-1)*6)+dir
    boundcount=boundcount+1
    unboundcount=unboundcount-1
    boundarray(boundcount)=((g-1)*maxlength*6)+((f-1)*6)+dir2


    totalenergy=totalenergy+deltaenergy
    do delbond=pick,unboundcount,1
        unboundarray(delbond)=unboundarray(delbond+2)
    end do
        unboundarray(unboundcount+1)=0
        unboundarray(unboundcount+2)=0
        bondformacccount=bondformacccount+1

end if



end subroutine bondformatt


subroutine bonddeleteatt
integer::pick,locator,m,l,dir,g,f,delbond,dir2,newlocator,maxl
logical:: formlogic
double precision:: deltaenergy,denergy
real::realpick
if(boundcount==0) return


deltaenergy = 0.0d0
formlogic = .false.
pick = int(ran2(seed)*boundcount)+1
realpick=real(pick)
if(int(2*((realpick/2)-floor(realpick/2)))==0) pick = pick-1


locator=boundarray(pick)
m=floor(real(locator-1)/(6*maxlength))+1
maxl=chlen(m)
newlocator=locator-(6*maxl*(m-1))
l=floor(real(newlocator-1)/6)+1
if(protcoords(m,l)%a111==0) then
bonddeleterejcount= bonddeleterejcount+1
return
end if


dir=newlocator-(6*(l-1))
!if(dir>6) write(6,*) 'DIR FAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


SELECT CASE (dir)

CASE (1)
g=protcoords(m,l)%am
f=protcoords(m,l)%al
dir2=2
CASE (2)
g=protcoords(m,l)%bm
f=protcoords(m,l)%bl
dir2=1
CASE (3)
g=protcoords(m,l)%cm
f=protcoords(m,l)%cl
dir2=4
CASE (4)
g=protcoords(m,l)%dm
f=protcoords(m,l)%dl
dir2=3
CASE(5)
g=protcoords(m,l)%em
f=protcoords(m,l)%el
dir2=6
CASE(6)
g=protcoords(m,l)%fm
f=protcoords(m,l)%fl
dir2=5
END SELECT


if(protcoords(g,f)%a111/=dir2) then
        write(6,*) 'FAILLLLLL - bond mismatch!' 
        write(6,*) boundarray(pick)
        write(6,*) m,l,g,f,dir,dir2,protcoords(m,l)%a111
        write(6,*) protcoords(m,l)
        write(6,*) protcoords(g,f)
        write(6,*) boundarray(1:100)
        stop 99
end if
if(g==m) then
bonddeleterejcount= bonddeleterejcount+1
return
end if



if(protcoords(g,f)%a111==0) then
bonddeleterejcount= bonddeleterejcount+1
return
end if

!if(interenergy(protcoords(m,l)%type,protcoords(g,f)%type) >=0) then
if(protcoords(m,l)%type==protcoords(g,f)%type)  then
bonddeleterejcount= bonddeleterejcount+1
return
end if



    deltaenergy = -1*interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
    denergy=(real(boundcount)/real(unboundcount+1))

    !if(deltaenergy <= 0.0) then
    !    formlogic=.true.
    !end if
    if((denergy*exp(-deltaenergy/kT)) > ran2(seed)) then
        formlogic=.true.
    else
        bonddeleterejcount=bonddeleterejcount+1
    return
    end if



        if(formlogic .eqv. .true.) then
    unboundcount=unboundcount+1
    boundcount=boundcount-1
    protcoords(m,l)%a111=0
    protcoords(m,l)%a11n=g
    protcoords(g,f)%a111=0
    protcoords(g,f)%a11n=m
        maxl=chlen(m)
    unboundarray(unboundcount)=((m-1)*maxlength*6)+((l-1)*6)+dir
    unboundcount=unboundcount+1
    boundcount=boundcount-1
    unboundarray(unboundcount)=((g-1)*maxlength*6)+((f-1)*6)+dir2

    !if((m==339).and.(l==3)) write(6,*) 'bondremove',unboundarray(unboundcount)
    !if((m==291).and.(l==3)) write(6,*) 'bondremove',unboundarray(unboundcount)

totalenergy=totalenergy+deltaenergy
    do delbond=pick,boundcount,1
        boundarray(delbond)=boundarray(delbond+2)
    end do
        boundarray(boundcount+1)=0
        boundarray(boundcount+2)=0
        bonddeleteacccount= bonddeleteacccount+1
        !if((m==79).or.(g==79)) write(6,*) 'unboundarray',unboundarray(unboundcount),m,l,dir 

     
end if




end subroutine bonddeleteatt

  subroutine entropy(deltaenergy,tempcoord,m,l)
    double precision,intent(inout) ::deltaenergy
    integer,intent(in)::m,l
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    double precision::oldent,newent,sep1,sep2,dx,dy,dz
    integer::maxl

    maxl = chlen(m)
    oldent = 0.0d0
    newent = 0.0d0
    if(l>1) then
          dx = min(abs(protcoords(m,l)%x - protcoords(m,l-1)%x), gridsize -abs(protcoords(m,l)%x - protcoords(m,l-1)%x))
          dy = min(abs(protcoords(m,l)%y - protcoords(m,l-1)%y), gridsize -abs(protcoords(m,l)%y - protcoords(m,l-1)%y))
          dz = min(abs(protcoords(m,l)%z - protcoords(m,l-1)%z), gridsize -  abs(protcoords(m,l)%z - protcoords(m,l-1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))
          
          dx = min(abs(tempcoord(l)%x - protcoords(m,l-1)%x), gridsize -abs(tempcoord(l)%x - protcoords(m,l-1)%x))
          dy = min(abs(tempcoord(l)%y - protcoords(m,l-1)%y), gridsize -abs(tempcoord(l)%y - protcoords(m,l-1)%y))
          dz = min(abs(tempcoord(l)%z - protcoords(m,l-1)%z), gridsize -  abs(tempcoord(l)%z - protcoords(m,l-1)%z))

          sep2 = sqrt((dx**2)+(dy**2)+(dz**2))

    
    oldent = oldent + entropycalc(sep1,m)
    newent = newent +entropycalc(sep2,m)
    end if

    if(l<maxl) then
          dx = min(abs(protcoords(m,l)%x - protcoords(m,l+1)%x), gridsize -abs(protcoords(m,l)%x - protcoords(m,l+1)%x))
          dy = min(abs(protcoords(m,l)%y - protcoords(m,l+1)%y), gridsize -abs(protcoords(m,l)%y - protcoords(m,l+1)%y))
          dz = min(abs(protcoords(m,l)%z - protcoords(m,l+1)%z), gridsize -  abs(protcoords(m,l)%z - protcoords(m,l+1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))
          
          dx = min(abs(tempcoord(l)%x - protcoords(m,l+1)%x), gridsize -abs(tempcoord(l)%x - protcoords(m,l+1)%x))
          dy = min(abs(tempcoord(l)%y - protcoords(m,l+1)%y), gridsize -abs(tempcoord(l)%y - protcoords(m,l+1)%y))
          dz = min(abs(tempcoord(l)%z - protcoords(m,l+1)%z), gridsize -  abs(tempcoord(l)%z - protcoords(m,l+1)%z))

          sep2 = sqrt((dx**2)+(dy**2)+(dz**2))

    
    oldent = oldent + entropycalc(sep1,m)
    newent = newent +entropycalc(sep2,m)
    end if

    deltaenergy = deltaenergy+(newent-oldent)
    

  end subroutine entropy

 double precision Function entropycalc(separ,m)
    double precision, intent(in) :: separ
        integer,intent(in)::m
double precision :: Na
Na = protcoords(m,1)%linker
 entropycalc = -1*((3*(separ**2)*kT)/(2*(Na**2)))
  end Function entropycalc


    subroutine entropyrep(deltaenergy,tempcoord,m,l)
    double precision,intent(inout) ::deltaenergy
    integer,intent(in)::m,l
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    double precision::oldent,newent,sep1,sep2,dx,dy,dz
    integer::maxl,f!,h

    maxl = chlen(m)
    if(l == 1) then
       !h = 2
       f = maxl
    else if(l==maxl) then
       f = 2
       !h = maxl-1
    end if
    oldent = 0.0d0
    newent = 0.0d0

          dx = min(abs(protcoords(m,f)%x - protcoords(m,f-1)%x), gridsize -abs(protcoords(m,f)%x - protcoords(m,f-1)%x))
          dy = min(abs(protcoords(m,f)%y - protcoords(m,f-1)%y), gridsize -abs(protcoords(m,f)%y - protcoords(m,f-1)%y))
          dz = min(abs(protcoords(m,f)%z - protcoords(m,f-1)%z), gridsize -  abs(protcoords(m,f)%z - protcoords(m,f-1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))
          
          dx = min(abs(tempcoord(l)%x - protcoords(m,l)%x), gridsize -abs(tempcoord(l)%x - protcoords(m,l)%x))
          dy = min(abs(tempcoord(l)%y - protcoords(m,l)%y), gridsize -abs(tempcoord(l)%y - protcoords(m,l)%y))
          dz = min(abs(tempcoord(l)%z - protcoords(m,l)%z), gridsize -  abs(tempcoord(l)%z - protcoords(m,l)%z))

          sep2 = sqrt((dx**2)+(dy**2)+(dz**2))

    
    oldent = entropycalc(sep1,m)
    newent = entropycalc(sep2,m)


    deltaenergy = deltaenergy+(newent-oldent)
    

  end subroutine entropyrep


    subroutine entropyinit(deltaenergy,m,l)
    double precision,intent(inout) ::deltaenergy
    integer,intent(in)::m,l
    double precision::sep1,dx,dy,dz
    integer::maxl

    maxl = chlen(m)
    if(l ==maxl) return
    
          dx = min(abs(protcoords(m,l)%x - protcoords(m,l+1)%x), gridsize -abs(protcoords(m,l)%x - protcoords(m,l+1)%x))
          dy = min(abs(protcoords(m,l)%y - protcoords(m,l+1)%y), gridsize -abs(protcoords(m,l)%y - protcoords(m,l+1)%y))
          dz = min(abs(protcoords(m,l)%z - protcoords(m,l+1)%z), gridsize -  abs(protcoords(m,l)%z - protcoords(m,l+1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))


    deltaenergy = deltaenergy + entropycalc(sep1,m)
    

  end subroutine entropyinit
  
 subroutine entropyinit1(deltaenergy,m,l)
    double precision,intent(inout) ::deltaenergy
    integer,intent(in)::m,l
    double precision::sep1,dx,dy,dz
    integer::maxl

    maxl = chlen(m)
    if(l ==maxl) return

          dx = min(abs(protcoords(m,l)%x - protcoords(m,l+1)%x), gridsize-abs(protcoords(m,l)%x - protcoords(m,l+1)%x))
          dy = min(abs(protcoords(m,l)%y - protcoords(m,l+1)%y), gridsize-abs(protcoords(m,l)%y - protcoords(m,l+1)%y))
          dz = min(abs(protcoords(m,l)%z - protcoords(m,l+1)%z), gridsize -abs(protcoords(m,l)%z - protcoords(m,l+1)%z))

          sep1 = sqrt((dx**2)+(dy**2)+(dz**2))


    deltaenergy = deltaenergy + (2*entropycalc(sep1,m))


  end subroutine entropyinit1


subroutine initialisearray
integer::m,l,maxl,axx,bxx,cxx



do axx=1,gridsize,1
do bxx=1,gridsize,1
do cxx=1,gridsize,1

cell(axx,bxx,cxx)%bm=0
cell(axx,bxx,cxx)%bl=0
end do
end do
end do

do m = 1,nprotein,1
maxl = chlen(m)
do l=1,maxl

        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bm = m
        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl = l


end do
end do



end subroutine initialisearray


subroutine updatearray(m,l,tempcoord)
integer,intent(in)::m,l
    type(prottemp),dimension(:),intent(in) :: tempcoord

        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bm = 0
        cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl = 0
        cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bm = m
        cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bl = l



end subroutine updatearray


subroutine updatearrayreptation(m,l,f,tempcoord)
integer,intent(in)::m,l,f
    type(prottemp),dimension(:),intent(in) :: tempcoord

        cell(protcoords(m,f)%x,protcoords(m,f)%y,protcoords(m,f)%z)%bm = 0
        cell(protcoords(m,f)%x,protcoords(m,f)%y,protcoords(m,f)%z)%bl = 0
        cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bm = m
        cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bl = l



end subroutine updatearrayreptation


  subroutine phase(mz,clsize,clnos,cb)
    !subroutine to find the radius of gyration of the largest cluster      !
    !also contains subroutine clustercount for finding the largest cluster !
    
    integer :: m,l,maxl,h,a,totcluspop,rogpop,maxrogpop,g
    double precision :: radofgy,maxradofgy
    type(centremass)::avepos,sumsq
    integer,intent(inout)::clsize,mz
    integer,dimension(:),allocatable::cllist
    integer,dimension(:,:),allocatable:: route
    integer,dimension(:),allocatable::rcounter,checklist
    integer,dimension(:,:),intent(inout)::cb
    integer,dimension(:),intent(inout)::clnos
    type(basicp) :: comcluster
    type(rprot) :: dpcomcluster
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
    call clustercom(mz,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge,dpcomcluster)
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

    if (restart .eqv. .true.) call rdf(dpcomcluster)

    maxradofgy = 0.0d0
    do g = 1,clsize
       m = cllist(g)
       maxl = chlen(m)
       do l = 1,maxl
          sumsq%x = min(abs(protcoords(m,l)%x-comcluster%x),gridsize-abs(protcoords(m,l)%x-comcluster%x))
          sumsq%y = min(abs(protcoords(m,l)%y-comcluster%y),gridsize-abs(protcoords(m,l)%y-comcluster%y))
          sumsq%z = min(abs(protcoords(m,l)%z-comcluster%z),gridsize-abs(protcoords(m,l)%z-comcluster%z))
          maxradofgy = maxradofgy + (sumsq%x**2) + (sumsq%y**2) + (sumsq%z**2)
       end do
       maxrogpop = maxrogpop + maxl
    end do
    maxradofgy = maxradofgy/maxrogpop


    write(77,*) time,radofgy,sqrt(radofgy),maxradofgy,sqrt(maxradofgy)


  end subroutine phase

  subroutine rdf(dpcomcluster)
    !subroutine to calculate the RDF of the beads from the centre of mass of the largest cluster  !
    type(rprot),intent(in) :: dpcomcluster
    integer :: m,l,maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta


    phasecount = phasecount+1
    do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          dx = min(abs(protcoords(m,l)%x - dpcomcluster%x), gridsize - abs(protcoords(m,l)%x - dpcomcluster%x))
          dy = min(abs(protcoords(m,l)%y - dpcomcluster%y), gridsize - abs(protcoords(m,l)%y - dpcomcluster%y))
          dz = min(abs(protcoords(m,l)%z - dpcomcluster%z), gridsize - abs(protcoords(m,l)%z - dpcomcluster%z))

          delta = sqrt((dx**2)+(dy**2)+(dz**2))

          lx = INT(delta/(1.0/binsize))


          summing(protcoords(m,l)%type,lx+1) = summing(protcoords(m,l)%type,lx+1)+1.0
       end do
    end do

  end subroutine rdf


  subroutine updatebondold(m,beadmin,beadmax,tempcoord)
    !Subroutine to give protcoords initial neighbourlist after           !
    !being generated in initial energy subroutine                        !
    integer,intent(in):: m,beadmin,beadmax
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f,bg,bf

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

    end do

  end subroutine updatebondold



  subroutine foundation
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z
    character(len = 10) :: BIN
    read(23,*) BIN
    !write(6,*) maxlength

    !write(67,*) ' '
    do m = 1,nprotein,1
       read(23,*) protcoords(m,1)%species,protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z,protcoords(m,1)%type,chlen(m)&
            ,protcoords(m,1)%linker

       maxl = chlen(m)
protcoords(m,1)%a111=0
protcoords(m,1)%a11n=0


       do l = 2,maxl,1
          !write(6,*) l
          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z,protcoords(m,l)%type
protcoords(m,l)%a111=0
protcoords(m,l)%a11n=0

       end do


    end do

    call energy(.TRUE.)
    call debug
    if(scalinginfo .eqv. .true.) call pdbsetupalt
  
        
        unboundcount = 0
boundcount = 0

        end subroutine foundation




  subroutine foundationrestart
    !read in initial positions of the protein
    integer::m,l,no1,no2,maxl,z
    character(len = 10) :: BIN
    read(23,*)lastruntime,totrmsbrute,suffclust,unboundcount,boundcount,totalenergy
    lastruntime = abs(lastruntime)
    if(lastruntime >= 1000000000) lastruntime = lastruntime/1000
    do m = 1,nprotein,1
       read(23,*) protcoords(m,1)%species, protcoords(m,1)%x, protcoords(m,1)%y, protcoords(m,1)%z,&
            protcoords(m,1)%type,chlen(m),protcoords(m,1)%linker,protcoords(m,1)%a111,protcoords(m,1)%a11n

       maxl = chlen(m)

       do l = 2,maxl,1
          !write(6,*) l
          read(23,*) protcoords(m,l)%species, protcoords(m,l)%x, protcoords(m,l)%y, protcoords(m,l)%z,&
               protcoords(m,l)%type,protcoords(m,l)%a111,protcoords(m,l)%a11n
       end do
    end do



    call energy(.TRUE.)
    call debug

        !unboundcount = 0
!boundcount = 0
call unboundcounter
    close(23)
  end subroutine foundationrestart

  subroutine clusterposition
    !subroutine to select the cluster move to be applied  !
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


    if((clt .eqv. .false.)) then
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

  end subroutine clusterposition



  subroutine rotate(m,l,b,tempcoords,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz)
    !redundant subroutine for rotation move !
    integer,intent(in):: m,l,b,xx,xy,xz,yx,yy,yz,zx,zy,zz,dx,dy,dz
    type(prottemp),dimension(:),intent(inout) :: tempcoords
    !write(6,*) dx,dy,dz
    tempcoords(b)%x = modulo(protcoords(m,l)%x + (xx*dx) + (xy*dy) + (xz*dz) -1,gridsize)+1
    tempcoords(b)%y = modulo(protcoords(m,l)%y + (yx*dx) + (yy*dy) + (yz*dz)-1,gridsize)+1
    tempcoords(b)%z = modulo(protcoords(m,l)%z + (zx*dx) + (zy*dy) + (zz*dz)-1,gridsize)+1

  end subroutine rotate

  subroutine pickmoves
    !Subroutine to select which moves to attemp   !
    integer:: scans,m,maxl,decide,l,a
        real::bondformdelete    
!double precision :: decide
    scans = 1
    !do scans = 1,(nprotein1*maxlength1 + nprotein2*maxlength2)
    decide = int(ran2(seed)*10)+1
    !if(time > 570000 .and. (time < 573000))write(6,*) 'decide',decide,time

    !if(protein(m,1)%species==3) then
    !call clusterposition
    !return
    !end if
    !write(6,*) 'start loop',decide
!call energy(.false.)
    if(decide <= 4) then
       do a=1,maxlength
       bondformdelete = ran2(seed)-0.5
        if(bondformdelete>0.0) then
        call bondformatt
        else if(bondformdelete<0.0) then
        call bonddeleteatt
        end if

       end do
    end if
    if((decide > 4) .and. (decide <= 6) )then  
           call reptation
    else if((decide > 6) .and. (decide <= 8)) then
             call positioning
    else if((decide > 8) ) then
            call clusterposition
     end if      

!write(6,*) 'decide',decide
!call energy(.false.)
!    write(6,*) 'decide',decide
!call debug



    !if((modulo(time,maxtime/100) == 0)) then
    !write(29,*) (time-equilib)*scans, real(totrmsbrute)/nprotein
    !end if

  end subroutine pickmoves


  logical Function Energydecision(denergy)
    !function to select whether to accept or reject a move !
    double precision, intent(in) :: denergy
    if(denergy <= 0.0) then
       Energydecision = .True.
    else if(exp(-denergy/kT) > ran2(seed)) then
       Energydecision = .True.
    else
       Energydecision = .False.
    end if
  end Function Energydecision



  subroutine positioning
    !displaces a bead by 0 or +/- 1 in x,y or z direction
    integer :: m,l,deltax1,deltax2,deltay1,deltay2,deltaz1,deltaz2,dx,dy,dz,bdir,str
    double precision :: deltaenergy
    logical :: rac,racc,adjver
    integer:: st,pr,maxl,maxback,decide
    !double precision::sumdebug
    type(prottemp),dimension(:),allocatable :: tempcoord
    !integer::fakex,fakey,fakez
    allocate(tempcoord(maxlength))
    !return    
    rac = .true.
    m = int(ran2(seed)*nprotein)+1

    !if(protcoords(m,1)%species == 3) then
    !call clusterposition
    !return 
    !end if

    maxl = chlen(m)
    l = int(ran2(seed)*maxl)+1

if(protcoords(m,l)%a111>0) then
       rac = .false.
        goto 37
    end if

    !if((m==79).and. (l==3)) write(6,*) protcoords(m,l)%a111 
    !if((m==428) .and. (l==2)) write(6,*) 'a111 prepositioning' ,protcoords(m,l)%a111
    
    dx = 0
    dy = 0
    dz = 0
    !write(6,*) 'positioning'


    decide = int(ran2(seed)*3)+1

    SELECT CASE (decide)

    CASE (1)
       dx = 1
    CASE (2)
       dx = 0
    CASE (3)
       dx = -1

    END SELECT

    decide = int(ran2(seed)*3)+1

    SELECT CASE (decide)

    CASE (1)
       dy = 1
    CASE (2)
       dy = 0
    CASE (3)
       dy = -1

    END SELECT

    decide = int(ran2(seed)*3)+1

    SELECT CASE(decide)

    CASE (1)
       dz = 1
    CASE (2)
       dz = 0
    CASE (3)
       dz = -1

    END SELECT


    tempcoord(l)%x = modulo(protcoords(m,l)%x+dx-1,gridsize)+1
    tempcoord(l)%y = modulo(protcoords(m,l)%y+dy-1,gridsize)+1
    tempcoord(l)%z = modulo(protcoords(m,l)%z+dz-1,gridsize)+1


    
    if(l>1) then
       dx = min(abs(tempcoord(l)%x - protcoords(m,l-1)%x),gridsize-abs(tempcoord(l)%x &
            - protcoords(m,l-1)%x)) 
       dy = min(abs(tempcoord(l)%y - protcoords(m,l-1)%y),gridsize- abs(tempcoord(l)%y &
            - protcoords(m,l-1)%y)) 
       dz = min(abs(tempcoord(l)%z - protcoords(m,l-1)%z),gridsize - abs(tempcoord(l)%z &
            - protcoords(m,l-1)%z)) 

       sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
       if(sumdebug > protcoords(m,1)%linker) then
          rac = .false.
          goto 37
       end if
    end if

    if(l<maxl) then

       dx = min(abs(tempcoord(l)%x - protcoords(m,l+1)%x),gridsize-abs(tempcoord(l)%x &
            - protcoords(m,l+1)%x)) 
       dy = min(abs(tempcoord(l)%y - protcoords(m,l+1)%y),gridsize- abs(tempcoord(l)%y &
            - protcoords(m,l+1)%y)) 
       dz = min(abs(tempcoord(l)%z - protcoords(m,l+1)%z),gridsize - abs(tempcoord(l)%z &
            - protcoords(m,l+1)%z)) 

       sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
       if(sumdebug > protcoords(m,1)%linker) then
          rac = .false.
          goto 37
       end if
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


if(cell(tempcoord(l)%x,tempcoord(l)%y,tempcoord(l)%z)%bm /=0) then
rac = .false.
goto 37
end if


    tempcoord(l)%am = 0
    tempcoord(l)%bm = 0
    tempcoord(l)%cm = 0
    tempcoord(l)%dm = 0
    tempcoord(l)%em = 0
    tempcoord(l)%fm = 0
    tempcoord(l)%al = 0
    tempcoord(l)%bl = 0
    tempcoord(l)%cl = 0
    tempcoord(l)%dl = 0
    tempcoord(l)%el = 0
    tempcoord(l)%fl = 0



    call removeenergyintra(m,l,deltaenergy,tempcoord)
    call removeenergyinter(m,l,deltaenergy,tempcoord)

    do str = 1,l-3,1
       call adjacent(m,l,str,tempcoord,adjver,bdir)
       if((adjver.eqv. .true.)) then
          call bondformintra(m,l,str,tempcoord,bdir,deltaenergy)
       end if
    end do
    do str = l+3,maxl,1
       call adjacent(m,l,str,tempcoord,adjver,bdir)
       if((adjver.eqv. .true.)) then
          call bondformintra(m,l,str,tempcoord,bdir,deltaenergy)
       end if
    end do



if(cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0)then
        if(cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=m)then
        
        pr=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

if(cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=0)then
     if(cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm/=m)then
 
          pr=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bm
                str=cell(modulo(tempcoord(l)%x-2,gridsize)+1,tempcoord(l)%y,tempcoord(l)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if



        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm/=0)then
        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm/=m)then

          pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y,gridsize)+1,tempcoord(l)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm/=0)then
        if(cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm/=m)then
                pr=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bm
                str=cell(tempcoord(l)%x,modulo(tempcoord(l)%y-2,gridsize)+1,tempcoord(l)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z,gridsize)+1)%bm/=0)then
        if(cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z,gridsize)+1)%bm/=m)then
      
          pr=cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z,gridsize)+1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z,gridsize)+1)%bl
        bdir = -3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if


        if(cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z-2,gridsize)+1)%bm/=0)then      
        if(cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z-2,gridsize)+1)%bm/=m)then
          pr=cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z-2,gridsize)+1)%bm
                str=cell(tempcoord(l)%x,tempcoord(l)%y,modulo(tempcoord(l)%z-2,gridsize)+1)%bl
        bdir = 3
                      call bondforminter(m,pr,l,str,tempcoord,bdir,deltaenergy)
end if
end if

call entropy(deltaenergy,tempcoord,m,l)

    if(Energydecision(deltaenergy) .eqv. .true.) then
       totalenergy = totalenergy + deltaenergy
        call updatearray(m,l,tempcoord)
         call updatepos(m,l,l,tempcoord,maxl)
       call updatebond(m,l,l,tempcoord)
         macc = macc + 1
call unboundcounter
    else
       rac = .false.
       goto 37
    end if


37  if (rac .eqv. .false.) then
       mrej = mrej + 1
    end if




    end subroutine positioning



  subroutine updatepos(chainnum,beadmin,beadmax,tempcoord,maxl)
    !updates bead positions
    integer,intent(in):: chainnum,beadmin,beadmax,maxl
    Type(prottemp),dimension(:),intent(in) :: tempcoord
    integer ::beadnum
    !moves beads to new positions and reassigns isobonding


    do beadnum = beadmin,beadmax,1
       protcoords(chainnum,beadnum)%x = tempcoord(beadnum)%x
       protcoords(chainnum,beadnum)%y = tempcoord(beadnum)%y
       protcoords(chainnum,beadnum)%z = tempcoord(beadnum)%z
    end do

    call comfind(chainnum,.TRUE.)


    !call rms(chainnum)
  end subroutine updatepos

  subroutine updatecluspos(tempcluscoord,tcom,cllist,clsize)
    !updates bead positions
    integer,dimension(:),intent(in) :: cllist
    integer,intent(in) :: clsize
    Type(prottemp),dimension(:,:),intent(in) :: tempcluscoord
    type(centremass),dimension(:),intent(in)::tcom
    integer ::b,l,a,maxback,ra
    !moves beads to new positions and reassigns isobonding

    do ra = 1,clsize
       a = cllist(ra)
       maxback = chlen(a)
       !call rms(a)
       do l = 1,maxback,1
          protcoords(a,l)%x = tempcluscoord(a,l)%x
          protcoords(a,l)%y = tempcluscoord(a,l)%y
          protcoords(a,l)%z = tempcluscoord(a,l)%z

       end do
       com(2,a)%x = com(1,a)%x
       com(2,a)%y = com(1,a)%y
       com(2,a)%z = com(1,a)%z
       com(1,a)%x = tcom(a)%x
       com(1,a)%y = tcom(a)%y
       com(1,a)%z = tcom(a)%z
       !call comfind(a,.TRUE.)
       !call rms(a)
       call tracking(a)
    end do

    !write(6,*) 'successful cluster move'
  end subroutine updatecluspos


  subroutine clusterassignmove(chainnum,cllist,clsize,cl)
    !Subroutine for assigning the chains connected to the    !
    ! chain of interest into a  clusters                     !  
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
    !subroutine used by clusterassignmove to check if a bead is bound to its nearest neighbour  !
    ! i.e. the bond vectors point towards each other                                            !

    integer,intent(in) :: bdir,m,l
    integer,intent(inout)::a,cluspop
    logical,dimension(:),intent(inout)::cl
    integer,dimension(:),intent(inout):: cllist
    integer :: f,g,bdir2,bdir3



SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
                g = protcoords(m,l)%am
                f = protcoords(m,l)%al

CASE (1)
bdir2=2
bdir3=1
                g = protcoords(m,l)%bm
                f = protcoords(m,l)%bl

CASE (-2)
bdir2=3
bdir3=4
                g = protcoords(m,l)%cm
                f = protcoords(m,l)%cl

CASE (2)
bdir2=4
bdir3=3
                g = protcoords(m,l)%dm
                f = protcoords(m,l)%dl

CASE(-3)
bdir2=5
bdir3=6
                g = protcoords(m,l)%em
                f = protcoords(m,l)%el

CASE(3)
bdir2=6
bdir3=5
                g = protcoords(m,l)%fm
                f = protcoords(m,l)%fl

END SELECT

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

if(protcoords(m,l)%a111/=0) then

                         if(int(interenergy(protcoords(g,f)%type,protcoords(m,l)%type))/=0) then

       if(cl(g) .eqv. .true.) goto 31
       a = a + 1
       cllist(a) = g
       cl(g) = .true.
       cluspop = cluspop + 1
    end if

end if


31  continue

  end subroutine clusteraid



  

  subroutine removeenergyintra(m,l,deltaenergy,tempcoord)
    !removes the energy for the bonds that are potentially going to be broken by the move
    integer,intent(in)::m,l
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout)::tempcoord
    integer :: g,f

    g = m
    if(protcoords(m,l)%am == m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%al
       tempcoord(l)%am = -1
    end if

    if(protcoords(m,l)%bm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%bl
       tempcoord(l)%bm = -1
    end if

    if(protcoords(m,l)%cm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
    end if

    if(protcoords(m,l)%dm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%dl
       tempcoord(l)%dm = -1
    end if

    if(protcoords(m,l)%em ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%el
       tempcoord(l)%em = -1
    end if

    if(protcoords(m,l)%fm ==m) then
       deltaenergy = deltaenergy - intraen
       f = protcoords(m,l)%fl
       tempcoord(l)%fm = -1
    end if

  end subroutine removeenergyintra


  subroutine removeenergyinter(m,l,deltaenergy,tempcoord)
    !removes energy for interchain bonds that will potentially be broken by the move
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
    end if

    if((protcoords(m,l)%bm /= 0)  .and. (protcoords(m,l)%bm /= m)) then
       g = protcoords(m,l)%bm
       f = protcoords(m,l)%bl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%bm = -1
    end if

    if((protcoords(m,l)%cm /= 0) .and. (protcoords(m,l)%cm /= m)) then
       g = protcoords(m,l)%cm
       f = protcoords(m,l)%cl
       tempcoord(l)%cm = -1
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
    end if

    if((protcoords(m,l)%dm /= 0) .and. (protcoords(m,l)%dm /= m)) then
       g = protcoords(m,l)%dm
       f = protcoords(m,l)%dl
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%dm = -1
    end if

    if((protcoords(m,l)%em /= 0) .and. (protcoords(m,l)%em /= m)) then
       g = protcoords(m,l)%em
       f = protcoords(m,l)%el
       !if((m == 5) .and. (l == 2) .and. (g == 10) .and. (f == 1)) write(6,*) 'successful bond rupture',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%em = -1
    end if

    if((protcoords(m,l)%fm /= 0) .and. (protcoords(m,l)%fm /= m)) then
       g = protcoords(m,l)%fm
       f = protcoords(m,l)%fl
       !if((m == 10) .and. (l == 1) .and. (g == 5) .and. (f == 2)) write(6,*) 'successful rupture of bond',time
       deltaenergy = deltaenergy - interen(protcoords(m,l)%type,protcoords(g,f)%type)
       tempcoord(l)%fm = -1
    end if

  end subroutine removeenergyinter

  subroutine bondformintra(m,l,st,tempcoord,bdir,deltaenergy)
    !Forms bond around the beads that have been moved (intrachain bonds)
    integer,intent(in)::m,l,st,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'start fail'
    if(bdir == -1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%am = m
       tempcoord(l)%al = st
    end if

    if(bdir == 1) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%bm = m
       tempcoord(l)%bl = st
    end if

    if(bdir == -2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%cm = m
       tempcoord(l)%cl = st
    end if

    if(bdir == 2) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%dm = m
       tempcoord(l)%dl = st
    end if

    if(bdir == -3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%em = m
       tempcoord(l)%el = st
    end if

    if(bdir == 3) then
       deltaenergy = deltaenergy + intraen
       tempcoord(l)%fm = m
       tempcoord(l)%fl = st
    end if
    !if(deltaenergy < 0.0) write(6,*) deltaenergy,intraenergy,'fail'
  end subroutine bondformintra

  subroutine bondforminter(m,g,l,st,tempcoord,bdir,deltaenergy)
    !Forms bonds from the beads that have been moved (inter chain)
    integer,intent(in)::m,l,st,g,bdir
    double precision,intent(inout)::deltaenergy
    type(prottemp),dimension(:),intent(inout):: tempcoord
    !check this and cluster bonds update

       deltaenergy = deltaenergy +interen(protcoords(m,l)%type,protcoords(g,st)%type)


    if(bdir == -1) then
       tempcoord(l)%am = g
       tempcoord(l)%al = st
    end if

    if(bdir == 1) then
       tempcoord(l)%bm = g
       tempcoord(l)%bl = st
    end if

    if(bdir == -2) then
       tempcoord(l)%cm = g
       tempcoord(l)%cl = st
    end if

    if(bdir == 2) then
       tempcoord(l)%dm = g
       tempcoord(l)%dl = st
    end if

    if(bdir == -3) then
       tempcoord(l)%em = g
       tempcoord(l)%el = st
    end if

    if(bdir == 3) then
       tempcoord(l)%fm = g
       tempcoord(l)%fl = st
    end if

    !if(modulo(time,1000)==0) write(6,*) 'interbond',bdir,deltaenergy,tbo(l)
  end subroutine bondforminter

  subroutine updatebondother(m,bead,tempcoord)
    !Subroutine to update bonds and nearest neighbour list for reptation move
    integer,intent(in):: m,bead
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f


    if(tempcoord(bead)%am == -1) then
       g = protcoords(m,bead)%am
       f = protcoords(m,bead)%al
       protcoords(g,f)%bl = 0
       protcoords(g,f)%bm = 0
    end if

    if(tempcoord(bead)%bm == -1) then
       g = protcoords(m,bead)%bm
       f = protcoords(m,bead)%bl
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


  end subroutine updatebondother




  subroutine updatebond(m,beadmin,beadmax,tempcoord)
    !updates bonds and nearest neighbour list 
    integer,intent(in):: m,beadmin,beadmax
    Type(prottemp),dimension(:),intent(inout) :: tempcoord
    integer ::l,g,f,bg,bf


    do l = beadmin,beadmax,1

       if(tempcoord(l)%am == -1) then
          g = protcoords(m,l)%am
          f = protcoords(m,l)%al
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
    !updates bonds and neearest neighbour after reptation (different to updatebondother)
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

    end do

  end subroutine updatebondrep



  subroutine reptation
    !perform reptation move
    double precision :: choose,deltaenergy
    integer :: g,g1,g2,g3,g4,g5,g6,m,l,dxx,dyx,dzx,direc
    integer:: st,pr,dx,dy,dz,chain2,beads2,bdir,maxl,maxback,str
    logical :: reptcont,adjver
    type(prottemp),dimension(:),allocatable :: tempcoord



    t = time + 1

    !choose = -0.5
    choose =  ran2(seed) - 0.5
    !write(6,*) 'reptation start',choose

    reptcont = .true.
    m =int(ran2(seed)*(nprotein))+1

    !if(protcoords(m,1)%species == 3) then
    !call clusterposition
    !return
    !end if


    maxl = chlen(m)

        do l=1,maxl
                if(protcoords(m,l)%a111>0) then
                rerej=rerej+1
                reject = reject+1
                return

                end if
        end do


    allocate(tempcoord(maxl))



    direc = int(ran2(seed)*5) +1
    if (choose > 0.0) then

19     continue

       dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
       dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
       dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))



       if(sqrt(real((dx**2) + (dy**2) + (dz**2))) > protcoords(m,1)%linker) goto 19




       tempcoord(1)%x = modulo(protcoords(m,1)%x + dx-1,gridsize)+1
       tempcoord(1)%y = modulo(protcoords(m,1)%y + dy-1,gridsize)+1
       tempcoord(1)%z = modulo(protcoords(m,1)%z + dz-1,gridsize)+1

       dx = min(abs(tempcoord(1)%x - protcoords(m,2)%x),gridsize-abs(tempcoord(1)%x &
            - protcoords(m,2)%x)) 
       dy = min(abs(tempcoord(1)%y - protcoords(m,2)%y),gridsize- abs(tempcoord(1)%y &
            - protcoords(m,2)%y)) 
       dz = min(abs(tempcoord(1)%z - protcoords(m,2)%z),gridsize - abs(tempcoord(1)%z &
            - protcoords(m,2)%z)) 

       sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
       if(sumdebug > real(protcoords(m,1)%linker)) then
          reptcont = .false.
          goto 83
       end if

       deltaenergy = 0.0d0

       do st = 1,maxl-1,1
          if(overlaps(m,1,st,tempcoord) .eqv. .false.) then
             reptcont = .false.
             goto 83
          end if
       end do



if(cell(tempcoord(1)%x,tempcoord(1)%y,tempcoord(1)%z)%bm /=0) then
reptcont = .false.
goto 83
end if



       tempcoord(maxl)%am = 0
       tempcoord(maxl)%bm = 0
       tempcoord(maxl)%cm = 0
       tempcoord(maxl)%dm = 0
       tempcoord(maxl)%em = 0
       tempcoord(maxl)%fm = 0
       tempcoord(maxl)%al = 0
       tempcoord(maxl)%bl = 0
       tempcoord(maxl)%cl = 0
       tempcoord(maxl)%dl = 0
       tempcoord(maxl)%el = 0
       tempcoord(maxl)%fl = 0

        call removeenergyintra(m,maxl,deltaenergy,tempcoord)
       call removeenergyinter(m,maxl,deltaenergy,tempcoord)

       tempcoord(1)%am = 0
       tempcoord(1)%bm = 0
       tempcoord(1)%cm = 0
       tempcoord(1)%dm = 0
       tempcoord(1)%em = 0
       tempcoord(1)%fm = 0
       tempcoord(1)%al = 0
       tempcoord(1)%bl = 0
       tempcoord(1)%cl = 0
       tempcoord(1)%dl = 0
       tempcoord(1)%el = 0
       tempcoord(1)%fl = 0

       do st = 3,maxl-1
          call adjacent(m,1,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,1,st,tempcoord,bdir,deltaenergy)
          end if
       end do



if(cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=0)then
if(cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=m)then
                pr=cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm  
              str=cell(modulo(tempcoord(1)%x,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if

if(cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=0)then
if(cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm/=m)then
      
          pr=cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bm
                str=cell(modulo(tempcoord(1)%x-2,gridsize)+1,tempcoord(1)%y,tempcoord(1)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if 
end if


        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bm/=0)then
        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bm/=m)then
                pr=cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bm
                str=cell(tempcoord(1)%x,modulo(tempcoord(1)%y,gridsize)+1,tempcoord(1)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bm/=0)then
        if(cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bm/=m)then
  
        pr=cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bm
                str=cell(tempcoord(1)%x,modulo(tempcoord(1)%y-2,gridsize)+1,tempcoord(1)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if
      
  if(cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z,gridsize)+1)%bm/=0)then
          if(cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z,gridsize)+1)%bm/=m)then        
        pr=cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z,gridsize)+1)%bm
                str=cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z,gridsize)+1)%bl
        bdir = -3
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z-2,gridsize)+1)%bm/=0)then
        if(cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z-2,gridsize)+1)%bm/=m)then
pr=cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z-2,gridsize)+1)%bm
                str=cell(tempcoord(1)%x,tempcoord(1)%y,modulo(tempcoord(1)%z-2,gridsize)+1)%bl
        bdir = 3
                      call bondforminter(m,pr,1,str,tempcoord,bdir,deltaenergy)
end if
end if


       !section to adjust for change in types


call entropyrep(deltaenergy,tempcoord,m,1)



       if(Energydecision(deltaenergy) .eqv. .true.) then
  
               !if(deltaenergy> 0.0) write(6,*) 'reptation energy change!!!!!!!!!!!!!!!!!'
               call updatearrayreptation(m,1,maxl,tempcoord)
          totalenergy = totalenergy + deltaenergy
          do l =maxl,2,-1
             protcoords(m,l)%x = protcoords(m,l-1)%x
             protcoords(m,l)%y = protcoords(m,l-1)%y
             protcoords(m,l)%z = protcoords(m,l-1)%z
               cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl =l 
          end do

          call updatebondother(m,maxl,tempcoord)
          call updatepos(m,1,1,tempcoord,maxl)
          call reptationadjust(m,maxl,2,tempcoord,-1)
          call updatebondrep(m,1,1,tempcoord,1)
         successful = successful + 1
          reptforward = reptforward + 1
          reacc = reacc +1      
        call unboundcounter       
     else
          reptcont = .false.
       end if
    else if (choose < 0.0) then


15     continue

       dx = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
       dy = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))
       dz = nint((ran2(seed)-0.5)*((2*protcoords(m,1)%linker)+1))



       if(sqrt(real((dx**2) + (dy**2) + (dz**2))) > protcoords(m,1)%linker) goto 15


       tempcoord(maxl)%x = modulo(protcoords(m,maxl)%x + dx-1,gridsize)+1
       tempcoord(maxl)%y = modulo(protcoords(m,maxl)%y + dy-1,gridsize)+1
       tempcoord(maxl)%z = modulo(protcoords(m,maxl)%z + dz-1,gridsize)+1

       dx = min(abs(tempcoord(maxl)%x - protcoords(m,maxl-1)%x),gridsize-abs(tempcoord(maxl)%x &
            - protcoords(m,maxl-1)%x)) 
       dy = min(abs(tempcoord(maxl)%y - protcoords(m,maxl-1)%y),gridsize- abs(tempcoord(maxl)%y &
            - protcoords(m,maxl-1)%y)) 
       dz = min(abs(tempcoord(maxl)%z - protcoords(m,maxl-1)%z),gridsize - abs(tempcoord(maxl)%z &
            - protcoords(m,maxl-1)%z)) 

       sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
       if(sumdebug > real(protcoords(m,1)%linker)) then
          reptcont = .false.
          goto 83
       end if

       deltaenergy = 0.0d0


       do st = 2,maxl,1
          if(overlaps(m,maxl,st,tempcoord) .eqv. .false.) then
             reptcont = .false.
             goto 83
          end if
       end do



if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm /=0) then
reptcont = .false.
goto 83
end if


       tempcoord(1)%am = 0
       tempcoord(1)%bm = 0
       tempcoord(1)%cm = 0
       tempcoord(1)%dm = 0
       tempcoord(1)%em = 0
       tempcoord(1)%fm = 0
       tempcoord(1)%al = 0
       tempcoord(1)%bl = 0
       tempcoord(1)%cl = 0
       tempcoord(1)%dl = 0
       tempcoord(1)%el = 0
       tempcoord(1)%fl = 0

       call removeenergyintra(m,1,deltaenergy,tempcoord)
       call removeenergyinter(m,1,deltaenergy,tempcoord)

       tempcoord(maxl)%am = 0
       tempcoord(maxl)%bm = 0
       tempcoord(maxl)%cm = 0
       tempcoord(maxl)%dm = 0
       tempcoord(maxl)%em = 0
       tempcoord(maxl)%fm = 0
       tempcoord(maxl)%al = 0
       tempcoord(maxl)%bl = 0
       tempcoord(maxl)%cl = 0
       tempcoord(maxl)%dl = 0
       tempcoord(maxl)%el = 0
       tempcoord(maxl)%fl = 0



       do st = 2,maxl-2
          call adjacent(m,maxl,st,tempcoord,adjver,bdir)
          if((adjver.eqv. .true.)) then
             call bondformintra(m,maxl,st,tempcoord,bdir,deltaenergy)
          end if
       end do


if(cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=0)then
if(cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=m)then

                pr=cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm
                str=cell(modulo(tempcoord(maxl)%x,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bl
        bdir = -1
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if

if(cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=0)then
if(cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm/=m)then
      
          pr=cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bm
                str=cell(modulo(tempcoord(maxl)%x-2,gridsize)+1,tempcoord(maxl)%y,tempcoord(maxl)%z)%bl
        bdir = 1
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if 
end if


        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bm/=0)then
        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bm/=m)then
          pr=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bm
                str=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y,gridsize)+1,tempcoord(maxl)%z)%bl
        bdir = -2
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bm/=0)then
        if(cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bm/=m)then
      
          pr=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bm
                str=cell(tempcoord(maxl)%x,modulo(tempcoord(maxl)%y-2,gridsize)+1,tempcoord(maxl)%z)%bl
        bdir = 2
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if

        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z,gridsize)+1)%bm/=0)then
        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z,gridsize)+1)%bm/=m)then
      
          pr=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z,gridsize)+1)%bm
                str=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z,gridsize)+1)%bl
        bdir = -3
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if


        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z-2,gridsize)+1)%bm/=0)then
        if(cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z-2,gridsize)+1)%bm/=m)then
      
          pr=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z-2,gridsize)+1)%bm
                str=cell(tempcoord(maxl)%x,tempcoord(maxl)%y,modulo(tempcoord(maxl)%z-2,gridsize)+1)%bl
        bdir = 3
                      call bondforminter(m,pr,maxl,str,tempcoord,bdir,deltaenergy)
end if
end if



call entropyrep(deltaenergy,tempcoord,m,maxl)


       if(Energydecision(deltaenergy) .eqv. .true.) then
          totalenergy = totalenergy + deltaenergy
          !if(deltaenergy>0.0) write(6,*) 'reptation energy fail!!!!!!!!!!!'
          call updatearrayreptation(m,maxl,1,tempcoord)
          do l =1,maxl-1,1
             protcoords(m,l)%x = protcoords(m,l+1)%x
             protcoords(m,l)%y = protcoords(m,l+1)%y
             protcoords(m,l)%z = protcoords(m,l+1)%z
               cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl =l

          end do

          call updatebondother(m,1,tempcoord)
          tempcoord(maxl-1)%x = protcoords(m,maxl-1)%x
          tempcoord(maxl-1)%y = protcoords(m,maxl-1)%y
          tempcoord(maxl-1)%z = protcoords(m,maxl-1)%z
          call updatepos(m,maxl-1,maxl,tempcoord,maxl)
          call reptationadjust(m,1,maxl-1,tempcoord,1)
          call updatebondrep(m,maxl,maxl,tempcoord,-1)
          successful = successful + 1
          reptbackward = reptbackward + 1
          reacc = reacc + 1     
                call unboundcounter  
       else
          reptcont = .false.
       end if
    end if
83  if (reptcont .eqv. .false.) then
       reject = reject + 1
       rerej = rerej +1
    end if


  end subroutine reptation



  subroutine reptationadjust(m,beadmin,beadmax,tempcoord,del)
    !moves bond neighbour list along the chain to account for reptation
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
             ! if((del == 1) .and. (protcoords(m,l+del)%al > l+2)) then
             ! protcoords(m,l)%al = protcoords(m,l+del)%al - del
             !else if((del == -1) .and. (protcoords(m,l+del)%al < l-2)) then
             !  protcoords(m,l)%al = protcoords(m,l+del)%al - del
             !else
             if(protcoords(m,l+del)%al /= buzz) then
                protcoords(m,l)%al = protcoords(m,l+del)%al - del
             else
                protcoords(m,l)%am = 0
             end if
      
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
                protcoords(m,l)%bl=0
             end if

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

    end do


  end subroutine reptationadjust





  logical Function overlaps (pr,ch1,ch2,tempcoords)
    !function to check the proposed position doesn't overlap with any other beads
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
    !subroutine to check the proposed cluster move doesn't result in any overlaps
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
    !calculates the energy of the system
    !init true if this is the first time energy is called (in this case
    !bonds/neighbour lists are assigned to every bead
    integer :: dx,dy,dz,m,l,f,g,bdir,delx,dely,delz,maxlengthss,maxl
    double precision :: initialenergy,olderenergy,dumoldenergy
    Type(prottemp),dimension(:),allocatable :: tempcoord
    !double precision,dimension(:),intent(inout)::chen
    logical :: adjver
    logical,intent(in)::init
        integer::bdir2,bdir3
    initialenergy = 0.0d0
    !totalenergy = 0.0d0
    allocate(tempcoord(maxlength))


    if(init .eqv. .true.) then
       do m = 1,nprotein
          maxl = chlen(m)
          do l = 1,maxl
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
             tempcoord(l)%al = 0
             tempcoord(l)%bl = 0
             tempcoord(l)%cl = 0
             tempcoord(l)%dl = 0
             tempcoord(l)%el = 0
             tempcoord(l)%fl = 0
     
        tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y
             tempcoord(l)%z = protcoords(m,l)%z
        call entropyinit1(initialenergy,m,l)
        end do
          do l = 1,maxl,1
             do f = 1,maxl
                if(abs(f-l)>2) then
                   call adjacent(m,l,f,tempcoord,adjver,bdir)
                   if((adjver .eqv. .true.)) then
                      call bondformintra(m,l,f,tempcoord,bdir,initialenergy)
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
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy)
                         !if(initialenergy < olderenergy) write(6,*) 'energyyyyy',&
                         !initialenergy,m,l,g,f
                      end if
                   end do
                end do
             end if
          end do
          call updatebondold(m,1,maxl,tempcoord)
       end do
    end if

    if(init .eqv. .false.) then
       do m=1,nprotein,1
          chen(m)= 0.0d0
       end do
       do m =1,nprotein,1
          maxl = chlen(m)
          do l =1,maxl,1
             tempcoord(l)%am = 0
             tempcoord(l)%bm = 0
             tempcoord(l)%cm = 0
             tempcoord(l)%dm = 0
             tempcoord(l)%em = 0
             tempcoord(l)%fm = 0
             tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y
             tempcoord(l)%z = protcoords(m,l)%z
        call entropyinit(initialenergy,m,l)
          end do
          do l = 1,maxl-3,1
             do f = l+3,maxl
                call adjacent(m,l,f,tempcoord,adjver,bdir)
                if((adjver .eqv. .true.)) then
                   dumoldenergy = initialenergy
                   call bondformintra(m,l,f,tempcoord,bdir,initialenergy)
                   chen(m) = chen(m) + (initialenergy-dumoldenergy)
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
SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
CASE (1)
bdir2=2
bdir3=1
CASE (-2)
bdir2=3
bdir3=4
CASE (2)
bdir2=4
bdir3=3
CASE(-3)
bdir2=5
bdir3=6
CASE(3)
bdir2=6
bdir3=5
END SELECT



                         dumoldenergy = initialenergy
                         if((protcoords(m,l)%a111==bdir2).and.(protcoords(g,f)%a111==bdir3)) then !.and. &
 !                       (protcoords(m,l)%a11n==g) .and.(protcoords(g,f)%a11n==m)) then
                                initialenergy =initialenergy+interenergy(protcoords(m,l)%type,protcoords(g,f)%type)
                           end if
                         call bondforminter(m,g,l,f,tempcoord,bdir,initialenergy)
                         chen(m) = chen(m) + (initialenergy-dumoldenergy)
  
                      end if
                   end do
                end do
             end do
          end if
       end do
    end if


    if((init .eqv. .true.).and.(restart .eqv. .true.)) return


    if(init .eqv. .true.) initialenergy = initialenergy/2
    if((time > 0) .and. (nint(initialenergy) /= nint(totalenergy))) then
        if (abs(initialenergy-totalenergy)<0.1) then
       write(6,*) 'energyfail',totalenergy,initialenergy,time
       finalfail = .true.
        end if
    end if

    if(init .eqv. .true.) totalenergy = initialenergy


  end subroutine energy


  subroutine adjacent(pr,ch1,ch2,tempory,adjver,bdir)
    !subroutine to check whether two beads are on adjacent lattice sites
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


  subroutine comfind(m,track)
    !find the centre of mass of each chain
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
    !linkerlength = 7

    !do a =1,maxl
    !write(6,*) 'bead',a,protcoords(m,a)%x,protcoords(m,a)%y,protcoords(m,a)%z
    !end do
    do l = 2,maxl,1
       dcomx = protcoords(m,l)%x-protcoords(m,l-1)%x
       if (dcomx >protcoords(m,1)%linker) dcomx = dcomx-gridsize
       if (dcomx <(-protcoords(m,1)%linker)) dcomx = gridsize+dcomx
       comx = comx + ((maxl+1-l)*dcomx)
       dcomy = protcoords(m,l)%y-protcoords(m,l-1)%y

       if (dcomy >protcoords(m,1)%linker) dcomy = dcomy-gridsize
       if (dcomy <(-protcoords(m,1)%linker)) dcomy = gridsize+dcomy

       comy = comy + ((maxl+1-l)*dcomy)
       dcomz = protcoords(m,l)%z-protcoords(m,l-1)%z
       if (dcomz >protcoords(m,1)%linker) dcomz = dcomz-gridsize
       if (dcomz <(-protcoords(m,1)%linker)) dcomz = gridsize+dcomz
       comz = comz + ((maxl+1-l)*dcomz)
    end do

    com(1,m)%x = modulo(protcoords(m,1)%x +(real(comx)/(maxl)) -1,real(gridsize))+1
    com(1,m)%y = modulo(protcoords(m,1)%y +(real(comy)/(maxl)) -1,real(gridsize))+1
    com(1,m)%z = modulo(protcoords(m,1)%z +(real(comz)/(maxl)) -1,real(gridsize))+1

    if(time > 0 .and. (track .eqv. .true.)) then
       call tracking(m)
    end if

    ! if(track .eqv. .true.) call rms(m)

  end subroutine comfind

  subroutine tracking(m)
    !find total displacement in x,y and z directions
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



  subroutine rms(m)
    !incorrect rms code - not called
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

  subroutine radiusofgyration
    !calculates the radius of gyration of the largest cluster
    integer::m,l,maxl,nobeads
    double precision :: rog,totrog!,avex,avey,avez
    totrog = 0.0d0
    nobeads = 0
    !avex = 0.0d0
    !avey = 0.0d0
    !avez = 0.0d0
    do m = 1,nprotein,1
       maxl = chlen(m)
       !do l = 1,maxl,1
       !avex = avex + protcoords(m,l)%x
       !avey = avey + protcoords(m,l)%y
       !avez = avez + protcoords(m,l)%z
       !end do
       !avex = avex/maxl
       !avey = avey/maxl
       !avez = avez/maxl

       do l = 1, maxl,1
          rog = (min(modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%x-com(1,m)%x,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%y-com(1,m)%y,gridsize*1.0))))**2+ &
               (min(modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0),(gridsize - &
               modulo(protcoords(m,l)%z-com(1,m)%z,gridsize*1.0))))**2
          totrog = totrog + rog
       end do
       nobeads = nobeads+maxl
    end do

    runningaveROG = runningaveROG + totrog
    polymerrog = polymerrog + (totrog/nprotein)
    if(noinfo .eqv. .false.) write(97,*) time-equilib, &
         totrog/(nobeads),SQRT(totrog/(nobeads))

  end subroutine radiusofgyration


  subroutine length
    !finds end to end separation of chains (only called if scaling info required) scalinginfo = TRUE
    integer :: ddx,ddy,ddz,dxsum,dysum,dzsum,m,l,maxl
    double precision :: normchainl
    t = time
    totchainlength = 0.0
    normchainl = 0.0d0

    do m = 1,nprotein,1

       maxl = chlen(m)
       dxsum = min(abs(protcoords(m,maxl)%x - protcoords(m,1)%x),gridsize-abs(protcoords(m,maxl)%x - protcoords(m,1)%x))
       dysum = min(abs(protcoords(m,maxl)%y - protcoords(m,1)%y),gridsize-abs(protcoords(m,maxl)%y - protcoords(m,1)%y))
       dzsum = min(abs(protcoords(m,maxl)%z - protcoords(m,1)%z),gridsize-abs(protcoords(m,maxl)%z - protcoords(m,1)%z))

       chainlength = ((dxsum**2) + (dysum**2) + (dzsum**2))
       totchainlength = totchainlength + chainlength
    end do
    !totchainlength = sum(chainlength)
    runningaveEtE = runningaveEtE + totchainlength  !sort this out
    normchainl = totchainlength/nprotein
    polymerl = polymerl + normchainl
    !actualchainlength = actualchainlength + avechainlength
    write(91,*) t,normchainl
  end subroutine length

  subroutine dataout
    !outputs system configuration
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit
    integer :: choice
    allocate(bondposit(nprotein,maxlength))
    !write(6,*) 'data'
    t = time
    !if(suffclust .eqv. .false.) choice = 67
    !if(suffclust .eqv. .true.) choice = 88
        choice = 88
    write(choice,*) ' '
    write(choice,*) 'timestep ', 'indexed'
    write(choice,*) 'pbc', real(4*gridsize),real(4*gridsize),real(4*gridsize)
    !write(67,*) ((nprotein1*maxlength1) + (nprotein2*maxlength2))
    acount = 0

    do m = 1,nprotein,1
       maxl = chlen(m)
       do l = 1,maxl,1
          write(choice,*) acount,real(4*protcoords(m,l)%x), real(4*protcoords(m,l)%y),real(4*protcoords(m,l)%z),&
                protcoords(m,l)%a111,protcoords(m,l)%a11n,totalenergy 
          acount = acount +1

       end do
    end do



  end subroutine dataout

  subroutine dataoutrestart
    !outputs final system configuration so restart can be carried out
    integer::m,l,maxl,f,g,maxback,acount
    type(rprot),dimension(:,:),allocatable :: bondposit

    open(59, file = 'checkpoint.xyz', action = 'write')
    write(59,*)backtime,totrmsbrute,suffclust,unboundcount,boundcount,totalenergy
    do m = 1,nprotein,1
       maxl = chlen(m)
       write(59,*) protcoords(m,1)%species,protcoords(m,1)%x,protcoords(m,1)%y,protcoords(m,1)%z,&
            protcoords(m,1)%type,maxl,protcoords(m,1)%linker,protcoords(m,1)%a111,protcoords(m,1)%a11n
       do l = 2,maxl,1
          write(59,*) protcoords(m,l)%species,protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z,protcoords(m,l)%type,&
                protcoords(m,l)%a111,protcoords(m,l)%a11n
       end do
    end do


  end subroutine dataoutrestart


  subroutine clusterenergy(tempcluscoord,clusen,m,clsize,cllist,cl)
    !calculates whether any bonds have been formed, and if they have rejects the move
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
                end if
             end do
          end do
23        continue
       end do
    end do
    !if(time> 1490 .and. (time < 1510)) write(6,*) 'clusterenergy =',energyofcluster
    !write(6,*) 'aaaaaa'
  end subroutine clusterenergy


  subroutine clusterenergycell(tempcluscoord,clusen,m,clsize,cllist,cl)
    !calculates whether any bonds have been formed, and if they have rejects the
    !move
    integer,dimension(:),intent(in) ::cllist
    logical,dimension(:),intent(in)::cl
    integer,intent(in) :: clsize,m
    double precision,intent(inout)::clusen
    type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer:: a,l,bdir,maxl,ra,pr,str
    logical :: clusmove,adjver
    allocate(tempcoord(maxlength))
    !error in this!!!

    !write(6,*)
    !'realsanity',m,l,tempcluscoord(472,4)%x,tempcluscoord(472,4)%y,tempcluscoord(472,4)%z
    clusen = 0.0
    !write(6,*) 'clsize',clsize
    do ra = 1,clsize
       a = cllist(ra)
       !write(6,*) 'a',a
       maxl = chlen(a)
          clusmove = .true.
          do l = 1,maxl,1



if(cell(modulo(tempcluscoord(a,l)%x,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm/=0)then
if(cell(modulo(tempcluscoord(a,l)%x,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm/=a)then
             
   pr=cell(modulo(tempcluscoord(a,l)%x,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm
                str=cell(modulo(tempcluscoord(a,l)%x,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bl
        bdir = -1
          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv. .true.) .and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if

end if
end if

if(cell(modulo(tempcluscoord(a,l)%x-2,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm/=0)then
if(cell(modulo(tempcluscoord(a,l)%x-2,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm/=a)then      
          pr=cell(modulo(tempcluscoord(a,l)%x-2,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bm
                str=cell(modulo(tempcluscoord(a,l)%x-2,gridsize)+1,tempcluscoord(a,l)%y,tempcluscoord(a,l)%z)%bl
        bdir = 1
          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv. .true.).and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if
end if
end if


        if(cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y,gridsize)+1,tempcluscoord(a,l)%z)%bm/=0)then
        if(cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y,gridsize)+1,tempcluscoord(a,l)%z)%bm/=a)then
      
          pr=cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y,gridsize)+1,tempcluscoord(a,l)%z)%bm
                str=cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y,gridsize)+1,tempcluscoord(a,l)%z)%bl
        bdir = -2
          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv..true.).and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if
end if
end if

        if(cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y-2,gridsize)+1,tempcluscoord(a,l)%z)%bm/=0)then
        if(cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y-2,gridsize)+1,tempcluscoord(a,l)%z)%bm/=a)then
      
          pr=cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y-2,gridsize)+1,tempcluscoord(a,l)%z)%bm
                str=cell(tempcluscoord(a,l)%x,modulo(tempcluscoord(a,l)%y-2,gridsize)+1,tempcluscoord(a,l)%z)%bl
        bdir = 2
          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv..true.).and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if

end if
end if

        if(cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z,gridsize)+1)%bm/=0)then
        if(cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z,gridsize)+1)%bm/=a)then
      
          pr=cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z,gridsize)+1)%bm
                str=cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z,gridsize)+1)%bl
        bdir = -3

          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv. .true.).and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if

end if
end if


        if(cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z-2,gridsize)+1)%bm/=0)then
        if(cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z-2,gridsize)+1)%bm/=a)then
                pr=cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z-2,gridsize)+1)%bm
                str=cell(tempcluscoord(a,l)%x,tempcluscoord(a,l)%y,modulo(tempcluscoord(a,l)%z-2,gridsize)+1)%bl
        bdir = 3
          if(cl(pr) .eqv. .false.) then

        call clusbondform(a,pr,l,str,tempcluscoord,bdir)
                   if((isbond .eqv. .true.).and.(interen(protcoords(a,l)%type,protcoords(pr,str)%type) < 0.0)) then
                      clusen = 100
                      return
                   end if
        end if
end if
end if

23        continue
       end do
    end do
    !if(time> 1490 .and. (time < 1510)) write(6,*) 'clusterenergy
    !=',energyofcluster
    !write(6,*) 'aaaaaa'
  end subroutine clusterenergycell


  subroutine clusbondform(m,g,l,f,tempcoord,bdir)
!checks the nearest neighbours of the beads moved in cluster move to check if bonds ahve formed
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



  subroutine clusrbond(m,axis,rotor,maxl)
    !rotate the bonds of the beads in a cluster during cluster rotation
    integer,intent(in):: m,axis,rotor,maxl
    logical :: bm
    integer :: sign,s


    do s = 1,maxl
       bm = .true.


       if(rotor == 0) then
          bm = .true.
          goto 87
       end if
       !write(6,*) bonddd(m,s),m,s
       if(axis == 1) then
          if(rotor == 1) then
             goto 87
          else if(rotor == -1) then
             goto 87
          end if
       else if(axis ==2 ) then
          if(rotor == 1) then
             goto 87
          else if(rotor == -1) then
             goto 87
          end if

       else if(axis ==3) then
          if(rotor == 1) then
             goto 87
          else if(rotor == -1) then
          end if

       end if


87     if(bm .eqv. .true.) continue
    end do
  end subroutine clusrbond

  subroutine clustertranslation(m,clsize,cb,cllist,cl)
    !translation move on a cluster which contains the randomly chosen bead
    integer,intent(inout)::m
    integer,intent(inout):: clsize
    logical,dimension(:),intent(in)::cl
    type(centremass),dimension(:),allocatable:: tcom
    integer,dimension(:),intent(inout) ::cllist
    integer,dimension(:,:),intent(inout) ::cb
    integer :: dx,dy,dz,vv,rvv,a,moved
    type(prottemp),dimension(:,:),allocatable :: tempcluscoord
    double precision :: clusen,direction
    logical :: moveallow
    allocate(tempcluscoord(nprotein,maxlength))
    allocate(tcom(nprotein))
    !return
    !write(6,*) 'translate', totalenergy
    direction = ran2(seed)
    cltatt = cltatt + 1
    moveallow = .true.
    if(direction <= 1.0/6) then
       dx = steplength
       dy = 0
       dz = 0
       moved =1
    else if((direction > 1.0/6) .and. (direction <= 1.0/3)) then
       dx = -(steplength)
       dy = 0
       dz = 0
       moved = 2
    else if((direction > 1.0/3) .and. (direction <= 1.0/2)) then
       dx = 0
       dy = steplength
       dz = 0
       moved = 3
    else if((direction > 1.0/2) .and. (direction <= 2.0/3)) then
       dx = 0
       dy = -(steplength)
       dz = 0
       moved = 4
    else  if((direction > 2.0/3) .and. (direction <= 5.0/6)) then
       dx = 0
       dy = 0
       dz = steplength
       moved = 5
    else if((direction > 5.0/6) .and. (direction <= 1.0/1)) then
       dx = 0
       dy = 0
       dz = -(steplength)
       moved = 6
    end if

    !write(6,*) 'info',dx,dy,dz
    call clusremove(cllist,clsize,tempcluscoord,cl)
    !write(6,*) 'corresponding info',tempcluscoord(464,6)%x,tempcluscoord(464,6)%y,tempcluscoord(464,6)%z

    if(steplength ==1) then
       if(clsize<(nprotein/4)) then
          do rvv = 1,clsize
             vv = cllist(rvv)

             call clustermove1(vv,dx,dy,dz,tempcluscoord,moveallow,cllist,cl,moved)

             if(moveallow .eqv. .false.) return
             tcom(vv)%x = modulo(com(1,vv)%x + dx -1,real(gridsize))+1
             tcom(vv)%y = modulo(com(1,vv)%y + dy -1,real(gridsize))+1
             tcom(vv)%z = modulo(com(1,vv)%z + dz -1,real(gridsize))+1
          end do

       else if(clsize>=(nprotein/4)) then
          do rvv = 1,clsize
             vv = cllist(rvv)
             call clustermove2(vv,dx,dy,dz,tempcluscoord,moveallow,cllist,cl)
             if(moveallow .eqv. .false.) return
             tcom(vv)%x = modulo(com(1,vv)%x + dx -1,real(gridsize))+1
             tcom(vv)%y = modulo(com(1,vv)%y + dy -1,real(gridsize))+1
             tcom(vv)%z = modulo(com(1,vv)%z + dz -1,real(gridsize))+1
          end do
       end if
    else if(steplength /= 1) then
       do rvv = 1,clsize
          vv = cllist(rvv)
          call clustermove2(vv,dx,dy,dz,tempcluscoord,moveallow,cllist,cl)
          if(moveallow .eqv. .false.) return
          tcom(vv)%x = modulo(com(1,vv)%x + dx -1,real(gridsize))+1
          tcom(vv)%y = modulo(com(1,vv)%y + dy -1,real(gridsize))+1
          tcom(vv)%z = modulo(com(1,vv)%z + dz -1,real(gridsize))+1
       end do
    end if


    clusen = 0.0


if((real(clsize)/nprotein)>3.0/4.0)then
    call clusterenergy(tempcluscoord,clusen,m,clsize,cllist,cl)
else
            call clusterenergycell(tempcluscoord,clusen,m,clsize,cllist,cl)
end if

    if(int(clusen) == 0) then
        call updatearrayclus(clsize,cllist,tempcluscoord)
       call updatecluspos(tempcluscoord,tcom,cllist,clsize)
       call updateclusbonds(clsize,cllist,cl,tempcluscoord)
        !call energy(.true.)
       cltacc = cltacc + 1
        call unboundcounter  
  end if


  end subroutine clustertranslation

  subroutine clusremove(cllist,clsize,tempcluscoord,cl)
    !removes neighbours that are not in the same cluster from the neighbour list
    !prior to cluster translation
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


subroutine updatearrayclus(clsize,cllist,tempcluscoord)
    integer :: m,mt,maxl,l
    integer,intent(inout)::clsize
    integer,dimension(:),intent(inout) ::cllist
    type(prottemp),dimension(:,:),intent(in):: tempcluscoord

    do mt = 1,clsize
       m = cllist(mt)
       maxl = chlen(m)
       do l = 1,maxl

cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bm=0
cell(protcoords(m,l)%x,protcoords(m,l)%y,protcoords(m,l)%z)%bl=0

        end do
   end do

    do mt = 1,clsize
       m = cllist(mt)
       maxl = chlen(m)
       do l = 1,maxl

cell(tempcluscoord(m,l)%x,tempcluscoord(m,l)%y,tempcluscoord(m,l)%z)%bm=m
cell(tempcluscoord(m,l)%x,tempcluscoord(m,l)%y,tempcluscoord(m,l)%z)%bl=l

        end do
   end do


end subroutine updatearrayclus

  subroutine updateclusbonds(clsize,cllist,cl,tempcluscoord)
    !Subroutine to update neighbour list after a cluster move
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

  subroutine clustermove1(chainno,delx,dely,delz,tempcluscoord,moveallow,cllist,cl,moved)
    !if the cluster move is only by 1 lattice spacing then overlap check is just checking that none of the
    !beads are neighbouring a bead that isn't in the cluster - in the direction of movement
    integer,intent(in) :: chainno,delx,dely,delz,moved
    logical,dimension(:),intent(in)::cl
    logical,intent(inout)::moveallow
    Type(prottemp),dimension(:,:),intent(inout) :: tempcluscoord
    integer,dimension(:),intent(inout) ::cllist
    integer :: b,pr,st,maxback,maxl
    logical :: clusmove

    clusmove = .true.
    moveallow = .true.
    maxl = chlen(chainno)


    call fastoverlap(chainno,maxl,cllist,clusmove,moved)

    if(clusmove .eqv. .false.) then
       moveallow = .false.
       return
    end if
    do b =1,maxl
       tempcluscoord(chainno,b)%x = modulo(protcoords(chainno,b)%x + delx-1,gridsize)+1
       tempcluscoord(chainno,b)%y = modulo(protcoords(chainno,b)%y + dely-1,gridsize)+1
       tempcluscoord(chainno,b)%z = modulo(protcoords(chainno,b)%z + delz-1,gridsize)+1


    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove1



  subroutine clustermove2(chainno,delx,dely,delz,tempcluscoord,moveallow,cllist,cl)
    !this subroutine is a full overlap check for a cluster step size of any length
    integer,intent(in) :: chainno,delx,dely,delz
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
39        continue
       end do

    end do
31  if(clusmove .eqv. .false.) moveallow = .false.
    !calculate the energy penalty and insure no overlaps

  end subroutine clustermove2


  subroutine fastoverlap(m,maxl,cllist,clusmove,moved)
    !check that none of the beads have neighbours that aren't in
    !the cluster in the direction of movement
    integer,intent(in)::m,maxl,moved
    integer,dimension(:),intent(inout) ::cllist
    logical,intent(inout)::clusmove
    integer::l

    if(moved == 1) then
       do l = 1,maxl
          if(protcoords(m,l)%am /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%am)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    else if(moved == 2) then
       do l = 1,maxl
          if(protcoords(m,l)%bm /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%bm)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    else if(moved == 3) then
       do l = 1,maxl
          if(protcoords(m,l)%cm /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%cm)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    else if(moved == 4) then
       do l = 1,maxl
          if(protcoords(m,l)%dm /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%dm)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    else if(moved == 5) then
       do l = 1,maxl
          if(protcoords(m,l)%em /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%em)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    else if(moved == 6) then
       do l = 1,maxl
          if(protcoords(m,l)%fm /=0) then
             if(ANY(cllist .eq. protcoords(m,l)%fm)) then
                continue
             else
                clusmove = .false.
                return
             end if
          end if
       end do
    end if



  end subroutine fastoverlap



  subroutine debug
    !checks system for overlaps, bbond discrepancies and linker breakages
    integer::m,l,f,g,dx,dy,dz,maxlengthss,maxl,zz


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
             dx = min(abs(protcoords(m,l)%x - protcoords(m,l-1)%x),gridsize-abs(protcoords(m,l)%x &
                  - protcoords(m,l-1)%x))
             dy = min(abs(protcoords(m,l)%y - protcoords(m,l-1)%y),gridsize-abs(protcoords(m,l)%y &
                  - protcoords(m,l-1)%y)) 
             dz = min(abs(protcoords(m,l)%z - protcoords(m,l-1)%z),gridsize-abs(protcoords(m,l)%z &
                  - protcoords(m,l-1)%z)) 

             sumdebug = sqrt(real((dx**2) + (dy**2) + (dz**2)))
             if(sumdebug > protcoords(m,1)%linker) then
                write(6,*) 'SPLIT',m,l,l-1,sumdebug
                finalfail = .true.
             end if
          end if

       end do
    end do





    do m =1,nprotein
       maxl = chlen(m)
       do l=1,maxl


          if(protcoords(m,l)%am /= 0) then
             g = protcoords(m,l)%am
             f = protcoords(m,l)%al
             if((m ==g) .and. (l==f)) then
                write(6,*) 'FAIL -bonded to itself',m,l,g,f,'a'
                finalfail = .true.
             end if
             if((protcoords(g,f)%bm /=m) .or. (protcoords(g,f)%bl /=l)) then
                write(6,*) 'FAIL-bond mismatch',m,l,g,f,'a',protcoords(g,f)%bm,&
                     protcoords(g,f)%bl
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
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
                     protcoords(g,f)%al
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
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
                     protcoords(g,f)%dl
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
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
                     protcoords(g,f)%cl
                finalfail = .true.
             end if
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
                     protcoords(g,f)%fl
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
                finalfail = .true.
             end if
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
                     protcoords(g,f)%el
                write(6,*) 'protein a',protcoords(m,l)
                write(6,*) 'protein b',protcoords(g,f)
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
!subroutine that allocates all chains to clusters and finds the largest cluster for analysis
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn,pr,str
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize
    integer,dimension(:,:),allocatable::cb,conn
    integer,dimension(:),allocatable::maxcluslist,concount,correlationone,correlationcluster
    double precision::maxclusenergy
!real :: time1,time2
    allocate(cb(nprotein,nprotein))
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(conn(nprotein,maxlength))
    allocate(concount(nprotein))
    allocate(correlationone(nprotein))
    allocate(correlationcluster(nprotein))
    correlationone(:) = 0
    correlationcluster(:) = 0
    !allocate(mashist(10))
    maxclusenergy = 0.0d0
    normalisedsize = 0
    clcount = 0
    clusterpop = 0
    clustcount = 0
    do m = 1,nprotein
       clnos(m) = m
       cllist(m) = 0
       histcl(m) = 0
       concount(m) =0
       do a =1,maxlength
          conn(m,a)=0
       end do
    end do

    a = 1

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do


!call CPU_TIME(time1)

    do m = 1,nprotein-1,1
       maxback = chlen(m)
          do l = 1,maxback
             tempcoord(l)%x = protcoords(m,l)%x
             tempcoord(l)%y = protcoords(m,l)%y

if(cell(modulo(protcoords(a,l)%x,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm/=0)then
if(cell(modulo(protcoords(a,l)%x,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm/=m)then
             
   pr=cell(modulo(protcoords(a,l)%x,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm
                str=cell(modulo(protcoords(a,l)%x,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bl
        bdir = -1
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)

end if
end if

if(cell(modulo(protcoords(a,l)%x-2,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm/=0)then
if(cell(modulo(protcoords(a,l)%x-2,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm/=m)then
      
          pr=cell(modulo(protcoords(a,l)%x-2,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bm
                str=cell(modulo(protcoords(a,l)%x-2,gridsize)+1,protcoords(a,l)%y,protcoords(a,l)%z)%bl
        bdir = 1
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)

end if
end if


        if(cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y,gridsize)+1,protcoords(a,l)%z)%bm/=0)then
        if(cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y,gridsize)+1,protcoords(a,l)%z)%bm/=m)then
      
          pr=cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y,gridsize)+1,protcoords(a,l)%z)%bm
                str=cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y,gridsize)+1,protcoords(a,l)%z)%bl
        bdir = -2
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)

end if
end if

        if(cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y-2,gridsize)+1,protcoords(a,l)%z)%bm/=0)then
        if(cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y-2,gridsize)+1,protcoords(a,l)%z)%bm/=m)then
                pr=cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y-2,gridsize)+1,protcoords(a,l)%z)%bm
                str=cell(protcoords(a,l)%x,modulo(protcoords(a,l)%y-2,gridsize)+1,protcoords(a,l)%z)%bl
        bdir = 2
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)


end if
end if

        if(cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z,gridsize)+1)%bm/=0)then
        if(cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z,gridsize)+1)%bm/=m)then
      
          pr=cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z,gridsize)+1)%bm
                str=cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z,gridsize)+1)%bl
        bdir = -3
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)

end if
end if

        if(cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z-2,gridsize)+1)%bm/=0)then
        if(cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z-2,gridsize)+1)%bm/=m)then
      
          pr=cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z-2,gridsize)+1)%bm
                str=cell(protcoords(a,l)%x,protcoords(a,l)%y,modulo(protcoords(a,l)%z-2,gridsize)+1)%bl
        bdir = 3
call clusterform(cb,m,l,pr,str,bdir,conn,correlationone,concount,clusyes,clnos)

end if
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
    allocate(maxcluslist(maxclus))

    z = 1
    atconn= 0
    do m = 1,nprotein
       if(protcoords(m,1)%species==1) atconn = atconn+concount(m)
       if(clnos(m) == content)then
          maxcluslist(z) = m
          z = z+1
correlationcluster(m) = 1
       end if
    end do

    cconn =0
    acconn =0

    do m = 1,maxclus   
       g = maxcluslist(m)
       cconn = cconn + concount(g)
       if(protcoords(g,1)%species ==1) acconn = acconn+concount(g)
       maxclusenergy = maxclusenergy + chen(g)
    end do

    tconn=sum(concount(:))

    if((maxclus >10) .and. (suffclust .eqv. .false.)) then
       suffclust = .true.
       !call pdbsetupalt
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

    ! write(6,*) 'tconn',tconn,'atconn',atconn,'acconn',acconn,'cconn',cconn,maxclus
    onepop = 0
    do m = 1,maxclus
       g = maxcluslist(m)
       if(protcoords(g,1)%species ==1) onepop = onepop+1
    end do

    if(restart .eqv. .true.) then
       do m = 1,nprotein
          if(histcl(m)/=0) write(38,*) time,histcl(m),histcl(m)**2
       end do
    end if

    !if(noinfo .eqv. .false.) then
    !if(clustcount /= 0) write(6,*) time,clustcount, 'cluster!'
    if(clustcount /= 0) write(82,*) time,clustcount, clusterpop, real(clusterpop)/clustcount &
         , real((normalisedsize) +(nprotein-clusterpop))/nprotein,maxclus,onepop,maxclus-onepop,maxclusenergy&
         ,sum(chen)-maxclusenergy,acconn,cconn-acconn,atconn-acconn,tconn-((atconn-acconn)+cconn) 

    !,clcount
    !if(time>10000) write(6,*) 'output',time, real((normalisedsize)+(nprotein-clusterpop))/nprotein,real((normalisedsize)+ &
    !(nprotein-clusterpop))
    if(clustcount == 0)  write(82,*) time,clustcount, clusterpop,0.0,0.0,0.0,0,0,0,0,0,0,0,0
    !end if

    write(886,*) correlationone(:)
    write(887,*) correlationcluster(:)

    !write(6,*) 'RELIABLE CLUSTER COUNT',clustcount,time2-time1
    call phase(content,maxclus,clnos,cb)
    deallocate(cllist)
    deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
  end subroutine clustercount

subroutine clusterform(cb,m,l,g,f,bdir,conn,correlationone,concount,clusyes,clnos)
integer,intent(in)::m,l,g,f,bdir
logical,intent(inout)::clusyes
integer::oldcl,zzz,bdir2,bdir3
    integer,dimension(:),intent(inout):: clnos
    integer,dimension(:,:),intent(inout)::cb,conn
    integer,dimension(:),intent(inout)::concount,correlationone




                        cb(m,g) = 1
                      cb(g,m) = 1
                      if(interen(protcoords(m,l)%type,protcoords(g,f)%type) <0.0) then

                         if(any(conn(m,:) .eq. g)) then
                            goto 63
                         else
                            concount(m) = concount(m) + 1
                            conn(m,concount(m)) = g
                         end if
63                       continue

                         if(any(conn(g,:) .eq. m)) then
                            goto 83
                         else
                            concount(g) = concount(g) + 1
                            conn(g,concount(g)) = m
                         end if
83                       continue

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
SELECT CASE (bdir)

CASE (-1)
bdir2=1
bdir3=2
CASE (1)
bdir2=2
bdir3=1
CASE (-2)
bdir2=3
bdir3=4
CASE (2)
bdir2=4
bdir3=3
CASE(-3)
bdir2=5
bdir3=6
CASE(3)
bdir2=6
bdir3=5
END SELECT



                         if((protcoords(m,l)%a111==bdir2).and.(protcoords(g,f)%a111==bdir3).and. &
                        (protcoords(m,l)%a11n==g).and.(protcoords(g,f)%a11n==m) .and. &
                        (int(interenergy(protcoords(g,f)%type,protcoords(m,l)%type))/= 0)) then

if(m == 1) correlationone(g) = 1

                      if(any(conn(m,:) .eq. g)) then
                         goto 65
                      else
                         concount(m) = concount(m) + 1
                         conn(m,concount(m)) = g
                      end if
65                    continue

                      if(any(conn(g,:) .eq. m)) then
                         goto 85
                      else
                         concount(g) = concount(g) + 1
                         conn(g,concount(g)) = m
                      end if
85                    continue

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

end subroutine clusterform

  subroutine strippedclustercount
!subroutine that allocates all chains to clusters and finds the largest cluster for analysis
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn,bdir2,bdir3
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl,bound
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize
    integer,dimension(:,:),allocatable::cb,conn
    integer,dimension(:),allocatable::maxcluslist,concount,correlationone,correlationcluster
    double precision::maxclusenergy
    !    real::time1,time2    
allocate(cb(nprotein,nprotein))
    allocate(clnos(nprotein))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(correlationone(nprotein))
    allocate(correlationcluster(nprotein))
    correlationone(:) = 0
    correlationcluster(:) = 0
    !allocate(mashist(10))
    maxclusenergy = 0.0d0
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

!call CPU_TIME(time1)

    do m = 1,nprotein-1
       maxback = chlen(m)
       do l = 1,maxback
          bound = .false.

if(protcoords(m,l)%a111/=0) then
bdir=protcoords(m,l)%a111

SELECT CASE (bdir)

CASE (1)
bdir2=1
bdir3=2 
                g = protcoords(m,l)%am
                f = protcoords(m,l)%al

CASE (2)
bdir2=2
bdir3=1
                g = protcoords(m,l)%bm
                f = protcoords(m,l)%bl

CASE (3)
bdir2=3
bdir3=4
                g = protcoords(m,l)%cm
                f = protcoords(m,l)%cl

CASE (4)
bdir2=4
bdir3=3
                g = protcoords(m,l)%dm
                f = protcoords(m,l)%dl

CASE(5)
bdir2=5
bdir3=6
                g = protcoords(m,l)%em
                f = protcoords(m,l)%el

CASE(6)
bdir2=6
bdir3=5
                g = protcoords(m,l)%fm
                f = protcoords(m,l)%fl

END SELECT


if(int(interenergy(protcoords(g,f)%type,protcoords(m,l)%type))/=0) then
bound = .true. 



          
          if((bound .eqv. .true.) .and. (g>m)) then
             if (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0) then
                if(m == 1) correlationone(g) = 1
                
                
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
          end if
      end if 
    
             
       end do
    end do
    
!call CPU_TIME(time2)

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
    allocate(maxcluslist(maxclus))

    z = 1
    
    do m = 1,nprotein
       if(clnos(m) == content)then
          maxcluslist(z) = m
          z = z+1
          correlationcluster(m) = 1
       end if
    end do


    write(886,*) correlationone(:)
    write(887,*) correlationcluster(:)

    !write(6,*) 'TEST cLUSTER COUNT',clustcount,time2-time1
   
    deallocate(cllist)
    deallocate(clnos)
    !deallocate(tempcoord)
    deallocate(histcl)
  end subroutine strippedclustercount






  subroutine clustercount2
!subroutine that allocates all chains to clusters and finds the largest cluster for analysis
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,a,content,onepop,z
    integer :: acconn,cconn,atconn,tconn,bdir2,bdir3
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl,bound
    type(prottemp),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize
    integer,dimension(:,:),allocatable::cb,conn
    integer,dimension(:),allocatable::maxcluslist,concount,correlationone,correlationcluster
    double precision::maxclusenergy
!real :: time1,time2
    allocate(cb(nprotein,nprotein))
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(conn(nprotein,maxlength))
    allocate(concount(nprotein))
    allocate(correlationone(nprotein))
    allocate(correlationcluster(nprotein))
    correlationone(:) = 0
    correlationcluster(:) = 0
    !allocate(mashist(10))
    maxclusenergy = 0.0d0
    normalisedsize = 0
    clcount = 0
    clusterpop = 0
    clustcount = 0
    do m = 1,nprotein
       clnos(m) = m
       cllist(m) = 0
       histcl(m) = 0
       concount(m) =0
       do a =1,maxlength
          conn(m,a)=0
       end do
    end do

    a = 1

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do


!call CPU_TIME(time1)


do m = 1,nprotein-1
       maxback = chlen(m)
       do l = 1,maxback
          bound = .false.

!!!!!!!!!!!!!!!!!!REPLACETHIS!!!!!!!!!!!!!!
if(protcoords(m,l)%a111/=0) then
bdir=protcoords(m,l)%a111

SELECT CASE (bdir)

CASE (1)
bdir2=1
bdir3=2
                g = protcoords(m,l)%am
                f = protcoords(m,l)%al

CASE (2)
bdir2=2
bdir3=1
                g = protcoords(m,l)%bm
                f = protcoords(m,l)%bl

CASE (3)
bdir2=3
bdir3=4
                g = protcoords(m,l)%cm
                f = protcoords(m,l)%cl

CASE (4)
bdir2=4
bdir3=3
                g = protcoords(m,l)%dm
                f = protcoords(m,l)%dl

CASE(5)
bdir2=5
bdir3=6
                g = protcoords(m,l)%em
                f = protcoords(m,l)%el

CASE(6)
bdir2=6
bdir3=5
                g = protcoords(m,l)%fm
                f = protcoords(m,l)%fl

END SELECT


bound = .true.
end if

          
          if((bound .eqv. .true.) .and. (g>m)) then

                        if (interenergy(protcoords(g,f)%type,protcoords(m,l)%type)  < 0.0) then

if(m == 1) correlationone(g) = 1
                      
                      if(any(conn(m,:) .eq. g)) then
                         goto 65
                      else
                         concount(m) = concount(m) + 1
                         conn(m,concount(m)) = g
                      end if
65                    continue

                      if(any(conn(g,:) .eq. m)) then
                         goto 85
                      else
                         concount(g) = concount(g) + 1
                         conn(g,concount(g)) = m
                      end if
85                    continue

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
        

!call CPU_TIME(time2)

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
    allocate(maxcluslist(maxclus))

    z = 1
    atconn= 0
    do m = 1,nprotein
       if(protcoords(m,1)%species==1) atconn = atconn+concount(m)
       if(clnos(m) == content)then
          maxcluslist(z) = m
          z = z+1
correlationcluster(m) = 1
       end if
    end do

    cconn =0
    acconn =0

    do m = 1,maxclus   
       g = maxcluslist(m)
       cconn = cconn + concount(g)
       if(protcoords(g,1)%species ==1) acconn = acconn+concount(g)
       maxclusenergy = maxclusenergy + chen(g)
    end do

    tconn=sum(concount(:))

    if((maxclus >10) .and. (suffclust .eqv. .false.)) then
       suffclust = .true.
       call pdbsetupalt
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

    ! write(6,*) 'tconn',tconn,'atconn',atconn,'acconn',acconn,'cconn',cconn,maxclus
    onepop = 0
    do m = 1,maxclus
       g = maxcluslist(m)
       if(protcoords(g,1)%species ==1) onepop = onepop+1
    end do

    if(restart .eqv. .true.) then
       do m = 1,nprotein
          if(histcl(m)/=0) write(38,*) time,histcl(m),histcl(m)**2
       end do
    end if

    !if(noinfo .eqv. .false.) then
    !if(clustcount /= 0) write(6,*) time,clustcount, 'cluster!'
    if(clustcount /= 0) write(82,*) time,clustcount, clusterpop, real(clusterpop)/clustcount &
         , real((normalisedsize) +(nprotein-clusterpop))/nprotein,maxclus,onepop,maxclus-onepop,maxclusenergy&
         ,sum(chen)-maxclusenergy,acconn,cconn-acconn,atconn-acconn,tconn-((atconn-acconn)+cconn) 

    !,clcount
    !if(time>10000) write(6,*) 'output',time, real((normalisedsize)+(nprotein-clusterpop))/nprotein,real((normalisedsize)+ &
    !(nprotein-clusterpop))
    if(clustcount == 0)  write(82,*) time,clustcount, clusterpop,0.0,0.0,0.0,0,0,0,0,0,0,0,0
    !end if

    write(886,*) correlationone(:)
    write(887,*) correlationcluster(:)

    !write(6,*) 'RELIABLE CLUSTER COUNT',clustcount,time2-time1
    call phase(content,maxclus,clnos,cb)
    deallocate(cllist)
    deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
  end subroutine clustercount2

  

  subroutine clustercom(m,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge,dpcomcluster)
    !finds the centre of mass of a cluster
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
    type(rprot),intent(inout) :: dpcomcluster

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
       !write(6,*) 'a'
       !if(clsize>1) write(6,*) 'component',b,com(1,b)%x,com(1,b)%y,com(1,b)%z
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
             if (delx > (gridsize)/2) delx = delx-gridsize
             if (delx < -(gridsize)/2) delx = delx + gridsize
             comx(b,f) =  delx
             comx(f,b) = -delx
             dely = com(1,f)%y - com(1,b)%y
             if (dely > (gridsize)/2) dely = dely-gridsize
             if (dely < -(gridsize)/2) dely = dely + gridsize
             comy(b,f) =  dely
             comy(f,b) = -dely
             delz = com(1,f)%z - com(1,b)%z
             if (delz > (gridsize)/2) delz = delz-gridsize
             if (delz < -(gridsize)/2) delz = delz + gridsize
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
       !write(6,*) 'checklist',checklist(dint)

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
             !write(6,*) 'rcounts',rcounter(b),rcounter(dummyint),a,cench
             checklist(zx) = a
             zx = zx +1
          end if
16        if(a == a) continue
       end do

    end do


    !write(6,*) 'backup 2',rcounter(4)
    !write(6,*) 'orig central chain',centralchain,toolarge
    call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
    !write(6,*) 'key'


    !write(6,*) 'backup 3',rcounter(4)

76  if(clsize == 1) then
       maxl = chlen(m)
       centralchain = m
       totcluspop = maxl
       ! write(6,*) 'individual chain',com(1,m)%x,com(1,m)%y,com(1,m)%z
    end if

    !if(clsize>1) write(6,*) 'totcluspop',totcluspop
    !write(6,*) centralchain
    !write(6,*) 'runtime', com(1,centralchain)%x
    !write(6,*) dummycom%x
    comcluster%x = INT(modulo((com(1,centralchain)%x) + (dummycom%x/totcluspop)-1,real(gridsize))+1)
    comcluster%y = INT(modulo((com(1,centralchain)%y) + (dummycom%y/totcluspop)-1,real(gridsize))+1)
    comcluster%z = INT(modulo((com(1,centralchain)%z) + (dummycom%z/totcluspop)-1,real(gridsize))+1)


    dpcomcluster%x = modulo((com(1,centralchain)%x) + (dummycom%x/totcluspop)-1,real(gridsize))+1
    dpcomcluster%y = modulo((com(1,centralchain)%y) + (dummycom%y/totcluspop)-1,real(gridsize))+1
    dpcomcluster%z = modulo((com(1,centralchain)%z) + (dummycom%z/totcluspop)-1,real(gridsize))+1
    write(3,*) time,dpcomcluster%x,dpcomcluster%y,dpcomcluster%z

53  if(compass .eqv. .true.) continue

  end subroutine clustercom


  subroutine comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
!finds the centre of mass of a large cluster
    integer,intent(in)::clsize
    integer,intent(inout),dimension(:,:)::route
    integer,intent(inout),dimension(:)::checklist,rcounter
    integer :: zna,a,b,rc,maxb
    double precision,intent(in),dimension(:,:) ::comx,comy,comz
    type(rprot),intent(inout) :: dummycom
    logical,intent(inout)::toolarge
    double precision :: bx,by,bz,maxx,maxy,maxz,minx,miny,minz
    !write(6,*) 'centralchain',checklist(1)


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

       ! write(6,*) 'dummycom', dummycom%x,dummycom%y,dummycom%z
       !write(6,*) 'x',comx(checklist(1),route(a,1)),checklist(1),route(a,1),&
       !  com(1,checklist(1))%x,com(1,route(a,1))%x
       ! write(6,*) 'y',comy(checklist(1),route(a,1)),checklist(1),route(a,1),&
       !      com(1,checklist(1))%y,com(1,route(a,1))%y
       !write(6,*) 'z',comz(checklist(1),route(a,1)),checklist(1),route(a,1),&
       !    com(1,checklist(1))%z,com(1,route(a,1))%z


       if(rcounter(a)>1) then
          do rc = 1,rcounter(a)-1,1
             !write(6,*) 'finalcoms', comz(centralchain,route(a,rc))
             !write(6,*) 'route', route(zna,rc)
             !write(6,*) 'total count =',rcounter(zna)
             bx  = bx + comx(route(a,rc),route(a,rc+1))
             by = by + comy(route(a,rc),route(a,rc+1))
             bz = bz + comz(route(a,rc),route(a,rc+1))
             !write(6,*) 'x',comx(route(a,rc),route(a,rc+1))
             !write(6,*) 'y',comy(route(a,rc),route(a,rc+1))
             !write(6,*) 'z',comz(route(a,rc),route(a,rc+1))
             !write(6,*) 'dummycom', dummycom%x,dummycom%y,dummycom%z
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
    !check to see calculated centre of mass is correct
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


  subroutine pdbsetup
    !sets up vtf file for trajectory
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
               'resid ',protcoords(m,l)%type
          write(67,*) 'atom', acount+1, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species)-1,'resid ', &
               protcoords(m,l)%type
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


  subroutine pdbsetupalt
    !new subroutine to setup vtf file for trajectory output
    integer::m,l,maxl,acount,dumres,infoout
    character(len=8) :: atnum
    acount = 0
    write(88,*) 'atom ', 'default ', 'radius ', 1.00000, 'name ','C'
    do m= 1,nprotein
       maxl = chlen(m)
       !if(protcoords(m,1)%species == 1) dumres = 1
       !if(protcoords(m,1)%species == 2) dumres = 2
       !if(dumres == 1) then
       !if(clnos(m) == clnos(g)) then
       !   infoout = 100
       !else
       !   infoout = protcoords(m,1)%species
       !end if


       do l = 1,maxl
          write(88,*) 'atom', acount, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species), &
               'resid ',protcoords(m,l)%type
          write(88,*) 'atom', acount+1, 'radius', 1.00000, 'name ',(2*protcoords(m,1)%species)-1,'resid ', &
               protcoords(m,l)%type
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
    !reads in data from setup file
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




    seed = 6                      !seed for random number generation
    maxlength = 20              !number of beads in chain
    nprotein = 1                 !number of chains
    gridsize= 100        !dimensions of the lattice
    maxtime= 10000                !length of simulation
    equilib = 5000              !equilibration time
    rept = 1                  !1 = reptation 0 = no reptations
    kT = 1.0
    film = .true.
    debugyes = .true.



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
          call get_integer(dum)   
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
          maxtime = maxtime
       CASE ('EQUIB_TIME')
          CALL get_integer(equilib)
          equilib = equilib
       CASE ('RIGHTANGLE')
          !CALL get_integer(right)

       CASE ('CRANK')
          !CALL get_integer(crank)

       CASE ('REPTATION')
          CALL get_integer(rept)

       CASE ('PIVOT')
          !CALL get_integer(piv)

       CASE ('KT')
          CALL get_dp(kT)

       CASE ('FILM')
          CALL get_logical(film)

       CASE ('DEBUGYES')
          CALL get_logical(debugyes)

          !CASE ('ISBOND')
          !CALL get_logical(isbond)

       CASE ('MAXPIV')
          !CALL get_integer(maxpiv)

       CASE ('CRANKLIMIT')
          !CALL get_integer(cranklimit)



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
          !call get_integer(linkerlength)

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



