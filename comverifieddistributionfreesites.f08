program move

  implicit none

  character(len=10)::BIN
  integer,dimension(:),allocatable::type,bulkneigh,clusneigh
  double precision,dimension(:),allocatable::coordx,coordy,coordz,en,sum1
  integer::totpoints,n,m,l,t,j,lx,timeofinterest,dummymax,count,e,maxlength,kx,t1,t2
  double precision::bondx,bondy,bondz,rx,ry,rz,readcomx,readcomy,readcomz
  integer,dimension(:),allocatable ::chlen,speciespop,specieslen,specieslink,totunb
  double precision::delta,r,dr,rho,intraen,totiso,sumall
  real,dimension(:,:),allocatable::summing,radial,allsum
integer::binsize,maxl,gridsize,unusedtime,maxtime,nprotein,mtype,nspecies,preinterest,time
integer::firstrun,ft,gt,x,halfbox
  character(len = 10) :: commandread,commandread2,commandread3,commandread4,commandread5
  character(len=10) :: u,v,w
  real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy
  logical ::isbond  
integer:: norma,normb,normc,normd
integer,dimension(:,:),allocatable :: runningclusterlist
integer,dimension(6)::totclhist,totblhist
  type rprot
     real :: x,y,z
     integer::bond,type
  end type rprot


  type rprot2
     real :: a1,a2,a3,a4,a5,a6,linker
     integer::x,y,z,bond,type
  end type rprot2


  type listt
     integer:: m,l
  end type listt


  type centremass
     double precision :: x,y,z
  end type centremass

type basicp
integer::x,y,z
end type

  type(centremass),dimension(:,:),allocatable :: com

  type(listt),dimension(:),allocatable:: abulk,bbulk,aclus,bclus,a2bulk,a2clus,b2clus,b2bulk
  type(rprot2),dimension(:,:),allocatable :: chains

open(211,file = 'connectivitytracker',action = 'write')
open(712,file= 'inputp.dat',action='write')
  open(21, file = 'm.vtf', action = 'read')
  open(29, file = 'RDF.dat', action = 'write')
  open(31, file = 'unnormRDF.dat', action = 'write')
  open(39, file='e.dat',action='read')
  open(38, file = 'freesitedist.dat', action = 'write')
  open(47, file='c.dat',action='read')
  open(68,file ='unbounddensity.dat',action='write')
open(168,file='densitybeadspershell.dat',action='write')  
open(196,file='availablesites.dat',action='write')
open(99,file = 'stability.dat',action = 'write')
close(99)
!open(196,file='availablesites.dat',action = 'write')

!write(168,*) 'FAIL'
!write(6,*) 'FAAAAAAAAAAAAAAAAAAAAAAIIIIIIIIIIIIIIIIIILLLLLLLLLLLLLLLLLLL'
isbond = .true.
  call read_setup
  write(6,*) 'mytype',mtype,'int 0',int(0.8),int(1.2),int(-0.2),int(-1.2)

  do ft = 1,mtype
     do gt = 1,mtype
        totiso = totiso + interen(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
  end do
  if(nint(totiso) == 0) isbond = .false.
isbond = .false.
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
  call get_command_argument(4,commandread4)
  read(commandread4, '(i4)') preinterest
  call get_command_argument(5,commandread5)
  read(commandread5, '(i4)') firstrun

  write(6,*) 'inputs', binsize,unusedtime,timeofinterest,preinterest
  if(specieslen(1)>=specieslen(2)) maxlength = specieslen(1)
  if(specieslen(1)<specieslen(2)) maxlength = specieslen(2)


halfbox = int((sqrt(3*(real(gridsize)**2)))/2)

  allocate(summing(4*mtype,((halfbox)*binsize)+2))
  allocate(radial(4*mtype,((halfbox)*binsize)+2))
allocate(allsum(mtype,((halfbox)*binsize)+2))
  allocate(com(2,nprotein))
  allocate(coordx(totpoints/2))
  allocate(coordy(totpoints/2))
  allocate(coordz(totpoints/2))
  allocate(chains(nprotein,maxlength))
  allocate(aclus(specieslen(1)*speciespop(1)))
  allocate(bclus(specieslen(2)*speciespop(2)))
  allocate(abulk(specieslen(1)*speciespop(1)))
  allocate(bbulk(specieslen(2)*speciespop(2)))
    allocate(a2clus(specieslen(1)*speciespop(1)))
  allocate(b2clus(specieslen(2)*speciespop(2)))
  allocate(a2bulk(specieslen(1)*speciespop(1)))
  allocate(b2bulk(specieslen(2)*speciespop(2)))
  allocate(en(timeofinterest))
  allocate(chlen(nprotein))
allocate(bulkneigh(mtype))
allocate(clusneigh(mtype))
allocate(totunb(4))
  allocate(type(totpoints/2))
  allocate(sum1(4*mtype))
allocate(runningclusterlist(nprotein,2))

clusneigh = 0
bulkneigh = 0
  summing(:,:) = 0
radial(:,:) = 0
allsum(:,:) = 0
totunb(:) = 0
  count = 1
totclhist(:) = 0
totblhist(:) = 0

  read(21,*) BIN,BIN,BIN,BIN,BIN,BIN
  do m = 1,nprotein
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,chains(m,1)%type
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
chains(m,1)%linker = specieslink(chains(m,1)%type)
     maxl = specieslen(chains(m,1)%type)
     chlen(m) = specieslen(chains(m,1)%type)
     count = count+1

     do l = 2,maxl
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,chains(m,l)%type
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        count = count+1
chains(m,l)%linker = specieslink(chains(m,l)%type)

     end do

  end do


  dummymax = totpoints - nprotein
  do m = 1,dummymax
     read(21,*)BIN,BIN,BIN
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !accounts for 2 extra lines in energy than in move
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do m = 1,preinterest
!  read(39,*) BIN,BIN
!read(47,*) BIN,BIN
read(47,*)BIN,BIN
end do

do m = 1,firstrun,1
     !read(21,*) BIN
     !read(21,*) BIN
     read(39,*) BIN,BIN
     read(47,*) BIN,BIN,BIN
     !do  l= 1,totpoints
        !read(21,*)BIN
     !end do
   end do
write(6,*) 'preinterest =',preinterest,firstrun
   
  do t = 1,unusedtime
     read(21,*) BIN
     read(21,*) BIN
     read(39,*) BIN,BIN
     read(47,*) BIN,BIN,BIN,BIN
     !read(21,*) BIN,BIN,BIN
     do m = 1,totpoints
        read(21,*)BIN
     end do
  end do


  
  do t = 1,timeofinterest,1

     read(21,*) BIN
     read(21,*) BIN
     read(39,*) t1,en(t)
     read(47,*)t2,readcomx,readcomy,readcomz
write(6,*) 'time',t1,t2
     if(t1 /= t2) stop 8
     write(6,*) 'energy',en(t),t,t1     
     do m = 1,nprotein,1
        read(21,*)  BIN,rx,ry,rz
        read(21,*) BIN,bondx,bondy,bondz

        chains(m,1)%x = rx/4
        chains(m,1)%y = ry/4
        chains(m,1)%z = rz/4


        if((chains(m,1)%x-(bondx/4.0) == -0.5) .or. (chains(m,1)%x-(bondx/4.0) == (gridsize-0.5))) chains(m,1)%bond = 1
        if((chains(m,1)%x-(bondx/4.0) == 0.5) .or. (chains(m,1)%x-(bondx/4.0) == -(gridsize-0.5))) chains(m,1)%bond = -1
        if((chains(m,1)%y-(bondy/4.0) == -0.5) .or. (chains(m,1)%y-(bondy/4.0) == (gridsize-0.5))) chains(m,1)%bond = 2
        if((chains(m,1)%y-(bondy/4.0) == 0.5) .or. (chains(m,1)%y-(bondy/4.0) == -(gridsize-0.5))) chains(m,1)%bond = -2
        if((chains(m,1)%z-(bondz/4.0) == -0.5) .or. (chains(m,1)%z-(bondz/4.0) == (gridsize-0.5))) chains(m,1)%bond = 3
        if((chains(m,1)%z-(bondz/4.0) == 0.5) .or. (chains(m,1)%z-(bondz/4.0) == -(gridsize-0.5))) chains(m,1)%bond = -3
        maxl = specieslen(chains(m,1)%type)

        if(chains(m,1)%bond ==0) write(6,*) 'no bond allocation',m,1
        do l = 2,maxl,1

           read(21,*)  BIN,rx,ry,rz
           read(21,*) BIN,bondx,bondy,bondz

           chains(m,l)%x = rx/4
           chains(m,l)%y = ry/4
           chains(m,l)%z = rz/4

           if((chains(m,l)%x-(bondx/4.0) == -0.5) .or. (chains(m,l)%x-(bondx/4.0) == (gridsize-0.5))) chains(m,l)%bond = 1
           if((chains(m,l)%x-(bondx/4.0) == 0.5) .or. (chains(m,l)%x-(bondx/4.0) == -(gridsize-0.5))) chains(m,l)%bond = -1
           if((chains(m,l)%y-(bondy/4.0) == -0.5) .or. (chains(m,l)%y-(bondy/4.0) == (gridsize-0.5))) chains(m,l)%bond = 2
           if((chains(m,l)%y-(bondy/4.0) == 0.5) .or. (chains(m,l)%y-(bondy/4.0) == -(gridsize-0.5))) chains(m,l)%bond = -2
           if((chains(m,l)%z-(bondz/4.0) == -0.5) .or. (chains(m,l)%z-(bondz/4.0) == (gridsize-0.5))) chains(m,l)%bond = 3
           if((chains(m,l)%z-(bondz/4.0) == 0.5) .or. (chains(m,l)%z-(bondz/4.0) == -(gridsize-0.5))) chains(m,l)%bond = -3
           !write(6,*) 'bond',chains(m,l)%bond
           if(chains(m,l)%bond ==0) write(6,*) 'no bond allocation',m,l
        end do
     end do
     write(6,*) 'finished'





     call clustercount
!call neighbour
  end do


  norma =  timeofinterest !sum(summing(1,:))
  normb = timeofinterest !sum(summing(2,:))
  normc = timeofinterest !sum(summing(3,:))
  normd =  timeofinterest !sum(summing(4,:))


  dr = 1.0/binsize
  rho = (totpoints/2)/real(gridsize**3)
write(29,*) 0,0,0,0,0
write(38,*) 0,0,0,0,0,0,0,0,0,0
write(68,*) 0,0,0,0,0,0,0,0

do l = 1,((gridsize/2)*binsize)+1,1
   do x =1,4*mtype,1
sum1(x) = 0.0d0
do kx = 1,l,1
sum1(x) = sum1(x) + radial(x,kx)
end do
end do


sumall = 0.0d0
do kx = 1,l,1
sumall = sumall + allsum(1,kx) + allsum(2,kx)
end do



     r = (real(l)/binsize) !-(1.0/(2*real(binsize)))  
     write(29,*)r,summing(1,l)/(norma*4*PI_8*(r**2)*rho*dr),summing(2,l)/(normb*4*PI_8*(r**2)*rho*dr)&
          ,summing(3,l)/(normc*4*PI_8*(r**2)*rho*dr),summing(4,l)/(normd*4*PI_8*(r**2)*rho*dr),&
          summing(5,l)/(norma*4*PI_8*(r**2)*rho*dr),summing(6,l)/(normb*4*PI_8*(r**2)*rho*dr)&
          ,summing(7,l)/(normc*4*PI_8*(r**2)*rho*dr),summing(8,l)/(normd*4*PI_8*(r**2)*rho*dr)
!write(31,*) r,sum1/norma,sum2/norma,sum3/norma,sum4/norma
     write(31,*) r,radial(1,l)/norma,radial(2,l)/normb
!summing(1,l)/norma,summing(2,l)/normb,summing(3,l)/normc,summing(4,l)/normd
     
        r = (real(l)/binsize)           
        write(38,*) r,sum1(1)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),sum1(2)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),&
          sum1(3)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),sum1(4)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),&
          sum1(5)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),sum1(6)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),&
          sum1(7)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),sum1(8)/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),&
        (sum1(3)+sum1(1))/(timeofinterest*(r**3)*(4.0/3.0)*PI_8),(sum1(4)+sum1(2))/(timeofinterest*(r**3)*(4.0/3.0)*PI_8)


write(68,*) r,radial(:,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)) ,&
radial(2,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)),&
radial(3,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)),&
radial(4,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0))

!write(168,*)r,allsum(:,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)) !,&
!allsum(2,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0))

!write(168,*)r,allsum(:,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)) ,&
! (allsum(1,l)+allsum(2,l))/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0))

write(168,*)r,allsum(:,l)/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)),&
 (allsum(1,l)+allsum(2,l))/(timeofinterest*((r**3)-((real(l-1)/binsize)**3))*PI_8*(4.0/3.0)),&
 sumall/(timeofinterest*(r**3)*PI_8*(4.0/3.0))


  end do

write(6,*) 'average unbound BULK and CLUSTER:',(real(totunb(:))/timeofinterest)-1
write(6,*) 'cluster neighbours=',clusneigh(:),real(clusneigh(1))/(totunb(1)-1),real(clusneigh(2))/(totunb(2)-1)
write(6,*) 'bulkneighbours=',bulkneigh(:),real(bulkneigh(1))/(totunb(3)-1),real(bulkneigh(2))/(totunb(4)-1)
write(6,*)'ratios=',(real(clusneigh(1))/(totunb(1)-1))/(real(bulkneigh(1))/(totunb(3)-1)) &
,(real(clusneigh(2))/(totunb(2)-1))/(real(bulkneigh(2))/(totunb(4)-1))
write(6,*) 'FINAL cluster neighbours:',real(totclhist)/timeofinterest
write(6,*) 'FINAL bulk neighbours:',real(totblhist)/timeofinterest

write(712,*) int((real(totunb(2))/timeofinterest)-1),int((real(totunb(4))/timeofinterest)-1),clusneigh(2),bulkneigh(2)


contains

subroutine neighbour(clnos,bound,content)
integer,intent(in)::content
integer,dimension(:),intent(in)::clnos
integer,dimension(:,:),intent(in)::bound
integer::m,l,maxl,xp,yp,zp
integer,dimension(gridsize,gridsize,gridsize)::cpoints,bpoints
integer,dimension(6):: clhist,blhist
cpoints(:,:,:) = 0
bpoints(:,:,:) = 0
clhist(:) = 0
blhist(:) = 0

do m = 1,nprotein,1
maxl = chlen(m)
if(chains(m,1)%type ==2) then
if(clnos(m) /= content) then
do l = 1,maxl
if (bound(m,l) == 0) then
if(chains(m,l)%a1 == 0)then
 bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(modulo(int(chains(m,l)%x)+1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) = &
bpoints(modulo(int(chains(m,l)%x)+1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) + 1
end if
if(chains(m,l)%a2 == 0) then
bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(modulo(int(chains(m,l)%x)-1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) = &
bpoints(modulo(int(chains(m,l)%x)-1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) + 1
end if
if(chains(m,l)%a3 == 0) then
bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(chains(m,l)%x,modulo(int(chains(m,l)%y)+1-1,gridsize)+1,chains(m,l)%z) = &
bpoints(chains(m,l)%x,modulo(int(chains(m,l)%y)+1-1,gridsize)+1,chains(m,l)%z) + 1
end if
if(chains(m,l)%a4 == 0) then
bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(chains(m,l)%x,modulo(int(chains(m,l)%y)-1-1,gridsize)+1,chains(m,l)%z) = &
bpoints(chains(m,l)%x,modulo(int(chains(m,l)%y)-1-1,gridsize)+1,chains(m,l)%z) + 1
end if
if(chains(m,l)%a5 == 0) then
bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(chains(m,l)%x,chains(m,l)%y,modulo(int(chains(m,l)%z)+1-1,gridsize)+1) = &
bpoints(chains(m,l)%x,chains(m,l)%y,modulo(int(chains(m,l)%z)+1-1,gridsize)+1) + 1
end if
if(chains(m,l)%a6 == 0) then
bulkneigh(chains(m,l)%type) = bulkneigh(chains(m,l)%type)+1
bpoints(chains(m,l)%x,chains(m,l)%y,modulo(int(chains(m,l)%z)-1-1,gridsize)+1) = &
bpoints(chains(m,l)%x,chains(m,l)%y,modulo(int(chains(m,l)%z)-1-1,gridsize)+1) + 1
end if
end if
end do
else if(clnos(m) == content) then
do l = 1,maxl
if(bound(m,l) == 0) then
if(chains(m,l)%a1 == 0)then
 clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(modulo(chains(m,l)%x+1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) = &
cpoints(modulo(chains(m,l)%x+1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) + 1
end if
if(chains(m,l)%a2 == 0) then
clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(modulo(chains(m,l)%x-1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) = &
cpoints(modulo(chains(m,l)%x-1-1,gridsize)+1,chains(m,l)%y,chains(m,l)%z) + 1
end if
if(chains(m,l)%a3 == 0) then
clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(chains(m,l)%x,modulo(chains(m,l)%y+1-1,gridsize)+1,chains(m,l)%z) = &
cpoints(chains(m,l)%x,modulo(chains(m,l)%y+1-1,gridsize)+1,chains(m,l)%z) + 1
end if
if(chains(m,l)%a4 == 0) then
clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(chains(m,l)%x,modulo(chains(m,l)%y-1-1,gridsize)+1,chains(m,l)%z) = &
cpoints(chains(m,l)%x,modulo(chains(m,l)%y-1-1,gridsize)+1,chains(m,l)%z) + 1
end if
if(chains(m,l)%a5 == 0) then
clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(chains(m,l)%x,chains(m,l)%y,modulo(chains(m,l)%z+1-1,gridsize)+1) = &
cpoints(chains(m,l)%x,chains(m,l)%y,modulo(chains(m,l)%z+1-1,gridsize)+1) + 1
end if
if(chains(m,l)%a6 == 0) then
clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
cpoints(chains(m,l)%x,chains(m,l)%y,modulo(chains(m,l)%z-1-1,gridsize)+1) = &
cpoints(chains(m,l)%x,chains(m,l)%y,modulo(chains(m,l)%z-1-1,gridsize)+1) + 1
end if




!if(chains(m,l)%a6 == 0) clusneigh(chains(m,l)%type) = clusneigh(chains(m,l)%type)+1
end if
end do
end if

end if
end do


do xp = 1,gridsize
do yp = 1,gridsize
do zp = 1,gridsize
if(cpoints(xp,yp,zp)== 1) clhist(1) = clhist(1)+1
if(cpoints(xp,yp,zp)== 2) clhist(2) = clhist(2)+1
if(cpoints(xp,yp,zp)== 3) clhist(3) = clhist(3)+1
if(cpoints(xp,yp,zp)== 4) clhist(4) = clhist(4)+1
if(cpoints(xp,yp,zp)== 5) clhist(5) = clhist(5)+1
if(cpoints(xp,yp,zp)== 6) clhist(6) = clhist(6)+1


if(bpoints(xp,yp,zp)== 1) blhist(1) = blhist(1)+1
if(bpoints(xp,yp,zp)== 2) blhist(2) = blhist(2)+1
if(bpoints(xp,yp,zp)== 3) blhist(3) = blhist(3)+1
if(bpoints(xp,yp,zp)== 4) blhist(4) = blhist(4)+1
if(bpoints(xp,yp,zp)== 5) blhist(5) = blhist(5)+1
if(bpoints(xp,yp,zp)== 6) blhist(6) = blhist(6)+1


end do
end do
end do

write(6,*) 'Cluster neighbours:',clhist(:)
write(6,*) 'Bulk neighbours',blhist(:)

do l = 1,6
totclhist(l) = totclhist(l) + clhist(l)
totblhist(l) = totblhist(l) + blhist(l)
end do
end subroutine neighbour


  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,content,onepop,z
    integer :: acconn,cconn,atconn,tconn,a,b,c,d,a2,b2,c2,d2
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(rprot),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer ::oldcl,dum1,dum2,zzz,maxclus,normalisedsize,dummy,unbcount,bcount,backupcount
    integer,dimension(:,:),allocatable::cb,conn
    integer,dimension(:),allocatable::maxcluslist,concount
    integer,dimension(:,:),allocatable::bound
    double precision::maxclusenergy,initialenergy,comx,comy,comz
         type(rprot) :: dpcomcluster



    allocate(cb(nprotein,nprotein))
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(conn(nprotein,maxlength))
    allocate(concount(nprotein))
    allocate(bound(nprotein,maxlength))
    !allocate(mashist(10))
    maxclusenergy = 0.0d0
    initialenergy = 0.0
    normalisedsize = 0
    clcount = 0
    clusterpop = 0
    clustcount = 0

runningclusterlist(:,1) = 0

do m= 1,nprotein,1
maxl= chlen(m)
do l =1,maxl
chains(m,l)%a1 = 0
chains(m,l)%a2 = 0
chains(m,l)%a3 = 0
chains(m,l)%a4 = 0
chains(m,l)%a5 = 0
chains(m,l)%a6 = 0


end do
end do
    do m = 1,nprotein
       clnos(m) = m
       cllist(m) = 0
       histcl(m) = 0
       concount(m) =0
       do a =1,maxlength
          bound(m,a) = 0
        conn(m,a) = 0
       end do
    

end do

    a = 1

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do


    do m = 1,nprotein-1,1
       maxback = chlen(m)
       do g = m+1,nprotein,1
         ! if((g/=m) .and. (clnos(m) /= clnos(g))) then
             maxlengthss = chlen(g)
             do l = 1,maxback
                tempcoord(l)%x = chains(m,l)%x
                tempcoord(l)%y =  chains(m,l)%y
                tempcoord(l)%z = chains(m,l)%z

                do f = l+1,maxback
                   call adjacent(m,l,f,tempcoord,adjver,bdir)              
                   if(adjver.eqv. .true.) then
                      if(bdir ==-1) then
                         chains(m,l)%a1 = 1
                         chains(m,f)%a2 = 1
                      else if(bdir ==1) then
                         chains(m,l)%a2 = 1
                         chains(m,f)%a1 = 1
                      else if(bdir ==-2) then
                         chains(m,l)%a3 = 1
                         chains(m,f)%a4 = 1
                      else if(bdir ==2) then
                         chains(m,l)%a4 = 1
                         chains(m,f)%a3 = 1
                      else if(bdir ==-3) then
                         chains(m,l)%a5 = 1
                         chains(m,f)%a6 = 1
                      else if(bdir ==3) then
                         chains(m,l)%a6 = 1
                         chains(m,f)%a5 = 1
                      end if
                      end if
                   end do
                   
                do f = 1,maxlengthss
                   call adjacent(g,l,f,tempcoord,adjver,bdir)              
                   if(adjver.eqv. .true.) then
                      if(bdir ==-1) then
                         chains(m,l)%a1 = 1
                         chains(g,f)%a2 = 1
                      else if(bdir ==1) then
                         chains(m,l)%a2 = 1
                         chains(g,f)%a1 = 1
                      else if(bdir ==-2) then
                         chains(m,l)%a3 = 1
                         chains(g,f)%a4 = 1
                      else if(bdir ==2) then
                         chains(m,l)%a4 = 1
                         chains(g,f)%a3 = 1
                      else if(bdir ==-3) then
                         chains(m,l)%a5 = 1
                         chains(g,f)%a6 = 1
                      else if(bdir ==3) then
                         chains(m,l)%a6 = 1
                         chains(g,f)%a5 = 1
                      end if
        
                if(isbond .eqv. .true.) then
                         cb(m,g) = 1
                         cb(g,m) = 1
                         !bound(m,l) = 1
                         !bound(g,f) = 1
                         if(interen(chains(m,l)%type,chains(g,f)%type) < 0.0) then
                     


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




       initialenergy = initialenergy + interen(chains(m,l)%type,chains(g,f)%type)    
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
                      if((bdir == chains(g,f)%bond) .and. (bdir == (-1*chains(m,l)%bond)) .and. &
                           (interenergy(chains(g,f)%type,chains(m,l)%type)  < 0.0)) then
                         initialenergy = initialenergy + interenergy(chains(m,l)%type,chains(g,f)%type) 



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


                         !write(6,*) 'energy update',interenergy(chains(m,l)%type,chains(g,f)%type),initialenergy
!if(bound(g,f) ==1) write(6,*) 'already bound!!!!!!!!!!!!!!!'
                         bound(m,l) = 1
                         bound(g,f) = 1

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
          !end if
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
    atconn= 0
 backupcount = 0
   do m = 1,nprotein
      if(chains(m,1)%type==1) atconn = atconn+concount(m)
       if(clnos(m) == content)then
          maxcluslist(z) = m
        runningclusterlist(z,1) = m
          z = z+1
backupcount = backupcount+1
       end if
    end do

   cconn =0
    acconn =0

    do m = 1,maxclus
       g = maxcluslist(m)
       cconn = cconn + concount(g)
       if(chains(g,1)%type ==1) acconn = acconn+concount(g)
       !maxclusenergy = maxclusenergy + chen(g)
    end do

    tconn=sum(concount(:))

   onepop = 0
    do m = 1,maxclus
       g = maxcluslist(m)
       if(chains(g,1)%type ==1) onepop = onepop+1
    end do



write(6,*) 'connectivity',acconn,cconn-acconn,atconn-acconn,tconn-((atconn-acconn)+cconn)
write(6,*) 'CONNECTIVITY RATIOS',real(acconn)/onepop,real(cconn-acconn)/(maxclus-onepop),&
real(atconn-acconn)/(speciespop(1) - onepop),real(tconn-((atconn-acconn)+cconn))/(speciespop(2)-(maxclus-onepop))
Write(6,*) 'CLUSTER SETUP',onepop,maxclus-onepop,speciespop(1) - onepop,speciespop(2) - (maxclus-onepop)
write(211,*) t,real(acconn)/onepop,real(cconn-acconn)/(maxclus-onepop),&
real(atconn-acconn)/(speciespop(1) - onepop),real(tconn-((atconn-acconn)+cconn))/(speciespop(2)-(maxclus-onepop))

call compare(backupcount)

call neighbour(clnos,bound,content)


    a = 1
    b = 1
    c = 1
    d = 1
    a2 = 1
    b2 = 1
    c2 = 1
    d2 = 1

    do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          if(bound(m,l) == 0) then
        if((chains(m,l)%a1 == 1) .and. (chains(m,l)%a2 == 1 ).and.  (chains(m,l)%a3 == 1) .and. &
        (chains(m,l)%a4 == 1) .and. (chains(m,l)%a5 == 1) .and. (chains(m,l)%a6 == 1)) then
                if(clnos(m) == content)then    
                        if(chains(m,l)%type ==1) then
                                a2clus(a2)%m = m
                                a2clus(a2)%l = l
                                a2 = a2+1
                        else if(chains(m,l)%type ==2) then
                                b2clus(b2)%m = m
                                b2clus(b2)%l = l
                                b2 = b2+1
                        end if
                else if(clnos(m) /= content) then
                        if(chains(m,l)%type ==1) then
                                a2bulk(c2)%m = m
                                a2bulk(c2)%l = l
                                c2 = c2+1
                        else if(chains(m,l)%type ==2) then
                                b2bulk(d2)%m = m
                                b2bulk(d2)%l = l
                                d2 = d2+1
                end if
             end if




   
   goto 119
else

             if(clnos(m) == content)then    
            if(chains(m,l)%type ==1) then
                   aclus(a)%m = m
                   aclus(a)%l = l
                   a = a+1
                else if(chains(m,l)%type ==2) then
                   bclus(b)%m = m
                   bclus(b)%l = l
                   b = b+1
                end if
             else if(clnos(m) /= content) then
                if(chains(m,l)%type ==1) then
                   abulk(c)%m = m
                   abulk(c)%l = l
                   c = c+1
                else if(chains(m,l)%type ==2) then
                   bbulk(d)%m = m
                   bbulk(d)%l = l
                   d = d+1
                end if
             end if
end if
119 continue
          end if
       end do
    end do
    write(6,*) 'aclus',a,'bclus',b,'abulk',c,'bbulk',d,en(t),maxclus
    write(6,*) 'a2clus',a2,'b2clus',b2,'a2bulk',c2,'b2bulk',d2,en(t),maxclus
    write(196,*) a,b,c,d,en(t),maxclus
   !write(96,*) 'FAIL' 

bcount = 0
    unbcount =0
      do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          if(bound(m,l)==1) bcount =bcount+1
          if(bound(m,l)==0) unbcount = unbcount+1

       end do
       end do

       dummy = a+b+c+d-4
       write(6,*) 'bound =',bcount,'unbound=',unbcount,'total=',dummy

totunb(1) = totunb(1) + a
totunb(2) = totunb(2) + b
totunb(3) = totunb(3) + c
totunb(4) = totunb(4) + d

do m =1,nprotein
call comfind(m)
end do

write(6,*) 'postcom time',t


    call phase(content,maxclus,clnos,cb,dpcomcluster)


comx = dpcomcluster%x
comy = dpcomcluster%y
comz = dpcomcluster%z

!comx = readcomx
!comy = readcomy
!comz = readcom
write(6,*) 'COM values', comx,comy,comz
    call paircorrelation(a,b,c,d,a2,b2,c2,d2,comx,comy,comz)
!write(6,*) 'energy',en(1),timeofinterest

call rdf(comx,comy,comz)

    if (nint(initialenergy) /= nint(en(t))) then
       write(6,*) 'energyfail!!!!!!!!!!!!!!',t,initialenergy,en(t)
    end if



    deallocate(cllist)
    deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
  end subroutine clustercount

subroutine compare(backupcount)
integer,intent(inout)::backupcount
integer :: a,leave,join


leave = 0
join = 0
do a = 1,nprotein,1
if(ANY(runningclusterlist(:,1) .eq. a)) then
        if(ANY(runningclusterlist(:,2) .eq. a)) then
        continue
        else
        join = join + 1
        end if
end if

if(ANY(runningclusterlist(:,2) .eq. a)) then
        if(ANY(runningclusterlist(:,1) .eq. a)) then
        continue
        else
        leave = leave + 1
        end if
end if

end do

open(99,file='stability.dat',access = 'append')
write(99,*) join,leave,join-leave
close(99)

runningclusterlist(:,2) = 0

do a = 1,backupcount,1
runningclusterlist(a,2) = runningclusterlist(a,1)
end do






end subroutine compare



 subroutine comfind(m)
    integer :: dcomx,dcomy,dcomz,maxl,comx,comy,comz
    integer :: l,a
    integer,intent(in)::m



       comx = 0
       comy = 0
       comz = 0
        maxl = chlen(m)
!linkerlength = 7

        !do a =1,maxl
           !write(6,*)
           !'bead',a,chains(m,a)%x,chains(m,a)%y,chains(m,a)%z
        !end do
        do l = 2,maxl,1
           dcomx = chains(m,l)%x-chains(m,l-1)%x
           if (dcomx >chains(m,1)%linker) dcomx = dcomx-gridsize
           if (dcomx <(-chains(m,1)%linker)) dcomx = gridsize+dcomx
           comx = comx + ((maxl+1-l)*dcomx)
           dcomy = chains(m,l)%y-chains(m,l-1)%y
           !if(abs(dcomy) > 7) write(6,*) 'over pbc',gridsize,linkerlength,dcomy
           if (dcomy >chains(m,1)%linker) dcomy = dcomy-gridsize
           if (dcomy <(-chains(m,1)%linker)) dcomy = gridsize+dcomy
           !if(abs(dcomy) > 7) write(6,*)
           !'faillllll',gridsize,linkerlength,dcomy
           comy = comy + ((maxl+1-l)*dcomy)
           dcomz = chains(m,l)%z-chains(m,l-1)%z
           if (dcomz >chains(m,1)%linker) dcomz = dcomz-gridsize
           if (dcomz <(-chains(m,1)%linker)) dcomz = gridsize+dcomz
           comz = comz + ((maxl+1-l)*dcomz)
        end do
        !write(6,*) 'deltas',comx,comy,comz
        com(1,m)%x = modulo(chains(m,1)%x+(real(comx)/(maxl))-1,real(gridsize))+1
        com(1,m)%y = modulo(chains(m,1)%y+(real(comy)/(maxl))-1,real(gridsize))+1
        com(1,m)%z = modulo(chains(m,1)%z+(real(comz)/(maxl))-1,real(gridsize))+1
        !end do
        !if((com(1,1)%y>1500) .and.(com(1,1)%y<2500))write(6,*)
        !'com',time,com(1,1)%y,chains(1,1)%y,comy
        !write(6,*) 'commmm 2'

  end subroutine comfind


  subroutine adjacent(pr,ch1,ch2,tempory,adjver,bdir)
    integer,intent(in):: pr,ch1,ch2
    integer,intent(inout)::bdir
    integer::dx,dy,dz,delx,dely,delz,signs
    logical,intent(inout)::adjver
    Type(rprot),dimension(:),intent(in) :: tempory
    dx = 4
    dy = 4
    dz = 4

    delx = tempory(ch1)%x - chains(pr,ch2)%x
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

    dely = tempory(ch1)%y - chains(pr,ch2)%y
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

    delz = tempory(ch1)%z - chains(pr,ch2)%z
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

subroutine phase(mz,clsize,clnos,cb,dpcomcluster)
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
    type(rprot),intent(inout) :: dpcomcluster
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
write(6,*) 'phase 1 time',t
    call clustercom(mz,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge,dpcomcluster)
write(6,*) 'phase2 time',t
  end subroutine phase



subroutine rdf(comx,comy,comz)
double precision,intent(in)::comx,comy,comz
integer :: m,l,maxl,lx
    double precision :: dx,dy,dz
    double precision :: delta



    do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          dx = min(abs(chains(m,l)%x - comx), gridsize -abs(chains(m,l)%x - comx))
          dy = min(abs(chains(m,l)%y - comy), gridsize -abs(chains(m,l)%y - comy))
          dz = min(abs(chains(m,l)%z - comz), gridsize -abs(chains(m,l)%z - comz))

   delta = sqrt((dx**2)+(dy**2)+(dz**2))

   lx = INT(delta/(1.0/binsize))


            allsum(chains(m,l)%type,lx+1) =allsum(chains(m,l)%type,lx+1)+1.0
       end do
       end do



end subroutine rdf



  subroutine clustercom(m,clsize,comcluster,cb,cllist,route,rcounter,totcluspop,checklist,toolarge,dpcomcluster)
          integer,intent(in)::m
          integer,intent(inout)::clsize,totcluspop
          integer,dimension(:),intent(inout) ::cllist,rcounter,checklist
          integer,dimension(:,:),intent(inout) ::cb,route
          integer ::a,b,base,maxl,maxb,ra,dint,centralchain,f,g,cench,rc,zx,maxltest
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

          !t = time+1

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
             !if(clsize>1) write(6,*)
             !'component',b,com(1,b)%x,com(1,b)%y,com(1,b)%z
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
16              if(a == a) continue
             end do

          end do


          !write(6,*) 'backup 2',rcounter(4)
          !write(6,*) 'orig central chain',centralchain,toolarge
          call comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
          !write(6,*) 'key'


          !write(6,*) 'backup 3',rcounter(4)

76        if(clsize == 1) then
             maxl = chlen(m)
             centralchain = m
             totcluspop = maxl
             ! write(6,*) 'individual chain',com(1,m)%x,com(1,m)%y,com(1,m)%z
          end if

          !if(clsize>1) write(6,*) 'totcluspop',totcluspop
          !write(6,*) centralchain
          !write(6,*) 'runtime', com(1,centralchain)%x
          !write(6,*) dummycom%x
          comcluster%x = INT(modulo((com(1,centralchain)%x)+(dummycom%x/totcluspop)-1,real(gridsize))+1)
          comcluster%y = INT(modulo((com(1,centralchain)%y)+(dummycom%y/totcluspop)-1,real(gridsize))+1)
          comcluster%z = INT(modulo((com(1,centralchain)%z)+(dummycom%z/totcluspop)-1,real(gridsize))+1)


          dpcomcluster%x = modulo((com(1,centralchain)%x)+(dummycom%x/totcluspop)-1,real(gridsize))+1
          dpcomcluster%y = modulo((com(1,centralchain)%y)+(dummycom%y/totcluspop)-1,real(gridsize))+1
          dpcomcluster%z = modulo((com(1,centralchain)%z)+(dummycom%z/totcluspop)-1,real(gridsize))+1
          write(3,*) time,dpcomcluster%x,dpcomcluster%y,dpcomcluster%z


if(int(dpcomcluster%x) /= int(readcomx)) write(6,*)'FAILLLLLL',dpcomcluster%x,readcomx
if(int(dpcomcluster%y) /= int(readcomy))write(6,*)'FAILLLLLL',dpcomcluster%y,readcomy
if(int(dpcomcluster%z) /= int(readcomz))write(6,*)'FAILLLLLL',dpcomcluster%z,readcomz




53        if(compass .eqv. .true.) continue

        end subroutine clustercom


        subroutine comlargecomp(clsize,dummycom,checklist,route,rcounter,comx,comy,comz,toolarge)
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
             !write(6,*) 'finalcoms',
             !comx(centralchain,route(zna,1)),rcounter(zna)

             bx = comx(checklist(1),route(a,1))
             by =comy(checklist(1),route(a,1))
             bz = comz(checklist(1),route(a,1))

             ! write(6,*) 'dummycom', dummycom%x,dummycom%y,dummycom%z
             !write(6,*)
             !'x',comx(checklist(1),route(a,1)),checklist(1),route(a,1),&
             !  com(1,checklist(1))%x,com(1,route(a,1))%x
             ! write(6,*)
             ! 'y',comy(checklist(1),route(a,1)),checklist(1),route(a,1),&
             !      com(1,checklist(1))%y,com(1,route(a,1))%y
             !write(6,*)
             !'z',comz(checklist(1),route(a,1)),checklist(1),route(a,1),&
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


  subroutine paircorrelation(a,b,c,d,a2,b2,c2,d2,comx,comy,comz)
    integer,intent(in)::a,b,c,d,a2,b2,c2,d2
    integer::k,h
    double precision::dx,dy,dz
double precision, intent(in)::comx,comy,comz







    if(a>1) then

       do h = 1,a-1
             dx = min(abs(comx-chains(aclus(h)%m,aclus(h)%l)%x),gridsize-&
                  abs(comx-chains(aclus(h)%m,aclus(h)%l)%x))
             dy = min(abs(comy-chains(aclus(h)%m,aclus(h)%l)%y),gridsize-&
                  abs(comy-chains(aclus(h)%m,aclus(h)%l)%y))
             dz = min(abs(comz-chains(aclus(h)%m,aclus(h)%l)%z),gridsize-&
                  abs(comz-chains(aclus(h)%m,aclus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(1,lx+1) = radial(1,lx+1)+(1.0)

             
          end do
       
       do h = 1,a-2
          do k = h+1,a-1

             if(k==h) goto 75

             dx = min(abs(chains(aclus(k)%m,aclus(k)%l)%x-chains(aclus(h)%m,aclus(h)%l)%x),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%x-chains(aclus(h)%m,aclus(h)%l)%x))
             dy = min(abs(chains(aclus(k)%m,aclus(k)%l)%y-chains(aclus(h)%m,aclus(h)%l)%y),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%y-chains(aclus(h)%m,aclus(h)%l)%y))
             dz = min(abs(chains(aclus(k)%m,aclus(k)%l)%z-chains(aclus(h)%m,aclus(h)%l)%z),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%z-chains(aclus(h)%m,aclus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             if(abs(lx)>(gridsize*binsize)) goto 75

             summing(1,lx+1) = summing(1,lx+1)+(1.0) !/((a-1)**2))

75           continue
          end do
       end do
    end if

    if(b>1) then


       do h = 1,b-1
             dx = min(abs(comx-chains(bclus(h)%m,bclus(h)%l)%x),gridsize-&
                  abs(comx-chains(bclus(h)%m,bclus(h)%l)%x))
             dy = min(abs(comy-chains(bclus(h)%m,bclus(h)%l)%y),gridsize-&
                  abs(comy-chains(bclus(h)%m,bclus(h)%l)%y))
             dz = min(abs(comz-chains(bclus(h)%m,bclus(h)%l)%z),gridsize-&
                  abs(comz-chains(bclus(h)%m,bclus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(2,lx+1) = radial(2,lx+1)+(1.0)

             
          end do
       
       do h = 1,b-2
          do k = h+1,b-1

             if(k==h) goto 95
             !write(6,*) bclus(k)%l,bclus(h)%l,h,k,b
             dx = min(abs(chains(bclus(k)%m,bclus(k)%l)%x-chains(bclus(h)%m,bclus(h)%l)%x),gridsize-&
                  abs(chains(bclus(k)%m,bclus(k)%l)%x-chains(bclus(h)%m,bclus(h)%l)%x))
             dy = min(abs(chains(bclus(k)%m,bclus(k)%l)%y-chains(bclus(h)%m,bclus(h)%l)%y),gridsize-&
                  abs(chains(bclus(k)%m,bclus(k)%l)%y-chains(bclus(h)%m,bclus(h)%l)%y))
             dz = min(abs(chains(bclus(k)%m,bclus(k)%l)%z-chains(bclus(h)%m,bclus(h)%l)%z),gridsize-&
                  abs(chains(bclus(k)%m,bclus(k)%l)%z-chains(bclus(h)%m,bclus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             if(abs(lx+1)>(gridsize*binsize)) goto 95

             summing(2,lx+1) = summing(2,lx+1)+(1.0) !/((b-1)**2))

95           continue
          end do
       end do
    end if

    if(c>1) then

       do h = 1,c-1
             dx = min(abs(comx-chains(abulk(h)%m,abulk(h)%l)%x),gridsize-&
                  abs(comx-chains(abulk(h)%m,abulk(h)%l)%x))
             dy = min(abs(comy-chains(abulk(h)%m,abulk(h)%l)%y),gridsize-&
                  abs(comy-chains(abulk(h)%m,abulk(h)%l)%y))
             dz = min(abs(comz-chains(abulk(h)%m,abulk(h)%l)%z),gridsize-&
                  abs(comz-chains(abulk(h)%m,abulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175
!if(lx>((gridsize/2)+1)) then
!write(6,*) 'too big',lx,dx,dy,dz,delta
!end if
             radial(3,lx+1) = radial(3,lx+1)+(1.0)

             
          end do

       
       do h = 1,c-2
          do k = h+1,c-1

             if(k==h) goto 45

             dx = min(abs(chains(abulk(k)%m,abulk(k)%l)%x-chains(abulk(h)%m,abulk(h)%l)%x),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%x-chains(abulk(h)%m,abulk(h)%l)%x))
             dy = min(abs(chains(abulk(k)%m,abulk(k)%l)%y-chains(abulk(h)%m,abulk(h)%l)%y),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%y-chains(abulk(h)%m,abulk(h)%l)%y))
             dz = min(abs(chains(abulk(k)%m,abulk(k)%l)%z-chains(abulk(h)%m,abulk(h)%l)%z),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%z-chains(abulk(h)%m,abulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
            ! if(lx ==0) lx =1
             if(abs(lx)>(gridsize*binsize)) goto 45

             summing(3,lx+1) = summing(3,lx+1)+(1.0) !/((c-1)**2))

45           continue
          end do
       end do
    end if

    if(d>1) then

         do h = 1,d-1
            dx = min(abs(comx-chains(bbulk(h)%m,bbulk(h)%l)%x),gridsize-&
                  abs(comx-chains(bbulk(h)%m,bbulk(h)%l)%x))
             dy = min(abs(comy-chains(bbulk(h)%m,bbulk(h)%l)%y),gridsize-&
                  abs(comy-chains(bbulk(h)%m,bbulk(h)%l)%y))
             dz = min(abs(comz-chains(bbulk(h)%m,bbulk(h)%l)%z),gridsize-&
                  abs(comz-chains(bbulk(h)%m,bbulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(4,lx+1) = radial(4,lx+1)+(1.0)

             
          end do
       
       do h = 1,d-2
          do k = h+1,d-1

             if(k==h) goto 55

             dx = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%x-chains(bbulk(h)%m,bbulk(h)%l)%x),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%x-chains(bbulk(h)%m,bbulk(h)%l)%x))
             dy = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%y-chains(bbulk(h)%m,bbulk(h)%l)%y),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%y-chains(bbulk(h)%m,bbulk(h)%l)%y))
             dz = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%z-chains(bbulk(h)%m,bbulk(h)%l)%z),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%z-chains(bbulk(h)%m,bbulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             
             if(abs(lx)>(gridsize*binsize))     goto 55
             
             summing(4,lx+1) = summing(4,lx+1)+(1.0) !/((d-1)**2))

55           continue
          end do
       end do
    end if



   if(a2>1) then

       do h = 1,a2-1
             dx = min(abs(comx-chains(a2clus(h)%m,a2clus(h)%l)%x),gridsize-&
                  abs(comx-chains(a2clus(h)%m,a2clus(h)%l)%x))
             dy = min(abs(comy-chains(a2clus(h)%m,a2clus(h)%l)%y),gridsize-&
                  abs(comy-chains(a2clus(h)%m,a2clus(h)%l)%y))
             dz = min(abs(comz-chains(a2clus(h)%m,a2clus(h)%l)%z),gridsize-&
                  abs(comz-chains(a2clus(h)%m,a2clus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(5,lx+1) = radial(5,lx+1)+(1.0)

             
          end do
       
       do h = 1,a2-2
          do k = h+1,a2-1


             dx = min(abs(chains(a2clus(k)%m,a2clus(k)%l)%x-chains(a2clus(h)%m,a2clus(h)%l)%x),gridsize-&
                  abs(chains(a2clus(k)%m,a2clus(k)%l)%x-chains(a2clus(h)%m,a2clus(h)%l)%x))
             dy = min(abs(chains(a2clus(k)%m,a2clus(k)%l)%y-chains(a2clus(h)%m,a2clus(h)%l)%y),gridsize-&
                  abs(chains(a2clus(k)%m,a2clus(k)%l)%y-chains(a2clus(h)%m,a2clus(h)%l)%y))
             dz = min(abs(chains(a2clus(k)%m,a2clus(k)%l)%z-chains(a2clus(h)%m,a2clus(h)%l)%z),gridsize-&
                  abs(chains(a2clus(k)%m,a2clus(k)%l)%z-chains(a2clus(h)%m,a2clus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1


             summing(5,lx+1) = summing(5,lx+1)+(1.0) !/((a-1)**2))

          end do
       end do
    end if

    if(b2>1) then


       do h = 1,b2-1
             dx = min(abs(comx-chains(b2clus(h)%m,b2clus(h)%l)%x),gridsize-&
                  abs(comx-chains(b2clus(h)%m,b2clus(h)%l)%x))
             dy = min(abs(comy-chains(b2clus(h)%m,b2clus(h)%l)%y),gridsize-&
                  abs(comy-chains(b2clus(h)%m,b2clus(h)%l)%y))
             dz = min(abs(comz-chains(b2clus(h)%m,b2clus(h)%l)%z),gridsize-&
                  abs(comz-chains(b2clus(h)%m,b2clus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(6,lx+1) = radial(6,lx+1)+(1.0)

             
          end do
       
       do h = 1,b2-2
          do k = h+1,b2-1

             dx = min(abs(chains(b2clus(k)%m,b2clus(k)%l)%x-chains(b2clus(h)%m,b2clus(h)%l)%x),gridsize-&
                  abs(chains(b2clus(k)%m,b2clus(k)%l)%x-chains(b2clus(h)%m,b2clus(h)%l)%x))
             dy = min(abs(chains(b2clus(k)%m,b2clus(k)%l)%y-chains(b2clus(h)%m,b2clus(h)%l)%y),gridsize-&
                  abs(chains(b2clus(k)%m,b2clus(k)%l)%y-chains(b2clus(h)%m,b2clus(h)%l)%y))
             dz = min(abs(chains(b2clus(k)%m,b2clus(k)%l)%z-chains(b2clus(h)%m,b2clus(h)%l)%z),gridsize-&
                  abs(chains(b2clus(k)%m,b2clus(k)%l)%z-chains(b2clus(h)%m,b2clus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1


             summing(6,lx+1) = summing(6,lx+1)+(1.0) !/((b-1)**2))

             end do
       end do
    end if

    if(c2>1) then

       do h = 1,c2-1
             dx = min(abs(comx-chains(a2bulk(h)%m,a2bulk(h)%l)%x),gridsize-&
                  abs(comx-chains(a2bulk(h)%m,a2bulk(h)%l)%x))
             dy = min(abs(comy-chains(a2bulk(h)%m,a2bulk(h)%l)%y),gridsize-&
                  abs(comy-chains(a2bulk(h)%m,a2bulk(h)%l)%y))
             dz = min(abs(comz-chains(a2bulk(h)%m,a2bulk(h)%l)%z),gridsize-&
                  abs(comz-chains(a2bulk(h)%m,a2bulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(7,lx+1) = radial(7,lx+1)+(1.0)

             
          end do

       
       do h = 1,c2-2
          do k = h+1,c2-1

    

             dx = min(abs(chains(a2bulk(k)%m,a2bulk(k)%l)%x-chains(a2bulk(h)%m,a2bulk(h)%l)%x),gridsize-&
                  abs(chains(a2bulk(k)%m,a2bulk(k)%l)%x-chains(a2bulk(h)%m,a2bulk(h)%l)%x))
             dy = min(abs(chains(a2bulk(k)%m,a2bulk(k)%l)%y-chains(a2bulk(h)%m,a2bulk(h)%l)%y),gridsize-&
                  abs(chains(a2bulk(k)%m,a2bulk(k)%l)%y-chains(a2bulk(h)%m,a2bulk(h)%l)%y))
             dz = min(abs(chains(a2bulk(k)%m,a2bulk(k)%l)%z-chains(a2bulk(h)%m,a2bulk(h)%l)%z),gridsize-&
                  abs(chains(a2bulk(k)%m,a2bulk(k)%l)%z-chains(a2bulk(h)%m,a2bulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1


             summing(7,lx+1) = summing(7,lx+1)+(1.0) !/((c-1)**2))

          end do
       end do
    end if

    if(d2>1) then

         do h = 1,d2-1
            dx = min(abs(comx-chains(b2bulk(h)%m,b2bulk(h)%l)%x),gridsize-&
                  abs(comx-chains(b2bulk(h)%m,b2bulk(h)%l)%x))
             dy = min(abs(comy-chains(b2bulk(h)%m,b2bulk(h)%l)%y),gridsize-&
                  abs(comy-chains(b2bulk(h)%m,b2bulk(h)%l)%y))
             dz = min(abs(comz-chains(b2bulk(h)%m,b2bulk(h)%l)%z),gridsize-&
                  abs(comz-chains(b2bulk(h)%m,b2bulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             !if(abs(lx)>(gridsize*binsize)) goto 175

             radial(8,lx+1) = radial(8,lx+1)+(1.0)

             
          end do
       
       do h = 1,d2-2
          do k = h+1,d2-1



             dx = min(abs(chains(b2bulk(k)%m,b2bulk(k)%l)%x-chains(b2bulk(h)%m,b2bulk(h)%l)%x),gridsize-&
                  abs(chains(b2bulk(k)%m,b2bulk(k)%l)%x-chains(b2bulk(h)%m,b2bulk(h)%l)%x))
             dy = min(abs(chains(b2bulk(k)%m,b2bulk(k)%l)%y-chains(b2bulk(h)%m,b2bulk(h)%l)%y),gridsize-&
                  abs(chains(b2bulk(k)%m,b2bulk(k)%l)%y-chains(b2bulk(h)%m,b2bulk(h)%l)%y))
             dz = min(abs(chains(b2bulk(k)%m,b2bulk(k)%l)%z-chains(b2bulk(h)%m,b2bulk(h)%l)%z),gridsize-&
                  abs(chains(b2bulk(k)%m,b2bulk(k)%l)%z-chains(b2bulk(h)%m,b2bulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             !if(lx ==0) lx =1
             

             
             summing(8,lx+1) = summing(8,lx+1)+(1.0) !/((d-1)**2))


          end do
       end do
    end if



    
  end subroutine paircorrelation




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
allocate(specieslink(nspecies))
       CASE ('TOTCHAINS')
          !CALL get_integer(nprotein)

          runtype = 1

       CASE ('NEWCHAIN')
          CALL get_integer(specieslen(runtype))
          CALL get_integer(speciespop(runtype))
          call get_integer(specieslink(runtype))   
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


          allocate(interen(mtype,mtype))
          allocate(interenergy(mtype,mtype))
          allocate(intraenergy(mtype,mtype))


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



