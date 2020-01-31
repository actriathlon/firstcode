program move

  implicit none

  character(len=10)::BIN
  integer,dimension(:),allocatable::type
  double precision,dimension(:),allocatable::coordx,coordy,coordz,en
  integer::totpoints,n,m,l,t,j,lx,timeofinterest,dummymax,count,e,maxlength
integer::totmonoc,totdic,totclusc,maxlink,totbulkc
  double precision::bondx,bondy,bondz,rx,ry,rz,rr,totmono,totdi,totclus,totbulk
  integer,dimension(:),allocatable :: chlen,speciespop,specieslen,linklen
  double precision::delta,r,dr,rho,intraen,totiso
  real,dimension(:,:),allocatable::summing
  integer::binsize,maxl,gridsize,unusedtime,maxtime,nprotein,mtype,nspecies,preinterest,time,ft,gt
  character(len = 10) :: commandread,commandread2,commandread3,commandread4
  character(len=10) :: u,v,w
  real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
  double precision,dimension(:,:),allocatable :: interenergy,interen,intraenergy
  logical ::isbond  
integer:: norma,normb,normc,normd

  type rprot
     real :: x,y,z
     integer::bond,type
  end type rprot

  type listt
     integer:: m,l
  end type listt

  type(listt),dimension(:),allocatable:: abulk,bbulk,aclus,bclus
  type(rprot),dimension(:),allocatable :: com
  type(rprot),dimension(:,:),allocatable :: chains



totmono=0
totdi =0
totclus = 0
totmonoc = 0
totdic = 0
totclusc = 0
totbulk = 0
totbulkc = 0
  open(21, file = 'movetagged.vtf', action = 'read')
  open(29, file = 'sepvals.dat', action = 'write')
  open(31, file = 'unnormRDF.dat', action = 'write')
  open(39, file='energy.dat',action='read')
open(49,file='histsep.dat',action='write')
  isbond = .true.
  call read_setup
  write(6,*) 'mytype',mtype,'int 0',int(0.2),int(1.2),int(-0.2),int(-1.2)

  do ft = 1,mtype
     do gt = 1,mtype
        totiso = totiso + interen(ft,gt)
        write(6,*) 'energies',ft,gt, 'dir',interenergy(ft,gt),'iso',interen(ft,gt)
     end do
  end do
  if(nint(totiso) == 0) isbond = .false.
!isbond = .false.
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


  write(6,*) 'inputs', binsize,unusedtime,timeofinterest,preinterest
  if(specieslen(1)>=specieslen(2)) maxlength = specieslen(1)
  if(specieslen(1)<specieslen(2)) maxlength = specieslen(2)

  allocate(summing(3,(maxlink*binsize)+1))
  allocate(com(totpoints/2))
  allocate(coordx(totpoints/2))
  allocate(coordy(totpoints/2))
  allocate(coordz(totpoints/2))
  allocate(chains(nprotein,maxlength))
  allocate(aclus(specieslen(1)*speciespop(1)))
  allocate(bclus(specieslen(2)*speciespop(2)))
  allocate(abulk(specieslen(1)*speciespop(1)))
  allocate(bbulk(specieslen(2)*speciespop(2)))
  allocate(en(timeofinterest))
  allocate(chlen(nprotein))
  allocate(type(totpoints/2))




  summing(:,:) = 0



  count = 1

  read(21,*) BIN,BIN,BIN,BIN,BIN,BIN
  do m = 1,nprotein
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,chains(m,1)%type
     read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN

     maxl = specieslen(chains(m,1)%type)
     chlen(m) = specieslen(chains(m,1)%type)
     count = count+1

     do l = 2,maxl
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,chains(m,l)%type
        read(21,*)BIN,BIN,BIN,BIN,BIN,BIN,BIN,BIN
        count = count+1
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
  read(39,*) BIN,BIN
   end do
write(6,*) 'preinterest =',preinterest
   
  do t = 1,unusedtime
     read(21,*) BIN
     read(21,*) BIN
     read(39,*) BIN,BIN
     !read(21,*) BIN,BIN,BIN
     do m = 1,totpoints
        read(21,*)BIN
     end do
  end do


  
  do t = 1,timeofinterest
time = t
     read(21,*) BIN
     read(21,*) BIN
     read(39,*) BIN,en(t)
!write(6,*) 'energy',en(t),t     
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

  end do


do lx = 1,maxlink*binsize,1
        rr = lx/real(binsize)
write(49,*) rr,summing(1,lx)/totmonoc,summing(2,lx)/totdic,summing(3,lx)/totclusc

!write(29,*) real(totmono)/totmonoc,real(totdi)/totdic,real(totclus)/totclusc
end do

write(29,*)real(totmono)/totmonoc,real(totdi)/totdic,real(totclus)/totclusc,totbulk/totbulkc


contains




  subroutine clustercount
    integer :: m,l,g,f,clcount,clusterpop,bdir,maxlengthss,maxback,clustcount,content,onepop,z
    integer ::acconn,cconn,atconn,tconn,a,b,c,d,sepmonocount,sepdicount,sepmaxcount,sss
    integer::sepx,sepy,sepz
    integer,dimension(:),allocatable:: clnos
    logical :: clusyes,adjver,newcl
    type(rprot),dimension(:),allocatable :: tempcoord
    integer,dimension(:),allocatable:: cllist,histcl!,mashist
    integer :: oldcl,dum1,dum2,zzz,maxclus,normalisedsize,dummy,unbcount,bcount
integer :: bulkcount
double precision :: bulkmax    
integer,dimension(:,:),allocatable::cb,conn
    integer,dimension(:),allocatable::maxcluslist,concount
    integer,dimension(:,:),allocatable::bound
 double precision::maxclusenergy,initialenergy,sepmono,sepdi,sepmax,totsep


    allocate(cb(nprotein,nprotein))
    allocate(clnos(nprotein))
    allocate(tempcoord(maxlength))
    allocate(cllist(nprotein))
    allocate(histcl(nprotein))
    allocate(conn(nprotein,maxlength))
    allocate(concount(nprotein))
    allocate(bound(nprotein,maxlength))
    !allocate(mashist(10))
sepmonocount = 0
sepdicount = 0
sepmaxcount = 0
sepmono = 0.0d0
sepdi = 0.0d0
sepmax = 0.0d0
bulkmax = 0.0d0
bulkcount = 0
    maxclusenergy = 0.0d0
    initialenergy = 0.0
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
          bound(m,a) = 0
       end do
    end do

    a = 1

    do m = 1,nprotein
       do g = 1,nprotein
          cb(m,g) = 0
       end do
    end do

        do m=1,nprotein,1
                maxl = chlen(m)
                do l = 1,maxl,1
                        tempcoord(l)%x = chains(m,l)%x
                        tempcoord(l)%y =  chains(m,l)%y
                        tempcoord(l)%z = chains(m,l)%z

                        do f = l+1,maxl,1
                                call adjacent(m,l,f,tempcoord,adjver,bdir)

                                if((adjver .eqv. .true.)) then
                                        initialenergy = initialenergy+intraen

                                        if((bdir == chains(m,f)%bond) .and. (bdir==(-1*chains(m,l)%bond))) then

                                                initialenergy = initialenergy + intraenergy(chains(m,l)%type,chains(m,f)%type)
                                        end if
                                end if
                        end do
                end do
        end do

        !write(6,*) 'INITIAL ENERGY', initialenergy

    do m = 1,nprotein-1,1
       maxback = chlen(m)
       do g = m+1,nprotein,1
         ! if((g/=m) .and. (clnos(m) /= clnos(g))) then
             maxlengthss = chlen(g)
             do l = 1,maxback,1
                        tempcoord(l)%x = chains(m,l)%x
                        tempcoord(l)%y =  chains(m,l)%y
                        tempcoord(l)%z = chains(m,l)%z

                do f = 1,maxlengthss,1
                   call adjacent(g,l,f,tempcoord,adjver,bdir)              
                   if(adjver.eqv. .true.) then
                      if(isbond .eqv. .true.) then
                         !cb(m,g) = 1
                         !cb(g,m) = 1
                         !bound(m,l) = 1
                         !bound(g,f) = 1
                            initialenergy = initialenergy + interen(chains(m,l)%type,chains(g,f)%type)    
                                 if(interen(chains(m,l)%type,chains(g,f)%type) < 0.0) then

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
                       cb(m,g) = 1
                       cb(g,m) = 1
 
                clusyes = .true.
                         end if
                      end if
                      if((bdir == chains(g,f)%bond) .and. (bdir == (-1*chains(m,l)%bond)) .and. &
                           (interenergy(chains(g,f)%type,chains(m,l)%type)  < 0.0)) then
                         initialenergy = initialenergy + interenergy(chains(m,l)%type,chains(g,f)%type) 
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
    do m = 1,nprotein
       if(clnos(m) == content)then
          maxcluslist(z) = m
          z = z+1
       end if
    end do


    a = 1
    b = 1
    c = 1
    d = 1

    do m = 1,nprotein
       maxl = chlen(m)
       do l = 1,maxl
          if(bound(m,l) == 0) then
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
       end do
    end do
    write(6,*) 'aclus',a,'bclus',b,'abulk',c,'bbulk',d,en(t),maxclus

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

    
    !call paircorrelation(a,b,c,d)
!write(6,*) 'energy',en(1),timeofinterest
    if (nint(initialenergy) /= nint(en(t))) then
       write(6,*) t,initialenergy,en(t)
    end if




  do m = 1,nprotein
       maxl = chlen(m)
        !if(histcl(m) == 1) then
        do l = 2,maxl
sepx =min(abs(chains(m,l)%x-chains(m,l-1)%x),gridsize-abs(chains(m,l)%x-chains(m,l-1)%x))
sepy=min(abs(chains(m,l)%y-chains(m,l-1)%y),gridsize-abs(chains(m,l)%y-chains(m,l-1)%y))
sepz=min(abs(chains(m,l)%z-chains(m,l-1)%z),gridsize-abs(chains(m,l)%z-chains(m,l-1)%z))

totsep = sqrt((real(sepx)**2)+(real(sepy)**2) +(real(sepz)**2))
sepmono = sepmono+totsep
sepmonocount = sepmonocount+1


sss = int(totsep*binsize)
summing(1,sss) = summing(1,sss)+1.0
end do
!end if
        if(histcl(m) == 2) then
        do l = 2,maxl

sepx=min(abs(chains(m,l)%x-chains(m,l-1)%x),gridsize-abs(chains(m,l)%x-chains(m,l-1)%x))
sepy=min(abs(chains(m,l)%y-chains(m,l-1)%y),gridsize-abs(chains(m,l)%y-chains(m,l-1)%y))
sepz=min(abs(chains(m,l)%z-chains(m,l-1)%z),gridsize-abs(chains(m,l)%z-chains(m,l-1)%z))

totsep = sqrt((real(sepx)**2)+(real(sepy)**2) +(real(sepz)**2))
sepdi = sepdi+totsep
sepdicount = sepdicount+1



sss = int(totsep*binsize)
summing(2,sss) = summing(2,sss)+1.0
end do
end if

        if(clnos(m) == content) then
        do l = 2,maxl


sepx=min(abs(chains(m,l)%x-chains(m,l-1)%x),gridsize-abs(chains(m,l)%x-chains(m,l-1)%x))
sepy=min(abs(chains(m,l)%y-chains(m,l-1)%y),gridsize-abs(chains(m,l)%y-chains(m,l-1)%y))
sepz=min(abs(chains(m,l)%z-chains(m,l-1)%z),gridsize-abs(chains(m,l)%z-chains(m,l-1)%z))

totsep = sqrt((real(sepx)**2)+(real(sepy)**2) +(real(sepz)**2))
sepmax = sepmax+totsep
sepmaxcount = sepmaxcount+1

sss = int(totsep*binsize)
!write(6,*) 'sss',sss,totsep,int(0.5)
summing(3,sss) = summing(3,sss)+1.0
end do

else if(clnos(m) /= content) then
        do l = 2,maxl


sepx=min(abs(chains(m,l)%x-chains(m,l-1)%x),gridsize-abs(chains(m,l)%x-chains(m,l-1)%x))
sepy=min(abs(chains(m,l)%y-chains(m,l-1)%y),gridsize-abs(chains(m,l)%y-chains(m,l-1)%y))
sepz=min(abs(chains(m,l)%z-chains(m,l-1)%z),gridsize-abs(chains(m,l)%z-chains(m,l-1)%z))

totsep = sqrt((real(sepx)**2)+(real(sepy)**2) +(real(sepz)**2))
bulkmax = bulkmax+totsep
bulkcount = bulkcount+1

!sss = int(totsep*binsize)
!write(6,*) 'sss',sss,totsep,int(0.5)
!summing(3,sss) = summing(3,sss)+1.0
end do



end if






end do
           open(113,file = 'beadspace.dat',access = 'append')
write(113,*) time,sepmono/sepmonocount,sepdi/sepdicount,sepmax/sepmaxcount
close(113)

totmono = totmono+ sepmono
totmonoc = totmonoc + sepmonocount
totdi = totdi + sepdi
totdic = totdic + sepdicount
totclus = totclus + sepmax
totclusc = totclusc + sepmaxcount
totbulk = totbulk + bulkmax
totbulkc = totbulkc + bulkcount
    deallocate(cllist)
    deallocate(clnos)
    deallocate(tempcoord)
    deallocate(histcl)
  end subroutine clustercount


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



  subroutine paircorrelation(a,b,c,d)
    integer,intent(in)::a,b,c,d
    integer::k,h
    double precision::dx,dy,dz

    if(a>1) then

       do h = 1,a-1
          do k = 1,a-1

             if(k==h) goto 75

             dx = min(abs(chains(aclus(k)%m,aclus(k)%l)%x-chains(aclus(h)%m,aclus(h)%l)%x),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%x-chains(aclus(h)%m,aclus(h)%l)%x))
             dy = min(abs(chains(aclus(k)%m,aclus(k)%l)%y-chains(aclus(h)%m,aclus(h)%l)%y),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%y-chains(aclus(h)%m,aclus(h)%l)%y))
             dz = min(abs(chains(aclus(k)%m,aclus(k)%l)%z-chains(aclus(h)%m,aclus(h)%l)%z),gridsize-&
                  abs(chains(aclus(k)%m,aclus(k)%l)%z-chains(aclus(h)%m,aclus(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             if(lx ==0) lx =1
             if(abs(lx)>(gridsize*binsize)) goto 75

             summing(1,lx) = summing(1,lx)+(1.0) !/((a-1)**2))

75           continue
          end do
       end do
    end if

    if(b>1) then
       do h = 1,b-1
          do k = 1,b-1

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
             if(lx ==0) lx =1
             if(abs(lx)>(gridsize*binsize)) goto 95

             summing(2,lx) = summing(2,lx)+(1.0) !/((b-1)**2))

95           continue
          end do
       end do
    end if

    if(c>1) then
       do h = 1,c-1
          do k = 1,c-1

             if(k==h) goto 45

             dx = min(abs(chains(abulk(k)%m,abulk(k)%l)%x-chains(abulk(h)%m,abulk(h)%l)%x),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%x-chains(abulk(h)%m,abulk(h)%l)%x))
             dy = min(abs(chains(abulk(k)%m,abulk(k)%l)%y-chains(abulk(h)%m,abulk(h)%l)%y),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%y-chains(abulk(h)%m,abulk(h)%l)%y))
             dz = min(abs(chains(abulk(k)%m,abulk(k)%l)%z-chains(abulk(h)%m,abulk(h)%l)%z),gridsize-&
                  abs(chains(abulk(k)%m,abulk(k)%l)%z-chains(abulk(h)%m,abulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             if(lx ==0) lx =1
             if(abs(lx)>(gridsize*binsize)) goto 45

             summing(3,lx) = summing(3,lx)+(1.0) !/((c-1)**2))

45           continue
          end do
       end do
    end if

    if(d>1) then
       do h = 1,d-1
          do k = 1,d-1

             if(k==h) goto 55

             dx = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%x-chains(bbulk(h)%m,bbulk(h)%l)%x),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%x-chains(bbulk(h)%m,bbulk(h)%l)%x))
             dy = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%y-chains(bbulk(h)%m,bbulk(h)%l)%y),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%y-chains(bbulk(h)%m,bbulk(h)%l)%y))
             dz = min(abs(chains(bbulk(k)%m,bbulk(k)%l)%z-chains(bbulk(h)%m,bbulk(h)%l)%z),gridsize-&
                  abs(chains(bbulk(k)%m,bbulk(k)%l)%z-chains(bbulk(h)%m,bbulk(h)%l)%z))


             delta = sqrt((dx**2)+(dy**2)+(dz**2))

             lx = INT(delta/(1.0/binsize))
             if(lx ==0) lx =1
             
             if(abs(lx)>(gridsize*binsize))     goto 55
             
             summing(4,lx) = summing(4,lx)+(1.0) !/((d-1)**2))

55           continue
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
allocate(linklen(nspecies))
       CASE ('TOTCHAINS')
          !CALL get_integer(nprotein)

          runtype = 1

       CASE ('NEWCHAIN')
          CALL get_integer(specieslen(runtype))
          CALL get_integer(speciespop(runtype))
          call get_integer(linklen(runtype))   
          do f = 1,specieslen(runtype)
             call get_integer(dummytype)
             if(dummytype>mtype) mtype = dummytype
          end do
          write(6,*) 'runtype',runtype,specieslen(runtype),speciespop(runtype)
          runtype = runtype + 1



maxlink = maxval(linklen)
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
          intraen = intraen
write(6,*) 'INTRA_ENERGY',intraen

       CASE ('INTERACTION')
          call get_integer(xl)
          call get_dp(intra)
          write(6,*) 'mtype',mtype
          do f = 1,mtype
             call get_dp(interenergy(xl,f))
             intraenergy(xl,f) = intra
             interenergy(xl,f) = -abs(interenergy(xl,f))
             call get_dp(interen(xl,f))
             !interen(xl,f) = -abs(interen(xl,f))
          end do

write(6,*) 'INTRA - dir',intraenergy(:,:)

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



