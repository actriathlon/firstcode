program thinning
  implicit none
  integer l,m,over,num,nprot,len,bead,p,q,acount,nosteps,freq
  double precision :: x,y,z,test,boxl
  character(len = 20):: dummy
      character(len=8) :: atnum

  nprot = 480   !number of proteins in the system
  len = 50     !2x the length of the chains
  len2 =24      !length of chain -1
  nosteps =726   ! number of output steps of the original file
  freq =10   !factor by which file is to be reduced by
  
  open(23, file = 'move.vtf', action = 'read')
       open(43, file = 'movethin.vtf', action = 'write')


       write(43,*) 'atom ', 'default ', 'radius ', 1.00000, 'name ','C'
       read(23,*) dummy,dummy,dummy,test,dummy,dummy



       
  do m = 1,nprot
     do l = 1,len
        read(23,*) dummy,bead,dummy,test,dummy,p,dummy,q
        write(43,*) 'atom',bead,'radius',1.00000000,'name',p,'resid',q
     end do
     end do
     read(23,*)
     write(43,*) ' '
     
     acount = 0
    do m= 1,nprot
       do l = 1,len2
          write(atnum,'(i7)') acount
          read(23,*)
          read(23,*)
          write(43,*) 'bond', adjustr(atnum)//':',acount+2
          write(43,*) 'bond', adjustr(atnum)//':',acount+1
          acount = acount+2
       end do
       read(23,*)
       write(atnum,'(i7)') acount
       write(43,*) 'bond', adjustr(atnum)//':',acount+1
       acount = acount +2
    end do

    
       
  do over = 1,nosteps
     if(modulo(over,freq) == 0) then
        read(23,*)
        write(43,*) ' '
        read(23,*) dummy
        write(43,*) 'timestep indexed'
        read(23,*) dummy,boxl,boxl,boxl
        !write(6,*) 'pbc',boxl,boxl,boxl
  do m = 1,nprot
     do l = 1,len
        read(23,*) num,x,y,z
        write(43,*) num,x,y,z
     end do
  end do
else
   read(23,*)
   !write(6,*) 
   read(23,*) 
   read(23,*)dummy
   !write(6,*) dummy
        do m = 1,nprot
     do l = 1,len
        read(23,*) num,x,y,z
     end do
  end do
  end if

end do
end program thinning
