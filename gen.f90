program main
  implicit none
  real*8 :: temp
  character(kind=1,len=3) :: atemp,aid
  character(kind=1,len=99) :: ainp(37),binp
  character(kind=1,len=10) :: cinp
  integer :: id,j
  real*8 :: xrand,yrand

  call random_seed()

 !temp=300.d0
  read*,temp
  print*,temp
  write(atemp,'(i3.3)') int(temp)
  call system("mkdir run-"//atemp)
  open(101,file='/home/chenjun/CO-Ru0001/3-xMD/2-xMD-HD/src-3layers-binary/INPUT.3layer',status='old',action='read')
  do j=1,37
     read(101,'(a99)') ainp(j)
  enddo
  close(101)

  do id=1,8
     write(aid,'(i3.3)') id

     call random_number(xrand); xrand=xrand*5.d0
     call random_number(yrand); yrand=yrand*5.d0

     binp=ainp(32)
     write(cinp,'(f10.5)') xrand
     binp(20:29)=cinp
     write(cinp,'(f10.5)') yrand
     binp(30:39)=cinp
     ainp(32)=binp

     binp=ainp(33)
     write(cinp,'(f10.5)') xrand
     binp(20:29)=cinp
     write(cinp,'(f10.5)') yrand
     binp(30:39)=cinp
     ainp(33)=binp

     binp=ainp(10)
     write(cinp,'(f10.2)') temp
     binp(1:10)=cinp
     ainp(10)=binp

     call system("mkdir run-"//atemp//"/run-"//aid)
     open(200+id,file="run-"//atemp//"/run-"//aid//"/INPUT",status='unknown')
     do j=1,37
        write(200+id,*) trim(ainp(j))
     enddo
     close(200+id)
     call system("ln -s /home/chenjun/CO-Ru0001/3-xMD/2-xMD-HD/src-3layers-binary/pes8.dat run-"//atemp//"/run-"//aid//"/pes8.dat")
     call system("ln -s /home/chenjun/CO-Ru0001/3-xMD/2-xMD-HD/src-3layers-binary/xmd.x run-"//atemp//"/run-"//aid//"/xmd_220405")
  enddo

end program
