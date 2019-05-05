module functions
contains

!----------------------------- ФУНКЦИЯ должна возвращать через своё имя
!                              то, что ей положено возвращать; в данном 
!function sxy(a,k)             ! случае матрицу, являющуюся содержимым
!integer a(:,:,:), k           ! k-ого листа книги (a).

                              ! Если провести компиляцию при включённой
                              ! опции -Wall (вывод предупреждений о
                              ! сомнительных моментах, которые, тем не менее,
                              ! всё-таки позволили построить исполнимый
                              ! код), то увидим:
!  gfortran -c -Wall -Wtabs functions.f95
!  functions.f95:24:
!
!  function syz(a,k)
!  1
!  Предупреждение: Return value of function 'syz' at (1) not set
!  functions.f95:15:
!
!  function sxz(a,k)
!  1
!  Предупреждение: Return value of function 'sxz' at (1) not set
!  functions.f95:6:
!
!  function sxy(a,k)
!  1
!  Предупреждение: Return value of function 'sxy' at (1) not set
!
! т.е. компилятор предупреждает, что несмотря на то, что описаные
! объекты являются по синтаксису функциями, тем не менее, они ничено не 
! возвращают через своё имя, что, вообще говоря, приводит к вопросу:
!
!            А зачем тогда они описывались как функции?
!
! Если хотим просто вывести нужную матрицу k-го листа, то для этого
! естественно было написать процедуру вывода:  subroutine wrt_sxy(a,k)
! причём так, чтобы выводимое было легко сопоставить с трёхмерным 
! массивом a (т.е. чтобы при сравнении A c выводимым было просто 
! сообразить, как оно расположено в A.
!
! При Вашей форме вывода нужно ещё помнить, что ФОРТРАН-матрицы по
! умолчанию выводятся по стольцам. Поэтому выведенное в строку
! на самом деле будет содержать столбец. А это ещё одна дополнительная
! нагрузка на наши "извилины": 
!
!   Результат должен быть НАГЛЯДЕН, т.е. обладал свойством восприниматься
!   человеком без особого напряга мыслительного процесса --- где
!   ошибиться мы прекрасно найдём и без этой дополнительной нагрузки ---
!   так зачем же ещё расширять поле возможных ошибок.
!
! Короче:  
!    function sxy1(a,k) result(s)
!    integer a(:,:,:), s(size(a,1),size(a,2))
!    s=a(:,:,k);
!    end function sxy1
!
!    subroutine wrt_sxy1(s)
!    integer s(:,:), nx,ny
!    character(3) sn
!    character(60) sf
!    nx=size(s,1); ny=size(s,2)
!    write(sf,'(i3)') ny
!    sf='('//sf//'i3)'; write(*,trim(sf)) ((s(i,j),j=1,ny),i=1,nx)
!    end subroutine wrt_sxy1

 
                              ! Вопрос о выводе этой матрицы решается
!write(*,*) a(:,:,k)          ! не функцией sxy, а программной единицей,
                              ! вызвавшей sxy. Нужно --- распечатает,
                              ! не нужно --- будет использовать в дальнейших
                              ! расчётах, если таковые имеются.
!end function sxy

function sxy1(a,k) result(s)
integer a(:,:,:), s(size(a,1),size(a,2))
s=a(:,:,k);
end function sxy1

subroutine wrt_sxy1(s)
integer s(:,:), nx,ny
character(3) sn
character(60) sf
nx=size(s,1); ny=size(s,2); !  write(*,*) ' nx=',nx,'    ny=',ny
write(sn,'(i3)') ny         !  write(*,*) ' sn=',sn 
sf='('//sn//'i7)';          !  write(*,*) ' sf=',trim(sf)
write(*,trim(sf)) ((s(i,j),j=1,ny),i=1,nx)
write(*,*) ' '
end subroutine wrt_sxy1


!----------------------------- ! Пожалуйста, по аналогии с sxy1 и wrt_sxy1
                               ! напишите соответствующие функции для 
function sxz1(a,k) result(s)              ! двух оставшихся сечений.
integer a(:,:,:), s(size(a,1), size(a,3))
s=a(:,k,:)
end function sxz1              ! П: Написал.

subroutine wrt_sxz1(s)
integer s(:,:), nx,ny
character(3) sn
character(60) sf
nx=size(s,1); ny=size(s,2);
write(sn,'(i3)') ny
sf='('//sn//'i7)';
write(*,trim(sf)) ((s(i,j),j=1,ny),i=1,nx)
write(*,*) ' '
end subroutine wrt_sxz1

!-----------------------------

function syz1(a,k) result(s)
integer a(:,:,:), s(size(a,2), size(a,3))
s=a(k,:,:)
end function syz1

subroutine wrt_syz1(s)
integer s(:,:), nx,ny
character(3) sn
character(60) sf
nx=size(s,1); ny=size(s,2);
write(sn,'(i3)') ny
sf='('//sn//'i7)';
write(*,trim(sf)) ((s(i,j),j=1,ny),i=1,nx)
write(*,*) ' '
end subroutine wrt_syz1

!-----------------------------

subroutine wrt_page(a); integer a(:,:,:), sh(3)
integer i, j, k
sh=shape(a); write(*,*) ' wrt_page(a):'
do k=1,sh(3)
  write(*,'(4i7)') ((a(i,j,k),j=1,sh(2)),i=1,sh(1)); write(*,*)
enddo
end subroutine wrt_page

end module
