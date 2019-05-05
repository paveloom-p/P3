module subprograms
implicit none

contains


subroutine calculation(i,top_era,era_incr) ! Перерасчёт координат на выбранную эпоху
implicit none

real(8) top_era
integer era_incr

! Формулы, указанные в руководстве Цветкова:
! alpha=aplha_0+(T-T_0)*mu_alpha ,
! delta=delta_0+(T-T_0)*mu_delta ,
! где (T-T_0) принимаю за era.

! Причем mu_alpha и mu_delta идут в размерности mas/yr,
! где mas - милли- угловая секунда.
! 
! 1 mas = 0.00000028 градуса

real(8) pmRA, pmDE, RAdeg, DEdeg
integer i, j, k

open(6, file='prmotions'); open(8, file='data_2')

do j=1, i-1, 1 
        read(6,*) pmRA, pmDE, RAdeg, DEdeg
                do k=era_incr, 1, -1
                write(8,*) "'", RAdeg+top_era/k*0.00000028*pmRA, DEdeg+top_era/k*0.00000028*pmDE, "'"
                write(4,*) '           ', RAdeg+top_era/k*0.00000028*pmRA, DEdeg+top_era/k*0.00000028*pmDE
                enddo
        write(8,*) "' '"
        write(4,*) ' '
enddo

close(8); close(6)

end subroutine calculation

!-------------------------------

subroutine plotting_data(i,order,era_incr) ! Подготовка данных к выводу
implicit none

integer i, j, k, era_incr
integer order(:)
character(100) a((i-1)*era_incr+(i-2)) ! Формальный массив для ввода/вывода data_2
character(100) b(i) ! Формальный массив для ввода/вывода data_1



open(5, file='data_1') ! Считывание data_1
        do j=1, i-1, 1
        read(5,*) b(j)
        enddo
close(5)

open(8, file='data_2') ! Считывание data_2
        do j=1, (i-1)*era_incr+(i-2), 1
        read(8,*) a(j)
        enddo
close(8)

open(9, file='plotting_data') ! Запись в файл для вывода
        do j=1, size(order), 1
                 write(9,*) b(order(j))        
        enddo

        write(9,'(/)')
         
        do k=1, era_incr, 1
        do j=1, size(order), 1
                 write(9,*) a(order(j)+(order(j)-1)*era_incr+(k-1)) 
        enddo
        write(9,*) ' '
        enddo
close(9)

end subroutine plotting_data

end module subprograms
