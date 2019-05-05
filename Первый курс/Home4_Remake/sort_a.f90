module sort; implicit none; contains
!...   ...   ...   ...   ...   ...
subroutine alter1(a,b,c,na,nb,nc)   ! �������� ��������� � ����
integer na, nb, nc, i               ! � ����� ��������-77      
integer a(na), b(nb), c(nc)
if (na<=nb) then; do i=1, na;   c(2*i-1)=a(i); c(2*i)  =b(i); enddo
                  do i=na+1,nb; c( na+i)=b(i); enddo
            else; do i=1, nb;   c(2*i-1)=a(i); c(2*i)  =b(i); enddo
                  do i=nb+1,na; c( nb+i)=a(i); enddo
endif
end  subroutine alter1
subroutine alter2(a,b,c,na,nb,nc)   ! ��������� � ����� ��������-77,
integer na, nb, nc, i               ! �� ��������� ��������-95:
integer a(na), b(nb), c(nc)         ! ��������� ����������, ������ �����.
if (na<=nb) then; c(1:2*na-1 :2)=a(1:na); c(2: 2*na  :2)=b(1:na)
                                          c(2*na+1:nc:1)=b(na+1:nb)
            else; c(1:2*nb-1:2)=a(1:nb);  c(2:2*nb:2)=b(1:nb)
                                          c( 2*nb+1:nc)=a(nb+1:na)
endif
end  subroutine alter2
subroutine alter3(a,b,c)             ! ��������� � ����� ��������-95 
integer a(:), b(:), c(:), na, nb, nc ! ����� �������, �������������  
na=size(a); nb=size(b); nc=size(c)   ! ����� ������� ������� size.   
if (na<=nb) then; c(1:2*na-1 :2)=a(1:na); c(2: 2*na  :2)=b(1:na)
                                          c(2*na+1:nc:1)=b(na+1:nb)
            else; c(1:2*nb-1:2)=a(1:nb);  c(2: 2*nb  :2)=b(1:nb)
                                          c( 2*nb+1:nc)=a(nb+1:na)
endif
end  subroutine alter3
function alter4(a,b) result(c)         ! ������� ��������-95 
integer a(:), b(:), c(size(a)+size(b)) ! ����� ���������� �����
integer na, nb, nc                     ! ��ϣ  ���   ������! 
na=size(a); nb=size(b); nc=size(c)
if (na<=nb) then; c(1:2*na-1 :2)=a(1:na); c(2: 2*na  :2)=b(1:na)
                                          c(2*na+1:nc:1)=b(na+1:nb)
            else; c(1:2*nb-1:2)=a(1:nb);  c(2: 2*nb  :2)=b(1:nb)
                                          c( 2*nb+1:nc)=a(nb+1:na)
endif
end function alter4
end module sort
