module functions
contains

!-----------------------------

function sxy(a,k)
integer a(:,:,:), k

write(*,*) a(:,:,k)

end function sxy

!-----------------------------

function sxz(a,k)
integer a(:,:,:), k

write(*,*) a(:,k,:)

end function sxz

!-----------------------------

function syz(a,k)
integer a(:,:,:), k

write(*,*) a(k,:,:)

end function syz

end module
