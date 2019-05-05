module quadglobal
use mpmodule
implicit none
integer ndebug, nerror, nquadl, ndigits1, ndigits2, nepsilon1, nepsilon2, &
  nwords1, nwords2
end module

subroutine f_main

!   David H. Bailey   8 April 2017

!   This is the ARPREC Fortran-90 version for a single processor.

!   This software is provided for research purposes only.  
!   Commercial usage requires license agreement.

!   This work was supported by the Director, Office of Science, Division
!   of Mathematical, Information, and Computational Sciences of the
!   U.S. Department of Energy under contract number DE-AC03-76SF00098.

!   This program demonstrates the 2D quadrature routine 'quadts2d', which employs
!   the error function.  The function quadts2d is suitable to integrate
!   a function that is continuous, infinitely differentiable and integrable on a
!   2D finite open interval.  It can also be used for certain integrals on
!   infinite intervals, by making a suitable change of variable -- see below.

!   The function(s) to be integrated is(are) defined in external function
!   subprogram(s) -- see the sample function subprograms below.  The name(s) of
!   the function subprogram(s) must be included in appropriate type and external
!   statements in the main program.

!   Inputs set in parameter statement below:
!   kdebug Debug level setting.  Default = 2.
!   ndp    Digits of precision.  May not exceed mpipl in file mpmod90.f.
!            In some cases, ndp must be significantly greater than the desired
!            tolerance in the result-- see the examples below.
!   neps   Log10 of the desired tolerance in the result (negative integer).
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!            (possibly doubles) the number of accurate digits in the result,
!            but also roughly quadruples the run time.  nq1 > 2.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!            default it is set to 8 * 2^nq1.  Increase nq2 if directed by a 
!            message produced in initqts.  Note that the dimension of the
!            wk and xk arrays starts with -1, so the length of these arrays is
!            (nq2+2) * 4 eight-byte words.

use mpmodule
use mpmodulex
use quadglobal
implicit none
integer i, kdebug, lev, ndp1, ndp2, neps1, neps2, nq1, nq2, n, n1
parameter (kdebug = 2, ndp1 = 120, ndp2 = 240, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 9, nq2 = 8 * 2 ** nq1)
double precision dplog10q, d1, d2, second, tm0, tm1, tm9
type (mp_real) c2, c3, cat, catalan, err, quadts2d, fun1, fun2, &
  fun3, fun4, fun5, fun6, fun7, fun8, t1, t2, t3, t4, wk(-1:nq2), xk(-1:nq2)
type (mp_realx) fun11, fun12, fun21, fun22, fun31, fun32, fun41, fun42, &
  fun51, fun52, fun61, fun62, fun71, fun72, fun81, fun82, y1, y2
external quadts2d, catalan, fun1, fun11, fun12, fun2, fun21, fun22, &
  fun3, fun31, fun32, fun4, fun41, fun42, fun5, fun51, fun52, fun6, &
  fun61, fun62, fun7, fun71, fun72, fun8, fun81, fun82, second
save

call mpinit (ndp1)
call mpinitx (ndp2)
call mpsetprec (ndp1)
call mpgetprecwords (nwords1)
call mpsetprec (ndp2)
call mpgetprecwords (nwords2)
call mpsetprecwords (nwords1)
ndigits1 = ndp1
ndigits2 = ndp2
nepsilon1 = neps1
nepsilon2 = neps2
mpoud = ndp1
ndebug = kdebug
nerror = 0
nquadl = nq1
write (6, 1) nquadl, ndigits1, ndigits2, nepsilon1, nepsilon2
1 format ('Quadts2d test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,' Epsilon2 =',i6)

!   Initialize quadrature tables wk and xk (weights and abscissas).

tm0 = second ()
call initqts (nq1, nq2, wk, xk)
tm1 = second ()
if (nerror > 0) stop
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)
cat = catalan ()
tm9 = tm1 - tm0
c2 = 2.d0
c3 = 3.d0

!   Begin quadrature tests.

write (6, 11)
11 format (/ &
  'Problem 1: Int_0^1 Int_0^1 sqrt(s^2+t^2) ds dt')
lev = 9
write (6, 3) lev
3 format ('Level =',i4)
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadts2d (fun1, fun11, fun12, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 4) tm1 - tm0
4 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call mpwrite (6, t1)
t2 = sqrt (c2) / 3.d0 - log (c2) / 6.d0 + log (c2 + sqrt (c2)) / 3.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
5 format ('Actual error =',f10.6,'x10^',i5)
tm9 = tm9 + (tm1 - tm0)

write (6, 12)
12 format (/&
  'Problem 2: Int_0^1 Int_0^1 sqrt(1+(s-t)^2) ds dt')
lev = 6
write (6, 3) lev
call mpsetprecwords (nwords2)
y1 = 0.d0
y2 = 1.d0
call mpsetprecwords (nwords1)
tm0 = second ()
t1 = quadts2d (fun2, fun21, fun22, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
t2 = -sqrt (c2) / 3.d0 - log (sqrt (c2) - 1.d0) / 2.d0 &
     + log (sqrt (c2) + 1.d0) / 2.d0 + c2 / 3.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 13)
13 format (/ &
  'Problem 3: Int_-1^1 Int_-1^1 1/sqrt(1+s^2+t^2) ds dt')
lev = 7
write (6, 3) lev
y1 = -1.d0
y2 = 1.d0
tm0 = second ()
t1 = quadts2d (fun3, fun31, fun32, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
t2 = 4.d0 * log (2.d0 + sqrt (c3)) - 2.d0 * mppic / 3.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 14)
14 format (/&
  'Problem 4: Int_0^pi Int_0^pi log (2-cos(s)-cos(t)) ds dt')
lev = 9
write (6, 3) lev
call mpsetprecwords (nwords2)
y1 = 0.d0
y2 = mppicx
call mpsetprecwords (nwords1)
tm0 = second ()
t1 = quadts2d (fun4, fun41, fun42, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
t2 = 4.d0 * mppic * cat - mppic**2 * log (c2)
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 15)
15 format (/&
  'Problem 5: Int_0^inf Int_0^inf sqrt(s^2+s*t+t^2) * exp(-s-t) ds dt')
lev = 9
write (6, 3) lev
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadts2d (fun5, fun51, fun52, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
t2 = 1.d0 + 0.75d0 * log (c3)
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 16)
16 format (/&
  'Problem 6: Int_0^1 Int_0^1 1/(sqrt((1-s)*(1-t))*(s+t)) ds dt')
lev = 9
write (6, 3) lev
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadts2d (fun6, fun61, fun62, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
t2 = 4.d0 * cat
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 17)
17 format (/ &
'Problem 7: Int_0^1 Int_0^t 1/sqrt(1+s^2+t^2) ds dt')
lev = 6
write (6, 3) lev
y1 = 0.d0
y2 = 1.d0
tm0 = second ()
t1 = quadts2d (fun7, fun71, fun72, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
t2 = -mppic / 12.d0 - log (c2) / 2.d0 + log (1.d0 + sqrt (c3))
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 18)
18 format (/ &
  'Problem 8: Int_0^pi Int_0^y cos(s)*sin(t)/exp[s+t] ds dt')
lev = 6
write (6, 3) lev
call mpsetprecwords (nwords2)
y1 = 0.d0
y2 = mppicx
call mpsetprecwords (nwords1)
tm0 = second ()
t1 = quadts2d (fun8, fun81, fun82, y1, y2, lev, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 4) tm1 - tm0
call mpwrite (6, t1)
t2 = 0.25d0 * (1.d0 + exp (-mppic))
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
tm9 = tm9 + (tm1 - tm0)

write (6, 6) tm9
6 format ('Total CPU time =',f12.6)
stop
end

function fun1 (s, t)

!   fun1 (s,t) = sqrt[s^2+t^2]

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun1, s1, t1
type (mp_realx) s, t

s1 = s
t1 = t
fun1 = sqrt (s**2 + t**2)
return
end

function fun11 (s)
use mpmodulex
implicit none
type (mp_realx) fun11, s
fun11 = 0.d0
return
end

function fun12 (s)
use mpmodulex
implicit none
type (mp_realx) fun12, s
fun12 = 1.d0
return
end

function fun2 (s, t)

!   fun2 (s,t) = Sqrt[1 + (s - t)^2]

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun2, s1, t1
type (mp_realx) s, t

s1 = s
t1 = t
fun2 = sqrt (1.d0 + (s1 - t1)**2)
return
end

function fun21 (s)
use mpmodulex
implicit none
type (mp_realx) fun21, s
fun21 = 0.d0
return
end

function fun22 (s)
use mpmodulex
implicit none
type (mp_realx) fun22, s
fun22 = 1.d0
return
end

function fun3 (s, t)

!   fun3 (s,t) = 1/sqrt[1+s^2+t^2]

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun3, s1, t1
type (mp_realx) s, t

s1 = s
t1 = t
fun3 = 1.d0 / sqrt (1.d0 + s**2 + t**2)
return
end

function fun31 (s)
use mpmodulex
implicit none
type (mp_realx) fun31, s
fun31 = -1.d0
return
end

function fun32 (s)
use mpmodulex
implicit none
type (mp_realx) fun32, s
fun32 = 1.d0
return
end

function fun4 (s, t)

!   fun4 (s,t) = log (2 - cos(s) - cos(t))
!               = log (2*sin(s/2)^2 + 2*sin(t/2)^2)

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun4, s1, t1
type (mp_realx) s, t

s1 = s
t1 = t
! fun4 = log (2.d0 - cos (s1) - cos (t1))
fun4 = log (2.d0 * sin (0.5d0 * s1) ** 2 + 2.d0 * sin (0.5d0 * t1) ** 2)
return
end

function fun41 (s)
use mpmodulex
implicit none
type (mp_realx) fun41, s
fun41 = 0.d0
return
end

function fun42 (s)
use mpmodulex
implicit none
type (mp_realx) fun42, s
fun42 = mppic
return
end

function fun5 (s, t)

!   fun5 (s,t) = sqrt ((1/s-1)^2 + (1/s-1)*(1/t-1) + (1/t-1)^2) 
!                / (s^2 * t^2 * exp(1/s + 1/t - 2))

use mpmodule
use mpmodulex
use quadglobal
implicit none
type (mp_real) fun5, s1, t1, s2, t2, sq
type (mp_realx) s, t

s1 = s
t1 = t
if (s1 > 0.d0 .and. t1 > 0.d0) then
  s2 = 1.d0 / s1 - 1.d0
  t2 = 1.d0 / t1 - 1.d0
  if (s2 + t2 < 4.d0 * ndigits1) then
    sq = sqrt (s2**2 + s2 * t2 + t2**2)
    fun5 = sq / (s1**2 * t1**2) * exp (-s2 - t2)
  else
    fun5 = 0.d0
  endif
else
  fun5 = 0.d0
endif

return
end

function fun51 (s)
use mpmodulex
implicit none
type (mp_realx) fun51, s
fun51 = 0.d0
return
end

function fun52 (s)
use mpmodulex
implicit none
type (mp_realx) fun52, s
fun52 = 1.d0
return
end

function fun6 (s, t)

!   fun6 (s,t) = 1/(sqrt((1-s)*(1-t)) * (s+t))

use mpmodule
use mpmodulex
use quadglobal
implicit none
type (mp_real) fun6, s1, t1, s2, t2
type (mp_realx) s, t

s1 = s
t1 = t
call mpsetprecwords (nwords2)
s2 = 1.d0 - s
t2 = 1.d0 - t
call mpsetprecwords (nwords1)
if (s2 > 0.d0 .and. t2 > 0.d0) then
  fun6 = 1.d0 / (sqrt (s2 * t2) * (s1 + t1))
else
  fun6 = 0.d0
endif
return
end

function fun61 (s)
use mpmodulex
implicit none
type (mp_realx) fun61, s
fun61 = 0.d0
return
end

function fun62 (s)
use mpmodulex
implicit none
type (mp_realx) fun62, s
fun62 = 1.d0
return
end

function fun7 (s, t)

!   fun7 (s,t) = 1/sqrt[1+s^2+t^2]

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun7, s1, t1
type (mp_realx) s, t

s1 = s
t1 = t
fun7 = 1.d0 / sqrt (1.d0 + s1**2 + t1**2)
return
end

function fun71 (s)
use mpmodulex
implicit none
type (mp_realx) fun71, s
fun71 = 0.d0
return
end

function fun72 (s)
use mpmodulex
implicit none
type (mp_realx) fun72, s
fun72 = s
return
end

function fun8 (s, t)

!   fun8 (s,t) = cos(s) * cos(t) / exp (s+t)

use mpmodule
use mpmodulex
implicit none
type (mp_real) fun8, s1, t1, s2, t2
type (mp_realx) s, t

s1 = s
t1 = t
fun8 = cos(s1) * sin(t1) / exp(s1+t1)
return
end

function fun81 (s)
use mpmodulex
implicit none
type (mp_realx) fun81, s
fun81 = 0.d0
return
end

function fun82 (s)
use mpmodulex
implicit none
type (mp_realx) fun82, s
fun82 = s
return
end

subroutine initqts (nq1, nq2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = tanh (pi/2*sinh(t)).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
!   for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
!   and (xk(i), i = -1 to n), where n = xk(-1).   The array x_k contains 
!   1 minus the abscissas; the wk array contains the weights at these abscissas.

!   David H Bailey    2004-09-10

use mpmodule
use quadglobal
implicit none
integer i, ierror, iprint, j, k, k1, nq1, nq2, ntab, ntabx
real*8 h
parameter (iprint = 1000)
type (mp_real) eps2, p2, spi, t1, t2, t3, t4, t5, u1, u2, &
  wk(-1:nq2), xk(-1:nq2)

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqts: Error function quadrature initialization')
endif

eps2 = mpreal (10.d0) ** nepsilon2
p2 = 0.5d0 * mppic
h = 0.5d0 ** nq1
wk(-1) = dble (nq1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = dble (k) * h

!   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wk(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xk(k) = 1.d0 / (t3 * t4)
  wk(k) = u1 / t4 ** 2
  if (wk(k) < eps2) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
nerror = 91
goto 130

100 continue

xk(-1) = dble (k)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqts: Error in quadrature initialization; code =',i5)

130  continue

return
end

function quadts2d (fun, fun1, fun2, y1, y2, lev, nq1, nq2, wk, xk)

!   This routine computes the 2-d integral of the function fun:
!       y2   fun2[x]
!     Int  Int  fun[x,y] dx dy
!       y1   fun1[x]
!   with up to nq1 iterations, with a target tolerance of 10^nepsilon1.
!   wk and xk are precomputed tables of abscissas and weights.  The function
!   fun is not evaluated at endpoints.  The array x_k contains 1 minus
!   the abscissas; the wk array contains the weights at these abscissas.

!   David H. Bailey     2004-09-10

use mpmodule
use mpmodulex
use quadglobal
implicit none
integer i, ierror, ip(0:100), izx, iz1, iz2, iz3, iz4, j, k, k1, k2, k3, &
  lev, n, nds, nq1, nq2, nqq1
parameter (izx = 2)
logical logx1, logx2, logy1, logy2
real*8 d1, d2, d3, d4, dwmx, dwi1, dwi2, dwi3, dwi4, dwimx, dwj1(0:nq2), &
  dwj2(0:nq2), dwj3(0:nq2), dwj4(0:nq2), dwjmx, dplog10q, eps1, eps2, err, h
type (mp_real) quadts2d, fun, tsum1, tsum2, tsum3, &
  s1, s2, s3, t1, t2, t3, t4, tw1, tw2, tw3, tw4, wij, wk(-1:nq2), xk(-1:nq2)
type (mp_realx) ax, bx, ay, by, f1, f2, fun1, fun2, xki, xx1, xx2, y1, y2, &
  ykj, yy1, yy2
integer jz1(0:nq2), jz2(0:nq2), jz3(0:nq2), jz4(0:nq2)
external fun, fun1, fun2, dplog10q

call mpsetprecwords (nwords2)
ay = 0.5d0 * (y2 - y1)
by = 0.5d0 * (y2 + y1)
call mpsetprecwords (nwords1)

if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadts2d: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

k = lev
h = 0.5d0 ** k
k1 = ip(nqq1-k)
k2 = ip(nqq1-k+1)
k3 = ip(nqq1-k+2)
dwmx = -99999.d0
dwimx = -99999.d0
dwjmx = -99999.d9
dwi1 = -99999.d0
dwi2 = -99999.d0
dwi3 = -99999.d0
dwi4 = -99999.d0
tsum1 = 0.d0
tsum2 = 0.d0
tsum3 = 0.d0

do i = 0, n, k1
  dwj1(i) = -99999.d0
  dwj2(i) = -99999.d0
  dwj3(i) = -99999.d0
  dwj4(i) = -99999.d0
  jz1(i) = 0
  jz2(i) = 0
  jz3(i) = 0
  jz4(i) = 0
enddo

!   Evaluate function at level k in x and y, avoiding unnecessary computation.

do j = 0, n, k1

  if (mod (j, 100) == 0) write (6, *) j

  call mpsetprecwords (nwords2)
  ykj = 1.d0 - xk(j)
  yy1 = - ay * ykj + by
  yy2 = ay * ykj + by
!    logy1 = yy1 > y1
!    logy2 = yy2 < y2
  logy1 = .true.
  logy2 = .true.
  f1 = fun1 (yy1)
  f2 = fun2 (yy1)
  ax = 0.5d0 * (f2 - f1)
  bx = 0.5d0 * (f2 + f1)
  call mpsetprecwords (nwords1)
  iz1 = 0
  iz2 = 0
  iz3 = 0
  iz4 = 0

  do i = 0, n, k1
    call mpsetprecwords (nwords2)
    xki = 1.d0 - xk(i)
    xx1 = - ax * xki + bx
    xx2 = ax * xki + bx
!        logx1 = xx1 > f1
!        logx2 = xx2 < f2
    logx1 = .true.
    logx2 = .true.
    call mpsetprecwords (nwords1)
    wij = wk(i) * wk(j)

    if (logx1 .and. logy1 .and. iz1 < izx .and. jz1(i) < izx) then
      t1 = fun (xx1, yy1)
      tw1 = t1 * wij
      dwi1 = dplog10q (tw1)
      dwj1(i) = dwi1
      if (dwi1 < nepsilon1) then
        iz1 = iz1 + 1
        jz1(i) = jz1(i) + 1
      else
        iz1 = 0
        jz1(i) = 0
      endif
    else
      t1 = 0.d0
      tw1 = 0.d0
    endif

    if (i > 0 .and. logx2 .and. logy1 .and. iz2 < izx .and. jz2(i) < izx) &
      then
      t2 = fun (xx2, yy1)
      tw2 = t2 * wij
      dwi2 = dplog10q (tw2)
      dwj2(i) = dwi2
      if (dwi2 < nepsilon1) then
        iz2 = iz2 + 1
        jz1(i) = jz1(i) + 1
      else
        iz2 = 0
        jz2(i) = 0
      endif
    else
      t2 = 0.d0
      tw2 = 0.d0
    endif

    t1 = ax * ay * (tw1 + tw2)
    tsum1 = tsum1 + t1
    if (mod (i, k2) == 0 .and. mod (j,k2) == 0) tsum2 = tsum2 + t1
    if (mod (i, k3) == 0 .and. mod (j,k3) == 0) tsum3 = tsum3 + t1
    dwmx = max (dwmx, dplog10q (tw1), dplog10q (tw2))
  enddo

  call mpsetprecwords (nwords2)
  f1 = fun1 (yy2)
  f2 = fun2 (yy2)
  ax = 0.5d0 * (f2 - f1)
  bx = 0.5d0 * (f2 + f1)
  call mpsetprecwords (nwords1)

  do i = 0, n, k1
    call mpsetprecwords(nwords2)
    xki = 1.d0 - xk(i)
    xx1 = - ax * xki + bx
    xx2 = ax * xki + bx
!        logx1 = xx1 > f1
!        logx2 = xx2 < f2
    logx1 = .true.
    logx2 = .true.
    call mpsetprecwords (nwords1)
    wij = wk(i) * wk(j)

    if (j > 0 .and. logx1 .and. logy2 .and. iz3 < izx .and. jz3(i) < izx) &
      then
      t3 = fun (xx1, yy2)
      tw3 = t3 * wij
      dwi3 = dplog10q (tw3)
      dwj3(i) = dwi3
      if (dwi3 < nepsilon1) then
        iz3 = iz3 + 1
        jz3(i) = jz3(i) + 1
      else
        iz3 = 0
        jz3(i) = 0
      endif
    else
      t3 = 0.d0
      tw3 = 0.d0
    endif

    if (i > 0 .and. j > 0 .and. logx2 .and. logy2 .and. iz4 < izx &
      .and. jz4(i) < izx) then
      t4 = fun (xx2, yy2)
      tw4 = t4 * wij
      dwi4 = dplog10q (tw4)
      dwj4(i) = dwi4
      if (dwi4 < nepsilon1) then
        iz4 = iz4 + 1
        jz4(i) = jz4(i) + 1
      else
        iz4 = 0
        jz4(i) = 0
      endif
    else
      t4 = 0.d0
      tw4 = 0.d0
    endif

    t1 = ax * ay * (tw3 + tw4)
    tsum1 = tsum1 + t1
    if (mod (i, k2) == 0 .and. mod (j,k2) == 0) tsum2 = tsum2 + t1
    if (mod (i, k3) == 0 .and. mod (j,k3) == 0) tsum3 = tsum3 + t1
    dwmx = max (dwmx, dplog10q (tw3), dplog10q (tw4))
  enddo

  dwimx = max (dwimx, dwi1, dwi2, dwi3, dwi4)
enddo

do i = 0, n, k1
  dwjmx = max (dwjmx, dwj1(i), dwj2(i), dwj3(i), dwj4(i))
enddo

!   Compute s1 = current integral approximation and err = error estimate.

s1 = h**2 * tsum1
s2 = 4.d0 * h**2 * tsum2
s3 = 16.d0 * h**2 * tsum3
eps1 = dwmx + nepsilon1
eps2 = max (dwimx, dwjmx)
d1 = dplog10q (abs (s1 - s2))
d2 = dplog10q (abs (s1 - s3))
d3 = eps1 - 1.d0
d4 = eps2 - 1.d0

if (k <= 2) then
  err = 0.d0
  elseif (d1 .eq. -99999.d0) then
  err = -99999.d0
else
  err = nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
endif

!   Output current integral approximation and error estimate, to 56 dp.

if (ndebug >= 2) then
  write (6, 2) nint (err)
2 format ('quadts2d: est error = 10^',i5,'; approx value =')
  call mpwrite (6, s1)
endif
if (k >= 3 .and. err < eps1) goto 140
if (k >= 3 .and. err < eps2) goto 120

write (6, 3) nint (err), nquadl
3 format ('quadts2d: Estimated error = 10^',i5/&
    'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
goto 140

120 continue

write (6, 4) nint (err), ndigits1
4 format ('quadts2d: Estimated error = 10^',i5/&
   'Increase working prec (Digits) for greater accuracy. Current Digits =',i4)
goto 140

130 continue

if (ierror > 0) nerror = ierror + 100
  write (6, 5) nerror
5 format ('quadts2d: Error in quadrature calculation; code =',i5)
s1 = 0.d0

140 continue

quadts2d = s1
return
end

function catalan ()
use mpmodule
use quadglobal
implicit none
integer k
real*8 dk, eps
type (mp_real) catalan, c1, c2, c4, c8, r16, t1, t2, t3
type (mp_real) x1, x2, x3, x4, x5, x6

c1 = 1.d0
c2 = 2.d0
c4 = 4.d0
c8 = 8.d0
r16 = 1.d0 / 16.d0
t1 = 0.d0
t2 = 1.d0
eps = mpreal (10.d0) ** nepsilon1

do k = 0, 10000000
  dk = k
  t3 = t2 * (c8 / (8.d0 * dk + 1.d0) ** 2 + c8 / (8.d0 * dk + 2.d0) ** 2 &
       + c4 / (8.d0 * dk + 3.d0) ** 2 - c2 / (8.d0 * dk + 5.d0) ** 2 &
       - c2 / (8.d0 * dk + 6.d0) ** 2 - c1 / (8.d0 * dk + 7.d0) ** 2)
  t1 = t1 + t3
  t2 = r16 * t2
  if (t3 < 1.d-5 * eps) goto 100
enddo

write (6, *) 'catalan: error - contact author'

100 continue

catalan = 1.d0 / 8.d0 * mppic * log (c2) + 1.d0 / 16.d0 * t1
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (mp_real) a

! call mpmdc (a%mpr, da, ia)
da = a
ia = 0
if (da .eq. 0.d0) then
  dplog10q = -99999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmdq (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

! call mpmdc (a%mpr, da, ia)
da = a
ia = 0
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end
