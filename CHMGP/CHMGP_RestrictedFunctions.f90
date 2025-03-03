SUBROUTINE TankBlockCycle(S, x, n, is1stReservoir, Q, totQ, E, infil, infilPar, a1, a2, t1, t2)
IMPLICIT NONE 
INTEGER :: n
REAL(8) :: S(n), x(n), Q(n), infil(n), E(n)
REAL(8) :: a1, a2
REAL(8) :: t1, t2
REAL(8) :: infilPar, totQ

REAL(8) :: Q1, Q2

INTEGER :: i

INTEGER :: is1stReservoir

do i = 2, n 
 
  if(is1stReservoir > 0) then 
    if(x(i).le.0) then
     S(i) = max(S(i-1) - E(i) - totQ, 0._8)
    else
     S(i) = max(S(i-1) + x(i) - totQ, 0._8)
    endif
  else
    S(i) = max(S(i-1) + x(i-1) - totQ, 0._8)
  endif
  
  Q1 = 0._8
  if(S(i) > t1) then 
    Q1 = a1 * (S(i) - t1) 
  endif  
  
  Q2 = 0._8
  if(S(i) > t2) then 
    Q2 = a2 * (S(i) - t2) 
  endif 
  
  infil(i) = infilPar * S(i)

  Q(i) = Q1 + Q2
  totQ = Q(i) + infil(i)

enddo

END SUBROUTINE TankBlockCycle


SUBROUTINE R4T(x, n, t1, a1, t2, a2, a23, t3, a3, a34, t4, a4, a45, a5, output)
IMPLICIT NONE

INTEGER :: n
REAL(8) :: x(n), output(n)
REAL(8) :: a1(n), a2(n), a23(n), a3(n), a34(n), a4(n), a45(n), a5(n)
REAL(8) :: t1(n), t2(n), t3(n), t4(n)
REAL(8) :: Q1(n), Q2(n), Q23(n), Q3(n), Q34(n), Q4(n), Q45(n), Q5(n) 
REAL(8) :: S1(n), S2(n), S3(n), S4(n)
INTEGER :: i

S1(1) = 1._8
S2(1) = 1._8
S3(1) = 1._8
S4(1) = 1._8

! Tank 1
if(S1(1) > t1(1)) then
  Q1(1) = a1(1) * (S1(1) - t1(1)) 
else 
  Q1(1) = 0._8
end if

if(S1(1) > t2(1)) then
  Q2(1) = a2(1) * (S1(1) - t2(1)) 
else 
  Q2(1) = 0._8
end if

Q23(1) = max(a23(1) * S1(1), 0._8)


! Tank 2
if(S2(1) > t3(1)) then
  Q3(1) = a3(1) * (S2(1) - t3(1)) 
else 
  Q3(1) = 0._8
end if

Q34(1) = max(a34(1) * S2(1), 0._8)


! Tank 3
if(S3(1) > t4(1)) then
  Q4(1) = a4(1) * (S3(1) - t4(1)) 
else 
  Q4(1) = 0._8
end if

Q45(1) = max( a45(1) * S3(1), 0._8) 


! Tank 4
Q5(1) = max( a5(1) * S4(1), 0._8)


! Main loop
do i = 2, n 

  ! Tank 1  
  S1(i) = max(S1(i-1) + x(i) - Q1(i-1) - Q2(i-1) - Q23(i-1), 0._8)  
  
  if(S1(i) > t1(i)) then
    Q1(i) = a1(i) * (S1(i) - t1(i)) 
  else 
    Q1(i) = 0._8
  end if
  
  if(S1(i) > t2(i)) then
    Q2(i) = a2(i) * (S1(i) - t2(i)) 
  else 
    Q2(i) = 0._8
  end if
  
  Q23(i) = max( a23(i) * S1(i) , 0._8)
  
  
  ! Tank 2
  S2(i) = max(S2(i-1) + Q23(i-1) - Q3(i-1) - Q34(i-1), 0._8) 
  
  if(S2(i) > t3(i)) then
    Q3(i) = a3(i) * (S2(i) - t3(i)) 
  else 
    Q3(i) = 0._8
  end if
  
  Q34(i) = max( a34(i) * S2(i) , 0._8)
  
  ! Tank 3
  S3(i) = max(S3(i-1) + Q34(i-1) - Q4(i-1) - Q45(i-1), 0._8)
  
  if(S3(i) > t4(i)) then
    Q4(i) = a4(i) * (S3(i) - t4(i)) 
  else 
    Q4(i) = 0._8
  end if
  
  Q45(i) = max( a45(i) * S3(i), 0._8) 
 
  ! Tank 4
  S4(i) = max(S4(i-1) + Q45(i-1) - Q5(i-1), 0._8) 
  
  Q5(i) = max( a5(i) * S4(i) , 0._8)
 
enddo

output = Q1 + Q2 + Q3 + Q4 + Q5

END SUBROUTINE R4T


SUBROUTINE R2T(x, n, t1, a1, t2, a2, a23, a3, output)
IMPLICIT NONE

INTEGER :: n
REAL(8) :: x(n), output(n)
REAL(8) :: a1(n), a2(n), a23(n), a3(n)
REAL(8) :: t1(n), t2(n)
REAL(8) :: Q1(n), Q2(n), Q23(n), Q3(n)
REAL(8) :: S1(n), S2(n)
INTEGER :: i

S1(1) = 1._8
S2(1) = 1._8

! Tank 1
if(S1(1) > t1(1)) then
  Q1(1) = a1(1) * (S1(1) - t1(1)) 
else 
  Q1(1) = 0._8
end if

if(S1(1) > t2(1)) then
  Q2(1) = a2(1) * (S1(1) - t2(1)) 
else 
  Q2(1) = 0._8
end if

Q23(1) = max(a23(1) * S1(1), 0._8)


! Tank 2
Q3(1) = max(a3(1) * S2(1), 0._8)

do i = 2, n 

 ! Tank 1  
  S1(i) = max(S1(i-1) + x(i) - Q1(i-1) - Q2(i-1) - Q23(i-1), 0._8)  
  
   if(S1(i) > t1(i)) then
    Q1(i) = a1(i) * (S1(i) - t1(i)) 
  else 
    Q1(i) = 0._8
  end if
  
  if(S1(i) > t2(i)) then
    Q2(i) = a2(i) * (S1(i) - t2(i)) 
  else 
    Q2(i) = 0._8
  end if
  
  Q23(i) = max( a23(i) * S1(i) , 0._8)
  
  ! Tank 2
  S2(i) = max(S2(i-1) + Q23(i-1) - Q3(i-1) , 0._8) 
  
  Q3(i) = max( a3(i) * S2(i) , 0._8)
  
enddo

output = Q1 + Q2 + Q3

END SUBROUTINE R2T
