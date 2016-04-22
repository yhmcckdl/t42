! 一维欧拉方程数值解 (2016-04-22)  --by Chen Chen 
! 激波管问题&激波-密度扰动波干扰问题求解程序
! Scheme：5阶WENO，2阶GVC
! FVS：Steger-Warming

 PROGRAM RiemanNumeriComput
 IMPLICIT NONE
 
 INTEGER i, j, k, iterN, points, scheme, prob
 REAL(8) dt, dx, time, t, length, time_end, time_step, n, x1, x2
 REAL(8) UL, DL, PL, UR, DR, PR, u1, d1, p1, u2, p2
 REAL(8),ALLOCATABLE,DIMENSION(:,:) ::  Fd, U, Q           ! U (D,DU,E), Q (D,U,P)

 COMMON /G1/ dt, dx, t, x1
 COMMON /G2/ points, scheme
 COMMON /G3/ UL, DL, PL, UR, DR, PR
 COMMON /G4/ u1, d1, p1, u2, p2


!===================激波管初始条件================
 UL = 0.0   ;  DL = 1.0    ;  PL = 1.0
 UR = 0.0   ;  DR = 0.125  ;  PR = 0.1
!============激波-密度扰动波干扰初始条件=========
 u1 = 2.629 ;  d1 = 3.857  ;  p1 = 10.333
 u2 = 0.0   ;  p2 = 1.0
!======================计算条件===================
 time    = 0.20  ;  x1     = 0.0      ;  x2 = 1.0
 points  = 2000  ;  dt     = 0.000005
 prob    = 1     ;  scheme = 2                     ! prob1 = sod   prob2 = homework 5-2   scheme1 = GVC2   scheme2 = WENO5
 length  = x2 - x1
 dx      = length / points
 iterN   = time / dt

 ALLOCATE(Fd(points,3))
 ALLOCATE(Q(points,3))
 ALLOCATE(U(points,3))

 CALL MESH_GENERATION(prob)     ! 生成网格     
 CALL READ_MESH(Q)              ! 读网格
 CALL CONSER_VARI_GET(Q, U)     ! 计算守恒变量
 
 t = 0.0
 
 DO i = 1, iterN
  CALL STEG_WARM(U, Fd)         ! Steger-Warming
  CALL RK_THREE(U, Fd)          ! 3阶Runge-Kutta
  t = t + dt

 CALL CPU_TIME(time_step)

 n = (t/time)*100.0
 WRITE(*,'(A25, F5.2, A2, A18,F8.5,A13,F6.2,A4)')'Calculation Progress : ',n,'%', '   At Time Step =', t, &
                                                '   CPU Time = ',time_step,'sec'
 END DO

 CALL OUTPUT(Q, U)              ! 输出文件
 CALL CPU_TIME(time_end)        ! CPU_TIME 返回本程序执行到此语句所用时间sec
 WRITE(*,'(A19,F6.2,A4)') 'Total CPU time =', time_end, 'sec'

 END


! ==============================Steger-Warming分裂==========================
! 功能：提供n时刻的守恒变量U矩阵，返回n时刻“U矩阵对应的”通量差分矩阵Fdj
! 特点：只与本时刻的值有关
! 步骤：1）计算原始变量；               2）1个点1个点处理，写出F1,F2；
!       3）调用差分格式计算F1d, F2d；   4）最后通量差分Fd = F1d + F2d
! ==========================================================================
 SUBROUTINE STEG_WARM(U, Fd)
 IMPLICIT NONE
 REAL(8),ALLOCATABLE,DIMENSION(:,:) :: Q, F1, F2, F1d, F2d
 REAL(8) :: LAMBDA(3), LAMBDA_P(3), LAMBDA_M(3)
 REAL(8) a, dt, dx, x1, t, P1, P2, P3, WP, M1, M2, M3, WM    ! 通量分裂中第三个分量的计算系数，见课件p37
 INTEGER DIREC, i, j, points, scheme, k
 COMMON /G1/ dt, dx, t, x1
 COMMON /G2/ points, scheme
 REAL(8) :: U(points,3), Fd(points,3)

 ALLOCATE(F1(points,3))
 ALLOCATE(F2(points,3))
 ALLOCATE(F1d(points,3))
 ALLOCATE(F2d(points,3))
 ALLOCATE(Q(points,3))


! 此循环计算F1 F2
! ==================
 DO i = 1, points
  ! 原始变量
  Q(i,1)     = U(i,1)                                                     ! 密度
  Q(i,2)     = U(i,2) / U(i,1)                                            ! 速度
  Q(i,3)     = 0.4 * ( U(i,3) - 0.5 * ( ( U(i,2) * U(i,2) ) / U(i,1) ) )  ! 压强 See Riemann book p89
  a          = SQRT( (1.4 * Q(i,3) ) / Q(i,1) )                           ! 声速
  ! 3个特征值
  LAMBDA(1)  = Q(i,2)
  LAMBDA(2)  = Q(i,2) - a
  LAMBDA(3)  = Q(i,2) + a
  ! 3个正负特征值
  DO j = 1, 3
   LAMBDA_P(j) = 0.5 * ( LAMBDA(j) + SQRT( ( LAMBDA(j) * LAMBDA(j) ) + 0.0000000001**2) )
   LAMBDA_M(j) = 0.5 * ( LAMBDA(j) - SQRT( ( LAMBDA(j) * LAMBDA(j) ) + 0.0000000001**2) )
  END DO
  ! F1 的3个分量
  F1(i,1)    = ( Q(i,1) / 2.8 ) * ( 2 * 0.4 * LAMBDA_P(1) + LAMBDA_P(2) + LAMBDA_P(3) )                   
  F1(i,2)    = ( Q(i,1) / 2.8 ) * ( 2 * 0.4 * LAMBDA_P(1) * Q(i,2) + LAMBDA_P(2) * LAMBDA(2) + LAMBDA_P(3) * LAMBDA(3) )
  P1         = 0.4 * LAMBDA_P(1) * Q(i,2) * Q(i,2)       ;  P2 = 0.5 * LAMBDA_P(2) * LAMBDA(2) * LAMBDA(2)
  P3         = 0.5 * LAMBDA_P(3) * LAMBDA(3) * LAMBDA(3) ;  WP = 2.0 * ( LAMBDA_P(2) + LAMBDA_P(3) ) * a * a
  
  F1(i,3)    = ( Q(i,1) / 2.8 ) * ( P1 + P2 + P3 + WP )
  ! F2 的3个分量
  F2(i,1)    = ( Q(i,1) / 2.8 ) * ( 2 * 0.4 * LAMBDA_M(1) + LAMBDA_M(2) + LAMBDA_M(3) )                   
  F2(i,2)    = ( Q(i,1) / 2.8 ) * ( 2 * 0.4 * LAMBDA_M(1) * Q(i,2) + LAMBDA_M(2) * LAMBDA(2) + LAMBDA_M(3) * LAMBDA(3) )
  M1         = 0.4 * LAMBDA_M(1) * Q(i,2) * Q(i,2)       ;  M2 = 0.5 * LAMBDA_M(2) * LAMBDA(2) * LAMBDA(2)
  M3         = 0.5 * LAMBDA_M(3) * LAMBDA(3) * LAMBDA(3) ;  WM = 2.0 * ( LAMBDA_M(2) + LAMBDA_M(3) ) * a * a
  
  F2(i,3)    = ( Q(i,1) / 2.8 ) * ( M1 + M2 + M3 + WM )
 
 END DO

! 计算差分F1 F2, F = F1 + F2   后差1 前差2
! ==================
 DO i = 1, 3            ! 按列计算差分
  IF(scheme == 1)THEN   ! GVC
    DIREC = 1
    CALL GVC2(DIREC, F1(:,i), F1d(:,i))
    DIREC = 2
    CALL GVC2(DIREC, F2(:,i), F2d(:,i))

  ELSE IF(scheme == 2)THEN  ! WENO
    DIREC = 1
    CALL WENO5(DIREC, F1(:,i), F1d(:,i))
    DIREC = 2
    CALL WENO5(DIREC, F2(:,i), F2d(:,i))
  END IF
   Fd(:,i) = F1d(:,i) + F2d(:,i)   ! 合并
!   WRITE(*,*)Fd(:,1) ; STOP
 END DO

 RETURN
 END SUBROUTINE STEG_WARM


! ============================2阶GVC格式==============================
! 功能：提供n时刻列向量Fk（维数为全局变量points）和迎风方向DIREC，
!       返回n时刻的导数列向量Fkd
! 特点：守恒型，提供边界条件处理
! 具体格式表达式见李新亮课件p35
! ====================================================================

 SUBROUTINE GVC2(DIREC, Fk, Fkd)
 IMPLICIT NONE
 REAL(8) JUDGE, dx, dt, t, x1, F2, F1  ! F2 = Fj + 0.5, F1 = Fj - 0.5
 INTEGER DIREC, points, scheme, i
 COMMON /G2/ points, scheme
 COMMON /G1/ dt, dx, t, x1

 REAL(8) :: Fk(points), Fkd(points)
 
 SELECT CASE(DIREC)
 CASE (1) ! 正通量

  Fkd(1) = 0.0 ; Fkd(2) = 0.0 ; Fkd(points-1) = 0.0 ; Fkd(points) = 0.0  ! 边界处理

  DO i = 3, (points-2)                                                   ! 只循环计算非边界点的导数值
   JUDGE = ABS( Fk(i) - Fk(i-1) ) - ABS( Fk(i+1) - Fk(i) )               ! 每次重构前先判断i点位于波前OR波后
   IF(JUDGE < 0.0) THEN
    F2 = 0.5 * ( 3.0 * Fk(i) - Fk(i-1) )
    F1 = 0.5 * ( 3.0 * Fk(i-1) - Fk(i-2) )                               ! 先分别重构每个点的左值和右值，再计算差分
   ELSE
    F2 = 0.5 * ( Fk(i) + Fk(i+1) )
    F1 = 0.5 * ( Fk(i-1) + Fk(i) )
   END IF
    Fkd(i) = ( F2 - F1 ) / dx                                            ! 守恒型格式
  END DO

 CASE (2) ! 负通量

  Fkd(1) = 0.0 ; Fkd(2) = 0.0 ; Fkd(points-1) = 0.0 ; Fkd(points) = 0.0

  DO i = 3, (points-2)
   JUDGE = ABS( Fk(i+1) - Fk(i+2) ) - ABS( Fk(i) - Fk(i+1) )
   IF(JUDGE < 0.0) THEN
    F2 = 0.5 * ( 3.0 * Fk(i+1) - Fk(i+2) )
    F1 = 0.5 * ( 3.0 * Fk(i) - Fk(i+1) )
   ELSE
    F2 = 0.5 * ( Fk(i) + Fk(i+1) )
    F1 = 0.5 * ( Fk(i-1) + Fk(i) )
   END IF
    Fkd(i) = ( F2 - F1 ) / dx
  END DO

 END SELECT

 RETURN
 END SUBROUTINE GVC2

! =============================5阶WENO格式============================
! 功能：提供n时刻列向量Fk（维数为全局变量points）和迎风方向DIREC，
!       返回n时刻的导数列向量Fkd
! 特点：守恒型，提供边界条件处理
! 具体格式表达式见李新亮课件p47
! ====================================================================
 SUBROUTINE WENO5(DIREC, Fk, Fkd)
 IMPLICIT NONE
 REAL(8) dx, dt, t, w1(points+1), w2(points+1), w3(points+1)
 REAL(8) IS1, IS2, IS3, a1, a2, a3, a, x1
 INTEGER DIREC, points, scheme, i

 COMMON /G2/ points, scheme
 COMMON /G1/ dt, dx, t, x1
 
 REAL(8) :: Fk(points), Fkd(points), F(points+5)
 REAL(8) :: Fk1(points), Fk2(points), Fk3(points), Fj(points+1)

 Fj = 0.0

 SELECT CASE(DIREC)
 CASE (1) ! 正通量
  ! 构造实际计算数组F，左侧3个虚网格，右侧2个虚网格
  F(1) = 0.0 ; F(2) = 0.0 ; F(3) = 0.0 ; F(points+4) = 0.0 ; F(points+5) = 0.0

  DO i = 1, points
   F(i+3) = Fk(i)
  END DO

  DO i = 1, (points+1)
   Fk1(i) = (1.0/3.0) * F(i) - (7.0/6.0) * F(i+1) + (11.0/6.0) * F(i+2)
   Fk2(i) = (-1.0/6.0) * F(i+1) + (5.0/6.0) * F(i+2) + (1.0/3.0) * F(i+3)
   Fk3(i) = (1.0/3.0) * F(i+2) + (5.0/6.0) * F(i+3) - (1.0/6.0) * F(i+4)
   IS1    = 0.25 * ( F(i) - 4.0 * F(i+1) + 3.0 * F(i+2) )**2.0 + (13.0/12.0) * ( F(i) - 2.0 * F(i+1) + F(i+2) )**2.0
   IS2    = 0.25 * ( F(i+1) - F(i+3) )**2.0 + (13.0/12.0) * ( F(i+1) - 2.0 * F(i+2) + F(i+3) )**2.0
   IS3    = 0.25 * ( 3.0 * F(i+2) - 4.0 * F(i+3) + F(i+4) )**2.0 + (13.0/12.0) * ( F(i+2) - 2.0 * F(i+3) + F(i+4) )**2.0
   a1     = 0.1 / ( ( 0.000001 + IS1 )**2.0 )
   a2     = 0.6 / ( ( 0.000001 + IS2 )**2.0 )
   a3     = 0.3 / ( ( 0.000001 + IS3 )**2.0 )
   a      = a1 + a2 + a3
   w1(i)  = a1 / a ; w2(i) = a2 / a ; w3(i) = a3 / a
  END DO

  DO i = 4, (points-1)
   Fj(i)  = w1(i) * Fk1(i) + w2(i) * Fk2(i) + w3(i) * Fk3(i)
  END DO

  DO i = 4, (points-2)
   Fkd(i)  = ( Fj(i+1) - Fj(i) )  / dx
  END DO

  ! 边界处理
  ! j = 1            ! j = 2 
  Fkd(1) = 0.0 ;     Fkd(2) = ( Fk3(3) - Fk3(2) ) / dx

  ! j = 3
  DO i = 3, 4
   IS2 = 0.25 * ( F(i+1) - F(i+3) )**2.0 + (13.0/12.0) * ( F(i+1) - 2.0 * F(i+2) + F(i+3) )**2.0
   IS3 = 0.25 * ( 3.0 * F(i+2) - 4.0 * F(i+3) + F(i+4) )**2.0 + (13.0/12.0) * ( F(i+2) - 2.0 * F(i+3) + F(i+4) )**2.0
   a2  = 0.6 / ( ( 0.000001 + IS2 )**2.0 )
   a3  = 0.3 / ( ( 0.000001 + IS3 )**2.0 )
   a   = a2 + a3
   w2(i) = a2 / a ; w3(i) = a3 / a
  END DO
  Fkd(3) = ( w2(4)*Fk2(4) + w3(4)*Fk3(4) - w2(3)*Fk2(3) - w3(3)*Fk3(3) ) / dx

  ! j = points-1
  DO i = (points-1), points
   IS1    = 0.25 * ( F(i) - 4.0 * F(i+1) + 3.0 * F(i+2) )**2.0 + (13.0/12.0) * ( F(i) - 2.0 * F(i+1) + F(i+2) )**2.0
   IS2    = 0.25 * ( F(i+1) - F(i+3) )**2.0 + (13.0/12.0) * ( F(i+1) - 2.0 * F(i+2) + F(i+3) )**2.0
   a1  = 0.6 / ( ( 0.000001 + IS1 )**2.0 )
   a2  = 0.3 / ( ( 0.000001 + IS2 )**2.0 )
   a   = a1 + a2
   w1(i) = a1 / a ; w2(i) = a2 / a
  END DO
  Fkd(points-1) = (w1(points)*Fk1(points)+w2(points)*Fk2(points)-w1(points-1)*Fk1(points-1)-w2(points-1)*Fk2(points-1))/dx

  ! j = points
  Fkd(points) = ( Fk1(points+1) - Fk1(points) ) / dx

 CASE (2) ! 负通量
  F(1) = 0.0 ; F(2) = 0.0 ; F(points+3) = 0.0 ; F(points+4) = 0.0 ; F(points+5) = 0.0
  DO i = 1, points
   F(i+2) = Fk(i)
  END DO

  DO i = 1, points+1
   Fk1(i) = (1.0/3.0) * F(i+4) - (7.0/6.0) * F(i+3) + (11.0/6.0) * F(i+2)
   Fk2(i) = (-1.0/6.0) * F(i+3) + (5.0/6.0) * F(i+2) + (1.0/3.0) * F(i+1)
   Fk3(i) = (1.0/3.0) * F(i+2) + (5.0/6.0) * F(i+1) - (1.0/6.0) * F(i) 
   IS1    = 0.25 * ( F(i+4) - 4.0 * F(i+3) + 3.0 * F(i+2) )**2.0 + (13.0/12.0) * ( F(i+4) - 2.0 * F(i+3) + F(i+2) )**2.0
   IS2    = 0.25 * ( F(i+3) - F(i+1) )**2.0 + (13.0/12.0) * ( F(i+3) - 2.0 * F(i+2) + F(i+1) )**2.0
   IS3    = 0.25 * ( 3.0 * F(i+2) - 4.0 * F(i+1) + F(i) )**2.0 + (13.0/12.0) * ( F(i+2) - 2.0 * F(i+1) + F(i) )**2.0
   a1     = 0.1 / ( ( 0.000001 + IS1 )**2.0 )
   a2     = 0.6 / ( ( 0.000001 + IS2 )**2.0 )
   a3     = 0.3 / ( ( 0.000001 + IS3 )**2.0 )
   a      = a1 + a2 + a3
   w1(i)  = a1 / a ; w2(i) = a2 / a ; w3(i) = a3 / a
  END DO

  DO i = 3, (points-2)
   Fj(i)  = w1(i) * Fk1(i) + w2(i) * Fk2(i) + w3(i) * Fk3(i)
  END DO

  DO i = 3, (points-3)
   Fkd(i)  = ( Fj(i+1) - Fj(i) )  / dx
  END DO

  ! 边界处理
  ! j = 1
  Fkd(1) = ( Fk1(2) - Fk1(1) ) / dx

  ! j = 2
  DO i = 2, 3
   IS1    = 0.25 * ( F(i) - 4.0 * F(i+1) + 3.0 * F(i+2) )**2.0 + (13.0/12.0) * ( F(i) - 2.0 * F(i+1) + F(i+2) )**2.0
   IS2    = 0.25 * ( F(i+1) - F(i+3) )**2.0 + (13.0/12.0) * ( F(i+1) - 2.0 * F(i+2) + F(i+3) )**2.0
   a1  = 0.6 / ( ( 0.000001 + IS1 )**2.0 )
   a2  = 0.3 / ( ( 0.000001 + IS2 )**2.0 )
   a   = a1 + a2
   w1(i) = a1 / a ; w2(i) = a2 / a
  END DO
  Fkd(2) = ( w1(3)*Fk1(3) + w2(3)*Fk2(3) - w1(2)*Fk1(2) - w2(2)*Fk2(2) ) / dx

  ! j = points        ! j = points-1
  Fkd(points) = 0.0 ; Fkd(points-1) = ( Fk3(points) - Fk3(points-1) ) / dx

  ! j = points-2
  DO i = (points-2), points-1
   IS2 = 0.25 * ( F(i+1) - F(i+3) )**2.0 + (13.0/12.0) * ( F(i+1) - 2.0 * F(i+2) + F(i+3) )**2.0
   IS3 = 0.25 * ( 3.0 * F(i+2) - 4.0 * F(i+3) + F(i+4) )**2.0 + (13.0/12.0) * ( F(i+2) - 2.0 * F(i+3) + F(i+4) )**2.0
   a2  = 0.6 / ( ( 0.000001 + IS2 )**2.0 )
   a3  = 0.3 / ( ( 0.000001 + IS3 )**2.0 )
   a   = a2 + a3
   w2(i) = a2 / a ; w3(i) = a3 / a
  END DO
  Fkd(points-2) = (w2(points-1)*Fk2(points-1)+w3(points-1)*Fk3(points-1)-w2(points-2)*Fk2(points-2)-w3(points-2)*Fk3(points-2))/dx

 END SELECT

 RETURN
 END SUBROUTINE WENO5


! ====================3阶Runge-Kutta方法==================
! 功能：提供 n 时刻的守恒变量矩阵 U 和通量差分矩阵 Fd，
!       然后“整体推进”到 n+dt 时刻的 U，以此更新 U
! ========================================================
 SUBROUTINE RK_THREE(u, fd)
 IMPLICIT NONE
 REAL(8),ALLOCATABLE,DIMENSION(:,:) :: u1, u2, u1d, u2d ,ud
 REAL(8) :: u(points,3), fd(points,3)
 REAL(8) dt, dx, t, x1
 INTEGER i, j, points, scheme

 COMMON /G1/ dt, dx, t, x1
 COMMON /G2/ points, scheme

 ALLOCATE(ud(points,3))
 ALLOCATE(u1(points,3))
 ALLOCATE(u2(points,3))
 ALLOCATE(u1d(points,3))
 ALLOCATE(u2d(points,3))

! u1矩阵
 DO i = 1, 3
  DO j = 1, points
   u1(j,i) = u(j,i) + dt * ( -fd(j,i) )
  END DO
 END DO

 CALL STEG_WARM(u1, fd)

! u2矩阵
 DO i = 1, 3
  DO j = 1, points
   u2(j,i) = 0.75 * u(j,i) + 0.25 * ( u1(j,i) + dt * ( -fd(j,i) ) )
  END DO
 END DO

 CALL STEG_WARM(u2, fd)

! 更新u矩阵
 DO i = 1, 3
  DO j = 1, points
   u(j,i) = (1.0/3.0) * u(j,i) + (2.0/3.0) * ( u2(j,i) + dt * ( -fd(j,i) ) )
  END DO
 END DO

 RETURN
 END SUBROUTINE RK_THREE



!=============================生成初始网格==============================
 SUBROUTINE MESH_GENERATION(prob)
 IMPLICIT NONE
 REAL(8) UL, DL, PL, UR, DR, PR, dt, dx, u1, d1, p1, u2, p2, x, x1, t
 INTEGER points, i, scheme, prob

 COMMON /G1/ dt, dx, t, x1
 COMMON /G3/ UL, DL, PL, UR, DR, PR
 COMMON /G4/ u1, d1, p1, u2, p2
 COMMON /G2/ points, scheme

 OPEN(66,FILE='input.dat')

 SELECT CASE(prob)
 CASE (1) ! 生成求解sod激波管初始网格
  DO i = 1, (points/2)
   WRITE(66,*) DL, UL, PL
  END DO
  DO i = 1, (points/2)
   WRITE(66,*) DR, UR, PR
  END DO
 CASE (2) ! 生成求解激波-密度扰动波干扰问题的初始网格
  DO i = 1, (points/10)
   WRITE(66,*) d1, u1, p1
  END DO
  x = 0.0
  DO i = 1, (points-points/10)
   WRITE(66,*) 1 + 0.3*SIN(40.0*x), u2, p2
   x = x + dx
  END DO
 END SELECT

 CLOSE(66)

 RETURN

 END SUBROUTINE MESH_GENERATION


!==============读网格================
 SUBROUTINE READ_MESH(Q)
 IMPLICIT NONE
 INTEGER points, i, scheme
 COMMON /G2/ points, scheme
 REAL(8) ::  Q(points,3)

 OPEN(33,FILE='input.dat')
 DO i = 1, points 
  READ(33,*) Q(i,:)
 END DO
 
 CLOSE(33)

 RETURN
 END SUBROUTINE READ_MESH


!=================输出文件=================
 SUBROUTINE OUTPUT(Q, U)
 IMPLICIT NONE
 INTEGER points, i, scheme
 REAL(8) x, dx, t, dt, x1
 COMMON /G2/ points, scheme
 REAL(8) ::  Q(points,3), U(points,3)
 COMMON /G1/ dt, dx, t, x1
 CALL PRIM_VARI_GET(Q, U)  ! 计算原始变量

 OPEN(44,FILE='output.dat')
 x = x1
 DO i = 1, points
  WRITE(44,*) x, Q(i,:)
  x = x + dx
 END DO

 CLOSE(44)

 RETURN
 END SUBROUTINE OUTPUT


!============================计算守恒变量==============================
 SUBROUTINE CONSER_VARI_GET(Q, U)
 IMPLICIT NONE
 INTEGER points, i, scheme
 COMMON /G2/ points, scheme
 REAL(8) ::  Q(points,3), U(points,3)

 DO i = 1, points
  U(i,1) = Q(i,1)
  U(i,2) = Q(i,2) * Q(i,1)
  U(i,3) = Q(i,1) * ( 0.5 * Q(i,2) * Q(i,2) + ( Q(i,3) / ( 0.4 * Q(i,1) ) ) )
 END DO

 RETURN
 END SUBROUTINE CONSER_VARI_GET


!=============================计算原始变量================================
 SUBROUTINE PRIM_VARI_GET(Q, U)
 IMPLICIT NONE
 INTEGER points, i, scheme
 COMMON /G2/ points, scheme
 REAL(8) ::  Q(points,3), U(points,3)
 
 DO i = 1, points
  Q(i,1)     = U(i,1)                                                     ! 密度
  Q(i,2)     = U(i,2) / U(i,1)                                            ! 速度
  Q(i,3)     = 0.4 * ( U(i,3) - 0.5 * ( ( U(i,2) * U(i,2) ) / U(i,1) ) )  ! 压强 See Riemann book p89

 END DO

 RETURN
 END SUBROUTINE PRIM_VARI_GET




! 编程总结：
! 1、所有模块统一思考 但最好是一个模块搞定并调试完后，再写另一个模块
! 2、子程序数组传递时，在子程序中的变量命名时，当用到COMMON变量时，要注意声明变量的顺序
! 3、复制粘贴时要注意修改必要的参数
! 4、数据结构问题：有些变量是否需要定义数组？ 函数之间传递什么数据结构效率较高？
! 5、遇到代编号i的问题，最好先用特殊的数字1,2,3进行验证，再推广
! 6、运行出现“已放弃 (核心已转储)”提示可能是因为数组越界，或者F(i)=0.0，但i未指定
! 7、遇到未知问题时可以先编写小程序进行测试



