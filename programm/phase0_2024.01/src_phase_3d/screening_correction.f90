
subroutine Screening_Potential
  use m_Crystal_Structure,  only : altv,rltv
  use m_Control_Parameters, only : sw_screening_correction, screening_alpha, gmaxp
  use m_PlaneWaveBasisSet,    only : gr_l, ngabc, kgp
  use m_Parallelization,    only : ista_kngp,iend_kngp
  use m_Screening

  implicit none
  integer i

  real(8) ga,gb,gc,g(3)
  real(8) Lmin, Ltmp(3)

  !! PHASEと接続する場合には以下の文をアンコメントする
  Screening%La = position(altv(1,1),altv(2,1),altv(3,1)) ! PHASEと接続する場合
  Screening%Lb = position(altv(1,2),altv(2,2),altv(3,2)) ! PHASEと接続する場合
  Screening%Lc = position(altv(1,3),altv(2,3),altv(3,3)) ! PHASEと接続する場合
  ! カットオフ波数
  Screening%Gmax  = gmaxp      ! PHASEと接続する場合

  !---- PHASE変数kgp,ngabcのコピーを用意する。
  allocate( Screening%ngabc(kgp,3) )

  do i=1, kgp
     Screening%ngabc(i,1) = ngabc(i,1)
     Screening%ngabc(i,2) = ngabc(i,2)
     Screening%ngabc(i,3) = ngabc(i,3)
  end do

  Screening%Nk = kgp ! カットオフ半径内の波数ベクトルの総数

  ! 遮蔽パラメタ。次元は長さの逆数。単位は原子単位。[1/a.u.length]
  Ltmp(1) = length( Screening%La )
  Ltmp(2) = length( Screening%Lb )
  Ltmp(3) = length( Screening%Lc )
  Lmin = minval(Ltmp)
  Screening%alpha = screening_alpha/Lmin

  write(6,*) "screening correction..."
  write(6,*) "altv:"
  write(6,*) altv
  write(6,*) "La Lb Lc:"
  write(6,*) Screening%La, Screening%Lb, Screening%Lc
  write(6,*) "Gmax  alpha:"
  write(6,*) Screening%Gmax, Screening%alpha, screening_alpha

  ! 遮蔽ポテンシャルの設定
  call setup_screening
  write(6,*) "Nk:"
  write(6,*) Screening%Nk, ista_kngp,iend_kngp

  ! 遮蔽ポテンシャルの値のテスト表示
!  do i=1, Screening%Nk
!
!     ga = ngabc(i,1)
!     gb = ngabc(i,2)
!     gc = ngabc(i,3)
!     g(1) = rltv(1,1)*ga+rltv(1,2)*gb+rltv(1,3)*gc
!     g(2) = rltv(2,1)*ga+rltv(2,2)*gb+rltv(2,3)*gc
!     g(3) = rltv(3,1)*ga+rltv(3,2)*gb+rltv(3,3)*gc
!
!     write(6,'(i10,2e20.10)') i, Screening%phik(i), gr_l(i)
!     write(12,'(i10,6f11.6)') i, &
!          g, sqrt(dot_product(g,g)), &
!          gr_l(i), Screening%phik(i)
!  end do

end subroutine Screening_Potential

