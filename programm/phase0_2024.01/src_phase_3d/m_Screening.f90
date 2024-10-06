!*********************************************************
!* 逆空間でのクーロン力を遮蔽するポテンシャルを作成するモジュール
!*
!* 2006年6月19日改良版
!*********************************************************

! 3次元座標型とその演算子の定義
module class_position
  use m_Const_Parameters, only : DP
  implicit none

  type position
     real(kind=DP) :: x, y, z
  end type

  interface operator(+)
     module procedure add
  end interface

  interface operator(-)
     module procedure sub
  end interface

  interface operator(*)
     module procedure smul
  end interface

  interface operator(/)
     module procedure sdiv
  end interface

  interface operator(*)
     module procedure iprod
  end interface

  interface operator(.IPROD.)
     module procedure iprod
  end interface

  interface operator(.OPROD.)
     module procedure oprod
  end interface

contains
  function add( r1, r2 ) result(r)
    type(position), intent(in) :: r1, r2
    type(position) :: r

    r%x = r1%x + r2%x
    r%y = r1%y + r2%y
    r%z = r1%z + r2%z
  end function

  function sub( r1, r2 ) result(r)
    type(position), intent(in) :: r1, r2
    type(position) :: r

    r%x = r1%x - r2%x
    r%y = r1%y - r2%y
    r%z = r1%z - r2%z
  end function

  function smul( r1, d ) result(r)
    type(position), intent(in) :: r1
    real(kind=DP), intent(in) :: d
    type(position) :: r

    r%x = r1%x * d
    r%y = r1%y * d
    r%z = r1%z * d
  end function

  function sdiv( r1, d ) result(r)
    type(position), intent(in) :: r1
    real(kind=DP), intent(in) :: d
    type(position) :: r

    r%x = r1%x / d
    r%y = r1%y / d
    r%z = r1%z / d
  end function

  function iprod( r1, r2 ) result(d)
    type(position), intent(in) :: r1, r2
    real(kind=DP) :: d

    d = r1%x*r2%x + r1%y*r2%y + r1%z*r2%z
  end function

  function oprod( r1, r2 ) result(r)
    type(position), intent(in) :: r1, r2
    type(position) :: r

    r%x = r1%y*r2%z - r1%z*r2%y
    r%y = r1%z*r2%x - r1%x*r2%z
    r%z = r1%x*r2%y - r1%y*r2%x
  end function

  function length( r ) result(d)
    type(position), intent(in) :: r
    real(kind=DP) :: d

    d = sqrt(r%x*r%x + r%y*r%y + r%z*r%z)
  end function

  function length2( r ) result(d)
    type(position), intent(in) :: r
    real(kind=DP) :: d

    d = (r%x*r%x + r%y*r%y + r%z*r%z)
  end function

end module

! 遮蔽ポテンシャルモジュール
module m_Screening
  use m_Const_Parameters, only : DP, PAI
!  use m_PlaneWaveBasisSet, only : ngabc, kgp
  use m_Screening_FFT  ! 簡略化したPHASEのFFTモジュールを利用する
  use class_position
  implicit none

  ! 数学定数の設定
  real(kind=DP), parameter :: M_PI       = PAI ! 円周率
  real(kind=DP), parameter :: M_2_SQRTPI = 1.12837916709551257390d0 ! 2/sqrt(pi)

  ! 本モジュールで使用する変数群
  type Screening_param
     type(position) La, Lb, Lc ! 単位胞の格子ベクトル

     !---- PHASE変数のコピーを自前で持つ
     real(kind=DP), pointer, dimension(:,:) ::  ngabc ! PHASE変数

     real(kind=DP) Gmax ! カットオフ波数
     real(kind=DP) alpha ! 遮蔽ポテンシャルの遮蔽パラメタ
     real(kind=DP), dimension(:), pointer :: phik ! 遮蔽ポテンシャル
     integer Nk    ! 波数ベクトルの総数
  end type

  type(Screening_param) :: Screening

contains

  ! 遮蔽ポテンシャルを設定するサブルーチン
  subroutine setup_screening
    real(kind=DP)  V             ! 単位胞の体積
    type(position) dKa, dKb, dKc ! 逆格子の刻み
    integer        Na0, Nb0, Nc0 ! 格子点数
    integer        Na, Nb, Nc ! 格子点数
    type(position) dLa, dLb, dLc ! 実格子の刻み
    real(kind=DP)  dV            ! 実格子の刻み
    real(kind=DP)  r             ! 実空間の距離
    real(kind=DP)  kk            ! 逆格子ベクトルのノルム

    integer ia, ib, ic     ! 3次元格子点のインデックス
    integer ja, jb, jc     ! 3次元格子点のインデックス
    integer i, j           ! 1次元インデックス

    integer, dimension(3) :: n_rGpv ! FFTモジュールの初期化用の変数

    real(kind=DP),  allocatable, dimension(:) :: phir ! 実格子での値
    type(position), allocatable, dimension(:) :: K    ! 逆格子ベクトル

    ! 体積の計算
    V  = Screening%La .IPROD. ( Screening%Lb .OPROD. Screening%Lc )

    ! 逆格子の刻み幅の設定
    dKa = (Screening%Lb .OPROD. Screening%Lc) * (2.0d0*M_PI/V)
    dKb = (Screening%Lc .OPROD. Screening%La) * (2.0d0*M_PI/V)
    dKc = (Screening%La .OPROD. Screening%Lb) * (2.0d0*M_PI/V)

    ! PHASE変数の n_rGpv (逆格子点の数の概算値) と同じ値を自分で計算する。
    n_rGpv(1) = abs(Screening%Gmax/length(dKa)) + 1
    n_rGpv(2) = abs(Screening%Gmax/length(dKb)) + 1
    n_rGpv(3) = abs(Screening%Gmax/length(dKc)) + 1

    ! PHASE変数の fft_box_size_CD (逆格子点の数) と同じ値を自分で計算する。
    call m_Screening_FFT_set_box_sizes( n_rGpv, 1 )
    call m_Screening_FFT_setup( 1, .false. )

    ! 格子点の数の入手
    Na0 = fft_box_size_CD(1,0)
    Nb0 = fft_box_size_CD(2,0)
    Nc0 = fft_box_size_CD(3,0)
    Na = fft_box_size_CD(1,1)
    Nb = fft_box_size_CD(2,1)
    Nc = fft_box_size_CD(3,1)
    ! 実格子の刻み幅の設定
    dLa = Screening%La*(1.0d0/Na)
    dLb = Screening%Lb*(1.0d0/Nb)
    dLc = Screening%Lc*(1.0d0/Nc)
    dV  = dLa .IPROD. ( dLb .OPROD. dLc) ! ひとつの刻みの体積

    ! 実空間での遮蔽ポテンシャルの長距離成分を保存する配列変数を確保する
    allocate( phir(Na0*Nb0*Nc0) )

    ! 実空間での遮蔽ポテンシャルの長距離成分を設定するループ
    do ic = 0, Nc-1  !実空間の格子点でのループ
       do ib = 0, Nb-1
          do ia = 0, Na-1
             ! 単位胞の8隅がポテンシャルの原点となるように座標を調整する
             if( 2*ia>Na ) then
                ja = ia-Na
             else
                ja = ia
             end if
             if( 2*ib>Nb ) then
                jb = ib-Nb
             else
                jb = ib
             end if
             if( 2*ic>Nc ) then
                jc = ic-Nc
             else
                jc = ic
             end if

             ! 最寄の隅からの距離
             r = length(dLa*dble(ja) + dLb*dble(jb) + dLc*dble(jc))

             ! 3次元の格子インデックスから1次元配列のインデックスへの変換
             i = (ic*Nb0+ib)*Na0 + (ia + 1)

             ! 実空間での遮蔽ポテンシャルの長距離成分の設定
             phir(i) = phirlong(r)

!             write(15,'(i5,2e15.5,6i5)') i,phir(i),r,ia,ib,ic,ja,jb,jc
!             write(17,'(5f10.6)') dLa%x*dble(ja), dLb%y*dble(jb), dLc%z*dble(jc), r, phir(i)

          end do
       end do
    end do

    ! 実空間での遮蔽ポテンシャルの長距離成分をフーリエ変換で逆空間に変換する
    call m_Screening_FFT_set_box_sizes( n_rGpv, 1 )
    call m_Screening_FFT_setup( 1, .false. )
    call m_Screening_FFT_alloc_CD_box
    call m_Screening_FFT_CD_inverse( 6, phir )
    call m_Screening_FFT_dealloc_CD_box


    ! 逆空間での遮蔽ポテンシャルを保存する配列変数を確保する
    allocate( K(Screening%Nk) )
    allocate( Screening%phik(Screening%Nk) )

    do j=1, Screening%Nk
       ja=Screening%ngabc(j,1)
       jb=Screening%ngabc(j,2)
       jc=Screening%ngabc(j,3)
       ! 逆格子点の波数ベクトルのノルムを計算する
       K(j) = dKa*dble(ja) + dKb*dble(jb) + dKc*dble(jc)
       kk = length2(K(j))

       ! 逆格子点のインデックスを少し変換する。
       if( jc < 0 ) then
          ic = jc + Nc
       else
          ic = jc
       end if
       if( jb < 0 ) then
          ib = jb + Nb
       else
          ib = jb
       end if
       if( ja < 0 ) then
          ia = -2*ja
       else
          ia = +2*ja
       end if
       ! 3次元の格子インデックスから1次元配列のインデックスへの変換
       i = (ic*Nb0+ib)*Na0 + (ia + 1)

       Screening%phik(j) = phir(i)*dV - phiklong(kk)

!       write(16,'(2i5,6e15.5,6i5)') j,i,Screening%phik(j), &
!            phir(i),dV,phir(i)*dV,phiklong(kk), sqrt(kk),&
!            ia,ib,ic,ja,jb,jc
!       write(18,'(6f10.6)') K(j)%x, K(j)%y, K(j)%z, kk, phir(i)*dV, phiklong(kk)
!       write(*,*) Screening%phik(j)
    end do
!    write (*,*) "Debug informations for the screening module"
!    write (*,*) "  Screening%Gmax = ", Screening%Gmax
!    write (*,*) "  n_rGpv() = ", n_rGpv(1), n_rGpv(2), n_rGpv(3)
!    write (*,*) "  N = ", Na, Nb, Nc
!    write (*,*) "  Screening%Nk = ", Screening%Nk
    ! 遮蔽ポテンシャルの値のテスト表示
!    write(*,*) "# i kx ky kz kk PhiK(i)"
!    do j=1, Screening%Nk
!       if( j>32 ) then
!          write(*,*) "omitted"
!          exit
!       end if
!       write(13,100) j, K(j)%x, K(j)%y, K(j)%z, length2(K(j)), Screening%phik(j)
!    end do

    deallocate( K )
    deallocate( phir )

100 format(i4,5f11.6)

    return
  end subroutine

  ! 実空間での遮蔽ポテンシャルの長距離成分を返す関数
  real(kind=DP) function phirlong( r )
    real(kind=DP) r
    real(kind=DP) derf ! 誤差関数

    if( r == 0.d0 ) then ! 原点は特別に対応する
       phirlong = Screening%alpha*M_2_SQRTPI
    else
       phirlong = derf(Screening%alpha*r)/r
    end if
  end function

  ! 逆空間での遮蔽ポテンシャルの長距離成分を返す関数
  real(kind=DP) function phiklong( kk )
    real(kind=DP) kk

    if( kk == 0.d0 ) then ! 原点は特別に対応する
       phiklong = -M_PI/(Screening%alpha*Screening%alpha)
    else
       phiklong = 4.0d0*M_PI/kk*exp(-kk/(4.0d0*Screening%alpha*Screening%alpha))
    end if
  end function

end module
