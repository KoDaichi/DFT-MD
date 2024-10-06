.. _install_chapter:

PHASE/0のインストール
=====================

PHASE/0をインストールする方法を説明します。個別の環境についてのインストール方法について `このGitHub <https://github.com/Materials-Science-Software-Consortium/phase0_install>`_ のページにも情報があります。

動作環境
--------

PHASE/0は、PCから最先端のスーパーコンピュータの様々な計算機環境で動作します。

PHASE/0プログラムはFortran（Fortran90およびFORTRAN77）とCで記述されています。これらのコンパイラが使える計算機システムが必要です。並列計算をする場合にはMPIライブラリがインストールされている必要があります。

必要（利用可能）なソフトウェア、ライブラリ

-  Fortran90コンパイラ、Cコンパイラ（必須）
-  MPIライブラリ（並列計算に必須）
-  行列演算ライブラリLAPACK, BLAS（オプション）
-  FFTライブラリFFTW（オプション）
-  Perl（オプション）・・・・PHASEツールで必要
-  Gnuplot（オプション）・・・・PHASEツールで必要
-  Python（オプション）・・・・PHASEツール，ASEで必要

PHASEの動作確認を行っている計算機環境を以下に示します。

PHASE/0が動作する主な計算機環境

===================== ================ =======================
計算機環境            コンパイラー     利用可能ライブラリ
===================== ================ =======================
Linux                 GNU Compiler     LAPACK, BLAS, ScaLAPACK,
                      Intel Compiler   MKL, FFTW3
macOS (Apple Silicon) GNU Compiler     LAPACK, BLAS, ScaLAPACK, FFTW3
NEC SX Series         Fortran90/SX     Mathkeisan(LAPACK)
                                       ASL(FFT)
Fujitsu FX1000, FX700 Fujitsu Compiler
===================== ================ =======================

-  MPIライブラリは、OpenMPI, IntelMPI, SGI MPTなどに対応しています。
-  MKLのFFTW3ラッパーを介してMKLのFFTを利用することもできます。
-  ASLのFFTW3ラッパーを介してASLのFFTを利用することもできます。

2次元版および3次元版について
----------------------------

PHASE/0は，バンドおよび\ **k**\ 点の2軸で並列化を施した2次元版と，バンド，\ **k**\ 点に加え\ **G**\ 点の並列軸を加えた3次元版があります。低並列時では2次元版の方が高速に動作することが多いですが，高並列時は3次元版に軍配があがります。ただし，一部機能は2次元版でしか利用することはできません。

インストール方法
----------------

Linux環境を例にしてインストール方法を説明します。ここでは、Linux環境にIntel Fortran compilerがインストールされていることを想定します。
別のコンパイラを利用する場合はコンパイラの選択のときに、そのコンパイラを選択してください。また、MPIライブラリとして、OpenMPIを使用することを仮定しています。プログラムモデルの選択のときに「Serial」を選択すればMPIライブラリを使用しないでプログラムを作成することもできます。

まず、PHASE/0パッケージ |PHASE020XX.YY|.tar.gzをPHASE/0をインストールするディレクトリに展開してください。

.. parsed-literal::

 $ tar zxf |PHASE020XX.YY|.tar.gz

ディレクトリ |PHASE020XX.YY| に移り、インストーラーを実行してください(2次元版はinstall.sh, 3次元版はinstall_3d.sh；以下は2次元版の例)

.. parsed-literal::

 $ cd |PHASE020XX.YY|
 $ ./install.sh

  === PHASE installer ===
  Do you want to install PHASE? (yes/no) [yes]

インストールするかどうか聞いてきますので、何も入力せずにEnterキーを押してください。

.. code-block:: text

 Supported platforms

 0) GNU Linux (IA32)
 1) GNU Linux (EM64T/AMD64)
 2) NEC SX Series
 x) Exit

 Enter number of your platform. [0]

対応する環境の一覧が表示されますので、「GNU Linux
(EM64T/AMD64)」に対応する1を入力して、Enterキーを押してください。

.. code-block:: text

 Supported compilers
 0) GNU compiler collection (gfortran)
 1) Intel Fortran compiler
 x) Exit
 Enter number of a desired compiler. [0]

対応するコンパイラーの一覧が表示されますので、「Intel Fortran compiler」に対応する1を入力して、Enterキーを押してください。

.. code-block:: text

 Supported programming-models
 0) Serial
 1) MPI parallel
 x) Exit
 Enter number of a desired programming-model. [0]

対応するプログラムモデルの一覧が表示されますので、「MPI
parallel」に対応する1を入力して、Enterキーを押してください。

.. code-block:: text

 Supported MPI libraries
 0) Open MPI
 1) SGI MPT
 2) Intel(R) MPI
 x) Exit
 Enter number of a desired MPI library. [0]

対応するMPIライブラリの一覧が表示されますので、「Open MPI」に対応する0を入力して、Enterキーを押してください。

.. code-block:: text

 Supported BLAS/LAPACK
 0) Netlib BLAS/LAPACK
 1) Intel Math Kernel Library (MKL)
 x) Exit
 Enter number of a desired library. [0]

対応するBLAS/LAPACKライブラリの一覧が表示されますので、「Netlib BLAS/LAPACK」に対応する0を入力して、Enterキーを押してください。

.. code-block:: text

 Supported FFT libraries
 0) Built-in FFT subroutnes
 1) FFTW3 library
 x) Exit
 Enter number of a desired library. [0]

対応するFFTライブラリの一覧が表示されますので、「Built-in FFT subroutnes」に対応する0を入力して、Enterキーを押してください。

.. code-block:: text

 Do you want to edit the makefile that has been generated? (yes/no/exit) [no]

作成されたMakefileを編集するかどうか聞いてきます。Makefileの確認や編集を行う必要がなければ、何も入力せずにEnterキーを押してください。

.. code-block:: text

 Do you want to make PHASE now? (yes/no) [yes]

PHASEのコンパイルとインストールを開始するかどうか聞いてきます。何も入力せずにEnterキーを押して、PHASEのコンパイルとインストールを始めてください。

.. code-block:: text

 PHASE was successfully installed.
 Do you want to check the installed programs? (yes/no) [no]

PHASEが正常にインストールされたことを告げるメッセージの後、プログラムのテスト計算を実行するかどうか聞いてくるので、必要があればyesを入力し、Enterキーを押してください。テスト計算をしないならば、noを入力してEnterキーを押してください。テスト計算を実行して以下のような出力が得られれば問題ありません。

.. code-block:: text

 Do you want to check the installed programs? (yes/no) [no]
 yes
 Checking total-energy calculation.
  Total energy : -7.897015156331 Hartree/cell
  Reference    : -7.897015156332 Hartree/cell
 Checking band-energy calculation.
  Valence band maximum : 0.233846 Hartree
  Reference            : 0.233846 Hartree

MPIプログラムの実行に用いるmpirunやmpiexecなどのコマンドを用いて実行します。

$HOME/|PHASE020XX.YY|/binを環境変数PATHに追加しておくと、PHASE/0のプログラムのパスを指定せずに実行でき便利です。

Bourne shell(ボーンシェル)系であれば、$HOME/.bashrcなどにPATHを記述します。

.. parsed-literal::

  export PATH=$HOME/|PHASE020XX.YY|/bin:$PATH

C shell(シーシェル)系であれば、$HOME/.cshrcにPATHを記述します。

.. parsed-literal::

  setenv PATH $HOME/|PHASE020XX.YY|/bin:$PATH

MPIライブラリのbinディレクトリにも必ずパスを通すようにしてください。

Bourne shell(ボーンシェル)系であれば、$HOME/.bashrcなどにPATHを記述します。

.. code-block:: shell

  export PATH=$HOME/openmpi/bin:$PATH

C shell(シーシェル)系であれば、$HOME/.cshrcにPATHを記述します。

.. code-block:: shell

  setenv PATH $HOME/openmpi/bin:$PATH

以下のようにして、PHASEを実行します。

.. code-block:: shell

  $ mpirun -np 2 phase ne=1 nk=2

**gfortranを用いる場合の注意点**

Fortranコンパイラーとしてバージョン10以上のgfortranを用いる場合、上述のinstall.shスクリプトを用いた方法は失敗します。\ ``src_phase/Makefile`` を編集する必要があります。

.. code-block:: make

   F90 = mpif90 -m64

この行に ``-fallow-argument-mismatch`` という文字列を追加したのちに ``make`` コマンドを実行してください。


.. _section_install_mpifftw:

Distributed-memory FFTWとリンクする方法（バージョン2022.01以降）
------------------------------------------------------------------
3次元版PHASE/0を `Distributed-memory FFTWライブラリー <https://fftw.org/doc/Distributed_002dmemory-FFTW-with-MPI.html>`_ とリンクし、MPI並列FFTの処理をこのライブラリーに任せることができます。MPI並列FFTの処理は3次元版のPHASE/0に組み込まれており、3次元FFTの2軸を面とみなし分割する仕組みとなっています。これに対しDistributed-memory FFTWは1軸を分割します。したがってPHASE/0内蔵の並列FFTの方がスケーラビリティが高いのですが、系の形状やカットオフエネルギー、G点並列数によっては1軸分割の方が通信が少なくなる場合があり、そのような場合に利用するとより高速な計算を実現することができます。

リンクする方法は環境などによって異なります。ここでは `Intel MKLに付属するdistributed-memory FFTW <https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/appendix-c-fftw-interface-to-onemkl/fftw3-interface-to-onemkl/mpi-fftw3-wrappers.html>`_ とリンクする方法を紹介します。

まずはラッパーをコンパイルします。ここではホームディレクトリーの下の ``mkl`` というディレクトリーにインストールする前提のコマンドを紹介します。

.. code-block:: shell

 cp -r $MKLROOT/interfaces/fftw3x_cdft .
 cd fftw3x_cdft
 make libintel64 MKLROOT=$MKLROOT INSTALL_DIR=$HOME/mkl

Makefileを以下のように編集します。

- ``CPPFLAGS = -D_USE_DATE_AND_TIME_ ...`` に ``-DMPI_FFTW`` を追加
- ``--start-group ... ${MKLHOME}/libmkl_sequential.a ...`` に ``${HOME}/mkl/libfftw3x_cdft_lp64.a ${MKLHOME}/libmkl_cdft_core.a`` を追加

環境変数CPATHにfftwのインクルードディレクトリーを追加します。

.. code-block:: shell

 export CPATH=$MKLROOT/include/fftw:$CPATH

この状態で ``make clean;make`` とするとdistributed-memory FFTWを利用することのできるバイナリーを得ることができます。

デフォルトの状態では内蔵の並列FFTを用います。Distributed-memory FFTWを使う場合は以下のような設定を施します。

.. code-block:: text

   control{
     sw_mpi_fftw = on
   }

.. _section_install_eigenexa:

EigenExaとリンクする方法（バージョン2024.01以降）
------------------------------------------------------------------
3次元版PHASE/0を `EigenExa <https://www.r-ccs.riken.jp/labs/lpnctrt/projects/eigenexa/>`_ とリンクし、分散並列による行列対角化をこのライブラリーに任せることができます。PHASE/0の標準的な分散並列行列対角化ライブラリーはScaLAPACKですが、計算機アーキテクチャーや問題規模、並列数によってはEigenExaの方が高速に動作する場合があります。

リンクする方法は環境によって異なります。ここでは、MPIとしてインテルMPI, コンパイラーとしてインテル Fortranコンパイラー・クラシックを用いた例を紹介します。

まずはEigenExaをコンパイルします。上記のサイトから最新版をダウンロードし、 ``README.md`` ファイルの指示に従ってコンパイルします。

.. code-block:: shell

  tar -zxvf EigeExa-2.12.tar.gz
  cd EigenExa-2.12
  ./bootstrap
  ./configure
  make install

インテルMPIの環境が設定されている場合、自動的にそれが選択されます。

Makefileに以下のような編集を施す必要があります。 

- ``CPPFLAGS = -D_USE_DATE_AND_TIME_ ...`` に ``-DUSE_EIGENLIB -DEIGEN_EXA_2_9`` を追加します
- ``LIBS`` から始まる行にEigenExaのライブラリーを次のように指定します ``LIBS = -L$(HOME)/EigenExa-2.12/lib -lEigenExa ...``
- OpenMPを有効にするため、コンパイラーオプションに ``-qopenmp`` を加えます。 ``F90 = mpiifort -qopenmp`` ``LINK = mpiifort -qopenmp``
- ESMを有効にしている場合、 ``libesm.a`` の次の行における ``-D__MPI__`` という文字列を ``-D__OPENMP`` に差し替えます

環境変数 ``CPATH`` にEigenExaのインクルードディレクトリーを追加します。

.. code-block:: shell

 export CPATH=$HOME/EigenExa-2.12/include:$CPATH

この状態で ``make clean;make`` とするとEigenExaを利用することのできる実行可能ファイルを得ることができます。

デフォルトの状態ではScaLAPACKを用います。EigenExaを使う場合は以下のような設定を施します。

.. code-block:: text

  wavefunction_solver{
    submat{
      eigen_exa{
        sw_eigen_exa = on
      }
    }
  }

また、実行する前に環境変数 ``LD_LIBRARY_PATH`` にEigenExaのライブラリーファイルが配置されているディレクトリーを指定する必要があります。

.. code-block:: shell

 export LD_LIBRARY_PATH=$HOME/EigenExa-2.12/lib:$LD_LIBRARY_PATH

さらに、環境変数 ``OMP_NUM_THREADS`` によってMPI並列数や利用可能なコア数を勘案しながら適切なスレッド並列数を指定します。

.. code-block:: shell

 export OMP_NUM_THREADS=2

なお、EigenExaはIntel oneAPI 2024年版以降に同梱されているインテルMPIとコンパイラーではコンパイルできなくなっているようです。2023年版以下のバージョンを用いるようにしてください。また、EigenExaとk点並列を同時に利用することはできません。

