!ソースコードを簡略化するためのマクロ群(fortranのソースコードでのみ使用)
!Copyright (C) 2014-2015 Ryosuke Tomita
        
#define nxyz_ 3
#define NameLen 20
#define NOP_MAX_CUB 48
#define NOP_MAX_HEX 24

#define CUBIC_C2 2:4,13:18
#define CUBIC_C3 5:12
#define CUBIC_C4 19:24
#define CUBIC_IE 25
#define CUBIC_IC2 26:28,37:42

#define HEX_C2 4,7:12
#define HEX_C3 3,5
#define HEX_C6 2,6
#define HEX_IE 13
#define HEX_IC2 16,19:24

#define READSELECT(nfin) \
  select case(readline(nfin));\
  case(-1); return;\
  case(0); exit;\
  end select
#define READTAG(tag) (index(linebuf, tag) > 0)
#define READLINEERR(nfin) (readline(nfin) < 1)
#define READERR(nfin, tag) (READLINEERR(nfin) .or. READTAG(tag))

#define COMPLEX_T complex(kind(0d0))
#define DOUBLE_T double precision
#define TAG_T character(len=2)

#define CMPLX_(a,b) dcmplx(a,b)
#define CEXP_(a) cdexp(a)
#define TEX_CAPTION_FOR_CHARTABLE(K) '\caption{Character table and relation between classes and operators for $'//K//'$}'
#define TEX_CAPTION_FOR_COMTABLE(K1,K2) '\caption{Compatibility table for $'//K1//'$ and $'//K2//'$}'
#define TEX_BEGINTABLE '\begin{table}[htbp]'
#define TEX_ADJUST '\centering'
#define TEX_ENDTABLE '\end{table}'
#define TEX_BEGINTABULAR(NC) '\begin{tabular}{*{',(NC),'}{|c}|}\hline'
#define TEX_ENDTABULAR '\end{tabular}'
#define TEX_COMPAT_NAME(K,IMR) (trim(K)//'_'//trim(IMR))
