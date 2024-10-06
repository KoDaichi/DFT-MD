/***
!=======================================================================
!
!  PROGRAM  PHASE/0 2014.03 ($Rev: 409 $)
!
!
!  AUTHOR(S): K. Mae   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================
***/
/*****************************************************************************
Header file for input-file parse program input_parse.c

last modification: 2002.12.26 by K. Mae
******************************************************************************/
/*#include	<time.h>*/
#include	<stdio.h>
/*#include	<fcntl.h>*/
/*#include	<memory.h>*/
#include	<stdlib.h>
#ifndef __APPLE__
#include	<malloc.h>
#else
#include        <malloc/malloc.h>
#endif
#include	<ctype.h>
#include	<string.h>
#include        <math.h>

#ifdef ABSOFT
#define openInputFile			OPENINPUTFILE
#define closeInputFile			CLOSEINPUTFILE
#define selectTop				SELECTTOP
#define selectBlock				SELECTBLOCK
#define selectParentBlock		SELECTPARENTBLOCK
#define getIntValue				GETINTVALUE
#define getRealValue			GETREALVALUE
#define getIntVectorValue		GETINTVECTORVALUE
#define getRealVectorValue		GETREALVECTORVALUE
#define getStringValue			GETSTRINGVALUE
#define selectFirstTableLine	SELECTFIRSTTABLELINE
#define selectNextTableLine		SELECTNEXTTABLELINE
#define getNumBlockUnits		GETNUMBLOCKUNITS
#define getBlockUnit			GETBLOCKUNIT
#define frac2real				FRAC2REAL
#define frac2intint				FRAC2INTINT
#elif INTEL
#define openInputFile			openinputfile_
#define closeInputFile			closeinputfile_
#define selectTop				selecttop_
#define selectBlock				selectblock_
#define selectParentBlock		selectparentblock_
#define getIntValue				getintvalue_
#define getRealValue			getrealvalue_
#define getIntVectorValue		getintvectorvalue_
#define getRealVectorValue		getrealvectorvalue_
#define getStringValue			getstringvalue_
#define selectFirstTableLine	selectfirsttableline_
#define selectNextTableLine		selectnexttableline_
#define getNumBlockUnits		getnumblockunits_
#define getBlockUnit			getblockunit_
#define frac2real				frac2real_
#define frac2intint				frac2intint_
#elif LOWERCASE
#define openInputFile                   openinputfile
#define closeInputFile                  closeinputfile
#define selectTop                       selecttop
#define selectBlock                     selectblock
#define selectParentBlock               selectparentblock
#define getIntValue                     getintvalue
#define getRealValue                    getrealvalue
#define getIntVectorValue               getintvectorvalue
#define getRealVectorValue              getrealvectorvalue
#define getStringValue                  getstringvalue
#define selectFirstTableLine            selectfirsttableline
#define selectNextTableLine             selectnexttableline
#define getNumBlockUnits                getnumblockunits
#define getBlockUnit                    getblockunit
#define frac2real                       frac2real
#define frac2intint                     frac2intint
#endif

#define MAXTAGLEN	256
#define MAXVALLEN	256
#define MAXBLKDEPTH	16
#define MAXNUMTAGS	128
#define MAXTBLTAGS	64
#define MAXNUMUNITS	16
#define MAXUNITLEN	20
#define MAXNUMDIVCOL	16
#define MAXNUMDIVLINE	16

#define cCR			'\n'
#define cCREQUIV	','
#define cSPACE		' '
#define cTAB		0x09
#define cEQ			'='
#define cEOF		0x1A
#define cCOMMENT	'!'
#define cCOMMAND	'#'
#define cBLKB		'{'
#define cBLKE		'}'
#define cDEFAULT	'*'
#define cQUOTE		'"'

#define sCOMMENT1	"!"
#define lenCOMMENT1	1
#define sCOMMENT2	"//"
#define lenCOMMENT2	2
#define sCOMMAND1	"!#"
#define lenCOMMAND1	2
#define sCOMMAND2	"#"
#define lenCOMMAND2	1


#define TABTAGSTR		"#tag"			/**allow !#tag      **/
#define TABDEFAULTSTR	"#default"		/**allow !#default  **/
#define TABUNITSTR		"#units"		/**allow !#units  **/
#define TABGROUPSTR		"#group"		/**allow !#group  **/
#define BOOLSTR1		"yes"
#define BOOLSTR0		"no"
#define BOOLSTR1a		"on"
#define BOOLSTR0a		"off"

#define eEOB			-999
#define eNOTBLOCK		-998
#define eSYNTAX			-1
#define eVALTYPE		-2
#define eALLOC			-3
#define eNOTABLETAG		-4
#define eNOTAG			-5
#define eINVALIDVAL		-6
#define eINVALIDUNIT	-7
#define eTABLE			-8
#define eBIGTABLE		-9

#define rEOB		999
#define rNORMAL		0
#define rNOTAG		1
#define rNOVAL		2
#define rNOLINE		3
#define rNOUNIT		4
#define rTOPBLOCK	5

#define TAGTYPE_EQ		1
#define TAGTYPE_BLOCK	2
#define TAGTYPE_GROUP	3
#define TAGTYPE_UNIT	4
#define TAGTYPE_TAG		5
#define TAGTYPE_DEFAULT	6

#define TAGTOP		"topoffile"

#define TOLOWER	1
#define TOUPPER	2


/*********************************
 variables
************************************/
/***** whole file *******/
char *rbuf;
int	rfilesize;
int wholesize;

/****** table info *****/
typedef struct table_info_t {
	int  numlines;
	int	 curline;
	int	 toppnt;
	int	 readpnt;
	int	 ntags;
	char tags[MAXTBLTAGS][MAXTAGLEN+1];
	char defaults[MAXTBLTAGS][MAXVALLEN+1];
} TABLE_INFO;

/***** current block info *******/
int blk_depth;
typedef struct blkinfo_t {
	char tag[MAXTAGLEN+1];
	int	 begin;
	int	 end;
	int	 nchildtags;
	char childtags[MAXNUMTAGS][MAXTAGLEN+1];
	int  childtype[MAXNUMTAGS];
	int	 readpnt[MAXNUMTAGS];
	int	 endpnt[MAXNUMTAGS];
	int  numtable;
	TABLE_INFO	*tblinfo;
	int  bigtablefg;
	int	 nunits;
	char units[MAXNUMUNITS][MAXUNITLEN+1];
	char group[MAXVALLEN+1];
} BLKINFO;
BLKINFO *blk_info, defaultinfo, wkinfo;


/***** big table info *******/
typedef struct bigtable_info_t {
	int  numdivcol;
	int  numdivline;
	int  numtotallines;
	int	 curdivl;
	int	 curline;
	int	 curend;
	char groupid[MAXNUMDIVLINE][MAXVALLEN+1];
	int  sumnumlines[MAXNUMDIVLINE];
	int  numlines[MAXNUMDIVLINE][MAXNUMDIVCOL];
	int	 toppnt[MAXNUMDIVLINE][MAXNUMDIVCOL];
	int	 endpnt[MAXNUMDIVLINE][MAXNUMDIVCOL];
	int	 readpnt[MAXNUMDIVLINE][MAXNUMDIVCOL];
	int	 ntags[MAXNUMDIVLINE][MAXNUMDIVCOL];
	char tags[MAXNUMDIVLINE][MAXNUMDIVCOL][MAXTBLTAGS][MAXTAGLEN+1];
	char defaults[MAXNUMDIVLINE][MAXNUMDIVCOL][MAXTBLTAGS][MAXVALLEN+1];
	int	 nunits[MAXNUMDIVLINE][MAXNUMDIVCOL];
	char units[MAXNUMDIVLINE][MAXNUMDIVCOL][MAXNUMUNITS][MAXUNITLEN+1];
} BIGTABLE_INFO;
BIGTABLE_INFO *bigtable;



/**** line end character detected by get_nextword_inline() ****/
char endchar;


