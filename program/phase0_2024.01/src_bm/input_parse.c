/***
!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 (rev.375)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  FUNCTION: errout, strupper, strlower, is_intstr, is_realstr, 
!            gotonextline_byCR, gotonextline, is_comment, 
!            is_command, is_separator, is_tabtagstr, is_tabunitsstr,
!         is_tabdefaultstr, is_tabgroupstr, skip_space_inline, skip_space
!          find_block_range, get_nextword_inline, get_nextword, 
!           get_tagno, get_tabletagno, get_readpnt, set_childtaginfo, 
!             set_tableinfo, set_blockunits, set_blockgroup, seek_val
!            get_defaultval, is_bigtable, set_bigtable, openInputFile,
!           closeInputFile, selectTop, selectBlock, getNumBlockUnits,
!           getBlockUnit, selectParentBlock, getIntValue, getRealValue,
!         getIntVectorValue, getRealVectorValue, getStringValue, 
!         selectFirstTableLine, selectNextTableLine, frac2real, frac2intint,
!
!  AUTHOR(S): K. Mae            July/04/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
***/
/******************************************************
CHASE-3PT input-file parse tool functions

last modification: 2002.12.26 by K. Mae
*******************************************************/
#include "input_parse.h"

int selectTop();
int selectParentBlock();

/****************************************************
Tool functions for string
*****************************************************/
void errout( char *tag, char *val )
{
	printf("Input file read error: '%s' has an invalid value '%s'.\n",tag, val );
}

int strupper( char *str )
{
int i;

	for( i = 0; str[i] != 0; i++ ){
		str[i] = toupper( str[i] );
	}
	return( i );
}

int strlower( char *str )
{
int i;

	for( i = 0; str[i] != 0; i++ ){
		str[i] = tolower( str[i] );
	}
	return( i );
}

int is_intstr( char *str )
{
int	i, ret = 1;
char c;
	
	for( i = 0; str[i] != 0; i++ ){
		c = str[i];
		if( !(c >= '0' && c <= '9' || (i==0 && (c=='+' || c=='-'))) ){
			ret = 0;
			break;
		}
	}
	return( ret );
}

int is_realstr( char *str )
{
int	i, ret = 1;
char c;
	
	for( i = 0; str[i] != 0; i++ ){
		c = str[i];
		if( !((c >= '0' && c <= '9') || c=='.' || c=='+' || c=='-'
		      || c=='e' || c=='E'|| c=='d'|| c=='D') ){
			ret = 0;
			break;
		}
	    if( c=='d'|| c=='D' ){
			str[i] = 'e';
		}
	}
	return( ret );
}


/***********************************************************************
 internal functions for getting input words
************************************************************************/
int gotonextline_byCR( char *buf, int p, int imax )
{
int i;

	for( i = p; i < imax; i++ ){
		if( buf[i] == cCR || buf[i] == cEOF ){
			i++;
			break;
		}
	}
	return( i );
}

int gotonextline( char *buf, int p, int imax )
{
int i;

	for( i = p; i < imax; i++ ){
		if( buf[i] == cCR || buf[i] == cCREQUIV || buf[i] == cEOF ){
			i++;
			break;
		}
	}
	return( i );
}

int is_comment( char *buf )
{
	/*if( *buf == cCOMMENT && *(buf+1) != cCOMMAND ){*/
	if( (strncmp(buf,sCOMMENT1,lenCOMMENT1) == 0 && strncmp(buf,sCOMMAND1,lenCOMMAND1) != 0)
	 || strncmp(buf,sCOMMENT2,lenCOMMENT2) == 0 ){
		return( 1 );
	}
	return( 0 );
}

int is_command( char *buf )
{
	/*if( *buf == cCOMMENT && *(buf+1) == cCOMMAND ){*/
	if( strncmp(buf,sCOMMAND1,lenCOMMAND1) == 0 || strncmp(buf,sCOMMAND2,lenCOMMAND2) == 0 ){
		return( 1 );
	}
	return( 0 );
}

int is_separator( char *buf )
{
char c;

	c = buf[0];
	/*if( c == cSPACE || c == cTAB || c == cCR || c == cCREQUIV || c == cEOF || is_comment(buf) ){*/
	if( c <= cSPACE || c == cCREQUIV || is_comment(buf) ){
		return( 1 );
	}
	return( 0 );
}

int is_tabtagstr( char *buf )
{
int s;

	s = (*buf==cCOMMENT)? 1: 0;
	if( strcmp( buf+s, TABTAGSTR ) == 0 ){
		return( 1 );
	}
	return( 0 );
}

int is_tabunitsstr( char *buf )
{
int s;

	s = (*buf==cCOMMENT)? 1: 0;
	if( strcmp( buf+s, TABUNITSTR ) == 0 ){
		return( 1 );
	}
	return( 0 );
}

int is_tabdefaultstr( char *buf )
{
int s;

	s = (*buf==cCOMMENT)? 1: 0;
	if( strcmp( buf+s, TABDEFAULTSTR ) == 0 ){
		return( 1 );
	}
	return( 0 );
}

int is_tabgroupstr( char *buf )
{
int s;

	s = (*buf==cCOMMENT)? 1: 0;
	if( strcmp( buf+s, TABGROUPSTR ) == 0 ){
		return( 1 );
	}
	return( 0 );
}


int skip_space_inline( char *buf, int p, int imax )
{
int i;

	i = p;
	while( i < imax ){
		/*if( buf[i] != cSPACE && buf[i] != cTAB ){*/
		if( buf[i] > cSPACE || buf[i] == cCR ){
			break;
		}
		i++;
	}
	if( is_comment( &buf[i] ) ){
		while( i < imax ){
			if( buf[i] == cCR /*|| buf[i] == cCREQUIV*/ ){
				break;
			}
			i++;
		}
	}
	return( i );
}

int skip_space( char *buf, int p, int imax, int *numcr )
{
int i, ncr;
char c;

	ncr = *numcr;
	i = p;
	while( i < imax ){
		c = buf[i];
		if( is_comment( &buf[i] ) ){
			if( (i = gotonextline_byCR( buf, i, imax )) == imax ){
				*numcr = ncr;
				return( imax );
			}
			ncr++;
			continue;
		}
		/*else if( c != cSPACE && c != cTAB && c != cCR && c != cCREQUIV  ){*/
		else if( c > cSPACE && c != cCREQUIV  ){
			break;
		}
		if( c == cCR ){
			ncr++;
		}
		i++;
	}
	*numcr = ncr;
	return( i );
}


int find_block_range( int *b, int *e, char *buf, int p, int imax )
{
int i, kakko = 0;
int lineno = 0; /*dummy*/

	if( (i = skip_space( buf, p, imax, &lineno )) == imax ){
		return( eEOB );
	}
	if( buf[i] != cBLKB ){
		return( eNOTBLOCK );
	}
	*b = i+1;
	kakko++;
	i++;
	while( i < imax ){
		if( buf[i] == cBLKE ){
			kakko--;
			if( kakko == 0 ){
				*e = i;
				break;
			}
		}
		else if( buf[i] == cBLKB ){
			kakko++;
		}
		else if( is_comment( &buf[i] ) ){
			if( (i = gotonextline_byCR( buf, i, imax )) == imax ){
				return( eEOB );
			}
			continue;
		}
		i++;
	}
	if( i >= imax ){
		i = imax;
	}
	return( i );
}


int get_nextword_inline( char *retstr, char *buf, int *p, int imax )
{
int i, j, quotefg = 0;
char c;

	i = *p;
	j = 0;
	if( (i = skip_space_inline( buf, i, imax )) == imax ){
		*p = i;
		return( rEOB );
	}
	while( i < imax ){
		c = buf[i];
		if( !quotefg ){
			if( is_separator( &buf[i] ) ){
				break;
			}
			else if( c == cBLKB || c == cBLKE || c == cEQ ){
				if( j == 0 ){
					retstr[j] = c;
					i++;
					j++;
				}
				break;
			}
		}
		if( c == cQUOTE ){
			quotefg = !quotefg;
		}
		if( j >= MAXVALLEN ){
			retstr[MAXVALLEN] = 0;
			printf("Input file read error: Tag or Value is too long. [%s]\n",retstr);
			return( eINVALIDVAL );
		}
		retstr[j] = c;
		i++;
		j++;
	}
	endchar = c;
	retstr[j] = 0;
	/*strlower( retstr );*/
	*p = i;

	if( i >= imax ){
		return( rEOB );
	}

	return( rNORMAL );
}


int get_nextword( char *retstr, char *buf, int *p, int imax, int *numcr )
{
int i, j, quotefg = 0;
char c;

	i = *p;
	j = 0;
	if( (i = skip_space( buf, i, imax, numcr )) == imax ){
		*p = i;
		return( rEOB );
	}
	while( i < imax ){
		c = buf[i];
		/*printf("i=%d c='%c' rbuf='%c'\n", i, c, rbuf[i]);*/
		if( !quotefg ){
			if( is_separator( &buf[i] ) ){
				break;
			}
			else if( c == cBLKB || c == cBLKE || c == cEQ ){
				if( j == 0 ){
					retstr[j] = c;
					i++;
					j++;
				}
				break;
			}
		}
		if( c == cQUOTE ){
			quotefg = !quotefg;
		}
		if( j >= MAXVALLEN ){
			retstr[MAXVALLEN] = 0;
			printf("Input file read error: Tag or Value is too long. [%s]\n",retstr);
			return( eINVALIDVAL );
		}
		retstr[j] = c;
		i++;
		j++;
	}
	retstr[j] = 0;
	/*strlower( retstr );*/
	*p = i;

	if( i >= imax ){
		return( rEOB );
	}

	return( rNORMAL );
}

int get_tagno( struct blkinfo_t *binfo, char *tag, int tagtype )
{
int i, n;
char wktag[MAXTAGLEN+1];

	strcpy( wktag, tag );
	strlower( wktag );
	n = binfo->nchildtags;
	for( i = 0; i < n; i++ ){
		if( binfo->childtype[i] != tagtype ){
			continue;
		}
		if( tagtype == TAGTYPE_GROUP   || tagtype == TAGTYPE_UNIT
		 || tagtype == TAGTYPE_DEFAULT || tagtype == TAGTYPE_TAG
		 || strcmp( wktag, binfo->childtags[i] ) == 0 ){
			break;
		}
	}
	if( i >= n ){
		return( eNOTAG );
	}

	return( i );
}


int get_tabletagno( struct blkinfo_t *binfo, char *tag )
{
int i, n;

	n = binfo->tblinfo->ntags;
	for( i = 0; i < n; i++ ){
		if( strcmp( tag, binfo->tblinfo->tags[i] ) == 0 ){
			break;
		}
	}
	if( i >= n ){
		return( eNOTAG );
	}

	return( i );
}


int get_readpnt( struct blkinfo_t *binfo, char *tag, int tagtype )
{
int i, n;

	n = binfo->nchildtags;
	for( i = 0; i < n; i++ ){
		if( binfo->childtype[i] != tagtype ){
			continue;
		}
		if( tagtype == TAGTYPE_GROUP   || tagtype == TAGTYPE_UNIT
		 || tagtype == TAGTYPE_DEFAULT || tagtype == TAGTYPE_TAG
		 || strcmp( tag, binfo->childtags[i] ) == 0 ){
			break;
		}
	}
	if( i >= n ){
		return( eNOTAG );
	}

	return( binfo->readpnt[i]);
}

int set_childtaginfo( struct blkinfo_t *binfo )
{
int i, wlen, b, e, ret, imax, n, tablefg;
char wbuf[MAXTAGLEN+1];

	i = binfo->begin;
	imax = binfo->end;
	binfo->nchildtags =	n = 0;

	tablefg = 0;
	binfo->tblinfo = NULL;

	while( i < imax ){
		get_nextword_inline( wbuf, rbuf, &i, imax );
		wlen = strlen( wbuf );
		if( wlen == 0 ){
			i = gotonextline( rbuf, i, imax );
		}
		else{
			if( tablefg && (!is_command(wbuf)) ){
				/**** skip table values *****/
				i = gotonextline( rbuf, i, imax );
				continue;
			}
			strlower(wbuf);
			strcpy( binfo->childtags[n], wbuf );
			/*** equation tag is read  *****/
			ret = find_block_range( &b, &e, rbuf, i, imax );
			if( ret == eNOTBLOCK ){
				/*** table command ****/
				if( is_command(wbuf) ){
					if( is_tabtagstr( wbuf ) ){
						binfo->childtype[n] = TAGTYPE_TAG;
						tablefg = 1;
					}
					else if( is_tabunitsstr( wbuf ) ){
						binfo->childtype[n] = TAGTYPE_UNIT;
					}
					else if( is_tabgroupstr( wbuf ) ){
						binfo->childtype[n] = TAGTYPE_GROUP;
					}
					else{
						return( eSYNTAX );
					}
					binfo->readpnt[n] = i;
				}
				/*** equation *****/
				else{
					get_nextword_inline( wbuf, rbuf, &i, imax );
					if( wbuf[0] != cEQ ){
						return( eSYNTAX );
					}
					binfo->childtype[n] = TAGTYPE_EQ;
					binfo->readpnt[n] = i;
				}
			}
			/****** error *********/
			else if( ret < 0 ){
				printf("Input file read error: Error in finding tags for '%s'.\n", binfo->tag);
				return( ret );
			}
			/*** block tag is read  *****/
			else{
				if( is_command(wbuf) ){
					binfo->childtype[n] = TAGTYPE_DEFAULT;
				}
				else{
					binfo->childtype[n] = TAGTYPE_BLOCK;
				}
				binfo->readpnt[n] = b;
				binfo->endpnt[n] = e;
				i = e;
			}
			i = gotonextline( rbuf, i, imax );
			n++;
		}
	}

	binfo->nchildtags = n;

	if( tablefg ){
		if( (binfo->tblinfo = (TABLE_INFO *)malloc(sizeof(TABLE_INFO)) ) == NULL ){
			printf("input file: get_childtaginfo, tableinfo malloc error\n");
			return( eALLOC );
		}
	}

/***/
#ifdef DEBUG
printf("tag='%s'\n",binfo->tag);
printf("begin=%d\n",binfo->begin);
printf("end=%d\n",binfo->end);
printf("nchildtags=%d\n",binfo->nchildtags);
for( i = 0; i < binfo->nchildtags;i++ ){
printf("childtags='%s' [%d]\n",binfo->childtags[i],binfo->childtype[i]);
}
scanf("%d",&i);
#endif
/****/
	return( rNORMAL );
}


int set_tableinfo( struct blkinfo_t *binfo )
{
int i, j, wlen, b, e, ret, imax, n, rp, lastno;
char wbuf[MAXTAGLEN+1];
/*char boolbuf[MAXTAGLEN+1];*/

	imax = binfo->end;
	if( (rp = get_readpnt( binfo, TABTAGSTR, TAGTYPE_TAG )) < 0 ){
		printf("input file: no table tags\n");
		return( eNOTABLETAG );
	}

	/**** set table tags ****/
	n = 0;
	while(1){
		get_nextword_inline( wbuf, rbuf, &rp, imax );
		wlen = strlen( wbuf );
		if( wlen == 0 ){
			break;
		}
		strcpy( binfo->tblinfo->tags[n], wbuf );
		binfo->tblinfo->defaults[n][0] = 0;
		n++;
	}
	binfo->tblinfo->ntags = n;

	/**** get read point of the 1st line ****/
/*printf("In set_Tableinfo\n");*/
	lastno = binfo->nchildtags-1;
/*printf("lastno=%d\n",lastno);*/
	if( binfo->childtype[lastno] == TAGTYPE_DEFAULT ){
		rp = binfo->endpnt[lastno];
	}
	else{
		rp = binfo->readpnt[lastno];
	}
/*printf("rp=%d\n",rp);*/
	rp = gotonextline( rbuf, rp, binfo->end );
/*printf("rp=%d after gptpnextline\n",rp);*/
	binfo->tblinfo->readpnt = binfo->tblinfo->toppnt = rp;

	/*** get the number of lines ****/
	if( rp == binfo->end ){
		n = 0;
	}
	else{
		n = 1;
	}
	while( 1 ){
		if( (rp = gotonextline( rbuf, rp, binfo->end )) == binfo->end ){
			break;
		}
		n++;
	}
	binfo->tblinfo->numlines = n;
/*printf("nimlines=%d\n",n);*/

	if( (i = get_tagno( binfo, TABDEFAULTSTR, TAGTYPE_DEFAULT )) < 0 ){
		/*** no default input *****/
		return( rNORMAL );
	}

	/**** get default inputs ****/
	strcpy( defaultinfo.tag, TABDEFAULTSTR );
	defaultinfo.begin = binfo->readpnt[i];
	defaultinfo.end = binfo->endpnt[i];
	set_childtaginfo( &defaultinfo );

	n = binfo->tblinfo->ntags;
	for( i = 0; i < defaultinfo.nchildtags; i++ ){
		rp = defaultinfo.readpnt[i];
		get_nextword_inline( wbuf, rbuf, &rp, defaultinfo.end );
		for( j = 0; j < n; j++ ){
			if( strcmp( defaultinfo.childtags[i], binfo->tblinfo->tags[j]) == 0 ){
				break;
			}
		}
		strcpy( binfo->tblinfo->tags[j], defaultinfo.childtags[i] );
		strcpy( binfo->tblinfo->defaults[j], wbuf );
		if( j >= n ){
			n++;
		}
	}

	binfo->tblinfo->ntags = n;


	return( rNORMAL );
}


int set_blockunits( struct blkinfo_t *binfo )
{
int i, j, wlen, b, e, ret, imax, n, rp;
char wbuf[MAXTAGLEN+1];

	imax = binfo->end;
	if( (rp = get_readpnt( binfo, TABUNITSTR, TAGTYPE_UNIT )) < 0 ){
		binfo->nunits = 0;
		binfo->units[0][0] = 0;
		return( rNORMAL );
	}

	/**** read units ****/
	n = 0;
	while(1){
		get_nextword_inline( wbuf, rbuf, &rp, imax );
		wlen = strlen( wbuf );
		if( wlen == 0 ){
			break;
		}
		else if( wlen > MAXUNITLEN ){
			binfo->units[n][0] = 0;
			binfo->nunits = n;
			printf("Input file read error: '%s' has an invalid unit '%s'.\n", binfo->tag, wbuf);
			return( eINVALIDVAL );
		}
		strlower( wbuf );
		strcpy( binfo->units[n], wbuf );
		n++;
	}
	binfo->units[n][0] = 0;
	binfo->nunits = n;

	return( rNORMAL );
}

int set_blockgroup( struct blkinfo_t *binfo )
{
int i, j, wlen, b, e, ret, imax, n, rp;
char wbuf[MAXVALLEN+1];

	imax = binfo->end;

	binfo->group[0] = 0;
	if( (rp = get_readpnt( binfo, TABGROUPSTR, TAGTYPE_GROUP )) < 0 ){
		return( rNORMAL );
	}

	ret = get_nextword_inline( wbuf, rbuf, &rp, imax );
	if( ret < 0 ){
		printf("Input file read error: '%s' has an invalid group id.\n", binfo->tag );
		return( eINVALIDVAL );
	}
	else if( (wlen = strlen( wbuf )) == 0 ){
		printf("Input file read error: '%s' has no group id value.\n", binfo->tag );
		return( eINVALIDVAL );
	}
	strlower( wbuf );
	strcpy( binfo->group, wbuf );

	return( rNORMAL );
}

int seek_val( char *tag, int *colno )
{
int i, j, k, tno, rp, n, cdivl;
char wbuf[MAXVALLEN+1];

	/*** normal input ***/
	if( blk_info[blk_depth].tblinfo == NULL && blk_info[blk_depth].bigtablefg == 0){
		if( (tno = get_tagno( &blk_info[blk_depth], tag, TAGTYPE_EQ )) < 0 ){
			return( eNOTAG );
		}
		rp = blk_info[blk_depth].readpnt[tno];
		*colno = 0;
	}
	/*** table input ***/
	else{
		if( blk_info[blk_depth].bigtablefg ){
			n = 0;
			cdivl = bigtable->curdivl;
			for( j = 0; j < bigtable->numdivcol; j++ ){
				for( i = 0; i < bigtable->ntags[cdivl][j]; i++ ){
					if( strcmp( tag, bigtable->tags[cdivl][j][i] ) == 0 ){
						break;
					}
					n++;
				}
				if( i < bigtable->ntags[cdivl][j] ){
					break;
				}
			}
			if( j >= bigtable->numdivcol ){
				return( eNOTAG );
			}
			rp = bigtable->readpnt[cdivl][j];
/*printf("dummy reading\n");*/
			for( k = 0; k < i; k++ ){
				get_nextword_inline( wbuf, rbuf, &rp, bigtable->endpnt[cdivl][j] );
/*printf("%s ", wbuf);*/
			}
			*colno = n;
			blk_info[blk_depth].end = bigtable->endpnt[cdivl][j];
/*printf("\n");*/
/*printf("divc=%d tagno=%d colno=%d readpnt=%d rp=%d\n",j,i,n,bigtable->readpnt[cdivl][j],rp);*/
		}
		else{
			if( (tno = get_tabletagno( &blk_info[blk_depth], tag )) < 0 ){
				return( eNOTAG );
			}
			rp = blk_info[blk_depth].tblinfo->readpnt;
			/**** seek by dummy reading ***/
			for( i = 0; i < tno; i++ ){
				get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );
			}
			*colno = tno;
		}
	}
	
	return( rp );
}

int get_defaultval( char *val, int colno )
{
int j, n, cdivl;

	if( blk_info[blk_depth].tblinfo == NULL && blk_info[blk_depth].bigtablefg == 0 ){
		return( rNOVAL );
	}
	else{
		if( blk_info[blk_depth].bigtablefg ){
			n = colno;
			cdivl = bigtable->curdivl;
			for( j = 0; j < bigtable->numdivcol; j++ ){
				if( n >= bigtable->ntags[cdivl][j] ){
					n -= bigtable->ntags[cdivl][j];
				}
				else{
					strcpy( val, bigtable->defaults[cdivl][j][n] );
					break;
				}
			}
		}
		else{
			strcpy( val, blk_info[blk_depth].tblinfo->defaults[colno] );
		}
		if( val[0] == 0 ){
			return( rNOVAL );
		}
	}
	return( rNORMAL );
}

int is_bigtable( struct blkinfo_t *binfo, char *tag )
{
int i, n, cnt;
char wktag[MAXTAGLEN+1];

	strcpy( wktag, tag );
	strlower( wktag );
	n = binfo->nchildtags;
	cnt = 0;
	for( i = 0; i < n; i++ ){
		if( strcmp( wktag, binfo->childtags[i] ) == 0 ){
			cnt++;
		}
	}
	if( cnt > 1 ){
		return( 1 );
	}
	return( 0 );
}

int set_bigtable( struct blkinfo_t *binfo, char *tag )
{
int i, j, k, kk, n, nn, cnt, ret, ndivc, ndivl, divc, divl;
char wktag[MAXTAGLEN+1];

	strcpy( wktag, tag );
	strlower( wktag );
	n = binfo->nchildtags;
	cnt = 0;
	ndivc = ndivl = 0;
	divc = divl = -1;
	bigtable->numtotallines; 
	for( i = 0; i < n; i++ ){
		if( strcmp( wktag, binfo->childtags[i] ) == 0 ){
			strcpy( wkinfo.tag, wktag );
			wkinfo.begin = binfo->readpnt[i];
			wkinfo.end = binfo->endpnt[i];
			wkinfo.bigtablefg = 0;
			if( (ret=set_childtaginfo( &wkinfo )) < 0 ){
				printf("In set_bigtable: set_childtaginfo=%d for %d th '%s' block.\n",ret,cnt,wktag);
				return( eSYNTAX);
			}
			if( wkinfo.tblinfo != NULL ){
				if( (ret=set_tableinfo( &wkinfo )) < 0 ){
					printf("In set_bigtable: set_tableinfo=%d for %d th '%s' block.\n",ret,cnt,wktag);
					free(wkinfo.tblinfo);
					return( eTABLE );
				}
			}
			else{
				printf("In set_bigtable: %d th '%s' block has no '%s' indication.\n",cnt,wktag,TABTAGSTR);
				return( eBIGTABLE);
			}
			if( (nn = get_tagno( &wkinfo, TABGROUPSTR, TAGTYPE_GROUP )) < 0 ){
				printf("In set_bigtable: Divided table '%s' must have '%s' indication.\n",wktag,TABGROUPSTR);
				free(wkinfo.tblinfo);
				return( eBIGTABLE );
			}
			if( (ret=set_blockgroup( &wkinfo )) < 0 ){
				printf("In set_bigtable: set_blockgroup=%d for %d th '%s' block.\n",ret,cnt,wktag);
				free(wkinfo.tblinfo);
				return( eINVALIDVAL);
			}
			if( (ret=set_blockunits( &wkinfo )) < 0 ){
				printf("In set_bigtable: set_blockunits=%d for %d th '%s' block.\n",ret,cnt,wktag);
				free(wkinfo.tblinfo);
				return( eINVALIDVAL);
			}
			for( j = 0; j < ndivl; j++ ){
				if( strcmp( wkinfo.group, bigtable->groupid[j] ) == 0 ){
					break;
				}
			}
			if( j < ndivl-1 ){		/*error*/
				printf("In set_bigtable: The same group id must appear successively for '%s'.\n", wktag);
				free(wkinfo.tblinfo);
				return( eBIGTABLE);
			}
			else if( j == ndivl-1 ){ /*add div col*/
				divc++;
				if( divc >= ndivc ){
					if( divl == 0 ){
						ndivc++;
					}
					else{
						printf("In set_bigtable: '%s' tables have different column divisions.\n",wktag);
						free(wkinfo.tblinfo);
						return( eBIGTABLE);
					}
				}
			}
			else{					/*add div line*/
				strcpy( bigtable->groupid[ndivl], wkinfo.group );
				ndivl++;
				divl++;
				if( cnt == 0 ){
					ndivc++;
				}
				divc = 0;
				bigtable->numtotallines += wkinfo.tblinfo->numlines;
				bigtable->sumnumlines[divl] = bigtable->numtotallines;
			}
			bigtable->numlines[divl][divc] = wkinfo.tblinfo->numlines;
			bigtable->toppnt[divl][divc] = wkinfo.tblinfo->toppnt;
			bigtable->endpnt[divl][divc] = wkinfo.end;
			bigtable->readpnt[divl][divc] = wkinfo.tblinfo->readpnt;
			bigtable->ntags[divl][divc] = nn = wkinfo.tblinfo->ntags;
			for( j = 0; j < nn; j++ ){
            	strcpy( bigtable->tags[divl][divc][j], wkinfo.tblinfo->tags[j] );
				strcpy( bigtable->defaults[divl][divc][j], wkinfo.tblinfo->defaults[j] );
			}
			bigtable->nunits[divl][divc] = nn = wkinfo.nunits;
			for( j = 0; j < nn; j++ ){
            	strcpy( bigtable->units[divl][divc][j], wkinfo.units[j] );
			}
			free(wkinfo.tblinfo);
			cnt++;
		}
	}
	bigtable->numdivcol = ndivc;
	bigtable->numdivline = ndivl;
	bigtable->curline = 0;
	bigtable->curdivl = 0;

	for( i = 0; i < ndivl; i++ ){
		for( j = 1; j < ndivc; j++ ){
			if( bigtable->numlines[i][j] != bigtable->numlines[i][0] ){
				printf("In set_bigtable: Column-divided tables '%s' have the different number of lines.\n",wktag);
				return( eBIGTABLE);
			}
		}
	}
	for( j = 0; j < ndivc; j++ ){
		for( i = 1; i < ndivl; i++ ){
			if( bigtable->ntags[i][j] != bigtable->ntags[0][j] ){
				printf("In set_bigtable: line-divided tables '%s' have the different number of columns.\n",wktag);
				return( eBIGTABLE);
			}
			for( k = 0; k < bigtable->ntags[i][j]; k++ ){
				for( kk = 0; kk < bigtable->ntags[0][j]; kk++ ){
					if( strcmp( bigtable->tags[i][j][k], bigtable->tags[0][j][kk] ) == 0 ){
						break;
					}
				}
				if( kk >= bigtable->ntags[0][j] ){
					printf("In set_bigtable: line-divided tables '%s' have different column tags.\n",wktag);
					return( eBIGTABLE);
				}
			}
		}
	}

	return( rNORMAL );
}

/**************************************************************************/


/**************************************************************************
 low-level interface functions
**************************************************************************/

/***************************************************
 open and read input file
****************************************************/
int openInputFile( char *fname )
{
FILE *fpr;
int i, j, wlen, linefg, inspacefg, tablefg, deffg, ret, eqfg;
int kakko, linetopfg, lineno, numcr, defkakkofg;
char wbuf[MAXTAGLEN+1];
char curtag[MAXTAGLEN+1];
char *wkrbuf;

/**
printf("openInputFile: fname=[%s]\n",fname);
scanf("%d",&i);
**/
	/*** open input file ****/
	if( (fpr = fopen( fname, "r" )) == NULL ){
		printf("input file: open error [%s]\n", fname);
		return( -1 );
	}

	/*** get file size ****/
	fseek( fpr, 0, SEEK_END );
	rfilesize = ftell( fpr );

#ifdef DEBUG
	printf("input file: file size=[%d]\n", rfilesize);
#endif

	/*** alloc temp read buffer ****/
	if( (wkrbuf=(char *)malloc( rfilesize+1 )) == NULL ){
		printf("input file: malloc error, size=[%d]\n", rfilesize);
		fclose( fpr );
		return( -1 );
	}

	/*** read file to the buffer ****/
	fseek( fpr, 0, SEEK_SET );
	wholesize = fread( wkrbuf, 1, rfilesize, fpr );
	if( ferror(fpr) ){
		printf("input file: fread error: wholesize=%d\n", wholesize );
		free( wkrbuf );
		fclose( fpr );
		return( -1 );
	}

	wkrbuf[wholesize] = 0;

#ifdef DEBUG
	printf("input file: wholesize=%d\n", wholesize );
	scanf("%d",&i);
	printf("%s\n", wkrbuf );
	scanf("%d",&i);
#endif

	fclose( fpr );

	/*** alloc read buffer ****/
	if( (rbuf=(char *)malloc( rfilesize+1 )) == NULL ){
		printf("input file: malloc error, size=[%d]\n", rfilesize);
		free( wkrbuf );
		return( -1 );
	}

	/*** reform to compact format ***/
	linefg = inspacefg = tablefg = deffg = defkakkofg = 0;
	kakko = 0;
	eqfg = 0;
	linetopfg = 1;
	lineno =1;
	i = j = 0;
	while( i < wholesize ){
		if( linefg ){
			ret = get_nextword_inline( wbuf, wkrbuf, &i, wholesize );
		}
		else{
			ret = get_nextword( wbuf, wkrbuf, &i, wholesize, &lineno );
		}
		if( ret < 0 ){
			free( wkrbuf );
			free( rbuf );
			return( -1 );
		}
		else if( ret == rEOB ){
			break;
		}
		/*printf("i=%d whole=%d word='%s'\n",i, wholesize, wbuf );*/
		wlen = strlen( wbuf );
		if( linefg && wlen == 0 ){
			rbuf[j] = cCR;
			j++;
			i = gotonextline( wkrbuf, i, wholesize );
			/***/
			if( deffg && !defkakkofg && endchar == cCR ){
				rbuf[j] = cBLKE;
				rbuf[j+1] = cCR;
				rbuf[j+2] = 0;
				j += 2;
				deffg = 0;
			}
			/***/
			linefg = 0;
			inspacefg = 0;
			linetopfg = 1;
			eqfg = 0;
			if( endchar == cCR ){
				lineno++;
			}
		}
		else if( wbuf[0] == cBLKB || wbuf[0] == cBLKE ){
			if( wbuf[0] == cBLKB ){
				if( linetopfg ){
					printf("Input file (Line=%d): No tag is found before '%c'.\n", lineno,wbuf[0]);
					free( wkrbuf );
					free( rbuf );
					return( -1 );
				}
				else if( eqfg ){
					printf("Input file (Line=%d): '%c' is found after '%c'.\n", lineno,wbuf[0],cEQ);
					free( wkrbuf );
					free( rbuf );
					return( -1 );
				}
				else if( is_tabdefaultstr(curtag) ){
					defkakkofg = 1;
					kakko++;
					continue;
				}
				else if( is_command(curtag) ){
					printf("Input file (Line=%d): '%c' is found in the command '%s'.\n", lineno,wbuf[0],curtag);
					free( wkrbuf );
					free( rbuf );
					return( -1 );
				}
				kakko++;
			}
			else{
				kakko--;
				if( rbuf[j-1] != cCR ){
					rbuf[j] = cCR;
					j++;
				}
				linefg = 0;
			}
			rbuf[j] = wbuf[0];
			rbuf[j+1] = cCR;
			rbuf[j+2] = 0;
			j += 2;
			tablefg = 0;
			inspacefg = 0;
			linetopfg = 1;
			eqfg = 0;
			defkakkofg = deffg = 0;
		}
		else if( wbuf[0] == cEQ ){
			if( linetopfg ){
				printf("Input file (Line=%d): No tag is found before '%c'.\n", lineno, wbuf[0]);
				free( wkrbuf );
				free( rbuf );
				return( -1 );
			}
			rbuf[j] = wbuf[0];
			rbuf[j+1] = 0;
			j++;
			linefg = 1;
			linetopfg = 0;
			eqfg = 1;
		}
		else if( is_command( wbuf ) ){
			strcpy( &rbuf[j], wbuf );
			strcpy( curtag, wbuf );
			j += wlen;
			linefg = 1;
			inspacefg = 1;
			if( is_tabtagstr( wbuf ) ){
				tablefg = 1;
			}
			else if( is_tabdefaultstr( wbuf ) ){
				/***/
				rbuf[j] = cBLKB;
				rbuf[j+1] = cCR;
				rbuf[j+2] = 0;
				j += 2;
				deffg = 1;
				inspacefg = 0;
				/****/
			}
			linetopfg = 0;
		}
		else if( wlen == 0){
			continue;
		}
		else{
			if( (!tablefg) && (!is_command(curtag)) && (!linetopfg) && (!eqfg) ){
				printf("Input file (Line=%d): The tag '%s' is followed by the value '%s' without '%c'.\n", lineno,curtag,wbuf,cEQ);
				free( wkrbuf );
				free( rbuf );
				return( -1 );
			}
			if( inspacefg ){
				rbuf[j] = cSPACE;
				j++;
			}
			strcpy( &rbuf[j], wbuf );
			if( linetopfg ){
				strcpy( curtag, wbuf );
			}
			j += wlen;
			if( tablefg && (!is_command(wbuf)) ){
				linefg = 1;
				inspacefg = 1;
			}
			else if( linefg ){
				inspacefg = 1;
			}
			linetopfg = 0;
		}
	}
	wholesize = j;
	free( wkrbuf );

	if( kakko != 0 ){
		printf("input file: { } do not match.\n");
		free( rbuf );
		return( -1 );
	}
	/**** alloc read buffer ****/
	if( (rbuf=(char *)realloc( rbuf, wholesize+1 )) == NULL ){
		printf("input file: realloc error, size=[%d]\n", wholesize);
		/*free( rbuf );*/
		return( -1 );
	}
#ifdef DEBUG
	printf("input file: new wholesize=%d\n", wholesize);
	rbuf[wholesize]=0;
	printf("%s\n", rbuf );
	scanf("%d",&i);
#endif

	/*** alloc blockinfo ******/
	if( (blk_info=(BLKINFO *)malloc( sizeof(BLKINFO)*MAXBLKDEPTH )) == NULL ){
		printf("input file: blockinfo malloc error, size=[%d]\n", sizeof(BLKINFO)*MAXBLKDEPTH);
		free( rbuf );
		return( -1 );
	}

	/*** alloc bigtable ******/
	if( (bigtable=(BIGTABLE_INFO *)malloc( sizeof(BIGTABLE_INFO) )) == NULL ){
		printf("input file: bigtable malloc error, size=[%d]\n", sizeof(BIGTABLE_INFO));
		free( rbuf );
		free( blk_info );
		return( -1 );
	}

	blk_depth = 0;
	strcpy( blk_info[blk_depth].tag, TAGTOP );
	blk_info[blk_depth].begin = 0;
	blk_info[blk_depth].end = wholesize;
	blk_info[blk_depth].bigtablefg = 0;

	if( set_childtaginfo( &blk_info[blk_depth] ) < 0 ){
		printf("input file: invalid tag found in the 0th block depth.\n");
		free( rbuf );
		free( blk_info );
		free( bigtable );
		return( -1 );
	}
	if( set_blockunits( &blk_info[blk_depth] ) < 0 ){
		printf("input file: invalid unit found in the beginning.\n");
		free( rbuf );
		free( blk_info );
		free( bigtable );
		return( -1 );
	}

	if( blk_info[blk_depth].tblinfo != NULL ){
		set_tableinfo( &blk_info[blk_depth] );
	}

/*printf("end open\n");*/
	return( rNORMAL );
}

/***************************************************
 dealloc read buffer
****************************************************/
int closeInputFile()
{
	selectTop();
	free(rbuf);
	free(blk_info);
	free( bigtable );
	return( rNORMAL );
}

/***************************************************
****************************************************/
int selectTop()
{
	for( ; blk_depth > 0; blk_depth-- ){
		if( blk_info[blk_depth].tblinfo != NULL ){
			free( blk_info[blk_depth].tblinfo );
			blk_info[blk_depth].tblinfo = NULL;
		}
	}

	return( rNORMAL );
}
/***************************************************
****************************************************/
int selectBlock( char *blocktag )
{
int i, dep, ret;
char wktag[MAXTAGLEN+1];
/**/
#ifdef DEBUG
printf("selectBlock: blocktag=[%s]\n",blocktag);
scanf("%d",&i);
#endif
/**/
	strcpy( wktag, blocktag );
	strlower( wktag );

	dep = blk_depth;

	if( is_bigtable( &blk_info[dep], wktag ) ){
		if( set_bigtable( &blk_info[dep], wktag ) < 0 ){
			return( eBIGTABLE );
		}
		blk_depth++;
		strcpy( blk_info[blk_depth].tag, wktag );
		blk_info[blk_depth].bigtablefg = 1;
		return( rNORMAL );
	}
	else{
		if( (i = get_tagno( &blk_info[dep], wktag, TAGTYPE_BLOCK )) < 0 ){
			return( rNOTAG );
		}
	}

	blk_depth++;
	strcpy( blk_info[blk_depth].tag, wktag );
	blk_info[blk_depth].begin = blk_info[dep].readpnt[i];
	blk_info[blk_depth].end = blk_info[dep].endpnt[i];
	blk_info[blk_depth].bigtablefg = 0;

	if( (ret=set_childtaginfo( &blk_info[blk_depth] )) < 0 ){
		printf("In selectBlock: set_childtaginfo=%d",ret);
		blk_depth--;
		return( eSYNTAX);
	}
	if( set_blockunits( &blk_info[blk_depth] ) < 0 ){
		printf("In selectBlock: set_blockunits=%d",ret);
		selectParentBlock();
		return( eINVALIDVAL);
	}

/*printf("In selectBlock: tblinfo=%x\n",blk_info[blk_depth].tblinfo);*/
	if( blk_info[blk_depth].tblinfo != NULL ){
		set_tableinfo( &blk_info[blk_depth] );
	}
#ifdef DEBUG
printf("selectBlock: end\n");
scanf("%d",&i);
#endif
	
	return( rNORMAL );
}

/***************************************************
 ****************************************************/
int getNumBlockUnits( int *numunits )
{
int n, j, cdivl;

	if( blk_info[blk_depth].bigtablefg ){
		n = 0;
		cdivl = bigtable->curdivl; 
		for( j = 0; j < bigtable->numdivcol; j++ ){
			n += bigtable->nunits[cdivl][j];
		}
		*numunits = n;
	}
	else{
		*numunits = blk_info[blk_depth].nunits;
	}
	return( rNORMAL );
}

/***************************************************
 ****************************************************/
int getBlockUnit( char *unit, int *no, int unitlen )
{
int i, j, n, len, ret, cdivl;

	n = *no;

	if( blk_info[blk_depth].bigtablefg ){
		cdivl = bigtable->curdivl; 
		for( j = 0; j < bigtable->numdivcol; j++ ){
			if( n >= bigtable->nunits[cdivl][j] ){
				n -= bigtable->nunits[cdivl][j];
			}
			else{
				len = strlen(bigtable->units[cdivl][j][n]);
				strncpy( unit, bigtable->units[cdivl][j][n], len );
				ret = rNORMAL;	
				break;
			}
		}
		if( j >= bigtable->numdivcol ){
			len = 0;
			ret = rNOUNIT;
		}
	}
	else{
		if( n >= blk_info[blk_depth].nunits ){
			len = 0;
			ret = rNOUNIT;
		}
		else{
			len = strlen(blk_info[blk_depth].units[n]);
			strncpy( unit, blk_info[blk_depth].units[n], len );
			ret = rNORMAL;	
		}
	}

	for( i = len; i < MAXUNITLEN; i++ ){
		unit[i] = cSPACE;
	}
	unitlen = MAXUNITLEN;
	return(ret);
}

/***************************************************
 ****************************************************/
int selectParentBlock()
{
	if( blk_depth == 0 ){
		return( rTOPBLOCK );
	}

	if( blk_info[blk_depth].tblinfo != NULL ){
		free( blk_info[blk_depth].tblinfo );
		blk_info[blk_depth].tblinfo = NULL;
	}
	blk_depth --;
	return( rNORMAL );
}

/***************************************************
get integer value at tag
****************************************************/
int getIntValue( char *tag, int *ret )
{
int rp, colno;
char wbuf[MAXVALLEN+1];
char wk[MAXTAGLEN+1];
int dmy;
/**/
#ifdef DEBUG
printf("getIntValue: tag=[%s]\n",tag);
scanf("%d",&dmy);
#endif
/**/
	strcpy( wk, tag );
	strlower( wk );

	if( (rp=seek_val( wk, &colno )) < 0 ){
		*ret = 0;
		return( rNOTAG );
	}

	get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );

	if( wbuf[0] == 0 || wbuf[0] == cDEFAULT ){
		if( get_defaultval( wbuf, colno ) != 0 ){
			*ret = 0;
			return( rNOVAL );
		}
	}

	strcpy( wk, wbuf );
	strlower( wk );
	if( strcmp( wk, BOOLSTR1 ) == 0 || strcmp( wk, BOOLSTR1a ) == 0){
		strcpy( wbuf, "1" );
	}
	else if( strcmp( wk, BOOLSTR0 ) == 0 || strcmp( wk, BOOLSTR0a ) == 0 ){
		strcpy( wbuf, "0" );
	}

	if( is_intstr( wbuf ) == 0 ){
		*ret = 0;
		/*errout( tag, wbuf );*/
		return( eINVALIDVAL );
	}
	*ret = atoi( wbuf );

	return( rNORMAL );
}

/***************************************************
get real value at tag
****************************************************/
int getRealValue( char *tag, double *ret, char *unit, int taglen, int unitlen )
{
int rp, colno;
double wkret;
char wbuf[MAXVALLEN+1];
char wktag[MAXTAGLEN+1];
int dmy, i, len, iret;
/**/
#ifdef DEBUG
printf("getRealValue: tag=[%s]\n",tag);
scanf("%d",&dmy);
#endif
/**/
	strcpy( wktag, tag );
	strlower( wktag );

	if( (rp=seek_val( wktag, &colno )) < 0 ){
		*ret = 0.0;
		return( rNOTAG );
	}

	get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );

	if( wbuf[0] == 0 || wbuf[0] == cDEFAULT ){
		if( get_defaultval( wbuf, colno ) != 0 ){
			*ret = 0.0;
			return( rNOVAL );
		}
	}

	if( is_realstr( wbuf ) == 0 ){
		*ret = 0.0;
		errout( tag, wbuf );
		return( eINVALIDVAL );
	}
	wkret = atof( wbuf );

	iret = rNORMAL;
	if( blk_info[blk_depth].tblinfo == NULL && blk_info[blk_depth].bigtablefg == 0 ){
		/**** read unit ****/
		get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );
		if( wbuf[0] == 0 ){	/** default unit ***/
			len = 0;
		}
		else{				/** given unit ***/
			len = strlen(wbuf);
			if( len > MAXUNITLEN ){
				len = 0;
				iret = eINVALIDUNIT;
			}
			else{
				strncpy( unit, wbuf, len );
			}
		}
	}
	else{
		len = 0;
	}
	for( i = len; i < MAXUNITLEN; i++ ){
		unit[i] = cSPACE;
	}
	unitlen = MAXUNITLEN;

	*ret = wkret;

	return( iret );
}

/***************************************************
get int vector value at tag
****************************************************/
int getIntVectorValue( char *tag, int *ret )
{
int i, j, rp, colno;
int wkret;
char wbuf[MAXVALLEN+1];
char wktag[MAXTAGLEN+1];
int dmy;
/**/
#ifdef DEBUG
printf("getRealVectorValue: tag=[%s]\n",tag);
scanf("%d",&dmy);
#endif
/**/
	strcpy( wktag, tag );
	strlower( wktag );

	/*** vector reading is NOT allowed for table ***/
	if( blk_info[blk_depth].tblinfo != NULL || blk_info[blk_depth].bigtablefg ){ 
		return( eTABLE );
	}

	if( (rp=seek_val( wktag, &colno )) < 0 ){
		return( rNOTAG );
	}

	for( i = 0; i < 3; i++ ){
		get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );

		if( wbuf[0] == 0 || wbuf[0] == cDEFAULT ){
			/*if( blk_info[blk_depth].tblinfo == NULL ){*/
				for( j = i; j < 3; j++ ){
					ret[j] = 0;
				}
				return( rNOVAL );
			/*}*/
			/********
			else{
				strcpy( wbuf, blk_info[blk_depth].tblinfo->defaults[colno] );
				if( wbuf[0] == 0 ){
					for( j = i; j < 3; j++ ){
						ret[j] = 0;
					}
					return( rNOVAL );
				}
			}
			********/
		}

		if( is_intstr( wbuf ) == 0 ){
			for( j = i; j < 3; j++ ){
				ret[j] = 0;
			}
			errout( tag, wbuf );
			return( eINVALIDVAL );
		}
		wkret = atoi( wbuf );
		ret[i] = wkret;
	}

	return( rNORMAL );
}

/***************************************************
get real vector value at tag
****************************************************/
int getRealVectorValue( char *tag, double *ret, char *unit, int taglen, int unitlen )
{
int i, j, rp, colno;
double wkret;
char wbuf[MAXVALLEN+1];
char wktag[MAXTAGLEN+1];
int dmy, iret, len;
/**/
#ifdef DEBUG
printf("getRealVectorValue: tag=[%s]\n",tag);
scanf("%d",&dmy);
#endif
/**/
	strcpy( wktag, tag );
	strlower( wktag );

	/*** vector reading is NOT allowed for table ***/
	if( blk_info[blk_depth].tblinfo != NULL || blk_info[blk_depth].bigtablefg ){ 
		return( eTABLE );
	}

	if( (rp=seek_val( wktag, &colno )) < 0 ){
		return( rNOTAG );
	}

	for( i = 0; i < 3; i++ ){
		get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );

		if( wbuf[0] == 0 || wbuf[0] == cDEFAULT ){
			/*if( blk_info[blk_depth].tblinfo == NULL ){*/
				for( j = i; j < 3; j++ ){
					ret[j] = 0.0;
				}
				return( rNOVAL );
			/*}*/
			/********
			else{
				strcpy( wbuf, blk_info[blk_depth].tblinfo->defaults[colno] );
				if( wbuf[0] == 0 ){
					for( j = i; j < 3; j++ ){
						ret[j] = 0.0;
					}
					return( rNOVAL );
				}
			}
			*********/
		}

		if( is_realstr( wbuf ) == 0 ){
			for( j = i; j < 3; j++ ){
				ret[j] = 0.0;
			}
			errout( tag, wbuf );
			return( eINVALIDVAL );
		}
		wkret = atof( wbuf );
		ret[i] = wkret;
	}

	iret = rNORMAL;
	/**** read unit ****/
	get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );
	if( wbuf[0] == 0 ){	/** default unit ***/
		len = 0;
	}
	else{				/** given unit ***/
		len = strlen(wbuf);
		if( len > MAXUNITLEN ){
			len = 0;
			iret = eINVALIDUNIT;
		}
		else{
			strncpy( unit, wbuf, len );
		}
	}

	for( i = len; i < MAXUNITLEN; i++ ){
		unit[i] = cSPACE;
	}
	unitlen = MAXUNITLEN;

	return( iret );
}

/***************************************************
get string value at tag
****************************************************/
int getStringValue( char *tag, char *ret, int *convfg, int taglen, int retlen )
{
int rp, colno, wlen;
char wbuf[MAXVALLEN+1];
char wktag[MAXTAGLEN+1];
int dmy, i;
/**/
#ifdef DEBUG
printf("getStringValue: tag=[%s]\n",tag);
scanf("%d",&dmy);
#endif
/**/
	
	strcpy( wktag, tag );
	strlower( wktag );

	if( (rp=seek_val( wktag, &colno )) < 0 ){
		ret[0] = 0;
		return( rNOTAG );
	}

	get_nextword_inline( wbuf, rbuf, &rp, blk_info[blk_depth].end );

	if( wbuf[0] == 0 || wbuf[0] == cDEFAULT ){
		if( get_defaultval( wbuf, colno ) != 0 ){
			ret[0] = 0;
			return( rNOVAL );
		}
	}
/**
printf("getStringValue: wbuf=[%s]\n",wbuf);
scanf("%d",&dmy);
**/
	if( wbuf[0] == cQUOTE ){
		wlen = strlen( wbuf );
		for( i = 1; i < wlen; i++ ){
			if( wbuf[i] == cQUOTE ){
				wbuf[i-1] = 0;
				break;
			}
			wbuf[i-1] = wbuf[i];
		}
	}

	if( *convfg == TOLOWER ){
		strlower( wbuf );
	}
	else if( *convfg == TOLOWER ){
		strupper( wbuf );
	}
	retlen = strlen(wbuf);
	strncpy( ret, wbuf, retlen );
	/*strcpy( retstring, wbuf );*/
	
/***/
	
	for( i = retlen; i < MAXVALLEN; i++ ){
		ret[i] = cSPACE;
	}
	retlen = MAXVALLEN;
/**/
/**
printf("getStringValue: ret=[%s]\n",ret);
scanf("%d",&dmy);
**/
	return( rNORMAL );
}


/***************************************************
seek the first line of a table
****************************************************/
int selectFirstTableLine()
{
int i, j, tno, rp, lastno;

	/*** check if this block has a table ***/
	if( blk_info[blk_depth].tblinfo == NULL && blk_info[blk_depth].bigtablefg == 0 ){
		printf("selectFirstTableLine: The block '%s' has no table.\n", blk_info[blk_depth].tag );
		return( eALLOC );
	}

	if( blk_info[blk_depth].bigtablefg ){
		for( i = 0; i < bigtable->numdivline; i++ ){
			for( j = 0; j < bigtable->numdivcol; j++ ){
				bigtable->readpnt[i][j] = bigtable->toppnt[i][j];
			}
		}
		bigtable->curline = 0;
		bigtable->curdivl = 0;
	}
	else{
	/*** table input ***/
	/**********
	lastno = blk_info[blk_depth].nchildtags-1;
	if( blk_info[blk_depth].childtype[lastno] == TAGTYPE_DEFAULT ){
		rp = blk_info[blk_depth].endpnt[lastno];
	}
	else{
		rp = blk_info[blk_depth].readpnt[lastno];
	}

	if( (rp = gotonextline( rbuf, rp, blk_info[blk_depth].end )) == blk_info[blk_depth].end ){
		return( rNOLINE );
	}
	blk_info[blk_depth].tblinfo->readpnt = rp;
	***********/
/*
printf("In selectFirstTableLine readpnt=%d toppnt=%d\n",blk_info[blk_depth].tblinfo->readpnt,blk_info[blk_depth].tblinfo->toppnt);
*/
		blk_info[blk_depth].tblinfo->readpnt = blk_info[blk_depth].tblinfo->toppnt;
		blk_info[blk_depth].tblinfo->curline = 0;
	}

	return( rNORMAL );
}

/***************************************************
seek the next line in a table
****************************************************/
int selectNextTableLine()
{
int rp, curdivl, j;

	/*** check if this block has a table ***/
	if( blk_info[blk_depth].tblinfo == NULL && blk_info[blk_depth].bigtablefg == 0 ){
		printf("selectNextTableLine: The block '%s' has no table.\n", blk_info[blk_depth].tag );
		return( eALLOC );
	}

	if( blk_info[blk_depth].bigtablefg ){
		if( bigtable->curline == bigtable->numtotallines-1 ){
			return( rNOLINE );
		}
		bigtable->curline++;
		curdivl = bigtable->curdivl;
		if( bigtable->curline >= bigtable->sumnumlines[curdivl] ){
			curdivl++;
			for( j = 0; j < bigtable->numdivcol; j++ ){
				bigtable->readpnt[curdivl][j] = bigtable->toppnt[curdivl][j];
			}
			bigtable->curdivl = curdivl;
		}
		else{
			for( j = 0; j < bigtable->numdivcol; j++ ){
				rp = bigtable->readpnt[curdivl][j];
				rp = gotonextline( rbuf, rp, bigtable->endpnt[curdivl][j] );
				bigtable->readpnt[curdivl][j] = rp;
			}
		}
	}
	else{
		rp = blk_info[blk_depth].tblinfo->readpnt;
		if( (rp = gotonextline( rbuf, rp, blk_info[blk_depth].end )) == blk_info[blk_depth].end ){
			return( rNOLINE );
		}
		blk_info[blk_depth].tblinfo->readpnt = rp;
		blk_info[blk_depth].tblinfo->curline++;
	}

	return( rNORMAL );
}



/************************************************************************************/

int frac2real( char *fracstr, double *ret )
{
int i, j, devfg, len;
char wbuf[MAXVALLEN+1];
double a, b;

	len = strlen( fracstr );
	/*printf("frac2real len=%d fracstr=[%s]\n",len,fracstr);*/

	for( i = 0; fracstr[i] == cSPACE; i++ );	/* slip space*/

	if( fracstr[i] == '0' || i >= len ){
		*ret = 0.0;
		return( rNORMAL );
	}

	/*** bunsi ***/
	j = 0;
	for( ; ((fracstr[i] >= '0' && fracstr[i] <= '9') || fracstr[i] == '-' || fracstr[i] == '.') && i < len; i++ ){
		wbuf[j] = fracstr[i];
		j++;
	}
	wbuf[j] = 0;
	a = atof( wbuf );

	for( ; fracstr[i] == cSPACE; i++ );		/* skip space */

	if( i >= len ){
		*ret = a;
		return( rNORMAL );
	}

	if( fracstr[i] == '/' ){
		i++;
	}
	else{
		*ret = 0.0;
		return( eINVALIDVAL );
	}

	for( ; fracstr[i] == cSPACE; i++ );	/* slip space*/

	if( fracstr[i] == '0' || i >= len ){
		*ret = 0.0;
		return( eINVALIDVAL );
	}

	/*** bunbo ****/
	j = 0;
	for( ; ((fracstr[i] >= '0' && fracstr[i] <= '9') || fracstr[i] == '-' || fracstr[i] == '.') && i < len; i++ ){
		wbuf[j] = fracstr[i];
		j++;
	}
	wbuf[j] = 0;
	b = atof( wbuf );
	*ret = a/b;
	/*printf("frac2real a=%lf b=%lf ret=%lf fracstr=[%s]\n",a,b,a/b,fracstr);*/
	return( rNORMAL );
}


int frac2intint( char *fracstr, int *retu, int *retd )
{
int i, j, devfg, len;
char wbuf[MAXVALLEN+1];
int a, b;

	len = strlen( fracstr );
	/*printf("frac2real len=%d fracstr=[%s]\n",len,fracstr);*/

	for( i = 0; fracstr[i] == cSPACE; i++ );	/* slip space*/

	if( fracstr[i] == '0' || i >= len ){
		*retu = 0;
		*retd = 1;
		return( rNORMAL );
	}

	/*** bunsi ***/
	j = 0;
	for( ; ((fracstr[i] >= '0' && fracstr[i] <= '9') || fracstr[i] == '-' /*|| fracstr[i] == '.'*/) && i < len; i++ ){
		wbuf[j] = fracstr[i];
		j++;
	}
	wbuf[j] = 0;
	a = atoi( wbuf );

	for( ; fracstr[i] == cSPACE; i++ );		/* skip space */

	if( i >= len ){
		*retu = a;
		*retd = 1;
		return( rNORMAL );
	}

	if( fracstr[i] == '/' ){
		i++;
	}
	else{
		*retu = 0;
		*retd = 1;
		return( eINVALIDVAL );
	}

	for( ; fracstr[i] == cSPACE; i++ );	/* slip space*/

	if( fracstr[i] == '0' || i >= len ){
		*retu = 0;
		*retd = 1;
		return( eINVALIDVAL );
	}

	/*** bunbo ****/
	j = 0;
	for( ; ((fracstr[i] >= '0' && fracstr[i] <= '9') || fracstr[i] == '-' /*|| fracstr[i] == '.'*/) && i < len; i++ ){
		wbuf[j] = fracstr[i];
		j++;
	}
	wbuf[j] = 0;
	b = atoi( wbuf );
	*retu = a;
	*retd = b;
	/*printf("frac2real a=%lf b=%lf ret=%lf fracstr=[%s]\n",a,b,a/b,fracstr);*/
	return( rNORMAL );
}


/************************************************************************************/

