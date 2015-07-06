//do we need to use the scope operator when we're inside the class, i.e. nbins is being used unqualified, and qualified within classes..
//#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include "ExprPredictor.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>


//#include "ExprPredictor.h"
#include <sstream>
//#include "SeqAnnotator.h"
 
bool isNt( int a )
{
    if ( a < 0 || a > 3 ) return false;
    else return true;	
}

int complement( int a )
{
    assert( a >= 0 && a < ALPHABET_SIZE );
            
    if ( a == 0 ) return 3;
    if ( a == 1 ) return 2;
    if ( a == 2 ) return 1;
    if ( a == 3 ) return 0;	
    if ( a == MISSING ) return MISSING;
    if ( a == GAP ) return GAP;	
}

int symbolToInt( char c )
{
    char upper = tolower( c );
    for ( int i = 0; i < ALPHABET_SIZE; i++ ) {
        if ( ALPHABET[ i ] == upper ) return i;	
    }
    
    return -1;
}

char strand2char( bool strand )
{
    if ( strand ) return '+';
    else return '-';	
}

bool char2strand( char c )
{
    assert( c == '+' || c == '-' );
    
    if ( c == '+' ) return true;
    else return false;
}

vector< double > createNtDistr( double gcContent )
{
    assert( gcContent >= 0 && gcContent <= 1.0 );
    
    vector< double > freqs( 4 );
    freqs[0] = ( 1.0 - gcContent ) / 2.0;
    freqs[1] = gcContent / 2.0;
    freqs[2] = freqs[1];
    freqs[3] = freqs[0];

    return freqs;
}

Sequence::Sequence( const string& str )
{
    for ( int i = 0; i < str.size(); i++ ) {
        int nt = symbolToInt( str[ i ] );	// could be a NNN or gap
        if ( nt >= 0 && nt < ALPHABET_SIZE ) {
            nts.push_back( nt );
        } else {
            cerr << "Illegal symbol: " << nt << " in " << str << endl;
            exit( 0 );	
        }       
    }
}

Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
    assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );	

    for ( int i = 0; i < length; i++ ) {
        if ( strand ) {	nts.push_back( other[ start + i ] ); }
        else { nts.push_back( complement( other[ start + length - 1 - i ] ) ); }
    }	
}

int Sequence::push_back( int nt )
{	
    assert( nt >= 0 && nt < ALPHABET_SIZE );
    nts.push_back( nt );
    
    return 0;
}

int Sequence::push_back( const Sequence& elem )
{
    for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ] );	
    return 0;
}

Sequence Sequence::compRevCompl() const
{
    return Sequence( *this, 0, size(), false );	
}

void Sequence::getNtCounts( vector< int >& counts ) const
{
    counts.clear();
    for ( int i = 0; i < NBASES; i++ ) {
        counts.push_back( 0 );
    }
    
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;	
    }
}

bool Sequence::containsMissing() const
{
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] == MISSING ) return true;	
    }	
    
    return false;
}

int Sequence::load( const string& file, string& name, int format )
{
    vector< Sequence > seqs;
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( seqs[ 0 ] );
    name = names[ 0 ];
    return rval;
}

int Sequence::load( const string& file, int format )
{
    string name;
    int rval = load( file, name, format );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Sequence& seq )
{
    // output the nts
    for ( int i = 0; i < seq.size(); i++ ) {
        os << ALPHABET[ seq[ i ] ];	
    }	
                    
    return os;
}

int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }
    seqs.clear();
    names.clear();
     
    // 	open the file
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

    string line;
    Sequence seq;
    
    // read sequences: FASTA format
    if ( format == FASTA ) {
        while ( getline( fin, line ) ) {
            // add the sequence and start a new sequence if the line starts with >
            //cout << line << endl;
            if ( line[ 0 ] == '>' ) { 	
                if ( seq.size() ) {
                    seqs.push_back( seq );
                    seq.clear();	
                }
                        
                stringstream ss( line.substr( 1 ) );
                string name; 
                ss >> name;
                names.push_back( name );
            } else { 
                // check if the line contains content
                int start = line.find_first_not_of( " \t\r" );
                int last = line.find_last_not_of( " \t\r" );
                if ( start == string::npos || last == string::npos ) continue;
                        
                // append the sequence	
                for ( int i = start; i <= last; i++ ) {
                    int nt = symbolToInt( line[ i ] );	// could be a NNN or gap
                    if ( nt >= 0 && nt < ALPHABET_SIZE ) {
                        seq.push_back( nt );
                    } else {
                        //cerr << "Illegal symbol: " << nt << " in " << file << endl;
                        return RET_ERROR;	
                    } 
                }
            }			
        }
            
        // add the last sequence
        if( seq.size() ) seqs.push_back( seq );
                        
        return 0;
    }	
}

int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );	
    return rval;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
    assert( seqs.size() == names.size() );
    
    // check if the format character is legal
    if ( format != FASTA ) { return RET_ERROR; }
            
    ofstream fout( file.c_str() );
    
    if ( format == FASTA ) {
        for ( int i = 0; i < seqs.size(); i++ ) {
            fout << ">" << names[ i ] << endl;
            fout << seqs[ i ] << endl;
        }
    }
    
    return 0;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
    // default name: integer starting from 1
    vector< string > names;
    for ( int i = 0; i < seqs.size(); i++ ) {
        char buffer[ 10 ];
        sprintf( buffer, "%i", i );
        names.push_back( string( buffer ) );	
    }	
    
    // print
    return writeSequences( file, seqs, names, format );
}

Matrix compWtmx( const Matrix& countMatrix, double pseudoCount )
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    
    int l = countMatrix.nRows();		// l: the length of motif		
    Matrix pwm( l, 4 );

//     // the sum of each position/column should be a const. (number of sequences)
//     double n = 0;		// number of sites used in the count matrix
//     for ( int j = 0; j < 4; j++ ) {
//         n += countMatrix( 0, j );
//     }
//     for ( int i = 1; i < l; i++ ) {
//         double count = 0;
//         for ( int j = 0; j < 4; j++ ) {
//             count += countMatrix( i, j );
//         }
//         if ( count != n ) { cout << "count matrix incorrect" << endl; exit( 1 ); }
//     }
    
    // the multinomial distribution at each column
    for ( int i = 0; i < l; i++ ) {
        double n = 0;       // total counts at this position 
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }	
    }

    return pwm;		
}

Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 )
{
    assert( background.size() == 4 );	
    
    init();
}

Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 )
{
    assert( background.size() == 4 );
    
    pwm = compWtmx( countMatrix, pseudoCount );
    init();
}

double Motif::LLR( const Sequence& elem ) const
{
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;
    
    double result = 0;
    for ( int i = 0; i < l; i++ ) {
        result += LLRMat( i, elem[ i ] ); // LLR matrix, M(i,b) = log( f_i(b) / p(b) ), where f_i(b) is the frequency of b at position i of PWM, and p(b) the frequency of b at the background	
    }
    
    return result;
}

double Motif::energy( const Sequence& elem ) const
{
	//return ( LLR( elem ) );
    return ( -LLR( elem ) +   maxLLR );	  /// this is delta G,  Gmaxsite - Gsite, expratio = exp( - dG ) = exp ( LLR- LLRmax), since LLRmax =10, and LLR < 10 it follows
// expratio = exp ( -x ) where x is positive or zero,  when x =0 we have q= maxBindwt*[c]*1.  if [c] =1; q=maxBindwt,  1/(1+100)
}

void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
    assert( rng != NULL );
    
    int l = pwm.nRows();
    Sequence sampleElem;
    for ( int i = 0; i < l; i++ ) {
        // nt. distribution at position i
        vector< double > distr = pwm.getRow( i );
        
        // sample nt. from this distribution	
        int nt = sampleMul( rng, distr );
        sampleElem.push_back( nt );
    }		
    
    if ( strand == 0 ) elem = sampleElem.compRevCompl();
    else elem = sampleElem;
}

int Motif::load( const string& file, const vector< double >& background, string& name )
{
    vector< Motif > motifs;
    vector< string > names;
    int rval = readMotifs( file, background, motifs, names );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( motifs[ 0 ] );
    name = names[ 0 ];
    return rval;				
}

int Motif::load( const string& file, const vector< double >& background )
{
    string name;
    int rval = load( file, background, name );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Motif& motif )
{
    os << motif.pwm;
    
    return os;
}

void Motif::init()
{
    int l = pwm.nRows();
    
    // compute the LLR matrix
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat( i, j ) = log( pwm( i, j ) / background[ j ] );
        }
    }
    
    // the strongest site
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite.push_back( b_max );	
    }
    
    // compute the LLR of the strongest site
    maxLLR = 0;
    for ( int i = 0; i < l; i++ ) {
        maxLLR += LLRMat( i, maxSite[ i ] );	
    }
}

Matrix countmatrixS( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       // return RET_ERROR;
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	//////////////see if you can apply maximum to seq, then that will allow for the largest sized matrix to be constructed.
//	int size = seq[1].size();

int length = tempmax ;  //maxstring(file.c_str());
//cout << " size " << size << endl;
//cout << " length  in countmatrix " << length << endl;



//ofstream otes("rcRead.txt");
	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
			//otes << ">" << seq_name[j] << endl << readseq <<endl;
			for( int i = 0;  i< readseq.size(); i++)  // this for loop is miss-counting by one? maybe overloaded size method?
			{
				if(readseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				//if ( readseq.nts.empty() ) continue;
			}
	
	


	} // for j
cout << " m "  << endl << m << endl;
	Matrix countMatrix = m.transpose();
	//counts = countMatrix;
cout << "countmat" << countMatrix << endl;
	 double pseudoCount =.1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		// l: the length of motif	
	//otes.close();
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
//cout << " length mal " << endl;
	   return pwm;
}


int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
    // 	open the file
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }
    motifs.clear(); 
    names.clear();

    string line;
    
    // read the motifs
    do {
        getline( fin, line );
        
        if ( line[ 0 ] != '>' ) continue;
        
        // read the names, length and pseudocount
        int MAX_SIZE = 100;
        char lineStr[ MAX_SIZE ];
        strcpy( lineStr, ( line.substr( 1 ) ).c_str() );
        char *name, *lengthStr, *pseudoCountStr;
        name = strtok( lineStr, " \t" );
        lengthStr = strtok( NULL, " \t" );
        pseudoCountStr = strtok( NULL, " \t" );
        int length;
        double pseudoCount;
        if ( lengthStr ) length = atoi( lengthStr );
        else { return RET_ERROR; }
        if ( pseudoCountStr ) pseudoCount = atof( pseudoCountStr );
        else pseudoCount = PSEUDO_COUNT;
        
        // read the count matrix
        Matrix countMat( length, NBASES );
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < NBASES; ++j ) {
                fin >> countMat( i, j );
            }	
        }
        
        // create the motif
        names.push_back( string( name ) );
        motifs.push_back( Motif( countMat, pseudoCount, background ) );	
    } while ( !fin.eof() );
                                    
    return 0;
}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
    vector< string > names;
    return readMotifs( file, background, motifs, names );	
}

ostream& operator<<( ostream& os, const Site& site )
{
    char strandChar = site.strand ? '+' : '-';
    os <<"\n"<< site.start + 1 << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;
    
    return os;
}

bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
	//if( a.start == b.start ) return false;
    
    return true;	
}
bool siteOverlap2( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
	if( a.start == b.start ) return false;
    
    return true;	
}
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy )
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    sites.clear();
    names.clear();

    SiteVec currVec;
    int nrecords = 0;       // number of ">" read so far
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;

        if ( line.substr( 0, 1 ) == ">" ) {
            stringstream ss( line.substr( 1 ) );
            string name; 
            ss >> name;
            names.push_back( name );
            nrecords++;
            if ( nrecords > 1 ) {
                sites.push_back( currVec );
                currVec.clear();
            }
        } else {
            int start;
            char strandChar;
            string factor;
            double energy = 0;
            stringstream ss( line );
            ss >> start >> strandChar >> factor;
            if ( readEnergy ) ss >> energy; 
            bool strand = strandChar == '+' ? 1 : 0;
            map<string, int>::const_iterator iter = factorIdxMap.find( factor );
            currVec.push_back( Site( start - 1, strand, iter->second , energy ) );
        }
    }

    sites.push_back( currVec );

    return 0;
}

int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy )
{
    vector< string > names;
    return readSites( file, factorIdxMap, sites, names, readEnergy );
}


int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
// 	cout << "start annotation:" << endl;
    sites.clear();
    
    // scan the sequence for the sites of all motifs
    for ( int i = 0; i < seq.size(); i++ ) {
        // test for each motifp.predictOcc( tsites, 5, cons, ff);
        for ( int k = 0; k < motifs.size(); k++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;
            
            // positive strand
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
//          cout << elem << Site( i, 1, k, exp( energy ) ) << endl;			
            if ( energy <= energyThrs[ k ] ) {  //1115 changed from <= to >=
                sites.push_back( Site( i, 1, k, energy ) );
            }	
            
            // negative strand
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
// 	    cout << rcElem << Site( i, 0, k, exp( energy ) ) << endl;			
            if ( energy <= energyThrs[ k ] ) {//1115 changed from <= to >= , 1129 changed back
                sites.push_back( Site( i, 0, k, energy ) );
            }				
        }	
    }
    
// 	cout << "end annotation" << endl;
    return sites.size();
}
bool siteSortPredicate(const Site& d1, const Site& d2)
{
  return d1.start < d2.start;
}
//////////////////////////////////// 8/12/11
// calculate Z in the context of dorsal, if twist coop site has higher Z than strong non-coop site Z, choose coop site.
/*
sitesp[0].push_back( siteMax( sitest[0] ));    
sitestoverlap( sitest[0], sitesp[0] ); 
vector< double  > sitest_tw_Z( 0 ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )

//for ( int k = 0; k < motifs.size(); k++ ) {
tsites.push_back( sitesp[0][0] );
SiteVec sitesorderedtemp = tsites;
	for( int i = 0; i < sitest[1].size(); i++){
		//sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
//cout << "sitest_tw_Z.size  " << sitest_tw_Z.size() << endl;
//cout << "maxZ " << maxZ( sitest_tw_Z ) << endl;
//cout << " sitest_sitest_allk_alli_Z[k] " << sitest_allk_alli_Z[k]<< endl;
int ind = maxZ( sitest_tw_Z );	
if ( cons[1] != 0 ) {
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
}
sitesp[2].push_back( siteMax( sitest[2] ));    
sitestoverlap( sitest[2], sitesp[2] ); 		        
//cout << " size of  sitest  " <<"\t" << sitest.size() << endl;
///////////////////////////////////////////////////////////////////
// initial maximum site has been created along with a temp V that does not contain overlapping sites with maximum site
//////////////////////////////////////////////////////////////////

// repeat siteMax up to 6 times for each motif. using both concentration and affinity and cooperation
int N = 1; // this is the number of sites
double previous;
double current;

do {
vector< double > ff(3); // focc passed by ref to predictocc.
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( previous > .5 )
{ break ; }
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
cout << " sitest_sitest_allk_alli_Z[k] " << sitest_allk_alli_Z[k]<< endl;
int ind = maxZ( sitest_allk_alli_Z[k] );	
cout << " ind " << ind << endl;
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  // for k
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( current > .5 )
{ break ; }
  // exit do loop once dl occ is almost 1.

//std::sort(sites.begin(), sites.end(), siteSortPredicate); 
} while( N < 6 );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
// sitesorderedtemp has sites ordered different than sitest, however the index of each Z,( in the vector sitest_allk..) is the same as the site index in sitest.
///////////////////////////////////////////////////////////////////


distribution constructor (sample space constructor).  We assume each nucleotide is a potential binding site.  If biological purpose of the nucleotide is to bind factor x, then we can see similarities with the pauli principle  which leads to the fermi dirac distribution. the nucleotide is bound or unbound and there can not exist more than one x protein bound.  The problem is different from the fermi distribution in the sense that there is more than one particle species.  The constraint we have is that each nucleotide is bound by factor x or y, or unbound.  To simply say that one uses two fd dists, (with four states xb xu, yb, yu) would indicate that yu is xb or xu.  However this notation is wrong, we want sx sy su,  (that is we are in the perspective of the sequence (not the perspective of the proteins))..  That is we say the the sequence is bound.. not the protein is bound..  sb not sx or sy..  and su not xu or yu..

For each nucleotide in the genome we do this.  Hence the sites are competing for occupancy.  To determine the occupancy of a given site in the single species picture (labeled by its chromsome coordiante) we simply create a big matrix (sample space constructor) of energy values.  then we simply apply the boltzman distribution to each state (if they were only one particle and there were a heat bath), however if there are more than one particle then we must place the constraint that only one particle can occupy one state at a time (pauli principle).  If the particles do not interact than all we have to do again is to create a list of their genetic energies and then fill the sites with lowest energy first (irrespective of their species identity).  Again applying the pauli principle.  ( So if particle x lowest energy site is agctttt and there are 1000 of those sites and particle x has 10 fold affinity for its lowest site than does particle y (which also has agctttt as its lowest energy site) , then if there  1000 particle x and 1000 particle y, particle x will occupy all of the agctttt sites, and particle y will have to occupy 1000 of its next lowest energy sites (which may be *gctttt or g*tttt or...gcttt* ) 

One of the aims of the modeling is to be able to say that particle x has 10 fold affinity for some site over particle y.  (or at least say xs > ys ).  

Once we have delabeled the particle species, our constraint becomes simply N the total particle number.  Then we can apply the fd stats, which yields all the occupied sites.  However if we wish to say which particle species is occuping a particular occupied site we have to put the labels back on.  The question then can be whether it is possible to simply keep the sites labeled.  (this is like following the trajectories of a particles and saying that each individual particle is different (it has a unique label).  But it's different in the sense that the different species have different affinities, whereas distinguishable particle problems usually assume each particle has the same energy for a given state (site).

One can calculate the chem p by saying N = sum(distribution)

But this chem p is not the same chem p, as the chem p from reactions.  There dG is per reaction, or per mole reaction (i.e. Na times per reaction), and dN  is the change is particle number, which if all the particles have the same label then dN is 1, and one can think of adding or subtracting one particle to the system. (the method is irrelevant (i.e. the path) since G is a state function) .  However if we put the labels back on the prob is not that simple sing dG =uxdgx + uydgy, but ux and uy and dgx and dgy are not state functions.  and one can't say if they remove one x, then dG = uxdgx.  (?? )

well what one can do is say that if there are 10 sites 5x and 4y all bound somewhere, then to determine which sites are occupied by who, one can simply go through the 3 to the 10 many body states and determine the ground state.  then say that state determines who is where.  Otherwise one says that P(x1) = P(x1|x2,x3,,x5,,y1,y2,y3,y4)P(x1,x2,..y4).  Technically (and importantly) the question which sites are occupied by who is ill-posed.  As in the statistical framework all sites are occupied by everyone, only in the limit of T -> 0 does the problem become deterministic (and even in this regime if two species have the same affinity there still well be residual entropy (uncertainty who is where), at finie T one can only give the probability the site is occupied by x or y, (i.e. one says the prob this site is  occupied by x or y not, one does not say the site <is> occupied.

With this being said, the ill-posed question actually is a very important question, as it gets at the idea of structure.  if one thinks of the proteins on the dna as a structure one gross feature of this is the position of the proteins on the dna.  This can possibly be inferred by knowledge of the probabilities of each many body state (proteins bound on the genome).  However, the occupancy is on the order of ns, so the prob of any particular site would be verry low (one should integrate over cell time, and ask was it likely this site was occupied during this time window). 

This is our aim.  What we do is use expression as our primary indicator of occupancy.  As for expression to occur Pol2 must have been occupying the gene of interest.
And since our proteins recruit (stick to) Pol2 we can infer that our proteins were there too.  

Hence using known targets of Dorsal we reconstruct a genome that is only specific to Dorsal.  Within those gene's regulatory sequence, we assume that fragments with an occupancy of .5 or more will cause nearby genes to be 'expressed'.  Since many factors can confound the underlying structure which caused the expression when one is in a tissue that has a gradient of dorsal and the target gene, we use the 'border' as the position where dorsal occupancy is approximately .5+-.15 the (border = target gene's expression is 1/2 max gene expression).  

Hence we have a conservation constraint (<N> = .5 ).  This allows us to construct structures which are governed by the Dorsal gradient and the conservation constraint.  Actually if it were simply Dorsal in the network this would be trivial, it just reduces to fd stats.  However, we have dts all interacting,  hence the occupancy of a given site is not only constrained by the Nd the number of Dorsal, but also by Ns and Nt, and by the interactions between these proteins.  Hence the structure or distribution is context dependent, leading to the idea of the 'grammar' of the constructs.

We can work in a framework where we constrain the structures to only 6 dorsal.  (i.e. at most 6 sites of Dorsal can define a fragment of DNA that activates the gene).  This is to say at most 6 Dorsal sites can cause the expression to be .5max.  (this is due to the fact that known constructs contain less than 6 sites) 

The interactions of Dorsal with Twist and Snail, can cause the Dorsal occupancy to be context dependent.  That is a nearby twist may coopertively bind Dorsal or a nearby sanil may antagonize Dorsal binding (or cause Dorsal's activating domain to be repressed) both mechs called short range repression (quenching).  This in turn indicates that which sites within a fragment are necessary for 1/2 expr depend on the gradient of snail and twist too, and snail and twist affinity for sites.  

Hence our model constructs the fragement by adding the lowest energy dorsal site first from the fragement distribution.  Then adding the lowest twist site based on all possible configurations of twist within the context of the dorsal site, then adding snail site in the context of dt (ds adding tw was not found to be distinguishable i.e they commute)
 
Then another dorsal site is added (assuming the occupancy of dorsal is not above 1/2 on the fragment) in the context of the dts, a routine is called to check that d2 does not overlap with d1 (although it may overlap with t1 and s1, and thereby compete with these two for binding).

Then another tw site is added (in the context of ddts).. Then another sn... up to a max of 6 sites..

Once all the structures of all fragments have been constructed the expression is determined by the sigmodal function 1/(1+e(<Nd>wd +<Ns>ws) .  We do not account for twist in expression since twist alone is not sufficient for expression.  Wd and Ws are parameters to be trained.  The parameters in the binding affinity are used in the construction of the grammar, thereby allowing the grammar to be influenced by the parameters.

To actually say one has a particular structure (as opposed to a set of structures), one should clearly just use the best structure.      


*/
    
int SeqAnnotator::annoty4( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  // try this without initializing
////////////////////////////////////////////////////////////////////////
// in neuroectoderm ( tw and dorsal must both be present).////////////////
/////////////////////////////////////////////////////////////////////////
//
//		Neuroectoderm
//
//////////////////////////////////////////////////////////////////////////
//cout << "neuro" << endl;
//if( cons[2] ==0 ) {   
//cout << "neuro" << endl;
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < 2; k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		  if( k == 0 ) {
	double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[0];
	tsites.clear();
	//for( int g = 0; g< sitest[k].size(); g++ ) {
//	cout << "s before wloop " << s.size() << endl;
	while( counter > 0 && s.size() > 0 ){
			//tsites.push_back(siteMax( sitest[0] ));
			counter=0;
			Site tempsite;
			 //double mine=100000;
		//	cout << "tempsite " << tempsite<< endl;
 			for ( int i = 0; i < s.size(); i++ ) {
	
				if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
 				}
			}
	//cout << " name " << names << endl;
	//cout << "s after floop " << s.size() << endl;
//cout << "tempsite " << tempsite<< endl;
	//cout << "counter after loop " << counter << endl;
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites ); 
			//cout << "s after floop " << s.size() << endl;
			/*
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  
			   return 1; 
			}
			else{ tsites.clear();}*/
          } // while
	// cout << " about to pocc" << endl;
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
//cout << " made it throug sort " << endl;
	  p.predictOcc( tsites , 5, cons, ff);
//cout << " made it throught predictocc " << endl;
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			//   cout << " about to return" << endl;
			//   cout << tsites << endl;
			  		
			   return 1; 
			}
			else{ tsites.clear();}
	  }  // if k
	  }  // for k


/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
///////////////////////////////////////////////
///////////////////////////////////////////////

vector< double  > sitest_tw_Z( 0 );
////////////////////
// the following codemakes no sense so we comment it out;
/////////////////// 8/12
 //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
/*
sitest[0].clear();
sitest[1].clear();
//for ( int k = 0; k < motifs.size(); k++ ) {
//cout<< "before sitemax " ;
 for ( int k = 0; k <2; k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
	            sitest[k].push_back( Site( i, 0, k, energy ) );

	          }
	  }
*/
int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
//cout << " aout to n++ :" << endl;
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( previous > .5 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal

	for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
		//sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
//cout << "sitest_tw_Z.size  " << sitest_tw_Z.size() << endl;
//cout << "maxZ " << maxZ( sitest_tw_Z ) << endl;
//cout << " sitest_sitest_allk_alli_Z[k] " << sitest_allk_alli_Z[k]<< endl;
int ind = maxZ( sitest_tw_Z );	
//cout << "ind  for twist, after second push " << ind << endl;
double maxZscore = maxZs( sitest_tw_Z );
//if ( cons[1] != 0 ) {
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();
////////////////////////////////////

///////////////////////// delete the following code when the above is uncommented
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );
///////////////////////////////delete the tsites pushback
//cout << "made it throuh " << endl;
return tsites.size();

}
//} //if cons2 == 0

///////////////////////////////////////////////////////////////////////////////////
//
//        mesoderm
//
///////////////////////////////////////////////////////////////////////////////////







/*

else {
//cout << "meso " << endl;
 typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < motifs.size(); k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		 if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		}
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		 if ( energy <= energyThrs[ k ] ) {	
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		}

	          }
	

 /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		  if( k == 0 ) {

tsites.clear();
			//tsites.push_back(siteMax( sitest[0] ));
tsites.push_back (siteMax( sitest[0] ));


			
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  
			   return 1; }
			else{ tsites.clear();}
		  }  //if k
	  } //for k


/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbestdd = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
//cout << " sitesp pushed k " << kbest ;
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
}
else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
}
sitesp[2].push_back (siteMax( sitest[2] ));
sitestoverlap( sitest[2], sitesp[2] );
//cout << " size of  sitest  " <<"\t" << sitest.size() << endl;
///////////////////////////////////////////////////////////////////
// initial maximum site has been created along with a temp V that does not contain overlapping sites with maximum site
//////////////////////////////////////////////////////////////////

// repeat siteMax up to 6 times for each motif. using both concentration and affinity and cooperation
int N = 1; // this is the number of sites
double previous;
double current;

do {
vector< double > ff(3); // focc passed by ref to predictocc.
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
/*
vector< vector< Site > > sitespatemp = sitesp;
double Zb = p.predictZ(tsites,cons);
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
//bool delt = DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore);
bool delt = DeltaZ(  sitespatemp ,cons ,  Zb,  p, previous);  // DeltaZ checks to see if Z is sensitive to removal of twist 
if( previous > .5 && delt)
{ 

sitesp[1].pop_back();  // remove last twist site

break ;
 }
*/
/*
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
int ind = maxZ( sitest_allk_alli_Z[k] );	
////////////////////////////////////////////
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  // for k
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( current > .5 )
{ break ; }
  // exit do loop once dl occ is almost 1.
if(N ==6) { 
break; 
} 


//std::sort(sites.begin(), sites.end(), siteSortPredicate); 
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
// sitesorderedtemp has sites ordered different than sitest, however the index of each Z,( in the vector sitest_allk..) is the same as the site index in sitest.

///////////////////////// delete the following code when the above is uncommented
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
///////////////////////////////delete the tsites pushback

return sitest.size();

}// if snail 0

}// else
*/
int SeqAnnotator::annotydorsalold(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  // try this without initializing
////////////////////////////////////////////////////////////////////////
// in neuroectoderm ( tw and dorsal must both be present).////////////////
/////////////////////////////////////////////////////////////////////////
//
//		Neuroectoderm
//
//////////////////////////////////////////////////////////////////////////
//cout << "neuro" << endl;
if( cons[2] ==0 ) {   
//cout << "neuro" << endl;
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < 2; k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		  
	  }  // for k

/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
//sitesp[1].push_back( sitest[1][ibest] );
//sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
///////////////////////////////////////////////
///////////////////////////////////////////////

vector< double  > sitest_tw_Z( 0 );
////////////////////

int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( previous > .2 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
	for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
		//sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
//cout << "sitest_tw_Z.size  " << sitest_tw_Z.size() << endl;
//cout << "maxZ " << maxZ( sitest_tw_Z ) << endl;
//cout << " sitest_sitest_allk_alli_Z[k] " << sitest_allk_alli_Z[k]<< endl;
int ind = maxZ( sitest_tw_Z );	
//cout << "ind  for twist, after second push " << ind << endl;
double maxZscore = maxZs( sitest_tw_Z );
//if ( cons[1] != 0 ) {
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();
////////////////////////////////////

///////////////////////// delete the following code when the above is uncommented
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .2 ) break;
} while ( abs(previous - current ) > .25 );
///////////////////////////////delete the tsites pushback

return tsites.size();


} //if c=0

///////////////////////////////////////////////////////////////////////////////////
//
//        mesoderm
//
///////////////////////////////////////////////////////////////////////////////////

else {

 typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < motifs.size(); k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		 if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		}
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		 if ( energy <= energyThrs[ k ] ) {	
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		}

	          }
	

 /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		 if( k == 0 ) {
	double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[0];
	tsites.clear();
	//cout << "s before wloop to check if s[0] is default site" << s[0] << endl;
	//for( int g = 0; g< sitest[k].size(); g++ ) {
	//cout << "s before wloop " << s.size() << endl;
	while( counter > 0 && s.size() > 0 ){
			//tsites.push_back(siteMax( sitest[0] ));
			counter=0;
			Site tempsite;
			 //double mine=100000;
	//		cout << "temsite " << tempsite << endl;
	//		cout << "counter before loop " << counter << endl;
 			for ( int i = 0; i < s.size(); i++ ) {
	//			cout << " s energy i " << s[i].energy << endl;
				if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
	//			cout << "tempsite " << tempsite << endl;
 				}
			}
	//cout << "s after floop " << s.size() << endl;
	//cout << "counter after loop " << counter << endl;
			if( counter != 0) {
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites );
			} 
	//		cout << "s after floop " << s.size() << endl;
			
          } // while
	 //cout << " about to pocc" << endl;
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
//cout << " made it throug sort " << endl;
	  p.predictOcc( tsites , 5, cons, ff);
//cout << " made it throught predictocc " << endl;
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
//			   cout << " about to return" << endl;
			   cout << tsites << endl;
			  s.clear();		
			   return 1; 
			}
			else{ s.clear(); tsites.clear();}
	  }  // if k
	  } //for k


/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[0][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbestdd = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
//cout << " sitesp pushed k " << kbest ;
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
}
else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
}
sitesp[2].push_back (siteMax( sitest[2] ));
sitestoverlap( sitest[2], sitesp[2] );
//cout << " size of  sitest  " <<"\t" << sitest.size() << endl;
///////////////////////////////////////////////////////////////////
// initial maximum site has been created along with a temp V that does not contain overlapping sites with maximum site
//////////////////////////////////////////////////////////////////

// repeat siteMax up to 6 times for each motif. using both concentration and affinity and cooperation
int N = 1; // this is the number of sites
double previous;
double current;

do {
vector< double > ff(3); // focc passed by ref to predictocc.
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}

N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
int ind = maxZ( sitest_allk_alli_Z[k] );	
////////////////////////////////////////////
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  // for k
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( current > .5 )
{ break ; }
  // exit do loop once dl occ is almost 1.
if(N ==6) { 
break; 
} 


//std::sort(sites.begin(), sites.end(), siteSortPredicate); 
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
// sitesorderedtemp has sites ordered different than sitest, however the index of each Z,( in the vector sitest_allk..) is the same as the site index in sitest.

///////////////////////// delete the following code when the above is uncommented
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
///////////////////////////////delete the tsites pushback

return sitest.size();

//}// if snail 0
}// else

//*/

}

int SeqAnnotator::annotydorsal(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
//cout << "meso " << endl;
//cout << " motifs.size " << motifs.size() <<  endl;
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  // try this without initializing
////////////////////////////////////////////////////////////////////////
// in neuroectoderm ( tw and dorsal must both be present).////////////////
/////////////////////////////////////////////////////////////////////////
//
//		Neuroectoderm
//
//////////////////////////////////////////////////////////////////////////
//cout << "neuro" << endl;
if( cons[2] ==0 ) {   
    typedef  vector< Site > sitestt;
    sitestt te;	
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    for ( int k = 0; k < 2; k++ ) {
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	    int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }
                   } // for i
		     if( k == 0 ) {
			for ( int i = 0; i < seq.size(); i++ ) {
			      	   int l = d2.length();
			 	    if ( i + l > seq.size() ) break ;
			 	    double energy;
				    Sequence elem( seq, i, l, 1 );
			 	    energy = d2.energy( elem );
				    if ( energy <= energyThrs[ 0 ] ) {	
			 	    sitest[0].push_back( Site( i, 1, 0, energy ) );
				    }
		     		    Sequence rcElem( seq, i, l, 0 );
		      	      	    energy = d2.energy( rcElem );
				    if ( energy <= energyThrs[ 0] ) {
				    sitest[0].push_back( Site( i, 0, 0, energy ) );
				    }
			  }  // for i
			  double mine =100;
			  int gbest;
			  int counter=1;
			  sitestt  s = sitest[k];
			  tsites.clear();
			  while( counter > 0 && s.size() > 0 ){
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					}
			  } // while
			  mine =100;
			  gbest = 100;
			  counter=1;
			  while( counter > 0 && s.size() > 0 ){  // in case d2 also has a set of motifs repeat while loop, making sure not to double count first motif of dl.
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					} 
			  } // while
			  s.clear();
			} // if k == 0
			else{
				double mine =100;
				int gbest;
				int counter=1;
				sitestt  s = sitest[k];
				while( counter > 0 && s.size() > 0 ){
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					
					} 
			  	} // while
				s.clear();
			} // else
	} //for k
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
	p.predictOcc( tsites , 5, cons, ff);
	if ( ff[0] >.15) {  
		return 1; 
	}
	else{ 
//for(int i=0; i<tsites.size(); i++){
//sitesp[0].push_back( tsites[i] );
//	sitestoverlap( sitest[0], sitesp[0] ); 
//}
tsites.clear();}
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
////////////////////// find way to run annotydorsal without twist motif, to check if twist is necessary for network
///  next 50ish lines up to 'do' are choosing best twdl form, it is necessary to do the initial twdl form outside of do loop, inside do, simply repeats procedure adding forms
	SiteVec sitesorderedtemp(0);
	double Ztemp=0;
	int ktempindex=0;
	int itempindex=0;
	double Zbest=0;
	int ibest =0;
	int  kbest=0;
/*
	for ( int k = 0; k < sitest[0].size(); k++ ) {
			ktempindex = k;
	     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
			sitesorderedtemp.push_back( sitest[0][k] );
			sitesorderedtemp.push_back( sitest[1][i] );
			std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
			Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
			sitesorderedtemp.clear();
			itempindex=i;
			if( Zbest < Ztemp ) {
				Zbest = Ztemp;
				ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
				kbest = ktempindex;
			}// if
			}// for i
	} // for k
	sitesp[0].push_back( sitest[0][kbest] );
	sitestoverlap( sitest[0], sitesp[0] ); 
	
	//cout << " sitesp pushed k " << kbest ;
	sitesp[1].push_back( sitest[1][ibest] );
	sitestoverlap( sitest[1], sitesp[1] ); 
*/
	vector< double  > sitest_tw_Z( 0 );
	vector< double  > sitest_tw_Z2( 0 );
	int N=1;
	double previous = 0;
	double current = 0;
//cout << " sitesp " << sitesp << endl;
	do {
		tsites.clear();
		for ( int k = 0; k < 2; k++ ) {
			for( int i = 0; i < sitesp[k].size(); i++){
				tsites.push_back( sitesp[k][i] );
			}
		}
		std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
		p.predictOcc( tsites, 5, cons, ff);
		previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
		if( previous > .15 )
		{ break ; }
		N++;
		Site aaa;
		aaa= siteMax( sitest[0] );
		sitesp[0].push_back( aaa);
		sitestoverlap( sitest[0], sitesp[0] );
		sitestoverlap( sitest[1], sitesp[1] );
		//tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal; otherwise place for loop over all dorsal sites here( changing alg from n+n to n*n, where n is the number of sites
		for ( int k = 0; k < sitest[0].size(); k++ ) {
			for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
				//sitesorderedtemp = tsites;
				sitesorderedtemp.push_back( sitest[0][k] );
				sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
				std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
				sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
				sitesorderedtemp.clear();
			}   // for i
		} // for k
		//for ( int k = 0; k < sitesp[0].size(); k++ ) {
			for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
				sitesorderedtemp.push_back(aaa) ; // = tsites;
				//sitesorderedtemp.push_back( sitesp[0][k] );
				sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
				std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
				sitest_tw_Z2.push_back( p.predictZ( sitesorderedtemp,  cons) );
				sitesorderedtemp.clear();
			}   // for i
		//} // for k
//cout << " sitesp 2 " << sitesp << endl;
		int ind = maxZ( sitest_tw_Z );	
		//cout << " sitest tw Z " << endl << sitest_tw_Z << endl;
		double maxZscore = maxZs( sitest_tw_Z );
		int ind2 = maxZ( sitest_tw_Z2 );	
		//cout << " sitest tw Z " << endl << sitest_tw_Z << endl;
		double maxZscore2 = maxZs( sitest_tw_Z2 );
		//if ( cons[1] != 0 ) {
	if( maxZscore > maxZscore2){
			sitesp[1].push_back( sitest[1][ ind ] ); 
		 	sitestoverlap( sitest[1], sitesp[1] );
	}
	else{
		sitesp[1].push_back( sitest[1][ ind2 ] ); 
		 	sitestoverlap( sitest[1], sitesp[1] );
	}
//cout << " sitesp 3" << sitesp << endl;
		if (N > 5) break;
		tsites.clear();
		///////////////////////// delete the following code when the above is uncommented
		for ( int k = 0; k < 2; k++ ) {
			for( int i = 0; i < sitesp[k].size(); i++){
				tsites.push_back( sitesp[k][i] );
			}
		}
		std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
		p.predictOcc( tsites, 5, cons, ff);
		current = ff[0]; 
		if( current > .15 ) break;
		} 
	while ( current < .15 );    //( abs(previous - current ) > .25 );
///////////////////////////////delete the tsites pushback
//cout << " sitesp " << sitesp << endl;
	return tsites.size();
} //if( cons[2] ==0 )
///////////////////////////////////////////////////////////////////////////////////
//
//        mesoderm
//
///////////////////////////////////////////////////////////////////////////////////
else {
typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < motifs.size(); k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq.size() ) break ;
		 	    double energy;
			    Sequence elem( seq, i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
			 if ( energy <= energyThrs[ k ] ) {	
		 	    sitest[k].push_back( Site( i, 1, k, energy ) );
			}
	     		    Sequence rcElem( seq, i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			 if ( energy <= energyThrs[ k ] ) {	
			    sitest[k].push_back( Site( i, 0, k, energy ) );
			}
	          }
if( k == 0 ) {
for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = d2.length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = d2.energy( elem );
		    if ( energy <= energyThrs[ 0 ] ) {	
         	    sitest[0].push_back( Site( i, 1, 0, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = d2.energy( rcElem );
		    if ( energy <= energyThrs[ 0] ) {
	            sitest[0].push_back( Site( i, 0, 0, energy ) );
		    }
	          }
 /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		// if( k == 0 ) {
//cout << "sitest[0].size() "<< endl;
//cout << sitest[0].size() << endl;
	double mine =100;
	int gbest;
	int counter=1;
	int count=0;
	sitestt  s = sitest[k];
	tsites.clear();
//cout << sitest[k] << endl;
//	cout << "s before wloop to check if s[0] is default site" << s[0] << endl;
	//for( int g = 0; g< sitest[k].size(); g++ ) {
	//cout << "s before wloop " << s.size() << endl;
	while( counter > 0 && s.size() > 0 ){
			//tsites.push_back(siteMax( sitest[0] ));
			counter=0;
			Site tempsite;
			 //double mine=100000;
	//		cout << "temsite " << tempsite << endl;
	//		cout << "counter before loop " << counter << endl;
 			for ( int i = 0; i < s.size(); i++ ) {
	//			cout << " s energy i " << s[i].energy << endl;
				if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
				count++;
			//	cout << "tempsite " << tempsite << endl;
 				}
			}
	//cout << "s after floop " << s.size() << endl;
	//cout << "counter after loop " << counter << endl;
			if( counter != 0) {
/*
			Site temp2 = tempsite;
			if( count > 1) {
				if( temp2.start == tempsite.start) { cout << " infin " << endl; 
cout << " sitekil " << temp2 << endl;

				vector< Site > kil;
				kil.push_back( temp2) ;
				sitestoverlap( s, kil );
				sitestoverlap( sitest[0], kil );
				continue;
 				}
			}
*/
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites );
			} 
		//	cout << "s after floop " << s.size() << endl;
			
          } // while
		s.clear();
		//cout << "seqsitest[0] meso " << sitest[0] << endl; 
	 //cout << " about to pocc" << endl;
	}  // if k
	else{
// cout << " size sitest[k] " << sitest[k].size() << endl;
	if(k==1) continue;
	if( sitest[k].size() == 0 ) {continue; }
		double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[k];
	//tsites.clear();
//	cout << "s before wloop to check if s[0] is default site" << s[0] << endl;
	//for( int g = 0; g< sitest[k].size(); g++ ) {
	//cout << "s before wloop " << s.size() << endl;
	while( counter > 0 && s.size() > 0 ){
			//tsites.push_back(siteMax( sitest[0] ));
			counter=0;
			Site tempsite;
			 //double mine=100000;
	//		cout << "temsite " << tempsite << endl;
			//cout << "counter before loop " << counter << endl;
 			for ( int i = 0; i < s.size(); i++ ) {
			//	cout << " s energy i " << s[i].energy << endl;
				if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
	//			cout << "tempsite " << tempsite << endl;
 				}
			}
	//cout << "s after floop " << s.size() << endl;
//	cout << "counter after loop " << counter << endl;
			if( counter != 0) {
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites );
			} 
	//		cout << "s after floop " << s.size() << endl;
			
          } // while
		s.clear();
	 //cout << " about to pocc" << endl;
	}  // else
	//  }  // if k
	  } //for k

std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
//cout << " made it throug sort " << endl;
	  p.predictOcc( tsites , 5, cons, ff);
//cout << " made it throught predictocc " << endl;
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
//			   cout << " about to return" << endl;
			  // cout << tsites << endl;
			 // s.clear();		
			   return 1; 
			}
			else{ tsites.clear();}
/////////////////////
//for (int i=0; i<motifs.size(); i++ ){
// for ( int k = 0; k < sitest[i].size(); k++ ) {
//cout << "inside newest" ;
//cout << sitest[i][k] << endl ;

//}}
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
} // for k
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
if( sitest[0].size() > 1 ) {
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbestdd = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
} //for k
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) { // cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
//cout << " sitesp pushed k " << kbest << endl;
//cout << " ibest " << ibest << endl;
//sitesp[1].push_back( sitest[1][ibest ] );   // it seems this is a problem.. 9 12 11
//sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
} // if zbes < zbe
else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
//cout << " sitesp pushed k in else 1 =" << kbest << endl ;

sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k in else" << kbest << endl;
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) { // cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
}
}// if sitest0 > 1
/*                      // this is causing double counting, why ?  9 12 11
if( sitest[0].size() ==1 ) {
sitesp[0].push_back( sitest[0][0] );


}
if ( sitest[2].size() != 0 ) {

sitesp[2].push_back (siteMax( sitest[2] ));
sitestoverlap( sitest[2], sitesp[2] );

}
*/
//cout << " size of  sitest  " <<"\t" << sitest.size() << endl;
///////////////////////////////////////////////////////////////////
// initial maximum site has been created along with a temp V that does not contain overlapping sites with maximum site
//////////////////////////////////////////////////////////////////

// repeat siteMax up to 6 times for each motif. using both concentration and affinity and cooperation
int N = 1; // this is the number of sites
double previous;
double current;

do {
vector< double > ff(3); // focc passed by ref to predictocc.
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
//cout << " tsites pushed in do " << endl;
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
/*  // unique was commented out, due to deletion of critical dorsal site in 1PEEt 9 12 11
vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
*/
if(sitest[k].size() > 0  ) {
int ind = maxZ( sitest_allk_alli_Z[k] );	
////////////////////////////////////////////
if ( cons[k] != 0 ) {
//cout << " checking snail cons " << endl;
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
//cout << " ending nail " << endl;
}
}
}  // for k
//cout << "maxZ  " << endl;
tsites.clear();
//cout << "maxZ2  " << endl;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
//cout << " k " << k << endl;
//cout << " i " << i << endl;
//cout << " sitesp[k][i] " <<  sitesp[k][i]  << endl; 
		tsites.push_back( sitesp[k][i] );
//cout << "maxZ3  " << endl;
	}
}
//cout<< " tsites.size() " << tsites.size() << endl;
//cout << "maxZ4 " << endl;
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
//cout << "maxZ5  " << endl;
vector< double > ff2(3);
//cout<< " tsites.size() " << tsites.size() << endl;
//cout << tsites << endl;
p.predictOcc( tsites, 5, cons, ff2);
//cout << "maxZ6 " << endl;

current = ff[0];  // using foccdl in the while condition is out of scope.  why?
//cout << "maxZ7 " << endl;
if( current > .5 )
{ break ; }
  // exit do loop once dl occ is almost 1.
if(N ==6) { 
//cout << "maxZ29  " << endl;
break; 
} 
 //cout << " about to iterate in do " << endl;

//std::sort(sites.begin(), sites.end(), siteSortPredicate); 
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
// sitesorderedtemp has sites ordered different than sitest, however the index of each Z,( in the vector sitest_allk..) is the same as the site index in sitest.

///////////////////////// delete the following code when the above is uncommented
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
///////////////////////////////delete the tsites pushback
//cout << " aobout to return annotydorsal " << endl;
return sitest.size();

//}// if snail 0
}// else

//*/

}


 int SeqAnnotator::annoty2( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  // try this without initializing
////////////////////////////////////////////////////////////////////////
// in neuroectoderm ( tw and dorsal must both be present).////////////////
/////////////////////////////////////////////////////////////////////////
//
//		Neuroectoderm
//
//////////////////////////////////////////////////////////////////////////
cout << "neuro" << endl;
if(cons[2] ==0 ) {   
//cout << "neuro" << endl;
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < 2; k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		  if( k == 0 ) {

tsites.clear();
			//tsites.push_back(siteMax( sitest[0] ));
tsites.push_back (siteMax( sitest[0] ));


			
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  
			   return 1; }
			else{ tsites.clear();}
		  }
	  }


/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
///////////////////////////////////////////////
///////////////////////////////////////////////

vector< double  > sitest_tw_Z( 0 );
////////////////////
// the following codemakes no sense so we comment it out;
/////////////////// 8/12
 //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
/*
sitest[0].clear();
sitest[1].clear();
//for ( int k = 0; k < motifs.size(); k++ ) {
//cout<< "before sitemax " ;
 for ( int k = 0; k <2; k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
	            sitest[k].push_back( Site( i, 0, k, energy ) );

	          }
	  }
*/
int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( previous > .5 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal

	for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
		//sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
//cout << "sitest_tw_Z.size  " << sitest_tw_Z.size() << endl;
//cout << "maxZ " << maxZ( sitest_tw_Z ) << endl;
//cout << " sitest_sitest_allk_alli_Z[k] " << sitest_allk_alli_Z[k]<< endl;
int ind = maxZ( sitest_tw_Z );	
//cout << "ind  for twist, after second push " << ind << endl;
double maxZscore = maxZs( sitest_tw_Z );
//if ( cons[1] != 0 ) {
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();
////////////////////////////////////

///////////////////////// delete the following code when the above is uncommented
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );
///////////////////////////////delete the tsites pushback

return tsites.size();


}

///////////////////////////////////////////////////////////////////////////////////
//
//        mesoderm
//
///////////////////////////////////////////////////////////////////////////////////









else {
//cout << "meso " << endl;
 typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	  for ( int k = 0; k < motifs.size(); k++ ) {
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		 if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		}
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		 if ( energy <= energyThrs[ k ] ) {	
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		}

	          }
	

 /// if coopertivity is not important than this should be triggered before coop is exercised in next loop.
		  if( k == 0 ) {

tsites.clear();
			//tsites.push_back(siteMax( sitest[0] ));
tsites.push_back (siteMax( sitest[0] ));


			
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  
			   return 1; }
			else{ tsites.clear();}
		  }  //if k
	  } //for k

/////////////////////}
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbest = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
//cout << " Ztemp " << Ztemp << endl;
//cout << " Zbest " << Zbest << endl;
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			kbestdd = ktempindex;
		}// if
		}// for i
//cout << " not inside Zbest " << endl;
}
//cout << " k and ibest " << kbest << '\t' << ibest << endl;
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
//cout << " sitesp pushed k " << kbest ;
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
//cout << " sitesp pushed i " << ibest ;
}
else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
//cout << " sitesp pushed k " << kbest ;
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
			  tsites= sitesp[0];
			   return 1; }
}
//sitesp[2].push_back (siteMax( sitest[2] ));  // commented 9 12 11
//sitestoverlap( sitest[2], sitesp[2] );
//cout << " size of  sitest  " <<"\t" << sitest.size() << endl;
///////////////////////////////////////////////////////////////////
// initial maximum site has been created along with a temp V that does not contain overlapping sites with maximum site
//////////////////////////////////////////////////////////////////

// repeat siteMax up to 6 times for each motif. using both concentration and affinity and cooperation
int N = 1; // this is the number of sites
double previous;
double current;

do {
vector< double > ff(3); // focc passed by ref to predictocc.
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
/*
vector< vector< Site > > sitespatemp = sitesp;
double Zb = p.predictZ(tsites,cons);
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
//bool delt = DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore);
bool delt = DeltaZ(  sitespatemp ,cons ,  Zb,  p, previous);  // DeltaZ checks to see if Z is sensitive to removal of twist 
if( previous > .5 && delt)
{ 

sitesp[1].pop_back();  // remove last twist site

break ;
 }
*/
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   // for i
vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
int ind = maxZ( sitest_allk_alli_Z[k] );	
////////////////////////////////////////////
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  // for k
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  // using foccdl in the while condition is out of scope.  why?
if( current > .5 )
{ break ; }
  // exit do loop once dl occ is almost 1.
if(N ==6) { 
break; 
} 


//std::sort(sites.begin(), sites.end(), siteSortPredicate); 
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
// sitesorderedtemp has sites ordered different than sitest, however the index of each Z,( in the vector sitest_allk..) is the same as the site index in sitest.

///////////////////////// delete the following code when the above is uncommented
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
///////////////////////////////delete the tsites pushback

return sitest.size();

}// if snail 0
}// else


/*
add .5 to all expressions (this includes dorsal twist and rho and other targets) except snail!  subtract .5 from all expression except snail,..  
add 1 to all expression, subtract one, go crazy and add 1.5. and subtract it..   How does this noise effect parameters..  this doesn't quite make sense since anything that was zero now becomes one?  it really only makes sense to add the derivative to all the profiles, so onlyu the borders have noise, this results in simply shifting the indices vector by one or two.. which means this vector must also find the transition region for dorsal (and hence twist ). ..

this is equivalent to shifting the transition indices up and down by one or two cells, the difference is that some profiles may have huge shifts if one simply adds .5 to every cell in the profile because of profiles like a constant .4 in the mesoderm, that transforms the

*/
int SeqAnnotator::annoty3( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites) 
{

//for (int i=0; i<nrow;i++) {
int ncol =e.nCols();
int i=ti;
	vector< double > reD;					 
	reD = e.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  

//cout << " mi " << names << " "  <<  mi << endl;
//if( ti <10 ) {
vector< double > cons = f.getCol( column );
tsites.clear();
vector< double > ff(3);  // try this without initializing

       typedef  vector< Site > sitestt;
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	{ int k =2;
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }
		  }
	  
double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[k];
sitestt  sp;
sp.clear();
Site sns;
for( int i =0; i < seqSites.size() ; i++ ) {
	sns = seqSites[i];
	if( sns.factorIdx == 2 ) {
	  sp.push_back( sns );
	}
}
sitestoverlap( sitest[2], sp );
sitestoverlap( s, sp );
sp.clear();
	//tsites.clear();
	cout << "s before wloop to check if s[0] is default site" << s[0] << endl;
	//for( int g = 0; g< sitest[k].size(); g++ ) {
	cout << "s before wloop " << s.size() << endl;
	while( counter > 0 && s.size() > 0 ){
			//tsites.push_back(siteMax( sitest[0] ));
			counter=0;
			Site tempsite;
			 //double mine=100000;
	//		cout << "temsite " << tempsite << endl;
	//		cout << "counter before loop " << counter << endl;
 			for ( int i = 0; i < s.size(); i++ ) {
	//			cout << " s energy i " << s[i].energy << endl;
				if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
				//cout << "tempsite " << tempsite << endl;
 				}
			}
	//cout << "s after floop " << s.size() << endl;
	//cout << "counter after loop " << counter << endl;
			if( counter != 0) {
			//tsites.push_back (tempsite);
//cout << "tempsite " << tempsite << endl;
			sp.push_back (tempsite);
			sitestoverlap( s, sp );
			} 
	//		cout << "s after floop " << s.size() << endl;
			
          } // while
		s.clear();
	cout << " sp.size " << sp.size() <<endl;
for( int i=0; i< seqSites.size(); i++){
tsites.push_back( seqSites[i] );
}
for( int i=0; i< sp.size(); i++){
tsites.push_back( sp[i] );
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
cout << " tsites  " << tsites << endl;
cout << " sp.sizee " << sp.size() << endl;
cout << " sp " << sp << endl <<" sitest "<< sitest <<endl;
//cout << " made it throug tsites "<< names << endl << tsites << endl;
//cout << " sp ::  " << sp << endl;
double	pre=p.predictExpr( tsites, 5, cons) ;
//next line made 3/10/2015 needs to be deleted if modifications are made to predictExpr (which should be returning dorsal occupancy)
//if(e(ti,column) -e(ti,column+2) <= 0 ) {cout << "this works " << endl;sp.clear(); return 1;} //|| e(ti,column)-e(ti,column+2)
// for mesoderm genes the if statement below should not get executed.., since mi will then be 1 or 0;
if(e(ti,column) -e(ti,mi) == 0 ) {sp.clear(); return 1;} //|| e(ti,column)-e(ti,column+2)
//if( abs(pre - e(ti,column) ) < .2 )
if( pre < .5 )
{sp.clear(); return 1; }
		
			else{   sitestoverlap( sitest[2],sp ) ; sp.clear(); } // tsites.clear();}	
 //p.predictOcc( tsites , 5, cons, ff);
cout << " made it throught predictocc " << endl;
//			if ( ff[0] <.5) {  //cout << " ff[0] " << ff[0] << endl;
//			   cout << " about to return" << endl;
			  // cout << tsites << endl;
			 // s.clear();		
			//   return 1; 
			//}
			//else{ }

	
	} // k
//if ( 1 ) { return 1; }
//sitestoverlap( sitest[k],sp ) ; 
//cout<< "\n"  << "tsites " << tsites << endl;
//cout << "\n" << "sitest snail " << "\n" << sitest[2] << endl;
/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
//SiteVec sitesorderedtemp = seqSites;
//double Ztemp=0;
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//  the following 50 or so lines of code are to initialize the do loop  ( this requires using seqSites, sitesp ),  sitesp holds the initial site for the do loop
//////////////////////////////////////////////////// // changed seqSites to tsites line 2780, 11/17/11
////////////////////////////////////////////////////
int ktempindex=0;
//int itempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
//int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;  // if this is seqSites, we need to erase tsites before this point, and not kill any of sitest[2]
		ktempindex = k;
     		  //for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[2][k] );
		//sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;    //e(ti,column) );

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			
		}// if
		
}
//if( !sitest[2].empty() ){
sitesp[2].clear();
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 
//}
//cout << " sitesp pushed i " << ibest ;
///////////////////////////////////////////////
///////////////////////////////////////////////

int N=1;
double previous = 0;
double current = 0;
if( !sitest[2].empty())  {
do {
//cout << " start do " << N << endl;
//tsites.clear();
//tsites = seqSites;
if( N == 1 ) {
 { int k = 2;
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
}// if N
//cout << "tsites " << endl << tsites << endl;
N++;
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
previous=p.predictExpr( tsites, 5, cons) ;
//cout << " tsites inside do " << "/n" << tsites << endl;
//if( abs(previous - e(ti,column) ) < .2 )
if( previous < .5 )
{ break ; }
if( previous < .2 )
{ break ; }


//if (N > 6) break;

//if ( 1 ) { return 1; }
if( !sitest[2].empty() ){
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;
		ktempindex = k;
     		  //for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[2][k] );
		//sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp-  0 );

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			
		}// if
		
}
tsites.push_back( sitest[2][kbest] ) ;

////////////////
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
//cout << "tsites after pushed kbest " << endl << tsites <<  endl;
//cout << "sitest after pushed kbest into tsites " << endl<< sitest[2] << endl;
current = p.predictExpr(tsites,5,  cons); 
if( current < .5 ) break;
//if( abs(current - e(ti,column) ) < .2 ) break;  // this and next line commented out 3/10/2015
//if( current  < .2 ) break;
sitesp[2].clear();
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] );
//cout << "sitest after being deleted by sitesp " << endl<< sitest[2] << endl;
//cout << "sitesp after deleteing kbest from sitest " << endl<< sitesp[2] << endl;
} // if sitest[2] not empty
else break;
/*
if( sitest[2].empty() ){
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 
}
{ int k = 2;
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate); 
*/
} while( N < 3  ); // !sitest[2].empty() );   //abs(previous - current ) > .25 );
} // if sitest2 empty
///////////////////////////////delete the tsites pushback
//}
//cout << " made it to the endd " << endl;
//cout << "sitest after while " << endl<< sitest[2] << endl;
//cout << "tsites after while " << endl<< tsites << endl;
// Delete1( tsites,seqSitesm1delete1 );
//Delete1( tsites,dd );  // this is for robustness analysis, checking to see how much the borders shift
return tsites.size();
}

int SeqAnnotator::annotyd( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites, vector< SiteVec >& dd) 
{

//if( ti <10 ) {
vector< double > cons = f.getCol( column );
tsites.clear();
vector< double > ff(3);  // try this without initializing

       typedef  vector< Site > sitestt;
       vector< sitestt > sitest( motifs.size() );    // sitestemp
       vector< sitestt > sitesp( motifs.size() );    // sitespermanent
	{ int k =2;
		//cout << sitest.size() << " sitest.size " << "\n " << endl; 
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		   // if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		   // }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		  //  if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		  //  }
		  }
	  }

//cout << sitest[2] << endl;
/////////////////////
//////////////////////8/12/11  if the seq is length n then choose 2 out of the n that yield lowest Z.
//SiteVec sitesorderedtemp = seqSites;
//double Ztemp=0;
int ktempindex=0;
//int itempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
//int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = seqSites;
		ktempindex = k;
     		  //for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[2][k] );
		//sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;    //e(ti,column) );

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			
		}// if
		
}
//cout<< " e and names "  << e(ti,column)<<'\t' << names << endl;
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 
//cout << " sitesp pushed i " << ibest ;
///////////////////////////////////////////////
///////////////////////////////////////////////

//cout << " abs .1 " <<abs( .1 ) << endl;
//cout << " abs -.1 " <<abs( -.1 ) << endl;

int N=1;
double previous = 0;
double current = 0;

do {

tsites.clear();
tsites = seqSites;

 { int k = 2;
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
previous=p.predictExpr( tsites, 5, cons) ;

if( abs(previous - e(ti,column) ) < .2 )
{ break ; }
if( previous < .2 )
{ break ; }
N++;

//if (N > 6) break;



for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;
		ktempindex = k;
     		  //for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[2][k] );
		//sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;// e(ti,column) );

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
			
		}// if
		
}
tsites.push_back( sitest[2][kbest] ) ;

////////////////
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
current = p.predictExpr(tsites,5,  cons); 
if( abs(current - e(ti,column) ) < .2 ) break;
if( current  < .2 ) break;
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 

} while( N < 3 || abs(previous - current ) > .25 );
///////////////////////////////delete the tsites pushback
//}
// Delete1( tsites,seqSitesm1delete1 );
Delete1( tsites,dd );  // this is for robustness analysis, checking to see how much the borders shift
return tsites.size();
}

bool SeqAnnotator::DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore)
{
vector< double > f(0);
vector< Site > tsitess(0);
tsitess.clear();
sitespa[1].pop_back();  // remove last twist site
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitespa[k].size(); i++){  
		tsitess.push_back( sitespa[k][i] );
	}
}
std::sort(tsitess.begin(), tsitess.end(), siteSortPredicate);
double Zafter = p.predictZ( tsitess, cons);
p.predictOcc(tsitess,5,cons,f);
double Oafter = f[0];
//current = ff[0];  // using foccdl in the while condition is out of scope.  why?
if ( (Zbefore - Zafter)/Zbefore < .1  && (Obefore - Oafter)/Obefore < .1) { sitespa.clear(); return true; }
sitespa.clear();
return false; 
}


/////////////////////////////////////////////////////

int SeqAnnotator::Delete1( vector< Site >& t,vector< vector< Site > >& tt )
{



int size = t.size();

	//int j = 0;
	//for ( int j = 0; j < sitest.size(); i++ ) {
	// for( ptr2 ; ptr2 < sitest.end(); ptr2++){
//vector< SiteVec > seqSitesm1delete;
int count=0;
for( int i = 0; i < size; i++){  
		SiteVec temp1;
		temp1=t;
		vector< Site >::iterator ptr2 = temp1.begin();
		ptr2 += i;
		//cout << "i*sizeof  " << ptr2 + i*sizeof(Site);
		//SiteVec temp;
		//temp2=t;
		//Site ts = *ptr2;
		temp1.erase(ptr2  );
		tt.push_back( temp1 );
		//temp1.insert(ptr2,ts);
	 	//ptr2++;
		
}
	
 


//tt= seqSitesm1delete;
return 1;
}
///////////////////////////////////////////////////////////////
int SeqAnnotator::maxZ( vector< double >& Z )
{
 int tempindex = 0 ;
 double maxZ=(*min_element(Z.begin(), Z.end() ) );
 for ( int i = 0; i < Z.size(); i++ ) {
	
	if ( maxZ < Z[i] ){
		maxZ = Z[i];
		tempindex = i;
 	}
  }
return tempindex;
}
double SeqAnnotator::maxZs( vector< double >& Z )
{
 int tempindex = 0 ;
 double maxZ=(*min_element(Z.begin(), Z.end() ) );
 for ( int i = 0; i < Z.size(); i++ ) {
	
	if ( maxZ < Z[i] ){
		maxZ = Z[i];
		tempindex = i;
 	}
  }
return maxZ;
}
void SeqAnnotator::sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp )
{
for ( int i = 0; i < sitesp.size(); i++ ) {
	
        vector< Site >::iterator ptr2 = sitest.begin();
	int j = 0;
	//for ( int j = 0; j < sitest.size(); i++ ) {
	// for( ptr2 ; ptr2 < sitest.end(); ptr2++){
	while( ptr2 != sitest.end() ){
		if (siteOverlap( sitesp[i], sitest[j], this->motifs) )  { // overlap returns true when they DO overlap
		//vector< Site >::iterator ptr2 = &sitest[j];   // overlap returns false when the DONT overlap
		sitest.erase(ptr2);
	//this->printPar2();
	//	cout << *ptr2 << " this is the value pointed to by ptr2 after the erase(ptr2) operation."<< endl;
		continue;  // don't increment j, if we erase the site at position j, because the sitest subvector shifts such that a new site fills site j.		
		}  // if j = 0, has a 10 bp overlapping site, then we don't want to increment j until the ptr2 loop executes 10 times, for each iteration j remains 0.
		j++;
		ptr2++;
	}
	
 }

}


Site SeqAnnotator::siteMax( vector< Site >& sitess)
{
 Site tempsite;
 double mine=100000;
 for ( int i = 0; i < sitess.size(); i++ ) {
	
	if ( mine > sitess[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
		mine = sitess[i].energy;
		tempsite = sitess[i];
 	}
  }
/*
 vector< Site >::iterator ptr2 = sites.begin();
	int j = 0;
	//for ( int j = 0; j < sitest.size(); i++ ) {
	// for( ptr2 ; ptr2 < sitest.end(); ptr2++){
	while( ptr2 < sites.end() ){
		//if (siteOverlap( sitesp[i], sitest[j], this->motifs) )  { 
		//vector< Site >::iterator ptr2 = &sitest[j];
		if( maxe < ptr2->energy ){
			maxe = ptr2->energy;
			temp
		sitest.erase(ptr2);
		cout << *ptr2 << " this is the value pointed to by ptr2 after the erase(ptr2) operation.";
		continue;  // don't increment j, if we erase the site at position j, because the sitest subvector shifts such that a new site fills site j.		
		}  // if j = 0, has a 10 bp overlapping site, then we don't want to increment j until the ptr2 loop executes 10 times, for each iteration j remains 0.
		j++;
		ptr2++;
	}
*/
 return tempsite;	
}
int SeqAnnotator::compEnergy( const Sequence& seq, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ ) {
        Sequence elem( seq, sites[i].start, motifs[sites[i].factorIdx].length(), sites[i].strand );
        sites[i].energy = motifs[sites[i].factorIdx].energy( elem );
        sites[i].wtRatio = exp( sites[i].energy-10 );  // 1115 changed - to + in exponential
    }

    return 0;
}


ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;
    if ( toupperStr ( modelOptionStr ) == "BINS" ) return BINS;
    cerr << "modelOptionStr is not a valid model option" << endl; 
    exit(1);
}

string getModelOptionStr( ModelType modelOption )
{
    if ( modelOption == LOGISTIC ) return "Logisitic";
    if ( modelOption == DIRECT ) return "Direct";
    if ( modelOption == QUENCHING ) return "Quenching";
    if ( modelOption == CHRMOD_UNLIMITED ) return "ChrMod_Unlimited";
    if ( modelOption == CHRMOD_LIMITED ) return "ChrMod_Limited";
    if ( modelOption == BINS ) return "Bins";
    return "Invalid";
}

string getIntOptionStr( FactorIntType intOption )
{
    if ( intOption == BINARY ) return "Binary";
    if ( intOption == GAUSSIAN ) return "Gaussian";
    if ( intOption == BINSF ) return "Binsf";
    return "Invalid";
}

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;

    cerr << "objOptionStr is not a valid option of objective function" << endl; 
    exit(1);
}

string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";

    return "Invalid";
}

string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}
double FactorIntFuncBinsf::compFactorInt( double normalInt, double dist, bool orientation )  const
{
     assert( dist >= 0 );
     return normalInt;
}

double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;	
}

double FactorIntFuncGaussian::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );

    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );    
}

double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) ); 
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;
}
//ExprPar::ExprPar( int _nFactors, IntMatrix& _coopMat, const vector< bool >& _repIndicators): coopMat( _coopMat ), repIndicators( _repIndicators ) {}
 
ExprPar::ExprPar( int _nFactors, vector< vector< vector<double> > > _theV , vector< vector< vector<double> > > _theVr): factorIntMat() , theV(_theV),theVr(_theVr)  // call to Matrix constructor, in the initialization list.  neet to assert that _theV has dimensions of nfactors..
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );       
/*
    for ( int i = 0; i < _nFactors; i++ ) {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
        repEffects.push_back( ExprPar::default_repression );
    }

    basalTxp = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo; 
*/
}
ExprPar::ExprPar( int _nFactors ) : factorIntMat()
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );       

   if (modelOption == BINS ) {theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));
		theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));}

   for ( int i = 0; i < _nFactors; i++ ) {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
       // repEffects.push_back( ExprPar::default_repression );
    }

/*    basalTxp = modelOption == LOGISTIC ? ExprPar::default_basal_Logistic : ExprPar::default_basal_Thermo; 
*/
}
	
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, double _basalTxp ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxp( _basalTxp )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() ); 	
   // assert( txpEffects.size() == maxBindingWts.size() && repEffects.size() == maxBindingWts.size() );
  // assert( txpEffects.size() == maxBindingWts.size() );
}

ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) : factorIntMat()
{	
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
//     assert( pars.size() == ( _nFactors * ( _nFactors + 1 ) / 2 + 2 * _nFactors + 2 ); 
    int counter = 0;
	
    // set maxBindingWts 
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    
//cout << "counter after maxb " <<  counter << '\t' ;
   // set the interaction matrix
if (modelOption == BINS ) {

theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j < i; j++ ) {         // 1109 this was changed from <= to <  due to abortion errors, this menas the offdiags are define below
            if ( coopMat( i, j ) ) {
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			  double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );  //524 restet interaction, was 1, same for factorIntMat
			  theV[i][j][k] = interaction;              
                 }
            }  //if construct
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theV[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i ; j < _nFactors; j++ ) {            // this changed from i+1 to i, as indicated in 1109 comment above for offdiags
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
            theV[i][j][k] =  theV[j][i][k] ;
	     //theV[i][j][k].push_back(theV[j][i][k]) ;
            }       
        }
    }    
//cout << " after thev " << counter << '\t' ;
/////////
theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
      //  for ( int j = 0; j <= i; j++ ) {  // 81611
 	for ( int j = 0; j < i; j++ ) { 
            if ( repIndicators[ i] ) {
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			 double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interactionr ), log( max_interactionr ) ) ) : exp( pars[counter++] );  //524 reste]t interaction, was 1, same for factorIntMat
			  theVr[i][j][k] = interaction;              
                 }
            }  //if construct
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theVr[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
            theVr[i][j][k] =  theVr[j][i][k] ;
	     //theV[i][j][k].push_back(theV[j][i][k]) ;
            }       
        }
    }       
//cout << " after thevr " << counter << '\t' ;
     // set the transcriptional effects
    for ( int i = 0; i < _nFactors; i++ ) {
// 2/24/2015
	double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
                txpEffects.push_back( effect );        
	/*
            double effect = searchOption == CONSTRAINED ?  inverse_infty_transform( pars[counter++],  min_effect_Thermo ,  max_effect_Thermo   ) :  pars[counter++] ;
            txpEffects.push_back( effect ); 
	*/
        
    }


 //double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
   //     basalTxp = basal;
//cout << " counter after txpeff  " << counter << endl;
}

}
ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators , bool a) : factorIntMat()
{	
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
//     assert( pars.size() == ( _nFactors * ( _nFactors + 1 ) / 2 + 2 * _nFactors + 2 ); 
    int counter = 0;
	
    // set maxBindingWts 
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) :  pars[counter++] ;
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    
//cout << "counter after maxb " <<  counter << '\t' ;
   // set the interaction matrix
if (modelOption == BINS ) {

theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j < i; j++ ) {         // 1109 this was changed from <= to <  due to abortion errors, this menas the offdiags are define below
            if ( coopMat( i, j ) ) {
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			  double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : pars[counter++] ;  //524 restet interaction, was 1, same for factorIntMat
			  theV[i][j][k] = interaction;              
                 }
            }  //if construct
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theV[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i ; j < _nFactors; j++ ) {            // this changed from i+1 to i, as indicated in 1109 comment above for offdiags
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
            theV[i][j][k] =  theV[j][i][k] ;
	     //theV[i][j][k].push_back(theV[j][i][k]) ;
            }       
        }
    }    
//cout << " after thev " << counter << '\t' ;
/////////
theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
      //  for ( int j = 0; j <= i; j++ ) {  // 81611
 	for ( int j = 0; j < i; j++ ) { 
            if ( repIndicators[ i] ) {
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			 double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interactionr ), log( max_interactionr ) ) ) : pars[counter++] ;  //524 reste]t interaction, was 1, same for factorIntMat
			  theVr[i][j][k] = interaction;              
                 }
            }  //if construct
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theVr[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
            theVr[i][j][k] =  theVr[j][i][k] ;
	     //theV[i][j][k].push_back(theV[j][i][k]) ;
            }       
        }
    }       
//cout << " after thevr " << counter << '\t' ;
     // set the transcriptional effects
    for ( int i = 0; i < _nFactors; i++ ) {
// 2/24/2015
	double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : pars[counter++] ;
                txpEffects.push_back( effect );        
	/*
            double effect = searchOption == CONSTRAINED ?  inverse_infty_transform( pars[counter++],  min_effect_Thermo ,  max_effect_Thermo   ) :  pars[counter++] ;
            txpEffects.push_back( effect ); 
	*/
        
    }


 //double basal = searchOption == CONSTRAINED ? inverse_infty_transform( pars[counter++], min_basal_Logistic, max_basal_Logistic ) : pars[counter++];
   //     basalTxp = basal;
//cout << " counter after txpeff  " << counter << endl;
}

}

void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
//cout << "actindicators" << actIndicators[0] << actIndicators[1]  << actIndicators[2] << endl;
//cout << "repindicators" << repIndicators[0] <<  repIndicators[1]  << repIndicators[2] << endl;

   		//cout << " in getfrepars  pars.size =  "  << pars.size() << endl;
    // write maxBindingWts
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }


	for ( int i = 0; i < nFactors(); i++ ) {
       	 for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
		for(int k=0; k< ExprPar::nbins ; k++){
			
               double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] ); 
	       pars.push_back( interaction);
		             
                 }

            }  //if construct
          }
        }
//cout << "in getfreepars size of par after coopmat " << pars.size() << endl;
        for ( int i = 0; i < nFactors(); i++ ) {
          for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
            if ( repIndicators[i] ) {
		for(int k=0; k< ExprPar::nbins ; k++){
         		 double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
	  	//	cout << "theVr" << "\t" << interaction ; 
			  pars.push_back( interaction);
		     //}    // 520          
                 }

             }  //if construct
           }
         }    
//cout << "in getfreepars size of par after theVr  " << pars.size() << endl;
	 for ( int i = 0; i < nFactors(); i++ ) {
//2/24/2015
 double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
                pars.push_back( effect );
/*
 double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
            pars.push_back( effect );
*/
	  }

//////
//cout << "in getfreepars size of par " << pars.size() << endl;
/////
}

void ExprPar::getFreePars3( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
//cout << "actindicators" << actIndicators[0] << actIndicators[1]  << actIndicators[2] << endl;
//cout << "repindicators" << repIndicators[0] <<  repIndicators[1]  << repIndicators[2] << endl;

   		//cout << " in getfrepars  pars.size =  "  << pars.size() << endl;
    // write maxBindingWts
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : maxBindingWts[i] ;
            pars.push_back( weight );
        }
    }


	for ( int i = 0; i < nFactors(); i++ ) {
       	 for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
		for(int k=0; k< ExprPar::nbins ; k++){
			
               double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) :  theV[i][j][k] ; 
	       pars.push_back( interaction);
		             
                 }

            }  //if construct
          }
        }
//cout << "in getfreepars size of par after coopmat " << pars.size() << endl;
        for ( int i = 0; i < nFactors(); i++ ) {
          for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
            if ( repIndicators[i] ) {
		for(int k=0; k< ExprPar::nbins ; k++){
         		 double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : theVr[i][j][k] ;
	  	//	cout << "theVr" << "\t" << interaction ; 
			  pars.push_back( interaction);
		     //}    // 520          
                 }

             }  //if construct
           }
         }    
//cout << "in getfreepars size of par after theVr  " << pars.size() << endl;
	 for ( int i = 0; i < nFactors(); i++ ) {
//2/24/2015
 double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : txpEffects[i] ;
                pars.push_back( effect );
/*
 double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
            pars.push_back( effect );
*/
	  }

//////
//cout << "in getfreepars size of par " << pars.size() << endl;
/////
}

void ExprPar::getFreePars2( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) 
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
//cout << "actindicators" << actIndicators[0] << actIndicators[1]  << actIndicators[2] << endl;
//cout << "repindicators" << repIndicators[0] <<  repIndicators[1]  << repIndicators[2] << endl;

   	//	cout << " in getfrepars  pars.size =  "  << pars.size() << endl;
    // write maxBindingWts
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }


	for ( int i = 0; i < nFactors(); i++ ) {
       	 for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
		for(int k=0; k< ExprPar::nbins ; k++){
			
               double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] ); 
	       pars.push_back( interaction);
		             
                 }

            }  //if construct
          }
        }
//cout << "in getfreepars size of par after coopmat " << pars.size() << endl;
        for ( int i = 0; i < nFactors(); i++ ) {
          for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
            if ( repIndicators[i] ) {
		for(int k=0; k< ExprPar::nbins ; k++){
         		 double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
	  	//	cout << "theVr" << "\t" << interaction ; 
			  pars.push_back( interaction);
		     //}    // 520          
                 }

             }  //if construct
           }
         }    
//cout << "in getfreepars size of par after theVr  " << pars.size() << endl;
	 for ( int i = 0; i < nFactors(); i++ ) {
 double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
            pars.push_back( effect );
	  }

//////
//cout << "in getfreepars size of par " << pars.size() << endl;
/////
}
void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const //, const vector< vector< vector<double > > > _theV ) const
{
//     os.setf( ios::fixed );
//     os.precision( 3 );
    
    // print the factor information
    for ( int i = 0; i < nFactors(); i++ ) {
        os << "motif"<<motifNames[i] << "\t" << "mBwt" << maxBindingWts[i] << "\t" <<"txpE"<< txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t" << repEffects[i];
        os << endl;
    }
/*
    // print the basal transcription
    os << "basal_transcription = " << basalTxp << endl;
  if(modelOption == LOGISTIC){  
    // print the cooperative interactions
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << motifNames[i] << "\t" << motifNames[j] << "\t" << factorIntMat( i, j ) << endl;
        }
    }
*/
if(modelOption== BINS){
     // print the binned parameters
for(int i =0; i < theV.size(); i++){
	for(int j =0; j < theV[i].size(); j++){
		for(int k =0; k < theV[i][j].size(); k++)  //  i don't think this should be k less then vij.size...
		{  if ( coopMat(i,j)){
			os << "theV" << "\t" << theV[i][j][k] ;
			}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		}
	}
}
}

}

int ExprPar::load( const string& file, IntMatrix& coopMat, const vector< bool >& repIndicators )
{
    // open the file
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl;	exit( 1 ); } 
//cout << " starting load " << endl;    
    // read the factor information
    vector< string > motifNames( nFactors() );
    for ( int i = 0; i < nFactors(); i++ ) {
        fin >> motifNames[i] >> maxBindingWts[i] >> txpEffects[i]; 
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) fin >> repEffects[i];
    }

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    /*
    // read the basal transcription
    string symbol, eqSign, value;
    fin >> symbol >> eqSign >> value;
    if ( symbol != "basal_transcription" || eqSign != "=" ) return RET_ERROR;
    basalTxp = atof( value.c_str() );
   */
theV = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
theVr = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
//cout << " made it throught maxB and txp " << endl;
    string value, sym;
	fin >> sym;
	fin >> sym;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            	if ( coopMat( i, j ) ) {
			for(int k=0; k< ExprPar::nbins; k++){
			   fin >> value;  
			  double interaction = atof( value.c_str() ) ; // for more complicated data use factorIdxMap[motifNames[i]] for i,j indices inside for loops
			  theV[i][j][k] = interaction;              
                	 }
            	}  //if construct
                else {
                     for(int k=0; k< ExprPar::nbins; k++){
			  theV[i][j][k] = 1;    
		     }
                 }      
         }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = i + 1; j < nFactors(); j++ ) {
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
            theV[i][j][k] =  theV[j][i][k] ;
            }       
        }
    }    

///////////////////////////////////
for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( repIndicators[ i] ) {
		fin >> sym;
	fin >> sym;
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			fin >> value;  
			  double interaction = atof( value.c_str() ) ; // for more complicated data use factorIdxMap[motifNames[i]] for i,j indices inside for loops
			  theVr[i][j][k] = interaction;              
                 }
            }  //if construct
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theVr[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = i + 1; j < nFactors(); j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
            theVr[i][j][k] =  theVr[j][i][k] ;
            }       
        }
    }       
   

    
  return 0;
}

void ExprPar::adjust()
{
    // adjust binding paramters//
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) maxBindingWts[i] *= 2.0;
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) maxBindingWts[i] /= 2.0;
    }
    // adjust theV matrix
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           //if(factorIntMat(i,j) != 1){
	   
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) ) { 
                        theV[i][j][k] *= 10.0; 
                       // for ( int i = 0; i < nFactors(); i++ ) {
                         //      for ( int j = i + 1; j < nFactors(); j++ ) {
	                   //           for ( int k = 0; k < ExprPar::ExprPar::nbins; k++ ) {
                                           //theV[i][j][k].push_back( theV[j][i][k] );
                                           theV[i][j][k] = theV[j][i][k] ;
                         
                     }
                                         //if ( factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
	                                  // if (factorIntMat(i,j) != 1) {
		                         //for(int k =0; k < ExprPar::nbins; k++){
                     if ( theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) { 
                                          //factorIntMat( i, j ) /= 2.0;
                                          //factorIntMat( j, i ) = factorIntMat( i, j ); 
                          theV[i][j][k] /= 2.0; 
                         //for ( int i = 0; i < nFactors(); i++ ) {
                           //      for ( int j = i + 1; j < nFactors(); j++ ) {
	                     //           for ( int k = 0; k < ExprPar::nbins; k++ ) {
                                            theV[i][j][k] = theV[j][i][k] ; 
                         
                     }
               }
            //}
        }
    }

for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           //if(factorIntMat(i,j) != 1){
	   
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) ) { 
                        theVr[i][j][k] *= 10.0; 
                       // for ( int i = 0; i < nFactors(); i++ ) {
                         //      for ( int j = i + 1; j < nFactors(); j++ ) {
	                   //           for ( int k = 0; k < ExprPar::ExprPar::nbins; k++ ) {
                                           //theV[i][j][k].push_back( theV[j][i][k] );
                                           theVr[i][j][k] = theVr[j][i][k] ;
                         
                     }
                                         //if ( factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
	                                  // if (factorIntMat(i,j) != 1) {
		                         //for(int k =0; k < ExprPar::nbins; k++){
                     if ( theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) { 
                                          //factorIntMat( i, j ) /= 2.0;
                                          //factorIntMat( j, i ) = factorIntMat( i, j ); 
                          theVr[i][j][k] /= 2.0; 
                         //for ( int i = 0; i < nFactors(); i++ ) {
                           //      for ( int j = i + 1; j < nFactors(); j++ ) {
	                     //           for ( int k = 0; k < ExprPar::nbins; k++ ) {
                                            theVr[i][j][k] = theVr[j][i][k] ; 
                         
                     }
               }
            //}
        }
    }

// adjust transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        
            if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) txpEffects[i] *= 2.0;
            if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) txpEffects[i] /= 2.0;
        
        
    }

}

ModelType ExprPar::modelOption = CHRMOD_UNLIMITED;  //   5/19
SearchType ExprPar::searchOption = UNCONSTRAINED;   //   5/19
int ExprPar::estBindingOption = 1;                  //   1. estimate binding parameters; 0. not estimate binding parameters
vector< int > ExprFunc::Bc(5,1);
vector< Sequence > ExprPredictor::seqsy;
vector< string > ExprPredictor::seqNmes;
Matrix ExprPredictor::factorExprData2;
Matrix ExprPredictor::exprData2;
double ExprPar::default_weight = 10;
double ExprPar::default_interaction = 1;
double ExprPar::default_effect_Logistic = 5;
double ExprPar::default_effect_Thermo = 5;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -9.0;
double ExprPar::default_basal_Thermo = 0.0001;
double ExprPar::min_weight = .0099999;		
double ExprPar::max_weight = 100;		
double ExprPar::min_interaction = .99999;	
double ExprPar::max_interaction = 100;
double ExprPar::min_interactionr = .001;	
double ExprPar::max_interactionr = 1.00001;
double ExprPar::min_effect_Logistic = 5;	
double ExprPar::max_effect_Logistic = 5;
// double ExprPar::min_effect_Direct = 0.01;
double ExprPar::min_effect_Thermo = 1; //.01;
double ExprPar::max_effect_Thermo = 100;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Logistic = -9.0;	
double ExprPar::max_basal_Logistic = -9.0;
double ExprPar::min_basal_Thermo =.001;	
double ExprPar::max_basal_Thermo = 0.001;
double ExprPar::delta = 0.0000001;  // was .1 1010
int ExprPar::nbins = 5;  //binwidth;
 // after the colon apparently all these Datamembers are being initialized??, strange that one can initialize maxContact with parenthesis..
//// c.r. page 402 as datamembers can be initialized in the initialization list using parenthesis (e.g. maxContact( x)  = maxContact = x...)

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par,  const vector < int >& _B, const vector < int >& _Br  ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), par( _par ),   B( _B ), Br( _Br )
{ 
    int nFactors = par.nFactors();
    //assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
    Bc = vector< int >(5,1);
}



void ExprFunc::predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs , vector< double >& fOcc)  // k is the current sequence

///////
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[ sites[i].factorIdx ] * sites[i].wtRatio );	
    }
// bindingWts = K*[c]  = [c]*exp( - (Gseq - Gcon) )  = the parameter absorbs ref info and concentration info.
//cout << bindingWts << "bindingwts " << endl;
    // Logistic model


if ( modelOption == BINS ) {


//////////////////////////101jc
// initialization
	vector< double > Y( n + 1 );
    Y[0] = 0;
// 	for ( int i = 1; i <= n; i++ ) Y.push_back( 1 );
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();
//cout << "z " << Z_bind << endl;
	for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
//		cout << "match" << match << endl;
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
//cout << Y_total << "  ytotal  " << endl;
              factorOcc.push_back( Y_total / Z_bind );
	  }
//cout << factorOcc << "  factorOcc  " << endl;
        double totalEffect = 0;
fOcc = factorOcc;

//cout << fOcc << "   fOcc   " << endl;
      for ( int i = 0; i < motifs.size(); i++ ) {
            double effect = factorOcc[i];  // par.txpEffects[i] * factorOcc[i]; // 103  // we need to make sure the order of par.txpEffects and factorOcc are the same for the factors, i.e i=i?
	
            totalEffect += effect;
	}
            // length correction
          //   totalEffect = totalEffect / (double)length;
        //}
        //} //return par.expRatio * logistic( log( par.basalTxp ) + totalEffect );
        //return   totalEffect;    // 103 //logistic( par.basalTxp + totalEffect );
   
}
}
double ExprFunc::predictExpr2( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

//////////////////////////101jc
// initialization
	vector< double > Y( n + 1 );
    Y[0] = 0;
// 	for ( int i = 1; i <= n; i++ ) Y.push_back( 1 );
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();
//cout << " Zbind  " << Z_bind << endl;
	for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		//cout << "match" << match << endl;
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
	//cout << Y_total << endl;
              factorOcc.push_back( Y_total / Z_bind );
	  }

        double totalEffect = 0;

     // for ( int i = 0; i < motifs.size(); i++ ) {
           //double effect = par.txpEffects[i] * factorOcc[i]; // 103  // we need to make sure the order of par.txpEffects and factorOcc are the same for the factors, i.e i=i?
	//totalEffect = factorOcc[0];
       //     totalEffect =par.txpEffects[0] * factorOcc[0] + par.txpEffects[1]*factorOcc[1]  -  par.txpEffects[2] * factorOcc[2] -5;
	//}
//cout << "factorOcc0 inside predictExpr " << endl;
	// totalEffect =par.txpEffects[0] * factorOcc[0]   -5;
	 totalEffect = factorOcc[0]  ;// -5;
//cout << "factorOcc0 " << factorOcc[0] <<endl;
            // length correction
          //   totalEffect = totalEffect / (double)length;
        //}
	//cout << " totalEffect " << totalEffect << endl;
	//cout << " exp " << exp(-totalEffect) << endl;
//cout << " 1/(exp + 1 ) " << 1/(1+exp(-totalEffect) )<< endl;
	if(totalEffect >= .5){  return  totalEffect; }
	else{ return 0;}
 //  1/(1+exp(-totalEffect) );                    //logistic(  totalEffect-1 );
        //return   totalEffect;    // 103 //logistic( par.basalTxp + totalEffect );
   
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

//////////////////////////101jc
// initialization
	vector< double > Y( n + 1 );
    Y[0] = 0;
// 	for ( int i = 1; i <= n; i++ ) Y.push_back( 1 );
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();
//cout << " Zbind  " << Z_bind << endl;
	for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		//cout << "match" << match << endl;
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
	//cout << Y_total << endl;
              factorOcc.push_back( Y_total / Z_bind );
	  }

        double totalEffect = 0;

     // for ( int i = 0; i < motifs.size(); i++ ) {
           //double effect = par.txpEffects[i] * factorOcc[i]; // 103  // we need to make sure the order of par.txpEffects and factorOcc are the same for the factors, i.e i=i?
	//totalEffect = factorOcc[0];
       //     totalEffect =par.txpEffects[0] * factorOcc[0] + par.txpEffects[1]*factorOcc[1]  -  par.txpEffects[2] * factorOcc[2] -5;
	//}
//cout << "factorOcc0 inside predictExpr " << endl;
	// totalEffect =par.txpEffects[0] * factorOcc[0]   -5;
	 totalEffect = factorOcc[0]  ;// -5;
//cout << "factorOcc0 " << factorOcc[0] <<endl;
            // length correction
          //   totalEffect = totalEffect / (double)length;
        //}
	//cout << " totalEffect " << totalEffect << endl;
	//cout << " exp " << exp(-totalEffect) << endl;
//cout << " 1/(exp + 1 ) " << 1/(1+exp(-totalEffect) )<< endl;
	 
       return  totalEffect;  //  1/(1+exp(-totalEffect) );                    //logistic(  totalEffect-1 );
        //return   totalEffect;    // 103 //logistic( par.basalTxp + totalEffect );
   
}
double ExprFunc::predictExpr3( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

//////////////////////////101jc
// initialization
	vector< double > Y( n + 1 );
    Y[0] = 0;
// 	for ( int i = 1; i <= n; i++ ) Y.push_back( 1 );
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();
//cout << " Zbind  " << Z_bind << endl;
	for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		//cout << "match" << match << endl;
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
	//cout << Y_total << endl;
              factorOcc.push_back( Y_total / Z_bind );
	  }

        double totalEffect = 0;

     // for ( int i = 0; i < motifs.size(); i++ ) {
           //double effect = par.txpEffects[i] * factorOcc[i]; // 103  // we need to make sure the order of par.txpEffects and factorOcc are the same for the factors, i.e i=i?
	//totalEffect = factorOcc[0];
            totalEffect =par.txpEffects[0] * factorOcc[0] + par.txpEffects[1]*factorOcc[1]  -  par.txpEffects[2] * factorOcc[2] -5;
	//}
//cout << "factorOcc0 inside predictExpr " << endl;
	// totalEffect =par.txpEffects[0] * factorOcc[0]   -5;
	// totalEffect = factorOcc[0]  ;// -5;
//cout << "factorOcc0 " << factorOcc[0] <<endl;
            // length correction
          //   totalEffect = totalEffect / (double)length;
        //}
	//cout << " totalEffect " << totalEffect << endl;
	//cout << "par.txpEffects[2] "<< par.txpEffects[2] << " factorOcc[2]"<< factorOcc[2]  <<endl;
//cout << "factorOcc0 " << factorOcc[0] <<endl;
//cout << "facotr occ " << factorOcc << endl;
	//cout << " exp " << exp(-totalEffect) << endl;
//cout << " 1/(exp + 1 ) " << 1/(1+exp(-totalEffect) )<< endl;
	 
       return  1/(1+exp(-totalEffect) );                    //logistic(  totalEffect-1 );
        //return   totalEffect;    // 103 //logistic( par.basalTxp + totalEffect );
   
}
double ExprFunc::predictZ( const SiteVec& _sites, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    // store the sequence 
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0 
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );

    }

//////////////////////////101jc
// initialization
	
    
	// compute the partition function of binding
	double Z_bind = compPartFunc();

	
       return Z_bind;
        //return   totalEffect;    // 103 //logistic( par.basalTxp + totalEffect );
   
}

    
double ExprFunc::compPartFunc()
{
    int n = sites.size() - 1;
//    cout << " sites.size = " << sites.size() << endl;
	// initialization
    Z = vector< double >( n + 1 );
    Z[0] = 1.0;
    Zt = vector< double >( n + 1 );
    Zt[0] = 1.0;
	// printPar( par ); 
	// recurrence
double sum=0;  //1109 
	for ( int i = 1; i <= n; i++ ) {
		sum = Zt[boundaries[i]];  //  declaration for sum was here inside for loop
	//cout << sum << "  sum " << endl;
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
			//cout << " compfactorint (ij)  " << compFactorInt( sites[ i ], sites[ j ] )*Z[j] << endl;
			//cout << " compfactorint (ij)  " << compFactorInt( sites[ i ], sites[ j ] ) << endl;	
		}
//cout << bindingWts << "  bindingWts " << endl;
		Z[i] = bindingWts[ i ] * sum;
//cout << "Z[i] = bindingwts[i]*sum = " << Z[i] << endl;
        Zt[i] = Z[i] + Zt[i - 1];
//cout << Zt[i] << "  zt " << endl;
	}
	
	// the partition function 
// 	double Z_bind = 1;
// 	for ( int i = 0; i < sites.size(); i++ ) {
// 		Z_bind += Z[ i ];	
// 	}
    double Z_bind = Zt[n];
//cout << " Zbind " << Z_bind << endl;
	return Z_bind;
}

ModelType ExprFunc::modelOption = QUENCHING;  //5/19





double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
// 	assert( !siteOverlap( a, b, motifs ) );
double dist = abs( a.start - b.start );
//int nbins;
//nbins = int(coopDistThr/binwidth);
double maxInt=1;
//vector< double> bords;
// nbins = mm.size() + 1 hence the number of if statements is mm+1
if ( ( repIndicators[a.factorIdx ] || repIndicators[ b.factorIdx ] ) && ( repIndicators[a.factorIdx ] && repIndicators[ b.factorIdx ] )   )
maxInt = 1;

else if ( repIndicators[ a.factorIdx ] || repIndicators[ b.factorIdx ]) {
	if (      dist <= Br[0] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 0 ]; } //both sides of matrix diagonal need  to have the same parameters
	if (      dist > Br[0] && dist <= Br[1] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 1 ]; }
	if (      dist > Br[1] && dist <= Br[2] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 2 ]; }
	if (      dist > Br[2] && dist <= Br[3] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 3 ]; }
	if (      dist > Br[3] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 4 ]; }

}
else {
	if (      dist <= B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::Bc[0]++; } 
	if (      dist > B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::Bc[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::Bc[2]++;}
	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::Bc[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::Bc[4]++;}

}

     bool orientation = ( a.strand == b.strand ); 
    return  intFunc->compFactorInt( maxInt, dist, orientation );	
}



bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
// 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}




ExprPredictor::ExprPredictor(const vector< vector< SiteVec > >& _seqSitesb, const vector< vector< double > >& _bindingData,  vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr,const vector< int >& _binBord ,const vector< int >& _binBordr, const vector < bool >& _indicator_bool, SeqAnnotator _anny, const string& _file,  vector< SiteVec >& _seqSitesbot,  vector< SiteVec >& _seqSitesm1,  vector< SiteVec >& _seqSitesm2 , vector< SiteVec >& _seqSitesf2 ,vector< SiteVec >& _seqSitesbotf2, vector< SiteVec >& _seqSitesm1f2, vector< SiteVec >& _seqSitesm2f2, vector< SiteVec >& _seqSitesf3 ,  vector< SiteVec >& _seqSitesbotf3,vector< SiteVec >& _seqSitesm1f3 , vector< SiteVec >& _seqSitesm2f3): seqSitesb( _seqSitesb ), bindingData( _bindingData), seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), intFunc( _intFunc ), coopMat( _coopMat ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), binBord(_binBord), binBordr(_binBordr),  indicator_bool ( _indicator_bool ), anny( _anny) , file( _file ) , seqSitesbot( _seqSitesbot ), seqSitesm1( _seqSitesm1 ), seqSitesm2( _seqSitesm2 ),seqSitesf2( _seqSitesf2 ) ,seqSitesbotf2(_seqSitesbotf2), seqSitesm1f2( _seqSitesm1f2) ,seqSitesm2f2( _seqSitesm2f2), seqSitesf3( _seqSitesf3) , seqSitesbotf3( _seqSitesbotf3),seqSitesm1f3( _seqSitesm1f3 ), seqSitesm2f3( _seqSitesm2f3), seqSitesm1d1(),d(),spaceweights(0)
{   



    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
    assert( coopMat.isSquare() && coopMat.isSymmetric() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() );
    assert( maxContact > 0 );
    assert( repIndicators.size() == nFactors() );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors() );
    assert( repressionDistThr >= 0 );

    // set the model option for ExprPar and ExprFunc
    ExprPar::modelOption = modelOption;
    ExprFunc::modelOption = modelOption;
    ExprPar::nbins = _binBord.size() +1;  //  add one since there exist a bin after the last border and before the first border
 //binBord = vector< int>(2,1);
    // set the values of the parameter range according to the model option
    if ( modelOption != LOGISTIC && modelOption != DIRECT && modelOption !=BINS ) {
        ExprPar::min_effect_Thermo = 0.99;
        ExprPar::min_interaction = 0.99;
    }

    // set the option of parameter estimation
    ExprPar::estBindingOption = estBindingOption;
}

double ExprPredictor::objFunc( const ExprPar& par ) const
{
    if ( objOption == SSE ) return compRMSE( par );	
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par ); 
}

double ExprPredictor::objFunc2(  ExprPar& par ) 
{
    return -compAvgCorr2( par );
}

double ExprPredictor::objFuncborder(  ExprPar& par ) 
{
    return compAvgCorrborder8( par );
}
double ExprPredictor::objFuncborder2(  ExprPar& par ) 
{
    return compAvgCorrborder2( par );  // changed to + sign 11/6/2011
}
int ExprPredictor::train( const ExprPar& par_init )
{
//if (!testPar(par_init) ) { 
   par_model = par_init;
//par_model.adjust() ;}
//cout << "entering train print parmodel " << endl;
//printPar(par_init);
//cout<< " adjust bureau " << endl;
 // if ( ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
 //  obj_model = objFuncborder( par_model );  //92511
       if ( nAlternations == 0 ) {obj_model = objFuncborder2( par_model );return 0;}
    // cout << "Random start " <<  1 << ":\tParameters = "; printPar( par_model );
      //  cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
    // alternate between two different methods
////////////////////
 // alternate between two different methods
    ExprPar par_result;
    double obj_result;
   for ( int i = 0; i < nAlternations; i++ ) {
//cout << "about to simplex " << endl;
    //    simplex_minimize( par_result, obj_result );
//if (!testPar(par_result) ) { par_result.adjust() ;}
  //      par_model = par_result; // par_result is declared inside train to keep scoping issues right, otherwise if simplex was a void method the objects created inside simplex would no longer be accessible outside of simplex, by passing the objects as parameters, it is a way to circumvent scoping issues...
        // par_model.adjust();
        gradient_minimize( par_result, obj_result );
//if (!testPar(par_result) ) { par_result.adjust() ;}
        par_model = par_result;
	
      //   par_model.adjust();
    }
	
   
 par_model = par_result; 
    obj_model = obj_result;
//printPar( par_model );
//cout << " the objective function ( - corrbor ) value : " << obj_model << endl;
printPar( par_model );
    return 0;	
}

int ExprPredictor::train4( const ExprPar& par_init )
{
 par_model = par_init;
cout << " inside train4 " << endl;
printPar( par_model );
    ExprPar par_result;
    double obj_result;

        gradient_minimize2( par_result, obj_result );

	
   
 par_model = par_result; 
    obj_model = obj_result;

printPar( par_model );
    return 0;	
}
int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
   
 train( par_init );
    // training with random starts
	ExprPar par_best = par_model;
	double obj_best = obj_model;
	for ( int i = 0; i < nRandStarts; i++ ) {
        	ExprPar par_curr = par_init; 

///////////////////////////////////////////
	//	cout << " print par before  randsamplepar " << endl;
		//printPar( par_curr );
		randSamplePar( rng, par_curr ); /////////////this has a bug 525
		//cout << "print par after randsample par " << endl;
///////////////////////////////////////////////
		//printPar( par_curr );
//testPar(par_curr);
//cout<< "print par after testpar : " << endl;
//printPar( par_curr );
		train( par_curr );
       		//cout << "Random start " << i + 1 << ":\tParameters = " << endl;
		 printPar( par_model );
       // 	cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
		if ( obj_model < obj_best ) {
			par_best = par_model;
			obj_best = obj_model;	
		}
	}    

    // training using the best parameters so far
    if ( nRandStarts ) train( par_best ); 
 //   cout << "Initial training:\tParameters = "; 
//printPar( par_best );
par_model = par_best;
obj_model = obj_best;
cout<< " about to print best par " << endl;
printPar( par_model );
cout <<" ending train rng " << endl;
    return 0;
}
/////////////////////////////////////////////////jacobc
int ExprPredictor::train3( const ExprPar& par_init, const gsl_rng* rng )
{
   par_model = par_init;

    if ( ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
    obj_model = objFunc( par_model );
       if ( nAlternations == 0 ) return 0;


/*



 par_model = par_init;
//cout << "entering train print parmodel " << endl;
//printPar(par_init);
//cout<< " adjust bureau " << endl;
 // if ( ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
 //  obj_model = objFuncborder( par_model );  //92511
       if ( nAlternations == 0 ) return 0;
    // cout << "Random start " <<  1 << ":\tParameters = "; printPar( par_model );
      //  cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
    // alternate between two different methods
////////////////////
 // alternate between two different methods
    ExprPar par_result;
    double obj_result;
   for ( int i = 0; i < nAlternations; i++ ) {
 simplex_minimize2( par_result, obj_result );
        par_model = par_result; // par_result is declared inside train to keep scoping issues right, otherwise if simplex was a void method the objects created inside simplex would no longer be accessible outside of simplex, by passing the objects as parameters, it is a way to circumvent scoping issues...
        // par_model.adjust();
        gradient_minimize2( par_result, obj_result );
        par_model = par_result;
	
      //   par_model.adjust();
    }
	
   
 par_model = par_result; 
    obj_model = obj_result;
      //   par_model.adjust();
   
  // compAvgCorr2( par_model );

*/
    return 0;
}
void ExprPredictor::printFile( ostream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
//    cout.precision( 3 ); 
//     cout.width( 8 ); 
    
    // print binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC || modelOption == DIRECT ) os << par.txpEffects[i] << "\t";
        else {
            if ( actIndicators[i] ) os << par.txpEffects[i] << "\t";
        }
    }
 // print the basal transcription
    os << par.basalTxp << "\t"; 
    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << par.factorIntMat( i, j ) << "\t";
        }
    }

   
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) os << par.repEffects[i] << "\t"; 
        }
    }
    os << jRMSE.getObj() << endl;
   
}
		

///////////////////////////////////////////////////jacobc424
void ExprPredictor::printFilePar_KfoldCV( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "format.tex",ios::app );
int j=0;
cout << seqSites.size() << endl;
vector< Site > tsites(0);

/*
fo << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            fo << par.maxBindingWts[i] << "\t"; 
        }        
    }
fo << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       fo << par.txpEffects[i] << "\t";
       
    }
*/
fo << endl<<  " cooperativity : " << endl;
 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			fo <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
fo << endl;
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
  
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {//cout << " j = " << j << endl;
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();



/////////////////
 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";

fo << ExprPredictor::seqNmes[m] << "&" << "cell"<< cell[m] << "&" ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}







fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();


}


int ExprPredictor::load( const string& fila )
{
   ofstream fout( fila.c_str() );
    vector< string > motifNames;
//cout << "size of motifNames " << motifNames.size() << endl;

string dorsal, twist, snail;
dorsal="dl";
twist="tw";
snail="sn";
        motifNames.push_back( dorsal );
	motifNames.push_back( twist );
	motifNames.push_back( snail );
         fout << setprecision(2) ;
    for ( int i = 0; i < nFactors(); i++ ) {
	
        fout << motifNames[i] << '\t'<< par_model.maxBindingWts[i]<< '\t' << par_model.txpEffects[i] << '\n';
   
    }

// factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }

   
 for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            	if ( coopMat( i, j ) ) {
			fout << motifNames[j]<< '\t' << motifNames[i];
			for(int k=0; k< ExprPar::nbins; k++){
			   fout << '\t' << par_model.theV[i][j][k];               
                	 }
			fout << '\n' ;
            	}  //if construct   
         }
    }

for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if ( repIndicators[ i] ) {
		fout << motifNames[j]<< '\t' << motifNames[i];
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			fout<< '\t'<< par_model.theVr[i][j][k]  ;             
                 }
		fout << '\n' ;
            }  //if construct    
        }
    }
fout.close();
return 0;
}	



void ExprPredictor::printFiled( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;         /// why is this cout here?
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "format.tex",ios::app );
int j=0;
cout << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else { cout << " j1 = " << j << endl;
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	
	foo << endl;
	foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
		//if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
			
			else { cout << " j = " << j << endl;
				

				if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();



/////////////////
 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
   // for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( cell[2*m] );
        double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< predicted  << "&" << "pcell "<< cell[2*m] << "&" ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
	// when i was a young boy my dad son when you grow you will be the saviour of the damned the black parade
// sometimes i get the feeling that she's watching over me, we'll carry on, we'll carry on..			
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 10 << "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

for( int e = 0; e < d[m].size() ; e++ ){  //d[m].size = seqSitesm1[m].size  // since seqSitesm1 is parsed into  all fragments with one missing site
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
	for( int i = 0; i < d[m][e].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < d[m][e][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( d[m][e][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( d[m][e][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( d[m][e][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

} // for e
fo << "\\\\&&\\\\&&\\\\" << endl;

}// for m
fo.close();


}
		


void ExprPredictor::printFile3( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) 
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "format.tex",ios::app );
int j=0;
//cout << "cout seqSites.size " << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;


/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 double bottomofdborder = .5;
 double topofdborder = .8;
 int nrow = factorExprData.nRows();    
 int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
//cout << " getrowi " << endl;				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					//cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests " << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << gsl_vector_get(transition_Indicesbs,2) << endl; 



// for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

//		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
//		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);

////////////////////////////////////////
ofstream hes( "hes.txt" );
	vector< double > pars;
        par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators ); 
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       
	//Hassan start:
	pars.clear();
	pars = free_pars;
	double step = .01;
	bool a=1;
	//Hassan end
ofstream occm("occmat.txt");
cout << AllData.size() << endl;
cout << "AllData[0].size()"<< AllData[0].size() << endl;
cout << "AllData[1].size()"<< AllData[1].size() << endl;
cout << "AllBorders.size() " << AllBorders.size() << endl;
//cout << "Allb[1].size()"<< AllBorders[1].size() << endl;
//cout << "Allb[1].size()"<< AllBorders<< endl;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
JM = gsl_matrix_alloc( AllData[0].size()* AllData.size(),pars.size() );  // number of sequences(*number of bords) by the number of parameters
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence ( for given bord, loop over all seqs)
// Jacobian						// there are 12 bords, unless seqSitesmbot and the other snail vector are added.
//cout <<"bords "  << bords <<endl;
//for(int i=0; i < seqSites.size(); i++ ) {
for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
	//hes<< "Jacobian of " << ExprPredictor::seqNmes[i]  <<endl ;
	for (int jj = 0; jj < pars.size() ; jj++ ) {
	//hes << " pars jj " << jj << endl;
			//ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
		
			vector< double > parsjuh =pars;
			vector< double > parsjdh= pars;
			
			////cout << " add step " << endl;
			//cout << "pars(0) + step " << pars[0] << endl;
			//cout << "pars 0 + step " << pars[0] + step << endl;
			parsjuh[jj] = pars[jj] + step;
			parsjdh[jj] = pars[jj] - step;
			
			
 		       // vector< double > concs = factorExprData.getCol( cell[2*i] );
			//vector< double > concs = factorExprData.getCol( AllBorders[bords][i] );
		//cout << "about to get border " << endl;
//cout << " factorExprData.getCol( AllBorders[i][bords] ) "  << factorExprData.getCol( AllBorders[i][bords] )<< endl;
			vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
	// k fixed j up		
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}
	 ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators , a);
	///// k fixed j down
	 all_pars.clear();
	 free_par_counter = 0;
	 fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}


	 ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );


	///// calculate fxx or fxy
//cout.setf( ios::fixed );
 //   cout.precision( 3 ); 
//     cout.width( 8 );
 vector< double > fOcc(0);
hes.setf( ios::fixed);
hes.width(3);
	count = bords* nSeqs();
			ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parjdh );
 		        double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( par_model );
			if( jj==0){
			 func->predictOcc( AllData[bords][i], 5, concs,  fOcc );
 ; 
			//double p=func->predictExpr(AllData[bords][i], 5, concs ); // AllData[bords][i] = seqSites[i]
			occm << ExprPredictor::seqNmes[i]<< '\t'<<AllBorders[i][bords] <<'\t' << fOcc << endl;
			} //if i ==0

 		        double p = func->predictExpr(AllData[bords][i], 5, concs ); // AllData[bords][i] = seqSites[i]
			//hes  <<  setprecision(2) << ( pjuh- pjdh )/(2*step) << '\t' ;
			Jacobian[jj] += ( pjuh- pjdh )/(2*step);
			gsl_matrix_set( JM,count+ i, jj, ( pjuh- pjdh )/(2*step) );  // i is seqs, jj is parameters
			//count++;  // using counter to indicate position in bords*seq
			parsjuh.clear();
			parsjdh.clear();
//cout << " nottttttt "  << endl;
	}  // for jj
	//hes << endl ; // jacobian for each sequence
//}// for i in jacobian
//hes << endl << endl;
//cout <<" Jacobian " << endl << endl;
// Hessian

//for(int i=0; i < seqSites.size(); i++ ) {
//for(int i=0; i < AllData[bords].size(); i++ ) {
	//hes<< "Hessian of " << ExprPredictor::seqNmes[i] <<bords <<endl ;
	for (int j = 0; j < pars.size() ; j++ ) {
		for (int k = 0; k < pars.size() ; k++ ) {
		//for (int k = 0; k <=j ; k++ ) {
			//ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
		//cout << "hes " << j << " " << k <<endl;
			vector< double > parsjuh =pars;
			vector< double > parsjdh= pars;
			vector< double > parskuh =pars;  // this and the next aren't being used
			vector< double > parskdh= pars;
			vector< double > parsjkdh= pars;
			vector< double > parsjkuh =pars;
			vector< double > parsjukdh= pars;
			vector< double > parsjdkuh =pars;
			parsjuh[j] = pars[j] + step;    // f_j    first partial
			parsjdh[j] = pars[j] - step;    // f_j
			parskuh[k] = pars[k] + step;    // f_k
			parskdh[k] = pars[k] - step;     // f_k
			parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
			parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(x+h,y+h)
			parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
			parsjkdh[k] = pars[k] - step;
			parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
			parsjukdh[k] = pars[k] - step;
			parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-d )
			parsjdkuh[k] = pars[k] + step;
 		      

//vector< double > concs = factorExprData.getCol( cell[2*i] );
				vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
	// k fixed j up		
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}
	 ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
	///// k fixed j down
	 all_pars.clear();
	 free_par_counter = 0;
	 fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}


	 ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );

//  j fixed k up
		all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parskuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
// j fixed k down
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parskdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );


	///// calculate fxx or fxy
//cout.setf( ios::fixed );
 //   cout.precision( 3 ); 
//     cout.width( 8 ); 
//hes.setf( ios::fixed);
//hes.width(3);
	if( j == k ) {
			ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parjdh );
 		        double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( par_model );
 		        double p = func->predictExpr(AllData[bords][i], 5, concs );
			//hes << ( pjuh -2*p + pjdh )/(step*step) << '\t' ;
			Hessian[j][k] += ( pjuh -2*p + pjdh )/(step*step) ;
			parsjuh.clear();
			parsjdh.clear();
			parskuh.clear();
			parskdh.clear();
			parsjkuh.clear();
			parsjkdh.clear();
			parsjdkuh.clear();
			parsjukdh.clear();
	}
	if ( j !=k ) {


		
//////////// both u or d
			all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parjkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
		///// step by 2
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parjkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
////////////  mixed
			all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parjdkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
		///// step by 2
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parjukdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );


		ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parjdh );
 		        double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( par_model );
 		        double p = func->predictExpr(AllData[bords][i], 5, concs );
			  func = this->createExprFunc( parkuh );
			double pkuh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parkdh );
 		        double pkdh = func->predictExpr(AllData[bords][i], 5, concs );
			func = this->createExprFunc( parjkuh );
			double pjkuh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parjkdh );
 		        double pjkdh = func->predictExpr(AllData[bords][i], 5, concs );
			func = this->createExprFunc( parjukdh );
			double pjukdh = func->predictExpr(AllData[bords][i], 5, concs );
			 func = this->createExprFunc( parjdkuh );
 		        double pjdkuh = func->predictExpr(AllData[bords][i], 5, concs );
			
		//	hes << ( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step) << '\t' ;
			Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
			parsjuh.clear();
			parsjdh.clear();
			parskuh.clear();
			parskdh.clear();
			parsjkuh.clear();
			parsjkdh.clear();




	}  // if j!= k 



		} // for k
	//	hes << endl;	
	} // for j
//hes<< endl  <<endl ;
} // for i
} // for bords
hes << "Jacobian full" << endl;
for(int i = 0; i < pars.size() ; i++ ) {
hes <<  setprecision(2)  << Jacobian[i] << '\t' ;

}
hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
	for(int k = 0; k < pars.size() ; k++ ) {
	hes <<  setprecision(2)  <<  Hessian[j][k] << '\t' ;
	}
hes << endl;
}
occm.close();
/////////////////////////////////
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );

const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
 	Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
        gsl_matrix_set_row(IHessian , j , Hvector);     
}

int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);

hes << " inverse Hes " << endl;

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
        gsl_matrix_get_row( Hvector ,inverse , j );
	for(int i = 0; i < pars.size() ; i++ ) {
	hes << gsl_vector_get(Hvector,i)  << '\t' ;
	}
hes << endl;
}

double epsrel = 0 ; //.01;
gsl_matrix *covar;
covar = gsl_matrix_alloc( pars.size(), pars.size() );
gsl_multifit_covar(JM, epsrel, covar);
hes << "covar = JtJ inverse; error = sqrt(JtJ^-1) ,pg 446 gsl man" << endl;   /// why does this say inverse?
for( int i =0 ; i < covar->size1; i++) {
	for( int j =0 ; j < covar->size2; j++) {
	hes << sqrt(gsl_matrix_get(covar,i,j)) << '\t';	
	}
hes << endl;
}
hes << endl;
gsl_matrix * JMt = gsl_matrix_alloc( JM->size2, JM->size1);
gsl_matrix_transpose_memcpy(JMt,JM);   //(dest,source)
 gsl_matrix * CV = gsl_matrix_alloc( JM->size2, JM->size2);

////////

int oint2 =1;
int* spoint2 =&oint2;
//gsl_vector_mul(JMt,JM);

gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
1.0, JMt, JM,
0.0, CV);    // dgemm simply multiples two matrices:  JMt*JM,.. namely J transpose time J, and stores in CV
gsl_linalg_LU_decomp (CV,  p, spoint2);
gsl_matrix * inverseCV=gsl_matrix_alloc( JM->size2, JM->size2);
gsl_linalg_LU_invert (CV, p, inverseCV);

hes << " inverse JtJ = Covar, manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(JM->size2 ); 
       gsl_matrix_get_row( Hvector ,inverseCV , j );
	for(int i = 0; i < JM->size2 ; i++ ) {
	hes << gsl_vector_get(Hvector,i)  << '\t' ;
	}
hes << endl;
//Jvector = vector2gsl(Jacobian); 
//gsl_vector_mul(Jvectort,Jvector);
//double iip = 1/ip;
}
hes << endl;
hes << "  RMSE*(JtJ)^(-1) , parameter uncertainties , manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(JM->size2 ); 
       gsl_matrix_get_row( Hvector ,inverseCV , j );
	for(int i = 0; i < JM->size2 ; i++ ) {
hes << jRMSE.getObj()* gsl_vector_get(Hvector,i) << '\t' ;
}
hes << endl;
}
hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;

vector<int> indices(0);
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence 		( for given bord, loop over all seqs)						// there are 12 bords, unless seqSitesmbot and 		the other 	snail vector are added.
	for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
	indices.push_back(i);
	}
}
//hes << "JM matrix, something like a design matrix " << endl ;

//for( int i =0 ; i < JM->size1; i++) {
/*
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence 		( for given bord, loop over all seqs)						// there are 12 bords, unless seqSitesmbot and 		the other 	snail vector are added.
	for(int i=0; i < AllData[bords].size(); i++ ) {

	hes << ExprPredictor::seqNmes[i] <<" border "<<AllBorders[i][bords] << endl ;
	for( int j =0 ; j < JM->size2; j++) {
	hes << gsl_matrix_get(JM,i*bords,j) << '\t';
	//hes << gsl_matrix_get(covar,i,j) << '\t';	
	}
hes << endl;
}
}
*/

hes.close();
cout << " Alldata.size " << AllData.size() << " AllBorders.size() "<<AllBorders.size()<< endl;
cout << " Alldata[0].size " << AllData[0].size() << " AllBorders[0].size() "<<AllBorders[0].size()<< endl;
cout << " Alldata[1].size " << AllData[1].size() << " AllBorders[1].size() "<<AllBorders[1].size()<< endl;
cout <<  " AllBorders[0] " << endl<<AllBorders[0]<< endl;
cout << "Alldata mat is mxn matrix while AllBorders mat is nxm mat." << endl;

/*
int fixedsize = fix_pars.size();  
int rows = Occupancy -> size1;
int columns = Occupancy -> size2;
//int rows = Occw -> size1;
//int columns = Occw -> size2;
int i,j,k;
k=0;
gsl_matrix *X = gsl_matrix_alloc(columns,columns);
gsl_matrix *V = gsl_matrix_alloc(columns,columns);
gsl_vector *S =gsl_vector_alloc(columns);
gsl_vector *xx =gsl_vector_alloc(columns);  // x is in the row space, therefore it is a vector in Rcolumns, with at most row dimensions.
gsl_vector *b =gsl_vector_alloc(rows);     // b is in the column space of A.
gsl_vector_set_all( b,0 );
gsl_vector *work=gsl_vector_alloc(columns);
int rsvd;		// A must have more rows than columns or the same num
int rsvds;               //(A,V,S,work), on output A is replaced by U  (A=USVt)
rsvd=gsl_linalg_SV_decomp(Occupancy,V,S,work);
//rsvds=gsl_linalg_SV_solve(Occupancy,V,S,b,xx);
//gsl_matrix_transpose(V);
//printf ("x = \n");
//gsl_vector_fprintf (stdout, xx, "%g");

for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,j) == 0) {
		 fix_pars.clear(); 
		for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			fix_pars.push_back(gsl_matrix_get(V,j,i));
		}
		break;
	
	}

}

//fix_pars.clear();
//for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,2) < .1 ) {
	//	 fix_pars.clear(); 
	//	for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			//fix_pars.push_back(gsl_matrix_get(V,2,j));  // this is causing some parameters to be out of range!!!
			cout << "gsl_matrix_get(V,2,j)   = " << gsl_matrix_get(V,2,j)  << endl;
			//par_model.txpEffects[j] = gsl_matrix_get(V,2,j);  // these parameters may be out of range and will fuck things up !!
		}
	//else break;
	
//}
for(int i =0; i < columns;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			cout << "gsl_vector_get(S,i)  , singular values  = " << gsl_vector_get(S,i)  << endl;
		}
gsl_matrix_free( Occupancy );
*/
///////////////////////////////////////////////////////////////////////
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {//cout << " j = " << j << endl;
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	
	foo << endl;
	foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
		//if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
			
			else {//cout << " j = " << j << endl;
				

				if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();



/////////////////
 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
   // for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( cell[2*m] );
        double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << setprecision(2) << predicted  << "&" << "cell "<< cell[2*m] << "&" ;
//cout << tsites << endl;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
	// when i was a young boy my dad son when you grow you will be the saviour of the damned the black parade
// sometimes i get the feeling that she's watching over me, we'll carry on, we'll carry on..			
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {  // tsitesbot was changed from tsites below
			if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
	for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
	for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}


fo << "\\\\&&\\\\&&\\\\" << endl;

}// for m
fo.close();


}
		



void ExprPredictor::printFile3b( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) 
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "format.tex",ios::app );
int j=0;
//cout << "cout seqSites.size " << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;


/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 double bottomofdborder = .5;
 double topofdborder = .8;
 int nrow = factorExprData.nRows();    
 int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
//cout << " getrowi " << endl;				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					//cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests " << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << gsl_vector_get(transition_Indicesbs,2) << endl; 



// for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

//		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
//		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);

////////////////////////////////////////
ofstream hes( "hes.txt" );
	vector< double > pars;
        par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators ); 
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       
	//Hassan start:
	pars.clear();
	pars = free_pars;
	double step = .01;
	bool a=1;
	//Hassan end
ofstream occm("occmat.txt");
cout << AllData.size() << endl;
cout << "AllData[0].size()"<< AllData[0].size() << endl;
cout << "AllData[1].size()"<< AllData[1].size() << endl;
cout << "AllBorders.size() " << AllBorders.size() << endl;
//cout << "Allb[1].size()"<< AllBorders[1].size() << endl;
//cout << "Allb[1].size()"<< AllBorders<< endl;

/*
int fixedsize = fix_pars.size();  
int rows = Occupancy -> size1;
int columns = Occupancy -> size2;
//int rows = Occw -> size1;
//int columns = Occw -> size2;
int i,j,k;
k=0;
gsl_matrix *X = gsl_matrix_alloc(columns,columns);
gsl_matrix *V = gsl_matrix_alloc(columns,columns);
gsl_vector *S =gsl_vector_alloc(columns);
gsl_vector *xx =gsl_vector_alloc(columns);  // x is in the row space, therefore it is a vector in Rcolumns, with at most row dimensions.
gsl_vector *b =gsl_vector_alloc(rows);     // b is in the column space of A.
gsl_vector_set_all( b,0 );
gsl_vector *work=gsl_vector_alloc(columns);
int rsvd;		// A must have more rows than columns or the same num
int rsvds;               //(A,V,S,work), on output A is replaced by U  (A=USVt)
rsvd=gsl_linalg_SV_decomp(Occupancy,V,S,work);
//rsvds=gsl_linalg_SV_solve(Occupancy,V,S,b,xx);
//gsl_matrix_transpose(V);
//printf ("x = \n");
//gsl_vector_fprintf (stdout, xx, "%g");

for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,j) == 0) {
		 fix_pars.clear(); 
		for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			fix_pars.push_back(gsl_matrix_get(V,j,i));
		}
		break;
	
	}

}

//fix_pars.clear();
//for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,2) < .1 ) {
	//	 fix_pars.clear(); 
	//	for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			//fix_pars.push_back(gsl_matrix_get(V,2,j));  // this is causing some parameters to be out of range!!!
			cout << "gsl_matrix_get(V,2,j)   = " << gsl_matrix_get(V,2,j)  << endl;
			//par_model.txpEffects[j] = gsl_matrix_get(V,2,j);  // these parameters may be out of range and will fuck things up !!
		}
	//else break;
	
//}
for(int i =0; i < columns;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			cout << "gsl_vector_get(S,i)  , singular values  = " << gsl_vector_get(S,i)  << endl;
		}
gsl_matrix_free( Occupancy );
*/
///////////////////////////////////////////////////////////////////////
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 
	


/////////////////
 j=0;
int jjj=10000;
	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }



		for( int ii = 0; ii < seqSites[m].size(); ii++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		
			if ( i== seqSites[m][ii].start && seqSites[m][ii].factorIdx != 2) {
				j=0;

		  if ( abs(i-jjj) < 6 ) {
		  	while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
					fo << "\\color{green}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;
				}	// i++; j++;}   
					jjj=i; break;
				
			}//if abs
		  else{
				if( seqSites[m][ii].factorIdx == 0 ) { 

					while( j < 9 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
					fo << "\\color{blue}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;}   
					jjj=i; break;
				}//if tsitesii
				if( seqSites[m][ii].factorIdx == 1 ) { 

				while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
					fo << "\\color{green}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++; jjj=i;}   

				break;
				 }  // although j > tsites[i++], hence they'll appear right after one another.
				if(seqSites[m][ii].factorIdx == 2 ) {
/*
				while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
					// fo << "\\color{red}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}";  i++; j++; jjj=i;}
					fo << "\\color{black}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}";  i++; j++; jjj=i;}
	  			break; 
*/				}

			}// else abs()
	// when i was a young boy my dad son when you grow you will be the saviour of the damned the black parade
// sometimes i get the feeling that she's watching over me, we'll carry on, we'll carry on..			
				
			 }// if
		}// for ii
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }

	}//for i


j=0;
fo << "\\\\&&\\\\";

j=0;
/*
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
*/
j=0;
/*
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
*/


j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSites[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSites[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSites[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSites[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSites[m][i].factorIdx == 2 ) { fo << "\\color{black}{e}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;



fo  << endl;

}// for m
fo.close();


}
		


void ExprPredictor::printFile5( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) 
{
 ofstream hes( "hes5.txt" );
	vector< double > pars;
        par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators ); 
	int pars_size = pars.size();
hes << " pars " << endl << pars << endl;
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       
	//Hassan start:
	pars.clear();
	pars = free_pars;
	double step = .0000000001;
	bool a = 1;
	//spaceweights=0;
	//Hassan end
//ofstream occm("occmat.txt");
//cout << AllData.size() << endl;
//cout << "AllData[1].size()"<< AllData[1].size() << endl;
//cout << AllBorders.size() << endl;
//cout << "Allb[1].size()"<< AllBorders[1].size() << endl;
//cout << "Allb[1].size()"<< AllBorders<< endl;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
///////////////////////////////////////////////////////
JM = gsl_matrix_alloc( nSeqs(),pars.size() );  // number of sequences(*number of bords) by the number of 
///////////////////

	for (int j = 0; j < pars.size() ; j++ ) {
		for (int k = 0; k < pars.size() ; k++ ) {
		
			vector< double > parsjuh =pars;
			vector< double > parsjdh= pars;
			vector< double > parskuh =pars;  // this and the next aren't being used
			vector< double > parskdh= pars;
			vector< double > parsjkdh= pars;
			vector< double > parsjkuh =pars;
			vector< double > parsjukdh= pars;
			vector< double > parsjdkuh =pars;
			parsjuh[j] = pars[j] + step;    // f_j    first partial
			parsjdh[j] = pars[j] - step;    // f_j
			parskuh[k] = pars[k] + step;    // f_k
			parskdh[k] = pars[k] - step;     // f_k
			parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
			parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(k+h,j+h)
			parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
			parsjkdh[k] = pars[k] - step;
			parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
			parsjukdh[k] = pars[k] - step;
			parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-h )
			parsjdkuh[k] = pars[k] + step;
 		      

	// k fixed j up		
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}
	 ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
	//printPar( parjuh );
	///// k fixed j down
	 all_pars.clear();
	 free_par_counter = 0;
	 fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}


	 ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );

//  j fixed k up
		all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parskuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
// j fixed k down
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parskdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );


	///// calculate fxx or fxy
//cout.setf( ios::fixed );
 //   cout.precision( 3 ); 
//     cout.width( 8 ); 
//hes.setf( ios::fixed);
//hes.width(3);
	if( j == k ) {
			//ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh =  this->compRMSE2(parjuh); //func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parjdh );
 		        double pjdh = this->compRMSE2(parjuh); // func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( par_model );
 		        double p =  this->compRMSE2(par_model);//func->predictExpr3(seqSitesm1[i], 5, concs );
			//hes << ( pjuh -2*p + pjdh )/(step*step) << '\t' ;
			Hessian[j][k] += ( pjuh -2*p + pjdh )/(step*step) ;
			parsjuh.clear();
			parsjdh.clear();
			parskuh.clear();
			parskdh.clear();
			parsjkuh.clear();
			parsjkdh.clear();
			parsjdkuh.clear();
			parsjukdh.clear();
	}
	if ( j !=k ) {


		
//////////// both u or d
			all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parjkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
		///// step by 2
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parjkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
////////////  mixed
			all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}
		 ExprPar parjdkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
		///// step by 2
		 all_pars.clear();
		 free_par_counter = 0;
		 fix_par_counter = 0;
		for( int index = 0; index < indicator_bool.size(); index ++ ){
			if(  indicator_bool[ index ]  ){
				all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
			}
			else{
				all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
			}
		}


		ExprPar parjukdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);


		//ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh =this->compRMSE2(parjuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parjdh );
 		        double pjdh = this->compRMSE2(parjdh);//func->predictExpr3(seqSitesm1[i], 5, concs );
			 //func = this->createExprFunc( par_model );
 		        double p =this->compRMSE2(par_model);// func->predictExpr3(seqSitesm1[i], 5, concs );
			//  func = this->createExprFunc( parkuh );
			double pkuh = this->compRMSE2(parkuh);//func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parkdh );
 		        double pkdh =this->compRMSE2(parkdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
		//	func = this->createExprFunc( parjkuh );
			double pjkuh =this->compRMSE2(parjkuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parjkdh );
 		        double pjkdh =this->compRMSE2(parjkdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
			//func = this->createExprFunc( parjukdh );
			double pjukdh =this->compRMSE2(parjukdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parjdkuh );
 		        double pjdkuh =this->compRMSE2(parjdkuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
			printPar(parjkuh );
			printPar(parjkdh );
			printPar(parjdkuh );
			printPar(parjkdh );
			hes << pjkuh << " " << pjukdh << " " << pjdkuh << " " << pjkdh << endl;
			
			
			hes << ( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step) << '\t' << "delta " << 4*step*step << endl ;
			Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
			parsjuh.clear();
			parsjdh.clear();
			parskuh.clear();
			parskdh.clear();
			parsjkuh.clear();
			parsjkdh.clear();




	}  // if j!= k 



		} // for k
		//hes << endl;	
	} // for j
//hes<< endl  <<endl ;
//} // for i

hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
	for(int k = 0; k < pars.size() ; k++ ) {
	hes <<  setprecision(6)  <<  Hessian[j][k] << '\t' ;
	}
hes << endl;
}
//occm.close();
/////////////////////////////////
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );

const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
 	Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
//cout <<" gslvector before scale " << gsl_vector_get(Hvector,j) << endl ;
	gsl_vector_scale (Hvector, .5);
//cout <<" gslvector after scale " << gsl_vector_get(Hvector,j) << endl ;
        gsl_matrix_set_row(IHessian , j , Hvector);     
}

int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);

hes << " inverse .5*Hes, hence var of paramters " << endl;

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
        gsl_matrix_get_row( Hvector ,inverse , j );
	for(int i = 0; i < pars.size() ; i++ ) {
	hes << gsl_vector_get(Hvector,i)  << '\t' ;
	}
hes << endl;
}

hes << " standard error of paramters (sqrt((.5*hess)^-1)) " << endl;

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
        gsl_matrix_get_row( Hvector ,inverse , j );
	for(int i = 0; i < pars.size() ; i++ ) {
	if(i==j ){
	hes <<" stderr par "<<i << " = " << sqrt(gsl_vector_get(Hvector,i)) << '\t';	
	}
	//hes << sqrt(gsl_vector_get(Hvector,i))  << '\t' ;
	}
hes << endl;
}



hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;



gsl_matrix * Hessian2 = gsl_matrix_alloc(pars.size(),pars.size() );

const size_t N2 = pars.size();
gsl_permutation * p2 = gsl_permutation_alloc (N2);

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector2 =gsl_vector_alloc(pars.size() ); 
 	Hvector2 = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized

        gsl_matrix_set_row(Hessian2 , j , Hvector2);     
}
//gsl_matrix_view m
//= gsl_matrix_view_array (Hessian2, pars.size(),pars.size());
gsl_vector *eval = gsl_vector_alloc (pars.size());
gsl_matrix *evec = gsl_matrix_alloc (pars.size(), pars.size());
gsl_eigen_symmv_workspace * w =
gsl_eigen_symmv_alloc (pars.size());
//gsl_eigen_symmv (&m.matrix, eval, evec, w);
gsl_eigen_symmv (Hessian2, eval, evec, w);
gsl_eigen_symmv_free (w);
gsl_eigen_symmv_sort (eval, evec,
GSL_EIGEN_SORT_ABS_ASC);
{
int i;
for (i = 0; i < pars.size(); i++)
{
double eval_i
= gsl_vector_get (eval, i);
gsl_vector_view evec_i
= gsl_matrix_column (evec, i);
printf ("eigenvalue = %g\n", eval_i);
printf ("eigenvector = \n");
gsl_vector_fprintf (stdout,
&evec_i.vector, "%g");
}
}
gsl_vector_free (eval);
gsl_matrix_free (evec);



hes.close();


}
void ExprPredictor::printFile4( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) 
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   os.close();
ofstream fo( "format.tex",ios::app );
int j=0;
//cout << "cout seqSites.size " << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;


/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 double bottomofdborder = .5;
 double topofdborder = .8;
 int nrow = factorExprData.nRows();    
 int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
//cout << " getrowi " << endl;				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					//cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests " << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << gsl_vector_get(transition_Indicesbs,2) << endl; 



// for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

//		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
//		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);

////////////////////////////////////////
ofstream hes( "hes4.txt" );
	vector< double > pars;
        par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators ); 
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       
	//Hassan start:
	pars.clear();
	pars = free_pars;
	double step = .001;
	bool a = 1;
	//spaceweights=0;
	//Hassan end
//ofstream occm("occmat.txt");
//cout << AllData.size() << endl;
//cout << "AllData[1].size()"<< AllData[1].size() << endl;
//cout << AllBorders.size() << endl;
//cout << "Allb[1].size()"<< AllBorders[1].size() << endl;
//cout << "Allb[1].size()"<< AllBorders<< endl;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
///////////////////////////////////////////////////////
JM = gsl_matrix_alloc( nSeqs(),pars.size() );  // number of sequences(*number of bords) by the number of 
///////////////////
for(int i=0; i < nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
	//hes<< "Jacobian of " << ExprPredictor::seqNmes[i]  <<endl ;
	for (int jj = 0; jj < pars.size() ; jj++ ) {
	
//hes<< endl  <<endl ;
} // for i
} // for bords
hes << "Jacobian full" << endl;
for(int i = 0; i < pars.size() ; i++ ) {
hes <<  setprecision(2)  << Jacobian[i] << '\t' ;

}
hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
	for(int k = 0; k < pars.size() ; k++ ) {
	hes <<  setprecision(6)  <<  Hessian[j][k] << '\t' ;
	}
hes << endl;
}
//occm.close();
/////////////////////////////////
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );

const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
 	Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
cout <<" gslvector before scale " << gsl_vector_get(Hvector,j) << endl ;
	gsl_vector_scale (Hvector, .5);
cout <<" gslvector after scale " << gsl_vector_get(Hvector,j) << endl ;
        gsl_matrix_set_row(IHessian , j , Hvector);     
}

int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);

hes << " inverse .5*Hes, hence var of paramters " << endl;

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
        gsl_matrix_get_row( Hvector ,inverse , j );
	for(int i = 0; i < pars.size() ; i++ ) {
	hes << gsl_vector_get(Hvector,i)  << '\t' ;
	}
hes << endl;
}

hes << " standard error of paramters (sqrt((.5*hess)^-1)) " << endl;

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() ); 
        gsl_matrix_get_row( Hvector ,inverse , j );
	for(int i = 0; i < pars.size() ; i++ ) {
	if(i==j ){
	hes <<" stderr par "<<i << " = " << sqrt(gsl_vector_get(Hvector,i)) << '\t';	
	}
	//hes << sqrt(gsl_vector_get(Hvector,i))  << '\t' ;
	}
hes << endl;
}




double epsrel =.001;
gsl_matrix *covar;
covar = gsl_matrix_alloc( pars.size(), pars.size() );
gsl_multifit_covar(JM, epsrel, covar);
hes << "covar = JtJ inverse, sqrt((JtJ)^-1) " << endl;   /// why does this say inverse?
for( int i =0 ; i < covar->size1; i++) {
	for( int j =0 ; j < covar->size2; j++) {
	if(i==j ){
	hes <<" stderr par "<<i << " = " << sqrt(gsl_matrix_get(covar,i,j)) << '\t';	
	}
	}
hes << endl;
}
hes << endl;
gsl_matrix * JMt = gsl_matrix_alloc( JM->size2, JM->size1);
gsl_matrix_transpose_memcpy(JMt,JM);   //(dest,source)
 gsl_matrix * CV = gsl_matrix_alloc( JM->size2, JM->size2);

////////

int oint2 =1;
int* spoint2 =&oint2;
//gsl_vector_mul(JMt,JM);

gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
1.0, JMt, JM,
0.0, CV);    // dgemm simply multiples two matrices:  JMt*JM,.. namely J transpose time J, and stores in CV
gsl_linalg_LU_decomp (CV,  p, spoint2);
gsl_matrix * inverseCV=gsl_matrix_alloc( JM->size2, JM->size2);
gsl_linalg_LU_invert (CV, p, inverseCV);
hes << " Jacobian " << endl;
	for(int j = 0; j < JM->size1 ; j++ ) {
		gsl_vector * Hvector =gsl_vector_alloc(JM->size2 ); 
	       gsl_matrix_get_row( Hvector ,JM , j );
		for(int i = 0; i < JM->size2 ; i++ ) {
		hes << gsl_vector_get(Hvector,i)  << '\t' ;
		}
	hes << endl;
	
	}
	hes << endl;
	hes << " inverse JtJ = Covar, manual " << endl;
	for(int j = 0; j < JM->size2 ; j++ ) {
		gsl_vector * Hvector =gsl_vector_alloc(JM->size2 ); 
	       gsl_matrix_get_row( Hvector ,inverseCV , j );
		for(int i = 0; i < JM->size2 ; i++ ) {
		hes << gsl_vector_get(Hvector,i)  << '\t' ;
		}
	hes << endl;
	
	}
	hes << endl;
	hes << "  MSE*(JtJ)^(-1) , parameter uncertainties , manual " << endl;
	hes << "  MSE=sum_i ( y_i - ym(x_i) )^2 / (n-p) " << endl;
	hes << "  where n= nSeqs()*nConds(), and p is the number of fit parameters " << endl;
	hes << "  it's possible n should not include all nConds(), since they are/seem dependent " << endl;
	for(int j = 0; j < JM->size2 ; j++ ) {
		gsl_vector * Hvector =gsl_vector_alloc(JM->size2 ); 
	       gsl_matrix_get_row( Hvector ,inverseCV , j );
		for(int i = 0; i < JM->size2 ; i++ ) {
			hes << jRMSE.getObj() /( nSeqs()*nConds() )* gsl_vector_get(Hvector,i) << '\t' ;
		}
		hes << endl;
	}
hes  << " chi^2jRMSE.getObj() " << jRMSE.getObj() << endl;
hes  << " hes << nSeqs()*nConds()" << nSeqs()*nConds() << endl;
hes  << " hes << chi^2/nSeqs()*nConds()" << jRMSE.getObj()/(nSeqs()*nConds()) << endl;



hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;



gsl_matrix * Hessian2 = gsl_matrix_alloc(pars.size(),pars.size() );

const size_t N2 = pars.size();
gsl_permutation * p2 = gsl_permutation_alloc (N2);

for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector2 =gsl_vector_alloc(pars.size() ); 
 	Hvector2 = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized

        gsl_matrix_set_row(Hessian2 , j , Hvector2);     
}
//gsl_matrix_view m
//= gsl_matrix_view_array (Hessian2, pars.size(),pars.size());
gsl_vector *eval = gsl_vector_alloc (pars.size());
gsl_matrix *evec = gsl_matrix_alloc (pars.size(), pars.size());
gsl_eigen_symmv_workspace * w =
gsl_eigen_symmv_alloc (pars.size());
//gsl_eigen_symmv (&m.matrix, eval, evec, w);
gsl_eigen_symmv (Hessian2, eval, evec, w);
gsl_eigen_symmv_free (w);
gsl_eigen_symmv_sort (eval, evec,
GSL_EIGEN_SORT_ABS_ASC);
{
int i;
for (i = 0; i < pars.size(); i++)
{
double eval_i
= gsl_vector_get (eval, i);
gsl_vector_view evec_i
= gsl_matrix_column (evec, i);
printf ("eigenvalue = %g\n", eval_i);
printf ("eigenvector = \n");
gsl_vector_fprintf (stdout,
&evec_i.vector, "%g");
}
}
gsl_vector_free (eval);
gsl_matrix_free (evec);



hes.close();

}


void ExprPredictor::printFile2( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
//os.width( 18 ); 
    
    // print binding weights
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
			//count = count + 1 ;//ouut << endl;
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "format.tex",ios::app );
int j=0;
//cout << "cout seqSites.size " << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;
/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 double bottomofdborder = .5;
 double topofdborder = .8;
 int nrow = factorExprData.nRows();    
 int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
//cout << " getrowi " << endl;				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					//cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests " << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << gsl_vector_get(transition_Indicesbs,2) << endl; 



// for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

//		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
//		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
//cout << "eh " << endl;
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 //cout << tsites << endl;
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {//cout << " j = " << j << endl;
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	
	foo << endl;
	foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
		//if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
			
			else {//cout << " j = " << j << endl;
				

				if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();


//cout << "heh" << endl;
/////////////////
 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
//cout << "heh" << endl;
//ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
   // for ( int j = 0; j < nConds(); j++ ) {
     //   vector< double > concs = factorExprData.getCol( cell[2*m] );
       // double predicted = func->predictExpr(tsites, 5, concs );
//cout << jRMSE.cell<< endl;
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
   // for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( cell[2*m] );
        double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << predicted  << "&" << "pcell "<< jRMSE.cell[2*m] << "&" ;
//fo << ExprPredictor::seqNmes[m]<< " "  << "&" << "cell "<< cell[2*m] << "&" ;  //<< setprecision(2) << predicted // cut out 111611
//cout << "heh" << endl;
//cout << cell[0] << endl;
//cout << "heh" << endl;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		//if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
	
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
//cout << " phe" << endl;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {  // tsitesbot was changed from tsites below
			if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
//cout << seqSitesf2 << endl;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
	for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
	for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
	
	//	if (tsitesbot.size() == 1 ) { break; }
		for (;;) {
			if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {//cout << " j = " << j << endl;
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
				if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
				if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }// else
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}


fo << "\\\\&&\\\\&&\\\\" << endl;

}// for m
fo.close();


}

int ExprPredictor::train()
{	
    // random number generator
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
    
    // training using the default initial values with random starts
    ExprPar par_default( nFactors() );
    train( par_default, rng ); 
    
    return 0;	
}

int ExprPredictor::predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs ) 
{
//printPar( par_model );
    targetExprs.clear();
    
    // create site representation of the target sequence
//     SiteVec targetSites;
//     SeqAnnotator ann( motifs, energyThrs );	
//     ann.annot( targetSeq, targetSites );
     //  ExprPar pp;     
    // predict the expression
    ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr3( targetSites, targetSeqLength, concs );
        targetExprs.push_back( predicted );
    }
    
    return 0; 
}



ModelType ExprPredictor::modelOption = LOGISTIC;
int ExprPredictor::estBindingOption = 1;    // 1. estimate binding parameters; 0. not estimate binding parameters
ObjType ExprPredictor::objOption = SSE;

double ExprPredictor::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
{
    vector< int > shifts; 
    for ( int s = -maxShift; s <= maxShift; s++ ) {
        shifts.push_back( s ); 
    }

    vector< double > cov; 
    vector< double > corr; 
    cross_corr( x, y, shifts, cov, corr ); 
    double result = 0, weightSum = 0; 
//     result = corr[maxShift]; 
    result = *max_element( corr.begin(), corr.end() );
//     for ( int i = 0; i < shifts.size(); i++ ) {
//         double weight = pow( shiftPenalty, abs( shifts[i] ) ); 
//         weightSum += weight; 
//         result += weight * corr[i]; 
//     }
//     result /= weightSum; 

    return result; 
}

int ExprPredictor::maxShift = 5; 
double ExprPredictor::shiftPenalty = 0.8; 
int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 0;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
int ExprPredictor::nSimplexIters = 3;
int ExprPredictor::nGradientIters = 6;

int ExprPredictor::randSamplePar( const gsl_rng* rng, ExprPar& par ) const
{
    int counter= 0 ;
    // sample binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
	    if( indicator_bool[counter] ) {
            double rand_weight = exp( gsl_ran_flat( rng, log( ExprPar::min_weight ), log( ExprPar::max_weight ) ) ); 
            par.maxBindingWts[i] = rand_weight;
	    counter++;
	    }
	    else{ counter++; } 
        }        
    }

///////////////////////////////// ////clean this up 525, or wherever the logistic was..
if(modelOption == BINS){
for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
		
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
			 // double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );  //524 restet interaction, was 1, same for factorIntMat
		if( indicator_bool[counter] ) {
		double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interaction ), log( ExprPar::max_interaction ) ) );
                par.theV[ i][ j][k] = rand_interaction;	   
		counter++;
		}       //for    
		else{counter++;}
                 } //if indicator bool
		
            }  //if construct
           
        }
    }
   
for ( int i = 0; i < nFactors(); i++ ) {
       // for ( int j = 0; j <= i; j++ ) { 81611
	for ( int j = 0; j < i; j++ ) {
            if (repIndicators[i] ) {
		
		for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
		if( indicator_bool[counter] ) {
			 // double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );  //524 restet interaction, was 1, same for factorIntMat
double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interactionr ), log( ExprPar::max_interactionr ) ) );
                par.theVr[ i][ j][k] = rand_interaction;
		counter++;
		}// for k
		else{counter++;}
		}// if indicator
			              
                
            }  //if construct
            
        }
    }
    
}
  
////////////////////////////
    // sample the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       if( indicator_bool[counter] ) {
            double rand_effect =  gsl_ran_flat( rng,  ExprPar::min_effect_Thermo , ExprPar::max_effect_Thermo  );
            par.txpEffects[i] = rand_effect;
	    counter++;
	    }
	else{counter++;}
       
    }
//cout << "counter = " << counter << '\t' << " indicator size " << indicator_bool.size() << endl;
   /*
    
    // sample the basal transcription
    double rand_basal;
    if ( modelOption == LOGISTIC ) 
        rand_basal = gsl_ran_flat( rng, ExprPar::min_basal_Logistic, ExprPar::max_basal_Logistic );
    else
        rand_basal = exp( gsl_ran_flat( rng, log( ExprPar::min_basal_Thermo ), log( ExprPar::max_basal_Thermo ) ) );
    par.basalTxp = rand_basal;
    */
    return 0;
}

bool ExprPredictor::testPar( const ExprPar& par ) const
{
    // test binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
                return false; 
        }        
    }

 // adjust theV matrix
   if( modelOption == BINS) {
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           if(coopMat(i,j)){
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( par.theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) { 
			return false; }
             }
           }
        }
     }
     for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           if(repIndicators[i]){
             for(int k =0; k < ExprPar::nbins; k++){
                  if ( par.theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) || par.theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) { 
			return false; }
             }
           }
        }
     }
    }
 // test the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
       
            if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        
    }

    return true;    
}
void ExprPredictor::printPar2(  ) 
{ 
cout <<"dltwrep" << endl;
}
/*
void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 3 ); 
//     cout.width( 8 ); 
    
    // print binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
		if(i == 0){
            		cout << " maxBindDl" << "\t"<< par.maxBindingWts[i] << "\t";
		}
		if(i==1){
			cout << " maxBindTw" << "\t"<< par.maxBindingWts[i] << "\t";
		}
		if(i==2){
			cout << " maxBindSn" << "\t"<< par.maxBindingWts[i] << "\t";
		} 
        }        
    }
    cout << endl;



//if(modelOption== BINS){
  for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if(coopMat(i,j))  { cout << endl;
			for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
				if( i == 1 && j ==0 ) { 
				cout <<"dltwact"<< k<<"\t" << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				if( i == 2 && j ==0 ) {  
				cout <<"dlsnact"<< k<<"\t" << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				
				if( i== 2 && j == 1) { 
				cout <<"twsnact"<< k<<"\t" << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
			}
			
  		}
  	}
  }
cout <<" whatttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt" <<endl;
for(int i =0; i < nFactors(); i++){
	for(int j =0; j <=i; j++){
		if(repIndicators[i])  { cout << endl;
			for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
				if( i == 1 && j ==0 ) { 
				cout <<"dltwrep"<< k<<"\t" << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				if( i == 2 && j == 0) {  
				cout <<"dlsnrep"<< k<<"\t" << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				
				if( i== 2 && j == 1) { 
				cout <<"twsnrep"<< k<<"\t" << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
			}
			
  		}
  	}
  }

//}
cout << endl;
   // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
      cout<<"txp"  << i<<"\t" << par.txpEffects[i] << "\t";
        
    }
cout << endl;

}
*/
void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 3 ); 
//     cout.width( 8 ); 
    
    // print binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
		if(i == 0){
            		cout << par.maxBindingWts[i] << "\t";
		}
		if(i==1){
			cout << par.maxBindingWts[i] << "\t";
		}
		if(i==2){
			cout << par.maxBindingWts[i] << "\t";
		} 
        }        
    }
   



//if(modelOption== BINS){
  for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if(coopMat(i,j))  { cout << endl;
			for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
				if( i == 1 && j ==0 ) { 
				cout << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				if( i == 2 && j ==0 ) {  
				cout  << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				
				if( i== 2 && j == 1) { 
				cout  << par.theV[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
			}
			
  		}
  	}
  }

for(int i =0; i < nFactors(); i++){
	for(int j =0; j <=i; j++){
		if(repIndicators[i])  { cout << endl;
			for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
				if( i == 1 && j ==0 ) { 
				cout << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				if( i == 2 && j == 0) {  
				cout << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
				
				if( i== 2 && j == 1) { 
				cout  << par.theVr[i][j][k]  <<"\t";}
				//count = count + 1 ;//ouut << endl;
			}
			
  		}
  	}
  }

//}

   // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
      cout << par.txpEffects[i] << "\t";
        
    }
cout << endl;

}

ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{	
    return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par,  binBord, binBordr );
// all the parameters are privates from ExprPredictor// except par, which is created inside gsl_obj_f
}
// motifs, intFunc.. are all private members of ExprPredictor... whereas par is passed as a reference parameter..
ExprFunc* ExprPredictor::createExprFunc2(  ExprPar& par )  // the heading of an accessor method is followed by const, to tell the compiler that such a method may not modify any of the instance variables.,  a declaration ,starting with const, (e.g. const ExprPar& Par) of a named variable or parameter (e.g. par) indicates that this variable can not be modified..  by declaring a parametere whose type is type& we have the type object refer directly to the object instead of making a copy of the object being passed,  this is called a reference parameter.. 
{	
    return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par, binBord , binBordr);
}// motifs, intFunc.. are all private members of ExprPredictor... whereas par is passed as a reference parameter..

double ExprPredictor::compRMSE( const ExprPar& par ) const
{

    // create the expression function
    ExprFunc* func = createExprFunc( par );
            double rss = 0;
    // error of each sequence
    double squaredErr = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }


        double beta;
        squaredErr += least_square( predictedExprs, observedExprs, beta );
    	
//for(int i = 0; i < 5; i++){
//cout << "Bincounts" <<endl; // i<< "   " << ExprFunc::Bc[i] << endl;
//}
      for ( int i = 0; i < nConds(); i++ ) {  // 10/04
        rss += ( predictedExprs[i] -  observedExprs[i] ) * ( predictedExprs[i] -  observedExprs[i] ); //beta*x is p, where x is now acting as a above, and beta is x above.
    }   
   }
double rmse =  rss/ ( nSeqs() * nConds() ) ; 
//cout << "RMSE " << rmse << endl;
//printPar( par );
 // double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
//cout << "RMSE " << rmse << endl;
    return rmse;
  //return rmse;
}
  double ExprPredictor::compf( const ExprPar& par ,  vector<double>&ff)  //gsl_vector * f ) 
{
    spaceweights=0;
    // create the expression function
    ExprFunc* func = createExprFunc( par );
	//printPar( par ) ;
            double rss = 0;
    // error of each sequence

    double genespaceweight=0; 
    for ( int i = 0; i < nSeqs(); i++ ) {
	//genespaceweight=genespaceweight+1;  // pseudocount for each gene
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
	vector< double > concsweight(0);
	for ( int j = 1; j < nConds()-1; j++ ) {
	 concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
	 genespaceweight=genespaceweight+concsweight[j];
	}

      for ( int j = 0; j < nConds(); j++ ) {  // 10/04
	
	// chi^2 = [ ( yi - ym)/sigma ] ^ 2    where sigma = sqrt( p(1-p) )
	double obser = observedExprs[j];
	double predic =predictedExprs[j];
	if(observedExprs[j]==0){ obser =.01;}
	//if(predictedExprs[j]==0){ predic =.001;}
	if(observedExprs[j]==1){ obser =.99;}
	//if(predictedExprs[j]==1){ predic =.999;}
	   // rss +=( predic -  obser) * ( predic -  obser )/ (obser*(1-obser )); 
	    rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic )); 
		//cout << " setting f " << " i " << i << " j " << endl;
		//gsl_vector_set( f,nConds()*i+j, predic-obser); 
		ff.push_back( (predic-obser)/ sqrt(obser));
      //  rss +=( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] )/ (observedExprs[j]*(1-observedExprs[j])+.1 ); 
	   //  cout <<" se/p(1-p)  "<< ( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] )/ (obser*(1-obser) ) << endl;
		//cout << " P[i](1-P[i]) " << (observedExprs[j]*(1-observedExprs[j]) ) << endl;
		//cout << " P(1-P) " << (obser*(1-obser) ) << endl;
		//cout << " obser " << obser << endl;
		//cout << " predic " << predic << endl;
	}
//for(int i = 0; i < 5; i++){
//cout << "Bincounts" <<endl; // i<< "   " << ExprFunc::Bc[i] << endl;
//}
      
   }

  for ( int i = 0; i < nSeqs(); i++ ) {
        for ( int j = 0; j < nConds(); j++ ) {
      
//cout <<" gslvector f "<<" i " <<i <<" j " << j << gsl_vector_get(f,nConds()*i+j) << endl ;
	}
    }
//cout <<" gslvector after f "  << endl ;  
	double rmse = rss;
    return rmse; // /(nSeqs()*nConds() );
  //return rmse;
}
double ExprPredictor::compRMSE2( const ExprPar& par ) 
{
    spaceweights=0;
    // create the expression function
    ExprFunc* func = createExprFunc( par );
	//printPar( par ) ;
            double rss = 0;
    // error of each sequence

    double genespaceweight=0; 
    for ( int i = 0; i < nSeqs(); i++ ) {
	//genespaceweight=genespaceweight+1;  // pseudocount for each gene
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
	vector< double > concsweight(0);
	for ( int j = 1; j < nConds()-1; j++ ) {
	 concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
	 genespaceweight=genespaceweight+concsweight[j];
	}

      for ( int j = 0; j < nConds(); j++ ) {  // 10/04
	
	// chi^2 = [ ( yi - ym)/sigma ] ^ 2    where sigma = sqrt( p(1-p) )
	double obser = observedExprs[j];
	double predic =predictedExprs[j];
	//if(observedExprs[j]==0){ obser =.001;}
	//if(predictedExprs[j]==0){ predic =.001;}
	//if(observedExprs[j]==1){ obser =.999;}
	//if(predictedExprs[j]==1){ predic =.999;}
	   // rss +=( predic -  obser) * ( predic -  obser )/ (obser*(1-obser )); 
	    rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic )); 
      //  rss +=( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] )/ (observedExprs[j]*(1-observedExprs[j])+.1 ); 
	   //  cout <<" se/p(1-p)  "<< ( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] )/ (obser*(1-obser) ) << endl;
		//cout << " P[i](1-P[i]) " << (observedExprs[j]*(1-observedExprs[j]) ) << endl;
		//cout << " P(1-P) " << (obser*(1-obser) ) << endl;
		//cout << " obser " << obser << endl;
		//cout << " predic " << predic << endl;
	}
//for(int i = 0; i < 5; i++){
//cout << "Bincounts" <<endl; // i<< "   " << ExprFunc::Bc[i] << endl;
//}
      
   }
spaceweights=genespaceweight;
double rmse =  rss;   //  / ( nSeqs() * nConds() )  ; 
//nSeqs() * nConds()  > genespaceweight, this is an attempt to correct for the correlated errors, by reducing the dof.

//cout << "RMSE " << rmse << endl;
//printPar( par );
 // double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
//cout << "RMSE " << rmse << endl;
    return rmse; // /(nSeqs()*nConds() );
  //return rmse;
}
/*
double ExprPredictor::compRMSE3( const ExprPar& par , int i) 
{

    // create the expression function
    ExprFunc* func = createExprFunc( par );
            double rss = 0;
//printPar( par ) ;
    // error of each sequence
    double squaredErr = 0;
    //for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
       }


    //    double beta;
    //    squaredErr += least_square( predictedExprs, observedExprs, beta );
    	
//for(int i = 0; i < 5; i++){
//cout << "Bincounts" <<endl; // i<< "   " << ExprFunc::Bc[i] << endl;
//}
vector< double > concsweight(0);
double genespaceweight=1; //1.0/nSeqs();  //setting to 1.0/nSeqs() is a pseudocount.
// j starts at 1 instead of j=0, due to exprData(i,j-1), and ends at nConds()-1, instead of nConds..
for ( int j = 1; j < nConds()-1; j++ ) {
  concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
}
 rss +=( predictedExprs[0] -  observedExprs[0] ) * ( predictedExprs[0] -  observedExprs[0] ); 
 rss +=( predictedExprs[nConds()] -  observedExprs[nConds()] ) * ( predictedExprs[nConds()] -  observedExprs[nConds()] ); 
      for ( int j = 1; j < nConds()-1; j++ ) {  // 10/04
	if(concsweight[j]< .05 ){ continue; }
        rss +=( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] ); 
	genespaceweight=genespaceweight+concsweight[j];
	
	}
     
// the chisquare degrees of freedom is spaceweights, which is the number of independent data points minus the number of fitted parameters.  The number of independent data points is determined by the number of independent errors (yi - ym(i)) where yi is observed expression at space point i, and ym is the model at point i.  For runs of errors that are consecutively negative or positive, we have correlated errors, which are not independent.  To determine runs of dependent errors we use the structure 'concsweight' which looks for differential expresssion.  Then based on the amount of differential expression we can estimate the amount of independence between consecutive data points.  Hence Concsweight has nothing to do with weighted least squares.  Concsweight is not the variance or std of the expression.  It's a structure designed to disentangle the dependency error structure along space.  For example a profile like 1111111111110000000000000 we want to reduce to just the two points 10 at the internal border.  This is because all the flanking expression in the profile has a heavy dependency on the adjacent flanking expression point in space.  These points should not be counted as a full degree of freedom (a data point).  Now, its possible that if there's equal amount of 1's as 0's the information needed for the two classes will be balanced (the fit will work fine), but we are still messing up the chi square statistic.  Of course, each gene has a dependency for the entire profile, this is partly based on the fact that seqsitesm1 was designed to cover repression based on the activation.  Regardless we will treat each gene independently, and for neuroectoderm genes repressed in mesoderm, they will be treated as if the mesoderm portion of the gene is a completely new gene (new degree of freedom).  While genes only differential expression is at the neuroectoderm dorsal ectoderm border, they'll be treated with less degrees of freedom (since they don't have the mesdoerm differential expression). 

// another way to do this is construct the Covariance data matrix using the expression profile.  Then the chisquare is not sum_ij (yij - ymij) over space i and gene j, rather chi= sum_x xT* sigma^-1 x, where x is the expression profile for gene x.  That is we assume each gene's profile is drawn from a multivariate normal (the dependencies from internal spatial points will be seen by offdiagonal elements inside sigma. 

// in the calculation of rss, we do use concsweight as a weigth like the variance in weighted least square, probably need to change calculation of spacewieght and genespaceweight to be updated by concsweight instead of 1 and .5.
double rmse = rss/genespaceweight; ///rss/ (  nConds() )  ; 
//cout << "RMSE " << rmse << endl;
//printPar( par );
 // double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
//cout << "RMSE " << rmse << endl;
    return   rmse;
  //return rmse;
}
*/
double ExprPredictor::compRMSE3( const ExprPar& par , int i) 
{

    // create the expression function
    ExprFunc* func = createExprFunc( par );
            double rss = 0;
//printPar( par ) ;
    // error of each sequence
    double squaredErr = 0;
    //for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
       }


    //    double beta;
    //    squaredErr += least_square( predictedExprs, observedExprs, beta );
    	
//for(int i = 0; i < 5; i++){
//cout << "Bincounts" <<endl; // i<< "   " << ExprFunc::Bc[i] << endl;
//}
vector< double > concsweight(0);
double genespaceweight=1; //1.0/nSeqs();  //setting to 1.0/nSeqs() is a pseudocount.
// j starts at 1 instead of j=0, due to exprData(i,j-1), and ends at nConds()-1, instead of nConds..
for ( int j = 1; j < nConds()-1; j++ ) {
  concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
}
// rss +=( predictedExprs[0] -  observedExprs[0] ) * ( predictedExprs[0] -  observedExprs[0] ); 
// rss +=( predictedExprs[nConds()] -  observedExprs[nConds()] ) * ( predictedExprs[nConds()] -  observedExprs[nConds()] ); 
      for ( int j = 0; j < nConds(); j++ ) {  // 10/04
	//if(concsweight[j]< .05 ){ continue; }
       // rss +=( predictedExprs[j] -  observedExprs[j] ) * ( predictedExprs[j] -  observedExprs[j] )/ (observedExprs[j]*(1-observedExprs[j]) ); 
	double obser = observedExprs[j];
	double predic =predictedExprs[j];
	//if(observedExprs[j]==0){ obser =.001;}
	//if(predictedExprs[j]==0){ predic =.001;}
	//if(observedExprs[j]==1){ obser =.999;}
	//if(predictedExprs[j]==1){ predic =.999;}
	//rss +=( predic -  obser) * ( predic -  obser )/ (obser*(1-obser )); 
	rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic )); 
	// cross entropy
	//rss += -observedExprs[j]*log( predictedExprs[j] ) -( observedExprs[j]*(1-observedExprs[j]) )*log( 1 - predictedExprs[j] ) ;

	//genespaceweight=genespaceweight+concsweight[j];
	
	}
     
// the chisquare degrees of freedom is spaceweights, which is the number of independent data points minus the number of fitted parameters.  The number of independent data points is determined by the number of independent errors (yi - ym(i)) where yi is observed expression at space point i, and ym is the model at point i.  For runs of errors that are consecutively negative or positive, we have correlated errors, which are not independent.  To determine runs of dependent errors we use the structure 'concsweight' which looks for differential expresssion.  Then based on the amount of differential expression we can estimate the amount of independence between consecutive data points.  Hence Concsweight has nothing to do with weighted least squares.  Concsweight is not the variance or std of the expression.  It's a structure designed to disentangle the dependency error structure along space.  For example a profile like 1111111111110000000000000 we want to reduce to just the two points 10 at the internal border.  This is because all the flanking expression in the profile has a heavy dependency on the adjacent flanking expression point in space.  These points should not be counted as a full degree of freedom (a data point).  Now, its possible that if there's equal amount of 1's as 0's the information needed for the two classes will be balanced (the fit will work fine), but we are still messing up the chi square statistic.  Of course, each gene has a dependency for the entire profile, this is partly based on the fact that seqsitesm1 was designed to cover repression based on the activation.  Regardless we will treat each gene independently, and for neuroectoderm genes repressed in mesoderm, they will be treated as if the mesoderm portion of the gene is a completely new gene (new degree of freedom).  While genes only differential expression is at the neuroectoderm dorsal ectoderm border, they'll be treated with less degrees of freedom (since they don't have the mesdoerm differential expression). 

// another way to do this is construct the Covariance data matrix using the expression profile.  Then the chisquare is not sum_ij (yij - ymij) over space i and gene j, rather chi= sum_x xT* sigma^-1 x, where x is the expression profile for gene x.  That is we assume each gene's profile is drawn from a multivariate normal (the dependencies from internal spatial points will be seen by offdiagonal elements inside sigma. 

// in the calculation of rss, we do use concsweight as a weigth like the variance in weighted least square, probably need to change calculation of spacewieght and genespaceweight to be updated by concsweight instead of 1 and .5.
double rmse = rss; ///rss/ (  nConds() )  ; 
//cout << "RMSE " << rmse << endl;
//printPar( par );
 // double rmse = sqrt( squaredErr / ( nSeqs() * nConds() ) ); 
//cout << "RMSE " << rmse << endl;
    return   rmse;//  /( nSeqs() * nConds() ) ;
  //return rmse;
}
double ExprPredictor::compRMSE4( const ExprPar& par , int i, int j) 
{

    // create the expression function
    ExprFunc* func = createExprFunc( par );
            double rss = 0;
    // error of each sequence
    double squaredErr = 0;
    //for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
      //  for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            //cout<< " predict "  << predicted << endl;
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
double obser=observed;
if(observed==0){ obser =.01;}
	//if(predicted==0){ predic =.001;}
	if(observed==1){ obser =.99;}
	//if(predictedExprs[j]==1){ predic =.999;}
	   // rss +=( predic -  obser) * ( predic -  obser )/ (obser*(1-obser )); 


    return    predicted/sqrt(obser) ; /// (obser*(1-obser ));//rmse;
  //Jacobian is w*(partial f / partial x)  where x is parameter and f is t value, w is weight (inverse variance)
}
////////////////////////////////////////////////////////////////////////
//compavgcorr takes the binding data of mulitiple infactors, and associates each infactor with its binding data, through correlation.,/
//all the correlations between infactor occupancy and binding data are stored in  vector corrs.
//the sum of all these correlations is used in the end.
////////////////////////////////////////////////////////////////////////

double ExprPredictor::compAvgCorr( const ExprPar& par ) const
{
	ExprFunc* func = createExprFunc( par );
	vector< double > fOcc;
	vector< double > corrs; 	
	for ( int i = 0; i < seqSitesb.size(); i++ ) {   //seqs is a vector<vector<Sequence>> a Sequence has a DNA sequence and DNA datalike length
		vector< double > predicted; 
		for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {			
			// predicted binding of the i-th factor to the j-th sequence in the i-th experiment // 914 jacobc i'th experiment is the same as i'th factor
			vector< double > concs = factorExprData.getCol( 0 );  //103

			func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
 			predicted.push_back( fOcc[i] );// fOcc[i] represents the Experiment number, so if data for two experiments i will be 0 and 1.
// hence we will have to put the sequences in the sequence data twice, once for twist and once for dorsal.
// notice how this chimes with corrs.push_back( correlation(predicted, bindingData[i ])) below, hence bindingData[i] will have both dl and twist.			
		}
		corrs.push_back( correlation( predicted, bindingData[ i ] ) );
	}	
	
	return mean( corrs ); 
}
////////////////////////////////////////////////
//works a little different than compAvgCorr, in that there is no reason to fill the predicted vector.
//other than that, the methods are the same: for each infactor the occupancy is calculated and placed in occstemp.
//occstemp is then copied into occs, i think vector's copy constructor should have no problem with the copying.
/////////////////////////////////////////////////////////////

void ExprPredictor::predict2(vector< vector< double > >& occs ) 
{
	ExprFunc* func = createExprFunc( par_model );
	vector< double > fOcc;
	vector< double > corrs;
	vector< vector< double > > occstemp(  seqSitesb.size(),vector< double >(seqSitesb[0].size()) ); // [0] is irrelevant since all experiments have same # of seqs	
	for ( int i = 0; i < seqSitesb.size(); i++ ) {   // number of experiments
		
		for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {		// number of sequences	
			vector< double > concs = factorExprData.getCol( 0 );  //103

			func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
 //predicted.push_back( fOcc[i] );// fOcc[i] represents the Experiment number, so if data for two experiments i will be 0 and 1.
// hence we will have to put the sequences in the sequence data twice, once for twist and once for dorsal.
// notice how this chimes with corrs.push_back( correlation(predicted, bindingData[i ])) below, hence bindingData[i] will have both dl and twist.
			//occs[i].push_back( fOcc[i] );
			occstemp[i][j] =  fOcc[i] ;
					
		}
	}	
	
occs = occstemp;
}
//////////////////////////////////////////////////////////////
// compavgcorr2 is used for the expression data, and is called in gsl_obj_f in an additive way.
// compavgcorr2 simply takes correlation between each pattern and the predicted pattern,
// then taking the average over all these, in the same additive way as above,hence each sequence is  equally weighted
// although one may say that there are hidden weights within the patterns, as some patterns may be more epensive
// due to irreproducibility, such as multimodals with sharp modes.
/////////////////////////////////////////////////////////


double ExprPredictor::compAvgCorr2(  ExprPar& par ) 
{
//cout << " exxxxxpre comavgcoorr2  " << endl;
    // create the expression function
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
////////////////////////////////////////////////////////
    // vector< double > pars;
//pars.clear();
//cout << pars.size() <<  "if no num to left pars is empty:" <<endl;  
   // par_model.adjust();
  //  par.getFreePars2( pars, coopMat, actIndicators, repIndicators ); 
Matrix f = factorExprData;
Matrix e = exprData;
/*
for ( int i = 0; i < nSeqs(); i++ ) {
           anny.annoty2( seqsy[ i], seqSites[ i ], f, e, *func , j);

        }
*/
//cout <<" seqsites "<< seqSites[0] << '\t';

//cout << " noottt hereeeee " << endl;
//////////////////////////////////////////////////////
    // Pearson correlation of each sequence
    double totalSim = 0;
    for ( int i = 0; i < 3; i++ ) {  // should be i < nSeqs()
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
          //  for ( int i = 0; i < nSeqs(); i++ ) {
if(j == 0 ) {
           anny.annoty2( seqsy[ i], seqSites[ i ], f, e, *func , j, ExprPredictor::seqNmes[i]);
}
else {anny.annoty2( seqsy[ i], seqSites[ i ], f, e, *func , j, " "); }
       // }
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
         //   cout <<" predicted= "<< predicted << '\t';
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
//cout << " size predicted ad observed " << predictedExprs.size() << "  " << observedExprs.size() << endl;
corrs.push_back( correlation( predictedExprs, observedExprs ) );
//        totalSim += corr( predictedExprs, observedExprs ); 
//         cout << "Sequence " << i << "\t" << corr( predictedExprs, observedExprs ) << endl;
    }	
par_model = par;
cout <<  endl;
if( corrs != corrs) {
cout << "corrs != corrs" << endl;
return .5;
}
//cout << " mean(corrs)  " << mean(corrs) << endl;
return mean( corrs );
  //  return totalSim / nSeqs();
}


double ExprPredictor::compAvgCorrborder8(  ExprPar& par ) 
{
double expmin = .15;
double expmax = .65;
double bottomofdborder = .15;
double topofdborder = .65;
/////////////////////////////////////////////////////
///////////////////////////////////////////////////
///// this may effect optimizer
vector< int > initint;
//cout << " size of iniint " << initint.size() << endl;
for ( int ab = 0; ab < nSeqs() ; ab++) {
AllBorders.push_back( initint );
}
par_model = par;
AllBorders.clear();
AllData.clear();
//cout << " about to create " << endl;

///////////////////////////////////////////////////////
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
///////////////////////**************************8///////////////////////
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	// create rng type
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
	//int label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
   // int count =0;
    while ( ss >> label ) {
	
	//for( int j=0; j< 10 ; j++ ) {
	
	//count++; 
	//ostringstream ssout;
	//ssout << count;
	//colLabels.push_back( ssout.str() );
	//}
	colLabels.push_back( label );
    }
    // read the data
//for(int i; i < 5; i++){
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
		//double cellrep = atof(val);    // the number of val that occur is the same as the num of label.
		//for( int j=0; j< 10 ; j++ ) {
		 	
/*
			rand_site_index = gsl_ran_flat( rng, -.15, .15 );
			if( val + rand_site_index > 1 ) {  // change val to cell_label ?
			 vals.push_back( 1 );
			}
	
			if( val + rand_site_index < 0 ) {
			 vals.push_back( 0 );
			}
	
			if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
			vals.push_back( val + rand_site_index );
			}
*/

	vals.push_back( val  );		
	//} // for j
			

	}// while ss >> val
	data.push_back( vals );
    }

//cout << " collabels size " << colLabels.size() << endl;
//cout << colLabels << endl;
//cout << " data size " << data.size() << endl;
//cout << "datacol size0 " << data[0].size() << endl;
/*
cout << "datacol size 1 " << data[1].size() << endl;
cout << data[1] << endl;
cout << "datacol size20 " << data[20].size() << endl;
cout << data[20] << endl;
*/
Matrix data1( data );
//data1.save( "matrixexpr.txt" );
//////////////////////////////////////////******************************888888/////////////////////
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
cout << " exprData nrows " << nrow << endl;
gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tit;
int tib;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
	tit = 0;
	tib =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
		gsl_vector_set(transition_Indicesb, i , tib);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tit =j;
					//cout << " tit " << tit << endl; 
					gsl_vector_set(transition_Indicest, i , tit); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tit+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tit + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tib=tit + counter;
									//cout << "tib  " << tib << endl;
									gsl_vector_set(transition_Indicesb, i , tib);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tit == 0) { 
				  tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tib = minindex;
				gsl_vector_set(transition_Indicest, i , tit); 
				gsl_vector_set(transition_Indicesb, i , tib); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else     
//cout << endl << endl;
}// for i

 vector< double > predictedExprs;
        vector< double > observedExprs;
    // Pearson correlation of each sequence
    double totalSim = 0;
double rms = 0;
//cout << "nrows " << nrow << endl;
cout << "transitionindicest "<< endl;
for (int i=0; i<nrow;i++) {
cout << endl << gsl_vector_get(transition_Indicest,i) << endl; 

}

cout << "transitionindicesb " << endl;
for (int i=0; i<nrow;i++) {
cout << endl << gsl_vector_get(transition_Indicesb,i) << endl; 

}

//////////////////////////////////////////////////////////////////////////////////////////////////
/* dorsal border of target genes (controlled by sequence and dorsal twist concentrations ) remeber transition indexes start at 0, spreadsheet starts at 1 */
//////////////////////////////////////////////////////////////////////////////////////////////////

int c = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
		c++;
		cell.push_back( gsl_vector_get(transition_Indicest,i) );
		cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         
 
            		double observed = .7;  // .7 was for debugging, needs to change
           	AllBorders[i].push_back( 0 );
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         	AllBorders[i].push_back( 0 );
           
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
	//cout << " before anny " << endl;
		   vector< double > concst = factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSites[ i ], f, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
		AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
// cout << " names ep seqnmes " << ExprPredictor::seqNmes[i] << endl;
        	    double predictedt = func->predictExpr( seqSites[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );

           // was 4
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );

		    vector< double > concs2 = factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbot[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
// was 2
           		AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
        	    double predicted2 = func->predictExpr( seqSitesbot[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
		    
       		}
	}  // for i
//////////////////////////////////////////
/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 bottomofdborder = .5;
 topofdborder = .8;
 nrow = factorExprData.nRows();    
 ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	cout << " for row " << i << " max index " << mi << endl;
	cout << " for row " << i << " max index " << mi << endl;
cout << "rowexpdata 0 p0int "<< endl << gsl_vector_get(rowexprData,0) << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					cout << " tit " << tits << " mi " << mi << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

cout << "transitionindicests "<< endl << gsl_vector_get(transition_Indicests,2) << endl; 
cout << "transitionindicesbs " << endl <<gsl_vector_get(transition_Indicesbs,2) << endl; 

//cout << "transitionindicests 0 p0int "<< endl << gsl_vector_get(transition_Indicests,2) rowexprData<< endl
 for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
		vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
		if( concheck[2]!=0){ seqSitesm1[ i ]=seqSites[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );


		double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
       //  AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );

 continue; } // don't add more snail sites, if i is not an nee.
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
        	    double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}


for ( int i = 0; i < nSeqs(); i++ ) {
//vector< vector< Site > > se;
//cout << " annoty3 "<< endl;
		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicesbs,2) );
        	 //  anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , 10, ExprPredictor::seqNmes[i],i, seqSites[i], se);
		vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
		if( concheck[2]!=0){seqSitesm2[ i ]=seqSitesbot[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) ); 

		 double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
        // AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );


continue; } // don't add more snail sites, if i is not an nee.
		   anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
           
        	    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}
/*  // uncomment if seqSitesm1bot and seqSitesm2bot are added to ExprPredictor structure, (these are the sites from neuro or meso at bottom of pattern..

 for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
		    anny.annoty3( seqsy[ i], seqSitesm1bot[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSitesbot[i]);
        	    double predictedt = func->predictExpr( seqSitesm1bot[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - 1 );
}
for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
		    anny.annoty3( seqsy[ i], seqSitesm2bot[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSitesbot[i]);
        	    double predictedt = func->predictExpr( seqSitesm2bot[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - 1 );
}
*/
//cout << " predicted Exprs " << predictedExprs << endl;
//cout << " observed Exprs " << observedExprs << endl;
//cout << " cell " <<  cell << endl;

//double cor = correlation( predictedExprs, observedExprs );

//cout << "cor = " << cor << endl;

/*
repeat whole procedure for shifted profiles of dorsal and twist up 2 cells down 2cells..  (actually multiply by 1.5 and .5 and add to oriigianal..
*/
Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
temp[i]= temp[i] + noise;
//temp.push_back(2);
if(temp[i] > 1 ) { temp[i] = 1 ;}
}
f2.setRow(0,temp);

//cout << "temp "<< temp << endl;
//cout << "f2 "<< f2 << endl;
 for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
               	
		//cell.push_back( gsl_vector_get(transition_Indicest,i) );
		//cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
     //    cout << " gsl_transindex = 0 in f2 " << endl;
           AllBorders[i].push_back( 0 );
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         AllBorders[i].push_back( 0 );
           
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
//cout << "inside annoty " << endl;
		//cout << "entering annoty42 " << endl;	

		   vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
//cout << "inside annoty42 " << endl;
           AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
        	    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );

		    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
// was 2
           		AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
        	    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
		    
       		}
	}  // for i
//////////////////////////////////////////
/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 bottomofdborder = .5;
 topofdborder = .8;
 nrow = factorExprData.nRows();    
 ncol =factorExprData.nCols();
//gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
//gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
//int tits;
//int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
//cout << " getrowi " << endl;				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
				//	cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests " << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << gsl_vector_get(transition_Indicesbs,2) << endl; 



 for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

		 vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicests,2) );  // here we retrieve the snail border..
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
			vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
		if( concheck[2]!=0){
			 seqSitesm1f2[ i ]=seqSitesf2[i]; AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
			 double predictedt = func->predictExpr( seqSitesm1f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
	        // AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
		continue; } // don't add more snail sites, if i is not an nee.
		    anny.annoty3( seqsy[ i], seqSitesm1f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
        	    double predictedt = func->predictExpr( seqSitesm1f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}


for ( int i = 0; i < nSeqs(); i++ ) {
//vector< vector< Site > > se;
		 vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicesbs,2) );
        	 //  anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , 10, ExprPredictor::seqNmes[i],i, seqSites[i], se);
		vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
	      if( concheck[2]!=0){ seqSitesm2f2[ i ]=seqSitesbotf2[i]; AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
		double predictedt = func->predictExpr( seqSitesm2f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
continue; } // don't add more snail sites, if i is not an nee.
		   anny.annoty3( seqsy[ i], seqSitesm2f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
           	AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
        	    double predictedt = func->predictExpr( seqSitesm2f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}
/*
 vector< SiteVec >& seqSites;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbot;
vector< SiteVec >& seqSitesm1;
vector< SiteVec >& seqSitesm2;

//////////////dorsal + noise
 vector< SiteVec >& seqSitesf2;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbotf2;
vector< SiteVec >& seqSitesm1f2;
vector< SiteVec >& seqSitesm2f2;
///////////////dorsal - noise
 vector< SiteVec >& seqSitesf3;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbotf3;
vector< SiteVec >& seqSitesm1f3;
vector< SiteVec >& seqSitesm2f3;
*/
//vector< vector< SiteVec > >AllData;
AllData.push_back(seqSites ) ;
AllData.push_back(seqSitesbot ) ;
AllData.push_back( seqSitesm1) ;
AllData.push_back(seqSitesm2 ) ;
AllData.push_back(seqSitesf2 ) ;
AllData.push_back(seqSitesbotf2 ) ;
AllData.push_back(seqSitesm1f2 ) ;
AllData.push_back(seqSitesm2f2 ) ;
/*
AllData.push_back( seqSitesf3 ) ;  // currently the f3 are not being used, what are they?
AllData.push_back( seqSitesbotf3);
AllData.push_back( seqSitesm1f3);
AllData.push_back(seqSitesm2f3 );
*/
//cout << " correlation( predictedExprs, observedExprs ) = " << correlation( predictedExprs, observedExprs ) << endl;
//return correlation( predictedExprs, observedExprs );
obj_model = sqrt(rms/(nSeqs()*8));  // this was set for testing of correlation matrix of parameters, 10.5.11
return  obj_model;  // rms;

}

double ExprPredictor::compAvgCorrborder(  ExprPar& par ) 
{
double expmin = .35;
double expmax = .65;
double bottomofdborder = .35;
double topofdborder = .65;
/////////////////////////////////////////////////////
///////////////////////////////////////////////////
///// this may effect optimizer
vector< int > initint;
//cout << " size of iniint " << initint.size() << endl;
for ( int ab = 0; ab < nSeqs() ; ab++) {
AllBorders.push_back( initint );
}
par_model = par;
AllBorders.clear();
AllData.clear();
//cout << " about to create " << endl;

///////////////////////////////////////////////////////
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
///////////////////////**************************8///////////////////////
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	// create rng type
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
	//int label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
   // int count =0;
    while ( ss >> label ) {
	
	//for( int j=0; j< 10 ; j++ ) {
	
	//count++; 
	//ostringstream ssout;
	//ssout << count;
	//colLabels.push_back( ssout.str() );
	//}
	colLabels.push_back( label );
    }
    // read the data
//for(int i; i < 5; i++){
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
		//double cellrep = atof(val);    // the number of val that occur is the same as the num of label.
		//for( int j=0; j< 10 ; j++ ) {
		 	
/*
			rand_site_index = gsl_ran_flat( rng, -.15, .15 );
			if( val + rand_site_index > 1 ) {  // change val to cell_label ?
			 vals.push_back( 1 );
			}
	
			if( val + rand_site_index < 0 ) {
			 vals.push_back( 0 );
			}
	
			if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
			vals.push_back( val + rand_site_index );
			}
*/

	vals.push_back( val  );		
	//} // for j
			

	}// while ss >> val
	data.push_back( vals );
    }

//cout << " collabels size " << colLabels.size() << endl;
//cout << colLabels << endl;
//cout << " data size " << data.size() << endl;
//cout << "datacol size0 " << data[0].size() << endl;
/*
cout << "datacol size 1 " << data[1].size() << endl;
cout << data[1] << endl;
cout << "datacol size20 " << data[20].size() << endl;
cout << data[20] << endl;
*/
Matrix data1( data );
//data1.save( "matrixexpr.txt" );
//////////////////////////////////////////******************************888888/////////////////////
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tit;
int tib;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
	tit = 0;
	tib =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
		gsl_vector_set(transition_Indicesb, i , tib);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tit =j;
					//cout << " tit " << tit << endl; 
					gsl_vector_set(transition_Indicest, i , tit); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tit+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tit + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tib=tit + counter;
									//cout << "tib  " << tib << endl;
									gsl_vector_set(transition_Indicesb, i , tib);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tit == 0) { 
				  tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tib = minindex;
				gsl_vector_set(transition_Indicest, i , tit); 
				gsl_vector_set(transition_Indicesb, i , tib); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else     
//cout << endl << endl;
}// for i

 vector< double > predictedExprs;
        vector< double > observedExprs;
    // Pearson correlation of each sequence
    double totalSim = 0;
double rms = 0;
//cout << "nrows " << nrow << endl;
//cout << "transitionindicest "<< endl;
for (int i=0; i<nrow;i++) {
//cout << endl << gsl_vector_get(transition_Indicest,i) << endl; 

}

//cout << "transitionindicesb " << endl;
for (int i=0; i<nrow;i++) {
//cout << endl << gsl_vector_get(transition_Indicesb,i) << endl; 

}

//////////////////////////////////////////////////////////////////////////////////////////////////
/* dorsal border of target genes (controlled by sequence and dorsal twist concentrations ) remeber transition indexes start at 0, spreadsheet starts at 1 */
//////////////////////////////////////////////////////////////////////////////////////////////////

int c = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
		c++;
		cell.push_back( gsl_vector_get(transition_Indicest,i) );
		cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         
 
            		double observed = .7;  // .7 was for debugging, needs to change
           	AllBorders[i].push_back( 0 );
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         	AllBorders[i].push_back( 0 );
           
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
	//cout << " before anny " << endl;
		   vector< double > concst = factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSites[ i ], f, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
		AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
// cout << " names ep seqnmes " << ExprPredictor::seqNmes[i] << endl;
        	    double predictedt = func->predictExpr( seqSites[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );

           // was 4
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );

		    vector< double > concs2 = factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbot[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
// was 2
           		AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
        	    double predicted2 = func->predictExpr( seqSitesbot[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
		    
       		}
	}  // for i
//////////////////////////////////////////
/*
mesoderm/neuroectoderm border (snail border)
*/// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
//////////////////////////////////////////
 bottomofdborder = .5;
 topofdborder = .8;
 nrow = factorExprData.nRows();    
 ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
//  int i = 2;  // get snail border
	vector< double > reD;	
				 
	reD = factorExprData.getRow(i);	
	//cout << " number of rows " << nrow << endl;
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
//cout << " about to hit exptrace mes " << endl;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
		gsl_vector_set(transition_Indicesbs, i , tibs);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tits =j;
					//cout << " tit " << tits << endl; 
					gsl_vector_set(transition_Indicests, i , tits); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									//cout << "tib  " << tibs << endl;
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else    

//cout << endl << endl;
}//for i

//cout << "transitionindicests "<< endl << gsl_vector_get(transition_Indicests,2) << endl; 
//cout << "transitionindicesbs " << endl <<gsl_vector_get(transition_Indicesbs,2) << endl; 



 for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.

		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
		vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
		if( concheck[2]!=0){ seqSitesm1[ i ]=seqSites[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );


		double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
        // AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
cout << " rms " << rms << endl; 
 continue; } // don't add more snail sites, if i is not an nee.
        	  // anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i],seqSitesm1d1[i], ddd);
	//  anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , 0, ExprPredictor::seqNmes[i],i, seqSites[i], seqSitesm1d1);
          //d.push_back(seqSitesm1d1);
		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
        	    double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}


for ( int i = 0; i < nSeqs(); i++ ) {
//vector< vector< Site > > se;
//cout << " annoty3 "<< endl;
		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicesbs,2) );
        	 //  anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , 10, ExprPredictor::seqNmes[i],i, seqSites[i], se);
		vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
		//cout << ExprPredictor::seqNmes[i] << " " << concheck << endl;
		if( concheck[2]!=0){seqSitesm2[ i ]=seqSitesbot[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) ); 

		 double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
       //  AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );


continue; } // don't add more snail sites, if i is not an nee.
		   anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
           
        	    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += ( predictedt - observedt )*( predictedt - observedt );
}

AllData.push_back(seqSites ) ;
AllData.push_back(seqSitesbot ) ;
AllData.push_back( seqSitesm1) ;
AllData.push_back(seqSitesm2 ) ;

/*
AllData.push_back( seqSitesf3 ) ;  // currently the f3 are not being used, what are they?
AllData.push_back( seqSitesbotf3);
AllData.push_back( seqSitesm1f3);
AllData.push_back(seqSitesm2f3 );
*/
//cout << " correlation( predictedExprs, observedExprs ) = " << correlation( predictedExprs, observedExprs ) << endl;
//return correlation( predictedExprs, observedExprs );
obj_model = sqrt(rms/(nSeqs()*4));  // this was set for testing of correlation matrix of parameters, 10.5.11
return  obj_model;  // rms;

}
/*
the uncertainty in the params can be thought of as 3 data points for linear, one could define 3 lines, m1,m2,m3, and take std of these three slopes.. which would be different than the range that m could take on.
for a hyperplane of parameters one could do the same and define take the std of the different parameters, for all the different potential planes...
however for a nonlinear object one should think of a particular parameter, with respect to the rmse, then one would see some sort of landscape..  one could collect all the local minima below a certain tolerance say rmse < .01 or something,,  than each parameter for each minima represents a potential value of the parameter, however it doesn't make sense to say  the parameter could have a std here, because say one has four minima,, a1,a2,a3,a4 are the different parameter values at the minima, if one says the std is roughly max(a1,,,a4) - min(a1,,a4),, the problem is all the values <a> could take on in a+-std, are not necessarily in the minima since a1 through a4 may be all defined for distinct values of b1-b4, and c1-c4,... where b and c are other parameters at the minima,, that is saying that the min(a) = f(b,c,..z)  where b,c,..z are the other parameters. 
*/









/*
the uncertainty in the params can be thought of as 3 data points for linear, one could define 3 lines, m1,m2,m3, and take std of these three slopes.. which would be different than the range that m could take on.
for a hyperplane of parameters one could do the same and define take the std of the different parameters, for all the different potential planes...
however for a nonlinear object one should think of a particular parameter, with respect to the rmse, then one would see some sort of landscape..  one could collect all the local minima below a certain tolerance say rmse < .01 or something,,  than each parameter for each minima represents a potential value of the parameter, however it doesn't make sense to say  the parameter could have a std here, because say one has four minima,, a1,a2,a3,a4 are the different parameter values at the minima, if one says the std is roughly max(a1,,,a4) - min(a1,,a4),, the problem is all the values <a> could take on in a+-std, are not necessarily in the minima since a1 through a4 may be all defined for distinct values of b1-b4, and c1-c4,... where b and c are other parameters at the minima,, that is saying that the min(a) = f(b,c,..z)  where b,c,..z are the other parameters. 
*/


// compAvgCorrborder2 was created for testing of output "the.txt", I think it can now be deleted..
double ExprPredictor::compAvgCorrborder2(  ExprPar& par ) 
{
double expmin = .35;
double expmax = .65;
double bottomofdborder = .35;
double topofdborder = .65;
/////////////////////////////////////////////////////
///////////////////////////////////////////////////
///// this may effect optimizer


par_model = par;


///////////////////////////////////////////////////////
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
///////////////////////**************************8///////////////////////
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	// create rng type
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
	//int label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
   // int count =0;
    while ( ss >> label ) {
	
	//for( int j=0; j< 10 ; j++ ) {
	
	//count++; 
	//ostringstream ssout;
	//ssout << count;
	//colLabels.push_back( ssout.str() );
	//}
	colLabels.push_back( label );
    }
    // read the data
//for(int i; i < 5; i++){
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
		//double cellrep = atof(val);    // the number of val that occur is the same as the num of label.
		//for( int j=0; j< 10 ; j++ ) {
		 	
/*
			rand_site_index = gsl_ran_flat( rng, -.15, .15 );
			if( val + rand_site_index > 1 ) {  // change val to cell_label ?
			 vals.push_back( 1 );
			}
	
			if( val + rand_site_index < 0 ) {
			 vals.push_back( 0 );
			}
	
			if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
			vals.push_back( val + rand_site_index );
			}
*/

	vals.push_back( val  );		
	//} // for j
			

	}// while ss >> val
	data.push_back( vals );
    }
/*
cout << " collabels size " << colLabels.size() << endl;
cout << colLabels << endl;
cout << " data size " << data.size() << endl;
cout << "datacol size0 " << data[0].size() << endl;
cout << "datacol size 1 " << data[1].size() << endl;
cout << data[1] << endl;
cout << "datacol size20 " << data[20].size() << endl;
cout << data[20] << endl;
*/
//cout << "datacol size0 " << data[0].size() << endl;
Matrix data1( data );
//data1.save( "matrixexpr.txt" );
//////////////////////////////////////////******************************888888/////////////////////
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tit;
int tib;
///////////////////////////////////////////////////////////// this loop first finds top of border, based on max exression
///////////////////////////////////////////////////////////// then the loop continues until expr passes the bottome border
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices) 
	
	double exptrace;
	double exppeek;
///////////////  set both indices to zero
	tit = 0;
	tib =0;
	if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
		gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
		gsl_vector_set(transition_Indicesb, i , tib);
		//continue;
		break;
	}
	else{

///////////////////////////////////////////// this is the main of the loop first setting the topborder
	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
			// the failure of the above condition indicates the trace is diminishing in value, so push index just before the topofborder flag fails.
			if(exppeek < topofdborder ) { 
				        tit =j;
					//cout << " tit " << tit << endl; 
					gsl_vector_set(transition_Indicest, i , tit); 
					//if ( gsl_vector_get(rowexprData,tit) > bottomofdborder ) {  // this should automatically be true
					//set the tib:
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tit+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tit + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tib=tit + counter;
									//cout << "tib  " << tib << endl;
									gsl_vector_set(transition_Indicesb, i , tib);
									break;	
								}//else
					} // for(;;)
						//} // if bottomofdboder
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tit == 0) { 
				  tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
				  tib = minindex;
				gsl_vector_set(transition_Indicest, i , tit); 
				gsl_vector_set(transition_Indicesb, i , tib); 
			} 	   
		 	                   // pushback the ti that is the transition index..  
			break;
			} // if < topofdborder
			// we push BEFORE the failure flag because we are training with assumption of dorsal occupancy, for sharp borders, we want expression ON as the index.
			// this condition may set tit to 0 for traces that are all zero, hence the minindex condition below
			////// the above condition sets the topindex, so we know need to set the bottom index, (note we have exhauted all the possible conditions on exptrace for the top
		}//  else
	}// for j
}// else     
//cout << endl << endl;
}// for i

 vector< double > predictedExprs;
        vector< double > observedExprs;
    // Pearson correlation of each sequence
    double totalSim = 0;
double rms = 0;


//////////////////////////////////////////////////////////////////////////////////////////////////
/* dorsal border of target genes (controlled by sequence and dorsal twist concentrations ) remeber transition indexes start at 0, spreadsheet starts at 1 */
//////////////////////////////////////////////////////////////////////////////////////////////////

Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
temp[i]= temp[i] + noise;
//temp.push_back(2);
if(temp[i] > 1 ) { temp[i] = 1 ;}
}
//f2.setRow(0,temp);

//cout << "temp "<< temp << endl;
//cout << "f2 "<< f2 << endl;
 for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
               	
		//cell.push_back( gsl_vector_get(transition_Indicest,i) );
		//cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
        cout << " gsl_transindex = 0 in f2 " << endl;
           
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            		predictedExprs.push_back( predicted );
         
           cout << " gslindex = 0 " << endl;
            		double observed = .7;  // .7 was for debugging, needs to change
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
//cout << "inside annoty " << endl;
		//cout << "entering annoty42 " << endl;	

		   vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
//cout << "inside annoty42 " << endl;
           
        	    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );

		    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
// was 2
           
        	    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     //rms += abs( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
		    
       		}
	}  // for i

return rms/nSeqs();

}
int Random(int n)
{
	return rand() % n ;
}

int ExprPredictor::createccdata()
{
double expmin = .35;
double expmax = .65;

int k=3;
/////////////////////////////////////////////////////
///////////////////////////////////////////////////
///// this may effect optimizer



///////////////////////////////////////////////////////
    
///////////////////////**************************8///////////////////////
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	// create rng type
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    
    // read the data
//for(int i; i < 5; i++){
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
	 rand_site_index = gsl_ran_flat( rng, -.15, .15 );
	if( val + rand_site_index > 1 ) {
	 vals.push_back( 1 );
	}
	
	if( val + rand_site_index < 0 ) {
	 vals.push_back( 0 );
	}
	
	if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
	vals.push_back( val + rand_site_index );
	}
	} // while ss >> val
        data.push_back( vals );
    }
Matrix data1( data );
//data1.save( "matrixexpr.txt" );
//////////////////////////////////////////******************************888888/////////////////////
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
if( m < expmin ) {
	gsl_vector_set(transition_Indices, i , ti);  // here ti = NULL
	//continue;
	break;
	}
else{
	 
        for (int j=mi; j< ncol; j++)  {  
		/*
	        double temp;
		temp = m/2;     
		*/	        
		if (expmax < gsl_vector_get(rowexprData,j) ) {
		//	cout << " expmax < e " << j <<'\t' ;
			continue;
	        }
		
		else {

			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..  
}// else     
}// for i

for (int i=0; i<nrow;i++) {
cout << "trans index : " << gsl_vector_get(transition_Indices,i) << endl;
}
/////////////////////////start with original seq and expr file (declassified)
ifstream seq_file("efastanotw10txt");  //("efastanotw6.txt");
	ifstream expr_file( "expre10.tab");//"expre6.tab");
assert (seq_file.is_open() &&  expr_file.is_open());
string header;
	//read the first line (header) in expression file
	getline(expr_file, header);
	
	string temp;
	vector <string> seq_name;
	vector <string> seq;
	vector <string> expr;
	
	seq.clear();
	expr.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name1 (temp, 1, temp.length() - 1);
		seq_name.push_back(name1);
		
		getline(seq_file, temp);
		seq.push_back(temp);
				
		getline(expr_file, temp);
		expr.push_back(temp);
	}
	
	int size = seq.size();
//cout << " size of seq.size  " << size << endl;
//cout << seq[0] << endl;
/////////////////////////create meso and neur files based on cell attribute of exprpredictor class (loop over nseqs, check cell, if cell < 8, append seq and exp to meso. etc..
		ofstream meso_seq_file ("mesoseq.txt");
		ofstream meso_expr_file ("mesoexpr.txt");
		
		ofstream neuro_seq_file("neuroseq.txt");
		ofstream neuro_expr_file("neuroexpr.txt");
		
		assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
		
		meso_expr_file << header << endl;
		neuro_expr_file << header << endl;
	for(int i = 0; i < size; i++){
		//ostringstream number;
		//number << i;
		if( gsl_vector_get(transition_Indices,i) < 10 ) {
				
				meso_seq_file << ">" << seq_name[i] << endl;
				meso_seq_file << seq[i] << endl;
				meso_expr_file << expr[i] << endl;	
		}
		else{
				neuro_seq_file << ">" << seq_name[i] << endl;
				neuro_seq_file << seq[i] << endl;
				neuro_expr_file << expr[i] << endl;
		
		}
	}
		neuro_seq_file.close();
		neuro_expr_file.close();
		
		meso_seq_file.close();
		meso_expr_file.close();
		seq_file.close();
		expr_file.close();


	//assert(argc == 4);
	ifstream seq_filem("mesoseq.txt");
	ifstream expr_filem("mesoexpr.txt");
	
	string headerm;
	//read the first line (header) in expression file
	getline(expr_filem, headerm);
	//cout << " headerm " << headerm.size() << endl;
	string tempm;
	vector <string> seq_namem;
	vector <string> seqm;
	vector <string> exprm;
	
	seqm.clear();
	exprm.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_filem.eof()){
		tempm = "";
		getline(seq_filem, tempm);
		if(tempm.length() == 0){
			break;
		}
		
		string namem (tempm, 1, tempm.length() - 1);
		seq_namem.push_back(namem);
		
		getline(seq_filem, tempm);
		seqm.push_back(tempm);
				
		getline(expr_filem, tempm);
		exprm.push_back(tempm);
	}
	
	int sizem = seqm.size();
	//cout << sizem << endl;
	//do cross-check for names
	for(int i = 0; i < seqm.size(); i++){
		//assert(same_name(seq_name[i], expr[i]));
	}
	
//cout << " here we test output" << endl;
	//cout << expr[1] << endl;
	//vector <int> indices (xin_indices, xin_indices + sizeof(xin_indices) / sizeof(int));

	vector <int> indicesm(sizem, 0);
	for(int i = 0; i < indicesm.size(); i++)
		indicesm[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesm.begin(), indicesm.end(), Random);		
	
	//for(int i = 0; i < indicesm.size(); i++)
		//cout << indicesm[i] << endl;

	int min_sizem = sizem/k;
	vector <int> partition_sizesm(k, min_sizem);
	
	int residuem = sizem - k * min_sizem;
	for(int i = 0; i < residuem; i++){
		partition_sizesm[i] ++;
	}	
	
	//cross check for partition sizes
	int summ = 0;
	for(int i = 0; i < k; i++){
		summ += partition_sizesm[i];
	}
	
	assert (summ == sizem);
	
	//starting index of each bin in the indices array
	int startm = 0;
	vector <int> startsm;
	startsm.clear();
	for(int i = 0; i < k; i++){
		startsm.push_back(startm);
		startm += partition_sizesm[i];
	}
	
	for(int i = 0; i < k; i++){
		//cout << startsm[i] << endl;
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numberm;
		numberm << i;
	//	cout << "about to cout i in meso" << i << endl;
		ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
		ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
		
		ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
		ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		train_expr_file << headerm << endl;
		test_expr_file << headerm << endl;
		
		for(int j = 0; j < sizem; j++){
			int indexm = indicesm[j];
			if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
				test_seq_file << ">" << seq_namem[indexm] << endl;
				test_seq_file << seqm[indexm] << endl;
				test_expr_file << exprm[indexm] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namem[indexm] << endl;
				train_seq_file << seqm[indexm] << endl;
				train_expr_file << exprm[indexm] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}

////////////////////neuro

	ifstream seq_filen("neuroseq.txt");
	ifstream expr_filen("neuroexpr.txt");

	string headern;
	//read the first line (header) in expression file
	getline(expr_filen, headern);
	
	string tempn;
	vector <string> seq_namen;
	vector <string> seqn;
	vector <string> exprn;
	
	seqn.clear();
	exprn.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_filen.eof()){
		tempn = "";
		getline(seq_filen, tempn);
		if(tempn.length() == 0){
			break;
		}
		
		string namen (tempn, 1, tempn.length() - 1);
		seq_namen.push_back(namen);
		
		getline(seq_filen, tempn);
		seqn.push_back(tempn);
				
		getline(expr_filen, tempn);
		exprn.push_back(tempn);
	}
	
	int sizen = seqn.size();
	//cout << "sizen  " << sizen << endl;
//cout << "seqn[0]  " << seqn[0] << endl;
	//do cross-check for names
	for(int i = 0; i < seqn.size(); i++){
		//assert(same_name(seq_name[i], expr[i]));
	}
	
//cout << " here we test output" << endl;
	//cout << expr[1] << endl;
	//vector <int> indices (xin_indices, xin_indices + sizeof(xin_indices) / sizeof(int));

	vector <int> indicesn(sizen, 0);
	for(int i = 0; i < indicesn.size(); i++)
		indicesn[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesn.begin(), indicesn.end(), Random);		
	
	//for(int i = 0; i < indicesn.size(); i++)
	//	cout << " indicesn " <<  indicesn[i] << endl;

	int min_sizen = sizen/k;
	vector <int> partition_sizesn(k, min_sizen);
	
	int residuen = sizen - k * min_sizen;
	for(int i = 0; i < residuen; i++){
		partition_sizesn[i] ++;
	}	
	
	//cross check for partition sizes
	int sumn = 0;
	for(int i = 0; i < k; i++){
		sumn += partition_sizesn[i];
	}
	
	assert (sumn == sizen);
	
	//starting index of each bin in the indices array
	int startn = 0;
	vector <int> startsn;
	startsn.clear();
	for(int i = 0; i < k; i++){
		startsn.push_back(startn);
		startn += partition_sizesn[i];
	}
	
	for(int i = 0; i < k; i++){
		//cout << startsn[i] << endl;
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numbern;
		numbern << i;
		//cout << "about to cout i in neuro " << i << endl;
		ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
		ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
		
		ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
		ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		//train_expr_file << header << endl;
		//test_expr_file << header << endl;
		
		for(int j = 0; j < sizen; j++){
			int indexn = indicesn[j];
			if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
				test_seq_file << ">" << seq_namen[indexn] << endl;
				test_seq_file << seqn[indexn] << endl;
				test_expr_file << exprn[indexn] << endl;				
			}
			else{
				train_seq_file<< ">" << seq_namen[indexn] << endl;
				train_seq_file << seqn[indexn] << endl;
				train_expr_file << exprn[indexn] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}
return 1;
}
int ExprPredictor::createccdata2()
{
double expmin = .35;
double expmax = .65;

int k=3;
/////////////////////////////////////////////////////
///////////////////////////////////////////////////
///// this may effect optimizer



///////////////////////////////////////////////////////
    
///////////////////////**************************8///////////////////////
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	// create rng type
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    // read the first row: row labels (ignore the first field)
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    
    // read the data
//for(int i; i < 5; i++){
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
	 rand_site_index = gsl_ran_flat( rng, -.15, .15 );
	if( val + rand_site_index > 1 ) {
	 vals.push_back( 1 );
	}
	
	if( val + rand_site_index < 0 ) {
	 vals.push_back( 0 );
	}
	
	if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
	vals.push_back( val + rand_site_index );
	}
	} // while ss >> val
        data.push_back( vals );
    }
Matrix data1( data );
//data1.save( "matrixexpr.txt" );
//////////////////////////////////////////******************************888888/////////////////////
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	//cout << " for row " << i << " max index " << mi << endl;
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
if( m < expmin ) {
	gsl_vector_set(transition_Indices, i , ti);  // here ti = NULL
	//continue;
	break;
	}
else{
	 
        for (int j=mi; j< ncol; j++)  {  
		/*
	        double temp;
		temp = m/2;     
		*/	        
		if (expmax < gsl_vector_get(rowexprData,j) ) {
		//	cout << " expmax < e " << j <<'\t' ;
			continue;
	        }
		
		else {

			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..  
}// else     
}// for i

for (int i=0; i<nrow;i++) {
//cout << "trans index : " << gsl_vector_get(transition_Indices,i) << endl;
}
/////////////////////////start with original seq and expr file (declassified)
//ifstream seq_file("efastanotw6.txt");
//	ifstream expr_file("expre6.tab");
//ifstream seq_file("efastanotw6b.txt");  //("efastanotw6.txt");
//	ifstream expr_file( "expre6.tab");//"expre6.tab");
ifstream seq_file("efastanotw10.txt");  //("efastanotw6.txt");
	ifstream expr_file( "expre10.tab");//"expre6.tab");
assert (seq_file.is_open() &&  expr_file.is_open());
string header;
	//read the first line (header) in expression file
	getline(expr_file, header);
	
	string temp;
	vector <string> seq_name;
	vector <string> seq;
	vector <string> expr;
	
	seq.clear();
	expr.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name1 (temp, 1, temp.length() - 1);
		seq_name.push_back(name1);
		
		getline(seq_file, temp);
		seq.push_back(temp);
				
		getline(expr_file, temp);
		expr.push_back(temp);
	}
	
	int size = seq.size();
//cout << " size of seq.size  " << size << endl;
//cout << seq[0] << endl;
/////////////////////////create meso and neur files based on cell attribute of exprpredictor class (loop over nseqs, check cell, if cell < 8, append seq and exp to meso. etc..
		ofstream meso_seq_file ("mesoseq.txt");
		ofstream meso_expr_file ("mesoexpr.txt");
		
		ofstream neuro_seq_file("neuroseq.txt");
		ofstream neuro_expr_file("neuroexpr.txt");
		
		assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
		
		meso_expr_file << header << endl;
		neuro_expr_file << header << endl;
//cout << " size " << size << endl;
	for(int i = 0; i < size; i++){
		//ostringstream number;
		//number << i;
		ostringstream numberm;
		numberm << i;
		
		//ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
		if( gsl_vector_get(transition_Indices,i) < 10 ) {
				
				meso_seq_file << ">" << seq_name[i] << endl;
				meso_seq_file << seq[i] << endl;
				meso_expr_file << expr[i] << endl;	
		}
		else{
				neuro_seq_file << ">" << seq_name[i] << endl;
				neuro_seq_file << seq[i] << endl;
				neuro_expr_file << expr[i] << endl;
		
		}
	}
		neuro_seq_file.close();
		neuro_expr_file.close();
		
		meso_seq_file.close();
		meso_expr_file.close();
		seq_file.close();
		expr_file.close();


	//assert(argc == 4);
	ifstream seq_filem("mesoseq.txt");
	ifstream expr_filem("mesoexpr.txt");
	
	string headerm;
	//read the first line (header) in expression file
	getline(expr_filem, headerm);
	//cout << " headerm " << headerm.size() << endl;
	string tempm;
	vector <string> seq_namem;
	vector <string> seqm;
	vector <string> exprm;
	
	seqm.clear();
	exprm.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_filem.eof()){
		tempm = "";
		getline(seq_filem, tempm);
		if(tempm.length() == 0){
			break;
		}
		
		string namem (tempm, 1, tempm.length() - 1);
		seq_namem.push_back(namem);
		
		getline(seq_filem, tempm);
		seqm.push_back(tempm);
				
		getline(expr_filem, tempm);
		exprm.push_back(tempm);
	}
	
	int sizem = seqm.size();
	//cout << sizem << endl;
	//do cross-check for names
	for(int i = 0; i < seqm.size(); i++){
		//assert(same_name(seq_name[i], expr[i]));
	}
	
//cout << " here we test output" << endl;
	//cout << expr[1] << endl;
	//vector <int> indices (xin_indices, xin_indices + sizeof(xin_indices) / sizeof(int));

	vector <int> indicesm(sizem, 0);
	for(int i = 0; i < indicesm.size(); i++)
		indicesm[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesm.begin(), indicesm.end(), Random);		
	
	//for(int i = 0; i < indicesm.size(); i++)
		//cout << indicesm[i] << endl;

	int min_sizem = sizem/k;
	vector <int> partition_sizesm(k, min_sizem);
	
	int residuem = sizem - k * min_sizem;
	for(int i = 0; i < residuem; i++){
		partition_sizesm[i] ++;
	}	
	
	//cross check for partition sizes
	int summ = 0;
	for(int i = 0; i < k; i++){
		summ += partition_sizesm[i];
	}
	
	assert (summ == sizem);
	
	//starting index of each bin in the indices array
	int startm = 0;
	vector <int> startsm;
	startsm.clear();
	for(int i = 0; i < k; i++){
		startsm.push_back(startm);
		startm += partition_sizesm[i];
	}
	
	for(int i = 0; i < k; i++){
		//cout << startsm[i] << endl;
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numberm;
		numberm << i;
		
		ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
		ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
		
		ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
		ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		train_expr_file << headerm << endl;
		test_expr_file << headerm << endl;
		
		for(int j = 0; j < sizem; j++){
			int indexm = indicesm[j];
			if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
				test_seq_file << ">" << seq_namem[indexm] << endl;
				test_seq_file << seqm[indexm] << endl;
				test_expr_file << exprm[indexm] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namem[indexm] << endl;
				train_seq_file << seqm[indexm] << endl;
				train_expr_file << exprm[indexm] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}

////////////////////neuro

	ifstream seq_filen("neuroseq.txt");
	ifstream expr_filen("neuroexpr.txt");

	string headern;
	//read the first line (header) in expression file
	getline(expr_filen, headern);
	
	string tempn;
	vector <string> seq_namen;
	vector <string> seqn;
	vector <string> exprn;
	
	seqn.clear();
	exprn.clear();
	
	//read the files and load the sequences in vector seq, expressions in vector expr
	while(!seq_filen.eof()){
		tempn = "";
		getline(seq_filen, tempn);
		if(tempn.length() == 0){
			break;
		}
		
		string namen (tempn, 1, tempn.length() - 1);
		seq_namen.push_back(namen);
		
		getline(seq_filen, tempn);
		seqn.push_back(tempn);
				
		getline(expr_filen, tempn);
		exprn.push_back(tempn);
	}
	
	int sizen = seqn.size();
	//cout << "sizen  " << sizen << endl;
//cout << "seqn[0]  " << seqn[0] << endl;
	//do cross-check for names
	for(int i = 0; i < seqn.size(); i++){
		//assert(same_name(seq_name[i], expr[i]));
	}
	
//cout << " here we test output" << endl;
	//cout << expr[1] << endl;
	//vector <int> indices (xin_indices, xin_indices + sizeof(xin_indices) / sizeof(int));

	vector <int> indicesn(sizen, 0);
	for(int i = 0; i < indicesn.size(); i++)
		indicesn[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesn.begin(), indicesn.end(), Random);		
	
	//for(int i = 0; i < indicesn.size(); i++)
		//cout << " indicesn " <<  indicesn[i] << endl;

	int min_sizen = sizen/k;
	vector <int> partition_sizesn(k, min_sizen);
	
	int residuen = sizen - k * min_sizen;
	for(int i = 0; i < residuen; i++){
		partition_sizesn[i] ++;
	}	
	
	//cross check for partition sizes
	int sumn = 0;
	for(int i = 0; i < k; i++){
		sumn += partition_sizesn[i];
	}
	
	assert (sumn == sizen);
	
	//starting index of each bin in the indices array
	int startn = 0;
	vector <int> startsn;
	startsn.clear();
	for(int i = 0; i < k; i++){
		startsn.push_back(startn);
		startn += partition_sizesn[i];
	}
	
	for(int i = 0; i < k; i++){
		//cout << startsn[i] << endl;
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numbern;
		numbern << i;
		
		ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
		ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
		
		ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
		ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		//train_expr_file << header << endl;
		//test_expr_file << header << endl;
		
		for(int j = 0; j < sizen; j++){
			int indexn = indicesn[j];
			if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
				test_seq_file << ">" << seq_namen[indexn] << endl;
				test_seq_file << seqn[indexn] << endl;
				test_expr_file << exprn[indexn] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namen[indexn] << endl;
				train_seq_file << seqn[indexn] << endl;
				train_expr_file << exprn[indexn] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}
return 1;
}
void ExprPredictor::compvar(vector< double >& vars) // added from Binsrep folder check READMEd there
{
	vector< double > predictedExprs;
        vector< double > stepobservedExprs;

    // create the expression function
    ExprFunc* func = createExprFunc( par_model );
            double rss = 0;
    // error of each sequence
double step = .001;
    vector< double > concs = factorExprData.getCol( 0 );

vector< double > concs2 = concs;
concs2[0] = concs[0] - step;
//cout <<" concs[0] = "<<    concs[0]   <<  endl;
    for ( int i = 0; i < nSeqs(); i++ ) {
        
       
         
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr( seqSites[ i ], 5, concs );
            predictedExprs.push_back( predicted );
           
            // observed expression for the i-th sequence at the j-th condition
            double stepobserved = func->predictExpr(seqSites[ i ], 5, concs2 );       //exprData( i, 0 );
            stepobservedExprs.push_back( stepobserved );
        
     }
double partialcu = .5 * ( concs2[0] + concs[0] );
assert( predictedExprs.size() == stepobservedExprs.size());
    for ( int i = 0; i <  predictedExprs.size(); i++ ) {
 //cout << "predictedExprs[i] = " <<    predictedExprs[i]  <<  endl;
 //cout << " stepobservedExprs[i] = " <<    stepobservedExprs[i] <<  endl;
//cout << " predictedExprs[i]  - stepobservedExprs[i] = " << predictedExprs[i]  -    stepobservedExprs[i] <<  endl;
//cout << " partialcu * (predictedExprs[i]  -  stepobservedExprs[i]  ) / step  = "  << partialcu * (predictedExprs[i]  -  stepobservedExprs[i]  ) / step << endl;
     vars.push_back(  partialcu * (predictedExprs[i]  -  stepobservedExprs[i]  ) / step )   ; 
    }
//cout << "concs[0] = " <<     concs[0]  <<  endl;
//cout << "concs2[0] = " <<     concs2[0]  <<  endl;
//cout << "partialcu = " <<     partialcu  <<  endl;
}
double ExprPredictor::compAvgCrossCorr( const ExprPar& par ) const
{
    // create the expression function
    ExprFunc* func = createExprFunc( par );
            
    // cross correlation similarity of each sequence
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            // predicted expression for the i-th sequence at the j-th condition
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            // observed expression for the i-th sequence at the j-th condition
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs ); 
    }	

    return totalSim / nSeqs();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ExprPredictor::compOccMat(const gsl_vector* v, void* params  ) 
{
int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
        for (int j=mi; j< ncol-mi; j++)  {  
	        double temp;
		temp = m/2;     
	        if (temp < gsl_vector_get(rowexprData,j) ) {
			continue;
	        }
		
		else {
			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..  
}     

ExprFunc* funcmat = createExprFunc2( par_model );
gsl_matrix * Occupancy = gsl_matrix_alloc(nrow,par_model.nFactors() );  //   this should be  done in the optimization do loop     
//gsl_matrix * Occw = gsl_matrix_alloc(nrow,4 );
for ( int i = 0; i < nSeqs(); i++ ) {
                                                                  // for loop i index = transition_Indices(i), at each transition getCol of transFactor exp;
        vector< double > concs = factorExprData.getCol( gsl_vector_get(transition_Indices,i) );
        gsl_vector * Occvector =gsl_vector_alloc(motifs.size()); // there are 4 w's
        vector <double > fOcc(3);
        funcmat->predictOcc( seqSites[ i ], seqLengths[i], concs, fOcc); // don't need TranI or i, since concs only has the transition
        Occvector = vector2gsl(fOcc);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
        gsl_matrix_set_row(Occupancy , i , Occvector);           
}


int fixedsize = fix_pars.size();  
int rows = Occupancy -> size1;
int columns = Occupancy -> size2;
//int rows = Occw -> size1;
//int columns = Occw -> size2;
int i,j,k;
k=0;
gsl_matrix *X = gsl_matrix_alloc(columns,columns);
gsl_matrix *V = gsl_matrix_alloc(columns,columns);
gsl_vector *S =gsl_vector_alloc(columns);
gsl_vector *xx =gsl_vector_alloc(columns);  // x is in the row space, therefore it is a vector in Rcolumns, with at most row dimensions.
gsl_vector *b =gsl_vector_alloc(rows);     // b is in the column space of A.
gsl_vector_set_all( b,0 );
gsl_vector *work=gsl_vector_alloc(columns);
int rsvd;		// A must have more rows than columns or the same num
int rsvds;               //(A,V,S,work), on output A is replaced by U  (A=USVt)
rsvd=gsl_linalg_SV_decomp(Occupancy,V,S,work);
//rsvds=gsl_linalg_SV_solve(Occupancy,V,S,b,xx);
//gsl_matrix_transpose(V);
//printf ("x = \n");
//gsl_vector_fprintf (stdout, xx, "%g");
/*
for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,j) == 0) {
		 fix_pars.clear(); 
		for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			fix_pars.push_back(gsl_matrix_get(V,j,i));
		}
		break;
	
	}

}
*/
//fix_pars.clear();
//for(int j =0; j< fixedsize; j++) {
	if(gsl_vector_get(S,2) < .1 ) {
	//	 fix_pars.clear(); 
	//	for(int i =0; i < fixedsize;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			//fix_pars.push_back(gsl_matrix_get(V,2,j));  // this is causing some parameters to be out of range!!!
			cout << "gsl_matrix_get(V,2,j)   = " << gsl_matrix_get(V,2,j)  << endl;
			//par_model.txpEffects[j] = gsl_matrix_get(V,2,j);  // these parameters may be out of range and will fuck things up !!
		}
	//else break;
	
//}
for(int i =0; i < columns;i++){  // the number of fix_pars must be aligned with the seq2e main file..
			cout << "gsl_vector_get(S,i)  , singular values  = " << gsl_vector_get(S,i)  << endl;
		}
gsl_matrix_free( Occupancy );
gsl_matrix_free( X );
gsl_matrix_free( V );
gsl_vector_free( xx );
gsl_vector_free( b );
gsl_vector_free( S );
//printf("rsvd = %d\n",rsvd);    	
//printf("rsvd = %d\n",rsvds);
   // return OccMatrix
}

//////////////////////////////////////////////////////////
   

int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )   // const
{

 cout << "Start minimization simplex" << endl;

    vector< double > pars;
  
printPar(par_model );
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
       cout << pars.size() <<  "if no num to left pars is empty:" <<endl; 
int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
//cout << endl;
//cout << "  pars " << endl << pars << endl;
printPar( par_model );
return 1;
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       

vector < double >free_parste(0) ;
free_parste= free_pars;
vector < double >fix_parste(0);
fix_parste=fix_pars;
//cout << pars.size() <<  "if n8 num to left pars is empty:" <<endl; 
//cout << free_pars.size() <<  "if n8 num to left freepars is empty:" <<endl;
//compOccMat( s-> x, (void*)this);
	pars.clear();
	int free_par_counterte = 0;
	int fix_par_counterte = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_parste[ free_par_counterte ++ ]);
		}
		else{
			pars.push_back( fix_parste[ fix_par_counterte ++ ]);
		}
	}
 par_result = ExprPar( pars, coopMat, actIndicators, repIndicators );
return 1;
//cout <<"about to inverse transform pars to exprpar " << endl;
//printPar( par_result);
	//Hassan start:
	
	
	pars.clear();
	pars = free_pars;
	//Hassan end
//cout << " pars.size " << pars.size() << endl;
//cout << " pars " << endl << pars << endl;
for (int i = 0; i < pars.size() ; i++ ) {
//cout << " pars"<< i << '\t'  << pars[i] << endl;
}
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f3;
    my_func.n = pars.size();
    my_func.params = (void*)this;    // params holds par_model and so does pars.., but these are essentially distinct objects, so one could change but the other not
    
    gsl_vector* x = vector2gsl( pars );

   const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );

    // create the minimizer
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );

    gsl_multimin_fminimizer_set( s, &my_func, x, ss ); // here is where the problem is  101
cout << "ssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
//count = 0;

    // iteration
    size_t iter = 0;
    int status;
    double size;	
printPar( par_model );
// gsl_matrix* svdmatrix = gsl_matrix_alloc(nSeqs(),fix_pars.size());
//compOccMat( s-> x, (void*)this);
    do {
        double f_prev = iter ? s->fval : 1.0E6;     // the function starts with some very large number
   
	iter++;

        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_fminimizer_iterate( s );
   //  cout <<"iter " << iter << endl;
        // check for error
        if ( status ) { cout<< " gsl_multimin_iterate positve: " << endl; break;}
//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
        ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
	//Hassan end
	  // check if the current values of parameters are valid
      
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) {cout << "testpar failed: " <<endl;  break;}
        
        // check for stopping condition
         double f_curr = s->fval;            // 525 this causes the parameters to be stable or unstable if uncommented
       //  double delta_f = abs( f_curr - f_prev ); 
        // if ( objOption == SSE && delta_f < min_delta_f_SSE ) { cout << "min f_sse " << endl ; break;}
        // if ( objOption == CORR && delta_f < min_delta_f_Corr ) {cout << "min_delta_f_corr positive : " << endl; break;}
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size,1e-1 ); //1e-6
 
 //status = gsl_multimin_test_size( size, 1e-1 );
//cout << "status   " << status  << endl;
 		if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }

        // print the current parameter and function values
      cout <<"iter" << "\t" << iter << "\t" << " f_curr " << f_curr<<endl;
  //  printPar( par_curr );
   //   printf( "\tf() = %8.5f size = %.3f\n", s->fval, size );
    } while (  status == GSL_CONTINUE && iter < nSimplexIters );
 
//Hassan start:
free_pars = gsl2vector( s-> x);
cout << pars.size() <<  "if n8 num to left pars is empty:" <<endl; 
printPar(par_model );
//cout << free_pars.size() <<  "if n8 num to left freepars is empty:" <<endl;
//compOccMat( s-> x, (void*)this);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
 par_result = ExprPar( pars, coopMat, actIndicators, repIndicators );
    obj_result = s->fval;
cout << "End minimization simplex" << endl;
    return 0;

}

int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )   // const
{
	cout << "Start gradient minimiz minimization" << endl;
printPar( par_model );
    // extract initial parameters
    vector< double > pars;
//par_model.adjust();
printPar( par_model );
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
       //Hassan start:
printPar( par_model );
cout<< "nopeeeeeeeee"<<endl;
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			//cout << "testing 1: " << pars[ index ] << endl;
			free_pars.push_back( pars[ index ]);
		}
		else{
			//cout << "testing 2: " << pars[ index ] << endl;
			fix_pars.push_back( pars[ index ] );
		}
	}
  
	
	pars.clear();
	pars = free_pars;
	//Hassan end     
    // set the objective function and its gradient
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f3;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    // set the initial values to be searched 
    gsl_vector* x = vector2gsl( pars ); 

	// CHECK POINT: evaluate gsl_obj_f() function
 	//cout << "binding at the initial value of parameters = " << gsl_obj_f( x, (void*)this ) << endl; 
		
	// choose the method of optimization and set its parameters
 	//const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
//    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs;
		const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_steepest_descent;
    // create the minimizer
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
//////////////////////////////////////////////
///////////////////////////////////////////////
/// if we only are concerned with 1 significant figure for each parameter, can we change the tol, and init_step to larger values?
//  originally they were:   double init_step = 0.02, tol = 0.1;
//  changed 81411
//////////////////////////////////////////////
    double init_step = .01, tol = .0001;
//cout << " before set " << endl;
//printPar( par_model );
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
    //   cout << "after set " << endl;
//printPar( par_model );     
    // iteration
    size_t iter = 0;
    int status;
    do {
        double f_prev = iter ? s->f : 1.0E6;     // the function starts with some very large number	
      // cout <<"iter " << iter << endl;
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
 //cout << " status " << status << endl;
//cout << " GSL_SUCCESS " << GSL_SUCCESS << endl;
//cout << " GSL_CONTINUE " << GSL_CONTINUE << endl;
	if (status ) {
//cout << "gsl_multimin_fdfminimizer_iterate( s ) positive " << endl;
//printf("%3d  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %10.5f\n", iter,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2),gsl_vector_get(s->x,3),gsl_vector_get(s->x,4),gsl_vector_get(s->x,5),gsl_vector_get(s->x,6),gsl_vector_get(s->x,7),gsl_vector_get(s->x,8),s->f);
 break; }
cout << " about to printPar( par_model ) " << endl;
     printPar( par_model );
//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
 ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
//Hassan end
cout << " about to printPar( par_curr ) " << endl;
     printPar( par_curr );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) { cout << "testpar failed: " <<endl;  printPar( par_curr );break;}
        double f_curr = s->f;
        double delta_f = abs( f_curr - f_prev ); 
cout << " deltaf " << delta_f << " f_curr " << f_curr << " iter " << iter << endl;
        //if ( objOption == SSE && delta_f < min_delta_f_SSE ) {cout << "min f_sse " << endl ; break;}
	 if ( objOption == SSE && delta_f < .0000001 ) {cout << "min f_sse " << endl ; break;}
        if ( objOption == CORR && delta_f < min_delta_f_Corr ) { cout << "min_delta_f_corr" << endl; break;}
        if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
  //      cout << " about to :status = gsl_multimin_test_gradient( s->gradient, .001 ) (if status=0 break)" << endl ;
        status = gsl_multimin_test_gradient( s->gradient, 5e-4 ); // was 5e-4, 81611
    //    cout << " status " << status << endl;
      // printPar( par_curr );

    } while ( status == GSL_CONTINUE && iter < nGradientIters );

//Hassan start:
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators );
	//Hassan end
    obj_result = s->f;
cout << " par _result " << endl;
printPar( par_result );
  // free the minimizer
    gsl_vector_free( x );    
    gsl_multimin_fdfminimizer_free( s );
    
    return 0;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
//printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
printf ("iter: %3u x = % 15.8f % 15.8f "
"|f(x)| = %g\n",
iter,
gsl_vector_get (s->x, 0),
gsl_vector_get (s->x, 1),
//gsl_vector_get (s->x, 2),
gsl_blas_dnrm2 (s->f));
}

int ExprPredictor::gradient_minimize2( ExprPar& par_result, double& obj_result )   // const
{
	cout << "Start gradient minimiz minimization" << endl;
    // extract initial parameters
    vector< double > pars;

    par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators ); 
    cout << "pars " << pars << endl;

	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			//cout << "testing 1: " << pars[ index ] << endl;
			free_pars.push_back( pars[ index ]);
		}
		else{
			//cout << "testing 2: " << pars[ index ] << endl;
			fix_pars.push_back( pars[ index ] );
		}
	}
  	pars.clear();
	pars = free_pars;
	
const gsl_multifit_fdfsolver_type *T;
gsl_multifit_fdfsolver *s;
int status;
unsigned int i, iter = 0;
const size_t n = nSeqs()*nConds();
const size_t p = pars.size();
cout << " n " << n << " p " << p << endl;
    gsl_multifit_function_fdf my_func;
    my_func.f = &gsl_obj_f3_fit;
cout << " aoub " << endl;	
   // my_func.df =NULL; // &gsl_obj_df_fit;
my_func.df = &gsl_obj_df_fit;
cout << " aee " << endl;
    my_func.fdf =&gsl_obj_fdf_fit;
cout << " noo " << endl;
    my_func.n = n ;  
    my_func.p = p ; 
    my_func.params = (void*)this;


    gsl_vector* xx = vector2gsl( pars ); 

T = gsl_multifit_fdfsolver_lmder;
s = gsl_multifit_fdfsolver_alloc (T, n, p);
cout << " ok a " << endl;
gsl_multifit_fdfsolver_set (s, &my_func, xx);
const gsl_matrix* a;//(s->J);
cout << " rint " << endl;
//print_state (iter, s);
bool aa=1;
do
{
iter++;
cout << "here5 " << endl;
status = gsl_multifit_fdfsolver_iterate(s);
printf ("status = %s\n", gsl_strerror (status));
//printf ("x = \n");
//gsl_vector_fprintf (stdout, x, "%g");
//printf("%3d  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %10.5f\n", iter,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2),gsl_vector_get(s->x,3),gsl_vector_get(s->x,4),gsl_vector_get(s->x,5),gsl_vector_get(s->x,6),gsl_vector_get(s->x,7),gsl_vector_get(s->x,8),s->f);
//cout << " iter, x, dx, f " << endl;
//printf("%3u  %.3f  %.3f  %g\n", iter,gsl_vector_get(s->x,0),gsl_vector_get(s->dx,0),s->f);
//print_state (iter, s);
//if (status)
//break;

gsl_vector* g=gsl_vector_alloc( s->J->size2 );
const gsl_matrix* a(s->J);
const gsl_vector* fa(s->f);
//gsl_multifit_fdfsolver_dif_df (s->x, &gsl_obj_fdf_fit, s->f, s->J);
//gsl_matrix_fprintf (stdout, a, "%g");
/*
for(int j = 0; j < a->size1 ; j++ ) {
		gsl_vector * Hvector =gsl_vector_alloc(a->size2 ); 
	       gsl_matrix_get_row( Hvector ,a , j );
		for(int i = 0; i < a->size2 ; i++ ) {
		cout << gsl_vector_get(Hvector,i)  << '\t' ;
		}
	cout << endl;
	
	}
cout << " ok " <<endl;
cout << " size J rows " << s->J->size1 << endl;
cout << "size J cols " << s->J->size2 << endl;
cout << " size f " << s->f->size << endl;
*/
		//cout << gsl_vector_get(s->x,0)  << '\t' ;
		//cout << gsl_vector_get(s->dx,0)  << '\t' ;
	cout << endl;
	
cout << " okinside grad minimize 3 " <<endl;
int aaa=0;
aaa =gsl_multifit_gradient(a,fa, g) ;
gsl_vector_fprintf (stdout, g, "%g");
cout << " grad above " << endl;
double epsabs=.001;
status= gsl_multifit_test_gradient (g, epsabs);
//printf ("status grad = %s\n", gsl_strerror (status));
status = gsl_multifit_test_delta (s->dx, s->x,1e-4, 1e-4);
gsl_vector_fprintf (stdout, s->dx, "%g");
gsl_vector_fprintf (stdout, s->x, "%g");
//printf ("status = %s\n", gsl_strerror (status));
cout << "here5b " << endl;
}
while (status == GSL_CONTINUE && iter < 15);
//while (status == GSL_CONTINUE && iter < 500);
//gsl_multifit_fdfsolver_dif_df (s->x, s, s->f, a);
free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators,aa );
	//Hassan end
    

gsl_matrix *covar = gsl_matrix_alloc (p, p);
gsl_multifit_covar (s->J, 0.0, covar);
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
{
double chi = gsl_blas_dnrm2(s->f);
double dof = n - p;
double c = GSL_MAX_DBL(1, chi / sqrt(dof));
obj_result = chi;
printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
printf("chisq = %g\n", chi);
printf ("par1 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
printf ("par2 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
//printf ("w_(tw) = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
//printf ("w_(sn) = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
}
printf ("status = %s\n", gsl_strerror (status));

cout << "covar = JtJ inverse, sqrt((JtJ)^-1) " << endl;   /// why does this say inverse?
for( int ii =0 ; ii < covar->size1; ii++) {
	for( int j =0 ; j < covar->size2; j++) {
	//if(i==j ){
	cout <<gsl_matrix_get(covar,ii,j) << '\t';	
	//}
	}
cout << endl;
}

gsl_multifit_fdfsolver_free (s);
  gsl_vector_free( xx );    

   

 
    
    return 0;
}



double gsl_obj_f( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
{ 

    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object. 
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
    if( (*predictor).clasvar == 0 ) {
double tempobj2 = predictor->objFuncborder2( par ); return tempobj2;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
}
else { 
//cout << "classvar " <<  (*predictor).clasvar << endl;
double obj = predictor->objFunc2( par );  return obj;	
}

    // call the ExprPredictor object to evaluate the objective function 
//    cout << "obj2= nnaaaa" << endl;
double tempobj2 = predictor->objFuncborder2( par ); // 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.  // 111611 changed objFuncborder2 to objFuncborder, and changed the negative sign in front of compavgcorrborder.
//double obj = predictor->objFunc( par );	

	
 //if (abs(tempobj2) > abs(obj2) ) {   par_result = par_model ; obj2 = tempobj2; }
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
	//double objsum = obj + tempobj2;    
	return tempobj2;
}
double gsl_obj_f3( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
{ 

    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object. 
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
  

if( (*predictor).clasvar == 1 ) {
double tempobj3 = predictor->objFuncborder( par ); return tempobj3;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
}
else { 
double tempobj3 = predictor->compRMSE2( par );   //objFuncborder2( par ); // 
	return tempobj3;	
}

}
double gsl_obj_f2( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
{ 

    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object. 
    //ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
    
    // call the ExprPredictor object to evaluate the objective function 
//    cout << "obj2= nnaaaa" << endl;
if( (*predictor).clasvar == 1 ) {
double tempobj2 = predictor->objFuncborder2( par ); return tempobj2;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
}
else { 
double obj = predictor->objFunc( par );  return obj;	
}
	
 //if (abs(tempobj2) > abs(obj2) ) {   par_result = par_model ; obj2 = tempobj2; }
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
//cout << "obj2= " << "\t" << "tempobj2= "  << abs(tempobj2)  << endl;
	//double objsum = obj + tempobj2;    
	
}
/*
double gsl_obj_f( const gsl_vector* v, void* params )
{ 
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;// apparently v gets passed as opposed to par_model, because gls_multimin_iter doesn't know of par_model so it only iterates on the v vector, which means v acts as a proxy for pars in ExprPar parameters(vector pars,....)
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
    ExprPar par( gsl2vector( v ), predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
//vector < double > x;
//x=gsl2vector(v);        
//for(int i=0; i < x.size(); i++){ cout << x[i] << endl;
// cout << "sizex " << x.size() << endl;
//}   
//cout <<  " hello " << endl; 

    // call the ExprPredictor object to evaluate the objective function 
    double obj = predictor->objFunc( par );
   //predictor->*this.printPar( par );	
//cout << "didn't make it " << endl;
    return obj;
}
*/
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-3;
    numeric_deriv( grad, gsl_obj_f3, v, params, step );	
}

void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f3( v, params ); 
    gsl_obj_df( v, params, grad );		
}

int gsl_obj_f3_fit( const gsl_vector* v, void* params, gsl_vector * f) //place compmatocc in here since we cant access v elsewhere
{ 

    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	vector<double> fix_pars(predictor->getfixpars());
	//vector<bool> free_pars(predictor->getfreepars());
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );


 
      
  //////////////////////
///////////////////////////
///////////////////////////
	vector<double>ff(0);
	predictor -> compf( par, ff);
	//f = gsl_vector_alloc(predictor->nConds()*i+j);
	//cout << " ff " << ff << endl;
	//f = vector2gsl(ff);
	for ( int i = 0; i < ff.size(); i++ ) gsl_vector_set( f, i, ff[ i ] );	
////////////////////////////////
/////////////////////////
	
	return GSL_SUCCESS;	
}



int gsl_obj_df_fit( const gsl_vector* v, void* params, gsl_matrix * J )
{
   //J =NULL;	  // see page  416 of gsl manual

 	ExprPredictor* predictor = (ExprPredictor*)params;
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	vector< double > fix_pars(predictor->getfixpars());
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	bool a = 1;
    ExprPar par22( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );

	vector< double > pars;
	pars=all_pars;
	//Hassan start:
	pars.clear();
	pars =  temp_free_pars;
	double step = .2;
	
vector< double > initilizer( pars.size() );
vector< vector< double > > init2( pars.size(), initilizer );
int count = 0; // used to indicate data point within bords*seq
///////////////////////////////////////////////////////
//J = gsl_matrix_alloc( predictor->nSeqs()*predictor->nConds(),pars.size() ); 
for(int i=0; i < predictor->nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
	for(int j=0; j < predictor->nConds(); j++ ) {
	for (int jj = 0; jj < pars.size() ; jj++ ) {

		
			vector< double > parsjuh =pars;
			vector< double > parsjdh= pars;
		
			parsjuh[jj] = pars[jj] + step;
			parsjdh[jj] = pars[jj] - step;
	// k fixed j up		
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	vector< bool > indicator_bool(predictor->getIndicatorbool(  ));
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}
	const	IntMatrix coopMat( predictor->getcoopmat());
	const vector< bool > actIndicators(predictor->getactIndicators());
	//actIndicators=predictor->getactIndicators() 
	const vector< bool > repIndicators(predictor->getrepIndicators()); 
	 ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);
	///// k fixed j down
	 all_pars.clear();
	 free_par_counter = 0;
	 fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}


	 ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );


	///// calculate fxx or fxy
//cout.setf( ios::fixed );
 //   cout.precision( 3 ); 
//     cout.width( 8 );
 vector< double > fOcc(0);
//hes.setf( ios::fixed);
//hes.width(3);
	count = predictor->nSeqs();
			//ExprFunc* func = this->createExprFunc( parjuh );
 		        double pjuh = predictor->compRMSE4(parjuh, i,j); // func->predictExpr3(seqSitesm1[i], 5, concs );
			// func = this->createExprFunc( parjdh );
 		        double pjdh =  predictor->compRMSE4(parjdh, i,j); //func->predictExpr3(seqSitesm1[i], 5, concs );
			//ExprFunc* func = this->createExprFunc( par_model );
			//if(i==0){
			 //func->predictOcc( seqSitesm1[i], 5, concs,  fOcc );
			//double p=func->predictExpr3(seqSitesm1[i], 5, concs ); // seqSitesm1[i] = seqSites[i]
			//occm << fOcc << endl;
			//} //if i ==0

 		        //double p = func->predictExpr3(seqSitesm1[i], 5, concs ); // seqSitesm1[i] = seqSites[i]
			//hes  <<  setprecision(2) << ( pjuh- pjdh )/(2*step) << '\t' ;
			//Jacobian[jj] += ( pjuh- pjdh )/(2.0*step);
			//nConds()*i+j
			gsl_matrix_set( J,predictor->nConds()*i+j, jj, ( pjuh- pjdh )/(2.0*step) );  // i is seqs, jj is parameters
			//count++;  // using counter to indicate position in bords*seq
			parsjuh.clear();
			parsjdh.clear();
//cout << " nottttttt "  << endl;
		} //for j nconds
	}  // for jj npars
   }// for i nseqs

//J=NULL;
return GSL_SUCCESS;
}

int gsl_obj_fdf_fit(  const gsl_vector* v, void* params, gsl_vector* f, gsl_matrix * J)
{
   // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    // parse the variables (parameters to be optimized)	
//     vector< double > expv;
//     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	vector<double> fix_pars(predictor->getfixpars());
	//vector<bool> free_pars(predictor->getfreepars());
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	bool a=1;
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() ,a);


 
      
  cout << " inside here " << endl;
vector<double >ff(0);
	predictor -> compf( par, ff);
	//cout << " ff " << endl; 
	//cout << ff << endl;
	//f = vector2gsl(ff);
	 for ( int i = 0; i < ff.size(); i++ ) gsl_vector_set( f, i, ff[ i ] );


/////////////////////////////////////////////////////////////////////////comput J
 
	vector< double > pars;
	pars=all_pars;
	//Hassan start:
	pars.clear();
	pars =  temp_free_pars;
	double step = .2;
	//bool a = 1;
vector< double > initilizer( pars.size() );
vector< vector< double > > init2( pars.size(), initilizer );
int count = 0; // used to indicate data point within bords*seq
///////////////////////////////////////////////////////
//J = gsl_matrix_alloc( predictor->nSeqs()*predictor->nConds(),pars.size() ); 
for(int i=0; i < predictor->nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
	for(int j=0; j < predictor->nConds(); j++ ) {
	for (int jj = 0; jj < pars.size() ; jj++ ) {

		
			vector< double > parsjuh =pars;
			vector< double > parsjdh= pars;
		
			parsjuh[jj] = pars[jj] + step;
			parsjdh[jj] = pars[jj] - step;
	// k fixed j up		
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	vector< bool > indicator_bool(predictor->getIndicatorbool(  ));
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}
	const	IntMatrix coopMat( predictor->getcoopmat());
	const vector< bool > actIndicators(predictor->getactIndicators());
	//actIndicators=predictor->getactIndicators() 
	const vector< bool > repIndicators(predictor->getrepIndicators()); 
	// ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators);
	ExprPar parjuh =ExprPar( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );
	//cout << " all pars parjuh " << all_pars << endl;
	///// k fixed j down
	 all_pars.clear();
	 free_par_counter = 0;
	 fix_par_counter = 0;
	for( int index = 0; index < indicator_bool.size(); index ++ ){
		if(  indicator_bool[ index ]  ){
			all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
		}
	}


	// ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators );
	ExprPar parjdh =ExprPar( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );
	//cout << " all pars parjdh " << all_pars << endl;

	///// calculate fxx or fxy
//cout.setf( ios::fixed );
 //   cout.precision( 3 ); 
//     cout.width( 8 );
 vector< double > fOcc(0);
//hes.setf( ios::fixed);
//hes.width(3);
	count = predictor->nSeqs();
			//ExprFunc* func = this->createExprFunc( parjuh );
			//cout << " i " << i << " j " << j << endl; 
 		        double pjuh = predictor->compRMSE4(parjuh, i,j); // func->predictExpr3(seqSitesm1[i], 5, concs );
			//cout << " pjuh " << pjuh << endl;
			// func = this->createExprFunc( parjdh );
 		        double pjdh =  predictor->compRMSE4(parjdh, i,j); //func->predictExpr3(seqSitesm1[i], 5, concs );
			//cout << " pduh " << pjdh << endl;
			//ExprFunc* func = this->createExprFunc( par_model );
			//if(i==0){
			 //func->predictOcc( seqSitesm1[i], 5, concs,  fOcc );
			//double p=func->predictExpr3(seqSitesm1[i], 5, concs ); // seqSitesm1[i] = seqSites[i]
			//occm << fOcc << endl;
			//} //if i ==0

 		        //double p = func->predictExpr3(seqSitesm1[i], 5, concs ); // seqSitesm1[i] = seqSites[i]
			//hes  <<  setprecision(2) << ( pjuh- pjdh )/(2*step) << '\t' ;
			//Jacobian[jj] += ( pjuh- pjdh )/(2.0*step);
			gsl_matrix_set( J,predictor->nConds()*i+j, jj, ( pjuh- pjdh )/(2.0*step) );  // i is seqs, jj is parameters
			//count++;  // using counter to indicate position in bords*seq
			parsjuh.clear();
			parsjdh.clear();
//cout << " nottttttt "  << endl;
		} //for j nconds
	}  // for jj npars
   }// for i nseqs

//J=NULL;

	 for ( int i = 0; i < predictor->nSeqs(); i++ ) {
        for ( int j = 0; j < predictor->nConds(); j++ ) {
//	 for ( int jj = 0; jj < pars.size(); jj++ ) {
	//size_t iii;
     //  for (  iii = 0; iii < ff.size(); iii++ ){
//cout<< setprecision(6) <<" gsl J "<<" predictor->nConds()*i+j " <<predictor->nConds()*i+j<< "   " << gsl_matrix_get( J,predictor->nConds()*i+j, jj)<< endl ;
//	}
    	}
        }
//cout <<" Jacobian after J "  << endl ;





	return GSL_SUCCESS;	
			
}
