#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H


#include <cctype>
#include <cstring>
#include <algorithm>
#include "Tools.h"
#include <stdio.h>


const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' };
const int NBASES = 4;	
const int ALPHABET_SIZE = 6;
const int MISSING = 4;		
const int GAP = 5;		

bool isNt( int a );
int complement( int a );
int symbolToInt( char c );
char strand2char( bool strand );
bool char2strand( char c );
vector< double > createNtDistr( double gcContent ); 

enum FileFormats { FASTA, PHYLIP };

class Sequence {
public:
    
    Sequence() : nts() {}     
    Sequence( const vector< int >& _nts ) : nts( _nts ) {}
    Sequence( const string& str );
    Sequence( const Sequence& other, int start, int length, bool strand = true ); 
    void copy( const Sequence& other ) { nts = other.nts; }
    Sequence( const Sequence& other ) { copy( other ); }
    
    
    Sequence& operator=( const Sequence& other ) { copy( other ); return *this; }

    
    bool operator==( const Sequence& other ) {
        if ( nts == other.nts ) return true;
        else return false;
    }
    
    
    int size() const { return nts.size(); }
    const int& operator[]( int pos ) const {
        assert( pos >= 0 && pos < size() );
        return nts[ pos ];	
    }
    int& operator[]( int pos ) {
        assert( pos >= 0 && pos < size() );
        return nts[ pos ];	
    }
    
    
    int push_back( int nt );
    int push_back( const Sequence& elem );	
    
    
    Sequence compRevCompl() const;
    
    
    void getNtCounts( vector< int >& counts ) const;
    bool containsMissing() const;
            
    
    void clear() { nts.clear(); }
    
    
    int load( const string& file, string& name, int format = FASTA );
    int load( const string& file, int format = FASTA );	
    
    
    friend ostream& operator<<( ostream& os, const Sequence& seq );
private:
    vector< int > nts;		
};

int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format = FASTA );
int readSequences( const string& file, vector< Sequence >& seqs, int format = FASTA );

int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format = FASTA );
int writeSequences( const string& file, const vector< Sequence >& seqs, int format = FASTA );


const double PSEUDO_COUNT = 0.25;

Matrix compWtmx( const Matrix& countMatrix, double pseudoCount );
Matrix countmatrixS( const string& file ); 
class Motif {
public:	
    
    Motif() : pwm(), background() {}	
    Motif( const Matrix& _pwm, const vector< double >& _background ); 
    Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background );		
    void copy( const Motif& other ) { pwm = other.pwm; background = other.background; LLRMat = other.LLRMat; maxSite = other.maxSite; maxLLR = other.maxLLR; }
    Motif( const Motif& other ) { copy( other ); }


    
    Motif& operator=( const Motif& other ) { copy( other ); return *this; }
    
    
    int length() const { return pwm.nRows(); }
    const Matrix& getPwm() const { return pwm; } 
	  Matrix& getPwm2() { return pwm; } 
    const vector< double >& getBackground() const { return background; }
    const Matrix& getLLRMat() const { return LLRMat; }
    const Sequence& getMaxSite() const { return maxSite; }
    double getMaxLLR() const { return maxLLR; }
    
    
    double LLR( const Sequence& elem ) const; 
    
    
    double energy( const Sequence& elem ) const;

    
    void sample( const gsl_rng* rng, Sequence& elem, bool strand = true ) const;
                    
    
    int load( const string& file, const vector< double >& background, string& name );
    int load( const string& file, const vector< double >& background );
    
    
    friend ostream& operator<<( ostream& os, const Motif& motif );	
private:
    Matrix pwm;	
    vector< double > background;	
    Matrix LLRMat;	
    Sequence maxSite;	
    double maxLLR;	
            
    
    void init();	
};

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names );
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs );


class Site {
public:
    
    Site() : start( 0 ), strand( true ), factorIdx( 0 ), energy( 0 ), wtRatio( 1 ) {}
    Site( int _start, bool _strand, int _factorIdx ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( 0 ), wtRatio( 1 ) {}
    Site( int _start, bool _strand, int _factorIdx, double _energy ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( _energy ) { wtRatio = exp(- energy ); }	
    void copy( const Site& other ) { start = other.start; strand = other.strand; factorIdx = other.factorIdx; energy = other.energy; wtRatio = other.wtRatio; }
    Site( const Site& other ) { copy( other ); }
    
    
    Site& operator=( const Site& other ) { copy( other ); return *this; }	
            
    friend ostream& operator<<( ostream& os, const Site& site );	
    
    int start;		
    bool strand;	
    int factorIdx;	
    double energy;	
    double wtRatio;	
};

bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs );


void sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp );
typedef vector< Site > SiteVec;		

int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy = false );
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy = false );


bool siteSortPredicate(const Site& d1, const Site& d2);
int Random(int n);




enum ModelType {
    LOGISTIC,   
    DIRECT,     
    QUENCHING,      
    CHRMOD_UNLIMITED,     
    CHRMOD_LIMITED,        
    BINS
};

ModelType getModelOption( const string& modelOptionStr );
string getModelOptionStr( ModelType modelOption ); 

enum FactorIntType {
    BINARY,     
    GAUSSIAN,    
    BINSF
}; 

string getIntOptionStr( FactorIntType intOption );

enum ObjType {
    SSE,    
    CORR,   
    CROSS_CORR  
};

ObjType getObjOption( const string& objOptionStr );
string getObjOptionStr( ObjType objOption );

enum SearchType {
    UNCONSTRAINED,  
    CONSTRAINED     
};

string getSearchOptionStr( SearchType searchOption );

class FactorIntFunc {
public:
    
    virtual double compFactorInt( double normalInt, double dist, bool orientation ) const = 0;

    
    virtual double getMaxDist() const = 0;    
};

class FactorIntFuncBinsf : public FactorIntFunc {
public:
    
FactorIntFuncBinsf(double _distThr,  int _nbins, double _orientationEffect = 1.0 ) : distThr( _distThr ),  nbins(_nbins),orientationEffect( _orientationEffect ) {}

    
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    
    double getMaxDist() const {
        return distThr;
    } 
private:
    int nbins;
    double distThr;		
    double orientationEffect;	
};


class FactorIntFuncBinary : public FactorIntFunc {
public:
    
    FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

    
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    
    double getMaxDist() const {
        return distThr;
    } 
    
private:
   
    double distThr;		
    double orientationEffect;	
};

class FactorIntFuncGaussian : public FactorIntFunc {
public: 
    
    FactorIntFuncGaussian( double _distThr, double _sigma ) : distThr( _distThr ), sigma( _sigma ) {
        assert( distThr > 0 && sigma > 0 );
    }

    
    double compFactorInt( double normalInt, double dist, bool orientation ) const; 

    
    double getMaxDist() const {
        return distThr;
    } 
private: 
    double distThr;     
    double sigma;       
};

class FactorIntFuncGeometric : public FactorIntFunc {
public:
    
    FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : distThr( _distThr ), spacingEffect( _spacingEffect ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

    
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    
    double getMaxDist() const {
        return distThr;
    } 
private:
    double distThr;		
    double spacingEffect;		
    double orientationEffect;	
};


class ExprPar {
public:
    
    ExprPar() : factorIntMat(), theV(), theVr() {}  
    ExprPar( int _nFactors );		
  
ExprPar( int _nFactors, vector< vector< vector<double> > > _theV , vector< vector< vector<double> > > _theVr);
    ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, double _basalTxp );
    ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators );	
	ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators , bool a);	
    void copy( const ExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; theV = other.theV; theVr = other.theVr;  }
    ExprPar( const ExprPar& other ) { copy( other ); }
	
    
    ExprPar& operator=( const ExprPar& other ) { copy( other ); return *this; }	
	
    
    int nFactors() const { return maxBindingWts.size(); }
	
    
    void getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const; 
	   void getFreePars2( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) ;
void getFreePars3( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const;
    
    void print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const;
    
 		
    
    int load( const string& file, IntMatrix& _coopMat , const vector< bool >& _repIndicators ); 

    
    void adjust(); 
    
    
    vector< double > maxBindingWts;			
    Matrix factorIntMat; 		
    vector< double > txpEffects;    
    vector< double > repEffects;    
    double basalTxp;        
    vector< vector< vector<double > > > theV ;
 vector< vector< vector<double > > > theVr ;
vector< vector< vector<double > > > getV()
{ return theV; }

    static ModelType modelOption;     
    static SearchType searchOption;    
    static int estBindingOption;    
    static int nbins;
    static double default_weight;	
    static double default_interaction;		
    static double default_effect_Logistic;   
    static double default_effect_Thermo;     
    static double default_repression;   
    static double default_basal_Logistic;       
    static double default_basal_Thermo;         
    static double min_weight;		
    static double max_weight;		
    static double min_interaction;	    
    static double max_interaction;	    
    static double min_interactionr;	    
    static double max_interactionr;	    
    static double min_effect_Logistic;   
    static double max_effect_Logistic;   
    static double min_effect_Thermo;    
    static double max_effect_Thermo;   
    static double min_repression;   
    static double max_repression;   
    static double min_basal_Logistic;    
    static double max_basal_Logistic;    
    static double min_basal_Thermo;   
    static double max_basal_Thermo;   
    static double delta;        
};


class ExprFunc {
public:
    
 

  
ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par,  const vector< int >& _B ,  const vector< int >& _Br);


void predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs, vector< double >& fOcc );

double predictZ( const SiteVec& _sites, const vector< double >& factorConcs );
    
    const vector< Motif >& getMotifs() const {
        return motifs;
    }
  
    
    double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs );
	double predictExpr2( const SiteVec& _sites, int length, const vector< double >& factorConcs );
double predictExpr3( const SiteVec& _sites, int length, const vector< double >& factorConcs );
   static vector<int > Bc;  
    static ModelType modelOption;     
void setsites( SiteVec& s) { sites = s ;}
private:
    
    const vector< Motif >& motifs; 
	 
  
   
    
    const FactorIntFunc* intFunc;   
    const vector< bool >& actIndicators;    
    int maxContact;     
    const vector< bool >& repIndicators;    
    const IntMatrix& repressionMat;    
    double repressionDistThr;   
   const vector<int >& B;     
const vector<int >& Br;
	
	vector< double > Z;
    vector< double > Zt;
				
	
	double compPartFunc();
vector< double > factorOcc;
    
    const ExprPar& par;
		    
    
    SiteVec sites;
    vector< int > boundaries;   

    
    vector< double > bindingWts; 
		
   
   
    
    
    double compFactorInt( const Site& a, const Site& b ) const;

    
    bool testRepression( const Site& a, const Site& b ) const;
};
class SeqAnnotator   {
public:
    
 SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs , Motif& _d2) : motifs( _motifs ), energyThrs( _energyThrs ), d2( _d2), B(),  seqSitesm1delete1(  ) { assert( motifs.size() == energyThrs.size() ); }

    SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs ) : motifs( _motifs ), energyThrs( _energyThrs ), B(),  seqSitesm1delete1(  ) { assert( motifs.size() == energyThrs.size() ); }
         SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs ,const  vector<int> _B ) : motifs( _motifs ), energyThrs( _energyThrs ), B(_B),seqSitesm1delete1( ) { assert( motifs.size() == energyThrs.size() ); }   
    
    int annot( const Sequence& seq, SiteVec& sites ) const;	

    int annoty( const Sequence& seq, SiteVec& sites, Matrix f, Matrix e, vector< double > p) ;

int annoty2( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; 
int annoty3( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites) ;
int annotyd( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites, vector< SiteVec >& dd) ;
int annoty4( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; 
int annotydorsal(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; 
int annotydorsalold(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ;
ExprPar SeqPar;
int maxZ( vector< double >& Z );
double maxZs( vector< double >& Z );
Site siteMax( vector< Site >& sites);
bool DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore);
void sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp );
    
    int compEnergy( const Sequence& seq, SiteVec& sites ) const;
int Delete1( vector< Site > &t ,vector< vector< Site > >& tt);
Motif d2;
private:
    vector< Motif > motifs;	
    vector< double > energyThrs;	
 const  vector<int> B;
vector< vector< Site > > seqSitesm1delete1;
};

class ExprPredictor {
public:
    
    ExprPredictor( const vector< vector< SiteVec > >& _seqSitesb, const vector< vector< double > >& _bindingData,  vector< SiteVec > & _seqSites, const vector< int >& _seqLengths,  Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const vector<int>& _binBord, const vector<int>& _binBordr , const vector < bool >& _indicator_bool, SeqAnnotator _anny, const string& _file ,  vector< SiteVec > & _seqSitesbot ,  vector< SiteVec >& _seqSitesm1,  vector< SiteVec >& _seqSitesm2, vector< SiteVec >& _seqSitesf2 ,vector< SiteVec >& _seqSitesbotf2, vector< SiteVec >& _seqSitesm1f2, vector< SiteVec >& _seqSitesm2f2, vector< SiteVec >& _seqSitesf3 ,  vector< SiteVec >& _seqSitesbotf3,vector< SiteVec >& _seqSitesm1f3 , vector< SiteVec >& _seqSitesm2f3);

    
    int nSeqs() const {
        return seqSites.size();
    }
	 const IntMatrix getcoopmat() {return coopMat; }
	const vector< bool > getactIndicators()  {return actIndicators; }
	const vector< bool > getrepIndicators()  {return repIndicators; }
	
	int nSeqsb() const {
        return seqSitesb.size();
    }
    int nFactors() const { 
        return motifs.size(); 
    }
    int nConds() const {
         return exprData.nCols();
    }
    const IntMatrix& getCoopMat() const {
        return coopMat;
    }
    const vector< bool >& getActIndicators() const {
        return actIndicators;
    }
    const vector< bool >& getRepIndicators() const {
        return repIndicators;
    }
    const IntMatrix& getRepressionMat() const {  
        return repressionMat;
    }
    const ExprPar& getPar() const { return par_model; }
 ExprPar& getPar2() { return par_model; }
    double getObj() const { return obj_model; }
	  void modifyIndicatorbool( vector< bool >&t ) { 
		indicator_bool=t;
	 } 
vector< bool > getIndicatorbool(  ) { 
		return indicator_bool;
	 } 
	 vector< double> getfixpars() {
         return fix_pars;
    }
	
        
    
    void printPar2(  );
const string& file;

    
   double compf( const ExprPar& par ,  vector< double >& f ) ;
    double compRMSE2( const ExprPar& par ) ;
   double compRMSE3( const ExprPar& par , int i) ;
	   double compRMSE4( const ExprPar& par , int i, int jc) ;
    double objFunc( const ExprPar& par ) const;
    double objFunc2(  ExprPar& par ) ;
	double objFuncborder(  ExprPar& par ) ;
double objFuncborder2(  ExprPar& par );
    
    int train( const ExprPar& par_init ); 	
    int train( const ExprPar& par_init, const gsl_rng* rng );   
    int train();	
	int train4( const ExprPar& par_init );
 int train3( const ExprPar& par_init, const gsl_rng* rng );   
      SeqAnnotator anny;      
    
    int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs ) ; 
 void predict2( vector< vector< double > >& occs ) ; 
    
void compOccMat(  const gsl_vector* x, void* params ) ;
void printFile( ostream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile2( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile3( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) ;
void printFile3b( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) ;
void printFile4( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) ;
void printFile5( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) ;
void printFiled( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFilePar_KfoldCV( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;

int createccdata();
int createccdata2();
int load( const string& fila );  
 
    
    static ModelType modelOption;     
    static int estBindingOption;    
    static ObjType objOption;       
static vector< Sequence > seqsy;
static vector< string > seqNmes;
 vector< int > cell ;
    
    static double exprSimCrossCorr( const vector< double >& x, const vector< double >& y ); 
    static int maxShift;    
    static double shiftPenalty;     

    
    static int nAlternations;   
    static int nRandStarts;     
    static double min_delta_f_SSE;      
    static double min_delta_f_Corr;     
    static double min_delta_f_CrossCorr;    
    static int nSimplexIters;       
    static int nGradientIters;      
   
void printBc( ) const;    
const vector< bool >& repIndicators;    
    
    void printPar( const ExprPar& par ) const;
 vector < bool > indicator_bool;
   vector < double > fix_pars;
    vector < double > free_pars;
 vector< SiteVec >& seqSites;		
 vector< SiteVec >& seqSitesbot;
vector< SiteVec >& seqSitesm1;
vector< SiteVec >& seqSitesm2;

 vector< SiteVec >& seqSitesf2;		
 vector< SiteVec >& seqSitesbotf2;
vector< SiteVec >& seqSitesm1f2;
vector< SiteVec >& seqSitesm2f2;
 vector< SiteVec >& seqSitesf3;		
 vector< SiteVec >& seqSitesbotf3;
vector< SiteVec >& seqSitesm1f3;
vector< SiteVec >& seqSitesm2f3;

int clasvar;
double spaceweights;
double getdof(  ){  return spaceweights;}
vector< vector< SiteVec > >AllData;  
vector< vector< int > >AllBorders;   
vector< double > Jacobian;
vector< vector< double > > Hessian;
static Matrix exprData2;		
    ExprFunc* createExprFunc2(  ExprPar& par ) ;
    static Matrix factorExprData2;
void compvar(vector< double >& vars);
vector< vector< Site > > seqSitesm1d1;
 vector< vector< vector< Site > > > d;
private:
 
    
    const vector< int >& seqLengths;           
  Matrix& exprData;		
    const vector< Motif >& motifs;		
   const Matrix& factorExprData;		
    const vector< vector< SiteVec > >& seqSitesb;		
   
  
const vector< vector< double > >& bindingData;		
    
    const FactorIntFunc* intFunc;   
    const IntMatrix& coopMat;       
    const vector< bool >& actIndicators;   
    int maxContact;     
    
    
    const IntMatrix& repressionMat;    
    double repressionDistThr;   
    
    
    ExprPar par_model;
    double obj_model;	
    double obj2;
  const  vector<int>& binBord;
const  vector<int>& binBordr;
   
    
    int randSamplePar( const gsl_rng* rng, ExprPar& par ) const; 

    
    bool testPar( const ExprPar& par ) const; 

  
    ExprFunc* createExprFunc( const ExprPar& par ) const;
 
     

    
    double compRMSE( const ExprPar& par ) const;		
    double compAvgCorr( const ExprPar& par ) const;     
 double compAvgCorr2( ExprPar& par ) ;
double compAvgCorrborder( ExprPar& par ) ;
    double compAvgCrossCorr( const ExprPar& par ) const;    
double compAvgCorrborder2(  ExprPar& par ) ;
double compAvgCorrborder8(  ExprPar& par ); 
    
    int simplex_minimize( ExprPar& par_result, double& obj_result ) ;
    int gradient_minimize( ExprPar& par_result, double& obj_result ) ;
    int gradient_minimize2( ExprPar& par_result, double& obj_result );
};

double gsl_obj_f( const gsl_vector* v, void* params );
double gsl_obj_f3( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad ); 
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad ); 
int gsl_obj_f3_fit( const gsl_vector* v, void* params,  gsl_vector * f );
int gsl_obj_df_fit( const gsl_vector* v, void* params, gsl_matrix * J); 
int gsl_obj_fdf_fit( const gsl_vector* v, void* params,  gsl_vector * f , gsl_matrix * J); 


#endif
