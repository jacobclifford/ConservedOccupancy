
#include "ExprPredictor.h"

//vector<int>ExprFunc::Bc(5,1);
int main( int argc, char* argv[] ) 
{
  
    // command line processing
    string seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, seqFileb, bIFile, dcFile, duFile, adamiFile, exprFile2;
    string outFile;     // output file
    double coopDistThr = 150;
    double factorIntSigma = 50.0;   // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 150;
    double energyThr = 5;
    int maxContact = 1;
bool free_fix_indicators[] = {0,0,0,
			    0,0,0,0,0,
///*
0,0,0,0,0,
0,0,0,0,0,
0,0,0 };
//	*/	//	1,1,1,1,1,
		//		1,1,1,1,1,
		//		1,1,1,1,1 };
	int nExps = 1;	// number of experiments
	vector <bool> indicator_bool ( free_fix_indicators, free_fix_indicators + sizeof( free_fix_indicators )/sizeof( bool ));
//cout << " indicator_bool " << indicator_bool[0] << endl;
//out << " size of indciator bool " << indicator_bool.size() << endl;
//cout << " sizeof(free_fix_indicators " << sizeof( free_fix_indicators )  << '\t' <<" sizeof(bool) " << sizeof( bool )   << endl;
//cout << "free_fix_indicators  = "<< free_fix_indicators[0] << endl;
    int binwidth; // = 30;
   // int nbins =  3  ;//int(coopDistThr/binwidth);  // this is different than the nbins private of ExprPar, or ExprPredictor
    ExprPredictor::nAlternations = 1;
	int nalt =1;
	int nrand2=1;
  //  vector< vector< vector<double> > >mm(3,vector< vector<double> >(3,vector<double>(3,1)));
    for ( int i = 1; i < argc; i++ ) {
	
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
	else if ( !strcmp( "-nrand", argv[ i ] ) )
            nrand2 = atoi( argv[++i] );
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];
	else if ( !strcmp( "-du", argv[ i ] ) )
            duFile = argv[ ++i ];  
	else if ( !strcmp( "-dc", argv[ i ] ) )
            dcFile = argv[ ++i ];
	else if ( !strcmp( "-sa", argv[ i ] ) )
            adamiFile = argv[ ++i ];            
        else if ( !strcmp( "-e", argv[ i ] ) )
            exprFile = argv[ ++i ];   
	else if ( !strcmp( "-e2", argv[ i ] ) )
            exprFile2 = argv[ ++i ];          
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
        else if ( !strcmp( "-f", argv[ i ] ) )
            factorExprFile = argv[ ++i ];    
        else if ( !strcmp( "-o", argv[ i ] ) )
            ExprPredictor::modelOption = getModelOption( argv[++i] );
        else if ( !strcmp( "-c", argv[ i ] ) )
            coopFile = argv[ ++i ];
        else if ( !strcmp( "-i", argv[ i ] ) )
            factorInfoFile = argv[ ++i ];            
        else if ( !strcmp( "-r", argv[ i ] ) )
            repressionFile = argv[ ++i ];  
        else if ( !strcmp( "-oo", argv[ i ] ) )
            ExprPredictor::objOption = getObjOption( argv[++i] );    
        else if ( !strcmp( "-mc", argv[i] ) )
            maxContact = atoi( argv[++i] );
        else if ( !strcmp( "-fo", argv[i] ) )
            outFile = argv[++i];    
        else if ( !strcmp( "-p", argv[i] ) )
            parFile = argv[++i]; 
        else if ( !strcmp( "-ct", argv[i] ) )
            coopDistThr = atof( argv[++i] ); 
else if ( !strcmp( "-binwt", argv[i] ) )
            binwidth = atof( argv[++i] ); 
        else if ( !strcmp( "-sigma", argv[i] ) )
            factorIntSigma = atof( argv[++i] );    
        else if ( !strcmp( "-rt", argv[i] ) )
            repressionDistThr = atof( argv[++i] );    
        else if ( !strcmp( "-et", argv[i] ) )
            energyThr = atof( argv[++i] );
        else if ( !strcmp( "-na", argv[i] ) )
            nalt = atoi( argv[++i] );  
  	else if ( !strcmp( "-n", argv[ i ] ) )
	    nExps = atoi( argv[ ++i ] );
	else if ( !strcmp( "-bs", argv[i] ) )
            seqFileb = argv[++i];    
        else if ( !strcmp( "-bi", argv[i] ) )
            bIFile = argv[++i]; 
    }
    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -et energyThr -na nAlternations -ct coopDistThr -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }           
vector< int > mmm(4);
mmm[0] = 30;
mmm[1] = 52;
mmm[2] = 155;
mmm[3] = 160;

vector< int > mmmr(4);
mmmr[0] = 50;
mmmr[1] = 100;
mmmr[2] = 150;
mmmr[3] = 200;
//mmm[4] = 50;


//     bool readSites = false;     // whether read sites (if true) or read sequences 
  // nbins = int(coopDistThr/binwidth); 
    // additional control parameters
    double gcContent = 0.4;
    FactorIntType intOption = BINSF;     // type of interaction function
    ExprPar::searchOption =  CONSTRAINED;      // search option: unconstrained; constrained. 
    ExprPar::estBindingOption = 1;
    int nbin = mmm.size() + 1;

//cout << "nbin = " << nbin << endl;
//cout << "mmm.size() "<< mmm.size() << endl;
    ExprPar::nbins = nbin;    // necessary in for loops for theV
   // ExprPredictor::nbins = nbins;
    //ExprFunc::coopDistThr=coopDistThr;
    //ExprFunc::nbins=nbins;
    //ExprPredictor::coopDistThr=coopDistThr;
     ExprPredictor::nAlternations=1;
    ExprPredictor::nRandStarts = 1;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 3;
    ExprPredictor::nGradientIters = 10;

    int rval;
    vector< vector< double > > data;    // buffer for reading matrix data
    vector< string > labels;    // buffer for reading the labels of matrix data
    string factor1, factor2;    // buffer for reading factor pairs

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
//////////////////////////////////////////////////7292011

rval = readSequences( seqFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );
/////////////////////////////////////////////////

    // read the expression data
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
    Matrix exprData( data ); 
    int nConds = exprData.nCols();
    //////////////////////////////////////////////73011 (this needs to be checked for shallow copy complications.
ExprPredictor::exprData2 = exprData;
//////////////////////////////////////////////////////////////
    // read the motifs

    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    rval = readMotifs( motifFile, background, motifs, motifNames ); 
    assert( rval != RET_ERROR );
    int nFactors = motifs.size();
//vector< vector< vector<double> > >mm(nFactors,vector< vector<double> >(nFactors,vector<double>(nFactors,1)));
///////////

Motif d;


d.copy(motifs[0]);
Matrix dtemp=d.getPwm2();

vector< double > e(4);
e[0] = 1;
e[1] = 1;
e[2] = 1;
e[3] = 10;
vector< vector< double > > d2s;
for(int i = 0 ; i < dtemp.nRows(); i++ ) {

	d2s.push_back( dtemp.getRow(i) ) ;

	//d2s.setRow( i,  dtemp.getRow(i) );
	if ( i == 4 ){
	d2s.push_back( e );

	//d2s.setRow( i,  e );
	}
}
Matrix d2sm( d2s );

//Motif d2( d2sm , background) ;

Matrix dcm =countmatrixS( dcFile );
Motif dcmm( dcm, background);


Matrix dum = countmatrixS( duFile );
Motif d2( dum, background);




////////////
    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
     
    // read the factor expression data
    labels.clear();
    data.clear();
    rval = readMatrix( factorExprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nFactors && condNames.size() == nConds );
    for ( int i = 0; i < nFactors; i++ ) assert( labels[i] == motifNames[i] );
    Matrix factorExprData( data ); 
    assert( factorExprData.nCols() == nConds ); 
    //////////////////////////////////////////////73011
ExprPredictor::factorExprData2 = factorExprData;
//////////////////////////////////////////////////////////////
    // site representation of the sequences
    vector< double > energyThrs( nFactors, energyThr ); 
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrs );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.annot( seqs[ i ], seqSites[ i ] );
            seqLengths[i] = seqs[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }
 vector< SiteVec > seqSitesbot( nSeqs );
seqSitesbot = seqSites;
 vector< SiteVec > seqSitesm1( nSeqs );
seqSitesm1 = seqSites;
 vector< SiteVec > seqSitesm2( nSeqs );
seqSitesm2 = seqSites;

cout << " asd" <<endl;
vector< SiteVec > seqSitesf2( nSeqs );
seqSitesf2 = seqSites;
vector< SiteVec > seqSitesbotf2(nSeqs);
seqSitesbotf2 = seqSites;
vector< SiteVec > seqSitesm1f2(nSeqs);
seqSitesm1f2 = seqSites; 
vector< SiteVec > seqSitesm2f2( nSeqs);
seqSitesm2f2 = seqSites;

vector< SiteVec >  seqSitesf3( nSeqs);
seqSitesf3 = seqSites;
vector< SiteVec >  seqSitesbotf3( nSeqs);
seqSitesbotf3 = seqSites;
vector< SiteVec > seqSitesm1f3( nSeqs );
seqSitesm1f3 = seqSites;
vector< SiteVec >  seqSitesm2f3( nSeqs);
seqSitesm2f3=seqSites;

vector< vector< Site > > seqSitesm1delete1;
seqSitesm1delete1 = seqSites;

///////////////////////////////////////////////////////////////////////////////////
SeqAnnotator anny( motifs, energyThrs , d2);
   
        
////////////////////////////////////////////////////////////////////////////////////
    // read the cooperativity matrix 
    IntMatrix coopMat( nFactors, nFactors, false );
    if ( !coopFile.empty() ) {
        ifstream fcoop( coopFile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the cooperativity file " << coopFile << endl;
            exit( 1 );
        }  
        while ( fcoop >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            coopMat( idx1, idx2 ) = true;
            coopMat( idx2, idx1 ) = true;
        }        
    } 

    // read the roles of factors
    vector< bool > actIndicators( nFactors, false );
    vector< bool > repIndicators( nFactors, false );
    if ( !factorInfoFile.empty() ) {
        ifstream finfo( factorInfoFile.c_str() );
        if ( !finfo ) {
            cerr << "Cannot open the factor information file " << factorInfoFile << endl;
            exit( 1 );
        }      
        string name;
        int i = 0, actRole, repRole;
        while ( finfo >> name >> actRole >> repRole ) {
            assert( name == motifNames[i] );
            if ( actRole ) actIndicators[i] = true;
            if ( repRole ) repIndicators[i] = true;
            i++;
        }
    }
    
    // read the repression matrix 
    IntMatrix repressionMat( nFactors, nFactors, false );
    if ( !repressionFile.empty() ) {
        ifstream frepr( repressionFile.c_str() );
        if ( !frepr ) {
            cerr << "Cannot open the repression file " << repressionFile << endl;
            exit( 1 );
        }        
        while ( frepr >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            repressionMat( idx1, idx2 ) = true;
        }        
    }
/*
/////////////////////////////////////////////////////////////////////////////////////////  ////////////  added data
 vector< vector< double > > datab;    // buffer for reading matrix data
    vector< string > labelsb;    // buffer for reading the labels of matrix data

 // read the sequences
    vector< Sequence > seqsb;
    vector< string > seqNamesb;
    rval = readSequences( seqFileb, seqsb, seqNamesb );
    assert( rval != RET_ERROR );
    int nSeqsb = seqsb.size();
//int nSeqsPerTF = allSeqs.size() / nExps;
    // read the chip data
    vector< string > condNamesb;  
    rval = readMatrix( bIFile, labelsb, condNamesb, datab );
    assert( rval != RET_ERROR );
    assert( labelsb.size() == nSeqsb );
    for ( int i = 0; i < nSeqsb; i++ ) assert( labelsb[i] == seqNamesb[i] );
    Matrix bData( datab ); 
    int nCondsb = bData.nCols();
    //////////////  11/10

 // site representation of the sequences
    vector< double > energyThrsb( nFactors, energyThr ); 
    vector< SiteVec > seqSitesb( nSeqsb );
    SeqAnnotator annb( motifs, energyThrsb );
 
        for ( int i = 0; i < nSeqsb; i++ ) {
            annb.annot( seqsb[ i ], seqSitesb[ i ] );
            
        }
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
*/
	// read the sequences
	vector< Sequence > allSeqs;
	vector< string > allNames;
	rval = readSequences( seqFileb, allSeqs, allNames );
    assert( rval != RET_ERROR );

	int nSeqsb = allSeqs.size() / nExps;
	vector< vector< Sequence > > seqs2( nExps );
	vector< vector< string > > names( nExps );	
	int counter = 0;
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
			seqs2[ i ].push_back( allSeqs[ counter ] );	
			names[ i ].push_back( allNames[ counter ] );
			counter++;
		}
	}

	// read the binding data: one row per experiment (binding of all sequences in that experiment)
	ifstream fdata( bIFile.c_str() );
    if ( !fdata ) {
        cerr << "Cannot find the binding data file " << bIFile << endl;
        exit( 1 );
    }
	vector< vector< double > > bindingData( nExps, vector<double>(nSeqsb) );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
			string name;
			fdata >> name;			
			if( name != names[ i ][ j ] ) { 
				cerr << "Error: " << names[ i ][ j ] << "\t" << name << endl;
				exit( 1 );
			}
			fdata >> bindingData[ i ][ j ];	
 
		}
	}	


vector< vector< SiteVec > > seqSitesb( seqs2.size() );  // shouldn't this be seqs2.size(), changed to 2 on 1129
	SeqAnnotator annb( motifs, energyThrs );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
// 			cout << "Exp " << i << "\tSeq " << j << endl;
			SiteVec sites;
			annb.annot( seqs2[ i ][ j ], sites );  // annot( sequence, sitevector)  // annot takes empty sitevector
			seqSitesb[ i ].push_back( sites );
			
		}
	}	
////////////////
  
    // create the expression predictor
    FactorIntFunc* intFunc; 
    if ( intOption == BINARY ) intFunc = new FactorIntFuncBinary( coopDistThr ); 
    else if ( intOption == GAUSSIAN ) intFunc = new FactorIntFuncGaussian( coopDistThr, factorIntSigma );
    else if ( intOption == BINSF ) intFunc = new FactorIntFuncBinsf( coopDistThr, nbin );
    else {
        cerr << "Interaction Function is invalid " << endl; exit( 1 ); 
    }

   ExprPredictor* predictor = new ExprPredictor(seqSitesb, bindingData, seqSites, seqLengths, exprData, motifs, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, mmm, mmmr, indicator_bool, anny, exprFile, seqSitesbot, seqSitesm1,seqSitesm2, seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, seqSitesbotf3,seqSitesm1f3, seqSitesm2f3 );  //520

    // read the initial parameter values
    ExprPar par_init( nFactors);  // 519
	if ( !parFile.empty() ) {
        rval = par_init.load( parFile, coopMat, repIndicators );
        if ( rval == RET_ERROR ) {
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        } 
    	}


	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)
 
 //   cout << "training values of parameters on variable sites" << endl;
(*predictor).clasvar = 1;
predictor->objFuncborder(  par_init );
//double obj = predictor->getObj();
//cout << " trained all objective = " << obj << endl;
//predictor->train( par_init, rng);
 double obj = predictor->getObj();
cout << " trained all objective  for clasvar 1= " << obj << endl;
predictor->printPar(predictor->getPar());
(*predictor).clasvar = 0;
   ExprPar::searchOption =  UNCONSTRAINED;
//ofstream to("ot.txt");
//predictor->printFile3(to,predictor->getPar(),*predictor);
//to.close();
   
//(*predictor).clasvar = 1;
cout << "finished first phase, training values of parameters on fixed sites" << endl;
bool free_fix_indicators2[] = {1,0,0,
			    1,0,0,0,0,
///*
0,0,0,0,0,
0,0,0,0,0,
0,0,0 };
vector <bool> indicator_bool2 ( free_fix_indicators2, free_fix_indicators2 + sizeof( free_fix_indicators2 )/sizeof( bool ));
predictor->modifyIndicatorbool( indicator_bool2);
//ExprPar::searchOption =  CONSTRAINED;
   ExprPar::searchOption =  UNCONSTRAINED;
 ExprPredictor::nRandStarts = nrand2;
 ExprPredictor::nAlternations=nalt;
    ExprPredictor::nSimplexIters = 20;
    ExprPredictor::nGradientIters = 300;
cout << " before train " << endl;
predictor->printPar(predictor->getPar());
predictor->train( predictor->getPar(), rng);
predictor->train4( predictor->getPar() );
cout << " after train " << endl;
predictor->printPar(predictor->getPar());
 //obj = predictor->getObj();
 obj = predictor->getObj();
cout << " trained all objective for clasvar0= " << obj << endl;
//predictor->objFuncborder( par_init  ); //predictor->getPar2() );
   ExprPar::searchOption =  UNCONSTRAINED;
ofstream to2("ot2.txt");
//predictor->printFile3b(to2,predictor->getPar(),*predictor);
to2.close();
ofstream to3("ot3.txt");
//predictor->printFile5(to3,predictor->getPar(),*predictor);
to3.close();
cout << " trained all objective = " << obj << endl;
cout << " 2*.1 = " << 2*.1 << endl;
cout << " degrees of freedom for the network " << predictor->getdof( ) << endl;
predictor->printPar(predictor->getPar());
ofstream fo("pars.txt");
//predictor->printFilePar_KfoldCV(fo, predictor->getPar(),*predictor);
fo.close();
   


/// use creatccdata2 with seq2expr using the same data (createcc asks for special files which have to be aligned with seq2expr) 9 13 11
//cout << "ccdata : " << endl;
//predictor->createccdata2();

///////////928
///*
///*
///*
string parFile2="pars2.txt";

seqs.clear();
seqNames.clear();
  
bool free_fix_indicators3[] = {1,0,0,
			    1,0,0,0,0,
///*
0,0,0,0,0,
0,0,0,0,0,
0,0,0 };
vector <bool> indicator_bool3 ( free_fix_indicators3, free_fix_indicators3 + sizeof( free_fix_indicators3 )/sizeof( bool ));
rval = readSequences(adamiFile, seqs, seqNames );
    assert( rval != RET_ERROR );
     nSeqs = seqs.size();
cout << " ok inside predictor2" << endl;
    // read the expression data
  data.clear();
  labels.clear(); 
    rval = readMatrix( exprFile2, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
	   vector< vector< double > > data2;    // buffer for reading matrix data
    vector< string > labels2; 
	//~exprData;
cout << " ok2 " << endl;
	Matrix exprData2( data ); 
cout << " ok3" << endl;
   nConds = exprData2.nCols();
    //////////////////////////////////////////////73011 (this needs to be checked for shallow copy complications.
//ExprPredictor::exprData2 = exprData;
//////////////////////////////////////////////////////////////
  seqSites.clear();
	vector< SiteVec > seqSitestes( nSeqs );
    seqSites=seqSitestes;
	    vector< int > seqLengthstes( nSeqs );
 seqLengths.clear();
	seqLengths=seqLengthstes;
    //seqLengths( nSeqs );
   // SeqAnnotator ann( motifs, energyThrs );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.annot( seqs[ i ], seqSites[ i ] );
            seqLengths[i] = seqs[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }
cout << " passed annot " << endl;
 //vector< SiteVec > seqSitesbot( nSeqs );
seqSitesbot.clear(); //vector< SiteVec > seqSitesm1( nSeqs );
seqSitesm1.clear(); 
seqSitesm2.clear();
seqSitesf2.clear();
seqSitesbotf2.clear();
seqSitesm1f2.clear();
seqSitesm2f2.clear();
seqSitesf3.clear();
seqSitesbotf3.clear();
seqSitesm1f3.clear();
seqSitesm2f3.clear();
seqSitesm1delete1.clear();

seqSitesbot=seqSites;
 //vector< SiteVec > seqSitesm1( nSeqs );
seqSitesm1 = seqSites;
 //vector< SiteVec > seqSitesm2( nSeqs );
seqSitesm2 = seqSites;

cout << " asd" <<endl;

seqSitesf2 = seqSites;

seqSitesbotf2 = seqSites;

seqSitesm1f2 = seqSites; 

seqSitesm2f2 = seqSites;


seqSitesf3 = seqSites;

seqSitesbotf3 = seqSites;

seqSitesm1f3 = seqSites;

seqSitesm2f3=seqSites;

seqSitesm1delete1 = seqSites;






 ExprPredictor* predictor2 = new ExprPredictor(seqSitesb, bindingData, seqSites, seqLengths, exprData2, motifs, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, mmm, mmmr, indicator_bool3, anny, exprFile2, seqSitesbot, seqSitesm1,seqSitesm2, seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, seqSitesbotf3,seqSitesm1f3, seqSitesm2f3 );  //520

rval = readSequences( adamiFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );
/////////////////////////////////////////////////
cout << " about to obj func " << endl;
 
predictor2->objFuncborder(  par_init );
cout << " about to train " << endl;
predictor2->train( par_init, rng);

predictor2->train4( predictor2->getPar() );
cout << " about to print after train4 " << endl;

ofstream to22("ot22.txt");
//predictor->printFile4(to22,predictor2->getPar(),*predictor2);
to22.close();
 
//predictor->load( parFile2 );
//*/
//predictor->printPar(predictor->getPar());
//ofstream fo("pars.txt");
//predictor->printFilePar_KfoldCV(fo, predictor->getPar(),*predictor);
//fo.close();
//*/

ofstream to("ot.txt");
predictor->printFile3b(to,predictor2->getPar(),*predictor);
to.close();


//ExprFunc* func1 = predictor->createExprFunc2( par_init );
//for(int i =0; i < ExprPredictor::seqsy.size() ; i++ ){
 //predictor->anny.annoty2( ExprPredictor::seqsy[ 0], seqSites[ 0 ], factorExprData, exprData, *func1 );
//	}
/*
par_init.load( parFile2, coopMat, repIndicators );
cout << "training parfile2 " << endl;
predictor->train( par_init);
*/
//fo << seqSites[0] << endl;
//fo << ExprPredictor::seqsy[0] << endl;
/*
ofstream fo( "expressout.txt", ios::app );

   fo << "Rows\t" << condNames << endl; 
 //   for ( int i = 0; i < nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor->predict( seqSites[0], seqLengths[0], targetExprs );
        vector< double > observedExprs = exprData.getRow( 0 );
  
        fo << seqNames[0]<< "\t" <<observedExprs << endl;      // observations
   //     fout << seqNames[i];
	fo << seqNames[0] << "\t" << targetExprs << endl;
//    }
//for ( int k = 0; k < motifs.size(); k++ ) {

int j=0;
//int i = seqSites[0][0].start;
for( int i = 0; i < seqSites[0].size() ; i++ ) {
while ( j < ExprPredictor::seqsy[0].size() ) {
	if (seqSites[0].empty() ) { break; }
	for ( ;; ) {
		if ( j < seqSites[0][i].start ) {
			fo << "b";
			j++;
		}
		else {
			if( seqSites[0][i].factorIdx == 0 ) { fo << "d"; j++  ; break; }
			if( seqSites[0][i].factorIdx == 1 ) { fo << "t"; j++ ; break; }
			if( seqSites[0][i].factorIdx == 2 ) { fo << "s"; j++ ; break; }
		} // else
	}
}
}
fo << endl;
fo.close();
*/
ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }
    fout << "Rows\t" << condNames << endl;
    for ( int i = 0; i < nSeqs; i++ ) {
      vector< double > observedExprs = exprData2.getRow( i );
	vector< double > targetExprsm1;
	vector< double > dorsalExprs = factorExprData.getRow( 2 );
        predictor2->predict(seqSitesm1[ i ], seqLengths[i], targetExprsm1 );
       // vector< double > observedExprs = exprData.getRow( i );
        
        // error
       // double beta; 
       // double error = sqrt( least_square( targetExprs, observedExprs, beta ) / nConds );
	fout << seqNames[i] << '\t';
        // print the results
        fout << observedExprs << endl;      // observations
   //     fout << seqNames[i];
	fout << seqNames[i];
        for ( int j = 0; j < nConds; j++ ) {
	if(j==0){ fout << "\t" <<  targetExprsm1[j]+.0001; continue;} 
	 fout << "\t" <<  targetExprsm1[j];       // predictions
	}
/*	fout << "dl";

        for ( int j = 0; j < nConds; j++ ) {
	fout << "\t" <<  factorExprData(1,j);    // factor dorsla
	}
*/
	fout << endl;
	fout << "Dl" << '\t';
        // print the results
        fout << dorsalExprs << endl;      // observations
	


}

/*
ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }
    fout << "Rows" ;
 for ( int j = 0; j < nConds-25; j++ )  fout << "\t" << condNames[j] ;
fout << endl;
    for ( int i = 0; i < nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor->predict( seqSitesf2[i], seqLengths[i], targetExprs );
        vector< double > observedExprs = exprData.getRow( i );
        
        // error
       // double beta; 
       // double error = sqrt( least_square( targetExprs, observedExprs, beta ) / nConds );
	fout << seqNames[i] ;
        // print the results
       for ( int j = 0; j < nConds-25; j++ )  fout<< "\t" <<setprecision(3) << observedExprs[j] ;      // observations
   //     fout << seqNames[i];
fout << endl;
	fout << seqNames[i];
        for ( int j = 0; j < nConds-25; j++ ) {fout << "\t"    << setprecision(3) << targetExprs[j]; }      // predictions
      // fout << '\t' << observedExprs ;  
fout << endl;
fout << endl;
	}
*/
/*
vector< double > targetExprsm1;
        predictor->predict(seqSitesm1[ i ], seqLengths[i], targetExprsm1 );
       // vector< double > observedExprs = exprData.getRow( i );
        
        // error
       // double beta; 
       // double error = sqrt( least_square( targetExprs, observedExprs, beta ) / nConds );
	fout << seqNames[i] << '\t';
        // print the results
        fout << observedExprs << endl;      // observations
   //     fout << seqNames[i];
	fout << seqNames[i];
        for ( int j = 0; j < nConds; j++ ) fout << "\t" <<  targetExprsm1[j];       // predictions
      // fout << '\t' << observedExprs ;  
fout << endl;


}
*/
/*
ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }

   // fout << "Rows\t" << condNames << endl;
vector< double > ttargetExprs;
vector< double > oobservedExprs;
    for ( int i = 0; i < nSeqs; i++ ) {
        vector< double > targetExprs;
        predictor->predict( seqSites[i], seqLengths[i], targetExprs );   // predictor's par_model contains best parameters
        vector< double > observedExprs = exprData.getRow( i );
        oobservedExprs.push_back(observedExprs[0]);
	ttargetExprs.push_back(targetExprs[0]);
        // error
        double beta; 
        double error = sqrt( least_square( targetExprs, observedExprs, beta ) / nConds );
// fout <<  observedExprs[i] << "\t"<< targetExprs[i] << endl;// beta * targetExprs[j];       // predictions
      //  fout << endl;

        // print the results
       fout  << observedExprs << "\t" << targetExprs << endl;// beta * targetExprs[j];       // predictions
        //fout << endl;
//fout << endl;
//fout << seqNames[i];
        //for ( int j = 0; j < nConds; j++ )
 //fout << "\t" << factorExprData.getRow(0);// beta * targetExprs[j];       // predictions
   //     fout << endl;
 
 
      } 
    cout <<  " correlation  " << correlation( ttargetExprs, oobservedExprs ) << endl;
*/
/*
vector< double > var;
 predictor->compvar( var);
//cout << var << endl;
ofstream ffout( "var.txt" );
 for ( int i = 0; i < nSeqs; i++ ) {
 ffout  << observedExprs[i] << "\t" << targetExprs[i] << "\t" << var[i] << endl;
}
 
       // print the agreement bewtween predictions and observations
     //   cout << seqNames[i] << "\t"; 
       // if ( ExprPredictor::objOption == SSE ) 
            //cout << error << endl;
       // else if ( ExprPredictor::objOption == CORR ) 
    //        cout << corr( targetExprs, observedExprs ) << endl;
        //else if ( ExprPredictor::objOption == CROSS_CORR )
         //   cout << ExprPredictor::exprSimCrossCorr( observedExprs, targetExprs ) << endl; 
    }
*/
    return 0;	
}

