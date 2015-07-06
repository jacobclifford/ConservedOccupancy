#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace std;
double THR = 0.65;

double pgp (const vector <double>& x, const vector<double>& y)
{
/*
Input: ex_vec (of zeroes and ones), pr_vec

Calculation:

unexprsum=sum(1-ex_vec) # bins of nonexpression
pr_un_expr = sum(ex_vec*pr_vec)/sum(ex_vec) # average prediction under
expression
pr_un_nexpr = sum((1-ex_vec)*(pr_vec))/unexprsum # average prediction
under non-expression
npr_un_nexpr = sum((1-ex_vec)*(1-pr_vec))/unexprsum # average area
above prediction under non-expression

if(pr_un_expr>= pr_un_nexpr+0.05){
       ret_val=pr_un_expr/2 - pr_un_nexpr + npr_un_nexpr/2
}else{
       ret_val=-0.5; ## minimum values above case can take
}

return(ret_val)
*/
	double unexprsum = 0;
	for(int i = 0; i < x.size(); i++){
		unexprsum += (1 - x[i]);
	}
	
	double sum_1 = 0;
	double sum_2 = 0;
	for(int i = 0; i < x.size(); i++){
		sum_1 += (x[i] * y[i]);
		sum_2 += x[i];
	}
	double pr_un_expr = sum_1/sum_2;
	sum_1 = 0;
	for(int i = 0; i < x.size(); i++){
		sum_1 += ((1 - x[i]) * y[i]);
	}
	double pr_un_nexpr = sum_1/unexprsum;
	sum_1 = 0;
	for(int i = 0; i < x.size(); i++){
		sum_1 += ((1 - x[i]) * (1 - y[i]));
	}
	double npr_un_nexpr = sum_1/unexprsum;
	double ret_val;
	if(pr_un_expr >= pr_un_nexpr + 0.05){
		ret_val = pr_un_expr / 2.0 - pr_un_nexpr + npr_un_nexpr / 2.0;
	}
	else{
		ret_val = -0.5;
	}
	return ret_val;
}

double least_square( const vector< double >& x, const vector< double >& y, double& beta)
{
    assert( x.size() == y.size() );
    int n = x.size();

    double numerator = 0, denom = 0;
    for ( int i = 0; i < n; i++ ) {
        numerator += x[i] * y[i];
        denom += x[i] * x[i];
    }
    beta = numerator / denom;
	beta = 1;
    double rss = 0;
    for ( int i = 0; i < n; i++ ) {
        rss +=  abs( y[i] - beta * x[i] );
    }

    return rss;
}

double least_square2( const vector< double >& x, const vector< double >& y, double& beta)
{
    assert( x.size() == y.size() );
    int n = x.size();

    double numerator = 0, denom = 0;
    for ( int i = 0; i < n; i++ ) {
        numerator += x[i] * y[i];
        denom += x[i] * x[i];
    }
    beta = numerator / denom;
	beta = 1;
    double rss = 0;
    for ( int i = 0; i < n; i++ ) {
        rss +=  ( y[i] - beta * x[i] )*( y[i] - beta * x[i] );
    }

    return rss;
}


double compute_cc (vector <double> &observed, vector <double> &predicted)
{
	assert(observed.size() == predicted.size());
	int n = observed.size();
	
	double sumx = 0;
	double sumy = 0;
	double sumxy = 0;
	double sumx2 = 0;
	double sumy2 = 0;
		
	for(int j = 0; j < n; j++){
		sumx += observed[j];
		sumx2 += observed[j] * observed[j];
		sumy += predicted[j];
		sumy2 += predicted[j] * predicted[j];
		sumxy += observed[j] * predicted[j];
	}
		
	double corr_coeff = (n * sumxy - sumx * sumy)/(sqrt(n * sumx2 - sumx * sumx) * sqrt(n * sumy2 - sumy * sumy));
	return corr_coeff;
}

string strip_space(string s)
{
	string t = "";
	int len = s.length();
	for(int i = 0; i < len; i++){
		if(!isspace(s[i])){
			t += s[i];
		}
	}
	return t;
}

int main(int argc, char* argv[])
{
/*
Usage:
name of the file containing data
*/
	
	ifstream info_file(argv[1]);  // DI.plot
	assert(info_file.is_open());
	
	vector <string> dat_files;
	vector <string> models;
	
	dat_files.clear();
	models.clear(); 
	
	while(!info_file.eof()){
		string s1, s2;
		info_file >> s1 >> s2;  // out.txt and Direct
		s1 = strip_space(s1);
		s2 = strip_space(s2);
		if( s1.length() == 0 || s2.length() == 0 ){
			break;
		}
		dat_files.push_back(s1);  // out.txt
		models.push_back(s2);     // Direct
	}
		
	ofstream plot_file (argv[2]);  // plotcmd file is created
	assert(plot_file.is_open());
	
	/*ofstream eps_file (argv[3]);
	assert(eps_file.is_open());*/
	cout << " heree " << endl;
	ofstream summary_file("summary.tex");
	assert(summary_file.is_open());	
// temp_file ( a )  is like an initization, so it's the same as temp_file = a.
	ifstream temp_file(dat_files[0].c_str());  // tmep_file is a ifstream object, this object points to dat_file0, which is Rows, it is pointer due to c_str()?
	string temp;
	getline(temp_file, temp);  // getline(stream, read) , where read is the variable that will be filled by the stream object, here holds observed data for first enhancer.
	temp_file.close();
	
	istringstream instr;
	instr.str(temp);  // initializes istringstream object to temp 
	instr >> temp;  // like cin, inputs characters from temp until whitespace is found
	assert (temp == "Rows");
	
	int start, end;
	instr >> start;
	while(!instr.eof()){   // for a string eof is newline?
		instr >> end;
	}
	
	int ncond = end - start + 1;
	
	int n_models = models.size();
	
	vector < double > sse (n_models, 0);
	
	vector < double > sum_cc (n_models, 0);
	
	vector < double > count_thr (n_models, 0);
	
	vector <double > sse3;
	vector <double > sse2;
	vector < double > cc_vector;
	vector < double > pgp_vector;
	
	vector < vector < double > > max_scale;
	max_scale.clear();
	
	cc_vector.clear();
	pgp_vector.clear();
	
	vector <string> seqs;
	seqs.clear();
	int nseq;
	for ( int index = 0; index < n_models; index++ ){
		ifstream input_file (dat_files[index].c_str());
		nseq = 0;
		string temp;
		getline( input_file, temp ); //reads the ROWS line
		max_scale.push_back( vector < double > () );
		while(!input_file.eof()){   // keeps looping  until out.txt is eof, so n_models can still be one.
			/*samee start*/
			//string temp;
			//getline( input_file, temp ); //reads the ROWS line
			/*samee end*/
			//cout << temp << endl;
			input_file >> temp;           // this holds the sequence label like rho2017
			//cout << temp << endl;
			string seq1 = strip_space( temp );
			if( seq1.length() == 0 ){
				break;
			}
			vector < double > obs;
			obs.clear();
			for( int i = 0; i < ncond; i++ ){
				double val;
				input_file >> val;     //  now everything after the label until nconds iterations of whitespace
				obs.push_back( val );
			}
			input_file >> temp;		// take label again, like rho2017
			//cout << temp << endl;
			string seq2 = strip_space( temp );
			assert( seq1 == seq2 );		// assert labels are the same
			seqs.push_back( seq1 );         // seqs holds the labels
			nseq ++;
			max_scale[ index ].push_back(0);
			vector <double> pre;
			pre.clear();
			for( int i = 0; i < ncond; i++ ){
				double val;
				input_file >> val;
				if( val > max_scale[ index ][ nseq - 1 ]){   //determining the range or scale for plotcmd file
					max_scale[ index ][ nseq - 1] = val;
				}
				pre.push_back(val);
			}
/////////////////////////////////////////////////////////////////dorsal
			input_file >> temp;           // this holds the sequence label like dl
			//cout << temp << endl;
			string seq3 = strip_space( temp );
			if( seq3.length() == 0 ){
				break;
			}
			//cout << "CC = " << cc << " for: " << seq1 << endl;
			vector <double> dl;
			dl.clear();
			for( int i = 0; i < ncond; i++ ){
				double val;
				input_file >> val;
				if( val > max_scale[ index ][ nseq - 1 ]){   //determining the range or scale for plotcmd file
					max_scale[ index ][ nseq - 1] = val;
				}
				dl.push_back(val);
			}
//////////////////////////////////////////////////////////////   snail
/*
			input_file >> temp;           // this holds the sequence label like sn
			//cout << temp << endl;
			string seq3 = strip_space( temp );
			if( seq3.length() == 0 ){
				break;
			}
			//cout << "CC = " << cc << " for: " << seq1 << endl;
			vector <double> snail;
			snail.clear();
			for( int i = 0; i < ncond; i++ ){
				double val;
				input_file >> val;
				if( val > max_scale[ index ][ nseq - 1 ]){   //determining the range or scale for plotcmd file
					max_scale[ index ][ nseq - 1] = val;
				}
				snail.push_back(val);
			}

*/
////////////////////////////////////////////////////////////// objectives
			double beta;
			double rss = least_square ( obs, pre, beta );
			sse [ index ] += rss;
			sse2.push_back( rss );

			double rss3 = least_square2 ( obs, pre, beta );
			sse3.push_back( rss3 );
		/*	
			for( int i = 0; i < pre.size(); i++ ){
				pre[ i ] *= beta;                      //  multiply by the amplitude// take this  out.
			}
		*/
			double cc = compute_cc ( obs, pre );
			if(cc > THR){
				count_thr[ index ] ++;
			}
			sum_cc [ index ] += cc;
			cc_vector.push_back( cc );
			
			
			double pgp_score = pgp ( obs, pre );
			pgp_vector.push_back( pgp_score );
			
			string output_file_name = seq1 + "_" + models[index] + ".tempdat";   //string concatenation op
			ofstream output_file ( output_file_name.c_str() );
			for( int i = 0; i < ncond; i++ ){
				output_file << i + start << "\t" << obs[ i ] << "\t" << pre[ i ] << "\t" << dl[i] <<endl;
			}
			output_file.close();
			
			
			temp = "";
			seq1 = "";
			seq2 = "";
			/*compute the p-value*/
			/*
			double cc_1_2 = compute_cc(pre_1, pre_2);
		
			double t = (cc_1 - cc_2) * sqrt(((ncond - 3)*(1 + cc_1_2))/(2 * (1 - cc_1 * cc_1 - cc_1_2 * cc_1_2 - cc_2 * cc_2 + 2 * cc_1 * cc_1_2 * cc_2)));
			*/

		}
		input_file.close();
			
		cout << nseq << " sequences read from " << dat_files[ index] << " for model: " << models[ index ] << endl;
		//exit(1);
	}

	vector < double > max_pre (nseq, 0);
	/*
	for( int index = 0; index < nseq; index ++ ){
		/*for( int i = 0; i < n_models; i++ ){
			string seq_dat_file_name = seqs[ index ] + "_" + models[ i ] + ".tempdat";
			ifstream seq_dat_file ( seq_dat_file_name.c_str() );
			double bin, obs_val, pre_val;
			for( int j = 0; j < ncond; j++ ){
				seq_dat_file >> bin >> obs_val >> pre_val;
				if( pre_val > max_pre[ index ] ){
					max_pre[ index ] = pre_val;
				}
			}
			seq_dat_file.close();
		}
		for( int i = 0; i < n_models; i++ ){
			string seq_dat_file_name = seqs[ index ] + "_" + models[ i ] + ".tempdat";
			ifstream seq_dat_file ( seq_dat_file_name.c_str() );
			ofstream seq_out_file ( (seqs[ index ] + "_" + models[ i ] + ".dat").c_str() );
			double bin, obs_val, pre_val;
			for( int j = 0; j < ncond; j++ ){
				seq_dat_file >> bin >> obs_val >> pre_val;
				seq_out_file << bin << "\t" << obs_val << "\t" << ( pre_val / max_scale[ i ][ index ] ) << endl;  // transform processed data to 0-1 range
			}
			seq_dat_file.close();
			seq_out_file.close();
		}
	}
	*/
	ifstream plot_priority_file ( argv[ 3 ] );
	vector < int > plot_priority;
	int plotable = 0;
	int first_plotable = -1;
	plot_priority.clear();
	for( int index = 0; index < nseq; index++ ){
		int plot_priority_val;
		plot_priority_file >> plot_priority_val;
		if( plot_priority_val ){
			plot_priority.push_back( plot_priority_val );
			plotable ++;
		}
		else{
			plot_priority.push_back( 0 );
		}
	}

	plot_priority_file.close();

	for( int priority = 1; priority <= 3; priority ++ ){
		for ( int index = 0; index < nseq; index ++ ){
			if ( plot_priority[ index ] == priority && first_plotable == -1 ){
				first_plotable = index;
				break;
			}
		}
		if( first_plotable > 0 ){
			break;
		}
	}


	int lmargin = 4;
	int bmargin = 3;
	int rmargin = 7;
	int tmargin = 2;


	double size1 = 1;
	double size2 = 1;
	
	double xmin = start;
	double xmax = end;
	
	double eps_x = 1;
	
	double ymin = 0;
	double ymax = 1;
	
	double eps_y = 0.025;
	
	string overall_font = "Helvetica";
	int overall_font_size = 30;
	
	string xlabel_font = "Helvetica";
	string ylabel_font = "Helvetica";
	int xlabel_font_size = 37;
	int ylabel_font_size = 37;
	
	string title_font = "Times-Italic";
	int title_font_size = 45;
	
	
	int col_count_in_fig = 3;
	
	
	int plotted = 0;


	//plot_file.precision( 3 );	
	plot_file.setf( ios::fixed, ios::floatfield );

	for( int priority = 1; priority <= 3; priority ++ ){
		for( int index = 0; index < nseq; index ++ ){
		
			if( plot_priority[ index ] == priority ){
			
				plot_file << "set term postscript eps enhanced color \"" << overall_font << "," << overall_font_size << "\"" << endl; 
				plot_file << "set lmargin " << lmargin << endl;
				plot_file << "set bmargin " << bmargin << endl;
				plot_file << "set rmargin " << rmargin << endl;
				plot_file << "set tmargin " << tmargin << endl;
				plot_file << "set size " << size1 << ", " << size2 << endl;
				plot_file << "set xrange [ " << start - eps_x << " : " << end + eps_x << " ]" << endl;
				plot_file << "set yrange [ " << ymin - eps_y << " : " << ymax + eps_y << " ]" << endl;
				
				string model_font_plot_color[] = { " rgb \"#FF0000\"",
									" rgb \"#228B22\"",
									" rgb \"#0000FF\""
									};

				
				
				int lw = 2;

				plot_file << "set style line 1 lt -1 lc " << "rgb \"red\" lw " << lw << endl;

				for ( int i = 1; i <= n_models; i++ ){
					plot_file << "set style line " << i + 1 << " lt -1 lc " << model_font_plot_color[ i ] << " lw " << lw << endl;
				}
	
				plot_file << "set output \"" << seqs[ index ] << ".eps\"" << endl;
				

				string title_font_color = " def";
				
				if( priority == 1 ){
					title_font_color = model_font_plot_color[ 1 ];
				}
				else if( priority == 3 ){
					title_font_color = model_font_plot_color[ 2 ];
				}

				plot_file << "set title \"" << seqs[ index ] << "\"" << " noenhanced font \"" << title_font << "," << title_font_size << "\" tc" << title_font_color << endl;
				
				string xlabel = ( plotted + col_count_in_fig >= plotable ) ? "DV z axis" : "";
				string ylabel = ( plotted % col_count_in_fig == 0 ) ? "Expr" : "";
				plot_file << "set xlabel \""<< xlabel << "\" font \"" << xlabel_font << "," << xlabel_font_size << "\""<< endl;
				plot_file << "set ylabel \""<< ylabel << "\" font \"" << ylabel_font << "," << ylabel_font_size << "\""<< endl;
			
			
				double legend_xloc = 1;
				double legend_yloc = ymax - 0.05;
				double legend_yloc_del = .1;

				for( int i = 0; i < n_models; i++ ){
					//plot_file << "set label \"" << setprecision( 3 ) <<  ceil( 1000 * cc_vector[ index + i * nseq ] )/1000  << "\" " << "at " << legend_xloc << ", " << legend_yloc << " tc " << model_font_plot_color[ i + 1] << endl;
					plot_file << "set label \""<< "CC=" << setprecision( 3 ) <<  ceil( 1000 * cc_vector[ index + i * nseq ] )/1000  << "\" " << "at " << legend_xloc << ", " << legend_yloc << " tc rgb \"black\"" << endl;
//					plot_file << "set label \""<< "AE=" << setprecision( 1 ) <<  sse2[index+i*nseq]  << "\" " << "at " << legend_xloc << ", " << legend_yloc-.2 << " tc rgb \"black\"" << endl;
plot_file << "set label \""<< "SE=" << setprecision( 1 ) <<  sse3[index+i*nseq]  << "\" " << "at " << legend_xloc << ", " << legend_yloc-.2 << " tc rgb \"black\"" << endl;
					//legend_yloc -= legend_yloc_del;	
			
				}		
			
			
			//plot_file << "set yrange [0:" << max_scale[ index ] << "]" << endl;
				int line_style = 1;
				if( index == first_plotable ){
					//plot_file << "set label \"Label\" at 75,0.8 textcolor lt 3" << endl;
					double legend_x = 10;
					plot_file << "set label \"observed.\" at " << legend_x << ", 0.95 tc rgb \"red\"" << endl;
					plot_file << "set label " << "\"" << "model" << "\""  << " at "<< legend_x  <<", 0.85 tc rgb \"#228B22\"" << endl;
					plot_file << "set label " << "\"" << "Dorsal"  << "\""  << " at "<< legend_x  <<", 0.75 tc rgb \"blue\"" << endl;
					//plot_file << "set label " << "\"" << "Snail"  << "\""  << " at "<< legend_x  <<", 0.65 tc rgb \"blue\"" << endl;
					if( n_models > 1 ){
						plot_file << "set label " << "\"" << models[ 1 ]  << "\""  << " at " << legend_x << ", 0.75 tc rgb \"#0000FF\"" << endl;
					}
					
				}
				
				for( int i = 0; i < n_models; i++){
					
					if( i == 0 ){
						plot_file << "plot \"" << seqs[ index ] + "_" + models[ 0 ] + ".tempdat" << "\" using 1:2 notitle with lines ls " << line_style ++;	
					}
					
					plot_file << ", \"" << seqs[ index ] + "_" + models[ i ] + ".tempdat" << "\" using 1:3 notitle with lines ls " << line_style ++;
					plot_file << ", \"" << seqs[ index ] + "_" + models[ i ] + ".tempdat" << "\" using 1:4 notitle with lines ls " << line_style ++;
					
					}	
				
				plot_file << endl;
				plot_file << "reset" << endl;		
				plotted ++;
		
			}
		}
	}

	plot_file.close();		


	//generate the format.tex file:

	ofstream tabular_format_file ( "format.tex" );
	
	plotted = 0;
	
	for ( int priority = 1; priority <= 3; priority ++ ){
		for ( int index = 0; index < nseq; index ++ ){
			if( plot_priority[ index ] == priority ){
				tabular_format_file << "\\resizebox{\\plotlen}{!}{\\includegraphics{" << seqs[ index ] << ".eps" << "}}";
				plotted ++;

				if ( ( plotted % col_count_in_fig ) == 0 ){
					tabular_format_file << " \\\\" << endl;
				}
				else{
					tabular_format_file << " &" << endl;
				}
			}
		}
	}
	
	assert ( plotted == plotable );
	for( ; ( plotted % col_count_in_fig ) != 0; ){
		plotted ++;
		if( plotted % col_count_in_fig == 0 )
			tabular_format_file << " \\\\" << endl;
		else{
			tabular_format_file << " &" << endl;
		}
	}
	tabular_format_file.close();

	cout << endl << endl;
	cout << "Generated Output for " << nseq << " Sequences at " << ncond << " conditions." << endl;
	cout << "Summary:" << endl;

	for( int index = 0; index < n_models; index ++ ){
		cout << models[ index ] << endl;
		summary_file << "\\textbf{" << models[ index ] << ":}" << endl << endl;
		double rmse = sqrt ( sse[ index ]/( ncond * nseq ) );
		double avg_cc = sum_cc[ index ]/nseq;
		cout << "Avg. CC: " << avg_cc << endl;
		cout << "RMSE: " << rmse << endl;
		cout << "Number of predictions above " << THR << ": " << count_thr[ index ] << endl;
	
		summary_file << "Avg. CC: " << avg_cc << endl << endl;
		summary_file << "RMSE: " << rmse << endl << endl;
		summary_file << "Number of predictions above " << THR << ": " << count_thr[ index ] << endl << endl;
	}

	cout << endl << endl;
	

	return 0;
}
