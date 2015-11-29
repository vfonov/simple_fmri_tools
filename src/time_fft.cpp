#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <vector>
#include <fftw3.h>
#include <getopt.h>
#include <complex>
#include <minc_io_exceptions.h>
#include <minc_1_simple.h>

using namespace  std;
using namespace  minc;
typedef std::complex<float> _complex;

void show_usage (const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <input> <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl
    << "--low <threshold> "<<endl
    << "--hi  <threshold> "<<endl
    << "--spectrum - calculate power spectrum instead of filtering"<<endl;
  
}

int main (int argc, char **argv)
{
  int verbose=1;
  int clobber=0;
  int low=0,hi=0;
  int spectrum=0;
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"spectrum", no_argument, &spectrum, 1},
    {"clobber", no_argument, &clobber, 1},
    {"low",  required_argument, 0, 'l'},
    {"hi",   required_argument, 0, 'h'},
		{0, 0, 0, 0}
		};
    
	int c;
	for (;;)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "l:h:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 'l':
			low = atoi(optarg);
			break;
		case 'h':
			hi = atoi(optarg);
			break;
		case 'v':
			cout << "Version: 1.0" << endl;
			return 0;
		case '?':
			/* getopt_long already printed an error message. */
		default:
			show_usage (argv[0]);
			return 1;
		}
	}
  
  if((argc - optind)<2)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input=argv[optind];
  std::string output=argv[optind+1];
  if (!clobber && !access (output.c_str(), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  
	try {
    minc_1_reader rdr;
    rdr.open(input.c_str());
    minc_1_writer wrt;
    wrt.open(output.c_str(),rdr.info(),3,rdr.datatype());

    if(rdr.info().size()==4 && rdr.info()[3].name=="time")
    {
      std::cout<<"Processing time dimension..."<<std::endl<<std::flush;
      fftwf_plan _plan,_iplan;
      //this is going to be ineffective way, i know!
      float *temp=new float[rdr.info()[3].length];
      unsigned int clength=rdr.info()[3].length;
      //while(clength<rdr.info[3].length) clength*=2;
      
      _complex* ctemp=new _complex[clength];
      _plan  = fftwf_plan_dft_1d(clength, reinterpret_cast<fftwf_complex*>(ctemp), reinterpret_cast<fftwf_complex*>(ctemp), FFTW_FORWARD,  FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
      _iplan = fftwf_plan_dft_1d(clength, reinterpret_cast<fftwf_complex*>(ctemp), reinterpret_cast<fftwf_complex*>(ctemp), FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
      
      minc_input_iterator<float> in(rdr);
      minc_output_iterator<float> out(wrt);
      in.begin();
      out.begin();
      
      for(int x=0;x<rdr.info()[0].length;x++)
        for(int y=0;y<rdr.info()[1].length;y++)
          for(int z=0;z<rdr.info()[2].length;z++)
          {
            //read the column
            for(int i=0;i<clength;i++)
            {
              ctemp[i]=_complex(0,0);
              if(i<rdr.info()[3].length)
              {
                ctemp[i]=_complex(in.value(),0);
                in.next();
              }
            }
            fftwf_execute(_plan);
            if(spectrum)
            { //calculate power spectrum
              for(int i=0;i<clength;i++)
                ctemp[i]=abs(ctemp[i])*abs(ctemp[i])/clength;
                
            } else {
              //do filtering here
              //simple cutoff for now
              for(int i=0;i<clength;i++)
              {
                ctemp[i]/=clength; //normalize
                if( ((i<=clength/2) && (i<low || i>=hi)) || 
                    ((i>clength/2)  && ((clength-i-1)<low || (clength-i-1)>=hi) )) //cutoff frequencies, don't forget about 'negative' freq
                  ctemp[i]=_complex(0,0);
              }
              fftwf_execute(_iplan);
            }
            //write out result
            for(int i=0;i<rdr.info()[3].length;i++)
            {
              out.value(ctemp[i].real());
              out.next();
            }
          }
      fftwf_destroy_plan(_plan);
      fftwf_destroy_plan(_iplan);
      std::cout<<"done!"<<std::endl;
      
    } else if(rdr.info().size()==3) {
      std::cerr<<"This file should have time dimension as fastest varying!"<<std::endl;
      return 1;
    }
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    cerr << "errno="<<errno<<std::endl;
    return 1;
  }
	return 0;
}
