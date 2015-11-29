//#include "minc_wrappers.h"
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
#include "minc_1_simple.h"

using namespace  std;
using namespace  minc;

typedef std::complex<float> _complex;

void normalize(std::vector<float>& array)
{
  if(array.empty())
    return;
  double mean=0.0;
  double sd=0.0;
  for(std::vector<float>::iterator i=array.begin();i!=array.end();i++)
  {
    mean+=*i;
    sd+=(*i)*(*i);
  }
  mean/=array.size();
  sd/=array.size();
  sd-=mean*mean;
  sd=sqrt(sd);
  for(std::vector<float>::iterator i=array.begin();i!=array.end();i++)
  {
    (*i)=((*i)-mean)/sd;
  }
}

class _fft
{
  protected:
    fftwf_plan _plan,_iplan;

  public:
    typedef std::vector<_complex> complex_array;
    complex_array buffer;
    _fft(int length):buffer(length)
    {
      _plan  = fftwf_plan_dft_1d(length, reinterpret_cast<fftwf_complex*>(&buffer[0]), reinterpret_cast<fftwf_complex*>(&buffer[0]), FFTW_FORWARD,  FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
      _iplan = fftwf_plan_dft_1d(length, reinterpret_cast<fftwf_complex*>(&buffer[0]), reinterpret_cast<fftwf_complex*>(&buffer[0]), FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
    }
    ~_fft()
    {
      fftwf_destroy_plan(_plan);
      fftwf_destroy_plan(_iplan);
    }
    
    void forward(void)
    {
      fftwf_execute(_plan);
    }
    
    void backward(void)
    {
      fftwf_execute(_iplan);
    }
    
    void put(const std::vector<float>& array)
    {
      //buffer.fill(_complex(0,0);
      complex_array::iterator j;
      std::vector<float>::const_iterator i;
      for(j=buffer.begin();j!=buffer.end();j++)
        (*j)=_complex(0,0);

      for(i=array.begin(),j=buffer.begin();i!=array.end() && j!=buffer.end();i++,j++)
      {
        (*j)=_complex(*i,0);
      }
    }
    
    void get(std::vector<float>& array) const
    {
      std::vector<float>::iterator i;
      complex_array::const_iterator j;
      for(i=array.begin();i!=array.end();i++)
        (*i)=0.0;
      for(i=array.begin(),j=buffer.begin();i!=array.end()&&j!=buffer.end();i++,j++)
      {
        *i=(*j).real();
      }
    }
    
    //square bandpass filter
    void filter(int low,int hi)
    {
      size_t length=buffer.size();
      for(int i=0;i<length;i++)
      {
        buffer[i]/=length; //normalize
        
        //cutoff frequencies, don't forget about 'negative' freq
        if( (( i<=length/2 ) && (  i<low || i >= hi)) ||  
            (( i>length/2  ) && ( (length-i-1)<low || (length-i-1) >= hi)) ) 
          buffer[i]=_complex(0,0);
      }
    }
    
    void normalize(void)
    {
      double mean=0.0;
      double sd=0.0;
      for(complex_array::iterator i=buffer.begin();i!=buffer.end();i++)
      {
        mean+=(*i).real();
        sd+=(*i).real()*(*i).real();
      }
      mean/=buffer.size();
      sd/=buffer.size();
      sd-=mean*mean;
      sd=sqrt(sd);
      for(complex_array::iterator i=buffer.begin();i!=buffer.end();i++)
      {
        (*i)=_complex(((*i).real()-mean)/sd,0);
      }
    }
};

void show_usage (const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <input1> <input_tamplate.txt>  <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl
    << "--low <Hz>"<<endl
    << "--hi <Hz>"<<endl
    << "--shift <f> z slice of the sample" <<endl;
  
}

int main (int argc, char **argv)
{
  int verbose=1;
  int clobber=0;
  int low=0,hi=0;
  double low_freq=0.0;
  double hi_freq=0.0;
  double sample_shift=0.0;
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"low",  required_argument, 0, 'l'},
    {"hi",   required_argument, 0, 'h'},
    {"shift",   required_argument, 0, 's'},
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
			//low = atoi(optarg);
      low_freq =atof(optarg);
			break;
		case 'h':
			//hi = atoi(optarg);
      hi_freq=atof(optarg);
			break;
		case 'v':
			cout << "Version: 1.0" << endl;
			return 0;
    case 's':
      sample_shift=atof(optarg);
			break;
      
		case '?':
			/* getopt_long already printed an error message. */
		default:
			show_usage (argv[0]);
			return 1;
		}
	}
  
  if((argc - optind)<3)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input=argv[optind];
  std::string input_tamp=argv[optind+1];
  std::string output=argv[optind+2];
  
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
    if(rdr.info().size()!=4 || rdr.info()[3].name!="time")
    {
      std::cerr<<"Input file "<< input.c_str() <<" should have time dimension as fastest varying!"<<std::endl;
      return 1;
    }
    
    if(rdr.info()[0].name!="zspace")
    {
      std::cerr<<"Input file "<< input.c_str() <<" should have zspace as slowest varying!"<<std::endl;
      return 1;
    }
    
    std::cout<<"Length:"<<rdr.info()[0].length<<std::endl;
    std::cout<<"Sample shift:"<<sample_shift<<std::endl;
    double  z_phase_shift=-1.0/rdr.info()[0].length;
    sample_shift=-sample_shift/rdr.info()[0].length;
    std::cout<<"Z space phase shift:"<<rdr.info()[3].step*z_phase_shift<<" sec/slice"<<std::endl;
    std::cout<<"sample phase shift:"<<rdr.info()[3].step*sample_shift<<" sec"<<std::endl;
    //z_phase_shift=0.0;
    
    std::cout<<"Processing time dimension..."<<std::endl<<std::flush;
    unsigned int clength=rdr.info()[3].length;
    hi= floor(rdr.info()[3].step*rdr.info()[3].length*hi_freq);
    low=floor(rdr.info()[3].step*rdr.info()[3].length*low_freq);
    std::cout<<"low="<<low<<" Hi="<<hi<<std::endl;
    std::vector<float> temp(clength);
    std::vector<float> sample(clength);
    
    //rdr_tamp.read(sample,clength);
    
    ifstream in_tamp(input_tamp.c_str());
    int i=0;
    while(in_tamp.good() && !in_tamp.eof() && i<clength )
    { 
      in_tamp>>sample[i];
      i++;
    }
    _fft filter(clength);
    filter.put(sample);
    filter.forward();
    filter.filter(low,hi);
    filter.backward();
    //filter.normalize();
    filter.get(sample);
    
    
    _fft convolve(clength*2);
    convolve.put(sample);
    convolve.normalize();
    ofstream out_tamp("filtered.txt");
    for( i=0;i<clength;i++)
      out_tamp<<convolve.buffer[i].real()<<"\n";
    
    convolve.forward();
    _fft::complex_array sample_fft(convolve.buffer); //save a copy
    for(int i=0;i<clength*2;i++)
      sample_fft[i]/=clength*2; //normalize
    
    minc_input_iterator<float> in(rdr);
    minc_output_iterator<float> out(wrt);
    in.begin();
    out.begin();

    for(int z=0;z<rdr.info()[0].length;z++)
      for(int y=0;y<rdr.info()[1].length;y++)
        for(int x=0;x<rdr.info()[2].length;x++)
        {
          //read the column
          for(int i=0;i<clength;i++)
          {
              temp[i]=in.value();
              in.next();
          }
          
          filter.put(temp);
          filter.forward();
          filter.filter(low,hi);
          filter.backward();
          //filter.normalize();
          filter.get(temp);
          
          convolve.put(temp);
          convolve.normalize();
          convolve.forward();
          
          //TODO: calculate phase shift here!
          
          for(int i=0;i<clength*2;i++)
          {
            _complex sample_phase_shift(cos(2*M_PI*sample_shift*i/(clength*2)),sin(2*M_PI*sample_shift*i/(clength*2)));
            _complex phase_shift(cos(2*M_PI*z_phase_shift*i/(clength*2)*z),sin(2*M_PI*z_phase_shift*i/(clength*2)*z));
            convolve.buffer[i]/=clength*2; //normalize
            convolve.buffer[i]*=phase_shift;
            convolve.buffer[i]*=std::conj(sample_fft[i]*sample_phase_shift);
          }
          convolve.backward();
          
          for(int i=0;i<clength;i++)
          {
            if(i<clength/2)
              out.value(convolve.buffer[clength*2-clength/2+i].real());
            else
              out.value(convolve.buffer[i-clength/2].real());
            out.next();
          }
        }

    std::cout<<"done!"<<std::endl;
    
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    
    cerr << "errno="<<errno<<std::endl;
    return 1;
  }
	return 0;
}
