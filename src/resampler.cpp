

#include <cmath>

#include "resampler.hh"

namespace resampler
{

  template <typename T>
    void nearest_neighbor(void *inraster, double in_x,
			  double in_y, unsigned long in_cols,
			  void *outraster,
			  unsigned long out_x, unsigned long out_y,
			  unsigned long out_cols)
  {
	  T* __restrict in = (T*)inraster;
	  T* __restrict out = (T*)outraster;
	  T value = in[(unsigned long)in_x + ((unsigned long)in_y * in_cols)];
	  
	  out[out_x + (out_y * out_cols)] = value;
	  return;
  }
  

  template <typename T>
  void bilinear(void *inraster, double in_x,
		double in_y, unsigned long in_cols,
		void *outraster,
		unsigned long out_x, unsigned long out_y,
		unsigned long out_cols)
  {
	  double x1, x2, y1, y2;
	  double ul, lr, ll, ur;
	  double val = 0;
	  double denom = 0;
	  T* __restrict in = (T*)inraster;
	  T* __restrict out = (T*)outraster;	  
	  
	  x1 = ceil(in_x);
	  x2 = floor(in_x);
	  y1 = floor(in_y);
	  y2 = ceil(in_y);
	  ul = in[(int)(x2 + (y1 * in_cols))];
	  lr = in[(int)(x1 + (y2 * in_cols))];
	  ll = in[(int)(x2 + (y2 * in_cols))];
	  ur = in[(int)(x1 + (y1 * in_cols))];
		  denom = (x2 - x1) * (y2 - y1);
		  
		  
	val += (ll / denom) * (x2 - in_x) * (y2 - in_y);
	val += (lr / denom) * (in_x - x1) * (y2 - in_y);
	val += (ul / denom) * (x2 - in_x) * (in_y - y1);
	val += (ur / denom) * (in_x - x1) * (in_y - y1);
	
		  out[out_x + (out_y * out_cols)] = (T)val;
	
	  
	  
	  return;
	  
  }
  

// This function is to coax the compiler to actually generate some
// code. Don't call it.

  void never_call_me()
  {
	  nearest_neighbor<unsigned char>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<signed char>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<int>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<unsigned int>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<long>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<unsigned long>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<float>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<double>(0,0,0,0,0,0,0,0);

	  
	  bilinear<unsigned char>(0,0,0,0,0,0,0,0);
	  bilinear<signed char>(0,0,0,0,0,0,0,0);
	  bilinear<int>(0,0,0,0,0,0,0,0);
	  bilinear<unsigned int>(0,0,0,0,0,0,0,0);
	  bilinear<long>(0,0,0,0,0,0,0,0);
	  bilinear<unsigned long>(0,0,0,0,0,0,0,0);
	  bilinear<float>(0,0,0,0,0,0,0,0);
	  bilinear<double>(0,0,0,0,0,0,0,0);

  }

}
