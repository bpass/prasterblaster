/* 
   Programmer: David Mattli
*/

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <vector>


#include "pRPL/prProcess.h"
#include "pRPL/neighborhood.h"
#include "pRPL/cellSpace.h"
#include "pRPL/layer.h"

//#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "resampler.hh"

using namespace pRPL;
using resampler::resampler_func;

class Reprojector
{
public:
	Reprojector(PRProcess prc, 
		    ProjectedRaster *_input, 
		    ProjectedRaster *_output);
	~Reprojector();
	void reproject();
	void parallelReproject();

private:
	double maxx, minx, maxy, miny;
	double inraster_ul_x, inraster_ul_y;
	double begin_time, end_time;
	int in_sub_rows, in_sub_cols;
	ProjectedRaster *input;
	ProjectedRaster *output;
	PRProcess prc;
	vector<unsigned char> inraster;
	vector<unsigned char> outraster;


	resampler_func resampler;
};

void FindMinBox(ProjectedRaster *input, Projection *outproj,
		double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);

void FindMinBox2(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize,
		 double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);

ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
