/*
  Programmer: David Mattli
*/

#include <cstdio>
#include <cstring>


#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "reprojector.hh"
#include "rasterreader.hh"

#include <gdal_priv.h>


double params[15] =  { 6370997.000000, 
		       0, 0, 0, 0, 
		       0, 0, 0, 0, 0, 
		       0, 0, 0, 0, 0};



int main(int argc, char *argv[]) 
{
	double ul_x, ul_y, lr_x, lr_y;
 	long int rows, cols;

	vector<unsigned char> *dat = 0;
	PRProcess prc(MPI_COMM_WORLD);
	Reprojector *re = 0;
	string output_filename = "/home/dmattli/Desktop/test/test.tif";

	prc.init(argc, argv);

	rows = cols = 0;
//	ProjectedRaster in("/home/dmattli/Desktop/mmr/glc_geographic_30sec.img");
	ProjectedRaster in("/home/dmattli/Desktop/mmr/veg_geographic_1deg.img");
	if (in.isReady() == true) {
	  printf("Image opened!\n");
	} else {
	  prc.abort();
	  return 1;
	}
	
	Projection *outproj;
	outproj = new Mollweide(params, METER, in.getDatum());
	
		/*
		  FindMinBox(&in, outproj, in.bitsPerPixel(), ul_x, ul_y, lr_x, lr_y);
		  FindMinBox(&in, outproj, in.bitsPerPixel(), ul_x, ul_y, lr_x, lr_y);
		  rows = (ul_y-lr_y) / in.getPixelSize();
		  cols = (lr_x-ul_x) / in.getPixelSize();


		  ProjectedRaster out("/home/dmattli/Desktop/test/test.tif", 
		  in.getRowCount(), in.getColCount(), 
		  in.getPixelType(), in.getPixelSize(), 
		  in.bandCount(), outproj, ul_x, ul_y);
		*/
	if (!prc.isMaster()) {
		output_filename = "";
	}
		ProjectedRaster out(output_filename,
				    &in,
				    outproj,
				    in.getPixelType(),
				    in.getPixelSize());

	
	if (!(in.isReady() && out.isReady())) {
		printf("Error in opening rasters\n");
		return 1;
	}
	
	re = new Reprojector(prc, &in, &out);
	printf("Beginnning reprojection...\n");
//	re->parallelReproject();
	
	// Cleanup
	//		delete re;
	//		delete outproj;	

	return 0;
}
