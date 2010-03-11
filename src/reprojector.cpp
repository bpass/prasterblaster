
#include <cstdio>

#include <mpi.h>

#include <gdal.h>
#include <gdal_priv.h>

#include "pRPL/prProcess.h"
#include "pRPL/neighborhood.h"
#include "pRPL/cellSpace.h"
#include "pRPL/layer.h"
#include "pRPL/basicTypes.h"

#include "reproject_transition.hh"

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "reprojector.hh"
#include "resampler.hh"

using namespace pRPL;

Reprojector::Reprojector(PRProcess _prc, ProjectedRaster *_input, ProjectedRaster *_output) 
	:
	prc(_prc), input(_input), output(_output)
{
	vector<CellCoord> coords;
	double ul_x, ul_y, lr_x, lr_y;
	int slice_rows = output->getRowCount() / prc.nPrcs();
	int in_beg_row = 0;

	maxx = maxy = 0;
	minx = miny = 1e+37;
	//	resampler = &resampler::nearest_neighbor<unsigned char>;
	
	for (int i=0; i<2; ++i) {
	  for (int j=0; j<2; ++j) {
	    coords.push_back(CellCoord(j, i));
	  }
	}

	ul_x = ul_y = lr_x = lr_y = 0.0;
	
	outraster.resize(output->getColCount() * slice_rows);
	FindMinBox2(output->ul_x, output->ul_y - (output->getPixelSize() * slice_rows * prc.id()),
		    output->getPixelSize(), 
		    slice_rows,
		    output->getColCount(), 
		    output->getProjection(), input->getProjection(),
		    input->getPixelSize(),
		    ul_x, ul_y, lr_x, lr_y);


	in_sub_rows = (ul_y - lr_y) / input->getPixelSize() + 1;
	in_sub_cols = (lr_x - ul_x) / input->getPixelSize() + 1;
	inraster.resize(in_sub_rows * in_sub_cols);
	inraster_ul_x = ul_x;
	inraster_ul_y = ul_y;
	in_beg_row = (input->ul_y - ul_y) / input->getPixelSize();

	if (prc.isMaster()) {
		outraster.resize(input->getRowCount() * input->getColCount());
	}
	
	// Read input raster
	input->readRaster(in_beg_row, in_sub_rows, &(inraster[0]));

	

	return;
};

Reprojector::~Reprojector()
{


return;
}

void Reprojector::parallelReproject()
{

	Transformer t;
	Coordinate temp, temp2;
	double in_ulx, in_uly, out_ulx, out_uly;
	double in_pixsize, out_pixsize;
	long in_rows, in_cols, out_rows, out_cols;
	vector<char> ot;
	Projection* outproj = output->getProjection();

	begin_time = MPI_Wtime();

	t.setInput(*output->getProjection());
	t.setOutput(*input->getProjection());

  
	in_rows = input->getRowCount();
	in_cols = input->getColCount();
	out_rows = output->getRowCount()/prc.nPrcs();
	out_cols = output->getColCount();
	in_pixsize = input->getPixelSize();
	out_pixsize = input->getPixelSize();
	in_ulx = input->ul_x;
	in_uly = inraster_ul_y;
	out_ulx = output->ul_x;
	out_uly = output->ul_y - (prc.id() * ((out_rows * out_pixsize)));
	temp.units = METER;


	for (int y = out_rows/2; y < out_rows; ++y) {
/*
		if (y % 100 == 0)
			printf("Node %d On row %d of %d\n", prc.id(), 
			       y + (prc.id() * (output->getRowCount()/prc.nPrcs())), 
			       output->getRowCount());
*/
		for (int x = 0; x < out_cols; ++x) {
			temp.x = ((double)x * out_pixsize) + out_ulx;
			temp.y = out_uly - ((double)y * out_pixsize);
			t.transform(&temp);
			outproj->inverse(temp.x, temp.y);
			outproj->forward(outproj->lon(),
					 outproj->lat(), 
					 &(temp2.x), &(temp2.y));
			temp.x = ((double)x * out_pixsize) + out_ulx;	
			temp.y = out_uly - ((double)y * out_pixsize);

			t.transform(&temp);

			// temp now contains coords to input projection
			temp.x -= inraster_ul_x;
			//temp.y += inraster_ul_y;
			temp.y = inraster_ul_y - temp.y;
			temp.x /= in_pixsize;
			temp.y /= out_pixsize;
			// temp is now scaled to input raster coords, now resample!
			//			printf("Pixsize! %f %f\n", in_pixsize, out_pixsize);
			//			printf("About to sample the point: (%f %f) to (%d %d)\n",
			//			       temp.x, temp.y, x, y);
/*
			if (((int)temp.x < 0) || ((int)temp.x > in_sub_cols)) {
				printf("Node %d: X out of range! %f %d\n", prc.id(),
				       temp.x, in_sub_cols);
				return;
			}

			if (((int)temp.y < 0) || ((int)temp.y > in_sub_rows)) {
				printf("Node %d: Y out of range! %f %d\n", prc.id(),
				       temp.y, in_sub_rows);
				return;
			}
				
*/	
			resampler::bilinear<unsigned char>(&(inraster[0]), 
								   temp.x, temp.y, in_cols, 
								   &(outraster[0]), x, y, out_cols);



		}
	}

	MPI_Gather(&(outraster[0]), out_rows*out_cols, MPI_BYTE, 
		   &(outraster[0]), out_rows*out_cols, MPI_BYTE,
		   0, prc.comm());

	if (prc.isMaster()) {
		printf("Writing output raster.\n");
		output->writeRaster(0, output->getRowCount(), &(outraster[0]));
	}
	

	if (outproj)
		delete outproj;


	end_time = MPI_Wtime();
	if (prc.isMaster()) {
	printf("\n\n__--^^ Reprojection took a total of %f seconds ^^--__\n\n", 
	       end_time - begin_time);
	}

	return;
}

void Reprojector::reproject()
{
	parallelReproject();

}

void FindMinBox2(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize,
		 double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y)
{
	double ul_x, ul_y;
	double pixsize; // in METER!
	Coordinate temp;
	Transformer t;
	double maxx, minx, maxy, miny;
	double ul_lon, ul_lat, lr_lon, lr_lat;

	maxx = maxy = 0;
	minx = miny = 1e+37;

	ul_x = in_ul_x;
	ul_y = in_ul_y;
	pixsize = out_pixsize;
  
	t.setInput(*inproj);
	t.setOutput(*outproj);

	// Find geographic corners of input
	inproj->inverse(ul_x, ul_y, &ul_lon, &ul_lat);
	inproj->inverse(ul_x+(cols*pixsize), 
					ul_y-(rows*pixsize), 
					&lr_lon, 
					&lr_lat);
	double delta_east = (lr_lon-ul_lon)/cols, 
		delta_north = (ul_lat-lr_lat)/rows;

	// Calculate minbox
	temp.x = ul_x;
	temp.y = ul_y;
	temp.units = METER;

	// Check top of map
	for (int x = 0; x < cols; ++x) {
		temp.x = (double)x*pixsize + ul_x;
		temp.y = ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (x*delta_east);
		temp.y = ul_lat;
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check bottom of map
	for (int x = 0; x < cols; ++x) {
		temp.x = (double)x*pixsize + ul_x;
		temp.y = (double)-rows*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (x * delta_east);
		temp.y = ul_lat - (rows*delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;

	}
 
	// Check Left side
	for (int y = 0; y < rows; ++y) {
		temp.x = ul_x;
		temp.y = (double)-y*(pixsize+1) + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon;
		temp.y = ul_lat - (y*delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check right side
	for (int y = 0; y < rows; ++y) {
		temp.x = (double)cols*(1+pixsize) + ul_x;
		temp.y = (double)-y*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (cols * delta_east);
		temp.y = ul_lat - (y * delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Set outputs
	_ul_x = minx;
	_ul_y = maxy;
	_lr_x = maxx;
	_lr_y = miny;
  


}

// Minbox finds number of rows and columns in output and upper-left
// corner of output
void FindMinBox(ProjectedRaster *input, Projection *outproj, double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y)
{
	double ul_x, ul_y;
	double pixsize; // in METER!
	Coordinate temp;
	Transformer t;
	int cols, rows;
	double maxx, minx, maxy, miny;
	double ul_lon, ul_lat, lr_lon, lr_lat;

	maxx = maxy = 0;
	minx = miny = 1e+37;

	ul_x = input->ul_x;
	ul_y = input->ul_y;
	cols = input->getColCount();
	rows = input->getRowCount();
	pixsize = out_pixsize;
  
	t.setInput(*input->getProjection());
	t.setOutput(*outproj);

	// Find geographic corners of input
	input->getProjection()->inverse(ul_x, ul_y, &ul_lon, &ul_lat);
	input->getProjection()->inverse(ul_x+(cols*pixsize), 
					ul_y-(rows*pixsize), 
					&lr_lon, 
					&lr_lat);
	double delta_east = (lr_lon-ul_lon)/cols, 
		delta_north = (ul_lat-lr_lat)/rows;

	// Calculate minbox
	temp.x = ul_x;
	temp.y = ul_y;
	temp.units = METER;

	// Check top of map
	for (int x = 0; x < cols; ++x) {
		temp.x = (double)x*pixsize + ul_x;
		temp.y = ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (x*delta_east);
		temp.y = ul_lat;
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check bottom of map
	for (int x = 0; x < cols; ++x) {
		temp.x = (double)x*pixsize + ul_x;
		temp.y = (double)-rows*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (x * delta_east);
		temp.y = ul_lat - (rows*delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;

	}
 
	// Check Left side
	for (int y = 0; y < rows; ++y) {
		temp.x = ul_x;
		temp.y = (double)-y*(pixsize+1) + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon;
		temp.y = ul_lat - (y*delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check right side
	for (int y = 0; y < rows; ++y) {
		temp.x = (double)cols*(1+pixsize) + ul_x;
		temp.y = (double)-y*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
		temp.x = ul_lon + (cols * delta_east);
		temp.y = ul_lat - (y * delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Set outputs
	_ul_x = minx;
	_ul_y = maxy;
	_lr_x = maxx;
	_lr_y = miny;
  

	return;
}

