
#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void encodeOneStep(const char* filename, const unsigned char* image, unsigned width, unsigned height) {
  /*Encode the image*/
  unsigned error = lodepng_encode32_file(filename, image, width, height);

  /*if there's an error, display it*/
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
}





int colour_(int *nwp, double fluxjac[][*nwp]) {

    printf("\n\nCall from the C program: Start colouring\n");



    int nw = *nwp;
    int i, j;
    double rf, gf, bf;
    double rf2, gf2, bf2;
	
    unsigned width = nw, height = nw;
    unsigned char* image = malloc(width * height * 4 );
    const char* filename = "flux-jacobian.png";

    int nl = 1, R;
    int nt = nw + 2*nl;
    int i0, j0, ii, jj;
    double rad, Cval, jmax, jmin, exp, coltol;



  
//  Note that arr[i][j] is same as *(*(arr+i)+j) 
    double *C[nt];
    for (i=0; i<nt; i++) 
         C[i] = (double *)malloc(nt * sizeof(double)); 
  
//    double C[nw+2*nl][nw+2*nl];
  







    exp = 0.05;




    R = nl;


//  initialise shading array with zeros
    for(j = 0; j < nw+2*nl; j++){
    for(i = 0; i < nw+2*nl; i++){
	
	C[i][j] = 0;
	
    }
    }


    jmax = 0.0;
    jmin = 0.0;
//  find max value in flux jacobian
    for(j = R; j < nw+nl; j++){
    for(i = R; i < nw+nl; i++){

    if (fluxjac[i-R][j-R] > jmax) {
	jmax = fluxjac[i-R][j-R];
    }

    if (fluxjac[i-R][j-R] < jmin) {
	jmin = fluxjac[i-R][j-R];
    }

    }
    }


	










coltol = 0.0;

//--------------------------------------------------------------------------------------------
//  initialise shading array with values from flux jacobian
    for(j = R; j < nw+nl; j++){
    for(i = R; i < nw+nl; i++){
	
//	C[i][j] = pow(fabs(fluxjac[i-R][j-R])/jmax,exp);




	if (fabs(fluxjac[i-R][j-R]) > coltol && fabs(fluxjac[i-R][j-R]) != 0.0) {
	    C[i][j] = 1.0; 
	} else {
	    C[i][j] = 0.0; 
	}
	


//	if (fabs(fluxjac[i-R][j-R]) == 0.0) {
//	    C[i][j] = 0.0; 
//	} else {
//	    C[i][j] = 1.0; 
//	}
	
    }
    }



//  if flux jacobian value is not zero, activate smearing
    for(j = R; j < nw+nl; j++){
    for(i = R; i < nw+nl; i++){
	

	i0 = i-R;
        j0 = j-R;	
	    

	if (fluxjac[i-R][j-R] != 0.0) {
	for (ii = 0; ii<=2*nl; ii++) {
	for (jj = 0; jj<=2*nl; jj++) {

	    rad = sqrt((jj-R)*(jj-R) + (ii-R)*(ii-R));
            Cval = 1.0/(rad*rad);
	    //C[i0+ii][j0+jj] = 1.0 - (1.0-C[i0+ii][j0+jj])*(1-Cval);

	}
	}
	}
	
    }
    }



//  start the looping; assign desired colour to relevant pixel
    for(j = 0; j < height; j++){
    for(i = 0; i < width; i++){
	
	rf = 255, gf = 255, bf = 255;
	rf2 = 255, gf2 = 255, bf2 = 255;

	Cval = C[i+R][j+R];

	rf = rf*(1-Cval) + Cval*0;
	gf = gf*(1-Cval);
	bf = bf*(1-Cval);


//        if (fluxjac[i-R][j-R] < 0.0) {
//
//	    rf = rf*(1-Cval);
//	    gf = gf*(1-Cval);
//	    bf = bf*(1-Cval)+Cval*0;
//
 //       }

	unsigned red = rf, green = gf, blue = bf;
	image[4 * width * j + 4 * i + 0] = red;
	image[4 * width * j + 4 * i + 1] = green;
	image[4 * width * j + 4 * i + 2] = blue;
	image[4 * width * j + 4 * i + 3] = 255;
	
    }
    }




//  run the formatting function
    encodeOneStep(filename, image, width, height);
    free(image);

    return 0;


}
