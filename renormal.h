#ifndef _Renormal_h_
#define _Renormal_h_

//----------------------------------------------------------------------//
//  Renormal: calculate fundamental matrix by renormalization           //
//                                                                      //
//  1999/2/12                                                           //
//  Programmed by Hitoshi Mishima                                       //
//  Department of Computer Science,                                     //
//  Gunma University Faculty of Engineering                             //
//----------------------------------------------------------------------//

int Renormal(ImagePosition* P1, ImagePosition* P2, int point_num, 
	     VectorN& eigen_val, Matrix* eigen_mat, double& J, 
	     int dim = 1, double min = 1.0e-8, int itmax = 50);

//  P1        : array of image points(camera1).         input
//  P2        : array of image points(camera2).         input
//  point_num : number of image points.                 input
//  eigen_val : eigenvalues of moment tensor            output
//  eigen_mat : eigenmatrices of moment tensor          output
//  J         : residual of renormalization             output
//  dim       : dimension of renormalization            input
//               (1 or 2, default: 1)
//  min       : constant for convergence test           input
//               (default: 1.0e-8)
//  itmax     : maximum iteration number                input
//               (>=1, default: 50)
//
//Return value
//  iteration number

#endif //_Renormal_h_
     
