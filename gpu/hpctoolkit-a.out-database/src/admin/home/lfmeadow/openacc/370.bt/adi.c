//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB BT code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"

#include <time.h>
void adi()
{
  struct timespec ts1, ts2;

  clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
  compute_rhs();
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
  double t0 = (ts2.tv_sec + 1e-9*ts2.tv_nsec) - (ts1.tv_sec + 1e-9*ts1.tv_nsec);

  x_solve();
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
  double t1 = (ts1.tv_sec + 1e-9*ts1.tv_nsec) - (ts2.tv_sec + 1e-9*ts2.tv_nsec);

  y_solve();
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
  double t2 = (ts2.tv_sec + 1e-9*ts2.tv_nsec) - (ts1.tv_sec + 1e-9*ts1.tv_nsec);

  z_solve();
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
  double t3 = (ts1.tv_sec + 1e-9*ts1.tv_nsec) - (ts2.tv_sec + 1e-9*ts2.tv_nsec);

  add();
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
  double t4 = (ts2.tv_sec + 1e-9*ts2.tv_nsec) - (ts1.tv_sec + 1e-9*ts1.tv_nsec);
  printf("adi: %.1f %.1f %.1f %.1f %.1f\n", 1000*t0, 1000*t1, 1000*t2, 1000*t3, 1000*t4);
}
