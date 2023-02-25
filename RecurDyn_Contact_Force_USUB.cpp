#include "stdafx.h"

#include "DllFunc.h"

#include "math.h"

 

ContactForce_API void __cdecl contact_force(double time, double upar[], int npar, double pen, double rvel[], int jflag, int iflag, double* result)

{

   using namespace rd_syscall;

   // Parameter Information

   //   time   : Simulation time of RD/Solver. (Input)

   //   upar   : Parameters defined by user. (Input)

   //   npar   : Number of user parameters. (Input)

   //   pen    : Penetration at the contact point. (Input)

   //   rvel    : Relative velocity between contact pair w.r.t. contact reference frame. (Input, Size : 3)

   //   jflag   : When RD/Solver evaluates a Jacobian, the flag is true. (Input)

   //   iflag   : When RD/Solver initializes arrays, the flag is true. (Input)

   //   result  : Returned contact force vector. (Output, Size : 6)

 

   // User Statement

   double dpen, K, C, M1, M2, M3, STV, DTV, SFC, DFC;

   double SFO, DFO, TVEL, FC[3], ZERO, MU;

   int i,ERRFLG;

 

   dpen = rvel[2];

   K = upar[0];

   C = upar[1];

   M1 = upar[2];

   M2 = upar[3];

   M3 = upar[4];

   STV = upar[5];

   DTV = upar[6];

   SFC = upar[7];

   DFC = upar[8];

   ZERO = 0.0;

 

   //---- COMPUTE CONTACT NORMAL FORCE

   if ( fabs(pen) < 1.0e-12 ) {

      SFO = 0.0;

   }

   else {

      SFO = K * pow(fabs(pen),M1);

   }

 

   if ( fabs(dpen) < 1.0e-12 ) {

      DFO = 0.0;

   }

   else {

      DFO = C * pow(fabs(pen),M3) * pow(fabs(dpen),M2);

   }

 

   if ( dpen > ZERO ) DFO = -DFO;

   FC[2] = SFO + DFO;

 

   //---- COMPUTE FRICTION FORCE

   TVEL = sqrt(rvel[0]*rvel[0] + rvel[1]*rvel[1]);

   FC[0] = 0.0;

   FC[1] = 0.0;

   if ( TVEL >= 1.0e-12 ) {

      if ( TVEL < STV ) {

         rd_havsin(TVEL,-STV,-SFC,STV,SFC,0,&MU,&ERRFLG);

      }

      else {

         rd_havsin(TVEL,STV,SFC,DTV,DFC,0,&MU,&ERRFLG);

      }

 

      double frictionForce = MU * fabs(FC[2]);

      if (TVEL == 0.0) {

         FC[0] = 0.0;

         FC[1] = 0.0;

      }

      else {

         FC[0] = -frictionForce*rvel[0]/TVEL;

         FC[1] = -frictionForce*rvel[1]/TVEL;

      }

   }

 

   for(i=0;i<3;i++) result[i] = FC[i];

}
