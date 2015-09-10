#include <inttypes.h>
#include <stdio.h>

#include <iostream>
#include <unistd.h>
#include <cstdlib>
// #include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>
#include <dirent.h>

#include <bh3defaultgrid.h>

#include <EXP2D_MatrixData.h>
#include <main.h>
#include <EXP2D_tools.h>
#include <EXP2D_binaryfile.h>
#include <EXP2D_rk4.hpp>
#include <EXP2D_runner.hpp>
#include <EXP2D_evaluation.h>
#include <plot_with_mgl.h>
#include <EXP2D_startgrids.h>


#include <TF2.h>
#include <TH2.h>
#include <TMath.h>

using namespace std;


// Fitting a 2-D histogram
// This tutorial illustrates :
//  - how to create a 2-d function
//  - fill a 2-d histogram randomly from this function
//  - fit the histogram
//  - display the fitted function on top of the histogram 
//
// This example can be executed via the interpreter or ACLIC
//   root > .x fit2.C
//   root > .x fit2.C++
//Author: Rene Brun
         
Double_t g2(Double_t *x, Double_t *par) { //This seems to be a Gaussian function. Why does the graph look weird?
   Double_t r1 = Double_t((x[0]-par[1])/par[2]);
   Double_t r2 = Double_t((x[1]-par[3])/par[4]);
   return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}   
Double_t fun2(Double_t *x, Double_t *par) { //I don't think I know what these are doing.
   Double_t *p1 = &par[0];
   Double_t *p2 = &par[5];
   Double_t *p3 = &par[10];
   //Double_t result = g2(x,p1) + g2(x,p2) + g2(x,p3);
   Double_t result = g2(x,p1) + g2(x,p2); //This is the formula, but why? Is it just adding each Gaussian, x and y? Removing the third term gets rid of the weird side bump.
   return result;
}

void fit2() {
   const Int_t npar = 15;
   Double_t f2params[npar] = 
      {100,-3,3,-3,3,160,0,0.8,0,0.9,40,4,0.7,4,0.7}; //what on earth are all these parameters for?
   TF2 *f2 = new TF2("f2",fun2,-10,10,-10,10, npar); //"fun2" seems to put together the Gaussian defined in g2.
   f2->SetParameters(f2params);

   //Create an histogram and fill it randomly with f2
   
   //TH2F *h2 = new TH2F("h2","from f2",40,-10,10,40,-10,10);
   //h2->SetOption("surf2z");
   Int_t nentries = 100000;
   //h2->FillRandom("f2",nentries);
  
   
   TH2F *h2 = new TH2F("h2","2DGaussian",40,-10,10,40,-10,10);
  h2->SetOption("surf2z");
  h2->GetXaxis()->SetTitle("My X Title");
  h2->GetYaxis()->SetTitle("My Y Title");
  int binx1 = h2->GetXaxis()->FindBin(-.3);
  int binx2 = h2->GetXaxis()->FindBin(.3);
  int biny1 = h2->GetYaxis()->FindBin(-.3);
  int biny2 = h2->GetYaxis()->FindBin(.3);
  TRandom3 foo(12345);
  for(int i=0;i<10000;i++)
    {
      h2->Fill(foo.Gaus(0,1),foo.Gaus(0,1));
    }
   
   //Fit h2 with original function f2
   Float_t ratio = 4*nentries/100000;
   f2params[ 0] *= ratio; //I have no idea what these are for. Scaling something to the number of entries?
   f2params[ 5] *= ratio;
   f2params[10] *= ratio;
   f2->SetParameters(f2params);
   h2->Fit("f2");
   f2->Draw("same, SURF");
}

int main(int argc, char *argv[])
{	
	fit2();

    return 0;
}