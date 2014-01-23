#include "Point.h"
#include "Castor.h"

#include <TH2D.h>
#include <TGeoManager.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPaveText.h>

#include <iostream>
using namespace std;

int main()
{

  CastorHalf farHalf(true);
  Point farMeasuredTop(-0.0859*1000.,0.3115*1000.,-15.8858*1000.);
  Point farMeasuredMid(-0.3195*1000.,0.0846*1000.,-14.8784*1000.);
  Point farMeasuredBottom(-0.3148*1000.,-0.0968*1000.,-15.8888*1000.);
  Point farMeasuredTopFeb(-0.0869*1000.,0.3108*1000.,-15.8865*1000.);
  Point farMeasuredMidFeb(-0.3211*1000.,0.0837*1000.,-14.8791*1000.);
  Point farMeasuredBottomFeb(-0.3148*1000.,-0.0972*1000.,-15.8895*1000.);
  farHalf.SetRelTarget1(farMeasuredTopFeb);
  farHalf.SetRelTarget2(farMeasuredMidFeb);
  farHalf.SetRelTarget3(farMeasuredBottomFeb);


  CastorHalf nearHalf(false);
  Point nearMeasuredTop(0.0963*1000., 0.3108*1000., -15.8900*1000.);
  Point nearMeasuredMid(0.3240*1000., 0.0842*1000., -14.8802*1000.);
  Point nearMeasuredBottom(0.3264*1000., -0.0999*1000., -15.8913*1000.);
  Point nearMeasuredTopFeb(0.0932*1000., 0.3107*1000., -15.8898*1000.);
  Point nearMeasuredMidFeb(0.3253*1000., 0.0850*1000., -14.8811*1000.);
  Point nearMeasuredBottomFeb(0.3233*1000., -0.0986*1000., -15.8913*1000.);
  nearHalf.SetRelTarget1(nearMeasuredTopFeb);
  nearHalf.SetRelTarget2(nearMeasuredMidFeb);
  nearHalf.SetRelTarget3(nearMeasuredBottomFeb);

  cout << endl << endl << " ---Far:---" << endl << endl;
  farHalf.Fit(0);
  cout << endl << endl << " ---NEAR:---" << endl << endl;
  nearHalf.Fit(0);
  
  TPaveText* txt_near = new TPaveText(0.6,0.7,0.8,0.85,"NDC b t l");
  txt_near->SetTextFont(42);
  txt_near->SetFillStyle(0);
  txt_near->SetTextColor(kBlack);
  txt_near->SetTextSize(0.033);
  txt_near->SetBorderSize(0);
  ostringstream txt_near_ss;
  txt_near_ss.precision(2);
  txt_near_ss.str(""); txt_near_ss<< "Near side"; txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "x=" << nearHalf.GetCenter()->GetX(); txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "y=" << nearHalf.GetCenter()->GetY(); txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "#theta=" << nearHalf.GetTheta()/TMath::DegToRad(); txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "#rho=" << nearHalf.GetRho()/TMath::DegToRad(); txt_near->AddText(txt_near_ss.str().c_str());

  TPaveText* txt_far = new TPaveText(0.2,0.7,0.4,0.85,"NDC b t l");
  txt_far->SetTextFont(42);
  txt_far->SetFillStyle(0);
  txt_far->SetTextColor(kBlack);
  txt_far->SetTextSize(0.033);
  txt_far->SetBorderSize(0);
  ostringstream txt_far_ss;
  txt_far_ss.precision(2);
  txt_far_ss.str(""); txt_far_ss<< "Far side"; txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "x=" << farHalf.GetCenter()->GetX(); txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "y=" << farHalf.GetCenter()->GetY(); txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "#theta=" << farHalf.GetTheta()/TMath::DegToRad(); txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "#rho=" << farHalf.GetRho()/TMath::DegToRad(); txt_far->AddText(txt_far_ss.str().c_str());


  TCanvas can1;
  can1.SetGrid();
  TH2D* frame = new TH2D("frame",";x [mm];y [mm]",100,-15,15,100,-15,15);
  frame->SetStats(kFALSE);
  frame->Draw("AXIS");
  farHalf.DrawXY(kRed);
  nearHalf.DrawXY(kBlue);
  can1.SaveAs("plotZoom.pdf");


  TCanvas can2;
  can2.SetGrid();
  TH2D* grid = new TH2D("frame",";x [mm];y [mm]",100,-350,350,100,-350,350);
  frame->SetStats(kFALSE);
  grid->Draw("AXIS");
  farHalf.DrawXY(kRed);
  nearHalf.DrawXY(kBlue);
  txt_far->Draw("SAME");
  txt_near->Draw("SAME");
  can2.SaveAs("plot.pdf");
  can2.SaveAs("plot.root");


  TCanvas can3;
  TGeoManager* geoM = new TGeoManager("world", "world");
  TGeoMaterial* mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium* med = new TGeoMedium("Vacuum",1,mat);
  TGeoVolume* geoWorld = geoM->MakeBox("world_volume", med, 2000,2000,2000);
  geoM->SetTopVolume(geoWorld);
  nearHalf.Draw(geoM);
  farHalf.Draw(geoM);
  geoM->Draw();
  geoM->CloseGeometry();
  geoM->Export("test.C");
  std:: cout << "bla" << std:: endl;
  //can3.SaveAs("plot.pdf");
  
  return 0;
}
