#include "Point.h"
#include "Castor.h"

#include <TH2D.h>
#include <TGeoManager.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLegend.h>

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
  gStyle->SetOptStat(0);

  CastorHalf farHalf(true);
  CastorHalf nearHalf(false);

  if(argc==2)
    {
      string filename(argv[1]);
      cout << "Opening file: " << filename << endl;
      ifstream file;
      file.open (filename.c_str());
      vector<Point> points(6);

      string xstr,ystr,zstr;
      int i = 0;
      while(file >> xstr >> ystr >> zstr)
        {
          double x = std::atof (xstr.c_str()) * 1000.;
          double y = std::atof (ystr.c_str()) * 1000.;
          double z = std::atof (zstr.c_str()) * 1000.;
          points[i].SetXYZ(x,y,z);
          cout << i << " " << x << " " << y << " " << z << endl;
          i++;
        }
      if(i!=6)
        {
          cerr << "Wrong number of targets in file" <<endl;
          exit(-1);
        }
      farHalf.SetRelTarget1(points[0]);
      farHalf.SetRelTarget2(points[1]);
      farHalf.SetRelTarget3(points[2]);
      nearHalf.SetRelTarget1(points[3]);
      nearHalf.SetRelTarget2(points[4]);
      nearHalf.SetRelTarget3(points[5]);

      file.close();
    }
  else
    {
      cout << "No input file given. Fitting default 2013 Jan measurment" << endl;
      Point farMeasuredTop(-0.0859*1000.,0.3115*1000.,-15.8858*1000.);
      Point farMeasuredMid(-0.3195*1000.,0.0846*1000.,-14.8784*1000.);
      Point farMeasuredBottom(-0.3148*1000.,-0.0968*1000.,-15.8888*1000.);

      Point farMeasuredTopFeb(-0.0860*1000.,0.3107*1000.,-15.8865*1000.);
      Point farMeasuredMidFeb(-0.3202*1000.,0.0836*1000.,-14.8792*1000.);
      Point farMeasuredBottomFeb(-0.3138*1000.,-0.0972*1000.,-15.8895*1000.);

      Point farMeasuredMid2014   (-0.3208*1000.,0.0847*1000.,-14.8775*1000.);
      Point farMeasuredTop2014   (-0.0909*1000.,0.3122*1000.,-15.8847*1000.);
      Point farMeasuredBottom2014(-0.3194*1000.,-0.0954*1000.,-15.8875*1000.);
      farHalf.SetRelTarget1(farMeasuredTop);
      farHalf.SetRelTarget2(farMeasuredMid);
      farHalf.SetRelTarget3(farMeasuredBottom);


      Point nearMeasuredTop(0.0963*1000., 0.3108*1000., -15.8900*1000.);
      Point nearMeasuredMid(0.3240*1000., 0.0842*1000., -14.8802*1000.);
      Point nearMeasuredBottom(0.3264*1000., -0.0999*1000., -15.8913*1000.);

      Point nearMeasuredTopFeb(0.0941*1000., 0.3107*1000., -15.8898*1000.);
      Point nearMeasuredMidFeb(0.3262*1000., 0.0850*1000., -14.8811*1000.);
      Point nearMeasuredBottomFeb(0.3242*1000., -0.0987*1000., -15.8914*1000.);

      Point nearMeasuredMid2014   (0.3180*1000., 0.0865*1000.,  -14.8794*1000.);
      Point nearMeasuredTop2014   (0.0873*1000., 0.3124*1000.,  -15.8880*1000.);
      Point nearMeasuredBottom2014(0.3207*1000., -0.0947*1000., -15.8901*1000.);
      nearHalf.SetRelTarget1(nearMeasuredTop);
      nearHalf.SetRelTarget2(nearMeasuredMid);
      nearHalf.SetRelTarget3(nearMeasuredBottom);
    }

  cout << endl << endl << " ---Far:---" << endl << endl;
  farHalf.Fit(0);
  cout << endl << endl << " ---NEAR:---" << endl << endl;
  nearHalf.Fit(0);

  TPaveText* txt_near = new TPaveText(0.60,0.6,0.75,0.75,"NDC");
  txt_near->SetTextFont(42);
  txt_near->SetFillColor(kWhite);
  txt_near->SetTextColor(kBlue);
  txt_near->SetLineColor(kBlue);
  txt_near->SetTextSize(0.033);
  txt_near->SetBorderSize(1);
  ostringstream txt_near_ss;
  txt_near_ss.precision(2);
  txt_near_ss.str(""); txt_near_ss<< "Near side"; txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "x=" << nearHalf.GetCenter()->GetX() << " mm"; txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "y=" << nearHalf.GetCenter()->GetY() << " mm"; txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "#theta=" << nearHalf.GetTheta()/TMath::DegToRad() << "#circ"; txt_near->AddText(txt_near_ss.str().c_str());
  txt_near_ss.str(""); txt_near_ss<< "#rho=" << nearHalf.GetRho()/TMath::DegToRad() << "#circ"; txt_near->AddText(txt_near_ss.str().c_str());

  TPaveText* txt_far = new TPaveText(0.25,0.6,0.4,0.75,"NDC");
  txt_far->SetTextFont(42);
  txt_far->SetFillColor(kWhite);
  txt_far->SetTextColor(kRed);
  txt_far->SetLineColor(kRed);
  txt_far->SetTextSize(0.033);
  txt_far->SetBorderSize(1);
  ostringstream txt_far_ss;
  txt_far_ss.precision(2);
  txt_far_ss.str(""); txt_far_ss<< "Far side"; txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "x=" << farHalf.GetCenter()->GetX() << " mm"; txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "y=" << farHalf.GetCenter()->GetY() << " mm"; txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "#theta=" << farHalf.GetTheta()/TMath::DegToRad() << "#circ"; txt_far->AddText(txt_far_ss.str().c_str());
  txt_far_ss.str(""); txt_far_ss<< "#rho=" << farHalf.GetRho()/TMath::DegToRad() << "#circ"; txt_far->AddText(txt_far_ss.str().c_str());


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
  grid->SetStats(kFALSE);
  grid->Draw("AXIS");
  farHalf.DrawXY(kRed);
  nearHalf.DrawXY(kBlue);
  txt_far->Draw("BR");
  txt_near->Draw("BR");
  TLegend leg(0.4,0.25,0.6,0.4);
  //histogram because root doesn't take attributes if new TMarker is initialised here
  TH1D* helper = new TH1D("a","a",0,0,1); helper->SetMarkerStyle(20); helper->SetMarkerSize(2); helper->SetMarkerColor(kRed);
  leg.AddEntry(helper,"shift","p");
  helper = new TH1D("a","a",0,0,1); helper->SetMarkerStyle(24); helper->SetMarkerSize(2); helper->SetMarkerColor(kRed);
  leg.AddEntry(helper,"measured target positions","p");
  helper = new TH1D("a","a",0,0,1); helper->SetMarkerStyle(5); helper->SetMarkerSize(2); helper->SetMarkerColor(kRed);
  leg.AddEntry(helper,"nominal target pos + shift","p");
  leg.SetFillColor(kWhite);
  leg.Draw();
  can2.SaveAs("plot.pdf");
  can2.SaveAs("plot.C");
  can2.SaveAs("plot.png");



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
  //can3.SaveAs("plot.pdf");

  cout << endl
       << "***************************RESULTS*****************************" << endl
       << "* " << "Far side" << endl
       << "* " << "x=" << farHalf.GetCenter()->GetX() << endl
       << "* " << "y=" << farHalf.GetCenter()->GetY() << endl
       << "* " << "theta=" << farHalf.GetTheta()/TMath::DegToRad() << " degree" << endl
       << "* " << "rho=" << farHalf.GetRho()/TMath::DegToRad() << " degree" << endl
       << "* " << "x (nonIP)=" << farHalf.GetCenterNonIP()->GetX() << endl
       << "* " << "y (nonIP)=" << farHalf.GetCenterNonIP()->GetY() << endl
       << "***************************************************************" << endl;
  cout << "* " << "Near side" << endl
       << "* " << "x=" << nearHalf.GetCenter()->GetX() << endl
       << "* " << "y=" << nearHalf.GetCenter()->GetY() << endl
       << "* " << "theta=" << nearHalf.GetTheta()/TMath::DegToRad() << " degree" << endl
       << "* " << "rho=" << nearHalf.GetRho()/TMath::DegToRad() << " degree" << endl
       << "* " << "x (nonIP)=" << nearHalf.GetCenterNonIP()->GetX() << endl
       << "* " << "y (nonIP)=" << nearHalf.GetCenterNonIP()->GetY() << endl
       << "***************************RESULTS*****************************" << endl << endl;

  return 0;
}
