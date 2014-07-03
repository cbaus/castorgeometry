#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TROOT.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

#define _frames 50
#define _do_animation 0
#define _do_fits 1

void animate_pos_mag()
{
  vector<string> sensors;
  sensors.push_back("IP_FAR_TOP");
  sensors.push_back("IP_FAR_BOTTOM");
  sensors.push_back("IP_NEAR_TOP");
  sensors.push_back("IP_NEAR_BOTTOM");

  vector<double> sensors_angle;
  sensors_angle.push_back((0-22.5)*TMath::DegToRad());
  sensors_angle.push_back((90-22.5)*TMath::DegToRad());
  sensors_angle.push_back((180+22.5)*TMath::DegToRad());
  sensors_angle.push_back((90+22.5)*TMath::DegToRad());

  vector<double> sensors_delta;

  TFile* file = TFile::Open("InputSensorNamesHFminus_v1-data_2013_time_fix_with_oracle_db_2_with_lumi_m.root");
  TTree* access;
  file->GetObject("SensorsTree",access);

  for(vector<string>::const_iterator sensor = sensors.begin(); sensor != sensors.end(); ++sensor)
    {
      cout << endl << endl << "--- Starting with Sensor " << (*sensor) << " ---" << endl << endl;
      TGraph *framem0 = 0;
      TGraph *framem1 = 0;
      TGraph *framem2 = 0;
      TGraph *framem3 = 0;
      TGraph *framem4 = 0;

      int entries = access->GetEntries("Timestamp<20130122000000000");
      int step = entries/_frames;
      
      ostringstream varname; varname << "CMS_CS_HF_MINUS_CASTOR_" << *sensor;

      

      if(_do_fits)
        {
          cout << "-- Fitting..." << endl;
          string comp1_start("20130108000000000");
          string comp1_end  ("20130109000000000");
          string comp2_start("20130117000000000");
          string comp2_end  ("20130118000000000");
          
          ostringstream comp1_selec, comp2_selec;
          comp1_selec << varname.str() << " < 100 && Timestamp >= " << comp1_start << " && Timestamp < " << comp1_end;
          comp2_selec << varname.str() << " < 100 && Timestamp >= " << comp2_start << " && Timestamp < " << comp2_end;

          TCanvas* c1 = new TCanvas;
          access->Draw((varname.str()+string(":Timestamp")).c_str(),comp1_selec.str().c_str());
          TGraph* comp1_hist = (TGraph*)(gROOT->FindObject("Graph")->Clone());
          
          TCanvas* c2 = new TCanvas;
          access->Draw((varname.str()+string(":Timestamp")).c_str(),comp2_selec.str().c_str());
          TGraph* comp2_hist = (TGraph*)gPad->GetPrimitive("Graph")->Clone(); //stupid root deletes graphs

          if(!comp1_hist || !comp1_hist->GetN()) throw runtime_error("comp1 empty");
          if(!comp2_hist || !comp2_hist->GetN()) throw runtime_error("comp2 empty");
          cout << "-Graphs drawn" << endl;

          TFitResultPtr comp1_fit = comp1_hist->Fit("pol0","S");
          TFitResultPtr comp2_fit = comp2_hist->Fit("pol0","S");

          cout << "- Fit convergencs status | comp1: " <<  int(comp1_fit) << " | comp2: " << int(comp2_fit) << endl;
          
          sensors_delta.push_back(comp2_fit->Parameter(0)-comp1_fit->Parameter(0));

          cout << "--...done fitting" << endl;

        }//if do fits



      if(_do_animation)
        {
          cout << "-- Animating..." << endl;
          access->Draw((varname.str()+string(">>hist")).c_str(),(varname.str()+string("<100 && Timestamp<20130122000000000")).c_str());
          TH1D* hist =  (TH1D*)gPad->GetPrimitive("hist");
          double ymax = hist->GetBinCenter(hist->GetNbinsX()-1)*1.2;
          double ymin = hist->GetBinCenter(0)*0.8;
      
          cout << endl << "Using for y-axis: " << ymin << "->" << ymax << endl;

          TH2D* histax = new TH2D("histax",";B [T];pos in mm",20,0,4,100,ymin,ymax);
          histax->SetStats(kFALSE);

          int namecount = 0;
          for (int frame=0; frame<entries; frame+=step)
            {
              if(framem3) framem4 = (TGraph*)framem3->Clone(); else framem4 = 0;
              if(framem2) framem3 = (TGraph*)framem2->Clone(); else framem3 = 0;
              if(framem1) framem2 = (TGraph*)framem1->Clone(); else framem2 = 0;
              if(framem0) framem1 = (TGraph*)framem0->Clone(); else framem1 = 0;

              ostringstream drawcmd,outfilename;
              outfilename << "pos_" << (*sensor) << "_" << setfill('0') << setw(5) << namecount++ << ".png";
              drawcmd << varname.str() << "<100 && Entry$ > " << frame << " && Entry$ <= " << (frame+step);

              cout << "Filename: " << outfilename.str() << endl;
              cout << "Entries: " << access->Draw((varname.str()+string(":MAGNET_FIELD_VALUE")).c_str(),drawcmd.str().c_str()) << endl;
      
              framem0 = (TGraph*)gPad->GetPrimitive("Graph");

              TCanvas* c = new TCanvas;
              histax->Draw();

              if(framem4)
                {
                  framem4->SetMarkerColor(kGray);
                  framem4->Draw("P");
                }
              if(framem3)
                {
                  framem3->SetMarkerColor(kGray+1);
                  framem3->Draw("P");
                }
              if(framem2)
                {
                  framem2->SetMarkerColor(kGray+2);
                  framem2->Draw("P");
                }
              if(framem1)
                {
                  framem1->SetMarkerColor(kGray+3);
                  framem1->Draw("P");
                }
              if(framem0 )framem0->Draw("P");
              else cout << "frame0m empty" << endl;

              c->SaveAs(outfilename.str().c_str());
            }//entry loop
          cout << "--...done animating (please run, e.g., \"convert pos_IP_NEAR_TOP_* pos_IP_NEAR_TOP.gif\"" << endl;
        }//if do animation
    }//sensor loop

  if(_do_fits)
    {
      assert(sensors.size() == sensors_delta.size());
      assert(sensors.size() == sensors_angle.size());
      for(int i = 0; i<int(sensors.size()); i++)
        {
          double r = sensors_delta[i];
          double alpha = sensors_angle[i];
          double x = r * cos(alpha);
          double y = r * sin(alpha);
          cout << "Sensor: " << sensors[i] << " Delta: " << sensors_delta[i] << " X: " << x << " Y: " << y << endl;
        }
    }
      
  
}//program
