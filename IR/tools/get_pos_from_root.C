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

void get_pos_from_root()
{
  vector<string> sensors;
  sensors.push_back("IP_FAR_TOP");
  sensors.push_back("IP_FAR_BOTTOM");
  sensors.push_back("IP_NEAR_TOP");
  sensors.push_back("IP_NEAR_BOTTOM");

  sensors.push_back("NONIP_FAR_TOP");
  sensors.push_back("NONIP_FAR_CENTER");
  sensors.push_back("NONIP_FAR_BOTTOM");
  sensors.push_back("NONIP_NEAR_TOP");
  sensors.push_back("NONIP_NEAR_CENTER");
  sensors.push_back("NONIP_NEAR_BOTTOM");

  vector<double> sensors_delta;

  TFile* file = TFile::Open("data/InputSensorNamesHFminus_v1-data_2013_time_fix_with_oracle_db_2_with_lumi_m.root");
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
          string comp_before_start("20130108000000000");
          string comp_before_end  ("20130109000000000");
          string comp_after_start("20130117000000000");
          string comp_after_end  ("20130118000000000");
          
          ostringstream comp_before_selec, comp_after_selec;
          comp_before_selec << varname.str() << " < 100 && Timestamp >= " << comp_before_start << " && Timestamp < " << comp_before_end;
          comp_after_selec << varname.str() << " < 100 && Timestamp >= " << comp_after_start << " && Timestamp < " << comp_after_end;

          TCanvas* c1 = new TCanvas;
          access->Draw((varname.str()+string(":Timestamp")).c_str(),comp_before_selec.str().c_str());
          TGraph* comp_before_hist = (TGraph*)(gROOT->FindObject("Graph")->Clone());
          
          TCanvas* c2 = new TCanvas;
          access->Draw((varname.str()+string(":Timestamp")).c_str(),comp_after_selec.str().c_str());
          TGraph* comp_after_hist = (TGraph*)gPad->GetPrimitive("Graph")->Clone(); //stupid root deletes graphs

          if(!comp_before_hist || !comp_before_hist->GetN()) throw runtime_error("comp_before empty");
          if(!comp_after_hist || !comp_after_hist->GetN()) throw runtime_error("comp_after empty");
          cout << "-Graphs drawn" << endl;

          TFitResultPtr comp_before_fit = comp_before_hist->Fit("pol0","SQ");
          TFitResultPtr comp_after_fit = comp_after_hist->Fit("pol0","SQ");

          cout << "\n-------------\n";
          //          cout << "- Fit convergencs status | before B field: " <<  int(comp_before_fit) << " | after B field: " << int(comp_after_fit) << endl;
          cout << "- Mean+-RMS before B field | after B field: " <<  comp_before_hist->GetMean(2) << "+-" << comp_before_hist->GetRMS(2) << " | " << comp_after_hist->GetMean(2) << "+-" << comp_after_hist->GetRMS(2) << endl;
          cout << "-------------\n\n\n";
          
          sensors_delta.push_back(comp_after_fit->Parameter(0)-comp_before_fit->Parameter(0));

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

}//program
