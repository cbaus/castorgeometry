#ifndef _CASTOR_H_
#define _CASTOR_H_

#include <iostream>
#include <iomanip>
#include <sstream>

#include "Point.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TGeoManager.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TVector.h"

const double castorLength = 15868. - 14390; //from begin to last rib


namespace myFit
{
  void FitFunction(int& /*npar*/, double* const /*grad*/,
                   double& chi2, double* const shift,
                   const int /*iFlag*/);
}

class CastorHalf {
 public:
 CastorHalf(bool farNotNear): f_theta(0)
    {
      f_center = new Point(0.,0.,0.);
      f_target1 = new Point(0.,0.,0.);
      f_target2 = new Point(0.,0.,0.);
      f_target3 = new Point(0.,0.,0.);
      f_far=farNotNear;

      if(f_far)
        {
          f_targetn1 = new Point(-90.96 , +317.21, castorLength);
          f_targetn2 = new Point(-317.21, +90.96,  castorLength-1010.);
          f_targetn3 = new Point(-317.21, -90.96,  castorLength);
        }
      else //if near
        {
          f_targetn1 = new Point(+92.90, +316.38, castorLength);
          f_targetn2 = new Point(+316.38, +92.90, castorLength-1010.);
          f_targetn3 = new Point(+316.38, -92.90, castorLength);
        }
    }
  ~CastorHalf()
    {
      delete f_center;
      delete f_target1;
      delete f_target2;
      delete f_target3;
    }
  void Draw(TGeoManager* geoM);
  void DrawXY(Color_t col);
  double GetTheta() {return f_theta;}
  double GetRho() {return f_rho;}
  Point* GetCenter() {return f_center;}
  Point* GetTarget1() {return new Point(*f_center+*f_target1);}
  Point* GetTarget2() {return new Point(*f_center+*f_target2);}
  Point* GetTarget3() {return new Point(*f_center+*f_target3);}
  void SetRelTarget1(Point p) {*f_target1 += p;}
  void SetRelTarget2(Point p) {*f_target2 += p;}
  void SetRelTarget3(Point p) {*f_target3 += p;}
  void Rotate(double theta) {f_theta = theta;}
  bool Fit(bool xyonly);
  void Translate(double x, double y, double z);
 private:
  Point* f_center;
  double f_theta;
  double f_rho;
  Point* f_target1;
  Point* f_target2;
  Point* f_target3;
  Point* f_targetn1;
  Point* f_targetn2;
  Point* f_targetn3;
  bool f_far;

};

void CastorHalf::Translate(double x, double y, double z)
{
  f_center->Translate(x,y,z);
  return;
}

void CastorHalf::Draw(TGeoManager* geoM)
{
  double angleBegin = 0;
  double angleEnd = 0;

  if(f_far)
    {
      angleBegin = 90;
      angleEnd = 270;
    }
  else
    {
      angleBegin = -90;
      angleEnd = 90;
    }

  std::ostringstream name;
  if (f_far)
    name << "geoFarHalf";
  else
    name << "geoNearHalf";
  TGeoMaterial* mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium* med = new TGeoMedium("Vacuum",1,mat);
  TGeoVolume* geoHalf = geoM->MakeTubs(name.str().c_str(), med, 200,290,castorLength/2.,angleBegin,angleEnd);

  TGeoTranslation t0(0,0,castorLength/2.); //shift because centre for turning around y-axis is in middle of CASTOR
  TGeoTranslation t0i(0,0,-castorLength/2.);
  TGeoTranslation t1(f_center->GetX(),f_center->GetY(),0);

  //std::cout << "Shifting by " << f_center->GetX() << " " << f_center->GetY() << std::endl;
  TGeoRotation r1;
  r1.RotateY(f_theta);
  TGeoHMatrix pos;
  pos = t1 * t0i * r1 * t0 ;
  TGeoHMatrix* posP = new TGeoHMatrix(pos);
  TGeoVolume* geoW = geoM->GetMasterVolume();
  //if(f_far)
  geoW->AddNode(geoHalf,1,posP);
  // else
  // geoW->AddNode(geoHalf,1);

  //Sensors
  std::ostringstream nametar;
  nametar.str("geoTarget1");
  if(f_far)
    nametar << "_far";
  else
    nametar << "_near";
  TGeoVolume* geoTarget1 = geoM->MakeSphere(nametar.str().c_str(), med, 0,5);
  TGeoTranslation tar1(f_targetn1->GetX(),f_targetn1->GetY(),f_targetn1->GetZ()-castorLength/2.);
  TGeoHMatrix* posT1 = new TGeoHMatrix(tar1);
  geoW->AddNode(geoTarget1,1,posT1);

  nametar.str("geoTarget2");
  if(f_far)
    nametar << "_far";
  else
    nametar << "_near";
  TGeoVolume* geoTarget2 = geoM->MakeSphere(nametar.str().c_str(), med, 0,5);
  TGeoTranslation tar2(f_targetn2->GetX(),f_targetn2->GetY(),f_targetn2->GetZ()-castorLength/2.);
  TGeoHMatrix* posT2 = new TGeoHMatrix(tar2);
  geoW->AddNode(geoTarget2,1,posT2);

  nametar.str("geoTarget3");
  if(f_far)
    nametar << "_far";
  else
    nametar << "_near";
  TGeoVolume* geoTarget3 = geoM->MakeSphere(nametar.str().c_str(), med, 0,5);
  TGeoTranslation tar3(f_targetn3->GetX(),f_targetn3->GetY(),f_targetn3->GetZ()-castorLength/2.);
  TGeoHMatrix* posT3 = new TGeoHMatrix(tar3);
  geoW->AddNode(geoTarget3,1,posT3);


}

void CastorHalf::DrawXY(Color_t col)
{
  TMarker* m = new TMarker(f_center->GetX(),f_center->GetY(),20);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");


  m = new TMarker(f_target1->GetX(),f_target1->GetY(),24);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  m = new TMarker(f_target2->GetX(),f_target2->GetY(),24);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  m = new TMarker(f_target3->GetX(),f_target3->GetY(),24);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  //nominal as crosses
  Point* targetn1 = new Point(f_targetn1); //copy for drawing. delete later
  Point* targetn2 = new Point(f_targetn2);
  Point* targetn3 = new Point(f_targetn3);
  targetn1->RotateX(f_theta);
  targetn1->RotateY(f_rho);
  targetn1->Translate(f_center->GetX(),f_center->GetY(),0);
  targetn2->RotateX(f_theta);
  targetn2->RotateY(f_rho);
  targetn2->Translate(f_center->GetX(),f_center->GetY(),0);
  targetn3->RotateX(f_theta);
  targetn3->RotateY(f_rho);
  targetn3->Translate(f_center->GetX(),f_center->GetY(),0);

  m = new TMarker(targetn1->GetX(),targetn1->GetY(),5);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  m = new TMarker(targetn2->GetX(),targetn2->GetY(),5);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  m = new TMarker(targetn3->GetX(),targetn3->GetY(),5);
  m->SetMarkerSize(2);
  m->SetMarkerColor(col);
  m->Draw("SAME");

  delete targetn1, targetn2, targetn3;


}
///Fit will move center to best possible position
bool CastorHalf::Fit(bool xyonly)
{
  const int param = xyonly?2:4;
  const int ndf = 9 - param;
  TMinuit theFitter(param);
  TVector targets(18);
  targets[0]=f_target1->GetX();
  targets[1]=f_target2->GetX();
  targets[2]=f_target3->GetX();
  targets[3]=f_target1->GetY();
  targets[4]=f_target2->GetY();
  targets[5]=f_target3->GetY();
  targets[6]=f_target1->GetZ();
  targets[7]=f_target2->GetZ();
  targets[8]=f_target3->GetZ();

  targets[ 9]=f_targetn1->GetX();
  targets[10]=f_targetn2->GetX();
  targets[11]=f_targetn3->GetX();
  targets[12]=f_targetn1->GetY();
  targets[13]=f_targetn2->GetY();
  targets[14]=f_targetn3->GetY();
  targets[15]=f_targetn1->GetZ();
  targets[16]=f_targetn2->GetZ();
  targets[17]=f_targetn3->GetZ();

  theFitter.SetFCN(myFit::FitFunction);
  theFitter.SetObjectFit(&targets);

  //INIT MISSING

  int ierflag;
  theFitter.mnparm(0,"x", 0,1,0,0,ierflag);
  theFitter.mnparm(1,"y", 0,1,0,0,ierflag);
  if(param == 4)
    {
      theFitter.mnparm(2,"theta", 0,0.1,-1,1,ierflag);
      theFitter.mnparm(3,"rho", 0,0.1,-1,1,ierflag);
    }

  double arglist[2]={100000, 0.1}; //calls, tolerance default 0.1
  theFitter.mnexcm("MIGRAD", arglist, 2, ierflag); //arglist has 2 params

  if (ierflag)
    {
      std::cerr << " MIGRAD failed. Error: " << ierflag << std::endl;
      return false;
    }

  theFitter.mnexcm("HESSE", arglist, 2, ierflag); //3rd option= 2 arg in list

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  theFitter.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  int status = icstat;
  double chi2 = amin;

  double x=0,y=0,theta=0,rho=0;
  double xe=0,ye=0,thetae=0,rhoe=0;
  theFitter.GetParameter(0,x,xe);
  theFitter.GetParameter(1,y,ye);
  if(param == 4)
    {
      theFitter.GetParameter(2,theta,thetae);
      theFitter.GetParameter(3,rho,rhoe);
    }

  Point shift(x,y,0);
  *f_center += shift;
  f_theta=theta;
  f_rho=rho;

  std::cout << "chi2/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf << std::endl;
  std::cout << "x = " << f_center->GetX() << std::endl;
  std::cout << "y = " << f_center->GetY() << std::endl;
  std::cout << "theta = " << theta / TMath::DegToRad() << "°" << std::endl;
  std::cout << "rho = " << rho / TMath::DegToRad() << "°" << std::endl << std::endl;

  /* TCanvas can; */
  /* ((TGraph*)theFitter.Contour(100))->Draw("AL"); */
  /* theFitter.SetErrorDef(2); */
  /* ((TGraph*)theFitter.Contour(100))->Draw("L"); */
  /* theFitter.SetErrorDef(3); */
  /* ((TGraph*)theFitter.Contour(100))->Draw("L"); */
  /* can.SaveAs("test.pdf"); */

  return true;
}
#endif //#ifndef _CASTOR_H_

void myFit::FitFunction(int& npar, double* const /*grad*/,
                        double& chi2, double* const shift,
                        const int /*iFlag*/)
{
  const double x = shift[0];
  const double y = shift[1];
  double theta = 0;
  double rho = 0;

  if(npar > 2)
    {
      theta = shift[2];
      rho = shift[3];
    }

  TVector* targets = (TVector*) gMinuit->GetObjectFit();

  Point* t1 = new Point((*targets)[0],(*targets)[3],(*targets)[15]);
  Point* t2 = new Point((*targets)[1],(*targets)[4],(*targets)[16]);
  Point* t3 = new Point((*targets)[2],(*targets)[5],(*targets)[17]);

  Point* tn1 = new Point((*targets)[ 9],(*targets)[12],(*targets)[15]); //FIXME
  Point* tn2 = new Point((*targets)[10],(*targets)[13],(*targets)[16]);
  Point* tn3 = new Point((*targets)[11],(*targets)[14],(*targets)[17]);

  if(npar>2)
    {
      tn1->RotateX(theta);
      tn2->RotateX(theta);
      tn3->RotateX(theta);
      tn1->RotateY(rho);
      tn2->RotateY(rho);
      tn3->RotateY(rho);
    }

  tn1->Translate(x,y,0);
  tn2->Translate(x,y,0);
  tn3->Translate(x,y,0);

  const double r1 = t1->GetDistance(tn1);
  const double r2 = t2->GetDistance(tn2);
  const double r3 = t3->GetDistance(tn3);

  const double r1e = t1->GetDistanceError(tn1);
  const double r2e = t2->GetDistanceError(tn2);
  const double r3e = t3->GetDistanceError(tn3);

  chi2 = 0;
  chi2 += pow(r1 / r1e,2);
  chi2 += pow(r2 / r2e,2); //mid point at castor front
  chi2 += pow(r3 / r3e,2);

  /* std::cerr << std::fixed << std::setprecision(8) << "x=" << x << "  y=" << y << " th=" << theta/TMath::DegToRad() << " rh=" << rho/TMath::DegToRad() */
  /*           << "   chi2=" << chi2 << "   r1=" << r1 << " r2=" << r2 << " r3=" << r3 << "   r1e=" << r1e << " r2e=" << r2e << std::endl */
  /*           << "  dx1=" << t1->GetX() << "/" << tn1->GetX() << " dx2=" << t2->GetX() << "/" << tn2->GetX() << "  dx3=" << t3->GetX() << "/" << tn3->GetX() << std::endl */
  /*           << "  dy1=" << t1->GetY() << "/" << tn1->GetY() << " dy2=" << t2->GetY() << "/" << tn2->GetY() << "  dy3=" << t3->GetY() << "/" << tn3->GetY() << std::endl */
  /*           << "  dz1=" << t1->GetZ() << "/" << tn1->GetZ() << " dz2=" << t2->GetZ() << "/" << tn2->GetZ() << "  dz3=" << t3->GetZ() << "/" << tn3->GetZ() << std::endl */
  /*           << std::endl; */

  delete t1,t2,t3;
  delete tn1,tn2,tn3;

  return;
}
