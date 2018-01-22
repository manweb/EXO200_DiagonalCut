const int nrBins = 200;
int binning = 1;
double binWidth = 5000.0*binning/nrBins;
double nSigma = 4.0;

double theta_ss = 0.1860;
double theta_ms = 0.2064;

double *peak_ss = new double[4];
double *peakRotated_ss = new double[4];
double *peakN_ss = new double[4];
double *peakNRotated_ss = new double[4];
double *cut_ss = new double[4];
double *cutRotated_ss = new double[4];
double *cutN_ss = new double[4];
double *cutNRotated_ss = new double[4];
double *peak_ms = new double[4];
double *peakRotated_ms = new double[4];
double *peakN_ms = new double[4];
double *peakNRotated_ms = new double[4];
double *cut_ms = new double[4];
double *cutRotated_ms = new double[4];
double *cutN_ms = new double[4];
double *cutNRotated_ms = new double[4];

void DiagonalCut_3()
{
  gStyle->SetPalette(100);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  Th228DiagonalCut(true);
  
  return;
  
  Co60DiagonalCut(false);
  Cs137DiagonalCut(false);
  
  TGraphErrors *grRotated_ss = new TGraphErrors(4,peakRotated_ss,cutRotated_ss,0,0);
  TGraphErrors *grRotated_ms = new TGraphErrors(4,peakRotated_ms,cutRotated_ms,0,0);
  TGraphErrors *gr_ss = new TGraphErrors(4,peak_ss,cut_ss,0,0);
  TGraphErrors *gr_ms = new TGraphErrors(4,peak_ms,cut_ms,0,0);
  
  TGraphErrors *grNRotated_ss = new TGraphErrors(4,peakRotated_ss,cutNRotated_ss,0,0);
  TGraphErrors *grNRotated_ms = new TGraphErrors(4,peakRotated_ms,cutNRotated_ms,0,0);
  TGraphErrors *grN_ss = new TGraphErrors(4,peakN_ss,cutN_ss,0,0);
  TGraphErrors *grN_ms = new TGraphErrors(4,peakN_ms,cutN_ms,0,0);
  
  grRotated_ss->SetTitle("Single Site");
  grRotated_ms->SetTitle("Multi Site");
  gr_ss->SetTitle("Single Site");
  gr_ms->SetTitle("Multi Site");
  
  grNRotated_ss->SetTitle("Single Site");
  grNRotated_ms->SetTitle("Multi Site");
  grN_ss->SetTitle("Single Site");
  grN_ms->SetTitle("Multi Site");
  
  grRotated_ss->SetMarkerStyle(20);
  grRotated_ss->SetMarkerSize(0.8);
  
  grRotated_ms->SetMarkerStyle(20);
  grRotated_ms->SetMarkerSize(0.8);
  
  gr_ss->SetMarkerStyle(20);
  gr_ss->SetMarkerSize(0.8);
  
  gr_ms->SetMarkerStyle(20);
  gr_ms->SetMarkerSize(0.8);
  
  grNRotated_ss->SetMarkerStyle(20);
  grNRotated_ss->SetMarkerSize(0.8);
  
  grNRotated_ms->SetMarkerStyle(20);
  grNRotated_ms->SetMarkerSize(0.8);
  
  grN_ss->SetMarkerStyle(20);
  grN_ss->SetMarkerSize(0.8);
  
  grN_ms->SetMarkerStyle(20);
  grN_ms->SetMarkerSize(0.8);
  
  grRotated_ss->GetXaxis()->SetLimits(0,4500);
  grRotated_ms->GetXaxis()->SetLimits(0,4500);
  gr_ss->GetXaxis()->SetLimits(0,3500);
  gr_ms->GetXaxis()->SetLimits(0,3500);
  
  grRotated_ss->GetYaxis()->SetRangeUser(0,12000);
  grRotated_ms->GetYaxis()->SetRangeUser(0,12000);
  gr_ss->GetYaxis()->SetRangeUser(0,12000);
  gr_ms->GetYaxis()->SetRangeUser(0,12000);
  
  TF1 *fit_ss = new TF1("fit_ss","pol2",0,5000);
  TF1 *fit_ms = new TF1("fit_ms","pol2",0,5000);
  TF1 *fitRotated_ss = new TF1("fitRotated_ss","pol2",0,5000);
  TF1 *fitRotated_ms = new TF1("fitRotated_ms","pol2",0,5000);
  
  TF1 *fitN_ss = new TF1("fitN_ss","pol2",0,5000);
  TF1 *fitN_ms = new TF1("fitN_ms","pol2",0,5000);
  TF1 *fitNRotated_ss = new TF1("fitNRotated_ss","pol2",0,5000);
  TF1 *fitNRotated_ms = new TF1("fitNRotated_ms","pol2",0,5000);
  
  fit_ss->SetLineWidth(2);
  fit_ss->SetLineColor(kRed);
  fit_ss->SetLineStyle(2);
  
  fit_ms->SetLineWidth(2);
  fit_ms->SetLineColor(kRed);
  fit_ms->SetLineStyle(2);
  
  fitRotated_ss->SetLineWidth(2);
  fitRotated_ss->SetLineColor(kRed);
  fitRotated_ss->SetLineStyle(2);
  
  fitRotated_ms->SetLineWidth(2);
  fitRotated_ms->SetLineColor(kRed);
  fitRotated_ms->SetLineStyle(2);
  
  fitN_ss->SetLineWidth(2);
  fitN_ss->SetLineColor(kRed);
  fitN_ss->SetLineStyle(2);
  
  fitN_ms->SetLineWidth(2);
  fitN_ms->SetLineColor(kRed);
  fitN_ms->SetLineStyle(2);
  
  fitNRotated_ss->SetLineWidth(2);
  fitNRotated_ss->SetLineColor(kRed);
  fitNRotated_ss->SetLineStyle(2);
  
  fitNRotated_ms->SetLineWidth(2);
  fitNRotated_ms->SetLineColor(kRed);
  fitNRotated_ms->SetLineStyle(2);
  
  gr_ss->Fit("fit_ss","r");
  gr_ms->Fit("fit_ms","r");
  grRotated_ss->Fit("fitRotated_ss","r");
  grRotated_ms->Fit("fitRotated_ms","r");
  
  grN_ss->Fit("fitN_ss","r");
  grN_ms->Fit("fitN_ms","r");
  grNRotated_ss->Fit("fitNRotated_ss","r");
  grNRotated_ms->Fit("fitNRotated_ms","r");
  
  TH2F *h_ss = new TH2F("h_ss","Single Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  TH2F *h_ms = new TH2F("h_ms","Multi Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  
  TH2F *hRotated_ss = new TH2F("hRotated_ss","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  TH2F *hRotated_ms = new TH2F("hRotated_ms","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  
  TChain *t = new TChain("t","tree");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2424_noReclustering_Fiducial.root");
  /*t->Add("../EnergyCalibration/analysis/alpha/V2/2426_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2431_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2432_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2433_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2434_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2447_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2448_noReclustering_Fiducial.root");*/
  
  t->Add("../EnergyCalibration/analysis/alpha/V2/2526_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2538_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2543_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2555_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2566_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2578_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2596_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2608_noReclustering_Fiducial.root");
  
  t->Draw("csc:epcrec >> h_ss","nsite == 1","goff");
  t->Draw("csc:epcrec >> h_ms","nsite > 1","goff");
  t->Draw("(csc*TMath::Cos(0.1860) - epcrec*TMath::Sin(0.1860)):(csc*TMath::Sin(0.1860) + epcrec*TMath::Cos(0.1860)) >> hRotated_ss","nsite == 1","goff");
  t->Draw("(csc*TMath::Cos(0.2064) - epcrec*TMath::Sin(0.2064)):(csc*TMath::Sin(0.2064) + epcrec*TMath::Cos(0.2064)) >> hRotated_ms","nsite > 1","goff");
  
  TCanvas *c1 = new TCanvas("c1","Diagonal cut",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  gr_ss->Draw("AZP");
  grN_ss->Draw("PZsame");
  c1->cd(2);
  gr_ms->Draw("AZP");
  grN_ms->Draw("ZPsame");
  
  TCanvas *c2 = new TCanvas("c2","Diagonal cut (rotated)",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  grRotated_ss->Draw("AZP");
  grNRotated_ss->Draw("ZPsame");
  c2->cd(2);
  grRotated_ms->Draw("AZP");
  grNRotated_ms->Draw("PZsame");
  
  TCanvas *c3 = new TCanvas("c3","2D spectrum",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  h_ss->Draw("colz");
  fit_ss->Draw("same");
  fitN_ss->Draw("same");
  c3->cd(2);
  h_ms->Draw("colz");
  fit_ms->Draw("same");
  fitN_ms->Draw("same");
  
  TCanvas *c4 = new TCanvas("c4","2D rotated spectrum",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  hRotated_ss->Draw("colz");
  fitRotated_ss->Draw("same");
  fitNRotated_ss->Draw("same");
  c4->cd(2);
  hRotated_ms->Draw("colz");
  fitRotated_ms->Draw("same");
  fitNRotated_ms->Draw("same");
  
  return;
}

void Th228DiagonalCut(bool disp)
{
  TH2F *h_ss = new TH2F("h_ss","Single Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  TH2F *h_ms = new TH2F("h_ms","Multi Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  
  TH2F *hRotated_ss = new TH2F("hRotated_ss","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  TH2F *hRotated_ms = new TH2F("hRotated_ms","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  
  TChain *t = new TChain("t","tree");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2424_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2426_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2431_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2432_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2433_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2434_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2447_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2448_noReclustering_Fiducial.root");
  
  t->Draw("csc:epcrec >> h_ss","nsite == 1","goff");
  t->Draw("csc:epcrec >> h_ms","nsite > 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ss",theta_ss,theta_ss,theta_ss,theta_ss),"nsite == 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ms",theta_ms,theta_ms,theta_ms,theta_ms),"nsite > 1","goff");
  
  int StartBin_ss = 4000.0/binWidth;
  int EndBin_ss = 4400.0/binWidth;
  int StartBin_ms = 4200.0/binWidth;
  int EndBin_ms = 4600.0/binWidth;
  
  TH1D *hPY_ss = hRotated_ss->ProjectionY("hPY_ss",StartBin_ss,EndBin_ss);
  TH1D *hPY_ms = hRotated_ms->ProjectionY("hPY_ms",StartBin_ms,EndBin_ms);
  TH1D *hPX_ss = hRotated_ss->ProjectionX("hPX_ss",int(7000.0/12000.0*nrBins),int(8500.0/12000.0*nrBins));
  TH1D *hPX_ms = hRotated_ms->ProjectionX("hPX_ms",int(7300.0/12000.0*nrBins),int(8800.0/12000.0*nrBins));
  
  TF1 *fitPY_ss = new TF1("fitPY_ss","gaus",6000,10000);
  fitPY_ss->SetParLimits(2,0,200);
  fitPY_ss->SetLineWidth(2);
  fitPY_ss->SetLineColor(kBlue);
  
  TF1 *fitPY_ms = new TF1("fitPY_ms","gaus",6000,10000);
  fitPY_ms->SetParLimits(2,0,200);
  fitPY_ms->SetLineWidth(2);
  fitPY_ms->SetLineColor(kBlue);
  
  TF1 *fitPX_ss = new TF1("fitPX_ss","gaus",4200,4400);
  fitPX_ss->SetParLimits(2,0,200);
  fitPX_ss->SetLineWidth(2);
  fitPX_ss->SetLineColor(kBlue);
  
  TF1 *fitPX_ms = new TF1("fitPX_ms","gaus",4200,4600);
  fitPX_ms->SetParLimits(2,0,200);
  fitPX_ms->SetLineWidth(2);
  fitPX_ms->SetLineColor(kBlue);
  
  hPY_ss->Fit("fitPY_ss","r");
  hPY_ms->Fit("fitPY_ms","r");
  hPX_ss->Fit("fitPX_ss","r");
  hPX_ms->Fit("fitPX_ms","r");
  
  double cutX_ss = 4258.49;
  double cutY_ss = fitPY_ss->GetParameter(1) + nSigma*fitPY_ss->GetParameter(2);
  double cutNY_ss = fitPY_ss->GetParameter(1) - nSigma*fitPY_ss->GetParameter(2);
  
  //peakRotated_ss[3] = fitPX_ss->GetParameter(1);
  peakRotated_ss[3] = cutX_ss;                      // <-- value from calibration
  cutRotated_ss[3] = cutY_ss;
  cutNRotated_ss[3] = cutNY_ss;
  peak_ss[3] = cutX_ss*TMath::Cos(theta_ms) - cutY_ss*TMath::Sin(theta_ms);
  peakN_ss[3] = cutX_ss*TMath::Cos(theta_ms) - cutNY_ss*TMath::Sin(theta_ms);
  cut_ss[3] = cutX_ss*TMath::Sin(theta_ms) + cutY_ss*TMath::Cos(theta_ms);
  cutN_ss[3] = cutX_ss*TMath::Sin(theta_ms) + cutNY_ss*TMath::Cos(theta_ms);
  
  double cutX_ms = 4402.39;
  double cutY_ms = fitPY_ms->GetParameter(1) + nSigma*fitPY_ms->GetParameter(2);
  double cutNY_ms = fitPY_ms->GetParameter(1) - nSigma*fitPY_ms->GetParameter(2);
  
  //peakRotated_ms[3] = fitPX_ms->GetParameter(1);
  peakRotated_ms[3] = cutX_ms;                      // <-- value from calibration
  cutRotated_ms[3] = cutY_ms;
  cutNRotated_ms[3] = cutNY_ms;
  peak_ms[3] = cutX_ms*TMath::Cos(theta_ms) - cutY_ms*TMath::Sin(theta_ms);
  peakN_ms[3] = cutX_ms*TMath::Cos(theta_ms) - cutNY_ms*TMath::Sin(theta_ms);
  cut_ms[3] = cutX_ms*TMath::Sin(theta_ms) + cutY_ms*TMath::Cos(theta_ms);
  cutN_ms[3] = cutX_ms*TMath::Sin(theta_ms) + cutNY_ms*TMath::Cos(theta_ms);
  
  if (!disp) {return;}
  
  TCanvas *c1 = new TCanvas("c1","2D spectrum",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h_ss->Draw("colz");
  c1->cd(2);
  h_ms->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","2D rotated spectrum",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  hRotated_ss->Draw("colz");
  c2->cd(2);
  hRotated_ms->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","Projection",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  hPY_ss->Draw("EZP");
  c3->cd(2);
  hPX_ss->Draw("EZP");
  
  TCanvas *c4 = new TCanvas("c4","Projection",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  hPY_ms->Draw("EZP");
  c4->cd(2);
  hPX_ms->Draw("EZP");
  
  return;
}

void Co60DiagonalCut(bool disp)
{
  TH2F *h_ss = new TH2F("h_ss","Single Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  TH2F *h_ms = new TH2F("h_ms","Multi Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  
  TH2F *hRotated_ss = new TH2F("hRotated_ss","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  TH2F *hRotated_ms = new TH2F("hRotated_ms","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  
  TChain *t = new TChain("t","tree");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2526_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2538_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2543_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2555_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2566_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2578_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2596_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2608_noReclustering_Fiducial.root");
  
  t->Draw("csc:epcrec >> h_ss","nsite == 1","goff");
  t->Draw("csc:epcrec >> h_ms","nsite > 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ss",theta_ss,theta_ss,theta_ss,theta_ss),"nsite == 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ms",theta_ms,theta_ms,theta_ms,theta_ms),"nsite > 1","goff");
  
  int StartBin1_ss = 1800.0/binWidth;
  int EndBin1_ss = 2000.0/binWidth;
  int StartBin1_ms = 1860.0/binWidth;
  int EndBin1_ms = 2120.0/binWidth;
  int StartBin2_ss = 2080.0/binWidth;
  int EndBin2_ss = 2280.0/binWidth;
  int StartBin2_ms = 2140.0/binWidth;
  int EndBin2_ms = 2380.0/binWidth;
  
  TH1D *hPY1_ss = hRotated_ss->ProjectionY("hPY1_ss",StartBin1_ss,EndBin1_ss);
  TH1D *hPY1_ms = hRotated_ms->ProjectionY("hPY1_ms",StartBin1_ms,EndBin1_ms);
  TH1D *hPY2_ss = hRotated_ss->ProjectionY("hPY2_ss",StartBin2_ss,EndBin2_ss);
  TH1D *hPY2_ms = hRotated_ms->ProjectionY("hPY2_ms",StartBin2_ms,EndBin2_ms);
  
  TF1 *fitPY1_ss = new TF1("fitPY1_ss","gaus",2000,6000);
  fitPY1_ss->SetParLimits(2,0,200);
  fitPY1_ss->SetLineWidth(2);
  fitPY1_ss->SetLineColor(kBlue);
  
  TF1 *fitPY1_ms = new TF1("fitPY1_ms","gaus",2000,6000);
  fitPY1_ms->SetParLimits(2,0,200);
  fitPY1_ms->SetLineWidth(2);
  fitPY1_ms->SetLineColor(kBlue);
  
  TF1 *fitPY2_ss = new TF1("fitPY2_ss","gaus",2000,6000);
  fitPY2_ss->SetParLimits(2,0,200);
  fitPY2_ss->SetLineWidth(2);
  fitPY2_ss->SetLineColor(kBlue);
  
  TF1 *fitPY2_ms = new TF1("fitPY2_ms","gaus",2000,6000);
  fitPY2_ms->SetParLimits(2,0,200);
  fitPY2_ms->SetLineWidth(2);
  fitPY2_ms->SetLineColor(kBlue);
  
  hPY1_ss->Fit("fitPY1_ss","r");
  hPY1_ms->Fit("fitPY1_ms","r");
  hPY2_ss->Fit("fitPY2_ss","r");
  hPY2_ms->Fit("fitPY2_ms","r");
  
  double cutX1_ss = 1917.52;
  double cutY1_ss = fitPY1_ss->GetParameter(1) + nSigma*fitPY1_ss->GetParameter(2);
  double cutNY1_ss = fitPY1_ss->GetParameter(1) - nSigma*fitPY1_ss->GetParameter(2);
  double cutX2_ss = 2179.38;
  double cutY2_ss = fitPY2_ss->GetParameter(1) + nSigma*fitPY2_ss->GetParameter(2);
  double cutNY2_ss = fitPY2_ss->GetParameter(1) - nSigma*fitPY2_ss->GetParameter(2);
  double cutX1_ms = 1984.67;
  double cutY1_ms = fitPY1_ms->GetParameter(1) + nSigma*fitPY1_ms->GetParameter(2);
  double cutNY1_ms = fitPY1_ms->GetParameter(1) - nSigma*fitPY1_ms->GetParameter(2);
  double cutX2_ms = 2255.21;
  double cutY2_ms = fitPY2_ms->GetParameter(1) + nSigma*fitPY2_ms->GetParameter(2);
  double cutNY2_ms = fitPY2_ms->GetParameter(1) - nSigma*fitPY2_ms->GetParameter(2);
  
  //peakRotated_ss[3] = fitPX_ss->GetParameter(1);
  peakRotated_ss[1] = cutX1_ss;                      // <-- value from calibration
  cutRotated_ss[1] = cutY1_ss;
  cutNRotated_ss[1] = cutNY1_ss;
  peak_ss[1] = cutX1_ss*TMath::Cos(theta_ms) - cutY1_ss*TMath::Sin(theta_ms);
  peakN_ss[1] = cutX1_ss*TMath::Cos(theta_ms) - cutNY1_ss*TMath::Sin(theta_ms);
  cut_ss[1] = cutX1_ss*TMath::Sin(theta_ms) + cutY1_ss*TMath::Cos(theta_ms);
  cutN_ss[1] = cutX1_ss*TMath::Sin(theta_ms) + cutNY1_ss*TMath::Cos(theta_ms);
  
  //peakRotated_ms[3] = fitPX_ms->GetParameter(1);
  peakRotated_ms[1] = cutX1_ms;                      // <-- value from calibration
  cutRotated_ms[1] = cutY1_ms;
  cutNRotated_ms[1] = cutNY1_ms;
  peak_ms[1] = cutX1_ms*TMath::Cos(theta_ms) - cutY1_ms*TMath::Sin(theta_ms);
  peakN_ms[1] = cutX1_ms*TMath::Cos(theta_ms) - cutNY1_ms*TMath::Sin(theta_ms);
  cut_ms[1] = cutX1_ms*TMath::Sin(theta_ms) + cutY1_ms*TMath::Cos(theta_ms);
  cutN_ms[1] = cutX1_ms*TMath::Sin(theta_ms) + cutNY1_ms*TMath::Cos(theta_ms);
  
  //peakRotated_ss[3] = fitPX_ss->GetParameter(1);
  peakRotated_ss[2] = cutX2_ss;                      // <-- value from calibration
  cutRotated_ss[2] = cutY2_ss;
  cutNRotated_ss[2] = cutNY2_ss;
  peak_ss[2] = cutX2_ss*TMath::Cos(theta_ms) - cutY2_ss*TMath::Sin(theta_ms);
  peakN_ss[2] = cutX2_ss*TMath::Cos(theta_ms) - cutNY2_ss*TMath::Sin(theta_ms);
  cut_ss[2] = cutX2_ss*TMath::Sin(theta_ms) + cutY2_ss*TMath::Cos(theta_ms);
  cutN_ss[2] = cutX2_ss*TMath::Sin(theta_ms) + cutNY2_ss*TMath::Cos(theta_ms);
  
  //peakRotated_ms[3] = fitPX_ms->GetParameter(1);
  peakRotated_ms[2] = cutX2_ms;                      // <-- value from calibration
  cutRotated_ms[2] = cutY2_ms;
  cutNRotated_ms[2] = cutNY2_ms;
  peak_ms[2] = cutX2_ms*TMath::Cos(theta_ms) - cutY2_ms*TMath::Sin(theta_ms);
  peakN_ms[2] = cutX2_ms*TMath::Cos(theta_ms) - cutNY2_ms*TMath::Sin(theta_ms);
  cut_ms[2] = cutX2_ms*TMath::Sin(theta_ms) + cutY2_ms*TMath::Cos(theta_ms);
  cutN_ms[2] = cutX2_ms*TMath::Sin(theta_ms) + cutNY2_ms*TMath::Cos(theta_ms);
  
  if (!disp) {return;}
  
  TCanvas *c1 = new TCanvas("c1","2D spectrum",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h_ss->Draw("colz");
  c1->cd(2);
  h_ms->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","2D rotated spectrum",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  hRotated_ss->Draw("colz");
  c2->cd(2);
  hRotated_ms->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","Projection",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  hPY1_ss->Draw("EZP");
  c3->cd(2);
  hPY2_ss->Draw("EZP");
  
  TCanvas *c4 = new TCanvas("c4","Projection",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  hPY1_ms->Draw("EZP");
  c4->cd(2);
  hPY2_ms->Draw("EZP");
  
  return;
}

void Cs137DiagonalCut(bool disp)
{
  TH2F *h_ss = new TH2F("h_ss","Single Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  TH2F *h_ms = new TH2F("h_ms","Multi Site Spectrum",nrBins,0,3500,nrBins,0,12000);
  
  TH2F *hRotated_ss = new TH2F("hRotated_ss","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  TH2F *hRotated_ms = new TH2F("hRotated_ms","Rotated Single Site Spectrum",nrBins,0,5000,nrBins,0,12000);
  
  TChain *t = new TChain("t","tree");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2450_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2469_noReclustering_Fiducial.root");
  t->Add("../EnergyCalibration/analysis/alpha/V2/2473_noReclustering_Fiducial.root");
  
  t->Draw("csc:epcrec >> h_ss","nsite == 1","goff");
  t->Draw("csc:epcrec >> h_ms","nsite > 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ss",theta_ss,theta_ss,theta_ss,theta_ss),"nsite == 1","goff");
  t->Draw(Form("(csc*TMath::Cos(%.4f) - epcrec*TMath::Sin(%.4f)):(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f)) >> hRotated_ms",theta_ms,theta_ms,theta_ms,theta_ms),"nsite > 1","goff");
  
  int StartBin_ss = 950.0/binWidth;
  int EndBin_ss = 1150.0/binWidth;
  int StartBin_ms = 960.0/binWidth;
  int EndBin_ms = 1260.0/binWidth;
  
  TH1D *hPY_ss = hRotated_ss->ProjectionY("hPY_ss",StartBin_ss,EndBin_ss);
  TH1D *hPY_ms = hRotated_ms->ProjectionY("hPY_ms",StartBin_ms,EndBin_ms);
  
  TF1 *fitPY_ss = new TF1("fitPY_ss","gaus",1000,4000);
  fitPY_ss->SetParLimits(2,0,200);
  fitPY_ss->SetLineWidth(2);
  fitPY_ss->SetLineColor(kBlue);
  
  TF1 *fitPY_ms = new TF1("fitPY_ms","gaus",1000,4000);
  fitPY_ms->SetParLimits(2,0,200);
  fitPY_ms->SetLineWidth(2);
  fitPY_ms->SetLineColor(kBlue);
  
  hPY_ss->Fit("fitPY_ss","r");
  hPY_ms->Fit("fitPY_ms","r");
  
  double cutX_ss = 1073.55;
  double cutY_ss = fitPY_ss->GetParameter(1) + nSigma*fitPY_ss->GetParameter(2);
  double cutNY_ss = fitPY_ss->GetParameter(1) - nSigma*fitPY_ss->GetParameter(2);
  
  //peakRotated_ss[3] = fitPX_ss->GetParameter(1);
  peakRotated_ss[0] = cutX_ss;                      // <-- value from calibration
  cutRotated_ss[0] = cutY_ss;
  cutNRotated_ss[0] = cutNY_ss;
  peak_ss[0] = cutX_ss*TMath::Cos(theta_ms) - cutY_ss*TMath::Sin(theta_ms);
  peakN_ss[0] = cutX_ss*TMath::Cos(theta_ms) - cutNY_ss*TMath::Sin(theta_ms);
  cut_ss[0] = cutX_ss*TMath::Sin(theta_ms) + cutY_ss*TMath::Cos(theta_ms);
  cutN_ss[0] = cutX_ss*TMath::Sin(theta_ms) + cutNY_ss*TMath::Cos(theta_ms);
  
  double cutX_ms = 1103.94;
  double cutY_ms = fitPY_ms->GetParameter(1) + nSigma*fitPY_ms->GetParameter(2);
  double cutNY_ms = fitPY_ms->GetParameter(1) - nSigma*fitPY_ms->GetParameter(2);
  
  //peakRotated_ms[3] = fitPX_ms->GetParameter(1);
  peakRotated_ms[0] = cutX_ms;                      // <-- value from calibration
  cutRotated_ms[0] = cutY_ms;
  cutNRotated_ms[0] = cutNY_ms;
  peak_ms[0] = cutX_ms*TMath::Cos(theta_ms) - cutY_ms*TMath::Sin(theta_ms);
  peakN_ms[0] = cutX_ms*TMath::Cos(theta_ms) - cutNY_ms*TMath::Sin(theta_ms);
  cut_ms[0] = cutX_ms*TMath::Sin(theta_ms) + cutY_ms*TMath::Cos(theta_ms);
  cutN_ms[0] = cutX_ms*TMath::Sin(theta_ms) + cutNY_ms*TMath::Cos(theta_ms);
  
  if (!disp) {return;}
  
  TCanvas *c1 = new TCanvas("c1","2D spectrum",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h_ss->Draw("colz");
  c1->cd(2);
  h_ms->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","2D rotated spectrum",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  hRotated_ss->Draw("colz");
  c2->cd(2);
  hRotated_ms->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","Projection",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  hPY_ss->Draw("EZP");
  //c3->cd(2);
  //hPX_ss->Draw("EZP");
  
  TCanvas *c4 = new TCanvas("c4","Projection",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  hPY_ms->Draw("EZP");
  //c4->cd(2);
  //hPX_ms->Draw("EZP");
  
  return;
}
