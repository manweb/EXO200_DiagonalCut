void SigmaVsEnergy()
{
  int n = 4;
  
  double energy[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double sigma1[4] = {0.059052033*661.657, 0.035700989*1173.228, 0.031326622*1332.492, 0.018434475*2614.511};
  double sigma1_err[4] = {0.002555952*661.657, 0.001297748*1173.228, 0.000805116*1332.492, 0.000840202*2614.511};
  double sigma2[4] = {0.15095*661.657, 0.11069*1173.228, 0.09840*1173.228, 0.06832*2614.511};
  double sigma2_err[4] = {sigma2[0]*0.1, sigma2[1]*0.1, sigma2[2]*0.1, sigma2[3]*0.1};
  
  TGraphErrors *gr1 = new TGraphErrors(n,energy,sigma1,0,sigma1_err);
  TGraphErrors *gr2 = new TGraphErrors(n,energy,sigma2,0,sigma2_err);
  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.8);
  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.8);
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AZP");
  
  TCanvas *c2 = new TCanvas();
  gr2->Draw("AZP");
}