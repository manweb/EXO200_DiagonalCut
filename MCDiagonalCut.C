using namespace RooFit;

int nsite;
bool flFiducial;
bool flHasMissingPositions;
bool flCut;

double ClusterCull(double *Energy_cc, double *Energy_sc, EXOScintillationCluster *sc, EXOEventData* ED);

double Charge_p0_ss = 68.263;
double Charge_p1_ss = 0.9369;
double Charge_p0_ms = 113.4;
double Charge_p1_ms = 0.9382;

double Charge_p0_err_ss = 2.4;
double Charge_p1_err_ss = 0.0012;
double Charge_p0_err_ms = 3.3;
double Charge_p1_err_ms = 0.0022;

double Scint_p0_ss = -143.2;
double Scint_p1_ss = 0.3194;     
double Scint_p0_ms = -143.2;         // <--- need multi site calibration
double Scint_p1_ms = 0.3194;         // <--- need multi site calibration

double Scint_p0_err_ss = 7.627;
double Scint_p1_err_ss = 0.001773;
double Scint_p0_err_ms = 7.627;      // <--- need multi site calibration
double Scint_p1_err_ms = 0.001773;   // <--- need multi site calibration

double Rotated_theta_ss = 0.1860;
double Rotated_theta_ms = 0.2064;

double Rotated_p0_ss = 3.81;
double Rotated_p1_ss = 0.6114;
double Rotated_p0_ms = 7.6951;
double Rotated_p1_ms = 0.59;

double Rotated_p0_err_ss = 2.0292;
double Rotated_p1_err_ss = 0.001133;
double Rotated_p0_err_ms = 2.3016;
double Rotated_p1_err_ms = 0.001212;

double Rotated_rho_ss = -0.895;
double Rotated_rho_ms = -0.890;

double p0_DiagonalCut_ss = 2717.72;
double p1_DiagonalCut_ss = 3.4844;

double p0_DiagonalCut_ms = 2894.05;
double p1_DiagonalCut_ms = 3.7024;

double MC_p0_ss = -1.988;
double MC_p1_ss = 1.026;
double MC_p0_ms = 3.993;
double MC_p1_ms = 1.02;

double ChargeThreshold = 0.0;

double Calib_p1_ss = p1_DiagonalCut_ss*Scint_p1_ss/Charge_p1_ss;
double Calib_p0_ss = Scint_p1_ss*(p0_DiagonalCut_ss - p1_DiagonalCut_ss*Charge_p0_ss/Charge_p1_ss);
double Calib_p1_ms = p1_DiagonalCut_ms*Scint_p1_ms/Charge_p1_ms;
double Calib_p0_ms = Scint_p1_ms*(p0_DiagonalCut_ms - p1_DiagonalCut_ms*Charge_p0_ms/Charge_p1_ms);

TF2 *f2_ss;
TF2 *f2_ms;

void MCDiagonalCut()
{
  gROOT->SetStyle("Plain");

  f2_ss = new TF2("f2_ss",gaussian,-1200,1200,-1200,1200,6);
  f2_ss->SetParameters(1.0,0.0,0.0,45.0,170.0,-0.5044);
  
  f2_ms = new TF2("f2_ms",gaussian,-1200,1200,-1200,1200,6);
  f2_ms->SetParameters(1.0,0.0,0.0,50.0,190.0,-0.6054);
  
  f2_ss->SetNpx(500);
  f2_ss->SetNpy(500);
  
  f2_ms->SetNpx(500);
  f2_ms->SetNpy(500);

  char *fname = "/nfs/slac/g/exo_data2/exo_data/data/MC/Phase2a/P2_SourceP4_px_228Th/P2_SourceP4_px_228Th_*.root";
  
  EXOEventData *ED = 0;
  
  TChain treeMC("tree");
  treeMC.Add(fname);
  
  treeMC.SetBranchAddress("EventBranch", &ED);
  
  TH1F *h_ss = new TH1F("h_ss","h_ss",700,100,8000);
  TH1F *h_ms = new TH1F("h_ms","h_ms",700,100,8000);
  TH1F *h_ss_culled = new TH1F("h_ss_culled","h_ss_culled",700,100,8000);
  TH1F *h_ms_culled = new TH1F("h_ms_culled","h_ms_culled",700,100,8000);
  
  h_ss_culled->SetLineColor(kRed);
  h_ms_culled->SetLineColor(kRed);
  
  TH2F *h2_ss = new TH2F("h2_ss","h2_ss",200,0,3500,200,0,3500);
  TH2F *h2_ms = new TH2F("h2_ms","h2_ms",200,0,3500,200,0,3500);

  int count = 0;
  UInt_t nEvents = treeMC.GetEntries();
  // loop over events
  for (UInt_t evtID = 0; evtID < nEvents; evtID++){
    treeMC.GetEntry(evtID);
    
    if (evtID%1000 == 0) {cout << evtID << " events processed" << endl;}
    
    // only take events where the number of scintillation clusters = 1
    int nsc = ED->GetNumScintillationClusters();
    if (nsc!=1) {continue;}
    
    // Get the scintillation cluster and apply charge clustering
    double Energy_cc = 0.0;
    double Energy_sc = 0.0;
    double Etrue = 0.0;
    double Etrue_culled = 0.0;
    EXOScintillationCluster *sc = ED->GetScintillationCluster(0);
    ClusterCull(&Energy_cc, &Energy_sc, &Etrue, &Etrue_culled, sc, ED);  // where the truth designates that this applies to MC
    
    // Notice that we are culling clusters before energy calibration
    // For the time being, the culling does not apply an energy threshold.
    // In the future, because the simulation is so close to the true 
    // energy scale, it probably doesn't matter.
    
    // We only proceed if...
    if (!flFiducial) {continue;}      // the culled event is fiducial
    if (nsite == 0) {continue;}        // and if there are any charge clusters observed.
    if (flHasMissingPositions) {continue;}    // and all charge clusters have position reconstructed

    if (nsite == 1) {
      h2_ss->Fill(Energy_sc,Energy_cc);   // Fill histogram with MC energy
      h_ss->Fill(Etrue_culled);
      if (!flCut) {h_ss_culled->Fill(Etrue_culled);}
      count++;
    }
    else {
      h2_ms->Fill(Energy_sc,Energy_cc);   // Fill histogram with MC energy
      h_ms->Fill(Etrue_culled);
      if (!flCut) {h_ms_culled->Fill(Etrue_culled);}
    }
    if (count > 38255) {break;}
  }//loop through all events
  
  TF1 *DiagCut_ss = new TF1("DiagCut_ss","[1]*x + [0]",0,3500);
  DiagCut_ss->SetParameters(Calib_p0_ss,Calib_p1_ss);
  DiagCut_ss->SetLineWidth(2);
  DiagCut_ss->SetLineColor(kRed);
  DiagCut_ss->SetLineStyle(2);
  
  TF1 *DiagCut_ms = new TF1("DiagCut_ms","[1]*x + [0]",0,3500);
  DiagCut_ms->SetParameters(Calib_p0_ss,Calib_p1_ss);
  DiagCut_ms->SetLineWidth(2);
  DiagCut_ms->SetLineColor(kRed);
  DiagCut_ms->SetLineStyle(2);
  
  TFile *fOut = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/DiagonalCut/MCDiagCut.root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h2_ss->Draw("colz");
  DiagCut_ss->Draw("same");
  c1->cd(2);
  h2_ms->Draw("colz");
  DiagCut_ms->Draw("same");

  c1->Write();
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  h_ss->Draw();
  h_ss_culled->Draw("same");
  c2->cd(2);
  h_ms->Draw();
  h_ms_culled->Draw("same");

  c2->Write();

  fOut->Close();
  
  return;
}

double ClusterCull(double *Energy_cc, double *Energy_sc, double *Et, double *Et_culled, EXOScintillationCluster *sc, EXOEventData* ED)
{
  // Initialize variables
  nsite = 0;
  flFiducial = true;
  flHasMissingPositions = false;
  flCut = false;
  Int_t ncl = sc->GetNumChargeClusters();
  
  double Etrue = 0.0;
  double Etrue_culled = 0.0;
  
  EXOMonteCarloData *mc = 0;
  
  mc = (EXOMonteCarloData*)ED->fMonteCarloData;
  Int_t nPixelatedCharge = mc->GetNumPixelatedChargeDeposits();
  EXOMCPixelatedChargeDeposit* cd = 0;
  Int_t nU = ED->GetNumUWireSignals();
  EXOChargeCluster* cc = 0;
  EXOUWireSignal* ws = 0;
  Double_t R = 0;
  
  // sum up all pixelated charge deposits on the u wires of all clusters, provided cluster is not in the dead layer
  for (Int_t i = 0; i < nPixelatedCharge; i++) {
    cd = (EXOMCPixelatedChargeDeposit*)mc->GetPixelatedChargeDeposit(i);
    R = TMath::Sqrt(pow(cd->GetPixelCenter().GetX(),2.)+pow(cd->GetPixelCenter().GetY(),2.));//let's hope these are mm.
    {
      Int_t chn = cd->fDepositChannel;
      for(Int_t l=0;l<nU;++l) {
        ws = ED->GetUWireSignal(l);
        if(ws->fChannel == chn) {
          if(TMath::Abs(ws->fTime - cd->fWireHitTime) > 3500)//may be >1 U-wire signals with the same chn
            continue;
          cc = ws->GetChargeCluster();
          if(!cc) 
            continue;//GetChargeCluster may return NULL
          if(cc->GetScintillationCluster() != sc)
            break;
          if ((cc->fRawEnergy * MC_p1_ss + MC_p0_ss) < ChargeThreshold)//should we apply threshold to individual U-wire energies instead? -I.O.
            break;
          //found a wire affected by i-th pixel that belongs to charge cluster above the threshold 
          //add energy on this wire due to this pixel
          if (R <= 168.0) {Etrue_culled += cd->fTotalEnergy*1000.0;}
          Etrue += cd->fTotalEnergy*1000.0;
          break;
        };
      };//loop through u wire signals in this event in the search for the one with the right channel number 
    };//loop through u wire signals affected by this cluster
  };//loop through pixels
  
  // first step is to loop over charge clusters and make sure each of them is inside the active 
  // volume and is above energy threshold
  for (int clID = 0; clID < ncl; clID++){
    cc = sc->GetChargeClusterAt(clID);
    
    R = TMath::Sqrt(cc->fX * cc->fX + cc->fY * cc->fY);
    
    // Check if all charge clusters have complete position reconstruction. If not, set a flag in order
    // to deal with these events later
    if (cc->fX == -999 || cc->fY == -999) {flHasMissingPositions = true;}

    if (R > 168.0 || TMath::Abs(cc->fZ) > 182.0 || TMath::Abs(cc->fZ) < 5.0) {flFiducial = false;}
    
    // increment charge cluster count
    nsite++;
  }
  
  double p0_ion = 35.67;
  double p1_ion = 0.00469;
  double p0_scint = 70.48;
  double p1_scint = 0.03911;
  double theta = -0.5044;
  
  double sigma_ion = Etrue*p1_ion + p0_ion;
  double sigma_scint = Etrue*p1_scint + p0_scint;
  
  //TF2 *f2_ss = new TF2("f2_ss",gaussian,-400,400,-400,400,6);
  //f2_ss->SetParameters(1.0,0.0,0.0,sigma_ion,sigma_scint,theta);
  //f2_ss->SetNpx(500);
  //f2_ss->SetNpy(500);
  
  double r1;
  double r2;
  
  if (nsite == 1) {f2_ss->GetRandom2(r1,r2);}
  else {f2_ms->GetRandom2(r1,r2);}
  
  *Energy_cc = Etrue + r1 - (Etrue - Etrue_culled);
  *Energy_sc = Etrue + r2;
  *Et = Etrue;
  *Et_culled = Etrue_culled;
  
  // apply diagonal cut
  if (nsite == 1 && *Energy_sc > *Energy_cc*Calib_p1_ss + Calib_p0_ss) {flCut = true;}
  if (nsite > 1 && *Energy_sc > *Energy_cc*Calib_p1_ms + Calib_p0_ms) {flCut = true;}
  
  return Etrue;
}

void test()
{
  TH2F *h2 = new TH2F("h2","h2",400,0,3500,400,0,3500);
  //TH2F *h2 = new TH2F("h2","h2",400,-5,5,400,-5,5);
  
  TH1F *h = new TH1F("h","h",400,0,5000);
  
  TF2 *f2_ss = new TF2("f2_ss",gaussian,-800,800,-800,800,6);
  f2_ss->SetParameters(1.0,0.0,0.0,71.2248,186.476,-0.5044);
  f2_ss->SetNpx(500);
  f2_ss->SetNpy(500);
  
  TRandom3 rand(0);
  
  double theta = 0.5044;
  
  //double rho = -1.0*TMath::Tan(theta);
  double rho = -0.99;
  //double rho = TMath::Sin(2*theta)/(4*71.2248**2) + TMath::Sin(2*theta)/(4*186.476**2);
  //cout << rho << endl;
  
  double sigma_x_true = 71.2248;
  double sigma_y_true = 186.476;
  double sigma_x_corr = sigma_x_true*TMath::Cos(theta);
  double sigma_y_corr = sigma_y_true*TMath::Cos(theta);
  
  for (int i = 0; i < 50000; i++) {
    /*double u1 = rand.Gaus();
    double u2 = rand.Gaus();
    double r1 = sigma_x_corr*u1;
    double r2 = sigma_y_corr*(u1*rho + TMath::Sqrt(1-rho**2)*u2);*/
    double r1;
    double r2;
    f2_ss->GetRandom2(r1,r2);
    h2->Fill(2615.0+r1,2615.0+r2);
    h->Fill((2615.0+r2)*TMath::Sin(theta) + (2615.0+r1)*TMath::Cos(theta));
    //h2->Fill(r1,r2);
  }
  
  /*TMatrixD parmat(2,2);
  parmat(0,0) = 71.2248**2;
  parmat(1,1) = 186.476**2;
  parmat(0,1) = 71.2248*186.476*rho;
  parmat(1,0) = parmat(0,1);
  
  TDecompChol chol(parmat);
  if(!chol.Decompose()) {
    cout<<"Cholesky decomposition failed for cor1mat. Poorly conditioned matrix."<<endl;
    exit(1);
  };
  
  TMatrixD cholT = chol.GetU();
  cholT.T();
  
  TVectorD parset(2);
  
  for (int i = 0; i < 50000; i++) {
    parset(0) = rand.Gaus();
    parset(1) = rand.Gaus();
  
    parset *= cholT;
    
    h2->Fill(2615.0+parset(0),2615.0+parset(1));
  }*/
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  h2->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  h->Draw();
}

void test2()
{
  TH2F *h2 = new TH2F("h2","h2",100,-5,5,100,-5,5);
  
  RooRealVar theta("theta","theta",-0.5044);
  RooRealVar sigma_x("sigma_x","sigma_x",0.5);
  RooRealVar sigma_y("sigma_y","sigma_y",0.8);
  RooRealVar x("x","x",-5,5);
  RooRealVar y("y","y",-5,5);
  
  RooGenericPdf pdf("pdf","pdf","exp(-((cos(@4)**2/(2*@1**2)+sin(@4)**2/(2*@3**2))*@0**2 + (sin(@4)**2/(2*@1**2)+cos(@4)**2/(2*@3**2))*@2**2 + 2*(-sin(2*@4)/(4*@1**2)+sin(2*@4)/(4*@3**2))*@0*@2))",RooArgList(x,sigma_x,y,sigma_y,theta));
  
  TH2F *h2_1 = (TH2F*)pdf.createHistogram("h2_1",x,Binning(100),YVar(y,Binning(100)));
  
  TRandom3 rand(0);
  
  for (int i = 0; i < 1000; i++) {
    //double r = rand.Gaus();
    //sigma_x.setVal(r);
    RooDataSet *data = (RooDataSet*)pdf.generate(RooArgList(x,y),1);
    double val_x = ((RooRealVar*)(data->get())->find("x"))->getVal();
    double val_y = ((RooRealVar*)(data->get())->find("y"))->getVal();
    h2->Fill(val_x,val_y);
    //cout << data->getObservables()->Print("x") << endl;
  }
  
  h2->Draw("colz");
  
  //cout << x.getVal() << endl;
  
  //((TTree*)data->tree())->Draw("x");
}

double gaussian(double *x, double *par) {
  double A = par[0];
  double x_0 = par[1];
  double y_0 = par[2];
  double sigma_x = par[3];
  double sigma_y = par[4];
  double a = TMath::Cos(par[5])**2 / (2*sigma_x**2) + TMath::Sin(par[5])**2 / (2*sigma_y**2);
  double b = -1.0 * TMath::Sin(2*par[5]) / (4*sigma_x**2) + TMath::Sin(2*par[5]) / (4*sigma_y**2);
  double c = TMath::Sin(par[5])**2 / (2*sigma_x**2) + TMath::Cos(par[5])**2 / (2*sigma_y**2);
  
  double val = A*TMath::Exp(-1.0*(a*(x[0] - x_0)**2 + 2*b*(x[0] - x_0)*(x[1] - y_0) + c*(x[1] - y_0)**2));
  
  return val;
}
