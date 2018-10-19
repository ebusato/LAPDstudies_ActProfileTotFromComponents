
class Beam {
public:
  double m_tp; // in sec
  double m_ts; // in sec
  double m_I; // nA of proton beam

  double getProtonsPerSecondAverage();
  double getProtonsPerSecondInSpill();
  double getSpillFraction();
};

double Beam::getSpillFraction()
{
  return m_ts/(m_ts+m_tp);
}

double Beam::getProtonsPerSecondAverage()
{
  return m_I*1e-9/1.6e-19;
}

double Beam::getProtonsPerSecondInSpill()
{
  return m_I*1e-9/1.6e-19*(m_tp+m_ts)/m_ts;
}

class Isotope {
public:
  double m_lambda; // sec^{-1}
  double m_rate; // # beta+ emitters per proton
  double m_rateprime; // # beta+ emitters per sec in spill
  Beam* p_beamPtr;
  
  Isotope(double lambda, double r, Beam* b);
  double NumberOfDecaysAfterActivation(double N0, double time);
  double NumberOfDecaysAfterDeactivation(double N0, double time);
  double TotNumber(double N0, double time);
};

Isotope::Isotope(double l, double r, Beam* b)
{
  m_lambda = l;
  m_rate = r;
  m_rateprime = m_rate*b->getProtonsPerSecondInSpill();
  cout << "rates: " << m_rate << "  " << m_rateprime << endl;
  p_beamPtr = b; 
}

// NO = initial number (before start of activation phase)
// time (in sec) = duration of activation period
double Isotope::NumberOfDecaysAfterActivation(double N0, double time)
{
  return N0+m_rateprime*p_beamPtr->getSpillFraction()*(time+(1/m_lambda)*(TMath::Exp(-1*m_lambda*time)-1));
}

// N0 = initial number (before start of deactivation phase)
// time (in sec) = duration of deactivation period
double Isotope::NumberOfDecaysAfterDeactivation(double N0, double time)
{
  return N0-N0*TMath::Exp(-1*m_lambda*time);
}

// time in sec
double Isotope::TotNumber(double N0, double time)
{
  return N0 + p_beamPtr->getProtonsPerSecondAverage()*m_rate*time;
}

void setGraphStyle(TGraph* g, int color)
{
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(2);
}

TH1* scaleHisto(TH1* h, int n)
{
  TString nameClone(h->GetName());
  nameClone+="_clone";
  TH1* hClone = (TH1*) h->Clone(nameClone.Data());
  hClone->Scale(n/hClone->Integral());
  return hClone;
}

void DrawStack(int idx, TH1* h_N12, std::vector<double> N12numberAt, TH1* h_C11, std::vector<double> C11numberAt,TH1* h_C10, std::vector<double> C10numberAt, TH1* h_O15, std::vector<double> O15numberAt, TH1* h_Gamma, std::vector<double> GammanumberAt)
{
  TH1* h_N12Clone = scaleHisto(h_N12, N12numberAt[idx]);
  TH1* h_C11Clone = scaleHisto(h_C11, C11numberAt[idx]);
  TH1* h_C10Clone = scaleHisto(h_C10, C10numberAt[idx]);
  TH1* h_O15Clone = scaleHisto(h_O15, O15numberAt[idx]);
  TH1* h_GammaClone = scaleHisto(h_Gamma, GammanumberAt[idx]);
  
  THStack* hStack = new THStack();
  hStack->Add(h_GammaClone);
  hStack->Add(h_C10Clone);
  hStack->Add(h_C11Clone);
  hStack->Add(h_N12Clone);
  hStack->Add(h_O15Clone);
  hStack->Draw("hist");
}

void ActProfileTotFromComponents()
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  double timeEndAct = 1; // sec (time when activation phase ends)
  double timeTot = 15*60; // sec
  int n = 2e5; // number of points
  
  Beam b;
  b.m_tp = 36e-9;
  b.m_ts = 4e-9;
  b.m_I = 8;
  
  Isotope N12(TMath::Log(2)/11e-3, 5.05425e-4, &b);
  Isotope C11(TMath::Log(2)/20.36/60, 5.0419e-3, &b);
  Isotope C10(TMath::Log(2)/19.29, 0.00879*0.0046, &b);
  Isotope O15(TMath::Log(2)/2.04/60, 2.8567e-3, &b);

  cout << "N12 tot number: " << N12.TotNumber(0, timeEndAct) << endl;
  cout << "C11 tot number: " << C11.TotNumber(0, timeEndAct) << endl;
  cout << "C10 tot number: " << C10.TotNumber(0, timeEndAct) << endl;
  cout << "O15 tot number: " << O15.TotNumber(0, timeEndAct) << endl;
  
  TGraph* gN12VsTime = new TGraph(n);
  setGraphStyle(gN12VsTime, 7);
  TGraph* gC11VsTime = new TGraph(n);
  setGraphStyle(gC11VsTime, kBlue-9);
  TGraph* gC10VsTime = new TGraph(n);
  setGraphStyle(gC10VsTime, kBlue);
  TGraph* gO15VsTime = new TGraph(n);
  setGraphStyle(gO15VsTime, 3);
  TGraph* gGammaVsTime = new TGraph(n);
  setGraphStyle(gGammaVsTime, kMagenta);
  
  std::vector<double> N12numberAt;
  std::vector<double> C11numberAt;
  std::vector<double> C10numberAt;
  std::vector<double> O15numberAt;
  std::vector<double> GammanumberAt;
  std::vector<double> times;
  
  // Init
  bool first = false;
  int i_firstOfDeactivation=0;
  for(double i=0; i<n; i++) {
    double timeLoc = i*timeTot/float(n-1);
    times.push_back(timeLoc);
    if(timeLoc > timeEndAct && !first) {
      i_firstOfDeactivation = i;
      first = true;
      cout << "First of deact:" << timeLoc << "  " << i_firstOfDeactivation << endl; 
    }
    N12numberAt.push_back(0);
    C11numberAt.push_back(0);
    C10numberAt.push_back(0);
    O15numberAt.push_back(0);
    GammanumberAt.push_back(0);
  }
  
  // Compute numbers for all times
  for(double i=0; i<times.size(); i++) {
    double t = times[i];

    if(i < i_firstOfDeactivation) {
      ////////////////////////////
      // Activation phase
      ////////////////////////////
      
      N12numberAt[i] = N12.NumberOfDecaysAfterActivation(0, t);
      C11numberAt[i] = C11.NumberOfDecaysAfterActivation(0, t);
      C10numberAt[i] = C10.NumberOfDecaysAfterActivation(0, t);
      O15numberAt[i] = O15.NumberOfDecaysAfterActivation(0, t);
      GammanumberAt[i] = 0.00879*0.0284*b.getProtonsPerSecondAverage()*t;

      //cout << t << "  sec: N12=" << N12numberAt[i] << "  C11=" << C11numberAt[i] << "  O15=" << O15numberAt[i] << endl;
    } else {
      ////////////////////////////
      // Deactivation phase
      ////////////////////////////

      N12numberAt[i] = N12numberAt[i_firstOfDeactivation-1] + N12.NumberOfDecaysAfterDeactivation(N12.TotNumber(0, timeEndAct) - N12numberAt[i_firstOfDeactivation-1], t-timeEndAct);
      C11numberAt[i] = C11numberAt[i_firstOfDeactivation-1] + C11.NumberOfDecaysAfterDeactivation(C11.TotNumber(0, timeEndAct) - C11numberAt[i_firstOfDeactivation-1], t-timeEndAct);
      C10numberAt[i] = C10numberAt[i_firstOfDeactivation-1] + C10.NumberOfDecaysAfterDeactivation(C10.TotNumber(0, timeEndAct) - C10numberAt[i_firstOfDeactivation-1], t-timeEndAct);
      O15numberAt[i] = O15numberAt[i_firstOfDeactivation-1] + O15.NumberOfDecaysAfterDeactivation(O15.TotNumber(0, timeEndAct) - O15numberAt[i_firstOfDeactivation-1], t-timeEndAct);
      GammanumberAt[i] = GammanumberAt[i_firstOfDeactivation-1];
      
      //cout << t << "  sec: N12=" << N12numberAt[i] << "  C11=" << C11numberAt[i] << "  O15=" << O15numberAt[i] << endl;
    }
    gN12VsTime->SetPoint(i, t, N12numberAt[i]);
    gC11VsTime->SetPoint(i, t, C11numberAt[i]);
    gC10VsTime->SetPoint(i, t, C10numberAt[i]);
    gO15VsTime->SetPoint(i, t, O15numberAt[i]);
    gGammaVsTime->SetPoint(i, t, GammanumberAt[i]);
  }

  TCanvas* c0 = new TCanvas("c0", "c0", 600, 600);
  TMultiGraph* multi = new TMultiGraph();
  multi->Add(gN12VsTime);
  multi->Add(gC11VsTime);
  multi->Add(gC10VsTime);
  multi->Add(gO15VsTime);
  multi->Add(gGammaVsTime);
  
  multi->Draw("ap");
  multi->GetYaxis()->SetRangeUser(0.001, 3e9);
  multi->Draw("ap");
  gPad->SetLogy();
  gPad->SetLogx();

  // Activity profile

  TFile* f = new TFile("data/Draw_Target_OnlyBeta_annihilation_option_Run_2_Setup_0.root");
  TH1* h_N12 = (TH1*) f->Get("hProfil_EachZ_1000070120_stack_6");
  TH1* h_C11 = (TH1*) f->Get("hProfil_EachZ_1000060110_stack_5");
  TH1* h_C10 = (TH1*) f->Get("hProfil_EachZ_1000060100_stack_4");
  TH1* h_O15 = (TH1*) f->Get("hProfil_EachZ_1000080150_stack_10");
  TH1* h_Gamma = (TH1*) f->Get("hProfil_EachZ_22_stack_1");

  TCanvas* c1 = new TCanvas("c1", "c1", 1800, 400);
  c1->Divide(4,1);
  c1->cd(1);
  int idx = 0.5/timeTot*(n-1);
  DrawStack(idx, h_N12, N12numberAt, h_C11, C11numberAt, h_C10, C10numberAt, h_O15, O15numberAt, h_Gamma, GammanumberAt);
  c1->cd(2);
  idx = 10/timeTot*(n-1);
  DrawStack(idx, h_N12, N12numberAt, h_C11, C11numberAt, h_C10, C10numberAt, h_O15, O15numberAt, h_Gamma, GammanumberAt);
  c1->cd(3);
  idx = 60/timeTot*(n-1);
  DrawStack(idx, h_N12, N12numberAt, h_C11, C11numberAt, h_C10, C10numberAt, h_O15, O15numberAt, h_Gamma, GammanumberAt);
  c1->cd(4);
  idx = 60*10/timeTot*(n-1);
  DrawStack(idx, h_N12, N12numberAt, h_C11, C11numberAt, h_C10, C10numberAt, h_O15, O15numberAt, h_Gamma, GammanumberAt);
  
}

