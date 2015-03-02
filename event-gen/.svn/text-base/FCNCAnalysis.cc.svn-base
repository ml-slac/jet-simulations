#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"

#include "FCNCAnalysis.h"
#include "FCNCTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// Constructor 
FCNCAnalysis::FCNCAnalysis(){
    if(fDebug) cout << "FCNCAnalysis::FCNCAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new FCNCTools();

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_largeR_ALTAS  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);

    if(fDebug) cout << "FCNCAnalysis::FCNCAnalysis End " << endl;
}

// Destructor 
FCNCAnalysis::~FCNCAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_largeR_ALTAS;
}

// Begin method
void FCNCAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for FCNC");
    
   DeclareBranches();
   ResetBranches();
   

   return;
}

// End
void FCNCAnalysis::End(){
    
    tT->Write();
    tF->Close();
    return;
}

// Analyze
void FCNCAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV){
    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    double METpx;
    double METpy;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    std::vector <fastjet::PseudoJet>           bhadrons;
    std::vector <fastjet::PseudoJet>           Bosons;
    std::vector <fastjet::PseudoJet>           chadrons;

    // Particle loop -----------------------------------------------------------
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

        fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
        p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));
	    
        if (tool->IsBHadron(pythia8->event[ip].id())){
	         bhadrons.push_back(p);
	    }
	    if (tool->IsCHadron(pythia8->event[ip].id())){
	        chadrons.push_back(p);
        }

        // truth bosons: only looking at these that decay to two particles, this excluded intermediate boson that are their own children
        //               or in final state
        int absID = abs(pythia8->event[ip].id());
        if( (  (absID == 23 || absID == 24 || absID == 25 || absID == 32 || absID ==  99999999) && pythia8->event.daughterList(ip).size() ==2 ) || 
            (  (absID == 23 || absID == 24 || absID == 25 || absID == 32 || absID ==  99999999) && pythia8->event[ip].isFinal())) { 
            
            if(fDebug){
                cout << "FCNCAnalysis::AnalyzeEvent " ;
                vector<int> daugthers = pythia8->event.daughterList(ip);
                cout << "boson " << pythia8->event[ip].id();
                if(daugthers.size() >0 ) cout << " with daughters " << pythia8->event[daugthers[0]].id() ;
                if(daugthers.size() >1 ) cout << " and " << pythia8->event[daugthers[1]].id();
                cout << endl;
            }

            Bosons.push_back(p);

            if(fTNBosonsFilled == MaxNBosons) {cout << "Warning: More than " << MaxNBosons << " bosons found!" << endl; continue;}
            fTBosonID [fTNBosonsFilled] = pythia8->event[ip].id();
            fTBosonPt [fTNBosonsFilled] = pythia8->event[ip].pT();
            fTBosonEta[fTNBosonsFilled] = pythia8->event[ip].eta();
            fTBosonPhi[fTNBosonsFilled] = pythia8->event[ip].phi();
            fTBosonM  [fTNBosonsFilled] = pythia8->event[ip].m();
            fTNBosonsFilled++;
        }

        // leptons
        bool particleIsSelLepton(false);
        // isolated electrons ---------------------------
        if (pythia8->event[ip].isFinal() && fabs(pythia8->event[ip].id())  ==11) {
            if (pythia8->event[ip].pT()  < 20)                             continue;
            if(! tool->IsIsolated(&pythia8->event[ip], pythia8, 0.2, 0.4)) continue;
            
            fTElesPt       [fTNElesFilled] = pythia8->event[ip].pT();
            fTElesEta      [fTNElesFilled] = pythia8->event[ip].eta();
            fTElesPhi      [fTNElesFilled] = pythia8->event[ip].phi();
            fTElesE        [fTNElesFilled] = pythia8->event[ip].m();
            fTElesCharge   [fTNElesFilled] = TMath::Sign(1,pythia8->event[ip].id());
            fTNElesFilled++;
            particleIsSelLepton = true;
        }

                      
        // isolated muons ---------------------------
        if (pythia8->event[ip].isFinal() && fabs(pythia8->event[ip].id())  ==13) {
            if (pythia8->event[ip].pT()  < 20)                             continue;
            if(! tool->IsIsolated(&pythia8->event[ip], pythia8, 0.2, 0.4)) continue;
            
            fTMuonsPt       [fTNMuonsFilled] = pythia8->event[ip].pT();
            fTMuonsEta      [fTNMuonsFilled] = pythia8->event[ip].eta();
            fTMuonsPhi      [fTNMuonsFilled] = pythia8->event[ip].phi();
            fTMuonsE        [fTNMuonsFilled] = pythia8->event[ip].m();
            fTMuonsCharge   [fTNMuonsFilled] = TMath::Sign(1,pythia8->event[ip].id());
            fTNMuonsFilled++;
            particleIsSelLepton = true;
        }
	
	// MET ----------------------------
	if ((fabs(pythia8->event[ip].id())  ==12) || (fabs(pythia8->event[ip].id())  ==14) || (fabs(pythia8->event[ip].id())  ==16)){
	  if (!pythia8->event[ip].isFinal() )      continue;
	  METpx+=pythia8->event[ip].px();
	  METpy+=pythia8->event[ip].py();
	}

	// particles for jets --------------
        if (!pythia8->event[ip].isFinal() )      continue;
        if (fabs(pythia8->event[ip].id())  ==11) continue;
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==13) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;
        if (pythia8->event[ip].pT()       < 0.5) continue;
        if (particleIsSelLepton)                 continue;

	    particlesForJets.push_back(p);

     } // end particle loop -----------------------------------------------

    // small R jets: ATLAS Style ------------------------------------------
    fastjet::ClusterSequence cs(particlesForJets, *m_jet_def);
    vector<fastjet::PseudoJet> myJets = fastjet::sorted_by_pt(cs.inclusive_jets(10.0)); //was 10
    for (unsigned int ij = 0; ij < myJets.size(); ij++){
        if(fTNJetsSmallRFilled == MaxNJetSmallR) {cout << "Warning: More than " << MaxNJetSmallR << " small R jets" << endl; continue;}
        fTJsmallPt       [fTNJetsSmallRFilled] = myJets[ij].pt();
        fTJsmallEta      [fTNJetsSmallRFilled] = myJets[ij].eta();
        fTJsmallPhi      [fTNJetsSmallRFilled] = myJets[ij].phi();
        fTJsmallM        [fTNJetsSmallRFilled] = myJets[ij].m();
        fTJsmallCharge   [fTNJetsSmallRFilled] = tool->JetCharge(myJets[ij],0.5);
        fTJsmallBtag     [fTNJetsSmallRFilled] = (tool->Btag(myJets[ij],bhadrons,chadrons,0.4,0.7,5,100)? 1:0);
        fTNJetsSmallRFilled++;
    }

    // large-R jets: Trimmed, ATLAS style ---------------------------------
    fastjet::Filter trimmer (fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05));
    fastjet::ClusterSequence csLargeR(particlesForJets, *m_jet_def_largeR_ALTAS);
    vector<fastjet::PseudoJet> myJetsLargeR_ATLAS = fastjet::sorted_by_pt(csLargeR.inclusive_jets(10.0)); //was 10
    for (unsigned int ij = 0; ij < myJetsLargeR_ATLAS.size(); ij++){
        fastjet::PseudoJet tj = trimmer(myJetsLargeR_ATLAS[ij]);
        if(fDebug){
        std::cout << ij << " trimmed:   charge and mass: " << tool->JetCharge(tj,                    0.5) << " " << tj.m()                     << " " << tool->Btag(tj                    ,bhadrons,chadrons,0.4,0.7,5,100) << std::endl;
        std::cout << ij << " untrimmed: charge and mass: " << tool->JetCharge(myJetsLargeR_ATLAS[ij],0.5) << " " << myJetsLargeR_ATLAS[ij].m() << " " << tool->Btag(myJetsLargeR_ATLAS[ij],bhadrons,chadrons,0.4,0.7,5,100) << std::endl;
        }
        if(fTNJetsLargeRFilled == MaxNJetLargeR) {cout << "Warning: More than " << MaxNJetLargeR << " large R jets" << endl; continue;}
        fTJlargeRPt         [fTNJetsLargeRFilled] = tj.pt();
        fTJlargeREta        [fTNJetsLargeRFilled] = tj.eta();
        fTJlargeRPhi        [fTNJetsLargeRFilled] = tj.phi();
        fTJlargeRM          [fTNJetsLargeRFilled] = tj.m();
        fTJlargeRMungroomed [fTNJetsLargeRFilled] = myJetsLargeR_ATLAS[ij].m();
        fTJlargeRCharge     [fTNJetsLargeRFilled] = tool->JetCharge(tj,0.5);
        fTJlargeRBtag       [fTNJetsLargeRFilled] = (tool->Btag(tj,bhadrons,chadrons,1,0.7,5,100)? 1:0);
        fTJlargeWplusMatch  [fTNJetsLargeRFilled] = (tool->BosonMatch(tj, Bosons, 1, 24)      ? 1:0);
        fTJlargeWminusMatch [fTNJetsLargeRFilled] = (tool->BosonMatch(tj, Bosons, 1, -24)     ? 1:0);
        fTJlargeZpMatch     [fTNJetsLargeRFilled] = (tool->BosonMatch(tj, Bosons, 1, 32)      ? 1:0);
        fTJlargeHpMatch     [fTNJetsLargeRFilled] = (tool->BosonMatch(tj, Bosons, 1, 99999999)? 1:0);
        fTNJetsLargeRFilled++;
    }



    // Event Selection: ATLAS Style ---------------------------------------                                                                                 
    
    if (fTNMuonsFilled+fTNElesFilled == 1){
      TLorentzVector mylepton;
      TVector2 MET = TVector2(METpx,METpy);
      if (fTNMuonsFilled==1){
	mylepton.SetPtEtaPhiE(fTMuonsPt[0],fTMuonsEta[0],fTMuonsPhi[0],fTMuonsE[0]);
      }
      else{
	mylepton.SetPtEtaPhiE(fTElesPt[0],fTElesEta[0],fTElesPhi[0],fTElesE[0]);
      }
      double MT = sqrt(2.0*mylepton.Pt()*MET.Mod()*(1.0 - cos( mylepton.Phi() - atan2(MET.Px(),MET.Py()))));
      int btags=0;
      for (int ijet=0;ijet<fTNJetsSmallRFilled;ijet++){
	btags+=fTJsmallBtag[ijet];
      }
      
      if (MET.Mod() > 20 && MET.Mod()+MT > 60 && btags>0){

	//Resolved Case:
	bool passresolve=false;
	double close=999.;
	int w1=-1;
	int w2=-1;
	for (int ijet=0;ijet<fTNJetsSmallRFilled; ijet++){
	  //four jets and two in the mass window
	  if (fTNJetsSmallRFilled<4) break;
	  TLorentzVector Jet1;
	  Jet1.SetPtEtaPhiM(fTJsmallPt[ijet],fTJsmallEta[ijet],fTJsmallPhi[ijet],fTJsmallM[ijet]);
	  for (int jjet=ijet+1;jjet<fTNJetsSmallRFilled; jjet++){
	    TLorentzVector Jet2;
	    Jet2.SetPtEtaPhiM(fTJsmallPt[jjet],fTJsmallEta[jjet],fTJsmallPhi[jjet],fTJsmallM[jjet]);
	    if (fabs((Jet1+Jet2).M()-80)<close){
	      close=(Jet1+Jet2).M();
	      w1=ijet;
	      w2=jjet;
	    }
	  }
	}
	if (close < 30){
	  passresolve=true;
	}

	//Boosted Case:
	int ifat=-1;
	bool passboost=false;
	for (int ijet=0;ijet<fTNJetsLargeRFilled; ijet++){
	  //at least one jet in the mass window
	  TLorentzVector fjet;
	  fjet.SetPtEtaPhiM(fTJlargeRPt[ijet],fTJlargeREta[ijet],fTJlargeRPhi[ijet],fTJlargeRM[ijet]);
	  if (fjet.Pt() > 150 && fabs(fjet.M()-80)<30){
	    ifat=ijet;
	    passboost=true;
	    break;
	  }
	}

	std::cout << "pass resolve, pass boost, mu numb, el numb " << passresolve << " " << passboost << " " << fTNMuonsFilled << " " << fTNElesFilled << std::endl;

	tT->Fill();
      }
    }

    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void FCNCAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   
   // smallR jets
   tT->Branch("NJetsFilledSmallR",         &fTNJetsSmallRFilled,       "NJetsFilledSmallR/I");
   tT->Branch("JsmallPt",                  &fTJsmallPt,                "JsmallPt[NJetsFilledSmallR]/F");
   tT->Branch("JsmallEta",                 &fTJsmallEta,               "JsmallEta[NJetsFilledSmallR]/F");
   tT->Branch("JsmallPhi",                 &fTJsmallPhi,               "JsmallPhi[NJetsFilledSmallR]/F");
   tT->Branch("JsmallM",                   &fTJsmallM,                 "JsmallM[NJetsFilledSmallR]/F");
   tT->Branch("JsmallCharge",              &fTJsmallCharge,            "JsmallCharge[NJetsFilledSmallR]/F");
   tT->Branch("JsmallBtag",                &fTJsmallBtag,              "JsmallBtag[NJetsFilledSmallR]/I");

   // largeR jets
   tT->Branch("NJetsFilledLargeR",         &fTNJetsLargeRFilled,        "NJetsFilledLargeR/I");
   tT->Branch("JlargeRPt",                 &fTJlargeRPt,                "JlargeRPt[NJetsFilledLargeR]/F");
   tT->Branch("JlargeREta",                &fTJlargeREta,               "JlargeREta[NJetsFilledLargeR]/F");
   tT->Branch("JlargeRPhi",                &fTJlargeRPhi,               "JlargeRPhi[NJetsFilledLargeR]/F");
   tT->Branch("JlargeRM",                  &fTJlargeRM,                 "JlargeRM[NJetsFilledLargeR]/F");
   tT->Branch("JlargeRMungroomed",         &fTJlargeRMungroomed,        "JlargeRMungroomed[NJetsFilledLargeR]/F");
   tT->Branch("JlargeRCharge",             &fTJlargeRCharge,            "JlargeRCharge[NJetsFilledLargeR]/F");
   tT->Branch("JlargeRBtag",               &fTJlargeRBtag,              "JlargeRBtag[NJetsFilledLargeR]/I");
   tT->Branch("JlargeWplusMatch",          &fTJlargeWplusMatch,         "JlargeWplusMatch[NJetsFilledLargeR]/I");
   tT->Branch("JlargeWminusMatch",         &fTJlargeWminusMatch,        "JlargeWminusMatch[NJetsFilledLargeR]/I");
   tT->Branch("JlargeZpMatch",             &fTJlargeZpMatch,            "JlargeZpMatch[NJetsFilledLargeR]/I");
   tT->Branch("JlargeHpMatch",             &fTJlargeHpMatch,            "JlargeHpMatch[NJetsFilledLargeR]/I");

   // Bosons
   tT->Branch("NBosonsFilled",             &fTNBosonsFilled,          "NBosonsFilled/I");
   tT->Branch("BosonID",                   &fTBosonID,                "BosonID[NBosonsFilled]/I");
   tT->Branch("BosonPt",                   &fTBosonPt,                "BosonPt[NBosonsFilled]/F");
   tT->Branch("BosonEta",                  &fTBosonEta,               "BosonEta[NBosonsFilled]/F");
   tT->Branch("BosonPhi",                  &fTBosonPhi,               "BosonPhi[NBosonsFilled]/F");
   tT->Branch("BosonM",                    &fTBosonM,                 "BosonM[NBosonsFilled]/F");

   // Eles
   tT->Branch("NElesFilled",              &fTNElesFilled,           "NElesFilled/I");
   tT->Branch("ElesPt",                   &fTElesPt,                "ElesPt[NElesFilled]/F");
   tT->Branch("ElesEta",                  &fTElesEta,               "ElesEta[NElesFilled]/F");
   tT->Branch("ElesPhi",                  &fTElesPhi,               "ElesPhi[NElesFilled]/F");
   tT->Branch("ElesE",                    &fTElesE,                 "ElesE[NElesFilled]/F");
   tT->Branch("ElesCharge",               &fTElesCharge,            "ElesCharge[NElesFilled]/I");

   // Muons
   tT->Branch("NMuonsFilled",              &fTNMuonsFilled,           "NMuonsFilled/I");
   tT->Branch("MuonsPt",                   &fTMuonsPt,                "MuonsPt[NMuonsFilled]/F");
   tT->Branch("MuonsEta",                  &fTMuonsEta,               "MuonsEta[NMuonsFilled]/F");
   tT->Branch("MuonsPhi",                  &fTMuonsPhi,               "MuonsPhi[NMuonsFilled]/F");
   tT->Branch("MuonsE",                    &fTMuonsE,                 "MuonsE[NMuonsFilled]/F");
   tT->Branch("MuonsCharge",               &fTMuonsCharge,            "MuonsCharge[NMuonsFilled]/I");

   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void FCNCAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNBosonsFilled               = 0;
      fTNJetsSmallRFilled           = 0;
      fTNJetsLargeRFilled           = 0;
      fTNElesFilled                 = 0;
      fTNMuonsFilled                = 0;
      for (int iP=0; iP < MaxNEles; ++iP){
          fTElesPt     [iP]= -999;
          fTElesEta    [iP]= -999;
          fTElesPhi    [iP]= -999;
          fTElesE      [iP]= -999;
          fTElesCharge [iP]= -999;
      }
      for (int iP=0; iP < MaxNMuons; ++iP){
          fTMuonsPt     [iP]= -999;
          fTMuonsEta    [iP]= -999;
          fTMuonsPhi    [iP]= -999;
          fTMuonsE      [iP]= -999;
          fTMuonsCharge [iP]= -999;
      }
      for (int iP=0; iP < MaxNBosons; ++iP){
          fTBosonPt  [iP]= -999;
          fTBosonEta [iP]= -999;
          fTBosonPhi [iP]= -999;
          fTBosonM   [iP]= -999;
          fTBosonID  [iP]= -999;
      }
      for (int iP=0; iP < MaxNJetSmallR; ++iP){
          fTJsmallPt      [iP]= -999;
          fTJsmallPhi     [iP]= -999;
          fTJsmallEta     [iP]= -999;
          fTJsmallM       [iP]= -999;
          fTJsmallCharge  [iP]= -999;
          fTJsmallBtag    [iP]= -999;
      }
      for (int iP=0; iP < MaxNJetLargeR; ++iP){
          fTJlargeRPt          [iP]= -999;
          fTJlargeRPhi         [iP]= -999;
          fTJlargeREta         [iP]= -999;
          fTJlargeRM           [iP]= -999;
          fTJlargeRMungroomed  [iP]= -999;
          fTJlargeRCharge      [iP]= -999;
          fTJlargeRBtag        [iP]= -999;
          fTJlargeWplusMatch   [iP]= -999;
          fTJlargeWminusMatch  [iP]= -999;
          fTJlargeZpMatch      [iP]= -999;
          fTJlargeHpMatch      [iP]= -999;
      }
}
