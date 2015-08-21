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
#include "TH2F.h"

#include "MIAnalysis.h"
#include "MITools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"

#include "Pythia8/Pythia.h"

// #include "Nsubjettiness.h"
#include "Njettiness.hh"
#include "Nsubjettiness.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

typedef pair<double, double> point;

double euclidean_distance(const point &x, const point &y)
{
    double d1 = x.first - y.first;
    double d2 = x.second - y.second;

    return sqrt(d1 * d1 + d2 * d2);
}

// Constructor 
MIAnalysis::MIAnalysis(int imagesize)
{
    imagesize *= imagesize;
    MaxN = imagesize;
    fTIntensity = new float[imagesize];
    // fTRotatedIntensity = new float[imagesize];
    // fTLocalDensity = new float[imagesize];
    // fTGlobalDensity = new float[imagesize];


    if(fDebug) cout << "MIAnalysis::MIAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new MITools();

    //model the detector as a 2D histogram   
    //                         xbins       y bins
    detector = new TH2F("", "", 100, -5, 5, 200, -10, 10);
    for(int i = 1; i <= 100; i++)
    {
        for (int j = 1; j <= 200; j++)
        {
            detector->SetBinContent(i,j,0);
        }
    }

    if(fDebug) cout << "MIAnalysis::MIAnalysis End " << endl;
}

// Destructor 
MIAnalysis::~MIAnalysis()
{
    delete tool;

    delete[] fTIntensity;
    // delete[] fTRotatedIntensity;
    // delete[] fTGlobalDensity;
    // delete[] fTLocalDensity;
}

// Begin method
void MIAnalysis::Begin()
{
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for MI");
   
   // for shit you want to do by hand
   DeclareBranches();
   ResetBranches();
   
   return;
}

// End
void MIAnalysis::End()
{
    tT->Write();
    tF->Close();
    return;
}

// Analyze
void MIAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV,
    int pixels, float range)
{

    if(fDebug) cout << "MIAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) cout << "MIAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet> particlesForJets;

    detector->Reset();
   
    // Particle loop ----------------------------------------------------------
    for (int ip=0; ip<pythia8->event.size(); ++ip){

        fastjet::PseudoJet p(pythia8->event[ip].px(),
                             pythia8->event[ip].py(), 
                             pythia8->event[ip].pz(),
                             pythia8->event[ip].e() );

        // particles for jets --------------
        if (!pythia8->event[ip].isFinal())       continue;

        //Skip neutrinos, PDGid = 12, 14, 16 
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;


        // find the particles rapidity and phi, then get the detector bins
        int ybin = detector->GetXaxis()->FindBin(p.rapidity());
        int phibin = detector->GetYaxis()->FindBin(p.phi());

        // do bin += value in the associated detector bin
        detector->SetBinContent(ybin, phibin, 
                                detector->GetBinContent(ybin, phibin) + p.e());
      
    }  
    // end particle loop -----------------------------------------------  

    //Now, we extract the energy from the calorimeter for processing by fastjet
    for (int i = 1; i <= detector->GetNbinsX(); i++)
    {
        for (int j = 1; j <= detector->GetNbinsY(); j++)
        {
            if (detector->GetBinContent(i, j) > 0)
            {
                double phi = detector->GetYaxis()->GetBinCenter(j);
                double eta = detector->GetXaxis()->GetBinCenter(i);
                double E = detector->GetBinContent(i, j);
                fastjet::PseudoJet p(0., 0., 0., 0.);

                //We measure E (not pT)!  And treat 'clusters' as massless.
                p.reset_PtYPhiM(E/cosh(eta), eta, phi, 0.); 
                particlesForJets.push_back(p);
            }
        }
    }

    fastjet::JetDefinition *m_jet_def = new fastjet::JetDefinition(
        fastjet::antikt_algorithm, 1.0);

    fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
        fastjet::SelectorPtFractionMin(1e-6));

    fastjet::ClusterSequence csLargeR(particlesForJets, *m_jet_def);

    vector<fastjet::PseudoJet> considered_jets = fastjet::sorted_by_pt(
        csLargeR.inclusive_jets(10.0));
    fastjet::PseudoJet leading_jet = trimmer(considered_jets[0]);

    //Now, let's make an image out of the leading jet. 
    vector<fastjet::PseudoJet> subjets = sorted_by_pt(leading_jet.pieces());

    if (subjets.size() < 2)
    {
        return;
    }

    // bool overall = true;
    // Dump kinematics of the leading jet.
    // For the leading SUBjet, put subjets[0].
    fTLeadingEta = leading_jet.eta();
    fTLeadingPhi = leading_jet.phi();
    fTLeadingPt = leading_jet.perp();
    fTLeadingM = leading_jet.m();
    
    
    vector<pair<double, double>  > consts_image;
    vector<pair<double, double>  > subjets_image; 

    for (int i = 0; i < 2; i++)
    {
        pair<double, double> subjet_hold;
        subjet_hold.first = subjets[i].eta();
        subjet_hold.second = subjets[i].phi();
        subjets_image.push_back(subjet_hold);
    }
    
    vector<fastjet::PseudoJet> sorted_consts = sorted_by_pt(leading_jet.constituents());

    for(int i = 0; i < sorted_consts.size(); i++)
    {
        pair<double, double> const_hold;
        const_hold.first = sorted_consts[i].eta();
        const_hold.second = sorted_consts[i].phi();
        consts_image.push_back(const_hold);
    }

    //Step 1: Center on the leading subjet
    for (int i =0; i < sorted_consts.size(); i++)
    {
        consts_image[i].first = consts_image[i].first-subjets_image[0].first;
        consts_image[i].second = consts_image[i].second-subjets_image[0].second;
    }
    for (int i =1; i >= 0; i--)
    {
        subjets_image[i].first = subjets_image[i].first - subjets_image[0].first;
        subjets_image[i].second = subjets_image[i].second- subjets_image[0].second;
    }

    //Step 2: Fill in the unrotated image
    //-------------------------------------------------------------------------   
    TH2F* orig_im = new TH2F("", "", pixels, -range, range, pixels, -range, range);

    for (int i = 0; i < sorted_consts.size(); i++)
    {
        orig_im->Fill(consts_image[i].first,consts_image[i].second,sorted_consts[i].e());
      //std::cout << i << "       " << consts_image[i].first  << " " << consts_image[i].second << std::endl;  
    }

    //Step 2b): fill in the density
    //-------------------------------------------------------------------------
    // double Rlocal = 0.5;
    // double Rglobal = 1.0;
    // TH2F* localdensity = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // TH2F* globaldensity = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // for (int i = 0; i < sorted_consts.size(); ++i)
    // {
    //     auto reference_px = consts_image[i];
    //     auto ref_pt = sorted_consts[i].e();

    //     double global_acc = 0;
    //     double local_acc = 0;

    //     for (int j = 0; j < sorted_consts.size(); ++j)
    //     {   
    //         if (j != i)
    //         {
    //             auto dR = euclidean_distance(reference_px, consts_image[j]);
    //             if (dR <= Rglobal)
    //             {
    //                 global_acc += (ref_pt / dR);
    //                 if (dR <= Rlocal)
    //                 {
    //                     local_acc += (ref_pt / dR);
    //                 }
    //             }
    //         }
    //     }
    //     localdensity->Fill(reference_px.first, reference_px.second, local_acc);
    //     globaldensity->Fill(reference_px.first, reference_px.second, global_acc);
    // }


    fTdeltaR = euclidean_distance(subjets_image[0], subjets_image[1]);

    //Step 3: Rotate so the subleading subjet is at -pi/2
    //-------------------------------------------------------------------------
    fTSubLeadingEta = subjets_image[1].first;
    fTSubLeadingPhi = subjets_image[1].second;
    double theta = atan(subjets_image[1].second/subjets_image[1].first)+2.*atan(1.); //atan(1) = pi/4
    
    if (-sin(theta)*subjets_image[1].first+cos(theta)*subjets_image[1].second > 0)
    {
        theta+=-4.*atan(1.);
    }

    fTRotationAngle = theta;

    for (int i=0; i<sorted_consts.size(); i++)
    {
        double x = consts_image[i].first;
        double y = consts_image[i].second;
        consts_image[i].first = cos(theta)*x + sin(theta)*y;
        consts_image[i].second = -sin(theta)*x + cos(theta)*y;
    }

    for (int i=0; i<2; i++)
    {
        double x = subjets_image[i].first;
        double y = subjets_image[i].second;
        subjets_image[i].first = cos(theta)*x + sin(theta)*y;
        subjets_image[i].second =  -sin(theta)*x + cos(theta)*y;
    }

    theta = atan(subjets_image[1].second/subjets_image[1].first);

    // //Step 4b): fill in the rotated image
    // //-------------------------------------------------------------------------
    // TH2F* rotatedimage = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // for (int i = 0; i < sorted_consts.size(); i++)
    // {
    //     rotatedimage->Fill(consts_image[i].first,consts_image[i].second,sorted_consts[i].e());
    // }
    

    //Step 5: Dump the images in the tree!
    //-------------------------------------------------------------------------
    int counter=0;
    for (int i=1; i<=orig_im->GetNbinsX(); i++)
    {
        for (int j=1; j<=orig_im->GetNbinsY(); j++)
        {
            // fTRotatedIntensity[counter] = rotatedimage->GetBinContent(i,j);
            fTIntensity[counter] = orig_im->GetBinContent(i,j);
            // fTLocalDensity[counter] = localdensity->GetBinContent(i, j);
            // fTGlobalDensity[counter] = globaldensity->GetBinContent(i, j);

            counter++;
        }
    }

    // Step 6: Fill in nsubjettiness (new)
    //----------------------------------------------------------------------------
    OnePass_WTA_KT_Axes axis_spec;
    NormalizedMeasure parameters(1.0, 1.0);

    // NormalizedMeasure parameters(1.0, 1.0);
    Nsubjettiness subjettiness_1(1, axis_spec, parameters);
    Nsubjettiness subjettiness_2(2, axis_spec, parameters);
    Nsubjettiness subjettiness_3(3, axis_spec, parameters);

    fTTau1 = (float) subjettiness_1.result(leading_jet);
    fTTau2 = (float) subjettiness_2.result(leading_jet);
    fTTau3 = (float) subjettiness_3.result(leading_jet);

    fTTau32 = (abs(fTTau2) < 1e-4 ? -10 : fTTau3 / fTTau2);
    fTTau21 = (abs(fTTau1) < 1e-4 ? -10 : fTTau2 / fTTau1);


    // // Step 7: Fill in nsubjettiness (old)
    // //----------------------------------------------------------------------------
    // OnePass_KT_Axes axis_spec_old;
    // NormalizedMeasure parameters_old(1.0, 1.0);

    // Nsubjettiness old_subjettiness_1(1, axis_spec_old, parameters_old);
    // Nsubjettiness old_subjettiness_2(2, axis_spec_old, parameters_old);
    // Nsubjettiness old_subjettiness_3(3, axis_spec_old, parameters_old);

    // fTTau1 = (float) old_subjettiness_1.result(leading_jet);
    // fTTau2 = (float) old_subjettiness_2.result(leading_jet);
    // fTTau3 = (float) old_subjettiness_3.result(leading_jet);

    // fTTau32old = (abs(fTTau2) < 1e-4 ? -10 : fTTau3 / fTTau2);
    // fTTau21old = (abs(fTTau1) < 1e-4 ? -10 : fTTau2 / fTTau1);

    tT->Fill();

    return;
}

// declate branches
void MIAnalysis::DeclareBranches()
{

    // Event Properties 
    tT->Branch("NFilled", &fTNFilled, "NFilled/I");

    tT->Branch("Intensity", *&fTIntensity, "Intensity[NFilled]/F");

    // tT->Branch("LocalDensity", *&fTLocalDensity, "LocalDensity[NFilled]/F");
    // tT->Branch("GlobalDensity", *&fTGlobalDensity, "GlobalDensity[NFilled]/F");

    // tT->Branch("RotatedIntensity", 
    //     *&fTRotatedIntensity, "RotatedIntensity[NFilled]/F");

    tT->Branch("SubLeadingEta", &fTSubLeadingEta, "SubLeadingEta/F");
    tT->Branch("SubLeadingPhi", &fTSubLeadingPhi, "SubLeadingPhi/F");

    tT->Branch("LeadingEta", &fTLeadingEta, "LeadingEta/F");
    tT->Branch("LeadingPhi", &fTLeadingPhi, "LeadingPhi/F");
    tT->Branch("LeadingPt", &fTLeadingPt, "LeadingPt/F");
    tT->Branch("LeadingM", &fTLeadingM, "LeadingM/F");
    tT->Branch("RotationAngle", &fTRotationAngle, "RotationAngle/F");

    tT->Branch("Tau1", &fTTau1, "Tau1/F");
    tT->Branch("Tau2", &fTTau2, "Tau2/F");
    tT->Branch("Tau3", &fTTau3, "Tau3/F");

    tT->Branch("DeltaR", &fTdeltaR, "DeltaR/F");

    tT->Branch("Tau32", &fTTau32, "Tau32/F");
    tT->Branch("Tau21", &fTTau21, "Tau21/F");
    // tT->Branch("Tau32old", &fTTau32old, "Tau32old/F");
    // tT->Branch("Tau21old", &fTTau21old, "Tau21old/F");
    return;
}

// resets vars
void MIAnalysis::ResetBranches(){
    // reset branches 
    fTNFilled = MaxN;
    fTSubLeadingPhi = -999;
    fTSubLeadingEta = -999;
    fTRotationAngle = -999;
  
    fTTau32 = -999;
    fTTau21 = -999;

    fTTau1 = 999;
    fTTau2 = 999;
    fTTau3 = 999;

    // fTTau32old = -999;
    // fTTau21old = -999;

    fTLeadingEta = -999;
    fTLeadingPhi = -999;
    fTLeadingPt = -999;
    fTLeadingM = -999;

    for (int iP=0; iP < MaxN; ++iP)
    {
        fTIntensity[iP]= -999;
        // fTRotatedIntensity[iP]= -999;
        // fTLocalDensity[iP]= -999;
        // fTGlobalDensity[iP]= -999;
    }
}
