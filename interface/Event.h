#ifndef EVENT_H_
#define EVENT_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TParameter.h"
#include "TEnv.h"
#include "TF1.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Analysis/goldilocks/interface/physicsObjects.h"
#include "Analysis/goldilocks/interface/configuration.h"
#include "Analysis/goldilocks/interface/truthMatching.h"
#include "Analysis/goldilocks/interface/deepLearning.h"
#include "Analysis/goldilocks/interface/ttbarReco.h"

#ifdef __CINT__
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<float> >+;
#endif



// Event Class
class Event {
  public:
    // Constructor
    Event( TTreeReader &myReader, configuration &cmaConfig);
    Event( const Event &obj);

    // Destructor
    virtual ~Event();

    // create hash tables in truth/reco TTree to match truth <-> reco events
    // uses configuration option matchTruthToReco to match truth to reco (reco loop)
    // OR match reco to truth (truth loop, for acceptance studies)
    void matchTruthWithReco();
    // check during looping over truth events, if reco event match is found
    bool isValidRecoEntry(){ return (m_entry > (long long)-1);}

    // Execute the event (load information and setup objects)
    virtual void execute(Long64_t entry);
    virtual void updateEntry(Long64_t entry);

    // Setup physics information
    void initialize_leptons();
    void initialize_neutrinos();
    void initialize_jets();
    void initialize_ljets();
    void initialize_eventWeights();
    void initialize_weights();
    void initialize_kinematics();
    void initialize_truth();

    virtual double getSystEventWeight(const std::string &syst, const int weightIndex=-1);

    // Clear stuff;
    virtual void finalize();
    virtual void clear();

    // Get physics object(s) information
    std::vector<Lepton> leptons() {return m_leptons;}
    std::vector<Neutrino> neutrinos() {return m_neutrinos;}
    std::vector<Ljet> ljets() {return m_ljets;}
    std::vector<Ljet> ljetsPUPPI() {return m_ljetsPUPPI;}
    std::vector<Jet> jets() {return m_jets;}

    virtual std::vector<int> btag_jets(const std::string &wkpt);
    virtual std::vector<int> btag_jets() {return m_btag_jets_default;}  // using configured b-tag WP

    virtual float HT() {return m_HT;}
    virtual float ST() {return m_ST;}
    virtual float met( const std::string& met_name );

    // Get truth physics information 
    std::vector<Lepton> truth_leptons() {return m_truth_leptons;}
    std::vector<Neutrino> truth_neutrinos() {return m_truth_neutrinos;}
    std::vector<Ljet> truth_ljets() {return m_truth_ljets;}
    std::vector<Jet>  truth_jets() {return m_truth_jets;}
    TruthTop m_truth_top;
    TruthTop m_truth_antitop;

    // metadata
//    virtual unsigned long long eventNumber();
//    virtual unsigned int runNumber();
    long long entry() {return m_entry;}
    virtual float eventNumber();
    virtual float runNumber();
    virtual int lumiblock();
    virtual float xsection() {return m_xsection;}
    virtual float kfactor()  {return m_kfactor;}
    virtual float sumOfWeights() {return m_sumOfWeights;}
    virtual std::string treeName() {return m_treeName;}

    // functions to calculate things
    void getBtaggedJets( Jet& jet );
    void getDNNInputs();      // return the DNN inputs to the user
    void getDNN();            // get the DNN output
    float DNN(){ return m_DNN;}

    void buildTtbar();
    void overlapRemoval(std::vector<int>& new_objects);
    std::vector<Top> ttbar() {return m_ttbar;}
    void deepLearningPrediction(Top& top);

    // Get weights
    virtual float nominal_weight() {return m_nominal_weight;}
    float weight_mc() {return 1.0;}
    float truth_weight_mc() {return 1.0;}
    float weight_jvt() {return 1.0;}
    float weight_pileup() {return 1.0;}
    float weight_lept_eff() {return 1.0;}
    float weight_btag() {return m_weight_btag_default;}
    float weight_btag(const std::string &wkpt);

    // Get weight systematics
    virtual std::map<std::string,float > weightSystematicsFloats();
    virtual std::map<std::string,std::vector<float> > weightSystematicsVectorFloats();
    virtual std::vector<std::string> listOfWeightSystematics() {return m_listOfWeightSystematics;}

  protected:

    // general information
    configuration *m_config;
    TTreeReader &m_ttree;
    TTreeReader m_truth_tree;
    std::string m_treeName;
    bool m_grid;
    bool m_isMC;
    long long m_entry;
    long long m_truth_entry;
    bool m_getDNN;
    bool m_buildNeutrinos;

    bool m_ee;
    bool m_mumu;
    bool m_emu;

    // event weight information
    double m_nominal_weight;
    double m_xsection;
    double m_kfactor;
    double m_sumOfWeights;
    double m_LUMI;
    std::map<int, float> m_mapXSection; // map DSID to XSection
    std::map<int, float> m_mapKFactor;  // map DSID to KFactor
    std::map<int, float> m_mapAMI;      // map DSID to sum of weights

    // physics object information
    std::vector<Lepton> m_leptons;
    std::vector<Neutrino> m_neutrinos;
    std::vector<Ljet> m_ljets;
    std::vector<Ljet> m_ljetsPUPPI;
    std::vector<Jet>  m_jets;

    // truth physics object information
    std::vector<Lepton> m_truth_leptons;
    std::vector<Neutrino> m_truth_neutrinos;
    std::vector<Ljet> m_truth_ljets;
    std::vector<Jet>  m_truth_jets;
    std::vector<Parton> m_truth_partons;
    std::vector<TruthTop> m_truth_tops;
    float m_truthDR = 0.6;  // truth-matching radius

    // b-tagged calo jets with various WP
    std::map<std::string, std::vector<int> > m_btag_jets;
    std::vector<int> m_btag_jets_default;
    float m_cMVAv2L;
    float m_cMVAv2M;
    float m_cMVAv2T;

    float m_HT;
    float m_ST;
    float m_met_met;
    float m_met_phi;

    lwt::LightweightNeuralNetwork* m_lwnn;
    std::map<std::string, double> m_dnnInputs;   // values for inputs to the DNN
    std::string m_dnnKey;
    float m_DNN;   // DNN score

    bool m_useLargeRJets;
    bool m_useJets;
    bool m_useLeptons;
    bool m_useNeutrinos;
    bool m_useTruth;

    ttbarReco* m_ttbarRecoTool;            // tool to perform ttbar reconstruction
    deepLearning* m_deepLearningTool;      // tool to perform deep learning
    truthMatching* m_truthMatchingTool;    // tool to perform truth-matching
    std::vector<Top> m_ttbar;              // container for ttbar system
    std::vector<int> m_ak4candidates;      // AK4 jets that aren't inside AK8 candidates and have matches
    std::vector<int> m_ak8candidates;      // AK8 jets that have matches
    std::vector<int> m_ak8candidatesPUPPI; // AK8 jets that have matches
    std::map<std::string,int> m_mapContainment;
    std::map<std::string,int> m_targetMap;

    // nominal b-tagging weight maps
    std::map<std::string, float> m_weight_trackjet_btag;
    std::map<std::string, float> m_weight_btag;
    float m_weight_btag_default;
    // Maps to keep track of weight systematics
    std::map<std::string,TTreeReaderValue<float> * > m_weightSystematicsFloats;
    std::map<std::string,TTreeReaderValue<std::vector<float>> * > m_weightSystematicsVectorFloats;
    std::vector<std::string> m_listOfWeightSystematics;

    // TTree variables [all possible ones]

    // *************
    // the following are from root files accessed 
    //    on 9 January 2018
    //    at root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Summer16_80X_Jan_2017_Ntp_v12X/
    TTreeReaderValue<double> * m_dnn_score;

    TTreeReaderValue<double> * m_avg_npv;
    TTreeReaderValue<double> * m_calomet;
    TTreeReaderValue<double> * m_calometphi;
    TTreeReaderValue<double> * m_dPhi0_CUT;
    TTreeReaderValue<double> * m_dPhi1_CUT;
    TTreeReaderValue<double> * m_dPhi2_CUT;
    TTreeReaderValue<double> * m_evtWeight;
    TTreeReaderValue<double> * m_genHT;
    TTreeReaderValue<double> * m_genmet;
    TTreeReaderValue<double> * m_genmetphi;
    TTreeReaderValue<double> * m_ht;
    TTreeReaderValue<double> * m_met;
    TTreeReaderValue<double> * m_metphi;
    TTreeReaderValue<double> * m_mht;
    TTreeReaderValue<double> * m_mhtphi;
    TTreeReaderValue<double> * m_mt2;
    TTreeReaderValue<double> * m_q;
    TTreeReaderValue<double> * m_stored_weight;
    TTreeReaderValue<double> * m_tru_npv;
    TTreeReaderValue<double> * m_x1;
    TTreeReaderValue<double> * m_x2;

    TTreeReaderValue<int> * m_globalTightHalo2016Filter;
    TTreeReaderValue<int> * m_goodVerticesFilter;
    TTreeReaderValue<int> * m_eeBadScFilter;
    TTreeReaderValue<int> * m_EcalDeadCellTriggerPrimitiveFilter;
    TTreeReaderValue<int> * m_nMuons_CUT;
    TTreeReaderValue<int> * m_nMuons;
    TTreeReaderValue<int> * m_nElectrons_CUT;
    TTreeReaderValue<int> * m_nElectrons;
    TTreeReaderValue<int> * m_nJets;
    TTreeReaderValue<int> * m_NJetsISR;
    TTreeReaderValue<int> * m_id1;
    TTreeReaderValue<int> * m_id2;
    TTreeReaderValue<int> * m_loose_nIsoTrks;
    TTreeReaderValue<int> * m_nIsoTrks_CUT;
    TTreeReaderValue<int> * m_nJets_CUT;
    TTreeReaderValue<int> * m_vtxSize;
    TTreeReaderValue<int> * m_npv;
    TTreeReaderValue<int> * m_nm1;
    TTreeReaderValue<int> * m_n0;
    TTreeReaderValue<int> * m_np1;
    TTreeReaderValue<unsigned int> * m_run;
    TTreeReaderValue<unsigned int> * m_lumi;
    TTreeReaderValue<unsigned long long> * m_event;
    TTreeReaderValue<unsigned int> * m_looseJetID;
    TTreeReaderValue<unsigned int> * m_tightJetID;
    TTreeReaderValue<unsigned int> * m_tightlepvetoJetID;
    TTreeReaderValue<unsigned int> * m_looseJetID_NoLep;
    TTreeReaderValue<unsigned int> * m_tightJetID_NoLep;
    TTreeReaderValue<unsigned int> * m_tightlepvetoJetID_NoLep;
    TTreeReaderValue<unsigned int> * m_BadChargedCandidateFilter;
    TTreeReaderValue<unsigned int> * m_BadPFMuonFilter;
    TTreeReaderValue<unsigned int> * m_HBHENoiseFilter;
    TTreeReaderValue<unsigned int> * m_HBHEIsoNoiseFilter;
    TTreeReaderValue<std::vector<double>> * m_muonsCharge;
    TTreeReaderValue<std::vector<double>> * m_muonsMtw;
    TTreeReaderValue<std::vector<double>> * m_muonsRelIso;
    TTreeReaderValue<std::vector<double>> * m_muonsMiniIso;
    TTreeReaderValue<std::vector<double>> * m_muonspfActivity;
    TTreeReaderValue<std::vector<double>> * m_specialFixMuonsCharge;
    TTreeReaderValue<std::vector<double>> * m_elesCharge;
    TTreeReaderValue<std::vector<double>> * m_elesMtw;
    TTreeReaderValue<std::vector<double>> * m_elesRelIso;
    TTreeReaderValue<std::vector<double>> * m_elesMiniIso;
    TTreeReaderValue<std::vector<double>> * m_elespfActivity;
    TTreeReaderValue<std::vector<double>> * m_recoJetsJecUnc;
    TTreeReaderValue<std::vector<double>> * m_recoJetsJecScaleRawToFull;
    TTreeReaderValue<std::vector<double>> * m_qgLikelihood;
    TTreeReaderValue<std::vector<double>> * m_qgPtD;
    TTreeReaderValue<std::vector<double>> * m_qgAxis2;
    TTreeReaderValue<std::vector<double>> * m_recoJetschargedHadronEnergyFraction;
    TTreeReaderValue<std::vector<double>> * m_recoJetschargedEmEnergyFraction;
    TTreeReaderValue<std::vector<double>> * m_recoJetsneutralEmEnergyFraction;
    TTreeReaderValue<std::vector<double>> * m_recoJetsmuonEnergyFraction;
    TTreeReaderValue<std::vector<double>> * m_recoJetsBtag_0;
    TTreeReaderValue<std::vector<double>> * m_recoJetsCharge_0;

    TTreeReaderValue<std::vector<double>> * m_deepFlavor_b;
    TTreeReaderValue<std::vector<double>> * m_deepFlavor_bb;
    TTreeReaderValue<std::vector<double>> * m_deepFlavor_lepb;
    TTreeReaderValue<std::vector<double>> * m_deepFlavor_c;
    TTreeReaderValue<std::vector<double>> * m_deepFlavor_uds;
    TTreeReaderValue<std::vector<double>> * m_deepFlavor_g;

    TTreeReaderValue<std::vector<double>> * m_tau1;
    TTreeReaderValue<std::vector<double>> * m_tau2;
    TTreeReaderValue<std::vector<double>> * m_tau3;
    TTreeReaderValue<std::vector<double>> * m_softDropMass;
    TTreeReaderValue<std::vector<double>> * m_ak8SubJetsBdisc;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8JetsDeepAK8;
    TTreeReaderValue<std::vector<double>> * m_puppitau1;
    TTreeReaderValue<std::vector<double>> * m_puppitau2;
    TTreeReaderValue<std::vector<double>> * m_puppitau3;
    TTreeReaderValue<std::vector<double>> * m_puppisoftDropMass;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_puppiSubJetsBdisc;
    TTreeReaderValue<std::vector<double>> * m_recoJetsJecUncLepCleaned;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_qgLikelihood;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_qgPtD;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_qgAxis2;
    TTreeReaderValue<std::vector<double>> * m_recoJetschargedHadronEnergyFractionLepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetsneutralEmEnergyFractionLepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetschargedEmEnergyFractionLepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetsmuonEnergyFractionLepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetsBtag_0_LepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetsCharge_0_LepCleaned;
    TTreeReaderValue<std::vector<double>> * m_recoJetsJecScaleRawToFull_LepCleaned;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_tau1;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_tau2;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_tau3;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_puppisoftDropMass;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_puppitau1;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_puppitau2;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_puppitau3;
    TTreeReaderValue<std::vector<double>> * m_prodJetsNoLep_puppiSubJetsBdisc;
    TTreeReaderValue<std::vector<double>> * m_W_emu_pfActivityVec;
    TTreeReaderValue<std::vector<double>> * m_W_tau_emu_pfActivityVec;
    TTreeReaderValue<std::vector<double>> * m_W_tau_prongs_pfActivityVec;
    TTreeReaderValue<std::vector<double>> * m_ScaleWeightsMiniAOD;
    TTreeReaderValue<std::vector<double>> * m_trksForIsoVeto_charge;
    TTreeReaderValue<std::vector<double>> * m_trksForIsoVeto_dz;
    TTreeReaderValue<std::vector<double>> * m_trksForIsoVeto_iso;
    TTreeReaderValue<std::vector<double>> * m_trksForIsoVeto_pfActivity;
    TTreeReaderValue<std::vector<double>> * m_loose_isoTrks_charge;
    TTreeReaderValue<std::vector<double>> * m_loose_isoTrks_dz;
    TTreeReaderValue<std::vector<double>> * m_loose_isoTrks_iso;
    TTreeReaderValue<std::vector<double>> * m_loose_isoTrks_mtw;
    TTreeReaderValue<std::vector<double>> * m_loose_isoTrks_pfActivity;
    TTreeReaderValue<std::vector<double>> * m_metMagUp;
    TTreeReaderValue<std::vector<double>> * m_metMagDown;
    TTreeReaderValue<std::vector<double>> * m_metPhiUp;
    TTreeReaderValue<std::vector<double>> * m_metPhiDown;
    TTreeReaderValue<std::vector<int>> * m_PassTrigger;
    TTreeReaderValue<std::vector<int>> * m_TriggerPrescales;
    TTreeReaderValue<std::vector<int>> * m_muonsFlagMedium;
    TTreeReaderValue<std::vector<int>> * m_muonsFlagTight;
    TTreeReaderValue<std::vector<int>> * m_specialFixtype;
    TTreeReaderValue<std::vector<int>> * m_elesFlagMedium;
    TTreeReaderValue<std::vector<int>> * m_elesFlagVeto;
    TTreeReaderValue<std::vector<int>> * m_recoJetsFlavor;
    TTreeReaderValue<std::vector<int>> * m_qgMult;
    TTreeReaderValue<std::vector<int>> * m_muMatchedJetIdx;
    TTreeReaderValue<std::vector<int>> * m_eleMatchedJetIdx;
    TTreeReaderValue<std::vector<int>> * m_looseisoTrksMatchedJetIdx;
    TTreeReaderValue<std::vector<int>> * m_trksForIsoVetoMatchedJetIdx;
    TTreeReaderValue<std::vector<int>> * m_prodJetsNoLep_qgMult;
    TTreeReaderValue<std::vector<int>> * m_genDecayIdxVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayPdgIdVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayMomIdxVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayMomRefVec;
    TTreeReaderValue<std::vector<int>> * m_W_emuVec;
    TTreeReaderValue<std::vector<int>> * m_W_tauVec;
    TTreeReaderValue<std::vector<int>> * m_W_tau_emuVec;
    TTreeReaderValue<std::vector<int>> * m_W_tau_prongsVec;
    TTreeReaderValue<std::vector<int>> * m_W_tau_nuVec;
    TTreeReaderValue<std::vector<int>> * m_selPDGid;
    TTreeReaderValue<std::vector<int>> * m_trksForIsoVeto_pdgId;
    TTreeReaderValue<std::vector<int>> * m_trksForIsoVeto_idx;
    TTreeReaderValue<std::vector<int>> * m_loose_isoTrks_pdgId;
    TTreeReaderValue<std::vector<int>> * m_loose_isoTrks_idx;
    TTreeReaderValue<std::vector<int>> * m_forVetoIsoTrksidx;
    TTreeReaderValue<std::vector<unsigned int>> * m_elesisEB;
    TTreeReaderValue<std::vector<std::string>> * m_TriggerNames;
    TTreeReaderValue<std::vector<std::string>> * m_genDecayStrVec;
    TTreeReaderValue<std::vector<std::string>> * m_ntpVersion;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_muonsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_specialFixMuonsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_elesLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_jetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_puppiJetsLVec;
    TTreeReaderValue<std::vector<std::vector<TLorentzVector>>> * m_puppiSubJetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_ak8JetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_ak8SubJetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_jetsLVecLepCleaned;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_prodJetsNoLep_puppiJetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_prodJetsNoLep_puppiSubJetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_genDecayLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_selGenParticle;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_genjetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_trksForIsoVetoLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_loose_isoTrksLVec;
};

#endif
