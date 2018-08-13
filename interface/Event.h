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
    bool isValidRecoEntry() const {return (m_entry > (long long)-1);}

    // Execute the event (load information and setup objects)
    virtual void execute(Long64_t entry);
    virtual void updateEntry(Long64_t entry);

    // Setup physics information
    void initialize_jets();
    void initialize_ljets();
    void initialize_weights();
    void initialize_truth();

    // Clear stuff;
    virtual void finalize();
    virtual void clear();

    // Get physics object(s) information
    std::vector<Ljet> ljets() const {return m_ljets;}
    std::vector<Ljet> ljetsPUPPI() const {return m_ljetsPUPPI;}
    std::vector<Jet> jets() const {return m_jets;}

    // Get truth physics information 
    TruthTop m_truth_top;
    TruthTop m_truth_antitop;

    // metadata
    long long entry() {return m_entry;}
    virtual std::string treeName() {return m_treeName;}

    // functions to calculate things
    float DNN(){ return m_DNN;}

    void buildTtbar();
    std::vector<Top> ttbar() {return m_ttbar;}
    void deepLearningPrediction(Top& top);

    // Get weights
    virtual float nominal_weight() const {return m_nominal_weight;}

  protected:

    // general information
    configuration *m_config;
    TTreeReader &m_ttree;
    TTreeReader m_truth_tree;
    std::string m_treeName;
    bool m_isMC;
    long long m_entry;
    long long m_truth_entry;
    bool m_DNNtraining;
    bool m_DNNinference;

    // event weight information
    double m_nominal_weight;

    // physics object information
    std::vector<Ljet> m_ljets;
    std::vector<Ljet> m_ljetsPUPPI;
    std::vector<Jet>  m_jets;

    // truth physics object information
    std::vector<Parton> m_truth_partons;
    std::vector<TruthTop> m_truth_tops;
    float m_truthDR = 0.6;  // truth-matching radius

    std::map<std::string, double> m_dnnInputs;   // values for inputs to the DNN
    std::string m_dnnKey;
    float m_DNN;   // DNN score

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


    // TTree variables [all possible ones]
    // *************
    // the following are from root files accessed 
    //    on 5 August 2018
    //    at root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Top_ntuple_V9
    TTreeReaderValue<double> * m_stored_weight;

    TTreeReaderValue<std::vector<int>> * m_PassTrigger;
    TTreeReaderValue<std::vector<std::string>> * m_TriggerNames;

    TTreeReaderValue<int> * m_globalTightHalo2016Filter;
    TTreeReaderValue<int> * m_goodVerticesFilter;
    TTreeReaderValue<int> * m_eeBadScFilter;
    TTreeReaderValue<int> * m_EcalDeadCellTriggerPrimitiveFilter;
    TTreeReaderValue<unsigned int> * m_BadChargedCandidateFilter;
    TTreeReaderValue<unsigned int> * m_BadPFMuonFilter;
    TTreeReaderValue<unsigned int> * m_HBHENoiseFilter;
    TTreeReaderValue<unsigned int> * m_HBHEIsoNoiseFilter;
    TTreeReaderValue<int> * m_duplicateMuonsFilter;
    TTreeReaderValue<int> * m_badMuonsFilter;
    TTreeReaderValue<int> * m_noBadMuonsFilter;

    TTreeReaderValue<int> * m_NJetsISR;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_ak4LVec;
    TTreeReaderValue<unsigned int> * m_ak4looseJetID;
    TTreeReaderValue<unsigned int> * m_ak4tightJetID;
    TTreeReaderValue<unsigned int> * m_ak4tightlepvetoJetID;
    TTreeReaderValue<unsigned int> * m_ak4looseJetID_NoLep;
    TTreeReaderValue<unsigned int> * m_ak4tightJetID_NoLep;
    TTreeReaderValue<unsigned int> * m_ak4tightlepvetoJetID_NoLep;
    TTreeReaderValue<std::vector<double>> * m_ak4qgLikelihood;
    TTreeReaderValue<std::vector<double>> * m_ak4qgPtD;
    TTreeReaderValue<std::vector<double>> * m_ak4qgAxis1;
    TTreeReaderValue<std::vector<double>> * m_ak4qgAxis2;
    TTreeReaderValue<std::vector<int>> * m_ak4qgMult;
    TTreeReaderValue<std::vector<int>> * m_ak4Flavor;
    TTreeReaderValue<std::vector<double>> * m_ak4Charge;
    TTreeReaderValue<std::vector<double>> * m_ak4deepCSV_b;
    TTreeReaderValue<std::vector<double>> * m_ak4deepCSV_bb;
    TTreeReaderValue<std::vector<double>> * m_ak4deepCSV_l;
    TTreeReaderValue<std::vector<double>> * m_ak4deepCSV_c;
    TTreeReaderValue<std::vector<double>> * m_ak4deepCSV_cc;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_b;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_bb;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_lepb;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_c;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_uds;
    TTreeReaderValue<std::vector<double>> * m_ak4deepFlavor_g;


    TTreeReaderValue<std::vector<TLorentzVector>> * m_ak8puppiJetsLVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_ak8DeepAK8LVec;

    TTreeReaderValue<std::vector<double>> * m_ak8puppiTau1;
    TTreeReaderValue<std::vector<double>> * m_ak8puppiTau2;
    TTreeReaderValue<std::vector<double>> * m_ak8puppiTau3;
    TTreeReaderValue<std::vector<double>> * m_ak8puppiSoftDropMass;
    TTreeReaderValue<std::vector<double>> * m_ak8SubJetsBdisc;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8JetsDeepAK8;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_puppiSubJetsBdisc;

    TTreeReaderValue<std::vector<std::vector<TLorentzVector>>> * m_ak8puppiSubJetsLVec;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8puppiSubJetsBdisc;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8puppiSubJetMult;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8puppiSubjetPtD;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8puppiSubjetAxis1;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8puppiSubjetAxis2;

    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8top;
    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8W;
    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8Z;
    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8Zbb;
    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8Hbb;
    TTreeReaderValue<std::vector<double>> * m_ak8DeepAK8H4q;
    TTreeReaderValue<std::vector<std::vector<double>>> * m_ak8DeepAK8;


    TTreeReaderValue<std::vector<int>> * m_selPDGid;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_genDecayLVec;
    TTreeReaderValue<std::vector<double>> * m_genMatched;
    TTreeReaderValue<std::vector<int>> * m_genDecayIdxVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayPdgIdVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayMomIdxVec;
    TTreeReaderValue<std::vector<int>> * m_genDecayMomRefVec;
    TTreeReaderValue<std::vector<TLorentzVector>> * m_selGenParticle;
};

#endif
