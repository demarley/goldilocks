#ifndef MINITREE_H_
#define MINITREE_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Analysis/goldilocks/interface/Event.h"
#include "Analysis/goldilocks/interface/physicsObjects.h"
#include "Analysis/goldilocks/interface/eventSelection.h"
#include "Analysis/goldilocks/interface/configuration.h"


class miniTree {
  public:
    // Default - so root can load based on a name;
    miniTree(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~miniTree();

    // Run once at the start of the job;
    virtual void initialize(TFile& outputFile);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveEvent(const std::map<std::string,double> features);

    // Clear stuff;
    virtual void finalize();


  protected:

    TTree * m_ttree;
    TTree * m_metadataTree;
    configuration * m_config;

    bool m_usePUPPI;

    /**** Training branches ****/
    // weights for inputs
    float m_xsection;
    float m_kfactor;
    float m_sumOfWeights;
    float m_weight;
    float m_nominal_weight;

    // Deep learning features
    unsigned int m_target;

    unsigned int m_nDeepAK8;
    std::vector<float> m_AK8_deepAK8;

    float m_AK8_pt;
    float m_AK8_SDmass;
    float m_AK8_tau1;
    float m_AK8_tau2;
    float m_AK8_tau3;
    float m_AK8_tau21;
    float m_AK8_tau32;
    float m_AK8_subjet0_bdisc;
    float m_AK8_subjet0_pTrel;
    float m_AK8_subjet0_charge;
    float m_AK8_subjet1_bdisc;
    float m_AK8_subjet1_pTrel;
    float m_AK8_subjet1_charge;

    float m_AK4_deepCSVb;
    float m_AK4_deepCSVbb;
    float m_AK4_deepCSVc;
    float m_AK4_deepCSVcc;
    float m_AK4_deepCSVl;
    float m_AK4_deepFlavorb;
    float m_AK4_deepFlavorbb;
    float m_AK4_deepFlavorc;
    float m_AK4_deepFlavoruds;
    float m_AK4_deepFlavorg;
    float m_AK4_deepFlavorlepb;

    float m_AK4_charge;

    float m_AK8AK4_mass;
    float m_AK8AK4_deltaR;

    /**** Metadata ****/
    // which sample has which target value
    // many ROOT files will be merged together to do the training!
    std::string m_name;
    unsigned int m_target_value;
    unsigned int m_nEvents;
};

#endif
