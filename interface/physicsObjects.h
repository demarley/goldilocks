#ifndef PHYSICSOBJECTS_H_
#define PHYSICSOBJECTS_H_

/* 
   Physics objects to be used in analyses
   This structure allows the Event class
   and other classes to access these objects
   without circular inclusion (which breaks!)
*/
#include "TLorentzVector.h"
#include <map>
#include <string>


// base object (consistent reference to TLorentzVector)
struct CmaBase {
    TLorentzVector p4;
};


// Truth information
struct Parton : CmaBase {
    int pdgId;
    int index;       // index in vector of truth partons
    int decayIdx;    // index in truth record
    int parent_ref;  // index in truth vector of parent
    int parent_idx;  // index in truth record of parent
    int top_index;   // index in truth_tops if this is a top
    int containment; // record value used to calculate containment

    // Heavy Object Booleans
    bool isTop;
    bool isW;
    bool isZ;
    bool isHiggs;
    // Lepton Booleans
    bool isLepton;
    bool isTau;
    bool isElectron;
    bool isMuon;
    bool isNeutrino;
    // Quark Booleans
    bool isQuark;
    bool isBottom;
    bool isLight;
};

struct TruthTop {
    // collect indices in truth_partons vector of top parton info
    bool isTop;
    bool isAntiTop;
    int Top;
    int W;
    int bottom;
    std::vector<int> Wdecays;   // for storing W daughters
    std::vector<int> daughters; // for storing non-W/bottom daughters

    bool isHadronic;  // W decays to quarks
    bool isLeptonic;  // W decays to leptons
};


// Reco information
struct Jet : CmaBase {
    // extra jet attributes
    float bdisc;
    std::map<std::string, bool> isbtagged;
    int true_flavor;
    float radius;
    double charge;
    double rho;      // jet energy density, 1 value per event (attaching to each jet for convenience)
    int index;       // index in vector of jets

    int truth_jet;   // index in vector of truth jets that is closest to this jet
    int containment; // level of containment for partons
    std::vector<int> truth_partons;  // vector containing partons that are truth-matched to jet
    int matchId;    // keep track of jets matched to top or anti-top

    float deepCSVb;
    float deepCSVbb;
    float deepCSVc;
    float deepCSVcc;
    float deepCSVl;

    float deepFlavorb;
    float deepFlavorbb;
    float deepFlavorc;
    float deepFlavoruds;
    float deepFlavorg;
    float deepFlavorlepb;
};

struct Ljet : Jet {
    // extra ljet attributes
    int isGood;
    float softDropMass;
    float tau1;
    float tau2;
    float tau3;
    float tau21;
    float tau32;

    double deepAK8top;
    double deepAK8W;
    double deepAK8Z;
    double deepAK8Zbb;
    double deepAK8Hbb;
    double deepAK8H4q;
    std::vector<double> deepAK8; // deepAK8 raw scores
    std::vector<Jet> subjets;    // soft-drop subjets
};



struct Top {
    // Define a top quark
    TLorentzVector p4;
    unsigned int target;        // for ML training
    std::map<std::string,double> dnn;

    bool isTop;
    bool isAntiTop;

    // contains all associated jets (hadronic or leptonic)
    std::vector<int> jets;  

    // hadronically-decaying top quark
    int ljet;               // large-R jet
};



// ------------------------ // 
// Struct to contain sample information (processing the input file)

struct Sample {
    std::string primaryDataset;
    float XSection;
    float KFactor;
    float sumOfWeights;
    unsigned int NEvents;
};

#endif
