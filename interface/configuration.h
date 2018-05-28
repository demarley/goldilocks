#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>

#include "Analysis/goldilocks/interface/tools.h"


class configuration {
  public:
    // Default - so root can load based on a name;
    configuration( const std::string &configFile );
    //configuration( const configuration& );
    configuration& operator=( const configuration& rhs );

    // Default - so we can clean up;
    virtual ~configuration();

    // Run once at the start of the job;
    virtual void initialize();
    std::string getConfigOption( std::string item );

    // Print configuration
    virtual void print();

    virtual bool isMC();              // must call "inspectFile(file)" or "isMC(file)" first!
    virtual bool isMC( TFile& file );
    bool isGridFile();
    bool isTtbar(){ return m_isTtbar;}
    bool isQCD(){ return m_isQCD;}

    // object declarations
    virtual bool useJets();
    virtual bool useNeutrinos();
    virtual bool useLeptons();
    virtual bool useLargeRJets();
    virtual bool useRCJets();
    virtual bool useTruth();

    std::string jet_btagWkpt();
    float cMVAv2L() {return m_cMVAv2L;}
    float cMVAv2M() {return m_cMVAv2M;}
    float cMVAv2T() {return m_cMVAv2T;}

    // functions about the TTree
    virtual bool isNominalTree();
    virtual bool isNominalTree( const std::string &tree_name );
    std::vector<std::string> treeNames();
    void setTreename(std::string treeName);
    std::string treename();

    // functions about the file
    virtual void inspectFile( TFile& file );
    std::vector<std::string> filesToProcess();
    void setFilename(std::string fileName);
    std::string filename();
    std::string primaryDataset() {return m_primaryDataset;}
    unsigned int NTotalEvents() {return m_NTotalEvents;}

    // return some values from config file
    std::string verboseLevel();
    std::string selection();
    std::string cutsfile();
    std::string outputFilePath();
    std::string customFileEnding();
    std::string configFileName();
    std::string getAbsolutePath();
    int nEventsToProcess();
    unsigned long long firstEvent();
    bool makeNewFile();
    bool makeHistograms();
    bool makeEfficiencies();

    // information for event weights
    std::string metadataFile();
    std::map<std::string,Sample> mapOfSamples(){return m_mapOfSamples;}
    Sample sample(){return m_mapOfSamples.at(m_primaryDataset);}
    virtual double LUMI();

    // weight systematics
    bool calcWeightSystematics();
    std::map<std::string,unsigned int> mapOfWeightVectorSystematics();
    std::vector<std::string> listOfWeightSystematics();
    std::string listOfWeightSystematicsFile();
    std::string listOfWeightVectorSystematicsFile();

    // DNN
    std::string dnnFile();
    bool getDNN();
    double minDNN();
    double maxDNN();
    std::string dnnKey();   // key for lwtnn

    // Reco/Truth event loops
    bool doRecoEventLoop();
    bool doTruthEventLoop();
    bool matchTruthToReco();
    void setMatchTruthToReco(bool truthToReco);
    std::map<std::string,int> mapOfPartonContainment() {return m_containmentMap;}
    std::map<int,std::string> mapOfPartonContainmentRev() {return m_containmentMapRev;}
    std::map<std::string,int> mapOfTargetValues() {return m_targetMap;}

    // misc. for dilepton ttbar
    bool buildNeutrinos();

    float beamEnergy() {return m_beamEnergy;}           // 13000.;
    double topQuarkMass() {return m_topQuarkMass;}      // 172.5
    double bQuarkMass() {return m_bQuarkMass;}          // 4.18
    double WMass() {return m_WMass;}                    // 80.2

    /// All analysis eras as needed
    enum Era{run2_13tev_25ns,     run2_13tev_2015_25ns, run2_13tev_2016_25ns, 
             run2_13tev_25ns_74X, undefined};
    Era convert(const std::string& era);       /// Convert an era from string to enum
    std::string convert(const Era& era);   /// Convert an era from enum to string
    double energyInTev(const Era& era) {return 13.;}     /// Return energy for given era in TeV

  protected:

    void check_btag_WP(const std::string &wkpt);

    std::map<std::string,std::string> m_map_config;
    const std::string m_configFile;

    bool m_isMC;
    bool m_isQCD;
    bool m_isTtbar;
    bool m_isGridFile;
    bool m_useTruth;

    // object declarations
    bool m_useJets;
    bool m_useLeptons;
    bool m_useLargeRJets;
    bool m_useRCJets;
    bool m_useNeutrinos;
    bool m_buildNeutrinos;

    // luminosity
    double m_LUMI      = 36074.56; // 2015+2016 luminosity
    double m_LUMI_2015 = 3212.96;
    double m_LUMI_2016 = 32861.6; // OflLumi-13TeV-008

    // return some values from config file
    std::string m_input_selection;
    std::string m_selection;
    std::string m_cutsfile;
    std::string m_treename;
    std::string m_filename;
    std::string m_primaryDataset;
    unsigned int m_NTotalEvents;
    std::string m_verboseLevel;
    int m_nEventsToProcess;
    unsigned long long m_firstEvent;
    std::string m_outputFilePath;
    std::string m_customFileEnding;
    bool m_makeNewFile;
    bool m_makeHistograms;
    bool m_makeEfficiencies;
    std::string m_cma_absPath;
    std::string m_metadataFile;
    bool m_getDNN;
    std::string m_dnnFile;
    std::string m_dnnKey;
    bool m_doRecoEventLoop;
    bool m_matchTruthToReco;

    std::string m_jet_btag_wkpt;   // "L","M","T"
    std::string m_tjet_btag_wkpt;
    std::string m_toptag_wkpt;

    // b-tagging (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco)
    // isBTagged = (jet.cMVAv2 > wkpt)
    std::vector<std::string> m_btag_WPs = {"L","M","T"};
    float m_cMVAv2L=-0.5884;
    float m_cMVAv2M=0.4432;
    float m_cMVAv2T=0.9432;

    std::vector<std::string> m_filesToProcess;
    std::vector<std::string> m_treeNames;

    bool m_calcWeightSystematics;
    std::map<std::string,unsigned int> m_mapOfWeightVectorSystematics;
    std::vector<std::string> m_listOfWeightSystematics;
    std::string m_listOfWeightSystematicsFile;
    std::string m_listOfWeightVectorSystematicsFile;

    std::map<std::string,Sample> m_mapOfSamples;  // map of Sample structs

    // OLD:
    std::map<std::string, float> m_XSection; // map DSID to XSection
    std::map<std::string, float> m_KFactor;  // map DSID to KFactor
    std::map<std::string, float> m_AMI;      // map DSID to sum of weights

    // -- Top Mass Variables -- //
    const double m_electronMass = 0.000511;
    const double m_muonMass     = 0.105658;
    const double m_bQuarkMass   = 4.8;
    const double m_WMass        = 80.4;
    const double m_topQuarkMass = 172.5;
    const float m_beamEnergy    = 13000.;
    const int SENTINEL    = -1000;
    const int NCHAN       = 4;
    const double m_sqrt_s = 13000;      // center-of-mass energy

    // Samples primary dataset names
    std::vector<std::string> m_qcdFiles   = {
      "QCD_HT1000to1500_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT1500to2000_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT2000toInf_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT200to300_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT300to500_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT500to700_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT700to1000_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",  
      "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"}; // possible qcd files
    std::vector<std::string> m_ttbarFiles = {
      "TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8",
      "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
      "TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"}; // possible Ttbar files

    // Degrees of 'containment' for parton matching to jets
    std::map<std::string,int> m_containmentMap = {
                    {"NONE",  0},
                    {"BONLY", 1},
                    {"QONLY", 2},
                    {"BQ",    3},   // B+Q
                    {"W",     4},   // Q+Q
                    {"FULL",  5}};  // (BQ)+Q || (W+B) || Q+Q+B
    std::map<int,std::string> m_containmentMapRev = {
                    {0,"NONE"},
                    {1,"BONLY"},
                    {2,"QONLY"},
                    {3,"BQ"},     // B+Q
                    {4,"W"},      // Q+Q
                    {5,"FULL"}};  // (BQ)+Q || (W+B) || Q+Q+B
    // map of target values used in training ML
    std::map<std::string,int> m_targetMap = {
      {"none",0},    // :: NONE = QCD (background)
      {"QB",1},      // :: QB-Q = Signal AK8(QB) + AK4(Q)
      {"W",2},       // :: QQ-B = Signal AK8(W)  + AK4(B)
      {"full",3},    // :: FULL = Signal AK8
      {"other",4} }; // :: Other = placeholder (resolved/Q-only/B-only/etc.)


    std::map<std::string,std::string> m_defaultConfigs = {
             {"useJets",               "true"},
             {"useLeptons",            "true"},
             {"useLargeRJets",         "true"},
             {"useRCJets",             "false"},
             {"useNeutrinos",          "true"},
             {"useTruth",              "false"},
             {"jet_btag_wkpt",         "70"},
             {"makeNewFile",           "true"},
             {"makeHistograms",        "true"},
             {"makeEfficiencies",      "true"},
             {"NEvents",               "-1"},
             {"firstEvent",            "0"},
             {"input_selection",       "grid"},
             {"selection",             "example"},
             {"output_path",           "./"},
             {"customFileEnding",      ""},
             {"calcWeightSystematics", "false"},
             {"weightSystematicsFile",       "config/weightSystematics.txt"},
             {"weightVectorSystematicsFile", "config/weightVectorSystematics.txt"},
             {"cutsfile",              "examples/config/cuts_example.txt"},
             {"inputfile",             "examples/config/miniSL_ALLfiles.txt"},
             {"treename",              "stopTreeMaker/AUX"},
             {"treenames",             "config/treenames_nominal.txt"},
             {"metadataFile",          "config/sampleMetaData.txt"},
             {"verboseLevel",          "INFO"},
             {"dnnFile",               "config/keras_ttbar_DNN.json"},
             {"dnnKey",                "dnn"},
             {"getDNN",                "false"},
             {"doRecoEventLoop",       "true"},
             {"buildNeutrinos",        "true"} };
};

#endif
