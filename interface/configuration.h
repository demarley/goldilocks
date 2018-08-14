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

    virtual bool isMC();             // must call "inspectFile(file)" or "isMC(file)" first!
    virtual bool isMC( TFile& file );
    bool isTtbar(){ return m_isTtbar;}
    bool isQCD(){ return m_isQCD;}

    // object declarations
    virtual bool useJets(){ return m_useJets;}
    virtual bool useLargeRJets(){ return m_useLargeRJets;}
    virtual bool usePUPPI(){ return m_usePUPPI;}
    virtual bool useTruth(){ return m_useTruth;}

    // functions about the TTree
    void setTreename(std::string treeName);
    std::string treename(){ return m_treename;}

    // functions about the file
    virtual void inspectFile( TFile& file );
    std::vector<std::string> filesToProcess(){ return m_filesToProcess;}
    void setFilename(std::string fileName);
    std::string filename(){ return m_filename;}
    std::string primaryDataset() {return m_primaryDataset;}
    unsigned int NTotalEvents() {return m_NTotalEvents;}

    // return some values from config file
    std::string verboseLevel(){ return m_verboseLevel;}
    std::string selection(){ return m_selection;}
    std::string cutsfile(){ return m_cutsfile;}
    std::string outputFilePath(){ return m_outputFilePath;}
    std::string customFileEnding(){ return m_customFileEnding;}
    std::string configFileName(){ return m_configFile;}
    std::string getAbsolutePath(){ return m_cma_absPath;}
    int nEventsToProcess(){ return m_nEventsToProcess;}
    unsigned long long firstEvent(){ return m_firstEvent;}
    bool makeNewFile(){ return m_makeNewFile;}
    bool makeHistograms(){ return m_makeHistograms;}

    // information for event weights
    std::string metadataFile(){ return m_metadataFile;}
    std::map<std::string,Sample> mapOfSamples(){return m_mapOfSamples;}
    Sample sample(){return m_mapOfSamples.at(m_primaryDataset);}

    // DNN
    std::string dnnFile(){ return m_dnnFile;}
    bool DNNtraining(){ return m_DNNtraining;}
    bool DNNinference(){ return m_DNNinference;}
    std::string dnnKey(){ return m_dnnKey;}   // key for lwtnn

    // Reco/Truth event loops
    std::map<std::string,int> mapOfPartonContainment() {return m_containmentMap;}
    std::map<int,std::string> mapOfPartonContainmentRev() {return m_containmentMapRev;}
    std::map<std::string,int> mapOfTargetValues() {return m_targetMap;}

  protected:

    void check_btag_WP(const std::string &wkpt);

    std::map<std::string,std::string> m_map_config;
    const std::string m_configFile;

    bool m_isMC;
    bool m_isQCD;
    bool m_isTtbar;
    bool m_isGridFile;
    bool m_useTruth;
    bool m_fileInspected;

    // object declarations
    bool m_useJets;
    bool m_useLargeRJets;
    bool m_usePUPPI;

    // return some values from config file
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
    std::string m_cma_absPath;
    std::string m_metadataFile;
    bool m_DNNtraining;
    bool m_DNNinference;
    std::string m_dnnFile;
    std::string m_dnnKey;

    std::vector<std::string> m_filesToProcess;
    std::map<std::string,Sample> m_mapOfSamples;  // map of Sample structs

    // OLD:
    std::map<std::string, float> m_XSection; // map DSID to XSection
    std::map<std::string, float> m_KFactor;  // map DSID to KFactor
    std::map<std::string, float> m_AMI;      // map DSID to sum of weights

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
      {"BQ",1},      // :: QB-Q = Signal AK8(QB) + AK4(Q)
      {"W",2},       // :: QQ-B = Signal AK8(W)  + AK4(B)
      {"ttbckg",3} };// :: QQ+X/QB+X : background
//      {"FULL",4} };  // :: FULL = Signal AK8


    std::map<std::string,std::string> m_defaultConfigs = {
             {"useJets",               "true"},
             {"useLargeRJets",         "true"},
             {"usePUPPI",              "false"},
             {"useTruth",              "false"},
             {"makeNewFile",           "true"},
             {"makeHistograms",        "true"},
             {"NEvents",               "-1"},
             {"firstEvent",            "0"},
             {"selection",             "none"},
             {"output_path",           "./"},
             {"customFileEnding",      ""},
             {"cutsfile",              "config/cuts_none.txt"},
             {"inputfile",             "config/inputfiles.txt"},
             {"treename",              "stopTreeMaker/AUX"},
             {"metadataFile",          "config/sampleMetaData.txt"},
             {"verboseLevel",          "INFO"},
             {"DNNtraining",           "false"},
             {"DNNinference",          "false"},
             {"dnnFile",               "config/keras_ttbar_DNN.json"},
             {"dnnKey",                "dnn"} };
};

#endif
