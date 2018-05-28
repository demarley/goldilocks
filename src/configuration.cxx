/*
Created:        --
Last Updated:   28 May 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Configuration class
  -- Read config file and use functions
     to return configurations later
*/
#include "Analysis/goldilocks/interface/configuration.h"
#include <stdlib.h>


configuration::configuration(const std::string &configFile) : 
  m_configFile(configFile),
  m_isMC(false),
  m_isQCD(false),
  m_isTtbar(false),
  m_isGridFile(false),
  m_useTruth(false),
  m_useJets(false),
  m_useLeptons(false),
  m_useLargeRJets(false),
  m_useRCJets(false),
  m_useNeutrinos(false),
  m_buildNeutrinos(true),
  m_input_selection("SetMe"),
  m_selection("SetMe"),
  m_cutsfile("SetMe"),
  m_treename("SetMe"),
  m_filename("SetMe"),
  m_primaryDataset("SetMe"),
  m_NTotalEvents(0),
  m_verboseLevel("SetMe"),
  m_nEventsToProcess(0),
  m_firstEvent(0),
  m_outputFilePath("SetMe"),
  m_customFileEnding("SetMe"),
  m_makeNewFile(false),
  m_makeHistograms(false),
  m_cma_absPath("SetMe"),
  m_metadataFile("SetMe"),
  m_getDNN(false),
  m_dnnFile("SetMe"),
  m_dnnKey("SetMe"),
  m_doRecoEventLoop(false),
  m_matchTruthToReco(true),
  m_jet_btag_wkpt("SetMe"),
  m_calcWeightSystematics(false),
  m_listOfWeightSystematicsFile("SetMe"),
  m_listOfWeightVectorSystematicsFile("SetMe"){
    m_XSection.clear();
    m_KFactor.clear();
    m_AMI.clear();
    m_map_config.clear();
  }

configuration::~configuration() {}

configuration &configuration::operator=(const configuration &rhs) { return *this; }

void configuration::initialize() {
    /* Initialize the configurations */
    std::vector<std::string> configurations; 
    cma::read_file( m_configFile, configurations ); // read config file into vector

    // fill map with values from configuration file
    for (const auto& config : configurations){
        // split config items by space
        std::istringstream cfg(config);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_map_config.insert( std::pair<std::string,std::string>(tokens.at(0),tokens.at(1)) );
    }

    // Protection against default settings missing in custom configuration
    // -- map of defaultConfigs defined in header (can't use 'verbose' tools, not defined yet!)
    for (const auto& defaultConfig : m_defaultConfigs){
        if ( m_map_config.find(defaultConfig.first) == m_map_config.end() ){ // item isn't in config file
            std::cout << " WARNING :: CONFIG : Configuration " << defaultConfig.first << " not defined" << std::endl;
            std::cout << " WARNING :: CONFIG : Setting value to default " << defaultConfig.second << std::endl;
            m_map_config[defaultConfig.first] = defaultConfig.second;
        }
    }


    // Set the verbosity level (the amount of output to the console)
    std::map<std::string,unsigned int> verboseMap = cma::verboseMap(); // load mapping of string to integer
    m_verboseLevel = getConfigOption("verboseLevel");
    if (verboseMap.find(m_verboseLevel)==verboseMap.end()){
        m_verboseLevel = "INFO";
        cma::setVerboseLevel(m_verboseLevel);

        cma::WARNING( "CONFIG : Verbose level selected, "+m_verboseLevel+", is not supported " );
        cma::WARNING( "CONFIG : Please select one of the following: " );
        for (const auto& dm : verboseMap)
            cma::WARNING( "CONFIG :          "+dm.first);
        cma::WARNING( "CONFIG : Continuing; setting verbose level to "+m_verboseLevel);
    }
    else{
        cma::setVerboseLevel(m_verboseLevel);
    }


    // Get the absolute path to goldilocks for loading
    char* cma_path = getenv("CYMINIANADIR");
    if (cma_path==NULL){
        cma::WARNING("CONFIG : environment variable " );
        cma::WARNING("CONFIG :    'CYMINIANADIR' " );
        cma::WARNING("CONFIG : is not set.  Using PWD to set path." );
        cma::WARNING("CONFIG : This may cause problems submitting batch jobs." );
        cma_path = getenv("PWD");
    }
    m_cma_absPath = cma_path;
    cma::DEBUG("CONFIG : path set to: "+m_cma_absPath );

    // Assign values
    m_nEventsToProcess = std::stoi(getConfigOption("NEvents"));
    m_firstEvent       = std::stoi(getConfigOption("firstEvent"));
    m_input_selection  = getConfigOption("input_selection"); // "grid", "pre", etc.
    m_selection        = getConfigOption("selection");

    // check that b-tag and top-tag WPs are recognized as one of supported values
    check_btag_WP(getConfigOption("jet_btag_wkpt"));

    m_jet_btag_wkpt    = getConfigOption("jet_btag_wkpt");
    m_outputFilePath   = getConfigOption("output_path");
    m_customFileEnding = getConfigOption("customFileEnding");
    m_cutsfile         = getConfigOption("cutsfile");
    m_useTruth         = cma::str2bool( getConfigOption("useTruth") );
    m_useJets          = cma::str2bool( getConfigOption("useJets") );
    m_useLeptons       = cma::str2bool( getConfigOption("useLeptons") );
    m_useLargeRJets    = cma::str2bool( getConfigOption("useLargeRJets") );
    m_useRCJets        = cma::str2bool( getConfigOption("useRCJets") );
    m_useNeutrinos     = cma::str2bool( getConfigOption("useNeutrinos") );
    m_makeNewFile      = cma::str2bool( getConfigOption("makeNewFile") );
    m_makeHistograms   = cma::str2bool( getConfigOption("makeHistograms") );
    m_makeEfficiencies = cma::str2bool( getConfigOption("makeEfficiencies") );
    m_dnnFile          = getConfigOption("dnnFile");
    m_dnnKey           = getConfigOption("dnnKey");
    m_getDNN           = cma::str2bool( getConfigOption("getDNN") );
    m_doRecoEventLoop  = cma::str2bool( getConfigOption("doRecoEventLoop") );
    m_matchTruthToReco = true;  // not needed in this analysis (so it's not a config option) but here in case we do later
    m_buildNeutrinos   = cma::str2bool( getConfigOption("buildNeutrinos") );
    m_metadataFile     = getConfigOption("metadataFile");
    m_calcWeightSystematics             = cma::str2bool( getConfigOption("calcWeightSystematics") );
    m_listOfWeightSystematicsFile       = getConfigOption("weightSystematicsFile");
    m_listOfWeightVectorSystematicsFile = getConfigOption("weightVectorSystematicsFile");

    cma::read_file( getConfigOption("inputfile"), m_filesToProcess );
    //cma::read_file( getConfigOption("treenames"), m_treeNames );
    m_treename = getConfigOption("treename");

    // Type of input file
    // if grid file, do some extra calculations
    m_isGridFile = (m_input_selection.compare("grid")==0);

    m_mapOfSamples.clear();
    cma::getSampleWeights( m_metadataFile,m_mapOfSamples );

    // systematics that are weights in the nominal tree
    m_listOfWeightSystematics.resize(0);
    cma::read_file(m_listOfWeightSystematicsFile,m_listOfWeightSystematics);

    std::vector<std::string> weightVectorSystematics; // vector btagging SF
    cma::read_file(m_listOfWeightVectorSystematicsFile,weightVectorSystematics);

    for (const auto& weightVectorSystematic : weightVectorSystematics){
        // split config items by space
        std::istringstream cfg(weightVectorSystematic);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_mapOfWeightVectorSystematics.insert( std::pair<std::string,unsigned int>( tokens.at(0),std::stoi(tokens.at(1)) ) );
    }

    return;
}


void configuration::print(){
    // -- Print the configuration
    std::cout << " ** CyMiniAna ** " << std::endl;
    std::cout << " --------------- " << std::endl;
    std::cout << " CONFIGURATION :: Printing configuration " << std::endl;
    std::cout << " " << std::endl;
    for (const auto& config : m_map_config){
        std::cout << " " << config.first << "\t\t\t" << config.second << std::endl;
    }
    std::cout << " --------------- " << std::endl;

    return;
}


std::string configuration::getConfigOption( std::string item ){
    /* Check that the item exists in the map & return it; otherwise throw exception  */
    std::string value("");

    try{
        value = m_map_config.at(item);
    }
    catch(const std::exception&){
        cma::ERROR("CONFIG : Option "+item+" does not exist in configuration.");
        cma::ERROR("CONFIG : This does not exist in the default configuration either.");
        cma::ERROR("CONFIG : Returing an empty string.");
    }

    return value;
}


void configuration::inspectFile( TFile& file ){
    /* Compare filenames to determine file type */
    m_isQCD   = false;
    m_isTtbar = false;
    m_primaryDataset = "";
    m_NTotalEvents = 0;

    // check if file is QCD
    for (const auto& qcd : m_qcdFiles){
        m_isQCD = (m_filename.find(qcd)!=std::string::npos);
        m_primaryDataset = qcd;
        if (m_isQCD) break;
    }

    // check if file is ttbar
    for (const auto& ttbar : m_ttbarFiles){
        m_isTtbar = (m_filename.find(ttbar)!=std::string::npos);
        m_primaryDataset = ttbar;
        if (m_isTtbar) break;
    }

    m_isMC = (m_isQCD || m_isTtbar);  // only considering ttbar vs qcd

    // get the metadata
    if (m_primaryDataset.size()>0) m_NTotalEvents = m_mapOfSamples.at(m_primaryDataset).NEvents;
    else{
        cma::WARNING("CONFIGURATION : Primary dataset name not found, checking the map");
        for (const auto& s : m_mapOfSamples){
            Sample samp = s.second;
            std::size_t found = m_filename.find(samp.primaryDataset);
            if (found!=std::string::npos){ 
                m_primaryDataset = samp.primaryDataset;
                m_NTotalEvents   = samp.NEvents;
                break;
            } // end if the filename contains this primary dataset
        } // end loop over map of samples (to access metadata info)
    } // end else


    // Protection against accessing truth information that may not exist
    if (!m_isMC && m_useTruth){
        cma::WARNING("CONFIGURATION : 'useTruth=true' but 'isMC=false'");
        cma::WARNING("CONFIGURATION : Setting 'useTruth' to false");
        m_useTruth = false;
    }

    return;
}


std::string configuration::verboseLevel(){
    /* Return the verbosity level */
    return m_verboseLevel;
}

std::string configuration::getAbsolutePath(){
    /* Return the absolute path of the goldilocks directory (batch jobs) */
    return m_cma_absPath;
}

void configuration::setTreename(std::string treeName){
    m_treename = treeName;
    return;
}

std::string configuration::treename(){
    return m_treename;
}

void configuration::setFilename(std::string fileName){
    m_filename = fileName;
    return;
}

std::string configuration::filename(){
    return m_filename;
}

bool configuration::isNominalTree(){
    return isNominalTree( m_treename );
}

bool configuration::isNominalTree( const std::string &tree_name ){
    /* Check if tree is a nominal one */
    bool isNominal(false);
    if (tree_name.compare("nominal")==0 || tree_name.compare("Nominal")==0)
        isNominal = true;
    else
        isNominal = false;

    return isNominal;
}

bool configuration::isMC(){
    return m_isMC;
}

bool configuration::isMC( TFile& file ){
    /* Check the sum of weights tree DSIDs (to determine Data || MC) */
    inspectFile( file );
    return m_isMC;
}

bool configuration::calcWeightSystematics(){
    return m_calcWeightSystematics;
}

std::map<std::string,unsigned int> configuration::mapOfWeightVectorSystematics(){
    return m_mapOfWeightVectorSystematics;
}

std::vector<std::string> configuration::listOfWeightSystematics(){
    return m_listOfWeightSystematics;
}

std::string configuration::listOfWeightSystematicsFile(){
    return m_listOfWeightSystematicsFile;
}

std::string configuration::listOfWeightVectorSystematicsFile(){
    return m_listOfWeightVectorSystematicsFile;
}

std::string configuration::metadataFile(){
    return m_metadataFile;
}

void configuration::check_btag_WP(const std::string &wkpt){
    if(! std::any_of(m_btag_WPs.begin(), m_btag_WPs.end(), [&](const std::string& s){return (s.compare(wkpt) == 0);} ) ) {
        cma::ERROR("CONFIG : Unknown b-tagging WP: "+wkpt+". Aborting!");
        cma::ERROR("CONFIG : Available b-tagging WPs: "+cma::vectorToStr(m_btag_WPs));
        exit(EXIT_FAILURE);
    }

    return;
}

bool configuration::useJets(){
    return m_useJets;
}

bool configuration::useLeptons(){
    return m_useLeptons;
}

bool configuration::useLargeRJets(){
    return m_useLargeRJets;
}

bool configuration::useRCJets(){
    return m_useRCJets;
}

bool configuration::useNeutrinos(){
    return m_useNeutrinos;
}

std::string configuration::configFileName(){
    return m_configFile;
}

std::string configuration::selection(){
    return m_selection;
}

std::string configuration::cutsfile(){
    return m_cutsfile;
}

std::string configuration::jet_btagWkpt(){
    return m_jet_btag_wkpt;
}

int configuration::nEventsToProcess(){
    return m_nEventsToProcess;
}

unsigned long long configuration::firstEvent(){
    return m_firstEvent;
}

std::string configuration::outputFilePath(){
    return m_outputFilePath;
}

std::string configuration::customFileEnding(){
    return m_customFileEnding;
}

double configuration::LUMI(){
    return m_LUMI;
}

bool configuration::useTruth(){
    return m_useTruth;
}

bool configuration::isGridFile(){
    return m_isGridFile;
}

std::vector<std::string> configuration::filesToProcess(){
    return m_filesToProcess;
}

std::vector<std::string> configuration::treeNames(){
    return m_treeNames;
}

bool configuration::makeNewFile(){
    return m_makeNewFile;
}

bool configuration::makeHistograms(){
    return m_makeHistograms;
}

bool configuration::makeEfficiencies(){
    return m_makeEfficiencies;
}


// values for DNN
std::string configuration::dnnFile(){
    return m_dnnFile;  // file with dnn values (*.json)
}
std::string configuration::dnnKey(){
    return m_dnnKey;   // key in dnnFile to access NN
}
// build the DNN
bool configuration::getDNN(){
    return m_getDNN;
}
bool configuration::doRecoEventLoop(){
    return m_doRecoEventLoop;
}

bool configuration::matchTruthToReco(){
    /* true -- match truth events to reco events (loop over reco events)
     * false -- match reco events to truth events (loop over truth events)
     */
    return m_matchTruthToReco;
}

void configuration::setMatchTruthToReco(bool truthToReco){
    m_matchTruthToReco = truthToReco;
}

bool configuration::buildNeutrinos(){
    return m_buildNeutrinos;
}

// the end
