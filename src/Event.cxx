/*
Created:        24 January 2018
Last Updated:   28 May     2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Event class
 Contains all the objects (& structs) with event information
*/
#include "Analysis/goldilocks/interface/Event.h"

// constructor
Event::Event( TTreeReader &myReader, configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_ttree(myReader),
  m_treeName("SetMe"),
  m_DNN(0.0),
  m_useLargeRJets(false),
  m_useJets(false),
  m_useLeptons(false),
  m_useNeutrinos(false),
  m_useTruth(false){
    m_isMC     = m_config->isMC();
    m_grid     = m_config->isGridFile();             // file directly from original analysis team
    m_treeName = m_ttree.GetTree()->GetName();       // for systematics
    m_getDNN   = m_config->getDNN();                 // build DNN

    m_mapContainment = m_config->mapOfPartonContainment();  // containment map (ints and strings)
    m_targetMap      = m_config->mapOfTargetValues();       // map of target values for NN training

    m_cMVAv2L = m_config->cMVAv2L();
    m_cMVAv2M = m_config->cMVAv2M();
    m_cMVAv2T = m_config->cMVAv2T();

    // Event Info
    m_run   = new TTreeReaderValue<unsigned int>(m_ttree,"run");
    m_lumi  = new TTreeReaderValue<unsigned int>(m_ttree,"lumi");
    m_event = new TTreeReaderValue<unsigned long long>(m_ttree,"event");
    m_avg_npv = new TTreeReaderValue<double>(m_ttree,"avg_npv");

    m_vtxSize = new TTreeReaderValue<int>(m_ttree,"vtxSize");
    m_npv = new TTreeReaderValue<int>(m_ttree,"npv");
    m_nm1 = new TTreeReaderValue<int>(m_ttree,"nm1");
    m_n0  = new TTreeReaderValue<int>(m_ttree,"n0");
    m_np1 = new TTreeReaderValue<int>(m_ttree,"np1");

    m_eeBadScFilter   = new TTreeReaderValue<int>(m_ttree,"eeBadScFilter");
    m_BadPFMuonFilter = new TTreeReaderValue<unsigned int>(m_ttree,"BadPFMuonFilter");
    m_HBHENoiseFilter = new TTreeReaderValue<unsigned int>(m_ttree,"HBHENoiseFilter");
    m_goodVerticesFilter = new TTreeReaderValue<int>(m_ttree,"goodVerticesFilter");
    m_HBHEIsoNoiseFilter = new TTreeReaderValue<unsigned int>(m_ttree,"HBHEIsoNoiseFilter");
    m_globalTightHalo2016Filter = new TTreeReaderValue<int>(m_ttree,"globalTightHalo2016Filter");
    m_BadChargedCandidateFilter = new TTreeReaderValue<unsigned int>(m_ttree,"BadChargedCandidateFilter");
    m_EcalDeadCellTriggerPrimitiveFilter = new TTreeReaderValue<int>(m_ttree,"EcalDeadCellTriggerPrimitiveFilter");

    // Kinematics
    m_met    = new TTreeReaderValue<double>(m_ttree,"met");
    m_metphi = new TTreeReaderValue<double>(m_ttree,"metphi");
    m_calomet    = new TTreeReaderValue<double>(m_ttree,"calomet");
    m_calometphi = new TTreeReaderValue<double>(m_ttree,"calometphi");

    m_PassTrigger  = new TTreeReaderValue<std::vector<int>>(m_ttree,"PassTrigger");
    m_TriggerNames = new TTreeReaderValue<std::vector<std::string>>(m_ttree,"TriggerNames");
    m_TriggerPrescales = new TTreeReaderValue<std::vector<int>>(m_ttree,"TriggerPrescales");

    // Tracks -- not used as of 12 January 2018
    m_nIsoTrks_CUT = new TTreeReaderValue<int>(m_ttree,"nIsoTrks_CUT");
    m_forVetoIsoTrksidx = new TTreeReaderValue<std::vector<int>>(m_ttree,"forVetoIsoTrksidx");

    m_trksForIsoVeto_pdgId  = new TTreeReaderValue<std::vector<int>>(m_ttree,"trksForIsoVeto_pdgId");
    m_trksForIsoVeto_idx    = new TTreeReaderValue<std::vector<int>>(m_ttree,"trksForIsoVeto_idx");
    m_trksForIsoVetoLVec    = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"trksForIsoVetoLVec");
    m_trksForIsoVeto_charge = new TTreeReaderValue<std::vector<double>>(m_ttree,"trksForIsoVeto_charge");
    m_trksForIsoVeto_dz     = new TTreeReaderValue<std::vector<double>>(m_ttree,"trksForIsoVeto_dz");
    m_trksForIsoVeto_iso    = new TTreeReaderValue<std::vector<double>>(m_ttree,"trksForIsoVeto_iso");
    m_trksForIsoVeto_pfActivity = new TTreeReaderValue<std::vector<double>>(m_ttree,"trksForIsoVeto_pfActivity");

    m_loose_nIsoTrks = new TTreeReaderValue<int>(m_ttree,"loose_nIsoTrks");
    m_loose_isoTrks_pdgId  = new TTreeReaderValue<std::vector<int>>(m_ttree,"loose_isoTrks_pdgId");
    m_loose_isoTrks_idx    = new TTreeReaderValue<std::vector<int>>(m_ttree,"loose_isoTrks_idx");
    m_loose_isoTrksLVec    = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"loose_isoTrksLVec");
    m_loose_isoTrks_charge = new TTreeReaderValue<std::vector<double>>(m_ttree,"loose_isoTrks_charge");
    m_loose_isoTrks_dz  = new TTreeReaderValue<std::vector<double>>(m_ttree,"loose_isoTrks_dz");
    m_loose_isoTrks_iso = new TTreeReaderValue<std::vector<double>>(m_ttree,"loose_isoTrks_iso");
    m_loose_isoTrks_mtw = new TTreeReaderValue<std::vector<double>>(m_ttree,"loose_isoTrks_mtw");
    m_loose_isoTrks_pfActivity = new TTreeReaderValue<std::vector<double>>(m_ttree,"loose_isoTrks_pfActivity");


    // JETS
    if (m_config->useLargeRJets()){
        m_useLargeRJets = true;

        m_puppiJetsLVec    = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"puppiAK8LVec");
        m_puppiSubJetsLVec = new TTreeReaderValue<std::vector<std::vector<TLorentzVector>>>(m_ttree,"puppiAK8SubjetLVec");
        m_ak8JetsLVec      = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"deepAK8LVec");
//        m_ak8SubJetsLVec  = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"ak8SubJetsLVec");

        m_qgLikelihood = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgLikelihood");
        m_qgPtD   = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgPtD");
        m_qgAxis2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgAxis2");
        m_qgMult  = new TTreeReaderValue<std::vector<int>>(m_ttree,"qgMult");
/*
        m_tau1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"tau1");
        m_tau2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"tau2");
        m_tau3 = new TTreeReaderValue<std::vector<double>>(m_ttree,"tau3");
        m_softDropMass    = new TTreeReaderValue<std::vector<double>>(m_ttree,"softDropMass");
        m_ak8SubJetsBdisc = new TTreeReaderValue<std::vector<double>>(m_ttree,"ak8SubJetsBdisc");
*/
        m_ak8JetsDeepAK8  = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"deepAK8raw");
        m_puppitau1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau1");
        m_puppitau2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau2");
        m_puppitau3 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau3");
        m_puppisoftDropMass = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8SoftDropMass");
        m_puppiSubJetsBdisc = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetBDisc");
    }

    if (m_config->useJets()){
        m_useJets = true;

        m_NJetsISR = new TTreeReaderValue<int>(m_ttree,"NJetsISR");
        m_jetsLVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"jetsLVec");
        m_recoJetsFlavor = new TTreeReaderValue<std::vector<int>>(m_ttree,"recoJetsFlavor");

        m_looseJetID = new TTreeReaderValue<unsigned int>(m_ttree,"looseJetID");
        m_tightJetID = new TTreeReaderValue<unsigned int>(m_ttree,"tightJetID");
        m_tightlepvetoJetID = new TTreeReaderValue<unsigned int>(m_ttree,"tightlepvetoJetID");
/*
        m_jetsLVecLepCleaned = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"jetsLVecLepCleaned");
        m_looseJetID_NoLep = new TTreeReaderValue<unsigned int>(m_ttree,"looseJetID_NoLep");
        m_tightJetID_NoLep = new TTreeReaderValue<unsigned int>(m_ttree,"tightJetID_NoLep");
        m_tightlepvetoJetID_NoLep = new TTreeReaderValue<unsigned int>(m_ttree,"tightlepvetoJetID_NoLep");
*/
        m_recoJetsJecUnc   = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsJecUnc");
        m_recoJetsBtag_0   = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsBtag_0");
        m_recoJetsCharge_0 = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsCharge_0");
        m_recoJetsJecScaleRawToFull  = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsJecScaleRawToFull");
        m_recoJetsmuonEnergyFraction = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsmuonEnergyFraction");
        m_recoJetschargedEmEnergyFraction = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetschargedEmEnergyFraction");
        m_recoJetsneutralEmEnergyFraction = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsneutralEmEnergyFraction");
        m_recoJetschargedHadronEnergyFraction = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetschargedHadronEnergyFraction");

        m_deepFlavor_b    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorb");
        m_deepFlavor_bb   = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorbb");
        m_deepFlavor_lepb = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorlepb");
        m_deepFlavor_c    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorc");
        m_deepFlavor_uds  = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavoruds");
        m_deepFlavor_g    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorg");
/*
        m_recoJetsJecUncLepCleaned    = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsJecUncLepCleaned");
        m_recoJetsBtag_0_LepCleaned   = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsBtag_0_LepCleaned");
        m_recoJetsCharge_0_LepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsCharge_0_LepCleaned");
        m_recoJetsJecScaleRawToFull_LepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsJecScaleRawToFull_LepCleaned");
        m_recoJetsmuonEnergyFractionLepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsmuonEnergyFractionLepCleaned");
        m_recoJetsneutralEmEnergyFractionLepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsneutralEmEnergyFractionLepCleaned");
        m_recoJetschargedEmEnergyFractionLepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetschargedEmEnergyFractionLepCleaned");
        m_recoJetschargedHadronEnergyFractionLepCleaned = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetschargedHadronEnergyFractionLepCleaned");
        m_prodJetsNoLep_puppiJetsLVec    = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"prodJetsNoLep_puppiJetsLVec");
        m_prodJetsNoLep_puppiSubJetsLVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"prodJetsNoLep_puppiSubJetsLVec");
        m_prodJetsNoLep_qgPtD   = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_qgPtD");
        m_prodJetsNoLep_qgMult  = new TTreeReaderValue<std::vector<int>>(m_ttree,"prodJetsNoLep_qgMult");
        m_prodJetsNoLep_qgAxis2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_qgAxis2");
        m_prodJetsNoLep_qgLikelihood = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_qgLikelihood");
        m_prodJetsNoLep_tau1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_tau1");
        m_prodJetsNoLep_tau2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_tau2");
        m_prodJetsNoLep_tau3 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_tau3");
        m_prodJetsNoLep_puppisoftDropMass = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_puppisoftDropMass");
        m_prodJetsNoLep_puppitau1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_puppitau1");
        m_prodJetsNoLep_puppitau2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_puppitau2");
        m_prodJetsNoLep_puppitau3 = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_puppitau3");
        m_prodJetsNoLep_puppiSubJetsBdisc = new TTreeReaderValue<std::vector<double>>(m_ttree,"prodJetsNoLep_puppiSubJetsBdisc");
*/
    }

    // LEPTONS
    m_useNeutrinos = (m_config->useNeutrinos());
    if (m_config->useLeptons()){
        m_useLeptons = true;

        m_nMuons       = new TTreeReaderValue<int>(m_ttree,"nMuons");
        m_nMuons_CUT   = new TTreeReaderValue<int>(m_ttree,"nMuons_CUT");
        m_muonsCharge  = new TTreeReaderValue<std::vector<double>>(m_ttree,"muonsCharge");
        m_muonsMtw     = new TTreeReaderValue<std::vector<double>>(m_ttree,"muonsMtw");
        m_muonsRelIso  = new TTreeReaderValue<std::vector<double>>(m_ttree,"muonsRelIso");
        m_muonsMiniIso = new TTreeReaderValue<std::vector<double>>(m_ttree,"muonsMiniIso");
        m_muonspfActivity = new TTreeReaderValue<std::vector<double>>(m_ttree,"muonspfActivity");
        m_muonsFlagMedium  = new TTreeReaderValue<std::vector<int>>(m_ttree,"muonsFlagMedium");
        m_muonsFlagTight = new TTreeReaderValue<std::vector<int>>(m_ttree,"muonsFlagTight");
        //m_muMatchedJetIdx  = new TTreeReaderValue<std::vector<int>>(m_ttree,"muMatchedJetIdx");
        m_muonsLVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"muonsLVec");

        m_nElectrons_CUT = new TTreeReaderValue<int>(m_ttree,"nElectrons_CUT");
        m_nElectrons = new TTreeReaderValue<int>(m_ttree,"nElectrons");
        m_elesCharge = new TTreeReaderValue<std::vector<double>>(m_ttree,"elesCharge");
        m_elesMtw = new TTreeReaderValue<std::vector<double>>(m_ttree,"elesMtw");
        m_elesRelIso = new TTreeReaderValue<std::vector<double>>(m_ttree,"elesRelIso");
        m_elesMiniIso = new TTreeReaderValue<std::vector<double>>(m_ttree,"elesMiniIso");
        m_elespfActivity = new TTreeReaderValue<std::vector<double>>(m_ttree,"elespfActivity");
        m_elesisEB = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"elesisEB");
        m_elesLVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"elesLVec");
        m_elesFlagMedium   = new TTreeReaderValue<std::vector<int>>(m_ttree,"elesFlagMedium");
        m_elesFlagVeto     = new TTreeReaderValue<std::vector<int>>(m_ttree,"elesFlagVeto");

        m_W_emuVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"W_emuVec");
        m_W_emu_pfActivityVec = new TTreeReaderValue<std::vector<double>>(m_ttree,"W_emu_pfActivityVec");
        m_W_tauVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"W_tauVec");
        m_W_tau_emuVec    = new TTreeReaderValue<std::vector<int>>(m_ttree,"W_tau_emuVec");
        m_W_tau_nuVec     = new TTreeReaderValue<std::vector<int>>(m_ttree,"W_tau_nuVec");
        m_W_tau_prongsVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"W_tau_prongsVec");
        m_W_tau_emu_pfActivityVec    = new TTreeReaderValue<std::vector<double>>(m_ttree,"W_tau_emu_pfActivityVec");
        m_W_tau_prongs_pfActivityVec = new TTreeReaderValue<std::vector<double>>(m_ttree,"W_tau_prongs_pfActivityVec");
    }


    // TRUTH
    m_useTruth = (m_config->useTruth());

    if (m_isMC){
/*
        m_genjetsLVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"genjetsLVec");
        m_genHT = new TTreeReaderValue<double>(m_ttree,"genHT");
        m_selGenParticle = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"selGenParticle");
        m_evtWeight = new TTreeReaderValue<double>(m_ttree,"evtWeight");
        m_q = new TTreeReaderValue<double>(m_ttree,"q");
        m_x1 = new TTreeReaderValue<double>(m_ttree,"x1");
        m_x2 = new TTreeReaderValue<double>(m_ttree,"x2");
        m_ScaleWeightsMiniAOD = new TTreeReaderValue<std::vector<double>>(m_ttree,"ScaleWeightsMiniAOD");
*/
        m_selPDGid = new TTreeReaderValue<std::vector<int>>(m_ttree,"selPDGid");
        m_genDecayLVec   = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"genDecayLVec");
        m_genDecayIdxVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayIdxVec");
        m_genDecayStrVec = new TTreeReaderValue<std::vector<std::string>>(m_ttree,"genDecayStrVec");
        m_genDecayPdgIdVec  = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayPdgIdVec");
        m_genDecayMomIdxVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayMomIdxVec");
        m_genDecayMomRefVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayMomRefVec");

        m_stored_weight = new TTreeReaderValue<double>(m_ttree,"stored_weight");

        m_genmet    = new TTreeReaderValue<double>(m_ttree,"genmet");
        m_genmetphi = new TTreeReaderValue<double>(m_ttree,"genmetphi");

        m_tru_npv = new TTreeReaderValue<double>(m_ttree,"tru_npv");
    } // end isMC


    // Truth matching tool
    m_truthMatchingTool = new truthMatching(cmaConfig);
    m_truthMatchingTool->initialize();

    m_ttbarRecoTool    = new ttbarReco(cmaConfig);
    m_deepLearningTool = new deepLearning(cmaConfig);

    // DNN material
    bool useDNN(false);
    if (!m_getDNN && useDNN)  // always false for now
        m_dnn_score = new TTreeReaderValue<double>(m_ttree,"dnn_score");
} // end constructor


Event::~Event() {}



void Event::initialize_eventWeights(){
    /* Create vectors of the systematics that are weights for the nominal events */
    std::map<std::string,unsigned int> mapWeightSystematics = m_config->mapOfWeightVectorSystematics();

    m_listOfWeightSystematics = m_config->listOfWeightSystematics();

    m_weightSystematicsFloats.clear();
    m_weightSystematicsVectorFloats.clear();

    // systematics from the nominal tree that are floats
    for (const auto& nom_syst : m_listOfWeightSystematics){
        if (!m_useLeptons && nom_syst.find("leptonSF")!=std::string::npos)
            continue;
        m_weightSystematicsFloats[nom_syst] = new TTreeReaderValue<float>(m_ttree,nom_syst.c_str());
    }

    // systematics from the nominal tree that are vectors
    for (const auto& syst : mapWeightSystematics)
        m_weightSystematicsVectorFloats[syst.first] = new TTreeReaderValue<std::vector<float>>(m_ttree,syst.first.c_str());

    return;
}


void Event::clear(){
    /* Clear many of the vectors/maps for each event -- SAFETY PRECAUTION */
    m_truth_ljets.clear();
    m_truth_jets.clear();
    m_truth_leptons.clear();
    m_truth_neutrinos.clear();
    m_truth_partons.clear();
    m_truth_tops.clear();
    
    m_ljets.clear();
    m_ljetsPUPPI.clear();
    m_jets.clear();
    m_leptons.clear();
    m_neutrinos.clear();

    m_btag_jets.clear();
    m_btag_jets_default.clear();
    m_weight_btag_default = 1.0;

    m_dnnInputs.clear();

    m_HT = 0;
    m_ST = 0;

    return;
}


void Event::updateEntry(Long64_t entry){
    /* Update the entry -> update all TTree variables */
    cma::DEBUG("EVENT : Update Entry "+std::to_string(entry) );

    m_entry = entry;
    cma::DEBUG("EVENT : Set entry for updating ");

    // make sure the entry exists
    // when looping over truth events, this condition is not always met
    if(isValidRecoEntry())
        m_ttree.SetEntry(m_entry);
    else
        cma::ERROR("EVENT : Invalid Reco entry "+std::to_string(m_entry)+"!");

    return;
}


void Event::execute(Long64_t entry){
    /* Get the values from the event */
    cma::DEBUG("EVENT : Execute event " );

    // Load data from root tree for this event
    updateEntry(entry);

    // Reset many event-level values
    clear();

    // Get the event weights (for cutflows/histograms)
    initialize_weights();
    cma::DEBUG("EVENT : Setup weights ");

    // Truth Information (before other physics objects, to do truth-matching)
    if (m_useTruth && m_isMC){
        initialize_truth();
        cma::DEBUG("EVENT : Setup truth information ");
    }

    // Jets
    if (m_useJets){
        initialize_jets();
        cma::DEBUG("EVENT : Setup small-R jets ");
    }

    // Large-R Jets
    if (m_useLargeRJets){
        initialize_ljets();
        cma::DEBUG("EVENT : Setup large-R jets ");
    }

    // Leptons
    if (m_useLeptons){
        initialize_leptons();
        cma::DEBUG("EVENT : Setup leptons ");
    }

    // Get some kinematics
    initialize_kinematics();
    cma::DEBUG("EVENT : Setup kinematic variables ");

    // Neutrinos
    if (m_useNeutrinos){
        // Need ALL other information from the event to do this
        initialize_neutrinos();
        cma::DEBUG("EVENT : Setup neutrinos ");
    }


    // ------------- //

    // Ttbar Reconstruction
    m_ttbarRecoTool->execute(m_jets,m_ljets);
    m_ttbar = m_ttbarRecoTool->tops();

    // DNN prediction for each Top object
    for (auto& top : m_ttbar)
        deepLearningPrediction(top);

    cma::DEBUG("EVENT : Setup Event ");

    return;
}


void Event::initialize_truth(){
    /* Setup struct of truth information */
    // Truth AK4 jets
    m_truth_jets.clear();

    unsigned int nGenJets = 0;  //(*m_genjetsLVec)->size();
    for (unsigned int i=0; i<nGenJets; i++){
        Jet jet;
        jet.p4     = (*m_genjetsLVec)->at(i);
        jet.index  = i;
        jet.radius = 0.4;
        m_truth_jets.push_back( jet );
    }

    // Truth AK8 jets
    m_truth_ljets.clear(); // not setup

    // Truth kinematics
    //MET,HT


    // Truth Partons
    m_truth_partons.clear();
    unsigned int nPartons( (*m_genDecayLVec)->size() );

    // Collect truth top information into one value
    unsigned int t_idx(0);  // keeping track of tops in m_truth_tops
    m_truth_tops.clear();

    // loop over truth partons
    for (unsigned int i=0; i<nPartons; i++){
        Parton parton;
        parton.p4 = (*m_genDecayLVec)->at(i);

        int pdgId = (*m_genDecayPdgIdVec)->at(i);
        unsigned int abs_pdgId = std::abs(pdgId);

        parton.pdgId = pdgId;

        // simple booleans for type
        parton.isTop = ( abs_pdgId==6 );
        parton.isW   = ( abs_pdgId==24 );
        parton.isLepton = ( abs_pdgId>=11 && abs_pdgId<=16 );
        parton.isQuark  = ( abs_pdgId<7 );

        if (parton.isLepton){
            parton.isTau  = ( abs_pdgId==15 );
            parton.isMuon = ( abs_pdgId==13 );
            parton.isElectron = ( abs_pdgId==11 );
            parton.isNeutrino = ( abs_pdgId==12 || abs_pdgId==14 || abs_pdgId==16 );
        }
        else if (parton.isQuark){
            parton.isLight  = ( abs_pdgId<5 );
            parton.isBottom = ( abs_pdgId==5 );
        }

        parton.index      = i;                              // index in vector of truth_partons
        parton.decayIdx   = (*m_genDecayIdxVec)->at(i);     // index in full truth record of parton
        parton.parent_ref = (*m_genDecayMomRefVec)->at(i);  // index in truth vector of parent
        parton.parent_idx = (*m_genDecayMomIdxVec)->at(i);  // index in full truth record of parent
        parton.top_index  = -1;                             // index in truth_tops vector
        parton.containment = 0;                             // value for determining matching

        // build truth top structs
        // in truth parton record, the top should arrive before its children
        TruthTop top;

        if (parton.isTop){
            top.Wdecays.clear();    // for storing W daughters
            top.daughters.clear();  // for storing non-W/bottom daughters

            top.Top       = parton.index;
            top.isTop     = (pdgId>0);
            top.isAntiTop = (pdgId<0);
            top.isHadronic = false;     // initialize
            top.isLeptonic = false;     // initialize
            parton.top_index = t_idx;

            m_truth_tops.push_back(top);   // store tops now, add information from children in the next iterations over partons
            t_idx++;
        }
        else if (parton.parent_ref>0){
            // check the parent! (ignore parent_ref<=0 -- doesn't exist)
            Parton parent = m_truth_partons.at(parton.parent_ref);

            // check if grandparent exists
            int top_index(-1);            // to refer to (grand)parent top quark in m_truth_tops
            if(parent.parent_ref>0) {
                Parton gparent = m_truth_partons.at(parent.parent_ref);  // parent of parent
                if (gparent.isTop) top_index = gparent.top_index;
            }

            // Parent is Top (W or b)
            if (parent.isTop){
                top = m_truth_tops.at(parent.top_index);
                if (parton.isW) top.W = parton.index;
                else if (parton.isBottom) {
                    top.bottom = parton.index;
                    parton.containment = m_mapContainment.at("BONLY");
                    if (top.isAntiTop) parton.containment*=-1;
                }
                else top.daughters.push_back( parton.index );        // non-W/bottom daughter
                m_truth_tops[parent.top_index] = top;                // update entry
            }
            // Parent is W
            else if (parent.isW && top_index>=0){
                top = m_truth_tops.at(top_index);
                top.Wdecays.push_back(parton.index);
                top.isHadronic = (parton.isQuark);
                top.isLeptonic = (parton.isLepton);

                parton.containment = m_mapContainment.at("QONLY");
                if (top.isAntiTop) parton.containment*=-1;

                m_truth_tops[top_index] = top;      // update entry
            }
        } // end else if

        // store for later access
        m_truth_partons.push_back( parton );
    } // end loop over truth partons

    m_truthMatchingTool->setTruthPartons(m_truth_partons);
    m_truthMatchingTool->setTruthTops(m_truth_tops);

    return;
}


void Event::initialize_jets(){
    /* Setup struct of jets (small-r) and relevant information 
     * b-tagging: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
        cMVAv2L -0.5884	
        cMVAv2M 0.4432	 	 
        cMVAv2T 0.9432
     */
    m_jets.clear();  // don't know a priori the number of jets that pass kinematics
    m_ak4candidates.clear();

    // book-keeping for the b-tagging working points. Set in 'getBtaggedJets()'
    m_btag_jets["L"].clear();
    m_btag_jets["M"].clear();
    m_btag_jets["T"].clear();

    unsigned int j_idx(0); // counting jets that pass kinematic cuts
    for (unsigned int i=0,size=(*m_jetsLVec)->size(); i<size; i++){
        Jet jet;
        jet.p4 = (*m_jetsLVec)->at(i);

        // kinematic cuts
        if (jet.p4.Pt() < 30 || std::abs( jet.p4.Eta() > 2.4 ) ) continue;

        // Other properties
        jet.bdisc = (*m_recoJetsBtag_0)->at(i);
        jet.charge = (*m_recoJetsCharge_0)->at(i);
        jet.true_flavor = (m_isMC) ? (*m_recoJetsFlavor)->at(i) : -1;
        jet.index  = j_idx;
        jet.radius = 0.4;
        jet.containment = 0;   // initialize in case this is data

        // b-tagging
        getBtaggedJets(jet);

        // truth matching
        if (m_useTruth){
            cma::DEBUG("EVENT : Truth match AK4 jets");

            // parton
            m_truthMatchingTool->matchJetToTruthTop(jet);
            if (jet.containment!=0) m_ak4candidates.push_back(jet.index);

            // truth jet
            m_truthMatchingTool->matchJetToTruthJet(jet,m_truth_jets);
        }

        m_jets.push_back(jet);
        j_idx++;
    }

    // store the configured b-tag WP information for fast+convenient access
    m_btag_jets_default = m_btag_jets.at(m_config->jet_btagWkpt());

    return;
}


void Event::initialize_ljets(){
    /* Setup struct of large-R jets and relevant information */
    unsigned int n_ljets  = (*m_ak8JetsLVec)->size();

    m_ak8candidates.clear();       // truth matching candidates
    m_ak8candidatesPUPPI.clear();  // truth matching candidates
    m_ljets.clear();               // AK8 CHS
    m_ljetsPUPPI.clear();          // AK8 PUPPI

    // -- CHS (DeepAK8)
    //    Only 4-vector and DeepAK8 scores available for this collection
    //    Use the 16 outputs of DeepAK8 for NN (rather than engineered features)
    unsigned int j_idx(0);  // counting ljets that pass kinematic cuts
    for (unsigned int i=0; i<n_ljets; i++){
        Ljet ljet;
        ljet.p4 = (*m_ak8JetsLVec)->at(i);

        ljet.isGood = (ljet.p4.Pt()>200. && fabs(ljet.p4.Eta())<2.4 && ljet.p4.M()>10.);  // no softdrop mass available
        // kinematic cuts
        if (!ljet.isGood) continue;

        ljet.index  = j_idx;
        ljet.radius = 0.8;
        ljet.containment = 0;   // set initial value (in case this is data)
        ljet.matchId     = -1;

        ljet.deepAK8.clear();
        for (const auto& x : (*m_ak8JetsDeepAK8)->at(i))
            ljet.deepAK8.push_back(x);

        // Truth-matching to jet
        ljet.truth_partons.clear();
        if (m_useTruth && m_config->isTtbar()) {
            cma::DEBUG("EVENT : Truth match AK8 subjets");  // match subjets (and then the AK8 jet) to truth tops

            // First AK8 subjet
            m_truthMatchingTool->matchJetToTruthTop(ljet);  // match to partons

            cma::DEBUG("EVENT : ++ Ljet containment = "+std::to_string(ljet.containment)+" for truth top "+std::to_string(ljet.matchId));

            // add partons to the AK8 jet
            for (const auto& p : ljet.truth_partons){
                cma::DEBUG("EVENT :    Ljet parton "+std::to_string(m_truth_partons.at(p).pdgId)+" with value "+std::to_string(m_truth_partons.at(p).containment));
                ljet.truth_partons.push_back(p);
            }

            if (ljet.containment!=0){
                m_ak8candidates.push_back( ljet.index );
                cma::DEBUG("EVENT : Ljet that is a candidate with cont. = "+std::to_string(ljet.containment));
            }
        } // end truth matching ljet to partons

        m_ljets.push_back(ljet);
        j_idx++;
    }


    // -- PUPPI
    j_idx   = 0;  // counting ljets that pass kinematic cuts
    n_ljets = (*m_puppiJetsLVec)->size();
    for (unsigned int i=0; i<n_ljets; i++){
        Ljet ljet;
        ljet.p4 = (*m_puppiJetsLVec)->at(i);
        ljet.softDropMass = (*m_puppitau3)->at(i*2+1);  //bug in (*m_puppisoftDropMass)->at(i); values stored in Tau3

        ljet.isGood = (ljet.p4.Pt()>200. && fabs(ljet.p4.Eta())<2.4 && ljet.softDropMass>20.);
        // kinematic cuts
        if (!ljet.isGood) continue;

        // soft drop subjets -- need to check all subjets (?)
        cma::DEBUG("EVENT : Loop over sub-ljets, "+std::to_string((*m_puppiSubJetsLVec)->size()));

        ljet.subjets.clear();
        // have to check all subjets because they may not match indices with the ljet
        unsigned int nsubjets = (*m_puppiSubJetsLVec)->at(i).size();
        for (unsigned int k=0; k<nsubjets; k++){
            Jet subjet;
            subjet.p4 = (*m_puppiSubJetsLVec)->at(i).at(k); //(*m_ak8SubJetsLVec)->at(j);
            subjet.radius = 0.4;

            if ( subjet.p4.DeltaR( ljet.p4 ) > subjet.radius) break; // only consider subjets within R<0.4; continue to next subjet j

            subjet.bdisc = (*m_puppiSubJetsBdisc)->at(i).at(k);    //(*m_ak8SubJetsBdisc)->at(j);
            subjet.index  = k;      // don't keep a vector of all matched subjets; maintain index in m_ak8SubJetsLVec
            subjet.containment = 0;
            ljet.subjets.push_back(subjet);
        }

        // soft drop quality cut
        if (ljet.subjets.size()<2) continue; // require at least two subjets

        // Other properties
        // bug in tau3/softdrop mass; SD mass stored in tau3 vector, too
        cma::DEBUG("EVENT : Other properties, "+std::to_string((*m_puppitau1)->size())+"; "+std::to_string((*m_ak8JetsLVec)->size()));
        ljet.tau1  = (*m_puppitau1)->at(i); //(*m_tau1)->at(i);
        ljet.tau2  = (*m_puppitau2)->at(i); //(*m_tau2)->at(i);
        ljet.tau3  = (*m_puppitau3)->at(i*2); //(*m_tau3)->at(i); bug stores softdropmass here, too
        ljet.tau21 = ljet.tau2 / ljet.tau1;
        ljet.tau32 = ljet.tau3 / ljet.tau2;
        ljet.index  = j_idx;
        ljet.radius = 0.8;
        ljet.containment = 0;   // set initial value (in case this is data)
        ljet.matchId     = -1;

        // Truth-matching to subjets
        ljet.truth_partons.clear();
        if (m_useTruth && m_config->isTtbar()) {
            cma::DEBUG("EVENT : Truth match AK8 subjets");  // match subjets (and then the AK8 jet) to truth tops

            // First AK8 subjet
            Jet sub0 = ljet.subjets.at(0);
            m_truthMatchingTool->matchJetToTruthTop(sub0);  // match to partons

            ljet.containment += sub0.containment;
            ljet.matchId      = sub0.matchId;       // set the truth parton to which the AK8 jet is matched
            cma::DEBUG("EVENT : ++ Subjet 0 containment = "+std::to_string(sub0.containment)+" for truth top "+std::to_string(sub0.matchId));

            // add partons to the AK8 jet
            for (const auto& p : sub0.truth_partons){
                cma::DEBUG("EVENT :    Subjet 0 parton "+std::to_string(m_truth_partons.at(p).pdgId)+" with value "+std::to_string(m_truth_partons.at(p).containment));
                ljet.truth_partons.push_back(p);
            }

            // Second AK8 subjet
            Jet sub1 = ljet.subjets.at(1);
            m_truthMatchingTool->matchJetToTruthTop(sub1);  // match to partons
            int sub_containment = sub1.containment;
            cma::DEBUG("EVENT : ++ Subjet 1 containment = "+std::to_string(sub_containment)+" for truth top "+std::to_string(sub1.matchId));

            // set the 'matchId'
            if (sub0.matchId>=0 && sub1.matchId>=0 && sub0.matchId != sub1.matchId){
                cma::WARNING("EVENT : *** AK8 Jet "+std::to_string(ljet.index)+" truth-matched to multiple top partons! "+std::to_string(sub0.matchId)+" & "+std::to_string(sub1.matchId));
                for (const auto& top : m_truth_tops){
                    std::string name = top.isTop ? "top" : "antitop";
                    Parton tquark    = m_truth_partons.at( top.Top );
                    cma::DEBUG("EVENT :     DeltaR (ljet, truth "+name+") = "+std::to_string(tquark.p4.DeltaR(ljet.p4)));
                    cma::DEBUG("EVENT :     DeltaR (sub0, truth "+name+") = "+std::to_string(tquark.p4.DeltaR(sub0.p4)));
                    cma::DEBUG("EVENT :     DeltaR (sub1, truth "+name+") = "+std::to_string(tquark.p4.DeltaR(sub1.p4)));
                }
                Parton tquark1 = m_truth_partons.at( m_truth_tops.at(0).Top );
                Parton tquark2 = m_truth_partons.at( m_truth_tops.at(1).Top );
                cma::DEBUG("EVENT :     DeltaR (truth top,truth antitop) = "+std::to_string(tquark1.p4.DeltaR(tquark2.p4)));
            }
            else if (sub1.matchId>=0 && sub0.matchId<0){
                ljet.matchId = sub1.matchId;
            }
            // else if (sub0.matchId>=0 && sub1.matchId<0): ljet initialized to sub0.matchId above
            // else (sub1.matchId == sub0.matchId) do nothing

            // add partons to ljet (if they don't already exist!)
            // e.g., one subjet has B and the other BQ, then the ljet should be BQ, not BBQ!
            for (const auto& p : sub1.truth_partons){
                int contValue = m_truth_partons.at(p).containment;
                cma::DEBUG("EVENT :    Subjet 1 parton "+std::to_string(m_truth_partons.at(p).pdgId)+" with value "+std::to_string(contValue));
                if (std::find(ljet.truth_partons.begin(), ljet.truth_partons.end(), p) == ljet.truth_partons.end())
                    ljet.truth_partons.push_back(p);
                else
                    sub_containment -= contValue;
            } 
            cma::DEBUG("EVENT :    Update ljet containment, add "+std::to_string(sub_containment));
            ljet.containment += sub_containment; // adjust the sub0 containment depending on overlaps

            if (ljet.containment!=0){
                m_ak8candidatesPUPPI.push_back( ljet.index );
                cma::DEBUG("EVENT : PUPPI Ljet that is a candidate with cont. = "+std::to_string(ljet.containment));
            }
        } // end truth matching ljet to partons

        m_ljetsPUPPI.push_back(ljet);
        j_idx++;
    }

    return;
}


void Event::initialize_leptons(){
    /* Setup struct of lepton and relevant information */
    m_leptons.clear();

    for (unsigned int i=0,size=(**m_nElectrons); i<size; i++){
        Lepton lep;
        lep.p4 = (*m_elesLVec)->at(i);
        lep.charge = (*m_elesCharge)->at(0);
        lep.isElectron = true;
        lep.isMuon     = false;
        lep.Iso        = (*m_elesMiniIso)->at(i);

        m_leptons.push_back(lep);
    }

    for (unsigned int i=0,size=(**m_nMuons); i<size; i++){
        Lepton lep;
        lep.p4 = (*m_muonsLVec)->at(i);
        lep.charge = (*m_muonsCharge)->at(0);
        lep.isElectron = true;
        lep.isMuon     = false;
        lep.Iso        = (*m_muonsMiniIso)->at(i);

        m_leptons.push_back(lep);
    }

    return;
}


void Event::initialize_neutrinos(){
    /* Build the neutrinos */
    m_neutrinos.clear();
    return;
}


void Event::initialize_kinematics(){
    /* Kinematic variables (HT, ST, MET) */
    m_met_met = *(*m_met);
    m_met_phi = *(*m_metphi);

    m_HT = 0.0;   // total transverse hadronic energy
    m_ST = 0.0;   // total transverse energy

    for (const auto& jet : m_jets)
        m_HT += jet.p4.Pt();

    m_ST += m_HT;

    m_ST = m_met_met;
    if (m_useLeptons){
        for (const auto& lep : m_leptons)
            m_ST += lep.p4.Pt(); 
    }

    return;
}


void Event::getBtaggedJets( Jet& jet ){
    /* Determine the b-tagging */
    jet.isbtagged["L"] = false;
    jet.isbtagged["M"] = false;
    jet.isbtagged["T"] = false;

    if (jet.bdisc > m_cMVAv2L){
        jet.isbtagged["L"] = true;
        m_btag_jets["L"].push_back(jet.index);  // 0 = index of this jet
        if (jet.bdisc > m_cMVAv2M){
            jet.isbtagged["M"] = true;
            m_btag_jets["M"].push_back(jet.index);
            if (jet.bdisc > m_cMVAv2T){
                jet.isbtagged["T"] = true;
                m_btag_jets["T"].push_back(jet.index);
            }
        }
    }

    return;
}


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;

    m_weight_btag.clear();
    if (m_isMC) m_nominal_weight = **m_stored_weight; //**m_evtWeight;

    return;
}

double Event::getSystEventWeight( const std::string &syst, const int weightIndex ){
    /* Calculate the event weight given some systematic
       -- only call for nominal events and systematic weights
       -- for non-nominal tree systematics, use the nominal event weight

       @param syst          Name of systematic (nominal or some weight systematic)
       @param weightIndex   Index of btagging SF; default to -1
    */
    double syst_event_weight(1.0);

    if (syst.compare("nominal")==0){
        // nominal event weight
        syst_event_weight  = m_nominal_weight;
    }
    else{
        // safety to catch something weird -- just return 1.0
        cma::WARNING("EVENT : Passed systematic variation, "+syst+", to Event::getSystEventWeight() ");
        cma::WARNING("EVENT : that is inconsistent with the goldilocks options of ");
        cma::WARNING("EVENT :     nominal, jvt, pileup, leptonSF, and bTagSF. ");
        cma::WARNING("EVENT : Returning a weight of 1.0. ");
        syst_event_weight = 1.0;
    }

    return syst_event_weight;
}


std::vector<int> Event::btag_jets(const std::string &wkpt){
    /* Small-R Jet b-tagging */
    std::string tmp_wkpt(wkpt);
    if(m_btag_jets.find(wkpt) == m_btag_jets.end()){
        cma::WARNING("EVENT : B-tagging working point "+wkpt+" does not exist.");
        cma::WARNING("EVENT : Return vector of b-tagged jets for default working point "+m_config->jet_btagWkpt());
        tmp_wkpt = m_config->jet_btagWkpt();
    }
    return m_btag_jets.at(tmp_wkpt);
}

float Event::met( const std::string& met_name ){
    // MET
    float met_value(0.0);
    if (met_name.compare("met")==0)
        met_value = m_met_met;
    else if (met_name.compare("phi")==0)
        met_value = m_met_phi;
    else{
        cma::WARNING("EVENT : Request for MET variable that is neither 'met' nor 'phi': "+met_name);
        cma::WARNING("EVENT : Returning 0.0");
    }

    return met_value;
}

void Event::deepLearningPrediction(Top& top){
    /* Return map of deep learning values */
    m_deepLearningTool->training(top,m_jets,m_ljets);
    top.dnn = m_deepLearningTool->features();

    if (m_getDNN){
        cma::DEBUG("EVENT : Calculate DNN ");
        m_deepLearningTool->inference(top,m_jets,m_ljets);     // decorate the top struct with DNN values
    }

    return;
}

float Event::weight_btag(const std::string &wkpt){
    std::string tmp_wkpt(wkpt);
    if(m_weight_btag.find(wkpt) == m_weight_btag.end()){
        cma::WARNING("EVENT : B-tagging working point "+wkpt+" does not exist");
        cma::WARNING("EVENT : Return calo-jet b-tag SF for default working point "+m_config->jet_btagWkpt());
        tmp_wkpt = m_config->jet_btagWkpt();
    }
    return m_weight_btag[tmp_wkpt];
}

// Get weight systematics
std::map<std::string,float> Event::weightSystematicsFloats(){
    /* systematics floats */
    std::map<std::string,float> tmp_weightSystematicsFloats;
    for (const auto& wsf : m_weightSystematicsFloats)
        tmp_weightSystematicsFloats[wsf.first] = **wsf.second;

    return tmp_weightSystematicsFloats;
}

std::map<std::string,std::vector<float> > Event::weightSystematicsVectorFloats(){
    /* weight systematics stored as vectors */
    std::map<std::string,std::vector<float> > tmp_weightSystematicsVectorFloats;
    for (const auto& wsf : m_weightSystematicsVectorFloats)
        tmp_weightSystematicsVectorFloats[wsf.first] = **wsf.second;

    return tmp_weightSystematicsVectorFloats;
}

float Event::eventNumber(){ return **m_event;}
float Event::runNumber(){   return **m_run;}
int Event::lumiblock(){     return **m_lumi;}



//
void Event::finalize(){
    /* delete variables */
    cma::DEBUG("EVENT : Finalize() ");
    delete m_BadChargedCandidateFilter;
    delete m_BadPFMuonFilter;
    delete m_EcalDeadCellTriggerPrimitiveFilter;
    delete m_HBHEIsoNoiseFilter;
    delete m_HBHENoiseFilter;
    delete m_PassTrigger;
    delete m_TriggerNames;
    delete m_TriggerPrescales;
    delete m_eeBadScFilter;

    delete m_run;
    delete m_event;
    delete m_lumi;

    if (m_useLeptons){
    cma::DEBUG("EVENT : Finalize -- Clear leptons");
    delete m_elesCharge;
    delete m_elesFlagMedium;
    delete m_elesFlagVeto;
    delete m_elesLVec;
    delete m_elesMiniIso;
    delete m_elesMtw;
    delete m_elesRelIso;
    delete m_elesisEB;
    delete m_elespfActivity;
    //delete m_muMatchedJetIdx;
    delete m_muonsCharge;
    delete m_muonsFlagMedium;
    delete m_muonsFlagTight;
    delete m_muonsLVec;
    delete m_muonsMiniIso;
    delete m_muonsMtw;
    delete m_muonsRelIso;
    delete m_muonspfActivity;
    delete m_n0;
    delete m_nElectrons;
    delete m_nElectrons_CUT;
    delete m_nIsoTrks_CUT;
    delete m_nMuons;
    delete m_nMuons_CUT;
    delete m_W_emuVec;
    delete m_W_emu_pfActivityVec;
    delete m_W_tauVec;
    delete m_W_tau_emuVec;
    delete m_W_tau_emu_pfActivityVec;
    delete m_W_tau_nuVec;
    delete m_W_tau_prongsVec;
    delete m_W_tau_prongs_pfActivityVec;
    }

    cma::DEBUG("EVENT : Finalize -- Clear misc 1");
    delete m_avg_npv;
    delete m_calomet;
    delete m_calometphi;
    delete m_forVetoIsoTrksidx;
    delete m_globalTightHalo2016Filter;
    delete m_goodVerticesFilter;

    delete m_trksForIsoVetoLVec;
    delete m_trksForIsoVeto_charge;
    delete m_trksForIsoVeto_dz;
    delete m_trksForIsoVeto_idx;
    delete m_trksForIsoVeto_iso;
    delete m_trksForIsoVeto_pdgId;
    delete m_trksForIsoVeto_pfActivity;
    delete m_vtxSize;
    delete m_loose_isoTrksLVec;
    delete m_loose_isoTrks_charge;
    delete m_loose_isoTrks_dz;
    delete m_loose_isoTrks_idx;
    delete m_loose_isoTrks_iso;
    delete m_loose_isoTrks_mtw;
    delete m_loose_isoTrks_pdgId;
    delete m_loose_isoTrks_pfActivity;
    delete m_loose_nIsoTrks;

    cma::DEBUG("EVENT : Finalize -- Clear Kinematics");
    delete m_met;
    delete m_nm1;
    delete m_np1;
    delete m_npv;

    if (m_useLargeRJets){
    cma::DEBUG("EVENT : Finalize -- Clear Ljets");
    delete m_puppiSubJetsBdisc;
    delete m_puppiSubJetsLVec;
    delete m_puppisoftDropMass;
    delete m_puppitau1;
    delete m_puppitau2;
    delete m_puppitau3;
/*
    delete m_softDropMass;
    delete m_tau1;
    delete m_tau2;
    delete m_tau3;
    delete m_ak8SubJetsBdisc;
    delete m_ak8SubJetsLVec;
*/
    delete m_ak8JetsLVec;
    delete m_qgAxis2;
    delete m_qgLikelihood;
    delete m_qgMult;
    delete m_qgPtD;
    }

    if (m_useJets){
    cma::DEBUG("EVENT : Finalize -- Clear Jets");
    delete m_jetsLVec;
    delete m_looseJetID;
    delete m_NJetsISR;
    delete m_puppiJetsLVec;
    delete m_recoJetsBtag_0;
    delete m_recoJetsCharge_0;
    delete m_recoJetsFlavor;
    delete m_recoJetsJecScaleRawToFull;
    delete m_recoJetsJecUnc;
    delete m_tightJetID;
    delete m_tightlepvetoJetID;
/*
    delete m_recoJetsBtag_0_LepCleaned;
    delete m_recoJetsCharge_0_LepCleaned;
    delete m_recoJetsJecScaleRawToFull_LepCleaned;
    delete m_jetsLVecLepCleaned;
    delete m_tightJetID_NoLep;
    delete m_tightlepvetoJetID_NoLep;
    delete m_looseJetID_NoLep;
    delete m_prodJetsNoLep_puppiSubJetsBdisc;
    delete m_prodJetsNoLep_puppiSubJetsLVec;
    delete m_prodJetsNoLep_puppisoftDropMass;
    delete m_prodJetsNoLep_puppitau1;
    delete m_prodJetsNoLep_puppitau2;
    delete m_prodJetsNoLep_puppitau3;
    delete m_prodJetsNoLep_tau1;
    delete m_prodJetsNoLep_tau2;
    delete m_prodJetsNoLep_tau3;
    delete m_prodJetsNoLep_puppiJetsLVec;
    delete m_prodJetsNoLep_qgAxis2;
    delete m_prodJetsNoLep_qgLikelihood;
    delete m_prodJetsNoLep_qgMult;
    delete m_prodJetsNoLep_qgPtD;
    delete m_recoJetsJecUncLepCleaned;
    delete m_recoJetschargedEmEnergyFraction;
    delete m_recoJetschargedEmEnergyFractionLepCleaned;
    delete m_recoJetschargedHadronEnergyFraction;
    delete m_recoJetschargedHadronEnergyFractionLepCleaned;
    delete m_recoJetsmuonEnergyFraction;
    delete m_recoJetsmuonEnergyFractionLepCleaned;
    delete m_recoJetsneutralEmEnergyFraction;
    delete m_recoJetsneutralEmEnergyFractionLepCleaned;
*/
    }

    if (m_isMC){
    cma::DEBUG("EVENT : Finalize -- Clear MC");
/*
    delete m_genjetsLVec;
    delete m_genHT;
    delete m_selGenParticle;
    delete m_evtWeight;
    delete m_q;
    delete m_x1;
    delete m_x2;
    delete m_ScaleWeightsMiniAOD;
*/
    delete m_selPDGid;
    delete m_tru_npv;
    delete m_stored_weight;
    delete m_genDecayIdxVec;
    delete m_genDecayLVec;
    delete m_genDecayMomIdxVec;
    delete m_genDecayMomRefVec;
    delete m_genDecayPdgIdVec;
    delete m_genDecayStrVec;
    delete m_genmet;
    delete m_genmetphi;
    } // end isMC

    return;
}

// THE END
