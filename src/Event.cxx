/*
Created:        24 January 2018
Last Updated:    5 August  2018

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
  m_DNN(-999.),
  m_useTruth(false){
    m_isMC     = m_config->isMC();
    m_treeName = m_ttree.GetTree()->GetName();              // for systematics
    m_DNNtraining  = m_config->DNNtraining();               // load DNN inputs for training
    m_DNNinference = m_config->DNNinference();              // load DNN inputs for inference

    m_mapContainment = m_config->mapOfPartonContainment();  // containment map (ints and strings)
    m_targetMap      = m_config->mapOfTargetValues();       // map of target values for NN training

    // ** LOAD BRANCHES FROM TTREE ** //
    // Event Info -- Filters and Triggers
    m_eeBadScFilter   = new TTreeReaderValue<int>(m_ttree,"eeBadScFilter");
    m_BadPFMuonFilter = new TTreeReaderValue<unsigned int>(m_ttree,"BadPFMuonFilter");
    m_noBadMuonsFilter = new TTreeReaderValue<int>(m_ttree,"noBadMuonsFilter");
    m_badMuonsFilter   = new TTreeReaderValue<int>(m_ttree,"badMuonsFilter");
    m_duplicateMuonsFilter = new TTreeReaderValue<int>(m_ttree,"duplicateMuonsFilter");
    m_HBHENoiseFilter = new TTreeReaderValue<unsigned int>(m_ttree,"HBHENoiseFilter");
    m_goodVerticesFilter = new TTreeReaderValue<int>(m_ttree,"goodVerticesFilter");
    m_HBHEIsoNoiseFilter = new TTreeReaderValue<unsigned int>(m_ttree,"HBHEIsoNoiseFilter");
    m_globalTightHalo2016Filter = new TTreeReaderValue<int>(m_ttree,"globalTightHalo2016Filter");
    m_BadChargedCandidateFilter = new TTreeReaderValue<unsigned int>(m_ttree,"BadChargedCandidateFilter");
    m_EcalDeadCellTriggerPrimitiveFilter = new TTreeReaderValue<int>(m_ttree,"EcalDeadCellTriggerPrimitiveFilter");

    m_PassTrigger  = new TTreeReaderValue<std::vector<int>>(m_ttree,"PassTrigger");
    m_TriggerNames = new TTreeReaderValue<std::vector<std::string>>(m_ttree,"TriggerNames");

    // AK8
    m_ak8puppiJetsLVec    = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"puppiAK8LVec");
    m_ak8puppiSubJetsLVec = new TTreeReaderValue<std::vector<std::vector<TLorentzVector>>>(m_ttree,"puppiAK8SubjetLVec");
    m_ak8puppiTau1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau1");
    m_ak8puppiTau2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau2");
    m_ak8puppiTau3 = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8Tau3");
    m_ak8puppiSoftDropMass = new TTreeReaderValue<std::vector<double>>(m_ttree,"puppiAK8SoftDropMass");
    m_ak8puppiSubJetsBdisc = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetBDisc");
    m_ak8puppiSubJetMult   = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetMult");
    m_ak8puppiSubjetPtD   = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetPtD");
    m_ak8puppiSubjetAxis1 = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetAxis1");
    m_ak8puppiSubjetAxis2 = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"puppiAK8SubjetAxis2");
    m_ak8DeepAK8LVec = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"deepAK8LVec");
    m_ak8DeepAK8top = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8btop");
    m_ak8DeepAK8W   = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8bW");
    m_ak8DeepAK8Z   = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8bZ");
    m_ak8DeepAK8Zbb = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8bZbb");
    m_ak8DeepAK8Hbb = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8bHbb");
    m_ak8DeepAK8H4q = new TTreeReaderValue<std::vector<double>>(m_ttree,"deepAK8bH4q");
    m_ak8DeepAK8    = new TTreeReaderValue<std::vector<std::vector<double>>>(m_ttree,"deepAK8raw");

    // AK4
    m_NJetsISR = new TTreeReaderValue<int>(m_ttree,"NJetsISR");
    m_ak4LVec  = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"jetsLVec");
    m_ak4Flavor = new TTreeReaderValue<std::vector<int>>(m_ttree,"recoJetsFlavor");
    m_ak4Charge = new TTreeReaderValue<std::vector<double>>(m_ttree,"recoJetsCharge_0");
    m_ak4looseJetID = new TTreeReaderValue<unsigned int>(m_ttree,"looseJetID");
    m_ak4tightJetID = new TTreeReaderValue<unsigned int>(m_ttree,"tightJetID");
    m_ak4tightlepvetoJetID = new TTreeReaderValue<unsigned int>(m_ttree,"tightlepvetoJetID");
    m_ak4deepCSV_b  = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepCSVb");
    m_ak4deepCSV_bb = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepCSVbb");
    m_ak4deepCSV_c  = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepCSVc");
    m_ak4deepCSV_cc = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepCSVcc");
    m_ak4deepCSV_l  = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepCSVl");
    m_ak4deepFlavor_b    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorb");
    m_ak4deepFlavor_bb   = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorbb");
    m_ak4deepFlavor_lepb = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorlepb");
    m_ak4deepFlavor_c    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorc");
    m_ak4deepFlavor_uds  = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavoruds");
    m_ak4deepFlavor_g    = new TTreeReaderValue<std::vector<double>>(m_ttree,"DeepFlavorg");
    m_ak4qgLikelihood = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgLikelihood");
    m_ak4qgPtD   = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgPtD");
    m_ak4qgAxis1 = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgAxis1");
    m_ak4qgAxis2 = new TTreeReaderValue<std::vector<double>>(m_ttree,"qgAxis2");
    m_ak4qgMult  = new TTreeReaderValue<std::vector<int>>(m_ttree,"qgMult");

    // TRUTH
    m_useTruth = (m_config->useTruth());
    if (m_isMC){
        m_selPDGid   = new TTreeReaderValue<std::vector<int>>(m_ttree,"selPDGid");
        m_genMatched = new TTreeReaderValue<std::vector<double>>(m_ttree,"genMatched");
        m_genDecayLVec   = new TTreeReaderValue<std::vector<TLorentzVector>>(m_ttree,"genDecayLVec");
        m_genDecayIdxVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayIdxVec");
        m_genDecayPdgIdVec  = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayPdgIdVec");
        m_genDecayMomIdxVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayMomIdxVec");
        m_genDecayMomRefVec = new TTreeReaderValue<std::vector<int>>(m_ttree,"genDecayMomRefVec");

        m_stored_weight = new TTreeReaderValue<double>(m_ttree,"stored_weight");
    } // end isMC

    // Truth matching tool
    m_truthMatchingTool = new truthMatching(cmaConfig);
    m_truthMatchingTool->initialize();

    m_ttbarRecoTool    = new ttbarReco(cmaConfig);
    m_deepLearningTool = new deepLearning(cmaConfig);

    // DNN material
} // end constructor

Event::~Event() {}

void Event::clear(){
    /* Clear many of the vectors/maps for each event -- SAFETY PRECAUTION */
    m_truth_partons.clear();
    m_truth_tops.clear();
    
    m_ljets.clear();
    m_ljetsPUPPI.clear();
    m_jets.clear();

    m_dnnInputs.clear();

    return;
}


void Event::updateEntry(Long64_t entry){
    /* Update the entry -> update all TTree variables */
    cma::DEBUG("EVENT : Update Entry "+std::to_string(entry) );

    m_entry = entry;

    // make sure the entry exists/is valid
    if(isValidRecoEntry())
        m_ttree.SetEntry(m_entry);
    else
        cma::ERROR("EVENT : Invalid Reco entry "+std::to_string(m_entry)+"!");

    cma::DEBUG("EVENT : Set entry for updating ");

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
    initialize_jets();
    cma::DEBUG("EVENT : Setup small-R jets ");

    // Large-R Jets
    initialize_ljets();
    cma::DEBUG("EVENT : Setup large-R jets ");

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
    /* Setup struct of jets (small-r) and relevant information */
    m_jets.clear();  // don't know a priori the number of jets that pass kinematics
    m_ak4candidates.clear();

    unsigned int j_idx(0); // counting jets that pass kinematic cuts
    for (unsigned int i=0,size=(*m_ak4LVec)->size(); i<size; i++){
        Jet jet;
        jet.p4 = (*m_ak4LVec)->at(i);

        // kinematic cuts
        if (jet.p4.Pt() < 30 || std::abs( jet.p4.Eta() > 2.4 ) ) continue;

        // Other properties
        jet.charge = (*m_ak4Charge)->at(i);
        jet.true_flavor = (m_isMC) ? (*m_ak4Flavor)->at(i) : -1;
        jet.index  = j_idx;
        jet.radius = 0.4;
        jet.containment = 0;   // initialize in case this is data

        jet.deepCSVb  = (*m_ak4deepCSV_b)->at(i);
        jet.deepCSVbb = (*m_ak4deepCSV_bb)->at(i);
        jet.deepCSVc  = (*m_ak4deepCSV_c)->at(i);
        jet.deepCSVcc = (*m_ak4deepCSV_cc)->at(i);
        jet.deepCSVl  = (*m_ak4deepCSV_l)->at(i);

        jet.deepFlavorb  = (*m_ak4deepFlavor_b)->at(i);
        jet.deepFlavorbb = (*m_ak4deepFlavor_bb)->at(i);
        jet.deepFlavorc  = (*m_ak4deepFlavor_c)->at(i);
        jet.deepFlavorg  = (*m_ak4deepFlavor_g)->at(i);
        jet.deepFlavoruds  = (*m_ak4deepFlavor_uds)->at(i);
        jet.deepFlavorlepb = (*m_ak4deepFlavor_lepb)->at(i);

        // truth matching
        if (m_useTruth){
            cma::DEBUG("EVENT : Truth match AK4 jets");

            // parton
            m_truthMatchingTool->matchJetToTruthTop(jet);
            if (jet.containment!=0) m_ak4candidates.push_back(jet.index);
        }

        m_jets.push_back(jet);
        j_idx++;
    }

    return;
}


void Event::initialize_ljets(){
    /* Setup struct of large-R jets and relevant information */
    unsigned int n_ljets  = (*m_ak8DeepAK8LVec)->size();

    m_ak8candidates.clear();       // truth matching candidates
    m_ak8candidatesPUPPI.clear();  // truth matching candidates
    m_ljets.clear();               // AK8 CHS
    m_ljetsPUPPI.clear();          // AK8 PUPPI

    cma::DEBUG("EVENT : DeepAK8s");
    // -- CHS (DeepAK8)
    //    Only 4-vector and DeepAK8 scores available for this collection
    //    Use the 16 outputs of DeepAK8 for NN (rather than engineered features)
    unsigned int j_idx(0);  // counting ljets that pass kinematic cuts
    for (unsigned int i=0; i<n_ljets; i++){
        Ljet ljet;
        ljet.p4 = (*m_ak8DeepAK8LVec)->at(i);

        // kinematic cuts
        ljet.isGood = (ljet.p4.Pt()>200. && fabs(ljet.p4.Eta())<2.4 && ljet.p4.M()>20.);  // no softdrop mass available
        if (!ljet.isGood) continue;

        ljet.index  = j_idx;
        ljet.radius = 0.8;
        ljet.containment = 0;   // set initial value (in case this is data)
        ljet.matchId     = -1;

        ljet.deepAK8.clear();
        for (const auto& x : (*m_ak8DeepAK8)->at(i))
            ljet.deepAK8.push_back(x);

        ljet.deepAK8top = (*m_ak8DeepAK8top)->at(i);
        ljet.deepAK8W   = (*m_ak8DeepAK8W)->at(i);
        ljet.deepAK8Z   = (*m_ak8DeepAK8Z)->at(i);
        ljet.deepAK8Zbb = (*m_ak8DeepAK8Zbb)->at(i);
        ljet.deepAK8Hbb = (*m_ak8DeepAK8Hbb)->at(i);
        ljet.deepAK8H4q = (*m_ak8DeepAK8H4q)->at(i);

        ljet.subjets.clear();     // no subjets for DeepAK8 jets
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
    cma::DEBUG("EVENT : PUPPI AK8s");
    j_idx   = 0;  // counting ljets that pass kinematic cuts
    n_ljets = (*m_ak8puppiJetsLVec)->size();
    for (unsigned int i=0; i<n_ljets; i++){
        Ljet ljet;
        ljet.p4 = (*m_ak8puppiJetsLVec)->at(i);

        ljet.tau1 = (*m_ak8puppiTau1)->at(i);
        ljet.tau2 = (*m_ak8puppiTau2)->at(i);
        ljet.tau3 = (*m_ak8puppiTau3)->at(i*2);
        ljet.softDropMass = (*m_ak8puppiTau3)->at(i*2+1);  // bug in (*m_ak8puppisoftDropMass)->at(i); values stored in Tau3

        ljet.isGood = (ljet.p4.Pt()>200. && fabs(ljet.p4.Eta())<2.4 && ljet.softDropMass>20.);
        // kinematic cuts
        if (!ljet.isGood) continue;

        // soft drop subjets -- need to check all subjets (?)
        cma::DEBUG("EVENT : Loop over sub-ljets, "+std::to_string((*m_ak8puppiSubJetsLVec)->size()));

        ljet.subjets.clear();
        unsigned int nsubjets = (*m_ak8puppiSubJetsLVec)->at(i).size();
        for (unsigned int k=0; k<nsubjets; k++){
            Jet subjet;
            subjet.p4 = (*m_ak8puppiSubJetsLVec)->at(i).at(k); //(*m_ak8SubJetsLVec)->at(j);
            subjet.radius = 0.4;
            subjet.bdisc = (*m_ak8puppiSubJetsBdisc)->at(i).at(k);    //(*m_ak8SubJetsBdisc)->at(j);
            subjet.index  = k;      // don't keep a vector of all matched subjets; maintain index in m_ak8SubJetsLVec
            subjet.containment = 0;
            ljet.subjets.push_back(subjet);
        }

        // soft drop quality cut
        if (ljet.subjets.size()<2) continue; // require at least two subjets

        // Other properties
        // bug in tau3/softdrop mass; SD mass stored in tau3 vector, too
        cma::DEBUG("EVENT : Other properties, "+std::to_string((*m_ak8puppiTau1)->size())+"; "+std::to_string((*m_ak8puppiJetsLVec)->size()));
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


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;
    if (m_isMC) m_nominal_weight = **m_stored_weight; //**m_evtWeight;

    return;
}


void Event::deepLearningPrediction(Top& top){
    /* Return map of deep learning values */
    m_deepLearningTool->training(top,m_jets,m_ljets);
    top.dnn = m_deepLearningTool->features();

    if (m_DNNinference){
        cma::DEBUG("EVENT : Calculate DNN ");
        m_deepLearningTool->inference(top,m_jets,m_ljets);     // decorate the top struct with DNN values
    }

    return;
}


// -- clean-up
void Event::finalize(){
    /* Delete variables used to access information from TTree */
    cma::DEBUG("EVENT : Finalize() ");
    delete m_PassTrigger;
    delete m_TriggerNames;

    delete m_BadChargedCandidateFilter;
    delete m_BadPFMuonFilter;
    delete m_EcalDeadCellTriggerPrimitiveFilter;
    delete m_HBHEIsoNoiseFilter;
    delete m_HBHENoiseFilter;
    delete m_eeBadScFilter;
    delete m_globalTightHalo2016Filter;
    delete m_goodVerticesFilter;

    cma::DEBUG("EVENT : Finalize -- Clear Ljets");

    cma::DEBUG("EVENT : Finalize -- Clear Jets");

    if (m_isMC){
      cma::DEBUG("EVENT : Finalize -- Clear MC");
      delete m_stored_weight;

      delete m_selPDGid;
      delete m_genDecayIdxVec;
      delete m_genDecayLVec;
      delete m_genDecayMomIdxVec;
      delete m_genDecayMomRefVec;
      delete m_genDecayPdgIdVec;
    } // end isMC

    return;
}

// THE END
