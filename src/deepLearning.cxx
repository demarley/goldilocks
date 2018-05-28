/*
Created:        19 February 2018
Last Updated:   28 May      2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for performing deep learning tasks
*/
#include "Analysis/goldilocks/interface/deepLearning.h"


deepLearning::deepLearning( configuration& cmaConfig ) :
  m_config(&cmaConfig),
  m_lwnn(nullptr){
    // Setup lwtnn
    std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
    lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
    m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
    m_dnnKey = m_config->dnnKey();
  }

deepLearning::~deepLearning() {
    delete m_lwnn;
}


void deepLearning::training(Top& top, const std::vector<Jet>& jets, const std::vector<Ljet>& ljets){
    /* Prepare inputs for training */
    m_jets  = jets;
    m_ljets = ljets;

    loadFeatures(top);

    return;
}

void deepLearning::inference(Top& top, const std::vector<Jet>& jets, const std::vector<Ljet>& ljets){
    /* Obtain results from LWTNN */
    m_jets  = jets;
    m_ljets = ljets;

    loadFeatures(top);
    m_discriminant = m_lwnn->compute(m_dnnInputs);
    top.dnn = m_discriminant;

    return;
}


void deepLearning::loadFeatures(const Top& top){
    /* Calculate DNN features */
    m_dnnInputs.clear();

    // feature calculations
    m_dnnInputs["target"] = top.target;

    Ljet ljet = m_ljets.at( top.ljet );

/* ORIGINAL (pre-DEEPAK8)
    m_dnnInputs["ljet_SDmass"] = ljet.softDropMass;
    m_dnnInputs["ljet_tau1"]   = ljet.tau1;
    m_dnnInputs["ljet_tau2"]   = ljet.tau2;
    m_dnnInputs["ljet_tau3"]   = ljet.tau3;
    m_dnnInputs["ljet_tau21"]  = ljet.tau21;
    m_dnnInputs["ljet_tau32"]  = ljet.tau32;

    m_dnnInputs["ljet_subjet0_bdisc"] = ljet.subjets.at(0).bdisc;
    m_dnnInputs["ljet_subjet0_pTrel"] = ljet.subjets.at(0).bdisc;

    m_dnnInputs["ljet_subjet1_bdisc"] = ljet.subjets.at(1).bdisc;
    m_dnnInputs["ljet_subjet1_pTrel"] = ljet.subjets.at(1).bdisc;
*/
    unsigned int i(0);
    for (const auto& x : ljet.deepAK8){
        std::string idx = std::to_string(i);
        m_dnnInputs["ljet_deepAK8_"+idx] = x;
        i++;
    }

    Jet jet = m_jets.at( top.jets.at(0) );
    m_dnnInputs["jet_bdisc"]  = jet.bdisc;
    m_dnnInputs["jet_charge"] = jet.charge;

    m_dnnInputs["ljet_jet_m"]      = (ljet.p4 + jet.p4).M();
    m_dnnInputs["ljet_jet_deltaR"] = ljet.p4.DeltaR( jet.p4 );

    m_dnnInputs["weight"] = (ljet.p4 + jet.p4).Pt();  // 1/ljet.p4.Pt() or something

    cma::DEBUG("EVENT : Set DNN input values ");

    return;
}

std::map<std::string,double> deepLearning::features(){
    /* return features */
    return m_dnnInputs;
}

std::map<std::string,double> deepLearning::predictions(){
    /* Return the full map to the user */
    return m_discriminant;
}

double deepLearning::prediction(){
    /* Return the score for the default key */
    return m_discriminant.at(m_dnnKey);
}

double deepLearning::prediction(const std::string& key){
    /* Just return the prediction (after execute!) */
    return m_discriminant.at(key);
}

// THE END //
