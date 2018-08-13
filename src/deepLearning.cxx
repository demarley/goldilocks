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
  m_lwnn(nullptr),
  m_dnnKey(""),
  m_usePuppi(false){
    m_usePuppi = m_config->usePUPPI();

    if (m_config->DNNinference()){
        // Setup lwtnn
        std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
        lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
        m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
        m_dnnKey = m_config->dnnKey();
    }
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

    if (m_usePuppi){
        m_dnnInputs["AK8_pt"] = ljet.p4.Pt();
        m_dnnInputs["AK8_SDmass"] = ljet.softDropMass;
        m_dnnInputs["AK8_tau1"]   = ljet.tau1;
        m_dnnInputs["AK8_tau2"]   = ljet.tau2;
        m_dnnInputs["AK8_tau3"]   = ljet.tau3;
        m_dnnInputs["AK8_tau21"]  = ljet.tau21;
        m_dnnInputs["AK8_tau32"]  = ljet.tau32;

        m_dnnInputs["AK8_subjet0_bdisc"] = ljet.subjets.at(0).bdisc;
        m_dnnInputs["AK8_subjet0_pTrel"] = ljet.subjets.at(0).p4.Pt() / ljet.p4.Pt();

        m_dnnInputs["AK8_subjet1_bdisc"] = ljet.subjets.at(1).bdisc;
        m_dnnInputs["AK8_subjet1_pTrel"] = ljet.subjets.at(1).p4.Pt() / ljet.p4.Pt();
    }

    else{
        unsigned int i(0);
        for (const auto& x : ljet.deepAK8){
            std::string idx = std::to_string(i);
            m_dnnInputs["AK8_deepAK8_"+idx] = x;
            i++;
        }
    }

    Jet jet = m_jets.at( top.jets.at(0) );
    m_dnnInputs["AK4_deepCSVb"]  = jet.deepCSVb;
    m_dnnInputs["AK4_deepCSVbb"] = jet.deepCSVbb;
    m_dnnInputs["AK4_deepCSVc"]  = jet.deepCSVc;
    m_dnnInputs["AK4_deepCSVcc"] = jet.deepCSVcc;
    m_dnnInputs["AK4_deepCSVl"]  = jet.deepCSVl;

    m_dnnInputs["AK4_deepFlavorb"]  = jet.deepFlavorb;
    m_dnnInputs["AK4_deepFlavorbb"] = jet.deepFlavorbb;
    m_dnnInputs["AK4_deepFlavorc"]  = jet.deepFlavorc;
    m_dnnInputs["AK4_deepFlavorg"]  = jet.deepFlavorg;
    m_dnnInputs["AK4_deepFlavoruds"]  = jet.deepFlavoruds;
    m_dnnInputs["AK4_deepFlavorlepb"] = jet.deepFlavorlepb;

    m_dnnInputs["AK4_charge"] = jet.charge;

    m_dnnInputs["AK8AK4_mass"]   = (ljet.p4 + jet.p4).M();
    m_dnnInputs["AK8AK4_deltaR"] = ljet.p4.DeltaR( jet.p4 );

    m_dnnInputs["weight"] = 1. / (ljet.p4 + jet.p4).Pt();  // 1/ljet.p4.Pt() or something

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
