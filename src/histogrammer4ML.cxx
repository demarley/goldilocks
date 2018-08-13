/*
Created:        29 August 2017
Last Updated:   28 May    2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make histograms for systematic uncertainties (& nominal) 
to go into plots || TRexFitter
*/
#include "Analysis/goldilocks/interface/histogrammer4ML.h"


histogrammer4ML::histogrammer4ML( configuration& cmaConfig, std::string name ) :
  histogrammer::histogrammer(cmaConfig,name),
  m_config(&cmaConfig),
  m_name(name),
  m_usePUPPI(false){
    m_usePUPPI = m_config->usePUPPI();
}

histogrammer4ML::~histogrammer4ML() {}


/**** INITIALIZE HISTOGRAMS ****/


void histogrammer4ML::initialize( TFile& outputFile ){
    /* Setup some values and book histograms */
    outputFile.cd();

    bookHists();

    return;
}

void histogrammer4ML::bookHists(){
    /* 
      Book histograms -- modify/inherit this function for analysis-specific hists 

      @param name   This is the string used to identify histograms for different systematics/event weights

      0 :: NONE = QCD (background)
      1 :: QB-Q = Signal AK8(QB) + AK4(Q)
      2 :: QQ-B = Signal AK8(W)  + AK4(B)
      3 :: QQ/QB-X = Signal AK8(QB/W) + AK4(other)
    */
    cma::DEBUG("HISTOGRAMMER : Init. histograms: "+m_name);

    // features for DNN
    for (const auto& target : m_targets){
        // AK4
        histogrammer::init_hist("AK4_deepCSVb_"+target+"_"+m_name,  200, -1,1);
        histogrammer::init_hist("AK4_deepCSVbb_"+target+"_"+m_name, 200, -1,1);
        histogrammer::init_hist("AK4_deepCSVc_"+target+"_"+m_name,  200, -1,1);
        histogrammer::init_hist("AK4_deepCSVcc_"+target+"_"+m_name, 200, -1,1);
        histogrammer::init_hist("AK4_deepCSVl_"+target+"_"+m_name,  200, -1,1);
        histogrammer::init_hist("AK4_deepFlavorb_"+target+"_"+m_name,   200, -1,1);
        histogrammer::init_hist("AK4_deepFlavorbb_"+target+"_"+m_name,  200, -1,1);
        histogrammer::init_hist("AK4_deepFlavorc_"+target+"_"+m_name,   200, -1,1);
        histogrammer::init_hist("AK4_deepFlavoruds_"+target+"_"+m_name, 200, -1,1);
        histogrammer::init_hist("AK4_deepFlavorg_"+target+"_"+m_name,   200, -1,1);
        histogrammer::init_hist("AK4_deepFlavorlepb_"+target+"_"+m_name,200, -1,1);

        histogrammer::init_hist("AK4_charge_"+target+"_"+m_name, 500, -5,5);

        // AK8 + AK4 system
        histogrammer::init_hist("AK8AK4_deltaR_"+target+"_"+m_name, 50,0,5);
        histogrammer::init_hist("AK8AK4_mass_"+target+"_"+m_name,  500,0,1000);
        histogrammer::init_hist("AK8AK4_mass_deltaR_"+target+"_"+m_name,    500,0,1000,50,0,5); // (AK8+AK4).M() vs DeltaR(AK8,AK4)

        if (m_usePUPPI){
            // High-level (PUPPI)
            histogrammer::init_hist("AK8_pt-"+target+"_"+m_name,     2000,  0.0, 2000.0);
            histogrammer::init_hist("AK8_SDmass-"+target+"_"+m_name,  500,  0.0,  500.0);
            histogrammer::init_hist("AK8_tau1-"+target+"_"+m_name,    200,  0.0,    2.0);
            histogrammer::init_hist("AK8_tau2-"+target+"_"+m_name,    200,  0.0,    2.0);
            histogrammer::init_hist("AK8_tau3-"+target+"_"+m_name,    200,  0.0,    2.0);
            histogrammer::init_hist("AK8_tau21-"+target+"_"+m_name,   100,  0.0,    1.0);
            histogrammer::init_hist("AK8_tau32-"+target+"_"+m_name,   100,  0.0,    1.0);
            histogrammer::init_hist("AK8_subjet0_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
            histogrammer::init_hist("AK8_subjet1_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
            histogrammer::init_hist("AK8_subjet0_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);
            histogrammer::init_hist("AK8_subjet1_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);

            // AK8 substructure vs DeltaR(AK8,AK4)
            histogrammer::init_hist("AK8SDmass_AK8AK4deltaR_"+target+"_"+m_name,500,0,1000,50,0,5); // AK8(SD mass) vs DeltaR(AK8,AK4)
            histogrammer::init_hist("AK8tau21_AK8AK4deltaR_"+target+"_"+m_name,  100,0,1,50,0,5);
            histogrammer::init_hist("AK8tau32_AK8AK4deltaR_"+target+"_"+m_name,  100,0,1,50,0,5);
        }
        else{
            // AK8 (DeepAK8)
            for (unsigned int i=0;i<16;i++)
                histogrammer::init_hist("AK8_deepAK8-"+std::to_string(i)+"_"+target+"_"+m_name,  100, 0,1);
        }
    }

    return;
}


/**** FILL HISTOGRAMS ****/
void histogrammer4ML::fill( const std::map<std::string,double> features, double weight ){
    /* Fill histograms -- 
       Fill information from single top object (inputs to deep learning)
    */
    std::string target = std::to_string( int(features.at("target")) );

    cma::DEBUG("HISTOGRAMMER : Fill histograms: "+m_name+"; target = "+target);

    // Features
    histogrammer::fill("AK4_deepCSVb_"+target+"_"+m_name,  	features.at("AK4_deepCSVb"),  weight);
    histogrammer::fill("AK4_deepCSVbb_"+target+"_"+m_name, 	features.at("AK4_deepCSVbb"), weight);
    histogrammer::fill("AK4_deepCSVc_"+target+"_"+m_name,  	features.at("AK4_deepCSVc"),  weight);
    histogrammer::fill("AK4_deepCSVcc_"+target+"_"+m_name, 	features.at("AK4_deepCSVcc"), weight);
    histogrammer::fill("AK4_deepCSVl_"+target+"_"+m_name,  	features.at("AK4_deepCSVl"),  weight);
    histogrammer::fill("AK4_deepFlavorb_"+target+"_"+m_name,    features.at("AK4_deepFlavorb"),    weight);
    histogrammer::fill("AK4_deepFlavorbb_"+target+"_"+m_name,   features.at("AK4_deepFlavorbb"),   weight);
    histogrammer::fill("AK4_deepFlavorc_"+target+"_"+m_name,    features.at("AK4_deepFlavorc"),    weight);
    histogrammer::fill("AK4_deepFlavoruds_"+target+"_"+m_name,  features.at("AK4_deepFlavoruds"),  weight);
    histogrammer::fill("AK4_deepFlavorg_"+target+"_"+m_name,    features.at("AK4_deepFlavorg"),    weight);
    histogrammer::fill("AK4_deepFlavorlepb_"+target+"_"+m_name, features.at("AK4_deepFlavorlepb"), weight);
    histogrammer::fill("AK4_charge_"+target+"_"+m_name, features.at("AK4_charge"), weight);

    if (m_usePUPPI){
        histogrammer::fill("AK8SDmass_AK8AK4deltaR_"+target+"_"+m_name, features.at("AK8_SDmass"), features.at("AK8AK4_deltaR"), weight );
        histogrammer::fill("AK8tau21_AK8AK4deltaR_"+target+"_"+m_name,  features.at("AK8_tau21"),  features.at("AK8AK4_deltaR"), weight );
        histogrammer::fill("AK8tau32_AK8AK4deltaR_"+target+"_"+m_name,  features.at("AK8_tau32"),  features.at("AK8AK4_deltaR"), weight );

        histogrammer::fill("AK8_SDmass-"+target+"_"+m_name, features.at("AK8_SDmass"), weight);
        histogrammer::fill("AK8_tau1-"+target+"_"+m_name,   features.at("AK8_tau1"),  weight);
        histogrammer::fill("AK8_tau2-"+target+"_"+m_name,   features.at("AK8_tau2"),  weight);
        histogrammer::fill("AK8_tau3-"+target+"_"+m_name,   features.at("AK8_tau3"),  weight);
        histogrammer::fill("AK8_tau21-"+target+"_"+m_name,  features.at("AK8_tau21"), weight);
        histogrammer::fill("AK8_tau32-"+target+"_"+m_name,  features.at("AK8_tau32"), weight);

        histogrammer::fill("AK8_subjet0_bdisc-"+target+"_"+m_name, features.at("AK8_subjet0_bdisc"), weight);
        histogrammer::fill("AK8_subjet0_pTrel-"+target+"_"+m_name, features.at("AK8_subjet0_pTrel"),  weight);
        histogrammer::fill("AK8_subjet1_bdisc-"+target+"_"+m_name, features.at("AK8_subjet1_bdisc"), weight);
        histogrammer::fill("AK8_subjet1_pTrel-"+target+"_"+m_name, features.at("AK8_subjet1_pTrel"), weight);
    }
    else{
        cma::DEBUG("HISTOGRAMMER : ljets deepAK8 ");
        for (unsigned int i=0;i<16;i++){
            std::string idx = std::to_string(i);
            histogrammer::fill("AK8_deepAK8-"+idx+"_"+target+"_"+m_name, features.at("AK8_deepAK8_"+idx), weight);
        }
    }

    cma::DEBUG("HISTOGRAMMER : AK8s-AK4s ");
    histogrammer::fill("AK8AK4_deltaR_"+target+"_"+m_name, features.at("AK8AK4_deltaR"), weight);
    histogrammer::fill("AK8AK4_mass_"+target+"_"+m_name,   features.at("AK8AK4_mass"),   weight);
    histogrammer::fill("AK8AK4_mass_deltaR_"+target+"_"+m_name, features.at("AK8AK4_mass"), features.at("AK8AK4_deltaR"), weight );

    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}

// THE END
