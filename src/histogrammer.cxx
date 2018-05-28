/*
Created:        29 August 2017
Last Updated:   29 May    2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make histograms!
*/
#include "Analysis/goldilocks/interface/histogrammer.h"


histogrammer::histogrammer( configuration& cmaConfig, std::string name ) :
  m_config(&cmaConfig),
  m_name(name),
  m_isMC(false),
  m_doSystWeights(false),
  m_useLargeRJets(false),
  m_useJets(false),
  m_useLeptons(false),
  m_putOverflowInLastBin(true),
  m_putUnderflowInFirstBin(true){
    m_map_histograms1D.clear();
    m_map_histograms2D.clear();
    m_map_histograms3D.clear();

    m_isMC = m_config->isMC();

    m_useLargeRJets = m_config->useLargeRJets();
    m_useJets       = m_config->useJets();
    m_useLeptons    = m_config->useLeptons();

    if (m_name.length()>0  && m_name.substr(m_name.length()-1,1).compare("_")!=0)
        m_name = m_name+"_"; // add '_' to end of string, if needed
  }

histogrammer::~histogrammer() {}


/**** INITIALIZE HISTOGRAMS ****/

// -- 1D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBins, const double x_min, const double x_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms1D["h_"+name] = new TH1D(("h_"+name).c_str(), ("h_"+name).c_str(),nBins,x_min,x_max);
    m_map_histograms1D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBins, const double *xbins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms1D["h_"+name] = new TH1D(("h_"+name).c_str(), ("h_"+name).c_str(),nBins,xbins);
    m_map_histograms1D["h_"+name]->Sumw2();

    return;
}
// -- 2D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms2D["h_"+name] = new TH2D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                            nBinsX,x_min,x_max,nBinsY,y_min,y_max);
    m_map_histograms2D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms2D["h_"+name] = new TH2D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                           nBinsX,xbins,nBinsY,ybins);
    m_map_histograms2D["h_"+name]->Sumw2();

    return;
}
// -- 3D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max,
                              const unsigned int nBinsZ, const double z_min, const double z_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms3D["h_"+name] = new TH3D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                            nBinsX,x_min,x_max,nBinsY,y_min,y_max,nBinsZ,z_min,z_max);
    m_map_histograms3D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins,
                              const unsigned int nBinsZ, const double *zbins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms3D["h_"+name] = new TH3D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                           nBinsX,xbins,nBinsY,ybins,nBinsZ,zbins);
    m_map_histograms3D["h_"+name]->Sumw2();

    return;
}


void histogrammer::initialize( TFile& outputFile, bool doSystWeights ){
    /* Setup some values and book histograms */
    m_doSystWeights = doSystWeights;
    outputFile.cd();

    m_mapContainment    = m_config->mapOfPartonContainment();
    m_mapContainmentRev = m_config->mapOfPartonContainmentRev();

    // loop over treenames (typically systematic uncertainties)
    for (auto& treename : m_config->treeNames() ){
        std::size_t found = treename.find("/");   // protection against directory structure
        if (found!=std::string::npos){
            treename = treename.substr(found+1);
        }

        bookHists( m_name+treename );
    }

    // weight systematics
    if (m_isMC && m_doSystWeights){
        // In ATLAS these were only necessary for the nominal tree.
        // Do not need to make them for every systematic variation!
        for (const auto& syst : m_config->listOfWeightSystematics()){
            bookHists( m_name+syst );
        } // end weight systematics

        // vector weight systematics
        for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
            for (unsigned int el=0;el<syst.second;++el){
                std::string weightIndex = std::to_string(el);
                bookHists( m_name+weightIndex+"_"+syst.first );
            } // end components of vector
        } // end vector weight systematics
    } // end if MC and save weight systematics

    return;
}

void histogrammer::bookHists( std::string name ){
    /* 
      Book histograms -- modify/inherit this function for analysis-specific hists 

      @param name   This is the string used to identify histograms for different systematics/event weights
    */
    m_names.resize(0); // append names to this to keep track of later

    cma::DEBUG("HISTOGRAMMER : Init. histograms: "+name);

    if (m_useLargeRJets){
        for (const auto& c : m_containments){
            std::string cname = c+"_"+name;
            init_hist("ljet_pt_"+cname,     2000,  0.0, 2000.0);
            init_hist("ljet_eta_"+cname,      50, -2.5,    2.5);
            init_hist("ljet_phi_"+cname,      64, -3.2,    3.2);
            init_hist("ljet_SDmass_"+cname,  500,  0.0,  500.0);
            init_hist("ljet_tau1_"+cname,    200,  0.0,    2.0);
            init_hist("ljet_tau2_"+cname,    200,  0.0,    2.0);
            init_hist("ljet_tau3_"+cname,    200,  0.0,    2.0);
            init_hist("ljet_tau21_"+cname,   100,  0.0,    1.0);
            init_hist("ljet_tau32_"+cname,   100,  0.0,    1.0);
            init_hist("ljet_subjet0_btag_"+cname, 100, 0.0, 1.0);
            init_hist("ljet_subjet1_btag_"+cname, 100, 0.0, 1.0);

            init_hist("ljet_pt_eta_"+cname,    200,  0.0, 2000.0,  50, -2.5, 2.5);  // pt vs eta (pt=x-axis)
            init_hist("ljet_pt_SDmass_"+cname, 200,  0.0, 2000.0,  50,  0, 500);    // pt vs SDmass (pt=x-axis)

            // plots with both AK4 and AK8 jets
            if (m_useJets){
                init_hist("deltaR_ljet_jet_"+cname, 500, 0.0, 5.0); // nearest AK4+AK8 (AK4 outside AK8)
                init_hist("ljet_jet_m_"+cname,  1000, 0.0, 1000.0); // nearest AK4+AK8 (AK4 outside AK8)
                init_hist("ljet_jet_pt_"+cname, 2000, 0.0, 2000.0); // nearest AK4+AK8 (AK4 outside AK8)
            }
        }
    }

    if (m_useLargeRJets && m_useJets){
        // AK8+AK4 system of interest
        std::vector<std::string> topAntiTop = {"top","antitop"};
        for (const auto& x : topAntiTop){
            // AK8 mass vs DeltaR(AK8,AK4)
            init_hist(x+"_MassvDR_ljet-W_jet-BONLY_"+name,  1000,0,1000,50,0,5);
            init_hist(x+"_MassvDR_ljet-BQ_jet-QONLY_"+name, 1000,0,1000,50,0,5);
            // AK8 Tau21 vs DeltaR(AK8,AK4)
            init_hist(x+"_Tau21vDR_ljet-W_jet-BONLY_"+name,  100,0,1,50,0,5);
            init_hist(x+"_Tau21vDR_ljet-BQ_jet-QONLY_"+name, 100,0,1,50,0,5);
            // AK8 Tau32 vs DeltaR(AK8,AK4)
            init_hist(x+"_Tau32vDR_ljet-W_jet-BONLY_"+name,  100,0,1,50,0,5);
            init_hist(x+"_Tau32vDR_ljet-BQ_jet-QONLY_"+name, 100,0,1,50,0,5);
            // DeltaR(AK8,AK4)
            init_hist(x+"_deltaR_ljet-W_jet-BONLY_"+name,  500, 0.0, 5.0);
            init_hist(x+"_mass_ljet-W_jet-BONLY_"+name,   1000, 0.0,1000.0);
            init_hist(x+"_deltaR_ljet-BQ_jet-QONLY_"+name, 500, 0.0, 5.0);
            init_hist(x+"_mass_ljet-BQ_jet-QONLY_"+name,  1000, 0.0,1000.0);
        }
    }

    if (m_useJets){
        init_hist("n_jets_"+name,   31, -0.5,  30.5);
        init_hist("n_btags_"+name,  11, -0.5,  10.5);

        init_hist("jet_pt_"+name,  500, 0.0,  500);
        init_hist("jet_eta_"+name,  50, -2.5, 2.5);
        init_hist("jet_phi_"+name,  64, -3.2, 3.2);
        init_hist("jet_bdisc_"+name, 200, -1,1);
    }


    if (m_useLeptons){
        init_hist("lep_pt_"+name,  500, 0.0, 2000);
        init_hist("lep_eta_"+name,  50, -2.5, 2.5);
        init_hist("lep_phi_"+name,  64, -3.2, 3.2);
    }

    // kinematics
    init_hist("met_met_"+name, 500,  0.0,  500);
    init_hist("met_phi_"+name, 6.4, -3.2,  3.2);
    init_hist("ht_"+name,     5000,  0.0, 5000);

    // DNN
    init_hist("dnn_"+name,  100, 0.0,   1.);

    return;
}




/**** FILL HISTOGRAMS ****/

void histogrammer::fill( const std::string& name, const double& value, const double& weight ){
    /* TH1D */
    TH1D* this_hist = m_map_histograms1D.at("h_"+name);

    this_hist->Fill(value,weight);

    return;
}

void histogrammer::fill( const std::string& name, 
                         const double& xvalue, const double& yvalue, const double& weight ){
    /* TH2D */
    TH2D* this_hist = m_map_histograms2D.at("h_"+name);

    this_hist->Fill(xvalue,yvalue,weight);

    return;
}

void histogrammer::fill( const std::string& name, 
                         const double& xvalue, const double& yvalue, const double& zvalue, const double& weight ){
    /* TH3D */
    TH3D* this_hist = m_map_histograms3D.at("h_"+name);

    this_hist->Fill(xvalue,yvalue,zvalue,weight);

    return;
}


void histogrammer::fill( Event& event ){
    /* Fill histograms -- fill histograms based on treename or systematic weights ("nominal" but different weight)
       This is the function to modify / inherit for analysis-specific purposes
    */
    std::string treeName = event.treeName();
    double event_weight  = event.nominal_weight();

    fill( m_name+treeName, event, event_weight );

    // if there are systematics stored as weights (e.g., b-tagging, pileup, etc.)
    // the following calls the fill() function with different event weights
    // to make histograms
    // In ATLAS, these weights only existed in the 'nominal' tree
    bool isNominal = m_config->isNominalTree( treeName );
    if (m_isMC && isNominal && m_doSystWeights){
        // weight systematics
        event_weight = 1.0;
        for (const auto& syst : m_config->listOfWeightSystematics()){
            event_weight = event.getSystEventWeight( syst );
            fill( m_name+syst, event, event_weight );
        } // end weight systematics

        // vector weight systematics
        event_weight = 1.0;
        for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
            for (unsigned int el=0;el<syst.second;++el){
                event_weight = event.getSystEventWeight( syst.first, el );
                std::string weightIndex = std::to_string(el);

                fill( m_name+weightIndex+"_"+syst.first, event, event_weight );
            } // end components of vector
        } // end vector weight systematics
    } // end if nominal and doSystWeights

    return;
}


void histogrammer::fill( const std::string& name, Event& event, double event_weight){
    /* Fill histograms -- just use information from the event and fill histogram
       This is the function to modify / inherit for analysis-specific purposes
    */
    cma::DEBUG("HISTOGRAMMER : Fill histograms: "+name+".");

    // Load some objects from Event
    std::vector<Jet> jets   = event.jets();
    std::vector<Ljet> ljets = event.ljets();
    std::vector<Lepton> leptons = event.leptons();
    std::vector<Top> tops = event.ttbar();  // reconstructed ttbar system

    cma::DEBUG("HISTOGRAMMER : event weight = "+std::to_string(event_weight) );
    int FULL  = m_mapContainment.at("FULL");
    int QB    = m_mapContainment.at("BQ");
    int W     = m_mapContainment.at("W");
    int QONLY = m_mapContainment.at("QONLY");
    int BONLY = m_mapContainment.at("BONLY");

    if (m_useLargeRJets){
        cma::DEBUG("HISTOGRAMMER : Fill large-R jets");
        for (const auto& ljet : ljets){
            if (std::abs(ljet.containment)>FULL){
                cma::WARNING("HISTOGRAMMER : event "+std::to_string(event.entry())+" - AK8 jet "+std::to_string(ljet.index)+"with odd containment: "+std::to_string(ljet.containment));
                continue;
            }
            std::string cname = m_mapContainmentRev[ std::abs(ljet.containment) ]+"_"+name;

            fill("ljet_pt_"+cname,  ljet.p4.Pt(),  event_weight);
            fill("ljet_eta_"+cname, ljet.p4.Eta(), event_weight);
            fill("ljet_phi_"+cname, ljet.p4.Phi(), event_weight);
            fill("ljet_SDmass_"+cname, ljet.softDropMass, event_weight);
            fill("ljet_tau1_"+cname,  ljet.tau1, event_weight);
            fill("ljet_tau2_"+cname,  ljet.tau2, event_weight);
            fill("ljet_tau3_"+cname,  ljet.tau3, event_weight);
            fill("ljet_tau21_"+cname, ljet.tau21, event_weight);
            fill("ljet_tau32_"+cname, ljet.tau32, event_weight);

            fill("ljet_pt_eta_"+cname,   ljet.p4.Pt(), ljet.p4.Eta(), event_weight);      // pt vs eta (pt=x-axis)
            fill("ljet_pt_SDmass_"+cname,ljet.p4.Pt(), ljet.softDropMass, event_weight);  // pt vs SDmass (pt=x-axis)

            fill("ljet_subjet0_btag_"+cname, ljet.subjets.at(0).bdisc, event_weight);
            fill("ljet_subjet1_btag_"+cname, ljet.subjets.at(1).bdisc, event_weight);

            // plots with both AK4 and AK8 jets
            if (m_useJets){
                double minDR(100.0); // initial value to check deltaR 
                int jet_idx(-1);
                unsigned int ii(0);

                for (const auto& jet : jets){
                    double DR = jet.p4.DeltaR(	ljet.p4	);
                    if ( DR > 0.8 && DR < minDR ){
                        minDR   = DR;
                        jet_idx = ii;
                    }
                    ii++;
                } // end loop over AK4

                if (jet_idx>=0){
                    TLorentzVector AK4_AK8 = (ljet.p4 + jets.at(jet_idx).p4);

                    fill("deltaR_ljet_jet_"+cname, minDR,    event_weight); // nearest AK4 (outside AK8)
                    fill("ljet_jet_m_"+cname,  AK4_AK8.M(),  event_weight);
                    fill("ljet_jet_pt_"+cname, AK4_AK8.Pt(), event_weight);
                }
            } // end if use AK4
        } // end loop over AK8
    } // end if use AK8


    if (m_useJets){
        cma::DEBUG("HISTOGRAMMER : Fill small-R jets");
        fill("n_btags_"+name, event.btag_jets().size(), event_weight );
        fill("n_jets_"+name, jets.size(), event_weight );

        for (const auto& jet : jets){
            fill("jet_pt_"+name,  jet.p4.Pt(),   event_weight);
            fill("jet_eta_"+name, jet.p4.Eta(),  event_weight);
            fill("jet_phi_"+name, jet.p4.Phi(),  event_weight);
            fill("jet_bdisc_"+name, jet.bdisc, event_weight);
        }
    }

    // make plots comparing AK8 with AK4 jets of certain containments
    if (m_useLargeRJets && m_useJets){

        cma::DEBUG("HISTOGRAMMER : Looping over reconstruced tops "+std::to_string(tops.size()));
        for (const auto& top : tops){
            cma::DEBUG("HISTOGRAMMER :  > Fatjet index = "+std::to_string(top.ljet));
            if (top.ljet<0) continue; // check AK8 exists in top

            Ljet ljet = ljets.at(top.ljet);
            cma::DEBUG("HISTOGRAMMER :  > Fatjet containment = "+std::to_string(ljet.containment));
            if (ljet.containment!=W && ljet.containment!=QB){
                cma::DEBUG("HISTOGRAMMER :  > Containment is not W or QB");
                continue;
            }

            Jet jet = jets.at(top.jets.at(0));  // get the AK4 associated with this

            std::string pid    = (top.isTop) ? "top" : "antitop";
            std::string ljetCt = m_mapContainmentRev[ljet.containment];
            std::string jetCt  = m_mapContainmentRev[jet.containment];

            cma::DEBUG("HISTOGRAMMER :  > Fill pid "+pid+"; ljet ct "+ljetCt+", jet ct "+jetCt);
            fill(pid+"_deltaR_ljet-"+ljetCt+"_jet-"+jetCt+"_"+name, jet.p4.DeltaR( ljet.p4 ), event_weight );
            fill(pid+"_mass_ljet-"+ljetCt+"_jet-"+jetCt+"_"+name,   (jet.p4+ljet.p4).M(), event_weight );
            fill(pid+"_MassvDR_ljet-"+ljetCt+"_jet-"+jetCt+"_"+name, (jet.p4+ljet.p4).M(), jet.p4.DeltaR( ljet.p4 ), event_weight );
            fill(pid+"_Tau21vDR_ljet-"+ljetCt+"_jet-"+jetCt+"_"+name, ljet.tau21, jet.p4.DeltaR( ljet.p4 ), event_weight );
            fill(pid+"_Tau32vDR_ljet-"+ljetCt+"_jet-"+jetCt+"_"+name, ljet.tau32, jet.p4.DeltaR( ljet.p4 ), event_weight );
        } // end loop over ttbar system
    } // end if useJets and useLargeRJets


    if (m_useLeptons){
        cma::DEBUG("HISTOGRAMMER : Fill leptons");
        for (const auto& lep : leptons){
            fill("lep_pt_"+name,  lep.p4.Pt(),  event_weight);
            fill("lep_eta_"+name, lep.p4.Eta(), event_weight);
            fill("lep_phi_"+name, lep.p4.Phi(), event_weight);
        }
    }


    // kinematics
    cma::DEBUG("HISTOGRAMMER : Fill kinematics");
    fill("met_met_"+name, event.met("met"), event_weight);
    fill("met_phi_"+name, event.met("phi"), event_weight);
    fill("ht_"+name,      event.HT(),       event_weight);

    // DNN
    cma::DEBUG("HISTOGRAMMER : Fill DNN");
    fill("dnn_"+name, event.DNN(), event_weight); // N/A

    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}





/**** OVER/UNDERFLOW ****/

void histogrammer::overUnderFlow(){
    /* Call overflow and underflow functions at once */
    overFlow();
    underFlow();
    return;
}


void histogrammer::overFlow() {
    /* Add overflow to last bin */
    if (!m_putOverflowInLastBin){
        cma::INFO("HISTOGRAMMER : Not putting overflow in last bin(s)");
        return;
    }
    else{
        // loop over 1D histograms
        for (const auto& hist : m_map_histograms1D){
            unsigned int numBins = hist.second->GetNbinsX();
            double overflowContent = hist.second->GetBinContent(numBins + 1);

            hist.second->AddBinContent(numBins,overflowContent); // add overflow to last bin
            hist.second->SetBinContent(numBins+1, 0);            // set overflow to 0
        }
        // loop over 2D histograms
        for (const auto& hist : m_map_histograms2D){
            unsigned int numBinsX = hist.second->GetXaxis()->GetNbins();
            unsigned int numBinsY = hist.second->GetYaxis()->GetNbins();

            // x-axis :: overflow in y
            for (unsigned int xx=1;xx<numBinsX+1;++xx){
                double overflowContent = hist.second->GetBinContent(xx,numBinsY+1);

                int lastBin     = hist.second->GetBin(xx,numBinsY);
                int overflowBin = hist.second->GetBin(xx,numBinsY+1);
                hist.second->AddBinContent(lastBin,overflowContent); // add overflow to last bin
                hist.second->SetBinContent(overflowBin,0);           // set overflow to 0
            }
            // y-axis :: overflow in x
            for (unsigned int yy=1;yy<numBinsY;++yy){
                double overflowContent = hist.second->GetBinContent(numBinsX+1,yy);

                int lastBin     = hist.second->GetBin(numBinsX,yy);
                int overflowBin = hist.second->GetBin(numBinsX+1,yy);
                hist.second->AddBinContent(lastBin,overflowContent); // add overflow to last bin
                hist.second->SetBinContent(overflowBin,0);           // set overflow to 0
            }
        } // end 2D histogram overflow
    } // end else put overflow in first bin

    return;
}

void histogrammer::underFlow() {
    /* Add underflow to first bin */
    if (!m_putUnderflowInFirstBin){
        cma::INFO("HISTOGRAMMER : Not putting underflow in first bin(s)");
        return;
    }
    else{
        // loop over 1D histograms
        for (const auto& hist : m_map_histograms1D){
            double underflowContent = hist.second->GetBinContent(0);

            hist.second->AddBinContent(1, underflowContent);  // add underflow to first bin
            hist.second->SetBinContent(0, 0);                 // set underflow to 0
        }
        // loop over 2D histograms
        for (const auto& hist : m_map_histograms2D){
            unsigned int numBinsX = hist.second->GetXaxis()->GetNbins();
            unsigned int numBinsY = hist.second->GetYaxis()->GetNbins();

            // x-axis :: underflow in y
            for (unsigned int xx=1;xx<numBinsX+1;++xx){
                double underflowContent = hist.second->GetBinContent(xx,numBinsY+1);

                int firstBin     = hist.second->GetBin(xx,1);
                int underflowBin = hist.second->GetBin(xx,0);
                hist.second->AddBinContent(firstBin,underflowContent); // add overflow to last bin
                hist.second->SetBinContent(underflowBin,0);            // set overflow to 0
            }
            // y-axis :: underflow in x
            for (unsigned int yy=1;yy<numBinsY;++yy){
                double underflowContent = hist.second->GetBinContent(0,yy);

                int firstBin     = hist.second->GetBin(1,yy);
                int underflowBin = hist.second->GetBin(0,yy);
                hist.second->AddBinContent(firstBin,underflowContent); // add overflow to last bin
                hist.second->SetBinContent(underflowBin,0);           // set overflow to 0
            }
        } // end 2D histogram underflow
    } // end else put underflow in first bin
}

// THE END