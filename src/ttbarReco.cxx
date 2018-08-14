/*
Created:        19 February 2018
Last Updated:   19 May      2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for performing deep learning tasks
*/
#include "Analysis/goldilocks/interface/ttbarReco.h"


ttbarReco::ttbarReco( configuration& cmaConfig ) :
  m_config(&cmaConfig){
    m_targetMap = m_config->mapOfTargetValues();
  }

ttbarReco::~ttbarReco() {}


void ttbarReco::execute(const std::vector<Jet>& jets, const std::vector<Ljet>& ljets){
    /* Build top quarks system
       Testing semi-resolved tagger; only interested in AK8(QB/W)+AK4
        0 :: NONE = QCD (background)
        1 :: QB-Q = Signal AK8(QB) + AK4(Q)
        2 :: QQ-B = Signal AK8(W)  + AK4(B)
        3 :: FULL = Signal AK8
        4 :: QQ+X/QB+X (Other = pseudo-signal)

        -- if (fully contained): top_cand.target = 'full'; no extra AK4 needed
        -- if (q/b-only) top_cand.target = 4; need two extra AK4 -> this is effectively the 'fully resolved' case
    */
    m_ttbar.clear();

    m_jets  = jets;
    m_ljets = ljets;


    // Reconstruct ttbar, define containment
    bool isQCD(m_config->isQCD());
    bool isTtbar(m_config->isTtbar());

    unsigned int target(4);

    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(m_ljets.size())+" ak8 candidates");
    for (const auto& ljet : ljets){
        // Overlap Removal (subtract AK4 if matched to AK8 subjet)
        std::vector<Jet> ak4candidates;
        overlapRemoval(ljet,ak4candidates,isTtbar);     // overlap removal of AK4 from AK8

        int absLjetCt = std::abs(ljet.containment);

        cma::DEBUG("TTBARRECO : Access AK8 jet matched to parton "+std::to_string(ljet.matchId)+"; containment = "+std::to_string(absLjetCt));

        for (const auto& jet : ak4candidates){

            Top top_cand;  // reconstructed top candidates
            top_cand.jets.clear();
            top_cand.ljet = ljet.index;

            if (isTtbar && (absLjetCt==BQ || absLjetCt==W)){
                cma::DEBUG("TTBARRECO : -- AK8 that is QB or W");

                int total_containment = jet.containment+ljet.containment;

                if (total_containment==FULL || total_containment==-FULL){
                    cma::DEBUG("TTBARRECO : AK8+AK4 top in ttbar");
                    target = (absLjetCt==BQ) ? m_targetMap.at("BQ") : m_targetMap.at("W");
                }
                else{
                    cma::DEBUG("TTBARRECO : AK8+AK4 background combo in ttbar");
                    target = m_targetMap.at("ttbckg");
                }
            } // end if ljet is QB or W
            else if (isQCD){
                // Just assign AK8+AK4 jets as background (target 0)
                target = m_targetMap.at("none");
            } // end if QCD
            else continue;     // only want specific ttbar scenarios & qcd

            top_cand.jets.push_back(jet.index);
            top_cand.target = target;
            m_ttbar.push_back( top_cand );
        } // end loop over AK4
    } // end loop over ak8 candidates

    cma::DEBUG("TTBARRECO : Ttbar built ");

    return;
}


void ttbarReco::overlapRemoval(const Ljet& ak8, std::vector<Jet>& new_objects, const bool isTtbar){
    /* 
      Remove AK4 jets that overlap with AK8 candidate 
      'match_truth' for checking the AK4 and AK8 match the same parton
    */
    new_objects.clear();

    for (const auto& jet : m_jets){
        unsigned int dRmatch(0);

        // Don't use an AK4 jet that matches to the same individual parton (q,q,b) as the AK8
        unsigned int n_parton_overlaps(0);
        if (isTtbar){
            for (const auto& tp : jet.truth_partons){
                if (std::find(ak8.truth_partons.begin(), ak8.truth_partons.end(), tp) != ak8.truth_partons.end()) n_parton_overlaps++;
            }
        }
        if (n_parton_overlaps>0) continue;


        if (ak8.subjets.size()<1)  // match directly to AK8 if there are no subjets
            dRmatch = cma::deltaRMatch(ak8.p4, jet.p4, ak8.radius);
        else{
            for (const auto& sj : ak8.subjets)
                dRmatch += cma::deltaRMatch(sj.p4, jet.p4, sj.radius);
        }

        if (dRmatch<1) new_objects.push_back( jet );
    }

    return;
}


// THE END //
