#ifndef TTBARRECO_H
#define TTBARRECO_H

#include <string>
#include <map>
#include <vector>

#include "Analysis/goldilocks/interface/tools.h"
#include "Analysis/goldilocks/interface/configuration.h"
#include "Analysis/goldilocks/interface/physicsObjects.h"


class ttbarReco {
  public:
    ttbarReco( configuration& cmaConfig );

    ~ttbarReco();

    std::vector<Top> tops();
    void execute(const std::vector<Jet>& jets, const std::vector<Ljet>& ljets);
    void overlapRemoval(const Ljet& ak8, std::vector<Jet>& new_objects, const bool match_truth);

  protected:

    configuration *m_config;

    std::vector<Top> m_ttbar;
    std::map<std::string,int> m_mapContainment;
    std::map<std::string,int> m_targetMap;

    std::vector<Jet> m_jets;
    std::vector<Ljet> m_ljets;
};

#endif
