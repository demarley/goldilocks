# goldilocks
The "not too boosted, not too resolved" top tagger

Tagging top quarks in the semi-merged regime.  
Built for CMS analyses considering a wide pT-range of top quark decays:  
A large-R jet (AK8) is paired with a small-R jet (AK4) to reconstruct the top quark.

The following describes steps necessary to setup the tagger in a CMSSW environment.  
_The current setup does not use the Producer/Analyzer model, thus porting this to a "Standalone" setup outside of CMSSW should be straightforward_

## Getting Started

The goldilocks framework has been developed in the CMSSW release `CMSSW_8_0_28_patch1`.
To begin, checkout the relevant packages, including the [hepPlotter]() submodule for goldilocks.

```
## setup CMSSW
cmsrel CMSSW_8_0_28_patch1
cd CMSSW_8_0_28_patch1/src/
cmsenv
git cms-init

## Lightweight NN (for running Keras models in C++)
mkdir lwtnn
git clone https://github.com/demarley/lwtnn.git -b CMSSW_8_0_X-compatible lwtnn/lwtnn

## This module
mkdir Analysis
git clone --recurse-submodules https://github.com/demarley/goldilocks.git Analysis/
```



# Comments or Questions
Contact the author
