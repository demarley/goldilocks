# Goldilocks
The _not too boosted, not too resolved_ top tagger

Tagging top quarks in the semi-merged regime (moderate pT).  
Built for CMS analyses considering a wide pT-range of top quark decays:  
A large-R jet (AK8) is paired with a small-R jet (AK4) to reconstruct the top quark.

The following describes steps necessary to setup the tagger in a CMSSW environment.  
_The current setup does not use the Producer/Analyzer model, thus porting this to a "Standalone" setup outside of CMSSW should be straightforward_

## Getting Started

The goldilocks framework has been developed in the CMSSW release `CMSSW_8_0_28_patch1`.
To begin, checkout the relevant packages, including the [hepPlotter](https://github.com/demarley/hepPlotter) module for goldilocks.

```
## setup CMSSW
cmsrel CMSSW_8_0_28_patch1
cd CMSSW_8_0_28_patch1/src/
cmsenv
git cms-init

## Lightweight NN (for running Keras models in C++)
mkdir lwtnn
git clone https://github.com/demarley/lwtnn.git -b CMSSW_8_0_X-compatible lwtnn/lwtnn

## Goldilocks and hepPlotter
mkdir Analysis
git clone https://github.com/demarley/hepPlotter.git Analysis/
git clone https://github.com/demarley/goldilocks.git Analysis/
```



# Comments or Questions
Please post an issue or submit a PR with any questions, comments, or suggestions.
