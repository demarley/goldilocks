# Goldilocks
The _not too boosted, not too resolved_ top tagger

Tagging hadronically-decaying top quarks (t->qqb) in the semi-merged regime (moderate transverse momentum).  
Built for CMS analyses considering a wide range of top quark decays:  
> A large-R jet (AK8) is paired with a small-R jet (AK4) to reconstruct the top quark.
> Goldilocks is applied to the AK8+AK4 system to classify it as one of the following categories:
> - Signal "W" (AK8='qq'; AK4='b')
> - Signal "QB" (AK8='qb'; AK4='q')
> - Top Background "tt bckg" (AK8='qq' or 'qb'; AK4=incorrect jet)
> - Multi-jet Background (AK8 and AK4 are chosen randomly from a non-top quark sample)


The following describes steps necessary to setup the tagger in a CMSSW environment.  
_The current setup does not use the Producer/Analyzer model, thus porting this to a "Standalone" setup outside of CMSSW should be straightforward_

## Getting Started

### CMSSW Environment
The goldilocks framework has been developed in the CMSSW release `CMSSW_8_0_28_patch1` to process flat ntuples and generate samples for training.

```
## setup CMSSW
cmsrel CMSSW_8_0_28_patch1
cd CMSSW_8_0_28_patch1/src/
cmsenv
git cms-init --upstream-only

## Lightweight NN (for running Keras models in C++)
mkdir lwtnn
git clone https://github.com/demarley/lwtnn.git -b CMSSW_8_0_X-compatible lwtnn/lwtnn

## Goldilocks
mkdir Analysis
git clone https://github.com/demarley/goldilocks.git Analysis/
```

### Python (Deep Learning) Environment

To perform the training, a python environment outside of CMSSW is used by the authors.
Two packages, 
[hepPlotter](https://github.com/demarley/hepPlotter) and [asimov](https://github.com/demarley/asimov), 
are necessary in addition to goldilocks:

```
# Checkout Asimov package (for performing the training)
git clone https://github.com/demarley/asimov.git

# Checkout hepPlotter package (for making plots)
git clone https://github.com/demarley/hepPlotter.git
```

To perform the training, please consult the wiki.


# Comments or Questions
Please post an issue or submit a PR with any questions, comments, or suggestions.
