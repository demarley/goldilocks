# Created: 10 January 2018
#
# Dan Marley
# daniel.edison.marley@cernSPAMNOT.ch
# Texas A&M University
# ---
# Tar the current CMSSW directory for running condor jobs


set eos_path = "/store/user/demarley/susy/semi-resolved-tagging/goldilocks/"

echo " Make tarball ..."

tar --exclude-caches-all --exclude-vcs -zcf $CMSSW_VERSION.tgz -C $CMSSW_BASE/.. $CMSSW_VERSION \
    --exclude=$CMSSW_VERSION/src/Analysis/goldilocks/data \
    --exclude=$CMSSW_VERSION/src/Analysis/goldilocks/batch \
    --exclude=$CMSSW_VERSION/src/Analysis/goldilocks/plots \
    --exclude=$CMSSW_VERSION/src/Analysis/goldilocks/training \
    --exclude=$CMSSW_VERSION/src/Analysis/goldilocks/inference \
    --exclude=$CMSSW_VERSION/lib/slc6_amd64_gcc530/.edmplugincache \
    --exclude=$CMSSW_VERSION/lib/slc6_amd64_gcc530/.poisonededmplugincache \
    --verbose


# eosrm "$eos_path"$CMSSW_VERSION.tgz
echo " Copy tarball to new location on EOS "$eos_path
xrdcp $CMSSW_VERSION.tgz root://cmseos.fnal.gov/"$eos_path"
