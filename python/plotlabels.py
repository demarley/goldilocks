"""
Extra labels for plots
"""
import os
import sys
from array import array

try:
    CMSSW_BASE = os.environ['CMSSW_BASE']
    from Analysis.hepPlotter.histogram1D import Histogram1D
    from Analysis.hepPlotter.histogram1D import Histogram2D
    import Analysis.hepPlotter.labels as hpl
    import Analysis.hepPlotter.tools as hpt
except KeyError:
    cwd = os.getcwd()
    print cwd
    hpd = cwd.rstrip("goldilocks")+"/hepPlotter/python/"
    if hpd not in sys.path:
        sys.path.insert(0,hpd)
    import tools as hpt


class Sample(object):
    """Class for organizing plotting information about physics samples"""
    def __init__(self,label='',color=''):
        self.label = label
        self.color = color

class Variable(object):
    """Class for organizing plotting information about variables"""
    def __init__(self,binning=[],label=''):
        self.binning = binning
        self.label   = label


def variable_labels():
    """Dictionaries that contain Variables objects."""
    _phi  = r'$\phi$'
    _eta  = r'$\eta$'
    _T    = r'$_{\text{T}}$ [GeV]'
    _mass = 'Mass [GeV]'
    bdisc_bins = array('d',[i*0.1 for i in range(11)])  # default value = -1

    variables = {}

    variables['AK8_C2']    = Variable(binning=hpt.hist1d(10,  0.,   0.6), label=r'C$_2^{\beta\text{=1}}$')
    variables['AK8_D2']    = Variable(binning=hpt.hist1d(20,  0.,   5.0), label=r'D$_2^{\beta\text{=1}}$')
    variables['AK8_d12']   = Variable(binning=hpt.hist1d(20,  0.,  125.), label=r'$\sqrt{\text{d}_{\text{12}}}$ [GeV]')
    variables['AK8_d23']   = Variable(binning=hpt.hist1d(12,  0.,   60.), label=r'$\sqrt{\text{d}_{\text{23}}}$ [GeV]')
    variables['AK8_eta']   = Variable(binning=hpt.hist1d(20, -3.,    3.), label=r'AK8 '+_eta)
    variables['AK8_phi']   = Variable(binning=hpt.hist1d(20, -2.,    2.), label=r'AK8 $\phi$')
    variables['AK8_m']     = Variable(binning=hpt.hist1d(40,  0.,  400.), label=r'AK8 '+_mass)
    variables['AK8_pt']    = Variable(binning=hpt.hist1d(14,200., 1500.), label=r'AK8 p'+_T)
    variables['AK8_tau1']  = Variable(binning=hpt.hist1d(10,  0.,   0.6), label=r'$\tau_{\text{1}}$')
    variables['AK8_tau2']  = Variable(binning=hpt.hist1d(10,  0.,   0.5), label=r'$\tau_{\text{2}}$')
    variables['AK8_tau21'] = Variable(binning=hpt.hist1d(11, 00.,   1.1), label=r'$\tau_{\text{21}}$')
    variables['AK8_tau3']  = Variable(binning=hpt.hist1d(10,  0.,   0.6), label=r'$\tau_{\text{3}}$')
    variables['AK8_tau32'] = Variable(binning=hpt.hist1d(11,  0.,   1.1), label=r'$\tau_{\text{32}}$')
    variables['AK8_softDropMass'] = Variable(binning=hpt.hist1d(40,0.,400.), label=r'AK8 '+_mass)
    variables['AK8_SDMass'] = variables['AK8_softDropMass']
    variables['AK8_subjet0_bdisc'] = Variable(binning=hpt.hist1d(10,0,1), label=r'AK8 Subjet 0 bDisc')
    variables['AK8_subjet0_pTrel'] = Variable(binning=hpt.hist1d(10,0,1), label=r'AK8 Subjet 0 p$_{\text{T}}^{\text{rel}}$')
    variables['AK8_subjet1_bdisc'] = Variable(binning=hpt.hist1d(10,0,1), label=r'AK8 Subjet 1 bDisc')
    variables['AK8_subjet1_pTrel'] = Variable(binning=hpt.hist1d(10,0,1), label=r'AK8 Subjet 1 p$_{\text{T}}^{\text{rel}}$')

    for i in range(16):
        variables['AK8_deepAK8_{0}'.format(i)] = Variable(binning=hpt.hist1d(10,0,1), label=r'DeepAK8[{0}]'.format(i))

    variables['AK4_deepCSVb']  = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepCSV(b)')
    variables['AK4_deepCSVbb'] = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepCSV(bb)')
    variables['AK4_deepCSVc']  = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepCSV(c)')
    variables['AK4_deepCSVcc'] = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepCSV(cc)')
    variables['AK4_deepCSVl']  = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepCSV(l)')

    variables['AK4_deepFlavorb']    = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(b)')
    variables['AK4_deepFlavorbb']   = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(bb)')
    variables['AK4_deepFlavorc']    = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(c)')
    variables['AK4_deepFlavoruds']  = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(uds)')
    variables['AK4_deepFlavorg']    = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(g)')
    variables['AK4_deepFlavorlepb'] = Variable(binning=hpt.hist1d(10,0,1),label=r'AK4 DeepFlavor(lepb)')

    variables['AK8AK4_mass']   = Variable(binning=hpt.hist1d(40,0.,400.), label=r'AK8+AK4 '+_mass)
    variables['AK8AK4_deltaR'] = Variable(binning=hpt.hist1d(10,0.,5.),   label=r'$\Delta$R(AK8,AK4)')

    return variables



def sample_labels():
    """Dictionaries that contain Samples objects.
       > The key values match those in config/sampleMetadata.txt.
         (The functions in util.py are setup to read the information this way.)
         If you want something unique, then you may need to specify 
         it manually in your plotting script
    """
    ## Sample information
    samples = {}

    samples['signal'] = Sample(label='Signal',color='b')
    samples['bckg']   = Sample(label='Bckg',color='r')

    ttbar = r't$\bar{\text{t}}$'
    samples['multijet'] = Sample(label=r'Multi-jet',  color='purple')
    samples['BQ']       = Sample(label=ttbar+' (QB)', color='red')
    samples['W']        = Sample(label=ttbar+' (W)',  color='blue')
    samples['ttbckg']   = Sample(label=ttbar+' bckg.',color='green')

    return samples
