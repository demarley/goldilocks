"""
Created:        11 November  2016
Last Updated:   16 February  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Base class for plotting deep learning

Designed for running on desktop at TAMU
with specific set of software installed
--> not guaranteed to work in CMSSW environment!

Does not use ROOT!
Instead, uses matplotlib to generate figures
"""
import os
import sys
import json
import util
from datetime import date
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='sans-serif')
from keras.utils.vis_utils import plot_model as keras_plot
from sklearn.metrics import roc_curve, auc

import hepPlotter.hepPlotterLabels as hpl
import hepPlotter.hepPlotterTools as hpt
from hepPlotter.hepPlotter import HepPlotter



class Target(object):
    """Class to contain information for targets used in training"""
    def __init__(self,name=""):
        self.name  = name     # Name of this target, e.g., 'signal'
        self.df    = None     # dataframe of this target's features
        self.color = 'k'
        self.label = ''
        self.target_value = -999
        self.binning = 1


class DeepLearningPlotter(object):
    """Plotting utilities for deep learning"""
    def __init__(self):
        """Give default values to member variables"""
        self.date = date.today().strftime('%d%b%Y')
        self.betterColors = hpt.betterColors()['linecolors']
        self.sample_labels   = hpl.sample_labels()
        self.variable_labels = hpl.variable_labels()

        self.msg_svc      = util.VERBOSE()
        self.filename     = ""
        self.output_dir   = ''
        self.image_format = 'png'
        self.process_label = ''      # if a single process is used for all training, set this

        self.classification = False  # 'binary','multi',False
        self.regression     = False  # True or False

        self.df = None
        self.targets = []
        self.baseline_target = None

        self.CMSlabelStatus = "Simulation Internal"


    def initialize(self,dataframe,target_names=[],target_values=[]):
        """
        Set parameters of class to make plots

        @param dataframe    The dataframe that contains physics information for training/testing
        """
        self.df = dataframe

        try:
            self.processlabel = self.sample_labels[self.filename].label   # process used in each plot
        except KeyError:
            self.processlabel = ''

        if self.classification:
            for i,(n,v) in enumerate(zip(target_names,target_values)):
                tmp    = Target(n)
                tmp.df = self.df.loc[self.df['target']==v]
                tmp.target_value = v
                tmp.label = self.sample_labels[n].label
                tmp.color = self.betterColors[i]
                self.targets.append(tmp)
        else: # regression
            try:
                tmp    = Target(target_names[0])
                tmp.df = self.df.loc[self.df['target']==target_values[0]]
                tmp.target_value = target_values[0]
            except TypeError:
                tmp    = Target(target_names)
                tmp.df = self.df.loc[self.df['target']==target_values]
                tmp.target_value = target_values

            tmp.label = self.sample_labels[tmp.name].label
            tmp.color = self.betterColors[i]
            self.targets.append(tmp)

        if self.baseline_target is None:
            for target in self.targets:
                if target.target_value==0: 
                    self.baseline_target = target
                    break

        return


    def features(self):
        """
        Plot the features
        For classification, compare different targets
        For regression, just plot the features        <- should do data/mc plots instead!
        """
        self.msg_svc.INFO("DL : Plotting features.")

        # pick one of the targets to compare against (critical for multi-classification)
        # -- compare with the 0 target (main background)
        ratio_denominator = {}
        ratio_numerator   = {}
        ratio_partners    = {}
        for t,target in enumerate(self.targets):
            ratio_denominator[target.name] = False
            ratio_numerator[target.name] = False
            ratio_partners[target.name] = []

            if target.target_value==self.baseline_target.target_value: 
                ratio_denominator[target.name] = True
                ratio_partners[target.name]    = [i.name for i in self.targets if i.target_value!=self.baseline_target.target_value ]
            else:
                ratio_numerator[target.name] = True
                ratio_partners[target.name] = [ i.name for i in self.targets if i.target_value==self.baseline_target.target_value ]

        # plot the features
        plt_feature = self.df.keys()
        for hi,feature in enumerate(plt_feature):

            if feature=='target': continue

            hist = HepPlotter("histogram",1)

            hist.normed  = True
            hist.stacked = False
            hist.logplot = False
            hist.binning = self.variable_labels[feature].binning
            hist.x_label = self.variable_labels[feature].label
            hist.y_label = "Events"
            hist.format  = self.image_format
            hist.saveAs  = self.output_dir+"/hist_"+feature+"_"+self.date
            hist.ratio_plot  = True
            hist.ratio_type  = 'significance'
            hist.y_ratio_label = r"S/$\sqrt{\text{B}}$"
            hist.CMSlabel    = 'top left'
            hist.CMSlabelStatus   = self.CMSlabelStatus
            hist.numLegendColumns = 1

            hist.initialize()

            # Add some extra text to the plot
            n_extra_text = 0
            if self.processlabel: 
                hist.extra_text.Add(self.processlabel,coords=[0.03,0.80]) # physics process that produces these features
                n_extra_text+=1

            base_data,_ = np.histogram(self.baseline_target.df[feature],bins=self.variable_labels[feature].binning,normed=True)
            for t,target in enumerate(self.targets):
                ratio_den = ratio_denominator[target.name]
                ratio_num = ratio_numerator[target.name]
                ratio_partner = ratio_partners[target.name]

                hist.Add(target.df[feature], name=target.name, draw='step',
                         linecolor=target.color, label=target.label,
                         ratio_partner=ratio_partner,ratio_den=ratio_den,ratio_num=ratio_num)
                if target.name!=self.baseline_target.name:
                    t_data,_ = np.histogram(target.df[feature],bins=self.variable_labels[feature].binning,normed=True)
                    separation = util.getSeparation(base_data,t_data)
                    hist.extra_text.Add("Sep({0},{1}) = {2:.2f}".format(self.baseline_target.name,target.name,separation),
                                        coords=[0.03,(0.80-0.08*n_extra_text)])
                    n_extra_text+=1

            p = hist.execute()
            hist.savefig()

        return


    def feature_correlations(self):
        """Plot correlations between features of the NN"""
        ## Correlation Matrices of Features (top/antitop) ##
        fontProperties = {'family':'sans-serif'}
        opts = {'cmap': plt.get_cmap("bwr"), 'vmin': -1, 'vmax': +1}

        for c,target in enumerate(self.targets):

            saveAs = "{0}/correlations_{1}_{2}".format(self.output_dir,target.name,self.date)

            allkeys = target.df.keys()
            keys = []
            for key in allkeys:
                if key!='target': keys.append(key)
            t_ = target.df[keys]
            corrmat = t_.corr()

            # Save correlation matrix to CSV file
            corrmat.to_csv("{0}.csv".format(saveAs))

            # Plot correlation matrix
            # -- Use matplotlib directly
            fig,ax = plt.subplots()

            heatmap1 = ax.pcolor(corrmat, **opts)
            cbar     = plt.colorbar(heatmap1, ax=ax)

            cbar.ax.set_yticklabels( [i.get_text().strip('$') for i in cbar.ax.get_yticklabels()], **fontProperties )

            labels = [self.variable_labels[feature].label for feature in corrmat.columns.values]
            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_xticklabels(labels, fontProperties, fontsize=12, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(labels, fontProperties, fontsize=12, minor=False)


            ## CMS/COM Energy Label + Signal name
            cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
            cms_stamp.coords = [0.02,1.00]
            cms_stamp.fontsize = 16
            cms_stamp.va = 'bottom'
            ax.text(0.02,1.00,cms_stamp.text,fontsize=cms_stamp.fontsize,
                    ha=cms_stamp.ha,va=cms_stamp.va,transform=ax.transAxes)

            energy_stamp    = hpl.EnergyStamp()
            energy_stamp.ha = 'right'
            energy_stamp.coords = [0.99,1.00]
            energy_stamp.fontsize = 16
            energy_stamp.va = 'bottom'
            ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                    fontsize=energy_stamp.fontsize,ha=energy_stamp.ha, va=energy_stamp.va, transform=ax.transAxes)

            ax.text(0.03,0.93,target.label,fontsize=16,ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format),
                        format=self.image_format,dpi=300,bbox_inches='tight')
            plt.close()

        return


    def prediction(self,train_data={},test_data={}):
        """Plot the training and testing predictions"""
        self.msg_svc.INFO("DL : Plotting DNN prediction. ")

        values = {0:np.array([1,0,0]), 1:np.array([0,1,0]), 2:np.array([0,0,1])}

        # Plot all k-fold cross-validation results
        for i,(train,trainY,test,testY) in enumerate(zip(train_data['X'],train_data['Y'],test_data['X'],test_data['Y'])):

            # Make a plot for each target value (e.g., QCD prediction to be QCD; Top prediction to be QCD, etc)
            for t_ind,tar in enumerate(self.targets):

                hist = HepPlotter("histogram",1)

                hist.ratio_plot    = True
                hist.ratio_type    = "ratio"
                hist.y_ratio_label = "Test/Train"
                hist.label_size    = 14
                hist.normed  = True  # compare shape differences (likely don't have the same event yield)
                hist.format  = self.image_format
                hist.saveAs  = "{0}/hist_DNN_prediction_kfold{1}_target{2}_{3}".format(self.output_dir,i,t_ind,self.date)
                hist.binning = [bb/10. for bb in range(11)]
                hist.stacked = False
                hist.logplot = False
                hist.x_label = "Prediction"
                hist.y_label = "Arbitrary Units"
                hist.CMSlabel = 'top left'
                hist.CMSlabelStatus   = self.CMSlabelStatus
                hist.numLegendColumns = 1

                hist.extra_text.Add("{0} Prediction".format(tar.name),coords=[0.03,0.80],fontsize=14)
                if self.processlabel: hist.extra_text.Add(self.processlabel,coords=[0.03,0.72],fontsize=14)

                hist.initialize()

                json_data = {}
                for t,target in enumerate(self.targets):

                    target_value = values[target.target_value]  # arrays for multiple predictions 

                    train_t = train[np.where((trainY==target_value).all(axis=1))][:,t_ind] # get the t_ind prediction for this target 
                    test_t  = test[np.where((testY==target_value).all(axis=1))][:,t_ind]

                    ## Training
                    hist.Add(train_t, 
                             name=target.name+'_train', linecolor=target.color,
                             linewidth=2, draw='step', label=target.label+" Train",
                             ratio_den=True,ratio_num=False,ratio_partner=target.name+'_test')
                    ## Testing
                    hist.Add(test_t, 
                             name=target.name+'_test', linecolor=target.color, color=target.color,
                             linewidth=0, draw='stepfilled', label=target.label+" Test", alpha=0.5,
                             ratio_den=False,ratio_num=True,ratio_partner=target.name+'_train')

                    ## Save data to JSON file
                    json_data[target.name+"_train"] = {}
                    json_data[target.name+"_test"]  = {}
                    d_tr,b_tr = np.histogram(train_t,bins=hist.binning)
                    d_te,b_te = np.histogram(test_t,bins=hist.binning)

                    json_data[target.name+"_train"]["binning"] = b_tr.tolist()
                    json_data[target.name+"_train"]["content"] = d_tr.tolist()
                    json_data[target.name+"_test"]["binning"] = b_te.tolist()
                    json_data[target.name+"_test"]["content"] = d_te.tolist()

#                separation = util.getSeparation(base_data,t_data)
#                hist.extra_text.Add("Sep({0},{1}) = {2:.2f}".format(self.baseline_target.name,target.name,separation),
#                                    coords=[0.03,(0.80-0.08*n_extra_text)])

                p = hist.execute()
                hist.savefig()

                # save results to JSON file (just histogram values & bins) to re-make plots
                with open("{0}.json".format(hist.saveAs), 'w') as outfile:
                    json.dump(json_data, outfile)

        return



    def ROC(self,fprs=[],tprs=[],accuracy={}):
        """Plot the ROC curve & save to text file"""
        self.msg_svc.INFO("DL : Plotting ROC curve.")

        saveAs = "{0}/roc_curve_{1}".format(self.output_dir,self.date)

        ## Use matplotlib directly
        fig,ax = plt.subplots()

        # Draw all of the ROC curves from the K-fold cross-validation
        ax.plot([0, 1], [0, 1], ls='--',label='No Discrimination',lw=2,c='gray')
        ax.axhline(y=1,lw=1,c='lightgray',ls='--')

        for ft,(fpr,tpr) in enumerate(zip(fprs,tprs)):
            # Plot ROC curve
            roc_auc = auc(fpr,tpr)
            ax.plot(fpr,tpr,label='K-fold {0} (AUC = {1:.2f})'.format(ft,roc_auc),lw=2)

            # save ROC curve to CSV file (to plot later)
            outfile_name = "{0}_{1}.csv".format(saveAs,ft)
            csv = [ "{0},{1}".format(fp,tp) for fp,tp in zip(fpr,tpr) ]
            util.to_csv(outfile_name,csv)

        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.5])

        ax.set_xlabel(r'$\epsilon$(anti-top)',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(["{0:.1f}".format(i) for i in ax.get_xticks()],fontsize=22)
        ax.set_ylabel(r'$\epsilon$(top)',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+["{0:.1f}".format(i) for i in ax.get_yticks()[1:-1]]+[''],fontsize=22)

        ## CMS/COM Energy Label
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        cms_stamp.coords = [0.03,0.97]
        cms_stamp.fontsize = 16
        ax.text(cms_stamp.coords[0],cms_stamp.coords[1],cms_stamp.text,fontsize=cms_stamp.fontsize,
                ha=cms_stamp.ha,va=cms_stamp.va,transform=ax.transAxes)

        energy_stamp    = hpl.EnergyStamp()
        energy_stamp.coords = [0.03,0.90]
        energy_stamp.fontsize = 16
        ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                fontsize=energy_stamp.fontsize,ha=energy_stamp.ha, va=energy_stamp.va, transform=ax.transAxes)

        text_args = {'ha':'left','va':'top','fontsize':18,'transform':ax.transAxes}
        if self.processlabel: ax.text(0.03,0.82,self.processlabel,**text_args)
        if accuracy: ax.text(0.03,0.75,r"Accuracy = {0:.2f}$\pm${1:.2f}".format(accuracy['mean'],accuracy['std']),**text_args)

        leg = ax.legend(loc=4,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig('{0}.{1}'.format(saveAs,self.image_format),
                    format=self.image_format,bbox_inches='tight',dpi=300)
        plt.close()

        return


    def plot_loss_history(self,history,ax=None,index=-1):
        """Draw history of model"""
        loss  = history.history['loss']
        x     = range(1,len(loss)+1)
        label = 'Loss {0}'.format(index) if index>=0 else 'Loss'
        ax.plot(x,loss,label=label)

        csv = [ "{0},{1}".format(i,j) for i,j in zip(x,loss) ]

        return csv


    def loss_history(self,history,kfold=0,val_loss=0.0):
        """Plot loss as a function of epoch for model"""
        self.msg_svc.INFO("DL : Plotting loss as a function of epoch number.")

        saveAs = "{0}/loss_epochs_{1}".format(self.output_dir,self.date)
        all_histories = type(history)==list

        # draw the loss curve
        fig,ax = plt.subplots()

        # also save the data to a CSV file
        if all_histories:
            for i,h in enumerate(history):
                csv = self.plot_loss_history(h,ax=ax,index=i)
                filename = "{0}_{1}.csv".format(saveAs,i)
                util.to_csv(filename,csv)
        else:
            csv = self.plot_loss_history(history,ax=ax)
            filename = "{0}.csv".format(saveAs)
            util.to_csv(filename,csv)

        ax.set_xlabel('Epoch',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(["{0:.1f}".format(i) for i in ax.get_xticks()],fontsize=22)
        ax.set_ylabel('Loss',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+["{0:.1f}".format(i) for i in ax.get_yticks()[1:-1]]+[''],fontsize=22)


        ## CMS/COM Energy Label
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        cms_stamp.coords = [0.03,0.97]
        cms_stamp.fontsize = 18
        ax.text(cms_stamp.coords[0],cms_stamp.coords[1],cms_stamp.text,fontsize=cms_stamp.fontsize,
                ha=cms_stamp.ha,va=cms_stamp.va,transform=ax.transAxes)

        energy_stamp    = hpl.EnergyStamp()
        energy_stamp.coords = [0.03,0.90]
        energy_stamp.fontsize = 18
        ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                fontsize=energy_stamp.fontsize,ha=energy_stamp.ha, va=energy_stamp.va, transform=ax.transAxes)

        text_args = {'ha':'left','va':'top','fontsize':18,'transform':ax.transAxes}
        text = "Validation Loss = {0}; {1} K-folds".format(val_loss,len(history)) if all_histories else "Validation Loss = {0}".format(val_loss)
        ax.text(0.03,0.76,text,**text_args)

        leg = ax.legend(loc=1,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig(self.output_dir+'/loss_epochs_{0}_{1}.{2}'.format(kfold,self.date,self.image_format),
                    format=self.image_format,bbox_inches='tight',dpi=200)
        plt.close()


        return


    def model(self,model,name):
        """Plot the model architecture to view later"""
        keras_plot(model,to_file='{0}/{1}_model.eps'.format(self.output_dir,name),show_shapes=True)
        return

## THE END ##
