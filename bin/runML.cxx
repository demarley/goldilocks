/*
Created:        17 February 2018
Last Updated:   28 May      2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Basic steering macro for running goldilocks
 - Make flat ntuples for machine learning
 - Make histograms to compare features
*/
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <map>
#include <fstream>
#include <string>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Analysis/goldilocks/interface/configuration.h"
#include "Analysis/goldilocks/interface/Event.h"
#include "Analysis/goldilocks/interface/eventSelection.h"
#include "Analysis/goldilocks/interface/flatTree4ML.h"
#include "Analysis/goldilocks/interface/tools.h"
#include "Analysis/goldilocks/interface/histogrammer4ML.h"


int main(int argc, char** argv) {
    /* Steering macro for goldilocks */
    if (argc < 2) {
        cma::HELP();
        return -1;
    }

    unsigned long long maxEntriesToRun(0);     // maximum number of entries in TTree
    unsigned int numberOfEventsToRun(0);       // number of events to run
    bool passEvent(false);                     // event passed selection

    // configuration
    configuration config(argv[1]);                         // configuration file
    config.initialize();

    int nEvents = config.nEventsToProcess();                      // requested number of events to run
    std::string outpathBase = config.outputFilePath();            // directory for output files
    unsigned long long firstEvent      = config.firstEvent();     // first event to begin running over
    std::vector<std::string> filenames = config.filesToProcess(); // list of files to process
    std::string selection(config.selection());                    // selection to apply
    std::string treename(config.treename());

    std::string customFileEnding( config.customFileEnding() );
    if (customFileEnding.length()>0  && customFileEnding.substr(0,1).compare("_")!=0){
        customFileEnding = "_"+customFileEnding; // add '_' to beginning of string, if needed
    }

    // event selection
    eventSelection evtSel( config );
    evtSel.initialize();
    unsigned int ncuts = evtSel.numberOfCuts();            // number of cuts in selection
    std::vector<std::string> cutNames = evtSel.cutNames(); // names of cuts


    // --------------- //
    // -- File loop -- //
    // --------------- //
    unsigned int numberOfFiles(filenames.size());
    unsigned int currentFileNumber(0);
    cma::INFO("RUNML : *** Starting file loop *** ");
    for (const auto& filename : filenames) {

        ++currentFileNumber;
        cma::INFO("RUNML :   Opening "+filename+"   ("+std::to_string(currentFileNumber)+"/"+std::to_string(numberOfFiles)+")");

        auto file = TFile::Open(filename.c_str());
        if (!file || file->IsZombie()){
            cma::WARNING("RUNML :  -- File: "+filename);
            cma::WARNING("RUNML :     does not exist or it is a Zombie. ");
            cma::WARNING("RUNML :     Continuing to next file. ");
            continue;
        }

        config.setFilename( filename );   // Use the filename to determine primary dataset and information about the sample
        config.inspectFile( *file );      // Determine information about the input file (metadata)
        Sample s = config.sample();       // load the Sample struct (xsection,kfactor,etc)

        std::vector<std::string> fileKeys;
        cma::getListOfKeys(file,fileKeys);      // keep track of ttrees in file


        // -- Output file -- //
        // CMS doesn't use 'mcChannelNumber', need to keep the same file names
        // therefore, make new directories for the different selections
        struct stat dirBuffer;
        std::string outpath = outpathBase+"/"+selection+customFileEnding;
        if ( !(stat((outpath).c_str(),&dirBuffer)==0 && S_ISDIR(dirBuffer.st_mode)) ){
            cma::DEBUG("RUNML : Creating directory for storing output: "+outpath);
            system( ("mkdir "+outpath).c_str() );  // make the directory so the files are grouped together
        }

        std::size_t pos   = filename.find_last_of(".");     // the last ".", i.e., ".root"
        std::size_t found = filename.find_last_of("/");     // the last "/"
        std::string outputFilename = filename.substr(found+1,pos-1-found); // betwee "/" and "."
        // hopefully this returns: "diboson_WW_361082" given something like:
        // "/some/path/to/file/diboson_WW_361082.root"

        std::string fullOutputFilename = outpath+"/"+outputFilename+".root";
        std::unique_ptr<TFile> outputFile(TFile::Open( fullOutputFilename.c_str(), "RECREATE"));
        cma::INFO("RUNML :   >> Saving to "+fullOutputFilename);


        histogrammer4ML histMaker(config,"ML");      // initialize histogrammer
        histMaker.initialize( *outputFile );

        // -- Cutflow histograms
        std::map<std::string, TH1D*>  h_cutflows;            // map of cutflow histograms (weights applied)
        std::map<std::string, TH1D*>  h_cutflows_unweighted; // map of cutflow histograms (raw # of events)


        // check that the ttree exists in this file before proceeding
        if (std::find(fileKeys.begin(), fileKeys.end(), treename) == fileKeys.end()){
            cma::INFO("RUNML : TTree "+treename+" is not present in this file, continuing to next TTree");
            continue;
        }

        // -- Cutflow histogram [initialize and label bins]
        h_cutflows[treename] = new TH1D( (treename+"_cutflow").c_str(),(treename+"_cutflow").c_str(),ncuts+1,0,ncuts+1);
        h_cutflows_unweighted[treename] = new TH1D( (treename+"_cutflow_unweighted").c_str(),(treename+"_cutflow_unweighted").c_str(),ncuts+1,0,ncuts+1);

        h_cutflows[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");
        h_cutflows_unweighted[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");

        for (unsigned int c=1;c<=ncuts;++c){
            h_cutflows[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
            h_cutflows_unweighted[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
        }

        // -- Load TTree to loop over
        cma::INFO("RUNML :      TTree "+treename);
        TTreeReader myReader(treename.c_str(), file);

        // -- Make new Tree in Root file
        flatTree4ML miniTTree(config);          // initialize TTree for new file
        miniTTree.initialize( *outputFile );

        // -- Number of Entries to Process -- //
        maxEntriesToRun = myReader.GetEntries(true);
        if (maxEntriesToRun<1) // skip files with no entries
            continue;

        if (nEvents < 0 || ((unsigned int)nEvents+firstEvent) > maxEntriesToRun)
            numberOfEventsToRun = maxEntriesToRun - firstEvent;
        else
            numberOfEventsToRun = nEvents;

        // ---------------- //
        // -- Event Loop -- //
        // ---------------- //
        Long64_t imod = 1;                     // print to the terminal
        Event event = Event(myReader, config);

        Long64_t eventCounter = 0;    // counting the events processed
        Long64_t entry = firstEvent;  // start at a different event!
        while (myReader.Next()) {

            if (eventCounter+1 > numberOfEventsToRun){
                cma::INFO("RUNML : Processed the desired number of events: "+std::to_string(eventCounter)+"/"+std::to_string(numberOfEventsToRun));
                break;
            }

            if (entry%imod==0){
                cma::INFO("RUNML :       Processing event "+std::to_string(entry) );
                if(imod<2e4) imod *=10;
            }

            // -- Build Event -- //
            cma::DEBUG("RUNML : Execute event");
            event.execute(entry);
            // now we have event object that has the event-level objects in it
            // pass this to the selection tools

            // -- Event Selection -- //
            cma::DEBUG("RUNML : Apply event selection");
            passEvent = evtSel.applySelection(event,*h_cutflows.at( treename.c_str() ),*h_cutflows_unweighted.at( treename.c_str() ));

            std::vector<Top> tops = event.ttbar();
            if (passEvent && tops.size()>0){
                cma::DEBUG("RUNML : Passed selection, now save information");

                // For ML, we are training on semi-boosted top quarks (AK8+AK4)
                // Only save features of the AK8+AK4 system to the output ntuple/histograms
                std::map<std::string,double> features2save;
                features2save["xsection"] = s.XSection;
                features2save["kfactor"]  = s.KFactor;
                features2save["sumOfWeights"] = s.sumOfWeights;
                features2save["nominal_weight"] = event.nominal_weight();

                for (const auto& top : tops){

                    for (const auto& x : top.dnn) features2save[x.first] = x.second;

                    miniTTree.saveEvent(features2save);
                    histMaker.fill(features2save);
                }
            }

            // iterate the entry and number of events processed
            ++entry;
            ++eventCounter;
        } // end event loop

        event.finalize();
        miniTTree.finalize();

        // put overflow/underflow content into the first and last bins
        histMaker.overUnderFlow();

        cma::INFO("RUNML :   END Running  "+filename);
        cma::INFO("RUNML :   >> Output at "+fullOutputFilename);

        outputFile->Write();
        outputFile->Close();

        // -- Clean-up stuff
        delete file;          // free up some memory 
        file = ((TFile *)0);  // (no errors for too many root files open)
    } // end file loop

    cma::INFO("RUNML : *** End of file loop *** ");
    cma::INFO("RUNML : Program finished. ");
}

// THE END
