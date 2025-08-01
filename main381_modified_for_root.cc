// main381.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: Higgs; electron-positron;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Simple example of Higgs pruduction at future e+e- colliders.

#include "Pythia8/Pythia.h"

//include ROOT functions to generate histograms
#include "TH1D.h"
#include "TFile.h"
#include "TApplication.h"
#include "TVirtualPad.h"
#include "TMath.h"
#include "TTree.h"


using namespace Pythia8;

//==========================================================================



int main() {
    
    // Number of events.
    int nEvent = 25000;
    
    // Generator. Incoming beams. (Switch off iniial-state photon radiation.)
    Pythia pythia;
    pythia.readString("Beams:idA = -11");
    pythia.readString("Beams:idB = 11");
    pythia.readString("Beams:eCM = 240.");
    //pythia.readString("PDF:lepton = off");

    // All Higgs production channels.
    pythia.readString("HiggsSM:all = on");
    pythia.readString("25:onMode = 0");
    pythia.readString("25:onifAny = 54");
    pythia.readString("23:onMode = 1");
    pythia.readString("51:all = S S 1 0 0 50.");
    pythia.readString("51:mWidth = 3.9466e-16");
    pythia.readString("51:tau0 = 500.");
    pythia.readString("25:addChannel = 1 1.0 100 51 51");
    pythia.readString("51:mayDecay = 1");
    pythia.readString("51:onMode = 0");
    pythia.readString("51:addChannel = 1 1.0 100 11 -11");
    //pythia.readString("23:onMode = off");
    //pythia.readString("23:onIfAny = 1 2 3 4 5");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;
    
    // Create file on which histogram(s) can be saved.
    TFile* outFile = new TFile("first_attempt_500mm.root", "RECREATE");

    TTree *t1 = new TTree("t1","t1");
    
    //Set up the TTree and define all the variables needed
    
    int pid;
    int MC_event;
    double x,y,z,t,energy,phi,theta,px,py,pz;
    std::vector<int> MotherList,SisterList,DaughterListRec;
    int list,status;
    bool isFinal;
    
    
    t1->Branch("energy",&energy,"energy/D");
    t1->Branch("x",&x,"x/D");
    t1->Branch("y",&y,"y/D");
    t1->Branch("z",&z,"z/D");
    t1->Branch("t",&t,"t/D");
    t1->Branch("pid",&pid,"pid/I");
    t1->Branch("phi",&phi,"phi/D");
    t1->Branch("theta",&theta,"theta/D");
    t1->Branch("px",&px,"px/D");
    t1->Branch("py",&py,"py/D");
    t1->Branch("pz",&pz,"pz/D");
    t1->Branch("MC_event",&MC_event,"MC_event/I");
    t1->Branch("MotherList",&MotherList);
    t1->Branch("list",&list);
    t1->Branch("status",&status);
    t1->Branch("isFinal",&isFinal);
    t1->Branch("SisterList",&SisterList);
    t1->Branch("DaughterListRec",&DaughterListRec);
    

    
    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;
        for (int ipt = 0; ipt < pythia.event.size();ipt++){
            pid = pythia.event[ipt].id();
            x = pythia.event[ipt].xProd();
            y = pythia.event[ipt].yProd();
            z = pythia.event[ipt].zProd();
            t = pythia.event[ipt].tProd();
            phi = pythia.event[ipt].phi();
            theta = pythia.event[ipt].theta();
            energy = pythia.event[ipt].e();
            px = pythia.event[ipt].px();
            py = pythia.event[ipt].py();
            pz = pythia.event[ipt].pz();
            list = ipt;
            SisterList = pythia.event[ipt].sisterList();
            status = pythia.event[ipt].status();
            MC_event = iEvent;
            MotherList = pythia.event[ipt].motherList();
            DaughterListRec = pythia.event[ipt].daughterListRecursive();
            isFinal = pythia.event[ipt].isFinal();
            t1->Fill();
            

        }
        
        
        
      }
    // Statistics on event generation.
    pythia.stat();

    // Write everything into a root file
    t1->Write();
    outFile->Close();
    
    delete outFile;
    return 0;
    
}
