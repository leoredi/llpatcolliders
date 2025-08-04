import ROOT as r
r.gROOT.SetBatch()

inputFile = r.TFile("main144.root")

inputTree = inputFile.Get("t")

for event in inputTree:
    print (event.particles[0].phi)

