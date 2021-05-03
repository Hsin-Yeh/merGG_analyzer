#!/usr/bin/env python

from __future__ import division

import ROOT
import os, sys, logging
ROOT.gROOT.SetBatch(True)
import argparse
import math
from array import array

def drawHistos(h, xaxis_name, plot_dir, plot_name, log_flag):
    c1 = ROOT.TCanvas()
    c1.SetGrid()
    h[0].GetYaxis().SetTitle("a.u.")
    h[0].GetXaxis().SetTitle("{}".format(xaxis_name))
    h[0].SetLineColor(209)
    h[1].SetLineColor(2)
    h[0].Draw("HIST")
    h[1].Draw("HISTsame")
    h[0].GetYaxis().SetRangeUser(0,max(h[0].GetMaximum(), h[1].GetMaximum()) * 1.3)
    legend = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend.AddEntry(h[0],"Prompt","l")
    legend.AddEntry(h[1],"Fake","l")
    legend.Draw("same")
    if(log_flag):
        h[0].GetYaxis().SetRangeUser(1e-5,max(h[0].GetMaximum(), h[1].GetMaximum())*1.3)
        c1.SetLogy()
    c1.SaveAs("{}/{}.png".format(plot_dir, plot_name))

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
def reducedPhotonID(ggChain):
    if ( abs(ggChain.Eta) < 1.4442 and ggChain.PFChIso < 25):# and ggChain.PFNeuIso < 24.032+0.01512*ggChain.Et+0.00002259*ggChain.Et*ggChain.Et):# and ggChain.PFPhoIso < 2.876 + 0.004017*ggChain.Et): 1.694
            return True
    elif ( 1.566 < abs(ggChain.Eta) < 2.5  and ggChain.HoverE < 0.0590 and ggChain.PFChIso < 2.089 and ggChain.PFNeuIso < 19.722+0.0117*ggChain.Et+0.000023*ggChain.Et*ggChain.Et and ggChain.PFPhoISO < 4.162 + 0.0037*ggChain.Et):
            return True
    else:
            return False



if __name__ == "__main__":

    logging.basicConfig(format='%(message)s', stream=sys.stderr, level=logging.INFO)

    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('--promptFiles','-p',nargs="+",help='input prompt filenames')
    parser.add_argument('--fakeFiles','-f',nargs="+",help='input fake filenames')
    parser.add_argument('--outputDirectory','-o',default="plots/merged",help='output plots directory')
    parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
    args = parser.parse_args()

    # Read infile
    ggChains=[]
    ggChains.append(ROOT.TChain("recoPhoton"))
    ggChains.append(ROOT.TChain("recoPhoton"))
    for promptFile in args.promptFiles:
        logging.info("Adding prompt file {}".format(promptFile))
        ggChains[0].Add(promptFile)
    for fakeFile in args.fakeFiles:
        logging.info("Adding fake file {}".format(fakeFile))
        ggChains[1].Add(fakeFile)


    h_Et=[]
    h_Eta=[]
    h_Phi=[]

    h_SigmaIEtaIEta=[]
    h_SigmaIEtaIPhi=[]
    h_SigmaIPhiIPhi=[]
    h_SigmaIEtaIEta_SigmaIPhiIPhi=[]
    h_SigmaEta=[]
    h_SigmaPhi=[]
    h_SigmaEta_SigmaPhi=[]
    h_Emax_E3x3=[]
    h_Emax_ESCraw=[]
    h_E5x5_ESCraw=[]
    h_E2_ESCraw=[]
    h_E2x2_E3x3=[]
    h_E2x5_ESCraw=[]
    h_E2_Emax=[]
    h_E1x3_ESCraw=[]
    h_E2x2_ESCraw=[]
    h_HoverE=[]

    #Isolation histograms
    h_PFPhoIso=[]
    h_PFChIso=[]
    h_PFNeuIso=[]
    h_PFPhoIso_afterSelection=[]
    h_PFChIso_afterSelection=[]
    h_PFNeuIso_afterSelection=[]


    for i in range(2):
        if ( i==0 ):
            name = "prompt"
        else:
            name = "fake"
        h_Et.append(ROOT.TH1D("h_{}_Et".format(name),"",100, 0, 200))
        h_Eta.append(ROOT.TH1D("h_{}_Eta".format(name),"",60,-1.5,1.5))
        h_Phi.append(ROOT.TH1D("h_{}_Phi".format(name),"",100,-5,5))
        h_SigmaIEtaIEta.append(ROOT.TH1D("h_{}_SigmaIEtaIEta".format(name),"",100,0,0.03))
        h_SigmaIEtaIPhi.append(ROOT.TH1D("h_{}_SigmaIEtaIPhi".format(name),"",100,-0.0002,0.0002))
        h_SigmaIPhiIPhi.append(ROOT.TH1D("h_{}_SigmaIPhiIPhi".format(name),"",100,0,0.03))
        h_SigmaIEtaIEta_SigmaIPhiIPhi.append(ROOT.TH1D("h_{}_SigmaIEtaIEta_SigmaIPhiIPhi".format(name),"",100,0,1.7))
        h_SigmaEta.append(ROOT.TH1D("h_{}_SigmaEta".format(name),"",100,0,0.025))
        h_SigmaPhi.append(ROOT.TH1D("h_{}_SigmaPhi".format(name),"",100,0,0.09))
        h_SigmaEta_SigmaPhi.append(ROOT.TH1D("h_{}_SigmaEta_SigmaPhi".format(name),"",100,0,1.5))
        h_Emax_E3x3.append(ROOT.TH1D("h_{}_Emax_E3x3".format(name),"",100,0,1))
        h_Emax_ESCraw.append(ROOT.TH1D("h_{}_Emax_ESCraw".format(name),"",100,0,1))
        h_E5x5_ESCraw.append(ROOT.TH1D("h_{}_E5x5_ESCraw".format(name),"",100,0,1.2))
        h_E2_ESCraw.append(ROOT.TH1D("h_{}_E2_ESCraw".format(name),"",100,0,0.5))
        h_E2x2_E3x3.append(ROOT.TH1D("h_{}_E2x2_E3x3".format(name),"",100,0,1))
        h_E2x5_ESCraw.append(ROOT.TH1D("h_{}_E2x5_ESCraw".format(name),"",100,0,1))
        h_E2_Emax.append(ROOT.TH1D("h_{}_E2_Emax".format(name),"",100,0,1))
        h_E1x3_ESCraw.append(ROOT.TH1D("h_{}_E1x3_ESCraw".format(name),"",100,0,1))
        h_E2x2_ESCraw.append(ROOT.TH1D("h_{}_E2x2_ESCraw".format(name),"",100,0,1))
        h_HoverE.append(ROOT.TH1D("h_{}_HoverE".format(name),"",100,0,0.3))
        h_PFPhoIso.append(ROOT.TH1D("h_{}_PFPhoIso".format(name),"",100,0,10))
        h_PFChIso.append(ROOT.TH1D("h_{}_PFChIso".format(name),"",100,0,30))
        h_PFNeuIso.append(ROOT.TH1D("h_{}_PFNeuIso".format(name),"",100,0,10))
        h_PFPhoIso_afterSelection.append(ROOT.TH1D("h_{}_PFPhoIso_afterSelection".format(name),"",100,0,10))
        h_PFChIso_afterSelection.append(ROOT.TH1D("h_{}_PFChIso_afterSelection".format(name),"",100,0,30))
        h_PFNeuIso_afterSelection.append(ROOT.TH1D("h_{}_PFNeuIso_afterSelection".format(name),"",100,0,10))



    # ==================== Loop ==================== #
    totalEvent=[]
    passedEvent=[]

    for i, ggChain in enumerate(ggChains):
        totalEvent.append(ggChain.GetEntries())
        passedEvent.append(0)
        for event_nm in range(ggChain.GetEntries()):
            if (event_nm % args.report == 0):
                logging.info("Processing event {}/{} ({:.1f}%) ".format(event_nm, ggChain.GetEntries(), 100*event_nm/ggChain.GetEntries()))

            # read parameters
            ggChain.GetEntry(event_nm)


            h_PFChIso[i].Fill(ggChain.PFChIso)
            h_PFNeuIso[i].Fill(ggChain.PFNeuIso)
            h_PFPhoIso[i].Fill(ggChain.PFPhoIso)
            if (ggChain.matchGenFlag == True):
                passedEvent[i]+=1
                # fill histograms
                h_Et[i].Fill(ggChain.Et)
                h_Eta[i].Fill(ggChain.Eta)
                h_Phi[i].Fill(ggChain.Phi)
                h_SigmaIEtaIEta[i].Fill(ggChain.SigmaIEtaIEtaFull5x5)
                h_SigmaIEtaIPhi[i].Fill(ggChain.SigmaIEtaIPhiFull5x5)
                h_SigmaIPhiIPhi[i].Fill(ggChain.SigmaIPhiIPhiFull5x5)
                # h_SigmaIEtaIEta_SigmaIPhiIPhi[i].Fill(ggChain.SigmaIEtaIEtaFull5x5/ggChain.SigmaIPhiIPhiFull5x5) #FIXME divided by zero
                h_SigmaEta[i].Fill(ggChain.SCEtaWidth)
                h_SigmaPhi[i].Fill(ggChain.SCPhiWidth)
                h_SigmaEta_SigmaPhi[i].Fill(ggChain.SCEtaWidth/ggChain.SCPhiWidth)
                h_Emax_E3x3[i].Fill(ggChain.SeedEnergy/ggChain.R9Full5x5)
                h_Emax_ESCraw[i].Fill(ggChain.SeedEnergy/ggChain.SCRawE)
                h_E5x5_ESCraw[i].Fill(ggChain.E5x5Full5x5/ggChain.SCRawE)
                # h_E2_ESCraw[i].Fill(ggChain.E2_ggChain.SCRawE)
                h_E2x2_E3x3[i].Fill(ggChain.E2x2Full5x5/ggChain.R9Full5x5)
                # h_E2x5_ESCraw[i].Fill(ggChain.E2x5_ggChain.SCRawE)
                # h_E2_Emax[i].Fill(ggChain.E2_ggChain.SeedEnergy)
                # h_E1x3_ESCraw[i].Fill(ggChain.E1x3_ggChain.SCRawE)
                h_E2x2_ESCraw[i].Fill(ggChain.E2x2Full5x5/ggChain.SCRawE)
                h_HoverE[i].Fill(ggChain.HoverE)
                h_PFChIso_afterSelection[i].Fill(ggChain.PFChIso)
                h_PFNeuIso_afterSelection[i].Fill(ggChain.PFNeuIso)
                h_PFPhoIso_afterSelection[i].Fill(ggChain.PFPhoIso)


    logging.info("====================")
    logging.info("Prompt Selection efficiency = {:.4f}".format(passedEvent[0] / totalEvent[0]))
    logging.info("Fake   Selection efficiency = {:.4f}".format(passedEvent[1] / totalEvent[1]))
    logging.info("====================")

    # ========== Normalization ========== #
    factor = 100.
    for i in range(2):
        h_Et[i].Scale(factor/h_Et[i].Integral())
        h_Eta[i].Scale(factor/h_Eta[i].Integral())
        h_Phi[i].Scale(factor/h_Phi[i].Integral())
        h_SigmaIEtaIEta[i].Scale(factor/h_SigmaIEtaIEta[i].Integral())
        h_SigmaIEtaIPhi[i].Scale(factor/h_SigmaIEtaIPhi[i].Integral())
        h_SigmaIPhiIPhi[i].Scale(factor/h_SigmaIPhiIPhi[i].Integral())
        # h_SigmaIEtaIEta_SigmaIPhiIPhi[i].Scale(factor/h_SigmaIEtaIEta_SigmaIPhiIPhi[i].Integral()) #FIXME divided by zero
        h_SigmaEta[i].Scale(factor/h_SigmaEta[i].Integral())
        h_SigmaPhi[i].Scale(factor/h_SigmaPhi[i].Integral())
        h_SigmaEta_SigmaPhi[i].Scale(factor/h_SigmaEta_SigmaPhi[i].Integral())
        h_Emax_E3x3[i].Scale(factor/h_Emax_E3x3[i].Integral())
        h_Emax_ESCraw[i].Scale(factor/h_Emax_ESCraw[i].Integral())
        h_E5x5_ESCraw[i].Scale(factor/h_E5x5_ESCraw[i].Integral())
        # h_E2_ESCraw[i].Scale(factor/h_E2_ESCraw[i].Integral())
        # h_E2x2_E3x3[i].Scale(factor/h_E2x2_E3x3[i].Integral())
        # h_E2x5_ESCraw[i].Scale(factor/h_E2x5_ESCraw[i].Integral())
        # h_E2_Emax[i].Scale(factor/h_E2_Emax[i].Integral())
        # h_E1x3_ESCraw[i].Scale(factor/h_E1x3_ESCraw[i].Integral())
        h_E2x2_ESCraw[i].Scale(factor/h_E2x2_ESCraw[i].Integral())
        h_HoverE[i].Fill(ggChain.HoverE)



    # ========== Plots ========== #
    plot_dir = args.outputDirectory

    drawHistos(h_Et,"E_{T}",plot_dir,"Et",False)
    drawHistos(h_Eta,"#eta",plot_dir,"Eta",False)
    drawHistos(h_Phi,"#phi",plot_dir,"Phi",False)
    drawHistos(h_SigmaIEtaIEta,"#sigma_{i#etai#eta}",plot_dir,"SigmaIEtaIEta",False)
    drawHistos(h_SigmaIPhiIPhi,"#sigma_{i#phii#phi}",plot_dir,"SigmaIPhiIPhi",False)
    drawHistos(h_SigmaIEtaIPhi,"#sigma_{i#etai#phi}",plot_dir,"SigmaIEtaIPhi",False)
    drawHistos(h_SigmaEta,"#sigma_{#eta}",plot_dir,"SigmaEta",False)
    drawHistos(h_SigmaPhi,"#sigma_{#phi}",plot_dir,"SigmaPhi",False)
    drawHistos(h_SigmaPhi,"#sigma_{#phi}",plot_dir,"SigmaPhi",False)
    drawHistos(h_SigmaEta_SigmaPhi,"#sigma_{#eta}/#sigma_{#phi}",plot_dir,"SigmaEta_SigmaPhi",False)
    drawHistos(h_Emax_E3x3,"E_{max}/E_{3x3}",plot_dir,"Emax_E3x3",False)
    drawHistos(h_Emax_ESCraw,"E_{max}/E_{SC}^{raw}",plot_dir,"Emax_ESCraw",False)
    drawHistos(h_E5x5_ESCraw,"E_{5x5}/E_{SC}^{raw}",plot_dir,"E5x5_ESCraw",False)
    drawHistos(h_E5x5_ESCraw,"E_{5x5}/E_{SC}^{raw}",plot_dir,"E5x5_ESCraw",False)
    drawHistos(h_E2x2_E3x3,"E_{2x2}/E_{3x3}",plot_dir,"E2x2_E3x3",False)
    drawHistos(h_HoverE,"HoverE",plot_dir,"HoverE",True)
    drawHistos(h_PFPhoIso,"PFPhoIso",plot_dir,"PFPhoIso",True)
    drawHistos(h_PFChIso,"PFChIso",plot_dir,"PFChIso",True)
    drawHistos(h_PFNeuIso,"PFNeuIso",plot_dir,"PFNeuIso",True)
    drawHistos(h_PFPhoIso_afterSelection,"PFPhoIso",plot_dir,"PFPhoIso_afterSelection",True)
    drawHistos(h_PFChIso_afterSelection,"PFChIso",plot_dir,"PFChIso_afterSelection",True)
    drawHistos(h_PFNeuIso_afterSelection,"PFNeuIso",plot_dir,"PFNeuIso_afterSelection",True)
