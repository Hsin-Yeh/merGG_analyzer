#!/usr/bin/env python

import ROOT
import os, sys, logging
ROOT.gROOT.SetBatch(True)
import argparse
import math
from array import array

if __name__ == "__main__":

    logging.basicConfig(format='%(message)s', stream=sys.stderr, level=logging.INFO)

    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('--promptFiles','-p',nargs="+",help='input prompt filenames')
    parser.add_argument('--fakeFiles','-f',nargs="+",help='input fake filenames')
    parser.add_argument('--outputDirectory','-o',default="plots/merged",help='output plots directory')
    parser.add_argument('--report','-r',default=5000,type=int,help='report every x events')
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
        h_HoverE.append(ROOT.TH1D("h_{}_HoverE".format(name),"",100,0,0.08))


    # ==================== Loop ==================== #
    for i, ggChain in enumerate(ggChains):
        for event_nm in range(ggChain.GetEntries()):
            if (event_nm % args.report == 0):
                logging.info("Processing event {} ".format(event_nm))

            # read parameters
            ggChain.GetEntry(event_nm)

            # fill histograms
            for j in range(2):
                h_Et[i].Fill(ggChain.Et[j])
                h_Eta[i].Fill(ggChain.Eta[j])
                h_Phi[i].Fill(ggChain.Phi[j])
                h_SigmaIEtaIEta[i].Fill(ggChain.SigmaIEtaIEtaFull5x5[j])
                h_SigmaIEtaIPhi[i].Fill(ggChain.SigmaIEtaIPhiFull5x5[j])
                h_SigmaIPhiIPhi[i].Fill(ggChain.SigmaIPhiIPhiFull5x5[j])
                # h_SigmaIEtaIEta_SigmaIPhiIPhi[i].Fill(ggChain.SigmaIEtaIEtaFull5x5/ggChain.SigmaIPhiIPhiFull5x5) #FIXME divided by zero
                h_SigmaEta[i].Fill(ggChain.SCEtaWidth[j])
                h_SigmaPhi[i].Fill(ggChain.SCPhiWidth[j])
                h_SigmaEta_SigmaPhi[i].Fill(ggChain.SCEtaWidth[j]/ggChain.SCPhiWidth[j])
                h_Emax_E3x3[i].Fill(ggChain.SeedEnergy[j]/ggChain.R9Full5x5[j])
                h_Emax_ESCraw[i].Fill(ggChain.SeedEnergy[j]/ggChain.SCRawE[j])
                h_E5x5_ESCraw[i].Fill(ggChain.E5x5Full5x5[j]/ggChain.SCRawE[j])
                # h_E2_ESCraw[i].Fill(ggChain.E2_ggChain.SCRawE)
                h_E2x2_E3x3[i].Fill(ggChain.E2x2Full5x5[j]/ggChain.R9Full5x5[j])
                # h_E2x5_ESCraw[i].Fill(ggChain.E2x5_ggChain.SCRawE)
                # h_E2_Emax[i].Fill(ggChain.E2_ggChain.SeedEnergy)
                # h_E1x3_ESCraw[i].Fill(ggChain.E1x3_ggChain.SCRawE)
                h_E2x2_ESCraw[i].Fill(ggChain.E2x2Full5x5[j]/ggChain.SCRawE[j])
                h_HoverE[i].Fill(ggChain.HoverE[j])


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
        h_E2x2_E3x3[i].Scale(factor/h_E2x2_E3x3[i].Integral())
        # h_E2x5_ESCraw[i].Scale(factor/h_E2x5_ESCraw[i].Integral())
        # h_E2_Emax[i].Scale(factor/h_E2_Emax[i].Integral())
        # h_E1x3_ESCraw[i].Scale(factor/h_E1x3_ESCraw[i].Integral())
        h_E2x2_ESCraw[i].Scale(factor/h_E2x2_ESCraw[i].Integral())
        h_HoverE[i].Fill(ggChain.HoverE)



    # ========== Plots ========== #
    c1 = ROOT.TCanvas()
    c1.SetGrid()
    plot_dir = args.outputDirectory

    h_Et[1].GetYaxis().SetTitle("a.u.")
    h_Et[1].GetXaxis().SetTitle("E_{T}")
    h_Et[0].SetLineColor(209)
    h_Et[1].SetLineColor(2)
    h_Et[1].Draw("HIST")
    h_Et[0].Draw("HISTsame")
    h_Et[1].GetYaxis().SetRangeUser(0,h_Et[1].GetMaximum()*1.3)
    legend_Et = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_Et.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_Et.AddEntry(h_Et[0],"Prompt","l")
    legend_Et.AddEntry(h_Et[1],"Fake","l")
    legend_Et.Draw("same")
    c1.Update()
    c1.SaveAs("{}/Et.png".format(plot_dir))

    h_Eta[1].GetYaxis().SetTitle("a.u.")
    h_Eta[1].GetXaxis().SetTitle("#eta")
    h_Eta[0].SetLineColor(209)
    h_Eta[1].SetLineColor(2)
    h_Eta[1].Draw("HIST")
    h_Eta[0].Draw("HISTsame")
    h_Eta[1].GetYaxis().SetRangeUser(0,h_Eta[1].GetMaximum()*1.3)
    legend_Eta = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_Eta.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_Eta.AddEntry(h_Eta[0],"Prompt","l")
    legend_Eta.AddEntry(h_Eta[1],"Fake","l")
    legend_Eta.Draw("same")
    c1.Update()
    c1.SaveAs("{}/Eta.png".format(plot_dir))

    h_Phi[1].GetYaxis().SetTitle("a.u.")
    h_Phi[1].GetXaxis().SetTitle("#phi")
    h_Phi[0].SetLineColor(209)
    h_Phi[1].SetLineColor(2)
    h_Phi[1].Draw("HIST")
    h_Phi[0].Draw("HISTsame")
    h_Phi[1].GetYaxis().SetRangeUser(0,h_Phi[1].GetMaximum()*1.3)
    legend_Phi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_Phi.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_Phi.AddEntry(h_Phi[0],"Prompt","l")
    legend_Phi.AddEntry(h_Phi[1],"Fake","l")
    legend_Phi.Draw("same")
    c1.Update()
    c1.SaveAs("{}/Phi.png".format(plot_dir))

    h_SigmaIEtaIEta[0].GetYaxis().SetTitle("a.u.")
    h_SigmaIEtaIEta[0].GetXaxis().SetTitle("#sigma_{i#etai#eta}")
    h_SigmaIEtaIEta[0].SetLineColor(209)
    h_SigmaIEtaIEta[1].SetLineColor(2)
    h_SigmaIEtaIEta[0].Draw("HIST")
    h_SigmaIEtaIEta[1].Draw("HISTsame")
    h_SigmaIEtaIEta[0].GetYaxis().SetRangeUser(0,h_SigmaIEtaIEta[0].GetMaximum()*1.3)
    legend_SigmaIEtaIEta = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaIEtaIEta.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaIEtaIEta.AddEntry(h_SigmaIEtaIEta[0],"Prompt","l")
    legend_SigmaIEtaIEta.AddEntry(h_SigmaIEtaIEta[1],"Fake","l")
    legend_SigmaIEtaIEta.Draw("same")
    c1.Update()
    c1.SaveAs("{}/SigmaIEtaIEta.png".format(plot_dir))

    h_SigmaIEtaIPhi[0].GetYaxis().SetTitle("a.u.")
    h_SigmaIEtaIPhi[0].GetXaxis().SetTitle("#sigma_{i#etai#phi}")
    h_SigmaIEtaIPhi[0].SetLineColor(209)
    h_SigmaIEtaIPhi[1].SetLineColor(2)
    h_SigmaIEtaIPhi[0].Draw("HIST")
    h_SigmaIEtaIPhi[1].Draw("HISTsame")
    h_SigmaIEtaIPhi[0].GetYaxis().SetRangeUser(0,h_SigmaIEtaIPhi[0].GetMaximum()*1.3)
    legend_SigmaIEtaIPhi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaIEtaIPhi.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaIEtaIPhi.AddEntry(h_SigmaIEtaIPhi[0],"Prompt","l")
    legend_SigmaIEtaIPhi.AddEntry(h_SigmaIEtaIPhi[1],"Fake","l")
    legend_SigmaIEtaIPhi.Draw("same")
    c1.Update()
    c1.SaveAs("{}/SigmaIEtaIPhi.png".format(plot_dir))

    h_SigmaIPhiIPhi[0].GetYaxis().SetTitle("a.u.")
    h_SigmaIPhiIPhi[0].GetXaxis().SetTitle("#sigma_{i#phii#phi}")
    h_SigmaIPhiIPhi[0].SetLineColor(209)
    h_SigmaIPhiIPhi[1].SetLineColor(2)
    h_SigmaIPhiIPhi[0].Draw("HIST")
    h_SigmaIPhiIPhi[1].Draw("HISTsame")
    h_SigmaIPhiIPhi[0].GetYaxis().SetRangeUser(0,h_SigmaIPhiIPhi[0].GetMaximum()*1.3)
    legend_SigmaIPhiIPhi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaIPhiIPhi.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaIPhiIPhi.AddEntry(h_SigmaIPhiIPhi[0],"Prompt","l")
    legend_SigmaIPhiIPhi.AddEntry(h_SigmaIPhiIPhi[1],"Fake","l")
    legend_SigmaIPhiIPhi.Draw("same")
    c1.Update()
    c1.SaveAs("{}/SigmaIPhiIPhi.png".format(plot_dir))

    # h_SigmaIEtaIEta_SigmaIPhiIPhi[1].GetYaxis().SetTitle("a.u.")
    # h_SigmaIEtaIEta_SigmaIPhiIPhi[1].GetXaxis().SetTitle("# SigmaIEtaIEta_SigmaIPhiIPhi")
    #     # h_SigmaIEtaIEta_SigmaIPhiIPhi[0].SetLineColor(2)
    # h_SigmaIEtaIEta_SigmaIPhiIPhi[1].SetLineColor(2)
    # h_SigmaIEtaIEta_SigmaIPhiIPhi[1].Draw("HIST")
    # h_SigmaIEtaIEta_SigmaIPhiIPhi[0].Draw("HISTsame")
    # legend_SigmaIEtaIEta_SigmaIPhiIPhi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_SigmaIEtaIEta_SigmaIPhiIPhi.AddEntry(h_SigmaIEtaIEta_SigmaIPhiIPhi[0],"Prompt","l")
    # legend_SigmaIEtaIEta_SigmaIPhiIPhi.AddEntry(h_SigmaIEtaIEta_SigmaIPhiIPhi[1],"Fake","l")
    # legend_SigmaIEtaIEta_SigmaIPhiIPhi.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/SigmaIEtaIEta_SigmaIPhiIPhi.png".format(plot_dir))

    h_SigmaEta[0].GetYaxis().SetTitle("a.u.")
    h_SigmaEta[0].GetXaxis().SetTitle("#sigma_{#eta}")
    h_SigmaEta[0].SetLineColor(209)
    h_SigmaEta[1].SetLineColor(2)
    h_SigmaEta[0].Draw("HIST")
    h_SigmaEta[1].Draw("HISTsame")
    h_SigmaEta[0].GetYaxis().SetRangeUser(0,h_SigmaEta[0].GetMaximum()*1.3)
    legend_SigmaEta = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaEta.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaEta.AddEntry(h_SigmaEta[0],"Prompt","l")
    legend_SigmaEta.AddEntry(h_SigmaEta[1],"Fake","l")
    legend_SigmaEta.Draw("same")
    c1.Update()
    c1.SaveAs("{}/SigmaEta.png".format(plot_dir))

    h_SigmaPhi[0].GetYaxis().SetTitle("a.u.")
    h_SigmaPhi[0].GetXaxis().SetTitle("#sigma_{#phi}")
    h_SigmaPhi[0].SetLineColor(209)
    h_SigmaPhi[1].SetLineColor(2)
    h_SigmaPhi[0].Draw("HIST")
    h_SigmaPhi[1].Draw("HISTsame")
    h_SigmaPhi[0].GetYaxis().SetRangeUser(0,h_SigmaPhi[0].GetMaximum()*1.3)
    legend_SigmaPhi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaPhi.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaPhi.AddEntry(h_SigmaPhi[0],"Prompt","l")
    legend_SigmaPhi.AddEntry(h_SigmaPhi[1],"Fake","l")
    legend_SigmaPhi.Draw("same")
    c1.Update()
    c1.SaveAs("{}/SigmaPhi.png".format(plot_dir))

    h_SigmaEta_SigmaPhi[1].GetYaxis().SetTitle("a.u.")
    h_SigmaEta_SigmaPhi[1].GetXaxis().SetTitle("#sigma_{#eta}/#sigma_{#phi}")
    h_SigmaEta_SigmaPhi[0].SetLineColor(209)
    h_SigmaEta_SigmaPhi[1].SetLineColor(2)
    h_SigmaEta_SigmaPhi[1].Draw("HIST")
    h_SigmaEta_SigmaPhi[0].Draw("HISTsame")
    h_SigmaEta_SigmaPhi[1].GetYaxis().SetRangeUser(0,h_SigmaEta_SigmaPhi[1].GetMaximum()*1.3)
    legend_SigmaEta_SigmaPhi = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_SigmaEta_SigmaPhi.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_SigmaEta_SigmaPhi.AddEntry(h_SigmaEta_SigmaPhi[0],"Prompt","l")
    legend_SigmaEta_SigmaPhi.AddEntry(h_SigmaEta_SigmaPhi[1],"Fake","l")
    legend_SigmaEta_SigmaPhi.Draw("same")
    c1.Update()
    c1.SaveAs("{}/h_SigmaEta_SigmaPhi.png".format(plot_dir))

    h_Emax_E3x3[1].GetYaxis().SetTitle("a.u.")
    h_Emax_E3x3[1].GetXaxis().SetTitle("E_{max}/E_{3x3}")
    h_Emax_E3x3[0].SetLineColor(209)
    h_Emax_E3x3[1].SetLineColor(2)
    h_Emax_E3x3[1].Draw("HIST")
    h_Emax_E3x3[0].Draw("HISTsame")
    h_Emax_E3x3[1].GetYaxis().SetRangeUser(0,h_Emax_E3x3[1].GetMaximum()*1.3)
    legend_Emax_E3x3 = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_Emax_E3x3.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_Emax_E3x3.AddEntry(h_Emax_E3x3[0],"Prompt","l")
    legend_Emax_E3x3.AddEntry(h_Emax_E3x3[1],"Fake","l")
    legend_Emax_E3x3.Draw("same")
    c1.Update()
    c1.SaveAs("{}/Emax_E3x3.png".format(plot_dir))

    h_Emax_ESCraw[1].GetYaxis().SetTitle("a.u.")
    h_Emax_ESCraw[1].GetXaxis().SetTitle("E_{max}/E_{SC}^{raw}")
    h_Emax_ESCraw[0].SetLineColor(209)
    h_Emax_ESCraw[1].SetLineColor(2)
    h_Emax_ESCraw[1].Draw("HIST")
    h_Emax_ESCraw[0].Draw("HISTsame")
    h_Emax_ESCraw[1].GetYaxis().SetRangeUser(0,h_Emax_ESCraw[1].GetMaximum()*1.3)
    legend_Emax_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_Emax_ESCraw.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_Emax_ESCraw.AddEntry(h_Emax_ESCraw[0],"Prompt","l")
    legend_Emax_ESCraw.AddEntry(h_Emax_ESCraw[1],"Fake","l")
    legend_Emax_ESCraw.Draw("same")
    c1.Update()
    c1.SaveAs("{}/Emax_ESCraw.png".format(plot_dir))

    h_E5x5_ESCraw[0].GetYaxis().SetTitle("a.u.")
    h_E5x5_ESCraw[0].GetXaxis().SetTitle("E_{5x5}/E_{SC}^{raw}")
    h_E5x5_ESCraw[0].SetLineColor(209)
    h_E5x5_ESCraw[1].SetLineColor(2)
    h_E5x5_ESCraw[0].Draw("HIST")
    h_E5x5_ESCraw[1].Draw("HISTsame")
    h_E5x5_ESCraw[0].GetYaxis().SetRangeUser(0,h_E5x5_ESCraw[0].GetMaximum()*1.3)
    legend_E5x5_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_E5x5_ESCraw.SetHeader("0<|#eta|<1.4442  E_{T}>10")
    legend_E5x5_ESCraw.AddEntry(h_E5x5_ESCraw[0],"Prompt","l")
    legend_E5x5_ESCraw.AddEntry(h_E5x5_ESCraw[1],"Fake","l")
    legend_E5x5_ESCraw.Draw("same")
    c1.Update()
    c1.SaveAs("{}/E5x5_ESCraw.png".format(plot_dir))

    # h_E2_ESCraw[1].GetYaxis().SetTitle("a.u.")
    # h_E2_ESCraw[1].GetXaxis().SetTitle("E2_ESCraw")
    # h_E2_ESCraw[1].SetLineColor(2)
    # h_E2_ESCraw[1].Draw("HIST")
    # h_E2_ESCraw[0].Draw("HISTsame")
    # legend_E2_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_E2_ESCraw.AddEntry(h_E2_ESCraw[0],"Prompt","l")
    # legend_E2_ESCraw.AddEntry(h_E2_ESCraw[1],"Fake","l")
    # legend_E2_ESCraw.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/E2_ESCraw.png".format(plot_dir))

    h_E2x2_E3x3[1].GetYaxis().SetTitle("a.u.")
    h_E2x2_E3x3[1].GetXaxis().SetTitle("E_{2x2}/E_{3x3}")
    h_E2x2_E3x3[0].SetLineColor(209)
    h_E2x2_E3x3[1].SetLineColor(2)
    h_E2x2_E3x3[1].Draw("HIST")
    h_E2x2_E3x3[0].Draw("HISTsame")
    h_E2x2_E3x3[1].GetYaxis().SetRangeUser(0,h_E2x2_E3x3[1].GetMaximum()*1.3)
    legend_E2x2_E3x3 = ROOT.TLegend(0.6,0.72,0.88,0.88)
    legend_E2x2_E3x3.SetHeader("0<#eta<1.4442  E_{T}>10")
    legend_E2x2_E3x3.AddEntry(h_E2x2_E3x3[0],"Prompt","l")
    legend_E2x2_E3x3.AddEntry(h_E2x2_E3x3[1],"Fake","l")
    legend_E2x2_E3x3.Draw("same")
    c1.Update()
    c1.SaveAs("{}/E2x2_E3x3.png".format(plot_dir))

    # h_E2x5_ESCraw[1].GetYaxis().SetTitle("a.u.")
    # h_E2x5_ESCraw[1].GetXaxis().SetTitle("# E2x5_ESCraw")
    # h_E2x5_ESCraw[1].SetLineColor(2)
    # h_E2x5_ESCraw[1].Draw("HIST")
    # h_E2x5_ESCraw[0].Draw("HISTsame")
    # legend_E2x5_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_E2x5_ESCraw.AddEntry(h_E2x5_ESCraw[0],"Prompt","l")
    # legend_E2x5_ESCraw.AddEntry(h_E2x5_ESCraw[1],"Fake","l")
    # legend_E2x5_ESCraw.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/h E2x5_ESCraw.png".format(plot_dir))

    # h_E2_Emax[1].GetYaxis().SetTitle("a.u.")
    # h_E2_Emax[1].GetXaxis().SetTitle("# E2_Emax")
    # h_E2_Emax[1].SetLineColor(2)
    # h_E2_Emax[1].Draw("HIST")
    # h_E2_Emax[0].Draw("HISTsame")
    # legend_E2_Emax = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_E2_Emax.AddEntry(h_E2_Emax[0],"Prompt","l")
    # legend_E2_Emax.AddEntry(h_E2_Emax[1],"Fake","l")
    # legend_E2_Emax.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/E2_Emax.png".format(plot_dir))

    # h_E1x3_ESCraw[1].GetYaxis().SetTitle("a.u.")
    # h_E1x3_ESCraw[1].GetXaxis().SetTitle("# E1x3_ESCraw")
    # h_E1x3_ESCraw[1].SetLineColor(2)
    # h_E1x3_ESCraw[1].Draw("HIST")
    # h_E1x3_ESCraw[0].Draw("HISTsame")
    # legend_E1x3_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_E1x3_ESCraw.AddEntry(h_1x3_ESCraw[0],"Prompt","l")
    # legend_E1x3_ESCraw.AddEntry(h_1x3_ESCraw[1],"Fake","l")
    # legend_E1x3_ESCraw.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/E1x3_ESCraw.png".format(plot_dir))

    # h_E2x2_ESCraw[1].GetYaxis().SetTitle("a.u.")
    # h_E2x2_ESCraw[1].GetXaxis().SetTitle("# E2x2_ESCraw")
    # h_E2x2_ESCraw[1].SetLineColor(2)
    # h_E2x2_ESCraw[1].Draw("HIST")
    # h_E2x2_ESCraw[0].Draw("HISTsame")
    # legend_E2x2_ESCraw = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_E2x2_ESCraw.AddEntry(h_E2x2_ESCraw[0],"Prompt","l")
    # legend_E2x2_ESCraw.AddEntry(h_E2x2_ESCraw[1],"Fake","l")
    # legend_E2x2_ESCraw.Draw("same")
    # c1.Update()
    # c1.SaveAs("{}/E2x2_ESCraw.png".format(plot_dir))

    # h_HoverE[1].GetYaxis().SetTitle("a.u.")
    # h_HoverE[1].GetXaxis().SetTitle("HoverE")
    # h_HoverE[1].SetLineColor(2)
    # h_HoverE[1].Draw("HIST")
    # h_HoverE[0].Draw("HISTsame")
    # h_HoverE[1].GetYaxis().SetRangeUser(0,h_HoverE[1].GetMaximum()*1.3)
    # legend_HoverE = ROOT.TLegend(0.6,0.72,0.88,0.88)
    # legend_HoverE.AddEntry(h_HoverE[0],"Prompt","l")
    # legend_HoverE.AddEntry(h_HoverE[1],"Fake","l")
    # legend_HoverE.Draw("same")
    # c1.Update()
    # c1.SetLogy()
    # c1.SaveAs("{}/HoverE.png".format(plot_dir))
