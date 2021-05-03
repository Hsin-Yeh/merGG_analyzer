#!/usr/bin/env python

import ROOT
import os, sys, logging
ROOT.gROOT.SetBatch(True)
import argparse
import math
from array import array

def color(i):
    if (i==0):
        return 1
    elif (i==1):
        return 2
    elif (i==2):
        return 4
    elif (i==3):
        return 7
    elif (i==4):
        return 8
    elif (i==5):
        return 9
    elif (i==6):
        return 50
    elif (i==7):
        return 51
    elif (i==8):
        return 28

def drawHistos(h, l, xtitle):
    for i, h in enumerate(h):
        if (i==0):
            h.GetYaxis().SetTitle("a.u.")
            h.GetXaxis().SetTitle(xtitle)
            h.Draw("HIST")
        else:
            h.Draw("HISTsame")
    l.Draw("same")



if __name__ == "__main__":

    logging.basicConfig(format='%(message)s', stream=sys.stderr, level=logging.INFO)

    parser = argparse.ArgumentParser(description='compare gen plots')
    parser.add_argument('in_filenames',nargs="+",help='input filenames')
    parser.add_argument('--report','-r',default=5000,type=int,help='report every x events')
    args = parser.parse_args()

    c1 = ROOT.TCanvas()
    h_Et=[]
    h_dR=[]
    l_Et = ROOT.TLegend(0.6,0.6,0.88,0.88)
    l_dR = ROOT.TLegend(0.6,0.6,0.88,0.88)
    inFiles=[]
    factor=100

    for i, in_filename in enumerate(args.in_filenames):
        print("Reading fle: {}".format(in_filename))
        inFiles.append(ROOT.TFile.Open(in_filename ,"READ"))
        a_mass_nominal = in_filename.split("/output")[0].split("/M")[1].replace("p",".")

        h_Et.append(inFiles[i].Get("h_pt_gen").Clone())
        h_Et[i].SetLineColor(color(i))
        h_Et[i].Scale(factor/h_Et[i].Integral())
        l_Et.AddEntry(h_Et[i],"M(a)={}GeV".format(a_mass_nominal),"l")

        h_dR.append(inFiles[i].Get("h_deltaR_gen").Clone())
        h_dR[i].Scale(factor/h_dR[i].Integral())
        h_dR[i].SetLineColor(color(i))
        l_dR.AddEntry(h_dR[i],"M(a)={}GeV".format(a_mass_nominal),"l")


    drawHistos(h_Et, l_Et, "pt")
    c1.SaveAs("genPt.png")

    drawHistos(h_dR, l_dR, "dR")
    c1.SaveAs("gendR.png")
