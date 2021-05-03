#!/usr/bin/env python

from __future__ import division

import ROOT
ROOT.gROOT.SetBatch(True)
import sys
import argparse
import math
from DataFormats.FWLite import Events, Handle
import EgammaUser.EgammaDAS2020.CoreTools as CoreTools


def isAncestor(a,p) :
        if a == p :
                return True
        for i in xrange(0,p.numberOfMothers()) :
                if isAncestor(a,p.mother(i)) :
                         return True
        return False

def deltaR(p1, p2):
        return math.sqrt( (p1.phi()-p2.phi())**2 + (p1.eta()-p2.eta())**2 )

def addMass(p1, p2):
        return (p1.p4() + p2.p4()).mass()

def deltaR_GenvsReco(eta1, phi1, eta2, phi2):
        return math.sqrt( (eta1-eta2)**2 + (phi1-phi2)**2 )

def addLVs(a, b):
        """add two Lorenz-vectors. a and b should be of the same type"""
        LV = type(a)

def reducedPhotonID(ggTree, ireco):
        HoverE = ggTree.phoHoverE[ireco];
        phoPFChIso = ggTree.phoPFChIso[ireco]
        phoPFNeuIso = ggTree.phoPFNeuIso[ireco]
        phoSigmaIEtaIEtaFull5x5 = ggTree.phoSigmaIEtaIEtaFull5x5[ireco]
        pt = ggTree.phoEt[ireco]
        eta = abs(ggTree.phoEta[ireco])

        if ( eta < 1.4442 and HoverE < 0.0597 and phoPFChIso < 1.295 and phoPFNeuIso < 10.910+0.0148*pt+0.000017*pt*pt ):
                return True
        elif ( 1.566 < eta < 2.5  and HoverE < 0.0481 and phoPFChIso < 1.011 and phoPFNeuIso < 5.931+0.0163*pt+0.000014*pt*pt ):
                return True
        else:
                return False


class RecoPhotons:
        nPho = 0
        phoE = None
        phoEt = None
        phoEta = None
        phoPhi = None





if __name__ == "__main__":

        #first we load in the FWLite libaries and enable them
        CoreTools.load_fwlitelibs()

        #now read the cmd line to find our input filename
        #local file "file:<filename>
        #remote file "root://
        #note unlike CMSSW proper its not smart enough to resolve to fall over to xrootd automatically
        #so it has to be specified
        parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
        parser.add_argument('in_filenames',nargs="+",help='input filenames')
        parser.add_argument('--prefix','-p',default='',help='file prefix')
        parser.add_argument('--out','-o',default="output.root",help='output filename')
        parser.add_argument('--report','-r',default=5000,type=int,help='report every x events')
        args = parser.parse_args()

        # example input : ~/mergedPhotonID/MC_sample/2016/HAHMHToAA_AToGG_MA-15GeV_TuneCUETP8M1_PSweights_13TeV-madgraph_pythia8/3EDB22FF-3BF6-EA11-A0E5-0242AC130002.root"
        year = (args.in_filenames[0].split("MC_sample/")[1].split("/")[0])
        a_mass_nominal = float(args.in_filenames[0].split("GeV")[0].split("-")[1].replace("p","."))
        print("Processing year: {}, m(a): {}".format(year, a_mass_nominal))

        print("Reading input miniAOD")
        in_filenames_with_prefix = ['{}{}'.format(args.prefix,x) for x in args.in_filenames]
        events = Events(in_filenames_with_prefix)

        handlePruned  = Handle ("std::vector<reco::GenParticle>")
        handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
        handlePhoton  = Handle ("std::vector<pat::Photon>")
        labelPruned = ("prunedGenParticles")
        labelPacked = ("packedGenParticles")
        labelPhoton = ("slimmedPhotons")

        fFile = ROOT.TFile("~/mergedPhotonID/MC_sample/2016/ggNtuple/1GeV/ggtree_mc.root","READ")
        ggTree = fFile.Get("ggNtuplizer/EventTree")

        h_deltaR_gen = ROOT.TH1D("h_deltaR_gen","",100,0,a_mass_nominal*0.3)
        h_deltaR_pat = ROOT.TH1D("h_deltaR_pat","",100,0,a_mass_nominal*0.3)
        h_deltaR_reco = ROOT.TH1D("h_deltaR_reco","",100,0,a_mass_nominal*0.3)
        h_deltaR_closePhotons_reco_reducedPhotonID = ROOT.TH1D("h_deltaR_closePhotons_reco_reducedPhotonID","",100,0,a_mass_nominal*0.3)
        h_deltaR_closePhotons_reco_LoosePhotonID = ROOT.TH1D("h_deltaR_closePhotons_reco_LoosePhotonID","",100,0,a_mass_nominal*0.3)
        h_mass_gen = ROOT.TH1D("h_mass_gen","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_mass_pat = ROOT.TH1D("h_mass_pat","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_mass_reco = ROOT.TH1D("h_mass_reco","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_mass_reco_reducedPhotonID = ROOT.TH1D("h_mass_reco_reducedPhotonID","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_mass_reco_LoosePhotonID = ROOT.TH1D("h_mass_reco_LoosePhotonID","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_mass_pat_merged = ROOT.TH1D("h_mass_pat_merged","",100,80,140)
        h_reco_phoCount_whilePatMerged = ROOT.TH2D("h_reco_phoCount_whilePatMerged","",5,0,5,5,0,5)
        h_pt_gen = ROOT.TH1D("h_pt_gen","GenPhoton pt",100, 0, 200)
        h_pt_reco = ROOT.TH1D("h_pt_reco","RecoPhoton pt",100, 0, 200)
        h_pt_reco_reducedPhotonID = ROOT.TH1D("h_pt_reco_reducedPhotonID","RecoPhoton pt",100, 0, 200)
        h_pt_reco_LoosePhotonID = ROOT.TH1D("h_pt_reco_LoosePhotonID","RecoPhoton pt",100, 0, 200)
        h_eta_gen = ROOT.TH1D("h_eta_gen","GenPhoton eta",60,-3,3)
        h_eta_reco = ROOT.TH1D("h_eta_reco","RecoPhoton eta",60,-3,3)
        h_eta_reco_reducedPhotonID = ROOT.TH1D("h_eta_reco_reducedPhotonID","RecoPhoton eta",60,-3,3)
        h_eta_reco_LoosePhotonID = ROOT.TH1D("h_eta_reco_LoosePhotonID","RecoPhoton eta",60,-3,3)
        h_phi_gen = ROOT.TH1D("h_phi_gen","GenPhoton phi",100,-5,5)
        h_phi_reco = ROOT.TH1D("h_phi_reco","RecoPhoton phi",100,-5,5)
        h_phi_reco_reducedPhotonID = ROOT.TH1D("h_phi_reco_reducedPhotonID","RecoPhoton phi",100,-5,5)
        h_phi_reco_LoosePhotonID = ROOT.TH1D("h_phi_reco_LoosePhotonID","RecoPhoton phi",100,-5,5)
        h_photonPair_energyRatio = ROOT.TH1D("h_photonPair_energyRatio","Photon pair energy ratio",101,0,1.01)
        h_photonPair_energyRatio_reco_reducedPhotonID = ROOT.TH1D("h_photonPair_energyRatio_reco_reducedPhotonID","Photon pair energy ratio",101,0,1.01)
        h_photonPair_energyRatio_reco_LoosePhotonID = ROOT.TH1D("h_photonPair_energyRatio_reco_LoosePhotonID","Photon pair energy ratio",101,0,1.01)
        h_pt_closePhotons_reco = ROOT.TH1D("h_pt_closePhotons_reco","Close RecoPhoton pt",100, 0, 200)
        h_pt_closePhotons_reco_reducedPhotonID = ROOT.TH1D("h_pt_closePhotons_reco_reducedPhotonID","Close RecoPhoton pt",100, 0, 200)
        h_pt_closePhotons_reco_LoosePhotonID = ROOT.TH1D("h_pt_closePhotons_reco_LoosePhotonID","Close RecoPhoton pt",100, 0, 200)

        pat4pho_event_count=0
        merged_event_count=0
        # loop over events
        for event_nm, event in enumerate(events):
                print ("------------------------------------------------------------")
                print ("Event={}".format(event_nm))
                event.getByLabel (labelPacked, handlePacked)
                event.getByLabel (labelPruned, handlePruned)
                event.getByLabel (labelPhoton, handlePhoton)
                # get the product
                packed = handlePacked.product()
                pruned = handlePruned.product()
                photon = handlePhoton.product()

                # get the reco photons HACK
                ggTree.GetEntry(event_nm)
                recoPhotons = RecoPhotons()
                recoPhotons.nPho = ggTree.nPho
                recoPhotons.phoE = ggTree.phoE
                recoPhotons.phoEt = ggTree.phoEt
                recoPhotons.phoEta = ggTree.phoEta

                genpho = []
                patmatch1_flag=False
                patmatch2_flag=False
                patClose_flag=False


                pho_GenvsPat_matchingID = []
                pho_GenvsReco_matchingID = []
                pho_GenvsPat_matching_nm = 0
                pho_GenvsReco_matching_nm = 0

                # Get the gen photon list
                for p in pruned :
                        if (p.pdgId() == 22 ):# and p.mother(0).pdgId()== 35 ):
                                print ("PdgId : {}   pt : {}  eta : {}   phi : {}  mother : {}" .format(p.pdgId(),p.pt(),p.eta(),p.phi(),p.mother(0).pdgId()))
                                genpho.append(p)
                                if (p.pt()>10):
                                        h_pt_gen.Fill(p.pt())
                                        h_eta_gen.Fill(p.eta())
                                        h_phi_gen.Fill(p.phi())

                # Fill kenemetic histograms of each photon pair ( total 2 photon pairs in h->aa->gggg )
                # for i in range(0, 2):
                #         a_mass = ( genpho[i*2].p4() + genpho[i*2+1].p4() ).mass()
                #         print("mass : {}  dR : {}".format(a_mass, deltaR(genpho[i*2], genpho[i*2+1])))
                #         h_deltaR_gen.Fill(deltaR(genpho[i*2], genpho[i*2+1]))
                #         h_mass_gen.Fill(a_mass)

                print ("    pat:photon")
                # Gen vs Pat matching
                for ipho, pho in enumerate(packed):
                        pho_GenvsPat_matchingID.append(-1)
                        if (pho.pdgId() == 22):
                                print ("    PdgId : {}   pt : {}  eta : {}   phi : {}  mass : {}" .format(pho.pdgId(),pho.pt(),pho.eta(),pho.phi(),pho.p4().mass()))
                                for igen, gen in enumerate(genpho):
                                        # pat photon matching to gen photon
                                        if ( deltaR(pho, gen) < 0.15 ):
                                                pho_GenvsPat_matchingID[ipho] = igen
                                                pho_GenvsPat_matching_nm+=1
                                                break

                # print("    {} , {}".format(pho_GenvsPat_matching_nm,pho_GenvsPat_matchingID))

                # for i, (ipho, ipho_GenvsPat_matchingID) in enumerate(zip(photon, pho_GenvsPat_matchingID)):
                #         for j, (jpho, jpho_GenvsPat_matchingID) in enumerate(zip(photon, pho_GenvsPat_matchingID)):
                #                 if ( ipho_GenvsPat_matchingID != -1 and jpho_GenvsPat_matchingID != -1 and i > j):
                #                         mass_photonPair = (ipho.p4(0)+jpho.p4(0)).mass()
                #                         delta_photonPair = deltaR(ipho, jpho)

                #                         # Matched photon pair
                #                         if (ipho_GenvsPat_matchingID == 0 and jpho_GenvsPat_matchingID == 1):
                #                                 print("    mass1 = {}, deltaR1 = {}".format(mass_photonPair , delta_photonPair))
                #                                 h_mass_pat.Fill(mass_photonPair)
                #                                 h_deltaR_pat.Fill(delta_photonPair)
                #                                 patmatch1_flag = True
                #                         elif (ipho_GenvsPat_matchingID == 2 and jpho_GenvsPat_matchingID == 3):
                #                                 print("    mass2 = {}, deltaR2 = {}".format(mass_photonPair, delta_photonPair))
                #                                 h_deltaR_pat.Fill(delta_photonPair)
                #                                 h_mass_pat.Fill(mass_photonPair)
                #                                 patmatch2_flag = True

                #                         # Very close photon pair ( still resolved )
                #                         elif (ipho_GenvsPat_matchingID == jpho_GenvsPat_matchingID):
                #                                 print("    very close photon pair!!!! mass = {}, deltaR = {}".format(mass_photonPair, delta_photonPair))
                #                                 h_deltaR_pat.Fill(delta_photonPair)
                #                                 h_mass_pat.Fill(mass_photonPair)
                #                                 patClose_flag = True


                #                         # two merged photon pairs
                #                         elif ( (ipho_GenvsPat_matchingID == 0 and jpho_GenvsPat_matchingID == 2) or (ipho_GenvsPat_matchingID == 2 or jpho_GenvsPat_matchingID == 0)):
                #                                 if (pho_GenvsPat_matching_nm == 2):
                #                                         print("    merged photon pair!!! mass = {}, deltaR = {}".format(mass_photonPair, delta_photonPair))
                #                                         h_mass_pat_merged.Fill(mass_photonPair)

        #         print ("        reco:photon")
        #         # Gen vs Reco matching
        #         for ireco in range(ggTree.nPho):
        #                 pho_GenvsReco_matchingID.append(-1) # Add index to the matching array
        #                 # Photon ID selection
        #                 print ("        pt : {}  eta : {}   phi : {}" .format(ggTree.phoEt[ireco],ggTree.phoEta[ireco],ggTree.phoPhi[ireco]))
        #                 # Gen vs Reco matching
        #                 for igen, gen in enumerate(genpho):
        #                         if ( deltaR_GenvsReco(gen.eta(), gen.phi(), ggTree.phoEta[ireco], ggTree.phoPhi[ireco]) < 0.15 ):
        #                                 pho_GenvsReco_matchingID[ireco] = igen
        #                                 pho_GenvsReco_matching_nm+=1
        #                                 h_pt_reco.Fill(ggTree.phoEt[ireco])
        #                                 h_eta_reco.Fill(ggTree.phoEta[ireco])
        #                                 h_phi_reco.Fill(ggTree.phoEta[ireco])
        #                                 if ( ggTree.phoIDbit[ireco]>>0&1 == 1 ):
        #                                         h_pt_reco_LoosePhotonID.Fill(ggTree.phoEt[ireco])
        #                                         h_eta_reco_LoosePhotonID.Fill(ggTree.phoEta[ireco])
        #                                         h_phi_reco_LoosePhotonID.Fill(ggTree.phoEta[ireco])
        #                                 if ( reducedPhotonID(ggTree, ireco) ):
        #                                         h_pt_reco_reducedPhotonID.Fill(ggTree.phoEt[ireco])
        #                                         h_eta_reco_reducedPhotonID.Fill(ggTree.phoEta[ireco])
        #                                         h_phi_reco_reducedPhotonID.Fill(ggTree.phoEta[ireco])
        #                                 break
        #         for pt, eta, phi, pid, status in zip(ggTree.mcPt, ggTree.mcEta, ggTree.mcPhi, ggTree.mcPID, ggTree.mcStatus):
        #                 print(pt, eta, phi, pid, status)
        #                 # print ( "    pt : {}  eta : {}  phi : {}  pdgID : {}  status: {}".format(ggTree.mcPt[ireco], ggTree.mcEta[ireco], ggTree.mcPhi[ireco]), ggTree.mcPID[ireco], ggTree.mcStatus[ireco] )
        #         # Reco photon pair
        #         for i, ipho_GenvsReco_matchingID in enumerate(pho_GenvsReco_matchingID):
        #                 for j, jpho_GenvsReco_matchingID in enumerate(pho_GenvsReco_matchingID):
        #                         if ( ipho_GenvsReco_matchingID != -1 and jpho_GenvsReco_matchingID != -1 and i > j ):
        #                                 p4_pho1 = ROOT.TLorentzVector()
        #                                 p4_pho2 = ROOT.TLorentzVector()
        #                                 p4_pho1.SetPtEtaPhiE( ggTree.phoEt[i], ggTree.phoEta[i], ggTree.phoPhi[i], ggTree.phoE[i])
        #                                 p4_pho2.SetPtEtaPhiE( ggTree.phoEt[j], ggTree.phoEta[j], ggTree.phoPhi[j], ggTree.phoE[j])
        #                                 mass_photonPair = (p4_pho1 + p4_pho2).M()
        #                                 delta_photonPair = math.sqrt(((ggTree.phoEta[i] - ggTree.phoEta[j])**2 + (ggTree.phoPhi[i]-ggTree.phoPhi[j])**2))
        #                                 # Matched photon pair
        #                                 if (ipho_GenvsReco_matchingID == 0 and jpho_GenvsReco_matchingID == 1):
        #                                         print("        mass1 = {}, deltaR1 = {}".format(mass_photonPair, delta_photonPair))
        #                                         h_mass_reco.Fill(mass_photonPair)
        #                                         h_deltaR_reco.Fill(delta_photonPair)
        #                                 elif (ipho_GenvsReco_matchingID == 2 and jpho_GenvsReco_matchingID == 3):
        #                                         print("        mass1 = {}, deltaR1 = {}".format(mass_photonPair, delta_photonPair))
        #                                         h_mass_reco.Fill(mass_photonPair)
        #                                         h_deltaR_reco.Fill(delta_photonPair)
        #                                 elif (ipho_GenvsReco_matchingID == jpho_GenvsReco_matchingID):
        #                                         print("        very close reco photon pair!!!! mass = {}, deltaR = {}".format(mass_photonPair, delta_photonPair))
        #                                         if ( ggTree.phoEt[i] < ggTree.phoEt[j] ):
        #                                                 photonPair_energyRatio = ggTree.phoEt[i]/ggTree.phoEt[j]
        #                                         else:
        #                                                 photonPair_energyRatio = ggTree.phoEt[j]/ggTree.phoEt[i]
        #                                         h_mass_reco.Fill(mass_photonPair)
        #                                         h_deltaR_reco.Fill(delta_photonPair)
        #                                         h_photonPair_energyRatio.Fill(photonPair_energyRatio)
        #                                         h_pt_closePhotons_reco.Fill(ggTree.phoEt[i])
        #                                         h_pt_closePhotons_reco.Fill(ggTree.phoEt[j])

        #                                 if ( ggTree.phoIDbit[i]>>0&1 == 1 and ggTree.phoIDbit[j]>>0&1 == 1):
        #                                         if (ipho_GenvsReco_matchingID == jpho_GenvsReco_matchingID):
        #                                                 if ( ggTree.phoEt[i] < ggTree.phoEt[j] ):
        #                                                         photonPair_energyRatio = ggTree.phoEt[i]/ggTree.phoEt[j]
        #                                                 else:
        #                                                         photonPair_energyRatio = ggTree.phoEt[j]/ggTree.phoEt[i]
        #                                                 h_mass_reco_LoosePhotonID.Fill(mass_photonPair)
        #                                                 h_deltaR_closePhotons_reco_LoosePhotonID.Fill(delta_photonPair)
        #                                                 h_photonPair_energyRatio_reco_LoosePhotonID.Fill(photonPair_energyRatio)
        #                                                 h_pt_closePhotons_reco_LoosePhotonID.Fill(ggTree.phoEt[i])
        #                                                 h_pt_closePhotons_reco_LoosePhotonID.Fill(ggTree.phoEt[j])

        #                                 if (reducedPhotonID(ggTree, i) and reducedPhotonID(ggTree, j)): # and ggTree.reducedPhotonIDbit[ireco]>>0&1 == 1 ):
        #                                         # Matched photon pair
        #                                         if (ipho_GenvsReco_matchingID == jpho_GenvsReco_matchingID):
        #                                                 if ( ggTree.phoEt[i] < ggTree.phoEt[j] ):
        #                                                         photonPair_energyRatio = ggTree.phoEt[i]/ggTree.phoEt[j]
        #                                                 else:
        #                                                         photonPair_energyRatio = ggTree.phoEt[j]/ggTree.phoEt[i]
        #                                                 h_mass_reco_reducedPhotonID.Fill(mass_photonPair)
        #                                                 h_deltaR_closePhotons_reco_reducedPhotonID.Fill(delta_photonPair)
        #                                                 h_photonPair_energyRatio_reco_reducedPhotonID.Fill(photonPair_energyRatio)
        #                                                 h_pt_closePhotons_reco_reducedPhotonID.Fill(ggTree.phoEt[i])
        #                                                 h_pt_closePhotons_reco_reducedPhotonID.Fill(ggTree.phoEt[j])

        # print ("Matching efficiency              : {}".format( h_pt_reco.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Matching + reduced Photon ID efficiency  : {}".format( h_pt_reco_reducedPhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Matching + Loose Photon ID efficiency    : {}".format( h_pt_reco_LoosePhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Close Matched Photons efficiency : {}".format( h_pt_closePhotons_reco.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Close Matched Photons + reduced Photon ID efficiency: {}".format( h_pt_closePhotons_reco_reducedPhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Close Matched Photons + Loose Photon ID efficiency  : {}".format( h_pt_closePhotons_reco_LoosePhotonID.GetEntries() / h_pt_gen.GetEntries() ))

        # c1 = ROOT.TCanvas()
        # h_deltaR_gen.SetTitle('Gen m(a)={}GeV'.format(a_mass_nominal))
        # h_deltaR_gen.GetXaxis().SetTitle('dR')
        # h_deltaR_gen.Draw()
        # c1.SaveAs("plots/Ma{}GeV/deltaR_gen.png".format(a_mass_nominal))

        # # h_deltaR_pat.SetTitle('Pat m(a)={}GeV'.format(a_mass_nominal))
        # # h_deltaR_pat.GetXaxis().SetTitle('dR')
        # # h_deltaR_pat.Draw()
        # # c1.Update()
        # # c1.SaveAs("plots/Ma{}GeV/deltaR_pat.png".format(a_mass_nominal))

        # h_deltaR_reco.SetTitle('Reco m(a)={}GeV'.format(a_mass_nominal))
        # h_deltaR_reco.GetXaxis().SetTitle('dR')
        # h_deltaR_reco.Draw()
        # h_deltaR_closePhotons_reco_reducedPhotonID.SetLineColor(2)
        # h_deltaR_closePhotons_reco_reducedPhotonID.Draw("same")
        # h_deltaR_closePhotons_reco_LoosePhotonID.SetLineColor(4)
        # h_deltaR_closePhotons_reco_LoosePhotonID.Draw("same")
        # legend_deltaR = ROOT.TLegend_DeltaR(0.55,0.55,0.9,0.9)
        # legend_deltaR.AddEntry(h_deltaR_reco,"Total reco")
        # legend_deltaR.AddEntry(h_deltaR_reco_reducedPhotonID,"reduced PhoID")
        # legend_deltaR.AddEntry(h_deltaR_reco_LoosePhotonID,"Loose PhoID")
        # legend_deltaR.AddEntry(h_deltaR_closePhotons_reco,"Close Photon no PhoID")
        # legend_deltaR.AddEntry(h_deltaR_closePhotons_reco_reducedPhotonID,"Close Photon reduced PhoID")
        # legend_deltaR.AddEntry(h_deltaR_closePhotons_reco_LoosePhotonID,"Close Photon Loose PhoID")

        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/deltaR_reco.png".format(a_mass_nominal))

        # h_mass_gen.SetTitle('Gen m(a)={}GeV'.format(a_mass_nominal))
        # h_mass_gen.GetXaxis().SetTitle('mass')
        # h_mass_gen.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/mass_gen.png".format(a_mass_nominal))

        # # h_mass_pat.SetTitle('Pat m(a)={}GeV'.format(a_mass_nominal))
        # # h_mass_pat.GetXaxis().SetTitle('mass')
        # # h_mass_pat.Draw()
        # # c1.Update()
        # # c1.SaveAs("plots/Ma{}GeV/mass_pat.png".format(a_mass_nominal))

        # h_mass_reco.SetTitle('Reco m(a)={}GeV'.format(a_mass_nominal))
        # h_mass_reco.GetXaxis().SetTitle('mass')
        # h_mass_reco.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/mass_reco.png".format(a_mass_nominal))

        # # h_mass_pat_merged.SetTitle('Merged Pat m(a)={}GeV'.format(a_mass_nominal))
        # # h_mass_pat_merged.GetXaxis().SetTitle('aa mass')
        # # h_mass_pat_merged.Draw()
        # # c1.Update()
        # # c1.SaveAs("plots/Ma{}GeV/mass_pat_merged.png".format(a_mass_nominal))

        # h_pt_gen.SetTitle('GenPhoton pt')
        # h_pt_gen.GetXaxis().SetTitle('pt')
        # h_pt_gen.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/pt_gen.png".format(a_mass_nominal))

        # h_eta_gen.SetTitle('Eta')
        # h_eta_gen.GetXaxis().SetTitle('eta')
        # h_eta_gen.Draw()
        # h_eta_reco.SetLineColor(2)
        # h_eta_reco.Draw("same")
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/eta.png".format(a_mass_nominal))


        # h_phi_gen.SetTitle('GenPhoton phi')
        # h_phi_gen.GetXaxis().SetTitle('phi')
        # h_phi_gen.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/phi_gen.png".format(a_mass_nominal))

        # h_pt_reco.SetTitle('RecoPhoton pt')
        # h_pt_reco.GetXaxis().SetTitle('pt')
        # h_pt_reco.SetStats(0)
        # h_pt_reco.Draw()
        # h_pt_reco_reducedPhotonID.SetLineColor(2)
        # h_pt_reco_reducedPhotonID.Draw("same")
        # h_pt_reco_LoosePhotonID.SetLineColor(4)
        # h_pt_reco_LoosePhotonID.Draw("same")
        # h_pt_closePhotons_reco.SetLineColor(6)
        # h_pt_closePhotons_reco.Draw("same")
        # h_pt_closePhotons_reco_reducedPhotonID.SetLineColor(8)
        # h_pt_closePhotons_reco_reducedPhotonID.Draw("same")
        # h_pt_closePhotons_reco_LoosePhotonID.SetLineColor(28)
        # h_pt_closePhotons_reco_LoosePhotonID.Draw("same")
        # legend_pt = ROOT.TLegend_Pt(0.55,0.55,0.9,0.9)
        # legend_pt.AddEntry(h_pt_reco,"Total reco")
        # legend_pt.AddEntry(h_pt_reco_reducedPhotonID,"reduced PhoID")
        # legend_pt.AddEntry(h_pt_reco_LoosePhotonID,"Loose PhoID")
        # legend_pt.AddEntry(h_pt_closePhotons_reco,"Close Photon no PhoID")
        # legend_pt.AddEntry(h_pt_closePhotons_reco_reducedPhotonID,"Close Photon reduced PhoID")
        # legend_pt.AddEntry(h_pt_closePhotons_reco_LoosePhotonID,"Close Photon Loose PhoID")
        # # legend_pt.SetTextSize(0.08)
        # legend_pt.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/pt.png".format(a_mass_nominal))


        # h_phi_reco.SetTitle('RecoPhoton phi')
        # h_phi_reco.GetXaxis().SetTitle('phi')
        # h_phi_reco.Draw()
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/phi_reco.png".format(a_mass_nominal))

        # h_photonPair_energyRatio.SetTitle('Photon pair energy ratio')
        # h_photonPair_energyRatio.GetXaxis().SetTitle('E_photon1/E_photon2')
        # h_photonPair_energyRatio.Draw()
        # h_photonPair_energyRatio_reco_reducedPhotonID.SetLineColor(2)
        # h_photonPair_energyRatio_reco_reducedPhotonID.Draw("same")
        # c1.Update()
        # c1.SaveAs("plots/Ma{}GeV/photon_pair_energy_ratio.png".format(a_mass_nominal))
