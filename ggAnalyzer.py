#!/usr/bin/env python

import ROOT
import os, sys, logging
ROOT.gROOT.SetBatch(True)
# ROOT.gROOT.Macro( os.path.expanduser( '~/.rootlogon.C' ))
import argparse
import math


def deltaR(pho1, pho2):
    return math.sqrt( (pho1.Eta-pho2.Eta)**2 + (pho1.Phi-pho2.Phi)**2 )

def reducedPhotonID(pho):
    if ( abs(pho.Eta) < 1.4442 and pho.HoverE < 0.0597 and pho.PFChIso < 1.295 and pho.PFNeuIso < 10.910+0.0148*pho.Et+0.000017*pho.Et*pho.Et ):
            return True
    elif ( 1.566 < abs(pho.Eta) < 2.5  and pho.HoverE < 0.0481 and pho.PFChIso < 1.011 and pho.PFNeuIso < 5.931+0.0163*pho.Et+0.000014*pho.Et*pho.Et ):
            return True
    else:
            return False

def isMergedPhoton(ipho, pho_list):
    pho = pho_list[ipho]
    merge = True
    for ip, p in enumerate(pho_list):
        if(ip != ipho):
            if (pho.matchGenID == p.matchGenID):
                merge = False
                break
    return merge


class recoPhoton():
    def __init__(self, ggTree, ireco):
        # From ggTree parameters
        self.E = ggTree.phoE[ireco]
        self.Et = ggTree.phoEt[ireco]
        self.Eta = ggTree.phoEta[ireco]
        self.Phi = ggTree.phoPhi[ireco]
        self.SigmaIEtaIEtaFull5x5 = ggTree.phoSigmaIEtaIEtaFull5x5[ireco]
        self.SigmaIEtaIPhiFull5x5 = ggTree.phoSigmaIEtaIPhiFull5x5[ireco]
        self.SigmaIPhiIPhiFull5x5 = ggTree.phoSigmaIPhiIPhiFull5x5[ireco]
        self.E2x2Full5x5 = ggTree.phoE2x2Full5x5[ireco]
        self.E5x5Full5x5 = ggTree.phoE5x5Full5x5[ireco]
        self.R9Full5x5 = ggTree.phoR9Full5x5[ireco]
        self.PFChWorstIso = ggTree.phoPFChWorstIso[ireco]
        self.IDMVA = ggTree.phoIDMVA[ireco]
        self.SCE = ggTree.phoSCE[ireco]
        self.SCRawE = ggTree.phoSCRawE[ireco]
        self.ESEnP1 = ggTree.phoESEnP1[ireco]
        self.ESEnP2 = ggTree.phoESEnP2[ireco]
        self.SCEta = ggTree.phoSCEta[ireco]
        self.SCPhi = ggTree.phoSCPhi[ireco]
        self.SCEtaWidth = ggTree.phoSCEtaWidth[ireco]
        self.SCPhiWidth = ggTree.phoSCPhiWidth[ireco]
        self.SCBrem = ggTree.phoSCBrem[ireco]
        self.hasPixelSeed = ggTree.phohasPixelSeed[ireco]
        self.EleVeto = ggTree.phoEleVeto[ireco]
        self.R9 = ggTree.phoR9[ireco]
        self.HoverE = ggTree.phoHoverE[ireco]
        self.PFChIso = ggTree.phoPFChIso[ireco]
        self.PFPhoIso = ggTree.phoPFPhoIso[ireco]
        self.PFNeuIso = ggTree.phoPFNeuIso[ireco]
        self.phoIDbit = ggTree.phoIDbit[ireco]
        self.SeedEnergy = ggTree.phoSeedEnergy[ireco]
        # self defined parameters
        self.matchGenID = -1
        self.P4 = ROOT.TLorentzVector()
        self.P4.SetPtEtaPhiE(self.Et, self.Eta, self.Phi, self.E)
        self.gen = None

    def setMatchGen(self, gen):
        self.gen = gen

class mcPhoton():
    def __init__(self, ggTree, imc):
        # from ggNtuple variable
        self.E = ggTree.mcE[imc]
        self.Et = ggTree.mcEt[imc]
        self.Eta = ggTree.mcEta[imc]
        self.Phi = ggTree.mcPhi[imc]
        self.PID = ggTree.mcPID[imc]
        self.MomPID = ggTree.mcMomPID[imc]
        # self defined parameters
        self.P4 = ROOT.TLorentzVector()
        self.P4.SetPtEtaPhiE( self.Et, self.Eta, self.Phi, self.E )
        self.genPair_ptRatio = -1
        self.genPair_deltaR = -1
        self.genPair_mass = -1
        self.genPair_Et = -1
        self.genPair_Eta = -1
        self.genPair_Phi = -1

    def setGenPair(self, genpho):
        self.genPair = genpho  # HACK I originally wanted to directly assign a genpair to a mcPhoton but I don't know how to initialize an empty class object
        self.genPair_ptRatio = min(genpho.Et, self.Et) / max(genpho.Et, self.Et)
        self.genPair_deltaR = deltaR(self, genpho)
        self.genPair_mass = (self.P4+genpho.P4).M()
        self.genPair_Et = genpho.Et
        self.genPair_Eta = genpho.Eta
        self.genPair_Phi = genpho.Phi

if __name__ == "__main__":

        logging.basicConfig(format='%(message)s', stream=sys.stderr, level=logging.DEBUG)

        parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
        parser.add_argument('in_filenames',nargs="+",help='input filenames')
        parser.add_argument('--prefix','-p',default='',help='file prefix')
        parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
        args = parser.parse_args()

        # Read infile
        ggChain = ROOT.TChain("ggNtuplizer/EventTree")
        for in_filename in args.in_filenames:
            logging.info("Adding fle: {}".format(in_filename))
            ggChain.Add(in_filename)

        # info
        year = (args.in_filenames[0].split("MC_Sample/")[1].split("/")[0])
        if ( args.in_filenames[0].find("ZPrimeToZHMuMuGG") != -1):
            a_mass_nominal = int(args.in_filenames[0].split(".root")[0].split("M_")[1]) #~/eos/mergedPhtonID/MC_Sample/2017/ZPrimeToZHMuMuGG/ggtree_mc_ZPrimeToZHMuMuGG_M_1000.root
            outputdir = "plots/ZPrimeToZHMuMuGG/{}/M{}".format(year, a_mass_nominal)
        elif (args.in_filenames[0].find("Haa") != -1):
            a_mass_nominal = float(args.in_filenames[0].split(".root")[0].split("ggtree_mc_")[1].replace("p",".")) # ~/eos/mergedPhtonID/MC_Sample/2017/Haa/ggtree_mc_0p1.root
            a_mass_nominal_text = args.in_filenames[0].split(".root")[0].split("ggtree_mc_")[1]
            outputdir = "plots/Haa/{}/M{}".format(year, a_mass_nominal_text)
        logging.info("Processing year: {}, m(a): {}".format(year, a_mass_nominal))
        logging.info("Output directory: {}".format(outputdir))

        # create outputdir if needed
        try:
            os.makedirs(outputdir)
        except OSError:
            logging.info("Creation of the directory %s failed" % outputdir)
        else:
            logging.info("Successfully created the directory %s" % outputdir)

        # outfile
        outfile = ROOT.TFile("{}/output.root".format(outputdir), 'recreate')

        # Define histos
        h_pt_gen = ROOT.TH1D("h_pt_gen","GenPhoton pt",100, 0, 3000)
        h_eta_gen = ROOT.TH1D("h_eta_gen","GenPhoton eta",60,-3,3)
        h_phi_gen = ROOT.TH1D("h_phi_gen","GenPhoton phi",100,-5,5)
        h_deltaR_gen = ROOT.TH1D("h_deltaR_gen","",100,0,0.8)
        h_mass_gen = ROOT.TH1D("h_mass_gen","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_ptRatio_gen = ROOT.TH1D("h_ptRatio_gen","",101,0,1.01)
        h_pt_deltaR_gen = ROOT.TH2D("h_pt_deltaR_gen","",100,0,200,100,0,a_mass_nominal*0.3)

        h_pt_reco = ROOT.TH1D("h_pt_reco","RecoPhoton pt",100, 0, 200)
        h_eta_reco = ROOT.TH1D("h_eta_reco","RecoPhoton eta",60,-3,3)
        h_phi_reco = ROOT.TH1D("h_phi_reco","RecoPhoton phi",100,-5,5)
        # h_pt_reco_reducedPhotonID = ROOT.TH1D("h_pt_reco_reducedPhotonID","RecoPhoton pt",100, 0, 200)
        # h_pt_reco_LoosePhotonID = ROOT.TH1D("h_pt_reco_LoosePhotonID","RecoPhoton pt",100, 0, 200)
        # h_eta_reco_reducedPhotonID = ROOT.TH1D("h_eta_reco_reducedPhotonID","RecoPhoton eta",60,-3,3)
        # h_eta_reco_LoosePhotonID = ROOT.TH1D("h_eta_reco_LoosePhotonID","RecoPhoton eta",60,-3,3)
        # h_phi_reco_reducedPhotonID = ROOT.TH1D("h_phi_reco_reducedPhotonID","RecoPhoton phi",100,-5,5)
        # h_phi_reco_LoosePhotonID = ROOT.TH1D("h_phi_reco_LoosePhotonID","RecoPhoton phi",100,-5,5)

        # ========== Resolved Photon pair histos ========== #
        h_resolved_mass_reco = ROOT.TH1D("h_resolved_mass_reco","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        h_resolved_deltaR_reco = ROOT.TH1D("h_resolved_deltaR_reco","",100,0,0.2)
        h_resolved_ptRatio_reco = ROOT.TH1D("h_resolved_ptRatio_reco","",101,0,1.01)
        h_resolved_ptRatio_gen = ROOT.TH1D("h_resolved_ptRatio_gen","",101,0,1.01)
        h_resolved_pt_reco = ROOT.TH1D("h_resolved_pt_reco","",100, 0, 200)
        h_resolved_eta_gen = ROOT.TH1D("h_resolved_eta_gen","",60,-3,3)
        h_resolved_phi_gen = ROOT.TH1D("h_resolved_phi_gen","",100,-5,5)
        h_resolved_pt_deltaR_gen = ROOT.TH2D("h_resolved_pt_deltaR_gen","",100,0,200,100,0,a_mass_nominal*0.3)
        h_resolved_ptRatio_deltaR_gen = ROOT.TH2D("h_resolved_ptRatio_deltaR_gen","",101,0,1.01,100,0,a_mass_nominal*0.3)

        # # Loose Photon ID
        # h_resolved_mass_reco_LoosePhotonID = ROOT.TH1D("h_resolved_mass_reco_LoosePhotonID","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        # h_resolved_deltaR_reco_LoosePhotonID = ROOT.TH1D("h_resolved_deltaR_reco_LoosePhotonID","",100,0,a_mass_nominal*0.3)
        # h_resolved_ptRatio_reco_LoosePhotonID = ROOT.TH1D("h_resolved_ptRatio_reco_LoosePhotonID","",101,0,1.01)
        # h_resolved_pt_reco_LoosePhotonID = ROOT.TH1D("h_resolved_pt_reco_LoosePhotonID","",100, 0, 200)

        # # Reduced Photon ID
        # h_resolved_mass_reco_reducedPhotonID = ROOT.TH1D("h_resolved_mass_reco_reducedPhotonID","",100,a_mass_nominal*0.6,a_mass_nominal*1.6)
        # h_resolved_deltaR_reco_reducedPhotonID = ROOT.TH1D("h_resolved_deltaR_reco_reducedPhotonID","",100,0,a_mass_nominal*0.3)
        # h_resolved_ptRatio_reco_reducedPhotonID = ROOT.TH1D("h_resolved_ptRatio_reco_reducedPhotonID","",101,0,1.01)
        # h_resolved_pt_reco_reducedPhotonID = ROOT.TH1D("h_resolved_pt_reco_reducedPhotonID","",100, 0, 200)

        # ========== Merged histos ========== #
        h_merged_pt_reco = ROOT.TH1D("h_merged_pt_reco","",100, 0, 200)
        h_merged_eta_reco = ROOT.TH1D("h_merged_eta_reco","",60,-3,3)
        h_merged_phi_reco = ROOT.TH1D("h_merged_phi_reco","",100,-5,5)
        h_merged_ptRatio_gen = ROOT.TH1D("h_merged_ptRatio_gen","",101,0,1.01)
        h_merged_pt_deltaR_gen = ROOT.TH2D("h_merged_pt_deltaR_gen","",100,0,200,100,0,a_mass_nominal*0.3)
        h_merged_ptRatio_deltaR_gen = ROOT.TH2D("h_merged_ptRatio_deltaR_gen","",101,0,1.01,100,0,a_mass_nominal*0.3)

        h_merged_dpt = ROOT.TH1D("h_merged_dpt","",100,0,10)




        for event_nm in range(ggChain.GetEntries()):
                if (event_nm % args.report == 0):
                    logging.info("Processing event {}/{} ({:.1f}%) ".format(event_nm, ggChain.GetEntries(), 100*event_nm/ggChain.GetEntries()))

                # read parameters
                ggChain.GetEntry(event_nm)


                genpho  = [] # Prompt gen photon
                recopho = [] # matched reco photon list

                for imc in range(ggChain.nMC):
                    if ( ggChain.mcPID[imc] == 22 and (ggChain.mcMomPID[imc] == 25 or ggChain.mcMomPID[imc] == 35) ): # ggChain.mcStatusFlag[imc]>>1&1 == 1 ):
                        photon = mcPhoton(ggChain, imc)
                        logging.debug ("PdgId : {}   pt : {}  eta : {}   phi : {}  mother : {}" .format( photon.PID, photon.Et, photon.Eta, photon.Phi, photon.MomPID ))
                        genpho.append(photon)

                for i in range(len(genpho)/2):
                    pho1 = genpho[i*2]
                    pho2 = genpho[i*2+1]
                    pho1.setGenPair(pho2)
                    pho2.setGenPair(pho1)
                    logging.debug("mass : {}  dR : {}".format(pho1.genPair_mass, pho1.genPair_deltaR))
                    if (pho1.Et>10 and pho2.Et>10 and abs(pho1.Eta) < 1.4442 and abs(pho2.Eta) < 1.4442):
                        h_pt_gen.Fill(pho1.Et)
                        h_eta_gen.Fill(pho1.Eta)
                        h_phi_gen.Fill(pho1.Phi)
                        h_pt_gen.Fill(pho2.Et)
                        h_eta_gen.Fill(pho2.Eta)
                        h_phi_gen.Fill(pho2.Phi)
                        h_deltaR_gen.Fill(pho1.genPair_deltaR)
                        h_mass_gen.Fill(pho1.genPair_mass)
                        h_ptRatio_gen.Fill(pho1.genPair_ptRatio)
                        h_pt_deltaR_gen.Fill(pho1.Et, pho1.genPair_deltaR)
                        h_pt_deltaR_gen.Fill(pho2.Et, pho2.genPair_deltaR)


                logging.debug ("        reco:photon")
                # Gen vs Reco matching
                for ireco in range(ggChain.nPho):
                        # Photon ID selection
                        logging.debug("        pt : {}  eta : {}   phi : {}" .format(ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                        photon = recoPhoton(ggChain,ireco)
                        # Gen vs Reco matching
                        for igen, gen in enumerate(genpho):
                                if ( deltaR(gen, photon) < 0.15 ):
                                        photon.matchGenID = igen
                                        photon.setMatchGen(gen)
                                        recopho.append(photon)
                                        h_pt_reco.Fill(photon.Et)
                                        h_eta_reco.Fill(photon.Eta)
                                        h_phi_reco.Fill(photon.Phi)
                                        # if ( photon.phoIDbit >>0&1 == 1 ):
                                        #         h_pt_reco_LoosePhotonID.Fill(photon.Et)
                                        #         h_eta_reco_LoosePhotonID.Fill(photon.Eta)
                                        #         h_phi_reco_LoosePhotonID.Fill(photon.Phi)
                                        # if ( reducedPhotonID(photon) ):
                                        #         h_pt_reco_reducedPhotonID.Fill(photon.Et)
                                        #         h_eta_reco_reducedPhotonID.Fill(photon.Eta)
                                        #         h_phi_reco_reducedPhotonID.Fill(photon.Phi)
                                        break

                # Reco photon pair
                for i, reco1 in enumerate(recopho):
                        # Separate to Merged and Resolved Catagories
                        if ( isMergedPhoton(i, recopho) ):
                                logging.debug("        merged photon pair!!!! i = {} pt = {}  eta = {}".format(i, reco1.Et, reco1.Eta))
                                h_merged_pt_reco.Fill(reco1.Et)
                                h_merged_eta_reco.Fill(reco1.Eta)
                                h_merged_phi_reco.Fill(reco1.Phi)
                                h_merged_pt_deltaR_gen.Fill(genpho[reco1.matchGenID].genPair_Et, genpho[reco1.matchGenID].genPair_deltaR)
                                h_merged_ptRatio_gen.Fill( reco1.gen.genPair_ptRatio )
                                h_merged_ptRatio_deltaR_gen.Fill(reco1.gen.genPair_ptRatio, reco1.gen.genPair_deltaR)

                                h_merged_dpt.Fill( (reco1.Et - genpho[reco1.matchGenID].Et)/ genpho[reco1.matchGenID].Et)
                        else:
                                for j in range(i+1, len(recopho)):
                                        reco2 = recopho[j]
                                        mass_photonPair = (reco1.P4 + reco2.P4).M()
                                        delta_photonPair = deltaR(reco1, reco2)

                                        # Matched photon pair
                                        if (reco1.matchGenID == 0 and reco2.matchGenID == 1):
                                                logging.debug("        mass1 = {}, deltaR1 = {}".format(mass_photonPair, delta_photonPair))
                                                h_resolved_mass_reco.Fill(mass_photonPair)
                                                h_resolved_deltaR_reco.Fill(delta_photonPair)
                                        elif (reco1.matchGenID == 2 and reco2.matchGenID == 3):
                                                logging.debug ("        mass2 = {}, deltaR2 = {}".format(mass_photonPair, delta_photonPair))
                                                h_resolved_mass_reco.Fill(mass_photonPair)
                                                h_resolved_deltaR_reco.Fill(delta_photonPair)
                                        elif (reco1.matchGenID == reco2.matchGenID):
                                                logging.debug("        very close reco photon pair!!!! i = {}, j = {}, mass = {}, deltaR = {}".format(i, j, mass_photonPair, delta_photonPair))
                                                h_resolved_mass_reco.Fill(mass_photonPair)
                                                h_resolved_deltaR_reco.Fill(delta_photonPair)
                                                h_resolved_ptRatio_reco.Fill(min(reco1.Et, reco2.Et)/max(reco1.Et, reco2.Et))
                                                h_resolved_pt_reco.Fill(reco1.Et)
                                                h_resolved_pt_reco.Fill(reco2.Et)
                                                h_resolved_ptRatio_gen.Fill( reco2.gen.genPair_ptRatio )
                                                h_resolved_pt_deltaR_gen.Fill(reco1.gen.Et, reco1.gen.genPair_deltaR)
                                                h_resolved_pt_deltaR_gen.Fill(reco2.gen.Et, reco2.gen.genPair_deltaR)
                                                h_resolved_ptRatio_deltaR_gen.Fill(reco1.gen.genPair_ptRatio, reco1.gen.genPair_deltaR)
                                                h_resolved_ptRatio_deltaR_gen.Fill(reco2.gen.genPair_ptRatio, reco2.gen.genPair_deltaR)
                                                h_resolved_eta_gen.Fill(reco1.gen.Eta)
                                                h_resolved_eta_gen.Fill(reco2.gen.Eta)
                                                h_resolved_phi_gen.Fill(reco1.gen.Phi)
                                                h_resolved_phi_gen.Fill(reco2.gen.Phi)

                                        # if ( reco1.phoIDbit>>0&1 == 1 and reco2.phoIDbit>>0&1 == 1):
                                        #         if (reco1.matchGenID == reco2.matchGenID):
                                        #                 h_resolved_mass_reco_LoosePhotonID.Fill(mass_photonPair)
                                        #                 h_resolved_deltaR_reco_LoosePhotonID.Fill(delta_photonPair)
                                        #                 h_resolved_ptRatio_reco_LoosePhotonID.Fill(min(reco1.Et, reco2.Et)/max(reco1.Et, reco2.Et))
                                        #                 h_resolved_pt_reco_LoosePhotonID.Fill(reco1.Et)
                                        #                 h_resolved_pt_reco_LoosePhotonID.Fill(reco2.Et)

                                        # if ( reducedPhotonID(reco1) and reducedPhotonID(reco2) ):
                                        #         if (reco1.matchGenID == reco2.matchGenID):
                                        #                 h_resolved_mass_reco_reducedPhotonID.Fill(mass_photonPair)
                                        #                 h_resolved_deltaR_reco_reducedPhotonID.Fill(delta_photonPair)
                                        #                 h_resolved_ptRatio_reco_reducedPhotonID.Fill(min(reco1.Et, reco2.Et)/max(reco1.Et, reco2.Et))
                                        #                 h_resolved_pt_reco_reducedPhotonID.Fill(reco1.Et)
                                        #                 h_resolved_pt_reco_reducedPhotonID.Fill(reco2.Et)



        print ("Matching efficiency                                  : {}".format( h_pt_reco.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Matching + reduced Photon ID efficiency              : {}".format( h_pt_reco_reducedPhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Matching + Loose Photon ID efficiency                : {}".format( h_pt_reco_LoosePhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        print ("Close Matched Photons efficiency                     : {}".format( h_resolved_pt_reco.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Close Matched Photons + reduced Photon ID efficiency : {}".format( h_resolved_pt_reco_reducedPhotonID.GetEntries() / h_pt_gen.GetEntries() ))
        # print ("Close Matched Photons + Loose Photon ID efficiency   : {}".format( h_resolved_pt_reco_LoosePhotonID.GetEntries() / h_pt_gen.GetEntries() ))

        h_pt_gen.Write()
        h_deltaR_gen.Write()

        # ========== Plots ========== #
        c1 = ROOT.TCanvas()

        # c1.SaveAs("{}/deltaR_gen.png".format(outputdir,a_mass_nominal))

        h_resolved_deltaR_reco.SetTitle('photon pair deltaR')
        h_resolved_deltaR_reco.GetXaxis().SetTitle('dR')
        h_resolved_deltaR_reco.GetYaxis().SetTitle('Events')
        h_resolved_deltaR_reco.Draw()
        h_resolved_deltaR_reco.SetLineColor(2)
        h_deltaR_gen.Draw("same")
        legend_deltaR = ROOT.TLegend(0.55,0.55,0.8,0.8)
        legend_deltaR.AddEntry(h_deltaR_gen,"gen","l")
        legend_deltaR.AddEntry(h_resolved_deltaR_reco,"reco","l")
        legend_deltaR.Draw("same")
        c1.Update()
        c1.SaveAs("{}/deltaR_reco.png".format(outputdir,a_mass_nominal))

        h_eta_gen.SetTitle('Eta')
        h_eta_gen.GetXaxis().SetTitle('Eta')
        h_eta_gen.GetYaxis().SetTitle('Events')
        h_eta_gen.Draw()
        h_eta_reco.SetLineColor(2)
        h_eta_reco.Draw('same')
        h_merged_eta_reco.SetLineColor(4)
        h_merged_eta_reco.Draw('same')
        h_resolved_eta_gen.SetLineColor(6)
        h_resolved_eta_gen.Draw('same')
        legend_eta = ROOT.TLegend(0.65,0.65,0.85,0.85)
        legend_eta.AddEntry(h_eta_gen,"gen","l")
        legend_eta.AddEntry(h_eta_reco,"reco","l")
        legend_eta.AddEntry(h_merged_eta_reco,"gen merged","l")
        legend_eta.AddEntry(h_resolved_eta_gen,"gen resolved","l")
        legend_eta.Draw('same')
        c1.Update()
        c1.SaveAs("{}/eta.png".format(outputdir,a_mass_nominal))

        h_phi_gen.SetTitle('Phi')
        h_phi_gen.GetXaxis().SetTitle('Phi')
        h_phi_gen.GetYaxis().SetTitle('Events')
        h_phi_gen.Draw()
        h_phi_reco.SetLineColor(2)
        h_phi_reco.Draw('same')
        h_merged_phi_reco.SetLineColor(4)
        h_merged_phi_reco.Draw('same')
        h_resolved_phi_gen.SetLineColor(6)
        h_resolved_phi_gen.Draw('same')
        legend_phi = ROOT.TLegend(0.65,0.65,0.85,0.85)
        legend_phi.AddEntry(h_phi_gen,"gen","l")
        legend_phi.AddEntry(h_phi_reco,"reco","l")
        legend_phi.AddEntry(h_merged_phi_reco,"gen merged","l")
        legend_phi.AddEntry(h_resolved_phi_gen,"gen resolved","l")
        legend_phi.Draw('same')
        c1.Update()
        c1.SaveAs("{}/phi.png".format(outputdir,a_mass_nominal))
        # h_mass_gen.SetTitle('Gen m(a)={}GeV'.format(outputdir,a_mass_nominal))
        # h_mass_gen.GetXaxis().SetTitle('mass')
        # h_mass_gen.Draw()
        # c1.Update()
        # c1.SaveAs("{}/mass_gen.png".format(outputdir,a_mass_nominal))

        # h_mass_reco.SetTitle('Reco m(a)={}GeV'.format(outputdir,a_mass_nominal))
        # h_mass_reco.GetXaxis().SetTitle('mass')
        # h_mass_reco.Draw()
        # c1.Update()
        # c1.SaveAs("{}/mass_reco.png".format(outputdir,a_mass_nominal))

        h_pt_reco.SetTitle('RecoPhoton pt')
        h_pt_reco.GetXaxis().SetTitle('pt')
        h_pt_reco.Draw()
        h_resolved_pt_reco.SetLineColor(6)
        h_resolved_pt_reco.Draw("same")
        h_merged_pt_reco.SetLineColor(28)
        h_merged_pt_reco.Draw("same")
        legend_pt = ROOT.TLegend(0.45,0.5,0.75,0.8)
        legend_pt.AddEntry(h_pt_reco,"Total reco","l")
        legend_pt.AddEntry(h_resolved_pt_reco,"Resolved","l")
        legend_pt.AddEntry(h_merged_pt_reco,"Merged","l")
        legend_pt.Draw("same")
        c1.Update()
        c1.SaveAs("{}/pt_reco.png".format(outputdir,a_mass_nominal))

        h_ptRatio_gen.SetTitle('photon pair, pt ratio')
        h_ptRatio_gen.GetXaxis().SetTitle('pt_min / pt_max')
        h_ptRatio_gen.Draw()
        h_merged_ptRatio_gen.SetLineColor(2)
        h_merged_ptRatio_gen.Draw("same")
        h_resolved_ptRatio_gen.SetLineColor(4)
        h_resolved_ptRatio_gen.Draw("same")
        # h_resolved_ptRatio_reco.SetLineColor(6)
        # h_resolved_ptRatio_reco.Draw("same")
        legend_ptRatio = ROOT.TLegend(0.45,0.5,0.75,0.8)
        legend_ptRatio.AddEntry(h_ptRatio_gen,"Total gen","l")
        legend_ptRatio.AddEntry(h_merged_ptRatio_gen,"merged gen","l")
        legend_ptRatio.AddEntry(h_resolved_ptRatio_gen,"resolved gen","l")
        # legend_ptRatio.AddEntry(h_resolved_ptRatio_reco,"resolved reco","l")
        legend_ptRatio.Draw("same")
        c1.Update()
        c1.SaveAs("{}/pt_ratio.png".format(outputdir,a_mass_nominal))

        h_efficiency_resolved_pt = h_resolved_pt_reco.Clone("h_efficiency_resolved_pt")
        h_efficiency_resolved_pt.Divide(h_pt_gen)
        h_efficiency_resolved_pt.SetTitle('Resolved efficiency vs pt')
        h_efficiency_resolved_pt.Draw("P")
        c1.Update()
        c1.SaveAs("{}/resoloved_efficiency_pt.png".format(outputdir,a_mass_nominal))

        h_efficiency_merged_pt = h_merged_pt_reco.Clone("h_efficiency_merged_pt")
        h_efficiency_merged_pt.Divide(h_pt_gen)
        h_efficiency_merged_pt.SetTitle('Merged efficiency vs pt')
        h_efficiency_merged_pt.Draw("P")
        c1.Update()
        c1.SaveAs("{}/merged_efficiency_pt.png".format(outputdir,a_mass_nominal))

        h_merged_pt_deltaR_gen.SetTitle('Merged deltaR vs Pt')
        h_merged_pt_deltaR_gen.GetXaxis().SetTitle('Pt')
        h_merged_pt_deltaR_gen.GetYaxis().SetTitle('dR')
        h_merged_pt_deltaR_gen.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/merged_dR_pt.png".format(outputdir,a_mass_nominal))

        h_resolved_pt_deltaR_gen.SetTitle('Resolved deltaR vs Pt')
        h_resolved_pt_deltaR_gen.GetXaxis().SetTitle('Pt')
        h_resolved_pt_deltaR_gen.GetYaxis().SetTitle('dR')
        h_resolved_pt_deltaR_gen.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/resolved_dR_pt.png".format(outputdir,a_mass_nominal))

        h_merged_ptRatio_deltaR_gen.SetTitle('Merged deltaR vs PtRatio')
        h_merged_ptRatio_deltaR_gen.GetXaxis().SetTitle('PtRatio')
        h_merged_ptRatio_deltaR_gen.GetYaxis().SetTitle('dR')
        h_merged_ptRatio_deltaR_gen.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/merged_dR_ptRatio.png".format(outputdir,a_mass_nominal))

        h_resolved_ptRatio_deltaR_gen.SetTitle('Resolved deltaR vs PtRatio')
        h_resolved_ptRatio_deltaR_gen.GetXaxis().SetTitle('PtRatio')
        h_resolved_ptRatio_deltaR_gen.GetYaxis().SetTitle('dR')
        h_resolved_ptRatio_deltaR_gen.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/resolved_dR_ptRatio.png".format(outputdir,a_mass_nominal))

        h_pt_deltaR_gen.SetTitle('deltaR vs Pt')
        h_pt_deltaR_gen.GetXaxis().SetTitle('Pt')
        h_pt_deltaR_gen.GetYaxis().SetTitle('dR')
        h_pt_deltaR_gen.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/dR_pt.png".format(outputdir,a_mass_nominal))

        h_efficiency_resolved_pt_deltaR = h_resolved_pt_deltaR_gen.Clone("h_efficiency_resolved_pt_deltaR")
        h_efficiency_resolved_pt_deltaR.Divide(h_pt_deltaR_gen)
        h_efficiency_resolved_pt_deltaR.SetTitle('Resolved efficiency vs pt&deltaR')
        h_efficiency_resolved_pt_deltaR.Draw("colz")
        c1.Update()
        c1.SaveAs("{}/resolved_efficiency_pt_dR.png".format(outputdir,a_mass_nominal))

        h_merged_dpt.GetXaxis().SetTitle('dpt')
        h_merged_dpt.Draw()
        c1.Update()
        c1.SaveAs("{}/h_merged_dpt.png".format(outputdir,a_mass_nominal))

        outfile.Write("",ROOT.TFile.kOverwrite)
        outfile.Close()
