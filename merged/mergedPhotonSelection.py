#!/usr/bin/env python

from __future__ import division

import ROOT
import os, sys, logging
ROOT.gROOT.SetBatch(True)
import argparse
import math
from array import array

def deltaR(pho1, pho2):
    return math.sqrt( (pho1.Eta-pho2.Eta)**2 + (pho1.Phi-pho2.Phi)**2 )

def deltaPt(gen, reco):
    return abs((reco.Et-gen.Et)/gen.Et)

def recoCommonSelection(reco):
    passed = False
    if ( abs(reco.Eta)<1.4442 and reco.Et > 10 and reco.PFChIso < 25):
        passed = True
    return passed

# Define reco & gen matching selection
def matchGenSelection(gen, reco):
    passed = False
    if ( deltaR(gen, reco) < 0.15 ):
        passed = True
    return passed

# Define selection to reco that does 'not' match to gen
def no_matchGenSelection(gen, reco):
    passed = False
    if ( deltaR(gen, reco) > 0.4 ):
        passed = True
    return passed


def passPromptMergedSelection(reco1, reco2):
    passed = True
    if ((reco1.matchGenID == reco2.matchGenID) or (reco1.matchGenID == 0 and reco2.matchGenID == 1) or (reco1.matchGenID == 1 and reco2.matchGenID == 0) or (reco1.matchGenID == 2 and reco2.matchGenID == 3) or (reco1.matchGenID == 3 and reco2.matchGenID == 2)):
        passed = False
    return passed

def passFakeMergedSelection(reco1, reco2):
    passed = False
    if (reco1.matchGenFlag == True):
        passed = True
    elif (reco1.matchGenFlag == False):
        if (reco2.matchGenFlag == True):
            passed = True
        elif (reco2.matchGenFlag == False and reco1.Et > reco2.Et):
            passed = True
    return passed


class recoPhoton():
    def __init__(self, ggChain, ireco):
        # From ggChain parameters
        self.E = ggChain.phoE[ireco]
        self.Et = ggChain.phoEt[ireco]
        self.Eta = ggChain.phoEta[ireco]
        self.Phi = ggChain.phoPhi[ireco]
        self.SigmaIEtaIEtaFull5x5 = ggChain.phoSigmaIEtaIEtaFull5x5[ireco]
        self.HoverE = ggChain.phoHoverE[ireco]
        self.PFChIso = ggChain.phoPFChIso[ireco]
        self.PFPhoIso = ggChain.phoPFPhoIso[ireco]
        self.PFNeuIso = ggChain.phoPFNeuIso[ireco]
        self.phoIDbit = ggChain.phoIDbit[ireco]
        # self defined parameters
        self.matchGenID = -1
        self.recoID = -1
        self.P4 = ROOT.TLorentzVector()
        self.P4.SetPtEtaPhiE(self.Et, self.Eta, self.Phi, self.E)
        self.gen = None
        self.matchGenFlag = False

    def setMatchGen(self, gen):
        self.gen = gen


class mcPhoton():
    def __init__(self, ggChain, imc):
        # from ggNtuple variable
        self.E = ggChain.mcE[imc]
        self.Et = ggChain.mcEt[imc]
        self.Eta = ggChain.mcEta[imc]
        self.Phi = ggChain.mcPhi[imc]
        self.PID = ggChain.mcPID[imc]
        self.MomPID = ggChain.mcMomPID[imc]


if __name__ == "__main__":

    logging.basicConfig(format='%(message)s', stream=sys.stderr, level=logging.INFO)

    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('in_filenames',nargs="+",help='input filenames')
    parser.add_argument('--prefix','-p',default='',help='file prefix')
    parser.add_argument('--out','-o',default="output.root",help='output filename')
    parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
    parser.add_argument('--savePrompt','-s',action="store_true",help='store prompt or fake photons')
    args = parser.parse_args()

    # Read infile
    ggChain = ROOT.TChain("ggNtuplizer/EventTree")
    for in_filename in args.in_filenames:
        print("Adding fle: {}".format(in_filename))
        ggChain.Add(in_filename)

    # Initialize output tree
    if (args.savePrompt==True):
        print("Selecting prompt photons to promptMerged{}.root".format(args.prefix))
        outfile = ROOT.TFile("promptMerged{}.root".format(args.prefix), 'recreate')
    elif (args.savePrompt==False):
        print("Selecting fake photons to fakeMerged{}.root".format(args.prefix))
        outfile = ROOT.TFile("fakeMerged{}.root".format(args.prefix), 'recreate')
    outtree = ROOT.TTree("recoPhoton", "Skimmed reco photon")

    E = array('d',[0])
    SigmaE = array('d',[0])
    Et = array('d',[0])
    Eta = array('d',[0])
    Phi = array('d',[0])
    CalibE = array('d',[0])
    CalibEt = array('d',[0])
    SCE = array('d',[0])
    SCRawE = array('d',[0])
    ESEnP1 = array('d',[0])
    ESEnP2 = array('d',[0])
    SCEta = array('d',[0])
    SCPhi = array('d',[0])
    SCEtaWidth = array('d',[0])
    SCPhiWidth = array('d',[0])
    SCBrem = array('d',[0])
    hasPixelSeed = array('d',[0])
    EleVeto = array('d',[0])
    R9 = array('d',[0])
    HoverE = array('d',[0])
    ESEffSigmaRR = array('d',[0])
    SigmaIEtaIEtaFull5x5 = array('d',[0])
    SigmaIEtaIPhiFull5x5 = array('d',[0])
    SigmaIPhiIPhiFull5x5 = array('d',[0])
    E2x2Full5x5 = array('d',[0])
    E5x5Full5x5 = array('d',[0])
    R9Full5x5 = array('d',[0])
    PFChIso = array('d',[0])
    PFChPVIso = array('d',[0])
    PFPhoIso = array('d',[0])
    PFNeuIso = array('d',[0])
    PFChWorstIso = array('d',[0])
    PFChWorstVetoIso = array('d',[0])
    EcalPFClusterIso = array('d',[0])
    HcalPFClusterIso = array('d',[0])
    IDMVA = array('d',[0])
    FiredSingleTrgs = array('d',[0])
    FiredDoubleTrgs = array('d',[0])
    FiredTripleTrgs = array('d',[0])
    FiredL1Trgs = array('d',[0])
    SeedTime = array('d',[0])
    SeedEnergy = array('d',[0])
    MIPTotEnergy = array('d',[0])
    MIPChi2 = array('d',[0])
    MIPSlope = array('d',[0])
    MIPIntercept = array('d',[0])
    MIPNhitCone = array('d',[0])
    MIPIsHalo = array('d',[0])
    IDbit = array('d',[0])
    Scale_stat_up = array('d',[0])
    Scale_stat_dn = array('d',[0])
    Scale_syst_up = array('d',[0])
    Scale_syst_dn = array('d',[0])
    Scale_gain_up = array('d',[0])
    Scale_gain_dn = array('d',[0])
    Resol_rho_up = array('d',[0])
    Resol_rho_dn = array('d',[0])
    Resol_phi_up = array('d',[0])
    Resol_phi_dn = array('d',[0])
    matchGenFlag = array('b',[0])
    outtree.Branch("E", E, 'E/D')
    outtree.Branch("SigmaE",     SigmaE, '    SigmaE/D')
    outtree.Branch("Et", Et, 'Et/D')
    outtree.Branch("Eta", Eta, 'Eta/D')
    outtree.Branch("Phi", Phi, 'Phi/D')
    outtree.Branch("CalibE", CalibE, 'CalibE/D')
    outtree.Branch("CalibEt", CalibEt, 'CalibEt/D')
    outtree.Branch("SCE", SCE, 'SCE/D')
    outtree.Branch("SCRawE", SCRawE, 'SCRawE/D')
    outtree.Branch("ESEnP1", ESEnP1, 'ESEnP1/D')
    outtree.Branch("ESEnP2", ESEnP2, 'ESEnP2/D')
    outtree.Branch("SCEta", SCEta, 'SCEta/D')
    outtree.Branch("SCPhi", SCPhi, 'SCPhi/D')
    outtree.Branch("SCEtaWidth", SCEtaWidth, 'SCEtaWidth/D')
    outtree.Branch("SCPhiWidth", SCPhiWidth, 'SCPhiWidth/D')
    outtree.Branch("SCBrem", SCBrem, 'SCBrem/D')
    outtree.Branch("hasPixelSeed", hasPixelSeed, 'hasPixelSeed/D')
    outtree.Branch("EleVeto", EleVeto, 'EleVeto/D')
    outtree.Branch("R9", R9, 'R9/D')
    outtree.Branch("HoverE", HoverE, 'HoverE/D')
    outtree.Branch("ESEffSigmaRR", ESEffSigmaRR, 'ESEffSigmaRR/D')
    outtree.Branch("SigmaIEtaIEtaFull5x5", SigmaIEtaIEtaFull5x5, 'SigmaIEtaIEtaFull5x5/D')
    outtree.Branch("SigmaIEtaIPhiFull5x5", SigmaIEtaIPhiFull5x5, 'SigmaIEtaIPhiFull5x5/D')
    outtree.Branch("SigmaIPhiIPhiFull5x5", SigmaIPhiIPhiFull5x5, 'SigmaIPhiIPhiFull5x5/D')
    outtree.Branch("E2x2Full5x5", E2x2Full5x5, 'E2x2Full5x5/D')
    outtree.Branch("E5x5Full5x5", E5x5Full5x5, 'E5x5Full5x5/D')
    outtree.Branch("R9Full5x5", R9Full5x5, 'R9Full5x5/D')
    outtree.Branch("PFChIso", PFChIso, 'PFChIso/D')
    outtree.Branch("PFChPVIso", PFChPVIso, 'PFChPVIso/D')
    outtree.Branch("PFPhoIso", PFPhoIso, 'PFPhoIso/D')
    outtree.Branch("PFNeuIso", PFNeuIso, 'PFNeuIso/D')
    outtree.Branch("PFChWorstIso", PFChWorstIso, 'PFChWorstIso/D')
    outtree.Branch("PFChWorstVetoIso", PFChWorstVetoIso, 'PFChWorstVetoIso/D')
    outtree.Branch("EcalPFClusterIso", EcalPFClusterIso, 'EcalPFClusterIso/D')
    outtree.Branch("HcalPFClusterIso", HcalPFClusterIso, 'HcalPFClusterIso/D')
    outtree.Branch("IDMVA", IDMVA, 'IDMVA/D')
    outtree.Branch("FiredSingleTrgs", FiredSingleTrgs, 'FiredSingleTrgs/D')
    outtree.Branch("FiredDoubleTrgs", FiredDoubleTrgs, 'FiredDoubleTrgs/D')
    outtree.Branch("FiredTripleTrgs", FiredTripleTrgs, 'FiredTripleTrgs/D')
    outtree.Branch("FiredL1Trgs", FiredL1Trgs, 'FiredL1Trgs/D')
    outtree.Branch("SeedTime", SeedTime, 'SeedTime/D')
    outtree.Branch("SeedEnergy", SeedEnergy, 'SeedEnergy/D')
    outtree.Branch("MIPTotEnergy", MIPTotEnergy, 'MIPTotEnergy/D')
    outtree.Branch("MIPChi2", MIPChi2, 'MIPChi2/D')
    outtree.Branch("MIPSlope", MIPSlope, 'MIPSlope/D')
    outtree.Branch("MIPIntercept", MIPIntercept, 'MIPIntercept/D')
    outtree.Branch("MIPNhitCone", MIPNhitCone, 'MIPNhitCone/D')
    outtree.Branch("MIPIsHalo", MIPIsHalo, 'MIPIsHalo/D')
    outtree.Branch("IDbit", IDbit, 'IDbit/D')
    outtree.Branch("Scale_stat_up", Scale_stat_up, 'Scale_stat_up/D')
    outtree.Branch("Scale_stat_dn", Scale_stat_dn, 'Scale_stat_dn/D')
    outtree.Branch("Scale_syst_up", Scale_syst_up, 'Scale_syst_up/D')
    outtree.Branch("Scale_syst_dn", Scale_syst_dn, 'Scale_syst_dn/D')
    outtree.Branch("Scale_gain_up", Scale_gain_up, 'Scale_gain_up/D')
    outtree.Branch("Scale_gain_dn", Scale_gain_dn, 'Scale_gain_dn/D')
    outtree.Branch("Resol_rho_up", Resol_rho_up, 'Resol_rho_up/D')
    outtree.Branch("Resol_rho_dn", Resol_rho_dn, 'Resol_rho_dn/D')
    outtree.Branch("Resol_phi_up", Resol_phi_up, 'Resol_phi_up/D')
    outtree.Branch("Resol_phi_dn", Resol_phi_dn, 'Resol_phi_dn/D')
    outtree.Branch("matchGenFlag", matchGenFlag, 'matchGenFlag/O')

    totalRecoPho = 0
    passedRecoPho = 0

    # ==================== Start Loop ==================== #
    for event_nm in range(ggChain.GetEntries()):
        if (event_nm % args.report == 0):
            logging.info("Processing event {}/{} ({:.1f}%) ".format(event_nm, ggChain.GetEntries(), 100*event_nm/ggChain.GetEntries()))

        # read parameters
        ggChain.GetEntry(event_nm)
        totalRecoPho+=ggChain.nPho

        genpho_list = []
        recopho_list = []

        logging.debug("==========Event {}==========".format(event_nm))
        logging.debug("gen:photon")
        for imc in range(ggChain.nMC):
            saveMC = False
            if ( args.savePrompt == True and ggChain.mcPID[imc] == 22 and (ggChain.mcMomPID[imc] == 25 or ggChain.mcMomPID[imc] == 35)):
                saveMC = True
            elif ( args.savePrompt == False and ggChain.mcPID[imc] == 22 and ggChain.mcStatusFlag[imc]>>1&1 == 1 ):
                saveMC = True
            if (  saveMC == True ):
                photon = mcPhoton(ggChain, imc)
                logging.debug ("PdgId : {}   pt : {}  eta : {}   phi : {}  mother : {}" .format( photon.PID, photon.Et, photon.Eta, photon.Phi, photon.MomPID ))
                genpho_list.append(photon)

        logging.debug ("    reco:photon")
        # Gen vs Reco matching
        for ireco in range(ggChain.nPho):
            logging.debug("    {}  pt : {}  eta : {}   phi : {}" .format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
            photon = recoPhoton(ggChain,ireco)
            fakeCandidate = True
            if(recoCommonSelection(photon) == False):  # NOTE Common Selection for both prompt and fake
                break
            for igen, gen in enumerate(genpho_list):
                if (matchGenSelection(gen, photon)):
                    if ( args.savePrompt == True ):
                        logging.debug("    prompt Candidate : {}  pt : {}  eta : {}   phi : {}" .format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                    elif ( args.savePrompt == False ):
                        logging.debug("    fake match gen Candidate : {}  pt : {}  eta : {}   phi : {}" .format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                    photon.matchGenID = igen
                    photon.recoID = ireco
                    photon.matchGenFlag = True
                    photon.setMatchGen(gen)
                    recopho_list.append(photon)
                    fakeCandidate = False
                    break
                if (args.savePrompt == False and no_matchGenSelection(gen, photon) == False):
                    fakeCandidate = False
            if (args.savePrompt == False and fakeCandidate == True ):
                logging.debug("    fake don't match gen Candidate : {}  pt : {}  eta : {}  phi : {}".format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                photon.recoID = ireco
                recopho_list.append(photon)

        for i, reco1 in enumerate(recopho_list):
            saveReco_flag = True
            for j, reco2 in enumerate(recopho_list):
                if ( i != j ):
                    if ( (args.savePrompt==True and passPromptMergedSelection(reco1,reco2)==False ) or (args.savePrompt==False and passFakeMergedSelection(reco1, reco2)==False ) ):
                        saveReco_flag = False
            if (saveReco_flag == True):
                passedRecoPho+=1
                logging.debug("        find merged photon, Et = {} Eta = {}".format(reco1.Et, reco1.Eta))
                # Fill tree
                E[0] = ggChain.phoE[reco1.recoID]
                SigmaE[0] = ggChain.phoSigmaE[reco1.recoID]
                Et[0] = ggChain.phoEt[reco1.recoID]
                Eta[0] = ggChain.phoEta[reco1.recoID]
                Phi[0] = ggChain.phoPhi[reco1.recoID]
                CalibE[0] = ggChain.phoCalibE[reco1.recoID]
                CalibEt[0] = ggChain.phoCalibEt[reco1.recoID]
                SCE[0] = ggChain.phoSCE[reco1.recoID]
                SCRawE[0] = ggChain.phoSCRawE[reco1.recoID]
                ESEnP1[0] = ggChain.phoESEnP1[reco1.recoID]
                ESEnP2[0] = ggChain.phoESEnP2[reco1.recoID]
                SCEta[0] = ggChain.phoSCEta[reco1.recoID]
                SCPhi[0] = ggChain.phoSCPhi[reco1.recoID]
                SCEtaWidth[0] = ggChain.phoSCEtaWidth[reco1.recoID]
                SCPhiWidth[0] = ggChain.phoSCPhiWidth[reco1.recoID]
                SCBrem[0] = ggChain.phoSCBrem[reco1.recoID]
                hasPixelSeed[0] = ggChain.phohasPixelSeed[reco1.recoID]
                EleVeto[0] = ggChain.phoEleVeto[reco1.recoID]
                R9[0] = ggChain.phoR9[reco1.recoID]
                HoverE[0] = ggChain.phoHoverE[reco1.recoID]
                ESEffSigmaRR[0] = ggChain.phoESEffSigmaRR[reco1.recoID]
                SigmaIEtaIEtaFull5x5[0] = ggChain.phoSigmaIEtaIEtaFull5x5[reco1.recoID]
                SigmaIEtaIPhiFull5x5[0] = ggChain.phoSigmaIEtaIPhiFull5x5[reco1.recoID]
                SigmaIPhiIPhiFull5x5[0] = ggChain.phoSigmaIPhiIPhiFull5x5[reco1.recoID]
                E2x2Full5x5[0] = ggChain.phoE2x2Full5x5[reco1.recoID]
                E5x5Full5x5[0] = ggChain.phoE5x5Full5x5[reco1.recoID]
                R9Full5x5[0] = ggChain.phoR9Full5x5[reco1.recoID]
                PFChIso[0] = ggChain.phoPFChIso[reco1.recoID]
                PFPhoIso[0] = ggChain.phoPFPhoIso[reco1.recoID]
                PFNeuIso[0] = ggChain.phoPFNeuIso[reco1.recoID]
                PFChWorstIso[0] = ggChain.phoPFChWorstIso[reco1.recoID]
                IDMVA[0] = ggChain.phoIDMVA[reco1.recoID]
                FiredSingleTrgs[0] = ggChain.phoFiredSingleTrgs[reco1.recoID]
                FiredDoubleTrgs[0] = ggChain.phoFiredDoubleTrgs[reco1.recoID]
                FiredTripleTrgs[0] = ggChain.phoFiredTripleTrgs[reco1.recoID]
                FiredL1Trgs[0] = ggChain.phoFiredL1Trgs[reco1.recoID]
                SeedTime[0] = ggChain.phoSeedTime[reco1.recoID]
                SeedEnergy[0] = ggChain.phoSeedEnergy[reco1.recoID]
                IDbit[0] = ggChain.phoIDbit[reco1.recoID]
                # Scale_stat_up[0] = ggChain.phoScale_stat_up[reco1.recoID]
                # Scale_stat_dn[0] = ggChain.phoScale_stat_dn[reco1.recoID]
                # Scale_syst_up[0] = ggChain.phoScale_syst_up[reco1.recoID]
                # Scale_syst_dn[0] = ggChain.phoScale_syst_dn[reco1.recoID]
                # Scale_gain_up[0] = ggChain.phoScale_gain_up[reco1.recoID]
                # Scale_gain_dn[0] = ggChain.phoScale_gain_dn[reco1.recoID]
                # Resol_rho_up[0] = ggChain.phoResol_rho_up[reco1.recoID]
                # Resol_rho_dn[0] = ggChain.phoResol_rho_dn[reco1.recoID]
                # Resol_phi_up[0] = ggChain.phoResol_phi_up[reco1.recoID]
                # Resol_phi_dn[0] = ggChain.phoResol_phi_dn[reco1.recoID]
                # self defined branch
                matchGenFlag[0] = reco1.matchGenFlag
                outtree.Fill()

    logging.info("Selection efficiency = {:.4f}".format(passedRecoPho / totalRecoPho))

    outfile.Write("",ROOT.TFile.kOverwrite)
    outfile.Close()
