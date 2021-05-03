#!/usr/bin/env python

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

def passPromptPairSelection(reco1, reco2):
    passed = False
    if ((reco1.matchGenID == reco2.matchGenID) or (reco1.matchGenID == 0 and reco2.matchGenID == 1) or (reco1.matchGenID == 1 and reco2.matchGenID == 0) or (reco1.matchGenID == 2 and reco2.matchGenID == 3) or (reco1.matchGenID == 3 and reco2.matchGenID == 2)):
        passed = True
    return passed

def passFakePairSelection(reco1, reco2):
    passed = False
    if ( deltaR(reco1, reco2) < 0.15 ):
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
    parser.add_argument('--report','-r',default=5000,type=int,help='report every x events')
    parser.add_argument('--savePrompt','-s',action="store_true",help='store prompt or fake photons')
    parser.add_argument('--isMerged','-m',action="store_true",help='store merge or resolved photons')
    args = parser.parse_args()

    # Read infile
    ggChain = ROOT.TChain("ggNtuplizer/EventTree")
    for in_filename in args.in_filenames:
        print("Adding {} file".format(in_filename))
        ggChain.Add(in_filename)

    # Initialize output tree
    if (args.savePrompt==True):
        print("Saving prompt photons")
        outfile = ROOT.TFile("prompt.root", 'recreate')
    elif (args.savePrompt==False):
        print("Saving fake photons")
        outfile = ROOT.TFile("fake.root", 'recreate')
    outtree = ROOT.TTree("recoPhoton", "Skimmed reco photon")

    E = array('d',2*[0])
    SigmaE = array('d',2*[0])
    Et = array('d',2*[0])
    Eta = array('d',2*[0])
    Phi = array('d',2*[0])
    CalibE = array('d',2*[0])
    CalibEt = array('d',2*[0])
    SCE = array('d',2*[0])
    SCRawE = array('d',2*[0])
    ESEnP1 = array('d',2*[0])
    ESEnP2 = array('d',2*[0])
    SCEta = array('d',2*[0])
    SCPhi = array('d',2*[0])
    SCEtaWidth = array('d',2*[0])
    SCPhiWidth = array('d',2*[0])
    SCBrem = array('d',2*[0])
    hasPixelSeed = array('d',2*[0])
    EleVeto = array('d',2*[0])
    R9 = array('d',2*[0])
    HoverE = array('d',2*[0])
    ESEffSigmaRR = array('d',2*[0])
    SigmaIEtaIEtaFull5x5 = array('d',2*[0])
    SigmaIEtaIPhiFull5x5 = array('d',2*[0])
    SigmaIPhiIPhiFull5x5 = array('d',2*[0])
    E2x2Full5x5 = array('d',2*[0])
    E5x5Full5x5 = array('d',2*[0])
    R9Full5x5 = array('d',2*[0])
    PFChIso = array('d',2*[0])
    PFChPVIso = array('d',2*[0])
    PFPhoIso = array('d',2*[0])
    PFNeuIso = array('d',2*[0])
    PFChWorstIso = array('d',2*[0])
    PFChWorstVetoIso = array('d',2*[0])
    EcalPFClusterIso = array('d',2*[0])
    HcalPFClusterIso = array('d',2*[0])
    IDMVA = array('d',2*[0])
    FiredSingleTrgs = array('d',2*[0])
    FiredDoubleTrgs = array('d',2*[0])
    FiredTripleTrgs = array('d',2*[0])
    FiredL1Trgs = array('d',2*[0])
    SeedTime = array('d',2*[0])
    SeedEnergy = array('d',2*[0])
    MIPTotEnergy = array('d',2*[0])
    MIPChi2 = array('d',2*[0])
    MIPSlope = array('d',2*[0])
    MIPIntercept = array('d',2*[0])
    MIPNhitCone = array('d',2*[0])
    MIPIsHalo = array('d',2*[0])
    xtalBits = array('d',2*[0])
    IDbit = array('d',2*[0])
    Scale_stat_up = array('d',2*[0])
    Scale_stat_dn = array('d',2*[0])
    Scale_syst_up = array('d',2*[0])
    Scale_syst_dn = array('d',2*[0])
    Scale_gain_up = array('d',2*[0])
    Scale_gain_dn = array('d',2*[0])
    Resol_rho_up = array('d',2*[0])
    Resol_rho_dn = array('d',2*[0])
    Resol_phi_up = array('d',2*[0])
    Resol_phi_dn = array('d',2*[0])
    outtree.Branch("E", E, 'E[2]/D')
    outtree.Branch("SigmaE",     SigmaE, '    SigmaE[2]/D')
    outtree.Branch("Et", Et, 'Et[2]/D')
    outtree.Branch("Eta", Eta, 'Eta[2]/D')
    outtree.Branch("Phi", Phi, 'Phi[2]/D')
    outtree.Branch("CalibE", CalibE, 'CalibE[2]/D')
    outtree.Branch("CalibEt", CalibEt, 'CalibEt[2]/D')
    outtree.Branch("SCE", SCE, 'SCE[2]/D')
    outtree.Branch("SCRawE", SCRawE, 'SCRawE[2]/D')
    outtree.Branch("ESEnP1", ESEnP1, 'ESEnP1[2]/D')
    outtree.Branch("ESEnP2", ESEnP2, 'ESEnP2[2]/D')
    outtree.Branch("SCEta", SCEta, 'SCEta[2]/D')
    outtree.Branch("SCPhi", SCPhi, 'SCPhi[2]/D')
    outtree.Branch("SCEtaWidth", SCEtaWidth, 'SCEtaWidth[2]/D')
    outtree.Branch("SCPhiWidth", SCPhiWidth, 'SCPhiWidth[2]/D')
    outtree.Branch("SCBrem", SCBrem, 'SCBrem[2]/D')
    outtree.Branch("hasPixelSeed", hasPixelSeed, 'hasPixelSeed[2]/D')
    outtree.Branch("EleVeto", EleVeto, 'EleVeto[2]/D')
    outtree.Branch("R9", R9, 'R9[2]/D')
    outtree.Branch("HoverE", HoverE, 'HoverE[2]/D')
    outtree.Branch("ESEffSigmaRR", ESEffSigmaRR, 'ESEffSigmaRR[2]/D')
    outtree.Branch("SigmaIEtaIEtaFull5x5", SigmaIEtaIEtaFull5x5, 'SigmaIEtaIEtaFull5x5[2]/D')
    outtree.Branch("SigmaIEtaIPhiFull5x5", SigmaIEtaIPhiFull5x5, 'SigmaIEtaIPhiFull5x5[2]/D')
    outtree.Branch("SigmaIPhiIPhiFull5x5", SigmaIPhiIPhiFull5x5, 'SigmaIPhiIPhiFull5x5[2]/D')
    outtree.Branch("E2x2Full5x5", E2x2Full5x5, 'E2x2Full5x5[2]/D')
    outtree.Branch("E5x5Full5x5", E5x5Full5x5, 'E5x5Full5x5[2]/D')
    outtree.Branch("R9Full5x5", R9Full5x5, 'R9Full5x5[2]/D')
    outtree.Branch("PFChIso", PFChIso, 'PFChIso[2]/D')
    outtree.Branch("PFChPVIso", PFChPVIso, 'PFChPVIso[2]/D')
    outtree.Branch("PFPhoIso", PFPhoIso, 'PFPhoIso[2]/D')
    outtree.Branch("PFNeuIso", PFNeuIso, 'PFNeuIso[2]/D')
    outtree.Branch("PFChWorstIso", PFChWorstIso, 'PFChWorstIso[2]/D')
    outtree.Branch("PFChWorstVetoIso", PFChWorstVetoIso, 'PFChWorstVetoIso[2]/D')
    outtree.Branch("EcalPFClusterIso", EcalPFClusterIso, 'EcalPFClusterIso[2]/D')
    outtree.Branch("HcalPFClusterIso", HcalPFClusterIso, 'HcalPFClusterIso[2]/D')
    outtree.Branch("IDMVA", IDMVA, 'IDMVA[2]/D')
    outtree.Branch("FiredSingleTrgs", FiredSingleTrgs, 'FiredSingleTrgs[2]/D')
    outtree.Branch("FiredDoubleTrgs", FiredDoubleTrgs, 'FiredDoubleTrgs[2]/D')
    outtree.Branch("FiredTripleTrgs", FiredTripleTrgs, 'FiredTripleTrgs[2]/D')
    outtree.Branch("FiredL1Trgs", FiredL1Trgs, 'FiredL1Trgs[2]/D')
    outtree.Branch("SeedTime", SeedTime, 'SeedTime[2]/D')
    outtree.Branch("SeedEnergy", SeedEnergy, 'SeedEnergy[2]/D')
    outtree.Branch("MIPTotEnergy", MIPTotEnergy, 'MIPTotEnergy[2]/D')
    outtree.Branch("MIPChi2", MIPChi2, 'MIPChi2[2]/D')
    outtree.Branch("MIPSlope", MIPSlope, 'MIPSlope[2]/D')
    outtree.Branch("MIPIntercept", MIPIntercept, 'MIPIntercept[2]/D')
    outtree.Branch("MIPNhitCone", MIPNhitCone, 'MIPNhitCone[2]/D')
    outtree.Branch("MIPIsHalo", MIPIsHalo, 'MIPIsHalo[2]/D')
    outtree.Branch("xtalBits", xtalBits, 'xtalBits[2]/D')
    outtree.Branch("IDbit", IDbit, 'IDbit[2]/D')
    outtree.Branch("Scale_stat_up", Scale_stat_up, 'Scale_stat_up[2]/D')
    outtree.Branch("Scale_stat_dn", Scale_stat_dn, 'Scale_stat_dn[2]/D')
    outtree.Branch("Scale_syst_up", Scale_syst_up, 'Scale_syst_up[2]/D')
    outtree.Branch("Scale_syst_dn", Scale_syst_dn, 'Scale_syst_dn[2]/D')
    outtree.Branch("Scale_gain_up", Scale_gain_up, 'Scale_gain_up[2]/D')
    outtree.Branch("Scale_gain_dn", Scale_gain_dn, 'Scale_gain_dn[2]/D')
    outtree.Branch("Resol_rho_up", Resol_rho_up, 'Resol_rho_up[2]/D')
    outtree.Branch("Resol_rho_dn", Resol_rho_dn, 'Resol_rho_dn[2]/D')
    outtree.Branch("Resol_phi_up", Resol_phi_up, 'Resol_phi_up[2]/D')
    outtree.Branch("Resol_phi_dn", Resol_phi_dn, 'Resol_phi_dn[2]/D')

    for event_nm in range(ggChain.GetEntries()):
        if (event_nm % args.report == 0):
            logging.info("Processing event {} ".format(event_nm))

        # read parameters
        ggChain.GetEntry(event_nm)

        genpho_list = []
        recopho_list = []

        logging.debug("==========Event {}==========".format(event_nm))
        logging.debug("gen:photon")
        for imc in range(ggChain.nMC):
            if ( ggChain.mcPID[imc] == 22 and ggChain.mcStatusFlag[imc]>>1&1 == 1 ):
                photon = mcPhoton(ggChain, imc)
                logging.debug ("PdgId : {}   pt : {}  eta : {}   phi : {}  mother : {}" .format( photon.PID, photon.Et, photon.Eta, photon.Phi, photon.MomPID ))
                genpho_list.append(photon)

        logging.debug ("    reco:photon")
        # Gen vs Reco matching
        for ireco in range(ggChain.nPho):
            photon = recoPhoton(ggChain,ireco)
            fakeCandidate = True
            for igen, gen in enumerate(genpho_list):
                if ( args.savePrompt == True ):
                    if ( deltaR(gen, photon) < 0.15 and deltaPt(gen, photon) < 0.1 and gen.Eta < 1.4442 and gen.Et > 10 and photon.Et > 10 ):
                        logging.debug("    prompt Candidate : {}  pt : {}  eta : {}   phi : {}" .format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                        photon.matchGenID = igen
                        photon.recoID = ireco
                        photon.setMatchGen(gen)
                        recopho_list.append(photon)
                elif ( args.savePrompt == False ):
                    if ( deltaR(gen, photon) < 0.4 or photon.Et < 10 or photon.Eta > 1.4442):
                        fakeCandidate = False
            if (args.savePrompt == False and fakeCandidate == True):
                logging.debug("    fake Candidate : {}  pt : {}  eta : {}   phi : {}" .format(ireco,ggChain.phoEt[ireco],ggChain.phoEta[ireco],ggChain.phoPhi[ireco]))
                photon.recoID = ireco
                recopho_list.append(photon)


        for i in range(len(recopho_list)):
            for j in range(i+1, len(recopho_list)):
                reco1 = recopho_list[i]
                reco2 = recopho_list[j]
                recoID_list = []
                if ( (args.savePrompt==True and passPromptPairSelection(reco1,reco2)) or (args.savePrompt==False and passFakePairSelection(reco1,reco2)) ):
                    # Place the smaller Et photon at the first index
                    logging.debug("        find photon pair, Et1 = {} Et2 = {} deltaR = {}".format(reco1.Et,reco2.Et,deltaR(reco1,reco2)))
                    if ( reco1.Et < reco2.Et ):
                        recoID_list.append(reco1.recoID)
                        recoID_list.append(reco2.recoID)
                    else:
                        recoID_list.append(reco2.recoID)
                        recoID_list.append(reco1.recoID)
                    # Fill tree
                    for i in range(2):
                        logging.debug("        save photon, Et = {}".format(ggChain.phoEt[recoID_list[i]]))
                        E[i] = ggChain.phoE[recoID_list[i]]
                        SigmaE[i] = ggChain.phoSigmaE[recoID_list[i]]
                        Et[i] = ggChain.phoEt[recoID_list[i]]
                        Eta[i] = ggChain.phoEta[recoID_list[i]]
                        Phi[i] = ggChain.phoPhi[recoID_list[i]]
                        CalibE[i] = ggChain.phoCalibE[recoID_list[i]]
                        CalibEt[i] = ggChain.phoCalibEt[recoID_list[i]]
                        SCE[i] = ggChain.phoSCE[recoID_list[i]]
                        SCRawE[i] = ggChain.phoSCRawE[recoID_list[i]]
                        ESEnP1[i] = ggChain.phoESEnP1[recoID_list[i]]
                        ESEnP2[i] = ggChain.phoESEnP2[recoID_list[i]]
                        SCEta[i] = ggChain.phoSCEta[recoID_list[i]]
                        SCPhi[i] = ggChain.phoSCPhi[recoID_list[i]]
                        SCEtaWidth[i] = ggChain.phoSCEtaWidth[recoID_list[i]]
                        SCPhiWidth[i] = ggChain.phoSCPhiWidth[recoID_list[i]]
                        SCBrem[i] = ggChain.phoSCBrem[recoID_list[i]]
                        hasPixelSeed[i] = ggChain.phohasPixelSeed[recoID_list[i]]
                        EleVeto[i] = ggChain.phoEleVeto[recoID_list[i]]
                        R9[i] = ggChain.phoR9[recoID_list[i]]
                        HoverE[i] = ggChain.phoHoverE[recoID_list[i]]
                        ESEffSigmaRR[i] = ggChain.phoESEffSigmaRR[recoID_list[i]]
                        SigmaIEtaIEtaFull5x5[i] = ggChain.phoSigmaIEtaIEtaFull5x5[recoID_list[i]]
                        SigmaIEtaIPhiFull5x5[i] = ggChain.phoSigmaIEtaIPhiFull5x5[recoID_list[i]]
                        SigmaIPhiIPhiFull5x5[i] = ggChain.phoSigmaIPhiIPhiFull5x5[recoID_list[i]]
                        E2x2Full5x5[i] = ggChain.phoE2x2Full5x5[recoID_list[i]]
                        E5x5Full5x5[i] = ggChain.phoE5x5Full5x5[recoID_list[i]]
                        R9Full5x5[i] = ggChain.phoR9Full5x5[recoID_list[i]]
                        PFChIso[i] = ggChain.phoPFChIso[recoID_list[i]]
                        PFChPVIso[i] = ggChain.phoPFChPVIso[recoID_list[i]]
                        PFPhoIso[i] = ggChain.phoPFPhoIso[recoID_list[i]]
                        PFNeuIso[i] = ggChain.phoPFNeuIso[recoID_list[i]]
                        PFChWorstIso[i] = ggChain.phoPFChWorstIso[recoID_list[i]]
                        PFChWorstVetoIso[i] = ggChain.phoPFChWorstVetoIso[recoID_list[i]]
                        EcalPFClusterIso[i] = ggChain.phoEcalPFClusterIso[recoID_list[i]]
                        HcalPFClusterIso[i] = ggChain.phoHcalPFClusterIso[recoID_list[i]]
                        IDMVA[i] = ggChain.phoIDMVA[recoID_list[i]]
                        FiredSingleTrgs[i] = ggChain.phoFiredSingleTrgs[recoID_list[i]]
                        FiredDoubleTrgs[i] = ggChain.phoFiredDoubleTrgs[recoID_list[i]]
                        FiredTripleTrgs[i] = ggChain.phoFiredTripleTrgs[recoID_list[i]]
                        FiredL1Trgs[i] = ggChain.phoFiredL1Trgs[recoID_list[i]]
                        SeedTime[i] = ggChain.phoSeedTime[recoID_list[i]]
                        SeedEnergy[i] = ggChain.phoSeedEnergy[recoID_list[i]]
                        MIPTotEnergy[i] = ggChain.phoMIPTotEnergy[recoID_list[i]]
                        MIPChi2[i] = ggChain.phoMIPChi2[recoID_list[i]]
                        MIPSlope[i] = ggChain.phoMIPSlope[recoID_list[i]]
                        MIPIntercept[i] = ggChain.phoMIPIntercept[recoID_list[i]]
                        MIPNhitCone[i] = ggChain.phoMIPNhitCone[recoID_list[i]]
                        MIPIsHalo[i] = ggChain.phoMIPIsHalo[recoID_list[i]]
                        IDbit[i] = ggChain.phoIDbit[recoID_list[i]]
                        # xtalBits = ggChain.phoxtalBits[recoID_list[i]]  # Index out of range FIXME
                        # Scale_stat_up = ggChain.phoScale_stat_up[recoID_list[i]]
                        # Scale_stat_dn = ggChain.phoScale_stat_dn[recoID_list[i]]
                        # Scale_syst_up = ggChain.phoScale_syst_up[recoID_list[i]]
                        # Scale_syst_dn = ggChain.phoScale_syst_dn[recoID_list[i]]
                        # Scale_gain_up = ggChain.phoScale_gain_up[recoID_list[i]]
                        # Scale_gain_dn = ggChain.phoScale_gain_dn[recoID_list[i]]
                        # Resol_rho_up = ggChain.phoResol_rho_up[recoID_list[i]]
                        # Resol_rho_dn = ggChain.phoResol_rho_dn[recoID_list[i]]
                        # Resol_phi_up = ggChain.phoResol_phi_up[recoID_list[i]]
                        # Resol_phi_dn = ggChain.phoResol_phi_dn[recoID_list[i]]
                        outtree.Fill()


    outfile.Write("",ROOT.TFile.kOverwrite)
    outfile.Close()
