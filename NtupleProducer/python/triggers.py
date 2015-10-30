import FWCore.ParameterSet.Config as cms

#TRIGGERLIST=cms.vstring()
TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### ===================== MC, Spring15 ============================
    ### e mu 
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"),
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), #ele filters
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), # mu filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1"),
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), # ele filters
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"), # muon filters
        channel = cms.int32(0)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu24_eta2p1_v1"),
        path1 = cms.vstring (""), # ele filters
        path2 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # muon filters
        channel = cms.int32(0)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu27_v1"),
        path1 = cms.vstring (""), # ele filters
        path2 = cms.vstring ("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"), # muon filters
        channel = cms.int32(0)
        ),

    ### e tauh
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1"),
        path1 = cms.vstring ("hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # e filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
        path1 = cms.vstring ("hltEle32WP75GsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WP75_Gsf"), # RE-MINIAOD V2
        path1 = cms.vstring ("hltSingleEle22WP75GsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),

  ### mu tauh
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"), #MINIAOD V1
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # mu filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_v1"), #MINIAOD V1
        path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu27_v1"), #MINIAOD V1
        path1 = cms.vstring ("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1"), #RE-MINIAOD V2
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),

### tauh tauh
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters (replicated)
        channel = cms.int32(3)
        ),

### ===================== DATA, Ott 2015 ============================
    ### e mu 
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"), #2015 RUN B/C
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), #ele filters
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), # mu filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2"), #2015 RUN B/C
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), # ele filters
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"), # muon filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3"), #2015 RUN D
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), #ele filters
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), # mu filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3"), #2015 RUN D
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), # ele filters
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"), # muon filters
        channel = cms.int32(0)
        ),

    ### e tauh
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1"), # DATA RUN2015B / C <-- PATH FOR C ?
        path1 = cms.vstring ("hltSingleEle22WPLooseGsfTrackIsoFilter" , "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"), # e filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIso" , "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v1"), # DATARUN2015 B
        path1 = cms.vstring ("hltEle32WP75GsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele23_WPLoose_Gsf_v"), # DATA RUN2015D
        path1 = cms.vstring ("hltEle23WPLooseGsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v2"), # DATA RUN2015 D FALLBACK
        path1 = cms.vstring ("hltEle22WPLooseL1IsoEG20erTau20erGsfTrackIsoFilter" , "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"), # e filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIso" , "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v2"), # DATARUN2015 D FALLBACK
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),


    ### mu tauh
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2"), #2015 RUN B/C
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # mu filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # tauh filters
        channel = cms.int32(2)
        ),
    # cms.PSet (
    #     HLT = cms.string("HLT_IsoMu24_eta2p1_IterTrk02_v2"),
    #     path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # mu filters
    #     path2 = cms.vstring (""), # tauh filters
    #     channel = cms.int32(2)
    #     ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_v2"), #2015 RUN B/C
        path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu18_v"), #2015 RUN D
        path1 = cms.vstring ("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_v1"), #2015 RUN D (FALLBACK)
        path1 = cms.vstring ("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),

    ### tauh tauh
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2"), # RUN 2015 B/C
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters (replicated)
        channel = cms.int32(3)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"), # RUN 2015 B/C
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"), # tauh filters
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"), # tauh filters (replicated)
        channel = cms.int32(3)
        )
    )

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST