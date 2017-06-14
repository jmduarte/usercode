#!/usr/bin/python
import ROOT
import itertools
import resource
from array import array
from math import sqrt, pi, log10, log, exp
# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
from utils import deltaR,SetVariable,DummyClass,productWithCheck,checkTriggerIndex

#ROOT.gROOT.LoadMacro("/scratch/sdonato/NtupleForPaolo/CMSSW_8_0_3_patch1/src/DataFormats/L1Trigger/interface/EtSumHelper.h")

Handle.productWithCheck = productWithCheck

maxJets         = 50
bunchCrossing   = 0
pt_min          = 20

def FillVector(source,variables,minPt=pt_min):
    variables.num[0] = 0
    for obj in source.productWithCheck():
        if obj.pt()<minPt: continue
        if variables.num[0]<len(variables.pt):
            for (name,var) in variables.__dict__.items():
                if name == "pt" :           var[variables.num[0]] = obj.pt()
                elif name == "eta" :        var[variables.num[0]] = obj.eta()
                elif name == "phi" :        var[variables.num[0]] = obj.phi()
                elif name == "mass" :       var[variables.num[0]] = obj.mass()
                elif name == "neHadEF" :    var[variables.num[0]] = obj.neutralHadronEnergyFraction()
                elif name == "neEmEF" :     var[variables.num[0]] = obj.neutralEmEnergyFraction()
                elif name == "chHadEF" :    var[variables.num[0]] = obj.chargedHadronEnergyFraction()
                elif name == "chEmEF" :     var[variables.num[0]] = obj.chargedEmEnergyFraction()
                elif name == "muEF" :       var[variables.num[0]] = obj.muonEnergyFraction()
                elif name == "mult" :       var[variables.num[0]] = obj.chargedMultiplicity()+obj.neutralMultiplicity();
                elif name == "neMult" :     var[variables.num[0]] = obj.neutralMultiplicity()
                elif name == "chMult" :     var[variables.num[0]] = obj.chargedMultiplicity()
            variables.num[0] += 1

def FillBtag(btags_source, jets, jet_btags):
    for i in range(jets.num[0]):
        jet_btags[i] = -20.
        dRmax = 0.3
        btags = btags_source.productWithCheck()
        for j in range(0,btags.size()):
            jetB = btags.key(j).get()
            dR = deltaR(jetB.eta(),jetB.phi(),jets.eta[i],jets.phi[i])
            if dR<dRmax:
                jet_btags[i] = max(0.,btags.value(j))
                dRmax = dR

def Matching(phi, eta, jets):
    index = -1
    for i in range(jets.num[0]):
        dRmax = 0.3
        dR = deltaR(eta,phi,jets.eta[i],jets.phi[i])
        if dR<dRmax:
            index = i
            dRmax = dR
    return index

def FillVectorAndBtag(jets,btags=0,Jets_csv=0):
    jet.num[0] = 0
    for jet in jets.productWithCheck():
        if jet.pt()<pt_min: continue
        if jet.num[0]<len(Jets_pt):                
            for (name,obj) in v.__dict__.keys():
                if name == "pt" :           obj[jet.num[0]] = jet.pt()
                elif name == "eta" :        obj[jet.num[0]] = jet.eta()
                elif name == "phi" :        obj[jet.num[0]] = jet.phi()
                elif name == "mass" :       obj[jet.num[0]] = jet.mass()
                elif name == "neHadEF" :    obj[jet.num[0]] = jet.neutralHadronEnergyFraction()
                elif name == "neEmEF" :     obj[jet.num[0]] = jet.neutralEmEnergyFraction()
                elif name == "chHadEF" :    obj[jet.num[0]] = jet.chargedHadronEnergyFraction()
                elif name == "chEmEF" :     obj[jet.num[0]] = jet.chargedEmEnergyFraction()
                elif name == "muEF" :       obj[jet.num[0]] = jet.muonEnergyFraction()
                elif name == "mult" :       obj[jet.num[0]] = jet.chargedMultiplicity()+jet.neutralMultiplicity();
                elif name == "neMult" :     obj[jet.num[0]] = jet.neutralMultiplicity()
                elif name == "chMult" :     obj[jet.num[0]] = jet.chargedMultiplicity()
            
            offlineCSV = -2.
            if not btags is 0:
                for j in range(0,btags.productWithCheck().size()):
                    jetB = btags.productWithCheck().key(j).get()
                    dR = deltaR(jetB.eta(),jetB.phi(),jet.eta(),jet.phi())
                    if dR<0.3:
                        offlineCSV = max(0.,btags.productWithCheck().value(j))
                        break
                Jets_csv[jet.num[0]] = offlineCSV
            jet.num[0] += 1

    #def BookVariable(tree,name="variable",type_='F'):
    #    var   = SetVariable(tree,name ,type_)
    #    return var

def getVertex(vertex_source):
    vertices = vertex_source.productWithCheck()
    if vertices.size()>0:
        return vertices.at(0).z()
    else:
        return -1000

def WithFallback(product,method="pt"):
    if product.size()>0:
        return getattr(product[0],method)()
    else:
        return -10

def BookVector(tree,name="vector",listMembers=[]):
    obj = DummyClass()
    obj.num   = SetVariable(tree,name+'_num' ,'I')
    for member in listMembers:
        if "match" in name:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'I',name+'_num',maxJets))
        else:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'F',name+'_num',maxJets))
    return obj

    ##########################################################################

def launchNtupleFromHLT(fileOutput,filesInput, secondaryFiles, maxEvents,preProcessing=True, firstEvent=0):
    bunchCrossing   = 12
    print "filesInput: ",filesInput
    print "fileOutput: ",fileOutput
    print "secondaryFiles: ",secondaryFiles
    print "maxEvents: ",maxEvents
    print "preProcessing: ",preProcessing
    print "firstEvent: ",firstEvent
    
    Signal = True
    if len(filesInput)>0 and ('QCD' in filesInput[0]):
        Signal = False
    print "Signal=",Signal
    
    ## Pre-processing
    if preProcessing:
        from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
        from PhysicsTools.HeppyCore.framework.config import MCComponent
        cmsRun_config = "hltForNtuples3_dump.py"
        preprocessor = CmsswPreprocessor(cmsRun_config)
        cfg = MCComponent("OutputHLT",filesInput, secondaryfiles=secondaryFiles)
        print "Run cmsswPreProcessing using:"
        print cfg.name
        print cfg.files
        print cfg.secondaryfiles
        print
        try:
            preprocessor.run(cfg,".",firstEvent,maxEvents)
        except:
            print "cmsswPreProcessing failed!"
            print "cat cmsRun_config.py"
            config = file(cmsRun_config)
            print config.read()
            print "cat cmsRun.log"
            log = file("cmsRun.log")
            print log.read()
            preprocessor.run(cfg,".",firstEvent,maxEvents)
            raise Exception("CMSSW preprocessor failed!")
    
    f = ROOT.TFile(fileOutput,"recreate")
    tree = ROOT.TTree("tree","tree")

    fwLiteInputs = ["cmsswPreProcessing.root"]
    if len(filesInput)==0: exit
    import os.path
    if not os.path.isfile(fwLiteInputs[0]):
        raise Exception( fwLiteInputs[0] + " does not exist.")
    events = Events (fwLiteInputs)

    ### list of input variables ###

    triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults::TEST")

    pileUp_source, pileUp_label = Handle("vector<PileupSummaryInfo>"), ("addPileupInfo")

    generator_source, generator_label       = Handle("GenEventInfoProduct"), ("generator")

    l1HT_source, l1HT_label                 = Handle("BXVector<l1t::EtSum>"), ("hltGtStage2Digis","EtSum")
    l1Jets_source, l1Jets_label               = Handle("BXVector<l1t::Jet>"), ("hltGtStage2Digis","Jet")

    offMet_source, offMet_label             = Handle("vector<reco::PFMET>"), ("pfMet")
    genMet_source, genMet_label             = Handle("vector<reco::GenMET>"), ("genMetTrue")

    pfMet_source, pfMet_label               = Handle("vector<reco::PFMET>"), ("hltPFMETProducer")
    pfMht_source, pfMht_label               = Handle("vector<reco::MET>"), ("hltPFMHTTightID")

    caloMet_source, caloMet_label           = Handle("vector<reco::CaloMET>"), ("hltMet")
    caloMht_source, caloMht_label           = Handle("vector<reco::MET>"), ("hltMht")
    caloMhtNoPU_source, caloMhtNoPU_label   = Handle("vector<reco::MET>"), ("hltMHTNoPU")

    caloJets_source, caloJets_label           = Handle("vector<reco::CaloJet>"), ("hltAK4CaloJetsCorrectedIDPassed")

    calobtag_source, calobtag_label         = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltCombinedSecondaryVertexBJetTagsCalo")

    caloPUid_source, caloPUid_label         = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltCaloJetFromPV")

    pfJets_source, pfJets_label               = Handle("vector<reco::PFJet>"), ("hltAK4PFJetsLooseIDCorrected")
    pfbtag_source, pfbtag_label             = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltCombinedSecondaryVertexBJetTagsPF")

    if Signal:
        genJets_source, genJets_label         = Handle("vector<reco::GenJet>"), ("ak4GenJetsNoNu")
        offJets_source, offJets_label         = Handle("vector<reco::PFJet>"), ("ak4PFJets")
        offbtag_source, offbtag_label         = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        
        FastPrimaryVertex_source, FastPrimaryVertex_label   = Handle("vector<reco::Vertex>"), ("hltFastPrimaryVertex")
        FPVPixelVertices_source, FPVPixelVertices_label     = Handle("vector<reco::Vertex>"), ("hltFastPVPixelVertices")
        PixelVertices_source, PixelVertices_label           = Handle("vector<reco::Vertex>"), ("hltPixelVertices")
        VerticesPF_source, VerticesPF_label                 = Handle("vector<reco::Vertex>"), ("hltVerticesPF")
        VerticesL3_source, VerticesL3_label                 = Handle("vector<reco::Vertex>"), ("hltVerticesL3")
        
        genParticles_source, genParticles_label             = Handle("vector<reco::GenParticle>"), ("genParticles")
    
    ### create output variables ###
    
    caloJets     = BookVector(tree,"caloJets",['pt','eta','phi','mass','matchOff','matchGen','puId','csv'])
    pfJets       = BookVector(tree,"pfJets",['pt','eta','phi','mass','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv'])
    
    l1Jets       = BookVector(tree,"l1Jets",['pt','eta','phi','matchOff','matchGen'])
    l1HT        = SetVariable(tree,'l1HT')
    caloMet       = BookVector(tree,"caloMet",['pt','phi'])
    caloMht       = BookVector(tree,"caloMht",['pt','phi'])
    caloMhtNoPU   = BookVector(tree,"caloMhtNoPU",['pt','phi'])
    pfMet         = BookVector(tree,"pfMet",['pt','phi'])
    pfMht         = BookVector(tree,"pfMht",['pt','phi'])
    l1Met       = SetVariable(tree,'l1Met')
    l1Met_phi   = SetVariable(tree,'l1Met_phi')
    l1Mht       = SetVariable(tree,'l1Mht')
    l1Mht_phi   = SetVariable(tree,'l1Mht')
    if Signal:
        genJets  = BookVector(tree,"genJets",['pt','eta','phi','mass','mcFlavour','mcPt'])
        offJets  = BookVector(tree,"offJets",['pt','eta','phi','mass','csv','matchGen'])
        offMet  = BookVector(tree,"offMet",['pt','phi'])
        genMet  = BookVector(tree,"genMet",['pt','phi'])
        FastPrimaryVertex   = SetVariable(tree,'FastPrimaryVertex')
        FPVPixelVertices    = SetVariable(tree,'FPVPixelVertices')
        PixelVertices       = SetVariable(tree,'PixelVertices')
        VerticesPF          = SetVariable(tree,'VerticesPF')
        VerticesL3          = SetVariable(tree,'VerticesL3')
        trueVertex          = SetVariable(tree,'trueVertex')

    evt         = SetVariable(tree,'evt')
    lumi        = SetVariable(tree,'lumi')
    run         = SetVariable(tree,'run')
    
    pu          = SetVariable(tree,'pu')
    ptHat       = SetVariable(tree,'ptHat')
    maxPUptHat  = SetVariable(tree,'maxPUptHat')
    
    f.cd()
    
    ##get trigger names
    events.to(0)
    for event in events: break
    event.getByLabel(triggerBitLabel, triggerBits)
    names = event.object().triggerNames(triggerBits.product())
    triggerNames = names.triggerNames()
    for name in triggerNames: name = name.split("_v")[0]
    nTriggers = len(triggerNames)
    triggerVars = {}
    for trigger in triggerNames:
        triggerVars[trigger]=array( 'i', [ 0 ] )
        tree.Branch( trigger, triggerVars[trigger], trigger+'/O' )

    ##event loop
    for iev,event in enumerate(events):
        if iev>maxEvents and maxEvents>=0: break
        event.getByLabel(triggerBitLabel, triggerBits)
        event.getByLabel(generator_label, generator_source)
        event.getByLabel(pileUp_label, pileUp_source)
        event.getByLabel(caloJets_label, caloJets_source)
        event.getByLabel(calobtag_label, calobtag_source)
        event.getByLabel(caloPUid_label, caloPUid_source)
        event.getByLabel(caloMet_label, caloMet_source)
        event.getByLabel(caloMht_label, caloMht_source)
        event.getByLabel(caloMhtNoPU_label, caloMhtNoPU_source)
        event.getByLabel(pfMet_label, pfMet_source)
        event.getByLabel(pfMht_label, pfMht_source)
        event.getByLabel(pfJets_label, pfJets_source)
        event.getByLabel(pfbtag_label, pfbtag_source)
        event.getByLabel(l1Jets_label, l1Jets_source)
        event.getByLabel(l1HT_label, l1HT_source)
        
        ####################################################
        

        run[0]          = event.eventAuxiliary().run()
        lumi[0]         = event.eventAuxiliary().luminosityBlock()
        evt[0]          = event.eventAuxiliary().event()
        
        if Signal:
            event.getByLabel(offMet_label, offMet_source)
            event.getByLabel(genMet_label, genMet_source)
            event.getByLabel(offJets_label, offJets_source)
            event.getByLabel(genJets_label, genJets_source)
            event.getByLabel(offbtag_label, offbtag_source)
            event.getByLabel(FastPrimaryVertex_label, FastPrimaryVertex_source)
            event.getByLabel(FPVPixelVertices_label, FPVPixelVertices_source)
            event.getByLabel(PixelVertices_label, PixelVertices_source)
            event.getByLabel(VerticesPF_label, VerticesPF_source)
            event.getByLabel(VerticesL3_label, VerticesL3_source)
            event.getByLabel(genParticles_label, genParticles_source)
            
            FastPrimaryVertex[0] = getVertex(FastPrimaryVertex_source)
            FPVPixelVertices[0] = getVertex(FPVPixelVertices_source)
            PixelVertices[0] = getVertex(PixelVertices_source)
            VerticesPF[0] = getVertex(VerticesPF_source)
            VerticesL3[0] = getVertex(VerticesL3_source)
            trueVertex[0] = genParticles_source.productWithCheck().at(2).vertex().z()
        
        FillVector(caloJets_source,caloJets)
        FillVector(pfJets_source,pfJets)
        FillVector(l1Jets_source,l1Jets)
        
        FillVector(caloMet_source,caloMet, 0)
        FillVector(caloMht_source,caloMht, 0)
        FillVector(caloMhtNoPU_source,caloMhtNoPU, 0)
        FillVector(pfMet_source,pfMet, 0)
        FillVector(pfMht_source,pfMht, 0)
        
        FillBtag(calobtag_source, caloJets, caloJets.csv)
        FillBtag(caloPUid_source, caloJets, caloJets.puId)
        FillBtag(pfbtag_source, pfJets, pfJets.csv)
        
        l1Met[0],l1Met_phi[0],l1Mht[0],l1Mht_phi[0],l1HT[0] = -1,-1,-1,-1,-1
        for et in l1HT_source.productWithCheck():
            if et.getType()==ROOT.l1t.EtSum.kMissingEt:
                (l1Met[0],l1Met_phi[0]) = (et.et(),et.phi())
            elif et.getType()==ROOT.l1t.EtSum.kMissingHt:
                (l1Mht[0],l1Mht_phi[0]) = (et.et(),et.phi())
            elif et.getType()==ROOT.l1t.EtSum.kTotalEt:
                pass
            elif et.getType()==ROOT.l1t.EtSum.kTotalHt:
                l1HT[0] = et.et()
        
        if Signal:
            FillVector(offMet_source,offMet, 0)
            FillVector(genMet_source,genMet, 0)
            
            FillVector(offJets_source,offJets,15)
            FillVector(genJets_source,genJets,15)
            
            FillBtag(offbtag_source, offJets, offJets.csv)
            
            for i in range(caloJets.num[0]):
                caloJets.matchOff[i] = Matching(caloJets.phi[i],caloJets.eta[i],offJets)
                caloJets.matchGen[i] = Matching(caloJets.phi[i],caloJets.eta[i],genJets)
            
            for i in range(pfJets.num[0]):
                pfJets.matchOff[i] = Matching(pfJets.phi[i],pfJets.eta[i],offJets)
                pfJets.matchGen[i] = Matching(pfJets.phi[i],pfJets.eta[i],genJets)
            
            for i in range(l1Jets.num[0]):
                l1Jets.matchOff[i] = Matching(l1Jets.phi[i],l1Jets.eta[i],offJets)
                l1Jets.matchGen[i] = Matching(l1Jets.phi[i],l1Jets.eta[i],genJets)
            
            for i in range(offJets.num[0]):
                offJets.matchGen[i] = Matching(offJets.phi[i],offJets.eta[i],genJets)
            
            for i in range(genJets.num[0]):
                genJets.mcFlavour[i] = -100
                genJets.mcPt[i] = -100
            
            for genParticle in genParticles_source.productWithCheck():
                if genParticle.pt()<5: continue
                if not (abs(genParticle.pdgId()) in [21,1,2,3,4,5,11,13,15]): continue
                if genParticle.mother().pt()>5 and (abs(genParticle.mother().pdgId()) in [21,1,2,3,4,5,11,13]): continue
                if evt[0]==7826939:
                    print "genParticle:"
                    print genParticle.pt(),genParticle.eta(),genParticle.phi(),genParticle.pdgId()
                    print "genJets:"
                for i in range(genJets.num[0]):
                    if genParticle.pt()<0.2*genJets.pt[i]: continue
                    if deltaR(genParticle.eta(),genParticle.phi(),genJets.eta[i],genJets.phi[i])<0.4:
                        if evt[0]==7826939:
                            print genJets.pt[i],genJets.eta[i],genJets.eta[i],genJets.mcFlavour[i]
                            print "not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1]):",not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1])
                        if abs(genParticle.pdgId())==5:
                            if genJets.mcFlavour[i]!=5 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==4 and not int(abs(genJets.mcFlavour[i])) in [5]:
                            if genJets.mcFlavour[i]!=4 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==3 and not int(abs(genJets.mcFlavour[i])) in [5,4]:
                            if genJets.mcFlavour[i]!=3 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==2 and not int(abs(genJets.mcFlavour[i])) in [5,4,3]:
                            if genJets.mcFlavour[i]!=2 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==1 and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2]:
                            if genJets.mcFlavour[i]!=1 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==21 and not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1]):
                            if genJets.mcFlavour[i]!=21 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [11,13] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21]:
                            if not (genJets.mcFlavour[i] in [11,13]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [15] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21,11,13]:
                            if not (genJets.mcFlavour[i] in [15]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [22] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21,11,13,15]:
                            if not (genJets.mcFlavour[i] in [22]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        if evt[0]==7826939:
                            print "newFlav:",genJets.mcFlavour[i]
                        
        ptHat[0]    = generator_source.product().qScale()
        
        if bunchCrossing>=pileUp_source.productWithCheck().size() or pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing()!=0:
            found=False
            for bunchCrossing in range(pileUp_source.productWithCheck().size()):
                if pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing() == 0 :
                    found=True;
                    break
            if not found:
                Exception("Check pileupSummaryInfos!")
            print "I'm using bunchCrossing=",bunchCrossing
        pu[0] = pileUp_source.productWithCheck().at(bunchCrossing).getTrueNumInteractions()
        
        maxPUptHat[0] = -1
        for ptHat_ in pileUp_source.productWithCheck().at(bunchCrossing).getPU_pT_hats():
            maxPUptHat[0] = max(maxPUptHat[0],ptHat_)
        
        names = event.object().triggerNames(triggerBits.product())
        for i,triggerName in enumerate(triggerNames):
            index = names.triggerIndex(triggerName)
#            print "index=",index
            if checkTriggerIndex(triggerName,index,names.triggerNames()):
                triggerVars[triggerName][0] = triggerBits.product().accept(index)
#                print "acc:",triggerBits.product().accept(index)
            else:
                triggerVars[triggerName][0] = 0
        
        if iev%10==1: print "Event: ",iev," done."
        tree.Fill()

    f.Write()
    f.Close()

if __name__ == "__main__":
    #filesInput = ["root://eoscms.cern.ch//eos/cms/store/mc/PhaseIFall16DR/GluGluToRSGravitonToHHTo4B_M-450_narrow_13TeV-madgraph/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/80000/F8161EEB-9810-E711-A85C-FA163E0B564E.root"]
    filesInput = ["root://eoscms.cern.ch//eos/cms/store/mc/PhaseIFall16DR/GluGluToRSGravitonToHHTo4B_M-450_narrow_13TeV-madgraph/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/80000/EAACC9D6-0F11-E711-A9B1-FA163EDAFEAB.root"]
#    filesInput = ["file:/mnt/t3nfs01/data01/shome/sdonato/QCD470_GEN-SIM-RAW_PhaseI_83X_FlatPU28to62.root"]
    secondaryFiles = []
    fileOutput = "tree.root"
    maxEvents = 100
    launchNtupleFromHLT(fileOutput,filesInput,secondaryFiles,maxEvents, preProcessing=False)
