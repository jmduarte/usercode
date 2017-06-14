dataset = '/GluGluToRSGravitonToHHTo4B_M-450_narrow_13TeV-madgraph/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/GEN-SIM-RAW'

from WMCore.Configuration import Configuration
config = Configuration()

name = 'HLT_Ntuple_Hbb_Signal_v0p1_localTest'
config.section_("General")
config.General.workArea = 'crab_'+name
config.General.transferLogs=True
#        config.General.requestName = name+"_"+dataset.replace('/',"_")
config.General.requestName = name+"_"+dataset.split('/')[1]

config.section_("JobType")
config.JobType.numCores = 4
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'crab_fake_pset.py'
config.JobType.scriptExe = 'crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
os.system("voms-proxy-info -path | xargs -i  cp {}  .")
config.JobType.inputFiles = [
                             'fwlite_config.py',
                             'hltForNtuples3_dump.py',
                             'script.py',
                             'python.tar.gz',
                             'utils.py',
                             'crab_config.py',
                             '.x509up_u636',
]

config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 ##FIXME: use 20

#        config.Data.totalUnits = -1
config.Data.totalUnits = 1

#        config.Data.totalUnits = 100*config.Data.unitsPerJob ##FIXME: use -1
#        config.Data.unitsPerJob = 1 #FIXME: use 20
#        config.Data.totalUnits = 1 #FIXME: use -1
config.Data.outLFNDirBase = '/store/user/woodson/' + name
config.Data.publication = True
#        config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
#        config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt' 
#'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.inputDataset = dataset
if "AODSIM" in dataset:
    config.Data.secondaryInputDataset = dataset.replace("AODSIM","GEN-SIM-RAW")

#        config.Data.publishDataName = config.General.requestName
config.Data.outputDatasetTag = config.General.requestName
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.blacklist = ['T0_*']
config.Site.storageSite = "T2_CH_CSCS"

print "submitting ",dataset

