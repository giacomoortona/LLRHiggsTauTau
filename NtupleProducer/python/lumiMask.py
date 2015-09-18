# lumi masks applied when running on data
# produce list with ConvertJSON.py in test/tools

#used by crab
JSONFILE = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-254349_13TeV_PromptReco_Collisions15_JSON_v2.txt"

#used in local productions
LUMIMASK = cms.untracked.VLuminosityBlockRange( *(
      '251244:85-251244:86',
      '251244:88-251244:93',
      '251244:96-251244:121',
      '251244:123-251244:156',
      '251244:158-251244:428',
      '251244:430-251244:442',
      '251251:1-251251:31',
      '251251:33-251251:97',
      '251251:99-251251:167',
      '251252:1-251252:283',
      '251252:285-251252:505',
      '251252:507-251252:554',
      '251561:1-251561:94',
      '251562:1-251562:439',
      '251562:443-251562:691',
      '251643:1-251643:216',
      '251643:222-251643:606',
      '251721:21-251721:36',
      '251721:123-251721:244',
      '251883:56-251883:56',
      '251883:58-251883:60',
      '251883:62-251883:144',
      '251883:156-251883:437',
))
