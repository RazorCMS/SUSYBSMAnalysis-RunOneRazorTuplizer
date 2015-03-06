SUSYBSMAnalysis-RunOneRazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 1 AOD

Instructions for compiling in CMSSW
--------------

    cmsrel CMSSW_5_3_26
    cd CMSSW_5_3_26/src
    cmsenv
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RunOneRazorTuplizer SUSYBSMAnalysis/RunOneRazorTuplizer
    scram b
    
Running the BASE ntuplizer
--------------

    cmsRun python/razorTuplizer_MC.py
    cmsRun python/razorTuplizer_data.py

Use the appropriate python config file for data or MC.    
Before running, check python/razorTuplizer.py to make sure that the correct global tag is defined. (process.GlobalTag.globaltag = ...)

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py
