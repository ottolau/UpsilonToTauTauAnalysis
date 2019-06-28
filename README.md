# UpsilonToTauTauAnalysis

##### To work with CMSSW_10_2_14 and head version, you do :
```
cmsrel CMSSW_10_2_14
cd CMSSW_10_2_14
cmsenv
git clone https://github.com/ottolau/UpsilonToTauTauAnalysis.git
scram b clean; scram b -j20
```

To run the analyzer, in UpsilonToTauTauAnalysis/UpsilonToTauTauAnalyzer/test/, you do:
```
cmsRun run_mc2018_102X.py
```

