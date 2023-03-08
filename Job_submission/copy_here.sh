#!/bin/bash
for i in {61..200}
do
gfal-copy gsiftp://dcache-cms-gridftp.desy.de/pnfs//desy.de/cms/tier2/store/group/upgrade/delphes_output/YR_Delphes/Delphes342pre15_hadd/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_200PU/TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_200PU_${i}.root  gsiftp://cms-gridftp.rcac.purdue.edu/store/user/abakshi/TTBar_Delphes/
done
