rm costMERGE.f
rm resMERGE.f
rm jMERGE.f
rm dudaMERGE.f


cat cost.f griddata.f > costMERGE.f
cat jRES.f griddata.f boundaries.f > resMERGE.f

cat boundaries.f cost.f griddata.f initialise.f nozzle_main.f resid.f stagedata.f timestep.f update.f > jMERGE.f
cat boundaries.f cost.f griddata.f initialise.f nozzle_main.f resid.f stagedata.f timestep.f update.f > dudaMERGE.f




#find . -maxdepth 1 -iname '*.f' -not -name 'wrapper.f' -exec cat {} +>allMERGE.f
##find . -maxdepth 1 -iname '*.f' -not -name 'wrapper.f' -o -not -name 'costMERGE_d-all.f' -exec cat {} +>allMERGE.f
#cp allMERGE.f ./jMERGE.f
#mv allMERGE.f ./dudaMERGE.f
