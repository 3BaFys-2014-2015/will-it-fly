#change path to correct path of the airfoils

PATH=/home/uauser/Shared/coord_seligFmt.tar/coord_seligFmt/*

source set_env.sh
for file in $PATH 
do
./build/bin/wif_demo_core_airfoil $file
done
source unset_env.sh
