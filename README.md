# Data processing pipeline for mcLARO T2&T2* mapping of liver before and after iron overload

## Step 1: save k-space data from scanarxiv (MATLAB)
Run `kdata  = myScanArchiveRead_T1wT2wQSM('Path_to_Scanarxiv', 128, 4);` to load and reorganize k-space data from scanarxiv.\
Then run `writecfl('Your_kdata_filename', kdata);` to save k-space data.

## Step 2: reconstruct multi-echo multi-contrast images from saved k-space data (Python)
You need to add
`export TOOLBOX_PATH=/data/Jiahao/src/bart`
`export PATH=${TOOLBOX_PATH}:${PATH}`
`export PYTHONPATH="${TOOLBOX_PATH}/python:$PYTHONPATH"`
to your .bashrc for bart recon path in Python.\
Then run `python ./recon.py Your_kdata_filename` to reconstruct multi-echo multi-contrast images.

## Step 3: dictionary matching for T2 mapping from multi-contrast images (MATLAB)
Run `T1T2_mapping.m` to generate `qMaps_all.mat` (a 4 dimensional array with T1 map, T2 map and proton density concatenated into the 4th dimension) from your reconstructed images. Make sure to change filename at line #4 before run this script.

------
You may run the same pipeline for both pre and post contrast injection data.\
I have run the pipeline for both data and saved them as "qMaps_all_pre.mat" and "qMaps_all_post.mat" in this repo. You can do it yourself in the future.
