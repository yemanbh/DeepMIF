rem To use this pipeline make sure you copied all your the code and the docker to a folder where you install deepmif
rem this are the input arguments
set docker_path=%1
set input_dir=%2
set output_dir=%3
set distance=%4
set cell_phenotypes=%5
set nuclear_markers=%6
set non_nuclear_markers=%7

rem set code_dir=C:\Users\yhagos\Dropbox (ICR)\Yeman\Projects\DeepMIFUI\v5_docker
rem load docker image
rem docker load -i %docker_path%

rem run deep mif pipeline docker; src files inside docker
docker run ^
   -v "%input_dir%:/input" ^
   -v "%output_dir%:/output" ^
    dmifdi:01 ^
    -distance %distance% -cell_phenotypes %cell_phenotypes% ^
    -nuclear_m %nuclear_markers% ^
    -non_nuclear_m %non_nuclear_markers%

