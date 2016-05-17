#! /bin/csh -f
# Thomas Reerink

# Run examples:
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_antarctica config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_antarctica
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_greenland  config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_greenland 
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_hemisphere config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_hemisphere
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_himalaya   config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_himalaya  
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/racmo2_clrun-to-im/config_oblimap_racmo2_clrun_to_im_greenland_10x10km config-files/oblimap/im-to-racmo2_clrun/config_oblimap_im_to_racmo2_clrun_greenland_10x10km


if($#argv == 0 || $#argv == 2) then

 set config_file_for_gcm_to_im = 'config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_greenland'
 set config_file_for_im_to_gcm = 'config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_greenland'
 if($#argv == 2) then
  set config_file_for_gcm_to_im = $1
  set config_file_for_im_to_gcm = $2
 endif


 ./src/oblimap_gcm_to_im_program ${config_file_for_gcm_to_im}
 
 # The post profiling:
 if(-e gmon.out) then
  gprof ./src/oblimap_gcm_to_im_program gmon.out > profiling-oblimap_gcm_to_im_program.txt
  rm -f gmon.out
 endif

 ./src/oblimap_im_to_gcm_program ${config_file_for_im_to_gcm}
 
 # The post profiling:
 if(-e gmon.out) then
  gprof ./src/oblimap_im_to_gcm_program gmon.out > profiling-oblimap_im_to_gcm_program.txt
  rm -f gmon.out
 endif


else
 echo ' This script runs without a argument, or requires two OPTIONAL arguments, e.g.:'
 echo '  ./oblimap-to-and-fro-mapping.csh'
 echo ' Or:'
 echo '  ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_greenland config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_greenland'
endif
