#! /bin/csh -f
# Thomas Reerink

# Run examples:
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/racmo2.3-greenland-to-im/config-oblimap-racmo2.3-to-im-greenland-10x10km config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-10x10km
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config_oblimap_ccsm_to_im_greenland                           config-files/oblimap/im-to-ccsm/config_oblimap_im_to_ccsm_greenland

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
