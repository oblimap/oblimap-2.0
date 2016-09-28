#! /bin/csh -f
# Thomas Reerink

if($#argv == 1 || $#argv == 3) then

  if($#argv == 1) then
    # Running and comparing with the benchmark results in the directory: oblimap-results/compare/
    # Run like:
    #  ./verify-oblimap-scan-search-block-size.csh config-files/oblimap/ec-earth-to-im/config_oblimap_ec-earth_to_im_greenland

    set config_filename   = $1

    set created_file_name = ` grep '_created_filename_config' ${config_filename} | sed -e 's/.*= .//' -e 's/.$//' `
    set compare_file_name = 'oblimap-results/compare/'` grep '_created_filename_config' ${config_filename} | sed -e 's/.*= .//' -e 's/.$//' -e 's/.*\///' `
    set command           = ` grep ${config_filename} ${config_filename} | sed -e 's/\\!//' `
    set diff_file_name    = 'ddiff.nc'

    time ./${command}

    echo ''
    echo ' The difference can be viewed with:'
    echo '  ncview' ${diff_file_name}
    echo ''

    ncdump -h ${created_file_name} > dump111.h ; ncdump -h ${compare_file_name} > dump222.h ; diff dump111.h dump222.h > diff111222.h ; rm -f dump111.h dump222.h; more diff111222.h; rm -f diff111222.h; diff ${created_file_name} ${compare_file_name}; ncdiff -O ${created_file_name} ${compare_file_name} ${diff_file_name}

  else if($#argv == 3) then
    # Verification of the size of the scan_search_block_size_config: Taking a scan_search_block_size_config = n in a first scanning mapping, and taking it n+1 in
    # a second scanning mapping allows us to proof if n is large enough. If the scanning mappings give identical results, then n was large enough, and maybe a 
    # smaller n will be large enough as well. If the scanning mappings differ, then n must be larger, so try again.
    # This script produces a temporary config file in which the scan_search_block_size_config is adjusted and the fast scanning mode is switched on.
    # Run like:
    #  ./verify-oblimap-scan-search-block-size.csh config-files/oblimap/ec-earth-to-im/config_oblimap_ec-earth_to_im_greenland 17 16

    set config_filename = $1

    set config_filename_1   = ${config_filename}'.tmp-'${2}
    set config_filename_2   = ${config_filename}'.tmp-'${3}
    set created_file_name_1 = ` grep '_created_filename_config' ${config_filename} | sed -e 's/.*= .//' -e 's/.$//' -e 's/.*\///' `'-'${2}
    set created_file_name_2 = ` grep '_created_filename_config' ${config_filename} | sed -e 's/.*= .//' -e 's/.$//' -e 's/.*\///' `'-'${3}

    # Copy the config file and edit the scan_search_block_size_config, level_of_automatic_oblimap_scanning_config and *_created_filename_config:
    sed ${config_filename} -e 's/txt/txt-'${2}'/' -e 's/.*scan_search_block_size_config.*/scan_search_block_size_config                             = '${2}'/' -e 's/.*level_of_automatic_oblimap_scanning_config.*/level_of_automatic_oblimap_scanning_config                = 2/' -e 's/_created_filename_config.*/_created_filename_config                                = "'${created_file_name_1}'" /' > ${config_filename_1}
    # And, for in case they were not present in the config file, adding at the end of the config file the scan_search_block_size_config and level_of_automatic_oblimap_scanning_config
    sed -i ${config_filename_1} -e 's/^\/$/level_of_automatic_oblimap_scanning_config                = 2/'      -e '$ a /'
    sed -i ${config_filename_1} -e 's/^\/$/scan_search_block_size_config                             = '${2}'/' -e '$ a /'
    # Grep the mapping run command from the config header, and exucute the mapping:
    set command = ` grep ${config_filename} ${config_filename} | sed -e 's/\\!//' -e 's/$/.tmp-'${2}'/' `
    time ./${command}

    # Copy the config file and edit the scan_search_block_size_config, level_of_automatic_oblimap_scanning_config and *_created_filename_config:
    sed ${config_filename} -e 's/txt/txt-'${3}'/' -e 's/.*scan_search_block_size_config.*/scan_search_block_size_config                             = '${3}'/' -e 's/.*level_of_automatic_oblimap_scanning_config.*/level_of_automatic_oblimap_scanning_config                = 2/' -e 's/_created_filename_config.*/_created_filename_config                                = "'${created_file_name_2}'" /' > ${config_filename_2}
    # And, for in case they were not present in the config file, adding at the end of the config file the scan_search_block_size_config and level_of_automatic_oblimap_scanning_config
    sed -i ${config_filename_2} -e 's/^\/$/level_of_automatic_oblimap_scanning_config                = 2/'      -e '$ a /'
    sed -i ${config_filename_2} -e 's/^\/$/scan_search_block_size_config                             = '${3}'/' -e '$ a /'
    # Grep the mapping run command from the config header, and exucute the mapping:
    set command = ` grep ${config_filename} ${config_filename} | sed -e 's/\\!//' -e 's/$/.tmp-'${3}'/' `
    time ./${command}
    
   #set diff_file_name = 'ddiff.nc'
    set diff_file_name = 'diff-scan-block-size-'${2}'-'${3}'.nc'
    
    echo ''
    echo ' The difference can be viewed with:'
    echo '  ncview' ${diff_file_name}
    echo ''
    
    ncdump -h ${created_file_name_1} > dump111.h ; ncdump -h ${created_file_name_2} > dump222.h ; diff dump111.h dump222.h > diff111222.h ; rm -f dump111.h dump222.h; more diff111222.h; rm -f diff111222.h; diff ${created_file_name_1} ${created_file_name_2}; ncdiff -O ${created_file_name_1} ${created_file_name_2} ${diff_file_name}

   #mv -f ${config_filename_1} 'config-file-test-'${2}
   #mv -f ${config_filename_2} 'config-file-test-'${3}
    rm -f ${config_filename_1} ${config_filename_2}
   #rm -f ${created_file_name_1} ${created_file_name_2}

  else
   echo ' invalid option'
  endif 

else
 echo ' Needs one argument:'
 echo '  ./verify-oblimap-scan-search-block-size.csh config-files/oblimap/ec-earth-to-im/config_oblimap_ec-earth_to_im_greenland'
 echo ' or needs three arguments:'
 echo '  ./verify-oblimap-scan-search-block-size.csh config-files/oblimap/ec-earth-to-im/config_oblimap_ec-earth_to_im_greenland 15 14'
 echo ' or for comparing the dynamic scan search block size with for example the full search mode:'
 echo '  ./verify-oblimap-scan-search-block-size.csh config-files/oblimap/ec-earth-to-im/config_oblimap_ec-earth_to_im_greenland -1 -3'
endif
