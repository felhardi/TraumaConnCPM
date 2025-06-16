import os
import cbig_network_correspondence as cnc

config_path = os.path.expanduser('/gpfs/milgram/project/gee_dylan/fah29/cbig/ICAconfig_sh.ini')
file_path = os.path.expanduser('/gpfs/milgram/project/gee_dylan/fah29/cbig/combined_posmask_sh_bin_MNI.nii.gz')

# Mattew Glasser2016 360-ROI with Ji2019 12 Cole-Anticevic networks
# HCP 25-node ICA maps
atlas_names_list = ["MG360J12", "HCPICA"]

ref_params = cnc.compute_overlap_with_atlases.DataParams(config_path, file_path)

cnc.compute_overlap_with_atlases.network_correspondence(
    ref_params,
    atlas_names_list,
    os.path.expanduser('/gpfs/milgram/project/gee_dylan/fah29/cbig/singleoverlap_results_sh')
)
