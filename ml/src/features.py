import numpy

'''
    Features, target, and ancillary variables for the cc1e0pi 
    machine learning-based selection.
'''

# features are complemented by binning information,
# can be useful for quick feature plotting, or for cleaning
# (even though common BDT models won't need it)
FEATURES_CC1E0PI = {
    # overall 
    'slice_vertex_x': numpy.arange(-360, 360+10, 10),
    'slice_vertex_y': numpy.arange(-200, 150+10, 10),
    'slice_vertex_z': numpy.arange(-900, 900+40, 40),
    # Pandora nu score
    'slice_nuscore_nufspfos'         : numpy.arange(0, 15+1, 1), 
    'slice_nuscore_nutothits'        : numpy.arange(0, 7000+200, 200), 
    'slice_nuscore_nuwgdirz'         : numpy.arange(-1, 1+0.02, 0.02), 
    'slice_nuscore_nusps'            : numpy.arange(0, 350+10, 10), 
    'slice_nuscore_nueigen'          : numpy.arange(0, 1+0.02, 0.02), 
    'slice_nuscore_crlongtrkdiry'    : numpy.arange(-1, 1+0.05, 0.05),
    'slice_nuscore_crlongtrkdef'     : numpy.arange(0, 0.4+0.01, 0.01), 
    'slice_nuscore_crlongtrkhitfrac' : numpy.arange(0, 1+0.02, 0.02), 
    'slice_nuscore_crmaxhits'        : numpy.arange(0, 1500+30, 30),
    # NuGraph2-enhanced tagging
    'slice_nugraph2_ng_vtx_hip_hits_ind1'  : numpy.arange(0, 30+1, 1), 
    'slice_nugraph2_ng_vtx_hip_hits_ind2'  : numpy.arange(0, 30+1, 1), 
    'slice_nugraph2_ng_vtx_hip_hits_coll'  : numpy.arange(0, 30+1, 1), 
    'slice_nugraph2_shr_hits_ind1'         : numpy.arange(0, 1000+10, 10), 
    'slice_nugraph2_shr_hits_ind2'         : numpy.arange(0, 1000+10, 10), 
    'slice_nugraph2_shr_hits_coll'         : numpy.arange(0, 1000+10, 10),
    'slice_nugraph2_unclust_shr_hits_ind1' : numpy.arange(0, 500+10, 10), 
    'slice_nugraph2_unclust_shr_hits_ind2' : numpy.arange(0, 500+10, 10), 
    'slice_nugraph2_unclust_shr_hits_coll' : numpy.arange(0, 500+10, 10),
    'slice_nugraph2_nshrpfps'              : numpy.arange(0, 10+1, 1), 
    'slice_nugraph2_nmippfps'              : numpy.arange(0, 10+1, 1), 
    'slice_nugraph2_nhippfps'              : numpy.arange(0, 10+1, 1),
    # NuGraph2-tagged leading shower
    'leading_shr_colldedx'    : numpy.arange(0, 10+0.1, 0.1), 
    'leading_shr_availdedx'   : numpy.arange(0, 10+0.1, 0.1), 
    'leading_shr_openangle'   : numpy.arange(0, 30+0.5, 0.5), 
    'leading_shr_convgap'     : numpy.arange(0, 10+0.2, 0.2),
    'leading_shr_hitshare_bp' : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_trackscore'  : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_shr_frac'    : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_hip_frac'    : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_mip_frac'    : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_mhl_frac'    : numpy.arange(0, 1+0.02, 0.02), 
    'leading_shr_dif_frac'    : numpy.arange(0, 1+0.02, 0.02),
    # NuGraph2-tagged leading proton
    'slice_n_sel_protons'  : numpy.arange(0, 7+1, 1), 
    'leading_p_trackscore' : numpy.arange(0, 1+0.02, 0.02), 
    'leading_p_shr_frac'   : numpy.arange(0, 1+0.02, 0.02), 
    'leading_p_hip_frac'   : numpy.arange(0, 1+0.02, 0.02), 
    'leading_p_mip_frac'   : numpy.arange(0, 1+0.02, 0.02), 
    'leading_p_mhl_frac'   : numpy.arange(0, 1+0.02, 0.02), 
    'leading_p_dif_frac'   : numpy.arange(0, 1+0.02, 0.02),
    # remaining particles
    'slice_n_other' : numpy.arange(0, 7+1, 1),
}

# we won't be working with multiple classes, at least initially
# but additional ground truth information is pretty useful
# for debugging and for physics
TARGET_CC1E0PI = [
    'slice_is_cc1e0pi', 'slice_is_cosmic', 'slice_is_nue', 'slice_is_numu',
    'slice_is_cc', 'slice_is_nc', 'slice_is_pi0', 'slice_is_oofv', 'slice_is_infv',
    'slice_is_qe', 'slice_is_mec', 'slice_is_res', 'slice_is_leading_shr_electron'
]

# more variables - make sure to don't let those slip into 
# the training (truth, energy)
MORE_VARS_CC1E0PI = [
    'run', 'event', 'truthindex', 'nnu', 'E', 'recoE', 'collshrE', 'trueshrE'
]