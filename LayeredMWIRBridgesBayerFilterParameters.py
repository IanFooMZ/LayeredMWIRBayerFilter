#
# Parameter file for the Bayer Filter CMOS optimization
#

import numpy as np

#
# Files
#
# project_name = 'layered_mwir_2d_lithography_bridges_rgb_10layers_3p5to5p5um_si_fixed_step_addcage_25x25x25um_f30um'
# project_name = 'layered_mwir_2d_lithography_bridges_rgb_10layers_3p5to5p5um_si_fixed_step_addcage_25x25x25um_f30um'
# project_name = 'layered_mwir_2d_lithography_bridges_rgb_10layers_3p5to5p5um_si_fixed_step_addcage_30x30x25um_f30um'
project_name = 'angInc_Cu450nm_0deg_src15_layered_mwir_2dlit_bridges_rgb_10layers_3p5to5p5um_si_fixed_step_addcage_30x30x25um_f30um'

# todo(gdrobert): consider the contrast here like you are doing in the cmos designs.. minimize energy into wrong quadrant

#
# Optical
#
background_index = 1.0
min_device_index = 1.0
max_device_index = 1.5

min_device_permittivity = min_device_index**2
max_device_permittivity = max_device_index**2

init_permittivity_0_1_scale = 0.5

focal_length_um = 30
focal_plane_center_lateral_um = 0
focal_plane_center_vertical_um = -focal_length_um

#
# Device
#
mesh_spacing_um = 0.25

device_size_lateral_um = 30#25
device_size_verical_um = 25

device_voxels_lateral = 1 + int(device_size_lateral_um / mesh_spacing_um)
# todo: fix!
# device_voxels_vertical = 1 + int(device_size_verical_um / mesh_spacing_um)
device_voxels_vertical = 2 + int(device_size_verical_um / mesh_spacing_um)

device_vertical_maximum_um = device_size_verical_um
device_vertical_minimum_um = 0

silicon_thickness_um = 3.5

num_sidewalls = 4
sidewall_thickness_um = 0.45
sidewall_material = 'Cu (Copper) - Palik'
sidewall_extend_focalplane = False
sidewall_x_positions_um = [device_size_lateral_um / 2 + sidewall_thickness_um / 2, 0, -device_size_lateral_um / 2 - sidewall_thickness_um / 2, 0]
sidewall_y_positions_um = [0, device_size_lateral_um / 2 + sidewall_thickness_um / 2, 0, -device_size_lateral_um / 2 - sidewall_thickness_um / 2]
sidewall_xspan_positions_um = [sidewall_thickness_um, device_size_lateral_um + sidewall_thickness_um * 2,
                               sidewall_thickness_um, device_size_lateral_um + sidewall_thickness_um * 2]
sidewall_yspan_positions_um = [device_size_lateral_um + sidewall_thickness_um * 2, sidewall_thickness_um, 
                               device_size_lateral_um + sidewall_thickness_um * 2, sidewall_thickness_um]


#
# Spectral
#
lambda_min_um = 3.5
lambda_max_um = 5.5
bandwidth_um = lambda_max_um - lambda_min_um

num_bands = 3
num_points_per_band = 10
num_design_frequency_points = num_bands * num_points_per_band

lambda_values_um = np.linspace(lambda_min_um, lambda_max_um, num_design_frequency_points)
max_intensity_by_wavelength = (device_size_lateral_um**2)**2 / (focal_length_um**2 * lambda_values_um**2)


num_dispersive_ranges = 1
dispersive_range_size_um = bandwidth_um / num_dispersive_ranges

dispersive_ranges_um = [
	[ lambda_min_um + idx * dispersive_range_size_um, lambda_min_um + ( idx + 1 ) * dispersive_range_size_um ] for idx in range( 0, num_dispersive_ranges )
]


#
# Fabrication Constraints
#
min_feature_size_um = 0.75
min_feature_size_voxels = min_feature_size_um / mesh_spacing_um
blur_half_width_voxels = int( np.ceil( (min_feature_size_voxels - 1) / 2. ) )

num_vertical_layers = 10

#
# FDTD
#
vertical_gap_size_um = 15
lateral_gap_size_um = 4

fdtd_region_size_vertical_um = 2 * vertical_gap_size_um + device_size_verical_um + focal_length_um
fdtd_region_size_lateral_um = 2 * lateral_gap_size_um + device_size_lateral_um
fdtd_region_maximum_vertical_um = device_size_verical_um + vertical_gap_size_um
fdtd_region_minimum_vertical_um = -focal_length_um - vertical_gap_size_um

fdtd_region_minimum_vertical_voxels = int( np.ceil(fdtd_region_size_vertical_um / mesh_spacing_um) )
fdtd_region_minimum_lateral_voxels = int( np.ceil(fdtd_region_size_lateral_um / mesh_spacing_um) )

fdtd_simulation_time_fs = 9000

#
# Forward Source
#
lateral_aperture_um = 1.1 * device_size_lateral_um
src_maximum_vertical_um = device_size_verical_um + vertical_gap_size_um * 2. / 3.
src_minimum_vertical_um = -focal_length_um - 0.5 * vertical_gap_size_um
# src_beam_rad = device_size_lateral_um/2
src_angle_incidence = 0 # degrees
src_phi_incidence = 0 # degrees

src_hgt_Si = 1.3
src_hgt_polymer = src_maximum_vertical_um + 1.5 - device_size_verical_um
n_Si = 3.43
n_polymer = max_device_index

assert ( src_maximum_vertical_um + 1 ) < ( fdtd_region_maximum_vertical_um - silicon_thickness_um ), "The source is in the silicon"


#
# Spectral and polarization selectivity information
#
# polarizations_focal_plane_map = [ ['x', 'y'], ['x', 'y'], ['x', 'y'], ['x', 'y'] ]
polarizations_focal_plane_map = [ ['x', 'y'], ['x'], ['x', 'y'], ['y'] ]
weight_focal_plane_map = [ 0.5, 1.0, 0.5, 1.0 ]
polarization_name_to_idx = { 'x':0, 'y':1, 'z':2 }
# We are assuming that the data is organized in order of increasing wavelength (i.e. - blue first, red last)
spectral_focal_plane_map = [
	[0, num_points_per_band],
	[num_points_per_band, 2 * num_points_per_band],
	[2 * num_points_per_band, 3 * num_points_per_band],
	[num_points_per_band, 2 * num_points_per_band]
]

#
# Adjoint sources
#
adjoint_vertical_um = -focal_length_um
num_focal_spots = 4
num_adjoint_sources = num_focal_spots
adjoint_x_positions_um = [device_size_lateral_um / 4., -device_size_lateral_um / 4., -device_size_lateral_um / 4., device_size_lateral_um / 4.]
adjoint_y_positions_um = [device_size_lateral_um / 4., device_size_lateral_um / 4., -device_size_lateral_um / 4., -device_size_lateral_um / 4.]


# dispersive_range_to_adjoint_src_map = [
	# [ 0 ],
	# [ 1, 3 ],
	# [ 2 ]
# ]

# adjoint_src_to_dispersive_range_map = [
	# 0,
	# 1,
	# 2,
	# 1
# ]

dispersive_range_to_adjoint_src_map = [
	[ 0, 1, 2, 3 ]
]

adjoint_src_to_dispersive_range_map = [
	0,
	0,
	0,
	0
]


#
# Optimization
#
num_epochs = 10
num_iterations_per_epoch = 30

use_fixed_step_size = False
fixed_step_size = 2.0

epoch_start_permittivity_change_max = 0.15
epoch_end_permittivity_change_max = 0.05
epoch_range_permittivity_change_max = epoch_start_permittivity_change_max - epoch_end_permittivity_change_max

epoch_start_permittivity_change_min = 0.05
epoch_end_permittivity_change_min = 0
epoch_range_permittivity_change_min = epoch_start_permittivity_change_min - epoch_end_permittivity_change_min

#
# Topology
#
topology_num_free_iterations_between_patches = 8

