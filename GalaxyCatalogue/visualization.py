import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import glob
import numpy as np
import unyt
from unyt import unyt_array
from unyt import Mpc, kpc
from typing import List, Union, Tuple, Dict

from swiftsimio import load, mask
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)
from swiftascmaps import evermore, evermore_r, red_tv, red_tv_r
from astropy.visualization import make_lupton_rgb

import velociraptor
from argumentparser import ArgumentParser
from simulation_data import SimInfo

def make_rotation_matrix(ang_momentum):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)
    edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum, axis="y")

    return face_on_rotation_matrix, edge_on_rotation_matrix


def get_projection_map(sim_info, npix, size, origin):

    size = unyt.unyt_array([size, size, size], 'Mpc')
    region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]
    region_defined = unyt.unyt_array([region[0][0], region[0][1], region[1][0], region[1][1]], 'Mpc')

    data_mask = mask(f"{sim_info.directory}/{sim_info.snapshot_name}")
    data_mask.constrain_spatial(region)
    data = load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=data_mask)

    data.stars.coordinates = data.stars.coordinates.to_physical()

    data.stars.smoothing_lengths = generate_smoothing_lengths(
        coordinates=data.stars.coordinates,
        boxsize=data.metadata.boxsize,
        kernel_gamma=2.5,
        neighbours=20,
        speedup_fac=2,
        dimension=3,
    )

    # data.stars.coordinates = data.stars.coordinates.to_physical()
    # data.stars.smoothing_lengths = data.stars.smoothing_lengths.to_physical()

    # Set here your vector if you want some rotation
    # rotation_vector = np.array([1, 0, 0]) #ad-hoc
    # face_on_rotation_matrix, edge_on_rotation_matrix = make_rotation_matrix(
    #     rotation_vector
    # )

    luminosities = [
        data.stars.luminosities.GAMA_i,
        data.stars.luminosities.GAMA_r,
        data.stars.luminosities.GAMA_g,
    ]
    rgb_image_face = np.zeros((npix, npix, len(luminosities)))

    for ilum in range(len(luminosities)):

        # Face on projection
        data.stars.usermass = luminosities[ilum]
        pixel_grid = project_pixel_grid(
            data=data.stars,
            resolution=npix,
            project="usermass",
            parallel=True,
            region=region_defined,
            # region=visualise_region,
            # rotation_center=unyt.unyt_array([x, y, z]),
            # rotation_matrix=face_on_rotation_matrix,
            boxsize=data.metadata.boxsize,
            backend="subsampled",
        )

        x_range = region[0][1] - region[0][0]
        y_range = region[1][1] - region[1][0]
        units = 1.0 / (x_range * y_range)
        units.convert_to_units(1.0 / (x_range.units * y_range.units))

        mass_map_face = unyt_array(pixel_grid, units=units)
        mass_map_face.convert_to_units(1.0 / kpc ** 2)
        try:
            mass_map_face[mass_map_face == 0.0] = mass_map_face[
                mass_map_face > 0.0
            ].min()
        except:
            mass_map_face[mass_map_face == 0.0] = 1.0e-10

        rgb_image_face[:, :, ilum] = mass_map_face.T

    stars_map = make_lupton_rgb(
        rgb_image_face[:, :, 0],
        rgb_image_face[:, :, 1],
        rgb_image_face[:, :, 2],
        Q=30,
        stretch=50,
    )

    data.dark_matter.smoothing_lengths = generate_smoothing_lengths(
        data.dark_matter.coordinates,
        data.metadata.boxsize,
        kernel_gamma=2.5,
        neighbours=57,
        speedup_fac=2,
        dimension=3,
    )

    # Project the dark matter mass
    dm_grid = project_pixel_grid(
        data=data.dark_matter,
        resolution=npix,
        project="masses",
        parallel=True,
        region=region_defined,
        # region=visualise_region,
        # rotation_center=unyt.unyt_array([x, y, z]),
        # rotation_matrix=face_on_rotation_matrix,
        boxsize=data.metadata.boxsize,
        backend="subsampled",
    )

    dm_map = unyt_array(dm_grid, units=units)
    dm_map.convert_to_units(1.0 / kpc ** 2)
    dm_map = dm_map.T

    return stars_map, dm_map, region_defined


def make_image(sim_info, npix, size, origin, halo_index):

    stars_map, dm_map, region = get_projection_map(
        sim_info, npix, size, origin
    )

    ######################
    fig = plt.figure(figsize=(6.0, 6.0))
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    ax = plt.subplot(1,1,1)
    ax.tick_params(labelleft=False, labelbottom=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    ax.imshow(stars_map, extent=region, interpolation='nearest')
    fig.savefig(sim_info.output_path+"/stars_light_map_halo_%i.png"%halo_index, dpi=500)
    plt.close()

    fig = plt.figure(figsize=(6.0, 6.0))
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    ax = plt.subplot(1,1,1)
    ax.tick_params(labelleft=False, labelbottom=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    ax.imshow(LogNorm()(dm_map), extent=region, interpolation='nearest', cmap=evermore_r)
    fig.savefig(sim_info.output_path+"/dm_mass_map_halo_%i.png"%halo_index, dpi=500)
    plt.close()


def plot_sample(sim_info, num_galaxies, npix, size):

    # We are impossing to plot central galaxies only
    sample = np.where(sim_info.halo_data.structure_type == 10)[0]

    # Do we need any data for Ben's parametric model?
    # Halo mass? Halo radius? Galaxy mass? Size? Morphology? Angular momentum?
    M200c = sim_info.halo_data.log10_halo_mass[sample]
    R200c = sim_info.halo_data.virial_radius[sample]

    for index in range(num_galaxies):

        x = sim_info.halo_data.xminpot[sample[index]]
        y = sim_info.halo_data.yminpot[sample[index]]
        z = sim_info.halo_data.zminpot[sample[index]]

        origin = unyt.unyt_array([x, y, z], 'Mpc') / sim_info.a  # to comoving

        make_image(sim_info, npix, size, origin, index)


if __name__ == "__main__":

    config_parameters = ArgumentParser()

    # Load all data and save it in SimInfo class
    sim_info = SimInfo(
        directory=config_parameters.directory_list[0],
        snapshot=config_parameters.snapshot_list[0],
        catalogue=config_parameters.catalogue_list[0],
        name=config_parameters.name_list[0],
        output=config_parameters.output_directory,
        simtype=config_parameters.sim_type[0],
    )

    # Some options:
    npix = 1024         # number of pixels
    size = 3            # [Mpc] size of image
    num_galaxies = 10   # number of galaxies to plot

    # This script will plot the most massive galaxies
    # in the box. If num_galaxies = 10, it will correspond
    # to making visualization of the 10 most massive galaxies
    # and their respective haloes
    plot_sample(sim_info, num_galaxies, npix, size)

