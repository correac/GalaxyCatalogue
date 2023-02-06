import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np
import unyt
from unyt import unyt_array
from unyt import Mpc, kpc, pc

from swiftsimio import load, mask
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)

from swiftascmaps import evermore, evermore_r, red_tv, red_tv_r, reputation, midnights
from astropy.visualization import make_lupton_rgb

from argumentparser import ArgumentParser
from simulation_data import SimInfo


def calculate_angular_momentum(part_pos, part_vel, part_mass, origin):

    pos = part_pos - origin
    # Compute distances
    distancesDATA = np.sqrt(np.sum(pos ** 2, axis=1))

    # Restrict particles to 30kpc region
    extract = distancesDATA < 30.0
    pos = pos[extract,:]
    vel = part_vel[extract,:]
    mass = part_mass[extract]

    distancesDATA = distancesDATA[extract]

    Mstar = np.sum(mass)
    dvVmass = np.sum(mass[:, np.newaxis] * vel, axis=0) / Mstar
    vel -= dvVmass # correct for bulk velocity

    # Compute momentum within 5 kpc region
    extract = distancesDATA < 5.0
    smomentum_inner_5kpc = np.cross(pos[extract,:], vel[extract, :])
    momentum_inner_5kpc = np.sum(mass[extract][:, np.newaxis] * smomentum_inner_5kpc, axis=0)

    return momentum_inner_5kpc

def make_rotation_matrix(ang_momentum):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)
    edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum, axis="y")

    return face_on_rotation_matrix, edge_on_rotation_matrix


def translate_region(origin, size, boxSize, region):

    for i in range(3):
        if (origin[i] + size[0] * 0.5 > boxSize[i]) or (origin[i] - size[0] * 0.5 < 0):
            region[i] = unyt.unyt_array([0, boxSize[i].value], 'Mpc')

    return region

def translate_coordinates(coordinates, origin, size, boxSize, region_defined):

    new_region = region_defined.copy()
    new_origin = origin.copy()

    for i in range(3):

        if (origin[i] + size[0] * 0.5 > boxSize[i]):
            select_region = np.where(coordinates[:,i] <= origin[i] + size[0] * 0.5 - boxSize[i])[0]
            coordinates[select_region,i] += boxSize[i]
            coordinates[:, i] -= (origin[i] + size[0] * 0.5 - boxSize[i])
            new_origin[i] -= (origin[i] + size[0] * 0.5 - boxSize[i])

            if i == 0:
                new_region[0:2] -= (origin[i] + size[0] * 0.5 - boxSize[i])
            if i == 1:
                new_region[2:4] -= (origin[i] + size[0] * 0.5 - boxSize[i])

        if (origin[i] - size[0] * 0.5 < 0):
            select_region = np.where(coordinates[:,i] >= origin[i] - size[0] * 0.5 + boxSize[i])[0]
            coordinates[select_region,i] -= boxSize[i]
            coordinates[:, i] -= (origin[i] - size[0] * 0.5)
            new_origin[i] -= (origin[i] - size[0] * 0.5)

            if i == 0:
                new_region[0:2] -= (origin[i] - size[0] * 0.5)
            if i == 1:
                new_region[2:4] -= (origin[i] - size[0] * 0.5)

    return coordinates, new_region, new_origin


def get_projection_map(sim_info, npix, size, origin, halo_index, rotation):

    rotation_vector = np.array([1, 0, 0])  # ad-hoc

    size = unyt.unyt_array([size, size, size], 'Mpc')
    region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]
    region_defined = unyt.unyt_array([region[0][0], region[0][1], region[1][0], region[1][1]], 'Mpc')
    region = translate_region(origin, size, sim_info.boxSize, region)

    data_mask = mask(f"{sim_info.directory}/{sim_info.snapshot_name}")
    data_mask.constrain_spatial(region)
    data = load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=data_mask)

    data.stars.coordinates = data.stars.coordinates.to_physical()

    data.stars.coordinates, region_defined, new_origin = \
        translate_coordinates(data.stars.coordinates, origin, size, sim_info.boxSize, region_defined)

    data.stars.smoothing_lengths = generate_smoothing_lengths(
        coordinates=data.stars.coordinates,
        boxsize=data.metadata.boxsize,
        kernel_gamma=2.5,
        neighbours=20,
        speedup_fac=2,
        dimension=3,
    )

    # Set here your vector if you want some rotation
    # if rotation flag = 1 or -1, we make rotation based on stellar angular momentum
    if np.abs(rotation) == 1:
        rotation_vector = calculate_angular_momentum(data.stars.coordinates, data.stars.velocities,
                                                     data.stars.masses, new_origin)

    # if rotation flag = 2 we make random rotation
    if rotation == 2:
        u = np.random.uniform(0,1,1)[0]
        theta = np.arccos(1. - 2 * u)
        phi = 2 * np.pi * np.random.uniform(0,1,1)[0]
        rotation_vector = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi),np.cos(theta)])

    # if rotation flag > 2 we rotate by given angle
    if rotation > 2:
        theta = rotation
        phi = 2 * np.pi * np.random.uniform(0,1,1)[0]
        rotation_vector = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    face_on_rotation_matrix, edge_on_rotation_matrix = make_rotation_matrix(rotation_vector)

    if rotation == -1:
        rotation_matrix = edge_on_rotation_matrix
    else:
        rotation_matrix = face_on_rotation_matrix

    # Calculate particles luminosity in i, r and g bands
    luminosities = [
        data.stars.luminosities.GAMA_i,
        data.stars.luminosities.GAMA_r,
        data.stars.luminosities.GAMA_g,
    ]
    rgb_image_face = np.zeros((npix, npix, len(luminosities)))

    for ilum in range(len(luminosities)):

        data.stars.usermass = luminosities[ilum]

        # Face on projection
        if rotation == 0:
            pixel_grid = project_pixel_grid(
                data=data.stars,
                resolution=npix,
                project="usermass",
                parallel=True,
                region=region_defined,
                boxsize=data.metadata.boxsize,
                backend="subsampled",
            )
        else:
            pixel_grid = project_pixel_grid(
                data=data.stars,
                resolution=npix,
                project="usermass",
                parallel=True,
                region=region_defined,
                rotation_center=new_origin,
                rotation_matrix=rotation_matrix,
                boxsize=data.metadata.boxsize,
                backend="subsampled",
            )

        x_range = region[0][1] - region[0][0]
        y_range = region[1][1] - region[1][0]
        units = 1.0 / (x_range * y_range)
        units.convert_to_units(1.0 / (x_range.units * y_range.units))

        mass_map_face = unyt_array(pixel_grid, units=units)
        if size[0] > 0.5:
            mass_map_face.convert_to_units(1.0 / kpc ** 2)
        else:
            mass_map_face.convert_to_units(1.0 / pc ** 2)
        try:
            mass_map_face[mass_map_face == 0.0] = mass_map_face[
                mass_map_face > 0.0
            ].min()
        except:
            mass_map_face[mass_map_face == 0.0] = 1.0e-10

        rgb_image_face[:, :, ilum] = mass_map_face.T

    # # Some plotting options based on aperture
    if size[0] > 1:
        Q = 30
        stretch = 50
    elif (size[0] <= 1) & (size[0] > 0.5):
        Q=50
        stretch=30
    elif (size[0] <= 0.5) & (size[0]>0.1):
        Q = 20
        stretch = 0.01
    else:
        Q=10
        stretch=0.5

    stars_map = make_lupton_rgb(
        rgb_image_face[:, :, 0],
        rgb_image_face[:, :, 1],
        rgb_image_face[:, :, 2],
        Q=Q,
        stretch=stretch,
    )

    data.dark_matter.coordinates = data.dark_matter.coordinates.to_physical()

    data.dark_matter.coordinates, _, _ = \
        translate_coordinates(data.dark_matter.coordinates, origin, size, sim_info.boxSize, region_defined)

    data.dark_matter.smoothing_lengths = generate_smoothing_lengths(
        data.dark_matter.coordinates,
        data.metadata.boxsize,
        kernel_gamma=2.5,
        neighbours=57,
        speedup_fac=2,
        dimension=3,
    )

    # Project the dark matter mass
    if rotation == 0:
        dm_grid = project_pixel_grid(
            data=data.dark_matter,
            resolution=npix,
            project="masses",
            parallel=True,
            region=region_defined,
            boxsize=data.metadata.boxsize,
            backend="subsampled",
        )
    else:
        dm_grid = project_pixel_grid(
            data=data.dark_matter,
            resolution=npix,
            project="masses",
            parallel=True,
            region=region_defined,
            rotation_center=new_origin,
            rotation_matrix=rotation_matrix,
            boxsize=data.metadata.boxsize,
            backend="subsampled",
        )

    dm_map = unyt_array(dm_grid, units=units)
    dm_map.convert_to_units(1.0 / kpc ** 2)
    dm_map = dm_map.T

    return stars_map, dm_map, region_defined


def make_image(sim_info, npix, size, origin, halo_index, rotation):

    stars_map, dm_map, region = get_projection_map(
        sim_info, npix, size, origin, halo_index, rotation
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
    filename = sim_info.output_path+"/stars_light_map_galaxy_%i"%halo_index
    filename += "_size_%.2fMpc"%size
    fig.savefig(filename+".png", dpi=500)
    plt.close()

    # Nienke select here a colormap
    colormap = evermore_r
    # colormap = midnights
    # colormap = 'binary_r'
    # colormap = 'viridis'
    # colormap = 'magma'


    fig = plt.figure(figsize=(6.0, 6.0))
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    ax = plt.subplot(1,1,1)
    ax.tick_params(labelleft=False, labelbottom=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    ax.imshow(LogNorm()(dm_map), extent=region, interpolation='nearest', cmap=colormap)
    filename = sim_info.output_path+"/dm_mass_map_galaxy_%i"%halo_index
    filename += "_size_%.2fMpc"%size
    fig.savefig(filename+".png", dpi=500)
    plt.close()

def make_output_data(sim_info, sample, size, num_galaxies):

    # Do we need any data for Ben's parametric model?
    # Halo mass? Halo radius? Galaxy mass? Size? Morphology? Angular momentum?
    M200c = sim_info.halo_data.log10_halo_mass[sample]
    R200c = sim_info.halo_data.virial_radius[sample]
    Mstellar = sim_info.halo_data.log10_stellar_mass[sample]
    num_structure = np.zeros(len(sample))

    num_haloes_all = len(sim_info.halo_data.log10_halo_mass)
    haloes_position = np.zeros((num_haloes_all, 3))
    haloes_position[:,0] = sim_info.halo_data.xminpot
    haloes_position[:,1] = sim_info.halo_data.yminpot
    haloes_position[:,2] = sim_info.halo_data.zminpot

    for index in range(num_galaxies):

        x = sim_info.halo_data.xminpot[sample[index]]
        y = sim_info.halo_data.yminpot[sample[index]]
        z = sim_info.halo_data.zminpot[sample[index]]

        origin = np.array([x, y, z])

        volume = unyt.unyt_array([size, size, size])
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(volume, origin)]
        region_defined = unyt.unyt_array([region[0][0], region[0][1], region[1][0], region[1][1]])

        haloes_new_position, _, new_origin = \
            translate_coordinates(haloes_position, origin, volume, sim_info.boxSize, region_defined)

        pos = haloes_new_position - new_origin
        distance = np.sqrt(np.sum(pos ** 2, axis=1))
        num_structure[index] = len(np.where(distance <= size)[0])


    # Output data
    filename = sim_info.output_path + "/galaxy_data_size_%.2fMpc.txt"%size

    with open(filename, 'w') as f:
        # Some header
        line = '# Galaxy No. - Stellar Mass - Virial Mass - Virial Radius - Num. Substructure'
        f.write(line)
        f.write('\n')

        for index in range(num_galaxies):
            line = "%i"%index+"  %.2f"%Mstellar[index]+"  %.2f"%M200c[index]+"  %.2f"%R200c[index]+"  %i"%num_structure[index]
            f.write(line)
            f.write('\n')


def plot_sample(sim_info, num_galaxies, npix, size, rotation):

    # We are imposing to plot central galaxies only
    sample = np.where(sim_info.halo_data.structure_type == 10)[0]

    # We output some information regarding the image we generate
    make_output_data(sim_info, sample, size, num_galaxies)

    for index in range(num_galaxies):

        x = sim_info.halo_data.xminpot[sample[index]]
        y = sim_info.halo_data.yminpot[sample[index]]
        z = sim_info.halo_data.zminpot[sample[index]]

        origin = unyt.unyt_array([x, y, z], 'Mpc') / sim_info.a  # to comoving

        make_image(sim_info, npix, size, origin, index, rotation)


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
    size = 1.0          # [Mpc] size of image
    num_galaxies = 1    # number of galaxies to plot
    rotation = 50       # Rotation flag. 0: no rotation,
                        #                1: rotation to make stellar disc face-on,
                        #               -1: rotation to make stellar disc edge-on,
                        #                2: random rotation,
                        #               >2: rotation in given angle

    # This script will plot the most massive galaxies
    # in the box. If num_galaxies = 10, it will correspond
    # to making visualization of the 10 most massive galaxies
    # and their respective haloes.
    plot_sample(sim_info, num_galaxies, npix, size, rotation)

