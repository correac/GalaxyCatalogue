import numpy as np
import unyt
import velociraptor
from swiftsimio import load
from typing import List, Union, Tuple, Dict


class HaloCatalogue:
    """
    General class containing halo properties
    """

    def __init__(
        self, path_to_catalogue: str, dm_particle_mass: float, simulation_type: str,
    ):
        """
        Parameters
        ----------
        path_to_catalogue: str
        Path to the catalogue with halo properties
        dm_particle_mass: unyt.array.unyt_quantity
        Minimum dark matter particle mass in units of Msun. Haloes that contain less than
        1000 dark mattter particles are disregarded
        """

        self.path_to_catalogue = path_to_catalogue

        # Load catalogue using velociraptor python library
        catalogue = velociraptor.load(self.path_to_catalogue)

        # Selecting haloes that contain at least 100 DM particles
        mask = np.where(
            catalogue.masses.mass_200crit.to("Msun").value >= unyt.unyt_quantity(100 * dm_particle_mass, "Msun")
        )[0]

        # Compute the number of haloes following the selection mask
        self.number_of_haloes = len(mask)

        # Structure type
        self.structure_type = catalogue.structure_type.structuretype[mask]

        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        self.log10_mvir = np.log10(
            catalogue.masses.mass_tot.to("Msun").value[mask]
        )

        self.log10_stellar_mass = np.log10(
            np.clip(catalogue.apertures.mass_star_50_kpc.to("Msun").value[mask], 1, 1e16)
        )

        self.concentration = catalogue.concentration.cnfw.value[mask]
        self.virial_radius = catalogue.radii.r_200crit.to("Mpc").value[mask]

        self.scale_radius = self.virial_radius / self.concentration

        # Ids of haloes satisfying the selection criterion
        self.halo_index = mask.copy()

        self.xminpot = catalogue.positions.xcminpot.to("Mpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("Mpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("Mpc").value[mask]

class SimInfo:

    def __init__(
        self,
        directory: str,
        snapshot: str,
        catalogue: str,
        name: Union[str, None],
        output: str,
        simtype: str,
    ):
        """
        Parameters
        ----------
        directory: str
        Run directory
        snapshot: str
        Name of the snapshot file
        catalogue: str
        Name of the catalogue file
        name:
        Name of the run
        galaxy_min_stellar_mass: unyt.array.unyt_quantity
        """

        self.directory = directory
        self.snapshot_name = snapshot
        self.catalogue_name = catalogue
        self.output_path = output
        self.simulation_type = simtype
        self.path_to_catalogue = f"{self.directory}/{self.catalogue_name}"

        # Load snapshot via swiftsimio
        self.snapshot = load(f"{self.directory}/{self.snapshot_name}")

        # Box size of the simulation in kpc
        self.boxSize = self.snapshot.metadata.boxsize.to("Mpc")

        # Cosmic scale factor
        self.a = self.snapshot.metadata.scale_factor

        self.hubble_time_Gyr = self.snapshot.metadata.time.to("Gyr").value

        self.Omega_m = self.snapshot.metadata.cosmology.Om0

        self.Omega_l = 1 - self.Omega_m

        self.h = self.snapshot.metadata.cosmology.H0.value / 100.

        self.z = self.snapshot.metadata.redshift

        self.rhocrit0 = 2.7754e11 * self.h ** 2  # Msun / Mpc^3

        # Conversion from internal units to kpc
        kpc = 3.08567758e21
        self.to_kpc_units = (
            self.snapshot.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            / kpc
        )

        # Conversion from internal units to Msun
        Msun = 1.9891e33
        self.to_Msun_units = (
            self.snapshot.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            / Msun
        )

        # Maximum softening
        self.softening = (
            self.snapshot.metadata.gravity_scheme[
                "Maximal physical DM softening length (Plummer equivalent) [internal units]"
            ][0] * self.to_kpc_units
        )

        # Maximum softening for baryons
        self.baryon_max_soft = (
            self.snapshot.metadata.gravity_scheme[
                "Maximal physical baryon softening length  [internal units]"
            ][0]
            * self.to_kpc_units
        )

        self.dm_particle_mass = 1e6  # Need to fix, table mass is empty

        # Object containing halo properties (from halo catalogue)
        self.halo_data = HaloCatalogue(
            path_to_catalogue=f"{self.directory}/{self.catalogue_name}",
            dm_particle_mass=self.dm_particle_mass,
            simulation_type=self.simulation_type
        )