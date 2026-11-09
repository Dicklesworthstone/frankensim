//! Source-owned rheology and fully implicit convection with physical RCR loads.

use super::{
    MaterialFluxError, MaterialFluxSystem, MaterialStokesWindkesselHistory,
    MaterialStokesWindkesselPort, MaterialStokesWindkesselSolution, physical_windkessel,
    with_kinematic_windkessel,
};
use crate::material::{
    MaterialFluxRheologyPoint, MaterialNavierStokesBoundary, MaterialNavierStokesRheologyStep,
    MaterialNavierStokesStep, rheology::RheologyField, with_kinematic_boundaries,
};
use crate::ns::boundary::NavierStokesControls;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Couple source-owned laminar convection and RCR/inertance at one BDF endpoint.
    ///
    /// The original source's frozen dynamic viscosity and density supply fluid
    /// stress, mass and convection. External storage uses the same physical
    /// reference normalization and paired fluid/compliance histories. None starts
    /// BDF1; older history and its interval select variable-step BDF2. Original
    /// natural/dissipative backflow and incoming reservoir momentum are retained.
    /// No external resistor or numerical BDF power is added to bulk fluid heat.
    ///
    /// # Errors
    /// Retains both dynamic laws and original nonlinear endpoint checks; no
    /// source, history, temperature or partial storage update is published.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_variable_step_windkessel_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesWindkesselSolution<MaterialNavierStokesStep>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let native = with_kinematic_windkessel(
            density,
            ports,
            loads,
            previous,
            older_history,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |loads, previous, older| {
                with_kinematic_boundaries(density, boundaries, |boundaries| {
                    self.system
                        .bdf2_navier_stokes_mixed_variable_step_windkessel_ports(
                            self.params,
                            &|point| force_density(point).map(|f| f / density),
                            boundaries,
                            ports,
                            loads,
                            previous,
                            older,
                            timestep_s,
                            controls,
                            absolute_pressure_tolerance_pa / density,
                            absolute_flow_tolerance_m2_s,
                        )
                })
                .map_err(MaterialFluxError::Boundary)
            },
        )?;
        physical_windkessel(native, density, |endpoint, _| {
            self.material_navier_stokes_step(endpoint)
        })
    }

    /// Resolve actual source rheology, laminar convection and hydraulic storage
    /// together, including Carreau-Yasuda sources at their accepted shear rate.
    ///
    /// The existing material resolver/cache owns every constitutive query and
    /// retains the exact accepted source per cell. Each nonlinear iteration
    /// rebuilds convection and viscosity before imposing the external port law;
    /// the final original equations, finite-resistance law and both unreduced
    /// dynamic laws must all pass. Pressure and the complete field relax together.
    /// Native velocity histories initialize both convection and canonical port
    /// flow histories; they are not replaced by the BDF mass-history blend.
    ///
    /// Density, thermodynamic conditions, composition and external C,L,R remain
    /// fixed through supplied histories. The reservoir reference/units and
    /// separate storage, physical losses and numerical BDF defect are identical
    /// to the linear physical Windkessel entry. Instantaneous bulk heat is not
    /// integrated temperature, and convective kinetic-energy transport remains
    /// separate from the external circuit account. No turbulence model is added.
    ///
    /// # Errors
    /// Preserves the original source error before solver refusals. Neither
    /// input history nor accepted device state is modified on failure; no
    /// candidate source set, flow or compliance state is returned prematurely.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<
        MaterialStokesWindkesselSolution<MaterialNavierStokesRheologyStep<G>>,
        MaterialFluxError,
    > {
        let density = self.transport.density_kg_m3();
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = with_kinematic_windkessel(
            density,
            ports,
            loads,
            previous,
            older_history,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |loads, previous, older| {
                with_kinematic_boundaries(density, boundaries, |boundaries| {
                    self.system
                        .bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports_checked(
                            self.params,
                            &|cell, _, rate| field.viscosity(cell, rate),
                            &|point| force_density(point).map(|f| f / density),
                            boundaries,
                            ports,
                            loads,
                            previous,
                            older,
                            timestep_s,
                            controls,
                            absolute_pressure_tolerance_pa / density,
                            absolute_flow_tolerance_m2_s,
                            &|| field.is_valid(),
                        )
                })
                .map_err(MaterialFluxError::Boundary)
            },
        );
        let (native, samples) = field.finish(result)?;
        physical_windkessel(native, density, |endpoint, _| {
            self.material_navier_stokes_rheology_step(endpoint, samples)
        })
    }
}
