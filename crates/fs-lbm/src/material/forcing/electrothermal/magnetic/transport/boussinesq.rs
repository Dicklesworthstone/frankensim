//! Density, resistivity and transport from one local material evaluation.
//!
//! This joins the existing Boussinesq projection to the existing electrical
//! source/heat/flow step. The raw material stays attached to each accepted
//! sample; reference density and heat capacity continue to own inertia and
//! sensible storage. No second buoyancy or electromagnetic law is introduced.

use super::*;
use crate::material::transport::boussinesq::{
    BackgroundBuoyancy, BoussinesqConfig, BoussinesqMaterial, candidate_with_buoyancy,
};

#[allow(clippy::too_many_arguments)]
fn advance<S: MovingSlab<D, Q> + BackgroundBuoyancy, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; D]],
    depth_m: f64, isolated: Option<f64>, config: BoussinesqConfig<D>, gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
    poll(gate)?;
    config.validate()?;
    super::advance::<S, BoussinesqMaterial<F>, D, Q>(
        slab, domain, plane, input, field, forces, isolated, gate,
        |cell, temperature, reference| {
            let response = resolve(cell, temperature, reference)?;
            // The original source is evaluated once. Its local density drives
            // buoyancy, its mu/k drive diffusion, and its paired resistivity
            // drives current at this SAME temperature and source context.
            ElectricalTransport::try_new(
                BoussinesqMaterial::try_new(response.material, reference, config)?,
                response.resistivity_ohm_m,
            )
        },
        |mapping, power, total_forces, samples| {
            let heat = mapping.heat_inputs(
                power, total_forces, mapping.lattice.layout()?.0, depth_m, input.limits.thermal,
            )?;
            candidate_with_buoyancy::<S, ElectricalTransport<BoussinesqMaterial<F>>, D, Q>(
                mapping, heat, total_forces, samples, config, gate,
                |sample| sample.material().resolved_transport().density_kg_m3(),
            )
        },
    )
}

impl MaterialSlab<ThermalLbm> {
    /// Advance current, Joule heat, magnetic force and material-driven buoyancy.
    ///
    /// A single local callback supplies electrical resistivity and the actual
    /// fluid source at the accepted temperature. The source's density supplies
    /// `(rho-rho0)*gravity`; its viscosity/conductivity supply local transport.
    /// The existing Boussinesq projection retains rho0/Cp0 for inertia/storage
    /// and checks the declared departures at BOTH initial and endpoint states.
    /// The slab's independent beta load is replaced for this step, not added.
    ///
    /// Field is tesla, external force is N/m³, depth is metres. As in the
    /// original 2-D magnetic step, the field must be out of plane. Electrical
    /// insulation remains an explicit `None`, not a missing-property fallback.
    /// Source receipts remain in `sample.material.material().material()`.
    /// Failure or cancellation preserves populations, loads and reference units.
    ///
    /// Rates/coefficients are explicit initial-state values for one step;
    /// endpoint evaluation checks source support, not midpoint convergence.
    /// Buoyancy work is NOT Joule heat and is not included in the reported
    /// electrical work identity. No gravitational/kinetic total-energy balance,
    /// induction, compressible work, species evolution or phase change is claimed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 2]], depth_m: f64, config: BoussinesqConfig<2>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        self.step_buoyant_material_motion(input, field_t, forces_n_m3, depth_m, None, config, gate, resolve)
    }

    /// Electrode-free counterpart: the same material buoyancy/transport with
    /// induced current, Joule heat and Lorentz force. `current_tolerance_a` is
    /// the original isolated current owner's local charge-balance tolerance.
    /// A nonempty electrode list is refused before material queries.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64,
        config: BoussinesqConfig<2>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        self.step_buoyant_material_motion(
            input, field_t, forces_n_m3, depth_m, Some(current_tolerance_a), config, gate, resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_buoyant_material_motion<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]],
        forces: &[[f64; 2]], depth_m: f64, isolated: Option<f64>,
        config: BoussinesqConfig<2>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) {
            return Err(invalid("invalid represented 2D cell volume").into());
        }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, F, 2, 5>(
            self, domain, plane, input, field, forces, depth_m, isolated, config, gate, resolve,
        )
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional material buoyancy/current/heat/flow coupling using
    /// actual dx³ volumes, all vector components and D3Q7 heat conversion.
    /// Shares sampling, accounting and atomic publication with the 2-D path.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 3]], config: BoussinesqConfig<3>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        self.step_buoyant_material_motion(input, field_t, forces_n_m3, None, config, gate, resolve)
    }

    /// Three-dimensional electrode-free counterpart with original isolated
    /// charge/current semantics and no fictitious electrical supply work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 3]], current_tolerance_a: f64,
        config: BoussinesqConfig<3>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        self.step_buoyant_material_motion(
            input, field_t, forces_n_m3, Some(current_tolerance_a), config, gate, resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_buoyant_material_motion<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]],
        forces: &[[f64; 3]], isolated: Option<f64>, config: BoussinesqConfig<3>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<BoussinesqMaterial<F>>, MaterialMagneticError> {
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, F, 3, 7>(
            self, domain, plane, input, field, forces, dx, isolated, config, gate, resolve,
        )
    }
}
