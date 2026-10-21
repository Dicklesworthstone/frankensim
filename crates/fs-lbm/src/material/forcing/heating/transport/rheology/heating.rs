//! Generalized-Newtonian stress power on the existing thermal/flow step.
//!
//! For the admitted incompressible law `stress = 2 mu D`, its irreversible
//! power is `stress:D = 2 mu D:D = mu gamma_dot^2`. The same initial material
//! evaluation supplies both collision viscosity and this nonnegative source.
//! No numerical kinetic-energy residual, pressure work, bulk-viscosity loss,
//! turbulence closure or elastic energy is relabelled as thermal power.

use super::super::super::Inputs;
use super::super::advance_with_sampler_and_sources;
use super::{
    FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError, MaterialSlab,
    MaterialTransportError, MaterialTransportStep, RheologyMaterial, RheologyPoint, RheologySolve,
    StressSlab, ThermalLbm, ThermalLbm3, TransportSample, poll, sample,
};
use fs_exec::CancelGate;

/// Explicit source partition for an accepted physical viscous-heating step.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousHeating {
    /// `mu * gamma_dot^2` [W/m3], one value per cell including zero wall entries.
    /// These are the applied initial-state values, held over one time step.
    pub power_w_m3: Vec<f64>,
    /// Integrated constitutive stress work converted into sensible heat [J].
    /// This is not an independently measured decrement of lattice kinetic energy.
    pub viscous_heat_j: f64,
    /// Separately authored signed volumetric heat integrated over the step [J].
    /// It must exclude the molecular viscous heating owned by this method.
    pub authored_heat_j: f64,
    /// Combined thermal source minus authored heat minus viscous heat [J].
    /// Checked against the same caller-specified thermal balance tolerance.
    pub source_partition_defect_j: f64,
}

/// Original material/thermal receipts plus the explicitly owned viscous heat.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousRheologyStep<F> {
    /// The same transport step and both source epochs as the unheated API.
    /// Its heat source is the sum of authored and physical viscous heat.
    pub transport: MaterialTransportStep<RheologyMaterial<F>>,
    /// Separate physical loss and external-source accounting.
    pub viscous: ViscousHeating,
}

#[allow(clippy::too_many_arguments)]
fn source_inputs<S, F, const D: usize>(
    mapping: &MaterialSlab<S>,
    authored: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    original: Inputs<D>,
    applied: &[TransportSample<RheologyMaterial<F>>],
    gate: &CancelGate,
) -> Result<(Inputs<D>, ViscousHeating), MaterialTransportError> {
    let units = mapping.force_units;
    let volume = (units.cell_size_m * units.cell_size_m) * depth_m;
    let mut power_w_m3 = vec![0.0; authored.len()];
    let mut combined = authored.to_vec();
    let mut viscous_heat_j = 0.0;
    for value in applied {
        poll(gate)?;
        let rate = value.material.point().shear_rate_per_s;
        let viscosity = value.transport.dynamic_viscosity_pa_s();
        // This multiplication order avoids overflowing rate^2 before a small
        // viscosity can scale it. Both values come from the cached query.
        let power = (viscosity * rate) * rate;
        let energy = (power * volume) * units.timestep_s;
        if !power.is_finite()
            || power < 0.0
            || (rate != 0.0 && power == 0.0)
            || !energy.is_finite()
            || (power != 0.0 && energy == 0.0)
        {
            return Err(MaterialTransportError::Invalid(
                "physical viscous heating is unrepresentable",
            ));
        }
        power_w_m3[value.cell] = power;
        // Preserve the existing zero-source path, including signed zero.
        if power != 0.0 {
            combined[value.cell] += power;
        }
        viscous_heat_j += energy;
    }
    if !viscous_heat_j.is_finite() {
        return Err(MaterialTransportError::Invalid(
            "summed physical viscous heat is unrepresentable",
        ));
    }
    let authored_heat_j = original.prescribed_j;
    let input = mapping.heat_inputs(&combined, forces, authored.len(), depth_m, limits)?;
    let source_partition_defect_j = (input.prescribed_j - authored_heat_j) - viscous_heat_j;
    if !source_partition_defect_j.is_finite()
        || source_partition_defect_j.abs() > limits.balance_tolerance_j
    {
        return Err(MaterialTransportError::Invalid(
            "authored and viscous heat partition exceeds the requested joule tolerance",
        ));
    }
    Ok((
        input,
        ViscousHeating {
            power_w_m3,
            viscous_heat_j,
            authored_heat_j,
            source_partition_defect_j,
        },
    ))
}

#[allow(clippy::too_many_arguments)]
fn advance<S: StressSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    solve: RheologySolve,
    gate: &CancelGate,
    mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
    solve.validate()?;
    let (transport, viscous) =
        advance_with_sampler_and_sources::<S, RheologyMaterial<F>, ViscousHeating, D, Q>(
            slab,
            power,
            forces,
            depth_m,
            limits,
            gate,
            |mapping, state, cells, plane| {
                sample::<S, F, D, Q>(
                    mapping,
                    state,
                    cells,
                    plane,
                    limits.temperature_range_k,
                    solve,
                    gate,
                    &mut resolve,
                )
            },
            |mapping, original, applied| {
                source_inputs(
                    mapping, power, forces, depth_m, limits, original, applied, gate,
                )
            },
        )?;
    Ok(ViscousRheologyStep { transport, viscous })
}

impl MaterialSlab<ThermalLbm> {
    /// Advance with material viscosity AND its physical viscous heat.
    ///
    /// `power_w_m3` supplies external heat only: this method adds the molecular
    /// `mu * gamma_dot^2` loss once, using the exact initial viscosity/shear
    /// samples that drive momentum. The original two source half-kicks advance
    /// sensible heat, and endpoint material queries see the resulting actual
    /// temperature. Repeated calls close explicit temperature/viscosity feedback.
    /// Source, conversion, temperature, heat-balance, endpoint and cancellation
    /// refusals all preserve the accepted slab.
    ///
    /// The incompressible generalized-Newtonian, fixed-density/Cp, low-Mach
    /// assumptions of `step_with_material_rheology` still apply. Loss power is
    /// frozen at the initial state: this is not an implicit thermal solve or a
    /// second-order accuracy claim. The thermal audit does not prove transient
    /// total-energy conservation; numerical momentum dissipation and pressure
    /// work are neither measured nor converted into heat here.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart, with the same physical source, feedback and atomicity.
    /// Cell volume is the retained `dx^3`; see the 2D method for model limits.
    pub fn step_with_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }
}
