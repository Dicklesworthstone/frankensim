//! One within-step material state for resistance, thermal diffusion and momentum.
//!
//! Reuses the magnetic midpoint iteration and the local material sampler. Each
//! trial starts at the accepted populations, with all three properties evaluated
//! together at the proposed midpoint. No additional field solver is introduced.

use super::*;
use super::super::midpoint::{
    self as driver, MagneticMidpointError, MagneticMidpointOptions, MagneticMidpointStep,
};
use crate::material::transport::sample_temperatures;

/// Existing source/coordinate controls plus independent local transport closure.
#[derive(Debug, Clone, Copy)]
pub struct MaterialMidpointOptions {
    /// Shared temperature, velocity, heat, impulse, work and iteration controls.
    pub coupling: MagneticMidpointOptions,
    /// Largest allowed |used-checked|/max(used,checked) for either positive
    /// viscosity or thermal conductivity at any fluid cell, in [0,1]. Zero
    /// requests exact coefficient agreement. Not a material uncertainty bound.
    pub transport_relative_tolerance: f64,
}
impl MaterialMidpointOptions {
    fn admit(self) -> Result<(), MaterialMidpointError> {
        self.coupling.admit()?;
        if !self.transport_relative_tolerance.is_finite()
            || !(0.0..=1.0).contains(&self.transport_relative_tolerance)
        {
            return Err(MaterialMidpointError::InvalidTransportTolerance);
        }
        Ok(())
    }
}

/// Relative differences at the actual old/new midpoint, not aggregate averages.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportResidual {
    /// Maximum cellwise relative difference in thermal conductivity.
    pub conductivity_relative: f64,
    /// Maximum cellwise relative difference in dynamic viscosity.
    pub viscosity_relative: f64,
}

/// Accepted fields, applied material objects and separately checked endpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialMidpointStep<F> {
    /// Existing midpoint source fields, evaluation coordinates and balances.
    /// Mechanical work is evaluated source power times dt, not kinetic change.
    pub step: MagneticMidpointStep,
    /// Objects used for BOTH the applied current solve and local heat/flow rates.
    pub applied: Vec<TransportSample<ElectricalTransport<F>>>,
    /// Endpoint support checks, not replacements for the applied material.
    pub endpoint: Vec<TransportSample<ElectricalTransport<F>>>,
    /// Fresh local coefficient checks at the actual old/new midpoint.
    pub transport_residual: TransportResidual,
}

/// Every failure leaves the accepted slab intact. Externally captured query
/// receipts remain caller-owned, including those from rejected trials.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialMidpointError {
    /// Transport tolerance is not finite or lies outside [0,1].
    InvalidTransportTolerance,
    /// Original material, geometry, electrical, thermal or cancellation error.
    Step(MaterialMagneticError),
    /// Original midpoint control or iteration failure.
    Iteration(MagneticMidpointError),
    /// No trial passed every source, coordinate and local transport gate.
    NonConvergence {
        /// Original iteration count and last coordinate/source discrepancies.
        iteration: MagneticMidpointError,
        /// Most recent coefficient recheck, or None when no trial reached it.
        last_checked_transport: Option<TransportResidual>,
    },
}
impl core::fmt::Display for MaterialMidpointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidTransportTolerance => write!(f, "material midpoint transport tolerance must be in [0,1]"),
            Self::Step(e) => write!(f, "{e}"),
            Self::Iteration(e) => write!(f, "{e}"),
            Self::NonConvergence { iteration, last_checked_transport } =>
                write!(f, "{iteration}; last local transport check {last_checked_transport:?}"),
        }
    }
}
impl std::error::Error for MaterialMidpointError {}
impl From<MagneticMidpointError> for MaterialMidpointError {
    fn from(e: MagneticMidpointError) -> Self { Self::Iteration(e) }
}
impl From<MaterialTransportError> for MaterialMidpointError {
    fn from(e: MaterialTransportError) -> Self { Self::Step(e.into()) }
}
impl From<ElectricHeatError> for MaterialMidpointError {
    fn from(e: ElectricHeatError) -> Self { Self::Step(e.into()) }
}
impl From<ChtError> for MaterialMidpointError {
    fn from(e: ChtError) -> Self { Self::Step(e.into()) }
}
impl From<ThermalSourceError> for MaterialMidpointError {
    fn from(e: ThermalSourceError) -> Self { Self::Step(e.into()) }
}

type Samples<F> = Vec<TransportSample<ElectricalTransport<F>>>;

fn transport_residual<F>(
    used: &[TransportSample<ElectricalTransport<F>>],
    checked: &[TransportSample<ElectricalTransport<F>>], gate: &CancelGate,
) -> Result<TransportResidual, MaterialMidpointError> {
    if used.len() != checked.len() {
        return Err(invalid("midpoint material sample count changed").into());
    }
    let mut residual = TransportResidual { conductivity_relative: 0.0, viscosity_relative: 0.0 };
    for (index, (a, b)) in used.iter().zip(checked).enumerate() {
        if index % 256 == 0 { poll(gate)?; }
        if a.cell != b.cell { return Err(invalid("midpoint material cell order changed").into()); }
        // Both snapshots contain admitted positive finite coefficients. Taking
        // their positive difference before dividing cannot overflow.
        let relative = |x: f64, y: f64| (x - y).abs() / x.max(y);
        residual.conductivity_relative = residual.conductivity_relative.max(relative(
            a.transport.thermal_conductivity_w_m_k(), b.transport.thermal_conductivity_w_m_k(),
        ));
        residual.viscosity_relative = residual.viscosity_relative.max(relative(
            a.transport.dynamic_viscosity_pa_s(), b.transport.dynamic_viscosity_pa_s(),
        ));
    }
    Ok(residual)
}

#[allow(clippy::too_many_arguments)]
fn advance<S: MovingSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field: &[[f64; 3]], external: &[[f64; D]],
    depth_m: f64, isolated: Option<f64>, options: MaterialMidpointOptions, gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
    poll(gate)?;
    if let Some(tolerance) = isolated
        && (!input.setup.electrodes.is_empty() || !tolerance.is_finite() || tolerance < 0.0)
    {
        return Err(invalid("isolated material midpoint needs no electrodes and a finite nonnegative current tolerance").into());
    }
    let cells = domain.cell_count();
    let transport = slab.force_units.transport;
    let dt = slab.force_units.timestep_s;
    let mut last_transport = None;
    let result = driver::solve_prepared::<S, Samples<F>, TransportResidual, MaterialMidpointError, D, Q>(
        slab, &domain, plane, input, field, external, options.coupling, gate,
        |mapping, temperature| Ok(sample_temperatures::<S, ElectricalTransport<F>, D, Q>(
            mapping, temperature, cells, plane, input.limits.thermal.temperature_range_k, gate, &mut resolve,
        )?),
        |temperature, velocity, applied| {
            let mut electrical = |cell: usize, _: f64, _: &FluidTransportState| {
                Ok(applied[cell - plane].material.resistivity_ohm_m())
            };
            Ok(source_at(&domain, plane, input, field, temperature, velocity,
                &transport, dt, isolated, gate, &mut electrical)?)
        },
        |old, power, force, applied| {
            let (lattice, heat) = candidate_with_samples::<S, ElectricalTransport<F>, D, Q>(
                old, power, force, depth_m, input.limits.thermal, applied, gate,
            )?;
            Ok((MaterialSlab {
                lattice, top_temperature_k: old.top_temperature_k,
                temperature_difference_k: old.temperature_difference_k,
                force_units: old.force_units,
            }, heat))
        },
        |used, checked| {
            let residual = transport_residual(used, checked, gate)?;
            last_transport = Some(residual);
            Ok((residual.conductivity_relative <= options.transport_relative_tolerance
                && residual.viscosity_relative <= options.transport_relative_tolerance, residual))
        },
    );
    match result {
        Ok(r) => Ok(MaterialMidpointStep {
            step: r.step, applied: r.applied, endpoint: r.endpoint, transport_residual: r.material_residual,
        }),
        Err(MaterialMidpointError::Iteration(iteration @ MagneticMidpointError::NonConvergence { .. })) =>
            Err(MaterialMidpointError::NonConvergence { iteration, last_checked_transport: last_transport }),
        Err(error) => Err(error),
    }
}

impl MaterialSlab<ThermalLbm> {
    /// Couple local k(T), mu(T), resistivity(T), current, heat and force at one
    /// 2D old/new temperature/velocity midpoint. Every trial uses the same old
    /// populations. Original concrete material objects supply all conversions
    /// once per query, including fluid electrode cells; wall reservoirs are not
    /// queried. Endpoint material admission precedes atomic publication.
    ///
    /// Sources use the existing current/contact/force owners. Bz [T], electrode
    /// drives, external forces [N/m3], contacts, density and Cp stay fixed over
    /// the call. Explicit depth [m] controls extensive energy/impulse, not 2D
    /// dynamics. Thermal and flow coefficients are re-evaluated each trial.
    /// Coordinate, source AND independent transport gates must all pass.
    ///
    /// No induction, latent heat, exact kinetic-plus-thermal conservation,
    /// second-order split integration, unconditional stability or calibrated
    /// material claim. Midpoint and endpoint queries do not prove continuous
    /// support between those points. Captured query side effects remain owned
    /// by the caller; cancellation is polled between queries/whole kernel steps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]],
        forces: &[[f64; 2]], depth_m: f64, options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        self.material_midpoint(input, field, forces, depth_m, None, options, gate, resolve)
    }

    /// Electrode-free counterpart: gauges retain their complete physical
    /// contributions. Requires an empty electrode list and current tolerance [A].
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]],
        forces: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64,
        options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        self.material_midpoint(input, field, forces, depth_m, Some(current_tolerance_a), options, gate, resolve)
    }

    #[allow(clippy::too_many_arguments)]
    fn material_midpoint<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 2]],
        depth_m: f64, isolated: Option<f64>, options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) {
            return Err(invalid("invalid represented 2D cell volume").into());
        }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, F, 2, 5>(self, domain, plane, input, field, forces, depth_m, isolated, options, gate, resolve)
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart using the SAME midpoint driver and local material sampler,
    /// actual D3Q7 scaling and retained cubic volumes. All vector components
    /// participate. The fixed-storage and source-work limitations remain.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 3]],
        options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        self.material_midpoint(input, field, forces, None, options, gate, resolve)
    }

    /// 3D isolated-conductor counterpart with a local current-continuity limit [A].
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 3]],
        current_tolerance_a: f64, options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        self.material_midpoint(input, field, forces, Some(current_tolerance_a), options, gate, resolve)
    }

    #[allow(clippy::too_many_arguments)]
    fn material_midpoint<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 3]],
        isolated: Option<f64>, options: MaterialMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMidpointStep<F>, MaterialMidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, F, 3, 7>(self, domain, plane, input, field, forces, dx, isolated, options, gate, resolve)
    }
}
