//! Bounded temperature/resistivity feedback over one existing thermal-flow step.
//!
//! Solve T_next = thermal_step(T_old, Joule((T_old + T_next)/2)). Each trial
//! starts from the same accepted populations, not the previous trial. Midpoint
//! source evaluation is not a second-order claim for the underlying split LBM.

use super::*;

/// Hard ceiling on nonlinear trials in one call (each uses at most two current solves).
pub const MAX_MIDPOINT_ITERATIONS: usize = 128;

/// Controls for the source/temperature fixed point; independent of current
/// solver tolerance and thermal energy-balance admission.
#[derive(Debug, Clone, Copy)]
pub struct MidpointOptions {
    /// Nonlinear trial budget, in `1..=MAX_MIDPOINT_ITERATIONS`.
    pub max_iterations: usize,
    /// Maximum cellwise midpoint-temperature consistency error [K].
    pub temperature_tolerance_k: f64,
    /// L1 sum of absolute cell Joule-power discrepancies times dt [J], after
    /// re-solving at the actual candidate midpoint. Opposing local errors
    /// cannot cancel as they could in a comparison of total source power.
    pub source_tolerance_j: f64,
    /// Explicit under-relaxation in `(0,1]`. It changes iteration, not the
    /// accepted equations, energy tolerance, dt or external force schedule.
    pub relaxation: f64,
}
impl MidpointOptions {
    fn admit(self) -> Result<(), MidpointError> {
        if self.max_iterations == 0 || self.max_iterations > MAX_MIDPOINT_ITERATIONS
            || !positive(self.temperature_tolerance_k)
            || !self.source_tolerance_j.is_finite() || self.source_tolerance_j < 0.0
            || !positive(self.relaxation) || self.relaxation > 1.0
        { Err(MidpointError::InvalidOptions) } else { Ok(()) }
    }
}

/// Accepted step, including source evaluation coordinates and consistency residuals.
#[derive(Debug, Clone, PartialEq)]
pub struct MidpointStep {
    /// Fields and balances generated using `evaluation_temperature_k`.
    pub step: ElectricHeatStep,
    /// Actual nonlinear trials; each restarts from the same pre-step state.
    pub iterations: usize,
    /// Actual current solves, including final consistency re-evaluations.
    pub electrical_solves: usize,
    /// Kelvin field used by the accepted source solve (includes wall cells).
    pub evaluation_temperature_k: Vec<f64>,
    /// Maximum difference from the actual old/new temperature midpoint [K].
    pub temperature_residual_k: f64,
    /// Sum of absolute per-cell source discrepancies at that actual midpoint [J].
    pub joule_residual_j: f64,
}

/// No failure publishes a partial time step. Nonconvergence is a numerical
/// refusal, not proof of physical thermal runaway or absence of a solution.
#[derive(Debug, Clone, PartialEq)]
pub enum MidpointError {
    /// Invalid nonlinear budget, relaxation or tolerance.
    InvalidOptions,
    /// Material, geometry, current, heat, energy or cancellation refusal.
    Step(ElectricHeatError),
    /// Budget exhausted without satisfying both source and temperature tests.
    NonConvergence {
        /// Number of attempted nonlinear trials.
        iterations: usize,
        /// Last midpoint-temperature mismatch [K].
        temperature_residual_k: f64,
        /// Last trial's source mismatch [J], or None when its temperature
        /// residual was too large to attempt a source consistency solve.
        joule_residual_j: Option<f64>,
    },
}
impl core::fmt::Display for MidpointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidOptions => write!(f, "invalid electric/thermal midpoint controls"),
            Self::Step(e) => write!(f, "{e}"),
            Self::NonConvergence { iterations, temperature_residual_k, joule_residual_j } =>
                write!(f, "electric/thermal midpoint did not converge in {iterations} trials: {temperature_residual_k} K, source {joule_residual_j:?} J"),
        }
    }
}
impl std::error::Error for MidpointError {}
impl From<ElectricHeatError> for MidpointError {
    fn from(e: ElectricHeatError) -> Self { Self::Step(e) }
}
impl From<ThermalSourceError> for MidpointError {
    fn from(e: ThermalSourceError) -> Self { Self::Step(e.into()) }
}

#[allow(clippy::too_many_arguments)] // state, source, geometry and controls have separate owners
fn solve<S: Slab<D, Q>, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, forces: &[[f64; D]], options: MidpointOptions,
    gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    mut step: impl FnMut(&mut MaterialSlab<S>, &[f64], &[[f64; D]], HeatSourceLimits) -> Result<HeatBalance, ThermalSourceError>,
) -> Result<MidpointStep, MidpointError> {
    let old = temperatures::<S, D, Q>(slab, input.limits.thermal, gate)?;
    let mut evaluation = old.clone();
    let mut last_temperature = 0.0_f64;
    let mut last_joule = None;
    let mut electrical_solves = 0;
    let dt = slab.force_units.timestep_s;
    for iteration in 1..=options.max_iterations {
        poll(gate)?;
        let source = solve_source(&domain, plane, &evaluation, &slab.force_units.transport,
            dt, input, gate, &mut resolve)?;
        electrical_solves += 1;
        let mut next = candidate::<S, D, Q>(slab);
        let heat = step(&mut next, &source.power_w_m3, forces, input.limits.thermal)?;
        let end = temperatures::<S, D, Q>(&next, input.limits.thermal, gate)?;
        let exact_midpoint: Vec<_> = old.iter().zip(&end).map(|(a,b)| 0.5*a + 0.5*b).collect();
        last_temperature = evaluation.iter().zip(&exact_midpoint)
            .fold(0.0_f64, |m,(a,b)| m.max((a-b).abs()));
        last_joule = None;
        if last_temperature <= options.temperature_tolerance_k {
            let checked = solve_source(&domain, plane, &exact_midpoint, &slab.force_units.transport,
                dt, input, gate, &mut resolve)?;
            electrical_solves += 1;
            let discrepancy: f64 = source.electric.joule_w.iter().zip(&checked.electric.joule_w)
                .map(|(used,actual)| (used-actual).abs() * dt).sum();
            if !discrepancy.is_finite() { return Err(invalid("midpoint Joule-source residual is unrepresentable").into()); }
            last_joule = Some(discrepancy);
            if discrepancy <= options.source_tolerance_j {
                let report = balance(source, heat, input.limits)?;
                poll(gate)?;
                *slab = next;
                return Ok(MidpointStep {
                    step: report, iterations: iteration, electrical_solves,
                    evaluation_temperature_k: evaluation,
                    temperature_residual_k: last_temperature, joule_residual_j: discrepancy,
                });
            }
        }
        // Endpoints and their midpoint are admitted Kelvin values. Relaxation
        // remains inside their convex interval, with no extrapolated material query.
        for (guess, exact) in evaluation.iter_mut().zip(exact_midpoint) {
            *guess += options.relaxation * (exact - *guess);
        }
    }
    Err(MidpointError::NonConvergence {
        iterations: options.max_iterations, temperature_residual_k: last_temperature,
        joule_residual_j: last_joule,
    })
}

impl MaterialSlab<ThermalLbm> {
    /// Advance one 2-D step with bounded midpoint electrical-material feedback.
    /// Geometry, source law, wall/force conventions and rollback are those of
    /// `step_electrothermal`. The callback must evaluate the same material law
    /// on repeated queries; externally retained trial receipts remain caller-owned.
    /// Contacts and drives are held fixed over this call. Bulk resistance is
    /// re-evaluated; the slab's density/Cp/thermal transport are not reinterpreted.
    pub fn step_electrothermal_midpoint(
        &mut self, input: ElectricHeatInput<'_>, forces_n_m3: &[[f64; 2]], depth_m: f64,
        options: MidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MidpointStep, MidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume").into()); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ]).map_err(ElectricHeatError::from)?;
        solve::<ThermalLbm, 2, 5>(self, domain, plane, input, forces_n_m3, options, gate, resolve,
            |s,p,f,limits| s.step_with_power_density(p,f,depth_m,limits))
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart of the 2-D midpoint source/temperature solve.
    /// Reuses the same fixed-point, final source recheck and energy admission.
    pub fn step_electrothermal_midpoint(
        &mut self, input: ElectricHeatInput<'_>, forces_n_m3: &[[f64; 3]],
        options: MidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MidpointStep, MidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx).map_err(ElectricHeatError::from)?;
        solve::<ThermalLbm3, 3, 7>(self, domain, plane, input, forces_n_m3, options, gate, resolve,
            |s,p,f,limits| s.step_with_power_density(p,f,limits))
    }
}
