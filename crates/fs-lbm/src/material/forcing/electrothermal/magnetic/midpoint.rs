//! Coupled current, Lorentz-force and Joule-source feedback over one flow step.
//!
//! Each trial restarts from the accepted populations. The existing magnetic
//! source owner is evaluated at the old/new temperature AND velocity midpoint.
//! This removes source lag from the declared fixed-point equation; it does not
//! change the underlying Guo velocity convention or certify total-energy,
//! temporal order, unconditional stability, or magnetic induction.

use super::*;
use super::super::midpoint::MAX_MIDPOINT_ITERATIONS;

/// Independent nonlinear, coordinate, heat, force and work admission controls.
#[derive(Debug, Clone, Copy)]
pub struct MagneticMidpointOptions {
    /// Trial budget in `1..=MAX_MIDPOINT_ITERATIONS`; at most two current solves
    /// per trial, including the independent final source-consistency solve.
    pub max_iterations: usize,
    /// Maximum cellwise midpoint-temperature mismatch [K].
    pub temperature_tolerance_k: f64,
    /// Maximum absolute component of midpoint-velocity mismatch [m/s].
    pub velocity_tolerance_m_s: f64,
    /// L1 sum of cell Joule-power differences times dt [J].
    pub joule_tolerance_j: f64,
    /// L1 sum of force-density component differences times cell volume and dt
    /// [N s]. Opposing forces cannot cancel in this criterion.
    pub impulse_tolerance_n_s: f64,
    /// L1 sum of local Lorentz-work differences [J], evaluated using each
    /// solve's own force AND velocity, not just its aggregate mechanical power.
    pub work_tolerance_j: f64,
    /// Under-relaxation in `(0,1]`, applied to both coordinate fields. Changes
    /// iteration only; does not alter the accepted equations or timestep.
    pub relaxation: f64,
}

impl MagneticMidpointOptions {
    pub(super) fn admit(self) -> Result<(), MagneticMidpointError> {
        let nonnegative = |x: f64| x.is_finite() && x >= 0.0;
        if self.max_iterations == 0 || self.max_iterations > MAX_MIDPOINT_ITERATIONS
            || !positive(self.temperature_tolerance_k) || !positive(self.velocity_tolerance_m_s)
            || !nonnegative(self.joule_tolerance_j) || !nonnegative(self.impulse_tolerance_n_s)
            || !nonnegative(self.work_tolerance_j)
            || !positive(self.relaxation) || self.relaxation > 1.0
        { return Err(MagneticMidpointError::InvalidOptions); }
        Ok(())
    }
}

/// Independently re-solved source differences at the actual candidate midpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagneticSourceResidual {
    /// Sum of absolute local Joule-energy differences [J].
    pub joule_j: f64,
    /// Sum of absolute local impulse-component differences [N s].
    pub impulse_n_s: f64,
    /// Sum of absolute local magnetic-work differences [J].
    pub mechanical_work_j: f64,
}

/// Accepted coupled step; all returned sources are the ones actually applied.
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticMidpointStep {
    /// Source-producing current, conductor-frame heat, force and instantaneous power fields.
    pub motion: MovingConductorSolution,
    /// Actual heat storage and independently measured thermal wall exchange.
    pub heat: HeatBalance,
    /// Electrical work at the source evaluation point times dt [J].
    pub electrode_energy_j: f64,
    /// Applied conductor-frame Joule energy [J].
    pub joule_energy_j: f64,
    /// Lorentz power at `evaluation_velocity_m_s` times dt [J]. This is NOT
    /// the measured kinetic-energy change of the split fluid solver.
    pub midpoint_mechanical_work_j: f64,
    /// Heat source minus applied Joule energy [J].
    pub transfer_defect_j: f64,
    /// Stored heat minus wall inputs plus midpoint mechanical work minus
    /// electrode work [J]; a source-accounting check, not total-energy closure.
    pub coupling_defect_j: f64,
    /// Number of tentative thermal/flow steps, each starting at the old state.
    pub iterations: usize,
    /// Current solves, including consistency re-evaluations.
    pub electrical_solves: usize,
    /// Kelvin coordinates used by the applied source, including wall cells.
    pub evaluation_temperature_k: Vec<f64>,
    /// Physical velocity used by the applied source, including zero wall rows.
    pub evaluation_velocity_m_s: Vec<[f64; 3]>,
    /// Maximum temperature-coordinate consistency error [K].
    pub temperature_residual_k: f64,
    /// Maximum velocity-component consistency error [m/s].
    pub velocity_residual_m_s: f64,
    /// L1 source differences from a fresh solve at the true candidate midpoint.
    pub source_residual: MagneticSourceResidual,
}

/// No error commits any candidate population, force field or material unit.
#[derive(Debug, Clone, PartialEq)]
pub enum MagneticMidpointError {
    /// Invalid nonlinear budget, relaxation or tolerance.
    InvalidOptions,
    /// An existing material, electrical, thermal, Mach or cancellation refusal.
    Step(ElectricHeatError),
    /// Numerical iteration exhausted its budget; not a claim of physical runaway.
    NonConvergence {
        /// Number of attempted thermal/flow trials.
        iterations: usize,
        /// Last temperature-coordinate mismatch [K].
        temperature_residual_k: f64,
        /// Last velocity-coordinate mismatch [m/s].
        velocity_residual_m_s: f64,
        /// Last source recheck, or None if coordinate gates prevented that check.
        source_residual: Option<MagneticSourceResidual>,
    },
}
impl core::fmt::Display for MagneticMidpointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidOptions => write!(f, "invalid magnetothermal midpoint controls"),
            Self::Step(e) => write!(f, "{e}"),
            Self::NonConvergence { iterations, temperature_residual_k, velocity_residual_m_s, source_residual } =>
                write!(f, "magnetothermal midpoint did not converge in {iterations} trials: {temperature_residual_k} K, {velocity_residual_m_s} m/s, sources {source_residual:?}"),
        }
    }
}
impl std::error::Error for MagneticMidpointError {}
impl From<ElectricHeatError> for MagneticMidpointError {
    fn from(error: ElectricHeatError) -> Self { Self::Step(error) }
}
impl From<ThermalSourceError> for MagneticMidpointError {
    fn from(error: ThermalSourceError) -> Self { Self::Step(error.into()) }
}

// Check local work, not a net mechanical-power difference in which a motor
// region and a braking region could cancel. No new force or heat model here.
fn source_residual(
    used: &MovingConductorSolution, checked: &MovingConductorSolution,
    used_velocity: &[[f64; 3]], checked_velocity: &[[f64; 3]],
    domain: &VoxelDomain, dt: f64, gate: &CancelGate,
) -> Result<MagneticSourceResidual, ElectricHeatError> {
    let mut result = MagneticSourceResidual { joule_j: 0.0, impulse_n_s: 0.0, mechanical_work_j: 0.0 };
    for c in 0..domain.cell_count() {
        if c % 256 == 0 { poll(gate)?; }
        result.joule_j += (used.electric.joule_w[c] - checked.electric.joule_w[c]).abs() * dt;
        let volume = domain.volume(c);
        let mut used_power = 0.0;
        let mut checked_power = 0.0;
        for a in 0..3 {
            let fa = used.lorentz_force_density_n_m3[c][a];
            let fb = checked.lorentz_force_density_n_m3[c][a];
            result.impulse_n_s += ((fa - fb).abs() * volume) * dt;
            used_power += fa * used_velocity[c][a];
            checked_power += fb * checked_velocity[c][a];
        }
        result.mechanical_work_j += ((used_power - checked_power).abs() * volume) * dt;
    }
    if ![result.joule_j, result.impulse_n_s, result.mechanical_work_j].iter().all(|v| v.is_finite()) {
        return Err(invalid("magnetic midpoint source residual is unrepresentable"));
    }
    Ok(result)
}

// The original frozen-transport API uses the same iteration as local-material
// coupling. Its empty prepared state adds no material queries or physics steps.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve<S: MovingSlab<D, Q>, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]], field: &[[f64; 3]], external: &[[f64; D]],
    options: MagneticMidpointOptions, isolated_current_tolerance_a: Option<f64>, gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    mut step: impl FnMut(&mut MaterialSlab<S>, &[f64], &[[f64; D]], HeatSourceLimits) -> Result<HeatBalance, ThermalSourceError>,
) -> Result<MagneticMidpointStep, MagneticMidpointError> {
    let transport = slab.force_units.transport;
    let dt = slab.force_units.timestep_s;
    let result = solve_prepared::<S, (), (), MagneticMidpointError, D, Q>(
        slab, &domain, plane, input, field, external, options, gate,
        |_, _| Ok(()),
        |temperature, velocity, _| Ok(source_at_with_open_faces(
            &domain, plane, input, open_faces, field, temperature, velocity,
            &transport, dt, isolated_current_tolerance_a, gate, &mut resolve,
        )?),
        |old, power, force, _| {
            let mut next = candidate::<S, D, Q>(old);
            let heat = step(&mut next, power, force, input.limits.thermal)?;
            Ok((next, heat))
        },
        |_, _| Ok((true, ())),
    )?;
    Ok(result.step)
}

// Internal data retained by the two actual physics adapters, not another
// material system. The driver never queries or converts a prepared object.
pub(super) struct PreparedMidpointStep<P, R> {
    pub step: MagneticMidpointStep,
    pub applied: P,
    pub endpoint: P,
    pub material_residual: R,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn solve_prepared<
    S: MovingSlab<D, Q>, P, R, E, const D: usize, const Q: usize,
>(
    slab: &mut MaterialSlab<S>, domain: &VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field: &[[f64; 3]], external: &[[f64; D]],
    options: MagneticMidpointOptions, gate: &CancelGate,
    mut prepare: impl FnMut(&MaterialSlab<S>, &[f64]) -> Result<P, E>,
    mut source_at_prepared: impl FnMut(&[f64], &[[f64; 3]], &P) -> Result<MagneticSource, E>,
    mut step_prepared: impl FnMut(&MaterialSlab<S>, &[f64], &[[f64; D]], &P) -> Result<(MaterialSlab<S>, HeatBalance), E>,
    mut check_prepared: impl FnMut(&P, &P) -> Result<(bool, R), E>,
) -> Result<PreparedMidpointStep<P, R>, E>
where E: From<ElectricHeatError> + From<MagneticMidpointError> {
    options.admit()?;
    let cells = domain.cell_count();
    admit_field::<D>(field, cells)?;
    let old_temperature = temperatures::<S, D, Q>(slab, input.limits.thermal, gate)?;
    let old_velocity = velocities::<S, D, Q>(slab, cells, plane, gate)?;
    let mut evaluation_temperature = old_temperature.clone();
    let mut evaluation_velocity = old_velocity.clone();
    let dt = slab.force_units.timestep_s;
    let mut electrical_solves = 0;
    let mut last_temperature = 0.0;
    let mut last_velocity = 0.0;
    let mut last_source = None;
    for iteration in 1..=options.max_iterations {
        poll(gate)?;
        let applied = prepare(slab, &evaluation_temperature)?;
        let source = source_at_prepared(&evaluation_temperature, &evaluation_velocity, &applied)?;
        electrical_solves += 1;
        let forces = combined_forces(external, &source.motion, gate)?;
        let (next, heat) = step_prepared(slab, &source.power, &forces, &applied)?;
        let end_temperature = temperatures::<S, D, Q>(&next, input.limits.thermal, gate)?;
        let end_velocity = velocities::<S, D, Q>(&next, cells, plane, gate)?;
        let mut midpoint_temperature = Vec::with_capacity(cells);
        let mut midpoint_velocity = Vec::with_capacity(cells);
        last_temperature = 0.0_f64;
        last_velocity = 0.0_f64;
        for c in 0..cells {
            if c % 256 == 0 { poll(gate)?; }
            let t = 0.5 * old_temperature[c] + 0.5 * end_temperature[c];
            let u = core::array::from_fn::<_, 3, _>(|a| 0.5 * old_velocity[c][a] + 0.5 * end_velocity[c][a]);
            last_temperature = last_temperature.max((t - evaluation_temperature[c]).abs());
            for a in 0..3 { last_velocity = last_velocity.max((u[a] - evaluation_velocity[c][a]).abs()); }
            midpoint_temperature.push(t);
            midpoint_velocity.push(u);
        }
        last_source = None;
        if last_temperature <= options.temperature_tolerance_k && last_velocity <= options.velocity_tolerance_m_s {
            let checked_material = prepare(slab, &midpoint_temperature)?;
            let checked = source_at_prepared(&midpoint_temperature, &midpoint_velocity, &checked_material)?;
            electrical_solves += 1;
            let residual = source_residual(&source.motion, &checked.motion, &evaluation_velocity,
                &midpoint_velocity, &domain, dt, gate)?;
            last_source = Some(residual);
            let (material_matches, material_residual) = check_prepared(&applied, &checked_material)?;
            if material_matches && residual.joule_j <= options.joule_tolerance_j
                && residual.impulse_n_s <= options.impulse_tolerance_n_s
                && residual.mechanical_work_j <= options.work_tolerance_j
            {
                let report = balance(source, heat, input)?;
                // Endpoint evaluation is a support check, never a replacement
                // for the material that actually supplied the accepted sources.
                let endpoint = prepare(&next, &end_temperature)?;
                poll(gate)?;
                *slab = next;
                return Ok(PreparedMidpointStep { applied, endpoint, material_residual, step: MagneticMidpointStep {
                    motion: report.motion, heat: report.heat,
                    electrode_energy_j: report.electrode_energy_j, joule_energy_j: report.joule_energy_j,
                    midpoint_mechanical_work_j: report.lagged_mechanical_work_j,
                    transfer_defect_j: report.transfer_defect_j, coupling_defect_j: report.coupling_defect_j,
                    iterations: iteration, electrical_solves,
                    evaluation_temperature_k: evaluation_temperature,
                    evaluation_velocity_m_s: evaluation_velocity,
                    temperature_residual_k: last_temperature, velocity_residual_m_s: last_velocity,
                    source_residual: residual,
                }});
            }
        }
        for c in 0..cells {
            if c % 256 == 0 { poll(gate)?; }
            evaluation_temperature[c] += options.relaxation * (midpoint_temperature[c] - evaluation_temperature[c]);
            for a in 0..3 {
                evaluation_velocity[c][a] += options.relaxation * (midpoint_velocity[c][a] - evaluation_velocity[c][a]);
            }
        }
    }
    Err(MagneticMidpointError::NonConvergence {
        iterations: options.max_iterations, temperature_residual_k: last_temperature,
        velocity_residual_m_s: last_velocity, source_residual: last_source,
    }.into())
}

impl MaterialSlab<ThermalLbm> {
    /// Advance one 2-D step with coupled velocity/temperature magnetic feedback.
    /// Current, force and heat use the old/new midpoint; every trial reuses the
    /// old populations. Only out-of-plane Bz is admitted. Depth controls the
    /// electrical geometry and extensive energy/impulse accounts, not motion.
    ///
    /// Drives, B, contacts and retained thermal coefficients are fixed during
    /// the call. The constitutive callback must represent the same law on every
    /// query; captured receipts are caller-owned. Failed convergence, source
    /// admission or cancellation does not replace the accepted slab. Trial Mach
    /// limits still apply: this is not an unconditionally stable solver.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 2]], depth_m: f64, options: MagneticMidpointOptions,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        self.step_magnetothermal_midpoint_with_open_faces(
            input, &[], field_t, external_n_m3, depth_m, options, gate, resolve,
        )
    }

    /// The same coupled step with explicitly insulating internal electrical faces.
    /// `open_faces` names adjacent full-grid cell pairs, including the wall-row
    /// indexing. It is fixed through every trial and final source recheck.
    /// Electrical separation does not cut fluid momentum or thermal transport;
    /// each wired conducting component still needs a voltage reference.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_open_faces(
        &mut self, input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]], field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 2]], depth_m: f64, options: MagneticMidpointOptions,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, external_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume").into()); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ]).map_err(ElectricHeatError::from)?;
        solve::<ThermalLbm, 2, 5>(self, domain, plane, input, open_faces, field_t, external_n_m3, options, None, gate, resolve,
            |s,p,f,limits| s.step_with_power_density(p,f,depth_m,limits))
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart of the coupled midpoint step. All vector
    /// components participate, with the existing shared current/heat/force and
    /// thermal owners. Work is midpoint source power times dt, not a measured
    /// kinetic-energy change. No induction or total-energy/order claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 3]], options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        self.step_magnetothermal_midpoint_with_open_faces(
            input, &[], field_t, external_n_m3, options, gate, resolve,
        )
    }

    /// Three-dimensional electrical-face topology in the same midpoint solve.
    /// Open faces carry no current, motional work, contact loss or Lorentz
    /// transfer. Material volume, fluid connectivity and heat transport remain
    /// those of the original slab; topology is held fixed over this step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_open_faces(
        &mut self, input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]], field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 3]], options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, external_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx).map_err(ElectricHeatError::from)?;
        solve::<ThermalLbm3, 3, 7>(self, domain, plane, input, open_faces, field_t, external_n_m3, options, None, gate, resolve,
            |s,p,f,limits| s.step_with_power_density(p,f,limits))
    }
}
