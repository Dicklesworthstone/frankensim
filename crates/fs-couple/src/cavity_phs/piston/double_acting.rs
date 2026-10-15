//! Two finite, material-resolved gas volumes driving ONE translating wall.
//!
//! This is the missing mechanical connection between closed chambers, not two
//! independently stepped pistons. The original gas midpoint law and piston
//! impulse solver are reused. Each gas keeps its own composition, source cards,
//! caloric state, transport, pressure and acoustic properties.

use super::{
    ClosedGasParcel, GasEnergySolve, GasHeatVolumeCandidate, GasHeatVolumeStep,
    GasPiston, PistonError, PistonInput, PistonSolve, ResolvedThermochemicalGasState,
    Trial, impulse, invalid, poll, positive,
};

/// Two opposed constant-area faces on the SAME rigid translating body.
///
/// Positive velocity expands the primary gas and compresses the opposing gas.
/// Only the primary owner holds mass and velocity. The opposing gas has no
/// fictitious independent inertia or motion coordinate. Different positive face
/// areas are allowed (e.g. a piston rod); any other load is explicit input.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleActingGasPiston {
    primary: GasPiston,
    opposing: ClosedGasParcel,
    opposing_area_m2: f64,
}

/// External energy and force over one common physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoublePistonInput {
    /// Positive elapsed time [s].
    pub duration_s: f64,
    /// Integrated heat into [primary, opposing] gas [J], not watts or reactions.
    pub heat_j: [f64; 2],
    /// Additional constant force [N], positive on primary expansion.
    /// Must exclude BOTH gas pressure loads, which this owner supplies.
    pub external_force_n: f64,
}

/// Original nonlinear controls plus the second material's own support bracket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoublePistonSolve {
    /// Primary volume bracket and gas controls; the energy tolerance now bounds
    /// BOTH gas energies plus the SINGLE wall kinetic energy and external work.
    pub piston: PistonSolve,
    /// Opposing gas temperature/energy/source work controls, never borrowed from
    /// the primary material. Both gases must admit every requested proposal.
    pub opposing_gas: GasEnergySolve,
    /// Absolute allowed roundoff in dV_opposing + A_opposing*dx [m3].
    /// This does not relax the independently checked joule/impulse balances.
    pub volume_tolerance_m3: f64,
}

/// Joint accepted state accounting. There is NO exterior-pressure supply term:
/// pressure work between these two finite gases and their wall is internal.
#[derive(Clone, Debug, PartialEq)]
pub struct DoublePistonStep {
    /// Original gas heat/work/entropy reports, [primary, opposing].
    pub gas: [GasHeatVolumeStep; 2],
    /// Source-bearing pressure/thermal/acoustic midpoint states in the same order.
    pub midpoint: [ResolvedThermochemicalGasState; 2],
    /// Common signed wall displacement [m].
    pub displacement_m: f64,
    /// Initial and accepted final wall velocity [m/s].
    pub velocity_m_s: [f64; 2],
    /// Independently measured change of the ONE wall's kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Prescribed extra force times the common displacement [J].
    pub external_force_work_j: f64,
    /// M*dv - dt*(A_primary*p_primary - A_opposing*p_opposing + F) [N s].
    pub momentum_residual_n_s: f64,
    /// dK - W_primary - W_opposing - W_external [J].
    pub mechanical_energy_residual_j: f64,
    /// dU_primary + dU_opposing + dK - Q_primary - Q_opposing - W_external [J].
    pub total_energy_residual_j: f64,
    /// Actual dV_opposing + A_opposing*dx [m3].
    pub volume_exchange_residual_m3: f64,
    /// Interior joint momentum iterations, excluding initial/bracket proposals.
    pub iterations: usize,
    /// Total gas proposals, including discarded trials: TWO per impulse trial.
    /// Each retains its own bounded caloric inversion and source resolution.
    pub material_trials: usize,
}

struct PairedTrial {
    primary: Trial,
    opposing: GasHeatVolumeCandidate,
    primary_input: PistonInput,
    volume_residual: f64,
}

impl DoubleActingGasPiston {
    /// Bind an existing accepted piston to a second closed material volume.
    /// No source is resampled, no gas state is reset, and no second body mass
    /// can be supplied. This constructs a new coupling, not an impulse-free
    /// connection of two formerly moving mechanical bodies.
    pub fn new(
        primary: GasPiston, opposing: ClosedGasParcel, opposing_area_m2: f64,
    ) -> Result<Self, PistonError> {
        if !positive(opposing_area_m2)
            || !positive(opposing.volume_m3() / opposing_area_m2)
            || !positive(opposing_area_m2 / primary.area_m2())
        {
            return Err(invalid("opposed gas needs a positive represented face area and cavity length"));
        }
        Ok(Self { primary, opposing, opposing_area_m2 })
    }

    /// Original primary dynamics owner. Immutable: a caller cannot advance only
    /// one side and leave the opposing volume or shared momentum behind.
    #[must_use]
    pub const fn primary(&self) -> &GasPiston { &self.primary }
    /// Original opposing material state, conserved gas mass and source context.
    #[must_use]
    pub const fn opposing(&self) -> &ClosedGasParcel { &self.opposing }
    /// Positive pressure-exposed area of the opposing face [m2].
    #[must_use]
    pub const fn opposing_area_m2(&self) -> f64 { self.opposing_area_m2 }

    fn opposing_volume(&self, primary_volume: f64) -> Result<f64, PistonError> {
        let displacement = (primary_volume - self.primary.gas().volume_m3()) / self.primary.area_m2();
        let volume = self.opposing.volume_m3() - self.opposing_area_m2 * displacement;
        if !displacement.is_finite() || !positive(volume)
            || !positive(volume / self.opposing_area_m2)
        {
            return Err(invalid("primary volume proposal collapses or overflows the opposing cavity"));
        }
        Ok(volume)
    }

    fn trial(
        &self, volume: f64, input: DoublePistonInput, controls: DoublePistonSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<impulse::Trial<PairedTrial>, PistonError> {
        poll(cancelled)?;
        let opposing_volume = self.opposing_volume(volume)?;
        let opposing = self.opposing.midpoint_heat_volume(
            opposing_volume, input.heat_j[1], controls.opposing_gas, &mut *cancelled,
        )?;
        let ratio = self.opposing_area_m2 / self.primary.area_m2();
        // The original one-sided trial already owns its gas/momentum equations.
        // Its exterior load is now a material proposal, NOT a fixed reservoir.
        let exterior = opposing.midpoint().state().pressure * ratio;
        if !positive(exterior) {
            return Err(invalid("opposing pressure load is unrepresentable"));
        }
        let primary_input = PistonInput {
            duration_s: input.duration_s, heat_j: input.heat_j[0],
            exterior_pressure_pa: exterior, external_force_n: input.external_force_n,
        };
        let primary = self.primary.trial(volume, primary_input, controls.piston, cancelled)?;
        // dV_opposing/dV_primary=-A_opposing/A_primary. Both material tangents
        // are of the SAME discrete heat-volume update, not -gamma*p/V.
        let tangent = primary.tangent - input.duration_s * self.opposing_area_m2 * ratio
            * opposing.step().pressure_final_volume_derivative_pa_m3;
        let volume_residual = (opposing.endpoint().volume_m3() - self.opposing.volume_m3())
            + self.opposing_area_m2 * primary.displacement;
        if !positive(tangent) || !volume_residual.is_finite()
            || volume_residual.abs() > controls.volume_tolerance_m3
        {
            return Err(invalid("shared-piston discrete tangent or volume exchange failed admission"));
        }
        Ok(impulse::Trial {
            residual: primary.residual, tangent,
            value: PairedTrial { primary, opposing, primary_input, volume_residual },
        })
    }

    /// Jointly solve both gas energies and the single wall momentum.
    ///
    /// dV_primary=A_primary*dx, dV_opposing=-A_opposing*dx,
    /// dx=dt*(v0+v1)/2, M*dv=dt*(A_primary*pm_primary-A_opposing*pm_opposing+F).
    /// Every proposal starts from the same two accepted material states.
    /// Heat, pressure, transport and acoustic response are recomputed by each
    /// original gas owner; different species or caloric data need no solver edit.
    ///
    /// Both geometric endpoint volumes must be positive BEFORE querying sources.
    /// At most 2*(maximum_iterations+3) gas proposals are made. Sources, geometry,
    /// balances, exhausted work and cancellation cannot publish either side alone.
    /// There are no valves, mass/species exchange, contact stops, crank linkages,
    /// combustion, wall heat capacity or acoustic propagation in this component.
    pub fn step(
        &mut self, input: DoublePistonInput, controls: DoublePistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonStep, PistonError> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.piston.volume_bracket_m3;
        if !positive(input.duration_s) || !input.heat_j.iter().all(|q| q.is_finite())
            || !input.external_force_n.is_finite() || !positive(controls.volume_tolerance_m3)
            || !positive(controls.piston.energy_tolerance_j)
            || !positive(lo) || !positive(hi) || lo >= hi
            || self.primary.gas().volume_m3() < lo || self.primary.gas().volume_m3() > hi
        {
            return Err(invalid("invalid shared-piston input, volume bracket or balance tolerances"));
        }
        self.opposing_volume(lo)?;
        self.opposing_volume(hi)?;
        let solved = impulse::solve(self.primary.gas().volume_m3(), controls.piston, |volume| {
            self.trial(volume, input, controls, &mut cancelled)
        })?;
        let PairedTrial { primary, opposing, primary_input, volume_residual } = solved.value;
        // Finish the original owner only on a private candidate. The accepted
        // pair remains untouched even if the final all-owner audit is refused.
        let opposing_midpoint = opposing.midpoint().clone();
        let mut next_primary = self.primary.clone();
        let primary_step = next_primary.finish(primary, primary_input, controls.piston,
            solved.iterations, solved.proposals, &mut cancelled)?;
        let (next_opposing, opposing_step) = opposing.into_parts();
        let force = self.primary.area_m2() * primary_step.midpoint.state().pressure
            - self.opposing_area_m2 * opposing_midpoint.state().pressure + input.external_force_n;
        let momentum = self.primary.mass_kg()
            * (primary_step.velocity_m_s[1] - primary_step.velocity_m_s[0]) - input.duration_s * force;
        let mechanical = (primary_step.kinetic_energy_change_j - primary_step.gas.work_by_gas_j)
            - opposing_step.work_by_gas_j - primary_step.external_force_work_j;
        let total = ((primary_step.gas.stored_energy_change_j + opposing_step.stored_energy_change_j)
            + primary_step.kinetic_energy_change_j) - input.heat_j[0] - input.heat_j[1]
            - primary_step.external_force_work_j;
        if ![force, momentum, mechanical, total].iter().all(|x| x.is_finite())
            || momentum.abs() > controls.piston.momentum_tolerance_n_s
            || mechanical.abs() > controls.piston.energy_tolerance_j
            || total.abs() > controls.piston.energy_tolerance_j
        {
            return Err(invalid("joint gas/piston momentum or independently measured total energy failed admission"));
        }
        let report = DoublePistonStep {
            gas: [primary_step.gas, opposing_step],
            midpoint: [primary_step.midpoint, opposing_midpoint],
            displacement_m: primary_step.displacement_m,
            velocity_m_s: primary_step.velocity_m_s,
            kinetic_energy_change_j: primary_step.kinetic_energy_change_j,
            external_force_work_j: primary_step.external_force_work_j,
            momentum_residual_n_s: momentum, mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total, volume_exchange_residual_m3: volume_residual,
            iterations: solved.iterations, material_trials: 2 * solved.proposals,
        };
        poll(&mut cancelled)?;
        self.primary = next_primary;
        self.opposing = next_opposing;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
