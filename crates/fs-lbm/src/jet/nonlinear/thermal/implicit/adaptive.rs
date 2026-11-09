//! Heat-aware adaptive coupled implicit integration.
//!
//! Each comparison uses one selected-method joint step and two joint half steps.
//! Only the fine path is retained; neither state nor physical heat is Richardson
//! extrapolated. Endpoint and heat differences are LOCAL indicators, not global
//! error certificates. Their relative scales use changes during this interval,
//! rather than a thermal datum, absolute Kelvin temperature or bulk translation.
//!
//! A separate accepted-prefix ledger enforces D <= D_abs + r*(H+D), where H is
//! physical viscous heat and D is backward-Euler kinetic damping. Numerical loss
//! is never heated. This accounting allowance does not bound all heat error:
//! transport, temporal quadrature and nonlinear residuals still matter.
//!
//! Reusable O(N) checkpoints preserve the whole-call transaction. Coarse and
//! rejected solves count against the same work budgets; their heat is discarded.

use fs_math::dd::Dd;

use super::*;

/// Absolute SI allowance plus a dimensionless relative allowance.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetAdaptiveTolerance {
    pub absolute: f64,
    pub relative: f64,
}
impl ThermalJetAdaptiveTolerance {
    fn valid(self) -> bool {
        self.absolute.is_finite()
            && self.absolute >= 0.0
            && self.relative.is_finite()
            && self.relative >= 0.0
    }
    fn ratio(self, error: f64, scale: f64) -> Result<f64, ThermalJetAdaptiveError> {
        Ok(ratio(error, finite(self.absolute + self.relative * scale)?))
    }
}

/// Local physical indicators and a separate whole-call damping allowance.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetAdaptiveControls {
    /// Minimum comparison interval [s]; each retained method step is half of this.
    /// A shorter final remainder may be tried, but not refined below itself.
    pub minimum_interval_s: f64,
    /// Consecutive halvings allowed at one accepted state, at most 52.
    pub maximum_refinements: usize,
    pub area_m2: ThermalJetAdaptiveTolerance,
    pub radius_m: ThermalJetAdaptiveTolerance,
    pub velocity_m_s: ThermalJetAdaptiveTolerance,
    pub temperature_k: ThermalJetAdaptiveTolerance,
    pub specific_energy_j_kg: ThermalJetAdaptiveTolerance,
    pub cell_heat_j: ThermalJetAdaptiveTolerance,
    /// Controls the sum of absolute cell-heat discrepancies, so opposite local
    /// errors cannot cancel in an aggregate heat comparison.
    pub total_heat_j: ThermalJetAdaptiveTolerance,
    pub absolute_numerical_damping_j: f64,
    /// Finite in [0,1): D <= absolute + relative*(H+D) at every accepted prefix.
    pub relative_numerical_damping: f64,
}
impl Default for ThermalJetAdaptiveControls {
    fn default() -> Self {
        let tolerance = |absolute| ThermalJetAdaptiveTolerance {
            absolute,
            relative: 1e-3,
        };
        Self {
            minimum_interval_s: 1e-12,
            maximum_refinements: 40,
            area_m2: tolerance(1e-10),
            radius_m: tolerance(1e-10),
            velocity_m_s: tolerance(1e-6),
            temperature_k: tolerance(1e-4),
            specific_energy_j_kg: tolerance(1e-4),
            cell_heat_j: tolerance(1e-8),
            total_heat_j: tolerance(1e-7),
            absolute_numerical_damping_j: 1e-8,
            relative_numerical_damping: 1e-3,
        }
    }
}
impl ThermalJetAdaptiveControls {
    fn validate(self) -> Result<(), ThermalJetAdaptiveError> {
        if !self.minimum_interval_s.is_finite()
            || self.minimum_interval_s <= 0.0
            || self.maximum_refinements > 52
            || ![
                self.area_m2,
                self.radius_m,
                self.velocity_m_s,
                self.temperature_k,
                self.specific_energy_j_kg,
                self.cell_heat_j,
                self.total_heat_j,
            ]
            .iter()
            .all(|t| t.valid())
            || !self.absolute_numerical_damping_j.is_finite()
            || self.absolute_numerical_damping_j < 0.0
            || !self.relative_numerical_damping.is_finite()
            || !(0.0..1.0).contains(&self.relative_numerical_damping)
        {
            return Err(ThermalJetAdaptiveError::Invalid(
                "finite physical tolerances, bounded refinement and damping fraction in [0,1) required",
            ));
        }
        Ok(())
    }
}

/// Local fine-path estimates |fine-coarse|/(2^p-1), with method order p.
/// The largest normalized indicator uses changes during the trial interval.
#[derive(Clone, Copy, Debug, Default)]
pub struct ThermalJetAdaptiveIndicators {
    pub area_m2: f64,
    pub radius_m: f64,
    pub velocity_m_s: f64,
    pub temperature_k: f64,
    pub specific_energy_j_kg: f64,
    pub cell_heat_j: f64,
    pub total_heat_j: f64,
    pub scaled_maximum: f64,
}
impl ThermalJetAdaptiveIndicators {
    fn for_method(mut self, method: ThermalJetImplicitMethod) -> Self {
        if method == ThermalJetImplicitMethod::StrangMidpoint {
            self.area_m2 /= 3.0;
            self.radius_m /= 3.0;
            self.velocity_m_s /= 3.0;
            self.temperature_k /= 3.0;
            self.specific_energy_j_kg /= 3.0;
            self.cell_heat_j /= 3.0;
            self.total_heat_j /= 3.0;
            self.scaled_maximum /= 3.0;
        }
        self
    }

    fn merge(&mut self, other: Self) {
        self.area_m2 = self.area_m2.max(other.area_m2);
        self.radius_m = self.radius_m.max(other.radius_m);
        self.velocity_m_s = self.velocity_m_s.max(other.velocity_m_s);
        self.temperature_k = self.temperature_k.max(other.temperature_k);
        self.specific_energy_j_kg = self.specific_energy_j_kg.max(other.specific_energy_j_kg);
        self.cell_heat_j = self.cell_heat_j.max(other.cell_heat_j);
        self.total_heat_j = self.total_heat_j.max(other.total_heat_j);
        self.scaled_maximum = self.scaled_maximum.max(other.scaled_maximum);
    }
}

#[derive(Debug)]
pub enum ThermalJetAdaptiveError {
    Invalid(&'static str),
    Implicit(ThermalJetImplicitError),
    Resolution {
        interval_s: f64,
        refinements: usize,
        local: Option<ThermalJetAdaptiveIndicators>,
        damping_ratio: f64,
    },
}
impl From<ThermalJetImplicitError> for ThermalJetAdaptiveError {
    fn from(value: ThermalJetImplicitError) -> Self {
        Self::Implicit(value)
    }
}
impl From<NonlinearJetError> for ThermalJetAdaptiveError {
    fn from(value: NonlinearJetError) -> Self {
        Self::Implicit(value.into())
    }
}
impl core::fmt::Display for ThermalJetAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thermoviscous jet: {self:?}")
    }
}
impl std::error::Error for ThermalJetAdaptiveError {}

#[derive(Clone, Copy, Debug)]
pub struct ThermalJetAdaptiveAdvance {
    /// Actual selected-method fine trajectory. Attempted substeps and nonlinear work
    /// include coarse/rejected trials; accepted substeps are the two half steps.
    pub evolution: ThermalJetImplicitAdvance,
    pub accepted_intervals: usize,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_local: ThermalJetAdaptiveIndicators,
    /// Includes the absolute-only reserve scheduling indicator when applicable.
    pub maximum_accepted_damping_ratio: f64,
    pub numerical_damping_allowance_j: f64,
    pub minimum_accepted_interval_s: f64,
    pub maximum_accepted_interval_s: f64,
}

struct Checkpoint {
    fields: Fields,
    energy: Vec<f64>,
    temperature: Vec<f64>,
}
impl Checkpoint {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            fields: Fields::zeros(n, gate)?,
            energy: zeros(n, gate)?,
            temperature: zeros(n, gate)?,
        })
    }
    fn save(
        &mut self,
        fields: &Fields,
        energy: &[f64],
        profile: &CarreauProfile,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        for i in 0..self.energy.len() {
            poll_at(gate, i)?;
            self.fields.a[i] = fields.a[i];
            self.fields.q[i] = fields.q[i];
            self.energy[i] = energy[i];
            self.temperature[i] = profile.temperature_k(i);
        }
        Ok(())
    }
    fn restore(
        &self,
        jet: &mut ThermalCarreauJet<'_>,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        for i in 0..self.energy.len() {
            poll_at(gate, i)?;
            jet.flow.candidate.a[i] = self.fields.a[i];
            jet.flow.candidate.q[i] = self.fields.q[i];
            jet.candidate[i] = self.energy[i];
        }
        Ok(())
    }
}

pub(super) struct Workspace {
    base: Checkpoint,
    coarse: Checkpoint,
    coarse_heat: Vec<f64>,
}
impl Workspace {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            base: Checkpoint::new(n, gate)?,
            coarse: Checkpoint::new(n, gate)?,
            coarse_heat: zeros(n, gate)?,
        })
    }
}

/// Dd keeps a small accepted increment visible beside a long accumulated history.
#[derive(Clone, Copy)]
struct Ledger {
    heat: Dd,
    damping: Dd,
}
impl Default for Ledger {
    fn default() -> Self {
        Self {
            heat: Dd::ZERO,
            damping: Dd::ZERO,
        }
    }
}
impl Ledger {
    fn append(self, heat: f64, damping: f64) -> Self {
        Self {
            heat: self.heat + Dd::from_f64(heat),
            damping: self.damping + Dd::from_f64(damping),
        }
    }
    fn remaining(self, c: ThermalJetAdaptiveControls) -> Result<f64, ThermalJetAdaptiveError> {
        let r = Dd::from_f64(c.relative_numerical_damping);
        let unused = Dd::from_f64(c.absolute_numerical_damping_j) + r * self.heat
            - (Dd::ONE - r) * self.damping;
        let value = finite(unused.to_f64())?;
        if unused.lt(Dd::ZERO) {
            return Err(ThermalJetAdaptiveError::Invalid(
                "accepted damping prefix exceeded its allowance",
            ));
        }
        Ok(value)
    }
    fn allowed(self, c: ThermalJetAdaptiveControls) -> Result<bool, ThermalJetAdaptiveError> {
        let h = finite(self.heat.to_f64())?;
        let d = finite(self.damping.to_f64())?;
        let r = Dd::from_f64(c.relative_numerical_damping);
        let excess = (Dd::ONE - r) * self.damping
            - r * self.heat
            - Dd::from_f64(c.absolute_numerical_damping_j);
        let allowance =
            finite(c.absolute_numerical_damping_j + c.relative_numerical_damping * (h + d))?;
        Ok(!Dd::ZERO.lt(excess) && d <= allowance)
    }
}

/// A positive sub-ULP cost must not disappear by first adding it to the history.
fn ratio(cost: f64, allowance: f64) -> f64 {
    if cost == 0.0 {
        0.0
    } else if allowance == 0.0 {
        f64::INFINITY
    } else {
        cost / allowance
    }
}

/// Centered kinetic energy without subtracting total and bulk kinetic energies.
fn centered_kinetic(
    fields: &Fields,
    mass_factor: f64,
    gate: &CancelGate,
) -> Result<f64, NonlinearJetError> {
    let n = fields.a.len();
    let first = Dd::from_f64(fields.q[0]) / Dd::from_f64(average(fields.a[0], fields.a[1]));
    let mut weights = Sum::default();
    let mut mean_change = Dd::ZERO;
    for i in 0..n {
        poll_at(gate, i)?;
        let w = positive(average(fields.a[i], fields.a[next(i, n)]))?;
        weights.add(w);
        mean_change = mean_change + Dd::from_f64(fields.q[i]) - Dd::from_f64(w) * first;
    }
    let mean_change = mean_change / Dd::from_f64(positive(weights.total())?);
    let mut energy = Sum::default();
    for i in 0..n {
        poll_at(gate, i)?;
        let w = average(fields.a[i], fields.a[next(i, n)]);
        let centered =
            finite((Dd::from_f64(fields.q[i]) / Dd::from_f64(w) - first - mean_change).to_f64())?;
        energy.add(finite(0.5 * mass_factor * w * centered * centered)?);
    }
    finite(energy.total())
}

fn component(
    start: Dd,
    coarse: Dd,
    fine: Dd,
    tolerance: ThermalJetAdaptiveTolerance,
    error_maximum: &mut f64,
    scaled_maximum: &mut f64,
) -> Result<(), ThermalJetAdaptiveError> {
    let error = finite((fine - coarse).to_f64())?.abs();
    let change = finite((fine - start).to_f64())?
        .abs()
        .max(finite((coarse - start).to_f64())?.abs());
    *error_maximum = error_maximum.max(error);
    *scaled_maximum = scaled_maximum.max(tolerance.ratio(error, change)?);
    Ok(())
}

impl Workspace {
    fn indicators(
        &self,
        jet: &ThermalCarreauJet<'_>,
        coarse_heat: f64,
        fine_heat: f64,
        c: ThermalJetAdaptiveControls,
        gate: &CancelGate,
    ) -> Result<ThermalJetAdaptiveIndicators, ThermalJetAdaptiveError> {
        let n = self.base.energy.len();
        let mass_factor = positive(
            core::f64::consts::PI * jet.flow.parameters.density_kg_m3 * jet.flow.spacing_m,
        )?;
        let mut result = ThermalJetAdaptiveIndicators::default();
        let mut heat_difference = Sum::default();
        let fine_local_heat = jet.flow.spatial.trial_viscous_heat_j()?;
        let dd = Dd::from_f64;
        for i in 0..n {
            poll_at(gate, i)?;
            let a0 = self.base.fields.a[i];
            let ac = self.coarse.fields.a[i];
            let af = jet.flow.stage.a[i];
            component(
                dd(a0),
                dd(ac),
                dd(af),
                c.area_m2,
                &mut result.area_m2,
                &mut result.scaled_maximum,
            )?;
            component(
                dd(fs_math::det::sqrt(a0)),
                dd(fs_math::det::sqrt(ac)),
                dd(fs_math::det::sqrt(af)),
                c.radius_m,
                &mut result.radius_m,
                &mut result.scaled_maximum,
            )?;
            let face = |fields: &Fields| average(fields.a[i], fields.a[next(i, n)]);
            component(
                dd(self.base.fields.q[i]) / dd(face(&self.base.fields)),
                dd(self.coarse.fields.q[i]) / dd(face(&self.coarse.fields)),
                dd(jet.flow.stage.q[i]) / dd(face(&jet.flow.stage)),
                c.velocity_m_s,
                &mut result.velocity_m_s,
                &mut result.scaled_maximum,
            )?;
            component(
                dd(self.base.temperature[i]),
                dd(self.coarse.temperature[i]),
                dd(jet.workspace.profile.temperature_k(i)),
                c.temperature_k,
                &mut result.temperature_k,
                &mut result.scaled_maximum,
            )?;
            component(
                dd(self.base.energy[i]) / dd(positive(mass_factor * a0)?),
                dd(self.coarse.energy[i]) / dd(positive(mass_factor * ac)?),
                dd(jet.stage[i]) / dd(positive(mass_factor * af)?),
                c.specific_energy_j_kg,
                &mut result.specific_energy_j_kg,
                &mut result.scaled_maximum,
            )?;
            let error = finite(fine_local_heat[i] - self.coarse_heat[i])?.abs();
            let scale = fine_local_heat[i].max(self.coarse_heat[i]);
            result.cell_heat_j = result.cell_heat_j.max(error);
            result.scaled_maximum = result
                .scaled_maximum
                .max(c.cell_heat_j.ratio(error, scale)?);
            heat_difference.add(error);
        }
        result.total_heat_j =
            finite(heat_difference.total())?.max(finite(fine_heat - coarse_heat)?.abs());
        result.scaled_maximum = result.scaled_maximum.max(
            c.total_heat_j
                .ratio(result.total_heat_j, fine_heat.max(coarse_heat))?,
        );
        Ok(result)
    }
}

struct Comparison {
    local: ThermalJetAdaptiveIndicators,
    heat: f64,
    upwind: f64,
    damping: f64,
}

impl ThermalCarreauJet<'_> {
    /// Adaptive first-order joint implicit evolution with physical-heat control.
    ///
    /// Only two-half-step trajectories are published. The exact kernel retains
    /// its transport/capillary CFL and original source bounds. The substep budget
    /// counts ALL coarse and fine trials, including discarded ones. Unsupported
    /// source states, failed solves, exhausted work/refinement, cancellation and
    /// final conservation refusal roll back the complete call.
    ///
    /// A native geometry stop publishes earlier complete comparison intervals;
    /// its unfinished coarse/fine comparison is discarded. Local tolerances do
    /// not certify global solution error. The separate prefix damping allowance
    /// is an actual energy account, not a bound on all physical heat error.
    ///
    /// With an absolute-only damping budget, a nonterminal trial reserves credit
    /// using current centered kinetic plus surface energy. This scheduling
    /// heuristic is released when centered kinetic energy fits remaining credit;
    /// every actual prefix still satisfies the original damping inequality.
    pub fn advance_implicit_adaptive(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        implicit: ThermalJetImplicitControls,
        adaptive: ThermalJetAdaptiveControls,
        gate: &CancelGate,
    ) -> Result<ThermalJetAdaptiveAdvance, ThermalJetAdaptiveError> {
        self.advance_implicit_adaptive_with_method(
            duration_s,
            budget,
            implicit,
            adaptive,
            ThermalJetImplicitMethod::ImexEuler,
            gate,
        )
    }

    /// Select first-order BE or second-order Strang/midpoint with the same
    /// source, work, heat and whole-call transaction controls.
    ///
    /// The fine-path local estimate divides the coarse/fine difference by
    /// 2^p-1 (one for BE, three for midpoint), with square/cube-root scheduling.
    /// Midpoint has no BE kinetic damping, but unresolved stiff modes may reverse
    /// or oscillate: local velocity, caloric and physical-heat controls remain
    /// essential. A zero damping account alone never certifies heat accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_implicit_adaptive_with_method(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        implicit: ThermalJetImplicitControls,
        adaptive: ThermalJetAdaptiveControls,
        method: ThermalJetImplicitMethod,
        gate: &CancelGate,
    ) -> Result<ThermalJetAdaptiveAdvance, ThermalJetAdaptiveError> {
        poll(gate)?;
        implicit.validate()?;
        adaptive.validate()?;
        let b = budget.flow;
        if !duration_s.is_finite()
            || duration_s < 0.0
            || !b.maximum_step_s.is_finite()
            || b.maximum_step_s <= 0.0
            || !b.cfl.is_finite()
            || b.cfl <= 0.0
            || b.cfl > 0.25
            || b.maximum_substeps == 0
            || !budget.energy_tolerance_j.is_finite()
            || budget.energy_tolerance_j < 0.0
        {
            return Err(ThermalJetAdaptiveError::Invalid(
                "duration or coupled step budget",
            ));
        }
        let mut work = match self.implicit_workspace.take() {
            Some(work) => work,
            None => super::Workspace::new(self.accepted.len(), gate)?,
        };
        if work.adaptive.is_none() {
            match Workspace::new(self.accepted.len(), gate) {
                Ok(scratch) => work.adaptive = Some(scratch),
                Err(error) => {
                    self.implicit_workspace = Some(work);
                    return Err(error.into());
                }
            }
        }
        let mut scratch = work.adaptive.take().expect("adaptive workspace prepared");
        let result = self.advance_adaptive_inner(
            duration_s,
            budget,
            implicit,
            adaptive,
            method,
            &mut work,
            &mut scratch,
            gate,
        );
        work.adaptive = Some(scratch);
        self.implicit_workspace = Some(work);
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn charged_implicit_step(
        &mut self,
        dt: f64,
        base: Evaluation,
        budget: ThermalJetStepBudget,
        implicit: ThermalJetImplicitControls,
        method: ThermalJetImplicitMethod,
        heat: TrialHeat,
        work: &mut super::Workspace,
        totals: &mut AdvanceTotals,
        gate: &CancelGate,
    ) -> Result<Option<TrialStep>, ThermalJetImplicitError> {
        poll(gate)?;
        if totals.attempted >= budget.flow.maximum_substeps {
            return Err(NonlinearJetError::SubstepBudget.into());
        }
        totals.attempted += 1;
        let result = self.try_implicit_step(
            dt,
            base,
            budget.flow.cfl,
            implicit,
            method,
            heat,
            work,
            &mut totals.evaluations,
            gate,
        )?;
        if let Some(step) = &result {
            totals.nonlinear += step.solved.nonlinear;
            totals.linear += step.solved.linear;
            totals.maximum_velocity = totals.maximum_velocity.max(step.solved.velocity_residual);
            totals.maximum_energy = totals.maximum_energy.max(step.solved.energy_residual);
        }
        Ok(result)
    }

    /// Leaves the fine endpoint in stage and its two-half physical heat in step.
    #[allow(clippy::too_many_arguments)]
    fn compare_implicit_interval(
        &mut self,
        dt: f64,
        base: Evaluation,
        budget: ThermalJetStepBudget,
        implicit: ThermalJetImplicitControls,
        adaptive: ThermalJetAdaptiveControls,
        method: ThermalJetImplicitMethod,
        work: &mut super::Workspace,
        scratch: &mut Workspace,
        totals: &mut AdvanceTotals,
        gate: &CancelGate,
    ) -> Result<Option<Comparison>, ThermalJetAdaptiveError> {
        let Some(coarse) = self.charged_implicit_step(
            dt,
            base,
            budget,
            implicit,
            method,
            TrialHeat::Reset,
            work,
            totals,
            gate,
        )?
        else {
            return Ok(None);
        };
        scratch
            .coarse
            .save(&self.flow.stage, &self.stage, &self.workspace.profile, gate)?;
        let cell_heat = self.flow.spatial.trial_viscous_heat_j()?;
        for (i, &value) in cell_heat.iter().enumerate() {
            poll_at(gate, i)?;
            scratch.coarse_heat[i] = value;
        }
        let coarse_heat = coarse.heat;
        scratch.base.restore(self, gate)?;
        let base = self.evaluate_transport(false, gate)?;
        let half = positive(0.5 * dt)?;
        let Some(first) = self.charged_implicit_step(
            half,
            base,
            budget,
            implicit,
            method,
            TrialHeat::Reset,
            work,
            totals,
            gate,
        )?
        else {
            return Ok(None);
        };
        core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
        core::mem::swap(&mut self.candidate, &mut self.stage);
        let middle = self.evaluate_transport(false, gate)?;
        // The original kernel checks the endpoint transport bound too.
        if half > stable_step(middle.rate_per_s, budget.flow.cfl) {
            return Ok(None);
        }
        let Some(second) = self.charged_implicit_step(
            half,
            middle,
            budget,
            implicit,
            method,
            TrialHeat::Add,
            work,
            totals,
            gate,
        )?
        else {
            return Ok(None);
        };
        let heat = finite(first.heat + second.heat)?;
        Ok(Some(Comparison {
            local: scratch
                .indicators(self, coarse_heat, heat, adaptive, gate)?
                .for_method(method),
            heat,
            upwind: finite(first.upwind + second.upwind)?,
            damping: finite(first.damping + second.damping)?,
        }))
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_adaptive_inner(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        implicit: ThermalJetImplicitControls,
        adaptive: ThermalJetAdaptiveControls,
        method: ThermalJetImplicitMethod,
        work: &mut super::Workspace,
        scratch: &mut Workspace,
        gate: &CancelGate,
    ) -> Result<ThermalJetAdaptiveAdvance, ThermalJetAdaptiveError> {
        let target = self.flow.time_s + duration_s;
        if !target.is_finite() || (duration_s > 0.0 && target <= self.flow.time_s) {
            return Err(NonlinearJetError::NonRepresentable.into());
        }
        for i in 0..self.accepted.len() {
            poll_at(gate, i)?;
            self.flow.candidate.a[i] = self.flow.accepted.a[i];
            self.flow.candidate.q[i] = self.flow.accepted.q[i];
            self.candidate[i] = self.accepted[i];
        }
        self.flow.spatial.reset_viscous_heat(gate)?;
        let b = budget.flow;
        let mut base = self.evaluate_transport(false, gate)?;
        let mut totals = AdvanceTotals::default();
        let mut ledger = Ledger::default();
        let mut upwind = Sum::default();
        let (mut elapsed, mut cap) = (0.0, b.maximum_step_s);
        let (mut accepted, mut attempted, mut rejected, mut refinements) = (0, 0, 0, 0);
        let mut maximum_local = ThermalJetAdaptiveIndicators::default();
        let mut maximum_damping = 0.0_f64;
        let (mut minimum_interval, mut maximum_interval) = (f64::INFINITY, 0.0_f64);
        let mut outcome = SlenderJetOutcome::Completed;
        let mass_factor = positive(
            core::f64::consts::PI * self.flow.parameters.density_kg_m3 * self.flow.spacing_m,
        )?;
        while elapsed < duration_s {
            poll(gate)?;
            let remaining_time = duration_s - elapsed;
            let dt = positive(
                remaining_time
                    .min(cap)
                    .min(stable_step(base.rate_per_s, b.cfl)),
            )?;
            let floor = adaptive.minimum_interval_s.min(remaining_time);
            if dt < floor {
                return Err(ThermalJetAdaptiveError::Resolution {
                    interval_s: dt,
                    refinements,
                    local: None,
                    damping_ratio: 0.0,
                });
            }
            let next_elapsed = if dt == remaining_time {
                duration_s
            } else {
                elapsed + dt
            };
            if next_elapsed <= elapsed
                || self.flow.time_s + next_elapsed <= self.flow.time_s + elapsed
            {
                return Err(NonlinearJetError::NonRepresentable.into());
            }
            scratch.base.save(
                &self.flow.candidate,
                &self.candidate,
                &self.workspace.profile,
                gate,
            )?;
            let centered = if adaptive.relative_numerical_damping == 0.0 {
                centered_kinetic(&scratch.base.fields, mass_factor, gate)?
            } else {
                0.0
            };
            attempted += 1;
            let comparison = match self.compare_implicit_interval(
                dt,
                base,
                budget,
                implicit,
                adaptive,
                method,
                work,
                scratch,
                &mut totals,
                gate,
            ) {
                Err(ThermalJetAdaptiveError::Implicit(ThermalJetImplicitError::Thermal(
                    ThermalJetError::Flow(NonlinearJetError::InitialDomain(violation)),
                ))) => {
                    scratch.base.restore(self, gate)?;
                    base = self.evaluate_transport(false, gate)?;
                    outcome = SlenderJetOutcome::DomainLimit(violation);
                    break;
                }
                Err(error) => return Err(error),
                Ok(value) => value,
            };
            let mut local = None;
            let mut damping_ratio = 0.0;
            if let Some(comparison) = comparison {
                local = Some(comparison.local);
                let remaining_credit = ledger.remaining(adaptive)?;
                let r = adaptive.relative_numerical_damping;
                let cost = finite((1.0 - r) * comparison.damping)?;
                let allowance = finite(remaining_credit + r * comparison.heat)?;
                damping_ratio = ratio(cost, allowance);
                // Absolute-only scheduling reserves credit for unresolved
                // centered kinetic/surface inventory BEFORE accepting this trial.
                // Releasing the reserve once K fits credit avoids an undamped
                // surface baseline trapping a physically quiet tail.
                if r == 0.0 && next_elapsed < duration_s && centered > remaining_credit {
                    let pool = positive(centered + base.diagnostics.surface_energy_j)?;
                    let fraction = (finite(comparison.heat + comparison.damping)? / pool).min(1.0);
                    damping_ratio =
                        damping_ratio.max(ratio(cost, finite(remaining_credit * fraction)?));
                }
                let next_ledger = ledger.append(comparison.heat, comparison.damping);
                let permitted = damping_ratio <= 1.0 && next_ledger.allowed(adaptive)?;
                if comparison.local.scaled_maximum <= 1.0 && permitted {
                    self.flow.spatial.accept_viscous_heat_step(gate)?;
                    ledger = next_ledger;
                    upwind.add(comparison.upwind);
                    totals.accepted += 2;
                    accepted += 1;
                    maximum_local.merge(comparison.local);
                    maximum_damping = maximum_damping.max(damping_ratio);
                    minimum_interval = minimum_interval.min(dt);
                    maximum_interval = maximum_interval.max(dt);
                    core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
                    core::mem::swap(&mut self.candidate, &mut self.stage);
                    base = self.evaluate_transport(false, gate)?;
                    elapsed = next_elapsed;
                    refinements = 0;
                    let indicator = comparison.local.scaled_maximum.max(damping_ratio);
                    let growth = if indicator == 0.0 {
                        2.0
                    } else {
                        let root = match method {
                            ThermalJetImplicitMethod::ImexEuler => fs_math::det::sqrt(indicator),
                            ThermalJetImplicitMethod::StrangMidpoint => {
                                fs_math::det::pow(indicator, 1.0 / 3.0)
                            }
                        };
                        (0.9 / root).min(2.0)
                    };
                    cap = positive((dt * growth).min(b.maximum_step_s))?;
                    continue;
                }
            }
            rejected += 1;
            scratch.base.restore(self, gate)?;
            base = self.evaluate_transport(false, gate)?;
            let smaller = 0.5 * dt;
            if refinements >= adaptive.maximum_refinements || smaller < floor || smaller == 0.0 {
                return Err(ThermalJetAdaptiveError::Resolution {
                    interval_s: dt,
                    refinements: refinements + 1,
                    local,
                    damping_ratio,
                });
            }
            refinements += 1;
            cap = smaller;
        }
        totals.heat = finite(ledger.heat.to_f64())?;
        totals.damping = finite(ledger.damping.to_f64())?;
        totals.upwind = finite(upwind.total())?;
        let allowance = finite(
            adaptive.absolute_numerical_damping_j
                + adaptive.relative_numerical_damping * (totals.heat + totals.damping),
        )?;
        let evolution = self.publish_implicit_result(
            method,
            duration_s,
            elapsed,
            outcome,
            base,
            totals,
            budget.energy_tolerance_j,
            gate,
        )?;
        Ok(ThermalJetAdaptiveAdvance {
            evolution,
            accepted_intervals: accepted,
            attempted_intervals: attempted,
            rejected_intervals: rejected,
            maximum_accepted_local: maximum_local,
            maximum_accepted_damping_ratio: maximum_damping,
            numerical_damping_allowance_j: allowance,
            minimum_accepted_interval_s: if accepted == 0 { 0.0 } else { minimum_interval },
            maximum_accepted_interval_s: maximum_interval,
        })
    }
}

#[cfg(test)]
mod tests;
