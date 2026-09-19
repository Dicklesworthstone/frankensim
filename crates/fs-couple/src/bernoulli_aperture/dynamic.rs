//! Stateful moving-slit junction with caller-supplied mechanics and contact.
//!
//! This is a public host for the existing massive-reed midpoint island, not a
//! new stepper. The caller supplies every pressure/body-flow sample and an
//! fs-dcontact obstacle (including one constructed from a matdb receipt).
//! No instrument preset, attack envelope, inferred lay stiffness or default
//! damping enters this path. Observer propagation and the external bore/plate
//! dynamics remain the caller's responsibility. Energy records describe this
//! local junction only; they are not physical-validation certificates.
//!
//! Storage and nonlinear trials currently allocate through the existing
//! owners. A bounded iteration count is not an allocation-free or hard-real-
//! time claim. Cancellation is observed between complete mechanical steps.

use super::BernoulliAperture;
use crate::acoustic_realize::AcousticRealizeError;
use crate::reed_bore::{FastSolveStats, ReedSolverMode, reed_pressure_face, reed_structural};
use crate::unilateral_contact::slit_contact_discrete_coefficients;
use fs_dcontact::{ContactStorage, Obstacle};
use fs_exec::CancelGate;
use fs_phs::Storage;
use fs_scenario::BeatingReed;

/// Explicit SI parameters for a moving slit and its characteristic load.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DynamicApertureSpec {
    /// Rest geometry and static closing pressure; together with stiffness these
    /// define the effective pressure area A = k H / P_c.
    pub aperture: BernoulliAperture,
    /// Moving effective mass [kg], strictly positive.
    pub mass_kg: f64,
    /// Restoring stiffness [N/m], strictly positive (no inferred fallback).
    pub stiffness_n_m: f64,
    /// Viscous damping ratio, nonnegative and explicitly supplied.
    pub damping_ratio: f64,
    /// Fluid density [kg/m^3].
    pub density_kg_m3: f64,
    /// Positive characteristic pressure/volume-flow impedance [Pa s/m^3].
    pub impedance_pa_s_m3: f64,
    /// Fixed mechanical step [s].
    pub time_step_s: f64,
    /// Total accepted-step budget, including steps preceding a resume.
    pub max_steps: u64,
}

/// Mechanical state; a negative opening represents penetration of the lay.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ApertureState {
    /// Opening coordinate [m].
    pub opening_m: f64,
    /// Opening velocity [m/s]; positive widens the slit.
    pub opening_velocity_m_s: f64,
}

/// Inputs held during one midpoint solve. Reverse pressure/flow is supported.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ApertureDrive {
    /// Upstream pressure [Pa], not an acoustic output or a pitch target.
    pub upstream_pressure_pa: f64,
    /// Incoming characteristic pressure wave [Pa].
    pub incoming_pressure_pa: f64,
    /// Additional body-supplied volume flow into the bore [m^3/s].
    pub body_flow_m3_s: f64,
}

/// One accepted physical step. Fluxes and work refer to its midpoint;
/// `state` and `stored_energy_j` refer to the end of the step.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ApertureFrame {
    /// One-based accepted-step ordinal (unchanged by block partitioning).
    pub step: u64,
    /// End time [s], derived from ordinal rather than accumulated additions.
    pub time_s: f64,
    /// End-of-step mechanical state.
    pub state: ApertureState,
    /// Nonnegative midpoint opening actually used in the Bernoulli solve [m].
    pub midpoint_opening_m: f64,
    /// Outgoing characteristic pressure [Pa].
    pub outgoing_pressure_pa: f64,
    /// Sum of incoming and outgoing characteristic pressure [Pa].
    pub bore_pressure_pa: f64,
    /// Signed Bernoulli jet flow [m^3/s].
    pub jet_flow_m3_s: f64,
    /// Signed flow due to moving pressure face, -A v_mid [m^3/s].
    pub swept_flow_m3_s: f64,
    /// Net characteristic flow (p_plus - p_minus) / Z [m^3/s].
    pub bore_flow_m3_s: f64,
    /// Jet + swept + supplied body flow - characteristic flow [m^3/s].
    pub flow_residual_m3_s: f64,
    /// Mechanical plus contact stored energy [J].
    pub stored_energy_j: f64,
    /// Change in stored energy over this step [J].
    pub storage_change_j: f64,
    /// Jet, structural and nonadhesive contact dissipation over this step [J].
    pub dissipated_energy_j: f64,
    /// Differential-pressure work on this junction [J]. This excludes any
    /// storage/loss in a downstream bore or an externally driven plate.
    pub pressure_work_j: f64,
}

impl ApertureFrame {
    /// Local energy-balance residual [J]; numerical diagnostic, not a proof.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j - self.pressure_work_j
    }
}

/// Why an otherwise valid block stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApertureTerminal {
    /// Every requested input was stepped.
    Complete,
    /// Cancellation was observed before the next step.
    Cancelled,
    /// The total accepted-step budget was reached.
    BudgetExhausted,
}

/// Only `completed` output slots were written. The suffix remains untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ApertureProgress {
    /// Accepted steps in this call, not over the model's whole lifetime.
    pub completed: usize,
    /// Completion/cancellation/budget state.
    pub terminal: ApertureTerminal,
}

/// A persistent moving slit whose contact law is never chosen implicitly.
pub struct DynamicAperture {
    spec: DynamicApertureSpec,
    reed: BeatingReed,
    lay: Obstacle,
    contact_storage: ContactStorage,
    state: ApertureState,
    accepted_steps: u64,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

// Ordinals remain exactly representable in the floating-point time projection.
const MAX_EXACT_STEPS: u64 = 1_u64 << 53;

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _x: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _x: &[f64], out: &mut [f64]) { out.fill(0.0); }
}

impl DynamicAperture {
    /// Admit explicit mechanics, initial state and a unit-opening obstacle.
    /// The obstacle's gap, weight, K, alpha, chi and provenance are preserved.
    /// For a no-contact comparison, supply an explicit zero-stiffness obstacle.
    ///
    /// # Errors
    /// Invalid or nonfinite mechanics, contact shape/law, time/budget or energy.
    pub fn new(
        spec: DynamicApertureSpec,
        state: ApertureState,
        lay: Obstacle,
    ) -> Result<Self, AcousticRealizeError> {
        if ![
            spec.aperture.rest_opening_m, spec.aperture.width_m,
            spec.aperture.closing_pressure_pa, spec.mass_kg, spec.stiffness_n_m,
            spec.damping_ratio, spec.density_kg_m3, spec.impedance_pa_s_m3,
            spec.time_step_s, state.opening_m, state.opening_velocity_m_s,
        ].iter().all(|x| x.is_finite())
            || spec.aperture.rest_opening_m <= 0.0 || spec.aperture.width_m <= 0.0
            || spec.aperture.closing_pressure_pa <= 0.0 || spec.mass_kg <= 0.0
            || spec.stiffness_n_m <= 0.0 || spec.damping_ratio < 0.0
            || spec.density_kg_m3 <= 0.0 || spec.impedance_pa_s_m3 <= 0.0
            || spec.time_step_s <= 0.0 || spec.max_steps == 0
            || spec.max_steps > MAX_EXACT_STEPS
            || !(spec.time_step_s * spec.max_steps as f64).is_finite()
        {
            return Err(invalid("dynamic aperture requires finite explicit mechanics and a positive step budget"));
        }
        // Check raw-parts escape values before the generic storage can index
        // them. Unit opening is a coordinate contract, not a material inference.
        if lay.n_points() != 1 || lay.collocation() != [-1.0]
            || lay.gaps().len() != 1 || lay.weights().len() != 1
            || !lay.gaps()[0].is_finite() || !lay.weights()[0].is_finite()
            || lay.weights()[0] < 0.0 || !lay.stiffness().is_finite()
            || lay.stiffness() < 0.0 || !lay.alpha().is_finite() || lay.alpha() < 1.0
            || !lay.internal_loss().is_finite() || lay.internal_loss() < 0.0
            || lay.provenance().trim().is_empty()
        {
            return Err(invalid("dynamic aperture requires a finite, provenance-labelled unit-opening contact law"));
        }
        slit_contact_discrete_coefficients(&lay, state.opening_m, state.opening_m)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let contact_storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![lay.clone()])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let reed = BeatingReed {
            rest_opening_m: spec.aperture.rest_opening_m,
            width_m: spec.aperture.width_m,
            closing_pressure_pa: spec.aperture.closing_pressure_pa,
            mass_kg: spec.mass_kg,
            stiffness_n_m: spec.stiffness_n_m,
            damping_ratio: spec.damping_ratio,
            // These fields do not enter the mechanical island. Every actual
            // input is an explicit ApertureDrive; no envelope is synthesized.
            blowing_pressure_pa: 0.0,
            attack_s: 0.0,
        };
        let model = Self { spec, reed, lay, contact_storage, state, accepted_steps: 0 };
        let (stiffness, damping) = reed_structural(reed);
        if !reed_pressure_face(reed).is_finite() || !stiffness.is_finite()
            || !damping.is_finite() || !model.energy_at(state).is_finite()
        {
            return Err(invalid("dynamic aperture derived mechanics or initial energy overflowed"));
        }
        Ok(model)
    }

    /// Current accepted mechanical state, including nonzero vibration on resume.
    #[must_use]
    pub const fn state(&self) -> ApertureState { self.state }

    /// Total accepted steps; failed trials do not advance this clock.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 { self.accepted_steps }

    /// Immutable explicit parameters, including the current total step budget.
    #[must_use]
    pub const fn spec(&self) -> &DynamicApertureSpec { &self.spec }

    /// Contact law supplied at admission; no synthetic source is substituted.
    #[must_use]
    pub fn contact_law(&self) -> &Obstacle { &self.lay }

    fn energy_at(&self, state: ApertureState) -> f64 {
        let displacement = state.opening_m - self.reed.rest_opening_m;
        0.5 * self.reed.mass_kg * state.opening_velocity_m_s * state.opening_velocity_m_s
            + 0.5 * self.reed.stiffness_n_m * displacement * displacement
            + self.contact_storage.hamiltonian(&[state.opening_m, 0.0])
    }

    /// Increase the total budget without resetting physical state or time.
    ///
    /// # Errors
    /// A non-increase, unrepresentable ordinal or overflowing time horizon.
    pub fn extend_step_budget(&mut self, new_total: u64) -> Result<(), AcousticRealizeError> {
        if new_total <= self.spec.max_steps || new_total > MAX_EXACT_STEPS
            || !(new_total as f64 * self.spec.time_step_s).is_finite()
        {
            return Err(invalid("dynamic aperture budget extension must increase a finite exact-step horizon"));
        }
        self.spec.max_steps = new_total;
        Ok(())
    }

    /// Advance one input transactionally through the existing midpoint solver.
    ///
    /// # Errors
    /// Budget, invalid input, numerical solve or nonfinite observation. Every
    /// failure leaves the old state and accepted-step count unchanged.
    pub fn step(&mut self, drive: ApertureDrive) -> Result<ApertureFrame, AcousticRealizeError> {
        if self.accepted_steps >= self.spec.max_steps {
            return Err(AcousticRealizeError::Reed { what: "dynamic aperture step budget exhausted" });
        }
        if ![drive.upstream_pressure_pa, drive.incoming_pressure_pa, drive.body_flow_m3_s]
            .iter().all(|x| x.is_finite())
        {
            return Err(invalid("dynamic aperture pressure and flow inputs must be finite"));
        }
        let old = self.state;
        let dt = self.spec.time_step_s;
        let (outgoing, opening, velocity) = crate::reed_bore::step_massive_reed(
            self.reed, self.spec.density_kg_m3, self.spec.impedance_pa_s_m3,
            drive.incoming_pressure_pa, drive.upstream_pressure_pa,
            old.opening_m, old.opening_velocity_m_s, dt, drive.body_flow_m3_s,
            Some(&self.lay), ReedSolverMode::Strict, &mut FastSolveStats::default(),
        )?;
        let state = ApertureState { opening_m: opening, opening_velocity_m_s: velocity };
        let midpoint_opening_m = f64::midpoint(old.opening_m, opening).max(0.0);
        let vm = f64::midpoint(old.opening_velocity_m_s, velocity);
        let bore = outgoing + drive.incoming_pressure_pa;
        let dp = drive.upstream_pressure_pa - bore;
        let jet = fs_phs::bernoulli_volume_flow(
            self.reed.width_m, midpoint_opening_m, dp, self.spec.density_kg_m3,
        );
        let swept = -reed_pressure_face(self.reed) * vm;
        let wave = (outgoing - drive.incoming_pressure_pa) / self.spec.impedance_pa_s_m3;
        let (elastic, contact_damping) = slit_contact_discrete_coefficients(
            &self.lay, old.opening_m, opening,
        ).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let force = (elastic - contact_damping * vm).max(0.0);
        let (_, damping) = reed_structural(self.reed);
        let energy = self.energy_at(state);
        let frame = ApertureFrame {
            step: self.accepted_steps + 1,
            time_s: (self.accepted_steps + 1) as f64 * dt,
            state,
            midpoint_opening_m,
            outgoing_pressure_pa: outgoing,
            bore_pressure_pa: bore,
            jet_flow_m3_s: jet,
            swept_flow_m3_s: swept,
            bore_flow_m3_s: wave,
            flow_residual_m3_s: jet + swept + drive.body_flow_m3_s - wave,
            stored_energy_j: energy,
            storage_change_j: energy - self.energy_at(old),
            dissipated_energy_j: dt * (dp * jet + damping * vm * vm + (elastic - force) * vm),
            pressure_work_j: dt * dp * (wave - drive.body_flow_m3_s),
        };
        if ![frame.time_s, frame.bore_pressure_pa, frame.jet_flow_m3_s,
            frame.swept_flow_m3_s, frame.bore_flow_m3_s, frame.flow_residual_m3_s,
            frame.stored_energy_j, frame.storage_change_j, frame.dissipated_energy_j,
            frame.pressure_work_j, frame.balance_residual_j()].iter().all(|v| v.is_finite())
        {
            return Err(AcousticRealizeError::Reed { what: "dynamic aperture observation left the finite set" });
        }
        // No fallible work follows publication of the accepted physical step.
        self.state = state;
        self.accepted_steps = frame.step;
        Ok(frame)
    }

    /// Step a caller-sized block with sample-boundary cancellation and resume.
    /// The output suffix is never touched after cancellation or exhaustion.
    /// Resume with the unused input/output suffix and a fresh cancellation gate;
    /// after budget exhaustion, explicitly extend the total step budget first.
    ///
    /// # Errors
    /// Shape mismatch (no progress), or a per-step refusal. On a numerical/input
    /// refusal, previously completed steps remain accepted; inspect
    /// `accepted_steps()` to recover that prefix, then retry the refused input.
    pub fn advance_block(
        &mut self,
        inputs: &[ApertureDrive],
        out: &mut [ApertureFrame],
        gate: &CancelGate,
    ) -> Result<ApertureProgress, AcousticRealizeError> {
        if inputs.len() != out.len() {
            return Err(invalid("dynamic aperture input/output block lengths must match"));
        }
        for (completed, (drive, slot)) in inputs.iter().zip(out.iter_mut()).enumerate() {
            if gate.is_requested() {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::Cancelled });
            }
            if self.accepted_steps >= self.spec.max_steps {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::BudgetExhausted });
            }
            *slot = self.step(*drive)?;
        }
        Ok(ApertureProgress { completed: inputs.len(), terminal: ApertureTerminal::Complete })
    }
}
