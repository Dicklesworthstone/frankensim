//! Reciprocal acoustic reaction on additional faces of the original moving body.
//!
//! Each independent, fixed-geometry DG domain supplies its exact held-flow mean
//! law p=H+Z*Q. Signed exposed area a maps the ONE wall velocity to Q=a*v and
//! pressure back to force -a*p. Eliminating the linear loads gives one affine
//! force, so the original two-gas impulse solve owns the only momentum search.
//!
//! These are additional acoustic domains coupled by mechanical faces, not the
//! two thermodynamic chamber inventories. No chamber volume is counted again,
//! no valve mass flow is inferred from acoustic volume flow, and no density
//! equality is required between independently sourced acoustic environments.

use super::{DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve, DoublePistonStep};
use crate::acoustic_realize::AcousticRealizeError;
use crate::cavity_phs::piston::PistonError;
use crate::dg_port::{DgPortLoad, DgTerminalFrame, sum_terminal_work};
use fs_duct::dg::DgStepOutput;

/// One distinct acoustic domain acted on by a face of the translating body.
/// Positive area injects flow on primary expansion; negative area withdraws it.
/// The original duct inlet remains the acoustic boundary; area is the physical
/// swept-volume coefficient of the attached face, not a replacement duct area.
#[derive(Clone, Copy, Debug)]
pub struct AcousticMotionFace<'a> {
    pub load: &'a DgPortLoad,
    pub signed_area_m2: f64,
}

/// Explicit work and admission limits, independent of the original gas solver.
/// No extra nonlinear iteration policy is needed: acoustic reaction contributes
/// its exact derivative to the existing shared-wall impulse equation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticMotionSolve {
    pub maximum_faces: usize,
    pub pressure_tolerance_pa: f64,
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
}

/// A refusal leaves every supplied acoustic and mechanical owner unchanged.
#[derive(Debug)]
pub enum AcousticMotionError<E> {
    Invalid(&'static str),
    Acoustic(AcousticRealizeError),
    Mechanical(E),
    Cancelled,
}
impl<E> From<AcousticRealizeError> for AcousticMotionError<E> {
    fn from(error: AcousticRealizeError) -> Self {
        match error {
            AcousticRealizeError::Cancelled => Self::Cancelled,
            other => Self::Acoustic(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for AcousticMotionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "acoustic/shared motion: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for AcousticMotionError<E> {}

/// Actual original load exchange, ordered exactly as the authored faces.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticMotionFaceStep {
    pub signed_area_m2: f64,
    pub inlet_flow_m3_s: f64,
    pub mean_pressure_pa: f64,
    pub inlet_work_j: f64,
    pub port_work_residual_j: f64,
    /// Includes original branch and finite terminal-channel storage and losses.
    pub duct: DgStepOutput,
    /// Exterior work is signed exchange, not unmodeled terminal dissipation.
    pub terminals: Vec<DgTerminalFrame>,
}

/// The actual mechanical/acoustic exchange, without duplicating gas storage.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticMotionStep {
    /// Sum of pressure forces on the original shared motion [N].
    pub force_n: f64,
    pub midpoint_velocity_m_s: f64,
    /// Native extra-force work minus the authored external-force work [J].
    pub mechanical_work_j: f64,
    /// Affine force minus the independently observed native pressure force [N].
    pub force_residual_n: f64,
    pub acoustic_storage_change_j: f64,
    pub acoustic_dissipated_j: f64,
    pub terminal_work_j: f64,
    /// Sum of original storage + dissipation + exterior work - inlet work [J].
    pub acoustic_energy_residual_j: f64,
    /// ACTUAL mechanical work plus all original DG inlet work [J].
    pub transfer_residual_j: f64,
    pub faces: Vec<AcousticMotionFaceStep>,
}

pub(super) struct PreparedAcousticMotion<'a> {
    faces: Vec<AcousticMotionFace<'a>>,
    responses: Vec<(f64, f64)>,
    duration_s: f64,
    controls: AcousticMotionSolve,
    pub(super) history_force_n: f64,
    pub(super) resistance_n_s_m: f64,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
pub(super) fn check<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticMotionError<E>> {
    if cancelled() {
        Err(AcousticMotionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Inspect every accepted history once. No propagator or physical state changes.
pub(super) fn prepare<'a, E>(
    faces: &[AcousticMotionFace<'a>],
    duration_s: f64,
    accepted_steps: u64,
    controls: AcousticMotionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedAcousticMotion<'a>, AcousticMotionError<E>> {
    check(cancelled)?;
    if faces.is_empty()
        || faces.len() > controls.maximum_faces
        || !positive(duration_s)
        || !positive(controls.pressure_tolerance_pa)
        || !positive(controls.velocity_tolerance_m_s)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(AcousticMotionError::Invalid(
            "acoustic motion needs a bounded nonempty face group, positive clock and physical tolerances",
        ));
    }
    let mut responses = Vec::with_capacity(faces.len());
    let mut history_force_n = 0.0;
    let mut resistance_n_s_m = 0.0;
    for (index, face) in faces.iter().enumerate() {
        check(cancelled)?;
        if !face.signed_area_m2.is_finite()
            || face.load.has_driven_mouth()
            || face.load.propagator().dt().to_bits() != duration_s.to_bits()
            || face.load.state().step_count() != accepted_steps
            || faces[..index]
                .iter()
                .any(|other| std::ptr::eq(other.load, face.load))
        {
            return Err(AcousticMotionError::Invalid(
                "each face requires a distinct original load on the mechanical clock, a finite area and no driven exterior mouth",
            ));
        }
        let relation = face.load.port_response()?;
        let h = relation.history_pressure_pa;
        let z = relation.impedance_pa_s_m3;
        if !h.is_finite() || !positive(z) {
            return Err(AcousticMotionError::Invalid(
                "original acoustic mean relation requires finite history and positive impedance",
            ));
        }
        // A zero-area face retains its own free acoustic continuation, while
        // contributing exactly zero force and avoiding signed-zero changes.
        if face.signed_area_m2 != 0.0 {
            history_force_n -= face.signed_area_m2 * h;
            resistance_n_s_m += (face.signed_area_m2 * z) * face.signed_area_m2;
        }
        if !history_force_n.is_finite() || !resistance_n_s_m.is_finite() {
            return Err(AcousticMotionError::Invalid(
                "aggregate acoustic pressure force or resistance is unrepresentable",
            ));
        }
        responses.push((h, z));
    }
    check(cancelled)?;
    Ok(PreparedAcousticMotion {
        faces: faces.to_vec(),
        responses,
        duration_s,
        controls,
        history_force_n,
        resistance_n_s_m,
    })
}

/// March each original complete load once, using the accepted native motion.
/// All accepted owners remain borrowed immutably; partial private proposals are
/// discarded on any late work failure or cancellation.
pub(super) fn finish<E>(
    prepared: PreparedAcousticMotion<'_>,
    mechanical: &DoublePistonStep,
    external_force_n: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<DgPortLoad>, AcousticMotionStep), AcousticMotionError<E>> {
    check(cancelled)?;
    let dt = prepared.duration_s;
    let controls = prepared.controls;
    let displacement = mechanical.displacement_m;
    let velocity = displacement / dt;
    let endpoint_mean = 0.5 * mechanical.velocity_m_s[0] + 0.5 * mechanical.velocity_m_s[1];
    let force = prepared.history_force_n - prepared.resistance_n_s_m * velocity;
    let external_work = external_force_n * displacement;
    let mechanical_work = mechanical.external_force_work_j - external_work;
    let predicted_work = force * displacement;
    if ![
        displacement,
        velocity,
        endpoint_mean,
        force,
        external_work,
        mechanical_work,
        predicted_work,
    ]
    .iter()
    .all(|value| value.is_finite())
        || (velocity - endpoint_mean).abs() > controls.velocity_tolerance_m_s
        || (mechanical_work - predicted_work).abs() > controls.energy_tolerance_j
    {
        return Err(AcousticMotionError::Invalid(
            "actual shared-wall velocity or native force work disagrees with its acoustic load",
        ));
    }
    let mut loads = Vec::with_capacity(prepared.faces.len());
    let mut reports = Vec::with_capacity(prepared.faces.len());
    let mut observed_force = 0.0;
    let mut force_allowance = 0.0;
    let mut inlet_work = 0.0;
    let mut storage_change = 0.0;
    let mut dissipated = 0.0;
    let mut terminal_work = 0.0;
    for (face, &(history, impedance)) in prepared.faces.iter().zip(&prepared.responses) {
        check(cancelled)?;
        let flow = face.signed_area_m2 * velocity;
        let predicted_pressure = history + impedance * flow;
        if !flow.is_finite() || !predicted_pressure.is_finite() {
            return Err(AcousticMotionError::Invalid(
                "actual shared-wall flow or acoustic pressure is unrepresentable",
            ));
        }
        let mut load = face.load.clone();
        let id = load.inlet();
        let trial = load.preview_step(flow)?;
        let pressure = trial
            .output()
            .port_outputs
            .iter()
            .find(|port| port.id == id)
            .ok_or(AcousticMotionError::Invalid(
                "native acoustic trial lost its inlet",
            ))?
            .conjugate_rate;
        let work = trial.inlet_work_j();
        let defect = work - pressure * (flow * dt);
        let delta = trial.output().ledger.stored_energy_j - trial.previous_stored_energy_j();
        let loss = trial.output().ledger.cumulative_dissipated_j
            - trial.previous_cumulative_dissipated_j();
        let terminals = trial.terminal_frames().to_vec();
        let exterior = sum_terminal_work(&terminals);
        if ![pressure, work, defect, delta, loss, exterior]
            .iter()
            .all(|value| value.is_finite())
            || (pressure - predicted_pressure).abs() > controls.pressure_tolerance_pa
            || defect.abs() > controls.energy_tolerance_j
            || ((delta + loss) + exterior - work).abs() > controls.energy_tolerance_j
        {
            return Err(AcousticMotionError::Invalid(
                "native acoustic pressure, inlet work or complete load energy failed admission",
            ));
        }
        observed_force -= face.signed_area_m2 * pressure;
        force_allowance += face.signed_area_m2.abs() * controls.pressure_tolerance_pa;
        inlet_work += work;
        storage_change += delta;
        dissipated += loss;
        terminal_work += exterior;
        check(cancelled)?;
        // Publication is only into this private clone. No accepted load has
        // been mutated if a subsequent face or the final combined audit fails.
        let duct = trial.accept();
        reports.push(AcousticMotionFaceStep {
            signed_area_m2: face.signed_area_m2,
            inlet_flow_m3_s: flow,
            mean_pressure_pa: pressure,
            inlet_work_j: work,
            port_work_residual_j: defect,
            duct,
            terminals,
        });
        loads.push(load);
    }
    let force_residual = force - observed_force;
    let acoustic_residual = ((storage_change + dissipated) + terminal_work) - inlet_work;
    let transfer = mechanical_work + inlet_work;
    if ![
        observed_force,
        force_allowance,
        inlet_work,
        storage_change,
        dissipated,
        terminal_work,
        force_residual,
        acoustic_residual,
        transfer,
    ]
    .iter()
    .all(|value| value.is_finite())
        || force_residual.abs() > force_allowance
        || acoustic_residual.abs() > controls.energy_tolerance_j
        || transfer.abs() > controls.energy_tolerance_j
    {
        return Err(AcousticMotionError::Invalid(
            "actual pressure reaction, acoustic storage or reciprocal work failed admission",
        ));
    }
    check(cancelled)?;
    Ok((
        loads,
        AcousticMotionStep {
            force_n: force,
            midpoint_velocity_m_s: velocity,
            mechanical_work_j: mechanical_work,
            force_residual_n: force_residual,
            acoustic_storage_change_j: storage_change,
            acoustic_dissipated_j: dissipated,
            terminal_work_j: terminal_work,
            acoustic_energy_residual_j: acoustic_residual,
            transfer_residual_j: transfer,
            faces: reports,
        },
    ))
}

/// One complete original gas/wall candidate and all original acoustic endpoints.
/// Native gas, branch, material-wall, finite-channel and rational-load histories
/// stay with their original owners; consume and publish all endpoints together.
#[derive(Clone, Debug)]
pub struct AcousticPistonCandidate {
    piston: DoubleActingGasPiston,
    loads: Vec<DgPortLoad>,
    step: DoublePistonStep,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}
impl AcousticPistonCandidate {
    pub const fn piston(&self) -> &DoubleActingGasPiston {
        &self.piston
    }
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    pub const fn piston_step(&self) -> &DoublePistonStep {
        &self.step
    }
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    pub fn into_parts(
        self,
    ) -> (
        DoubleActingGasPiston,
        Vec<DgPortLoad>,
        DoublePistonStep,
        AcousticMotionStep,
        f64,
    ) {
        (
            self.piston,
            self.loads,
            self.step,
            self.motion,
            self.total_energy_residual_j,
        )
    }
}

impl DoubleActingGasPiston {
    /// Couple the original two independent finite gases and their ONE moving
    /// body to additional acoustic faces. Pressure is the original DG acoustic
    /// perturbation; chamber absolute-pressure work remains with this piston.
    /// Any uncompensated mean pressure on the additional faces belongs in the
    /// authored external force or a balanced reference configuration.
    ///
    /// Every load uses the same interval and accepted step count. Different
    /// acoustic gases, geometry, branches, wall memories and terminal loads are
    /// retained. A repeated reference to one load is refused instead of owning
    /// one acoustic history twice. Zero area preserves mechanical decoupling.
    ///
    /// The exact affine mean reaction enters the existing impulse tangent; each
    /// original DG owner is stepped once after momentum convergence. The held
    /// acoustic-flow approximation uses the native midpoint displacement and
    /// requires time refinement for moving-boundary accuracy. Fixed acoustic
    /// geometry, linear perturbation acoustics and additional distinct domains
    /// remain explicit; there is no valve jet or changing-chamber acoustic model.
    pub fn acoustic_candidate(
        &self,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        faces: &[AcousticMotionFace<'_>],
        motion: AcousticMotionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticPistonCandidate, AcousticMotionError<PistonError>> {
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            motion,
            &mut cancelled,
        )?;
        let force = if prepared.history_force_n == 0.0 {
            input.external_force_n
        } else {
            input.external_force_n + prepared.history_force_n
        };
        let mut next = self.clone();
        let step = next
            .step_with_midpoint_resistance(
                DoublePistonInput {
                    external_force_n: force,
                    ..input
                },
                controls,
                prepared.resistance_n_s_m,
                &mut cancelled,
            )
            .map_err(|error| match error {
                PistonError::Cancelled => AcousticMotionError::Cancelled,
                other => AcousticMotionError::Mechanical(other),
            })?;
        let (loads, exchange) = finish(prepared, &step, input.external_force_n, &mut cancelled)?;
        let total = step.total_energy_residual_j
            + exchange.acoustic_energy_residual_j
            + exchange.transfer_residual_j;
        if !total.is_finite() || total.abs() > motion.energy_tolerance_j {
            return Err(AcousticMotionError::Invalid(
                "gas, shared kinetic and complete acoustic energy do not close",
            ));
        }
        check(&mut cancelled)?;
        Ok(AcousticPistonCandidate {
            piston: next,
            loads,
            step,
            motion: exchange,
            total_energy_residual_j: total,
        })
    }
}

#[cfg(test)]
pub(super) mod tests;
