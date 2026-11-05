//! Resolved capacitor geometry with finite native circuit charge evolution.
//!
//! For `H(q,s) = q²/(2 C(s))`, put `a(s) = 1/C(s)`. The symmetric
//! coordinate-increment discrete gradient satisfies the exact identity
//! `delta H = (q0+q1)(a0+a1) delta q/4
//!            + (q0²+q1²)(a1-a0)/4`.
//! Thus the original circuit advances with `1/C_eff = (a0+a1)/2`, while
//! the original guided body receives `F = -(q0²+q1²) delta a/(4 delta s)`.
//! The same displacement closes both subproblems. No independent RC, rigid-body
//! or electrostatic integrator is introduced, and resistor heat is never
//! inferred from an energy residual.

use super::{
    CapacitanceSample, ElectricalControl, Error, MotionControls, PlanarCapacitorMotion,
    PlanarCapacitorPath, Result, ScalarError, check, feedback, finite, positive,
};
use crate::material_circuit::{
    MaterialCircuit, MaterialCircuitSolve, ResistorHeat, contacts::ContactHeat,
};
use fs_ga::{Se3, Vec3, Wrench};
use fs_math::dd::Dd;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, PrescribedBaseImpulseResponse, apply_prescribed_base_impulses,
};
use fs_mbd_material::ResolvedBody;
use fs_phs::{StepRecord, circuit::Branch};

/// Joint electrical, material and moving-field history. The complete original
/// circuit remains attached, including unrelated capacitors, mutual inductance,
/// resistor sources and electrical contacts. Only the named capacitor changes
/// with the resolved geometry; its charge is a native circuit coordinate.
#[derive(Clone)]
pub struct PlanarCapacitorCircuitMotion {
    motion: PlanarCapacitorMotion,
    circuit: MaterialCircuit,
    branch: usize,
    charge_slot: usize,
    electrical_energy_j: f64,
}

/// External sources held over one common physical interval.
#[derive(Clone, Copy, Debug)]
pub struct CircuitMotionInput<'a> {
    /// Native circuit port order. These inputs carry source work, not heat.
    pub electrical_ports: &'a [f64],
    /// Force along the positive electrode's translation axis [N].
    pub external_force_n: f64,
    pub duration_s: f64,
}

/// Native responses and independent physical endpoint energy accounting.
pub struct CircuitMotionStep {
    /// Native circuit subproblem. Its `delta_h` uses C_eff, not the changing
    /// physical geometry; `electrical_energy_j` is the physical storage pair.
    pub circuit: StepRecord,
    pub mechanics: PrescribedBaseImpulseResponse,
    /// Complete original electrical storage at each actual field geometry [J].
    /// This already includes `field_energy_j`; do not count that subset twice.
    pub electrical_energy_j: [f64; 2],
    pub field_energy_j: [f64; 2],
    pub electrical_supply_j: f64,
    pub external_mechanical_work_j: f64,
    pub field_work_on_body_j: f64,
    pub discrete_field_force_n: f64,
    pub displacement_m: f64,
    pub kinematic_residual_m: f64,
    /// Original material bulk loss at the same circuit midpoint efforts.
    pub resistors: Vec<ResistorHeat>,
    /// Original series-interface loss, separate from bulk resistor heat.
    pub contacts: Vec<ContactHeat>,
    pub heat_j: f64,
    pub electrical_energy_residual_j: f64,
    pub mechanical_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub field_solves: usize,
}

/// One indivisible field, body and circuit proposal; no accepted owner mutates.
pub struct CircuitMotionCandidate {
    state: PlanarCapacitorCircuitMotion,
    step: CircuitMotionStep,
}
impl CircuitMotionCandidate {
    pub const fn state(&self) -> &PlanarCapacitorCircuitMotion {
        &self.state
    }
    pub const fn step(&self) -> &CircuitMotionStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarCapacitorCircuitMotion, CircuitMotionStep) {
        (self.state, self.step)
    }
}

impl PlanarCapacitorCircuitMotion {
    /// Construct a new joint state. `branch` is the original graph index of
    /// the two-electrode field. Its first node corresponds to the path's moving
    /// positive terminal; its second corresponds to the fixed negative terminal.
    /// Its numeric C placeholder is replaced by the actual initial field solve,
    /// preserving supplied charge and every other source/storage binding.
    /// This initializes an assembly; it is not a live parameter-change API.
    ///
    /// The active branch must have no competing material reduction, nonlinear
    /// charge law or capacitance block. Other branches retain their laws. Ideal
    /// constraints must admit the original charge at the initial, intermediate
    /// and final capacitances; incompatible clamps refuse without projection.
    /// Use `PlanarCapacitorMotion` for directly held ideal electrode voltage.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        path: PlanarCapacitorPath,
        body: ResolvedBody,
        circuit: MaterialCircuit,
        branch: usize,
        coordinate_m: f64,
        velocity_m_s: f64,
        initial_inputs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if !circuit
            .graph()
            .branches
            .get(branch)
            .is_some_and(|entry| matches!(entry.2, Branch::Capacitor { .. }))
        {
            return Err(Error::Invalid(
                "resolved moving field needs an original capacitor branch",
            ));
        }
        let charge_slot = circuit.graph().branches[..branch]
            .iter()
            .filter(|entry| matches!(entry.2, Branch::Capacitor { .. }))
            .count();
        let charge_c = circuit.charges_c()[charge_slot];
        let motion = PlanarCapacitorMotion::new(
            path,
            body,
            ElectricalControl::IsolatedCharge { charge_c },
            coordinate_m,
            velocity_m_s,
            cancelled,
        )?;
        let (circuit, electrical_energy_j) = circuit.with_resolved_capacitance(
            branch,
            motion.capacitance().capacitance_f(),
            initial_inputs,
            cancelled,
        )?;
        check(cancelled)?;
        Ok(Self {
            motion,
            circuit,
            branch,
            charge_slot,
            electrical_energy_j,
        })
    }

    pub fn path(&self) -> &PlanarCapacitorPath {
        self.motion.path()
    }
    pub fn body(&self) -> &ResolvedBody {
        self.motion.body()
    }
    pub fn mechanical_model(&self) -> &ArticulatedModel {
        self.motion.mechanical_model()
    }
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn branch(&self) -> usize {
        self.branch
    }
    pub const fn coordinate_m(&self) -> f64 {
        self.motion.coordinate_m()
    }
    pub const fn velocity_m_s(&self) -> f64 {
        self.motion.velocity_m_s()
    }
    pub const fn capacitance(&self) -> CapacitanceSample {
        self.motion.capacitance()
    }
    pub const fn charge_c(&self) -> f64 {
        self.motion.charge_c()
    }
    pub const fn voltage_v(&self) -> f64 {
        self.motion.voltage_v()
    }
    pub const fn field_force_n(&self) -> f64 {
        self.motion.field_force_n()
    }
    pub const fn field_energy_j(&self) -> f64 {
        self.motion.field_energy_j()
    }
    pub const fn electrical_energy_j(&self) -> f64 {
        self.electrical_energy_j
    }

    /// Advance charge/flux, native body momentum and resolved geometry together.
    /// Every nonzero displacement trial re-solves the original dielectric mesh;
    /// its C and analytic shape derivative enter the same symmetric work rule.
    /// The native circuit owns KCL/KVL, source power and material/contact loss.
    /// The native articulated impulse owner owns body velocity and kinetic work.
    ///
    /// Material snapshots, topology, the other storage laws and mechanical
    /// inertia stay fixed. Heat is reported to explicit receivers, not silently
    /// deposited in an invented reservoir. A thermal/material update requires
    /// a separate joint owner. All original mesh/path restrictions apply.
    ///
    /// The work split is symmetric; timestep and mesh refinement remain
    /// necessary. Energy admission alone does not bound either error. Unresolved C
    /// differences use the endpoint analytic derivatives, with the independent
    /// physical work checks still enforced. Failure/cancellation changes nothing.
    pub fn step_candidate(
        &self,
        input: CircuitMotionInput<'_>,
        motion_controls: MotionControls,
        circuit_controls: MaterialCircuitSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CircuitMotionCandidate> {
        check(cancelled)?;
        if !positive(input.duration_s)
            || !input.external_force_n.is_finite()
            || !positive(motion_controls.energy_tolerance_j)
            || !positive(circuit_controls.energy_tolerance_j)
            || !positive(circuit_controls.power_tolerance_w)
        {
            return Err(Error::Invalid(
                "positive duration and work/power limits and finite force required",
            ));
        }
        let old = self.motion.electrical;
        let c0 = self.capacitance().capacitance_f();
        let [axis_x, axis_y] = self.path().translation_axis();
        let axis = Vec3::new(axis_x, axis_y, 0.0);
        let mut field_solves = 0;
        let evaluate = |displacement: f64, _trials: usize, cancelled: &mut _| -> Result<(_, f64)> {
            let coordinate = finite(self.coordinate_m() + displacement)?;
            if displacement != 0.0 && coordinate == self.coordinate_m() {
                return Err(Error::Invalid(
                    "electrode displacement disappeared in coordinate rounding",
                ));
            }
            let sample = if displacement == 0.0 {
                self.capacitance()
            } else {
                if coordinate != 0.0 {
                    field_solves += 1;
                }
                self.path().sample(coordinate, cancelled)?
            };
            let c1 = sample.capacitance_f();
            // Harmonic mean without overflowing C0+C1 or their product.
            let (small, large) = (c0.min(c1), c0.max(c1));
            let effective = finite(small / (0.5 + 0.5 * (small / large)))?;
            let (trial, _) = self.circuit.with_resolved_capacitance(
                self.branch,
                effective,
                input.electrical_ports,
                cancelled,
            )?;
            let electrical_candidate = trial.step_with_contacts_candidate(
                trial.resistors().to_vec(),
                trial.contacts().to_vec(),
                input.electrical_ports,
                input.duration_s,
                circuit_controls,
                &mut *cancelled,
            )?;
            let (electrical_candidate, electrical_step) = electrical_candidate.into_parts();
            let charge_c = electrical_candidate.charges_c()[self.charge_slot];
            let control = ElectricalControl::IsolatedCharge { charge_c };
            let electrical = control.at(sample)?;
            self.path().admit_amplitude(sample, electrical.voltage_v)?;
            let dc = finite(c1 - c0)?;
            let derivative = if displacement == 0.0 {
                self.capacitance().derivative_f_per_m()
            } else if dc.abs() <= 16.0 * f64::EPSILON * c0.max(c1) {
                finite(
                    0.5 * self.capacitance().derivative_f_per_m()
                        + 0.5 * sample.derivative_f_per_m(),
                )?
            } else {
                finite(dc / displacement)?
            };
            // Symmetric geometry force includes BOTH charge endpoints; using
            // Q0*Q1 here would lose work whenever circuit charge changes.
            let force = finite(
                (0.25 * (old.charge_c / c0)) * (derivative * (old.charge_c / c1))
                    + (0.25 * (charge_c / c1)) * (derivative * (charge_c / c0)),
            )?;
            let impulse = finite(input.duration_s * finite(force + input.external_force_n)?)?;
            let midpoint = finite(self.coordinate_m() + 0.5 * displacement)?;
            let mechanics = apply_prescribed_base_impulses(
                self.mechanical_model(),
                BaseState::stationary(Se3::identity()),
                &[midpoint],
                &[self.velocity_m_s()],
                &[0.0],
                &[Wrench::new(Vec3::default(), axis.scale(impulse))],
            )?;
            let velocity = mechanics.joint_velocity_after[0];
            let residual = finite(
                displacement
                    - input.duration_s * finite(0.5 * self.velocity_m_s() + 0.5 * velocity)?,
            )?;
            let field_work = finite(force * displacement)?;
            let external_work = finite(input.external_force_n * displacement)?;
            let (circuit, next_energy) = electrical_candidate.with_resolved_capacitance(
                self.branch,
                c1,
                input.electrical_ports,
                cancelled,
            )?;
            let change = Dd::from_f64(next_energy) - Dd::from_f64(self.electrical_energy_j);
            let electrical_residual = finite(
                (change + Dd::from_f64(field_work) + Dd::from_f64(electrical_step.heat_j)
                    - Dd::from_f64(electrical_step.native.supplied))
                .to_f64(),
            )?;
            let mechanical_residual = finite(
                (Dd::from_f64(mechanics.kinetic_energy_change)
                    - Dd::from_f64(field_work)
                    - Dd::from_f64(external_work))
                .to_f64(),
            )?;
            let total_residual = finite(
                (change
                    + Dd::from_f64(mechanics.kinetic_energy_change)
                    + Dd::from_f64(electrical_step.heat_j)
                    - Dd::from_f64(electrical_step.native.supplied)
                    - Dd::from_f64(external_work))
                .to_f64(),
            )?;
            let mut motion = self.motion.clone();
            motion.sample = sample;
            motion.control = control;
            motion.electrical = electrical;
            motion.velocity_m_s = velocity;
            let state = Self {
                motion,
                circuit,
                branch: self.branch,
                charge_slot: self.charge_slot,
                electrical_energy_j: next_energy,
            };
            let step = CircuitMotionStep {
                electrical_energy_j: [self.electrical_energy_j, next_energy],
                field_energy_j: [old.energy_j, electrical.energy_j],
                electrical_supply_j: electrical_step.native.supplied,
                circuit: electrical_step.native,
                mechanics,
                external_mechanical_work_j: external_work,
                field_work_on_body_j: field_work,
                discrete_field_force_n: force,
                displacement_m: displacement,
                kinematic_residual_m: residual,
                resistors: electrical_step.resistors,
                contacts: electrical_step.contacts,
                heat_j: electrical_step.heat_j,
                electrical_energy_residual_j: electrical_residual,
                mechanical_energy_residual_j: mechanical_residual,
                total_energy_residual_j: total_residual,
                field_solves,
            };
            Ok((CircuitMotionCandidate { state, step }, residual))
        };
        let candidate = feedback::solve_scalar(
            motion_controls.displacement_bracket_m,
            motion_controls.displacement_tolerance_m,
            motion_controls.maximum_iterations,
            evaluate,
            cancelled,
        )
        .map_err(|error| match error {
            ScalarError::Evaluation(error) => error,
            ScalarError::Invalid(reason) => Error::Invalid(reason),
            ScalarError::Cancelled => Error::Cancelled,
            ScalarError::OutsideBracket { lower, upper } => Error::OutsideBracket {
                lower_residual_m: lower,
                upper_residual_m: upper,
            },
            ScalarError::NotConverged {
                iterations,
                residual,
            } => Error::NotConverged {
                iterations,
                residual_m: residual,
            },
        })?;
        let step = &candidate.step;
        if [
            step.electrical_energy_residual_j,
            step.mechanical_energy_residual_j,
            step.total_energy_residual_j,
            step.mechanics.kinetic_work_residual,
        ]
        .iter()
        .any(|residual| residual.abs() > motion_controls.energy_tolerance_j)
        {
            return Err(Error::Balance {
                field_residual_j: step.electrical_energy_residual_j,
                mechanical_residual_j: step.mechanical_energy_residual_j,
                total_residual_j: step.total_energy_residual_j,
            });
        }
        check(cancelled)?;
        Ok(candidate)
    }
}
