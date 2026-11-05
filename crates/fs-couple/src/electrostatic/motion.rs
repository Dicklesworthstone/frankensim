//! Reciprocal guided motion driven by freshly solved finite-element fields.
//!
//! One affine mesh coordinate resolves a two-electrode capacitance C(s), and
//! the original field's virtual work supplies its analytic derivative. The
//! existing articulated impulse owner supplies body inertia, momentum and
//! kinetic work; the existing bounded scalar feedback owner finds one common
//! displacement. No second field, rigid-body or root solver is implemented.
//!
//! Charge is either held on an isolated electrode, or voltage is held by an
//! ideal source whose signed V dQ work is retained. The corresponding field
//! work is -d[Q²/(2C)] or +d[C V²/2]. The optional `circuit` owner supplies
//! finite-circuit charge evolution on this same field path. Breakdown, charge
//! transport inside the dielectric and radiation remain separate owners.
//! Energy conventions: <https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_theory.05.79.html>.

mod path;
pub use path::{CapacitanceSample, PlanarCapacitorPath};
/// Finite material RLC supply coupled to the same resolved moving field.
#[cfg(feature = "resolved-electrostatic-circuit-motion")]
pub mod circuit;

use std::sync::Arc;

use crate::cavity_phs::piston::thermal::feedback::{self, ScalarError};
use fs_em::electrostatic::ElectrostaticError;
use fs_ga::{Se3, Vec3, Wrench};
use fs_math::dd::Dd;
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointModel, Link, PrescribedBaseImpulseResponse,
    apply_prescribed_base_impulses,
};
use fs_mbd_material::ResolvedBody;

#[derive(Debug)]
pub enum Error {
    Field(ElectrostaticError),
    Mechanics(ArticulatedError),
    #[cfg(feature = "resolved-electrostatic-circuit-motion")]
    Circuit(crate::material_circuit::MaterialCircuitError),
    Invalid(&'static str),
    Cancelled,
    OutsideBracket {
        lower_residual_m: f64,
        upper_residual_m: f64,
    },
    NotConverged {
        iterations: usize,
        residual_m: f64,
    },
    Balance {
        field_residual_j: f64,
        mechanical_residual_j: f64,
        total_residual_j: f64,
    },
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved electrostatic motion: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<ElectrostaticError> for Error {
    fn from(error: ElectrostaticError) -> Self {
        match error {
            ElectrostaticError::Cancelled => Self::Cancelled,
            error => Self::Field(error),
        }
    }
}
impl From<ArticulatedError> for Error {
    fn from(error: ArticulatedError) -> Self {
        Self::Mechanics(error)
    }
}
#[cfg(feature = "resolved-electrostatic-circuit-motion")]
impl From<crate::material_circuit::MaterialCircuitError> for Error {
    fn from(error: crate::material_circuit::MaterialCircuitError) -> Self {
        match error {
            crate::material_circuit::MaterialCircuitError::Cancelled => Self::Cancelled,
            error => Self::Circuit(error),
        }
    }
}
type Result<T> = core::result::Result<T, Error>;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "unrepresentable electrostatic motion arithmetic",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Fixed electrical control over the entire retained history. Switching this
/// control requires a separate charge/source-work transaction, not replacement
/// of this enum in an already accepted state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElectricalControl {
    IsolatedCharge { charge_c: f64 },
    HeldVoltage { voltage_v: f64 },
}

#[derive(Clone, Copy, Debug)]
pub struct MotionControls {
    /// Proposed displacement relative to the accepted coordinate [m]. Must
    /// contain zero; all evaluated coordinates must stay in the mesh path.
    pub displacement_bracket_m: [f64; 2],
    pub displacement_tolerance_m: f64,
    /// Independent field/source, mechanical and combined work allowances [J].
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct MotionInput {
    pub duration_s: f64,
    /// Held external force in the positive electrode's translation direction.
    pub external_force_n: f64,
}

#[derive(Clone, Copy, Debug)]
struct ElectricalState {
    voltage_v: f64,
    charge_c: f64,
    energy_j: f64,
    force_n: f64,
}
impl ElectricalControl {
    fn at(self, sample: CapacitanceSample) -> Result<ElectricalState> {
        let capacitance = sample.capacitance_f();
        let (voltage_v, charge_c) = match self {
            Self::IsolatedCharge { charge_c } => (finite(charge_c / capacitance)?, charge_c),
            Self::HeldVoltage { voltage_v } => (voltage_v, finite(voltage_v * capacitance)?),
        };
        finite(voltage_v)?;
        finite(charge_c)?;
        let energy_j = finite((0.5 * charge_c) * voltage_v)?;
        let force_n = finite((0.5 * voltage_v) * (voltage_v * sample.derivative_f_per_m()))?;
        if energy_j < 0.0 || ((voltage_v != 0.0 || charge_c != 0.0) && energy_j == 0.0) {
            return Err(Error::Invalid(
                "nonzero electrical storage disappeared in SI arithmetic",
            ));
        }
        Ok(ElectricalState {
            voltage_v,
            charge_c,
            energy_j,
            force_n,
        })
    }
}

/// Accepted field geometry and one native, material-derived translating body.
/// The opposite electrode and guide are prescribed. The supplied body must
/// represent all inertia moving rigidly with the positive electrode; interior
/// mesh motion supplies no additional invented mass. Deforming material inertia
/// or moving supports require a richer mechanical mapping.
#[derive(Clone)]
pub struct PlanarCapacitorMotion {
    path: Arc<PlanarCapacitorPath>,
    body: Arc<ResolvedBody>,
    model: Arc<ArticulatedModel>,
    control: ElectricalControl,
    sample: CapacitanceSample,
    electrical: ElectricalState,
    velocity_m_s: f64,
}

/// Independent endpoint storage and actual native momentum/work observations.
/// The native inertia map applies the generalized axial force at the link
/// origin. It describes this guided coordinate; its support wrench is not a
/// reconstruction of the electrode's full electric traction or torque.
pub struct MotionStep {
    pub native: PrescribedBaseImpulseResponse,
    pub field_energy_j: [f64; 2],
    pub field_energy_change_j: f64,
    pub electrical_supply_j: f64,
    pub external_mechanical_work_j: f64,
    pub field_work_on_body_j: f64,
    pub discrete_field_force_n: f64,
    pub displacement_m: f64,
    pub kinematic_residual_m: f64,
    pub field_energy_residual_j: f64,
    pub mechanical_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    /// New coordinate field solves; the accepted old response is reused.
    pub field_solves: usize,
}
pub struct MotionCandidate {
    state: PlanarCapacitorMotion,
    step: MotionStep,
}
impl MotionCandidate {
    pub const fn state(&self) -> &PlanarCapacitorMotion {
        &self.state
    }
    pub const fn step(&self) -> &MotionStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarCapacitorMotion, MotionStep) {
        (self.state, self.step)
    }
}

impl PlanarCapacitorMotion {
    pub fn new(
        path: PlanarCapacitorPath,
        body: ResolvedBody,
        control: ElectricalControl,
        coordinate_m: f64,
        velocity_m_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        finite(velocity_m_s)?;
        let sample = path.sample(coordinate_m, cancelled)?;
        let electrical = control.at(sample)?;
        path.admit_amplitude(sample, electrical.voltage_v)?;
        let [x, y] = path.translation_axis();
        let model = ArticulatedModel::new(vec![Link::new(
            "resolved-electrostatic-electrode",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(x, y, 0.0), None)?,
            body.inertia(),
        )])?;
        check(cancelled)?;
        Ok(Self {
            path: Arc::new(path),
            body: Arc::new(body),
            model: Arc::new(model),
            control,
            sample,
            electrical,
            velocity_m_s,
        })
    }
    pub fn path(&self) -> &PlanarCapacitorPath {
        &self.path
    }
    pub fn body(&self) -> &ResolvedBody {
        &self.body
    }
    pub fn mechanical_model(&self) -> &ArticulatedModel {
        &self.model
    }
    pub const fn control(&self) -> ElectricalControl {
        self.control
    }
    pub const fn coordinate_m(&self) -> f64 {
        self.sample.coordinate_m()
    }
    pub const fn velocity_m_s(&self) -> f64 {
        self.velocity_m_s
    }
    pub const fn capacitance(&self) -> CapacitanceSample {
        self.sample
    }
    pub const fn voltage_v(&self) -> f64 {
        self.electrical.voltage_v
    }
    pub const fn charge_c(&self) -> f64 {
        self.electrical.charge_c
    }
    pub const fn field_energy_j(&self) -> f64 {
        self.electrical.energy_j
    }
    pub const fn field_force_n(&self) -> f64 {
        self.electrical.force_n
    }

    /// Propose one implicit midpoint displacement with a secant electric force.
    /// Each nonzero trial rebuilds and solves the ORIGINAL dielectric field at
    /// its proposed geometry. Force and source work come from that same C(s).
    /// Native articulated momentum and independent endpoint field energies must
    /// close three work checks. Refusal/cancellation mutates no accepted state.
    ///
    /// The secant is the discrete energy gradient, not a finite-difference
    /// substitute for the instantaneous force API. When capacitance differences
    /// are unresolved at floating-point precision, the endpoint analytic shape
    /// derivatives supply a symmetric local force; all work gates still apply.
    /// No numerical discrepancy is converted to heat or projected out of state.
    /// One step performs at most maximum_iterations+1 new field solves, each
    /// with the original assembly and Krylov budgets. No wall-clock claim is made.
    pub fn step_candidate(
        &self,
        input: MotionInput,
        controls: MotionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MotionCandidate> {
        check(cancelled)?;
        if !positive(input.duration_s)
            || !input.external_force_n.is_finite()
            || !positive(controls.energy_tolerance_j)
        {
            return Err(Error::Invalid(
                "positive duration and work allowance and finite external force required",
            ));
        }
        let [axis_x, axis_y] = self.path.translation_axis();
        let axis = Vec3::new(axis_x, axis_y, 0.0);
        let old = self.electrical;
        let mut field_solves = 0;
        let evaluate = |displacement: f64,
                        _trials: usize,
                        cancelled: &mut _|
         -> Result<(MotionCandidate, f64)> {
            let coordinate = finite(self.coordinate_m() + displacement)?;
            if displacement != 0.0 && coordinate == self.coordinate_m() {
                return Err(Error::Invalid(
                    "nonzero electrode displacement disappeared in coordinate rounding",
                ));
            }
            let sample = if displacement == 0.0 {
                self.sample
            } else {
                if coordinate != 0.0 {
                    field_solves += 1;
                }
                self.path.sample(coordinate, cancelled)?
            };
            let electrical = self.control.at(sample)?;
            self.path.admit_amplitude(sample, electrical.voltage_v)?;
            let dc = finite(sample.capacitance_f() - self.sample.capacitance_f())?;
            let force = if displacement == 0.0 {
                old.force_n
            } else if dc.abs()
                <= 16.0 * f64::EPSILON * self.sample.capacitance_f().max(sample.capacitance_f())
            {
                finite(0.5 * old.force_n + 0.5 * electrical.force_n)?
            } else {
                // Both controls give W_field=V0*V1*(C1-C0)/2. For held
                // voltage this is +delta U; for fixed charge it is -delta U.
                finite((0.5 * old.voltage_v) * (electrical.voltage_v * (dc / displacement)))?
            };
            let impulse = finite(input.duration_s * finite(force + input.external_force_n)?)?;
            let midpoint = finite(self.coordinate_m() + 0.5 * displacement)?;
            let native = apply_prescribed_base_impulses(
                &self.model,
                BaseState::stationary(Se3::identity()),
                &[midpoint],
                &[self.velocity_m_s],
                &[0.0],
                &[Wrench::new(Vec3::default(), axis.scale(impulse))],
            )?;
            let velocity = native.joint_velocity_after[0];
            let midpoint_velocity = finite(0.5 * self.velocity_m_s + 0.5 * velocity)?;
            let residual = finite(displacement - input.duration_s * midpoint_velocity)?;
            let field_work = finite(force * displacement)?;
            let external_work = finite(input.external_force_n * displacement)?;
            let field_change =
                finite((Dd::from_f64(electrical.energy_j) - Dd::from_f64(old.energy_j)).to_f64())?;
            let electrical_supply = match self.control {
                ElectricalControl::IsolatedCharge { .. } => 0.0,
                ElectricalControl::HeldVoltage { voltage_v } => finite(
                    (Dd::from_f64(voltage_v)
                        * (Dd::from_f64(electrical.charge_c) - Dd::from_f64(old.charge_c)))
                    .to_f64(),
                )?,
            };
            let field_residual = finite(
                (Dd::from_f64(field_change) + Dd::from_f64(field_work)
                    - Dd::from_f64(electrical_supply))
                .to_f64(),
            )?;
            let mechanical_residual = finite(
                (Dd::from_f64(native.kinetic_energy_change)
                    - Dd::from_f64(field_work)
                    - Dd::from_f64(external_work))
                .to_f64(),
            )?;
            let total_residual = finite(
                (Dd::from_f64(field_change) + Dd::from_f64(native.kinetic_energy_change)
                    - Dd::from_f64(electrical_supply)
                    - Dd::from_f64(external_work))
                .to_f64(),
            )?;
            let mut state = self.clone();
            state.sample = sample;
            state.electrical = electrical;
            state.velocity_m_s = velocity;
            let step = MotionStep {
                native,
                field_energy_j: [old.energy_j, electrical.energy_j],
                field_energy_change_j: field_change,
                electrical_supply_j: electrical_supply,
                external_mechanical_work_j: external_work,
                field_work_on_body_j: field_work,
                discrete_field_force_n: force,
                displacement_m: displacement,
                kinematic_residual_m: residual,
                field_energy_residual_j: field_residual,
                mechanical_energy_residual_j: mechanical_residual,
                total_energy_residual_j: total_residual,
                field_solves,
            };
            Ok((MotionCandidate { state, step }, residual))
        };
        let candidate = feedback::solve_scalar(
            controls.displacement_bracket_m,
            controls.displacement_tolerance_m,
            controls.maximum_iterations,
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
            step.field_energy_residual_j,
            step.mechanical_energy_residual_j,
            step.total_energy_residual_j,
            step.native.kinetic_work_residual,
        ]
        .iter()
        .any(|x| x.abs() > controls.energy_tolerance_j)
        {
            return Err(Error::Balance {
                field_residual_j: step.field_energy_residual_j,
                mechanical_residual_j: step.mechanical_energy_residual_j,
                total_residual_j: step.total_energy_residual_j,
            });
        }
        check(cancelled)?;
        Ok(candidate)
    }
}
