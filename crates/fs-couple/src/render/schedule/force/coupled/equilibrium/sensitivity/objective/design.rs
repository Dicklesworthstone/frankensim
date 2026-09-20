//! Multi-load-case inverse design over the existing equilibrium and adjoint owners.
//! Candidates are rebuilt from an immutable zero-state template. No optimizer,
//! finite-difference gradient, warm-up dynamics or material identification lives
//! here. A result contains every case or none; spent work is never refunded.
use super::*;
use crate::render::schedule::force::{coupled::contact::multiple::MultiContactConfig, project};
use crate::render::RenderError;
use fs_dcontact::Obstacle;

/// A field assigned a physical value by a design variable. Indices refer to
/// construction order, not the flattened modal layout. Damping is deliberately
/// absent: static observations cannot identify it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DesignField {
    /// Spring stiffness [N/m].
    SpringStiffness(usize),
    /// Stress-free spring extension [m].
    SpringRest(usize),
    /// Original normal-law coefficient [N/m^alpha].
    ContactStiffness(usize),
    /// Original contact gap [m].
    ContactGap(usize),
    /// Original quadrature weight, with the source law's declared convention.
    ContactWeight(usize),
    /// Physical force [N] on one actuator of one load case.
    ActuatorForce { case: usize, actuator: usize },
}

/// One dimensionless decision x, assigning p = reference + scale*x to every
/// listed field. Sharing requires identical parameter kinds (and identical
/// exponents for contact stiffness). These are domain bounds, not a box-KKT
/// solver. Reference values need not equal the original template coefficients.
#[derive(Clone, Debug)]
pub struct DesignVariable {
    /// Unique nonempty name, at most 128 bytes.
    pub name: String,
    /// Physical value at x=0.
    pub reference: f64,
    /// Positive physical units per unit decision coordinate.
    pub scale: f64,
    /// Finite lower bound on the physical value.
    pub minimum: f64,
    /// Finite upper bound, strictly greater than minimum.
    pub maximum: f64,
    /// Nonempty set of distinct fields assigned this value.
    pub fields: Vec<DesignField>,
}

/// A signed, additive physical load. Its shape maps N to N/sqrt(kg).
#[derive(Clone, Debug)]
pub struct DesignLoad {
    /// Explicit mass-normalized physical actuator map.
    pub attachment: ModalAttachment,
    /// Held physical force [N].
    pub force_n: f64,
}

/// Independent stationary experiment, not the continuation of another case.
#[derive(Clone, Debug)]
pub struct DesignLoadCase {
    /// Unique nonempty name, at most 128 bytes.
    pub name: String,
    /// Additive loads; multiple loads may act on the same attachment.
    pub loads: Vec<DesignLoad>,
    /// Existing physical displacement objectives, with their own scales/weights.
    pub targets: Vec<DisplacementTarget>,
}

/// Existing physics/derivative budgets plus bounded design-family dimensions.
#[derive(Clone, Copy, Debug)]
pub struct DesignBudget {
    /// Original mechanical admission and balance limits.
    pub coupling: ModalCouplingConfig,
    /// Original joint contact preload limits.
    pub contact: MultiContactConfig,
    /// Original primal/adjoint and activity-margin limits.
    pub sensitivity: SensitivityBudget,
    /// In 1..=64.
    pub max_cases: usize,
    /// In 1..=128.
    pub max_variables: usize,
    /// Total variable-to-field bindings, at most 1024.
    pub max_bindings: usize,
    /// Targets plus actuators per case, at most 1024.
    pub max_ports_per_case: usize,
}

/// A structural, physical or resource refusal; never a substitute objective.
#[derive(Debug)]
pub enum DesignError {
    /// Incompatible or nonfinite problem data.
    Invalid { what: &'static str },
    /// Numerical physical value outside the declared finite domain.
    OutsideBounds { variable: usize, value: f64, minimum: f64, maximum: f64 },
    /// Cancellation observed before returning a complete result.
    Cancelled,
    /// A cumulative allowance has been exhausted.
    Budget { what: &'static str },
    /// The original force projection error, not a fabricated bad objective.
    Projection { case: usize, source: RenderError },
    /// A complete-case solve, activity-margin or adjoint refusal.
    Case { case: usize, source: ModalCouplingError },
}
impl core::fmt::Display for DesignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "equilibrium design: {what}"),
            Self::OutsideBounds { variable, value, minimum, maximum } => write!(f,
                "design variable {variable}: {value} outside [{minimum}, {maximum}]"),
            Self::Cancelled => write!(f, "equilibrium design cancelled"),
            Self::Budget { what } => write!(f, "equilibrium design budget: {what}"),
            Self::Projection { case, source } => write!(f, "load case {case} projection: {source}"),
            Self::Case { case, source } => write!(f, "load case {case}: {source}"),
        }
    }
}
impl std::error::Error for DesignError {}

/// Attempts, including failures. Primal/adjoint internals retain their original
/// operation budgets; these counts describe calls, not exact flops or elapsed time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DesignWork {
    /// Objective callback attempts, including invalid candidate points.
    pub evaluations: usize,
    /// Attempted primal-plus-adjoint load cases.
    pub case_solves: usize,
}

/// Caller-owned cumulative limits. Invalid candidates consume an evaluation;
/// a case consumes its solve allowance before physical work starts. A cancelled
/// pre-call consumes neither. Raising limits never erases accumulated work.
#[derive(Clone, Debug)]
pub struct DesignControl {
    maximum_evaluations: usize,
    maximum_case_solves: usize,
    work: DesignWork,
}
impl DesignControl {
    #[must_use]
    /// Start explicit cumulative limits (zero deliberately permits no work).
    pub fn new(maximum_evaluations: usize, maximum_case_solves: usize) -> Self {
        Self { maximum_evaluations, maximum_case_solves, work: DesignWork::default() }
    }
    #[must_use]
    /// Inspect spent work without changing limits.
    pub const fn work(&self) -> DesignWork { self.work }
    /// Raise cumulative allowances without refunding failed attempts.
    pub fn extend(&mut self, maximum_evaluations: usize, maximum_case_solves: usize) -> Result<(), DesignError> {
        if maximum_evaluations < self.maximum_evaluations || maximum_case_solves < self.maximum_case_solves {
            return Err(bad("cumulative design limits may only be extended"));
        }
        self.maximum_evaluations = maximum_evaluations;
        self.maximum_case_solves = maximum_case_solves;
        Ok(())
    }
}

/// Completed case evidence, in input case order.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignCaseResult {
    /// Dimensionless summed displacement objective.
    pub value: f64,
    /// Predictions in physical metres, in target order.
    pub observations_m: Vec<f64>,
    /// Original force-balance and contact-margin check.
    pub equilibrium: EquilibriumLinearizationReport,
    /// Full-coordinate relative adjoint residual.
    pub adjoint_relative_residual: f64,
}
/// Gradients are with respect to the dimensionless decisions, including both
/// physical scaling and summation over shared fields and independent cases.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignEvaluation {
    /// Dimensionless summed displacement objective.
    pub value: f64,
    /// One total derivative per scaled decision coordinate.
    pub gradient: Vec<f64>,
    /// Decoded parameter values in variable declaration order.
    pub physical_parameters: Vec<f64>,
    /// Every completed experiment, never a partial family.
    pub cases: Vec<DesignCaseResult>,
}

/// Immutable, reusable physical problem, with no retained speculative states.
/// Models must have zero initial Q/V. Modal bases, masses, observation maps,
/// exponents and damping stay fixed; only explicitly bound fields are changed.
pub struct EquilibriumDesign {
    models: Vec<ModalAcousticTimeModel>,
    springs: Vec<ModalConnection>,
    contacts: Vec<(ModalContact, ModalContactConfig)>,
    cases: Vec<DesignLoadCase>,
    variables: Vec<DesignVariable>,
    budget: DesignBudget,
    modes: usize,
    offsets: Vec<usize>,
}
impl EquilibriumDesign {
    /// Admit the immutable model, distinct field bindings and complete cases.
    pub fn new(models: Vec<ModalAcousticTimeModel>, springs: Vec<ModalConnection>,
        contacts: Vec<(ModalContact, ModalContactConfig)>, cases: Vec<DesignLoadCase>,
        variables: Vec<DesignVariable>, budget: DesignBudget, gate: &CancelGate) -> Result<Self, DesignError>
    {
        checkpoint(gate)?;
        if !(1..=64).contains(&budget.max_cases) || cases.is_empty() || cases.len() > budget.max_cases
            || !(1..=128).contains(&budget.max_variables) || variables.is_empty() || variables.len() > budget.max_variables
            || budget.max_bindings > 1024 || budget.max_ports_per_case > 1024
            || budget.contact.max_contacts > 32 || contacts.len() > budget.contact.max_contacts
            || springs.len() > budget.coupling.max_connections {
            return Err(bad("case, variable, binding, port or contact admission exceeded"));
        }
        if models.iter().flat_map(|m| m.states()).any(|s|
            s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(bad("design templates require zero Q/V; no supplied vibration may be discarded"));
        }
        // Reuse structural admission before any candidate is attempted.
        let template = CoupledModalSystem::new(models, springs.clone(), budget.coupling, gate)
            .map_err(|source| case_error(0, source))?;
        for (contact, config) in &contacts {
            contact_column(&template, contact, *config).map_err(|source| case_error(0, source))?;
        }
        let mut seen = Vec::new();
        for (i, variable) in variables.iter().enumerate() {
            checkpoint(gate)?;
            if !name_ok(&variable.name) || variables[..i].iter().any(|v| v.name == variable.name)
                || [variable.reference, variable.scale, variable.minimum, variable.maximum].iter().any(|x| !x.is_finite())
                || variable.scale <= 0.0 || variable.minimum >= variable.maximum
                || variable.reference < variable.minimum || variable.reference > variable.maximum || variable.fields.is_empty() {
                return Err(bad("variable needs a unique name, positive scale, ordered finite bounds and fields"));
            }
            for field in &variable.fields {
                if seen.len() == budget.max_bindings || seen.contains(field) {
                    return Err(bad("duplicate field assignment or binding budget exceeded"));
                }
                validate_field(*field, &springs, &contacts, &cases, variable.minimum)?;
                if !same_quantity(variable.fields[0], *field, &contacts) {
                    return Err(bad("shared variable fields must have the same physical parameter kind and exponent"));
                }
                seen.push(*field);
            }
        }
        for (i, case) in cases.iter().enumerate() {
            checkpoint(gate)?;
            if !name_ok(&case.name) || cases[..i].iter().any(|c| c.name == case.name) || case.targets.is_empty()
                || case.loads.len().checked_add(case.targets.len()).is_none_or(|n| n > budget.max_ports_per_case) {
                return Err(bad("case needs a unique name and a bounded nonempty target family"));
            }
            for load in &case.loads {
                check_map(&template, &load.attachment)?;
                if !load.force_n.is_finite() { return Err(bad("physical load must be finite")); }
            }
            for target in &case.targets {
                check_map(&template, &target.attachment)?;
                if !target.target_m.is_finite() || !target.scale_m.is_finite() || target.scale_m <= 0.0
                    || !target.weight.is_finite() || target.weight < 0.0 { return Err(bad("invalid physical displacement objective")); }
            }
        }
        checkpoint(gate)?;
        Ok(Self { modes: template.mode_count(), offsets: template.offsets.clone(), models: template.models,
            springs, contacts, cases, variables, budget })
    }
    #[must_use]
    /// Decision definitions and their physical scaling.
    pub fn variables(&self) -> &[DesignVariable] { &self.variables }
    #[must_use]
    /// Original independent experiments, never mutated by a candidate.
    pub fn load_cases(&self) -> &[DesignLoadCase] { &self.cases }

    /// Decode the complete decision vector before cloning/solving any physics.
    pub fn physical_parameters(&self, point: &[f64]) -> Result<Vec<f64>, DesignError> {
        if point.len() != self.variables.len() || point.iter().any(|x| !x.is_finite()) {
            return Err(bad("one finite decision coordinate is required per variable"));
        }
        point.iter().zip(&self.variables).enumerate().map(|(i, (x, v))| {
            let value = v.reference + v.scale * x;
            if !value.is_finite() || value < v.minimum || value > v.maximum {
                Err(DesignError::OutsideBounds { variable: i, value, minimum: v.minimum, maximum: v.maximum })
            } else { Ok(value) }
        }).collect()
    }

    /// Re-solve each independent case at the trial parameters, then use one
    /// existing adjoint per case. No partial objective/gradient escapes. A case
    /// changing contact activity is re-admitted from its own solved state; a
    /// near-switch derivative still refuses under the original margin rule.
    pub fn evaluate(&self, point: &[f64], control: &mut DesignControl, gate: &CancelGate)
        -> Result<DesignEvaluation, DesignError>
    {
        checkpoint(gate)?;
        if control.work.evaluations >= control.maximum_evaluations { return Err(DesignError::Budget { what: "evaluations" }); }
        control.work.evaluations += 1;
        let parameters = self.physical_parameters(point)?;
        if self.cases.len() > control.maximum_case_solves.saturating_sub(control.work.case_solves) {
            return Err(DesignError::Budget { what: "complete load-case family" });
        }
        let mut springs = self.springs.clone();
        let mut cases = self.cases.clone();
        let mut coefficients: Vec<_> = self.contacts.iter().map(|(c, _)|
            [c.law.stiffness(), c.law.gaps()[0], c.law.weights()[0]]).collect();
        for (variable, &value) in self.variables.iter().zip(&parameters) {
            for field in &variable.fields {
                match *field {
                    DesignField::SpringStiffness(i) => springs[i].stiffness_n_m = value,
                    DesignField::SpringRest(i) => springs[i].rest_extension_m = value,
                    DesignField::ContactStiffness(i) => coefficients[i][0] = value,
                    DesignField::ContactGap(i) => coefficients[i][1] = value,
                    DesignField::ContactWeight(i) => coefficients[i][2] = value,
                    DesignField::ActuatorForce { case, actuator } => cases[case].loads[actuator].force_n = value,
                }
            }
        }
        let mut contacts = self.contacts.clone();
        for ((c, _), values) in contacts.iter_mut().zip(coefficients) {
            let old = &c.law;
            if values != [old.stiffness(), old.gaps()[0], old.weights()[0]] {
                // Modified coefficients are predictions, not the original receipt.
                c.law = Obstacle::new(vec![-1.0], 1, 1, vec![values[1]], vec![values[2]], values[0], old.alpha(),
                    format!("inverse-design candidate from {}", old.provenance()))
                    .and_then(|o| o.with_internal_loss(old.internal_loss()))
                    .map_err(|source| case_error(0, ModalCouplingError::ContactLaw(source)))?;
            }
        }
        let mut result = DesignEvaluation { value: 0.0, gradient: vec![0.0; self.variables.len()],
            physical_parameters: parameters, cases: Vec::with_capacity(cases.len()) };
        for (index, case) in cases.iter().enumerate() {
            checkpoint(gate)?;
            control.work.case_solves += 1;
            let mut network = CoupledModalSystem::new(self.models.clone(), springs.clone(), self.budget.coupling, gate)
                .map_err(|source| case_error(index, source))?;
            let mut columns = Vec::with_capacity(case.loads.len());
            let mut forces = Vec::with_capacity(case.loads.len());
            let mut actuators = Vec::with_capacity(case.loads.len());
            for load in &case.loads {
                let mut column = vec![0.0; self.modes];
                let start = self.offsets[load.attachment.component];
                column[start..start + load.attachment.shapes.len()].copy_from_slice(&load.attachment.shapes);
                columns.push(column); forces.push(load.force_n); actuators.push(load.attachment.clone());
            }
            // Reuse the same physical-to-modal force projection as rendering.
            let external = if columns.is_empty() { vec![0.0; self.modes] } else {
                project(&columns, &forces).map_err(|source| DesignError::Projection { case: index, source })?
            };
            if contacts.is_empty() { network.initialize_static_equilibrium(&external, gate) }
            else { network.initialize_contact_equilibrium(&external, &contacts, self.budget.contact, gate) }
                .map_err(|source| case_error(index, source))?;
            let linear = EquilibriumLinearization::new(&network, &external, &contacts, self.budget.sensitivity, gate)
                .map_err(|source| case_error(index, source))?;
            let objective = linear.displacement_objective(&case.targets, &actuators, self.budget.max_ports_per_case, gate)
                .map_err(|source| case_error(index, source))?;
            result.value = number(result.value + objective.value)?;
            for (j, variable) in self.variables.iter().enumerate() {
                let mut physical = 0.0;
                for &field in &variable.fields {
                    let pullback = &objective.residual_pullback;
                    let derivative = match field {
                        DesignField::SpringStiffness(i) => -pullback.springs[i].stiffness,
                        DesignField::SpringRest(i) => -pullback.springs[i].rest_extension,
                        DesignField::ContactStiffness(i) => -pullback.contacts[i].stiffness,
                        DesignField::ContactGap(i) => -pullback.contacts[i].gap,
                        DesignField::ContactWeight(i) => -pullback.contacts[i].weight,
                        DesignField::ActuatorForce { case, actuator } => if case == index {
                            objective.physical_force_gradient[actuator]
                        } else { 0.0 },
                    };
                    physical = number(physical + derivative)?;
                }
                result.gradient[j] = number(result.gradient[j] + number(physical * variable.scale)?)?;
            }
            result.cases.push(DesignCaseResult { value: objective.value, observations_m: objective.observations_m,
                equilibrium: linear.report(), adjoint_relative_residual: objective.adjoint.relative_residual });
        }
        checkpoint(gate)?;
        Ok(result)
    }
}

fn check_map(network: &CoupledModalSystem, map: &ModalAttachment) -> Result<(), DesignError> {
    let model = network.models.get(map.component).ok_or_else(|| bad("design load/target names an unknown component"))?;
    if map.shapes.len() != model.modes().len() || map.shapes.iter().any(|x| !x.is_finite()) {
        return Err(bad("design load/target must match the complete finite modal basis"));
    }
    Ok(())
}
fn validate_field(field: DesignField, springs: &[ModalConnection], contacts: &[(ModalContact, ModalContactConfig)],
    cases: &[DesignLoadCase], minimum: f64) -> Result<(), DesignError>
{
    let valid = match field {
        DesignField::SpringStiffness(i) => i < springs.len() && minimum >= 0.0,
        DesignField::SpringRest(i) => i < springs.len(),
        DesignField::ContactStiffness(i) | DesignField::ContactWeight(i) => i < contacts.len() && minimum >= 0.0,
        DesignField::ContactGap(i) => i < contacts.len(),
        DesignField::ActuatorForce { case, actuator } => cases.get(case).is_some_and(|c| actuator < c.loads.len()),
    };
    if valid { Ok(()) } else { Err(bad("unknown design field or negative stiffness/weight domain")) }
}
fn same_quantity(a: DesignField, b: DesignField, contacts: &[(ModalContact, ModalContactConfig)]) -> bool {
    match (a, b) {
        (DesignField::ContactStiffness(i), DesignField::ContactStiffness(j)) => contacts[i].0.law.alpha() == contacts[j].0.law.alpha(),
        _ => std::mem::discriminant(&a) == std::mem::discriminant(&b),
    }
}
fn name_ok(name: &str) -> bool { !name.trim().is_empty() && name.len() <= 128 }
fn number(x: f64) -> Result<f64, DesignError> { if x.is_finite() { Ok(x) } else { Err(bad("design arithmetic is not finite")) } }
fn checkpoint(gate: &CancelGate) -> Result<(), DesignError> { if gate.is_requested() { Err(DesignError::Cancelled) } else { Ok(()) } }
fn bad(what: &'static str) -> DesignError { DesignError::Invalid { what } }
fn case_error(case: usize, source: ModalCouplingError) -> DesignError {
    if matches!(&source, ModalCouplingError::Cancelled) { DesignError::Cancelled } else { DesignError::Case { case, source } }
}
