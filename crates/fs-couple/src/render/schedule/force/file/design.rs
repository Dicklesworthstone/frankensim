//! Authored mechanical models plus independent experiments -> static inverse design.
//! Reuses the performance decoder, source port maps, design-field admission,
//! nonlinear equilibrium and adjoints. No second model parser or optimizer.
use super::{
    MAX_PORT_WEIGHTS, MAX_PROJECTION_TERMS, MODAL_COUPLED_PERFORMANCE_SCHEMA,
    MODAL_CONTACT_PERFORMANCE_SCHEMA, MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA,
    ModalForceVoice, ModalPerformanceError, ModalPerformanceInfo, ParsedPerformance,
    Reader, Row, input,
};
use crate::render::schedule::force::{ForceInitialization, coupled::ModalAttachment};
use crate::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use crate::render::schedule::force::coupled::equilibrium::sensitivity::{
    SensitivityBudget, objective::DisplacementTarget,
};
use crate::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignBudget, DesignError, DesignField, DesignLoad, DesignLoadCase,
    DesignVariable, EquilibriumDesign,
    constraints::{ResponseConstraint, ResponseQuantity, ConstraintSense},
};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;

/// Design-file schema. Model files retain their existing v2/v3/v4 schemas.
pub const EQUILIBRIUM_DESIGN_SCHEMA: &str = "frankensim-equilibrium-design-v1";
/// Domain of the exact design bytes, separate from the source model identity.
pub const EQUILIBRIUM_DESIGN_HASH_DOMAIN: &str = "org.frankensim.fs-couple.equilibrium-design-input.v1";
/// Apply this byte ceiling while reading, before allocating an unbounded file.
pub const MAX_EQUILIBRIUM_DESIGN_BYTES: usize = 1024 * 1024;

/// Distinguish the offending input and preserve the original physical refusal.
#[derive(Debug)]
pub enum DesignFileError {
    /// Existing model decoder or model-template restriction.
    Model(ModalPerformanceError),
    /// Design syntax, counts, or referenced port.
    Input(ModalPerformanceError),
    /// Existing physical design admission or cancellation.
    Design(DesignError),
}
impl core::fmt::Display for DesignFileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Model(e) => write!(f, "design model: {e}"),
            Self::Input(ModalPerformanceError::Input { line, what }) => write!(f, "design input line {line}: {what}"),
            Self::Input(e) => write!(f, "design input: {e}"),
            Self::Design(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for DesignFileError {}
impl From<ModalPerformanceError> for DesignFileError {
    fn from(e: ModalPerformanceError) -> Self { Self::Input(e) }
}
impl From<DesignError> for DesignFileError {
    fn from(e: DesignError) -> Self { Self::Design(e) }
}

/// Complete immutable problem and the exact two source identities. Loading does
/// not optimize, spend an objective evaluation, or certify an equilibrium.
pub struct EquilibriumDesignFile {
    problem: EquilibriumDesign,
    model_info: ModalPerformanceInfo,
    design_hash: ContentHash,
}
impl EquilibriumDesignFile {
    /// Load a zero-state v2/v3/v4 model and explicit load cases/design bindings.
    ///
    /// A template must use retain-state/free-mass voices, zero Q/V, zero initial
    /// port forces and no events. Preloads, vibration or friction are not silently
    /// discarded. Model clocks/pressure settings keep their parsing admission but
    /// do not define this stationary displacement objective. Source per-component,
    /// connection and contact limits remain unchanged. Joint preload limits may
    /// tighten, never enlarge, an existing v4 contact-set budget.
    ///
    /// Shape copies across loads, targets and constraints are capped at 65536 coefficients,
    /// independently of the source's shape cap. Counts are checked before their
    /// allocations. Displacement rows name explicitly declared source ports.
    /// An optional constraint_limits/constraints section follows all variables;
    /// absent sections preserve the original unconstrained design semantics.
    pub fn from_bytes(model_bytes: &[u8], design_bytes: &[u8], gate: &CancelGate)
        -> Result<Self, DesignFileError>
    {
        poll(gate)?;
        if design_bytes.len() > MAX_EQUILIBRIUM_DESIGN_BYTES {
            return Err(input(1, "design input exceeds 1 MiB").into());
        }
        let parsed = ParsedPerformance::from_bytes(model_bytes, 1).map_err(DesignFileError::Model)?;
        poll(gate)?;
        if !matches!(parsed.info.schema, MODAL_COUPLED_PERFORMANCE_SCHEMA |
            MODAL_CONTACT_PERFORMANCE_SCHEMA | MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA) {
            return Err(model_error("static design requires a v2/v3/v4 model with explicit coupling limits; friction is unsupported"));
        }
        if !parsed.events.is_empty() || parsed.voices.iter().any(|voice|
            voice.initialization != ForceInitialization::RetainState
            || voice.initial_force_n.iter().any(|f| *f != 0.0)
            || voice.model.states().iter().any(|s|
                s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0)) {
            return Err(model_error("design templates require retained zero states, zero port loads and no scheduled events"));
        }
        let (connections, coupling) = parsed.coupled.ok_or_else(|| model_error("missing model coupling limits"))?;
        let text = std::str::from_utf8(design_bytes).map_err(|_| input(1, "design input must be UTF-8"))?;
        let mut reader = Reader { lines: text.lines(), line: 0 };
        reader.row(EQUILIBRIUM_DESIGN_SCHEMA)?.finish()?;
        let mut row = reader.row("preload_limits")?;
        let preload = MultiContactConfig {
            max_contacts: row.count(32)?, max_sweeps: row.count(128)?,
            max_setup_terms: row.count(MAX_PROJECTION_TERMS)?,
        };
        row.finish()?;
        if preload.max_sweeps == 0 || preload.max_contacts < parsed.info.contacts {
            return Err(input(reader.line, "preload needs positive sweeps and capacity for every model contact").into());
        }
        if let Some((_, original)) = &parsed.multiple {
            if preload.max_contacts > original.max_contacts || preload.max_sweeps > original.max_sweeps
                || preload.max_setup_terms > original.max_setup_terms {
                return Err(input(reader.line, "preload limits must not enlarge the model contact-set limits").into());
            }
        }
        let mut row = reader.row("sensitivity_limits")?;
        let sensitivity = SensitivityBudget {
            max_contacts: row.count(32)?, max_setup_terms: row.count(MAX_PROJECTION_TERMS)?,
            max_query_terms: row.count(MAX_PROJECTION_TERMS)?, minimum_contact_margin_m: row.scalar()?,
        };
        row.finish()?;
        if sensitivity.minimum_contact_margin_m < 0.0 || sensitivity.max_contacts < parsed.info.contacts {
            return Err(input(reader.line, "sensitivity needs a nonnegative activity margin and every model contact").into());
        }
        let mut row = reader.row("design_limits")?;
        let budget = DesignBudget { coupling, contact: preload, sensitivity,
            max_cases: row.count(64)?, max_variables: row.count(128)?,
            max_bindings: row.count(1024)?, max_ports_per_case: row.count(1024)?,
        };
        row.finish()?;
        let mut row = reader.row("cases")?;
        let case_count = row.count(budget.max_cases)?;
        row.finish()?;
        if case_count == 0 { return Err(input(reader.line, "at least one complete load case is required").into()); }
        let mut cases = Vec::with_capacity(case_count);
        let mut remaining_shapes = MAX_PORT_WEIGHTS;
        for _ in 0..case_count {
            poll(gate)?;
            let mut row = reader.row("case")?;
            let name = name(&mut row)?;
            let load_count = row.count(budget.max_ports_per_case)?;
            let target_count = row.count(budget.max_ports_per_case - load_count)?;
            row.finish()?;
            if target_count == 0 { return Err(input(reader.line, "each case needs displacement targets").into()); }
            let mut loads = Vec::with_capacity(load_count);
            for _ in 0..load_count {
                poll(gate)?;
                let mut row = reader.row("load")?;
                let attachment = port(&mut row, &parsed.voices, &mut remaining_shapes)?;
                let force_n = row.scalar()?;
                row.finish()?;
                loads.push(DesignLoad { attachment, force_n });
            }
            let mut targets = Vec::with_capacity(target_count);
            for _ in 0..target_count {
                poll(gate)?;
                let mut row = reader.row("target")?;
                let attachment = port(&mut row, &parsed.voices, &mut remaining_shapes)?;
                let target_m = row.scalar()?;
                let scale_m = row.scalar()?;
                let weight = row.scalar()?;
                row.finish()?;
                targets.push(DisplacementTarget { attachment, target_m, scale_m, weight });
            }
            cases.push(DesignLoadCase { name, loads, targets });
        }
        let mut row = reader.row("variables")?;
        let count = row.count(budget.max_variables)?;
        row.finish()?;
        let mut variables = Vec::with_capacity(count);
        let mut bindings_left = budget.max_bindings;
        for _ in 0..count {
            poll(gate)?;
            let mut row = reader.row("variable")?;
            let name = name(&mut row)?;
            let reference = row.scalar()?;
            let scale = row.scalar()?;
            let minimum = row.scalar()?;
            let maximum = row.scalar()?;
            let fields_count = row.count(bindings_left)?;
            row.finish()?;
            bindings_left -= fields_count;
            let mut fields = Vec::with_capacity(fields_count);
            for _ in 0..fields_count {
                let mut row = reader.row("bind")?;
                let field = match row.word()? {
                    "spring-stiffness" => DesignField::SpringStiffness(row.parse()?),
                    "spring-rest" => DesignField::SpringRest(row.parse()?),
                    "contact-stiffness" => DesignField::ContactStiffness(row.parse()?),
                    "contact-gap" => DesignField::ContactGap(row.parse()?),
                    "contact-weight" => DesignField::ContactWeight(row.parse()?),
                    "actuator-force" => DesignField::ActuatorForce { case: row.parse()?, actuator: row.parse()? },
                    _ => return Err(input(row.line, "unknown static design field; no ignored bindings").into()),
                };
                row.finish()?;
                fields.push(field);
            }
            variables.push(DesignVariable { name, reference, scale, minimum, maximum, fields });
        }
        let (constraints, maximum_constraints) = if reader.lines.clone().next().is_some() {
            read_constraints(&mut reader, &parsed.voices, &mut remaining_shapes, gate)?
        } else { (Vec::new(), 0) };
        if reader.lines.next().is_some() { return Err(input(reader.line + 1, "unexpected trailing design record").into()); }
        let contacts = match parsed.multiple {
            Some((contacts, _)) => contacts,
            None => parsed.contact.into_iter().collect(),
        };
        let models = parsed.voices.into_iter().map(|voice| voice.model).collect();
        let problem = EquilibriumDesign::new(models, connections, contacts, cases, variables, budget, gate)?
            .with_constraints(constraints, maximum_constraints, gate)?;
        poll(gate)?;
        Ok(Self { problem, model_info: parsed.info,
            design_hash: hash_domain(EQUILIBRIUM_DESIGN_HASH_DOMAIN, design_bytes) })
    }

    /// Borrow the actual model-independent physics objective for an optimizer.
    #[must_use]
    pub fn problem(&self) -> &EquilibriumDesign { &self.problem }
    /// Move it to a caller-owned study without rereading either input.
    #[must_use]
    pub fn into_problem(self) -> EquilibriumDesign { self.problem }
    /// Source model description and its original versioned input hash.
    #[must_use]
    pub const fn model_info(&self) -> ModalPerformanceInfo { self.model_info }
    /// Exact design bytes, not a claim about the origin of observations.
    #[must_use]
    pub const fn design_hash(&self) -> ContentHash { self.design_hash }
}

fn port(row: &mut Row<'_>, voices: &[ModalForceVoice], remaining: &mut usize)
    -> Result<ModalAttachment, ModalPerformanceError>
{
    let component: usize = row.parse()?;
    let index: usize = row.parse()?;
    let shapes = voices.get(component).and_then(|v| v.port_shapes.get(index))
        .ok_or_else(|| input(row.line, "load/target must name an existing model component and port"))?;
    *remaining = remaining.checked_sub(shapes.len())
        .ok_or_else(|| input(row.line, "total copied load/target shape coefficients exceed 65536"))?;
    Ok(ModalAttachment { component, shapes: shapes.clone() })
}
fn name(row: &mut Row<'_>) -> Result<String, ModalPerformanceError> {
    let text = row.word()?;
    if text.len() > 128 { return Err(input(row.line, "case/variable name exceeds 128 bytes")); }
    Ok(text.to_owned())
}
fn model_error(what: &'static str) -> DesignFileError { DesignFileError::Model(input(1, what)) }
fn poll(gate: &CancelGate) -> Result<(), DesignFileError> {
    if gate.is_requested() { Err(DesignError::Cancelled.into()) } else { Ok(()) }
}

// Additive, strict section: no alternate mechanical parser or implicit units.
fn read_constraints(reader: &mut Reader<'_>, voices: &[ModalForceVoice], remaining_shapes: &mut usize,
    gate: &CancelGate) -> Result<(Vec<ResponseConstraint>, usize), DesignFileError>
{
    let mut row = reader.row("constraint_limits")?;
    let maximum = row.count(64)?;
    row.finish()?;
    let mut row = reader.row("constraints")?;
    let count = row.count(maximum)?;
    row.finish()?;
    let mut constraints = Vec::with_capacity(count);
    for _ in 0..count {
        poll(gate)?;
        let mut row = reader.row("constraint")?;
        let name = name(&mut row)?;
        let case = row.parse()?;
        let quantity = match row.word()? {
            "displacement" => ResponseQuantity::Displacement(port(&mut row, voices, remaining_shapes)?),
            "spring-force" => ResponseQuantity::SpringForce(row.parse()?),
            "contact-force" => ResponseQuantity::ContactForce(row.parse()?),
            "contact-penetration" => ResponseQuantity::ContactPenetration(row.parse()?),
            _ => return Err(input(row.line, "unknown physical constraint response").into()),
        };
        let sense = match row.word()? {
            "at-most" => ConstraintSense::AtMost,
            "at-least" => ConstraintSense::AtLeast,
            "equal" => ConstraintSense::Equal,
            _ => return Err(input(row.line, "expected at-most, at-least or equal").into()),
        };
        let bound = row.scalar()?;
        let scale = row.scalar()?;
        row.finish()?;
        constraints.push(ResponseConstraint { name, case, quantity, sense, bound, scale });
    }
    Ok((constraints, maximum))
}

#[cfg(test)]
#[path = "design/constraints_tests.rs"]
mod constraints_tests;
