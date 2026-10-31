//! Memoryless scalar conduction at the fields solved by the existing circuit.
//!
//! Material curves supply conductivity or resistivity; this module supplies no
//! second current law or circuit integrator. Each bounded constitutive trial
//! restarts from the same accepted charge/flux. Only the original electrical
//! solver advances time and accounts for source work and Joule heat.

use fs_matdb::MaterialCard;
use fs_material::conductor::{ConductorError, ResolvedConductor};
use fs_qty::{Dims, QuantitySpec};

use super::contacts::{ContactHeat, SeriesElectricalContact};
use super::{
    MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitError, MaterialCircuitSolve,
    MaterialCircuitStep, MaterialResistor, ResistorHeat, poll,
};

pub mod contact;
use contact::{ContactConductivitySample, ContactResponse};

/// Work cap on complete circuit/source evaluations in one proposed time step.
pub const MAX_CONDUCTIVITY_ITERATIONS: usize = 128;

/// Physical coordinate to update in the original material query.
/// Signed coordinates follow the circuit branch's first-to-second orientation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConductionCoordinate {
    /// Bulk voltage divided by the original conductor length [V/m].
    ElectricField,
    /// Absolute bulk voltage divided by conductor length [V/m].
    ElectricFieldMagnitude,
    /// Branch current divided by the original cross-sectional area [A/m2].
    CurrentDensity,
    /// Absolute branch current divided by cross-sectional area [A/m2].
    CurrentDensityMagnitude,
}
impl ConductionCoordinate {
    pub(crate) fn dimensions(self) -> Dims {
        match self {
            Self::ElectricField | Self::ElectricFieldMagnitude => Dims([1, 1, -3, 0, -1, 0]),
            Self::CurrentDensity | Self::CurrentDensityMagnitude => Dims([-2, 0, 0, 0, 1, 0]),
        }
    }
    fn value(self, point: ConductionPoint) -> f64 {
        match self {
            Self::ElectricField => point.electric_field_v_m,
            Self::ElectricFieldMagnitude => point.electric_field_v_m.abs(),
            Self::CurrentDensity => point.current_density_a_m2,
            Self::CurrentDensityMagnitude => point.current_density_a_m2.abs(),
        }
    }
}

/// Bind one explicitly typed axis already present in a conductor's query.
#[derive(Clone, Copy, Debug)]
pub struct ConductionAxis<'a> {
    /// Original material coordinate name; no conventional spelling is inferred.
    pub name: &'a str,
    /// Signed or magnitude interpretation declared by the source owner.
    pub coordinate: ConductionCoordinate,
}

/// Same immutable material card for an original circuit resistor branch.
/// Every resistor needs one binding; field/current axes cannot remain frozen.
pub struct ConductorResponse<'a> {
    /// Original graph branch index, not its resistor-only position.
    pub branch: usize,
    /// Must match the complete source bundle already retained by that branch.
    pub card: &'a MaterialCard,
    /// At most one electric-field and one current-density coordinate. Empty for
    /// a source with no electrical operating coordinates (e.g. a thermal law).
    pub axes: &'a [ConductionAxis<'a>],
}

/// Bounded constitutive fixed-point controls, independent of native DAE work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductivitySolve {
    /// Complete native circuit trials in `1..=MAX_CONDUCTIVITY_ITERATIONS`.
    pub maximum_iterations: usize,
    /// Absolute applied/evaluated resistance tolerance [ohm], nonnegative.
    pub absolute_resistance_tolerance_ohm: f64,
    /// Relative resistance tolerance in `[0, 1)`; at least one tolerance is > 0.
    pub relative_resistance_tolerance: f64,
}
impl ConductivitySolve {
    pub(crate) fn validate(self) -> Result<(), MaterialCircuitError> {
        if !(1..=MAX_CONDUCTIVITY_ITERATIONS).contains(&self.maximum_iterations)
            || !self.absolute_resistance_tolerance_ohm.is_finite()
            || self.absolute_resistance_tolerance_ohm < 0.0
            || !self.relative_resistance_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_resistance_tolerance)
            || (self.absolute_resistance_tolerance_ohm == 0.0
                && self.relative_resistance_tolerance == 0.0)
        {
            return Err(invalid("invalid bounded conductivity solve controls"));
        }
        Ok(())
    }
}

/// Actual bulk branch quantities from the native algebraic midpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionPoint {
    /// Original circuit branch index.
    pub branch: usize,
    /// Signed first-to-second bulk voltage [V].
    pub voltage_v: f64,
    /// Signed first-to-second current [A].
    pub current_a: f64,
    /// Signed `voltage_v / length_m` [V/m].
    pub electric_field_v_m: f64,
    /// Signed `current_a / area_m2` [A/m2].
    pub current_density_a_m2: f64,
}

/// Original complete bundle re-evaluated at the actual solved branch fields.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductivitySample {
    /// Actual coordinates supplied to the original material query.
    pub point: ConductionPoint,
    /// Same geometry, source basis, card and selected claims at those coordinates.
    pub conductor: ResolvedConductor,
    /// Signed evaluated resistance minus the resistance used by the native step.
    pub resistance_defect_ohm: f64,
    /// Signed `V - I R_evaluated` [V].
    pub voltage_defect_v: f64,
    /// Signed `V I - I^2 R_evaluated` [W].
    pub power_defect_w: f64,
}

/// Native evolution plus the separately admitted constitutive residual.
pub struct MaterialConductivityStep {
    /// Authoritative native storage, source work and Joule heat. Its retained
    /// materials are the coefficients actually used to compute this step.
    pub electrical: MaterialCircuitStep,
    /// Source responses at that step's solved fields, in original branch order.
    /// Their difference from applied coefficients is bounded, never hidden.
    pub evaluated: Vec<ConductivitySample>,
    /// Complete interface bundles evaluated at each actual contact drop/current.
    /// Empty when contacts have no electrical response binding.
    pub evaluated_contacts: Vec<ContactConductivitySample>,
    /// All circuit trials, including discarded constitutive iterates.
    pub iterations: usize,
    /// Operating-field bundle refreshes; excludes the one seed-admission replay
    /// per branch before the first circuit trial.
    pub source_evaluations: usize,
}

/// Uncommitted electrical endpoint; no accepted state or thermal store changes.
pub struct MaterialConductivityCandidate {
    circuit: MaterialCircuit,
    step: MaterialConductivityStep,
}
impl MaterialConductivityCandidate {
    /// Proposed charge/flux state with the applied source coefficients retained.
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    /// Native energy account and actual-field source evaluations.
    pub const fn step(&self) -> &MaterialConductivityStep {
        &self.step
    }
    /// Publish alongside any thermal owner receiving the native branch heats.
    pub fn into_parts(self) -> (MaterialCircuit, MaterialConductivityStep) {
        (self.circuit, self.step)
    }
}

fn invalid(message: &'static str) -> MaterialCircuitError {
    MaterialCircuitError::Invalid(message)
}

fn admit<'a, 'b>(
    circuit: &MaterialCircuit,
    seeds: &[MaterialResistor],
    sources: &'a [ConductorResponse<'b>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<&'a ConductorResponse<'b>>, MaterialCircuitError> {
    if circuit.resistors().is_empty()
        || sources.len() != circuit.resistors().len()
        || seeds.len() != circuit.resistors().len()
    {
        return Err(invalid(
            "every original resistor needs one conductivity response",
        ));
    }
    let mut ordered: Vec<_> = sources.iter().collect();
    ordered.sort_by_key(|binding| binding.branch);
    for ((binding, original), seed) in ordered.iter().zip(circuit.resistors()).zip(seeds) {
        poll(cancelled)?;
        let material = seed.conductor.material();
        if binding.branch != original.branch
            || binding.card.content_hash() != material.card_identity()
            || seed.branch != original.branch
            || seed.conductor.source_basis() != original.conductor.source_basis()
            || seed.conductor.source_key() != original.conductor.source_key()
            || seed.conductor.length_m().to_bits() != original.conductor.length_m().to_bits()
            || seed.conductor.area_m2().to_bits() != original.conductor.area_m2().to_bits()
            || binding.axes.len() > 2
        {
            return Err(invalid(
                "conductivity response must retain branch, card and at most two field axes",
            ));
        }
        let point = material.query();
        // A thermal owner may update the context, but not replace the original
        // source selection, requirements, card, geometry or constitutive basis.
        let replay = original
            .conductor
            .material()
            .at_query_point(binding.card, &point)
            .map_err(ConductorError::MaterialState)?;
        // The initial source may have been admitted with SingleClaimOnly or a
        // preference policy; refresh records an explicit pin. Those honest
        // policy receipts differ, while the selected claims and full property
        // requirements must remain identical.
        if replay.properties().len() != material.properties().len()
            || replay
                .properties()
                .iter()
                .zip(material.properties())
                .any(|(expected, seed)| {
                    expected.requirement() != seed.requirement()
                        || expected.answer().receipt.selected != seed.answer().receipt.selected
                        || expected.value_si().to_bits() != seed.value_si().to_bits()
                })
        {
            return Err(invalid(
                "refreshed electrical seed changed its original source bundle",
            ));
        }
        for (index, axis) in binding.axes.iter().enumerate() {
            if !point.axes().contains_key(axis.name)
                || point.axis_quantities().get(axis.name)
                    != Some(&QuantitySpec::dimensional(axis.coordinate.dimensions()))
                || binding.axes[..index].iter().any(|other| {
                    other.name == axis.name
                        || other.coordinate.dimensions() == axis.coordinate.dimensions()
                })
            {
                return Err(invalid(
                    "conductivity axes require unique existing SI field/current-density coordinates",
                ));
            }
        }
        for (name, quantity) in point.axis_quantities() {
            if [
                ConductionCoordinate::ElectricField.dimensions(),
                ConductionCoordinate::CurrentDensity.dimensions(),
            ]
            .contains(&quantity.dims())
                && !binding.axes.iter().any(|axis| axis.name == name.as_str())
            {
                return Err(invalid(
                    "an electrical operating coordinate cannot remain frozen",
                ));
            }
        }
        // A legacy untyped query must not hide the independent field variable
        // recorded by its original constitutive response.
        for property in material.properties() {
            let claim = binding
                .card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("conductivity bundle lost its selected claim"))?;
            for (name, dims) in claim.value.axes() {
                if [
                    ConductionCoordinate::ElectricField.dimensions(),
                    ConductionCoordinate::CurrentDensity.dimensions(),
                ]
                .contains(&dims)
                    && !binding.axes.iter().any(|axis| axis.name == name)
                {
                    return Err(invalid(
                        "a source law's electrical coordinate cannot remain frozen",
                    ));
                }
            }
        }
    }
    Ok(ordered)
}

// The outer owner may refresh thermal context, but contact identity, geometry,
// coordinate meanings and every selected property remain the original bundle.
fn admit_contacts(
    circuit: &MaterialCircuit,
    contacts: &[SeriesElectricalContact],
    electrical_response: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaterialCircuitError> {
    if contacts.len() != circuit.contacts().len() {
        return Err(invalid("every original contact needs one refreshed source"));
    }
    for contact in contacts {
        poll(cancelled)?;
        let original = circuit.contacts().iter().find(|old| old.id() == contact.id())
            .ok_or_else(|| invalid("unbound conductivity contact"))?;
        let before = original.source().resolved();
        let after = contact.source().resolved();
        if contact.branch() != original.branch()
            || contact.area_m2().to_bits() != original.area_m2().to_bits()
            || before.card_identity() != after.card_identity()
            || before.surface_materials() != after.surface_materials()
            || before.axis_quantities() != after.axis_quantities()
            || before.query_point().len() != after.query_point().len()
            || before.query_point().iter().zip(after.query_point())
                .any(|((a, _), (b, _))| a != b)
            || before.properties().len() != after.properties().len()
            || before.properties().iter().zip(after.properties()).any(|(a, b)| {
                a.requirement() != b.requirement()
                    || a.answer().receipt.selected != b.answer().receipt.selected
            })
        {
            return Err(invalid("conductivity contact refresh changed its original source or geometry"));
        }
        if !electrical_response && (after.properties().iter().any(|property|
            property.source_axes().iter().any(|(_, dims)| contact::electrical(*dims)))
            || after.axis_quantities().values().any(|quantity| contact::electrical(quantity.dims()))) {
            return Err(invalid("contact electrical operating coordinates require their own nonlinear closure"));
        }
    }
    Ok(())
}

fn refresh(
    original: &ResolvedConductor,
    binding: &ConductorResponse<'_>,
    point: ConductionPoint,
) -> Result<ResolvedConductor, MaterialCircuitError> {
    let mut query = original.material().query();
    for axis in binding.axes {
        query = query
            .with(axis.name, axis.coordinate.value(point))
            .map_err(|_| invalid("solved electrical coordinate is unrepresentable"))?;
    }
    let material = original
        .material()
        .at_query_point(binding.card, &query)
        .map_err(ConductorError::MaterialState)?;
    Ok(original.with_material(material)?)
}

// The constitutive iteration is shared by the original fixed and moving-field
// descriptors. Each driver returns its OWN native candidate and loss account;
// this seam cannot manufacture a circuit, force, heat or accepted history.
pub(super) trait ConductivityDriver {
    fn resistor_heats(&self) -> &[ResistorHeat];
    fn contact_heats(&self) -> &[ContactHeat] { &[] }
}
impl ConductivityDriver for MaterialCircuitCandidate {
    fn resistor_heats(&self) -> &[ResistorHeat] { &self.step().resistors }
    fn contact_heats(&self) -> &[ContactHeat] { &self.step().contacts }
}
pub(super) struct ConductivityTrial<D> {
    pub candidate: D,
    pub evaluated: Vec<ConductivitySample>,
    pub evaluated_contacts: Vec<ContactConductivitySample>,
    pub iterations: usize,
    pub source_evaluations: usize,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_native<D: ConductivityDriver, C: FnMut() -> bool>(
    circuit: &MaterialCircuit,
    mut seeds: Vec<MaterialResistor>,
    sources: &[ConductorResponse<'_>],
    mut contacts: Vec<SeriesElectricalContact>,
    contact_sources: Option<&[ContactResponse<'_>]>,
    controls: MaterialCircuitSolve,
    solve: ConductivitySolve,
    mut native: impl FnMut(Vec<MaterialResistor>, Vec<SeriesElectricalContact>, &mut C)
        -> Result<D, MaterialCircuitError>,
    mut cancelled: C,
) -> Result<ConductivityTrial<D>, MaterialCircuitError> {
    poll(&mut cancelled)?;
    solve.validate()?;
    if seeds.len() != circuit.resistors().len() {
        return Err(invalid(
            "every original resistor needs one refreshed material seed",
        ));
    }
    seeds.sort_by_key(|seed| seed.branch);
    let bindings = admit(circuit, &seeds, sources, &mut cancelled)?;
    super::contacts::admit(&mut contacts)?;
    admit_contacts(circuit, &contacts, contact_sources.is_some(), &mut cancelled)?;
    let contact_bindings = contact::admit(circuit, &contacts, contact_sources, &mut cancelled)?;
    let mut applied = seeds.clone();
    let mut applied_contacts = contacts.clone();
    let mut maximum_resistance_defect = 0.0_f64;
    let mut maximum_power_defect = 0.0_f64;
    for iterations in 1..=solve.maximum_iterations {
        poll(&mut cancelled)?;
        let candidate = native(applied, applied_contacts, &mut cancelled)?;
        if candidate.resistor_heats().len() != seeds.len()
            || candidate.contact_heats().len() != contacts.len() {
            return Err(invalid("native conductivity trial lost a bulk or contact loss"));
        }
        let mut evaluated = Vec::with_capacity(bindings.len());
        let mut next = Vec::with_capacity(bindings.len());
        let mut converged = true;
        maximum_resistance_defect = 0.0;
        maximum_power_defect = 0.0;
        for ((heat, original), binding) in
            candidate.resistor_heats().iter().zip(&seeds).zip(&bindings)
        {
            poll(&mut cancelled)?;
            let point = ConductionPoint {
                branch: original.branch,
                voltage_v: heat.voltage_v,
                current_a: heat.current_a,
                electric_field_v_m: heat.voltage_v / original.conductor.length_m(),
                current_density_a_m2: heat.current_a / original.conductor.area_m2(),
            };
            if !point.electric_field_v_m.is_finite()
                || !point.current_density_a_m2.is_finite()
                || point.voltage_v != 0.0 && point.electric_field_v_m == 0.0
                || point.current_a != 0.0 && point.current_density_a_m2 == 0.0
            {
                return Err(invalid(
                    "solved conductor field/current density is unrepresentable",
                ));
            }
            let conductor = refresh(&original.conductor, binding, point)?;
            let used = heat.material.conductor.resistance_ohm();
            let actual = conductor.resistance_ohm();
            let resistance_defect_ohm = actual - used;
            let voltage_defect_v = point.voltage_v - point.current_a * actual;
            let power_defect_w =
                point.voltage_v * point.current_a - conductor.joule_power_w(point.current_a)?;
            let tolerance = solve.absolute_resistance_tolerance_ohm
                + solve.relative_resistance_tolerance * actual.max(used);
            if [resistance_defect_ohm, voltage_defect_v, power_defect_w, tolerance]
                .iter().any(|value| !value.is_finite())
            {
                return Err(invalid("constitutive residual or tolerance is unrepresentable"));
            }
            maximum_resistance_defect = maximum_resistance_defect.max(resistance_defect_ohm.abs());
            maximum_power_defect = maximum_power_defect.max(power_defect_w.abs());
            converged &= resistance_defect_ohm.abs() <= tolerance
                && power_defect_w.abs() <= controls.power_tolerance_w;
            next.push(MaterialResistor { branch: original.branch, conductor: conductor.clone() });
            evaluated.push(ConductivitySample { point, conductor, resistance_defect_ohm,
                voltage_defect_v, power_defect_w });
        }
        let mut evaluated_contacts = Vec::with_capacity(contact_bindings.len());
        let mut next_contacts = if contact_bindings.is_empty() { contacts.clone() }
            else { Vec::with_capacity(contacts.len()) };
        for ((heat, original), binding) in candidate.contact_heats().iter()
            .zip(&contacts).zip(&contact_bindings) {
            poll(&mut cancelled)?;
            let sample = contact::evaluate(original, binding, heat)?;
            let actual = sample.contact.resistance_ohm();
            let used = heat.contact.resistance_ohm();
            let tolerance = solve.absolute_resistance_tolerance_ohm
                + solve.relative_resistance_tolerance * actual.max(used);
            if !tolerance.is_finite() {
                return Err(invalid("contact constitutive tolerance is unrepresentable"));
            }
            maximum_resistance_defect = maximum_resistance_defect.max(sample.resistance_defect_ohm.abs());
            maximum_power_defect = maximum_power_defect.max(sample.power_defect_w.abs());
            converged &= sample.resistance_defect_ohm.abs() <= tolerance
                && sample.power_defect_w.abs() <= controls.power_tolerance_w;
            next_contacts.push(sample.contact.clone());
            evaluated_contacts.push(sample);
        }
        poll(&mut cancelled)?;
        if converged {
            return Ok(ConductivityTrial { candidate, evaluated, evaluated_contacts, iterations,
                source_evaluations: iterations * (bindings.len() + contact_bindings.len()) });
        }
        applied = next;
        applied_contacts = next_contacts;
    }
    Err(MaterialCircuitError::ConductivityConvergence {
        iterations: solve.maximum_iterations,
        resistance_defect_ohm: maximum_resistance_defect,
        power_defect_w: maximum_power_defect,
    })
}

impl MaterialCircuit {
    /// Resolve conductivity/resistivity from the actual solved electric field
    /// and/or current density, retaining the original complete source bundle.
    /// The supplied axes must already exist with their SI quantity descriptors;
    /// temperature, pressure, composition and all other coordinates stay fixed.
    /// Signed versus magnitude conventions are explicit on each axis binding.
    ///
    /// Each trial uses the SAME original charge/flux and native Kirchhoff step.
    /// At its algebraic midpoint, `E=V/L` and `J=I/A` refresh the original pinned
    /// claims. Only those material responses supply the next resistance; there
    /// is no fit, damping, clipping, artificial resistor or source fallback.
    ///
    /// Acceptance requires every resistance residual to meet `solve` AND every
    /// constitutive power residual to meet `controls.power_tolerance_w`. Native
    /// Joule heat remains counted once; evaluated material heat is a residual,
    /// never an extra source. Both applied and actual-field evaluated sources
    /// are retained, so finite nonlinear tolerance is inspectable.
    ///
    /// This is an instantaneous scalar, fixed-geometry, fixed-C/L midpoint
    /// model. It does not add carrier/charge transport, dielectric memory, arcs,
    /// hysteresis, Hall response or a convergence/time-accuracy guarantee.
    /// Contact-bound circuits require their own jointly refreshed contact law.
    /// Source/domain refusal, nonconvergence and cancellation preserve `self`.
    pub fn conductivity_candidate(
        &self,
        sources: &[ConductorResponse<'_>],
        inputs: &[f64],
        dt: f64,
        controls: MaterialCircuitSolve,
        solve: ConductivitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialConductivityCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        self.conductivity_with_materials_candidate(
            self.resistors().to_vec(),
            sources,
            inputs,
            dt,
            controls,
            solve,
            cancelled,
        )
    }

    /// The same nonlinear electrical solve at explicitly refreshed material
    /// context, for an outer thermal/pressure/material-state transaction.
    ///
    /// Every seed must be an exact-source refresh of its original complete
    /// bundle, with unchanged branch, geometry and conductivity/resistivity
    /// basis. For example, a thermal owner can provide its candidate T while
    /// this owner closes E/J at that SAME T. A trial never spends the original
    /// electrical storage, and field coordinates cannot be held at the seed.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_materials_candidate(
        &self,
        seeds: Vec<MaterialResistor>,
        sources: &[ConductorResponse<'_>],
        inputs: &[f64],
        dt: f64,
        controls: MaterialCircuitSolve,
        solve: ConductivitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialConductivityCandidate, MaterialCircuitError> {
        let trial = solve_native(self, seeds, sources, Vec::new(), None, controls, solve,
            |applied, _, check| self.step_candidate(applied, inputs, dt, controls, check), cancelled)?;
        let (circuit, electrical) = trial.candidate.into_parts();
        Ok(MaterialConductivityCandidate { circuit, step: MaterialConductivityStep {
            electrical, evaluated: trial.evaluated, evaluated_contacts: trial.evaluated_contacts,
            iterations: trial.iterations,
            source_evaluations: trial.source_evaluations } })
    }

    /// Close nonlinear BULK conduction with explicit same-context series contacts.
    ///
    /// The original descriptor solves the sum of bulk and contact drops. Its
    /// bulk-only voltage feeds `E=V_bulk/L`; contact loss is retained separately
    /// and never enters the bulk constitutive power residual or heat twice.
    /// An outer thermal owner supplies fresh contact and bulk sources at each
    /// proposed temperature/phase state. Each inner trial restarts from the
    /// same accepted charge/flux and holds those contact snapshots fixed.
    ///
    /// Complete contact source selections, ordered surfaces, coordinate meanings,
    /// branch identities and areas must remain unchanged. Electrical-field,
    /// voltage/current dependent contacts need their own nonlinear closure and
    /// refuse; no interface storage, arcing or geometry evolution is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_contacts_candidate(
        &self,
        seeds: Vec<MaterialResistor>,
        contacts: Vec<SeriesElectricalContact>,
        sources: &[ConductorResponse<'_>],
        inputs: &[f64],
        dt: f64,
        controls: MaterialCircuitSolve,
        solve: ConductivitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialConductivityCandidate, MaterialCircuitError> {
        self.conductivity_with_contact_response_candidate(seeds, contacts, sources,
            None, inputs, dt, controls, solve, cancelled)
    }

    /// Jointly close bulk E/J and every contact's own voltage/current response.
    ///
    /// Each interface binding declares its ordered-surface direction and typed
    /// voltage, current or current-density axes. Complete pinned interface
    /// bundles are refreshed at the SAME native midpoint as bulk materials.
    /// Both coefficient and constitutive-power residuals must converge. Native
    /// contact heat is deposited once; evaluated-source defects are not heat.
    ///
    /// `None` holds electrically independent contacts fixed. `Some` requires
    /// exactly one original-card binding per contact, including fixed laws.
    /// The zero-thickness contact owns no electric field, capacitance, carrier
    /// storage, breakdown, switching history or geometrical evolution. Strong
    /// nonlinearity can refuse the bounded fixed-point iteration atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn conductivity_with_contact_response_candidate(
        &self,
        seeds: Vec<MaterialResistor>,
        contacts: Vec<SeriesElectricalContact>,
        sources: &[ConductorResponse<'_>],
        contact_sources: Option<&[ContactResponse<'_>]>,
        inputs: &[f64],
        dt: f64,
        controls: MaterialCircuitSolve,
        solve: ConductivitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialConductivityCandidate, MaterialCircuitError> {
        let trial = solve_native(self, seeds, sources, contacts, contact_sources, controls, solve,
            |applied, interfaces, check| self.step_with_contacts_candidate(
                applied, interfaces, inputs, dt, controls, check), cancelled)?;
        let (circuit, electrical) = trial.candidate.into_parts();
        Ok(MaterialConductivityCandidate { circuit, step: MaterialConductivityStep {
            electrical, evaluated: trial.evaluated, evaluated_contacts: trial.evaluated_contacts,
            iterations: trial.iterations,
            source_evaluations: trial.source_evaluations } })
    }
}
