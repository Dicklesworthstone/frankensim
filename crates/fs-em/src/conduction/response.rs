//! Scalar source response at the fields recovered by the original spatial solve.
//!
//! Each bounded iteration assembles the SAME FEEC conductivity operator, solves
//! its original voltage/current traces, and refreshes complete pinned sources at
//! actual |E| and |J|. Applied/evaluated conductivity, current and power must all
//! agree. Optional relaxation requeries original sources at damped E/J
//! coordinates; it never blends conductivities or adds iteration error as heat.
//! This Picard closure need not converge for every positive source curve and
//! proves neither uniqueness nor stability.

use super::*;
use fs_material::conductor::{ConductorError, response::ScalarConductorResponse};
use fs_material::electrical_interface::response::ScalarElectricalInterfaceResponse;
use std::sync::Arc;

/// Original local contact source laws and uniform-trace response gates.
pub mod contact;
mod power;
use contact::{
    ContactResponseControls, ContactResponseError, ContactResponseSample, ContactResponseSources,
};

pub const MAX_RESPONSE_ITERATIONS: usize = 128;

/// Safeguarded adaptive relaxation of original electrical source coordinates.
/// Every trial weight lies in `[minimum_weight, maximum_weight]` within
/// `(0, 1]`. Trials requery the source at interpolated E/J or contact coordinates;
/// they never interpolate conductivity, resistance, temperature or source pins.
/// Positive resistance alone does not guarantee convergence or uniqueness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveConductionRelaxation {
    /// Weight used before a secant history is available.
    pub initial_weight: f64,
    /// Strictly positive lower bound for the numerical weight.
    pub minimum_weight: f64,
    /// Upper bound, at most one, for convex source-coordinate trials.
    pub maximum_weight: f64,
}
impl AdaptiveConductionRelaxation {
    pub(crate) fn prepare(
        self,
    ) -> Result<fs_ports::BoundedAitkenRelaxation, fs_ports::AitkenRelaxationError> {
        fs_ports::BoundedAitkenRelaxation::new(
            self.initial_weight,
            self.minimum_weight,
            self.maximum_weight,
        )
    }
}

/// Constitutive gates in addition to the native algebraic/current/work gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionResponseControls {
    pub maximum_iterations: usize,
    pub absolute_conductivity_tolerance_s_m: f64,
    pub relative_conductivity_tolerance: f64,
    pub current_density_tolerance_a_m2: f64,
}
impl ConductionResponseControls {
    fn validate(self) -> Result<(), ConductionResponseError> {
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&self.maximum_iterations)
            || !self.absolute_conductivity_tolerance_s_m.is_finite()
            || self.absolute_conductivity_tolerance_s_m < 0.0
            || !self.relative_conductivity_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_conductivity_tolerance)
            || (self.absolute_conductivity_tolerance_s_m == 0.0
                && self.relative_conductivity_tolerance == 0.0)
            || !positive(self.current_density_tolerance_a_m2)
        {
            return Err(ConductionResponseError::Invalid(
                "invalid bounded spatial response controls",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum ConductionResponseError {
    Field(ConductionError),
    Material(ConductorError),
    Contact(ContactResponseError),
    Invalid(&'static str),
    NotConverged {
        iterations: usize,
        maximum_conductivity_defect_s_m: f64,
        maximum_current_density_defect_a_m2: f64,
        absolute_power_defect_w: f64,
        maximum_contact_resistance_defect_ohm_m2: f64,
        maximum_contact_voltage_defect_v: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ConductionResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial conductor response: {self:?}")
    }
}
impl std::error::Error for ConductionResponseError {}
impl From<ConductionError> for ConductionResponseError {
    fn from(error: ConductionError) -> Self {
        match error {
            ConductionError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<ConductorError> for ConductionResponseError {
    fn from(error: ConductorError) -> Self {
        Self::Material(error)
    }
}
impl From<ContactResponseError> for ConductionResponseError {
    fn from(error: ContactResponseError) -> Self {
        Self::Contact(error)
    }
}

/// Complete continuing source and numerical policy for a frozen field handoff.
/// Bulk laws, contact laws and relaxation travel together, so a coupled owner
/// cannot silently freeze a nonlinear interface or lose its selected solver policy.
#[derive(Clone, Debug)]
pub struct ConductionResponseSources {
    pub(crate) bulk: Arc<[ScalarConductorResponse]>,
    pub(crate) relaxation: f64,
    pub(crate) adaptive_relaxation: Option<AdaptiveConductionRelaxation>,
    pub(crate) contacts: Option<ContactResponseSources>,
}
impl ConductionResponseSources {
    pub fn bulk(&self) -> &[ScalarConductorResponse] {
        &self.bulk
    }
    pub fn relaxation(&self) -> f64 {
        self.relaxation
    }
    pub fn adaptive_relaxation(&self) -> Option<AdaptiveConductionRelaxation> {
        self.adaptive_relaxation
    }
    pub fn contacts(&self) -> Option<&ContactResponseSources> {
        self.contacts.as_ref()
    }
}
impl From<Vec<ScalarConductorResponse>> for ConductionResponseSources {
    fn from(bulk: Vec<ScalarConductorResponse>) -> Self {
        Self {
            bulk: bulk.into(),
            relaxation: 1.0,
            adaptive_relaxation: None,
            contacts: None,
        }
    }
}
impl From<Arc<[ScalarConductorResponse]>> for ConductionResponseSources {
    fn from(bulk: Arc<[ScalarConductorResponse]>) -> Self {
        Self {
            bulk,
            relaxation: 1.0,
            adaptive_relaxation: None,
            contacts: None,
        }
    }
}

/// Same geometry and original sources; no independent circuit or thermal state.
/// A frozen linearization is private outside an explicit coupled-solve callback.
/// Its actual-amplitude solve is repeated inside the source closure; accepted
/// nonlinear fields cannot be scaled to a different excitation or obtained by
/// silently freezing an operating point.
#[derive(Clone)]
pub struct PlanarConductionResponse {
    native: PlanarConduction,
    sources: Arc<[ScalarConductorResponse]>,
    relaxation: f64,
    adaptive_relaxation: Option<AdaptiveConductionRelaxation>,
    contact_sources: Option<ContactResponseSources>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConductionResponseSample {
    pub triangle: usize,
    /// Complete original material evaluated at the actual reported E and J.
    pub conductor: ResolvedConductorState,
    /// Evaluated minus applied conductivity [S/m].
    pub conductivity_defect_s_m: f64,
    /// Native J minus the source-evaluated J(E) [A/m2], in mesh axes.
    pub current_density_defect_a_m2: [f64; 2],
    /// Native local E.J minus source-evaluated power over this cell [W].
    pub power_defect_w: f64,
}

/// Original native field observations and separate constitutive residuals.
/// Applied sources remain in `candidate.problem().conductors()`; evaluated
/// sources below agree to the declared gates, without post-solve substitution.
#[derive(Clone, Debug)]
pub struct ConductionResponseReport {
    /// Absolute reporting view; physical fields use the retained differences.
    pub voltage_v: Vec<f64>,
    pub relative_voltage_v: Vec<f64>,
    pub reference_voltage_v: Vec<f64>,
    pub electric_field_v_m: Vec<[f64; 2]>,
    pub current_density_a_m2: Vec<[f64; 2]>,
    pub cell_power_w: Vec<f64>,
    /// Native surface losses in canonical contact order, excluded from bulk
    /// `cell_power_w` and included exactly once in `power_w`.
    pub contact_response: Vec<ConductionContactResponse>,
    pub power_w: f64,
    pub nodal_boundary_current_a: Vec<Option<f64>>,
    pub nodal_applied_current_a: Vec<Option<f64>>,
    pub boundary_current_a: f64,
    pub maximum_free_current_residual_a: f64,
    pub power_balance_residual_w: f64,
    pub linear: Option<SolveReport>,
    pub evaluated: Vec<ConductionResponseSample>,
    /// Actual source queries and constitutive residuals for uniform contacts.
    pub evaluated_contacts: Vec<ContactResponseSample>,
    /// Every native spatial solve, including discarded constitutive trials.
    pub field_solves: usize,
    /// Complete actual-field AND relaxed-coordinate source queries, excluding
    /// initial binding admission. Full steps and fixed sources add no queries.
    pub source_evaluations: usize,
}

pub struct ConductionResponseCandidate {
    problem: PlanarConductionResponse,
    report: ConductionResponseReport,
}
impl ConductionResponseCandidate {
    pub fn problem(&self) -> &PlanarConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &ConductionResponseReport {
        &self.report
    }
    pub fn into_parts(self) -> (PlanarConductionResponse, ConductionResponseReport) {
        (self.problem, self.report)
    }
}

#[derive(Clone, Copy)]
enum ResponseDrive<'a> {
    Traces {
        voltage: &'a [(usize, f64)],
        currents: &'a [ConductionCurrentTrace],
    },
    Current {
        positive: &'a [usize],
        negative: &'a [usize],
        current_a: f64,
    },
}

impl PlanarConductionResponse {
    /// Bind one source law per original triangle. Every source is a positive
    /// memoryless scalar conductor with explicitly owned magnitude coordinates.
    /// Scalar fixed-state and nonlinear regions may coexist. Ordinary tensor,
    /// field-storage and moving-conductor paths retain their own admission.
    pub fn new(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        Self::new_with_phase_contacts(mesh, sources, Vec::new(), Vec::new(), limits, cancelled)
    }

    /// Nonlinear bulk response and original matched-trace surface conduction.
    /// Contact laws remain static at explicitly owned external/phase states;
    /// only declared bulk E/J coordinates are closed by this response owner.
    pub fn new_with_phase_contacts(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        Self::from_bindings(
            mesh,
            sources,
            contacts,
            liquid_mass_fraction_axes,
            None,
            limits,
            cancelled,
        )
    }

    /// Close each contact's original R'' at its actual A-to-B voltage jump and
    /// normal current density. The source and trace vectors correspond in INPUT
    /// order and are canonicalized together. Both bulk and contact response use
    /// the existing native field iteration, current balance and Joule accounting.
    ///
    /// Electrically varying contacts require uniform jump/current along each trace.
    /// Fixed contacts retain the original nonuniform P1 response.
    /// A nonuniform trace refuses; it is never replaced by a face-average law.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contact_response(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<ConductionContact>,
        contact_laws: Vec<ScalarElectricalInterfaceResponse>,
        contact_controls: ContactResponseControls,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        poll(cancelled)?;
        if contacts.len() != contact_laws.len() || contacts.len() > limits.max_triplets / 16 {
            return Err(ConductionResponseError::Invalid(
                "complete bounded interface response coverage required",
            ));
        }
        let mut ordered: Vec<_> = contacts.into_iter().zip(contact_laws).collect();
        ordered.sort_by_key(|(c, _)| {
            let mut a = c.first_vertices;
            let mut b = c.second_vertices;
            a.sort_unstable();
            b.sort_unstable();
            (a, b)
        });
        let (contacts, laws): (Vec<_>, Vec<_>) = ordered.into_iter().unzip();
        let axes = laws
            .iter()
            .map(|law| law.liquid_mass_fraction_axes().clone())
            .collect();
        let bindings = ContactResponseSources::new(laws, contact_controls)?;
        Self::from_bindings(
            mesh,
            sources,
            contacts,
            axes,
            Some(bindings),
            limits,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_bindings(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        contact_sources: Option<ContactResponseSources>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        poll(cancelled)?;
        if sources.len() != mesh.faces().len() || sources.len() > limits.max_cells {
            return Err(ConductionResponseError::Invalid(
                "one bounded source response per triangle",
            ));
        }
        let mut conductors = Vec::with_capacity(sources.len());
        for source in &sources {
            poll(cancelled)?;
            source.admit(source.original())?;
            conductors.push(source.original().clone());
        }
        let native = PlanarConduction::assemble_with_admitted_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            |_| true,
            |index, source, axes| {
                contact::admit_source(index, source, axes, contact_sources.as_ref())
            },
        )?;
        Ok(Self {
            native,
            sources: sources.into(),
            relaxation: 1.0,
            adaptive_relaxation: None,
            contact_sources,
        })
    }

    /// Use a fixed coordinate-relaxation weight in `(0, 1]` for subsequent
    /// nonlinear solves. One is the original undamped iteration. A smaller
    /// weight can stabilize alternating positive-resistance laws; it is not a
    /// guarantee of convergence and does not weaken any physical residual gate.
    ///
    /// Trials requery the pinned source at interpolated owned E/J magnitudes.
    /// They never interpolate conductivity, temperature or phase. The existing
    /// iteration limit bounds both field solves and extra source queries.
    pub fn with_relaxation(mut self, weight: f64) -> Result<Self, ConductionResponseError> {
        if !(positive(weight) && weight <= 1.0) {
            return Err(ConductionResponseError::Invalid(
                "source-coordinate relaxation must be in (0, 1]",
            ));
        }
        self.relaxation = weight;
        self.adaptive_relaxation = None;
        Ok(self)
    }

    /// Adapt the weight from successive owned source-coordinate defects using
    /// the shared bounded Aitken kernel. Original bulk/contact laws are requeried
    /// at each trial and all actual-field current/power gates remain binding.
    /// History is local to one solve and is discarded on refusal or cancellation.
    /// Calling `with_relaxation` again selects a fixed weight.
    pub fn with_adaptive_relaxation(
        mut self,
        controls: AdaptiveConductionRelaxation,
    ) -> Result<Self, ConductionResponseError> {
        controls.prepare().map_err(|_| {
            ConductionResponseError::Invalid(
                "adaptive source relaxation requires 0 < minimum <= initial <= maximum <= 1",
            )
        })?;
        self.relaxation = controls.initial_weight;
        self.adaptive_relaxation = Some(controls);
        Ok(self)
    }

    pub fn adaptive_relaxation(&self) -> Option<AdaptiveConductionRelaxation> {
        self.adaptive_relaxation
    }

    /// Fixed weight, or the initial adaptive weight, retained across candidates.
    pub fn relaxation(&self) -> f64 {
        self.relaxation
    }

    pub fn mesh(&self) -> &TriComplex2 {
        self.native.mesh()
    }
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.native.conductors()
    }
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        &self.sources
    }
    pub fn contacts(&self) -> &[ConductionContact] {
        self.native.contacts()
    }
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        self.native.contact_liquid_mass_fraction_axes()
    }

    pub fn contact_response_sources(&self) -> Option<&ContactResponseSources> {
        self.contact_sources.as_ref()
    }

    fn admit_contacts(
        &self,
        contacts: &[ConductionContact],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ConductionResponseError> {
        if contacts.len() != self.contacts().len() {
            return Err(ConductionResponseError::Invalid(
                "complete original surface contact coverage required",
            ));
        }
        for (next, original) in contacts.iter().zip(self.contacts()) {
            poll(cancelled)?;
            let before = original.source.resolved();
            let after = next.source.resolved();
            if next.first_vertices != original.first_vertices
                || next.second_vertices != original.second_vertices
                || next.source.surfaces() != original.source.surfaces()
                || before.card_identity() != after.card_identity()
                || before.axis_quantities() != after.axis_quantities()
                || before.query_point().len() != after.query_point().len()
                || before
                    .query_point()
                    .iter()
                    .zip(after.query_point())
                    .any(|((a, _), (b, _))| a != b)
                || before.properties().len() != after.properties().len()
                || before
                    .properties()
                    .iter()
                    .zip(after.properties())
                    .any(|(a, b)| {
                        a.requirement() != b.requirement()
                            || a.source_axes() != b.source_axes()
                            || a.answer().receipt.selected != b.answer().receipt.selected
                    })
            {
                return Err(ConductionResponseError::Invalid(
                    "surface geometry, material and selected sources must remain pinned during response closure",
                ));
            }
        }
        Ok(())
    }

    /// Transfer this operator AND its complete source laws to a continuing
    /// coupled owner. The operator alone is only a frozen linearization and
    /// cannot be advanced as a nonlinear material. Rejoin the same parts with
    /// [`Self::from_frozen_parts`] before every source-response solve.
    /// The source bundle retains contact closure and numerical relaxation too.
    pub fn into_frozen_parts(self) -> (PlanarConduction, ConductionResponseSources) {
        (
            self.native,
            ConductionResponseSources {
                bulk: self.sources,
                relaxation: self.relaxation,
                adaptive_relaxation: self.adaptive_relaxation,
                contacts: self.contact_sources,
            },
        )
    }

    /// Rejoin an external owner's retained spatial operator and original laws
    /// without assembling another matrix or selecting another source. Every
    /// complete conductor is admitted by its original law. Already admitted
    /// native surface contacts and their phase-coordinate ownership are retained.
    pub fn from_frozen_parts(
        native: PlanarConduction,
        sources: ConductionResponseSources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        poll(cancelled)?;
        if sources.bulk.len() != native.conductors().len()
            || sources
                .contacts
                .as_ref()
                .is_some_and(|c| c.laws().len() != native.contacts().len())
        {
            return Err(ConductionResponseError::Invalid(
                "complete bulk response laws required for the retained operator",
            ));
        }
        for (source, conductor) in sources.bulk.iter().zip(native.conductors()) {
            poll(cancelled)?;
            source.admit(conductor)?;
        }
        for (index, (contact, axes)) in native
            .contacts()
            .iter()
            .zip(native.contact_liquid_mass_fraction_axes())
            .enumerate()
        {
            poll(cancelled)?;
            contact::admit_source(index, &contact.source, axes, sources.contacts.as_ref())?;
        }
        Ok(Self {
            native,
            sources: sources.bulk,
            relaxation: sources.relaxation,
            adaptive_relaxation: sources.adaptive_relaxation,
            contact_sources: sources.contacts,
        })
    }

    /// Solve the original voltage/current traces at this fixed external context.
    /// Prescribed current density uses actual native boundary loads, never a
    /// scaled unit solution. Equipotential total-current electrodes use
    /// [`Self::solve_current_candidate`] and its actual terminal constraint.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_candidate(
        &self,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_with_materials_candidate(
            self.conductors().to_vec(),
            voltage,
            currents,
            assembly,
            field,
            controls,
            cancelled,
        )
    }

    /// Close E/J at an explicitly refreshed original material context. An outer
    /// finite-heat owner supplies each proposed T; all other source roles/pins,
    /// geometry and basis stay original. Every rejected trial leaves this owner
    /// unchanged. Actual E.J and terminal work are authoritative; nonlinear
    /// coefficient/current/power defects remain separate numerical residuals.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_with_materials_and_contacts_candidate(
            seeds,
            self.contacts().to_vec(),
            voltage,
            currents,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same actual-amplitude response solve with explicitly refreshed
    /// original interface sources. Trace pairing, source cards, selected claims
    /// and phase-coordinate roles remain pinned; no source policy is rerun.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_driven_with_materials_candidate(
            seeds,
            contacts,
            ResponseDrive::Traces { voltage, currents },
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Close source response under an ideal equipotential total-current drive.
    /// Every frozen FEEC trial enforces the authored current on the original
    /// electrode geometry, then evaluates source data at its actual E and J.
    /// No unit-field source query or post-convergence amplitude scaling occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_current_candidate(
        &self,
        positive: &[usize],
        negative: &[usize],
        current_a: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_current_with_materials_candidate(
            self.conductors().to_vec(),
            positive,
            negative,
            current_a,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same total-current closure at an explicitly refreshed original
    /// thermal/phase context. Source roles, pins and geometry remain original.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_current_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        positive: &[usize],
        negative: &[usize],
        current_a: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_current_with_materials_and_contacts_candidate(
            seeds,
            self.contacts().to_vec(),
            positive,
            negative,
            current_a,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Ideal total-current closure with the same explicitly refreshed bulk and
    /// surface sources. Surface voltage jumps affect every native field trial.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_current_with_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        positive: &[usize],
        negative: &[usize],
        current_a: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_driven_with_materials_candidate(
            seeds,
            contacts,
            ResponseDrive::Current {
                positive,
                negative,
                current_a,
            },
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    fn solve_driven_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        drive: ResponseDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_coupled_with_materials_and_contacts_candidate(
            seeds,
            contacts,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                let field = match drive {
                    ResponseDrive::Traces { voltage, currents } => native
                        .solve_with_current_traces(
                            voltage,
                            currents,
                            field_controls,
                            &mut || cancelled(),
                        )?,
                    ResponseDrive::Current {
                        positive,
                        negative,
                        current_a,
                    } => native.solve_current(
                        positive,
                        negative,
                        current_a,
                        field_controls,
                        &mut || cancelled(),
                    )?,
                };
                Ok::<_, ConductionResponseError>((field, ()))
            },
        )
        .map(|(candidate, ())| candidate)
    }

    /// Close material response at the actual field of an external physical
    /// owner, such as a continuing Kirchhoff circuit. Each callback receives
    /// the current frozen FEEC operator and must return a native field from
    /// that EXACT operator plus its proposed driver state. The callback must
    /// start every trial from its same accepted physical history; a constitutive
    /// iteration does not advance time. Only the driver state paired with the
    /// converged source/field candidate is returned.
    ///
    /// Source evaluation, current/work defects, iteration limits and source
    /// continuity are identical to the voltage- and current-driven paths.
    /// Source data is evaluated at actual amplitude after the callback, never
    /// at a normalized field used internally by a frozen linear solver. The
    /// callback owns its circuit/port balances and any external work budget.
    /// Returned `field_solves` counts callbacks, one native field solve each.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_coupled_with_materials_candidate<D, E>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
        solve: impl for<'p> FnMut(
            &'p PlanarConduction,
            &mut dyn FnMut() -> bool,
        ) -> Result<(ConductionField<'p>, D), E>,
    ) -> Result<(ConductionResponseCandidate, D), E>
    where
        E: From<ConductionResponseError>,
    {
        self.solve_coupled_with_materials_and_contacts_candidate(
            seeds,
            self.contacts().to_vec(),
            assembly,
            field_controls,
            controls,
            cancelled,
            solve,
        )
    }

    /// Couple the same nonlinear bulk closure and refreshed original contacts
    /// to a native circuit/physical driver. Interface laws are held at this
    /// trial's declared thermal/phase context; an outer material owner evolves
    /// those states. Every constitutive trial contains the surface operator,
    /// so its current, voltage jump and Joule heat participate in native work.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_coupled_with_materials_and_contacts_candidate<D, E>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
        mut solve: impl for<'p> FnMut(
            &'p PlanarConduction,
            &mut dyn FnMut() -> bool,
        ) -> Result<(ConductionField<'p>, D), E>,
    ) -> Result<(ConductionResponseCandidate, D), E>
    where
        E: From<ConductionResponseError>,
    {
        poll(cancelled).map_err(ConductionResponseError::from)?;
        controls.validate()?;
        self.native
            .validate_controls(0, field_controls)
            .map_err(ConductionResponseError::from)?;
        if seeds.len() != self.sources.len() || seeds.len() > assembly.max_cells {
            return Err(ConductionResponseError::Invalid(
                "complete bounded source response seeds required",
            )
            .into());
        }
        self.admit_contacts(&contacts, cancelled)?;
        for (source, seed) in self.sources.iter().zip(&seeds) {
            poll(cancelled).map_err(ConductionResponseError::from)?;
            source.admit(seed).map_err(ConductionResponseError::from)?;
        }
        let mut applied = seeds.clone();
        let mut applied_contacts = contacts.clone();
        let mut adaptive = self
            .adaptive_relaxation
            .map(AdaptiveConductionRelaxation::prepare)
            .transpose()
            .map_err(|_| ConductionResponseError::Invalid("invalid adaptive source relaxation"))?;
        let mut operating_coordinates = Vec::new();
        let mut relaxation_queries = 0usize;
        let mut contact_queries = 0usize;
        let (mut maximum_contact_resistance, mut maximum_contact_voltage) = (0.0_f64, 0.0_f64);
        let (mut maximum_sigma, mut maximum_current, mut absolute_power) =
            (0.0_f64, 0.0_f64, 0.0_f64);
        for iteration in 1..=controls.maximum_iterations {
            poll(cancelled).map_err(ConductionResponseError::from)?;
            let native = PlanarConduction::assemble_with_admitted_contacts(
                self.mesh().clone(),
                applied,
                applied_contacts,
                self.contact_liquid_mass_fraction_axes().to_vec(),
                assembly,
                cancelled,
                |_| true,
                |index, source, axes| {
                    contact::admit_source(index, source, axes, self.contact_sources.as_ref())
                },
            )
            .map_err(ConductionResponseError::from)?;
            let (field, driver) = solve(&native, cancelled)?;
            if !std::ptr::eq(field.problem(), &native) {
                return Err(ConductionResponseError::Invalid(
                    "coupled response field belongs to another frozen operator",
                )
                .into());
            }
            if field.boundary_current_a().abs() > field_controls.current_tolerance_a
                || field.maximum_free_current_residual_a() > field_controls.current_tolerance_a
                || field.power_balance_residual_w().abs() > field_controls.power_tolerance_w
            {
                return Err(ConductionResponseError::Field(ConductionError::Balance {
                    current_residual_a: field.boundary_current_a(),
                    maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                    power_residual_w: field.power_balance_residual_w(),
                })
                .into());
            }
            let mut evaluated = Vec::with_capacity(seeds.len());
            let mut next = Vec::with_capacity(seeds.len());
            let mut evaluated_contacts = Vec::new();
            let mut next_contacts = native.contacts().to_vec();
            let mut converged = true;
            maximum_sigma = 0.0;
            maximum_current = 0.0;
            absolute_power = 0.0;
            maximum_contact_resistance = 0.0;
            maximum_contact_voltage = 0.0;
            for (triangle, (source, seed)) in self.sources.iter().zip(&seeds).enumerate() {
                poll(cancelled).map_err(ConductionResponseError::from)?;
                let e = field.electric_field_v_m()[triangle];
                let j = field.current_density_a_m2()[triangle];
                let conductor = source
                    .at_fields(seed, e, j)
                    .map_err(ConductionResponseError::from)?;
                let sigma = finite(conductor.resistivity_ohm_m().recip())
                    .map_err(ConductionResponseError::from)?;
                let used = native.conductivity_s_m[triangle];
                let conductivity_defect_s_m =
                    finite(sigma - used).map_err(ConductionResponseError::from)?;
                let evaluated_j = conductor.conductivity().current_density_a_m2(e);
                let defect = [
                    finite(j[0] - evaluated_j[0]).map_err(ConductionResponseError::from)?,
                    finite(j[1] - evaluated_j[1]).map_err(ConductionResponseError::from)?,
                ];
                let current_defect =
                    finite(defect[0].hypot(defect[1])).map_err(ConductionResponseError::from)?;
                let volume =
                    self.mesh()
                        .face_measure(triangle)
                        .ok_or(ConductionResponseError::Invalid(
                            "native conductor cell volume missing",
                        ))?;
                let power_defect_w = finite(volume * e[0].mul_add(defect[0], e[1] * defect[1]))
                    .map_err(ConductionResponseError::from)?;
                let tolerance = finite(
                    controls.absolute_conductivity_tolerance_s_m
                        + controls.relative_conductivity_tolerance * sigma.max(used),
                )
                .map_err(ConductionResponseError::from)?;
                maximum_sigma = maximum_sigma.max(conductivity_defect_s_m.abs());
                maximum_current = maximum_current.max(current_defect);
                absolute_power = finite(absolute_power + power_defect_w.abs())
                    .map_err(ConductionResponseError::from)?;
                converged &= conductivity_defect_s_m.abs() <= tolerance
                    && current_defect <= controls.current_density_tolerance_a_m2;
                next.push(conductor.clone());
                evaluated.push(ConductionResponseSample {
                    triangle,
                    conductor,
                    conductivity_defect_s_m,
                    current_density_defect_a_m2: defect,
                    power_defect_w,
                });
            }
            if let Some(bindings) = &self.contact_sources {
                for (index, law) in bindings.laws().iter().enumerate() {
                    poll(cancelled).map_err(ConductionResponseError::from)?;
                    let heat = &field.contact_response()[index];
                    let (sample, contact_converged) = contact::evaluate(
                        index,
                        law,
                        &contacts[index].source,
                        &native.contacts()[index].source,
                        heat.voltage_jump_v,
                        heat.current_density_a_m2,
                        heat.power_w,
                        bindings.controls(),
                    )
                    .map_err(ConductionResponseError::from)?;
                    contact_queries =
                        contact_queries
                            .checked_add(2)
                            .ok_or(ConductionResponseError::Invalid(
                                "contact source query count overflows",
                            ))?;
                    maximum_contact_resistance =
                        maximum_contact_resistance.max(sample.maximum_resistance_defect_ohm_m2);
                    maximum_contact_voltage =
                        maximum_contact_voltage.max(sample.maximum_voltage_defect_v);
                    maximum_current =
                        maximum_current.max(sample.maximum_current_density_defect_a_m2);
                    absolute_power = finite(absolute_power + sample.absolute_power_defect_w)
                        .map_err(ConductionResponseError::from)?;
                    converged &= contact_converged;
                    next_contacts[index].source = sample.evaluated[0].clone();
                    evaluated_contacts.push(sample);
                }
            }
            poll(cancelled).map_err(ConductionResponseError::from)?;
            if converged && absolute_power <= field_controls.power_tolerance_w {
                let report = ConductionResponseReport {
                    voltage_v: field.voltage_v().to_vec(),
                    relative_voltage_v: field.relative_voltage_v().to_vec(),
                    reference_voltage_v: field.reference_voltage_v().to_vec(),
                    electric_field_v_m: field.electric_field_v_m().to_vec(),
                    current_density_a_m2: field.current_density_a_m2().to_vec(),
                    cell_power_w: field.cell_power_w().to_vec(),
                    contact_response: field.contact_response().to_vec(),
                    power_w: field.power_w(),
                    nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
                    nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
                    boundary_current_a: field.boundary_current_a(),
                    maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                    power_balance_residual_w: field.power_balance_residual_w(),
                    linear: field.report().cloned(),
                    evaluated,
                    evaluated_contacts,
                    field_solves: iteration,
                    source_evaluations: iteration
                        .checked_mul(seeds.len())
                        .and_then(|count| count.checked_add(relaxation_queries))
                        .and_then(|count| count.checked_add(contact_queries))
                        .ok_or(ConductionResponseError::Invalid(
                            "source evaluation count overflows",
                        ))?,
                };
                return Ok((
                    ConductionResponseCandidate {
                        problem: Self {
                            native,
                            sources: self.sources.clone(),
                            relaxation: self.relaxation,
                            adaptive_relaxation: self.adaptive_relaxation,
                            contact_sources: self.contact_sources.clone(),
                        },
                        report,
                    },
                    driver,
                ));
            }
            let relaxation = if iteration < controls.maximum_iterations {
                if let Some(adaptive) = &mut adaptive {
                    operating_coordinates.clear();
                    for ((law, current), following) in
                        self.sources.iter().zip(native.conductors()).zip(&next)
                    {
                        poll(cancelled).map_err(ConductionResponseError::from)?;
                        for (axis, _) in law.axes() {
                            operating_coordinates.push(
                                current
                                    .material()
                                    .state_coordinate(axis)
                                    .zip(following.material().state_coordinate(axis))
                                    .ok_or(ConductionResponseError::Invalid(
                                        "missing owned source coordinate",
                                    ))?,
                            );
                        }
                    }
                    if let Some(bindings) = &self.contact_sources {
                        for (index, law) in bindings.laws().iter().enumerate() {
                            poll(cancelled).map_err(ConductionResponseError::from)?;
                            for (axis, _) in law.axes() {
                                operating_coordinates.push(
                                    native.contacts()[index]
                                        .source
                                        .resolved()
                                        .state_coordinate(axis)
                                        .zip(
                                            next_contacts[index]
                                                .source
                                                .resolved()
                                                .state_coordinate(axis),
                                        )
                                        .ok_or(ConductionResponseError::Invalid(
                                            "missing owned contact coordinate",
                                        ))?,
                                );
                            }
                        }
                    }
                    adaptive.next_weight(&operating_coordinates).map_err(|_| {
                        ConductionResponseError::Invalid(
                            "invalid adaptive source-coordinate history",
                        )
                    })?
                } else {
                    self.relaxation
                }
            } else {
                1.0
            };
            if iteration < controls.maximum_iterations && relaxation < 1.0 {
                for ((source, current), following) in
                    self.sources.iter().zip(native.conductors()).zip(&mut next)
                {
                    poll(cancelled).map_err(ConductionResponseError::from)?;
                    if !source.axes().is_empty() {
                        *following = source
                            .relax_operating_point(current, following, relaxation)
                            .map_err(ConductionResponseError::from)?;
                        relaxation_queries = relaxation_queries.checked_add(1).ok_or(
                            ConductionResponseError::Invalid("source evaluation count overflows"),
                        )?;
                    }
                }
                if let Some(bindings) = &self.contact_sources {
                    for (index, law) in bindings.laws().iter().enumerate() {
                        poll(cancelled).map_err(ConductionResponseError::from)?;
                        if !law.axes().is_empty() {
                            next_contacts[index].source = law
                                .relax_operating_point(
                                    &native.contacts()[index].source,
                                    &next_contacts[index].source,
                                    relaxation,
                                )
                                .map_err(ContactResponseError::from)
                                .map_err(ConductionResponseError::from)?;
                            relaxation_queries = relaxation_queries.checked_add(1).ok_or(
                                ConductionResponseError::Invalid(
                                    "source evaluation count overflows",
                                ),
                            )?;
                        }
                    }
                }
            }
            applied = next;
            applied_contacts = next_contacts;
        }
        Err(ConductionResponseError::NotConverged {
            iterations: controls.maximum_iterations,
            maximum_conductivity_defect_s_m: maximum_sigma,
            maximum_current_density_defect_a_m2: maximum_current,
            absolute_power_defect_w: absolute_power,
            maximum_contact_resistance_defect_ohm_m2: maximum_contact_resistance,
            maximum_contact_voltage_defect_v: maximum_contact_voltage,
        }
        .into())
    }
}
