//! Scalar source closure at the actual tetrahedral electric field and current.
//!
//! Each bounded Picard iteration uses the original native bulk/contact operator.
//! Only potentials are combined inside a frozen terminal response; the entire
//! response is rebuilt after a conductivity change. Nonlinear fields cannot be
//! scaled to another drive. Optional coordinate relaxation requeries the same
//! original source; no coefficient blending or numerical heat is introduced.
//! Positive source values do not guarantee convergence or uniqueness.

use super::*;
use crate::conduction::response::contact::{
    self, ContactResponseControls, ContactResponseError, ContactResponseSample,
    ContactResponseSources,
};
pub use crate::conduction::response::{AdaptiveConductionRelaxation, ConductionResponseSources};
use fs_material::conductor::{ConductorError, response::ScalarConductorResponse};
use fs_material::electrical_interface::response::ScalarElectricalInterfaceResponse;

/// Constitutive and total native-solve limits, separate from field tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionResponseControls {
    pub maximum_iterations: usize,
    /// Boundary drives cost one solve per iteration; N electrodes cost N-1.
    pub maximum_field_solves: usize,
    pub absolute_conductivity_tolerance_s_m: f64,
    pub relative_conductivity_tolerance: f64,
    pub current_density_tolerance_a_m2: f64,
}

#[derive(Debug)]
pub enum SpatialConductionResponseError {
    Field(SpatialConductionError),
    Material(ConductorError),
    Contact(ContactResponseError),
    Invalid(&'static str),
    Cancelled,
    FieldSolveBudget {
        completed: usize,
        required: usize,
        limit: usize,
    },
    NotConverged {
        iterations: usize,
        maximum_conductivity_defect_s_m: f64,
        maximum_current_density_defect_a_m2: f64,
        absolute_power_defect_w: f64,
        maximum_contact_resistance_defect_ohm_m2: f64,
        maximum_contact_voltage_defect_v: f64,
    },
}
impl From<SpatialConductionError> for SpatialConductionResponseError {
    fn from(error: SpatialConductionError) -> Self {
        match error {
            SpatialConductionError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<ConductorError> for SpatialConductionResponseError {
    fn from(error: ConductorError) -> Self {
        Self::Material(error)
    }
}
impl From<ContactResponseError> for SpatialConductionResponseError {
    fn from(error: ContactResponseError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for SpatialConductionResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral conductor source response: {self:?}")
    }
}
impl std::error::Error for SpatialConductionResponseError {}
type Result<T, E = SpatialConductionResponseError> = core::result::Result<T, E>;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    super::check(cancelled).map_err(Into::into)
}
fn finite(value: f64) -> Result<f64> {
    super::finite(value).map_err(Into::into)
}

/// Actual-amplitude drive held throughout this stationary constitutive closure.
#[derive(Clone, Copy, Debug)]
pub enum SpatialConductionResponseDrive<'a> {
    Boundary(&'a SpatialConductionBoundary),
    /// Held positive power between two equipotential face sets. Every source
    /// trial solves its own resistance and closes actual E/J. Supply limits
    /// apply only after material closure, never to intermediate coefficients.
    Power {
        positive_faces: &'a [usize],
        negative_faces: &'a [usize],
        source: crate::conduction::ConductionPowerSource,
        reference_voltage_v: f64,
    },
    /// Signed current INTO the selected electrodes, including zero for floating
    /// electrodes. Other electrodes prescribe voltage in original terminal order;
    /// currents follow the supplied current-terminal order.
    MixedTerminals {
        terminals: &'a [&'a [usize]],
        current_terminals: &'a [usize],
        voltages_v: &'a [f64],
        currents_a: &'a [f64],
        response: SpatialConductanceControls,
        mixed: SpatialMixedConductanceControls,
    },
}
impl<'a> From<&'a SpatialConductionBoundary> for SpatialConductionResponseDrive<'a> {
    fn from(boundary: &'a SpatialConductionBoundary) -> Self {
        Self::Boundary(boundary)
    }
}

/// Same immutable geometry and original source laws; no second current model.
/// The frozen operator stays private so a nonlinear operating point cannot be
/// accidentally exposed as a reusable linear resistance.
#[derive(Clone)]
pub struct SpatialConductionResponse {
    native: SpatialConduction,
    sources: Arc<[ScalarConductorResponse]>,
    relaxation: f64,
    adaptive_relaxation: Option<AdaptiveConductionRelaxation>,
    contact_sources: Option<ContactResponseSources>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialConductionResponseSample {
    pub tetrahedron: usize,
    /// Complete original source evaluated at the actual reported E and J.
    pub conductor: ResolvedConductorState,
    /// Evaluated minus applied conductivity.
    pub conductivity_defect_s_m: f64,
    /// Native J minus source-evaluated sigma(E,J) E, in original mesh axes.
    pub current_density_defect_a_m2: [f64; 3],
    pub power_defect_w: f64,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionResponseReport {
    /// Original native field, including separately retained contact Joule heat.
    pub field: SpatialConductionField,
    pub evaluated: Vec<SpatialConductionResponseSample>,
    pub evaluated_contacts: Vec<ContactResponseSample>,
    pub iterations: usize,
    /// Includes every basis solve in discarded constitutive iterations.
    pub field_solves: usize,
    /// Actual-field and damped-coordinate source queries, excluding initial
    /// binding. Full steps and fixed sources need no extra damping queries.
    pub source_evaluations: usize,
}

pub struct SpatialConductionResponseCandidate {
    problem: SpatialConductionResponse,
    report: SpatialConductionResponseReport,
}
impl SpatialConductionResponseCandidate {
    pub fn problem(&self) -> &SpatialConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &SpatialConductionResponseReport {
        &self.report
    }
    pub fn into_parts(self) -> (SpatialConductionResponse, SpatialConductionResponseReport) {
        (self.problem, self.report)
    }
}

impl SpatialConductionResponse {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_contacts(
            positions_m,
            tetrahedra,
            sources,
            Vec::new(),
            limits,
            cancelled,
        )
    }

    /// One complete scalar source per tetrahedron. The original contact laws
    /// retain fixed external states; this owner closes bulk E/J coordinates only.
    /// Temperature/phase contexts are prescribed here, not advanced thermally.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contacts(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<SpatialConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::from_bindings(
            positions_m,
            tetrahedra,
            sources,
            contacts,
            None,
            limits,
            cancelled,
        )
    }

    /// Original local R'' laws on matched triangular interfaces. Electrically
    /// varying contacts require uniform solved voltage jump/current density; varying traces
    /// refuse instead of acquiring an averaged or fitted constitutive law.
    /// Source bindings follow INPUT contacts and are canonicalized with them.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contact_response(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<SpatialConductionContact>,
        contact_laws: Vec<ScalarElectricalInterfaceResponse>,
        contact_controls: ContactResponseControls,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if contacts.len() != contact_laws.len() || contacts.len() > limits.max_triplets / 36 {
            return Err(SpatialConductionResponseError::Invalid(
                "complete bounded 3D contact source coverage required",
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
        let bindings = ContactResponseSources::new(laws, contact_controls)?;
        Self::from_bindings(
            positions_m,
            tetrahedra,
            sources,
            contacts,
            Some(bindings),
            limits,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_bindings(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<SpatialConductionContact>,
        contact_sources: Option<ContactResponseSources>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if sources.len() != tetrahedra.len()
            || sources.len() > limits.max_cells
            || contacts.len() > limits.max_triplets / 36
        {
            return Err(SpatialConductionResponseError::Invalid(
                "complete bounded source/contact coverage required",
            ));
        }
        let mut conductors = Vec::with_capacity(sources.len());
        for source in &sources {
            check(cancelled)?;
            source.admit(source.original())?;
            conductors.push(source.original().clone());
        }
        let geometry = Arc::new(
            prepare_geometry(positions_m, tetrahedra, limits, cancelled)
                .map_err(SpatialConductionError::from)?,
        );
        // The original laws above own EVERY E/J coordinate. This private
        // assembly seam never relaxes admission for ordinary Ohmic snapshots.
        let contact_axes = contact_sources.as_ref().map_or_else(
            || vec![[None, None]; contacts.len()],
            |bindings| {
                bindings
                    .laws()
                    .iter()
                    .map(|law| law.liquid_mass_fraction_axes().clone())
                    .collect()
            },
        );
        let native =
            SpatialConduction::assemble_admitted(geometry, conductors, limits, cancelled, |_| {
                true
            })?
            .bind_admitted_contacts(
                contacts,
                contact_axes,
                limits,
                cancelled,
                |i, s, a| {
                    contact::admit_source(i, s, a, contact_sources.as_ref()).map_err(|_| {
                        SpatialConductionError::Invalid(
                            "3D contact response source admission refused",
                        )
                    })
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

    /// Select a fixed source-coordinate relaxation weight in `(0, 1]`.
    /// One preserves the original undamped Picard path. Smaller weights can
    /// stabilize alternating resistance responses by requerying original laws
    /// at interpolated owned E/J magnitudes. Temperature/phase, source pins and
    /// all actual-field current/energy acceptance gates remain unchanged.
    /// No global convergence or uniqueness claim follows from relaxation.
    pub fn with_relaxation(mut self, weight: f64) -> Result<Self> {
        if !(super::terminals::positive(weight) && weight <= 1.0) {
            return Err(SpatialConductionResponseError::Invalid(
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
    ) -> Result<Self> {
        controls.prepare().map_err(|_| {
            SpatialConductionResponseError::Invalid(
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

    pub fn mesh(&self) -> &TetComplex {
        self.native.mesh()
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        self.native.positions_m()
    }
    pub fn exterior_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.native.exterior_faces()
    }
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.native.conductors()
    }
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        &self.sources
    }
    pub fn contacts(&self) -> &[SpatialConductionContact] {
        self.native.contacts()
    }
    pub fn contact_response_sources(&self) -> Option<&ContactResponseSources> {
        self.contact_sources.as_ref()
    }
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        self.native.contact_liquid_mass_fraction_axes()
    }
    pub fn contact_cells(&self) -> &[[usize; 2]] {
        self.native.contact_cells()
    }
    pub fn cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.native.cell_volume_m3(cell)
    }

    /// Transfer the admitted native operator together with every original bulk,
    /// contact and relaxation policy to a continuing coupled owner. A frozen
    /// operator alone does not represent a nonlinear operating point.
    pub fn into_frozen_parts(self) -> (SpatialConduction, ConductionResponseSources) {
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

    /// Rejoin an original frozen operator and its complete source policy without
    /// rebuilding geometry or silently dropping nonlinear contact ownership.
    pub fn from_frozen_parts(
        native: SpatialConduction,
        sources: ConductionResponseSources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if sources.bulk.len() != native.conductors().len()
            || sources
                .contacts
                .as_ref()
                .is_some_and(|c| c.laws().len() != native.contacts().len())
        {
            return Err(SpatialConductionResponseError::Invalid(
                "complete original bulk/contact response laws required",
            ));
        }
        for (source, conductor) in sources.bulk.iter().zip(native.conductors()) {
            check(cancelled)?;
            source.admit(conductor)?;
        }
        for (index, (contact, axes)) in native
            .contacts()
            .iter()
            .zip(native.contact_liquid_mass_fraction_axes())
            .enumerate()
        {
            check(cancelled)?;
            contact::admit_source(index, &contact.source, axes, sources.contacts.as_ref())
                .map_err(|_| {
                    SpatialConductionError::Invalid("3D contact response source admission refused")
                })?;
        }
        Ok(Self {
            native,
            sources: sources.bulk,
            relaxation: sources.relaxation,
            adaptive_relaxation: sources.adaptive_relaxation,
            contact_sources: sources.contacts,
        })
    }

    /// Close the original sources against actual E/J at the authored drive.
    /// Native current, Hodge and supplied-work gates remain active. Separate
    /// sigma, current-density and sum-ABSOLUTE-local-power gates admit the
    /// constitutive response; errors in different cells cannot cancel.
    ///
    /// Geometry is reused. Every iteration rebuilds native conductivity/contact
    /// assembly and terminal response, charging all N-1 basis solves before work.
    /// Source bounds, nonconvergence and cancellation publish no candidate.
    /// This is a bounded Picard closure, not a claim of global convergence,
    /// uniqueness, nonuniform nonlinear contact quadrature or a transient thermal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_candidate<'d>(
        &self,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCandidate> {
        self.solve_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            drive,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Close the same original sources on explicitly refreshed external context,
    /// e.g. a finite-enthalpy thermal trial. Geometry, contact pairing, original
    /// source pins and numerical policy are retained; only accepted candidates
    /// may replace an enclosing physical owner's state.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_materials_and_contacts_candidate<'d>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<SpatialConductionContact>,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCandidate> {
        check(cancelled)?;
        let drive = drive.into();
        if let SpatialConductionResponseDrive::Power { source, .. } = drive {
            source
                .validate()
                .map_err(SpatialConductionResponseError::Invalid)?;
        }
        let solves_per_iteration = match drive {
            SpatialConductionResponseDrive::Boundary(_)
            | SpatialConductionResponseDrive::Power { .. } => 1,
            SpatialConductionResponseDrive::MixedTerminals { terminals, .. } => {
                terminals.len().checked_sub(1).filter(|&n| n > 0).ok_or(
                    SpatialConductionResponseError::Invalid(
                        "at least two mixed-drive electrodes required",
                    ),
                )?
            }
        };
        let (candidate, power_resistance_ohm) = self.solve_coupled_candidate(
            seeds,
            contacts,
            solves_per_iteration,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                let mut checkpoint = || cancelled();
                let cancelled = &mut checkpoint;
                let mut power_resistance_ohm = None;
                let field = match drive {
                    SpatialConductionResponseDrive::Boundary(boundary) => {
                        native.solve(boundary, field_controls, cancelled)?
                    }
                    SpatialConductionResponseDrive::Power {
                        positive_faces,
                        negative_faces,
                        source,
                        reference_voltage_v,
                    } => {
                        let response = native.resistance_at_reference_voltage(
                            positive_faces,
                            negative_faces,
                            reference_voltage_v,
                            field_controls,
                            cancelled,
                        )?;
                        power_resistance_ohm = Some(response.resistance_ohm());
                        response.field_at_power(
                            crate::conduction::ConductionPowerSource {
                                maximum_voltage_v: f64::MAX,
                                maximum_current_a: f64::MAX,
                                ..source
                            },
                            field_controls.current_tolerance_a,
                            field_controls.power_tolerance_w,
                            cancelled,
                        )?
                    }
                    SpatialConductionResponseDrive::MixedTerminals {
                        terminals,
                        current_terminals,
                        voltages_v,
                        currents_a,
                        response,
                        mixed,
                    } => {
                        let matrix = native.conductance_matrix(
                            terminals,
                            response,
                            field_controls,
                            cancelled,
                        )?;
                        let mixed =
                            matrix.prepare_current_control(current_terminals, mixed, cancelled)?;
                        mixed.field_at_controls(
                            voltages_v,
                            currents_a,
                            field_controls.current_tolerance_a,
                            field_controls.power_tolerance_w,
                            cancelled,
                        )?
                    }
                };
                Ok::<_, SpatialConductionResponseError>((field, power_resistance_ohm))
            },
        )?;
        if let SpatialConductionResponseDrive::Power { source, .. } = drive {
            source
                .operating_point(power_resistance_ohm.ok_or(
                    SpatialConductionResponseError::Invalid("lost power-drive resistance"),
                )?)
                .map_err(SpatialConductionResponseError::Invalid)?;
        }
        Ok(candidate)
    }

    /// Close a continuing external circuit on every frozen native N-electrode
    /// response. The callback supplies actual terminal voltages and a private
    /// driver proposal; the original matrix recovers the field itself, so an
    /// unrelated or hand-constructed field cannot satisfy constitutive gates.
    ///
    /// Every N-1 basis extraction is charged before callback work. The caller
    /// owns circuit-step and outer thermal budgets and must keep all trials on
    /// one accepted storage history. No nonlinear terminal matrix is cached.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_coupled_terminal_candidate<D, E>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<SpatialConductionContact>,
        terminals: &[&[usize]],
        matrix_controls: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
        mut solve: impl FnMut(
            &SpatialConduction,
            &SpatialConductanceMatrix<'_>,
            &mut dyn FnMut() -> bool,
        ) -> Result<(Vec<f64>, D), E>,
    ) -> Result<(SpatialConductionResponseCandidate, D), E>
    where
        E: From<SpatialConductionResponseError>,
    {
        let solves_per_iteration = terminals.len().checked_sub(1).filter(|&n| n > 0).ok_or(
            SpatialConductionResponseError::Invalid("at least two coupled electrodes required"),
        )?;
        self.solve_coupled_candidate(
            seeds,
            contacts,
            solves_per_iteration,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                let matrix = native
                    .conductance_matrix(terminals, matrix_controls, field_controls, &mut || {
                        cancelled()
                    })
                    .map_err(SpatialConductionResponseError::from)?;
                let (voltages, driver) = solve(native, &matrix, cancelled)?;
                let field = matrix
                    .field_at_voltages(
                        &voltages,
                        field_controls.current_tolerance_a,
                        field_controls.power_tolerance_w,
                        &mut || cancelled(),
                    )
                    .map_err(SpatialConductionResponseError::from)?;
                Ok((field, driver))
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_coupled_candidate<D, E>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<SpatialConductionContact>,
        solves_per_iteration: usize,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
        mut solve: impl FnMut(
            &SpatialConduction,
            &mut dyn FnMut() -> bool,
        ) -> Result<(SpatialConductionField, D), E>,
    ) -> Result<(SpatialConductionResponseCandidate, D), E>
    where
        E: From<SpatialConductionResponseError>,
    {
        check(cancelled)?;
        if !(1..=crate::conduction::response::MAX_RESPONSE_ITERATIONS)
            .contains(&controls.maximum_iterations)
            || controls.maximum_field_solves == 0
            || !controls.absolute_conductivity_tolerance_s_m.is_finite()
            || controls.absolute_conductivity_tolerance_s_m < 0.0
            || !controls.relative_conductivity_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_conductivity_tolerance)
            || (controls.absolute_conductivity_tolerance_s_m == 0.0
                && controls.relative_conductivity_tolerance == 0.0)
            || !super::terminals::positive(controls.current_density_tolerance_a_m2)
        {
            return Err(SpatialConductionResponseError::Invalid(
                "invalid bounded source response controls",
            )
            .into());
        }
        if seeds.len() != self.sources.len()
            || seeds.len() > assembly.max_cells
            || contacts.len() != self.contacts().len()
        {
            return Err(SpatialConductionResponseError::Invalid(
                "complete original bulk/contact coverage required",
            )
            .into());
        }
        for (law, seed) in self.sources.iter().zip(&seeds) {
            check(cancelled)?;
            law.admit(seed)
                .map_err(SpatialConductionResponseError::from)?;
        }
        for (original, next) in self.contacts().iter().zip(&contacts) {
            check(cancelled)?;
            let old = original.source.resolved();
            let new = next.source.resolved();
            if original.first_vertices != next.first_vertices
                || original.second_vertices != next.second_vertices
                || original.source.surfaces() != next.source.surfaces()
                || old.card_identity() != new.card_identity()
                || old.axis_quantities() != new.axis_quantities()
                || old
                    .query_point()
                    .iter()
                    .map(|(n, _)| n)
                    .ne(new.query_point().iter().map(|(n, _)| n))
                || old.properties().len() != new.properties().len()
                || old.properties().iter().zip(new.properties()).any(|(a, b)| {
                    a.requirement() != b.requirement()
                        || a.source_axes() != b.source_axes()
                        || a.answer().receipt.selected != b.answer().receipt.selected
                })
            {
                return Err(SpatialConductionResponseError::Invalid(
                    "contact refresh must preserve original geometry and complete source roles",
                )
                .into());
            }
        }
        let mut applied = seeds.clone();
        let mut applied_contacts = contacts.clone();
        let mut field_solves = 0usize;
        let mut adaptive = self
            .adaptive_relaxation
            .map(AdaptiveConductionRelaxation::prepare)
            .transpose()
            .map_err(|_| {
                SpatialConductionResponseError::Invalid("invalid adaptive source relaxation")
            })?;
        let mut operating_coordinates = Vec::new();
        let mut relaxation_queries = 0usize;
        let mut contact_queries = 0usize;
        let (mut maximum_contact_resistance, mut maximum_contact_voltage) = (0.0_f64, 0.0_f64);
        let (mut maximum_sigma, mut maximum_current, mut absolute_power) =
            (0.0_f64, 0.0_f64, 0.0_f64);
        for iteration in 1..=controls.maximum_iterations {
            check(cancelled)?;
            let needed = field_solves.checked_add(solves_per_iteration);
            if needed.is_none_or(|n| n > controls.maximum_field_solves) {
                return Err(SpatialConductionResponseError::FieldSolveBudget {
                    completed: field_solves,
                    required: solves_per_iteration,
                    limit: controls.maximum_field_solves,
                }
                .into());
            }
            // All applied bundles come from their original admitted source
            // (construction or the preceding actual-field query).
            let native = SpatialConduction::assemble_admitted(
                Arc::clone(&self.native.geometry),
                applied,
                assembly,
                cancelled,
                |_| true,
            )
            .map_err(SpatialConductionResponseError::from)?
            .bind_admitted_contacts(
                applied_contacts,
                self.native.contact_liquid_mass_fraction_axes.clone(),
                assembly,
                cancelled,
                |i, s, a| {
                    contact::admit_source(i, s, a, self.contact_sources.as_ref()).map_err(|_| {
                        SpatialConductionError::Invalid(
                            "3D contact response source admission refused",
                        )
                    })
                },
            )
            .map_err(SpatialConductionResponseError::from)?;
            let (field, driver) = solve(&native, cancelled)?;
            field_solves += solves_per_iteration;
            let mut evaluated = Vec::with_capacity(seeds.len());
            let mut next = Vec::with_capacity(seeds.len());
            let mut next_contacts = native.contacts().to_vec();
            let mut evaluated_contacts = Vec::new();
            let mut converged = true;
            maximum_sigma = 0.0;
            maximum_current = 0.0;
            absolute_power = 0.0;
            maximum_contact_resistance = 0.0;
            maximum_contact_voltage = 0.0;
            let mut power_correction = 0.0;
            for (tetrahedron, (source, seed)) in self.sources.iter().zip(&seeds).enumerate() {
                check(cancelled)?;
                let e = field.electric_field_v_m[tetrahedron];
                let j = field.current_density_a_m2[tetrahedron];
                let conductor = source
                    .at_spatial_fields(seed, e, j)
                    .map_err(SpatialConductionResponseError::from)?;
                let sigma = finite(conductor.resistivity_ohm_m().recip())?;
                let used = native.conductors[tetrahedron].resistivity_ohm_m().recip();
                let conductivity_defect_s_m = finite(sigma - used)?;
                let defect = [
                    finite(j[0] - sigma * e[0])?,
                    finite(j[1] - sigma * e[1])?,
                    finite(j[2] - sigma * e[2])?,
                ];
                let current_defect = finite(defect[0].hypot(defect[1]).hypot(defect[2]))?;
                let volume = native.cell_volume_m3(tetrahedron).ok_or(
                    SpatialConductionResponseError::Invalid("original cell volume absent"),
                )?;
                let power_defect_w = finite(volume * dot(e, defect))?;
                let tolerance = finite(
                    controls.absolute_conductivity_tolerance_s_m
                        + controls.relative_conductivity_tolerance * sigma.max(used),
                )?;
                maximum_sigma = maximum_sigma.max(conductivity_defect_s_m.abs());
                maximum_current = maximum_current.max(current_defect);
                let adjusted = finite(power_defect_w.abs() - power_correction)?;
                let sum = finite(absolute_power + adjusted)?;
                power_correction = finite((sum - absolute_power) - adjusted)?;
                absolute_power = sum;
                converged &= conductivity_defect_s_m.abs() <= tolerance
                    && current_defect <= controls.current_density_tolerance_a_m2;
                next.push(conductor.clone());
                evaluated.push(SpatialConductionResponseSample {
                    tetrahedron,
                    conductor,
                    conductivity_defect_s_m,
                    current_density_defect_a_m2: defect,
                    power_defect_w,
                });
            }
            if let Some(bindings) = &self.contact_sources {
                for (index, law) in bindings.laws().iter().enumerate() {
                    check(cancelled)?;
                    let heat = &field.contacts[index];
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
                    .map_err(SpatialConductionResponseError::from)?;
                    contact_queries = contact_queries.checked_add(3).ok_or(
                        SpatialConductionResponseError::Invalid(
                            "contact source query count overflows",
                        ),
                    )?;
                    maximum_contact_resistance =
                        maximum_contact_resistance.max(sample.maximum_resistance_defect_ohm_m2);
                    maximum_contact_voltage =
                        maximum_contact_voltage.max(sample.maximum_voltage_defect_v);
                    maximum_current =
                        maximum_current.max(sample.maximum_current_density_defect_a_m2);
                    let adjusted = finite(sample.absolute_power_defect_w - power_correction)?;
                    let sum = finite(absolute_power + adjusted)?;
                    power_correction = finite((sum - absolute_power) - adjusted)?;
                    absolute_power = sum;
                    converged &= contact_converged;
                    next_contacts[index].source = sample.evaluated[0].clone();
                    evaluated_contacts.push(sample);
                }
            }
            check(cancelled)?;
            if converged && absolute_power <= field_controls.power_tolerance_w {
                let report = SpatialConductionResponseReport {
                    field,
                    evaluated,
                    evaluated_contacts,
                    iterations: iteration,
                    field_solves,
                    source_evaluations: iteration
                        .checked_mul(seeds.len())
                        .and_then(|count| count.checked_add(relaxation_queries))
                        .and_then(|count| count.checked_add(contact_queries))
                        .ok_or(SpatialConductionResponseError::Invalid(
                            "source evaluation count overflows",
                        ))?,
                };
                return Ok((
                    SpatialConductionResponseCandidate {
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
                        check(cancelled)?;
                        for (axis, _) in law.axes() {
                            operating_coordinates.push(
                                current
                                    .material()
                                    .state_coordinate(axis)
                                    .zip(following.material().state_coordinate(axis))
                                    .ok_or(SpatialConductionResponseError::Invalid(
                                        "missing owned source coordinate",
                                    ))?,
                            );
                        }
                    }
                    if let Some(bindings) = &self.contact_sources {
                        for (index, law) in bindings.laws().iter().enumerate() {
                            check(cancelled)?;
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
                                        .ok_or(SpatialConductionResponseError::Invalid(
                                            "missing owned contact coordinate",
                                        ))?,
                                );
                            }
                        }
                    }
                    adaptive.next_weight(&operating_coordinates).map_err(|_| {
                        SpatialConductionResponseError::Invalid(
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
                    check(cancelled)?;
                    if !source.axes().is_empty() {
                        *following = source
                            .relax_operating_point(current, following, relaxation)
                            .map_err(SpatialConductionResponseError::from)?;
                        relaxation_queries = relaxation_queries.checked_add(1).ok_or(
                            SpatialConductionResponseError::Invalid(
                                "source evaluation count overflows",
                            ),
                        )?;
                    }
                }
                if let Some(bindings) = &self.contact_sources {
                    for (index, law) in bindings.laws().iter().enumerate() {
                        check(cancelled)?;
                        if !law.axes().is_empty() {
                            next_contacts[index].source = law
                                .relax_operating_point(
                                    &native.contacts()[index].source,
                                    &next_contacts[index].source,
                                    relaxation,
                                )
                                .map_err(ContactResponseError::from)
                                .map_err(SpatialConductionResponseError::from)?;
                            relaxation_queries = relaxation_queries.checked_add(1).ok_or(
                                SpatialConductionResponseError::Invalid(
                                    "source query count overflows",
                                ),
                            )?;
                        }
                    }
                }
            }
            applied = next;
            applied_contacts = next_contacts;
        }
        Err(SpatialConductionResponseError::NotConverged {
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
