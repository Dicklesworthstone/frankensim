//! Original nonlinear scalar sources at the actual meridional E and J.
//!
//! The private shared response engine closes complete conductivity/resistivity
//! laws at each cell's recovered field, using the exact swept volume and native
//! radial contact moments. Voltage/current traces, total-current electrodes and
//! regulated power all solve at actual amplitude inside the source iteration.
//! Accepted nonlinear fields cannot be rescaled into another operating point.
//!
//! Static/phase interfaces retain their full P1 surface response. Electrically
//! varying interfaces require uniform endpoint jump and normal current density;
//! nonuniform traces refuse rather than substituting a surface-average law.
//! Source, algebraic, current and work residuals have the existing independent
//! gates. Relaxed source coordinates never become an extra heat contribution.
//!
//! Geometry is fixed, fields are local [radial, axial], and every physical
//! quantity refers to the authored angular span. This scalar memoryless model
//! adds no azimuthal current, Hall response, storage, motion or hysteresis.
//! Positive source values alone prove neither convergence nor uniqueness.

use super::admit_terminal_traces;
use crate::conduction::{
    ConductionContact, ConductionControls, ConductionCurrentTrace, ConductionPowerSource, poll,
    response::{
        AdaptiveConductionRelaxation, ConductionResponseCandidate, ConductionResponseControls,
        ConductionResponseError, ConductionResponseReport, PlanarConductionResponse,
        contact::{ContactResponseControls, ContactResponseSources},
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::{ResolvedConductorState, response::ScalarConductorResponse},
    electrical_interface::response::ScalarElectricalInterfaceResponse,
};
use fs_rep_mesh::TriComplex2;

/// Source-bound nonlinear transport on the original (r,z) mesh.
/// The frozen planar implementation is private and never escapes as a scalable
/// resistance, moving conductor or independently accepted Ohmic operator.
#[derive(Clone)]
pub struct AxisymmetricConductionResponse {
    core: PlanarConductionResponse,
}

/// Applied native field and evaluated source residuals at one accepted drive.
/// Both report and continuing model retain the actual swept geometry. A rejected
/// trial leaves its input model unchanged; continuation consumes the two together.
pub struct AxisymmetricConductionResponseCandidate {
    problem: AxisymmetricConductionResponse,
    report: ConductionResponseReport,
}
impl AxisymmetricConductionResponseCandidate {
    pub fn problem(&self) -> &AxisymmetricConductionResponse {
        &self.problem
    }

    /// Cell components are [radial, axial]; watts refer to the angular span.
    /// Surface losses are separate from bulk powers and included once in total
    /// power. Applied and evaluated sources remain separate observations.
    pub fn report(&self) -> &ConductionResponseReport {
        &self.report
    }

    pub fn into_parts(self) -> (AxisymmetricConductionResponse, ConductionResponseReport) {
        (self.problem, self.report)
    }

    fn from_core(candidate: ConductionResponseCandidate) -> Self {
        let (core, report) = candidate.into_parts();
        Self {
            problem: AxisymmetricConductionResponse { core },
            report,
        }
    }
}

fn admit_metric(
    mesh: &TriComplex2,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ConductionResponseError> {
    poll(cancelled)?;
    if !mesh.metric().is_axisymmetric() {
        return Err(ConductionResponseError::Invalid(
            "axisymmetric source response requires an explicit radial metric",
        ));
    }
    Ok(())
}

impl AxisymmetricConductionResponse {
    /// Bind one original positive scalar law per meridional triangle.
    /// Fixed and explicitly E/J-dependent regions may coexist. Ordinary Ohmic
    /// admission remains separate and cannot freeze these operating axes.
    pub fn new(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        Self::new_with_phase_contacts(mesh, sources, Vec::new(), Vec::new(), limits, cancelled)
    }

    /// Add original matched swept interfaces at explicitly supplied external
    /// temperature/phase states. Electrical contact axes need the separate
    /// `new_with_contact_response` constructor; no contact law is inferred.
    pub fn new_with_phase_contacts(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        admit_metric(&mesh, cancelled)?;
        Ok(Self {
            core: PlanarConductionResponse::from_metric_bindings(
                mesh,
                sources,
                contacts,
                liquid_mass_fraction_axes,
                None,
                limits,
                cancelled,
            )?,
        })
    }

    /// Close each original interface law at its actual A-to-B voltage jump
    /// and normal current density. Trace/source pairs are supplied in the same
    /// INPUT order and canonicalized together. Electrically varying traces
    /// must be uniform within their declared dimensional gates.
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
        admit_metric(&mesh, cancelled)?;
        Ok(Self {
            core: PlanarConductionResponse::from_metric_contact_response(
                mesh,
                sources,
                contacts,
                contact_laws,
                contact_controls,
                limits,
                cancelled,
            )?,
        })
    }

    /// Requery original sources at relaxed owned electrical coordinates.
    /// The weight lies in (0,1]; all native and constitutive gates remain binding.
    pub fn with_relaxation(self, weight: f64) -> Result<Self, ConductionResponseError> {
        Ok(Self {
            core: self.core.with_relaxation(weight)?,
        })
    }

    /// Use the original bounded Aitken source-coordinate policy. Its history
    /// belongs to one solve; only the configured policy continues on acceptance.
    pub fn with_adaptive_relaxation(
        self,
        controls: AdaptiveConductionRelaxation,
    ) -> Result<Self, ConductionResponseError> {
        Ok(Self {
            core: self.core.with_adaptive_relaxation(controls)?,
        })
    }

    pub fn relaxation(&self) -> f64 {
        self.core.relaxation()
    }

    pub fn adaptive_relaxation(&self) -> Option<AdaptiveConductionRelaxation> {
        self.core.adaptive_relaxation()
    }

    pub fn mesh(&self) -> &TriComplex2 {
        self.core.mesh()
    }

    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.core.conductors()
    }

    pub fn sources(&self) -> &[ScalarConductorResponse] {
        self.core.sources()
    }

    /// Original source order after native canonical trace admission.
    pub fn contacts(&self) -> &[ConductionContact] {
        self.core.contacts()
    }

    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        self.core.contact_liquid_mass_fraction_axes()
    }

    pub fn contact_response_sources(&self) -> Option<&ContactResponseSources> {
        self.core.contact_response_sources()
    }

    /// Close original voltage and inward-current-density traces at actual
    /// field amplitude. Current loads use the exact radial surface moments;
    /// every unconstrained physical boundary is insulating.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_candidate(
        &self,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        self.solve_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            voltage,
            currents,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same source closure with explicitly refreshed original bulk states.
    /// An outer owner supplies temperature/phase; all source identities,
    /// coordinates, material roles and meridional geometry remain original.
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
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
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

    /// Refresh original bulk and interface thermal/phase context together,
    /// then close their electrical coordinates at the same native field.
    /// Bulk E.J and exact surface heat remain the physical loss; conductivity,
    /// current and absolute constitutive work defects remain numerical gates.
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
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        self.core
            .solve_with_materials_and_contacts_candidate(
                seeds,
                contacts,
                voltage,
                currents,
                assembly,
                field_controls,
                controls,
                cancelled,
            )
            .map(AxisymmetricConductionResponseCandidate::from_core)
    }

    /// Hold an equipotential terminal's signed TOTAL current. Each frozen
    /// native solve imposes that constraint before evaluating its actual E/J.
    /// Both terminals need their own positive-area physical surface.
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
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        self.solve_current_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
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

    /// Total-current closure with refreshed original bulk states.
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
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
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

    /// Total-current closure with refreshed original bulk and contact states.
    /// Every trial retains the exact swept geometry and source/contact roles.
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
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        admit_terminal_traces(self.mesh(), &[positive, negative], cancelled)?;
        self.core
            .solve_current_with_materials_and_contacts_candidate(
                seeds,
                contacts,
                positive,
                negative,
                current_a,
                assembly,
                field_controls,
                controls,
                cancelled,
            )
            .map(AxisymmetricConductionResponseCandidate::from_core)
    }

    /// Hold terminal power while closing bulk and interface response at the
    /// actual operating point. Compliance is checked on the converged field;
    /// an accepted nonlinear response is never scaled afterwards.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_power_candidate(
        &self,
        positive: &[usize],
        negative: &[usize],
        source: ConductionPowerSource,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        self.solve_power_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            positive,
            negative,
            source,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Regulated-power closure at refreshed original bulk/contact context.
    /// Physical source domains, independent current/work checks, constitutive
    /// iteration limits and cancellation apply to every native trial.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_power_with_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        positive: &[usize],
        negative: &[usize],
        source: ConductionPowerSource,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseCandidate, ConductionResponseError> {
        admit_terminal_traces(self.mesh(), &[positive, negative], cancelled)?;
        self.core
            .solve_power_with_materials_and_contacts_candidate(
                seeds,
                contacts,
                positive,
                negative,
                source,
                assembly,
                field_controls,
                controls,
                cancelled,
            )
            .map(AxisymmetricConductionResponseCandidate::from_core)
    }
}
