//! Regional material geometry in the existing shared mechanical/acoustic bank.
//!
//! Every triangle retains its resolved density, section stiffness, thickness and
//! material orientation. The existing DKT/modal owner reduces it once; the SAME
//! body/cavity image reacts on every string. No averaged wood/metal, preset body
//! frequencies, independent receiver filters or new physical stepper.
//!
//! The chart is a perfectly bonded flat linear plate, not a carved shell or a
//! joint/slip law. Loss is either explicitly authored or selected from each
//! region's isotropic Kelvin-Voigt source. The explicit tensor constructor retains
//! full nonproportional viscous coupling; elastic data alone do not imply loss.
//! The original constructors use a common bridge footprint; the spatial variants
//! retain one physical footprint per string. Both share the same body eigenbasis.
//! The compact pressure law has no delay or broadband exterior-fluid back-loading.

/// Independent spatial body contacts through the same material assembly.
pub mod spatial;
/// Material strain histories coupled to every string through the original body.
pub mod relaxation;
pub use spatial::PlateBridgeFootprints;

use super::*;
use crate::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
use super::plate::SharedPlateBridge;
use super::super::plate::{MaterialBodyCavity, reduce_radiating_body, validate_body_medium};
use crate::thin_plate::{PlateChartRadiation, ResolvedPlateChart};

/// Frozen regional plate sources and the actual bridge/observer geometry.
#[derive(Clone, Debug)]
pub struct SharedMaterialPlateChartBody {
    /// Complete mesh, regional material identities, section/mass prescriptions.
    /// A selected regional bending loss is retained and used in the body modes.
    pub specimen: ResolvedPlateChart,
    /// Boundary support, prestress, eigenvalue window, retained modes and nodal
    /// unit-total-force footprint. Authored damping must be zero with source loss.
    pub reduction: PlateChartRadiation,
    /// Exactly one interior bridge law per declared string. All act through the
    /// same nodal footprint unless explicit spatial footprints are selected.
    pub bridges: Vec<SharedPlateBridge>,
    /// Maximum full plate DOFs (three per node), checked before FE assembly.
    pub max_plate_dofs: usize,
    /// Maximum triangles, checked before FE assembly and the modal solve.
    pub max_plate_elements: usize,
    /// Original body modal time-domain admission.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Compact baffled half-space observer distance [m].
    pub listener_m: f64,
}

impl ReciprocalStringBank {
    /// Compile regional solid sources, arbitrary admitted plate mesh, all strings
    /// and an optional common cavity into the existing two-way friction network.
    ///
    /// Replacing one named region changes the actual element mass/stiffness,
    /// certified modes, signed bridge participation and radiating areas together.
    /// Selected material viscosity contributes actual regional curvature loss,
    /// including conservative prestress dilution, to the same retained modes.
    /// The exterior gas supplies string drag and pressure density. An admitted
    /// distinct cavity gas supplies internal compression/inertia. Material and
    /// geometry receipts remain attached through scheduling.
    /// A new realization starts at rest; this does not morph a vibrating state.
    ///
    /// Geometry assembly/eigensolving is synchronous in its original owner. The
    /// explicit mesh budgets bound admission; cancellation is checked before and
    /// after that call, not falsely promised inside the eigensolver.
    ///
    /// # Errors
    /// Incomplete bridges, excessive mesh/network dimensions, invalid controls,
    /// detached force footprint, stale cavity gas, material/modal refusal or
    /// cancellation. No partially assembled bank is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_chart(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_material_plate_chart_prepared(cx, sources, gas, body, cavity,
            None, law, options, gate)
    }

    /// Select the shared ordinary or prepared owner before original string
    /// histories are restored. The same retained plate basis, bridges and
    /// source budgets are used by both owners.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_chart_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, cavity: Option<MaterialBodyCavity>,
        prepared: impl Into<ModalPreparation>, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate_chart(cx, sources, gas, body, cavity, prepared.into(),
            law, options, gate, None, None)
    }

    /// Keep every regional viscous cross term in the shared string/body/bow solve.
    /// The viscosity must be bound to THIS exact chart. Material tensors are
    /// projected on the original elastic modes alongside bridge/radiation maps,
    /// then realized as passive dashpot columns in the original connection solver.
    /// No diagonal modal loss, material-name dispatch or extra time solver.
    ///
    /// Modes remain a truncated real basis, not complex modes of the damped
    /// continuum. Full modal C is retained; mesh/modal/time refinement is still
    /// necessary. Strong loss outside the retained spatial basis is not certified.
    /// Source-bound zero viscosity retains the old zero-loss arithmetic exactly.
    /// Original material bands, network force/energy budgets and cancellation apply.
    #[allow(clippy::too_many_arguments)]
    pub fn from_viscous_material_plate_chart(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, viscosity: PlateViscosity,
        limits: PlateViscosityLimits, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_viscous_material_plate_chart_prepared(cx, sources, gas, body,
            viscosity, limits, cavity, None, law, options, gate)
    }

    /// Preserve the full viscous matrix with preparation selected before any
    /// original string material histories enter the shared connected state.
    #[allow(clippy::too_many_arguments)]
    pub fn from_viscous_material_plate_chart_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, viscosity: PlateViscosity,
        limits: PlateViscosityLimits, cavity: Option<MaterialBodyCavity>,
        prepared: impl Into<ModalPreparation>, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate_chart(cx, sources, gas, body, cavity, prepared.into(),
            law, options, gate, Some((viscosity, limits)), None)
    }

    #[allow(clippy::too_many_arguments)]
    fn build_material_plate_chart(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, cavity: Option<MaterialBodyCavity>,
        prepared: ModalPreparation, law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
        tensor: Option<(PlateViscosity, PlateViscosityLimits)>,
        footprints: Option<PlateBridgeFootprints>,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let chart = body.specimen.chart();
        let dofs = chart.mesh.nodes.len().checked_mul(3)
            .ok_or_else(|| invalid("regional plate DOF count overflow"))?;
        if dofs == 0 || dofs > body.max_plate_dofs
            || chart.mesh.tris.is_empty() || chart.mesh.tris.len() > body.max_plate_elements {
            return Err(invalid("regional plate exceeds its explicit mesh assembly budget"));
        }
        let r = &body.reduction;
        if r.n_modes == 0 || !r.damping_ratio.is_finite() || r.damping_ratio < 0.0
            || !r.eigenvalue_window.0.is_finite() || !r.eigenvalue_window.1.is_finite()
            || r.eigenvalue_window.0 < 0.0 || r.eigenvalue_window.1 <= r.eigenvalue_window.0
            || r.unit_force_weights.len() != chart.mesh.nodes.len()
            || r.unit_force_weights.iter().any(|x| !x.is_finite()) {
            return Err(invalid("regional plate requires valid modal controls and a complete finite nodal force footprint"));
        }
        let string_modes = super::plate::material_body_preflight(
            &sources, &body.bridges, r.n_modes, cavity.as_ref(), &options, gate)?;
        if let Some(ports) = &footprints {
            if ports.weights.len() != sources.len() {
                return Err(invalid("spatial body contacts require exactly one footprint per material string"));
            }
        }
        validate_body_medium(&gas, options.sample_rate_hz, body.listener_m, cavity.as_ref(), gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        if let Some((viscosity, limits)) = &tensor {
            viscosity.preflight(&body.specimen, r.n_modes, *limits, gate)
                .map_err(ReciprocalStringError::Material)?;
            let channels = if viscosity.is_zero() { 0 } else { r.n_modes };
            let links = sources.len() + channels + cavity.as_ref().map_or(0, |c|
                1 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0));
            if links > options.network.max_connections
                || links > crate::render::schedule::force::coupled::MAX_MODAL_CONNECTIONS {
                return Err(invalid("complete material damping, bridge and cavity channels exceed the original connection budget"));
            }
        }
        // Force maps, full material damping and pressure maps must come from
        // ONE eigenbasis, including mode signs and any degenerate eigenspace.
        // Retain the original complete DKT vectors, including rotations, so
        // later source-bound memory never repeats an eigensolve or reconstructs
        // curvature from a transverse nodal trace. The original mesh/mode
        // budgets bound this extra retained image.
        let basis_terms = dofs.checked_mul(r.n_modes)
            .ok_or_else(|| invalid("retained body constitutive basis size overflow"))?;
        let projection_terms = footprints.as_ref().map_or(basis_terms,
            |p| p.max_projection_terms);
        let (radiators, ports, damping, basis) = body.specimen.material_radiators_with_basis_and_force_ports(
            r, footprints.as_ref().map_or(&[][..], |p| p.weights.as_slice()),
            projection_terms, tensor.as_ref().map(|(v,l)| (v,*l,gate)), gate)
            .map_err(ReciprocalStringError::Material)?;
        let viscous = tensor.map(|(source,_)| crate::plate_viscosity::ViscousPlateImage {
            source, radiators: radiators.clone(), damping: damping.expect("selected original tensor is projected"),
        });
        poll(gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let physical_radiators = viscous.as_ref().map_or(radiators.as_slice(), |v|v.radiators());
        let mut image = reduce_radiating_body(physical_radiators, &gas, options.sample_rate_hz,
            body.modal_budget, body.listener_m, cavity.as_ref(), None, gate)?;
        if let Some(v) = &viscous { image.body_connections.extend(v.damping().connections(0)); }
        let mut bank = super::plate::finish_material_body_with_ports(cx, sources, &body.bridges,
            string_modes, gas, image, footprints.as_ref().map(|_| ports.as_slice()), prepared, law, options, gate)?;
        bank.source.material_chart = Some(body);
        bank.source.material_chart_basis = Some(basis);
        bank.source.material_chart_viscosity = viscous;
        bank.source.material_chart_bridges = footprints;
        bank.source.material_cavity = cavity;
        poll(gate)?;
        Ok(bank)
    }
}
