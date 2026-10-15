//! Regional material geometry in the existing shared mechanical/acoustic bank.
//!
//! Every triangle retains its resolved density, section stiffness, thickness and
//! material orientation. The existing DKT/modal owner reduces it once; the SAME
//! body/cavity image reacts on every string. No averaged wood/metal, preset body
//! frequencies, independent receiver filters or new physical stepper.
//!
//! The chart is a perfectly bonded flat linear plate, not a carved shell or a
//! joint/slip law. Its damping remains the explicitly authored chart loss model;
//! resolving elastic source properties does not certify material-derived loss.
//! All strings use one authored physical bridge footprint. The compact pressure
//! law has no propagation delay or broadband exterior-fluid back-loading.

use super::*;
use super::plate::SharedPlateBridge;
use super::super::plate::{MaterialBodyCavity, reduce_radiating_body, validate_body_medium};
use crate::thin_plate::{PlateChartRadiation, ResolvedPlateChart, certified_chart_radiators};

/// Frozen regional plate sources and the actual bridge/observer geometry.
#[derive(Clone, Debug)]
pub struct SharedMaterialPlateChartBody {
    /// Complete mesh, regional material identities, section/mass prescriptions.
    pub specimen: ResolvedPlateChart,
    /// Boundary support, prestress, eigenvalue window, retained modes, explicit
    /// loss and nodal unit-total-force footprint. Not an authored mode bank.
    pub reduction: PlateChartRadiation,
    /// Exactly one interior bridge law per declared string. All act through the
    /// same nodal footprint, not through unrelated modal attachment vectors.
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
    /// The same gas supplies string drag, cavity compression/inertia and pressure
    /// density. Material and geometry receipts remain attached through scheduling.
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
        validate_body_medium(&gas, options.sample_rate_hz, body.listener_m, cavity.as_ref(), gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let radiators = certified_chart_radiators(chart, r).map_err(ReciprocalStringError::Material)?;
        poll(gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let image = reduce_radiating_body(&radiators, &gas, options.sample_rate_hz,
            body.modal_budget, body.listener_m, cavity.as_ref(), gate)?;
        let mut bank = super::plate::finish_material_body(cx, sources, &body.bridges,
            string_modes, gas, image, law, options, gate)?;
        bank.source.material_chart = Some(body);
        bank.source.material_cavity = cavity;
        poll(gate)?;
        Ok(bank)
    }
}
