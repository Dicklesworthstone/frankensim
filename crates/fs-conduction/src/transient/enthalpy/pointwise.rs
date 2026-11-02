//! One enthalpy boundary slot for ambient or reflected pointwise radiation.
//! Both operators use the original nonsymmetric Newton-FGMRES owner.

use crate::{
    ConductionError, PointwiseAmbientRadiation, PointwiseEnclosureAssembly,
    PointwiseEnclosureRadiation, PointwiseEnclosureReport, PointwiseRadiationAssembly,
    PointwiseRadiationReport, MaterialPointwiseAmbientRadiation,
    MaterialPointwiseEnclosureRadiation, MaterialPointwiseEnclosureAssembly,
};
use fs_exec::Cx;

#[derive(Clone, Copy)]
pub(super) enum PointwiseBoundary<'a> {
    Ambient(&'a PointwiseAmbientRadiation<'a>),
    MaterialAmbient(&'a MaterialPointwiseAmbientRadiation<'a>),
    Enclosure(&'a PointwiseEnclosureRadiation<'a>),
    MaterialEnclosure(&'a MaterialPointwiseEnclosureRadiation<'a>),
}

#[derive(Debug)]
pub(super) enum PointwiseBoundaryAssembly {
    Ambient(PointwiseRadiationAssembly),
    Enclosure(PointwiseEnclosureAssembly),
    MaterialEnclosure(MaterialPointwiseEnclosureAssembly),
}

pub(super) enum PointwiseBoundaryReport {
    Ambient(PointwiseRadiationReport),
    Enclosure(PointwiseEnclosureReport),
}

impl PointwiseBoundary<'_> {
    pub(super) fn assemble(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        specific_h: &[f64],
    ) -> Result<PointwiseBoundaryAssembly, ConductionError> {
        match self {
            Self::Ambient(owner) => owner
                .assemble(cx, temperature)
                .map(PointwiseBoundaryAssembly::Ambient),
            Self::MaterialAmbient(owner) => owner
                .assemble_at_temperature(cx, specific_h, temperature)
                .map(PointwiseBoundaryAssembly::Ambient),
            Self::Enclosure(owner) => owner
                .assemble(cx, temperature)
                .map(PointwiseBoundaryAssembly::Enclosure),
            Self::MaterialEnclosure(owner) => owner
                .assemble_at_temperature(cx, specific_h, temperature)
                .map(PointwiseBoundaryAssembly::MaterialEnclosure),
        }
    }

    pub(super) fn load_and_report(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        specific_h: &[f64],
    ) -> Result<(Vec<f64>, PointwiseBoundaryReport), ConductionError> {
        match self {
            Self::Ambient(owner) => owner
                .load_and_report(cx, temperature)
                .map(|(load, report)| (load, PointwiseBoundaryReport::Ambient(report))),
            Self::MaterialAmbient(owner) => owner
                .load_and_report(cx, specific_h, temperature)
                .map(|(load, report)| (load, PointwiseBoundaryReport::Ambient(report))),
            Self::Enclosure(owner) => owner
                .load_and_report(cx, temperature)
                .map(|(load, report)| (load, PointwiseBoundaryReport::Enclosure(report))),
            Self::MaterialEnclosure(owner) => owner
                .load_and_report(cx, specific_h, temperature)
                .map(|(load, report)| (load, PointwiseBoundaryReport::Enclosure(report))),
        }
    }
}

impl PointwiseBoundaryAssembly {
    pub(super) fn admit_adjoint(&self, cx: &Cx<'_>) -> Result<(), ConductionError> {
        match self {
            Self::Ambient(assembly) => assembly.admit_adjoint(cx),
            Self::Enclosure(assembly) => assembly.admit_adjoint(cx),
            Self::MaterialEnclosure(assembly) => assembly.admit_adjoint(cx),
        }
    }

    pub(super) fn parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<crate::radiation::pointwise::PointwiseRadiationParameterGradient, ConductionError> {
        match self {
            Self::Ambient(assembly) => assembly.parameter_pullback(cx, outward_load_weights),
            Self::Enclosure(assembly) => assembly.parameter_pullback(cx, outward_load_weights),
            Self::MaterialEnclosure(assembly) => assembly.parameter_pullback(cx, outward_load_weights),
        }
    }

    pub(super) fn add_scaled_tangent_transpose(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        column_scale: Option<&[f64]>,
        scale: f64,
        include_row: impl Fn(usize) -> bool,
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        match self {
            Self::Ambient(assembly) => assembly.add_scaled_tangent_transpose(
                cx, direction, column_scale, scale, include_row, output,
            ),
            Self::Enclosure(assembly) => assembly.add_scaled_tangent_transpose(
                cx, direction, column_scale, scale, include_row, output,
            ),
            // Material columns already include the original T'(h) and epsilon'(h).
            Self::MaterialEnclosure(assembly) => assembly.add_scaled_tangent_transpose(
                cx, direction, scale, include_row, output,
            ),
        }
    }

    // Use the positive local emission diagonal as the enclosure Picard
    // preconditioner; all nonsymmetric reflection stays in the true action.
    pub(super) fn preconditioner_diagonal(&self, vertex: usize) -> f64 {
        match self {
            Self::Ambient(assembly) => assembly.temperature_diagonal_w_k[vertex],
            Self::Enclosure(assembly) => assembly.emission_diagonal()[vertex],
            Self::MaterialEnclosure(assembly) => assembly.emission_temperature_diagonal_w_k[vertex],
        }
    }

    pub(super) fn add_scaled_tangent(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        column_scale: Option<&[f64]>,
        scale: f64,
        include_row: impl Fn(usize) -> bool,
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        match self {
            Self::Ambient(assembly) => {
                assembly.add_scaled_tangent(cx, direction, column_scale, scale, include_row, output)
            }
            Self::Enclosure(assembly) => {
                assembly.add_scaled_tangent(cx, direction, column_scale, scale, include_row, output)
            }
            Self::MaterialEnclosure(assembly) => {
                assembly.add_scaled_tangent(cx, direction, scale, include_row, output)
            }
        }
    }
}

impl PointwiseBoundaryReport {
    pub(super) fn outward_heat_w(&self) -> f64 {
        match self {
            Self::Ambient(report) => report.outward_heat_w,
            Self::Enclosure(report) => report.outward_heat_w,
        }
    }

    pub(super) fn into_reports(
        self,
    ) -> (
        Option<PointwiseRadiationReport>,
        Option<PointwiseEnclosureReport>,
    ) {
        match self {
            Self::Ambient(report) => (Some(report), None),
            Self::Enclosure(report) => (None, Some(report)),
        }
    }
}

