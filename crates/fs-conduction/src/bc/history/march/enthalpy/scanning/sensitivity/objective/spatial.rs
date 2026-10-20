//! Mesh-bound spatial readouts of the original nodal thermal fields.
//!
//! A readout is a nonnegative normalized sparse row. Interior P1 probes use
//! barycentric weights; region means use actual tetrahedral volumes; declared
//! sensor or reference-mass weights use `weighted_mean`. Temperature and phase
//! fraction are sampled BEFORE averaging, never inferred from mean enthalpy.
//! Geometry and weights are fixed; sensor dynamics and optical response are
//! not modeled. Volume means of mass fraction are NOT bulk mass fractions:
//! use reference-mass weights for that quantity.

use super::{ScanSensitivityError, ScanObservable, ScanObservation, ScanStep,
    ScanObjectiveBudget, ScanObjectiveEvaluation, ScanObjectiveSampling,
    ScanningEnthalpyRun, EnthalpyStepConfig, LinearConfig, finite, invalid, poll, zeros};
use super::responses::{ScanResponse, ScanResponseBudget, ScanResponseEvaluation};
use super::responses::peaks::{ScanPeakResponse, ScanPeakEvaluation};
use super::checkpoint::{ScanCheckpointBudget, ScanCheckpointWork};
use crate::ConductionMesh;
use fs_exec::Cx;
use std::collections::BTreeMap;

/// One immutable spatial mean tied to the exact mesh used for construction.
/// Cloning copies only its bounded sparse row, not the physical mesh.
#[derive(Debug, Clone)]
pub struct ScanReadout<'m> {
    mesh: &'m ConductionMesh,
    weights: Vec<(usize, f64)>,
}

impl PartialEq for ScanReadout<'_> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.mesh, other.mesh) && self.weights == other.weights
    }
}

impl<'m> ScanReadout<'m> {
    /// Exact original vertex readout, without interpolation or normalization drift.
    pub fn vertex(cx: &Cx<'_>, mesh: &'m ConductionMesh, vertex: usize)
        -> Result<Self, ScanSensitivityError>
    {
        Self::weighted_mean(cx, mesh, &[(vertex, 1.0)], 1)
    }

    /// A P1 point in one explicit tetrahedron. Nonnegative homogeneous
    /// barycentric weights are normalized; `[1;4]` means its centroid. Negative
    /// coordinates (extrapolation) and the all-zero tuple refuse. No nearest-node
    /// snapping or search across coincident material interfaces is performed.
    pub fn in_tetrahedron(cx: &Cx<'_>, mesh: &'m ConductionMesh, element: usize,
        barycentric: [f64; 4]) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let tet = mesh.complex().tets.get(element)
            .ok_or_else(|| invalid("spatial readout names an absent tetrahedron"))?;
        let mut weights: Vec<_> = tet.iter().zip(barycentric)
            .map(|(&v, w)| (v as usize, w)).collect();
        weights.sort_by_key(|&(v, _)| v);
        Self::weighted_mean(cx, mesh, &weights, 4)
    }

    /// Volume mean of a P1 field over selected complete tetrahedra. This is the
    /// exact element integral of the represented linear field, evaluated in
    /// binary64, not a node-count average or an error bound on the continuum.
    /// Elements must be nonempty, strictly increasing, and fit `max_elements`.
    pub fn volume_mean(cx: &Cx<'_>, mesh: &'m ConductionMesh,
        elements: &[usize], max_elements: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        if elements.is_empty() || elements.len() > max_elements {
            return Err(invalid("volume readout needs a nonempty bounded element selection"));
        }
        let mut weights = BTreeMap::<usize, f64>::new();
        let mut previous = None;
        for (i, &element) in elements.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if element >= mesh.element_count() || previous.is_some_and(|old| old >= element) {
                return Err(invalid("volume readout elements must be distinct increasing IDs inside the mesh"));
            }
            previous = Some(element);
            let weight = finite(mesh.element_volume(element) / 4.0)?;
            if weight <= 0.0 { return Err(invalid("volume readout element measure is not representable")); }
            for &vertex in &mesh.complex().tets[element] {
                let entry = weights.entry(vertex as usize).or_default();
                *entry = finite(*entry + weight)?;
            }
        }
        let weights: Vec<_> = weights.into_iter().collect();
        Self::weighted_mean(cx, mesh, &weights, mesh.vertex_count())
    }

    /// A declared finite sensor footprint or reference-mass mean. Entries must
    /// have distinct increasing vertex IDs, finite nonnegative weights and at
    /// least one positive weight. Common weight units cancel on normalization.
    /// Zero entries are removed, not allowed to seed unrelated material corners.
    ///
    /// Admission precedes allocation. Normalization scales by the largest weight
    /// to avoid an overflowing raw sum. A positive coefficient that underflows
    /// to zero refuses instead of silently discarding a sensor contribution.
    pub fn weighted_mean(cx: &Cx<'_>, mesh: &'m ConductionMesh,
        weights: &[(usize, f64)], max_weights: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        if weights.is_empty() || weights.len() > max_weights || weights.len() > mesh.vertex_count() {
            return Err(invalid("spatial readout needs a nonempty bounded sparse row"));
        }
        let mut maximum = 0.0_f64;
        let mut previous = None;
        for (i, &(vertex, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if vertex >= mesh.vertex_count() || previous.is_some_and(|old| old >= vertex)
                || !weight.is_finite() || weight < 0.0
            { return Err(invalid("readout weights need distinct increasing mesh vertices and finite nonnegative values")); }
            previous = Some(vertex);
            maximum = maximum.max(weight);
        }
        if maximum == 0.0 { return Err(invalid("spatial readout has zero total weight")); }
        let mut sum = 0.0;
        for (i, &(_, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            let scaled = weight / maximum;
            if weight > 0.0 && scaled == 0.0 {
                return Err(invalid("positive readout weight underflows during normalization"));
            }
            sum = finite(sum + scaled)?;
        }
        let mut normalized = Vec::new();
        normalized.try_reserve_exact(weights.len())
            .map_err(|_| invalid("spatial readout allocation refused"))?;
        for (i, &(vertex, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if weight == 0.0 { continue; }
            let coefficient = finite((weight / maximum) / sum)?;
            if coefficient == 0.0 { return Err(invalid("positive normalized readout coefficient underflows")); }
            normalized.push((vertex, coefficient));
        }
        poll(cx, weights.len())?;
        Ok(Self { mesh, weights: normalized })
    }

    /// Actual normalized coefficients in ascending original vertex order.
    pub fn weights(&self) -> &[(usize, f64)] { &self.weights }

    /// Immutable mesh identity. Equal vertex counts do not establish identity.
    pub fn mesh(&self) -> &'m ConductionMesh { self.mesh }

    /// Sample a field already expressed in its measured physical coordinate.
    /// Unsampled vertices are not read. A refusal returns no partial scalar.
    pub fn evaluate(&self, cx: &Cx<'_>, field: &[f64]) -> Result<f64, ScanSensitivityError> {
        poll(cx, 0)?;
        if field.len() != self.mesh.vertex_count() { return Err(invalid("readout field must match its exact mesh shape")); }
        let mut value = 0.0;
        for (i, &(vertex, weight)) in self.weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            value = finite(value + finite(weight * finite(field[vertex])?)?)?;
        }
        poll(cx, self.weights.len())?;
        Ok(value)
    }

    /// Transpose the SAME spatial row into the original nodal coordinates.
    /// Thermal-chart conversion is a separate step, after this scatter.
    pub fn pullback(&self, cx: &Cx<'_>, seed: f64) -> Result<Vec<f64>, ScanSensitivityError> {
        poll(cx, 0)?;
        finite(seed)?;
        let mut out = zeros(self.mesh.vertex_count())?;
        self.add_seed(cx, &mut out, seed)?;
        Ok(out)
    }

    fn add_seed(&self, cx: &Cx<'_>, out: &mut [f64], seed: f64) -> Result<(), ScanSensitivityError> {
        for (i, &(vertex, weight)) in self.weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            out[vertex] = finite(out[vertex] + finite(weight * seed)?)?;
        }
        poll(cx, self.weights.len())
    }
}

/// One physical measurement of an interior probe or finite spatial footprint.
/// Its single residual is formed AFTER the weighted readout, not once per node.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanSpatialObservation<'m> {
    /// Zero-based fixed-grid endpoint.
    pub endpoint: usize,
    /// Immutable measurement geometry, tied to the original mesh.
    pub readout: ScanReadout<'m>,
    /// Quantity sampled at each node before applying the spatial row.
    pub observable: ScanObservable,
    /// Target in K, J/kg or dimensionless fraction, according to `observable`.
    pub target: f64,
    /// Positive residual scale in the target's physical units.
    pub scale: f64,
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::ScanObservation {}
    impl Sealed for super::ScanSpatialObservation<'_> {}
}

/// The two supported observation layouts share one physical objective engine.
/// Sealed to preserve the relation between each readout and its exact transpose.
/// Every stored sparse row has an explicit constructor cap; at most
/// `max_observations * mesh.vertex_count()` coefficients can be represented.
pub trait ScanObservationSpec: sealed::Sealed {
    /// Fixed-grid endpoint.
    fn endpoint(&self) -> usize;
    /// Original physical coordinate, before any spatial reduction.
    fn observable(&self) -> ScanObservable;
    /// Physical fitting target.
    fn target(&self) -> f64;
    /// Positive physical residual scale.
    fn scale(&self) -> f64;
    /// Check that the readout addresses this exact mesh before solving physics.
    #[doc(hidden)]
    fn validate_readout(&self, mesh: &ConductionMesh) -> Result<(), ScanSensitivityError>;
    /// Sample the original endpoint without changing it.
    #[doc(hidden)]
    fn prediction(&self, cx: &Cx<'_>, report: &ScanStep) -> Result<f64, ScanSensitivityError>;
    /// Add the spatial transpose to an internal nodal seed. A numerical refusal
    /// may leave a partial seed; the objective owner discards that entire trial.
    #[doc(hidden)]
    fn add_seed(&self, cx: &Cx<'_>, seed: &mut [f64], value: f64) -> Result<(), ScanSensitivityError>;
}

impl ScanObservationSpec for ScanObservation {
    fn endpoint(&self) -> usize { self.endpoint }
    fn observable(&self) -> ScanObservable { self.observable }
    fn target(&self) -> f64 { self.target }
    fn scale(&self) -> f64 { self.scale }
    fn validate_readout(&self, mesh: &ConductionMesh) -> Result<(), ScanSensitivityError> {
        if self.vertex >= mesh.vertex_count() { return Err(invalid("observation vertex is outside its mesh")); }
        Ok(())
    }
    fn prediction(&self, cx: &Cx<'_>, report: &ScanStep) -> Result<f64, ScanSensitivityError> {
        poll(cx, 0)?;
        field(report, self.observable).get(self.vertex).copied()
            .ok_or_else(|| invalid("observation vertex is outside its field")).and_then(finite)
    }
    fn add_seed(&self, cx: &Cx<'_>, seed: &mut [f64], value: f64) -> Result<(), ScanSensitivityError> {
        poll(cx, 0)?;
        let target = seed.get_mut(self.vertex).ok_or_else(|| invalid("observation seed has the wrong shape"))?;
        *target = finite(*target + finite(value)?)?;
        Ok(())
    }
}

impl ScanObservationSpec for ScanSpatialObservation<'_> {
    fn endpoint(&self) -> usize { self.endpoint }
    fn observable(&self) -> ScanObservable { self.observable }
    fn target(&self) -> f64 { self.target }
    fn scale(&self) -> f64 { self.scale }
    fn validate_readout(&self, mesh: &ConductionMesh) -> Result<(), ScanSensitivityError> {
        if !std::ptr::eq(self.readout.mesh, mesh) {
            return Err(invalid("spatial observation belongs to a different mesh; matching vertex counts are insufficient"));
        }
        Ok(())
    }
    fn prediction(&self, cx: &Cx<'_>, report: &ScanStep) -> Result<f64, ScanSensitivityError> {
        self.readout.evaluate(cx, field(report, self.observable))
    }
    fn add_seed(&self, cx: &Cx<'_>, seed: &mut [f64], value: f64) -> Result<(), ScanSensitivityError> {
        poll(cx, 0)?;
        if seed.len() != self.readout.mesh.vertex_count() { return Err(invalid("spatial observation seed has the wrong shape")); }
        finite(value)?;
        self.readout.add_seed(cx, seed, value)
    }
}

fn field(report: &ScanStep, observable: ScanObservable) -> &[f64] {
    match observable {
        ScanObservable::Temperature => &report.solution().temperature,
        ScanObservable::SpecificEnthalpy => &report.solution().specific_enthalpy_j_kg,
        ScanObservable::LiquidMassFraction => &report.solution().liquid_mass_fraction,
    }
}

/// Fit spatial or vertex observations with optional original point/peak limits
/// and checkpointing. The same forward and reverse owners serve both layouts;
/// constraints do not add forward trajectories, and replay checks the actual
/// weighted predictions. Readouts are fixed, not differentiated mesh geometry.
///
/// Existing numerical/material/phase/energy/cancellation refusals propagate.
/// This returns the original physical report and optional measured replay work.
/// Point and peak constraints retain their existing nodal interpretation.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_observations<O: ScanObservationSpec>(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling, step: EnthalpyStepConfig,
    transpose: LinearConfig, budget: ScanObjectiveBudget, responses: &[ScanResponse],
    peaks: &[ScanPeakResponse], response_budget: ScanResponseBudget, max_peak_samples: usize,
    checkpoints: Option<ScanCheckpointBudget>,
) -> Result<(ScanObjectiveEvaluation, Vec<ScanResponseEvaluation>, Vec<ScanPeakEvaluation>,
    Option<ScanCheckpointWork>), ScanSensitivityError>
{
    if let Some(checkpoints) = checkpoints {
        let result = super::checkpoint::evaluate_observations_checkpointed(cx, run, endpoints,
            observations, sampling, step, transpose, budget, responses, peaks,
            response_budget, max_peak_samples, checkpoints)?;
        Ok((result.objective, result.responses, result.peaks, Some(result.work)))
    } else {
        let (objective, points, peaks) = super::responses::peaks::evaluate_observations_with_peaks(
            cx, run, endpoints, observations, sampling, step, transpose, budget, responses,
            peaks, response_budget, max_peak_samples)?;
        Ok((objective, points, peaks, None))
    }
}
