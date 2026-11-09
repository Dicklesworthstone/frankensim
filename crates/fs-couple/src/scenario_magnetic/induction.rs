//! Authored magnetic histories on the existing accepted induction timeline.
//!
//! This is an adapter to the native A-phi candidate, not another integrator.
//! Each admitted interval is held or affine in WORLD coordinates, starts at
//! the accepted physical clock and matches its retained prescribed trace.

use super::{
    MagneticPotentialBoundary3, ScenarioMagneticError, native_error, poll, potential_sample,
    vector_potential_boundary_at,
};
use fs_em::spatial::vector_potential::induction::{
    SpatialInduction, SpatialInductionCandidate, SpatialInductionControls,
};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_scenario::payload::{OutsideDomainPolicy, Payload, SampleSource, TableInterpolation};
use fs_scenario::{BcValue, BoundaryCondition, FrameMotion, FrameTree, WORLD};

fn affine_row(
    row: &BoundaryCondition,
    frames: &FrameTree,
    start: f64,
    end: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    // The original sampler has already validated both endpoint frames, source
    // semantics, quantity, component count and polar-vector parity.
    let Some(BcValue::Typed(Payload::Vector(value))) = &row.value else {
        return Err(E::Invalid(
            "induction requires a typed magnetic vector potential",
        ));
    };
    let mut frame = value.meta().frame();
    let mut depth = 0usize;
    while frame != WORLD {
        poll(cancelled)?;
        depth += 1;
        if depth > 256 {
            return Err(E::Invalid(
                "magnetic frame ancestry exceeds the admitted bound",
            ));
        }
        let entry = frames
            .frames
            .iter()
            .find(|entry| entry.id == frame)
            .ok_or(E::Invalid("magnetic frame ancestry is unavailable"))?;
        if !matches!(entry.motion, FrameMotion::Fixed { .. }) {
            return Err(E::Invalid(
                "affine induction traces require fixed frame ancestry",
            ));
        }
        frame = entry.parent;
    }
    match value.source() {
        SampleSource::Fixed(_) => {}
        SampleSource::Distribution { .. } => {
            return Err(E::Invalid(
                "induction cannot infer a realization of a distribution",
            ));
        }
        SampleSource::Table {
            times,
            interpolation,
            outside,
            ..
        } => {
            if *outside == OutsideDomainPolicy::Periodic {
                return Err(E::Invalid(
                    "periodic magnetic histories require an explicit unwrapped interval",
                ));
            }
            if times.len() == 1 && *outside == OutsideDomainPolicy::Clamp {
                return Ok(());
            }
            let next = times.partition_point(|time| time.value <= start);
            if times.get(next).is_some_and(|time| time.value < end) {
                return Err(E::Invalid(
                    "split induction steps at every interior magnetic-history knot",
                ));
            }
            if *interpolation == TableInterpolation::StepLeft
                && potential_sample(row, frames, start, cancelled)?
                    != potential_sample(row, frames, end, cancelled)?
            {
                return Err(E::Invalid(
                    "a discontinuous magnetic trace cannot be replaced by an affine ramp",
                ));
            }
        }
    }
    poll(cancelled)
}

/// Advance canonical magnetic forcing through the native Faraday candidate.
///
/// Samples the original boundary at `state.time_s()` and that clock plus
/// `duration_s`; a separate caller clock cannot desynchronize the physics.
/// Fixed values and nonperiodic linear table segments are admitted. Held
/// StepLeft segments are admitted only when no jump is crossed, including at
/// the endpoint. Interior knots, clamp transitions, periodic tables and dynamic
/// frame ancestry refuse rather than being silently converted into ramps.
/// Callers split intervals at authored knots; this function adds no scheduler.
///
/// Every sampled start edge must agree with the accepted A cochain within the
/// caller's positive `controls.flux_tolerance_wb`. This is a trace-matching
/// allowance in webers, not an error bound on the continuous solution. Native
/// preparation owns the complete fixed support, shared edges, current/charge
/// compatibility, relative gauge and physical energy balance.
///
/// Impressed currents and held mesh-vertex velocities retain their original
/// native meaning. The returned unpublished candidate carries the actual
/// Joule heat, Lorentz work, magnetic storage and boundary work. Rejection or
/// cancellation leaves the accepted field, clock and history unchanged.
#[allow(clippy::too_many_arguments)]
pub fn step_induction3<'m, C: SpatialConductorMaterial>(
    state: &SpatialInduction<'m, C>,
    rows: &[MagneticPotentialBoundary3<'_>],
    frames: &FrameTree,
    endpoint_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    duration_s: f64,
    controls: SpatialInductionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialInductionCandidate<'m, C>, ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    let start = state.time_s();
    let end = start + duration_s;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !start.is_finite()
        || !end.is_finite()
        || end <= start
        || !controls.flux_tolerance_wb.is_finite()
        || controls.flux_tolerance_wb <= 0.0
    {
        return Err(E::Invalid(
            "finite advancing induction time and positive trace tolerance required",
        ));
    }
    let initial = vector_potential_boundary_at(state.problem(), rows, frames, start, cancelled)?;
    let endpoint = vector_potential_boundary_at(state.problem(), rows, frames, end, cancelled)?;
    for row in rows {
        poll(cancelled)?;
        match row {
            MagneticPotentialBoundary3::Uniform { row, .. } => {
                affine_row(row, frames, start, end, cancelled)?
            }
            MagneticPotentialBoundary3::EdgeAverages { values, .. } => {
                for value in values {
                    affine_row(value, frames, start, end, cancelled)?;
                }
            }
        }
    }
    for trace in &initial.potential {
        let edges = state
            .problem()
            .face_edges(trace.face)
            .ok_or(E::Invalid("original magnetic face edges are unavailable"))?;
        for (edge, expected) in edges.into_iter().zip(trace.edge_integrals_wb) {
            poll(cancelled)?;
            let difference = state.potential_wb()[edge] - expected;
            if !difference.is_finite() || difference.abs() > controls.flux_tolerance_wb {
                return Err(E::Invalid(
                    "magnetic history does not start at the accepted boundary trace",
                ));
            }
        }
    }
    let candidate = state
        .step_with_boundary_candidate(
            endpoint_current_a_m2,
            velocity_m_s,
            &endpoint,
            duration_s,
            controls,
            cancelled,
        )
        .map_err(native_error)?;
    poll(cancelled)?;
    Ok(candidate)
}
