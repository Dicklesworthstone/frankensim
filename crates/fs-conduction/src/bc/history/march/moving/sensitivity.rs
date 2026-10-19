//! Derivatives of the SAME normalized P1 absorbed-power source as the primal.
//!
//! Controls are instantaneous world center, physical half widths, and absorbed
//! watts at a fixed time. The normalization changes when the footprint moves or
//! widens, including at a clipped domain boundary. Its derivative is retained:
//! moving a source redistributes its fixed watts rather than creating heat.
//! This is not a derivative of mesh geometry, schedule times, event selection,
//! the support-count admission rule, or incident optical absorption.

/// Chain sampled physical controls into authored scalar schedule values.
pub mod schedule;
mod power;
/// Authored path and footprint edits for joint physical inverse design.
pub mod geometry;

use super::{MovingHeatSample, MovingHeatSource, invalid, poll, zeros};
use crate::{ConductionError, ConductionMesh, ScalarField, COMPONENT_POWER_DIMS};
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_scenario::signal::TimeSignal;

/// A physical perturbation of one sampled source, with all histories held fixed.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MovingHeatDirection {
    /// World center perturbation, metres.
    pub center_m: [f64; 3],
    /// World-axis half-width perturbation, metres.
    pub half_width_m: [f64; 3],
    /// Absorbed-power perturbation, watts.
    pub absorbed_power_w: f64,
}

/// Derivative of a supplied scalar objective with respect to sampled controls.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MovingHeatGradient {
    /// Objective derivative per metre of world center translation.
    pub center_m: [f64; 3],
    /// Objective derivative per metre of physical half width.
    pub half_width_m: [f64; 3],
    /// Objective derivative per absorbed watt.
    pub absorbed_power_w: f64,
}

/// Immutable source/mesh binding with O(vertices) storage, not a dense Jacobian.
/// The unit-power projection is admitted even when the primal power is zero;
/// an off source outside the mesh has a valid primal but no admitted power
/// derivative. At zero power only feasible nonnegative-power directions apply.
#[derive(Debug)]
pub struct MovingHeatLinearization<'m> {
    sample: MovingHeatSample,
    mesh: &'m ConductionMesh,
    half_width_m: [f64; 3],
    unit_source: ScalarField,
    volumes: Vec<f64>,
}

impl MovingHeatSource {
    /// Reuse the original projection and bind its instantaneous derivatives.
    /// Does not differentiate the sampled time, path coefficients or knots.
    /// The compact squared weight is C1 where it enters/leaves support; its
    /// first derivative is zero on that boundary. No C2 claim is made there.
    ///
    /// # Errors
    /// Original source refusal, unresolved unit-power footprint, cancellation,
    /// or unrepresentable unit projection/geometry. Work is two source samples
    /// plus one volume pass, under the source's original mesh-size caps.
    pub fn linearize_at<'m>(&self, cx: &Cx<'_>, mesh: &'m ConductionMesh, time_s: f64)
        -> Result<MovingHeatLinearization<'m>, ConductionError>
    { self.linearize_sample(cx, mesh, time_s, false) }

    /// Same derivative binding with exact end-left held path/power values.
    /// The chosen side of the event is fixed, not differentiated or averaged.
    pub fn linearize_before<'m>(&self, cx: &Cx<'_>, mesh: &'m ConductionMesh, time_s: f64)
        -> Result<MovingHeatLinearization<'m>, ConductionError>
    { self.linearize_sample(cx, mesh, time_s, true) }

    fn linearize_sample<'m>(&self, cx: &Cx<'_>, mesh: &'m ConductionMesh,
        time_s: f64, before: bool) -> Result<MovingHeatLinearization<'m>, ConductionError>
    {
        let sample = self.sample(cx, mesh, time_s, before)?;
        // The original owner computes the exact same normalized spatial shape.
        // Construct constants instead of cloning a possibly large time table.
        let unit = Self {
            center_m: sample.center_m.map(|v| TimeSignal::Constant(QtyAny::new(v, super::LENGTH))),
            half_width_m: self.half_width_m,
            absorbed_power_w: TimeSignal::Constant(QtyAny::new(1.0, COMPONENT_POWER_DIMS)),
            budget: self.budget,
        }.sample(cx, mesh, 0.0, false)?;
        let mut volumes = zeros(mesh.vertex_count())?;
        for (e, tet) in mesh.complex().tets.iter().enumerate() {
            if e % 256 == 0 { poll(cx, e)?; }
            let quarter = mesh.element_volume(e) / 4.0;
            for &v in tet { volumes[v as usize] = finite(volumes[v as usize] + quarter)?; }
        }
        poll(cx, mesh.vertex_count())?;
        Ok(MovingHeatLinearization { sample, mesh, half_width_m: self.half_width_m,
            unit_source: unit.source, volumes })
    }
}

impl MovingHeatLinearization<'_> {
    /// Original primal sample; no independently reimplemented forward model.
    #[must_use]
    pub const fn sample(&self) -> &MovingHeatSample { &self.sample }

    /// Contract a P1 source-density objective seed with all seven controls.
    /// `seed[i]` is dJ/dq_i, NOT a nodal heat-load seed. The enthalpy adjoint's
    /// `source_density` is already in these coordinates; applying the P1 mass
    /// matrix again here would count the spatial integration twice.
    ///
    /// # Errors
    /// Wrong seed shape, nonfinite data/arithmetic, or cancellation. No partial
    /// gradient is returned. Numerical/material admission decisions stay fixed.
    pub fn pullback(&self, cx: &Cx<'_>, seed: &[f64])
        -> Result<MovingHeatGradient, ConductionError>
    {
        poll(cx, 0)?;
        if seed.len() != self.volumes.len() {
            return Err(invalid("moving-source derivative needs one seed per P1 source-density node"));
        }
        let mut result = MovingHeatGradient::default();
        for (i, &value) in seed.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            finite(value)?;
            result.absorbed_power_w = finite(result.absorbed_power_w + value * self.unit_source.at(i))?;
        }
        for (i, &value) in seed.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            // q_i=P*w_i/S. The subtracted term is the pullback through S;
            // it remains present when the footprint is clipped by the domain.
            let centered = finite(value - self.volumes[i] * result.absorbed_power_w)?;
            let derivatives = self.numerator_derivatives(i)?;
            for axis in 0..3 {
                let factor = finite(self.sample.absorbed_power_w * centered)?;
                result.center_m[axis] = finite(result.center_m[axis] + factor * derivatives[axis])?;
                result.half_width_m[axis] = finite(result.half_width_m[axis] + factor * derivatives[axis + 3])?;
            }
        }
        poll(cx, seed.len())?;
        Ok(result)
    }

    /// Apply the source-density Jacobian without constructing its seven columns.
    /// Return d(q_i) in nodal order. The derivative of total deposited watts is
    /// exactly the power direction in exact arithmetic; center/width changes
    /// contribute zero to that total, including on clipped footprints.
    pub fn apply(&self, cx: &Cx<'_>, direction: MovingHeatDirection)
        -> Result<Vec<f64>, ConductionError>
    {
        poll(cx, 0)?;
        for value in direction.center_m.into_iter().chain(direction.half_width_m)
            .chain(std::iter::once(direction.absorbed_power_w)) { finite(value)?; }
        let mut numerator = zeros(self.volumes.len())?;
        let mut normalization = 0.0;
        for (i, out) in numerator.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            let derivatives = self.numerator_derivatives(i)?;
            for axis in 0..3 {
                *out = finite(*out + derivatives[axis] * direction.center_m[axis])?;
                *out = finite(*out + derivatives[axis + 3] * direction.half_width_m[axis])?;
            }
            normalization = finite(normalization + self.volumes[i] * *out)?;
        }
        for (i, out) in numerator.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            let shape = self.unit_source.at(i);
            let shape_change = finite(*out - shape * normalization)?;
            *out = finite(direction.absorbed_power_w * shape + self.sample.absorbed_power_w * shape_change)?;
        }
        poll(cx, numerator.len())?;
        Ok(numerator)
    }

    // Return dw_i/S, BEFORE the quotient correction; no derivative of the
    // numerically stabilizing max-weight factor is needed because it cancels.
    fn numerator_derivatives(&self, i: usize) -> Result<[f64; 6], ConductionError> {
        let shape = self.unit_source.at(i);
        if shape == 0.0 { return Ok([0.0; 6]); }
        let mut r = [0.0; 3];
        let mut r2 = 0.0;
        for axis in 0..3 {
            r[axis] = (self.mesh.positions()[i][axis] - self.sample.center_m[axis]) / self.half_width_m[axis];
            r2 += r[axis] * r[axis];
        }
        let remainder = 1.0 - r2;
        if !remainder.is_finite() || remainder <= 0.0 {
            return Err(invalid("moving-source derivative lost the admitted positive footprint support"));
        }
        // shape/(1-r²) = (1-r²)/S avoids squaring the small normalization.
        let coefficient = finite(4.0 * (shape / remainder))?;
        let mut result = [0.0; 6];
        for axis in 0..3 {
            result[axis] = finite((coefficient * r[axis]) / self.half_width_m[axis])?;
            result[axis + 3] = finite((coefficient * r[axis] * r[axis]) / self.half_width_m[axis])?;
        }
        Ok(result)
    }
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(invalid("moving-source derivative arithmetic is not representable")) }
}
