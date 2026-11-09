//! A finite turbulent inlet that preserves the caller's prescribed bulk flux.
//!
//! The area-averaged normal velocity is the linear functional L*a of the
//! original four unit-variance Gaussian quadratures. Projecting
//! a_c = a - L^T (L*a)/(L*L^T) gives their exact one-time Gaussian distribution
//! conditioned on zero flux. Its covariance is I-L^T*L/(L*L^T); no energy is
//! restored. See Hoffman & Ribak (1991), DOI 10.1086/186160, for Gaussian
//! constrained realizations. The existing exact face filter defines L, not
//! point samples on a particular mesh.
//!
//! Every polarization remains transverse, so the resulting continuous field
//! stays solenoidal. In the existing anelastic construction rho*u=rho_ref*v,
//! the same constraint also makes the turbulent mass flux through the patch
//! zero, up to representable floating-point roundoff. Mean inflow remains
//! caller-owned. No wall correction is applied.
//!
//! The current phase, amplitudes and clock come from the original source;
//! refresh never draws noise or advances it. This is an instantaneous linear
//! projection of its history, not conditioning the complete OU path on all
//! future constraints. Modal independence and the unconditioned Mann covariance
//! do not hold for the projected frame. Restore the original source checkpoint
//! and rebuild this derived frame to recover the same samples.

use super::{
    FlowSample, MannBoxFilter, MannSpectralModel, Refusal, SpectralMode, det, dot, filter,
    finite_vector, refuse, sample_modes,
};
use crate::{MAX_MODES, fetch_mass_consistent::ExponentialDensityProfile};

/// Conservative face averages on a complete, nonuniform inlet subdivision.
pub mod finite_volume;

/// Fixed, nondegenerate parallelogram in the original Mann physical frame.
/// The ordered half edges define its normal by their right-handed cross product.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannInflowPatch {
    center_m: [f64; 3],
    half_edges_m: [[f64; 3]; 2],
    normal: [f64; 3],
    area_m2: f64,
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use a matching Mann model, a finite nondegenerate inlet patch and a sufficient numerical-storage budget",
    )
}

fn numerical() -> Refusal {
    error(
        "mann-inflow-range",
        "inlet geometry, constraint or projected field is not representable",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error("mann-inflow-cancelled", "inlet refresh cancelled"))
    } else {
        Ok(())
    }
}

fn storage<T>(count: usize) -> Result<Vec<T>, Refusal> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| error("mann-inflow-allocation", "unable to reserve inlet storage"))?;
    Ok(result)
}

fn zeros(count: usize) -> Result<Vec<f64>, Refusal> {
    let mut values = storage(count)?;
    values.resize(count, 0.0);
    Ok(values)
}

impl MannInflowPatch {
    pub fn new(center_m: [f64; 3], half_edges_m: [[f64; 3]; 2]) -> Result<Self, Refusal> {
        finite_vector(center_m)?;
        for edge in half_edges_m {
            finite_vector(edge)?;
        }
        let length = half_edges_m.map(|edge| edge.into_iter().fold(0.0, det::hypot));
        if !length.into_iter().all(|x| x.is_finite() && x > 0.0) {
            return Err(numerical());
        }
        let a = half_edges_m[0].map(|x| x / length[0]);
        let b = half_edges_m[1].map(|x| x / length[1]);
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let magnitude = cross.into_iter().fold(0.0, det::hypot);
        let area_m2 = (4.0 * magnitude * length[0]) * length[1];
        if !area_m2.is_finite() || area_m2 <= 0.0 || magnitude <= 0.0 {
            return Err(numerical());
        }
        Ok(Self {
            center_m,
            half_edges_m,
            normal: cross.map(|x| x / magnitude),
            area_m2,
        })
    }

    pub const fn center_m(&self) -> [f64; 3] {
        self.center_m
    }
    pub const fn half_edges_m(&self) -> [[f64; 3]; 2] {
        self.half_edges_m
    }
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    fn footprint(&self) -> [[f64; 3]; 3] {
        [self.half_edges_m[0], self.half_edges_m[1], [0.0; 3]]
    }
}

/// Reusable derived realization with one exact area-integral flux constraint.
/// Refresh accepts any original temporal source sharing the immutable model;
/// it does not own a separate stochastic trajectory. Sampling and successful
/// refresh allocate nothing. Cancelled/refused refresh retains the prior frame.
#[derive(Clone, Debug)]
pub struct MannInflowFrame {
    model: MannSpectralModel,
    patch: MannInflowPatch,
    filter: MannBoxFilter,
    normal_coefficients: Vec<[f64; 2]>,
    active_constraint: bool,
    time_s: f64,
    amplitudes: [Vec<f64>; 4],
    trial: [Vec<f64>; 4],
}

impl PartialEq for MannInflowFrame {
    fn eq(&self, other: &Self) -> bool {
        self.model == other.model
            && self.patch == other.patch
            && self.filter == other.filter
            && self.normal_coefficients == other.normal_coefficients
            && self.active_constraint == other.active_constraint
            && self.time_s == other.time_s
            && self.amplitudes == other.amplitudes
    }
}

impl MannInflowFrame {
    /// Numerical capacity for geometry, filter, two coefficients per mode,
    /// and all eight accepted/trial amplitude arrays. Fixed metadata excluded.
    pub fn required_workspace_bytes(mode_count: usize) -> Result<usize, Refusal> {
        if !(1..=MAX_MODES).contains(&mode_count) {
            return Err(numerical());
        }
        Ok(mode_count * (core::mem::size_of::<SpectralMode>() + 11 * core::mem::size_of::<f64>()))
    }

    /// Actual retained numerical capacities, including the copied spectrum.
    pub fn workspace_bytes(&self) -> usize {
        self.model.modes.capacity() * core::mem::size_of::<SpectralMode>()
            + self.filter.workspace_bytes()
            + self.normal_coefficients.capacity() * core::mem::size_of::<[f64; 2]>()
            + self
                .amplitudes
                .iter()
                .chain(&self.trial)
                .map(|v| v.capacity() * core::mem::size_of::<f64>())
                .sum::<usize>()
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn patch(&self) -> MannInflowPatch {
        self.patch
    }

    /// Same finite Fourier synthesis and analytic spatial gradient, evaluated
    /// with the accepted projected quadratures. No new boundary at z=0.
    pub fn sample(&self, point: [f64; 3]) -> Result<FlowSample, Refusal> {
        sample_modes(&self.model, self.time_s, self.paths(), point, None)
    }

    /// Exact face/cell average of the projected field, using the existing filter.
    pub fn sample_filtered(
        &self,
        point: [f64; 3],
        filter: &MannBoxFilter,
    ) -> Result<FlowSample, Refusal> {
        sample_modes(
            &self.model,
            self.time_s,
            self.paths(),
            point,
            filter.gains(&self.model)?,
        )
    }

    /// Original anelastic density transformation of this constrained field.
    /// Its normal turbulent mass-flux integral vanishes through the same patch
    /// when the entire patch lies inside the supplied density profile's domain.
    pub fn sample_anelastic(
        &self,
        point: [f64; 3],
        density: &ExponentialDensityProfile,
    ) -> Result<crate::AnelasticFlowSample, Refusal> {
        density.compose_anelastic(None, self.sample(point)?, point[2])
    }

    /// Exact represented spatial momentum/stress, including the modal
    /// correlations introduced by the constraint; no ensemble replacement.
    pub fn sample_filtered_moments(
        &self,
        point: [f64; 3],
        observation: &filter::moments::MannBoxMoments,
        cancelled: impl FnMut() -> bool,
    ) -> Result<filter::moments::MannFilteredMoments, Refusal> {
        filter::moments::sample(
            &self.model,
            self.time_s,
            self.paths(),
            point,
            observation,
            cancelled,
        )
    }

    /// Exact Favre mass/momentum of the same constrained anelastic realization.
    pub fn sample_anelastic_filtered_moments(
        &self,
        point: [f64; 3],
        observation: &filter::anelastic::MannAnelasticBoxMoments,
        cancelled: impl FnMut() -> bool,
    ) -> Result<filter::anelastic::MannFavreMoments, Refusal> {
        filter::anelastic::sample(
            &self.model,
            self.time_s,
            self.paths(),
            point,
            observation,
            cancelled,
        )
    }

    fn paths(&self) -> [&[f64]; 4] {
        core::array::from_fn(|i| self.amplitudes[i].as_slice())
    }
}

pub(super) fn prepare(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    patch: MannInflowPatch,
    maximum_workspace_bytes: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MannInflowFrame, Refusal> {
    poll(&mut cancelled)?;
    let n = model.mode_count();
    if MannInflowFrame::required_workspace_bytes(n)? > maximum_workspace_bytes {
        return Err(error(
            "mann-inflow-budget",
            "inlet storage exceeds its byte budget",
        ));
    }
    let filter = MannBoxFilter::prepare(
        model,
        patch.footprint(),
        maximum_workspace_bytes,
        &mut cancelled,
    )?;
    let gains = filter.gains(model)?.expect("nondegenerate face");
    let mut coefficients = storage(n)?;
    let mut norm = 0.0;
    for (m, mode) in model.modes.iter().enumerate() {
        poll(&mut cancelled)?;
        let mut pair = [0.0; 2];
        for p in 0..2 {
            let normal = dot(patch.normal, mode.polarization[p]);
            pair[p] = normal * gains[m];
            if !pair[p].is_finite() || (normal != 0.0 && gains[m] != 0.0 && pair[p] == 0.0) {
                return Err(numerical());
            }
            norm = det::hypot(norm, pair[p]);
        }
        coefficients.push(pair);
    }
    if !norm.is_finite() {
        return Err(numerical());
    }
    if norm != 0.0 {
        for pair in &mut coefficients {
            for value in pair {
                let normalized = *value / norm;
                if !normalized.is_finite() || (*value != 0.0 && normalized == 0.0) {
                    return Err(numerical());
                }
                *value = normalized;
            }
        }
    }
    let mut modes = storage(n)?;
    modes.extend(model.modes.iter().cloned());
    let retained = MannSpectralModel {
        geometry_seed: model.geometry_seed,
        parameters: model.parameters,
        quadrature: model.quadrature,
        quadrature_rule: model.quadrature_rule,
        advection_mps: model.advection_mps,
        modes,
    };
    let mut frame = MannInflowFrame {
        model: retained,
        patch,
        filter,
        normal_coefficients: coefficients,
        active_constraint: norm != 0.0,
        time_s: time,
        amplitudes: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
        trial: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
    };
    if frame.workspace_bytes() > maximum_workspace_bytes {
        return Err(error(
            "mann-inflow-budget",
            "allocated inlet capacities exceed the byte budget",
        ));
    }
    refresh(model, time, amplitudes, &mut frame, &mut cancelled)?;
    Ok(frame)
}

pub(super) fn refresh(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    frame: &mut MannInflowFrame,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), Refusal> {
    poll(&mut cancelled)?;
    frame.filter.gains(model)?;
    if !time.is_finite() || amplitudes.iter().any(|p| p.len() != model.mode_count()) {
        return Err(numerical());
    }
    let mut projection = 0.0;
    let mut compensation = 0.0;
    let mut denominator = 0.0;
    if frame.active_constraint {
        let position =
            core::array::from_fn(|i| frame.patch.center_m[i] - model.advection_mps[i] * time);
        finite_vector(position)?;
        for (m, mode) in model.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            let phase = dot(mode.k, position);
            if !phase.is_finite() {
                return Err(numerical());
            }
            let (s, c) = det::sin_cos(phase);
            let [a, b] = frame.normal_coefficients[m];
            let weights = [a * c, a * s, b * c, b * s];
            for p in 0..4 {
                frame.trial[p][m] = weights[p];
                let term = weights[p] * amplitudes[p][m] - compensation;
                let sum = projection + term;
                compensation = (sum - projection) - term;
                projection = sum;
                denominator += weights[p] * weights[p];
            }
        }
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(numerical());
        }
        projection /= denominator;
        if !projection.is_finite() {
            return Err(numerical());
        }
    }
    for m in 0..model.mode_count() {
        poll(&mut cancelled)?;
        for p in 0..4 {
            let value = if frame.active_constraint {
                amplitudes[p][m] - projection * frame.trial[p][m]
            } else {
                amplitudes[p][m]
            };
            if !value.is_finite() {
                return Err(numerical());
            }
            frame.trial[p][m] = value;
        }
    }
    poll(&mut cancelled)?;
    core::mem::swap(&mut frame.amplitudes, &mut frame.trial);
    frame.time_s = time;
    Ok(())
}

#[cfg(test)]
mod tests;
