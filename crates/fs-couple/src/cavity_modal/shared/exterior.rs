//! Geometry-bound exterior reaction for a common infinite rigid baffle.
//!
//! Original signed surface traces and the actual circular neck opening share
//! the existing Rayleigh kernel, positive-residue fitter and coupled runtime.
//! This is a convergent area quadrature, not a closed-body exterior BEM solve.
use super::{CavityNeckLoss, ModalCouplingError, RadiationModel, SharedCavity};
use crate::render::plate::impact::radiation::fit;
use fs_bem::helmholtz::{
    HelmholtzError, Medium,
    rayleigh_modal::{AreaSample, PreparedSurface, MAX_MODAL_RAYLEIGH_CELLS},
};
use fs_exec::CancelGate;
use fs_math::{det, eft::two_sum};

/// Circular opening on the SAME baffle plane as the supplied body cells.
/// Radius and gas come from the original cavity, never from an acoustic fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularBaffledMouth {
    /// Opening centre in the common plane [m].
    pub centre_m: [f64; 2],
    /// Uniform radial subdivisions of the actual neck radius.
    pub radial_rings: usize,
    /// Ring i has this many times (2*i+1) azimuthal sectors. Every cell
    /// consequently has the same area, with its exact sector-area centroid.
    pub first_ring_sectors: usize,
}

/// Explicit spatial, frequency and preparation work for the existing owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaffledExteriorOptions {
    /// Sampled angular-frequency band [rad/s], not an out-of-band guarantee.
    pub band_rad_s: [f64; 2],
    /// Odd 33..=257 rows: even rows train, odd rows independently check.
    pub frequency_samples: usize,
    /// Total body and mouth cells, without dropping any original field.
    pub maximum_cells: usize,
    /// Ceiling on frequencies*cells*fields*(cells+fields).
    pub maximum_sweep_terms: u64,
    /// Complete passive pole count; the original runtime also checks capacity.
    pub maximum_poles: usize,
}

/// Numerical resolution and original fitter errors in the common mechanical
/// mass-normalized basis. These are sampled estimates, not continuum bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct BaffledExteriorReport {
    pub options: BaffledExteriorOptions,
    pub mouth: CircularBaffledMouth,
    pub surface_cells: usize,
    pub physical_fields: usize,
    pub acoustic_poles: usize,
    pub sweep_terms: u64,
    pub minimum_equivalent_disc_ppw: f64,
    pub peak_complex_error: f64,
    pub rms_complex_error: f64,
    pub peak_resistance_error: f64,
    pub rms_resistance_error: f64,
}

#[derive(Debug)]
pub enum BaffledExteriorError {
    Invalid(&'static str),
    Surface(HelmholtzError),
    Fit(String),
    Coupling(ModalCouplingError),
    Cancelled,
}
impl core::fmt::Display for BaffledExteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shared cavity baffled exterior: {self:?}")
    }
}
impl std::error::Error for BaffledExteriorError {}

fn invalid(what: &'static str) -> BaffledExteriorError {
    BaffledExteriorError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), BaffledExteriorError> {
    if gate.is_requested() { Err(BaffledExteriorError::Cancelled) } else { Ok(()) }
}

impl SharedCavity {
    /// Compute and attach the full surface-mouth impedance from the original
    /// signed traces. Each body cell carries ALL original surface modes in
    /// flattened component order, in m/s per mass-normalized modal velocity.
    /// Its integrated trace must agree with the cavity's signed volume shape.
    ///
    /// All cells lie in one common plane with the same outward normal. Body
    /// quadrature must exclude the opening. The circular mouth uses its actual
    /// radius and uniform volume flow; only an ExteriorMouth with a declared
    /// flanged exterior is admitted. No exterior end correction is counted twice.
    ///
    /// Spatial resolution and the sampled passive fit are separate controls.
    /// Existing fit error limits are unchanged. The fitted band covers original
    /// positive surface frequencies, the reference neck frequency and the
    /// original conservative upper bound; coupled/out-of-band excitation still
    /// requires an appropriately chosen band. No microphone or extra integrator
    /// is added, and every signed mutual reaction reaches the same runtime.
    pub fn with_baffled_exterior(
        self, body_cells: &[AreaSample], mouth: CircularBaffledMouth,
        options: BaffledExteriorOptions, gate: &CancelGate,
    ) -> Result<(Self, BaffledExteriorReport), BaffledExteriorError> {
        let (model, report) = self.prepare_baffled_exterior(body_cells, mouth, options, gate)?;
        let source_modes: Vec<_> = (0..=self.surface_modes).collect();
        let loaded = self.with_exterior_radiation(&model, &source_modes, options.maximum_poles, gate)
            .map_err(BaffledExteriorError::Coupling)?;
        poll(gate)?;
        Ok((loaded, report))
    }

    /// Pure preparation for this exact geometry, gas and original modal basis.
    /// The returned model uses original surface velocities and PHYSICAL neck
    /// volume flow, as required by with_exterior_radiation. The paired attaching
    /// method preserves this binding automatically; reusing it on a changed
    /// geometry, gas or modal basis is not supported.
    pub fn prepare_baffled_exterior(
        &self, body_cells: &[AreaSample], mouth: CircularBaffledMouth,
        options: BaffledExteriorOptions, gate: &CancelGate,
    ) -> Result<(RadiationModel, BaffledExteriorReport), BaffledExteriorError> {
        poll(gate)?;
        let [low, high] = options.band_rad_s;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || high <= low
            || !(33..=257).contains(&options.frequency_samples) || options.frequency_samples % 2 == 0
            || !(1..=MAX_MODAL_RAYLEIGH_CELLS).contains(&options.maximum_cells)
            || options.maximum_sweep_terms == 0 || !(1..=256).contains(&options.maximum_poles)
            || 8.0 * high * self.system.sample_period_s() >= core::f64::consts::PI
        {
            return Err(invalid("baffled exterior requires a resolved positive band and explicit cell/work/pole limits"));
        }
        if self.source.loss() != CavityNeckLoss::ExteriorMouth
            || self.source.geometry().outer_end != fs_scenario::NeckEnd::Flanged
            || self.system.samples_rendered() != 0 || self.radiation_neck_port.is_some()
        {
            return Err(invalid("baffled exterior requires one cold, internally loaded neck with a flanged exterior"));
        }
        let reference = self.source.isolated_omega_rad_s();
        if reference < low || reference > high
            || self.system.linear_frequency_upper_bound_rad_s() > high
            || self.system.components().iter().flat_map(|m| m.modes()).take(self.surface_modes)
                .any(|m| m.angular_frequency_rad_s > 0.0
                    && (m.angular_frequency_rad_s < low || m.angular_frequency_rad_s > high))
        {
            return Err(invalid("exterior band must cover the original structural and cavity reference frequencies"));
        }
        let (samples, sweep_terms) = self.baffled_samples(body_cells, mouth, options, gate)?;
        let surface = PreparedSurface::new(&samples).map_err(BaffledExteriorError::Surface)?;
        let gas = self.source.gas();
        let medium = Medium { density: gas.density, sound_speed: gas.sound_speed };
        let ppw = surface.equivalent_disc_ppw(high, medium).map_err(BaffledExteriorError::Surface)?;
        let extent = surface.maximum_distance_m();
        let spacing = if extent > 0.0 {
            (core::f64::consts::TAU * gas.sound_speed / (6.0 * extent)).min(high - low)
        } else { high - low };
        let mut omega: Vec<_> = (0..options.frequency_samples).map(|i|
            low + (high - low) * i as f64 / (options.frequency_samples - 1) as f64).collect();
        omega[0] = low;
        omega[options.frequency_samples - 1] = high;
        let mut data = Vec::with_capacity(omega.len());
        for &frequency in &omega {
            poll(gate)?;
            data.push(surface.power_consistent_impedance(frequency, medium)
                .map_err(BaffledExteriorError::Surface)?);
        }
        poll(gate)?;
        // Fit ALL fields in the same mechanical units. Fitting a mixed
        // modal-velocity/volume-flow matrix would let the tiny mouth area alone
        // dominate the error objective and erase weak structural radiation.
        let fitted = fit::fit_resolved(&omega, &data, self.surface_modes + 1, spacing)
            .map_err(BaffledExteriorError::Fit)?;
        poll(gate)?;
        if fitted.model.poles.len() > options.maximum_poles {
            return Err(invalid("complete exterior fit exceeds its pole budget; no radiation field was dropped"));
        }
        let report = BaffledExteriorReport {
            options, mouth, surface_cells: samples.len(), physical_fields: self.surface_modes + 1,
            acoustic_poles: fitted.model.poles.len(), sweep_terms, minimum_equivalent_disc_ppw: ppw,
            peak_complex_error: fitted.peak_error, rms_complex_error: fitted.rms_error,
            peak_resistance_error: fitted.resistance_peak_error,
            rms_resistance_error: fitted.resistance_rms_error,
        };
        let mut model = fitted.model;
        let neck_projection = self.source.neck_volume_shape_m2_per_sqrt_kg();
        for pole in &mut model.poles {
            let original = pole.coupling[self.surface_modes];
            let physical = original / neck_projection;
            if !physical.is_finite() || (original != 0.0 && physical == 0.0) {
                return Err(invalid("physical neck radiation pullback is not representable"));
            }
            pole.coupling[self.surface_modes] = physical;
        }
        Ok((model, report))
    }

    fn baffled_samples(
        &self, body_cells: &[AreaSample], mouth: CircularBaffledMouth,
        options: BaffledExteriorOptions, gate: &CancelGate,
    ) -> Result<(Vec<AreaSample>, u64), BaffledExteriorError> {
        if mouth.centre_m.iter().any(|x| !x.is_finite()) || mouth.radial_rings == 0
            || mouth.first_ring_sectors < 3 || self.surface_modes + 1 > 128 || body_cells.is_empty()
        {
            return Err(invalid("complete body fields and a finite, resolved circular opening are required"));
        }
        let mouth_cells = mouth.radial_rings.checked_mul(mouth.radial_rings)
            .and_then(|n| n.checked_mul(mouth.first_ring_sectors))
            .ok_or_else(|| invalid("mouth quadrature count overflow"))?;
        let cells = body_cells.len().checked_add(mouth_cells)
            .ok_or_else(|| invalid("exterior cell count overflow"))?;
        let fields = (self.surface_modes + 1) as u64;
        let count = cells as u64;
        let sweep_terms = count.checked_mul(fields)
            .and_then(|v| v.checked_mul(count + fields))
            .and_then(|v| v.checked_mul(options.frequency_samples as u64))
            .ok_or_else(|| invalid("exterior sweep work overflow"))?;
        if cells > options.maximum_cells || sweep_terms > options.maximum_sweep_terms {
            return Err(invalid("complete body and mouth exceed the declared quadrature or sweep budget"));
        }
        let radius = self.source.geometry().neck_radius_m;
        let mut volume = vec![0.0; self.surface_modes];
        let mut correction = volume.clone();
        let mut absolute = volume.clone();
        let mut samples = Vec::with_capacity(cells);
        for cell in body_cells {
            poll(gate)?;
            if !cell.area_m2.is_finite() || cell.area_m2 <= 0.0
                || cell.position_m.iter().any(|x| !x.is_finite())
                || cell.mode_shape.len() != self.surface_modes || cell.mode_shape.iter().any(|x| !x.is_finite())
                || (cell.position_m[0] - mouth.centre_m[0]).hypot(cell.position_m[1] - mouth.centre_m[1]) < radius
            {
                return Err(invalid("body quadrature must retain finite original fields outside the actual mouth"));
            }
            for (i, &shape) in cell.mode_shape.iter().enumerate() {
                let weight = cell.area_m2 * shape;
                let (sum, error) = two_sum(volume[i], weight);
                volume[i] = sum;
                correction[i] += error;
                absolute[i] += weight.abs();
            }
            let mut sample = cell.clone();
            sample.mode_shape.push(0.0);
            samples.push(sample);
        }
        for i in 0..self.surface_modes {
            let actual = volume[i] + correction[i];
            let expected = self.column[i] * self.source.neck_area_m2();
            if !actual.is_finite() || !absolute[i].is_finite()
                || (actual - expected).abs() > 128.0 * f64::EPSILON * absolute[i].max(expected.abs())
            {
                return Err(invalid("exterior trace must integrate to the same signed cavity volume projection"));
            }
        }
        let area = self.source.neck_area_m2() / mouth_cells as f64;
        let neck_shape = det::sqrt(self.source.neck_mass_kg()).recip();
        for ring in 0..mouth.radial_rings {
            poll(gate)?;
            let sectors = mouth.first_ring_sectors * (2 * ring + 1);
            let inner = ring as f64 / mouth.radial_rings as f64;
            let outer = (ring + 1) as f64 / mouth.radial_rings as f64;
            let angle = core::f64::consts::PI / sectors as f64;
            // Exact centroid of an annular sector. Evaluate the radial ratio
            // after factoring the powers, avoiding subtraction of close cubes.
            let radial = radius * (2.0 / 3.0) * (outer * outer + outer * inner + inner * inner)
                / (outer + inner) * det::sin(angle) / angle;
            for sector in 0..sectors {
                let theta = core::f64::consts::TAU * (sector as f64 + 0.5) / sectors as f64;
                let mut mode_shape = vec![0.0; self.surface_modes + 1];
                mode_shape[self.surface_modes] = neck_shape;
                samples.push(AreaSample {
                    position_m: [mouth.centre_m[0] + radial * det::cos(theta),
                        mouth.centre_m[1] + radial * det::sin(theta)],
                    area_m2: area, mode_shape,
                });
            }
        }
        Ok((samples, sweep_terms))
    }
}

#[cfg(test)]
mod tests;
