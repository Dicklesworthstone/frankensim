//! Causal molecular propagation of original complete stress histories.
//!
//! The ORIGINAL caloric reservoirs define ell(s)=sqrt(rho/K(s)). Extract the
//! physical wavefront r/c_inf and prepare four bounded, dimensionless symbols
//! `F_m(s)=exp[-r(s ell(s)-s/c_inf)] (c_inf ell(s))^m`, m=0..3.
//! The original Q, Qdot, Qddot and integral then supply every near/far pressure
//! and velocity term. This does not filter pressure alone or differentiate
//! recorded noise inside the propagation operator.
//!
//! BDF2 convolution quadrature uses delta(z)=3/2-2z+z²/2, evaluated at actual
//! right-half-plane Laplace points. Contour inversion prepares causal weights;
//! a separately zero-padded FFT performs LINEAR convolution. Retaining the
//! first N weights is complete for the first N outputs, not a tail model.
//! Reference: Lubich (1988), Numer. Math. 52, 129–145,
//! <https://doi.org/10.1007/BF01398686>. Uncorrected second-order time accuracy
//! requires smooth compatible onset; finite acceleration jumps are admitted
//! without that order claim. Time, contour, original source interpolation and
//! finite-cell errors remain independent refinement axes.

use super::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, invalid, nonfinite, poll,
};
use crate::AeroacError;
use crate::curle3d::quadrupole::recorded::QuadrupoleSampleClock;
use crate::material::MaterialAeroacError;
use fs_fft::{C64 as FftComplex, Fft};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, dd::Dd, det};

mod observe;
/// Complete volume/force/stress composition using this same molecular cache.
pub mod sources;

const SCOPE: &str = "Initially resting prescribed stationary 3D stress sources in one frozen \
homogeneous molecularly relaxing gas, including complete pressure/vector velocity and exact \
ideal images. Frozen-speed wavefronts are explicit; finite-horizon BDF2 convolution quadrature \
has separate timestep, contour and source-interpolation errors. No classical viscothermal \
loss, moving medium/source, finite impedance, general scattering, feedback, automatic \
source/bandwidth validity, or measured-SPL qualification.";

/// One finite, uniform causal observation window, beginning at source onset.
#[derive(Clone, Copy, Debug)]
pub struct MolecularTransientConfig {
    pub timestep_s: f64,
    pub samples: usize,
    /// Power of two, at least twice samples. Independently refine this contour.
    pub transform_len: usize,
    /// rho=exp(-exponent/transform_len), in [8,24]. Aliasing and roundoff
    /// compete; this is a numerical contour choice, not an error certificate.
    pub contour_exponent: f64,
}

/// Explicit numerical payload and work bounds before allocating a cache.
#[derive(Clone, Copy, Debug)]
pub struct MolecularTransientBudget {
    pub maximum_paths: usize,
    pub maximum_samples: usize,
    pub maximum_transform_len: usize,
    /// Peak numerical backing arrays, cached kernels, path geometry and result.
    /// Caller histories, original material metadata and allocator overhead excluded.
    pub maximum_workspace_bytes: usize,
    /// Conservative FFT stage-element count for preparation plus ONE evaluation.
    pub maximum_transform_work: usize,
}

impl Default for MolecularTransientBudget {
    fn default() -> Self {
        Self {
            maximum_paths: 64,
            maximum_samples: 32768,
            maximum_transform_len: 131072,
            maximum_workspace_bytes: 256 * 1024 * 1024,
            maximum_transform_work: 1 << 29,
        }
    }
}

/// Immutable original material, path geometry and four CQ kernels per path.
/// Warm batches use O(paths*L log L + paths*N) fixed-order work, where L is the
/// transform length and N the output count, with no material resampling.
pub struct PreparedMolecularQuadrupoleObserver {
    molecular: RelaxationState,
    geometry: PreparedTransientQuadrupoleObserver,
    config: MolecularTransientConfig,
    fft: Fft,
    // None is the exact empty/zero-reservoir limit: all four filters are identity.
    kernels: Option<Vec<[Vec<FftComplex>; 4]>>,
    workspace_bytes: usize,
}

/// One finite prefix of the received physical field; later sound is unspecified.
#[derive(Clone, Debug)]
pub struct MolecularQuadrupoleRecord {
    clock: QuadrupoleSampleClock,
    fields: Vec<RetardedQuadrupoleField>,
}

impl MolecularQuadrupoleRecord {
    pub const fn clock(&self) -> QuadrupoleSampleClock {
        self.clock
    }
    pub fn fields(&self) -> &[RetardedQuadrupoleField] {
        &self.fields
    }
    /// Original integer clock in DD, retaining sub-ULP path delays at large epochs.
    pub fn time_at_clock_dd(&self, index: usize) -> Option<Dd> {
        (index < self.fields.len()).then(|| {
            Dd::from_f64(self.clock.epoch_s())
                + Dd::from_f64(self.clock.timestep_s()) * Dd::from_f64(index as f64)
        })
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

impl PreparedMolecularQuadrupoleObserver {
    /// Prepare direct paths from the original molecular gas, not an independently
    /// authored propagation speed. Each source retains its physical tensor order.
    pub fn new(
        molecular: RelaxationState,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        config: MolecularTransientConfig,
        budget: MolecularTransientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialAeroacError> {
        poll(&mut cancelled)?;
        let frozen = FluidAcousticState::try_new(
            molecular.gas().density,
            molecular.frozen_sound_speed_m_s(),
        )
        .map_err(MaterialAeroacError::Material)?;
        let geometry = PreparedTransientQuadrupoleObserver::new_cancellable(
            frozen,
            observer_m,
            source_positions_m,
            budget.maximum_paths,
            &mut cancelled,
        )?;
        Self::from_paths(molecular, geometry, config, budget, cancelled)
    }

    /// Bind existing direct/ideal-wall/wedge paths. Their density and speed MUST
    /// equal this material's density and frozen wavefront speed. Every image
    /// keeps its original sign, tensor transformation and separate physical delay.
    pub fn from_paths(
        molecular: RelaxationState,
        geometry: PreparedTransientQuadrupoleObserver,
        config: MolecularTransientConfig,
        budget: MolecularTransientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialAeroacError> {
        poll(&mut cancelled)?;
        let c = molecular.frozen_sound_speed_m_s();
        if geometry.acoustic().density_kg_m3() != molecular.gas().density
            || geometry.acoustic().sound_speed_m_s() != c
        {
            return Err(
                invalid("molecular paths must retain the original frozen-speed medium").into(),
            );
        }
        let paths = geometry.path_count();
        let n = config.transform_len;
        if !config.timestep_s.is_finite()
            || config.timestep_s <= 0.0
            || !config.timestep_s.recip().is_finite()
            || config.samples < 2
            || config.samples > budget.maximum_samples
            || n < 4
            || !n.is_power_of_two()
            || n > budget.maximum_transform_len
            || config.samples > n / 2
            || paths > budget.maximum_paths
            || !config.contour_exponent.is_finite()
            || !(8.0..=24.0).contains(&config.contour_exponent)
        {
            return Err(
                invalid("molecular transient clock, horizon, contour or extent budget").into(),
            );
        }
        let work = paths
            .checked_mul(n)
            .and_then(|v| v.checked_mul(n.ilog2() as usize))
            .and_then(|v| v.checked_mul(28));
        if work.is_none_or(|v| v > budget.maximum_transform_work) {
            return Err(invalid("molecular transient transform work budget").into());
        }
        // 4 cached complex kernels/path; conservative 8 extra complex buffers
        // cover construction, material contour, FFT tables/subtables and work.
        // 16 real projected signals + 4 DD sums + result fields per output.
        let workspace_bytes = paths
            .checked_mul(4)
            .and_then(|v| v.checked_add(8))
            .and_then(|v| v.checked_mul(n))
            .and_then(|v| v.checked_mul(core::mem::size_of::<FftComplex>()))
            .and_then(|v| {
                v.checked_add(config.samples.checked_mul(
                    16 * core::mem::size_of::<f64>()
                        + 4 * core::mem::size_of::<Dd>()
                        + core::mem::size_of::<RetardedQuadrupoleField>(),
                )?)
            })
            .and_then(|v| {
                v.checked_add(paths.checked_mul(core::mem::size_of::<super::Transfer>())?)
            })
            .and_then(|v| {
                v.checked_add(paths.checked_mul(core::mem::size_of::<[Vec<FftComplex>; 4]>())?)
            })
            .filter(|&v| v <= budget.maximum_workspace_bytes)
            .ok_or_else(|| invalid("molecular transient numerical workspace byte budget"))?;
        let fft = Fft::new(n);
        let active = molecular
            .processes()
            .iter()
            .any(|p| p.heat_capacity_j_kg_k != 0.0);
        let kernels = if active && paths != 0 {
            let rho = det::exp(-config.contour_exponent / n as f64);
            let mut contour = Vec::with_capacity(n / 2 + 1);
            for bin in 0..=n / 2 {
                if bin.is_multiple_of(64) {
                    poll(&mut cancelled)?;
                }
                let z = if bin == 0 {
                    C64::from_re(rho)
                } else if bin == n / 2 {
                    C64::from_re(-rho)
                } else {
                    let (sin, cos) = det::sin_cos(-core::f64::consts::TAU * bin as f64 / n as f64);
                    C64::new(cos, sin).scale(rho)
                };
                let s = (C64::from_re(1.5) - z.scale(2.0) + (z * z).scale(0.5))
                    .scale(config.timestep_s.recip());
                let response = molecular
                    .laplace_response(s)
                    .map_err(MaterialAeroacError::Material)?;
                contour.push((
                    response.slowness_s_m.scale(c),
                    response.excess_propagation_rate_per_m,
                ));
            }
            let mut scratch = vec![FftComplex::default(); n];
            let mut kernels = Vec::with_capacity(paths);
            for path in &geometry.transfers {
                poll(&mut cancelled)?;
                let range = path.delay_s * c;
                let mut rows: [Vec<FftComplex>; 4] =
                    core::array::from_fn(|_| vec![FftComplex::default(); n]);
                for (bin, &(slowness, excess)) in contour.iter().enumerate() {
                    if bin.is_multiple_of(64) {
                        poll(&mut cancelled)?;
                    }
                    let exponent = excess.scale(range);
                    if !finite(exponent) {
                        return Err(invalid("molecular path exponent overflow").into());
                    }
                    let (sin, cos) = det::sin_cos(-exponent.im);
                    let mut value = C64::new(cos, sin).scale(det::exp(-exponent.re));
                    for row in &mut rows {
                        if !finite(value) {
                            return Err(invalid("molecular radial symbol overflow").into());
                        }
                        row[bin] = FftComplex::new(value.re, value.im);
                        if bin != 0 && bin != n / 2 {
                            row[n - bin] = row[bin].conj();
                        }
                        value = value * slowness;
                    }
                }
                for row in &mut rows {
                    poll(&mut cancelled)?;
                    fft.inverse(row, &mut scratch);
                    let mut weight = 1.0;
                    for (index, value) in row.iter_mut().enumerate() {
                        if index.is_multiple_of(256) {
                            poll(&mut cancelled)?;
                        }
                        if index < config.samples {
                            let real = value.re / weight;
                            if !real.is_finite() {
                                return Err(nonfinite("molecular CQ coefficient").into());
                            }
                            *value = FftComplex::new(real, 0.0);
                            weight *= rho;
                        } else {
                            *value = FftComplex::default();
                        }
                    }
                    fft.forward(row, &mut scratch);
                    if row.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
                        return Err(nonfinite("molecular convolution kernel").into());
                    }
                }
                kernels.push(rows);
            }
            Some(kernels)
        } else {
            None
        };
        poll(&mut cancelled)?;
        Ok(Self {
            molecular,
            geometry,
            config,
            fft,
            kernels,
            workspace_bytes,
        })
    }

    pub const fn molecular_state(&self) -> &RelaxationState {
        &self.molecular
    }
    pub const fn paths(&self) -> &PreparedTransientQuadrupoleObserver {
        &self.geometry
    }
    pub const fn config(&self) -> MolecularTransientConfig {
        self.config
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub const fn scope(&self) -> &'static str {
        SCOPE
    }
}

#[cfg(test)]
mod tests;
