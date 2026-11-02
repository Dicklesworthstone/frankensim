//! Loaded resonance poles of the SAME accepted assembly and continuous load.
//!
//! The retained transient generator is scaled and passed to the shared L1
//! nonsymmetric spectral service. No MIMO realization, structural projection,
//! pole pruning, surrogate damping or independent equation assembly occurs.

use super::PreparedAssemblyMimoTransient;
use fs_exec::Cx;
use fs_math::{c64::C64, det};
use fs_spectral::dense::{
    DenseSpectrumError, DenseSpectrumOptions, eigenpair_row_residuals, nonsymmetric_spectrum,
};
use fs_vfit::passivity::CertificateClass;

/// Explicit L1 analysis, physical-unit conversion and output budgets.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyResonanceControls {
    pub numerical: DenseSpectrumOptions,
    /// Positive [rad/s] reference, independent of the transient sample clock.
    /// Eigenvalues are solved for s/reference and then returned in physical units.
    pub reference_angular_frequency_rad_s: f64,
    /// Additional retained/input/output matrix and state payload, excluding the
    /// separately bounded L1 workspace, borrowed assembly and load source.
    pub maximum_payload_bytes: usize,
    /// Conservative coordinate conversion and original-row audit work: 32*d^3.
    pub maximum_back_transform_terms: u64,
    /// Maximum original [q,p,load] row-norm residual, in (0,1). Both this and
    /// the independently scaled L1 residual must pass for every returned pole.
    pub maximum_original_relative_residual: f64,
}

#[derive(Debug)]
pub enum AssemblyResonanceError {
    Invalid(&'static str),
    Numerical(DenseSpectrumError),
    OriginalResidual {
        index: usize,
        kinematic: f64,
        momentum: f64,
        load: f64,
    },
    Cancelled,
}
impl core::fmt::Display for AssemblyResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "loaded assembly resonance refused: {self:?}")
    }
}
impl std::error::Error for AssemblyResonanceError {}

/// One numerical pole, including relaxing, unstable and port-dark modes.
/// It is not necessarily a peak in a particular physical driving response.
#[derive(Clone, Debug)]
pub struct AssemblyResonancePole {
    /// x(t)=v*exp(s*t); Re(s)<0 means decay, with no clipping or loss floor.
    pub laplace_pole_per_s: C64,
    /// omega=i*s under exp(-i*omega*t); a decaying pole has Im(omega)<0.
    /// Both conjugate Laplace poles remain present, with opposite Re(omega).
    pub angular_frequency_rad_s: C64,
    /// Full original [q,p,load] chart, p=M*v, including every real MIMO memory
    /// coordinate (pairs use the original Re/-Im chart). Maximum component is
    /// one with fixed phase; this numerical normalization supplies no energy.
    pub state: Vec<C64>,
    pub scaled_relative_residual: f64,
    pub kinematic_relative_residual: f64,
    pub momentum_relative_residual: f64,
    pub load_relative_residual: f64,
    /// Only |Im(s)| lies in the declared real-frequency source band. This is
    /// not a physical certification of analytic continuation off that axis.
    pub oscillation_frequency_in_source_band: bool,
}
impl AssemblyResonancePole {
    /// Positive oscillation frequency of this pole, NOT a forced-response peak.
    pub fn frequency_hz(&self) -> f64 {
        self.laplace_pole_per_s.im.abs() / core::f64::consts::TAU
    }
    /// Signed decay rate; a growing pole returns a negative number.
    pub fn decay_rate_per_s(&self) -> f64 {
        -self.laplace_pole_per_s.re
    }
    /// Cycle-decay convention |Im(s)|/(-2 Re(s)), only for decaying oscillation.
    /// No finite Q is invented for an undamped, growing or purely relaxing pole.
    /// This does not infer a -3 dB linewidth for a coupled/nonnormal transfer.
    pub fn quality_factor(&self) -> Option<f64> {
        let s = self.laplace_pole_per_s;
        if s.re >= 0. || s.im == 0. {
            return None;
        }
        let q = 0.5 * (s.im.abs() / (-s.re));
        q.is_finite().then_some(q)
    }
}

/// Complete numerical state spectrum, in native canonical Laplace (re,im)
/// order. Its length is exactly 2*retained_coordinates + original_load_states.
/// Neither weak coupling nor band membership removes a root or its multiplicity.
#[derive(Clone, Debug)]
pub struct AssemblyResonanceSpectrum {
    pub poles: Vec<AssemblyResonancePole>,
    pub retained_coordinates: usize,
    pub original_load_states: usize,
    pub source_band_rad_s: [f64; 2],
    pub continuous_source_class: CertificateClass,
    pub sampled_source_class: CertificateClass,
    /// Numerical basis resolution in the declared scaled coordinates. Neither
    /// a small backward residual nor this value is an eigenvalue enclosure.
    pub eigenbasis_reciprocal_condition: f64,
}

impl PreparedAssemblyMimoTransient<'_, '_, '_, '_> {
    /// Analyze the exact complete generator already used by this transient.
    /// This neither advances nor resets structural/load history, and its poles
    /// do not depend on the chosen transient sample clock. The existing source
    /// band and report classes are retained, without a new passivity claim.
    ///
    /// The base coordinate scaling uses sqrt(diag(K)),sqrt(diag(M^-1)); load
    /// coordinates keep the existing exact owner's numerical scales. Every
    /// resulting eigenstate is returned to ORIGINAL coordinates and its rows
    /// are independently checked against the retained unscaled generator.
    ///
    /// The shared L1 owner refuses an unresolved/defective eigenvector basis.
    /// Success concerns this complete finite generator, not omitted structural
    /// modes, continuum completeness, source accuracy outside its band, or an
    /// interval bound on small decay rates. Physical forced peaks/antiresonances
    /// still use the existing source-bound harmonic response.
    pub fn resonance_spectrum(
        &self,
        controls: AssemblyResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<AssemblyResonanceSpectrum, AssemblyResonanceError> {
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        let n = self.structure.vibration().mode_count();
        let mut scales = Vec::with_capacity(2 * n + self.joint.extra_scales.len());
        for i in 0..n {
            scales.push(det::sqrt(self.structure.stiffness_matrix()[i * n + i]));
        }
        for i in 0..n {
            scales.push(det::sqrt(self.structure.inverse_mass_matrix()[i * n + i]));
        }
        scales.extend_from_slice(&self.joint.extra_scales);
        let poles = analyze_generator(
            &self.joint.generator_per_s,
            n,
            &scales,
            self.source.limits().band_rad_s,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(AssemblyResonanceSpectrum {
            poles: poles.0,
            retained_coordinates: n,
            original_load_states: scales.len() - 2 * n,
            source_band_rad_s: self.source.limits().band_rad_s,
            continuous_source_class: self.source.reports().0.class,
            sampled_source_class: self.source.reports().1.class,
            eigenbasis_reciprocal_condition: poles.1,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AssemblyResonanceError> {
    if cancelled() {
        Err(AssemblyResonanceError::Cancelled)
    } else {
        Ok(())
    }
}
fn numerical(error: DenseSpectrumError) -> AssemblyResonanceError {
    if error == DenseSpectrumError::Cancelled {
        AssemblyResonanceError::Cancelled
    } else {
        AssemblyResonanceError::Numerical(error)
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn analyze_generator(
    generator: &[f64],
    modes: usize,
    scales: &[f64],
    band: [f64; 2],
    controls: AssemblyResonanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<AssemblyResonancePole>, f64), AssemblyResonanceError> {
    poll(cancelled)?;
    let d = scales.len();
    let size = d as u128;
    let rate = controls.reference_angular_frequency_rad_s;
    if modes == 0
        || 2 * modes > d
        || d.checked_mul(d) != Some(generator.len())
        || d > controls.numerical.maximum_dimension
        || !rate.is_finite()
        || rate <= 0.
        || !controls.maximum_original_relative_residual.is_finite()
        || controls.maximum_original_relative_residual <= 0.
        || controls.maximum_original_relative_residual >= 1.
        || 64 * size * size + 512 * size > controls.maximum_payload_bytes as u128
        || 32 * size.pow(3) > controls.maximum_back_transform_terms as u128
        || scales.iter().any(|v| !v.is_finite() || *v <= 0.)
    {
        return Err(AssemblyResonanceError::Invalid(
            "complete spectrum, valid scales/residuals and sufficient conversion budgets required",
        ));
    }
    let mut matrix = vec![C64::ZERO; d * d];
    let mut original = vec![C64::ZERO; d * d];
    for i in 0..d {
        poll(cancelled)?;
        for j in 0..d {
            let a = generator[i * d + j];
            let value = ((a / rate) * scales[i]) / scales[j];
            if !a.is_finite() || !value.is_finite() || (a != 0. && value == 0.) {
                return Err(AssemblyResonanceError::Invalid(
                    "original generator scaling is unrepresentable",
                ));
            }
            matrix[i * d + j] = C64::from_re(value);
            original[i * d + j] = C64::from_re(a);
        }
    }
    let spectrum = nonsymmetric_spectrum(&matrix, d, controls.numerical, &mut *cancelled)
        .map_err(numerical)?;
    let mut result = Vec::with_capacity(d);
    for (index, ((&scaled, vector), &scaled_residual)) in spectrum
        .values
        .iter()
        .zip(&spectrum.right_vectors)
        .zip(&spectrum.relative_residuals)
        .enumerate()
    {
        poll(cancelled)?;
        let value = scaled.scale(rate);
        if !finite(value)
            || (scaled.re != 0. && value.re == 0.)
            || (scaled.im != 0. && value.im == 0.)
        {
            return Err(AssemblyResonanceError::Invalid(
                "physical pole rate is unrepresentable",
            ));
        }
        let mut state = Vec::with_capacity(d);
        for (&v, &scale) in vector.iter().zip(scales) {
            let x = C64::new(v.re / scale, v.im / scale);
            if !finite(x) || (v.re != 0. && x.re == 0.) || (v.im != 0. && x.im == 0.) {
                return Err(AssemblyResonanceError::Invalid(
                    "original eigenstate is unrepresentable",
                ));
            }
            state.push(x);
        }
        let pivot = (0..d)
            .max_by(|&a, &b| state[a].abs().total_cmp(&state[b].abs()))
            .unwrap();
        let magnitude = state[pivot].abs();
        if !magnitude.is_finite() || magnitude == 0. {
            return Err(AssemblyResonanceError::Invalid(
                "original eigenstate has no finite amplitude",
            ));
        }
        let phase = C64::new(state[pivot].re / magnitude, state[pivot].im / magnitude);
        for x in &mut state {
            *x = C64::new(x.re / magnitude, x.im / magnitude) / phase;
            if !finite(*x) {
                return Err(AssemblyResonanceError::Invalid(
                    "original eigenstate phase normalization overflow",
                ));
            }
        }
        let rows = eigenpair_row_residuals(&original, d, value, &state, &mut *cancelled)
            .map_err(numerical)?;
        let worst = |range: core::ops::Range<usize>| rows[range].iter().copied().fold(0., f64::max);
        let (kinematic, momentum, load) = (
            worst(0..modes),
            worst(modes..2 * modes),
            worst(2 * modes..d),
        );
        if kinematic.max(momentum).max(load) > controls.maximum_original_relative_residual {
            return Err(AssemblyResonanceError::OriginalResidual {
                index,
                kinematic,
                momentum,
                load,
            });
        }
        let frequency = value.im.abs();
        result.push(AssemblyResonancePole {
            laplace_pole_per_s: value,
            angular_frequency_rad_s: C64::new(-value.im, value.re),
            state,
            scaled_relative_residual: scaled_residual,
            kinematic_relative_residual: kinematic,
            momentum_relative_residual: momentum,
            load_relative_residual: load,
            oscillation_frequency_in_source_band: frequency >= band[0] && frequency <= band[1],
        });
    }
    poll(cancelled)?;
    Ok((result, spectrum.eigenbasis_reciprocal_condition))
}

#[cfg(test)]
mod tests;
