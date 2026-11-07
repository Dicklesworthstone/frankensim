//! Analytic frequency dispersion of the original full tensor-Prony law.

use super::{Mandel6, TensorProny, ViscoError};
use crate::visco::prony_spectrum::PronySpectrumWeights;
use fs_math::dd::Dd;

mod plane_stress;
mod normal_stress;
pub use normal_stress::NormalStressHarmonic;
pub(crate) use plane_stress::projection::Sum as HarmonicWorkSum;
mod strain;
pub use plane_stress::PlaneStressHarmonic;
mod response;
pub use response::{TensorPronyHarmonic, TensorPronyHarmonicBranch, TensorPronyHarmonicPoint};

/// Physical tensor derivatives with respect to angular frequency [rad/s].
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyFrequencyDerivative {
    /// dC'/domega [Pa s], retaining all anisotropic/mutual entries.
    pub storage_pa_s: Mandel6,
    /// dC''/domega [Pa s], including the original Kelvin viscosity.
    pub loss_pa_s: Mandel6,
    /// d(C''/omega)/domega [Pa s²], evaluated directly from each branch.
    /// The continuous DC limit is zero. No subtraction of two nearly equal
    /// loss/omega terms is used; the Kelvin contribution is exactly zero.
    pub effective_viscosity_pa_s2: Mandel6,
}

pub(super) fn admit_frequency(law: &TensorProny, omega: f64) -> Result<(), ViscoError> {
    if !omega.is_finite() {
        return Err(ViscoError::Parameters {
            what: "angular frequency omega must be finite",
        });
    }
    if omega < law.band_rad_s.0 || omega > law.band_rad_s.1 {
        return Err(ViscoError::OutOfBand {
            omega,
            band: law.band_rad_s,
        });
    }
    Ok(())
}

// A modulus-only consumer must not request a potentially unrepresentable
// effective-viscosity slope. Both queries share the exact branch aggregation.
pub(super) fn modulus_frequency_derivative(
    law: &TensorProny,
    omega: f64,
) -> Result<(Mandel6, Mandel6), ViscoError> {
    let [storage, loss] = derivative_tensors(law, omega, PronySpectrumWeights::modulus_derivative)?;
    Ok((storage, loss))
}

fn derivative_tensors<const N: usize>(
    law: &TensorProny,
    omega: f64,
    branch_values: impl Fn(&PronySpectrumWeights, f64) -> Option<[f64; N]>,
) -> Result<[Mandel6; N], ViscoError> {
    admit_frequency(law, omega)?;
    let invalid = || ViscoError::Parameters {
        what: "tensor-Prony physical frequency derivative is not representable",
    };
    // Storage/loss lead each supported spectral query; viscosity is optional.
    let mut sums = [[[Dd::ZERO; 6]; 6]; N];
    if let Some(kelvin) = &law.kv_pa_s {
        for i in 0..6 {
            for j in 0..=i {
                sums[1][i][j] = Dd::from_f64(kelvin[i][j]);
            }
        }
    }
    for branch in &law.branches {
        let weights = PronySpectrumWeights::new(omega, branch.tau_s).ok_or_else(invalid)?;
        for i in 0..6 {
            for j in 0..=i {
                let values = branch_values(&weights, branch.c_pa[i][j]).ok_or_else(invalid)?;
                for k in 0..N {
                    sums[k][i][j] = sums[k][i][j] + Dd::from_f64(values[k]);
                }
            }
        }
    }
    let values: [Mandel6; N] = sums.map(|matrix| {
        core::array::from_fn(|i| core::array::from_fn(|j| matrix[i.max(j)][i.min(j)].to_f64()))
    });
    if values
        .iter()
        .flatten()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(invalid());
    }
    Ok(values)
}

impl TensorProny {
    /// Admit the maximum physical Mandel strain over an entire harmonic cycle.
    ///
    /// `strain` is a peak phasor, not an RMS amplitude. Both harmonic sign
    /// conventions describe the same ellipse and have the same peak norm.
    /// The largest eigenvalue of its scaled 2x2 Gram matrix gives the exact
    /// cycle extremum without sampling phases or summing component maxima.
    /// The original tensor material's strain limit is retained. This admits
    /// one spatial point; it does not certify unsampled spatial variation.
    pub fn validate_harmonic_mandel_strain(
        &self,
        strain: &[fs_math::c64::C64; 6],
    ) -> Result<(), ViscoError> {
        admit_harmonic_strain(strain, self.strain_limit)
    }

    /// Differentiate the original storage, loss and effective harmonic viscosity.
    ///
    /// The same constitutive frequency band as `complex_modulus` applies. The
    /// algebraic law supplies one-sided slopes at admitted band endpoints.
    /// Full tensor coefficients, branch rates and Kelvin viscosity remain in
    /// their original frame. No modal fit or frequency perturbation is used.
    /// Products retain binary scale until their physical coefficient is applied;
    /// compensated sums retain signed cancellation among relaxation slopes.
    ///
    /// This is frequency-domain dispersion only. C''/omega is not a replacement
    /// for the material's hereditary state during time integration. A finite
    /// primal query remains available if the requested tangent is unrepresentable.
    pub fn complex_modulus_frequency_derivative(
        &self,
        omega: f64,
    ) -> Result<TensorPronyFrequencyDerivative, ViscoError> {
        let [storage_pa_s, loss_pa_s, effective_viscosity_pa_s2] =
            derivative_tensors(self, omega, PronySpectrumWeights::derivative)?;
        Ok(TensorPronyFrequencyDerivative {
            storage_pa_s,
            loss_pa_s,
            effective_viscosity_pa_s2,
        })
    }
}

fn admit_harmonic_strain(strain: &[fs_math::c64::C64; 6], limit: f64) -> Result<(), ViscoError> {
    match strain::peak_norm(strain) {
        Some(peak) if peak <= limit => Ok(()),
        _ => Err(ViscoError::Parameters {
            what: "harmonic strain exceeds the original tensor material small-strain domain",
        }),
    }
}
