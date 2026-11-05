//! Finite-volume momentum and kinetic energy of a represented Mann realization.
//!
//! A spatial average of velocity does not determine the average of its square.
//! For centered cosine/sine modal vectors C_m/S_m and uniform box gain G(k),
//! the central second moment uses
//! A_mn = (G(k_m-k_n)+G(k_m+k_n))/2 - G(k_m)G(k_n),
//! B_mn = (G(k_m-k_n)-G(k_m+k_n))/2.
//! Mixed sine/cosine terms vanish by symmetry. All mode pairs are retained;
//! this is the instantaneous spatial stress, not an ensemble closure.
//! The exact product-to-sum observation is limited by floating-point roundoff,
//! especially when the footprint collapses and its stress approaches zero.

use super::{MannBoxFilter, error, poll, sinc};
use crate::mann_field::{MannSpectralModel, finite_vector, sample_modes};
use crate::{FlowSample, MAX_MODES, Refusal};
use fs_math::det;

/// Spatial mean and central second moment of one accepted realization.
#[derive(Clone, Debug, PartialEq)]
pub struct MannFilteredMoments {
    /// The original exact box-averaged velocity and analytic gradient.
    pub mean: FlowSample,
    /// `overline(u_i u_j) - overline(u_i) overline(u_j)` [m²/s²].
    /// Symmetric, positive semidefinite in exact arithmetic. Small roundoff
    /// eigenvalues are not clipped or replaced by an isotropic approximation.
    pub subfilter_stress_m2_s2: [[f64; 3]; 3],
}

impl MannFilteredMoments {
    /// Exact represented convective momentum flux per unit density [m²/s²].
    /// Pressure, density fluctuations and any mean flow are separate inputs.
    #[must_use]
    pub fn momentum_flux_per_density(&self) -> [[f64; 3]; 3] {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                self.mean.u[i] * self.mean.u[j] + self.subfilter_stress_m2_s2[i][j]
            })
        })
    }

    /// Kinetic energy per unit mass contained in the spatial fluctuations.
    #[must_use]
    pub fn subfilter_kinetic_energy(&self) -> f64 {
        0.5 * (0..3)
            .map(|i| self.subfilter_stress_m2_s2[i][i])
            .sum::<f64>()
    }
}

/// Prepared exact quadratic observation, bound to one spectrum and footprint.
/// Preparation stores two coefficients per unordered pair of modes. Sampling
/// needs O(modes²) work, no heap allocations, and 12 KiB bounded stack scratch.
/// Use the original `MannBoxFilter` when only the linear mean is needed.
#[derive(Clone, Debug, PartialEq)]
pub struct MannBoxMoments {
    filter: MannBoxFilter,
    pair_coefficients: Vec<[f64; 2]>,
}

impl MannBoxMoments {
    /// Prepare mean plus quadratic observations with explicit storage and
    /// mode-pair work limits. A point observation stores no coefficients and
    /// returns an identically zero subfilter stress. Cancelled preparation
    /// publishes no observation and never changes a temporal field.
    pub fn prepare(
        model: &MannSpectralModel,
        half_edges_m: [[f64; 3]; 3],
        maximum_workspace_bytes: usize,
        maximum_mode_pairs: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        let filter =
            MannBoxFilter::prepare(model, half_edges_m, maximum_workspace_bytes, &mut cancelled)?;
        let n = if filter.gains.is_empty() {
            0
        } else {
            model.mode_count()
        };
        let pairs = n * (n + 1) / 2;
        let required = filter.workspace_bytes() + pairs * core::mem::size_of::<[f64; 2]>();
        if pairs > maximum_mode_pairs || required > maximum_workspace_bytes {
            return Err(error(
                "mann-moments-budget",
                "quadratic observation exceeds its pair or storage budget",
            ));
        }
        let mut coefficients = Vec::new();
        coefficients.try_reserve_exact(pairs).map_err(|_| {
            error(
                "mann-filter-allocation",
                "unable to reserve quadratic observation coefficients",
            )
        })?;
        if filter.workspace_bytes() + coefficients.capacity() * core::mem::size_of::<[f64; 2]>()
            > maximum_workspace_bytes
        {
            return Err(error(
                "mann-moments-budget",
                "allocated quadratic coefficient capacity exceeds the byte budget",
            ));
        }
        let gain = |k: [f64; 3]| -> Result<f64, Refusal> {
            let mut value = 1.0;
            for edge in half_edges_m {
                let phase = super::dot(k, edge);
                if !phase.is_finite() {
                    return Err(error(
                        "mann-filter-invalid",
                        "quadratic filter phase is not representable",
                    ));
                }
                value *= sinc(phase);
            }
            Ok(value)
        };
        for m in 0..n {
            for j in m..n {
                poll(&mut cancelled)?;
                let minus = gain(core::array::from_fn(|i| {
                    model.modes[m].k[i] - model.modes[j].k[i]
                }))?;
                let plus = gain(core::array::from_fn(|i| {
                    model.modes[m].k[i] + model.modes[j].k[i]
                }))?;
                coefficients.push([
                    0.5 * (minus + plus) - filter.gains[m] * filter.gains[j],
                    0.5 * (minus - plus),
                ]);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            filter,
            pair_coefficients: coefficients,
        })
    }

    /// Allocated numerical storage including the original linear filter.
    #[must_use]
    pub fn workspace_bytes(&self) -> usize {
        self.filter.workspace_bytes()
            + self.pair_coefficients.capacity() * core::mem::size_of::<[f64; 2]>()
    }

    /// Unordered mode pairs evaluated by each non-point observation.
    #[must_use]
    pub fn mode_pairs(&self) -> usize {
        self.pair_coefficients.len()
    }
}

pub(crate) fn sample(
    model: &MannSpectralModel,
    time_s: f64,
    amplitudes: [&[f64]; 4],
    point_m: [f64; 3],
    observation: &MannBoxMoments,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MannFilteredMoments, Refusal> {
    poll(&mut cancelled)?;
    finite_vector(point_m)?;
    let gains = observation.filter.gains(model)?;
    let mean = sample_modes(model, time_s, amplitudes, point_m, gains)?;
    let mut stress = [[0.0; 3]; 3];
    if gains.is_some() {
        let position = core::array::from_fn(|i| point_m[i] - model.advection_mps[i] * time_s);
        let mut vectors = [[[0.0; 3]; 2]; MAX_MODES];
        for (m, mode) in model.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            let phase = super::dot(mode.k, position);
            let (s, c) = (det::sin(phase), det::cos(phase));
            for i in 0..3 {
                vectors[m][0][i] = mode.polarization[0][i]
                    * (amplitudes[0][m] * c + amplitudes[1][m] * s)
                    + mode.polarization[1][i] * (amplitudes[2][m] * c + amplitudes[3][m] * s);
                vectors[m][1][i] = mode.polarization[0][i]
                    * (-amplitudes[0][m] * s + amplitudes[1][m] * c)
                    + mode.polarization[1][i] * (-amplitudes[2][m] * s + amplitudes[3][m] * c);
            }
        }
        let mut pairs = observation.pair_coefficients.iter();
        for m in 0..model.mode_count() {
            for j in m..model.mode_count() {
                poll(&mut cancelled)?;
                let coefficients = pairs.next().expect("prepared pair count");
                for a in 0..3 {
                    for b in a..3 {
                        for q in 0..2 {
                            let mut product = vectors[m][q][a] * vectors[j][q][b];
                            if m != j {
                                product += vectors[j][q][a] * vectors[m][q][b];
                            }
                            stress[a][b] += coefficients[q] * product;
                        }
                    }
                }
            }
        }
        for a in 0..3 {
            for b in 0..a {
                stress[a][b] = stress[b][a];
            }
        }
    }
    if !stress.iter().flatten().all(|v| v.is_finite()) {
        return Err(error(
            "mann-filter-invalid",
            "quadratic spatial moment is not finite",
        ));
    }
    poll(&mut cancelled)?;
    Ok(MannFilteredMoments {
        mean,
        subfilter_stress_m2_s2: stress,
    })
}

#[cfg(test)]
mod tests;
