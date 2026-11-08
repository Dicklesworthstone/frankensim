//! Finite-volume momentum and kinetic energy of a represented Mann realization.
//!
//! A spatial average of velocity does not determine the average of its square.
//! For centered cosine/sine modal vectors C_m/S_m and uniform box gain G(k),
//! the central second moment uses
//! A_mn = (G(k_m-k_n)+G(k_m+k_n))/2 - G(k_m)G(k_n),
//! B_mn = (G(k_m-k_n)-G(k_m+k_n))/2.
//! Mixed sine/cosine terms vanish by symmetry. All mode pairs are retained;
//! this is the instantaneous spatial stress, not an ensemble closure.
//! Small-footprint coefficients are integrated as central moments directly,
//! retaining the quadratic sine and quartic cosine terms without subtracting
//! rounded unit gains. Sampling retains the same prepared pair coefficients.

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

// For s uniform on [-1,1], return Cov(cos(a s),cos(b s)), E[sin(a s)sin(b s)]
// and E[cos(a s)] E[cos(b s)]. Keeping the first quantity CENTRAL is essential:
// it is O(a² b²), below the rounding error of unit gains on small footprints.
fn edge_moments(a: f64, b: f64) -> Result<[f64; 3], Refusal> {
    if !a.is_finite() || !b.is_finite() {
        return Err(error(
            "mann-filter-invalid",
            "quadratic filter phase is not representable",
        ));
    }
    let product = sinc(a) * sinc(b);
    if a == 0.0 || b == 0.0 {
        return Ok([0.0, 0.0, product]);
    }
    if a.abs().max(b.abs()) <= 0.5 {
        // Taylor coefficients times exact uniform moments. In particular,
        // Cov(s^(2m),s^(2n)) = 4mn/[(2m+2n+1)(2m+1)(2n+1)].
        // Eight terms per factor leave < 3e-20 absolute omitted sine terms
        // at the branch boundary; small phases retain their own scale.
        const TERMS: usize = 8;
        let (a2, b2) = (a * a, b * b);
        let mut cosine = [[0.0; 2]; TERMS];
        let mut sine = [[0.0; 2]; TERMS];
        let mut even = [-0.5 * a2, -0.5 * b2];
        let mut odd = [a, b];
        for n in 0..TERMS {
            cosine[n] = even;
            sine[n] = odd;
            let degree = (2 * (n + 1)) as f64;
            even[0] *= -a2 / ((degree + 1.0) * (degree + 2.0));
            even[1] *= -b2 / ((degree + 1.0) * (degree + 2.0));
            odd[0] *= -a2 / (degree * (degree + 1.0));
            odd[1] *= -b2 / (degree * (degree + 1.0));
        }
        let mut central = 0.0;
        let mut odd_product = 0.0;
        for m in 0..TERMS {
            for n in 0..TERMS {
                let dm = (2 * (m + 1)) as f64;
                let dn = (2 * (n + 1)) as f64;
                central += cosine[m][0] * cosine[n][1] * (dm * dn)
                    / ((dm + dn + 1.0) * (dm + 1.0) * (dn + 1.0));
                odd_product += sine[m][0] * sine[n][1] / (dm + dn - 1.0);
            }
        }
        Ok([central, odd_product, product])
    } else {
        let (minus, plus) = (a - b, a + b);
        if !minus.is_finite() || !plus.is_finite() {
            return Err(error(
                "mann-filter-invalid",
                "quadratic filter phase is not representable",
            ));
        }
        let (minus, plus) = (sinc(minus), sinc(plus));
        Ok([
            0.5 * (minus + plus) - product,
            0.5 * (minus - plus),
            product,
        ])
    }
}

fn pair_moments(
    first: [f64; 3],
    second: [f64; 3],
    half_edges: [[f64; 3]; 3],
) -> Result<[f64; 2], Refusal> {
    let mut central = 0.0;
    let mut odd_product = 0.0;
    let mut product = 1.0;
    for edge in half_edges {
        let [a, b, p] = edge_moments(super::dot(first, edge), super::dot(second, edge))?;
        // Angle addition over independent symmetric box coordinates. Expand
        // the central term before evaluating: (A+P)(a+p)+Bb-Pp would lose
        // the very stress that this observation is supposed to retain.
        let next_central = central * p + product * a + central * a + odd_product * b;
        odd_product = odd_product * (a + p) + (central + product) * b;
        central = next_central;
        product *= p;
    }
    Ok([central, odd_product])
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
        for m in 0..n {
            for j in m..n {
                poll(&mut cancelled)?;
                coefficients.push(pair_moments(
                    model.modes[m].k,
                    model.modes[j].k,
                    half_edges_m,
                )?);
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
