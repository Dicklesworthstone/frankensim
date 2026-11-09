//! Positive, component-balanced importance quadrature of the original tensor.
//!
//! Coordinates are s=log(|k|L), mu=kz/|k| and phi in [0,pi/2]. A bounded
//! table approximates each diagonal energy density only to CHOOSE the proposal
//! distribution. Normalizing those three proposal densities balances their
//! sampling effort; it never rescales any physical tensor or realization.
//! A 5% uniform-density mixture keeps positive support throughout the band.
//!
//! Stratified/shifted unit points invert its three conditional CDFs. If their
//! joint density is p(s,mu,phi), one ky-reflected mode pair has weight
//! 2 |k|^3/(M p) per mode, where M is the number of pairs. Reflection supplies
//! the other ky half-plane and real sine/cosine quadratures supply the -k
//! conjugate; together these cover the full sphere exactly in the integration
//! measure. Every mode is factored from the unchanged Phi(k). There is no
//! empirical covariance matching or self-normalized importance weighting.

use super::*;

const RADIAL: usize = 32;
const POLAR: usize = 32;
const AZIMUTH: usize = 16;
const CELLS: usize = RADIAL * POLAR * AZIMUTH;
const UNIFORM_FRACTION: f64 = 0.05;

fn error(message: &str) -> Refusal {
    refuse(
        "mann-importance-invalid",
        message.into(),
        "use an even mode count in 2..=256 and a representable finite kL band",
    )
}

fn storage<T>(count: usize) -> Result<Vec<T>, Refusal> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| {
        error("unable to reserve bounded importance-quadrature preparation storage")
    })?;
    Ok(values)
}

// Private preparation can be reused by numerical batteries. The public
// constructor drops it once the actual mode geometry has been generated.
pub(super) struct Prepared {
    parameters: MannParams,
    quadrature: MannQuadrature,
    cells: Vec<f64>,
    rows: Vec<f64>,
    radial_cdf: [f64; RADIAL],
    mu_edges: [f64; POLAR + 1],
    log_min: f64,
    delta_log: f64,
    delta_phi: f64,
}

impl Prepared {
    pub(super) fn new(parameters: MannParams, quadrature: MannQuadrature) -> Result<Self, Refusal> {
        admit_model_inputs(&parameters, quadrature, [0.0; 3])?;
        if quadrature.n_modes < 2 || !quadrature.n_modes.is_multiple_of(2) {
            return Err(error(
                "balanced quadrature requires independent reflected mode pairs",
            ));
        }
        let log_min = det::ln(quadrature.min_kl);
        let delta_log = (det::ln(quadrature.max_kl) - log_min) / RADIAL as f64;
        let delta_phi = core::f64::consts::FRAC_PI_2 / AZIMUTH as f64;
        if !delta_log.is_finite() || delta_log <= 0.0 {
            return Err(error("logarithmic importance band is not representable"));
        }
        // Nonuniform mu cells resolve the poles. Interpolation and the measure
        // are in mu itself, retaining uniform solid angle in the isotropic limit.
        let mu_edges: [f64; POLAR + 1] =
            core::array::from_fn(|i| -det::cos(core::f64::consts::PI * i as f64 / POLAR as f64));
        let unit = MannParams {
            alpha_eps23: 1.0,
            length_m: 1.0,
            gamma: parameters.gamma,
        };
        let mut raw: Vec<[f64; 3]> = storage(CELLS)?;
        let mut totals = [0.0; 3];
        for radial in 0..RADIAL {
            let radius = det::exp(log_min + (radial as f64 + 0.5) * delta_log);
            for polar in 0..POLAR {
                let mu = 0.5 * (mu_edges[polar] + mu_edges[polar + 1]);
                let width = mu_edges[polar + 1] - mu_edges[polar];
                let planar = det::sqrt((1.0 - mu * mu).max(0.0));
                for azimuth in 0..AZIMUTH {
                    let (sine, cosine) = det::sin_cos((azimuth as f64 + 0.5) * delta_phi);
                    let tensor = mann_tensor(
                        &unit,
                        [
                            radius * planar * cosine,
                            radius * planar * sine,
                            radius * mu,
                        ],
                    )?;
                    let mass: [f64; 3] =
                        core::array::from_fn(|i| tensor[i][i] * radius * radius * radius * width);
                    if mass.iter().any(|&v| !v.is_finite() || v < 0.0) {
                        return Err(error("proposal component energy is not representable"));
                    }
                    for i in 0..3 {
                        totals[i] += mass[i];
                    }
                    raw.push(mass);
                }
            }
        }
        if totals.iter().any(|&v| !v.is_finite() || v <= 0.0) {
            return Err(error(
                "proposal component integral is not positive and finite",
            ));
        }
        let mut cells = storage(CELLS)?;
        for (index, mass) in raw.iter().enumerate() {
            let polar = (index / AZIMUTH) % POLAR;
            let width = mu_edges[polar + 1] - mu_edges[polar];
            let balanced: f64 = (0..3).map(|i| mass[i] / totals[i]).sum();
            // Each normalized component integrates to one in this proposal;
            // the uniform part has the same total mass, three. It is not a
            // floor on the physical spectrum or a restored spectral tail.
            let uniform = 3.0 * width / (2 * RADIAL * AZIMUTH) as f64;
            cells.push((1.0 - UNIFORM_FRACTION) * balanced + UNIFORM_FRACTION * uniform);
        }
        drop(raw);
        let mut rows = storage(RADIAL * POLAR)?;
        let mut radial_cdf = [0.0; RADIAL];
        let mut cumulative = 0.0;
        for (radial, radial_total) in radial_cdf.iter_mut().enumerate() {
            let mut row_total = 0.0;
            for polar in 0..POLAR {
                let mut cell_total = 0.0;
                for azimuth in 0..AZIMUTH {
                    let index = (radial * POLAR + polar) * AZIMUTH + azimuth;
                    let next = cell_total + cells[index];
                    if !next.is_finite() || next <= cell_total {
                        return Err(error("proposal cell probability is not representable"));
                    }
                    cell_total = next;
                    cells[index] = next;
                }
                let next = row_total + cell_total;
                if !next.is_finite() || next <= row_total {
                    return Err(error("proposal polar probability is not representable"));
                }
                row_total = next;
                rows.push(next);
            }
            let next = cumulative + row_total;
            if !next.is_finite() || next <= cumulative {
                return Err(error("proposal radial probability is not representable"));
            }
            cumulative = next;
            *radial_total = next;
        }
        Ok(Self {
            parameters,
            quadrature,
            cells,
            rows,
            radial_cdf,
            mu_edges,
            log_min,
            delta_log,
            delta_phi,
        })
    }

    // Returns physical k and its full-sphere measure per reflected mode.
    fn draw(&self, unit: [f64; 3]) -> Result<([f64; 3], f64), Refusal> {
        let (radial, r_fraction, p_r) = invert(&self.radial_cdf, unit[0])?;
        let (polar, mu_fraction, p_mu) =
            invert(&self.rows[radial * POLAR..(radial + 1) * POLAR], unit[1])?;
        let row = (radial * POLAR + polar) * AZIMUTH;
        let (azimuth, phi_fraction, p_phi) = invert(&self.cells[row..row + AZIMUTH], unit[2])?;
        let radius = det::exp(self.log_min + (radial as f64 + r_fraction) * self.delta_log)
            / self.parameters.length_m;
        let width = self.mu_edges[polar + 1] - self.mu_edges[polar];
        let mu = self.mu_edges[polar] + mu_fraction * width;
        let planar = det::sqrt((1.0 - mu * mu).max(0.0));
        let (sine, cosine) = det::sin_cos((azimuth as f64 + phi_fraction) * self.delta_phi);
        let k = [
            radius * planar * cosine,
            radius * planar * sine,
            radius * mu,
        ];
        // Include all three conditional probabilities rather than assuming
        // floating-point prefix differences telescope exactly between rows.
        let jacobian = self.delta_log * width * self.delta_phi / p_r / p_mu / p_phi;
        let weight =
            2.0 * jacobian / (self.quadrature.n_modes / 2) as f64 * radius * radius * radius;
        if !k.iter().all(|v| v.is_finite()) || !weight.is_finite() || weight <= 0.0 {
            return Err(error(
                "importance wavevector or positive volume weight is not representable",
            ));
        }
        Ok((k, weight))
    }

    pub(super) fn build(
        &self,
        geometry_seed: u64,
        advection_mps: [f64; 3],
    ) -> Result<MannSpectralModel, Refusal> {
        admit_model_inputs(&self.parameters, self.quadrature, advection_mps)?;
        let mut stream = StreamKey {
            seed: geometry_seed,
            kernel: MANN_FIELD_KERNEL,
            tile: 0,
        }
        .stream();
        let shifts = [stream.next_f64(), stream.next_f64(), stream.next_f64()];
        let mut modes = storage(self.quadrature.n_modes)?;
        let pairs = self.quadrature.n_modes / 2;
        for index in 0..pairs {
            let unit = [
                (index as f64 + shifts[0]) / pairs as f64,
                shifted_unit(radical_inverse(index + 1, 2), shifts[1]),
                shifted_unit(radical_inverse(index + 1, 3), shifts[2]),
            ];
            let (k, weight) = self.draw(unit)?;
            let phi = mann_tensor(&self.parameters, k)?;
            let polarization = polarizations(k, phi, weight)?;
            modes.push(SpectralMode {
                k,
                weight,
                polarization,
            });
            // Phi(kx,-ky,kz) = R Phi(k) R^T, R=diag(1,-1,1).
            // The mode index is distinct, so both OU owners automatically
            // retain independent stochastic lanes for the reflected partner.
            modes.push(SpectralMode {
                k: [k[0], -k[1], k[2]],
                weight,
                polarization: polarization.map(|p| [p[0], -p[1], p[2]]),
            });
        }
        Ok(MannSpectralModel {
            geometry_seed,
            parameters: self.parameters,
            quadrature: self.quadrature,
            quadrature_rule: MannQuadratureRule::BalancedEnergyV1,
            advection_mps,
            modes,
        })
    }
}

// Piecewise-uniform inverse CDF with its actual conditional cell probability.
fn invert(cdf: &[f64], unit: f64) -> Result<(usize, f64, f64), Refusal> {
    if !unit.is_finite() || !(0.0..=1.0).contains(&unit) {
        return Err(error("importance coordinate is outside the unit interval"));
    }
    let total = cdf[cdf.len() - 1];
    let position = unit * total;
    let index = cdf
        .partition_point(|&value| value <= position)
        .min(cdf.len() - 1);
    let left = if index == 0 { 0.0 } else { cdf[index - 1] };
    let mass = cdf[index] - left;
    let fraction = (position - left) / mass;
    let probability = mass / total;
    if !fraction.is_finite()
        || !(0.0..=1.0).contains(&fraction)
        || !probability.is_finite()
        || probability <= 0.0
    {
        return Err(error(
            "conditional importance probability is not representable",
        ));
    }
    Ok((index, fraction, probability))
}

#[cfg(test)]
mod tests;
