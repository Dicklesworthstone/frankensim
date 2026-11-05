//! Exact driven frequency tangents of the complete mechanical/acoustic/material
//! pencil. T*y'=f'-T'*y uses the SAME assembled operator and LU as T*y=f.
//! No frequency difference, acoustic inverse, damping fit or new eigenproblem.

use super::*;
use fs_modal::harmonic::resonance::HarmonicTransferSample;

/// Original response limits plus the complete additional tangent payload/work.
#[derive(Clone, Copy, Debug)]
pub struct BemHarmonicFrequencyOptions {
    pub response: BemHarmonicOptions,
    /// cases*(32*full_dimension+512) + 128*full_dimension bytes, including
    /// returned tangent fields and temporary original-coordinate contractions.
    pub max_additional_response_bytes: usize,
    /// 256*cases*full_dimension^2 conservative contraction/solve/audit units.
    /// The original response limits already charge the one assembly and LU.
    pub max_additional_work: u64,
}

/// Signed derivatives of cycle-average powers per rad/s. Geometry, medium,
/// physical maps and material parameters are fixed; their dispersion remains
/// in the differentiated states. Boundary work remains separate from heat.
#[derive(Clone, Copy, Debug)]
pub struct BemHarmonicPowerDerivative {
    pub input_w_per_rad_s: f64,
    pub kelvin_w_per_rad_s: f64,
    pub relaxation_w_per_rad_s: f64,
    pub radiated_w_per_rad_s: f64,
    pub balance_w_per_rad_s: f64,
    pub apparent_scale_w_per_rad_s: f64,
    pub relative_balance: f64,
}

/// Force-scaled response and exact derivative on the same physical branch.
/// Every derivative is per rad/s and keeps the original coordinate ordering.
/// The residuals audit T*y'+T'*y=f' in the unscaled original equations.
#[derive(Clone, Debug)]
pub struct BemHarmonicFrequencyResponse {
    pub response: BemHarmonicResponse,
    pub modal_displacement_derivative: Vec<C64>,
    pub modal_velocity_derivative: Vec<C64>,
    pub surface_pressure_derivative_pa: Vec<C64>,
    pub surface_normal_velocity_derivative_m_s: Vec<C64>,
    pub elastic_memory_derivative_sqrt_j: Vec<C64>,
    pub mechanical_derivative_residual: f64,
    pub boundary_derivative_residual: f64,
    pub memory_derivative_residual: f64,
    pub power_derivative: BemHarmonicPowerDerivative,
}

impl BemHarmonicFrequencyResponse {
    /// Fixed real displacement receiver in the retained physical mode basis.
    pub fn displacement_transfer(
        &self,
        weights: &[f64],
    ) -> Result<HarmonicTransferSample, BemModesError> {
        transfer(
            weights,
            &self.response.modal_displacement,
            &self.modal_displacement_derivative,
        )
    }
    /// Fixed real velocity receiver, including the derivative of -i*omega.
    pub fn velocity_transfer(
        &self,
        weights: &[f64],
    ) -> Result<HarmonicTransferSample, BemModesError> {
        transfer(
            weights,
            &self.response.modal_velocity,
            &self.modal_velocity_derivative,
        )
    }
    /// Fixed real surface-pressure receiver in original panel order. Area
    /// weighting or sensor normalization belongs in the supplied weights.
    pub fn surface_pressure_transfer(
        &self,
        weights: &[f64],
    ) -> Result<HarmonicTransferSample, BemModesError> {
        transfer(
            weights,
            &self.response.surface_pressure_pa,
            &self.surface_pressure_derivative_pa,
        )
    }
    /// Fixed real outward surface-velocity receiver in original panel order.
    pub fn surface_velocity_transfer(
        &self,
        weights: &[f64],
    ) -> Result<HarmonicTransferSample, BemModesError> {
        transfer(
            weights,
            &self.response.surface_normal_velocity_m_s,
            &self.surface_normal_velocity_derivative_m_s,
        )
    }
}

fn transfer(
    weights: &[f64],
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, BemModesError> {
    HarmonicTransferSample::project_real(weights, values, derivatives)
        .map_err(|_| bad("complete finite physical receiver and representable transfer required"))
}

impl BemLoadedModes {
    /// Differentiate complete driven responses with one assembly and one LU.
    /// Optional source slopes are in the same original structural coordinates;
    /// None holds force phasors fixed. Geometry, medium, basis and material
    /// parameters stay fixed, while acoustic and Maxwell dispersion remain.
    ///
    /// Original response and differentiated mechanical/boundary/memory rows
    /// must independently pass their declared numerical residual limits. The
    /// physical source, Kelvin, Maxwell and surface-power derivatives must
    /// close independently. Any refusal/cancellation discards the whole batch.
    pub fn harmonic_response_batch_with_frequency_derivative(
        &self,
        omega: f64,
        forces: &[&[C64]],
        force_derivatives: Option<&[&[C64]]>,
        options: &BemHarmonicFrequencyOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<BemHarmonicFrequencyResponse>, BemModesError> {
        poll(&mut cancelled)?;
        let (m, p, n) = (self.mode_count(), self.panel_count(), self.dimension());
        if force_derivatives.is_some_and(|slopes| {
            slopes.len() != forces.len()
                || slopes
                    .iter()
                    .any(|f| f.len() != m || f.iter().any(|&z| !finite(z)))
        }) {
            return Err(bad(
                "one complete finite force derivative per harmonic load case required",
            ));
        }
        let bytes = forces.len() as u128 * (32 * n as u128 + 512) + 128 * n as u128;
        let work = 256 * forces.len() as u128 * (n as u128).pow(2);
        if bytes > options.max_additional_response_bytes as u128
            || work > options.max_additional_work as u128
        {
            return Err(bad(
                "harmonic frequency tangent exceeds its additional payload or work budget",
            ));
        }
        let prepared =
            self.prepare_harmonic_batch(omega, forces, &options.response, &mut cancelled)?;
        let zero = vec![C64::ZERO; m];
        let mut results = Vec::with_capacity(forces.len());
        for (case, &force) in forces.iter().enumerate() {
            poll(&mut cancelled)?;
            let source_slope = force_derivatives.map_or(zero.as_slice(), |values| values[case]);
            let state = prepared.solve_rhs(force, &mut cancelled)?;
            let mut rhs = vec![C64::ZERO; n];
            for i in 0..n {
                poll(&mut cancelled)?;
                let mut sum = source_slope.get(i).copied().unwrap_or(C64::ZERO);
                let mut correction = C64::ZERO;
                for j in 0..n {
                    let term = -prepared.raw.value.derivative[i * n + j] * state[j] - correction;
                    let next = sum + term;
                    correction = (next - sum) - term;
                    sum = next;
                }
                rhs[i] = sum;
            }
            let mut tangent = prepared.solve_rhs(&rhs, &mut cancelled)?;
            let residuals = tangent_residuals(
                &prepared.raw,
                &state,
                &tangent,
                source_slope,
                m,
                p,
                self.settings.physical_residual_tolerance,
                &mut cancelled,
            )?;
            let velocity: Vec<_> = state[..m]
                .iter()
                .zip(&tangent[..m])
                .map(|(&x, &dx)| C64::new(0.0, -1.0) * x + C64::new(0.0, -omega) * dx)
                .collect();
            let normal = self.harmonic_normal_velocity(&velocity, &mut cancelled)?;
            let response = self.harmonic_response_from_state(
                omega,
                force,
                state,
                &prepared.raw,
                &options.response,
                &mut cancelled,
            )?;
            let power_derivative = self.harmonic_power_derivative(
                force,
                source_slope,
                &response,
                &velocity,
                &tangent[m..m + p],
                &normal,
                &tangent[m + p..],
                options.response.power_relative_tolerance,
                &mut cancelled,
            )?;
            let memory = tangent.split_off(m + p);
            let pressure = tangent.split_off(m);
            results.push(BemHarmonicFrequencyResponse {
                response,
                modal_displacement_derivative: tangent,
                modal_velocity_derivative: velocity,
                surface_pressure_derivative_pa: pressure,
                surface_normal_velocity_derivative_m_s: normal,
                elastic_memory_derivative_sqrt_j: memory,
                mechanical_derivative_residual: residuals[0],
                boundary_derivative_residual: residuals[1],
                memory_derivative_residual: residuals[2],
                power_derivative,
            });
        }
        poll(&mut cancelled)?;
        Ok(results)
    }

    #[allow(clippy::too_many_arguments)] // original physical work and its tangent
    fn harmonic_power_derivative(
        &self,
        force: &[C64],
        force_slope: &[C64],
        response: &BemHarmonicResponse,
        velocity_slope: &[C64],
        pressure_slope: &[C64],
        normal_slope: &[C64],
        memory_slope: &[C64],
        tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<BemHarmonicPowerDerivative, BemModesError> {
        // A common numerical amplitude scale for BOTH orders preserves every
        // product, even if the force crosses zero while its slope is nonzero.
        // This scaling is not a norm or a comparison of unlike physical units.
        let gain = force
            .iter()
            .chain(force_slope)
            .map(|z| z.abs())
            .fold(0.0, f64::max);
        let gain = if gain == 0.0 { 1.0 } else { gain };
        let normalize = |values: &[C64]| -> Result<Vec<C64>, BemModesError> {
            values
                .iter()
                .map(|&z| {
                    let value = divide(z, gain);
                    if !finite(value)
                        || (z.re != 0.0 && value.re == 0.0)
                        || (z.im != 0.0 && value.im == 0.0)
                    {
                        return Err(bad(
                            "harmonic power derivative normalization is not representable",
                        ));
                    }
                    Ok(value)
                })
                .collect()
        };
        let f = normalize(force)?;
        let df = normalize(force_slope)?;
        let v = normalize(&response.modal_velocity)?;
        let dv = normalize(velocity_slope)?;
        let p = normalize(&response.surface_pressure_pa)?;
        let dp = normalize(pressure_slope)?;
        let u = normalize(&response.surface_normal_velocity_m_s)?;
        let du = normalize(normal_slope)?;
        let e = normalize(&response.elastic_memory_sqrt_j)?;
        let de = normalize(memory_slope)?;
        let mut sums = [0.0_f64; 5];
        let mut corrections = [0.0_f64; 5];
        let mut add = |channel: usize, value: f64| {
            let term = value - corrections[channel];
            let next = sums[channel] + term;
            corrections[channel] = (next - sums[channel]) - term;
            sums[channel] = next;
        };
        let m = self.mode_count();
        for i in 0..m {
            poll(cancelled)?;
            add(0, 0.5 * (df[i] * v[i].conj() + f[i] * dv[i].conj()).re);
            add(
                4,
                0.5 * (df[i].abs() * v[i].abs() + f[i].abs() * dv[i].abs()),
            );
            for j in 0..m {
                let c = self.viscosity[i * m + j];
                add(1, 0.5 * (dv[i].conj() * v[j] + v[i].conj() * dv[j]).re * c);
                add(
                    4,
                    0.5 * (dv[i].abs() * v[j].abs() + v[i].abs() * dv[j].abs()) * c.abs(),
                );
            }
        }
        for ((&e, &de), &tau) in e.iter().zip(&de).zip(&self.memory.times) {
            poll(cancelled)?;
            add(2, (e.conj() * de).re / tau);
            add(4, e.abs() * de.abs() / tau);
        }
        for i in 0..self.panel_count() {
            poll(cancelled)?;
            let area = self.surface.areas()[i];
            add(
                3,
                0.5 * (dp[i] * u[i].conj() + p[i] * du[i].conj()).re * area,
            );
            add(
                4,
                0.5 * (dp[i].abs() * u[i].abs() + p[i].abs() * du[i].abs()) * area,
            );
        }
        let balance = sums[0] - sums[1] - sums[2] - sums[3];
        let relative = if sums[4] == 0.0 {
            balance.abs()
        } else {
            balance.abs() / sums[4]
        };
        if sums.iter().any(|x| !x.is_finite()) || !relative.is_finite() {
            return Err(bad("harmonic power derivative overflow"));
        }
        if relative > tolerance {
            return Err(BemModesError::PowerBalance {
                relative,
                tolerance,
            });
        }
        let restore = |value: f64| -> Result<f64, BemModesError> {
            let physical = (value * gain) * gain;
            if !physical.is_finite() || (value != 0.0 && physical == 0.0) {
                return Err(bad(
                    "physical harmonic power derivative is not representable",
                ));
            }
            Ok(physical)
        };
        Ok(BemHarmonicPowerDerivative {
            input_w_per_rad_s: restore(sums[0])?,
            kelvin_w_per_rad_s: restore(sums[1])?,
            relaxation_w_per_rad_s: restore(sums[2])?,
            radiated_w_per_rad_s: restore(sums[3])?,
            balance_w_per_rad_s: restore(balance)?,
            apparent_scale_w_per_rad_s: restore(sums[4])?,
            relative_balance: relative,
        })
    }
}

#[allow(clippy::too_many_arguments)] // two original states and physical row blocks
fn tangent_residuals(
    raw: &PhysicalPencil,
    state: &[C64],
    tangent: &[C64],
    source_slope: &[C64],
    modes: usize,
    panels: usize,
    tolerance: f64,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<[f64; 3], BemModesError> {
    let n = state.len();
    let norm = state.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let tangent_norm = tangent.iter().map(|v| v.abs()).fold(0.0, f64::max);
    // Same original-row rounding policy as the response audit, now for two
    // complex contractions. This is numerical evidence, not an interval proof.
    let units = (2 * n + 4) as f64;
    let gamma = units * (f64::EPSILON / 2.0) / (1.0 - units * (f64::EPSILON / 2.0));
    let mut residuals = [0.0_f64; 3];
    for i in 0..n {
        poll(cancelled)?;
        let source = source_slope.get(i).copied().unwrap_or(C64::ZERO);
        let mut sum = -source;
        let mut correction = C64::ZERO;
        let mut denominator = source.abs();
        let mut rounding_scale = source.abs();
        for j in 0..n {
            for term in [
                raw.value.matrix[i * n + j] * tangent[j],
                raw.value.derivative[i * n + j] * state[j],
            ] {
                let term = term - correction;
                let next = sum + term;
                correction = (next - sum) - term;
                sum = next;
            }
            let value_scale = raw.envelope[i * n + j];
            let derivative_scale = raw.value.derivative[i * n + j].abs();
            denominator += value_scale * tangent[j].abs() + derivative_scale * state[j].abs();
            rounding_scale += value_scale * tangent_norm + derivative_scale * norm;
        }
        if !finite(sum) || !denominator.is_finite() || !rounding_scale.is_finite() {
            return Err(bad(
                "original harmonic tangent residual is not representable",
            ));
        }
        let relative = if denominator == 0.0 {
            if sum.abs() == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            sum.abs() / denominator
        };
        let floor = gamma * rounding_scale;
        let reported = if relative <= tolerance {
            relative
        } else if sum.abs() <= floor {
            sum.abs() * tolerance / floor
        } else {
            relative
        };
        let block = if i < modes {
            0
        } else if i < modes + panels {
            1
        } else {
            2
        };
        residuals[block] = residuals[block].max(reported);
    }
    if residuals.iter().any(|&r| !r.is_finite() || r > tolerance) {
        return Err(BemModesError::PhysicalResidual {
            mechanical: residuals[0],
            boundary: residuals[1],
            memory: residuals[2],
        });
    }
    poll(cancelled)?;
    Ok(residuals)
}

#[cfg(test)]
#[path = "../harmonic_frequency_tests.rs"]
mod tests;
