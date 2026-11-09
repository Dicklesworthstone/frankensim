//! Conservative observations of the original exponential-density Mann field.
//!
//! With rho=rho_c exp(-kappa*dz) and u=(rho_ref/rho)v, Favre filtering gives
//! mean(rho*u)=rho_ref*mean(v), u_tilde=mean(rho*u)/mean(rho), and the positive
//! central stress tau=mean(rho*u*u)/mean(rho)-u_tilde*u_tilde. No density at a
//! representative point may replace either of those integrals.
//!
//! Write f=exp(kappa*dz/2)v, h=exp(-kappa*dz/2), D=mean(h*h). The stress is
//! (rho_ref/rho_c)^2/D times the Gram matrix of f after its projection on h.
//! For centered covariances C, that matrix is
//! C_fg + (C_hh*m_f*m_g-m_h*(m_f*C_gh+m_g*C_fh)-C_fh*C_gh)/D.
//! Every term is small on a small footprint; no nearly equal raw moments are
//! subtracted. The existing bounded complex-exponential covariance kernel
//! supplies these moments. All features share exp(-kappa*vertical_radius/2)
//! normalization, which cancels analytically without a growing exponential.
//!
//! These are spatial observations of the same prescribed outer-flow source,
//! not a wall model, buoyancy evolution or a new SGS constitutive closure.

use super::MannBoxFilter;
use super::moments::MannBoxMoments;
use crate::fetch_mass_consistent::ExponentialDensityProfile;
use crate::mann_field::wall::filter::{moments::pair_covariance, scaled_edge_gain};
use crate::mann_field::{MannSpectralModel, dot, finite_vector, sample_modes};
use crate::{FlowSample, MAX_MODES, Refusal, refuse};
use fs_math::det;

type Complex = [f64; 2];
type Exponents = [Complex; 3];
type Covariance = [[f64; 2]; 2];

/// Exact represented mass, momentum and kinetic energy in a fixed footprint.
#[derive(Clone, Debug, PartialEq)]
pub struct MannFavreMoments {
    /// Uniform spatial average of physical density [kg/m³].
    pub mean_density_kg_m3: f64,
    /// Uniform average of rho*u [kg/(m² s)], independent of density contrast
    /// at fixed reference density and original solenoidal realization.
    pub mean_mass_flux_kg_m2_s: [f64; 3],
    /// Favre velocity [m/s] and its derivative as the fixed footprint translates
    /// [1/s]. This is a density-weighted mean, not a uniform velocity average.
    pub favre_mean: FlowSample,
    /// mean(rho*(u-u_tilde)*(u-u_tilde))/mean(rho) [m²/s²].
    /// Positive semidefinite in exact arithmetic; roundoff is not clipped.
    pub favre_subfilter_stress_m2_s2: [[f64; 3]; 3],
    /// Uniform average of rho*u_i*u_j [Pa]. No pressure or viscous stress.
    pub convective_momentum_flux_pa: [[f64; 3]; 3],
    /// Uniform average of rho*|u|²/2 [J/m³]; no thermal energy or heating.
    pub kinetic_energy_density_j_m3: f64,
}

#[derive(Clone, Debug, PartialEq)]
enum Observation {
    Constant(MannBoxMoments),
    Stratified {
        filter: MannBoxFilter,
        density_factor: f64,
        coefficients: Vec<Covariance>,
    },
}

/// Source-bound exact Favre moments of the existing anelastic perturbation.
/// Boxes, skew parallelepipeds and collapsed faces/lines are admitted. At each
/// query the entire footprint must lie inside the retained density interval.
/// The original Mann frame's +z must denote the profile's physical height.
#[derive(Clone, Debug, PartialEq)]
pub struct MannAnelasticBoxMoments {
    density: ExponentialDensityProfile,
    vertical_radius_m: f64,
    observation: Observation,
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use the matching Mann source, a finite footprint inside the density domain, and sufficient explicit storage and mode-pair budgets",
    )
}

fn numerical() -> Refusal {
    error(
        "mann-anelastic-moments-range",
        "density-weighted moment, gradient or physical flux is not representable",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-anelastic-moments-cancelled",
            "density-weighted observation cancelled",
        ))
    } else {
        Ok(())
    }
}

fn feature_mean(exponents: Exponents) -> Complex {
    let mut mean = [1.0, 0.0];
    for [real, imaginary] in exponents {
        let gain = scaled_edge_gain(real, imaginary);
        mean = [
            mean[0] * gain[0] - mean[1] * gain[1],
            mean[0] * gain[1] + mean[1] * gain[0],
        ];
    }
    mean
}

impl MannAnelasticBoxMoments {
    /// Prepare exact products with no sampling approximation. Constant density
    /// or zero vertical extent reuses the original MannBoxMoments. Otherwise
    /// each stored unordered mode pair uses four f64 coefficients, plus the
    /// original linear-filter gains. Preparation additionally evaluates n+1
    /// density-feature pairs and uses 20 KiB of fixed feature scratch, plus
    /// the existing kernel's bounded local temporaries. Sampling costs O(n²),
    /// uses 12 KiB of modal-vector scratch, and allocates no storage.
    /// The pair limit counts stored mode pairs; the byte limit counts actual
    /// retained numerical capacities, excluding fixed metadata and stack.
    pub fn prepare(
        model: &MannSpectralModel,
        density: ExponentialDensityProfile,
        half_edges_m: [[f64; 3]; 3],
        maximum_workspace_bytes: usize,
        maximum_mode_pairs: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(&mut cancelled)?;
        let radius = half_edges_m.iter().map(|e| e[2].abs()).sum::<f64>();
        if !radius.is_finite() || !half_edges_m.iter().flatten().all(|v| v.is_finite()) {
            return Err(numerical());
        }
        let kappa = density.inverse_scale_height_per_m();
        let observation = if kappa == 0.0 || radius == 0.0 {
            Observation::Constant(MannBoxMoments::prepare(
                model,
                half_edges_m,
                maximum_workspace_bytes,
                maximum_mode_pairs,
                &mut cancelled,
            )?)
        } else {
            let filter = MannBoxFilter::prepare(
                model,
                half_edges_m,
                maximum_workspace_bytes,
                &mut cancelled,
            )?;
            let n = model.mode_count();
            let pairs = n * (n + 1) / 2;
            let bytes = filter.workspace_bytes() + pairs * core::mem::size_of::<Covariance>();
            if pairs > maximum_mode_pairs || bytes > maximum_workspace_bytes {
                return Err(error(
                    "mann-anelastic-moments-budget",
                    "density-weighted observation exceeds its pair or byte budget",
                ));
            }
            let mut coefficients = Vec::new();
            coefficients.try_reserve_exact(pairs).map_err(|_| {
                error(
                    "mann-anelastic-moments-allocation",
                    "unable to reserve moments",
                )
            })?;
            if filter.workspace_bytes()
                + coefficients.capacity() * core::mem::size_of::<Covariance>()
                > maximum_workspace_bytes
            {
                return Err(error(
                    "mann-anelastic-moments-budget",
                    "allocated moment capacities exceed the byte budget",
                ));
            }
            let h = half_edges_m.map(|edge| [-0.5 * kappa * edge[2], 0.0]);
            if !h.iter().flatten().all(|v| v.is_finite()) {
                return Err(numerical());
            }
            let mean_h = feature_mean(h)[0];
            let variance_h = pair_covariance(h, h)?[0][0];
            // D includes the common bounded-feature normalization.
            let density_factor = mean_h * mean_h + variance_h;
            if !density_factor.is_finite() || density_factor <= 0.0 {
                return Err(numerical());
            }
            let mut exponents = [[[0.0; 2]; 3]; MAX_MODES];
            let mut means = [[0.0; 2]; MAX_MODES];
            let mut cross_h = [[0.0; 2]; MAX_MODES];
            for (m, mode) in model.modes.iter().enumerate() {
                poll(&mut cancelled)?;
                exponents[m] = core::array::from_fn(|j| [-h[j][0], dot(mode.k, half_edges_m[j])]);
                means[m] = feature_mean(exponents[m]);
                let cross = pair_covariance(exponents[m], h)?;
                cross_h[m] = [cross[0][0], cross[1][0]];
            }
            for m in 0..n {
                for j in m..n {
                    poll(&mut cancelled)?;
                    let central = pair_covariance(exponents[m], exponents[j])?;
                    let coefficient = core::array::from_fn(|a| {
                        core::array::from_fn(|b| {
                            let correction = (variance_h * means[m][a] * means[j][b]
                                - mean_h
                                    * (means[m][a] * cross_h[j][b] + means[j][b] * cross_h[m][a])
                                - cross_h[m][a] * cross_h[j][b])
                                / density_factor;
                            (central[a][b] + correction) / density_factor
                        })
                    });
                    if !coefficient.iter().flatten().all(|v| v.is_finite()) {
                        return Err(numerical());
                    }
                    coefficients.push(coefficient);
                }
            }
            Observation::Stratified {
                filter,
                density_factor,
                coefficients,
            }
        };
        poll(&mut cancelled)?;
        Ok(Self {
            density,
            vertical_radius_m: radius,
            observation,
        })
    }

    pub fn workspace_bytes(&self) -> usize {
        match &self.observation {
            Observation::Constant(original) => original.workspace_bytes(),
            Observation::Stratified {
                filter,
                coefficients,
                ..
            } => {
                filter.workspace_bytes()
                    + coefficients.capacity() * core::mem::size_of::<Covariance>()
            }
        }
    }

    /// Stored unordered mode pairs evaluated by each non-point observation.
    pub fn mode_pairs(&self) -> usize {
        match &self.observation {
            Observation::Constant(original) => original.mode_pairs(),
            Observation::Stratified { coefficients, .. } => coefficients.len(),
        }
    }
}

fn finish(
    density: f64,
    mass_flux: [f64; 3],
    mean: FlowSample,
    stress: [[f64; 3]; 3],
) -> Result<MannFavreMoments, Refusal> {
    let mut momentum = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            momentum[i][j] = mass_flux[i] * mean.u[j] + density * stress[i][j];
            momentum[j][i] = momentum[i][j];
        }
    }
    let energy = (0..3).map(|i| 0.5 * momentum[i][i]).sum::<f64>();
    if !density.is_finite()
        || density <= 0.0
        || !energy.is_finite()
        || mass_flux
            .iter()
            .chain(mean.u.iter())
            .chain(mean.grad.iter().flatten())
            .chain(stress.iter().flatten())
            .chain(momentum.iter().flatten())
            .any(|v| !v.is_finite())
    {
        return Err(numerical());
    }
    Ok(MannFavreMoments {
        mean_density_kg_m3: density,
        mean_mass_flux_kg_m2_s: mass_flux,
        favre_mean: mean,
        favre_subfilter_stress_m2_s2: stress,
        convective_momentum_flux_pa: momentum,
        kinetic_energy_density_j_m3: energy,
    })
}

pub(in crate::mann_field) fn sample(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    point: [f64; 3],
    observation: &MannAnelasticBoxMoments,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MannFavreMoments, Refusal> {
    poll(&mut cancelled)?;
    finite_vector(point)?;
    let density = observation.density;
    let rho_low = density.density_kg_m3(point[2] - observation.vertical_radius_m)?;
    density.density_kg_m3(point[2] + observation.vertical_radius_m)?;
    let rho_center = density.density_kg_m3(point[2])?;
    let kappa = density.inverse_scale_height_per_m();
    let result = match &observation.observation {
        Observation::Constant(original) => {
            let old =
                super::moments::sample(model, time, amplitudes, point, original, &mut cancelled)?;
            let flow = density.compose_anelastic(None, old.mean, point[2])?;
            let mut stress = old.subfilter_stress_m2_s2;
            if kappa != 0.0 {
                let scale = density.reference_density_kg_m3() / rho_center;
                for value in stress.iter_mut().flatten() {
                    *value = scale * (*value * scale);
                }
            }
            finish(rho_center, flow.mass_flux_kg_m2_s, flow.flow, stress)?
        }
        Observation::Stratified {
            filter,
            density_factor,
            coefficients,
        } => {
            let original = sample_modes(model, time, amplitudes, point, filter.gains(model)?)?;
            // rho_low contains the opposite exponential normalization. This
            // stays within endpoint densities without ever forming exp(kappa*r).
            let rho_bar = rho_low * density_factor;
            let mean_scale = density.reference_density_kg_m3() / rho_bar;
            let local_scale = density.reference_density_kg_m3() / rho_center;
            if ![rho_bar, mean_scale, local_scale]
                .into_iter()
                .all(|v| v.is_finite() && v > 0.0)
            {
                return Err(numerical());
            }
            let mass_flux = original.u.map(|v| density.reference_density_kg_m3() * v);
            let mut mean = original;
            for (velocity, gradient) in mean.u.iter_mut().zip(&mut mean.grad) {
                *velocity *= mean_scale;
                for derivative in &mut *gradient {
                    *derivative *= mean_scale;
                }
                gradient[2] += kappa * *velocity;
            }
            let position = core::array::from_fn(|i| point[i] - model.advection_mps[i] * time);
            let mut vectors = [[[0.0; 3]; 2]; MAX_MODES];
            for (m, mode) in model.modes.iter().enumerate() {
                poll(&mut cancelled)?;
                let (s, c) = det::sin_cos(dot(mode.k, position));
                for i in 0..3 {
                    vectors[m][0][i] = local_scale
                        * (mode.polarization[0][i] * (amplitudes[0][m] * c + amplitudes[1][m] * s)
                            + mode.polarization[1][i]
                                * (amplitudes[2][m] * c + amplitudes[3][m] * s));
                    vectors[m][1][i] = local_scale
                        * (mode.polarization[0][i]
                            * (-amplitudes[0][m] * s + amplitudes[1][m] * c)
                            + mode.polarization[1][i]
                                * (-amplitudes[2][m] * s + amplitudes[3][m] * c));
                }
            }
            let mut stress = [[0.0; 3]; 3];
            let mut pairs = coefficients.iter();
            for m in 0..model.mode_count() {
                for n in m..model.mode_count() {
                    poll(&mut cancelled)?;
                    let covariance = pairs.next().expect("prepared mode-pair count");
                    for i in 0..3 {
                        for j in i..3 {
                            for a in 0..2 {
                                for b in 0..2 {
                                    let mut product = vectors[m][a][i] * vectors[n][b][j];
                                    if m != n {
                                        product += vectors[n][b][i] * vectors[m][a][j];
                                    }
                                    stress[i][j] += covariance[a][b] * product;
                                }
                            }
                        }
                    }
                }
            }
            for i in 0..3 {
                for j in 0..i {
                    stress[i][j] = stress[j][i];
                }
            }
            finish(rho_bar, mass_flux, mean, stress)?
        }
    };
    poll(&mut cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
