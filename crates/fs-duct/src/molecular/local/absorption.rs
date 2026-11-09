//! Positive losses of the SAME original joint molecular/wall pressure law.
//!
//! In exp(-i omega t), Cp has nonnegative imaginary part. The original
//! radial temperature mode has D_n=Cp+i*d_n/omega, d_n=Cp_eq*lambda_n and
//! mean-square weight w_n. Its Fourier and reservoir availability losses per
//! length are A/(2*rho*T0) w_n [d_n, omega*Im(Cp)] |p/D_n|^2.
//! They add to the original shunt resistance without subtracting two powers.
//! The original positive series resistance supplies direct viscous heat.
//!
//! These are cycle-mean losses for peak p/U phasors. Fourier/reservoir loss
//! does not independently specify caloric deposition or wall/fluid transport.

use super::*;
use fs_math::dd::Dd;

/// One original circular constitutive sample, retained at field preparation.
/// The material, wall measure, Cp and wave are from the same native query.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MolecularCylinder<'source> {
    pub(crate) state: &'source RelaxationState,
    pub(crate) wall: Option<&'source ZkDiffusive>,
    pub(crate) radius_m: f64,
    pub(crate) x_m: f64,
    pub(crate) omega: f64,
    pub(crate) cp: C64,
    pub(crate) wave: SegmentWave,
    // Same queried mechanical surface law used in the retained composed wave.
    // Positive real part is wall work, not an additional gas heat mechanism.
    pub(crate) mechanical_admittance_per_metre: C64,
}

/// Geometry and the same frozen radial sample of one native spherical slice.
/// The acoustic specific coefficients are constant on this slice; their
/// volume-flow normalization follows its actual varying cross-section.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MolecularCone<'source> {
    pub(crate) inlet_radius_m: f64,
    pub(crate) outlet_radius_m: f64,
    pub(crate) length_m: f64,
    pub(crate) original: MolecularCylinder<'source>,
}

impl MolecularCylinder<'_> {
    /// Diagonal per-length weights for [|U|^2, |p|^2, |p|^2], in order
    /// viscous, Fourier and molecular. The shared interval owner applies 1/2.
    /// No acoustic or caloric material is queried again.
    pub(crate) fn loss_weights(&self) -> Result<[f64; 3], DuctError> {
        let gas = self.state.gas();
        let area = core::f64::consts::PI * self.radius_m * self.radius_m;
        let mut result = [0.0; 3];
        if let Some(wall) = self.wall {
            let (series, equilibrium_shunt) = wall.response_at(self.omega).map_err(|reason| {
                DuctError::WallDiffusionUnresolved {
                    reason,
                    radius_m: self.radius_m,
                    x_m: self.x_m,
                }
            })?;
            result[0] = series.re * (gas.density / area);
            if self.state.processes().is_empty() {
                // Preserve the original equilibrium branch and its exact
                // sound-speed normalization; no molecular channel is added.
                result[1] = equilibrium_shunt.re
                    * (area / (gas.density * gas.sound_speed * gas.sound_speed));
            } else {
                let scale = (area / gas.density) / gas.temperature;
                let mut losses = [Dd::ZERO; 2];
                for term in wall.shunt_terms() {
                    let (weight, denominator, wall_capacity) =
                        thermal_mode(self.state, term, self.omega, self.cp);
                    let root = denominator
                        .recip()
                        .scale(det::sqrt(weight))
                        .scale(det::sqrt(scale))
                        .scale(det::sqrt(self.omega));
                    for (sum, capacity) in losses.iter_mut().zip([wall_capacity, self.cp.im]) {
                        *sum = *sum + positive_square(root, capacity)?;
                    }
                }
                result[1] = losses[0].to_f64();
                result[2] = losses[1].to_f64();
            }
        } else if self.cp.im != 0.0 {
            // The native bare modulus is p*Cp/Cv. Retain its ACTUAL pressure
            // normalization rather than replacing p by approximately rho*R*T.
            let scale = (area / gas.pressure) * gas.specific_gas_constant;
            let root = self
                .cp
                .recip()
                .scale(det::sqrt(scale))
                .scale(det::sqrt(self.omega));
            result[2] = positive_square(root, self.cp.im)?.to_f64();
        }
        if result.into_iter().any(|v| !v.is_finite() || v < 0.0)
            || (self.cp.im > 0.0 && result[2] == 0.0)
            || (self.wall.is_some() && result[0] == 0.0)
        {
            return Err(bad("original molecular cylinder loss is unrepresentable"));
        }
        Ok(result)
    }
}

fn positive_square(root: C64, capacity: f64) -> Result<Dd, DuctError> {
    if capacity == 0.0 {
        return Ok(Dd::ZERO);
    }
    if !capacity.is_finite() || capacity < 0.0 || !finite(root) {
        return Err(bad(
            "original molecular temperature loss form is unrepresentable",
        ));
    }
    let z = root.scale(det::sqrt(capacity));
    let re = Dd::from_f64(z.re);
    let im = Dd::from_f64(z.im);
    let power = re * re + im * im;
    if !finite(z) || !power.to_f64().is_finite() || power.to_f64() <= 0.0 {
        return Err(bad(
            "positive molecular temperature loss is unrepresentable",
        ));
    }
    Ok(power)
}
