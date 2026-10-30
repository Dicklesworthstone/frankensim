//! Bounded isothermal/isobaric vapor--liquid flash for Peng--Robinson.
//!
//! Wilson K-values only initialize and bracket the solve. Accepted endpoints
//! close Rachford--Rice material balance and component fugacity equality.

use fs_math::det;
use fs_qty::{Pressure, Temperature};

use super::{
    PengRobinsonErrorV1, PengRobinsonMixtureV1, PengRobinsonPhaseRootV1,
    PengRobinsonPhaseV1,
};

pub const MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1: usize = 256;
pub const MAX_RACHFORD_RICE_ITERATIONS_V1: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PengRobinsonFlashOptionsV1 {
    pub maximum_iterations: usize,
    pub rachford_rice_iterations: usize,
    pub log_k_tolerance: f64,
}
impl Default for PengRobinsonFlashOptionsV1 {
    fn default() -> Self {
        Self { maximum_iterations: 128, rachford_rice_iterations: 80, log_k_tolerance: 1.0e-10 }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum PengRobinsonFlashErrorV1 {
    Eos(PengRobinsonErrorV1),
    InvalidOptions,
    NoTwoPhaseBracket { f_liquid: f64, f_vapor: f64 },
    RachfordRiceUnresolved,
    FugacityIterationStalled { iterations: usize, residual: f64 },
    Unrepresentable(&'static str),
}
impl core::fmt::Display for PengRobinsonFlashErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Peng-Robinson flash refusal: {self:?}")
    }
}
impl std::error::Error for PengRobinsonFlashErrorV1 {}
impl From<PengRobinsonErrorV1> for PengRobinsonFlashErrorV1 {
    fn from(value: PengRobinsonErrorV1) -> Self { Self::Eos(value) }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PengRobinsonFlashV1 {
    temperature: Temperature,
    pressure: Pressure,
    vapor_fraction_mol: f64,
    liquid_mole_fractions: Vec<f64>,
    vapor_mole_fractions: Vec<f64>,
    liquid: PengRobinsonPhaseRootV1,
    vapor: PengRobinsonPhaseRootV1,
    iterations: usize,
    log_k_residual: f64,
    material_balance_residual: f64,
}
impl PengRobinsonFlashV1 {
    #[must_use] pub const fn temperature(&self) -> Temperature { self.temperature }
    #[must_use] pub const fn pressure(&self) -> Pressure { self.pressure }
    #[must_use] pub const fn vapor_fraction_mol(&self) -> f64 { self.vapor_fraction_mol }
    #[must_use] pub fn liquid_mole_fractions(&self) -> &[f64] { &self.liquid_mole_fractions }
    #[must_use] pub fn vapor_mole_fractions(&self) -> &[f64] { &self.vapor_mole_fractions }
    #[must_use] pub const fn liquid(&self) -> &PengRobinsonPhaseRootV1 { &self.liquid }
    #[must_use] pub const fn vapor(&self) -> &PengRobinsonPhaseRootV1 { &self.vapor }
    #[must_use] pub const fn iterations(&self) -> usize { self.iterations }
    #[must_use] pub const fn log_k_residual(&self) -> f64 { self.log_k_residual }
    #[must_use] pub const fn material_balance_residual(&self) -> f64 { self.material_balance_residual }
}

impl PengRobinsonMixtureV1 {
    /// Attempt a two-phase PT flash. Missing Wilson/Rachford--Rice bracket
    /// refuses; it is not silently reinterpreted as a stability proof.
    pub fn two_phase_flash_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
    ) -> Result<PengRobinsonFlashV1, PengRobinsonFlashErrorV1> {
        if options.maximum_iterations == 0
            || options.maximum_iterations > MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1
            || options.rachford_rice_iterations == 0
            || options.rachford_rice_iterations > MAX_RACHFORD_RICE_ITERATIONS_V1
            || !options.log_k_tolerance.is_finite()
            || options.log_k_tolerance <= 0.0
        {
            return Err(PengRobinsonFlashErrorV1::InvalidOptions);
        }
        self.evaluate_pt(temperature, pressure)?;
        let t = temperature.value();
        let p = pressure.value();
        let z = self.mole_fractions();
        let mut log_k = Vec::with_capacity(z.len());
        for component in self.components() {
            let exponent = 5.373 * (1.0 + component.acentric_factor)
                * (1.0 - component.critical_temperature.value() / t);
            let k = component.critical_pressure.value() / p * det::exp(exponent);
            if !k.is_finite() || k <= 0.0 {
                return Err(PengRobinsonFlashErrorV1::Unrepresentable("Wilson K-value"));
            }
            log_k.push(det::ln(k));
        }

        let components = self.components().to_vec();
        let kij = self.binary_interactions().to_vec();
        let mut last_residual = f64::INFINITY;
        for iteration in 0..options.maximum_iterations {
            let mut k = Vec::with_capacity(log_k.len());
            for value in &log_k {
                let next = det::exp(*value);
                if !next.is_finite() || next <= 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable("iterated K-value"));
                }
                k.push(next);
            }
            let beta = rachford_rice(z, &k, options.rachford_rice_iterations)?;
            let mut x = Vec::with_capacity(z.len());
            let mut y = Vec::with_capacity(z.len());
            for i in 0..z.len() {
                let denominator = 1.0 + beta * (k[i] - 1.0);
                if !denominator.is_finite() || denominator <= 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable("phase denominator"));
                }
                let xi = z[i] / denominator;
                let yi = k[i] * xi;
                if !xi.is_finite() || xi < 0.0 || !yi.is_finite() || yi < 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable("phase composition"));
                }
                x.push(xi);
                y.push(yi);
            }
            normalize(&mut x)?;
            normalize(&mut y)?;

            let liquid_model = PengRobinsonMixtureV1::new(components.clone(), x.clone(), kij.clone())?;
            let vapor_model = PengRobinsonMixtureV1::new(components.clone(), y.clone(), kij.clone())?;
            let liquid_eval = liquid_model.evaluate_pt(temperature, pressure)?;
            let vapor_eval = vapor_model.evaluate_pt(temperature, pressure)?;
            let liquid = liquid_eval.select(PengRobinsonPhaseV1::Liquid);
            let vapor = vapor_eval.select(PengRobinsonPhaseV1::Vapor);
            let lp = liquid.log_fugacity_coefficients();
            let vp = vapor.log_fugacity_coefficients();

            let mut target = Vec::with_capacity(z.len());
            let mut residual = 0.0_f64;
            for i in 0..z.len() {
                let value = lp[i] - vp[i];
                if !value.is_finite() {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable("fugacity K-value"));
                }
                residual = residual.max((value - log_k[i]).abs());
                target.push(value);
            }
            last_residual = residual;
            if residual <= options.log_k_tolerance {
                let mut material = 0.0_f64;
                for i in 0..z.len() {
                    material = material.max(
                        ((1.0 - beta) * x[i] + beta * y[i] - z[i]).abs()
                    );
                }
                if !material.is_finite() || material > 2048.0 * f64::EPSILON * z.len() as f64 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable("material balance"));
                }
                return Ok(PengRobinsonFlashV1 {
                    temperature, pressure, vapor_fraction_mol: beta,
                    liquid_mole_fractions: x, vapor_mole_fractions: y,
                    liquid: liquid.clone(), vapor: vapor.clone(),
                    iterations: iteration + 1, log_k_residual: residual,
                    material_balance_residual: material,
                });
            }
            for i in 0..log_k.len() {
                log_k[i] = 0.5 * (log_k[i] + target[i]);
            }
        }
        Err(PengRobinsonFlashErrorV1::FugacityIterationStalled {
            iterations: options.maximum_iterations,
            residual: last_residual,
        })
    }
}

fn normalize(values: &mut [f64]) -> Result<(), PengRobinsonFlashErrorV1> {
    let sum: f64 = values.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(PengRobinsonFlashErrorV1::Unrepresentable("phase composition sum"));
    }
    for value in values { *value /= sum; }
    Ok(())
}

fn rr_value(z: &[f64], k: &[f64], beta: f64) -> Result<f64, PengRobinsonFlashErrorV1> {
    let mut value = 0.0;
    for i in 0..z.len() {
        let denominator = 1.0 + beta * (k[i] - 1.0);
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(PengRobinsonFlashErrorV1::Unrepresentable("Rachford-Rice denominator"));
        }
        value += z[i] * (k[i] - 1.0) / denominator;
    }
    if !value.is_finite() {
        return Err(PengRobinsonFlashErrorV1::Unrepresentable("Rachford-Rice residual"));
    }
    Ok(value)
}

fn rachford_rice(
    z: &[f64], k: &[f64], iterations: usize,
) -> Result<f64, PengRobinsonFlashErrorV1> {
    let f0 = rr_value(z, k, 0.0)?;
    let f1 = rr_value(z, k, 1.0)?;
    if f0 <= 0.0 || f1 >= 0.0 {
        return Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket {
            f_liquid: f0, f_vapor: f1,
        });
    }
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..iterations {
        let mid = 0.5 * (lo + hi);
        if rr_value(z, k, mid)? > 0.0 { lo = mid; } else { hi = mid; }
    }
    let beta = 0.5 * (lo + hi);
    let residual = rr_value(z, k, beta)?;
    if !beta.is_finite() || beta <= 0.0 || beta >= 1.0
        || residual.abs() > 64.0 * f64::EPSILON * z.len() as f64
    {
        return Err(PengRobinsonFlashErrorV1::RachfordRiceUnresolved);
    }
    Ok(beta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::MolarMass;
    use crate::eos::peng_robinson::PengRobinsonComponentV1;

    fn methane() -> PengRobinsonComponentV1 {
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.016_042_46),
        }
    }
    fn ethane() -> PengRobinsonComponentV1 {
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(305.322),
            critical_pressure: Pressure::new(4.8722e6),
            acentric_factor: 0.0995,
            molar_mass: MolarMass::new(0.030_069_04),
        }
    }

    #[test]
    fn methane_ethane_flash_closes_material_and_fugacity() {
        let model = PengRobinsonMixtureV1::new(
            vec![methane(), ethane()], vec![0.5, 0.5], vec![],
        ).unwrap();
        let flash = model.two_phase_flash_pt(
            Temperature::new(220.0), Pressure::new(3.0e6), Default::default(),
        ).unwrap();
        assert!((flash.vapor_fraction_mol() - 0.302_585_915_664_206_8).abs() < 2e-8);
        for (actual, expected) in flash.liquid_mole_fractions().iter().zip([
            0.364_782_477_310_312_44, 0.635_217_522_689_687_6,
        ]) {
            assert!((*actual - expected).abs() < 2e-8);
        }
        for (actual, expected) in flash.vapor_mole_fractions().iter().zip([
            0.811_655_631_974_075_9, 0.188_344_368_025_924_16,
        ]) {
            assert!((*actual - expected).abs() < 2e-8);
        }
        assert!(flash.log_k_residual() <= 1.0e-10);
        assert!(flash.material_balance_residual() < 1.0e-12);
        let lf = flash.liquid().log_fugacity_coefficients();
        let vf = flash.vapor().log_fugacity_coefficients();
        for i in 0..2 {
            let lhs = det::ln(flash.liquid_mole_fractions()[i]) + lf[i];
            let rhs = det::ln(flash.vapor_mole_fractions()[i]) + vf[i];
            assert!((lhs - rhs).abs() < 2e-10);
        }
    }

    #[test]
    fn flash_refuses_absent_bracket_and_invalid_work() {
        let model = PengRobinsonMixtureV1::new(
            vec![methane(), ethane()], vec![0.5, 0.5], vec![],
        ).unwrap();
        assert!(matches!(
            model.two_phase_flash_pt(
                Temperature::new(400.0), Pressure::new(1.0e5), Default::default(),
            ),
            Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket { .. })
        ));
        assert!(matches!(
            model.two_phase_flash_pt(
                Temperature::new(220.0), Pressure::new(3.0e6),
                PengRobinsonFlashOptionsV1 { maximum_iterations: 0, ..Default::default() },
            ),
            Err(PengRobinsonFlashErrorV1::InvalidOptions)
        ));
    }
}
