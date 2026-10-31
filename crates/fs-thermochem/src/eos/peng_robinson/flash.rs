//! Bounded isothermal/isobaric vapor--liquid flash for Peng--Robinson.
//!
//! Wilson K-values initialize the solve. If that estimate misses a phase
//! split, bounded tangent-plane searches supply incipient-phase compositions.
//! Single-phase iterates are allowed, but published endpoints must have a
//! nontrivial split, material/fugacity closure, and no Gibbs-energy increase.

use fs_math::det;
use fs_qty::{MolarMass, Pressure, Temperature};

use super::{
    PengRobinsonErrorV1, PengRobinsonMixtureV1, PengRobinsonPhaseRootV1, PengRobinsonPhaseV1,
};

pub const MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1: usize = 256;
pub const MAX_RACHFORD_RICE_ITERATIONS_V1: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PengRobinsonFlashOptionsV1 {
    /// Per-search limit. There are at most two split attempts and two
    /// incipient-phase searches, or at most `6 * maximum_iterations + 2`
    /// evaluations of the bounded cubic EOS in the complete call.
    pub maximum_iterations: usize,
    pub rachford_rice_iterations: usize,
    pub log_k_tolerance: f64,
}
impl Default for PengRobinsonFlashOptionsV1 {
    fn default() -> Self {
        Self {
            maximum_iterations: 128,
            rachford_rice_iterations: 80,
            log_k_tolerance: 1.0e-10,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum PengRobinsonFlashErrorV1 {
    Eos(PengRobinsonErrorV1),
    InvalidOptions,
    NoTwoPhaseBracket { f_liquid: f64, f_vapor: f64 },
    RachfordRiceUnresolved,
    FugacityIterationStalled { iterations: usize, residual: f64 },
    StabilityIterationStalled { iterations: usize, residual: f64 },
    NoSeparatedPhases,
    GibbsEnergyIncreased { dimensionless_change: f64 },
    Cancelled,
    Unrepresentable(&'static str),
}
impl core::fmt::Display for PengRobinsonFlashErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Peng-Robinson flash refusal: {self:?}")
    }
}
impl std::error::Error for PengRobinsonFlashErrorV1 {}
impl From<PengRobinsonErrorV1> for PengRobinsonFlashErrorV1 {
    fn from(value: PengRobinsonErrorV1) -> Self {
        Self::Eos(value)
    }
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
    gibbs_change_over_rt: f64,
}
impl PengRobinsonFlashV1 {
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    #[must_use]
    pub const fn vapor_fraction_mol(&self) -> f64 {
        self.vapor_fraction_mol
    }
    #[must_use]
    pub fn liquid_mole_fractions(&self) -> &[f64] {
        &self.liquid_mole_fractions
    }
    #[must_use]
    pub fn vapor_mole_fractions(&self) -> &[f64] {
        &self.vapor_mole_fractions
    }
    #[must_use]
    pub const fn liquid(&self) -> &PengRobinsonPhaseRootV1 {
        &self.liquid
    }
    #[must_use]
    pub const fn vapor(&self) -> &PengRobinsonPhaseRootV1 {
        &self.vapor
    }
    #[must_use]
    /// Iterations in the accepted split attempt. This excludes a failed
    /// initial attempt and the optional two incipient-phase searches; the
    /// complete work ceiling is documented by `maximum_iterations`.
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    #[must_use]
    pub const fn log_k_residual(&self) -> f64 {
        self.log_k_residual
    }
    #[must_use]
    pub const fn material_balance_residual(&self) -> f64 {
        self.material_balance_residual
    }
    /// Molar Gibbs change relative to the minimum-Gibbs feed root, divided by
    /// RT. Standard chemical references cancel by component material balance.
    #[must_use]
    pub const fn gibbs_change_over_rt(&self) -> f64 {
        self.gibbs_change_over_rt
    }
}

impl PengRobinsonMixtureV1 {
    /// Attempt a two-phase PT flash, with tangent-plane initialization when
    /// Wilson's approximation misses the real-EOS split. Absence of a split
    /// refuses; it is not silently reinterpreted as a global stability proof.
    pub fn two_phase_flash_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
    ) -> Result<PengRobinsonFlashV1, PengRobinsonFlashErrorV1> {
        self.two_phase_flash_pt_cancellable(temperature, pressure, options, || false)
    }

    /// The same equilibrium closure, polling before each phase evaluation.
    pub fn two_phase_flash_pt_cancellable(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
        mut cancelled: impl FnMut() -> bool,
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
        if cancelled() {
            return Err(PengRobinsonFlashErrorV1::Cancelled);
        }
        let feed = self.evaluate_pt(temperature, pressure)?;
        let feed_root = feed.select(PengRobinsonPhaseV1::Stable);
        let log_k = wilson_log_k(self, temperature, pressure)?;
        let first = self.flash_from_seed(
            temperature,
            pressure,
            options,
            log_k,
            feed_root,
            &mut cancelled,
        );
        if !matches!(
            first,
            Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket { .. })
                | Err(PengRobinsonFlashErrorV1::NoSeparatedPhases)
                | Err(PengRobinsonFlashErrorV1::FugacityIterationStalled { .. })
        ) {
            return first;
        }
        let stability = self.stability_pt_cancellable(
            temperature,
            pressure,
            super::stability::PengRobinsonStabilityOptionsV1 {
                maximum_iterations: options.maximum_iterations,
                stationarity_tolerance: options.log_k_tolerance,
                tangent_plane_tolerance: options.log_k_tolerance,
            },
            &mut cancelled,
        )?;
        let Some(witness) = stability.instability_witness() else {
            return first;
        };
        let vapor_like =
            witness.root().molar_volume_m3_per_mol() > feed_root.molar_volume_m3_per_mol();
        let sign = if vapor_like { 1.0 } else { -1.0 };
        let log_k = self
            .mole_fractions()
            .iter()
            .zip(witness.mole_fractions())
            .map(|(&z, &w)| {
                if z > 0.0 {
                    sign * (det::ln(w) - det::ln(z))
                } else {
                    0.0
                }
            })
            .collect();
        self.flash_from_seed(
            temperature,
            pressure,
            options,
            log_k,
            feed_root,
            &mut cancelled,
        )
    }

    fn flash_from_seed(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
        mut log_k: Vec<f64>,
        feed: &PengRobinsonPhaseRootV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PengRobinsonFlashV1, PengRobinsonFlashErrorV1> {
        let z = self.mole_fractions();
        let mut liquid_model = self.clone();
        let mut vapor_model = self.clone();
        let mut k = vec![0.0; z.len()];
        let mut x = vec![0.0; z.len()];
        let mut y = vec![0.0; z.len()];
        let mut target = vec![0.0; z.len()];
        let mut last_residual = f64::INFINITY;
        for iteration in 0..options.maximum_iterations {
            if cancelled() {
                return Err(PengRobinsonFlashErrorV1::Cancelled);
            }
            for (i, value) in log_k.iter().enumerate() {
                let next = if z[i] == 0.0 { 1.0 } else { det::exp(*value) };
                if !next.is_finite() || next <= 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "iterated K-value",
                    ));
                }
                k[i] = next;
            }
            // Endpoint iterates describe an incipient phase, not a published
            // two-phase equilibrium. This permits the real EOS to correct a
            // Wilson estimate outside its approximate bubble/dew bracket.
            let beta = match rachford_rice(z, &k, options.rachford_rice_iterations) {
                Ok(beta) => beta,
                Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket { f_liquid, .. }) => {
                    if f_liquid <= 0.0 { 0.0 } else { 1.0 }
                }
                Err(error) => return Err(error),
            };
            for i in 0..z.len() {
                let denominator = 1.0 + beta * (k[i] - 1.0);
                if !denominator.is_finite() || denominator <= 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "phase denominator",
                    ));
                }
                let xi = z[i] / denominator;
                let yi = k[i] * xi;
                if !xi.is_finite() || xi < 0.0 || !yi.is_finite() || yi < 0.0 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "phase composition",
                    ));
                }
                x[i] = xi;
                y[i] = yi;
            }
            normalize(&mut x)?;
            normalize(&mut y)?;

            set_composition(&mut liquid_model, &x)?;
            set_composition(&mut vapor_model, &y)?;
            let liquid_eval = liquid_model.evaluate_pt(temperature, pressure)?;
            if cancelled() {
                return Err(PengRobinsonFlashErrorV1::Cancelled);
            }
            let vapor_eval = vapor_model.evaluate_pt(temperature, pressure)?;
            let liquid = liquid_eval.select(PengRobinsonPhaseV1::Liquid);
            let vapor = vapor_eval.select(PengRobinsonPhaseV1::Vapor);
            let lp = liquid.log_fugacity_coefficients();
            let vp = vapor.log_fugacity_coefficients();

            let mut residual = 0.0_f64;
            for i in 0..z.len() {
                if z[i] == 0.0 {
                    target[i] = 0.0;
                    continue;
                }
                let value = lp[i] - vp[i];
                if !value.is_finite() {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "fugacity K-value",
                    ));
                }
                residual = residual.max((value - log_k[i]).abs());
                target[i] = value;
            }
            last_residual = residual;
            if residual <= options.log_k_tolerance {
                if beta <= 0.0 || beta >= 1.0 {
                    return Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket {
                        f_liquid: rr_value(z, &k, 0.0)?,
                        f_vapor: rr_value(z, &k, 1.0)?,
                    });
                }
                let separation = x
                    .iter()
                    .zip(&y)
                    .map(|(x, y)| (x - y).abs())
                    .fold(0.0_f64, f64::max);
                // At a pure-fluid/azeotropic composition, PT and feed
                // composition cannot determine phase amount. Distinct cubic
                // volumes alone must not legitimize a roundoff-selected beta.
                if separation <= 32.0 * options.log_k_tolerance {
                    return Err(PengRobinsonFlashErrorV1::NoSeparatedPhases);
                }
                let mut material = 0.0_f64;
                let mut gibbs_change = 0.0;
                for i in 0..z.len() {
                    material = material.max(((1.0 - beta) * x[i] + beta * y[i] - z[i]).abs());
                    if z[i] > 0.0 {
                        if x[i] <= 0.0 || y[i] <= 0.0 {
                            return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                                "phase component underflow",
                            ));
                        }
                        let reference = det::ln(z[i]) + feed.log_fugacity_coefficients()[i];
                        gibbs_change += (1.0 - beta) * x[i] * (det::ln(x[i]) + lp[i] - reference)
                            + beta * y[i] * (det::ln(y[i]) + vp[i] - reference);
                    }
                }
                if !material.is_finite() || material > 2048.0 * f64::EPSILON * z.len() as f64 {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "material balance",
                    ));
                }
                if !gibbs_change.is_finite() {
                    return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                        "flash Gibbs change",
                    ));
                }
                if gibbs_change > options.log_k_tolerance {
                    return Err(PengRobinsonFlashErrorV1::GibbsEnergyIncreased {
                        dimensionless_change: gibbs_change,
                    });
                }
                if cancelled() {
                    return Err(PengRobinsonFlashErrorV1::Cancelled);
                }
                return Ok(PengRobinsonFlashV1 {
                    temperature,
                    pressure,
                    vapor_fraction_mol: beta,
                    liquid_mole_fractions: x,
                    vapor_mole_fractions: y,
                    liquid: liquid.clone(),
                    vapor: vapor.clone(),
                    iterations: iteration + 1,
                    log_k_residual: residual,
                    material_balance_residual: material,
                    gibbs_change_over_rt: gibbs_change,
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

pub(super) fn wilson_log_k(
    model: &PengRobinsonMixtureV1,
    temperature: Temperature,
    pressure: Pressure,
) -> Result<Vec<f64>, PengRobinsonFlashErrorV1> {
    model
        .components()
        .iter()
        .map(|component| {
            let value = det::ln(component.critical_pressure.value()) - det::ln(pressure.value())
                + 5.373
                    * (1.0 + component.acentric_factor)
                    * (1.0 - component.critical_temperature.value() / temperature.value());
            if value.is_finite() {
                Ok(value)
            } else {
                Err(PengRobinsonFlashErrorV1::Unrepresentable(
                    "Wilson log K-value",
                ))
            }
        })
        .collect()
}

/// Reuse already-admitted immutable component and interaction data. Only the
/// trial composition and its physical molar mass change; no O(n²) parameter
/// cloning or readmission is needed inside the nonlinear iterations.
pub(super) fn set_composition(
    model: &mut PengRobinsonMixtureV1,
    fractions: &[f64],
) -> Result<(), PengRobinsonFlashErrorV1> {
    let mut mass = 0.0;
    for (i, &x) in fractions.iter().enumerate() {
        if !x.is_finite() || x < 0.0 {
            return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                "trial mole fraction",
            ));
        }
        mass += x * model.components[i].molar_mass.value();
    }
    if !mass.is_finite() || mass <= 0.0 {
        return Err(PengRobinsonFlashErrorV1::Unrepresentable(
            "trial molar mass",
        ));
    }
    model.mole_fractions.copy_from_slice(fractions);
    model.molar_mass = MolarMass::new(mass);
    Ok(())
}

fn normalize(values: &mut [f64]) -> Result<(), PengRobinsonFlashErrorV1> {
    let sum: f64 = values.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(PengRobinsonFlashErrorV1::Unrepresentable(
            "phase composition sum",
        ));
    }
    for value in values {
        *value /= sum;
    }
    Ok(())
}

fn rr_value(z: &[f64], k: &[f64], beta: f64) -> Result<f64, PengRobinsonFlashErrorV1> {
    let mut value = 0.0;
    for i in 0..z.len() {
        let denominator = 1.0 + beta * (k[i] - 1.0);
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(PengRobinsonFlashErrorV1::Unrepresentable(
                "Rachford-Rice denominator",
            ));
        }
        value += z[i] * (k[i] - 1.0) / denominator;
    }
    if !value.is_finite() {
        return Err(PengRobinsonFlashErrorV1::Unrepresentable(
            "Rachford-Rice residual",
        ));
    }
    Ok(value)
}

fn rachford_rice(z: &[f64], k: &[f64], iterations: usize) -> Result<f64, PengRobinsonFlashErrorV1> {
    let f0 = rr_value(z, k, 0.0)?;
    let f1 = rr_value(z, k, 1.0)?;
    if f0 <= 0.0 || f1 >= 0.0 {
        return Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket {
            f_liquid: f0,
            f_vapor: f1,
        });
    }
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..iterations {
        let mid = 0.5 * (lo + hi);
        if rr_value(z, k, mid)? > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let beta = 0.5 * (lo + hi);
    let residual = rr_value(z, k, beta)?;
    if !beta.is_finite()
        || beta <= 0.0
        || beta >= 1.0
        || residual.abs() > 64.0 * f64::EPSILON * z.len() as f64
    {
        return Err(PengRobinsonFlashErrorV1::RachfordRiceUnresolved);
    }
    Ok(beta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eos::peng_robinson::PengRobinsonComponentV1;
    use fs_qty::MolarMass;

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
        let model =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.5, 0.5], vec![]).unwrap();
        let flash = model
            .two_phase_flash_pt(
                Temperature::new(220.0),
                Pressure::new(3.0e6),
                Default::default(),
            )
            .unwrap();
        assert!((flash.vapor_fraction_mol() - 0.302_585_915_664_206_8).abs() < 2e-8);
        for (actual, expected) in flash
            .liquid_mole_fractions()
            .iter()
            .zip([0.364_782_477_310_312_44, 0.635_217_522_689_687_6])
        {
            assert!((*actual - expected).abs() < 2e-8);
        }
        for (actual, expected) in flash
            .vapor_mole_fractions()
            .iter()
            .zip([0.811_655_631_974_075_9, 0.188_344_368_025_924_16])
        {
            assert!((*actual - expected).abs() < 2e-8);
        }
        assert!(flash.log_k_residual() <= 1.0e-10);
        assert!(flash.material_balance_residual() < 1.0e-12);
        assert!(flash.gibbs_change_over_rt() < 0.0);
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
        let model =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.5, 0.5], vec![]).unwrap();
        assert!(matches!(
            model.two_phase_flash_pt(
                Temperature::new(400.0),
                Pressure::new(1.0e5),
                Default::default(),
            ),
            Err(PengRobinsonFlashErrorV1::NoTwoPhaseBracket { .. })
        ));
        assert!(matches!(
            model.two_phase_flash_pt(
                Temperature::new(220.0),
                Pressure::new(3.0e6),
                PengRobinsonFlashOptionsV1 {
                    maximum_iterations: 0,
                    ..Default::default()
                },
            ),
            Err(PengRobinsonFlashErrorV1::InvalidOptions)
        ));
    }

    #[test]
    fn g1_flash_recovers_split_outside_wilson_dew_bracket() {
        let model =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.028, 0.972], vec![])
                .unwrap();
        let t = Temperature::new(240.0);
        let p = Pressure::new(1.0e6);
        let k: Vec<_> = wilson_log_k(&model, t, p)
            .unwrap()
            .into_iter()
            .map(det::exp)
            .collect();
        assert!(rr_value(model.mole_fractions(), &k, 1.0).unwrap() > 0.0);
        let flash = model.two_phase_flash_pt(t, p, Default::default()).unwrap();
        // Independent binary calculation: solve ln(x_i phi_liq_i) =
        // ln(y_i phi_vap_i), then use the lever rule. No Wilson data in oracle.
        assert!((flash.liquid_mole_fractions()[0] - 0.003_716_176_809_946_948).abs() < 2e-10);
        assert!((flash.vapor_mole_fractions()[0] - 0.029_432_582_824_918_83).abs() < 2e-10);
        assert!((flash.vapor_fraction_mol() - 0.944_293_039_078_447).abs() < 2e-8);
        assert!(flash.gibbs_change_over_rt() < 0.0);
        assert!(flash.material_balance_residual() < 1e-12);
        for i in 0..2 {
            let liquid_potential = det::ln(flash.liquid_mole_fractions()[i])
                + flash.liquid().log_fugacity_coefficients()[i];
            let vapor_potential = det::ln(flash.vapor_mole_fractions()[i])
                + flash.vapor().log_fugacity_coefficients()[i];
            assert!((liquid_potential - vapor_potential).abs() < 2e-10);
        }
    }

    #[test]
    fn g0_tangent_plane_witness_detects_instability_of_feed_root() {
        let model =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.5, 0.5], vec![]).unwrap();
        let result = model
            .stability_pt(
                Temperature::new(220.0),
                Pressure::new(3e6),
                Default::default(),
            )
            .unwrap();
        let witness = result.instability_witness().unwrap();
        assert!(witness.tangent_plane_distance_over_rt() < -1e-3);
        let mut direct = 0.0;
        for i in 0..2 {
            let w = witness.mole_fractions()[i];
            direct += w
                * (det::ln(w / model.mole_fractions()[i])
                    + witness.root().log_fugacity_coefficients()[i]
                    - result.feed().log_fugacity_coefficients()[i]);
        }
        assert!((direct - witness.tangent_plane_distance_over_rt()).abs() < 2e-14);
        assert!(witness.stationarity_residual() <= 1e-10);

        let gas = model
            .stability_pt(
                Temperature::new(400.0),
                Pressure::new(1e5),
                Default::default(),
            )
            .unwrap();
        assert!(gas.instability_witness().is_none());
        for trial in [gas.vapor_seed(), gas.liquid_seed()] {
            assert!(trial.tangent_plane_distance_over_rt() >= -1e-10);
            assert!(trial.stationarity_residual() <= 1e-10);
        }
    }

    #[test]
    fn g3_phase_equilibrium_preserves_absences_and_component_order() {
        let plain =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.5, 0.5], vec![]).unwrap();
        let absent = PengRobinsonMixtureV1::new(
            vec![ethane(), methane(), methane()],
            vec![0.5, 0.5, 0.0],
            vec![],
        )
        .unwrap();
        let t = Temperature::new(220.0);
        let p = Pressure::new(3e6);
        let a = plain.two_phase_flash_pt(t, p, Default::default()).unwrap();
        let b = absent.two_phase_flash_pt(t, p, Default::default()).unwrap();
        assert!((a.vapor_fraction_mol() - b.vapor_fraction_mol()).abs() < 2e-12);
        assert!((a.liquid_mole_fractions()[0] - b.liquid_mole_fractions()[1]).abs() < 2e-12);
        assert_eq!(b.liquid_mole_fractions()[2], 0.0);
        assert_eq!(b.vapor_mole_fractions()[2], 0.0);
        let stability = absent.stability_pt(t, p, Default::default()).unwrap();
        assert!(stability.instability_witness().is_some());
        assert_eq!(stability.vapor_seed().mole_fractions()[2], 0.0);
        assert_eq!(stability.liquid_seed().mole_fractions()[2], 0.0);
    }

    #[test]
    fn g4_phase_searches_honor_cancellation_and_work_limits() {
        use super::super::stability::PengRobinsonStabilityOptionsV1;
        let model =
            PengRobinsonMixtureV1::new(vec![methane(), ethane()], vec![0.5, 0.5], vec![]).unwrap();
        let t = Temperature::new(220.0);
        let p = Pressure::new(3e6);
        let mut polls = 0;
        assert_eq!(
            model.two_phase_flash_pt_cancellable(t, p, Default::default(), || {
                polls += 1;
                polls == 3
            }),
            Err(PengRobinsonFlashErrorV1::Cancelled)
        );
        assert_eq!(polls, 3);
        assert_eq!(
            model.stability_pt_cancellable(t, p, Default::default(), || true),
            Err(PengRobinsonFlashErrorV1::Cancelled)
        );
        assert!(matches!(
            model.stability_pt(
                t,
                p,
                PengRobinsonStabilityOptionsV1 {
                    maximum_iterations: 1,
                    ..Default::default()
                }
            ),
            Err(PengRobinsonFlashErrorV1::StabilityIterationStalled { iterations: 1, .. })
        ));
    }
}
