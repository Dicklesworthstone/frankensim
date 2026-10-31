//! Pure-fluid coexistence from equality of the original PR fugacities.
//!
//! At fixed T, `d(ln(phi_l)-ln(phi_v))/d(ln p) = Z_l-Z_v`. Safeguarded
//! Newton iteration closes coexistence without a fitted vapor-pressure law.
//! Single-phase bracket endpoints identify the liquid or vapor side; they
//! never become a fictitious second phase. PT alone leaves phase amount unset.
//!
//! The volume lever rule and the original departure functions supply phase
//! amount and latent calorics. Their original ideal reference must be added
//! ONCE by its source owner; identical phases have no ideal-mixing entropy.
//! Peng and Robinson (1976), <https://doi.org/10.1021/i160057a011>.
//! This considers fluid phases only, without solids, interfacial energy,
//! nucleation, chemical reactions or critical-region regularization.

use super::{
    PengRobinsonComponentV1, PengRobinsonErrorV1, PengRobinsonEvaluationV1, PengRobinsonMixtureV1,
    PengRobinsonPhaseRootV1,
};
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use crate::eos::IdealGasMolarVolumeQuantityV1 as MolarVolume;
use fs_math::det;
use fs_qty::{Pressure, Temperature};

// For x=v/b>1 the minimum of (x²+2x-1)²/[2(x+1)(x-1)²]
// lies at x³-3x²-3x-3=0. These separate the pure-PR spinodal branches.
const CRITICAL_VOLUME_OVER_B: f64 = 3.951_373_035_591_441;
const CRITICAL_ATTRACTION_OVER_RTB: f64 = 5.877_359_948_604_404;

/// Whole-call work and independent dimensionless convergence controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PengRobinsonSaturationOptionsV1 {
    /// Total EOS evaluations, including both bracket endpoints: 2 through 512.
    pub maximum_evaluations: usize,
    /// Maximum absolute liquid-minus-vapor logarithmic fugacity difference.
    pub log_fugacity_tolerance: f64,
    /// Maximum magnitude of the final Newton correction in log pressure.
    pub relative_pressure_tolerance: f64,
}
impl Default for PengRobinsonSaturationOptionsV1 {
    fn default() -> Self {
        Self {
            maximum_evaluations: 128,
            log_fugacity_tolerance: 1e-10,
            relative_pressure_tolerance: 1e-10,
        }
    }
}

/// No output candidate is published on domain, source or convergence refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum PengRobinsonSaturationErrorV1 {
    /// Refusal from the original EOS owner.
    Eos(PengRobinsonErrorV1),
    /// Invalid controls or missing physical constraint.
    Invalid(&'static str),
    /// Supplied pressure bounds do not bracket coexistence.
    NotBracketed,
    /// The two phases or their pressure sensitivity are unresolved.
    CriticalRegionUnresolved,
    /// Whole-call EOS budget exhausted before both gates closed.
    NoConvergence {
        /// Total EOS evaluations, including both bracket endpoints.
        evaluations: usize,
        /// Last two-phase logarithmic fugacity difference, if encountered.
        log_fugacity_residual: Option<f64>,
    },
    /// A finite physical quantity could not be represented.
    Unrepresentable(&'static str),
    /// Cancellation at an EOS boundary or before result publication.
    Cancelled,
}
impl From<PengRobinsonErrorV1> for PengRobinsonSaturationErrorV1 {
    fn from(value: PengRobinsonErrorV1) -> Self {
        Self::Eos(value)
    }
}
impl core::fmt::Display for PengRobinsonSaturationErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pure Peng-Robinson coexistence: {self:?}")
    }
}
impl std::error::Error for PengRobinsonSaturationErrorV1 {}
use PengRobinsonSaturationErrorV1 as E;

/// Distinct mechanically stable roots at common T, p and chemical potential.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonSaturationV1 {
    component: PengRobinsonComponentV1,
    temperature: Temperature,
    pressure: Pressure,
    liquid: PengRobinsonPhaseRootV1,
    vapor: PengRobinsonPhaseRootV1,
    log_fugacity_residual: f64,
    relative_pressure_correction: f64,
    evaluations: usize,
}
impl PengRobinsonSaturationV1 {
    /// Original pure-component parameters, including molar mass.
    #[must_use]
    pub const fn component(&self) -> &PengRobinsonComponentV1 {
        &self.component
    }
    /// Common phase temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Coexistence pressure computed from the original EOS.
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Saturated liquid root.
    #[must_use]
    pub const fn liquid(&self) -> &PengRobinsonPhaseRootV1 {
        &self.liquid
    }
    /// Saturated vapor root.
    #[must_use]
    pub const fn vapor(&self) -> &PengRobinsonPhaseRootV1 {
        &self.vapor
    }
    /// Signed liquid-minus-vapor logarithmic fugacity residual.
    #[must_use]
    pub const fn log_fugacity_residual(&self) -> f64 {
        self.log_fugacity_residual
    }
    /// Observed magnitude of the final Newton correction in log pressure.
    #[must_use]
    pub const fn relative_pressure_correction(&self) -> f64 {
        self.relative_pressure_correction
    }
    /// Total EOS evaluations, including bracket endpoints.
    #[must_use]
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    /// Vapor-minus-liquid enthalpy [J/mol]. The common ideal reference cancels.
    #[must_use]
    pub fn latent_enthalpy_j_mol(&self) -> f64 {
        self.vapor.residual_enthalpy_j_mol() - self.liquid.residual_enthalpy_j_mol()
    }
    /// Vapor-minus-liquid internal energy [J/mol], without a second heat debit.
    #[must_use]
    pub fn latent_internal_energy_j_mol(&self) -> f64 {
        self.vapor.residual_internal_energy_j_mol() - self.liquid.residual_internal_energy_j_mol()
    }
    /// Clapeyron slope `Delta h / (T Delta v)` [Pa/K].
    #[must_use]
    pub fn saturation_pressure_slope_pa_k(&self) -> f64 {
        self.latent_enthalpy_j_mol()
            / (self.temperature.value()
                * (self.vapor.molar_volume_m3_per_mol() - self.liquid.molar_volume_m3_per_mol()))
    }
    /// Resolve phase amount from an explicit mean molar volume. Endpoints are
    /// admitted; outside the coexistence dome refuses without clipping quality.
    pub fn at_molar_volume(
        &self,
        volume: MolarVolume,
    ) -> Result<PengRobinsonSaturatedMixtureV1, E> {
        let v = volume.value();
        let vl = self.liquid.molar_volume_m3_per_mol();
        let vv = self.vapor.molar_volume_m3_per_mol();
        if !v.is_finite() || v < vl || v > vv {
            return Err(E::Invalid(
                "molar volume must lie between the two saturation roots",
            ));
        }
        let result = self.at_vapor_fraction((v - vl) / (vv - vl))?;
        let residual = (result.molar_volume().value() - v).abs() / v;
        if !residual.is_finite() || residual > 128.0 * f64::EPSILON {
            return Err(E::Unrepresentable("two-phase volume balance"));
        }
        Ok(result)
    }
    /// Explicit vapor amount fraction. For one species mole and mass quality
    /// coincide; this is not volume fraction, and PT does not choose its value.
    pub fn at_vapor_fraction(
        &self,
        vapor_fraction: f64,
    ) -> Result<PengRobinsonSaturatedMixtureV1, E> {
        if !vapor_fraction.is_finite() || !(0.0..=1.0).contains(&vapor_fraction) {
            return Err(E::Invalid("vapor amount fraction must lie in [0,1]"));
        }
        Ok(PengRobinsonSaturatedMixtureV1 {
            saturation: self.clone(),
            vapor_fraction,
        })
    }
}

/// Two macroscopic equilibrium phases with amount-weighted extensive calorics.
/// No surface excess, slip, transport, heat-capacity or sound-speed law is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonSaturatedMixtureV1 {
    saturation: PengRobinsonSaturationV1,
    vapor_fraction: f64,
}
impl PengRobinsonSaturatedMixtureV1 {
    /// Original coexistence solution and component binding.
    #[must_use]
    pub const fn saturation(&self) -> &PengRobinsonSaturationV1 {
        &self.saturation
    }
    /// Vapor mole fraction, equal to mass quality for this pure component.
    #[must_use]
    pub const fn vapor_fraction(&self) -> f64 {
        self.vapor_fraction
    }
    fn mix(&self, liquid: f64, vapor: f64) -> f64 {
        (1.0 - self.vapor_fraction) * liquid + self.vapor_fraction * vapor
    }
    /// Additive molar volume [m³/mol].
    #[must_use]
    pub fn molar_volume(&self) -> MolarVolume {
        MolarVolume::new(self.mix(
            self.saturation.liquid.molar_volume_m3_per_mol(),
            self.saturation.vapor.molar_volume_m3_per_mol(),
        ))
    }
    /// Conserved mass divided by total phase volume, not averaged density.
    #[must_use]
    pub fn mass_density_kg_m3(&self) -> f64 {
        self.saturation.component.molar_mass.value() / self.molar_volume().value()
    }
    /// Enthalpy departure [J/mol]. Add the original ideal reference once.
    #[must_use]
    pub fn residual_enthalpy_j_mol(&self) -> f64 {
        self.mix(
            self.saturation.liquid.residual_enthalpy_j_mol(),
            self.saturation.vapor.residual_enthalpy_j_mol(),
        )
    }
    /// Internal-energy departure [J/mol], already including phase latent energy.
    #[must_use]
    pub fn residual_internal_energy_j_mol(&self) -> f64 {
        self.mix(
            self.saturation.liquid.residual_internal_energy_j_mol(),
            self.saturation.vapor.residual_internal_energy_j_mol(),
        )
    }
    /// Entropy departure from ideal pure gas at the same T,p [J/(mol K)].
    /// No ideal-mixing entropy is added between identical species.
    #[must_use]
    pub fn residual_entropy_j_mol_k(&self) -> f64 {
        let entropy = |root: &PengRobinsonPhaseRootV1| {
            root.residual_enthalpy_j_mol() / self.saturation.temperature.value()
                - R * root.log_fugacity_coefficients()[0]
        };
        self.mix(
            entropy(&self.saturation.liquid),
            entropy(&self.saturation.vapor),
        )
    }
}

struct Trial {
    evaluated: PengRobinsonEvaluationV1,
    gap: Option<f64>,
    side: f64,
}
impl Trial {
    fn new(evaluated: PengRobinsonEvaluationV1) -> Self {
        let roots = evaluated.phase_roots();
        let liquid = &roots[0];
        let vapor = &roots[roots.len() - 1];
        let gap = (roots.len() > 1)
            .then(|| liquid.log_fugacity_coefficients()[0] - vapor.log_fugacity_coefficients()[0]);
        let side = gap.unwrap_or_else(|| {
            if liquid.molar_volume_m3_per_mol() / liquid.covolume_m3_per_mol
                > CRITICAL_VOLUME_OVER_B
            {
                1.0
            } else {
                -1.0
            }
        });
        Self {
            evaluated,
            gap,
            side,
        }
    }
    fn log_pressure(&self) -> f64 {
        det::ln(self.evaluated.pressure().value())
    }
    fn correction(&self) -> Option<f64> {
        let roots = self.evaluated.phase_roots();
        self.gap.map(|gap| {
            gap / (roots[0].compressibility_factor()
                - roots[roots.len() - 1].compressibility_factor())
        })
    }
    fn accept(
        &self,
        component: PengRobinsonComponentV1,
        options: PengRobinsonSaturationOptionsV1,
        evaluations: usize,
    ) -> Result<Option<PengRobinsonSaturationV1>, E> {
        let Some(gap) = self.gap else {
            return Ok(None);
        };
        let correction = self.correction().unwrap_or(f64::NAN).abs();
        if gap.abs() > options.log_fugacity_tolerance {
            return Ok(None);
        }
        if !correction.is_finite() {
            return Err(E::CriticalRegionUnresolved);
        }
        if correction > options.relative_pressure_tolerance {
            return Ok(None);
        }
        let roots = self.evaluated.phase_roots();
        let liquid = &roots[0];
        let vapor = &roots[roots.len() - 1];
        let vl = liquid.molar_volume_m3_per_mol();
        let vv = vapor.molar_volume_m3_per_mol();
        if (vv - vl) / vv <= 64.0 * det::sqrt(f64::EPSILON)
            || [liquid, vapor].into_iter().any(|root| {
                let dpdv = root.pressure_derivatives().1;
                !dpdv.is_finite() || dpdv >= 0.0
            })
        {
            return Err(E::CriticalRegionUnresolved);
        }
        let result = PengRobinsonSaturationV1 {
            component,
            temperature: self.evaluated.temperature(),
            pressure: self.evaluated.pressure(),
            liquid: liquid.clone(),
            vapor: vapor.clone(),
            log_fugacity_residual: gap,
            relative_pressure_correction: correction,
            evaluations,
        };
        if [
            result.latent_enthalpy_j_mol(),
            result.latent_internal_energy_j_mol(),
            result.saturation_pressure_slope_pa_k(),
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
        {
            return Err(E::Unrepresentable(
                "positive latent energy and Clapeyron slope",
            ));
        }
        Ok(Some(result))
    }
}

fn evaluate_trial(
    model: &PengRobinsonMixtureV1,
    temperature: Temperature,
    pressure: Pressure,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial, E> {
    if cancelled() {
        return Err(E::Cancelled);
    }
    *evaluations += 1;
    let trial = Trial::new(model.evaluate_pt(temperature, pressure)?);
    if !trial.side.is_finite() {
        return Err(E::Unrepresentable("logarithmic fugacity difference"));
    }
    Ok(trial)
}

impl PengRobinsonMixtureV1 {
    /// Solve pure-fluid coexistence inside explicit pressure bounds at one T.
    pub fn pure_saturation_at_temperature(
        &self,
        temperature: Temperature,
        pressure_bracket: [Pressure; 2],
        options: PengRobinsonSaturationOptionsV1,
    ) -> Result<PengRobinsonSaturationV1, E> {
        self.pure_saturation_at_temperature_cancellable(
            temperature,
            pressure_bracket,
            options,
            || false,
        )
    }
    /// Cancellable bounded coexistence solve. Both tolerances must close;
    /// their observed residuals are numerical tests, not rigorous enclosures.
    /// Cancellation is checked before every EOS call and result publication.
    pub fn pure_saturation_at_temperature_cancellable(
        &self,
        temperature: Temperature,
        pressure_bracket: [Pressure; 2],
        options: PengRobinsonSaturationOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PengRobinsonSaturationV1, E> {
        if cancelled() {
            return Err(E::Cancelled);
        }
        if self.components.len() != 1 || self.mole_fractions != [1.0] {
            return Err(E::Invalid(
                "pure coexistence requires exactly one component with fraction one",
            ));
        }
        let component = self.components[0];
        let t = temperature.value();
        if !t.is_finite() || t <= 0.0 || t >= component.critical_temperature.value() {
            return Err(E::Invalid(
                "pure coexistence requires 0 < T < critical temperature",
            ));
        }
        let [pl, pu] = pressure_bracket.map(|pressure| pressure.value());
        if !pl.is_finite() || !pu.is_finite() || pl <= 0.0 || pl >= pu {
            return Err(E::Invalid(
                "positive finite ordered pressure bracket required",
            ));
        }
        if !(2..=512).contains(&options.maximum_evaluations)
            || [
                options.log_fugacity_tolerance,
                options.relative_pressure_tolerance,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0 || x >= 1.0)
        {
            return Err(E::Invalid("invalid coexistence work limit or tolerance"));
        }
        let mut evaluations = 0;
        let mut lower = evaluate_trial(
            self,
            temperature,
            pressure_bracket[0],
            &mut evaluations,
            &mut cancelled,
        )?;
        let attraction = lower.evaluated.dimensionless_a() / lower.evaluated.dimensionless_b();
        if !attraction.is_finite()
            || attraction <= CRITICAL_ATTRACTION_OVER_RTB * (1.0 + 128.0 * f64::EPSILON)
        {
            return Err(E::CriticalRegionUnresolved);
        }
        let mut upper = evaluate_trial(
            self,
            temperature,
            pressure_bracket[1],
            &mut evaluations,
            &mut cancelled,
        )?;
        for endpoint in [&lower, &upper] {
            if let Some(result) = endpoint.accept(component, options, evaluations)? {
                if cancelled() {
                    return Err(E::Cancelled);
                }
                return Ok(result);
            }
        }
        if lower.side < 0.0 || upper.side > 0.0 {
            return Err(E::NotBracketed);
        }
        let mut last_gap = upper.gap.or(lower.gap);
        while evaluations < options.maximum_evaluations {
            let lo = lower.log_pressure();
            let hi = upper.log_pressure();
            let best = if lower.gap.map_or(f64::INFINITY, f64::abs)
                <= upper.gap.map_or(f64::INFINITY, f64::abs)
            {
                &lower
            } else {
                &upper
            };
            let newton = best.log_pressure() - best.correction().unwrap_or(f64::NAN);
            let q = if newton.is_finite() && newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if q <= lo || q >= hi {
                break;
            }
            let trial = evaluate_trial(
                self,
                temperature,
                Pressure::new(det::exp(q)),
                &mut evaluations,
                &mut cancelled,
            )?;
            last_gap = trial.gap.or(last_gap);
            if let Some(result) = trial.accept(component, options, evaluations)? {
                if cancelled() {
                    return Err(E::Cancelled);
                }
                return Ok(result);
            }
            if trial.side >= 0.0 {
                lower = trial;
            } else {
                upper = trial;
            }
        }
        Err(E::NoConvergence {
            evaluations,
            log_fugacity_residual: last_gap,
        })
    }
}

#[cfg(test)]
mod tests;
