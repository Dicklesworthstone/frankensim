//! Voce isotropic saturation on the native associative J2 return map.
//!
//! The yield radius is `sigma_y + H alpha + Q (1-exp(-b alpha))`. The native
//! linear Prager backstress and the complete original `J2State` are retained.
//! A bounded scalar Newton solve determines the plastic increment; the original
//! tensor correction uses the saturation SECANT, while its algorithmic tangent
//! uses the ENDPOINT derivative. No tensor flow rule or tensor tangent is copied.
//!
//! The selected hardening is energetic: its additional recoverable storage is
//! `Q [alpha-(1-exp(-b alpha))/b]`. Only `sigma_y delta_alpha` is physical
//! dissipation. The convex hardening Bregman remainder is backward-Euler
//! numerical loss, never plastic heat. This is rate-independent small-strain
//! isotropic saturation, not nonlinear backstress recovery, damage or fatigue.

use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use super::{J2Plasticity, J2State, Update};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};
use fs_math::det;

/// Source-bound, fixed-condition Voce parameters and retained material receipts.
pub mod source;

/// Linear isotropic/Prager J2 with an independently selected Voce increment.
/// The zero-Q limit delegates bit-for-bit to the original native law.
#[derive(Clone, Debug, PartialEq)]
pub struct VoceJ2 {
    backbone: J2Plasticity,
    saturation_pa: f64,
    rate_per_strain: f64,
}

impl VoceJ2 {
    /// Admit Q >= 0 [Pa] and b > 0 [per equivalent plastic strain]. Q is the
    /// added saturation strength, not an absolute final yield stress. The
    /// original linear hardening remains active independently of saturation.
    pub fn new(
        backbone: J2Plasticity,
        saturation_pa: f64,
        rate_per_strain: f64,
    ) -> Result<Self, MaterialError> {
        backbone
            .validate_overstress_history(&[0.0; 6], &backbone.initial_state())
            .map_err(|e| MaterialError::Parameters {
                what: e.to_string(),
            })?;
        let (_, mu) = backbone.elastic.lame();
        let slope = saturation_pa * rate_per_strain;
        if !saturation_pa.is_finite()
            || saturation_pa < 0.0
            || !rate_per_strain.is_finite()
            || rate_per_strain <= 0.0
            || !slope.is_finite()
            || (saturation_pa > 0.0 && slope == 0.0)
            || !(backbone.yield_stress + saturation_pa).is_finite()
            || !(3.0 * mu + backbone.combined_hardening() + slope).is_finite()
        {
            return Err(MaterialError::Parameters {
                what: "Voce hardening needs finite Q >= 0, b > 0 and representable strength/slope"
                    .into(),
            });
        }
        Ok(Self {
            backbone,
            saturation_pa,
            rate_per_strain,
        })
    }

    /// Original elastic, initial-yield and linear hardening parameters.
    pub const fn backbone(&self) -> &J2Plasticity {
        &self.backbone
    }
    /// Additional asymptotic isotropic strength [Pa].
    pub const fn saturation_pa(&self) -> f64 {
        self.saturation_pa
    }
    /// Approach rate with respect to equivalent plastic strain.
    pub const fn rate_per_strain(&self) -> f64 {
        self.rate_per_strain
    }

    pub(super) fn saturation_strength(&self, alpha: f64) -> f64 {
        self.saturation_pa * -det::expm1(-self.rate_per_strain * alpha)
    }

    /// Stable incremental hardening slope, including the zero-increment limit.
    /// Shared by the rate-independent and combined overstress returns.
    pub(super) fn saturation_secant_pa(&self, alpha: f64, delta: f64) -> f64 {
        let x = self.rate_per_strain * delta;
        let ratio = if x == 0.0 { 1.0 } else { -det::expm1(-x) / x };
        self.saturation_pa * self.rate_per_strain
            * det::exp(-self.rate_per_strain * alpha) * ratio
    }

    /// Isotropic hardening force H alpha + Q (1-exp(-b alpha)) [Pa].
    pub fn isotropic_hardening_pa(&self, alpha: f64) -> f64 {
        if self.saturation_pa == 0.0 {
            return self.backbone.hardening * alpha;
        }
        self.backbone.hardening * alpha + self.saturation_strength(alpha)
    }

    /// Endpoint derivative H + Q b exp(-b alpha) [Pa]. The independent
    /// kinematic modulus is not an isotropic yield-radius derivative.
    pub fn isotropic_hardening_derivative_pa(&self, alpha: f64) -> f64 {
        if self.saturation_pa == 0.0 {
            return self.backbone.hardening;
        }
        self.backbone.hardening
            + self.saturation_pa * self.rate_per_strain * det::exp(-self.rate_per_strain * alpha)
    }

    /// Native effective-stress yield residual at the actual retained history.
    pub fn yield_function_with_state(&self, stress: &Voigt, state: &J2State) -> f64 {
        let original = self.backbone.yield_function_with_state(stress, state);
        if self.saturation_pa == 0.0 {
            original
        } else {
            original - self.saturation_strength(state.alpha)
        }
    }

    /// Admit a committed history at this law's actual temperature and yield
    /// radius. A thermal trial may subsequently soften that radius and yield.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &J2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), EnergyError> {
        if self.saturation_pa == 0.0 {
            return self
                .backbone
                .validate_committed_history(strain, state, yield_tolerance_pa);
        }
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa <= 0.0 {
            return Err(EnergyError::Invalid(
                "yield allowance must be finite and positive",
            ));
        }
        self.backbone.validate_overstress_history(strain, state)?;
        let elastic = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
        let stress = self.backbone.elastic.stress(&elastic, &());
        let residual = self.yield_function_with_state(&stress, state);
        if !residual.is_finite() {
            return Err(EnergyError::Invalid(
                "committed Voce stress is unrepresentable",
            ));
        }
        if residual > yield_tolerance_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: yield_tolerance_pa,
            });
        }
        Ok(())
    }

    /// Shared scalar consistency solve. The tensor owner supplies 3 mu; the
    /// exact axial restriction supplies E. This is an increment, not a new
    /// accumulated history. NaN identifies an unrepresentable raw response.
    pub(super) fn plastic_increment(
        &self,
        effective_trial_pa: f64,
        alpha: f64,
        elastic_correction_modulus_pa: f64,
    ) -> f64 {
        // Match the native shifted predictor's arithmetic grouping exactly.
        // A one-ulp disagreement at yield must not turn an elastic scalar
        // result into a failed tensor return.
        let radius = (self.backbone.yield_stress + self.saturation_strength(alpha))
            + self.backbone.hardening * alpha;
        let excess = effective_trial_pa - radius;
        if !excess.is_finite() || !alpha.is_finite() || alpha < 0.0 {
            return f64::NAN;
        }
        if excess <= 0.0 {
            return 0.0;
        }
        let a = elastic_correction_modulus_pa + self.backbone.combined_hardening();
        if !a.is_finite() || a <= 0.0 {
            return f64::NAN;
        }
        if self.saturation_pa == 0.0 {
            return excess / a;
        }
        let remaining = self.saturation_pa * det::exp(-self.rate_per_strain * alpha);
        if remaining == 0.0 {
            return excess / a;
        }
        let mut lower = 0.0;
        let mut upper = excess / a;
        let slope =
            self.saturation_pa * self.rate_per_strain * det::exp(-self.rate_per_strain * alpha);
        let mut delta = excess / (a + slope);
        if !upper.is_finite() || upper <= 0.0 || !delta.is_finite() || delta <= 0.0 {
            return f64::NAN;
        }
        for _ in 0..64 {
            let x = self.rate_per_strain * delta;
            let increment = remaining * -det::expm1(-x);
            let residual = a * delta + increment - excess;
            if !residual.is_finite() {
                return f64::NAN;
            }
            if residual.abs() <= 4.0 * f64::EPSILON * excess {
                return delta;
            }
            if residual < 0.0 {
                lower = delta;
            } else {
                upper = delta;
            }
            let derivative = a + slope * det::exp(-x);
            let newton = delta - residual / derivative;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == delta || next == lower || next == upper {
                return next;
            }
            delta = next;
        }
        f64::NAN
    }

    fn return_parameters(&self, strain: &Voigt, state: &J2State) -> (J2Plasticity, f64) {
        let mut shifted = self.backbone.clone();
        shifted.yield_stress += self.saturation_strength(state.alpha);
        let trial = shifted.radial_return(strain, state, 0.0);
        if !trial.plastic {
            return (shifted, 0.0);
        }
        let (_, mu) = self.backbone.elastic.lame();
        let delta = self.plastic_increment(trial.q_trial, state.alpha, 3.0 * mu);
        if !delta.is_finite() || delta <= 0.0 {
            return (shifted, f64::NAN);
        }
        let secant = self.saturation_secant_pa(state.alpha, delta);
        (shifted, secant)
    }

    fn returned(&self, strain: &Voigt, state: &J2State) -> Update {
        let (shifted, secant) = self.return_parameters(strain, state);
        shifted.radial_return(strain, state, secant)
    }

    // Integral_0^alpha Q (1-exp(-b p)) dp. The series preserves small-strain
    // storage when subtracting alpha and (1-exp(-b alpha))/b would cancel.
    pub(super) fn saturation_storage(&self, alpha: f64) -> f64 {
        self.saturation_storage_scaled(alpha, self.saturation_pa)
    }

    // d psi_Voce/dQ at fixed b and retained alpha, including Q=0. Do not
    // divide stored energy by Q: an affine thermal path can end at zero Q.
    pub(super) fn saturation_storage_per_pa(&self, alpha: f64) -> f64 {
        self.saturation_storage_scaled(alpha, 1.0)
    }

    fn saturation_storage_scaled(&self, alpha: f64, saturation_pa: f64) -> f64 {
        let x = self.rate_per_strain * alpha;
        if x < 0.1 {
            let ratio = -1.0 / 362880.0 + x / 3628800.0;
            let ratio = 1.0 / 40320.0 + x * ratio;
            let ratio = -1.0 / 5040.0 + x * ratio;
            let ratio = 1.0 / 720.0 + x * ratio;
            let ratio = -1.0 / 120.0 + x * ratio;
            let ratio = 1.0 / 24.0 + x * ratio;
            let ratio = -1.0 / 6.0 + x * ratio;
            let ratio = 0.5 + x * ratio;
            (saturation_pa * self.rate_per_strain) * alpha * (alpha * ratio)
        } else {
            saturation_pa * (alpha * (1.0 + det::expm1(-x) / x))
        }
    }

    // R(alpha_1) delta_alpha - Integral_alpha_0^alpha_1 R(p) dp >= 0.
    // Evaluating this convex remainder directly avoids subtracting two almost
    // equal hardening work/storage terms on small increments.
    pub(super) fn saturation_numerical_loss(&self, alpha: f64, delta: f64) -> f64 {
        if delta == 0.0 {
            return 0.0;
        }
        let x = self.rate_per_strain * delta;
        let ratio = if x < 0.1 {
            let ratio = -1.0 / 45360.0 + x / 403200.0;
            let ratio = 1.0 / 5760.0 + x * ratio;
            let ratio = -1.0 / 840.0 + x * ratio;
            let ratio = 1.0 / 144.0 + x * ratio;
            let ratio = -1.0 / 30.0 + x * ratio;
            let ratio = 1.0 / 8.0 + x * ratio;
            let ratio = -1.0 / 3.0 + x * ratio;
            x * (0.5 + x * ratio)
        } else {
            -det::expm1(-x) / x - det::exp(-x)
        };
        (self.saturation_pa * det::exp(-self.rate_per_strain * alpha)) * delta * ratio
    }

    /// Original native radial-return work account with the independent Voce
    /// storage and convex integration loss. Physical plastic heat remains
    /// sigma_y delta_alpha. Every strain/history gate and linear hardening
    /// account stays owned by the original J2 energy implementation.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        if self.saturation_pa == 0.0 {
            return self
                .backbone
                .advance_with_energy(old_strain, old_state, strain, tolerances);
        }
        let (shifted, secant) = self.return_parameters(strain, old_state);
        if !secant.is_finite() {
            return Err(EnergyError::Invalid("Voce return is not representable"));
        }
        // The temporary offset/secant only obtains the SAME native tensor
        // update and linear energy account. It is not physical viscosity.
        let mut next = shifted
            .advance_with_overstress_energy(old_strain, old_state, strain, tolerances, secant)?;
        let delta = next.state.alpha - old_state.alpha;
        let residual = self.yield_function_with_state(&next.stress_pa, &next.state);
        if !residual.is_finite() {
            return Err(EnergyError::Invalid("nonfinite Voce yield residual"));
        }
        if residual > tolerances.yield_pa || delta > 0.0 && residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        let initial = self.saturation_storage(old_state.alpha);
        let final_state = self.saturation_storage(next.state.alpha);
        let numerical = self.saturation_numerical_loss(old_state.alpha, delta);
        if [initial, final_state, numerical]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(EnergyError::Invalid(
                "Voce storage or numerical loss is not representable",
            ));
        }
        next.energy.initial.hardening_j_m3 += initial;
        next.energy.final_state.hardening_j_m3 += final_state;
        next.energy.physical_dissipation_j_m3 = self.backbone.yield_stress * delta;
        next.energy.numerical_dissipation_j_m3 += numerical;
        next.energy.residual_j_m3 = next.energy.mechanical_work_j_m3
            - next.energy.stored_change_j_m3()
            - next.energy.physical_dissipation_j_m3
            - next.energy.numerical_dissipation_j_m3;
        if !next.energy.initial.total_j_m3().is_finite()
            || !next.energy.final_state.total_j_m3().is_finite()
            || !next.energy.physical_dissipation_j_m3.is_finite()
            || !next.energy.numerical_dissipation_j_m3.is_finite()
            || !next.energy.residual_j_m3.is_finite()
        {
            return Err(EnergyError::Invalid(
                "Voce energy account is not representable",
            ));
        }
        if next.energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual {
                residual_j_m3: next.energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(next)
    }
}

impl SmallStrainLaw for VoceJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        if self.saturation_pa == 0.0 {
            return self.backbone.stress(strain, state);
        }
        self.returned(strain, state).stress
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        if self.saturation_pa == 0.0 {
            return self.backbone.tangent(strain, state);
        }
        let returned = self.returned(strain, state);
        let slope = self.saturation_pa
            * self.rate_per_strain
            * det::exp(-self.rate_per_strain * (state.alpha + returned.dgamma));
        self.backbone.tangent_from_return(&returned, slope)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        if self.saturation_pa == 0.0 {
            return self.backbone.update_state(strain, state);
        }
        let (shifted, secant) = self.return_parameters(strain, state);
        shifted.update_state_with_overstress(strain, state, secant)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "J2 yield radius sigma_y + H alpha + Q(1-exp(-b alpha)); linear Prager backstress",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        if self.saturation_pa == 0.0 {
            return self.backbone.card();
        }
        ModelCard::new(
            "material.j2-voce-plasticity",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain".into(),
                "associative rate-independent J2".into(),
                "energetic linear plus Voce isotropic hardening".into(),
                "independent linear Prager kinematic hardening".into(),
            ],
            self.backbone.card().validity,
            vec![
                "no nonlinear backstress recovery, damage or fatigue-life calibration".into(),
                "no finite-strain, creep or thermal-coefficient evolution".into(),
            ],
            0.05,
        )
    }
}
