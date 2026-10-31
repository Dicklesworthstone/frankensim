//! Joint Voce hardening and associative Perzyna kinetics.
//!
//! For A = 3 mu + H_iso + H_k and the excess above the OLD Voce surface,
//! backward Euler solves A d + Q exp(-b alpha) (1-exp(-b d)) + f_v(d) = excess,
//! f_v(d) = s0 [eta d/(dt s0)]^(1/n). This is one monotone scalar return,
//! not sequential rate-independent and viscous corrections. The native tensor
//! return uses the summed secants; its tangent uses the endpoint derivatives.
//! Only (sigma_y + f_v) d dissipates physically. Voce energy and its convex
//! integration remainder retain their original storage/numerical owners.
//!
//! The associative Perzyna definition follows the same source as the native
//! power-law owner: https://neml.readthedocs.io/en/latest/vp_flow/perzyna.html .
//! Small-strain, energetic isotropic/Prager hardening only; no sub-yield creep,
//! nonlinear backstress recovery, finite-strain flow or material calibration.

use super::{EnergyError, EnergyIncrement, EnergyTolerances, PowerLawOverstressJ2};
use crate::plastic::voce::VoceJ2;
use crate::plastic::{J2Plasticity, J2State};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};
use fs_math::det;

/// A jointly integrated saturation law and explicitly timed physical kinetics.
/// Q=0 delegates exactly to the original selected linear or power-law owner.
#[derive(Clone, Debug, PartialEq)]
pub struct VoceOverstressJ2 {
    backbone: VoceJ2,
    kinetics: PowerLawOverstressJ2,
    linear_selection: bool,
    log_kinetic_coefficient: f64,
}

impl VoceOverstressJ2 {
    /// Select alpha_dot = max(f,0)/eta, eta [Pa s], at actual physical dt [s].
    pub fn new_linear(
        backbone: VoceJ2,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        // At n=1 the normalization cancels identically and is never queried
        // as a material property or exposed as a selected stress scale.
        let mut law = Self::new_power_law(backbone, viscosity_pa_s, 1.0, 1.0, duration_s)?;
        law.linear_selection = true;
        Ok(law)
    }

    /// Select alpha_dot = (s0/eta) max(f/s0,0)^n with n>=1 and positive s0.
    pub fn new_power_law(
        backbone: VoceJ2,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let kinetics = PowerLawOverstressJ2::new(
            backbone.backbone().clone(),
            viscosity_pa_s,
            stress_scale_pa,
            exponent,
            duration_s,
        )?;
        let (_, mu) = backbone.backbone().elastic.lame();
        let denominator = 3.0 * mu
            + backbone.backbone().combined_hardening()
            + backbone.saturation_pa() * backbone.rate_per_strain()
            + kinetics.linear.overstress_modulus_pa;
        if !denominator.is_finite() {
            return Err(MaterialError::Parameters {
                what: "combined Voce and overstress slope is not representable".into(),
            });
        }
        let log_kinetic_coefficient = det::ln(stress_scale_pa)
            + (det::ln(viscosity_pa_s) - det::ln(duration_s) - det::ln(stress_scale_pa)) / exponent;
        Ok(Self {
            backbone,
            kinetics,
            linear_selection: false,
            log_kinetic_coefficient,
        })
    }

    pub const fn backbone(&self) -> &VoceJ2 {
        &self.backbone
    }
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.kinetics.viscosity_pa_s()
    }
    pub const fn duration_s(&self) -> f64 {
        self.kinetics.duration_s()
    }
    pub const fn exponent(&self) -> f64 {
        self.kinetics.exponent()
    }
    /// No stress normalization is selected by the linear constructor.
    pub const fn stress_scale_pa(&self) -> Option<f64> {
        if self.linear_selection {
            None
        } else {
            Some(self.kinetics.stress_scale_pa())
        }
    }

    pub(crate) fn endpoint_overstress(&self, delta: f64) -> f64 {
        if delta == 0.0 {
            return 0.0;
        }
        if self.exponent() == 1.0 {
            return self.kinetics.linear.overstress_modulus_pa * delta;
        }
        det::exp(self.log_kinetic_coefficient + det::ln(delta) / self.exponent())
    }

    // Shared inverse kinetics for the native multi-branch consistency owner.
    // Logarithmic iterates avoid forming an overflowing stress power.
    pub(crate) fn overstress_from_log_increment(&self, log_delta: f64) -> f64 {
        if self.exponent() == 1.0 {
            self.kinetics.linear.overstress_modulus_pa * det::exp(log_delta)
        } else {
            det::exp(self.log_kinetic_coefficient + log_delta / self.exponent())
        }
    }

    pub(crate) const fn linear_overstress_modulus_pa(&self) -> f64 {
        self.kinetics.linear.overstress_modulus_pa
    }

    fn plastic_increment(&self, q_trial: f64, alpha: f64) -> f64 {
        let native = self.backbone.backbone();
        let (_, mu) = native.elastic.lame();
        if self.exponent() == 1.0 {
            return self.backbone.plastic_increment(
                q_trial,
                alpha,
                3.0 * mu + self.kinetics.linear.overstress_modulus_pa,
            );
        }
        let excess = q_trial
            - ((native.yield_stress + self.backbone.saturation_strength(alpha))
                + native.hardening * alpha);
        let a = 3.0 * mu + native.combined_hardening();
        if !excess.is_finite() || !alpha.is_finite() || alpha < 0.0 {
            return f64::NAN;
        }
        if excess <= 0.0 {
            return 0.0;
        }
        let upper_delta = excess / a;
        if !upper_delta.is_finite() || upper_delta <= 0.0 {
            return f64::NAN;
        }
        // A logarithmic bracket resolves extremely small rate increments
        // without powers of stress or an unbounded material-dependent loop.
        let mut lower = det::ln(f64::from_bits(1));
        let mut upper = det::ln(upper_delta);
        let evaluate = |log_delta: f64| {
            let delta = det::exp(log_delta);
            let viscous = det::exp(self.log_kinetic_coefficient + log_delta / self.exponent());
            let secant = self.backbone.saturation_secant_pa(alpha, delta);
            let residual = (a + secant) * delta + viscous - excess;
            let slope = self.backbone.saturation_pa()
                * self.backbone.rate_per_strain()
                * det::exp(-self.backbone.rate_per_strain() * (alpha + delta));
            (
                delta,
                residual,
                (a + slope) * delta + viscous / self.exponent(),
            )
        };
        if lower > upper || evaluate(lower).1 > 0.0 {
            return f64::NAN;
        }
        let mut position = upper;
        for _ in 0..96 {
            let (delta, residual, derivative) = evaluate(position);
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return f64::NAN;
            }
            if residual.abs() <= 8.0 * f64::EPSILON * excess {
                return delta;
            }
            if residual < 0.0 {
                lower = position;
            } else {
                upper = position;
            }
            let newton = position - residual / derivative;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == position || next == lower || next == upper {
                return det::exp(next);
            }
            position = next;
        }
        f64::NAN
    }

    fn return_parameters(&self, strain: &Voigt, state: &J2State) -> (J2Plasticity, f64) {
        let mut shifted = self.backbone.backbone().clone();
        shifted.yield_stress += self.backbone.saturation_strength(state.alpha);
        let trial = shifted.radial_return(strain, state, 0.0);
        if !trial.plastic {
            return (shifted, 0.0);
        }
        let delta = self.plastic_increment(trial.q_trial, state.alpha);
        if !delta.is_finite() || delta <= 0.0 {
            return (shifted, f64::NAN);
        }
        let viscous = if self.exponent() == 1.0 {
            self.kinetics.linear.overstress_modulus_pa
        } else {
            self.endpoint_overstress(delta) / delta
        };
        let secant = self.backbone.saturation_secant_pa(state.alpha, delta) + viscous;
        (shifted, secant)
    }

    /// Propose native history and independently check the joint rate equation
    /// and first law. No accepted history is modified on any return path.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        if self.backbone.saturation_pa() == 0.0 {
            return self
                .kinetics
                .advance_with_energy(old_strain, old_state, strain, tolerances);
        }
        let (shifted, secant) = self.return_parameters(strain, old_state);
        if !secant.is_finite() || secant < 0.0 {
            return Err(EnergyError::Invalid(
                "combined Voce-overstress return is not representable",
            ));
        }
        let mut next = shifted
            .advance_with_overstress_energy(old_strain, old_state, strain, tolerances, secant)?;
        let delta = next.state.alpha - old_state.alpha;
        let viscous = self.endpoint_overstress(delta);
        let residual = self
            .backbone
            .yield_function_with_state(&next.stress_pa, &next.state)
            - viscous;
        if !viscous.is_finite() || !residual.is_finite() {
            return Err(EnergyError::Invalid(
                "combined Voce-overstress endpoint is not representable",
            ));
        }
        if residual > tolerances.yield_pa || delta > 0.0 && residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        let initial = self.backbone.saturation_storage(old_state.alpha);
        let final_state = self.backbone.saturation_storage(next.state.alpha);
        let numerical = self
            .backbone
            .saturation_numerical_loss(old_state.alpha, delta);
        if [initial, final_state, numerical]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(EnergyError::Invalid(
                "combined Voce storage/loss is not representable",
            ));
        }
        next.energy.initial.hardening_j_m3 += initial;
        next.energy.final_state.hardening_j_m3 += final_state;
        next.energy.physical_dissipation_j_m3 =
            (self.backbone.backbone().yield_stress + viscous) * delta;
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
                "combined Voce-overstress energy is not representable",
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

impl SmallStrainLaw for VoceOverstressJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        if self.backbone.saturation_pa() == 0.0 {
            return self.kinetics.stress(strain, state);
        }
        let (shifted, secant) = self.return_parameters(strain, state);
        shifted.radial_return(strain, state, secant).stress
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        if self.backbone.saturation_pa() == 0.0 {
            return self.kinetics.tangent(strain, state);
        }
        let (shifted, secant) = self.return_parameters(strain, state);
        let returned = shifted.radial_return(strain, state, secant);
        let delta = returned.dgamma;
        let viscous_slope = if self.exponent() == 1.0 {
            self.kinetics.linear.overstress_modulus_pa
        } else if delta > 0.0 {
            (self.endpoint_overstress(delta) / delta) / self.exponent()
        } else {
            0.0
        };
        let saturation_slope = self.backbone.saturation_pa()
            * self.backbone.rate_per_strain()
            * det::exp(-self.backbone.rate_per_strain() * (state.alpha + delta));
        shifted.tangent_from_return(&returned, saturation_slope + viscous_slope)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        if self.backbone.saturation_pa() == 0.0 {
            return self.kinetics.update_state(strain, state);
        }
        let (shifted, secant) = self.return_parameters(strain, state);
        shifted.update_state_with_overstress(strain, state, secant)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        if self.backbone.saturation_pa() == 0.0 {
            return if self.linear_selection {
                self.kinetics.linear.admissibility()
            } else {
                self.kinetics.admissibility()
            };
        }
        MaterialAdmissibility {
            failure_envelope: "small-strain associative Voce/Prager J2 overstress, n >= 1",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        if self.backbone.saturation_pa() == 0.0 {
            return if self.linear_selection {
                self.kinetics.linear.card()
            } else {
                self.kinetics.card()
            };
        }
        ModelCard::new("material.j2-voce-overstress", "0.1.0", Ambition::Solid,
            vec!["small strain".into(), "joint associative Perzyna/Voce return, n >= 1".into(),
                "energetic isotropic and Prager hardening".into(), "backward Euler at explicit physical dt".into()],
            self.backbone.card().validity,
            vec!["no sub-yield creep, nonlinear backstress recovery, finite-strain flow or empirical calibration".into()], 0.05)
    }
}
