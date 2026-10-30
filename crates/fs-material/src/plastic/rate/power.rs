//! Power-law Perzyna kinetics on the existing J2 return map.
//!
//! For f = q(stress-backstress)-sigma_y-H_iso*alpha, the selected law is
//! alpha_dot = (s0/eta) max(f/s0,0)^n, eta [Pa s], s0 [Pa], n >= 1.
//! The normalization retains the original eta when n=1. Associative definition:
//! https://neml.readthedocs.io/en/latest/vp_flow/perzyna.html
//!
//! Backward Euler solves A*dalpha+s0*(eta*dalpha/(dt*s0))^(1/n)=f_trial,
//! A=3mu+H_iso+H_k. The positive secant f/dalpha supplies the original radial
//! return; its tangent instead uses df/dalpha=f/(n*dalpha). Only the endpoint
//! work (sigma_y+f)*dalpha heats. Viscosity stores no hardening energy.
//! This small-strain threshold law does not describe sub-yield creep,
//! nonlinear backstress saturation, finite-strain flow, or damage.

use super::{LinearOverstressJ2, J2Plasticity, J2State};
use super::{EnergyError, EnergyIncrement, EnergyTolerances};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Positive stress normalization [Pa], independently selected from material data.
pub const OVERSTRESS_STRESS_SCALE_PROPERTY: &str = "plastic_overstress_stress_scale";
/// Dimensionless power exponent, at least one throughout its source path.
pub const OVERSTRESS_EXPONENT_PROPERTY: &str = "plastic_overstress_exponent";

/// Nonlinear kinetics with the original J2 history and explicit physical time.
#[derive(Clone, Debug, PartialEq)]
pub struct PowerLawOverstressJ2 {
    linear: LinearOverstressJ2,
    stress_scale_pa: f64,
    exponent: f64,
}
impl PowerLawOverstressJ2 {
    /// Admit positive finite eta, s0, dt and finite n >= 1. The original
    /// small-strain energetic law retains the same coefficient admission.
    pub fn new(backbone: J2Plasticity, viscosity_pa_s: f64, stress_scale_pa: f64,
        exponent: f64, duration_s: f64) -> Result<Self, MaterialError> {
        if !stress_scale_pa.is_finite() || stress_scale_pa <= 0.0
            || !exponent.is_finite() || exponent < 1.0 {
            return Err(MaterialError::Parameters {
                what: "power-law overstress requires positive finite stress scale and exponent >= 1".into(),
            });
        }
        Ok(Self { linear: LinearOverstressJ2::new(backbone, viscosity_pa_s, duration_s)?,
            stress_scale_pa, exponent })
    }
    pub const fn backbone(&self) -> &J2Plasticity { self.linear.backbone() }
    pub const fn viscosity_pa_s(&self) -> f64 { self.linear.viscosity_pa_s() }
    pub const fn duration_s(&self) -> f64 { self.linear.duration_s() }
    pub const fn stress_scale_pa(&self) -> f64 { self.stress_scale_pa }
    pub const fn exponent(&self) -> f64 { self.exponent }

    // Solve log(m), m=f/dalpha, without overflowing f^n:
    // log(m)/n -(1-1/n)*log(1+A/m) = log(eta/dt)/n
    //                                  +(1-1/n)*log(s0/f_trial).
    // The left side is strictly increasing. A fixed 96-step bisection is
    // bounded independently of the material. Unrepresentable roots return NaN
    // so the existing mechanical/energy gates refuse the entire candidate.
    fn secant_modulus(&self, strain: &Voigt, state: &J2State) -> f64 {
        if self.exponent == 1.0 { return self.linear.overstress_modulus_pa; }
        let backbone = self.backbone();
        let trial = backbone.radial_return(strain, state, 0.0);
        if !trial.plastic { return 0.0; }
        let f_trial = trial.q_trial - (backbone.yield_stress + backbone.hardening * state.alpha);
        let (_, mu) = backbone.elastic.lame();
        let a = 3.0 * mu + backbone.combined_hardening();
        let inverse_n = 1.0 / self.exponent;
        let complement = 1.0 - inverse_n;
        let log_a = det::ln(a);
        let rhs = (det::ln(self.viscosity_pa_s()) - det::ln(self.duration_s())) * inverse_n
            + complement * (det::ln(self.stress_scale_pa) - det::ln(f_trial));
        let residual = |log_m: f64| log_m * inverse_n
            - complement * softplus(log_a - log_m) - rhs;
        let mut lower = det::ln(f64::from_bits(1));
        let mut upper = det::ln(f64::MAX);
        if !rhs.is_finite() || residual(lower) > 0.0 || residual(upper) < 0.0 {
            return f64::NAN;
        }
        for _ in 0..96 {
            let middle = 0.5 * lower + 0.5 * upper;
            if middle == lower || middle == upper { break; }
            if residual(middle) < 0.0 { lower = middle; } else { upper = middle; }
        }
        let modulus = det::exp(0.5 * lower + 0.5 * upper);
        if modulus > 0.0 && (a + modulus).is_finite() { modulus } else { f64::NAN }
    }

    fn endpoint_overstress(&self, delta_alpha: f64) -> f64 {
        if delta_alpha == 0.0 { return 0.0; }
        det::exp(det::ln(self.stress_scale_pa)
            + (det::ln(self.viscosity_pa_s()) + det::ln(delta_alpha)
                - det::ln(self.duration_s()) - det::ln(self.stress_scale_pa)) / self.exponent)
    }

    /// Propose the original native history and independently audit the nonlinear
    /// flow equation and physical heat, without mutating accepted state.
    pub fn advance_with_energy(&self, old_strain: &Voigt, old_state: &J2State,
        strain: &Voigt, tolerances: EnergyTolerances) -> Result<EnergyIncrement, EnergyError> {
        if self.exponent == 1.0 {
            return self.linear.advance_with_energy(old_strain, old_state, strain, tolerances);
        }
        let modulus = self.secant_modulus(strain, old_state);
        if !modulus.is_finite() {
            return Err(EnergyError::Invalid("power-law overstress return is not representable"));
        }
        let mut next = self.backbone().advance_with_overstress_energy(old_strain, old_state,
            strain, tolerances, modulus)?;
        let delta = next.state.alpha - old_state.alpha;
        let overstress = self.endpoint_overstress(delta);
        let residual = self.backbone().yield_function_with_state(&next.stress_pa, &next.state) - overstress;
        if !overstress.is_finite() || !residual.is_finite() {
            return Err(EnergyError::Invalid("power-law endpoint rate is not representable"));
        }
        if residual > tolerances.yield_pa || delta > 0.0 && residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual { residual_pa: residual, tolerance_pa: tolerances.yield_pa });
        }
        // The secant is only an implementation device. Deposit the independently
        // evaluated physical law and retain every original storage/numerical
        // account, then recheck the complete first-law identity.
        next.energy.physical_dissipation_j_m3 = self.backbone().yield_stress * delta + overstress * delta;
        next.energy.residual_j_m3 = next.energy.mechanical_work_j_m3 - next.energy.stored_change_j_m3()
            - next.energy.physical_dissipation_j_m3 - next.energy.numerical_dissipation_j_m3;
        if !next.energy.residual_j_m3.is_finite() || next.energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual { residual_j_m3: next.energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3 });
        }
        Ok(next)
    }
}

fn softplus(x: f64) -> f64 {
    if x > 0.0 { x + det::ln(1.0 + det::exp(-x)) }
    else { det::ln(1.0 + det::exp(x)) }
}

impl SmallStrainLaw for PowerLawOverstressJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State { self.linear.initial_state() }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        if self.exponent == 1.0 { return self.linear.stress(strain, state); }
        self.backbone().radial_return(strain, state, self.secant_modulus(strain, state)).stress
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        if self.exponent == 1.0 { return self.linear.tangent(strain, state); }
        let modulus = self.secant_modulus(strain, state);
        let update = self.backbone().radial_return(strain, state, modulus);
        self.backbone().tangent_from_return(&update, modulus / self.exponent)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        if self.exponent == 1.0 { return self.linear.update_state(strain, state); }
        self.backbone().update_state_with_overstress(strain, state, self.secant_modulus(strain, state))
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility { failure_envelope: "small-strain associative J2 power-law overstress, n >= 1",
            ..self.backbone().admissibility() }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new("material.j2-power-law-overstress", "0.1.0", Ambition::Solid,
            vec!["small strain".into(), "associative power-law overstress n >= 1".into(),
                "energetic isotropic/Prager hardening".into(), "backward Euler at explicit physical dt".into()],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0,
                self.backbone().elastic.strain_limit.min(0.05)),
            vec!["no sub-yield creep, finite-strain flow, damage or measured qualification".into()], 0.05)
    }
}
