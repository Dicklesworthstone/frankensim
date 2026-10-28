//! Linear associative J2 overstress on the original radial-return owner.
//!
//! With `f = q(stress-backstress) - sigma_y - H_iso alpha`, the declared Perzyna n=1 law is
//! `alpha_dot = max(f, 0)/eta`, eta in Pa s. Backward Euler gives
//! `delta_alpha = max(f_trial,0)/(3 mu + H_iso + H_k + eta/dt)`. The endpoint may
//! legitimately remain outside the rate-independent yield surface. The
//! original consistent tangent uses the same denominator. No artificial
//! hardening is stored: only the selected physical hardening terms enter free energy, while
//! `sigma_y delta_alpha + eta delta_alpha^2/dt` is physical dissipation.
//!
//! This is a bounded small-strain linear-overstress model, not arbitrary
//! rate dependence, finite-strain creep, damage, or a measured material law.

use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use super::{J2Plasticity, J2State};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_qty::Dims;

/// Positive linear-overstress viscosity in coherent Pa s, from material data.
pub const OVERSTRESS_VISCOSITY_PROPERTY: &str = "plastic_overstress_viscosity";
/// Pa s, independent of the dimensionally identical fluid shear-viscosity role.
pub const OVERSTRESS_VISCOSITY_DIMS: Dims = Dims([-1, 1, -1, 0, 0, 0]);

/// Native J2 coefficients plus one explicitly selected physical time interval.
/// A fresh value must be formed for each step; eta is not an added hardening law.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearOverstressJ2 {
    backbone: J2Plasticity,
    viscosity_pa_s: f64,
    duration_s: f64,
    overstress_modulus_pa: f64,
}

impl LinearOverstressJ2 {
    /// Admit a finite positive viscosity and actual physical step duration.
    pub fn new(
        backbone: J2Plasticity,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let modulus = viscosity_pa_s / duration_s;
        let (_, mu) = backbone.elastic.lame();
        if [
            viscosity_pa_s,
            duration_s,
            modulus,
            3.0 * mu + backbone.combined_hardening() + modulus,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(MaterialError::Parameters {
                what:
                    "linear overstress needs representable positive eta, dt and return denominator"
                        .into(),
            });
        }
        // Public backbone fields can have been modified; reuse its energy admission.
        backbone
            .validate_overstress_history(&[0.0; 6], &backbone.initial_state())
            .map_err(|e| MaterialError::Parameters {
                what: e.to_string(),
            })?;
        Ok(Self {
            backbone,
            viscosity_pa_s,
            duration_s,
            overstress_modulus_pa: modulus,
        })
    }

    /// The physical elastic, yield and energetic-hardening coefficients.
    pub const fn backbone(&self) -> &J2Plasticity {
        &self.backbone
    }
    /// Original source viscosity in Pa s.
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.viscosity_pa_s
    }
    /// The actual time represented by this backward-Euler return map.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Native return mapping and independently admitted work/heat balance.
    /// Accepted history may have positive overstress; it must remain a valid
    /// finite deviatoric plastic state within the same small-strain domain.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.backbone.advance_with_overstress_energy(
            old_strain,
            old_state,
            strain,
            tolerances,
            self.overstress_modulus_pa,
        )
    }
}

impl SmallStrainLaw for LinearOverstressJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.backbone
            .radial_return(strain, state, self.overstress_modulus_pa)
            .stress
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        self.backbone
            .tangent_with_overstress(strain, state, self.overstress_modulus_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        self.backbone
            .update_state_with_overstress(strain, state, self.overstress_modulus_pa)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "small-strain linear J2 overstress; positive f = eta alpha_dot is admitted",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-linear-overstress",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain".into(),
                "associative linear overstress".into(),
                if self.backbone.kinematic_hardening == 0.0 { "energetic linear isotropic hardening" }
                    else { "energetic linear isotropic and Prager hardening" }.into(),
                "backward Euler at explicit physical dt".into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                self.backbone.elastic.strain_limit.min(0.05),
            ),
            vec![
                "no nonlinear viscosity or rate exponent".into(),
                "no finite-strain creep or damage".into(),
            ],
            0.05,
        )
    }
}
