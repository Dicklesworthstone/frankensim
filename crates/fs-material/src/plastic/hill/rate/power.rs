//! Nonlinear associative overstress on the original five-mode Hill return.
//!
//! With f=q(sigma-X)-Y-R(alpha), the selected kinetic law is
//! alpha_dot=(s0/eta)*max(f/s0,0)^n, eta [Pa s], s0 [Pa], n>=1.
//! Alpha keeps the original Hill work-conjugate normalization. Backward Euler
//! adds f_rate=s0*(eta*delta_alpha/(dt*s0))^(1/n) to the endpoint radius.
//! Its tangent uses f_rate/(n*delta_alpha), not a secant or finite difference.
//!
//! The concave resistance contributes f_rate*(1-1/n) to y-dalpha*y' in the
//! native rank-one tangent. One monotone scalar solve therefore retains the
//! anisotropic direction, Voce saturation and Prager hardening together.
//! Powers use dimensionless logarithms and analytic brackets; n=1 calls the
//! existing linear path exactly. Point evaluation allocates no memory.
//!
//! Physical heat is (Y+f_rate)*delta_alpha. No kinetic coefficient enters
//! stored energy or numerical dissipation. Retained overstress is admitted;
//! each successful hold represents another explicitly supplied interval.
//! Reference: <https://neml.readthedocs.io/en/latest/vp_flow/perzyna.html>.

use super::super::oriented::OrientedHill48;
use super::super::{
    EnergyTolerances, Hill48, HillEnergyIncrement, HillError, HillResponse, HillState,
    PowerRateReturn, RateReturn, log_ratio, positive,
};
use super::{LinearOverstressHill48, validate_history};
use crate::plastic::energy::StoredEnergy;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// A power-law kinetic selection with fixed physical coefficients, original
/// material axes and one actual backward-Euler interval.
#[derive(Clone, Debug, PartialEq)]
pub struct PowerLawOverstressHill48 {
    linear: LinearOverstressHill48,
    stress_scale_pa: f64,
    exponent: f64,
    prepared: PowerRateReturn,
}

impl PowerLawOverstressHill48 {
    /// Select eta [Pa s], s0 [Pa], n>=1 and physical dt [s] in aligned axes.
    /// The original physical Hill, Voce and Prager coefficients are unchanged.
    pub fn new(
        backbone: Hill48,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, HillError> {
        Self::new_oriented(
            OrientedHill48::new(backbone, IDENTITY)?,
            viscosity_pa_s,
            stress_scale_pa,
            exponent,
            duration_s,
        )
    }

    /// Preserve the original immutable material frame and local history.
    pub fn new_oriented(
        backbone: OrientedHill48,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, HillError> {
        if !positive(stress_scale_pa) || !exponent.is_finite() || exponent < 1.0 {
            return Err(HillError::Parameters(
                "Hill power-law overstress needs positive finite s0 and finite n>=1",
            ));
        }
        let linear = LinearOverstressHill48::new_oriented(backbone, viscosity_pa_s, duration_s)?;
        let modulus_pa = linear.prepared.modulus_pa;
        let prepared = PowerRateReturn {
            modulus_pa,
            stress_scale_pa,
            inverse_exponent: 1.0 / exponent,
            log_modulus_ratio: log_ratio(modulus_pa, linear.local_law().return_scale),
        };
        Ok(Self {
            linear,
            stress_scale_pa,
            exponent,
            prepared,
        })
    }

    /// Original physical law and fixed material orientation.
    pub const fn backbone(&self) -> &OrientedHill48 {
        self.linear.backbone()
    }
    /// Original elastic, yield and energetic hardening in material axes.
    pub const fn local_law(&self) -> &Hill48 {
        self.linear.local_law()
    }
    /// Independently selected viscosity against the Hill multiplier [Pa s].
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.linear.viscosity_pa_s()
    }
    /// Independent kinetic stress normalization [Pa].
    pub const fn stress_scale_pa(&self) -> f64 {
        self.stress_scale_pa
    }
    /// Dimensionless forward-flow exponent, at least one.
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }
    /// Actual duration represented by one accepted candidate [s].
    pub const fn duration_s(&self) -> f64 {
        self.linear.duration_s()
    }

    /// Actual excess above the static yield surface for a proposed increment.
    pub fn endpoint_overstress_pa(&self, increment: f64) -> Result<f64, HillError> {
        if self.exponent == 1.0 {
            return self.linear.endpoint_overstress_pa(increment);
        }
        if !increment.is_finite() || increment < 0.0 {
            return Err(HillError::State(
                "Hill kinetic increment must be finite and nonnegative",
            ));
        }
        let value = self.prepared.overstress(increment);
        if !value.is_finite() {
            return Err(HillError::Numerical(
                "Hill nonlinear overstress is not representable",
            ));
        }
        Ok(value)
    }

    /// Structural, finite, isochoric and reachable history admission; the
    /// static surface does not restrict an accepted finite-rate endpoint.
    pub fn validate_history(&self, strain: &Voigt, state: &HillState) -> Result<(), HillError> {
        validate_history(self.backbone(), strain, state)
    }

    /// Accepted-history consumer seam. Kinetic consistency belongs to an
    /// explicit interval, not a duration-free static-yield projection.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &HillState,
        tolerance_pa: f64,
    ) -> Result<(), HillError> {
        self.linear
            .validate_committed_history(strain, state, tolerance_pa)
    }

    /// Spatial stress and exact tangent from one frozen-history candidate.
    /// yield_residual_pa remains the actual static overstress.
    pub fn response(&self, strain: &Voigt, state: &HillState) -> Result<HillResponse, HillError> {
        if self.exponent == 1.0 {
            return self.linear.response(strain, state);
        }
        self.backbone()
            .spatial_response(self.local_law().response_with_rate(
                &self.backbone().local_strain(strain),
                state,
                Some(RateReturn::Power(&self.prepared)),
            )?)
    }

    /// Original recoverable energy with no kinetic storage.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &HillState,
    ) -> Result<StoredEnergy, HillError> {
        self.backbone().stored_energy(strain, state)
    }

    /// Check nonlinear flow, retained state, work, physical heat and the
    /// original material's distinct integration losses without committing.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        self.advance_with_energy_from(self.backbone(), old_strain, old_state, strain, tolerances)
    }

    /// Endpoint-frozen coefficient seam for the original thermal source owner.
    /// Previous history is admitted at its actual old coefficients and frame.
    pub(crate) fn advance_with_energy_from(
        &self,
        previous: &OrientedHill48,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        if self.exponent == 1.0 {
            return self
                .linear
                .advance_with_energy_from(previous, old_strain, old_state, strain, tolerances);
        }
        if self.backbone().material_to_spatial() != previous.material_to_spatial() {
            return Err(HillError::Parameters(
                "thermal Hill history requires the original material frame",
            ));
        }
        let mut increment = self.local_law().advance_with_rate_energy_from(
            previous.law(),
            &previous.local_strain(old_strain),
            old_state,
            &self.backbone().local_strain(strain),
            tolerances,
            Some(RateReturn::Power(&self.prepared)),
        )?;
        increment.strain = *strain;
        increment.response = self.backbone().spatial_response(increment.response)?;
        Ok(increment)
    }
}

impl SmallStrainLaw for PowerLawOverstressHill48 {
    type State = HillState;
    fn initial_state(&self) -> HillState {
        self.linear.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &HillState) -> Voigt {
        if self.exponent == 1.0 {
            return self.linear.stress(strain, state);
        }
        let local = self.local_law().returned_with_rate(
            &self.backbone().local_strain(strain),
            state,
            Some(RateReturn::Power(&self.prepared)),
        );
        let stress = local.map_or([f64::NAN; 6], |trial| {
            self.backbone().spatial_stress(&trial.stress_pa)
        });
        if stress.iter().all(|value| value.is_finite()) {
            stress
        } else {
            [f64::NAN; 6]
        }
    }
    fn tangent(&self, strain: &Voigt, state: &HillState) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &HillState) -> HillState {
        if self.exponent == 1.0 {
            return self.linear.update_state(strain, state);
        }
        self.local_law()
            .returned_with_rate(
                &self.backbone().local_strain(strain),
                state,
                Some(RateReturn::Power(&self.prepared)),
            )
            .map_or_else(|_| state.clone(), |trial| trial.state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "small-strain associative Hill power-law overstress; n>=1 and positive static overstress admitted",
            ..self.backbone().admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.hill48-power-law-overstress",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain and rotation; fixed material axes".into(),
                "associative anisotropic power-law overstress n>=1 at explicit dt".into(),
                "native Hill alpha with energetic linear/Voce/Prager hardening".into(),
                "exact coupled tangent and distinct physical/integration dissipation".into(),
            ],
            self.backbone().card().validity,
            vec![
                "no sub-yield creep, dynamic backstress recovery or evolving texture".into(),
                "no finite-strain forming, damage or measured-alloy calibration".into(),
            ],
            0.05,
        )
    }
}
