//! Associative linear overstress for the original anisotropic Hill return.
//!
//! With `f=q(sigma-X)-Y-R(alpha)`, the explicitly timed law is
//! `alpha_dot=max(f,0)/eta`, eta [Pa s]. Alpha is the native Hill
//! WORK-CONJUGATE multiplier, not generally von Mises equivalent strain.
//! Backward Euler adds `(eta/dt)*delta_alpha` only to the incremental
//! radius in the original five-mode scalar solve. The physical Y/H/Q/b/C
//! coefficients and stored energy remain unchanged. The consistent tangent
//! differentiates that same coupled solve, including nonlinear Voce saturation.
//!
//! Valid endpoints retain positive static overstress. Holding the deformation
//! for another actual interval can therefore relax stress and advance history.
//! Repeated candidate evaluations against the same incoming history do not
//! consume time; commit only after the surrounding physical step succeeds.
//!
//! Physical dissipation is `Y*delta_alpha + (eta/dt)*delta_alpha^2`.
//! Original elastic/hardening backward-Euler losses remain separate; no
//! viscosity is stored as hardening and no numerical loss is relabeled as heat.
//! The associative Perzyna family is described by
//! <https://neml.readthedocs.io/en/latest/vp_flow/perzyna.html>.
//! This explicit n=1 law supplies no sub-yield creep, nonlinear rate exponent,
//! dynamic recovery, finite deformation, evolving texture or fitted alloy data.

mod power;
pub use power::PowerLawOverstressHill48;

use super::oriented::OrientedHill48;
use super::{
    EnergyTolerances, Hill48, HillEnergyIncrement, HillError, HillResponse, HillState,
    LinearRateReturn, RateReturn, positive,
};
use crate::plastic::energy::StoredEnergy;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// A fixed physical interval with independently selected viscosity and material
/// axes. Local plastic history stays in the original immutable material frame.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearOverstressHill48 {
    backbone: OrientedHill48,
    viscosity_pa_s: f64,
    duration_s: f64,
    prepared: LinearRateReturn,
}

/// Structural history admission without a static-yield projection or an
/// invented time interval. Thermal callers use the coefficients at actual old T.
pub(super) fn validate_history(
    backbone: &OrientedHill48,
    strain: &Voigt,
    state: &HillState,
) -> Result<(), HillError> {
    backbone
        .law()
        .admit(&backbone.local_strain(strain), state)?;
    backbone.stored_energy(strain, state)?;
    Ok(())
}

impl LinearOverstressHill48 {
    /// Select finite positive eta [Pa s] and actual duration [s] in aligned axes.
    pub fn new(backbone: Hill48, viscosity_pa_s: f64, duration_s: f64) -> Result<Self, HillError> {
        Self::new_oriented(
            OrientedHill48::new(backbone, IDENTITY)?,
            viscosity_pa_s,
            duration_s,
        )
    }

    /// Preserve the original fixed material orientation and all physical
    /// hardening while preparing the kinetic denominator once per interval.
    pub fn new_oriented(
        backbone: OrientedHill48,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, HillError> {
        let modulus_pa = viscosity_pa_s / duration_s;
        if !positive(viscosity_pa_s) || !positive(duration_s) || !positive(modulus_pa) {
            return Err(HillError::Parameters(
                "Hill overstress needs finite positive eta, dt and representable eta/dt",
            ));
        }
        let law = backbone.law();
        let incremental_hardening = law.hardening + modulus_pa;
        let returns = law
            .modes
            .eigenvalues
            .map(|value| incremental_hardening + law.return_twice_shear * value);
        let scale = returns.iter().copied().fold(0.0_f64, f64::max)
            + law.saturation_pa * law.saturation_rate_per_strain;
        let ratios = returns.map(|value| value / scale);
        if !positive(scale) || ratios.iter().any(|value| !positive(*value)) {
            return Err(HillError::Parameters(
                "Hill viscous return moduli or ratios are not representable",
            ));
        }
        Ok(Self {
            backbone,
            viscosity_pa_s,
            duration_s,
            prepared: LinearRateReturn {
                modulus_pa,
                scale,
                ratios,
            },
        })
    }

    /// Original physical law and its immutable material orientation.
    pub const fn backbone(&self) -> &OrientedHill48 {
        &self.backbone
    }

    /// Original physical elastic, yield and hardening coefficients in local axes.
    pub const fn local_law(&self) -> &Hill48 {
        self.backbone.law()
    }

    /// Explicit viscosity against the native Hill work-conjugate multiplier.
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.viscosity_pa_s
    }

    /// The actual duration represented by one accepted backward-Euler update.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Endpoint excess above the STATIC surface for a proposed plastic increment.
    pub fn endpoint_overstress_pa(&self, increment: f64) -> Result<f64, HillError> {
        if !increment.is_finite() || increment < 0.0 {
            return Err(HillError::State(
                "Hill kinetic increment must be finite and nonnegative",
            ));
        }
        let value = self.prepared.modulus_pa * increment;
        if !value.is_finite() {
            return Err(HillError::Numerical(
                "Hill endpoint overstress is not representable",
            ));
        }
        Ok(value)
    }

    /// Admit the original finite, isochoric, reachable history. Positive static
    /// overstress is allowed and is not projected before a subsequent interval.
    pub fn validate_history(&self, strain: &Voigt, state: &HillState) -> Result<(), HillError> {
        validate_history(&self.backbone, strain, state)
    }

    /// Compatibility with accepted-history consumers. The stress allowance must
    /// be positive, but static yield is not an admissibility condition for this
    /// rate law; kinetic consistency belongs to a particular physical increment.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &HillState,
        tolerance_pa: f64,
    ) -> Result<(), HillError> {
        if !positive(tolerance_pa) {
            return Err(HillError::Parameters(
                "Hill history allowance must be finite and positive",
            ));
        }
        self.validate_history(strain, state)
    }

    /// One candidate with spatial stress/tangent and unchanged local history.
    /// `yield_residual_pa` reports the actual STATIC overstress, not the
    /// mismatch after subtracting eta*delta_alpha/dt.
    pub fn response(&self, strain: &Voigt, state: &HillState) -> Result<HillResponse, HillError> {
        self.backbone
            .spatial_response(self.local_law().response_with_rate(
                &self.backbone.local_strain(strain),
                state,
                Some(RateReturn::Linear(&self.prepared)),
            )?)
    }

    /// Original recoverable energy; no viscosity-dependent storage is added.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &HillState,
    ) -> Result<StoredEnergy, HillError> {
        self.backbone.stored_energy(strain, state)
    }

    /// Checked endpoint work, physical yield/viscous dissipation and the
    /// original elastic/isotropic/kinematic integration losses.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        self.advance_with_energy_from(&self.backbone, old_strain, old_state, strain, tolerances)
    }

    /// Endpoint-frozen coefficient work for the source-bound thermal owner.
    /// Previous history is admitted at its original coefficients; reciprocal
    /// thermal work and entropy remain that owner's responsibility.
    pub(crate) fn advance_with_energy_from(
        &self,
        previous: &OrientedHill48,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        if self.backbone.material_to_spatial() != previous.material_to_spatial() {
            return Err(HillError::Parameters(
                "thermal Hill history requires the original material frame",
            ));
        }
        let mut increment = self.local_law().advance_with_rate_energy_from(
            previous.law(),
            &previous.local_strain(old_strain),
            old_state,
            &self.backbone.local_strain(strain),
            tolerances,
            Some(RateReturn::Linear(&self.prepared)),
        )?;
        increment.strain = *strain;
        increment.response = self.backbone.spatial_response(increment.response)?;
        Ok(increment)
    }
}

impl SmallStrainLaw for LinearOverstressHill48 {
    type State = HillState;

    fn initial_state(&self) -> HillState {
        self.local_law().initial_state()
    }

    fn stress(&self, strain: &Voigt, state: &HillState) -> Voigt {
        // Line searches use the same return without building or rotating a tangent.
        let local = self.local_law().returned_with_rate(
            &self.backbone.local_strain(strain),
            state,
            Some(RateReturn::Linear(&self.prepared)),
        );
        let stress = local.map_or([f64::NAN; 6], |trial| {
            self.backbone.spatial_stress(&trial.stress_pa)
        });
        if stress.iter().all(|value| value.is_finite()) {
            stress
        } else {
            [f64::NAN; 6]
        }
    }

    fn tangent(&self, strain: &Voigt, state: &HillState) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |response| response.tangent_pa)
    }

    fn update_state(&self, strain: &Voigt, state: &HillState) -> HillState {
        self.local_law()
            .returned_with_rate(
                &self.backbone.local_strain(strain),
                state,
                Some(RateReturn::Linear(&self.prepared)),
            )
            .map_or_else(|_| state.clone(), |trial| trial.state)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "small-strain associative Hill linear overstress; positive static overstress is admitted",
            ..self.backbone.admissibility()
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.hill48-linear-overstress",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain and rotation; fixed material axes".into(),
                "associative anisotropic linear overstress at explicit physical dt".into(),
                "native Hill work-conjugate alpha and energetic linear/Voce/Prager hardening"
                    .into(),
                "exact coupled tangent; separate physical and integration dissipation".into(),
            ],
            self.backbone.card().validity,
            vec![
                "no sub-yield creep, nonlinear rate exponent or backstress recovery".into(),
                "no finite-strain forming, evolving texture, damage or measured-alloy calibration"
                    .into(),
            ],
            0.05,
        )
    }
}
