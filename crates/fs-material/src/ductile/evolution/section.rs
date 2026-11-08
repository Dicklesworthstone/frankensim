//! Progressive ductile work through the original plane-stress section.
//!
//! The section's existing reduction owns transverse equilibrium and its full
//! nonsymmetric tangent. This audit replays the single old-to-candidate native
//! return at each accepted thickness point; it never advances the candidate
//! again. Plastic dissipation, damage release and endpoint numerical loss keep
//! their original meanings. Neither a work deficit nor the work-scale diagnostic
//! is converted into heat, a crack surface or removed material.
//!
//! All work below is per original shell area. The material's immutable band
//! length is independent of thickness quadrature; this adapter neither selects
//! a crack direction nor substitutes shell thickness for that length.

use core::fmt::Debug;

use super::{
    DifferentiableDuctileInitiationLaw, DuctileEvolutionError, DuctileJ2Conditions, DuctileJ2Law,
    DuctileJ2State,
};
use crate::plastic::section::{LayeredPlaneStressSection, LayeredPlaneStressState, SectionError};

/// Physical storage and separate loss accounts per original section area, J/m2.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DuctileSectionEnergy {
    pub initial_stored_j_m2: f64,
    pub final_stored_j_m2: f64,
    /// Sum of pointwise storage changes, avoiding subtraction of large totals.
    pub stored_change_j_m2: f64,
    pub plastic_dissipation_j_m2: f64,
    pub damage_release_j_m2: f64,
    pub numerical_dissipation_j_m2: f64,
    pub material_work_j_m2: f64,
    pub generalized_work_j_m2: f64,
    /// Signed work of the small remaining transverse equilibrium residual.
    pub transverse_residual_work_j_m2: f64,
    /// Reference-thickness integral of `Gp * delta(Phi(x)) / h`, where
    /// Phi is the selected integrity profile's primitive. Only exponential
    /// attenuation makes delta(Phi) equal delta(damage). This is a diagnostic,
    /// NOT an additional loss term or work per newly created crack area.
    pub degraded_equivalent_plastic_work_j_m2: f64,
    /// Recovered `[Nx, Ny, Nxy, Mx, My, Mxy]`, in N/m and N respectively.
    pub resultants: [f64; 6],
    /// Generalized plus transverse work minus storage and the three losses.
    pub balance_residual_j_m2: f64,
}

/// A refused audit does not alter either endpoint or consume physical time.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileSectionEnergyError {
    Section(SectionError),
    Material {
        point: usize,
        source: DuctileEvolutionError,
    },
    State {
        point: usize,
        what: &'static str,
    },
    Work {
        residual_j_m2: f64,
        tolerance_j_m2: f64,
    },
}

impl core::fmt::Display for DuctileSectionEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ductile section energy refused: {self:?}")
    }
}
impl std::error::Error for DuctileSectionEnergyError {}
impl From<SectionError> for DuctileSectionEnergyError {
    fn from(value: SectionError) -> Self {
        Self::Section(value)
    }
}

impl<L> DuctileJ2Law<'_, L> {
    /// Actual interval and held temperature of this prepared constitutive law.
    /// Reading these does not retime the native return or rebind its history.
    pub const fn conditions(&self) -> DuctileJ2Conditions {
        self.conditions
    }
}

impl<'m, L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq>
    LayeredPlaneStressSection<DuctileJ2Law<'m, L>>
{
    /// Audit one complete accepted section increment using its exact original
    /// material owners, conditions and native energy allowances. All thickness
    /// points must describe this same single physical interval; a whole load
    /// path cannot be collapsed into one audit. Coefficients and band lengths
    /// remain bound to each point's immutable material owner.
    pub fn audit_ductile_energy(
        &self,
        old: &LayeredPlaneStressState<DuctileJ2State<'m, L>>,
        new: &LayeredPlaneStressState<DuctileJ2State<'m, L>>,
        tolerance_j_m2: f64,
    ) -> Result<DuctileSectionEnergy, DuctileSectionEnergyError> {
        self.validate_state(old)?;
        self.validate_state(new)?;
        if !tolerance_j_m2.is_finite() || tolerance_j_m2 <= 0.0 {
            return Err(SectionError::Invalid(
                "positive finite ductile section work allowance required",
            )
            .into());
        }
        let duration_s = self.layers()[0].law.conditions.duration_s;
        if self
            .layers()
            .iter()
            .any(|layer| layer.law.conditions.duration_s != duration_s)
        {
            return Err(SectionError::Invalid(
                "ductile section layers must share one physical interval",
            )
            .into());
        }
        let mut account = DuctileSectionEnergy::default();
        for (index, point) in self.points().iter().enumerate() {
            let law = &self.layers()[point.layer].law;
            let before = &old.points[index];
            let after = &new.points[index];
            if before.strain != *before.history.strain()
                || after.strain != *after.history.strain()
                || !law.owner.owns_state(&before.history)
                || !law.owner.owns_state(&after.history)
            {
                return Err(DuctileSectionEnergyError::State {
                    point: index,
                    what: "full strain and history must retain their original ductile owner",
                });
            }
            let trial = law
                .owner
                .trial(&after.strain, &before.history, law.conditions)
                .map_err(|source| DuctileSectionEnergyError::Material {
                    point: index,
                    source,
                })?;
            if trial.state != after.history {
                return Err(DuctileSectionEnergyError::State {
                    point: index,
                    what: "accepted history differs from its single original ductile return",
                });
            }
            let stress = trial.stress_pa;
            if stress.iter().any(|x| !x.is_finite())
                || [stress[2], stress[4], stress[5]]
                    .iter()
                    .any(|x| x.abs() > self.reduction_config().stress_tolerance_pa)
            {
                return Err(DuctileSectionEnergyError::State {
                    point: index,
                    what: "accepted transverse stress exceeds the original plane-stress allowance",
                });
            }
            let energy = trial.energy;
            let weight = point.weight_m;
            account.initial_stored_j_m2 += weight * energy.initial_stored_j_m3;
            account.final_stored_j_m2 += weight * energy.final_stored_j_m3;
            account.stored_change_j_m2 +=
                weight * (energy.final_stored_j_m3 - energy.initial_stored_j_m3);
            account.plastic_dissipation_j_m2 += weight * energy.plastic_dissipation_j_m3;
            account.damage_release_j_m2 += weight * energy.damage_release_j_m3;
            account.numerical_dissipation_j_m2 += weight * energy.numerical_dissipation_j_m3;
            account.material_work_j_m2 += weight * energy.mechanical_work_j_m3;
            account.degraded_equivalent_plastic_work_j_m2 += weight
                * (energy.degraded_equivalent_plastic_work_increment_j_m2
                    / law.owner.settings.characteristic_length_m);
            let transverse = stress[2] * (after.strain[2] - before.strain[2])
                + 2.0 * stress[4] * (after.strain[4] - before.strain[4])
                + 2.0 * stress[5] * (after.strain[5] - before.strain[5]);
            account.transverse_residual_work_j_m2 += weight * transverse;
            for (i, component) in [0, 1, 3].into_iter().enumerate() {
                account.resultants[i] += weight * stress[component];
                account.resultants[i + 3] += (weight * point.z_m) * stress[component];
            }
        }
        account.generalized_work_j_m2 = (0..6)
            .map(|i| {
                account.resultants[i] * (new.generalized_strain[i] - old.generalized_strain[i])
            })
            .sum();
        let work = account.generalized_work_j_m2 + account.transverse_residual_work_j_m2;
        account.balance_residual_j_m2 = work
            - account.stored_change_j_m2
            - account.plastic_dissipation_j_m2
            - account.damage_release_j_m2
            - account.numerical_dissipation_j_m2;
        if [
            account.initial_stored_j_m2,
            account.final_stored_j_m2,
            account.plastic_dissipation_j_m2,
            account.damage_release_j_m2,
            account.numerical_dissipation_j_m2,
            account.degraded_equivalent_plastic_work_j_m2,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || account.resultants.iter().any(|x| !x.is_finite())
        {
            return Err(SectionError::Invalid(
                "ductile section storage and losses must remain finite and nonnegative",
            )
            .into());
        }
        for residual in [
            work - account.material_work_j_m2,
            account.balance_residual_j_m2,
        ] {
            if !residual.is_finite() || residual.abs() > tolerance_j_m2 {
                return Err(DuctileSectionEnergyError::Work {
                    residual_j_m2: residual,
                    tolerance_j_m2,
                });
            }
        }
        Ok(account)
    }
}
