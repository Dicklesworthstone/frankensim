//! Reference-thickness work accounting for native energetic material laws.
//!
//! Audit the exact old/new histories of one already evaluated section. Native
//! `advance_with_energy` remains the sole constitutive energy owner; this
//! adapter neither derives heat from a work deficit nor advances an accepted
//! history. Transverse equilibrium residual work is retained with its sign.
//! Coefficients and their energetic interpretation must be unchanged over the
//! interval. This account does not transport heat or calibrate its partition.

use super::{LayeredPlaneStressSection, LayeredPlaneStressState, SectionError};
use crate::creep::plastic::source::ResolvedNortonJ2;
use crate::creep::plastic::{NortonJ2, NortonJ2State};
use crate::creep::recovery::RecoveringNortonCreep;
use crate::creep::recovery::plastic::RecoveringNortonJ2;
use crate::creep::recovery::plastic::source::ResolvedRecoveringNortonJ2;
use crate::creep::recovery::source::ResolvedRecoveringNorton;
use crate::creep::source::ResolvedNortonCreep;
use crate::creep::{CreepEnergy, CreepLaw, CreepTolerances, NortonCreep};
use crate::plastic::J2Plasticity;
use crate::plastic::chaboche::ChabocheJ2;
use crate::plastic::energy::{EnergyError, EnergyTolerances, IncrementEnergy, StoredEnergy};
use crate::plastic::hill::Hill48;
use crate::plastic::hill::oriented::OrientedHill48;
use crate::plastic::hill::source::ResolvedHill48;
use crate::plastic::rate::{LinearOverstressJ2, PowerLawOverstressJ2, VoceOverstressJ2};
use crate::plastic::tabulated::TabulatedJ2;
use crate::plastic::voce::VoceJ2;
use crate::{MaterialError, SmallStrainLaw, Voigt};

/// Explicit opt-in to a native mechanical energy account. Implementations must
/// use exactly the same update and physical duration as `SmallStrainLaw`.
pub trait NativeSectionEnergyLaw: SmallStrainLaw {
    /// Original physical interval of a timed native return. `None` is valid
    /// only for a rate-independent law. Motion and material must advance over
    /// the same interval; this metadata does not rescale or rebuild the law.
    fn physical_duration_s(&self) -> Option<f64> {
        None
    }

    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &Self::State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError>;
}

macro_rules! native_energy {
    ($law:ty, $validate:expr) => {
        native_energy!($law, $validate, |_| None);
    };
    ($law:ty, $validate:expr, $duration:expr) => {
        impl NativeSectionEnergyLaw for $law {
            fn physical_duration_s(&self) -> Option<f64> {
                let duration: fn(&$law) -> Option<f64> = $duration;
                duration(self)
            }

            fn native_energy_increment(
                &self,
                old_strain: &Voigt,
                old_state: &Self::State,
                strain: &Voigt,
                tolerances: EnergyTolerances,
            ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
                let validate: fn(&$law, &Voigt, &Self::State, f64) -> Result<(), EnergyError> =
                    $validate;
                validate(self, old_strain, old_state, tolerances.yield_pa).map_err(|e| {
                    MaterialError::State {
                        what: e.to_string(),
                    }
                })?;
                let value = self
                    .advance_with_energy(old_strain, old_state, strain, tolerances)
                    .map_err(|e| MaterialError::State {
                        what: e.to_string(),
                    })?;
                Ok((value.state, value.stress_pa, value.energy))
            }
        }
    };
}
native_energy!(J2Plasticity, J2Plasticity::validate_committed_history);
native_energy!(VoceJ2, VoceJ2::validate_committed_history);
native_energy!(TabulatedJ2, TabulatedJ2::validate_committed_history);
native_energy!(
    LinearOverstressJ2,
    |law, strain, state, _| law.backbone().validate_overstress_history(strain, state),
    |law| Some(law.duration_s())
);
native_energy!(
    PowerLawOverstressJ2,
    |law, strain, state, _| law.backbone().validate_overstress_history(strain, state),
    |law| Some(law.duration_s())
);
native_energy!(
    VoceOverstressJ2,
    |law, strain, state, _| law
        .backbone()
        .backbone()
        .validate_overstress_history(strain, state),
    |law| Some(law.duration_s())
);

// Hill's checked account owns incoming-history admission and the complete
// local-axis history. The oriented owner transforms work-conjugate stress and
// strain while retaining that history; no J2 projection or second return map
// belongs at the section boundary.
macro_rules! native_hill_energy {
    ($law:ty) => {
        impl NativeSectionEnergyLaw for $law {
            fn native_energy_increment(
                &self,
                old_strain: &Voigt,
                old_state: &Self::State,
                strain: &Voigt,
                tolerances: EnergyTolerances,
            ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
                let value = self
                    .advance_with_energy(old_strain, old_state, strain, tolerances)
                    .map_err(|e| MaterialError::State {
                        what: e.to_string(),
                    })?;
                Ok((value.response.state, value.response.stress_pa, value.energy))
            }
        }
    };
}
native_hill_energy!(Hill48);
native_hill_energy!(OrientedHill48);

// These explicit isothermal implementations keep thermal expansion/entropy
// out of the mechanical-only account. Primary strain hardening is kinetic:
// neither permanent creep nor resistance becomes recoverable hardening energy.
macro_rules! native_creep_energy {
    ($law:ty) => {
        impl NativeSectionEnergyLaw for $law {
            fn physical_duration_s(&self) -> Option<f64> {
                Some(self.kinetics().duration_s())
            }

            fn native_energy_increment(
                &self,
                old_strain: &Voigt,
                old_state: &Self::State,
                strain: &Voigt,
                tolerances: EnergyTolerances,
            ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
                let value = self
                    .creep_increment(old_strain, old_state, strain, creep_tolerances(tolerances))
                    .map_err(|e| MaterialError::State {
                        what: e.to_string(),
                    })?;
                let elastic_storage = |elastic_j_m3| StoredEnergy {
                    elastic_j_m3,
                    hardening_j_m3: 0.0,
                    kinematic_hardening_j_m3: 0.0,
                };
                let energy = creep_energy(
                    value.energy,
                    elastic_storage(value.energy.initial_stored_j_m3),
                    elastic_storage(value.energy.final_stored_j_m3),
                );
                Ok((value.state, value.stress_pa, energy))
            }
        }
    };
}
native_creep_energy!(NortonCreep);
native_creep_energy!(RecoveringNortonCreep);

impl NativeSectionEnergyLaw for NortonJ2 {
    fn physical_duration_s(&self) -> Option<f64> {
        Some(self.creep().duration_s())
    }

    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &Self::State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
        let value = self
            .advance_with_energy(old_strain, old_state, strain, creep_tolerances(tolerances))
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })?;
        let energy = creep_energy(
            value.energy,
            joint_creep_storage(self, old_strain, old_state)?,
            joint_creep_storage(self, strain, &value.state)?,
        );
        Ok((value.state, value.stress_pa, energy))
    }
}

impl NativeSectionEnergyLaw for RecoveringNortonJ2 {
    fn physical_duration_s(&self) -> Option<f64> {
        Some(self.kinetics().duration_s())
    }

    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &Self::State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
        let value = self
            .creep_increment(old_strain, old_state, strain, creep_tolerances(tolerances))
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })?;
        let energy = creep_energy(
            value.energy,
            joint_creep_storage(self.norton_j2(), old_strain, &old_state.native)?,
            joint_creep_storage(self.norton_j2(), strain, &value.state.native)?,
        );
        Ok((value.state, value.stress_pa, energy))
    }
}

// Query the original energetic J2 owner with the shared elastic strain.
// Subtracting H alpha_p^2/2 from a large total would lose this partition.
// Recovering kinetic resistance contributes no additional stored energy.
fn joint_creep_storage(
    law: &NortonJ2,
    total: &Voigt,
    state: &NortonJ2State,
) -> Result<StoredEnergy, MaterialError> {
    let elastic = core::array::from_fn(|i| {
        total[i] - state.creep.creep_strain[i] - state.plastic.plastic_strain[i]
    });
    law.plastic()
        .energy_of(&elastic, state.plastic.alpha, &state.plastic.plastic_strain)
        .map_err(|e| MaterialError::State {
            what: e.to_string(),
        })
}

fn creep_tolerances(value: EnergyTolerances) -> CreepTolerances {
    CreepTolerances {
        stress_pa: value.yield_pa,
        energy_j_m3: value.energy_j_m3,
    }
}

fn creep_energy(
    native: CreepEnergy,
    initial: StoredEnergy,
    final_state: StoredEnergy,
) -> IncrementEnergy {
    IncrementEnergy {
        initial,
        final_state,
        mechanical_work_j_m3: native.mechanical_work_j_m3,
        physical_dissipation_j_m3: native.physical_dissipation_j_m3,
        numerical_dissipation_j_m3: native.numerical_dissipation_j_m3,
        residual_j_m3: native.residual_j_m3,
    }
}

macro_rules! native_source_energy {
    ($law:ty) => {
        impl NativeSectionEnergyLaw for $law {
            fn physical_duration_s(&self) -> Option<f64> {
                self.law().physical_duration_s()
            }

            fn native_energy_increment(
                &self,
                old_strain: &Voigt,
                old_state: &Self::State,
                strain: &Voigt,
                tolerances: EnergyTolerances,
            ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
                self.law()
                    .native_energy_increment(old_strain, old_state, strain, tolerances)
            }
        }
    };
}
native_source_energy!(ResolvedNortonCreep);
native_source_energy!(ResolvedNortonJ2);
native_source_energy!(ResolvedRecoveringNorton);
native_source_energy!(ResolvedRecoveringNortonJ2);
native_source_energy!(ResolvedHill48);

impl<const N: usize> NativeSectionEnergyLaw for ChabocheJ2<N> {
    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &Self::State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
        self.validate_history(old_strain, old_state)
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })?;
        let elastic = core::array::from_fn(|i| old_strain[i] - old_state.plastic.plastic_strain[i]);
        let stress = self.backbone().backbone().elastic.stress(&elastic, &());
        let yield_residual = self.yield_function_with_state(&stress, old_state);
        if !yield_residual.is_finite() || yield_residual > tolerances.yield_pa {
            return Err(MaterialError::State {
                what: "initial Chaboche section history exceeds its yield allowance".into(),
            });
        }
        let value = self
            .advance_with_energy(old_strain, old_state, strain, tolerances)
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })?;
        // Native physical dissipation already contains dynamic recovery.
        Ok((value.state, value.stress_pa, value.energy))
    }
}

impl<L: NativeSectionEnergyLaw + ?Sized> NativeSectionEnergyLaw for &L {
    fn physical_duration_s(&self) -> Option<f64> {
        (**self).physical_duration_s()
    }

    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &Self::State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(Self::State, Voigt, IncrementEnergy), MaterialError> {
        (**self).native_energy_increment(old_strain, old_state, strain, tolerances)
    }
}

/// Recoverable native energy integrated through reference thickness, J/m².
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SectionStoredEnergy {
    pub elastic_j_m2: f64,
    pub hardening_j_m2: f64,
    pub kinematic_hardening_j_m2: f64,
}
impl SectionStoredEnergy {
    pub fn total_j_m2(self) -> f64 {
        self.elastic_j_m2 + self.hardening_j_m2 + self.kinematic_hardening_j_m2
    }
    fn add(&mut self, value: StoredEnergy, thickness: f64) {
        self.elastic_j_m2 += thickness * value.elastic_j_m3;
        self.hardening_j_m2 += thickness * value.hardening_j_m3;
        self.kinematic_hardening_j_m2 += thickness * value.kinematic_hardening_j_m3;
    }
}

/// Complete endpoint work in J/m², with the actual constrained residual work.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SectionEnergyAudit {
    pub initial: SectionStoredEnergy,
    pub final_state: SectionStoredEnergy,
    pub physical_dissipation_j_m2: f64,
    pub numerical_dissipation_j_m2: f64,
    pub material_work_j_m2: f64,
    pub generalized_work_j_m2: f64,
    /// sigma_zz delta_ezz + 2 sigma_yz delta_eyz + 2 sigma_zx delta_ezx.
    /// Signed local constraint residual work, never another heat source.
    pub transverse_residual_work_j_m2: f64,
    /// Generalized plus transverse work minus storage change and both losses.
    pub balance_residual_j_m2: f64,
    pub resultants: [f64; 6],
}
impl SectionEnergyAudit {
    pub fn stored_change_j_m2(&self) -> f64 {
        (self.final_state.elastic_j_m2 - self.initial.elastic_j_m2)
            + (self.final_state.hardening_j_m2 - self.initial.hardening_j_m2)
            + (self.final_state.kinematic_hardening_j_m2 - self.initial.kinematic_hardening_j_m2)
    }
}

/// Original state/native-energy refusal or an independently bounded work gap.
#[derive(Clone, Debug, PartialEq)]
pub enum SectionEnergyError {
    Section(SectionError),
    Native {
        point: usize,
        source: MaterialError,
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
impl core::fmt::Display for SectionEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native section energy refused: {self:?}")
    }
}
impl std::error::Error for SectionEnergyError {}
impl From<SectionError> for SectionEnergyError {
    fn from(error: SectionError) -> Self {
        Self::Section(error)
    }
}

impl<L: NativeSectionEnergyLaw> LayeredPlaneStressSection<L> {
    /// Audit exactly one old-to-new accepted return, without modifying either
    /// state. Multiple physical intervals cannot be collapsed into one audit.
    /// Timed laws replay their original single interval from the OLD history.
    pub fn audit_energy(
        &self,
        old: &LayeredPlaneStressState<L::State>,
        new: &LayeredPlaneStressState<L::State>,
        material: EnergyTolerances,
        tolerance_j_m2: f64,
    ) -> Result<SectionEnergyAudit, SectionEnergyError> {
        self.validate_state(old)?;
        self.validate_state(new)?;
        if [tolerance_j_m2, material.energy_j_m3, material.yield_pa]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(
                SectionError::Invalid("positive finite section work allowance required").into(),
            );
        }
        let mut account = SectionEnergyAudit::default();
        for (index, point) in self.points.iter().enumerate() {
            let before = &old.points[index];
            let after = &new.points[index];
            let (history, stress, energy) = self.layers[point.layer]
                .law
                .native_energy_increment(&before.strain, &before.history, &after.strain, material)
                .map_err(|source| SectionEnergyError::Native {
                    point: index,
                    source,
                })?;
            if history != after.history {
                return Err(SectionEnergyError::State {
                    point: index,
                    what: "accepted history does not equal the original native increment",
                });
            }
            if [stress[2], stress[4], stress[5]]
                .iter()
                .any(|x| !x.is_finite() || x.abs() > self.reduction.stress_tolerance_pa)
            {
                return Err(SectionEnergyError::State {
                    point: index,
                    what: "accepted transverse stress exceeds the plane-stress allowance",
                });
            }
            let weight = point.weight_m;
            account.initial.add(energy.initial, weight);
            account.final_state.add(energy.final_state, weight);
            account.physical_dissipation_j_m2 += weight * energy.physical_dissipation_j_m3;
            account.numerical_dissipation_j_m2 += weight * energy.numerical_dissipation_j_m3;
            account.material_work_j_m2 += weight * energy.mechanical_work_j_m3;
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
            - account.stored_change_j_m2()
            - account.physical_dissipation_j_m2
            - account.numerical_dissipation_j_m2;
        if [
            account.initial.total_j_m2(),
            account.final_state.total_j_m2(),
            account.physical_dissipation_j_m2,
            account.numerical_dissipation_j_m2,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(SectionError::Invalid(
                "section storage and losses must be finite and nonnegative",
            )
            .into());
        }
        for residual in [
            work - account.material_work_j_m2,
            account.balance_residual_j_m2,
        ] {
            if !residual.is_finite() || residual.abs() > tolerance_j_m2 {
                return Err(SectionEnergyError::Work {
                    residual_j_m2: residual,
                    tolerance_j_m2,
                });
            }
        }
        Ok(account)
    }
}
