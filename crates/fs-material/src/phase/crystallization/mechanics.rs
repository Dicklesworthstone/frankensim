//! Explicit stress-free locking and density-driven crystallization shrinkage.
//!
//! The caller chooses the locking temperature, crystal MASS fraction and total
//! strain; no universal percolation threshold or stiffness/fraction rule is
//! inferred. Subsequent free volume comes from the SAME original material pair.
//! Its isotropic linear change is `(v/v_lock)^(1/3)-1`, added to the locked
//! infinitesimal strain. This includes thermal density changes already present
//! in those sources; a second thermal expansion law must not be added.
//!
//! A frozen snapshot delegates stress and tangent to a caller-supplied original
//! stateless elastic law at elastic strain. It is directly consumable by native
//! plane/tetrahedral/hexahedral solvers, with no alternative mechanics solver.
//! Locking, current phase applicability and the explicit small-strain bounds
//! need independent material justification. Complete network unlocking has a
//! separate explicit fixed-strain transition in [`unlocking`]; no universal
//! melting criterion, plasticity, anisotropic crystallite growth, evolving
//! stiffness or pressure feedback is inferred. Updating snapshots supplies
//! one-way mechanics: the associated
//! elastic/parameter work is NOT included in the thermal crystallization heat
//! balance, and this adapter makes no closed coupled-energy claim.

use super::{CrystallizationCaloricPair, CrystallizationVolumeState};
use crate::phase::PhaseStateError;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_blake3::ContentHash;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Explicit loss of elastic support and its heat/energy-transfer account.
pub mod unlocking;

#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationMechanicsError {
    Phase(PhaseStateError),
    Invalid(&'static str),
    PairMismatch,
}
impl From<PhaseStateError> for CrystallizationMechanicsError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl core::fmt::Display for CrystallizationMechanicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization mechanics: {self:?}")
    }
}
impl std::error::Error for CrystallizationMechanicsError {}

/// Explicit component limits in (0,1), including sqrt(2)-weighted tensor shear
/// in Mandel. The caller must justify the actual infinitesimal-strain regime;
/// the upper admission boundary is not permission to treat order-one strain as linear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationStrainLimits {
    pub maximum_absolute_total_mandel_strain: f64,
    pub maximum_absolute_elastic_mandel_strain: f64,
    pub maximum_absolute_free_isotropic_strain: f64,
}

/// One caller-chosen stress-free reference. Its sources and history are immutable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationLock {
    volume: CrystallizationVolumeState,
    total_strain: Voigt,
}
impl CrystallizationLock {
    pub fn new(
        pair: &CrystallizationCaloricPair<'_>,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
        total_strain: Voigt,
    ) -> Result<Self, CrystallizationMechanicsError> {
        if total_strain.iter().any(|e| !e.is_finite()) {
            return Err(CrystallizationMechanicsError::Invalid(
                "nonfinite locking total strain",
            ));
        }
        Ok(Self {
            volume: pair.volume_state(temperature_k, crystalline_mass_fraction)?,
            total_strain,
        })
    }
    pub const fn caloric_pair_identity(&self) -> ContentHash {
        self.volume.caloric_pair_identity
    }
    pub const fn volume_state(&self) -> &CrystallizationVolumeState {
        &self.volume
    }
    pub const fn total_strain(&self) -> &Voigt {
        &self.total_strain
    }

    /// Prepare one fixed thermal/kinetic mechanical snapshot. The borrowed law
    /// must be elastic, valid at the current state, and stress free at zero
    /// elastic strain. Existing prestress or thermal eigenstrain is refused.
    pub fn freeze<'a, L: SmallStrainLaw<State = ()> + ?Sized>(
        &self,
        pair: &CrystallizationCaloricPair<'_>,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
        elastic: &'a L,
        limits: CrystallizationStrainLimits,
    ) -> Result<FrozenCrystallizationElastic<'a, L>, CrystallizationMechanicsError> {
        let invalid = CrystallizationMechanicsError::Invalid;
        if pair.identity() != self.volume.caloric_pair_identity {
            return Err(CrystallizationMechanicsError::PairMismatch);
        }
        if [
            limits.maximum_absolute_total_mandel_strain,
            limits.maximum_absolute_elastic_mandel_strain,
            limits.maximum_absolute_free_isotropic_strain,
        ]
        .into_iter()
        .any(|b| !b.is_finite() || b <= 0.0 || b >= 1.0)
        {
            return Err(invalid("finite positive strain limits below one required"));
        }
        if crystalline_mass_fraction < self.volume.crystalline_mass_fraction {
            return Err(invalid(
                "loss of crystallinity requires an explicit melting/unlocking law",
            ));
        }
        let declaration = elastic.admissibility();
        if !declaration.has_stored_energy
            || !declaration.tangent_symmetric
            || elastic
                .stress(&[0.0; 6], &())
                .iter()
                .any(|s| !s.is_finite() || *s != 0.0)
            || elastic
                .tangent(&[0.0; 6], &())
                .iter()
                .flatten()
                .any(|c| !c.is_finite())
        {
            return Err(invalid(
                "original stateless elastic law must have zero prestress and finite tangent",
            ));
        }
        let volume = pair.volume_state(temperature_k, crystalline_mass_fraction)?;
        let reference_v = self.volume.specific_volume_m3_kg;
        let ratio = volume.specific_volume_m3_kg / reference_v;
        let stretch = fs_math::det::pow(ratio, 1.0 / 3.0);
        // (r^(1/3)-1)=(r-1)/(r^(2/3)+r^(1/3)+1) retains tiny changes.
        let free_isotropic_strain = ((volume.specific_volume_m3_kg - reference_v) / reference_v)
            / (stretch * stretch + stretch + 1.0);
        if !stretch.is_finite()
            || stretch <= 0.0
            || !free_isotropic_strain.is_finite()
            || free_isotropic_strain.abs() > limits.maximum_absolute_free_isotropic_strain
        {
            return Err(invalid(
                "free density-driven strain leaves its declared small-strain limit",
            ));
        }
        let free_strain = std::array::from_fn(|i| {
            self.total_strain[i] + if i < 3 { free_isotropic_strain } else { 0.0 }
        });
        if !within(
            &self.total_strain,
            limits.maximum_absolute_total_mandel_strain,
        ) || !within(&free_strain, limits.maximum_absolute_total_mandel_strain)
        {
            return Err(invalid(
                "locked or free total strain leaves its declared limit",
            ));
        }
        Ok(FrozenCrystallizationElastic {
            elastic,
            lock: *self,
            volume,
            limits,
            free_strain,
            free_isotropic_strain,
        })
    }
}

/// Native constitutive snapshot. Temperature/history changes require a new snapshot.
pub struct FrozenCrystallizationElastic<'a, L: SmallStrainLaw<State = ()> + ?Sized> {
    elastic: &'a L,
    lock: CrystallizationLock,
    volume: CrystallizationVolumeState,
    limits: CrystallizationStrainLimits,
    free_strain: Voigt,
    free_isotropic_strain: f64,
}
impl<L: SmallStrainLaw<State = ()> + ?Sized> FrozenCrystallizationElastic<'_, L> {
    pub const fn source(&self) -> &L {
        self.elastic
    }
    pub const fn lock(&self) -> &CrystallizationLock {
        &self.lock
    }
    pub const fn volume_state(&self) -> &CrystallizationVolumeState {
        &self.volume
    }
    pub const fn free_strain(&self) -> &Voigt {
        &self.free_strain
    }
    pub const fn free_isotropic_strain(&self) -> f64 {
        self.free_isotropic_strain
    }
    /// Total and elastic strains both retain their separate supplied limits.
    pub fn elastic_strain(&self, total: &Voigt) -> Result<Voigt, CrystallizationMechanicsError> {
        let elastic = std::array::from_fn(|i| total[i] - self.free_strain[i]);
        if !within(total, self.limits.maximum_absolute_total_mandel_strain)
            || !within(&elastic, self.limits.maximum_absolute_elastic_mandel_strain)
        {
            return Err(CrystallizationMechanicsError::Invalid(
                "total or elastic strain leaves the admitted mechanical domain",
            ));
        }
        Ok(elastic)
    }
}
impl<L: SmallStrainLaw<State = ()> + ?Sized> SmallStrainLaw
    for FrozenCrystallizationElastic<'_, L>
{
    type State = ();
    fn initial_state(&self) {}
    fn stress(&self, total: &Voigt, (): &()) -> Voigt {
        match self.elastic_strain(total) {
            Ok(strain) => self.elastic.stress(&strain, &()),
            Err(_) => [f64::NAN; 6],
        }
    }
    fn tangent(&self, total: &Voigt, (): &()) -> Tangent6 {
        match self.elastic_strain(total) {
            Ok(strain) => self.elastic.tangent(&strain, &()),
            Err(_) => [[f64::NAN; 6]; 6],
        }
    }
    fn update_state(&self, _total: &Voigt, (): &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        let mut declaration = self.elastic.admissibility();
        declaration.failure_envelope = "original elastic domain plus declared total/elastic/free strain bounds; frozen history";
        declaration
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.frozen-crystallization-elastic",
            "0.1.0",
            Ambition::Frontier,
            vec![
                "explicit stress-free lock; original supplied elastic law".into(),
                "isotropic additive-specific-volume shrinkage; bounded infinitesimal strain".into(),
                format!(
                    "original caloric pair {:?}; T={} K; crystal mass fraction={}",
                    self.volume.caloric_pair_identity,
                    self.volume.temperature_k,
                    self.volume.crystalline_mass_fraction
                ),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                core::f64::consts::FRAC_1_SQRT_2 * self.limits.maximum_absolute_total_mandel_strain,
            ),
            vec![
                "one-way frozen snapshot; no automatic melting/reset or thermal feedback".into(),
                "locking and isotropic shrinkage have no calibrated model-discrepancy bound".into(),
            ],
            f64::INFINITY,
        )
    }
}

fn within(strain: &Voigt, limit: f64) -> bool {
    crate::tensor::to_mandel(strain)
        .iter()
        .all(|e| e.is_finite() && e.abs() <= limit)
}
