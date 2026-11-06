//! Source-bound mass and mean forces for thermal tensile motion.
//!
//! Density joins the same original fracture and instantaneous-expansion source
//! bundle. An inertial consumer binds its chosen initial density snapshot once
//! to the original reference geometry and retains that mass inventory. Later
//! density queries do not authorize a mass change or a momentum impulse.
//!
//! The borrowed force law maps the total endpoint strain to mechanical strain
//! and uses the original analytic path-mean stress. After explicit coefficient
//! transfer, the constitutive path has fixed new coefficients and affine
//! mechanical strain: the difference of total displacement strain and free
//! thermal strain along their common path. Native endpoint history is retained.
//!
//! Mean stress is conjugate to both paths. Its contraction with total strain
//! change supplies mechanical work, and minus its trace times free-strain
//! change supplies eigenstrain work. Their sum is the original native material
//! work. This adapter cannot infer either split without the OLD thermal source;
//! the motion owner retains that source and performs the complete accounting.
//! No time integrator, heat conversion, extra damage law or second-order thermal
//! trajectory claim is supplied by the force adapter.

use super::{
    BrittleSourceError, ResolvedThermalBrittle, resolve_with_requirements,
    thermal_brittle_requirements,
};
use crate::brittle::path::BrittlePathStep;
use crate::brittle::{BrittleCrackBandSettings, BrittleState};
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Density;

/// Original fracture/elastic/expansion requirements plus positive density.
///
/// Static thermal consumers keep their existing requirement set; only a
/// consumer that binds reference mass needs the additional property.
pub fn thermal_brittle_dynamic_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = thermal_brittle_requirements()?;
    requirements.push(ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    Ok(requirements)
}

/// Resolve one original source bundle for thermally driven tensile motion.
///
/// Reference selection, integrated expansion and current-state refresh share
/// the static thermal resolver. All nine selected claims remain attached and
/// cannot change silently during a temperature successor. Density is positive
/// at reference and current states under the original source/domain admission.
///
/// The returned density value belongs to the CURRENT material snapshot. A
/// receiving motion owner explicitly binds its initial snapshot to reference
/// geometry once; this function does not infer a current/reference volume map.
pub fn resolve_thermal_brittle_dynamic_material(
    card: &MaterialCard,
    reference_query: &QueryPoint,
    current_query: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    let resolved = resolve_with_requirements(
        card,
        reference_query,
        current_query,
        selection,
        settings,
        &thermal_brittle_dynamic_requirements()?,
    )?;
    resolved.band().reference_density_kg_m3()?;
    Ok(resolved)
}

/// TOTAL-strain force view with the original native MECHANICAL-strain history.
///
/// Force stress is the exact native path mean. The accepted history and its
/// physical stress remain the endpoint response. A missing classical mean
/// derivative remains missing; an explicitly admitted native solver fallback
/// may differentiate this same force response.
#[derive(Clone, Copy, Debug)]
pub struct ThermalBrittleMeanSmallStrain<'a> {
    owner: &'a ResolvedThermalBrittle,
}

impl ResolvedThermalBrittle {
    /// Borrow the exact mean-force view for a coupled displacement/eigenstrain
    /// path. First transfer coefficients, then solve the thermal deformation;
    /// a prior restrained thermal update would commit a different physical path.
    #[must_use]
    pub fn mean_small_strain_law(&self) -> ThermalBrittleMeanSmallStrain<'_> {
        ThermalBrittleMeanSmallStrain { owner: self }
    }
}

impl<'a> ThermalBrittleMeanSmallStrain<'a> {
    /// Original endpoint, mean stress and optional derivative with respect to
    /// the NEW total tensor strain at fixed new temperature and old history.
    ///
    /// The additive eigenstrain has zero displacement derivative, so the mean
    /// derivative is exactly the original mechanical one. Returned path work is
    /// MECHANICAL-strain work, not total-displacement work alone. Total and
    /// mechanical strain bounds, owner identity and native work gates all apply.
    pub fn trial(
        &self,
        total: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<(BrittlePathStep<'a>, Option<Tangent6>), BrittleSourceError> {
        let mechanical = self.owner.mechanical_strain(total)?;
        Ok(self
            .owner
            .band()
            .law()
            .trial_with_mean_stress_tangent(&mechanical, state)?)
    }
}

impl<'a> SmallStrainLaw for ThermalBrittleMeanSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State {
        self.owner.small_strain_law().initial_state()
    }

    fn stress(&self, total: &Voigt, state: &Self::State) -> Voigt {
        // Force assembly must not pay for an unused mean-stress derivative.
        self.owner
            .mechanical_strain(total)
            .and_then(|mechanical| {
                self.owner
                    .band()
                    .law()
                    .trial_with_mean_stress(&mechanical, state)
                    .map_err(BrittleSourceError::from)
            })
            .map(|path| path.mean_stress_pa)
            .unwrap_or([f64::NAN; 6])
    }

    fn tangent(&self, total: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(total, state)
            .ok()
            .and_then(|(_, tangent)| tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }

    fn update_state(&self, total: &Voigt, state: &Self::State) -> Self::State {
        // Native endpoint evolution (including its refusal marker), never an
        // averaged history or a second material update.
        self.owner.small_strain_law().update_state(total, state)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }

    fn card(&self) -> ModelCard {
        let mut card = self.owner.small_strain_law().card();
        card.name = "material.brittle-source-thermal-path-mean".into();
        card.assumptions.push(
            "exact native mean stress along affine mechanical strain under fixed new coefficients; original endpoint history retained".into(),
        );
        card.assumptions.push(
            "total displacement and thermal eigenstrain evolve together after coefficient-only history transfer; their work ports sum to native mechanical-strain work".into(),
        );
        card.known_failures.push(
            "native mean tangent may be absent at nonsmooth contact/damage paths; the solver must explicitly admit its own fallback for this same force law".into(),
        );
        card.known_failures.push(
            "mean force stress differs from endpoint stress; no mass matrix, time integration or higher-order thermal-source trajectory is supplied".into(),
        );
        card.assumptions.sort_unstable();
        card.known_failures.sort_unstable();
        card
    }
}
