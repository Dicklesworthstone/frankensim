//! Same-card measured Cauchy hardening and explicit Perzyna kinetics.
//!
//! All coefficients are resolved at one fixed original source context. The
//! physical interval selects an integration step, not another material or
//! another source query. Rebinding duration shares the original observations,
//! receipts and immutable curve arrays while retaining caller-owned history.
//! Changing temperature requires the corresponding thermodynamic family.

use super::TabulatedOverstressJ2;
use crate::plastic::J2State;
use crate::plastic::energy::{
    EnergyError, EnergyIncrement, EnergyTolerances, IncrementEnergy, StoredEnergy,
};
use crate::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use crate::plastic::section::energy::NativeSectionEnergyLaw;
use crate::plastic::tabulated::source::{
    HardeningCurveRequest, ResolvedTabulatedHardening, ResolvedTabulatedJ2,
    TabulatedJ2SourceRequest, TabulatedKinematicHardening, TabulatedSourceError,
    YieldStressMeasure, resolve_tabulated_hardening, tabulated_j2_requirements,
};
use crate::plastic::tabulated::{TabulatedError, TabulatedJ2Response};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};
use std::sync::Arc;

/// Explicit physical rate law. Extra properties never select a different law.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TabulatedOverstressKinetics {
    /// alpha_dot = max(f,0)/eta; only the original positive viscosity is required.
    Linear,
    /// alpha_dot = (s0/eta)*max(f/s0,0)^n; s0>0 and n>=1 are independent data.
    PowerLaw,
}

/// Original hardening/kinematic interpretation and one actual physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TabulatedOverstressSourceRequest {
    pub plasticity: TabulatedJ2SourceRequest,
    pub kinetics: TabulatedOverstressKinetics,
    pub duration_s: f64,
}

/// Complete scalar requirements, composable with density or other source roles.
/// No fluid-viscosity alias, inferred stress scale or fitted exponent is used.
pub fn tabulated_overstress_requirements(
    kinematic: TabulatedKinematicHardening,
    kinetics: TabulatedOverstressKinetics,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = tabulated_j2_requirements(kinematic)?;
    requirements.push(ScalarPropertyRequirement::try_new(
        OVERSTRESS_VISCOSITY_PROPERTY,
        OVERSTRESS_VISCOSITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    if kinetics == TabulatedOverstressKinetics::PowerLaw {
        requirements.push(ScalarPropertyRequirement::try_new(
            OVERSTRESS_STRESS_SCALE_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        requirements.push(ScalarPropertyRequirement::try_new(
            OVERSTRESS_EXPONENT_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    Ok(requirements)
}

/// A native physical step sharing the complete original source and evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTabulatedOverstressJ2 {
    source: Arc<ResolvedTabulatedJ2>,
    kinetics: TabulatedOverstressKinetics,
    law: TabulatedOverstressJ2,
}

impl ResolvedTabulatedOverstressJ2 {
    /// Adopt a complete original bundle, retaining any additional source roles.
    /// The native hardening owner already checks exact card/context agreement.
    pub fn try_from_source(
        source: ResolvedTabulatedJ2,
        kinetics: TabulatedOverstressKinetics,
        duration_s: f64,
    ) -> Result<Self, TabulatedSourceError> {
        let material = source.material();
        let eta = coefficient(
            material,
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
        )?;
        let law = match kinetics {
            TabulatedOverstressKinetics::Linear => {
                TabulatedOverstressJ2::new_linear(source.law().clone(), eta, duration_s)?
            }
            TabulatedOverstressKinetics::PowerLaw => TabulatedOverstressJ2::new_power_law(
                source.law().clone(),
                eta,
                coefficient(material, OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS)?,
                coefficient(material, OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE)?,
                duration_s,
            )?,
        };
        Ok(Self {
            source: Arc::new(source),
            kinetics,
            law,
        })
    }

    /// Select another actual interval without resolving, copying or changing data.
    /// Every trial must still start from the same caller-owned accepted history.
    pub fn for_duration(&self, duration_s: f64) -> Result<Self, TabulatedSourceError> {
        let law = match self.kinetics {
            TabulatedOverstressKinetics::Linear => TabulatedOverstressJ2::new_linear(
                self.law.backbone().clone(),
                self.law.viscosity_pa_s(),
                duration_s,
            )?,
            TabulatedOverstressKinetics::PowerLaw => TabulatedOverstressJ2::new_power_law(
                self.law.backbone().clone(),
                self.law.viscosity_pa_s(),
                self.law
                    .stress_scale_pa()
                    .expect("explicit power-law source selection"),
                self.law.exponent(),
                duration_s,
            )?,
        };
        Ok(Self {
            source: self.source.clone(),
            kinetics: self.kinetics,
            law,
        })
    }

    pub fn source(&self) -> &ResolvedTabulatedJ2 {
        &self.source
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source.material()
    }
    pub fn curve(&self) -> &ResolvedTabulatedHardening {
        self.source.curve()
    }
    pub const fn kinetics(&self) -> TabulatedOverstressKinetics {
        self.kinetics
    }
    pub const fn law(&self) -> &TabulatedOverstressJ2 {
        &self.law
    }
    pub const fn duration_s(&self) -> f64 {
        self.law.duration_s()
    }

    pub fn response(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<TabulatedJ2Response, TabulatedError> {
        self.law.response(strain, state)
    }
    pub fn validate_history(&self, strain: &Voigt, state: &J2State) -> Result<(), EnergyError> {
        self.law.validate_history(strain, state)
    }
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<StoredEnergy, EnergyError> {
        self.law.stored_energy(strain, state)
    }
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.law
            .advance_with_energy(old_strain, old_state, strain, tolerances)
    }
}

impl SmallStrainLaw for ResolvedTabulatedOverstressJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

impl NativeSectionEnergyLaw for ResolvedTabulatedOverstressJ2 {
    fn physical_duration_s(&self) -> Option<f64> {
        Some(self.duration_s())
    }
    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(J2State, Voigt, IncrementEnergy), MaterialError> {
        self.law
            .native_energy_increment(old_strain, old_state, strain, tolerances)
    }
}

/// Resolve all scalar coefficients and the original complete Cauchy hardening
/// curve at one fixed physical context. Caller selection is preserved separately
/// for scalar roles and the full hardening observations.
pub fn resolve_tabulated_overstress_j2(
    card: &MaterialCard,
    fixed_context: &QueryPoint,
    request: &TabulatedOverstressSourceRequest,
    material_selection: MaterialPropertySelection,
    curve_selection: ClaimSelection,
) -> Result<ResolvedTabulatedOverstressJ2, TabulatedSourceError> {
    let curve = resolve_tabulated_hardening(
        card,
        fixed_context,
        &HardeningCurveRequest {
            stress_measure: YieldStressMeasure::Cauchy,
            maximum_points: request.plasticity.maximum_points,
        },
        curve_selection,
    )?;
    let material = resolve_material_state_point(
        card,
        fixed_context,
        &tabulated_overstress_requirements(
            request.plasticity.kinematic_hardening,
            request.kinetics,
        )?,
        material_selection,
    )?;
    let source = ResolvedTabulatedJ2::try_from_material(
        material,
        curve,
        request.plasticity.kinematic_hardening,
        request.plasticity.strain_limit,
    )?;
    ResolvedTabulatedOverstressJ2::try_from_source(source, request.kinetics, request.duration_s)
}

fn coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
) -> Result<f64, TabulatedSourceError> {
    let invalid = || {
        TabulatedSourceError::Material(MaterialError::Parameters {
            what: format!(
                "tabulated flow coefficient {name} requires its exact finite SI quantity at fixed external conditions"
            ),
        })
    };
    let property = material.property(name).ok_or_else(invalid)?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
        || property.angular_frequency_band_rad_s()?.is_some()
    {
        return Err(invalid());
    }
    Ok(property.value_si())
}
