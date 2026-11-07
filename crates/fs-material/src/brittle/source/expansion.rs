//! Source-bound isotropic thermal eigenstrain for the original brittle owner.
//!
//! The shared material source integrates the selected instantaneous alpha(T)
//! exactly over its admitted cells. Total tensor strain is split as
//! epsilon_mechanical = epsilon_total - integral(alpha dT) I. Only the mechanical
//! strain enters the original crack-band/cohesive law; no damage law is added.
//! Tensor-Voigt shear entries and their work weights remain unchanged.
//!
//! Coefficient transfer and eigenstrain evolution are distinct operations.
//! A solid consumer first transfers coefficients at held mechanical strain,
//! then solves the new eigenstrain and displacement together using the borrowed
//! total-strain law. Free expansion must not be preceded by a fictitious fully
//! restrained trial. The optional held-total transfer below instead represents
//! an explicitly restrained material experiment. Both orderings retain source
//! claims, reference temperature and irreversible native history.
//!
//! Coefficient work, eigenstrain work and fracture loss are separate signed
//! energy ports. None is inferred to be heat or entropy. Endpoint coefficient
//! changes are ordered before the straight mechanical-strain path; callers
//! resolve intermediate source changes, reversals and fracture events.

/// Same-source reference density and native path-mean thermal force response.
pub mod dynamics;

use super::{
    BrittleSourceError, ResolvedBrittleCrackBand, brittle_crack_band_requirements, coefficient,
};
use crate::brittle::temperature::BrittleTemperatureStep;
use crate::brittle::{
    BrittleCrackBandSettings, BrittleError, BrittleState, BrittleStep, finite, positive,
};
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    MaterialStatePointError, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_matdb::{MaterialCard, QueryPoint};

/// Complete original fracture/elastic requirements plus instantaneous alpha(T).
/// The expansion coefficient may be negative or zero and must have SI units 1/K.
pub fn thermal_brittle_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = brittle_crack_band_requirements()?;
    requirements.push(ScalarPropertyRequirement::try_new(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        ScalarAdmissibility::Finite,
    )?);
    Ok(requirements)
}

/// One original brittle source and its selected, reference-bound free expansion.
/// Histories remain the native mechanical-strain [BrittleState]; the reference
/// is immutable even when the current temperature changes.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalBrittle {
    band: ResolvedBrittleCrackBand,
    expansion: IntegratedIsotropicThermalExpansion,
}

/// Resolve all fracture/elastic/expansion claims together at the reference and
/// keep those exact claims through the current temperature. Only canonical T
/// changes; other coordinates and their meanings must remain identical.
/// Source gaps along the expansion path refuse even when both endpoints exist.
pub fn resolve_thermal_brittle(
    card: &MaterialCard,
    reference_query: &QueryPoint,
    current_query: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    resolve_with_requirements(
        card,
        reference_query,
        current_query,
        selection,
        settings,
        &thermal_brittle_requirements()?,
    )
}

/// One source-selection, reference-binding and integrated-expansion owner for
/// static and inertial consumers. Extra requirements remain in the complete
/// source bundle and are refreshed with their original selected claims.
fn resolve_with_requirements(
    card: &MaterialCard,
    reference_query: &QueryPoint,
    current_query: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
    requirements: &[ScalarPropertyRequirement],
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    if reference_query.axis_quantities() != current_query.axis_quantities() {
        return Err(BrittleSourceError::Invalid(
            "thermal expansion must retain the reference coordinate conventions",
        ));
    }
    let material = resolve_material_state_point(card, reference_query, requirements, selection)?;
    coefficient(
        &material,
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
    )?;
    let selected = material
        .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
        .ok_or(BrittleSourceError::Invalid(
            "missing thermal-expansion claim",
        ))?
        .answer()
        .receipt
        .selected;
    let expansion = integrate_isotropic_thermal_expansion(
        card,
        reference_query,
        current_query,
        MaterialPropertySelection::PinnedByProperty(vec![(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
            selected,
        )]),
    )?;
    let reference = ResolvedBrittleCrackBand::try_from_material(material, settings)?;
    let temperature =
        expansion
            .current()
            .resolved()
            .state_coordinate("T")
            .ok_or(BrittleSourceError::Invalid(
                "missing current absolute temperature",
            ))?;
    // This native refresh retains the exact prepared plane/modal projection.
    let band = reference.at_temperature(card, "T", temperature)?;
    let candidate = ResolvedThermalBrittle { band, expansion };
    candidate.check_binding()?;
    Ok(candidate)
}

impl ResolvedThermalBrittle {
    /// Complete source bundle and original prepared crack-band owner.
    #[must_use]
    pub const fn band(&self) -> &ResolvedBrittleCrackBand {
        &self.band
    }

    /// Prepare a discovered plane and geometric width while retaining this
    /// exact thermal source, expansion integral and immutable reference.
    ///
    /// The receiver may change only crack-plane geometry. The complete original
    /// material bundle (including extra multiphysics receipts), coefficients and
    /// strain/work controls must agree. This is source preparation, not history
    /// transfer: any elastic handoff must separately check stress and storage,
    /// and an existing damaged state cannot be attached to the new owner.
    pub fn with_resolved_band(
        &self,
        band: ResolvedBrittleCrackBand,
    ) -> Result<Self, BrittleSourceError> {
        if band.material() != self.band.material()
            || band.law().elastic() != self.band.law().elastic()
            || band.law().parameters() != self.band.law().parameters()
        {
            return Err(BrittleSourceError::Invalid(
                "thermal crack-plane preparation must retain the complete original source",
            ));
        }
        let original = self.band.law().settings();
        let target = band.law().settings();
        if target.maximum_strain != original.maximum_strain
            || target.energy_tolerance_j_m3 != original.energy_tolerance_j_m3
        {
            return Err(BrittleSourceError::Invalid(
                "thermal crack-plane preparation must retain strain and work controls",
            ));
        }
        let candidate = Self {
            band,
            expansion: self.expansion.clone(),
        };
        candidate.check_binding()?;
        Ok(candidate)
    }

    /// Exact selected expansion curve, reference/current receipts and integral.
    #[must_use]
    pub const fn expansion(&self) -> &IntegratedIsotropicThermalExpansion {
        &self.expansion
    }

    /// Signed integral of the instantaneous coefficient from the fixed reference.
    #[must_use]
    pub fn free_linear_strain(&self) -> f64 {
        self.expansion.free_linear_strain()
    }

    /// Convert TOTAL tensor strain to the native mechanical strain. Both strain
    /// norms remain inside the original small-strain domain; shear is unchanged.
    pub fn mechanical_strain(&self, total: &Voigt) -> Result<Voigt, BrittleSourceError> {
        self.band.law().check_strain(total)?;
        let free = self.free_linear_strain();
        let mechanical = core::array::from_fn(|i| total[i] - if i < 3 { free } else { 0.0 });
        self.band.law().check_strain(&mechanical)?;
        Ok(mechanical)
    }

    /// Recover total strain from an accepted native mechanical history.
    /// Coefficient-only preparation is an intermediate state: its retained
    /// mechanical strain must still be evolved by the coupled displacement solve.
    pub fn total_strain(&self, state: &BrittleState<'_>) -> Result<Voigt, BrittleSourceError> {
        if !core::ptr::eq(state.owner, self.band.law()) {
            return Err(BrittleError::StateMismatch.into());
        }
        self.band.law().check_strain(state.strain())?;
        let free = self.free_linear_strain();
        let total = core::array::from_fn(|i| state.strain()[i] + if i < 3 { free } else { 0.0 });
        self.band.law().check_strain(&total)?;
        Ok(total)
    }

    /// Disposable original material trial at a TOTAL tensor-strain endpoint.
    pub fn trial_total<'a>(
        &'a self,
        total: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<BrittleStep<'a>, BrittleSourceError> {
        let mechanical = self.mechanical_strain(total)?;
        Ok(self.band.law().trial(&mechanical, state)?)
    }

    /// Borrowed assembly law; its state is the native mechanical history.
    /// The virgin state has zero mechanical strain and stress-free total strain
    /// free_linear_strain() * I, which need not equal zero displacement.
    #[must_use]
    pub fn small_strain_law(&self) -> ThermalBrittleSmallStrain<'_> {
        ThermalBrittleSmallStrain { owner: self }
    }

    /// Refresh all original claims and reintegrate the SAME expansion claim from
    /// the SAME reference. This checks the complete reference-to-current path.
    pub fn at_temperature(
        &self,
        card: &MaterialCard,
        temperature_k: f64,
    ) -> Result<Self, BrittleSourceError> {
        let band = self.band.at_temperature(card, "T", temperature_k)?;
        let expansion = integrate_isotropic_thermal_expansion(
            card,
            &self.expansion.reference().resolved().query(),
            &band.material().query(),
            MaterialPropertySelection::PinnedByProperty(vec![(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
                self.expansion.selected_claim(),
            )]),
        )?;
        let next = Self { band, expansion };
        next.check_binding()?;
        self.require_temperature_successor(&next)?;
        Ok(next)
    }

    /// Require the original card, selected claims, reference, held context,
    /// coordinate conventions, plane, width and native projection.
    pub fn require_temperature_successor(&self, next: &Self) -> Result<(), BrittleSourceError> {
        self.band.require_temperature_successor(&next.band, "T")?;
        if self.expansion.reference().resolved().identity()
            != next.expansion.reference().resolved().identity()
            || self.expansion.selected_claim() != next.expansion.selected_claim()
        {
            return Err(BrittleSourceError::Invalid(
                "thermal successor changed the original expansion reference or selected claim",
            ));
        }
        Ok(())
    }

    /// Prepare coefficient changes at held MECHANICAL strain, preserving native
    /// damage and accumulated fracture work. The solid consumer then solves
    /// displacement and thermal eigenstrain together using the next total-strain
    /// law. This preparation alone is not a held-total-strain thermal experiment.
    pub fn transfer_coefficients<'b>(
        &self,
        state: &BrittleState<'_>,
        next: &'b Self,
        energy_tolerance_j_m3: f64,
    ) -> Result<BrittleTemperatureStep<'b>, BrittleSourceError> {
        self.require_temperature_successor(next)?;
        self.total_strain(state)?;
        self.band
            .transfer_temperature(state, &next.band, "T", energy_tolerance_j_m3)
    }

    /// Explicitly restrained thermal experiment at fixed TOTAL strain.
    ///
    /// First transfer coefficients at held mechanical strain/damage, allowing
    /// native irreversible envelope relaxation. Then move mechanical strain by
    /// -(new_free - old_free) * I under the new coefficients. The original
    /// analytic mean stress supplies the eigenstrain work. This ordering is not
    /// the simultaneous free-expansion equilibrium operation of a solid solver.
    /// Every work balance is checked before the new history is returned.
    pub fn transfer_temperature<'b>(
        &self,
        state: &BrittleState<'_>,
        next: &'b Self,
        energy_tolerance_j_m3: f64,
    ) -> Result<BrittleThermalStep<'b>, BrittleSourceError> {
        positive(energy_tolerance_j_m3)?;
        let coefficient = self.transfer_coefficients(state, next, energy_tolerance_j_m3)?;
        let delta_free = finite(next.free_linear_strain() - self.free_linear_strain())?;
        let mechanical =
            core::array::from_fn(|i| state.strain()[i] - if i < 3 { delta_free } else { 0.0 });
        let path = next
            .band
            .law()
            .trial_with_mean_stress(&mechanical, &coefficient.state)?;
        let eigenstrain_work_j_m3 = finite(
            -delta_free
                * (path.mean_stress_pa[0] + path.mean_stress_pa[1] + path.mean_stress_pa[2]),
        )?;
        let eigenstrain_work_defect_j_m3 = finite(eigenstrain_work_j_m3 - path.path_work_j_m3)?;
        check_work(eigenstrain_work_defect_j_m3, energy_tolerance_j_m3)?;
        let damage_dissipated_j_m3 =
            finite(coefficient.damage_dissipated_j_m3 + path.endpoint.energy.dissipated)?;
        let stored_change_j_m3 =
            finite(path.endpoint.state.stored_energy_j_m3() - state.stored_energy_j_m3())?;
        let defect_j_m3 = finite(
            coefficient.parameter_work_j_m3 + eigenstrain_work_j_m3
                - stored_change_j_m3
                - damage_dissipated_j_m3,
        )?;
        check_work(defect_j_m3, energy_tolerance_j_m3)?;
        if path.endpoint.state.damage() < state.damage()
            || path.endpoint.state.dissipated_j_m3() < state.dissipated_j_m3()
            || damage_dissipated_j_m3 < 0.0
        {
            return Err(BrittleSourceError::Invalid(
                "thermal fracture history cannot heal",
            ));
        }
        Ok(BrittleThermalStep {
            state: path.endpoint.state,
            stress_pa: path.endpoint.stress_pa,
            mean_stress_pa: path.mean_stress_pa,
            coefficient_parameter_work_j_m3: coefficient.parameter_work_j_m3,
            eigenstrain_work_j_m3,
            coefficient_damage_dissipated_j_m3: coefficient.damage_dissipated_j_m3,
            eigenstrain_damage_dissipated_j_m3: path.endpoint.energy.dissipated,
            damage_dissipated_j_m3,
            stored_change_j_m3,
            defect_j_m3,
            eigenstrain_work_defect_j_m3,
        })
    }

    fn check_binding(&self) -> Result<(), BrittleSourceError> {
        let material = self.band.material();
        let current = self.expansion.current().resolved();
        coefficient(
            material,
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
        )?;
        if material.card_identity() != current.card_identity()
            || material.query_point() != current.query_point()
            || material.axis_quantities() != current.axis_quantities()
            || material
                .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
                .is_none_or(|property| {
                    property.answer().receipt.selected != self.expansion.selected_claim()
                })
        {
            return Err(BrittleSourceError::Invalid(
                "brittle coefficients and expansion must share the exact source state",
            ));
        }
        let free = self.free_linear_strain();
        self.band
            .law()
            .check_strain(&[free, free, free, 0.0, 0.0, 0.0])?;
        Ok(())
    }
}

/// Checked restrained thermal update. All energies are densities in J/m3.
/// The mean stress applies only to the eigenstrain stage under new coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleThermalStep<'a> {
    pub state: BrittleState<'a>,
    pub stress_pa: Voigt,
    pub mean_stress_pa: Voigt,
    pub coefficient_parameter_work_j_m3: f64,
    pub eigenstrain_work_j_m3: f64,
    pub coefficient_damage_dissipated_j_m3: f64,
    pub eigenstrain_damage_dissipated_j_m3: f64,
    pub damage_dissipated_j_m3: f64,
    pub stored_change_j_m3: f64,
    pub defect_j_m3: f64,
    /// Physical eigenstrain port minus the native mechanical-path work.
    pub eigenstrain_work_defect_j_m3: f64,
}

/// Total-strain view of the native mechanical history, with no new state owner.
#[derive(Clone, Copy, Debug)]
pub struct ThermalBrittleSmallStrain<'a> {
    owner: &'a ResolvedThermalBrittle,
}

impl<'a> SmallStrainLaw for ThermalBrittleSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State {
        self.owner.band.law().initial_state()
    }
    fn stress(&self, total: &Voigt, state: &Self::State) -> Voigt {
        self.owner
            .trial_total(total, state)
            .map(|step| step.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, total: &Voigt, state: &Self::State) -> Tangent6 {
        self.owner
            .trial_total(total, state)
            .map(|step| step.branch_tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, total: &Voigt, state: &Self::State) -> Self::State {
        self.owner
            .trial_total(total, state)
            .map(|step| step.state)
            .unwrap_or(BrittleState {
                owner: self.owner.band.law(),
                strain: [f64::NAN; 6],
                point: self.owner.band.law().cohesive.bonded_point(),
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.band.small_strain_law().admissibility()
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.brittle-source-isotropic-thermal-expansion", "0.1.0", Ambition::Solid,
            vec![
                "source-integrated isotropic eigenstrain at one immutable temperature reference".into(),
                "total and mechanical tensor strains satisfy the native small-strain bound".into(),
                "native fixed-plane crack-band fracture, irreversible history and tensor-Voigt shear".into(),
                "initial native state is stress-free at total strain equal to thermal eigenstrain".into(),
                "explicit coefficient-first preparation precedes the coupled displacement solve".into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude", 0.0, self.owner.band.law().settings().maximum_strain,
            ),
            vec![
                "no evolving crack plane, crushing, anisotropic expansion or finite thermal stretch".into(),
                "no heat/entropy inference or unresolved thermal-cycle integration".into(),
                "source support and fracture-energy matching are not material qualification".into(),
            ],
            f64::INFINITY,
        )
    }
}

fn check_work(defect: f64, tolerance: f64) -> Result<(), BrittleSourceError> {
    if defect.abs() > tolerance {
        return Err(BrittleError::EnergyDefect(defect).into());
    }
    Ok(())
}
