//! Shared Maxwellian electron-impact coefficients consumed by the existing
//! concentration kinetics owner. This owner still has one common temperature;
//! its thermodynamic reverse law must not be evaluated at a separate Te.

use super::*;
use crate::electron_collision::{
    AVOGADRO_PER_MOL, ELECTRON_MASS_KG, ElectronCollisionKindV1,
    ElectronCrossSectionV1, MaxwellianRateControlsV1, MaxwellianRateWeightingV1,
};

/// Explicitly selected forward-law family, retaining its original source.
#[derive(Clone, Debug, PartialEq)]
pub enum ForwardRateV1 {
    /// Authored modified Arrhenius coefficient.
    Arrhenius(ArrheniusRateV1),
    /// Source-bound binary electron-impact single-ionization coefficient.
    ElectronImpact(Box<ElectronImpactRateV1>),
}
impl From<ArrheniusRateV1> for ForwardRateV1 {
    fn from(value: ArrheniusRateV1) -> Self { Self::Arrhenius(value) }
}
impl From<ElectronImpactRateV1> for ForwardRateV1 {
    fn from(value: ElectronImpactRateV1) -> Self { Self::ElectronImpact(Box::new(value)) }
}
impl ForwardRateV1 {
    pub(super) fn disabled(&self) -> bool {
        matches!(self, Self::Arrhenius(a) if a.reference_rate().value() == 0.0)
    }
    pub(super) fn log_rate(&self, temperature: f64, cancelled: &mut impl FnMut() -> bool) -> Result<(f64, f64), KineticsErrorV1> {
        match self {
            Self::Arrhenius(a) => a.log_rate(temperature),
            Self::ElectronImpact(impact) => {
                let p = impact.cross_section.evaluate_properties(
                    Temperature::new(temperature), MaxwellianRateWeightingV1::Event,
                    impact.controls, cancelled,
                ).map_err(|error| match error {
                    crate::electron_collision::ElectronCollisionErrorV1::Cancelled => KineticsErrorV1::Cancelled,
                    error => KineticsErrorV1::ElectronCollision(error),
                })?;
                // q = NA <sigma v> c_e c_target, with concentrations in mol/m³.
                Ok((p.log_coefficient_m3_s() + fs_math::det::ln(AVOGADRO_PER_MOL),
                    p.logarithmic_temperature_derivative_per_k()))
            }
        }
    }
}

/// Binary electron-impact ionization using the actual Maxwellian integral.
/// The target identity and threshold come from its retained cross section.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronImpactRateV1 {
    electron: SpeciesId,
    cross_section: ElectronCrossSectionV1,
    controls: MaxwellianRateControlsV1,
    maximum_target_thermal_speed_ratio: f64,
}
impl ElectronImpactRateV1 {
    /// Bind an explicitly named electron, ionization source and integration
    /// controls. The stationary-target approximation is admitted subsequently
    /// against source masses: `sqrt(me/m_target) <= maximum_ratio` at the
    /// common temperature. The limit is a caller choice in `(0,1)`, not an
    /// experimental accuracy certificate.
    pub fn new(
        electron: SpeciesId,
        cross_section: ElectronCrossSectionV1,
        controls: MaxwellianRateControlsV1,
        maximum_target_thermal_speed_ratio: f64,
    ) -> Result<Self, KineticsErrorV1> {
        controls.validate().map_err(KineticsErrorV1::ElectronCollision)?;
        if cross_section.kind() != ElectronCollisionKindV1::Ionization
            || electron == *cross_section.target_species()
            || !positive(maximum_target_thermal_speed_ratio)
            || maximum_target_thermal_speed_ratio >= 1.0
        {
            return Err(KineticsErrorV1::Invalid("binary electron-impact source or thermal-speed limit"));
        }
        Ok(Self { electron, cross_section, controls, maximum_target_thermal_speed_ratio })
    }
    /// Explicit electron identity, checked against the source-bank charge/nuclei.
    pub const fn electron_species(&self) -> &SpeciesId { &self.electron }
    /// Original ionization cross-section source, including its energy threshold.
    pub const fn cross_section(&self) -> &ElectronCrossSectionV1 { &self.cross_section }
    /// Original Maxwellian work and error-estimate controls.
    pub const fn controls(&self) -> MaxwellianRateControlsV1 { self.controls }
    /// Maximum source target/electron thermal-speed ratio.
    pub const fn maximum_target_thermal_speed_ratio(&self) -> f64 { self.maximum_target_thermal_speed_ratio }

    pub(super) fn admit(&self, reaction: &ElementaryReactionV1, bank: &IdealGasEquilibriumV1) -> Result<(), KineticsErrorV1> {
        let species = bank.elemental_matrix().species();
        let electron = species.binary_search(&self.electron).map_err(|_| KineticsErrorV1::Invalid("electron absent from source bank"))?;
        let target = species.binary_search(self.cross_section.target_species()).map_err(|_| KineticsErrorV1::Invalid("impact target absent from source bank"))?;
        let masses = bank.models();
        let electron_mass = masses[electron].molar_mass().value() / AVOGADRO_PER_MOL;
        if bank.charge_vector().get(electron) != Some(-1)
            || (0..bank.elemental_matrix().row_count()).any(|i| bank.elemental_matrix().get(i, electron) != Some(0))
            || ((electron_mass - ELECTRON_MASS_KG) / ELECTRON_MASS_KG).abs() > 16.0 * f64::EPSILON
            || bank.charge_vector().get(target) != Some(0)
            || fs_math::det::sqrt(masses[electron].molar_mass().value() / masses[target].molar_mass().value()) > self.maximum_target_thermal_speed_ratio
        {
            return Err(KineticsErrorV1::Invalid("electron/neutral-target identity, physical mass or stationary-target approximation"));
        }
        let on_side = |side: &[(SpeciesId, u32)], id: &SpeciesId, order| side.iter().any(|(s, n)| s == id && *n == order);
        if reaction.pressure_dependence != PressureDependenceV1::Independent
            || reaction.reactants.len() != 2 || reaction.products.len() != 2
            || !on_side(&reaction.reactants, &self.electron, 1)
            || !on_side(&reaction.reactants, self.cross_section.target_species(), 1)
            || !on_side(&reaction.products, &self.electron, 2)
            || !reaction.products.iter().any(|(id, order)| *order == 1 && bank.charge_vector().charge(id) == Some(1))
        {
            return Err(KineticsErrorV1::Invalid("cross-section law requires e + neutral -> 2e + singly charged ion without pressure multiplier"));
        }
        Ok(())
    }
}
