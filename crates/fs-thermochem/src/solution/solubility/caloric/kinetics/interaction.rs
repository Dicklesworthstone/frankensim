//! Constitutive interaction data for the original finite-bath thermal solver.
//!
//! The shared solver retains all phase stocks and calls the original kinetic
//! laws itself. An interaction source supplies only its conditional saturation
//! and the additional enthalpy of those explicit states. It cannot replace a
//! rate law, move an inventory or implement another time integrator.

use super::{SolubilityCaloricErrorV1, SolubilityCaloricModelV1};
use crate::solution::solubility::kinetics::{SolubilityKineticStateV1, SolubilityKineticsV1};
use fs_qty::{AmountConcentration, Temperature, Volume};

/// Immutable interaction and caloric source for one original ordered mixture.
///
/// Implementations must derive conditional saturation and enthalpy from the
/// SAME thermodynamic potential, including any temperature-dependent pure-solid
/// calibration. `admit_binding` must establish positive frozen-inventory heat
/// capacity throughout the offered domain, not merely composition convexity at
/// one temperature. Unsupported sources and states must refuse explicitly.
///
/// The solver advances phases at each stage's initial temperature, then solves
/// the heat balance with those final phase stocks held fixed. The existing
/// full/two-half controller measures this first-order phase/heat splitting
/// error. The ideal endpoint-temperature trajectory remains a separate native
/// limit when no interaction source is selected.
///
/// All callbacks must be allocation-free, deterministic and bounded by the
/// declared component/source limits. Cancellation is polled by the original
/// solver between components and scalar temperature trials. A source is bound
/// by its shared Arc identity and travels with thermal states and outcomes.
pub trait SolubilityThermalInteractionV1: core::fmt::Debug + Send + Sync {
    /// Temperature support for BOTH interactions and their caloric derivative.
    fn temperature_domain(&self) -> [Temperature; 2];

    /// Qualify the original phase calorics, ordered pair laws and supported
    /// thermal branch. A common species name alone is not source compatibility.
    fn admit_binding(
        &self,
        caloric: &SolubilityCaloricModelV1,
        kinetics: &[SolubilityKineticsV1],
    ) -> Result<(), SolubilityCaloricErrorV1>;

    /// Conditional interfacial concentration of one component while all other
    /// current dissolved concentrations are frozen. Earlier declared-order
    /// coordinates may already have advanced during this stage.
    fn conditional_saturation(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        phases: &[SolubilityKineticStateV1],
        component: usize,
    ) -> Result<AmountConcentration, SolubilityCaloricErrorV1>;

    /// Direct additional enthalpy increment [J] for explicit initial and final
    /// states at the SAME conserved solvent volume and component inventories.
    /// The shared owner separately accounts for original sensible and ideal
    /// phase-transfer enthalpy exactly once. Implementations should evaluate
    /// differences stably rather than subtract large absolute energy datums.
    fn enthalpy_increment_j(
        &self,
        solvent_volume: Volume,
        from_temperature: Temperature,
        from_phases: &[SolubilityKineticStateV1],
        to_temperature: Temperature,
        to_phases: &[SolubilityKineticStateV1],
    ) -> Result<f64, SolubilityCaloricErrorV1>;
}
