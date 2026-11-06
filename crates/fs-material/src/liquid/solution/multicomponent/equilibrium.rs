//! Original liquid-property sources for every accepted equilibrium-candidate phase.
//!
//! The thermochemical owner selects a constrained, seeded local candidate. This
//! adapter retains that candidate and resolves each physical phase at its own
//! finite component inventory. Only extensive mass and phase volumes are
//! aggregated. No effective viscosity, conductivity, heat capacity, morphology,
//! phase-count guarantee or global Gibbs-minimum claim is inferred.

use std::sync::Arc;

use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::multicomponent::MAX_REGULAR_SOLUTION_SPECIES_V1;
use fs_thermochem::solution::multicomponent::equilibrium::MulticomponentEquilibriumCandidateV1;

use super::{
    MulticomponentLiquidComponent, MulticomponentLiquidError, MulticomponentLiquidSolutionLaw,
    ResolvedMulticomponentLiquidCell, Result, SourceOwner, poll, positive, positive_value,
    sum_nonnegative,
};

const RELATIVE_ACCOUNTING_TOLERANCE: f64 = 64.0 * f64::EPSILON;

impl MulticomponentLiquidSolutionLaw<'_> {
    /// Materialize every positive phase of the original thermodynamic candidate.
    ///
    /// The candidate must retain the SAME complete regular-solution model as
    /// this source binding: species, both pair matrices and temperature domain.
    /// Its original finite feed and phase inventories are checked without a
    /// trace floor, rescaling, repartitioning or a second equilibrium search.
    /// Each phase reuses `homogeneous`, including its original source selections
    /// and strict dependent-coordinate accuracy. Source support is required at
    /// the resulting phases, not at an artificial homogeneous feed composition.
    ///
    /// Any source, inventory, representation or cancellation failure discards
    /// the complete unpublished result. The law and its sources remain unchanged.
    pub fn resolve_equilibrium(
        &self,
        equilibrium: MulticomponentEquilibriumCandidateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidEquilibrium> {
        poll(&mut cancelled)?;
        let feed = equilibrium.feed();
        let amounts = equilibrium.component_moles();
        let count = self.species().len();
        if feed.model() != self.thermodynamics()
            || amounts.len() != count
            || feed.mole_fractions().len() != count
            || equilibrium.phases().is_empty()
        {
            return Err(MulticomponentLiquidError::Invalid(
                "equilibrium requires the original complete thermodynamic owner and finite feed",
            ));
        }
        let total_moles = positive_value(
            sum_nonnegative(amounts.iter().copied(), "equilibrium feed amount")?,
            "positive equilibrium feed amount",
        )?;
        if !relative_match(equilibrium.total_moles(), total_moles)
            || amounts.iter().zip(feed.mole_fractions()).any(|(&n, &x)| {
                let fraction = n / total_moles;
                (n > 0.0 && fraction == 0.0) || !relative_match(fraction, x)
            })
        {
            return Err(MulticomponentLiquidError::Arithmetic(
                "equilibrium feed composition and original inventory disagree",
            ));
        }
        let temperature = feed.temperature();
        let mut phases = Vec::with_capacity(equilibrium.phases().len());
        for phase in equilibrium.phases() {
            poll(&mut cancelled)?;
            let accepted = phase.thermodynamics();
            if accepted.model() != feed.model()
                || accepted.temperature() != temperature
                || phase.component_moles().len() != count
                || accepted.mole_fractions().len() != count
            {
                return Err(MulticomponentLiquidError::Invalid(
                    "every equilibrium phase requires the original species, temperature and model",
                ));
            }
            let phase_moles = positive_value(
                sum_nonnegative(
                    phase.component_moles().iter().copied(),
                    "equilibrium phase amount",
                )?,
                "positive equilibrium phase amount",
            )?;
            if !relative_match(phase.total_moles(), phase_moles)
                || !relative_match(phase.mole_fraction_of_feed(), phase_moles / total_moles)
            {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "equilibrium phase amount and original inventory disagree",
                ));
            }
            let cell = self.homogeneous(temperature, phase.component_moles(), &mut cancelled)?;
            if cell
                .mole_fractions()
                .iter()
                .zip(accepted.mole_fractions())
                .any(|(&resolved, &selected)| !relative_match(resolved, selected))
            {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "source-resolved phase composition differs from the accepted equilibrium phase",
                ));
            }
            phases.push(cell);
        }

        let mut component_mass = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for (i, (&amount, component)) in amounts.iter().zip(self.components()).enumerate() {
            poll(&mut cancelled)?;
            let recovered = sum_nonnegative(
                phases.iter().map(|phase| phase.component_moles()[i]),
                "recovered equilibrium component amount",
            )?;
            if !relative_match(recovered, amount) {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "source-resolved phases do not conserve an original component amount",
                ));
            }
            let mass = amount * component.molar_mass.value();
            if !mass.is_finite() || (amount > 0.0 && !positive(mass)) {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "equilibrium component mass",
                ));
            }
            component_mass[i] = mass;
        }
        let total_mass_kg = positive_value(
            sum_nonnegative(component_mass[..count].iter().copied(), "equilibrium mass")?,
            "positive equilibrium mass",
        )?;
        let recovered_mass = sum_nonnegative(
            phases.iter().map(ResolvedMulticomponentLiquidCell::mass_kg),
            "recovered equilibrium mass",
        )?;
        if !relative_match(recovered_mass, total_mass_kg) {
            return Err(MulticomponentLiquidError::Arithmetic(
                "source-resolved phase masses do not conserve original feed mass",
            ));
        }
        let total_volume_m3 = positive_value(
            sum_nonnegative(
                phases
                    .iter()
                    .map(ResolvedMulticomponentLiquidCell::volume_m3),
                "total equilibrium phase volume",
            )?,
            "positive total equilibrium phase volume",
        )?;
        poll(&mut cancelled)?;
        Ok(ResolvedMulticomponentLiquidEquilibrium {
            equilibrium,
            phases: phases.into_boxed_slice(),
            owner: Arc::clone(&self.owner),
            total_mass_kg,
            total_volume_m3,
        })
    }
}

/// The original local candidate with separately source-resolved physical phases.
///
/// Phase order is the original thermochemical order, not density/stratification
/// order. Every phase retains its component amounts and full liquid receipts.
/// The candidate's seeded local-search scope and qualifications are unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMulticomponentLiquidEquilibrium {
    equilibrium: MulticomponentEquilibriumCandidateV1,
    phases: Box<[ResolvedMulticomponentLiquidCell]>,
    owner: Arc<SourceOwner>,
    total_mass_kg: f64,
    total_volume_m3: f64,
}

impl ResolvedMulticomponentLiquidEquilibrium {
    pub fn species(&self) -> &[SpeciesId] {
        self.owner.thermodynamics.species()
    }
    pub fn components(&self) -> &[MulticomponentLiquidComponent] {
        &self.owner.components
    }
    /// Each phase's actual density, viscosity, heat capacity and conductivity.
    pub fn phases(&self) -> &[ResolvedMulticomponentLiquidCell] {
        &self.phases
    }
    /// Original accepted candidate, including its numerical/local-search report.
    pub const fn equilibrium_candidate(&self) -> &MulticomponentEquilibriumCandidateV1 {
        &self.equilibrium
    }
    pub fn temperature(&self) -> Temperature {
        self.equilibrium.feed().temperature()
    }
    /// Original finite feed, in canonical `species()` order, preserved exactly.
    pub fn total_component_moles(&self) -> &[f64] {
        self.equilibrium.component_moles()
    }
    pub fn total_moles(&self) -> f64 {
        self.equilibrium.total_moles()
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass_kg
    }
    /// Sum of actual phase masses divided by their own source-resolved density.
    /// This adds phase volumes; it assumes no pure-component mixing rule.
    pub const fn total_volume_m3(&self) -> f64 {
        self.total_volume_m3
    }
}

fn relative_match(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && actual >= 0.0
        && expected >= 0.0
        && if expected == 0.0 {
            actual == 0.0
        } else {
            actual > 0.0 && (actual - expected).abs() <= RELATIVE_ACCOUNTING_TOLERANCE * expected
        }
}
