//! Charge-driven productive reaction amounts with a finite reagent inventory.
//!
//! The caller declares the net reaction, deposited product, electron count,
//! molar mass and current efficiency. The existing Faraday owner computes the
//! product amount; the exact written-equation coefficient sets the extent.
//! The original conserved reaction then admits every reagent against its stock.
//! No stock is clipped and failed calls leave the input unchanged.
//!
//! Unallocated charge has NO implied material pathway. In particular, at less
//! than unit efficiency the modeled oxygen, acid and water changes are only
//! those of the declared productive reaction. They are not a complete cell
//! balance or a prediction of hydrogen evolution, redox shuttles or corrosion.
//! No voltage, energy, transport limit, kinetics or efficiency is inferred.
//! Faraday and copper-cell context: Beukes and Badenhorst (2009),
//! <https://www.saimm.co.za/Journal/v109n06p343.pdf>.

/// Reversible potentials bound to the original equation and electron basis.
pub mod potential;

use super::{
    ConservedReaction, ReactionBalanceError, faraday_yield, nonnegative, product, representable,
};
use crate::{MassAmountBasis, SpeciesId};
use fs_qty::{Amount, Current, Dimensionless, ElectricCharge, MolarMass, Time};

type Result<T> = core::result::Result<T, ReactionBalanceError>;

/// An original conserved equation and an explicit Faraday product association.
#[derive(Debug, Clone)]
pub struct ElectrolysisReaction {
    reaction: ConservedReaction,
    product: SpeciesId,
    molar_mass: MolarMass,
    electrons_per_product: u32,
    product_coefficient: i128,
}

/// Atomic result for the productive pathway, with unused chemistry explicit.
#[derive(Debug, Clone)]
pub struct ElectrolysisAdvance {
    /// Original reaction, product association and electron/mass conventions.
    pub model: ElectrolysisReaction,
    /// Total integrated external charge [C].
    pub passed_charge: ElectricCharge,
    /// Caller-declared fraction assigned to the productive pathway, in [0,1].
    pub current_efficiency: Dimensionless,
    /// Passed charge times current efficiency [C].
    pub productive_charge: ElectricCharge,
    /// Charge outside the productive pathway [C]; its side chemistry is absent.
    pub unallocated_charge: ElectricCharge,
    /// Admitted amount of the integer reaction as written [mol].
    pub extent: Amount,
    /// Product formed by that extent, including its original molar mass [kg/mol].
    pub product_yield: MassAmountBasis,
    /// Canonical inventory after ONLY the declared productive reaction [mol].
    pub modeled_final_inventory: Vec<(SpeciesId, Amount)>,
}

impl ElectrolysisReaction {
    /// Bind one positively formed species to its electron transfer per formula.
    ///
    /// For `2 Cu2+ + 2 H2O -> 2 Cu + O2 + 4 H+`, select Cu with TWO electrons
    /// per deposited atom. The written equation therefore transfers FOUR.
    /// Atom/charge conservation does not infer or authenticate this supplied
    /// electron association, molar mass, current efficiency or electrode route.
    pub fn new(
        reaction: ConservedReaction,
        product: SpeciesId,
        molar_mass: MolarMass,
        electrons_per_product: u32,
    ) -> Result<Self> {
        if electrons_per_product == 0 {
            return Err(ReactionBalanceError::Invalid(
                "electron count per deposited product must be positive",
            ));
        }
        MassAmountBasis::from_amount(Amount::new(0.0), molar_mass)?;
        let product_coefficient = reaction
            .terms()
            .iter()
            .find(|term| term.species == product && term.coefficient > 0)
            .ok_or_else(|| ReactionBalanceError::Species {
                species: product.clone(),
                reason: "Faraday product must be a formed species in the original equation",
            })?
            .coefficient;
        Ok(Self {
            reaction,
            product,
            molar_mass,
            electrons_per_product,
            product_coefficient,
        })
    }

    /// Original exactly conserved net reaction, including its certificate.
    pub const fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original deposited-product species, retaining its phase identity.
    pub const fn product(&self) -> &SpeciesId {
        &self.product
    }
    /// Caller-declared mass per mole of deposited product [kg/mol].
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }
    /// Electron transfer per deposited product formula unit.
    pub const fn electrons_per_product(&self) -> u32 {
        self.electrons_per_product
    }
    /// Electron transfer per integer equation, retaining its coefficient basis.
    pub const fn electrons_per_equation(&self) -> u128 {
        self.product_coefficient as u128 * self.electrons_per_product as u128
    }

    /// Advance from nonnegative integrated charge and explicit efficiency.
    ///
    /// Supply every active species exactly once, including zero products; use
    /// the original reaction's inventory conventions. Excess charge refuses
    /// when its productive extent exceeds any reagent, rather than silently
    /// reducing efficiency. Zero charge/efficiency leaves the inventory intact.
    /// The extent ratio is rounded down one step only if its exact product
    /// would exceed the original Faraday amount. This does not repair stock.
    pub fn advance_charge(
        &self,
        initial: &[(SpeciesId, Amount)],
        charge: ElectricCharge,
        current_efficiency: Dimensionless,
    ) -> Result<ElectrolysisAdvance> {
        let faraday = faraday_yield(
            charge,
            self.molar_mass,
            self.electrons_per_product,
            current_efficiency.value(),
        )?;
        let coefficient = self.product_coefficient as f64;
        let mut extent = faraday.amount().value() / coefficient;
        if coefficient.mul_add(extent, -faraday.amount().value()) > 0.0 {
            extent = extent.next_down();
        }
        representable(
            extent,
            faraday.amount().value() > 0.0,
            "electrolysis extent",
        )?;
        let extent = Amount::new(extent);
        let inventory = self.reaction.advance_inventory(initial, extent)?;
        let deposited = product(extent.value(), coefficient, "electrolysis product")?;
        let product_yield = MassAmountBasis::from_amount(Amount::new(deposited), self.molar_mass)?;
        representable(
            product_yield.mass().value(),
            deposited > 0.0,
            "electrolysis mass",
        )?;
        let productive = product(
            charge.value(),
            current_efficiency.value(),
            "productive electrolysis charge",
        )?;
        let unallocated = product(
            charge.value(),
            1.0 - current_efficiency.value(),
            "unallocated electrolysis charge",
        )?;
        Ok(ElectrolysisAdvance {
            model: self.clone(),
            passed_charge: charge,
            current_efficiency,
            productive_charge: ElectricCharge::new(productive),
            unallocated_charge: ElectricCharge::new(unallocated),
            extent,
            product_yield,
            modeled_final_inventory: inventory,
        })
    }

    /// Integrate a constant nonnegative current [A] over elapsed time [s].
    /// For variable current, callers supply independently integrated charge to
    /// [`Self::advance_charge`]. This helper implies no waveform or controller.
    pub fn advance_current(
        &self,
        initial: &[(SpeciesId, Amount)],
        current: Current,
        elapsed: Time,
        current_efficiency: Dimensionless,
    ) -> Result<ElectrolysisAdvance> {
        nonnegative(current.value(), "current must be finite and nonnegative")?;
        nonnegative(
            elapsed.value(),
            "elapsed time must be finite and nonnegative",
        )?;
        let charge = ElectricCharge::new(product(
            current.value(),
            elapsed.value(),
            "constant-current integrated charge",
        )?);
        self.advance_charge(initial, charge, current_efficiency)
    }
}

#[cfg(test)]
#[path = "electrolysis/tests.rs"]
mod tests;
