//! Finite isothermal mixing in one original binary liquid model.
//!
//! At fixed temperature and pressure, conserved pure-component enthalpy and
//! entropy datums cancel. Only the equilibrium mixing contributions remain.
//! The required heat is therefore delta H_mix, including both coexistence
//! phases. A prescribed-temperature reservoir takes the opposite heat account;
//! its entropy change is -Q/T and total generated entropy is -delta G_mix/T.

use super::{
    BinaryCoexistenceOptionsV1, BinaryLiquidError, BinaryLiquidSolutionLaw,
    ResolvedBinaryLiquidEquilibrium, poll,
};
use std::sync::Arc;

/// One finite, constant-pressure mixing transaction at a common temperature.
///
/// The returned heat duty must be exchanged with the isothermal owner. It is
/// not an extra heat-capacity term, latent heat, or an adiabatic temperature
/// prediction. Bulk equilibrium supplies no mixing time, interfacial energy,
/// stirring work, or kinetic-energy dissipation.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryLiquidMixing {
    equilibrium: ResolvedBinaryLiquidEquilibrium,
    heat_into_liquid_j: f64,
    gibbs_change_j: f64,
    liquid_entropy_change_j_per_k: f64,
    entropy_generated_j_per_k: f64,
}

impl BinaryLiquidMixing {
    /// Final conserved species, equilibrium phases and their original sources.
    pub const fn equilibrium(&self) -> &ResolvedBinaryLiquidEquilibrium {
        &self.equilibrium
    }

    /// Take the final inventory for another mixing, transport or buoyancy step.
    pub fn into_equilibrium(self) -> ResolvedBinaryLiquidEquilibrium {
        self.equilibrium
    }

    /// Heat required INTO the liquid to retain its prescribed temperature [J].
    /// Positive is endothermic mixing; negative releases heat to the reservoir.
    pub const fn heat_into_liquid_j(&self) -> f64 {
        self.heat_into_liquid_j
    }

    /// Final minus initial Gibbs mixing energy [J], nonpositive analytically.
    /// The raw floating-point difference is retained, including near-zero roundoff.
    pub const fn gibbs_change_j(&self) -> f64 {
        self.gibbs_change_j
    }

    /// Entropy change of the liquids alone [J/K]; it need not be positive.
    pub const fn liquid_entropy_change_j_per_k(&self) -> f64 {
        self.liquid_entropy_change_j_per_k
    }

    /// Entropy change of the prescribed-temperature heat reservoir [J/K].
    pub fn reservoir_entropy_change_j_per_k(&self) -> f64 {
        -self.heat_into_liquid_j / self.equilibrium.temperature().value()
    }

    /// Combined liquid/reservoir entropy increase, `-delta G/T` [J/K].
    /// This equals `delta S_liquid-Q/T` to floating-point roundoff. No clamp
    /// turns a small signed roundoff residual into a claimed physical source.
    pub const fn entropy_generated_j_per_k(&self) -> f64 {
        self.entropy_generated_j_per_k
    }
}

impl BinaryLiquidSolutionLaw<'_> {
    /// Combine two finite equilibrium inventories at their SAME temperature.
    ///
    /// Both inputs must retain this exact immutable source binding: pair law,
    /// ordered species and molar masses, selected material claims, original
    /// query context and source roles. Identically reconstructed bindings are
    /// accepted. Merely sharing species names or material labels is insufficient.
    ///
    /// Every species is conserved separately. The original equilibrium owner
    /// determines whether the combined inventory is homogeneous or phase-split,
    /// then refreshes the original transport claims at each phase composition.
    /// Actual phase enthalpies determine heat duty; evaluating only the feed's
    /// homogeneous mixing enthalpy would incorrectly heat a demixed liquid.
    ///
    /// Pure-component caloric datums cancel at common T and fixed pressure.
    /// The source's mixture Cp is therefore neither integrated nor counted as
    /// a second mixing contribution. Different temperatures require a separate
    /// sensible-enthalpy balance and refuse here.
    ///
    /// Work is bounded by one original equilibrium solve. All inputs remain
    /// immutable; cancellation, arithmetic or source failure publishes no result.
    pub fn mix_isothermal(
        &self,
        left: &ResolvedBinaryLiquidEquilibrium,
        right: &ResolvedBinaryLiquidEquilibrium,
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidMixing, BinaryLiquidError> {
        poll(&mut cancelled)?;
        let owns = |parcel: &ResolvedBinaryLiquidEquilibrium| {
            Arc::ptr_eq(&self.owner, &parcel.owner) || self.owner.as_ref() == parcel.owner.as_ref()
        };
        if !owns(left) || !owns(right) {
            return Err(BinaryLiquidError::Invalid(
                "isothermal mixing requires the original thermodynamic and material binding",
            ));
        }
        if left.temperature != right.temperature {
            return Err(BinaryLiquidError::Invalid(
                "isothermal mixing requires equal absolute temperatures",
            ));
        }
        let components = [
            combined_component(left.component_moles[0], right.component_moles[0])?,
            combined_component(left.component_moles[1], right.component_moles[1])?,
        ];
        let equilibrium =
            self.equilibrate(left.temperature, components, options, &mut cancelled)?;
        let moles = [
            left.component_moles.iter().sum::<f64>(),
            right.component_moles.iter().sum::<f64>(),
        ];
        let initial = [
            left.thermodynamic_equilibrium(),
            right.thermodynamic_equilibrium(),
        ];
        let final_state = equilibrium.thermodynamic_equilibrium();
        let heat_into_liquid_j = change(
            moles,
            initial.map(|state| state.mixing_enthalpy_j_per_mol()),
            final_state.mixing_enthalpy_j_per_mol(),
            "mixing heat duty",
        )?;
        let gibbs_change_j = change(
            moles,
            initial.map(|state| state.mixing_gibbs_j_per_mol()),
            final_state.mixing_gibbs_j_per_mol(),
            "mixing Gibbs change",
        )?;
        let liquid_entropy_change_j_per_k = change(
            moles,
            initial.map(|state| state.mixing_entropy_j_per_mol_k()),
            final_state.mixing_entropy_j_per_mol_k(),
            "mixing liquid entropy change",
        )?;
        let entropy_generated_j_per_k = checked(
            -gibbs_change_j / left.temperature.value(),
            "mixing entropy generation",
        )?;
        checked(
            -heat_into_liquid_j / left.temperature.value(),
            "mixing reservoir entropy change",
        )?;
        poll(&mut cancelled)?;
        Ok(BinaryLiquidMixing {
            equilibrium,
            heat_into_liquid_j,
            gibbs_change_j,
            liquid_entropy_change_j_per_k,
            entropy_generated_j_per_k,
        })
    }
}

// Cancel the common final molar value before multiplying by inventory. This
// gives exact zero for identical intensive states and avoids subtracting two
// large extensive energies. Species and molar balance hold to roundoff.
fn change(
    moles: [f64; 2],
    initial: [f64; 2],
    final_value: f64,
    name: &'static str,
) -> Result<f64, BinaryLiquidError> {
    let left = checked(moles[0] * (final_value - initial[0]), name)?;
    let right = checked(moles[1] * (final_value - initial[1]), name)?;
    checked(left + right, name)
}

fn checked(value: f64, name: &'static str) -> Result<f64, BinaryLiquidError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BinaryLiquidError::Arithmetic(name))
    }
}

fn combined_component(left: f64, right: f64) -> Result<f64, BinaryLiquidError> {
    let total = checked(left + right, "combined component inventory")?;
    if left > 0.0 && right > 0.0 && (total == left || total == right) {
        return Err(BinaryLiquidError::Arithmetic(
            "positive incoming component lost in addition",
        ));
    }
    Ok(total)
}
