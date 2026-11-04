//! Original liquid sources for homogeneous cells in a diffuse-interface field.
//!
//! A spatial concentration field already represents phase separation. Each
//! local cell therefore retains the original HOMOGENEOUS thermodynamic branch,
//! including a metastable or spinodal composition, when refreshing transport.
//! Calling bulk equilibrium independently in every cell would replace the
//! field's free energy and destroy its interfacial dynamics.
//!
//! The square-gradient model requires common fixed partial molar volumes. Its
//! explicit site concentration must agree with each actual source density and
//! the original molar masses. An incompatible density source refuses; this
//! adapter does not replace that source with an additive-density approximation.

use super::{
    BinaryLiquidError, BinaryLiquidSolutionLaw, ResolvedBinaryLiquidPhase, poll, positive_value,
};
use fs_qty::Temperature;
use fs_thermochem::solution::BinarySolutionStateV1;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;

/// A local homogeneous branch and the unchanged original source-resolved phase.
/// This type deliberately makes no claim of local bulk phase equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryLiquidCell {
    thermodynamics: BinarySolutionStateV1,
    phase: ResolvedBinaryLiquidPhase,
}
impl ResolvedBinaryLiquidCell {
    pub const fn thermodynamics(&self) -> &BinarySolutionStateV1 {
        &self.thermodynamics
    }
    /// Existing fluid/thermal transport consumer with the actual local sources.
    pub const fn phase(&self) -> &ResolvedBinaryLiquidPhase {
        &self.phase
    }
}

impl BinaryLiquidSolutionLaw<'_> {
    /// Refresh a finite HOMOGENEOUS local cell without a coexistence solve.
    /// Both original component amounts and every selected source claim survive.
    /// Phase stability is reported by the original homogeneous law; the caller
    /// supplies the spatial gradient model, history and interfacial mechanics.
    pub fn homogeneous_cell(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryLiquidCell, BinaryLiquidError> {
        poll(&mut cancelled)?;
        if component_moles.iter().any(|n| !n.is_finite() || *n < 0.0) {
            return Err(BinaryLiquidError::Invalid(
                "finite nonnegative component moles required",
            ));
        }
        let total = positive_value(component_moles[0] + component_moles[1], "local cell amount")?;
        let thermodynamics = self.thermodynamics.evaluate_composition(
            temperature,
            [component_moles[0] / total, component_moles[1] / total],
        )?;
        let phase = self.resolve_homogeneous_phase(
            temperature,
            thermodynamics.mole_fractions(),
            component_moles,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(ResolvedBinaryLiquidCell {
            thermodynamics,
            phase,
        })
    }
}

/// The exact gradient thermodynamics and original mixture source binding.
pub struct BinaryLiquidSquareGradientLaw<'a, 'card> {
    liquid: &'a BinaryLiquidSolutionLaw<'card>,
    gradient: BinarySquareGradientV1,
    maximum_relative_volume_mismatch: f64,
}
impl<'a, 'card> BinaryLiquidSquareGradientLaw<'a, 'card> {
    /// Admit an explicit same-pair field owner and a numerical source-volume
    /// allowance in [0,1e-6]. Kappa and equal partial molar volumes require
    /// physical qualification by the caller; agreement at one composition is
    /// not evidence for an entire mixture. Every later cell is checked again.
    pub fn try_new(
        liquid: &'a BinaryLiquidSolutionLaw<'card>,
        gradient: BinarySquareGradientV1,
        maximum_relative_volume_mismatch: f64,
    ) -> Result<Self, BinaryLiquidError> {
        if liquid.thermodynamics() != gradient.thermodynamics()
            || !maximum_relative_volume_mismatch.is_finite()
            || !(0.0..=1e-6).contains(&maximum_relative_volume_mismatch)
        {
            return Err(BinaryLiquidError::Invalid(
                "same original pair and bounded relative source-volume allowance required",
            ));
        }
        Ok(Self {
            liquid,
            gradient,
            maximum_relative_volume_mismatch,
        })
    }
    pub const fn gradient_model(&self) -> &BinarySquareGradientV1 {
        &self.gradient
    }
    pub const fn liquid_law(&self) -> &BinaryLiquidSolutionLaw<'card> {
        self.liquid
    }

    /// Materialize one actual field-cell inventory. No bulk split or property
    /// averaging occurs. Source failure or cancellation publishes no result.
    pub fn resolve_cell(
        &self,
        component_moles: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryLiquidCell, BinaryLiquidError> {
        let cell = self.liquid.homogeneous_cell(
            self.gradient.temperature(),
            component_moles,
            &mut cancelled,
        )?;
        let expected = positive_value(
            (component_moles[0] + component_moles[1]) / self.gradient.molar_concentration_mol_m3(),
            "fixed-site cell volume",
        )?;
        let actual = cell.phase.volume_m3();
        let mismatch = (actual - expected).abs() / actual.max(expected);
        if !mismatch.is_finite() || mismatch > self.maximum_relative_volume_mismatch {
            return Err(BinaryLiquidError::Invalid(
                "original density violates the fixed equal-partial-molar-volume field model",
            ));
        }
        poll(&mut cancelled)?;
        Ok(cell)
    }
}
