//! Affine temperature evolution of the original bounded Kirchhoff yield table.
//!
//! Plastic-strain knots stay fixed. Each absolute yield ordinate evolves as
//! `Y_i(T)=Y_i(T_ref)+s_i*(T-T_ref)`. Admitting the complete table at both ends
//! of the family's temperature interval proves positivity and nondecreasing
//! resistance throughout it. No temperature or plastic-strain extrapolation,
//! interpolation between incompatible strain coordinates, or heat fraction is
//! inferred. Thermal softening changes resistance without erasing history.
//!
//! The energetic convention is the original integral of `Y(alpha,T)-Y(0,T)`.
//! Its temperature derivative is prepared exactly on the original strain grid.
//! Sharing that immutable preparation keeps frozen endpoints cheap and retains
//! entropy even where the entire endpoint curve becomes a yield plateau.
//! Fixed-history free energy is affine in temperature, so this addition has no
//! mechanical heat-capacity term. The existing elastic/expansion bound remains.

use super::{AffineThermalFiniteJ2, FiniteJ2ThermalError, finite};
use crate::plastic::finite::FiniteJ2Error;
use crate::plastic::tabulated::thermal::AffineTabulatedHardening;
use crate::plastic::tabulated::{TabulatedError, TabulatedHardening};

/// One reference curve and explicitly supplied derivatives on the original
/// Kirchhoff/original-volume basis. Shared table arithmetic carries no stress
/// conversion; this native finite wrapper retains its original interpretation.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2ThermalTabulatedHardening {
    affine: AffineTabulatedHardening,
}

impl FiniteJ2ThermalTabulatedHardening {
    /// Prepare the signed derivative of hardening storage. The receiving finite
    /// family admits its temperature span, initial yield and native history bounds.
    pub fn new(
        reference: TabulatedHardening,
        yield_stress_slopes_pa_k: Vec<f64>,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let affine = AffineTabulatedHardening::new(reference, yield_stress_slopes_pa_k).map_err(
            |error| match error {
                TabulatedError::Parameters(reason) => FiniteJ2ThermalError::Parameters(reason),
                error => native_error(error),
            },
        )?;
        Ok(Self { affine })
    }

    /// Exact immutable curve at the receiving family's reference temperature.
    pub fn reference(&self) -> &TabulatedHardening {
        self.affine.reference()
    }

    /// Explicit derivatives of the original absolute yield ordinates [Pa/K].
    pub fn yield_stress_slopes_pa_k(&self) -> &[f64] {
        self.affine.yield_stress_slopes_pa_k()
    }

    /// Fixed-history derivative of the original Kirchhoff hardening storage.
    /// Its negative contributes to entropy; no material-point allocation occurs.
    pub fn hardening_energy_temperature_derivative_j_m3_k(
        &self,
        alpha: f64,
    ) -> Result<f64, FiniteJ2ThermalError> {
        self.affine
            .hardening_energy_temperature_derivative_j_m3_k(alpha)
            .map_err(native_error)
    }

    /// Prepare an admitted finite endpoint outside constitutive iterations.
    pub(super) fn frozen_table(
        &self,
        temperature_change_k: f64,
    ) -> Result<TabulatedHardening, FiniteJ2ThermalError> {
        self.affine
            .at_temperature_change(temperature_change_k)
            .map_err(native_error)
    }
}

fn native_error(error: TabulatedError) -> FiniteJ2ThermalError {
    match error {
        TabulatedError::Numerical(reason) => FiniteJ2ThermalError::Numerical(reason),
        error => FiniteJ2Error::Tabulated(error).into(),
    }
}

impl AffineThermalFiniteJ2 {
    /// Select one bounded table whose yield ordinates evolve affinely with
    /// temperature. The table replaces the complete isotropic yield radius:
    /// original H and dH/dT must both be zero, and initial yield/value derivative
    /// must match its first knot. The complete supplied strain/temperature
    /// rectangle is admitted. Different selected curves cannot replace history.
    ///
    /// Perzyna and simultaneous Norton/primary creep retain this selection.
    /// Pure Norton creep has no yield hardening and refuses; Voce and tabulated
    /// hardening are exclusive. This does not infer annealing or texture change.
    pub fn with_tabulated_hardening(
        mut self,
        hardening: FiniteJ2ThermalTabulatedHardening,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let p = self.parameters;
        if self.norton_creep
            || self.voce_hardening.is_some()
            || p.reference.hardening_modulus_pa != 0.0
            || p.slopes.hardening_modulus_pa_k != 0.0
            || p.reference.yield_stress_pa != hardening.reference().initial_yield_stress_pa()
            || p.slopes.yield_stress_pa_k != hardening.yield_stress_slopes_pa_k()[0]
            || self.limits.maximum_accumulated_plastic_strain
                > hardening.reference().maximum_plastic_strain()
            || self
                .tabulated_hardening
                .as_ref()
                .is_some_and(|old| old != &hardening)
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "thermal table needs matching initial yield and derivative, zero linear hardening, full history support, a plastic family and no competing or replacement selection",
            ));
        }
        finite(
            p.temperature_range_k[1] * hardening.affine.maximum_absolute_storage_slope_j_m3_k(),
        )?;
        self.tabulated_hardening = Some(hardening);
        // Affine ordinate values and differences achieve extrema at these
        // endpoints. Native admission also checks every return stiffness.
        for temperature in p.temperature_range_k {
            self.at_temperature_for_inspection(temperature)?;
        }
        Ok(self)
    }

    /// Exact original curve and temperature derivatives, shared by endpoints.
    pub fn tabulated_hardening(&self) -> Option<&FiniteJ2ThermalTabulatedHardening> {
        self.tabulated_hardening.as_ref()
    }
}
