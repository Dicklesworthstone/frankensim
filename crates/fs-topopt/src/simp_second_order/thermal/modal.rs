//! Raw-density frequency design on the original thermal mechanical realization.
//! Stiffness and prestress use the original SIMP map; consistent inertia uses
//! the same realization's physical occupancy times its source reference density.
//! A robust goal therefore uses eroded occupancy for inertia, while its separate
//! volume/mass inventory remains dilated. The inventory bound is not a different
//! physical body silently spliced into the vibration pencil.
//!
//! These are exact first derivatives of the retained discrete initial-stress
//! approximation. No modal Hessian, continuum stability certificate, or hidden
//! low-density mass interpolation is claimed.

use super::SimpThermalComplianceGoal;
use crate::control::{SolveControl, SolveWork};
use crate::elasticity::second_order::thermal::modal::{
    ThermalFrequencyError, ThermalFrequencyGoal, ThermalFrequencyOptions,
};
use fs_spectral::generalized::GeneralizedError;

/// One source-bound prestressed spectrum and its complete raw-design gradient.
pub struct SimpThermalFrequencyGoal<'g, 'a> {
    source: &'g SimpThermalComplianceGoal<'a>,
    physical: ThermalFrequencyGoal<'g, 'a>,
    gradient: Vec<f64>,
    work: SolveWork,
}
impl<'g, 'a> SimpThermalFrequencyGoal<'g, 'a> {
    pub fn value(&self) -> f64 {
        self.physical.value()
    }
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }
    pub fn source(&self) -> &'g SimpThermalComplianceGoal<'a> {
        self.source
    }
    pub fn physical(&self) -> &ThermalFrequencyGoal<'g, 'a> {
        &self.physical
    }
    /// Cumulative caller work after both modal/adjoint work and raw pullback.
    pub fn work(&self) -> SolveWork {
        self.work
    }
}

impl<'a> SimpThermalComplianceGoal<'a> {
    /// Reuse the original thermal equilibria and design map for frequency design.
    /// One selected case supplies physical prestress independently of its
    /// compliance weight. The native aggregate uses one equilibrium adjoint;
    /// the complete modulus/inertia pullback adds ONE transpose-filter solve.
    ///
    /// Consistent inertia uses source density times actual projected occupancy,
    /// with no artificial mass floor or rho^6 void-mode suppression. All the
    /// native positivity, source-domain, residual and cancellation gates apply.
    pub fn prestressed_frequency_goal<'g>(
        &'g self,
        case_index: usize,
        options: ThermalFrequencyOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<SimpThermalFrequencyGoal<'g, 'a>, ThermalFrequencyError> {
        let physical = self.physical.prestressed_frequency_goal(
            case_index,
            &self.map.projected,
            options,
            control,
        )?;
        let gradient = self
            .map
            .pullback_stiffness_and_mass_gradient(
                physical.modulus_gradient(),
                physical.mass_fraction_gradient(),
                control,
            )
            .map_err(GeneralizedError::Operator)?;
        control
            .checkpoint("simp-thermal-frequency-publish")
            .map_err(|error| GeneralizedError::Operator(error.into()))?;
        Ok(SimpThermalFrequencyGoal {
            source: self,
            physical,
            gradient,
            work: control.work(),
        })
    }
}
