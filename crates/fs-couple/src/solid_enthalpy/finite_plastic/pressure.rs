//! Current-area pressure through retained finite plasticity and heat transport.
//!
//! The original pressure force and follower tangent participate in each native
//! mechanical proposal. Thermal feedback always restarts from accepted history;
//! physical yield/creep/viscous heat and reversible exchange enter the original
//! caloric cells once. Pressure path work, its geometric endpoint correction and
//! numerical kinetic loss remain separate from material heat.
//!
//! The supplied thermal network still owns its reference conductances; no moving
//! heat-transfer surface, exterior-fluid momentum or time-accuracy claim follows.

use super::feedback::{self, MechanicalEnergy, MechanicalProposal};
use super::{
    FiniteThermoplasticCandidate, FiniteThermoplasticConfig, FiniteThermoplasticControls,
    FiniteThermoplasticError, FiniteThermoplasticMaterial, FiniteThermoplasticMaterialFamily,
    FiniteThermoplasticSolid, HyperTetSettings, finite,
};
use crate::iqn_ils::driver::CouplingMethod;
use fs_exec::Cx;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticPressureWork, FiniteThermoplasticTetEnergy,
};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// Native material account and independent original-pressure work. The combined
/// balance uses the complete external-work defect; constitutive defect and
/// pressure geometry correction remain individually available.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteThermoplasticPressureEnergy {
    pub material: FiniteThermoplasticTetEnergy,
    pub pressure: FiniteThermoplasticPressureWork,
}

impl MechanicalEnergy for FiniteThermoplasticPressureEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.material.internal_energy_change_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.pressure.external_work_defect_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.material.dead_load_endpoint_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.material.support_endpoint_work_j
    }
    fn follower_pressure_work_j(&self) -> f64 {
        self.pressure.work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        0.0
    }
}

/// Complete thermal, plastic/momentum and observer candidate with pressure work.
pub type FiniteThermoplasticPressureCandidate<'a, 't, O, M = FiniteThermoplasticMaterial<'a>> =
    FiniteThermoplasticCandidate<'a, 't, O, FiniteThermoplasticPressureEnergy, M>;

impl<'a, M: FiniteThermoplasticMaterialFamily> FiniteThermoplasticSolid<'a, M> {
    /// Hold the original pressure field over the actual coupled interval.
    ///
    /// Uniform pressure retains exact -p delta(V) work; per-face pressure retains
    /// straight-displacement-segment work, and spatial pressure variants retain
    /// their original conservative potentials and wet-face integration. Original
    /// source tables, plastic/creep histories and reference inertia are preserved.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure<'t, 'p, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteThermoplasticPressureCandidate<'a, 't, O, M>,
        FiniteThermoplasticError<O::Error>,
    > {
        self.advance_with_pressure_and_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            surface,
            pressure,
            settings,
            mechanical_controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Pressure-driven finite forming with the original bounded Picard or IQN-ILS
    /// enthalpy iteration. All temperature, cell heat and energy gates still
    /// decide acceptance; acceleration contributes no physical state or heat.
    ///
    /// # Errors
    /// Original source/geometry/pressure refusal, thermal or mechanical solve,
    /// combined energy gate, iteration limit, observer failure or cancellation.
    /// Accepted solid, transport and observer histories are unchanged on failure.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure_and_coupling<'t, 'p, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteThermoplasticPressureCandidate<'a, 't, O, M>,
        FiniteThermoplasticError<O::Error>,
    > {
        let pressure = pressure.into();
        let next = feedback::advance_with_coupling(
            self.binding(),
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            method,
            cx,
            |temperatures, duration_s| {
                let step = self
                    .solid
                    .advance_with_pressure(
                        temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        surface,
                        pressure,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Solid)?;
                let heat_j = step
                    .element_energy
                    .iter()
                    .enumerate()
                    .map(|(i, energy)| {
                        finite(
                            self.solid
                                .reference_volume_m3(i)
                                .expect("retained tetrahedron index")
                                * energy.caloric_heat_j_m3(),
                        )
                    })
                    .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?;
                let pressure = step.pressure.ok_or(FiniteThermoplasticError::Invalid(
                    "native thermoplastic pressure proposal omitted its boundary-work account",
                ))?;
                Ok(MechanicalProposal {
                    state: step.solid,
                    equilibrium: step.equilibrium,
                    energy: FiniteThermoplasticPressureEnergy {
                        material: step.energy,
                        pressure,
                    },
                    inertia: step.inertia,
                    heat_j,
                })
            },
        )?;
        Ok(FiniteThermoplasticCandidate {
            thermal: next.thermal,
            solid: Self {
                solid: next.state,
                cells: next.cells,
                time_s: target_s,
            },
            observer: next.observer,
            report: next.report,
        })
    }
}
