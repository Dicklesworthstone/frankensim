//! Follower pressure through the original porous thermomechanical feedback.
//!
//! The same current-area load is evaluated inside every original solid solve.
//! Temperature proposals never advance plastic or motion history. Physical
//! compaction heat and reciprocal exchange enter the native caloric timeline;
//! pressure quadrature correction and numerical kinetic loss never do.
//!
//! Thermal transport still uses its supplied reference conductances. This
//! module does not infer deformed heat-transfer areas or a surrounding fluid.

use super::{
    FiniteThermalPorousCandidate, FiniteThermalPorousSolid, FiniteThermoplasticConfig,
    FiniteThermoplasticControls, FiniteThermoplasticError, HyperTetSettings, MechanicalProposal,
    feedback, finite,
};
use crate::iqn_ils::driver::CouplingMethod;
use crate::solid_enthalpy::finite_plastic::feedback::MechanicalEnergy;
use fs_exec::Cx;
use fs_solid::hyper3::plastic::porous::FinitePorousPressureEnergy;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// Original material storage/heat plus independently measured pressure work.
/// The material endpoint defect remains available separately. The complete
/// external-work defect additionally contains the pressure path correction.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteThermalPorousPressureEnergy {
    pub material: FiniteThermoplasticTetEnergy,
    pub pressure: FinitePorousPressureEnergy,
}

impl MechanicalEnergy for FiniteThermalPorousPressureEnergy {
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

/// Atomic thermal, retained solid and observer proposal with pressure accounts.
pub type FiniteThermalPorousPressureCandidate<'a, 't, O> =
    FiniteThermalPorousCandidate<'a, 't, O, FiniteThermalPorousPressureEnergy>;

impl<'a> FiniteThermalPorousSolid<'a> {
    /// Hold a pressure field on the current surface through native heat feedback.
    ///
    /// Uniform pressure uses exact -p delta(V) work; a per-face field uses the
    /// native straight-displacement-segment work. Spatial pressure variants
    /// retain their original conservative potentials and wet-face integration.
    /// Optional inertia keeps its original reference mass and accepted velocity.
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
    ) -> Result<FiniteThermalPorousPressureCandidate<'a, 't, O>, FiniteThermoplasticError<O::Error>>
    {
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

    /// Pressure-driven compaction with the existing bounded Picard or IQN-ILS
    /// iteration. Acceptance still checks original per-cell heat, temperature
    /// and both coupled energy residuals; acceleration adds no physical state.
    /// Every refusal retains the complete accepted thermal/solid/observer triple.
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
    ) -> Result<FiniteThermalPorousPressureCandidate<'a, 't, O>, FiniteThermoplasticError<O::Error>>
    {
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
                    .map_err(FiniteThermoplasticError::Porous)?;
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
                    "native porous pressure proposal omitted its boundary-work account",
                ))?;
                Ok(MechanicalProposal {
                    state: step.solid,
                    equilibrium: step.equilibrium,
                    energy: FiniteThermalPorousPressureEnergy {
                        material: step.energy,
                        pressure,
                    },
                    inertia: step.inertia,
                    heat_j,
                })
            },
        )?;
        Ok(FiniteThermalPorousCandidate {
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
