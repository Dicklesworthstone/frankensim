//! Sealed gas, retained thermal compaction and original caloric transport.
//!
//! Every feedback proposal solves the same accepted solid history against the
//! original cavity reference sources. Gas pressure, material temperature and
//! permanent volume change meet in the native current-geometry Newton solve.
//! Only the original material's signed physical caloric heat enters its cells.
//!
//! Gas storage and explicit exterior pressure work remain separate accounts.
//! Gas work cancels internally in the combined balance, while each cavity's
//! endpoint-force correction remains a numerical work defect. No cavity energy
//! is deposited a second time as solid heat.
//!
//! These are resolved, homogeneous, isentropic cavities. Gas-wall heat transfer,
//! gas transport and pore-fluid exchange are not inferred. The solid keeps its
//! supplied reference thermal network, source inventory and optional inertia.

/// Explicit gas-wall heat exchange with retained actual-volume gas epochs.
pub mod heat;

use super::{
    FiniteThermalPorousCandidate, FiniteThermalPorousSolid, FiniteThermoplasticConfig,
    FiniteThermoplasticControls, FiniteThermoplasticError, HyperTetSettings, MechanicalProposal,
    feedback, finite,
};
use crate::iqn_ils::driver::CouplingMethod;
use crate::solid_enthalpy::finite_plastic::feedback::MechanicalEnergy;
use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::dd::Dd;
use fs_solid::hyper3::HyperTetSolution;
use fs_solid::hyper3::plastic::porous::FinitePorousPressureEnergy;
use fs_solid::hyper3::plastic::porous::cavity::FinitePorousCavityEnergy;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy;
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, AdiabaticCavityGas, TetCavityGasLoad};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// Native solid equilibrium and gas endpoints at the same accepted temperature.
/// Cavities retain their caller order and independent original inventories.
#[derive(Debug)]
pub struct FiniteThermalPorousCavityEquilibrium<G = ClosedGasParcel> {
    pub solid: HyperTetSolution,
    pub cavities: Vec<TetCavityGasLoad<G>>,
}

/// Original material account plus distinct gas and external-work receipts.
/// The material endpoint defect remains material-only. The complete defect
/// additionally contains held-pressure and every independent cavity correction.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteThermalPorousCavitiesEnergy {
    pub material: FiniteThermoplasticTetEnergy,
    pub pressure: Option<FinitePorousPressureEnergy>,
    pub cavity_energy: Vec<FinitePorousCavityEnergy>,
    /// Compensated sum of the original gas owners' direct energy increments [J].
    pub gas_stored_energy_change_j: f64,
    /// Compensated sum of explicit exterior pressure work on cavity walls [J].
    pub exterior_work_j: f64,
    /// Complete signed external endpoint-work defect, never a heat source [J].
    pub external_work_defect_j: f64,
}

impl MechanicalEnergy for FiniteThermalPorousCavitiesEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.material.internal_energy_change_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.external_work_defect_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.material.dead_load_endpoint_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.material.support_endpoint_work_j
    }
    fn follower_pressure_work_j(&self) -> f64 {
        self.pressure
            .as_ref()
            .map_or(0.0, |pressure| pressure.work_j)
    }
    fn cavity_gas_storage_change_j(&self) -> f64 {
        self.gas_stored_energy_change_j
    }
    fn cavity_exterior_work_j(&self) -> f64 {
        self.exterior_work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        0.0
    }
}

/// Atomic thermal, retained solid, original observer and sealed-gas endpoint.
/// Reuse the original reference while following one isentrope. To begin a new
/// physical heat state at a deformed cavity, explicitly bind a CavityGasEpoch;
/// the heat module retains these physical gas/wall histories as one system.
pub type FiniteThermalPorousCavitiesCandidate<'a, 't, O, G = ClosedGasParcel> =
    FiniteThermalPorousCandidate<
        'a,
        't,
        O,
        FiniteThermalPorousCavitiesEnergy,
        FiniteThermalPorousCavityEquilibrium<G>,
    >;

impl<'a> FiniteThermalPorousSolid<'a> {
    /// Solve compaction and caloric heat feedback with independent sealed gases.
    ///
    /// Optional held pressure is additional to each cavity's explicit exterior
    /// pressure; the caller selects which faces and surroundings each represents.
    /// Every source, geometry, material, equilibrium and coupled energy refusal
    /// retains the original thermal/solid/observer state and gas references.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_cavities<'t, O: DepositionObserver, G: AdiabaticCavityGas>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteThermalPorousCavitiesCandidate<'a, 't, O, G>,
        FiniteThermoplasticError<O::Error>,
    > {
        self.advance_with_cavities_and_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            cavities,
            pressure,
            settings,
            mechanical_controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Use the original bounded Picard or IQN-ILS enthalpy iteration.
    ///
    /// Every proposal and the final actual-temperature verification start from
    /// the same accepted plastic, motion and source histories. Acceptance checks
    /// cellwise heat, temperature, native transport and the complete solid/gas
    /// first law. Acceleration changes no physical state or constitutive law.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_cavities_and_coupling<'t, O: DepositionObserver, G: AdiabaticCavityGas>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteThermalPorousCavitiesCandidate<'a, 't, O, G>,
        FiniteThermoplasticError<O::Error>,
    > {
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
                let coupled = self
                    .solid
                    .advance_with_cavities(
                        temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        cavities,
                        pressure,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Porous)?;
                let step = coupled.step;
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
                let (mut gas_storage, mut exterior_work) = (Dd::ZERO, Dd::ZERO);
                for energy in &coupled.cavity_energy {
                    super::poll(cx)?;
                    gas_storage = gas_storage + Dd::from_f64(energy.gas_stored_energy_change_j);
                    exterior_work = exterior_work + Dd::from_f64(energy.exterior_work_j);
                }
                Ok(MechanicalProposal {
                    state: step.solid,
                    equilibrium: FiniteThermalPorousCavityEquilibrium {
                        solid: step.equilibrium,
                        cavities: coupled.cavities,
                    },
                    energy: FiniteThermalPorousCavitiesEnergy {
                        material: step.energy,
                        pressure: step.pressure,
                        cavity_energy: coupled.cavity_energy,
                        gas_stored_energy_change_j: finite(gas_storage.to_f64())?,
                        exterior_work_j: finite(exterior_work.to_f64())?,
                        external_work_defect_j: coupled.external_work_defect_j,
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
