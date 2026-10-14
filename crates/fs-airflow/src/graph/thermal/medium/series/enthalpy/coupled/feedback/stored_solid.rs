//! Joint solid/coolant storage using the original spatial enthalpy solver.
//!
//! Both histories use the same backward-Euler interval. Every nonlinear trial
//! starts from the immutable accepted solid enthalpy, never a previous trial.
//! The actual FEM Robin flux enters the existing coolant/reservoir calculation;
//! body heat is not also injected into the tank. Only a fully closed step
//! publishes either history. Pump/pipe heat omissions keep their original gates.
//!
//! Solid reference mass/geometry is fixed; reservoir pressure is fixed and its
//! volume may vary. This is not moving-material transport, pipe residence time,
//! phase-dependent mechanics, or a complete mechanical-energy certificate.

use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ElementMaterials, ScalarField, ThermalBoundary, ThermalInterfaces,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
};
use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;

use crate::AirflowError;
use crate::conjugate::SolidRegionState;
use super::super::LiquidSolidResponse;
use super::pressure::parallel::ParallelLiquidNetwork;
use super::pressure::parallel::solid::CoupledParallelLiquidSolution;
use super::pressure::parallel::pump::{PumpHeadCurve, PumpOperatingPoint};
use super::pressure::parallel::pump::solid::PumpCouplingControls;
use super::recirculation::{
    LiquidReservoir, RecirculatingLiquidStep, RecirculationError, ReservoirStepControls,
};

/// Accepted histories on their original storage owners and ONE physical clock.
/// No mutable slice or independent solid clock is exposed. The spatial storage
/// owner retains its original chart and invariant nodal masses across steps.
#[derive(Debug)]
pub struct StoredCoolingState<'a, 's, 'm, 'c> {
    reservoir: LiquidReservoir<'a>,
    storage: &'s EnthalpyBackwardEuler<'m, 'c>,
    solid_specific_enthalpy_j_kg: Vec<f64>,
}

impl<'a, 's, 'm, 'c> StoredCoolingState<'a, 's, 'm, 'c> {
    /// Bind a uniformly charted spatial solid to an existing coolant reservoir.
    /// The two initial histories are explicitly declared at the reservoir time.
    /// Subsequent steps cannot exchange the storage owner or advance just one.
    pub fn try_new(
        cx: &Cx<'_>, reservoir: LiquidReservoir<'a>,
        storage: &'s EnthalpyBackwardEuler<'m, 'c>, solid_specific_enthalpy_j_kg: Vec<f64>,
    ) -> Result<Self, StoredCoolingError> {
        poll(cx)?;
        if solid_specific_enthalpy_j_kg.len() != storage.reference_nodal_masses_kg().len() {
            return Err(StoredCoolingError::Invalid("solid history must match the original spatial storage"));
        }
        for (vertex, &h) in solid_specific_enthalpy_j_kg.iter().enumerate() {
            if vertex % 256 == 0 { poll(cx)?; }
            storage.phase_curve().state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
        }
        poll(cx)?;
        Ok(Self { reservoir, storage, solid_specific_enthalpy_j_kg })
    }

    /// Accepted coolant state and shared physical time.
    pub const fn reservoir(&self) -> &LiquidReservoir<'a> { &self.reservoir }
    /// Accepted mesh-ordered solid specific enthalpy, J/kg.
    pub fn solid_specific_enthalpy_j_kg(&self) -> &[f64] { &self.solid_specific_enthalpy_j_kg }
    /// Original spatial chart/mass owner, never replaced by a thermal trial.
    pub const fn storage(&self) -> &'s EnthalpyBackwardEuler<'m, 'c> { self.storage }
}

/// Original spatial inputs. Only the boundary partition is rebuilt from the
/// current coolant rows. Conductivity may be tensor/temperature dependent and
/// explicitly element assigned; reference storage is the state's original one.
#[derive(Clone, Copy, Debug)]
pub struct StoredSolidProblem<'m> {
    pub mesh: &'m ConductionMesh,
    pub material: &'m ConductivityModel,
    pub element_materials: Option<&'m ElementMaterials>,
    /// Endpoint volumetric body heating, W/m3, owned by the solid only.
    pub source: &'m ScalarField,
    /// Original matching thermal interfaces, not an external energy source.
    pub interfaces: Option<&'m ThermalInterfaces>,
}

/// All existing inner controls plus an independently checked combined account.
#[derive(Clone, Copy, Debug)]
pub struct StoredCoolingControls {
    pub pump: PumpCouplingControls,
    pub reservoir: ReservoirStepControls,
    pub solid: EnthalpyStepConfig,
    /// Absolute allowance for BOTH aggregate interface mismatch and combined
    /// storage minus external input, J. This never relaxes either owner's gates.
    pub joint_energy_tolerance_j: f64,
}

/// Original typed physics errors survive the existing callback boundary.
#[derive(Debug)]
pub enum StoredCoolingError {
    Invalid(&'static str),
    Circuit(RecirculationError),
    Solid(EnthalpyError),
    Boundary(ConductionError),
    EnergyBalance { account: &'static str, residual_j: f64, tolerance_j: f64 },
}
impl core::fmt::Display for StoredCoolingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "stored solid/coolant step: {self:?}")
    }
}
impl std::error::Error for StoredCoolingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Circuit(e) => Some(e), Self::Solid(e) => Some(e),
            Self::Boundary(e) => Some(e), _ => None,
        }
    }
}
impl From<RecirculationError> for StoredCoolingError {
    fn from(e: RecirculationError) -> Self { Self::Circuit(e) }
}
impl From<EnthalpyError> for StoredCoolingError {
    fn from(e: EnthalpyError) -> Self { Self::Solid(e) }
}
impl From<ConductionError> for StoredCoolingError {
    fn from(e: ConductionError) -> Self { Self::Boundary(e) }
}

/// Actual original pump/common-solid result, with the original FEM solution.
pub type StoredCoolingCircuit<'p, 'a> =
    PumpOperatingPoint<'p, CoupledParallelLiquidSolution<'a, EnthalpyStepSolution>>;

#[derive(Debug)]
pub struct StoredCoolingStep<'p, 'a> {
    pub liquid: RecirculatingLiquidStep<StoredCoolingCircuit<'p, 'a>>,
    /// dt * (body generation - prescribed outward Neumann flux), J.
    /// Coolant Robin transfer and internal contacts are NOT external sources.
    pub solid_external_heat_j: f64,
    /// dt * solid Robin outflow minus actual net coolant return energy, J.
    pub interface_energy_residual_j: f64,
    /// Sum of actual solid and reservoir storage minus both external inputs, J.
    pub combined_energy_residual_j: f64,
}

impl<'n, 'a, 'card> ParallelLiquidNetwork<'n, 'a, 'card> {
    /// Advance the actual P1 solid enthalpy field AND the recirculating coolant.
    ///
    /// `boundary` only maps current rows to the original mesh. The numerical
    /// solve is always `EnthalpyBackwardEuler::advance`, using one unchanged
    /// history and exactly `end_time_s - state.reservoir().time_s()`. All Robin
    /// regions must belong to this coolant network in branch-major order. Other
    /// faces may be adiabatic or prescribed Neumann; Dirichlet remains refused.
    /// Boundary callbacks must not write externally or mutate accepted state.
    ///
    /// Existing local nodal, solid-energy, regional watt, pump, source and
    /// reservoir-joule gates remain active. The combined check is additional.
    /// Failed or cancelled calls change neither history; returned fields are
    /// exactly the accepted endpoint, not a re-solve after publication.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stored_solid_with_pump<'p, F>(
        &self, cx: &Cx<'_>, state: &mut StoredCoolingState<'a, '_, '_, '_>,
        end_time_s: f64, reservoir_external_heat_j: f64, initial_walls_k: &[f64],
        pump: &'p PumpHeadCurve, problem: StoredSolidProblem<'_>,
        controls: StoredCoolingControls, mut boundary: F,
    ) -> Result<StoredCoolingStep<'p, 'a>, StoredCoolingError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<ThermalBoundary, ConductionError>,
    {
        poll(cx)?;
        let dt = end_time_s - state.reservoir.time_s();
        let tolerance = controls.joint_energy_tolerance_j;
        if !end_time_s.is_finite() || !dt.is_finite() || dt <= 0.0
            || !tolerance.is_finite() || tolerance <= 0.0
        { return Err(StoredCoolingError::Invalid("finite forward time and positive joint joule tolerance required")); }
        // The existing reservoir API commits on success; operate on a disposable
        // copy until the solid and combined accounts have also passed.
        let mut reservoir = state.reservoir.clone();
        let storage = state.storage;
        let old = &state.solid_specific_enthalpy_j_kg;
        let mut original_error = None;
        let attempted = self.advance_reservoir_with_pump_coupled_iqn(
            cx, &mut reservoir, end_time_s, reservoir_external_heat_j, initial_walls_k,
            pump, &controls.pump, controls.reservoir, |cx, rows| {
                let trial = (|| -> Result<_, StoredCoolingError> {
                    let mapped = boundary(cx, rows)?;
                    let solution = storage.advance(cx, ConductionProblem {
                        mesh: problem.mesh, boundary: &mapped, material: problem.material,
                        element_materials: problem.element_materials, source: problem.source,
                    }, problem.interfaces, old, dt, controls.solid)?;
                    let regions = solution.robin_fluxes.iter()
                        .map(SolidRegionState::from_robin_flux).collect();
                    Ok(LiquidSolidResponse { regions, state: solution })
                })();
                trial.map_err(|error| {
                    original_error = Some(error);
                    AirflowError::InvalidConjugateInput {
                        field: "spatial stored-solid callback refused; original error retained", value_bits: 0,
                    }
                })
            },
        );
        let liquid = match attempted {
            Ok(result) => result,
            Err(error) => return Err(original_error.unwrap_or_else(|| error.into())),
        };
        let solid = &liquid.circuit.state.solid.state;
        let solid_external_heat_j = dt * (solid.source_w - solid.neumann_out_w);
        let interface_energy_residual_j = dt * solid.robin_out_w - liquid.returned_heat_j;
        // Difference coordinates avoid dependence on either chart's arbitrary
        // absolute enthalpy origin. Do not add interface transfer twice.
        let combined_energy_residual_j = (solid.stored_energy_change_j - solid_external_heat_j)
            + (liquid.stored_enthalpy_change_j - liquid.external_heat_j);
        for (account, residual) in [
            ("solid/coolant interface", interface_energy_residual_j),
            ("combined solid/reservoir storage", combined_energy_residual_j),
        ] {
            if !residual.is_finite() || residual.abs() > tolerance {
                return Err(StoredCoolingError::EnergyBalance {
                    account, residual_j: residual, tolerance_j: tolerance,
                });
            }
        }
        if !solid_external_heat_j.is_finite() || solid.specific_enthalpy_j_kg.len() != old.len() {
            return Err(StoredCoolingError::Invalid("unrepresentable solid input or changed history shape"));
        }
        poll(cx)?;
        // No fallible work, allocation or callback after publication begins.
        state.solid_specific_enthalpy_j_kg.copy_from_slice(&solid.specific_enthalpy_j_kg);
        state.reservoir = reservoir;
        Ok(StoredCoolingStep { liquid, solid_external_heat_j,
            interface_energy_residual_j, combined_energy_residual_j })
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), StoredCoolingError> {
    super::checkpoint(cx).map_err(RecirculationError::from).map_err(StoredCoolingError::from)
}
