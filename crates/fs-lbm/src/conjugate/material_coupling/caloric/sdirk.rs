//! Physical-stage current feedback on the original caloric SDIRK2 owner.

use super::{
    CaloricElectrothermalConfig, CaloricElectrothermalReport, CoupledStep, HeatStep, admit,
    advance_with_heat, refusal,
};
use crate::conjugate::radiation::{RadiationConfig, Radiators};
use crate::conjugate::tensor_conduction::{advance_radiating_sdirk2_prepared, caloric_stage_clock};
use crate::conjugate::transient::TransientRecord;
use crate::conjugate::{
    ChtError, ElectricSetup, ElectricSolution, ElectricalContact, SolidMaterial,
    TENSOR_SDIRK2_GAMMA, TensorCaloricEndpoint, TensorCaloricRadiationReport,
    TensorCaloricSdirk2Step, TensorCaloricStorage, TensorConductivityLaw, ThermalSetup, Voxel,
    VoxelDomain, advance_caloric_tensor_conduction_sdirk2, poll,
};
use fs_exec::CancelGate;

/// Both original physical stages and the fields supplying their Joule sources.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricElectrothermalSdirk2Step {
    /// Original native stages and integrated heat ledger; Joule heat is included.
    pub thermal: TensorCaloricSdirk2Step,
    /// Applied fields at the corresponding two physical stage times.
    pub applied_electric: [ElectricSolution; 2],
    /// Current refreshed at EACH returned stage's original energy and Kelvin.
    /// These fields supplied no additional heat. Index one is the final state.
    pub endpoint_electric: [ElectricSolution; 2],
    /// Positive tableau-weighted electrical heat and defects; extrema range
    /// over both stages. Current solves equal 2 * (passes + 1).
    pub report: CaloricElectrothermalReport,
}

impl HeatStep<2> for TensorCaloricSdirk2Step {
    fn endpoints(&self) -> [&TensorCaloricEndpoint; 2] {
        self.stages.each_ref()
    }

    fn record(&self) -> &TransientRecord {
        &self.record
    }

    fn into_energies(self) -> [Vec<f64>; 2] {
        self.stages.map(|stage| stage.specific_energy_j_kg)
    }
}

impl CoupledStep<TensorCaloricSdirk2Step, 2> {
    fn sdirk2(self) -> CaloricElectrothermalSdirk2Step {
        CaloricElectrothermalSdirk2Step {
            thermal: self.thermal,
            applied_electric: self.applied_electric,
            endpoint_electric: self.endpoint_electric,
            report: self.report,
        }
    }
}

/// Advance spatial caloric/current feedback at both original SDIRK2 stages.
///
/// Electric drives, electrical contacts, authored heat/boundaries and explicit
/// prescribed energies are sampled at start + gamma*dt and start + dt, with
/// gamma = TENSOR_SDIRK2_GAMMA. The resistivity callback receives that physical
/// stage time, cell, occupancy, original specific energy and temperature. Its
/// source-domain, insulation and original-energy semantics match the BE API.
///
/// Every outer proposal invokes the original TWO-STAGE heat owner with the
/// SAME physical opening energy in both rows. Each row has its own applied
/// electrical source; the native owner alone carries actual first-stage heat
/// into stage two. A discarded proposal advances no inventory or physical time.
/// Both stage energies and both refreshed local Joule fields must converge.
///
/// Electrical energy, thermostat heat and noncancelling local source defects
/// use dt * [(1-gamma)*stage1 + gamma*stage2]. Source combination and native
/// tableau/reduction rounding enter the transfer defect. Each stage's own
/// prescribed mask selects thermostat Joule heat. Native storage and reservoir
/// corrections remain with the original heat ledger, including changing masks.
///
/// Native nonlinear/Krylov/backtracking limits cover both stages per proposal;
/// maximum_passes bounds proposals. Source-law budgets persist across them.
/// Smooth forcing and fixed constraints permit second-order time accuracy when
/// coupling/solver tolerances are resolved. Split forcing jumps and sharp phase
/// corners; no second-order claim is made across an unresolved event. This
/// remains fixed-inventory, stationary DC coupling without electromagnetic
/// storage. Any stage refusal or cancellation returns no partial interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_electrothermal_sdirk2(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    electric: [&ElectricSetup; 2],
    contacts: [&[ElectricalContact]; 2],
    thermal: [&ThermalSetup; 2],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: [&[(usize, f64)]; 2],
    start_time_s: f64,
    time_step_s: f64,
    config: &CaloricElectrothermalConfig,
    resistivity: impl FnMut(f64, usize, Voxel, f64, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<CaloricElectrothermalSdirk2Step, ChtError> {
    poll(gate)?;
    let (_, stage_time_s) = caloric_stage_clock(start_time_s, time_step_s)
        .ok_or_else(|| refusal("SDIRK2 needs strictly ordered finite physical stage times"))?;
    advance_with_heat(
        domain,
        solids,
        electric,
        contacts,
        thermal,
        conductivity_laws,
        storage,
        start_time_s,
        time_step_s,
        stage_time_s,
        [1.0 - TENSOR_SDIRK2_GAMMA, TENSOR_SDIRK2_GAMMA],
        config,
        resistivity,
        |heated| {
            advance_caloric_tensor_conduction_sdirk2(
                domain,
                solids,
                heated,
                conductivity_laws,
                storage,
                prescribed_specific_energy_j_kg,
                start_time_s,
                time_step_s,
                &config.thermal,
                gate,
            )
        },
        gate,
    )
    .map(CoupledStep::sdirk2)
}

/// Advance stage-correct caloric/current feedback with physical gray radiation.
///
/// This uses the same stage source iteration and electrical heat accounting as
/// advance_caloric_tensor_electrothermal_sdirk2. The original radiating owner
/// converges emission/irradiation at BOTH stages, carries physical first-stage
/// radiation once and returns the actual stage reports. Integrated escaped
/// heat is dt * [(1-gamma)*report[0].radiated_w + gamma*report[1].radiated_w];
/// thermal.record.radiated_w remains the final endpoint power.
///
/// Ray geometry is built once for the complete coupled interval. Geometry,
/// emissivities and surroundings stay fixed; participating-media radiation and
/// moving surfaces are outside the original model. Native radiative and thermal
/// work ceilings cover both stages per proposal; maximum_passes bounds the
/// number of proposals. The same authored radiation must not also be supplied
/// through a cell sink or source. Failure returns no stage or interval prefix.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_electrothermal_radiating_sdirk2(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    electric: [&ElectricSetup; 2],
    contacts: [&[ElectricalContact]; 2],
    thermal: [&ThermalSetup; 2],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: [&[(usize, f64)]; 2],
    start_time_s: f64,
    time_step_s: f64,
    config: &CaloricElectrothermalConfig,
    radiation: &RadiationConfig,
    resistivity: impl FnMut(f64, usize, Voxel, f64, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<
    (
        CaloricElectrothermalSdirk2Step,
        [TensorCaloricRadiationReport; 2],
    ),
    ChtError,
> {
    poll(gate)?;
    admit(config, start_time_s, time_step_s, domain.cell_count())?;
    domain.check_materials(solids.len())?;
    let (_, stage_time_s) = caloric_stage_clock(start_time_s, time_step_s).ok_or_else(|| {
        refusal("radiating SDIRK2 needs strictly ordered finite physical stage times")
    })?;
    config
        .maximum_passes
        .checked_mul(radiation.max_iterations)
        .ok_or_else(|| refusal("radiating thermal work bound is unrepresentable"))?;
    let radiators = Radiators::build(domain, solids, radiation, gate)?;
    let mut accepted_radiation = None;
    let step = advance_with_heat(
        domain,
        solids,
        electric,
        contacts,
        thermal,
        conductivity_laws,
        storage,
        start_time_s,
        time_step_s,
        stage_time_s,
        [1.0 - TENSOR_SDIRK2_GAMMA, TENSOR_SDIRK2_GAMMA],
        config,
        resistivity,
        |heated| {
            let (step, reports) = advance_radiating_sdirk2_prepared(
                domain,
                solids,
                heated,
                conductivity_laws,
                storage,
                prescribed_specific_energy_j_kg,
                start_time_s,
                time_step_s,
                &config.thermal,
                radiation,
                Some(&radiators),
                gate,
            )?;
            accepted_radiation = Some(reports);
            Ok(step)
        },
        gate,
    )?;
    let reports = accepted_radiation
        .ok_or_else(|| refusal("accepted radiating stages are missing their reports"))?;
    Ok((step.sdirk2(), reports))
}
