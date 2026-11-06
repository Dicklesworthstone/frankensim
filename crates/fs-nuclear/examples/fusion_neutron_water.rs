//! D-T burn -> original neutron transport -> finite IF97 water receiver.
//!
//! Run: cargo run -p fs-nuclear --features thermal --example fusion_neutron_water
//! Neutron coefficients and deposition fractions below are authored illustrative
//! inputs, not evaluated water cross sections or a blanket-design prediction.
//! Nuclear capture Q/activation and neutron angular transport are not modeled.
//! The water receiver is isobaric and uses reference geometry, not a rigid
//! sealed volume. Its storage is enthalpy, so the combined energy account
//! includes the receiver's isobaric pressure work through that enthalpy change.

#[path = "../../fs-gas/examples/dt_fusion/source.rs"]
mod plasma_source;

use fs_gas::fusion::{DtFusionHeating, DtFusionSolve, advance_dt_fusion};
use fs_material::phase::water::{IsobaricWater, WaterChartControls};
use fs_nuclear::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind, LeakageBoundary,
};
use fs_nuclear::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryEntry, FusionNeutronHistoryLimits,
};
use fs_nuclear::fusion_neutrons::{NeutronKineticBudget, NeutronKineticModel, NeutronKineticState};
use fs_nuclear::irradiation::NeutronGroups;
use fs_nuclear::thermal::{NeutronKineticCellSource, NeutronKineticHeating};
use fs_qty::{Power, Pressure, Temperature, Time};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::fusion::{KELVIN_PER_KEV, NuclearSpeciesV1, SPEED_OF_LIGHT_M_S};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use std::error::Error;
use std::sync::Arc;

fn main() -> Result<(), Box<dyn Error>> {
    let duration_s = 1e-3;
    let (plasma, initial) = plasma_source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let mut controls = DtFusionSolve::new(
        Temperature::new(0.2 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1e4),
        Temperature::new(20.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    controls.relative_tolerance = 1e-8;
    controls.record_neutron_history = true;
    let burn = advance_dt_fusion(
        &initial,
        &plasma,
        DtFusionHeating {
            alpha_electron_fraction: 0.8,
            electron_power: Power::new(0.0),
            heavy_power: Power::new(0.0),
        },
        Time::new(duration_s),
        controls,
        || false,
    )?;
    let emission = burn.neutron_emission()?;

    let water = IsobaricWater::new(
        Pressure::new(101_325.0),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 500.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65_535,
            },
        },
        || false,
    )?;
    let water_initial = water.state_at_temperature(Temperature::new(300.0), || false)?;
    let mass_kg = 0.01;
    let volume_m3 = mass_kg * water_initial.specific_volume_m3_kg();
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg,
        specific_h_j_kg: water_initial.specific_enthalpy_j_kg(),
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.0, 0.0, 0.02],
            volume_m3,
            surface_area_m2: 0.003,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    let thermal_controls = StepConfig {
        max_dt_s: duration_s / 8.0,
        energy_tolerance_j: 1e-6,
    };
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(1), &water)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 4096,
        },
    )?;
    let (mut thermal, _) = thermal.advance_to(0.0, thermal_controls, || Ok(()))?;
    let initial_material_energy_j = thermal.energy_j();

    let mev = 1.602_176_634e-13;
    let energies = [0.1 * mev, 14.0 * mev, 15.0 * mev];
    let neutron_rest_j = NuclearSpeciesV1::Neutron.mass_kg() * SPEED_OF_LIGHT_M_S.powi(2);
    let speeds: Vec<_> = energies
        .iter()
        .map(|&energy| {
            let x = energy / neutron_rest_j;
            SPEED_OF_LIGHT_M_S * (x * (2.0 + x)).sqrt() / (1.0 + x)
        })
        .collect();
    let mut reactions = Vec::new();
    let mut fractions = Vec::new();
    for group in 1..3 {
        reactions.push(DiffusionReaction {
            data_id: format!("illustrative-fast-to-slow-{group}"),
            cell: 0,
            incident_group: group,
            macroscopic_cross_section_m_inv: 0.5,
            kind: DiffusionReactionKind::Scatter { outgoing_group: 0 },
        });
        fractions.push(1.0); // Explicit local recoil-energy approximation.
    }
    for group in 0..3 {
        reactions.push(DiffusionReaction {
            data_id: format!("illustrative-capture-{group}"),
            cell: 0,
            incident_group: group,
            macroscopic_cross_section_m_inv: 0.1,
            kind: DiffusionReactionKind::Capture,
        });
        fractions.push(0.2); // Remaining INCIDENT kinetic energy is outgoing radiation.
    }
    let diffusion = Arc::new(DiffusionModel::new(
        DiffusionData {
            data_id: "illustrative-one-cell-neutron-exposure".into(),
            groups: Arc::new(NeutronGroups::new(
                "illustrative-ascending-groups".into(),
                vec![0.0, mev, 14.5 * mev, 20.0 * mev],
            )?),
            speeds_m_s: speeds,
            cell_volumes_m3: vec![volume_m3],
            delayed_groups: vec![],
            reactions,
            links: vec![],
            leakage: vec![LeakageBoundary {
                cell: 0,
                conductance_m2: vec![0.2 * volume_m3; 3],
            }],
        },
        DiffusionLimits::default(),
    )?);
    let model = Arc::new(NeutronKineticModel::new(
        "illustrative-kinetic-energy-partition-v1".into(),
        diffusion,
        energies.to_vec(),
        fractions,
    )?);
    let history = Arc::new(FusionNeutronHistory::project(
        Arc::clone(&model),
        burn.neutron_history()
            .expect("requested accepted history")
            .iter()
            .map(|&interval| FusionNeutronHistoryEntry { cell: 0, interval })
            .collect(),
        FusionNeutronHistoryLimits::default(),
    )?);
    let neutrons = Arc::new(NeutronKineticState::new(model, vec![0.0; 3])?);
    let sources = NeutronKineticCellSource::bind_with_history(
        neutrons,
        Time::new(0.0),
        vec![(CellId(1), MaterialId(1))],
        Arc::clone(&history),
    )?;
    let mut heating = NeutronKineticHeating::new(
        &thermal,
        sources,
        NeutronKineticBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e6,
            },
            energy_tolerance_j: 1e-6,
        },
    )?;

    println!(
        "Constructed D-T burn with {} accepted neutron source intervals.",
        history.entries().len()
    );
    println!("Authored multigroup diffusion coefficients; no evaluated blanket/activation claim.");
    println!(
        "time_s,water_k,stored_neutron_j,deposited_j,outgoing_reaction_radiation_j,leaked_neutron_j,total_energy_residual_j"
    );
    for target in [0.5 * duration_s, duration_s, 2.0 * duration_s] {
        let (next_thermal, next_heating, _) =
            heating.advance_to(&thermal, target, thermal_controls, || Ok(()))?;
        thermal = next_thermal;
        heating = next_heating;
        let retained = heating.source(CellId(1)).expect("bound neutron source");
        let storage = retained.state().kinetic_energy_j()?;
        let exports = heating.cumulative_released()?;
        let emitted: f64 = history
            .entries()
            .iter()
            .map(|entry| {
                let interval = &entry.interval;
                let duration = interval.emission().duration().value();
                let elapsed = (target - interval.start().value()).clamp(0.0, duration);
                interval.emission().total_kinetic_energy().value() * (elapsed / duration)
            })
            .sum();
        let residual = storage + exports.total_j() - emitted;
        println!(
            "{target:.6e},{:.9},{storage:.9e},{:.9},{:.9},{:.9},{residual:.6e}",
            thermal.temperature_k(CellId(1))?,
            exports.deposited_j,
            exports.radiation_j,
            exports.escaped_j
        );
        assert!(residual.abs() < 1e-5);
        assert!(
            (thermal.energy_j() - initial_material_energy_j - exports.deposited_j).abs() < 1e-6
        );
    }
    let exports = heating.cumulative_released()?;
    let stored_neutron_j = heating
        .source(CellId(1))
        .expect("bound source")
        .state()
        .kinetic_energy_j()?;
    let plasma_change_j = burn.state().heavy_internal_energy().value()
        + burn.state().electron_internal_energy().value()
        - initial.heavy_internal_energy().value()
        - initial.electron_internal_energy().value();
    let whole_residual = plasma_change_j + thermal.energy_j() - initial_material_energy_j
        + stored_neutron_j
        + exports.radiation_j
        + exports.escaped_j
        - burn.balance().nuclear_energy_released_j;
    println!(
        "Plasma internal energy + neutron kinetic energy + receiver enthalpy residual \
         (including isobaric pressure work): {whole_residual:.6e} J"
    );
    assert!(whole_residual.abs() < 1e-5);
    assert!(thermal.temperature_k(CellId(1))? > 300.0);
    Ok(())
}
