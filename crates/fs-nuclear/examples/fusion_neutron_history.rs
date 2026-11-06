//! Two evolving D-T sources -> accepted histories -> spatial neutron transport.
//! Run: cargo run -p fs-nuclear --example fusion_neutron_history
//!
//! Neutron coefficients/heat fractions are explicitly illustrative, not evaluated
//! blanket data. This two-fast-group model includes capture, reciprocal diffusion
//! and leakage; 20% of captured INCIDENT kinetic energy is local heat, 80% radiation.
//! It supplies no moderation, capture Q, activation, or photon/material heat solve.

#[path = "../../fs-gas/examples/dt_fusion/source.rs"]
mod plasma_source;

use fs_gas::fusion::{
    DtFusionHeating, DtFusionSolve, DtFusionTransport, FusionEnergyConfinement,
    advance_dt_fusion_with_transport,
};
use fs_nuclear::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionLink, DiffusionModel,
    DiffusionReaction, DiffusionReactionKind, LeakageBoundary,
};
use fs_nuclear::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryEntry, FusionNeutronHistoryLimits,
    NeutronKineticHistoryBudget,
};
use fs_nuclear::fusion_neutrons::{NeutronKineticBudget, NeutronKineticModel, NeutronKineticState};
use fs_nuclear::irradiation::NeutronGroups;
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::fusion::{KELVIN_PER_KEV, NuclearSpeciesV1, SPEED_OF_LIGHT_M_S};
use std::{error::Error, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
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
    controls.record_neutron_history = true;
    controls.relative_tolerance = 1e-5;
    controls.absolute_amount_tolerance_mol = 1e-20;
    controls.initial_step_s = 1e-3;
    controls.maximum_step_s = 2e-3;
    let mut entries = Vec::new();
    let mut emitted = [0.0; 2];
    let mut plasma_remainder_j = 0.0;
    println!("Two 1-litre D-T plasmas: nD=nT=5e19/m3; prompt alpha heat is 80% electrons.");
    println!("Prescribed thermal losses cool each burn on its own accepted clock.");
    let cases = [(0.012, 10.0, 0.008), (0.020, 14.0, 0.012)];
    for (cell, (duration, ion_kev, tau)) in cases.into_iter().enumerate() {
        let (plasma, initial) =
            plasma_source::plasma([5e19, 5e19, 0.0], 0.001, [ion_kev, 0.5 * ion_kev]);
        let thermal = DtFusionTransport {
            energy_confinement: Some(FusionEnergyConfinement::new(
                Time::new(tau),
                Time::new(0.75 * tau),
                Temperature::new(2.0 * KELVIN_PER_KEV),
                Temperature::new(KELVIN_PER_KEV),
            )?),
            ..DtFusionTransport::default()
        };
        let burn = advance_dt_fusion_with_transport(
            &initial,
            &plasma,
            DtFusionHeating {
                alpha_electron_fraction: 0.8,
                electron_power: Power::new(0.0),
                heavy_power: Power::new(0.0),
            },
            &thermal,
            Time::new(duration),
            controls,
            || false,
        )?;
        let source = burn.neutron_emission()?;
        let b = burn.balance();
        let delta = burn.state().heavy_internal_energy().value()
            - initial.heavy_internal_energy().value()
            + burn.state().electron_internal_energy().value()
            - initial.electron_internal_energy().value();
        let other = b.heavy_transport_energy_out_j + b.electron_transport_energy_out_j;
        let residual =
            delta + other + source.total_kinetic_energy().value() - b.nuclear_energy_released_j;
        assert!(residual.abs() < 1e-7);
        plasma_remainder_j += delta + other - b.nuclear_energy_released_j;
        emitted[0] += source.total_particle_count();
        emitted[1] += source.total_kinetic_energy().value();
        let intervals = burn.neutron_history().expect("recording enabled");
        let rate = |i: usize| {
            intervals[i].emission().total_particle_count()
                / intervals[i].emission().duration().value()
        };
        let (first, last, count) = (rate(0), rate(intervals.len() - 1), intervals.len());
        println!("cell {cell}: {count} intervals, {first:.6e} -> {last:.6e} n/s");
        println!("  Independent plasma energy residual: {residual:.3e} J");
        for &interval in intervals {
            entries.push(FusionNeutronHistoryEntry { cell, interval });
        }
    }

    let mev = 1.602_176_634e-13;
    let energies = [14.0 * mev, 15.0 * mev];
    let rest = NuclearSpeciesV1::Neutron.mass_kg() * SPEED_OF_LIGHT_M_S.powi(2);
    let speeds = energies.map(|energy| {
        let x = energy / rest;
        SPEED_OF_LIGHT_M_S * (x * (2.0 + x)).sqrt() / (1.0 + x)
    });
    let groups = Arc::new(NeutronGroups::new(
        "illustrative-DT-groups".into(),
        vec![0.0, 14.5 * mev, 20.0 * mev],
    )?);
    let mut reactions = Vec::new();
    for cell in 0..2 {
        for group in 0..2 {
            reactions.push(DiffusionReaction {
                data_id: format!("illustrative-capture-{cell}-{group}"),
                cell,
                incident_group: group,
                macroscopic_cross_section_m_inv: 0.1,
                kind: DiffusionReactionKind::Capture,
            });
        }
    }
    let diffusion = Arc::new(DiffusionModel::new(
        DiffusionData {
            data_id: "illustrative-two-cell-diffusion".into(),
            groups: Arc::clone(&groups),
            speeds_m_s: speeds.to_vec(),
            cell_volumes_m3: vec![0.5, 0.75],
            delayed_groups: vec![],
            reactions,
            links: vec![DiffusionLink {
                cell_a: 0,
                cell_b: 1,
                conductance_m2: vec![0.02; 2],
            }],
            leakage: (0..2)
                .map(|cell| LeakageBoundary {
                    cell,
                    conductance_m2: vec![0.01; 2],
                })
                .collect(),
        },
        DiffusionLimits::default(),
    )?);
    let model = Arc::new(NeutronKineticModel::new(
        "illustrative-kinetic-heat-partition".into(),
        Arc::clone(&diffusion),
        energies.to_vec(),
        vec![0.2; 4],
    )?);
    let history = FusionNeutronHistory::project(
        Arc::clone(&model),
        entries,
        FusionNeutronHistoryLimits::default(),
    )?;
    assert!(Arc::ptr_eq(history.model(), &model) && Arc::ptr_eq(diffusion.groups(), &groups));
    println!(
        "Original groups: {}; speeds {speeds:?} m/s",
        groups.data_id()
    );
    println!("Capture 0.1/m; face/leak conductances 0.02/0.01m2: illustrative inputs.");
    let budget = NeutronKineticHistoryBudget {
        transport: NeutronKineticBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-7 * emitted[0].max(1.0),
            },
            energy_tolerance_j: 1e-7 * emitted[1].max(1.0),
        },
        // Fast-neutron transport has a much shorter clock than the plasma burn.
        maximum_step_s: 2e-5,
        maximum_substeps: 16_384,
    };
    let initial = NeutronKineticState::new(Arc::clone(&model), vec![0.0; 4])?;
    let end = Time::new(0.020_020); // Both burns have ended; transport continues for 20 microseconds.
    let (neutrons, report) =
        initial.trial_with_history(Time::new(0.0), end, &history, budget, || Ok(()))?;
    let cells = neutrons.cumulative_cell_released();
    let heat = cells.iter().map(|e| e.deposited_j).sum::<f64>();
    let radiation = cells.iter().map(|e| e.radiation_j).sum::<f64>();
    let leakage = cells.iter().map(|e| e.escaped_j).sum::<f64>();
    let stored = neutrons.kinetic_energy_j()?;
    let source_residual = report.injected_kinetic_energy_j - emitted[1];
    let count_residual = report.injected_neutrons - emitted[0];
    let transport_residual = stored + heat + radiation + leakage - report.injected_kinetic_energy_j;
    let combined = plasma_remainder_j + stored + heat + radiation + leakage;
    assert!(count_residual.abs() <= budget.transport.diffusion.population_tolerance);
    assert!(source_residual.abs() <= budget.transport.energy_tolerance_j);
    assert!(transport_residual.abs() <= budget.transport.energy_tolerance_j);
    assert!(combined.abs() <= 1e-7 + budget.transport.energy_tolerance_j);
    assert!(history.rates_at(end)?.iter().all(|&rate| rate == 0.0));
    assert!(Arc::ptr_eq(neutrons.model(), &model));
    assert!(Arc::ptr_eq(neutrons.neutrons().model(), &diffusion));
    println!("time_s,stored_J,heat_J,radiation_J,leakage_J,source_residual_J,transport_residual_J");
    println!(
        "{:.6},{stored:.9e},{heat:.9e},{radiation:.9e},{leakage:.9e},{source_residual:.3e},{transport_residual:.3e}",
        end.value()
    );
    println!(
        "Source count residual: {count_residual:.3e}; plasma + neutron + exported-energy residual: {combined:.3e} J"
    );
    Ok(())
}
