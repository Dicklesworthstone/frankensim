//! G1 piecewise-source absorption and G3 restart/clock/spatial composition.
//! Coefficients and joule-scale groups are manufactured numerical fixtures.

use super::*;
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind,
};
use crate::irradiation::NeutronGroups;
use fs_qty::{Amount, Energy};
use fs_thermochem::{
    electron_collision::AVOGADRO_PER_MOL, fusion::neutrons::FusionNeutronEmissionV1,
};

fn model(cells: usize) -> Arc<NeutronKineticModel> {
    let reactions = (0..cells)
        .flat_map(|cell| {
            (0..2).map(move |group| DiffusionReaction {
                data_id: format!("manufactured-capture-{cell}-{group}"),
                cell,
                incident_group: group,
                macroscopic_cross_section_m_inv: 0.7,
                kind: DiffusionReactionKind::Capture,
            })
        })
        .collect::<Vec<_>>();
    let diffusion = DiffusionModel::new(
        DiffusionData {
            data_id: "manufactured-history".into(),
            groups: Arc::new(
                NeutronGroups::new("manufactured-energy".into(), vec![0.0, 2.5, 8.0]).unwrap(),
            ),
            speeds_m_s: vec![1.0, 1.0],
            cell_volumes_m3: vec![1.0; cells],
            delayed_groups: vec![],
            reactions,
            links: vec![],
            leakage: vec![],
        },
        DiffusionLimits::default(),
    )
    .unwrap();
    Arc::new(
        NeutronKineticModel::new(
            "manufactured-kinetic-account".into(),
            Arc::new(diffusion),
            vec![1.0, 4.0],
            vec![0.75; cells * 2],
        )
        .unwrap(),
    )
}

fn entry(cell: usize, start: f64, end: f64, count: f64, energy: f64) -> FusionNeutronHistoryEntry {
    let emission = FusionNeutronEmissionV1::new(
        Time::new(end - start),
        [Amount::new(count / AVOGADRO_PER_MOL), Amount::new(0.0)],
        [Energy::new(count * energy), Energy::new(0.0)],
    )
    .unwrap();
    FusionNeutronHistoryEntry {
        cell,
        interval: FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), emission).unwrap(),
    }
}

fn budget() -> NeutronKineticHistoryBudget {
    NeutronKineticHistoryBudget {
        transport: NeutronKineticBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-7,
            },
            energy_tolerance_j: 1e-7,
        },
        maximum_step_s: 10.0,
        maximum_substeps: 100,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1.0),
        "actual={actual:e} expected={expected:e}"
    );
}

fn fixtures() -> Vec<FusionNeutronHistoryEntry> {
    vec![
        entry(1, 0.5, 1.5, 80.0, 2.5),
        entry(0, 1.0, 2.0, 50.0, 3.0),
        entry(0, 0.0, 1.0, 100.0, 2.0),
    ]
}

#[test]
fn prepared_history_retains_original_multicell_intervals_and_right_limits() {
    let entries = fixtures();
    let history = FusionNeutronHistory::project(
        model(2),
        entries.clone(),
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    assert_eq!(history.entries(), entries);
    for (time, rate) in [
        (0.0, 100.0),
        (0.5, 180.0),
        (1.0, 130.0),
        (1.5, 50.0),
        (2.0, 0.0),
        (5.0, 0.0),
    ] {
        close(
            history.rates_at(Time::new(time)).unwrap().iter().sum(),
            rate,
        );
    }
    assert_eq!(history.next_after(Time::new(0.5)), Some(Time::new(1.0)));
    assert_eq!(history.next_after(Time::new(2.0)), None);
    let (number, energy) = history
        .interval_moments(Time::new(0.75), Time::new(0.25))
        .unwrap();
    close(number, 45.0);
    close(energy, 100.0);
    assert!(
        history
            .interval_moments(Time::new(0.75), Time::new(0.5))
            .is_err()
    );
}

#[test]
fn continuous_history_matches_analytic_capture_and_continues_after_shutoff() {
    let owner = model(2);
    let entries = fixtures();
    let history = FusionNeutronHistory::project(
        owner.clone(),
        entries.clone(),
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![0.0; 4]).unwrap();
    let (final_state, report) = initial
        .trial_with_history(
            Time::new(0.0),
            Time::new(4.0),
            &history,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let mut stored = 0.0;
    let mut injected = 0.0;
    let mut populations = [0.0; 2];
    for entry in &entries {
        let emission = entry.interval.emission();
        let duration = emission.duration().value();
        let rate = emission.total_particle_count() / duration;
        let number = rate / 0.7
            * -(-0.7 * duration).exp_m1()
            * (-0.7 * (4.0 - entry.interval.end().value())).exp();
        populations[entry.cell] += number;
        let mean = emission.total_kinetic_energy().value() / emission.total_particle_count();
        stored += number * mean;
        injected += emission.total_kinetic_energy().value();
    }
    for (cell, expected) in populations.into_iter().enumerate() {
        close(
            final_state.neutrons().neutron_counts()[cell * 2..cell * 2 + 2]
                .iter()
                .sum(),
            expected,
        );
    }
    close(report.injected_neutrons, 230.0);
    close(report.injected_kinetic_energy_j, injected);
    close(report.final_kinetic_energy_j, stored);
    close(report.released.deposited_j, 0.75 * (injected - stored));
    close(report.released.radiation_j, 0.25 * (injected - stored));
    assert_eq!(report.released.escaped_j, 0.0);
    assert!(report.energy_balance_residual_j.abs() < 1e-8);
    assert_eq!(report.substeps, 5);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0; 4]);
}

#[test]
fn restart_inside_an_active_interval_preserves_source_and_energy() {
    let owner = model(2);
    let history = FusionNeutronHistory::project(
        owner.clone(),
        fixtures(),
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![0.0; 4]).unwrap();
    let (whole, report) = initial
        .trial_with_history(
            Time::new(0.0),
            Time::new(4.0),
            &history,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let (middle, first) = initial
        .trial_with_history(Time::new(0.0), Time::new(1.25), &history, budget(), || {
            Ok(())
        })
        .unwrap();
    let (split, second) = middle
        .trial_with_history(Time::new(1.25), Time::new(4.0), &history, budget(), || {
            Ok(())
        })
        .unwrap();
    for (a, b) in whole
        .neutrons()
        .neutron_counts()
        .iter()
        .zip(split.neutrons().neutron_counts())
    {
        close(*a, *b);
    }
    close(
        first.injected_neutrons + second.injected_neutrons,
        report.injected_neutrons,
    );
    close(
        first.injected_kinetic_energy_j + second.injected_kinetic_energy_j,
        report.injected_kinetic_energy_j,
    );
    close(
        first.released.deposited_j + second.released.deposited_j,
        report.released.deposited_j,
    );
    for (a, b) in whole
        .cumulative_cell_released()
        .iter()
        .zip(split.cumulative_cell_released())
    {
        close(a.deposited_j, b.deposited_j);
    }
}

#[test]
fn late_births_retain_more_neutrons_than_one_whole_burn_average() {
    let owner = model(1);
    let late = FusionNeutronHistory::project(
        owner.clone(),
        vec![entry(0, 0.0, 0.5, 10.0, 2.0), entry(0, 0.5, 1.0, 90.0, 2.0)],
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let flat = FusionNeutronHistory::project(
        owner.clone(),
        vec![entry(0, 0.0, 1.0, 100.0, 2.0)],
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![0.0; 2]).unwrap();
    let (a, ar) = initial
        .trial_with_history(Time::new(0.0), Time::new(1.0), &late, budget(), || Ok(()))
        .unwrap();
    let (b, br) = initial
        .trial_with_history(Time::new(0.0), Time::new(1.0), &flat, budget(), || Ok(()))
        .unwrap();
    close(ar.injected_neutrons, br.injected_neutrons);
    close(ar.injected_kinetic_energy_j, br.injected_kinetic_energy_j);
    assert!(a.kinetic_energy_j().unwrap() > b.kinetic_energy_j().unwrap() + 1.0);
    assert!(ar.released.deposited_j < br.released.deposited_j);
}

#[test]
fn overlapping_domains_storage_and_cancellation_are_atomic_refusals() {
    let owner = model(2);
    assert!(
        FusionNeutronHistory::project(
            owner.clone(),
            vec![entry(0, 0.0, 1.0, 1.0, 2.0), entry(0, 0.5, 2.0, 1.0, 2.0)],
            FusionNeutronHistoryLimits::default()
        )
        .is_err()
    );
    assert!(
        FusionNeutronHistory::project(
            owner.clone(),
            vec![entry(2, 0.0, 1.0, 1.0, 2.0)],
            FusionNeutronHistoryLimits::default()
        )
        .is_err()
    );
    assert!(
        FusionNeutronHistory::project(
            owner.clone(),
            vec![entry(0, 0.0, 1.0, 1.0, 20.0)],
            FusionNeutronHistoryLimits::default()
        )
        .is_err()
    );
    assert!(
        FusionNeutronHistory::project(
            owner.clone(),
            fixtures(),
            FusionNeutronHistoryLimits {
                maximum_intervals: 2,
                maximum_rate_values: 1000
            }
        )
        .is_err()
    );
    assert!(
        FusionNeutronHistory::project(
            owner.clone(),
            fixtures(),
            FusionNeutronHistoryLimits {
                maximum_intervals: 3,
                maximum_rate_values: 10
            }
        )
        .is_err()
    );
    let history = FusionNeutronHistory::project(
        owner.clone(),
        fixtures(),
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![0.0; 4]).unwrap();
    let mut short = budget();
    short.maximum_substeps = 1;
    assert!(
        initial
            .trial_with_history(Time::new(0.0), Time::new(4.0), &history, short, || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    assert!(
        initial
            .trial_with_history(Time::new(0.0), Time::new(4.0), &history, budget(), || {
                polls += 1;
                if polls > 8 {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0; 4]);
    assert_eq!(
        initial.cumulative_cell_released(),
        vec![EnergyChannels::default(); 2]
    );
}

#[test]
fn empty_forcing_still_evolves_resident_neutrons_and_admits_idle_calls() {
    let owner = model(1);
    let history =
        FusionNeutronHistory::project(owner.clone(), vec![], FusionNeutronHistoryLimits::default())
            .unwrap();
    assert!(history.rates_at(Time::new(f64::NAN)).is_err());
    let initial = NeutronKineticState::new(owner, vec![4.0, 6.0]).unwrap();
    let (idle, zero) = initial
        .trial_with_history(
            Time::new(2.0),
            Time::new(2.0),
            &history,
            budget(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        idle.neutrons().neutron_counts(),
        initial.neutrons().neutron_counts()
    );
    assert_eq!(zero.substeps, 0);
    let (later, report) = initial
        .trial_with_history(
            Time::new(2.0),
            Time::new(3.0),
            &history,
            budget(),
            || Ok(()),
        )
        .unwrap();
    close(later.kinetic_energy_j().unwrap(), 28.0 * (-0.7_f64).exp());
    assert_eq!(report.injected_neutrons, 0.0);
    assert!(report.released.deposited_j > 0.0);
}

#[test]
fn bounded_transport_steps_resolve_stiff_spatial_cycles_and_source_shutoff() {
    use crate::diffusion::{DiffusionLink, LeakageBoundary};

    // Symmetric sources and material data keep both populations equal, so the
    // reciprocal spatial transfers cancel in the independent scalar solution:
    // dN/dt = source - lambda*N. The long source interval spans 1024 removal
    // lifetimes; the source-off tail spans two more. All times are exact binary
    // fractions, avoiding a tiny final remainder in this stiffness fixture.
    let lambda = 8192.0;
    let duration = 0.125;
    let step = 1.0 / lambda;
    let end = duration + 2.0 * step;
    let entries = (0..2)
        .map(|cell| entry(cell, 0.0, duration, lambda * duration, 4.0))
        .collect::<Vec<_>>();
    let emission = *entries[0].interval.emission();
    // Use the retained channel mean exactly: converting particles to moles and
    // back need not reproduce the authored 4 J representative bit for bit.
    let energy = emission.channels()[0]
        .mean_particle_energy()
        .unwrap()
        .value();
    let source_count = emission.total_particle_count();
    let resident_count = 4.0;
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-stiff-spatial-cycle".into(),
                groups: Arc::new(
                    NeutronGroups::new("manufactured-stiff-energy".into(), vec![0.0, 8.0]).unwrap(),
                ),
                speeds_m_s: vec![lambda],
                cell_volumes_m3: vec![1.0; 2],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|cell| DiffusionReaction {
                        data_id: format!("manufactured-stiff-capture-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 0.75,
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![DiffusionLink {
                    cell_a: 0,
                    cell_b: 1,
                    conductance_m2: vec![1.0],
                }],
                leakage: (0..2)
                    .map(|cell| LeakageBoundary {
                        cell,
                        conductance_m2: vec![0.25],
                    })
                    .collect(),
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let owner = Arc::new(
        NeutronKineticModel::new(
            "manufactured-stiff-kinetic-account".into(),
            diffusion,
            vec![energy],
            vec![0.5; 2],
        )
        .unwrap(),
    );
    let history = FusionNeutronHistory::project(
        owner.clone(),
        entries,
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![resident_count; 2]).unwrap();
    let mut limits = budget();
    limits.maximum_substeps = 2048;
    // A source-boundary-only advance would require unrepresentable surviving
    // coefficients in the original positive propagator, which must refuse.
    assert!(
        initial
            .trial_with_history(Time::new(0.0), Time::new(end), &history, limits, || Ok(()))
            .is_err()
    );
    limits.maximum_step_s = step;
    let (final_state, report) = initial
        .trial_with_history(Time::new(0.0), Time::new(end), &history, limits, || Ok(()))
        .unwrap();
    let remaining = (resident_count * (-lambda * duration).exp()
        - source_count / (lambda * duration) * (-lambda * duration).exp_m1())
        * (-2.0_f64).exp();
    assert_eq!(report.substeps, 1026);
    close(report.injected_neutrons, 2.0 * source_count);
    close(
        report.injected_kinetic_energy_j,
        2.0 * emission.total_kinetic_energy().value(),
    );
    for (cell, &population) in final_state.neutrons().neutron_counts().iter().enumerate() {
        close(population, remaining);
        let removed_energy = energy * (resident_count + source_count - remaining);
        close(
            report.cell_released[cell].deposited_j,
            0.375 * removed_energy,
        );
        close(
            report.cell_released[cell].radiation_j,
            0.375 * removed_energy,
        );
        close(report.cell_released[cell].escaped_j, 0.25 * removed_energy);
    }
    close(report.final_kinetic_energy_j, 2.0 * energy * remaining);
    assert!(report.energy_roundoff_j <= limits.transport.energy_tolerance_j);
    assert!(report.population_roundoff <= limits.transport.diffusion.population_tolerance);
    assert_eq!(history.rates_at(Time::new(end)).unwrap(), &[0.0; 2]);
    assert_eq!(initial.neutrons().neutron_counts(), &[resident_count; 2]);
}
