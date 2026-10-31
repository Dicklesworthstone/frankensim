use super::*;
use fs_gas::material::transport::parcel::adaptive::{
    ParcelMolecularAdvanceControls, ParcelMolecularAdvanceError,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::DiffusionNetworkError;

fn accuracy() -> ParcelMolecularAdvanceControls {
    ParcelMolecularAdvanceControls {
        maximum_step_s: 20.0,
        minimum_step_s: 1e-6,
        maximum_evaluations: 5000,
        species_tolerance_kg: 1e-7,
        temperature_tolerance_k: 1e-2,
        relative_tolerance: 0.0,
    }
}
fn equal_density(
    bank: &ThermochemicalGasModel,
    x: [f64; 2],
    t: [f64; 2],
) -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([
        (1, parcel(bank, x[0], t[0], 1e5 * t[0] / 400.0)),
        (2, parcel(bank, x[1], t[1], 1e5 * t[1] / 400.0)),
    ])
}
fn endpoints(parcels: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels.keys().map(|&id| (id, endpoint())).collect()
}
fn differences(parcels: &BTreeMap<u64, ClosedGasParcel>) -> [f64; 2] {
    [
        mass(&parcels[&2], "A") / parcels[&2].mass_kg()
            - mass(&parcels[&1], "A") / parcels[&1].mass_kg(),
        parcels[&2].state().state().temperature - parcels[&1].state().state().temperature,
    ]
}

// Independent exact two-mode solution. Equal masses/densities and equal Cv
// conserve mean x and T; consequently the actual face T,p,x, material lambda
// and authored H difference are constant even while both endpoints evolve.
// δx'=-a(δx+s δT), δT'=-b δT+(H/Cv)δx'.
fn analytic_modes(time: f64, initial: [f64; 2]) -> [f64; 2] {
    let mean_t = 400.0;
    let density = 1e5 * MASS / (R * mean_t);
    let cv = 3.0 * R / MASS;
    let h = 100000.0;
    let a = 2.0 * link().area_m2 * DIFFUSIVITY / (1e-3 * link().distance_m);
    let b = 2.0 * link().area_m2 * CONDUCTIVITY / (density * 1e-3 * cv * link().distance_m);
    let s = 0.25 * MASS * h / (R * mean_t * mean_t);
    let k = [[a, a * s], [a * h / cv, b + a * h / cv * s]];
    let mean = 0.5 * (k[0][0] + k[1][1]);
    let radius = (mean * mean - a * b).sqrt();
    let c = (radius * time).cosh();
    let sh = (radius * time).sinh() / radius;
    let decay = (-mean * time).exp();
    [
        decay * (c * initial[0] - sh * ((k[0][0] - mean) * initial[0] + k[0][1] * initial[1])),
        decay * (c * initial[1] - sh * (k[1][0] * initial[0] + (k[1][1] - mean) * initial[1])),
    ]
}

#[test]
fn adaptive_soret_dufour_relaxation_converges_in_both_independent_physical_scales() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let original = equal_density(&bank, [0.7, 0.3], [390.0, 410.0]);
    let exact = analytic_modes(40.0, differences(&original));
    let mut errors = Vec::new();
    let mut work = Vec::new();
    // Tighten either finite species or finite thermal accuracy independently.
    for (species, temperature) in [(1e-6, 1.0), (1e-9, 1.0), (1.0, 1e-4), (1e-9, 1e-4)] {
        let (_, mut op) = fixture(MASS, 0.0, 0.0);
        let next = op
            .advance_candidate(
                &original,
                &[link()],
                GasTransportPaths::default(),
                [0.0, 40.0],
                &endpoints(&original),
                controls(),
                ParcelMolecularAdvanceControls {
                    species_tolerance_kg: species,
                    temperature_tolerance_k: temperature,
                    ..accuracy()
                },
                || false,
            )
            .unwrap();
        let actual = differences(next.parcels());
        errors.push([(actual[0] - exact[0]).abs(), (actual[1] - exact[1]).abs()]);
        let report = next.report();
        work.push(report.evaluations);
        assert_eq!(report.interval_s, [0.0, 40.0]);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.evaluations >= 3 * report.accepted_intervals);
        assert!(
            report
                .maximum_accepted_error_ratios
                .iter()
                .all(|&r| r <= 1.0)
        );
        assert!(
            report.window.energy_change_j.abs() <= controls().physical.exchange.energy_tolerance_j
        );
        assert!(
            report.window.maximum_species_mass_residual_kg
                <= controls().physical.exchange.mass_tolerance_kg
        );
        assert!(report.window.entropy_change_j_k > 0.0);
        // Gross accepted A transfer is its actual finite inventory change;
        // no rejected coarse or fine probe may enter this account.
        let gross: f64 = report
            .transfers
            .iter()
            .filter(|row| row.species == id("A"))
            .map(|row| {
                if row.donor == 1 {
                    row.mass_kg
                } else {
                    -row.mass_kg
                }
            })
            .sum();
        close(
            gross,
            mass(&original[&1], "A") - mass(&next.parcels()[&1], "A"),
            2e-17,
        );
    }
    assert!(
        errors[1][0] < 0.2 * errors[0][0],
        "species errors: {errors:?}"
    );
    assert!(
        errors[2][1] < 0.2 * errors[0][1],
        "thermal errors: {errors:?}"
    );
    assert!(errors[3][0] < 0.2 * errors[0][0] && errors[3][1] < 0.2 * errors[0][1]);
    assert!(work[1] > work[0] && work[2] > work[0]);
}

#[test]
fn adaptive_native_paths_replay_only_accepted_half_steps_and_keep_gross_heat() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 400.0, 1e5)),
        (2, parcel(&bank, 0.2, 400.0, 1e5)),
        (3, parcel(&bank, 0.5, 400.0, 8e4)),
    ]);
    let solvers = endpoints(&original);
    let nozzles = [SnapshotNozzleLink {
        id: 91,
        parcels: [1, 3],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 1.0,
        },
        nozzle_solve: endpoint(),
    }];
    let paths = GasTransportPaths {
        nozzles: &nozzles,
        ..GasTransportPaths::default()
    };
    let next = op
        .advance_candidate(
            &original,
            &[link()],
            paths,
            [0.0, 0.3],
            &solvers,
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 0.1,
                species_tolerance_kg: 1.0,
                temperature_tolerance_k: 1.0,
                ..accuracy()
            },
            || false,
        )
        .unwrap();
    let report = next.report();
    assert_eq!(report.accepted_intervals, 3);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.evaluations, 9);
    let (_, mut replay) = fixture(MASS, 0.0, 0.0);
    let mut state = original.clone();
    let mut streams: BTreeMap<(u64, u64, SpeciesId), [f64; 2]> = BTreeMap::new();
    let mut heat: BTreeMap<[u64; 2], [f64; 2]> = BTreeMap::new();
    let mut last = None;
    for (start, end) in [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3)] {
        let middle = start + 0.5 * (end - start);
        for dt in [middle - start, end - middle] {
            let next = replay
                .step_candidate(&state, &[link()], paths, dt, &solvers, controls(), || false)
                .unwrap();
            for row in &next.step().transport.exchange.fluxes {
                let sum = streams
                    .entry((
                        row.transfer.donor,
                        row.transfer.receiver,
                        row.transfer.species.clone(),
                    ))
                    .or_default();
                sum[0] += row.transfer.mass_kg;
                sum[1] += row.carried_enthalpy_j;
            }
            for row in &next.step().transport.exchange.heat {
                heat.entry(row.parcels).or_default()
                    [usize::from(row.heat_left_to_right_j < 0.0)] += row.heat_left_to_right_j.abs();
            }
            let (parcels, step) = next.into_parts();
            state = parcels;
            last = Some(step);
        }
    }
    assert_eq!(next.parcels(), &state);
    assert_eq!(&report.last_substep, last.as_ref().unwrap());
    assert_eq!(report.transfers.len(), streams.len());
    for row in &report.transfers {
        let sum = streams[&(row.donor, row.receiver, row.species.clone())];
        close(row.mass_kg, sum[0], 1e-18);
        close(row.carried_enthalpy_j, sum[1], 1e-10);
    }
    for row in &report.heat {
        let sum = heat[&row.parcels];
        close(row.left_to_right_j, sum[0], 1e-12);
        close(row.right_to_left_j, sum[1], 1e-12);
    }
    assert!(report.last_substep.transport.nozzles[0].transferred_mass_kg > 0.0);
}

#[test]
fn adaptive_entropy_and_withdrawal_refusals_refine_without_reclassifying_sources() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = equal_density(&bank, [0.5, 0.5], [390.0, 410.0]);
    let solvers = endpoints(&original);
    assert!(matches!(
        op.step_candidate(
            &original,
            &[link()],
            GasTransportPaths::default(),
            200.0,
            &solvers,
            controls(),
            || false
        ),
        Err(ParcelMolecularError::Transport(
            GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease { .. })
        ))
    ));
    let next = op
        .advance_candidate(
            &original,
            &[link()],
            GasTransportPaths::default(),
            [0.0, 200.0],
            &solvers,
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 200.0,
                species_tolerance_kg: 1.0,
                temperature_tolerance_k: 100.0,
                ..accuracy()
            },
            || false,
        )
        .unwrap();
    assert!(next.report().rejected_intervals > 0);
    assert!(next.report().window.entropy_change_j_k > 0.0);
    let gradient = equal_density(&bank, [0.8, 0.2], [400.0, 400.0]);
    let mut limited = controls();
    limited.physical.exchange.maximum_species_fraction = 1e-3;
    let next = op
        .advance_candidate(
            &gradient,
            &[link()],
            GasTransportPaths::default(),
            [0.0, 0.5],
            &endpoints(&gradient),
            limited,
            ParcelMolecularAdvanceControls {
                maximum_step_s: 0.5,
                species_tolerance_kg: 1.0,
                temperature_tolerance_k: 100.0,
                ..accuracy()
            },
            || false,
        )
        .unwrap();
    assert!(next.report().rejected_intervals > 0);
    assert_eq!(next.report().interval_s, [0.0, 0.5]);
}

#[test]
fn adaptive_stiff_fourier_heat_refines_the_caloric_bracket_with_no_species_transfer() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = equal_density(&bank, [1.0, 1.0], [350.0, 450.0]);
    let mut solvers = endpoints(&original);
    for solve in solvers.values_mut() {
        solve.maximum_temperature_k = 600.0;
    }
    let face = ParcelMolecularLink {
        area_m2: 1.0,
        ..link()
    };
    let refused = op.step_candidate(
        &original,
        &[face],
        GasTransportPaths::default(),
        10.0,
        &solvers,
        controls(),
        || false,
    );
    assert!(matches!(
        refused,
        Err(ParcelMolecularError::Transport(
            GasTransportError::Exchange(SpeciesNetworkError::Parcel {
                error: fs_material::thermochemical_gas::energy::GasEnergyError::OutsideBracket { .. },
                ..
            })
        ))
    ));
    let next = op
        .advance_candidate(
            &original,
            &[face],
            GasTransportPaths::default(),
            [0.0, 10.0],
            &solvers,
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 10.0,
                temperature_tolerance_k: 0.01,
                ..accuracy()
            },
            || false,
        )
        .unwrap();
    let capacity = original[&1].mass_kg() * 3.0 * R / MASS;
    let rate = 2.0 * face.area_m2 * CONDUCTIVITY / (face.distance_m * capacity);
    let difference = 100.0 * (-rate * 10.0).exp();
    close(
        next.parcels()[&1].state().state().temperature,
        400.0 - 0.5 * difference,
        0.02,
    );
    close(
        next.parcels()[&2].state().state().temperature,
        400.0 + 0.5 * difference,
        0.02,
    );
    assert!(next.report().transfers.is_empty());
    assert!(next.report().rejected_intervals > 0);
    assert!(next.report().window.entropy_change_j_k > 0.0);
    for id in [1, 2] {
        close(
            mass(&next.parcels()[&id], "A"),
            mass(&original[&id], "A"),
            1e-18,
        );
        assert_eq!(mass(&next.parcels()[&id], "B"), 0.0);
    }
}

#[test]
fn adaptive_late_source_and_work_failures_publish_neither_prefix_nor_clock() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 649.95, 1e5)),
        (2, parcel(&bank, 0.2, 649.95, 1e5)),
    ]);
    let mut state = original.clone();
    let mut time = 7.0;
    let result = op.advance_to(
        &mut state,
        &mut time,
        7.25,
        &[link()],
        GasTransportPaths::default(),
        &endpoints(&original),
        controls(),
        ParcelMolecularAdvanceControls {
            maximum_step_s: 0.025,
            species_tolerance_kg: 1.0,
            temperature_tolerance_k: 1.0,
            ..accuracy()
        },
        || false,
    );
    assert!(
        matches!(result, Err(ParcelMolecularAdvanceError::Physical { time_s, error:
        ParcelMolecularError::Face(DissipationError::Material(_)), .. }) if time_s > 7.0),
        "{result:?}"
    );
    assert_eq!(state, original);
    assert_eq!(time, 7.0);
    let original = equal_density(&bank, [0.7, 0.3], [390.0, 410.0]);
    let mut state = original.clone();
    let mut time = 0.0;
    let result = op.advance_to(
        &mut state,
        &mut time,
        0.03,
        &[link()],
        GasTransportPaths::default(),
        &endpoints(&original),
        controls(),
        ParcelMolecularAdvanceControls {
            maximum_step_s: 0.01,
            maximum_evaluations: 4,
            species_tolerance_kg: 1.0,
            temperature_tolerance_k: 1.0,
            ..accuracy()
        },
        || false,
    );
    assert!(
        matches!(result, Err(ParcelMolecularAdvanceError::WorkLimit {
        time_s, accepted_intervals: 1, evaluations: 4 }) if time_s > 0.0)
    );
    assert_eq!(state, original);
    assert_eq!(time, 0.0);
}

#[test]
fn adaptive_late_cancellation_and_unresolved_clocks_remain_atomic() {
    let (bank, mut first) = fixture(MASS, 0.0, 0.0);
    let original = equal_density(&bank, [0.7, 0.3], [390.0, 410.0]);
    let solvers = endpoints(&original);
    let control = ParcelMolecularAdvanceControls {
        maximum_step_s: 0.01,
        species_tolerance_kg: 1.0,
        temperature_tolerance_k: 1.0,
        ..accuracy()
    };
    let mut first_window_polls = 0;
    first
        .advance_candidate(
            &original,
            &[link()],
            GasTransportPaths::default(),
            [0.0, 0.01],
            &solvers,
            controls(),
            control,
            || {
                first_window_polls += 1;
                false
            },
        )
        .unwrap();
    let (_, mut op) = fixture(MASS, 0.0, 0.0);
    let mut state = original.clone();
    let mut time = 0.0;
    let mut polls = 0;
    // The same first interval is complete before this checkpoint; no accepted
    // prefix of the longer requested window may be exposed on cancellation.
    let result = op.advance_to(
        &mut state,
        &mut time,
        1.0,
        &[link()],
        GasTransportPaths::default(),
        &solvers,
        controls(),
        control,
        || {
            polls += 1;
            polls > first_window_polls
        },
    );
    assert!(matches!(
        result,
        Err(ParcelMolecularAdvanceError::Cancelled)
    ));
    assert_eq!(state, original);
    assert_eq!(time, 0.0);
    let clock = 2.0_f64.powi(50);
    let mut time = clock;
    let result = op.advance_to(
        &mut state,
        &mut time,
        clock.next_up(),
        &[link()],
        GasTransportPaths::default(),
        &solvers,
        controls(),
        ParcelMolecularAdvanceControls {
            maximum_step_s: 1.0,
            ..control
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(ParcelMolecularAdvanceError::TimeResolution { .. })
    ));
    assert_eq!(state, original);
    assert_eq!(time.to_bits(), clock.to_bits());
    let equilibrium = equal_density(&bank, [0.5, 0.5], [400.0, 400.0]);
    let mut state = equilibrium.clone();
    let report = op
        .advance_to(
            &mut state,
            &mut time,
            clock + 0.75,
            &[link()],
            GasTransportPaths::default(),
            &endpoints(&equilibrium),
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 0.75,
                ..control
            },
            || false,
        )
        .unwrap();
    assert_eq!(time.to_bits(), (clock + 0.75).to_bits());
    assert_eq!(state, equilibrium);
    assert_eq!(report.minimum_accepted_substep_s, 0.25);
    assert_eq!(report.maximum_accepted_substep_s, 0.5);
    // Decimal endpoint rounding must not manufacture an unresolvable one-ULP
    // remainder merely because an interior span rounds just above its cap.
    let mut time = 0.0;
    let report = op
        .advance_to(
            &mut state,
            &mut time,
            0.4,
            &[link()],
            GasTransportPaths::default(),
            &endpoints(&equilibrium),
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 0.1,
                ..control
            },
            || false,
        )
        .unwrap();
    assert_eq!(time.to_bits(), 0.4_f64.to_bits());
    assert_eq!(report.accepted_intervals, 4);
    assert_eq!(state, equilibrium);
    // Three rounded ordinary spans would otherwise stop one representable
    // value before the target. The final span itself must include that tail.
    let mut time = 0.0;
    let target = 0.005_f64;
    let report = op
        .advance_to(
            &mut state,
            &mut time,
            target,
            &[link()],
            GasTransportPaths::default(),
            &endpoints(&equilibrium),
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: (target / 3.0).next_down(),
                ..control
            },
            || false,
        )
        .unwrap();
    assert_eq!(time.to_bits(), target.to_bits());
    assert_eq!(report.accepted_intervals, 3);
    assert_eq!(state, equilibrium);
}

#[test]
fn adaptive_reference_shifts_preserve_the_finite_trajectory_and_window_audit() {
    let run = |formation, shift| {
        let (bank, mut op) = fixture(MASS, formation, shift);
        let original = equal_density(&bank, [0.7, 0.3], [390.0, 410.0]);
        op.advance_candidate(
            &original,
            &[link()],
            GasTransportPaths::default(),
            [0.0, 20.0],
            &endpoints(&original),
            controls(),
            ParcelMolecularAdvanceControls {
                maximum_step_s: 1.0,
                ..accuracy()
            },
            || false,
        )
        .unwrap()
    };
    let baseline = run(0.0, 0.0);
    let shifted = run(1e14, 2.0_f64.powi(40));
    for id in [1, 2] {
        close(
            baseline.parcels()[&id].state().state().temperature,
            shifted.parcels()[&id].state().state().temperature,
            1e-8,
        );
        close(
            mass(&baseline.parcels()[&id], "A"),
            mass(&shifted.parcels()[&id], "A"),
            1e-16,
        );
    }
    close(
        baseline.report().window.energy_change_j,
        shifted.report().window.energy_change_j,
        1e-8,
    );
    close(
        baseline.report().window.entropy_change_j_k,
        shifted.report().window.entropy_change_j_k,
        1e-10,
    );
}
