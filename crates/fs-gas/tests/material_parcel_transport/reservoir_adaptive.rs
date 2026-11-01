//! Independent transient and gross-stream checks for open molecular networks.

use super::*;
use fs_gas::material::transport::parcel::adaptive::{
    ParcelMolecularAdvanceControls, ParcelMolecularAdvanceError,
};
use fs_gas::material::transport::parcel::reservoir::{
    MolecularReservoir, ReservoirMolecularLink, adaptive::ReservoirMolecularAdvanceControls,
};

const CAPACITY: f64 = 0.75;
const CONDUCTANCE: f64 = 0.05;

fn accuracy(temperature: f64, heat: f64) -> ReservoirMolecularAdvanceControls {
    ReservoirMolecularAdvanceControls {
        state: ParcelMolecularAdvanceControls {
            maximum_step_s: 2.0,
            minimum_step_s: 1e-6,
            maximum_evaluations: 10_000,
            species_tolerance_kg: 1e-7,
            temperature_tolerance_k: temperature,
            relative_tolerance: 0.0,
        },
        species_transfer_tolerance_kg: 1e-8,
        heat_transfer_tolerance_j: heat,
    }
}

fn bath(bank: &ThermochemicalGasModel, x: f64, t: f64, p: f64) -> MolecularReservoir {
    MolecularReservoir::new(
        bank.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        point(t, p),
        || false,
    )
    .unwrap()
}

fn pure_pair(bank: &ThermochemicalGasModel) -> BTreeMap<u64, ClosedGasParcel> {
    // Constant Cv=3R/M and p/T=100000/400 give C=3*p_ref*V/T_ref
    // =0.75 J/K in both original 1e-3 m3 rigid cells.
    BTreeMap::from([
        (1, parcel(bank, 1.0, 350.0, 1e5 * 350.0 / 400.0)),
        (2, parcel(bank, 1.0, 450.0, 1e5 * 450.0 / 400.0)),
    ])
}

fn endpoints(gas: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gas.keys().map(|&id| (id, endpoint())).collect()
}

fn internal() -> ParcelMolecularLink {
    ParcelMolecularLink {
        id: 7,
        parcels: [1, 2],
        area_m2: 0.1,
        distance_m: 0.1,
    }
}

fn boundary(reservoir: &MolecularReservoir) -> ReservoirMolecularLink<'_> {
    ReservoirMolecularLink {
        id: 21,
        parcel: 1,
        reservoir,
        area_m2: 0.1,
        distance_m: 0.1,
    }
}

fn temperatures(gas: &BTreeMap<u64, ClosedGasParcel>) -> [f64; 2] {
    [
        gas[&1].state().state().temperature,
        gas[&2].state().state().temperature,
    ]
}

// Independent continuum reference, with no native trajectory or caloric
// inversion: u_i=T_i-500; u'=-K*u, K=G/C * [[2,-1],[-1,1]].
// Cp=4R follows the authored NASA coefficient, Cv=3R in each rigid cell.
// Pure A has no species diffusion, so Soret/Dufour contributions vanish.
fn thermal_reference(time: f64) -> ([f64; 2], f64) {
    let slow = CONDUCTANCE / CAPACITY * (3.0 - 5.0_f64.sqrt()) / 2.0;
    let fast = CONDUCTANCE / CAPACITY * (3.0 + 5.0_f64.sqrt()) / 2.0;
    let initial = [-150.0, -50.0];
    let k_initial = [
        CONDUCTANCE / CAPACITY * (2.0 * initial[0] - initial[1]),
        CONDUCTANCE / CAPACITY * (initial[1] - initial[0]),
    ];
    let mut slow_modes = [0.0; 2];
    let mut fast_modes = [0.0; 2];
    let mut temperature = [0.0; 2];
    for i in 0..2 {
        slow_modes[i] = (fast * initial[i] - k_initial[i]) / (fast - slow);
        fast_modes[i] = (k_initial[i] - slow * initial[i]) / (fast - slow);
        temperature[i] =
            500.0 + slow_modes[i] * (-slow * time).exp() + fast_modes[i] * (-fast * time).exp();
    }
    // Integrate G*(500-T1) independently, not by reusing a native ledger or
    // subtracting its endpoint energies.
    let heat = -CONDUCTANCE
        * (slow_modes[0] * -(-slow * time).exp_m1() / slow
            + fast_modes[0] * -(-fast * time).exp_m1() / fast);
    (temperature, heat)
}

#[test]
fn adaptive_open_internal_heat_converges_to_two_exact_modes_and_bath_integral() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let original = pure_pair(&bank);
    let saved = original.clone();
    let held = bath(&bank, 1.0, 500.0, 125_000.0);
    let saved_bath = held.clone();
    let expected = thermal_reference(10.0);
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1.0, 0.01] {
        let (_, mut op) = fixture(MASS, 0.0, 0.0);
        let next = op
            .advance_network_candidate(
                &original,
                &[internal()],
                &[boundary(&held)],
                [0.0, 10.0],
                &endpoints(&original),
                controls(),
                accuracy(tolerance, 1e6),
                || false,
            )
            .unwrap();
        let actual = temperatures(next.parcels());
        errors.push(
            (actual[0] - expected.0[0])
                .abs()
                .max((actual[1] - expected.0[1]).abs()),
        );
        let report = next.report();
        work.push(report.evaluations);
        assert_eq!(report.interval_s, [0.0, 10.0]);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.evaluations >= 3 * report.accepted_intervals);
        assert!(
            report
                .maximum_accepted_error_ratios
                .iter()
                .all(|&v| v <= 1.0)
        );
        assert!(report.transfers.is_empty());
        for row in &report.reservoir_species {
            assert_eq!(row.mass_into_gas_kg, 0.0);
            assert_eq!(row.mass_out_of_gas_kg, 0.0);
            assert_eq!(row.net_mass_into_gas_kg, 0.0);
        }
        assert_eq!(report.reservoir_heat.len(), 1);
        let external = &report.reservoir_heat[0];
        assert_eq!((external.id, external.parcel), (21, 1));
        assert!(external.heat_into_gas_j > 0.0);
        assert_eq!(external.heat_out_of_gas_j, 0.0);
        close(
            external.net_heat_into_gas_j,
            CAPACITY * ((actual[0] - 350.0) + (actual[1] - 450.0)),
            controls().physical.exchange.energy_tolerance_j,
        );
        assert_eq!(report.heat.len(), 1);
        let internal_heat = &report.heat[0];
        assert_eq!(internal_heat.parcels, [1, 2]);
        // Cell 2 initially heats cell 1; later the bath-warmed cell 1 heats
        // cell 2. Both physical directions must survive the integrated ledger.
        assert!(internal_heat.left_to_right_j > 1e-3);
        assert!(internal_heat.right_to_left_j > 1e-3);
        close(
            internal_heat.left_to_right_j - internal_heat.right_to_left_j,
            CAPACITY * (actual[1] - 450.0),
            controls().physical.exchange.energy_tolerance_j,
        );
        let entropy = CAPACITY * (actual[0] / 350.0).ln() + CAPACITY * (actual[1] / 450.0).ln()
            - external.net_heat_into_gas_j / 500.0;
        close(report.window.entropy_production_j_k, entropy, 2e-11);
        assert!(report.window.entropy_production_j_k > 0.0);
        assert!(
            report.window.energy_residual_j.abs()
                <= controls().physical.exchange.energy_tolerance_j
        );
        for id in [1, 2] {
            close(
                next.parcels()[&id].mass_kg(),
                original[&id].mass_kg(),
                1e-18,
            );
            assert_eq!(mass(&next.parcels()[&id], "B"), 0.0);
            close(original[&id].mass_kg() * 3.0 * R / MASS, CAPACITY, 2e-15);
        }
        if tolerance == 0.01 {
            close(external.net_heat_into_gas_j, expected.1, 0.3);
        }
    }
    assert!(
        errors[1] < 0.2 * errors[0],
        "temperature errors: {errors:?}"
    );
    assert!(errors[1] < 0.2, "temperature errors: {errors:?}");
    assert!(work[1] > work[0]);
    assert_eq!(original, saved);
    assert_eq!(held, saved_bath);
}

#[test]
fn bath_heat_stream_accuracy_is_independent_of_loose_state_accuracy() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let original = pure_pair(&bank);
    let held = bath(&bank, 1.0, 500.0, 125_000.0);
    let expected_heat = thermal_reference(10.0).1;
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1e6, 1e-3] {
        let (_, mut op) = fixture(MASS, 0.0, 0.0);
        let next = op
            .advance_network_candidate(
                &original,
                &[internal()],
                &[boundary(&held)],
                [0.0, 10.0],
                &endpoints(&original),
                controls(),
                accuracy(1e6, tolerance),
                || false,
            )
            .unwrap();
        errors.push((next.report().reservoir_heat[0].net_heat_into_gas_j - expected_heat).abs());
        work.push(next.report().evaluations);
        assert!(next.report().maximum_accepted_error_ratios[3] <= 1.0);
        if tolerance == 1e-3 {
            assert!(next.report().rejected_intervals > 0);
        }
    }
    assert!(
        errors[1] < 0.2 * errors[0],
        "heat-integral errors: {errors:?}"
    );
    assert!(work[1] > work[0]);

    let (_, mut op) = fixture(MASS, 0.0, 0.0);
    let mut at_floor = accuracy(1e6, 1e-8);
    at_floor.state.minimum_step_s = 2.0;
    let refusal = op.advance_network_candidate(
        &original,
        &[internal()],
        &[boundary(&held)],
        [0.0, 2.0],
        &endpoints(&original),
        controls(),
        at_floor,
        || false,
    );
    assert!(
        matches!(
            refusal,
            Err(ParcelMolecularAdvanceError::StreamAccuracy { error_ratios, .. })
                if error_ratios[1] > 1.0
        ),
        "{refusal:?}"
    );
}

#[test]
fn balanced_baths_retain_exact_gross_species_heat_and_mixing_entropy() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
    let rich = bath(&bank, 0.9, 400.0, 1e5);
    let lean = bath(&bank, 0.1, 400.0, 1e5);
    let links = [
        ReservoirMolecularLink {
            id: 21,
            area_m2: 0.01,
            ..boundary(&rich)
        },
        ReservoirMolecularLink {
            id: 23,
            area_m2: 0.01,
            ..boundary(&lean)
        },
    ];
    let mut accuracy = accuracy(1e-5, 1e-8);
    accuracy.state.maximum_step_s = 20.0;
    accuracy.species_transfer_tolerance_kg = 1e-12;
    let next = op
        .advance_network_candidate(
            &original,
            &[],
            &links,
            [0.0, 40.0],
            &endpoints(&original),
            controls(),
            accuracy,
            || false,
        )
        .unwrap();
    let density = 1e5 * MASS / (R * 400.0);
    let expected_mass = density * DIFFUSIVITY * 0.01 / 0.1 * 0.4 * 40.0;
    let report = next.report();
    assert_eq!(report.accepted_intervals, 2);
    assert_eq!(report.accepted_substeps, 4);
    assert_eq!(report.evaluations, 6);
    assert_eq!(report.reservoir_species.len(), 4);
    for (boundary_id, species, sign) in [
        (21, "A", 1.0),
        (21, "B", -1.0),
        (23, "A", -1.0),
        (23, "B", 1.0),
    ] {
        let row = report
            .reservoir_species
            .iter()
            .find(|row| row.id == boundary_id && row.species == id(species))
            .unwrap();
        assert_eq!(row.parcel, 1);
        close(row.net_mass_into_gas_kg, sign * expected_mass, 2e-17);
        close(
            row.mass_into_gas_kg,
            if sign > 0.0 { expected_mass } else { 0.0 },
            2e-17,
        );
        close(
            row.mass_out_of_gas_kg,
            if sign < 0.0 { expected_mass } else { 0.0 },
            2e-17,
        );
    }
    for (boundary_id, sign) in [(21, 1.0), (23, -1.0)] {
        let row = report
            .reservoir_heat
            .iter()
            .find(|row| row.id == boundary_id)
            .unwrap();
        close(
            row.net_heat_into_gas_j,
            sign * 100_000.0 * expected_mass,
            2e-9,
        );
        assert!(
            if sign > 0.0 {
                row.heat_into_gas_j
            } else {
                row.heat_out_of_gas_j
            } > 0.0
        );
    }
    for species in ["A", "B"] {
        close(
            mass(&next.parcels()[&1], species),
            mass(&original[&1], species),
            2e-17,
        );
    }
    close(next.parcels()[&1].state().state().temperature, 400.0, 2e-8);
    let entropy = 2.0 * expected_mass * R / MASS * 9.0_f64.ln();
    close(report.window.entropy_production_j_k, entropy, 2e-9);
    assert!(
        report.window.energy_residual_j.abs() <= controls().physical.exchange.energy_tolerance_j
    );
}

#[test]
fn active_internal_and_bath_face_order_preserves_the_complete_candidate() {
    let (bank, mut first) = fixture(MASS, 0.0, 0.0);
    let original = pure_pair(&bank);
    let held = bath(&bank, 1.0, 500.0, 125_000.0);
    // Separate exposed patches add to the same G=0.05 W/K reference.
    let mut internal_links = [
        ParcelMolecularLink {
            id: 9,
            area_m2: 0.06,
            ..internal()
        },
        ParcelMolecularLink {
            id: 7,
            area_m2: 0.04,
            ..internal()
        },
    ];
    let mut bath_links = [
        ReservoirMolecularLink {
            id: 23,
            area_m2: 0.07,
            ..boundary(&held)
        },
        ReservoirMolecularLink {
            id: 21,
            area_m2: 0.03,
            ..boundary(&held)
        },
    ];
    let accepted = first
        .advance_network_candidate(
            &original,
            &internal_links,
            &bath_links,
            [0.0, 1.0],
            &endpoints(&original),
            controls(),
            accuracy(0.01, 1e-3),
            || false,
        )
        .unwrap();
    internal_links.reverse();
    bath_links.reverse();
    let (_, mut second) = fixture(MASS, 0.0, 0.0);
    let reversed = second
        .advance_network_candidate(
            &original,
            &internal_links,
            &bath_links,
            [0.0, 1.0],
            &endpoints(&original),
            controls(),
            accuracy(0.01, 1e-3),
            || false,
        )
        .unwrap();
    assert_eq!(accepted, reversed);
    assert_eq!(accepted.report().last_substep.internal_links.len(), 2);
    assert_eq!(accepted.report().last_substep.links.len(), 2);
    assert!(
        accepted
            .report()
            .reservoir_heat
            .iter()
            .all(|row| row.heat_into_gas_j > 0.0)
    );
}

#[test]
fn open_network_work_and_late_source_failures_publish_no_state_or_clock_prefix() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = pure_pair(&bank);
    let held = bath(&bank, 1.0, 500.0, 125_000.0);
    let mut limited = accuracy(1.0, 1.0);
    limited.state.maximum_step_s = 0.02;
    limited.state.maximum_evaluations = 4;
    let mut state = original.clone();
    let mut time = 7.0;
    let result = op.advance_network_to(
        &mut state,
        &mut time,
        7.06,
        &[internal()],
        &[boundary(&held)],
        &endpoints(&original),
        controls(),
        limited,
        || false,
    );
    assert!(
        matches!(result, Err(ParcelMolecularAdvanceError::WorkLimit {
        time_s, accepted_intervals: 1, evaluations: 4,
    }) if time_s > 7.0),
        "{result:?}"
    );
    assert_eq!(state, original);
    assert_eq!(time, 7.0);

    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 649.95, 1e5)),
        (2, parcel(&bank, 0.2, 649.95, 1e5)),
    ]);
    let held = bath(&bank, 0.8, 649.95, 1e5);
    let bath_link = ReservoirMolecularLink {
        area_m2: 1e-8,
        ..boundary(&held)
    };
    let mut loose = accuracy(1.0, 1.0);
    loose.state.maximum_step_s = 0.025;
    loose.state.species_tolerance_kg = 1.0;
    loose.species_transfer_tolerance_kg = 1.0;
    let mut state = original.clone();
    let mut time = 7.0;
    let result = op.advance_network_to(
        &mut state,
        &mut time,
        7.25,
        &[link()],
        &[bath_link],
        &endpoints(&original),
        controls(),
        loose,
        || false,
    );
    assert!(
        matches!(result, Err(ParcelMolecularAdvanceError::Physical {
        time_s, error: ParcelMolecularError::Face(DissipationError::Material(_)), ..
    }) if time_s > 7.0),
        "{result:?}"
    );
    assert_eq!(state, original);
    assert_eq!(time, 7.0);
}

#[test]
fn open_network_late_cancellation_and_represented_clocks_remain_atomic() {
    let (bank, mut first) = fixture(MASS, 0.0, 0.0);
    let original = pure_pair(&bank);
    let held = bath(&bank, 1.0, 500.0, 125_000.0);
    let solvers = endpoints(&original);
    let mut control = accuracy(1.0, 1.0);
    control.state.maximum_step_s = 0.02;
    let mut first_interval_polls = 0;
    first
        .advance_network_candidate(
            &original,
            &[internal()],
            &[boundary(&held)],
            [0.0, 0.02],
            &solvers,
            controls(),
            control,
            || {
                first_interval_polls += 1;
                false
            },
        )
        .unwrap();
    let (_, mut op) = fixture(MASS, 0.0, 0.0);
    let mut state = original.clone();
    let mut time = 0.0;
    let mut polls = 0;
    let result = op.advance_network_to(
        &mut state,
        &mut time,
        1.0,
        &[internal()],
        &[boundary(&held)],
        &solvers,
        controls(),
        control,
        || {
            polls += 1;
            polls > first_interval_polls
        },
    );
    assert!(matches!(
        result,
        Err(ParcelMolecularAdvanceError::Cancelled)
    ));
    assert_eq!(state, original);
    assert_eq!(time, 0.0);

    let clock = 2.0_f64.powi(50);
    time = clock;
    control.state.maximum_step_s = 1.0;
    let result = op.advance_network_to(
        &mut state,
        &mut time,
        clock.next_up(),
        &[internal()],
        &[boundary(&held)],
        &solvers,
        controls(),
        control,
        || false,
    );
    assert!(matches!(
        result,
        Err(ParcelMolecularAdvanceError::TimeResolution { .. })
    ));
    assert_eq!(state, original);
    assert_eq!(time.to_bits(), clock.to_bits());

    let equilibrium = BTreeMap::from([
        (1, parcel(&bank, 1.0, 400.0, 1e5)),
        (2, parcel(&bank, 1.0, 400.0, 1e5)),
    ]);
    let held = bath(&bank, 1.0, 400.0, 1e5);
    let mut state = equilibrium.clone();
    control.state.maximum_step_s = 0.75;
    let report = op
        .advance_network_to(
            &mut state,
            &mut time,
            clock + 0.75,
            &[internal()],
            &[boundary(&held)],
            &endpoints(&equilibrium),
            controls(),
            control,
            || false,
        )
        .unwrap();
    assert_eq!(time.to_bits(), (clock + 0.75).to_bits());
    assert_eq!(state, equilibrium);
    assert_eq!(report.minimum_accepted_substep_s, 0.25);
    assert_eq!(report.maximum_accepted_substep_s, 0.5);

    let target = 0.005_f64;
    time = 0.0;
    control.state.maximum_step_s = (target / 3.0).next_down();
    let report = op
        .advance_network_to(
            &mut state,
            &mut time,
            target,
            &[internal()],
            &[boundary(&held)],
            &endpoints(&equilibrium),
            controls(),
            control,
            || false,
        )
        .unwrap();
    assert_eq!(time.to_bits(), target.to_bits());
    assert_eq!(report.accepted_intervals, 3);
    assert_eq!(state, equilibrium);
}
