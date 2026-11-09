use super::*;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn model() -> SingletTripletDecay {
    SingletTripletDecay::try_new([4e-19, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap()
}

fn limits() -> PumpColumnLimits {
    PumpColumnLimits {
        relative_tolerance: 2e-5,
        maximum_cell_evaluations: 500_000,
        ..PumpColumnLimits::default()
    }
}

#[test]
fn g1_population_residence_retains_the_real_stable_ground_state() {
    let model =
        SingletTripletDecay::try_new([4e-19, 3e-19], [2.0, 0.0], [0.0; 2], [0.0; 2]).unwrap();
    let plan = SaturableSingletTripletPlan::prepare_with_population_integrals(
        model,
        [0.0; 2],
        [0.0; 2],
        0.5,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let residence = plan
        .integrate_populations_with_cancel([2.0, 3.0, 4.0], || false)
        .unwrap();
    let singlet = 1.5 * (1.0 - (-1.0_f64).exp());
    close(residence[0], 2.5 - singlet, 2e-13);
    close(residence[1], singlet, 2e-13);
    close(residence[2], 2.0, 2e-13);
    close(residence.iter().sum(), 4.5, 2e-13);
}

#[test]
fn g1_passive_host_absorption_keeps_independent_spectral_heat_and_transmission() {
    let depths = [[0.2, 0.7], [0.4, 0.1]];
    let cells = depths.map(|passive_optical_depth| SaturablePumpCell {
        model: model(),
        optical_depth_per_ground: [[0.0; 2]; 2],
        passive_optical_depth,
    });
    let photons = [1e19, 2e19];
    let energy = [5e-19, 4e-19];
    let column = SaturablePumpColumn::try_new(&cells, energy, 2, || false).unwrap();
    let step = column
        .advance(&[[0.0; 3]; 2], photons, 0.3, limits(), || false)
        .unwrap();
    assert_eq!(step.substeps, 1);
    for lane in 0..2 {
        let incoming = photons[lane] * 0.3;
        let outgoing = incoming * (-(depths[0][lane] + depths[1][lane])).exp();
        close(step.transmitted_photons[lane], outgoing, 2e-13);
        close(
            step.passive_absorbed_photons[lane],
            incoming - outgoing,
            2e-13,
        );
        close(
            step.passive_absorbed_energy_j[lane],
            (incoming - outgoing) * energy[lane],
            2e-13,
        );
        assert_eq!(step.absorbed_photons[lane], 0.0);
    }
    assert!(step.residual_j.abs() < 1e-12);
    assert!(step.cells.iter().all(|cell| cell.passive_heat_j > 0.0));
}

#[test]
fn g1_one_shot_saturable_slab_matches_independent_fluence_solution() {
    // No decay, no stimulated emission: saturation fluence is one photon per
    // cross section. exp(a*G)-1 decays as exp(-a*incident_fluence).
    let count = 1e19;
    let depth_per_ground = 2.0 / count;
    let column = SaturablePumpColumn::try_new(
        &[SaturablePumpCell {
            model: model(),
            optical_depth_per_ground: [[depth_per_ground, 0.0]],
            passive_optical_depth: [0.0],
        }],
        [5e-19],
        1,
        || false,
    )
    .unwrap();
    let incident = 5.0 / depth_per_ground;
    let step = column
        .advance(&[[count, 0.0, 0.0]], [incident], 1.0, limits(), || false)
        .unwrap();
    // Independent 80-digit evaluation of ln(1+exp(-2)*(exp(5)-1)).
    let outgoing_fluence = 3.042_148_271_140_494_6;
    close(
        step.transmitted_photons[0] * depth_per_ground,
        outgoing_fluence,
        2e-5,
    );
    close(
        step.absorbed_photons[0] * depth_per_ground,
        5.0 - outgoing_fluence,
        2e-5,
    );
    close(
        step.cells[0].kinetics.populations[0] * depth_per_ground,
        outgoing_fluence - 3.0,
        2e-5,
    );
    assert!(step.absorbed_photons[0] < count);
    close(
        step.cells[0].kinetics.populations.iter().sum(),
        count,
        2e-12,
    );
    assert!(step.substeps > 2);
    assert!(step.residual_j.abs() < 1e-10);
}

#[test]
fn g3_bleaching_front_uses_upstream_transmission_and_native_time_refinement() {
    let count = 1e19;
    let coefficient = 1.0 / count;
    let cells = [SaturablePumpCell {
        model: model(),
        optical_depth_per_ground: [[coefficient, 0.0]],
        passive_optical_depth: [0.0],
    }; 2];
    let column = SaturablePumpColumn::try_new(&cells, [5e-19], 2, || false).unwrap();
    let initial = [[count, 0.0, 0.0]; 2];
    let rate = [5.0 / coefficient];
    let whole = column
        .advance(&initial, rate, 1.0, limits(), || false)
        .unwrap();
    let mid = column
        .advance(&initial, rate, 0.5, limits(), || false)
        .unwrap();
    let populations: Vec<_> = mid
        .cells
        .iter()
        .map(|cell| cell.kinetics.populations)
        .collect();
    let end = column
        .advance(&populations, rate, 0.5, limits(), || false)
        .unwrap();
    let front_fluence = (1.0 + (-1.0_f64).exp() * 5.0_f64.exp_m1()).ln();
    close(
        whole.cells[0].transmitted_photons[0] * coefficient,
        front_fluence,
        3e-5,
    );
    close(
        whole.transmitted_photons[0] * coefficient,
        3.042_148_271_140_494_6,
        3e-5,
    );
    assert!(whole.cells[0].kinetics.populations[0] < whole.cells[1].kinetics.populations[0]);
    close(
        whole.transmitted_photons[0],
        mid.transmitted_photons[0] + end.transmitted_photons[0],
        5e-5,
    );
    for index in 0..2 {
        close(
            whole.cells[index].kinetics.populations[0],
            end.cells[index].kinetics.populations[0],
            5e-5,
        );
        close(
            whole.cells[index].absorbed_photons[0],
            mid.cells[index].absorbed_photons[0] + end.cells[index].absorbed_photons[0],
            5e-5,
        );
    }
}

#[test]
fn g3_spectral_competition_keeps_target_events_separate_from_photon_energies() {
    let model =
        SingletTripletDecay::try_new([4e-19, 3e-19], [3.0, 0.2], [0.7, 0.1], [1.0, 0.05]).unwrap();
    let cell = SaturablePumpCell {
        model,
        optical_depth_per_ground: [[0.7e-19, 0.1e-19], [0.2e-19, 0.8e-19]],
        passive_optical_depth: [0.2, 0.4],
    };
    let column = SaturablePumpColumn::try_new(&[cell], [5e-19, 4e-19], 1, || false).unwrap();
    let step = column
        .advance(&[[1e19, 0.0, 0.0]], [2e19, 1e19], 0.4, limits(), || false)
        .unwrap();
    let k = &step.cells[0].kinetics;
    close(k.populations.iter().sum(), 1e19, 2e-12);
    close(
        k.pump_events.iter().sum(),
        step.absorbed_photons.iter().sum(),
        2e-12,
    );
    close(
        k.absorbed_pump_energy_j.iter().sum(),
        step.absorbed_energy_j.iter().sum(),
        2e-12,
    );
    assert!(k.pump_events.into_iter().all(|value| value > 0.0));
    assert!(k.crossing_events.into_iter().all(|value| value > 0.0));
    for lane in 0..2 {
        close(
            step.incident_photons[lane],
            step.transmitted_photons[lane]
                + step.absorbed_photons[lane]
                + step.passive_absorbed_photons[lane],
            2e-12,
        );
    }
    assert!(step.residual_j.abs() < 1e-10);
}

#[test]
fn g0_column_admission_work_and_cancellation_do_not_modify_populations() {
    let cell = SaturablePumpCell {
        model: model(),
        optical_depth_per_ground: [[0.2, 0.0]],
        passive_optical_depth: [0.0],
    };
    assert!(SaturablePumpColumn::try_new(&[cell], [3e-19], 1, || false).is_err());
    assert!(SaturablePumpColumn::try_new(&[cell], [5e-19], 0, || false).is_err());
    let column = SaturablePumpColumn::try_new(&[cell], [5e-19], 1, || false).unwrap();
    let populations = [[10.0, 0.0, 0.0]];
    let before = populations;
    let zero = column
        .advance(&populations, [25.0], 0.0, limits(), || false)
        .unwrap();
    assert_eq!(zero.cells[0].kinetics.populations, populations[0]);
    assert_eq!(zero.incident_photons, [0.0]);
    for limited in [
        PumpColumnLimits {
            maximum_root_iterations: 2,
            ..limits()
        },
        PumpColumnLimits {
            maximum_substeps: 1,
            ..limits()
        },
        PumpColumnLimits {
            maximum_cell_evaluations: 2,
            ..limits()
        },
    ] {
        assert!(matches!(
            column.advance(&populations, [25.0], 1.0, limited, || false),
            Err(SingletTripletError::Evolution(DecayError::Budget { .. }))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        column
            .advance(&populations, [25.0], 1.0, limits(), || {
                polls += 1;
                polls == 25
            })
            .unwrap_err(),
        SingletTripletError::Evolution(DecayError::Cancelled)
    );
    assert_eq!(populations, before);
}

#[test]
fn g3_a_weak_high_energy_pump_keeps_its_heat_below_the_mean_photon_energy_ulp() {
    let column = SaturablePumpColumn::try_new(
        &[SaturablePumpCell {
            model: model(),
            optical_depth_per_ground: [[1e-19, 0.0]; 2],
            passive_optical_depth: [0.0; 2],
        }],
        [4e-19, 5e-19],
        1,
        || false,
    )
    .unwrap();
    let step = column
        .advance(&[[1e19, 0.0, 0.0]], [1e19, 10.0], 0.5, limits(), || false)
        .unwrap();
    let heat = step.cells[0].kinetics.pump_relaxation_heat_j[0];
    assert!(heat > 0.0);
    close(heat, step.absorbed_photons[1] * (5e-19 - 4e-19), 2e-12);
    // Averaging these photon energies first rounds to the resonant value.
    assert_eq!(4e-19 + 1e-18 * (5e-19 - 4e-19), 4e-19);
}
