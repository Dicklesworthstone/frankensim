//! G0/G1/G3/G4/G5: reaction competition, Bateman source convolution, stable
//! limits, source/parent conservation, cancellation and checkpoint replay.
//! All data are explicitly synthetic and verify equations, not an isotope.

use super::*;
use crate::DecayHeatGroup;

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn budget() -> IrradiationBudget {
    IrradiationBudget {
        max_neutron_groups: 8,
        source: FissionBudget {
            max_decay_groups: 24,
            energy_tolerance_j: 1e-8,
            fuel_tolerance_nuclei: 1e-9,
        },
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (tolerance {tolerance:.4e})"
    );
}

fn fixture(
    k_f: f64,
    k_c: f64,
    lambda: f64,
    parent: f64,
    stored: f64,
) -> (IrradiationState, GroupFlux) {
    let groups =
        Arc::new(NeutronGroups::new("synthetic-one-group-joule".into(), vec![0.0, 1.0]).unwrap());
    let fission_model = Arc::new(
        FissionModel::new(
            "synthetic-fission-yield".into(),
            energy(3.0, 2.0, 1.0),
            vec![DecayHeatGroup {
                decay_constant_s_inv: lambda,
                energy_per_fission: energy(1.0, 0.5, 0.5),
            }],
        )
        .unwrap(),
    );
    let model = Arc::new(
        IrradiationModel::new(IrradiationData {
            data_id: "synthetic-microscopic-reactions".into(),
            nuclide: NuclideId::new(92, 235, 0).unwrap(),
            groups: Arc::clone(&groups),
            cross_sections: vec![GroupCrossSections {
                fission_m2: k_f,
                capture_m2: k_c,
            }],
            temperature_min_k: 250.0,
            temperature_max_k: 1000.0,
            fission_model,
            capture_energy_per_event: energy(1.0, 2.0, 0.0),
        })
        .unwrap(),
    );
    (
        IrradiationState::new(model, parent, vec![stored]).unwrap(),
        GroupFlux::new(groups, vec![1.0]).unwrap(),
    )
}

#[test]
fn g1_multigroup_rates_and_competing_capture_use_the_same_parent_depletion() {
    let (single, _) = fixture(0.0, 0.0, 0.2, 100.0, 0.0);
    let groups = Arc::new(
        NeutronGroups::new("synthetic-two-group-joule".into(), vec![0.0, 1e-19, 1e-12]).unwrap(),
    );
    let mut data = single.model.data.clone();
    data.groups = Arc::clone(&groups);
    data.cross_sections = vec![
        GroupCrossSections {
            fission_m2: 1e-22,
            capture_m2: 2e-22,
        },
        GroupCrossSections {
            fission_m2: 3e-22,
            capture_m2: 1e-22,
        },
    ];
    let mut state = IrradiationState::new(
        Arc::new(IrradiationModel::new(data).unwrap()),
        100.0,
        vec![0.0],
    )
    .unwrap();
    let flux = GroupFlux::new(groups, vec![2e20, 1e20]).unwrap();
    let constants = state
        .model()
        .rate_constants(&flux, 300.0, budget(), || Ok(()))
        .unwrap();
    close(constants.fission_s_inv, 0.05, 1e-17);
    close(constants.capture_s_inv, 0.05, 1e-17);
    let power = state.power_w(&flux, 300.0, budget(), || Ok(())).unwrap();
    close(power.deposited_j, 20.0, 1e-14);
    close(power.radiation_j, 20.0, 1e-14);
    close(power.escaped_j, 5.0, 1e-14);
    let report = state
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let expected_parent = 100.0 * (-0.2_f64).exp();
    let fissions = 50.0 * -(-0.2_f64).exp_m1();
    let stored = 10.0 * ((-0.2_f64).exp() - (-0.4_f64).exp()) / 0.1;
    let delayed_release = 2.0 * fissions - stored;
    close(state.remaining_parent_nuclei(), expected_parent, 2e-14);
    close(report.fissions, fissions, 3e-15);
    close(report.captures, fissions, 3e-15);
    close(state.stored_delayed_energy_j(), stored, 3e-14);
    close(
        report.released.deposited_j,
        4.0 * fissions + 0.5 * delayed_release,
        3e-14,
    );
    close(
        report.released.radiation_j,
        4.0 * fissions + 0.25 * delayed_release,
        3e-14,
    );
    close(
        report.released.escaped_j,
        fissions + 0.25 * delayed_release,
        1e-14,
    );
    close(report.fission_energy_j, 8.0 * fissions, 2e-14);
    close(report.capture_energy_j, 3.0 * fissions, 1e-14);
    close(
        report.generated_energy_j,
        report.released.total_j() + stored,
        5e-14,
    );
    let rates = state
        .reaction_rates(&flux, 300.0, budget(), || Ok(()))
        .unwrap();
    close(rates.fissions_per_s, expected_parent * 0.05, 1e-15);
    close(rates.captures_per_s, expected_parent * 0.05, 1e-15);
}

#[test]
fn g1_equal_depletion_and_decay_rates_use_the_bateman_repeated_root_limit() {
    let (mut state, flux) = fixture(0.125, 0.125, 0.25, 120.0, 7.0);
    let report = state
        .advance(3.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let survival = (-0.75_f64).exp();
    let stored = (7.0 + 2.0 * 0.125 * 120.0 * 3.0) * survival;
    close(state.stored_delayed_energy_j(), stored, 2e-14);
    close(state.remaining_parent_nuclei(), 120.0 * survival, 1e-14);
    close(report.fissions, 60.0 * (1.0 - survival), 1e-14);
    close(report.captures, report.fissions, 0.0);
    close(
        report.generated_energy_j + 7.0,
        report.released.total_j() + stored,
        6e-14,
    );
}

#[test]
fn g1_pure_capture_does_not_create_fissions_or_fission_decay_energy() {
    let (mut state, flux) = fixture(0.0, 0.3, 0.2, 100.0, 9.0);
    let report = state
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let captures = 100.0 * -(-0.6_f64).exp_m1();
    let old_release = 9.0 * -(-0.4_f64).exp_m1();
    assert_eq!(report.fissions, 0.0);
    assert_eq!(report.fission_energy_j, 0.0);
    close(report.captures, captures, 1e-14);
    close(
        state.stored_delayed_energy_j(),
        9.0 * (-0.4_f64).exp(),
        1e-15,
    );
    close(
        report.released.deposited_j,
        captures + old_release * 0.5,
        1e-14,
    );
    close(
        report.released.radiation_j,
        2.0 * captures + old_release * 0.25,
        2e-14,
    );
    close(report.released.escaped_j, old_release * 0.25, 1e-15);
}

#[test]
fn g1_pure_fission_and_shutdown_retain_the_irradiation_history() {
    let (mut state, flux) = fixture(0.15, 0.0, 0.2, 100.0, 0.0);
    let burn = state
        .advance(4.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(burn.captures, 0.0);
    close(burn.fissions, 100.0 * -(-0.6_f64).exp_m1(), 1e-14);
    let before = state.clone();
    let zero = GroupFlux::new(Arc::clone(flux.groups()), vec![0.0]).unwrap();
    let report = state
        .advance(5.0, &zero, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.fissions, 0.0);
    assert_eq!(report.captures, 0.0);
    assert_eq!(report.generated_energy_j, 0.0);
    assert_eq!(
        state.remaining_parent_nuclei(),
        before.remaining_parent_nuclei()
    );
    assert_eq!(state.cumulative_fissions(), before.cumulative_fissions());
    assert_eq!(state.cumulative_captures(), before.cumulative_captures());
    close(
        state.stored_delayed_energy_j(),
        before.stored_delayed_energy_j() * (-1.0_f64).exp(),
        2e-14,
    );
    close(
        report.released.total_j(),
        before.stored_delayed_energy_j() * (1.0 - (-1.0_f64).exp()),
        3e-14,
    );
    assert!(report.released.deposited_j > 0.0);
    let power = state.power_w(&zero, 300.0, budget(), || Ok(())).unwrap();
    assert_eq!(power, state.decay_power_w().unwrap());
}

#[test]
fn g3_partition_and_checkpoint_replay_preserve_competition_and_energy() {
    let (original, flux) = fixture(0.1, 0.2, 0.3, 200.0, 17.0);
    let (whole, _) = original
        .trial(8.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let mut split = original;
    for _ in 0..32 {
        split
            .advance(0.125, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
    }
    let saved = split.clone();
    for _ in 0..32 {
        split
            .advance(0.125, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
    }
    let mut replay = saved;
    for _ in 0..32 {
        replay
            .advance(0.125, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
    }
    assert_eq!(split, replay);
    close(
        split.remaining_parent_nuclei(),
        whole.remaining_parent_nuclei(),
        3e-14,
    );
    close(
        split.cumulative_fissions(),
        whole.cumulative_fissions(),
        3e-14,
    );
    close(
        split.cumulative_captures(),
        whole.cumulative_captures(),
        6e-14,
    );
    close(
        split.stored_delayed_energy_j(),
        whole.stored_delayed_energy_j(),
        5e-14,
    );
    close(
        split.cumulative_released().deposited_j,
        whole.cumulative_released().deposited_j,
        2e-13,
    );
    close(
        split.cumulative_released().radiation_j,
        whole.cumulative_released().radiation_j,
        2e-13,
    );
    close(
        split.cumulative_released().escaped_j,
        whole.cumulative_released().escaped_j,
        5e-14,
    );
}

#[test]
fn g3_nearly_coincident_rates_are_continuous_on_both_sides() {
    let mut stores = Vec::new();
    for lambda in [0.25 * (1.0 - 1e-12), 0.25, 0.25 * (1.0 + 1e-12)] {
        let (mut state, flux) = fixture(0.125, 0.125, lambda, 120.0, 7.0);
        state
            .advance(3.0, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
        stores.push(state.stored_delayed_energy_j());
    }
    assert!(stores[0] > stores[1] && stores[1] > stores[2]);
    close(stores[0], stores[1], 3e-11);
    close(stores[2], stores[1], 3e-11);
}

#[test]
fn g1_tiny_new_decay_heat_survives_both_small_and_large_parent_exposure() {
    for k in [1e-18, 0.1, 10.0] {
        let (base, _) = fixture(k, 0.0, 1e-18, 100.0, 0.0);
        let mut data = base.model.data.clone();
        data.fission_model = Arc::new(
            FissionModel::new(
                "synthetic-delayed-only".into(),
                EnergyChannels::default(),
                vec![DecayHeatGroup {
                    decay_constant_s_inv: 1e-18,
                    energy_per_fission: energy(2.0, 0.0, 0.0),
                }],
            )
            .unwrap(),
        );
        let flux = GroupFlux::new(Arc::clone(&data.groups), vec![1.0]).unwrap();
        let mut state = IrradiationState::new(
            Arc::new(IrradiationModel::new(data).unwrap()),
            100.0,
            vec![0.0],
        )
        .unwrap();
        let report = state
            .advance(1.0, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
        // First-order expansion in lambda only: integral_0^1
        // 2*N0*k*exp(-k*s)*lambda*(1-s) ds. Higher lambda terms <1e-35.
        let expected = if k < 1e-10 {
            100.0 * k * 1e-18
        } else {
            200.0 * 1e-18 * (1.0 - (1.0 - (-k).exp()) / k)
        };
        assert!(report.released.deposited_j > 0.0);
        close(report.released.deposited_j, expected, expected * 3e-14);
        assert_eq!(report.released.radiation_j, 0.0);
        assert_eq!(report.released.escaped_j, 0.0);
    }
}

#[test]
fn g1_stiff_depletion_retains_a_representable_parent_and_delayed_tail() {
    let (base, _) = fixture(1000.0, 0.0, 1000.0, 1e300, 0.0);
    let mut data = base.model.data.clone();
    data.fission_model = Arc::new(
        FissionModel::new(
            "synthetic-stiff-delayed-only".into(),
            EnergyChannels::default(),
            vec![DecayHeatGroup {
                decay_constant_s_inv: 1000.0,
                energy_per_fission: energy(1.0, 0.0, 0.0),
            }],
        )
        .unwrap(),
    );
    let flux = GroupFlux::new(Arc::clone(&data.groups), vec![1.0]).unwrap();
    let mut state = IrradiationState::new(
        Arc::new(IrradiationModel::new(data).unwrap()),
        1e300,
        vec![0.0],
    )
    .unwrap();
    let wide = IrradiationBudget {
        source: FissionBudget {
            energy_tolerance_j: 1e286,
            fuel_tolerance_nuclei: 1e286,
            ..budget().source
        },
        ..budget()
    };
    let report = state.advance(1.0, &flux, 300.0, wide, || Ok(())).unwrap();
    let parent = (1e300_f64.ln() - 1000.0).exp();
    let stored = (1e300_f64.ln() + 1000.0_f64.ln() - 1000.0).exp();
    assert!(state.remaining_parent_nuclei() > 0.0);
    assert!(state.stored_delayed_energy_j() > 0.0);
    close(state.remaining_parent_nuclei(), parent, parent * 7e-14);
    close(state.stored_delayed_energy_j(), stored, stored * 9e-14);
    close(report.fissions, 1e300, 0.0);
    close(report.generated_energy_j, 1e300, 0.0);
}

#[test]
fn g1_very_large_exposure_does_not_erase_slow_delayed_heat_by_ratio_underflow() {
    let (base, _) = fixture(1e308, 0.0, 1e-20, 1e20, 0.0);
    let mut data = base.model.data.clone();
    data.fission_model = Arc::new(
        FissionModel::new(
            "synthetic-impulsive-burn-slow-decay".into(),
            EnergyChannels::default(),
            vec![DecayHeatGroup {
                decay_constant_s_inv: 1e-20,
                energy_per_fission: energy(1.0, 0.0, 0.0),
            }],
        )
        .unwrap(),
    );
    let flux = GroupFlux::new(Arc::clone(&data.groups), vec![1.0]).unwrap();
    let mut state = IrradiationState::new(
        Arc::new(IrradiationModel::new(data).unwrap()),
        1e20,
        vec![0.0],
    )
    .unwrap();
    // The 1 J release is too small to subtract from a 1e20 J reservoir in
    // f64; that defect must be charged, while the positive heat is retained.
    let admitted = IrradiationBudget {
        source: FissionBudget {
            energy_tolerance_j: 10.0,
            ..budget().source
        },
        ..budget()
    };
    let report = state
        .advance(1.0, &flux, 300.0, admitted, || Ok(()))
        .unwrap();
    assert_eq!(report.fissions, 1e20);
    assert_eq!(state.remaining_parent_nuclei(), 0.0);
    close(report.released.deposited_j, 1.0, 1e-15);
    assert!(report.energy_roundoff_j >= 1.0);
}

#[test]
fn g0_repeated_sub_ulp_fissions_and_captures_are_not_lost() {
    let (base, _) = fixture(5e-17, 5e-17, 1.0, 1e16, 0.0);
    let mut data = base.model.data.clone();
    data.fission_model = Arc::new(
        FissionModel::new(
            "synthetic-prompt-only".into(),
            energy(1.0, 0.0, 0.0),
            vec![],
        )
        .unwrap(),
    );
    data.capture_energy_per_event = energy(1.0, 0.0, 0.0);
    let flux = GroupFlux::new(Arc::clone(&data.groups), vec![1.0]).unwrap();
    let mut state =
        IrradiationState::new(Arc::new(IrradiationModel::new(data).unwrap()), 1e16, vec![])
            .unwrap();
    for _ in 0..64 {
        state
            .advance(1.0, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
    }
    assert_eq!(state.remaining_parent_nuclei(), 1e16 - 64.0);
    close(state.cumulative_fissions(), 32.0, 2e-13);
    close(state.cumulative_captures(), 32.0, 2e-13);
    close(state.cumulative_released().deposited_j, 64.0, 4e-13);
}

#[test]
fn g0_group_data_temperature_and_budget_admission_refuse_ambiguous_inputs() {
    assert!(NuclideId::new(0, 1, 0).is_err());
    assert!(NuclideId::new(92, 91, 0).is_err());
    assert_ne!(
        NuclideId::new(92, 235, 0).unwrap(),
        NuclideId::new(92, 235, 1).unwrap()
    );
    for edges in [
        vec![],
        vec![0.0],
        vec![1.0, 0.0],
        vec![0.0, 0.0],
        vec![0.0, f64::NAN],
    ] {
        assert!(NeutronGroups::new("bad".into(), edges).is_err());
    }
    let (state, flux) = fixture(0.1, 0.2, 0.3, 10.0, 0.0);
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(GroupFlux::new(Arc::clone(flux.groups()), vec![bad]).is_err());
        let mut data = state.model.data.clone();
        data.cross_sections[0].fission_m2 = bad;
        assert!(IrradiationModel::new(data).is_err());
    }
    assert!(GroupFlux::new(Arc::clone(flux.groups()), vec![]).is_err());
    for (name, edges) in [
        ("different-name", vec![0.0, 1.0]),
        ("synthetic-one-group-joule", vec![0.0, 2.0]),
    ] {
        let wrong = GroupFlux::new(
            Arc::new(NeutronGroups::new(name.into(), edges).unwrap()),
            vec![0.0],
        )
        .unwrap();
        assert!(
            state
                .trial(0.0, &wrong, 300.0, budget(), || Ok(()))
                .is_err()
        );
    }
    for bad_temperature in [249.0, 1001.0, f64::NAN] {
        assert!(
            state
                .trial(0.0, &flux, bad_temperature, budget(), || Ok(()))
                .is_err()
        );
    }
    assert!(
        state
            .trial(
                1.0,
                &flux,
                300.0,
                IrradiationBudget {
                    max_neutron_groups: 0,
                    ..budget()
                },
                || Ok(())
            )
            .is_err()
    );
    assert!(
        state
            .trial(
                1.0,
                &flux,
                300.0,
                IrradiationBudget {
                    source: FissionBudget {
                        max_decay_groups: 0,
                        ..budget().source
                    },
                    ..budget()
                },
                || Ok(())
            )
            .is_err()
    );
    let (same, empty) = state.trial(0.0, &flux, 300.0, budget(), || Ok(())).unwrap();
    assert_eq!(same, state);
    assert_eq!(empty, IrradiationStep::default());
}

#[test]
fn g4_late_cancellation_or_accounting_refusal_publishes_no_irradiation() {
    let (mut state, flux) = fixture(0.1, 0.2, 0.3, 100.0, 7.0);
    let original = state.clone();
    let mut polls = 0;
    state
        .trial(1.0, &flux, 300.0, budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut current = 0;
    let result = state.advance(1.0, &flux, 300.0, budget(), || {
        current += 1;
        if current == polls {
            Err(FissionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(result, Err(FissionError::Cancelled));
    assert_eq!(state, original);
    for duration in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            state
                .advance(duration, &flux, 300.0, budget(), || Ok(()))
                .is_err()
        );
        assert_eq!(state, original);
    }
    let (mut stored, flux) = fixture(0.0, 0.0, 1e-20, 0.0, 1e20);
    let original = stored.clone();
    let strict = IrradiationBudget {
        source: FissionBudget {
            energy_tolerance_j: 0.1,
            ..budget().source
        },
        ..budget()
    };
    assert_eq!(
        stored.advance(1.0, &flux, 300.0, strict, || Ok(())),
        Err(FissionError::Budget("energy accounting roundoff"))
    );
    assert_eq!(stored, original);
}
