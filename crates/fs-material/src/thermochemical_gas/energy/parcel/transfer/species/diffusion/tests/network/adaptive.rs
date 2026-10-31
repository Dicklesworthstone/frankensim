use super::*;
use network::adaptive::{
    DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceError, advance_candidate,
};

fn accuracy(relative: f64) -> DiffusionNetworkAdvanceControls {
    DiffusionNetworkAdvanceControls {
        maximum_step_s: 25.0,
        minimum_step_s: 1e-6,
        maximum_evaluations: 10000,
        absolute_species_tolerance_kg: 1e-13,
        relative_species_tolerance: relative,
        temperature_tolerance_k: 1e-6,
    }
}
fn composition(parcel: &ClosedGasParcel) -> f64 {
    parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[0]
}
fn species_mass(parcel: &ClosedGasParcel, name: &str) -> f64 {
    let mixture = parcel.model().thermodynamics();
    let amount = parcel.mass_kg()
        / parcel
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
    mixture
        .components()
        .iter()
        .zip(mixture.mole_composition().fractions())
        .find(|(model, _)| model.species().as_str() == name)
        .map_or(0.0, |(model, x)| amount * x * model.molar_mass().value())
}

#[test]
fn g1_adaptive_molecular_window_reaches_requested_time_and_refines_to_binary_exponential() {
    let d = card(2e-5);
    let point = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &point)];
    let links = [link(7, [0, 1], &pairs)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.25))]);
    let bounds = endpoints(&initial);
    let saved = initial.clone();
    assert!(matches!(
        network::step_candidate(&initial, &links, 25.0, &bounds, solve(), || false),
        Err(DiffusionNetworkError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    let mut errors = Vec::new();
    for tolerance in [4e-4, 2.5e-5] {
        let next = advance_candidate(
            &initial,
            &links,
            25.0,
            &bounds,
            solve(),
            accuracy(tolerance),
            || false,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.duration_s, 25.0);
        assert!(report.accepted_intervals > 1 && report.rejected_intervals > 0);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.evaluations > 3 * report.accepted_intervals);
        assert!(
            report.maximum_species_error_ratio <= 1.0
                && report.maximum_temperature_error_ratio <= 1.0
        );
        assert!(report.total_mass_residual_kg.abs() <= solve().exchange.mass_tolerance_kg);
        assert!(report.maximum_species_mass_residual_kg <= solve().exchange.mass_tolerance_kg);
        assert!(report.total_energy_residual_j.abs() <= solve().exchange.energy_tolerance_j);
        assert!(report.entropy_change_j_k > 0.0);
        let gap = composition(&next.parcels()[&0]) - composition(&next.parcels()[&1]);
        errors.push((gap - 0.5 * (-1.0_f64).exp()).abs());
        let moved = species_mass(&initial[&0], "A") - species_mass(&next.parcels()[&0], "A");
        let account = report
            .transfers
            .iter()
            .find(|row| row.donor == 0 && row.receiver == 1 && row.species.as_str() == "A")
            .unwrap();
        near(account.mass_kg, moved, 1e-14);
        // The initial 25 s rejected probe moves much more than the accepted
        // relaxation. Its nominal flux must not appear in the physical account.
        assert!(account.mass_kg < 2e-4);
    }
    assert!(
        errors[0] < 0.003 && errors[1] < 0.32 * errors[0],
        "{errors:?}"
    );
    assert_eq!(initial, saved);
}

#[test]
fn g1_g3_adaptive_branch_mode_preserves_simultaneous_transfer_and_link_order() {
    let d = card(2e-5);
    let point = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &point)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.25)), (2, gas(0.25))]);
    let bounds = endpoints(&initial);
    let first = advance_candidate(
        &initial,
        &[link(3, [0, 1], &pairs), link(8, [0, 2], &pairs)],
        10.0,
        &bounds,
        solve(),
        accuracy(1e-4),
        || false,
    )
    .unwrap();
    let reordered = advance_candidate(
        &initial,
        &[link(999, [0, 2], &pairs), link(2, [0, 1], &pairs)],
        10.0,
        &bounds,
        solve(),
        accuracy(1e-4),
        || false,
    )
    .unwrap();
    assert_eq!(first.parcels(), reordered.parcels());
    assert_eq!(first.report().transfers, reordered.report().transfers);
    assert_eq!(first.report().evaluations, reordered.report().evaluations);
    assert_eq!(first.parcels()[&1], first.parcels()[&2]);
    let exact_center = 5.0 / 12.0 + (-0.6_f64).exp() / 3.0;
    let exact_outer = 5.0 / 12.0 - (-0.6_f64).exp() / 6.0;
    near(composition(&first.parcels()[&0]), exact_center, 0.001);
    near(composition(&first.parcels()[&1]), exact_outer, 0.001);
}

fn formation_gas(x: f64, offsets: [f64; 2]) -> ClosedGasParcel {
    use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};
    let original = gas(x);
    let components = original
        .model()
        .thermodynamics()
        .components()
        .iter()
        .enumerate()
        .map(|(i, model)| {
            let regions = model
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a7").unwrap().value += offsets[i];
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                model.species().clone(),
                model.molar_mass(),
                model.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    let mixture = FrozenIdealGasMixtureModelV1::new(
        components,
        Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
    )
    .unwrap();
    let model = ThermochemicalGasModel::new(
        mixture,
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(400.0, 100000.0), 0.001).unwrap()
}

#[test]
fn g3_adaptive_accuracy_and_energy_accounts_ignore_species_formation_datums() {
    let d = card(2e-5);
    let point = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &point)];
    let links = [link(1, [0, 1], &pairs)];
    let ordinary = BTreeMap::from([
        (0, formation_gas(0.75, [0.0; 2])),
        (1, formation_gas(0.25, [0.0; 2])),
    ]);
    let shifted = BTreeMap::from([
        (0, formation_gas(0.75, [1e18, -3e18])),
        (1, formation_gas(0.25, [1e18, -3e18])),
    ]);
    let run = |states: &BTreeMap<u64, ClosedGasParcel>| {
        advance_candidate(
            states,
            &links,
            10.0,
            &endpoints(states),
            solve(),
            accuracy(1e-4),
            || false,
        )
        .unwrap()
    };
    let a = run(&ordinary);
    let b = run(&shifted);
    for id in ordinary.keys() {
        assert_eq!(
            a.parcels()[id].mass_kg().to_bits(),
            b.parcels()[id].mass_kg().to_bits()
        );
        assert_eq!(
            a.parcels()[id].state().state(),
            b.parcels()[id].state().state()
        );
        assert_eq!(
            composition(&a.parcels()[id]).to_bits(),
            composition(&b.parcels()[id]).to_bits()
        );
    }
    assert_eq!(a.report().evaluations, b.report().evaluations);
    assert_eq!(a.report().accepted_intervals, b.report().accepted_intervals);
    assert_eq!(
        a.report().maximum_species_error_ratio.to_bits(),
        b.report().maximum_species_error_ratio.to_bits()
    );
    assert_eq!(
        a.report().total_energy_residual_j.to_bits(),
        b.report().total_energy_residual_j.to_bits()
    );
    assert_ne!(
        a.report().transfers[0].carried_enthalpy_j,
        b.report().transfers[0].carried_enthalpy_j
    );
}

#[test]
fn g0_g4_adaptive_window_refuses_accuracy_work_missing_sources_and_late_cancel_atomically() {
    let d = card(2e-5);
    let point = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &point)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.25))]);
    let saved = initial.clone();
    let bounds = endpoints(&initial);
    let links = [link(1, [0, 1], &pairs)];
    let strict = DiffusionNetworkAdvanceControls {
        minimum_step_s: 1.0,
        maximum_step_s: 1.0,
        absolute_species_tolerance_kg: 1e-20,
        relative_species_tolerance: 0.0,
        ..accuracy(1e-4)
    };
    assert!(matches!(
        advance_candidate(&initial, &links, 1.0, &bounds, solve(), strict, || false),
        Err(DiffusionNetworkAdvanceError::Accuracy { .. })
    ));
    let one_attempt = DiffusionNetworkAdvanceControls {
        maximum_evaluations: 3,
        ..accuracy(1e-4)
    };
    assert!(matches!(
        advance_candidate(
            &initial,
            &links,
            10.0,
            &bounds,
            solve(),
            one_attempt,
            || false
        ),
        Err(DiffusionNetworkAdvanceError::WorkLimit { evaluations: 3, .. })
    ));
    assert!(matches!(
        advance_candidate(
            &initial,
            &links,
            10.0,
            &bounds,
            solve(),
            DiffusionNetworkAdvanceControls {
                maximum_step_s: 0.5,
                maximum_evaluations: 4,
                ..accuracy(1e-4)
            },
            || false
        ),
        Err(DiffusionNetworkAdvanceError::WorkLimit {
            elapsed_s: 0.5,
            accepted_intervals: 1,
            evaluations: 4,
        })
    ));
    let missing = [link(1, [0, 1], &[])];
    assert!(matches!(
        advance_candidate(
            &initial,
            &missing,
            10.0,
            &bounds,
            solve(),
            accuracy(1e-4),
            || false
        ),
        Err(DiffusionNetworkAdvanceError::Step {
            error: DiffusionNetworkError::Face { .. },
            ..
        })
    ));
    // The caller context is inside the source's range, but the actual face is
    // 400 K. A smaller timestep must not manufacture material applicability.
    let mut claim = d.claims_for(BINARY_DIFFUSIVITY_PROPERTY)[0].1.clone();
    claim.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        200.0,
        300.0,
    );
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    let bounded = MaterialCard::assemble(d.id().clone(), claims, vec![]).unwrap();
    let stale = super::point(250.0, 100000.0);
    let unsupported = [pair("A", "B", &bounded, &stale)];
    assert!(matches!(
        advance_candidate(
            &initial,
            &[link(1, [0, 1], &unsupported)],
            10.0,
            &bounds,
            solve(),
            one_attempt,
            || false
        ),
        Err(DiffusionNetworkAdvanceError::Step {
            elapsed_s: 0.0,
            step_s: 10.0,
            error: DiffusionNetworkError::Face { .. },
        })
    ));
    let mut polls = 0;
    advance_candidate(
        &initial,
        &links,
        10.0,
        &bounds,
        solve(),
        accuracy(1e-4),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut called = 0;
    assert_eq!(
        advance_candidate(
            &initial,
            &links,
            10.0,
            &bounds,
            solve(),
            accuracy(1e-4),
            || {
                called += 1;
                called == polls - 1
            }
        ),
        Err(DiffusionNetworkAdvanceError::Cancelled)
    );
    assert_eq!(initial, saved);
}

#[test]
fn g0_zero_transfer_windows_keep_states_and_finish_the_exact_short_remainder() {
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.25))]);
    let closed = DiffusionNetworkLink {
        geometry: GasDiffusionGeometry {
            area_m2: 0.0,
            ..geometry()
        },
        ..link(1, [0, 1], &[])
    };
    let controls = DiffusionNetworkAdvanceControls {
        minimum_step_s: 1.5,
        maximum_step_s: 2.0,
        ..accuracy(1e-4)
    };
    let next = advance_candidate(
        &initial,
        &[closed],
        5.0,
        &endpoints(&initial),
        solve(),
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(next.parcels(), &initial);
    assert_eq!(next.report().duration_s, 5.0);
    assert_eq!(next.report().accepted_intervals, 3);
    assert_eq!(next.report().accepted_substeps, 6);
    assert_eq!(next.report().evaluations, 9);
    assert_eq!(next.report().minimum_accepted_substep_s, 0.5);
    assert!(next.report().transfers.is_empty());
    assert_eq!(next.report().entropy_change_j_k, 0.0);
    assert_eq!(next.report().total_energy_residual_j, 0.0);
}
