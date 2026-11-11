mod adaptive;

use super::*;
use network::binary::{
    BinaryDiffusionNetworkControls, BinaryDiffusionNetworkError, advance_candidate,
};

fn resolution(substeps: usize, edges: usize) -> BinaryDiffusionNetworkControls {
    BinaryDiffusionNetworkControls {
        substeps,
        maximum_pair_evaluations: substeps * (2 * edges - 1),
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
    let moles = parcel.mass_kg()
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
        .map_or(0.0, |(model, x)| moles * x * model.molar_mass().value())
}

// Independent NASA fixture integral: Cp_A/R=3.5+0.001*T, Cp_B/R=2.5.
fn internal_energy(parcel: &ClosedGasParcel) -> f64 {
    let x = composition(parcel);
    let t = parcel.state().state().temperature;
    let moles = parcel.mass_kg()
        / parcel
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * moles
        * (x * (2.5 * t + 0.0005 * t * t) + (1.0 - x) * 1.5 * t)
}

#[test]
fn g1_one_edge_reaches_unequal_volume_equilibrium_beyond_explicit_stability() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let right = ClosedGasParcel::new(gas(0.25).model().clone(), context.clone(), 0.004).unwrap();
    let initial = BTreeMap::from([(0, gas(0.75)), (1, right)]);
    let bounds = endpoints(&initial);
    let links = [link(7, [1, 0], &pairs)];
    let mut physical = solve();
    physical.exchange.maximum_species_fraction = 0.9;
    assert!(network::step_candidate(&initial, &links, 1e8, &bounds, physical, || false,).is_err());
    let next = advance_candidate(
        &initial,
        &links,
        1e8,
        &bounds,
        physical,
        resolution(1, 1),
        || false,
    )
    .unwrap();
    let mut pair_controls = controls();
    pair_controls.exchange.maximum_species_fraction = 0.9;
    let direct = initial[&0]
        .diffuse_binary_exponential_with(
            &initial[&1],
            geometry(),
            1e8,
            &pairs,
            pair_controls,
            || false,
        )
        .unwrap();
    assert_eq!(&next.parcels()[&0], direct.exchange().left());
    assert_eq!(&next.parcels()[&1], direct.exchange().right());
    // Exact volume-weighted equilibrium, independent of molar masses.
    near(composition(&next.parcels()[&0]), 0.35, 3e-13);
    near(composition(&next.parcels()[&1]), 0.35, 3e-13);
    let report = next.report();
    assert_eq!(report.duration_s, 1e8);
    assert_eq!(report.pair_evaluations, 1);
    assert_eq!(report.last_sweep.len(), 1);
    assert_eq!(report.last_sweep[0].parcels, [0, 1]);
    assert_eq!(report.last_sweep[0].relaxation.relaxation_fraction, 1.0);
    assert!(report.entropy_change_j_k > 0.0);
    assert!(report.endpoint_energy_change_j.abs() <= physical.exchange.energy_tolerance_j);
    let moved = species_mass(&initial[&0], "A") - species_mass(&next.parcels()[&0], "A");
    let accounted = report
        .transfers
        .iter()
        .find(|row| row.donor == 0 && row.species.as_str() == "A")
        .unwrap();
    near(accounted.mass_kg, moved, 1e-14);
}

#[test]
fn g1_symmetric_chain_refines_and_accounts_for_every_finite_inventory() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let links = [link(9, [1, 2], &pairs), link(2, [0, 1], &pairs)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.5)), (2, gas(0.25))]);
    let bounds = endpoints(&initial);
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let next = advance_candidate(
            &initial,
            &links,
            10.0,
            &bounds,
            solve(),
            resolution(count, 2),
            || false,
        )
        .unwrap();
        let exact = [
            0.5 + 0.25 * (-0.2_f64).exp(),
            0.5,
            0.5 - 0.25 * (-0.2_f64).exp(),
        ];
        errors.push(
            (0..3)
                .map(|i| (composition(&next.parcels()[&(i as u64)]) - exact[i]).abs())
                .fold(0.0_f64, f64::max),
        );
        let report = next.report();
        assert_eq!(report.pair_evaluations, 3 * count);
        assert_eq!(report.last_sweep.len(), 3);
        assert_eq!(
            report
                .last_sweep
                .iter()
                .map(|stage| stage.parcels)
                .collect::<Vec<_>>(),
            vec![[0, 1], [1, 2], [0, 1]],
        );
        for (stage, fraction) in report.last_sweep.iter().zip([0.5, 1.0, 0.5]) {
            near(
                stage.relaxation.duration_s,
                fraction * 10.0 / count as f64,
                1e-14,
            );
        }
        for name in ["A", "B"] {
            let before: f64 = initial.values().map(|p| species_mass(p, name)).sum();
            let after: f64 = next.parcels().values().map(|p| species_mass(p, name)).sum();
            near(after, before, solve().exchange.mass_tolerance_kg);
            for &id in initial.keys() {
                let net: f64 = report
                    .transfers
                    .iter()
                    .filter(|row| row.species.as_str() == name)
                    .map(|row| {
                        if row.receiver == id {
                            row.mass_kg
                        } else if row.donor == id {
                            -row.mass_kg
                        } else {
                            0.0
                        }
                    })
                    .sum();
                near(
                    species_mass(&next.parcels()[&id], name),
                    species_mass(&initial[&id], name) + net,
                    1e-14,
                );
            }
        }
        assert!(report.total_mass_residual_kg.abs() <= solve().exchange.mass_tolerance_kg);
        assert!(report.maximum_species_mass_residual_kg <= solve().exchange.mass_tolerance_kg);
        assert!(report.total_energy_residual_j.abs() <= solve().exchange.energy_tolerance_j);
        assert!(report.endpoint_energy_change_j.abs() <= solve().exchange.energy_tolerance_j);
        assert!(report.entropy_change_j_k > 0.0);
        near(
            next.parcels().values().map(internal_energy).sum(),
            initial.values().map(internal_energy).sum(),
            solve().exchange.energy_tolerance_j,
        );
    }
    assert!(
        errors
            .windows(2)
            .all(|e| e[0] / e[1] > 3.8 && e[0] / e[1] < 4.2),
        "{errors:?}"
    );
}

fn temperature_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 1e-5), (1800.0, 9e-5)],
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            // The face query sends an absolute temperature. The claim declares
            // the same quantity over its knots, as the sibling card() does.
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                200.0,
                1800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic temperature-dependent pair".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "binary pair".into(),
            phase: "gas".into(),
            process: "source refresh".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn g3_cycle_rechecks_actual_sources_and_ignores_link_order_ids_and_orientation() {
    let d = temperature_card();
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let initial: BTreeMap<_, _> = [(0, 0.75, 390.0), (1, 0.5, 400.0), (2, 0.25, 410.0)]
        .into_iter()
        .map(|(id, x, t)| {
            (
                id,
                ClosedGasParcel::new(gas(x).model().clone(), point(t, 100000.0), 0.001).unwrap(),
            )
        })
        .collect();
    let bounds = endpoints(&initial);
    let make = |id, endpoints| {
        let mut edge = link(id, endpoints, &pairs);
        edge.solve.maximum_temperature_difference_k = 50.0;
        edge.solve.maximum_pressure_difference_pa = 50_000.0;
        edge
    };
    let links = [make(3, [0, 1]), make(8, [1, 2]), make(2, [0, 2])];
    let first = advance_candidate(
        &initial,
        &links,
        1.0,
        &bounds,
        solve(),
        resolution(2, 3),
        || false,
    )
    .unwrap();
    let reordered = [make(999, [2, 0]), make(1, [2, 1]), make(22, [1, 0])];
    let second = advance_candidate(
        &initial,
        &reordered,
        1.0,
        &bounds,
        solve(),
        resolution(2, 3),
        || false,
    )
    .unwrap();
    assert_eq!(first.parcels(), second.parcels());
    assert_eq!(first.report().transfers, second.report().transfers);
    assert_eq!(first.report().pair_evaluations, 10);
    assert_eq!(first.report().last_sweep.len(), 5);
    let mut refreshed = false;
    for stage in &first.report().last_sweep {
        let flux = &stage.initial_flux;
        let sample = &flux.pairs[0];
        near(
            sample.diffusivity_m2_s,
            5e-8 * flux.face_temperature_k,
            1e-17,
        );
        assert!(
            sample
                .material
                .query_point()
                .contains(&("temperature".into(), flux.face_temperature_k),)
        );
        let initial_mean = 0.5 * initial[&stage.parcels[0]].state().state().temperature
            + 0.5 * initial[&stage.parcels[1]].state().state().temperature;
        refreshed |= (flux.face_temperature_k - initial_mean).abs() > 1e-6;
    }
    assert!(refreshed);
    // Independent exponential amounts and quadratic NASA internal-energy
    // inversion, with D(T)=5e-8*T at every actual split stage.
    for (id, x, t) in [
        (0, 0.7358020822132906, 390.23034486598084),
        (1, 0.5000093712636019, 399.96902617255586),
        (2, 0.2649164105895538, 409.48818128683985),
    ] {
        near(composition(&first.parcels()[&id]), x, 3e-12);
        near(first.parcels()[&id].state().state().temperature, t, 2e-7);
    }
    near(
        first.parcels().values().map(internal_energy).sum(),
        initial.values().map(internal_energy).sum(),
        solve().exchange.energy_tolerance_j,
    );
}

#[test]
fn g4_work_source_structure_and_late_cancel_refuse_without_partial_publication() {
    let d = card(2e-5);
    let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &context)];
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.5)), (2, gas(0.25))]);
    let saved = initial.clone();
    let bounds = endpoints(&initial);
    let links = [link(3, [0, 1], &pairs), link(8, [1, 2], &pairs)];
    let budget = BinaryDiffusionNetworkControls {
        maximum_pair_evaluations: 5,
        ..resolution(2, 2)
    };
    assert!(matches!(
        advance_candidate(&initial, &links, 0.1, &bounds, solve(), budget, || false,),
        Err(BinaryDiffusionNetworkError::WorkLimit {
            required_pair_evaluations: 6,
            maximum_pair_evaluations: 5,
        })
    ));
    let missing = [link(3, [0, 1], &pairs), link(8, [1, 2], &[])];
    assert!(matches!(
        advance_candidate(
            &initial,
            &missing,
            0.1,
            &bounds,
            solve(),
            resolution(2, 2),
            || false,
        ),
        Err(BinaryDiffusionNetworkError::Face { id: 8, .. })
    ));
    let duplicate = [link(3, [0, 1], &pairs), link(8, [1, 0], &pairs)];
    assert!(matches!(
        advance_candidate(
            &initial,
            &duplicate,
            0.1,
            &bounds,
            solve(),
            resolution(2, 2),
            || false,
        ),
        Err(BinaryDiffusionNetworkError::Invalid(_))
    ));
    let mut count = 0;
    advance_candidate(
        &initial,
        &links,
        0.1,
        &bounds,
        solve(),
        resolution(2, 2),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let mut calls = 0;
    assert_eq!(
        advance_candidate(
            &initial,
            &links,
            0.1,
            &bounds,
            solve(),
            resolution(2, 2),
            || {
                calls += 1;
                calls == count
            },
        ),
        Err(BinaryDiffusionNetworkError::Cancelled)
    );
    assert_eq!(initial, saved);
}

#[test]
fn g0_closed_faces_keep_original_states_without_pair_defaults() {
    let initial = BTreeMap::from([(0, gas(0.75)), (1, gas(0.25))]);
    let closed = DiffusionNetworkLink {
        geometry: GasDiffusionGeometry {
            area_m2: 0.0,
            ..geometry()
        },
        ..link(1, [0, 1], &[])
    };
    let next = advance_candidate(
        &initial,
        &[closed],
        100.0,
        &endpoints(&initial),
        solve(),
        resolution(2, 1),
        || false,
    )
    .unwrap();
    assert_eq!(next.parcels(), &initial);
    assert!(next.report().transfers.is_empty());
    assert_eq!(next.report().entropy_change_j_k, 0.0);
    assert!(next.report().last_sweep[0].initial_flux.pairs.is_empty());
}
