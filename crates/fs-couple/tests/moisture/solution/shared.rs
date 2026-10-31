//! One aqueous inventory with independently evolving, original specimens.
//! These analytical/manufactured fixtures are authored, not calibration data.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::atmosphere::DryBasisSorption;
use fs_couple::moisture::cartesian::{MoistureCartesian, MoistureCartesianGeometry};

fn participant<'a>(
    id: u64,
    slab: &'a MoistureSlab,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> MoistureParticipant<'a> {
    MoistureParticipant {
        id,
        slab,
        faces,
        transport: config(),
    }
}
fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}
fn grid(
    length: f64,
    cells: [usize; 3],
    d: f64,
    initial: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureCartesian {
    MoistureCartesian::new(
        card([d; 2], 0.3),
        point(0.1),
        ["moisture_diffusivity"; 3],
        MaterialPropertySelection::SingleClaimOnly,
        MoistureCartesianGeometry {
            lengths_m: [length; 3],
            cells,
            frame: ContentHash([4; 32]),
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        initial,
        cx,
    )
    .unwrap()
}
fn spatial<'a>(
    id: u64,
    moisture: &'a MoistureCartesian,
    faces: &'a [Option<DryBasisSorption>; 6],
) -> CartesianMoistureParticipant<'a> {
    CartesianMoistureParticipant {
        id,
        moisture,
        faces,
        transport: config(),
        nonlinear: nonlinear(),
    }
}

#[test]
fn shared_solution_has_one_manufactured_activity_and_keeps_opposing_local_transfers() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let activity = 0.95_f64;
        let slope = [0.2, 0.25];
        let alpha = [0.08, 4.0 * 2.0 * 0.03 / 0.7_f64.powi(2)];
        let dry = [0.5, 0.35];
        let weight = [
            dry[0] * alpha[0] / (1.0 + alpha[0]),
            dry[1] * alpha[1] / (1.0 + alpha[1]),
        ];
        let fa = [
            Some(atmosphere::sorption(slope[0])),
            Some(atmosphere::sorption(slope[0])),
        ];
        let fb = [
            Some(atmosphere::sorption(slope[1])),
            Some(atmosphere::sorption(slope[1])),
        ];
        for balanced in [false, true] {
            // In the second case nonzero local exchanges cancel. Net reservoir
            // equilibrium must not suppress either specimen's physical response.
            let old = [
                0.17,
                if balanced {
                    slope[1] * activity + weight[0] / weight[1] * (slope[0] * activity - 0.17)
                } else {
                    0.26
                },
            ];
            let expected: [f64; 2] = std::array::from_fn(|i| {
                (old[i] + alpha[i] * slope[i] * activity) / (1.0 + alpha[i])
            });
            let change: [f64; 2] = std::array::from_fn(|i| dry[i] * (expected[i] - old[i]));
            let initial_volume = final_volume + (change[0] + change[1]) / WATER_DENSITY;
            let amount = -activity.ln() * final_volume / WATER_VOLUME;
            let bath = reservoir(initial_volume, vec![0.6 * amount, 0.4 * amount]);
            let a = slab(card([0.01; 2], 0.3), 1.0, vec![old[0]], cx);
            let b = slab(card([0.03; 2], 0.3), 0.7, vec![old[1]], cx);
            let endpoint = bath
                .advance_shared(
                    &[participant(20, &a, &fa), participant(5, &b, &fb)],
                    2.0,
                    2,
                    controls(initial_volume),
                    cx,
                )
                .unwrap();
            let reversed = bath
                .advance_shared(
                    &[participant(5, &b, &fb), participant(20, &a, &fa)],
                    2.0,
                    2,
                    controls(initial_volume),
                    cx,
                )
                .unwrap();
            assert_eq!(
                endpoint.solids().iter().map(|p| p.id).collect::<Vec<_>>(),
                vec![5, 20]
            );
            assert_eq!(
                endpoint.reservoir().solution(),
                reversed.reservoir().solution()
            );
            near(
                endpoint.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            near(
                endpoint.reservoir().solution().water_activity().unwrap(),
                activity,
                1e-12,
            );
            for (id, i) in [(20, 0), (5, 1)] {
                let body = endpoint.solids().iter().find(|p| p.id == id).unwrap();
                let repeat = reversed.solids().iter().find(|p| p.id == id).unwrap();
                near(body.moisture.moisture()[0], expected[i], 2e-13);
                near(
                    body.transport.boundary_water_kg.iter().sum(),
                    change[i],
                    2e-13,
                );
                assert_eq!(body.moisture.time_s(), 2.0);
                assert_eq!(body.moisture.moisture(), repeat.moisture.moisture());
                assert_eq!(
                    body.transport.boundary_water_kg,
                    repeat.transport.boundary_water_kg
                );
                for source in body.transport.surface_sorption.iter().flatten() {
                    near(source.query().axes()["relative_humidity"], activity, 1e-12);
                }
            }
            assert!(change[0] > 1e-4 && change[1] < -1e-4);
            near(
                endpoint.step().water_to_solids_kg,
                change[0] + change[1],
                2e-13,
            );
            near(endpoint.step().water_balance_residual_kg, 0.0, 1e-13);
            near(endpoint.step().direct_water_balance_residual_kg, 0.0, 1e-13);
            assert_eq!(
                endpoint.reservoir().solution().solute_moles(),
                bath.solution().solute_moles()
            );
            if balanced {
                assert_eq!(endpoint.step().iterations, 0);
                assert_eq!(
                    endpoint.reservoir().solution().solvent_volume_m3(),
                    initial_volume
                );
                near(endpoint.step().water_to_solids_kg, 0.0, 1e-13);
            }
            assert_eq!(bath.solution().solvent_volume_m3(), initial_volume);
            assert_eq!(a.moisture(), &[old[0]]);
            assert_eq!(b.moisture(), &[old[1]]);
            assert_eq!(a.time_s(), 0.0);
            assert_eq!(b.time_s(), 0.0);
        }
    });
}

#[test]
fn shared_cartesian_solution_uses_nonideal_activity_and_advances_sealed_companions() {
    with_cx(|cx| {
        let final_volume = 1e-4;
        let activity = 0.95_f64;
        let dt = 0.005;
        let length = [0.1_f64, 0.2];
        let diffusion = [0.01, 0.02];
        let dry = [0.5, 4.0];
        let old = [0.1, 0.26];
        let slope = [0.25, 0.2];
        let alpha: [f64; 2] = std::array::from_fn(|i| 12.0 * dt * diffusion[i] / length[i].powi(2));
        let expected: [f64; 2] =
            std::array::from_fn(|i| (old[i] + alpha[i] * slope[i] * activity) / (1.0 + alpha[i]));
        let water: f64 = (0..2).map(|i| dry[i] * (expected[i] - old[i])).sum();
        let initial_volume = final_volume + water / WATER_DENSITY;
        let interaction = 0.002;
        // Independent inversion of -ln(a_w)/v_water = c + B*c²/2.
        let osmotic_concentration = -activity.ln() / WATER_VOLUME;
        let concentration = 2.0 * osmotic_concentration
            / (1.0 + (1.0 + 2.0 * interaction * osmotic_concentration).sqrt());
        let amount = concentration * final_volume;
        let model = Arc::new(
            VirialSolutionModel::new(293.15, WATER_VOLUME, 1.0, vec![3000.0], vec![interaction])
                .unwrap(),
        );
        let bath = FiniteMoistureSolution::new(
            AqueousSolution::nonideal(model.clone(), initial_volume, vec![amount]).unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        let a = grid(length[0], [1; 3], diffusion[0], vec![old[0]], cx);
        let b = grid(length[1], [1; 3], diffusion[1], vec![old[1]], cx);
        let sealed_profile: Vec<_> = (0..8)
            .map(|i| if i % 2 == 0 { 0.08 } else { 0.12 })
            .collect();
        let c = grid(0.02, [2; 3], 0.01, sealed_profile.clone(), cx);
        let fa = std::array::from_fn(|_| Some(atmosphere::sorption(slope[0])));
        let fb = std::array::from_fn(|_| Some(atmosphere::sorption(slope[1])));
        let sealed = std::array::from_fn(|_| None);
        let endpoint = bath
            .advance_shared_cartesian(
                &[
                    spatial(30, &c, &sealed),
                    spatial(20, &b, &fb),
                    spatial(10, &a, &fa),
                ],
                dt,
                10,
                controls(initial_volume),
                cx,
            )
            .unwrap();
        near(
            endpoint.reservoir().solution().water_activity().unwrap(),
            activity,
            1e-12,
        );
        near(
            endpoint.reservoir().solution().solvent_volume_m3(),
            final_volume,
            2e-16,
        );
        for (i, body) in endpoint.solids()[..2].iter().enumerate() {
            near(body.moisture.moisture()[0], expected[i], 2e-13);
            near(
                body.transport.boundary_water_kg.iter().sum(),
                dry[i] * (expected[i] - old[i]),
                2e-13,
            );
            near(
                body.transport.transport_start[0].query().axes()["moisture_content"],
                old[i],
                0.0,
            );
            near(
                body.transport.transport_end[0].query().axes()["moisture_content"],
                body.moisture.moisture()[0],
                0.0,
            );
            for source in body.transport.surface_sorption.iter().flatten() {
                near(source.query().axes()["relative_humidity"], activity, 1e-12);
            }
        }
        let companion = &endpoint.solids()[2];
        assert_eq!(companion.id, 30);
        assert_eq!(companion.transport.boundary_water_kg, [0.0; 6]);
        assert!(companion.moisture.moisture()[0] > 0.08);
        assert!(companion.moisture.moisture()[1] < 0.12);
        near(companion.moisture.water_mass_kg(), c.water_mass_kg(), 1e-13);
        assert!(endpoint.solids().iter().all(|p| p.moisture.time_s() == dt));
        near(endpoint.step().water_to_solids_kg, water, 2e-13);
        near(endpoint.step().water_balance_residual_kg, 0.0, 1e-13);
        near(endpoint.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(endpoint.reservoir().solution().solute_moles(), &[amount]);
        assert!(std::ptr::eq(
            endpoint.reservoir().solution().nonideal_model().unwrap(),
            model.as_ref()
        ));
        // The sealed companion still consumes eight real cells of the budget.
        assert!(
            bath.advance_shared_cartesian(
                &[
                    spatial(10, &a, &fa),
                    spatial(20, &b, &fb),
                    spatial(30, &c, &sealed)
                ],
                dt,
                9,
                controls(initial_volume),
                cx
            )
            .is_err()
        );
        assert_eq!(c.moisture(), sealed_profile);
        assert_eq!(c.time_s(), 0.0);
    });
}

#[test]
fn pure_shared_solution_debits_total_uptake_once_and_preserves_single_specimen_api() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = reservoir(volume, vec![0.0]);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let b = slab(card([0.02; 2], 0.3), 0.5, vec![0.08], cx);
        let fa = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let fb = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let expected = [(0.1 + 0.02 * 0.2) / 1.02, (0.08 + 0.16 * 0.25) / 1.16];
        let water = 0.5 * (expected[0] - 0.1) + 0.25 * (expected[1] - 0.08);
        let endpoint = bath
            .advance_shared(
                &[participant(1, &a, &fa), participant(2, &b, &fb)],
                0.5,
                2,
                controls(volume),
                cx,
            )
            .unwrap();
        assert_eq!(endpoint.step().water_activity, [1.0; 2]);
        assert_eq!(endpoint.step().iterations, 0);
        near(
            endpoint.solids()[0].moisture.moisture()[0],
            expected[0],
            1e-14,
        );
        near(
            endpoint.solids()[1].moisture.moisture()[0],
            expected[1],
            1e-14,
        );
        near(
            bath.solvent_mass_kg() - endpoint.reservoir().solvent_mass_kg(),
            water,
            1e-13,
        );
        near(endpoint.step().water_to_solids_kg, water, 1e-14);
        let one = bath
            .advance_shared(&[participant(8, &a, &fa)], 0.5, 1, controls(volume), cx)
            .unwrap();
        let legacy = bath
            .advance(&a, 0.5, &fa, config(), controls(volume), cx)
            .unwrap();
        assert_eq!(one.reservoir().solution(), legacy.reservoir().solution());
        assert_eq!(
            one.solids()[0].moisture.moisture(),
            legacy.moisture().moisture()
        );
        assert_eq!(
            one.solids()[0].transport.boundary_water_kg,
            legacy.step().transport.boundary_water_kg
        );
    });
}

fn limited_sorption() -> DryBasisSorption {
    let key = "limited_shared_solution_sorption";
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(
                key,
                qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ),
            value: PropertyValue::Curve {
                abscissa: "relative_humidity".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 0.0), (1.0, 0.25)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.0,
                    0.9,
                )
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    293.15,
                    293.15,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored limited sorption source".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored sorption".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    DryBasisSorption::new(card, point, key, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn shared_solution_late_source_work_clock_and_cancel_refusals_publish_nothing() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = reservoir(volume, vec![0.5]);
        let saved = bath.solution().clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let b = slab(card([0.02; 2], 0.3), 1.0, vec![0.08], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let unsupported = [Some(limited_sorption()), Some(limited_sorption())];
        assert!(bath.solution().water_activity().unwrap() > 0.9);
        assert!(
            bath.advance(&a, 2.0, &faces, config(), controls(volume), cx)
                .unwrap()
                .moisture()
                .moisture()[0]
                > a.moisture()[0]
        );
        assert!(matches!(
            bath.advance_shared(
                &[participant(1, &a, &faces), participant(2, &b, &unsupported)],
                2.0,
                2,
                controls(volume),
                cx
            ),
            Err(SolutionMoistureError::Moisture(MoistureError::Material(_)))
        ));
        let mut short = controls(volume);
        short.maximum_iterations = 1;
        assert!(matches!(
            bath.advance_shared(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                2.0,
                2,
                short,
                cx
            ),
            Err(SolutionMoistureError::NotConverged { iterations: 1, .. })
        ));
        assert!(
            bath.advance_shared(
                &[participant(1, &a, &faces), participant(1, &b, &faces)],
                2.0,
                2,
                controls(volume),
                cx
            )
            .is_err()
        );
        let later = a
            .advance(
                0.1,
                &[MoistureFace::Sealed, MoistureFace::Sealed],
                config(),
                cx,
            )
            .unwrap()
            .0;
        assert!(
            bath.advance_shared(
                &[participant(1, &later, &faces), participant(2, &b, &faces)],
                2.0,
                2,
                controls(volume),
                cx
            )
            .is_err()
        );
        assert_eq!(bath.solution(), &saved);
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.08]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let bath = reservoir(1e-4, vec![0.02]);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], &cx);
        let b = slab(card([0.02; 2], 0.3), 1.0, vec![0.08], &cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        gate.request();
        assert!(matches!(
            bath.advance_shared(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                2.0,
                2,
                controls(1e-4),
                &cx
            ),
            Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
        assert_eq!(bath.solution().solvent_volume_m3(), 1e-4);
        assert_eq!(bath.solution().solute_moles(), &[0.02]);
    });
}
