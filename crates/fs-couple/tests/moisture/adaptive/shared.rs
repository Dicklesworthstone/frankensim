//! Genuine specimen fields share temporal accuracy, not a fictitious geometry.
//! Independent continuous-time references use authored verification coefficients.
use super::*;
use fs_couple::moisture::{
    atmosphere::{CartesianMoistureParticipant, DryBasisSorption, MoistureParticipant},
    solution::SolutionMoistureError,
};

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
fn adaptive_shared_solution_resolves_opposing_fields_when_the_bath_is_stationary() {
    with_cx(|cx| {
        let target = 10.0;
        let volume = 1e-3;
        let bath = solution::reservoir(volume, vec![0.0]);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.3], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let fields = [participant(12, &a, &faces), participant(4, &b, &faces)];
        let loose = bath
            .advance_shared_adaptive_to(
                &fields,
                target,
                2,
                solution::controls(volume),
                time_controls(target, 1e-3),
                cx,
            )
            .unwrap();
        let temporal = time_controls(target, 1e-6);
        let fine = bath
            .advance_shared_adaptive_to(
                &fields,
                target,
                2,
                solution::controls(volume),
                temporal,
                cx,
            )
            .unwrap();
        // Both stores have k=4D/L². Their equal and opposite continuous-time
        // trajectories leave pure solvent activity AND water inventory fixed.
        let expected = [
            0.2 + 0.1 * (-0.04 * target).exp(),
            0.2 - 0.1 * (-0.04 * target).exp(),
        ];
        let error = |candidate: &SharedSolutionMoistureAdvanceCandidate| {
            candidate
                .solids()
                .iter()
                .zip(expected)
                .map(|(p, m)| (p.moisture.moisture()[0] - m).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine) < 0.25 * error(&loose));
        assert!(error(&fine) < 1e-4);
        let report = fine.report();
        assert!(report.summary.rejected_intervals > 0);
        assert!(report.summary.evaluations > loose.report().summary.evaluations);
        assert_eq!(report.summary.interval_s, [0.0, target]);
        assert_eq!(report.summary.activity, Some([1.0; 2]));
        near(
            fine.reservoir().solvent_mass_kg(),
            bath.solvent_mass_kg(),
            temporal.water_balance_tolerance_kg,
        );
        near(
            report.summary.water_change_kg,
            0.0,
            temporal.water_balance_tolerance_kg,
        );
        assert_eq!(
            report.specimens.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![4, 12]
        );
        for (i, (field, account)) in fine.solids().iter().zip(&report.specimens).enumerate() {
            assert_eq!(field.id, account.id);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(field.transport.interval_s[1], target);
            assert!(
                account.maximum_accepted_error_ratio > 0.0
                    && account.maximum_accepted_error_ratio <= 1.0
            );
            near(
                account.boundary_water_kg.iter().sum(),
                account.water_change_kg,
                temporal.water_balance_tolerance_kg,
            );
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 1);
            assert!(account.water_change_kg.abs() > 0.01);
        }
        let reversed = bath
            .advance_shared_adaptive_to(
                &[participant(4, &b, &faces), participant(12, &a, &faces)],
                target,
                2,
                solution::controls(volume),
                temporal,
                cx,
            )
            .unwrap();
        assert_eq!(fine.reservoir().solution(), reversed.reservoir().solution());
        assert_eq!(
            report.summary.evaluations,
            reversed.report().summary.evaluations
        );
        for (left, right) in fine.solids().iter().zip(reversed.solids()) {
            assert_eq!(left.id, right.id);
            assert_eq!(left.moisture.moisture(), right.moisture.moisture());
        }
        // The same two physical rates and dry inventories on real Cartesian
        // geometries exercise the shared six-face entry point independently.
        let ga = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.1], cx);
        let gb = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.3], cx);
        let sides = [faces[0].clone(), faces[1].clone(), None, None, None, None];
        let spatial_end = bath
            .advance_shared_cartesian_adaptive_to(
                &[spatial(12, &ga, &sides), spatial(4, &gb, &sides)],
                target,
                2,
                solution::controls(volume),
                temporal,
                cx,
            )
            .unwrap();
        for ((field, account), expected) in spatial_end
            .solids()
            .iter()
            .zip(&spatial_end.report().specimens)
            .zip(expected)
        {
            near(field.moisture.moisture()[0], expected, 1e-4);
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
            assert_eq!(field.moisture.geometry().lengths_m, [1.0, 0.02, 0.05]);
            near(
                field.transport.transport_end[0].query().axes()["moisture_content"],
                field.moisture.moisture()[0],
                0.0,
            );
        }
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.3]);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
    });
}

#[test]
fn adaptive_shared_vapor_matches_mean_and_difference_modes_with_one_gas_debit() {
    with_cx(|cx| {
        let target = 2.0;
        let r0 = 0.8;
        let air = atmosphere::atmosphere(r0, 1.0);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.20], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let temporal = time_controls(target, 2e-5);
        let endpoint = air
            .advance_shared_adaptive_to(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                target,
                2,
                atmosphere::controls(),
                temporal,
                cx,
            )
            .unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / r0;
        let dry = a.dry_mass_kg() + b.dry_mass_kg();
        let mean0 = 0.125;
        let k = 0.04;
        let slope = 0.25;
        let scale = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * r0 + slope * dry * mean0 / capacity) / scale;
        let mean = equilibrium + (mean0 - equilibrium) * (-k * scale * target).exp();
        let difference = 0.075 * (-k * target).exp();
        let expected = [mean - difference, mean + difference];
        let expected_rh = r0 - dry * (mean - mean0) / capacity;
        near(
            endpoint.reservoir().relative_humidity().unwrap(),
            expected_rh,
            2e-3,
        );
        assert!(endpoint.report().summary.rejected_intervals > 0);
        for (i, (field, account)) in endpoint
            .solids()
            .iter()
            .zip(&endpoint.report().specimens)
            .enumerate()
        {
            near(field.moisture.moisture()[0], expected[i], 5e-5);
            assert_eq!(field.moisture.time_s(), target);
            assert_eq!(account.water_change_kg.is_sign_positive(), i == 0);
            near(
                account.boundary_water_kg.iter().sum(),
                account.water_change_kg,
                temporal.water_balance_tolerance_kg,
            );
            for source in field.transport.surface_sorption.iter().flatten() {
                near(
                    source.query().axes()["relative_humidity"],
                    endpoint.reservoir().relative_humidity().unwrap(),
                    1e-11,
                );
            }
        }
        let report = &endpoint.report().summary;
        near(
            endpoint.reservoir().vapor_mass_kg().unwrap() - air.vapor_mass_kg().unwrap()
                + report.water_change_kg,
            0.0,
            temporal.water_balance_tolerance_kg,
        );
        near(
            atmosphere::carrier_mass(endpoint.reservoir()),
            atmosphere::carrier_mass(&air),
            1e-12,
        );
        assert!(
            report.maximum_carrier_mass_residual_kg.unwrap() <= temporal.water_balance_tolerance_kg
        );
        assert!(
            report.total_mass_residual_kg.unwrap().abs() <= temporal.water_balance_tolerance_kg
        );
        assert_eq!(report.last_substep.interval_s[1], target);
        let ga = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.05], cx);
        let gb = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.20], cx);
        let sides = [faces[0].clone(), faces[1].clone(), None, None, None, None];
        let spatial_end = air
            .advance_shared_cartesian_adaptive_to(
                &[spatial(2, &gb, &sides), spatial(1, &ga, &sides)],
                target,
                2,
                atmosphere::controls(),
                temporal,
                cx,
            )
            .unwrap();
        near(
            spatial_end.reservoir().relative_humidity().unwrap(),
            expected_rh,
            2e-3,
        );
        for (i, (field, account)) in spatial_end
            .solids()
            .iter()
            .zip(&spatial_end.report().specimens)
            .enumerate()
        {
            near(field.moisture.moisture()[0], expected[i], 5e-5);
            assert_eq!(&account.boundary_water_kg[2..], &[0.0; 4]);
            near(
                account.water_change_kg,
                field.moisture.dry_mass_kg() * (field.moisture.moisture()[0] - [0.05, 0.20][i]),
                1e-14,
            );
        }
        assert_eq!(a.moisture(), &[0.05]);
        assert_eq!(b.moisture(), &[0.20]);
        near(air.relative_humidity().unwrap(), r0, 1e-14);
    });
}

fn source_limited_slab(cx: &Cx<'_>) -> MoistureSlab {
    let base = card([0.01; 2], 0.3);
    let mut claim = base.claims_for("moisture_diffusivity")[0].1.clone();
    claim.key = PropertyKey::new("limited_diffusivity", MOISTURE_DIFFUSIVITY_DIMS);
    claim.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
            293.15,
        )
        .with_quantity(
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            0.0,
            0.15,
        );
    let mut claims = base.claims().clone();
    claims.insert_claim(claim).unwrap();
    let material =
        MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
    MoistureSlab::new(
        material,
        point(0.1),
        "limited_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        geometry(1.0),
        vec![0.1],
        cx,
    )
    .unwrap()
}

#[test]
fn shared_adaptive_late_work_domain_and_source_refusals_discard_the_entire_window() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = solution::reservoir(volume, vec![0.0]);
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let b = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.25)),
            Some(atmosphere::sorption(0.25)),
        ];
        let temporal = time_controls(10.0, 1.0);
        // The first two fine steps leave V=5.4166...e-5 m³ inside the
        // declared [5e-5,1.5e-4] bracket. The next interval exhausts it.
        assert!(matches!(
            bath.advance_shared_adaptive_to(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                20.0,
                2,
                solution::controls(volume),
                temporal,
                cx
            ),
            Err(MoistureAdvanceError::Step {
                time_s,
                evaluations: 4,
                error: SolutionMoistureError::Unbracketed { .. },
                ..
            }) if time_s == 10.0
        ));
        let mut short = temporal;
        short.maximum_evaluations = 3;
        assert!(matches!(
            bath.advance_shared_adaptive_to(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                20.0,
                2,
                solution::controls(volume),
                short,
                cx
            ),
            Err(MoistureAdvanceError::WorkLimit {
                time_s,
                accepted_intervals: 1,
                evaluations: 3
            }) if time_s == 10.0
        ));
        let unsupported = source_limited_slab(cx);
        assert!(matches!(
            bath.advance_shared_adaptive_to(
                &[
                    participant(1, &a, &faces),
                    participant(2, &unsupported, &faces)
                ],
                20.0,
                2,
                solution::controls(volume),
                temporal,
                cx
            ),
            Err(MoistureAdvanceError::Step {
                evaluations: 1,
                error: SolutionMoistureError::Moisture(MoistureError::Material(_)),
                ..
            })
        ));
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
        assert_eq!(bath.solution().solute_moles(), &[0.0]);
        assert_eq!(a.moisture(), &[0.1]);
        assert_eq!(b.moisture(), &[0.1]);
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
        let air = atmosphere::atmosphere(0.8, 1.0);
        let a = grid([1.0; 3], [1; 3], 0.01, vec![0.1], &cx);
        let b = grid([1.0; 3], [1; 3], 0.01, vec![0.2], &cx);
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(0.25)));
        gate.request();
        assert!(matches!(
            air.advance_shared_cartesian_adaptive_to(
                &[spatial(1, &a, &faces), spatial(2, &b, &faces)],
                1.0,
                2,
                atmosphere::controls(),
                time_controls(1.0, 1e-3),
                &cx
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
        near(air.relative_humidity().unwrap(), 0.8, 1e-14);
    });
}
