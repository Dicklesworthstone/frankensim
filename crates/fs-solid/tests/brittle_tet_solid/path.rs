//! G1/G3: adaptive loading uses the actual source-bound P1 solve and history.
use super::*;

use fs_solid::brittle::solid::path::{
    BrittleTetPathControls, BrittleTetPathError, BrittleTetRefinementCause,
};

fn path_controls() -> BrittleTetPathControls {
    BrittleTetPathControls {
        initial_fraction: 1.0,
        minimum_fraction: 1.0 / 65_536.0,
        maximum_fraction: 1.0,
        maximum_accepted_steps: 512,
        maximum_trials: 4096,
        maximum_damage_increment: 0.08,
        initial_equilibrium_tolerance_n: 1e-6,
        work: BrittleTetControls {
            maximum_endpoint_work_defect_j: 1e-3,
            ..controls()
        },
    }
}

#[test]
fn adaptive_native_path_resolves_softening_and_complete_surface_work() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "adaptive fracture");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let band = &virgin.band_geometry()[0];
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * band.band_width_m());
        let strain = 1.2 * failure;
        let supports = prescribed(0.03, strain);
        let settings = path_controls();
        let one_shot = virgin.advance(
            &supports,
            &[],
            IncrementSettings::default(),
            settings.work,
            cx,
        );
        assert!(matches!(
            one_shot,
            Err(BrittleTetError::EnergyGate {
                account: "signed endpoint integration defect",
                ..
            })
        ));
        let step = virgin
            .advance_path(&supports, &[], IncrementSettings::default(), settings, cx)
            .unwrap();
        assert!(step.rejected_trials > 0);
        assert_eq!(step.trials, step.records.len() + step.rejected_trials);
        assert_eq!(step.records.last().unwrap().fraction, 1.0);
        assert!(step.records.len() > 4);
        assert!(step.solid.solution().element_states[0].failed());
        close(step.solid.solution().element_strain[0][0], strain, 2e-12);
        let fracture = FRACTURE_WORK * band.projected_area_m2();
        close(step.solid.fracture_dissipation_j(), fracture, 2e-9);
        close(step.energy.fracture_dissipation_j, fracture, 2e-9);
        close(step.energy.material_work_j, fracture, 2e-9);
        assert!(step.energy.stored_after_j.abs() < 1e-10);
        assert!(step.energy.material_balance_error_j.abs() < 1e-9);
        assert!(step.energy.equilibrium_work_residual_j.abs() < 1e-9);
        let mut previous_fraction = 0.0;
        let mut previous_damage = 0.0;
        let mut absolute_defect = 0.0;
        let mut signed_defect = 0.0;
        for record in &step.records {
            assert!(record.fraction > previous_fraction);
            assert!(record.maximum_damage >= previous_damage);
            assert!(record.maximum_damage_increment <= settings.maximum_damage_increment);
            close(
                record.maximum_damage - previous_damage,
                record.maximum_damage_increment,
                2e-11,
            );
            let epsilon = record.fraction * strain;
            let stress = if epsilon <= onset {
                YOUNG * epsilon
            } else if epsilon < failure {
                STRENGTH * (failure - epsilon) / (failure - onset)
            } else {
                0.0
            };
            // Independent uniaxial endpoint force times prescribed motion.
            let expected_work = band.reference_volume_m3()
                * stress
                * ((record.fraction - previous_fraction) * strain);
            assert!((record.energy.support_endpoint_work_j - expected_work).abs() < 1e-9);
            assert!(
                record.energy.endpoint_work_defect_j.abs()
                    <= settings.work.maximum_endpoint_work_defect_j
                        * (record.fraction - previous_fraction)
            );
            absolute_defect += record.energy.endpoint_work_defect_j.abs();
            signed_defect += record.energy.endpoint_work_defect_j;
            previous_fraction = record.fraction;
            previous_damage = record.maximum_damage;
        }
        close(step.absolute_endpoint_work_defect_j, absolute_defect, 2e-12);
        close(step.energy.endpoint_work_defect_j, signed_defect, 2e-10);
        assert!(
            step.absolute_endpoint_work_defect_j <= settings.work.maximum_endpoint_work_defect_j
        );
        assert!(step.absolute_endpoint_work_defect_j >= step.energy.endpoint_work_defect_j.abs());
        assert_eq!(virgin.fracture_dissipation_j(), 0.0);
        assert_eq!(virgin.solution().element_states[0].damage(), 0.0);
        assert!(std::ptr::eq(
            step.solid.materials()[0].law(),
            materials[0].law()
        ));
        let replay = virgin
            .advance_path(&supports, &[], IncrementSettings::default(), settings, cx)
            .unwrap();
        assert_eq!(step.solid.solution(), replay.solid.solution());
        assert_eq!(step.records, replay.records);
        assert_eq!(step.energy, replay.energy);
        assert_eq!(step.trials, replay.trials);
    });
}

#[test]
fn retained_damage_survives_unloading_reloading_and_compressive_reclosure() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "adaptive reversal");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let peak = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let damage = peak.solution().element_states[0].damage();
        let loss = peak.fracture_dissipation_j();
        let stress = peak.solution().element_stress[0][0];
        let unloaded = peak
            .advance_path(
                &prescribed(0.03, 0.0005),
                &[],
                IncrementSettings::default(),
                path_controls(),
                cx,
            )
            .unwrap();
        close(
            unloaded.solid.solution().element_stress[0][0],
            0.25 * stress,
            2e-9,
        );
        let reloaded = unloaded
            .solid
            .advance_path(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                path_controls(),
                cx,
            )
            .unwrap();
        close(reloaded.solid.solution().element_stress[0][0], stress, 2e-9);
        let closed = reloaded
            .solid
            .advance_path(
                &prescribed(0.03, -0.0002),
                &[],
                IncrementSettings::default(),
                path_controls(),
                cx,
            )
            .unwrap();
        close(
            closed.solid.solution().element_stress[0][0],
            -0.0002 * YOUNG,
            2e-9,
        );
        for step in [&unloaded, &reloaded, &closed] {
            close(
                step.solid.solution().element_states[0].damage(),
                damage,
                2e-10,
            );
            close(step.solid.fracture_dissipation_j(), loss, 2e-10);
            assert!(step.energy.fracture_dissipation_j.abs() < 1e-10);
            assert!(
                step.records
                    .iter()
                    .all(|r| r.maximum_damage_increment < 1e-12)
            );
            assert!(step.energy.material_balance_error_j.abs() < 1e-9);
        }
    });
}

#[test]
fn explicit_force_endpoints_keep_initial_equilibrium_and_physical_work() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "adaptive force");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let volume = virgin.band_geometry()[0].reference_volume_m3();
        let force = 0.5 * STRENGTH * volume / 0.03;
        let fixed: Vec<_> = prescribed(0.03, 0.0)
            .into_iter()
            .filter(|(dof, _)| *dof != 3)
            .collect();
        let fine = BrittleTetPathControls {
            work: BrittleTetControls {
                maximum_endpoint_work_defect_j: 1e-5,
                ..controls()
            },
            ..path_controls()
        };
        let step = virgin
            .advance_path(
                &fixed,
                &[(3, [0.0, 0.5 * force]), (3, [0.0, 0.5 * force])],
                IncrementSettings::default(),
                fine,
                cx,
            )
            .unwrap();
        close(
            step.solid.solution().element_strain[0][0],
            0.5 * STRENGTH / YOUNG,
            2e-9,
        );
        close(
            step.energy.material_work_j,
            0.5 * volume * (0.5 * STRENGTH).powi(2) / YOUNG,
            2e-9,
        );
        assert_eq!(step.energy.support_endpoint_work_j, 0.0);
        assert!(step.energy.dead_load_endpoint_work_j > step.energy.material_work_j);
        assert!(step.rejected_trials > 0);
        assert!(matches!(
            step.solid.advance_path(
                &fixed,
                &[(3, [0.0, force])],
                IncrementSettings::default(),
                fine,
                cx,
            ),
            Err(BrittleTetPathError::InitialEquilibrium { .. })
        ));
        let unloaded = step
            .solid
            .advance_path(
                &fixed,
                &[(3, [force, 0.5 * force])],
                IncrementSettings::default(),
                fine,
                cx,
            )
            .unwrap();
        close(
            unloaded.solid.solution().element_strain[0][0],
            0.25 * STRENGTH / YOUNG,
            2e-9,
        );
        assert!(unloaded.energy.material_work_j < 0.0);
        assert_eq!(unloaded.solid.fracture_dissipation_j(), 0.0);
    });
}

#[test]
fn rejected_whole_intervals_leave_the_original_solid_available_for_retry() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "adaptive refusal");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let snapshot = virgin.solution().clone();
        let target = prescribed(0.03, 0.5 * STRENGTH / YOUNG);
        let fine = BrittleTetPathControls {
            work: BrittleTetControls {
                maximum_endpoint_work_defect_j: 1e-5,
                ..controls()
            },
            ..path_controls()
        };
        assert!(matches!(
            virgin.advance_path(
                &target,
                &[],
                IncrementSettings::default(),
                BrittleTetPathControls {
                    maximum_accepted_steps: 1,
                    ..fine
                },
                cx,
            ),
            Err(BrittleTetPathError::WorkLimit {
                accepted_steps: 1,
                ..
            })
        ));
        assert!(matches!(
            virgin.advance_path(
                &target,
                &[],
                IncrementSettings::default(),
                BrittleTetPathControls {
                    maximum_trials: 1,
                    ..fine
                },
                cx,
            ),
            Err(BrittleTetPathError::WorkLimit {
                accepted_steps: 0,
                trials: 1
            })
        ));
        assert!(matches!(
            virgin.advance_path(
                &target,
                &[],
                IncrementSettings::default(),
                BrittleTetPathControls {
                    minimum_fraction: 1.0,
                    ..fine
                },
                cx,
            ),
            Err(BrittleTetPathError::MinimumInterval {
                cause: BrittleTetRefinementCause::Endpoint(BrittleTetError::EnergyGate {
                    account: "signed endpoint integration defect",
                    ..
                }),
                ..
            })
        ));
        assert!(matches!(
            virgin.advance_path(
                &target,
                &[],
                IncrementSettings::default(),
                BrittleTetPathControls {
                    minimum_fraction: 1.0,
                    maximum_damage_increment: 0.001,
                    work: controls(),
                    ..path_controls()
                },
                cx,
            ),
            Ok(_)
        )); // This endpoint stays elastic, so the damage gate is inactive.
        let damaging = prescribed(0.03, 0.002);
        assert!(matches!(
            virgin.advance_path(
                &damaging,
                &[],
                IncrementSettings::default(),
                BrittleTetPathControls {
                    minimum_fraction: 1.0,
                    maximum_damage_increment: 0.001,
                    work: controls(),
                    ..path_controls()
                },
                cx,
            ),
            Err(BrittleTetPathError::MinimumInterval {
                cause: BrittleTetRefinementCause::DamageIncrement { .. },
                ..
            })
        ));
        assert!(matches!(
            virgin.advance_path(
                &target,
                &[],
                IncrementSettings {
                    max_newton_iterations: 0,
                    ..IncrementSettings::default()
                },
                fine,
                cx,
            ),
            Err(BrittleTetPathError::Endpoint(BrittleTetError::Solve(
                SmallStrainTetError::InvalidInput { .. }
            )))
        ));
        assert!(matches!(
            virgin.advance_path(
                &[(0, 0.0), (0, 0.0)],
                &[],
                IncrementSettings::default(),
                fine,
                cx,
            ),
            Err(BrittleTetPathError::Endpoint(BrittleTetError::Solve(
                SmallStrainTetError::InvalidInput { .. }
            )))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                virgin.advance_path(&target, &[], IncrementSettings::default(), fine, cancelled,),
                Err(BrittleTetPathError::Endpoint(BrittleTetError::Cancelled))
            ));
        });
        assert_eq!(virgin.solution(), &snapshot);
        assert_eq!(virgin.fracture_dissipation_j(), 0.0);
        let resumed = virgin
            .advance_path(&target, &[], IncrementSettings::default(), fine, cx)
            .unwrap();
        close(
            resumed.solid.solution().element_strain[0][0],
            0.5 * STRENGTH / YOUNG,
            2e-9,
        );
        assert!(
            resumed.absolute_endpoint_work_defect_j <= fine.work.maximum_endpoint_work_defect_j
        );
    });
}
