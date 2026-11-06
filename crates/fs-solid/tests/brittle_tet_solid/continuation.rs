//! G1/G3 native arclength fracture: analytical traction and total Gc*A,
//! source/history continuity, physical reversals and transactional refusal.

use super::*;
use fs_solid::brittle::solid::continuation::{
    BrittleTetArcBudget, BrittleTetArcControls, BrittleTetArcDirection, BrittleTetContinuation,
    BrittleTetContinuationError, BrittleTetPathScales,
};
use fs_solid::continuation::ArcSettings;

fn arc_settings() -> ArcSettings {
    ArcSettings {
        ds: 0.2,
        ds_min: 1e-6,
        ds_max: 1.0,
        max_corrector_iters: 24,
        tol: 3e-10,
    }
}

fn arc_controls() -> BrittleTetArcControls {
    BrittleTetArcControls {
        work: controls(),
        maximum_damage_increment: 1.0,
    }
}

fn arc_budget() -> BrittleTetArcBudget {
    BrittleTetArcBudget {
        maximum_steps: 128,
        maximum_operator_evaluations: 256,
    }
}

fn fixed_supports() -> Vec<(usize, f64)> {
    prescribed(0.03, 0.0)
        .into_iter()
        .filter(|(dof, _)| *dof != 3)
        .collect()
}

fn scales(solid: &BrittleTetSolid<'_>) -> BrittleTetPathScales {
    BrittleTetPathScales {
        displacement_m: 0.03 * STRENGTH / YOUNG,
        force_n: STRENGTH * solid.band_geometry()[0].reference_volume_m3() / 0.03,
    }
}

#[test]
fn native_arc_crosses_tensile_peak_and_preserves_complete_fracture_energy() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "arclength fracture");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let units = scales(&virgin);
        let force = units.force_n;
        // Duplicate entries and a force on an eliminated support stay physical.
        let loads = [(3, 0.25 * force), (3, 0.75 * force), (0, 2.0 * force)];
        let fixed = fixed_supports();
        let settings = arc_settings();
        let mut path = BrittleTetContinuation::try_from_solid(
            &virgin,
            &fixed,
            &loads,
            0.0,
            units,
            arc_budget(),
            &settings,
            cx,
        )
        .unwrap();
        assert_eq!(path.free_dofs(), &[3, 7, 11]);
        assert_eq!(path.reference_force_n()[3], force);
        let band = &virgin.band_geometry()[0];
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * band.band_width_m());
        let mut previous_strain = 0.0;
        let mut previous_damage = 0.0;
        let mut previous_loss = 0.0;
        let mut exact_work = 0.0;
        let mut load_turns = 0;
        for index in 0..90 {
            let step = path.advance(&settings, arc_controls(), None, cx).unwrap();
            if index == 24 {
                let replay = path.advance(&settings, arc_controls(), None, cx).unwrap();
                assert_eq!(
                    step.continuation.solid().solution(),
                    replay.continuation.solid().solution()
                );
                assert_eq!(step.continuation.trace(), replay.continuation.trace());
                assert_eq!(step.energy, replay.energy);
            }
            assert!(step.force_residual_norm_n <= settings.tol * force);
            assert!(step.operator_evaluations <= arc_budget().maximum_operator_evaluations);
            assert!(step.energy.material_balance_error_j.abs() < 1e-9);
            assert!(step.energy.equilibrium_work_residual_j.abs() < 1e-9);
            assert_eq!(step.energy.support_endpoint_work_j, 0.0);
            let solid = step.continuation.solid();
            let solution = solid.solution();
            let strain = solution.element_strain[0][0];
            let expected_stress = if strain < onset {
                YOUNG * strain
            } else if strain < failure {
                STRENGTH * (failure - strain) / (failure - onset)
            } else {
                0.0
            };
            assert!(strain > previous_strain);
            assert!((solution.element_stress[0][0] - expected_stress).abs() < 3e-8 * STRENGTH);
            assert!((step.continuation.load_factor() - expected_stress / STRENGTH).abs() < 3e-8);
            for component in &solution.element_stress[0][1..] {
                assert!(component.abs() < 3e-8 * STRENGTH);
            }
            let reaction_x: f64 = fixed
                .iter()
                .zip(&solution.increments[0].reactions_n)
                .filter(|((dof, _), _)| dof % 3 == 0)
                .map(|(_, reaction)| reaction)
                .sum();
            assert!(
                (reaction_x + 3.0 * step.continuation.load_factor() * force).abs() < 1e-8 * force
            );
            let point = &solution.element_states[0];
            assert!(point.damage() >= previous_damage);
            assert!(solid.fracture_dissipation_j() >= previous_loss);
            assert!(std::ptr::eq(solid.materials()[0].law(), materials[0].law()));
            previous_strain = strain;
            previous_damage = point.damage();
            previous_loss = solid.fracture_dissipation_j();
            exact_work += step.energy.material_work_j;
            load_turns += usize::from(step.load_direction_changed);
            path = step.continuation;
        }
        assert!(path.solid().solution().element_states[0].failed());
        assert!(previous_strain > failure);
        assert!(load_turns >= 1);
        close(
            path.solid().fracture_dissipation_j(),
            FRACTURE_WORK * band.projected_area_m2(),
            3e-9,
        );
        close(
            exact_work,
            path.solid().stored_energy_j() + path.solid().fracture_dissipation_j(),
            3e-9,
        );
        assert_eq!(path.accepted_steps(), 90);
        assert_eq!(path.trace().len(), 91);
        assert_eq!(virgin.fracture_dissipation_j(), 0.0);
    });
}

#[test]
fn explicit_reversal_unloads_the_original_damaged_solid_without_healing() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "arclength reversal");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let damaged = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let units = scales(&damaged);
        let lambda = damaged.solution().element_stress[0][0] / STRENGTH;
        let path = BrittleTetContinuation::try_from_solid(
            &damaged,
            &fixed_supports(),
            &[(3, units.force_n)],
            lambda,
            units,
            arc_budget(),
            &arc_settings(),
            cx,
        )
        .unwrap();
        let mut direction = [[0.0; 3]; 4];
        direction[1][0] = -units.displacement_m;
        let step = path
            .advance(
                &arc_settings(),
                arc_controls(),
                Some(BrittleTetArcDirection {
                    displacement_m: &direction,
                    load_factor: -0.1,
                }),
                cx,
            )
            .unwrap();
        let before = &damaged.solution().element_states[0];
        let after = &step.continuation.solid().solution().element_states[0];
        assert!(after.strain()[0] < before.strain()[0]);
        assert_eq!(after.damage(), before.damage());
        assert_eq!(after.dissipated_j_m3(), before.dissipated_j_m3());
        close(
            step.continuation.solid().solution().element_stress[0][0],
            damaged.solution().element_stress[0][0] * after.strain()[0] / before.strain()[0],
            3e-9,
        );
        assert!(step.energy.material_work_j < 0.0);
        assert_eq!(step.energy.fracture_dissipation_j, 0.0);
        assert_eq!(path.solid().solution(), damaged.solution());
    });
}

#[test]
fn work_damage_step_and_cancellation_refusals_preserve_the_same_checkpoint() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "arclength refusals");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let units = scales(&virgin);
        let fixed = fixed_supports();
        let loads = [(3, units.force_n)];
        let settings = arc_settings();
        let path = BrittleTetContinuation::try_from_solid(
            &virgin,
            &fixed,
            &loads,
            0.0,
            units,
            arc_budget(),
            &settings,
            cx,
        )
        .unwrap();
        let strict = BrittleTetArcControls {
            work: BrittleTetControls {
                maximum_endpoint_work_defect_j: 0.0,
                ..controls()
            },
            ..arc_controls()
        };
        assert!(matches!(
            path.advance(&settings, strict, None, cx),
            Err(BrittleTetContinuationError::Solid(
                BrittleTetError::EnergyGate { .. }
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                path.advance(&settings, arc_controls(), None, cancelled),
                Err(BrittleTetContinuationError::Solid(
                    BrittleTetError::Cancelled
                ))
            ));
        });
        let narrow = BrittleTetContinuation::try_from_solid(
            &virgin,
            &fixed,
            &loads,
            0.0,
            units,
            BrittleTetArcBudget {
                maximum_operator_evaluations: 1,
                ..arc_budget()
            },
            &settings,
            cx,
        )
        .unwrap();
        assert!(matches!(
            narrow.advance(&settings, arc_controls(), None, cx),
            Err(BrittleTetContinuationError::OperatorBudget { evaluations: 1 })
        ));
        let once = BrittleTetContinuation::try_from_solid(
            &virgin,
            &fixed,
            &loads,
            0.0,
            units,
            BrittleTetArcBudget {
                maximum_steps: 1,
                ..arc_budget()
            },
            &settings,
            cx,
        )
        .unwrap()
        .advance(&settings, arc_controls(), None, cx)
        .unwrap()
        .continuation;
        assert!(matches!(
            once.advance(&settings, arc_controls(), None, cx),
            Err(BrittleTetContinuationError::StepBudget)
        ));
        assert_eq!(path.solid().solution(), virgin.solution());
        assert_eq!(path.accepted_steps(), 0);
        assert_eq!(path.trace(), &[(0.0, 0.0)]);

        let damaged = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let damaged_path = BrittleTetContinuation::try_from_solid(
            &damaged,
            &fixed,
            &loads,
            damaged.solution().element_stress[0][0] / STRENGTH,
            units,
            arc_budget(),
            &settings,
            cx,
        )
        .unwrap();
        let mut direction = [[0.0; 3]; 4];
        direction[1][0] = units.displacement_m;
        let no_growth = BrittleTetArcControls {
            maximum_damage_increment: 0.0,
            ..arc_controls()
        };
        assert!(matches!(
            damaged_path.advance(
                &settings,
                no_growth,
                Some(BrittleTetArcDirection {
                    displacement_m: &direction,
                    load_factor: -0.01,
                }),
                cx
            ),
            Err(BrittleTetContinuationError::DamageIncrement { .. })
        ));
        assert_eq!(damaged_path.solid().solution(), damaged.solution());
    });
}

#[test]
fn initial_equilibrium_supports_loads_and_physical_scales_are_admitted() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "arclength initial admission");
        let materials = [material(&source, &geometry, 0)];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let units = scales(&virgin);
        let fixed = fixed_supports();
        let loads = [(3, units.force_n)];
        let settings = arc_settings();
        assert!(matches!(
            BrittleTetContinuation::try_from_solid(
                &virgin,
                &fixed,
                &loads,
                0.5,
                units,
                arc_budget(),
                &settings,
                cx,
            ),
            Err(BrittleTetContinuationError::InitialEquilibrium { .. })
        ));
        let mut duplicate = fixed.clone();
        duplicate.push(fixed[0]);
        assert!(
            BrittleTetContinuation::try_from_solid(
                &virgin,
                &duplicate,
                &loads,
                0.0,
                units,
                arc_budget(),
                &settings,
                cx,
            )
            .is_err()
        );
        assert!(
            BrittleTetContinuation::try_from_solid(
                &virgin,
                &fixed,
                &[(0, units.force_n)],
                0.0,
                units,
                arc_budget(),
                &settings,
                cx,
            )
            .is_err()
        );
        assert!(
            BrittleTetContinuation::try_from_solid(
                &virgin,
                &fixed,
                &loads,
                0.0,
                BrittleTetPathScales {
                    displacement_m: 0.0,
                    ..units
                },
                arc_budget(),
                &settings,
                cx,
            )
            .is_err()
        );
        let path = BrittleTetContinuation::try_from_solid(
            &virgin,
            &fixed,
            &loads,
            0.0,
            units,
            arc_budget(),
            &settings,
            cx,
        )
        .unwrap();
        let mut moving_support = [[0.0; 3]; 4];
        moving_support[0][0] = units.displacement_m;
        assert!(
            path.advance(
                &settings,
                arc_controls(),
                Some(BrittleTetArcDirection {
                    displacement_m: &moving_support,
                    load_factor: 1.0,
                }),
                cx
            )
            .is_err()
        );
    });
}
