//! G1: consistent-mass inertia, discovered tensile planes and original replay.
//! Manufactured elastic/softening references are implementation checks only.

use super::*;
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::Density;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicControls, BrittleTetDynamicError};
use fs_solid::brittle::solid::initiation::dynamics::{
    TensileElasticTetMotion, TensileElasticTetMotionStep, TensileTetDynamicError,
    TensileTetMotionProbe,
};

fn dynamic_source(density: f64) -> MaterialCard {
    let base = source_card(STRENGTH, "native inertial tensile initiation", true);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            value: PropertyValue::Scalar {
                value: density,
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic reference density".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn dynamic_material(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
) -> ResolvedBrittleCrackBand {
    let geometry = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_brittle_dynamic_material(
        source,
        &QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        geometry.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn dynamic_controls() -> BrittleTetDynamicControls {
    BrittleTetDynamicControls {
        work: controls(),
        kinetic_work_tolerance_j: 1e-10,
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn elastic_step<'a>(probe: TensileTetMotionProbe<'a>) -> TensileElasticTetMotionStep<'a> {
    match probe {
        TensileTetMotionProbe::Elastic(step) => step,
        TensileTetMotionProbe::Onset(_) => panic!("expected a strictly elastic native motion step"),
    }
}

#[test]
fn g1_free_elastic_motion_uses_original_consistent_mass_and_numerical_accounts() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let density = 2000.0;
        let card = dynamic_source(density);
        let materials = [dynamic_material(&card, &geometry)];
        let solid = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let v = [[0.2, -0.1, 0.05]; 4];
        let original =
            TensileElasticTetMotion::try_from_solid(&solid, &v, 1.25, 2, initiation_settings(), cx)
                .unwrap();
        let mass = density * geometry.element_volumes_m3()[0];
        let h = 1e-6;
        let force = 2.0;
        let loads: Vec<_> = (0..4).map(|node| (3 * node, force / 4.0)).collect();
        let first = elastic_step(
            original
                .probe(
                    h,
                    &[],
                    &loads,
                    IncrementSettings::default(),
                    dynamic_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let dv = h * force / mass;
        let expected = [v[0][0] + dv, v[0][1], v[0][2]];
        for (velocity, displacement) in first
            .motion
            .velocity_mps()
            .iter()
            .zip(&first.motion.solution().displacement_m)
        {
            for axis in 0..3 {
                close(velocity[axis], expected[axis], 2e-8);
                close(displacement[axis], h * expected[axis], 2e-8);
            }
        }
        close(first.motion.reference_mass_kg(), mass, 1e-13);
        assert_eq!(first.motion.reference_density_kg_m3(), &[density]);
        close(
            first.energy.inertia.numerical_loss_j,
            0.5 * mass * dv * dv,
            2e-7,
        );
        close(
            first.energy.dead_load_endpoint_work_j,
            force * h * expected[0],
            2e-8,
        );
        assert_eq!(first.energy.fracture_dissipation_j, 0.0);
        assert!(first.energy.total_balance_error_j.abs() < 1e-10);
        assert_eq!(first.motion.accepted_steps(), 1);
        assert_eq!(first.motion.remaining_steps(), 1);
        assert_eq!(original.velocity_mps(), &v);
        assert_eq!(original.time_s(), 1.25);
        let coast = elastic_step(
            first
                .motion
                .probe(
                    h,
                    &[],
                    &[],
                    IncrementSettings::default(),
                    dynamic_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        close(
            coast.motion.kinetic_energy_j(),
            first.motion.kinetic_energy_j(),
            2e-8,
        );
        assert_eq!(coast.motion.accepted_steps(), 2);
        assert_eq!(coast.motion.remaining_steps(), 0);
        assert!(matches!(
            coast.motion.probe(
                h,
                &[],
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                None,
                cx
            ),
            Err(TensileTetDynamicError::Motion(
                BrittleTetDynamicError::StepBudget { .. }
            ))
        ));
    });
}

#[test]
fn g1_inertial_tensile_plane_handoff_preserves_consumed_time_mass_momentum_and_original_request() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let density = 2000.0;
        let card = dynamic_source(density);
        let materials = [dynamic_material(&card, &geometry)];
        let solid = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let original = TensileElasticTetMotion::try_from_solid(
            &solid,
            &[[0.0; 3]; 4],
            2.0,
            2,
            initiation_settings(),
            cx,
        )
        .unwrap();
        let h = 1e-5;
        let length = 0.02;
        let volume = geometry.element_volumes_m3()[0];
        let area = volume / length;
        let modal_mass = density * volume / 10.0;
        let stiffness = YOUNG * volume / (length * length);
        let onset_strain = STRENGTH / YOUNG;
        let first_u = 0.4 * onset_strain * length;
        let first_force = (stiffness + modal_mass / (h * h)) * first_u;
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 7)
            .map(|dof| (dof, 0.0))
            .collect();
        let settings = IncrementSettings {
            max_newton_iterations: 12,
            ..IncrementSettings::default()
        };
        let accepted = elastic_step(
            original
                .probe(
                    h,
                    &supports,
                    &[(7, first_force)],
                    settings,
                    dynamic_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let motion = &accepted.motion;
        close(motion.solution().displacement_m[2][1], first_u, 2e-9);
        close(motion.velocity_mps()[2][1], first_u / h, 2e-9);
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 1);
        let before = motion.solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let force = 1.1 * STRENGTH * area;
        // Preserve duplicate free loads and an applied force on a support.
        let forces = [(7, 0.6 * force), (0, 2.0), (7, 0.4 * force)];
        let prediction = match motion
            .probe(
                h,
                &supports,
                &forces,
                settings,
                dynamic_controls(),
                None,
                cx,
            )
            .unwrap()
        {
            TensileTetMotionProbe::Onset(onset) => onset,
            TensileTetMotionProbe::Elastic(_) => panic!("expected inertial tensile crossing"),
        };
        assert_eq!(prediction.initial_solution(), &before);
        assert_eq!(prediction.requested_dt_s(), h);
        assert_eq!(
            prediction.requested_interval_s(),
            [motion.time_s(), motion.time_s() + h]
        );
        assert_eq!(prediction.prescribed_end_m(), &supports);
        assert_eq!(prediction.nodal_forces_end_n(), &forces);
        let normal = prediction.events()[0].as_ref().unwrap().native().normal();
        close(normal[1].abs(), 1.0, 2e-12);
        assert!(normal[0].abs() < 1e-12 && normal[2].abs() < 1e-12);
        assert_eq!(materials[0].law().normal(), [1.0, 0.0, 0.0]);
        let predictor_u = (force + modal_mass / (h * h) * (first_u + h * velocity[2][1]))
            / (stiffness + modal_mass / (h * h));
        close(
            prediction.elastic_predictor().displacement_m[2][1],
            predictor_u,
            2e-9,
        );
        assert!(predictor_u / length > onset_strain);

        // A smaller physical interval is solved again, not interpolated along
        // the full-step predictor. It resolves an admitted elastic endpoint.
        let short_h = 0.1 * h;
        let shorter = elastic_step(
            motion
                .probe(
                    short_h,
                    &supports,
                    &forces,
                    settings,
                    dynamic_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let short_u = (force
            + modal_mass / (short_h * short_h) * (first_u + short_h * velocity[2][1]))
            / (stiffness + modal_mass / (short_h * short_h));
        close(
            shorter.motion.solution().displacement_m[2][1],
            short_u,
            2e-9,
        );
        assert!((short_u - (first_u + 0.1 * (predictor_u - first_u))).abs() > 1e-8);
        assert_eq!(motion.solution(), &before);

        let bands = prediction.resolve_bands(cx).unwrap();
        let handoff = prediction.handoff(&bands, 1e-9, cx).unwrap();
        assert_eq!(
            handoff.motion.solid().solution().displacement_m,
            before.displacement_m
        );
        assert_eq!(handoff.motion.velocity_mps(), velocity.as_slice());
        assert_eq!(
            handoff.motion.reference_density_kg_m3(),
            motion.reference_density_kg_m3()
        );
        assert_eq!(
            handoff.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        assert_eq!(handoff.motion.kinetic_energy_j(), motion.kinetic_energy_j());
        assert_eq!(handoff.motion.time_s(), motion.time_s());
        assert_eq!(handoff.motion.accepted_steps(), 1);
        assert_eq!(handoff.motion.remaining_steps(), 1);
        assert_eq!(handoff.motion.solid().fracture_dissipation_j(), 0.0);
        close(handoff.elastic_storage_j, motion.stored_energy_j(), 2e-12);
        assert!(handoff.storage_defect_j.abs() < 1e-12);

        // Changed reference matter cannot enter through a matching elastic law.
        let foreign = [resolve_brittle_dynamic_material(
            &dynamic_source(density * 2.0),
            &QueryPoint::new()
                .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            bands[0].law().settings(),
        )
        .unwrap()];
        assert!(prediction.handoff(&foreign, 1e-9, cx).is_err());
        let strict = BrittleTetDynamicControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..dynamic_controls()
        };
        assert!(prediction.replay(&bands, 1e-9, strict, cx).is_err());
        with_cx(true, |cancelled| {
            assert!(
                prediction
                    .replay(&bands, 1e-9, dynamic_controls(), cancelled)
                    .is_err()
            );
        });
        assert_eq!(motion.solution(), &before);
        assert_eq!(motion.velocity_mps(), velocity.as_slice());
        assert_eq!(motion.accepted_steps(), 1);

        let replay = prediction
            .replay(&bands, 1e-9, dynamic_controls(), cx)
            .unwrap();
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * bands[0].law().band_width_m());
        let softening_stiffness = -area * STRENGTH / (length * (failure - onset_strain));
        let intercept = area * STRENGTH * failure / (failure - onset_strain);
        let expected_u = (force + modal_mass / (h * h) * (first_u + h * velocity[2][1])
            - intercept)
            / (modal_mass / (h * h) + softening_stiffness);
        assert!(expected_u / length > onset_strain && expected_u / length < failure);
        close(
            replay.motion.solid().solution().displacement_m[2][1],
            expected_u,
            2e-8,
        );
        close(
            replay.motion.velocity_mps()[2][1],
            (expected_u - first_u) / h,
            2e-8,
        );
        assert!(replay.motion.solid().solution().element_states[0].damage() > 0.0);
        assert!(replay.energy.fracture_dissipation_j > 0.0);
        assert!(replay.energy.inertia.numerical_loss_j > 0.0);
        assert!(replay.energy.total_balance_error_j.abs() < 1e-9);
        assert_eq!(replay.motion.accepted_steps(), 2);
        assert_eq!(replay.motion.remaining_steps(), 0);
        assert_eq!(replay.motion.time_s(), motion.time_s() + h);
        assert_eq!(
            replay.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        let repeated = prediction
            .replay(&bands, 1e-9, dynamic_controls(), cx)
            .unwrap();
        assert_eq!(
            repeated.motion.solid().solution(),
            replay.motion.solid().solution()
        );
        assert_eq!(repeated.motion.velocity_mps(), replay.motion.velocity_mps());
        assert_eq!(motion.solution(), &before);
    });
}
