use crate::cohesive_contact::objective::source::SourceBoundObjectiveCohesivePoint;
use crate::cohesive_contact::objective::{ObjectiveCohesiveConfig, ObjectiveCohesivePoint};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_mbd::{UnitQuaternion, Vec3};
use fs_qty::{Dims, Pressure, Temperature};

const AREA: f64 = 0.002;

fn config() -> ObjectiveCohesiveConfig {
    ObjectiveCohesiveConfig {
        maximum_cayley_increment: 8.0,
        kinematic_tolerance_m: 2e-14,
        material_energy_tolerance_j_m2: 1e-8,
        spatial_work_tolerance_j: 1e-10,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    for (a, b) in [a.x, a.y, a.z].into_iter().zip([b.x, b.y, b.z]) {
        close(a, b, tolerance);
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn fixture(history: &str) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (COHESIVE_NORMAL_STIFFNESS_PROPERTY, stiffness, 1e9, 5e8),
        (COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY, stiffness, 4e8, 3e8),
        (COHESIVE_NORMAL_STRENGTH_PROPERTY, Pressure::DIMS, 2e5, 8e4),
        (
            COHESIVE_TANGENTIAL_STRENGTH_PROPERTY,
            Pressure::DIMS,
            1.6e5,
            8e4,
        ),
        (
            COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            50.0,
        ),
        (
            COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            400.0,
            200.0,
        ),
        (COHESIVE_BK_EXPONENT_PROPERTY, Dims::NONE, 2.0, 1.5),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, 2e9, 3e9),
        (
            COHESIVE_NORMAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            2e5,
            4e5,
        ),
        (
            COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            3e5,
            6e5,
        ),
        // Both families are present so a caller cannot silently change the
        // original envelope merely because another family resolves successfully.
        (COHESIVE_STIFFNESS_PROPERTY, stiffness, 1e9, 5e8),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5, 8e4),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            50.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, cold), (340.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal mixed-mode interface; no measured calibration"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured objective thermal fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit original bond".into(),
            third_body: None,
            environment: "prescribed interface temperature".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(
    card: &InterfaceSystemCard,
    temperature: f64,
    envelope: CohesiveEnvelope,
) -> ResolvedBilinearCohesive {
    resolve_viscous_bilinear_cohesive(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        envelope,
        0.02,
        5.0,
    )
    .unwrap()
}

fn bound(
    source: &ResolvedBilinearCohesive,
    jump: [f64; 3],
    frame: UnitQuaternion,
) -> SourceBoundObjectiveCohesivePoint<'_> {
    let material = source
        .law()
        .bonded_point()
        .trial_over_time(jump, 0.02, config().material_energy_tolerance_j_m2)
        .unwrap()
        .equilibrium
        .point;
    let point = ObjectiveCohesivePoint::new(material, frame, Vec3::new(0.003, -0.002, 0.004), AREA)
        .unwrap();
    SourceBoundObjectiveCohesivePoint::new(
        point,
        source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
    )
    .unwrap()
}

#[test]
fn fixed_jump_temperature_preserves_exact_geometry_mode_and_irreversible_history() {
    let card = fixture("original thermal bond");
    let cold = resolve(&card, 280.0, CohesiveEnvelope::FrozenOnsetMixedMode);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), 1.3).unwrap();
    let loaded = bound(&cold, [4e-4, 5e-4, -2e-4], frame);
    let separation = frame.rotate_body_to_world(
        loaded
            .numerical()
            .reference_offset_m()
            .add(vector([-2e-4, 5e-4, -2e-4])),
    );
    let separation = Vec3::new(separation.x.next_up(), separation.y, separation.z);
    let point = loaded
        .trial(frame, separation, 0.02, config())
        .unwrap()
        .point();
    let saved = point;
    let old = point.numerical().material();
    let original_mode = old.mixed_mode_history().unwrap();
    let step = point.transfer_temperature(&hot, "T", config()).unwrap();
    let next = step.point();
    assert_eq!(next.numerical().frame(), point.numerical().frame());
    assert_eq!(
        next.numerical().reference_offset_m(),
        point.numerical().reference_offset_m()
    );
    assert_eq!(next.numerical().separation_world_m(), separation);
    assert_eq!(next.numerical().area_m2(), AREA);
    assert_eq!(next.numerical().material().jump_m(), old.jump_m());
    assert!(core::ptr::eq(next.source(), &hot));
    assert!(core::ptr::eq(next.numerical().material().law(), hot.law()));
    let mode = next.numerical().material().mixed_mode_history().unwrap();
    assert_eq!(mode.onset_direction, original_mode.onset_direction);
    assert_eq!(
        mode.energy_shear_fraction,
        original_mode.energy_shear_fraction
    );
    assert!(next.numerical().material().damage() >= old.damage());
    let g = old.jump_m();
    let parameter_work = 0.5
        * AREA
        * ((1.0 - old.damage())
            * ((5e8 - 1e9) * g[0].max(0.0).powi(2) + (3e8 - 4e8) * (g[1] * g[1] + g[2] * g[2]))
            + (3e9 - 2e9) * g[0].min(0.0).powi(2));
    close(step.parameter_work_j, parameter_work, 2e-12);
    close(
        step.parameter_work_j,
        AREA * step.material.parameter_work_j_m2,
        2e-14,
    );
    close(
        step.damage_dissipated_j,
        AREA * step.material.damage_dissipated_j_m2,
        2e-14,
    );
    close(
        step.stored_change_j,
        AREA * (next.numerical().material().stored_energy_j_m2() - old.stored_energy_j_m2()),
        2e-14,
    );
    close(
        step.parameter_work_j - step.stored_change_j - step.damage_dissipated_j,
        0.0,
        2e-12,
    );
    close(step.energy_residual_j, 0.0, 2e-12);
    let cooled = next
        .transfer_temperature(&cold, "T", config())
        .unwrap()
        .point();
    assert!(cooled.numerical().material().damage() >= next.numerical().material().damage());
    close(
        cooled.numerical().material().dissipated_j_m2(),
        old.dissipated_j_m2(),
        2e-12,
    );
    assert_eq!(
        cooled
            .numerical()
            .material()
            .mixed_mode_history()
            .unwrap()
            .onset_direction,
        original_mode.onset_direction
    );
    assert_eq!(point, saved);
}

#[test]
fn rotating_thermal_mixed_load_rupture_and_reclosure_close_actual_world_work() {
    let card = fixture("original thermal bond");
    let cold = resolve(&card, 280.0, CohesiveEnvelope::FrozenOnsetMixedMode);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.45).unwrap();
    let mut point = bound(&cold, [4e-4, 5e-4, -2e-4], frame);
    let initial_storage = AREA * point.numerical().material().stored_energy_j_m2();
    let initial_fracture = AREA * point.numerical().material().dissipated_j_m2();
    let original_mode = point.numerical().material().mixed_mode_history().unwrap();
    let (mut work, mut parameter_work, mut fracture, mut viscosity) = (0.0, 0.0, 0.0, 0.0);
    for (index, (source, jump)) in [
        (&hot, [8e-4, 7e-4, -6e-4]),
        (&hot, [0.004, 0.008, -0.003]),
        (&cold, [-2e-4, 0.001, -2e-4]),
    ]
    .into_iter()
    .enumerate()
    {
        let frame = point
            .numerical()
            .frame()
            .right_exp(Vec3::new(0.4, -0.3, 0.2))
            .unwrap();
        let separation =
            frame.rotate_body_to_world(point.numerical().reference_offset_m().add(vector(jump)));
        let old = point;
        let thermal = point.transfer_temperature(source, "T", config()).unwrap();
        let manual = thermal
            .point()
            .trial(frame, separation, 0.02, config())
            .unwrap();
        let step = point
            .trial_with_temperature(source, "T", frame, separation, 0.02, config())
            .unwrap();
        assert_eq!(step.numerical(), manual.numerical());
        assert_eq!(step.temperature().material, thermal.material);
        let next = step.point();
        assert!(core::ptr::eq(next.source(), source));
        close(
            step.total_stored_change_j,
            AREA * (next.numerical().material().stored_energy_j_m2()
                - old.numerical().material().stored_energy_j_m2()),
            2e-12,
        );
        close(
            step.total_fracture_dissipated_j,
            AREA * (next.numerical().material().dissipated_j_m2()
                - old.numerical().material().dissipated_j_m2()),
            2e-12,
        );
        close(
            step.viscous_dissipated_j,
            AREA * step.numerical().material.viscous_dissipated_j_m2,
            2e-14,
        );
        let wrench = step.numerical().discrete_wrench;
        close_vec(
            wrench.force_on_a_n.add(wrench.force_on_b_n),
            Vec3::ZERO,
            1e-15,
        );
        // Actual attachment translation and the work-conjugate frame rotation.
        let a0 = Vec3::new(0.3, -0.7, 1.2);
        let a1 = a0.add(Vec3::new(0.01, -0.02, 0.03));
        let b0 = a0.add(old.numerical().separation_world_m());
        let b1 = a1.add(next.numerical().separation_world_m());
        let body_work = wrench.force_on_a_n.dot(a1.sub(a0))
            + wrench.force_on_b_n.dot(b1.sub(b0))
            + wrench
                .frame_couple_n_m
                .dot(step.numerical().cayley_increment_world);
        let moment = a0
            .scale(0.5)
            .add(a1.scale(0.5))
            .cross(wrench.force_on_a_n)
            .add(b0.scale(0.5).add(b1.scale(0.5)).cross(wrench.force_on_b_n))
            .add(wrench.frame_couple_n_m);
        close_vec(moment, Vec3::ZERO, 3e-12);
        close(-body_work, step.numerical().spatial_work_j, 3e-12);
        close(
            -body_work + step.temperature().parameter_work_j,
            step.total_stored_change_j
                + step.total_fracture_dissipated_j
                + step.viscous_dissipated_j,
            3e-12,
        );
        close(step.energy_residual_j, 0.0, 2e-12);
        let mode = next.numerical().material().mixed_mode_history().unwrap();
        assert_eq!(mode.onset_direction, original_mode.onset_direction);
        assert_eq!(
            mode.energy_shear_fraction,
            original_mode.energy_shear_fraction
        );
        if index == 0 {
            assert!(step.temperature().damage_dissipated_j > 0.0);
            assert!(step.total_fracture_dissipated_j > step.temperature().damage_dissipated_j);
            assert!(step.viscous_dissipated_j > 0.0);
        } else {
            assert!(next.numerical().material().failed());
        }
        if index == 2 {
            assert_eq!(step.viscous_dissipated_j, 0.0);
            assert!(next.numerical().material().traction_pa()[0] < 0.0);
            assert_eq!(next.numerical().material().traction_pa()[1], 0.0);
        }
        work += step.numerical().spatial_work_j;
        parameter_work += step.temperature().parameter_work_j;
        fracture += step.total_fracture_dissipated_j;
        viscosity += step.viscous_dissipated_j;
        point = next;
    }
    close(
        fracture,
        AREA * point.numerical().material().dissipated_j_m2() - initial_fracture,
        2e-12,
    );
    close(
        work + parameter_work,
        AREA * point.numerical().material().stored_energy_j_m2() - initial_storage
            + fracture
            + viscosity,
        3e-12,
    );
}

#[test]
fn thermal_mixed_world_trial_covaries_under_a_common_observer_rotation() {
    let card = fixture("original thermal bond");
    let cold = resolve(&card, 280.0, CohesiveEnvelope::FrozenOnsetMixedMode);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let point = bound(&cold, [4e-4, 5e-4, -2e-4], UnitQuaternion::IDENTITY);
    let observer = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -3.0, 2.0), 2.4).unwrap();
    let rotated = point.transport_frame(observer).unwrap();
    let rotation = Vec3::new(0.7, -0.4, 0.2);
    let a_frame = UnitQuaternion::IDENTITY.right_exp(rotation).unwrap();
    let b_frame = observer.right_exp(rotation).unwrap();
    let local_separation = point
        .numerical()
        .reference_offset_m()
        .add(vector([7e-4, -5e-4, 6e-4]));
    let a = point
        .trial_with_temperature(
            &hot,
            "T",
            a_frame,
            a_frame.rotate_body_to_world(local_separation),
            0.02,
            config(),
        )
        .unwrap();
    let b = rotated
        .trial_with_temperature(
            &hot,
            "T",
            b_frame,
            b_frame.rotate_body_to_world(local_separation),
            0.02,
            config(),
        )
        .unwrap();
    assert_eq!(a.temperature().material, b.temperature().material);
    assert_eq!(
        a.point().numerical().material().mixed_mode_history(),
        b.point().numerical().material().mixed_mode_history()
    );
    close(
        a.point().numerical().material().damage(),
        b.point().numerical().material().damage(),
        2e-12,
    );
    for (left, right) in [
        (a.numerical().discrete_wrench, b.numerical().discrete_wrench),
        (a.numerical().endpoint_wrench, b.numerical().endpoint_wrench),
    ] {
        close_vec(
            right.force_on_a_n,
            observer.rotate_body_to_world(left.force_on_a_n),
            2e-12,
        );
        close_vec(
            right.force_on_b_n,
            observer.rotate_body_to_world(left.force_on_b_n),
            2e-12,
        );
        close_vec(
            right.frame_couple_n_m,
            observer.rotate_body_to_world(left.frame_couple_n_m),
            2e-12,
        );
    }
    close_vec(
        b.numerical().cayley_increment_world,
        observer.rotate_body_to_world(a.numerical().cayley_increment_world),
        2e-14,
    );
    for (a, b) in [
        (a.numerical().spatial_work_j, b.numerical().spatial_work_j),
        (
            a.temperature().parameter_work_j,
            b.temperature().parameter_work_j,
        ),
        (a.total_stored_change_j, b.total_stored_change_j),
        (a.total_fracture_dissipated_j, b.total_fracture_dissipated_j),
        (a.viscous_dissipated_j, b.viscous_dissipated_j),
    ] {
        close(a, b, 2e-12);
    }
    close(a.energy_residual_j, 0.0, 2e-12);
    close(b.energy_residual_j, 0.0, 2e-12);
}

#[test]
fn wrong_source_model_and_late_trial_refusals_leave_the_original_point_retryable() {
    let card = fixture("original thermal bond");
    let cold = resolve(&card, 280.0, CohesiveEnvelope::FrozenOnsetMixedMode);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let other_card = fixture("different bond preparation and original source");
    let foreign = resolve(&other_card, 340.0, CohesiveEnvelope::FrozenOnsetMixedMode);
    let wrong_model = resolve(&card, 340.0, CohesiveEnvelope::CommonMode);
    let point = bound(&cold, [4e-4, 5e-4, -2e-4], UnitQuaternion::IDENTITY);
    let saved = point;
    let frame = UnitQuaternion::IDENTITY
        .right_exp(Vec3::new(0.4, -0.3, 0.2))
        .unwrap();
    let separation = frame.rotate_body_to_world(
        point
            .numerical()
            .reference_offset_m()
            .add(vector([8e-4, 7e-4, -6e-4])),
    );
    let accepted = point
        .trial_with_temperature(&hot, "T", frame, separation, 0.02, config())
        .unwrap();
    for source in [&foreign, &wrong_model] {
        assert!(point.transfer_temperature(source, "T", config()).is_err());
        assert!(
            point
                .trial_with_temperature(source, "T", frame, separation, 0.02, config())
                .is_err()
        );
    }
    for bad in [
        ObjectiveCohesiveConfig {
            maximum_cayley_increment: f64::NAN,
            ..config()
        },
        ObjectiveCohesiveConfig {
            kinematic_tolerance_m: -1.0,
            ..config()
        },
        ObjectiveCohesiveConfig {
            material_energy_tolerance_j_m2: 0.0,
            ..config()
        },
        ObjectiveCohesiveConfig {
            spatial_work_tolerance_j: 0.0,
            ..config()
        },
    ] {
        assert!(point.transfer_temperature(&hot, "T", bad).is_err());
        assert!(
            point
                .trial_with_temperature(&hot, "T", frame, separation, 0.02, bad)
                .is_err()
        );
    }
    assert!(
        point
            .transfer_temperature(&hot, "unowned_temperature", config())
            .is_err()
    );
    assert!(
        point
            .trial_with_temperature(&hot, "T", frame, separation, 0.0, config())
            .is_err()
    );
    assert!(
        point
            .trial_with_temperature(
                &hot,
                "T",
                frame,
                Vec3::new(f64::NAN, 0.0, 0.0),
                0.02,
                config()
            )
            .is_err()
    );
    let rotation_bound = ObjectiveCohesiveConfig {
        maximum_cayley_increment: 0.1,
        ..config()
    };
    // The source/temperature stage itself releases real damage, but the later
    // finite-rotation refusal must not publish that otherwise valid candidate.
    assert!(
        point
            .transfer_temperature(&hot, "T", rotation_bound)
            .unwrap()
            .damage_dissipated_j
            > 0.0
    );
    assert!(
        point
            .trial_with_temperature(&hot, "T", frame, separation, 0.02, rotation_bound)
            .is_err()
    );
    assert_eq!(point, saved);
    let retry = point
        .trial_with_temperature(&hot, "T", frame, separation, 0.02, config())
        .unwrap();
    assert_eq!(retry.point(), accepted.point());
    assert_eq!(retry.numerical(), accepted.numerical());
    assert_eq!(
        retry.total_fracture_dissipated_j,
        accepted.total_fracture_dissipated_j
    );
    assert_eq!(
        retry.temperature().parameter_work_j,
        accepted.temperature().parameter_work_j
    );
}
