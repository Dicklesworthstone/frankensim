use super::*;
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurfaceImpulseCell3};
use fs_mbd::{Pose, UnitQuaternion};

fn controls() -> LiquidRigidDragControls {
    LiquidRigidDragControls {
        maximum_cells: 64,
        maximum_stencil_cells: 16,
        maximum_speed_m_s: 0.12,
        maximum_momentum_defect_kg_m_s: 1e-11,
        maximum_angular_momentum_defect_kg_m2_s: 1e-11,
        maximum_work_defect_j: 1e-11,
    }
}
fn units() -> LiquidMomentumUnits3 {
    LiquidMomentumUnits3 {
        cell_size_m: 0.25,
        velocity_unit_m_s: 1.0,
        mass_unit_kg: 2.0,
        grid_origin_world_m: Vec3::new(0.2, -0.3, 0.1),
    }
}
fn surface() -> FreeSurface3 {
    FreeSurface3::new(
        4,
        4,
        4,
        0.9,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, _, _| x < 3,
    )
    .unwrap()
}
fn near(x: f64, y: f64, tolerance: f64) {
    assert!(
        (x - y).abs() <= tolerance,
        "{x:e} != {y:e} (tol={tolerance:e})"
    );
}
fn near_vec(x: Vec3, y: Vec3, tolerance: f64) {
    for (a, b) in a(x).into_iter().zip(a(y)) {
        near(a, b, tolerance);
    }
}
fn samples(
    liquid: &mut FreeSurface3,
    indices: &[usize],
    units: LiquidMomentumUnits3,
) -> Vec<FreeSurfaceImpulseCell3> {
    liquid
        .prepare_impulse_stencil(
            indices,
            controls().native(units, units.admit().unwrap()).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
        .cells()
        .collect()
}
fn body() -> (RigidBodyState, MassProperties) {
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.4, -0.2), 0.43).unwrap();
    let properties = MassProperties::from_principal_frame(
        3.2,
        Vec3::ZERO,
        Vec3::new(0.14, 0.21, 0.26),
        principal,
    )
    .unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(-0.3, 1.0, 0.2), 0.71).unwrap();
    let pose = Pose::new(Vec3::new(0.71, 0.13, 0.52), orientation).unwrap();
    let omega = orientation.rotate_world_to_body(Vec3::new(0.012, -0.017, 0.02));
    let body = RigidBodyState::new(
        pose,
        Vec3::new(0.03, -0.01, 0.01).scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap();
    (body, properties)
}
fn total(
    cells: &[FreeSurfaceImpulseCell3],
    body: RigidBodyState,
    properties: MassProperties,
    u: LiquidMomentumUnits3,
) -> (Vec3, Vec3, f64) {
    let mut p = body.linear_momentum_world();
    let mut l = body
        .pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
        .add(body.pose().position_world().cross(p));
    let mut energy = body.kinetic_energy(properties).unwrap();
    for c in cells {
        let mass = c.mass * u.mass_unit_kg;
        let velocity = v(c.velocity).scale(u.velocity_unit_m_s);
        let pos = u
            .grid_origin_world_m
            .add(v(c.position).scale(u.cell_size_m));
        p = p.add(velocity.scale(mass));
        l = l.add(pos.cross(velocity.scale(mass)));
        energy += 0.5 * mass * velocity.dot(velocity);
    }
    (p, l, energy)
}

#[test]
fn partial_fill_recoils_with_retained_mass_against_finite_body() {
    let mut liquid = surface();
    let u = units();
    let ids = [2];
    let sample = samples(&mut liquid, &ids, u)[0];
    near(sample.mass, 0.5, 2e-16);
    let m = sample.mass * u.mass_unit_kg;
    let properties = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap();
    let center = u
        .grid_origin_world_m
        .add(v(sample.position).scale(u.cell_size_m));
    let pose = Pose::new(center, UnitQuaternion::IDENTITY).unwrap();
    let mut rigid = RigidBodyState::new(pose, Vec3::new(0.18, 0.0, 0.0), Vec3::ZERO).unwrap();
    let initial = rigid;
    let velocity = initial.center_of_mass_velocity_world(properties).unwrap().x;
    let rate = 3.0;
    let duration = 0.2;
    let j = 3.0 * m / (3.0 + m) * velocity * -(-rate * duration * (1.0 + m / 3.0)).exp_m1();
    let report = apply_rigid_drag_native(
        &mut liquid,
        &mut rigid,
        properties,
        &ids,
        u,
        rate,
        duration,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    near(report.liquid_impulse_kg_m_s.x, j, 2e-17);
    near(
        rigid.linear_momentum_world().x,
        initial.linear_momentum_world().x - j,
        5e-17,
    );
    let end = samples(&mut liquid, &ids, u)[0];
    near(end.velocity[0] * u.velocity_unit_m_s, j / m, 3e-17);
    assert_eq!(rigid.pose(), initial.pose());
    assert_eq!(liquid.grid().steps(), 0);
    let expected_loss = 0.5 * 3.0 * m / (3.0 + m)
        * velocity
        * velocity
        * (1.0 - (-2.0 * rate * duration * (1.0 + m / 3.0)).exp());
    near(report.drag_dissipation_j, expected_loss, 5e-18);
}

#[test]
fn rotated_full_inertia_closes_linear_angular_and_kinetic_accounts() {
    let ids = [5, 6, 17, 18, 21, 22, 33, 34, 37, 38];
    let mut liquid = surface();
    let u = units();
    let (mut rigid, properties) = body();
    let initial = rigid;
    let original = samples(&mut liquid, &ids, u);
    let before = total(&original, rigid, properties, u);
    let ledger = liquid.ledger_mass();
    let epoch = liquid.geometry_epoch();
    let clock = liquid.grid().steps();
    let original_fill: Vec<_> = ids
        .iter()
        .map(|i| {
            liquid
                .sample_at(i % 4, (i / 4) % 4, i / 16)
                .unwrap()
                .fill_fraction
                .to_bits()
        })
        .collect();
    let report = apply_rigid_drag_native(
        &mut liquid,
        &mut rigid,
        properties,
        &ids,
        u,
        5.0,
        0.2,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let end = samples(&mut liquid, &ids, u);
    let after = total(&end, rigid, properties, u);
    near_vec(before.0, after.0, 1e-14);
    near_vec(before.1, after.1, 1e-14);
    near(before.2 - after.2, report.drag_dissipation_j, 3e-17);
    assert!(report.drag_dissipation_j > 0.0);
    assert!(
        report
            .liquid_angular_impulse_kg_m2_s
            .dot(report.liquid_angular_impulse_kg_m2_s)
            > 1e-8
    );
    near(report.kinetic_balance_residual_j, 0.0, 2e-17);
    assert_eq!(initial.pose(), rigid.pose());
    assert_eq!(ledger.to_bits(), liquid.ledger_mass().to_bits());
    assert_eq!(
        (epoch, clock),
        (liquid.geometry_epoch(), liquid.grid().steps())
    );
    for ((i, fill), (a, b)) in ids.iter().zip(original_fill).zip(original.iter().zip(end)) {
        assert_eq!(
            fill,
            liquid
                .sample_at(i % 4, (i / 4) % 4, i / 16)
                .unwrap()
                .fill_fraction
                .to_bits()
        );
        assert_eq!(a.mass.to_bits(), b.mass.to_bits());
    }
}

#[test]
fn native_exact_source_composes_and_stiff_limit_has_no_residual_slip() {
    let ids = [5, 6, 17, 18, 21, 22, 33, 34, 37, 38];
    let u = units();
    let (original, properties) = body();
    let mut full_liquid = surface();
    let mut twice_liquid = surface();
    let mut full = original;
    let mut twice = original;
    let gate = CancelGate::new_clock_free();
    let report = apply_rigid_drag_native(
        &mut full_liquid,
        &mut full,
        properties,
        &ids,
        u,
        7.0,
        0.3,
        controls(),
        &gate,
    )
    .unwrap();
    let first = apply_rigid_drag_native(
        &mut twice_liquid,
        &mut twice,
        properties,
        &ids,
        u,
        7.0,
        0.1,
        controls(),
        &gate,
    )
    .unwrap();
    let second = apply_rigid_drag_native(
        &mut twice_liquid,
        &mut twice,
        properties,
        &ids,
        u,
        7.0,
        0.2,
        controls(),
        &gate,
    )
    .unwrap();
    near_vec(
        full.linear_momentum_world(),
        twice.linear_momentum_world(),
        5e-16,
    );
    near_vec(
        full.angular_momentum_body(),
        twice.angular_momentum_body(),
        5e-16,
    );
    for (a, b) in samples(&mut full_liquid, &ids, u)
        .iter()
        .zip(samples(&mut twice_liquid, &ids, u))
    {
        near_vec(v(a.velocity), v(b.velocity), 5e-16);
    }
    near(
        report.drag_dissipation_j,
        first.drag_dissipation_j + second.drag_dissipation_j,
        4e-17,
    );
    let mut stiff_liquid = surface();
    let mut stiff = original;
    let report = apply_rigid_drag_native(
        &mut stiff_liquid,
        &mut stiff,
        properties,
        &ids,
        u,
        1e300,
        1e100,
        controls(),
        &gate,
    )
    .unwrap();
    assert!(report.drag_dissipation_j > 0.0);
    let com = stiff.center_of_mass_velocity_world(properties).unwrap();
    let omega = stiff.pose().orientation().rotate_body_to_world(
        properties
            .angular_velocity_body_checked(stiff.angular_momentum_body())
            .unwrap(),
    );
    for cell in samples(&mut stiff_liquid, &ids, u) {
        let arm = u
            .grid_origin_world_m
            .add(v(cell.position).scale(u.cell_size_m))
            .sub(stiff.pose().position_world());
        near_vec(
            v(cell.velocity).scale(u.velocity_unit_m_s),
            com.add(omega.cross(arm)),
            2e-16,
        );
    }
}

#[test]
fn common_uniform_velocity_changes_neither_exchange_nor_dissipation() {
    let ids = [5, 6, 17, 18, 21, 22, 33, 34, 37, 38];
    let u = units();
    let (original, properties) = body();
    let mut base_liquid = surface();
    let mut boosted_liquid = surface();
    let mut base = original;
    let boost = Vec3::new(0.025, -0.018, 0.007);
    let mut boosted = RigidBodyState::new(
        original.pose(),
        original
            .linear_momentum_world()
            .add(boost.scale(properties.mass())),
        original.angular_momentum_body(),
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    let stage = boosted_liquid
        .prepare_impulse_stencil(
            &ids,
            controls().native(u, u.admit().unwrap()).unwrap(),
            &gate,
        )
        .unwrap();
    let impulses: Vec<_> = stage
        .cells()
        .map(|c| a(boost.scale(c.mass / u.velocity_unit_m_s)))
        .collect();
    let _ = stage
        .with_impulses(&impulses, &gate)
        .unwrap()
        .commit(&gate)
        .unwrap();
    let a = apply_rigid_drag_native(
        &mut base_liquid,
        &mut base,
        properties,
        &ids,
        u,
        4.0,
        0.2,
        controls(),
        &gate,
    )
    .unwrap();
    let b = apply_rigid_drag_native(
        &mut boosted_liquid,
        &mut boosted,
        properties,
        &ids,
        u,
        4.0,
        0.2,
        controls(),
        &gate,
    )
    .unwrap();
    near_vec(a.liquid_impulse_kg_m_s, b.liquid_impulse_kg_m_s, 2e-16);
    near_vec(
        a.liquid_angular_impulse_kg_m2_s,
        b.liquid_angular_impulse_kg_m2_s,
        1e-16,
    );
    near(a.drag_dissipation_j, b.drag_dissipation_j, 3e-17);
}

#[test]
fn zero_source_is_identity_and_late_refusal_preserves_both_owners() {
    let ids = [5, 6, 17, 18, 21, 22, 33, 34, 37, 38];
    let u = units();
    let (mut rigid, properties) = body();
    let original = rigid;
    let mut liquid = surface();
    let bits = liquid.grid().state_bits();
    let gate = CancelGate::new_clock_free();
    let report = apply_rigid_drag_native(
        &mut liquid,
        &mut rigid,
        properties,
        &ids,
        u,
        0.0,
        0.2,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(report.drag_dissipation_j, 0.0);
    assert_eq!(original, rigid);
    assert_eq!(bits, liquid.grid().state_bits());
    let mut strict = controls();
    strict.maximum_angular_momentum_defect_kg_m2_s = 0.0;
    let error = apply_rigid_drag_native(
        &mut liquid,
        &mut rigid,
        properties,
        &ids,
        u,
        3.7,
        0.23,
        strict,
        &gate,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            LiquidMomentumError::Balance {
                account: "angular momentum and exact source response",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(original, rigid);
    assert_eq!(bits, liquid.grid().state_bits());
    gate.request();
    assert!(matches!(
        apply_rigid_drag_native(
            &mut liquid,
            &mut rigid,
            properties,
            &ids,
            u,
            3.7,
            0.23,
            controls(),
            &gate
        ),
        Err(LiquidMomentumError::Cancelled)
    ));
    assert_eq!(original, rigid);
    assert_eq!(bits, liquid.grid().state_bits());
}
