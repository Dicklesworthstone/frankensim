//! Prescribed angular rates retain SI units and the original accepted clock.
use super::*;

fn spinning_model<'a>(
    liquid: &'a ResolvedLiquidState,
    solid: &'a ResolvedBody,
    rate: f64,
) -> ResolvedImmersedPlanes<'a, ResolvedLiquidState, 1> {
    let mut state = initial()[0];
    state.velocity_m_s = [0.0; 3];
    let mut surface = plane();
    surface.point_m[0] = state.center_m[0] - RADIUS - 0.005 * RADIUS;
    surface.tangential_gaps_m = Some([0.01 * RADIUS, 0.0001 * RADIUS]);
    Model::spheres(
        liquid,
        domain(),
        [SubmergedSphere::new(solid, RADIUS).unwrap()],
        [state],
        [0.0; 3],
        film(),
        &CancelGate::new(),
    )
    .unwrap()
    .with_spinning_prescribed_planes(&[surface], &[rate], &CancelGate::new())
    .unwrap()
}

#[test]
fn physical_plane_spin_and_local_velocity_use_the_original_time_and_length_scales() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let rate = 0.05;
    let model = spinning_model(&liquid, &solid, rate);
    assert_eq!(model.source_plane_axial_spins_rad_s(), &[rate]);
    assert_eq!(model.planes_lattice()[0].axial_spin(), rate * DT);
    close(model.plane_axial_spin_rad_s(0).unwrap(), rate, 2e-17);
    let source = model.plane_si(0).unwrap();
    let point = [source.point_m[0], 0.014, -0.017];
    let actual = model.plane_velocity_at_si(0, point).unwrap();
    let omega = vector(source.normal).scale(rate);
    let expected =
        vector(source.velocity_m_s) + omega.cross(vector(point) - vector(source.point_m));
    for (a, e) in actual.into_iter().zip(array(expected)) {
        close(a, e, 2e-18);
    }
    assert!(model.plane_axial_spin_rad_s(1).is_err());
    assert!(model.plane_velocity_at_si(1, point).is_err());
    assert!(model.plane_velocity_at_si(0, [f64::NAN; 3]).is_err());
    assert!(model.plane_velocity_at_si(0, [0.0, 1e6, 0.0]).is_err());
    let stationary = super::model(&liquid, &solid);
    assert_eq!(stationary.source_plane_axial_spins_rad_s(), &[0.0]);
    assert_eq!(stationary.plane_axial_spin_rad_s(0).unwrap(), 0.0);
}

#[test]
fn rotating_plane_reaction_retains_si_work_and_common_clock_through_material_rebind() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let solid = body(800.0);
    let rate = 0.005;
    let mut actual = spinning_model(&first, &solid, rate);
    let mut reference = spinning_model(&first, &solid, rate);
    let origin = [0.011, -0.013, 0.009];
    let gate = CancelGate::new();
    let control = DryBodyCouplingConfig3::default();
    let before = actual.assembly().state_si(0).unwrap();
    let fluid_before = full_fluid_momentum(actual.assembly(), origin);
    let initial_plane = actual.plane_si(0).unwrap();
    let source = actual.source_planes_si()[0];
    let report = actual
        .step_forces_si([[0.0; 6]], origin, control, &gate)
        .unwrap();
    assert_eq!(
        report,
        reference
            .step_forces_si([[0.0; 6]], origin, control, &gate)
            .unwrap()
    );
    check_plane_balance(&actual, before, fluid_before, origin, &report);
    let assembly = actual.assembly();
    let reaction = &report.planes[0];
    let impulse = vector(reaction.impulse.map(|p| assembly.impulse_ns(p).unwrap()));
    let angular = vector(
        reaction
            .angular_impulse
            .map(|l| assembly.angular_impulse_nms(l).unwrap()),
    );
    let midpoint_axis =
        vector(initial_plane.point_m) + vector(initial_plane.velocity_m_s).scale(0.5 * DT);
    let axis_torque = angular - (midpoint_axis - vector(origin)).cross(impulse);
    let work = impulse.dot(vector(initial_plane.velocity_m_s))
        + vector(initial_plane.normal).scale(rate).dot(axis_torque);
    close(assembly.energy_j(reaction.work).unwrap(), work, 1e-19);
    assert!(reaction.work < 0.0);
    assert!(
        report
            .motion
            .native
            .coupling
            .lubrication
            .active_tangential_plane_pairs
            > 0
    );
    assert_eq!(actual.source_planes_si(), &[source]);
    assert_eq!(actual.source_plane_axial_spins_rad_s(), &[rate]);
    close(actual.plane_axial_spin_rad_s(0).unwrap(), rate, 2e-18);
    let populations = fields(actual.assembly());
    let bodies = *actual.assembly().bodies_lattice();
    let planes = actual.planes_lattice().to_vec();
    actual.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert_eq!(actual.planes_lattice(), planes);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert_eq!(fields(actual.assembly()), populations);
    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        actual
            .step_forces_si([[0.0; 6]], origin, control, &stopped)
            .is_err()
    );
    assert_eq!(actual.planes_lattice(), planes);
    assert_eq!(actual.assembly().accepted_steps(), 1);
    assert_eq!(
        actual
            .step_forces_si([[0.0; 6]], origin, control, &gate)
            .unwrap(),
        reference
            .step_forces_si([[0.0; 6]], origin, control, &gate)
            .unwrap()
    );
    assert_eq!(fields(actual.assembly()), fields(reference.assembly()));
    close(actual.plane_axial_spin_rad_s(0).unwrap(), rate, 2e-18);
    assert_eq!(actual.assembly().accepted_steps(), 2);
}

#[test]
fn zero_spin_preserves_prescribed_plane_results_and_invalid_rates_refuse_attachment() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let build = || {
        let mut state = initial()[0];
        state.velocity_m_s = [0.0; 3];
        Model::spheres(
            &liquid,
            domain(),
            [SubmergedSphere::new(&solid, RADIUS).unwrap()],
            [state],
            [0.0; 3],
            film(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let gate = CancelGate::new();
    assert!(
        build()
            .with_spinning_prescribed_planes(&[plane()], &[], &gate)
            .is_err()
    );
    for rate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            build()
                .with_spinning_prescribed_planes(&[plane()], &[rate], &gate)
                .is_err()
        );
    }
    let mut actual = build()
        .with_spinning_prescribed_planes(&[plane()], &[0.0], &gate)
        .unwrap();
    let mut original = build().with_prescribed_planes(&[plane()], &gate).unwrap();
    assert_eq!(actual.planes_lattice(), original.planes_lattice());
    let control = DryBodyCouplingConfig3::default();
    assert_eq!(
        actual
            .step_forces_si([[0.0; 6]], [0.0; 3], control, &gate)
            .unwrap(),
        original
            .step_forces_si([[0.0; 6]], [0.0; 3], control, &gate)
            .unwrap()
    );
    assert_eq!(fields(actual.assembly()), fields(original.assembly()));
    assert_eq!(
        actual.assembly().bodies_lattice(),
        original.assembly().bodies_lattice()
    );
}
