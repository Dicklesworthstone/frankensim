use super::*;
use fs_mbd::tangential_friction::adhesive::jkr_torsional::geometry::RigidJkrTorsionalGeometryError;

fn spin(f: &mut Fixture, i: usize, omega: Vec3) {
    f.bodies[i] = RigidBodyState::new(
        f.bodies[i].pose(),
        f.bodies[i].linear_momentum_world(),
        f.properties[i]
            .angular_momentum_body_checked(omega)
            .unwrap(),
    )
    .unwrap();
}

fn geometry_controls(source: &SourceJkrMixedMode) -> RigidJkrTorsionalControls {
    let mut c = controls(source);
    // Actual COMs near +/-ARM and rotated markers are represented as f64.
    // The geometry gate must include their world-coordinate rounding scale;
    // the constitutive solve and physical work tolerances remain unchanged.
    let representation = 256.0 * f64::EPSILON * ARM;
    c.geometry_tolerance_m = c.geometry_tolerance_m.max(representation);
    c.displacement_tolerance_m = c.displacement_tolerance_m.max(representation);
    c
}

fn total_energy(f: &Fixture, source: &SourceJkrMixedMode, state: RigidJkrTorsionalState) -> f64 {
    let q = state.coordinates();
    source
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap()
        .reversible_energy_j
        + (0..2)
            .map(|i| f.bodies[i].kinetic_energy(f.properties[i]).unwrap())
            .sum::<f64>()
}

#[test]
fn accepted_axial_travel_evolves_real_geometry_without_resetting_twist_or_clock() {
    let source = source();
    let (_, _, time, twist) = scales(&source);
    let c = geometry_controls(&source);
    let mut f = axial();
    spin(&mut f, 0, Vec3::new(0.0, 0.0, 0.05 * twist / time));
    let mut pair = PreparedRigidJkrTorsionalPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&source),
        c,
        &CancelGate::new(),
    )
    .unwrap();
    let initial_twist = 0.05 * twist;
    let mut state = pair.initialize(initial_twist, &CancelGate::new()).unwrap();
    let initial_energy = total_energy(&f, &source, state);
    let initial_radius = source
        .evaluate_torsion_displacement_controlled(state.coordinates()[0], state.coordinates()[1])
        .unwrap()
        .contact_radius_m;
    let initial_bodies = f.bodies;
    let dt = 0.002 * time;
    for step in 0..96 {
        let candidate = pair.advance(state, dt, &CancelGate::new()).unwrap();
        if step == 0 {
            // The accepted normal travel must actually reach the native poses
            // before another phase can use the new contact geometry.
            let held = PreparedRigidJkrTorsionalPair::new(
                &candidate.bodies_after,
                &f.properties,
                f.contact(&source),
                c,
                &CancelGate::new(),
            )
            .unwrap();
            assert!(
                held.advance(candidate.state_after, dt, &CancelGate::new())
                    .is_err()
            );
        }
        let next = pair
            .continue_axial_geometry(candidate.state_after, &CancelGate::new())
            .unwrap();
        assert_eq!(
            next.state().coordinates(),
            candidate.state_after.coordinates()
        );
        assert_eq!(next.state().time_s(), candidate.state_after.time_s());
        assert_eq!(next.state().steps(), candidate.state_after.steps());
        for i in 0..2 {
            assert_eq!(
                next.bodies()[i].linear_momentum_world(),
                candidate.bodies_after[i].linear_momentum_world()
            );
            assert_eq!(
                next.bodies()[i].angular_momentum_body(),
                candidate.bodies_after[i].angular_momentum_body()
            );
        }
        close(
            next.report().audit.kinetic_energy_residual_j,
            0.0,
            c.energy_tolerance_j,
        );
        f.bodies = *next.bodies();
        let (new_pair, new_state, _) = next.into_parts();
        pair = new_pair;
        state = new_state;
        close(
            total_energy(&f, &source, state),
            initial_energy,
            20.0 * c.energy_tolerance_j,
        );
    }
    assert_eq!(state.steps(), 96);
    close(state.time_s(), 96.0 * dt, 64.0 * f64::EPSILON * 96.0 * dt);
    assert_ne!(f.bodies[0].pose(), initial_bodies[0].pose());
    let q = state.coordinates();
    let final_contact = source
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap();
    assert!((final_contact.contact_radius_m - initial_radius).abs() > 1e-9 * initial_radius);
    let markers: [Vec3; 2] =
        core::array::from_fn(|i| f.bodies[i].pose().point_world_from_body(f.arms[i]).unwrap());
    close(
        f.contact(&source).reference_approach_m - markers[1].sub(markers[0]).z,
        q[0],
        c.geometry_tolerance_m,
    );
    let angle: [f64; 2] = f.bodies.map(|body| {
        let [w, _, _, z] = body.pose().orientation().components();
        2.0 * z.atan2(w)
    });
    close(
        initial_twist + angle[0] - angle[1],
        q[1],
        20.0 * c.twist_tolerance_rad,
    );
}

#[test]
fn cancelled_geometry_handoff_preserves_the_accepted_native_and_source_endpoints() {
    let source = source();
    let (_, _, time, twist) = scales(&source);
    let f = axial();
    let pair = PreparedRigidJkrTorsionalPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&source),
        geometry_controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    let initial = pair.initialize(0.05 * twist, &CancelGate::new()).unwrap();
    let candidate = pair
        .advance(initial, 0.002 * time, &CancelGate::new())
        .unwrap();
    let before_bodies = candidate.bodies_after;
    let before_state = candidate.state_after;
    let gate = CancelGate::new();
    gate.request();
    assert!(pair.continue_axial_geometry(before_state, &gate).is_err());
    assert_eq!(candidate.bodies_after, before_bodies);
    assert_eq!(candidate.state_after, before_state);
    let first = pair
        .continue_axial_geometry(before_state, &CancelGate::new())
        .unwrap();
    let retry = pair
        .continue_axial_geometry(before_state, &CancelGate::new())
        .unwrap();
    assert_eq!(first.bodies(), retry.bodies());
    assert_eq!(first.state(), retry.state());
    assert_eq!(first.report(), retry.report());
}

#[test]
fn loose_audit_tolerances_cannot_authorize_noncommuting_finite_rotation() {
    let source = source();
    let (_, _, time, twist) = scales(&source);
    let mut f = axial();
    spin(&mut f, 0, Vec3::new(1.0, 0.0, 0.0));
    spin(&mut f, 1, Vec3::new(1.0, 0.0, 0.0));
    let c = RigidJkrTorsionalControls {
        velocity_tolerance_m_s: 1e6,
        angular_velocity_tolerance_rad_s: 1e6,
        angular_momentum_tolerance_nms: 1e6,
        ..controls(&source)
    };
    let pair = PreparedRigidJkrTorsionalPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&source),
        c,
        &CancelGate::new(),
    )
    .unwrap();
    let state = pair.initialize(0.05 * twist, &CancelGate::new()).unwrap();
    let candidate = pair
        .advance(state, 0.002 * time, &CancelGate::new())
        .unwrap();
    assert!(matches!(
        pair.continue_axial_geometry(candidate.state_after, &CancelGate::new()),
        Err(RigidJkrTorsionalGeometryError::NoncommutingMotion),
    ));
    assert_eq!(candidate.bodies_after[0].pose(), f.bodies[0].pose());
    assert_eq!(state.steps(), 0);
}
