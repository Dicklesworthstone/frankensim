//! Moving finite film carriers share the sampled sphere/fluid transition.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, DryBodyCouplingConfig3, E3,
    LubricationPlane3, MovingFiniteLubricationPlane3, MovingLubricationSurface3,
    MovingPartialSaturationBody3, PartialSaturationError3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
const RADIUS: f64 = 0.7;
const LOAD: [f64; 6] = [2e-5, -1e-5, 3e-5, 1e-6, -2e-6, 1e-6];
const CARRIER_LOAD: [f64; 6] = [-3e-5, 2e-5, -1e-5, -2e-6, 1e-6, 3e-6];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn array(v: Vec3) -> [f64; 3] { [v.x, v.y, v.z] }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a.into_iter().zip(b).map(|(x,y)| x*y).sum() }
fn grid() -> BoundaryGrid3 { BoundaryGrid3::new(8, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic()) }
fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center, mass, velocity,
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap(),
        inertia: InertiaTensor3::new(Mat3 { m: [0.12, 0.009, -0.004, 0.009, 0.17, 0.006, -0.004, 0.006, 0.21].map(|x| mass*x) }).unwrap(),
        angular_momentum_body: Vec3::new(2e-6, -3e-6, 4e-6),
    }
}
fn sphere() -> MovingPartialSaturationBody3 { body([3.5, 4.5, 4.5], 3.0, [-8e-5, 4e-5, -2e-5]) }
fn carrier() -> MovingFiniteLubricationPlane3 {
    let body = body([2.6, 4.1, 4.8], 5.0, [2e-5, -1e-5, 3e-5]);
    MovingFiniteLubricationPlane3::new(
        body,
        array(body.orientation.inverse().rotate(Vec3::new(0.19, 0.4, -0.3)).unwrap()),
        array(body.orientation.inverse().rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap()),
        0.03, 0.001,
    ).unwrap().with_tangential(0.013, 0.001).unwrap()
}
fn prescribed() -> LubricationPlane3 {
    LubricationPlane3::new([4.21, 4.3, 4.6], [-1.0, 0.0, 0.0], [1e-5, -2e-5, 1e-5], 0.03, 0.001)
        .unwrap().with_tangential(0.013, 0.001).unwrap().with_axial_spin(3e-5).unwrap()
}
fn surfaces() -> [MovingLubricationSurface3; 2] {
    [MovingLubricationSurface3::Finite(carrier()), MovingLubricationSurface3::Prescribed(prescribed())]
}
fn finite(surface: MovingLubricationSurface3) -> MovingFiniteLubricationPlane3 {
    let MovingLubricationSurface3::Finite(p) = surface else { panic!("finite surface") }; p
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.11, 0.1, 0.001, 2).unwrap().with_resistance_mode_limit(10).unwrap()
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..512).map(|i| g.populations(i % 8, (i / 8) % 8, i / 64)).collect()
}
fn momentum(b: MovingPartialSaturationBody3, origin: [f64; 3]) -> [f64; 6] {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| b.center[i] - origin[i]), p);
    core::array::from_fn(|i| if i < 3 { p[i] } else { spin[i - 3] + orbital[i - 3] })
}
fn applied(b: MovingPartialSaturationBody3, j: [f64; 6], origin: [f64; 3]) -> [f64; 6] {
    let orbital = cross(core::array::from_fn(|i| b.center[i] + 0.5*b.velocity[i] - origin[i]), core::array::from_fn(|i| j[i]));
    core::array::from_fn(|i| if i < 3 { j[i] } else { j[i] + orbital[i - 3] })
}
fn energy(b: MovingPartialSaturationBody3) -> f64 {
    0.5*b.mass*dot(b.velocity,b.velocity) + 0.5*dot(array(b.orientation.rotate(b.angular_momentum_body).unwrap()), b.angular_velocity_world().unwrap())
}

#[test]
fn moving_finite_carrier_and_spinning_reservoir_close_the_shared_fluid_tick() {
    let mut g = grid();
    let mut fields = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut spheres = [sphere()];
    let mut surfaces = surfaces();
    for _ in 0..3 {
        let old_sphere = spheres[0];
        let old_carrier = finite(surfaces[0]);
        let before = snapshot(&g);
        let report = g.step_with_lubricated_spheres_and_surfaces(
            &mut fields, [RADIUS], &mut spheres, [LOAD], model(), &mut surfaces,
            &[CARRIER_LOAD, [0.0; 6]], ORIGIN, &CancelGate::new(),
        ).unwrap();
        let carrier = finite(surfaces[0]);
        let account = report.surfaces[0].finite.unwrap();
        assert!(report.surfaces[1].finite.is_none());
        assert_eq!(report.motion.coupling.lubrication.active_plane_pairs, 2);
        assert_eq!(report.motion.coupling.lubrication.active_resistance_modes, 10);
        assert!(report.motion.coupling.bodies[0].hydrodynamic.coupled_cells > 0);
        assert_ne!(carrier.normal(), old_carrier.normal());
        assert_ne!(carrier.point(), old_carrier.point());
        assert_eq!(carrier.body().inertia, old_carrier.body().inertia);
        let (q_mid, _) = old_carrier.body().inertia.momentum_step(
            old_carrier.body().orientation, old_carrier.body().angular_momentum_body,
            Vec3::new(0.0, 0.0, 0.0), 0.5,
        ).unwrap();
        assert_eq!(account.midpoint_orientation, q_mid);
        let expected = q_mid.matrix().compose(carrier.body().inertia.matrix()).compose(q_mid.matrix().transpose());
        let solved = report.surfaces[0].film.finite_body.unwrap();
        for i in 0..3 { for j in 0..3 {
            close(solved.body_after.inertia_world[i][j], expected.m[3*i+j], 2e-16);
        }}
        let after = snapshot(&g);
        let mut change = [0.0; 6];
        for cell in 0..512 {
            let dp: [f64; 3] = core::array::from_fn(|i| (0..Q3).map(|q| {
                let e = [E3[q].0, E3[q].1, E3[q].2];
                f64::from(e[i]) * (after[cell][q] - before[cell][q])
            }).sum());
            let dl = cross([ (cell%8) as f64 + 0.5-ORIGIN[0], ((cell/8)%8) as f64+0.5-ORIGIN[1], (cell/64) as f64+0.5-ORIGIN[2] ], dp);
            for i in 0..3 { change[i] += dp[i]; change[i+3] += dl[i]; }
        }
        for (before, after, load) in [(old_sphere, spheres[0], LOAD), (old_carrier.body(), carrier.body(), CARRIER_LOAD)] {
            let p0 = momentum(before, ORIGIN); let p1 = momentum(after, ORIGIN);
            let external = applied(before, load, ORIGIN);
            for i in 0..6 { change[i] += p1[i] - p0[i] - external[i]; }
        }
        for i in 0..3 {
            close(change[i] + report.surfaces[1].film.impulse[i], 0.0, 2e-12);
            close(change[i+3] + report.surfaces[1].film.angular_impulse[i] - report.motion.coupling.periodic.angular_momentum_change[i], 0.0, 8e-12);
        }
        let p0 = momentum(old_carrier.body(), ORIGIN); let p1 = momentum(carrier.body(), ORIGIN);
        for i in 0..6 {
            close(account.momentum_change[i], p1[i]-p0[i], 3e-18);
            close(account.momentum_residual[i], 0.0, 3e-17);
        }
        close(account.kinetic_energy_change, energy(carrier.body()) - energy(old_carrier.body()), 1e-22);
        close(account.work_residual, 0.0, 5e-21);
        assert!(report.motion.coupling.lubrication.dissipated_energy >= 0.0);
        close(report.motion.coupling.lubrication.work_residual, 0.0, 5e-21);
    }
}

#[test]
fn full_tensor_surface_dynamics_are_independent_of_the_reporting_origin() {
    let origins = [ORIGIN, [8.0, -5.0, 2.0]];
    let mut results = Vec::new();
    for origin in origins {
        let mut g = grid(); let mut fields = [BodyFrameSampler3::new(&g, 4).unwrap()];
        let mut bodies = [sphere()]; let mut walls = surfaces();
        let report = g.step_with_lubricated_spheres_and_surfaces(&mut fields, [RADIUS], &mut bodies,
            [LOAD], model(), &mut walls, &[CARRIER_LOAD, [0.0; 6]], origin, &CancelGate::new()).unwrap();
        results.push((bodies, walls, snapshot(&g), report));
    }
    assert_eq!(results[0].0, results[1].0);
    assert_eq!(results[0].1, results[1].1);
    assert_eq!(results[0].2, results[1].2);
    let shift = core::array::from_fn(|i| origins[1][i] - origins[0][i]);
    for k in 0..2 {
        let a = results[0].3.surfaces[k]; let b = results[1].3.surfaces[k];
        let delta = cross(shift, a.film.impulse);
        for i in 0..3 { close(b.film.angular_impulse[i], a.film.angular_impulse[i]-delta[i], 2e-18); }
        close(a.film.work, b.film.work, 2e-22);
    }
    let a = results[0].3.surfaces[0].finite.unwrap(); let b = results[1].3.surfaces[0].finite.unwrap();
    for (x,y) in [(a.momentum_change,b.momentum_change),(a.external_impulse,b.external_impulse)] {
        let delta = cross(shift, core::array::from_fn(|i| x[i]));
        for i in 0..3 { close(y[i+3], x[i+3]-delta[i], 2e-18); }
    }
}

#[test]
fn prescribed_and_empty_surface_paths_preserve_existing_moving_and_dry_arithmetic() {
    for dry in [false, true] { for planes in [Vec::new(), vec![prescribed()]] {
        let mut a = grid(); let mut b = grid();
        let mut fa = [BodyFrameSampler3::new(&a, 4).unwrap()]; let mut fb = [BodyFrameSampler3::new(&b, 4).unwrap()];
        let mut ba = [sphere()]; let mut bb = ba;
        let mut walls: Vec<_> = planes.iter().copied().map(MovingLubricationSurface3::Prescribed).collect();
        let loads = vec![[0.0; 6]; planes.len()]; let gate = CancelGate::new();
        if dry {
            let old = a.step_with_dry_lubricated_spheres_and_planes(&mut fa, [RADIUS], &mut ba, [LOAD], model(), &planes, ORIGIN, DryBodyCouplingConfig3::default(), &gate).unwrap();
            let new = b.step_with_dry_lubricated_spheres_and_surfaces(&mut fb, [RADIUS], &mut bb, [LOAD], model(), &mut walls, &loads, ORIGIN, DryBodyCouplingConfig3::default(), &gate).unwrap();
            assert_eq!(old.motion, new.motion);
            assert_eq!(old.planes, new.surfaces.iter().map(|s| s.film).collect::<Vec<_>>());
        } else {
            let old = a.step_with_lubricated_spheres_and_planes(&mut fa, [RADIUS], &mut ba, [LOAD], model(), &planes, ORIGIN, &gate).unwrap();
            let new = b.step_with_lubricated_spheres_and_surfaces(&mut fb, [RADIUS], &mut bb, [LOAD], model(), &mut walls, &loads, ORIGIN, &gate).unwrap();
            assert_eq!(old.motion, new.motion);
            assert_eq!(old.planes, new.surfaces.iter().map(|s| s.film).collect::<Vec<_>>());
        }
        assert_eq!(ba,bb); assert_eq!(snapshot(&a),snapshot(&b));
        for (&old,&new) in planes.iter().zip(&walls) {
            let MovingLubricationSurface3::Prescribed(new) = new else { panic!("prescribed") };
            assert_eq!(new.point(), core::array::from_fn(|i| old.point()[i]+old.velocity()[i]));
        }
    }}
}

#[test]
fn dry_candidate_refusals_preserve_all_owners_and_retry_matches_fresh() {
    let mut g = grid(); let original = snapshot(&g);
    let mut fields = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()]; let initial = bodies; let mut walls = surfaces(); let original_walls = walls;
    let gate = CancelGate::new();
    let run = |g: &mut BoundaryGrid3, fields: &mut [BodyFrameSampler3; 1], bodies: &mut [MovingPartialSaturationBody3; 1], walls: &mut [MovingLubricationSurface3], loads: &[[f64; 6]], config, gate: &CancelGate| {
        g.step_with_dry_lubricated_spheres_and_surfaces(fields, [RADIUS], bodies, [LOAD], model(), walls, loads, ORIGIN, config, gate)
    };
    for loads in [vec![], vec![CARRIER_LOAD], vec![CARRIER_LOAD, LOAD], vec![CARRIER_LOAD, [0.0; 6], [0.0; 6]]] {
        assert!(run(&mut g,&mut fields,&mut bodies,&mut walls,&loads,DryBodyCouplingConfig3::default(),&gate).is_err());
        assert_eq!(bodies,initial); assert_eq!(walls,original_walls); assert_eq!(snapshot(&g),original);
    }
    let loads = [CARRIER_LOAD, [0.0; 6]];
    for count in [2, 65] {
        let mut duplicate = vec![MovingLubricationSurface3::Prescribed(prescribed()); count];
        let unchanged = duplicate.clone();
        assert!(run(&mut g, &mut fields, &mut bodies, &mut duplicate,
            &vec![[0.0; 6]; count], DryBodyCouplingConfig3::default(), &gate).is_err());
        assert_eq!(duplicate, unchanged);
        assert_eq!(bodies, initial);
        assert_eq!(snapshot(&g), original);
    }
    let limited = DryBodyCouplingConfig3 { maximum_evaluations: 1, ..DryBodyCouplingConfig3::default() };
    assert!(run(&mut g,&mut fields,&mut bodies,&mut walls,&loads,limited,&gate).is_err());
    assert_eq!(bodies,initial); assert_eq!(walls,original_walls); assert_eq!(snapshot(&g),original);
    let cancelled = CancelGate::new(); cancelled.request();
    assert_eq!(run(&mut g,&mut fields,&mut bodies,&mut walls,&loads,DryBodyCouplingConfig3::default(),&cancelled),Err(PartialSaturationError3::Cancelled));
    assert_eq!(bodies,initial); assert_eq!(walls,original_walls); assert_eq!(snapshot(&g),original);
    let accepted = run(&mut g,&mut fields,&mut bodies,&mut walls,&loads,DryBodyCouplingConfig3::default(),&gate).unwrap();
    assert!(accepted.motion.iterations > 0);
    let carrier_report = accepted.surfaces[0].finite.unwrap();
    let carrier_before = momentum(finite(original_walls[0]).body(), ORIGIN);
    let carrier_after = momentum(finite(walls[0]).body(), ORIGIN);
    for i in 0..6 {
        close(carrier_report.momentum_change[i], carrier_after[i] - carrier_before[i], 3e-18);
        close(carrier_report.momentum_residual[i], 0.0, 3e-17);
        close(accepted.motion.combined_momentum_residual[i], 0.0, 3e-10);
    }
    close(carrier_report.work_residual, 0.0, 5e-21);
    let mut fresh = grid(); let mut fields2 = [BodyFrameSampler3::new(&fresh,4).unwrap()];
    let mut bodies2 = initial; let mut walls2 = original_walls;
    let repeated = run(&mut fresh,&mut fields2,&mut bodies2,&mut walls2,&loads,DryBodyCouplingConfig3::default(),&gate).unwrap();
    assert_eq!(accepted,repeated); assert_eq!(bodies,bodies2); assert_eq!(walls,walls2); assert_eq!(snapshot(&g),snapshot(&fresh));
}

#[test]
fn finite_carrier_crossing_in_the_second_half_refuses_before_fluid_publication() {
    let mut g = grid();
    let mut fields = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let initial = bodies;
    let p = carrier();
    let local_point = p.body().orientation.inverse().rotate(Vec3::new(
        p.point()[0] - p.body().center[0], p.point()[1] - p.body().center[1], p.point()[2] - p.body().center[2],
    )).unwrap();
    let local_normal = p.body().orientation.inverse().rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
    // This inactive film permits a bounded external load to move the finite
    // plane through the sphere during the SECOND half-drift. Both the initial
    // and collision-midpoint gaps are valid, so this exercises late admission.
    let p = MovingFiniteLubricationPlane3::new(p.body(), array(local_point), array(local_normal), 0.001, 0.0001).unwrap();
    let mut walls = [MovingLubricationSurface3::Finite(p)];
    let original_walls = walls;
    let original_fluid = snapshot(&g);
    let error = g.step_with_lubricated_spheres_and_surfaces(
        &mut fields, [RADIUS], &mut bodies, [LOAD], model(), &mut walls,
        &[[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]], ORIGIN, &CancelGate::new(),
    ).unwrap_err();
    assert!(matches!(error, PartialSaturationError3::Invalid { reason } if reason.contains("half-drift cannot exclude contact")), "{error:?}");
    assert_eq!(bodies, initial);
    assert_eq!(walls, original_walls);
    assert_eq!(snapshot(&g), original_fluid);
    let accepted = g.step_with_lubricated_spheres_and_surfaces(
        &mut fields, [RADIUS], &mut bodies, [LOAD], model(), &mut walls,
        &[CARRIER_LOAD], ORIGIN, &CancelGate::new(),
    ).unwrap();
    let mut fresh = grid();
    let mut fresh_fields = [BodyFrameSampler3::new(&fresh, 4).unwrap()];
    let mut fresh_bodies = initial;
    let mut fresh_walls = original_walls;
    let repeated = fresh.step_with_lubricated_spheres_and_surfaces(
        &mut fresh_fields, [RADIUS], &mut fresh_bodies, [LOAD], model(), &mut fresh_walls,
        &[CARRIER_LOAD], ORIGIN, &CancelGate::new(),
    ).unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(bodies, fresh_bodies);
    assert_eq!(walls, fresh_walls);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}

#[test]
fn clear_endpoints_cannot_hide_a_rotating_collision_or_crossing_duplicate() {
    let mut g = grid();
    let original_fluid = snapshot(&g);
    let mut fields = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let initial = bodies;
    let body = MovingPartialSaturationBody3 {
        center: [2.79, 4.5, 4.5],
        orientation: So3::identity(),
        mass: 5.0,
        inertia: InertiaTensor3::new(Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        }).unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 4.0 * core::f64::consts::PI),
    };
    let (half_orientation, _) = body.inertia.momentum_step(body.orientation,
        body.angular_momentum_body, Vec3::new(0.0, 0.0, 0.0), 0.5).unwrap();
    let half_normal = array(half_orientation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap());
    let half_arm = core::array::from_fn(|i| initial[0].center[i] + 0.5 * initial[0].velocity[i] - body.center[i]);
    assert!(dot(half_normal, half_arm) - RADIUS > 0.009);
    let surface = MovingFiniteLubricationPlane3::new(body, [0.0; 3], [1.0, 0.0, 0.0], 0.03, 0.001).unwrap();
    let mut walls = [MovingLubricationSurface3::Finite(surface)];
    let original_walls = walls;
    let error = g.step_with_lubricated_spheres_and_surfaces(&mut fields, [RADIUS],
        &mut bodies, [LOAD], model(), &mut walls, &[[0.0; 6]], ORIGIN, &CancelGate::new()).unwrap_err();
    assert!(matches!(error, PartialSaturationError3::Invalid { reason } if reason.contains("half-drift cannot exclude contact")), "{error:?}");
    assert_eq!(bodies, initial);
    assert_eq!(walls, original_walls);
    assert_eq!(snapshot(&g), original_fluid);

    let plane = |x, v| MovingLubricationSurface3::Prescribed(
        LubricationPlane3::new([x, 0.0, 0.0], [1.0, 0.0, 0.0], [v, 0.0, 0.0], 0.03, 0.001).unwrap()
    );
    let mut walls = [plane(0.0, 0.04), plane(0.01, -0.04)];
    let original_walls = walls;
    let error = g.step_with_lubricated_spheres_and_surfaces(&mut fields, [RADIUS],
        &mut bodies, [LOAD], model(), &mut walls, &[[0.0; 6]; 2], ORIGIN, &CancelGate::new()).unwrap_err();
    assert!(matches!(error, PartialSaturationError3::Invalid { reason } if reason.contains("cannot be proven distinct")), "{error:?}");
    assert_eq!(bodies, initial);
    assert_eq!(walls, original_walls);
    assert_eq!(snapshot(&g), original_fluid);
}
