//! G1 simultaneous analytic recoil; G3 rigid-frame work; G4 value rollback.

use super::super::finite_planes::FinitePlaneInput;
use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        inertia_world: [
            [0.12, 0.009, -0.004],
            [0.009, 0.17, 0.006],
            [-0.004, 0.006, 0.21],
        ]
        .map(|row| row.map(|x| mass * x)),
        velocity,
        angular_velocity: [0.0; 3],
    }
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}

fn energy(body: PartialSaturationBody3) -> f64 {
    0.5 * body.mass * dot(body.velocity, body.velocity)
        + 0.5
            * dot(
                body.angular_velocity,
                body.inertia_world.map(|r| dot(r, body.angular_velocity)),
            )
}

fn momentum(body: PartialSaturationBody3, origin: [f64; 3]) -> [f64; 6] {
    let p = body.velocity.map(|v| body.mass * v);
    let orbital = cross(difference(body.center, origin), p);
    let spin = body.inertia_world.map(|r| dot(r, body.angular_velocity));
    core::array::from_fn(|i| {
        if i < 3 {
            p[i]
        } else {
            orbital[i - 3] + spin[i - 3]
        }
    })
}

fn world_load(body: PartialSaturationBody3, load: [f64; 6], origin: [f64; 3]) -> [f64; 6] {
    let orbital = cross(
        difference(body.center, origin),
        core::array::from_fn(|i| load[i]),
    );
    core::array::from_fn(|i| {
        if i < 3 {
            load[i]
        } else {
            load[i] + orbital[i - 3]
        }
    })
}

fn work(before: PartialSaturationBody3, after: PartialSaturationBody3, load: [f64; 6]) -> f64 {
    (0..3)
        .map(|i| {
            0.5 * load[i] * (before.velocity[i] + after.velocity[i])
                + 0.5 * load[i + 3] * (before.angular_velocity[i] + after.angular_velocity[i])
        })
        .sum()
}

fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.11, 0.2, 0.001, 2)
        .unwrap()
        .with_resistance_mode_limit(10)
        .unwrap()
}

fn evaluate(
    initial: PartialSaturationBody3,
    response: Response,
    surfaces: &mut [LubricationSurface3],
    loads: &[[f64; 6]],
    origin: [f64; 3],
    gate: &CancelGate,
) -> Result<(
    PartialSaturationBody3,
    SphereLubricationExchange3<1>,
    Vec<PlaneLubricationExchange3>,
)> {
    let native = prepare(surfaces, loads, gate)?;
    let mut proposed = [initial];
    let (film, reports) = solve(
        model(),
        [0.5],
        &native,
        &[initial],
        &[response],
        &mut proposed,
        origin,
        gate,
    )?;
    accept(surfaces, &reports);
    Ok((proposed[0], film, reports))
}

fn finite(surface: LubricationSurface3) -> PartialSaturationBody3 {
    let LubricationSurface3::Finite(plane) = surface else {
        panic!("finite test carrier")
    };
    plane.body()
}

#[test]
fn mixed_normal_films_match_independent_two_inertia_equations_with_fluid_response() {
    let initial = body([3.0, 4.0, 4.0], 2.0, [-0.01, 0.0, 0.0]);
    let wall = body([2.25, 4.0, 4.0], 3.0, [0.015, 0.0, 0.0]);
    let prescribed = LubricationPlane3::new(
        [3.75, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-0.004, 0.0, 0.0],
        0.75,
        0.02,
    )
    .unwrap();
    let mut surfaces = [
        LubricationSurface3::Finite(
            FiniteLubricationPlane3::new([2.25, 0.0, 0.0], [1.0, 0.0, 0.0], wall, 0.75, 0.02)
                .unwrap(),
        ),
        LubricationSurface3::Prescribed(prescribed),
    ];
    let load = [0.004, 0.0, 0.0, 0.0, 0.0, 0.0];
    let wall_load = [-0.003, 0.0, 0.0, 0.0, 0.0, 0.0];
    let fluid_weight = 0.8;
    let fluid_velocity = 0.002;
    let mut response = Response::new(load);
    response.add(
        fluid_weight,
        [0.0; 3],
        [fluid_velocity - initial.velocity[0], 0.0, 0.0],
    );
    let zeta = 6.0 * core::f64::consts::PI * 0.11 * 0.25 * (4.0 - 1.0 / 0.75);
    // Eliminate the two film tractions from physical body and carrier balances.
    // This 2x2 velocity system is independent of the production modal Schur solve.
    let a = initial.mass + 0.5 * fluid_weight + zeta;
    let b = -0.5 * zeta;
    let d = wall.mass + 0.5 * zeta;
    let rhs_s = load[0] + fluid_weight * (fluid_velocity - initial.velocity[0])
        - zeta * (initial.velocity[0] - wall.velocity[0])
        - zeta * (initial.velocity[0] - prescribed.velocity()[0]);
    let rhs_p = wall_load[0] + zeta * (initial.velocity[0] - wall.velocity[0]);
    let determinant = a * d - b * b;
    let dv = (d * rhs_s - b * rhs_p) / determinant;
    let du = (a * rhs_p - b * rhs_s) / determinant;
    let (after, film, reports) = evaluate(
        initial,
        response,
        &mut surfaces,
        &[wall_load, [0.0; 6]],
        [-2.0, 0.7, 1.2],
        &CancelGate::new(),
    )
    .unwrap();
    let wall_after = finite(surfaces[0]);
    close(after.velocity[0], initial.velocity[0] + dv, 4e-17);
    close(wall_after.velocity[0], wall.velocity[0] + du, 4e-17);
    assert_eq!(surfaces[1], LubricationSurface3::Prescribed(prescribed));
    assert!(reports[0].finite_body.is_some());
    assert!(reports[1].finite_body.is_none());
    assert_eq!(reports[0].finite_body.unwrap().body_after, wall_after);
    assert_eq!(film.active_plane_pairs, 2);
    let vmid = 0.5 * (initial.velocity[0] + after.velocity[0]);
    let umid = 0.5 * (wall.velocity[0] + wall_after.velocity[0]);
    let left = -zeta * (vmid - umid);
    let right = -zeta * (vmid - prescribed.velocity()[0]);
    close(reports[0].impulse[0], -left, 4e-17);
    close(reports[1].impulse[0], -right, 4e-17);
    close(
        film.dissipated_energy,
        zeta * ((vmid - umid).powi(2) + (vmid - prescribed.velocity()[0]).powi(2)),
        3e-18,
    );
    let hydro = fluid_weight * (fluid_velocity - vmid);
    close(
        initial.mass * dv + wall.mass * du - load[0] - wall_load[0] - hydro + reports[1].impulse[0],
        0.0,
        5e-17,
    );
    close(
        energy(after) - energy(initial) + energy(wall_after)
            - energy(wall)
            - work(initial, after, load)
            - work(wall, wall_after, wall_load)
            - hydro * vmid
            + reports[1].work
            + film.dissipated_energy,
        0.0,
        4e-18,
    );
    close(film.work_residual, 0.0, 3e-18);
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [v[2], v[0], v[1]]
}
fn rotate_body(b: PartialSaturationBody3) -> PartialSaturationBody3 {
    let order = [2, 0, 1];
    PartialSaturationBody3 {
        center: rotate(b.center),
        velocity: rotate(b.velocity),
        angular_velocity: rotate(b.angular_velocity),
        inertia_world: core::array::from_fn(|i| {
            core::array::from_fn(|j| b.inertia_world[order[i]][order[j]])
        }),
        ..b
    }
}
fn rotate_load(j: [f64; 6]) -> [f64; 6] {
    [j[2], j[0], j[1], j[5], j[3], j[4]]
}

#[test]
fn mixed_spinning_surface_and_full_finite_recoil_keep_world_momentum_and_work_objective() {
    let mut initial = body([3.0, 4.0, 4.0], 2.0, [-0.012, 0.007, -0.003]);
    initial.angular_velocity = [0.003, -0.004, 0.002];
    let mut wall = body([2.2, 4.3, 3.7], 4.0, [-0.001, 0.003, 0.002]);
    wall.angular_velocity = [-0.003, 0.004, 0.005];
    let load = [0.0004, -0.0002, 0.0001, 0.00003, -0.00001, 0.00002];
    let wall_load = [0.0002, -0.0001, 0.0003, 0.00001, 0.00002, -0.00003];
    let arm = [0.2, -0.3, 0.1];
    let fluid_velocity = [0.001, -0.002, 0.003];
    let origin = [-2.0, 0.7, 1.2];
    let mut reference = None;
    for rotated in [false, true] {
        let transform = |v| if rotated { rotate(v) } else { v };
        let sphere = if rotated {
            rotate_body(initial)
        } else {
            initial
        };
        let carrier = if rotated { rotate_body(wall) } else { wall };
        let j = if rotated { rotate_load(load) } else { load };
        let jp = if rotated {
            rotate_load(wall_load)
        } else {
            wall_load
        };
        let r = transform(arm);
        let u = transform(fluid_velocity);
        let o = transform(origin);
        let prescribed = LubricationPlane3::new(
            transform([3.505, 3.8, 4.2]),
            transform([-1.0, 0.0, 0.0]),
            transform([0.0, 0.003, -0.002]),
            0.012,
            0.001,
        )
        .unwrap()
        .with_tangential(0.01, 0.001)
        .unwrap()
        .with_axial_spin(0.003)
        .unwrap();
        let mut surfaces = [
            LubricationSurface3::Finite(
                FiniteLubricationPlane3::new(
                    transform([2.495, 0.0, 0.0]),
                    transform([1.0, 0.0, 0.0]),
                    carrier,
                    0.012,
                    0.001,
                )
                .unwrap()
                .with_tangential(0.01, 0.001)
                .unwrap(),
            ),
            LubricationSurface3::Prescribed(prescribed),
        ];
        let mut response = Response::new(j);
        let slip = difference(
            u,
            point_velocity(sphere.velocity, sphere.angular_velocity, r),
        );
        response.add(0.4, r, slip);
        let (after, film, reports) = evaluate(
            sphere,
            response,
            &mut surfaces,
            &[jp, [0.0; 6]],
            o,
            &CancelGate::new(),
        )
        .unwrap();
        let carrier_after = finite(surfaces[0]);
        assert_eq!(film.active_resistance_modes, 10);
        assert!(film.dissipated_energy > 0.0);
        assert!(
            reports[0]
                .finite_body
                .unwrap()
                .spin_impulse
                .iter()
                .any(|v| v.abs() > 1e-5)
        );
        let vmid = core::array::from_fn(|i| 0.5 * (sphere.velocity[i] + after.velocity[i]));
        let wmid = core::array::from_fn(|i| {
            0.5 * (sphere.angular_velocity[i] + after.angular_velocity[i])
        });
        let at_point = point_velocity(vmid, wmid, r);
        let hydro = difference(u, at_point).map(|v| 0.4 * v);
        let hydro_spin = cross(r, hydro);
        let hydro_load = core::array::from_fn(|i| if i < 3 { hydro[i] } else { hydro_spin[i - 3] });
        let before_s = momentum(sphere, o);
        let after_s = momentum(after, o);
        let before_p = momentum(carrier, o);
        let after_p = momentum(carrier_after, o);
        let ext_s = world_load(sphere, j, o);
        let ext_p = world_load(carrier, jp, o);
        let h = world_load(sphere, hydro_load, o);
        for i in 0..6 {
            let prescribed_reaction = if i < 3 {
                reports[1].impulse[i]
            } else {
                reports[1].angular_impulse[i - 3]
            };
            close(
                after_s[i] - before_s[i] + after_p[i] - before_p[i] - ext_s[i] - ext_p[i] - h[i]
                    + prescribed_reaction,
                0.0,
                2e-15,
            );
        }
        close(
            energy(after) - energy(sphere) + energy(carrier_after)
                - energy(carrier)
                - work(sphere, after, j)
                - work(carrier, carrier_after, jp)
                - dot(hydro, at_point)
                + reports[1].work
                + film.dissipated_energy,
            0.0,
            2e-17,
        );
        close(film.work_residual, 0.0, 2e-17);
        close(reports[0].finite_body.unwrap().work_residual, 0.0, 2e-17);
        assert_eq!(surfaces[1], LubricationSurface3::Prescribed(prescribed));
        if let Some((old_after, old_carrier, old_loss, old_work)) = reference {
            let expected_sphere = rotate_body(old_after);
            let expected_carrier = rotate_body(old_carrier);
            for i in 0..3 {
                close(after.velocity[i], expected_sphere.velocity[i], 2e-16);
                close(
                    after.angular_velocity[i],
                    expected_sphere.angular_velocity[i],
                    2e-16,
                );
                close(
                    carrier_after.velocity[i],
                    expected_carrier.velocity[i],
                    2e-16,
                );
                close(
                    carrier_after.angular_velocity[i],
                    expected_carrier.angular_velocity[i],
                    2e-16,
                );
            }
            close(film.dissipated_energy, old_loss, 2e-17);
            close(reports[1].work, old_work, 2e-17);
        } else {
            reference = Some((
                after,
                carrier_after,
                film.dissipated_energy,
                reports[1].work,
            ));
        }
    }
}

#[test]
fn homogeneous_and_empty_surface_sets_retain_the_original_response_bits() {
    let initial = body([3.0, 4.0, 4.0], 2.0, [-0.01, 0.002, -0.003]);
    let wall = body([2.25, 4.0, 4.0], 3.0, [0.015, -0.001, 0.002]);
    let geometry =
        LubricationPlane3::new([2.25, 0.0, 0.0], [1.0, 0.0, 0.0], wall.velocity, 0.75, 0.02)
            .unwrap();
    let carrier =
        FiniteLubricationPlane3::new(geometry.point(), geometry.normal(), wall, 0.75, 0.02)
            .unwrap();
    let load = [0.003, 0.001, -0.002, 0.0002, -0.0001, 0.0003];
    let mut finite_geometry = geometry;
    finite_geometry.recoil = Some(FinitePlaneInput {
        body: wall,
        external: load,
        carrier_anchor: None,
    });
    let response = Response::new([0.0; 6]);
    for (mut surfaces, external, legacy) in [
        (
            vec![LubricationSurface3::Prescribed(geometry)],
            vec![[0.0; 6]],
            vec![geometry],
        ),
        (
            vec![LubricationSurface3::Finite(carrier)],
            vec![load],
            vec![finite_geometry],
        ),
        (Vec::new(), Vec::new(), Vec::new()),
    ] {
        let gate = CancelGate::new();
        let native = prepare(&surfaces, &external, &gate).unwrap();
        assert_eq!(native, legacy);
        let mut expected = [initial];
        let expected_report = solve(
            model(),
            [0.5],
            &legacy,
            &[initial],
            &[response],
            &mut expected,
            [0.0; 3],
            &gate,
        )
        .unwrap();
        let (after, film, reports) =
            evaluate(initial, response, &mut surfaces, &external, [0.0; 3], &gate).unwrap();
        assert_eq!(after, expected[0]);
        assert_eq!((film, reports), expected_report);
    }
}

#[test]
fn invalid_mixed_loads_and_cancelled_or_refused_solve_leave_carriers_unchanged() {
    let initial = body([3.0, 4.0, 4.0], 2.0, [0.0; 3]);
    let wall = body([2.25, 4.0, 4.0], 3.0, [0.0; 3]);
    let surfaces = [
        LubricationSurface3::Finite(
            FiniteLubricationPlane3::new([2.25, 0.0, 0.0], [1.0, 0.0, 0.0], wall, 0.75, 0.02)
                .unwrap(),
        ),
        LubricationSurface3::Prescribed(
            LubricationPlane3::new([3.75, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3], 0.75, 0.02)
                .unwrap(),
        ),
    ];
    let gate = CancelGate::new();
    assert!(prepare(&surfaces, &[[0.0; 6]], &gate).is_err());
    assert!(prepare(&[surfaces[0]; 65], &[[0.0; 6]; 65], &gate).is_err());
    for loads in [
        [[f64::NAN; 6], [0.0; 6]],
        [[0.0; 6], [0.0, 0.0, 0.0, 0.0, 1e-30, 0.0]],
        [[1000.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0; 6]],
    ] {
        let mut trial = surfaces;
        assert!(
            evaluate(
                initial,
                Response::default(),
                &mut trial,
                &loads,
                [0.0; 3],
                &gate
            )
            .is_err()
        );
        assert_eq!(trial, surfaces);
    }
    gate.request();
    let mut trial = surfaces;
    assert!(
        evaluate(
            initial,
            Response::default(),
            &mut trial,
            &[[0.0; 6]; 2],
            [0.0; 3],
            &gate
        )
        .is_err()
    );
    assert_eq!(trial, surfaces);
    let fresh = CancelGate::new();
    let retried = evaluate(
        initial,
        Response::default(),
        &mut trial,
        &[[0.0; 6]; 2],
        [0.0; 3],
        &fresh,
    )
    .unwrap();
    let mut untouched = surfaces;
    let direct = evaluate(
        initial,
        Response::default(),
        &mut untouched,
        &[[0.0; 6]; 2],
        [0.0; 3],
        &fresh,
    )
    .unwrap();
    assert_eq!(retried, direct);
    assert_eq!(trial, untouched);
}
