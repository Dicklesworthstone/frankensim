use super::super::pipe::{FlowInterval, PipeLoadGeometry};
use super::*;
use fs_qty::Area;

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a} != {b}");
}
fn vector_near(a: [f64; 3], b: [f64; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        near(a, b);
    }
}
fn state() -> JunctionState {
    JunctionState {
        node: 1,
        position_m: [0.0; 3],
        about_m: [0.0; 3],
        duration_s: 1.0,
        density_kg_m3: 1.0,
        pressure_pa: 4.0,
        exterior_pressure_pa: 0.0,
    }
}
fn ports(final_q: [f64; 3], areas: [f64; 3]) -> [PortInterval; 3] {
    [
        PortInterval {
            outward_normal: [-1.0, 0.0, 0.0],
            area_m2: areas[0],
            outflows_m3_s: [0.0, -final_q[0]],
        },
        PortInterval {
            outward_normal: [0.0, 1.0, 0.0],
            area_m2: areas[1],
            outflows_m3_s: [0.0, final_q[1]],
        },
        PortInterval {
            outward_normal: [1.0, 0.0, 0.0],
            area_m2: areas[2],
            outflows_m3_s: [0.0, final_q[2]],
        },
    ]
}
fn geometry(points: Vec<[f64; 3]>, area: f64) -> PipeLoadGeometry {
    PipeLoadGeometry::try_new(points, Area::new(area), 8, || false).unwrap()
}
fn pipe_interval(q: [f64; 2], p: [f64; 2], exposed_ports: [bool; 2]) -> FlowInterval {
    FlowInterval {
        duration_s: 1.0,
        volume_flows_m3_s: q,
        density_kg_m3: 1.0,
        pressures_pa: p,
        exterior_pressure_pa: 0.0,
        exposed_ports,
        gravity_m_s2: [0.0; 3],
    }
}

/// Original descriptor incidence/storage step, with one zero-storage pressure
/// multiplier and three native branch impulses. No replacement time integrator.
fn native_startup() -> ([f64; 3], f64) {
    use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, step_descriptor};
    let mut j = vec![0.0; 16];
    for (branch, sign) in [(0, -1.0), (1, 1.0), (2, 1.0)] {
        j[branch * 4 + 3] = sign;
        j[3 * 4 + branch] = -sign;
    }
    let mut r = vec![0.0; 16];
    let mut g = vec![0.0; 12];
    for (branch, resistance) in [1.0, 2.0, 2.0].into_iter().enumerate() {
        r[branch * 4 + branch] = resistance;
        g[branch * 3 + branch] = 1.0;
    }
    let h = QuadraticStorage::new(vec![1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5], 3).unwrap();
    let system = DescriptorPortHamiltonian::new(4, 3, 3, j, r, g, Box::new(h)).unwrap();
    let result = step_descriptor(&system, &[0.0; 4], &[8.0, 0.0, 0.0], 1.0).unwrap();
    near(result.balance_residual(), 0.0);
    let final_q = [result.x[0], result.x[1] / 2.0, result.x[2] / 2.0];
    for (q, expected) in final_q.into_iter().zip([8.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0]) {
        near(q, expected);
    }
    near(0.5 * result.x[3], 4.0);
    (final_q, 0.5 * result.x[3])
}

#[test]
fn native_startup_turning_and_bent_pipe_loads_drive_real_mechanical_momentum() {
    use fs_ga::{Mat3, Se3, Vec3, Wrench};
    use fs_mbd::articulated::{
        ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
        apply_free_floating_impulses,
    };
    let (q, pressure) = native_startup();
    let mut junction_state = state();
    junction_state.pressure_pa = pressure;
    let paths = [
        geometry(vec![[-2.0, 0.0, 0.0], [0.0; 3]], 2.0),
        geometry(vec![[0.0; 3], [0.0, 2.0, 0.0]], 1.0),
        geometry(vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], 1.0),
    ];
    let actual_ports = std::array::from_fn::<_, 3, _>(|branch| {
        let (position, outward_normal) = paths[branch].junction_port(branch != 0);
        assert_eq!(position, [0.0; 3]);
        PortInterval {
            outward_normal,
            area_m2: paths[branch].area_m2,
            outflows_m3_s: [0.0, if branch == 0 { -q[branch] } else { q[branch] }],
        }
    });
    let junction = evaluate(junction_state, &actual_ports, || false).unwrap();
    vector_near(junction.load.force_n, [124.0 / 27.0, -124.0 / 27.0, 0.0]);
    near(junction.kinetic_energy_transport_j, 0.0);
    let pressures = [[8.0, pressure], [pressure, 0.0], [pressure, 0.0]];
    let mut components = Vec::new();
    let mut external = Vec::new();
    for branch in 0..3 {
        components.push(
            paths[branch]
                .evaluate(
                    pipe_interval([0.0, q[branch]], pressures[branch], [true; 2]),
                    [0.0; 3],
                    |_| 1.0,
                    || false,
                )
                .unwrap(),
        );
        external.push(
            paths[branch]
                .evaluate(
                    pipe_interval(
                        [0.0, q[branch]],
                        pressures[branch],
                        if branch == 0 {
                            [true, false]
                        } else {
                            [false, true]
                        },
                    ),
                    [0.0; 3],
                    |_| 1.0,
                    || false,
                )
                .unwrap(),
        );
    }
    components.push(junction.load);
    let summed = sum_loads(1.0, [0.0; 3], &components).unwrap();
    let total = sum_loads(1.0, [0.0; 3], &external).unwrap();
    vector_near(summed.impulse_ns, total.impulse_ns);
    vector_near(summed.angular_impulse_nms, total.angular_impulse_nms);
    // Independent whole-volume balance: exposed endpoint traction/flux minus
    // retained fluid momentum; the off-axis bend retains a genuine spin impulse.
    vector_near(total.impulse_ns, [284.0 / 27.0, -140.0 / 27.0, 0.0]);
    vector_near(total.angular_impulse_nms, [0.0, 0.0, -52.0 / 27.0]);
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "hydraulic support",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    let response = apply_free_floating_impulses(
        &model,
        FreeFloatingBaseState::stationary(Se3::identity()),
        &[],
        &[],
        &[],
        &[Wrench::new(
            v(total.angular_impulse_nms),
            v(total.impulse_ns),
        )],
    )
    .unwrap();
    let array = |v: Vec3| [v.x, v.y, v.z];
    vector_near(
        array(response.linear_momentum_change_world),
        total.impulse_ns,
    );
    vector_near(
        array(response.angular_momentum_change_world),
        total.angular_impulse_nms,
    );
    near(response.impulse_work, 42224.0 / 1215.0);
    near(response.kinetic_work_residual, 0.0);
}

#[test]
fn unequal_areas_expose_signed_energy_transport_and_reversal_flux() {
    let (q, _) = native_startup();
    let original = ports(q, [2.0, 0.5, 1.5]);
    let forward = evaluate(state(), &original, || false).unwrap();
    near(forward.kinetic_energy_transport_j, 176.0 / 243.0);
    vector_near(forward.load.force_n, [226.0 / 81.0, -86.0 / 27.0, 0.0]);
    let backward = original.map(|mut port| {
        port.outflows_m3_s = [-port.outflows_m3_s[1], 0.0];
        port
    });
    let reverse = evaluate(state(), &backward, || false).unwrap();
    near(reverse.kinetic_energy_transport_j, -176.0 / 243.0);
    vector_near(reverse.load.impulse_ns, forward.load.impulse_ns);
    let crossing = original.map(|mut port| {
        port.outflows_m3_s = [-port.outflows_m3_s[1], port.outflows_m3_s[1]];
        port
    });
    let mut unpressurized = state();
    unpressurized.pressure_pa = 0.0;
    let crossed = evaluate(unpressurized, &crossing, || false).unwrap();
    near(crossed.kinetic_energy_transport_j, 0.0);
    vector_near(crossed.load.impulse_ns, [64.0 / 81.0, -32.0 / 27.0, 0.0]);
    // Signed cubic transport and quadratic impulse are additive on the SAME
    // affine path; this does not assert equivalence of separately solved steps.
    let a = original.map(|mut p| {
        p.outflows_m3_s[1] *= 0.5;
        p
    });
    let b = original.map(|mut p| {
        p.outflows_m3_s[0] = 0.5 * p.outflows_m3_s[1];
        p
    });
    let mut half = state();
    half.duration_s = 0.5;
    let a = evaluate(half, &a, || false).unwrap();
    let b = evaluate(half, &b, || false).unwrap();
    near(
        a.kinetic_energy_transport_j + b.kinetic_energy_transport_j,
        forward.kinetic_energy_transport_j,
    );
    vector_near(
        super::super::pipe::add(a.load.impulse_ns, b.load.impulse_ns),
        forward.load.impulse_ns,
    );
}

#[test]
fn closed_loop_spinup_has_only_the_retained_fluid_reaction() {
    let paths = [
        geometry(vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], 1.0),
        geometry(vec![[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0; 3]], 1.0),
    ];
    // No cut is external. Even enormous internal static pressure must not
    // contaminate the finite net torque -rho*dQ/dt*integral(r cross dr).
    let contributions: Vec<_> = paths
        .iter()
        .map(|pipe| {
            pipe.evaluate(
                pipe_interval([0.0, 1.0], [1e30; 2], [false; 2]),
                [0.0; 3],
                |_| 1.0,
                || false,
            )
            .unwrap()
        })
        .collect();
    let total = sum_loads(1.0, [0.0; 3], &contributions).unwrap();
    vector_near(total.force_n, [0.0; 3]);
    vector_near(total.torque_nm, [0.0, 0.0, -2.0]);
}

#[test]
fn pressure_datum_arbitrary_origin_and_incomplete_stream_refusal() {
    let (q, _) = native_startup();
    let ports = ports(q, [2.0, 1.0, 1.0]);
    let mut moved = state();
    moved.about_m = [2.0, -3.0, 4.0];
    moved.pressure_pa += 1e5;
    moved.exterior_pressure_pa += 1e5;
    let load = evaluate(moved, &ports, || false).unwrap();
    vector_near(
        load.load.torque_nm,
        [-496.0 / 27.0, -496.0 / 27.0, -124.0 / 27.0],
    );
    assert!(evaluate(state(), &ports[..2], || false).is_err());
    let unrepresented_reversal = [PortInterval {
        outward_normal: [1.0, 0.0, 0.0],
        area_m2: 1.0,
        outflows_m3_s: [-1.0, 1.0],
    }];
    assert!(evaluate(state(), &unrepresented_reversal, || false).is_err());
    let mut invalid = ports;
    invalid[0].area_m2 = 0.0;
    assert!(evaluate(state(), &invalid, || false).is_err());
    assert_eq!(
        evaluate(state(), &ports, || true),
        Err(PipeLoadError::Cancelled)
    );
}

#[test]
fn total_pressure_reducer_closes_work_and_drives_native_mechanical_momentum() {
    use fs_ga::{Mat3, Se3, Vec3, Wrench};
    use fs_mbd::articulated::{
        ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
        apply_free_floating_impulses,
    };
    let (q, pressure) = native_startup();
    let paths = [
        geometry(vec![[-2.0, 0.0, 0.0], [0.0; 3]], 2.0),
        geometry(vec![[0.0; 3], [0.0, 1.0, 0.0]], 0.5),
        geometry(vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], 1.0),
    ];
    // The same original L=(1,2,2) now has genuinely different port speeds.
    let actual_ports = ports(q, [2.0, 0.5, 1.0]);
    let mut node = state();
    node.pressure_pa = pressure;
    let junction = evaluate_in_chart(
        node,
        &actual_ports,
        super::super::super::PressureChart::Total,
        || false,
    )
    .unwrap();
    vector_near(junction.load.force_n, [116.0 / 27.0, -70.0 / 27.0, 0.0]);
    near(junction.kinetic_energy_transport_j, 8.0 / 9.0);
    // Independent cubic quadrature of pressure work on the point volume.
    let supplied: f64 = actual_ports
        .iter()
        .map(|port| {
            let f = |q: f64| -(pressure - 0.5 * (q / port.area_m2).powi(2)) * q;
            let q = port.outflows_m3_s[1];
            (f(0.0) + 4.0 * f(0.5 * q) + f(q)) / 6.0
        })
        .sum();
    near(supplied, junction.kinetic_energy_transport_j);
    let total_pressures = [[8.0, pressure], [pressure, 0.0], [pressure, 0.0]];
    let mut components = Vec::new();
    let mut external = Vec::new();
    for branch in 0..3 {
        // Independent integral of the affine dynamic pressure for traction.
        let speed = q[branch] / paths[branch].area_m2;
        let head = speed * speed / 6.0;
        let pressures = total_pressures[branch].map(|p| p - head);
        let flow = [0.0, q[branch]];
        components.push(
            paths[branch]
                .evaluate(
                    pipe_interval(flow, pressures, [true; 2]),
                    [0.0; 3],
                    |_| 1.0,
                    || false,
                )
                .unwrap(),
        );
        external.push(
            paths[branch]
                .evaluate(
                    pipe_interval(
                        flow,
                        pressures,
                        if branch == 0 {
                            [true, false]
                        } else {
                            [false, true]
                        },
                    ),
                    [0.0; 3],
                    |_| 1.0,
                    || false,
                )
                .unwrap(),
        );
    }
    components.push(junction.load);
    let summed = sum_loads(1.0, [0.0; 3], &components).unwrap();
    let total = sum_loads(1.0, [0.0; 3], &external).unwrap();
    vector_near(summed.impulse_ns, total.impulse_ns);
    vector_near(summed.angular_impulse_nms, total.angular_impulse_nms);
    vector_near(total.impulse_ns, [268.0 / 27.0, -96.0 / 27.0, 0.0]);
    vector_near(total.angular_impulse_nms, [0.0, 0.0, -44.0 / 27.0]);
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "total-pressure support",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    let result = apply_free_floating_impulses(
        &model,
        FreeFloatingBaseState::stationary(Se3::identity()),
        &[],
        &[],
        &[],
        &[Wrench::new(
            v(total.angular_impulse_nms),
            v(total.impulse_ns),
        )],
    )
    .unwrap();
    let a = |v: Vec3| [v.x, v.y, v.z];
    vector_near(a(result.linear_momentum_change_world), total.impulse_ns);
    vector_near(
        a(result.angular_momentum_change_world),
        total.angular_impulse_nms,
    );
    near(result.impulse_work, 102268.0 / 3645.0);
    near(result.kinetic_work_residual, 0.0);
    assert_eq!(
        evaluate_in_chart(
            node,
            &actual_ports,
            super::super::super::PressureChart::Total,
            || true
        ),
        Err(PipeLoadError::Cancelled)
    );
}
