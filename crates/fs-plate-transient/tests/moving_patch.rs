//! G0/G3/G4: moving finite-area loads through the real spatial dynamics.
#[path = "../examples/moving_plate.rs"]
#[allow(dead_code)]
mod moving;

use fs_plate::{
    PlateMesh, PlateModel,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePointStencil},
};
use fs_plate_transient::{MovingPlateLoad, MovingPlateLoadBudget, PlateLoad, PlateLoadKeyframe};
use fs_time::galpha::{GeneralizedAlpha, SecondOrderState, galpha_step};

fn footprint() -> [PlatePatchPoint; 4] {
    [
        [-0.025, -0.02],
        [0.025, -0.02],
        [-0.025, 0.02],
        [0.025, 0.02],
    ]
    .map(|point| PlatePatchPoint {
        point,
        area_m2: 0.0005,
    })
}
fn keyframes() -> [PlateLoadKeyframe; 3] {
    [
        PlateLoadKeyframe {
            time_s: 0.0,
            center_m: [0.12, 0.27],
            load: [1.0, 0.002, -0.003],
        },
        PlateLoadKeyframe {
            time_s: moving::STEP * 6.0,
            center_m: [0.5, 0.32],
            load: [1.7, -0.004, 0.002],
        },
        PlateLoadKeyframe {
            time_s: moving::STEP * moving::STEPS as f64,
            center_m: [0.88, 0.27],
            load: [0.6, 0.003, 0.001],
        },
    ]
}
fn budget(mesh: &PlateMesh) -> MovingPlateLoadBudget {
    MovingPlateLoadBudget {
        max_keyframes: 3,
        patch: PlatePatchBudget {
            max_points: 4,
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:.16e} != {b:.16e}"
    );
}

// Independently assemble the specified quadrature through the existing point
// transfer, with explicit interpolation of this three-knot test path.
fn point_forcing(
    mesh: &PlateMesh,
    plate: &PlateModel,
    knots: &[PlateLoadKeyframe; 3],
    time: f64,
) -> Vec<f64> {
    let i = usize::from(time > knots[1].time_s);
    let t = (time - knots[i].time_s) / (knots[i + 1].time_s - knots[i].time_s);
    let center: [f64; 2] =
        std::array::from_fn(|j| (1.0 - t) * knots[i].center_m[j] + t * knots[i + 1].center_m[j]);
    let total: [f64; 3] =
        std::array::from_fn(|j| (1.0 - t) * knots[i].load[j] + t * knots[i + 1].load[j]);
    let mut result = vec![0.0; plate.free];
    for point in footprint() {
        PlatePointStencil::locate(
            mesh,
            plate,
            [center[0] + point.point[0], center[1] + point.point[1]],
            budget(mesh).patch.location,
        )
        .unwrap()
        .add_load(total.map(|v| v * 0.25), &mut result)
        .unwrap();
    }
    result
}

#[test]
fn moving_patch_interpolation_and_all_knot_load_derivatives_preserve_work() {
    let (mesh, plate) = moving::build_plate().unwrap();
    let points = footprint();
    let knots = keyframes();
    let load = MovingPlateLoad::new(&mesh, &plate, &points, &knots, budget(&mesh)).unwrap();
    assert_eq!(load.parameter_count(), 9);
    assert_eq!(load.horizon_s(), [knots[0].time_s, knots[2].time_s]);
    let q: Vec<_> = (0..plate.free).map(|i| (i % 5) as f64 * 0.01).collect();
    let v: Vec<_> = (0..plate.free)
        .map(|i| (i % 7) as f64 * 0.03 - 0.08)
        .collect();
    for time in [0.0, 0.004, knots[1].time_s, 0.019, knots[2].time_s] {
        let mut force = vec![99.0; plate.free];
        load.forcing(time, &mut force).unwrap();
        let expected = point_forcing(&mesh, &plate, &knots, time);
        for (&a, &b) in force.iter().zip(&expected) {
            close(a, b, 1e-13);
        }
        let (sample, power) = load.sample(time, &q, &v).unwrap();
        let (_, physical) = load.stencil(time).unwrap();
        close(dot(&force, &v), power, 1e-13);
        close(
            dot(&force, &q),
            physical[0] * sample.displacement + physical[1] * sample.slopes[1]
                - physical[2] * sample.slopes[0],
            1e-13,
        );
        let mut gradient = vec![99.0; 9];
        load.forcing_vjp(time, &v, &mut gradient).unwrap();
        for parameter in 0..9 {
            let mut plus = knots;
            let mut minus = knots;
            let delta = 1e-5;
            plus[parameter / 3].load[parameter % 3] += delta;
            minus[parameter / 3].load[parameter % 3] -= delta;
            let fd = (dot(&point_forcing(&mesh, &plate, &plus, time), &v)
                - dot(&point_forcing(&mesh, &plate, &minus, time), &v))
                / (2.0 * delta);
            close(gradient[parameter], fd, 1e-10);
        }
    }
    // Translating observation over a displaced, motionless plate supplies no power.
    let zero = vec![0.0; plate.free];
    assert_eq!(load.sample(0.004, &q, &zero).unwrap().1, 0.0);
    assert_eq!(load.sample(0.019, &q, &zero).unwrap().1, 0.0);
}

#[test]
fn refusals_are_atomic_and_repeated_or_reversed_time_queries_are_identical() {
    let (mesh, plate) = moving::build_plate().unwrap();
    let points = footprint();
    let knots = keyframes();
    let load = MovingPlateLoad::new(&mesh, &plate, &points, &knots, budget(&mesh)).unwrap();
    let mut force = vec![7.0; plate.free];
    let mut gradient = vec![8.0; 9];
    for time in [-1.0, 1.0, f64::NAN] {
        let before = force.clone();
        assert!(load.forcing(time, &mut force).is_err());
        assert_eq!(force, before);
        let before = gradient.clone();
        assert!(load.forcing_vjp(time, &force, &mut gradient).is_err());
        assert_eq!(gradient, before);
    }
    load.forcing(0.004, &mut force).unwrap();
    let first = force.clone();
    load.forcing(0.019, &mut force).unwrap();
    load.forcing(0.004, &mut force).unwrap();
    assert_eq!(force, first);
    let mut invalid = knots;
    invalid[1].time_s = invalid[0].time_s;
    assert!(MovingPlateLoad::new(&mesh, &plate, &points, &invalid, budget(&mesh)).is_err());
    assert!(MovingPlateLoad::new(&mesh, &plate, &[], &knots, budget(&mesh)).is_err());
    assert!(
        MovingPlateLoad::new(
            &mesh,
            &plate,
            &points,
            &knots,
            MovingPlateLoadBudget {
                max_keyframes: 2,
                ..budget(&mesh)
            }
        )
        .is_err()
    );

    // Both knot footprints fit, but the translating footprint crosses a hole.
    let mut cut_mesh = mesh.clone();
    cut_mesh.tris.retain(|tri| {
        let x = tri.iter().map(|&node| cut_mesh.nodes[node].0).sum::<f64>() / 3.0;
        !(0.25..0.75).contains(&x)
    });
    let crossing = [knots[0], knots[2]];
    let load =
        MovingPlateLoad::new(&cut_mesh, &plate, &points, &crossing, budget(&cut_mesh)).unwrap();
    let before = force.clone();
    assert!(load.forcing(knots[1].time_s, &mut force).is_err());
    assert_eq!(force, before);
    let mut gradient = vec![8.0; 6];
    let before = gradient.clone();
    assert!(
        load.forcing_vjp(knots[1].time_s, &force, &mut gradient)
            .is_err()
    );
    assert_eq!(gradient, before);
}

// Independent dense time stepping and point-quadrature forcing: no moving
// patch adapter, matrix-free PlateDynamics action or time-adjoint is used.
fn dense_run(
    mesh: &PlateMesh,
    plate: &PlateModel,
    knots: &[PlateLoadKeyframe; 3],
) -> (SecondOrderState, f64) {
    let n = plate.free;
    let m: Vec<_> = plate.m.to_dense().iter().map(|v| 0.9 * v).collect();
    let k: Vec<_> = plate.k.to_dense().iter().map(|v| 1.1 * v).collect();
    let c: Vec<_> = m.iter().zip(&k).map(|(m, k)| 3.0 * m + 2e-4 * k).collect();
    let dense = GeneralizedAlpha::new(&m, &c, &k, n, moving::STEP, moving::RHO);
    let force = point_forcing(mesh, plate, knots, 0.0);
    let mut state = SecondOrderState::new(
        0.0,
        &vec![0.0; n],
        &vec![0.0; n],
        &(0..n).map(|i| force[i] / m[i * n + i]).collect::<Vec<_>>(),
    );
    let weights = moving::sensor_weights(mesh, plate).unwrap();
    let alpha_f = moving::RHO / (1.0 + moving::RHO);
    let mut losses = vec![0.0];
    for step in 1..=moving::STEPS {
        let force = point_forcing(
            mesh,
            plate,
            knots,
            moving::STEP.mul_add(1.0 - alpha_f, state.t),
        );
        galpha_step(&dense, &mut state.q, &mut state.v, &mut state.a, &force);
        state.t += moving::STEP;
        state.steps += 1;
        if moving::SAMPLES.contains(&step) {
            losses.push(0.5 * (dot(&weights, &state.q) / moving::DEFLECTION_SCALE_M).powi(2));
        }
    }
    (state, losses.into_iter().rev().sum())
}

#[test]
fn actual_patch_dynamics_and_initialization_adjoint_match_dense_quadrature() {
    let (mesh, plate) = moving::build_plate().unwrap();
    let points = footprint();
    let knots = keyframes();
    let load = MovingPlateLoad::new(&mesh, &plate, &points, &knots, budget(&mesh)).unwrap();
    let actual = moving::evaluate_load(
        &mesh,
        &plate,
        moving::dynamics_parameters([1.1, 0.9, 3.0, 2e-4, 0.0]),
        &load,
    )
    .unwrap();
    assert_eq!(actual.gradient.len(), 13);
    let (expected, loss) = dense_run(&mesh, &plate, &knots);
    close(actual.loss, loss, 2e-10);
    for (a, b) in [
        (&actual.final_state.q, &expected.q),
        (&actual.final_state.v, &expected.v),
        (&actual.final_state.a, &expected.a),
    ] {
        let scale = b.iter().fold(1e-8_f64, |m, v| m.max(v.abs()));
        for (&a, &b) in a.iter().zip(b) {
            assert!((a - b).abs() < 2e-8 * scale);
        }
    }
    for parameter in 0..9 {
        let mut plus = knots;
        let mut minus = knots;
        let delta = if parameter % 3 == 0 { 2e-5 } else { 2e-7 };
        plus[parameter / 3].load[parameter % 3] += delta;
        minus[parameter / 3].load[parameter % 3] -= delta;
        let fd = (dense_run(&mesh, &plate, &plus).1 - dense_run(&mesh, &plate, &minus).1)
            / (2.0 * delta);
        assert!(
            (actual.gradient[4 + parameter] - fd).abs() < 5e-5 * (1e-8 + fd.abs()),
            "knot load {parameter}: adjoint {}, dense FD {fd}",
            actual.gradient[4 + parameter]
        );
    }
}
