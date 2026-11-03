use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn reference([x, y, z]: [f64; 3]) -> [f64; 3] {
    [
        2.0 + x + 0.3 * y * z,
        -1.0 + 1.2 * y + 0.2 * x * z,
        0.5 + z + 0.4 * x * y,
    ]
}
fn displacement([x, y, z]: [f64; 3]) -> [f64; 3] {
    [0.1 * x * x, -0.15 * y * y, 0.08 * z * z]
}
fn velocity([x, y, z]: [f64; 3]) -> [f64; 3] {
    [
        0.4 + x * y + 0.2 * z * z,
        -0.3 + 0.7 * y * z + 0.1 * x * x,
        0.9 - 0.5 * x * z,
    ]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-14,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn curved_quadratic_trace_and_transpose_preserve_impulse_torque_and_work() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        // Independently tabulated natural nodes in the public canonical edge
        // order. The polynomial fields make a P1 or clamped-weight trace fail.
        let natural = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let nodes = natural.map(reference);
        let u = natural.map(displacement);
        let v = natural.map(velocity);
        let cells = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]];
        let point = Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric: [0.4, 0.2, 0.15, 0.25],
        };
        let p = [0.2, 0.15, 0.25];
        let expected_x: [f64; 3] =
            std::array::from_fn(|axis| reference(p)[axis] + displacement(p)[axis]);
        let expected_v = velocity(p);
        let trace = kinematics(point, &cells, &nodes, &u, &v, &cx).unwrap();
        for axis in 0..3 {
            close(trace.position_m[axis], expected_x[axis]);
            close(trace.velocity_mps[axis], expected_v[axis]);
        }
        let impulse = [1.3, -0.7, 0.45];
        let entries = distribution(point, &cells, nodes.len(), impulse, &cx).unwrap();
        assert!((0..4).all(|i| entries[3 * i].1 < 0.0));
        assert!((4..10).all(|i| entries[3 * i].1 > 0.0));
        let origin = [-0.3, 0.6, 1.2];
        let mut resultant = [0.0; 3];
        let mut torque = [0.0; 3];
        let mut work = 0.0;
        for (dof, load) in entries {
            let node = dof / 3;
            let axis = dof % 3;
            resultant[axis] += load;
            work += load * v[node][axis];
            let arm = std::array::from_fn(|k| nodes[node][k] + u[node][k] - origin[k]);
            let mut force = [0.0; 3];
            force[axis] = load;
            let moment = cross(arm, force);
            for k in 0..3 {
                torque[k] += moment[k];
            }
        }
        let expected_torque = cross(std::array::from_fn(|k| expected_x[k] - origin[k]), impulse);
        for axis in 0..3 {
            close(resultant[axis], impulse[axis]);
            close(torque[axis], expected_torque[axis]);
        }
        close(
            work,
            (0..3).map(|axis| impulse[axis] * expected_v[axis]).sum(),
        );
        for barycentric in [[-0.1, 0.2, 0.4, 0.5], [0.1; 4], [f64::NAN, 0.0, 0.0, 1.0]] {
            let invalid = Tet10ImpulsePoint {
                barycentric,
                ..point
            };
            assert!(kinematics(invalid, &cells, &nodes, &u, &v, &cx).is_err());
            assert!(distribution(invalid, &cells, nodes.len(), impulse, &cx).is_err());
        }
    });
}
