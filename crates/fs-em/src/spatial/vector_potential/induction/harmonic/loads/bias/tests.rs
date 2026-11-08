use super::*;
use fs_ga::{Se3, Vec3};
use fs_math::det;

#[test]
fn bias_carrier_preserves_peak_force_torque_and_all_p1_virtual_work() {
    let positions = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 4.0],
    ];
    let volume = 4.0;
    let bias = [0.7, -0.3, 1.1];
    let current = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            C64::new(0.4 + i as f64 - 0.7 * j as f64, -0.3 + 0.5 * (i + j) as f64)
        })
    });
    let density = force_density(current, bias);
    let limits = SpatialLorentzLoadLimits {
        max_vertices: 4,
        max_cells: 1,
    };
    let loads = density.map(|d| {
        integrate(
            &positions,
            [CellLoad {
                vertices: [0, 1, 2, 3],
                volume_m3: volume,
                force_density_n_m3: &d,
            }]
            .into_iter(),
            limits,
            &mut || false,
        )
        .unwrap()
    });
    let velocity = positions.map(|p| [0.2 + p[0] + 0.4 * p[2], -0.3 + p[1], 0.8 - 0.5 * p[0]]);
    let wrenches = loads
        .each_ref()
        .map(|load| load.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array());
    let powers = loads
        .each_ref()
        .map(|load| load.power_w(&velocity, &mut || false).unwrap());
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    for phase in [0.0, 0.19, 1.2, 2.7, 4.8] {
        let (sine, cosine) = det::sin_cos(phase);
        let mut nodal = [[0.0; 3]; 4];
        let mut wrench = [0.0; 6];
        let mut power = 0.0;
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|i| if point == i { high } else { low });
            let j: [f64; 3] = std::array::from_fn(|a| {
                (0..4)
                    .map(|i| shape[i] * (current[i][a].re * cosine + current[i][a].im * sine))
                    .sum()
            });
            let p: [f64; 3] =
                std::array::from_fn(|a| (0..4).map(|i| shape[i] * positions[i][a]).sum());
            let v: [f64; 3] =
                std::array::from_fn(|a| (0..4).map(|i| shape[i] * velocity[i][a]).sum());
            let force = Vec3::new(j[0], j[1], j[2]).cross(Vec3::new(bias[0], bias[1], bias[2]));
            let moment = Vec3::new(p[0], p[1], p[2]).cross(force);
            let f = [force.x, force.y, force.z];
            let m = [moment.x, moment.y, moment.z];
            for a in 0..3 {
                wrench[a] += volume / 4.0 * m[a];
                wrench[3 + a] += volume / 4.0 * f[a];
                power += volume / 4.0 * f[a] * v[a];
                for i in 0..4 {
                    nodal[i][a] += volume / 4.0 * shape[i] * f[a];
                }
            }
        }
        for i in 0..4 {
            for a in 0..3 {
                let actual = cosine * loads[0].nodal_forces_n()[i][a]
                    + sine * loads[1].nodal_forces_n()[i][a];
                assert!((actual - nodal[i][a]).abs() < 2e-13);
            }
        }
        for a in 0..6 {
            assert!((cosine * wrenches[0][a] + sine * wrenches[1][a] - wrench[a]).abs() < 2e-13);
        }
        assert!((cosine * powers[0] + sine * powers[1] - power).abs() < 2e-13);
    }
}

#[test]
fn bias_force_is_linear_in_both_sources_and_rotates_at_the_fundamental() {
    let current = [[C64::new(1.2, -0.7), C64::new(-0.2, 0.8), C64::new(0.6, 1.1)]; 4];
    let b = [0.4, -0.9, 0.2];
    let base = force_density(current, b);
    let phase = C64::new(0.6, 0.8).scale(1.7);
    let changed = force_density(current.map(|j| j.map(|v| v * phase)), b.map(|v| -2.3 * v));
    for vertex in 0..4 {
        for axis in 0..3 {
            let expected =
                (C64::new(base[0][vertex][axis], base[1][vertex][axis]) * phase).scale(-2.3);
            assert!((changed[0][vertex][axis] - expected.re).abs() < 4e-15);
            assert!((changed[1][vertex][axis] - expected.im).abs() < 4e-15);
        }
    }
    assert_eq!(force_density(current, [0.0; 3]), [[[0.0; 3]; 4]; 2]);
}
