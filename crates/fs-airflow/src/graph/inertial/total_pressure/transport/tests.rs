use super::*;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn cubic_pressure_work_matches_independent_quadrature_during_reversal() {
    let (rho, area, dt) = (2.0, 1.0, 0.5);
    let work = branch_work(dt, rho, area, [20.0, 17.0], [-1.0, 3.0]).unwrap();
    // Simpson is exact for this independently evaluated cubic p(t) Q(t).
    let integrate = |p: f64| {
        let f = |q: f64| (p - 0.5 * rho * (q / area).powi(2)) * q;
        dt * (f(-1.0) + 4.0 * f(1.0) + f(3.0)) / 6.0
    };
    near(
        work.static_pressure_work_into_branch_j[0],
        integrate(20.0),
        2e-14,
    );
    near(
        work.static_pressure_work_into_branch_j[1],
        -integrate(17.0),
        2e-14,
    );
    near(work.average_static_pressures_pa[0], 20.0 - 7.0 / 3.0, 2e-14);
    near(work.kinetic_energy_into_branch_j[0], 2.5, 2e-14);
    // Multiplying two interval averages would be the wrong work.
    assert!(
        (work.static_pressure_work_into_branch_j[0] - work.average_static_pressures_pa[0] * dt)
            .abs()
            > 1.0
    );
    for port in 0..2 {
        near(
            work.static_pressure_work_into_branch_j[port] + work.kinetic_energy_into_branch_j[port],
            work.total_pressure_work_into_branch_j[port],
            2e-14,
        );
    }
}

#[test]
fn actual_native_descriptor_closes_unequal_area_streams_and_total_energy() {
    use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, step_descriptor};
    let system = DescriptorPortHamiltonian::new(
        3,
        2,
        2,
        vec![0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0],
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        Box::new(QuadraticStorage::new(vec![0.5, 0.0, 0.0, 1.0 / 3.0], 2).unwrap()),
    )
    .unwrap();
    for inlet in [4.0, -30.0] {
        let dt = 0.2;
        let step = step_descriptor(&system, &[0.8, 1.2, 0.0], &[inlet, 0.0], dt).unwrap();
        let average = (5.0 * 0.4 + 0.5 * dt * inlet) / (5.0 + dt);
        let final_flow = 2.0 * average - 0.4;
        near(step.y[0], average, 2e-13);
        near(step.y[1], average, 2e-13);
        near(step.x[0] / 2.0, final_flow, 2e-13);
        near(step.x[1] / 3.0, final_flow, 2e-13);
        let node_total = 0.5 * step.x[2];
        let first = branch_work(dt, 2.0, 1.0, [inlet, node_total], [0.4, step.x[0] / 2.0]).unwrap();
        let second = branch_work(dt, 2.0, 0.5, [node_total, 0.0], [0.4, step.x[1] / 3.0]).unwrap();
        let static_out = first.static_pressure_work_into_branch_j[1]
            + second.static_pressure_work_into_branch_j[0];
        let kinetic_out =
            first.kinetic_energy_into_branch_j[1] + second.kinetic_energy_into_branch_j[0];
        near(static_out + kinetic_out, 0.0, 2e-12);
        // The reducer's different speeds require different STATIC port pressures.
        assert!(second.average_static_pressures_pa[0] < first.average_static_pressures_pa[1]);
        let external_pressure = first.static_pressure_work_into_branch_j[0]
            + second.static_pressure_work_into_branch_j[1];
        let external_kinetic =
            first.kinetic_energy_into_branch_j[0] + second.kinetic_energy_into_branch_j[1];
        near(
            step.delta_h + step.dissipated - external_pressure - external_kinetic,
            0.0,
            2e-12,
        );
        near(external_pressure + external_kinetic, step.supplied, 2e-12);
    }
}

#[test]
fn source_bounds_include_zero_crossing_and_same_branch_pressure_drop() {
    let bounds = static_pressure_bounds([50.0, 47.0], [-2.0, 2.0], 3.0, 0.5).unwrap();
    assert_eq!(bounds, [[26.0, 23.0], [50.0, 47.0]]);
    for i in 0..=256 {
        let q = -2.0 + i as f64 / 64.0;
        let head = 0.5 * 3.0 * (q / 0.5).powi(2);
        for port in 0..2 {
            let p = [50.0, 47.0][port] - head;
            assert!(p >= bounds[0][port] && p <= bounds[1][port]);
        }
    }
    for endpoint in bounds {
        assert_eq!(endpoint[0] - endpoint[1], 3.0);
    }
    let same_sign = static_pressure_bounds([50.0, 47.0], [-1.0, -2.0], 3.0, 0.5).unwrap();
    assert_eq!(same_sign, [[26.0, 23.0], [44.0, 41.0]]);
}

#[test]
fn affine_partition_gauge_and_sign_preserve_physical_stream_work() {
    let full = branch_work(0.6, 3.0, 0.7, [8.0, 6.0], [2.0, -1.0]).unwrap();
    let a = branch_work(0.2, 3.0, 0.7, [8.0, 6.0], [2.0, 1.0]).unwrap();
    let b = branch_work(0.4, 3.0, 0.7, [8.0, 6.0], [1.0, -1.0]).unwrap();
    let backwards = branch_work(0.6, 3.0, 0.7, [8.0, 6.0], [-2.0, 1.0]).unwrap();
    let shifted = branch_work(0.6, 3.0, 0.7, [10008.0, 10006.0], [2.0, -1.0]).unwrap();
    for port in 0..2 {
        near(
            a.static_pressure_work_into_branch_j[port] + b.static_pressure_work_into_branch_j[port],
            full.static_pressure_work_into_branch_j[port],
            3e-14,
        );
        near(
            a.kinetic_energy_into_branch_j[port] + b.kinetic_energy_into_branch_j[port],
            full.kinetic_energy_into_branch_j[port],
            3e-14,
        );
        assert_eq!(
            backwards.average_static_pressures_pa[port],
            full.average_static_pressures_pa[port]
        );
        near(
            backwards.static_pressure_work_into_branch_j[port],
            -full.static_pressure_work_into_branch_j[port],
            3e-14,
        );
    }
    near(
        shifted
            .static_pressure_work_into_branch_j
            .iter()
            .sum::<f64>(),
        full.static_pressure_work_into_branch_j.iter().sum::<f64>(),
        1e-12,
    );
    assert!(branch_work(0.0, 3.0, 0.7, [8.0; 2], [1.0; 2]).is_err());
    assert!(branch_work(1.0, 3.0, 0.7, [f64::NAN; 2], [1.0; 2]).is_err());
    assert!(branch_work(1.0, 3.0, 0.7, [8.0; 2], [f64::MAX; 2]).is_err());
}
