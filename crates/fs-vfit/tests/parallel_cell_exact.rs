//! Physical endpoint, work and decay regressions for exact held-flow loads.

use core::f64::consts::TAU;
use fs_math::c64::C64;
use fs_vfit::relaxation::{
    ExactParallelCellLoad, ParallelCell, ParallelCellLoad, ParallelCellSpec, ParallelLcCell,
};

fn cell(r: f64, l: f64, c: f64) -> ParallelCell {
    ParallelCell {
        resistance_pa_s_m3: r,
        inertance_pa_s2_m3: l,
        compliance_m3_pa: c,
    }
}

fn load(cells: &[ParallelCell], dt: f64) -> ExactParallelCellLoad {
    ExactParallelCellLoad::new(ParallelCellSpec::new(cells).unwrap(), dt).unwrap()
}

fn close(found: f64, expected: f64, scale: f64) {
    assert!(
        (found - expected).abs() <= 3e-11 * scale,
        "found {found:.17e}, expected {expected:.17e}, scale {scale:.17e}"
    );
}

#[test]
fn a_stiff_rl_cell_relaxes_physically_instead_of_retaining_an_alternating_tail() {
    let cell = cell(10_000.0, 0.2, 0.0);
    let dt = 0.02; // 1000 physical decay times in one interval.
    let spec = ParallelCellSpec::new(&[cell]).unwrap();
    let mut exact = ExactParallelCellLoad::new(spec, dt).unwrap();
    let charged = exact.step(0.3).unwrap();
    close(charged.inductor_flows_m3_s()[0], 0.3, 0.3);
    close(charged.stored_energy_j, 0.009, 0.009);
    close(charged.dissipated_energy_j, 0.009, 0.009);
    close(charged.supplied_work_j, 0.018, 0.018);
    assert_eq!(charged.pressures_pa(), &[0.0]);
    let released = exact.step(0.0).unwrap();
    assert_eq!(released.inductor_flows_m3_s(), &[0.0]);
    assert_eq!(released.pressures_pa(), &[-0.0]);
    close(
        released.dissipated_energy_j,
        charged.stored_energy_j,
        charged.stored_energy_j,
    );
    assert_eq!(released.supplied_work_j, 0.0);
    close(released.balance_residual_j(), 0.0, charged.stored_energy_j);

    let mut midpoint = ParallelCellLoad::new(spec, dt).unwrap();
    let before = midpoint.step(0.3).unwrap().inductor_flows_m3_s()[0];
    let after = midpoint.step(0.0).unwrap().inductor_flows_m3_s()[0];
    assert!(
        after < -0.99 * before,
        "fixture must expose the unresolved midpoint tail"
    );
}

// Independent closed form of L i' = p, C p' = U-i-p/R from zero storage.
fn analytic_step(cell: ParallelCell, u: f64, t: f64) -> (f64, f64) {
    let omega2 = 1.0 / (cell.inertance_pa_s2_m3 * cell.compliance_m3_pa);
    let alpha = 0.5 / (cell.resistance_pa_s_m3 * cell.compliance_m3_pa);
    let discriminant = alpha * alpha - omega2;
    let (retention, impulse) = if discriminant == 0.0 {
        let e = (-alpha * t).exp();
        (e * (1.0 + alpha * t), e * t)
    } else if discriminant < 0.0 {
        let w = (-discriminant).sqrt();
        let e = (-alpha * t).exp();
        (
            e * ((w * t).cos() + alpha / w * (w * t).sin()),
            e * (w * t).sin() / w,
        )
    } else {
        let root = discriminant.sqrt();
        let (r1, r2) = (-alpha + root, -alpha - root);
        let (e1, e2) = ((r1 * t).exp(), (r2 * t).exp());
        ((r1 * e2 - r2 * e1) / (r1 - r2), (e1 - e2) / (r1 - r2))
    };
    (u * (1.0 - retention), u * impulse / cell.compliance_m3_pa)
}

#[test]
fn rlc_endpoints_and_independent_resistor_heat_match_all_damping_regimes() {
    for (cell, dt) in [
        (cell(1.0, 1.0, 1.0), 0.6),
        (cell(0.5, 1.0, 1.0), 0.6),
        (cell(0.2, 1.0, 1.0), 0.6),
        (cell(3.2, 0.0037, 0.00012), 0.0002),
        (cell(3.2e6, 3.7e4, 1.2e-10), 0.0002),
    ] {
        let u = 0.7;
        let frame = load(&[cell], dt).step(u).unwrap();
        let (i, p) = analytic_step(cell, u, dt);
        close(frame.inductor_flows_m3_s()[0], i, i.abs());
        close(frame.pressures_pa()[0], p, p.abs());
        close(
            frame.pressure_pa,
            cell.inertance_pa_s2_m3 * i / dt,
            frame.pressure_pa.abs(),
        );
        // Simpson integration of the analytic resistor pressure is independent
        // of both the modal heat kernel and work-minus-storage subtraction.
        let intervals = 2000;
        let mut integral = 0.0;
        for n in 0..=intervals {
            let (_, pressure) = analytic_step(cell, u, dt * f64::from(n) / f64::from(intervals));
            let weight = if n == 0 || n == intervals {
                1.0
            } else if n % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * (pressure / cell.resistance_pa_s_m3) * pressure;
        }
        let heat = integral * dt / (3.0 * f64::from(intervals));
        close(frame.dissipated_energy_j, heat, heat.abs());
        close(frame.balance_residual_j(), 0.0, frame.supplied_work_j.abs());
        assert!(frame.dissipated_energy_j >= 0.0);
    }
}

#[test]
fn a_held_flow_can_be_subdivided_without_changing_physical_decay_or_heat() {
    let cells = [
        cell(5000.0, 0.2, 0.0),
        cell(0.5, 1.0, 1.0),
        cell(3.0, 0.2, 0.1),
    ];
    let dt = 0.2;
    let mut coarse = load(&cells, dt);
    let mut fine = load(&cells, dt / 100.0);
    for u in [0.2, -0.4, 0.1, 0.0, 0.0] {
        let (direct, history) = coarse.affine_pressure().unwrap();
        let frame = coarse.step(u).unwrap();
        close(frame.pressure_pa, direct * u + history, 1.0);
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..100 {
            let small = fine.step(u).unwrap();
            heat += small.dissipated_energy_j;
            work += small.supplied_work_j;
        }
        for (&a, &b) in coarse
            .inductor_flows_m3_s()
            .iter()
            .zip(fine.inductor_flows_m3_s())
        {
            close(a, b, 1.0);
        }
        for (&a, &b) in coarse.pressures_pa().iter().zip(fine.pressures_pa()) {
            close(a, b, 1.0);
        }
        close(
            frame.dissipated_energy_j,
            heat,
            frame.stored_energy_j + heat,
        );
        close(
            frame.supplied_work_j,
            work,
            frame.stored_energy_j + heat + work.abs(),
        );
    }
}

#[test]
fn rest_static_flow_and_absent_resistors_have_their_physical_limits() {
    let mut cells = load(&[cell(1.0, 1.0, 0.0), cell(0.5, 1.0, 1.0)], 0.5);
    let rest = cells.step(0.0).unwrap();
    assert_eq!(rest.stored_energy_j, 0.0);
    assert_eq!(rest.dissipated_energy_j, 0.0);
    assert_eq!(cells.sampled_impedance(0.0).unwrap(), C64::ZERO);
    for _ in 0..300 {
        cells.step(0.25).unwrap();
    }
    let steady = cells.step(0.25).unwrap();
    for &i in steady.inductor_flows_m3_s() {
        close(i, 0.25, 0.25);
    }
    close(steady.pressure_pa, 0.0, 0.25);
    assert!(steady.dissipated_energy_j < 1e-27);
    let spec = ParallelCellSpec::lossless_lc(&[ParallelLcCell {
        inertance_pa_s2_m3: 1.0,
        compliance_m3_pa: 1.0,
    }])
    .unwrap();
    let frame = ExactParallelCellLoad::new(spec, 0.7)
        .unwrap()
        .step(0.4)
        .unwrap();
    close(
        frame.inductor_flows_m3_s()[0],
        0.4 * (1.0 - 0.7_f64.cos()),
        0.4,
    );
    close(frame.pressures_pa()[0], 0.4 * 0.7_f64.sin(), 0.4);
    assert_eq!(frame.dissipated_energy_j, 0.0);
    close(frame.balance_residual_j(), 0.0, frame.supplied_work_j);
    assert!(
        ParallelCellSpec::new(&[cell(0.0, 1.0, 1.0)]).is_err(),
        "a shorted cell is outside the positive-R declaration"
    );
}

#[test]
fn sampled_impedance_matches_the_measured_mean_pressure_transfer() {
    let cells = [
        cell(0.6, 0.3, 0.0),
        cell(0.9, 0.05, 0.4),
        cell(2.5, 0.01, 0.02),
    ];
    let dt = 0.02;
    let count = 1024;
    for bin in [1, 7, 43, 377] {
        let mut runtime = load(&cells, dt);
        let angle = TAU * f64::from(bin) / f64::from(count);
        for n in 0..8192 {
            runtime.step((angle * f64::from(n)).cos()).unwrap();
        }
        let mut measured = C64::ZERO;
        for n in 8192..8192 + count {
            let phase = angle * f64::from(n);
            let pressure = runtime.step(phase.cos()).unwrap().pressure_pa;
            measured = measured
                + C64::new(phase.cos(), -phase.sin()).scale(2.0 * pressure / f64::from(count));
        }
        let expected = runtime
            .sampled_impedance(f64::from(bin) / (f64::from(count) * dt))
            .unwrap();
        close(measured.re, expected.re, expected.abs());
        close(measured.im, expected.im, expected.abs());
        assert!(expected.re >= 0.0);
    }
}

#[test]
fn sampled_queries_keep_small_resistance_and_refuse_unresolved_lossless_poles() {
    let tiny = load(&[cell(1e-308, 1e308, 1.0)], 1.0);
    let direct = tiny.affine_pressure().unwrap().0;
    assert!(direct > 0.0);
    let sampled = tiny.sampled_impedance(0.25).unwrap();
    close(sampled.re / direct, 1.0, 1.0);
    close(sampled.im / direct, 0.0, 1.0);
    let spec = ParallelCellSpec::lossless_lc(&[ParallelLcCell {
        inertance_pa_s2_m3: 1.0,
        compliance_m3_pa: 1.0,
    }])
    .unwrap();
    let lossless = ExactParallelCellLoad::new(spec, 0.7).unwrap();
    assert!(lossless.sampled_impedance(1.0 / TAU).is_err());
    let away = lossless.sampled_impedance(0.1).unwrap();
    close(away.re, 0.0, away.abs());
}

#[test]
fn failed_exact_trials_preserve_every_accepted_cell() {
    let mut runtime = load(&[cell(2.0, 0.1, 0.0), cell(0.5, 1.0, 1.0)], 0.1);
    runtime.step(0.7).unwrap();
    let mut untouched = runtime.clone();
    assert!(runtime.step(f64::NAN).is_err());
    assert!(runtime.step(f64::MAX).is_err());
    assert_eq!(runtime.step(-0.1).unwrap(), untouched.step(-0.1).unwrap());
    assert!(runtime.sampled_impedance(-1.0).is_err());
    assert!(runtime.sampled_impedance(6.0).is_err());
}
