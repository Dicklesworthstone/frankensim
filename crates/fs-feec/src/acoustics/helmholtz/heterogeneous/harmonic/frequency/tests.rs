//! G1/G3: independent resonator/compliance formulas, a resolved pressure-soft
//! duct, coherent delayed sources, exact constraints and factor reuse.

use super::*;
use fs_la::eigen_complex::lu_complex_with_cancel;
use fs_math::det;
use fs_sparse::Coo;

struct Fixture {
    omega: f64,
    matrices: [Csr; 4], // K, M, G, B.
    free: Vec<usize>,
    scales: Vec<f64>,
    factor: LuComplex,
}

impl Fixture {
    fn new(omega: f64, matrices: [Csr; 4], free: Vec<usize>) -> Self {
        let n = free.len();
        let mut dense = vec![C64::ZERO; n * n];
        for (i, &row) in free.iter().enumerate() {
            for (matrix, multiplier, rotate) in [
                (&matrices[0], 1.0, false),
                (&matrices[1], -omega * omega, false),
                (&matrices[3], -omega, false),
                (&matrices[2], omega, true),
            ] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    if let Ok(j) = free.binary_search(&column) {
                        if rotate {
                            dense[i * n + j].im += multiplier * value;
                        } else {
                            dense[i * n + j].re += multiplier * value;
                        }
                    }
                }
            }
        }
        let scales: Vec<_> = dense
            .chunks_exact(n)
            .map(|row| 1.0 / det::sqrt(row.iter().map(|z| z.abs()).fold(0.0, f64::max)))
            .collect();
        for i in 0..n {
            for j in 0..n {
                dense[i * n + j] = dense[i * n + j].scale(scales[i]).scale(scales[j]);
            }
        }
        let factor = lu_complex_with_cancel(&dense, n, &mut || false).unwrap();
        Self {
            omega,
            matrices,
            free,
            scales,
            factor,
        }
    }

    fn pencil(&self) -> PressurePencil<'_> {
        PressurePencil {
            omega: self.omega,
            stiffness: &self.matrices[0],
            mass: &self.matrices[1],
            conductance: &self.matrices[2],
            susceptance: &self.matrices[3],
            factor: &self.factor,
            scales: &self.scales,
            free_dofs: &self.free,
            residual_tolerance: 5e-13,
        }
    }
}

fn matrix(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[i * n + j] != 0.0 {
                coo.push(i, j, values[i * n + j]);
            }
        }
    }
    coo.assemble()
}

fn scalar(omega: f64, coefficients: [f64; 4]) -> Fixture {
    Fixture::new(omega, coefficients.map(|v| matrix(&[v], 1)), vec![0])
}

fn close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "actual {actual:?}, expected {expected:?}, relative error {:e}",
        (actual - expected).abs() / expected.abs().max(1e-30)
    );
}

#[test]
fn high_q_resonator_keeps_complex_frequency_slope_at_the_peak() {
    let drive = C64::new(0.8, -0.3);
    let loss = 0.002;
    for omega in [19.8, 19.999, 20.0, 20.001, 20.2] {
        let f = scalar(omega, [400.0, 1.0, loss, 0.0]);
        let pencil = f.pencil();
        let primary = pencil.solve(&[drive], &mut || false).unwrap();
        let derivative = pencil
            .frequency_derivative(&primary.pressure, &[drive], &[C64::ZERO], &mut || false)
            .unwrap();
        // Independent scalar transfer function and quotient derivative.
        let a = C64::new(400.0 - omega * omega, omega * loss);
        let a_prime = C64::new(-2.0 * omega, loss);
        let numerator = C64::new(0.0, omega) * drive;
        let expected = numerator / a;
        let tangent = (C64::new(0.0, 1.0) * drive * a - numerator * a_prime) / (a * a);
        close(primary.pressure[0], expected, 1e-10);
        close(derivative.pressure[0], tangent, 2e-10);
        assert!(derivative.backward_error < 5e-13);
        if omega == 20.0 {
            close(primary.pressure[0], drive.scale(1.0 / loss), 1e-14);
            close(
                derivative.pressure[0],
                drive * C64::new(0.0, -2.0 / (loss * loss)),
                1e-14,
            );
        }
    }
}

#[test]
fn sealed_compliance_retains_its_low_frequency_phase_and_slope() {
    let drive = C64::new(0.002, -0.003);
    let slope = C64::new(-0.001, 0.004);
    let compliance = 0.8;
    for omega in [1e-8, 0.3, 30.0] {
        let f = scalar(omega, [0.0, compliance, 0.0, 0.0]);
        let pencil = f.pencil();
        let primary = pencil.solve(&[drive], &mut || false).unwrap();
        let derivative = pencil
            .frequency_derivative(&primary.pressure, &[drive], &[slope], &mut || false)
            .unwrap();
        close(
            primary.pressure[0],
            drive * C64::new(0.0, -1.0 / (omega * compliance)),
            2e-14,
        );
        let expected = drive * C64::new(0.0, 1.0 / (omega * omega * compliance))
            + slope * C64::new(0.0, -1.0 / (omega * compliance));
        close(derivative.pressure[0], expected, 2e-14);
    }
}

fn duct(cells: usize, omega: f64) -> Fixture {
    let (length, rho, c) = (1.2, 1.3, 2.0);
    let n = cells + 1;
    let dx = length / cells as f64;
    let mut k = Coo::new(n, n);
    let mut m = Coo::new(n, n);
    for cell in 0..cells {
        for i in 0..2 {
            for j in 0..2 {
                k.push(
                    cell + i,
                    cell + j,
                    if i == j { 1.0 } else { -1.0 } / (rho * dx),
                );
                m.push(
                    cell + i,
                    cell + j,
                    dx * if i == j { 2.0 } else { 1.0 } / (6.0 * rho * c * c),
                );
            }
        }
    }
    Fixture::new(
        omega,
        [
            k.assemble(),
            m.assemble(),
            Coo::new(n, n).assemble(),
            Coo::new(n, n).assemble(),
        ],
        (0..cells).collect(),
    ) // Exact pressure release at x=L.
}

#[test]
fn pressure_release_duct_tangent_converges_to_analytic_standing_wave() {
    let (length, rho, c, theta) = (1.2, 1.3, 2.0, 1.4);
    let omega = theta * c / length;
    let velocity = C64::new(0.02, -0.01);
    let exact = velocity * C64::new(0.0, rho * c * det::sin(theta) / det::cos(theta));
    let tangent = velocity * C64::new(0.0, rho * length / det::cos(theta).powi(2));
    let mut errors = Vec::new();
    for cells in [8, 16, 32] {
        let f = duct(cells, omega);
        let mut load = vec![C64::ZERO; cells + 1];
        load[0] = velocity;
        let pencil = f.pencil();
        let primary = pencil.solve(&load, &mut || false).unwrap();
        let derivative = pencil
            .frequency_derivative(
                &primary.pressure,
                &load,
                &vec![C64::ZERO; cells + 1],
                &mut || false,
            )
            .unwrap();
        assert_eq!(primary.pressure[cells], C64::ZERO);
        assert_eq!(derivative.pressure[cells], C64::ZERO);
        errors.push((
            (primary.pressure[0] - exact).abs() / exact.abs(),
            (derivative.pressure[0] - tangent).abs() / tangent.abs(),
        ));
    }
    for pair in errors.windows(2) {
        assert!(pair[1].0 < 0.3 * pair[0].0, "pressure errors {errors:?}");
        assert!(pair[1].1 < 0.3 * pair[0].1, "tangent errors {errors:?}");
    }
    assert!(errors[2].1 < 0.003, "tangent errors {errors:?}");
}

#[test]
fn coherent_delayed_sources_and_complex_boundaries_match_frequency_sweep() {
    let matrices = [
        matrix(&[5.0, -1.0, -1.0, 7.0], 2),
        matrix(&[1.0, 0.2, 0.2, 0.7], 2),
        matrix(&[0.3, -0.1, -0.1, 0.4], 2),
        matrix(&[0.2, 0.05, 0.05, -0.1], 2),
    ];
    let source = |w: f64| {
        let a = C64::new(det::cos(0.7 * w), -det::sin(0.7 * w));
        let b = C64::new(det::cos(1.3 * w), -det::sin(1.3 * w)).scale(0.8);
        let da = a * C64::new(0.0, -0.7);
        let db = b * C64::new(0.0, -1.3);
        ([a + b, a - b], [da + db, da - db])
    };
    let omega = 2.4;
    let f = Fixture::new(omega, matrices.clone(), vec![0, 1]);
    let (load, slopes) = source(omega);
    let primary = f.pencil().solve(&load, &mut || false).unwrap();
    let derivative = f
        .pencil()
        .frequency_derivative(&primary.pressure, &load, &slopes, &mut || false)
        .unwrap();
    let h = 1e-4;
    let fields: Vec<_> = [-2.0, -1.0, 1.0, 2.0]
        .iter()
        .map(|offset| {
            let w = omega + offset * h;
            let shifted = Fixture::new(w, matrices.clone(), vec![0, 1]);
            shifted
                .pencil()
                .solve(&source(w).0, &mut || false)
                .unwrap()
                .pressure
        })
        .collect();
    for i in 0..2 {
        let fd = (fields[0][i] - fields[1][i].scale(8.0) + fields[2][i].scale(8.0) - fields[3][i])
            .scale(1.0 / (12.0 * h));
        close(derivative.pressure[i], fd, 3e-9);
    }
}

#[test]
fn zero_primary_drive_keeps_its_nonzero_source_derivative() {
    let f = scalar(1.4, [3.0, 0.7, 0.2, -0.1]);
    let pencil = f.pencil();
    let slope = [C64::new(1.0, -2.0)];
    let primary = pencil.solve(&[C64::ZERO], &mut || false).unwrap();
    let derivative = pencil
        .frequency_derivative(&primary.pressure, &[C64::ZERO], &slope, &mut || false)
        .unwrap();
    let forced = pencil.solve(&slope, &mut || false).unwrap();
    assert_eq!(primary.pressure, vec![C64::ZERO]);
    close(derivative.pressure[0], forced.pressure[0], 2e-14);
}

#[test]
fn physical_equation_rescaling_preserves_pressure_and_frequency_tangent() {
    let mut results = Vec::new();
    for scale in [1e-120, 1.0, 1e120] {
        let f = scalar(3.0, [12.0, 1.0, 0.3, -0.2].map(|v| v * scale));
        let load = [C64::new(0.8, -0.3).scale(scale)];
        let slope = [C64::new(0.04, 0.07).scale(scale)];
        let primary = f.pencil().solve(&load, &mut || false).unwrap();
        let derivative = f
            .pencil()
            .frequency_derivative(&primary.pressure, &load, &slope, &mut || false)
            .unwrap();
        results.push((primary.pressure[0], derivative.pressure[0]));
    }
    for (p, dp) in &results {
        close(*p, results[1].0, 2e-14);
        close(*dp, results[1].1, 3e-14);
    }
}

#[test]
fn cancellation_and_invalid_derivatives_leave_the_retained_factor_reusable() {
    let f = scalar(2.0, [3.0, 1.0, 0.2, 0.1]);
    let pencil = f.pencil();
    let load = [C64::ONE];
    let primary = pencil.solve(&load, &mut || false).unwrap();
    let baseline = pencil
        .frequency_derivative(&primary.pressure, &load, &[C64::ZERO], &mut || false)
        .unwrap();
    for stop in [1, 4, 7] {
        let mut calls = 0;
        let result =
            pencil.frequency_derivative(&primary.pressure, &load, &[C64::ZERO], &mut || {
                calls += 1;
                calls >= stop
            });
        assert!(matches!(result, Err(PressureError::Cancelled)));
    }
    assert!(matches!(
        pencil.frequency_derivative(&primary.pressure, &load, &[], &mut || false),
        Err(PressureError::Invalid(_))
    ));
    assert!(matches!(
        pencil.frequency_derivative(
            &primary.pressure,
            &load,
            &[C64::new(f64::NAN, 0.0)],
            &mut || false
        ),
        Err(PressureError::Invalid(_))
    ));
    let replay = pencil
        .frequency_derivative(&primary.pressure, &load, &[C64::ZERO], &mut || false)
        .unwrap();
    assert_eq!(baseline.pressure, replay.pressure);
    assert_eq!(
        baseline.backward_error.to_bits(),
        replay.backward_error.to_bits()
    );
}
