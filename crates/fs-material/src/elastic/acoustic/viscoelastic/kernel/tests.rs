use super::*;

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "actual={actual:?} expected={expected:?}"
    );
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:e} expected={expected:e}"
    );
}

// Independent four-unknown displacement/traction continuity system. This does
// not use the production two-wave determinant or its scattering coefficients.
fn boundary_solution(slab: &Slab, input: [C64; 2]) -> [C64; 4] {
    let one = C64::ONE;
    let zero = C64::ZERO;
    let [l, r] = slab.exterior.map(|z| C64::from_re(z).recip());
    let z = slab.impedance.recip();
    let g = propagation(slab.k, slab.length).unwrap();
    let mut a = [
        [-one, zero, one, g, input[0]],
        [l, zero, z, -(g * z), input[0] * l],
        [zero, -one, g, one, input[1]],
        [zero, -r, g * z, -z, -(input[1] * r)],
    ];
    for col in 0..4 {
        let pivot = (col..4)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap();
        a.swap(col, pivot);
        let p = a[col][col];
        for j in col..5 {
            a[col][j] = a[col][j] / p;
        }
        for row in 0..4 {
            if row == col {
                continue;
            }
            let scale = a[row][col];
            for j in col..5 {
                a[row][j] = a[row][j] - scale * a[col][j];
            }
        }
    }
    core::array::from_fn(|i| a[i][4])
}

#[test]
fn g1_independent_boundary_equations_and_spatial_constitutive_heat() {
    let moduli = [C64::new(80.0, 20.0), C64::new(30.0, 8.0)];
    let input = [C64::new(1.1, 0.2), C64::new(-0.4, 0.7)];
    for thickness in [0.01, 0.3, 1.2] {
        let slab = Slab::new(2.0, moduli, 50.0, thickness, [20.0, 35.0]).unwrap();
        let driven = slab.drive(input, 0.3).unwrap();
        assert!(driven.residual.abs() < 1e-13 * driven.incident_power.into_iter().sum::<f64>());
        let [left, right, a, b] = boundary_solution(&slab, input);
        near(driven.outgoing[0], left, 2e-13);
        near(driven.outgoing[1], right, 2e-13);
        let front = slab.field(driven.inward, driven.parity, 0.0).unwrap();
        let back = slab.field(driven.inward, driven.parity, thickness).unwrap();
        near(front.pressure, input[0] + left, 2e-13);
        near(front.velocity, (input[0] - left).scale(1.0 / 20.0), 2e-13);
        near(back.pressure, input[1] + right, 2e-13);
        near(back.velocity, (right - input[1]).scale(1.0 / 35.0), 2e-13);
        let n = 16_384;
        let mut strain_integral = 0.0;
        for i in 0..n {
            let x = thickness * (i as f64 + 0.5) / n as f64;
            let pressure = a * propagation(slab.k, x).unwrap()
                + b * propagation(slab.k, thickness - x).unwrap();
            let strain = -pressure / (moduli[0] + moduli[1].scale(4.0 / 3.0));
            strain_integral += strain.norm_sq() * thickness / n as f64;
        }
        close(
            driven.heat[0],
            0.3 * 0.5 * 50.0 * moduli[0].im * strain_integral,
            3e-8,
        );
        close(
            driven.heat[1],
            0.3 * (2.0 / 3.0) * 50.0 * moduli[1].im * strain_integral,
            3e-8,
        );
        let field = slab
            .field(driven.inward, driven.parity, 0.37 * thickness)
            .unwrap();
        close(
            field.heat_density[0],
            0.5 * 50.0 * moduli[0].im * field.strain.norm_sq(),
            2e-13,
        );
        close(
            field.heat_density[1],
            (2.0 / 3.0) * 50.0 * moduli[1].im * field.strain.norm_sq(),
            2e-13,
        );
    }
}

#[test]
fn g1_lossless_matched_delay_and_zero_thickness_interface() {
    let slab = Slab::new(
        4.0,
        [C64::from_re(60.0), C64::from_re(30.0)],
        20.0,
        0.4,
        [20.0; 2],
    )
    .unwrap();
    assert_eq!(slab.scattering[0], C64::ZERO);
    near(
        slab.scattering[2],
        propagation(C64::from_re(4.0), 0.4).unwrap(),
        2e-15,
    );
    let driven = slab.drive([C64::ONE, C64::ZERO], 1.0).unwrap();
    assert_eq!(driven.heat, [0.0; 2]);
    close(driven.outgoing_power[1], driven.incident_power[0], 2e-14);
    let identity = Slab::new(
        4.0,
        [C64::new(60.0, 50.0), C64::new(30.0, 20.0)],
        20.0,
        0.0,
        [20.0, 35.0],
    )
    .unwrap();
    let driven = identity.drive([C64::ONE, C64::ZERO], 1.0).unwrap();
    near(driven.outgoing[0], C64::from_re(15.0 / 55.0), 2e-15);
    near(driven.outgoing[1], C64::from_re(70.0 / 55.0), 2e-15);
    assert_eq!(driven.heat, [0.0; 2]);
    assert_eq!(driven.strain_bound, 0.0);
    assert!(identity.field(driven.inward, driven.parity, 0.0).is_err());
}

#[test]
fn g1_positive_weak_loss_and_antisymmetric_thin_heat_below_roundoff() {
    let slab = Slab::new(
        4.0,
        [C64::new(60.0, 1e-24), C64::from_re(30.0)],
        20.0,
        0.8,
        [20.0; 2],
    )
    .unwrap();
    let driven = slab.drive([C64::ONE, C64::ZERO], 1.0).unwrap();
    assert!(driven.heat[0] > 0.0);
    close(driven.heat[0], 0.5 * 20.0 * 1e-24 * 0.8 / 10000.0, 3e-14);
    assert!(driven.heat[0] < f64::EPSILON * driven.incident_power[0]);
    for (length, pressure) in [(1e-8, 1.0), (1e-170, 1e150)] {
        let slab = Slab::new(
            4.0,
            [C64::new(60.0, 1e-12), C64::from_re(30.0)],
            20.0,
            length,
            [20.0; 2],
        )
        .unwrap();
        let driven = slab
            .drive([C64::from_re(pressure), C64::from_re(-pressure)], 1.0)
            .unwrap();
        // Independent thin odd-wave limit: integral |p|²=P² k² d³/3.
        let p_length = pressure * length;
        let expected = ((p_length * p_length) * length) * (20.0 * 1e-12 * 16.0 / 60000.0);
        assert!(driven.heat[0] > 0.0);
        close(driven.heat[0], expected, 3e-12);
    }
}

#[test]
fn g1_opaque_slab_uses_half_space_heat_without_growing_waves() {
    let slab = Slab::new(
        2.0,
        [C64::new(80.0, 20.0), C64::new(30.0, 8.0)],
        1000.0,
        1000.0,
        [20.0, 35.0],
    )
    .unwrap();
    let driven = slab.drive([C64::ONE, C64::ZERO], 0.7).unwrap();
    assert_eq!(driven.outgoing[1], C64::ZERO);
    let reference = (slab.impedance - C64::from_re(20.0)) / (slab.impedance + C64::from_re(20.0));
    near(driven.outgoing[0], reference, 2e-15);
    close(
        driven.heat[0] + driven.heat[1],
        0.7 * (1.0 - reference.norm_sq()) / 40.0,
        2e-13,
    );
    close(
        driven.heat[0] / driven.heat[1],
        20.0 / ((4.0 / 3.0) * 8.0),
        2e-13,
    );
    assert_eq!(
        slab.field(driven.inward, driven.parity, 1000.0)
            .unwrap()
            .heat_density,
        [0.0; 2]
    );
}

#[test]
fn g3_coherent_subdivision_and_drive_area_scaling() {
    let moduli = [C64::new(80.0, 20.0), C64::new(30.0, 8.0)];
    let whole = Slab::new(2.0, moduli, 50.0, 0.8, [20.0; 2]).unwrap();
    let half = Slab::new(2.0, moduli, 50.0, 0.4, [20.0; 2]).unwrap();
    let input = [C64::new(1.0, 0.3), C64::new(-0.2, 0.1)];
    let driven = whole.drive(input, 1.0).unwrap();
    let s = half.scattering;
    let middle_right = (s[2] * input[0] + s[3] * s[1] * input[1]) / (C64::ONE - s[3] * s[0]);
    let middle_left = s[0] * middle_right + s[1] * input[1];
    let left = half.drive([input[0], middle_left], 1.0).unwrap();
    let right = half.drive([middle_right, input[1]], 1.0).unwrap();
    near(driven.outgoing[0], left.outgoing[0], 2e-13);
    near(driven.outgoing[1], right.outgoing[1], 2e-13);
    for channel in 0..2 {
        close(
            driven.heat[channel],
            left.heat[channel] + right.heat[channel],
            2e-13,
        );
        close(
            whole.drive(input, 4.0).unwrap().heat[channel],
            4.0 * driven.heat[channel],
            2e-13,
        );
        close(
            whole.drive(input.map(|p| p.scale(2.0)), 1.0).unwrap().heat[channel],
            4.0 * driven.heat[channel],
            2e-13,
        );
    }
    let symmetric = whole.drive([C64::ONE, C64::ONE], 1.0).unwrap();
    let antisymmetric = whole.drive([C64::ONE, -C64::ONE], 1.0).unwrap();
    assert_eq!(symmetric.incident_power, antisymmetric.incident_power);
    assert!((symmetric.heat[0] / antisymmetric.heat[0] - 1.0).abs() > 0.1);
    assert_eq!(whole.drive([C64::ZERO; 2], 1.0).unwrap().heat, [0.0; 2]);
}

#[test]
fn g0_invalid_physical_and_observation_inputs_refuse() {
    let m = [C64::new(80.0, 20.0), C64::new(30.0, 8.0)];
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(Slab::new(invalid, m, 50.0, 0.8, [20.0; 2]).is_err());
        assert!(Slab::new(2.0, m, invalid, 0.8, [20.0; 2]).is_err());
        assert!(Slab::new(2.0, m, 50.0, 0.8, [invalid, 20.0]).is_err());
    }
    assert!(Slab::new(2.0, [C64::new(80.0, -1.0), m[1]], 50.0, 0.8, [20.0; 2]).is_err());
    let slab = Slab::new(2.0, m, 50.0, 0.8, [20.0; 2]).unwrap();
    assert!(slab.drive([C64::ONE; 2], 0.0).is_err());
    assert!(
        slab.drive([C64::new(f64::NAN, 0.0), C64::ZERO], 1.0)
            .is_err()
    );
    let drive = slab.drive([C64::ONE, C64::ZERO], 1.0).unwrap();
    for depth in [-1.0, 0.9, f64::NAN] {
        assert!(slab.field(drive.inward, drive.parity, depth).is_err());
    }
}
