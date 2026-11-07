use super::super::{det, propagation};
use super::*;

fn near(a: C64, b: C64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-280), "{a:?} != {b:?}");
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-280), "{a:e} != {b:e}");
}
fn sample(length: f64, rho: f64, bulk: C64, shear: C64, exterior: [f64; 2]) -> Layer {
    Layer {
        wave: Slab::new(rho, [bulk, shear], 50.0, length, exterior).unwrap(),
        maximum_strain: 10.0,
    }
}
fn heterogeneous(reference: [f64; 2]) -> Vec<Layer> {
    vec![
        sample(
            0.3,
            2.0,
            C64::new(80.0, 20.0),
            C64::new(30.0, 8.0),
            [20.0, reference[0]],
        ),
        sample(
            0.17,
            7.0,
            C64::new(170.0, 3.0),
            C64::new(45.0, 27.0),
            reference,
        ),
        sample(
            0.41,
            3.0,
            C64::new(95.0, 7.0),
            C64::new(20.0, 1.0),
            [reference[1], 35.0],
        ),
    ]
}

// Independent dense system in the PHYSICAL inward waves of each solid.
// Its equations contain no internal reference impedances or scattering/star
// products: the 2(n+1) rows directly impose traction and velocity continuity.
fn boundary_solution(layers: &[Layer], input: [C64; 2]) -> Vec<C64> {
    let n = layers.len();
    let m = 2 * n + 2;
    let mut matrix = vec![vec![C64::ZERO; m + 1]; m];
    let z = layers[0].wave.impedance.recip();
    let g = propagation(layers[0].wave.k, layers[0].wave.length).unwrap();
    let ext = 1.0 / layers[0].wave.exterior[0];
    matrix[0][0] = -C64::ONE;
    matrix[0][2] = C64::ONE;
    matrix[0][3] = g;
    matrix[0][m] = input[0];
    matrix[1][0] = C64::from_re(ext);
    matrix[1][2] = z;
    matrix[1][3] = -(g * z);
    matrix[1][m] = input[0].scale(ext);
    for i in 0..n - 1 {
        let row = 2 + 2 * i;
        let a = &layers[i].wave;
        let b = &layers[i + 1].wave;
        let ga = propagation(a.k, a.length).unwrap();
        let gb = propagation(b.k, b.length).unwrap();
        let p = 2 + 2 * i;
        matrix[row][p] = ga;
        matrix[row][p + 1] = C64::ONE;
        matrix[row][p + 2] = -C64::ONE;
        matrix[row][p + 3] = -gb;
        matrix[row + 1][p] = ga / a.impedance;
        matrix[row + 1][p + 1] = -a.impedance.recip();
        matrix[row + 1][p + 2] = -b.impedance.recip();
        matrix[row + 1][p + 3] = gb / b.impedance;
    }
    let slab = &layers[n - 1].wave;
    let g = propagation(slab.k, slab.length).unwrap();
    let z = slab.impedance.recip();
    let ext = 1.0 / slab.exterior[1];
    let row = 2 * n;
    matrix[row][1] = -C64::ONE;
    matrix[row][row] = g;
    matrix[row][row + 1] = C64::ONE;
    matrix[row][m] = input[1];
    matrix[row + 1][1] = C64::from_re(-ext);
    matrix[row + 1][row] = g * z;
    matrix[row + 1][row + 1] = -z;
    matrix[row + 1][m] = input[1].scale(-ext);
    for col in 0..m {
        let pivot = (col..m)
            .max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs()))
            .unwrap();
        matrix.swap(col, pivot);
        let diagonal = matrix[col][col];
        for j in col..=m {
            matrix[col][j] = matrix[col][j] / diagonal;
        }
        for i in 0..m {
            if i == col {
                continue;
            }
            let multiplier = matrix[i][col];
            for j in col..=m {
                matrix[i][j] = matrix[i][j] - multiplier * matrix[col][j];
            }
        }
    }
    matrix.iter().map(|row| row[m]).collect()
}

#[test]
fn g1_heterogeneous_boundary_equations_and_constitutive_layer_heat() {
    let original = heterogeneous([28.0, 52.0]);
    let input = [C64::new(1.1, 0.2), C64::new(-0.4, 0.7)];
    let solution = boundary_solution(&original, input);
    let stack = Stack::new(original.clone(), &mut || false).unwrap();
    let result = stack.drive(input, 0.3, &mut || false).unwrap();
    near(result.outgoing_pressure[0], solution[0], 3e-13);
    near(result.outgoing_pressure[1], solution[1], 3e-13);
    for (i, layer) in original.iter().enumerate() {
        let slab = layer.wave;
        let local = stack.layer_drive(i, input, 0.3).unwrap();
        let (a, b) = (solution[2 + 2 * i], solution[3 + 2 * i]);
        let n = 8192;
        let mut integral = 0.0;
        for j in 0..=n {
            let x = slab.length * j as f64 / n as f64;
            let p = a * propagation(slab.k, x).unwrap()
                + b * propagation(slab.k, slab.length - x).unwrap();
            let strain = p / slab.modulus;
            let weight = if j == 0 || j == n {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * strain.norm_sq() * slab.length / (3.0 * n as f64);
        }
        close(
            local.heat[0],
            0.3 * 0.5 * slab.omega * slab.moduli[0].im * integral,
            2e-12,
        );
        close(
            local.heat[1],
            0.3 * (2.0 / 3.0) * slab.omega * slab.moduli[1].im * integral,
            2e-12,
        );
        assert!(local.heat.into_iter().all(|h| h > 0.0));
        if i + 1 < original.len() {
            let left = stack.field(i, input, slab.length).unwrap();
            let right = stack.field(i + 1, input, 0.0).unwrap();
            near(left.pressure, right.pressure, 3e-13);
            near(left.velocity, right.velocity, 3e-13);
            assert!((left.strain - right.strain).abs() > 1e-5);
        }
    }
    close(
        result.incident.iter().sum::<f64>(),
        result.outgoing.iter().sum::<f64>() + result.heat.iter().sum::<f64>(),
        3e-13,
    );
}

#[test]
fn g3_reference_coordinates_and_physical_subdivision_leave_solution_unchanged() {
    let original = heterogeneous([28.0, 52.0]);
    let stack = Stack::new(original.clone(), &mut || false).unwrap();
    let changed = Stack::new(heterogeneous([3.0, 210.0]), &mut || false).unwrap();
    let split = original[1].wave;
    let mut subdivided = vec![original[0]];
    for (length, exterior) in [(0.07, [28.0, 60.0]), (0.10, [60.0, 52.0])] {
        subdivided.push(sample(
            length,
            7.0,
            split.moduli[0],
            split.moduli[1],
            exterior,
        ));
    }
    subdivided.push(original[2]);
    let subdivided = Stack::new(subdivided, &mut || false).unwrap();
    let input = [C64::new(0.8, -0.3), C64::new(0.1, 0.6)];
    let original_power = stack.drive(input, 0.7, &mut || false).unwrap();
    for other in [&changed, &subdivided] {
        let power = other.drive(input, 0.7, &mut || false).unwrap();
        for i in 0..2 {
            near(
                power.outgoing_pressure[i],
                original_power.outgoing_pressure[i],
                3e-13,
            );
            close(power.heat[i], original_power.heat[i], 3e-13);
        }
    }
    let heat = stack.layer_drive(1, input, 0.7).unwrap().heat;
    let a = subdivided.layer_drive(1, input, 0.7).unwrap().heat;
    let b = subdivided.layer_drive(2, input, 0.7).unwrap().heat;
    for i in 0..2 {
        close(heat[i], a[i] + b[i], 3e-13);
    }
    let before = stack.field(1, input, 0.12).unwrap();
    let after = subdivided.field(2, input, 0.05).unwrap();
    near(before.pressure, after.pressure, 3e-13);
    near(before.velocity, after.velocity, 3e-13);
}

#[test]
fn g1_single_layer_and_zero_thickness_limits_preserve_original_owner() {
    let layer = sample(
        0.4,
        2.0,
        C64::new(80.0, 20.0),
        C64::new(30.0, 8.0),
        [20.0, 35.0],
    );
    let input = [C64::new(0.3, -0.9), C64::new(0.1, 0.2)];
    let original = layer.wave.drive(input, 0.7).unwrap();
    let stack = Stack::new(vec![layer], &mut || false).unwrap();
    let result = stack.drive(input, 0.7, &mut || false).unwrap();
    assert_eq!(stack.layer(0).unwrap().moduli, layer.wave.moduli);
    assert_eq!(stack.scattering, layer.wave.scattering);
    for i in 0..2 {
        near(result.outgoing_pressure[i], original.outgoing[i], 3e-15);
        close(result.heat[i], original.heat[i], 3e-15);
    }
    let zero = vec![
        sample(
            0.0,
            7.0,
            C64::new(75.0, 8.0),
            C64::new(31.0, 7.0),
            [20.0, 49.0],
        ),
        sample(
            0.0,
            3.0,
            C64::new(14.0, 2.0),
            C64::new(6.0, 1.0),
            [49.0, 35.0],
        ),
    ];
    let stack = Stack::new(zero, &mut || false).unwrap();
    let result = stack
        .drive([C64::ONE, C64::ZERO], 1.0, &mut || false)
        .unwrap();
    near(
        result.outgoing_pressure[0],
        C64::from_re(15.0 / 55.0),
        3e-15,
    );
    near(
        result.outgoing_pressure[1],
        C64::from_re(70.0 / 55.0),
        3e-15,
    );
    assert_eq!(result.heat, [0.0; 2]);
    assert!(stack.field(0, input, 0.0).is_err());
}

#[test]
fn g1_opaque_shielding_and_weak_thin_standing_wave_work_are_positive() {
    let mut layers = heterogeneous([28.0, 52.0]);
    layers[0] = sample(
        1e6,
        2.0,
        C64::new(80.0, 20.0),
        C64::new(30.0, 8.0),
        [20.0, 28.0],
    );
    let stack = Stack::new(layers, &mut || false).unwrap();
    let input = [C64::ONE, C64::ZERO];
    let power = stack.drive(input, 1.0, &mut || false).unwrap();
    assert_eq!(power.outgoing[1], 0.0);
    assert_eq!(stack.layer_drive(1, input, 1.0).unwrap().heat, [0.0; 2]);
    assert_eq!(stack.layer_drive(2, input, 1.0).unwrap().heat, [0.0; 2]);
    close(
        power.heat.iter().sum::<f64>(),
        power.incident[0] - power.outgoing[0],
        2e-14,
    );
    for length in [1e-8, 1e-110] {
        let input = [C64::from_re(1e100), C64::from_re(-1e100)];
        let mut whole = sample(
            length,
            4.0,
            C64::new(60.0, 1e-24),
            C64::from_re(30.0),
            [20.0; 2],
        );
        whole.maximum_strain = 1e110;
        let mut half = sample(
            length / 2.0,
            4.0,
            C64::new(60.0, 1e-24),
            C64::from_re(30.0),
            [20.0; 2],
        );
        half.maximum_strain = 1e110;
        let direct = whole.wave.drive(input, 1.0).unwrap();
        let split = Stack::new(vec![half; 2], &mut || false).unwrap();
        let result = split.drive(input, 1.0, &mut || false).unwrap();
        assert!(result.heat[0] > 0.0);
        assert_eq!(result.heat[1], 0.0);
        close(result.heat[0], direct.heat[0], 5e-12);
        if length < 1e-100 {
            // Here 1-cos(k*d) is below representability, while the first-order
            // imaginary field and its dimensional heat remain representable.
            let original_field = whole
                .wave
                .field(direct.inward, direct.parity, length / 4.0)
                .unwrap();
            let subdivided_field = split.field(0, input, length / 4.0).unwrap();
            near(subdivided_field.pressure, original_field.pressure, 5e-13);
            near(subdivided_field.velocity, original_field.velocity, 5e-13);
            close(
                subdivided_field.heat_density[0],
                original_field.heat_density[0],
                5e-12,
            );
        }
    }
}

#[test]
fn g3_coherent_phase_area_and_reciprocity() {
    let stack = Stack::new(heterogeneous([28.0, 52.0]), &mut || false).unwrap();
    near(
        stack.scattering[1].scale(35.0),
        stack.scattering[2].scale(20.0),
        3e-14,
    );
    let input = [C64::new(0.8, 0.1), C64::new(0.3, -0.2)];
    let first = stack.drive(input, 0.4, &mut || false).unwrap();
    let (s, c) = det::sin_cos(0.71);
    let phase = C64::new(c, s);
    let rotated = stack
        .drive(input.map(|p| p * phase), 0.4, &mut || false)
        .unwrap();
    let larger = stack.drive(input, 2.8, &mut || false).unwrap();
    let opposed = stack
        .drive([input[0], -input[1]], 0.4, &mut || false)
        .unwrap();
    for i in 0..2 {
        near(
            rotated.outgoing_pressure[i],
            first.outgoing_pressure[i] * phase,
            3e-14,
        );
        close(rotated.heat[i], first.heat[i], 3e-14);
        close(larger.heat[i], 7.0 * first.heat[i], 3e-14);
    }
    assert!((opposed.heat[0] - first.heat[0]).abs() > 1e-5);
}

#[test]
fn g0_g4_support_input_join_and_cancellation_refusals_publish_nothing() {
    let original = heterogeneous([28.0, 52.0]);
    assert!(Stack::new(vec![], &mut || false).is_err());
    let mut bad = original.clone();
    bad[1].wave.exterior[0] = 4.0;
    assert!(Stack::new(bad, &mut || false).is_err());
    let mut bad = original.clone();
    bad[1].wave.omega = 51.0;
    assert!(Stack::new(bad, &mut || false).is_err());
    let mut weak = original.clone();
    weak[1].maximum_strain = 1e-10;
    let weak = Stack::new(weak, &mut || false).unwrap();
    assert_eq!(
        weak.drive([C64::ONE, C64::ZERO], 1.0, &mut || false),
        Err(Error::OutsideSourceSupport)
    );
    for after in [1, 4, 7, 10] {
        let mut calls = 0;
        assert!(matches!(
            Stack::new(original.clone(), &mut || {
                calls += 1;
                calls == after
            }),
            Err(Error::Cancelled)
        ));
    }
    let stack = Stack::new(original, &mut || false).unwrap();
    for after in [1, 2, 4, 5] {
        let mut calls = 0;
        assert_eq!(
            stack.drive([C64::ONE, C64::ZERO], 1.0, &mut || {
                calls += 1;
                calls == after
            }),
            Err(Error::Cancelled)
        );
    }
    for area in [-1.0, 0.0, f64::INFINITY, f64::NAN] {
        assert!(
            stack
                .drive([C64::ONE, C64::ZERO], area, &mut || false)
                .is_err()
        );
    }
    assert!(
        stack
            .drive([C64::new(f64::NAN, 0.0), C64::ZERO], 1.0, &mut || false)
            .is_err()
    );
    assert!(stack.layer_drive(3, [C64::ONE, C64::ZERO], 1.0).is_err());
    assert!(
        stack
            .drive([C64::ONE, C64::ZERO], 1.0, &mut || false)
            .is_ok()
    );
}
