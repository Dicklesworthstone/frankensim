use super::*;
use fs_phs::linear_exact::exponential;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn advance(matrix: &[f64], x: &[f64]) -> Vec<f64> {
    matrix
        .chunks_exact(x.len())
        .map(|row| dot(row, x))
        .collect()
}
fn step_matrix(generator: &[f64], inlet: &[f64], held: f64, dt: f64) -> Vec<f64> {
    let n = inlet.len();
    let d = n + 1;
    let mut augmented = vec![0.0; d * d];
    for i in 0..n {
        for j in 0..n {
            augmented[i * d + j] = dt * generator[i * n + j];
        }
        augmented[i * d + n] = dt * inlet[i] * held;
    }
    exponential(&augmented, d).unwrap()
}

#[test]
fn pressure_release_retains_cavity_resonance_and_the_original_outward_flow() {
    // Two compliant volumes joined by an inertance, in energy coordinates.
    // Releasing the second pressure leaves the FIRST volume's exact oscillator.
    let (omega, outlet, overlap, held) = (730.0_f64, 510.0, 0.4, 0.2);
    let source = [0.0, 0.0, -omega, 0.0, 0.0, outlet, omega, -outlet, 0.0];
    let retained = eliminate(
        &source,
        &[vec![0.0, -1.0, 0.0]],
        &[vec![0.0, 1.0, 0.0]],
        &[1.0, overlap, 0.0],
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(retained.generator.len(), 9); // No retained state is discarded.
    close(retained.inlet_flow_factors[0], overlap, 1e-14);
    let duration = 0.011;
    for panels in [11, 22, 44] {
        let dt = duration / panels as f64;
        let transition = step_matrix(&retained.generator, &retained.inlet_column, held, dt);
        let mut x = vec![1.0, 0.0, 0.0, 1.0];
        for step in 1..=panels {
            x = advance(&transition, &x);
            let time = dt * step as f64;
            let (sin, cos) = (omega * time).sin_cos();
            let pressure = cos + held / omega * sin;
            let velocity = sin + held / omega * (1.0 - cos);
            close(x[0], pressure, 3e-12);
            close(x[1], 0.0, 1e-14);
            close(x[2], velocity, 3e-12);
            let outward =
                dot(&retained.flow_columns[0], &x[..3]) + retained.inlet_flow_factors[0] * held;
            close(outward, outlet * velocity + overlap * held, 3e-12);
            let authored_work = held * (sin / omega + held / (omega * omega) * (1.0 - cos));
            close(
                0.5 * (x[0] * x[0] + x[2] * x[2]) - 0.5,
                authored_work,
                3e-12,
            );
        }
    }
}

#[test]
fn strictly_proper_relaxation_preserves_loaded_frequency_decay_and_passive_work() {
    // Original acoustic oscillator plus Z(s)=k/(s+a). The authored terminal
    // has E=k*z^2/2 and heat=a*k*z^2; no generic fitted storage is assumed.
    let (omega, k, a) = (700.0_f64, 3.0_f64, 80.0_f64);
    let retained = eliminate(
        &[0.0, -omega, 0.0, omega, 0.0, 0.0, 0.0, 0.0, -a],
        &[vec![-1.0, 0.0, 1.0]],
        &[vec![1.0, 0.0, -k]],
        &[1.0, 0.0, 0.0],
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let decay = a / (2.0 * (1.0 + k));
    let natural = omega * (k / (1.0 + k)).sqrt();
    let beta = (natural * natural - decay * decay).sqrt();
    let duration = 5.0 * core::f64::consts::TAU / beta;
    let initial = vec![1.0, 0.0, 1.0 / k, 1.0];
    let initial_energy = 0.5 * (1.0 + 1.0 / k);
    for panels in [80, 160, 320] {
        let dt = duration / panels as f64;
        let transition = step_matrix(&retained.generator, &retained.inlet_column, 0.0, dt);
        let mut x = initial.clone();
        for step in 1..=panels {
            x = advance(&transition, &x);
            let time = step as f64 * dt;
            let (sin, cos) = (beta * time).sin_cos();
            let envelope = (-decay * time).exp();
            let pressure = envelope * (cos - decay / beta * sin);
            let velocity = envelope * omega / beta * sin;
            close(x[0], pressure, 2e-11);
            close(x[1], velocity, 2e-11);
            close(x[2], pressure / k, 2e-11);
            assert!(0.5 * (x[0] * x[0] + x[1] * x[1] + k * x[2] * x[2]) <= initial_energy + 1e-11);
        }
        close(x[0], (-decay * duration).exp(), 2e-11);
    }
    // Independent Simpson integration of ORIGINAL pressure*flow and physical
    // terminal heat, along the actual retained exact trajectory.
    let panels = 4096;
    let dt = duration / panels as f64;
    let transition = step_matrix(&retained.generator, &retained.inlet_column, 0.0, dt);
    let mut x = initial.clone();
    let mut work = 0.0;
    let mut heat = 0.0;
    for step in 0..=panels {
        let weight = if step == 0 || step == panels {
            1.0
        } else if step % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let outward = dot(&retained.flow_columns[0], &x[..3]);
        work += weight * x[0] * outward;
        heat += weight * a * k * x[2] * x[2];
        if step != panels {
            x = advance(&transition, &x);
        }
    }
    work *= dt / 3.0;
    heat *= dt / 3.0;
    let terminal_change = 0.5 * k * (x[2] * x[2] - initial[2] * initial[2]);
    let final_energy = 0.5 * (x[0] * x[0] + x[1] * x[1] + k * x[2] * x[2]);
    close(work, terminal_change + heat, 2e-9);
    close(final_energy + heat, initial_energy, 2e-9);
    assert!(heat > 0.0);
    eprintln!(
        "loaded omega={natural:.9}, decay={decay:.9}, Q={:.9}",
        natural / (2.0 * decay)
    );
}

#[test]
fn mutual_strictly_proper_terminal_matches_independent_physical_reduction() {
    // Reciprocal passive K/(s+a), including a signed off-diagonal residue.
    // Both original terminal states stay present. The reference instead uses
    // (I+K) z' = -w*v-a*z+b*u and v'=w^T K z, derived physically.
    let (k11, k12, k22, a, w1, w2, held) = (
        3.0_f64, -0.7_f64, 2.0_f64, 40.0_f64, 310.0_f64, -190.0_f64, 0.15_f64,
    );
    let b = [1.0, 0.3];
    let source = [
        0.0, 0.0, -w1, 0.0, 0.0, 0.0, 0.0, -w2, 0.0, 0.0, w1, w2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -a,
        0.0, 0.0, 0.0, 0.0, 0.0, -a,
    ];
    let inputs = [
        vec![-1.0, 0.0, 0.0, 1.0, 0.0],
        vec![0.0, -1.0, 0.0, 0.0, 1.0],
    ];
    let h = [
        vec![1.0, 0.0, 0.0, -k11, -k12],
        vec![0.0, 1.0, 0.0, -k12, -k22],
    ];
    let inlet = [b[0], b[1], 0.0, 0.0, 0.0];
    let retained = eliminate(&source, &inputs, &h, &inlet, &CancelGate::new_clock_free()).unwrap();
    let denominator = (1.0 + k11) * (1.0 + k22) - k12 * k12;
    let r = [
        (1.0 + k22) / denominator,
        -k12 / denominator,
        -k12 / denominator,
        (1.0 + k11) / denominator,
    ];
    let reference = [
        -a * r[0],
        -a * r[1],
        -r[0] * w1 - r[1] * w2,
        -a * r[2],
        -a * r[3],
        -r[2] * w1 - r[3] * w2,
        w1 * k11 + w2 * k12,
        w1 * k12 + w2 * k22,
        0.0,
    ];
    let reference_inlet = [r[0] * b[0] + r[1] * b[1], r[2] * b[0] + r[3] * b[1], 0.0];
    for (actual, expected) in retained.inlet_flow_factors.iter().zip(reference_inlet) {
        close(*actual, expected, 2e-14);
    }
    let dt = 0.0001;
    let transition = step_matrix(&retained.generator, &retained.inlet_column, held, dt);
    let physical = step_matrix(&reference, &reference_inlet, held, dt);
    let (z1, z2, velocity) = (0.1, -0.2, 0.05);
    let mut x = vec![
        k11 * z1 + k12 * z2,
        k12 * z1 + k22 * z2,
        velocity,
        z1,
        z2,
        1.0,
    ];
    let mut reference_state = vec![z1, z2, velocity, 1.0];
    for _ in 0..400 {
        x = advance(&transition, &x);
        reference_state = advance(&physical, &reference_state);
        close(
            x[0],
            k11 * reference_state[0] + k12 * reference_state[1],
            3e-11,
        );
        close(
            x[1],
            k12 * reference_state[0] + k22 * reference_state[1],
            3e-11,
        );
        close(x[2], reference_state[2], 3e-11);
        close(x[3], reference_state[0], 3e-11);
        close(x[4], reference_state[1], 3e-11);
        for row in &h {
            close(dot(row, &x[..5]), 0.0, 3e-11);
        }
    }
    // Re-expressing pressure constraints in different units cannot change
    // the physical generator or its source-flow correction.
    let scaled = [
        h[0].iter().map(|v| 1e9 * v).collect(),
        h[1].iter().map(|v| 1e-7 * v).collect(),
    ];
    let reexpressed = eliminate(
        &source,
        &inputs,
        &scaled,
        &inlet,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for (actual, expected) in retained.generator.iter().zip(&reexpressed.generator) {
        close(*actual, *expected, 1e-13);
    }
}

#[test]
fn singular_hidden_constraints_and_cancelled_preparation_refuse() {
    let gate = CancelGate::new_clock_free();
    // The pressure multiplier cannot affect this pressure derivative:
    // a second differentiation / a higher-index DAE would be necessary.
    assert!(
        eliminate(
            &[0.0, 1.0, -1.0, 0.0],
            &[vec![0.0, 1.0]],
            &[vec![1.0, 0.0]],
            &[1.0, 0.0],
            &gate
        )
        .is_err()
    );
    // The scalar passive constraint must oppose outward flow.
    assert!(eliminate(&[0.0], &[vec![1.0]], &[vec![1.0]], &[0.0], &gate).is_err());
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(matches!(
        eliminate(&[0.0], &[vec![-1.0]], &[vec![1.0]], &[0.0], &cancel),
        Err(AcousticRealizeError::Cancelled)
    ));
}
