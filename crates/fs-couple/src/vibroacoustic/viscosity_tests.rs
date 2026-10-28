use super::*;
use fs_exec::CancelGate;

fn viscosity(matrix: Vec<f64>, modes: usize, resistance: f64) -> ModalPlateViscosity {
    ModalPlateViscosity::from_psd_modal_matrix(
        matrix, modes, resistance, &CancelGate::new_clock_free(),
    ).expect("passive material or relative dashpot")
}

fn model(
    omegas: Vec<f64>, coupling: Vec<f64>, eta_s: f64, eta_a: f64,
    radiation: Option<Vec<C64>>,
) -> VibroacousticModel {
    let n = omegas.len();
    VibroacousticModel::try_new(
        &StructuralModes { omegas, shapes: vec![vec![1.0]; n], loss_factor: eta_s },
        &CavityModes {
            omegas: vec![4.4], lambdas: vec![0.7], interface: vec![vec![1.0]],
            loss_factor: eta_a, rho0: 1.2, c0: 3.0,
        },
        coupling, radiation,
    ).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-14),
        "actual={actual:.17e}, expected={expected:.17e}",
    );
}
fn close_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-14),
        "actual={actual:?}, expected={expected:?}",
    );
}

#[test]
fn relative_dashpot_preserves_common_motion_and_damps_only_relative_motion() {
    let c = 1.7;
    let omega = 1.2;
    let m = model(vec![3.0, 3.0], vec![0.0, 0.0], 0.0, 0.0, None)
        .with_structural_viscosity(viscosity(vec![c, -c, -c, c], 2, 11.0)).unwrap();
    let common = m.frf(omega, &[C64::ONE, C64::ONE]).unwrap();
    let undamped = C64::from_re(1.0 / (9.0 - omega * omega));
    for &b in &common.b { close_complex(b, undamped, 2e-14); }
    assert!(common.power.structural < 1e-30, "{:?}", common.power);
    let relative = m.frf(omega, &[C64::ONE, -C64::ONE]).unwrap();
    let exact = C64::ONE / C64::new(9.0 - omega * omega, -2.0 * c * omega);
    close_complex(relative.b[0], exact, 2e-14);
    close_complex(relative.b[1], -exact, 2e-14);
    let expected_power = 2.0 * c * omega * omega * exact.norm_sq();
    close(relative.power.structural, expected_power, 2e-14);
    close(relative.power.input, expected_power, 2e-14);
    assert!(relative.power.balance_residual.abs() < 2e-14);
}

#[test]
fn full_viscosity_matches_independent_cavity_elimination_and_cycle_work() {
    let ws = [2.0, 3.6];
    let coupling = [0.25, -0.16];
    let c = [0.7, -0.49, -0.49, 0.343];
    let z = [
        C64::new(0.22, -0.12), C64::new(-0.06, 0.04),
        C64::new(-0.06, 0.04), C64::new(0.14, -0.08),
    ];
    let forces = [C64::new(1.0, 0.3), C64::new(-0.4, 0.2)];
    let m = model(ws.to_vec(), coupling.to_vec(), 0.0, 0.035, Some(z.to_vec()))
        .with_structural_viscosity(viscosity(c.to_vec(), 2, 3.0)).unwrap();
    for omega in [0.7, 1.9, 2.6, 4.3, 6.0] {
        // Eliminate the cavity analytically, then invert the remaining
        // complex 2x2 system by its determinant, independently of production LU.
        let da = C64::new(4.4 * 4.4 - omega * omega, -0.035 * 4.4 * 4.4);
        let gamma = 1.2 * 3.0 * 3.0 / 0.7 * omega * omega;
        let mut d = [C64::ZERO; 4];
        for r in 0..2 { for s in 0..2 {
            d[2*r+s] = C64::new(
                if r == s { ws[r]*ws[r] - omega*omega } else { 0.0 },
                -omega*c[2*r+s],
            ) - C64::new(0.0, omega)*z[2*r+s]
                - C64::from_re(gamma*coupling[r]*coupling[s])/da;
        }}
        let determinant = d[0]*d[3] - d[1]*d[2];
        let b = [
            (d[3]*forces[0] - d[1]*forces[1])/determinant,
            (d[0]*forces[1] - d[2]*forces[0])/determinant,
        ];
        let a = (b[0].scale(coupling[0]) + b[1].scale(coupling[1])).scale(gamma)/da;
        let actual = m.frf(omega, &forces).unwrap();
        close_complex(actual.b[0], b[0], 4e-13);
        close_complex(actual.b[1], b[1], 4e-13);
        close_complex(actual.a[0], a, 4e-13);

        // Integrate instantaneous physical dashpot work over a full period.
        // This uses the original signed C, independently of its factorization.
        let velocity = [C64::new(0.0, -omega)*b[0], C64::new(0.0, -omega)*b[1]];
        let mut cycle_power = 0.0;
        for sample in 0..64 {
            let phase = core::f64::consts::TAU * sample as f64 / 64.0;
            let (sine, cosine) = phase.sin_cos();
            let v = velocity.map(|v|v.re*cosine + v.im*sine);
            let instantaneous = v[0]*(c[0]*v[0] + c[1]*v[1])
                + v[1]*(c[2]*v[0] + c[3]*v[1]);
            assert!(instantaneous >= -1e-15);
            cycle_power += instantaneous / 64.0;
        }
        close(actual.power.structural, cycle_power, 4e-13);
        assert!(actual.power.structural > 0.0);
        assert!(actual.power.cavity > 0.0);
        assert!(actual.power.radiated > 0.0);
        assert!(actual.power.balance_residual.abs() < 4e-13, "{:?}", actual.power);
        assert!(actual.power.interface_residual.abs() < 4e-12, "{:?}", actual.power);

        let drive0 = m.frf(omega, &[C64::ONE, C64::ZERO]).unwrap();
        let drive1 = m.frf(omega, &[C64::ZERO, C64::ONE]).unwrap();
        close_complex(drive0.b[1], drive1.b[0], 4e-13);
    }
}

#[test]
fn signed_modal_basis_changes_preserve_response_and_all_power_owners() {
    let c = vec![0.8, -0.3, -0.3, 0.6];
    let z = vec![C64::new(0.4, -0.2), C64::new(0.1, 0.03),
                 C64::new(0.1, 0.03), C64::new(0.5, -0.1)];
    let force = [C64::new(0.9, 0.2), C64::new(-0.3, 0.5)];
    let m = model(vec![2.0, 3.6], vec![0.25, -0.16], 0.0, 0.035, Some(z.clone()))
        .with_structural_viscosity(viscosity(c.clone(), 2, 1.0)).unwrap();
    let signs = [-1.0, 1.0];
    let transformed_c = (0..4).map(|i|c[i]*signs[i/2]*signs[i%2]).collect();
    let transformed_z = (0..4).map(|i|z[i].scale(signs[i/2]*signs[i%2])).collect();
    let flipped = model(vec![2.0, 3.6], vec![-0.25, -0.16], 0.0, 0.035, Some(transformed_z))
        .with_structural_viscosity(viscosity(transformed_c, 2, 1e120)).unwrap();
    let actual = m.frf(2.3, &force).unwrap();
    let changed = flipped.frf(2.3, &[-force[0], force[1]]).unwrap();
    close_complex(changed.b[0], -actual.b[0], 2e-13);
    close_complex(changed.b[1], actual.b[1], 2e-13);
    close_complex(changed.a[0], actual.a[0], 2e-13);
    for (a, b) in [
        (changed.power.input, actual.power.input),
        (changed.power.structural, actual.power.structural),
        (changed.power.interface, actual.power.interface),
        (changed.power.cavity, actual.power.cavity),
        (changed.power.radiated, actual.power.radiated),
    ] { close(a, b, 2e-13); }
}

#[test]
fn truncation_uses_the_original_principal_matrix_and_channel_prefix() {
    let full = model(vec![2.0, 3.6, 5.7], vec![0.25, -0.16, 0.3], 0.0, 0.035, None)
        .with_structural_viscosity(viscosity(
            vec![2.0, -1.0, 0.2, -1.0, 1.5, -0.1, 0.2, -0.1, 1.0], 3, 7.0,
        )).unwrap();
    let reduced = model(vec![2.0, 3.6], vec![0.25, -0.16], 0.0, 0.035, None)
        .with_structural_viscosity(viscosity(vec![2.0, -1.0, -1.0, 1.5], 2, 0.02)).unwrap();
    let force = [C64::new(0.8, 0.2), C64::new(-0.1, 0.4), C64::new(3.0, -2.0)];
    for omega in [1.1, 2.7, 4.1] {
        let retained = full.frf_truncated(omega, &force, 2, 1).unwrap();
        let expected = reduced.frf(omega, &force[..2]).unwrap();
        for (&a, &b) in retained.b.iter().zip(&expected.b).chain(retained.a.iter().zip(&expected.a)) {
            close_complex(a, b, 2e-13);
        }
        close(retained.power.structural, expected.power.structural, 2e-13);
        assert!(retained.power.balance_residual.abs() < 2e-13);
    }
    // The convergence driver must use the same damping in both solves.
    let checked = full.frf_with_convergence(2.7, &force, 10.0).unwrap();
    let half = full.frf_truncated(2.7, &force, 1, 1).unwrap();
    close(checked.truncation_delta,
        (checked.response.power.input - half.power.input).abs()/checked.response.power.input.abs(),
        2e-14);
}

#[test]
fn structural_loss_ownership_and_basis_size_are_explicit() {
    let scalar = model(vec![2.0, 3.6], vec![0.25, -0.16], 0.02, 0.035, None);
    assert!(matches!(
        scalar.with_structural_viscosity(viscosity(vec![0.2, 0.0, 0.0, 0.3], 2, 1.0)),
        Err(VibroError::BadParameter { .. }),
    ));
    let unloaded = model(vec![2.0, 3.6], vec![0.25, -0.16], 0.0, 0.035, None);
    assert!(matches!(
        unloaded.clone().with_structural_viscosity(viscosity(vec![0.2], 1, 1.0)),
        Err(VibroError::ShapeMismatch { .. }),
    ));
    let zero = unloaded.clone()
        .with_structural_viscosity(viscosity(vec![0.0; 4], 2, 1.0)).unwrap();
    let force = [C64::new(0.8, 0.3), C64::new(0.1, -0.1)];
    let original = unloaded.frf(2.7, &force).unwrap();
    let explicit_zero = zero.frf(2.7, &force).unwrap();
    for (&a, &b) in original.b.iter().zip(&explicit_zero.b).chain(original.a.iter().zip(&explicit_zero.a)) {
        assert_eq!((a.re.to_bits(), a.im.to_bits()), (b.re.to_bits(), b.im.to_bits()));
    }
    assert_eq!(explicit_zero.power.structural, 0.0);
}
