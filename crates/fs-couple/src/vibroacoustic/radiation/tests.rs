//! G1 pulsating-sphere loaded motion; G0/G3 coherent work and passive directions.
use super::*;
use crate::vibroacoustic::{
    CavityModes, StructuralModes, VibroacousticModel, radiation_sample, with_radiation_cx,
};
use fs_modal::harmonic::HarmonicResponseError;

fn model(frequencies: &[f64], loading: ModalRadiationLoading) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: frequencies.to_vec(),
            shapes: vec![vec![1.0]; frequencies.len()],
            loss_factor: 0.0,
        },
        &CavityModes {
            omegas: vec![0.0],
            lambdas: vec![1.0],
            interface: vec![vec![1.0]],
            loss_factor: 0.0,
            rho0: 1.2,
            c0: 340.0,
        },
        vec![0.0; frequencies.len()],
        Some(loading),
    )
    .unwrap()
}

fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

#[test]
fn physical_sphere_reaction_changes_motion_and_closes_radiated_work() {
    // Exact spherical outgoing p(r)=C exp(i k r)/r gives Z_surface/(rho c)
    // = -i ka/(1-i ka). Integrating force over the surface gives this modal
    // load for a uniform radial displacement with unit mechanical mass.
    let (rho, sound, radius, dry_omega) = (1000.0, 1500.0, 0.05, 1000.0_f64);
    let area = 4.0 * core::f64::consts::PI * radius * radius;
    let force = C64::ONE;
    for omega in [100.0, 500.0, 1100.0, 3000.0] {
        let ka = omega * radius / sound;
        let z = C64::new(ka * ka, -ka).scale(rho * sound * area / (1.0 + ka * ka));
        let coupled = model(&[dry_omega], radiation_sample(omega, 1, vec![z]));
        let response = coupled.frf(omega, &[force]).unwrap();
        let exact = force
            / (C64::from_re(dry_omega * dry_omega - omega * omega)
                - C64::new(0.0, omega) * z);
        close(response.b[0], exact);
        let velocity = C64::new(0.0, -omega) * exact;
        let radiated = 0.5 * z.re * velocity.norm_sq();
        assert!((response.power.radiated / radiated - 1.0).abs() < 3e-13);
        assert!(response.power.balance_residual.abs() < 3e-12);
        assert!(response.power.radiated > 0.0);
        // The original constructor accepted a matrix at any frequency and
        // silently froze its added mass/resistance. All solve paths now refuse.
        assert!(coupled.frf(omega * 1.01, &[force]).is_err());
        assert!(
            coupled
                .frf_truncated(omega * 1.01, &[force], 1, 1)
                .is_err()
        );
        assert!(
            coupled
                .frf_with_convergence(omega * 1.01, &[force], 1.0)
                .is_err()
        );
    }
}

#[test]
fn complete_hermitian_admission_rejects_active_mutual_directions() {
    with_radiation_cx(|cx, _| {
        for matrix in [
            vec![
                C64::ONE,
                C64::from_re(2.0),
                C64::from_re(2.0),
                C64::ONE,
            ],
            vec![
                C64::ONE,
                C64::new(0.0, 2.0),
                C64::new(0.0, -2.0),
                C64::ONE,
            ],
            // A huge reactive impedance is not a tolerance for negative work.
            vec![
                C64::new(-1.0, 1e200),
                C64::ZERO,
                C64::ZERO,
                C64::ONE,
            ],
        ] {
            assert!(matches!(
                ModalRadiationLoading::try_new(100.0, 2, matrix, cx),
                Err(VibroError::Radiation(
                    HarmonicResponseError::NotPositiveSemidefinite { .. }
                )),
            ));
        }
        for matrix in [
            vec![],
            vec![C64::ONE; 3],
            vec![C64::new(f64::NAN, 0.0); 4],
        ] {
            assert!(ModalRadiationLoading::try_new(100.0, 2, matrix, cx).is_err());
        }
        assert!(ModalRadiationLoading::try_new(0.0, 1, vec![C64::ONE], cx).is_err());
    });
}

#[test]
fn reactive_exchange_cannot_pollute_lossless_or_weakly_radiating_work() {
    // Symmetric imaginary impedance is exactly lossless. Evaluate the work
    // before the coupled solve so a huge reactance cannot make its fields tiny
    // and hide cancellation of fictitious positive/negative radiation.
    let reactive = [
        C64::new(0.0, 3e150),
        C64::new(0.0, -7e150),
        C64::new(0.0, -7e150),
        C64::new(0.0, 5e150),
    ];
    let x = [C64::new(0.7, -0.4), C64::new(-0.2, 0.9)];
    let omega = 1.7;
    assert_eq!(outward_power(&reactive, 2, omega, &x), 0.0);
    let mut passive = reactive;
    let resistance = [2.0, -1.0, -1.0, 0.5];
    for (z, real) in passive.iter_mut().zip(resistance) {
        z.re = real;
    }
    let loading = radiation_sample(omega, 2, passive.to_vec());
    assert_eq!(loading.at_omega(omega).unwrap(), passive);
    let velocity_difference = (x[0] - x[1].scale(0.5)).scale(omega);
    let expected = velocity_difference.norm_sq();
    assert!((outward_power(&passive, 2, omega, &x) / expected - 1.0).abs() < 2e-14);
    let null_motion = [C64::new(1.0, 0.5), C64::new(2.0, 1.0)];
    assert_eq!(outward_power(&passive, 2, 1.0, &null_motion), 0.0);
    let retained = outward_power(&passive, 2, omega, &x[..1]);
    assert!((retained / (omega * omega * x[0].norm_sq()) - 1.0).abs() < 2e-14);
}

#[test]
fn mutual_radiation_keeps_coherent_cancellation_and_the_original_basis() {
    let omega = 2.0;
    let resistance = 0.7;
    let z = vec![
        C64::from_re(resistance),
        C64::from_re(-resistance),
        C64::from_re(-resistance),
        C64::from_re(resistance),
    ];
    let coupled = model(&[3.0, 3.0], radiation_sample(omega, 2, z));
    let common = coupled.frf(omega, &[C64::ONE, C64::ONE]).unwrap();
    close(common.b[0], C64::from_re(0.2));
    close(common.b[1], common.b[0]);
    assert!(common.power.radiated.abs() < 1e-28);
    let differential = coupled.frf(omega, &[C64::ONE, -C64::ONE]).unwrap();
    let exact = C64::new(5.0, -2.0 * omega * resistance).recip();
    close(differential.b[0], exact);
    close(differential.b[1], -exact);
    assert!(differential.power.radiated > 0.0);
    assert!(differential.power.balance_residual.abs() < 1e-13);
    with_radiation_cx(|cx, gate| {
        gate.request();
        assert!(matches!(
            ModalRadiationLoading::try_new(omega, 1, vec![C64::ONE], cx),
            Err(VibroError::Cancelled),
        ));
    });
}
