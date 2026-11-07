//! The original material law supplies every memory time and stiffness.

use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_modal::relaxation_modes::{MaxwellDampingBranch, maxwell_eigenvalues};

#[test]
fn material_loss_target_drives_an_axial_resonator_without_sampling_modal_zeta() {
    // The material owner generates this Prony source from a synthetic loss
    // target. In a uniform axial member, strain=q/L and energy=V E strain²/2.
    let source = GeneralizedMaxwell::matching_loss(2e6, 50.0, 0.2).unwrap();
    let (length, area, mass, viscous_drag) = (0.5, 0.002, 2.0, 0.3);
    let volume = area * length;
    let relaxed_stiffness = source.e_inf * area / length;
    let branches: Vec<_> = source
        .terms
        .iter()
        .filter(|(delta, _)| *delta > 0.0)
        .map(|&(delta, tau)| MaxwellDampingBranch {
            strain: vec![1.0 / length],
            stiffness: volume * delta,
            relaxation_time_s: tau,
        })
        .collect();
    let poles =
        maxwell_eigenvalues(&[mass], &[viscous_drag], &[relaxed_stiffness], 1, &branches).unwrap();
    assert_eq!(poles.len(), 3);
    let (delta, tau) = source.terms[0];
    let coefficients = [
        mass * tau,
        mass + viscous_drag * tau,
        viscous_drag + (source.e_inf + delta) * area / length * tau,
        relaxed_stiffness,
    ];
    for s in poles {
        let mut residual = C64::ZERO;
        let mut scale = 0.0;
        for a in coefficients {
            residual = residual * s + C64::from_re(a);
            scale = scale * s.abs() + a.abs();
        }
        assert!(residual.abs() < 2e-11 * scale);
        assert!(s.re < 0.0);
    }
    // The retained material variables also reproduce the same owner's full
    // frequency-dependent axial stiffness, including both asymptotic moduli.
    for omega in [0.0, 5.0, 50.0, 500.0, 5e4] {
        let (storage, loss) = source.modulus(omega);
        let mut reconstructed = C64::from_re(relaxed_stiffness);
        for b in &branches {
            let iw_tau = C64::new(0.0, omega * b.relaxation_time_s);
            reconstructed = reconstructed
                + (iw_tau / (C64::ONE + iw_tau)).scale(b.stiffness * b.strain[0] * b.strain[0]);
        }
        let expected = C64::new(storage, loss).scale(area / length);
        assert!((reconstructed - expected).abs() < 3e-12 * expected.abs());
    }
}
