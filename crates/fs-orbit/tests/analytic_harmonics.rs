//! G1/G3: analytic AFT tangents retain the original physical periodic solve.

use fs_math::{c64::C64, det};
use fs_orbit::{
    HbAnchor, HbBudget, HbScales, IslandStateJacobian, OrbitError, OrbitProblem,
    solve_hb_scaled_seeded,
};
use std::cell::Cell;

#[derive(Clone, Copy)]
enum Tangent {
    Absent,
    Exact,
    Refuse,
    Incomplete,
}

struct Cascade {
    tangent: Tangent,
    island_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}
impl Cascade {
    fn new(tangent: Tangent) -> Self {
        Self {
            tangent,
            island_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }
    fn derivative(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        if matches!(self.tangent, Tangent::Refuse) {
            return false;
        }
        out[0] = -1.0;
        out[1] = 0.0;
        out[2] = 3.0 * x[0] * x[0];
        if !matches!(self.tangent, Tangent::Incomplete) {
            out[3] = -1.0;
        }
        true
    }
}
impl OrbitProblem for Cascade {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        self.island_calls.set(self.island_calls.get() + 1);
        out[0] = det::cos(t) - x[0];
        out[1] = x[0].powi(3) - x[1];
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        if matches!(self.tangent, Tangent::Absent) {
            None
        } else {
            Some(Self::derivative)
        }
    }
}

fn budget() -> HbBudget {
    HbBudget {
        harmonics: 3,
        tolerance: 1.0e-11,
        ..HbBudget::default()
    }
}
fn scales() -> HbScales {
    HbScales {
        states: vec![1.0; 2],
        equations: vec![1.0; 2],
        angular_frequency: 1.0,
    }
}

#[test]
fn exact_cubic_tangent_preserves_generated_harmonics_with_fewer_source_evaluations() {
    let exact = Cascade::new(Tangent::Exact);
    let probed = Cascade::new(Tangent::Absent);
    for source in [&exact, &probed] {
        let result = solve_hb_scaled_seeded(
            source,
            HbAnchor::Forced { omega: 1.0 },
            vec![0.0; 14],
            &budget(),
            &scales(),
        )
        .unwrap();
        for (state, harmonic, expected) in [
            (0, 1, C64::new(0.25, -0.25)),
            (1, 1, C64::new(0.0, -3.0 / 32.0)),
            (1, 3, C64::new(-1.0 / 80.0, 1.0 / 160.0)),
        ] {
            assert!((result.orbit.coeffs[state][harmonic] - expected).abs() < 1.0e-10);
        }
        assert!(result.scaled_residual < budget().tolerance);
    }
    assert!(exact.tangent_calls.get() > 0);
    assert!(exact.tangent_calls.get() < exact.island_calls.get());
    assert!(
        3 * exact.island_calls.get() < probed.island_calls.get(),
        "exact {}, finite difference {}",
        exact.island_calls.get(),
        probed.island_calls.get()
    );
    assert_eq!(probed.tangent_calls.get(), 0);
}

#[test]
fn declared_tangent_refusal_or_unwritten_entries_never_fall_back_to_probes() {
    for mode in [Tangent::Refuse, Tangent::Incomplete] {
        let source = Cascade::new(mode);
        assert_eq!(
            solve_hb_scaled_seeded(
                &source,
                HbAnchor::Forced { omega: 1.0 },
                vec![0.0; 14],
                &budget(),
                &scales()
            )
            .unwrap_err(),
            OrbitError::TangentRefused { sample: 0 }
        );
        assert_eq!(
            source.island_calls.get(),
            16,
            "only the original base residual was evaluated"
        );
        assert_eq!(source.tangent_calls.get(), 1);
    }
}

#[test]
fn exact_backbone_tangent_preserves_physical_scales_frequency_and_deflation() {
    struct Backbone {
        calls: Cell<usize>,
    }
    impl OrbitProblem for Backbone {
        fn dim(&self) -> usize {
            2
        }
        fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
            out[0] = x[1];
            out[1] = -4.0 * x[0];
        }
        fn port(&self, s: C64) -> Vec<C64> {
            vec![s.scale(1.0e18), C64::ZERO, C64::ZERO, s.scale(1.0e-18)]
        }
        fn autonomous(&self) -> bool {
            true
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|source, _, _, out| {
                source.calls.set(source.calls.get() + 1);
                out.copy_from_slice(&[0.0, 1.0, -4.0, 0.0]);
                true
            })
        }
    }
    let scales = HbScales {
        states: vec![1.0e-12, 1.0e6],
        equations: vec![1.0e6, 1.0e-12],
        angular_frequency: 2.0,
    };
    for deflate_equilibrium in [false, true] {
        for phase in [0.0, 0.7] {
            let source = Backbone {
                calls: Cell::new(0),
            };
            let x = C64::new(det::cos(phase), det::sin(phase)).scale(0.4);
            let v = C64::new(0.0, 2.0) * x;
            let seed = vec![
                0.0,
                0.0,
                x.re * 1.0e-12,
                v.re * 1.0e6,
                x.im * 1.0e-12,
                v.im * 1.0e6,
                2.6,
            ];
            let result = solve_hb_scaled_seeded(
                &source,
                HbAnchor::Backbone {
                    amplitude: 0.4e-12,
                    omega_guess: 2.6,
                },
                seed,
                &HbBudget {
                    harmonics: 1,
                    tolerance: 1.0e-10,
                    deflate_equilibrium,
                    ..HbBudget::default()
                },
                &scales,
            )
            .unwrap();
            assert!((result.orbit.first_harmonic_amplitude(0) / 0.4e-12 - 1.0).abs() < 1.0e-9);
            assert!((result.orbit.omega - 2.0).abs() < 1.0e-9);
            assert!(
                (result.orbit.coeffs[1][1].scale(1.0e-6)
                    - C64::new(0.0, 2.0) * result.orbit.coeffs[0][1].scale(1.0e12))
                .abs()
                    < 1.0e-9
            );
            assert!(
                source.calls.get() > 0,
                "deflation must retain the exact source tangent"
            );
        }
    }
}
