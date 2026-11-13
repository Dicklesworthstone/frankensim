//! G1/G3: actual radiation, nonlinear storage/memory, mutual loads and work.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::phs_orbit::{
    PeriodicPortInput, PhsOrbitError, audit_forced_orbit,
    impedance::{
        PeriodicPortImpedance, PhsImpedanceEvidence, PortImpedanceError, PortImpedanceLimits,
    },
    solve_forced_harmonic_balance_scaled_seeded,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_orbit::{HbBudget, HbRefinementBudget, HbScales, OrbitError, OrbitProblem};
use fs_phs::{
    ContinuousStorageTangent, PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage,
};
use std::{cell::Cell, rc::Rc};

const LIMITS: PortImpedanceLimits = PortImpedanceLimits {
    maximum_states: 16,
    maximum_ports: 8,
    maximum_harmonics: 31,
    maximum_hb_unknowns: 1024,
    maximum_products: 20_000_000,
    maximum_workspace_bytes: 4 << 20,
};
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 27,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn hb(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 60,
        tolerance: 1e-11,
        deflate_equilibrium: false,
    }
}
fn scales(dimension: usize, omega: f64) -> HbScales {
    HbScales {
        states: vec![1.0; dimension],
        equations: vec![1.0; dimension],
        angular_frequency: omega,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_actual_baffled_piston_changes_the_orbit_and_its_physical_work() {
    with_cx(&CancelGate::new(), |cx| {
        let (mass, stiffness, damping, omega, radius) = (0.01, 2500.0, 0.03, 500.0, 0.05);
        let area = std::f64::consts::PI * radius * radius;
        let source = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let force = C64::new(0.2, -0.07); // physical exp(-i omega t) peak phasor
        let input = PeriodicPortInput::new(
            omega,
            vec![vec![C64::from_re(0.04), force.conj().scale(0.5)]],
        )
        .unwrap();
        // The exact physical piston owner supplies every nonzero frequency.
        // Its finite DC limit is zero; pressure/volume-flow converts to force/
        // velocity with area^2, keeping the physical port work unchanged.
        let mut matrices = vec![vec![C64::ZERO]];
        for harmonic in 1..=3 {
            let (r, x) =
                fs_phs::baffled_piston_impedance(1.2, 343.0, radius, harmonic as f64 * omega)
                    .unwrap();
            matrices.push(vec![C64::new(r, x).scale(area * area)]);
        }
        let z = matrices[1][0];
        assert!(z.re > 0.0 && z.im < 0.0);
        let load = PeriodicPortImpedance::new(omega, 1, matrices, LIMITS, cx).unwrap();
        let prepared = load.prepare(&source, &input, cx).unwrap();
        assert!(prepared.time_domain_rhs().is_none());
        assert!(prepared.time_domain_state_jacobian().is_none());
        let result = prepared
            .solve_scaled_seeded(vec![0.0; 3 * 7], &hb(3), &scales(3, omega))
            .unwrap();
        let q = force
            / (C64::new(stiffness - mass * omega * omega, -damping * omega)
                + C64::new(0.0, -omega) * z);
        let v = C64::new(0.0, -omega) * q;
        complex(result.state.coeffs[0][1].conj().scale(2.0), q, 2e-10);
        complex(
            result.state.coeffs[1][1].conj().scale(2.0),
            v.scale(mass),
            2e-10,
        );
        complex(result.flow_coefficients[0][1].conj().scale(2.0), v, 2e-9);
        complex(
            result.effort_coefficients[0][1].conj().scale(2.0),
            z * v,
            2e-9,
        );
        close(result.state.coeffs[0][0].re, 0.04 / stiffness, 2e-12);
        let unloaded = force / C64::new(stiffness - mass * omega * omega, -damping * omega);
        assert!(
            q.abs() < 0.5 * unloaded.abs(),
            "radiation inertia must change the resonant response"
        );
        let a = prepared.audit(&result, 128).unwrap();
        let b = prepared.audit(&result, 257).unwrap();
        close(
            a.mean_boundary_absorbed_power_w,
            0.5 * z.re * v.norm_sq(),
            2e-10,
        );
        close(
            a.mean_native_resistance_loss_w,
            0.5 * damping * v.norm_sq(),
            2e-10,
        );
        close(a.mean_source_power_w, 0.5 * (force.conj() * v).re, 2e-10);
        close(
            a.mean_boundary_absorbed_power_w,
            b.mean_boundary_absorbed_power_w,
            2e-12,
        );
        assert!(a.cycle_power_residual_w.abs() < 1e-10);
        assert!(a.physical_flow_rms_defect[0] < 1e-10);
    });
}

// Retain the actual fs-phs nonlinear/Maxwell storage while exposing its
// work-conjugate strain-rate / restoring-effort port. No replacement law.
struct NativeStorage(PortHamiltonian);
impl Storage for NativeStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0.storage().hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.0.storage().gradient(x, out);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        self.0.storage().continuous_tangent()
    }
}
fn nonlinear_memory(extra_resistance: f64) -> PortHamiltonian {
    let native = fs_phs::duffing_oscillator(2.0, 4.0, 7.0, 0.3)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 3.0,
            relaxation_time_s: 0.4,
        }])
        .unwrap();
    let (j, r, _) = native.structure();
    let j = j.to_vec();
    let mut r = r.to_vec();
    r[0] += extra_resistance;
    PortHamiltonian::new(
        3,
        1,
        j,
        r,
        vec![1.0, 0.0, 0.0],
        Box::new(NativeStorage(native)),
    )
    .unwrap()
}

#[test]
fn g1_nonlinear_flow_and_maxwell_states_match_the_same_physical_resistance() {
    with_cx(&CancelGate::new(), |cx| {
        let (omega, resistance) = (1.3, 0.2);
        let source = nonlinear_memory(0.0);
        let reference = nonlinear_memory(resistance);
        let input =
            PeriodicPortInput::new(omega, vec![vec![C64::from_re(0.1), C64::new(0.2, -0.04)]])
                .unwrap();
        let load = PeriodicPortImpedance::new(
            omega,
            1,
            vec![vec![C64::from_re(resistance)]; 16],
            LIMITS,
            cx,
        )
        .unwrap();
        let prepared = load.prepare(&source, &input, cx).unwrap();
        let callback = prepared
            .island_state_jacobian()
            .expect("actual continuous nonlinear memory tangent");
        let mut tangent = vec![f64::NAN; 16];
        assert!(callback(
            &prepared,
            0.0,
            &[0.2, -0.1, 0.05, 0.3],
            &mut tangent
        ));
        // The nonlinear ORIGINAL output is grad_q H; y columns are zero in
        // the island because the complete boundary reaction lives in P(s).
        close(tangent[3 * 4], 4.0 + 21.0 * 0.2 * 0.2 + 3.0, 1e-13);
        close(tangent[3 * 4 + 2], -3.0_f64.sqrt(), 1e-13);
        assert!((0..4).all(|row| tangent[row * 4 + 3] == 0.0));
        let refinement = HbRefinementBudget {
            max_harmonics: 15,
            residual_tolerances: vec![2e-7; 4],
            state_tolerances: vec![2e-7; 4],
            frequency_tolerance: 1e-10,
        };
        let result = prepared
            .solve_adaptive_scaled_seeded(vec![0.0; 4 * 7], &hb(3), &refinement, &scales(4, omega))
            .unwrap();
        let PhsImpedanceEvidence::Adaptive(checks) = &result.evidence else {
            panic!("adaptive evidence");
        };
        assert!(checks.len() >= 2);
        let harmonics = result.state.coeffs[0].len() - 1;
        let expected = solve_forced_harmonic_balance_scaled_seeded(
            &reference,
            &input,
            vec![0.0; 3 * (2 * harmonics + 1)],
            &hb(harmonics),
            &scales(3, omega),
        )
        .unwrap();
        assert_eq!(
            result.state.coeffs.len(),
            3,
            "original Maxwell coordinate must remain"
        );
        for coordinate in 0..3 {
            for h in 0..=harmonics {
                complex(
                    result.state.coeffs[coordinate][h],
                    expected.orbit.coeffs[coordinate][h],
                    2e-9,
                );
            }
        }
        for h in 0..=harmonics {
            complex(
                result.state.coeffs[2][h],
                result.state.coeffs[0][h].scale(3.0_f64.sqrt())
                    / C64::new(1.0, h as f64 * omega * 0.4),
                2e-9,
            );
        }
        assert!(
            result.flow_coefficients[0][3].abs() > 1e-5,
            "nonlinear output harmonics must load the source"
        );
        let a = prepared.audit(&result, 257).unwrap();
        let b = prepared.audit(&result, 513).unwrap();
        let reference_work = audit_forced_orbit(&reference, &input, &expected.orbit, 513).unwrap();
        assert!(a.mean_native_resistance_loss_w > 0.0 && a.mean_boundary_absorbed_power_w > 0.0);
        close(
            a.mean_native_resistance_loss_w + a.mean_boundary_absorbed_power_w,
            reference_work.mean_resistance_loss_w,
            2e-9,
        );
        close(
            a.mean_boundary_absorbed_power_w,
            a.retained_boundary_absorbed_power_w,
            2e-9,
        );
        close(a.mean_source_power_w, b.mean_source_power_w, 2e-10);
        assert!(a.cycle_power_residual_w.abs() < 2e-9);
        assert!(a.physical_flow_rms_defect[0] < 2e-7);
    });
}

#[test]
fn g3_full_mutual_nonreciprocal_load_preserves_dc_work_and_unit_changes() {
    with_cx(&CancelGate::new(), |cx| {
        let omega = 1.7;
        let dc = [1.2, -0.4];
        let force = [C64::new(0.5, -0.2), C64::new(-0.1, 0.3)];
        // Two translational momenta: H=p1^2/(2m1)+p2^2/(2m2).
        // The load contains dissipative mutual work, an antisymmetric lossless
        // port coupling, and reciprocal reactive inertia. No entry is dropped.
        let matrix = |w: f64| {
            vec![
                C64::new(0.7, -0.1 * w),
                C64::new(1.1, -0.04 * w),
                C64::new(-0.5, -0.04 * w),
                C64::new(0.5, -0.2 * w),
            ]
        };
        let z0 = matrix(0.0);
        let z = matrix(omega);
        let mut recorded: Option<(Vec<C64>, f64)> = None;
        for transformer in [1.0, 1e-3, 1e3] {
            let system = PortHamiltonian::new(
                2,
                2,
                vec![0.0; 4],
                vec![0.2, 0.0, 0.0, 0.4],
                vec![transformer, 0.0, 0.0, transformer],
                Box::new(QuadraticStorage::new(vec![0.5, 0.0, 0.0, 1.0], 2).unwrap()),
            )
            .unwrap();
            let input = PeriodicPortInput::new(
                omega,
                (0..2)
                    .map(|a| {
                        vec![
                            C64::from_re(dc[a] / transformer),
                            force[a].conj().scale(0.5 / transformer),
                        ]
                    })
                    .collect(),
            )
            .unwrap();
            let matrices = (0..=3)
                .map(|h| {
                    matrix(h as f64 * omega)
                        .into_iter()
                        .map(|value| value.scale(1.0 / (transformer * transformer)))
                        .collect()
                })
                .collect();
            let load = PeriodicPortImpedance::new(omega, 2, matrices, LIMITS, cx).unwrap();
            let prepared = load.prepare(&system, &input, cx).unwrap();
            let units = HbScales {
                states: vec![1.0, 1.0, transformer, transformer],
                equations: vec![1.0, 1.0, transformer, transformer],
                angular_frequency: omega,
            };
            let result = prepared
                .solve_scaled_seeded(vec![0.0; 4 * 7], &hb(3), &units)
                .unwrap();
            let solve = |a: C64, b: C64, c: C64, d: C64, f: [C64; 2]| {
                let determinant = a * d - b * c;
                [
                    (d * f[0] - b * f[1]) / determinant,
                    (a * f[1] - c * f[0]) / determinant,
                ]
            };
            let vdc = solve(
                z0[0] + C64::from_re(0.2),
                z0[1],
                z0[2],
                z0[3] + C64::from_re(0.4),
                [C64::from_re(dc[0]), C64::from_re(dc[1])],
            );
            let v = solve(
                z[0] + C64::new(0.2, -2.0 * omega),
                z[1],
                z[2],
                z[3] + C64::new(0.4, -omega),
                force,
            );
            for a in 0..2 {
                complex(
                    result.flow_coefficients[a][0].scale(1.0 / transformer),
                    vdc[a],
                    2e-10,
                );
                complex(
                    result.flow_coefficients[a][1]
                        .conj()
                        .scale(2.0 / transformer),
                    v[a],
                    2e-10,
                );
            }
            let audit = prepared.audit(&result, 129).unwrap();
            let expected_source = dc[0] * vdc[0].re
                + dc[1] * vdc[1].re
                + 0.5 * (force[0].conj() * v[0] + force[1].conj() * v[1]).re;
            close(audit.mean_source_power_w, expected_source, 2e-10);
            assert!(audit.mean_boundary_absorbed_power_w > 0.0);
            assert!(audit.cycle_power_residual_w.abs() < 2e-10);
            if let Some((state, work)) = &recorded {
                for (actual, expected) in result.state.coeffs.iter().flatten().zip(state) {
                    complex(*actual, *expected, 2e-10);
                }
                close(audit.mean_boundary_absorbed_power_w, *work, 2e-10);
            } else {
                recorded = Some((
                    result
                        .state
                        .coeffs
                        .iter()
                        .flatten()
                        .copied()
                        .collect::<Vec<_>>(),
                    audit.mean_boundary_absorbed_power_w,
                ));
            }
        }
    });
}

struct CountedStorage(Rc<Cell<usize>>);
impl Storage for CountedStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * x[0] * x[0]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.0.set(self.0.get() + 1);
        out[0] = x[0];
    }
}
#[test]
fn g1_missing_bands_nonpassive_or_complex_dc_refuse_before_source_evaluation() {
    with_cx(&CancelGate::new(), |cx| {
        let calls = Rc::new(Cell::new(0));
        let source = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![1.0],
            vec![1.0],
            Box::new(CountedStorage(calls.clone())),
        )
        .unwrap();
        let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO, C64::ONE]]).unwrap();
        let table =
            PeriodicPortImpedance::new(2.0, 1, vec![vec![C64::ONE]; 2], LIMITS, cx).unwrap();
        let prepared = table.prepare(&source, &input, cx).unwrap();
        assert!(matches!(
            prepared.solve_scaled_seeded(vec![0.0; 10], &hb(2), &scales(2, 2.0)),
            Err(PortImpedanceError::Invalid(_))
        ));
        let refinement = HbRefinementBudget {
            max_harmonics: 3,
            residual_tolerances: vec![1e-6; 2],
            state_tolerances: vec![1e-6; 2],
            frequency_tolerance: 1e-8,
        };
        assert!(matches!(
            prepared.solve_adaptive_scaled_seeded(
                vec![0.0; 6],
                &hb(1),
                &refinement,
                &scales(2, 2.0)
            ),
            Err(PortImpedanceError::Invalid(_))
        ));
        assert_eq!(calls.get(), 0);
        assert!(matches!(
            PeriodicPortImpedance::new(
                2.0,
                1,
                vec![vec![C64::new(1.0, 1e-100)], vec![C64::ONE]],
                LIMITS,
                cx
            ),
            Err(PortImpedanceError::Invalid(_))
        ));
        // Positive individual resistances cannot hide a negative mutual direction.
        assert!(matches!(
            PeriodicPortImpedance::new(
                2.0,
                2,
                vec![vec![C64::ONE, C64::from_re(2.0), C64::from_re(2.0), C64::ONE]; 2],
                LIMITS,
                cx
            ),
            Err(PortImpedanceError::Admission { harmonic: 0, .. })
        ));
        assert!(matches!(
            PeriodicPortImpedance::new(2.0, 1, vec![vec![C64::ONE], vec![]], LIMITS, cx),
            Err(PortImpedanceError::Invalid(_))
        ));
    });
}

#[test]
fn g1_declared_tangent_refusal_and_cancel_cannot_become_an_unloaded_solution() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let source = fs_phs::mass_spring_damper(2.0, 4.0, 0.3).unwrap();
        let input = PeriodicPortInput::new(1.3, vec![vec![C64::ZERO, C64::ONE]]).unwrap();
        let table =
            PeriodicPortImpedance::new(1.3, 1, vec![vec![C64::ONE]; 2], LIMITS, cx).unwrap();
        let refused = |_: &[f64], _: &[f64], _: &mut [f64]| false;
        let prepared = table
            .prepare(&source, &input, cx)
            .unwrap()
            .with_storage_hessian(&refused);
        assert!(matches!(
            prepared.solve_scaled_seeded(vec![0.0; 9], &hb(1), &scales(3, 1.3)),
            Err(PortImpedanceError::Regime(PhsOrbitError::Orbit(
                OrbitError::TangentRefused { .. }
            )))
        ));
        gate.request();
        assert!(matches!(
            prepared.solve_scaled_seeded(vec![0.0; 9], &hb(1), &scales(3, 1.3)),
            Err(PortImpedanceError::Cancelled)
        ));
    });
}
