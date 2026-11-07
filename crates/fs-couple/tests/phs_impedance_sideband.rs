//! G1/G3: loaded original-state sidebands, signed dispersive tangents and DC.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::phs_orbit::{
    PeriodicPortInput,
    impedance::{
        PeriodicPortImpedance, PortImpedanceError, PortImpedanceLimits,
        sideband::{
            ImpedanceSidebandLimits, PhsImpedanceSidebandSource, PortImpedanceFrequencySample,
            SidebandPortImpedance,
        },
    },
    sideband::{PhsPortSidebandControls, PhsPortSidebandError, PhsPortSidebandSource},
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_orbit::{
    HbBudget, HbScales,
    sideband::{
        SidebandControls, SidebandError, SidebandScales,
        sensitivity::SidebandCarrierDerivativeControls,
    },
};
use fs_phs::{
    ContinuousStorageTangent, PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage,
};

const LOAD_LIMITS: PortImpedanceLimits = PortImpedanceLimits {
    maximum_states: 16,
    maximum_ports: 8,
    maximum_harmonics: 31,
    maximum_hb_unknowns: 1024,
    maximum_products: 50_000_000,
    maximum_workspace_bytes: 8 << 20,
};
const PROBE_LIMITS: ImpedanceSidebandLimits = ImpedanceSidebandLimits {
    maximum_frequency_samples: 128,
    maximum_carriers: 16,
    maximum_products: 50_000_000,
    maximum_workspace_bytes: 8 << 20,
};
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 81,
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
fn hb_scales(dimension: usize, omega: f64) -> HbScales {
    HbScales {
        states: vec![1.0; dimension],
        equations: vec![1.0; dimension],
        angular_frequency: omega,
    }
}
fn scales(dimension: usize) -> SidebandScales {
    SidebandScales {
        states: vec![1.0; dimension],
        equations: vec![1.0; dimension],
    }
}
fn controls(ports: usize, harmonics: usize) -> PhsPortSidebandControls {
    PhsPortSidebandControls {
        sideband: SidebandControls {
            harmonics,
            samples_per_period: 128,
            maximum_unknowns: 512,
            maximum_payload_bytes: 64 << 20,
            maximum_work: 1_000_000_000,
            maximum_base_residual: 1e-7,
            maximum_tangent_sampling_defect: 1e-8,
            maximum_relative_residual: 1e-10,
            maximum_omitted_residual: 1e-7,
        },
        output_harmonics: harmonics,
        initial_samples: 128,
        max_samples: 1024,
        max_hessian_actions: 50_000,
        max_workspace_bytes: 16 << 20,
        max_observation_work: 1_000_000_000,
        envelope_tolerances: vec![1e-8; ports],
    }
}
fn derivative_controls() -> SidebandCarrierDerivativeControls {
    SidebandCarrierDerivativeControls {
        maximum_payload_bytes: 16 << 20,
        maximum_work: 500_000_000,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-7,
    }
}
fn complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}
fn samples(
    omega: f64,
    carriers: &[f64],
    harmonics: usize,
    mut data: impl FnMut(f64) -> (Vec<C64>, Option<Vec<C64>>),
) -> Vec<PortImpedanceFrequencySample> {
    let mut frequencies: Vec<_> = carriers
        .iter()
        .flat_map(|&carrier| {
            (-(harmonics as isize)..=harmonics as isize)
                .map(move |k| (k as f64).mul_add(omega, carrier).abs())
        })
        .collect();
    frequencies.sort_by(f64::total_cmp);
    frequencies.dedup();
    frequencies
        .into_iter()
        .map(|angular_frequency| {
            let (impedance, angular_derivative) = data(angular_frequency);
            PortImpedanceFrequencySample {
                angular_frequency,
                impedance,
                angular_derivative,
            }
        })
        .collect()
}
// Physical resistance, one relaxation branch and inertance, exp(-iwt).
fn dispersive(w: f64) -> (C64, C64) {
    let pole = C64::new(1.0, -0.4 * w);
    (
        C64::new(0.3, -0.2 * w) + C64::from_re(0.7) / pole,
        C64::new(0.0, -0.2) + C64::new(0.0, 0.28) / (pole * pole),
    )
}

#[test]
fn dispersive_load_and_same_lu_derivative_match_both_signs_and_dc_crossings() {
    with_cx(|cx| {
        let (mass, stiffness, damping, omega) = (2.0, 5.0, 0.2, 2.0);
        let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let input =
            PeriodicPortInput::new(omega, vec![vec![C64::from_re(0.1), C64::new(0.12, -0.03)]])
                .unwrap();
        let load = PeriodicPortImpedance::new(
            omega,
            1,
            (0..=1)
                .map(|h| vec![dispersive(h as f64 * omega).0])
                .collect(),
            LOAD_LIMITS,
            cx,
        )
        .unwrap();
        let mut loaded = load.prepare(&system, &input, cx).unwrap();
        let base = loaded
            .solve_scaled_seeded(vec![0.0; 9], &hb(1), &hb_scales(3, omega))
            .unwrap();
        let carriers = vec![-3.0, -2.0, 0.0, 2.0, 3.0];
        let table = SidebandPortImpedance::new(
            &load,
            carriers.clone(),
            1,
            samples(omega, &carriers, 1, |w| {
                let (z, dz) = dispersive(w);
                (vec![z], Some(vec![dz]))
            }),
            PROBE_LIMITS,
            cx,
        )
        .unwrap();
        let mut source =
            PhsImpedanceSidebandSource::new(&mut loaded, &base, &table, &controls(1, 1)).unwrap();
        let mut tangent = source.prepare(&scales(3), || false).unwrap();
        assert!(tangent.base_residual() < 1e-9);
        let drive = [C64::new(0.2, 0.1), C64::new(1.0, -0.2), C64::new(-0.3, 0.2)];
        let drive_slope = [
            C64::new(0.03, -0.01),
            C64::new(-0.02, 0.04),
            C64::new(0.01, 0.02),
        ];
        for carrier in carriers {
            let mut factor = tangent.factor(carrier, || false).unwrap();
            let pair = factor
                .solve_with_carrier_derivative(
                    &drive,
                    Some(&drive_slope),
                    derivative_controls(),
                    &[1e-8],
                    || false,
                )
                .unwrap();
            let again = factor.solve(&drive, || false).unwrap();
            for band in 0..3 {
                let w = ((band as isize - 1) as f64).mul_add(omega, carrier);
                // Direct analytic expression at SIGNED w; no table convention reuse.
                let pole = C64::new(1.0, 0.4 * w);
                let z = C64::new(0.3, 0.2 * w) + C64::from_re(0.7) / pole;
                let dz = C64::new(0.0, 0.2) + C64::new(0.0, -0.28) / (pole * pole);
                let iw = C64::new(0.0, w);
                let denominator = C64::new(stiffness - mass * w * w, damping * w) + iw * z;
                let derivative =
                    C64::new(-2.0 * mass * w, damping) + C64::new(0.0, 1.0) * z + iw * dz;
                let q = drive[band] / denominator;
                let dq = (drive_slope[band] - derivative * q) / denominator;
                let velocity = iw * q;
                let dv = C64::new(0.0, 1.0) * q + iw * dq;
                complex(
                    pair.response.response.state_response.coefficients[2 * band],
                    q,
                    2e-9,
                );
                complex(
                    pair.response.response.state_response.coefficients[2 * band + 1],
                    velocity.scale(mass),
                    2e-9,
                );
                complex(
                    pair.response.state_angular_derivative.coefficients[2 * band],
                    dq,
                    2e-9,
                );
                complex(
                    pair.response.state_angular_derivative.coefficients[2 * band + 1],
                    dv.scale(mass),
                    2e-9,
                );
                complex(
                    pair.response.response.port_coefficients[band],
                    velocity,
                    2e-9,
                );
                complex(pair.response.port_angular_derivative[band], dv, 2e-9);
                complex(pair.flow_coefficients[band], velocity, 2e-9);
                complex(pair.flow_angular_derivative[band], dv, 2e-9);
                complex(pair.effort_coefficients[band], z * velocity, 2e-9);
                complex(
                    pair.effort_angular_derivative[band],
                    z * dv + dz * velocity,
                    2e-9,
                );
                complex(
                    again.flow_coefficients[band],
                    pair.flow_coefficients[band],
                    2e-11,
                );
            }
            assert!(pair.retained_flow_mismatch[0] < 1e-10);
            assert!(pair.derivative_flow_mismatch[0] < 1e-10);
        }
    });
}

#[test]
fn complete_mutual_nonreciprocal_load_preserves_conjugate_port_units() {
    with_cx(|cx| {
        let omega = 2.0;
        let carrier = 1.25;
        let z = [
            C64::from_re(0.6),
            C64::from_re(0.28),
            C64::from_re(-0.12),
            C64::from_re(0.9),
        ];
        let forcing = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
        let a = [
            z[0] + C64::new(0.2, carrier),
            z[1] - C64::ONE,
            z[2] + C64::ONE,
            z[3] + C64::new(0.3, carrier),
        ];
        let determinant = a[0] * a[3] - a[1] * a[2];
        let expected = [
            (a[3] * forcing[0] - a[1] * forcing[1]) / determinant,
            (a[0] * forcing[1] - a[2] * forcing[0]) / determinant,
        ];
        for units in [[1.0, 1.0], [10.0, 0.1]] {
            let system = PortHamiltonian::new(
                2,
                2,
                vec![0.0, 1.0, -1.0, 0.0],
                vec![0.2, 0.0, 0.0, 0.3],
                vec![units[0], 0.0, 0.0, units[1]],
                Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
            )
            .unwrap();
            let transformed: Vec<_> = (0..4)
                .map(|i| z[i].scale(1.0 / (units[i / 2] * units[i % 2])))
                .collect();
            let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO]; 2]).unwrap();
            let load =
                PeriodicPortImpedance::new(omega, 2, vec![transformed.clone(); 2], LOAD_LIMITS, cx)
                    .unwrap();
            let mut loaded = load.prepare(&system, &input, cx).unwrap();
            let base = loaded
                .solve_scaled_seeded(vec![0.0; 12], &hb(1), &hb_scales(4, omega))
                .unwrap();
            let table = SidebandPortImpedance::new(
                &load,
                vec![carrier],
                0,
                samples(omega, &[carrier], 0, |_| (transformed.clone(), None)),
                PROBE_LIMITS,
                cx,
            )
            .unwrap();
            let mut source =
                PhsImpedanceSidebandSource::new(&mut loaded, &base, &table, &controls(2, 0))
                    .unwrap();
            let mut tangent = source.prepare(&scales(4), || false).unwrap();
            let mut factor = tangent.factor(carrier, || false).unwrap();
            let result = factor
                .solve(
                    &[
                        forcing[0].scale(1.0 / units[0]),
                        forcing[1].scale(1.0 / units[1]),
                    ],
                    || false,
                )
                .unwrap();
            for p in 0..2 {
                complex(
                    result.response.state_response.coefficients[p],
                    expected[p],
                    2e-10,
                );
                complex(
                    result.response.port_coefficients[p],
                    expected[p].scale(units[p]),
                    2e-10,
                );
                complex(
                    result.flow_coefficients[p],
                    expected[p].scale(units[p]),
                    2e-10,
                );
                complex(
                    result.effort_coefficients[p],
                    (z[p * 2] * expected[0] + z[p * 2 + 1] * expected[1]).scale(1.0 / units[p]),
                    2e-10,
                );
            }
        }
    });
}

// Delegate to the actual nonlinear Duffing/Maxwell storage and native Hessian.
// Change only the original work-conjugate port, exposing nonlinear restoring flow.
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
fn nonlinear_memory(resistance: f64) -> PortHamiltonian {
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
    r[0] += resistance;
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
fn nonlinear_loaded_memory_mixes_sidebands_like_the_same_native_resistance() {
    with_cx(|cx| {
        let (omega, carrier, resistance, h) = (1.3, 0.7, 0.2, 5);
        let system = nonlinear_memory(0.0);
        let reference = nonlinear_memory(resistance);
        let input = PeriodicPortInput::new(
            omega,
            vec![vec![C64::from_re(0.01), C64::new(0.04, -0.008)]],
        )
        .unwrap();
        let load = PeriodicPortImpedance::new(
            omega,
            1,
            vec![vec![C64::from_re(resistance)]; 12],
            LOAD_LIMITS,
            cx,
        )
        .unwrap();
        let mut loaded = load.prepare(&system, &input, cx).unwrap();
        let base = loaded
            .solve_scaled_seeded(vec![0.0; 4 * 23], &hb(11), &hb_scales(4, omega))
            .unwrap();
        let table = SidebandPortImpedance::new(
            &load,
            vec![carrier],
            h,
            samples(omega, &[carrier], h, |_| {
                (vec![C64::from_re(resistance)], None)
            }),
            PROBE_LIMITS,
            cx,
        )
        .unwrap();
        let mut c = controls(1, h);
        c.output_harmonics = 12;
        c.sideband.maximum_base_residual = 2e-6;
        c.sideband.maximum_omitted_residual = 2e-5;
        c.envelope_tolerances[0] = 2e-5;
        let mut source = PhsImpedanceSidebandSource::new(&mut loaded, &base, &table, &c).unwrap();
        let mut tangent = source.prepare(&scales(4), || false).unwrap();
        let mut factor = tangent.factor(carrier, || false).unwrap();
        let mut forcing = vec![C64::ZERO; 2 * h + 1];
        forcing[h] = C64::new(0.8, -0.2);
        let actual = factor.solve(&forcing, || false).unwrap();
        let mut reference_source =
            PhsPortSidebandSource::forced(&reference, &input, &base.state, &c).unwrap();
        let mut reference_tangent = reference_source.prepare(&scales(3), || false).unwrap();
        let mut reference_factor = reference_tangent.factor(carrier, || false).unwrap();
        let expected = reference_factor.solve(&forcing, || false).unwrap();
        assert_eq!(
            actual.response.state_dimension, 3,
            "native Maxwell coordinate remains present"
        );
        for (&a, &b) in actual
            .response
            .state_response
            .coefficients
            .iter()
            .zip(&expected.state_response.coefficients)
        {
            complex(a, b, 3e-8);
        }
        for (&a, &b) in actual
            .response
            .port_coefficients
            .iter()
            .zip(&expected.port_coefficients)
        {
            complex(a, b, 3e-8);
        }
        assert!(actual.response.state_response.coefficients[3 * h + 2].abs() > 1e-4);
        assert!(
            actual.response.state_response.coefficients[3 * (h + 2)].abs() > 1e-8,
            "nonlinear tangent must create shifted response"
        );
        assert!(
            (actual.response.port_coefficients[c.output_harmonics]
                - actual.response.state_response.coefficients[3 * h])
                .abs()
                > 1e-3,
            "original physical output is the nonlinear storage action, not a state copy"
        );
        assert!(actual.retained_flow_mismatch[0] < 2e-5);
    });
}

#[test]
fn admission_refuses_uncovered_loads_conflicting_base_values_and_false_dc_slopes() {
    with_cx(|cx| {
        let load =
            PeriodicPortImpedance::new(2.0, 1, vec![vec![C64::ONE]; 2], LOAD_LIMITS, cx).unwrap();
        assert!(SidebandPortImpedance::new(&load, vec![0.5], 1, vec![], PROBE_LIMITS, cx).is_err());
        let sample = |w, z, dz| PortImpedanceFrequencySample {
            angular_frequency: w,
            impedance: vec![z],
            angular_derivative: dz,
        };
        assert!(
            SidebandPortImpedance::new(
                &load,
                vec![2.0],
                0,
                vec![sample(2.0, C64::from_re(1.1), None)],
                PROBE_LIMITS,
                cx
            )
            .is_err()
        );
        assert!(
            SidebandPortImpedance::new(
                &load,
                vec![0.0],
                0,
                vec![sample(0.0, C64::new(1.0, 0.01), None)],
                PROBE_LIMITS,
                cx
            )
            .is_err()
        );
        assert!(
            SidebandPortImpedance::new(
                &load,
                vec![0.0],
                0,
                vec![sample(0.0, C64::ONE, Some(vec![C64::from_re(0.01)]))],
                PROBE_LIMITS,
                cx
            )
            .is_err()
        );
        assert!(matches!(
            SidebandPortImpedance::new(
                &load,
                vec![0.5],
                0,
                vec![sample(0.5, C64::from_re(-1.0), None)],
                PROBE_LIMITS,
                cx
            ),
            Err(PortImpedanceError::FrequencyAdmission { .. })
        ));
        let mut tiny = PROBE_LIMITS;
        tiny.maximum_workspace_bytes = 1;
        assert!(matches!(
            SidebandPortImpedance::new(
                &load,
                vec![0.5],
                0,
                vec![sample(0.5, C64::ONE, None)],
                tiny,
                cx
            ),
            Err(PortImpedanceError::Budget(_))
        ));
    });
}

#[test]
fn loaded_base_is_rechecked_and_missing_derivative_or_cancellation_never_falls_back() {
    with_cx(|cx| {
        let system = fs_phs::mass_spring_damper(2.0, 5.0, 0.2).unwrap();
        let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO]]).unwrap();
        let load =
            PeriodicPortImpedance::new(2.0, 1, vec![vec![C64::ONE]; 2], LOAD_LIMITS, cx).unwrap();
        let mut loaded = load.prepare(&system, &input, cx).unwrap();
        let base = loaded
            .solve_scaled_seeded(vec![0.0; 9], &hb(1), &hb_scales(3, 2.0))
            .unwrap();
        let table = SidebandPortImpedance::new(
            &load,
            vec![0.5],
            0,
            samples(2.0, &[0.5], 0, |_| (vec![C64::ONE], None)),
            PROBE_LIMITS,
            cx,
        )
        .unwrap();
        {
            let mut bad = base.clone();
            bad.flow_coefficients[0][0] = C64::from_re(0.1);
            let mut source =
                PhsImpedanceSidebandSource::new(&mut loaded, &bad, &table, &controls(1, 0))
                    .unwrap();
            assert!(matches!(
                source.prepare(&scales(3), || false),
                Err(PortImpedanceError::Sideband(
                    PhsPortSidebandError::Sideband(SidebandError::BaseUnresolved { .. })
                ))
            ));
        }
        let mut source =
            PhsImpedanceSidebandSource::new(&mut loaded, &base, &table, &controls(1, 0)).unwrap();
        let mut tangent = source.prepare(&scales(3), || false).unwrap();
        assert!(tangent.factor(0.6, || false).is_err());
        let mut factor = tangent.factor(0.5, || false).unwrap();
        assert!(matches!(
            factor.solve_with_carrier_derivative(
                &[C64::ONE],
                None,
                derivative_controls(),
                &[1e-8],
                || false
            ),
            Err(PortImpedanceError::Invalid(_))
        ));
        assert!(matches!(
            factor.solve(&[C64::ONE], || true),
            Err(PortImpedanceError::Cancelled)
        ));
        assert!(factor.solve(&[C64::ONE], || false).is_ok());
    });
}
