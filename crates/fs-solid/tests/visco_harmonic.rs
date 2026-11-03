//! G1/G3 dispersive forced resonance on the original tensor-Prony solid.
//! These analytic/manufactured specimens are software checks, not material data.

use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, HarmonicResponseError};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TET: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 873,
                    kernel_id: 73,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 100,
        maximum_payload_bytes: 4 * 1024 * 1024,
        maximum_relative_residual: 2e-12,
    }
}
fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn rank_one(vector: [f64; 6], factor: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| factor * vector[i] * vector[j]))
}
fn solid(law: TensorProny, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.0,
        material_state_identity: [73; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &TET,
        },
        TetViscoField::Uniform(&material),
        &fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-15),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-15),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn scalar_law(loss: bool) -> TensorProny {
    let branch = if loss {
        vec![PronyBranch::new(diagonal([18.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 0.1).unwrap()]
    } else {
        Vec::new()
    };
    let law = TensorProny::new_mandel(diagonal([6.0; 6]), branch, 0.02).unwrap();
    if loss {
        law.with_kelvin_viscosity_mandel(diagonal([1.2, 0.0, 0.0, 0.0, 0.0, 0.0]))
            .unwrap()
    } else {
        law
    }
}

#[test]
fn g1_hereditary_resonance_matches_analytic_complex_modulus_and_heat() {
    with_cx(|cx, _| {
        // Unit tetra: volume=1/6, B_xx=1, M_33=rho*V/10=1 kg.
        // K_eq=1, one branch K_j=3, tau=.1, and Kelvin C=.2.
        let source = solid(scalar_law(true), &[3], cx);
        for w in [0.5, 1.0, 1.1, 2.0, 5.0, 10.0] {
            let x = 0.1 * w;
            let storage = 1.0 + 3.0 * x * x / (1.0 + x * x);
            let loss = 0.2 * w + 3.0 * x / (1.0 + x * x);
            let dynamic = C64::new(storage - w * w, -loss);
            let prepared = source.prepare_harmonic_response(w, controls(), cx).unwrap();
            assert_eq!(prepared.free_dofs(), &[3]);
            assert_eq!(prepared.omega_rad_s(), w);
            for f in [C64::new(0.001, 0.0), C64::new(-0.0004, 0.0007)] {
                let expected = f / dynamic;
                let result = prepared.solve(&[f], cx).unwrap();
                near(result.displacement_m[0], expected);
                near(
                    result.velocity_m_s[0],
                    C64::new(w * expected.im, -w * expected.re),
                );
                let heat = 0.5 * w * loss * expected.norm_sq();
                close(result.average_material_heat_w, heat, 2e-11);
                close(result.average_input_power_w, heat, 2e-11);
                assert!(result.power_balance_residual_w.abs() <= 2e-11 * heat);
                assert!(result.relative_residual < 2e-12);
                let replay = prepared.solve(&[f], cx).unwrap();
                assert_eq!(result.displacement_m, replay.displacement_m);
                assert_eq!(
                    result.average_material_heat_w.to_bits(),
                    replay.average_material_heat_w.to_bits()
                );
            }
        }
    });
}

#[test]
fn g1_full_nonproportional_tensor_and_consistent_mass_preserve_forced_phasors() {
    with_cx(|cx, _| {
        let mut eq = diagonal([6.0, 9.0, 12.0, 15.0, 18.0, 21.0]);
        eq[0][3] = 1.2;
        eq[3][0] = 1.2;
        let law = TensorProny::new_mandel(
            eq,
            vec![PronyBranch::new(rank_one([1.0, 0.0, 0.0, 2.0, 0.0, 0.0], 12.0), 0.3).unwrap()],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(rank_one([1.0, 0.0, 0.0, -1.0, 0.0, 0.0], 0.6))
        .unwrap();
        let source = solid(law, &[3, 6], cx);
        let s = core::f64::consts::FRAC_1_SQRT_2;
        // Two x-displacement rows on distinct nodes retain M_01=.5.
        let mass = [1.0, 0.5, 0.5, 1.0];
        let equilibrium = [1.0, 0.2 * s, 0.2 * s, 1.25];
        let branch = [2.0, 4.0 * s, 4.0 * s, 4.0];
        let viscosity = [0.1, -0.1 * s, -0.1 * s, 0.05];
        let q = [C64::new(0.0004, -0.0002), C64::new(-0.0001, 0.0003)];
        for w in [0.4, 1.5, 4.0, 10.0] {
            let x = 0.3 * w;
            let k: [f64; 4] =
                core::array::from_fn(|i| equilibrium[i] + branch[i] * x * x / (1.0 + x * x));
            let loss: [f64; 4] =
                core::array::from_fn(|i| viscosity[i] * w + branch[i] * x / (1.0 + x * x));
            let force: Vec<_> = (0..2)
                .map(|i| {
                    (0..2).fold(C64::ZERO, |sum, j| {
                        sum + C64::new(k[2 * i + j] - w * w * mass[2 * i + j], -loss[2 * i + j])
                            * q[j]
                    })
                })
                .collect();
            let result = source
                .prepare_harmonic_response(w, controls(), cx)
                .unwrap()
                .solve(&force, cx)
                .unwrap();
            for (&actual, &expected) in result.displacement_m.iter().zip(&q) {
                near(actual, expected);
            }
            let cross = q[0].re * q[1].re + q[0].im * q[1].im;
            let expected_heat = 0.5
                * w
                * (loss[0] * q[0].norm_sq() + 2.0 * loss[1] * cross + loss[3] * q[1].norm_sq());
            close(result.average_material_heat_w, expected_heat, 2e-10);
            close(result.average_input_power_w, expected_heat, 2e-10);
        }
    });
}

#[test]
fn g0_full_cycle_strain_admission_catches_between_phase_violation() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(true), &[3], cx);
        let w = 2.0;
        let x = 0.1 * w;
        let dynamic = C64::new(
            1.0 + 3.0 * x * x / (1.0 + x * x) - w * w,
            -(0.2 * w + 3.0 * x / (1.0 + x * x)),
        );
        let prepared = source.prepare_harmonic_response(w, controls(), cx).unwrap();
        // Each quadrature-phase strain .015 passes the .02 limit separately,
        // while the true cycle maximum sqrt(2)*.015 exceeds it.
        let q = C64::new(0.015, 0.015);
        assert!(matches!(
            prepared.solve(&[dynamic * q], cx),
            Err(ViscoHarmonicError::Element(ViscoElementError::Material(_)))
        ));
        near(
            prepared
                .solve(&[dynamic * q.scale(0.5)], cx)
                .unwrap()
                .displacement_m[0],
            q.scale(0.5),
        );
    });
}

#[test]
fn g0_source_band_resource_and_unrepresentable_frequency_refuse() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(true).with_band((0.5, 10.0)).unwrap(), &[3], cx);
        assert!(matches!(
            source.prepare_harmonic_response(11.0, controls(), cx),
            Err(ViscoHarmonicError::Element(ViscoElementError::Material(_)))
        ));
        for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                source.prepare_harmonic_response(frequency, controls(), cx),
                Err(ViscoHarmonicError::Harmonic(
                    HarmonicResponseError::Invalid(_)
                ))
            ));
        }
        let mut budget = controls();
        budget.maximum_coordinates = 0;
        assert!(matches!(
            source.prepare_harmonic_response(1.0, budget, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        budget = controls();
        // One coordinate: 24 source + 128 matrix + 128 vector bytes.
        budget.maximum_payload_bytes = 279;
        assert!(matches!(
            source.prepare_harmonic_response(1.0, budget, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Budget { requested: 280, .. }
            ))
        ));
        budget.maximum_payload_bytes = 280;
        assert!(source.prepare_harmonic_response(1.0, budget, cx).is_ok());
    });
}

#[test]
fn g0_g4_lossless_response_remains_lossless_and_cancellation_refuses() {
    with_cx(|cx, gate| {
        let source = solid(scalar_law(false), &[3], cx);
        let prepared = source
            .prepare_harmonic_response(2.0, controls(), cx)
            .unwrap();
        let result = prepared.solve(&[C64::new(0.003, 0.0)], cx).unwrap();
        near(result.displacement_m[0], C64::new(-0.001, 0.0));
        assert_eq!(result.average_material_heat_w, 0.0);
        assert_eq!(result.average_input_power_w, 0.0);
        gate.request();
        assert!(matches!(
            prepared.solve(&[C64::new(0.003, 0.0)], cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Cancelled
            ))
        ));
        assert!(
            source
                .prepare_harmonic_response(2.0, controls(), cx)
                .is_err()
        );
    });
}
