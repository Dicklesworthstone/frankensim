//! G1/G3/G4: continuous all-mode forcing, physical power and source-phase rollback.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolidControls, CompleteSpectralSolidError, DampedSpectralSolid,
    HarmonicDampedSpectralSolid,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 35,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 128,
        maximum_preparation_bytes: 8 * 1024 * 1024,
        maximum_preparation_multiply_terms: 100_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-10,
        maximum_steps: 1000,
    }
}
fn matrix(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            coo.push(i, j, values[i * n + j]);
        }
    }
    coo.assemble()
}
fn pencil(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    TetElasticAssembly {
        stiffness: matrix(k, n),
        mass: matrix(m, n),
        free_dofs: (0..n).map(|i| 2 * i + 1).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} +/- {tolerance:e}"
    );
}
fn bits(s: &HarmonicDampedSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([
            s.time_s(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
            s.cumulative_work_j(),
            s.cumulative_dissipated_j(),
        ])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn all_128_coupled_modes_keep_physical_phase_and_independent_power() {
    with_cx(false, |cx| {
        let n = 128;
        let shape = [[1.0, 0.5], [0.25, 1.0]];
        let inverse = [[8.0 / 7.0, -4.0 / 7.0], [-2.0 / 7.0, 8.0 / 7.0]];
        let omega: Vec<_> = (0..n)
            .map(|i| {
                if i == 0 {
                    0.0
                } else {
                    2.0 + 0.75 * (i % 7) as f64
                }
            })
            .collect();
        let rates: Vec<_> = omega
            .iter()
            .enumerate()
            .map(|(i, w)| {
                if i == 0 {
                    0.8
                } else {
                    [0.1, 2.0, 4.0][i % 3] * w
                }
            })
            .collect();
        let lambda: Vec<_> = omega.iter().map(|w| w * w).collect();
        let (mut k, mut m, mut c, mut basis) = (
            vec![0.0; n * n],
            vec![0.0; n * n],
            vec![0.0; n * n],
            vec![0.0; n * n],
        );
        for first in (0..n).step_by(2) {
            for a in 0..2 {
                for b in 0..2 {
                    basis[(first + a) * n + first + b] = shape[a][b];
                    for j in 0..2 {
                        let product = inverse[j][a] * inverse[j][b];
                        m[(first + a) * n + first + b] += product;
                        k[(first + a) * n + first + b] += lambda[first + j] * product;
                        c[(first + a) * n + first + b] += rates[first + j] * product;
                    }
                }
            }
        }
        let transform = |values: &[f64], map: [[f64; 2]; 2]| -> Vec<f64> {
            (0..n)
                .map(|i| map[i % 2][0] * values[i / 2 * 2] + map[i % 2][1] * values[i / 2 * 2 + 1])
                .collect()
        };
        let dual = [
            [inverse[0][0], inverse[1][0]],
            [inverse[0][1], inverse[1][1]],
        ];
        let fc: Vec<_> = (0..n).map(|i| 0.03 - 0.001 * (i % 10) as f64).collect();
        let fs: Vec<_> = (0..n).map(|i| 0.01 * (i % 3) as f64 - 0.01).collect();
        let drive = 3.7;
        let (mut a, mut b) = (vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            let detuning = lambda[i] - drive * drive;
            let loss = rates[i] * drive;
            let denominator = detuning * detuning + loss * loss;
            a[i] = (detuning * fc[i] - loss * fs[i]) / denominator;
            b[i] = (loss * fc[i] + detuning * fs[i]) / denominator;
        }
        let initial_v: Vec<_> = b.iter().map(|value| drive * value).collect();
        let mut solid = HarmonicDampedSpectralSolid::prepare_with_basis(
            &pencil(&k, &m, n),
            &matrix(&c, n),
            &lambda,
            &basis,
            drive,
            0.23,
            &transform(&a, shape),
            &transform(&initial_v, shape),
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = (transform(&fc, dual), transform(&fs, dual));
        let (initial_cosine, initial_sine) = (cosine.clone(), sine.clone());
        for interval in 0..20 {
            let start = interval as f64 * 0.23;
            let end = (interval + 1) as f64 * 0.23;
            let frame = solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
            let (sin, cos) = (drive * end).sin_cos();
            let q: Vec<_> = (0..n).map(|i| a[i] * cos + b[i] * sin).collect();
            let v: Vec<_> = (0..n).map(|i| drive * (b[i] * cos - a[i] * sin)).collect();
            for (actual, expected) in solid.displacement_m().iter().zip(transform(&q, shape)) {
                near(*actual, expected, 1e-12);
            }
            for (actual, expected) in solid.velocity_mps().iter().zip(transform(&v, shape)) {
                near(*actual, expected, 3e-12);
            }
            // Integrals of cos², sin² and cos*sin give independent physical power.
            let cc = 0.5 * (end - start)
                + ((2.0 * drive * end).sin() - (2.0 * drive * start).sin()) / (4.0 * drive);
            let ss = end - start - cc;
            let cs = ((2.0 * drive * start).cos() - (2.0 * drive * end).cos()) / (4.0 * drive);
            let (mut work, mut heat) = (0.0, 0.0);
            for i in 0..n {
                let (vc, vs) = (drive * b[i], -drive * a[i]);
                work += fc[i] * vc * cc + fs[i] * vs * ss + (fc[i] * vs + fs[i] * vc) * cs;
                heat += rates[i] * (vc * vc * cc + vs * vs * ss + 2.0 * vc * vs * cs);
                near(
                    cosine[i],
                    initial_cosine[i] * cos + initial_sine[i] * sin,
                    1e-14,
                );
                near(
                    sine[i],
                    initial_sine[i] * cos - initial_cosine[i] * sin,
                    1e-14,
                );
            }
            near(frame.energy.supplied_work_j, work, 1e-11);
            near(frame.energy.dissipated_j, heat, 1e-11);
            assert!(frame.energy.dissipated_j > 0.0);
            assert!(frame.energy.cumulative_balance_error_j.abs() < 1e-11);
        }
        assert_eq!(solid.free_dofs().len(), 128);
    });
}

#[test]
fn undamped_resonance_grows_continuously_and_near_resonance_is_regular() {
    with_cx(false, |cx| {
        let assembly = pencil(&[16.0], &[4.0], 1);
        let mut solid = HarmonicDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &matrix(&[0.0], 1),
            &[4.0],
            &[0.5],
            2.0,
            0.13,
            &[0.0],
            &[0.0],
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([4.0], [0.0]);
        for interval in 0..60 {
            let start = interval as f64 * 0.13;
            let end = (interval + 1) as f64 * 0.13;
            let frame = solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
            near(
                solid.displacement_m()[0],
                0.25 * end * (2.0 * end).sin(),
                3e-12,
            );
            near(
                solid.velocity_mps()[0],
                0.25 * (2.0 * end).sin() + 0.5 * end * (2.0 * end).cos(),
                3e-12,
            );
            // Primitive of 4*cos(2t)*v(t), independently of endpoint energy.
            let work = 0.5 * (end * end - start * start)
                + 0.25 * (end * (4.0 * end).sin() - start * (4.0 * start).sin())
                - ((4.0 * end).cos() - (4.0 * start).cos()) / 16.0;
            near(frame.energy.supplied_work_j, work, 2e-11);
            assert_eq!(frame.energy.dissipated_j, 0.0);
        }
        for frequency in [2.0 - 1e-12, 2.0 + 1e-12] {
            let mut near_resonance = HarmonicDampedSpectralSolid::prepare_with_basis(
                &assembly,
                &matrix(&[0.0], 1),
                &[4.0],
                &[0.5],
                frequency,
                0.37,
                &[0.0],
                &[0.0],
                controls(),
                cx,
            )
            .unwrap();
            near_resonance.advance(&[4.0], &[0.0], cx).unwrap();
            near(
                near_resonance.displacement_m()[0],
                0.25 * 0.37 * (2.0_f64 * 0.37).sin(),
                1e-12,
            );
        }
    });
}

#[test]
fn damped_resonance_keeps_its_initial_transient_and_viscous_heat() {
    with_cx(false, |cx| {
        let response = |t: f64| {
            let damped = 3.99_f64.sqrt();
            let decay = (-0.1 * t).exp();
            let q = 2.5 * (2.0 * t).sin() - 5.0 / damped * decay * (damped * t).sin();
            let v = 5.0 * (2.0 * t).cos()
                - 5.0 * decay * ((damped * t).cos() - 0.1 / damped * (damped * t).sin());
            (q, v)
        };
        let mut solid = HarmonicDampedSpectralSolid::prepare_with_basis(
            &pencil(&[16.0], &[4.0], 1),
            &matrix(&[0.8], 1),
            &[4.0],
            &[0.5],
            2.0,
            0.11,
            &[0.0],
            &[0.0],
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([4.0], [0.0]);
        for interval in 0..25 {
            let start = interval as f64 * 0.11;
            let end = (interval + 1) as f64 * 0.11;
            let frame = solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
            let expected = response(end);
            near(solid.displacement_m()[0], expected.0, 2e-13);
            near(solid.velocity_mps()[0], expected.1, 4e-13);
            let (mut work, mut heat) = (0.0, 0.0);
            for node in 0..=256 {
                let t = start + (end - start) * node as f64 / 256.0;
                let v = response(t).1;
                let weight = if node == 0 || node == 256 {
                    1.0
                } else if node % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                work += weight * 4.0 * (2.0 * t).cos() * v * (end - start) / 768.0;
                heat += weight * 0.8 * v * v * (end - start) / 768.0;
            }
            near(frame.energy.supplied_work_j, work, 2e-12);
            near(frame.energy.dissipated_j, heat, 2e-12);
        }
    });
}

#[test]
fn dc_rigid_drag_matches_held_arithmetic_and_preserves_source_bits() {
    with_cx(false, |cx| {
        let assembly = pencil(&[0.0], &[4.0], 1);
        let c = matrix(&[12.0], 1);
        let mut plain = DampedSpectralSolid::prepare_with_basis(
            &assembly,
            &c,
            &[0.0],
            &[0.5],
            0.2,
            &[0.7],
            &[0.2],
            controls(),
            cx,
        )
        .unwrap();
        let mut harmonic = HarmonicDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &c,
            &[0.0],
            &[0.5],
            0.0,
            0.2,
            &[0.7],
            &[0.2],
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([8.0], [-17.0]);
        for _ in 0..10 {
            let expected = plain.advance(&cosine, cx).unwrap();
            let actual = harmonic
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
            assert_eq!(actual.energy, expected.energy);
            assert_eq!(
                harmonic.displacement_m()[0].to_bits(),
                plain.displacement_m()[0].to_bits()
            );
            assert_eq!(
                harmonic.velocity_mps()[0].to_bits(),
                plain.velocity_mps()[0].to_bits()
            );
            assert_eq!(cosine, [8.0]);
            assert_eq!(sine, [-17.0]);
        }
    });
}

#[test]
fn refusal_preserves_motion_and_phase_and_impulses_retain_harmonic_heat() {
    with_cx(false, |cx| {
        let mut limited = controls();
        limited.maximum_steps = 1;
        let mut solid = HarmonicDampedSpectralSolid::prepare_with_basis(
            &pencil(&[4.0], &[1.0], 1),
            &matrix(&[0.5], 1),
            &[4.0],
            &[1.0],
            2.0,
            0.2,
            &[0.1],
            &[0.2],
            limited,
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([0.8], [0.3]);
        let before = bits(&solid);
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.advance_continuous(&mut cosine, &mut sine, cancelled),
                Err(CompleteSpectralSolidError::Cancelled)
            ))
        });
        assert_eq!(bits(&solid), before);
        assert_eq!((cosine, sine), ([0.8], [0.3]));
        let (mut huge_cosine, mut huge_sine) = ([f64::MAX], [f64::MAX]);
        assert!(
            solid
                .advance_continuous(&mut huge_cosine, &mut huge_sine, cx)
                .is_err()
        );
        assert_eq!((huge_cosine, huge_sine), ([f64::MAX], [f64::MAX]));
        assert_eq!(bits(&solid), before);
        solid
            .advance_continuous(&mut cosine, &mut sine, cx)
            .unwrap();
        let heat = solid.cumulative_dissipated_j();
        assert!(heat > 0.0);
        let kick = solid.apply_momentum_impulse(&[0.2], 1e-13, cx).unwrap();
        assert_eq!(kick.energy.cumulative_dissipated_j, heat);
        assert_eq!(solid.steps(), 1);
        let (before, sources) = (bits(&solid), (cosine, sine));
        assert!(
            solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .is_err()
        );
        assert_eq!(bits(&solid), before);
        assert_eq!((cosine, sine), sources);
    });
}
