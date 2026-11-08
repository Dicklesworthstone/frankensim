//! G1/G3/G4: complete simultaneous-tone motion, cross heat and affine sources.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::time::ModalToneBudget;
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolidControls, CompleteSpectralSolidError, MultiToneDampedSpectralSolid,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

const SHAPE: [[f64; 2]; 2] = [[1.0, 0.5], [0.25, 1.0]];
const INVERSE: [[f64; 2]; 2] = [[8.0 / 7.0, -4.0 / 7.0], [-2.0 / 7.0, 8.0 / 7.0]];
const DUAL: [[f64; 2]; 2] = [[8.0 / 7.0, -2.0 / 7.0], [-4.0 / 7.0, 8.0 / 7.0]];

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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} +/- {tolerance:e}"
    );
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
fn coupled(lambda: &[f64], rates: &[f64]) -> (TetElasticAssembly, Csr, Vec<f64>) {
    let n = lambda.len();
    let (mut k, mut m, mut c, mut basis) = (
        vec![0.0; n * n],
        vec![0.0; n * n],
        vec![0.0; n * n],
        vec![0.0; n * n],
    );
    for first in (0..n).step_by(2) {
        for a in 0..2 {
            for b in 0..2 {
                basis[(first + a) * n + first + b] = SHAPE[a][b];
                for j in 0..2 {
                    let product = INVERSE[j][a] * INVERSE[j][b];
                    m[(first + a) * n + first + b] += product;
                    k[(first + a) * n + first + b] += lambda[first + j] * product;
                    c[(first + a) * n + first + b] += rates[first + j] * product;
                }
            }
        }
    }
    (pencil(&k, &m, n), matrix(&c, n), basis)
}
fn mix(values: &[f64], width: usize, map: [[f64; 2]; 2]) -> Vec<f64> {
    let mut out = vec![0.0; values.len()];
    for row in 0..values.len() / width {
        for j in 0..width {
            out[row * width + j] = map[row % 2][0] * values[(row / 2 * 2) * width + j]
                + map[row % 2][1] * values[(row / 2 * 2 + 1) * width + j];
        }
    }
    out
}
fn bits(s: &MultiToneDampedSpectralSolid) -> Vec<u64> {
    let s = s.solid();
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
fn trig_integral(w: f64, start: f64, end: f64, sine: bool) -> f64 {
    let z = 0.5 * w * (end - start);
    let sinc = if z.abs() < 1e-5 {
        1.0 - z * z / 6.0 + z * z * z * z / 120.0
    } else {
        z.sin() / z
    };
    let phase = w * (0.5 * start + 0.5 * end);
    (end - start) * sinc * if sine { phase.sin() } else { phase.cos() }
}
fn sin_sin(a: f64, b: f64, start: f64, end: f64) -> f64 {
    0.5 * (trig_integral(a - b, start, end, false) - trig_integral(a + b, start, end, false))
}
fn cos_sin(a: f64, b: f64, start: f64, end: f64) -> f64 {
    0.5 * (trig_integral(b + a, start, end, true) + trig_integral(b - a, start, end, true))
}

#[test]
fn all_128_modes_retain_cross_frequency_heat_and_spatial_force_phase() {
    with_cx(false, |cx| {
        let n = 128;
        let frequencies = [1.7, 1.9, 11.3];
        let k = frequencies.len();
        let omega: Vec<_> = (0..n)
            .map(|i| {
                if i == 0 {
                    0.0
                } else {
                    2.0 + 0.75 * (i % 7) as f64
                }
            })
            .collect();
        let lambda: Vec<_> = omega.iter().map(|w| w * w).collect();
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
        let (assembly, viscosity, basis) = coupled(&lambda, &rates);
        let a: Vec<_> = (0..n * k)
            .map(|p| 0.001 * (1 + (p / k + 2 * (p % k)) % 5) as f64 / (p % k + 1) as f64)
            .collect();
        let q0: Vec<_> = (0..n).map(|i| (0..k).map(|j| a[i * k + j]).sum()).collect();
        let mut fc = vec![0.0; n * k];
        let mut fs = vec![0.0; n * k];
        for i in 0..n {
            for j in 0..k {
                fc[i * k + j] = a[i * k + j] * (lambda[i] - frequencies[j] * frequencies[j]);
                fs[i * k + j] = -rates[i] * frequencies[j] * a[i * k + j];
            }
        }
        let (mut cosine, mut sine) = (mix(&fc, k, DUAL), mix(&fs, k, DUAL));
        let (initial_cosine, initial_sine) = (cosine.clone(), sine.clone());
        let mut solid = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &viscosity,
            &lambda,
            &basis,
            &frequencies,
            0.37,
            &mix(&q0, 1, SHAPE),
            &vec![0.0; n],
            false,
            controls(),
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let mut cross_signal = 0.0_f64;
        for interval in 0..12 {
            let start = interval as f64 * 0.37;
            let end = (interval + 1) as f64 * 0.37;
            let frame = solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
            let mut q = vec![0.0; n];
            let mut v = vec![0.0; n];
            let (mut work, mut heat, mut isolated_heat) = (0.0, 0.0, 0.0);
            for i in 0..n {
                for j in 0..k {
                    q[i] += a[i * k + j] * (frequencies[j] * end).cos();
                    v[i] -= a[i * k + j] * frequencies[j] * (frequencies[j] * end).sin();
                    for l in 0..k {
                        let vv = a[i * k + l] * frequencies[l];
                        work -= vv
                            * (fc[i * k + j] * cos_sin(frequencies[j], frequencies[l], start, end)
                                + fs[i * k + j]
                                    * sin_sin(frequencies[j], frequencies[l], start, end));
                        let loss = rates[i]
                            * a[i * k + j]
                            * frequencies[j]
                            * vv
                            * sin_sin(frequencies[j], frequencies[l], start, end);
                        heat += loss;
                        if j == l {
                            isolated_heat += loss;
                        }
                    }
                    let (sn, cs) = (frequencies[j] * end).sin_cos();
                    let p = i * k + j;
                    near(
                        cosine[p],
                        initial_cosine[p] * cs + initial_sine[p] * sn,
                        2e-13,
                    );
                    near(
                        sine[p],
                        initial_sine[p] * cs - initial_cosine[p] * sn,
                        2e-13,
                    );
                }
            }
            for (actual, expected) in solid.solid().displacement_m().iter().zip(mix(&q, 1, SHAPE)) {
                near(*actual, expected, 2e-12);
            }
            for (actual, expected) in solid.solid().velocity_mps().iter().zip(mix(&v, 1, SHAPE)) {
                near(*actual, expected, 4e-12);
            }
            near(frame.energy.supplied_work_j, work, 2e-11);
            near(frame.energy.dissipated_j, heat, 2e-11);
            cross_signal = cross_signal.max((heat - isolated_heat).abs());
        }
        assert!(
            cross_signal > 1e-6,
            "fixture must distinguish combined heat from isolated-tone heat"
        );
        assert_eq!(solid.solid().free_dofs().len(), 128);
    });
}

#[test]
fn affine_quadratures_follow_the_manufactured_original_coordinate_trajectory() {
    with_cx(false, |cx| {
        let lambda = [0.0, 9.0];
        let rates = [0.5, 6.0];
        let frequencies = [0.0, 2.0, 5.0];
        let a = [0.1, 0.03, -0.02, 0.04, 0.07, 0.01];
        let b = [0.02, -0.01, 0.005, -0.03, 0.02, -0.006];
        let (assembly, viscosity, basis) = coupled(&lambda, &rates);
        let state = |t: f64| {
            let (mut q, mut v) = (vec![0.0; 2], vec![0.0; 2]);
            for i in 0..2 {
                for j in 0..3 {
                    let (sn, cs) = (frequencies[j] * t).sin_cos();
                    q[i] += (a[3 * i + j] + b[3 * i + j] * t) * cs;
                    v[i] +=
                        b[3 * i + j] * cs - (a[3 * i + j] + b[3 * i + j] * t) * frequencies[j] * sn;
                }
            }
            (mix(&q, 1, SHAPE), mix(&v, 1, SHAPE))
        };
        // t sets the envelope; origin sets the COMMON carrier basis of its endpoints.
        let source = |t: f64, origin: f64| {
            let (mut c, mut s) = (vec![0.0; 6], vec![0.0; 6]);
            for i in 0..2 {
                for j in 0..3 {
                    let amplitude = a[3 * i + j] + b[3 * i + j] * t;
                    let fc = (lambda[i] - frequencies[j] * frequencies[j]) * amplitude
                        + rates[i] * b[3 * i + j];
                    let fs = -frequencies[j] * (2.0 * b[3 * i + j] + rates[i] * amplitude);
                    let (sn, cs) = (frequencies[j] * origin).sin_cos();
                    c[3 * i + j] = fc * cs + fs * sn;
                    s[3 * i + j] = fs * cs - fc * sn;
                }
            }
            (mix(&c, 3, DUAL), mix(&s, 3, DUAL))
        };
        let (q0, v0) = state(0.0);
        let mut solid = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &viscosity,
            &lambda,
            &basis,
            &frequencies,
            0.23,
            &q0,
            &v0,
            true,
            controls(),
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let (mut c, mut s) = source(0.0, 0.0);
        for interval in 0..8 {
            let start = interval as f64 * 0.23;
            let end = (interval + 1) as f64 * 0.23;
            let (ce, se) = source(end, start);
            let frame = solid
                .advance_affine_continuous(&mut c, &mut s, &ce, &se, cx)
                .unwrap();
            let (q, v) = state(end);
            let (expected_c, expected_s) = source(end, end);
            for i in 0..2 {
                near(solid.solid().displacement_m()[i], q[i], 3e-13);
                near(solid.solid().velocity_mps()[i], v[i], 5e-13);
            }
            for i in 0..6 {
                near(c[i], expected_c[i], 3e-13);
                near(s[i], expected_s[i], 3e-13);
            }
            // Integrate original physical F dot v and v^T C v independently.
            let (mut work, mut heat) = (0.0, 0.0);
            for node in 0..=256 {
                let t = start + (end - start) * node as f64 / 256.0;
                let v = state(t).1;
                let f = source(t, t).0;
                let weight = if node == 0 || node == 256 {
                    1.0
                } else if node % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for i in 0..2 {
                    work +=
                        weight * (0..3).map(|j| f[3 * i + j]).sum::<f64>() * v[i] * (end - start)
                            / 768.0;
                    let (columns, values) = viscosity.row(i);
                    let resistance: f64 = columns
                        .iter()
                        .zip(values)
                        .map(|(&j, &value)| value * v[j])
                        .sum();
                    heat += weight * v[i] * resistance * (end - start) / 768.0;
                }
            }
            near(frame.energy.supplied_work_j, work, 3e-11);
            near(frame.energy.dissipated_j, heat, 3e-11);
        }
    });
}

#[test]
fn duplicate_tones_cancel_and_refusals_preserve_all_source_endpoints() {
    with_cx(false, |cx| {
        let assembly = pencil(&[4.0], &[1.0], 1);
        let c = matrix(&[0.5], 1);
        let mut cancel = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &c,
            &[4.0],
            &[1.0],
            &[3.0, 3.0],
            0.2,
            &[0.0],
            &[0.0],
            false,
            controls(),
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([1.0, -1.0], [0.2, -0.2]);
        let frame = cancel
            .advance_continuous(&mut cosine, &mut sine, cx)
            .unwrap();
        assert_eq!(frame.energy.dissipated_j, 0.0);
        assert_eq!(frame.energy.supplied_work_j, 0.0);
        assert_eq!(cancel.solid().displacement_m(), [0.0]);
        assert_eq!(cancel.solid().velocity_mps(), [0.0]);
        let before = bits(&cancel);
        let (mut huge_c, mut huge_s) = ([f64::MAX, -f64::MAX], [f64::MAX, -f64::MAX]);
        assert!(
            cancel
                .advance_continuous(&mut huge_c, &mut huge_s, cx)
                .is_err()
        );
        assert_eq!(bits(&cancel), before);
        assert_eq!(
            (huge_c, huge_s),
            ([f64::MAX, -f64::MAX], [f64::MAX, -f64::MAX])
        );
        let mut solid = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly,
            &c,
            &[4.0],
            &[1.0],
            &[0.0, 3.0],
            0.2,
            &[0.1],
            &[0.2],
            true,
            controls(),
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let (mut c0, mut s0) = ([0.5, 0.6], [0.0, -0.2]);
        let (c1, s1) = ([0.1, 0.9], [0.0, 0.3]);
        let before = bits(&solid);
        assert!(solid.advance(&c0, &s0, cx).is_err());
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.advance_affine_continuous(&mut c0, &mut s0, &c1, &s1, cancelled),
                Err(CompleteSpectralSolidError::Cancelled)
            ))
        });
        assert_eq!((c0, s0), ([0.5, 0.6], [0.0, -0.2]));
        assert_eq!(bits(&solid), before);
        solid
            .advance_affine_continuous(&mut c0, &mut s0, &c1, &s1, cx)
            .unwrap();
        let heat = solid.solid().cumulative_dissipated_j();
        assert!(heat > 0.0);
        let kick = solid.apply_momentum_impulse(&[0.2], 1e-13, cx).unwrap();
        assert_eq!(kick.energy.cumulative_dissipated_j, heat);
        let (before, phases) = (bits(&solid), (c0, s0));
        assert!(
            solid
                .advance_affine_continuous(&mut c0, &mut s0, &[0.1, f64::NAN], &s1, cx)
                .is_err()
        );
        assert_eq!(bits(&solid), before);
        assert_eq!((c0, s0), phases);
        let mut tight = controls();
        tight.maximum_preparation_bytes = 96 + 4096 + 32 * 2;
        assert!(
            MultiToneDampedSpectralSolid::prepare_with_basis(
                &assembly,
                &c,
                &[4.0],
                &[1.0],
                &[0.0, 3.0],
                0.2,
                &[0.1],
                &[0.2],
                true,
                tight,
                ModalToneBudget::realtime_reference(),
                cx
            )
            .is_err()
        );
    });
}
