//! Physical oracles for the exact-linear, non-iterative SOS integrator.
use fs_nlmodal::sav::{SavConfig, SavError, SavStepper};
use fs_nlmodal::{KcStringParams, SosModalStorage, StressChannel, kirchhoff_carrier_string};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::mass_spring_damper;

fn config(offset_j: f64, count: usize) -> SavConfig {
    SavConfig {
        offset_j,
        max_modes: count,
        max_omega_dt: 2.8,
    }
}

#[test]
fn g1_zero_nonlinearity_is_the_same_exact_damped_forced_trajectory() {
    let omegas: Vec<f64> = [1000.0, 5000.0, 10000.0]
        .map(|f| std::f64::consts::TAU * f)
        .into();
    let damping = [0.0, 3.0, 80.0];
    let storage = SosModalStorage {
        omegas,
        channels: vec![],
    };
    let dt = 1.0 / 48000.0;
    let cfg = config(0.01, 3);
    let mut stepper = SavStepper::new(&storage, &damping, cfg, dt).unwrap();
    let initial = [1.0e-5, 0.03, -2.0e-6, 0.05, 7.0e-6, -0.07];
    let mut state = stepper.initial_state(&initial).unwrap();
    let mut exact: Vec<_> = storage
        .omegas
        .iter()
        .zip(damping)
        .map(|(&omega, c)| {
            PreparedExactZoh::new_owned(
                mass_spring_damper(1.0, omega * omega, c).unwrap(),
                dt,
                LinearSplitBudget {
                    max_modes: 2,
                    max_omega_dt: cfg.max_omega_dt,
                },
            )
            .unwrap()
        })
        .collect();
    let mut reference = initial;
    for k in 0..100 {
        let force = [1.0 + (k as f64 * 0.13).sin(), -0.2, 0.4];
        let record = stepper.step(&mut state, &force).unwrap();
        let mut expected_heat = 0.0;
        let mut expected_work = 0.0;
        for i in 0..3 {
            let mut next = [0.0; 2];
            let mut mean = [0.0];
            let r = exact[i]
                .step_into(
                    &reference[2 * i..2 * i + 2],
                    &[force[i]],
                    &mut next,
                    &mut mean,
                )
                .unwrap();
            reference[2 * i..2 * i + 2].copy_from_slice(&next);
            expected_heat += r.dissipated;
            expected_work += dt * force[i] * mean[0];
        }
        for (actual, expected) in state.coordinates().iter().zip(reference) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert_eq!(record.rank_one_solves, 0);
        assert_eq!(record.auxiliary_energy_defect_j, 0.0);
        assert_eq!(record.dissipated_j.to_bits(), expected_heat.to_bits());
        assert!((record.port_work_j - expected_work).abs() < 1e-18);
    }
}

#[test]
fn g1_linear_partials_have_unwarped_phase_at_audio_rate() {
    for frequency in [1000.0, 5000.0, 10000.0] {
        let omega = std::f64::consts::TAU * frequency;
        let storage = SosModalStorage {
            omegas: vec![omega],
            channels: vec![],
        };
        let dt = 1.0 / 48000.0;
        let mut runtime = SavStepper::new(&storage, &[0.0], config(1.0, 1), dt).unwrap();
        let mut q = runtime.initial_state(&[1.0, 0.0]).unwrap();
        let mut p = runtime.initial_state(&[0.0, 1.0]).unwrap();
        runtime.step(&mut q, &[0.0]).unwrap();
        runtime.step(&mut p, &[0.0]).unwrap();
        let phase = (0.5 * (q.coordinates()[0] + p.coordinates()[1])).acos();
        let cents = 1200.0 * (phase / (omega * dt)).log2();
        assert!(cents.abs() < 1e-8, "{frequency} Hz: {cents} cents");
    }
}

fn coupled_storage() -> SosModalStorage {
    SosModalStorage {
        omegas: vec![7.0, 11.0],
        channels: vec![
            StressChannel {
                coefficient: 600.0,
                coupling: vec![1.0, 0.2, 0.2, -0.3],
            },
            StressChannel {
                coefficient: 80.0,
                coupling: vec![0.0, 0.7, 0.7, 0.4],
            },
        ],
    }
}

// Independent component expansion of the full original force, with no SAV
// state, exact-linear cache, or production Storage::gradient call.
fn rhs(storage: &SosModalStorage, x: &[f64; 4], damping: &[f64; 2], force: &[f64; 2]) -> [f64; 4] {
    let mut out = [
        x[1],
        force[0] - storage.omegas[0] * storage.omegas[0] * x[0] - damping[0] * x[1],
        x[3],
        force[1] - storage.omegas[1] * storage.omegas[1] * x[2] - damping[1] * x[3],
    ];
    for c in &storage.channels {
        let e = &c.coupling;
        let s = e[0] * x[0] * x[0] + (e[1] + e[2]) * x[0] * x[2] + e[3] * x[2] * x[2];
        out[1] -= 0.5 * c.coefficient * s * (2.0 * e[0] * x[0] + (e[1] + e[2]) * x[2]);
        out[3] -= 0.5 * c.coefficient * s * ((e[1] + e[2]) * x[0] + 2.0 * e[3] * x[2]);
    }
    out
}

fn rk4(storage: &SosModalStorage, mut x: [f64; 4], duration: f64, steps: usize) -> ([f64; 4], f64) {
    let h = duration / steps as f64;
    let mut maximum_second_energy: f64 = 0.0;
    for _ in 0..steps {
        let f = |x: &[f64; 4]| rhs(storage, x, &[0.0; 2], &[0.0; 2]);
        let add = |a: [f64; 4], b: [f64; 4], c: f64| std::array::from_fn(|i| a[i] + c * b[i]);
        let k1 = f(&x);
        let k2 = f(&add(x, k1, h * 0.5));
        let k3 = f(&add(x, k2, h * 0.5));
        let k4 = f(&add(x, k3, h));
        x = std::array::from_fn(|i| x[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
        maximum_second_energy =
            maximum_second_energy.max(0.5 * (x[3] * x[3] + 121.0 * x[2] * x[2]));
    }
    (x, maximum_second_energy)
}

#[test]
fn g1_signed_multichannel_transfer_and_auxiliary_defect_converge_to_physical_motion() {
    let storage = coupled_storage();
    let initial = [0.2, 0.0, 0.0, 0.0];
    let duration = 2.0;
    let (reference, transferred) = rk4(&storage, initial, duration, 65536);
    assert!(
        transferred > 0.01,
        "nonlinear channel must excite the initially silent mode"
    );
    let mut errors = Vec::new();
    let mut defects = Vec::new();
    for steps in [100, 200, 400, 800] {
        let mut runtime =
            SavStepper::new(&storage, &[0.0; 2], config(1.0, 2), duration / steps as f64).unwrap();
        let mut state = runtime.initial_state(&initial).unwrap();
        let mut max_defect: f64 = 0.0;
        let mut initial_modified = None;
        for _ in 0..steps {
            let record = runtime.step(&mut state, &[0.0; 2]).unwrap();
            assert_eq!(record.rank_one_solves, 1);
            let h0 = *initial_modified.get_or_insert(record.modified_energy_j);
            assert!((record.modified_energy_j - h0).abs() < 2e-11 * h0);
            max_defect = max_defect.max(record.auxiliary_energy_defect_j.abs());
        }
        let x = state.coordinates();
        let error = ((7.0 * (x[0] - reference[0])).powi(2)
            + (x[1] - reference[1]).powi(2)
            + (11.0 * (x[2] - reference[2])).powi(2)
            + (x[3] - reference[3]).powi(2))
        .sqrt();
        errors.push(error);
        defects.push(max_defect);
    }
    eprintln!("multichannel state errors={errors:?}; maximum physical-energy defects={defects:?}");
    for values in [&errors, &defects] {
        for pair in values.windows(2).skip(1) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (1.8..2.25).contains(&order),
                "second-order physical convergence: {values:?}, p={order}"
            );
        }
    }
    assert!(errors[3] < 5e-4);
    assert!(defects[3] < 1e-4);
}

fn elliptic_duffing_frequency(omega: f64, ratio: f64) -> f64 {
    let m = ratio / (2.0 * (1.0 + ratio));
    let (mut a, mut b) = (1.0, (1.0 - m).sqrt());
    for _ in 0..20 {
        let next = 0.5 * (a + b);
        b = (a * b).sqrt();
        a = next;
    }
    omega * (1.0 + ratio).sqrt() * a / std::f64::consts::TAU
}

// Cubic Hermite interpolation uses the physical endpoint velocities. Its
// crossing error is higher order than the stepper and does not mask its ladder.
fn crossing_fraction(a: &[f64], b: &[f64], h: f64) -> f64 {
    let (mut low, mut high) = (0.0, 1.0);
    for _ in 0..48 {
        let s = 0.5 * (low + high);
        let s2 = s * s;
        let s3 = s2 * s;
        let q = (2.0 * s3 - 3.0 * s2 + 1.0) * a[0]
            + (s3 - 2.0 * s2 + s) * h * a[1]
            + (-2.0 * s3 + 3.0 * s2) * b[0]
            + (s3 - s2) * h * b[1];
        if q > 0.0 {
            low = s;
        } else {
            high = s;
        }
    }
    0.5 * (low + high)
}

#[test]
fn g1_5khz_duffing_backbone_refines_to_the_independent_elliptic_frequency() {
    let omega = std::f64::consts::TAU * 5000.0;
    let amplitude = 1e-5;
    for ratio in [0.01, 0.1, 0.5] {
        let storage = SosModalStorage {
            omegas: vec![omega],
            channels: vec![StressChannel {
                coefficient: ratio * omega * omega / (amplitude * amplitude),
                coupling: vec![1.0],
            }],
        };
        let reference = elliptic_duffing_frequency(omega, ratio);
        let mut errors = Vec::new();
        for substeps in [1, 2, 4, 8, 16, 32, 64, 128] {
            let dt = 1.0 / (48000.0 * substeps as f64);
            let mut runtime = SavStepper::new(&storage, &[0.0], config(0.05, 1), dt).unwrap();
            let mut state = runtime.initial_state(&[amplitude, 0.0]).unwrap();
            let mut crossings = Vec::new();
            for k in 0..(1000 * substeps) {
                let previous = [state.coordinates()[0], state.coordinates()[1]];
                runtime.step(&mut state, &[0.0]).unwrap();
                if previous[0] > 0.0 && state.coordinates()[0] <= 0.0 {
                    crossings.push(
                        (k as f64 + crossing_fraction(&previous, state.coordinates(), dt)) * dt,
                    );
                    if crossings.len() == 65 {
                        break;
                    }
                }
            }
            assert_eq!(crossings.len(), 65);
            let measured = 64.0 / (crossings[64] - crossings[0]);
            errors.push((1200.0 * (measured / reference).log2()).abs());
        }
        eprintln!("Duffing ratio={ratio}: cents errors={errors:?}");
        // Strong nonlinearity has substantial higher-order error on the coarse
        // rungs. Measure formal order on three resolved rungs, without accepting
        // a coarse result merely because its modified energy is conserved.
        for pair in errors.windows(2).skip(5) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (1.8..2.2).contains(&order),
                "Duffing temporal order {order}: {errors:?}"
            );
        }
        assert!(errors[4] < 0.1, "16-substep physical pitch: {errors:?}");
    }
}

#[test]
fn g1_held_force_keeps_the_original_static_nonlinear_compliance() {
    let storage = coupled_storage();
    let equilibrium = [0.13, 0.0, -0.07, 0.0];
    let restoring = rhs(&storage, &equilibrium, &[0.0; 2], &[0.0; 2]);
    let force = [-restoring[1], -restoring[3]];
    for dt in [0.05, 0.025, 0.0125] {
        let mut runtime = SavStepper::new(&storage, &[0.2, 0.5], config(0.2, 2), dt).unwrap();
        let mut state = runtime.initial_state(&equilibrium).unwrap();
        for _ in 0..100 {
            let r = runtime.step(&mut state, &force).unwrap();
            assert!(r.dissipated_j < 1e-24);
        }
        for (actual, expected) in state.coordinates().iter().zip(equilibrium) {
            assert!(
                (actual - expected).abs() < 2e-13,
                "static compliance changed: {:?}",
                state.coordinates()
            );
        }
    }
}

#[test]
fn g0_forced_damped_modified_energy_ledger_and_checkpoint_replay() {
    let storage = coupled_storage();
    let mut runtime = SavStepper::new(&storage, &[0.2, 0.3], config(1.0, 2), 0.001).unwrap();
    let mut state = runtime.initial_state(&[0.2, 0.1, -0.04, 0.05]).unwrap();
    let mut checkpoint = state.clone();
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut h0 = None;
    let mut last = 0.0;
    for k in 0..2000 {
        let force = [0.15 * (k as f64 * 0.01).sin(), -0.03];
        let record = runtime.step(&mut state, &force).unwrap();
        h0.get_or_insert(record.modified_energy_j - record.modified_energy_change_j);
        assert!(record.dissipated_j >= 0.0);
        assert_eq!(record.rank_one_solves, 1);
        work += record.port_work_j;
        heat += record.dissipated_j;
        last = record.modified_energy_j;
    }
    assert!((last - h0.unwrap() + heat - work).abs() < 2e-10);
    assert!(heat > 0.1);
    for k in 0..2000 {
        runtime
            .step(&mut checkpoint, &[0.15 * (k as f64 * 0.01).sin(), -0.03])
            .unwrap();
    }
    assert_eq!(state, checkpoint);
}

#[test]
fn g0_weak_nonlinear_energy_is_not_erased_by_a_large_auxiliary_offset() {
    let storage = SosModalStorage {
        omegas: vec![100.0],
        channels: vec![StressChannel {
            coefficient: 1.0e12,
            coupling: vec![1.0],
        }],
    };
    let mut runtime = SavStepper::new(&storage, &[0.0], config(1.0e8, 1), 1.0e-4).unwrap();
    let mut state = runtime.initial_state(&[1.0e-9, 0.0]).unwrap();
    let record = runtime.step(&mut state, &[0.0]).unwrap();
    // V0=2.5e-25 J is far below one ulp of C0. Storing r alone would erase
    // the entire nonlinear energy; the shifted auxiliary retains its update.
    assert!(record.auxiliary_energy_defect_j.abs() < 1.0e-27);
    assert!(record.physical_energy_j > 0.0);
    assert_eq!(record.rank_one_solves, 1);
}

#[test]
fn g0_kirchhoff_carrier_retains_128_modes_with_one_scalar_solve() {
    let storage = kirchhoff_carrier_string(
        &KcStringParams {
            length: 1.0,
            tension: 100.0,
            lin_density: 0.01,
            ea: 20000.0,
        },
        128,
    )
    .unwrap();
    let mut runtime =
        SavStepper::new(&storage, &[0.1; 128], config(0.1, 128), 1.0 / 48000.0).unwrap();
    let mut initial = vec![0.0; 256];
    for k in 0..128 {
        initial[2 * k] = 1e-5 / ((k + 1) * (k + 1)) as f64;
    }
    let mut state = runtime.initial_state(&initial).unwrap();
    let record = runtime.step(&mut state, &[0.0; 128]).unwrap();
    assert_eq!(state.coordinates().len(), 256);
    assert_eq!(record.rank_one_solves, 1);
    assert!(record.physical_energy_j > 0.0);
    assert!(
        state
            .coordinates()
            .chunks_exact(2)
            .all(|pair| pair[0] != 0.0)
    );
}

#[test]
fn g0_invalid_bands_capacities_offsets_and_cancelled_steps_are_transactional() {
    let storage = coupled_storage();
    for offset in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(SavStepper::new(&storage, &[0.0; 2], config(offset, 2), 0.01).is_err());
    }
    assert!(SavStepper::new(&storage, &[0.0; 2], config(1.0, 1), 0.01).is_err());
    let cfg = SavConfig {
        max_omega_dt: 0.11,
        ..config(1.0, 2)
    };
    assert!(SavStepper::new(&storage, &[0.0; 2], cfg, 0.01).is_ok());
    let above = f64::from_bits(0.01_f64.to_bits() + 1);
    assert!(matches!(
        SavStepper::new(&storage, &[0.0; 2], cfg, above),
        Err(SavError::UnresolvedMode { mode: 1, .. })
    ));
    assert!(
        SavStepper::new(
            &storage,
            &[0.0; 2],
            SavConfig {
                max_omega_dt: std::f64::consts::PI,
                ..cfg
            },
            0.01
        )
        .is_err()
    );
    for damping in [[-1.0, 0.0], [f64::NAN, 0.0], [0.0, f64::INFINITY]] {
        assert!(SavStepper::new(&storage, &damping, config(1.0, 2), 0.01).is_err());
    }
    let mut runtime = SavStepper::new(&storage, &[0.0; 2], config(1.0, 2), 0.01).unwrap();
    let mut state = runtime.initial_state(&[0.2, 0.0, 0.03, 0.0]).unwrap();
    let unchanged = state.clone();
    assert!(runtime.step(&mut state, &[f64::NAN, 0.0]).is_err());
    assert_eq!(state, unchanged);
    let mut calls = 0;
    let result = runtime.step_controlled(&mut state, &[0.1, 0.0], || {
        calls += 1;
        calls == 7
    });
    assert!(matches!(result, Err(SavError::Cancelled)));
    assert_eq!(state, unchanged);
    assert!(runtime.step(&mut state, &[f64::MAX; 2]).is_err());
    assert_eq!(state, unchanged);
    let mut retry = unchanged;
    runtime.step(&mut state, &[0.1, 0.0]).unwrap();
    runtime.step(&mut retry, &[0.1, 0.0]).unwrap();
    assert_eq!(state, retry);
}
