//! Independent Duffing-Maxwell evolution and the original film/contact owner.
use super::*;
use fs_couple::render::plate::impact::membrane::{
    MembraneRelaxationBudget, MembraneRelaxationSpectrum,
};
use fs_couple::render::plate::impact::relaxation::InitialMemory;

const RATIO: f64 = 0.7;
const TAU: f64 = 0.00015;
fn budget() -> MembraneRelaxationBudget {
    MembraneRelaxationBudget {
        maximum_source_arms: 8,
        maximum_memory_states: 64,
        maximum_setup_terms: 2_000_000,
    }
}
fn spectrum() -> MembraneRelaxationSpectrum {
    MembraneRelaxationSpectrum {
        branches: vec![(RATIO, TAU)],
        omega_band_rad_s: (0.0, 1e7),
    }
}
fn loaded(p: MembranePotential, amplitude: f64, initial: InitialMemory) -> MembranePotential {
    p.with_proportional_relaxation(spectrum(), vec![2.0 * amplitude], initial, budget())
        .unwrap()
}

// Independent continuous q,v,m,heat system. m is viscous q^2, not the
// production QR coordinate. beta is the ACTUAL film's quartic coefficient.
fn rhs(x: [f64; 4], k: f64, beta: f64) -> [f64; 4] {
    let [q, v, m, _] = x;
    let gap = q * q - m;
    [
        v,
        -k * q - beta * q * q * q - RATIO * beta * q * gap,
        gap / TAU,
        0.5 * RATIO * beta * gap * gap / TAU,
    ]
}
fn continuous(k: f64, beta: f64, q: f64, relaxed: bool, time: f64) -> [f64; 4] {
    let mut x = [q, 0.0, if relaxed { q * q } else { 0.0 }, 0.0];
    let dt = time / 16384.0;
    for _ in 0..16384 {
        let a = rhs(x, k, beta);
        let b = rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]), k, beta);
        let c = rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]), k, beta);
        let d = rhs(core::array::from_fn(|i| x[i] + dt * c[i]), k, beta);
        for i in 0..4 {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn nonlinear_film_memory_matches_independent_dynamics_and_heat_on_three_clocks() {
    let p = reduced(0.25);
    let k = p.reduction().omegas()[0].powi(2);
    let beta = 4.0 * p.reduction().stretching_energy(&[1.0]);
    let amplitude = 0.08 / p.reduction().maximum_slope(&[1.0]);
    let duration = 192.0 * config().dt_s;
    for relaxed in [false, true] {
        let expected = continuous(k, beta, amplitude, relaxed, duration);
        let mut errors = Vec::new();
        for steps in [192, 384, 768] {
            let law = loaded(
                p.clone(),
                amplitude,
                if relaxed {
                    InitialMemory::Relaxed
                } else {
                    InitialMemory::Unrelaxed
                },
            );
            assert_eq!(law.relaxation_state_count(), 1);
            let mut cfg = config();
            cfg.dt_s = duration / f64::from(steps);
            let system =
                ImpactSystem::new(vec![body(law, amplitude, 0.0)], vec![], vec![], vec![], cfg)
                    .unwrap();
            let energy = system.stored_energy_j();
            let equilibrium = 0.5 * k * amplitude * amplitude + 0.25 * beta * amplitude.powi(4);
            let initial = equilibrium
                + if relaxed {
                    0.0
                } else {
                    0.25 * RATIO * beta * amplitude.powi(4)
                };
            assert!((energy - initial).abs() < 2e-12 * initial);
            let mut system = system.prepare_analytic().unwrap();
            let gate = CancelGate::new_clock_free();
            let mut heat = 0.0;
            for _ in 0..steps {
                let sample = system.step(&[0.0], &gate).unwrap();
                assert!(sample.dissipated_energy_j >= 0.0);
                heat += sample.dissipated_energy_j;
            }
            assert!((system.stored_energy_j() + heat - energy).abs() < 1e-8 * energy + 1e-11);
            let state = system.state();
            let error = ((state[0] - expected[0]) / amplitude)
                .abs()
                .max(((state[1] - expected[1]) / (k.sqrt() * amplitude)).abs())
                .max((heat - expected[3]).abs() / energy);
            errors.push(error);
            assert!(
                heat > 0.0 && system.relaxation_observation().unwrap().dissipated_power_w >= 0.0
            );
        }
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "second-order film history: {errors:?}"
        );
        assert!(
            errors[2] < 2e-4,
            "resolved physical film state/heat: {errors:?}"
        );
    }
}

#[test]
fn complete_mixed_pair_energy_survives_prepared_state_and_acoustic_composition() {
    let (d, _, rim) = film();
    let n = d.model.free;
    let mut k = vec![0.0; n * n];
    let mut m = k.clone();
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = d.model.k.get(i, j);
            m[i * n + j] = d.model.m.get(i, j);
        }
    }
    let modes: Vec<_> = fs_modal::eigh_gen_dense(&k, &m, n)
        .unwrap()
        .into_iter()
        .take(2)
        .collect();
    let p = MembranePotential::from_pencil(
        &d.mesh,
        &d.section,
        &d.model,
        &modes,
        &rim,
        MembraneReductionBudget {
            max_modes: 4,
            max_nodes: 100,
            max_facet_pairs: 1000,
            max_solve_entries: 100000,
            relative_tolerance: 1e-6,
        },
        0.25,
    )
    .unwrap();
    let q = [
        0.04 / p.reduction().maximum_slope(&[1.0, 0.0]),
        -0.03 / p.reduction().maximum_slope(&[0.0, 1.0]),
    ];
    let extra = RATIO * p.reduction().stretching_energy(&q);
    let equilibrium = p.reduction().potential(&q);
    let law = p
        .with_proportional_relaxation(
            spectrum(),
            q.iter().map(|q| 2.0 * q.abs()).collect(),
            InitialMemory::Unrelaxed,
            budget(),
        )
        .unwrap();
    assert_eq!(law.relaxation_state_count(), 3);
    let make = || {
        ImpactSystem::new(
            vec![ImpactBody {
                potential: BodyPotential::Membrane(law.clone()),
                initial: q
                    .iter()
                    .map(|q| ModalAcousticState {
                        displacement_m_sqrt_kg: *q,
                        velocity_m_sqrt_kg_per_s: 0.0,
                    })
                    .collect(),
                damping_per_s: vec![0.0; 2],
            }],
            vec![],
            vec![],
            vec![],
            config(),
        )
        .unwrap()
    };
    let system = make();
    assert!((system.stored_energy_j() - equilibrium - extra).abs() < 2e-11 * (equilibrium + extra));
    use fs_couple::render::plate::impact::radiation::{Model, Pole};
    let air = Model {
        ports: 2,
        poles: vec![Pole {
            omega: 2000.0,
            zeta: 0.1,
            coupling: vec![10.0, -5.0],
        }],
    };
    let mut prepared = system
        .with_radiation_load(&air, &[0, 1], 1)
        .unwrap()
        .prepare_analytic()
        .unwrap();
    let mut reference = make()
        .with_radiation_load(&air, &[0, 1], 1)
        .unwrap()
        .prepare()
        .unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..32 {
        prepared.step(&[0.0; 2], &gate).unwrap();
        reference.step(&[0.0; 2], &gate).unwrap();
        for (a, b) in prepared.state().iter().zip(reference.state()) {
            assert!((a - b).abs() < 2e-7 * a.abs().max(b.abs()).max(1e-5));
        }
    }
    let state = prepared.state().to_vec();
    let memory = prepared.relaxation_observation().unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(prepared.step(&[0.0; 2], &cancelled).is_err());
    assert_eq!(prepared.state(), state);
    assert_eq!(prepared.relaxation_observation().unwrap(), memory);
}

#[test]
fn film_source_history_clock_and_duplicate_loss_refuse_without_fallback() {
    let p = reduced(0.25);
    let q = 0.04 / p.reduction().maximum_slope(&[1.0]);
    let mut small = budget();
    small.maximum_memory_states = 0;
    assert!(
        p.clone()
            .with_proportional_relaxation(spectrum(), vec![q], InitialMemory::Relaxed, small)
            .is_err()
    );
    let mut small = budget();
    small.maximum_setup_terms = 1;
    assert!(
        p.clone()
            .with_proportional_relaxation(spectrum(), vec![q], InitialMemory::Relaxed, small)
            .is_err()
    );
    assert!(
        p.clone()
            .with_proportional_relaxation(spectrum(), vec![], InitialMemory::Relaxed, budget())
            .is_err()
    );
    let law = loaded(p.clone(), q, InitialMemory::Relaxed);
    assert!(
        law.clone()
            .with_proportional_relaxation(spectrum(), vec![q], InitialMemory::Relaxed, budget())
            .is_err()
    );
    let mut duplicate = body(law, q, 0.0);
    duplicate.damping_per_s[0] = 1.0;
    assert!(ImpactSystem::new(vec![duplicate], vec![], vec![], vec![], config()).is_err());
    let mut fast = spectrum();
    fast.branches[0].1 = config().dt_s / 3.0;
    let law = p
        .clone()
        .with_proportional_relaxation(fast, vec![2.0 * q], InitialMemory::Relaxed, budget())
        .unwrap();
    assert!(ImpactSystem::new(vec![body(law, q, 0.0)], vec![], vec![], vec![], config()).is_err());
    let mut narrow = spectrum();
    narrow.omega_band_rad_s.1 = 1.0;
    let law = p
        .clone()
        .with_proportional_relaxation(narrow, vec![2.0 * q], InitialMemory::Relaxed, budget())
        .unwrap();
    assert!(ImpactSystem::new(vec![body(law, q, 0.0)], vec![], vec![], vec![], config()).is_err());
    let law = p
        .with_proportional_relaxation(spectrum(), vec![q], InitialMemory::Relaxed, budget())
        .unwrap();
    assert!(
        ImpactSystem::new(
            vec![body(law, 1.01 * q, 0.0)],
            vec![],
            vec![],
            vec![],
            config()
        )
        .is_err()
    );
}

#[test]
fn exponential_film_retains_original_history_heat_and_refinement() {
    use fs_couple::render::plate::impact::ExponentialImpactSystem;
    let p = reduced(0.25);
    let k = p.reduction().linear_stiffness_matrix()[0];
    let beta = 4.0 * p.reduction().stretching_energy(&[1.0]);
    let amplitude = 0.06 / p.reduction().maximum_slope(&[1.0]);
    let duration = 192.0 * config().dt_s;
    let gate = CancelGate::new_clock_free();
    for relaxed in [false, true] {
        let reference = continuous(k, beta, amplitude, relaxed, duration);
        let mut errors = Vec::new();
        for steps in [24, 48, 96] {
            let law = loaded(
                p.clone(),
                amplitude,
                if relaxed {
                    InitialMemory::Relaxed
                } else {
                    InitialMemory::Unrelaxed
                },
            );
            let mut cfg = config();
            cfg.dt_s = duration / f64::from(steps);
            let mut system = ExponentialImpactSystem::new(
                vec![body(law, amplitude, 0.0)],
                vec![],
                vec![],
                cfg,
                &gate,
            )
            .unwrap();
            assert_eq!(system.state().len(), 3);
            assert_eq!(system.relaxation_branch_count(), 1);
            let energy = system.stored_energy_j();
            let mut heat = 0.0;
            for _ in 0..steps {
                let frame = system.step(&[0.0], &gate).unwrap();
                assert!(frame.dissipated_energy_j >= 0.0);
                assert!(
                    frame.balance_residual_j.abs()
                        <= cfg.energy_absolute_tolerance_j + cfg.energy_relative_tolerance * energy
                );
                heat += frame.dissipated_energy_j;
            }
            assert!((system.stored_energy_j() + heat - energy).abs() < 1e-8 * energy + 1e-11);
            let x = system.state();
            let error = ((x[0] - reference[0]) / amplitude)
                .abs()
                .max(((x[1] - reference[1]) / (k.sqrt() * amplitude)).abs())
                .max((heat - reference[3]).abs() / energy);
            assert!(heat > 0.0);
            errors.push(error);
        }
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "exponential nonlinear film state/heat: {errors:?}"
        );
        assert!(errors[2] < 3e-4, "resolved film state and heat: {errors:?}");
    }
}

fn two_mode_film() -> MembranePotential {
    let (disk, _, rim) = film();
    let n = disk.model.free;
    let mut k = vec![0.0; n * n];
    let mut m = k.clone();
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = disk.model.k.get(i, j);
            m[i * n + j] = disk.model.m.get(i, j);
        }
    }
    let modes: Vec<_> = fs_modal::eigh_gen_dense(&k, &m, n)
        .unwrap()
        .into_iter()
        .take(2)
        .collect();
    MembranePotential::from_pencil(
        &disk.mesh,
        &disk.section,
        &disk.model,
        &modes,
        &rim,
        MembraneReductionBudget {
            max_modes: 4,
            max_nodes: 100,
            max_facet_pairs: 1000,
            max_solve_entries: 100000,
            relative_tolerance: 1e-6,
        },
        0.25,
    )
    .unwrap()
}

#[test]
fn exponential_mixed_film_memory_and_radiation_keep_complete_state_and_atomic_acceptance() {
    use fs_couple::render::plate::impact::ExponentialImpactSystem;
    use fs_couple::render::plate::impact::radiation::{Model, Pole};
    let p = two_mode_film();
    let q = [
        0.04 / p.reduction().maximum_slope(&[1.0, 0.0]),
        -0.025 / p.reduction().maximum_slope(&[0.0, 1.0]),
    ];
    let expected = p.reduction().potential(&q) + RATIO * p.reduction().stretching_energy(&q);
    let law = p
        .with_proportional_relaxation(
            spectrum(),
            q.iter().map(|q| 2.0 * q.abs()).collect(),
            InitialMemory::Unrelaxed,
            budget(),
        )
        .unwrap();
    let physical = ImpactBody {
        potential: BodyPotential::Membrane(law),
        initial: q
            .iter()
            .map(|q| ModalAcousticState {
                displacement_m_sqrt_kg: *q,
                velocity_m_sqrt_kg_per_s: 0.0,
            })
            .collect(),
        damping_per_s: vec![0.0; 2],
    };
    let gate = CancelGate::new_clock_free();
    let mut system = ExponentialImpactSystem::new(vec![physical], vec![], vec![], config(), &gate)
        .unwrap()
        .with_radiation_load(
            &Model {
                ports: 2,
                poles: vec![Pole {
                    omega: 2000.0,
                    zeta: 0.1,
                    coupling: vec![10.0, -5.0],
                }],
            },
            &[1, 0],
            1,
            &gate,
        )
        .unwrap();
    assert_eq!(system.relaxation_branch_count(), 3);
    assert_eq!(system.state().len(), 4 + 3 + 2);
    assert!((system.stored_energy_j() - expected).abs() < 3e-11 * expected);
    let mut heat = 0.0;
    for _ in 0..32 {
        let frame = system.step(&[0.0; 2], &gate).unwrap();
        assert!(frame.dissipated_energy_j >= 0.0);
        heat += frame.dissipated_energy_j;
    }
    assert!(heat > 0.0);
    assert!((system.stored_energy_j() + heat - expected).abs() < 1e-8 * expected + 1e-11);
    let state = system.state().to_vec();
    let history = system.relaxation_observation();
    let sample = system.samples();
    let membrane = system.membrane_observation(0).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.step(&[0.0; 2], &cancelled),
        Err(ImpactError::Cancelled)
    ));
    system.set_iteration_limit(0).unwrap();
    assert!(system.step(&[0.0; 2], &gate).is_err());
    assert_eq!(system.state(), state);
    assert_eq!(system.relaxation_observation(), history);
    assert_eq!(system.samples(), sample);
    assert_eq!(system.membrane_observation(0).unwrap(), membrane);
    system.set_iteration_limit(50).unwrap();
    system.step(&[0.0; 2], &gate).unwrap();
    assert_eq!(system.samples(), sample + 1);
}

#[test]
fn exponential_film_preserves_source_domains_and_counts_every_history_before_preparation() {
    use fs_couple::render::plate::impact::ExponentialImpactSystem;
    let gate = CancelGate::new_clock_free();
    let p = reduced(0.25);
    let q = 0.04 / p.reduction().maximum_slope(&[1.0]);
    assert!(matches!(
        ExponentialImpactSystem::new(
            vec![body(p.clone(), q, 0.0)],
            vec![],
            vec![],
            config(),
            &gate
        ),
        Err(ImpactError::Invalid(
            "exponential nonlinear films require their original declared strain/history envelope"
        ))
    ));
    let bounded_elastic = p
        .clone()
        .with_proportional_relaxation(
            MembraneRelaxationSpectrum {
                branches: vec![],
                omega_band_rad_s: (0.0, 1e7),
            },
            vec![2.0 * q],
            InitialMemory::Relaxed,
            budget(),
        )
        .unwrap();
    let mut lossless = ExponentialImpactSystem::new(
        vec![body(bounded_elastic, q, 0.0)],
        vec![],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    assert_eq!(lossless.state().len(), 2);
    assert_eq!(lossless.relaxation_branch_count(), 0);
    assert_eq!(lossless.stored_energy_j(), p.reduction().potential(&[q]));
    assert_eq!(
        lossless.step(&[0.0], &gate).unwrap().dissipated_energy_j,
        0.0
    );

    for (band, tau) in [(1.0, TAU), (1e7, config().dt_s / 3.0)] {
        let law = p
            .clone()
            .with_proportional_relaxation(
                MembraneRelaxationSpectrum {
                    branches: vec![(RATIO, tau)],
                    omega_band_rad_s: (0.0, band),
                },
                vec![2.0 * q],
                InitialMemory::Relaxed,
                budget(),
            )
            .unwrap();
        assert!(
            ExponentialImpactSystem::new(vec![body(law, q, 0.0)], vec![], vec![], config(), &gate)
                .is_err()
        );
    }
    let law = p
        .with_proportional_relaxation(spectrum(), vec![q], InitialMemory::Relaxed, budget())
        .unwrap();
    assert!(
        ExponentialImpactSystem::new(
            vec![body(law, 1.01 * q, 0.0)],
            vec![],
            vec![],
            config(),
            &gate
        )
        .is_err()
    );
    let p = two_mode_film();
    let law = p
        .with_proportional_relaxation(
            MembraneRelaxationSpectrum {
                branches: vec![(0.1, TAU); 8],
                omega_band_rad_s: (0.0, 1e7),
            },
            vec![1e-4; 2],
            InitialMemory::Unrelaxed,
            budget(),
        )
        .unwrap();
    assert_eq!(law.relaxation_state_count(), 24);
    let physical = ImpactBody {
        potential: BodyPotential::Membrane(law),
        initial: vec![
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 0.0
            };
            2
        ],
        damping_per_s: vec![0.0; 2],
    };
    // 10 mechanical modes plus 120 genuine material coordinates exceed the
    // unchanged complete exact-flow cache. No arm or mode may be discarded.
    assert!(matches!(
        ExponentialImpactSystem::new(vec![physical; 5], vec![], vec![], config(), &gate),
        Err(ImpactError::Invalid(
            "complete exponential impact state exceeds the original dense owner capacity"
        ))
    ));
}
