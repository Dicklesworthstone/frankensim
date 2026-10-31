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
