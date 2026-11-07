//! Omitted-mode compliance among fixed obstacles (pexfp.299, MagentaBadger
//! 7592 gate 3): the joint solve sees each support in series with it, the
//! ledger stores f G f / 2, and the step-held quasi-static response books its
//! declared numerical loss dR G dR / 2. The same holds beside bound nonlinear
//! storage, in its steps and its static preload.

use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;

// obstacle(): the synthetic support's linear (alpha 1) stiffness [N/m].
const SUPPORT_N_M: f64 = 2000.0;
// 2 K g = 1 for the opposed pair: its net stiffness 2K halves.
const PAIR_G: f64 = 2.5e-4;

fn moving(velocity: f64) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(100.0, velocity)], vec![], config(1000.0), &CancelGate::new_clock_free()).unwrap()
}
fn pair(compliance: Option<f64>, velocity: f64) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let contacts = [1.0, -1.0].map(|sign| (obstacle(sign, -1e-4), normal_limits())).to_vec();
    let joint = multi(friction(moving(velocity))).with_obstacles(contacts, limits(), &gate).unwrap();
    match compliance {
        Some(g) => joint.with_obstacle_residual_compliance(vec![g, -g, -g, g], &gate).unwrap(),
        None => joint,
    }
}
fn single(compliance: Option<f64>) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let joint = multi(friction(moving(0.0)))
        .with_obstacles(vec![(obstacle(1.0, -1e-4), normal_limits())], limits(), &gate).unwrap();
    match compliance {
        Some(g) => joint.with_obstacle_residual_compliance(vec![g], &gate).unwrap(),
        None => joint,
    }
}

// Damped angular frequency [rad/s] from the modal velocity's zero crossings.
fn frequency(joint: &mut MultiFrictionModalSystem, steps: u32) -> f64 {
    let dt = joint.sample_period_s();
    let mut crossings = Vec::new();
    let mut last = joint.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
    for step in 1..=steps {
        let frame = joint.step(&[0.0], &[released()]).unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let v = joint.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
        if last != 0.0 && v.signum() != last.signum() { crossings.push(dt*(f64::from(step)-v/(v-last))); }
        last = v;
    }
    assert!(crossings.len() >= 4, "only {} crossings", crossings.len());
    core::f64::consts::PI*(crossings.len()-1) as f64/(crossings[crossings.len()-1]-crossings[0])
}

#[test]
fn supports_act_in_series_with_the_omitted_modes() {
    // Pair: net stiffness 2K/(1 + 2Kg) on a 100 rad/s mode; single: K/(1 + Kg).
    // The 0.02 modal damping moves these by under 2e-4.
    for (mut joint, expected, alone) in [
        (pair(Some(PAIR_G), 0.004), 1e4+2.0*SUPPORT_N_M/(1.0+2.0*SUPPORT_N_M*PAIR_G), false),
        (pair(None, 0.004), 1e4+2.0*SUPPORT_N_M, false),
        (single(Some(2.0*PAIR_G)), 1e4+SUPPORT_N_M/(1.0+2.0*SUPPORT_N_M*PAIR_G), true),
        (single(None), 1e4+SUPPORT_N_M, true),
    ] {
        let paired = joint.obstacle_pairs();
        let bound = joint.obstacle_residual_compliance().is_some();
        let measured = frequency(&mut joint, 12_000);
        assert!((measured/expected.sqrt()-1.0).abs() < 1e-3,
            "single {alone} residual {bound}: {measured} rad/s against {}", expected.sqrt());
        assert_eq!(paired, if bound && !alone { vec![(0, 1)] } else { vec![] });
        for n in joint.last_obstacle_frames().unwrap() { assert!(n.normal_force_n > 0.0, "a support released"); }
    }
}

#[test]
fn the_ledger_needs_the_omitted_modes_storage_and_books_their_numerical_loss() {
    let mut joint = pair(Some(PAIR_G), 0.004);
    let (mut stored, mut net) = (0.0, 0.0);
    let mut unbooked = 0;
    let steps = 4096;
    for _ in 0..steps {
        let f = joint.step(&[0.0], &[released()]).unwrap().clone();
        let n = joint.last_obstacle_frames().unwrap();
        // G = g [[1, -1], [-1, 1]] sees only the net force.
        let now = n[0].normal_force_n-n[1].normal_force_n;
        assert!((f.obstacle_residual_energy_j-0.5*PAIR_G*now*now).abs() <= 1e-18);
        let change = now-net;
        assert!((f.obstacle_residual_loss_j-0.5*PAIR_G*change*change).abs() <= 1e-18);
        assert!(f.obstacle_residual_loss_j >= 0.0);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        // The same step's identity with f G f / 2 dropped from both ends.
        let dropped = f.energy_residual_j-(f.obstacle_residual_energy_j-stored);
        if dropped.abs() > f.energy_tolerance_j { unbooked += 1; }
        (stored, net) = (f.obstacle_residual_energy_j, now);
    }
    assert!(unbooked > 3*steps/4, "without f G f / 2 the identity still closed on {} of {steps} steps", steps-unbooked);
}

#[test]
fn residual_compliance_refuses_what_it_does_not_own() {
    let gate = CancelGate::new_clock_free();
    assert!(multi(friction(moving(0.0))).with_obstacle_residual_compliance(vec![1e-4], &gate).is_err());
    let pinch = || vec![PAIR_G, -PAIR_G, -PAIR_G, PAIR_G];
    for bad in [vec![1e-4], vec![1e-4, 0.0, 0.0], vec![-1e-4, 0.0, 0.0, 1e-4], vec![1e-4, 2e-4, 2e-4, 1e-4],
        vec![1e-4, 1e-5, 1.000_000_1e-5, 1e-4], vec![f64::NAN, 0.0, 0.0, 1e-4]] {
        assert!(pair(None, 0.0).with_obstacle_residual_compliance(bad, &gate).is_err());
    }
    assert!(pair(Some(PAIR_G), 0.0).with_obstacle_residual_compliance(pinch(), &gate).is_err());
    let mut played = pair(None, 0.0);
    played.step(&[0.0], &[released()]).unwrap();
    assert!(played.with_obstacle_residual_compliance(pinch(), &gate).is_err());
}

// Two Hertz stops on a string-like component: retained modes 300 and 800 rad/s,
// omitted 900 and 1400 rad/s. Omitted shapes are the last two of each row.
const STOP_SHAPES: [[f64; 4]; 2] = [[1.0, 0.6, 0.5, -0.4], [0.7, -0.9, 0.3, 0.6]];
const STOP_OMEGAS: [f64; 4] = [300.0, 800.0, 900.0, 1400.0];
const STOP_GAPS: [f64; 2] = [2e-5, 1e-5];
// A retained-mode load: the dense comparison gives the omitted modes no external force.
const STOP_LOAD: [f64; 2] = [40.0, 10.0];

fn modes(count: usize) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48000, STOP_OMEGAS[..count].iter().map(|&omega| ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO,
    }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn stiff_stop(shapes: &[f64], gap: f64, stiffness: f64) -> ModalContact {
    ModalContact { left: ModalAttachment { component: 0, shapes: shapes.to_vec() },
        right: ModalAttachment { component: 0, shapes: vec![0.0; shapes.len()] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], stiffness, 1.5,
            "synthetic Hertz support".into()).unwrap() }
}
fn stop(shapes: &[f64], gap: f64) -> ModalContact { stiff_stop(shapes, gap, 2e7) }
fn joint_over(model: ModalAcousticTimeModel, coupling: ModalCouplingConfig) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut port = vec![0.0; model.modes().len()];
    port[0] = 1.0;
    let network = CoupledModalSystem::new(vec![model], vec![], coupling, &gate).unwrap();
    multi(FrictionModalSystem::new(network, ModalFrictionPort {
        left: ModalAttachment { component: 0, shapes: port }, right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } },
        ModalFrictionConfig { max_iterations: 96, maximum_force_n: 50.0, maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-12 },
        &gate).unwrap())
}
fn retained_joint() -> MultiFrictionModalSystem { joint_over(modes(2), config(1000.0)) }
fn omitted_compliance() -> Vec<f64> {
    let mut g = vec![0.0; 4];
    for i in 0..2 { for j in 0..2 { for m in 2..4 {
        g[i*2+j] += STOP_SHAPES[i][m]*STOP_SHAPES[j][m]/(STOP_OMEGAS[m]*STOP_OMEGAS[m]);
    } } }
    g
}
fn displacements(joint: &MultiFrictionModalSystem) -> Vec<f64> {
    joint.components()[0].states().iter().map(|s| s.displacement_m_sqrt_kg).collect()
}
fn accepted_reactions(joint: &MultiFrictionModalSystem) -> Vec<f64> {
    joint.obstacles.as_ref().unwrap().accepted.iter().map(|f| f.normal_force_n).collect()
}

#[test]
fn released_preload_matches_the_dense_solve_with_the_omitted_modes_kept() {
    use crate::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
    let gate = CancelGate::new_clock_free();
    // The original owner on all four modes is the reference.
    let dense_contacts: Vec<_> = STOP_SHAPES.iter().zip(STOP_GAPS)
        .map(|(shapes, gap)| (stop(shapes, gap), normal_limits())).collect();
    let mut dense = CoupledModalSystem::new(vec![modes(4)], vec![], config(1000.0), &gate).unwrap();
    let dense_energy = dense.initialize_contact_equilibrium(&[STOP_LOAD[0], STOP_LOAD[1], 0.0, 0.0], &dense_contacts,
        MultiContactConfig { max_contacts: 8, max_sweeps: 128, max_setup_terms: 100000 }, &gate).unwrap();
    let dense_q: Vec<f64> = dense.components()[0].states().iter().map(|s| s.displacement_m_sqrt_kg).collect();
    let dense_reactions: Vec<f64> = STOP_SHAPES.iter().zip(STOP_GAPS).map(|(shapes, gap)| {
        let x: f64 = shapes.iter().zip(&dense_q).map(|(a, b)| a*b).sum();
        2e7*(x-gap).powf(1.5)
    }).collect();

    // Two modes plus G, with part of each retained column carried by the correction map.
    let correction = [0.2, -0.1];
    let contacts: Vec<_> = STOP_SHAPES.iter().zip(STOP_GAPS)
        .map(|(shapes, gap)| (stop(&[shapes[0]-correction[0], shapes[1]-correction[1]], gap), normal_limits())).collect();
    let corrections = vec![vec![ModalAttachment { component: 0, shapes: correction.to_vec() }]; 2];
    let mut joint = retained_joint().with_obstacle_residual_ports_and_corrections(contacts, corrections,
        omitted_compliance(), limits(), &gate).unwrap();
    let energy = joint.initialize_released_equilibrium(&STOP_LOAD, &gate).unwrap();
    let (q, reactions) = (displacements(&joint), accepted_reactions(&joint));

    for k in 0..2 { assert!((q[k]-dense_q[k]).abs() <= 1e-9*dense_q[k].abs(), "q[{k}]: {} against {}", q[k], dense_q[k]); }
    for i in 0..2 {
        assert!((reactions[i]-dense_reactions[i]).abs() <= 1e-8*dense_reactions[i], "R[{i}]: {} against {}", reactions[i], dense_reactions[i]);
    }
    assert!((energy-dense_energy).abs() <= 1e-9*dense_energy, "energy {energy} against {dense_energy}");
    assert!((energy-joint.total_energy_j().unwrap()).abs() <= 1e-15, "the ledger reads a different start");
    // Force balance on the retained modes, with the complete columns.
    for k in 0..2 {
        let load: f64 = STOP_SHAPES.iter().zip(&reactions).map(|(shapes, r)| shapes[k]*r).sum();
        let omega2 = STOP_OMEGAS[k]*STOP_OMEGAS[k];
        assert!((omega2*q[k]-(STOP_LOAD[k]-load)).abs() <= 1e-7, "mode {k} is not in balance");
    }
    assert!(joint.last_obstacle_frames().is_none() && joint.samples_rendered() == 0, "a report was fabricated");

    // Without G the same stops settle measurably elsewhere, so this test can fail.
    let plain_contacts: Vec<_> = STOP_SHAPES.iter().zip(STOP_GAPS)
        .map(|(shapes, gap)| (stop(&shapes[..2], gap), normal_limits())).collect();
    let mut plain = retained_joint().with_obstacles(plain_contacts, limits(), &gate).unwrap();
    plain.initialize_released_equilibrium(&STOP_LOAD, &gate).unwrap();
    assert!((displacements(&plain)[0]-dense_q[0]).abs() > 0.01*dense_q[0].abs());

    // The held load keeps it: no first-step jump in the omitted modes' lag.
    for _ in 0..64 {
        let frame = joint.step(&STOP_LOAD, &[released()]).unwrap().clone();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.obstacle_residual_loss_j <= 1e-18, "the settled omitted modes jumped: {}", frame.obstacle_residual_loss_j);
        for (now, then) in joint.last_obstacle_frames().unwrap().iter().zip(&reactions) {
            assert!((now.normal_force_n-then).abs() <= 1e-6*then);
        }
        for state in joint.components()[0].states() { assert!(state.velocity_m_sqrt_kg_per_s.abs() <= 1e-8); }
    }
    let release = joint.step(&[0.0, 0.0], &[released()]).unwrap();
    assert!(release.energy_residual_j.abs() <= release.energy_tolerance_j);
}

#[test]
fn released_preload_settles_each_support_in_series_with_the_omitted_modes() {
    let gate = CancelGate::new_clock_free();
    let k = SUPPORT_N_M;
    // Opposed pair, closed by 1e-4 m each: G sees the net force, whose stiffness 2K halves.
    let mut joint = pair(Some(PAIR_G), 0.0);
    let energy = joint.initialize_released_equilibrium(&[0.1], &gate).unwrap();
    let q = 0.1/(1e4+2.0*k/(1.0+2.0*k*PAIR_G));
    let net = 2.0*k*q/(1.0+2.0*k*PAIR_G);
    assert!((displacements(&joint)[0]-q).abs() <= 1e-9*q);
    let reactions = accepted_reactions(&joint);
    assert!((reactions[0]-reactions[1]-net).abs() <= 1e-9*net);
    assert!((reactions[0]-k*(q-PAIR_G*net+1e-4)).abs() <= 1e-9*reactions[0]);
    assert!((energy-joint.total_energy_j().unwrap()).abs() <= 1e-15);
    // One support: x = q - g R, closed by 1e-4 m.
    let g = 2.0*PAIR_G;
    let mut joint = single(Some(g));
    joint.initialize_released_equilibrium(&[0.2], &gate).unwrap();
    let series = k/(1.0+k*g);
    let q = (0.2-series*1e-4)/(1e4+series);
    assert!((displacements(&joint)[0]-q).abs() <= 1e-9*q.abs());
    assert!((accepted_reactions(&joint)[0]-series*(q+1e-4)).abs() <= 1e-9*series*(q+1e-4));
}

#[test]
fn a_stop_stiff_enough_to_miss_on_rounding_settles_on_the_measured_difference() {
    // A finger in small: a stop so stiff it closes by nanometres takes nearly all of a
    // load on its own column, so its closure is a tiny difference of far larger sums.
    // The tolerance is the bowed cards': 1e-9 of the reaction.
    const STIFFNESS: f64 = 1e13;
    const OMITTED: f64 = 2.5e-5;
    let gate = CancelGate::new_clock_free();
    let shape: Vec<f64> = (0..64).map(|m| (1.3*f64::from(m)+0.4).sin()).collect();
    let joint = |compliance: Option<f64>| {
        let model = ModalAcousticTimeModel::try_new(48000, (0..64).map(|m| ModalAcousticMode {
            angular_frequency_rad_s: 40.0+6.0*f64::from(m), damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO,
        }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
        let joint = joint_over(model, ModalCouplingConfig { max_modes: 64, max_setup_terms: 1_000_000, ..config(1000.0) });
        let mut tight = normal_limits();
        tight.force_relative_tolerance = 1e-9;
        let joint = joint.with_obstacles(vec![(stiff_stop(&shape, 0.0, STIFFNESS), tight)], limits(), &gate).unwrap();
        match compliance {
            Some(g) => joint.with_obstacle_residual_compliance(vec![g], &gate).unwrap(),
            None => joint,
        }
    };
    // Without the correction both paths miss on some of these loads (35db3c6df counted
    // misses without omitted modes, and its one-word mutation missed with them).
    for i in 0..256 {
        let scale = 2.0+f64::from(i)/128.0;
        let load: Vec<f64> = shape.iter().map(|s| s*scale).collect();
        for omitted in [None, Some(OMITTED)] {
            let mut settled = joint(omitted);
            let energy = settled.initialize_released_equilibrium(&load, &gate)
                .unwrap_or_else(|error| panic!("load scale {scale}, omitted {omitted:?}: {error}"));
            // Without omitted modes the solve is the network's own contact equilibrium, whose
            // gate checks this law on the actual rounded states before it returns; the obstacle
            // frames publish no reaction until a sample, so settling is what this path shows.
            let Some(g) = omitted else { continue };
            // The law on the accepted states, by this test's own sum.
            let reaction = accepted_reactions(&settled)[0];
            let closure = shape.iter().zip(displacements(&settled)).map(|(a, b)| a*b).sum::<f64>()-g*reaction;
            assert!((reaction-STIFFNESS*closure.powf(1.5)).abs() <= 1e-12+1e-9*reaction, "load scale {scale}");
            assert!((energy-settled.total_energy_j().unwrap()).abs() <= 1e-15, "the ledger reads a different start");
        }
    }
}

#[test]
fn a_zero_compliance_is_the_original_preload_to_the_bit() {
    let gate = CancelGate::new_clock_free();
    let (mut bound, mut plain) = (single(Some(0.0)), single(None));
    let energy = bound.initialize_released_equilibrium(&[0.1], &gate).unwrap();
    assert_eq!(energy.to_bits(), plain.initialize_released_equilibrium(&[0.1], &gate).unwrap().to_bits());
    assert_eq!(states(&bound), states(&plain));
}

#[test]
fn a_refused_residual_preload_leaves_the_held_reactions_alone() {
    let gate = CancelGate::new_clock_free();
    let mut small = normal_limits();
    small.maximum_force_n = 1e-3;
    let network = CoupledModalSystem::new(vec![model(100.0, 0.0)], vec![], config(1000.0), &gate).unwrap();
    let mut joint = multi(friction(network)).with_obstacles(vec![(obstacle(1.0, -1e-4), small)], limits(), &gate).unwrap()
        .with_obstacle_residual_compliance(vec![1e-4], &gate).unwrap();
    let before = joint.total_energy_j().unwrap();
    assert!(joint.initialize_released_equilibrium(&[0.1], &gate).is_err());
    assert_eq!(accepted_reactions(&joint), vec![0.0]);
    assert_eq!(joint.total_energy_j().unwrap().to_bits(), before.to_bits());
    assert_eq!(displacements(&joint), vec![0.0]);
}

// Bound nonlinear storage beside the condensed supports: a quartic hardening of
// the stopped mode, bound before the supports, is the cheap owner of the
// nonlinear path. Its outer balance carries f G f / 2 and dR G dR / 2, and its
// static preload settles each support in series with G.
const HARDENING: f64 = 1e12;
struct Hardening;
impl Storage for Hardening {
    fn hamiltonian(&self, q: &[f64]) -> f64 { 0.25*HARDENING*q[0].powi(4) }
    fn gradient(&self, q: &[f64], out: &mut [f64]) { out[0] = HARDENING*q[0].powi(3); }
}
impl super::super::super::super::nonlinear::ModalPotential for Hardening {
    fn coordinate_count(&self) -> usize { 1 }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() != 1 || !q[0].is_finite() { return Err(invalid("hardening test coordinate")); }
        Ok(())
    }
}
fn hardened(signs: &[f64], compliance: Option<Vec<f64>>, velocity: f64) -> MultiFrictionModalSystem {
    use super::super::super::super::nonlinear::{ModalPotentialBinding, NonlinearModalConfig};
    let gate = CancelGate::new_clock_free();
    let core = friction(moving(velocity)).with_nonlinear_potentials(vec![ModalPotentialBinding {
        component: 0, coordinate_start: 0, potential: Box::new(Hardening) }], NonlinearModalConfig {
        max_iterations: 64, maximum_generalized_force: 1e5, force_absolute_tolerance: 1e-11,
        force_relative_tolerance: 1e-10 }, &gate).unwrap();
    let contacts = signs.iter().map(|&sign| (obstacle(sign, -1e-4), normal_limits())).collect();
    let joint = multi(core).with_obstacles(contacts, limits(), &gate).unwrap();
    match compliance {
        Some(g) => joint.with_obstacle_residual_compliance(g, &gate).unwrap(),
        None => joint,
    }
}

fn translation_limits() -> ObstacleMotionLimits {
    ObstacleMotionLimits { maximum_abs_offset_m: 0.01, maximum_speed_m_s: 1.0 }
}

// Two original bodies, with a prescribed translation added to their relative
// contact coordinate. Undamped oscillators give an independent analytic step.
fn translated_residual(motion_first: bool, maximum_pressure: f64) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let models = [300.0, 500.0].into_iter().map(|omega| {
        ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
            angular_frequency_rad_s: omega, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }], ModalAcousticTimeBudget::audible_reference()).unwrap()
    }).collect();
    let network = CoupledModalSystem::new(models, vec![], config(maximum_pressure), &gate).unwrap();
    let core = friction(network).with_acceleration_pressure(vec![1.0, 0.0]).unwrap();
    let mut contact = obstacle(1.0, 0.0);
    contact.right = attachment(1, 0.8);
    let joint = multi(core).with_obstacles(vec![(contact, normal_limits())], limits(), &gate).unwrap();
    if motion_first {
        joint.with_obstacle_motion(translation_limits(), &gate).unwrap()
            .with_obstacle_residual_compliance(vec![PAIR_G], &gate).unwrap()
    } else {
        joint.with_obstacle_residual_compliance(vec![PAIR_G], &gate).unwrap()
            .with_obstacle_motion(translation_limits(), &gate).unwrap()
    }
}

#[test]
fn nonlinear_storage_books_the_omitted_modes_storage_and_numerical_loss() {
    let mut joint = hardened(&[1.0, -1.0], Some(vec![PAIR_G, -PAIR_G, -PAIR_G, PAIR_G]), 0.004);
    assert_eq!(joint.obstacle_pairs(), vec![(0, 1)]);
    let (mut stored, mut net) = (0.0, 0.0);
    let mut unbooked = 0;
    let steps = 4096;
    for _ in 0..steps {
        let f = joint.step(&[0.0], &[released()]).unwrap().clone();
        assert!(f.nonlinear.is_some(), "the step left the nonlinear owner");
        let n = joint.last_obstacle_frames().unwrap();
        let now = n[0].normal_force_n-n[1].normal_force_n;
        assert!((f.obstacle_residual_energy_j-0.5*PAIR_G*now*now).abs() <= 1e-18);
        let change = now-net;
        assert!((f.obstacle_residual_loss_j-0.5*PAIR_G*change*change).abs() <= 1e-18);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        // The same step's identity with f G f / 2 dropped from both ends.
        let dropped = f.energy_residual_j-(f.obstacle_residual_energy_j-stored);
        if dropped.abs() > f.energy_tolerance_j { unbooked += 1; }
        (stored, net) = (f.obstacle_residual_energy_j, now);
    }
    assert!(unbooked > 3*steps/4, "without f G f / 2 the identity still closed on {} of {steps} steps", steps-unbooked);
}

#[test]
fn nonlinear_storage_with_a_zero_compliance_is_the_original_preload_to_the_bit() {
    let gate = CancelGate::new_clock_free();
    let (mut bound, mut plain) = (hardened(&[1.0], Some(vec![0.0]), 0.0), hardened(&[1.0], None, 0.0));
    let energy = bound.initialize_released_equilibrium(&[0.5], &gate).unwrap();
    assert_eq!(energy.to_bits(), plain.initialize_released_equilibrium(&[0.5], &gate).unwrap().to_bits());
    assert_eq!(states(&bound), states(&plain));
    assert_eq!(bound.nonlinear_preload_diagnostics(), plain.nonlinear_preload_diagnostics());
}

#[test]
fn nonlinear_storage_preload_settles_the_support_in_series_with_the_omitted_modes() {
    let gate = CancelGate::new_clock_free();
    let (k, load) = (SUPPORT_N_M, 0.5);
    // k g = 1: the support's series stiffness halves.
    let g = 2.0*PAIR_G;
    let mut joint = hardened(&[1.0], Some(vec![g]), 0.0);
    let energy = joint.initialize_released_equilibrium(&[load], &gate).unwrap();
    let report = *joint.nonlinear_preload_diagnostics().unwrap();
    assert!(report.force_residual <= report.force_tolerance);
    let (q, r) = (displacements(&joint)[0], accepted_reactions(&joint)[0]);
    // The support's law at the closure its own reaction implies: x = q - g R, closed by 1e-4 m.
    let penetration = q-g*r+1e-4;
    assert!((r-k*penetration).abs() <= 1e-9*r, "R {r} against k x {}", k*penetration);
    // The retained mode balances the load against its linear, hardening and support forces.
    let balance = 1e4*q+HARDENING*q.powi(3)+r;
    assert!((balance-load).abs() <= 1e-9*load, "balance {balance} against {load}");
    assert!(HARDENING*q.powi(3) > 1e-3*load, "the hardening must matter here");
    // Stored: the mode, its hardening, the support and the omitted modes' f G f / 2.
    let expected = 0.5*1e4*q*q+0.25*HARDENING*q.powi(4)+0.5*k*penetration*penetration+0.5*g*r*r;
    assert!((energy-expected).abs() <= 1e-9*expected, "energy {energy} against {expected}");
    assert!((energy-joint.total_energy_j().unwrap()).abs() <= 1e-15, "the ledger reads a different start");
    assert!(joint.last_obstacle_frames().is_none() && joint.samples_rendered() == 0, "a report was fabricated");

    // Without G the same support settles measurably elsewhere, so this test can fail.
    let mut plain = hardened(&[1.0], None, 0.0);
    plain.initialize_released_equilibrium(&[load], &gate).unwrap();
    assert!((displacements(&plain)[0]-q).abs() > 0.01*q.abs());

    // The held load keeps it: no first-step jump in the omitted modes' lag.
    for _ in 0..64 {
        let frame = joint.step(&[load], &[released()]).unwrap().clone();
        assert!(frame.nonlinear.is_some(), "the step left the nonlinear owner");
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.obstacle_residual_loss_j <= 1e-18, "the settled omitted modes jumped: {}", frame.obstacle_residual_loss_j);
        let now = joint.last_obstacle_frames().unwrap()[0].normal_force_n;
        assert!((now-r).abs() <= 1e-6*r, "R {now} against the preload's {r}");
        for state in joint.components()[0].states() { assert!(state.velocity_m_sqrt_kg_per_s.abs() <= 1e-8); }
    }
}

#[test]
fn nonlinear_storage_refuses_moving_supports_beside_residual_compliance() {
    let gate = CancelGate::new_clock_free();
    let travel = ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0 };
    let moving = hardened(&[1.0], None, 0.0).with_obstacle_motion(travel, &gate).unwrap();
    assert!(moving.with_obstacle_residual_compliance(vec![1e-4], &gate).is_err());
    // The reverse binding order refuses too: the compliance first, then the translation.
    let compliant = hardened(&[1.0], Some(vec![1e-4]), 0.0);
    assert!(compliant.with_obstacle_motion(travel, &gate).is_err());
    let mut played = hardened(&[1.0], None, 0.0);
    played.step(&[0.0], &[released()]).unwrap();
    assert!(played.with_obstacle_residual_compliance(vec![1e-4], &gate).is_err());
    assert!(hardened(&[1.0], None, 0.0).with_obstacle_residual_compliance(vec![1e-4], &gate).is_ok());
}

#[test]
fn g1_translated_residual_contact_matches_analytic_force_motion_and_energy() {
    let mut joint = translated_residual(false, 1000.0);
    let mut reversed = translated_residual(true, 1000.0);
    let omega: [f64; 2] = [300.0, 500.0];
    let normal = [1.0, -0.8];
    let dt = joint.sample_period_s();
    // Stable closed-form unit-force displacement: (1-cos(w dt))/w^2.
    let response = omega.map(|w| 2.0*(0.5*w*dt).sin().powi(2)/(w*w));
    let compliance = normal.iter().zip(response).map(|(n, d)| n*n*d).sum::<f64>();
    let (mut q, mut v) = ([0.0; 2], [0.0; 2]);
    let (mut old_offset, mut old_reaction, mut old_energy) = (0.0, 0.0, 0.0);
    for command in [Some(-1e-5), Some(-2e-5), None, Some(-3e-5), Some(-2.5e-5), Some(-2e-5)] {
        let offset = command.unwrap_or(old_offset);
        let free_q: [f64; 2] = std::array::from_fn(|i| q[i]*(omega[i]*dt).cos()
            + v[i]*(omega[i]*dt).sin()/omega[i]);
        let free_v: [f64; 2] = std::array::from_fn(|i| v[i]*(omega[i]*dt).cos()
            - omega[i]*q[i]*(omega[i]*dt).sin());
        let old_x = normal[0]*q[0]+normal[1]*q[1]-PAIR_G*old_reaction-old_offset;
        let free_x = normal[0]*free_q[0]+normal[1]*free_q[1]-offset;
        // U=k*x^2/2 while compressed. Its secant gives R=k*(x0+x1)/2,
        // with x1=free_x-(N^T D N+G)R. No contact root or native report is
        // used to manufacture this force or either body's endpoint.
        let reaction = SUPPORT_N_M*(old_x+free_x)/(2.0+SUPPORT_N_M*(compliance+PAIR_G));
        q = std::array::from_fn(|i| free_q[i]-response[i]*normal[i]*reaction);
        v = std::array::from_fn(|i| free_v[i]-(omega[i]*dt).sin()/omega[i]*normal[i]*reaction);
        let x = normal[0]*q[0]+normal[1]*q[1]-PAIR_G*reaction-offset;
        assert!(old_x >= 0.0 && x > 0.0, "analytic linear-contact branch must stay compressed");
        let work = -reaction*(offset-old_offset);
        let tail_energy = 0.5*PAIR_G*reaction*reaction;
        let tail_loss = 0.5*PAIR_G*(reaction-old_reaction).powi(2);
        let contact_energy = 0.5*SUPPORT_N_M*x*x;
        let mechanical_energy = (0..2).map(|i| 0.5*(v[i]*v[i]+omega[i]*omega[i]*q[i]*q[i])).sum::<f64>();
        let energy = mechanical_energy+contact_energy+tail_energy;
        assert!((energy-old_energy+tail_loss-work).abs() < 2e-16);
        let step = |system: &mut MultiFrictionModalSystem| match command {
            Some(end) => system.step_with_obstacle_offsets(&[0.0, 0.0], &[released()], &[end]).unwrap().clone(),
            None => system.step(&[0.0, 0.0], &[released()]).unwrap().clone(),
        };
        let frame = step(&mut joint);
        assert_eq!(frame, step(&mut reversed), "binding order must not change the physical path");
        let stop = joint.last_obstacle_frames().unwrap()[0];
        assert!((stop.normal_force_n-reaction).abs() < 2e-12);
        for i in 0..2 {
            let state = joint.components()[i].states()[0];
            assert!((state.displacement_m_sqrt_kg-q[i]).abs() < 2e-15);
            assert!((state.velocity_m_sqrt_kg_per_s-v[i]).abs() < 2e-11);
        }
        assert!((frame.obstacle_work_j-work).abs() < 2e-16);
        assert!((stop.stored_energy_j-contact_energy).abs() < 2e-16);
        assert!((frame.obstacle_residual_energy_j-tail_energy).abs() < 2e-16);
        assert!((frame.obstacle_residual_loss_j-tail_loss).abs() < 2e-16);
        assert_eq!(frame.obstacle_dissipation_j, 0.0, "omitted-mode loss is not material heat");
        assert_eq!(stop.support_offset_m, offset);
        assert_eq!(joint.obstacle_offsets_m().unwrap(), &[offset]);
        assert!((joint.total_energy_j().unwrap()-energy).abs() < 2e-16);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        (old_offset, old_reaction, old_energy) = (offset, reaction, energy);
    }
    assert!(joint.components()[0].states()[0].velocity_m_sqrt_kg_per_s < 0.0);
    assert!(joint.components()[1].states()[0].velocity_m_sqrt_kg_per_s > 0.0);
}

#[test]
fn g0_translated_opposed_stops_keep_the_full_residual_matrix_and_each_actuator_work() {
    let gate = CancelGate::new_clock_free();
    let mut joint = pair(Some(PAIR_G), 0.0)
        .with_obstacle_motion(translation_limits(), &gate).unwrap();
    assert_eq!(joint.obstacle_pairs(), vec![(0, 1)]);
    let mut before = joint.total_energy_j().unwrap();
    let mut offsets = [0.0; 2];
    let mut previous_net = 0.0;
    for endpoint in [[-1e-5, -1e-5], [-2e-5, -1e-5], [-1.5e-5, 0.0], [-1e-5, 1e-5]] {
        let frame = joint.step_with_obstacle_offsets(&[0.0], &[released()], &endpoint).unwrap().clone();
        let stops = joint.last_obstacle_frames().unwrap();
        let net = stops[0].normal_force_n-stops[1].normal_force_n;
        let work = -(0..2).map(|i| stops[i].normal_force_n*(endpoint[i]-offsets[i])).sum::<f64>();
        let tail = 0.5*PAIR_G*net*net;
        let numerical = 0.5*PAIR_G*(net-previous_net).powi(2);
        let state = joint.components()[0].states()[0];
        let q = state.displacement_m_sqrt_kg;
        let energy = 0.5*(state.velocity_m_sqrt_kg_per_s.powi(2)+10000.0*q*q)
            + 0.5*SUPPORT_N_M*((q-endpoint[0]-PAIR_G*net+1e-4).powi(2)
                + (-q-endpoint[1]+PAIR_G*net+1e-4).powi(2)) + tail;
        assert!(stops.iter().all(|stop| stop.normal_force_n > 0.0));
        for i in 0..2 {
            assert_eq!(stops[i].support_offset_m, endpoint[i]);
            assert_eq!(stops[i].actuator_work_j, -stops[i].normal_force_n*(endpoint[i]-offsets[i]));
        }
        assert!((frame.obstacle_work_j-work).abs() < 1e-18);
        assert!((frame.obstacle_residual_energy_j-tail).abs() < 1e-18);
        assert!((frame.obstacle_residual_loss_j-numerical).abs() < 1e-18);
        assert!((joint.total_energy_j().unwrap()-energy).abs() < 1e-18);
        assert!((energy-before+frame.network_dissipation_j+numerical-work).abs() < 1e-13);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        (before, offsets, previous_net) = (energy, endpoint, net);
    }
}

#[test]
fn g4_refused_residual_motion_keeps_reactions_offsets_clock_and_retry_state() {
    let mut joint = translated_residual(false, 0.01);
    let mut reference = translated_residual(false, 0.01);
    for system in [&mut joint, &mut reference] {
        system.step_with_obstacle_offsets(&[0.0, 0.0], &[released()], &[-1e-6]).unwrap();
    }
    let before = states(&joint);
    let frame = joint.last_frame().unwrap().clone();
    let stops = joint.last_obstacle_frames().unwrap().to_vec();
    let energy = joint.total_energy_j().unwrap();
    // Motion is within its speed/offset limits. The candidate reaches the late
    // acceleration-pressure gate, after solving nonzero residual reactions.
    assert!(matches!(joint.step_with_obstacle_offsets(&[0.0, 0.0], &[released()], &[-2e-5]),
        Err(ModalFrictionError::Network(ModalCouplingError::Budget { .. }))));
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(joint.step_with_obstacle_offsets_under_gate(&[0.0, 0.0], &[released()], &[-2e-5], &cancelled).is_err());
    assert_eq!(states(&joint), before);
    assert_eq!(joint.last_frame(), Some(&frame));
    assert_eq!(joint.last_obstacle_frames().unwrap(), stops);
    assert_eq!(joint.obstacle_offsets_m().unwrap(), &[-1e-6]);
    assert_eq!(joint.total_energy_j().unwrap(), energy);
    assert_eq!(joint.samples_rendered(), 1);
    // Ordinary stepping holds the ACCEPTED offset and lag, not rejected scratch.
    assert_eq!(joint.step(&[0.0, 0.0], &[released()]).unwrap(),
        reference.step(&[0.0, 0.0], &[released()]).unwrap());
    assert_eq!(states(&joint), states(&reference));
    assert_eq!(joint.last_obstacle_frames(), reference.last_obstacle_frames());
    assert_eq!(joint.obstacle_offsets_m().unwrap(), &[-1e-6]);
}

#[test]
fn g5_zero_translation_preserves_the_existing_residual_trajectory_to_the_bit() {
    let gate = CancelGate::new_clock_free();
    let mut fixed = pair(Some(PAIR_G), 0.004);
    let mut translated = pair(Some(PAIR_G), 0.004)
        .with_obstacle_motion(translation_limits(), &gate).unwrap();
    for _ in 0..64 {
        assert_eq!(fixed.step(&[0.0], &[released()]).unwrap(),
            translated.step_with_obstacle_offsets(&[0.0], &[released()], &[0.0, 0.0]).unwrap());
        assert_eq!(states(&fixed), states(&translated));
        for (a, b) in fixed.components().iter().zip(translated.components()) {
            for (a, b) in a.states().iter().zip(b.states()) {
                assert_eq!(a.displacement_m_sqrt_kg.to_bits(), b.displacement_m_sqrt_kg.to_bits());
                assert_eq!(a.velocity_m_sqrt_kg_per_s.to_bits(), b.velocity_m_sqrt_kg_per_s.to_bits());
            }
        }
        assert_eq!(fixed.last_obstacle_frames(), translated.last_obstacle_frames());
    }
}
