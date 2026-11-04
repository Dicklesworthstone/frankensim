//! Constant-loss connections on the prepared owner: each coupled mode carries the anchored loss factor times its
//! strain-energy share (the module notes), with the spring and every ledger unchanged.
use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::ConstantLossConnection;
use fs_math::c64::C64;
use fs_ports::PortedModule;

const TAU: f64 = core::f64::consts::TAU;
/// The maple pack's anchor, 2 pi 200 rad/s, and its shear loss factor.
const REFERENCE: f64 = TAU * 200.0;
const ETA: f64 = 0.0162;

fn config(modes: usize, links: usize) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: modes, max_connections: links,
        max_setup_terms: modes * (links + 1) * (links + 1),
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e8,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

fn budget(modes: usize) -> PreparedModalBudget {
    PreparedModalBudget { max_modes: modes, max_dense_bytes: 128 * 1024 * 1024, max_multiply_terms: 100_000_000_000 }
}

fn model(rate: u32, modes: &[(f64, f64)], states: &[(f64, f64)]) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(rate,
        modes.iter().map(|&(w, z)| ModalAcousticMode { angular_frequency_rad_s: w, damping_ratio: z,
            pressure_per_modal_velocity: C64::ZERO }).collect(),
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&states.iter().map(|&(q, v)| ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }).collect::<Vec<_>>()).unwrap();
    model
}

fn spectrum_shape(omega: f64) -> f64 {
    ConstantLossSpectrum::fit().shape(omega)
}

fn damping_column(system: &CoupledModalSystem, link: usize) -> Option<Vec<f64>> {
    system.exact.as_ref().unwrap().damping_column(link).map(<[f64]>::to_vec)
}

/// One mode tied to ground by a joint (k, c = eta k / omega_ref) with the joint holding half the strain energy.
fn single(rate: u32, omega: f64, anchored: bool) -> CoupledModalSystem {
    let (b, k) = (0.8, 0.5 * omega * omega / 0.64);
    let free = det::sqrt(omega * omega - k * b * b);
    let links = vec![ModalColumnConnection { column: vec![b], stiffness_n_m: k,
        damping_n_s_m: ETA * k / REFERENCE, rest_extension_m: 0.0 }];
    let anchors = [ConstantLossConnection { connection: 0, reference_omega_rad_s: REFERENCE }];
    CoupledModalSystem::new_with_columns_inner(vec![model(rate, &[(free, 0.0)], &[(1e-4, 0.0)])], links,
        config(1, 1), Some(budget(1)), if anchored { &anchors[..] } else { &[] }, &CancelGate::new_clock_free()).unwrap()
}

// The receipt of the law on one mode: its damping column is s b with s^2 = omega_ref G(omega) / omega, so its loss
// factor c d^2 / omega is eta x share x G(omega) at 200 Hz, 1 kHz and 5 kHz, against the Kelvin-Voigt eta x share x
// omega / omega_ref.
#[test]
fn a_single_anchored_joint_gives_its_mode_the_constant_loss_factor_at_every_frequency() {
    for hz in [200.0, 1000.0, 5000.0] {
        let omega = TAU * hz;
        let system = single(192_000, omega, true);
        let link = &system.connections[0];
        let share = link.stiffness_n_m * system.columns[0][0].powi(2) / (omega * omega);
        let d = damping_column(&system, 0).unwrap();
        let loss = link.damping_n_s_m * d[0] * d[0] / omega;
        let expected = ETA * share * spectrum_shape(omega);
        eprintln!("{hz} Hz: share {share:.6}, loss factor {loss:.6e}, expected {expected:.6e}, Kelvin-Voigt {:.6e}",
            ETA * share * omega / REFERENCE);
        assert!((loss / expected - 1.0).abs() < 1e-10, "{hz} Hz: {loss:e} != {expected:e}");
        assert!(damping_column(&single(192_000, omega, false), 0).is_none());
    }
}

// The free decay of that mode at 1 kHz over 50 periods, read in total energy at whole periods: exp(-2 zeta omega t) with
// 2 zeta = eta x share x G, to 1e-4. The transition is exact, so this measures the law and not a time discretization.
#[test]
fn the_anchored_mode_decays_at_its_constant_loss_rate() {
    let omega = TAU * 1000.0;
    let mut system = single(48_000, omega, true);
    let zeta = 0.5 * ETA * 0.5 * spectrum_shape(omega);
    let start = system.total_energy_j().unwrap();
    for period in 1..=50 {
        for _ in 0..48 {
            system.step(&[0.0]).unwrap();
        }
        let ratio = system.total_energy_j().unwrap() / start;
        let expected = (-2.0 * zeta * omega * period as f64 / 1000.0).exp();
        assert!((ratio / expected - 1.0).abs() < 1e-4, "period {period}: {ratio:e} != {expected:e}");
    }
}

fn three(rate: u32, anchored: bool) -> CoupledModalSystem {
    let modes = [(TAU * 180.0, 0.002), (TAU * 950.0, 0.004), (TAU * 3100.0, 0.003)];
    let states = [(2e-5, 0.01), (-1e-5, -0.02), (4e-6, 0.03)];
    let links = vec![
        ModalColumnConnection { column: vec![1.0, -0.6, 0.3], stiffness_n_m: 4e6,
            damping_n_s_m: ETA * 4e6 / REFERENCE, rest_extension_m: 1e-7 },
        ModalColumnConnection { column: vec![-0.2, 0.7, 0.9], stiffness_n_m: 2e7,
            damping_n_s_m: 30.0, rest_extension_m: -1e-7 },
    ];
    let anchors = [ConstantLossConnection { connection: 0, reference_omega_rad_s: REFERENCE }];
    CoupledModalSystem::new_with_columns_inner(vec![model(rate, &modes, &states)], links, config(3, 2),
        Some(budget(3)), if anchored { &anchors[..] } else { &[] }, &CancelGate::new_clock_free()).unwrap()
}

// Three modes and two noncommuting joints, one anchored. An independent Jacobi eigensolve of the coupled stiffness
// gives every coupled mode a: the anchored dashpot's modal loss factor c (v_a.d)^2 / omega_a is eta x share_a x
// G(omega_a), share_a = k (v_a.b)^2 / omega_a^2, and its modal coupling c (v_a.d)(v_b.d) is
// c s_a s_b (v_a.b)(v_b.b). The Kelvin-Voigt joint keeps its own column.
#[test]
fn every_coupled_mode_carries_the_anchored_loss_factor_times_its_strain_energy_share() {
    let system = three(192_000, true);
    assert!(damping_column(&system, 1).is_none());
    let d = damping_column(&system, 0).unwrap();
    let (link, b) = (&system.connections[0], &system.columns[0]);
    let n = 3;
    let mut k = vec![0.0; n * n];
    for (i, mode) in system.models[0].modes().iter().enumerate() {
        k[i * n + i] = mode.angular_frequency_rad_s.powi(2);
    }
    for (column, law) in system.columns.iter().zip(&system.connections) {
        for i in 0..n {
            for j in 0..n {
                k[i * n + j] += law.stiffness_n_m * column[i] * column[j];
            }
        }
    }
    let (values, vectors) = fs_la::eigen::jacobi_eigh(&k, n);
    let project = |a: usize, x: &[f64]| (0..n).map(|i| vectors[i * n + a] * x[i]).sum::<f64>();
    let s = |a: usize| det::sqrt(REFERENCE * spectrum_shape(det::sqrt(values[a])) / det::sqrt(values[a]));
    for a in 0..n {
        let omega = det::sqrt(values[a]);
        let share = link.stiffness_n_m * project(a, b).powi(2) / values[a];
        let loss = link.damping_n_s_m * project(a, &d).powi(2) / omega;
        let expected = ETA * share * spectrum_shape(omega);
        eprintln!("coupled mode {a}: {:.2} Hz, share {share:.4e}, loss factor {loss:.6e}, expected {expected:.6e}",
            omega / TAU);
        assert!((loss / expected - 1.0).abs() < 1e-9, "mode {a}: {loss:e} != {expected:e}");
        for c in 0..n {
            let coupling = link.damping_n_s_m * project(a, &d) * project(c, &d);
            let expected = link.damping_n_s_m * s(a) * s(c) * project(a, b) * project(c, b);
            let scale = link.damping_n_s_m * s(a) * s(c) * b.iter().map(|x| x * x).sum::<f64>();
            assert!((coupling - expected).abs() <= 1e-9 * scale, "({a}, {c}): {coupling:e} != {expected:e}");
        }
    }
}

// The prepared energy identity and the per-port ledger hold with an anchored joint: every accepted step's port losses
// sum to its connection loss (the ledger refuses otherwise), its audit is clean, and the anchored network loses less
// to the joint than the Kelvin-Voigt one above the anchor.
#[test]
fn the_anchored_network_keeps_its_energy_identity_and_its_per_port_ledger() {
    let mut anchored = three(48_000, true);
    let mut viscous = three(48_000, false);
    let (mut lost, mut lost_viscous) = (0.0, 0.0);
    for tick in 0..240 {
        let external = [0.002, if tick % 2 == 0 { -0.003 } else { 0.001 }, 0.0005];
        let frame = anchored.step(&external).unwrap().clone();
        lost += frame.connection_dissipation_j;
        lost_viscous += viscous.step(&external).unwrap().connection_dissipation_j;
        let transfers = anchored.last_window_transfers().unwrap();
        assert_eq!(transfers.len(), 2);
        assert!(anchored.last_window_transfer_error().is_none());
        let fs_ports::ModuleWindowAudit::Available(audit) = anchored.audit_last_window().unwrap() else {
            panic!("anchored prepared physical account unavailable")
        };
        assert!(audit.violations().is_empty());
    }
    eprintln!("joint loss over 240 samples: anchored {lost:e} J, Kelvin-Voigt {lost_viscous:e} J");
    assert!(lost > 0.0 && lost < lost_viscous);
}

// The cold binder is the constructor: binding the anchor to the built prepared owner gives the same steps bit for bit.
// It refuses a started, an unprepared or an already anchored owner, and an invalid list.
#[test]
fn the_cold_binder_equals_the_constructor_and_refuses_what_it_cannot_honour() {
    let gate = CancelGate::new_clock_free();
    let anchor = ConstantLossConnection { connection: 0, reference_omega_rad_s: REFERENCE };
    let mut built = three(48_000, true);
    let mut bound = three(48_000, false).with_constant_loss_connections(&[anchor], &gate).unwrap();
    assert_eq!(bound.constant_loss_connections(), vec![anchor]);
    assert_eq!(built.constant_loss_connections(), vec![anchor]);
    for tick in 0..32 {
        let external = [0.001 * tick as f64, -0.002, 0.0];
        assert_eq!(built.step(&external).unwrap().clone(), bound.step(&external).unwrap().clone());
        assert_eq!(flatten(&built.models), flatten(&bound.models));
    }
    assert!(built.with_constant_loss_connections(&[anchor], &gate).is_err());
    assert!(three(48_000, true).with_constant_loss_connections(&[anchor], &gate).is_err());
    for invalid in [
        ConstantLossConnection { connection: 2, reference_omega_rad_s: REFERENCE },
        ConstantLossConnection { connection: 0, reference_omega_rad_s: 0.0 },
        ConstantLossConnection { connection: 0, reference_omega_rad_s: f64::NAN },
    ] {
        assert!(three(48_000, false).with_constant_loss_connections(&[invalid], &gate).is_err());
    }
    assert!(three(48_000, false).with_constant_loss_connections(&[anchor, anchor], &gate).is_err());
    assert!(three(48_000, false).with_constant_loss_connections(&[], &gate).is_err());
    let unprepared = CoupledModalSystem::new_with_columns(
        vec![model(48_000, &[(TAU * 180.0, 0.002)], &[(0.0, 0.0)])],
        vec![ModalColumnConnection { column: vec![1.0], stiffness_n_m: 1e6, damping_n_s_m: 1.0,
            rest_extension_m: 0.0 }],
        config(1, 1), &gate).unwrap();
    assert!(unprepared.with_constant_loss_connections(&[anchor], &gate).is_err());
    assert!(CoupledModalSystem::new_with_columns_inner(
        vec![model(48_000, &[(TAU * 180.0, 0.002)], &[(0.0, 0.0)])],
        vec![ModalColumnConnection { column: vec![1.0], stiffness_n_m: 1e6, damping_n_s_m: 1.0,
            rest_extension_m: 0.0 }],
        config(1, 1), None, &[anchor], &gate).is_err());
    // A spring without a dashpot has nothing to anchor.
    assert!(CoupledModalSystem::new_with_columns_inner(
        vec![model(48_000, &[(TAU * 180.0, 0.002)], &[(0.0, 0.0)])],
        vec![ModalColumnConnection { column: vec![1.0], stiffness_n_m: 1e6, damping_n_s_m: 0.0,
            rest_extension_m: 0.0 }],
        config(1, 1), Some(budget(1)), &[anchor], &gate).is_err());
}
