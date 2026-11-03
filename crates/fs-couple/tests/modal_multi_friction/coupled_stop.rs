//! Known answers for a stop bank coupled to a friction port (MagentaBadger 12431).
//!
//! The oracle is the bare `CoupledModalSystem::step` under chosen held forces: a
//! port traction `t*` and a stop reaction `R*`. It is independent of the friction
//! assembly, the stops and the rescue, as in
//! `both_sticking_ports_recover_independently_manufactured_joint_reactions`. From
//! its endpoint the test sets the stop law (a touching unit stop with exponent
//! 1.5 and no loss, so from rest its reaction at closure p is exactly
//! `k p^1.5 / 2.5`, with `k` chosen to make that `R*`) and the drive (the end
//! speed `u = b . v_end`, the stick target). The bank must then reach `t*`, `R*`
//! and the known regime.
//!
//! | row | port | stop | known answer |
//! |---|---|---|---|
//! | A | stuck, prescribed N = 1 N | on the port's own body | t* 0.30 N < 0.8 N; R* 0.10 N |
//! | B | stuck, N = the stop's reaction | on the linked body, pressed by 0.8 N | t* 0.20 N < 0.8 R* = 0.40 N; R* 0.50 N |
//! | C | as B, surface 0.05 m/s faster than u | as B | sliding: t* = 0.4 R* = 0.20 N; R* 0.50 N |
//!
//! Each row runs with the full sweep (128 sweeps), which checks the oracle, and
//! with one sweep. One sweep solves the stop at zero traction, then the port,
//! whose traction then moves the stop through `nt`, so the step goes to the
//! active-set rescue. Before the condensed rescue a coupled stop declined there.
//!
//! The bound on t (MagentaBadger 12431: derived, not chosen). Write the stop's
//! closure `x = g + nt t - nn R`, its law `R = phi(x)`, `D = phi'(x*)`, and the
//! port's slip on the endpoint line `s = f - a t + T R` (`a`, `T`, `nt`, `nn` are
//! the network's own one-step responses, read below from the oracle with unit
//! loads). An accepted step has a stop residual `rho = R - phi(x)` with
//! `|rho| <= tau_R`, the stop's force tolerance. Linearizing about the oracle's
//! `(t*, R*)`:
//!   dR = (D nt dt + rho) / (1 + D nn).
//! A stuck row is accepted with `|s| <= eps_s`, the sticking tolerance, and
//! `s(t*, R*) = 0`, so `s = J dt + T rho / (1 + D nn)` with the condensed
//! derivative `J = -a + T D nt / (1 + D nn)`:
//!   |dt| <= (eps_s + |T| tau_R / (1 + D nn)) / |J|.
//! A sliding Coulomb row is accepted with `|t - mu_k R| <= tau_f`, the port's
//! force tolerance, so `dt (1 - mu_k D nt / (1 + D nn)) = sigma + mu_k rho / (1 + D nn)`:
//!   |dt| <= (tau_f + mu_k tau_R / (1 + D nn)) / |1 - mu_k D nt / (1 + D nn)|.
//! Then `|dR| <= (|D nt| |dt| + tau_R) / (1 + D nn)`. The oracle's own arithmetic
//! (its endpoint, `k` from `p`, the law's `pow`) adds rounding of a few ulps; the
//! bound allows `16 eps (t* + R*)` for it, about 3e-15 N. The second-order term
//! of the linearization is below 1e-20 N here. A measured error above this bound
//! is a finding, not a reason to move it.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::normal_load::ContactNormalLoad;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;

const MU_S: f64 = 0.8;
const MU_K: f64 = 0.4;
const ALPHA: f64 = 1.5;
const DT: f64 = 1.0 / 48000.0;

#[derive(Clone, Copy)]
struct Row {
    name: &'static str,
    stop_body: usize,
    contact_loaded: bool,
    traction: f64,
    reaction: f64,
    press: f64,
    sliding: bool,
}
const A: Row = Row { name: "A", stop_body: 0, contact_loaded: false, traction: 0.30, reaction: 0.10, press: 0.0, sliding: false };
const B: Row = Row { name: "B", stop_body: 1, contact_loaded: true, traction: 0.20, reaction: 0.50, press: 0.8, sliding: false };
const C: Row = Row { name: "C", stop_body: 1, contact_loaded: true, traction: MU_K * 0.50, reaction: 0.50, press: 0.8, sliding: true };

fn port_config() -> ModalFrictionConfig { config() }
fn stop_config() -> ModalContactConfig {
    ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0, maximum_penetration_m: 1e-3,
        force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-12 }
}
fn shape(body: usize) -> f64 { attachment(body).shapes[0] }
fn column(body: usize) -> [f64; 2] { let mut c = [0.0; 2]; c[body] = shape(body); c }
fn dot2(x: [f64; 2], y: [f64; 2]) -> f64 { x[0]*y[0] + x[1]*y[1] }

/// The bare network's endpoint from rest under modal forces: (q_end, v_end).
fn oracle(forces: [f64; 2]) -> ([f64; 2], [f64; 2]) {
    let mut network = network(true, 0.0);
    network.step(&forces).unwrap();
    let state = |i: usize| network.components()[i].states()[0];
    ([state(0).displacement_m_sqrt_kg, state(1).displacement_m_sqrt_kg],
        [state(0).velocity_m_sqrt_kg_per_s, state(1).velocity_m_sqrt_kg_per_s])
}

/// The row's manufactured stop stiffness, drive and external load.
fn manufacture(row: Row) -> (f64, FrictionDrive, [f64; 2], f64) {
    let (b, n) = (column(0), column(row.stop_body));
    let mut forces = [b[0]*row.traction, 0.0];
    forces[row.stop_body] -= n[row.stop_body]*row.reaction;
    let external = [0.0, shape(1)*row.press];
    forces[1] += external[1];
    let (q, v) = oracle(forces);
    let closure = dot2(n, q);
    assert!(closure > 0.0, "row {}: the oracle must close the stop", row.name);
    let stiffness = (ALPHA + 1.0)*row.reaction/closure.powf(ALPHA);
    let surface = dot2(b, v) + if row.sliding { 0.05 } else { 0.0 };
    let normal = if row.contact_loaded { 0.0 } else { 1.0 };
    (stiffness, drive(surface, normal), external, closure)
}

fn bank(row: Row, stiffness: f64, sweeps: usize) -> MultiFrictionModalSystem {
    let gate = CancelGate::new();
    let stop = ModalContact { left: attachment(row.stop_body),
        right: ModalAttachment { component: row.stop_body, shapes: vec![0.0] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], stiffness, ALPHA,
            "manufactured known-answer stop, not a measured material".into()).unwrap() };
    let bank = FrictionModalSystem::new(network(true, 0.0), ModalFrictionPort { left: attachment(0), right: None,
        law: FrictionLaw::Coulomb { static_mu: MU_S, kinetic_mu: MU_K } }, port_config(), &gate).unwrap()
        .with_additional_ports(vec![], MultiFrictionConfig { max_ports: 4, max_sweeps: sweeps, max_setup_terms: 100_000 }, &gate)
        .unwrap()
        .with_obstacles(vec![(stop, stop_config())], FrictionObstacleConfig { max_contacts: 4, max_setup_terms: 100_000 }, &gate)
        .unwrap();
    if row.contact_loaded {
        bank.with_contact_normal_loads(vec![ContactNormalLoad { friction_port: 0, obstacle: 0 }], &gate).unwrap()
    } else { bank }
}

/// The derived bounds on |t - t*| and |R - R*| (module comment).
fn bounds(row: Row, closure: f64) -> (f64, f64) {
    let (b, n) = (column(0), column(row.stop_body));
    let (qb, vb) = oracle(b);
    let (qn, vn) = oracle(n);
    let a = dot2(b, vb)*DT;
    let t_map = dot2(b, vn)*DT;
    let nt = dot2(n, qb);
    let nn = dot2(n, qn);
    let d = ALPHA*row.reaction/closure;
    let tau_r = stop_config().force_absolute_tolerance_n + stop_config().force_relative_tolerance*row.reaction;
    let tau_f = port_config().force_absolute_tolerance_n + port_config().force_relative_tolerance*row.traction;
    let eps_s = port_config().sticking_tolerance_m;
    let spread = 1.0 + d*nn;
    let rounding = 16.0*f64::EPSILON*(row.traction + row.reaction);
    let dt_bound = if row.sliding {
        (tau_f + MU_K*tau_r/spread)/(1.0 - MU_K*d*nt/spread).abs()
    } else {
        (eps_s + t_map.abs()*tau_r/spread)/(-a + t_map*d*nt/spread).abs()
    } + rounding;
    (dt_bound, ((d*nt).abs()*dt_bound + tau_r)/spread + rounding)
}

fn known_answer(row: Row, sweeps: usize) {
    let (stiffness, drive, external, closure) = manufacture(row);
    let mut bank = bank(row, stiffness, sweeps);
    let frame = bank.step(&external, &[drive])
        .unwrap_or_else(|error| panic!("row {} at {sweeps} sweeps refused: {error}", row.name)).clone();
    let port = frame.ports[0];
    let stop = bank.last_obstacle_frames().unwrap()[0];
    let (t_bound, r_bound) = bounds(row, closure);
    let (t_error, r_error) = ((port.traction_n - row.traction).abs(), (stop.normal_force_n - row.reaction).abs());
    eprintln!("KAT row {} sweeps {sweeps}: |t-t*| {t_error:.3e} N <= {t_bound:.3e}; |R-R*| {r_error:.3e} N <= {r_bound:.3e}; \
        sweeps {} rescues {} regime {:?}", row.name, frame.sweeps, frame.active_set_rescues, port.regime);
    let regime = if row.sliding { ModalFrictionRegime::Sliding } else { ModalFrictionRegime::Sticking };
    assert_eq!(port.regime, regime, "row {}", row.name);
    if sweeps == 1 {
        assert_eq!((frame.sweeps, frame.active_set_rescues), (1, 1), "row {}: one sweep goes to the rescue", row.name);
    } else {
        assert_eq!(frame.active_set_rescues, 0, "row {}: the full sweep needs no rescue", row.name);
    }
    assert!(t_error <= t_bound, "row {} at {sweeps} sweeps: |t-t*| {t_error:e} N exceeds the derived {t_bound:e} N", row.name);
    assert!(r_error <= r_bound, "row {} at {sweeps} sweeps: |R-R*| {r_error:e} N exceeds the derived {r_bound:e} N", row.name);
    if row.contact_loaded { assert_eq!(port.normal_force_n.to_bits(), stop.normal_force_n.to_bits()); }
    if !row.sliding {
        assert!(port.end_slip_speed_m_s.unwrap().abs()*DT <= port_config().sticking_tolerance_m, "row {}", row.name);
    } else {
        assert!(port.end_slip_speed_m_s.unwrap() > 0.0, "row {}: the surface outruns the contact", row.name);
    }
    assert!(stop.constitutive_residual_n.abs() <= stop.force_tolerance_n, "row {}", row.name);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j, "row {}", row.name);
}

#[test]
fn row_a_full_sweep_reaches_the_oracle() { known_answer(A, 128); }
#[test]
fn row_b_full_sweep_reaches_the_oracle() { known_answer(B, 128); }
#[test]
fn row_c_full_sweep_reaches_the_oracle() { known_answer(C, 128); }
#[test]
fn row_a_one_sweep_is_rescued_to_the_oracle() { known_answer(A, 1); }
#[test]
fn row_b_one_sweep_is_rescued_to_the_oracle() { known_answer(B, 1); }
#[test]
fn row_c_one_sweep_is_rescued_to_the_oracle() { known_answer(C, 1); }

#[test]
fn residual_stop_retains_three_participant_reactions_and_independent_work_balance() {
    let gate = CancelGate::new_clock_free();
    let omega: [f64; 3] = [800.0, 1100.0, 600.0];
    // String, physical contact receiver, then its distinct bridge body.
    let normal = [100.0, -60.0, 35.0];
    let tail = 2e-5;
    let stiffness = 50_000.0;
    let network = CoupledModalSystem::new(omega.iter().map(|w| image(*w, 0.0, 0.0)).collect(),
        vec![], budgets(), &gate).unwrap();
    let contact = ModalContact {
        left: ModalAttachment { component: 0, shapes: vec![normal[0]] },
        right: ModalAttachment { component: 1, shapes: vec![-normal[1]] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], stiffness, 1.0,
            "manufactured linear contact for reciprocal residual test".into()).unwrap(),
    };
    let mut bank = FrictionModalSystem::new(network, ModalFrictionPort {
        left: attachment(0), right: None, law: law(),
    }, config(), &gate).unwrap().with_additional_ports(vec![], joint_config(), &gate).unwrap()
        .with_obstacle_residual_ports_and_corrections(vec![(contact.clone(), stop_config())],
            vec![ModalAttachment { component: 2, shapes: vec![normal[2]] }], vec![tail],
            FrictionObstacleConfig { max_contacts: 4, max_setup_terms: 100_000 }, &gate).unwrap();
    assert_eq!(bank.obstacle(0).unwrap().right, contact.right);
    assert_eq!(bank.obstacle_normal_column(0).unwrap(), normal);

    let external = [1.0, 0.0, 0.0];
    let response: Vec<_> = omega.iter().map(|w| (1.0-(w*DT).cos())/(w*w)).collect();
    let free = normal[0]*response[0]*external[0];
    let compliance: f64 = normal.iter().zip(&response).map(|(n, d)| n*n*d).sum();
    // From rest, the linear contact's discrete-gradient force is k*x_end/2,
    // x_end = N^T D external - (N^T D N + G)*R. This is independent of
    // the normal solver and tests the bridge AND the physical receiver.
    let reaction = stiffness*free/(2.0+stiffness*(compliance+tail));
    let closure = 2.0*reaction/stiffness;
    let frame = bank.step(&external, &[drive(0.0, 0.0)]).unwrap().clone();
    let stop = bank.last_obstacle_frames().unwrap()[0];
    assert!((stop.normal_force_n-reaction).abs() < 2e-12);
    let mut mechanical_energy = 0.0;
    let mut external_work = 0.0;
    for i in 0..3 {
        let force = external[i]-normal[i]*reaction;
        let q = response[i]*force;
        let v = (omega[i]*DT).sin()/omega[i]*force;
        let actual = bank.components()[i].states()[0];
        assert!((actual.displacement_m_sqrt_kg-q).abs() < 2e-15, "coordinate {i}");
        assert!((actual.velocity_m_sqrt_kg_per_s-v).abs() < 2e-11, "reaction {i}");
        mechanical_energy += 0.5*(v*v+omega[i]*omega[i]*q*q);
        external_work += external[i]*q;
    }
    assert!(bank.components()[1].states()[0].velocity_m_sqrt_kg_per_s > 0.0);
    assert!(bank.components()[2].states()[0].velocity_m_sqrt_kg_per_s < 0.0);
    let contact_energy = 0.5*stiffness*closure*closure;
    let residual_energy = 0.5*tail*reaction*reaction;
    assert!((stop.stored_energy_j-contact_energy).abs() < 1e-16);
    assert!((stop.residual_energy_j-residual_energy).abs() < 1e-16);
    assert_eq!(stop.residual_loss_j, stop.residual_energy_j);
    assert!((mechanical_energy+contact_energy+2.0*residual_energy-external_work).abs() < 1e-16);
    assert!((bank.total_energy_j().unwrap()-mechanical_energy-contact_energy-residual_energy).abs() < 1e-16);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}

#[test]
fn residual_stop_refuses_incomplete_or_invalid_retained_corrections() {
    let gate = CancelGate::new_clock_free();
    for case in 0..4 {
        let contact = ModalContact {
            left: attachment(0), right: attachment(1),
            law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
                "synthetic correction admission".into()).unwrap(),
        };
        let mut corrections = vec![attachment(1)];
        match case {
            0 => corrections.clear(),
            1 => corrections[0].component = 2,
            2 => corrections[0].shapes.push(1.0),
            _ => corrections[0].shapes[0] = f64::NAN,
        }
        assert!(joint(false).with_obstacle_residual_ports_and_corrections(
            vec![(contact, stop_config())], corrections, vec![2e-5],
            FrictionObstacleConfig { max_contacts: 4, max_setup_terms: 100_000 }, &gate).is_err());
    }
}
