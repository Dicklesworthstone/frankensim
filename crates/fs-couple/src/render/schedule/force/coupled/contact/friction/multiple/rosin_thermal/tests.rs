//! Smith & Woodhouse thermal rosin in the joint owner, on the paper's own rig:
//! the section 5.1 Perspex rod and wedge at 3 N, a = 0.5 mm, delta = 10 um,
//! Table 1 rosin, and the Eq. (2) fit to Fig. 4 (fs-tribo's g1 constants).
//! A heavy single mode under each port keeps the slip speed at the drive speed.

use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_math::c64::C64;
use fs_math::det;
use fs_tribo::rosin_thermal::{RosinContact, SteadySlidingCurve, SteadySlidingFit, ThermalSolid};

const RATE: u32 = 48_000;
const RIG_N: f64 = 3.0;
const RADIUS_M: f64 = 0.5e-3;
const PERSPEX: ThermalSolid = ThermalSolid {
    conductivity_w_per_m_k: 0.2,
    density_kg_per_m3: 1200.0,
    specific_heat_j_per_kg_k: 1460.0,
};
const FIT: SteadySlidingFit = SteadySlidingFit { amplitudes: [0.4, 0.45], decay_speeds_mps: [0.01, 0.1], asymptote: 0.35 };

/// Eq. (2) written out independently of the law.
fn eq2(v: f64) -> f64 {
    0.4 * det::exp(-v / 0.01) + 0.45 * det::exp(-v / 0.1) + 0.35
}

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 4,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    }
}
fn config() -> ModalFrictionConfig {
    ModalFrictionConfig {
        max_iterations: 96,
        maximum_force_n: 50.0,
        maximum_slip_speed_m_s: 100.0,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    }
}
fn model(hz: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        RATE,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: core::f64::consts::TAU * hz,
            damping_ratio: 0.05,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
/// One 100 kg mode under each of `cells` ports, every port on `law`.
fn system(cells: usize, law: FrictionLaw) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    #[allow(clippy::cast_precision_loss)]
    let models = (0..cells).map(|k| model(180.0 + 20.0 * k as f64)).collect();
    let network = CoupledModalSystem::new(models, vec![], budgets(), &gate).unwrap();
    let port = |component| ModalFrictionPort {
        left: ModalAttachment { component, shapes: vec![0.1] },
        right: None,
        law: law.clone(),
    };
    FrictionModalSystem::new(network, port(0), config(), &gate)
        .unwrap()
        .with_additional_ports(
            (1..cells).map(|c| (port(c), config())).collect(),
            MultiFrictionConfig { max_ports: 4, max_sweeps: 128, max_setup_terms: 1_000_000 },
            &gate,
        )
        .unwrap()
}
fn coulomb() -> FrictionLaw {
    FrictionLaw::Coulomb { static_mu: 1.2, kinetic_mu: 1.2 }
}
/// The rig as the WHOLE contact, on `system`'s own clock.
fn params(system: &MultiFrictionModalSystem, history_capacity: usize, speed_max: f64) -> SmithWoodhouseParams {
    let rig = RosinContact {
        radius_m: RADIUS_M,
        reference_normal_n: RIG_N,
        layer_thickness_m: 10e-6,
        stationary: PERSPEX,
        moving: PERSPEX,
    };
    SmithWoodhouseParams {
        rosin_density_kg_per_m3: 1010.0,
        rosin_specific_heat_j_per_kg_k: 250.0,
        steady_fit: SteadySlidingCurve::Exponential(FIT),
        calibration: rig,
        calibration_speed_max_mps: speed_max,
        contact: rig,
        step_s: system.sample_period_s(),
        history_capacity,
    }
}
fn sliding(normal: f64, slip: f64) -> FrictionPortFrame {
    FrictionPortFrame {
        normal_force_n: normal,
        slip_distance_m: slip,
        regime: if slip == 0.0 { ModalFrictionRegime::Sticking } else { ModalFrictionRegime::Sliding },
        ..FrictionPortFrame::default()
    }
}
fn frame(ports: Vec<FrictionPortFrame>) -> MultiFrictionFrame {
    MultiFrictionFrame { ports, ..MultiFrictionFrame::default() }
}
fn released() -> FrictionPortFrame {
    FrictionPortFrame { regime: ModalFrictionRegime::Released, ..FrictionPortFrame::default() }
}
/// One sample through the joint owner of a one-port `sys`: the step-start
/// coefficient, the step with the bow at `normal` and `v`, then the heat.
fn stroke(sys: &mut MultiFrictionModalSystem, rosin: &mut RosinThermalPorts, normal: f64, v: f64)
    -> (MultiFrictionFrame, Result<(), ModalFrictionError>)
{
    rosin.apply(sys).unwrap();
    let f = sys.step(&vec![0.0], &vec![FrictionDrive { speed_m_s: v, normal_force_n: normal }]).unwrap().clone();
    let accepted = rosin.accept(&f);
    (f, accepted)
}
/// Land at the first sample, slide `loaded` samples at the rig load, then lift
/// for `off` samples. Returns the rosin and the temperature at the lift.
fn land_and_lift(sys: &mut MultiFrictionModalSystem, loaded: usize, off: usize) -> RosinThermalPorts {
    let mut rosin = RosinThermalPorts::new(params(sys, loaded, 0.3), vec![0], sys).unwrap();
    for _ in 0..loaded { stroke(sys, &mut rosin, RIG_N, 0.1).1.unwrap(); }
    for _ in 0..off { stroke(sys, &mut rosin, 0.0, 0.1).1.unwrap(); }
    rosin
}

#[test]
fn bowed_port_at_constant_speed_returns_the_calibration_curve() {
    // Steady sliding through the joint owner. The stationary half-space keeps
    // drawing K T/sqrt(pi D t), so mu(t) = eq2(v) + c t^{-1/2} + O(t^{-1}), as
    // fs-tribo's g1 steady-sliding test: the deviation halves over each 4x in
    // time to within sqrt(t_c/t), t_c = 8a/(3 pi v), and one Richardson step in
    // t^{-1/2} removes at least 90 % of the deviation left at 1 s. Every step's
    // traction is mu N from the step-start temperature, and every temperature is
    // an independent law's, driven with the frames' own slips, to the bit.
    for (cells, v) in [(1_usize, 0.05), (1, 0.1), (1, 0.2), (2, 0.1)] {
        let mut sys = system(cells, coulomb());
        let marks = [RATE as usize / 16, RATE as usize / 4, RATE as usize];
        let mut rosin = RosinThermalPorts::new(params(&sys, marks[2], 0.3), (0..cells).collect(), &mut sys).unwrap();
        let single = SmithWoodhouse::new(params(&sys, marks[2], 0.3)).unwrap();
        let mut independent = single.initial_state();
        let normal = RIG_N / f64::from(u32::try_from(cells).unwrap());
        let drives = vec![FrictionDrive { speed_m_s: v, normal_force_n: normal }; cells];
        let area = core::f64::consts::PI * RADIUS_M * RADIUS_M;
        let mut mu = Vec::with_capacity(marks[2]);
        for _ in 0..marks[2] {
            rosin.apply(&mut sys).unwrap();
            let expected = area / RIG_N * single.yield_stress_pa(independent.temperature_rise_k()).unwrap();
            let f = sys.step(&vec![0.0; cells], &drives).unwrap().clone();
            rosin.accept(&f).unwrap();
            for (i, port) in f.ports.iter().enumerate() {
                assert_eq!(port.regime, ModalFrictionRegime::Sliding, "cell {i} slides throughout");
                assert!((port.traction_n.abs() - rosin.coefficients()[i] * normal).abs() <= 1e-9 * normal,
                    "cell {i}: |F| = mu N");
            }
            assert_eq!(rosin.coefficients()[0].to_bits(), expected.to_bits(), "step-start coefficient");
            let input = TangentialStepInput { normal_n: RIG_N, slip_increment_m: [f.ports[0].slip_distance_m, 0.0],
                dt_s: single.params().step_s, contact_temperature_k: None, nominal_pressure_pa: None,
                compliance_m_per_n: 0.0,
                solve: TangentialSolveConfig { max_iterations: 1, maximum_force_n: f64::MAX,
                    force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-12 } };
            independent = single.advance(&independent, &input).unwrap().state;
            if cells == 1 {
                assert_eq!(rosin.temperatures_k()[0].to_bits(), independent.temperature_rise_k().to_bits());
            }
            mu.push(f.ports[0].traction_n.abs() / normal);
        }
        let transit = 8.0 * RADIUS_M / (3.0 * core::f64::consts::PI * v);
        #[allow(clippy::cast_precision_loss)]
        let times: Vec<f64> = marks.iter().map(|&n| n as f64 / f64::from(RATE)).collect();
        let dev: Vec<f64> = marks.iter().map(|&n| mu[n - 1] - eq2(v)).collect();
        let (r1, r2) = (det::sqrt(times[1]), det::sqrt(times[2]));
        let extrapolated = (r2 * mu[marks[2] - 1] - r1 * mu[marks[1] - 1]) / (r2 - r1);
        let ratio = dev[2] / dev[1];
        println!("{{\"test\":\"rosin_constant_speed\",\"cells\":{cells},\"v_mps\":{v},\"mu\":[{},{},{}],\"eq2\":{},\"ratio\":{ratio},\"extrapolated\":{extrapolated},\"heat_s\":{}}}",
            mu[marks[0] - 1], mu[marks[1] - 1], mu[marks[2] - 1], eq2(v), rosin.heat_seconds());
        assert!(dev[0] > dev[1] && dev[1] > dev[2] && dev[2] > 0.0, "approach from above: {dev:?}");
        assert!((ratio - 0.5).abs() <= det::sqrt(transit / times[1]), "t^-1/2 law: {ratio}");
        assert!((extrapolated - eq2(v)).abs() <= 0.1 * dev[2], "Richardson limit {extrapolated} vs Eq. (2) {}", eq2(v));
    }
}

#[test]
fn one_and_three_cells_give_the_same_temperature_under_uniform_slip() {
    // The cells cut the contact along the string, and a hair point crosses all
    // of it in every cell: refining the ribbon must not change the physics.
    // Fig. 18's forced sinusoid (4 Hz, 0.08 m/s peak) passes through stick.
    let mut one = system(1, coulomb());
    let mut three = system(3, coulomb());
    let steps = RATE as usize / 4;
    let mut a = RosinThermalPorts::new(params(&one, steps, 0.3), vec![0], &mut one).unwrap();
    let mut b = RosinThermalPorts::new(params(&three, steps, 0.3), vec![0, 1, 2], &mut three).unwrap();
    let h = one.sample_period_s();
    let mut stuck = 0;
    for k in 0..steps {
        #[allow(clippy::cast_precision_loss)]
        let t = (k as f64 + 0.5) * h;
        let slip = if (k / 600) % 3 == 2 { 0.0 } else { 0.08 * det::sin(core::f64::consts::TAU * 4.0 * t) * h };
        if slip == 0.0 { stuck += 1; }
        a.apply(&mut one).unwrap();
        b.apply(&mut three).unwrap();
        for c in b.coefficients() {
            assert_eq!(c.to_bits(), a.coefficients()[0].to_bits(), "step {k}: coefficient");
        }
        a.accept(&frame(vec![sliding(RIG_N, slip)])).unwrap();
        b.accept(&frame(vec![sliding(1.0, slip); 3])).unwrap();
        for t3 in b.temperatures_k() {
            assert_eq!(t3.to_bits(), a.temperatures_k()[0].to_bits(), "step {k}: temperature");
        }
    }
    assert!(stuck > 0 && a.temperatures_k()[0] > 0.0);
    assert_eq!((a.accepted_steps(), a.heat_steps()), (steps, vec![steps]));
}

#[test]
fn a_landing_at_step_k_equals_the_same_contact_loaded_from_k_with_the_prefix_released() {
    // The history's t = 0 is the cell's first loaded step. A: the joint owner
    // with the bow released for k samples, then loaded. B: a fresh history fed
    // A's own loaded frames from its first step. Every coefficient and
    // temperature must agree to the bit. A's capacity is exactly its loaded
    // steps, so a released prefix that took history would refuse.
    let (k, n, v) = (240, 960, 0.1);
    let mut sys_a = system(1, coulomb());
    let mut sys_b = system(1, coulomb());
    let mut a = RosinThermalPorts::new(params(&sys_a, n, 0.3), vec![0], &mut sys_a).unwrap();
    let mut b = RosinThermalPorts::new(params(&sys_b, n, 0.3), vec![0], &mut sys_b).unwrap();
    let ambient = a.law().initial_state().temperature_rise_k();
    for j in 0..k {
        let (f, accepted) = stroke(&mut sys_a, &mut a, 0.0, v);
        accepted.unwrap();
        assert_eq!((f.ports[0].regime, f.ports[0].traction_n), (ModalFrictionRegime::Released, 0.0), "prefix step {j}");
    }
    assert_eq!(a.phases(), [RosinCellPhase::BeforeLanding]);
    assert_eq!((a.accepted_steps(), a.heat_steps()), (k, vec![0]));
    assert_eq!(a.temperatures_k()[0].to_bits(), ambient.to_bits());
    for j in 0..n {
        let (f, accepted) = stroke(&mut sys_a, &mut a, RIG_N, v);
        accepted.unwrap();
        assert_eq!(f.ports[0].regime, ModalFrictionRegime::Sliding, "loaded step {j}");
        b.apply(&mut sys_b).unwrap();
        b.accept(&f).unwrap();
        assert_eq!(a.coefficients()[0].to_bits(), b.coefficients()[0].to_bits(), "loaded step {j}: coefficient");
        assert_eq!(a.temperatures_k()[0].to_bits(), b.temperatures_k()[0].to_bits(), "loaded step {j}: temperature");
    }
    assert_eq!((a.phases(), b.phases()), ([RosinCellPhase::InContact { landed: k }].as_slice(),
        [RosinCellPhase::InContact { landed: 0 }].as_slice()));
    assert_eq!((a.accepted_steps(), a.heat_steps(), b.heat_steps()), (k + n, vec![n], vec![n]));
    assert!(a.temperatures_k()[0] > ambient);
    // Cells land on their own: cell 1, loaded one step after cell 0, starts
    // its history at its own first contact.
    let mut two = system(2, coulomb());
    let mut ribbon = RosinThermalPorts::new(params(&two, 4, 0.3), vec![0, 1], &mut two).unwrap();
    let mut fresh = RosinThermalPorts::new(params(&two, 4, 0.3), vec![0, 1], &mut two).unwrap();
    let h = two.sample_period_s();
    ribbon.accept(&frame(vec![sliding(1.5, 0.2 * h), released()])).unwrap();
    assert_eq!(ribbon.phases(), [RosinCellPhase::InContact { landed: 0 }, RosinCellPhase::BeforeLanding]);
    ribbon.accept(&frame(vec![sliding(1.5, 0.2 * h), sliding(1.5, 0.1 * h)])).unwrap();
    fresh.accept(&frame(vec![sliding(1.5, 0.1 * h); 2])).unwrap();
    assert_eq!(ribbon.phases(), [RosinCellPhase::InContact { landed: 0 }, RosinCellPhase::InContact { landed: 1 }]);
    assert_eq!(ribbon.heat_steps(), vec![2, 1]);
    assert_eq!(ribbon.temperatures_k()[1].to_bits(), fresh.temperatures_k()[1].to_bits());
}

#[test]
fn after_the_lift_the_force_is_zero_and_the_history_is_frozen() {
    // Lifted, the cell adds nothing: the owner releases the port, its traction
    // is zero, and its temperature and heat steps stay as they were at the
    // lift, to the bit, however long the bow is off. Every step off the
    // string takes the frozen temperature's coefficient, pi a^2 k_y(T)/N_ref.
    let (loaded, off) = (480, 720);
    let mut sys = system(1, coulomb());
    let mut rosin = land_and_lift(&mut sys, loaded, 0);
    let temperature = rosin.temperatures_k()[0];
    assert!(temperature > rosin.law().initial_state().temperature_rise_k());
    let area = core::f64::consts::PI * RADIUS_M * RADIUS_M;
    let coefficient = area / RIG_N * rosin.law().yield_stress_pa(temperature).unwrap();
    for j in 0..off {
        let (f, accepted) = stroke(&mut sys, &mut rosin, 0.0, 0.1);
        accepted.unwrap();
        assert_eq!((f.ports[0].regime, f.ports[0].traction_n), (ModalFrictionRegime::Released, 0.0), "off step {j}");
        assert_eq!(rosin.temperatures_k()[0].to_bits(), temperature.to_bits(), "off step {j}: temperature");
        assert_eq!(rosin.coefficients()[0].to_bits(), coefficient.to_bits(), "off step {j}: coefficient");
    }
    assert_eq!(rosin.phases(), [RosinCellPhase::Lifted { landed: 0, lifted: loaded }]);
    assert_eq!((rosin.accepted_steps(), rosin.heat_steps()), (loaded + off, vec![loaded]));
}

#[test]
fn a_second_landing_refuses_by_name_and_commits_nothing() {
    // Off-string cooling is unsourced, so the law cannot take the bow back.
    let mut sys = system(1, coulomb());
    let mut rosin = land_and_lift(&mut sys, 240, 240);
    let before = (rosin.temperatures_k(), rosin.phases().to_vec(), rosin.accepted_steps(), rosin.heat_steps());
    let (_, accepted) = stroke(&mut sys, &mut rosin, RIG_N, 0.1);
    let text = accepted.unwrap_err().to_string();
    assert!(text.contains("lands again after its lift"), "{text}");
    assert_eq!((rosin.temperatures_k(), rosin.phases().to_vec(), rosin.accepted_steps(), rosin.heat_steps()), before);
    assert_eq!(rosin.phases(), [RosinCellPhase::Lifted { landed: 0, lifted: 240 }]);
}

#[test]
fn refusals_are_named_and_leave_the_histories_unchanged() {
    let message = |e: ModalFrictionError| e.to_string();
    // A non-Coulomb port keeps its own law.
    let mut stribeck = system(1, FrictionLaw::Stribeck { static_mu: 1.2, kinetic_mu: 0.35,
        characteristic_speed: 0.1, viscous_per_speed: 0.0 });
    let p = params(&stribeck, 8, 0.3);
    assert!(message(RosinThermalPorts::new(p, vec![0], &mut stribeck).unwrap_err()).contains("plain Coulomb port"));
    // Ports must exist and be distinct.
    let mut two = system(2, coulomb());
    for ports in [vec![], vec![0, 0], vec![2]] {
        let p = params(&two, 8, 0.3);
        assert!(message(RosinThermalPorts::new(p, ports, &mut two).unwrap_err()).contains("distinct existing ports"));
    }
    // The law's step is the network's sample period.
    let mut p = params(&two, 8, 0.3);
    p.step_s *= 2.0;
    assert!(message(RosinThermalPorts::new(p, vec![0, 1], &mut two).unwrap_err()).contains("differs from the network sample period"));
    // A normal load that is negative or not finite, and an unloaded cell that
    // still carries traction, refuse; nothing commits, the loaded cell included.
    let mut rosin = RosinThermalPorts::new(params(&two, 8, 0.3), vec![0, 1], &mut two).unwrap();
    let h = two.sample_period_s();
    rosin.accept(&frame(vec![sliding(1.5, 0.2 * h); 2])).unwrap();
    let before = (rosin.temperatures_k(), rosin.phases().to_vec(), rosin.accepted_steps(), rosin.heat_steps());
    let dragging = FrictionPortFrame { traction_n: 1e-3, ..released() };
    for (bad, text) in [(sliding(-1.0, 0.2 * h), "finite nonnegative normal load"),
                        (sliding(f64::NAN, 0.2 * h), "finite nonnegative normal load"),
                        (dragging, "unloaded rosin thermal cell carries traction")] {
        let e = rosin.accept(&frame(vec![sliding(1.5, 0.2 * h), bad])).unwrap_err();
        assert!(message(e).contains(text), "{text}");
        assert_eq!((rosin.temperatures_k(), rosin.phases().to_vec(), rosin.accepted_steps(), rosin.heat_steps()), before);
    }
    // The history capacity.
    let mut capped = RosinThermalPorts::new(params(&two, 2, 0.3), vec![0, 1], &mut two).unwrap();
    for _ in 0..2 { capped.accept(&frame(vec![sliding(1.5, 0.2 * h); 2])).unwrap(); }
    let e = capped.accept(&frame(vec![sliding(1.5, 0.2 * h); 2])).unwrap_err();
    assert!(message(e).contains("smith_woodhouse_history_steps"));
    // An endpoint above the calibration refuses before publishing either
    // cell. Cell 0 remains inside the calibrated range and stages first;
    // cell 1 slides at 0.3 m/s beyond this 2 mm/s calibration's heat range.
    let mut hot = RosinThermalPorts::new(params(&two, 4_000, 0.002), vec![0, 1], &mut two).unwrap();
    let limit = hot.law().temperature_limit_k();
    let snapshot = |r: &RosinThermalPorts| (
        r.states.clone(), r.phases().to_vec(), r.accepted_steps(),
        r.coefficients().to_vec(), r.heat_seconds(),
    );
    let mut refused = None;
    for k in 0..4_000 {
        let before = snapshot(&hot);
        if let Err(e) = hot.accept(&frame(vec![sliding(1.5, 0.001 * h), sliding(1.5, 0.3 * h)])) {
            assert_eq!(snapshot(&hot), before, "refused step {k} must commit no cell");
            refused = Some((k, message(e)));
            break;
        }
        hot.apply(&mut two).unwrap();
    }
    let (k, text) = refused.expect("the temperature limit is reached");
    assert!(text.contains("smith_woodhouse_temperature_rise"), "{text}");
    assert!(hot.temperatures_k().iter().all(|&t| t <= limit), "step {k}");
    assert_eq!((hot.accepted_steps(), hot.heat_steps()), (k, vec![k; 2]));
    // The last valid checkpoint still defines both coefficients and can
    // continue loaded cooling, without resetting either thermal history.
    let before = snapshot(&hot);
    hot.apply(&mut two).unwrap();
    assert_eq!(snapshot(&hot), before);
    hot.accept(&frame(vec![sliding(1.5, 0.0); 2])).unwrap();
    assert_eq!(hot.heat_steps(), vec![k + 1; 2]);
}
