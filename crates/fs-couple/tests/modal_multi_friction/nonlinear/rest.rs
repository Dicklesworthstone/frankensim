//! A reaction iteration whose joint sweep stops at a different sweep from one trial to the next.
//!
//! Three bowed ports hold the three modes of one component. Its first mode carries the manufactured
//! quartic stretch energy of the cone-edge tests and starts displaced. The ports have the bowed
//! card's tolerances. The joint sweep starts each trial from zero tractions and stops at its first
//! sweep within the 1e-12 m sticking tolerance, so its answer is a step function of the held reaction.
//! Two speeds of the grid 0.01 + 0.19 i / 400 m/s (i = 59 and 60), each one step from the same start,
//! stop at different sweeps. Bisection between them by the stop reaches a speed whose step neither
//! the first pass nor the held pass accepts (on 8fa3491d19 it is refused with the connection residual).
//! The pass that runs the sweep to rest must accept it, as its trial count past twice the budget
//! discloses, with the energy closed.
use super::*;
use super::cone_edge::{MU_S, StretchQuartic};

const BUDGET: usize = 64;
const NORMAL_N: f64 = 10.0;
const START_M_SQRT_KG: f64 = 1e-4;
const SHAPES: [[f64; 3]; 3] = [[20.0, 4.0, -2.0], [18.0, -5.0, 3.0], [17.0, 2.0, 6.0]];

fn bow_config() -> ModalFrictionConfig {
    ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0, maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 }
}
fn bow_port(shapes: [f64; 3]) -> ModalFrictionPort {
    ModalFrictionPort { left: ModalAttachment { component: 0, shapes: shapes.to_vec() }, right: None,
        law: FrictionLaw::ShiftedPower { static_mu: MU_S, kinetic_mu: 0.0, characteristic_speed: 0.011, exponent: 0.4 } }
}
fn displaced() -> MultiFrictionModalSystem {
    let gate = CancelGate::new();
    let mut model = ModalAcousticTimeModel::try_new(48_000, [800.0, 1100.0, 1700.0].iter().map(|&angular_frequency_rad_s|
        ModalAcousticMode { angular_frequency_rad_s, damping_ratio: 0.002, pressure_per_modal_velocity: C64::ZERO })
        .collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    let q = START_M_SQRT_KG;
    model.restore_states(&[q, 0.3 * q, -0.2 * q].map(|displacement_m_sqrt_kg|
        ModalAcousticState { displacement_m_sqrt_kg, velocity_m_sqrt_kg_per_s: 0.0 })).unwrap();
    FrictionModalSystem::new(CoupledModalSystem::new(vec![model], vec![], budgets(), &gate).unwrap(),
        bow_port(SHAPES[0]), bow_config(), &gate).unwrap()
        .with_nonlinear_potentials(vec![ModalPotentialBinding { component: 0, coordinate_start: 0,
            potential: Box::new(StretchQuartic(1e13)) }], NonlinearModalConfig { max_iterations: BUDGET,
            maximum_generalized_force: 1e5, force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10 }, &gate)
        .unwrap()
        .with_additional_ports(SHAPES[1..].iter().map(|&shapes| (bow_port(shapes), bow_config())).collect(),
            joint_config(), &gate).unwrap()
}
// One step from the displaced start: the frame, and the energy before and after it.
fn first_step(speed: f64) -> (MultiFrictionFrame, f64, f64) {
    let mut system = displaced();
    let before = system.total_energy_j().unwrap();
    let frame = system.step(&[0.0; 3], &[0; 3].map(|_| drive(speed, NORMAL_N)))
        .unwrap_or_else(|error| panic!("one step at {speed:e} m/s refused: {error}")).clone();
    (frame, before, system.total_energy_j().unwrap())
}

#[test]
fn reaction_trials_between_two_sweep_stops_are_accepted_with_the_sweep_at_rest() {
    let trials = |frame: &MultiFrictionFrame| frame.nonlinear.map_or(0, |n| n.nonlinear_iterations);
    let grid = |i: u32| 0.01 + 0.19 * f64::from(i) / 400.0;
    let (mut low, mut high) = (grid(59), grid(60));
    let stop = first_step(low).0.sweeps;
    assert_ne!(stop, first_step(high).0.sweeps, "the two speeds' steps stop at different sweeps");
    let (frame, before, after) = loop {
        let mid = 0.5 * (low + high);
        assert!(low < mid && mid < high, "no step between {low:e} and {high:e} m/s needed the sweep at rest");
        let probe = first_step(mid);
        if trials(&probe.0) > 2 * BUDGET {
            break probe;
        }
        if probe.0.sweeps == stop { low = mid } else { high = mid }
    };
    let work = frame.external_work_j + frame.drive_work_j;
    let loss = frame.network_dissipation_j + frame.friction_dissipation_j + frame.declared_slip_work_j;
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert!((after - before + loss - work).abs() <= 1e-12 + 1e-8 * work.abs().max(loss));
}
