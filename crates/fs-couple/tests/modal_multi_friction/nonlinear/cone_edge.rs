//! G1: a port at the static cone's edge under a stiffening stretch reaction.
//!
//! One mode carries a manufactured quartic stretch energy and one bowed port
//! with a steeply weakening shifted-power law (the bowed card's law with a
//! smaller weakening speed, so that leaving the cone jumps to the high-slip
//! root; no low sliding branch exists). The normal load is set from the
//! converged sticking traction T* of the same step, just under T*/mu_s:
//! - sticking drags the mode further, the stretch reaction grows, and the
//!   held no-slip force T* then lies just outside the cone;
//! - the high-slip root leaves the mode nearly still, the reaction falls, and
//!   the no-slip force then lies inside the cone again.
//! Sticking is therefore inadmissible under its own reaction, while sliding on
//! the high branch is consistent under its own. A reaction iteration that
//! re-decides every trial from the step's accepted memory alternates between
//! the two and never settles; the step must instead slide, with energy closed.
use super::*;

const MU_S: f64 = 1.2147108804996467;
const BOW_M_S: f64 = 0.1;

struct StretchQuartic(f64);
impl Storage for StretchQuartic {
    fn hamiltonian(&self, q: &[f64]) -> f64 { 0.25 * self.0 * q[0].powi(4) }
    fn gradient(&self, q: &[f64], out: &mut [f64]) { out[0] = self.0 * q[0].powi(3); }
}
impl ModalPotential for StretchQuartic {
    fn coordinate_count(&self) -> usize { 1 }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() != 1 || !q[0].is_finite() || q[0].abs() > 1e-3 {
            return Err(ModalCouplingError::Invalid { what: "manufactured finite stretch quartic domain" });
        }
        Ok(())
    }
}

fn weakening() -> FrictionLaw {
    FrictionLaw::ShiftedPower { static_mu: MU_S, kinetic_mu: 0.0, characteristic_speed: 1e-4, exponent: 0.4 }
}
fn port(shapes: Vec<f64>) -> ModalFrictionPort {
    ModalFrictionPort { left: ModalAttachment { component: 0, shapes }, right: None, law: weakening() }
}
fn single() -> FrictionModalSystem { bowed(vec![20.0, 0.0]) }
fn bowed(shapes: Vec<f64>) -> FrictionModalSystem {
    let gate = CancelGate::new();
    let model = ModalAcousticTimeModel::try_new(48_000, [800.0, 1100.0].iter().map(|&angular_frequency_rad_s|
        ModalAcousticMode { angular_frequency_rad_s, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    FrictionModalSystem::new(CoupledModalSystem::new(vec![model], vec![], budgets(), &gate).unwrap(),
        port(shapes), config(), &gate).unwrap()
        .with_nonlinear_potentials(vec![ModalPotentialBinding { component: 0, coordinate_start: 0,
            potential: Box::new(StretchQuartic(1.7e20)) }], NonlinearModalConfig {
            max_iterations: 64, maximum_generalized_force: 1e5,
            force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10,
        }, &gate).unwrap()
}
// The second port touches only the other mode and stays lifted, so the bowed
// port reaches the same edge through the joint sweep (the bow-ribbon route).
fn joint() -> MultiFrictionModalSystem {
    single().with_additional_ports(vec![(port(vec![0.0, 10.0]), config())], joint_config(), &CancelGate::new()).unwrap()
}
// Two bowed ports with nearly the same shapes on both modes, as neighbouring
// stations of a bow ribbon have. Their joint sweep cannot converge in its 128
// sweeps, so every reaction trial goes to the active-set rescue.
fn ribbon() -> MultiFrictionModalSystem {
    bowed(vec![20.0, 1.0]).with_additional_ports(vec![(port(vec![19.0, -1.0]), config())], joint_config(),
        &CancelGate::new()).unwrap()
}
fn edge_load(sticking_traction_n: f64) -> f64 { sticking_traction_n / MU_S * (1.0 - 1e-6) }

#[test]
fn g1_cone_edge_port_slides_where_its_own_stretch_reaction_refuses_sticking() {
    let mut probe = single();
    let stuck = *probe.step(&[0.0, 0.0], drive(BOW_M_S, 75.0)).unwrap();
    assert_eq!(stuck.regime, ModalFrictionRegime::Sticking);
    let normal = edge_load(stuck.traction_n.abs());
    let mut system = single();
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0, 0.0], drive(BOW_M_S, normal))
        .unwrap_or_else(|error| panic!("cone-edge step refused at {normal} N: {error}"));
    assert_eq!(frame.regime, ModalFrictionRegime::Sliding);
    assert!(frame.traction_n.abs() < MU_S * normal);
    let after = system.total_energy_j().unwrap();
    let work = frame.external_work_j + frame.drive_work_j;
    let loss = frame.network_dissipation_j + frame.friction_dissipation_j + frame.declared_slip_work_j;
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert!((after - before + loss - work).abs() <= 1e-12 + 1e-8 * work.abs().max(loss));
}

#[test]
fn g1_cone_edge_port_slides_through_the_joint_sweep() {
    let mut probe = joint();
    let stuck = probe.step(&[0.0, 0.0], &[drive(BOW_M_S, 75.0), drive(0.0, 0.0)]).unwrap().clone();
    assert_eq!(stuck.ports[0].regime, ModalFrictionRegime::Sticking);
    let normal = edge_load(stuck.ports[0].traction_n.abs());
    let mut system = joint();
    let before = system.total_energy_j().unwrap();
    let frame = system.step(&[0.0, 0.0], &[drive(BOW_M_S, normal), drive(0.0, 0.0)])
        .unwrap_or_else(|error| panic!("cone-edge joint step refused at {normal} N: {error}")).clone();
    assert_eq!(frame.ports[0].regime, ModalFrictionRegime::Sliding);
    assert!(frame.ports[0].traction_n.abs() < MU_S * normal);
    let after = system.total_energy_j().unwrap();
    let work = frame.external_work_j + frame.drive_work_j;
    let loss = frame.network_dissipation_j + frame.friction_dissipation_j + frame.declared_slip_work_j;
    assert!((after - before + loss - work).abs() <= 1e-12 + 1e-8 * work.abs().max(loss));
}

// Each port is loaded at its own cone edge. The rescue keeps no branch memory
// and starts its walk at all-stuck on every reaction trial. A partition that
// holds under the reaction of another partition's response and not under its
// own makes the walk alternate with the trials, until the held pass keeps the
// partition it last resolved.
#[test]
fn g1_cone_edge_ribbon_slides_through_the_active_set_rescue() {
    let drives = |normal: [f64; 2]| [drive(BOW_M_S, normal[0]), drive(BOW_M_S, normal[1])];
    let mut probe = ribbon();
    let stuck = probe.step(&[0.0, 0.0], &drives([75.0, 75.0])).unwrap().clone();
    assert_eq!((stuck.sweeps, stuck.active_set_rescues), (128, 1), "the probe goes through the rescue");
    assert!(stuck.ports.iter().all(|port| port.regime == ModalFrictionRegime::Sticking));
    let normal = [edge_load(stuck.ports[0].traction_n.abs()), edge_load(stuck.ports[1].traction_n.abs())];
    let mut system = ribbon();
    let before = system.total_energy_j().unwrap();
    let frame = system.step(&[0.0, 0.0], &drives(normal))
        .unwrap_or_else(|error| panic!("cone-edge ribbon step refused at {normal:?} N: {error}")).clone();
    assert_eq!((frame.sweeps, frame.active_set_rescues), (128, 1));
    assert!(frame.ports.iter().any(|port| port.regime == ModalFrictionRegime::Sliding));
    for (port, normal) in frame.ports.iter().zip(normal) {
        assert!(port.traction_n.abs() <= MU_S * normal * (1.0 + 1e-9));
    }
    let after = system.total_energy_j().unwrap();
    let work = frame.external_work_j + frame.drive_work_j;
    let loss = frame.network_dissipation_j + frame.friction_dissipation_j + frame.declared_slip_work_j;
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert!((after - before + loss - work).abs() <= 1e-12 + 1e-8 * work.abs().max(loss));
}
