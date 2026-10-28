//! Independent shared-quartic energy and endpoint-acceleration oracle.
use super::*;
use fs_couple::render::schedule::force::coupled::ModalCouplingError;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::{
    ModalPotential, ModalPotentialBinding, NonlinearModalConfig,
};
use fs_phs::Storage;

#[path = "nonlinear/obstacles.rs"]
mod obstacles;

#[path = "nonlinear/preload.rs"]
mod preload;

struct Quartic(f64);
impl Storage for Quartic {
    fn hamiltonian(&self, q: &[f64]) -> f64 { 0.25 * self.0 * (q[0]*q[0] + q[1]*q[1]).powi(2) }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        let stiffness = self.0 * (q[0]*q[0] + q[1]*q[1]);
        for (force, q) in out.iter_mut().zip(q) { *force = stiffness * q; }
    }
}
impl ModalPotential for Quartic {
    fn coordinate_count(&self) -> usize { 2 }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() != 2 || q.iter().any(|q| !q.is_finite()) {
            return Err(ModalCouplingError::Invalid { what: "finite two-coordinate quartic fixture required" });
        }
        Ok(())
    }
}

#[test]
fn simultaneous_bow_pressure_uses_endpoint_gradient_and_original_external_work() {
    let gate = CancelGate::new();
    let omega = [800.0, 1100.0];
    let mut model = ModalAcousticTimeModel::try_new(48_000, omega.iter().map(|&angular_frequency_rad_s|
        ModalAcousticMode { angular_frequency_rad_s, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[
        ModalAcousticState { displacement_m_sqrt_kg: 1e-4, velocity_m_sqrt_kg_per_s: 0.0 },
        ModalAcousticState { displacement_m_sqrt_kg: -1.5e-4, velocity_m_sqrt_kg_per_s: 0.0 },
    ]).unwrap();
    let port = |shapes| ModalFrictionPort { left: ModalAttachment { component: 0, shapes },
        right: None, law: law() };
    let core = FrictionModalSystem::new(
        CoupledModalSystem::new(vec![model], vec![], budgets(), &gate).unwrap(),
        port(vec![20.0, 0.0]), config(), &gate).unwrap()
        .with_acceleration_pressure(vec![1.0, 2.0]).unwrap()
        .with_nonlinear_potentials(vec![ModalPotentialBinding { component: 0, coordinate_start: 0,
            potential: Box::new(Quartic(1e12)) }], NonlinearModalConfig {
            max_iterations: 64, maximum_generalized_force: 1e5,
            force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10,
        }, &gate).unwrap();
    let mut system = core.with_additional_ports(vec![(port(vec![0.0, 10.0]), config())], joint_config(), &gate).unwrap();
    let initial = system.total_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for sample in 0..96 {
        let old: Vec<_> = system.components()[0].states().iter().map(|s| s.displacement_m_sqrt_kg).collect();
        let drives = if sample < 80 { [drive(0.3, 1.0), drive(-0.2, 0.7)] }
            else { [drive(0.0, 0.0); 2] };
        let external = [0.02, -0.03];
        let frame = system.step(&external, &drives).unwrap().clone();
        let q: Vec<_> = system.components()[0].states().iter().map(|s| s.displacement_m_sqrt_kg).collect();
        let stiffness = 1e12 * (q[0]*q[0] + q[1]*q[1]);
        let pressure = [1.0, 2.0].into_iter().enumerate().map(|(i, weight)|
            weight * (external[i] + [20.0, 10.0][i]*frame.ports[i].traction_n
                - omega[i]*omega[i]*q[i] - stiffness*q[i])).sum::<f64>();
        assert!((frame.observer_pressure_pa - pressure).abs() < 1e-10);
        let external_work = external.iter().zip(q.iter().zip(&old)).map(|(f, (q, old))| f*(q-old)).sum::<f64>();
        assert!((frame.external_work_j - external_work).abs() < 1e-18);
        assert!((frame.nonlinear.unwrap().nonlinear_energy_j - Quartic(1e12).hamiltonian(&q)).abs() < 1e-15);
        work += frame.external_work_j + frame.drive_work_j;
        loss += frame.network_dissipation_j + frame.friction_dissipation_j;
        assert!((system.total_energy_j().unwrap() - initial + loss - work).abs() < 1e-9);
    }
}
