//! Independent contact energy, actuator work and endpoint-force checks with DG storage.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::normal_load::ContactNormalLoad;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;
use fs_dcontact::Obstacle;

const DT: f64 = 1.0 / 48_000.0;
const OMEGA: [f64; 2] = [100.0, 150.0];
const DAMPING: f64 = 0.01;
const BETA: f64 = 1e12;
const STIFFNESS: f64 = 2000.0;
const CONTACT_LOSS: f64 = 0.2;
const BOWS: [[f64; 2]; 2] = [[1.0, 0.0], [0.6, 0.8]];
const PRESSURE: [f64; 2] = [1.0, -0.4];
const INITIAL: [f64; 2] = [2e-5, 1e-4];

struct ConfinedQuartic(f64);
impl Storage for ConfinedQuartic {
    fn hamiltonian(&self, q: &[f64]) -> f64 { Quartic(BETA).hamiltonian(q) }
    fn gradient(&self, q: &[f64], out: &mut [f64]) { Quartic(BETA).gradient(q, out); }
}
impl ModalPotential for ConfinedQuartic {
    fn coordinate_count(&self) -> usize { 2 }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        Quartic(BETA).validate_coordinates(q)?;
        if q[0].abs() > self.0 {
            return Err(ModalCouplingError::Invalid { what: "quartic fixture displacement domain" });
        }
        Ok(())
    }
}

fn fixture(initial: [f64; 2], maximum_first: f64, contact_loaded: bool) -> MultiFrictionModalSystem {
    let gate = CancelGate::new();
    let mut model = ModalAcousticTimeModel::try_new(48_000, OMEGA.into_iter().map(|omega|
        ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: DAMPING,
            pressure_per_modal_velocity: C64::ZERO }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&initial.map(|q| ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: 0.0,
    })).unwrap();
    let row = |shapes: [f64; 2]| ModalAttachment { component: 0, shapes: shapes.to_vec() };
    let port = |shapes| ModalFrictionPort { left: row(shapes), right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } };
    let core = FrictionModalSystem::new(
        CoupledModalSystem::new(vec![model], vec![], budgets(), &gate).unwrap(),
        port(BOWS[0]), config(), &gate).unwrap()
        .with_acceleration_pressure(PRESSURE.to_vec()).unwrap()
        .with_nonlinear_potentials(vec![ModalPotentialBinding { component: 0, coordinate_start: 0,
            potential: Box::new(ConfinedQuartic(maximum_first)) }], NonlinearModalConfig {
            max_iterations: 64, maximum_generalized_force: 1e5,
            force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10,
        }, &gate).unwrap();
    let normal = ModalContact { left: row([0.0, 1.0]), right: row([0.0, 0.0]),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], STIFFNESS, 1.0,
            "synthetic shared-quartic contact oracle".into()).unwrap()
            .with_internal_loss(CONTACT_LOSS).unwrap() };
    let system = core.with_additional_ports(vec![(port(BOWS[1]), config())], joint_config(), &gate).unwrap()
        .with_obstacles(vec![(normal, ModalContactConfig {
            max_iterations: 96, maximum_force_n: 50.0, maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10,
        })], FrictionObstacleConfig { max_contacts: 1, max_setup_terms: 10_000 }, &gate).unwrap();
    if contact_loaded {
        system.with_contact_normal_loads(vec![ContactNormalLoad { friction_port: 0, obstacle: 0 }], &gate).unwrap()
    } else { system }
}
fn moving(system: MultiFrictionModalSystem) -> MultiFrictionModalSystem {
    system.with_obstacle_motion(ObstacleMotionLimits {
        maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0,
    }, &CancelGate::new()).unwrap()
}
fn drives(contact_loaded: bool) -> [FrictionDrive; 2] {
    [drive(0.2, if contact_loaded { 0.0 } else { 0.3 }), drive(-0.12, 0.2)]
}
fn contact_energy(closure: f64) -> f64 { 0.5 * STIFFNESS * closure.max(0.0).powi(2) }

#[test]
fn quartic_fixed_and_moving_stops_close_physical_work_and_endpoint_pressure() {
    for translating in [false, true] {
        // The fixed stop is independent of either prescribed bow. The moving
        // stop instead supplies bow zero's normal load from the current solve.
        let mut system = fixture(INITIAL, 1.0, translating);
        if translating { system = moving(system); }
        let initial = system.total_energy_j().unwrap();
        let (mut supply, mut loss, mut old_offset) = (0.0, 0.0, 0.0);
        let (mut released, mut recontacted) = (false, false);
        let mut contact_loss = 0.0;
        for sample in 0..192 {
            let before = system.components()[0].states().to_vec();
            let offset = if !translating || sample < 32 { 0.0 }
                else if sample < 96 { 4e-6 * (sample - 31) as f64 }
                else if sample < 128 { 256e-6 }
                else { 4e-6 * (191 - sample) as f64 };
            let external = [0.03, 0.1];
            let inputs = drives(translating);
            let frame = if translating {
                system.step_with_obstacle_offsets(&external, &inputs, &[offset])
            } else { system.step(&external, &inputs) }.unwrap().clone();
            let normal = system.last_obstacle_frames().unwrap()[0];
            let after = system.components()[0].states();
            let q = [after[0].displacement_m_sqrt_kg, after[1].displacement_m_sqrt_kg];
            let delta = [q[0] - before[0].displacement_m_sqrt_kg, q[1] - before[1].displacement_m_sqrt_kg];
            let old_closure = before[1].displacement_m_sqrt_kg - old_offset;
            let closure = q[1] - offset;
            let change = closure - old_closure;
            let energy_change = contact_energy(closure) - contact_energy(old_closure);
            // Analytic integral of a unilateral linear spring, including a
            // crossing of its zero-force boundary; no contact-owner calls.
            let elastic = if closure > 0.0 && old_closure > 0.0 {
                0.5 * STIFFNESS * (closure + old_closure)
            } else if change == 0.0 { STIFFNESS * closure.max(0.0) }
            else { energy_change / change };
            let reaction = (elastic * (1.0 + CONTACT_LOSS * change / DT)).max(0.0);
            assert!((normal.normal_force_n - reaction).abs() < 1e-9);
            assert!((frame.obstacle_energy_j - contact_energy(closure)).abs() < 1e-15);
            assert!((frame.obstacle_dissipation_j - (normal.normal_force_n * change - energy_change)).abs() < 1e-12);
            let support_work = -normal.normal_force_n * (offset - old_offset);
            assert!((frame.obstacle_work_j - support_work).abs() < 1e-18);
            assert_eq!(normal.actuator_work_j, frame.obstacle_work_j);
            assert_eq!(normal.support_offset_m, offset);
            if translating {
                assert_eq!(frame.ports[0].normal_force_n, normal.normal_force_n);
                if normal.normal_force_n == 0.0 {
                    released = true;
                    assert_eq!(frame.ports[0].traction_n, 0.0);
                    assert_eq!(frame.ports[0].regime, ModalFrictionRegime::Released);
                } else if sample >= 128 { recontacted = true; }
            }
            let mut bow_work = 0.0;
            for (index, port) in frame.ports.iter().enumerate() {
                let surface_distance = inputs[index].speed_m_s * DT;
                let slip = surface_distance - BOWS[index][0] * delta[0] - BOWS[index][1] * delta[1];
                bow_work += port.traction_n * surface_distance;
                assert!((port.dissipation_j - port.traction_n * slip).abs() < 1e-12);
            }
            let external_work = external[0] * delta[0] + external[1] * delta[1];
            assert!((frame.external_work_j - external_work).abs() < 1e-17);
            assert!((frame.drive_work_j - bow_work).abs() < 1e-17);
            let quartic = 0.25 * BETA * (q[0]*q[0] + q[1]*q[1]).powi(2);
            let linear = (0..2).map(|i| 0.5 * (after[i].velocity_m_sqrt_kg_per_s.powi(2)
                + OMEGA[i].powi(2) * q[i].powi(2))).sum::<f64>();
            assert!((frame.network_energy_j - linear - quartic).abs() < 1e-14);
            let diagnostic = frame.nonlinear.unwrap();
            assert!((diagnostic.nonlinear_energy_j - quartic).abs() < 1e-15);
            assert!(diagnostic.reaction_residual <= diagnostic.reaction_tolerance);
            let pressure = (0..2).map(|i| PRESSURE[i] * (external[i]
                + BOWS[0][i] * frame.ports[0].traction_n + BOWS[1][i] * frame.ports[1].traction_n
                - if i == 1 { normal.normal_force_n } else { 0.0 }
                - 2.0 * DAMPING * OMEGA[i] * after[i].velocity_m_sqrt_kg_per_s
                - OMEGA[i].powi(2) * q[i] - BETA * (q[0]*q[0] + q[1]*q[1]) * q[i])).sum::<f64>();
            assert!((frame.observer_pressure_pa - pressure).abs() < 1e-9);
            supply += external_work + bow_work + support_work;
            loss += frame.network_dissipation_j + frame.friction_dissipation_j + frame.obstacle_dissipation_j;
            contact_loss += frame.obstacle_dissipation_j;
            assert!((system.total_energy_j().unwrap() - initial + loss - supply).abs() < 1e-9);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            old_offset = offset;
        }
        assert!(contact_loss > 0.0);
        if translating { assert!(released && recontacted, "normal load must release and return with the contact"); }
    }
}

#[test]
fn rejected_constitutive_trial_and_cancel_preserve_joint_state_then_retry_exactly() {
    let mut system = moving(fixture(INITIAL, 3e-5, true));
    let mut reference = moving(fixture(INITIAL, 3e-5, true));
    let external = [0.03, 0.1];
    let inputs = drives(true);
    for sample in 1..=8 {
        let offsets = [-1e-6 * f64::from(sample)];
        assert_eq!(system.step_with_obstacle_offsets(&external, &inputs, &offsets).unwrap(),
            reference.step_with_obstacle_offsets(&external, &inputs, &offsets).unwrap());
    }
    let states = all_states(&system);
    let accepted = system.last_frame().unwrap().clone();
    let normals = system.last_obstacle_frames().unwrap().to_vec();
    let offsets = system.obstacle_offsets_m().unwrap().to_vec();
    let energy = system.total_energy_j().unwrap();
    // This finite force passes the ordinary mechanical budgets but its staged
    // displacement leaves the material domain after normal/friction prediction.
    let error = system.step_with_obstacle_offsets(&[1e5, 0.1], &inputs, &[-1e-5]).unwrap_err();
    assert!(error.to_string().contains("quartic fixture displacement domain"), "{error}");
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(system.step_with_obstacle_offsets_under_gate(&external, &inputs, &[-1e-5], &cancelled).is_err());
    assert_eq!(all_states(&system), states);
    assert_eq!(system.last_frame(), Some(&accepted));
    assert_eq!(system.last_obstacle_frames().unwrap(), normals);
    assert_eq!(system.obstacle_offsets_m().unwrap(), offsets);
    assert_eq!(system.total_energy_j().unwrap(), energy);
    assert_eq!(system.samples_rendered(), 8);
    for sample in 9..=24 {
        let offsets = [-1e-6 * f64::from(sample)];
        assert_eq!(system.step_with_obstacle_offsets(&external, &inputs, &offsets).unwrap(),
            reference.step_with_obstacle_offsets(&external, &inputs, &offsets).unwrap());
        assert_eq!(all_states(&system), all_states(&reference));
        assert_eq!(system.last_obstacle_frames(), reference.last_obstacle_frames());
    }
}

#[test]
fn static_preload_solves_quartic_and_normal_storage_together_then_releases() {
    let mut system = fixture([0.0; 2], 1.0, true);
    let energy = system.initialize_released_equilibrium(&[0.0, 0.2], &CancelGate::new()).unwrap();
    // Independent scalar force equation: (omega^2+k)q + beta q^3 = F.
    let (mut low, mut high) = (0.0, 0.2 / (OMEGA[1].powi(2) + STIFFNESS));
    for _ in 0..96 {
        let q = f64::midpoint(low, high);
        if (OMEGA[1].powi(2) + STIFFNESS) * q + BETA * q.powi(3) > 0.2 { high = q; }
        else { low = q; }
    }
    let q = f64::midpoint(low, high);
    assert_eq!(all_states(&system)[0].displacement_m_sqrt_kg, 0.0);
    assert!((all_states(&system)[1].displacement_m_sqrt_kg - q).abs() < 1e-10*q);
    let expected = 0.5*(OMEGA[1].powi(2) + STIFFNESS)*q*q + 0.25*BETA*q.powi(4);
    assert!((energy-expected).abs() < 1e-10*expected);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
    assert!(system.last_obstacle_frames().is_none());
    assert_eq!(system.total_energy_j().unwrap(), energy);
    let released = [drive(0.0, 0.0); 2];
    for _ in 0..8 {
        let frame = system.step(&[0.0, 0.2], &released).unwrap();
        assert!(frame.observer_pressure_pa.abs() < 1e-8, "held contact equilibrium radiates");
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    let normal = system.last_obstacle_frames().unwrap()[0];
    assert!((normal.normal_force_n - STIFFNESS*q).abs() < 1e-9);
    let frame = system.step(&[0.0; 2], &released).unwrap();
    assert!(frame.observer_pressure_pa.abs() > 0.05);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}
