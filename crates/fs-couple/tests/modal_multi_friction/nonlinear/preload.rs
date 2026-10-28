//! G1/G3: independently solved force-loaded nonlinear storage and release.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionError;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::{
    MappedModalPotentialBinding, ModalCoordinate,
};

struct LocalQuartic {
    coefficient: f64,
    maximum: f64,
}
impl Storage for LocalQuartic {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        0.25 * self.coefficient * q[0].powi(4)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        out[0] = self.coefficient * q[0].powi(3);
    }
}
impl ModalPotential for LocalQuartic {
    fn coordinate_count(&self) -> usize {
        1
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() != 1 || !q[0].is_finite() || q[0].abs() > self.maximum {
            Err(ModalCouplingError::Invalid {
                what: "manufactured finite-strain quartic domain",
            })
        } else {
            Ok(())
        }
    }
}
fn settings(iterations: usize) -> NonlinearModalConfig {
    NonlinearModalConfig {
        max_iterations: iterations,
        maximum_generalized_force: 1e5,
        force_absolute_tolerance: 1e-11,
        force_relative_tolerance: 1e-10,
    }
}
fn base(memory: bool, gauge: bool) -> FrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let point = |component| ModalAttachment {
        component,
        shapes: vec![1.0],
    };
    let mut models = vec![image(100.0, 0.0, 0.0), image(150.0, 0.0, 0.0)];
    if gauge {
        models.push(
            ModalAcousticTimeModel::try_new_with_free_coordinates(
                48000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        );
    }
    let mut network = CoupledModalSystem::new(
        models,
        vec![ModalConnection {
            left: point(0),
            right: point(1),
            stiffness_n_m: 4000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
        budgets(),
        &gate,
    )
    .unwrap();
    if memory {
        let mut projection = vec![0.0; 2 * network.mode_count()];
        projection[0] = 1.0;
        projection[2] = -0.5;
        network = network
            .with_material_relaxation(
                vec![fs_phs::RelaxationBranch {
                    projection,
                    stiffness: 5000.0,
                    relaxation_time_s: 0.01,
                }],
                None,
                &gate,
            )
            .unwrap();
    }
    let count = network.mode_count();
    let mut observer = vec![0.0; count];
    observer[0] = 1.0;
    FrictionModalSystem::new(
        network,
        ModalFrictionPort {
            left: point(0),
            right: None,
            law: law(),
        },
        config(),
        &gate,
    )
    .unwrap()
    .with_acceleration_pressure(observer)
    .unwrap()
}
fn bind(
    core: FrictionModalSystem,
    coefficient: f64,
    maximum: f64,
    iterations: usize,
) -> FrictionModalSystem {
    core.with_nonlinear_potentials(
        vec![ModalPotentialBinding {
            component: 0,
            coordinate_start: 0,
            potential: Box::new(LocalQuartic {
                coefficient,
                maximum,
            }),
        }],
        settings(iterations),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn cubic_root(force: f64, coefficient: f64) -> f64 {
    let linear = 14000.0 - 4000.0 * 4000.0 / 26500.0;
    let (mut low, mut high) = (0.0, force / linear);
    for _ in 0..96 {
        let q = f64::midpoint(low, high);
        if linear * q + coefficient * q.powi(3) > force {
            high = q;
        } else {
            low = q;
        }
    }
    f64::midpoint(low, high)
}

#[test]
fn strong_nonlinear_preload_matches_cubic_receiver_equilibrium_then_releases() {
    let gate = CancelGate::new_clock_free();
    for coefficient in [1e10, 1e12, 1e14] {
        let mut system = bind(base(false, false), coefficient, 0.003, 64);
        let stored = system
            .initialize_released_equilibrium(&[100.0, 0.0], &gate)
            .unwrap();
        let q = cubic_root(100.0, coefficient);
        let receiver = 4000.0 * q / 26500.0;
        assert!((states(&system)[0].displacement_m_sqrt_kg - q).abs() < 1e-10 * q);
        assert!((states(&system)[1].displacement_m_sqrt_kg - receiver).abs() < 1e-10 * q);
        let expected = 0.5 * 10000.0 * q * q
            + 0.5 * 22500.0 * receiver * receiver
            + 0.5 * 4000.0 * (q - receiver).powi(2)
            + 0.25 * coefficient * q.powi(4);
        assert!((stored - expected).abs() < 1e-9 * expected);
        assert_eq!(stored, system.total_energy_j().unwrap());
        assert_eq!(system.samples_rendered(), 0);
        assert!(system.last_frame().is_none());
        let report = *system.nonlinear_preload_diagnostics().unwrap();
        assert!(report.iterations <= 64 && report.force_residual <= report.force_tolerance);
        assert!(report.nonlinear_energy_j > 0.0);
        if coefficient == 1e14 {
            assert!(
                coefficient * q.powi(3) > 90.0,
                "nonlinear load must dominate this test"
            );
        }
        let held = *system.step(&[100.0, 0.0], drive(0.0, 0.0)).unwrap();
        assert!(
            held.observer_pressure_pa.abs() < 1e-7,
            "a held equilibrium must not launch a false pluck"
        );
        let release = *system.step(&[0.0, 0.0], drive(0.0, 0.0)).unwrap();
        assert!(release.observer_pressure_pa.abs() > 90.0);
        assert!(release.energy_residual_j.abs() <= release.energy_tolerance_j);
        assert_eq!(system.nonlinear_preload_diagnostics(), Some(&report));
    }
}

#[test]
fn zero_added_energy_keeps_linear_limit_and_relaxed_material_history() {
    let gate = CancelGate::new_clock_free();
    let mut linear = base(true, true);
    let mut nonlinear = bind(base(true, true), 0.0, 1.0, 64);
    let a = linear
        .initialize_released_equilibrium(&[1.0, 0.0, 0.0], &gate)
        .unwrap();
    let b = nonlinear
        .initialize_released_equilibrium(&[1.0, 0.0, 0.0], &gate)
        .unwrap();
    assert!((a - b).abs() < 1e-13 * a);
    for (a, b) in states(&linear).iter().zip(states(&nonlinear)) {
        assert!((a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg).abs() < 1e-15);
        assert_eq!(b.velocity_m_sqrt_kg_per_s, 0.0);
    }
    let states = states(&nonlinear);
    let expected = 5000.0_f64.sqrt()
        * (states[0].displacement_m_sqrt_kg - 0.5 * states[1].displacement_m_sqrt_kg);
    assert!((nonlinear.material_memory()[0] - expected).abs() < 1e-15);
    assert_eq!(nonlinear.material_energy_j().unwrap(), 0.0);
    assert_eq!(
        states[2].displacement_m_sqrt_kg, 0.0,
        "an unloaded fluid position gauge stays a gauge"
    );
    for _ in 0..32 {
        let a = *linear.step(&[0.0; 3], drive(0.0, 0.0)).unwrap();
        let b = *nonlinear.step(&[0.0; 3], drive(0.0, 0.0)).unwrap();
        assert!((a.observer_pressure_pa - b.observer_pressure_pa).abs() < 1e-10);
    }
}

struct CoupledStorage;
impl Storage for CoupledStorage {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        0.5 * 10000.0 * (q[0] - q[1]).powi(2)
            + 0.5 * 22500.0 * q[1] * q[1]
            + 0.25 * 1e12 * q[0].powi(4)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        out[0] = 10000.0 * (q[0] - q[1]) + 1e12 * q[0].powi(3);
        out[1] = 10000.0 * (q[1] - q[0]) + 22500.0 * q[1];
    }
}
impl ModalPotential for CoupledStorage {
    fn coordinate_count(&self) -> usize {
        2
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() == 2 && q.iter().all(|v| v.is_finite()) {
            Ok(())
        } else {
            Err(ModalCouplingError::Invalid {
                what: "finite mapped strain fixture",
            })
        }
    }
}
#[test]
fn mapped_full_strain_replaces_linear_stiffness_once_in_multi_contact_owner() {
    let gate = CancelGate::new_clock_free();
    let core = base(false, false)
        .with_mapped_nonlinear_potentials(
            vec![MappedModalPotentialBinding {
                coordinates: vec![
                    ModalCoordinate {
                        component: 0,
                        mode: 0,
                    },
                    ModalCoordinate {
                        component: 1,
                        mode: 0,
                    },
                ],
                replaces_linear: vec![0, 1],
                potential: Box::new(CoupledStorage),
            }],
            settings(128),
            &gate,
        )
        .unwrap();
    let mut system = core
        .with_additional_ports(vec![], joint_config(), &gate)
        .unwrap();
    let q0 = 0.0002_f64;
    let q1 = 14000.0 * q0 / 36500.0;
    let force = 14000.0 * (q0 - q1) + 1e12 * q0.powi(3);
    let energy = system
        .initialize_released_equilibrium(&[force, 0.0], &gate)
        .unwrap();
    let expected = CoupledStorage.hamiltonian(&[q0, q1]) + 0.5 * 4000.0 * (q0 - q1).powi(2);
    assert!((energy - expected).abs() < 1e-9 * expected);
    assert!((all_states(&system)[0].displacement_m_sqrt_kg - q0).abs() < 1e-10 * q0);
    assert!((all_states(&system)[1].displacement_m_sqrt_kg - q1).abs() < 1e-10 * q0);
    let report = system.nonlinear_preload_diagnostics().unwrap();
    assert!(report.force_residual <= report.force_tolerance);
}

#[test]
fn failed_domain_work_budget_or_cancellation_keeps_all_accepted_state() {
    let gate = CancelGate::new_clock_free();
    for (maximum, iterations) in [(0.01, 1), (1e-9, 64)] {
        let mut system = bind(base(true, false), 1e14, maximum, iterations);
        let before = states(&system);
        let memory = system.material_memory().to_vec();
        assert!(matches!(
            system.initialize_released_equilibrium(&[100.0, 0.0], &gate),
            Err(ModalFrictionError::NonlinearEquilibrium { .. })
        ));
        assert_eq!(states(&system), before);
        assert_eq!(system.material_memory(), memory);
        assert_eq!(system.samples_rendered(), 0);
        assert!(system.last_frame().is_none());
        assert!(system.nonlinear_preload_diagnostics().is_none());
    }
    let mut system = bind(base(false, false), 1e14, 0.01, 64);
    let before = states(&system);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.initialize_released_equilibrium(&[1.0, 0.0], &cancelled),
        Err(ModalFrictionError::Network(ModalCouplingError::Cancelled))
    ));
    assert_eq!(states(&system), before);
    system
        .initialize_released_equilibrium(&[1.0, 0.0], &gate)
        .unwrap();
    let loaded = states(&system);
    let report = *system.nonlinear_preload_diagnostics().unwrap();
    assert!(
        system
            .initialize_released_equilibrium(&[0.0, 0.0], &gate)
            .is_err()
    );
    assert_eq!(states(&system), loaded);
    assert_eq!(system.nonlinear_preload_diagnostics(), Some(&report));
}
