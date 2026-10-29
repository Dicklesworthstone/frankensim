//! Source-domain admission uses the same staged physical state as all forces.
use super::*;
use crate::tensor_material::PlaneStressPoint;
use fs_material::visco_tensor::{TensorProjectionBudget, TensorProny, TensorPronyBranch};

#[test]
fn g1_zero_coupling_history_is_owned_once_across_disconnected_islands() {
    let gate = CancelGate::new_clock_free();
    let mut system = CoupledModalSystem::new(
        vec![
            model(48_000, 100.0, 0.0, 0.0, 0.0),
            model(48_000, 200.0, 0.0, 0.0, 0.0),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![
            arm(vec![0.0; 4], 1.0, 0.001),
            arm(vec![0.0, 0.0, 1.0, 0.0], 1.0, 0.002),
        ],
        Some(vec![0.02, 0.0]),
        &gate,
    )
    .unwrap();
    let mut heat = 0.0;
    for tick in 1..=8 {
        let frame = system.step(&[0.0, 0.0]).unwrap().clone();
        let expected = 0.02 * (-f64::from(tick) / 48.0).exp();
        close(system.material_memory()[0], expected);
        assert_eq!(system.material_memory()[1], 0.0);
        assert!(
            states(&system)
                .iter()
                .all(|s| *s == ModalAcousticState::default())
        );
        heat += frame.material_dissipation_j;
        close(frame.material_energy_j + heat, 0.5 * 0.02_f64.powi(2));
    }
}

fn tensor_network() -> CoupledModalSystem {
    tensor_network_coordinates(false)
}
fn tensor_network_coordinates(shared:bool) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut factor = [[0.0; 6]; 6];
    factor[0] = [10.0, 0.0, 4.0, 1.0, 0.0, 0.0];
    let law = TensorProny::from_tensor_shear(
        fs_material::IsotropicElastic::new(1000.0, 0.25, 0.01)
            .unwrap()
            .stiffness(),
        vec![TensorPronyBranch::from_gram_rows(factor, 0.001).unwrap()],
        1e-6,
    )
    .unwrap();
    let image = PlaneStressPoint::new(
        &law,
        &[[1.0, 0.0, 0.0]],
        1.0,
        TensorProjectionBudget {
            max_coordinates: 4,
            max_scalar_terms: 1_000_000,
        },
        &gate,
    )
    .unwrap()
    .relaxation_image(1_000_000, &gate)
    .unwrap();
    let times:Vec<_>=image.relaxation_rates_per_s().iter().map(|rate|rate.recip()).collect();
    let shared=shared.then(||fs_material::visco_tensor::coordinates::RelaxationCoordinates::new(
        image.coupling_rows(),&times,1,
        fs_material::visco_tensor::coordinates::RelaxationCoordinateBudget {
            maximum_states:6,maximum_scalar_terms:1_000_000,
        },||false).unwrap());
    let (rows,times)=match &shared {
        Some(shared)=>(shared.coupling_rows(),shared.relaxation_times_s()),
        None=>(image.coupling_rows(),times.as_slice()),
    };
    let count=rows.len();
    let branches = rows
        .iter()
        .zip(times)
        .map(|(row, time)| arm(vec![row[0], 0.0], 1.0, *time))
        .collect();
    let omega = image.equilibrium_stiffness()[0].sqrt();
    mechanical(48_000, omega, 0.0)
        .with_material_relaxation(branches, None, &gate)
        .unwrap()
        .with_plane_stress_material_points(
            vec![ModalPlaneStressPoint {
                modes: 0..1,
                memories: 0..count,
                memory_map: shared.as_ref().map(|s|s.memory_map(0..6).unwrap()),
                strain_scale: 1.0,
                image,
            }],
            &gate,
        )
        .unwrap()
}

#[test]
fn g0_tensor_strain_rejects_static_dynamic_and_cancelled_candidates_atomically() {
    for shared in [false,true] {
    let gate = CancelGate::new_clock_free();
    let mut system = tensor_network_coordinates(shared);
    let before = states(&system);
    let history = system.material_memory().to_vec();
    let error = system
        .initialize_static_equilibrium(&[0.01], &gate)
        .unwrap_err();
    assert!(error.to_string().contains("surface strain"));
    assert_eq!(states(&system), before);
    assert_eq!(system.material_memory(), history);
    assert!(system.last_frame().is_none());
    let error = system.step(&[50_000.0]).unwrap_err();
    assert!(error.to_string().contains("surface strain"));
    assert_eq!(states(&system), before);
    assert_eq!(system.material_memory(), history);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(system.step_under_gate(&[1e-4], &cancelled).is_err());
    let mut clean = tensor_network_coordinates(shared);
    assert_eq!(system.step(&[1e-4]).unwrap(), clean.step(&[1e-4]).unwrap());
    assert_eq!(states(&system), states(&clean));
    assert_eq!(system.material_memory(), clean.material_memory());
    }
}

#[test]
fn g0_tensor_source_domain_admits_the_contact_endpoint_not_its_free_trial() {
    use super::super::super::contact::{ContactModalSystem, ModalContactConfig};
    let gate = CancelGate::new_clock_free();
    let mut network = tensor_network();
    // A numerical proposal can leave the source domain before contact has
    // acted. It changes scratch only and must not abort the surrounding solve.
    network.stage_inner(&[50_000.0], Some(&gate)).unwrap();
    assert!(
        network
            .validate_staged_material_points(Some(&gate))
            .is_err()
    );
    assert_eq!(states(&network), [ModalAcousticState::default()]);
    assert_eq!(network.material_memory(), [0.0; 6]);
    let contact = |network, stiffness| {
        ContactModalSystem::new_with_column(
            network,
            vec![1.0],
            fs_dcontact::Obstacle::new(
                vec![-1.0],
                1,
                1,
                vec![0.0],
                vec![1.0],
                stiffness,
                1.0,
                "manufactured linear stop for source-domain admission".into(),
            )
            .unwrap(),
            ModalContactConfig {
                max_iterations: 128,
                maximum_force_n: 1e6,
                maximum_penetration_m: 1e-3,
                force_absolute_tolerance_n: 1e-9,
                force_relative_tolerance: 1e-12,
            },
            &gate,
        )
        .unwrap()
    };
    let mut stopped = contact(network, 1e12);
    let frame = *stopped.step(&[50_000.0]).unwrap();
    assert!(frame.normal_force_n > 40_000.0);
    assert!(stopped.components()[0].states()[0].displacement_m_sqrt_kg < 1e-6);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert_eq!(stopped.samples_rendered(), 1);
    // A softer stop leaves the actual material beyond its source domain.
    // That endpoint must refuse without publishing contact, mechanics or clock.
    let mut insufficient = contact(tensor_network(), 1e4);
    assert!(
        insufficient
            .step(&[50_000.0])
            .unwrap_err()
            .to_string()
            .contains("surface strain")
    );
    assert_eq!(
        insufficient.components()[0].states(),
        [ModalAcousticState::default()]
    );
    assert!(insufficient.last_frame().is_none());
    assert_eq!(insufficient.samples_rendered(), 0);
}

#[test]
fn g3_tensor_domain_survives_extra_coordinates_nonlinear_preload_and_cold_restart() {
    use super::super::super::contact::friction::nonlinear::{
        MappedModalPotentialBinding, ModalCoordinate, ModalPotential, NonlinearModalConfig,
    };
    struct Quartic;
    impl fs_phs::Storage for Quartic {
        fn hamiltonian(&self, q: &[f64]) -> f64 {
            0.25e8 * q[0].powi(4)
        }
        fn gradient(&self, q: &[f64], g: &mut [f64]) {
            g[0] = 1e8 * q[0].powi(3);
        }
    }
    impl ModalPotential for Quartic {
        fn coordinate_count(&self) -> usize {
            1
        }
        fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
            if q.len() != 1 || !q[0].is_finite() {
                return Err(invalid("quartic coordinate"));
            }
            Ok(())
        }
    }
    let gate = CancelGate::new_clock_free();
    let mut owner = friction(tensor_network())
        .with_internal_components(vec![model(48_000, 200.0, 0.0, 0.0, 0.0)], &gate)
        .unwrap()
        .with_mapped_nonlinear_potentials(
            vec![MappedModalPotentialBinding {
                coordinates: vec![ModalCoordinate {
                    component: 0,
                    mode: 0,
                }],
                replaces_linear: vec![],
                potential: Box::new(Quartic),
            }],
            NonlinearModalConfig {
                max_iterations: 128,
                maximum_generalized_force: 1e6,
                force_absolute_tolerance: 1e-11,
                force_relative_tolerance: 1e-9,
            },
            &gate,
        )
        .unwrap()
        .with_additional_ports(
            vec![],
            super::super::super::contact::friction::multiple::MultiFrictionConfig {
                max_ports: 1,
                max_sweeps: 128,
                max_setup_terms: 10_000,
            },
            &gate,
        )
        .unwrap();
    let before = owner.components()[0].states().to_vec();
    assert!(
        owner
            .initialize_released_equilibrium(&[0.01, 0.0], &gate)
            .unwrap_err()
            .to_string()
            .contains("surface strain")
    );
    assert_eq!(owner.components()[0].states(), before);
    assert_eq!(owner.material_memory(), [0.0; 6]);
    assert!(
        owner
            .restore_component_states(
                0,
                &[ModalAcousticState {
                    displacement_m_sqrt_kg: 1e-3,
                    velocity_m_sqrt_kg_per_s: 0.0,
                }]
            )
            .unwrap_err()
            .to_string()
            .contains("surface strain")
    );
    assert_eq!(owner.components()[0].states(), before);
    owner
        .initialize_released_equilibrium(&[1e-4, 0.0], &gate)
        .unwrap();
    let settled = owner.components()[0].states()[0];
    let released = FrictionDrive {
        speed_m_s: 0.0,
        normal_force_n: 0.0,
    };
    owner.step(&[1e-4, 0.0], &[released]).unwrap();
    close(
        owner.components()[0].states()[0].displacement_m_sqrt_kg,
        settled.displacement_m_sqrt_kg,
    );
    let memory = owner.material_memory().to_vec();
    let accepted = owner.components()[0].states().to_vec();
    assert!(
        owner
            .step(&[50_000.0, 0.0], &[released])
            .unwrap_err()
            .to_string()
            .contains("surface strain")
    );
    assert_eq!(owner.material_memory(), memory);
    assert_eq!(owner.components()[0].states(), accepted);
    assert_eq!(owner.samples_rendered(), 1);
}
