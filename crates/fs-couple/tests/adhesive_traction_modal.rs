//! Original spatial loads on real modal, constitutive-memory and air states.
#[path = "../../fs-contact/tests/normal_work_derivative/fixtures.rs"]
mod fixtures;

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_contact::normal_patch::{
    maugis::{MaugisEnvelope, SourceMaugisSphereContact},
    traction::{
        AdhesiveTractionProfile,
        sectors::{PolarTransferBudget, RadialMomentBudget},
    },
};
use fs_couple::{
    modal_acoustic_time::{
        ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    },
    render::{
        plate::impact::radiation::{Model, Pole},
        schedule::force::coupled::{
            CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
            ModalCouplingError, adhesive_traction::*,
        },
    },
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_tribo::InputAuthority;

fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn network(rich: bool) -> CoupledModalSystem {
    let models = [500.0, 900.0]
        .into_iter()
        .enumerate()
        .map(|(i, omega)| {
            let mut m = ModalAcousticTimeModel::try_new(
                48_000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: omega,
                    damping_ratio: if rich { 0.02 } else { 0.0 },
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            if rich {
                m.restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: [1e-6, -2e-6][i],
                    velocity_m_sqrt_kg_per_s: [0.01, -0.005][i],
                }])
                .unwrap();
            }
            m
        })
        .collect();
    let links = if rich {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![0.3],
            },
            stiffness_n_m: 4e4,
            damping_n_s_m: 20.0,
            rest_extension_m: 1e-6,
        }]
    } else {
        vec![]
    };
    let native = CoupledModalSystem::new(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-14,
            energy_relative_tolerance: 1e-9,
        },
        &gate(),
    )
    .unwrap();
    if !rich {
        return native;
    }
    native
        .with_material_relaxation(
            vec![RelaxationBranch {
                projection: vec![1.0, 0.0, -0.4, 0.0],
                stiffness: 3e4,
                relaxation_time_s: 0.002,
            }],
            Some(vec![0.00007]),
            &gate(),
        )
        .unwrap()
        .with_radiation_load(
            &Model {
                ports: 2,
                poles: vec![Pole {
                    omega: 700.0,
                    zeta: 0.1,
                    coupling: vec![300.0, -200.0],
                }],
            },
            &[0, 1],
            1,
            &gate(),
        )
        .unwrap()
}
fn config(profile: &AdhesiveTractionProfile) -> AdhesiveTractionModalConfig {
    let c = profile.interaction_radius_m();
    let ring = profile.annulus(0.0, c).unwrap();
    let f = ring.hertz_force_n - ring.adhesive_force_n;
    AdhesiveTractionModalConfig {
        maximum_trace_terms: 256,
        transfer: PolarTransferBudget {
            radial_moment: RadialMomentBudget {
                absolute_tolerance_n_m: 1e-12 * f * c,
                relative_tolerance: 1e-10,
                maximum_evaluations: 65_536,
            },
            force_absolute_tolerance_n: 1e-11 * f,
            moment_absolute_tolerance_n_m: 1e-11 * f * c,
        },
        work_absolute_tolerance_j: 1e-24,
        work_relative_tolerance: 2e-10,
        maximum_quadrature_work_estimate_j: 1e-16,
    }
}
fn trace(radii_m: Vec<f64>, a: f64, bending_only: bool) -> PolarModalTrace {
    let mut coefficients = Vec::new();
    for _ in 1..radii_m.len() {
        for sign in [1.0, -1.0, -1.0, 1.0] {
            coefficients.extend(if bending_only {
                // Continuous piecewise-affine |x|/a and its opposite-body trace.
                [
                    AffineNormalModalTrace {
                        value: 0.0,
                        gradient_per_m: [sign / a, 0.0],
                    },
                    AffineNormalModalTrace {
                        value: 0.0,
                        gradient_per_m: [-sign / a, 0.0],
                    },
                ]
            } else {
                [
                    AffineNormalModalTrace {
                        value: 1.0,
                        gradient_per_m: [0.3 * sign / a, 0.0],
                    },
                    AffineNormalModalTrace {
                        value: -0.7,
                        gradient_per_m: [0.0, 0.4 / a],
                    },
                ]
            });
        }
    }
    PolarModalTrace {
        radii_m,
        angles_rad: vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
        coefficients,
    }
}
fn states(n: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    n.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * a.abs().max(b.abs()) + 1e-22,
        "{a:e} != {b:e}"
    );
}
fn same_native(a: &CoupledModalSystem, b: &CoupledModalSystem) {
    for (x, y) in states(a).iter().zip(states(b)) {
        close(x.displacement_m_sqrt_kg, y.displacement_m_sqrt_kg);
        close(x.velocity_m_sqrt_kg_per_s, y.velocity_m_sqrt_kg_per_s);
    }
    for (x, y) in a.material_memory().iter().zip(b.material_memory()) {
        close(*x, *y);
    }
    for (x, y) in a.radiation_state().iter().zip(b.radiation_state()) {
        close(*x, *y);
    }
}
fn jkr_profile() -> AdhesiveTractionProfile {
    AdhesiveTractionProfile::jkr_attached(&fixtures::jkr(), 0.0).unwrap()
}
fn maugis() -> SourceMaugisSphereContact {
    let original = fixtures::jkr();
    let normal = original.normal_state().clone();
    let r = original.effective_radius_m();
    let a0 = (PI * normal.adhesion_energy_j_per_m2() * r * r
        / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0))
        .cbrt();
    let d0 = a0 * a0 / r;
    SourceMaugisSphereContact::new(
        normal,
        r,
        2.0 * d0 / (PI * 0.2),
        MaugisEnvelope {
            maximum_indentation_m: original.envelope().maximum_indentation_m,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}

#[test]
fn zero_net_load_still_drives_opposite_flexible_bodies_through_patch_moments() {
    let source = fixtures::jkr();
    let e = source.normal_state().elastic().reduced_modulus_pa();
    let w = source.normal_state().adhesion_energy_j_per_m2();
    let r = source.effective_radius_m();
    let a0 = (4.5 * PI * w * r * r / e).cbrt();
    let profile = AdhesiveTractionProfile::jkr_attached(&source, a0 * a0 / (3.0 * r)).unwrap();
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    assert!(whole.normal_force_n.abs() < 1e-11 * whole.hertz_force_n);
    // Independent closed-form integral of |x|p/a, with original H/JKR forces.
    let bending = 3.0 * whole.hertz_force_n / 8.0 + whole.adhesive_force_n / 2.0;
    assert!(bending.abs() > 0.1 * whole.hertz_force_n);
    let mut reference = network(false);
    let mut owner = AdhesiveTractionModalSystem::new(
        network(false),
        trace(vec![0.0, a], a, true),
        config(&profile),
        &gate(),
    )
    .unwrap();
    let frame = *owner.step(&profile, &[0.0; 2], &gate()).unwrap();
    reference.step(&[bending, -bending]).unwrap();
    same_native(owner.network(), &reference);
    assert!(states(owner.network())[0].velocity_m_sqrt_kg_per_s < 0.0);
    assert!(states(owner.network())[1].velocity_m_sqrt_kg_per_s > 0.0);
    assert!(frame.spatial_traction_work_j > 0.0);
    assert_eq!(frame.quadrature_work_estimate_j, 0.0);
    assert_eq!(frame.transfer.radial_evaluations, 0);
    assert_eq!(frame.transfer.model_identity, source.identity());
    assert_eq!(frame.transfer.authority, InputAuthority::SyntheticFixture);
    close(
        frame.traction_work_on_network_j,
        reference.last_frame().unwrap().external_work_j,
    );
}

#[test]
fn prescribed_spatial_load_retains_every_native_state_loss_and_continuing_clock() {
    let profile = jkr_profile();
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    let bending = 3.0 * whole.hertz_force_n / 8.0 + whole.adhesive_force_n / 2.0;
    let force = [
        whole.normal_force_n + 0.3 * bending,
        -0.7 * whole.normal_force_n,
    ];
    let mut reference = network(true);
    let mut continuing = network(true);
    reference.step(&[0.0; 2]).unwrap();
    continuing.step(&[0.0; 2]).unwrap();
    let mut owner = AdhesiveTractionModalSystem::new(
        continuing,
        trace(vec![0.0, 0.4 * a, a], a, false),
        config(&profile),
        &gate(),
    )
    .unwrap();
    for ordinal in 2..10 {
        let external = [0.001, -0.0003];
        let frame = *owner.step(&profile, &external, &gate()).unwrap();
        let expected = reference
            .step(&[external[0] + force[0], external[1] + force[1]])
            .unwrap()
            .clone();
        assert_eq!(frame.sample, ordinal);
        same_native(owner.network(), &reference);
        let native = owner.network().last_frame().unwrap();
        for (x, y) in [
            (
                native.material_dissipation_j,
                expected.material_dissipation_j,
            ),
            (
                native.radiation_dissipation_j,
                expected.radiation_dissipation_j,
            ),
            (
                native.connection_dissipation_j,
                expected.connection_dissipation_j,
            ),
            (
                native.component_dissipation_j,
                expected.component_dissipation_j,
            ),
        ] {
            assert!(x > 0.0);
            close(x, y);
        }
        close(
            frame.traction_work_on_network_j + frame.other_external_work_j,
            native.external_work_j,
        );
        assert!(frame.work_residual_j.abs() <= frame.work_tolerance_j);
        assert!(native.energy_residual_j.abs() <= native.energy_tolerance_j);
    }
}

#[test]
fn original_maugis_contact_noncontact_and_detached_fields_share_the_native_receiver() {
    let source = maugis();
    let contact =
        AdhesiveTractionProfile::maugis(&source, &source.evaluate_contact_energy(0.0).unwrap())
            .unwrap();
    let limits = source.noncontact_limits().unwrap().unwrap();
    let delta = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let noncontact =
        AdhesiveTractionProfile::maugis(&source, &source.evaluate_noncontact(delta).unwrap())
            .unwrap();
    let detached = AdhesiveTractionProfile::maugis(
        &source,
        &source
            .evaluate_detached(-2.0 * source.interaction_range_m())
            .unwrap(),
    )
    .unwrap();
    assert_eq!(noncontact.contact_radius_m(), 0.0);
    let a = contact.contact_radius_m();
    let c = 1.2
        * contact
            .interaction_radius_m()
            .max(noncontact.interaction_radius_m());
    let map = trace(vec![0.0, 0.4 * a, a, c], a, false);
    let controls = config(&contact);
    let mut owner =
        AdhesiveTractionModalSystem::new(network(true), map.clone(), controls, &gate()).unwrap();
    let first = *owner.step(&contact, &[0.0; 2], &gate()).unwrap();
    assert!(first.transfer.radial_evaluations > 0 && first.quadrature_work_estimate_j > 0.0);
    let second = *owner.step(&noncontact, &[0.0; 2], &gate()).unwrap();
    assert_eq!(second.transfer.radial_evaluations, 0);
    assert_eq!(second.quadrature_work_estimate_j, 0.0);
    let last = *owner.step(&detached, &[0.0; 2], &gate()).unwrap();
    assert_eq!(last.traction_work_on_network_j, 0.0);
    assert_eq!(last.spatial_traction_work_j, 0.0);
    let mut exact_work = controls;
    exact_work.maximum_quadrature_work_estimate_j = 0.0;
    let mut refused =
        AdhesiveTractionModalSystem::new(network(true), map, exact_work, &gate()).unwrap();
    assert!(matches!(
        refused.step(&contact, &[0.0; 2], &gate()),
        Err(AdhesiveTractionModalError::QuadratureWork { .. })
    ));
    assert_eq!(refused.network().samples_rendered(), 0);
}

#[test]
fn receiving_refusal_late_cancellation_and_invalid_support_publish_nothing() {
    let profile = jkr_profile();
    let a = profile.contact_radius_m();
    let mut owner = AdhesiveTractionModalSystem::new(
        network(true),
        trace(vec![0.0, a], a, false),
        config(&profile),
        &gate(),
    )
    .unwrap();
    owner.step(&profile, &[0.0; 2], &gate()).unwrap();
    let accepted = states(owner.network());
    let memory = owner.network().material_memory().to_vec();
    let radiation = owner.network().radiation_state().to_vec();
    let native_frame = owner.network().last_frame().unwrap().clone();
    let frame = *owner.last_frame().unwrap();
    let refused: Result<(), AdhesiveTractionModalError> =
        owner.step_with_admission(&profile, &[0.0; 2], &gate(), |_, _| {
            Err(AdhesiveTractionModalError::Invalid(
                "receiving owner refused",
            ))
        });
    assert!(refused.is_err());
    let cancelled = gate();
    let refused: Result<(), AdhesiveTractionModalError> =
        owner.step_with_admission(&profile, &[0.0; 2], &cancelled, |_, _| {
            cancelled.request();
            Ok(())
        });
    assert!(matches!(
        refused,
        Err(AdhesiveTractionModalError::Native(
            ModalCouplingError::Cancelled
        ))
    ));
    let large = AdhesiveTractionProfile::jkr_attached(&fixtures::jkr(), 8e-6).unwrap();
    assert!(owner.step(&large, &[0.0; 2], &gate()).is_err());
    assert_eq!(states(owner.network()), accepted);
    assert_eq!(owner.network().material_memory(), memory);
    assert_eq!(owner.network().radiation_state(), radiation);
    assert_eq!(owner.network().last_frame().unwrap(), &native_frame);
    assert_eq!(owner.last_frame().unwrap(), &frame);
    owner.step(&profile, &[0.0; 2], &gate()).unwrap();
    assert_eq!(owner.network().samples_rendered(), 2);
}

#[test]
fn contact_sized_work_gate_catches_a_small_load_erased_by_background_force() {
    let profile = jkr_profile();
    let a = profile.contact_radius_m();
    let mut controls = config(&profile);
    controls.work_absolute_tolerance_j = 0.0;
    controls.work_relative_tolerance = 1e-12;
    let map = PolarModalTrace {
        radii_m: vec![0.0, a],
        angles_rad: vec![0.0, TAU],
        coefficients: vec![
            AffineNormalModalTrace {
                value: 1e-20,
                gradient_per_m: [0.0; 2],
            },
            AffineNormalModalTrace::default(),
        ],
    };
    let mut owner =
        AdhesiveTractionModalSystem::new(network(true), map, controls, &gate()).unwrap();
    let before = states(owner.network());
    assert!(matches!(
        owner.step(&profile, &[1.0, 0.0], &gate()),
        Err(AdhesiveTractionModalError::Work { .. })
    ));
    assert_eq!(states(owner.network()), before);
    assert_eq!(owner.network().samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
}

#[test]
fn trace_shape_and_storage_limits_are_admitted_before_mechanical_work() {
    let profile = jkr_profile();
    let a = profile.contact_radius_m();
    let mut controls = config(&profile);
    controls.maximum_trace_terms = 7;
    assert!(
        AdhesiveTractionModalSystem::new(
            network(false),
            trace(vec![0.0, a], a, false),
            controls,
            &gate()
        )
        .is_err()
    );
    let mut map = trace(vec![0.0, a], a, false);
    map.coefficients[7].gradient_per_m[0] = f64::NAN;
    assert!(
        AdhesiveTractionModalSystem::new(network(false), map, config(&profile), &gate()).is_err()
    );
    let mut map = trace(vec![0.0, a], a, false);
    map.coefficients.pop();
    assert!(
        AdhesiveTractionModalSystem::new(network(false), map, config(&profile), &gate()).is_err()
    );
}
