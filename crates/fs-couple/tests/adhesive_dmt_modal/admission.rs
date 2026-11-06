//! Small original contact work and atomic receiving owners on real mechanics.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_material::state_point::MaterialPropertySelection;

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod thermal;

fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}

// Unit mass-normalized free coordinates isolate a common external acceleration
// from the same original source's relative normal and tangential coordinates.
fn unit_columns() -> [Vec<f64>; 3] {
    core::array::from_fn(|i| {
        let mut b = vec![0.0; 6];
        b[i] = -1.0;
        b[3 + i] = 1.0;
        b
    })
}
fn unit_network(common_speed: f64, relative: [f64; 3]) -> CoupledModalSystem {
    let models = [-0.5, 0.5]
        .into_iter()
        .map(|sign| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000_000,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    };
                    3
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&relative.map(|v| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: common_speed + sign * v,
                }))
                .unwrap();
            model
        })
        .collect();
    CoupledModalSystem::new(
        models,
        vec![],
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 0,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e12,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-26,
            energy_relative_tolerance: 1e-9,
        },
        &gate(),
    )
    .unwrap()
}
fn unit_owner<'a>(
    normal: &'a fs_contact::normal_patch::dmt::SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    common_speed: f64,
    relative: [f64; 3],
    approach: f64,
) -> AdhesiveDmtModalSystem<'a> {
    let materials = normal
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    AdhesiveDmtModalSystem::new(
        unit_network(common_speed, relative),
        normal,
        shear,
        [&materials[0], &materials[1]],
        [approach, 0.0, 0.0],
        unit_columns(),
        config(),
        &gate(),
    )
    .unwrap()
}

#[test]
fn mechanical_receiving_refusal_and_late_cancellation_retain_original_native_history() {
    let (normal, shear) = sources(0.01, 1e7, 6e6);
    let initial = [0.2 * Z, 0.0, 0.0];
    let mut owner = make(&normal, &shear, [0.0, 1e-3, -0.5e-3], initial, true, RATE);
    let mut replay = make(&normal, &shear, [0.0, 1e-3, -0.5e-3], initial, true, RATE);
    owner.step(&[0.0; 6], [0.0; 3]).unwrap();
    replay.step(&[0.0; 6], [0.0; 3]).unwrap();
    let native = states(&owner);
    let source_state = owner.state();
    let last = *owner.last_frame().unwrap();
    let energy = owner.total_energy_j().unwrap().to_bits();
    let refused: Result<(), AdhesiveDmtModalError> =
        owner.step_with_admission(&[0.0; 6], [0.0; 3], &gate(), |frame, candidate| {
            assert_eq!(frame.sample, 2);
            assert_eq!(candidate.len(), 2);
            assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
            Err(AdhesiveDmtModalError::Network(
                ModalCouplingError::Cancelled,
            ))
        });
    assert!(refused.is_err());
    let cancelled = gate();
    let refused: Result<(), AdhesiveDmtModalError> =
        owner.step_with_admission(&[0.0; 6], [0.0; 3], &cancelled, |_, _| {
            cancelled.request();
            Ok(())
        });
    assert!(refused.is_err());
    assert_eq!(states(&owner), native);
    assert_eq!(owner.state(), source_state);
    assert_eq!(owner.last_frame(), Some(&last));
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.samples_rendered(), 1);
    let accepted: Result<_, AdhesiveDmtModalError> =
        owner.step_with_admission(&[0.0; 6], [0.0; 3], &gate(), |frame, _| Ok(*frame));
    assert_eq!(
        accepted.unwrap(),
        *replay.step(&[0.0; 6], [0.0; 3]).unwrap()
    );
    assert_eq!(states(&owner), states(&replay));
    assert_eq!(owner.state(), replay.state());
}

#[test]
fn original_contact_work_closes_independently_in_contact_gap_and_separation() {
    let (normal, shear) = sources(0.01, 1e7, 6e6);
    for approach in [0.2 * Z, -0.4 * Z, -2.0 * Z] {
        let mut owner = make(
            &normal,
            &shear,
            [0.0, 1e-3, -0.5e-3],
            [approach, 0.0, 0.0],
            true,
            RATE,
        );
        let before = normal.evaluate(approach).unwrap().reversible_energy_j;
        let prescribed = [0.0, 1e-4, -1e-4];
        let dt = owner.sample_period_s();
        let frame = *owner.step(&[0.0; 6], prescribed).unwrap();
        let actuator: f64 = (0..3)
            .map(|i| frame.mean_force_n[i] * prescribed[i] * dt)
            .sum();
        close(frame.prescribed_work_j, actuator, 1e-30, 1e-13);
        let residual = frame.surface_energy_j - before
            + frame.friction_dissipation_j
            + frame.contact_work_on_network_j
            - actuator;
        assert!(residual.abs() <= frame.contact_work_tolerance_j);
        assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
        assert!(frame.contact_work_tolerance_j < frame.energy_tolerance_j);
    }
}

#[test]
fn background_force_cannot_erase_adhesion_and_still_receive_source_work_credit() {
    let (normal, shear) = sources(0.01, 1e7, 6e6);
    let mut owner = unit_owner(&normal, &shear, 0.0, [1e-3, 0.0, 0.0], -0.5 * Z);
    let native = states(&owner);
    let state = owner.state();
    // Both native bodies share the large drive, so only the original small
    // relative velocity performs contact work. Each adhesive force correction
    // is smaller than half an ulp of its corresponding background input.
    let external = [1e9, 0.0, 0.0, 1e9, 0.0, 0.0];
    let mut called = false;
    let refused: Result<(), AdhesiveDmtModalError> =
        owner.step_with_admission(&external, [0.0; 3], &gate(), |_, _| {
            called = true;
            Ok(())
        });
    assert!(matches!(
        refused,
        Err(AdhesiveDmtModalError::ContactWork { .. })
    ));
    assert!(!called);
    assert_eq!(states(&owner), native);
    assert_eq!(owner.state(), state);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    owner.step(&[0.0; 6], [0.0; 3]).unwrap();
}

#[test]
fn thermal_work_is_resolved_below_native_energy_roundoff_and_publication_is_atomic() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let normal = fixture.dmt(Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let approach = 0.2 * Z;
    let refresh = normal
        .refresh_at_held_indentation(
            approach,
            fixture.cards(),
            AdhesiveTemperatures {
                surface_a_k: 310.0,
                surface_b_k: 290.0,
                interface_k: 315.0,
            },
        )
        .unwrap();
    let next_shear = shear.rebind_dmt_temperature(&normal, &refresh).unwrap();
    let mut still = unit_owner(&normal, &shear, 0.0, [0.0; 3], approach);
    let mut moving = unit_owner(&normal, &shear, 1e4, [0.0; 3], approach);
    let before = states(&moving);
    let original_state = moving.state();
    let energy = moving.total_energy_j().unwrap().to_bits();
    let refused: Result<(), AdhesiveDmtModalError> =
        moving.apply_source_refresh_with_admission(&refresh, &next_shear, &gate(), |_, state| {
            assert_eq!(state.position_m(), original_state.position_m());
            assert_eq!(state.regime(), original_state.regime());
            Err(AdhesiveDmtModalError::Network(
                ModalCouplingError::Cancelled,
            ))
        });
    assert!(refused.is_err());
    let cancelled = gate();
    let refused: Result<(), AdhesiveDmtModalError> =
        moving.apply_source_refresh_with_admission(&refresh, &next_shear, &cancelled, |_, _| {
            cancelled.request();
            Ok(())
        });
    assert!(refused.is_err());
    assert_eq!(states(&moving), before);
    assert_eq!(moving.state(), original_state);
    assert_eq!(moving.total_energy_j().unwrap().to_bits(), energy);
    assert!(core::ptr::eq(moving.normal(), &normal));
    assert!(core::ptr::eq(moving.shear(), &shear));
    let low = still
        .apply_source_refresh(&refresh, &next_shear, &gate())
        .unwrap();
    let high = moving
        .apply_source_refresh(&refresh, &next_shear, &gate())
        .unwrap();
    assert_eq!(high.energy_before_j, high.energy_after_j);
    assert_ne!(high.parameter_work_j, 0.0);
    assert_eq!(high.parameter_work_j, low.parameter_work_j);
    assert_eq!(high.energy_residual_j, low.energy_residual_j);
    assert_eq!(high.energy_tolerance_j, low.energy_tolerance_j);
    assert!(high.energy_residual_j.abs() <= high.energy_tolerance_j);
    assert!(high.energy_tolerance_j < 1e-8 * high.parameter_work_j.abs());
    assert_eq!(states(&moving), before);
    assert_eq!(moving.samples_rendered(), 0);
    assert!(moving.last_frame().is_none());
    assert!(
        moving
            .apply_source_refresh(&refresh, &next_shear, &gate())
            .is_err()
    );
}
