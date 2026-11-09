//! G1/G3/G4: actual body recoil, finite snap heat and one accepted lifecycle.

#[path = "native/dissipation.rs"]
mod dissipation;

use super::{CRITICAL_SHIFT, binding, controls as heat_controls, curves, events, source, thermal};
use fs_contact::normal_patch::adhesive::{JkrBranch, rolling::events::SourceJkrRollingEvents};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{JkrRollingSnapHeatError, native::*},
};
use fs_exec::CancelGate;
use fs_mbd::{
    MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3,
    tangential_friction::{
        RigidContactBinding,
        adhesive::rolling::{RigidAdhesiveRollingControls, events::*},
    },
};
use fs_thermal::enthalpy::{CellId, DepositionError};
use fs_tribo::partial_slip::TangentFrame;

const DT: f64 = 1e-4;

fn controls() -> RigidJkrRollingSnapHeatControls {
    RigidJkrRollingSnapHeatControls {
        mechanical: RigidAdhesiveRollingControls {
            maximum_bodies: 4,
            maximum_step_s: 0.01,
            maximum_angle_increment_rad: 1e-4,
            maximum_iterations: 96,
            angle_tolerance_rad: 1e-18,
            geometry_tolerance_m: 1e-15,
            angular_velocity_tolerance_rad_s: 1e-10,
            linear_momentum_tolerance_ns: 1e-25,
            angular_momentum_tolerance_nms: 1e-25,
            energy_tolerance_j: 1e-20,
        },
        heat: heat_controls(),
    }
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(
        &self,
        source: &'a SourceJkrRollingEvents,
    ) -> RigidAdhesiveRollingEventContact<'a> {
        RigidAdhesiveRollingEventContact {
            binding: self.binding,
            frame: self.frame,
            source,
        }
    }
}

fn fixture(
    source: &SourceJkrRollingEvents,
    detached: bool,
) -> (Fixture, RigidAdhesiveRollingEventState) {
    let properties = [
        MassProperties::from_principal_frame(
            2e-5,
            Vec3::ZERO,
            Vec3::new(0.8e-9, 1e-9, 1.2e-9),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3e-5,
            Vec3::ZERO,
            Vec3::new(1.5e-9, 2e-9, 2.5e-9),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
    ];
    let arms = [Vec3::new(0.0, 0.0, 0.01), Vec3::new(0.0, 0.0, -0.01)];
    let rate = 3.0 * CRITICAL_SHIFT / source.normal_source().effective_radius_m() / DT;
    let mut f = Fixture {
        bodies: core::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(arms[i].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                properties[i]
                    .angular_momentum_body_checked(Vec3::new(
                        if i == 0 { rate } else { 0.0 },
                        0.0,
                        0.0,
                    ))
                    .unwrap(),
            )
            .unwrap()
        }),
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    };
    let contact = f.contact(source);
    let gate = CancelGate::new();
    let mut state = if detached {
        contact.initial_detached_state(
            &f.bodies,
            &f.properties,
            -1e-6,
            controls().mechanical,
            &gate,
        )
    } else {
        contact.initial_attached_state(&f.bodies, &f.properties, 0.0, controls().mechanical, &gate)
    }
    .unwrap();
    if !detached {
        // The preload comes from actual finite-body recoil, never history adoption.
        apply_rigid_adhesive_rolling_events(
            &mut f.bodies,
            &f.properties,
            &mut state,
            contact,
            0.0,
            DT,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        assert!(state.source_state().rolling_storage_j() > 0.0);
    }
    (f, state)
}

#[test]
fn native_recoil_separation_and_recontact_publish_the_same_source_and_actual_heat() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let mut heat = binding(
        events.clone(),
        state.source_state().clone(),
        &caloric,
        &thermal,
    );
    let gate = CancelGate::new();
    for (phase, delta) in [0.0, normal.minimum_attached_indentation_m(), 0.0]
        .into_iter()
        .enumerate()
    {
        let before = f.bodies;
        let cells = thermal.cells().clone();
        let old_storage = state.source_state().normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j();
        let expected = prepare_rigid_adhesive_rolling_events(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&events),
            delta,
            DT,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let report = heat
            .apply_native(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                delta,
                DT,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(f.bodies, expected.bodies_after);
        assert_eq!(state, expected.state_after);
        assert_eq!(report.mechanics, expected.report);
        assert_eq!(heat.state(), state.source_state());
        assert_eq!(report.heat.final_state, *state.source_state());
        let kinetic: f64 = (0..2)
            .map(|i| {
                before[i]
                    .kinetic_energy_change_to(f.properties[i], f.bodies[i])
                    .unwrap()
            })
            .sum();
        let contact = state.source_state().normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j()
            - old_storage;
        let actual_heat: f64 = [CellId(0), CellId(1)]
            .into_iter()
            .map(|id| {
                let cell = &thermal.cells()[&id];
                cell.mass_kg * (cell.specific_h_j_kg - cells[&id].specific_h_j_kg)
            })
            .sum();
        let actual_residual = kinetic
            + contact
            + actual_heat
            + report.heat.work.other_work_j()
            + report.mechanics.irreversible_rolling_loss_j
            - report.mechanics.applied_work_j;
        assert!(actual_residual.abs() <= 2e-20);
        assert!(
            report.combined_energy_residual_j.abs()
                <= controls().mechanical.energy_tolerance_j
                    + controls().heat.transfer_energy_tolerance_j
        );
        assert_eq!(thermal.time_s(), 0.0);
        if phase == 0 {
            assert_ne!(f.bodies, before);
            assert!(report.mechanics.irreversible_rolling_loss_j > 0.0);
            assert_eq!(report.heat.work.dissipated_work_j(), 0.0);
            assert_eq!(thermal.cells(), &cells);
        } else if phase == 1 {
            assert_eq!(f.bodies, before);
            assert!(report.mechanics.mechanical.is_none());
            assert_eq!(
                state.source_state().normal_response().branch,
                JkrBranch::Detached
            );
            assert!(report.heat.rolling_release.dissipated_work_j() > 0.0);
            assert!(actual_heat > 0.0);
        } else {
            assert_ne!(f.bodies, before);
            let mechanical = report.mechanics.mechanical.as_ref().unwrap();
            assert_eq!(mechanical.contact.rolling().storage_before_j, 0.0);
            assert_eq!(report.heat.rolling_release.dissipated_work_j(), 0.0);
            assert!(report.heat.normal_release.dissipated_work_j() > 0.0);
        }
    }
}

#[test]
fn caloric_refusal_final_cancellation_and_stale_replay_preserve_all_four_owners() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, true);
    let original_bodies = f.bodies;
    let original_state = state.clone();
    let caloric = curves(&normal, 1.0);
    let mut narrow = curves(&normal, 1e-5);
    narrow[0] = caloric[0].clone();
    let mut refused_thermal = thermal(&narrow, 1.0);
    let cells = refused_thermal.cells().clone();
    let mut refused_heat = binding(
        events.clone(),
        state.source_state().clone(),
        &narrow,
        &refused_thermal,
    );
    let stamp = refused_heat.stamp();
    let gate = CancelGate::new();
    let mechanical = prepare_rigid_adhesive_rolling_events(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&events),
        0.0,
        DT,
        controls().mechanical,
        &gate,
    )
    .unwrap();
    assert_ne!(mechanical.bodies_after, original_bodies);
    assert!(matches!(
        refused_heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut refused_thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(()),
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
                DepositionError::MaterialDomain
            ))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(refused_heat.stamp(), stamp);
    assert_eq!(refused_thermal.cells(), &cells);

    let mut thermal = thermal(&caloric, 1.0);
    let original_cells = thermal.cells().clone();
    let mut heat = binding(
        events.clone(),
        state.source_state().clone(),
        &caloric,
        &thermal,
    );
    let stamp = heat.stamp();
    let mismatch = f
        .contact(&events)
        .initial_attached_state(&f.bodies, &f.properties, 0.0, controls().mechanical, &gate)
        .unwrap();
    assert!(matches!(
        heat.propose_native(
            &f.bodies,
            &f.properties,
            &mismatch,
            f.binding,
            f.frame,
            &thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrRollingSnapHeatError::HistoryMismatch)
    ));
    let mut calls = 0;
    let expected = heat
        .propose_native(
            &f.bodies,
            &f.properties,
            &state,
            f.binding,
            f.frame,
            &thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || {
                polls += 1;
                if polls == calls + 1 {
                    gate.request();
                }
                Ok(())
            },
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(DepositionError::Cancelled))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), &original_cells);
    let gate = CancelGate::new();
    let retry = heat
        .apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(&retry, expected.report());
    assert_eq!(&f.bodies, expected.bodies_after());
    assert_eq!(&state, expected.native_state());
    assert_eq!(heat.stamp(), expected.history().stamp());
    assert_eq!(thermal.cells(), expected.thermal().cells());
    let accepted_bodies = f.bodies;
    let accepted_state = state.clone();
    let accepted_stamp = heat.stamp();
    let accepted_cells = thermal.cells().clone();
    assert!(matches!(
        heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp,
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrRollingSnapHeatError::HistoryMismatch)
    ));
    assert_eq!(f.bodies, accepted_bodies);
    assert_eq!(state, accepted_state);
    assert_eq!(heat.stamp(), accepted_stamp);
    assert_eq!(thermal.cells(), &accepted_cells);
}
