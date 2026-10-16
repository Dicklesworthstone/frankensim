//! Reusable-buffer transitions retain the owned diagnostic API's exact states,
//! energy/work/loss records, arbitrary-duration fallback and atomic refusals.

use fs_couple::modal_acoustic_time::{
    ModalAcousticFrame, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeError, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_couple::render::{ModalStringVoice, RenderContext, RenderError, RenderVoice};
use fs_math::c64::C64;

fn model(count: usize, damping: f64, budget: ModalAcousticTimeBudget) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        (0..count)
            .map(|index| ModalAcousticMode {
                angular_frequency_rad_s: core::f64::consts::TAU * (220.0 + 55.0 * index as f64),
                damping_ratio: damping,
                pressure_per_modal_velocity: C64::new(0.75, -0.125),
            })
            .collect(),
        budget,
    )
    .unwrap()
}

fn assert_frame_bits(actual: &ModalAcousticFrame, expected: &ModalAcousticFrame) {
    for (a, b) in [
        (actual.observer_pressure_pa, expected.observer_pressure_pa),
        (actual.total_modal_energy_j, expected.total_modal_energy_j),
        (actual.input_work_j, expected.input_work_j),
        (actual.viscous_dissipation_j, expected.viscous_dissipation_j),
        (actual.dissipation_roundoff_tolerance_j, expected.dissipation_roundoff_tolerance_j),
    ] {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(actual.modal_energy_j.len(), expected.modal_energy_j.len());
    for (a, b) in actual.modal_energy_j.iter().zip(&expected.modal_energy_j) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(actual.modal_viscous_losses.len(), expected.modal_viscous_losses.len());
    for (a, b) in actual.modal_viscous_losses.iter().zip(&expected.modal_viscous_losses) {
        assert_eq!(a.energy_j.to_bits(), b.energy_j.to_bits());
        assert_eq!(a.roundoff_tolerance_j.to_bits(), b.roundoff_tolerance_j.to_bits());
    }
}

#[test]
fn workspace_and_owned_frames_agree_in_every_damping_and_duration_regime() {
    for damping in [0.0, 0.03, 1.0 - 2e-8, 1.0, 1.0 + 2e-8, 2.0] {
        for fraction in [1.0, 0.5, 0.37] {
            let mut owned = model(3, damping, ModalAcousticTimeBudget::audible_reference());
            let mut buffered = owned.clone();
            let mut workspace = ModalAcousticWorkspace::new(&buffered);
            let dt = fraction * owned.sample_period_s();
            for step in 0..256 {
                let forces = if step < 80 { [0.5, -0.25, 0.125] } else { [0.0; 3] };
                let expected = owned.step_duration(&forces, dt).unwrap();
                let actual = buffered.step_duration_into(&forces, dt, &mut workspace).unwrap();
                assert_frame_bits(actual, &expected);
                assert_eq!(buffered.states(), owned.states());
            }
        }
    }
}

#[test]
fn borrowed_diagnostics_and_model_state_keep_their_allocated_buffers() {
    let mut model = model(32, 0.03, ModalAcousticTimeBudget::audible_reference());
    let mut workspace = ModalAcousticWorkspace::new(&model);
    let state_buffer = model.states().as_ptr();
    let first = model.step_into(&[0.5; 32], &mut workspace).unwrap();
    let energy_buffer = first.modal_energy_j.as_ptr();
    let loss_buffer = first.modal_viscous_losses.as_ptr();
    let energy_capacity = first.modal_energy_j.capacity();
    let loss_capacity = first.modal_viscous_losses.capacity();
    for _ in 0..1024 {
        let frame = model.step_into(&[0.0; 32], &mut workspace).unwrap();
        assert_eq!(frame.modal_energy_j.as_ptr(), energy_buffer);
        assert_eq!(frame.modal_viscous_losses.as_ptr(), loss_buffer);
        assert_eq!(frame.modal_energy_j.capacity(), energy_capacity);
        assert_eq!(frame.modal_viscous_losses.capacity(), loss_capacity);
        assert_eq!(model.states().as_ptr(), state_buffer);
    }
}

#[test]
fn workspace_size_and_malformed_input_refusals_preserve_state_and_allow_retry() {
    let mut model = model(2, 0.03, ModalAcousticTimeBudget::audible_reference());
    let other = model.clone();
    let small = ModalAcousticTimeModel::try_new(
        48_000, vec![model.modes()[0]], ModalAcousticTimeBudget::audible_reference(),
    ).unwrap();
    let mut wrong = ModalAcousticWorkspace::new(&small);
    assert!(matches!(model.step_into(&[0.5; 2], &mut wrong),
        Err(ModalAcousticTimeError::InvalidInput { .. })));
    assert_eq!(model.states(), other.states());
    let mut workspace = ModalAcousticWorkspace::new(&model);
    for forces in [&[][..], &[f64::NAN, 0.0], &[f64::INFINITY, 0.0]] {
        assert!(model.step_into(forces, &mut workspace).is_err());
        assert_eq!(model.states(), other.states());
    }
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.step_duration_into(&[0.0; 2], duration, &mut workspace).is_err());
        assert_eq!(model.states(), other.states());
    }
    let mut expected = other;
    assert_frame_bits(model.step_into(&[0.5; 2], &mut workspace).unwrap(),
        &expected.step(&[0.5; 2]).unwrap());
    assert_eq!(model.states(), expected.states());
}

#[test]
fn late_pressure_and_energy_budget_refusals_do_not_publish_candidate_states() {
    for tight in [
        ModalAcousticTimeBudget { maximum_abs_pressure_pa: 1e-12,
            ..ModalAcousticTimeBudget::audible_reference() },
        ModalAcousticTimeBudget { maximum_total_energy_j: 1e-18,
            ..ModalAcousticTimeBudget::audible_reference() },
        ModalAcousticTimeBudget { maximum_abs_velocity_m_sqrt_kg_per_s: 1e-12,
            ..ModalAcousticTimeBudget::audible_reference() },
        ModalAcousticTimeBudget { maximum_abs_displacement_m_sqrt_kg: 1e-16,
            ..ModalAcousticTimeBudget::audible_reference() },
    ] {
        let mut buffered = model(2, 0.03, tight);
        let mut owned = buffered.clone();
        let mut workspace = ModalAcousticWorkspace::new(&buffered);
        let expected = owned.step(&[1.0, 0.5]).unwrap_err();
        let actual = buffered.step_into(&[1.0, 0.5], &mut workspace).unwrap_err();
        assert_eq!(actual, expected);
        assert_eq!(buffered.states(), &[ModalAcousticState::default(); 2]);
        // Failed scratch is overwritten; it never contaminates the next sample.
        assert_eq!(buffered.step_into(&[0.0; 2], &mut workspace).unwrap().observer_pressure_pa, 0.0);
    }
}

#[test]
fn borrowed_step_reuses_restored_history_not_workspace_history() {
    let mut buffered = model(2, 0.03, ModalAcousticTimeBudget::audible_reference());
    let mut workspace = ModalAcousticWorkspace::new(&buffered);
    buffered.step_into(&[0.5; 2], &mut workspace).unwrap();
    let states = [ModalAcousticState {
        displacement_m_sqrt_kg: 1e-5,
        velocity_m_sqrt_kg_per_s: -0.03,
    }; 2];
    buffered.restore_states(&states).unwrap();
    let mut owned = buffered.clone();
    assert_frame_bits(buffered.step_into(&[0.0; 2], &mut workspace).unwrap(),
        &owned.step(&[0.0; 2]).unwrap());
    assert_eq!(buffered.states(), owned.states());
}

#[test]
fn modal_render_keeps_context_poisoning_after_a_physics_refusal() {
    let tight = ModalAcousticTimeBudget {
        maximum_abs_pressure_pa: 1e-12,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    let voice = ModalStringVoice::new(model(2, 0.03, tight), vec![1.0; 2]).unwrap();
    let mut context = RenderContext::new(vec![RenderVoice::ModalString(voice)], 8);
    assert!(matches!(context.block(&mut [0.0; 8]), Err(RenderError::Modal(_))));
    let mut untouched = [91.0; 8];
    assert!(matches!(context.block(&mut untouched), Err(RenderError::Poisoned)));
    assert_eq!(untouched, [91.0; 8]);
    assert_eq!(context.samples_rendered(), 0);
}
