//! The shared force compiler must preserve the original projector and runtime.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use fs_couple::render::schedule::force::program::{PhysicalForceInput, ProjectedForceProgram};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::{ControlDelta, RenderError};
use fs_math::c64::C64;

fn config() -> ForceRenderConfig {
    ForceRenderConfig { sample_rate_hz: 48_000, max_block: 128,
        max_events: 32, max_controls: 128, max_projection_terms: 4096 }
}
fn input() -> PhysicalForceInput {
    PhysicalForceInput::new(vec![vec![2.0, -1.0], vec![0.5, 3.0]], vec![0.5, -0.25]).unwrap()
}
fn event(sample: u64, port: usize, force_n: f64) -> ModalForceEvent {
    ModalForceEvent { sample, voice: 0, port, force_n }
}
fn events() -> Vec<ModalForceEvent> {
    vec![event(19, 0, 0.0), event(7, 0, 9.0), event(7, 1, 0.25),
        event(7, 0, 1.0), event(51, 1, 0.0)]
}
fn model() -> ModalAcousticTimeModel {
    // Synthetic numerical observation; not an instrument radiation model.
    ModalAcousticTimeModel::try_new(48_000, [500.0, 800.0].map(|omega| ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.02,
        pressure_per_modal_velocity: C64::new(0.0, 1.0),
    }).to_vec(), ModalAcousticTimeBudget::audible_reference()).unwrap()
}

#[test]
fn all_same_time_changes_sum_before_projection_and_release_only_the_named_port() {
    let p = ProjectedForceProgram::compile(vec![input()], events(), config()).unwrap();
    assert_eq!(p.initial_forces(), &[vec![0.875, -1.25]]);
    let mut held = p.initial_forces()[0].clone(); let mut next = 0;
    for sample in 0..80 {
        while next < p.events().len() && p.events()[next].sample == sample {
            let e = p.events()[next]; held[e.mode] = e.force_n_per_sqrt_kg; next += 1;
        }
        let (a, b) = if sample < 7 {(0.5, -0.25)} else if sample < 19 {(1.0, 0.25)}
            else if sample < 51 {(0.0, 0.25)} else {(0.0, 0.0)};
        assert_eq!(held, vec![2.0*a+0.5*b, -a+3.0*b]);
    }
}

#[test]
fn original_renderer_uses_identical_compiled_controls_and_direct_motion() {
    let p = ProjectedForceProgram::compile(vec![input()], events(), config()).unwrap();
    let mut direct = model(); let mut force = p.initial_forces()[0].clone(); let mut next = 0;
    let expected: Vec<_> = (0..96).map(|sample| {
        while next < p.events().len() && p.events()[next].sample == sample {
            let e = p.events()[next]; force[e.mode] = e.force_n_per_sqrt_kg; next += 1;
        }
        direct.step(&force).unwrap().observer_pressure_pa.to_bits()
    }).collect();
    for partition in [1, 37, 96] {
        let voice = ModalForceVoice::new(model(), vec![vec![2.0, -1.0], vec![0.5, 3.0]],
            vec![0.5, -0.25], ForceInitialization::RetainState).unwrap();
        let mut renderer = ScheduledRenderer::from_modal_forces(vec![voice], events(), config()).unwrap();
        assert_eq!(renderer.pending_controls().len(), p.events().len());
        for (c, e) in renderer.pending_controls().iter().zip(p.events()) {
            assert_eq!(c.sample, e.sample);
            assert!(matches!(c.delta, ControlDelta::SetModalForce { voice: 0, mode, force_n_per_sqrt_kg }
                if mode == e.mode && force_n_per_sqrt_kg.to_bits() == e.force_n_per_sqrt_kg.to_bits()));
        }
        let mut out = vec![0.0; 96];
        for block in out.chunks_mut(partition) { renderer.block(block).unwrap(); }
        assert_eq!(out.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), expected);
    }
}

#[test]
fn force_work_is_conjugate_and_source_images_do_not_overwrite_each_other() {
    let second = PhysicalForceInput::new(vec![vec![-4.0]], vec![2.0]).unwrap();
    let p = ProjectedForceProgram::compile(vec![input(), second], vec![], config()).unwrap();
    assert_eq!(p.initial_forces()[1], [-8.0]);
    let delta = [0.003, -0.002];
    let mechanical = p.initial_forces()[0].iter().zip(delta).map(|(g, d)| g*d).sum::<f64>();
    let physical = 0.5*(2.0*delta[0]-delta[1])-0.25*(0.5*delta[0]+3.0*delta[1]);
    assert!((mechanical-physical).abs() < 1e-17);
}

#[test]
fn initial_preload_force_is_not_replaced_by_a_sample_zero_release() {
    let p = ProjectedForceProgram::compile(vec![input()], vec![event(0,0,0.0),event(0,1,0.0)], config()).unwrap();
    assert_eq!(p.initial_forces()[0], [0.875,-1.25]);
    assert_eq!(p.events().len(), 2);
    assert!(p.events().iter().all(|e| e.sample == 0 && e.force_n_per_sqrt_kg == 0.0));
}

#[test]
fn every_authored_input_and_all_three_original_budgets_are_enforced() {
    for columns in [vec![], vec![vec![]], vec![vec![f64::NAN]], vec![vec![1.0],vec![1.0,2.0]]] {
        assert!(PhysicalForceInput::new(columns, vec![0.0]).is_err());
    }
    for bad in [event(u64::MAX,0,0.0), event(0,2,0.0),event(0,0,f64::INFINITY)] {
        assert!(ProjectedForceProgram::compile(vec![input()], vec![bad], config()).is_err());
    }
    let wrong = ModalForceEvent { voice: 1, ..event(0,0,0.0) };
    assert!(matches!(ProjectedForceProgram::compile(vec![input()],vec![wrong],config()),
        Err(RenderError::UnknownVoice { index: 1 })));
    for which in 0..3 {
        let mut c = config();
        match which {0=>c.max_events=0,1=>c.max_controls=0,_=>c.max_projection_terms=3}
        assert!(matches!(ProjectedForceProgram::compile(vec![input()],events(),c),Err(RenderError::Sizing{..})));
    }
    let overflowing = PhysicalForceInput::new(vec![vec![f64::MAX]], vec![2.0]).unwrap();
    assert!(ProjectedForceProgram::compile(vec![overflowing], vec![], config()).is_err());
}
