//! G0/G3: complete-frame acceleration against the actual CPU instrument and receiver.
use super::*;
use crate::grand_piano::audio::AudioStream;
use crate::grand_piano::engine::external_step_test_stream;
use std::cell::{Cell, RefCell};
use std::rc::Rc;

const EMPTY_SCORE: &str = "sample,event,key,value\n";

fn stream(substeps: usize, score: &str, frames: usize) -> AudioStream {
    external_step_test_stream(substeps, score, frames as u64)
}

fn piano(instrument: &mut Instrument) -> Piano {
    let (bank, dampers) = instrument.gpu_images().expect("substeps must be admitted");
    let mut program = Program::new(&bank, dampers.as_ref()).unwrap();
    let (lifted, sustain) = instrument.damper_controls();
    program
        .set_dampers(dampers.as_ref(), &lifted, sustain)
        .unwrap();
    Piano {
        backend: Box::new(Cpu(Emulator::new(&program))),
        required_backend: None,
        program,
        damper_key: Some((lifted, sustain.to_bits())),
        positions: Vec::new(),
        check: 0,
        stats: Stats::default(),
        started: Instant::now(),
    }
}

fn near(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).fold(0.0_f64, |m, x| m.max(x.abs()));
    let error = a
        .iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(a.iter().chain(b).all(|x| x.is_finite()));
    assert!(
        error <= 1e-12 + 1e-9 * scale,
        "error {error:e}, scale {scale:e}"
    );
}

fn state_near(a: &Instrument, b: &Instrument) {
    near(&a.bank.q, &b.bank.q);
    near(&a.bank.v, &b.bank.v);
    let ledger = |p: &Instrument| {
        let e = p.accounting;
        [
            e.input_work_j,
            e.modal_loss_j,
            e.damper_loss_j,
            e.felt_loss_j,
            e.felt_relaxation_loss_j,
            e.shank_loss_j,
            e.catch_loss_j,
            p.energy_j(),
        ]
    };
    near(&ledger(a), &ledger(b));
    for p in [a, b] {
        assert!(
            (p.accounting.input_work_j - p.accounting.dissipated_j() - p.energy_j()).abs() < 1e-7
        );
        assert!(p.accounting.max_balance_error_j < 1e-7);
    }
}

struct CountFrames {
    piano: Piano,
    frames: Rc<Cell<usize>>,
}
impl Accelerator for CountFrames {
    fn advance(
        &mut self,
        instrument: &mut Instrument,
        frames: usize,
        traces: &mut [f64],
    ) -> Result<usize, String> {
        let done = self.piano.advance(instrument, frames, traces)?;
        assert!(done <= frames);
        self.frames.set(self.frames.get() + done);
        Ok(done)
    }
}

fn accelerate(stream: &mut AudioStream) -> Rc<Cell<usize>> {
    let mut piano = piano(stream.instrument_mut());
    piano.check = 64; // Still mechanics steps, including at 4x and 16x.
    let frames = Rc::new(Cell::new(0));
    stream.set_test_accelerator(Box::new(CountFrames {
        piano,
        frames: frames.clone(),
    }));
    frames
}

#[test]
fn multirate_stream_matches_cpu_pressure_controls_and_final_energy() {
    let score = "sample,event,key,value\n0,note_on,69,2\n129,sustain,0,1\n257,note_off,69,0\n384,sustain,0,0.5\n511,sostenuto,0,1\n639,sustain,0,0\n767,silent_key_down,69,0\n899,note_off,69,0\n1025,sostenuto,0,0\n";
    for substeps in [1, 4, 16] {
        let frames = 1203;
        let mut actual = stream(substeps, score, frames);
        let mut reference = stream(substeps, score, frames);
        let accelerated = accelerate(&mut actual);
        let mut got = vec![0.0; frames];
        let mut expected = vec![0.0; frames];
        reference.render_block(&mut expected).unwrap();
        let mut start = 0;
        for size in [1, 7, 64, 0, 127, 13, 256].into_iter().cycle() {
            let end = (start + size).min(frames);
            actual.render_block(&mut got[start..end]).unwrap();
            start = end;
            if start == frames {
                break;
            }
        }
        assert!(
            accelerated.get() > 0,
            "must exercise the accelerator at {substeps}x"
        );
        assert!(expected.iter().any(|x| x.abs() > 1e-10));
        near(&got, &expected);
        assert_eq!(actual.sample_position(), frames as u64);
        state_near(actual.instrument(), reference.instrument());
    }
}

#[test]
fn events_at_accelerated_boundary_dispatch_before_that_frame_and_horizon_is_exact() {
    for frames in [1, 2, 3, 17, 513] {
        let event = frames - 1;
        let score = format!("sample,event,key,value\n{event},note_on,69,2\n");
        let mut actual = stream(4, &score, frames);
        let mut reference = stream(4, &score, frames);
        let accelerated = accelerate(&mut actual);
        let mut prefix = vec![0.0; event];
        actual.render_block(&mut prefix).unwrap();
        reference.render_block(&mut vec![0.0; event]).unwrap();
        assert_eq!(actual.instrument().accounting.input_work_j, 0.0);
        assert_eq!(actual.sample_position(), event as u64);
        actual.render_block(&mut [0.0]).unwrap();
        reference.render_block(&mut [0.0]).unwrap();
        assert!(actual.instrument().accounting.input_work_j > 0.0);
        assert_eq!(actual.sample_position(), frames as u64);
        state_near(actual.instrument(), reference.instrument());
        if frames > 2 {
            assert!(accelerated.get() > 0);
        }
    }
}

#[test]
fn short_requests_and_nonmultiple_check_caps_never_commit_partial_frames() {
    for substeps in [1, 4, 16] {
        for frames in [0, 1, 2, 3, 17] {
            for cap in [substeps - 1, substeps, 2 * substeps + 1] {
                let mut actual = stream(substeps, EMPTY_SCORE, 100);
                let mut reference = stream(substeps, EMPTY_SCORE, 100);
                actual.instrument_mut().note_on(69, 2.0).unwrap();
                reference.instrument_mut().note_on(69, 2.0).unwrap();
                let mut accelerator = piano(actual.instrument_mut());
                accelerator.check = cap;
                let mut traces = vec![f64::NAN; frames * actual.instrument().board_trace_len()];
                let done = accelerator
                    .advance(actual.instrument_mut(), frames, &mut traces)
                    .unwrap();
                let expected = if cap == 0 {
                    frames
                } else {
                    frames.min(cap / substeps)
                };
                assert_eq!(done, expected);
                let len = reference.instrument().board_trace_len();
                for frame in 0..done {
                    let mut trace = vec![0.0; len];
                    reference
                        .instrument_mut()
                        .step_with_board_trace(&mut trace)
                        .unwrap();
                    near(&traces[frame * len..(frame + 1) * len], &trace);
                }
                assert!(traces[done * len..].iter().all(|x| x.is_nan()));
                state_near(actual.instrument(), reference.instrument());
            }
        }
    }
}

#[test]
fn planner_prefix_rounds_down_to_whole_output_frames() {
    let substeps = 4;
    let mut selected = None;
    for velocity in [2.0, 3.0, 5.0, 7.1, 8.0] {
        let mut actual = stream(substeps, EMPTY_SCORE, 1000);
        actual.instrument_mut().note_on(69, velocity).unwrap();
        let planned =
            actual
                .instrument_mut()
                .contact_free_steps(SEGMENT, PLAN_CLEARANCE_M, PLAN_MARGIN);
        if planned >= MIN_SEGMENT && !planned.is_multiple_of(substeps) {
            selected = Some((actual, velocity, planned));
            break;
        }
    }
    let (mut actual, velocity, planned) =
        selected.expect("fixture must end inside an output frame");
    let mut reference = stream(substeps, EMPTY_SCORE, 1000);
    reference.instrument_mut().note_on(69, velocity).unwrap();
    let mut accelerator = piano(actual.instrument_mut());
    let mut traces = vec![f64::NAN; (SEGMENT / substeps) * actual.instrument().board_trace_len()];
    let done = accelerator
        .advance(actual.instrument_mut(), SEGMENT / substeps, &mut traces)
        .unwrap();
    assert_eq!(done, planned / substeps);
    assert_eq!(accelerator.stats.steps, done * substeps);
    for _ in 0..done {
        reference.instrument_mut().step().unwrap();
    }
    state_near(actual.instrument(), reference.instrument());
    // The first CPU frame after the shortened stretch must start at the same jack/contact state.
    for _ in 0..32 {
        actual.instrument_mut().step().unwrap();
        reference.instrument_mut().step().unwrap();
    }
    state_near(actual.instrument(), reference.instrument());
}

/// Run the real emulator, changing only a reported position (contact refusal)
/// or the reported loss (ledger refusal) to exercise the production rollback.
struct FaultBackend {
    cpu: Cpu,
    refuse_step: Option<usize>,
    bad_loss: bool,
    backend_error: bool,
    calls: Rc<RefCell<Vec<usize>>>,
}
impl Backend for FaultBackend {
    fn name(&self) -> String {
        "CPU emulator with injected refusal".into()
    }
    fn load_dampers(&mut self, p: &Program) -> Result<(), String> {
        self.cpu.load_dampers(p)
    }
    fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> Result<(f64, f64), String> {
        let mut losses = self.cpu.run(p, q, v, steps, trace, positions)?;
        let first = self.calls.borrow().is_empty();
        self.calls.borrow_mut().push(steps);
        if first {
            if let Some(step) = self.refuse_step {
                positions[step * p.contacts.len()] = -1.0;
            }
            if self.bad_loss {
                losses.0 += 1.0;
            }
            if self.backend_error {
                return Err("injected failure after real emulator work".into());
            }
        }
        Ok(losses)
    }
}

#[test]
fn contact_refusal_at_every_substep_restores_whole_frame_and_only_reruns_complete_prefix() {
    for frame in [0, 3] {
        for offset in 0..4 {
            let mut actual = stream(4, EMPTY_SCORE, 1000);
            let mut reference = stream(4, EMPTY_SCORE, 1000);
            actual.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
            reference.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
            let before = format!("{:?}", actual.instrument().checkpoint_external_hammers());
            let mut accelerator = piano(actual.instrument_mut());
            let calls = Rc::new(RefCell::new(Vec::new()));
            accelerator.backend = Box::new(FaultBackend {
                cpu: Cpu(Emulator::new(&accelerator.program)),
                refuse_step: Some(frame * 4 + offset),
                bad_loss: false,
                backend_error: false,
                calls: calls.clone(),
            });
            let mut traces = vec![0.0; 8 * actual.instrument().board_trace_len()];
            let done = accelerator
                .advance(actual.instrument_mut(), 8, &mut traces)
                .unwrap();
            assert_eq!(done, frame);
            assert_eq!(accelerator.stats.rollbacks, 1);
            assert_eq!(accelerator.stats.steps, frame * 4);
            assert_eq!(
                *calls.borrow(),
                if frame == 0 {
                    vec![32]
                } else {
                    vec![32, frame * 4]
                }
            );
            if frame == 0 {
                assert_eq!(
                    before,
                    format!("{:?}", actual.instrument().checkpoint_external_hammers())
                );
            }
            for _ in 0..done {
                reference.instrument_mut().step().unwrap();
            }
            state_near(actual.instrument(), reference.instrument());
            // Repeating the rejected frame on the CPU must not double-count jack time or felt/work.
            for _ in 0..64 {
                actual.instrument_mut().step().unwrap();
                reference.instrument_mut().step().unwrap();
            }
            state_near(actual.instrument(), reference.instrument());
        }
    }
}

#[test]
fn ledger_refusal_restores_the_entire_unpublished_stretch() {
    let mut actual = stream(4, EMPTY_SCORE, 1000);
    let mut reference = stream(4, EMPTY_SCORE, 1000);
    actual.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
    reference.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
    let before = format!("{:?}", actual.instrument().checkpoint_external_hammers());
    let q = actual.instrument().bank.q.clone();
    let v = actual.instrument().bank.v.clone();
    let modal = actual.instrument().bank.last_modal_loss_j;
    let mut accelerator = piano(actual.instrument_mut());
    accelerator.backend = Box::new(FaultBackend {
        cpu: Cpu(Emulator::new(&accelerator.program)),
        refuse_step: None,
        bad_loss: true,
        backend_error: false,
        calls: Rc::new(RefCell::new(Vec::new())),
    });
    let mut traces = vec![0.0; 8 * actual.instrument().board_trace_len()];
    assert!(
        accelerator
            .advance(actual.instrument_mut(), 8, &mut traces)
            .unwrap_err()
            .contains("energy ledger")
    );
    assert_eq!(q, actual.instrument().bank.q);
    assert_eq!(v, actual.instrument().bank.v);
    assert_eq!(modal, actual.instrument().bank.last_modal_loss_j);
    assert_eq!(
        before,
        format!("{:?}", actual.instrument().checkpoint_external_hammers())
    );
    for _ in 0..64 {
        actual.instrument_mut().step().unwrap();
        reference.instrument_mut().step().unwrap();
    }
    state_near(actual.instrument(), reference.instrument());
}

/// G1/G3: a device can fail after partially writing the returned state. Execute
/// actual emulator work before injecting that error, then compare whole-frame
/// state/accounting with the independent CPU instrument.
#[test]
fn required_gpu_failure_restores_the_frame_without_cpu_fallback() {
    for required in ["cuda", "metal"] {
        let mut actual = stream(4, EMPTY_SCORE, 1000);
        let mut reference = stream(4, EMPTY_SCORE, 1000);
        actual.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
        reference.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
        let before = format!("{:?}", actual.instrument().checkpoint_external_hammers());
        let q = actual.instrument().bank.q.clone();
        let v = actual.instrument().bank.v.clone();
        let mut accelerator = piano(actual.instrument_mut());
        accelerator.required_backend = Some(required.into());
        let calls = Rc::new(RefCell::new(Vec::new()));
        accelerator.backend = Box::new(FaultBackend {
            cpu: Cpu(Emulator::new(&accelerator.program)),
            refuse_step: None,
            bad_loss: false,
            backend_error: true,
            calls: calls.clone(),
        });
        let mut traces = vec![0.0; 8 * actual.instrument().board_trace_len()];
        let error = accelerator.advance(actual.instrument_mut(), 8, &mut traces).unwrap_err();
        assert!(error.contains(&format!("required {required} backend failed")), "{error}");
        assert_eq!(*calls.borrow(), vec![32]);
        assert_eq!(accelerator.backend.name(), "CPU emulator with injected refusal");
        assert_eq!(accelerator.stats.steps, 0);
        assert_eq!(q, actual.instrument().bank.q);
        assert_eq!(v, actual.instrument().bank.v);
        assert_eq!(before, format!("{:?}", actual.instrument().checkpoint_external_hammers()));
        state_near(actual.instrument(), reference.instrument());
        for _ in 0..64 {
            actual.instrument_mut().step().unwrap();
            reference.instrument_mut().step().unwrap();
        }
        state_near(actual.instrument(), reference.instrument());
    }
}

#[test]
fn automatic_backend_failure_keeps_the_real_cpu_trajectory() {
    let mut actual = stream(4, EMPTY_SCORE, 1000);
    let mut reference = stream(4, EMPTY_SCORE, 1000);
    actual.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
    reference.instrument_mut().jack_on(69, 10.0, 0.007).unwrap();
    let mut accelerator = piano(actual.instrument_mut());
    let calls = Rc::new(RefCell::new(Vec::new()));
    accelerator.backend = Box::new(FaultBackend {
        cpu: Cpu(Emulator::new(&accelerator.program)),
        refuse_step: None,
        bad_loss: false,
        backend_error: true,
        calls: calls.clone(),
    });
    let len = actual.instrument().board_trace_len();
    let mut traces = vec![0.0; 8 * len];
    assert_eq!(accelerator.advance(actual.instrument_mut(), 8, &mut traces).unwrap(), 8);
    assert_eq!(*calls.borrow(), vec![32]);
    for frame in 0..8 {
        let mut trace = vec![0.0; len];
        reference.instrument_mut().step_with_board_trace(&mut trace).unwrap();
        near(&traces[frame * len..(frame + 1) * len], &trace);
    }
    state_near(actual.instrument(), reference.instrument());
}
