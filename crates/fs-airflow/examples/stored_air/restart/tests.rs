use super::*;
use super::super::{distance, evolution_for, history, mesh, parse, with_cx};

fn tick(run: &mut AdaptiveEvolution<State>, model: &Model<'_>, cx: &fs_exec::Cx<'_>, n: usize) -> AdaptiveReport {
    run.advance(n, &mut |old, interval, x| model.trial(cx, old, interval, x),
        &mut |_, a, b, _| Ok(distance(a, b, 1e-4)), &mut || false).unwrap()
}
fn location(name: &str) -> PathBuf {
    static SERIAL: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let serial = SERIAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    std::env::temp_dir().join(format!("fs-stored-air-{name}-{}-{now}-{serial}", std::process::id()))
}

#[test]
fn serialized_fem_restarts_match_uninterrupted_fields_and_cumulative_work() {
    with_cx(|cx| {
        let mesh = mesh(1).unwrap(); let inputs = Inputs::default();
        let model = Model::new(cx, &mesh, inputs).unwrap();
        let identity = binding(&model, inputs, &[7;32]);
        let make = || evolution_for(&mesh, &model.duty, 0.5);
        let mut full = make(); let mut expected_work = Progress::default();
        expected_work.record(&tick(&mut full, &model, cx, 4096)).unwrap();
        assert!(full.is_complete());
        let mut split = make(); let mut work = Progress::default();
        for _ in 0..4096 {
            let report = tick(&mut split, &model, cx, 1); work.record(&report).unwrap();
            let checkpoint = bytes(&split, work, &identity).unwrap();
            let mut fresh = make();
            assert_eq!(restore(&mut fresh, &checkpoint, &identity).unwrap(), work);
            assert_eq!(fresh, split); split = fresh;
            if report.complete { break; }
        }
        assert_eq!(split, full); assert_eq!(work, expected_work);
        assert_eq!(split.state().source_input_j.to_bits(), full.state().source_input_j.to_bits());
    });
}

#[test]
fn changed_physics_estimator_forcing_mesh_and_executable_refuse_restart() {
    with_cx(|cx| {
        use fs_conduction::duty::{DutyCycle, DutySegment};
        let grid = mesh(1).unwrap(); let inputs = Inputs::default();
        let model = Model::new(cx, &grid, inputs).unwrap();
        let original = evolution_for(&grid, &model.duty, 0.5);
        let identity = binding(&model, inputs, &[7;32]);
        let saved = bytes(&original, Progress::default(), &identity).unwrap();
        let mut different = Vec::new();
        let other_inputs = Inputs { air_capacity_j_k: 0.04, ..inputs };
        let other = Model::new(cx, &grid, other_inputs).unwrap();
        different.push(binding(&other, other_inputs, &[7;32]));
        let other_inputs = Inputs { ventilation_w_k: 0.0, ..inputs };
        let other = Model::new(cx, &grid, other_inputs).unwrap();
        different.push(binding(&other, other_inputs, &[7;32]));
        different.push(binding(&model, Inputs { tolerance_k: 2e-4, ..inputs }, &[7;32]));
        different.push(binding(&model, inputs, &[8;32]));
        let duty = DutyCycle::new(vec![DutySegment::constant(5.0, 0.5).unwrap()]).unwrap();
        let other = Model::with_duty(cx, &grid, inputs, duty).unwrap();
        different.push(binding(&other, inputs, &[7;32]));
        let other_grid = mesh(2).unwrap(); let other = Model::new(cx, &other_grid, inputs).unwrap();
        different.push(binding(&other, inputs, &[7;32]));
        for id in different {
            let mut run = original.clone();
            assert!(restore(&mut run, &saved, &id).is_err()); assert_eq!(run, original);
        }
        // Additional work is explicitly allowed; it must not alter the model.
        assert_eq!(binding(&model, Inputs { attempts: 12, ..inputs }, &[7;32]), identity);
    });
}

#[test]
fn field_codec_refuses_truncation_wrong_mesh_and_nonphysical_values() {
    let state = State { solid_k: vec![300.0, 301.0], air_k: 300.5, net_input_j: 0.2, source_input_j: 0.3 };
    let payload = encode(&state, Progress::default()).unwrap();
    assert_eq!(decode(&payload, 2).unwrap().0, state);
    assert!(decode(&payload, 3).is_err());
    for n in 0..payload.len() { assert!(decode(&payload[..n], 2).is_err()); }
    for (offset, value) in [(40, f64::NAN), (40, -1.0), (56, 0.0), (64, f64::INFINITY), (72, -1.0)] {
        let mut bad = payload.clone(); bad[offset..offset+8].copy_from_slice(&value.to_le_bytes());
        assert!(decode(&bad, 2).is_err());
    }
}

#[test]
fn directory_resume_ignores_partial_writes_and_never_replaces_prior_generations() {
    with_cx(|cx| {
        let grid = mesh(1).unwrap(); let inputs = Inputs::default();
        let model = Model::new(cx, &grid, inputs).unwrap();
        let identity = binding(&model, inputs, &[7;32]);
        let mut run = evolution_for(&grid, &model.duty, 0.5);
        let directory = location("journal"); let mut writer = Writer::new(&directory).unwrap();
        let old = writer.publish(&run, Progress::default(), &identity).unwrap();
        let old_bytes = fs::read(&old).unwrap();
        assert!(Writer::new(&directory).is_err());
        let mut work = Progress::default(); work.record(&tick(&mut run, &model, cx, 1)).unwrap();
        let last = writer.publish(&run, work, &identity).unwrap();
        fs::write(directory.join("checkpoint-99999999.pending"), b"interrupted write").unwrap();
        assert_eq!(read(&directory).unwrap(), fs::read(&last).unwrap());
        assert_eq!(read(&old).unwrap(), old_bytes);
        // A corrupt completed generation is refused, not silently replaced by
        // an older state that would conceal lost work.
        fs::write(&last, b"corrupt").unwrap();
        let mut fresh = evolution_for(&grid, &model.duty, 0.5); let before = fresh.clone();
        assert!(restore(&mut fresh, &read(&directory).unwrap(), &identity).is_err());
        assert_eq!(fresh, before); assert_eq!(read(&old).unwrap(), old_bytes);
    });
}

#[test]
fn contradictory_work_counts_and_oversized_files_are_refused() {
    with_cx(|cx| {
        let grid = mesh(1).unwrap(); let inputs = Inputs::default();
        let model = Model::new(cx, &grid, inputs).unwrap();
        let identity = binding(&model, inputs, &[7;32]);
        let mut run = evolution_for(&grid, &model.duty, 0.5); let before = run.clone();
        let bad = bytes(&run, Progress { attempts: 0, evaluations: 1, rejected: 0 }, &identity).unwrap();
        assert!(restore(&mut run, &bad, &identity).is_err()); assert_eq!(run, before);
        let path = location("oversize"); File::create(&path).unwrap().set_len((LIMIT+1) as u64).unwrap();
        assert!(read(&path).is_err());
    });
}

#[test]
fn restart_options_are_explicit_and_duplicate_flags_refuse() {
    let args = ["--resume", "prior", "--checkpoint-dir", "next", "--attempts", "1"];
    let options = parse(&args.map(String::from)).unwrap();
    assert_eq!(options.resume.as_deref(), Some(Path::new("prior")));
    assert_eq!(options.checkpoint_dir.as_deref(), Some(Path::new("next")));
    assert_eq!(options.inputs.attempts, 1);
    assert!(parse(&["--resume", "a", "--resume", "b"].map(String::from)).is_err());
    assert!(parse(&["--checkpoint-dir"].map(String::from)).is_err());
    assert!(history::default_pulse().window_s() > 0.0);
}
