use super::*;

fn score(body: &str) -> Result<Program, Error> {
    Program::parse(&format!("{}\n{body}", score::HEADER))
}
fn image(p: &Program) -> Vec<(f64, f64)> {
    p.knots.iter().map(|k| (k.time_s, k.force_n)).collect()
}
const HIT: &str = "shape,hit,0,0\nshape,hit,0.125,4\nshape,hit,0.25,0\n";

#[test]
fn tempo_mapped_rolls_preserve_physical_gesture_duration_and_legacy_csv_arithmetic() {
    let p = score(&format!("tempo,0,120\ntempo,2,60\n{HIT}roll,0,1,4,hit,1\n")).unwrap();
    let old = Program::parse("0,0\n0.125,4\n0.25,0\n0.5,0\n0.625,4\n0.75,0\n\
        1,0\n1.125,4\n1.25,0\n2,0\n2.125,4\n2.25,0\n").unwrap();
    assert_eq!(image(&p), image(&old));
    for count in [1, 3, 17, 1000] {
        for i in 0..count {
            let a = f64::from(i) * 2.5 / f64::from(count);
            let b = f64::from(i + 1) * 2.5 / f64::from(count);
            assert_eq!(p.average(a, b).to_bits(), old.average(a, b).to_bits());
        }
    }
    assert!(p.admit(0.125, 17, 1.0).is_err(), "the final gesture tail must fit");
    assert!(p.admit(0.125, 18, 1.0).is_ok());
}

#[test]
fn overlapping_signed_gestures_add_at_every_breakpoint_without_a_force_tail() {
    let p = score("tempo,0,60\nshape,p,0,0\nshape,p,0.25,4\nshape,p,0.5,0\n\
        shape,p,0.75,-2\nshape,p,1,0\nstroke,0.25,p,0.5\nstroke,0,p,1\n").unwrap();
    let one = Program::parse("0,0\n0.25,4\n0.5,0\n0.75,-2\n1,0").unwrap();
    for count in [1, 3, 19, 1024] {
        let mut impulse = 0.0;
        for i in 0..count {
            let a = f64::from(i) * 1.5 / f64::from(count);
            let b = f64::from(i + 1) * 1.5 / f64::from(count);
            let expected = one.average(a, b) + 0.5 * one.average(a - 0.25, b - 0.25);
            assert!((p.average(a, b) - expected).abs() < 1e-12);
            impulse += p.average(a, b) * (b - a);
        }
        assert!((impulse - 0.75).abs() < 1e-12);
    }
    assert_eq!(p.knots.first().unwrap().force_n, 0.0);
    assert_eq!(p.knots.last().unwrap().force_n, 0.0);
    assert_eq!(p.average(1.25, 10.0), 0.0);
    // Simultaneous authored gestures are summed, not deduplicated as notes.
    let p = score(&format!("tempo,0,60\n{HIT}stroke,0,hit,1\nstroke,0,hit,0.5")).unwrap();
    assert_eq!(image(&p), vec![(0.0, 0.0), (0.125, 6.0), (0.25, 0.0)]);
}

#[test]
fn undefined_nonfinite_discontinuous_and_ambiguous_scores_refuse() {
    let body = format!("tempo,0,60\n{HIT}stroke,0,hit,1\n");
    for bad in [
        body.replace("tempo,0,60\n", ""), body.replace("tempo,0,60", "tempo,1,60"),
        body.replace("tempo,0,60", "tempo,0,0"), body.replace("tempo,0,60", "tempo,0,NaN"),
        body.replace("tempo,0,60", "tempo,0,1e-320"),
        body.replace("tempo,0,60", "tempo,0,60\ntempo,0,120"),
        body.replace("shape,hit,0,0", "shape,hit,0,1"),
        body.replace("shape,hit,0,0", "shape,hit,0.01,0"),
        body.replace("shape,hit,0.25,0", "shape,hit,0.25,1"),
        body.replace("shape,hit,0.25,0", "shape,hit,0.125,0"),
        body.replace("shape,hit,0.125,4", "shape,hit,0.125,inf"),
        body.replace("stroke,0,hit,1", "stroke,0,missing,1"),
        body.replace("stroke,0,hit,1", "stroke,0,hit,-1"),
        body.replace("stroke,0,hit,1", "stroke,NaN,hit,1"),
        body.replace("stroke,0,hit,1", "roll,0,0,4,hit,1"),
        body.replace("stroke,0,hit,1", "roll,0,1,0,hit,1"),
        body.replace("stroke,0,hit,1", "roll,0,1,4.5,hit,1"),
        body.replace("stroke,0,hit,1", "stroke,1e300,hit,1"),
        body.replace("stroke,0,hit,1", "roll,1e20,1,2,hit,1"),
        body.replace("stroke,0,hit,1", "stroke,0,hit,1e308"),
        body.replace("stroke,0,hit,1", "unknown,0,hit,1"),
    ] { assert!(score(&bad).is_err(), "accepted {bad}"); }
    assert!(score("").is_err());
    assert!(Program::parse("frankensim-stick-score-v2\n0,0\n1,0").is_err());
    assert!(score(&format!("{body}\n{}", score::HEADER)).is_err());
    // Even unplayed force cards are admitted; a typo is not hidden by a rest.
    assert!(score(&format!("{body}shape,unused,0,1\nshape,unused,1,0")).is_err());
}

#[test]
fn expansion_overlap_and_combined_force_remain_bounded_without_truncation() {
    assert!(score(&format!("tempo,0,60\n{HIT}roll,0,1,40000,hit,1")).is_err());
    assert!(score(&format!("tempo,0,60\n{HIT}roll,0,1,22000,hit,1")).is_err());
    let error = score("tempo,0,60\nshape,p,0,0\nshape,p,0.5,1\nshape,p,1,0\n\
        roll,0,0.0001,2048,p,1").unwrap_err();
    assert!(error.to_string().contains("overlap-compilation"));
    let p = score("tempo,0,60\nshape,p,0,0\nshape,p,0.5,1e5\nshape,p,1,0\n\
        stroke,0,p,1\nstroke,0,p,1").unwrap();
    assert!(p.admit(0.001, 1000, 1e6).is_err());
}

#[test]
fn file_selection_keeps_both_hand_ports_and_existing_launch_controls() {
    use std::io::Write;
    let path = std::env::temp_dir().join(format!("frankensim-stick-score-{}-{}.score",
        std::process::id(), std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    let text = format!("# explicit force score\n{}\ntempo,0,60\n{HIT}stroke,0,hit,1", score::HEADER);
    let mut f = std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    f.write_all(text.as_bytes()).unwrap(); drop(f);
    let mut args = vec!["drum".into(), "--stick-force-file".into(), path.to_str().unwrap().into(),
        "--second-stick-force-file".into(), path.to_str().unwrap().into(),
        "--strike-speed-m-s".into(), "0".into()];
    let first = option(&mut args).unwrap().unwrap();
    let second = second_option(&mut args).unwrap().unwrap();
    assert_eq!(image(&first), image(&second));
    assert_eq!(args, vec!["drum", "--strike-speed-m-s", "0"]);
    let mut d = StickDrive::new_inputs(vec![
        Input { program: first, coordinate: 0, tip_weight: 2.0 },
        Input { program: second, coordinate: 2, tip_weight: 0.5 },
    ], 0.125, 2, 3).unwrap();
    let expected = d.forces(&[1.0, 3.0, -2.0]).unwrap().to_vec();
    assert_eq!(expected, [5.0, 3.0, -1.0]);
    assert!(d.forces(&[0.0]).is_err());
    assert_eq!(d.forces(&[1.0, 3.0, -2.0]).unwrap(), expected.as_slice());
    assert_eq!(d.accepted, 0);
}

#[test]
fn rhythmic_two_hand_programs_reproduce_direct_forces_in_the_actual_contact_solver() {
    use super::super::super::{drum_with_sticks, Stroke};
    use fs_exec::CancelGate;
    let first = Stroke { speed_m_s: 4.0, position_m: Some([0.06, 0.01]) };
    let second = Stroke { speed_m_s: 2.5, position_m: Some([-0.05, 0.02]) };
    let make = || drum_with_sticks(256, 2e-6, false, true, None, false,
        first, false, None, None, Some(second)).unwrap();
    let mut played = make(); let mut direct = make();
    let port = played.second_stick.unwrap(); let n = played.force.len();
    let physical = "0,0\n0.000128,0\n0.000256,2\n0.000512,0";
    let shape = "shape,p,0,0\nshape,p,0.000128,0\nshape,p,0.000256,2\nshape,p,0.000512,0\n";
    let inputs = vec![
        Input { program: score(&format!("tempo,0,120\n{shape}stroke,0,p,1")).unwrap(),
            coordinate: 0, tip_weight: played.stick_weight },
        Input { program: score(&format!("tempo,0,120\n{shape}stroke,0,p,0.5")).unwrap(),
            coordinate: port.coordinate, tip_weight: port.weight },
    ];
    played.system = played.system.with_stick_drives(inputs, 2e-6, 256, n).unwrap();
    let mut stages = StickDrive::new_inputs(vec![
        Input { program: Program::parse(physical).unwrap(), coordinate: 0, tip_weight: direct.stick_weight },
        Input { program: Program::parse(&physical.replace(",2", ",1")).unwrap(),
            coordinate: port.coordinate, tip_weight: port.weight },
    ], 2e-6, 256, n).unwrap();
    let gate = CancelGate::new_clock_free(); let mut work = 0.0;
    for tick in 0..256 {
        if tick == 80 {
            let state = played.system.state().to_vec();
            let mut invalid = vec![0.0; n]; invalid[port.coordinate] = 1e7;
            assert!(played.system.step(&invalid, &gate).is_err());
            assert_eq!(played.system.state(), state);
        }
        let a = played.system.step(&played.force, &gate).unwrap();
        let b = direct.system.step(stages.forces(&direct.force).unwrap(), &gate).unwrap();
        stages.accept();
        assert_eq!(played.system.state(), direct.system.state());
        assert_eq!(a.supplied_work_j, b.supplied_work_j);
        assert_eq!(a.dissipated_energy_j, b.dissipated_energy_j);
        assert!(a.balance_residual_j.abs() < 1e-7);
        work += a.supplied_work_j.abs();
    }
    assert!(work > 0.0);
    assert!((played.system.state()[2 * port.coordinate + 1] * port.weight - second.speed_m_s).abs() > 1e-8);
}
