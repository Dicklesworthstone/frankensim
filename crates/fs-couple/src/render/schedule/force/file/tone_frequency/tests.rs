use super::*;

fn records() -> [PhysicalTone; 2] {
    [
        PhysicalTone {
            voice: 1,
            port: 0,
            omega: 200.0,
            cosine: 3.0,
            sine: 4.0,
        },
        PhysicalTone {
            voice: 0,
            port: 1,
            omega: 200.0,
            cosine: 0.6,
            sine: 0.8,
        },
    ]
}

fn read(text: &str, storage: usize) -> Result<Option<ToneFrequencyPlan>, String> {
    ToneFrequencyPlan::read_optional(
        &mut Reader {
            lines: text.lines(),
            line: 0,
        },
        &records(),
        16,
        1000,
        0.9,
        None,
        &mut { storage },
    )
}

#[test]
fn original_rows_and_accepted_cursor_preserve_the_integrated_phase_through_silence() {
    let mut plan = read(
        "frequency_sweeps 1\nsweep 1 4 1\nomega 0 100\nomega 4 300\nomega 8 0\nomega 12 200\n",
        8192,
    )
    .unwrap()
    .unwrap();
    assert_eq!(plan.count(), 1);
    assert!(!plan.contains(0));
    assert!(plan.contains(1));
    assert!(plan.step(0, 0, [1.0; 2]).unwrap().is_none());
    assert!(plan.step(1, 1, [1.0; 2]).is_err());
    assert!(plan.step(2, 0, [1.0; 2]).is_err());
    let omega = |sample: u64| {
        if sample <= 4 {
            100.0 + 50.0 * sample as f64
        } else if sample <= 8 {
            300.0 - 75.0 * (sample - 4) as f64
        } else if sample <= 12 {
            50.0 * (sample - 8) as f64
        } else {
            200.0
        }
    };
    let mut phase = 0.0;
    let mut c = 0.6;
    let mut s = 0.8;
    for sample in 0..16 {
        let gain = if (4..10).contains(&sample) {
            [0.0; 2]
        } else {
            [1.0; 2]
        };
        let step = plan.step(1, sample, gain).unwrap().unwrap();
        let endpoint = step.endpoints(c, s).unwrap();
        assert_eq!(
            endpoint,
            plan.step(1, sample, gain)
                .unwrap()
                .unwrap()
                .endpoints(c, s)
                .unwrap()
        );
        if gain[0] == 0.0 {
            assert_eq!(endpoint, [0.0; 4]);
        }
        (c, s) = step.advance_force_phase(c, s).unwrap();
        phase += 0.001 * (omega(sample) + omega(sample + 1)) * 0.5;
        let (sn, cs) = phase.sin_cos();
        assert!((c - (0.6 * cs + 0.8 * sn)).abs() < 2.0e-14);
        assert!((s - (0.8 * cs - 0.6 * sn)).abs() < 2.0e-14);
        // No candidate query changed the physical clock or knot cursor.
        assert!(plan.step(1, sample + 1, gain).is_err());
        plan.commit(sample + 1);
    }
    assert!(plan.step(1, 16, [1.0; 2]).is_err());
    let bound = plan.maximum_force_error_bound_n();
    let curve = plan.take(1).unwrap();
    assert_eq!(curve.force_error_bound_n(), bound);
    assert!(!plan.contains(1));
}

#[test]
fn force_error_admission_uses_original_newtons_and_all_gain_segments() {
    let curve = "frequency_sweeps 1\nsweep 1 2 1\nomega 0 100\nomega 16 400\n";
    let small = read(curve, 8192).unwrap().unwrap();
    let large = read(&curve.replace("sweep 1", "sweep 0"), 8192)
        .unwrap()
        .unwrap();
    assert!(
        (large.maximum_force_error_bound_n() - 5.0 * small.maximum_force_error_bound_n()).abs()
            < 1.0e-15
    );
    let allowance = 1.1 * small.maximum_force_error_bound_n();
    let text = curve.replace("sweep 1 2 1", &format!("sweep 1 2 {allowance}"));
    assert!(read(&text, 8192).is_ok());
    let mut reader = Reader {
        lines: "envelopes 1\nenvelope 1 3\ngain 0 0\ngain 8 2\ngain 16 -1\n".lines(),
        line: 0,
    };
    let mut envelopes = ToneEnvelopePlan::read_optional(&mut reader, 2, 16, &mut 8192)
        .unwrap()
        .unwrap();
    assert_eq!(envelopes.bounds(0).unwrap(), [1.0, 0.0]);
    assert_eq!(envelopes.bounds(1).unwrap(), [2.0, 0.375]);
    envelopes.commit(12);
    assert_eq!(envelopes.bounds(1).unwrap(), [2.0, 0.375]);
    assert!(
        ToneFrequencyPlan::read_optional(
            &mut Reader {
                lines: text.lines(),
                line: 0
            },
            &records(),
            16,
            1000,
            0.9,
            Some(&envelopes),
            &mut 8192,
        )
        .is_err()
    );
    let plan = read(curve, 8192).unwrap().unwrap();
    // A caller cannot replace the admitted gains with an unbounded ramp.
    assert!(plan.step(1, 0, [0.0, 1.0e9]).is_err());
    assert!(plan.step(1, 0, [1.0; 2]).is_ok());
    assert!(plan.step(1, 0, [f64::NAN, 1.0]).is_err());
    assert!(plan.step(1, 0, [1.0; 2]).is_ok());
}

#[test]
fn malformed_unbudgeted_and_guard_violating_sweeps_refuse_before_rendering() {
    let valid = "frequency_sweeps 1\nsweep 1 2 1\nomega 0 100\nomega 16 300\n";
    assert!(read(valid, 1).is_err());
    for invalid in [
        valid.replace("frequency_sweeps 1", "frequency_sweeps 3"),
        valid.replace("sweep 1 2 1", "sweep 2 2 1"),
        valid.replace("sweep 1 2 1", "sweep 1 0 1"),
        valid.replace("sweep 1 2 1", "sweep 1 65537 1"),
        valid.replace("sweep 1 2 1", "sweep 1 2 0"),
        valid.replace("sweep 1 2 1", "sweep 1 2 NaN"),
        valid.replace("sweep 1 2 1", "sweep 1 2 1e-15"),
        valid.replace("omega 0 100", "omega 1 100"),
        valid.replace("omega 16 300", "omega 0 300"),
        valid.replace("omega 16 300", "omega 17 300"),
        valid.replace("omega 16 300", "omega 16 -1"),
        valid.replace("omega 16 300", "omega 16 NaN"),
        valid.replace("omega 16 300", "omega 16 3000"),
        valid.replace("omega 16 300", "omega 16 300 extra"),
        "frequency_sweeps 2\nsweep 1 1 1\nomega 0 100\nsweep 1 1 1\nomega 0 100\n".into(),
    ] {
        assert!(read(&invalid, 8192).is_err(), "admitted {invalid}");
    }
    assert!(read("frequency_sweeps 0\n", 0).unwrap().is_none());
    let mut reader = Reader {
        lines: "broadband 1\n".lines(),
        line: 0,
    };
    assert!(
        ToneFrequencyPlan::read_optional(&mut reader, &records(), 16, 1000, 0.9, None, &mut 0,)
            .unwrap()
            .is_none()
    );
    assert_eq!(reader.lines.next(), Some("broadband 1"));
    // Carrier-equal controls retain exact affine amplitude reconstruction.
    let plan = read("frequency_sweeps 1\nsweep 0 1 1e-30\nomega 0 200\n", 8192)
        .unwrap()
        .unwrap();
    assert_eq!(plan.maximum_force_error_bound_n(), 0.0);
    assert_eq!(
        plan.step(0, 0, [1.0; 2])
            .unwrap()
            .unwrap()
            .endpoints(3.0, 4.0)
            .unwrap(),
        [3.0, 4.0, 3.0, 4.0]
    );
}
