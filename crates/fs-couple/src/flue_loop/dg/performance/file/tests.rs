use super::*;

fn read(text: &str) -> Result<Option<CircularOutletReceiver>, FlueDgInputError> {
    let mut rows = Rows { lines: text.lines(), line: 0 };
    read_observation(&mut rows)
}

#[test]
fn observation_syntax_retains_the_explicit_receiver_and_one_record_boundary() {
    assert!(read("observation inlet-pressure").unwrap().is_none());
    let mut rows = Rows {
        lines: "observation baffled -0.2 0.1 0.3 8 32 2000\nnext-owner-record".lines(),
        line: 0,
    };
    let receiver = read_observation(&mut rows).unwrap().unwrap();
    assert_eq!(receiver.position_m, [-0.2, 0.1, 0.3]);
    assert_eq!(receiver.radial_rings, 8);
    assert_eq!(receiver.angular_points, 32);
    assert_eq!(receiver.maximum_frequency_hz, 2000.0);
    assert_eq!(rows.lines.next(), Some("next-owner-record"));
}

#[test]
fn missing_unknown_nonfinite_and_extra_observation_fields_refuse() {
    for text in [
        "",
        "observation",
        "receiver inlet-pressure",
        "observation exterior",
        "observation inlet-pressure extra",
        "observation baffled",
        "observation baffled 0 0 0.2 8 32",
        "observation baffled 0 0 0.2 8 32 2000 extra",
        "observation baffled NaN 0 0.2 8 32 2000",
        "observation baffled 0 0 inf 8 32 2000",
        "observation baffled 0 0 0.2 -1 32 2000",
        "observation baffled 0 0 0.2 8 1.5 2000",
        "observation baffled 0 0 0.2 8 32 NaN",
    ] {
        assert!(read(text).is_err(), "admitted {text}");
    }
}

#[test]
fn mouth_records_reuse_explicit_geometry_and_preserve_the_following_record() {
    let mut rows = Rows {
        lines: "mouth_radiation -0.03 0.01 16 8\nmouth_reaction compact-baffled\nnext-owner-record".lines(),
        line: 0,
    };
    let (mouth, reaction) = read_mouth(&mut rows, true).unwrap();
    let mouth = mouth.unwrap();
    assert_eq!(mouth.center_m, [-0.03, 0.01]);
    assert_eq!(mouth.cells, [16, 8]);
    assert_eq!(reaction, MouthReaction::CompactBaffled);
    assert_eq!(rows.lines.next(), Some("next-owner-record"));

    let mut rows = Rows { lines: "next-owner-record".lines(), line: 0 };
    assert_eq!(read_mouth(&mut rows, false).unwrap(), (None, MouthReaction::None));
    assert_eq!(rows.lines.next(), Some("next-owner-record"));
}

#[test]
fn independent_acoustic_aperture_uses_the_shared_parser_and_retains_record_boundaries() {
    let mut rows = Rows {
        lines: concat!("mouth_aperture rectangle 0.012 0.008\n",
            "mouth_radiation -0.03 0.01 16 8\n",
            "mouth_reaction mutual-baffled\nnext-owner-record").lines(),
        line: 0,
    };
    let (mouth, reaction) = read_mouth(&mut rows, true).unwrap();
    let mouth = mouth.unwrap();
    assert_eq!(mouth.aperture_size_m, Some([0.012, 0.008]));
    assert_eq!(mouth.center_m, [-0.03, 0.01]);
    assert_eq!(mouth.cells, [16, 8]);
    assert_eq!(reaction, MouthReaction::MutualBaffled);
    assert_eq!(rows.lines.next(), Some("next-owner-record"));
}

#[test]
fn mouth_records_refuse_missing_exterior_missing_rectangle_and_ambiguous_fields() {
    for (text, exterior) in [
        ("mouth_radiation -0.03 0 16 8", false),
        ("mouth_aperture rectangle 0.012 0.008", true),
        ("mouth_aperture circle 0.012 0.008\nmouth_radiation -0.03 0 16 8", true),
        ("mouth_aperture rectangle 0 0.008\nmouth_radiation -0.03 0 16 8", true),
        ("mouth_aperture rectangle NaN 0.008\nmouth_radiation -0.03 0 16 8", true),
        ("mouth_aperture rectangle 0.012 0.008 extra\nmouth_radiation -0.03 0 16 8", true),
        ("mouth_aperture rectangle 0.012 0.008\nmouth_radiation -0.03 0 16 8", false),
        ("mouth_reaction compact-baffled", true),
        ("mouth_reaction mutual-baffled", true),
        ("mouth_radiation -0.03 0 16", true),
        ("mouth_radiation -0.03 0 16 8 extra", true),
        ("mouth_radiation NaN 0 16 8", true),
        ("mouth_radiation -0.03 0 -1 8", true),
        ("mouth_radiation -0.03 0 16 8\nmouth_reaction unknown", true),
        ("mouth_radiation -0.03 0 16 8\nmouth_reaction compact-baffled extra", true),
        ("mouth_radiation -0.03 0 16 8\nmouth_reaction mutual-baffled extra", true),
    ] {
        let mut rows = Rows { lines: text.lines(), line: 0 };
        assert!(read_mouth(&mut rows, exterior).is_err(), "admitted {text}");
    }
}


#[test]
fn shared_mouth_selector_distinguishes_joint_feedback_from_compact_and_observation_only() {
    for (token, expected) in [
        ("none", MouthReaction::None),
        ("compact-baffled", MouthReaction::CompactBaffled),
        ("mutual-baffled", MouthReaction::MutualBaffled),
    ] {
        let text = format!("mouth_radiation -0.003 0 8 8\nmouth_reaction {token}\nnext-owner-record");
        let mut rows = Rows { lines: text.lines(), line: 0 };
        let (mouth, reaction) = read_mouth(&mut rows, true).unwrap();
        assert_eq!(mouth.unwrap().cells, [8, 8]);
        assert_eq!(reaction, expected);
        assert_eq!(rows.lines.next(), Some("next-owner-record"));
    }
}
