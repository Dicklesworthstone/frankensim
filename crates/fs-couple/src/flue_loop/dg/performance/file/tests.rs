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
