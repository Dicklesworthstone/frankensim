//! Front-door admission for physical Rayleigh backreaction. Kernel and passive
//! energy-exchange tests live with their numerical owners.
use super::*;

fn options(extra: &[&str]) -> Result<Options, String> {
    let mut args: Vec<String> = ["--preset", "steinway-d", "--render",
        "rayleigh-admission-test-unwritten.wav"].into_iter().map(str::to_owned).collect();
    args.extend(extra.iter().map(|s| (*s).to_owned()));
    Options::parse(&args)
}

#[test]
fn load_is_explicit_and_front_receiver_can_share_one_sided_load_with_cavity() {
    assert!(options(&[]).unwrap().rayleigh_loading.is_none());
    let selected = options(&["--rayleigh-loading", "1,20,2000,65",
        "--cavity", "declared-cavity.fspc", "--microphone-right", "1.2,1,1"]).unwrap();
    let load = selected.rayleigh_loading.unwrap();
    assert_eq!(load.sides, 1);
    assert_eq!(load.samples, 65);
    assert_eq!(load.min_hz, 20.0);
    assert_eq!(load.max_hz, 2000.0);
    assert_eq!(selected.microphone_right, Some([1.2, 1.0, 1.0]));
    assert!(options(&["--rayleigh-loading", "2,20,2000,65"]).is_ok());
}

#[test]
fn load_rejects_missing_geometry_diagnostics_duplicate_sides_and_audio_band() {
    for extra in [
        vec!["--rayleigh-loading", "2,20,2000,65", "--cavity", "declared-cavity.fspc"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--diagnostic-volume"],
        vec!["--rayleigh-loading", "1,20,3600,65", "--sample-rate", "8000"],
        vec!["--rayleigh-loading", "1,20,3000,65", "--substeps", "1"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--rayleigh-loading", "2,20,2000,65"],
    ] {
        assert!(options(&extra).is_err(), "accepted incompatible {extra:?}");
    }
    for args in [
        vec!["--rayleigh-loading", "1,20,2000,65"],
        vec!["--preset", "steinway-d", "--rayleigh-loading", "1,20,2000,65"],
        vec!["--render", "rayleigh-admission-test-unwritten.wav", "--rayleigh-loading", "1,20,2000,65"],
        vec!["--board", "modal-board.csv", "--render", "rayleigh-admission-test-unwritten.wav", "--rayleigh-loading", "1,20,2000,65"],
    ] {
        assert!(Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>()).is_err());
    }
}

#[test]
fn primary_radiation_ports_are_explicit_bounded_and_require_the_load() {
    assert_eq!(options(&[]).unwrap().radiation_ports, None);
    for ports in ["1", "32"] {
        let selected = options(&["--rayleigh-loading", "1,20,2000,65",
            "--radiation-ports", ports]).unwrap();
        assert_eq!(selected.radiation_ports, Some(ports.parse().unwrap()));
    }
    for extra in [
        vec!["--radiation-ports", "1"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--radiation-ports", "0"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--radiation-ports", "33"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--radiation-ports", "1", "--radiation-ports", "2"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--radiation-ports", "NaN"],
        vec!["--rayleigh-loading", "1,20,2000,65", "--radiation-ports", "2", "--diagnostic-volume"],
    ] {
        assert!(options(&extra).is_err(), "accepted incompatible {extra:?}");
    }
}
