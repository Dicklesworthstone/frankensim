//! Characterisation of the committed duduk card's plain-bore stand-in (bead rc-1109, duduk first render).
//!
//! This test RECORDS a known miss; it is not a tolerance of the instrument and not a claim that the
//! engine reproduces the duduk. The card (data/instruments/geometry/duduk-a-simon-maugeais-dalmont2024)
//! carries the reed region as a plain 12 mm bore because the reed interior is unpublished. Against the
//! three first peaks printed in the CC BY legend of Maugeais and Dalmont 2024, Figure 7 (G-1 180.8 Hz,
//! G 392.2 Hz, C+1 576.8 Hz), the passive input-impedance peak of this stand-in lies +28.9, +19.6 and
//! +53.8 cents sharp (fs-duct Bessel losses, dry air at 293.15 K and 101325 Pa, unflanged open end, all
//! ten 9.0 mm holes without undercut). For comparison, the Figure 6 round interior (22.07 cm3) puts these
//! peaks 250 to 450 cents flat and the interior derived by perimeter conservation (8.83 cm3) 88 to 348
//! cents sharp; the plain bore (11.65 cm3) is the one the card uses.
//!
//! Error model of the pinned values: the peak is refined on a 0.01 Hz grid (0.03 to 0.06 cents), the
//! recorded cents carry 0.05 cents of rounding, and the same lines were identical on two hosts to the
//! printed digits; the pin is therefore 0.3 cents. A failure means the engine or the card changed: read
//! the card's unknown u-reed-interior before moving a number.

use fs_duct::{input_impedance, Duct, HoleState, LossModel, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use std::path::PathBuf;

fn card_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/instruments/geometry/duduk-a-simon-maugeais-dalmont2024")
}

/// Rows of a committed specimen table, after its schema line and its header line.
fn rows(name: &str) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(card_dir().join(name)).expect("read the committed card table");
    text.lines()
        .skip(2)
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn dimension(part: &str, feature: &str, quantity: &str) -> f64 {
    rows("dimensions.tsv")
        .into_iter()
        .find(|r| r[0] == part && r[1] == feature && r[2] == quantity)
        .unwrap_or_else(|| panic!("card row {part} {feature} {quantity}"))[3]
        .parse()
        .expect("card value")
}

/// The card's duct for one fingering, from the reed tip: a plain bore split at each hole, as the lowering does.
fn duct(fingering: &str) -> Duct {
    let open: Vec<usize> = rows("fingerings.tsv")
        .into_iter()
        .find(|r| r[0] == fingering)
        .unwrap_or_else(|| panic!("fingering {fingering}"))[4]
        .split(',')
        .filter_map(|n| n.parse().ok())
        .collect();
    let d0 = dimension("bore", "station-000", "diameter");
    let d1 = dimension("bore", "station-001", "diameter");
    assert_eq!(d0, d1, "the card is no longer a plain bore: this characterisation does not apply");
    let radius = d0 / 2.0;
    let end = dimension("bore", "station-001", "axial-position") - dimension("bore", "station-000", "axial-position");
    let mut segments = Vec::new();
    let mut at = dimension("bore", "station-000", "axial-position");
    for number in 1..=10usize {
        let name = format!("hole-{number:02}");
        let x = dimension("tone-holes", &name, "axial-position");
        segments.push(Segment::Cylinder { radius, length: x - at });
        segments.push(Segment::ToneHole {
            hole_radius: dimension("tone-holes", &name, "hole-diameter") / 2.0,
            chimney_height: dimension("tone-holes", &name, "chimney-height"),
            bore_radius: radius,
            state: if open.contains(&number) { HoleState::Open } else { HoleState::Closed },
        });
        at = x;
    }
    segments.push(Segment::Cylinder { radius, length: end - at });
    Duct { segments }
}

/// First peak of |Z_in| within [0.85, 1.25] x the legend frequency, refined to 0.01 Hz.
fn first_peak(duct: &Duct, legend_hz: f64) -> f64 {
    let state = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("gas");
    let magnitude = |f: f64| -> f64 {
        input_impedance(duct, &state, std::f64::consts::TAU * f, LossModel::Bessel, Termination::UnflangedOpen)
            .expect("impedance")
            .impedance
            .abs()
    };
    let (lo, hi, step) = (0.85 * legend_hz, 1.25 * legend_hz, 0.25);
    let n = ((hi - lo) / step) as usize;
    let coarse: Vec<f64> = (0..=n).map(|i| magnitude(lo + step * i as f64)).collect();
    let top = (1..n)
        .filter(|&i| coarse[i] > coarse[i - 1] && coarse[i] >= coarse[i + 1])
        .max_by(|&a, &b| coarse[a].total_cmp(&coarse[b]))
        .expect("an interior maximum inside the window");
    let centre = lo + step * top as f64;
    let (mut best_f, mut best) = (centre, coarse[top]);
    let mut f = centre - 1.5 * step;
    while f <= centre + 1.5 * step {
        let m = magnitude(f);
        if m > best {
            best = m;
            best_f = f;
        }
        f += 0.01;
    }
    best_f
}

#[test]
fn plain_bore_stand_in_lies_sharp_of_the_figure_7_legend_peaks_by_the_recorded_cents() {
    // (fingering id, legend peak of Figure 7 [Hz], recorded residual of the stand-in [cents])
    let recorded = [("gm1", 180.8, 28.9), ("g", 392.2, 19.6), ("cp1", 576.8, 53.8)];
    for (fingering, legend, expected_cents) in recorded {
        let peak = first_peak(&duct(fingering), legend);
        let cents = 1200.0 * (peak / legend).log2();
        assert!(
            (cents - expected_cents).abs() < 0.3,
            "{fingering}: first peak {peak:.2} Hz is {cents:.2} cents from the legend {legend} Hz, recorded {expected_cents}"
        );
    }
}
