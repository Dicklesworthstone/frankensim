//! The Model D board's loss is DATA (pexfp.373). A `damping,material` board
//! takes every stiffness share's loss factor from a loss pack, projects the
//! shares' strain energies once on the source modes, and damps each mode by
//! its own energy-weighted loss. A mode outside a strict claim's measured
//! frequencies refuses by name; nothing falls back to the authored 0.015.
//! Swapping the panel's along-grain claim moves each mode by its energy
//! fraction times the change; drier air moves the panel's along-grain loss
//! along the Sitka moisture curve, and every mode's loss with it.
use super::super::board_geometry::{BoardGeometry, ritz::RitzOptions,
    loss::{CLUSTER_PAIR_RATIO, DIAGONAL_TOLERANCE, MaterialLossReceipt, ShareLoss}};
use super::*;
use fs_matdb::PropertyValue;

const KEYS: [u8; 88] = {
    let mut keys = [0u8; 88];
    let mut i = 0;
    while i < 88 {
        keys[i] = 21 + i as u8;
        i += 1;
    }
    keys
};
const OPEN: LossPolicy = LossPolicy { extrapolate_frequency: true, sugar_pine_loss_unknown: true };
const AUTHORED_RATIO: f64 = 0.015;

fn humid(temperature_k: f64, relative_humidity: f64) -> BoardEnvironment {
    BoardEnvironment::from_gas(&GasState::try_new_moist_air(temperature_k, 101_325.0, relative_humidity).unwrap())
        .unwrap()
}
fn prepared(text: &str, upper_hz: f64) -> (super::super::board_geometry::PreparedBoard, MaterialLossReceipt) {
    let (board, receipt) = BoardGeometry::read(text).unwrap().prepare_with_loss_receipt(&KEYS, upper_hz).unwrap();
    (board, receipt.expect("a damping,material board returns its receipt"))
}
fn eta_of(receipt: &MaterialLossReceipt, share: &str) -> f64 {
    match &receipt.shares.iter().find(|(n, _)| n == share).unwrap().1 {
        ShareLoss::Claim { eta, .. } => *eta,
        ShareLoss::Unknown { .. } => 0.0,
    }
}
fn index_of(receipt: &MaterialLossReceipt, share: &str) -> usize {
    receipt.shares.iter().position(|(n, _)| n == share).unwrap()
}
/// The Sitka moisture curve read straight from its pack: linear between the
/// moisture rows that bracket `mc` (a fraction).
fn sitka_curve(mc: f64) -> f64 {
    let pack = load_pack(SITKA_LOSS_MOISTURE).unwrap();
    let claims = pack.claims().claims_for("loss_factor_longitudinal");
    let PropertyValue::Curve { knots, .. } = &claims[0].1.value else { panic!("curve") };
    let i = knots.iter().position(|k| k.0 >= mc).unwrap();
    let (a, b) = (knots[i - 1], knots[i]);
    a.1 + (b.1 - a.1) * (mc - a.0) / (b.0 - a.0)
}

// Planted negative (DONE-WHEN): strict claims measured at other frequencies
// refuse the board by name at its first mode; sugar pine without its flag
// refuses before any board exists; nothing falls back to 0.015.
#[test]
fn a_mode_outside_a_strict_claims_frequencies_refuses_by_name_and_nothing_falls_back() {
    let materials = reference_materials().unwrap();
    let strict = LossPolicy { extrapolate_frequency: false, sugar_pine_loss_unknown: true };
    let preset = build_with_material_loss(8, true, &materials, strict).unwrap();
    assert!(preset.geometry.contains("\ndamping,material\n") && !preset.geometry.contains("damping,0.015"));
    let error = BoardGeometry::read(&preset.geometry).unwrap().prepare(&KEYS, 600.0).unwrap_err();
    assert!(error.contains("board mode 0 (") && error.contains("lies outside the 62.36-68.34 Hz")
        && error.contains("'panel-l'") && error.contains("declare extrapolate-frequency")
        && error.contains("no authored damping replaces it"), "{error}");
    let error = model_d_loss_rows(&materials.environment,
        LossPolicy { extrapolate_frequency: true, sugar_pine_loss_unknown: false }).unwrap_err();
    assert!(error.contains("sugar pine") && error.contains("Wegst 2006") && error.contains("sugar_pine_loss_unknown"),
        "{error}");
    // A material board must declare every share, and loss rows need damping,material.
    let rows = model_d_loss_rows(&materials.environment, OPEN).unwrap();
    let text = build_with_material_loss(8, true, &materials, OPEN).unwrap().geometry;
    let without = text.lines().filter(|l| !l.starts_with("loss-share,bridges-torsion")).collect::<Vec<_>>().join("\n");
    assert!(BoardGeometry::read(&without).unwrap_err().contains("'bridges-torsion' is undeclared"));
    let authored = build_with_materials(8, true, &materials).unwrap().geometry;
    let mixed = authored.replacen("damping,0.015\n", &format!("damping,0.015\n{rows}"), 1);
    assert!(BoardGeometry::read(&mixed).unwrap_err().contains("loss rows belong to a damping,material board"));
}

// G1 (DONE-WHEN regression): each mode's damping is its energy-weighted
// material loss, eta_a = sum_c eta_c f_c,a over shares that partition its
// strain energy; the receipt names the unknown sugar pine and Poisson shares
// as a lower bound, every extrapolated claim, the measured coupling, and the
// disagreement with the authored 0.015 per mode and per band.
#[test]
fn each_mode_is_damped_by_its_energy_weighted_material_loss_and_the_receipt_says_what_is_missing() {
    let materials = reference_materials().unwrap();
    let text = build_with_material_loss(8, true, &materials, OPEN).unwrap().geometry;
    let (board, receipt) = prepared(&text, 600.0);
    assert_eq!(board.modes.len(), receipt.modes.len());
    assert!(board.provenance.contains("damping from material loss shares"));
    let unknown = ["panel-12", "ribs-bending", "ribs-torsion", "cutoff-bar-bending", "cutoff-bar-torsion"];
    let eta_l = sitka_curve(materials.environment.moisture_content);
    assert!((eta_of(&receipt, "panel-l") / eta_l - 1.0).abs() < 1e-12, "{} vs {eta_l}", eta_of(&receipt, "panel-l"));
    for (a, (mode, m)) in board.modes.iter().zip(&receipt.modes).enumerate() {
        // The shares partition K, so the fractions sum to phi^T K phi / lambda: one up to the
        // board's own certificates, its mass-normalization admission (1e-7) and the mode's
        // certified eigenvalue interval.
        let (lo, hi) = board.frequency_intervals_hz[a];
        let bound = 1e-7 + (hi * hi - lo * lo) / (mode.frequency_hz * mode.frequency_hz);
        let sum: f64 = m.fractions.iter().sum();
        assert!((sum - 1.0).abs() <= bound, "mode {a}: fractions sum {sum}, certified bound {bound:e}");
        let predicted: f64 = receipt.shares.iter().zip(&m.fractions).map(|((name, _), f)| eta_of(&receipt, name) * f).sum();
        assert!((m.loss_factor / predicted - 1.0).abs() < 1e-12, "mode {a}");
        assert!((mode.damping_ratio / (m.loss_factor / 2.0) - 1.0).abs() < 1e-12, "mode {a}");
        let missing: f64 = unknown.iter().map(|s| m.fractions[index_of(&receipt, s)]).sum();
        assert!((m.unknown_fraction - missing).abs() < 1e-15, "mode {a}");
        assert!(m.extrapolated.contains(&index_of(&receipt, "panel-l")), "mode {a}: panel-l is measured at 62-68 Hz");
        // Below Spycher's first point the sensitivity is the loss itself.
        if m.frequency_hz <= receipt.sensitivity[0].from_hz {
            assert!((m.sensitivity[0] / m.loss_factor - 1.0).abs() < 1e-12, "mode {a}");
        }
        assert!(m.max_coupling.is_finite() && m.neighbours.iter().flatten().all(|c| c.ratio.is_finite()));
    }
    assert!(receipt.caveats.iter().any(|c| c.contains("'ribs-bending' loss unknown") && c.contains("lower bound")));
    assert!(receipt.caveats.iter().any(|c| c.contains("'panel-l' loss factor") && c.contains("outside its measured")));
    println!("{}", receipt.tsv());
    for band in receipt.band_summary(&[0.0, 200.0, 400.0, 600.0]) {
        println!("band {:.0}-{:.0} Hz: {} modes, damping ratio mean {:.5} (min {:.5}, max {:.5}) = {:.3} x the authored {AUTHORED_RATIO}",
            band.band_hz.0, band.band_hz.1, band.modes, band.mean_ratio, band.min_ratio, band.max_ratio,
            band.mean_ratio / AUTHORED_RATIO);
    }
    for caveat in &receipt.caveats { println!("caveat: {caveat}"); }
}

// G1 (DONE-WHEN positive): replacing the panel's along-grain claim (Qiu 2026
// at the board's EMC) by Ono 1996's (board 1, 1049.4 Hz) moves every mode's
// loss factor by exactly the change times that mode's panel-L energy fraction;
// the elastic board and its modes are unchanged.
#[test]
fn swapping_the_panel_along_grain_claim_moves_each_mode_by_its_energy_fraction() {
    let materials = reference_materials().unwrap();
    let text = build_with_material_loss(8, true, &materials, OPEN).unwrap().geometry;
    let ono = load_pack(SITKA_LOSS).unwrap();
    let replacement = pinned_share("panel-l", &ono, "loss_factor_longitudinal", 1049.4, "extrapolate-frequency",
        "board 1 of two, swapped in for the test").unwrap();
    let swapped: String = text.lines().map(|line| if line.starts_with("loss-share,panel-l,") {
        replacement.trim_end().to_string() } else { line.to_string() }).collect::<Vec<_>>().join("\n");
    assert_ne!(swapped, text);
    let ((before_board, before), (after_board, after)) = (prepared(&text, 600.0), prepared(&swapped, 600.0));
    assert_eq!(before_board.frequency_intervals_hz, after_board.frequency_intervals_hz);
    let (old, new) = (eta_of(&before, "panel-l"), eta_of(&after, "panel-l"));
    assert!((new - 0.00667).abs() < 1e-12 && new < old, "{old} -> {new}");
    let l = index_of(&before, "panel-l");
    for (a, (x, y)) in before.modes.iter().zip(&after.modes).enumerate() {
        assert_eq!(x.fractions, y.fractions, "mode {a}: same stiffness, same split");
        let predicted = (new - old) * x.fractions[l];
        assert!(((y.loss_factor - x.loss_factor) - predicted).abs() <= 1e-12 * predicted.abs() + 1e-16, "mode {a}");
        assert!(y.loss_factor < x.loss_factor, "mode {a}");
        println!("mode {a} {:.1} Hz: f_L {:.4}, eta {:.5} -> {:.5} (predicted {:.5})", x.frequency_hz, x.fractions[l],
            x.loss_factor, y.loss_factor, x.loss_factor + predicted);
    }
}

// G1 (DONE-WHEN positive): RH 0.65 -> 0.40 at 293.15 K takes the wood from
// EMC 0.120 to 0.077. The panel's along-grain loss moves along the Sitka
// moisture curve (its pack, read here), and every mode's loss with it by that
// change times its panel-L fraction; the other shares have no moisture law and
// are held. The board's own moduli also change, which shifts the fractions;
// that part is computed and printed beside the moisture-law part, not hidden.
#[test]
fn drier_air_moves_each_mode_as_the_sitka_moisture_curve_predicts() {
    let (wet, dry) = (humid(293.15, 0.65), humid(293.15, 0.40));
    let board = |environment: &BoardEnvironment| {
        let materials = model_d_materials(environment).unwrap();
        prepared(&build_with_material_loss(8, true, &materials, OPEN).unwrap().geometry, 450.0).1
    };
    let (before, after) = (board(&wet), board(&dry));
    let (old, new) = (eta_of(&before, "panel-l"), eta_of(&after, "panel-l"));
    assert!((old / sitka_curve(wet.moisture_content) - 1.0).abs() < 1e-12);
    assert!((new / sitka_curve(dry.moisture_content) - 1.0).abs() < 1e-12);
    assert!(new < old, "drier wood loses less along the grain: {old} -> {new}");
    for share in ["panel-r", "panel-lr", "bridges-bending", "bridges-torsion"] {
        assert_eq!(eta_of(&before, share), eta_of(&after, share), "{share} has no moisture law");
    }
    let l = index_of(&before, "panel-l");
    let common = before.modes.len().min(after.modes.len()).min(10);
    assert!(common >= 5, "{common} common modes");
    for a in 0..common {
        let (x, y) = (&before.modes[a], &after.modes[a]);
        let law = (new - old) * y.fractions[l];
        // The rest of the change comes only from the fraction shift.
        let shift: f64 = before.shares.iter().enumerate().map(|(c, (name, _))|
            eta_of(&before, name) * (y.fractions[c] - x.fractions[c])).sum();
        assert!(((y.loss_factor - x.loss_factor) - (law + shift)).abs() < 1e-12 * x.loss_factor, "mode {a}");
        assert!(y.loss_factor < x.loss_factor, "mode {a}");
        println!("mode {a}: {:.1} -> {:.1} Hz, eta {:.5} -> {:.5}: moisture law {law:+.2e}, fraction shift {shift:+.2e}",
            x.frequency_hz, y.frequency_hz, x.loss_factor, y.loss_factor);
    }
}

// G1 (DONE-WHEN regression): the shipped resolution and slice. Every source
// mode of the 2400 Hz slice and every reduced coordinate is printed against
// the authored 0.015, with the per-band summary and the near-pair coupling the
// diagonal neglects. Ruling 6431: each coupled cluster's exact eigenvalues
// measure the diagonal's error per mode against the declared 5 %; the result
// is printed, never tuned. The shipped .bin is untouched by this test.
#[test]
fn the_shipped_slice_reports_material_damping_against_the_authored_ratio() {
    let materials = reference_materials().unwrap();
    let preset = build_resolved_with_material_loss(16, 0.03, true, &materials, OPEN).unwrap();
    let options = RitzOptions::parse("pod,32,16,16,8000").unwrap();
    let (board, receipt) = BoardGeometry::read(&preset.geometry).unwrap()
        .prepare_reduced_with_loss_receipt(&KEYS, 2400.0, false, false, false, false, &options).unwrap();
    let receipt = receipt.unwrap();
    assert!(board.physical_damping.is_some());
    println!("{}", receipt.tsv());
    for band in receipt.band_summary(&[0.0, 200.0, 400.0, 800.0, 1600.0, 2400.0]) {
        println!("band {:.0}-{:.0} Hz: {} source modes, damping ratio mean {:.5} (min {:.5}, max {:.5}) = {:.3} x the authored {AUTHORED_RATIO}",
            band.band_hz.0, band.band_hz.1, band.modes, band.mean_ratio, band.min_ratio, band.max_ratio,
            band.mean_ratio / AUTHORED_RATIO);
    }
    let worst = receipt.modes.iter().enumerate().flat_map(|(a, m)| m.neighbours.iter().flatten().map(move |c| (a, c.mode, c.ratio)))
        .max_by(|x, y| x.2.total_cmp(&y.2)).unwrap();
    println!("largest near-pair coupling |H_ab|/|w_a^2 - w_b^2|: modes {} and {}: {:.3e}", worst.0, worst.1, worst.2);
    for (k, cluster) in receipt.clusters.iter().enumerate() {
        assert!(cluster.len() >= 2 && cluster.iter().all(|&a| receipt.modes[a].cluster == Some(k)), "cluster {k}");
        for &a in cluster {
            let m = &receipt.modes[a];
            println!("cluster {k} mode {a} {:.2} Hz: diagonal eta {:.5}, exact {:.5} ({:+.2}%), overlap {:.3}, \
                left out {:.1e}", m.frequency_hz, m.loss_factor, m.exact_loss_factor, 100.0 * m.diagonal_error(),
                m.overlap, m.neglected);
        }
    }
    assert!(receipt.max_outside_ratio < CLUSTER_PAIR_RATIO, "{}", receipt.max_outside_ratio);
    for (a, m) in receipt.modes.iter().enumerate() {
        assert!(m.exact_loss_factor > 0.0 && m.exact_loss_factor.is_finite(), "mode {a}");
        if m.cluster.is_none() { assert_eq!(m.exact_loss_factor, m.loss_factor, "mode {a}"); }
    }
    let (above, largest, at) = receipt.diagonal_check();
    println!("diagonal check: {} clusters, largest |eta_exact/eta_diagonal - 1| {:.3}% at mode {at} ({:.1} Hz); \
        {} modes above the declared {DIAGONAL_TOLERANCE}: {above:?}", receipt.clusters.len(), 100.0 * largest,
        receipt.modes[at].frequency_hz, above.len());
    assert!(receipt.caveats.iter().any(|c| c.starts_with("diagonal loss factor against the exact eigenvalues")));
    assert!(receipt.caveats.iter().any(|c| c.starts_with("off-diagonal loss H_ab takes sqrt") && c.contains("exact here")));
    for (i, mode) in board.modes.iter().enumerate() {
        println!("coordinate {i:2}: {:9.3} Hz, damping ratio {:.5} = {:.3} x {AUTHORED_RATIO}", mode.frequency_hz,
            mode.damping_ratio, mode.damping_ratio / AUTHORED_RATIO);
    }
    for caveat in &receipt.caveats { println!("caveat: {caveat}"); }
    assert_eq!(receipt.modes.len(), board.reduction.as_ref().unwrap().source_modes);
    assert!(receipt.modes.iter().all(|m| m.loss_factor > 0.0 && m.loss_factor.is_finite()));
}
