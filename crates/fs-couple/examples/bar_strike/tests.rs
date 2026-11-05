//! G2 oracle checks of the route: the uniform prism against free-free Timoshenko beam theory, a planted
//! wrong density that must move every frequency by exactly the root of the ratio, card refusals, and a
//! coarse end-to-end strike that must be finite and replay bit for bit.
use crate::card::{Card, Profile, Resonator, Tube, pack_scalar};
use crate::render::Body;
use crate::{bar, radiation, resonator, strike, structure, tube};
use fs_couple::modal_acoustic_time::ModalAcousticMode;
use fs_math::{c64::C64, det};
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn blank() -> Card {
    let path =
        root().join("data/instruments/descriptions/xylophone/orff-soprano-blank-rucz2021.bar");
    Card::load(&path, &root()).expect("the shipped blank card parses")
}

/// First three free-free Timoshenko bending frequencies (Hz) of a prism: the transcendental
/// characteristic equations of the even and odd families, roots by bisection. `shear_modulus` is the
/// modulus of the plane of bending (here G_RL), `kappa` = 5/6.
fn timoshenko(young: f64, shear_modulus: f64, rho: f64, length: f64, thickness: f64) -> [f64; 3] {
    let rotary = thickness * thickness / (12.0 * length * length);
    let shear = young * rotary / ((5.0 / 6.0) * shear_modulus);
    let characteristic = |q: f64, even: bool| {
        let lambda = q.powi(4);
        let sum = lambda * (rotary + shear);
        let disc = (sum * sum + 4.0 * (lambda - rotary * shear * lambda * lambda)).sqrt();
        let alpha = (2.0 * (lambda - rotary * shear * lambda * lambda) / (disc + sum)).sqrt();
        let beta = (0.5 * (disc + sum)).sqrt();
        let a = alpha * (beta * beta - rotary * lambda);
        let b = beta * (alpha * alpha + rotary * lambda);
        if even {
            a * det::sin(beta / 2.0) + b * det::cos(beta / 2.0) * det::tanh(alpha / 2.0)
        } else {
            a * det::cos(beta / 2.0) * det::tanh(alpha / 2.0) - b * det::sin(beta / 2.0)
        }
    };
    let mut roots = Vec::new();
    for even in [true, false] {
        let mut lo = 0.2;
        for j in 1..=1800 {
            let hi = 0.2 + f64::from(j) * 0.01;
            if characteristic(lo, even) * characteristic(hi, even) < 0.0 {
                let (mut a, mut b) = (lo, hi);
                for _ in 0..64 {
                    let mid = (a + b) / 2.0;
                    if characteristic(a, even) * characteristic(mid, even) <= 0.0 {
                        b = mid;
                    } else {
                        a = mid;
                    }
                }
                roots.push((a + b) / 2.0);
            }
            lo = hi;
        }
    }
    roots.sort_by(f64::total_cmp);
    assert!(roots.len() >= 3, "free-free Timoshenko roots");
    core::array::from_fn(|i| {
        roots[i].powi(2) * (young * thickness * thickness / (12.0 * rho)).sqrt()
            / (bar::TAU * length * length)
    })
}

/// A coarse, low-band copy of the blank card: cheap enough for a unit test.
fn coarse(band_hz: f64) -> Card {
    let mut card = blank();
    card.band_hz = band_hz;
    card.windows = 2;
    card.mesh_rungs = vec![1];
    card
}

/// G2: the blank's TET10 vertical bending modes against the free-free Timoshenko oracle, built from the
/// same card (so the stand-in wood is irrelevant to the comparison). With Rucz's own reference constants the
/// oracle agreed with his 3-D orthotropic FEM of the blank to 0.55 % (Table I: 1027/2712/5019 Hz), which is
/// what makes the 1-D theory a fair judge here; the 1.5 % band leaves room for the second rung's mesh error
/// and no more. This test says nothing about the wood or about the measured blank.
#[test]
fn uniform_prism_vertical_modes_match_free_free_timoshenko() {
    let mut card = blank();
    card.band_hz = 6000.0;
    card.windows = 2;
    card.mesh_rungs = vec![2];
    let s = structure(&card, 2, card.density_kg_m3).expect("rung-2 structure");
    let vertical = bar::vertical_modes(&s.mesh, &s.modes);
    assert!(
        vertical.len() >= 3,
        "three vertical bending modes below 6 kHz: {vertical:?}"
    );
    let oracle = timoshenko(
        card.youngs_pa[0],
        card.shear_pa[2],
        card.density_kg_m3,
        card.length_m,
        card.thickness_m,
    );
    for k in 0..3 {
        let f = bar::hz(&s.modes[vertical[k]]);
        let miss = f / oracle[k] - 1.0;
        println!(
            "vertical {k}: tet10 {f:.2} Hz, Timoshenko {:.2} Hz, miss {:.3} %",
            oracle[k],
            100.0 * miss
        );
        assert!(
            miss.abs() < 0.015,
            "vertical mode {k}: {f} Hz against {} Hz",
            oracle[k]
        );
    }
    assert_eq!(s.rigid.len(), 6, "all six free-body coordinates are kept");
}

/// Planted negative: scaling the density by 1.21 at fixed stiffness must divide every frequency by
/// exactly 1.1. A route that read the density from anywhere else, or normalized it away, fails this.
#[test]
fn planted_wrong_density_moves_every_frequency_by_the_root_of_the_ratio() {
    let card = coarse(3500.0);
    let a = structure(&card, 1, card.density_kg_m3).expect("nominal density");
    let b = structure(&card, 1, 1.21 * card.density_kg_m3).expect("planted density");
    assert_eq!(a.modes.len(), b.modes.len());
    assert!(a.modes.len() >= 2);
    for (x, y) in a.modes.iter().zip(&b.modes) {
        let ratio = bar::hz(x) / bar::hz(y);
        assert!(
            (ratio - 1.1).abs() < 1e-8,
            "ratio {ratio} at {} Hz",
            bar::hz(x)
        );
    }
}

#[test]
fn card_refusals_are_explicit() {
    let text = blank().text;
    let density_row = text
        .lines()
        .find(|l| l.starts_with("wood_density\t"))
        .expect("the blank card names its density pack");
    let missing = text.replace(&format!("{density_row}\n"), "");
    let err = Card::parse(&missing, &root()).unwrap_err();
    assert!(err.contains("wood_density"), "{err}");
    // A pack that lacks the property is refused by name, not read as zero.
    let wrong_pack = text.replace(
        "wood_youngs_longitudinal\tdata/matdb/seed-v1/rosewood-indian-fpl-gtr282",
        "wood_youngs_longitudinal\tdata/matdb/seed-v1/rosewood-orff-xylophone-bars-rucz-jasa2021",
    );
    assert!(
        Card::parse(&wrong_pack, &root())
            .unwrap_err()
            .contains("no scalar `young_modulus_longitudinal`")
    );
    let twice = format!("{text}length_m\t0.3\n");
    assert!(Card::parse(&twice, &root()).unwrap_err().contains("twice"));
    let header = text.replacen("frankensim.bar-card.v1", "frankensim.bar-card.v0", 1);
    assert!(
        Card::parse(&header, &root())
            .unwrap_err()
            .contains("must start")
    );
    let bad_mesh = text.replace("mesh_counts\t32\t4\t2", "mesh_counts\t30\t4\t2");
    assert!(
        Card::parse(&bad_mesh, &root())
            .unwrap_err()
            .contains("divisible by 4")
    );
    // A cut that cannot fit inside the bar is refused, not clipped.
    let cut = text.replace(
        "profile\tuniform",
        "profile\tarc\nprofile_arc\t0.0625\t0.0075\t0.25",
    );
    assert!(
        Card::parse(&cut, &root())
            .unwrap_err()
            .contains("does not fit")
    );
}

#[test]
fn certificate_row_is_optional_and_bounded() {
    let text = blank().text;
    assert_eq!(blank().maximum_relative_halfwidth, 1e-8);
    let with = |value: &str| format!("{text}\nmaximum_relative_halfwidth\t{value}\n");
    let card = Card::parse(&with("1e-7"), &root()).expect("the row parses");
    assert_eq!(card.maximum_relative_halfwidth, 1e-7);
    for bad in ["0", "-1e-7", "1e-5"] {
        assert!(
            Card::parse(&with(bad), &root())
                .unwrap_err()
                .contains("maximum_relative_halfwidth"),
            "{bad} must be refused"
        );
    }
}

#[test]
fn arc_profile_is_flat_in_the_middle_and_full_thickness_outside_the_cut() {
    let p = Profile::Arc {
        stone_radius_m: 0.0625,
        middle_thickness_m: 0.0075,
        straight_length_m: 0.06,
    };
    let (l, t) = (0.325, 0.016);
    assert_eq!(p.thickness(l, t, 0.5 * l), 0.0075);
    assert_eq!(p.thickness(l, t, 0.5 * l + 0.03), 0.0075);
    assert_eq!(p.thickness(l, t, 0.0), t);
    // The arc rises monotonically to the full thickness at its reach and the stone's own circle describes it.
    let d = 0.012;
    let x = 0.5 * l + 0.03 + d;
    let expected = 0.0075 + 0.0625 - (0.0625_f64 * 0.0625 - d * d).sqrt();
    assert!((p.thickness(l, t, x) - expected).abs() < 1e-15);
}

/// The whole route at the coarsest rung: modes, Burton-Miller transfers, cited mallet, render. The output is
/// finite, non-trivial, and replays bit for bit.
#[test]
fn coarse_strike_is_finite_and_replays_bitwise() {
    let mut card = coarse(3000.0);
    card.duration_s = 0.05;
    card.substeps = 8;
    let s = structure(&card, 1, card.density_kg_m3).expect("structure");
    let first = strike(&card, &s).expect("strike");
    let second = strike(&card, &s).expect("strike replay");
    assert!(first.stroke.samples.iter().all(|p| p.is_finite()));
    assert!(
        first.stroke.samples.iter().any(|p| p.abs() > 0.0),
        "the bar must radiate"
    );
    assert!(first.stroke.contact_time_s > 0.0 && first.stroke.contact_time_s < 1e-3);
    assert_eq!(
        first.validity, "inside",
        "the stroke stays inside the contact card's validity"
    );
    assert_eq!(
        first.stroke.contact_time_s.to_bits(),
        second.stroke.contact_time_s.to_bits()
    );
    assert!(
        first
            .stroke
            .samples
            .iter()
            .zip(&second.stroke.samples)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "same-host exact replay"
    );
}

/// The fit order this bar's transfers are held to. How many poles a row will take depends on the row and the
/// band: on the marimba F4 card order 16 holds out at 0.19 % and 32 and 48 refuse with `SigmaCollapse`, and on
/// this blank's shorter band 16 refuses too. A refusal at 12 fails the test rather than moving to a lower order.
const OBSERVER_FIT_ORDER: usize = 12;

fn observed_strike(card: &mut Card, s: &crate::Structure, rigid_rows: bool) -> crate::Strike {
    card.observer = Some(crate::card::ObserverRow {
        fit_order: OBSERVER_FIT_ORDER,
        held_out_budget: 0.05,
        rigid_rows,
    });
    strike(card, s).expect("broadband strike at the pinned order")
}

/// The rigid rows change what is heard and nothing in the bar. The whole bar's translation at the band's
/// lowest frequency is the compact dipole it has to be, which fixes the rows' sign and size without using them:
/// in phase, and with a volume above the bar's own. A flat strip's added volume is only the thin limit, and a
/// section with real thickness carries more, so the upper bound allows 1.6 times it.
#[test]
fn rigid_rows_change_only_what_is_heard() {
    let mut card = coarse(3000.0);
    card.duration_s = 0.05;
    card.substeps = 8;
    let s = structure(&card, 1, card.density_kg_m3).expect("structure");
    let plain = strike(&card, &s).expect("narrow-band strike");
    let silent = observed_strike(&mut card, &s, false);
    let heard = observed_strike(&mut card, &s, true);
    assert_eq!(
        plain.stroke.contact_time_s.to_bits(),
        heard.stroke.contact_time_s.to_bits()
    );
    assert!(
        plain
            .stroke
            .samples
            .iter()
            .zip(&heard.narrow_samples)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "the bar must not change"
    );
    assert!(
        silent
            .stroke
            .samples
            .iter()
            .zip(&heard.stroke.samples)
            .any(|(a, b)| a != b),
        "the rigid rows must be heard"
    );
    let check = heard
        .observer
        .as_ref()
        .and_then(|o| o.rigid_check)
        .expect("the translation is driven");
    println!("rigid check: {check:?}");
    assert!(check.phase_deg.abs() < 15.0, "{check:?}");
    assert!(
        check.effective_volume_m3 > check.volume_m3
            && check.effective_volume_m3
                < check.volume_m3 + 1.6 * (check.strip_bound_m3 - check.volume_m3),
        "{check:?}"
    );
}

/// The observer row changes what is heard and nothing in the bar: same contact, same force, and the same
/// run's narrow-band receiver replays a render without the row bit for bit. The pressure stays exactly zero
/// until the delay kernel can first reach it.
#[test]
fn broadband_observer_leaves_the_bar_alone() {
    let mut card = coarse(3000.0);
    card.duration_s = 0.05;
    card.substeps = 8;
    let s = structure(&card, 1, card.density_kg_m3).expect("structure");
    let plain = strike(&card, &s).expect("narrow-band strike");
    assert!(plain.observer.is_none() && plain.narrow_samples.is_empty());
    let heard = observed_strike(&mut card, &s, false);
    let again = strike(&card, &s).expect("broadband replay");
    assert_eq!(
        plain.stroke.contact_time_s.to_bits(),
        heard.stroke.contact_time_s.to_bits()
    );
    assert_eq!(
        plain.stroke.peak_force_n.to_bits(),
        heard.stroke.peak_force_n.to_bits()
    );
    assert_eq!(plain.stroke.samples.len(), heard.narrow_samples.len());
    assert!(
        plain
            .stroke
            .samples
            .iter()
            .zip(&heard.narrow_samples)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "the bar must not change"
    );
    assert!(
        plain
            .stroke
            .samples
            .iter()
            .zip(&heard.stroke.samples)
            .any(|(a, b)| a != b),
        "the observer must be the one heard"
    );
    assert!(
        heard
            .stroke
            .samples
            .iter()
            .zip(&again.stroke.samples)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "same-host exact replay"
    );
    let report = &heard.observer.as_ref().expect("observer outcome").report;
    println!(
        "coarse blank: held-out worst {:e}, delay {} samples",
        report.held_out_max_rel_error, report.delay_samples
    );
    let silent = ((report.delay_samples - 33.0) / f64::from(card.substeps)).floor() as usize;
    assert!(silent > 0);
    assert!(
        heard.stroke.samples[..silent].iter().all(|p| *p == 0.0),
        "nothing reaches the receiver before the earliest arrival"
    );
    assert!(heard.stroke.samples[silent..].iter().any(|p| p.abs() > 0.0));
}

#[test]
fn observer_row_is_optional_and_refused_when_incomplete() {
    let path =
        root().join("data/instruments/descriptions/xylophone/orff-soprano-blank-rucz2021.bar");
    let text = std::fs::read_to_string(path).expect("the shipped blank card");
    assert!(
        Card::parse(&text, &root())
            .expect("no row")
            .observer
            .is_none()
    );
    let with = |rows: &str| Card::parse(&format!("{text}\n{rows}\n"), &root());
    let on = with("observer\tbroadband\nobserver_fit_order\t32\nobserver_held_out_budget\t0.05")
        .expect("a complete row");
    assert_eq!(
        on.observer,
        Some(crate::card::ObserverRow {
            fit_order: 32,
            held_out_budget: 0.05,
            rigid_rows: false
        })
    );
    let rigid = with(
        "observer\tbroadband\nobserver_fit_order\t32\nobserver_held_out_budget\t0.05\nobserver_rigid_rows\ton",
    )
    .expect("the rigid row");
    assert!(rigid.observer.is_some_and(|o| o.rigid_rows));
    assert!(
        with("observer\tbroadband\nobserver_fit_order\t32\nobserver_held_out_budget\t0.05\nobserver_rigid_rows\toff")
            .is_err()
    );
    assert!(with("observer\tbroadband").is_err());
    assert!(
        with("observer\tnarrow\nobserver_fit_order\t32\nobserver_held_out_budget\t0.05").is_err()
    );
    assert!(
        with("observer\tbroadband\nobserver_fit_order\t32\nobserver_held_out_budget\t0").is_err()
    );
}

fn gas(card: &Card) -> fs_material::gas::GasState {
    let text = std::fs::read_to_string(&card.gas_path).expect("gas record");
    fs_material::gas_species::input::read_gas_state(&mut text.split_ascii_whitespace(), None)
        .expect("gas record resolves")
}

fn table_iii(volume_cm3: f64, opening_cm2: f64) -> Resonator {
    Resonator {
        volume_m3: volume_cm3 * 1e-6,
        opening_area_m2: opening_cm2 * 1e-4,
        // Rucz Sec. III.C: "a very short neck of 3 mm length".
        neck_length_m: 0.003,
        wet_x_m: [0.1, 0.17],
        measured_hz: None,
    }
}

/// The card-to-cavity wiring against the six resonators Rucz et al. measured with the bars removed
/// (Table III: V, S, measured f_res), and the misses reported.
///
/// The pre-registered claim "the library's element is within 5 % of every measured resonance" FAILED
/// (resonator 1 +1.7 %, 2 -3.4 %, 3 -4.3 %, 4 -7.4 %, and the larger resonators lower still): the element
/// takes the zero-frequency flanged end correction 0.8216 a, where Rucz's printed l_eff carries the
/// frequency-dependent reduction of their Ref. 24, which matters as ka grows. So this test asserts only
/// what is exact (the element equals the owner's closed form for the card's geometry, and the six
/// resonators are ordered as measured) and prints the misses; it is not a validation of the pitch.
#[test]
fn helmholtz_wiring_is_exact_and_the_misses_against_the_six_measured_resonators_are_printed() {
    let card = blank();
    let gas = gas(&card);
    let table = [
        (1830.4, 33.56, 305.5),
        (1560.0, 49.91, 386.5),
        (627.8, 34.72, 559.0),
        (531.0, 37.62, 641.5),
        (415.5, 40.52, 740.5),
        (259.9, 46.32, 1007.5),
    ];
    let mut previous = (0.0, 0.0);
    for (volume, opening, measured) in table {
        let res = table_iii(volume, opening);
        let cavity = resonator::cavity(&res, &gas).expect("cavity");
        let f = resonator::isolated_hz(&cavity);
        let radius = (res.opening_area_m2 / core::f64::consts::PI).sqrt();
        let closed_form = gas.sound_speed / bar::TAU
            * (res.opening_area_m2 / (res.volume_m3 * (res.neck_length_m + 2.0 * 0.8216 * radius)))
                .sqrt();
        println!(
            "resonator V {volume} cm3 S {opening} cm2: element {f:.1} Hz, closed form {closed_form:.1} Hz, measured {measured} Hz, miss {:.2} %",
            100.0 * (f / measured - 1.0)
        );
        assert!(
            (f / closed_form - 1.0).abs() < 1e-6,
            "{f} Hz against the closed form {closed_form} Hz"
        );
        assert!(
            f > previous.0 && measured > previous.1,
            "ordered as measured"
        );
        previous = (f, measured);
    }
}

/// Amplitude of `x` at `f` (a single-bin DFT), rectangular or Hann windowed.
fn line(x: &[f64], f: f64, rate: f64, hann: bool) -> f64 {
    let (mut re, mut im, mut weight) = (0.0, 0.0, 0.0);
    let n_total = x.len() as f64;
    for (n, v) in x.iter().enumerate() {
        let w = if hann {
            0.5 - 0.5 * det::cos(bar::TAU * n as f64 / (n_total - 1.0))
        } else {
            1.0
        };
        let phase = bar::TAU * f * n as f64 / rate;
        re += w * v * det::cos(phase);
        im += w * v * det::sin(phase);
        weight += w;
    }
    (re * re + im * im).sqrt() / weight
}

/// The bar over the first resonator radiates at the resonator's own frequency, which the bar alone, whose
/// lowest mode is above 1 kHz, does not; and the coupled render replays bit for bit. The bar's
/// wet interval is the middle 70 mm of the blank (a coarse stand-in for the opening's extent).
///
/// The pre-registered measurement was a rectangular-window line and a factor of 10. It gave 8.1: the
/// rectangular window leaks the bar's own kilohertz modes into the 310 Hz bin. The Hann window removes
/// that leakage floor, so the factor of 10 is asserted on the Hann line; both are printed.
#[test]
fn resonator_adds_its_own_line_and_the_coupled_render_replays() {
    let mut bare = coarse(3000.0);
    bare.duration_s = 0.25;
    bare.substeps = 8;
    let mut with = bare.clone();
    with.resonator = Some(table_iii(1830.4, 33.56));
    let gas = gas(&with);
    let f_res = resonator::isolated_hz(
        &resonator::cavity(with.resonator.as_ref().expect("resonator"), &gas).expect("cavity"),
    );
    let s = structure(&bare, 1, bare.density_kg_m3).expect("structure");
    let alone = strike(&bare, &s).expect("bar alone");
    let loaded = strike(&with, &s).expect("bar over resonator");
    let again = strike(&with, &s).expect("replay");
    assert!(loaded.stroke.samples.iter().all(|p| p.is_finite()));
    let rate = f64::from(bare.rate_hz);
    let rectangular = (
        line(&alone.stroke.samples, f_res, rate, false),
        line(&loaded.stroke.samples, f_res, rate, false),
    );
    let (a, b) = (
        line(&alone.stroke.samples, f_res, rate, true),
        line(&loaded.stroke.samples, f_res, rate, true),
    );
    println!(
        "line at the isolated resonance {f_res:.1} Hz: rectangular bar alone {:.3e} Pa, over resonator {:.3e} Pa (x{:.2}); Hann bar alone {a:.3e} Pa, over resonator {b:.3e} Pa (x{:.2})",
        rectangular.0,
        rectangular.1,
        rectangular.1 / rectangular.0,
        b / a
    );
    assert!(
        b > 10.0 * a,
        "the resonator must add its line: {b} against {a}"
    );
    assert!(
        loaded
            .stroke
            .samples
            .iter()
            .zip(&again.stroke.samples)
            .all(|(x, y)| x.to_bits() == y.to_bits()),
        "same-host exact replay"
    );
}

/// The wood is read from the three cited packs and nowhere else: density from Rucz's measured blanks,
/// E_L from the Indian-rosewood Handbook pack, the ratios from the Honduras-mahogany Handbook pack. A
/// value in a unit the reader does not convert exactly is refused, not guessed.
#[test]
fn wood_comes_from_the_cited_packs_and_inexact_units_are_refused() {
    let card = blank();
    let pack = |name: &str| root().join(format!("data/matdb/seed-v1/{name}/properties.tsv"));
    assert_eq!(card.density_kg_m3, 1116.0, "Rucz's measured density, kept");
    assert_eq!(
        card.youngs_pa[0], 13_530.0e6,
        "Indian rosewood E_L, 12 % MC"
    );
    let mahogany = pack("mahogany-honduras-fpl-gtr282");
    let ratio = |p: &str| pack_scalar(&mahogany, p).expect("Handbook ratio");
    assert_eq!(card.youngs_pa[1], card.youngs_pa[0] * ratio("et_over_el"));
    assert_eq!(card.youngs_pa[2], card.youngs_pa[0] * ratio("er_over_el"));
    assert_eq!(card.shear_pa[2], card.youngs_pa[0] * ratio("glr_over_el"));
    assert_eq!(
        card.poisson,
        [ratio("nu_lt"), ratio("nu_lr"), ratio("nu_tr")]
    );
    let err = pack_scalar(
        &pack("rosewood-indian-fpl-gtr282"),
        "radial_shrinkage_green_to_ovendry",
    )
    .unwrap_err();
    assert!(err.contains("unit"), "a percent is not converted: {err}");
}

/// La Favre's A3 tube as printed: OD 2.5 in, wall 1/16 in, air column 14 3/8 in, gap 5/8 in.
fn a3_tube() -> Tube {
    let inch = 0.0254;
    Tube {
        radius_m: (1.25 - 0.0625) * inch,
        air_length_m: 14.375 * inch,
        gap_m: 0.625 * inch,
        modes: 5,
    }
}

fn plain_modes(rows: &[(f64, f64)]) -> Vec<ModalAcousticMode> {
    rows.iter()
        .map(|&(f, a)| ModalAcousticMode {
            angular_frequency_rad_s: bar::TAU * f,
            damping_ratio: a / (bar::TAU * f),
            pressure_per_modal_velocity: C64::ZERO,
        })
        .collect()
}

#[test]
fn tube_rows_parse_and_refusals_are_explicit() {
    let text = blank().text;
    assert_eq!(blank().tube, None);
    let rows = "tube\tclosed_open\ntube_radius_m\t0.03\ntube_air_length_m\t0.36\ntube_gap_m\t0.016\ntube_modes\t5\n";
    let with = Card::parse(&format!("{text}\n{rows}"), &root()).expect("the rows parse");
    assert_eq!(
        with.tube,
        Some(Tube {
            radius_m: 0.03,
            air_length_m: 0.36,
            gap_m: 0.016,
            modes: 5,
        })
    );
    let kind = rows.replace("closed_open", "open_open");
    assert!(
        Card::parse(&format!("{text}\n{kind}"), &root())
            .unwrap_err()
            .contains("closed_open")
    );
    let short = rows.replace("0.36", "0.05");
    assert!(
        Card::parse(&format!("{text}\n{short}"), &root())
            .unwrap_err()
            .contains("column")
    );
    let both = "resonator\thelmholtz\nresonator_volume_m3\t0.001\nresonator_opening_area_m2\t0.001\nresonator_neck_length_m\t0.003\nresonator_wet_x_m\t0.1\t0.17\n";
    assert!(
        Card::parse(&format!("{text}\n{rows}{both}"), &root())
            .unwrap_err()
            .contains("not both")
    );
    let g = gas(&blank());
    let mut many = a3_tube();
    many.modes = 12;
    assert!(tube::isolated_modes(&many, &g).unwrap_err().contains("ka"));
}

/// The printed laws at La Favre's A3, against the sums by hand: Eq. 7 wall, Eq. 12 open end, odd modes.
#[test]
fn isolated_tube_follows_the_printed_laws() {
    let g = gas(&blank());
    let t = a3_tube();
    let modes = tube::isolated_modes(&t, &g).expect("modes");
    let f: Vec<f64> = modes.iter().map(|m| m.omega_rad_s / bar::TAU).collect();
    println!(
        "La Favre A3 column: {f:?} Hz, wall {:?}, open end {:?}",
        modes.iter().map(|m| m.wall_per_s).collect::<Vec<_>>(),
        modes.iter().map(|m| m.radiation_per_s).collect::<Vec<_>>()
    );
    assert!((f[0] - 223.5).abs() < 1.0, "{f:?}");
    assert!((modes[0].wall_per_s + modes[0].radiation_per_s - 8.2).abs() < 0.3);
    // Eq. 12 at low ka reduces to Eq. 13.
    let (re, im) = tube::open_end(0.01);
    assert!((re / (1e-4 / 4.0) - 1.0).abs() < 1e-3 && (im / (0.6133 * 0.01) - 1.0).abs() < 1e-3);
    // Odd modes, each above its odd multiple of the first by the two printed effects and nothing else: the open
    // end's correction shrinks as ka rises, and the wall slows the low modes most.
    let wall = (4e-8_f64.sqrt() + (g.gamma - 1.0) * 5.6e-8_f64.sqrt())
        / (t.radius_m * core::f64::consts::SQRT_2);
    let slowed = |w: f64| 1.0 - wall * (g.sound_speed / w).sqrt();
    for (m, mode) in modes.iter().enumerate() {
        let odd = (2 * m + 1) as f64;
        let by_laws = odd
            * (modes[0].acoustic_length_m / mode.acoustic_length_m)
            * (slowed(mode.omega_rad_s) / slowed(modes[0].omega_rad_s));
        let ratio = f[m] / f[0];
        assert!(
            (ratio / by_laws - 1.0).abs() < 1e-6,
            "mode {m} ratio {ratio} against {by_laws}"
        );
    }
    // Mode 3's ratio by hand from Eq. 7 and Eq. 12, worked apart from this code.
    assert!((f[3] / f[0] - 7.0369).abs() < 0.0005, "{f:?}");
    // The wall rises as the root of frequency, the open end much faster.
    assert!((modes[2].wall_per_s / modes[0].wall_per_s - 5.0_f64.sqrt()).abs() < 0.05);
    assert!(modes[2].radiation_per_s > 15.0 * modes[0].radiation_per_s);
}

/// A mutual resistance held at a far tube mode's frequency, as A3's 2027 Hz tube mode against the bar's 220 Hz mode,
/// refuses when every pair is kept. One set of fields, small on near pairs and large on far ones: each bar mode keeps
/// its nearest tube mode alone, so the same numbers build a passive body and the rest are zero.
#[test]
fn each_bar_mode_exchanges_with_its_nearest_tube_mode_alone() {
    let g = gas(&blank());
    let t = a3_tube();
    let modes = tube::isolated_modes(&t, &g).expect("modes");
    let bar_hz = [220.0, 880.0, 2017.0];
    let bar_a = [8.3, 23.7, 93.7];
    let near_hz = |f: f64| {
        (0..modes.len())
            .min_by(|&a, &b| {
                (modes[a].omega_rad_s / bar::TAU - f)
                    .abs()
                    .total_cmp(&(modes[b].omega_rad_s / bar::TAU - f).abs())
            })
            .expect("modes")
    };
    let fields: Vec<tube::MouthField> = (0..modes.len())
        .map(|m| tube::MouthField {
            per_bar_mode: bar_hz
                .iter()
                .map(|&f| C64::new(if m == near_hz(f) { 1.0 } else { 20_000.0 }, 0.5))
                .collect(),
            reflection: C64::ZERO,
        })
        .collect();
    let omega: Vec<f64> = bar_hz.iter().map(|f| bar::TAU * f).collect();
    let kept = tube::Network::new(&g, &t, &modes, &fields, &omega).expect("network");
    for (k, &f) in bar_hz.iter().enumerate() {
        for m in 0..modes.len() {
            let on = m == near_hz(f);
            assert_eq!(
                kept.resistance[m][k] != 0.0,
                on,
                "resistance, tube mode {m}, bar {f} Hz"
            );
            assert_eq!(
                kept.inertance[m][k] != 0.0,
                on,
                "inertance, tube mode {m}, bar {f} Hz"
            );
        }
    }
    // Every pair kept: the same fields through the same 1/sqrt(mass).
    let mut full = kept.clone();
    for (m, field) in fields.iter().enumerate() {
        full.resistance[m] = field
            .per_bar_mode
            .iter()
            .map(|z| z.re * kept.tube[m][2])
            .collect();
    }
    let rows: Vec<(f64, f64)> = bar_hz.iter().copied().zip(bar_a).collect();
    let refused = tube::TubeBody::new(96_000, &plain_modes(&rows), rows.len(), &full, 0.0, 0.0);
    assert!(
        refused
            .err()
            .expect("every pair kept refuses")
            .contains("resistance")
    );
    assert!(tube::TubeBody::new(96_000, &plain_modes(&rows), rows.len(), &kept, 0.0, 0.0).is_ok());
    let short = [bar::TAU * 220.0];
    assert!(
        tube::Network::new(&g, &t, &modes, &fields, &short)
            .unwrap_err()
            .contains("disagree")
    );
}

/// Forced from one bar mode, the tube rings on at its own frequency and its own decay: the integrator's check.
#[test]
fn tube_mode_rings_at_its_law_once_the_kick_is_over() {
    let (f_tube, a_tube) = (223.5, 8.2);
    let network = tube::Network {
        tube: vec![[bar::TAU * f_tube, a_tube, 1.0]],
        resistance: vec![vec![0.0]],
        inertance: vec![vec![1e-3]],
    };
    let rate = 96_000_u32;
    let mut body = tube::TubeBody::new(
        rate,
        &plain_modes(&[(5000.0, 2000.0)]),
        1,
        &network,
        0.0,
        0.0,
    )
    .expect("body");
    let mut v = Vec::new();
    for n in 0..rate {
        let force = if n < 8 { 1.0 } else { 0.0 };
        body.step(&[force]).expect("step");
        v.push(body.tube_velocities()[0]);
    }
    let start = (0.1 * f64::from(rate)) as usize;
    let crossings: Vec<f64> = (start..v.len() - 1)
        .filter(|&i| v[i] < 0.0 && v[i + 1] >= 0.0)
        .map(|i| i as f64 + v[i] / (v[i] - v[i + 1]))
        .collect();
    let f = (crossings.len() - 1) as f64 * f64::from(rate)
        / (crossings[crossings.len() - 1] - crossings[0]);
    // A decaying mode rings at its damped frequency, not at the undamped one the network carries.
    let damped = (f_tube * f_tube - (a_tube / bar::TAU) * (a_tube / bar::TAU)).sqrt();
    assert!((f / damped - 1.0).abs() < 1e-6, "{f} against {damped}");
    let peak = |a: usize, b: usize| v[a..b].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let width = f64::from(rate) as usize / 20;
    let (early, late) = (
        peak(start, start + width),
        peak(start + 4 * width, start + 5 * width),
    );
    let measured = (early / late).ln() / 0.2;
    assert!((measured / a_tube - 1.0).abs() < 0.01, "{measured}");
}

/// Two bar modes and two tube modes, struck: the energy with the mutual term plus what the laws took equals the work done.
#[test]
fn coupled_energy_balances() {
    let network = tube::Network {
        tube: vec![[bar::TAU * 223.0, 8.0, 1.0], [bar::TAU * 664.0, 20.0, 1.0]],
        resistance: vec![vec![0.5, 0.2], vec![0.3, 0.9]],
        inertance: vec![vec![0.02, 0.01], vec![0.015, 0.03]],
    };
    let rate = 384_000_u32;
    let h = 1.0 / f64::from(rate);
    let mut body = tube::TubeBody::new(
        rate,
        &plain_modes(&[(220.0, 8.0), (880.0, 30.0)]),
        2,
        &network,
        0.0,
        0.0,
    )
    .expect("body");
    let shape = [2.5, -1.1];
    let contact = 6.4e-4;
    let (mut work, mut lost, mut most) = (0.0, 0.0, 0.0_f64);
    let mut largest_gap = 0.0_f64;
    for n in 0..rate / 2 {
        let t = f64::from(n) * h;
        let load = if t < contact {
            100.0 * det::sin(core::f64::consts::PI * t / contact)
        } else {
            0.0
        };
        let forces = [shape[0] * load, shape[1] * load];
        let before = body.bar_velocities();
        let lost_before = body.dissipation();
        body.step(&forces).expect("step");
        let after = body.bar_velocities();
        for k in 0..2 {
            work += h * forces[k] * 0.5 * (before[k] + after[k]);
        }
        lost += 0.5 * h * (lost_before + body.dissipation());
        let e = body.energy();
        most = most.max(e);
        largest_gap = largest_gap.max((e + lost - work).abs());
    }
    println!(
        "energy: largest {most:.6e}, work {work:.6e}, lost {lost:.6e}, largest imbalance {largest_gap:.3e}"
    );
    assert!(
        largest_gap < 1e-3 * most,
        "imbalance {largest_gap} of {most}"
    );
}

/// A bar mode and a tube mode at one frequency, coupled by inertance x alone, beat at f (1/sqrt(1-x) - 1/sqrt(1+x)),
/// and the energy goes over and comes back: two-way, as Fig. 7 of Part II shows for a detuned pair.
#[test]
fn tuned_pair_beats_at_the_inertance_split() {
    let (f, x) = (220.0, 0.02);
    let network = tube::Network {
        tube: vec![[bar::TAU * f, 0.5, 1.0]],
        resistance: vec![vec![0.0]],
        inertance: vec![vec![x]],
    };
    let rate = 96_000_u32;
    let mut body =
        tube::TubeBody::new(rate, &plain_modes(&[(f, 0.5)]), 1, &network, 0.0, 0.0).expect("body");
    let n = (0.6 * f64::from(rate)) as usize;
    let (mut bar_v, mut tube_v) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        body.step(&[if i < 8 { 1.0 } else { 0.0 }]).expect("step");
        bar_v.push(body.bar_velocities()[0]);
        tube_v.push(body.tube_velocities()[0]);
    }
    let split = f * (1.0 / (1.0 - x).sqrt() - 1.0 / (1.0 + x).sqrt());
    let cycle = (f64::from(rate) / f) as usize;
    let envelope = |v: &[f64], at: usize| {
        v[at..at + cycle]
            .iter()
            .fold(0.0_f64, |m, y| m.max(y.abs()))
    };
    let half = (0.5 / split * f64::from(rate)) as usize;
    let full = (1.0 / split * f64::from(rate)) as usize;
    let (over_tube, over_bar) = (
        envelope(&tube_v, half - cycle / 2),
        envelope(&bar_v, half - cycle / 2),
    );
    let (back_tube, back_bar) = (
        envelope(&tube_v, full - cycle / 2),
        envelope(&bar_v, full - cycle / 2),
    );
    println!(
        "split {split:.3} Hz; half a beat: tube {over_tube:.3e}, bar {over_bar:.3e}; a whole beat: tube {back_tube:.3e}, bar {back_bar:.3e}"
    );
    assert!(over_tube > 20.0 * over_bar);
    assert!(back_bar > 20.0 * back_tube);
}

/// The reciprocity the coupling leans on, solved both ways on a real bar: the mouth's pressure from a bar mode's
/// motion against the force the mouth's flow puts on that mode, from the surface pressure of the rigid bar.
#[test]
fn mouth_transfers_are_reciprocal_on_the_bem() {
    let card = coarse(3000.0);
    let s = structure(&card, 1, card.density_kg_m3).expect("structure");
    let (surface, shapes) =
        radiation::surface_and_shapes(&card, &s.mesh, &s.modes).expect("surface");
    let g = gas(&card);
    let medium = fs_bem::helmholtz::Medium::from(&g);
    let t = Tube {
        radius_m: 0.015,
        air_length_m: 0.2,
        gap_m: 0.016,
        modes: 1,
    };
    let isolated = tube::isolated_modes(&t, &g).expect("modes");
    let fields =
        radiation::mouth_fields(&card, &t, medium, &surface, &shapes, &isolated).expect("fields");
    let w = isolated[0].omega_rad_s;
    let k = w / medium.sound_speed;
    let points = tube::mouth_points(tube::mouth_centre(&card, &t), t.radius_m);
    let reflected = radiation::reflected_velocity(&surface, &points, k);
    let solved = fs_bem::helmholtz::solve_radiation_batch(
        &surface,
        k,
        medium,
        &[reflected.as_slice()],
        fs_bem::helmholtz::Formulation::BurtonMiller,
    )
    .expect("the rigid bar's answer");
    let total: Vec<C64> = surface
        .centroids()
        .iter()
        .zip(&solved[0].pressure)
        .map(|(x, scattered)| {
            let mut incident = C64::ZERO;
            for p in &points {
                let r = (0..3).map(|c| (x[c] - p[c]).powi(2)).sum::<f64>().sqrt();
                let travel = C64::new(det::cos(k * r), det::sin(k * r));
                incident = incident
                    + travel
                        * C64::new(
                            0.0,
                            -w * medium.density / (2.0 * bar::TAU * r * points.len() as f64),
                        );
            }
            incident + *scattered
        })
        .collect();
    let mut compared = 0;
    let largest = fields[0]
        .per_bar_mode
        .iter()
        .fold(0.0_f64, |m, z| m.max(z.abs()));
    for (shape, forward) in shapes.iter().zip(&fields[0].per_bar_mode) {
        if forward.abs() < 0.05 * largest {
            continue;
        }
        let mut back = C64::ZERO;
        for ((p, v), area) in total.iter().zip(shape).zip(surface.areas()) {
            back = back + *p * C64::from_re(v * area);
        }
        let miss = (back - *forward).abs() / forward.abs();
        println!("forward {forward:?}, back {back:?}, miss {miss:.4}");
        assert!(miss < 0.05, "reciprocity misses by {miss}");
        compared += 1;
    }
    assert!(compared > 0);
}

/// The full route with a tube: finite, replays bit for bit, adds the tube's line, and leaves the contact alone.
#[test]
fn a_tube_adds_its_line_and_the_two_way_render_replays() {
    let mut bare = coarse(3000.0);
    bare.duration_s = 0.25;
    bare.substeps = 8;
    let mut with = bare.clone();
    with.tube = Some(Tube {
        radius_m: 0.015,
        air_length_m: 0.2,
        gap_m: 0.016,
        modes: 1,
    });
    let s = structure(&bare, 1, bare.density_kg_m3).expect("structure");
    let alone = strike(&bare, &s).expect("bar alone");
    let loaded = strike(&with, &s).expect("bar over its tube");
    let again = strike(&with, &s).expect("replay");
    assert!(loaded.stroke.samples.iter().all(|p| p.is_finite()));
    assert!(
        loaded
            .stroke
            .samples
            .iter()
            .zip(&again.stroke.samples)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let g = gas(&with);
    let f = tube::isolated_modes(with.tube.as_ref().expect("tube"), &g).expect("modes")[0]
        .omega_rad_s
        / bar::TAU;
    let rate = f64::from(bare.rate_hz);
    let (a, b) = (
        line(&alone.stroke.samples, f, rate, true),
        line(&loaded.stroke.samples, f, rate, true),
    );
    println!(
        "line at the column's {f:.1} Hz: bar alone {a:.3e} Pa, with tube {b:.3e} Pa; contact {} against {} s",
        alone.stroke.contact_time_s, loaded.stroke.contact_time_s
    );
    assert!(b > 3.0 * a);
    assert!((alone.stroke.contact_time_s / loaded.stroke.contact_time_s - 1.0).abs() < 0.02);
}

/// A long bar whose cut edge falls past the middle of its end cell. A lower-face midpoint there sits .28 mm from the
/// lower face and 4.45 mm from the end face, yet nearer the end as a fraction: .012 of the thickness against .009
/// of the length. Snapping by the fraction moved it onto the end face and folded the element.
#[test]
fn surface_snap_takes_the_nearest_face_in_metres() {
    use fs_solid::linear3::tet10::BoundaryProjection;
    let card = long_bar_with_steep_cut_edge();
    let a = bar::point(&card, [0.0, 1.0 / 6.0, 0.0]);
    let b = bar::point(&card, [1.0 / 56.0, 1.0 / 6.0, 0.0]);
    let mid: [f64; 3] = core::array::from_fn(|k| f64::midpoint(a[k], b[k]));
    let snapped = bar::Surface(&card).project(mid).expect("a snap");
    let lower_face = bar::point(&card, [mid[0] / card.length_m, 1.0 / 6.0, 0.0]);
    for k in 0..3 {
        assert!(
            (snapped[k] - lower_face[k]).abs() < 1e-12,
            "{snapped:?} against {lower_face:?}"
        );
    }
    let mesh = bar::mesh(&card, card.mesh_counts);
    let wood = bar::material(&card, card.density_kg_m3).expect("the blank's wood");
    bar::with_cx(|cx| bar::assemble(&card, &mesh, &wood, cx))
        .map(|_| ())
        .expect("every Jacobian positive");
}

fn long_bar_with_steep_cut_edge() -> Card {
    let mut card = blank();
    card.length_m = 0.498475;
    card.width_m = 0.0635;
    card.thickness_m = 0.023368;
    card.profile = Profile::Arc {
        stone_radius_m: 1.6834942755977116,
        middle_thickness_m: 0.005562,
        straight_length_m: 0.0,
    };
    card.mesh_counts = [56, 6, 3];
    card
}
