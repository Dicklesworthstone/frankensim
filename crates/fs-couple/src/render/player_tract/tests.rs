use super::*;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::schedule::reed::ReedPerformance;
use crate::upstream_airway::AirwaySection;
use fs_material::gas::GasSpec;

/// Kaburagi, Kato, Fukuda and Taguchi 2021, Figs. 4 and 5 as read into the data file (its header names the reading).
const AREAS: &str = include_str!(
    "../../../../../data/instruments/descriptions/players/kaburagi2021-saxophonists-tract-area-functions.tsv"
);
const VISCOUS_INPUT: &str = include_str!("../../../examples/reed-viscous-slit.performance");
/// The held pressure of [`card`] [Pa].
const HELD_PA: f64 = 2800.0;
/// Samples before the hold window and in it (0.2 s each at 48 kHz).
const HALF: usize = 9600;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

/// Participant A's long-tube low-note tract (the posture the E3 and F3 pages take), valve side first: the data file's
/// `A_long_low` column reversed, each section Table 1's 16.0 cm / 40 long, then the glottis of Sec. 4.1 (0.22 cm^2,
/// 0.25 cm).
fn long_low() -> Vec<(f64, f64)> {
    let mut rows = AREAS.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<&str> = rows.next().unwrap().split('\t').collect();
    let column = header
        .iter()
        .position(|name| *name == "A_long_low")
        .unwrap();
    let mut areas: Vec<f64> = rows
        .map(|line| {
            line.split('\t')
                .nth(column)
                .unwrap()
                .parse::<f64>()
                .unwrap()
                * 1.0e-4
        })
        .collect();
    assert_eq!(areas.len(), 40);
    areas.reverse();
    let mut sections: Vec<(f64, f64)> =
        areas.into_iter().map(|area| (0.160 / 40.0, area)).collect();
    sections.push((0.0025, 0.22e-4));
    sections
}

fn chain(sections: &[(f64, f64)], generations: usize) -> AirwayChain {
    let tubes = sections
        .iter()
        .map(|&(l, a)| AirwaySection::from_area(l, a, 1, None).unwrap())
        .collect();
    AirwayChain::new(tubes, fs_duct::LossModel::Bessel)
        .unwrap()
        .ahead_of(&AirwayChain::weibel_saab_2021_b(generations, None).unwrap())
        .unwrap()
}

fn tract_rows(pressure: &str, sections: &[(f64, f64)], generations: usize) -> String {
    let mut rows = format!("player_tract {pressure}\n");
    for (length_m, area_m2) in sections {
        rows.push_str(&format!("player_tract_section {length_m:e} {area_m2:e}\n"));
    }
    rows.push_str(&format!(
        "player_airway weibel-1963 saab-2021-table-1-B {generations}\n"
    ));
    rows
}

/// The slit fixture as a held note: 0.4 s, `HELD_PA` reached in 10 ms and held to the end, with `tract` before
/// `compile_limits`.
fn card(tract: Option<&str>) -> String {
    let mut text = VISCOUS_INPUT.replace("audio 48000 4801 20000", "audio 48000 19200 20000");
    let events = text.find("events\t5").unwrap();
    text.truncate(events);
    text.push_str("events\t1\nevent\t0e0\t1e-2\tpressure-pa\t2.8e3\n");
    match tract {
        Some(rows) => text.replace("compile_limits", &format!("{rows}compile_limits")),
        None => text,
    }
}

/// Render `text`; return the Hann-weighted means of the mouth pressure and of the flow into the tract over 0.2-0.4 s,
/// the plain mean of the mouth pressure there, the output, and the tract's `R_dc` and margin.
fn hold(text: &str) -> (f64, f64, f64, Vec<f64>, f64, Option<(f64, f64)>) {
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), HALF).unwrap();
    let mut out = vec![0.0; 2 * HALF];
    performance.block(&mut out[..HALF]).unwrap();
    let (mut weight, mut mouth, mut flow, mut plain) = (0.0, 0.0, 0.0, 0.0);
    for n in 0..HALF {
        performance.block(&mut out[HALF + n..HALF + n + 1]).unwrap();
        let tract = performance.player_tract().unwrap();
        let w = 0.5 - 0.5 * (core::f64::consts::TAU * (n as f64 + 0.5) / HALF as f64).cos();
        weight += w;
        mouth += w * tract.mouth_pressure_pa();
        flow += w * tract.flow_m3_s();
        plain += tract.mouth_pressure_pa();
    }
    let tract = performance.player_tract().unwrap();
    let (r_dc, margin) = (tract.dc_resistance(), tract.passivity_margin());
    assert!(out.iter().all(|p| p.is_finite()));
    (
        mouth / weight,
        flow / weight,
        plain / HALF as f64,
        out,
        r_dc,
        margin,
    )
}

// The fold is the reed law with the tract in series: solved on the folded boundary and unfolded, the bore's wave
// satisfies the original massless junction at the drop p_oc - Z (U - u_body) - p_bore, and the massive path's load
// coordinate is the bore's plus Z with the same incoming wave.
#[test]
fn folded_junction_is_the_reed_law_with_the_tract_in_series() {
    let gas = air();
    let airway =
        UpstreamAirway::realize(&chain(&long_low(), 3), &gas, 48_000, &LineLadder::default())
            .unwrap();
    let reed = BeatingReed {
        rest_opening_m: 0.0004,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: HELD_PA,
        attack_s: 0.0,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    };
    let zc = 2.0e6;
    for &(d, h, u_body, series) in &[
        (0.0, 0.0, 0.0, 0.0),
        (0.3, -150.0, 2.0e-6, 0.0),
        (-0.4, 420.0, -5.0e-6, 3.0),
        (0.9, 900.0, 0.0, -2.0),
    ] {
        let mut tract = PlayerTract {
            airway: airway.clone(),
            pressure: TractPressure::Lung,
            dc_resistance: 0.0,
            pole: 0.0,
            flow_lp: 0.0,
            margin: None,
            mouth_pa: 0.0,
            flow_m3_s: 0.0,
            gas,
            ledger: AirwayLedger::default(),
        };
        let boundary = ScatteringResponse {
            direct: d,
            history: h,
        };
        let (folded, p_junction, step) = tract.fold(HELD_PA, series, boundary, zc, u_body).unwrap();
        let q = crate::reed_bore::solve_reed_wave(
            reed,
            gas.density,
            zc,
            folded.direct,
            folded.history,
            p_junction,
            0.0,
            u_body,
        )
        .unwrap();
        let a = step.unfold(q);
        let p_minus = d * a + h;
        let p_bore = a + p_minus;
        let u = (a - p_minus) / zc;
        let mouth = tract.commit(step, u - u_body).unwrap();
        assert!(
            (mouth - (step.open_circuit_pa - step.impedance * (u - u_body))).abs()
                <= 1e-9 * HELD_PA
        );
        let jet = crate::reed_bore::aperture_of(reed).volume_flow(mouth - p_bore, gas.density);
        let scale = jet.abs().max(u.abs()).max(1e-9);
        assert!(
            (jet + u_body - u).abs() <= 1e-8 * scale,
            "d {d} h {h}: jet {jet} body {u_body} wave {u}"
        );
        let load = zc * (1.0 + folded.direct) / (1.0 - folded.direct);
        let expected = zc * (1.0 + d) / (1.0 - d) + step.impedance;
        assert!(
            (load - expected).abs() <= 1e-12 * expected,
            "load {load} against {expected}"
        );
        let incoming = folded.history / (1.0 - folded.direct);
        assert!((incoming - h / (1.0 - d)).abs() <= 1e-12 * (1.0 + (h / (1.0 - d)).abs()));
    }
}

// MB 18593: with a mouth pressure the hold mean is the card's within 1e-3; with a lung pressure it is the card's less
// R_dc times the mean flow (the airway's DC law), so the planted test sees the difference the subtraction makes.
#[test]
fn mouth_hold_mean_is_the_card_pressure_and_lung_hold_mean_is_the_dc_law() {
    let sections = long_low();
    let (mouth, flow, plain, _, r_dc, margin) =
        hold(&card(Some(&tract_rows("mouth 0.02", &sections, 17))));
    eprintln!(
        "mouth 20 ms: Hann mean {mouth} Pa, plain mean {plain} Pa, flow {flow} m3/s, R_dc {r_dc}, margin {margin:?}"
    );
    assert!(
        (mouth - HELD_PA).abs() <= 1e-3 * HELD_PA,
        "mouth hold mean {mouth}"
    );
    assert!(margin.is_some_and(|(m, _)| m >= 0.0));
    let (lung, lung_flow, lung_plain, _, lung_r_dc, lung_margin) =
        hold(&card(Some(&tract_rows("lung", &sections, 17))));
    eprintln!(
        "lung: Hann mean {lung} Pa, plain mean {lung_plain} Pa, flow {lung_flow} m3/s, R_dc*U {}",
        lung_r_dc * lung_flow
    );
    assert!(lung_margin.is_none() && lung_r_dc == r_dc);
    assert!(
        (lung - (HELD_PA - lung_r_dc * lung_flow)).abs() <= 1e-3 * HELD_PA,
        "lung hold mean {lung}"
    );
}

// MB 18593: tau is a numerical regulation constant; across 5 to 50 ms the hold mean moves by less than 1e-3.
#[test]
fn hold_mean_moves_under_1e_3_between_tau_5_and_50_ms() {
    let sections = long_low();
    let (short, ..) = hold(&card(Some(&tract_rows("mouth 0.005", &sections, 17))));
    let (long, ..) = hold(&card(Some(&tract_rows("mouth 0.05", &sections, 17))));
    eprintln!("tau 5 ms {short} Pa, tau 50 ms {long} Pa");
    assert!(
        (short - long).abs() <= 1e-3 * HELD_PA,
        "{short} against {long}"
    );
}

// MB 18593's planted case: a flat subtraction, H = 1, on participant A's long-low tract is an active load and refuses by
// name; the one-pole of 20 ms on the same bins is admitted.
#[test]
fn flat_subtraction_refuses_on_the_long_low_tract_and_the_one_pole_is_admitted() {
    let gas = air();
    let tract = chain(&long_low(), 17);
    let airway = UpstreamAirway::realize(&tract, &gas, 48_000, &LineLadder::default()).unwrap();
    let n = airway.receipt().reflectance.n_fft;
    let flat = admitted_margin(&tract, &gas, 48_000, n, |_| C64::ONE).unwrap_err();
    eprintln!("flat: {flat}");
    assert!(flat.to_string().contains("player tract refused"), "{flat}");
    let dt = 1.0 / 48_000.0;
    let pole = 1.0 - fs_math::det::exp(-dt / 0.02);
    let (margin, hz) = admitted_margin(&tract, &gas, 48_000, n, |omega| {
        one_pole_response(pole, omega, dt)
    })
    .unwrap();
    eprintln!("one-pole 20 ms on {n} bins: margin {margin} Pa s/m^3 at {hz} Hz");
    assert!(margin >= 0.0);
}

#[test]
fn malformed_tract_rows_refuse_before_rendering() {
    let section = "player_tract_section 4e-3 2e-4\n";
    let airway = "player_airway weibel-1963 saab-2021-table-1-B 3\n";
    for rows in [
        format!("player_tract\n{section}"),
        format!("player_tract mouth\n{section}"),
        format!("player_tract mouth 0\n{section}"),
        format!("player_tract mouth -0.02\n{section}"),
        format!("player_tract nose 0.02\n{section}"),
        format!("player_tract lung extra\n{section}"),
        "player_tract lung\n".to_string(),
        format!("{section}player_tract lung\n"),
        "player_tract lung\nplayer_tract_section 0 2e-4\n".to_string(),
        "player_tract lung\nplayer_tract_section 4e-3\n".to_string(),
        format!("player_tract lung\n{section}player_tract_wall 15 -1 0\n"),
        format!("player_tract lung\n{section}player_airway weibel-1963 saab-2021-table-1-B 0\n"),
        format!("player_tract lung\n{section}player_airway weibel-1963 saab-2021-table-1-B 25\n"),
        format!("player_tract lung\n{section}player_airway weibel-1963 saab-2021-table-1-A 3\n"),
        format!("player_tract lung\n{section}{airway}{airway}"),
        format!("player_tract lung\n{section}player_tract lung\n{section}"),
    ] {
        assert!(
            ReedPerformance::from_bytes(card(Some(&rows)).as_bytes(), HALF).is_err(),
            "{rows}"
        );
    }
    let good = format!("player_tract lung\n{section}player_tract_wall 15 1.48e6 1.6e4\n{airway}");
    assert!(
        ReedPerformance::from_bytes(card(Some(&good)).as_bytes(), HALF)
            .unwrap()
            .player_tract()
            .is_some()
    );
    assert!(
        ReedPerformance::from_bytes(card(None).as_bytes(), HALF)
            .unwrap()
            .player_tract()
            .is_none()
    );
}

#[test]
fn the_tract_changes_the_sound() {
    let mut plain = ReedPerformance::from_bytes(card(None).as_bytes(), HALF).unwrap();
    let mut without = vec![0.0; 2 * HALF];
    for part in without.chunks_mut(HALF) {
        plain.block(part).unwrap();
    }
    let (.., with, _, _) = hold(&card(Some(&tract_rows("mouth 0.02", &long_low(), 17))));
    let energy: f64 = without.iter().map(|p| p * p).sum();
    let changed: f64 = without
        .iter()
        .zip(&with)
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    eprintln!("relative change {}", changed / energy);
    assert!(energy > 0.0 && changed > 1e-8 * energy);
}
