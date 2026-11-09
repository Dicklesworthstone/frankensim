//! The `upstream_airway` row on posted flue cards. Without the row three cards, each cut to a 2048-sample window with
//! its first gesture event only, render the pressure bits recorded on the code before the row existed (FNV-1a 64 of
//! the little-endian f64 bits, printed by a probe on 2b00d383c before this change): recorder c2 held (no windway law, the
//! Auvray drive, breath noise), recorder f1 (the fully developed slit) and the flute D4 forced Bickley jet with breath
//! noise. With the row, the two recorder windows take Story's (2008) /u/ tract, lips first, ahead of all 24 printed
//! generations of Weibel's tree with the voice cards' tract wall, and render other bits from a port that is never the
//! ideal one. Without the row's `two-way` token they keep the one-way windows recorded on the one-way tree; with it the
//! windway draws air back in the attack (MB 18587), a declared stand-in coefficient.
use fs_couple::flue_loop::airway::{
    FlueAirwayClosures, FlueAirwayFrame, closed_end_jump_pa, sampled_step_undershoot, source_step_height,
};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const RECORDER_C2_HELD: &str = include_str!("../../../data/instruments/descriptions/recorder/recorder-c2-held.performance");
const RECORDER_F1: &str = include_str!("../../../data/instruments/descriptions/recorder/recorder-f1.performance");
const FLUTE_D4_FORCED: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4-forced-jet-momentum-jet-noise-re700.performance");
const STORY_2008: &str = include_str!("../../../data/instruments/descriptions/voice/story2008-table-I-area-functions.tsv");
const WINDOW: usize = 2048;
/// Digests of the cut cards without the row, recorded at 2b00d383c before this change.
const PLAIN: [(&str, u64); 3] = [("recorder-c2-held", 0xb097_b1f1_a1ca_e2e0), ("recorder-f1", 0xb9b5_3fb0_9bbd_bae2), ("flute-d4-forced", 0x6863_ee01_5d19_cf12)];
/// All the generations Saab et al. (2021) print; the tree's DC resistance changes by under 1e-3 after the 20th.
const GENERATIONS: usize = 24;
/// The voice cards' `tract_wall` (surface density, stiffness per area, resistance).
const WALL: &str = "airway_wall 15 1.48e+06 1.6e4";

fn gate() -> CancelGate { CancelGate::new_clock_free() }

/// The card with `WINDOW` output samples at full scale 1 and its first gesture event only.
fn cut(text: &str) -> String {
    let start = text.find("\noutput ").expect("an output row") + 1;
    let end = start + text[start..].find('\n').expect("output row end");
    let fields: Vec<&str> = text[start..end].split_ascii_whitespace().collect();
    assert_eq!(fields.len(), 4, "output row: {fields:?}");
    let text = format!("{}output {WINDOW} 1 {}{}", &text[..start], fields[3], &text[end..]);
    let at = text.find("\nevents\t").expect("an events row");
    let (head, tail) = text.split_at(at);
    let lines: Vec<&str> = tail[1..].split('\n').collect();
    assert!(lines[1].starts_with("event\t"), "{}", lines[1]);
    let newline = if tail.ends_with('\n') { "\n" } else { "" };
    format!("{head}\nevents\t1\n{}{newline}", lines[1])
}

fn digest(pressure: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for p in pressure {
        for byte in p.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    }
    h
}

fn card(name: &str) -> String {
    match name {
        "recorder-c2-held" => cut(RECORDER_C2_HELD),
        "recorder-f1" => cut(RECORDER_F1),
        "flute-d4-forced" => cut(FLUTE_D4_FORCED),
        _ => unreachable!("{name}"),
    }
}

/// Story's /u/ (the table's last column), lips first: `airway_section` rows of the section length and area in SI.
fn story_u_rows() -> Vec<String> {
    let mut areas = Vec::new();
    let mut delta_cm = None;
    for line in STORY_2008.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0].starts_with("Δ") {
            delta_cm = Some(fields.last().unwrap().parse::<f64>().unwrap());
        } else if let Ok(section) = fields[0].parse::<usize>() {
            assert_eq!(section, areas.len() + 1);
            areas.push(fields.last().unwrap().parse::<f64>().unwrap());
        }
    }
    assert_eq!(areas.len(), 44);
    let length_m = delta_cm.expect("the section length row") * 1.0e-2;
    areas.iter().rev().map(|cm2| format!("airway_section {length_m} {}", cm2 * 1.0e-4)).collect()
}

/// The cut card with the airway row after `jet_noise`, or after `jet_limits` without one.
fn with_row(text: &str) -> String {
    let rows = story_u_rows();
    let record = format!("upstream_airway weibel-1963 saab-2021-table-1-B {GENERATIONS} {} story-2008-table-I-u\n{}\n{WALL}\n",
        rows.len(), rows.join("\n"));
    let anchor = if text.contains("\njet_noise ") { "\njet_noise " } else { "\njet_limits " };
    let start = text.find(anchor).expect("an anchor row") + 1;
    let end = start + text[start..].find('\n').expect("anchor row end") + 1;
    format!("{}{record}{}", &text[..end], &text[end..])
}

/// The cut card with the airway row ending in the `two-way` token (MB 18587).
fn with_row_two_way(text: &str) -> String {
    let one_way = with_row(text);
    assert_eq!(one_way.matches(" story-2008-table-I-u\n").count(), 1);
    one_way.replacen(" story-2008-table-I-u\n", " story-2008-table-I-u two-way\n", 1)
}

/// The rendered window's digest, the voice's last airway frame and the samples the windway passed no flow.
fn render(text: &str) -> (u64, Option<FlueAirwayFrame>, Option<FlueAirwayClosures>) {
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate()).unwrap_or_else(|e| panic!("{e}"));
    let record = file.upstream_airway().cloned();
    assert_eq!(record.is_some(), file.renderer().voice().upstream_airway().is_some());
    let mut renderer = file.into_renderer();
    let mut pressure = vec![0.0_f64; WINDOW];
    for block in pressure.chunks_mut(256) { renderer.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    (digest(&pressure), renderer.voice().upstream_airway_frame().copied(), renderer.voice().upstream_airway_closures())
}

fn plain(name: &str) -> u64 { PLAIN.iter().find(|(n, _)| *n == name).expect("a recorded digest").1 }

#[test]
fn cards_without_the_row_keep_their_bits() {
    for (name, want) in PLAIN {
        let (got, frame, closed) = render(&card(name));
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-plain\",\"card\":\"{name}\",\"digest\":\"{got:016x}\"}}");
        assert_eq!((frame, closed), (None, None), "{name}");
        assert_eq!(got, want, "{name}: without the row the bits moved ({got:016x} against {want:016x})");
    }
}

/// The last voice sample (the card's 288 kHz clock) at which the attack may close the windway: half the window, whose
/// 2048 output samples at 48 kHz are 12288 voice samples.
const ATTACK_CLOSURES_END: u64 = 6144;

/// The attack check (MB 18574) on `c`, the no-flow samples of a window with the row: no no-flow sample after
/// `ATTACK_CLOSURES_END`. `Err` names the last one.
fn attack_closures_end(c: &FlueAirwayClosures) -> Result<(), String> {
    match c.last_sample {
        Some(last) if last > ATTACK_CLOSURES_END => {
            Err(format!("a no-flow sample at {last}, after the attack's end {ATTACK_CLOSURES_END}: {c:?}"))
        }
        _ => Ok(()),
    }
}

/// With the row the window's bits move, and its no-flow samples are judged (MB 18430/18441/18574). Before the source
/// FIR's arrival `tau` the least open-circuit pressure is at least `-B` times the closed-end jump, `B` the deepest
/// undershoot before its arrival of the ideal lossless band-limited step sampled at this chain's own `tau`
/// (`sampled_step_undershoot`; a causal chain has nothing before its arrival, and a band-limited realization with loss
/// and a rolloff rings no more than the lossless brickwall one at the same fraction; over every fraction it is at most
/// 0.1412 for `tau` of 20 samples or more). The jump is taken from the chain itself: the realized step height (the
/// source FIR's running sum at its first maximum after the arrival) times the drive of the largest lung pressure before
/// the arrival (`closed_end_jump_pa`); twice that lung pressure is printed beside it, being the jump only where the
/// step height is one.
///
/// From the arrival on, THIS ASSERTION CHANGED (MB 18574). It was: no closure below -1e-6 of the lung pressure. That
/// assumed no reverse drive, and Story's measured /u/ tract produces one: the lung front reaches the closed windway
/// amplified (a step height near 18), and once the windway opens its flow rings the tract until the mouth-side pressure
/// falls below ambient. The windway is one-way, so air drawn back has no law here and it closes: on recorder c2, 269
/// no-flow samples from the arrival with the least open-circuit pressure near -164 Pa, all within the first 4.39 ms.
/// The old assertion encoded a false expectation, not a physics property. Now: no no-flow sample after
/// `ATTACK_CLOSURES_END` (`attack_closures_end`, shown to bite by `the_attack_closure_check_bites`), and the count,
/// runs, span and least are printed. The -1e-6 bound stays on the library's held-step tests (`ua_flue_001`,
/// `ua_flue_002`).
#[test]
fn the_story_u_row_reaches_the_window() {
    for name in ["recorder-c2-held", "recorder-f1"] {
        let text = with_row(&card(name));
        let (got, frame, closed) = render(&text);
        let (frame, closed) = (frame.expect("a frame"), closed.expect("closures"));
        let file = FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate()).unwrap_or_else(|e| panic!("{e}"));
        let record = file.upstream_airway().expect("the record").clone();
        assert_eq!((record.generations, record.duct_sections, record.source.as_str()), (24, 44, "story-2008-table-I-u"));
        assert_eq!(record.chain.sections().len(), 44 + 24);
        let voice = file.renderer().voice();
        let airway = voice.upstream_airway().expect("a realized airway");
        let receipt = airway.receipt().clone();
        assert_eq!(airway.sample_rate_hz(), 288_000, "{name}: ATTACK_CLOSURES_END counts the 288 kHz clock");
        let tau = closed.arrival_samples;
        assert_eq!(voice.upstream_airway_closures().expect("closures").arrival_samples, tau, "{name}: the arrival");
        let (b, step) = (sampled_step_undershoot(tau), source_step_height(airway, tau));
        let jump = closed_end_jump_pa(airway, tau, closed.lung_before_arrival_pa);
        drop(file);
        let index = |i: Option<u64>| i.map_or_else(|| "null".to_owned(), |i| i.to_string());
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-row\",\"card\":\"{name}\",\"digest\":\"{got:016x}\",\"zs\":{:e},\"r_dc\":{:e},\"r0\":{:e},\"series_impedance\":{:e},\"reflectance_taps\":{},\"source_taps\":{},\"lung_pa\":{:e},\"open_circuit_pa\":{:e},\"foot_pa\":{:e},\"flow_m3_s\":{:e},\"valve_pa\":{:e},\"no_flow_samples\":{},\"runs\":{},\"first_sample\":{},\"last_sample\":{},\"attack_closures_end\":{ATTACK_CLOSURES_END},\"least_open_circuit_pa\":{:e},\"arrival_samples\":{tau:e},\"undershoot_b\":{b:e},\"step_height\":{step:e},\"lung_before_arrival_pa\":{:e},\"jump_chain_pa\":{jump:e},\"jump_2_lung_pa\":{:e},\"bound_before_pa\":{:e},\"least_before_arrival_pa\":{:e},\"samples_from_arrival\":{},\"first_sample_from_arrival\":{},\"least_from_arrival_pa\":{:e}}}",
            receipt.port_impedance, receipt.dc_resistance, receipt.same_sample_reflection, receipt.series_impedance,
            receipt.reflectance.n_fft, receipt.source_taps, frame.lung_pa, frame.port.open_circuit_pa, frame.foot_pa,
            frame.flow_m3_s, frame.valve_pa, closed.samples, closed.runs, index(closed.first_sample),
            index(closed.last_sample), closed.least_open_circuit_pa, closed.lung_before_arrival_pa,
            2.0 * closed.lung_before_arrival_pa, -b * jump, closed.least_before_arrival_pa, closed.samples_from_arrival,
            index(closed.first_sample_from_arrival), closed.least_from_arrival_pa);
        assert_ne!(got, plain(name), "{name}: the row did not reach the window");
        assert!(tau >= 20.0 && b <= 0.1412, "{name}: B {b:e} at tau {tau:e}");
        assert!(closed.least_before_arrival_pa >= -b * jump, "{name}: before the arrival, past B of the chain's jump {closed:?}");
        if let Err(e) = attack_closures_end(&closed) { panic!("{name}: {e}"); }
        assert!(frame.port.impedance > 0.0 && frame.foot_pa > 0.0 && frame.foot_pa < frame.port.open_circuit_pa);
        assert!((frame.valve_pa - frame.foot_pa).abs() <= 1e-9 * frame.port.open_circuit_pa);
    }
}

/// The attack check bites (MB 18574): over a synthetic open-circuit trace of the window's 12288 voice samples, closed
/// (not positive) on samples 0 to 390 and on every third sample from 391 to 1264 like the c2 attack, the check passes;
/// with one more closure at 8000 it fails, and naming 8000. A closure at `ATTACK_CLOSURES_END` itself passes, one sample
/// later fails.
#[test]
fn the_attack_closure_check_bites() {
    // The trace's closures, tallied into the record's fields the way the voice counts them on its own clock.
    let tally = |extra: Option<u64>| {
        let mut c = FlueAirwayClosures { arrival_samples: 390.938_989_024_689_2, ..FlueAirwayClosures::default() };
        for n in 0..12_288_u64 {
            let open_circuit_pa = if n <= 390 || (n <= 1264 && n % 3 == 0) || Some(n) == extra { -1.0 } else { 1.0 };
            if open_circuit_pa > 0.0 { continue; }
            if c.last_sample.is_none_or(|last| last + 1 != n) { c.runs += 1; }
            c.samples += 1;
            c.first_sample.get_or_insert(n);
            c.last_sample = Some(n);
            c.least_open_circuit_pa = c.least_open_circuit_pa.min(open_circuit_pa);
            if f64::from(u32::try_from(n).unwrap()) >= c.arrival_samples {
                c.samples_from_arrival += 1;
                c.first_sample_from_arrival.get_or_insert(n);
                c.least_from_arrival_pa = c.least_from_arrival_pa.min(open_circuit_pa);
            }
        }
        c
    };
    let attack = tally(None);
    assert_eq!((attack.first_sample, attack.last_sample, attack.first_sample_from_arrival), (Some(0), Some(1263), Some(393)));
    assert!(attack_closures_end(&attack).is_ok(), "{attack:?}");
    let late = attack_closures_end(&tally(Some(8000))).expect_err("a closure at 8000 passed the attack check");
    assert!(late.contains("at 8000"), "{late}");
    assert!(attack_closures_end(&tally(Some(ATTACK_CLOSURES_END))).is_ok());
    assert!(attack_closures_end(&tally(Some(ATTACK_CLOSURES_END + 1))).is_err());
    println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-attack-check\",\"end\":{ATTACK_CLOSURES_END},\"attack_samples\":{},\"attack_runs\":{},\"late\":\"{}\"}}",
        attack.samples, attack.runs, late.split(':').next().unwrap_or(""));
}

/// The one-way row's window digests, recorded on the one-way tree (fa17, export 2b00d383c) before the two-way windway
/// existed: without the `two-way` token the row keeps them (MB 18587, test d).
const ONE_WAY: [(&str, u64); 2] = [("recorder-c2-held", 0x1ce0_c8b0_47d5_56d3), ("recorder-f1", 0xee87_e8e8_0a88_a49c)];

/// Without the `two-way` token the Story /u/ row renders the one-way windows bit for bit (MB 18587, test d).
#[test]
fn the_row_without_the_token_keeps_its_one_way_bits() {
    for (name, want) in ONE_WAY {
        let (got, _, closed) = render(&with_row(&card(name)));
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-one-way\",\"card\":\"{name}\",\"digest\":\"{got:016x}\"}}");
        assert_eq!(closed.expect("closures").least_flow_m3_s, 0.0, "{name}: one way, no flow drawn back");
        assert_eq!(got, want, "{name}: without the token the bits moved ({got:016x} against {want:016x})");
    }
}

/// With the `two-way` token the Story /u/ row draws air back (MB 18587, test e). On both recorder windows the attack's
/// negative open-circuit pressures pass reverse flow, counted as before and none after `ATTACK_CLOSURES_END`; the
/// most negative flux is negative; the window's bits differ from the one-way window's; the last frame passes forward
/// flow, finite. The line prints the counts, the span, the least open-circuit pressure and the least flux beside the
/// one-way values recorded in `ONE_WAY`.
#[test]
fn the_story_u_row_two_way_draws_air_back() {
    for (name, one_way) in ONE_WAY {
        let text = with_row_two_way(&card(name));
        let file = FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate()).unwrap_or_else(|e| panic!("{e}"));
        assert!(file.upstream_airway().expect("the record").two_way && file.renderer().voice().upstream_airway_two_way());
        drop(file);
        let (got, frame, closed) = render(&text);
        let (frame, c) = (frame.expect("a frame"), closed.expect("closures"));
        let index = |i: Option<u64>| i.map_or_else(|| "null".to_owned(), |i| i.to_string());
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-two-way\",\"card\":\"{name}\",\"digest\":\"{got:016x}\",\"one_way_digest\":\"{one_way:016x}\",\"reverse_or_no_flow_samples\":{},\"runs\":{},\"first_sample\":{},\"last_sample\":{},\"least_open_circuit_pa\":{:e},\"least_flow_m3_s\":{:e},\"samples_from_arrival\":{},\"least_before_arrival_pa\":{:e},\"least_from_arrival_pa\":{:e},\"lung_pa\":{:e},\"foot_pa\":{:e},\"flow_m3_s\":{:e}}}",
            c.samples, c.runs, index(c.first_sample), index(c.last_sample), c.least_open_circuit_pa, c.least_flow_m3_s,
            c.samples_from_arrival, c.least_before_arrival_pa, c.least_from_arrival_pa, frame.lung_pa, frame.foot_pa,
            frame.flow_m3_s);
        assert!(c.samples > 0 && c.least_flow_m3_s < 0.0, "{name}: no air drawn back {c:?}");
        if let Err(e) = attack_closures_end(&c) { panic!("{name}: {e}"); }
        assert_ne!(got, one_way, "{name}: the token did not reach the window");
        assert!(frame.flow_m3_s > 0.0 && frame.foot_pa > 0.0 && frame.foot_pa < frame.port.open_circuit_pa);
    }
}

/// The render binary's sidecar states the voice's own no-flow samples and the airway's receipt for the slit card with
/// the row, one way and two ways, the same as the library's render of the same text, and carries no airway record
/// without the row.
#[test]
fn the_render_sidecar_states_the_no_flow_samples() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-flue-airway-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    for (row, two) in [(true, false), (true, true), (false, false)] {
        let text = match (row, two) {
            (true, true) => with_row_two_way(&card("recorder-f1")),
            (true, false) => with_row(&card("recorder-f1")),
            _ => card("recorder-f1"),
        };
        let input = dir.join(format!("recorder-f1-{row}-{two}.performance"));
        std::fs::write(&input, &text).unwrap();
        let out = dir.join(format!("recorder-f1-{row}-{two}.wav"));
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("wind").arg(&input).arg(&out).args(["--block", "256"]).output().unwrap();
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
        let sidecar = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        let (_, _, closed) = render(&text);
        let Some(c) = closed else {
            assert!(!row && !sidecar.contains("upstream_airway"), "{sidecar}");
            continue;
        };
        let index = |i: Option<u64>| i.map_or_else(|| "null".to_owned(), |i| i.to_string());
        let want = format!("\"no_flow\":{{\"samples\":{},\"runs\":{},\"first_sample\":{},\"last_sample\":{},\"least_open_circuit_pa\":{:e},\"arrival_samples\":{:e},\"lung_before_arrival_pa\":{:e},\"least_before_arrival_pa\":{:e},\"samples_from_arrival\":{},\"first_sample_from_arrival\":{},\"least_from_arrival_pa\":{:e},\"least_flow_m3_s\":{:e},",
            c.samples, c.runs, index(c.first_sample), index(c.last_sample), c.least_open_circuit_pa, c.arrival_samples,
            c.lung_before_arrival_pa, c.least_before_arrival_pa, c.samples_from_arrival, index(c.first_sample_from_arrival),
            c.least_from_arrival_pa, c.least_flow_m3_s);
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-airway-sidecar\",\"card\":\"recorder-f1\",\"two_way\":{two},{}}}}}",
            &want[..want.len() - 1]);
        assert!(row && sidecar.contains(&format!(
            "\"upstream_airway\":{{\"source\":\"story-2008-table-I-u\",\"generations\":24,\"duct_sections\":44,\"sections\":68,\"wall\":true,\"two_way\":{two},")),
            "{sidecar}");
        assert!(sidecar.contains(&want), "{want} not in {sidecar}");
    }
}

#[test]
fn the_row_refuses_what_it_does_not_cover() {
    let text = with_row(&card("recorder-f1"));
    for (from, to) in [
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story", "upstream_airway weibel-1963 saab-2021-table-1-B 25 44 story"),
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story", "upstream_airway weibel-1963 saab-2021-table-1-A 24 44 story"),
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u", "upstream_airway weibel-1963 saab-2021-table-1-B 24 44"),
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story", "upstream_airway weibel-1963 saab-2021-table-1-B 24 45 story"),
        ("airway_wall 15 ", "airway_wall -15 "),
        // The `two-way` token (MB 18587): only as the last token, after an explicit source token.
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u\n", "upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u twoway\n"),
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u\n", "upstream_airway weibel-1963 saab-2021-table-1-B 24 44 two-way\n"),
        ("upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u\n", "upstream_airway weibel-1963 saab-2021-table-1-B 24 44 story-2008-table-I-u two-way two-way\n"),
    ] {
        assert_eq!(text.matches(from).count(), 1, "{from}");
        assert!(FlueFilePerformance::from_bytes(text.replacen(from, to, 1).as_bytes(), 256, &gate()).is_err(), "{to}");
    }
    // The card centreline speed is held on the held pressure: the slit card takes it without the row and refuses the
    // row by name with it.
    let (jet, centred) = ("\njet 90 1e-6 3 0\n", "\njet 90 1e-6 3 0\njet_speed card-centreline\n");
    assert_eq!(text.matches(jet).count(), 1);
    assert!(FlueFilePerformance::from_bytes(card("recorder-f1").replacen(jet, centred, 1).as_bytes(), 256, &gate()).is_ok());
    let error = FlueFilePerformance::from_bytes(text.replacen(jet, centred, 1).as_bytes(), 256, &gate())
        .err().expect("a refusal").to_string();
    assert!(error.contains(fs_couple::flue_loop::airway::AIRWAY_REFUSES_CARD_CENTRELINE), "{error}");
}
