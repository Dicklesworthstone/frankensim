//! MB 20929 (b), T1 to T4: the lumped closed hole (`duct_node closed-hole`)
//! against its closed form, against the one-step closed waveguide chimney it
//! stands in for, the waveguide path it leaves alone, and its refusals.
//! MB 23551, T5 to T8: its printed lengths behind `matching` and `series`.
//! MB 25086, T9 and T10: `series` is the paper's arm, a series mass only.
use super::*;
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::bernoulli_aperture::viscothermal::{
    ViscothermalTimeOptions, with_regional_viscothermal_sections_options,
};
use core::f64::consts::{PI, TAU};
use fs_material::gas::GasSpec;
use fs_vfit::waveguide::network::WaveguideNetwork;

const RATE: u32 = 288_000;
const BORE: f64 = 0.0061;
const HOLE: f64 = 0.00247;
const POLICY: FractionalTubePolicy = FractionalTubePolicy {
    maximum_frequency_hz: 3000.0,
    maximum_summed_phase_error_rad: 0.05,
};
// Each a whole number of samples per period at RATE.
const FREQUENCIES: [f64; 5] = [50.0, 500.0, 1000.0, 2000.0, 3000.0];

fn dt() -> f64 {
    1.0 / f64::from(RATE)
}

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

/// One step in the test gas, `c dt`, as the fractional owner computes its cell.
fn step_m() -> f64 {
    air().sound_speed * dt()
}

fn read(text: &str) -> Result<Input, FlueInputError> {
    let mut records = Vec::new();
    let mut rows = Rows { lines: text.lines(), line: 0 };
    let input = Input::read(&mut rows, air(), None, &mut records, &CancelGate::new_clock_free())?;
    assert!(rows.lines.next().is_none());
    Ok(input)
}

/// Inlet, a bore section, the hole's junction, a bore section to a matched end.
fn lumped(height_m: f64) -> String {
    format!(
        "network 3 2\nduct_node inlet\nduct_node closed-hole {:e} {:e}\nduct_node reflection 0\n\
         duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\n",
        HOLE, height_m, BORE, BORE
    )
}

/// The same bore with the hole as a closed waveguide chimney of `height_m`.
fn chimney(height_m: f64) -> String {
    format!(
        "network 4 3\nduct_node inlet\nduct_node junction\nduct_node reflection 0\nduct_node reflection 1\n\
         duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\nduct_section 1 3 {:e} {:e} 1e-9\n",
        BORE, BORE, height_m, HOLE
    )
}

fn network(text: &str) -> WaveguideNetwork {
    let mut input = read(text).unwrap();
    input.check_hole_steps(RATE).unwrap();
    input.graph.max_wave_memory_bytes = 8 << 20;
    let (network, _, _) = input
        .graph
        .realize_fractional_waveguide(&input.gases, dt(), POLICY, &CancelGate::new_clock_free())
        .unwrap();
    network
}

/// Drive the inlet with a unit sine at `hz`, let it settle, and return the
/// steady phasors of what `observe` reads after each step, projected on whole
/// periods. The inlet absorbs what returns and the far end is matched, so
/// nothing but the drive persists.
fn phasors<const N: usize>(
    network: &mut WaveguideNetwork,
    hz: f64,
    observe: impl Fn(&WaveguideNetwork) -> [f64; N],
) -> [[f64; 2]; N] {
    let period = f64::from(RATE) / hz;
    assert_eq!(period.fract(), 0.0);
    let period = period as usize;
    let span = 16_384usize.div_ceil(period) * period;
    let omega = TAU * hz / f64::from(RATE);
    let mut sums = [[0.0; 2]; N];
    for n in 0..2 * span {
        let phase = omega * (n % period) as f64;
        network.step(phase.sin()).unwrap();
        if n >= span {
            for (sum, value) in sums.iter_mut().zip(observe(&*network)) {
                sum[0] += value * phase.cos();
                sum[1] -= value * phase.sin();
            }
        }
    }
    sums
}

fn div(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let d = b[0] * b[0] + b[1] * b[1];
    [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]
}

/// `|a - b| / |b|`.
fn gap(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]) / b[0].hypot(b[1])
}

/// The closed form at the bilinear frequency: `j (w L - 1 / (w C))`,
/// `w = (2 / dt) tan(omega dt / 2)`.
fn bilinear(hz: f64, inertance: f64, compliance: f64) -> [f64; 2] {
    let w = 2.0 / dt() * (PI * hz / f64::from(RATE)).tan();
    [0.0, w * inertance - 1.0 / (w * compliance)]
}

fn law(height_m: f64) -> (f64, f64) {
    let gas = air();
    let area = PI * HOLE * HOLE;
    (
        gas.density * height_m / (3.0 * area),
        area * height_m / (gas.density * gas.sound_speed * gas.sound_speed),
    )
}

#[test]
fn t1_the_closed_hole_is_the_closed_form_shunt_and_the_engine_plays_its_bilinear_impedance() {
    let height_m = 0.3e-3;
    let input = read(&lumped(height_m)).unwrap();
    let (inertance, compliance) = law(height_m);
    let near = |a: f64, b: f64| (a - b).abs() <= 4.0 * f64::EPSILON * b.abs();
    let [(hole, _)] = input.closed.as_slice() else { panic!("one closed hole") };
    assert_eq!((hole.node, hole.hole_radius_m, hole.height_m), (1, HOLE, height_m));
    assert!(near(hole.inertance_pa_s2_m3, inertance), "{} {inertance}", hole.inertance_pa_s2_m3);
    assert!(near(hole.compliance_m3_pa, compliance), "{} {compliance}", hole.compliance_m3_pa);
    let NetworkNode::Shunt { load } = input.graph.nodes[1] else { panic!("the hole is a shunt") };
    let base = load.base();
    assert_eq!(base.resistance_pa_s_m3, 0.0);
    assert_eq!(base.inertance_pa_s2_m3, hole.inertance_pa_s2_m3);
    assert_eq!(base.compliance_m3_pa, Some(hole.compliance_m3_pa));
    assert_eq!(input.graph.sections.len(), 2);
    let mut network = network(&lumped(height_m));
    for hz in FREQUENCIES {
        let [p, u] = phasors(&mut network, hz, |n| {
            let f = n.node_frame(1).unwrap();
            [f.load_pressure_pa, f.load_flow_m3_s]
        });
        let engine = div(p, u);
        let closed = bilinear(hz, hole.inertance_pa_s2_m3, hole.compliance_m3_pa);
        let error = gap(engine, closed);
        println!("T1 {hz} Hz: engine shunt Z {:.9e} j, closed form {:.9e} j Pa s/m^3 (engine real part {:.3e}), relative difference {error:.3e}",
            engine[1], closed[1], engine[0]);
        assert!(error <= 1e-9, "{hz} Hz: {error:e}");
    }
}

#[test]
fn t2_at_one_step_the_lumped_hole_meets_the_closed_waveguide_chimney() {
    let gas = air();
    let step = step_m();
    let below = step * (1.0 - 1e-9);
    // The chimney at exactly one step is one sample each way; the lumped hole sits just below.
    let mut lump = network(&lumped(below));
    let mut wave = network(&chimney(step));
    let z_hole = gas.density * gas.sound_speed / (PI * HOLE * HOLE);
    let (inertance, compliance) = law(below);
    println!("T2 crossover: one step c dt = {step:.9e} m (rate {RATE} Hz), lumped hole at {below:.9e} m, waveguide chimney at {step:.9e} m");
    for hz in FREQUENCIES {
        let x = TAU * hz / f64::from(RATE);
        let [p, u, junction_lump] = phasors(&mut lump, hz, |n| {
            let f = n.node_frame(1).unwrap();
            [f.load_pressure_pa, f.load_flow_m3_s, f.pressure_pa]
        });
        let [junction_wave] = phasors(&mut wave, hz, |n| [n.node_frame(1).unwrap().pressure_pa]);
        let z_lump = div(p, u);
        // A one-sample closed line: Z_h (1 + z^-2) / (1 - z^-2) = -j Z_h cot(omega dt).
        let z_wave = [0.0, -z_hole / x.tan()];
        let tau = (0.5 * x).tan();
        let analytic = tau * tau / 3.0 / (1.0 - tau * tau);
        let branch = gap(z_lump, z_wave);
        let junction = gap(junction_lump, junction_wave);
        let closed = gap(z_lump, bilinear(hz, inertance, compliance));
        println!("T2 {hz} Hz: branch |Z_lump - Z_chimney| / |Z_chimney| = {branch:.6e} (analytic (tau^2/3)/(1-tau^2) = {analytic:.6e}); \
                  junction pressure |p_lump - p_chimney| / |p_chimney| = {junction:.6e}; lumped engine vs its closed form {closed:.3e}");
        assert!(closed <= 1e-9, "{hz} Hz: {closed:e}");
        // The lumped hole's height is 1e-9 below the step, which moves L and C by 1e-9.
        assert!((branch - analytic).abs() <= 1e-6 * analytic + 4e-9, "{hz} Hz: {branch:e} {analytic:e}");
        // With matched bore lines and a lossless branch, |dp/p| = |dY| / |2 Y_bore + Y_lump| <= |dZ| / |Z_chimney|.
        assert!(junction <= branch * (1.0 + 1e-6) + 1e-12, "{hz} Hz: {junction:e} {branch:e}");
    }
}

#[test]
fn t3_a_closed_hole_of_one_step_or_more_is_a_duct_section_and_an_open_hole_below_one_step_refuses() {
    let step = step_m();
    for height_m in [step, step * (1.0 + 1e-9), 2.0 * step] {
        let error = read(&lumped(height_m)).unwrap().check_hole_steps(RATE).unwrap_err();
        assert!(error.to_string().contains("closed-hole of one step or more"), "{error}");
    }
    read(&lumped(step * (1.0 - 1e-9))).unwrap().check_hole_steps(RATE).unwrap();
    // A card without the node keeps its chimney section and its closed end, and no record.
    let input = read(&chimney(step)).unwrap();
    assert!(input.closed.is_empty());
    assert_eq!(input.graph.nodes, vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 0.0 }, NetworkNode::Termination { reflection: 1.0 }]);
    assert_eq!(input.graph.sections.len(), 3);
    input.check_hole_steps(RATE).unwrap();
    // An open side hole: 0.3 mm plus its inner length is above one step at 288 kHz, below it at 96 kHz.
    let open = format!(
        "network 4 3\nduct_node inlet\nduct_node junction\nduct_node reflection 0\nduct_node reflection -1\n\
         duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\nduct_section 1 3 3e-4 {:e} 1e-9 inner {:e}\n",
        BORE, BORE, HOLE, BORE
    );
    let input = read(&open).unwrap();
    println!("T3 open hole: graph chimney {:.6e} m, one step {:.6e} m at {RATE} Hz and {:.6e} m at 96000 Hz",
        input.inner[0].network_length_m, step, air().sound_speed / 96_000.0);
    input.check_hole_steps(RATE).unwrap();
    let error = input.check_hole_steps(96_000).unwrap_err();
    assert!(error.to_string().contains("open side hole shorter than one step"), "{error}");
}

#[test]
fn t4_zero_or_negative_length_or_radius_refuses_on_the_node_and_on_the_section() {
    for (b, t) in [(0.0, 3e-4), (HOLE, 0.0), (-HOLE, 3e-4), (HOLE, -3e-4)] {
        let text = format!(
            "network 3 2\nduct_node inlet\nduct_node closed-hole {b:e} {t:e}\nduct_node reflection 0\n\
             duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\n",
            BORE, BORE
        );
        let error = read(&text).err().expect("refused");
        assert!(error.to_string().contains("closed-hole needs a positive hole radius and height"), "{error}");
    }
    for (t, b) in [(0.0, HOLE), (3e-4, 0.0)] {
        let text = format!(
            "network 4 3\nduct_node inlet\nduct_node junction\nduct_node reflection 0\nduct_node reflection 1\n\
             duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\nduct_section 1 3 {t:e} {b:e} 1e-9\n",
            BORE, BORE
        );
        let error = read(&text).err().expect("refused");
        assert!(error.to_string().contains("positive geometry"), "{error}");
    }
}

#[test]
fn the_closed_hole_keeps_its_shunt_through_the_shared_loss_subdivision() {
    let gate = CancelGate::new_clock_free();
    let mut input = read(&lumped(3e-4)).unwrap();
    input.check_hole_steps(RATE).unwrap();
    // 16 cells: each 6.25 mm cell within 0.5 rad at 3 kHz, each quarter above one step.
    let mut rows = Rows {
        lines: "viscothermal sections 2\nloss_section 0 16 50 3000\nloss_section 1 16 50 3000\n".lines(),
        line: 0,
    };
    let (selections, sampling) = input.read_losses(&mut rows).unwrap();
    assert!(rows.lines.next().is_none());
    input.graph.max_wave_memory_bytes = 8 << 20;
    let shunt = input.graph.nodes[1];
    let (graph, gases, _) = with_regional_viscothermal_sections_options(
        input.graph,
        input.gases,
        dt(),
        &selections,
        ViscothermalTimeOptions { sampling, fractional: Some(POLICY) },
        &gate,
    )
    .unwrap();
    // The shunt stays once, at its own node: subdivision appends its new nodes, and the receipt names node 1.
    assert_eq!(graph.nodes.iter().filter(|&&node| node == shunt).count(), 1);
    assert_eq!(graph.nodes[1], shunt, "the closed hole's node index survives the loss subdivision");
    let (mut network, _, _) = graph.realize_fractional_waveguide(&gases, dt(), POLICY, &gate).unwrap();
    let (mut heat, mut work) = (0.0, 0.0);
    for sample in 0..4096 {
        let frame = network.step((TAU * f64::from(sample) / 576.0).sin()).unwrap();
        assert!(frame.interior_loss_j >= 0.0);
        heat += frame.terminal_loss_j + frame.interior_loss_j;
        work += frame.inlet_work_j;
    }
    assert!(heat > 0.0);
    assert!((network.stored_energy_j() + heat - work).abs() < 3e-10 * work.abs());
}

#[test]
fn the_closed_hole_refuses_a_lone_section_a_wider_hole_and_a_side_hole_on_its_junction() {
    let lone = format!(
        "network 2 1\nduct_node inlet\nduct_node closed-hole {:e} 3e-4\nduct_section 0 1 1e-1 {:e} 1e-9\n",
        HOLE, BORE
    );
    let error = read(&lone).err().expect("refused");
    assert!(error.to_string().contains("at least two sections"), "{error}");
    let wide = lumped(3e-4).replace(&format!("closed-hole {:e}", HOLE), "closed-hole 7e-3");
    let error = read(&wide).err().expect("refused");
    assert!(error.to_string().contains("no wider than its junction's sections"), "{error}");
    let shared = format!(
        "network 4 3\nduct_node inlet\nduct_node closed-hole {:e} 3e-4\nduct_node reflection 0\nduct_node reflection -1\n\
         duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\nduct_section 1 3 3e-4 {:e} 1e-9 inner {:e}\n",
        HOLE, BORE, BORE, HOLE, BORE
    );
    let error = read(&shared).err().expect("refused");
    assert!(error.to_string().contains("cannot share its junction"), "{error}");
}

// MB 23551: the printed closed-hole lengths behind `matching` and `series`, the open hole's words.

/// Eq. 8 of Lefebvre and Scavone 2012, by hand: `t_m = b d (1 + 0.207 d^3) / 8`, `d = b / a`.
fn printed_matching(b: f64, a: f64) -> f64 {
    let d = b / a;
    b * d * (1.0 + 0.207 * d * d * d) / 8.0
}

/// Eq. 34, by hand: the closed series length `t_a = (-0.12 - 0.17 tanh(2.4 t / b)) b d^2`.
fn printed_closed_series(b: f64, a: f64, t: f64) -> f64 {
    let d = b / a;
    (-0.12 - 0.17 * (2.4 * t / b).tanh()) * b * d * d
}

/// The lumped card with `words` and the bore radius a after them; the second bore section's radius is `r2`.
fn named(height_m: f64, words: &str, r2: f64) -> String {
    format!(
        "network 3 2\nduct_node inlet\nduct_node closed-hole {:e} {:e} {words} {:e}\nduct_node reflection 0\n\
         duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\n",
        HOLE, height_m, BORE, BORE, r2
    )
}

#[test]
fn t5_matching_takes_the_shunt_at_t_plus_t_m_and_the_engine_plays_it() {
    let height_m = 0.3e-3;
    let input = read(&named(height_m, "matching", BORE)).unwrap();
    let [(hole, _)] = input.closed.as_slice() else { panic!("one closed hole") };
    let t_m = printed_matching(HOLE, BORE);
    let got = hole.matching_length_m.expect("matching");
    assert!((got - t_m).abs() <= 4.0 * f64::EPSILON * t_m, "fs-duct's t_m {got:e} against Eq. 8 {t_m:e}");
    assert_eq!((hole.height_m, hole.bore_radius_m, hole.series_length_m, hole.series_bore),
        (height_m, Some(BORE), None, None));
    assert_eq!(hole.lumped_height_m(), height_m + got);
    let (inertance, compliance) = law(hole.lumped_height_m());
    let near = |a: f64, b: f64| (a - b).abs() <= 4.0 * f64::EPSILON * b.abs();
    assert!(near(hole.inertance_pa_s2_m3, inertance) && near(hole.compliance_m3_pa, compliance));
    // The hole's volume grows by t_m / t; matching alone moves no bore section.
    let plain = read(&lumped(height_m)).unwrap();
    let ratio = plain.closed[0].0.compliance_m3_pa / hole.compliance_m3_pa;
    println!("T5 t_m {got:.6e} m ({:.1} % of t), C without the word / with it {ratio:.6}", 100.0 * got / height_m);
    assert!((ratio - height_m / (height_m + got)).abs() <= 1e-12, "{ratio}");
    assert!(input.graph.sections.iter().all(|s| s.length_m == 1e-1));
    let mut network = network(&named(height_m, "matching", BORE));
    for hz in FREQUENCIES {
        let [p, u] = phasors(&mut network, hz, |n| {
            let f = n.node_frame(1).unwrap();
            [f.load_pressure_pa, f.load_flow_m3_s]
        });
        let error = gap(div(p, u), bilinear(hz, hole.inertance_pa_s2_m3, hole.compliance_m3_pa));
        println!("T5 {hz} Hz: engine shunt against the closed form at t + t_m, relative difference {error:.3e}");
        assert!(error <= 1e-9, "{hz} Hz: {error:e}");
    }
}

#[test]
fn t6_with_matching_the_lumped_hole_meets_a_closed_chimney_t_plus_t_m_long_at_one_step() {
    let gas = air();
    let step = step_m();
    let t_m = printed_matching(HOLE, BORE);
    // The admission reads t + t_m: a hole just under t_m short of one step passes without the word and refuses with it.
    let short = step - t_m * (1.0 - 1e-6);
    read(&lumped(short)).unwrap().check_hole_steps(RATE).unwrap();
    let error = read(&named(short, "matching", BORE)).unwrap().check_hole_steps(RATE).unwrap_err();
    assert!(error.to_string().contains("closed-hole of one step or more"), "{error}");
    // t + t_m just below one step against the closed chimney written t + t_m long at one step: t_m on both paths.
    let below = step * (1.0 - 1e-9) - t_m;
    let input = read(&named(below, "matching", BORE)).unwrap();
    let hole = input.closed[0].0;
    let mut lump = network(&named(below, "matching", BORE));
    let mut wave = network(&chimney(step));
    let z_hole = gas.density * gas.sound_speed / (PI * HOLE * HOLE);
    println!("T6 crossover: one step {step:.9e} m, lumped t {below:.9e} + t_m {:.9e} = {:.9e} m, chimney {step:.9e} m",
        hole.matching_length_m.unwrap(), hole.lumped_height_m());
    for hz in FREQUENCIES {
        let x = TAU * hz / f64::from(RATE);
        let [p, u, junction_lump] = phasors(&mut lump, hz, |n| {
            let f = n.node_frame(1).unwrap();
            [f.load_pressure_pa, f.load_flow_m3_s, f.pressure_pa]
        });
        let [junction_wave] = phasors(&mut wave, hz, |n| [n.node_frame(1).unwrap().pressure_pa]);
        let z_lump = div(p, u);
        let z_wave = [0.0, -z_hole / x.tan()];
        let tau = (0.5 * x).tan();
        let analytic = tau * tau / 3.0 / (1.0 - tau * tau);
        let branch = gap(z_lump, z_wave);
        let junction = gap(junction_lump, junction_wave);
        println!("T6 {hz} Hz: branch {branch:.6e} (analytic {analytic:.6e}), junction pressure {junction:.6e}");
        assert!((branch - analytic).abs() <= 1e-6 * analytic + 4e-9, "{hz} Hz: {branch:e} {analytic:e}");
        assert!(junction <= branch * (1.0 + 1e-6) + 1e-12, "{hz} Hz: {junction:e} {branch:e}");
    }
}

#[test]
fn t7_series_takes_eq_34_off_the_two_bore_sections_by_the_open_holes_arithmetic() {
    let height_m = 0.3e-3;
    let r2 = 0.9 * BORE;
    for words in ["series", "matching series"] {
        let input = read(&named(height_m, words, r2)).unwrap();
        let [(hole, _)] = input.closed.as_slice() else { panic!("one closed hole") };
        let t_a = printed_closed_series(HOLE, BORE, height_m);
        let got = hole.series_length_m.expect("series");
        assert!(got < 0.0 && (got - t_a).abs() <= 8.0 * f64::EPSILON * t_a.abs(), "fs-duct's t_a {got:e} against Eq. 34 {t_a:e}");
        // The open `series` arithmetic: (|t_a| / 2) (r_i / a)^2 off each section, the same series mass in its own area.
        let share = |r: f64| -got / 2.0 * (r / BORE) * (r / BORE);
        assert_eq!(hole.series_bore, Some([(0, share(BORE)), (1, share(r2))]));
        assert_eq!(input.graph.sections[0].length_m, 1e-1 - share(BORE));
        assert_eq!(input.graph.sections[1].length_m, 1e-1 - share(r2));
        let mass = share(BORE) / (BORE * BORE) + share(r2) / (r2 * r2);
        assert!((mass - got.abs() / (BORE * BORE)).abs() <= 1e-12 * mass, "{mass:e}");
        assert_eq!(hole.matching_length_m.is_some(), words.starts_with("matching"));
        if hole.matching_length_m.is_none() {
            assert_eq!(hole.lumped_height_m().to_bits(), height_m.to_bits(), "series alone leaves the shunt at t");
        }
        println!("T7 {words}: t_a {got:.6e} m, bore sections shorter by {:.6e} and {:.6e} m", share(BORE), share(r2));
    }
}

#[test]
fn t8_the_words_refuse_by_name_and_without_them_nothing_moves() {
    let height_m = 0.3e-3;
    let gas = air();
    for (text, refusal) in [
        (lumped(height_m).replace(&format!("{:e}\nduct_node reflection", height_m),
            &format!("{:e} {:e}\nduct_node reflection", height_m, BORE)), "takes a bore radius only after matching or series"),
        (lumped(height_m).replace(&format!("{:e}\nduct_node reflection", height_m),
            &format!("{:e} matching\nduct_node reflection", height_m)), "needs the main bore's radius after the words"),
        (named(height_m, "matching", BORE).replace(&format!("matching {:e}", BORE), &format!("matching {:e}", HOLE)),
            "narrower than the bore"),
        (named(height_m, "series", BORE).replace(&format!("series {:e}", BORE), &format!("series {:e}", -BORE)),
            "narrower than the bore"),
        (named(height_m, "series", BORE).replace(&format!("series {:e}", BORE), "series inf"), "physical scalars must be finite"),
        (named(height_m, "series", BORE).replace("duct_section 1 2 1e-1", "duct_section 1 2 1e-5"),
            "longer than its share of t_a"),
        (format!(
            "network 4 3\nduct_node inlet\nduct_node closed-hole {:e} {height_m:e} series {:e}\nduct_node reflection 0\nduct_node reflection 0\n\
             duct_section 0 1 1e-1 {:e} 1e-9\nduct_section 1 2 1e-1 {:e} 1e-9\nduct_section 1 3 1e-1 {:e} 1e-9\n",
            HOLE, BORE, BORE, BORE, BORE), "exactly two bore sections"),
        (format!("network 2 1\nduct_node inlet\nduct_node closed-hole {:e} 3e-4 matching {:e}\nduct_section 0 1 1e-1 {:e} 1e-9\n",
            HOLE, BORE, BORE), "at least two sections"),
    ] {
        let error = read(&text).err().unwrap_or_else(|| panic!("accepted: {text}"));
        assert!(error.to_string().contains(refusal), "refused as `{error}`, expected `{refusal}`");
    }
    // Without the words: no record, the landed law's L and C bit for bit, the bore untouched.
    let input = read(&lumped(height_m)).unwrap();
    let hole = input.closed[0].0;
    assert_eq!((hole.bore_radius_m, hole.matching_length_m, hole.series_length_m, hole.series_bore), (None, None, None, None));
    assert_eq!(hole.series_volume_m3, None);
    assert_eq!(hole.lumped_height_m().to_bits(), height_m.to_bits());
    let (inertance, compliance) = ClosedHoleLump::law(HOLE, height_m, gas.density, gas.sound_speed).unwrap();
    assert_eq!((hole.inertance_pa_s2_m3.to_bits(), hole.compliance_m3_pa.to_bits()), (inertance.to_bits(), compliance.to_bits()));
    assert!(input.graph.sections.iter().all(|s| s.length_m.to_bits() == 1e-1f64.to_bits()));
}

#[test]
fn t9_series_gives_the_lost_bore_volume_back_to_the_hole_so_the_junction_loses_mass_only() {
    let height_m = 0.3e-3;
    let r2 = 0.9 * BORE;
    let gas = air();
    let stiffness = gas.density * gas.sound_speed * gas.sound_speed;
    let near = |a: f64, b: f64| (a - b).abs() <= 4.0 * f64::EPSILON * b.abs();
    // The junction's volume: both bore sections and the hole's C rho c^2. Its series mass: rho l / S over the sections.
    let volume = |input: &Input| -> f64 {
        input.graph.sections.iter().map(|s| PI * s.radius_m * s.radius_m * s.length_m).sum::<f64>()
            + input.closed[0].0.compliance_m3_pa * stiffness
    };
    let mass = |input: &Input| -> f64 {
        input.graph.sections.iter().map(|s| gas.density * s.length_m / (PI * s.radius_m * s.radius_m)).sum()
    };
    let plain = read(&named(height_m, "matching", r2)).unwrap();
    let input = read(&named(height_m, "matching series", r2)).unwrap();
    let hole = input.closed[0].0;
    let t_a = hole.series_length_m.expect("series");
    let given = hole.series_volume_m3.expect("series");
    let [(i, x), (j, y)] = hole.series_bore.expect("series");
    let sections = &input.graph.sections;
    let lost = PI * (sections[i].radius_m * sections[i].radius_m * x + sections[j].radius_m * sections[j].radius_m * y);
    assert!(near(given, lost), "V_a {given:e} against the sections' lost volume {lost:e}");
    let (v0, v1) = (volume(&plain), volume(&input));
    assert!((v1 - v0).abs() <= 16.0 * f64::EPSILON * v0, "junction volume {v0:e} without series, {v1:e} with it");
    let (dm, want) = (mass(&plain) - mass(&input), gas.density * t_a.abs() / (PI * BORE * BORE));
    assert!((dm - want).abs() <= 1e-10 * want, "series mass off the junction {dm:e} against rho |t_a| / S_a {want:e}");
    // The hole's own law at t + t_m is untouched but for V_a in C.
    let (inertance, compliance) = law(hole.lumped_height_m());
    assert_eq!(hole.inertance_pa_s2_m3.to_bits(), plain.closed[0].0.inertance_pa_s2_m3.to_bits());
    assert!(near(hole.inertance_pa_s2_m3, inertance) && near(hole.compliance_m3_pa, compliance + given / stiffness));
    println!("T9 V_a {given:.6e} m^3 against the hole's {:.6e}; junction volume moved {:.3e} (relative), mass off {dm:.6e}",
        compliance * stiffness, (v1 - v0) / v0);
    // The engine plays the shunt with V_a in C.
    let mut network = network(&named(height_m, "matching series", r2));
    for hz in FREQUENCIES {
        let [p, u] = phasors(&mut network, hz, |n| {
            let f = n.node_frame(1).unwrap();
            [f.load_pressure_pa, f.load_flow_m3_s]
        });
        let error = gap(div(p, u), bilinear(hz, hole.inertance_pa_s2_m3, hole.compliance_m3_pa));
        println!("T9 {hz} Hz: engine shunt against the closed form with V_a, relative difference {error:.3e}");
        assert!(error <= 1e-9, "{hz} Hz: {error:e}");
    }
}

type Cx = [f64; 2];
type Two = [[Cx; 2]; 2];

fn cmul(a: Cx, b: Cx) -> Cx {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn chain(a: Two, b: Two) -> Two {
    let entry = |r: usize, c: usize| {
        let (p, q) = (cmul(a[r][0], b[0][c]), cmul(a[r][1], b[1][c]));
        [p[0] + q[0], p[1] + q[1]]
    };
    [[entry(0, 0), entry(0, 1)], [entry(1, 0), entry(1, 1)]]
}

/// A lossless bore section of `length` (negative: taken off) and `area`: `[[cos, j Z sin], [j sin / Z, cos]]`.
fn section(length: f64, k: f64, area: f64, rho_c: f64) -> Two {
    let (z, (sin, cos)) = (rho_c / area, (k * length).sin_cos());
    [[[cos, 0.0], [0.0, z * sin]], [[0.0, sin / z], [cos, 0.0]]]
}

/// The hole's series L-C to zero pressure as a shunt: admittance `1 / j (w L - 1 / (w C))`.
fn hole_shunt(w: f64, inertance: f64, compliance: f64) -> Two {
    let x = w * inertance - 1.0 / (w * compliance);
    [[[1.0, 0.0], [0.0, 0.0]], [[0.0, -1.0 / x], [1.0, 0.0]]]
}

fn arm(z: Cx) -> Two {
    [[[1.0, 0.0], z], [[0.0, 0.0], [1.0, 0.0]]]
}

#[test]
fn t10_the_closed_junction_with_series_is_the_papers_arm_to_second_order() {
    let height_m = 0.3e-3;
    let r2 = 0.9 * BORE;
    let gas = air();
    let rho_c = gas.density * gas.sound_speed;
    let input = read(&named(height_m, "matching series", r2)).unwrap();
    let hole = input.closed[0].0;
    let t_a = hole.series_length_m.expect("series");
    let [(i, x), (j, y)] = hole.series_bore.expect("series");
    let area = |s: usize| PI * input.graph.sections[s].radius_m * input.graph.sections[s].radius_m;
    // C_0: the hole at t + t_m without V_a, the shunt of the paper's T (its -Z_a / 4 left out, as the row leaves it).
    let (_, c0) = law(hole.lumped_height_m());
    let (l, c) = (hole.inertance_pa_s2_m3, hole.compliance_m3_pa);
    // Eq. 34's shares and V_a from t_a and the sections' radii, not the engine's record, for the A and D accounting below.
    let share = |s: usize| -t_a / 2.0 * (input.graph.sections[s].radius_m / BORE).powi(2);
    let (x_p, y_p) = (share(i), share(j));
    let c_p = c0 + (area(i) * x_p + area(j) * y_p) / (gas.density * gas.sound_speed * gas.sound_speed);
    for hz in FREQUENCIES {
        let w = TAU * hz;
        let k = w / gas.sound_speed;
        // The paper's T: Z_a / 2, the shunt, Z_a / 2, with Z_a = j k Z_c t_a in the bore's area at the hole.
        let half = [0.0, k * rho_c / (PI * BORE * BORE) * t_a / 2.0];
        let paper = chain(chain(arm(half), hole_shunt(w, l, c0)), arm(half));
        let realised = |compliance: f64| chain(chain(section(-x, k, area(i), rho_c), hole_shunt(w, l, compliance)),
            section(-y, k, area(j), rho_c));
        let fixed = realised(c);
        let entry = |m: Two, r: usize, s: usize| {
            let d = [m[r][s][0] - paper[r][s][0], m[r][s][1] - paper[r][s][1]];
            d[0].hypot(d[1]) / if r == s { 1.0 } else { paper[r][s][0].hypot(paper[r][s][1]) }
        };
        // Named tolerance, second order in the junction's lengths: e2 = (k (|t_a| + t + t_m))^2 on A, B and D; on the
        // admittance C, V_a sits behind the hole column's mass L, not beside it, which leaves w^2 L (C^2 - C_0^2) / C_0,
        // plus the shortened sections' own (k |t_a|)^2.
        let e2 = (k * (t_a.abs() + hole.lumped_height_m())).powi(2);
        let behind = w * w * l * (c * c - c0 * c0) / c0;
        let bound_c = 1.1 * behind + 4.0 * (k * t_a).powi(2);
        let [b, cc] = [entry(fixed, 0, 1), entry(fixed, 1, 0)];
        let without = entry(realised(c0), 1, 0);
        // A and D are 1 + O((k l)^2) in the paper's T and in the realisation alike, so a bound on either's distance from
        // the paper at e2 cannot fail: the bare shunt, no arm at all, meets it. They are accounted for instead, signed, from
        // Eq. 34's shares and V_a: A = 1 + (Z_a / 2) Y(C_0 + V_a / rho c^2) - (k^2 / 2)(x^2 + y^2 + 2 (S_j / S_i) x y), D
        // the same with S_i / S_j: the arm's half-mass on the hole's admittance, less each shortened section's mass on the
        // volume the other loses. What is left is fourth order, under e2 of the paper's own |A - 1|; V_a left out, or the
        // shares doubled or swapped, leave 0.3 to 0.9 of it.
        let half_y = cmul(half, [0.0, -1.0 / (w * l - 1.0 / (w * c_p))]);
        let own = |ratio: f64| 0.5 * k * k * (x_p * x_p + y_p * y_p + 2.0 * ratio * x_p * y_p);
        let off = |m: Cx, less: f64| (m[0] - 1.0 - half_y[0] + less).hypot(m[1] - half_y[1]);
        let (a, d) = (off(fixed[0][0], own(area(j) / area(i))), off(fixed[1][1], own(area(i) / area(j))));
        let depart = (paper[0][0][0] - 1.0).hypot(paper[0][0][1]);
        let bound_ad = e2 * depart + 16.0 * f64::EPSILON;
        println!("T10 {hz} Hz: A {:.3e} D {:.3e} of the paper's |A - 1| {depart:.3e} (e2 {e2:.3e}); B {b:.3e}; C {cc:.3e} \
                  (V_a behind L {behind:.3e}, bound {bound_c:.3e}); without V_a, C {without:.3e}", a / depart, d / depart);
        assert!(a <= bound_ad && d <= bound_ad, "{hz} Hz: A {a:e} D {d:e} off the accounting, against {bound_ad:e}");
        assert!(b <= e2, "{hz} Hz: B {b:e} against {e2:e}");
        assert!(cc <= bound_c, "{hz} Hz: C {cc:e} against {bound_c:e}");
        // The shortening alone, without V_a, misses the arm's shunt admittance by most of the hole.
        assert!(without >= 0.5, "{hz} Hz: without V_a the admittance is off by only {without:e}");
    }
}
