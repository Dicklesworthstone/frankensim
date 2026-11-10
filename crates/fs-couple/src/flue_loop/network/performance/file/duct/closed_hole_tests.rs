//! MB 20929 (b), T1 to T4: the lumped closed hole (`duct_node closed-hole`)
//! against its closed form, against the one-step closed waveguide chimney it
//! stands in for, the waveguide path it leaves alone, and its refusals.
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
