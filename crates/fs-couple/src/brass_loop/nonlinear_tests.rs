use super::*;
use fs_duct::modal::mm_input_impedance;
use fs_duct::{LossModel, Segment, Termination};
use fs_gas::dynamics::CaloricDomain;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{
    ThermochemicalPhase, ThermochemicalSources, read_frozen_thermochemical_model, read_gas_record,
};
use fs_math::det;
use fs_qty::Temperature;

#[path = "../render/schedule/reed/field_input/tests/support.rs"]
mod gas_test_support;

const INLET_M2: f64 = 2.2e-4;
const OPENING_M2: f64 = 0.012 * 8.0e-4;

const RECORD: &str = "thermochemical 293.15 101325 mole 1 specimen-gas 1";

// A synthetic specimen gas through the production pack, model and caloric
// owners, with the material state every linear owner reads from that source.
fn specimen_with_state() -> (FrozenGas, GasState) {
    let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let (model, _) =
        read_frozen_thermochemical_model(&mut RECORD.split_ascii_whitespace(), &sources).unwrap();
    let state = *read_gas_record(&mut RECORD.split_ascii_whitespace(), None, Some(&sources))
        .unwrap()
        .state();
    let gas = FrozenGas::new(
        model.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
        || false,
    )
    .unwrap();
    (gas, state)
}

fn specimen() -> FrozenGas {
    specimen_with_state().0
}

fn state(gas: &FrozenGas, pressure: f64, velocity: f64) -> Primitive {
    let density = pressure / (gas.gas_constant() * 293.15);
    gas.primitive(gas.conserved(density, velocity, pressure).unwrap(), || {
        false
    })
    .unwrap()
}

#[test]
fn a_closed_lip_is_the_interiors_own_wall() {
    let gas = specimen();
    let interior = state(&gas, 101_375.0, 0.3);
    let mouth = state(&gas, 104_325.0, 0.0);
    let closed = lip_face(&gas, interior, mouth, 0.0, INLET_M2).unwrap();
    let impedance = interior.density * interior.sound_speed;
    assert!(closed.flux.mass.abs() <= 1e-12 * interior.density * interior.sound_speed);
    assert_eq!(closed.mass_flow_kg_s, closed.flux.mass * INLET_M2);
    // Gas leaving the closed end draws an expansion of rho c u.
    let acoustic = interior.pressure - impedance * 0.3;
    assert!(
        (closed.pressure_pa - acoustic).abs() <= 2e-3 * impedance * 0.3,
        "{} {acoustic}",
        closed.pressure_pa
    );
}

#[test]
fn an_open_lip_meters_the_mouth_through_continuity_with_the_bores_wave() {
    let gas = specimen();
    let interior = state(&gas, 101_325.0, 0.0);
    let mouth = state(&gas, 104_325.0, 0.0);
    let open = lip_face(&gas, interior, mouth, OPENING_M2, INLET_M2).unwrap();
    assert!(open.pressure_pa > interior.pressure && open.pressure_pa < mouth.pressure);
    let (jet, choked) = gas.nozzle_exit(mouth, open.pressure_pa, || false).unwrap();
    assert!(!choked);
    let metered = jet.density * jet.velocity * OPENING_M2;
    assert!(
        (open.mass_flow_kg_s - metered).abs() <= 1e-9 * metered,
        "{} {metered}",
        open.mass_flow_kg_s
    );
    // Entering gas carries exactly the mouth's stagnation enthalpy.
    assert_eq!(
        open.flux.energy,
        open.flux.mass * gas.stagnation_enthalpy(mouth.temperature).unwrap()
    );
    // The cup face rises by the interior's simple wave, rho c times the face speed.
    let face_density = open.pressure_pa / (gas.gas_constant() * mouth.temperature);
    let rise = open.pressure_pa - interior.pressure;
    let acoustic = interior.density * interior.sound_speed * open.flux.mass / face_density;
    println!(
        "open lip: face rise {rise:.3} Pa, rho c u {acoustic:.3} Pa, jet {:.3} m/s, flow {:.4e} kg/s",
        jet.velocity, open.mass_flow_kg_s
    );
    assert!((rise - acoustic).abs() <= 1e-2 * rise);
}

#[test]
fn cup_pressure_above_the_mouth_drives_gas_back_through_the_lip() {
    let gas = specimen();
    let mouth = state(&gas, 101_325.0, 0.0);
    let pushing = state(&gas, 103_325.0, -0.5);
    let back = lip_face(&gas, pushing, mouth, OPENING_M2, INLET_M2).unwrap();
    assert!(back.mass_flow_kg_s < 0.0 && back.pressure_pa > mouth.pressure);
    let wave = gas
        .end_wave(pushing, -1.0, back.pressure_pa, || false)
        .unwrap();
    let donor = Primitive {
        velocity: 0.0,
        ..wave.primitive
    };
    let (jet, _) = gas.nozzle_exit(donor, mouth.pressure, || false).unwrap();
    let metered = -jet.density * jet.velocity * OPENING_M2;
    assert!(
        (back.mass_flow_kg_s - metered).abs() <= 1e-9 * metered.abs(),
        "{} {metered}",
        back.mass_flow_kg_s
    );
    // Leaving gas carries the interior star state's own enthalpy flux.
    let expected = gas.face_flux(wave.primitive, wave.signal_speed).unwrap();
    assert!((back.flux.energy - expected.energy).abs() <= 1e-9 * expected.energy.abs());
}

// A prescribed end volume flow: the interior's own end wave supplies the pressure.
fn flow_port(
    gas: &FrozenGas,
    stage: &PortStage,
    velocity: f64,
    guess: &mut f64,
) -> Result<(FaceFlux, f64), GasError> {
    let inside = stage.left.1;
    let mut pressure = *guess;
    for _ in 0..64 {
        let wave = gas.end_wave(inside, -1.0, pressure, || false)?;
        let next = pressure - (wave.primitive.velocity - velocity) / wave.velocity_per_pressure;
        let done = (next - pressure).abs() <= 1e-13 * pressure;
        pressure = next;
        if done {
            break;
        }
    }
    let wave = gas.end_wave(inside, -1.0, pressure, || false)?;
    *guess = pressure;
    let face = Primitive {
        velocity,
        ..wave.primitive
    };
    Ok((gas.face_flux(face, wave.signal_speed)?, pressure))
}

// |sum x_n e^{-i omega n dt}|.
fn dft_magnitude(signal: &[f64], omega_dt: f64) -> f64 {
    let (s, c) = (det::sin(omega_dt), det::cos(omega_dt));
    let (mut re, mut im, mut pr, mut pi) = (0.0, 0.0, 1.0, 0.0);
    for &x in signal {
        re += x * pr;
        im -= x * pi;
        let next = pr * c - pi * s;
        pi = pr * s + pi * c;
        pr = next;
    }
    re.hypot(im)
}

fn decibels(x: f64) -> f64 {
    20.0 * det::ln(x) / core::f64::consts::LN_10
}

// Peak frequency, level [dB] and half-power quality factor of a resonance of
// `level` near `near` [Hz].
fn resonance(level: &impl Fn(f64) -> f64, near: f64) -> (f64, f64, f64) {
    let cents = |c: f64| near * det::exp(c / 1200.0 * core::f64::consts::LN_2);
    let grid: Vec<f64> = (-300..=300).map(|c| cents(f64::from(c) * 0.5)).collect();
    let values: Vec<f64> = grid.iter().map(|&f| level(f)).collect();
    let best = (1..grid.len() - 1)
        .max_by(|&a, &b| values[a].total_cmp(&values[b]))
        .unwrap();
    let (a, b, c) = (values[best - 1], values[best], values[best + 1]);
    let offset = 0.5 * (a - c) / (a - 2.0 * b + c);
    let peak = grid[best] * det::exp(offset * 0.5 / 1200.0 * core::f64::consts::LN_2);
    let top = level(peak).max(b);
    let half = |outer: f64| {
        let (mut inside, mut outside) = (peak, outer);
        assert!(
            level(outside) < top - 3.0103,
            "resonance near {near} Hz is not isolated"
        );
        for _ in 0..60 {
            let mid = inside.midpoint(outside);
            if level(mid) > top - 3.0103 {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        inside.midpoint(outside)
    };
    let (low, high) = (half(peak * 0.85), half(peak * 1.15));
    (peak, top, peak / (high - low))
}

// Condition 1 of the GO for this bore (mail 9115). Tolerances fixed before
// the first run: each of the first six input-impedance peaks within 10 cents,
// its level within 1.5 dB and its half-power Q within 10 % of the multimodal
// image the linear MmLineBank is realized from (exact Bessel wall loss, the
// same gas, a flanged end). The ends differ in kind and the difference is
// part of what is measured: a rigid baffled piston here, Zorumski's
// multimodal flanged pipe there.
#[test]
fn a_small_amplitude_nonlinear_bore_reproduces_the_linear_banks_input_impedance() {
    let (gas, state) = specimen_with_state();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.006,
                length: 0.30,
            },
            Segment::Cone {
                inlet_radius: 0.006,
                outlet_radius: 0.012,
                length: 0.25,
            },
            Segment::Cone {
                inlet_radius: 0.012,
                outlet_radius: 0.035,
                length: 0.12,
            },
        ],
    };
    let step = 1.0 / 512_000.0;
    let mut bore = NonlinearBore::new(
        &duct,
        gas,
        state.temperature,
        state.pressure,
        Some(WallGas {
            density_kg_m3: state.density,
            dynamic_viscosity_pa_s: state.dynamic_viscosity,
            gamma: state.gamma,
            prandtl: state.prandtl,
        }),
        BoreSpec {
            maximum_cell_m: 2.4e-3,
            inlet_volume_m3: 0.0,
            step_s: step,
            cfl: 0.45,
            wall: Some(WallMemory {
                // 64 poles (the fs-phs cap): the 35 mm end is admitted on
                // 20-8000 Hz at 6e-3; 24 resolve only up to about 6 mm.
                poles: 64,
                band_hz: [20.0, 8000.0],
                maximum_relative_error: 1e-2,
            }),
        },
    )
    .unwrap();
    let room = bore.room();
    assert!((room.sound_speed / state.sound_speed - 1.0).abs() < 1e-6);
    println!("linear limit: wall {:?}", bore.wall_receipt());
    let source = bore.gas().clone();
    let inlet = bore.inlet_area_m2();
    let (width, peak_flow) = (1.0e-4, 1.0e-8);
    let flow = |t: f64| {
        if t < width {
            peak_flow * 0.5 * (1.0 - det::cos(core::f64::consts::TAU * t / width))
        } else {
            0.0
        }
    };
    let samples = (0.3 / step) as usize;
    let (mut pressure, mut volume) = (Vec::with_capacity(samples), Vec::with_capacity(samples));
    let mut guess = room.pressure;
    for n in 0..samples {
        let mut face = room.pressure;
        bore.step(|stage| {
            let (flux, p) = flow_port(&source, stage, flow(stage.time) / inlet, &mut guess)?;
            if stage.stage == 0 {
                face = p;
            }
            Ok(Some(flux))
        })
        .unwrap();
        pressure.push(face - room.pressure);
        volume.push(flow(n as f64 * step));
    }
    let linear = |f: f64| {
        let response = mm_input_impedance(
            &duct,
            &state,
            core::f64::consts::TAU * f,
            LossModel::Bessel,
            Termination::FlangedOpen,
            3,
            1,
        )
        .unwrap();
        decibels(response.plane_impedance.abs())
    };
    let nonlinear = |f: f64| {
        let omega_dt = core::f64::consts::TAU * f * step;
        decibels(dft_magnitude(&pressure, omega_dt) / dft_magnitude(&volume, omega_dt))
    };
    let grid: Vec<f64> = (40..2500).map(f64::from).collect();
    let levels: Vec<f64> = grid.iter().map(|&f| linear(f)).collect();
    let near: Vec<f64> = (1..grid.len() - 1)
        .filter(|&i| levels[i] > levels[i - 1] && levels[i] >= levels[i + 1])
        .map(|i| grid[i])
        .take(6)
        .collect();
    assert_eq!(near.len(), 6, "linear peaks {near:?}");
    for (k, &f) in near.iter().enumerate() {
        let (f_lin, db_lin, q_lin) = resonance(&linear, f);
        let (f_nl, db_nl, q_nl) = resonance(&nonlinear, f_lin);
        let cents = 1200.0 * det::ln(f_nl / f_lin) / core::f64::consts::LN_2;
        println!(
            "linear limit peak {}: bank {f_lin:8.2} Hz {db_lin:6.2} dB Q {q_lin:6.2} | euler {f_nl:8.2} Hz {db_nl:6.2} dB Q {q_nl:6.2} | {cents:+6.2} cents {:+5.2} dB Q {:+5.1} %",
            k + 1,
            db_nl - db_lin,
            100.0 * (q_nl / q_lin - 1.0)
        );
        assert!(cents.abs() <= 10.0);
        assert!((db_nl - db_lin).abs() <= 1.5);
        assert!((q_nl / q_lin - 1.0).abs() <= 0.10);
    }
}

// NOZ#653 of the R5C Bb4 card: the face pressure p sat one ulp above the mouth pressure (a lip just opening at equilibrium) and the primitive recovered
// from the end wave's conserved variables carried a pressure two ulp below p, hence below the mouth pressure. The outflow donor's pressure is p by
// construction, so the metered flow must not depend on how far the recovered pressure fell short of it.
#[test]
fn the_outflow_donor_is_the_wave_pressure_not_its_rounded_recovery() {
    let gas = specimen();
    let mouth = state(&gas, 101_325.0, 0.0);
    let p = mouth.pressure.next_up();
    assert!(p > mouth.pressure);
    let mut star = state(&gas, p, 0.0);
    star.pressure = p.next_down().next_down();
    assert!(star.pressure < mouth.pressure, "the planted star must sit below the mouth pressure");
    let metered_low = metered(&gas, star, mouth, OPENING_M2, p).expect("outflow at equilibrium must not be refused");
    assert!(metered_low.is_finite() && metered_low <= 0.0, "{metered_low}");
    // Independent of the recovery's rounding: the same star with its pressure exactly p, and one ulp above, meter the same flow bit for bit.
    for recovered in [p, p.next_up()] {
        let exact = Primitive { pressure: recovered, ..star };
        let metered_exact = metered(&gas, exact, mouth, OPENING_M2, p).unwrap();
        assert_eq!(metered_low.to_bits(), metered_exact.to_bits());
    }
    // A real outflow, far above the mouth pressure, is still metered and still outward.
    let high = mouth.pressure * 1.02;
    let outward = metered(&gas, state(&gas, high, 0.0), mouth, OPENING_M2, high).unwrap();
    assert!(outward < 0.0, "{outward}");
}

// The planted negative for the test above: the nozzle's precondition is untouched. A backpressure above a stationary donor by a real margin, or a
// donor that is not exactly stationary, is still refused with the same message; equality is still accepted.
#[test]
fn the_nozzle_still_refuses_a_backpressure_above_its_stationary_donor() {
    let gas = specimen();
    let donor = state(&gas, 101_325.0, 0.0);
    assert!(gas.nozzle_exit(donor, donor.pressure, || false).is_ok());
    let above = gas.nozzle_exit(donor, donor.pressure * (1.0 + 1e-6), || false).unwrap_err();
    assert!(format!("{above:?}").contains("stationary nozzle donor or backpressure"), "{above:?}");
    let moving = Primitive { velocity: 1e-9, ..donor };
    let moving = gas.nozzle_exit(moving, donor.pressure * (1.0 - 1e-6), || false).unwrap_err();
    assert!(format!("{moving:?}").contains("stationary nozzle donor or backpressure"), "{moving:?}");
}

// The opt-in lip jet noise on the nonlinear bore. The probe lip (authored, 12 mm wide, 0.5 mm rest gap, 1e-4 m^2 face, tension 1.6)
// on the lip tests' three-segment bore, the specimen gas (mu 1.8e-5 Pa s, air-like), the real cards' clock (32 kHz output, 16
// mechanics per output, one gas step each: 512 kHz) and the wall memory the bore tests use.
fn jet_lip(card: Option<JetNoiseCard>) -> NonlinearLip {
    let mut lip = match card {
        Some(card) => probe_lip().with_jet_noise(card).expect("a usable card"),
        None => probe_lip(),
    };
    lip.apply(BrassControl::SetLipTension(1.6)).expect("tension");
    lip
}

// The same probe lip with the opt-in constriction noise (the cup radius is the card's), or without it.
fn constriction_lip(card: Option<ConstrictionNoiseCard>) -> NonlinearLip {
    let mut lip = match card {
        Some(card) => probe_lip().with_constriction_noise(card).expect("a usable card"),
        None => probe_lip(),
    };
    lip.apply(BrassControl::SetLipTension(1.6)).expect("tension");
    lip
}

fn probe_lip() -> NonlinearLip {
    let (gas, state) = specimen_with_state();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder { radius: 0.006, length: 0.30 },
            Segment::Cone { inlet_radius: 0.006, outlet_radius: 0.012, length: 0.25 },
            Segment::Cone { inlet_radius: 0.012, outlet_radius: 0.035, length: 0.12 },
        ],
    };
    let lip = LipIsland {
        mass_kg: 1.8e-4,
        stiffness_n_m: 350.0,
        damping_n_s_m: 0.05,
        width_m: 0.012,
        rest_gap_m: 5.0e-4,
        face_area_m2: 1.0e-4,
        provenance: "AUTHORED probe values (Estimate)".to_string(),
    };
    let clock = NonlinearClock {
        maximum_cell_m: 2.4e-3,
        output_rate_hz: 32_000,
        mechanics_per_output: 16,
        gas_steps_per_mechanics: 1,
        wall: WallMemory { poles: 64, band_hz: [20.0, 8000.0], maximum_relative_error: 1e-2 },
    };
    let wall_gas = WallGas {
        density_kg_m3: state.density,
        dynamic_viscosity_pa_s: state.dynamic_viscosity,
        gamma: state.gamma,
        prandtl: state.prandtl,
    };
    NonlinearLip::new(&duct, gas, state.temperature, state.pressure, wall_gas, lip, 1.5e-6, clock).expect("nonlinear lip")
}

// Eight blocks of 32 output samples (4096 gas steps) of cup overpressure at a steady `blow` [Pa].
fn jet_run(lip: &mut NonlinearLip, blow: f64) -> Vec<f64> {
    lip.apply(BrassControl::SetBlowingPressure(blow)).expect("blow");
    let (mut out, mut cup, mut flows) = (Vec::new(), vec![0.0f64; 32], vec![0.0f64; 32 * 16]);
    for _ in 0..8 {
        lip.step_block(&mut cup, &mut flows).expect("block");
        out.extend_from_slice(&cup);
    }
    out
}

fn jet_card(onset: f64, seed: u64) -> JetNoiseCard {
    JetNoiseCard { seed, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: onset }
}

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

// Gauntlet G0 and G3. A lip without the jet noise, and one whose gate can never open, are bit-identical; with the gate open (the lip at
// its rest gap passes Re_h about 2300 at 3 kPa) the cup pressure changes, stays finite, and is seeded.
#[test]
fn lip_jet_noise_on_the_nonlinear_bore_shut_gate_is_bit_identical_and_an_open_gate_changes_the_pressure_seeded() {
    let plain = jet_run(&mut jet_lip(None), 3000.0);
    assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
    let shut = jet_run(&mut jet_lip(Some(jet_card(1.0e12, 7))), 3000.0);
    assert!(same_bits(&plain, &shut), "a gate that never opens adds exactly zero");
    let open = jet_run(&mut jet_lip(Some(jet_card(1700.0, 7))), 3000.0);
    assert!(open.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&plain, &open), "an open gate must reach the cup pressure");
    assert!(same_bits(&open, &jet_run(&mut jet_lip(Some(jet_card(1700.0, 7))), 3000.0)), "a declared seed repeats bit for bit");
    assert!(!same_bits(&open, &jet_run(&mut jet_lip(Some(jet_card(1700.0, 8))), 3000.0)), "another seed gives another noise");
}

// Gauntlet G3 with a planted negative. At 300 Pa the jet is about 22 m/s and Re_h = mdot / (w mu) stays far under 1700: the row with the
// 1700 gate is bit-identical to no row, while the same row with NO gate speaks (so the gate is what shut the first). With no mouth pressure
// there is no forward flow and no source, whatever the gate; reverse flow is silent at the unit level (the source needs q > 0).
#[test]
fn lip_jet_noise_on_the_nonlinear_bore_is_gated_on_the_jet_and_silent_without_forward_flow() {
    let plain = jet_run(&mut jet_lip(None), 300.0);
    let gated = jet_run(&mut jet_lip(Some(jet_card(1700.0, 7))), 300.0);
    assert!(same_bits(&plain, &gated), "below the gate the source stays shut");
    let ungated = jet_run(&mut jet_lip(Some(jet_card(0.0, 7))), 300.0);
    assert!(!same_bits(&plain, &ungated), "an ungated source speaks at the same pressure, so the gate is what shut the first");
    let still = jet_run(&mut jet_lip(None), 0.0);
    let none = jet_run(&mut jet_lip(Some(jet_card(0.0, 7))), 0.0);
    assert!(same_bits(&still, &none), "no mouth pressure, no forward flow: no source even with no gate");
}

// The probe bore's first station is a 6 mm radius cylinder (`probe_lip`): that is the cup-rim radius these tests use, the way a real
// card's `constriction_noise` radius is the radius of its own first bore station.
const PROBE_CUP_RADIUS_M: f64 = 6.0e-3;

fn constriction_card(seed: u64, duct_radius_m: f64) -> ConstrictionNoiseCard {
    ConstrictionNoiseCard { seed, duct_radius_m }
}

// Gauntlet G0 and G3 for the constriction noise on the nonlinear bore. A lip without it is the lip built before it existed; with no mouth
// pressure there is no forward flow and no source; the planted negative (a 12 m "duct", whose cutoff 100 / r = 8.3 Hz is under the bank's
// lowest band, so nothing is driven at any pressure) is bit-identical to no row at the pressure where the real radius speaks; the real
// radius changes the cup pressure, keeps it finite, repeats bit for bit under a declared seed, and another seed gives another noise.
#[test]
fn lip_constriction_noise_on_the_nonlinear_bore_is_absent_silent_without_a_jet_and_speaks_seeded_on_the_printed_law() {
    let plain = jet_run(&mut constriction_lip(None), 3000.0);
    assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
    let still = jet_run(&mut constriction_lip(None), 0.0);
    let none = jet_run(&mut constriction_lip(Some(constriction_card(7, PROBE_CUP_RADIUS_M))), 0.0);
    assert!(same_bits(&still, &none), "no mouth pressure, no forward flow: no source");
    let cut_off = jet_run(&mut constriction_lip(Some(constriction_card(7, 12.0))), 3000.0);
    assert!(same_bits(&plain, &cut_off), "no band below the duct cutoff, so the law supplies nothing");
    let open = jet_run(&mut constriction_lip(Some(constriction_card(7, PROBE_CUP_RADIUS_M))), 3000.0);
    assert!(open.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&plain, &open), "the printed law drives bands at this pressure, so it must reach the cup pressure");
    assert!(
        same_bits(&open, &jet_run(&mut constriction_lip(Some(constriction_card(7, PROBE_CUP_RADIUS_M))), 3000.0)),
        "a declared seed repeats bit for bit"
    );
    assert!(
        !same_bits(&open, &jet_run(&mut constriction_lip(Some(constriction_card(8, PROBE_CUP_RADIUS_M))), 3000.0)),
        "another seed gives another noise"
    );
}

// Gauntlet G3. The two opt-in sources share the lip's series drive and the mouth reservoir without displacing each other: the pair gives
// a cup pressure different from either alone, finite, and repeatable. The jet noise alone is the jet noise of the lip built before the
// constriction noise existed (the refactored drive adds the same terms in the same order when only one is present).
#[test]
fn lip_constriction_noise_and_jet_noise_on_the_nonlinear_bore_add_in_the_same_series_drive() {
    let both = |seed: u64| {
        let mut lip = probe_lip()
            .with_jet_noise(jet_card(1700.0, 7))
            .expect("jet card")
            .with_constriction_noise(constriction_card(seed, PROBE_CUP_RADIUS_M))
            .expect("constriction card");
        lip.apply(BrassControl::SetLipTension(1.6)).expect("tension");
        jet_run(&mut lip, 3000.0)
    };
    let jet_only = jet_run(&mut jet_lip(Some(jet_card(1700.0, 7))), 3000.0);
    let constriction_only = jet_run(&mut constriction_lip(Some(constriction_card(7, PROBE_CUP_RADIUS_M))), 3000.0);
    let pair = both(7);
    assert!(pair.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&pair, &jet_only), "the constriction noise reaches the cup beside the jet noise");
    assert!(!same_bits(&pair, &constriction_only), "the jet noise reaches the cup beside the constriction noise");
    assert!(same_bits(&pair, &both(7)), "the pair repeats bit for bit");
    assert!(!same_bits(&pair, &both(8)), "another constriction seed gives another pair");
}
