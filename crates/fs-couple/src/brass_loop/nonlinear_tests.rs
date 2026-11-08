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

pub(super) fn specimen() -> FrozenGas {
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
    let wave = gas
        .end_wave(interior, -1.0, open.pressure_pa, || false)
        .unwrap();
    let rise = open.pressure_pa - interior.pressure;
    let acoustic = interior.density * interior.sound_speed * wave.primitive.velocity;
    println!(
        "open lip: face rise {rise:.3} Pa, rho c u {acoustic:.3} Pa, jet {:.3} m/s, flow {:.4e} kg/s",
        jet.velocity, open.mass_flow_kg_s
    );
    assert!((rise - acoustic).abs() <= 1e-2 * rise);
}

#[test]
fn an_open_lip_closes_static_pressure_mass_and_total_enthalpy_together() {
    let gas = specimen();
    let interior = state(&gas, 101_325.0, 0.0);
    // This specimen's Cp/R=3.5+0.001*T is known independently of the
    // production caloric inversion. Include substantial pressure and area
    // ratios so retaining the mouth temperature cannot pass the energy check.
    for (mouth_pressure, opening) in [
        (104_325.0, OPENING_M2),
        (136_325.0, 0.75 * INLET_M2),
        (202_650.0, 0.5 * INLET_M2),
    ] {
        let mouth = state(&gas, mouth_pressure, 0.0);
        let open = lip_face(&gas, interior, mouth, opening, INLET_M2).unwrap();
        let wave = gas
            .end_wave(interior, -1.0, open.pressure_pa, || false)
            .unwrap();
        let velocity = wave.primitive.velocity;
        assert!(velocity > 0.0);
        let density = open.flux.mass / velocity;
        let temperature = open.pressure_pa / (gas.gas_constant() * density);
        let cooling = mouth.temperature - temperature;
        let kinetic = 0.5 * velocity * velocity;
        let enthalpy_drop =
            gas.gas_constant() * cooling * (3.5 + 0.001 * mouth.temperature - 0.0005 * cooling);
        assert!(temperature >= gas.domain().minimum_temperature && cooling > 0.0);
        assert!(
            (enthalpy_drop - kinetic).abs() <= 2e-8 + 2e-9 * kinetic,
            "p0={mouth_pressure}, kinetic={kinetic}, caloric drop={enthalpy_drop}"
        );
        assert!(
            (open.flux.momentum - (open.flux.mass * velocity + open.pressure_pa)).abs()
                <= 1e-10 * open.pressure_pa
        );
        let (jet, _) = gas.nozzle_exit(mouth, open.pressure_pa, || false).unwrap();
        let metered = jet.density * jet.velocity * opening;
        assert!((open.mass_flow_kg_s / metered - 1.0).abs() <= 2e-9);
        // Adiabatic mixing converts the nozzle jet's excess kinetic energy
        // back into heat; the slower face still lies below stagnation T.
        assert!(velocity < jet.velocity && temperature > jet.temperature);
        assert_eq!(
            open.flux.energy,
            open.flux.mass * gas.stagnation_enthalpy(mouth.temperature).unwrap()
        );
        let face_energy =
            open.flux.mass * (gas.stagnation_enthalpy(temperature).unwrap() + kinetic);
        assert!((open.flux.energy - face_energy).abs() <= 2e-9 * open.flux.energy.abs());
        let cp_r = 3.5 + 0.001 * temperature;
        let sound_speed = (cp_r / (cp_r - 1.0) * gas.gas_constant() * temperature).sqrt();
        assert!(open.flux.signal_speed >= velocity + sound_speed - 1e-10 * sound_speed);
    }
}

#[test]
fn a_cold_auxiliary_pressure_does_not_refuse_supported_choked_lip_inflow() {
    let gas = specimen();
    let interior = state(&gas, 101_325.0, 0.0);
    let mouth = state(&gas, 6.0 * interior.pressure, 0.0);
    let upper = gas
        .end_wave(interior, -1.0, mouth.pressure, || false)
        .unwrap();
    assert!(matches!(
        gas.stagnation_to_static_state(
            mouth.temperature,
            mouth.pressure,
            upper.primitive.velocity,
            || false
        ),
        Err(GasError::CaloricPath(_))
    ));
    let open = lip_face(&gas, interior, mouth, OPENING_M2, INLET_M2).unwrap();
    assert!(open.pressure_pa > interior.pressure && open.pressure_pa < mouth.pressure);
    let (jet, choked) = gas.nozzle_exit(mouth, open.pressure_pa, || false).unwrap();
    assert!(choked);
    let metered = jet.density * jet.velocity * OPENING_M2;
    assert!((open.mass_flow_kg_s / metered - 1.0).abs() < 2e-9);
    let wave = gas
        .end_wave(interior, -1.0, open.pressure_pa, || false)
        .unwrap();
    let density = open.flux.mass / wave.primitive.velocity;
    let temperature = open.pressure_pa / (gas.gas_constant() * density);
    assert!(temperature >= gas.domain().minimum_temperature && temperature < mouth.temperature);
    let total_enthalpy = gas.stagnation_enthalpy(temperature).unwrap()
        + 0.5 * wave.primitive.velocity * wave.primitive.velocity;
    assert!((open.flux.energy / open.flux.mass - total_enthalpy).abs() < 1e-8 * total_enthalpy);
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

#[test]
fn nonlinear_lip_tension_work_is_signed_atomic_and_separate_from_gas_exchange() {
    let mut lip = probe_lip();
    lip.apply(BrassControl::SetLipTension(1.5)).unwrap();
    assert_eq!(lip.lip_tension_work_j(), [0.0; 2]);
    lip.lip_state = vec![2e-4, 1e-5];
    let retained = lip.lip_state.clone();
    let gas_before = lip.bore.solver().states().to_vec();
    let time_before = lip.bore.solver().time();
    let energy_before = lip.lip_sys.hamiltonian(&retained);
    // Two opposite actuator transfers before any gas or mechanical advancement.
    lip.apply(BrassControl::SetLipTension(1.75)).unwrap();
    assert!(lip.lip_tension_work_j()[0] > 0.0);
    lip.apply(BrassControl::SetLipTension(1.25)).unwrap();
    let expected = 0.5
        * lip.lip.stiffness_n_m
        * (1.25_f64.powi(2) - 1.5_f64.powi(2))
        * retained[0]
        * retained[0];
    let work = lip.lip_tension_work_j();
    assert!(work[0] < 0.0);
    assert!((work[0] - expected).abs() < 2e-14 * expected.abs());
    assert!(
        (lip.lip_sys.hamiltonian(&retained) - energy_before - expected).abs()
            < 2e-14 * expected.abs()
    );
    assert_eq!(work[1], 0.0);
    assert_eq!(lip.lip_state, retained);
    assert_eq!(lip.bore.solver().states(), gas_before);
    assert_eq!(lip.bore.solver().time(), time_before);
    let account = lip.tension_work;
    let energy = lip.lip_sys.hamiltonian(&retained);
    for invalid in [f64::NAN, 0.0, -1.0, f64::INFINITY, 1e300] {
        assert!(lip.apply(BrassControl::SetLipTension(invalid)).is_err());
        assert_eq!(lip.tension, 1.25);
        assert_eq!(lip.lip_state, retained);
        assert_eq!(lip.lip_sys.hamiltonian(&retained), energy);
        assert_eq!(lip.tension_work, account);
    }
    lip.apply(BrassControl::SetLipTension(1.25)).unwrap();
    assert_eq!(lip.tension_work, account);
    assert!(lip.step_block(&mut [0.0], &mut []).is_err());
    assert_eq!(lip.tension_work, account);
    let mut pressure = [0.0];
    let mut flow = vec![0.0; lip.clock.mechanics_per_output as usize];
    lip.step_block(&mut pressure, &mut flow).unwrap();
    let diag = lip.diagnostics().last().unwrap();
    assert_eq!(diag.lip_tension_work_j, work);
    assert_eq!(diag.blow_work_j, 0.0);
    lip.step_block(&mut pressure, &mut flow).unwrap();
    assert_eq!(
        lip.diagnostics().last().unwrap().lip_tension_work_j,
        [0.0; 2]
    );
    assert_eq!(lip.lip_tension_work_j(), work);
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


fn finite_spec(volume_m3: f64, feed_area_m2: f64) -> FiniteMouthSpec {
    FiniteMouthSpec {
        volume_m3,
        feed_area_m2,
        initial_pressure_pa: 106_325.0,
        initial_temperature_k: 293.15,
    }
}

#[test]
fn finite_mouth_pressure_drives_the_actual_lip_without_a_control_reset() {
    let mut lip = probe_lip().with_finite_mouth(finite_spec(2e-5, 0.0)).unwrap();
    let before = lip.finite_mouth_state().unwrap();
    // A changed gesture drives the upstream feed. It cannot create or erase
    // stored gas, including when that feed is sealed.
    lip.apply(BrassControl::SetBlowingPressure(35_000.0)).unwrap();
    assert_eq!(lip.finite_mouth_state().unwrap(), before);
    lip.apply(BrassControl::SetBlowingPressure(0.0)).unwrap();
    assert_eq!(lip.finite_mouth_state().unwrap(), before);
    let dt = lip.clock.gas_step_s();
    let (_, _, face, _, drive, _) = lip.gas_step(dt).unwrap();
    assert!((drive - 5000.0).abs() < 1e-8, "actual stored overpressure drives the lip");
    assert!(face.mass_flow_kg_s > 0.0 && face.pressure_pa > lip.room.pressure);
    let report = *lip.finite_mouth_report().unwrap();
    assert_eq!(report.before, before);
    assert_eq!(report.after, lip.finite_mouth_state().unwrap());
    assert_eq!(report.feed_mass_kg, 0.0);
    assert_eq!(report.feed_energy_j, 0.0);
    assert!(report.bore_mass_kg > 0.0);
    assert!(report.after.pressure_pa < before.pressure_pa);
    assert!(report.after.temperature_k < before.temperature_k);
    assert!(report.mass_residual_kg.abs() < 1e-12 * before.mass_kg);
    assert!(report.energy_residual_j.abs() < 1e-12 * before.energy_j.abs());
    for _ in 0..31 {
        lip.gas_step(dt).unwrap();
    }
    assert!(lip.finite_mouth_state().unwrap().pressure_pa < report.after.pressure_pa);
}

#[test]
fn refused_finite_mouth_stage_keeps_the_lip_and_storage_at_the_last_accepted_state() {
    let mut lip = probe_lip().with_finite_mouth(finite_spec(1e-12, 0.0)).unwrap();
    let before = lip.finite_mouth_state().unwrap();
    let lip_state = lip.lip_state.clone();
    let cells = lip.bore.solver().states().to_vec();
    let clock = lip.bore.solver().time();
    let face = lip.face_pa;
    let dt = lip.clock.gas_step_s();
    let refusal = lip.gas_step(dt).unwrap_err();
    assert!(refusal.contains("CflExceeded"));
    assert_eq!(lip.finite_mouth_state().unwrap(), before);
    assert!(lip.finite_mouth_report().is_none());
    assert_eq!(lip.lip_state, lip_state);
    assert_eq!(lip.face_pa, face);
    assert_eq!(lip.bore.solver().states(), cells);
    assert_eq!(lip.bore.solver().time(), clock);
}

#[test]
fn finite_mouth_is_attached_once_before_stepping_and_admits_noise_rows_in_either_order() {
    let spec = finite_spec(2e-5, 1e-6);
    assert!(probe_lip().with_finite_mouth(spec).unwrap().with_finite_mouth(spec).is_err());
    let mut advanced = probe_lip();
    advanced.apply(BrassControl::SetBlowingPressure(3000.0)).unwrap();
    let dt = advanced.clock.gas_step_s();
    advanced.gas_step(dt).unwrap();
    assert!(advanced.with_finite_mouth(spec).is_err());
    assert!(probe_lip().with_finite_mouth(spec).unwrap().with_jet_noise(jet_card(0.0, 7)).is_ok());
    assert!(probe_lip().with_jet_noise(jet_card(0.0, 7)).unwrap().with_finite_mouth(spec).is_ok());
    assert!(probe_lip().with_finite_mouth(spec).unwrap()
        .with_constriction_noise(constriction_card(7, 0.006)).is_ok());
    assert!(probe_lip().with_constriction_noise(constriction_card(7, 0.006)).unwrap()
        .with_finite_mouth(spec).is_ok());
}

// The finite mouth with the opt-in noise rows. `noise_first` puts the noise rows on the lip before the finite mouth is attached.
fn finite_noise_lip(spec: FiniteMouthSpec, jet: Option<JetNoiseCard>, constriction: Option<ConstrictionNoiseCard>, noise_first: bool) -> NonlinearLip {
    let attach = |mut lip: NonlinearLip| {
        if let Some(card) = jet {
            lip = lip.with_jet_noise(card).expect("jet card");
        }
        if let Some(card) = constriction {
            lip = lip.with_constriction_noise(card).expect("constriction card");
        }
        lip
    };
    let mut lip = if noise_first {
        attach(probe_lip()).with_finite_mouth(spec).expect("finite mouth")
    } else {
        attach(probe_lip().with_finite_mouth(spec).expect("finite mouth"))
    };
    lip.apply(BrassControl::SetLipTension(1.6)).expect("tension");
    lip
}

// The closure's order test: the noise rows and the finite mouth give the same lip whichever is attached first, for the jet noise, the
// constriction noise and both; the noise reaches the cup pressure (it differs from the same finite mouth without it) and stays finite.
#[test]
fn finite_mouth_with_series_noise_is_the_same_lip_in_either_builder_order_and_the_noise_reaches_the_cup() {
    let spec = finite_spec(2e-5, 1e-5);
    let plain = jet_run(&mut finite_noise_lip(spec, None, None, false), 3000.0);
    assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
    let cases: [(Option<JetNoiseCard>, Option<ConstrictionNoiseCard>); 3] = [
        (Some(jet_card(1700.0, 7)), None),
        (None, Some(constriction_card(7, PROBE_CUP_RADIUS_M))),
        (Some(jet_card(1700.0, 7)), Some(constriction_card(7, PROBE_CUP_RADIUS_M))),
    ];
    for (jet, constriction) in cases {
        let after = jet_run(&mut finite_noise_lip(spec, jet, constriction, false), 3000.0);
        let before = jet_run(&mut finite_noise_lip(spec, jet, constriction, true), 3000.0);
        assert!(after.iter().all(|p| p.is_finite()));
        assert!(same_bits(&after, &before), "either builder order gives the same lip");
        assert!(!same_bits(&plain, &after), "the noise must reach the cup pressure through the finite mouth");
    }
}

// Gates that never open add exactly zero: the finite mouth with a shut jet gate, with a constriction cutoff above the band, or both
// is bit-identical to the finite mouth without a noise row (the offset is exactly 0.0, and the step is the one built before the offset).
#[test]
fn finite_mouth_with_silent_noise_rows_is_bit_identical_to_the_finite_mouth_alone() {
    let spec = finite_spec(2e-5, 1e-5);
    let plain = jet_run(&mut finite_noise_lip(spec, None, None, false), 3000.0);
    for (jet, constriction) in [
        (Some(jet_card(1.0e12, 7)), None),
        (None, Some(constriction_card(7, 12.0))),
        (Some(jet_card(1.0e12, 7)), Some(constriction_card(7, 12.0))),
    ] {
        let silent = jet_run(&mut finite_noise_lip(spec, jet, constriction, false), 3000.0);
        assert!(same_bits(&plain, &silent), "a silent noise row adds exactly zero to the finite mouth");
    }
}

// The large-volume limit, the physics check of the closure. A 1 m3 plenum that starts at the feed pressure and the room temperature does
// not move during the 8 ms of this run, so it is the original held reservoir: the finite mouth with the noise must give the cup pressure
// of the ideal mouth with the same noise (same seeds, same blow) to a part in 1e4 of its rms. The check can only see a leak of the noise into
// the wrong place if the noise is itself well above that tolerance: its precondition is ten times the tolerance, 1e-3 of the rms of the signal.
// Every case prints its measured ratios before any assertion, so a miss is read as a number and not as a bare failure.
#[test]
fn finite_mouth_with_series_noise_in_the_large_volume_limit_is_the_ideal_mouth_with_the_same_noise() {
    let room = probe_lip().room;
    let blow = 3000.0;
    let spec = FiniteMouthSpec {
        volume_m3: 1.0,
        feed_area_m2: 1e-2,
        initial_pressure_pa: room.pressure + blow,
        initial_temperature_k: room.temperature,
    };
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
    let mut report = Vec::new();
    for (name, jet, constriction) in [
        ("jet", Some(jet_card(1700.0, 7)), None),
        ("constriction", None, Some(constriction_card(7, PROBE_CUP_RADIUS_M))),
        ("both", Some(jet_card(1700.0, 7)), Some(constriction_card(7, PROBE_CUP_RADIUS_M))),
    ] {
        let mut ideal = probe_lip();
        if let Some(card) = jet {
            ideal = ideal.with_jet_noise(card).expect("jet card");
        }
        if let Some(card) = constriction {
            ideal = ideal.with_constriction_noise(card).expect("constriction card");
        }
        ideal.apply(BrassControl::SetLipTension(1.6)).expect("tension");
        let ideal = jet_run(&mut ideal, blow);
        let finite = jet_run(&mut finite_noise_lip(spec, jet, constriction, false), blow);
        let quiet = jet_run(&mut jet_lip(None), blow);
        let difference: Vec<f64> = ideal.iter().zip(&finite).map(|(a, b)| a - b).collect();
        let noise: Vec<f64> = ideal.iter().zip(&quiet).map(|(a, b)| a - b).collect();
        report.push((name, rms(&noise) / rms(&ideal), rms(&difference) / rms(&ideal), rms(&ideal)));
    }
    for (name, noise_ratio, difference_ratio, ideal_rms) in &report {
        eprintln!("T8 {name}: noise / ideal rms {noise_ratio:.3e}, (finite - ideal) / ideal rms {difference_ratio:.3e}, ideal rms {ideal_rms:.3} Pa");
    }
    for (name, noise_ratio, difference_ratio, _) in &report {
        assert!(*noise_ratio > 1e-3, "{name}: the noise is only {noise_ratio:.3e} of the signal, too small for the limit to mean anything");
        assert!(*difference_ratio <= 1e-4, "{name}: finite mouth differs from the ideal mouth by {difference_ratio:.3e} of the rms");
    }
}

// MB 3915: the linear loop's two opt-in lip rows on the nonlinear bore. The probe lip at its default tension 1.0 (the only tension the
// second mode runs at) with the second mode (Cullen 2000 Table A.2 embouchure 7: 244.5 Hz, Q 11, 1/mu -0.073 m^2/kg, on the probe lip's
// 1e-4 m^2 face), the flow-area exponent, or neither.
fn rows_lip(q: Option<f64>, mode2: Option<LipMode2>) -> NonlinearLip {
    let mut lip = probe_lip();
    if let Some(mode) = mode2 {
        lip = lip.with_lip_mode2(mode).expect("a usable second mode");
    }
    if let Some(q) = q {
        lip = lip.with_lip_area_exponent(q).expect("a usable exponent");
    }
    lip
}

fn embouchure_7_mode2() -> LipMode2 {
    LipMode2::from_response(244.5, 11.0, -0.073, 1.0e-4).expect("embouchure 7")
}

fn rms(x: &[f64]) -> f64 {
    (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt()
}

// G0 and G3. q = 1 (the row set to exactly 1) steps exactly the plain lip, bit for bit, on the nonlinear bore; q = 1.4 reaches the cup
// pressure by more than 1 % of its rms, stays finite and repeats bit for bit.
#[test]
fn lip_area_exponent_one_on_the_nonlinear_bore_is_bit_identical_and_q_reaches_the_cup_pressure() {
    let plain = jet_run(&mut rows_lip(None, None), 3000.0);
    assert!(plain.iter().any(|p| p.abs() > 1.0), "the probe cup must carry pressure for the comparison to mean anything");
    let unit = jet_run(&mut rows_lip(Some(1.0), None), 3000.0);
    assert!(same_bits(&plain, &unit), "q = 1 must be the plain w h on the nonlinear bore");
    let steeper = jet_run(&mut rows_lip(Some(1.4), None), 3000.0);
    assert!(steeper.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&plain, &steeper), "q = 1.4 must reach the cup pressure");
    assert!((rms(&steeper) - rms(&plain)).abs() > 0.01 * rms(&plain), "{} against {}", rms(&steeper), rms(&plain));
    assert!(same_bits(&steeper, &jet_run(&mut rows_lip(Some(1.4), None), 3000.0)), "the same card repeats bit for bit");
}

// The aperture is the linear loop's: `w h (h/H0)^(q-1)`, the rest gap's aperture unchanged, zero when closed, and `w h` bit for bit
// without the row; the exponent is refused by name outside [1, 4].
#[test]
fn lip_area_exponent_on_the_nonlinear_bore_scales_the_aperture_about_the_rest_gap_as_the_linear_loop_does() {
    let lip = rows_lip(Some(2.0), None);
    let (w, h0) = (lip.lip.width_m, lip.lip.rest_gap_m);
    assert_eq!(lip.flow_width_height(h0).to_bits(), (w * h0).to_bits());
    let wide = 2.0 * h0;
    assert!((lip.flow_width_height(wide) - w * wide * 2.0).abs() <= 1e-15 * w * wide * 2.0);
    let narrow = 0.5 * h0;
    assert!((lip.flow_width_height(narrow) - w * narrow * 0.5).abs() <= 1e-15 * w * narrow);
    assert_eq!(lip.flow_width_height(0.0), 0.0);
    let q14 = rows_lip(Some(1.4), None);
    let plain = rows_lip(None, None);
    for gap in [0.3 * h0, h0, 1.7 * h0, 3.0 * h0] {
        let expected = w * gap * (gap / h0).powf(1.4_f64 - 1.0);
        assert_eq!(q14.flow_width_height(gap).to_bits(), expected.to_bits(), "gap {gap}");
        assert_eq!(plain.flow_width_height(gap).to_bits(), (w * gap).to_bits(), "no row: the plain w h at gap {gap}");
    }
    for q in [0.0, 0.99, 4.01, -1.0, f64::NAN, f64::INFINITY] {
        let refusal = probe_lip().with_lip_area_exponent(q).err().expect("outside [1, 4]");
        assert!(refusal.contains("lip_area_exponent"), "{q}: {refusal}");
    }
    assert!(probe_lip().with_lip_area_exponent(4.0).is_ok());
}

// G3. The second lip mode reaches the cup pressure by more than 1 % of its rms, stays finite, repeats bit for bit and commits its own
// state; both rows together differ from each alone.
#[test]
fn the_second_lip_mode_on_the_nonlinear_bore_reaches_the_cup_pressure_and_commits_its_state() {
    let plain = jet_run(&mut rows_lip(None, None), 3000.0);
    let two = jet_run(&mut rows_lip(None, Some(embouchure_7_mode2())), 3000.0);
    assert!(two.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&plain, &two), "a second mode must reach the cup pressure");
    assert!((rms(&two) - rms(&plain)).abs() > 0.01 * rms(&plain), "{} against {}", rms(&two), rms(&plain));
    assert!(same_bits(&two, &jet_run(&mut rows_lip(None, Some(embouchure_7_mode2())), 3000.0)), "the same card repeats bit for bit");
    let mut lip = rows_lip(None, Some(embouchure_7_mode2()));
    jet_run(&mut lip, 3000.0);
    assert!(lip.lip2.as_ref().is_some_and(|second| second.state.iter().any(|x| *x != 0.0)), "the second mode's state is committed");
    let both = jet_run(&mut rows_lip(Some(1.4), Some(embouchure_7_mode2())), 3000.0);
    assert!(both.iter().all(|p| p.is_finite()));
    assert!(!same_bits(&both, &two) && !same_bits(&both, &jet_run(&mut rows_lip(Some(1.4), None), 3000.0)));
}

// The second mode has no printed pre-stress law: it refuses any tension but 1.0 by name, whether the tension comes after it or before.
#[test]
fn the_second_lip_mode_on_the_nonlinear_bore_runs_only_at_unit_tension_and_names_the_refusal() {
    let mut lip = rows_lip(None, Some(embouchure_7_mode2()));
    lip.apply(BrassControl::SetLipTension(1.0)).expect("tension 1 is the printed state");
    let refusal = lip.apply(BrassControl::SetLipTension(1.2)).expect_err("another tension has no printed law");
    assert!(refusal.contains("lip_response_mode2"), "{refusal}");
    let mut one_mode = rows_lip(Some(1.4), None);
    one_mode.apply(BrassControl::SetLipTension(1.2)).expect("the one-mode lip scales as before");
    let mut tensioned = probe_lip();
    tensioned.apply(BrassControl::SetLipTension(1.3)).expect("tension");
    let refusal = tensioned.with_lip_mode2(embouchure_7_mode2()).err().expect("added after a tension change");
    assert!(refusal.contains("lip_response_mode2"), "{refusal}");
}
