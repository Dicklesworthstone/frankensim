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
