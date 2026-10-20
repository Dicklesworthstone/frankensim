//! G1/G3 gas swaps: independent ideal-gas/monopole oracles and actual live PCM.
//! These qualify the reduced model, not an experimentally measured instrument.

use super::*;

fn gas(record: &str) -> GasState {
    let mut fields = record.split_ascii_whitespace();
    let state = fs_material::gas_species::input::read_gas_state(&mut fields, None).unwrap();
    assert!(fields.next().is_none());
    state
}

fn near(actual: f64, expected: f64) {
    let tolerance = 2e-11 * expected.abs().max(1e-16);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

// A single rigid piston coordinate makes the mechanical oracle analytic and
// independent of grid eigenvalue estimation or the instrument's fitted losses.
fn piston(gas: GasState) -> Instrument {
    let mut ins = Instrument::with_gas(96_000.0, gas).unwrap();
    let mut body = Body::point("test piston", 0.002);
    body.h = 0.02;
    body.area = body.h * body.h;
    body.radiates = true;
    body.w[0] = 1e-4;
    body.wp[0] = body.w[0];
    body.sync();
    let head = ins.add(body);
    ins.cavities.push(Cavity {
        heads: vec![head],
        volume: 0.001,
        port: None,
    });
    ins.mics.push(Mic::new([0.0, 0.0, 1.0]));
    ins
}

#[test]
fn pressure_temperature_and_composition_reach_force_delay_and_pressure() {
    let air = gas("300 100000 0");
    let dense = gas("300 200000 0");
    let hot = gas("320 100000 0");
    let carbon_dioxide = gas("mixture 300 100000 mole 1 CO2 1");
    near(dense.density, 2.0 * air.density);
    near(dense.sound_speed, air.sound_speed);
    assert!(hot.sound_speed > air.sound_speed);
    assert!(carbon_dioxide.sound_speed < air.sound_speed);

    for medium in [air, dense, hot, carbon_dioxide] {
        let mut ins = piston(medium);
        let image = ins.image();
        assert_eq!(image.gas, medium);
        // Uniform adiabatic gas: F = -(gamma p A^2/V) x.
        let area = 0.02_f64.powi(2);
        let stiffness = medium.gamma * medium.pressure * area.powi(2) / 0.001;
        near(ins.admission_rate2().0, stiffness / 0.002);
        let tap = image.mics[0].taps[0];
        near(
            tap.2 as f64 + tap.4 / (tap.3 + tap.4),
            ins.rate / medium.sound_speed,
        );
        near(
            tap.3 + tap.4,
            medium.density * area * ins.rate.powi(2) / (2.0 * PI),
        );
        let dt = ins.k();
        ins.step(&mut [0.0]);
        near(
            ins.bodies[0].w[0],
            1e-4 * (1.0 - dt * dt * stiffness / 0.002),
        );
    }
}

#[test]
fn declared_medium_reaches_drum_added_mass_port_and_radiation() {
    let air = gas("300 100000 0");
    let dense = gas("300 200000 0");
    let carbon_dioxide = gas("mixture 300 100000 mole 1 CO2 1");
    let mut spec = kit::kick();
    spec.grid = 0.035;
    spec.oversample = 8;
    spec.port = Some(0.025);
    let mut a = kit::drum_with_gas(&spec, &[], air).unwrap().instrument;
    let mut b = kit::drum_with_gas(&spec, &[], dense).unwrap().instrument;
    let mut c = kit::drum_with_gas(&spec, &[], carbon_dioxide)
        .unwrap()
        .instrument;
    assert_eq!(a.bodies[0].mass, b.bodies[0].mass);
    assert_eq!(a.bodies[0].stiff.val, c.bodies[0].stiff.val);
    near(
        a.bodies[0].air_mass,
        8.0 * air.density * spec.radius.powi(3) / 3.0,
    );
    near(b.bodies[0].air_mass, 2.0 * a.bodies[0].air_mass);
    near(b.rads[0].resistance, 2.0 * a.rads[0].resistance);
    near(
        b.cavities[0].port.as_ref().unwrap().resistance,
        2.0 * a.cavities[0].port.as_ref().unwrap().resistance,
    );
    near(
        c.rads[0].resistance / a.rads[0].resistance,
        carbon_dioxide.density * air.sound_speed / (air.density * carbon_dioxide.sound_speed),
    );
    assert!(b.admission_rate2().0 > a.admission_rate2().0);
    assert_eq!(a.image().gas, air);
    assert_eq!(b.image().gas, dense);
    assert_eq!(c.image().gas, carbon_dioxide);
}

#[test]
fn live_composition_changes_motion_and_pcm_and_survives_clone_and_blocks() {
    let mut air = live::Live::with_gas_record("tom-12", 8.0, "300 100000 0").unwrap();
    let mut co2 =
        live::Live::with_gas_record("tom-12", 8.0, "mixture 300 100000 mole 1 CO2 1").unwrap();
    air.strike(0.25 * air.radius, 0.0, 2.0);
    co2.strike(0.25 * co2.radius, 0.0, 2.0);
    let mut split = co2.clone();
    let (mut a, mut b, mut blocks) = ([0.0; 1024], [0.0; 1024], [0.0; 1024]);
    air.render(&mut a, 1.0);
    co2.render(&mut b, 1.0);
    split.render(&mut blocks[..319], 1.0);
    split.render(&mut blocks[319..], 1.0);
    assert!(a.iter().chain(&b).all(|x| x.is_finite()));
    assert!(a.iter().any(|x| *x != 0.0));
    assert_ne!(a, b, "a declared gas swap must reach rendered sound");
    assert_eq!(
        b, blocks,
        "the retained medium must survive clones and block boundaries"
    );
    assert_ne!(air.ins.bodies[air.target].w, co2.ins.bodies[co2.target].w);
    assert_eq!(co2.ins.gas(), split.ins.gas());
}

#[test]
fn incoherent_and_incomplete_gas_declarations_refuse_before_construction() {
    let mut broken = standard_gas();
    broken.density *= 2.0;
    assert!(Instrument::with_gas(48_000.0, broken).is_err());
    assert!(live::Live::with_gas("tom-12", 8.0, broken).is_err());
    for record in [
        "300 0 0",
        "300 100000 0 trailing",
        "mixture 300 100000 mole 2 N2 0.2 O2 0.3",
    ] {
        assert!(
            live::Live::with_gas_record("tom-12", 8.0, record).is_err(),
            "{record}"
        );
    }
}
