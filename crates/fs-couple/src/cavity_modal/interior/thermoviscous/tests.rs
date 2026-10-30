use super::*;
use fs_duct::compose::{ComposeOptions, Composition, compose_chain};
use fs_duct::{Segment, SegmentWave};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::acoustics::thermoviscous::{ThermoviscousBoundary, ThermoviscousLimits};
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::GasSpec;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_rep_mesh::TetComplex;

mod fluid;

const LENGTH: f64 = 0.015;
const HEIGHT: f64 = 0.002;
const WIDTH: f64 = 0.004;
const BAND: [f64; 2] = [20.0, 2000.0];

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn source(gas: GasState, relaxing: bool, rate_scale: f64) -> RelaxationState {
    // Synthetic reservoir isolates the caloric equations; this is not a
    // measured gas absorption law or an accuracy claim about terrestrial air.
    let data = if relaxing {
        format!(
            "declared 1 synthetic 1 vibration 1 900 1 power {} 293.15 0 200 500",
            101325.0 * 0.001 * rate_scale
        )
    } else {
        "declared 1 synthetic 0".to_owned()
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

fn points(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * HEIGHT, p[2] * WIDTH])
        .collect()
}

fn build<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    gas: &GasState,
    isothermal: bool,
) -> ThermoviscousFe<'c> {
    let geometry = HelmholtzFe::new(
        mesh,
        points,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let end = |x: f64| {
        faces
            .iter()
            .copied()
            .filter(|face| {
                face.iter()
                    .all(|&v| (points[v as usize][0] - x).abs() < 1e-14)
            })
            .collect::<Vec<_>>()
    };
    let left = end(0.0);
    let right = end(LENGTH);
    let walls: Vec<_> = faces
        .into_iter()
        .filter(|f| !left.contains(f) && !right.contains(f))
        .collect();
    ThermoviscousFe::new(
        mesh,
        points,
        source_medium(gas, 0.0).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: if isothermal { walls } else { Vec::new() },
            ports: vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
}

fn close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:?} != {b:?}, scale={scale:e}"
    );
}

#[test]
fn original_gas_and_molecular_reservoirs_reach_resolved_work_and_fields() {
    let (mesh, xyz) = kuhn_cube(2);
    let xyz = points(&xyz);
    let gas = gas();
    let f = 300.0;
    let drive = [C64::new(1.0, 0.2), C64::new(0.8, -0.1)];
    for isothermal in [false, true] {
        let real = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, isothermal), gas, BAND)
            .unwrap();
        let empty = ThermoviscousInterior::from_relaxation(
            build(&mesh, &xyz, &gas, isothermal),
            source(gas, false, 1.0),
            BAND,
        )
        .unwrap();
        let molecular = ThermoviscousInterior::from_relaxation(
            build(&mesh, &xyz, &gas, isothermal),
            source(gas, true, 1.0),
            BAND,
        )
        .unwrap();
        let changed = ThermoviscousInterior::from_relaxation(
            build(&mesh, &xyz, &gas, isothermal),
            source(gas, true, 0.1),
            BAND,
        )
        .unwrap();
        let base = real.response(f, &drive, &mut || false).unwrap();
        let parity = empty.response(f, &drive, &mut || false).unwrap();
        for (a, b) in base.inward_flow_m3_s.iter().zip(parity.inward_flow_m3_s) {
            close(*a, b, a.abs().max(1e-30), 2e-10);
        }
        let result = molecular.response(f, &drive, &mut || false).unwrap();
        let rate = changed.response(f, &drive, &mut || false).unwrap();
        assert!(result.power.viscous_w > 0.0);
        assert!(result.power.conductive_w > 0.0);
        assert!(result.power.caloric_w > 0.0);
        assert!(result.input_power_w() > 0.0);
        assert!(result.power.balance_residual_w().abs() < 2e-9 * result.input_power_w());
        let effect = result
            .inward_flow_m3_s
            .iter()
            .zip(&rate.inward_flow_m3_s)
            .map(|(a, b)| (*a - *b).abs())
            .fold(0.0_f64, f64::max)
            / result
                .inward_flow_m3_s
                .iter()
                .map(|z| z.abs())
                .fold(0.0_f64, f64::max);
        println!(
            "resolved molecular isothermal={isothermal} rate_effect={effect:e} work={:?}",
            result.power
        );
        assert!(effect > 1e-5);
        assert_eq!(molecular.gas(), &gas);
        assert_eq!(molecular.relaxation().unwrap().gas(), &gas);
        assert_eq!(molecular.band_hz(), BAND);
    }
}

#[test]
fn resolved_molecular_feature_closes_with_original_duct_and_complex_load() {
    let gas = gas();
    let (mesh, xyz) = kuhn_cube(2);
    let xyz = points(&xyz);
    let model = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let reference = gas.characteristic_impedance / (HEIGHT * WIDTH);
    let options = ComposeOptions {
        rule: Composition::Scattering,
        maximum_condition: 1e10,
        maximum_blocks: 2,
        preserve_well_conditioned_bits: false,
    };
    for f in [80.0, 300.0, 1100.0] {
        let block = model.duct_two_port(f, reference, &mut || false).unwrap();
        let cp = model
            .relaxation()
            .unwrap()
            .caloric_response(core::f64::consts::TAU * f)
            .unwrap()
            .cp_j_kg_k
            .conj();
        let plus = model
            .inner()
            .scattering_with_heat_capacity(f, cp, &[reference; 2], &mut || false)
            .unwrap();
        for (a, b) in block.scattering.iter().zip(plus) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), (-b.im).to_bits());
        }
        let tail = TwoPort::from_wave(
            &Segment::Cylinder {
                radius: (HEIGHT * WIDTH / core::f64::consts::PI).sqrt(),
                length: 0.11,
            },
            &SegmentWave {
                wavenumber: C64::from_re(core::f64::consts::TAU * f / gas.sound_speed),
                characteristic_impedance: C64::from_re(reference),
                specific_impedance: C64::from_re(gas.characteristic_impedance),
                shear_number: f64::INFINITY,
            },
            reference,
        )
        .unwrap();
        let load = Some(C64::new(0.8 * reference, 0.3 * reference));
        let terminal = compose_chain(&[tail], load, options).unwrap();
        let joined = compose_chain(&[block, tail], load, options).unwrap();
        // Independently close the traction/flow equations using the original
        // FE admittance, rather than converting its scattering block again.
        let y: Vec<_> = model
            .admittance(f, &mut || false)
            .unwrap()
            .iter()
            .map(|z| z.conj())
            .collect();
        let p2 = -terminal.z_in * y[2] / (C64::ONE + terminal.z_in * y[3]);
        let q1 = y[0] + y[1] * p2;
        let tail_flow = -(y[2] + y[3] * p2);
        close(joined.z_in, q1.recip(), reference, 1e-9);
        close(
            joined.u_out_over_u_in,
            terminal.u_out_over_u_in * tail_flow / q1,
            1.0,
            1e-9,
        );
        close(
            joined.p_out_over_u_in,
            terminal.p_out_over_u_in * tail_flow / q1,
            reference,
            1e-9,
        );
        assert!(joined.z_in.re > 0.0);
    }
}

#[test]
fn source_gas_mismatch_band_and_cancellation_refuse_before_results() {
    let gas = gas();
    let (mesh, xyz) = kuhn_cube(1);
    let xyz = points(&xyz);
    let humid = GasState::try_new_moist_air(305.0, 101325.0, 0.8).unwrap();
    assert!(ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), humid, BAND).is_err());
    let corrupt = GasState {
        specific_heat_cp: gas.specific_heat_cp * 1.01,
        ..gas
    };
    assert!(source_medium(&corrupt, 0.0).is_err());
    assert!(source_medium(&gas, -1.0).is_err());
    assert!(
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, [0.0, 20.0]).is_err()
    );
    let model = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    assert!(model.admittance(BAND[0] * 0.99, &mut || false).is_err());
    assert!(
        model
            .duct_two_port(BAND[1] * 1.01, 1e7, &mut || false)
            .is_err()
    );
    assert!(model.response(100.0, &[C64::ONE; 2], &mut || true).is_err());
    assert!(model.duct_two_port(100.0, 0.0, &mut || false).is_err());
}
