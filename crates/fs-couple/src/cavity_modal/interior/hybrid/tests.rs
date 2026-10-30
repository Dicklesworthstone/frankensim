use super::*;
use crate::cavity_modal::interior::thermoviscous::{ThermoviscousInterior, source_medium};
use fs_duct::compose::{ComposeOptions, Composition, compose_chain};
use fs_duct::{Segment, SegmentWave};
use fs_feec::acoustics::helmholtz::boundary_layer::{
    BoundaryLayerGeometry, BoundaryLayerTransport,
};
use fs_feec::acoustics::helmholtz::multiport::{BoundaryMode, InteriorMultiport};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, HelmholtzModes, Medium};
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits,
};
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxationState, RelaxingMedium, input::read_relaxation};
use fs_math::c64::C64;
use fs_rep_mesh::TetComplex;

const H_LENGTH: f64 = 0.04;
const L_LENGTH: f64 = 0.015;
const HEIGHT: f64 = 0.002;
const WIDTH: f64 = 0.004;
const BAND: [f64; 2] = [20.0, 2000.0];

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn positions(points: &[[f64; 3]], resolved: bool) -> Vec<[f64; 3]> {
    let (offset, length) = if resolved {
        (H_LENGTH, L_LENGTH)
    } else {
        (0.0, H_LENGTH)
    };
    points
        .iter()
        .map(|p| [offset + p[0] * length, p[1] * HEIGHT, p[2] * WIDTH])
        .collect()
}

fn medium(gas: &GasState) -> Medium {
    Medium {
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    }
}

fn ports(fe: &HelmholtzFe<'_>, start: f64, length: f64) -> Vec<BoundaryMode> {
    [start, start + length]
        .map(|x| {
            BoundaryMode::uniform(
                &fe.boundary_faces()
                    .into_iter()
                    .filter(|face| {
                        face.iter()
                            .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-14)
                    })
                    .collect::<Vec<_>>(),
            )
        })
        .to_vec()
}

fn pressure_model(
    mesh: &TetComplex,
    points: &[[f64; 3]],
    medium: Medium,
    extra_coordinate: bool,
    transport: Option<BoundaryLayerTransport>,
) -> InteriorMultiport {
    let fe = HelmholtzFe::new(mesh, points, 2, medium).unwrap();
    let n = fe.ndof();
    // The complete original pressure FE basis; no modal truncation or static
    // fit is used in this source/connection integration fixture.
    let basis = HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes: (0..n)
            .map(|i| {
                let mut v = vec![0.0; n];
                v[i] = 1.0;
                v
            })
            .collect(),
        below_low: 0,
        below_high: n,
    };
    let mut ports = ports(&fe, 0.0, H_LENGTH);
    if extra_coordinate {
        ports.push(ports[0].clone());
    }
    match transport {
        None => InteriorMultiport::new(&fe, &basis, &ports, false, &|| false),
        Some(transport) => {
            let walls = fe
                .boundary_faces()
                .into_iter()
                .filter(|face| {
                    !ports
                        .iter()
                        .any(|p| p.traces.iter().any(|(f, _)| f == face))
                })
                .collect::<Vec<_>>();
            InteriorMultiport::new_with_boundary_layers(
                &fe,
                &basis,
                &ports,
                false,
                &walls,
                transport,
                BoundaryLayerGeometry {
                    minimum_curvature_radius_m: f64::INFINITY,
                    minimum_half_gap_m: HEIGHT / 2.0,
                    maximum_layer_ratio: 0.05,
                },
                &|| false,
            )
        }
    }
    .unwrap()
}

fn resolved_model<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    gas: &GasState,
) -> ThermoviscousFe<'c> {
    let fe = HelmholtzFe::new(mesh, points, 2, medium(gas)).unwrap();
    let ports = ports(&fe, H_LENGTH, L_LENGTH);
    let walls = fe
        .boundary_faces()
        .into_iter()
        .filter(|face| {
            !ports
                .iter()
                .any(|p| p.traces.iter().any(|(f, _)| f == face))
        })
        .collect::<Vec<_>>();
    ThermoviscousFe::new(
        mesh,
        points,
        source_medium(gas, 0.0).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: walls,
            ports,
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
}

fn close(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 2e-9 * scale,
        "{a:?} != {b:?}, normalized error={:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn source_hybrid_fields_work_and_original_duct_composition_match_independent_closure() {
    let gas = gas();
    let (hm, cube) = kuhn_cube(2);
    let hp = positions(&cube, false);
    let (lm, cube) = kuhn_cube(2);
    let lp = positions(&cube, true);
    let h = pressure_model(&hm, &hp, medium(&gas), false, None);
    let source =
        ThermoviscousInterior::from_gas(resolved_model(&lm, &lp, &gas), gas, BAND).unwrap();
    let hybrid = HybridInterior::new(&h, &source, &[1], &[0], &|| false).unwrap();
    assert!(core::ptr::eq(hybrid.source(), &source));
    assert_eq!(hybrid.source().gas(), &gas);
    assert_eq!(hybrid.inner().external_port_count(), 2);
    let reference = gas.characteristic_impedance / (HEIGHT * WIDTH);
    let options = ComposeOptions {
        rule: Composition::Scattering,
        maximum_condition: 1e10,
        maximum_blocks: 2,
        preserve_well_conditioned_bits: false,
    };
    for f in [80.0, 300.0, 1100.0] {
        let prepared = hybrid
            .prepare(f, &[reference; 2], &[reference; 2], &|| false)
            .unwrap();
        let field = prepared
            .response(&[C64::new(1.0, 0.3), C64::new(-0.2, 0.1)], &|| false)
            .unwrap();
        assert!(!field.helmholtz.amplitudes.is_empty());
        assert!(!field.thermoviscous.velocity_m_s.is_empty());
        assert!(!field.thermoviscous.temperature_k.is_empty());
        assert!(field.power.thermoviscous.viscous_w > 0.0);
        assert!(field.power.thermoviscous.conductive_w > 0.0);
        assert!(field.power.external_input_w > 0.0);
        assert!(
            field.power.balance_residual_w().abs() < 2e-8 * field.power.external_input_w,
            "{:?}",
            field.power
        );
        let interface_scale =
            field.power.helmholtz_interface_w.abs() + field.power.thermoviscous_interface_w.abs();
        assert!(field.power.interface_residual_w().abs() < 2e-10 * interface_scale);
        let plus = prepared.scattering(&|| false).unwrap();
        let minus = hybrid
            .duct_scattering(f, &[reference; 2], &[reference; 2], &|| false)
            .unwrap();
        for (a, b) in minus.iter().zip(&plus) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), (-b.im).to_bits());
        }
        let block = hybrid.duct_two_port(f, reference, &|| false).unwrap();
        assert_eq!(block.scattering.as_slice(), minus.as_slice());
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
        // Independent closure in original H impedance / L admittance
        // coordinates. No hybrid scattering or transformer is used here.
        // Uniform traces give P_H1=P_L0 and Q_H1=-Q_L0. The attached duct
        // gives P_L1=-Z_tail*Q_L1, with inward flows on both FE owners.
        let z = h
            .impedance(f, &|| false)
            .unwrap()
            .iter()
            .map(|v| v.conj())
            .collect::<Vec<_>>();
        let y = source
            .admittance(f, &mut || false)
            .unwrap()
            .iter()
            .map(|v| v.conj())
            .collect::<Vec<_>>();
        let p1_over_p0 = -terminal.z_in * y[2] / (C64::ONE + terminal.z_in * y[3]);
        let y_loaded = y[0] + y[1] * p1_over_p0;
        let p0_over_input_flow = z[2] / (C64::ONE + z[3] * y_loaded);
        let expected_impedance = z[0] - z[1] * y_loaded * p0_over_input_flow;
        let tail_flow_over_input_flow = -(y[2] + y[3] * p1_over_p0) * p0_over_input_flow;
        close(joined.z_in, expected_impedance, reference);
        close(
            joined.u_out_over_u_in,
            terminal.u_out_over_u_in * tail_flow_over_input_flow,
            1.0,
        );
        close(
            joined.p_out_over_u_in,
            terminal.p_out_over_u_in * tail_flow_over_input_flow,
            reference,
        );
        assert!(joined.z_in.re > 0.0);
        eprintln!(
            "source hybrid f={f}Hz original H/L/duct impedance error={:e}",
            (joined.z_in - expected_impedance).abs() / reference
        );
    }
}

fn molecular(gas: GasState) -> RelaxationState {
    let data = "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

#[test]
fn source_hybrid_refuses_medium_molecular_band_reference_and_coordinate_mismatches() {
    let gas = gas();
    let (hm, cube) = kuhn_cube(1);
    let hp = positions(&cube, false);
    let (lm, cube) = kuhn_cube(1);
    let lp = positions(&cube, true);
    let h = pressure_model(&hm, &hp, medium(&gas), false, None);
    let source =
        ThermoviscousInterior::from_gas(resolved_model(&lm, &lp, &gas), gas, BAND).unwrap();
    let wrong_medium = Medium {
        density_kg_m3: gas.density * (1.0 + 1e-12),
        ..medium(&gas)
    };
    let wrong_h = pressure_model(&hm, &hp, wrong_medium, false, None);
    assert!(HybridInterior::new(&wrong_h, &source, &[1], &[0], &|| false).is_err());
    let transport = BoundaryLayerTransport {
        dynamic_viscosity_pa_s: gas.dynamic_viscosity,
        thermal_conductivity_w_m_k: gas.thermal_conductivity,
        cp_j_kg_k: gas.specific_heat_cp,
        heat_capacity_ratio: gas.gamma,
    };
    let same_wall = pressure_model(&hm, &hp, medium(&gas), false, Some(transport));
    assert!(HybridInterior::new(&same_wall, &source, &[1], &[0], &|| false).is_ok());
    let wrong_wall = pressure_model(
        &hm,
        &hp,
        medium(&gas),
        false,
        Some(BoundaryLayerTransport {
            thermal_conductivity_w_m_k: gas.thermal_conductivity * 1.01,
            ..transport
        }),
    );
    assert!(HybridInterior::new(&wrong_wall, &source, &[1], &[0], &|| false).is_err());
    let relaxing = ThermoviscousInterior::from_relaxation(
        resolved_model(&lm, &lp, &gas),
        molecular(gas),
        BAND,
    )
    .unwrap();
    // The current H owner has equilibrium real compression. Refuse a
    // molecular source rather than silently applying its calorics to L only.
    assert!(HybridInterior::new(&h, &relaxing, &[1], &[0], &|| false).is_err());
    assert!(HybridInterior::new(&h, &source, &[1], &[0], &|| true).is_err());
    let hybrid = HybridInterior::new(&h, &source, &[1], &[0], &|| false).unwrap();
    let reference = gas.characteristic_impedance / (HEIGHT * WIDTH);
    for f in [BAND[0] * 0.99, BAND[1] * 1.01, f64::NAN] {
        assert!(
            hybrid
                .prepare(f, &[reference; 2], &[reference; 2], &|| false)
                .is_err()
        );
    }
    assert!(
        hybrid
            .prepare(300.0, &[reference], &[reference; 2], &|| false)
            .is_err()
    );
    assert!(
        hybrid
            .prepare(300.0, &[reference; 2], &[0.0, reference], &|| false)
            .is_err()
    );
    assert!(
        hybrid
            .prepare(300.0, &[reference; 2], &[reference; 2], &|| true)
            .is_err()
    );
    assert!(
        hybrid
            .duct_scattering(300.0, &[reference; 2], &[reference; 2], &|| true)
            .is_err()
    );
    let h_three = pressure_model(&hm, &hp, medium(&gas), true, None);
    let three = HybridInterior::new(&h_three, &source, &[1], &[0], &|| false).unwrap();
    assert_eq!(three.inner().external_port_count(), 3);
    assert!(three.duct_two_port(300.0, reference, &|| false).is_err());
}
