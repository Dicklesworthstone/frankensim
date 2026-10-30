use super::*;
use fs_duct::compose::{ComposeOptions, Composition, compose_chain};
use fs_duct::{Segment, SegmentWave};
use fs_feec::acoustics::helmholtz::Medium;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ThermalBacking, WallThermalTransport,
};
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::GasSpec;

const LENGTH: f64 = 0.06;
const AREA: f64 = 0.04 * 0.04;
const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
    minimum_curvature_radius_m: f64::INFINITY,
    minimum_half_gap_m: 0.02,
    maximum_layer_ratio: 0.05,
};

fn complete(fe: &HelmholtzFe<'_>) -> HelmholtzModes {
    let n = fe.ndof();
    HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes: (0..n)
            .map(|i| {
                let mut p = vec![0.0; n];
                p[i] = 1.0;
                p
            })
            .collect(),
        below_low: 0,
        below_high: n,
    }
}
fn end(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-14)
        })
        .collect()
}
fn medium(gas: &GasState) -> Medium {
    Medium {
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    }
}
fn analytic(gas: &GasState, length: f64, frequency: f64, reference: f64) -> TwoPort {
    TwoPort::from_wave(
        &Segment::Cylinder {
            radius: (AREA / core::f64::consts::PI).sqrt(),
            length,
        },
        &SegmentWave {
            wavenumber: C64::from_re(core::f64::consts::TAU * frequency / gas.sound_speed),
            characteristic_impedance: C64::from_re(reference),
            specific_impedance: C64::from_re(gas.density * gas.sound_speed),
            shear_number: f64::INFINITY,
        },
        reference,
    )
    .unwrap()
}
fn options() -> ComposeOptions {
    ComposeOptions {
        rule: Composition::Scattering,
        maximum_condition: 1e8,
        maximum_blocks: 2,
        preserve_well_conditioned_bits: false,
    }
}

#[test]
fn spatial_fe_section_composes_with_original_duct_including_complex_phase() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let reference = gas.density * gas.sound_speed / AREA;
    let tail_length = 0.17;
    let mut errors = Vec::new();
    for cells in [1, 2, 3] {
        let (mesh, points) = kuhn_cube(cells);
        let points: Vec<_> = points
            .iter()
            .map(|p| [p[0] * LENGTH, p[1] * 0.04, p[2] * 0.04])
            .collect();
        let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
        let modes = complete(&fe);
        let ports = [
            BoundaryMode::uniform(&end(&fe, 0.0)),
            BoundaryMode::uniform(&end(&fe, LENGTH)),
        ];
        let interior = interior_multiport(&fe, &modes, &ports, false, &gas, &|| false).unwrap();
        let mut worst = 0.0_f64;
        for frequency in [160.0, 480.0, 900.0] {
            let block = duct_two_port(&interior, frequency, reference, &|| false).unwrap();
            let plus = interior
                .scattering(frequency, &[reference; 2], &|| false)
                .unwrap();
            for (minus, original) in block.scattering.iter().zip(plus) {
                assert_eq!(minus.re.to_bits(), original.re.to_bits());
                assert_eq!(minus.im.to_bits(), (-original.im).to_bits());
            }
            // The exp(-i omega t) duct's forward delay is exp(+i*k*L).
            // This phase assertion fails if the FE matrix is copied unchanged.
            assert!(block.scattering[2].im > 0.0);
            let load = Some(C64::new(0.7 * reference, 0.3 * reference));
            let actual = compose_chain(
                &[block, analytic(&gas, tail_length, frequency, reference)],
                load,
                options(),
            )
            .unwrap();
            let expected = compose_chain(
                &[analytic(&gas, LENGTH + tail_length, frequency, reference)],
                load,
                options(),
            )
            .unwrap();
            for (a, b, scale) in [
                (actual.z_in, expected.z_in, reference),
                (actual.p_out_over_u_in, expected.p_out_over_u_in, reference),
                (actual.u_out_over_u_in, expected.u_out_over_u_in, 1.0),
            ] {
                worst = worst.max((a - b).abs() / scale);
            }
            for (a, b) in actual
                .scattering
                .scattering
                .iter()
                .zip(expected.scattering.scattering)
            {
                worst = worst.max((*a - b).abs());
            }
        }
        println!("FE to duct composition cells={cells} complex error={worst:e}");
        errors.push(worst);
    }
    assert!(errors.windows(2).all(|e| e[1] < e[0]), "{errors:?}");
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn source_gas_and_conductive_wall_reach_the_same_passive_duct_block() {
    let gases = [
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        GasState::try_new_moist_air(310.15, 101325.0, 0.8).unwrap(),
    ];
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * 0.04, p[2] * 0.04])
        .collect();
    let mut responses = Vec::new();
    for gas in gases {
        let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
        let modes = complete(&fe);
        let left = end(&fe, 0.0);
        let right = end(&fe, LENGTH);
        let ports = [BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)];
        let walls: Vec<_> = fe
            .boundary_faces()
            .into_iter()
            .filter(|f| !left.contains(f) && !right.contains(f))
            .collect();
        let build = || {
            thermoviscous_multiport(&fe, &modes, &ports, false, &walls, &gas, GEOMETRY, &|| {
                false
            })
            .unwrap()
        };
        let isothermal = build();
        let conductive = build()
            .with_thermal_wall(ThermalWall::Conductive(ConductiveThermalWall {
                transport: WallThermalTransport {
                    conductivity_w_m_k: 0.01,
                    volumetric_heat_capacity_j_m3_k: 1000.0,
                },
                depth: ThermalWallDepth::Slab {
                    thickness_m: 1e-4,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: 0.1,
            }))
            .unwrap();
        let reference = gas.density * gas.sound_speed / AREA;
        let f = 600.0;
        let a = duct_two_port(&isothermal, f, reference, &|| false).unwrap();
        let b = duct_two_port(&conductive, f, reference, &|| false).unwrap();
        assert_ne!(a.scattering, b.scattering);
        for block in [a, b] {
            for column in 0..2 {
                let returned = block.scattering[column].abs().powi(2)
                    + block.scattering[2 + column].abs().powi(2);
                assert!(returned < 1.0 && returned > 0.0);
            }
        }
        let response = conductive
            .prepare(f, &|| false)
            .unwrap()
            .response(
                &[C64::new(1e-6, 0.3e-6), C64::new(-0.4e-6, 0.2e-6)],
                &|| false,
            )
            .unwrap();
        let loss = conductive
            .boundary_layer_power(f, &response.amplitudes)
            .unwrap();
        let total = loss.viscous_w + loss.thermal_w;
        assert!(total > 0.0);
        assert!((response.input_power_w() - total).abs() < 5e-10 * total);
        let refs = [reference, 2.0 * reference];
        let generic = duct_scattering(&conductive, f, &refs, &|| false).unwrap();
        let original = conductive.scattering(f, &refs, &|| false).unwrap();
        assert_eq!(
            generic,
            original.into_iter().map(C64::conj).collect::<Vec<_>>()
        );
        responses.push(b.scattering);
    }
    assert_ne!(responses[0], responses[1]);
}

#[test]
fn source_mismatch_mode_truncation_and_cancel_refuse_before_export() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let (mesh, points) = kuhn_cube(1);
    let fe = HelmholtzFe::new(&mesh, &points, 1, medium(&gas)).unwrap();
    let modes = complete(&fe);
    let ports = [BoundaryMode::uniform(&end(&fe, 0.0))];
    let one = interior_multiport(&fe, &modes, &ports, false, &gas, &|| false).unwrap();
    assert!(duct_two_port(&one, 300.0, 100.0, &|| false).is_err());
    assert!(duct_scattering(&one, 300.0, &[100.0], &|| true).is_err());
    assert!(duct_scattering(&one, 300.0, &[0.0], &|| false).is_err());
    let foreign = HelmholtzFe::new(
        &mesh,
        &points,
        1,
        Medium {
            sound_speed_m_s: f64::from_bits(gas.sound_speed.to_bits() + 1),
            ..medium(&gas)
        },
    )
    .unwrap();
    assert!(interior_multiport(&foreign, &modes, &ports, false, &gas, &|| false).is_err());
    let corrupt = GasState {
        density: gas.density * 1.01,
        ..gas
    };
    assert!(interior_multiport(&fe, &modes, &ports, false, &corrupt, &|| false).is_err());
}
