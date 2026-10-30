//! Regional sources must affect the original native pressure/flow operator.
//! Declared liquid tuples below are manufactured plumbing oracles, not measurements.

use fs_duct::compose::ComposeOptions;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgFluidRegion, DgMesh, DgSidePort, DgWallResolution, PortId,
};
use fs_duct::{ConeSlices, Duct, LossModel, Segment, WallModel};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{If97Region, If97State};
use fs_thermochem::water::transport::If97TransportState;

fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn liquid(rho: f64, speed: f64, viscosity: f64, cv: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(rho, speed).unwrap(),
        FluidTransportState::try_new(rho, 4200.0, 0.6, viscosity).unwrap(),
        cv,
    )
    .unwrap()
}
fn wall_resolution() -> DgWallResolution {
    DgWallResolution {
        poles_per_field: 16,
        quadrature_points: 3,
        band_rad_s: [50.0, 1000.0].map(|f| core::f64::consts::TAU * f),
        maximum_relative_error: 1e-4,
    }
}
fn water(temperature: f64) -> If97TransportState {
    let state = If97State::at_pressure_temperature(
        Pressure::new(3e6),
        Temperature::new(temperature),
        If97Region::Region1,
    )
    .unwrap();
    If97TransportState::from_state(If97FluidState::Gibbs(state)).unwrap()
}
fn lossless_input(load: C64, impedance: f64, phase: f64) -> C64 {
    let (sin, cos) = phase.sin_cos();
    (load.scale(cos) + C64::new(0.0, -impedance * sin))
        / (C64::from_re(cos) + load * C64::new(0.0, -sin / impedance))
}

#[test]
fn original_if97_regions_reflect_at_their_interface_and_match_the_outlet_material() {
    let (radius, length, split) = (0.008, 0.7, 0.3);
    let regions = [
        DgFluidRegion {
            end_m: split,
            material: water(300.0),
        },
        DgFluidRegion {
            end_m: length,
            material: water(450.0),
        },
    ];
    let acoustic: Vec<_> = regions
        .iter()
        .map(|region| region.material.acoustic_state().unwrap())
        .collect();
    let dg = DgDuct::from_material_regions(
        &cylinder(radius, length),
        &regions,
        DgMesh {
            elements: 6,
            order: 4,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Matched,
    )
    .unwrap();
    assert_eq!(dg.material_regions(), &regions);
    assert_eq!(dg.acoustic_regions(), acoustic.as_slice());
    assert!(dg.thermoviscous_regions().is_none());
    assert_eq!(
        dg.node_positions()
            .iter()
            .filter(|x| (**x - split).abs() < 1e-14)
            .count(),
        2
    );
    let area = core::f64::consts::PI * radius * radius;
    let z0 = acoustic[0].characteristic_impedance_pa_s_m() / area;
    let z1 = acoustic[1].characteristic_impedance_pa_s_m() / area;
    let omega = core::f64::consts::TAU * 600.0;
    let flow = 1e-9;
    let response = dg
        .harmonic_response(
            omega,
            &[(PortId(0), C64::from_re(flow))],
            20_000_000,
            1e-10,
            &mut || false,
        )
        .unwrap();
    let actual = response.ports[0].pressure_pa.scale(1.0 / flow);
    let phase = omega * split / acoustic[0].sound_speed_m_s();
    // The second region is terminated in its OWN matched impedance. Its
    // upstream interface reflection is R exp(2 i k_0 L_0), for exp(-i omega t).
    let expected = lossless_input(C64::from_re(z1), z0, phase);
    assert!((actual - expected).abs() / expected.abs() < 2e-5);
    let reflection = (actual - C64::from_re(z0)) / (actual + C64::from_re(z0));
    let coefficient = (z1 - z0) / (z1 + z0);
    let expected_reflection = C64::new(
        coefficient * (2.0 * phase).cos(),
        coefficient * (2.0 * phase).sin(),
    );
    assert!((reflection - expected_reflection).abs() < 2e-5);
    assert!(reflection.abs() > 1e-3);
    near(
        response.matched_absorption_w,
        response.input_power_w,
        response.input_power_w,
        2e-8,
    );
}

#[test]
fn interior_transport_changes_the_native_response_without_changing_the_ends() {
    let (radius, length) = (0.00025, 0.06);
    let baseline = liquid(1000.0, 1500.0, 0.001, 4200.0);
    let mut resistances = Vec::new();
    let mut identities = Vec::new();
    for (viscosity, cv, memory_dofs) in [
        (0.001, 4200.0, 144),
        (0.004, 4200.0, 144),
        (0.001, 4100.0, 192),
    ] {
        let middle = liquid(1000.0, 1500.0, viscosity, cv);
        let regions = [
            DgFluidRegion {
                end_m: 0.02,
                material: baseline,
            },
            DgFluidRegion {
                end_m: 0.04,
                material: middle,
            },
            DgFluidRegion {
                end_m: length,
                material: baseline,
            },
        ];
        let dg = DgDuct::from_material_regions(
            &cylinder(radius, length),
            &regions,
            DgMesh {
                elements: 3,
                order: 2,
            },
            DgBoundary::Port(PortId(0)),
            DgBoundary::IdealOpen,
        )
        .unwrap()
        .with_spatial_viscothermal(wall_resolution())
        .unwrap();
        assert_eq!(dg.material_regions(), &regions);
        assert_eq!(
            dg.thermoviscous_regions().unwrap(),
            &[baseline, middle, baseline]
        );
        // Every region has viscous histories. Only the final case has
        // Cp > Cv in its middle region, adding exactly its 48 thermal arms.
        assert_eq!(
            dg.spatial_viscothermal_admission().unwrap().memory_dofs,
            memory_dofs
        );
        assert_eq!(dg.total_dofs(), 2 * 3 * 3 + memory_dofs);
        let omega = core::f64::consts::TAU * 500.0;
        let flow = 1e-12;
        let response = dg
            .harmonic_response(
                omega,
                &[(PortId(0), C64::from_re(flow))],
                50_000_000,
                1e-10,
                &mut || false,
            )
            .unwrap();
        let actual = response.ports[0].pressure_pa.scale(1.0 / flow);
        // Independent frequency-domain exact Bessel chain, retaining each
        // original material. No DG diffusion history or masses enter it.
        let mut exact = C64::ZERO;
        for index in (0..regions.len()).rev() {
            let start = if index == 0 {
                0.0
            } else {
                regions[index - 1].end_m
            };
            exact = fs_duct::fluid::input_impedance_load_with_loss(
                &cylinder(radius, regions[index].end_m - start),
                &regions[index].material,
                omega,
                LossModel::Bessel,
                Some(exact),
                WallModel::Rigid,
                ConeSlices::Auto,
                ComposeOptions::default(),
                || false,
            )
            .unwrap()
            .z_in;
        }
        assert!(actual.re > 0.0);
        assert!((actual - exact).abs() / exact.abs() < 8e-4);
        near(
            response.wall_dissipation_w,
            response.input_power_w,
            response.input_power_w,
            2e-8,
        );
        resistances.push(actual.re);
        identities.push(dg.operator_identity());
    }
    assert!(resistances[1] > resistances[0]);
    assert_ne!(identities[0], identities[1]);
    assert_ne!(identities[0], identities[2]);
}

#[test]
fn regional_compressibility_reaches_the_original_compliant_wall_work_map() {
    let (radius, length, split, pressure) = (0.01, 1.0, 0.4, 100.0);
    let regions = [
        DgFluidRegion {
            end_m: split,
            material: liquid(1000.0, 1500.0, 0.001, 4200.0),
        },
        DgFluidRegion {
            end_m: length,
            material: liquid(700.0, 900.0, 0.002, 4200.0),
        },
    ];
    let pin = WallPin {
        surface_density: 20.0,
        stiffness_per_area: 2e12,
        resistance: 0.0,
    };
    let dg = DgDuct::from_material_regions(
        &cylinder(radius, length),
        &regions,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: PortId(7),
        downstream_element: 1,
    }])
    .unwrap()
    .with_compliant_wall(pin)
    .unwrap();
    assert_eq!(dg.continuous_operator().ports()[0].0, PortId(7));
    let m = dg.field_dofs();
    let x = dg.to_energy_coords(&vec![pressure; m], &vec![0.0; m]);
    let native = dg.continuous_operator();
    let n = native.dimension();
    let generator = native.generator_per_s();
    let derivative = |x: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| generator[i * n + j] * x[j]).sum())
            .collect()
    };
    let second = derivative(&derivative(&x));
    let (pressure_acceleration, _) = dg.from_energy_coords(&second);
    // Initially p is constant, q=0 and the wall is at rest. Radial mechanics
    // give w''=p/sigma; local continuity gives p''=-(2 K_s/r) p/sigma.
    // Inlet-only bulk would fail in every node of the second liquid region.
    for (region, values) in regions.iter().zip(pressure_acceleration.chunks(3)) {
        let expected = -2.0 * region.material.acoustic().adiabatic_bulk_modulus_pa() / radius
            * pressure
            / pin.surface_density;
        for value in values {
            near(*value, expected, expected.abs(), 2e-10);
        }
    }
    let mut state = dg.project_initial_state(|_| pressure, |_| 0.0);
    let area = core::f64::consts::PI * radius * radius;
    let expected_storage = 0.5
        * area
        * pressure
        * pressure
        * (split / regions[0].material.acoustic().adiabatic_bulk_modulus_pa()
            + (length - split) / regions[1].material.acoustic().adiabatic_bulk_modulus_pa());
    near(
        state.stored_energy_j(),
        expected_storage,
        expected_storage,
        2e-12,
    );
    let propagator = dg.propagator(2e-7, 1e-11).unwrap();
    for _ in 0..8 {
        let report = propagator.step(&mut state, &[]).unwrap();
        assert!(report.ledger.balance_checked);
        near(
            state.stored_energy_j(),
            expected_storage,
            expected_storage,
            2e-9,
        );
    }
    assert!(dg.compliant_wall_state(&state).unwrap().stored_energy_j > 0.0);
}

#[test]
fn uniform_fluid_in_a_taper_uses_the_same_spatial_gas_wall_kernel() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150_000.0).unwrap();
    let geometry = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.0005,
            outlet_radius: 0.0008,
            length: 0.03,
        }],
    };
    let mesh = DgMesh {
        elements: 2,
        order: 2,
    };
    let old = DgDuct::new(
        &geometry,
        &gas,
        mesh,
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_spatial_viscothermal(wall_resolution())
    .unwrap();
    let fluid = DgDuct::from_material(
        &geometry,
        &gas,
        mesh,
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_spatial_viscothermal(wall_resolution())
    .unwrap();
    assert_eq!(fluid.material(), &gas);
    assert_eq!(
        fluid.thermoviscous_material(),
        Some(gas.thermoviscous_state().unwrap())
    );
    assert_eq!(old.total_dofs(), fluid.total_dofs());
    let old_operator = old.continuous_operator();
    let fluid_operator = fluid.continuous_operator();
    let scale = old_operator
        .generator_per_s()
        .iter()
        .map(|x| x.abs())
        .fold(0.0, f64::max);
    for (old, fluid) in old_operator
        .generator_per_s()
        .iter()
        .zip(fluid_operator.generator_per_s())
    {
        near(*fluid, *old, scale, 4e-13);
    }
}

#[derive(Clone)]
struct SplitView {
    acoustic: FluidAcousticState,
    thermal: ThermoviscousFluidState,
}
impl FluidAcoustics for SplitView {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}
impl FluidTransport for SplitView {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.thermal.transport())
    }
}
impl ThermoviscousFluid for SplitView {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        Ok(self.thermal)
    }
}

#[test]
fn unresolved_regions_mesh_and_inconsistent_source_views_refuse() {
    let geometry = cylinder(0.001, 1.0);
    let source = liquid(1000.0, 1500.0, 0.001, 4100.0);
    let mesh = DgMesh {
        elements: 2,
        order: 2,
    };
    let make = |regions: &[DgFluidRegion<ThermoviscousFluidState>]| {
        DgDuct::from_material_regions(
            &geometry,
            regions,
            mesh,
            DgBoundary::Closed,
            DgBoundary::Closed,
        )
    };
    assert!(make(&[]).is_err());
    for ends in [[0.5, 0.5], [0.5, 0.9], [f64::NAN, 1.0], [0.0, 1.0]] {
        assert!(
            make(&ends.map(|end_m| DgFluidRegion {
                end_m,
                material: source
            }))
            .is_err()
        );
    }
    let regions = [
        DgFluidRegion {
            end_m: 0.25,
            material: source,
        },
        DgFluidRegion {
            end_m: 1.0,
            material: source,
        },
    ];
    let split_geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.001,
                length: 0.5,
            },
            Segment::Cylinder {
                radius: 0.002,
                length: 0.5,
            },
        ],
    };
    assert!(
        DgDuct::from_material_regions(
            &split_geometry,
            &regions,
            mesh,
            DgBoundary::Closed,
            DgBoundary::Closed,
        )
        .is_err()
    ); // Three combined geometry/material spans require three elements.
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150_000.0).unwrap();
    let mut stale = gas;
    stale.density *= 2.0;
    assert!(
        DgDuct::from_material_regions(
            &geometry,
            &[
                DgFluidRegion {
                    end_m: 0.5,
                    material: gas
                },
                DgFluidRegion {
                    end_m: 1.0,
                    material: stale
                },
            ],
            mesh,
            DgBoundary::Closed,
            DgBoundary::Closed,
        )
        .is_err()
    );
    let inconsistent = [DgFluidRegion {
        end_m: 1.0,
        material: SplitView {
            acoustic: FluidAcousticState::try_new(1000.0, 1200.0).unwrap(),
            thermal: source,
        },
    }];
    let bare = DgDuct::from_material_regions(
        &geometry,
        &inconsistent,
        mesh,
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    assert!(bare.with_spatial_viscothermal(wall_resolution()).is_err());
}
