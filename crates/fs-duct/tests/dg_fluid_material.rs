//! Original fluid evaluations reach the native time-domain pressure/flow owner.
//! Manufactured transport below tests constitutive plumbing, not measurements.

use fs_duct::compose::ComposeOptions;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgSidePort, PortId, PortInput};
use fs_duct::{ConeSlices, Duct, LossModel, RadiusProfile, RiccatiSteps, Segment, WallModel};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{If97Region, If97State};
use fs_thermochem::water::transport::If97TransportState;
use std::sync::Arc;

fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn mesh() -> DgMesh {
    DgMesh {
        elements: 2,
        order: 3,
    }
}
fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn liquid(mu: f64, cv: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap(),
        cv,
    )
    .unwrap()
}

#[test]
fn original_if97_pressure_state_reaches_native_modes_and_characteristic_work() {
    let length = 1.0;
    let radius = 0.00025;
    let geometry = cylinder(radius, length);
    let mut impedances = Vec::new();
    let mut lossy_impedances = Vec::new();
    for pressure in [3e6, 40e6] {
        let water = If97State::at_pressure_temperature(
            Pressure::new(pressure),
            Temperature::new(300.0),
            If97Region::Region1,
        )
        .unwrap();
        let water = If97TransportState::from_state(If97FluidState::Gibbs(water)).unwrap();
        let acoustic = water.acoustic_state().unwrap();
        let dg = DgDuct::from_material(
            &geometry,
            &water,
            DgMesh {
                elements: 4,
                order: 4,
            },
            DgBoundary::Closed,
            DgBoundary::Closed,
        )
        .unwrap();
        assert_eq!(dg.material(), &water);
        assert_eq!(dg.acoustic_material(), acoustic);
        let amplitude = 500.0;
        let mut state = dg.project_initial_state(
            |x| amplitude * (core::f64::consts::PI * x / length).cos(),
            |_| 0.0,
        );
        let initial = state.stored_energy_j();
        let expected_storage =
            core::f64::consts::PI * radius * radius * length * amplitude * amplitude
                / (4.0 * acoustic.adiabatic_bulk_modulus_pa());
        near(initial, expected_storage, expected_storage, 2e-5);
        let dt = length / (8.0 * acoustic.sound_speed_m_s());
        let propagator = dg.propagator(dt, 1e-10).unwrap();
        for step in 1..=8 {
            let report = propagator.step(&mut state, &[]).unwrap();
            assert!(report.ledger.balance_checked);
            let factor = (core::f64::consts::PI * step as f64 / 8.0).cos();
            let error = dg
                .l2_pressure_error(&state, |x| {
                    amplitude * factor * (core::f64::consts::PI * x / length).cos()
                })
                .unwrap();
            assert!(
                error / amplitude < 1e-4,
                "native standing-wave error {error}"
            );
            near(state.stored_energy_j(), initial, initial, 3e-10);
        }
        let open = DgDuct::from_material(
            &geometry,
            &water,
            DgMesh {
                elements: 4,
                order: 4,
            },
            DgBoundary::Port(PortId(0)),
            DgBoundary::IdealOpen,
        )
        .unwrap();
        let omega = 400.0;
        let flow = 1e-9;
        let response = open
            .harmonic_response(
                omega,
                &[(PortId(0), C64::from_re(flow))],
                20_000_000,
                1e-10,
                &mut || false,
            )
            .unwrap();
        let actual = response.ports[0].pressure_pa.scale(1.0 / flow);
        let impedance =
            acoustic.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * radius * radius);
        let expected = C64::new(
            0.0,
            -impedance * (omega * length / acoustic.sound_speed_m_s()).tan(),
        );
        assert!((actual - expected).abs() / expected.abs() < 1e-6);
        impedances.push(actual.im);

        // Pressure-dependent IAPWS viscosity/conductivity and the SAME IF97
        // Cp/Cv now enter the original DG diffusion memories, not just rho/c.
        let lossy = DgDuct::from_material(
            &geometry,
            &water,
            DgMesh {
                elements: 2,
                order: 2,
            },
            DgBoundary::Port(PortId(0)),
            DgBoundary::IdealOpen,
        )
        .unwrap()
        .with_viscothermal(16, [20.0, 200.0], 1e-4)
        .unwrap();
        let tuple = water.thermoviscous_state().unwrap();
        assert_eq!(lossy.material().thermoviscous_state().unwrap(), tuple);
        near(
            lossy.viscothermal_model().unwrap().series_drag_per_s(),
            8.0 * water.dynamic_viscosity().value() / (acoustic.density_kg_m3() * radius * radius),
            lossy.viscothermal_model().unwrap().series_drag_per_s(),
            3e-15,
        );
        let response = lossy
            .harmonic_response(
                omega,
                &[(PortId(0), C64::from_re(flow))],
                50_000_000,
                1e-10,
                &mut || false,
            )
            .unwrap();
        let actual = response.ports[0].pressure_pa.scale(1.0 / flow);
        let exact = fs_duct::fluid::input_impedance_load_with_loss(
            &geometry,
            &water,
            omega,
            LossModel::Bessel,
            Some(C64::ZERO),
            WallModel::Rigid,
            ConeSlices::Auto,
            ComposeOptions::default(),
            || false,
        )
        .unwrap()
        .z_in;
        assert!((actual - exact).abs() / exact.abs() < 3e-4);
        assert!(actual.re > 0.0);
        near(
            response.wall_dissipation_w,
            response.input_power_w,
            response.input_power_w,
            2e-8,
        );
        lossy_impedances.push(actual);
    }
    assert_ne!(impedances[0], impedances[1]);
    assert_ne!(lossy_impedances[0], lossy_impedances[1]);
}

#[test]
fn legacy_gas_and_generic_owner_share_geometry_ports_and_exact_time_propagation() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150_000.0).unwrap();
    for segment in [
        Segment::Cylinder {
            radius: 0.01,
            length: 0.3,
        },
        Segment::Cone {
            inlet_radius: 0.01,
            outlet_radius: 0.02,
            length: 0.3,
        },
        Segment::Profile {
            profile: Arc::new(RadiusProfile::exponential(0.01, 1.0, 0.3).unwrap()),
            steps: RiccatiSteps(16),
        },
    ] {
        let geometry = Duct {
            segments: vec![segment],
        };
        for outlet in [DgBoundary::Closed, DgBoundary::Matched] {
            let sides = [DgSidePort {
                id: PortId(1),
                downstream_element: 1,
            }];
            let old = DgDuct::new(&geometry, &gas, mesh(), DgBoundary::Port(PortId(0)), outlet)
                .unwrap()
                .with_side_ports(&sides)
                .unwrap();
            let fluid =
                DgDuct::from_material(&geometry, &gas, mesh(), DgBoundary::Port(PortId(0)), outlet)
                    .unwrap()
                    .with_side_ports(&sides)
                    .unwrap();
            assert_eq!(old.operator_identity(), fluid.operator_identity());
            let a = old.propagator(1e-5, 1e-11).unwrap();
            let b = fluid.propagator(1e-5, 1e-11).unwrap();
            let mut sa = old.zero_state();
            let mut sb = fluid.zero_state();
            for step in 0..6 {
                let input = [
                    PortInput {
                        id: PortId(0),
                        effort: 1e-6,
                    },
                    PortInput {
                        id: PortId(1),
                        effort: if step < 3 { -2e-7 } else { 0.0 },
                    },
                ];
                assert_eq!(
                    a.step(&mut sa, &input).unwrap(),
                    b.step(&mut sb, &input).unwrap()
                );
            }
            assert_eq!(sa, sb);
        }
    }
}

#[test]
fn liquid_viscosity_and_caloric_response_drive_native_diffusion_and_pressure() {
    let geometry = cylinder(0.00025, 0.05);
    let mut resistances = Vec::new();
    for (mu, cv) in [(0.001, 4100.0), (0.002, 4100.0), (0.001, 4200.0)] {
        let material = liquid(mu, cv);
        let bare = DgDuct::from_material(
            &geometry,
            &material,
            DgMesh {
                elements: 2,
                order: 2,
            },
            DgBoundary::Port(PortId(0)),
            DgBoundary::IdealOpen,
        )
        .unwrap();
        let dg = bare
            .clone()
            .with_viscothermal(16, [50.0, 2000.0], 1e-4)
            .unwrap();
        assert_eq!(dg.material().thermoviscous_state().unwrap(), material);
        assert!(dg.has_viscothermal_memory());
        assert_eq!(
            dg.viscothermal_model().unwrap().shunt_terms().is_empty(),
            cv == 4200.0
        );
        for frequency in [100.0, 500.0, 1000.0] {
            let omega = core::f64::consts::TAU * frequency;
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
            let exact = fs_duct::fluid::input_impedance_load_with_loss(
                &geometry,
                &material,
                omega,
                LossModel::Bessel,
                Some(C64::ZERO),
                WallModel::Rigid,
                ConeSlices::Auto,
                ComposeOptions::default(),
                || false,
            )
            .unwrap()
            .z_in;
            assert!(actual.re > 0.0);
            assert!((actual - exact).abs() / exact.abs() < 3e-4);
            near(
                response.wall_dissipation_w,
                response.input_power_w,
                response.input_power_w,
                2e-8,
            );
            if frequency == 500.0 {
                resistances.push(actual.re);
            }
        }
        let propagator = dg.propagator(1e-5, 1e-10).unwrap();
        let mut state =
            dg.project_initial_state(|x| 10.0 * (core::f64::consts::PI * x / 0.05).sin(), |_| 0.0);
        let initial = state.stored_energy_j();
        let mut previous = initial;
        for _ in 0..16 {
            let result = propagator.step(&mut state, &[]).unwrap();
            assert!(!result.ledger.balance_checked);
            assert!(state.stored_energy_j() <= previous + 2e-10 * initial);
            assert!(dg.wall_dissipated_power_w(&state).unwrap() >= 0.0);
            previous = state.stored_energy_j();
        }
        assert!(previous < initial);
        let mut foreign = bare.zero_state();
        let retained = foreign.clone();
        assert!(propagator.step(&mut foreign, &[]).is_err());
        assert_eq!(foreign, retained);
    }
    assert!(resistances[1] > resistances[0]);
}

#[test]
fn liquid_pressure_and_compliant_wall_share_the_native_physical_energy() {
    let (radius, length, p0) = (0.01, 0.3, 100.0);
    let material = liquid(0.001, 4100.0);
    let pin = WallPin {
        surface_density: 20.0,
        stiffness_per_area: 2e12,
        resistance: 2e6,
    };
    let dg = DgDuct::from_material(
        &cylinder(radius, length),
        &material,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_compliant_wall(pin)
    .unwrap();
    let bulk = material.acoustic().adiabatic_bulk_modulus_pa();
    let beta = 2.0 * bulk / radius;
    let stiffness = pin.stiffness_per_area + beta;
    let damping = pin.resistance / (2.0 * pin.surface_density);
    let omega = (stiffness / pin.surface_density - damping * damping).sqrt();
    let dt = 1e-7;
    let propagator = dg.propagator(dt, 1e-11).unwrap();
    let mut state = dg.project_initial_state(|_| p0, |_| 0.0);
    let initial = state.stored_energy_j();
    for step in 1..=24 {
        propagator.step(&mut state, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        let t = step as f64 * dt;
        let (sin, cos) = (omega * t).sin_cos();
        let w = p0 / stiffness * (1.0 - (-damping * t).exp() * (cos + damping / omega * sin));
        let v = p0 / (pin.surface_density * omega) * (-damping * t).exp() * sin;
        let wall = dg.compliant_wall_state(&state).unwrap();
        for ((pressure, displacement), velocity) in state
            .pressure()
            .iter()
            .zip(&wall.displacement_m)
            .zip(&wall.velocity_m_s)
        {
            near(*pressure, p0 - beta * w, p0, 3e-9);
            near(*displacement, w, p0 / stiffness, 3e-9);
            near(*velocity, v, p0 / (pin.surface_density * omega), 3e-9);
        }
        let area = core::f64::consts::TAU * radius * length;
        near(
            wall.stored_energy_j,
            0.5 * area * (pin.surface_density * v * v + pin.stiffness_per_area * w * w),
            initial,
            3e-9,
        );
        near(
            dg.wall_dissipated_power_w(&state).unwrap(),
            area * pin.resistance * v * v,
            initial / dt,
            3e-9,
        );
    }
}
