//! G1/G3: Newtonian fluid transport and actual caloric response reach the
//! native pressure/flow solve. Liquid coefficients here are manufactured
//! numerical fixtures, not a claim of measured material qualification.

use fs_duct::compose::{ChainResult, ComposeOptions};
use fs_duct::fluid::{self, FluidDuctError};
use fs_duct::{
    ConeSlices, Duct, DuctError, LossModel, RadiusProfile, RiccatiSteps, Segment, WallModel,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas_species::input::read_gas_record;
use fs_math::c64::C64;
use fs_phs::WallPin;
use std::sync::Arc;

fn options() -> ComposeOptions {
    ComposeOptions {
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}

fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}

fn liquid(viscosity: f64, conductivity: f64, cv: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, conductivity, viscosity).unwrap(),
        cv,
    )
    .unwrap()
}

fn near(actual: C64, expected: C64, scale: f64, tolerance: f64) {
    let error = (actual - expected).abs() / scale.max(expected.abs());
    assert!(error < tolerance, "{actual:?} != {expected:?}: {error}");
}

fn solve<M: ThermoviscousFluid + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    loss: LossModel,
    load: Option<C64>,
    wall: WallModel<'_>,
    slices: ConeSlices,
) -> ChainResult {
    fluid::input_impedance_load_with_loss(
        duct,
        material,
        omega,
        loss,
        load,
        wall,
        slices,
        options(),
        || false,
    )
    .unwrap()
}

#[test]
fn gas_mixtures_keep_original_bessel_and_wide_tube_propagation() {
    for text in [
        "mixture 293.15 100000 mole 2 N2 0.79 O2 0.21",
        "mixture 320 200000 mole 1 He 1",
        "mixture 300 100000 mole 1 CO2 1",
    ] {
        let mut words = text.split_ascii_whitespace();
        let record = read_gas_record(&mut words, None, None).unwrap();
        let gas = record.state();
        for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
            for omega in [100.0, 1000.0, 10000.0] {
                let original = fs_duct::segment_wave(gas, 0.02, omega, loss).unwrap();
                let generic = fluid::segment_wave_with_loss(gas, 0.02, omega, loss).unwrap();
                near(generic.wavenumber, original.wavenumber, 1e-20, 3e-13);
                near(
                    generic.characteristic_impedance,
                    original.characteristic_impedance,
                    1.0,
                    3e-13,
                );
                near(
                    generic.specific_impedance,
                    original.specific_impedance,
                    1.0,
                    3e-13,
                );
            }
            let duct = Duct {
                segments: vec![
                    Segment::Cylinder {
                        radius: 0.02,
                        length: 0.1,
                    },
                    Segment::Cone {
                        inlet_radius: 0.02,
                        outlet_radius: 0.03,
                        length: 0.2,
                    },
                    Segment::Profile {
                        profile: Arc::new(RadiusProfile::exponential(0.03, 1.0, 0.2).unwrap()),
                        steps: RiccatiSteps(24),
                    },
                ],
            };
            let wall = WallPin {
                surface_density: 0.2,
                stiffness_per_area: 4e7,
                resistance: 300.0,
            };
            for wall in [WallModel::Rigid, WallModel::Pin(&wall)] {
                let load = C64::new(2e5, -1e4);
                let original =
                    fs_duct::input_impedance_load_with_wall(&duct, gas, 700.0, loss, load, wall)
                        .unwrap();
                let generic = solve(&duct, gas, 700.0, loss, Some(load), wall, ConeSlices::Auto);
                near(generic.z_in, original.impedance, 1.0, 2e-11);
                near(
                    generic.u_out_over_u_in,
                    original.u_mouth_over_u_in,
                    1.0,
                    2e-11,
                );
                near(
                    generic.p_out_over_u_in,
                    original.p_mouth_over_u_in,
                    1.0,
                    2e-11,
                );
            }
        }
    }
}

#[test]
fn narrow_liquid_pipe_recovers_poiseuille_resistance_and_actual_isothermal_compliance() {
    let (radius, length, omega) = (1e-4_f64, 0.02, 0.01);
    let duct = cylinder(radius, length);
    let volume = core::f64::consts::PI * radius * radius * length;
    let mut closed_pressures = Vec::new();
    for cv in [4100.0, 4000.0] {
        let material = liquid(1e-3, 0.6, cv);
        let open = solve(
            &duct,
            &material,
            omega,
            LossModel::Bessel,
            Some(C64::ZERO),
            WallModel::Rigid,
            ConeSlices::Auto,
        );
        // Independent steady Hagen–Poiseuille limit at small Womersley number.
        let resistance = 8.0 * 1e-3 * length / (core::f64::consts::PI * radius.powi(4));
        assert!((open.z_in.re / resistance - 1.0).abs() < 2e-5);
        assert!(open.z_in.re > 0.0);
        // Isothermal compressibility is gamma/K_s. Keeping K_s and rho fixed
        // while changing only Cv must change the low-frequency closed response.
        let compliance = volume * (4200.0 / cv) / (1000.0 * 1500.0 * 1500.0);
        let closed = solve(
            &duct,
            &material,
            omega,
            LossModel::Bessel,
            None,
            WallModel::Rigid,
            ConeSlices::Auto,
        );
        assert!((closed.z_in.im * omega * compliance - 1.0).abs() < 2e-5);
        closed_pressures.push(closed.z_in.im);
        assert!(
            fluid::segment_wave_with_loss(&material, radius, omega, LossModel::WideTube).is_err()
        );
    }
    assert!((closed_pressures[1] / closed_pressures[0] - 4000.0 / 4100.0).abs() < 2e-5);
}

#[test]
fn viscosity_and_heat_diffusion_change_native_transmission_and_absorbed_power() {
    let (radius, length, omega) = (0.001, 3.0, 1e5);
    let duct = cylinder(radius, length);
    let mut transmitted = Vec::new();
    for viscosity in [1e-3, 4e-3] {
        let material = liquid(viscosity, 0.6, 4100.0);
        let wave =
            fluid::segment_wave_with_loss(&material, radius, omega, LossModel::Bessel).unwrap();
        let load = wave.characteristic_impedance;
        let result = solve(
            &duct,
            &material,
            omega,
            LossModel::Bessel,
            Some(load),
            WallModel::Rigid,
            ConeSlices::Auto,
        );
        let expected_flow = C64::new(
            (wave.wavenumber.re * length).cos(),
            (wave.wavenumber.re * length).sin(),
        )
        .scale((-wave.wavenumber.im * length).exp());
        near(result.z_in, load, 1.0, 3e-12);
        near(result.u_out_over_u_in, expected_flow, 1.0, 3e-12);
        let absorption = result.z_in.re - load.re * result.u_out_over_u_in.abs().powi(2);
        let expected_absorption = load.re * (1.0 - (-2.0 * wave.wavenumber.im * length).exp());
        assert!((absorption - expected_absorption).abs() < 4e-12 * load.re);
        assert!(absorption > 0.0);
        transmitted.push(result.u_out_over_u_in.abs());
    }
    assert!(transmitted[1] < 0.7 * transmitted[0]);

    // Independent Kirchhoff thin-boundary-layer attenuation, with the actual
    // liquid Cp/Cv and Pr. The full Bessel response tends to this first order.
    let material = liquid(1e-3, 0.6, 4100.0);
    let omega = 1e6;
    let nu: f64 = 1e-6;
    let pr: f64 = 1e-3 * 4200.0 / 0.6;
    let alpha =
        (nu * omega / 2.0).sqrt() / (radius * 1500.0) * (1.0 + (4200.0 / 4100.0 - 1.0) / pr.sqrt());
    let exact = fluid::segment_wave_with_loss(&material, radius, omega, LossModel::Bessel).unwrap();
    let wide =
        fluid::segment_wave_with_loss(&material, radius, omega, LossModel::WideTube).unwrap();
    assert!((exact.wavenumber.im / alpha - 1.0).abs() < 0.003);
    assert!((wide.wavenumber.im / alpha - 1.0).abs() < 3e-15);
}

#[test]
fn liquid_diffusion_uses_local_geometry_and_retains_wall_work() {
    let material = liquid(1e-3, 0.6, 4100.0);
    let (radius, length, omega) = (0.001, 0.3, 5000.0);
    let direct = cylinder(radius, length);
    let smooth = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(RadiusProfile::exponential(radius, 0.0, length).unwrap()),
            steps: RiccatiSteps(32),
        }],
    };
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: radius,
            outlet_radius: radius,
            length,
        }],
    };
    let load = C64::from_re(1000.0 * 1500.0 / (core::f64::consts::PI * radius * radius));
    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 2e12,
        resistance: 1e7,
    };
    for wall in [WallModel::Rigid, WallModel::Pin(&wall)] {
        let expected = solve(
            &direct,
            &material,
            omega,
            LossModel::Bessel,
            Some(load),
            wall,
            ConeSlices::Auto,
        );
        for geometry in [&smooth, &cone] {
            let result = solve(
                geometry,
                &material,
                omega,
                LossModel::Bessel,
                Some(load),
                wall,
                ConeSlices::Fixed(16),
            );
            near(result.z_in, expected.z_in, load.re, 2e-10);
            near(result.u_out_over_u_in, expected.u_out_over_u_in, 1.0, 2e-10);
        }
        assert!(expected.z_in.re - load.re * expected.u_out_over_u_in.abs().powi(2) > 0.0);
    }
    let tapered = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: radius,
            outlet_radius: 2.0 * radius,
            length,
        }],
    };
    let chain = fluid::duct_chain_with_loss(
        &tapered,
        &material,
        omega,
        LossModel::Bessel,
        WallModel::Rigid,
        ConeSlices::Fixed(16),
        options(),
        || false,
    )
    .unwrap();
    assert_eq!(chain.blocks.len(), 16);
    let first_midpoint = radius + radius / 32.0;
    let expected_shear = first_midpoint * (1000.0 * omega / 1e-3).sqrt();
    assert!((chain.min_shear_number / expected_shear - 1.0).abs() < 3e-15);
    let coarse = solve(
        &tapered,
        &material,
        omega,
        LossModel::Bessel,
        Some(load),
        WallModel::Rigid,
        ConeSlices::Fixed(8),
    );
    let medium = solve(
        &tapered,
        &material,
        omega,
        LossModel::Bessel,
        Some(load),
        WallModel::Rigid,
        ConeSlices::Fixed(16),
    );
    let fine = solve(
        &tapered,
        &material,
        omega,
        LossModel::Bessel,
        Some(load),
        WallModel::Rigid,
        ConeSlices::Fixed(64),
    );
    assert!((medium.z_in - fine.z_in).abs() < 0.4 * (coarse.z_in - fine.z_in).abs());
}

#[test]
fn wide_tube_thermal_admission_and_zero_thermal_coupling_are_physical() {
    // Low-Pr liquid: thin viscous layer does not imply a thin thermal layer.
    let conducting = liquid(1e-3, 42000.0, 4100.0);
    assert!(matches!(
        fluid::segment_wave_with_loss(&conducting, 0.001, 1e4, LossModel::WideTube),
        Err(FluidDuctError::Duct(DuctError::BadParameter { what })) if what.contains("thermal shear")
    ));
    let bessel = fluid::segment_wave_with_loss(&conducting, 0.001, 1e4, LossModel::Bessel).unwrap();
    assert!(bessel.wavenumber.im > 0.0);
    // A large caloric ratio also defeats a first-order thermal expansion even
    // when the thermal layer is thin. Full Bessel keeps that coupling intact.
    let large_caloric_ratio = liquid(1e-3, 0.6, 4.2);
    assert!(matches!(
        fluid::segment_wave_with_loss(&large_caloric_ratio, 0.001, 1e4, LossModel::WideTube),
        Err(FluidDuctError::Duct(DuctError::BadParameter { what })) if what.contains("coupled thermal")
    ));
    assert!(
        fluid::segment_wave_with_loss(&large_caloric_ratio, 0.001, 1e4, LossModel::Bessel).is_ok()
    );
    // Cp==Cv means no compressional thermal coupling. A thermal argument below
    // the Bessel kernel's own domain must not invalidate a response independent
    // of that unused diffusion channel.
    let no_thermal = liquid(1e-3, 4.2e20, 4200.0);
    let ordinary_conductivity = liquid(1e-3, 0.6, 4200.0);
    for loss in [LossModel::WideTube, LossModel::Bessel] {
        let a = fluid::segment_wave_with_loss(&no_thermal, 0.001, 1e4, loss).unwrap();
        let b = fluid::segment_wave_with_loss(&ordinary_conductivity, 0.001, 1e4, loss).unwrap();
        assert_eq!(a.wavenumber, b.wavenumber);
        assert_eq!(a.characteristic_impedance, b.characteristic_impedance);
    }
    let mut gas = fs_material::gas::GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        300.0,
        100_000.0,
    )
    .unwrap();
    gas.prandtl *= 1.01;
    assert!(matches!(
        fluid::segment_wave_with_loss(&gas, 0.001, 1e4, LossModel::Bessel),
        Err(FluidDuctError::Material(_))
    ));
}
