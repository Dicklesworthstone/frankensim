//! G1/G3: real material/EOS substitutions reach the native pipe solve.
//! Analytic plane-wave and compliant-tube references are independent of the
//! wave-basis/scattering implementation. Supplied wall numbers are numerical
//! fixtures, not measured material data.

use fs_duct::compose::{self, ChainResult, ComposeOptions};
use fs_duct::fluid::{self, FluidDuctError};
use fs_duct::{
    ConeSlices, Duct, DuctError, LossModel, RadiusProfile, RiccatiSteps, Segment, WallModel,
};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::read_gas_record;
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};
use std::sync::Arc;

#[path = "fluid_acoustics/thermoviscous.rs"]
mod thermoviscous;

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

fn water(temperature: f64, pressure: f64) -> If97State {
    If97State::at_pressure_temperature(
        Pressure::new(pressure),
        Temperature::new(temperature),
        If97Region::Region1,
    )
    .unwrap()
}

fn near(a: C64, b: C64, scale: f64, tolerance: f64) {
    let error = (a - b).abs() / scale.max(b.abs());
    assert!(error < tolerance, "{a:?} != {b:?}: relative error {error}");
}

fn solve<M: FluidAcoustics + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    load: Option<C64>,
    wall: WallModel<'_>,
) -> ChainResult {
    fluid::input_impedance_load(
        duct,
        material,
        omega,
        load,
        wall,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap()
}

#[test]
fn gas_and_liquid_substitutions_drive_pipe_pressure_flow_and_travel_phase() {
    let gas = |text: &str| {
        let mut words = text.split_ascii_whitespace();
        let record = read_gas_record(&mut words, None, None).unwrap();
        assert!(words.next().is_none());
        *record.state()
    };
    let air = gas("mixture 293.15 100000 mole 2 N2 0.79 O2 0.21");
    let helium = gas("mixture 293.15 100000 mole 1 He 1");
    let cold = water(293.15, 100_000.0);
    let warm = water(333.15, 100_000.0);
    let compressed = water(293.15, 10_000_000.0);
    let materials: [&dyn FluidAcoustics; 5] = [&air, &helium, &cold, &warm, &compressed];
    let (radius, length) = (0.01, 0.7);
    let area = core::f64::consts::PI * radius * radius;
    let duct = cylinder(radius, length);
    let mut phases = Vec::new();
    for material in materials {
        let state = material.acoustic_state().unwrap();
        let zc = state.density_kg_m3() * state.sound_speed_m_s() / area;
        for omega in [71.0, 491.0, 1301.0] {
            let phase = omega * length / state.sound_speed_m_s();
            let matched = solve(
                &duct,
                material,
                omega,
                Some(C64::from_re(zc)),
                WallModel::Rigid,
            );
            let transmission = C64::new(phase.cos(), phase.sin());
            near(matched.z_in, C64::from_re(zc), zc, 2e-12);
            near(matched.u_out_over_u_in, transmission, 1.0, 2e-12);
            near(matched.p_out_over_u_in, transmission.scale(zc), zc, 2e-12);
            let outlet_power = zc * matched.u_out_over_u_in.abs().powi(2);
            assert!((matched.z_in.re - outlet_power).abs() < 2e-12 * zc);

            let closed = solve(&duct, material, omega, None, WallModel::Rigid);
            near(closed.z_in, C64::new(0.0, zc / phase.tan()), zc, 3e-12);
            near(closed.u_out_over_u_in, C64::ZERO, 1.0, 3e-12);
            near(
                closed.p_out_over_u_in,
                C64::new(0.0, zc / phase.sin()),
                zc,
                3e-12,
            );
            let open = solve(&duct, material, omega, Some(C64::ZERO), WallModel::Rigid);
            near(open.z_in, C64::new(0.0, -zc * phase.tan()), zc, 3e-12);
        }
        phases.push(491.0 * length / state.sound_speed_m_s());
    }
    // Same geometry and drive: composition, temperature and actual liquid
    // compression change propagation, without a material-name code path.
    assert!((phases[0] / phases[1] - 1.0).abs() > 0.5);
    assert!((phases[2] / phases[3] - 1.0).abs() > 0.01);
    assert!((phases[2] / phases[4] - 1.0).abs() > 0.001);
}

#[test]
fn gas_pressure_changes_impedance_without_fabricating_a_sound_speed_change() {
    let low = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let high = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200_000.0).unwrap();
    let duct = cylinder(0.012, 0.43);
    let a = solve(&duct, &low, 500.0, None, WallModel::Rigid);
    let b = solve(&duct, &high, 500.0, None, WallModel::Rigid);
    near(b.z_in, a.z_in.scale(2.0), a.z_in.abs(), 2e-12);
    for radius in [0.003, 0.05] {
        for omega in [10.0, 4000.0] {
            let original = fs_duct::segment_wave(&low, radius, omega, LossModel::Lossless).unwrap();
            let shared = fluid::segment_wave(&low, radius, omega).unwrap();
            for (a, b) in [
                (original.wavenumber, shared.wavenumber),
                (
                    original.characteristic_impedance,
                    shared.characteristic_impedance,
                ),
                (original.specific_impedance, shared.specific_impedance),
            ] {
                assert_eq!(
                    [a.re.to_bits(), a.im.to_bits()],
                    [b.re.to_bits(), b.im.to_bits()]
                );
            }
        }
    }
}

#[test]
fn liquid_bulk_compressibility_and_wall_compliance_share_one_wave_equation() {
    let water = water(300.0, 3_000_000.0);
    let state = water.acoustic_state().unwrap();
    let (radius, length, omega) = (0.01, 0.7, 500.0);
    let area = core::f64::consts::PI * radius * radius;
    let duct = cylinder(radius, length);
    let smooth = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(RadiusProfile::exponential(radius, 0.0, length).unwrap()),
            steps: RiccatiSteps(24),
        }],
    };
    let mut pressures = Vec::new();
    for stiffness in [2e11, 8e11] {
        let wall = WallPin {
            surface_density: 0.2,
            stiffness_per_area: stiffness,
            resistance: 0.0,
        };
        // Independent per-length inertia and compliance, including the wall's
        // own surface inertia; K_s is adiabatic, not the isothermal modulus.
        let inertance = state.density_kg_m3() / area;
        let compliance = area / state.adiabatic_bulk_modulus_pa()
            + core::f64::consts::TAU * radius / (stiffness - omega * omega * wall.surface_density);
        let k = omega * (inertance * compliance).sqrt();
        let zc = (inertance / compliance).sqrt();
        let expected = C64::new(0.0, zc / (k * length).tan());
        let response = solve(&duct, &water, omega, None, WallModel::Pin(&wall));
        near(response.z_in, expected, zc, 3e-12);
        near(
            solve(&smooth, &water, omega, None, WallModel::Pin(&wall)).z_in,
            expected,
            zc,
            3e-11,
        );
        pressures.push(response.z_in.im);
    }
    assert!((pressures[0] / pressures[1] - 1.0).abs() > 0.2);

    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 2e11,
        resistance: 1e8,
    };
    let load = C64::from_re(state.characteristic_impedance_pa_s_m() / area);
    let damped = solve(&duct, &water, omega, Some(load), WallModel::Pin(&wall));
    let absorbed = damped.z_in.re - load.re * damped.u_out_over_u_in.abs().powi(2);
    assert!(
        absorbed > 1e-5 * load.re,
        "wall dissipation did not reach the native port"
    );
}

#[test]
fn cylinder_cone_and_smooth_gas_chains_retain_existing_native_response() {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 100_000.0).unwrap();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.1,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.016,
                length: 0.2,
            },
            Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(0.016, 1.0, 0.2).unwrap()),
                steps: RiccatiSteps(24),
            },
        ],
    };
    let pin = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 4e7,
        resistance: 300.0,
    };
    for wall in [WallModel::Rigid, WallModel::Pin(&pin)] {
        for omega in [400.0, 1700.0] {
            let load = C64::new(1e6, -2e5);
            let original = fs_duct::input_impedance_load_with_wall(
                &duct,
                &air,
                omega,
                LossModel::Lossless,
                load,
                wall,
            )
            .unwrap();
            let shared = solve(&duct, &air, omega, Some(load), wall);
            near(shared.z_in, original.impedance, 1.0, 2e-12);
            near(
                shared.u_out_over_u_in,
                original.u_mouth_over_u_in,
                1.0,
                2e-12,
            );
            near(
                shared.p_out_over_u_in,
                original.p_mouth_over_u_in,
                1.0,
                2e-12,
            );
        }
    }
}

#[test]
fn source_failure_unsupported_geometry_work_caps_and_cancellation_are_explicit() {
    struct Unsupported;
    impl FluidAcoustics for Unsupported {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Err(MaterialError::State {
                what: "source phase no longer supported".into(),
            })
        }
    }
    let duct = cylinder(0.01, 0.7);
    assert_eq!(
        fluid::duct_chain(
            &duct,
            &Unsupported,
            500.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false
        )
        .unwrap_err(),
        FluidDuctError::Material(MaterialError::State {
            what: "source phase no longer supported".into()
        }),
    );
    assert!(matches!(
        fluid::duct_chain(
            &duct,
            &Unsupported,
            500.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        })),
    ));
    let mut invalid = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    invalid.sound_speed *= 1.01;
    assert!(matches!(
        fluid::segment_wave(&invalid, 0.01, 500.0),
        Err(FluidDuctError::Material(_))
    ));
    let water = water(300.0, 3_000_000.0);
    assert!(
        fluid::input_impedance_load(
            &duct,
            &water,
            500.0,
            Some(C64::from_re(-1.0)),
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
        .is_err()
    );
    let profile = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(RadiusProfile::exponential(0.01, 1.0, 0.7).unwrap()),
            steps: RiccatiSteps(24),
        }],
    };
    let bounded = ComposeOptions {
        maximum_blocks: 23,
        ..options()
    };
    assert!(
        fluid::duct_chain(
            &profile,
            &water,
            500.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            bounded,
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        fluid::duct_chain(
            &profile,
            &water,
            500.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || {
                polls += 1;
                polls == 5
            }
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. })),
    ));
    let holes = Duct {
        segments: vec![Segment::ToneHole {
            hole_radius: 0.003,
            chimney_height: 0.01,
            bore_radius: 0.01,
            state: fs_duct::HoleState::Open,
        }],
    };
    assert!(
        fluid::duct_chain(
            &holes,
            &water,
            500.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false
        )
        .is_err()
    );
    // A failed attempt did not alter any reusable source or geometry.
    let prepared = fluid::duct_chain(
        &duct,
        &water,
        500.0,
        WallModel::Rigid,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap();
    let reflection = compose::compose_chain_reflection(&prepared.blocks, None, options()).unwrap();
    assert!((reflection.reflection.abs() - 1.0).abs() < 2e-12);
}
