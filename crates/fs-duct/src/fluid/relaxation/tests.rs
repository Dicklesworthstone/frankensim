use super::*;
use crate::compose::{self, ChainResult, Composition};
use crate::{DuctError, Segment};
use fs_material::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_math::c64::C64;
use fs_qty::Pressure;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn chart(pressure_pa: f64, lo: f64, hi: f64) -> RegionalWater {
    RegionalWater::new(
        Pressure::new(pressure_pa),
        RegionalWaterControls {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap()
}

fn options(rule: Composition) -> ComposeOptions {
    ComposeOptions {
        rule,
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}

fn cylinder(radius: f64, length: f64, pieces: usize) -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius,
                length: length / pieces as f64,
            };
            pieces
        ],
    }
}

fn solve(
    medium: &RelaxingWaterDuct<'_>,
    duct: &Duct,
    load: C64,
    options: ComposeOptions,
) -> ChainResult {
    let prepared = medium
        .duct_chain(duct, WallModel::Rigid, ConeSlices::Auto, options, || false)
        .unwrap();
    compose::compose_chain(&prepared.blocks, Some(load), options).unwrap()
}

fn near(actual: C64, expected: C64, scale: f64, tolerance: f64) {
    let error = (actual - expected).abs() / scale.max(expected.abs());
    assert!(
        error < tolerance,
        "{actual:?} != {expected:?}, error={error}"
    );
}

#[test]
fn g1_original_phase_density_and_relaxation_drive_matched_pipe_power_and_phase() {
    let (radius, length, tau) = (0.002, 0.15, 0.003);
    let area = core::f64::consts::PI * radius * radius;
    let duct = cylinder(radius, length, 1);
    for (pressure, lo, hi) in [(1e5, 300.0, 500.0), (20e6, 625.0, 645.0)] {
        let law = chart(pressure, lo, hi);
        for quality in [0.2, 0.6] {
            let thermal = law.state_at_quality(quality, || false).unwrap();
            let material = RelaxingRegionalWater::try_new(thermal, tau).unwrap();
            for omega in [0.01 / tau, 1.0 / tau, 100.0 / tau] {
                let medium = RelaxingWaterDuct::try_new(&material, omega).unwrap();
                assert!(core::ptr::eq(medium.material(), &material));
                assert_eq!(*medium.material().equilibrium().source(), thermal);
                assert_eq!(*medium.material().frozen().source(), thermal);
                assert_eq!(medium.angular_frequency_rad_s(), omega);
                assert_eq!(
                    medium.response().density_kg_m3().to_bits(),
                    thermal.bulk_density_kg_m3().to_bits()
                );
                let response = medium.response();
                let wave = medium.segment_wave(radius).unwrap();
                assert_eq!(wave.wavenumber, response.wavenumber_rad_m());
                assert!(wave.wavenumber.im > 0.0);
                let impedance = response.specific_impedance_pa_s_m().scale(1.0 / area);
                near(
                    wave.characteristic_impedance,
                    impedance,
                    impedance.abs(),
                    2e-14,
                );
                near(
                    (wave.specific_impedance * wave.wavenumber).scale(1.0 / omega),
                    C64::from_re(thermal.bulk_density_kg_m3()),
                    thermal.bulk_density_kg_m3(),
                    3e-14,
                );

                // Independent matched traveling wave, exp(-i omega t).
                let attenuation = (-wave.wavenumber.im * length).exp();
                let phase = wave.wavenumber.re * length;
                let transmission = C64::new(phase.cos(), phase.sin()).scale(attenuation);
                for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
                    let result = solve(&medium, &duct, impedance, options(rule));
                    near(result.z_in, impedance, impedance.abs(), 2e-11);
                    near(result.u_out_over_u_in, transmission, 1.0, 2e-11);
                    near(
                        result.p_out_over_u_in,
                        impedance * transmission,
                        impedance.abs(),
                        2e-11,
                    );
                    // Actual pressure/flow powers for unit inlet volume flow.
                    let inlet_power = 0.5 * result.z_in.re;
                    let outlet_power =
                        0.5 * (result.p_out_over_u_in * result.u_out_over_u_in.conj()).re;
                    let expected_loss = 0.5 * impedance.re * (1.0 - attenuation * attenuation);
                    assert!(inlet_power > 0.0 && outlet_power >= -1e-12 * inlet_power);
                    assert!(outlet_power <= inlet_power * (1.0 + 2e-11));
                    assert!(
                        ((inlet_power - outlet_power) - expected_loss).abs() < 3e-11 * inlet_power
                    );
                }
            }
        }
    }
}

#[test]
fn g1_native_chain_split_and_opaque_relaxing_line_preserve_the_physical_port() {
    let law = chart(1e5, 300.0, 500.0);
    let material =
        RelaxingRegionalWater::try_new(law.state_at_quality(0.3, || false).unwrap(), 0.003)
            .unwrap();
    let medium = RelaxingWaterDuct::try_new(&material, 700.0).unwrap();
    let radius = 0.004;
    let wave = medium.segment_wave(radius).unwrap();
    let load = wave.characteristic_impedance.scale(1.7);
    let opts = options(Composition::Scattering);
    let whole = solve(&medium, &cylinder(radius, 0.35, 1), load, opts);
    let split = solve(&medium, &cylinder(radius, 0.35, 32), load, opts);
    near(split.z_in, whole.z_in, load.abs(), 2e-11);
    near(split.u_out_over_u_in, whole.u_out_over_u_in, 1.0, 2e-11);
    near(
        split.p_out_over_u_in,
        whole.p_out_over_u_in,
        load.abs(),
        2e-11,
    );

    // 900 attenuation lengths would overflow a growing transfer basis.
    // The native scattering owner must retain the finite matched inlet.
    let opaque_length = 900.0 / wave.wavenumber.im;
    let opaque = solve(
        &medium,
        &cylinder(radius, opaque_length, 1),
        wave.characteristic_impedance,
        opts,
    );
    near(
        opaque.z_in,
        wave.characteristic_impedance,
        load.abs(),
        2e-11,
    );
    near(opaque.u_out_over_u_in, C64::ZERO, 1.0, 1e-300);
    near(opaque.p_out_over_u_in, C64::ZERO, load.abs(), 1e-300);
}

#[test]
fn g3_original_latent_heat_step_changes_density_and_pipe_response_at_fixed_temperature() {
    let law = chart(20e6, 625.0, 645.0);
    let initial = law.state_at_quality(0.2, || false).unwrap();
    let target = law.state_at_quality(0.6, || false).unwrap();
    let mass = 0.3;
    let heat = mass * (target.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
    let heated = law
        .advance_heat(initial, mass, heat, 2e-5, || false)
        .unwrap();
    assert_eq!(heated.state.temperature_k(), initial.temperature_k());
    assert!(heated.state.bulk_density_kg_m3() < initial.bulk_density_kg_m3());
    let before = RelaxingRegionalWater::try_new(initial, 0.002).unwrap();
    let after = RelaxingRegionalWater::try_new(heated.state, 0.002).unwrap();
    let a = RelaxingWaterDuct::try_new(&before, 500.0).unwrap();
    let b = RelaxingWaterDuct::try_new(&after, 500.0).unwrap();
    assert_eq!(*b.material().equilibrium().source(), heated.state);
    assert_eq!(
        b.response().density_kg_m3().to_bits(),
        heated.state.bulk_density_kg_m3().to_bits()
    );
    let wa = a.segment_wave(0.003).unwrap();
    let wb = b.segment_wave(0.003).unwrap();
    assert!((wb.characteristic_impedance / wa.characteristic_impedance - C64::ONE).abs() > 0.1);
    assert!((wb.wavenumber / wa.wavenumber - C64::ONE).abs() > 0.01);
    let duct = cylinder(0.003, 0.4, 3);
    let load = C64::from_re(wa.characteristic_impedance.re);
    let opts = options(Composition::Scattering);
    let pa = solve(&a, &duct, load, opts).z_in;
    let pb = solve(&b, &duct, load, opts).z_in;
    assert!((pb - pa).abs() > 0.01 * pa.abs());
}

#[test]
fn g0_frequency_geometry_work_bounds_and_cancellation_remain_native_refusals() {
    let law = chart(1e5, 300.0, 500.0);
    let material =
        RelaxingRegionalWater::try_new(law.state_at_quality(0.4, || false).unwrap(), 0.003)
            .unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            RelaxingWaterDuct::try_new(&material, omega),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
    }
    let medium = RelaxingWaterDuct::try_new(&material, 700.0).unwrap();
    assert!(matches!(
        medium.segment_wave(0.0),
        Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
    ));
    let opts = options(Composition::Scattering);
    assert!(matches!(
        medium.duct_chain(
            &Duct { segments: vec![] },
            WallModel::Rigid,
            ConeSlices::Auto,
            opts,
            || false
        ),
        Err(FluidDuctError::Duct(DuctError::EmptyDuct))
    ));
    let duct = cylinder(0.004, 0.3, 3);
    let bounded = ComposeOptions {
        maximum_blocks: 2,
        ..opts
    };
    assert!(
        medium
            .duct_chain(&duct, WallModel::Rigid, ConeSlices::Auto, bounded, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        medium.duct_chain(&duct, WallModel::Rigid, ConeSlices::Auto, opts, || {
            polls += 1;
            polls >= 3
        }),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let wave = medium.segment_wave(0.004).unwrap();
    let recovered = solve(&medium, &duct, wave.characteristic_impedance, opts);
    near(
        recovered.z_in,
        wave.characteristic_impedance,
        wave.characteristic_impedance.abs(),
        2e-11,
    );
}
