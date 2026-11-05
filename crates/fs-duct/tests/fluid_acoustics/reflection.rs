//! Finite material reflection at pressure antinodes and physical interfaces.
use super::*;
use fs_duct::compose::Composition;
use fs_duct::fluid::regional;
use fs_duct::network::NetworkLoad;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

fn phase(theta: f64) -> C64 {
    C64::new(theta.cos(), theta.sin())
}

#[test]
fn gas_and_liquid_reflections_remain_finite_at_impedance_poles() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let liquid = water(300.0, 2_000_000.0);
    let geometry = cylinder(0.012, 0.4);
    for material in [&gas as &dyn FluidAcoustics, &liquid] {
        let c = material.acoustic_state().unwrap().sound_speed_m_s();
        for (load, theta, sign) in [
            (None, core::f64::consts::PI, 1.0),
            (Some(C64::ZERO), core::f64::consts::FRAC_PI_2, -1.0),
        ] {
            let pole = theta * c / 0.4;
            for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
                let opts = ComposeOptions { rule, ..options() };
                assert!(matches!(
                    fluid::input_impedance_load(
                        &geometry,
                        material,
                        pole,
                        load,
                        WallModel::Rigid,
                        ConeSlices::Auto,
                        opts,
                        || false,
                    ),
                    Err(FluidDuctError::Duct(DuctError::Singular))
                ));
                for offset in [0.0, -1e-5, 1e-5] {
                    let actual = fluid::input_reflection_load(
                        &geometry,
                        material,
                        pole * (1.0 + offset),
                        load,
                        WallModel::Rigid,
                        ConeSlices::Auto,
                        opts,
                        || false,
                    )
                    .unwrap();
                    // Two propagation traversals with exp(-i omega t).
                    near(
                        actual.reflection,
                        phase(2.0 * theta * (1.0 + offset)).scale(sign),
                        1.0,
                        2e-12,
                    );
                    assert!((actual.reflection.abs() - 1.0).abs() < 2e-12);
                    if offset == 0.0 {
                        near(
                            actual
                                .referenced_to(3.7 * actual.reference_impedance)
                                .unwrap(),
                            C64::ONE,
                            1.0,
                            2e-12,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn regional_reflection_retains_interface_contrast_delay_and_power() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 2_000_000.0).unwrap();
    let liquid = water(300.0, 2_000_000.0);
    let radius = 0.012;
    let omega = 900.0;
    let geometry = Duct {
        segments: [0.2, 0.37]
            .map(|length| Segment::Cylinder { radius, length })
            .to_vec(),
    };
    for materials in [
        [&gas as &dyn FluidAcoustics, &liquid],
        [&liquid as &dyn FluidAcoustics, &gas],
    ] {
        let states = materials.map(|m| m.acoustic_state().unwrap());
        let area = core::f64::consts::PI * radius * radius;
        let [z1, z2] = states.map(|s| s.characteristic_impedance_pa_s_m() / area);
        let interface = (z2 - z1) / (z2 + z1);
        let expected = phase(2.0 * omega * 0.2 / states[0].sound_speed_m_s()).scale(interface);
        for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
            let actual = regional::input_reflection_load(
                &geometry,
                &materials,
                omega,
                Some(C64::from_re(z2)),
                WallModel::Rigid,
                ConeSlices::Auto,
                ComposeOptions { rule, ..options() },
                || false,
            )
            .unwrap();
            near(actual.reflection, expected, 1.0, 2e-11);
            assert_eq!(actual.reference_impedance, z1);
            // The second line is matched in its own material. All power not
            // reflected at this lossless interface reaches that termination.
            let transmitted_power = 4.0 * (z1 / z2) / (1.0 + z1 / z2).powi(2);
            let reflected_power = actual.reflection.abs().powi(2);
            assert!((1.0 - reflected_power - transmitted_power).abs() < 3e-11);
        }
    }
}

fn numerical_liquid() -> ThermoviscousFluidState {
    // Explicit numerical coefficients, not a claimed measurement of water.
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap(),
        3900.0,
    )
    .unwrap()
}

#[test]
fn lossy_reflection_preserves_wall_diffusion_and_region_owners() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let liquid = numerical_liquid();
    let geometry = Duct {
        segments: [0.2, 0.3]
            .map(|length| Segment::Cylinder {
                radius: 0.01,
                length,
            })
            .to_vec(),
    };
    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 4e7,
        resistance: 300.0,
    };
    for material in [&gas as &dyn ThermoviscousFluid, &liquid] {
        for boundary in [WallModel::Rigid, WallModel::Pin(&wall)] {
            for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
                let opts = ComposeOptions { rule, ..options() };
                let actual = fluid::input_reflection_load_with_loss(
                    &geometry,
                    material,
                    900.0,
                    LossModel::Bessel,
                    None,
                    boundary,
                    ConeSlices::Auto,
                    opts,
                    || false,
                )
                .unwrap();
                let response = fluid::input_impedance_load_with_loss(
                    &geometry,
                    material,
                    900.0,
                    LossModel::Bessel,
                    None,
                    boundary,
                    ConeSlices::Auto,
                    opts,
                    || false,
                )
                .unwrap();
                let reference = C64::from_re(actual.reference_impedance);
                near(
                    actual.reflection,
                    (response.z_in - reference) / (response.z_in + reference),
                    1.0,
                    2e-12,
                );
                assert!(actual.reflection.abs() < 1.0);
                let regions = regional::input_reflection_load_with_loss(
                    &geometry,
                    &[material, material],
                    900.0,
                    LossModel::Bessel,
                    None,
                    boundary,
                    ConeSlices::Auto,
                    opts,
                    || false,
                )
                .unwrap();
                near(regions.reflection, actual.reflection, 1.0, 2e-12);
            }
        }
    }
}

#[test]
fn material_boundary_reflection_uses_one_coherent_snapshot() {
    struct Once {
        state: ThermoviscousFluidState,
        calls: Cell<usize>,
    }
    impl FluidAcoustics for Once {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            panic!("boundary must use the coherent thermoviscous snapshot")
        }
    }
    impl FluidTransport for Once {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            panic!("boundary must use the coherent thermoviscous snapshot")
        }
    }
    impl ThermoviscousFluid for Once {
        fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            assert_eq!(self.calls.get(), 1);
            Ok(self.state)
        }
    }
    let material = Once {
        state: numerical_liquid(),
        calls: Cell::new(0),
    };
    let geometry = cylinder(0.01, 0.4);
    let boundary = NetworkLoad::Radiation {
        radius: 0.01,
        termination: fs_duct::Termination::UnflangedOpen,
    };
    let actual = fluid::input_reflection_boundary(
        &geometry,
        &material,
        900.0,
        LossModel::Bessel,
        &boundary,
        WallModel::Rigid,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap();
    assert_eq!(material.calls.get(), 1);
    let response = fluid::input_impedance_boundary(
        &geometry,
        &material.state,
        900.0,
        LossModel::Bessel,
        &boundary,
        WallModel::Rigid,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap();
    let reference = C64::from_re(actual.reference_impedance);
    near(
        actual.reflection,
        (response.z_in - reference) / (response.z_in + reference),
        1.0,
        2e-12,
    );
    assert!(actual.reflection.abs() < 1.0);
}

#[test]
fn controlled_reflection_cancels_during_the_terminated_sweep() {
    let blocks = [compose::TwoPort::through(10.0).unwrap(); 3];
    let mut polls = 0;
    let result = compose::compose_chain_reflection_controlled(
        &blocks,
        None,
        ComposeOptions {
            rule: Composition::ImpedanceRecursion,
            ..options()
        },
        || {
            polls += 1;
            polls > blocks.len()
        },
    );
    assert!(matches!(result, Err(DuctError::Cancelled { .. })));
    assert_eq!(polls, blocks.len() + 1);
    let actual = compose::compose_chain_reflection(&blocks, None, options()).unwrap();
    assert_eq!(actual.reflection, C64::ONE);
}
