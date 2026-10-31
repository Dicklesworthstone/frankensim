//! G1/G3: regional materials retain original operators and interface work.
use super::*;
use fs_duct::compose::Composition;
use fs_duct::fluid::regional;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use std::cell::Cell;

#[test]
fn uniform_regions_retain_original_geometry_wall_and_diffusion_response() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.012,
                length: 0.1,
            },
            Segment::Cone {
                inlet_radius: 0.012,
                outlet_radius: 0.016,
                length: 0.2,
            },
            Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(0.016, 0.5, 0.2).unwrap()),
                steps: RiccatiSteps(24),
            },
        ],
    };
    let materials: [&dyn ThermoviscousFluid; 3] = [&gas, &gas, &gas];
    let pin = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 4e7,
        resistance: 300.0,
    };
    for wall in [WallModel::Rigid, WallModel::Pin(&pin)] {
        for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
            let original = fluid::duct_chain_with_loss(
                &geometry,
                &gas,
                700.0,
                loss,
                wall,
                ConeSlices::Fixed(12),
                options(),
                || false,
            )
            .unwrap();
            let bound = regional::duct_chain_with_loss(
                &geometry,
                &materials,
                700.0,
                loss,
                wall,
                ConeSlices::Fixed(12),
                options(),
                || false,
            )
            .unwrap();
            assert_eq!(bound.reference_impedance, original.reference_impedance);
            assert_eq!(bound.min_shear_number, original.min_shear_number);
            assert_eq!(bound.blocks.len(), original.blocks.len());
            for (actual, expected) in bound.blocks.iter().zip(&original.blocks) {
                assert_eq!(actual.scattering, expected.scattering);
                assert_eq!(actual.reference_impedance, expected.reference_impedance);
            }
            for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
                let opts = ComposeOptions { rule, ..options() };
                let load = Some(C64::new(1e6, -2e5));
                let expected = fluid::input_impedance_load_with_loss(
                    &geometry,
                    &gas,
                    700.0,
                    loss,
                    load,
                    wall,
                    ConeSlices::Fixed(12),
                    opts,
                    || false,
                )
                .unwrap();
                let actual = regional::input_impedance_load_with_loss(
                    &geometry,
                    &materials,
                    700.0,
                    loss,
                    load,
                    wall,
                    ConeSlices::Fixed(12),
                    opts,
                    || false,
                )
                .unwrap();
                assert_eq!(actual.z_in, expected.z_in);
                assert_eq!(actual.u_out_over_u_in, expected.u_out_over_u_in);
                assert_eq!(actual.p_out_over_u_in, expected.p_out_over_u_in);
            }
        }
    }
}

// Independent plane-wave ABCD oracle under exp(-i omega t). It uses only the
// original rho/c values and trig, never the native wave/scattering builders.
fn analytic_line(material: &dyn FluidAcoustics, radius: f64, length: f64, omega: f64) -> [C64; 4] {
    let state = material.acoustic_state().unwrap();
    let z = state.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * radius * radius);
    let phase = omega * length / state.sound_speed_m_s();
    [
        C64::from_re(phase.cos()),
        C64::new(0.0, -z * phase.sin()),
        C64::new(0.0, -phase.sin() / z),
        C64::from_re(phase.cos()),
    ]
}

#[test]
fn if97_and_mixed_source_interfaces_preserve_pressure_flow_and_power() {
    let cold = water(300.0, 2_000_000.0);
    let warm = water(350.0, 2_000_000.0);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 2_000_000.0).unwrap();
    let (radius, omega) = (0.012, 900.0);
    let geometry = Duct {
        segments: [0.2, 0.37]
            .map(|length| Segment::Cylinder { radius, length })
            .to_vec(),
    };
    let combinations: [[&dyn FluidAcoustics; 2]; 2] = [[&cold, &warm], [&gas, &warm]];
    let mut impedances = Vec::new();
    for materials in combinations {
        let left = analytic_line(materials[0], radius, 0.2, omega);
        let right = analytic_line(materials[1], radius, 0.37, omega);
        let product = [
            left[0] * right[0] + left[1] * right[2],
            left[0] * right[1] + left[1] * right[3],
            left[2] * right[0] + left[3] * right[2],
            left[2] * right[1] + left[3] * right[3],
        ];
        let z = materials[1]
            .acoustic_state()
            .unwrap()
            .characteristic_impedance_pa_s_m()
            / (core::f64::consts::PI * radius * radius);
        let load = C64::from_re(z);
        let flow = (product[2] * load + product[3]).recip();
        let input = (product[0] * load + product[1]) * flow;
        for rule in [Composition::ImpedanceRecursion, Composition::Scattering] {
            let actual = regional::input_impedance_load(
                &geometry,
                &materials,
                omega,
                Some(load),
                WallModel::Rigid,
                ConeSlices::Auto,
                ComposeOptions { rule, ..options() },
                || false,
            )
            .unwrap();
            near(actual.z_in, input, input.abs(), 2e-11);
            near(actual.u_out_over_u_in, flow, 1.0, 2e-11);
            near(actual.p_out_over_u_in, load * flow, load.abs(), 2e-11);
            let output_power = z * actual.u_out_over_u_in.abs().powi(2);
            assert!((actual.z_in.re - output_power).abs() < 2e-10 * actual.z_in.abs());
            impedances.push(actual.z_in);
        }
    }
    assert!((impedances[0] - impedances[2]).abs() > 0.1 * impedances[0].abs());
}

#[test]
fn regional_sources_are_evaluated_once_and_failed_admission_publishes_nothing() {
    struct Once {
        state: FluidAcousticState,
        calls: Cell<usize>,
    }
    impl FluidAcoustics for Once {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            let calls = self.calls.get();
            self.calls.set(calls + 1);
            if calls == 0 {
                Ok(self.state)
            } else {
                Err(MaterialError::State {
                    what: "source queried twice".into(),
                })
            }
        }
    }
    let materials = [1000.0, 950.0].map(|rho| Once {
        state: FluidAcousticState::try_new(rho, 1500.0).unwrap(),
        calls: Cell::new(0),
    });
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.2
            };
            2
        ],
    };
    assert!(
        regional::duct_chain(
            &geometry,
            &materials[..1],
            700.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        regional::duct_chain(
            &geometry,
            &materials,
            700.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        }))
    ));
    assert!(materials.iter().all(|m| m.calls.get() == 0));
    regional::input_impedance_load(
        &geometry,
        &materials,
        700.0,
        Some(C64::from_re(2e9)),
        WallModel::Rigid,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap();
    assert!(materials.iter().all(|m| m.calls.get() == 1));
    assert!(matches!(
        regional::duct_chain(
            &geometry,
            &materials,
            700.0,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false
        ),
        Err(FluidDuctError::Material(_))
    ));
}

#[test]
fn regional_diffusion_respects_one_global_slice_budget() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let geometry = Duct {
        segments: vec![
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.012,
                length: 0.2,
            },
            Segment::Cone {
                inlet_radius: 0.012,
                outlet_radius: 0.014,
                length: 0.2,
            },
        ],
    };
    let solve = |maximum_blocks| {
        regional::duct_chain_with_loss(
            &geometry,
            &[gas, gas],
            700.0,
            LossModel::Bessel,
            WallModel::Rigid,
            ConeSlices::Fixed(8),
            ComposeOptions {
                maximum_blocks,
                ..options()
            },
            || false,
        )
    };
    assert!(
        matches!(solve(15), Err(FluidDuctError::Duct(DuctError::BadParameter { what }))
        if what.contains("block work cap"))
    );
    assert_eq!(solve(16).unwrap().blocks.len(), 16);
    let mut stale = gas;
    stale.sound_speed *= 1.01;
    assert!(matches!(
        regional::duct_chain_with_loss(
            &geometry,
            &[gas, stale],
            700.0,
            LossModel::Bessel,
            WallModel::Rigid,
            ConeSlices::Fixed(8),
            options(),
            || false,
        ),
        Err(FluidDuctError::Material(_))
    ));
}
