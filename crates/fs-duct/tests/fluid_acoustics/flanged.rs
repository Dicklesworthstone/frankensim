//! G1/G3: the original radiation operator serves gas and liquid materials.

use super::*;
use fs_duct::flanged::{
    flanged_modal_impedance, flanged_modal_impedance_from_material,
    flanged_pipe_impedance_from_material,
};
use fs_duct::network::NetworkLoad;
use fs_duct::{FLANGED_MAX_KA, FLANGED_MAX_MODES, FLANGED_MODES, Termination};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

#[test]
fn gas_radiation_bits_and_liquid_piston_asymptote_retain_the_original_owner() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let liquid = water(300.0, 3_000_000.0);
    let acoustic = liquid.acoustic_state().unwrap();
    let radius = 0.012;
    let area = core::f64::consts::PI * radius * radius;
    for modes in [1, 4, FLANGED_MODES] {
        for ka in [0.01, 0.8, 4.2] {
            let omega = ka * gas.sound_speed / radius;
            let legacy = fs_duct::flanged_pipe_impedance(&gas, radius, omega, modes).unwrap();
            let shared = flanged_pipe_impedance_from_material(&gas, radius, omega, modes).unwrap();
            assert_eq!(
                (legacy.re.to_bits(), legacy.im.to_bits()),
                (shared.re.to_bits(), shared.im.to_bits())
            );
            assert_eq!(
                flanged_modal_impedance(&gas, radius, omega, modes).unwrap(),
                flanged_modal_impedance_from_material(&gas, radius, omega, modes).unwrap()
            );
            let water_omega = omega * acoustic.sound_speed_m_s() / gas.sound_speed;
            let z =
                flanged_pipe_impedance_from_material(&liquid, radius, water_omega, modes).unwrap();
            let ratio =
                acoustic.characteristic_impedance_pa_s_m() / (gas.density * gas.sound_speed);
            near(z, legacy.scale(ratio), 1.0, 3e-11);
            assert!(z.re >= 0.0);
        }
    }
    // Independent leading terms of the infinite-baffle circular piston:
    // R/(rho*c/A) = (ka)^2/2; X/(rho*c/A) = -8*ka/(3*pi).
    let ka = 1e-3;
    let omega = ka * acoustic.sound_speed_m_s() / radius;
    let z = flanged_pipe_impedance_from_material(&liquid, radius, omega, 1).unwrap();
    let z0 = acoustic.characteristic_impedance_pa_s_m() / area;
    assert!((z.re / (0.5 * ka * ka * z0) - 1.0).abs() < 3e-7);
    assert!((z.im / (-8.0 * ka * z0 / (3.0 * core::f64::consts::PI)) - 1.0).abs() < 3e-7);
}

#[test]
fn liquid_multimodal_radiation_preserves_reciprocity_and_positive_power() {
    let liquid = water(330.0, 2_000_000.0);
    let acoustic = liquid.acoustic_state().unwrap();
    let radius = 0.03;
    let modes = 4;
    let flows = [
        C64::new(0.4, -0.2),
        C64::new(-0.3, 0.7),
        C64::new(0.5, 0.1),
        C64::new(-0.1, -0.4),
    ];
    for ka in [0.5, 4.5, 9.0] {
        let z = flanged_modal_impedance_from_material(
            &liquid,
            radius,
            ka * acoustic.sound_speed_m_s() / radius,
            modes,
        )
        .unwrap();
        let mut power = 0.0;
        for i in 0..modes {
            let mut pressure = C64::ZERO;
            for j in 0..modes {
                assert_eq!(z[i * modes + j], z[j * modes + i]);
                pressure = pressure + z[i * modes + j] * flows[j];
            }
            power += 0.5 * (flows[i].conj() * pressure).re;
        }
        assert!(power > 0.0, "ka={ka}, radiated power={power}");
    }
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    samples: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("the boundary must reuse its single complete material snapshot")
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("the boundary must reuse its single complete material snapshot")
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn radiating_liquid_pipe_uses_one_material_state_and_closes_pressure_flow_power() {
    // Authored transport fixture; the source EOS independently supplies rho/c.
    let acoustic = water(300.0, 3_000_000.0).acoustic_state().unwrap();
    let transport =
        FluidTransportState::try_new(acoustic.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let state = ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap();
    let source = CountedFluid {
        state,
        samples: Cell::new(0),
    };
    let (radius, length) = (0.012, 0.17);
    let duct = cylinder(radius, length);
    let omega = 0.6 * acoustic.sound_speed_m_s() / radius;
    let load = NetworkLoad::Radiation {
        radius,
        termination: Termination::FlangedOpen,
    };
    let z_load =
        flanged_pipe_impedance_from_material(&state, radius, omega, FLANGED_MODES).unwrap();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let result = fluid::input_impedance_boundary(
            &duct,
            &source,
            omega,
            loss,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
        .unwrap();
        let expected = fluid::input_impedance_load_with_loss(
            &duct,
            &state,
            omega,
            loss,
            Some(z_load),
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
        .unwrap();
        assert_eq!(result.z_in, expected.z_in);
        near(
            result.p_out_over_u_in,
            z_load * result.u_out_over_u_in,
            1.0,
            2e-14,
        );
        let transmitted = z_load.re * result.u_out_over_u_in.abs().powi(2);
        if loss == LossModel::Lossless {
            let phase = omega * length / acoustic.sound_speed_m_s();
            let zc = acoustic.characteristic_impedance_pa_s_m()
                / (core::f64::consts::PI * radius * radius);
            let c = C64::from_re(phase.cos());
            let is = C64::new(0.0, -phase.sin());
            let analytic = (z_load * c + is.scale(zc)) / (c + is * z_load.scale(1.0 / zc));
            near(result.z_in, analytic, zc, 3e-14);
            assert!((result.z_in.re / transmitted - 1.0).abs() < 3e-14);
        } else {
            assert!(
                result.z_in.re > transmitted,
                "wall diffusion must dissipate additional power"
            );
        }
    }
    assert_eq!(source.samples.get(), 2);
    assert!(matches!(
        fluid::input_impedance_boundary(
            &duct,
            &source,
            omega,
            LossModel::Lossless,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        }))
    ));
    assert_eq!(source.samples.get(), 2);
}

#[test]
fn original_source_refusals_and_radiation_bounds_survive_generic_admission() {
    struct Refused;
    impl FluidAcoustics for Refused {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Err(MaterialError::State {
                what: "original fluid phase refused".into(),
            })
        }
    }
    for error in [
        flanged_pipe_impedance_from_material(&Refused, 0.01, 1000.0, 4).unwrap_err(),
        flanged_modal_impedance_from_material(&Refused, 0.01, 1000.0, 4).unwrap_err(),
        fs_duct::compact_radiation_impedance_from_material(
            Termination::UnflangedOpen,
            &Refused,
            0.01,
            1000.0,
        )
        .unwrap_err(),
    ] {
        assert_eq!(
            error,
            FluidDuctError::Material(MaterialError::State {
                what: "original fluid phase refused".into(),
            })
        );
    }
    let liquid = water(300.0, 3_000_000.0);
    let sound = liquid.acoustic_state().unwrap().sound_speed_m_s();
    let omega = (FLANGED_MAX_KA + 1.0) * sound / 0.01;
    assert!(matches!(
        flanged_pipe_impedance_from_material(&liquid, 0.01, omega, 4),
        Err(FluidDuctError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));
    for modes in [0, FLANGED_MAX_MODES + 1] {
        assert!(flanged_modal_impedance_from_material(&liquid, 0.01, 1000.0, modes).is_err());
    }
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    gas.sound_speed *= 1.1;
    assert!(matches!(
        flanged_pipe_impedance_from_material(&gas, 0.01, 1000.0, 4),
        Err(FluidDuctError::Material(_))
    ));
    assert!(fs_duct::flanged_pipe_impedance(&gas, 0.01, 1000.0, 4).is_err());
}

#[test]
fn original_unflanged_models_share_material_coefficients_and_independent_liquid_references() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let liquid = water(300.0, 3_000_000.0);
    let acoustic = liquid.acoustic_state().unwrap();
    let radius = 0.015;
    let ka = 0.5;
    let gas_omega = ka * gas.sound_speed / radius;
    let omega = ka * acoustic.sound_speed_m_s() / radius;
    let z0 = acoustic.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * radius * radius);
    for termination in [
        Termination::UnflangedOpen,
        Termination::UnflangedThinWall {
            quadrature: fs_duct::LsQuadrature { nodes: 256 },
        },
        Termination::FlangedOpen,
    ] {
        let original =
            fs_duct::compact_radiation_impedance(termination, &gas, radius, gas_omega).unwrap();
        let shared = fs_duct::compact_radiation_impedance_from_material(
            termination,
            &gas,
            radius,
            gas_omega,
        )
        .unwrap();
        assert_eq!(
            (original.re.to_bits(), original.im.to_bits()),
            (shared.re.to_bits(), shared.im.to_bits())
        );
        let z =
            fs_duct::compact_radiation_impedance_from_material(termination, &liquid, radius, omega)
                .unwrap();
        assert!(z.re > 0.0 && z.im < 0.0);
        if termination == Termination::UnflangedOpen {
            near(
                z,
                C64::new(ka * ka / 4.0, -0.6133 * ka).scale(z0),
                z0,
                3e-15,
            );
        } else if matches!(termination, Termination::UnflangedThinWall { .. }) {
            // Independent 35-digit integration of original LS kernels, retained
            // by the numerical owner's reference table: ka=.5, |R| and l/a.
            let reflection = (z - C64::from_re(z0)) / (z + C64::from_re(z0));
            let r_abs = 0.8964413976209895;
            let end_correction = 0.5818364010713508;
            let phase = 2.0 * ka * end_correction;
            let reference = C64::new(-r_abs * phase.cos(), -r_abs * phase.sin());
            near(reflection, reference, 1.0, 1e-3);
        }
    }
}

#[test]
fn every_native_radiation_boundary_reuses_one_fluid_snapshot_and_original_limits() {
    let acoustic = water(300.0, 3_000_000.0).acoustic_state().unwrap();
    let transport =
        FluidTransportState::try_new(acoustic.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let source = CountedFluid {
        state: ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap(),
        samples: Cell::new(0),
    };
    let radius = 0.01;
    let duct = cylinder(radius, 0.2);
    let omega = 0.5 * acoustic.sound_speed_m_s() / radius;
    let thin_wall = Termination::UnflangedThinWall {
        quadrature: fs_duct::LsQuadrature { nodes: 128 },
    };
    for termination in [
        Termination::UnflangedOpen,
        thin_wall,
        Termination::FlangedOpen,
    ] {
        let load = NetworkLoad::Radiation {
            radius,
            termination,
        };
        let result = fluid::input_impedance_boundary(
            &duct,
            &source,
            omega,
            LossModel::Lossless,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
        .unwrap();
        let z = fs_duct::compact_radiation_impedance_from_material(
            termination,
            &acoustic,
            radius,
            omega,
        )
        .unwrap();
        near(
            result.p_out_over_u_in,
            z * result.u_out_over_u_in,
            1.0,
            3e-14,
        );
        let power = z.re * result.u_out_over_u_in.abs().powi(2);
        assert!((result.z_in.re / power - 1.0).abs() < 3e-14);
    }
    assert_eq!(source.samples.get(), 3);
    for (load, impedance) in [
        (NetworkLoad::Closed, None),
        (NetworkLoad::PressureRelease, Some(C64::ZERO)),
    ] {
        let result = fluid::input_impedance_boundary(
            &duct,
            &source,
            omega,
            LossModel::Lossless,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
        .unwrap();
        let expected = solve(&duct, &acoustic, omega, impedance, WallModel::Rigid);
        assert_eq!(result.z_in, expected.z_in);
    }
    for termination in [Termination::Closed, Termination::IdealOpen] {
        assert!(
            fs_duct::compact_radiation_impedance_from_material(
                termination,
                &acoustic,
                radius,
                omega
            )
            .is_err()
        );
    }
    let fit_above_band = 1.01 * acoustic.sound_speed_m_s() / radius;
    assert!(matches!(
        fs_duct::compact_radiation_impedance_from_material(
            Termination::UnflangedOpen,
            &acoustic,
            radius,
            fit_above_band
        ),
        Err(FluidDuctError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));
    let above_cutoff = 4.0 * acoustic.sound_speed_m_s() / radius;
    assert!(matches!(
        fs_duct::compact_radiation_impedance_from_material(
            thin_wall,
            &acoustic,
            radius,
            above_cutoff
        ),
        Err(FluidDuctError::Duct(
            DuctError::LevineSchwingerCutoff { .. }
        ))
    ));
    let bad_rule = Termination::UnflangedThinWall {
        quadrature: fs_duct::LsQuadrature { nodes: 15 },
    };
    assert!(
        fs_duct::compact_radiation_impedance_from_material(bad_rule, &acoustic, radius, omega)
            .is_err()
    );
    for (r, w) in [(0.0, omega), (radius, f64::NAN), (radius, -omega)] {
        assert!(
            fs_duct::compact_radiation_impedance_from_material(
                Termination::UnflangedOpen,
                &acoustic,
                r,
                w
            )
            .is_err()
        );
    }
}
