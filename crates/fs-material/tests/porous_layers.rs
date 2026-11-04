//! Physics controls for the source-bound multilayer reflection owner.
use core::cell::Cell;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{
    AngularPowerBudget, AngularPowerError, DiffuseIncidence, IncidentPower, Layer,
    LayerAbsorptionBudget, LayerAbsorptionError,
    MAX_ANGULAR_POWER_SAMPLES, MAX_LAYERS, RigidBackedLayers, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;

fn fluid() -> ThermoviscousFluidState {
    // Explicit analytic fixture, not measured ambient data.
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}
fn foam() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn single_layer_reproduces_original_owner_in_actual_gas_states() {
    let law = foam();
    for temperature in [293.15, 323.15] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap();
        for frequency in [150.0, 500.0, 1000.0, 1600.0] {
            let omega = core::f64::consts::TAU * frequency;
            let stack = RigidBackedLayers::new(
                &gas,
                omega,
                &[Layer::Porous {
                    law: &law,
                    thickness_m: 0.06,
                }],
            )
            .unwrap();
            let sample = law.effective_fluid(&gas, omega).unwrap();
            assert_eq!(stack.material_response(0), Some(&sample));
            let zs = sample.surface_impedance_rigid(0.06).unwrap();
            let z0 = C64::from_re(stack.fluid().acoustic().characteristic_impedance_pa_s_m());
            let response = stack.reflection(1.0).unwrap();
            near(response.pressure_reflection, (zs - z0) / (zs + z0), 3e-14);
            near(
                response.normalized_surface_admittance().unwrap(),
                z0 / zs,
                3e-14,
            );
        }
    }
}

// Independent conventional complex transfer matrices. These intentionally use
// growing cos/sin, only in a well-conditioned oracle, not the production E form.
fn cos_sin(z: C64) -> (C64, C64) {
    (
        C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh()),
        C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
    )
}
fn matrix_response(stack: &RigidBackedLayers<'_>, cosine: f64, load: C64) -> (C64, C64) {
    let z0 = stack.fluid().acoustic().characteristic_impedance_pa_s_m();
    let kt = stack.exterior_wavenumber_rad_m() * (1.0 - cosine * cosine).sqrt();
    let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    for i in 0..stack.layer_count() {
        let (kn, normal_inertia, d) = match stack.layer(i).unwrap() {
            Layer::FluidGap { thickness_m } => (
                C64::from_re(stack.exterior_wavenumber_rad_m() * cosine),
                C64::from_re(stack.exterior_wavenumber_rad_m() * cosine),
                thickness_m,
            ),
            Layer::FluidLayer { thickness_m, .. } => {
                let state = stack.fluid_response(i).unwrap();
                let k = stack.angular_frequency_rad_s() / state.sound_speed_m_s();
                let normal_squared = k * k - kt * kt;
                let kn = if normal_squared >= 0.0 {
                    C64::from_re(normal_squared.sqrt())
                } else {
                    C64::new(0.0, -(-normal_squared).sqrt())
                };
                (
                    kn,
                    C64::from_re(
                        stack.angular_frequency_rad_s() * state.density_kg_m3() * cosine / z0,
                    ),
                    thickness_m,
                )
            }
            Layer::Porous { thickness_m, .. } => {
                let (zc, k) = stack.material_response(i).unwrap().characteristic();
                let mut kn = (k * k - C64::from_re(kt * kt)).sqrt();
                if kn.im > 0.0 {
                    kn = -kn;
                }
                (kn, (zc * k).scale(cosine / z0), thickness_m)
            }
        };
        let (cos, sin) = cos_sin(kn.scale(d));
        let j_sin = C64::new(0.0, 1.0) * sin;
        let layer = if kn == C64::ZERO {
            [
                [C64::ONE, C64::new(0.0, d) * normal_inertia],
                [C64::ZERO, C64::ONE],
            ]
        } else {
            let ratio = normal_inertia / kn;
            [[cos, ratio * j_sin], [j_sin / ratio, cos]]
        };
        let mut product = [[C64::ZERO; 2]; 2];
        for row in 0..2 {
            for col in 0..2 {
                product[row][col] = m[row][0] * layer[0][col] + m[row][1] * layer[1][col];
            }
        }
        m = product;
    }
    let pressure = m[0][0] + m[0][1] * load;
    let velocity = m[1][0] + m[1][1] * load;
    (
        (pressure - velocity) / (pressure + velocity),
        C64::from_re(2.0) / (pressure + velocity),
    )
}
fn matrix_reflection(stack: &RigidBackedLayers<'_>, cosine: f64) -> C64 {
    matrix_response(stack, cosine, C64::ZERO).0
}

#[test]
fn ordered_porous_gap_stack_matches_independent_transfer_matrices() {
    let fluid = fluid();
    let first = foam();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &first,
            thickness_m: 0.03,
        },
        Layer::FluidGap { thickness_m: 0.05 },
        Layer::Porous {
            law: &second,
            thickness_m: 0.02,
        },
    ];
    for frequency in [150.0, 500.0, 1000.0, 1600.0] {
        let omega = core::f64::consts::TAU * frequency;
        let stack = RigidBackedLayers::new(&fluid, omega, &layers).unwrap();
        let reverse =
            RigidBackedLayers::new(&fluid, omega, &[layers[2], layers[1], layers[0]]).unwrap();
        for cosine in [1.0, 0.85, 0.5, 0.05] {
            let actual = stack.reflection(cosine).unwrap();
            let expected = matrix_reflection(&stack, cosine);
            near(actual.pressure_reflection, expected, 4e-14);
            assert!((actual.absorbed_fraction - (1.0 - expected.norm_sq())).abs() < 4e-14);
            assert!(actual.absorbed_fraction >= 0.0);
        }
        assert!(
            (stack.reflection(1.0).unwrap().pressure_reflection
                - reverse.reflection(1.0).unwrap().pressure_reflection)
                .abs()
                > 1e-4
        );
    }
}

#[test]
fn subdivision_and_zero_thickness_preserve_coherent_reflection() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU * 800.0;
    let whole = RigidBackedLayers::new(
        &fluid,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.06,
        }],
    )
    .unwrap();
    let split = RigidBackedLayers::new(
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.01,
            },
            Layer::FluidGap { thickness_m: 0.0 },
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
            Layer::Porous {
                law: &law,
                thickness_m: 0.03,
            },
        ],
    )
    .unwrap();
    for cosine in [1.0, 0.5, 0.05] {
        near(
            split.reflection(cosine).unwrap().pressure_reflection,
            whole.reflection(cosine).unwrap().pressure_reflection,
            4e-14,
        );
    }
}

#[test]
fn lossless_gaps_preserve_round_trip_phase_poles_and_zero_work() {
    let fluid = fluid();
    // k0=1 rad/m makes quarter/half-wave controls explicit.
    let omega = fluid.acoustic().sound_speed_m_s();
    for d in [
        0.0,
        0.1,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI,
    ] {
        let stack =
            RigidBackedLayers::new(&fluid, omega, &[Layer::FluidGap { thickness_m: d }]).unwrap();
        for cosine in [1.0, 0.5, 1e-200] {
            let response = stack.reflection(cosine).unwrap();
            let phase = 2.0 * cosine * d;
            near(
                response.pressure_reflection,
                C64::new(phase.cos(), -phase.sin()),
                3e-15,
            );
            assert_eq!(response.absorbed_fraction, 0.0);
            if let Ok(beta) = response.normalized_surface_admittance() {
                assert_eq!(beta.re, 0.0);
            }
        }
    }
    let empty = RigidBackedLayers::new(&fluid, omega, &[]).unwrap();
    assert_eq!(
        empty.reflection(0.25).unwrap().pressure_reflection,
        C64::ONE
    );
    assert_eq!(
        empty
            .reflection(0.25)
            .unwrap()
            .normalized_surface_admittance()
            .unwrap(),
        C64::ZERO
    );
}

#[test]
fn attenuating_thick_limit_and_extreme_grazing_keep_finite_positive_work() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU * 500.0;
    // Conventional cosh matrices overflow; the attenuating half-space limit is finite.
    let thick = RigidBackedLayers::new(
        &fluid,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: 1e308,
        }],
    )
    .unwrap();
    let (zc, k) = thick.material_response(0).unwrap().characteristic();
    let z0 = fluid.acoustic().characteristic_impedance_pa_s_m();
    for cosine in [1.0, 0.5, 0.01] {
        let kt = thick.exterior_wavenumber_rad_m() * (1.0 - cosine * cosine).sqrt();
        let mut kn = (k * k - C64::from_re(kt * kt)).sqrt();
        if kn.im > 0.0 {
            kn = -kn;
        }
        let ratio = (zc * k / kn).scale(cosine / z0);
        let actual = thick.reflection(cosine).unwrap();
        near(
            actual.pressure_reflection,
            (ratio - C64::ONE) / (ratio + C64::ONE),
            3e-14,
        );
        assert!(actual.absorbed_fraction > 0.0);
    }
    let stack = RigidBackedLayers::new(
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
            Layer::Porous {
                law: &law,
                thickness_m: 0.04,
            },
            Layer::FluidGap { thickness_m: 0.1 },
        ],
    )
    .unwrap();
    let baseline = stack.reflection(1e-100).unwrap();
    for cosine in [1e-200, 1e-300] {
        let response = stack.reflection(cosine).unwrap();
        assert!(response.absorbed_fraction > 0.0);
        assert!(
            (response.absorbed_fraction / cosine - baseline.absorbed_fraction / 1e-100).abs()
                < 2e-13
        );
        assert!(
            (response.pressure_reflection.im / cosine - baseline.pressure_reflection.im / 1e-100)
                .abs()
                < 2e-13
        );
    }
}

#[test]
fn weak_absorption_survives_unit_reflected_magnitude_roundoff() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU * 500.0;
    let d = 1e-100;
    let stack = RigidBackedLayers::new(
        &fluid,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: d,
        }],
    )
    .unwrap();
    let response = stack.reflection(1.0).unwrap();
    let bulk = stack
        .material_response(0)
        .unwrap()
        .dynamic_bulk_modulus_pa();
    let beta = C64::new(
        0.0,
        omega * d * fluid.acoustic().characteristic_impedance_pa_s_m(),
    ) / bulk;
    assert_eq!(response.pressure_reflection.norm_sq(), 1.0);
    assert!(response.absorbed_fraction > 0.0);
    assert!((response.absorbed_fraction / (4.0 * beta.re) - 1.0).abs() < 2e-13);
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    samples: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn all_layers_share_one_actual_snapshot_and_keep_source_domain() {
    let fluid = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let law = foam();
    let layers = [
        Layer::Porous {
            law: &law,
            thickness_m: 0.03,
        },
        Layer::FluidGap { thickness_m: 0.05 },
        Layer::Porous {
            law: &law,
            thickness_m: 0.01,
        },
    ];
    let stack = RigidBackedLayers::new(&fluid, core::f64::consts::TAU * 500.0, &layers).unwrap();
    for cosine in [1.0, 0.5, 0.01] {
        stack.reflection(cosine).unwrap();
    }
    assert_eq!(fluid.samples.get(), 1);
    assert_eq!(stack.material_response(0).unwrap().fluid(), &fluid.state);
    assert_eq!(stack.material_response(2).unwrap().fluid(), &fluid.state);
    assert!(stack.material_response(1).is_none());
    assert!(
        matches!(stack.layer(0), Some(Layer::Porous { law: source, .. }) if core::ptr::eq(source, &law))
    );
    assert!(stack.layer(3).is_none());
    for cosine in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(stack.reflection(cosine).is_err());
    }
    assert!(
        RigidBackedLayers::new(
            &fluid,
            core::f64::consts::TAU * 149.0,
            &[Layer::Porous {
                law: &law,
                thickness_m: 0.0
            }]
        )
        .is_err()
    );
    assert!(RigidBackedLayers::new(&fluid, 1.0, &[Layer::FluidGap { thickness_m: -0.1 }]).is_err());
    assert!(
        RigidBackedLayers::new(
            &fluid,
            1.0,
            &[Layer::FluidGap { thickness_m: 0.0 }; MAX_LAYERS + 1]
        )
        .is_err()
    );
    for omega in [0.0, f64::NAN, f64::INFINITY] {
        assert!(RigidBackedLayers::new(&fluid, omega, &[]).is_err());
    }
}

#[test]
fn receiving_interface_obeys_pressure_velocity_refraction_and_critical_limits() {
    let fluid = fluid();
    let receiver = FluidAcousticState::try_new(1.2, 425.0).unwrap();
    let stack =
        TransmittingLayers::new(&fluid, &receiver, core::f64::consts::TAU * 500.0, &[]).unwrap();
    for cosine in [1.0, 0.8, 0.6, 0.5, 0.001] {
        let response = stack.scattering(cosine).unwrap();
        let q = response.receiving_normal_wavenumber_rad_m
            / C64::from_re(stack.prepared_layers().exterior_wavenumber_rad_m());
        let load = q.scale(1.0 / cosine);
        near(
            response.pressure_transmission,
            C64::ONE + response.pressure_reflection,
            2e-14,
        );
        near(
            C64::ONE - response.pressure_reflection,
            load * response.pressure_transmission,
            2e-14,
        );
        assert_eq!(response.absorbed_fraction, 0.0);
        assert!((response.reflected_fraction + response.transmitted_fraction - 1.0).abs() < 2e-14);
        // Binary64 0.6 lies just below the exact critical cosine 3/5.
        // Its nonzero evanescent normal component must not be rounded away.
        if cosine <= 0.6 {
            assert!(response.receiving_normal_wavenumber_rad_m.im < 0.0);
            assert_eq!(response.receiving_normal_wavenumber_rad_m.re, 0.0);
            assert_eq!(response.transmitted_fraction, 0.0);
            assert!(response.pressure_reflection.im > 0.0);
        } else {
            let expected_q = (0.8_f64 * 0.8 - (1.0 - cosine * cosine)).sqrt();
            assert!((q.re - expected_q).abs() < 1e-15);
            assert_eq!(q.im, 0.0);
        }
    }
}

#[test]
fn original_sound_speed_contrast_preserves_grazing_transmission_and_layer_phase() {
    let exterior = fluid();
    let faster = FluidAcousticState::try_new(1.2, f64::from_bits(340.0_f64.to_bits() + 1)).unwrap();
    let cosine = 2.107_342_425_544_701_4e-8;
    // Exact rational arithmetic on the original binary64 speeds and cosine,
    // followed by an 85-decimal-digit sqrt, gives this positive k_n/k0.
    // Squaring the rounded speed ratio instead misclassifies it as evanescent.
    let q = 1.047_454_809_181_105_1e-8;
    let omega = 340.0; // k0=1 rad/m; thickness below fixes the normal phase.
    let interface = TransmittingLayers::new(&exterior, &faster, omega, &[]).unwrap();
    let result = interface.scattering(cosine).unwrap();
    let load = q / cosine;
    let expected_r = (1.0 - load) / (1.0 + load);
    let expected_t = 2.0 / (1.0 + load);
    near(result.pressure_reflection, C64::from_re(expected_r), 3e-15);
    near(
        result.pressure_transmission,
        C64::from_re(expected_t),
        3e-15,
    );
    assert!((result.transmitted_fraction - load * expected_t * expected_t).abs() < 3e-15);
    assert!(result.transmitted_fraction > 0.88);
    assert_eq!(result.receiving_normal_wavenumber_rad_m.im, 0.0);
    assert_eq!(result.absorbed_fraction, 0.0);

    let phase = 0.75_f64;
    let thickness = phase / q;
    let whole_layers = [Layer::FluidLayer {
        fluid: &faster,
        thickness_m: thickness,
    }];
    let split_layers = [Layer::FluidLayer {
        fluid: &faster,
        thickness_m: thickness * 0.5,
    }; 2];
    let normal_impedance_ratio = cosine / q;
    // Independent pressure/velocity transfer with the original normal mode.
    let p = C64::new(phase.cos(), normal_impedance_ratio * phase.sin());
    let v = C64::new(phase.cos(), phase.sin() / normal_impedance_ratio);
    let expected_r = (p - v) / (p + v);
    let expected_t = C64::from_re(2.0) / (p + v);
    let expected_rigid = (C64::from_re(phase.cos())
        - C64::new(0.0, phase.sin() / normal_impedance_ratio))
        / (C64::from_re(phase.cos()) + C64::new(0.0, phase.sin() / normal_impedance_ratio));
    for layers in [&whole_layers[..], &split_layers[..]] {
        let stack = TransmittingLayers::new(&exterior, &exterior, omega, layers).unwrap();
        let response = stack.scattering(cosine).unwrap();
        near(response.pressure_reflection, expected_r, 5e-15);
        near(response.pressure_transmission, expected_t, 5e-15);
        assert_eq!(response.absorbed_fraction, 0.0);
        near(
            stack
                .prepared_layers()
                .reflection(cosine)
                .unwrap()
                .pressure_reflection,
            expected_rigid,
            5e-15,
        );
    }
}

#[test]
fn critical_side_from_original_speeds_keeps_evanescent_pressure_and_zero_power() {
    let exterior = fluid();
    let faster = FluidAcousticState::try_new(1.2, 500.0).unwrap();
    let cosine = 0.733_212_111_192_934_4;
    // Independent exact-input rational/sqrt oracle. The old rounded-ratio
    // discriminant was positive, assigning spurious outgoing normal power.
    let q = C64::new(0.0, -3.952_211_632_119_969_6e-9);
    let load = q.scale(1.0 / cosine);
    let stack = TransmittingLayers::new(&exterior, &faster, 340.0, &[]).unwrap();
    let result = stack.scattering(cosine).unwrap();
    near(result.receiving_normal_wavenumber_rad_m, q, 5e-23);
    near(
        result.pressure_reflection,
        (C64::ONE - load) / (C64::ONE + load),
        3e-15,
    );
    near(
        result.pressure_transmission,
        C64::from_re(2.0) / (C64::ONE + load),
        3e-15,
    );
    assert_eq!(result.transmitted_fraction, 0.0);
    assert_eq!(result.absorbed_fraction, 0.0);
    assert!(result.pressure_reflection.im > 0.0);
    assert!(
        stack
            .scattering(cosine.next_up())
            .unwrap()
            .transmitted_fraction
            > 0.0
    );
}

#[test]
fn heterogeneous_fluid_and_porous_layers_match_independent_refracted_matrices() {
    let exterior = fluid();
    let faster = FluidAcousticState::try_new(2.3, 425.0).unwrap();
    let slower = FluidAcousticState::try_new(0.6, 280.0).unwrap();
    let receiver = FluidAcousticState::try_new(1.8, 380.0).unwrap();
    let foam = foam();
    let layers = [
        Layer::FluidLayer {
            fluid: &faster,
            thickness_m: 0.025,
        },
        Layer::Porous {
            law: &foam,
            thickness_m: 0.018,
        },
        Layer::FluidGap { thickness_m: 0.007 },
        Layer::FluidLayer {
            fluid: &slower,
            thickness_m: 0.04,
        },
    ];
    for frequency in [150.0, 500.0, 1600.0] {
        let stack = TransmittingLayers::new(
            &exterior,
            &receiver,
            core::f64::consts::TAU * frequency,
            &layers,
        )
        .unwrap();
        assert_eq!(stack.prepared_layers().fluid_response(0), Some(&faster));
        assert_eq!(stack.prepared_layers().fluid_response(3), Some(&slower));
        assert!(stack.prepared_layers().fluid_response(1).is_none());
        for cosine in [1.0, 0.8, 0.5, 0.05] {
            let q2 = (340.0_f64 / 380.0).powi(2) - (1.0 - cosine * cosine);
            let q = if q2 >= 0.0 {
                C64::from_re(q2.sqrt())
            } else {
                C64::new(0.0, -(-q2).sqrt())
            };
            let load = q.scale((1.2 / 1.8) / cosine);
            let (expected_r, expected_t) = matrix_response(stack.prepared_layers(), cosine, load);
            let actual = stack.scattering(cosine).unwrap();
            near(actual.pressure_reflection, expected_r, 8e-14);
            near(actual.pressure_transmission, expected_t, 8e-14);
            assert!((actual.transmitted_fraction - load.re * expected_t.norm_sq()).abs() < 8e-14);
            near(
                stack
                    .prepared_layers()
                    .reflection(cosine)
                    .unwrap()
                    .pressure_reflection,
                matrix_reflection(stack.prepared_layers(), cosine),
                8e-14,
            );
            assert!(
                (actual.reflected_fraction
                    + actual.transmitted_fraction
                    + actual.absorbed_fraction
                    - 1.0)
                    .abs()
                    < 8e-14
            );
        }
    }
}

#[test]
fn finite_fluid_transmits_evanescent_sound_without_turning_it_into_heat() {
    let exterior = fluid();
    let faster = FluidAcousticState::try_new(2.4, 425.0).unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    let cosine = 0.5;
    // Independent closed form for a symmetric fluid slab. Its interior k_n
    // is imaginary, but the coupled two interfaces transmit finite power.
    let beta = (1.0_f64 - cosine * cosine - 0.8 * 0.8).sqrt();
    let a = (2.4 / 1.2) * cosine;
    let mut previous = 1.0;
    for thickness in [0.005, 0.05, 0.5] {
        let stack = TransmittingLayers::new(
            &exterior,
            &exterior,
            omega,
            &[Layer::FluidLayer {
                fluid: &faster,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let actual = stack.scattering(cosine).unwrap();
        let x = (omega / 340.0) * beta * thickness;
        let denominator = C64::new(2.0 * x.cosh(), (a / beta - beta / a) * x.sinh());
        near(
            actual.pressure_transmission,
            C64::from_re(2.0) / denominator,
            2e-14,
        );
        near(
            actual.pressure_reflection,
            C64::new(0.0, (a / beta + beta / a) * x.sinh()) / denominator,
            2e-14,
        );
        assert!(actual.transmitted_fraction > 0.0 && actual.transmitted_fraction < previous);
        assert_eq!(actual.absorbed_fraction, 0.0);
        previous = actual.transmitted_fraction;
    }
    // A thick evanescent layer tends to the single-interface phase, retaining
    // zero heat even after its transmitted amplitude underflows to zero.
    let thick = TransmittingLayers::new(
        &exterior,
        &exterior,
        omega,
        &[Layer::FluidLayer {
            fluid: &faster,
            thickness_m: 1000.0,
        }],
    )
    .unwrap();
    let actual = thick.scattering(cosine).unwrap();
    near(
        actual.pressure_reflection,
        C64::new(1.0, beta / a) / C64::new(1.0, -beta / a),
        2e-14,
    );
    assert_eq!(actual.pressure_transmission, C64::ZERO);
    assert_eq!(actual.transmitted_fraction, 0.0);
    assert_eq!(actual.absorbed_fraction, 0.0);
}

#[test]
fn critical_fluid_layer_has_the_finite_pressure_velocity_limit() {
    let exterior = fluid();
    let faster = FluidAcousticState::try_new(2.4, 425.0).unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    let thickness = 0.07;
    let stack = TransmittingLayers::new(
        &exterior,
        &exterior,
        omega,
        &[Layer::FluidLayer {
            fluid: &faster,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    // At sin(theta)=c0/c1=0.8, dp/dz=-i*omega*rho*v_n and dv_n/dz=0.
    // The finite transfer is [1, i*omega*rho*d/Z_ref; 0, 1], not an
    // infinite normal impedance or a zero-throughput half-space.
    let inertia = omega * 2.4 * thickness * 0.6 / (1.2 * 340.0);
    let denominator = C64::new(2.0, inertia);
    let expected_t = C64::from_re(2.0) / denominator;
    let expected_r = C64::new(0.0, inertia) / denominator;
    let exact = stack.scattering(0.6).unwrap();
    near(exact.pressure_transmission, expected_t, 3e-14);
    near(exact.pressure_reflection, expected_r, 3e-14);
    near(
        stack
            .prepared_layers()
            .reflection(0.6)
            .unwrap()
            .pressure_reflection,
        C64::ONE,
        1e-15,
    );
    for cosine in [0.6 - 1e-10, 0.6 + 1e-10] {
        let adjacent = stack.scattering(cosine).unwrap();
        near(adjacent.pressure_transmission, expected_t, 2e-9);
        near(adjacent.pressure_reflection, expected_r, 2e-9);
        assert_eq!(adjacent.absorbed_fraction, 0.0);
    }
}

#[test]
fn fluid_matching_layer_and_original_source_snapshot_are_retained() {
    struct Source {
        samples: Cell<usize>,
        state: FluidAcousticState,
    }
    impl FluidAcoustics for Source {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            self.samples.set(self.samples.get() + 1);
            Ok(self.state)
        }
    }
    let exterior = fluid();
    let receiver = FluidAcousticState::try_new(2.4, 425.0).unwrap();
    let matching = Source {
        samples: Cell::new(0),
        state: FluidAcousticState::try_new(
            ((1.2_f64 * 340.0) * (2.4 * 425.0)).sqrt() / 425.0,
            425.0,
        )
        .unwrap(),
    };
    let omega = core::f64::consts::TAU * 500.0;
    let eighth_wave = 425.0 / (8.0 * 500.0);
    let stack = TransmittingLayers::new(
        &exterior,
        &receiver,
        omega,
        &[
            Layer::FluidLayer {
                fluid: &matching,
                thickness_m: eighth_wave,
            },
            Layer::FluidLayer {
                fluid: &matching,
                thickness_m: eighth_wave,
            },
        ],
    )
    .unwrap();
    assert_eq!(matching.samples.get(), 1);
    let original: &dyn FluidAcoustics = &matching;
    assert!(
        matches!(stack.prepared_layers().layer(0), Some(Layer::FluidLayer { fluid, .. }) if core::ptr::eq(fluid, original))
    );
    for _ in 0..3 {
        let response = stack.scattering(1.0).unwrap();
        assert!(response.reflected_fraction < 1e-28);
        assert!((response.transmitted_fraction - 1.0).abs() < 2e-14);
        assert_eq!(response.absorbed_fraction, 0.0);
    }
    assert_eq!(matching.samples.get(), 1);

    struct RefusingSource;
    impl FluidAcoustics for RefusingSource {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Err(MaterialError::State {
                what: "original source refuses this state".into(),
            })
        }
    }
    // Even zero thickness does not bypass the original source's admission.
    assert!(
        RigidBackedLayers::new(
            &exterior,
            omega,
            &[Layer::FluidLayer {
                fluid: &RefusingSource,
                thickness_m: 0.0
            }]
        )
        .is_err()
    );
}

#[test]
fn transmitting_porous_gap_stack_matches_full_independent_matrices() {
    let fluid = fluid();
    let first = foam();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &first,
            thickness_m: 0.03,
        },
        Layer::FluidGap { thickness_m: 0.05 },
        Layer::Porous {
            law: &second,
            thickness_m: 0.02,
        },
    ];
    for (density, speed) in [(1.2, 340.0), (2.4, 425.0), (0.8, 280.0)] {
        let receiver = FluidAcousticState::try_new(density, speed).unwrap();
        for frequency in [150.0, 500.0, 1600.0] {
            let stack = TransmittingLayers::new(
                &fluid,
                &receiver,
                core::f64::consts::TAU * frequency,
                &layers,
            )
            .unwrap();
            for cosine in [1.0, 0.9, 0.5, 0.05] {
                let q2 = (340.0 / speed) * (340.0 / speed) - (1.0 - cosine * cosine);
                let q = if q2 >= 0.0 {
                    C64::from_re(q2.sqrt())
                } else {
                    C64::new(0.0, -(-q2).sqrt())
                };
                let load = q.scale((1.2 / density) / cosine);
                let (expected_r, expected_t) =
                    matrix_response(stack.prepared_layers(), cosine, load);
                let actual = stack.scattering(cosine).unwrap();
                near(actual.pressure_reflection, expected_r, 5e-14);
                near(actual.pressure_transmission, expected_t, 5e-14);
                assert!(
                    (actual.transmitted_fraction - load.re * expected_t.norm_sq()).abs() < 5e-14
                );
                assert!(
                    (actual.reflected_fraction
                        + actual.transmitted_fraction
                        + actual.absorbed_fraction
                        - 1.0)
                        .abs()
                        < 5e-14
                );
                assert!(actual.absorbed_fraction >= 0.0);
            }
        }
    }
}

#[test]
fn transmitting_reciprocity_and_subdivision_retain_complex_phase() {
    let fluid = fluid();
    let first = foam();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    let front = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &first,
                thickness_m: 0.03,
            },
            Layer::FluidGap { thickness_m: 0.05 },
            Layer::Porous {
                law: &second,
                thickness_m: 0.02,
            },
        ],
    )
    .unwrap();
    let back = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &second,
                thickness_m: 0.02,
            },
            Layer::FluidGap { thickness_m: 0.05 },
            Layer::Porous {
                law: &first,
                thickness_m: 0.03,
            },
        ],
    )
    .unwrap();
    let split = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &first,
                thickness_m: 0.01,
            },
            Layer::Porous {
                law: &first,
                thickness_m: 0.02,
            },
            Layer::FluidGap { thickness_m: 0.05 },
            Layer::Porous {
                law: &second,
                thickness_m: 0.02,
            },
        ],
    )
    .unwrap();
    for cosine in [1.0, 0.9, 0.5] {
        let a = front.scattering(cosine).unwrap();
        let b = back.scattering(cosine).unwrap();
        let c = split.scattering(cosine).unwrap();
        near(a.pressure_transmission, b.pressure_transmission, 5e-14);
        assert!((a.pressure_reflection - b.pressure_reflection).abs() > 1e-4);
        near(a.pressure_transmission, c.pressure_transmission, 5e-14);
        near(a.pressure_reflection, c.pressure_reflection, 5e-14);
    }
}

#[test]
fn sixty_four_matched_gaps_transmit_delay_without_reflection_or_loss() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        core::f64::consts::TAU * 500.0,
        &[Layer::FluidGap { thickness_m: 0.003 }; MAX_LAYERS],
    )
    .unwrap();
    for cosine in [1.0, 0.5, 1e-200, 1e-300] {
        let response = stack.scattering(cosine).unwrap();
        let phase = stack.prepared_layers().exterior_wavenumber_rad_m()
            * cosine
            * 0.003
            * (MAX_LAYERS as f64);
        near(response.pressure_reflection, C64::ZERO, 3e-14);
        near(
            response.pressure_transmission,
            C64::new(phase.cos(), -phase.sin()),
            3e-14,
        );
        assert_eq!(response.absorbed_fraction, 0.0);
        assert!((response.transmitted_fraction - 1.0).abs() < 5e-14);
    }
}

#[test]
fn transmitting_thick_absorber_has_no_fictitious_leak_or_receiver_dependence() {
    let fluid = fluid();
    let law = foam();
    let receiver = FluidAcousticState::try_new(2.4, 425.0).unwrap();
    let stack = TransmittingLayers::new(
        &fluid,
        &receiver,
        core::f64::consts::TAU * 500.0,
        &[Layer::Porous {
            law: &law,
            thickness_m: 1e308,
        }],
    )
    .unwrap();
    for cosine in [1.0, 0.5, 0.01, 1e-200] {
        let response = stack.scattering(cosine).unwrap();
        let rigid = stack.prepared_layers().reflection(cosine).unwrap();
        near(
            response.pressure_reflection,
            rigid.pressure_reflection,
            5e-14,
        );
        assert_eq!(response.pressure_transmission, C64::ZERO);
        assert_eq!(response.transmitted_fraction, 0.0);
        assert!(response.absorbed_fraction > 0.0);
        assert!((response.absorbed_fraction - rigid.absorbed_fraction).abs() < 5e-14);
    }
}

#[test]
fn tiny_pressure_transmission_keeps_weighted_power_and_small_normal_wavenumber() {
    let fluid = fluid();
    let receiver = FluidAcousticState::try_new(1e-200, 340.0).unwrap();
    let stack = TransmittingLayers::new(&fluid, &receiver, 1000.0, &[]).unwrap();
    let response = stack.scattering(1.0).unwrap();
    assert!(response.pressure_transmission.re > 0.0);
    assert_eq!(response.pressure_transmission.norm_sq(), 0.0);
    assert!(response.transmitted_fraction > 0.0);
    assert!((response.transmitted_fraction / (4.0e-200 / 1.2) - 1.0).abs() < 4e-13);
    // Matched specific impedance, with a speed ratio whose square underflows.
    let fast = FluidAcousticState::try_new(1.2e-200, 340e200).unwrap();
    let normal = TransmittingLayers::new(&fluid, &fast, 1000.0, &[])
        .unwrap()
        .scattering(1.0)
        .unwrap();
    assert!(normal.receiving_normal_wavenumber_rad_m.re > 0.0);
    near(normal.pressure_reflection, C64::ZERO, 3e-14);
    near(normal.pressure_transmission, C64::ONE, 3e-14);
    assert!((normal.transmitted_fraction - 1.0).abs() < 3e-14);
    // Equal speeds make normal admittance contrast independent of angle, even
    // if multiplying density contrast by cosine first creates a subnormal.
    let dense = FluidAcousticState::try_new(2.4e123, 340.0).unwrap();
    let matched_speed = TransmittingLayers::new(&fluid, &dense, 1000.0, &[]).unwrap();
    let normal_power = matched_speed.scattering(1.0).unwrap().transmitted_fraction;
    let grazing_power = matched_speed
        .scattering(1e-200)
        .unwrap()
        .transmitted_fraction;
    assert!(grazing_power > 0.0);
    assert!((grazing_power / normal_power - 1.0).abs() < 3e-13);
}

struct CountedReceiver {
    state: FluidAcousticState,
    samples: Cell<usize>,
}
impl FluidAcoustics for CountedReceiver {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn one_original_source_has_one_state_across_incident_layers_and_receiver() {
    // A second query deliberately advances to a different coherent state.
    // Composing repeated calls would create a fictitious refracting interface.
    struct EvolvingFluid {
        states: [ThermoviscousFluidState; 2],
        samples: Cell<usize>,
    }
    impl EvolvingFluid {
        fn next(&self) -> ThermoviscousFluidState {
            let index = self.samples.get();
            self.samples.set(index + 1);
            self.states[index.min(1)]
        }
    }
    impl FluidAcoustics for EvolvingFluid {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Ok(self.next().acoustic())
        }
    }
    impl FluidTransport for EvolvingFluid {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            Ok(self.next().transport())
        }
    }
    impl ThermoviscousFluid for EvolvingFluid {
        fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
            Ok(self.next())
        }
    }
    let first = fluid();
    let second = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(2.4, 425.0).unwrap(),
        FluidTransportState::try_new(2.4, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    for with_layers in [false, true] {
        let source = EvolvingFluid {
            states: [first, second],
            samples: Cell::new(0),
        };
        let layers = [Layer::FluidLayer {
            fluid: &source,
            thickness_m: 0.01,
        }; 2];
        let stack = TransmittingLayers::new(
            &source,
            &source,
            omega,
            if with_layers { &layers } else { &[] },
        )
        .unwrap();
        assert_eq!(source.samples.get(), 1);
        assert_eq!(stack.receiving_state(), &first.acoustic());
        for cosine in [1.0, 0.6, 0.01] {
            let response = stack.scattering(cosine).unwrap();
            assert!(response.reflected_fraction < 1e-28);
            assert!((response.transmitted_fraction - 1.0).abs() < 2e-14);
        }
        assert_eq!(source.samples.get(), 1);
    }
    // The incident source is distinct; a finite layer and the receiver still
    // share their first evaluation. Equal-valued but distinct owners do not.
    let interior = EvolvingFluid {
        states: [second, first],
        samples: Cell::new(0),
    };
    let distinct = EvolvingFluid {
        states: [second, first],
        samples: Cell::new(0),
    };
    let layers = [Layer::FluidLayer {
        fluid: &interior,
        thickness_m: 0.01,
    }];
    let shared = TransmittingLayers::new(&first, &interior, omega, &layers).unwrap();
    assert_eq!(interior.samples.get(), 1);
    assert_eq!(shared.receiving_state(), &second.acoustic());
    let separate = TransmittingLayers::new(&first, &distinct, omega, &layers).unwrap();
    assert_eq!(interior.samples.get(), 2);
    assert_eq!(distinct.samples.get(), 1);
    assert_eq!(
        separate.prepared_layers().fluid_response(0),
        Some(&first.acoustic())
    );
    assert_eq!(separate.receiving_state(), &second.acoustic());
    shared.scattering(0.8).unwrap();
    separate.scattering(0.8).unwrap();
    assert_eq!(interior.samples.get(), 2);
    assert_eq!(distinct.samples.get(), 1);
}

#[test]
fn transmission_prepares_both_original_sources_once_and_preserves_admission() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let receiving = CountedReceiver {
        state: FluidAcousticState::try_new(2.4, 425.0).unwrap(),
        samples: Cell::new(0),
    };
    let law = foam();
    let layers = [Layer::Porous {
        law: &law,
        thickness_m: 0.03,
    }];
    let stack = TransmittingLayers::new(
        &incident,
        &receiving,
        core::f64::consts::TAU * 500.0,
        &layers,
    )
    .unwrap();
    for cosine in [1.0, 0.6, 0.5] {
        stack.scattering(cosine).unwrap();
    }
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiving.samples.get(), 1);
    assert_eq!(stack.receiving_state(), &receiving.state);
    assert_eq!(
        stack
            .prepared_layers()
            .material_response(0)
            .unwrap()
            .fluid(),
        &incident.state
    );
    for cosine in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(stack.scattering(cosine).is_err());
    }
    assert!(
        TransmittingLayers::new(
            &incident,
            &receiving,
            core::f64::consts::TAU * 149.0,
            &layers
        )
        .is_err()
    );
    assert!(
        TransmittingLayers::new(
            &incident,
            &receiving,
            1000.0,
            &[Layer::FluidGap { thickness_m: 0.01 }; MAX_LAYERS + 1]
        )
        .is_err()
    );
}

#[test]
fn angular_power_uses_incident_watts_without_averaging_complex_pressure_or_cosines_twice() {
    let incident = fluid();
    let law = foam();
    let layers = [Layer::Porous {
        law: &law,
        thickness_m: 0.04,
    }];
    let rigid = RigidBackedLayers::new(&incident, core::f64::consts::TAU * 500.0, &layers).unwrap();
    let transmitting = TransmittingLayers::new(
        &incident,
        &incident,
        core::f64::consts::TAU * 500.0,
        &layers,
    )
    .unwrap();
    let samples = [
        IncidentPower {
            cosine: 1.0,
            power_w: 2.0,
        },
        IncidentPower {
            cosine: 0.25,
            power_w: 6.0,
        },
        IncidentPower {
            cosine: 0.5,
            power_w: 0.0,
        },
    ];
    let budget = AngularPowerBudget {
        max_samples: samples.len(),
    };
    let first = rigid.reflection(1.0).unwrap();
    let second = rigid.reflection(0.25).unwrap();
    let result = rigid.angular_power(&samples, budget, || false).unwrap();
    assert_eq!(result.incident_power_w, 8.0);
    assert_eq!(result.evaluated_directions, 2);
    assert_eq!(result.transmitted_power_w, 0.0);
    assert!(
        (result.absorbed_power_w
            - (2.0 * first.absorbed_fraction + 6.0 * second.absorbed_fraction))
            .abs()
            < 1e-14
    );
    let coherent =
        (first.pressure_reflection.scale(0.25) + second.pressure_reflection.scale(0.75)).norm_sq();
    assert!((result.reflected_fraction - coherent).abs() > 1e-4);
    assert!(result.balance_residual_w.abs() < 1e-13);

    let result = transmitting
        .angular_power(&samples, budget, || false)
        .unwrap();
    let first = transmitting.scattering(1.0).unwrap();
    let second = transmitting.scattering(0.25).unwrap();
    for (actual, expected) in [
        (
            result.reflected_power_w,
            2.0 * first.reflected_fraction + 6.0 * second.reflected_fraction,
        ),
        (
            result.transmitted_power_w,
            2.0 * first.transmitted_fraction + 6.0 * second.transmitted_fraction,
        ),
        (
            result.absorbed_power_w,
            2.0 * first.absorbed_fraction + 6.0 * second.absorbed_fraction,
        ),
    ] {
        assert!((actual - expected).abs() < 2e-14);
    }
    assert!(result.absorbed_power_w > 0.0 && result.transmitted_power_w > 0.0);
    assert!(result.balance_residual_w.abs() < 1e-13);
}

#[test]
fn diffuse_porous_absorption_refines_to_an_independent_polar_angle_integral() {
    let incident = fluid();
    let law = foam();
    let layers = [
        Layer::Porous {
            law: &law,
            thickness_m: 0.04,
        },
        Layer::FluidGap { thickness_m: 0.02 },
    ];
    let stack = RigidBackedLayers::new(&incident, core::f64::consts::TAU * 500.0, &layers).unwrap();
    // Independent growing-transfer-matrix oracle in THETA rather than the
    // production mu quadrature; coefficients include internal refraction.
    for minimum_cosine in [0.0_f64, 0.4] {
        let maximum_angle = minimum_cosine.acos();
        let intervals = 16_384;
        let width = maximum_angle / intervals as f64;
        let mut reference = 0.0;
        for i in 0..intervals {
            let theta = (i as f64 + 0.5) * width;
            let reflection = matrix_reflection(&stack, theta.cos());
            reference += (1.0 - reflection.norm_sq()) * theta.cos() * theta.sin();
        }
        reference *= 2.0 * width / (1.0 - minimum_cosine * minimum_cosine);
        let mut errors = Vec::new();
        for samples in [16, 32, 64] {
            let result = stack
                .diffuse_power(
                    7.0,
                    DiffuseIncidence {
                        minimum_cosine,
                        polar_samples: samples,
                    },
                    AngularPowerBudget {
                        max_samples: samples,
                    },
                    || false,
                )
                .unwrap();
            assert_eq!(result.incident_power_w, 7.0);
            assert_eq!(result.evaluated_directions, samples);
            assert!(result.balance_residual_w.abs() < 5e-14);
            errors.push((result.absorbed_fraction - reference).abs());
        }
        assert!(
            errors[1] < 0.3 * errors[0],
            "mu_min={minimum_cosine}: {errors:?}"
        );
        assert!(
            errors[2] < 0.3 * errors[1],
            "mu_min={minimum_cosine}: {errors:?}"
        );
        assert!(errors[2] < 5e-5, "mu_min={minimum_cosine}: {errors:?}");
        if minimum_cosine == 0.0 {
            assert!((reference - stack.reflection(1.0).unwrap().absorbed_fraction).abs() > 0.02);
        }
    }
}

#[test]
fn diffuse_limits_keep_matched_transmission_critical_refraction_and_zero_absorption_distinct() {
    let incident = fluid();
    let budget = AngularPowerBudget { max_samples: 4096 };
    let diffuse = DiffuseIncidence {
        minimum_cosine: 0.0,
        polar_samples: 4096,
    };
    let rigid = RigidBackedLayers::new(&incident, 1000.0, &[]).unwrap();
    let reflection = rigid.diffuse_power(3.0, diffuse, budget, || false).unwrap();
    assert_eq!(reflection.reflected_power_w, 3.0);
    assert_eq!(reflection.absorbed_power_w, 0.0);
    assert_eq!(reflection.transmitted_power_w, 0.0);
    let matched = TransmittingLayers::new(&incident, &incident, 1000.0, &[]).unwrap();
    let result = matched
        .diffuse_power(3.0, diffuse, budget, || false)
        .unwrap();
    assert!(result.reflected_power_w < 1e-25);
    assert!((result.transmitted_power_w - 3.0).abs() < 3e-14);
    assert_eq!(result.absorbed_power_w, 0.0);

    let receiver = FluidAcousticState::try_new(1.2, 680.0).unwrap();
    let faster = TransmittingLayers::new(&incident, &receiver, 1000.0, &[]).unwrap();
    let evanescent = faster.scattering(0.5).unwrap();
    assert_eq!(evanescent.transmitted_fraction, 0.0);
    assert!(evanescent.pressure_transmission.abs() > 0.0);
    let result = faster
        .diffuse_power(3.0, diffuse, budget, || false)
        .unwrap();
    // Only sin(theta)<c_inc/c_recv=1/2 carries transmitted power. Its diffuse
    // incident-power share is 1/4, an independent exact upper bound.
    assert!(result.transmitted_fraction > 0.0 && result.transmitted_fraction < 0.25);
    assert_eq!(result.absorbed_power_w, 0.0);
    assert!(result.balance_residual_w.abs() < 5e-14);
    let zero = faster
        .diffuse_power(0.0, diffuse, budget, || false)
        .unwrap();
    assert_eq!(zero.incident_power_w, 0.0);
    assert_eq!(zero.reflected_power_w, 0.0);
    assert_eq!(zero.transmitted_power_w, 0.0);
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert_eq!(zero.transmitted_fraction, result.transmitted_fraction);
}

#[test]
fn angular_queries_retain_sources_and_refuse_invalid_weights_budgets_and_cancellation() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let receiver = CountedReceiver {
        state: FluidAcousticState::try_new(2.4, 340.0).unwrap(),
        samples: Cell::new(0),
    };
    let law = foam();
    let stack = TransmittingLayers::new(
        &incident,
        &receiver,
        core::f64::consts::TAU * 500.0,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.03,
        }],
    )
    .unwrap();
    let budget = AngularPowerBudget { max_samples: 32 };
    let diffuse = DiffuseIncidence {
        minimum_cosine: 0.0,
        polar_samples: 32,
    };
    let baseline = stack.diffuse_power(2.0, diffuse, budget, || false).unwrap();
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
    assert_eq!(
        stack.diffuse_power(2.0, diffuse, budget, || true),
        Err(AngularPowerError::Cancelled),
    );
    let mut checkpoints = 0;
    assert_eq!(
        stack.diffuse_power(2.0, diffuse, budget, || {
            checkpoints += 1;
            checkpoints == 10
        }),
        Err(AngularPowerError::Cancelled),
    );
    assert_eq!(
        stack.diffuse_power(2.0, diffuse, budget, || false).unwrap(),
        baseline
    );
    for maximum in [0, 31, MAX_ANGULAR_POWER_SAMPLES + 1] {
        assert!(
            stack
                .diffuse_power(
                    2.0,
                    diffuse,
                    AngularPowerBudget {
                        max_samples: maximum
                    },
                    || false
                )
                .is_err()
        );
    }
    for minimum in [-0.1, 1.0, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .diffuse_power(
                    2.0,
                    DiffuseIncidence {
                        minimum_cosine: minimum,
                        ..diffuse
                    },
                    budget,
                    || false
                )
                .is_err()
        );
    }
    for power in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .angular_power(
                    &[IncidentPower {
                        cosine: 1.0,
                        power_w: power
                    }],
                    budget,
                    || false
                )
                .is_err()
        );
    }
    for cosine in [0.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .angular_power(
                    &[IncidentPower {
                        cosine,
                        power_w: 0.0
                    }],
                    budget,
                    || false
                )
                .is_err()
        );
    }
    assert!(stack.angular_power(&[], budget, || false).is_err());
    assert!(
        stack
            .angular_power(
                &[IncidentPower {
                    cosine: 1.0,
                    power_w: 0.0
                }],
                budget,
                || false
            )
            .is_err()
    );
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
}

// Conventional, unscaled transfer of pressure and Z_ref*v_n. This oracle is
// deliberately limited to well-conditioned layers; the production report uses
// neither these growing functions nor spatial quadrature.
fn layer_state_oracle(kn: C64, inertia: C64, d: f64, right: [C64; 2]) -> [C64; 2] {
    if kn == C64::ZERO {
        return [right[0] + C64::new(0.0, d) * inertia * right[1], right[1]];
    }
    let (cos, sin) = cos_sin(kn.scale(d));
    let i_sin = C64::new(0.0, 1.0) * sin;
    let ratio = inertia / kn;
    [
        cos * right[0] + ratio * i_sin * right[1],
        i_sin / ratio * right[0] + cos * right[1],
    ]
}

fn layer_loss_oracle(
    stack: &RigidBackedLayers<'_>,
    cosine: f64,
    load: C64,
) -> (Vec<f64>, Vec<[f64; 2]>) {
    let omega = stack.angular_frequency_rad_s();
    let zref = stack.fluid().acoustic().characteristic_impedance_pa_s_m() / cosine;
    let kt = stack.exterior_wavenumber_rad_m() * (1.0 - cosine * cosine).sqrt();
    let parameters: Vec<_> = (0..stack.layer_count())
        .map(|i| {
            let (k, density, d) = match stack.layer(i).unwrap() {
                Layer::Porous { thickness_m, .. } => {
                    let response = stack.material_response(i).unwrap();
                    (
                        response.characteristic().1,
                        response.dynamic_density_kg_m3(),
                        thickness_m,
                    )
                }
                Layer::FluidGap { thickness_m } => (
                    C64::from_re(stack.exterior_wavenumber_rad_m()),
                    C64::from_re(stack.fluid().acoustic().density_kg_m3()),
                    thickness_m,
                ),
                Layer::FluidLayer { thickness_m, .. } => {
                    let fluid = stack.fluid_response(i).unwrap();
                    (
                        C64::from_re(omega / fluid.sound_speed_m_s()),
                        C64::from_re(fluid.density_kg_m3()),
                        thickness_m,
                    )
                }
            };
            let squared = k * k - C64::from_re(kt * kt);
            let mut kn = squared.sqrt();
            if kn.im > 0.0 {
                kn = -kn;
            }
            (kn, density.scale(omega / zref), d)
        })
        .collect();
    let mut states = vec![[C64::ZERO; 2]; stack.layer_count() + 1];
    states[stack.layer_count()] = [C64::ONE, load];
    for i in (0..stack.layer_count()).rev() {
        let (kn, inertia, d) = parameters[i];
        states[i] = layer_state_oracle(kn, inertia, d, states[i + 1]);
    }
    let amplitude = C64::from_re(2.0) / (states[0][0] + states[0][1]);
    let fluxes: Vec<_> = states
        .iter()
        .map(|s| {
            let p = amplitude * s[0];
            let q = amplitude * s[1];
            (p * q.conj()).re
        })
        .collect();
    let losses = fluxes.windows(2).map(|f| f[0] - f[1]).collect();
    let integrals = (0..stack.layer_count())
        .map(|i| {
            let Some(response) = stack.material_response(i) else {
                return [0.0; 2];
            };
            let density = response.dynamic_density_kg_m3();
            let compliance = response.dynamic_bulk_modulus_pa().recip();
            let (kn, inertia, d) = parameters[i];
            let panels = 1024;
            let mut values = [0.0; 2];
            for j in 0..=panels {
                let remaining = d * (1.0 - j as f64 / panels as f64);
                let state = layer_state_oracle(kn, inertia, remaining, states[i + 1]);
                let p = amplitude * state[0];
                let normal_velocity = (amplitude * state[1]).scale(1.0 / zref);
                let tangential_velocity = if kt == 0.0 {
                    C64::ZERO
                } else {
                    p / density.scale(omega / kt)
                };
                let viscous = omega
                    * -density.im
                    * zref
                    * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
                let thermal = omega * -compliance.im * zref * p.norm_sq();
                let weight = if j == 0 || j == panels {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                values[0] += weight * viscous;
                values[1] += weight * thermal;
            }
            values.map(|v| v * d / (3.0 * panels as f64))
        })
        .collect();
    (losses, integrals)
}

#[test]
fn g1_layer_absorption_matches_interface_work_and_independent_spatial_integrals() {
    let fluid = fluid();
    let first = foam();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let middle = FluidAcousticState::try_new(0.9, 425.0).unwrap();
    let layers = [
        Layer::Porous {
            law: &first,
            thickness_m: 0.015,
        },
        Layer::FluidGap { thickness_m: 0.01 },
        Layer::FluidLayer {
            fluid: &middle,
            thickness_m: 0.02,
        },
        Layer::Porous {
            law: &second,
            thickness_m: 0.025,
        },
    ];
    let budget = LayerAbsorptionBudget {
        max_layers: layers.len(),
    };
    for frequency in [400.0, 1100.0] {
        let stack =
            TransmittingLayers::new(&fluid, &fluid, core::f64::consts::TAU * frequency, &layers)
                .unwrap();
        for cosine in [1.0, 0.7, 0.2] {
            for rigid in [false, true] {
                let report = if rigid {
                    stack
                        .prepared_layers()
                        .layer_absorption(cosine, 7.0, budget, || false)
                        .unwrap()
                } else {
                    stack
                        .layer_absorption(cosine, 7.0, budget, || false)
                        .unwrap()
                };
                let (work, integrated) = layer_loss_oracle(
                    stack.prepared_layers(),
                    cosine,
                    if rigid { C64::ZERO } else { C64::ONE },
                );
                for (i, row) in report.layers().iter().enumerate() {
                    assert!((row.absorbed_fraction() - work[i]).abs() < 8e-13);
                    assert!((row.viscous_absorbed_fraction - integrated[i][0]).abs() < 2e-11);
                    assert!((row.thermal_absorbed_fraction - integrated[i][1]).abs() < 2e-11);
                    assert_eq!(row.absorbed_power_w, 7.0 * row.absorbed_fraction());
                    assert!(row.viscous_absorbed_fraction >= 0.0);
                    assert!(row.thermal_absorbed_fraction >= 0.0);
                }
                assert_eq!(report.layers()[1].absorbed_fraction(), 0.0);
                assert_eq!(report.layers()[2].absorbed_fraction(), 0.0);
                assert!(report.balance_residual_w.abs() < 7e-12);
                if rigid {
                    let response = stack.prepared_layers().reflection(cosine).unwrap();
                    assert_eq!(
                        report.reflected_fraction,
                        response.pressure_reflection.norm_sq()
                    );
                    assert_eq!(report.transmitted_fraction, 0.0);
                    assert!((report.absorbed_fraction - response.absorbed_fraction).abs() < 8e-13);
                } else {
                    let response = stack.scattering(cosine).unwrap();
                    assert_eq!(report.reflected_fraction, response.reflected_fraction);
                    assert_eq!(report.transmitted_fraction, response.transmitted_fraction);
                    assert!((report.absorbed_fraction - response.absorbed_fraction).abs() < 8e-13);
                }
            }
        }
    }
}

#[test]
fn g3_layer_absorption_subdivision_preserves_mechanisms_and_declared_zero_layers() {
    let fluid = fluid();
    let law = foam();
    let omega = core::f64::consts::TAU * 800.0;
    let whole = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.06,
        }],
    )
    .unwrap();
    let split = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.01,
            },
            Layer::FluidGap { thickness_m: 0.0 },
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
            Layer::Porous {
                law: &law,
                thickness_m: 0.0,
            },
            Layer::Porous {
                law: &law,
                thickness_m: 0.03,
            },
        ],
    )
    .unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 5 };
    for cosine in [1.0, 0.5, 0.01] {
        let a = whole
            .layer_absorption(cosine, 3.0, budget, || false)
            .unwrap();
        let b = split
            .layer_absorption(cosine, 3.0, budget, || false)
            .unwrap();
        let viscous: f64 = b.layers().iter().map(|l| l.viscous_absorbed_fraction).sum();
        let thermal: f64 = b.layers().iter().map(|l| l.thermal_absorbed_fraction).sum();
        assert!((a.layers()[0].viscous_absorbed_fraction - viscous).abs() < 2e-13);
        assert!((a.layers()[0].thermal_absorbed_fraction - thermal).abs() < 2e-13);
        assert_eq!(b.layers().len(), 5);
        assert_eq!(b.layers()[1].absorbed_power_w, 0.0);
        assert_eq!(b.layers()[3].absorbed_power_w, 0.0);
        assert!((a.absorbed_power_w - b.absorbed_power_w).abs() < 1e-12);
    }
}

#[test]
fn g1_local_loss_retains_sub_roundoff_transmission_and_thin_rigid_velocity_nodes() {
    let base = fluid();
    let fluid =
        ThermoviscousFluidState::from_heat_capacity_ratio(base.acoustic(), base.transport(), 1.0)
            .unwrap();
    let law = foam();
    let omega = core::f64::consts::TAU * 1000.0;
    let response = law.effective_fluid(&fluid, omega).unwrap();
    let density = response.dynamic_density_kg_m3();
    let bulk = response.dynamic_bulk_modulus_pa();
    assert_eq!(bulk.im, 0.0);
    let z0 = fluid.acoustic().characteristic_impedance_pa_s_m();
    let budget = LayerAbsorptionBudget { max_layers: 1 };
    let d = 1e-20;
    let transmitted = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[Layer::Porous {
            law: &law,
            thickness_m: d,
        }],
    )
    .unwrap();
    for cosine in [1.0, 0.4] {
        let report = transmitted
            .layer_absorption(cosine, 2.0, budget, || false)
            .unwrap();
        let zref = z0 / cosine;
        let kt = (omega / fluid.acoustic().sound_speed_m_s()) * (1.0 - cosine * cosine).sqrt();
        // Direct zero-thickness field limit: p=1, v_n=1/Z_ref and v_t=k_t/(omega*rho).
        let expected = omega
            * -density.im
            * zref
            * d
            * (1.0 / (zref * zref) + kt * kt / (omega * omega * density.norm_sq()));
        assert!(expected > 0.0 && expected < f64::EPSILON);
        assert!((report.absorbed_fraction / expected - 1.0).abs() < 2e-12);
        assert_eq!(report.layers()[0].thermal_absorbed_fraction, 0.0);
    }
    for d in [1e-8, 1e-12, 1e-20] {
        let rigid = RigidBackedLayers::new(
            &fluid,
            omega,
            &[Layer::Porous {
                law: &law,
                thickness_m: d,
            }],
        )
        .unwrap();
        let report = rigid.layer_absorption(1.0, 2.0, budget, || false).unwrap();
        // Rigid back: p≈2 and v_n≈2 i omega (d-z)/K, giving a true d³ loss.
        let expected = 4.0 * omega.powi(3) * -density.im * z0 * d.powi(3) / (3.0 * bulk.norm_sq());
        assert!(report.absorbed_fraction > 0.0);
        assert!((report.absorbed_fraction / expected - 1.0).abs() < 2e-11);
        assert_eq!(report.layers()[0].thermal_absorbed_fraction, 0.0);
    }
}

#[test]
fn g0_layer_absorption_keeps_evanescent_lossless_and_opaque_limits_distinct() {
    let fluid = fluid();
    let receiver = FluidAcousticState::try_new(1.1, 425.0).unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    let budget = LayerAbsorptionBudget { max_layers: 2 };
    let lossless = TransmittingLayers::new(
        &fluid,
        &fluid,
        omega,
        &[
            Layer::FluidLayer {
                fluid: &receiver,
                thickness_m: 0.07,
            },
            Layer::FluidGap { thickness_m: 0.02 },
        ],
    )
    .unwrap();
    for cosine in [1.0, 0.8, 0.6, 0.2] {
        let report = lossless
            .layer_absorption(cosine, 4.0, budget, || false)
            .unwrap();
        assert_eq!(report.absorbed_power_w, 0.0);
        assert!(report.layers().iter().all(|p| p.absorbed_fraction() == 0.0));
        assert!(report.transmitted_power_w > 0.0);
    }
    let law = foam();
    let opaque = TransmittingLayers::new(
        &fluid,
        &receiver,
        omega,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 1e5,
            },
            Layer::Porous {
                law: &law,
                thickness_m: 0.02,
            },
        ],
    )
    .unwrap();
    let report = opaque.layer_absorption(0.2, 4.0, budget, || false).unwrap();
    assert_eq!(report.transmitted_power_w, 0.0);
    assert!(report.layers()[0].absorbed_power_w > 0.0);
    assert_eq!(report.layers()[1].absorbed_power_w, 0.0);
    assert!(
        (report.absorbed_fraction - opaque.scattering(0.2).unwrap().absorbed_fraction).abs()
            < 2e-13
    );
    let zero = opaque.layer_absorption(0.2, 0.0, budget, || false).unwrap();
    assert_eq!(zero.absorbed_fraction, report.absorbed_fraction);
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert!(zero.layers().iter().all(|l| l.absorbed_power_w == 0.0));
    let grazing = opaque
        .layer_absorption(1e-80, 4.0, budget, || false)
        .unwrap();
    let extreme = opaque
        .layer_absorption(1e-200, 4.0, budget, || false)
        .unwrap();
    for (a, b) in [
        (
            grazing.layers()[0].viscous_absorbed_fraction,
            extreme.layers()[0].viscous_absorbed_fraction,
        ),
        (
            grazing.layers()[0].thermal_absorbed_fraction,
            extreme.layers()[0].thermal_absorbed_fraction,
        ),
    ] {
        assert!(b > 0.0);
        assert!(((a / 1e-80) / (b / 1e-200) - 1.0).abs() < 3e-11);
    }
    let empty = RigidBackedLayers::new(&fluid, omega, &[]).unwrap();
    let empty = empty
        .layer_absorption(1.0, 3.0, LayerAbsorptionBudget { max_layers: 0 }, || false)
        .unwrap();
    assert!(empty.layers().is_empty());
    assert_eq!(empty.reflected_power_w, 3.0);
    assert_eq!(empty.absorbed_power_w, 0.0);
}

#[test]
fn g4_layer_absorption_retains_sources_and_refuses_work_and_cancellation_atomically() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let receiver = CountedReceiver {
        state: incident.state.acoustic(),
        samples: Cell::new(0),
    };
    let law = foam();
    let stack = TransmittingLayers::new(
        &incident,
        &receiver,
        core::f64::consts::TAU * 500.0,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: 0.01,
            },
            Layer::FluidGap { thickness_m: 0.02 },
            Layer::Porous {
                law: &law,
                thickness_m: 0.03,
            },
        ],
    )
    .unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 3 };
    let baseline = stack.layer_absorption(0.7, 5.0, budget, || false).unwrap();
    for stop in [1, 3, 6, 8] {
        let mut polls = 0;
        assert_eq!(
            stack.layer_absorption(0.7, 5.0, budget, || {
                polls += 1;
                polls == stop
            }),
            Err(LayerAbsorptionError::Cancelled)
        );
        assert_eq!(
            stack.layer_absorption(0.7, 5.0, budget, || false).unwrap(),
            baseline
        );
    }
    for max_layers in [0, 2, MAX_LAYERS + 1] {
        assert!(
            stack
                .layer_absorption(0.7, 5.0, LayerAbsorptionBudget { max_layers }, || false)
                .is_err()
        );
    }
    for cosine in [0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .layer_absorption(cosine, 5.0, budget, || false)
                .is_err()
        );
    }
    for power in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .layer_absorption(0.7, power, budget, || false)
                .is_err()
        );
    }
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
}
