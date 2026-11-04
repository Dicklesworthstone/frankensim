//! Exact finite-plate controls from original fluid and elastic material states.
use fs_blake3::ContentHash;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::layers::{
    ElasticPlateLayers, ElasticTransmittingLayers, Layer, LayerAbsorptionBudget,
    LayerAbsorptionError, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;

fn fluid(density: f64, speed: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(density, speed).unwrap(),
        FluidTransportState::try_new(density, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}

fn solid(density: f64, cp: f64, cs: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        density,
        density * (cp * cp - (4.0 / 3.0) * cs * cs),
        density * cs * cs,
        ContentHash([24; 32]),
        0.001,
    )
    .unwrap()
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

fn near_real(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

// Independent six boundary equations. Unknowns are R, forward P/S velocity at
// the front, backward P/S velocity at the back, and T. No parity reduction,
// equivalent impedance, production bounded-wave function or Wronskian is used.
fn boundary_oracle(
    fluid: [f64; 2],
    solid: [f64; 3],
    receiver: [f64; 2],
    omega: f64,
    thickness: f64,
    cosine: f64,
) -> [C64; 6] {
    let [rf, cf] = fluid;
    let [rs, cp, cs] = solid;
    let [rr, cr] = receiver;
    let sine = ((1.0 - cosine) * (1.0 + cosine)).sqrt();
    let outgoing = |squared: f64| {
        if squared >= 0.0 {
            C64::from_re(squared.sqrt())
        } else {
            C64::new(0.0, -(-squared).sqrt())
        }
    };
    let p = outgoing(1.0 - (cp * sine / cf).powi(2));
    let s = outgoing(1.0 - (cs * sine / cf).powi(2));
    let ss = cs * sine / cf;
    let a = 1.0 - 2.0 * ss * ss;
    let zp = rs * cp / (rf * cf);
    let zs = rs * cs / (rf * cf);
    let exp = |z: C64| C64::new(z.im.cos(), z.im.sin()).scale(z.re.exp());
    let ep = exp(C64::new(0.0, -omega * thickness / cp) * p);
    let es = exp(C64::new(0.0, -omega * thickness / cs) * s);
    let yr = outgoing((cf / cr).powi(2) - sine * sine).scale(rf / rr);
    let real = C64::from_re;
    let z = C64::ZERO;
    let mut matrix = [
        [
            real(1.0),
            real(-zp * a),
            s.scale(2.0 * zs * ss),
            ep.scale(-zp * a),
            (s * es).scale(2.0 * zs * ss),
            z,
            real(-1.0),
        ],
        [
            real(-cosine),
            -p,
            real(ss),
            p * ep,
            es.scale(-ss),
            z,
            real(-cosine),
        ],
        [
            z,
            p.scale(2.0 * ss),
            real(a),
            (p * ep).scale(-2.0 * ss),
            es.scale(-a),
            z,
            z,
        ],
        [
            z,
            ep.scale(zp * a),
            (s * es).scale(-2.0 * zs * ss),
            real(zp * a),
            s.scale(-2.0 * zs * ss),
            real(-1.0),
            z,
        ],
        [z, p * ep, es.scale(-ss), -p, real(ss), -yr, z],
        [
            z,
            (p * ep).scale(2.0 * ss),
            es.scale(a),
            p.scale(-2.0 * ss),
            real(-a),
            z,
            z,
        ],
    ];
    for column in 0..6 {
        let pivot = (column..6)
            .max_by(|i, j| {
                matrix[*i][column]
                    .abs()
                    .total_cmp(&matrix[*j][column].abs())
            })
            .unwrap();
        matrix.swap(column, pivot);
        let divisor = matrix[column][column];
        assert!(divisor.abs() > 0.0);
        for entry in column..7 {
            matrix[column][entry] = matrix[column][entry] / divisor;
        }
        for row in 0..6 {
            if row == column {
                continue;
            }
            let multiplier = matrix[row][column];
            for entry in column..7 {
                matrix[row][entry] = matrix[row][entry] - multiplier * matrix[column][entry];
            }
        }
    }
    core::array::from_fn(|index| matrix[index][6])
}

#[test]
fn normal_plate_is_the_exact_original_longitudinal_fluid_layer() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = fluid(800.0, 1200.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let omega = 12000.0;
    let cp = plate.directional_waves([0.0, 0.0, 1.0]).unwrap().modes[2].phase_speed_m_s;
    let p_limit = FluidAcousticState::try_new(2700.0, cp).unwrap();
    // Includes compression half-wave transparency and inactive shear resonances.
    for thickness in [
        1e-12,
        0.001,
        0.1,
        core::f64::consts::PI * cp / omega,
        core::f64::consts::TAU * cp / omega,
    ] {
        let stack =
            ElasticPlateLayers::new(&incident, &plate, &receiver, omega, thickness, &[]).unwrap();
        assert!(core::ptr::eq(stack.plate_state(), &plate));
        assert_eq!(stack.thickness_m(), thickness);
        let actual = stack.scattering(1.0).unwrap();
        let reference = TransmittingLayers::new(
            &incident,
            &receiver,
            omega,
            &[Layer::FluidLayer {
                fluid: &p_limit,
                thickness_m: thickness,
            }],
        )
        .unwrap()
        .scattering(1.0)
        .unwrap();
        near(
            actual.pressure_reflection,
            reference.pressure_reflection,
            4e-15,
        );
        near(
            actual.pressure_transmission,
            reference.pressure_transmission,
            4e-15,
        );
        near_real(
            actual.transmitted_fraction,
            reference.transmitted_fraction,
            4e-15,
        );
        assert_eq!(actual.absorbed_fraction, 0.0);
    }
    let transparent = ElasticPlateLayers::new(
        &incident,
        &plate,
        &incident,
        omega,
        core::f64::consts::PI * cp / omega,
        &[],
    )
    .unwrap()
    .scattering(1.0)
    .unwrap();
    near(transparent.pressure_transmission, -C64::ONE, 1e-14);
    near_real(transparent.transmitted_fraction, 1.0, 5e-15);
}

#[test]
fn zero_plate_recovers_the_existing_fluid_receiver_through_entrance_layers() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = fluid(750.0, 1800.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let layers = [Layer::FluidGap { thickness_m: 0.017 }];
    let actual =
        ElasticPlateLayers::new(&incident, &plate, &receiver, 9000.0, 0.0, &layers).unwrap();
    let expected = TransmittingLayers::new(&incident, &receiver, 9000.0, &layers).unwrap();
    for cosine in [1.0, 0.97, 0.8, 0.6, 0.3, 0.001] {
        let a = actual.scattering(cosine).unwrap();
        let b = expected.scattering(cosine).unwrap();
        assert_eq!(a.pressure_reflection, b.pressure_reflection);
        assert_eq!(a.pressure_transmission, b.pressure_transmission);
        assert_eq!(a.transmitted_fraction, b.transmitted_fraction);
        assert_eq!(a.absorbed_fraction, 0.0);
    }
}

#[test]
fn oblique_plate_matches_all_six_independent_boundary_conditions() {
    let incident = fluid(1000.0, 1500.0);
    for (rho_s, cp, cs) in [(2700.0, 6000.0, 3000.0), (2200.0, 4000.0, 3200.0)] {
        let plate = solid(rho_s, cp, cs);
        for (rho_r, cr) in [(1000.0, 1500.0), (650.0, 1900.0), (1300.0, 900.0)] {
            let receiver = fluid(rho_r, cr);
            for thickness in [1e-8, 0.007, 0.2, 1.4] {
                let stack =
                    ElasticPlateLayers::new(&incident, &plate, &receiver, 16000.0, thickness, &[])
                        .unwrap();
                let [cp, cs] = stack.plate_phase_speeds_m_s();
                for sine in [0.03_f64, 0.15, 0.3, 0.48, 0.65, 0.9] {
                    let cosine = (1.0 - sine * sine).sqrt();
                    let actual = stack.scattering(cosine).unwrap();
                    let expected = boundary_oracle(
                        [1000.0, 1500.0],
                        [rho_s, cp, cs],
                        [rho_r, cr],
                        16000.0,
                        thickness,
                        cosine,
                    );
                    near(actual.pressure_reflection, expected[0], 7e-11);
                    near(actual.pressure_transmission, expected[5], 7e-11);
                    near(
                        actual.plate_front_pressure,
                        C64::ONE + actual.pressure_reflection,
                        8e-14,
                    );
                    near_real(
                        actual.reflected_fraction + actual.transmitted_fraction,
                        1.0,
                        2e-13,
                    );
                    assert_eq!(actual.absorbed_fraction, 0.0);
                }
            }
        }
    }
}

#[test]
fn all_three_critical_transitions_preserve_outgoing_branches_and_power() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = fluid(900.0, 2200.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let stack = ElasticPlateLayers::new(&incident, &plate, &receiver, 8000.0, 0.04, &[]).unwrap();
    let [cp, cs] = stack.plate_phase_speeds_m_s();
    for (mode, speed) in [cp, cs, 2200.0].into_iter().enumerate() {
        let critical = (1.0 - (1500.0 / speed).powi(2)).sqrt();
        let mut previous = None;
        for offset in [-16_i64, -1, 0, 1, 16] {
            let cosine = f64::from_bits((critical.to_bits() as i64 + offset) as u64);
            let actual = stack.scattering(cosine).unwrap();
            let kn = [
                actual.longitudinal_normal_wavenumber_rad_m,
                actual.shear_normal_wavenumber_rad_m,
                actual.receiving_normal_wavenumber_rad_m,
            ][mode];
            if offset == -16 {
                assert!(kn.im < 0.0 && kn.re == 0.0);
            }
            if offset == 16 {
                assert!(kn.re > 0.0 && kn.im == 0.0);
            }
            if mode == 2 && kn.re == 0.0 {
                assert_eq!(actual.transmitted_fraction, 0.0);
            }
            near_real(
                actual.reflected_fraction + actual.transmitted_fraction,
                1.0,
                1e-13,
            );
            if let Some(prior) = previous {
                near(actual.pressure_transmission, prior, 5e-6);
            }
            previous = Some(actual.pressure_transmission);
        }
    }
}

#[test]
fn finite_evanescent_plate_tunnels_without_artificial_absorption_or_growth() {
    let incident = fluid(1000.0, 1500.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let cosine = 0.6;
    let omega = 15000.0;
    let halfspace = ElasticTransmittingLayers::new(&incident, &plate, omega, &[])
        .unwrap()
        .scattering(cosine)
        .unwrap();
    let mut previous = 1.0;
    for thickness in [0.0001, 0.1, 1.0, 10.0, 100.0, 100000.0] {
        let actual = ElasticPlateLayers::new(&incident, &plate, &incident, omega, thickness, &[])
            .unwrap()
            .scattering(cosine)
            .unwrap();
        assert!(actual.longitudinal_normal_wavenumber_rad_m.im < 0.0);
        assert!(actual.shear_normal_wavenumber_rad_m.im < 0.0);
        assert!(actual.transmitted_fraction <= previous + 2e-14);
        assert_eq!(actual.absorbed_fraction, 0.0);
        near_real(
            actual.reflected_fraction + actual.transmitted_fraction,
            1.0,
            1e-13,
        );
        if thickness == 0.0001 {
            assert!(actual.transmitted_fraction > 0.999);
        }
        if thickness == 10.0 {
            assert!(actual.transmitted_fraction > 0.0 && actual.transmitted_fraction < 1e-45);
            assert!(actual.pressure_transmission.abs() > 0.0);
        }
        if thickness >= 100.0 {
            near(
                actual.pressure_reflection,
                halfspace.pressure_reflection,
                8e-14,
            );
            assert_eq!(actual.transmitted_fraction, 0.0);
        }
        previous = actual.transmitted_fraction;
    }
}

#[test]
fn air_loaded_plate_retains_the_narrow_flexural_tunneling_resonance() {
    let air = fluid(1.2, 343.0);
    let plate = solid(7800.0, 5900.0, 3200.0);
    let omega = 120046.33170621305;
    let mut peak = 0.0;
    for frequency_ratio in [0.99, 1.0, 1.01] {
        let stack =
            ElasticPlateLayers::new(&air, &plate, &air, omega * frequency_ratio, 0.001, &[])
                .unwrap();
        let actual = stack.scattering(0.6).unwrap();
        let [cp, cs] = stack.plate_phase_speeds_m_s();
        let expected = boundary_oracle(
            [1.2, 343.0],
            [7800.0, cp, cs],
            [1.2, 343.0],
            omega * frequency_ratio,
            0.001,
            0.6,
        );
        // The independent six-wave basis is poorly conditioned at this
        // high-contrast subwavelength resonance; retain a physical tolerance.
        near(actual.pressure_reflection, expected[0], 3e-8);
        near(actual.pressure_transmission, expected[5], 3e-8);
        assert!(actual.longitudinal_normal_wavenumber_rad_m.im < 0.0);
        assert!(actual.shear_normal_wavenumber_rad_m.im < 0.0);
        near_real(
            actual.reflected_fraction + actual.transmitted_fraction,
            1.0,
            2e-13,
        );
        if frequency_ratio == 1.0 {
            peak = actual.transmitted_fraction;
            assert!(peak > 0.99999);
        } else {
            assert!(actual.transmitted_fraction < 0.006);
        }
    }
    assert!(peak > 0.99999);
}

#[test]
fn reciprocal_fluid_loading_and_entrance_gap_preserve_plate_phase_and_power() {
    let front = fluid(1000.0, 1500.0);
    let back = fluid(800.0, 1200.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let omega = 23000.0;
    let cosine = 0.91_f64;
    let reverse_cosine = (1.0 - (1200.0_f64 / 1500.0).powi(2) * (1.0 - cosine * cosine)).sqrt();
    let forward = ElasticPlateLayers::new(&front, &plate, &back, omega, 0.13, &[])
        .unwrap()
        .scattering(cosine)
        .unwrap();
    let reverse = ElasticPlateLayers::new(&back, &plate, &front, omega, 0.13, &[])
        .unwrap()
        .scattering(reverse_cosine)
        .unwrap();
    near_real(
        forward.transmitted_fraction,
        reverse.transmitted_fraction,
        2e-13,
    );
    let load = (1000.0 * 1500.0 * reverse_cosine) / (800.0 * 1200.0 * cosine);
    near(
        forward.pressure_transmission.scale(load),
        reverse.pressure_transmission,
        2e-12,
    );
    let gap = 0.031;
    let delayed = ElasticPlateLayers::new(
        &front,
        &plate,
        &back,
        omega,
        0.13,
        &[Layer::FluidGap { thickness_m: gap }],
    )
    .unwrap()
    .scattering(cosine)
    .unwrap();
    let phase = -omega * gap * cosine / 1500.0;
    let delay = C64::new(phase.cos(), phase.sin());
    near(
        delayed.pressure_reflection,
        forward.pressure_reflection * delay * delay,
        2e-14,
    );
    near(
        delayed.pressure_transmission,
        forward.pressure_transmission * delay,
        2e-14,
    );
    near_real(
        delayed.transmitted_fraction,
        forward.transmitted_fraction,
        2e-14,
    );
    near(
        delayed.plate_front_pressure,
        forward.plate_front_pressure * delay,
        2e-14,
    );
}

#[test]
fn porous_entrance_loss_remains_separate_from_elastic_mode_conversion() {
    let incident = fluid(1.2, 343.0);
    let plate = solid(7800.0, 5900.0, 3200.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.025,
        },
        Layer::FluidGap { thickness_m: 0.01 },
    ];
    let stack =
        ElasticPlateLayers::new(&incident, &plate, &incident, 8000.0, 0.00015, &layers).unwrap();
    for cosine in [1.0, 0.99, 0.8, 0.4] {
        let actual = stack.scattering(cosine).unwrap();
        assert!(actual.absorbed_fraction > 0.0);
        assert!(actual.transmitted_fraction > 0.0);
        near_real(
            actual.reflected_fraction + actual.transmitted_fraction + actual.absorbed_fraction,
            1.0,
            8e-14,
        );
    }
}

#[test]
fn plate_rejects_unadmitted_geometry_angles_and_anisotropic_sources() {
    let incident = fluid(1000.0, 1500.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    for thickness in [-0.1, f64::NAN, f64::INFINITY] {
        assert!(
            ElasticPlateLayers::new(&incident, &plate, &incident, 8000.0, thickness, &[]).is_err()
        );
    }
    let generic = ElasticAcousticState::try_new_mandel(
        *plate.stiffness_mandel_pa(),
        plate.density_kg_m3(),
        ContentHash([25; 32]),
        plate.maximum_absolute_mandel_strain(),
    )
    .unwrap();
    assert!(ElasticPlateLayers::new(&incident, &generic, &incident, 8000.0, 0.1, &[]).is_err());
    let stack = ElasticPlateLayers::new(&incident, &plate, &incident, 8000.0, 0.1, &[]).unwrap();
    for cosine in [0.0, -0.5, 1.01, f64::NAN, f64::INFINITY] {
        assert!(stack.scattering(cosine).is_err());
    }
}

// One original porous entrance layer, followed by an incident-fluid gap and
// the finite plate. The independent six-wave boundary solution supplies the
// load. Conventional complex propagation and spatial quadrature then evaluate
// the local viscous/thermal loss; production face states and analytic integrals
// are not used. Fields are normalized to one peak incident pascal.
fn independent_entrance_heat(
    stack: &ElasticPlateLayers<'_>,
    cosine: f64,
    samples: usize,
) -> [f64; 2] {
    let prepared = stack.prepared_layers();
    let incident = prepared.fluid().acoustic();
    let receiver = stack.receiving_state();
    let [cp, cs] = stack.plate_phase_speeds_m_s();
    let omega = prepared.angular_frequency_rad_s();
    let bare = boundary_oracle(
        [incident.density_kg_m3(), incident.sound_speed_m_s()],
        [stack.plate_state().density_kg_m3(), cp, cs],
        [receiver.density_kg_m3(), receiver.sound_speed_m_s()],
        omega,
        stack.thickness_m(),
        cosine,
    );
    let Layer::Porous { thickness_m: d, .. } = prepared.layer(0).unwrap() else {
        panic!("oracle requires the original porous entrance layer");
    };
    let Layer::FluidGap { thickness_m: gap } = prepared.layer(1).unwrap() else {
        panic!("oracle requires the original incident-fluid gap");
    };
    let phi = omega * gap * cosine / incident.sound_speed_m_s();
    let jsin = C64::new(0.0, phi.sin());
    let p = (C64::ONE + bare[0]).scale(phi.cos()) + jsin * (C64::ONE - bare[0]);
    let q = (C64::ONE - bare[0]).scale(phi.cos()) + jsin * (C64::ONE + bare[0]);
    let response = prepared.material_response(0).unwrap();
    let rho = response.dynamic_density_kg_m3();
    let compliance = response.dynamic_bulk_modulus_pa().recip();
    let k = response.characteristic().1;
    let zref = incident.characteristic_impedance_pa_s_m() / cosine;
    let kt = omega / incident.sound_speed_m_s() * (1.0 - cosine * cosine).sqrt();
    let mut kn = (k * k - C64::from_re(kt * kt)).sqrt();
    if kn.im > 0.0 {
        kn = -kn;
    }
    let ratio = rho.scale(omega / zref) / kn;
    let cos_sin = |z: C64| {
        (
            C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh()),
            C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
        )
    };
    let (cos, sin) = cos_sin(kn.scale(d));
    let jsin = C64::new(0.0, 1.0) * sin;
    let p_front = cos * p + ratio * jsin * q;
    let q_front = cos * q + jsin * p / ratio;
    let amplitude = C64::from_re(2.0) / (p_front + q_front);
    let (p, q) = (p * amplitude, q * amplitude);
    let mut heat = [0.0; 2];
    for sample in 0..samples {
        let z = (sample as f64 + 0.5) * d / samples as f64;
        let (cos, sin) = cos_sin(kn.scale(z));
        let jsin = C64::new(0.0, 1.0) * sin;
        let pressure = cos * p + ratio * jsin * q;
        let normal_velocity = (cos * q + jsin * p / ratio).scale(1.0 / zref);
        let tangential_velocity = (pressure / rho).scale(kt / omega);
        heat[0] +=
            omega * zref * -rho.im * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
        heat[1] += omega * zref * -compliance.im * pressure.norm_sq();
    }
    heat.map(|value| value * d / samples as f64)
}

#[test]
fn finite_plate_heat_matches_independent_boundary_and_local_loss() {
    let incident = fluid(1.2, 340.0);
    let receiver = fluid(1.0, 300.0);
    let plate = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.025,
        },
        Layer::FluidGap { thickness_m: 0.012 },
    ];
    let budget = LayerAbsorptionBudget { max_layers: 2 };
    for thickness in [0.003, 0.035, 0.3] {
        let stack =
            ElasticPlateLayers::new(&incident, &plate, &receiver, 4000.0, thickness, &layers)
                .unwrap();
        for sine in [0.2_f64, 0.65, 0.985] {
            let cosine = (1.0 - sine * sine).sqrt();
            let response = stack.scattering(cosine).unwrap();
            let actual = stack
                .layer_absorption(cosine, 5.0, budget, || false)
                .unwrap();
            let expected = independent_entrance_heat(&stack, cosine, 4096);
            for (observed, expected) in [
                actual.layers()[0].viscous_absorbed_fraction,
                actual.layers()[0].thermal_absorbed_fraction,
            ]
            .into_iter()
            .zip(expected)
            {
                assert!(observed > 0.0);
                assert!(
                    (observed / expected - 1.0).abs() < 2e-7,
                    "d={thickness}, sin(theta)={sine}, observed={observed}, expected={expected}"
                );
            }
            assert_eq!(actual.layers().len(), 2);
            assert_eq!(actual.layers()[1].absorbed_power_w, 0.0);
            assert_eq!(actual.reflected_fraction, response.reflected_fraction);
            assert_eq!(actual.transmitted_fraction, response.transmitted_fraction);
            assert!(actual.transmitted_fraction > 0.0);
            near_real(actual.absorbed_fraction, response.absorbed_fraction, 3e-13);
            assert!(actual.balance_residual_w.abs() < 1e-12);
            if sine > 0.98 {
                assert!(response.longitudinal_normal_wavenumber_rad_m.im < 0.0);
                assert!(response.shear_normal_wavenumber_rad_m.im < 0.0);
            }
        }
    }
}

#[test]
fn zero_and_normal_plate_limits_keep_original_layer_heating() {
    let incident = fluid(1.2, 340.0);
    let receiver = fluid(1.0, 300.0);
    let plate = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.025,
        },
        Layer::FluidGap { thickness_m: 0.012 },
    ];
    let omega = 4000.0;
    let budget = LayerAbsorptionBudget { max_layers: 2 };
    let zero = ElasticPlateLayers::new(&incident, &plate, &receiver, omega, 0.0, &layers).unwrap();
    let original = TransmittingLayers::new(&incident, &receiver, omega, &layers).unwrap();
    for cosine in [1.0, 0.8, 0.2] {
        let actual = zero
            .layer_absorption(cosine, 5.0, budget, || false)
            .unwrap();
        let expected = original
            .layer_absorption(cosine, 5.0, budget, || false)
            .unwrap();
        assert_eq!(actual, expected);
    }
    let cp = zero.plate_phase_speeds_m_s()[0];
    let longitudinal = FluidAcousticState::try_new(plate.density_kg_m3(), cp).unwrap();
    let mut heat = [0.0; 2];
    for (index, thickness) in [0.5, 1.0]
        .map(|multiple| multiple * core::f64::consts::PI * cp / omega)
        .into_iter()
        .enumerate()
    {
        let elastic =
            ElasticPlateLayers::new(&incident, &plate, &receiver, omega, thickness, &layers)
                .unwrap();
        let reference = TransmittingLayers::new(
            &incident,
            &receiver,
            omega,
            &[
                layers[0],
                layers[1],
                Layer::FluidLayer {
                    fluid: &longitudinal,
                    thickness_m: thickness,
                },
            ],
        )
        .unwrap();
        let actual = elastic
            .layer_absorption(1.0, 5.0, budget, || false)
            .unwrap();
        let expected = reference
            .layer_absorption(1.0, 5.0, LayerAbsorptionBudget { max_layers: 3 }, || false)
            .unwrap();
        near_real(
            actual.layers()[0].viscous_absorbed_fraction,
            expected.layers()[0].viscous_absorbed_fraction,
            4e-14,
        );
        near_real(
            actual.layers()[0].thermal_absorbed_fraction,
            expected.layers()[0].thermal_absorbed_fraction,
            4e-14,
        );
        near_real(
            actual.transmitted_fraction,
            expected.transmitted_fraction,
            4e-14,
        );
        assert_eq!(expected.layers()[2].absorbed_power_w, 0.0);
        heat[index] = actual.absorbed_fraction;
    }
    // The lossless plate's quarter/half-wave loading changes heat in the
    // actual porous material; no plate absorption is needed to obtain it.
    assert!((heat[0] - heat[1]).abs() > 1e-3, "{heat:?}");
}

#[test]
fn thin_porous_heat_keeps_evanescent_plate_tunneling_below_roundoff() {
    let incident = fluid(1.2, 340.0);
    let plate = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let omega = 4000.0;
    let d = 1e-20;
    let cosine = 0.2_f64;
    let stack = ElasticPlateLayers::new(
        &incident,
        &plate,
        &incident,
        omega,
        0.02,
        &[Layer::Porous {
            law: &foam,
            thickness_m: d,
        }],
    )
    .unwrap();
    let actual = stack
        .layer_absorption(cosine, 5.0, LayerAbsorptionBudget { max_layers: 1 }, || {
            false
        })
        .unwrap();
    let [cp, cs] = stack.plate_phase_speeds_m_s();
    let bare = boundary_oracle(
        [1.2, 340.0],
        [30.0, cp, cs],
        [1.2, 340.0],
        omega,
        stack.thickness_m(),
        cosine,
    );
    let response = stack.prepared_layers().material_response(0).unwrap();
    let rho = response.dynamic_density_kg_m3();
    let zref = incident.acoustic().characteristic_impedance_pa_s_m() / cosine;
    let pressure = C64::ONE + bare[0];
    let normal_velocity = (C64::ONE - bare[0]).scale(1.0 / zref);
    let tangential_velocity = (pressure / rho).scale((1.0 - cosine * cosine).sqrt() / 340.0);
    // First-order heat is thickness times the local loss at the original
    // plate's exact front trace. This retains the oblique tangential channel.
    let viscous =
        omega * d * zref * -rho.im * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
    let thermal =
        omega * d * zref * -response.dynamic_bulk_modulus_pa().recip().im * pressure.norm_sq();
    assert!((actual.layers()[0].viscous_absorbed_fraction / viscous - 1.0).abs() < 1e-11);
    assert!((actual.layers()[0].thermal_absorbed_fraction / thermal - 1.0).abs() < 1e-11);
    assert!(actual.absorbed_fraction > 0.0 && actual.absorbed_fraction < f64::EPSILON);
    assert!(actual.absorbed_power_w > 0.0);
    assert!(actual.transmitted_fraction > 0.0);
    let waves = stack.scattering(cosine).unwrap();
    assert!(waves.longitudinal_normal_wavenumber_rad_m.im < 0.0);
    assert!(waves.shear_normal_wavenumber_rad_m.im < 0.0);
}

#[test]
fn finite_plate_heat_obeys_entrance_budget_and_atomic_cancellation() {
    let incident = fluid(1.2, 340.0);
    let plate = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.025,
        },
        Layer::FluidGap { thickness_m: 0.012 },
    ];
    let stack =
        ElasticPlateLayers::new(&incident, &plate, &incident, 4000.0, 0.02, &layers).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 2 };
    let expected = stack.layer_absorption(0.8, 5.0, budget, || false).unwrap();
    for cancel_on in 1..=2 * layers.len() + 2 {
        let mut calls = 0;
        assert_eq!(
            stack.layer_absorption(0.8, 5.0, budget, || {
                calls += 1;
                calls == cancel_on
            }),
            Err(LayerAbsorptionError::Cancelled)
        );
    }
    assert_eq!(
        stack.layer_absorption(0.8, 5.0, budget, || false).unwrap(),
        expected
    );
    assert!(
        stack
            .layer_absorption(0.8, 5.0, LayerAbsorptionBudget { max_layers: 1 }, || false)
            .is_err()
    );
    assert!(
        stack
            .layer_absorption(0.8, f64::NAN, budget, || false)
            .is_err()
    );
    let zero = stack.layer_absorption(0.8, 0.0, budget, || false).unwrap();
    assert_eq!(zero.absorbed_fraction, expected.absorbed_fraction);
    assert_eq!(zero.absorbed_power_w, 0.0);
    let lossless =
        ElasticPlateLayers::new(&incident, &plate, &incident, 4000.0, 0.02, &[]).unwrap();
    let actual = lossless
        .layer_absorption(0.2, 5.0, LayerAbsorptionBudget { max_layers: 0 }, || false)
        .unwrap();
    assert!(actual.layers().is_empty());
    assert_eq!(actual.absorbed_power_w, 0.0);
    assert!(actual.transmitted_power_w > 0.0);
}
