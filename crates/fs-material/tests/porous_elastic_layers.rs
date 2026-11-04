//! Analytic fluid/solid refraction controls, using original material states.
use fs_blake3::ContentHash;
use fs_material::elastic::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::layers::{
    ElasticTransmittingLayers, Layer, LayerAbsorptionBudget, LayerAbsorptionError,
    TransmittingLayers,
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
        ContentHash([23; 32]),
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

// Independent full boundary system, without the production shear null vector:
// continuity of normal stress, continuity of normal velocity, zero shear stress.
// Unknowns are pressure R and P/S velocity relative to incident velocity.
fn boundary_oracle(rho_f: f64, cf: f64, rho_s: f64, cp: f64, cs: f64, ci: f64) -> [C64; 3] {
    let sine = (1.0 - ci * ci).sqrt();
    let sp = cp * sine / cf;
    let ss = cs * sine / cf;
    let outgoing = |s: f64| {
        if s <= 1.0 {
            C64::from_re((1.0 - s * s).sqrt())
        } else {
            C64::new(0.0, -(s * s - 1.0).sqrt())
        }
    };
    let cp_cos = outgoing(sp);
    let cs_cos = outgoing(ss);
    let a = 1.0 - 2.0 * ss * ss;
    let zp = rho_s * cp / (rho_f * cf);
    let zs = rho_s * cs / (rho_f * cf);
    let mut system = [
        [
            C64::ONE,
            C64::from_re(-zp * a),
            cs_cos.scale(2.0 * zs * ss),
            -C64::ONE,
        ],
        [
            C64::from_re(-ci),
            -cp_cos,
            C64::from_re(ss),
            C64::from_re(-ci),
        ],
        [
            C64::ZERO,
            cp_cos.scale(2.0 * ss),
            C64::from_re(a),
            C64::ZERO,
        ],
    ];
    for column in 0..3 {
        let pivot = (column..3)
            .max_by(|i, j| {
                system[*i][column]
                    .abs()
                    .total_cmp(&system[*j][column].abs())
            })
            .unwrap();
        system.swap(column, pivot);
        let divisor = system[column][column];
        for entry in column..4 {
            system[column][entry] = system[column][entry] / divisor;
        }
        for row in 0..3 {
            if row == column {
                continue;
            }
            let multiplier = system[row][column];
            for entry in column..4 {
                system[row][entry] = system[row][entry] - multiplier * system[column][entry];
            }
        }
    }
    [system[0][3], system[1][3], system[2][3]]
}

#[test]
fn normal_incidence_is_the_original_longitudinal_impedance_limit() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = solid(2700.0, 6000.0, 3000.0);
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 8000.0, &[]).unwrap();
    assert!(core::ptr::eq(stack.receiving_state(), &receiver));
    let response = stack.scattering(1.0).unwrap();
    let zp = receiver.density_kg_m3() * stack.receiving_phase_speeds_m_s()[0];
    let zf = incident.acoustic().characteristic_impedance_pa_s_m();
    let reflection = (zp - zf) / (zp + zf);
    near(
        response.pressure_reflection,
        C64::from_re(reflection),
        3e-15,
    );
    near(
        response.surface_pressure_transmission,
        C64::from_re(2.0 * zp / (zp + zf)),
        3e-15,
    );
    near(
        response.longitudinal_velocity_transmission,
        C64::from_re(2.0 * zf / (zp + zf)),
        3e-15,
    );
    assert_eq!(response.shear_velocity_transmission, C64::ZERO);
    assert_eq!(response.shear_transmitted_fraction, 0.0);
    near_real(
        response.longitudinal_transmitted_fraction,
        4.0 * zf * zp / (zf + zp).powi(2),
        4e-15,
    );
    assert_eq!(response.absorbed_fraction, 0.0);
}

#[test]
fn oblique_modes_satisfy_the_independent_three_boundary_equations() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = solid(2700.0, 6000.0, 3000.0);
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 8000.0, &[]).unwrap();
    let [cp, cs] = stack.receiving_phase_speeds_m_s();
    // Both modes propagate, only SV propagates, then both are evanescent.
    for sine in [0.08_f64, 0.2, 0.35, 0.47, 0.6, 0.95] {
        let cosine = (1.0 - sine * sine).sqrt();
        let actual = stack.scattering(cosine).unwrap();
        let expected = boundary_oracle(1000.0, 1500.0, 2700.0, cp, cs, cosine);
        near(actual.pressure_reflection, expected[0], 8e-14);
        near(
            actual.longitudinal_velocity_transmission,
            expected[1],
            8e-14,
        );
        near(actual.shear_velocity_transmission, expected[2], 8e-14);
        near(
            actual.surface_pressure_transmission,
            C64::ONE + expected[0],
            8e-14,
        );
        near_real(
            actual.reflected_fraction
                + actual.longitudinal_transmitted_fraction
                + actual.shear_transmitted_fraction,
            1.0,
            2e-14,
        );
        if sine > 0.25 {
            assert_eq!(actual.longitudinal_transmitted_fraction, 0.0);
            assert!(actual.longitudinal_normal_wavenumber_rad_m.im < 0.0);
            assert!(actual.longitudinal_velocity_transmission.abs() > 0.0);
        }
        if sine > 0.5 {
            assert_eq!(actual.shear_transmitted_fraction, 0.0);
            assert!(actual.shear_normal_wavenumber_rad_m.im < 0.0);
            assert!(actual.shear_velocity_transmission.abs() > 0.0);
        } else {
            assert!(actual.shear_transmitted_fraction > 0.0);
        }
    }
}

#[test]
fn both_critical_transitions_keep_outgoing_fields_and_close_power() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = solid(2700.0, 6000.0, 3000.0);
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 8000.0, &[]).unwrap();
    for (mode, speed) in stack.receiving_phase_speeds_m_s().into_iter().enumerate() {
        let ratio = 1500.0 / speed;
        let critical = (1.0 - ratio * ratio).sqrt();
        for offset in [-8_i64, 8] {
            let cosine = f64::from_bits((critical.to_bits() as i64 + offset) as u64);
            let response = stack.scattering(cosine).unwrap();
            let kn = [
                response.longitudinal_normal_wavenumber_rad_m,
                response.shear_normal_wavenumber_rad_m,
            ][mode];
            let power = [
                response.longitudinal_transmitted_fraction,
                response.shear_transmitted_fraction,
            ][mode];
            if offset > 0 {
                assert!(kn.re > 0.0 && kn.im == 0.0);
                assert!(power > 0.0);
            } else {
                assert!(kn.re == 0.0 && kn.im < 0.0);
                assert_eq!(power, 0.0);
            }
            near_real(
                response.reflected_fraction
                    + response.longitudinal_transmitted_fraction
                    + response.shear_transmitted_fraction,
                1.0,
                3e-14,
            );
        }
    }
}

#[test]
fn forty_five_degree_shear_and_rayleigh_pressure_zero_are_not_poles() {
    let incident = fluid(1.0, 0.5);
    // An auxetic solid allows theta_S=45deg while P still propagates.
    let auxetic = solid(2.0, 1.25, 1.0);
    let stack = ElasticTransmittingLayers::new(&incident, &auxetic, 1.0, &[]).unwrap();
    let shear_speed = stack.receiving_phase_speeds_m_s()[1];
    let sine = 0.5 * core::f64::consts::FRAC_1_SQRT_2 / shear_speed;
    let response = stack.scattering((1.0 - sine * sine).sqrt()).unwrap();
    assert!(response.longitudinal_transmitted_fraction < 1e-27);
    assert!(response.shear_transmitted_fraction > 0.1);
    near_real(
        response.reflected_fraction + response.shear_transmitted_fraction,
        1.0,
        4e-15,
    );

    // For lambda=mu (nu=1/4), the traction-free Rayleigh root is
    // (c_R/c_S)^2=2-2/sqrt(3). Its surface normal stress is zero.
    let law = IsotropicElastic::new(5.0, 0.25, 0.001).unwrap();
    let receiver = ElasticAcousticState::from_isotropic(&law, 2.0, ContentHash([24; 32])).unwrap();
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 1.0, &[]).unwrap();
    let rayleigh = stack.receiving_phase_speeds_m_s()[1] * (2.0 - 2.0 / 3.0_f64.sqrt()).sqrt();
    let sine = 0.5 / rayleigh;
    let response = stack.scattering((1.0 - sine * sine).sqrt()).unwrap();
    near(response.pressure_reflection, -C64::ONE, 3e-14);
    near(response.surface_pressure_transmission, C64::ZERO, 3e-14);
    assert!(response.longitudinal_velocity_transmission.abs() > 0.0);
    assert!(response.shear_velocity_transmission.abs() > 0.0);
    assert_eq!(response.longitudinal_transmitted_fraction, 0.0);
    assert_eq!(response.shear_transmitted_fraction, 0.0);
}

#[test]
fn finite_layers_share_the_original_transfer_and_loss_owner() {
    let incident = fluid(1.2, 340.0);
    let receiver = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let film = FluidAcousticState::try_new(2.0, 500.0).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.025,
        },
        Layer::FluidLayer {
            fluid: &film,
            thickness_m: 0.012,
        },
        Layer::FluidGap { thickness_m: 0.02 },
    ];
    let elastic = ElasticTransmittingLayers::new(&incident, &receiver, 4000.0, &layers).unwrap();
    let receiving_fluid = FluidAcousticState::try_new(
        receiver.density_kg_m3(),
        elastic.receiving_phase_speeds_m_s()[0],
    )
    .unwrap();
    let fluid_stack =
        TransmittingLayers::new(&incident, &receiving_fluid, 4000.0, &layers).unwrap();
    let a = elastic.scattering(1.0).unwrap();
    let b = fluid_stack.scattering(1.0).unwrap();
    near(a.pressure_reflection, b.pressure_reflection, 2e-14);
    near(
        a.surface_pressure_transmission,
        b.pressure_transmission,
        2e-14,
    );
    near_real(
        a.longitudinal_transmitted_fraction,
        b.transmitted_fraction,
        2e-14,
    );
    near_real(a.absorbed_fraction, b.absorbed_fraction, 2e-14);
    assert!(a.absorbed_fraction > 0.1);
    for cosine in [0.95, 0.8, 0.5, 0.1] {
        let response = elastic.scattering(cosine).unwrap();
        assert!(response.absorbed_fraction > 0.0);
        near_real(
            response.reflected_fraction
                + response.longitudinal_transmitted_fraction
                + response.shear_transmitted_fraction
                + response.absorbed_fraction,
            1.0,
            3e-14,
        );
    }
}

#[test]
fn transparent_incident_gap_preserves_modal_power_and_advances_phase() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = solid(2700.0, 6000.0, 3000.0);
    let omega = 8000.0;
    let distance = 0.07;
    let bare = ElasticTransmittingLayers::new(&incident, &receiver, omega, &[]).unwrap();
    let gap = ElasticTransmittingLayers::new(
        &incident,
        &receiver,
        omega,
        &[Layer::FluidGap {
            thickness_m: distance,
        }],
    )
    .unwrap();
    for cosine in [0.99, 0.95, 0.8] {
        let a = bare.scattering(cosine).unwrap();
        let b = gap.scattering(cosine).unwrap();
        let angle = -omega * distance * cosine / 1500.0;
        let phase = C64::new(angle.cos(), angle.sin());
        near(
            b.pressure_reflection,
            a.pressure_reflection * phase * phase,
            6e-15,
        );
        near(
            b.surface_pressure_transmission,
            a.surface_pressure_transmission * phase,
            6e-15,
        );
        near(
            b.longitudinal_velocity_transmission,
            a.longitudinal_velocity_transmission * phase,
            6e-15,
        );
        near(
            b.shear_velocity_transmission,
            a.shear_velocity_transmission * phase,
            6e-15,
        );
        near_real(
            b.longitudinal_transmitted_fraction,
            a.longitudinal_transmitted_fraction,
            6e-15,
        );
        near_real(
            b.shear_transmitted_fraction,
            a.shear_transmitted_fraction,
            6e-15,
        );
    }
}

#[test]
fn general_elastic_states_and_invalid_incidence_are_not_silently_projected() {
    let incident = fluid(1000.0, 1500.0);
    let receiver = solid(2700.0, 6000.0, 3000.0);
    assert!(receiver.is_explicitly_isotropic());
    let general = ElasticAcousticState::try_new_mandel(
        *receiver.stiffness_mandel_pa(),
        receiver.density_kg_m3(),
        receiver.material_state_identity(),
        receiver.maximum_absolute_mandel_strain(),
    )
    .unwrap();
    assert!(!general.is_explicitly_isotropic());
    assert!(ElasticTransmittingLayers::new(&incident, &general, 8000.0, &[]).is_err());
    let cubic =
        ElasticAcousticState::from_cubic(2700.0, 90e9, 20e9, 35e9, ContentHash([25; 32]), 0.001)
            .unwrap();
    assert!(ElasticTransmittingLayers::new(&incident, &cubic, 8000.0, &[]).is_err());
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 8000.0, &[]).unwrap();
    for cosine in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(stack.scattering(cosine).is_err());
    }
}

// Independent P/S boundary equations followed by conventional transfer matrices
// and spatial midpoint quadrature of the local dissipated power. Growing
// trigonometric functions are confined to this moderate-thickness oracle.
// Neither the production layer-face observer nor its analytic loss integral
// participates. The exterior incident pressure is one peak pascal.
fn independent_layer_heat(
    stack: &ElasticTransmittingLayers<'_>,
    cosine: f64,
    samples: usize,
) -> Vec<[f64; 2]> {
    let prepared = stack.prepared_layers();
    let fluid = prepared.fluid().acoustic();
    let omega = prepared.angular_frequency_rad_s();
    let zref = fluid.characteristic_impedance_pa_s_m() / cosine;
    let kt = omega / fluid.sound_speed_m_s() * (1.0 - cosine * cosine).sqrt();
    let wave = |index| {
        let (rho, k, d) = match prepared.layer(index).unwrap() {
            Layer::Porous { thickness_m, .. } => {
                let response = prepared.material_response(index).unwrap();
                (
                    response.dynamic_density_kg_m3(),
                    response.characteristic().1,
                    thickness_m,
                )
            }
            Layer::FluidLayer { thickness_m, .. } => {
                let response = prepared.fluid_response(index).unwrap();
                (
                    C64::from_re(response.density_kg_m3()),
                    C64::from_re(omega / response.sound_speed_m_s()),
                    thickness_m,
                )
            }
            Layer::FluidGap { thickness_m } => (
                C64::from_re(fluid.density_kg_m3()),
                C64::from_re(omega / fluid.sound_speed_m_s()),
                thickness_m,
            ),
        };
        let mut kn = (k * k - C64::from_re(kt * kt)).sqrt();
        if kn.im > 0.0 {
            kn = -kn;
        }
        (kn, rho.scale(omega / zref) / kn, d)
    };
    let cos_sin = |z: C64| {
        (
            C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh()),
            C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
        )
    };
    let mut matrix = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    for index in 0..prepared.layer_count() {
        let (kn, ratio, d) = wave(index);
        let (cos, sin) = cos_sin(kn.scale(d));
        let jsin = C64::new(0.0, 1.0) * sin;
        let layer = [[cos, ratio * jsin], [jsin / ratio, cos]];
        let mut product = [[C64::ZERO; 2]; 2];
        for row in 0..2 {
            for col in 0..2 {
                product[row][col] = matrix[row][0] * layer[0][col] + matrix[row][1] * layer[1][col];
            }
        }
        matrix = product;
    }
    let [cp, cs] = stack.receiving_phase_speeds_m_s();
    let bare = boundary_oracle(
        fluid.density_kg_m3(),
        fluid.sound_speed_m_s(),
        stack.receiving_state().density_kg_m3(),
        cp,
        cs,
        cosine,
    );
    let rear = [C64::ONE + bare[0], C64::ONE - bare[0]];
    let entrance = [
        matrix[0][0] * rear[0] + matrix[0][1] * rear[1],
        matrix[1][0] * rear[0] + matrix[1][1] * rear[1],
    ];
    let amplitude = C64::from_re(2.0) / (entrance[0] + entrance[1]);
    let mut p = entrance[0] * amplitude;
    let mut q = entrance[1] * amplitude; // q = zref * normal particle velocity
    let mut result = Vec::with_capacity(prepared.layer_count());
    for index in 0..prepared.layer_count() {
        let (kn, ratio, d) = wave(index);
        let mut heat = [0.0; 2];
        if let Some(response) = prepared.material_response(index) {
            let rho = response.dynamic_density_kg_m3();
            let compliance = response.dynamic_bulk_modulus_pa().recip();
            for sample in 0..samples {
                let z = (sample as f64 + 0.5) * d / samples as f64;
                let (cos, sin) = cos_sin(kn.scale(z));
                let jsin = C64::new(0.0, 1.0) * sin;
                let pressure = cos * p - ratio * jsin * q;
                let normal_velocity = (cos * q - jsin * p / ratio).scale(1.0 / zref);
                let tangential_velocity = (pressure / rho).scale(kt / omega);
                // Divide omega/2 times the local quadratic loss by the
                // incident intensity 1/(2*zref). Both velocity components
                // dissipate through the ORIGINAL dynamic mass density.
                heat[0] += omega
                    * zref
                    * -rho.im
                    * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
                heat[1] += omega * zref * -compliance.im * pressure.norm_sq();
            }
            heat = heat.map(|value| value * d / samples as f64);
        }
        result.push(heat);
        let (cos, sin) = cos_sin(kn.scale(d));
        let jsin = C64::new(0.0, 1.0) * sin;
        (p, q) = (cos * p - ratio * jsin * q, cos * q - jsin * p / ratio);
    }
    result
}

#[test]
fn elastic_layer_heat_matches_independent_fields_across_mode_conversion() {
    let incident = fluid(1.2, 340.0);
    let receiver = solid(30.0, 700.0, 350.0);
    let first = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let second = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &first,
            thickness_m: 0.025,
        },
        Layer::FluidGap { thickness_m: 0.012 },
        Layer::Porous {
            law: &second,
            thickness_m: 0.014,
        },
    ];
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 4000.0, &layers).unwrap();
    // Both P/S radiate, then only S radiates, then both are evanescent. The
    // latter load still changes layer interference and dissipation.
    for sine in [0.2_f64, 0.65, 0.985] {
        let cosine = (1.0 - sine * sine).sqrt();
        let response = stack.scattering(cosine).unwrap();
        let report = stack
            .layer_absorption(cosine, 4.5, LayerAbsorptionBudget { max_layers: 3 }, || {
                false
            })
            .unwrap();
        let oracle = independent_layer_heat(&stack, cosine, 4096);
        assert_eq!(report.layers().len(), layers.len());
        assert_eq!(report.reflected_fraction, response.reflected_fraction);
        assert_eq!(
            report.transmitted_fraction,
            response.longitudinal_transmitted_fraction + response.shear_transmitted_fraction
        );
        for (layer, expected) in report.layers().iter().zip(oracle) {
            for (actual, expected) in [
                layer.viscous_absorbed_fraction,
                layer.thermal_absorbed_fraction,
            ]
            .into_iter()
            .zip(expected)
            {
                if expected == 0.0 {
                    assert_eq!(actual, 0.0);
                } else {
                    assert!(actual > 0.0);
                    assert!(
                        (actual / expected - 1.0).abs() < 2e-7,
                        "sin(theta)={sine}, actual={actual}, oracle={expected}"
                    );
                }
            }
            near_real(
                layer.absorbed_power_w,
                4.5 * layer.absorbed_fraction(),
                1e-15,
            );
        }
        near_real(report.absorbed_fraction, response.absorbed_fraction, 3e-14);
        assert!(report.balance_residual_w.abs() < 5e-14);
        if sine > 0.5 {
            assert_eq!(response.longitudinal_transmitted_fraction, 0.0);
        }
        if sine > 0.98 {
            assert_eq!(report.transmitted_fraction, 0.0);
        } else {
            assert!(response.shear_transmitted_fraction > 0.0);
        }
    }
}

#[test]
fn elastic_layer_heat_retains_sub_roundoff_loss_and_normal_fluid_limit() {
    let incident = fluid(1.2, 340.0);
    let receiver = solid(1.2, 340.0, 170.0);
    let receiving_fluid = FluidAcousticState::try_new(1.2, 340.0).unwrap();
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let omega = 4000.0;
    let d = 1e-20;
    let layers = [Layer::Porous {
        law: &foam,
        thickness_m: d,
    }];
    let elastic = ElasticTransmittingLayers::new(&incident, &receiver, omega, &layers).unwrap();
    let fluid_stack = TransmittingLayers::new(&incident, &receiving_fluid, omega, &layers).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 1 };
    let report = elastic
        .layer_absorption(1.0, 7.0, budget, || false)
        .unwrap();
    let original = fluid_stack
        .layer_absorption(1.0, 7.0, budget, || false)
        .unwrap();
    let response = elastic.prepared_layers().material_response(0).unwrap();
    let z0 = incident.acoustic().characteristic_impedance_pa_s_m();
    // Matched normal transmission: p=1 and v=1/z0 at d=0. Integrating
    // the local loss to first order in d gives these independent limits.
    let expected_viscous = omega * d * -response.dynamic_density_kg_m3().im / z0;
    let expected_thermal = omega * d * z0 * -response.dynamic_bulk_modulus_pa().recip().im;
    let row = report.layers()[0];
    assert!(report.absorbed_fraction > 0.0 && report.absorbed_fraction < f64::EPSILON);
    assert!((row.viscous_absorbed_fraction / expected_viscous - 1.0).abs() < 2e-13);
    assert!((row.thermal_absorbed_fraction / expected_thermal - 1.0).abs() < 2e-13);
    assert!((report.absorbed_fraction / original.absorbed_fraction - 1.0).abs() < 2e-13);
    assert!(report.absorbed_power_w > 0.0);
    assert_eq!(
        elastic.scattering(1.0).unwrap().shear_transmitted_fraction,
        0.0
    );
}

#[test]
fn elastic_layer_heat_preserves_subdivision_and_opaque_limits() {
    let incident = fluid(1.2, 340.0);
    let receiver = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let whole = ElasticTransmittingLayers::new(
        &incident,
        &receiver,
        4000.0,
        &[Layer::Porous {
            law: &foam,
            thickness_m: 0.04,
        }],
    )
    .unwrap();
    let split = ElasticTransmittingLayers::new(
        &incident,
        &receiver,
        4000.0,
        &[
            Layer::Porous {
                law: &foam,
                thickness_m: 0.017,
            },
            Layer::FluidGap { thickness_m: 0.0 },
            Layer::Porous {
                law: &foam,
                thickness_m: 0.023,
            },
        ],
    )
    .unwrap();
    for cosine in [1.0, 0.9, 0.5, 0.1] {
        let a = whole
            .layer_absorption(cosine, 3.0, LayerAbsorptionBudget { max_layers: 1 }, || {
                false
            })
            .unwrap();
        let b = split
            .layer_absorption(cosine, 3.0, LayerAbsorptionBudget { max_layers: 3 }, || {
                false
            })
            .unwrap();
        near_real(a.absorbed_fraction, b.absorbed_fraction, 4e-14);
        near_real(a.transmitted_fraction, b.transmitted_fraction, 4e-14);
        near_real(
            a.layers()[0].viscous_absorbed_fraction,
            b.layers().iter().map(|v| v.viscous_absorbed_fraction).sum(),
            4e-14,
        );
        near_real(
            a.layers()[0].thermal_absorbed_fraction,
            b.layers().iter().map(|v| v.thermal_absorbed_fraction).sum(),
            4e-14,
        );
        assert_eq!(b.layers()[1].absorbed_fraction(), 0.0);
    }
    let opaque = ElasticTransmittingLayers::new(
        &incident,
        &receiver,
        4000.0,
        &[Layer::Porous {
            law: &foam,
            thickness_m: 1000.0,
        }],
    )
    .unwrap();
    let report = opaque
        .layer_absorption(0.8, 3.0, LayerAbsorptionBudget { max_layers: 1 }, || false)
        .unwrap();
    assert_eq!(report.transmitted_power_w, 0.0);
    assert!(report.layers()[0].viscous_absorbed_fraction > 0.0);
    assert!(report.layers()[0].thermal_absorbed_fraction > 0.0);
    near_real(
        report.absorbed_fraction + report.reflected_fraction,
        1.0,
        3e-14,
    );
}

#[test]
fn elastic_layer_heat_keeps_lossless_rows_zero_and_refuses_cancelled_work() {
    let incident = fluid(1.2, 340.0);
    let receiver = solid(30.0, 700.0, 350.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let layers = [
        Layer::Porous {
            law: &foam,
            thickness_m: 0.0,
        },
        Layer::FluidGap { thickness_m: 0.04 },
    ];
    let stack = ElasticTransmittingLayers::new(&incident, &receiver, 4000.0, &layers).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 2 };
    let report = stack.layer_absorption(0.8, 4.0, budget, || false).unwrap();
    assert_eq!(report.absorbed_fraction, 0.0);
    assert!(report.layers().iter().all(|v| v.absorbed_fraction() == 0.0));
    for cancel_on in 1..=2 * layers.len() + 2 {
        let mut calls = 0;
        assert_eq!(
            stack.layer_absorption(0.8, 4.0, budget, || {
                calls += 1;
                calls == cancel_on
            }),
            Err(LayerAbsorptionError::Cancelled)
        );
        assert_eq!(
            stack.layer_absorption(0.8, 4.0, budget, || false).unwrap(),
            report
        );
    }
    let zero = stack.layer_absorption(0.8, 0.0, budget, || false).unwrap();
    assert_eq!(zero.reflected_fraction, report.reflected_fraction);
    assert_eq!(zero.transmitted_fraction, report.transmitted_fraction);
    assert_eq!(
        [
            zero.reflected_power_w,
            zero.transmitted_power_w,
            zero.absorbed_power_w
        ],
        [0.0; 3]
    );
    assert!(
        stack
            .layer_absorption(0.8, 4.0, LayerAbsorptionBudget { max_layers: 1 }, || false)
            .is_err()
    );
    assert!(
        stack
            .layer_absorption(0.8, 4.0, LayerAbsorptionBudget { max_layers: 65 }, || false)
            .is_err()
    );
    for cosine in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .layer_absorption(cosine, 4.0, budget, || false)
                .is_err()
        );
    }
    for watts in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .layer_absorption(0.8, watts, budget, || false)
                .is_err()
        );
    }
    let empty = ElasticTransmittingLayers::new(&incident, &receiver, 4000.0, &[]).unwrap();
    let report = empty
        .layer_absorption(0.8, 4.0, LayerAbsorptionBudget { max_layers: 0 }, || false)
        .unwrap();
    assert!(report.layers().is_empty());
    assert_eq!(report.absorbed_power_w, 0.0);
}
