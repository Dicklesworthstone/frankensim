//! Exact finite-plate controls from original fluid and elastic material states.
use fs_blake3::ContentHash;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::layers::{
    ElasticPlateLayers, ElasticTransmittingLayers, Layer, TransmittingLayers,
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
