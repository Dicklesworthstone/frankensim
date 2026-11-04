//! G1/G3: original elastic states, coherent cavities and conserved normal work.
use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::layers::{
    ElasticPlateLayers, ElasticPlateStack, ElasticPlateStackLayer as Section,
    ElasticTransmittingLayers, Layer, MAX_LAYERS, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;

fn fluid(rho: f64, c: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(rho, c).unwrap(),
        FluidTransportState::try_new(rho, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}

fn solid(rho: f64, cp: f64, cs: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 * cs * cs / 3.0),
        rho * cs * cs,
        ContentHash([29; 32]),
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
        "{actual:e} != {expected:e}"
    );
}

fn outgoing(square: f64) -> C64 {
    if square >= 0.0 {
        C64::from_re(square.sqrt())
    } else {
        C64::new(0.0, -(-square).sqrt())
    }
}

fn phase(kn: C64, thickness: f64) -> C64 {
    C64::new((-kn.re * thickness).cos(), (-kn.re * thickness).sin())
        .scale((kn.im * thickness).exp())
}

// Independent global wave-amplitude solve. Twelve unknowns are entrance R,
// four P/S amplitudes in each solid, two fluid-gap amplitudes, and exit T.
// The 12 original continuity/shear-free equations do not use the production
// parity reduction, transfer pair, determinant factorization or composition.
fn double_plate_oracle(
    incident: [f64; 2],
    plates: [[f64; 4]; 2],
    gap: [f64; 3],
    receiver: [f64; 2],
    omega: f64,
    cosine: f64,
) -> (C64, C64) {
    let [rho, ci] = incident;
    let sine = ((1.0 - cosine) * (1.0 + cosine)).sqrt();
    let face = |[rs, cp, cs, thickness]: [f64; 4]| {
        let p = outgoing(1.0 - (cp * sine / ci).powi(2));
        let s = outgoing(1.0 - (cs * sine / ci).powi(2));
        let ss = cs * sine / ci;
        let a = 1.0 - 2.0 * ss * ss;
        let zp = rs * cp / (rho * ci);
        let zs = rs * cs / (rho * ci);
        let ep = phase(p.scale(omega / cp), thickness);
        let es = phase(s.scale(omega / cs), thickness);
        let pressure = [C64::from_re(zp * a), s.scale(-2.0 * zs * ss)];
        let velocity = [p, C64::from_re(-ss)];
        let shear = [p.scale(2.0 * ss), C64::from_re(a)];
        let front = [
            [pressure[0], pressure[1], pressure[0] * ep, pressure[1] * es],
            [
                velocity[0],
                velocity[1],
                -velocity[0] * ep,
                -velocity[1] * es,
            ],
            [shear[0], shear[1], -shear[0] * ep, -shear[1] * es],
        ];
        let back = [
            [pressure[0] * ep, pressure[1] * es, pressure[0], pressure[1]],
            [
                velocity[0] * ep,
                velocity[1] * es,
                -velocity[0],
                -velocity[1],
            ],
            [shear[0] * ep, shear[1] * es, -shear[0], -shear[1]],
        ];
        (front, back)
    };
    let (front_a, back_a) = face(plates[0]);
    let (front_b, back_b) = face(plates[1]);
    let qg = outgoing((ci / gap[1]).powi(2) - sine * sine);
    let eg = phase(qg.scale(omega / ci), gap[2]);
    let yg = qg.scale(rho / gap[0]);
    let yr = outgoing((ci / receiver[1]).powi(2) - sine * sine).scale(rho / receiver[0]);
    let mut matrix = [[C64::ZERO; 13]; 12];
    for column in 0..4 {
        for quantity in 0..3 {
            matrix[quantity][1 + column] = if quantity == 2 {
                front_a[quantity][column]
            } else {
                -front_a[quantity][column]
            };
            matrix[3 + quantity][1 + column] = back_a[quantity][column];
            matrix[6 + quantity][7 + column] = if quantity == 2 {
                front_b[quantity][column]
            } else {
                -front_b[quantity][column]
            };
            matrix[9 + quantity][7 + column] = back_b[quantity][column];
        }
    }
    matrix[0][0] = C64::ONE;
    matrix[0][12] = -C64::ONE;
    matrix[1][0] = C64::from_re(-cosine);
    matrix[1][12] = C64::from_re(-cosine);
    matrix[3][5] = -C64::ONE;
    matrix[3][6] = -eg;
    matrix[4][5] = -yg;
    matrix[4][6] = yg * eg;
    matrix[6][5] = eg;
    matrix[6][6] = C64::ONE;
    matrix[7][5] = yg * eg;
    matrix[7][6] = -yg;
    matrix[9][11] = -C64::ONE;
    matrix[10][11] = -yr;
    for column in 0..12 {
        let pivot = (column..12)
            .max_by(|&i, &j| matrix[i][column].abs().total_cmp(&matrix[j][column].abs()))
            .unwrap();
        matrix.swap(column, pivot);
        let divisor = matrix[column][column];
        assert!(divisor.abs() > 0.0);
        for entry in column..13 {
            matrix[column][entry] = matrix[column][entry] / divisor;
        }
        for row in 0..12 {
            if row == column {
                continue;
            }
            let multiplier = matrix[row][column];
            for entry in column..13 {
                matrix[row][entry] = matrix[row][entry] - multiplier * matrix[column][entry];
            }
        }
    }
    (matrix[0][12], matrix[11][12])
}

#[test]
fn single_plate_and_zero_plate_retain_the_existing_owners() {
    let outside = fluid(1000.0, 1500.0);
    let receiver = fluid(800.0, 1800.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let gap = Layer::FluidGap { thickness_m: 0.031 };
    for thickness in [0.0, 1e-8, 0.007, 0.2, 1.4] {
        let layers = [
            Section::Fluid(gap),
            Section::Plate {
                state: &plate,
                thickness_m: thickness,
            },
        ];
        let stack = ElasticPlateStack::new(&outside, &receiver, 16000.0, &layers).unwrap();
        let original =
            ElasticPlateLayers::new(&outside, &plate, &receiver, 16000.0, thickness, &[gap])
                .unwrap();
        assert_eq!(stack.layers().len(), 2);
        match stack.layers().nth(1).unwrap() {
            Section::Plate { state, .. } => assert!(core::ptr::eq(state, &plate)),
            _ => panic!("original plate was lost"),
        }
        for cosine in [1.0, 0.99, 0.94, 0.8, 0.55, 0.03] {
            let a = stack.scattering(cosine).unwrap();
            let b = original.scattering(cosine).unwrap();
            near(a.pressure_reflection, b.pressure_reflection, 5e-13);
            near(a.pressure_transmission, b.pressure_transmission, 5e-13);
            near_real(a.transmitted_fraction, b.transmitted_fraction, 5e-13);
            assert_eq!(a.absorbed_fraction, 0.0);
        }
    }
    for gaps in [vec![], vec![gap]] {
        let sections: Vec<_> = gaps.iter().copied().map(Section::Fluid).collect();
        let stack = ElasticPlateStack::new(&outside, &receiver, 16000.0, &sections).unwrap();
        let original = TransmittingLayers::new(&outside, &receiver, 16000.0, &gaps).unwrap();
        for cosine in [1.0, 0.8, 0.3] {
            let a = stack.scattering(cosine).unwrap();
            let b = original.scattering(cosine).unwrap();
            assert_eq!(a.pressure_reflection, b.pressure_reflection);
            assert_eq!(a.pressure_transmission, b.pressure_transmission);
            assert_eq!(a.transmitted_fraction, b.transmitted_fraction);
        }
    }
}

#[test]
fn double_plate_solves_twelve_original_boundary_conditions() {
    let outside = fluid(1000.0, 1500.0);
    let second = solid(3600.0, 5100.0, 2700.0);
    for (rho, cp, cs) in [(2700.0, 6000.0, 3000.0), (2200.0, 4000.0, 3200.0)] {
        let first = solid(rho, cp, cs);
        for gap_speed in [1100.0, 2400.0] {
            let gap = fluid(900.0, gap_speed);
            for receiver_speed in [1200.0, 1800.0] {
                let receiver = fluid(800.0, receiver_speed);
                for thickness in [0.007, 0.09] {
                    for gap_thickness in [0.01, 0.07] {
                        let layers = [
                            Section::Plate {
                                state: &first,
                                thickness_m: thickness,
                            },
                            Section::Fluid(Layer::FluidLayer {
                                fluid: &gap,
                                thickness_m: gap_thickness,
                            }),
                            Section::Plate {
                                state: &second,
                                thickness_m: 0.013,
                            },
                        ];
                        let stack =
                            ElasticPlateStack::new(&outside, &receiver, 17000.0, &layers).unwrap();
                        for sine in [0.05_f64, 0.3, 0.6, 0.85] {
                            let cosine = (1.0 - sine * sine).sqrt();
                            let actual = stack.scattering(cosine).unwrap();
                            let (r, t) = double_plate_oracle(
                                [1000.0, 1500.0],
                                [[rho, cp, cs, thickness], [3600.0, 5100.0, 2700.0, 0.013]],
                                [900.0, gap_speed, gap_thickness],
                                [800.0, receiver_speed],
                                17000.0,
                                cosine,
                            );
                            near(actual.pressure_reflection, r, 4e-10);
                            near(actual.pressure_transmission, t, 4e-10);
                            near_real(
                                actual.reflected_fraction + actual.transmitted_fraction,
                                1.0,
                                3e-13,
                            );
                            assert_eq!(actual.absorbed_fraction, 0.0);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn fluid_cavity_resonance_transmits_through_two_reflecting_evanescent_plates() {
    let outside = fluid(1.2, 343.0);
    let plate = solid(7800.0, 5900.0, 3200.0);
    let omega = 118845.86838915091;
    // From the independent single-plate boundary phase: delta=arg(R)+pi.
    // Each individual plate transmits <0.006, but coherent cavity resonance
    // transmits almost all power. Multiplying power transmissions would fail.
    let gap = 0.015474056334315662;
    for ratio in [0.99, 1.0, 1.01] {
        let layers = [
            Section::Plate {
                state: &plate,
                thickness_m: 0.001,
            },
            Section::Fluid(Layer::FluidGap {
                thickness_m: gap * ratio,
            }),
            Section::Plate {
                state: &plate,
                thickness_m: 0.001,
            },
        ];
        let actual = ElasticPlateStack::new(&outside, &outside, omega, &layers)
            .unwrap()
            .scattering(0.6)
            .unwrap();
        let (r, t) = double_plate_oracle(
            [1.2, 343.0],
            [[7800.0, 5900.0, 3200.0, 0.001]; 2],
            [1.2, 343.0, gap * ratio],
            [1.2, 343.0],
            omega,
            0.6,
        );
        near(actual.pressure_reflection, r, 2e-7);
        near(actual.pressure_transmission, t, 2e-7);
        near_real(
            actual.reflected_fraction + actual.transmitted_fraction,
            1.0,
            8e-13,
        );
        if ratio == 1.0 {
            assert!(actual.transmitted_fraction > 0.99999);
        } else {
            assert!(actual.transmitted_fraction < 0.02);
        }
    }
    let single = ElasticPlateLayers::new(&outside, &plate, &outside, omega, 0.001, &[])
        .unwrap()
        .scattering(0.6)
        .unwrap();
    assert!(single.transmitted_fraction < 0.006);
    assert!(single.longitudinal_normal_wavenumber_rad_m.im < 0.0);
    assert!(single.shear_normal_wavenumber_rad_m.im < 0.0);
}

#[test]
fn reversing_unequal_ports_and_plate_order_preserves_reciprocal_power_and_phase() {
    let front = fluid(1000.0, 1500.0);
    let back = fluid(800.0, 1200.0);
    let gap = fluid(900.0, 1600.0);
    let a = solid(2700.0, 6000.0, 3000.0);
    let b = solid(2200.0, 4000.0, 3200.0);
    let layers = [
        Section::Plate {
            state: &a,
            thickness_m: 0.024,
        },
        Section::Fluid(Layer::FluidLayer {
            fluid: &gap,
            thickness_m: 0.041,
        }),
        Section::Plate {
            state: &b,
            thickness_m: 0.07,
        },
    ];
    let reverse = [layers[2], layers[1], layers[0]];
    let cosine = 0.88_f64;
    let reverse_cosine = (1.0 - (1200.0_f64 / 1500.0).powi(2) * (1.0 - cosine * cosine)).sqrt();
    let f = ElasticPlateStack::new(&front, &back, 14000.0, &layers)
        .unwrap()
        .scattering(cosine)
        .unwrap();
    let r = ElasticPlateStack::new(&back, &front, 14000.0, &reverse)
        .unwrap()
        .scattering(reverse_cosine)
        .unwrap();
    let load = (1000.0 * 1500.0) / (800.0 * 1200.0) * reverse_cosine / cosine;
    near_real(f.transmitted_fraction, r.transmitted_fraction, 2e-12);
    near(
        f.pressure_transmission.scale(load),
        r.pressure_transmission,
        2e-12,
    );
}

#[test]
fn critical_gaps_and_thick_plates_keep_tunneling_distinct_from_heat() {
    let outside = fluid(1000.0, 1500.0);
    let fast_gap = fluid(700.0, 2500.0);
    let first = solid(2700.0, 6000.0, 3000.0);
    let second = solid(2200.0, 4000.0, 3200.0);
    let layers = [
        Section::Plate {
            state: &first,
            thickness_m: 0.007,
        },
        Section::Fluid(Layer::FluidLayer {
            fluid: &fast_gap,
            thickness_m: 0.01,
        }),
        Section::Plate {
            state: &second,
            thickness_m: 0.013,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, 16000.0, &layers).unwrap();
    for critical_sine in [1500.0_f64 / 6000.0, 1500.0 / 3000.0, 1500.0 / 2500.0] {
        for factor in [1.0 - 1e-7, 1.0, 1.0 + 1e-7] {
            let cosine = (1.0 - (critical_sine * factor).powi(2)).sqrt();
            let result = stack.scattering(cosine).unwrap();
            near_real(
                result.reflected_fraction + result.transmitted_fraction,
                1.0,
                7e-13,
            );
            assert_eq!(result.absorbed_fraction, 0.0);
            assert!(result.transmitted_fraction > 0.0);
        }
    }
    let thick = [
        Section::Plate {
            state: &first,
            thickness_m: 100.0,
        },
        layers[1],
        layers[2],
    ];
    let opaque = ElasticPlateStack::new(&outside, &outside, 16000.0, &thick)
        .unwrap()
        .scattering(0.7)
        .unwrap();
    let half = ElasticTransmittingLayers::new(&outside, &first, 16000.0, &[])
        .unwrap()
        .scattering(0.7)
        .unwrap();
    assert_eq!(opaque.transmitted_fraction, 0.0);
    assert_eq!(opaque.absorbed_fraction, 0.0);
    near(opaque.pressure_reflection, half.pressure_reflection, 5e-13);
}

#[test]
fn porous_cavity_loss_uses_the_original_material_and_transfer() {
    let outside = fluid(1.2, 343.0);
    let a = solid(7800.0, 5900.0, 3200.0);
    let b = solid(2700.0, 6000.0, 3000.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let porous = Layer::Porous {
        law: &foam,
        thickness_m: 0.025,
    };
    let layers = [
        Section::Plate {
            state: &a,
            thickness_m: 0.00015,
        },
        Section::Fluid(porous),
        Section::Plate {
            state: &b,
            thickness_m: 0.0002,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, 8000.0, &layers).unwrap();
    for cosine in [1.0, 0.99, 0.8, 0.4] {
        let result = stack.scattering(cosine).unwrap();
        assert!(result.absorbed_fraction > 0.0 && result.transmitted_fraction > 0.0);
        near_real(
            result.reflected_fraction + result.transmitted_fraction + result.absorbed_fraction,
            1.0,
            2e-13,
        );
    }
    let p_a = FluidAcousticState::try_new(7800.0, 5900.0).unwrap();
    let p_b = FluidAcousticState::try_new(2700.0, 6000.0).unwrap();
    let reference = TransmittingLayers::new(
        &outside,
        &outside,
        8000.0,
        &[
            Layer::FluidLayer {
                fluid: &p_a,
                thickness_m: 0.00015,
            },
            porous,
            Layer::FluidLayer {
                fluid: &p_b,
                thickness_m: 0.0002,
            },
        ],
    )
    .unwrap()
    .scattering(1.0)
    .unwrap();
    let actual = stack.scattering(1.0).unwrap();
    near(
        actual.pressure_reflection,
        reference.pressure_reflection,
        3e-14,
    );
    near(
        actual.pressure_transmission,
        reference.pressure_transmission,
        3e-14,
    );
    near_real(actual.absorbed_fraction, reference.absorbed_fraction, 3e-14);
}

#[test]
fn maximum_stack_retains_exact_normal_delay_and_distinct_fluid_snapshots() {
    let outside = fluid(1.2, 343.0);
    let plate = solid(1.2, 343.0, 150.0);
    let layers: Vec<_> = (0..MAX_LAYERS)
        .map(|index| {
            if index % 2 == 0 {
                Section::Plate {
                    state: &plate,
                    thickness_m: 0.002,
                }
            } else {
                Section::Fluid(Layer::FluidGap {
                    thickness_m: 0.0013,
                })
            }
        })
        .collect();
    let result = ElasticPlateStack::new(&outside, &outside, 4000.0, &layers)
        .unwrap()
        .scattering(1.0)
        .unwrap();
    near(result.pressure_reflection, C64::ZERO, 3e-14);
    near(
        result.pressure_transmission,
        phase(C64::from_re(4000.0 / 343.0), 32.0 * 0.0033),
        3e-14,
    );
    near_real(result.transmitted_fraction, 1.0, 3e-14);
    struct CountingFluid {
        calls: std::cell::Cell<usize>,
        state: FluidAcousticState,
    }
    impl FluidAcoustics for CountingFluid {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }
    let other = CountingFluid {
        calls: std::cell::Cell::new(0),
        state: FluidAcousticState::try_new(1.5, 320.0).unwrap(),
    };
    let shared = [
        Section::Fluid(Layer::FluidLayer {
            fluid: &other,
            thickness_m: 0.01,
        }),
        Section::Plate {
            state: &plate,
            thickness_m: 0.002,
        },
        Section::Fluid(Layer::FluidLayer {
            fluid: &other,
            thickness_m: 0.02,
        }),
    ];
    let stack = ElasticPlateStack::new(&outside, &other, 4000.0, &shared).unwrap();
    assert_eq!(other.calls.get(), 1);
    for cosine in [1.0, 0.9, 0.7] {
        stack.scattering(cosine).unwrap();
    }
    assert_eq!(other.calls.get(), 1);
    assert!(core::ptr::eq(
        stack.receiving_source(),
        &other as &dyn FluidAcoustics
    ));
}

#[test]
fn unseparated_solids_and_unsupported_sources_are_refused() {
    let outside = fluid(1000.0, 1500.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let p = Section::Plate {
        state: &plate,
        thickness_m: 0.001,
    };
    for layers in [
        vec![p, p],
        vec![p, Section::Fluid(Layer::FluidGap { thickness_m: 0.0 }), p],
    ] {
        assert!(ElasticPlateStack::new(&outside, &outside, 8000.0, &layers).is_err());
    }
    let generic = ElasticAcousticState::try_new_mandel(
        *plate.stiffness_mandel_pa(),
        plate.density_kg_m3(),
        ContentHash([30; 32]),
        0.001,
    )
    .unwrap();
    assert!(
        ElasticPlateStack::new(
            &outside,
            &outside,
            8000.0,
            &[Section::Plate {
                state: &generic,
                thickness_m: 0.001
            }]
        )
        .is_err()
    );
    for thickness in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ElasticPlateStack::new(
                &outside,
                &outside,
                8000.0,
                &[Section::Plate {
                    state: &plate,
                    thickness_m: thickness
                }]
            )
            .is_err()
        );
    }
    assert!(
        ElasticPlateStack::new(
            &outside,
            &outside,
            8000.0,
            &vec![Section::Fluid(Layer::FluidGap { thickness_m: 0.01 }); MAX_LAYERS + 1]
        )
        .is_err()
    );
    let stack = ElasticPlateStack::new(&outside, &outside, 8000.0, &[p]).unwrap();
    for cosine in [0.0, -0.1, 1.001, f64::NAN] {
        assert!(stack.scattering(cosine).is_err());
    }
}
