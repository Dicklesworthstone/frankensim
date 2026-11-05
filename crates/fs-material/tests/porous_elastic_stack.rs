//! G1/G3: original elastic states, coherent cavities and conserved normal work.
use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::layers::{
    AcousticSpectrumError, AngularSpectrumBudget, ElasticPlateLayers, ElasticPlateStack,
    ElasticPlateStackLayer as Section, ElasticTransmittingLayers, IncidentPlaneWave, Layer,
    LayerAbsorptionBudget, LayerAbsorptionError, MAX_LAYERS, TransmittingLayers,
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
    let qg = outgoing((ci / gap[1]).powi(2) - sine * sine);
    let solution = double_plate_gap_solution(
        incident,
        plates,
        (qg.scale(omega / ci), qg.scale(rho / gap[0]), gap[2]),
        receiver,
        omega,
        cosine,
    );
    (solution[0], solution[11])
}

// The same independent twelve equations admit a complex equivalent-fluid gap;
// its returned amplitudes also determine local heat by spatial quadrature.
// Gap data: normal wavenumber, normal admittance times rho_inc*c_inc, thickness.
fn double_plate_gap_solution(
    incident: [f64; 2],
    plates: [[f64; 4]; 2],
    gap: (C64, C64, f64),
    receiver: [f64; 2],
    omega: f64,
    cosine: f64,
) -> [C64; 12] {
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
    let eg = phase(gap.0, gap.2);
    let yg = gap.1;
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
    std::array::from_fn(|i| matrix[i][12])
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

#[test]
fn interleaved_heat_matches_independent_plate_boundary_fields_and_spatial_integrals() {
    // Liquid/solid contrast keeps the independent twelve-amplitude system
    // well conditioned even across P/S critical angles. The separate air/metal
    // controls below exercise thin plates and strongly mismatched impedances.
    let outside = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4180.0, 0.6, 1e-3).unwrap(),
        1.01,
    )
    .unwrap();
    let first = solid(2700.0, 6000.0, 3000.0);
    let second = solid(3600.0, 5100.0, 2700.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let omega = 8000.0;
    let thickness = 0.025;
    let sections = [
        Section::Plate {
            state: &first,
            thickness_m: 0.007,
        },
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: thickness,
        }),
        Section::Plate {
            state: &second,
            thickness_m: 0.013,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, omega, &sections).unwrap();
    let material = foam.effective_fluid(&outside, omega).unwrap();
    let rho = material.dynamic_density_kg_m3();
    let compliance = material.dynamic_bulk_modulus_pa().recip();
    let k = material.characteristic().1;
    let z0 = 1000.0 * 1500.0;
    for cosine in [1.0_f64, 0.99, 0.97, 0.9, 0.8, 0.4] {
        let kt = omega / 1500.0 * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
        let kn = (k * k - C64::from_re(kt * kt)).sqrt();
        let normal_admittance = kn / rho.scale(omega);
        let oracle = double_plate_gap_solution(
            [1000.0, 1500.0],
            [
                [2700.0, 6000.0, 3000.0, 0.007],
                [3600.0, 5100.0, 2700.0, 0.013],
            ],
            (kn, normal_admittance.scale(z0), thickness),
            [1000.0, 1500.0],
            omega,
            cosine,
        );
        let spectrum = stack
            .prepare_angular_spectrum(
                &[IncidentPlaneWave {
                    direction: [((1.0 - cosine) * (1.0 + cosine)).sqrt(), 0.0, cosine],
                    pressure_pa: C64::ONE,
                }],
                AngularSpectrumBudget { max_components: 1 },
                &mut || false,
            )
            .unwrap();
        let local = spectrum.prepare_layer_field(1, &mut || false).unwrap();
        let entrance = spectrum.incident_field([0.0; 3], &mut || false).unwrap();
        near(entrance.pressure_pa, C64::ONE + oracle[0], 3e-9);
        near(
            entrance.particle_velocity_m_s[2],
            (C64::ONE - oracle[0]).scale(cosine / z0),
            3e-9 / z0,
        );
        near(
            spectrum
                .transmitted_field([0.0; 3], &mut || false)
                .unwrap()
                .pressure_pa,
            oracle[11],
            3e-9,
        );
        let mut integrals = [0.0; 2];
        let mut sampled_heat = [0.0; 2];
        let cells = 2048;
        for i in 0..=cells {
            let x = thickness * i as f64 / cells as f64;
            let forward = oracle[5] * phase(kn, x);
            let backward = oracle[6] * phase(kn, thickness - x);
            let pressure = forward + backward;
            let normal_velocity = (forward - backward) * normal_admittance;
            let tangential_velocity = pressure * (C64::from_re(kt / omega) / rho);
            let weight = if i == 0 || i == cells {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integrals[0] += weight * (normal_velocity.norm_sq() + tangential_velocity.norm_sq());
            integrals[1] += weight * pressure.norm_sq();
            let sample = local.field_with_heat([0.0, 0.0, x], &mut || false).unwrap();
            if i % (cells / 2) == 0 {
                near(sample.field.pressure_pa, pressure, 3e-9);
                near(
                    sample.field.particle_velocity_m_s[0],
                    tangential_velocity,
                    3e-9 / z0,
                );
                near(
                    sample.field.particle_velocity_m_s[2],
                    normal_velocity,
                    3e-9 / z0,
                );
            }
            sampled_heat[0] += weight * sample.viscous_heating_w_m3;
            sampled_heat[1] += weight * sample.thermal_heating_w_m3;
        }
        let scale = omega * z0 / cosine * thickness / (3.0 * cells as f64);
        let expected = [
            -rho.im * scale * integrals[0],
            -compliance.im * scale * integrals[1],
        ];
        let scattering = stack.scattering(cosine).unwrap();
        let report = stack
            .layer_absorption(cosine, 2.3, LayerAbsorptionBudget { max_layers: 3 }, || {
                false
            })
            .unwrap();
        assert_eq!(report.layers().len(), 3);
        assert_eq!(report.layers()[0].absorbed_fraction(), 0.0);
        assert_eq!(report.layers()[2].absorbed_fraction(), 0.0);
        let porous = report.layers()[1];
        let flux_scale = 2.0 * z0 / cosine * thickness / (3.0 * cells as f64);
        near_real(
            sampled_heat[0] * flux_scale,
            porous.viscous_absorbed_fraction,
            3e-9 * porous.viscous_absorbed_fraction,
        );
        near_real(
            sampled_heat[1] * flux_scale,
            porous.thermal_absorbed_fraction,
            3e-9 * porous.thermal_absorbed_fraction,
        );
        near_real(
            porous.viscous_absorbed_fraction,
            expected[0],
            3e-9 * expected[0],
        );
        near_real(
            porous.thermal_absorbed_fraction,
            expected[1],
            3e-9 * expected[1],
        );
        near(scattering.pressure_reflection, oracle[0], 2e-9);
        near(scattering.pressure_transmission, oracle[11], 2e-9);
        near_real(
            report.absorbed_fraction,
            1.0 - oracle[0].norm_sq() - oracle[11].norm_sq(),
            3e-10,
        );
        assert_eq!(report.reflected_fraction, scattering.reflected_fraction);
        assert_eq!(report.transmitted_fraction, scattering.transmitted_fraction);
        near_real(
            report.absorbed_power_w,
            2.3 * porous.absorbed_fraction(),
            1e-14,
        );
        assert!(report.balance_residual_w.abs() < 2e-13);
    }
}

#[test]
fn porous_heat_before_between_and_after_plates_retains_the_normal_longitudinal_limit() {
    let outside = fluid(1.2, 343.0);
    let first = solid(7800.0, 5900.0, 3200.0);
    let second = solid(2700.0, 6000.0, 3000.0);
    let first_fluid = FluidAcousticState::try_new(7800.0, 5900.0).unwrap();
    let second_fluid = FluidAcousticState::try_new(2700.0, 6000.0).unwrap();
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let fleece = PorousFluidLaw::try_new(PorousCard::pinta_fleece_foam02(), 1.81e-5).unwrap();
    let a = Layer::Porous {
        law: &foam,
        thickness_m: 0.01,
    };
    let b = Layer::Porous {
        law: &fleece,
        thickness_m: 0.025,
    };
    let c = Layer::Porous {
        law: &foam,
        thickness_m: 0.007,
    };
    let sections = [
        Section::Fluid(a),
        Section::Plate {
            state: &first,
            thickness_m: 0.00015,
        },
        Section::Fluid(b),
        Section::Plate {
            state: &second,
            thickness_m: 0.0002,
        },
        Section::Fluid(c),
    ];
    let reference_layers = [
        a,
        Layer::FluidLayer {
            fluid: &first_fluid,
            thickness_m: 0.00015,
        },
        b,
        Layer::FluidLayer {
            fluid: &second_fluid,
            thickness_m: 0.0002,
        },
        c,
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, 8000.0, &sections).unwrap();
    let reference = TransmittingLayers::new(&outside, &outside, 8000.0, &reference_layers).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 5 };
    let report = stack.layer_absorption(1.0, 1.7, budget, || false).unwrap();
    let expected = reference
        .layer_absorption(1.0, 1.7, budget, || false)
        .unwrap();
    assert_eq!(stack.layer_count(), 5);
    for (i, (actual, expected)) in report.layers().iter().zip(expected.layers()).enumerate() {
        near_real(
            actual.viscous_absorbed_fraction,
            expected.viscous_absorbed_fraction,
            3e-13,
        );
        near_real(
            actual.thermal_absorbed_fraction,
            expected.thermal_absorbed_fraction,
            3e-13,
        );
        near_real(actual.absorbed_power_w, expected.absorbed_power_w, 5e-13);
        if i % 2 == 0 {
            assert!(actual.absorbed_power_w > 0.0);
        } else {
            assert_eq!(actual.absorbed_power_w, 0.0);
        }
    }
    let zero = stack.layer_absorption(1.0, 0.0, budget, || false).unwrap();
    assert_eq!(zero.absorbed_fraction, report.absorbed_fraction);
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert!(zero.layers().iter().all(|row| row.absorbed_power_w == 0.0));
    let waves = [
        IncidentPlaneWave {
            direction: [0.0, 0.0, 1.0],
            pressure_pa: C64::new(0.6, -0.2),
        },
        IncidentPlaneWave {
            direction: [0.0, 0.0, 1.0],
            pressure_pa: C64::new(0.1, 0.3),
        },
    ];
    let budget = AngularSpectrumBudget { max_components: 2 };
    let actual = stack
        .prepare_angular_spectrum(&waves, budget, &mut || false)
        .unwrap();
    let expected = reference
        .prepare_angular_spectrum(&waves, budget, &mut || false)
        .unwrap();
    for (index, thickness) in [(0, 0.01), (2, 0.025), (4, 0.007)] {
        let actual = actual.prepare_layer_field(index, &mut || false).unwrap();
        let expected = expected.prepare_layer_field(index, &mut || false).unwrap();
        for z in [0.0, 0.5 * thickness, thickness] {
            let a = actual
                .field_with_heat([0.0, 0.0, z], &mut || false)
                .unwrap();
            let b = expected
                .field_with_heat([0.0, 0.0, z], &mut || false)
                .unwrap();
            near(a.field.pressure_pa, b.field.pressure_pa, 2e-12);
            near(
                a.field.particle_velocity_m_s[2],
                b.field.particle_velocity_m_s[2],
                2e-14,
            );
            near_real(
                a.viscous_heating_w_m3,
                b.viscous_heating_w_m3,
                2e-10 * b.viscous_heating_w_m3,
            );
            near_real(
                a.thermal_heating_w_m3,
                b.thermal_heating_w_m3,
                2e-10 * b.thermal_heating_w_m3,
            );
        }
    }
}

#[test]
fn interleaved_heat_keeps_weak_loss_and_opaque_plate_shielding() {
    let outside = fluid(1.2, 343.0);
    let matched = solid(1.2, 343.0, 150.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let thickness = 1e-22;
    let sections = [
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: thickness,
        }),
        Section::Plate {
            state: &matched,
            thickness_m: 0.01,
        },
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: 2.0 * thickness,
        }),
    ];
    let omega = 8000.0;
    let stack = ElasticPlateStack::new(&outside, &outside, omega, &sections).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 3 };
    let report = stack.layer_absorption(1.0, 1.0, budget, || false).unwrap();
    let material = foam.effective_fluid(&outside, omega).unwrap();
    let z0 = outside.acoustic().characteristic_impedance_pa_s_m();
    let viscous_per_m = -omega * material.dynamic_density_kg_m3().im / z0;
    let thermal_per_m = -omega * material.dynamic_bulk_modulus_pa().recip().im * z0;
    for (i, d) in [(0, thickness), (2, 2.0 * thickness)] {
        let row = report.layers()[i];
        near_real(
            row.viscous_absorbed_fraction,
            viscous_per_m * d,
            2e-13 * viscous_per_m * d,
        );
        near_real(
            row.thermal_absorbed_fraction,
            thermal_per_m * d,
            2e-13 * thermal_per_m * d,
        );
    }
    assert!(report.absorbed_fraction > 0.0 && report.absorbed_fraction < f64::EPSILON);
    let opaque = solid(7800.0, 5900.0, 3200.0);
    let sections = [
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: 0.01,
        }),
        Section::Plate {
            state: &opaque,
            thickness_m: 100.0,
        },
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: 0.01,
        }),
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, omega, &sections).unwrap();
    let report = stack.layer_absorption(0.7, 1.0, budget, || false).unwrap();
    assert!(report.layers()[0].absorbed_power_w > 0.0);
    assert_eq!(report.layers()[1].absorbed_power_w, 0.0);
    assert_eq!(report.layers()[2].absorbed_power_w, 0.0);
    assert_eq!(report.transmitted_power_w, 0.0);
    let spectrum = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.51_f64.sqrt(), 0.0, 0.7],
                pressure_pa: C64::ONE,
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let front = spectrum.prepare_layer_field(0, &mut || false).unwrap();
    let back = spectrum.prepare_layer_field(2, &mut || false).unwrap();
    let front = front
        .field_with_heat([0.0, 0.0, 0.005], &mut || false)
        .unwrap();
    let back = back
        .field_with_heat([0.0, 0.0, 0.005], &mut || false)
        .unwrap();
    assert!(front.viscous_heating_w_m3 > 0.0 && front.thermal_heating_w_m3 > 0.0);
    assert_eq!(back.viscous_heating_w_m3, 0.0);
    assert_eq!(back.thermal_heating_w_m3, 0.0);
    assert_eq!(back.field.pressure_pa, C64::ZERO);
}

#[test]
fn interleaved_absorption_counts_lossless_sections_and_cancels_both_passes() {
    let outside = fluid(1.2, 343.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let sections = [
        Section::Plate {
            state: &plate,
            thickness_m: 0.002,
        },
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: 0.015,
        }),
        Section::Fluid(Layer::FluidGap { thickness_m: 0.0 }),
        Section::Plate {
            state: &plate,
            thickness_m: 0.001,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, 8000.0, &sections).unwrap();
    let budget = LayerAbsorptionBudget { max_layers: 4 };
    let report = stack.layer_absorption(0.8, 1.0, budget, || false).unwrap();
    assert_eq!(report.layers()[2].absorbed_power_w, 0.0);
    assert!(
        stack
            .layer_absorption(0.8, 1.0, LayerAbsorptionBudget { max_layers: 3 }, || false)
            .is_err()
    );
    for stop in 1..=2 * sections.len() + 2 {
        let mut polls = 0;
        assert!(matches!(
            stack.layer_absorption(0.8, 1.0, budget, || {
                polls += 1;
                polls == stop
            }),
            Err(LayerAbsorptionError::Cancelled)
        ));
    }
    assert_eq!(
        stack.layer_absorption(0.8, 1.0, budget, || false).unwrap(),
        report
    );
    for watts in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            stack
                .layer_absorption(0.8, watts, budget, || false)
                .is_err()
        );
    }
    let empty = ElasticPlateStack::new(&outside, &outside, 8000.0, &[]).unwrap();
    let empty = empty
        .layer_absorption(0.8, 1.0, LayerAbsorptionBudget { max_layers: 0 }, || false)
        .unwrap();
    assert!(empty.layers().is_empty());
    assert_eq!(empty.absorbed_power_w, 0.0);
    assert_eq!(empty.transmitted_power_w, 1.0);
}

#[test]
fn coherent_plate_cavity_heat_matches_complex_boundary_superposition_and_nodes() {
    let outside = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4180.0, 0.6, 1e-3).unwrap(),
        1.01,
    )
    .unwrap();
    let first = solid(2700.0, 6000.0, 3000.0);
    let second = solid(3600.0, 5100.0, 2700.0);
    let foam = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let omega = 8000.0;
    let thickness = 0.025;
    let sections = [
        Section::Plate {
            state: &first,
            thickness_m: 0.007,
        },
        Section::Fluid(Layer::Porous {
            law: &foam,
            thickness_m: thickness,
        }),
        Section::Plate {
            state: &second,
            thickness_m: 0.013,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, omega, &sections).unwrap();
    let material = foam.effective_fluid(&outside, omega).unwrap();
    let rho = material.dynamic_density_kg_m3();
    let compliance = material.dynamic_bulk_modulus_pa().recip();
    let k = material.characteristic().1;
    let waves = [
        (0.99_f64, 0.0_f64, C64::new(0.8, 0.3)),
        (0.9, 0.7, C64::new(-0.4, 0.6)),
        (0.4, 2.0, C64::new(0.5, -0.8)),
    ]
    .map(|(cosine, azimuth, pressure_pa)| {
        let sine = ((1.0 - cosine) * (1.0 + cosine)).sqrt();
        IncidentPlaneWave {
            direction: [sine * azimuth.cos(), sine * azimuth.sin(), cosine],
            pressure_pa,
        }
    });
    let oracle = waves.map(|wave| {
        let tangent = [wave.direction[0], wave.direction[1]].map(|d| omega / 1500.0 * d);
        let kn = (k * k - C64::from_re(tangent[0].powi(2) + tangent[1].powi(2))).sqrt();
        let yn = kn / rho.scale(omega);
        let solution = double_plate_gap_solution(
            [1000.0, 1500.0],
            [
                [2700.0, 6000.0, 3000.0, 0.007],
                [3600.0, 5100.0, 2700.0, 0.013],
            ],
            (kn, yn.scale(1.5e6), thickness),
            [1000.0, 1500.0],
            omega,
            wave.direction[2],
        );
        (wave.pressure_pa, tangent, kn, yn, solution[5], solution[6])
    });
    let spectrum = stack
        .prepare_angular_spectrum(
            &waves,
            AngularSpectrumBudget { max_components: 3 },
            &mut || false,
        )
        .unwrap();
    let field = spectrum.prepare_layer_field(1, &mut || false).unwrap();
    assert!(core::ptr::eq(spectrum.stack(), &stack));
    assert!(core::ptr::eq(field.spectrum(), &spectrum));
    assert_eq!(field.layer_index(), 1);
    assert_eq!(spectrum.component_count(), 3);
    let mut interference_changes_heat = false;
    for position in [
        [0.0, 0.0, 0.0],
        [0.013, -0.024, 0.007],
        [-0.03, 0.01, thickness],
    ] {
        let mut pressure = C64::ZERO;
        let mut velocity = [C64::ZERO; 3];
        let mut incoherent_heat = 0.0;
        for (amplitude, kt, kn, yn, forward, backward) in oracle {
            let lateral = phase(C64::from_re(kt[0] * position[0] + kt[1] * position[1]), 1.0);
            let a = amplitude * lateral * forward * phase(kn, position[2]);
            let b = amplitude * lateral * backward * phase(kn, thickness - position[2]);
            let p = a + b;
            let v = [
                p * (C64::from_re(kt[0] / omega) / rho),
                p * (C64::from_re(kt[1] / omega) / rho),
                (a - b) * yn,
            ];
            pressure = pressure + p;
            for axis in 0..3 {
                velocity[axis] = velocity[axis] + v[axis];
            }
            incoherent_heat += -0.5
                * omega
                * (rho.im * v.iter().map(|v| v.norm_sq()).sum::<f64>()
                    + compliance.im * p.norm_sq());
        }
        let viscous = -0.5 * omega * rho.im * velocity.iter().map(|v| v.norm_sq()).sum::<f64>();
        let thermal = -0.5 * omega * compliance.im * pressure.norm_sq();
        let sample = field.field_with_heat(position, &mut || false).unwrap();
        near(sample.field.pressure_pa, pressure, 5e-9);
        for axis in 0..3 {
            near(
                sample.field.particle_velocity_m_s[axis],
                velocity[axis],
                5e-9 / 1.5e6,
            );
        }
        near_real(sample.viscous_heating_w_m3, viscous, 3e-8 * viscous);
        near_real(sample.thermal_heating_w_m3, thermal, 3e-8 * thermal);
        interference_changes_heat |=
            (viscous + thermal - incoherent_heat).abs() > 0.05 * incoherent_heat;
    }
    assert!(
        interference_changes_heat,
        "coherent fields must not become a power average"
    );

    // Opposite-phase, equal-angle sources cancel p and v_n at x=0 while their
    // transverse velocities reinforce. Identical directions cancel everything.
    let a = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::ONE,
    };
    for opposite_direction in [false, true] {
        let b = IncidentPlaneWave {
            direction: [if opposite_direction { -0.6 } else { 0.6 }, 0.0, 0.8],
            pressure_pa: -C64::ONE,
        };
        let spectrum = stack
            .prepare_angular_spectrum(
                &[a, b],
                AngularSpectrumBudget { max_components: 2 },
                &mut || false,
            )
            .unwrap();
        let local = spectrum.prepare_layer_field(1, &mut || false).unwrap();
        let sample = local
            .field_with_heat([0.0, 0.0, 0.01], &mut || false)
            .unwrap();
        assert_eq!(sample.field.pressure_pa, C64::ZERO);
        assert_eq!(sample.field.particle_velocity_m_s[2], C64::ZERO);
        assert_eq!(sample.thermal_heating_w_m3, 0.0);
        if opposite_direction {
            assert!(sample.viscous_heating_w_m3 > 0.0);
        } else {
            assert_eq!(sample.viscous_heating_w_m3, 0.0);
        }
    }
}

#[test]
fn elastic_spectrum_observes_lossless_fluid_identities_and_refuses_atomically() {
    let outside = fluid(1.2, 343.0);
    let plate = solid(2700.0, 6000.0, 3000.0);
    let sections = [
        Section::Plate {
            state: &plate,
            thickness_m: 0.002,
        },
        Section::Fluid(Layer::FluidGap { thickness_m: 0.015 }),
        Section::Fluid(Layer::FluidGap { thickness_m: 0.0 }),
        Section::Plate {
            state: &plate,
            thickness_m: 0.001,
        },
    ];
    let stack = ElasticPlateStack::new(&outside, &outside, 8000.0, &sections).unwrap();
    let waves = [IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::ONE,
    }];
    let budget = AngularSpectrumBudget { max_components: 1 };
    assert!(
        stack
            .prepare_angular_spectrum(
                &waves,
                AngularSpectrumBudget { max_components: 0 },
                &mut || false
            )
            .is_err()
    );
    let mut preparation_polls = 0;
    let spectrum = stack
        .prepare_angular_spectrum(&waves, budget, &mut || {
            preparation_polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=preparation_polls {
        let mut polls = 0;
        assert!(matches!(
            stack.prepare_angular_spectrum(&waves, budget, &mut || {
                polls += 1;
                polls == stop
            }),
            Err(AcousticSpectrumError::Cancelled)
        ));
    }
    assert!(spectrum.prepare_layer_field(0, &mut || false).is_err());
    assert!(
        spectrum
            .prepare_layer_field(sections.len(), &mut || false)
            .is_err()
    );
    let mut preparation_polls = 0;
    let local = spectrum
        .prepare_layer_field(1, &mut || {
            preparation_polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=preparation_polls {
        let mut polls = 0;
        assert!(matches!(
            spectrum.prepare_layer_field(1, &mut || {
                polls += 1;
                polls == stop
            }),
            Err(AcousticSpectrumError::Cancelled)
        ));
    }
    let mut observation_polls = 0;
    let sample = local
        .field_with_heat([0.0, 0.0, 0.007], &mut || {
            observation_polls += 1;
            false
        })
        .unwrap();
    assert_eq!(sample.viscous_heating_w_m3, 0.0);
    assert_eq!(sample.thermal_heating_w_m3, 0.0);
    for stop in 1..=observation_polls {
        let mut polls = 0;
        assert!(matches!(
            local.field_with_heat([0.0, 0.0, 0.007], &mut || {
                polls += 1;
                polls == stop
            }),
            Err(AcousticSpectrumError::Cancelled)
        ));
    }
    assert_eq!(
        local
            .field_with_heat([0.0, 0.0, 0.007], &mut || false)
            .unwrap(),
        sample
    );
    let identity = spectrum.prepare_layer_field(2, &mut || false).unwrap();
    near(
        identity.field([0.0; 3], &mut || false).unwrap().pressure_pa,
        local
            .field([0.0, 0.0, 0.015], &mut || false)
            .unwrap()
            .pressure_pa,
        1e-12,
    );
    assert!(identity.field_with_heat([0.0; 3], &mut || false).is_err());
    assert!(identity.field([0.0, 0.0, 1e-6], &mut || false).is_err());
}
