//! G1/G3/G4: original caloric dispersion, stress radiation and physical records.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
};
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::{
    CoherenceBudget, CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactSourceObserver,
};
use fs_aeroac::material::MaterialAeroacError;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::{PI, TAU};

const TAU_R: f64 = 0.00023;

fn medium(relaxing: bool) -> RelaxationState {
    medium_with_rate(relaxing, TAU_R)
}

fn medium_with_rate(relaxing: bool, tau: f64) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("manufactured-relaxing-gas").unwrap();
    let processes = if relaxing {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * tau,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }]
    } else {
        vec![]
    };
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap()
}

fn tensor() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.7, 0.1),
            C64::new(-0.3, 0.2),
            C64::new(0.9, -0.4),
            C64::new(0.23, -0.15),
            C64::new(-0.11, 0.3),
            C64::new(0.4, 0.2),
        ],
    }
}

fn owner(gas: &RelaxationState, f: f64, x: [f64; 3], y: &[[f64; 3]]) -> PreparedQuadrupoleObserver {
    PreparedQuadrupoleObserver::new_molecular(gas.clone(), f, x, y, y.len()).unwrap()
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-18),
        "{actual:?} != {expected:?}"
    );
}

// Eliminate the manufactured caloric reservoir independently in real arithmetic.
// K = p (cp0 - i omega tau cp_inf)/(cv0 - i omega tau cv_inf).
fn wave(gas: &RelaxationState, omega: f64) -> C64 {
    let s = gas.gas();
    let cv0 = s.specific_heat_cp - s.specific_gas_constant;
    let cvi = cv0 - s.specific_gas_constant;
    let cpi = cvi + s.specific_gas_constant;
    let y = omega * TAU_R;
    let den = cv0 * cv0 + y * y * cvi * cvi;
    let re = s.pressure * (s.specific_heat_cp * cv0 + y * y * cpi * cvi) / den;
    let im = s.pressure * y * (s.specific_heat_cp * cvi - cpi * cv0) / den;
    let factor = s.density * omega * omega / (re * re + im * im);
    let k2re = re * factor;
    let k2im = -im * factor;
    let real = ((k2re.hypot(k2im) + k2re) * 0.5).sqrt();
    C64::new(real, k2im / (2.0 * real))
}

fn green(k: C64, r: f64) -> C64 {
    let (sin, cos) = (k.re * r).sin_cos();
    C64::new(cos, sin).scale((-k.im * r).exp() / (4.0 * PI * r))
}

// Independent radial differentiation and full tensor contraction.
fn pressure(k: C64, x: [f64; 3], source: QuadrupolePhasor) -> C64 {
    let r = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let n = x.map(|v| v / r);
    let q = source.stress_moment_n_m;
    let radial = C64::new(-k.im - 1.0 / r, k.re);
    let g = green(k, r);
    let gp = g * radial;
    let gpp = g * (radial * radial + C64::from_re(1.0 / (r * r)));
    let projection = q[0].scale(n[0] * n[0])
        + q[1].scale(n[1] * n[1])
        + q[2].scale(n[2] * n[2])
        + q[3].scale(2.0 * n[0] * n[1])
        + q[4].scale(2.0 * n[0] * n[2])
        + q[5].scale(2.0 * n[1] * n[2]);
    (gpp - gp.scale(1.0 / r)) * projection + gp.scale(1.0 / r) * (q[0] + q[1] + q[2])
}

#[test]
fn g1_material_dispersion_reaches_the_full_green_hessian_and_vector_velocity() {
    let gas = medium(true);
    for f in [20.0, 800.0, 12000.0] {
        let omega = TAU * f;
        let k = wave(&gas, omega);
        for x in [[0.021, -0.013, 0.017], [0.6, 0.8, -0.3]] {
            let observer = owner(&gas, f, x, &[[0.0; 3]]);
            close(observer.wavenumber_per_m(), k, 2e-14);
            assert_eq!(observer.molecular_state(), Some(&gas));
            let field = observer.evaluate(&[tensor()]).unwrap();
            close(field.pressure_pa, pressure(k, x, tensor()), 3e-12);
            let r = x.iter().map(|v| v * v).sum::<f64>().sqrt();
            let scale = r.min(1.0 / k.abs());
            for axis in 0..3 {
                let errors = [0.002, 0.001, 0.0005].map(|step| {
                    let h = step * scale;
                    let mut plus = x;
                    let mut minus = x;
                    plus[axis] += h;
                    minus[axis] -= h;
                    let gradient =
                        (pressure(k, plus, tensor()) - pressure(k, minus, tensor())).scale(0.5 / h);
                    let velocity = C64::new(gradient.im, -gradient.re)
                        .scale(1.0 / (gas.gas().density * omega));
                    (velocity - field.particle_velocity_m_s[axis]).abs()
                });
                assert!(errors[1] < 0.27 * errors[0], "{f} {axis} {errors:?}");
                assert!(errors[2] < 0.3 * errors[1], "{f} {axis} {errors:?}");
                assert!(errors[2] < 2e-6 * field.particle_velocity_m_s[axis].abs());
            }
        }
    }
}

#[test]
fn g1_isotropic_stress_has_material_phase_and_exact_caloric_power_loss() {
    let gas = medium(true);
    let f = 800.0;
    let k = wave(&gas, TAU * f);
    let q = C64::new(0.003, -0.002);
    let source = QuadrupolePhasor {
        stress_moment_n_m: [q, q, q, C64::ZERO, C64::ZERO, C64::ZERO],
    };
    let radii = [0.12, 0.45, 1.3];
    let mut powers = [0.0; 3];
    for (index, r) in radii.into_iter().enumerate() {
        let observer = owner(&gas, f, [r, 0.0, 0.0], &[[0.0; 3]]);
        let field = observer.evaluate(&[source]).unwrap();
        let expected_p = -(k * k) * q * green(k, r);
        close(field.pressure_pa, expected_p, 2e-13);
        close(
            field.particle_velocity_m_s[0],
            expected_p * (k + C64::new(0.0, 1.0 / r)).scale(1.0 / (gas.gas().density * TAU * f)),
            2e-13,
        );
        assert_eq!(observer.source_delay_s(0), Some(r / gas.gas().sound_speed));
        close(
            C64::from_re(observer.source_phase_delay_s(0).unwrap()),
            C64::from_re(r * k.re / (TAU * f)),
            2e-14,
        );
        assert!(observer.source_phase_delay_s(0).unwrap() < observer.source_delay_s(0).unwrap());
        powers[index] = 4.0 * PI * r * r * field.active_intensity_w_m2[0];
    }
    let k2 = k * k;
    for i in 0..2 {
        // Integral over the shell of omega Im(1/K) |p|²/2. The trace's
        // p=-k² q G supplies its independent source amplitude.
        let radial_integral =
            ((-2.0 * k.im * radii[i]).exp() - (-2.0 * k.im * radii[i + 1]).exp()) / (2.0 * k.im);
        let dissipation = k2.im / (2.0 * gas.gas().density * TAU * f) * (k2 * q).norm_sq()
            / (4.0 * PI)
            * radial_integral;
        close(
            C64::from_re(powers[i] - powers[i + 1]),
            C64::from_re(dissipation),
            2e-12,
        );
        assert!(powers[i] > powers[i + 1] && powers[i + 1] > 0.0);
    }
}

#[test]
fn g3_empty_reservoirs_recover_existing_lossless_tensor_fields() {
    let gas = medium(false);
    for f in [20.0, 800.0, 12000.0] {
        let x = [0.6, -0.8, 0.3];
        let molecular = owner(&gas, f, x, &[[0.0; 3]]);
        let lossless = PreparedQuadrupoleObserver::new(
            gas.gas().acoustic_state().unwrap(),
            f,
            x,
            &[[0.0; 3]],
            1,
        )
        .unwrap();
        assert_eq!(molecular.wavenumber_per_m().im, 0.0);
        assert_eq!(lossless.source_delay_s(0), lossless.source_phase_delay_s(0));
        let actual = molecular.evaluate(&[tensor()]).unwrap();
        let expected = lossless.evaluate(&[tensor()]).unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 3e-12);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                3e-12,
            );
        }
    }
}

#[test]
fn g1_material_wavenumber_survives_tensor_images_and_wedge_boundaries() {
    let gas = medium(true);
    let frequency = 800.0;
    let positions = [[0.2, 0.3, 0.1]];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundary).unwrap();
        let x = [0.5, 0.0, 0.4];
        let direct = owner(&gas, frequency, x, &positions);
        let on_wall = direct.clone().with_orthogonal_planes(&[plane], 2).unwrap();
        for component in 0..6 {
            let mut q = QuadrupolePhasor::default();
            q.stress_moment_n_m[component] = C64::ONE;
            let original = direct.evaluate(&[q]).unwrap();
            let field = on_wall.evaluate(&[q]).unwrap();
            if boundary == IdealPlaneBoundary::Rigid {
                assert!(
                    field.particle_velocity_m_s[1].abs()
                        <= 1e-12
                            * original
                                .particle_velocity_m_s
                                .iter()
                                .map(|v| v.abs())
                                .fold(0.0, f64::max)
                );
            } else {
                assert!(field.pressure_pa.abs() <= 1e-12 * original.pressure_pa.abs());
            }
        }
        assert_eq!(on_wall.molecular_state(), Some(&gas));
        assert!(on_wall.scope().contains("relaxing gas"));
        // An off-wall independent source-image sum catches a lossless image
        // accidentally folded into an otherwise molecular direct path.
        let x = [0.5, 0.4, 0.4];
        let q = tensor();
        let mut image_q = q;
        image_q.stress_moment_n_m[3] = -image_q.stress_moment_n_m[3];
        image_q.stress_moment_n_m[5] = -image_q.stress_moment_n_m[5];
        if boundary == IdealPlaneBoundary::PressureRelease {
            image_q.stress_moment_n_m = image_q.stress_moment_n_m.map(|v| -v);
        }
        let image = [positions[0][0], -positions[0][1], positions[0][2]];
        let expected = owner(&gas, frequency, x, &[positions[0], image])
            .evaluate(&[q, image_q])
            .unwrap();
        let actual = owner(&gas, frequency, x, &positions)
            .with_orthogonal_planes(&[plane], 2)
            .unwrap()
            .evaluate(&[q])
            .unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 2e-12);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                2e-12,
            );
        }
    }
    let angle = PI / 3.0;
    let walls = [
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new(
            [0.0; 3],
            [angle.sin(), -angle.cos(), 0.0],
            IdealPlaneBoundary::Rigid,
        )
        .unwrap(),
    ];
    let source_position = [[0.3, 0.1, 0.0]];
    for x in [[0.6, 0.0, 0.2], [0.6 * angle.cos(), 0.6 * angle.sin(), 0.2]] {
        let observer = owner(&gas, frequency, x, &source_position)
            .with_wedge(walls, 3, 6)
            .unwrap();
        let field = observer.evaluate(&[tensor()]).unwrap();
        let normal = if x[1] == 0.0 {
            [0.0, 1.0, 0.0]
        } else {
            walls[1].inward_normal()
        };
        let normal_velocity = (0..3).fold(C64::ZERO, |sum, axis| {
            sum + field.particle_velocity_m_s[axis].scale(normal[axis])
        });
        let scale = field
            .particle_velocity_m_s
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        assert!(normal_velocity.abs() < 1e-11 * scale);
        assert!(observer.scope().contains("relaxing gas"));
    }
}

#[test]
fn g1_recorded_stress_keeps_original_coherence_units_and_molecular_propagation() {
    let gas = medium(true);
    let positions = [[0.0; 3], [0.04, -0.02, 0.01]];
    let dt = 1.0 / 8192.0;
    let clock = QuadrupoleSampleClock::new(0.0, dt, 0.0).unwrap();
    let mut record = QuadrupoleWelchStream::new(
        &positions,
        clock,
        64,
        WelchConfig {
            segment_samples: 64,
            hop_samples: 64,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        QuadrupoleRecordingBudget {
            maximum_sources: 2,
            welch: WelchBudget {
                maximum_channels: 12,
                maximum_samples_per_channel: 64,
                maximum_transforms: 12,
                maximum_workspace_bytes: 1_000_000,
            },
        },
        || false,
    )
    .unwrap();
    let a = tensor();
    let b = QuadrupolePhasor {
        stress_moment_n_m: a.stress_moment_n_m.map(|v| v * C64::new(-0.6, 0.3)),
    };
    let bin = 6;
    for sample in 0..64 {
        let phase = TAU * bin as f64 * sample as f64 / 64.0;
        let moments = [a, b].map(|q| {
            q.stress_moment_n_m
                .map(|v| v.re * phase.cos() + v.im * phase.sin())
        });
        record
            .push_sample(&positions, clock.time_at(sample).unwrap(), &moments, || {
                false
            })
            .unwrap();
    }
    let spectrum = record.finish().unwrap();
    let observer = spectrum
        .prepare_molecular_observer(gas.clone(), [0.3, 0.2, 0.1], bin, || false)
        .unwrap();
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    let actual = spectrum
        .evaluate_bin(&observer, bin, work, || false)
        .unwrap();
    let field = observer.evaluate(&[a, b]).unwrap();
    close(
        C64::from_re(actual.mean_square_pressure_pa2),
        C64::from_re(field.pressure_pa.norm_sq() * 0.5),
        3e-12,
    );
    for axis in 0..3 {
        close(
            C64::from_re(actual.active_intensity_w_m2[axis]),
            C64::from_re(field.active_intensity_w_m2[axis]),
            3e-12,
        );
        close(
            C64::from_re(actual.reactive_intensity_w_m2[axis]),
            C64::from_re(field.reactive_intensity_w_m2[axis]),
            3e-12,
        );
    }
    assert_eq!(actual.scope, observer.scope());
    assert_eq!(observer.source_position_m(1), Some(positions[1]));
    assert!(matches!(
        spectrum.prepare_molecular_observer(gas.clone(), [1.0; 3], 0, || false),
        Err(MaterialAeroacError::Aeroac(AeroacError::Welch(_)))
    ));
    assert!(matches!(
        spectrum.prepare_molecular_observer(gas, [1.0; 3], bin, || true),
        Err(MaterialAeroacError::Aeroac(AeroacError::Cancelled { .. }))
    ));
}

#[test]
fn g1_mixed_recorded_mechanisms_share_one_material_and_keep_cross_terms() {
    let gas = medium(true);
    let dt = 1.0 / 8192.0;
    let bin = 6;
    let f = bin as f64 / (64.0 * dt);
    let receiver = [0.3, 0.2, 0.1];
    let compact_position = [[0.0, 0.02, 0.0]];
    let stress_position = [[0.04, 0.03, 0.01]];
    let source = CompactSourcePhasor {
        volume_flow_m3_s: C64::new(0.003, 0.001),
        force_on_fluid_n: [C64::new(0.1, -0.08), C64::from_re(0.07), C64::from_re(0.02)],
    };
    let q = tensor();
    let amplitudes: [C64; 10] = core::array::from_fn(|channel| match channel {
        0 => source.volume_flow_m3_s,
        1..=3 => source.force_on_fluid_n[channel - 1],
        _ => q.stress_moment_n_m[channel - 4],
    });
    let samples: Vec<_> = (0..64)
        .map(|sample| {
            let phase = TAU * bin as f64 * sample as f64 / 64.0;
            amplitudes.map(|a| a.re * phase.cos() + a.im * phase.sin())
        })
        .collect();
    let spectrum = WelchCrossSpectrum::estimate_interleaved(
        &samples,
        dt,
        WelchConfig {
            segment_samples: 64,
            hop_samples: 64,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 10,
            maximum_samples_per_channel: 64,
            maximum_transforms: 10,
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    )
    .unwrap();
    let compact = PreparedCompactSourceObserver::new_molecular(
        gas.clone(),
        f,
        receiver,
        &compact_position,
        1,
    )
    .unwrap();
    let quadrupole = owner(&gas, f, receiver, &stress_position);
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    let wall =
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    for reflected in [false, true] {
        let mut compact = compact.clone();
        let mut quadrupole = quadrupole.clone();
        if reflected {
            compact = compact.with_orthogonal_planes(&[wall], 2).unwrap();
            quadrupole = quadrupole.with_orthogonal_planes(&[wall], 2).unwrap();
        }
        let a = compact.evaluate(&[source]).unwrap();
        let b = quadrupole.evaluate(&[q]).unwrap();
        let p = a.pressure_pa + b.pressure_pa;
        let actual = compact
            .evaluate_welch_bin_with_quadrupoles(&quadrupole, &spectrum, bin, work, || false)
            .unwrap();
        close(
            C64::from_re(actual.mean_square_pressure_pa2),
            C64::from_re(0.5 * p.norm_sq()),
            3e-12,
        );
        let incoherent = 0.5 * (a.pressure_pa.norm_sq() + b.pressure_pa.norm_sq());
        assert!((actual.mean_square_pressure_pa2 - incoherent).abs() > 1e-3 * incoherent);
        for axis in 0..3 {
            let u = a.particle_velocity_m_s[axis] + b.particle_velocity_m_s[axis];
            let intensity = p.scale(0.5) * u.conj();
            close(
                C64::new(
                    actual.active_intensity_w_m2[axis],
                    actual.reactive_intensity_w_m2[axis],
                ),
                intensity,
                3e-12,
            );
        }
        assert!(actual.scope.contains("relaxing gas"));
    }
    let lossless = PreparedCompactSourceObserver::new(
        gas.gas().acoustic_state().unwrap(),
        f,
        receiver,
        &compact_position,
        1,
    )
    .unwrap();
    assert!(
        lossless
            .evaluate_welch_bin_with_quadrupoles(&quadrupole, &spectrum, bin, work, || false)
            .is_err()
    );
    let different_rates = owner(
        &medium_with_rate(true, 2.0 * TAU_R),
        f,
        receiver,
        &stress_position,
    );
    assert_eq!(different_rates.acoustic(), quadrupole.acoustic());
    assert!(
        compact
            .evaluate_welch_bin_with_quadrupoles(&different_rates, &spectrum, bin, work, || false)
            .is_err()
    );
}

#[test]
fn g4_material_and_preparation_refusals_remain_typed_and_transactional() {
    let gas = medium(true);
    assert!(matches!(
        PreparedQuadrupoleObserver::new_molecular(gas.clone(), 0.0, [1.0; 3], &[], 0),
        Err(MaterialAeroacError::Material(_))
    ));
    assert!(
        PreparedQuadrupoleObserver::new_molecular(gas.clone(), 800.0, [0.0; 3], &[[0.0; 3]], 1)
            .is_err()
    );
    assert!(
        PreparedQuadrupoleObserver::new_molecular(gas.clone(), 800.0, [1.0; 3], &[[0.0; 3]], 0)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        PreparedQuadrupoleObserver::new_molecular_cancellable(
            gas,
            800.0,
            [1.0; 3],
            &[[0.0; 3]; 2],
            2,
            || {
                polls += 1;
                polls >= 4
            }
        ),
        Err(MaterialAeroacError::Aeroac(AeroacError::Cancelled { .. }))
    ));
    assert_eq!(polls, 4);
}
