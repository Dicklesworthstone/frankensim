//! Independent continuous/analytic authority and separate live spatial axes.
use core::f64::consts::TAU;
use fs_duct::modal::time_domain::{
    ModalEnd, ModalTimeError, assemble_lossless, assemble_lossless_regional,
};
use fs_duct::modal::{ModeCount, MouthLoad, mm_profile_impedance, modal_reflection_from_impedance};
use fs_duct::{DuctError, LossModel, RadiusProfile, RiccatiSteps, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, PoleOptions};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn harmonic_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 16 * 1024 * 1024,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2.0e-10,
    }
}

fn reflection(z: &[C64], reference: f64, n: usize) -> Vec<C64> {
    modal_reflection_from_impedance(z, &vec![C64::from_re(reference); n], n).unwrap()
}

fn distance(a: &[C64], b: &[C64]) -> f64 {
    det::sqrt(a.iter().zip(b).map(|(&x, &y)| (x - y).abs().powi(2)).sum())
}

// Exercise the actual canonical G/Q/J system through the independent shared
// harmonic owner, not an algebraic duplicate of the assembler's formulas.
fn spatial_reflection(
    profile: &RadiusProfile,
    gas: &GasState,
    n: usize,
    elements: usize,
    end: ModalEnd,
    omega: f64,
) -> Vec<C64> {
    let image =
        assemble_lossless(profile, gas, ModeCount(n), elements, end, 124, || false).unwrap();
    let system = image.system();
    assert_eq!(system.state_dim(), 2 * image.source().coordinates());
    assert_eq!(system.port_dim(), n);
    assert!(system.structure().1.iter().all(|&r| r == 0.0));
    let harmonic = HarmonicModel::from_ode(system, harmonic_budget(), || false).unwrap();
    let factor = harmonic.factor(omega, || false).unwrap();
    let mut impedance = vec![C64::ZERO; n * n];
    for column in 0..n {
        let mut input = vec![C64::ZERO; n];
        input[column] = C64::ONE;
        let response = factor.response(&input, || false).unwrap();
        for row in 0..n {
            impedance[row * n + column] = response.output[row];
        }
        assert_eq!(response.dissipated_power_w, 0.0);
    }
    for row in 0..n {
        for column in 0..n {
            assert!(
                (impedance[row * n + column] - impedance[column * n + row]).abs()
                    < 1.0e-8 * image.source().inlet_reference_impedance()
            );
        }
    }
    reflection(&impedance, image.source().inlet_reference_impedance(), n)
}

fn authority(
    profile: &RadiusProfile,
    gas: &GasState,
    n: usize,
    end: ModalEnd,
    omega: f64,
) -> Vec<C64> {
    let termination = match end {
        ModalEnd::Rigid => Termination::Closed,
        ModalEnd::PressureRelease => Termination::IdealOpen,
        ModalEnd::Ports => panic!("ideal-end authority requires an ideal end"),
    };
    let exact = mm_profile_impedance(
        profile,
        gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Analytic(termination),
        ModeCount(n),
        RiccatiSteps(512),
    )
    .unwrap();
    let z0 =
        gas.density * gas.sound_speed / (core::f64::consts::PI * profile.inlet_radius().powi(2));
    reflection(&exact.impedance_matrix, z0, n)
}

#[test]
fn uniform_all_mode_ideal_ends_follow_propagating_and_evanescent_formulas() {
    let gas = air();
    let (radius, length) = (0.05, 0.04);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let image =
            assemble_lossless(&profile, &gas, ModeCount(3), 16, end, 124, || false).unwrap();
        assert_eq!(
            image.source().nodes_m().len(),
            16 + usize::from(end == ModalEnd::Rigid)
        );
        let z0 = image.source().inlet_reference_impedance();
        for hz in [700.0, 5000.0] {
            let k0 = TAU * hz / gas.sound_speed;
            let mut exact = vec![C64::ZERO; 9];
            for (mode, &gamma) in image.source().radial_eigenvalues().iter().enumerate() {
                let squared = k0 * k0 - (gamma / radius).powi(2);
                let imaginary = if squared > 0.0 {
                    let k = det::sqrt(squared);
                    let tangent = det::tan(k * length);
                    z0 * k0 / k
                        * match end {
                            ModalEnd::Rigid => 1.0 / tangent,
                            ModalEnd::PressureRelease => -tangent,
                            ModalEnd::Ports => unreachable!(),
                        }
                } else {
                    let alpha = det::sqrt(-squared);
                    let e = det::exp(-2.0 * alpha * length);
                    -z0 * k0 / alpha
                        * match end {
                            ModalEnd::Rigid => (1.0 + e) / (1.0 - e),
                            ModalEnd::PressureRelease => (1.0 - e) / (1.0 + e),
                            ModalEnd::Ports => unreachable!(),
                        }
                };
                exact[mode * 3 + mode] = C64::new(0.0, imaginary);
            }
            let actual = spatial_reflection(&profile, &gas, 3, 16, end, TAU * hz);
            let expected = reflection(&exact, z0, 3);
            assert!(distance(&actual, &expected) < 0.025, "end={end:?}, hz={hz}");
            // No accidental coupling survives a truly constant radius.
            for i in 0..3 {
                for j in 0..3 {
                    if i != j {
                        assert!(actual[i * 3 + j].abs() < 1.0e-10);
                    }
                }
            }
        }
    }
}

#[test]
fn axial_ladder_converges_every_matrix_entry_to_original_continuous_riccati() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.12).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let expected = authority(&profile, &gas, 2, end, TAU * 600.0);
        let errors: Vec<_> = [4, 8, 16]
            .into_iter()
            .map(|elements| {
                distance(
                    &spatial_reflection(&profile, &gas, 2, elements, end, TAU * 600.0),
                    &expected,
                )
            })
            .collect();
        assert!(
            errors[0] > 0.03,
            "a live spatial defect is required: {errors:?}"
        );
        assert!(
            errors[1] < 0.6 * errors[0] && errors[2] < 0.4 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.01, "{errors:?}");
        assert!(
            expected[1].abs() > 1.0e-3,
            "the continuous cross-mode reflection must be live"
        );
    }
}

#[test]
fn separate_radial_ladder_changes_and_converges_plane_reflection() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 10.0, 0.12).unwrap();
    let expected = authority(&profile, &gas, 4, ModalEnd::Rigid, TAU * 600.0)[0];
    let values: Vec<_> = [1, 2, 4]
        .into_iter()
        .map(|n| spatial_reflection(&profile, &gas, n, 14, ModalEnd::Rigid, TAU * 600.0)[0])
        .collect();
    let errors: Vec<_> = values.iter().map(|&v| (v - expected).abs()).collect();
    assert!(
        errors[0] > 0.008,
        "a live omitted-mode defect is required: {errors:?}"
    );
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.003, "{errors:?}");
    assert!((values[1] - values[0]).abs() > 0.004);
    assert!((values[2] - values[1]).abs() > 0.001);
}

#[test]
fn source_preserves_all_coordinates_and_refuses_capacity_or_cancellation() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.12).unwrap();
    // All 136 states exceed the shared full exact-flow capacity; never keep
    // only the low eigenmodes, even if the caller supplies a larger ceiling.
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(4),
            16,
            ModalEnd::Rigid,
            1000,
            || false
        ),
        Err(ModalTimeError::Invalid(_))
    ));
    assert!(matches!(
        assemble_lossless(&profile, &gas, ModeCount(2), 4, ModalEnd::Rigid, 19, || {
            false
        }),
        Err(ModalTimeError::Invalid(_))
    ));
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(0),
            4,
            ModalEnd::Rigid,
            124,
            || false
        ),
        Err(ModalTimeError::Duct(DuctError::BadParameter { .. }))
    ));
    let mut calls = 0;
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(2),
            8,
            ModalEnd::Rigid,
            124,
            || {
                calls += 1;
                calls == 3
            }
        ),
        Err(ModalTimeError::Duct(DuctError::Cancelled {
            completed_steps: 1
        }))
    ));
}

#[test]
fn both_face_ports_converge_to_the_full_analytic_cylinder_impedance() {
    let gas = air();
    let (radius, length, modes) = (0.05, 0.04, 3);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    for hz in [700.0, 5000.0] {
        let omega = TAU * hz;
        let k0 = omega / gas.sound_speed;
        let mut errors = Vec::new();
        let mut reflection_error = 0.0;
        for elements in [4, 8, 16] {
            let image = assemble_lossless(
                &profile,
                &gas,
                ModeCount(modes),
                elements,
                ModalEnd::Ports,
                124,
                || false,
            )
            .unwrap();
            let source = image.source();
            let ports = 2 * modes;
            assert_eq!(source.port_count(), ports);
            assert_eq!(source.nodes_m().len(), elements + 1);
            assert_eq!(image.system().port_dim(), ports);
            let z0 = source.inlet_reference_impedance();
            assert_eq!(source.mouth_reference_impedance(), z0);
            let mut expected = vec![C64::ZERO; ports * ports];
            for (mode, &gamma) in source.radial_eigenvalues().iter().enumerate() {
                let squared = k0 * k0 - (gamma / radius).powi(2);
                // P'=i*rho*omega*U/S; the two boundary flows both point IN.
                // Hence Z_00=Z_LL=i*Zc*cot(kL), Z_0L=Z_L0=i*Zc*csc(kL).
                let (diagonal, transfer) = if squared > 0.0 {
                    let k = det::sqrt(squared);
                    let scale = z0 * k0 / k / det::sin(k * length);
                    (scale * det::cos(k * length), scale)
                } else {
                    let alpha = det::sqrt(-squared);
                    let e = det::exp(-alpha * length);
                    let scale = -z0 * k0 / alpha / (1.0 - e * e);
                    (scale * (1.0 + e * e), scale * 2.0 * e)
                };
                expected[mode * ports + mode] = C64::new(0.0, diagonal);
                expected[(modes + mode) * ports + modes + mode] = C64::new(0.0, diagonal);
                expected[mode * ports + modes + mode] = C64::new(0.0, transfer);
                expected[(modes + mode) * ports + mode] = C64::new(0.0, transfer);
            }
            let harmonic =
                HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
            let factor = harmonic.factor(omega, || false).unwrap();
            let mut actual = vec![C64::ZERO; ports * ports];
            for column in 0..ports {
                let mut input = vec![C64::ZERO; ports];
                input[column] = C64::ONE;
                let response = factor.response(&input, || false).unwrap();
                for row in 0..ports {
                    actual[row * ports + column] = response.output[row];
                    if row % modes != column % modes {
                        assert!(response.output[row].abs() < 1.0e-10 * z0);
                    }
                }
            }
            errors.push(distance(&actual, &expected) / z0);
            reflection_error = distance(
                &reflection(&actual, z0, ports),
                &reflection(&expected, z0, ports),
            );
        }
        assert!(errors[0] > 1.0e-4, "live axial defect required: {errors:?}");
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.35 * errors[1],
            "hz={hz}, errors={errors:?}"
        );
        // Keep raw impedance convergence above, including its physical pole
        // sensitivity. The bounded scattering matrix is the reflection QoI.
        assert!(
            reflection_error < 0.025,
            "hz={hz}, reflection_error={reflection_error:e}, impedance_errors={errors:?}"
        );
    }
}

fn warm_air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 600.0, 101_325.0).unwrap()
}

fn two_regions(cold: GasState, warm: GasState, elements: usize) -> Vec<GasState> {
    (0..elements)
        .map(|i| if i < elements / 2 { cold } else { warm })
        .collect()
}

// Independent constant-section pressure/axial-flow transfer, regular also at
// transverse cutoff. No production propagation or assembly helper is called.
fn exact_region(gas: GasState, radius: f64, length: f64, gamma: f64, omega: f64) -> [C64; 4] {
    let k_squared = (omega / gas.sound_speed).powi(2) - (gamma / radius).powi(2);
    let (cosine, sin_over_k) = if k_squared > 0.0 {
        let k = k_squared.sqrt();
        ((k * length).cos(), (k * length).sin() / k)
    } else if k_squared < 0.0 {
        let alpha = (-k_squared).sqrt();
        ((alpha * length).cosh(), (alpha * length).sinh() / alpha)
    } else {
        (1.0, length)
    };
    let inertia = gas.density * omega / (core::f64::consts::PI * radius * radius);
    [
        C64::from_re(cosine),
        C64::new(0.0, -inertia * sin_over_k),
        C64::new(0.0, -k_squared * sin_over_k / inertia),
        C64::from_re(cosine),
    ]
}

fn two_region_transfer(
    cold: GasState,
    warm: GasState,
    radius: f64,
    length: f64,
    gamma: f64,
    omega: f64,
) -> [C64; 4] {
    let a = exact_region(cold, radius, 0.5 * length, gamma, omega);
    let b = exact_region(warm, radius, 0.5 * length, gamma, omega);
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

#[test]
fn repeated_regional_gas_preserves_the_original_full_operator_bit_for_bit() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.12).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease, ModalEnd::Ports] {
        let uniform =
            assemble_lossless(&profile, &gas, ModeCount(2), 8, end, 124, || false).unwrap();
        let regions = vec![gas; 8];
        let regional =
            assemble_lossless_regional(&profile, &regions, ModeCount(2), end, 124, || false)
                .unwrap();
        assert_eq!(regional.source().regional_gases(), Some(regions.as_slice()));
        assert_eq!(regional.source().gas_at_element(7), Some(&gas));
        assert_eq!(regional.source().gas_at_element(8), None);
        assert_eq!(uniform.source().gas_at_element(7), Some(&gas));
        assert_eq!(uniform.source().regional_gases(), None);
        assert_eq!(uniform.source().mass(), regional.source().mass());
        assert_eq!(uniform.source().stiffness(), regional.source().stiffness());
        assert_eq!(
            uniform.source().inverse_mass(),
            regional.source().inverse_mass()
        );
        assert_eq!(uniform.system().structure(), regional.system().structure());
        assert_eq!(
            uniform.source().inlet_reference_impedance(),
            regional.source().inlet_reference_impedance()
        );
        assert_eq!(
            uniform.source().mouth_reference_impedance(),
            regional.source().mouth_reference_impedance()
        );
    }
}

#[test]
fn regional_density_and_cutoff_match_independent_two_medium_pressure_flow_continuity() {
    let (cold, warm) = (air(), warm_air());
    let (radius, length, omega) = (0.05, 0.035, TAU * 5000.0);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let mut errors = Vec::new();
        for elements in [4, 8, 16] {
            let gases = two_regions(cold, warm, elements);
            let image =
                assemble_lossless_regional(&profile, &gases, ModeCount(2), end, 124, || false)
                    .unwrap();
            let source = image.source();
            assert_eq!(source.gas_at_element(elements / 2), Some(&warm));
            assert_eq!(source.gas(), &cold);
            let area = core::f64::consts::PI * radius * radius;
            assert_eq!(
                source.inlet_reference_impedance(),
                cold.density * cold.sound_speed / area
            );
            assert_eq!(
                source.mouth_reference_impedance(),
                warm.density * warm.sound_speed / area
            );
            let gamma = source.radial_eigenvalues()[1];
            assert!(omega > cold.sound_speed * gamma / radius);
            assert!(omega < warm.sound_speed * gamma / radius);
            let mut expected = vec![C64::ZERO; 4];
            for (mode, &gamma) in source.radial_eigenvalues().iter().enumerate() {
                let t = two_region_transfer(cold, warm, radius, length, gamma, omega);
                expected[mode * 2 + mode] = if end == ModalEnd::Rigid {
                    t[0] / t[2]
                } else {
                    t[1] / t[3]
                };
            }
            let model =
                HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
            let factor = model.factor(omega, || false).unwrap();
            let mut actual = vec![C64::ZERO; 4];
            for mode in 0..2 {
                let mut input = vec![C64::ZERO; 2];
                input[mode] = C64::ONE;
                let response = factor.response(&input, || false).unwrap();
                for row in 0..2 {
                    actual[row * 2 + mode] = response.output[row];
                }
                assert!(response.mean_stored_energy_j > 0.0);
                assert_eq!(response.dissipated_power_w, 0.0);
                assert!(
                    response.power_residual_w.abs() < 1e-10 * source.inlet_reference_impedance()
                );
            }
            assert!(actual[1].abs() + actual[2].abs() < 1e-10 * source.inlet_reference_impedance());
            errors.push(distance(
                &reflection(&actual, source.inlet_reference_impedance(), 2),
                &reflection(&expected, source.inlet_reference_impedance(), 2),
            ));
        }
        assert!(
            errors[0] > 1e-4,
            "regional transfer must have a live spatial defect: {errors:?}"
        );
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
            "{end:?}: {errors:?}"
        );
        assert!(errors[2] < 0.025, "{end:?}: {errors:?}");
    }
}

#[test]
fn regional_closed_open_resonance_converges_to_the_physical_interface_root() {
    let (cold, warm) = (air(), warm_air());
    let (radius, length) = (0.02, 0.4);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    // U_in=0 and P_out=0 give D(omega)=0 for the original two-media
    // transfer matrix; this is an independently bracketed physical root.
    let (mut lo, mut hi) = (1.0, 0.9 * core::f64::consts::PI * cold.sound_speed / length);
    assert!(two_region_transfer(cold, warm, radius, length, 0.0, hi)[3].re < 0.0);
    for _ in 0..60 {
        let mid = f64::midpoint(lo, hi);
        if two_region_transfer(cold, warm, radius, length, 0.0, mid)[3].re > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let expected = f64::midpoint(lo, hi);
    let cold_only = core::f64::consts::PI * cold.sound_speed / (2.0 * length);
    assert!((expected / cold_only - 1.0).abs() > 0.1);
    let mut errors = Vec::new();
    for elements in [8, 16, 32] {
        let image = assemble_lossless_regional(
            &profile,
            &two_regions(cold, warm, elements),
            ModeCount(1),
            ModalEnd::PressureRelease,
            124,
            || false,
        )
        .unwrap();
        let model = HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
        let seed_frequency = 0.97 * expected;
        let state = model
            .factor(seed_frequency, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap()
            .state;
        let pole = model
            .resonance_pole(
                C64::new(0.0, -seed_frequency),
                &state,
                &PoleOptions {
                    spectral_scale: expected,
                    step_tolerance: 1e-11,
                    residual_tolerance: 1e-11,
                    max_iterations: 40,
                    max_dimension: 124,
                    max_dense_bytes: 16 * 1024 * 1024,
                    max_work: 40_000_000,
                },
                || false,
            )
            .unwrap();
        assert!(pole.laplace_pole_per_s.re.abs() < 1e-9 * expected);
        assert!(pole.stored_energy_j > 0.0);
        assert_eq!(pole.dissipated_power_w, 0.0);
        errors.push((pole.angular_frequency_rad_s.re / expected - 1.0).abs());
    }
    assert!(errors[0] > 1e-4, "{errors:?}");
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn regional_admission_refuses_incoherent_or_unbalanced_gas_and_cancellation() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.02, 0.0, 0.4).unwrap();
    let assemble = |gases: &[GasState]| {
        assemble_lossless_regional(&profile, gases, ModeCount(1), ModalEnd::Rigid, 124, || {
            false
        })
    };
    assert!(matches!(assemble(&[]), Err(ModalTimeError::Invalid(_))));
    let other_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 600.0, 100_000.0).unwrap();
    assert!(matches!(
        assemble(&[gas, other_pressure]),
        Err(ModalTimeError::Invalid(_))
    ));
    let mut incoherent = gas;
    incoherent.density *= 2.0;
    assert!(matches!(
        assemble(&[gas, incoherent]),
        Err(ModalTimeError::Duct(_))
    ));
    assert!(matches!(
        assemble_lossless_regional(
            &profile,
            &[gas, gas],
            ModeCount(1),
            ModalEnd::Rigid,
            124,
            || true
        ),
        Err(ModalTimeError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(matches!(
        assemble_lossless_regional(
            &profile,
            &[gas, gas],
            ModeCount(1),
            ModalEnd::Rigid,
            2,
            || false
        ),
        Err(ModalTimeError::Invalid(_))
    ));
}
