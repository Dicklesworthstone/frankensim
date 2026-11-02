use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

fn material(pressure: f64, time_scale: f64, relaxing: bool) -> RelaxationState {
    // Manufactured caloric data, not a measured-air correlation. Pressure
    // changes density and the declared collision time p*tau together.
    let data = if relaxing {
        format!(
            "declared 1 synthetic 1 vibration 1 900 1 power {} 293.15 0 200 500",
            101_325.0 * 0.001 * time_scale
        )
    } else {
        "declared 1 synthetic 0".into()
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

fn observer(surface: &SpherePanels, state: RelaxationState) -> MolecularBemObservation<'_> {
    observer_with_formulation(
        surface,
        state,
        MolecularObservationFormulation::TriangleCbie,
        [1.0, 3_000.0],
    )
}

fn observer_with_formulation(
    surface: &SpherePanels,
    state: RelaxationState,
    formulation: MolecularObservationFormulation,
    band: [f64; 2],
) -> MolecularBemObservation<'_> {
    let n = surface.areas().len();
    MolecularBemObservation::new_with_formulation(
        surface,
        state,
        &[SourceId::StructuralMode("radial velocity [m/s]".into())],
        &[vec![1.0; n]],
        &[
            ReceiverSpec {
                id: "near".into(),
                position_m: [0.0, 0.0, 1.0],
            },
            ReceiverSpec {
                id: "far".into(),
                position_m: [0.0, 0.0, 3.0],
            },
        ],
        0.01,
        band,
        formulation,
        ComplexOptions {
            max_panels: n,
            max_bytes: 32 * n * n,
            max_work: formulation.kernel_terms_per_pair() * (n * n) as u64,
        },
    )
    .unwrap()
}

#[test]
fn fixed_burton_miller_band_crosses_cbie_limit_and_keeps_physical_sphere_response() {
    // The same original BM operator is used below and through the sphere's
    // first interior Dirichlet resonance. No material loss hides that resonance.
    let radius = 0.1;
    let surface = SpherePanels::icosphere(radius, 2).unwrap();
    let state = material(101_325.0, 1.0, false);
    let gas = *state.gas();
    let band = [1.0, 3.25 * gas.sound_speed / radius];
    let model = observer_with_formulation(
        &surface,
        state,
        MolecularObservationFormulation::CentroidBurtonMiller,
        band,
    );
    assert_eq!(
        model.formulation(),
        MolecularObservationFormulation::CentroidBurtonMiller
    );
    assert!(band[1] / gas.sound_speed > model.plain_cbie_limit_rad_per_m());
    let area: f64 = surface.areas().iter().sum();
    for ka in [0.2, 1.0, core::f64::consts::PI] {
        let omega = ka * gas.sound_speed / radius;
        let k = model.source().relaxing_wavenumber(omega).unwrap();
        let q = C64::new(0.0, gas.density * omega);
        let exact_surface = q / C64::new(-1.0 / radius - k.im, k.re);
        let phase = k.scale(1.0 - radius);
        let exact_receiver = exact_surface
            * C64::new(det::cos(phase.re), det::sin(phase.re)).scale(radius * det::exp(-phase.im));
        let reaction = model.impedance_row(omega, 0, || false).unwrap()[0].scale(1.0 / area);
        let field = model.frequency_sample(omega, || false).unwrap();
        assert!(reaction.re > 0.0);
        assert!(
            (reaction - exact_surface).abs() < 0.08 * exact_surface.abs(),
            "ka={ka}: {reaction:?} vs {exact_surface:?}"
        );
        assert!(
            (field[0] - exact_receiver).abs() < 0.12 * exact_receiver.abs(),
            "ka={ka}: {:?} vs {exact_receiver:?}",
            field[0]
        );
    }
}

#[test]
fn formulation_changes_cache_identity_without_changing_material_or_earliest_arrival() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let state = material(101_325.0, 1.0, true);
    let triangle = observer(&surface, state.clone());
    let combined = observer_with_formulation(
        &surface,
        state,
        MolecularObservationFormulation::CentroidBurtonMiller,
        [1.0, 3_000.0],
    );
    assert_eq!(triangle.source(), combined.source());
    assert_eq!(triangle.geometric_delays_s(), combined.geometric_delays_s());
    assert_ne!(triangle.operator_hash(), combined.operator_hash());
    assert!(matches!(
        combined.frequency_sample(500.0, || true),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    assert!(combined.frequency_sample(3_001.0, || false).is_err());
}

#[test]
fn gas_reservoirs_change_both_reaction_and_physical_receiver_phase() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let a = observer(&surface, material(101_325.0, 0.25, true));
    let b = observer(&surface, material(101_325.0, 2.0, true));
    assert_eq!(a.source().gas(), b.source().gas());
    assert_ne!(
        a.operator_hash(),
        b.operator_hash(),
        "equal rho/c cannot hide changed molecular data"
    );
    let w = 600.0;
    let pa = a.frequency_sample(w, || false).unwrap();
    let pb = b.frequency_sample(w, || false).unwrap();
    let za = a.impedance_row(w, 0, || false).unwrap()[0];
    let zb = b.impedance_row(w, 0, || false).unwrap()[0];
    assert!((za - zb).abs() > 1e-3 * za.abs());
    assert!((pa[4] - pb[4]).abs() > 1e-3 * pa[4].abs());
    assert!(za.re > 0.0 && zb.re > 0.0);
    // Far/near transfer of a pulsating sphere: geometric spreading times the
    // SAME complex k. No extra r/c phase or post-sum absorption is permitted.
    let k = a.source().relaxing_wavenumber(w).unwrap();
    let phase = k.scale(2.0);
    let expected =
        C64::new(det::cos(phase.re), det::sin(phase.re)).scale(det::exp(-phase.im) / 3.0);
    assert!(((pa[4] / pa[0]) - expected).abs() < 0.01 * expected.abs());
    let pressure = observer(&surface, material(202_650.0, 0.25, true));
    assert_ne!(pressure.operator_hash(), a.operator_hash());
    assert!((pressure.impedance_row(w, 0, || false).unwrap()[0] - za).abs() > 0.5 * za.abs());
}

#[test]
fn no_reservoirs_recover_original_triangle_boundary_and_receivers() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let state = material(101_325.0, 1.0, false);
    let gas = *state.gas();
    let model = observer(&surface, state);
    let omega = 500.0;
    let fields = vec![C64::ONE; surface.areas().len()];
    let solution = helmholtz::solve_radiation(
        &surface,
        omega / gas.sound_speed,
        Medium::from(&gas),
        &fields,
        helmholtz::Formulation::PlainCbie,
    )
    .unwrap();
    let original = laplace::observe_frequency(
        &surface,
        &solution,
        &[[0.0, 0.0, 1.0], [0.0, 0.0, 3.0]],
        0.01,
    )
    .unwrap();
    let actual = model.frequency_sample(omega, || false).unwrap();
    for (index, field) in original.iter().enumerate() {
        let expected = [
            field.pressure_pa,
            field.gradient_pa_per_m[0],
            field.gradient_pa_per_m[1],
            field.gradient_pa_per_m[2],
        ];
        for (a, b) in actual[4 * index..4 * index + 4].iter().zip(expected) {
            assert!((*a - b).abs() < 1e-12 * field.pressure_pa.abs().max(1.0));
        }
    }
    assert!(
        model
            .frequency_sample(0.0, || false)
            .unwrap()
            .iter()
            .all(|&v| v == C64::ZERO)
    );
}

#[test]
fn molecular_front_delay_caching_band_and_cancellation_remain_source_bound() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let model = observer(&surface, material(101_325.0, 1.0, true));
    let geometry = Geometry::new(&surface, &[[0.0, 0.0, 1.0], [0.0, 0.0, 3.0]], 0.01).unwrap();
    for (&delay, &distance) in model
        .geometric_delays_s()
        .iter()
        .zip(geometry.clearances_m())
    {
        let expected = distance / model.source().frozen_sound_speed_m_s();
        assert!((delay - expected).abs() <= 4.0 * f64::EPSILON * expected);
        assert!(delay < distance / model.source().gas().sound_speed);
    }
    assert!(matches!(
        model.frequency_sample(500.0, || true),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    assert!(model.frequency_sample(3_001.0, || false).is_err());
    let n = surface.areas().len();
    let too_high = MolecularBemObservation::new(
        &surface,
        material(101_325.0, 1.0, true),
        &[SourceId::StructuralMode("body".into())],
        &[vec![1.0; n]],
        &[ReceiverSpec {
            id: "near".into(),
            position_m: [0.0, 0.0, 1.0],
        }],
        0.01,
        [1.0, 100_000.0],
        ComplexOptions {
            max_panels: n,
            max_bytes: 32 * n * n,
            max_work: 512 * (n * n) as u64,
        },
    );
    assert!(
        too_high.is_err(),
        "do not cross a representation resonance or switch quadrature inside a bake"
    );
}

/// G1/G3: spherical propagation in a relaxing medium keeps the actual
/// pressure gradient and its reactive velocity; p/(rho*c) is not momentum.
#[test]
fn harmonic_molecular_fields_preserve_momentum_dispersion_and_absorbed_power() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let model = observer(&surface, material(101_325.0, 0.25, true));
    let omega = 600.0;
    let gas = model.source().gas();
    let k = model.source().relaxing_wavenumber(omega).unwrap();
    assert!(k.im > 0.0);
    let samples = model.frequency_sample(omega, || false).unwrap();
    let prepared = model.prepare_harmonic(omega, || false).unwrap();
    assert_eq!(prepared.producer_hash(), model.operator_hash());
    assert_eq!(prepared.sources(), &model.boundary.sources);
    assert_eq!(prepared.receivers(), &model.boundary.receivers);
    assert!(matches!(
        prepared.producer_uncertainty(),
        ProducerUncertainty::Unresolved { .. }
    ));
    let mut fields = [HarmonicField::default(); 2];
    prepared
        .evaluate_coherent(&[C64::ONE], &mut fields, || false)
        .unwrap();
    let mut statistics = [HarmonicFieldStatistics::default(); 2];
    prepared
        .evaluate_independent_sources(
            &[1.0],
            HarmonicWorkBudget {
                maximum_products: 100,
            },
            &mut statistics,
            || false,
        )
        .unwrap();
    for (receiver, radius) in [1.0, 3.0].into_iter().enumerate() {
        let p = fields[receiver].pressure_pa;
        assert_eq!(p, samples[4 * receiver]);
        for axis in 0..3 {
            let expected =
                samples[4 * receiver + axis + 1] * C64::new(0.0, -1.0 / (gas.density * omega));
            assert_eq!(fields[receiver].particle_velocity_m_s[axis], expected);
        }
        // Differentiate p=A exp(ikr)/r independently: radial velocity is
        // p*(k+i/r)/(rho*omega), including spherical and molecular reactance.
        let radial = p * (k + C64::new(0.0, 1.0 / radius)).scale(1.0 / (gas.density * omega));
        let actual = fields[receiver].particle_velocity_m_s[2];
        assert!((actual - radial).abs() < 0.01 * radial.abs());
        let exact_intensity = p * radial.conj().scale(0.5);
        let intensity = statistics[receiver].complex_intensity_w_m2[2];
        assert!((intensity - exact_intensity).abs() < 0.01 * exact_intensity.abs());
        assert!(intensity.re > 0.0 && intensity.im < 0.0);
        assert!(
            (actual - p.scale(1.0 / (gas.density * gas.sound_speed))).abs() > 0.01 * actual.abs()
        );
    }
    // Outgoing power decreases with radius by exp(-2 Im(k) delta_r).
    // Multiplying pressure power by 1/(rho*c) or adding absorption again
    // would break these jointly observed pressure/velocity relations.
    let power_ratio = 9.0 * statistics[1].complex_intensity_w_m2[2].re
        / statistics[0].complex_intensity_w_m2[2].re;
    let expected = det::exp(-4.0 * k.im);
    assert!((power_ratio - expected).abs() < 0.02 * expected);
    assert!(power_ratio < 1.0);
}

/// G0/G3: structural velocity and volume flow keep their distinct units but
/// interfere before averaging, including pressure/velocity microphone terms.
#[test]
fn molecular_harmonics_preserve_mixed_source_interference_and_microphone_phase() {
    use crate::near_field::FirstOrder;

    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let area: f64 = surface.areas().iter().sum();
    let base = observer(&surface, material(101_325.0, 1.0, true));
    let model = MolecularBemObservation::new(
        &surface,
        base.source().clone(),
        &[
            SourceId::StructuralMode("surface velocity".into()),
            SourceId::ApertureVolumeFlow("volume flow".into()),
        ],
        &[vec![1.0; n], vec![1.0 / area; n]],
        &base.boundary.receivers,
        0.01,
        [1.0, 3_000.0],
        ComplexOptions {
            max_panels: n,
            max_bytes: 32 * n * n,
            max_work: 512 * (n * n) as u64,
        },
    )
    .unwrap();
    let prepared = model.prepare_harmonic(600.0, || false).unwrap();
    let source = C64::new(0.3, 0.4);
    let mut field = [HarmonicField::default(); 2];
    prepared
        .evaluate_coherent(&[source, C64::ZERO], &mut field, || false)
        .unwrap();
    let patterns = [
        FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap(),
        FirstOrder::new(0.0, [0.0, 0.0, -1.0]).unwrap(),
    ];
    let medium = Medium::from(model.source().gas());
    let expected = core::array::from_fn::<_, 2, _>(|i| {
        // Existing microphone calibration uses real equilibrium rho*c. The
        // dispersive velocity itself still came from the momentum equation.
        field[i].pressure_pa.scale(patterns[i].pressure_fraction())
            + field[i].particle_velocity_m_s[2].scale(
                (1.0 - patterns[i].pressure_fraction()) * medium.density * medium.sound_speed,
            )
    });
    let budget = HarmonicWorkBudget {
        maximum_products: 200,
    };
    let mut array = prepared.microphones(&patterns).unwrap();
    let cross = array
        .evaluate_coherence_modes(&[source, C64::ZERO], 1, budget, || false)
        .unwrap();
    for i in 0..2 {
        for j in 0..2 {
            let reference = expected[i] * expected[j].conj().scale(0.5);
            assert!((cross[i * 2 + j] - reference).abs() < 1e-12 * reference.abs());
        }
    }
    assert!(cross[1].im.abs() > 0.01 * cross[1].abs());
    let mut cancelled = [HarmonicField::default(); 2];
    prepared
        .evaluate_coherent(&[source, source.scale(-area)], &mut cancelled, || false)
        .unwrap();
    for i in 0..2 {
        assert!(cancelled[i].pressure_pa.abs() < 1e-12 * field[i].pressure_pa.abs());
        for axis in 0..3 {
            assert!(
                cancelled[i].particle_velocity_m_s[axis].abs()
                    < 1e-12 * field[i].particle_velocity_m_s[2].abs()
            );
        }
    }
}

/// G0/G3: an empty reservoir bank recovers the original physical observation;
/// zero frequency, out-of-band queries and cancellation never publish a bank.
#[test]
fn molecular_harmonic_lossless_limit_and_admission() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let model = observer_with_formulation(
        &surface,
        material(101_325.0, 1.0, false),
        MolecularObservationFormulation::CentroidBurtonMiller,
        [1.0, 3_000.0],
    );
    let omega = 500.0;
    let prepared = model.prepare_harmonic(omega, || false).unwrap();
    let baseline = model.boundary.prepare_harmonic(omega, || false).unwrap();
    let mut actual = [HarmonicField::default(); 2];
    let mut expected = [HarmonicField::default(); 2];
    prepared
        .evaluate_coherent(&[C64::ONE], &mut actual, || false)
        .unwrap();
    baseline
        .evaluate_coherent(&[C64::ONE], &mut expected, || false)
        .unwrap();
    for i in 0..2 {
        assert!(
            (actual[i].pressure_pa - expected[i].pressure_pa).abs()
                < 1e-12 * expected[i].pressure_pa.abs()
        );
        for axis in 0..3 {
            assert!(
                (actual[i].particle_velocity_m_s[axis] - expected[i].particle_velocity_m_s[axis])
                    .abs()
                    < 1e-12 * expected[i].particle_velocity_m_s[2].abs()
            );
        }
    }
    for frequency in [0.0, -1.0, 0.5, 3_001.0, f64::NAN, f64::INFINITY] {
        assert!(model.prepare_harmonic(frequency, || false).is_err());
    }
    assert!(matches!(
        model.prepare_harmonic(omega, || true),
        Err(MolecularObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
}

fn patch_observer(surface: &SpherePanels, state: RelaxationState) -> MolecularBemObservation<'_> {
    patch_observer_at_band(surface, state, [1.0, 3_000.0])
}

fn patch_observer_at_band(
    surface: &SpherePanels,
    state: RelaxationState,
    band: [f64; 2],
) -> MolecularBemObservation<'_> {
    let n = surface.areas().len();
    let mut structure = vec![0.0; n];
    let mut mouth = vec![0.0; n];
    structure[..n / 4].fill(1.0);
    let area: f64 = surface.areas()[n / 4..n / 2].iter().sum();
    mouth[n / 4..n / 2].fill(1.0 / area);
    MolecularBemObservation::new(
        surface,
        state,
        &[
            SourceId::StructuralMode("prescribed patch velocity".into()),
            SourceId::ApertureVolumeFlow("prescribed mouth flow".into()),
        ],
        &[structure, mouth],
        &[ReceiverSpec {
            id: "receiver".into(),
            position_m: [0.2, -0.4, 1.0],
        }],
        0.01,
        band,
        ComplexOptions {
            max_panels: n,
            max_bytes: 32 * n * n,
            max_work: 512 * (n * n) as u64,
        },
    )
    .unwrap()
}

fn patch_wall(surface: &SpherePanels, medium: Medium, resistance: f64) -> LocalBoundary {
    let n = surface.areas().len();
    LocalBoundary::try_new(
        medium,
        vec![
            LocalWall::Rigid,
            LocalWall::Mechanical(fs_phs::WallPin {
                surface_density: 0.1,
                stiffness_per_area: 4.0e5,
                resistance,
            }),
        ],
        (0..n).map(|i| usize::from(i >= n / 2)).collect(),
    )
    .unwrap()
}

/// G0/G3: reconstruct pressure feedback through the ORIGINAL rigid mobility,
/// independently of the Robin assembly. The same source basis then supplies
/// receiver momentum, work-conjugate reaction and passive wall loss.
#[test]
fn molecular_wall_observation_matches_independent_pressure_feedback_and_work() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let state = material(101_325.0, 0.25, true);
    let gas = *state.gas();
    let model = patch_observer(&surface, state)
        .with_local_boundary(patch_wall(&surface, Medium::from(&gas), 200.0))
        .unwrap();
    let omega = 600.0;
    let mut basis = vec![vec![C64::ZERO; n]; n];
    for (i, field) in basis.iter_mut().enumerate() {
        field[i] = C64::ONE;
    }
    let refs: Vec<_> = basis.iter().map(Vec::as_slice).collect();
    let rigid = model
        .exterior
        .solve_batch(
            &surface,
            omega,
            &refs,
            helmholtz::Formulation::PlainCbie,
            &model.limits,
            || false,
        )
        .unwrap();
    // p = Z (v_drive - Y p). This uses the mechanical coefficients directly,
    // not the normalized beta or original complex Robin A/B transformation.
    let wall_y = C64::ONE / C64::new(200.0, 4.0e5 / omega - 0.1 * omega);
    let mut feedback = vec![C64::ZERO; n * n];
    for i in 0..n {
        feedback[i * n + i] = C64::ONE;
        for (j, solution) in rigid.iter().enumerate().skip(n / 2) {
            feedback[i * n + j] = feedback[i * n + j] + solution.pressure_pa()[i] * wall_y;
        }
    }
    let factor = fs_la::eigen_complex::lu_complex(&feedback, n).unwrap();
    let solved = model.boundary_solutions(omega, &mut || false).unwrap();
    let harmonic = model.prepare_harmonic(omega, || false).unwrap();
    let amplitudes = [C64::new(0.7, -0.3), C64::new(-0.003, 0.002)];
    let mut combined_pressure = vec![C64::ZERO; n];
    let mut combined_velocity = vec![C64::ZERO; n];
    let mut combined_drive = vec![C64::ZERO; n];
    let mut expected_receiver = C64::ZERO;
    let mut expected_velocity = [C64::ZERO; 3];
    for (source, drive) in model.boundary.fields.iter().enumerate() {
        let mut expected: Vec<_> = (0..n)
            .map(|i| {
                rigid
                    .iter()
                    .zip(drive)
                    .fold(C64::ZERO, |sum, (solution, &v)| {
                        sum + solution.pressure_pa()[i] * v
                    })
            })
            .collect();
        factor.solve(&mut expected);
        let solution = &solved[source];
        let radiation = &solution.radiation;
        for i in 0..n {
            assert!(
                (radiation.pressure_pa()[i] - expected[i]).abs()
                    < 2e-11 * expected.iter().map(|z| z.abs()).fold(1.0, f64::max)
            );
            let velocity = drive[i]
                - if i < n / 2 {
                    C64::ZERO
                } else {
                    expected[i] * wall_y
                };
            assert!(
                (radiation.velocity_m_s()[i] - velocity).abs() < 2e-11 * velocity.abs().max(1.0)
            );
            if i < n / 2 {
                assert_eq!(radiation.velocity_m_s()[i], drive[i]);
            }
            combined_pressure[i] = combined_pressure[i] + expected[i] * amplitudes[source];
            combined_velocity[i] = combined_velocity[i] + velocity * amplitudes[source];
            combined_drive[i] = combined_drive[i] + drive[i] * amplitudes[source];
        }
        let independent_wall_work: f64 = (n / 2..n)
            .map(|i| 0.5 * 200.0 * (expected[i] * wall_y).abs().powi(2) * surface.areas()[i])
            .sum();
        assert!(independent_wall_work > 0.0);
        assert!(
            (solution.absorbed_power_w - independent_wall_work).abs()
                < 2e-11 * independent_wall_work
        );
        assert!(
            (solution.driving_power_w - solution.absorbed_power_w - radiation.boundary_work_w())
                .abs()
                < 2e-12 * solution.driving_power_w.abs().max(1.0)
        );
        let row = model.impedance_row(omega, source, || false).unwrap();
        for (j, radiating) in solved.iter().enumerate() {
            let expected_reaction = radiating
                .radiation
                .pressure_pa()
                .iter()
                .zip(drive)
                .zip(surface.areas())
                .fold(C64::ZERO, |sum, ((&p, &v), &area)| {
                    sum + p * v.conj().scale(area)
                });
            assert_eq!(row[j], expected_reaction);
        }
        let receiver = radiation
            .observe(&surface, &model.boundary.points, 0.01, || false)
            .unwrap()[0];
        expected_receiver = expected_receiver + receiver.pressure_pa * amplitudes[source];
        for axis in 0..3 {
            expected_velocity[axis] = expected_velocity[axis]
                + receiver.gradient_pa_per_m[axis]
                    * amplitudes[source]
                    * C64::new(0.0, -1.0 / (gas.density * omega));
        }
    }
    let mut actual = [HarmonicField::default(); 1];
    harmonic
        .evaluate_coherent(&amplitudes, &mut actual, || false)
        .unwrap();
    assert!((actual[0].pressure_pa - expected_receiver).abs() < 1e-12 * expected_receiver.abs());
    for axis in 0..3 {
        assert!(
            (actual[0].particle_velocity_m_s[axis] - expected_velocity[axis]).abs()
                < 1e-12 * expected_velocity[axis].abs().max(1e-10)
        );
    }
    let work = |velocity: &[C64]| -> f64 {
        combined_pressure
            .iter()
            .zip(velocity)
            .zip(surface.areas())
            .map(|((&p, &v), &area)| 0.5 * (p * v.conj()).re * area)
            .sum()
    };
    let mixed_wall_work: f64 = (n / 2..n)
        .map(|i| 0.5 * 200.0 * combined_velocity[i].abs().powi(2) * surface.areas()[i])
        .sum();
    assert!(
        (work(&combined_drive) - work(&combined_velocity) - mixed_wall_work).abs()
            < 2e-12 * work(&combined_drive).abs().max(1.0)
    );
}

#[test]
fn molecular_wall_identity_rigid_limit_and_physical_source_admission() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let state = material(101_325.0, 1.0, true);
    let medium = Medium::from(state.gas());
    let rigid = patch_observer(&surface, state.clone());
    let explicit_rigid = patch_observer(&surface, state.clone())
        .with_local_boundary(
            LocalBoundary::try_new(
                medium,
                vec![LocalWall::Rigid],
                vec![0; surface.areas().len()],
            )
            .unwrap(),
        )
        .unwrap();
    assert_eq!(
        rigid.frequency_sample(600.0, || false).unwrap(),
        explicit_rigid.frequency_sample(600.0, || false).unwrap()
    );
    let wall_a = patch_observer(&surface, state.clone())
        .with_local_boundary(patch_wall(&surface, medium, 200.0))
        .unwrap();
    let wall_b = patch_observer(&surface, state.clone())
        .with_local_boundary(patch_wall(&surface, medium, 800.0))
        .unwrap();
    assert_ne!(wall_a.operator_hash(), wall_b.operator_hash());
    assert_ne!(wall_a.operator_hash(), rigid.operator_hash());
    assert_eq!(wall_a.geometric_delays_s(), rigid.geometric_delays_s());
    assert_eq!(wall_a.source(), rigid.source());
    let a = wall_a.frequency_sample(600.0, || false).unwrap();
    let b = wall_b.frequency_sample(600.0, || false).unwrap();
    assert!((a[0] - b[0]).abs() > 1e-3 * a[0].abs());
    assert!(wall_a.frequency_sample(600.0, || true).is_err());
    assert!(wall_a.prepare_harmonic(3_001.0, || false).is_err());
    // Uniform radial motion overlaps the mechanical panels, even though that
    // wall has finite passive impedance. It must not be relabeled as a drive.
    assert!(
        observer(&surface, state.clone())
            .with_local_boundary(patch_wall(&surface, medium, 200.0))
            .is_err()
    );
    assert!(
        patch_observer(&surface, state.clone())
            .with_local_boundary(patch_wall(
                &surface,
                Medium {
                    density: medium.density * 1.1,
                    sound_speed: medium.sound_speed,
                },
                200.0
            ))
            .is_err()
    );
    // A lossless wall resonance inside the declared band is not regularized
    // into a finite admittance or replaced with a rigid panel.
    assert!(
        patch_observer(&surface, state)
            .with_local_boundary(patch_wall(&surface, medium, 0.0))
            .is_err()
    );
}

#[test]
fn molecular_porous_observation_retains_original_material_and_frequency_support() {
    use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
    use fs_material::porous::{PorousCard, PorousFluidLaw};
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let state = material(101_325.0, 0.25, true);
    let gas = *state.gas();
    let medium = Medium::from(&gas);
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let wall = LocalBoundary::try_new(
        medium,
        vec![
            LocalWall::Rigid,
            LocalWall::PorousRigid {
                law: law.clone(),
                fluid: gas.thermoviscous_state().unwrap(),
                thickness_m: 0.03,
            },
        ],
        (0..n).map(|i| usize::from(i >= n / 2)).collect(),
    )
    .unwrap();
    assert!(
        patch_observer(&surface, state.clone())
            .with_local_boundary(wall.clone())
            .is_err()
    );
    let band = [
        core::f64::consts::TAU * 200.0,
        core::f64::consts::TAU * 400.0,
    ];
    let model = patch_observer_at_band(&surface, state.clone(), band)
        .with_local_boundary(wall.clone())
        .unwrap();
    assert_eq!(
        model.local_boundary().unwrap().content_hash(),
        wall.content_hash()
    );
    assert_eq!(model.source(), &state);
    for frequency_hz in [220.0, 370.0] {
        let omega = core::f64::consts::TAU * frequency_hz;
        let native =
            crate::helmholtz::bc::PorousRigidBoundary::new(&law, &gas, omega, 0.03).unwrap();
        let coefficients = wall.normalized_admittance(omega, || false).unwrap();
        assert!(coefficients[..n / 2].iter().all(|&beta| beta == C64::ZERO));
        assert!(
            coefficients[n / 2..]
                .iter()
                .all(|&beta| beta == native.normalized_admittance())
        );
        let solutions = model.boundary_solutions(omega, &mut || false).unwrap();
        assert!(
            solutions
                .iter()
                .all(|solution| solution.absorbed_power_w > 0.0)
        );
        let observed = model.frequency_sample(omega, || false).unwrap();
        let rigid = patch_observer_at_band(&surface, state.clone(), band)
            .frequency_sample(omega, || false)
            .unwrap();
        assert!((observed[0] - rigid[0]).abs() > 1e-3 * rigid[0].abs());
    }
    // The material's positive lower support is not silently extended by the
    // surrounding gas's analytically zero DC radiator response.
    assert!(model.frequency_sample(0.0, || false).is_err());
    assert!(model.frequency_sample(band[0] * 0.5, || false).is_err());
}
