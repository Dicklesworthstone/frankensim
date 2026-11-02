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
