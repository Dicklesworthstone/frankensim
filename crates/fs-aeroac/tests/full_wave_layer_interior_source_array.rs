//! G1/G3: physical layer-source arrays radiate coherently into the incident exterior.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactSourceLimits, PreparedFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::array::{
    FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE, PreparedInteriorSourceFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::{
    PreparedInteriorSourceFullWaveForce, PreparedInteriorSourceFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CompactFieldStatistics, CompactSourcePhasor, PreparedCompactSourceObserver,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn layer(states: [FluidAcousticState; 3], thickness: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        states[0],
        states[1],
        states[2],
        [0.0; 3],
        [0.0, 0.0, 1.0],
        thickness,
    )
    .unwrap()
}
fn limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 1_048_576,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 36,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-8,
            absolute_velocity_tolerance_m_s_per_n: 1e-11,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 36,
        },
    }
}
fn sources() -> [CompactSourcePhasor; 2] {
    [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(0.0001, -0.0003),
            force_on_fluid_n: [
                C64::new(0.2, 0.1),
                C64::new(-0.1, 0.05),
                C64::new(0.1, -0.2),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-0.0002, 0.0001),
            force_on_fluid_n: [
                C64::new(-0.04, 0.03),
                C64::new(0.06, -0.02),
                C64::new(0.08, 0.09),
            ],
        },
    ]
}
fn scaled(source: CompactSourcePhasor, weight: C64) -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s * weight,
        force_on_fluid_n: source.force_on_fluid_n.map(|v| v * weight),
    }
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn check_moments(actual: &CompactFieldStatistics, fields: &[[C64; 4]]) {
    let p2 = fields
        .iter()
        .fold(C64::ZERO, |sum, f| sum + f[0] * f[0].conj())
        .scale(0.5);
    close(C64::from_re(actual.mean_square_pressure_pa2), p2);
    for i in 0..3 {
        let intensity = fields
            .iter()
            .fold(C64::ZERO, |sum, f| sum + f[0] * f[i + 1].conj())
            .scale(0.5);
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            intensity,
        );
        for j in 0..3 {
            let covariance = fields
                .iter()
                .fold(C64::ZERO, |sum, f| sum + f[i + 1] * f[j + 1].conj())
                .scale(0.5);
            close(actual.velocity_covariance_m2_s2[i][j], covariance);
        }
    }
}

#[test]
fn homogeneous_and_invisible_second_face_keep_physical_fluids_and_original_clock() {
    let state = fluid(1.2, 300.0);
    let slab = layer([state; 3], 0.08);
    let positions = [[0.0, 0.0, 0.025], [-0.01, 0.005, 0.04]];
    let receiver = [0.04, -0.012, -0.03];
    let mut budget = limits();
    budget.maximum_total_evaluations = 1;
    budget.monopole.maximum_evaluations = 1;
    budget.force.maximum_evaluations = 1;
    let array = PreparedInteriorSourceFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let original =
        PreparedCompactSourceObserver::new(state, 700.0, receiver, &positions, 2).unwrap();
    let actual = array.evaluate(&sources()).unwrap();
    let expected = original.evaluate(&sources()).unwrap();
    close(actual.pressure_pa, expected.pressure_pa);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            expected.particle_velocity_m_s[i],
        );
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            C64::new(
                expected.active_intensity_w_m2[i],
                expected.reactive_intensity_w_m2[i],
            ),
        );
    }
    let a = array.evaluate_independent_sources(&sources()).unwrap();
    let b = original.evaluate_independent_sources(&sources()).unwrap();
    close(
        C64::from_re(a.mean_square_pressure_pa2),
        C64::from_re(b.mean_square_pressure_pa2),
    );
    for i in 0..3 {
        close(
            C64::new(a.active_intensity_w_m2[i], a.reactive_intensity_w_m2[i]),
            C64::new(b.active_intensity_w_m2[i], b.reactive_intensity_w_m2[i]),
        );
        for j in 0..3 {
            close(
                a.velocity_covariance_m2_s2[i][j],
                b.velocity_covariance_m2_s2[i][j],
            );
        }
    }
    assert_eq!(array.evaluations(), 0);
    assert_eq!(array.layer().layer_state(), state);
    assert_eq!(array.layer().thickness_m(), 0.08);
    assert_eq!(array.frequency_hz(), 700.0);
    assert_eq!(array.observer_position_m(), receiver);
    assert_eq!(array.source_count(), 2);
    assert_eq!(array.source_position_m(1), Some(positions[1]));
    assert_eq!(array.source_position_m(2), None);
    assert!(array.source_diagnostics(2).is_none());
    for i in 0..2 {
        let d = array.source_diagnostics(i).unwrap();
        assert_eq!(d.monopole.evaluations, 0);
        assert_eq!(d.force.evaluations, 0);
    }
    assert_eq!(actual.scope, FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE);
    assert_eq!(a.scope, FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE);
    let line = array.evaluate_harmonic_line(&sources()).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);

    // The physical source fluid is the layer; the exterior receiver retains
    // the incident density in velocity and force transmission.
    let source_fluid = fluid(2.4, 300.0);
    let invisible_exit = PreparedInteriorSourceFullWaveCompactSourceObserver::new(
        layer([state, source_fluid, source_fluid], 0.08),
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let reversed =
        PlanarFluidInterface::new(source_fluid, state, [0.0; 3], [0.0, 0.0, -1.0]).unwrap();
    let original = PreparedFullWaveCompactSourceObserver::new(
        reversed,
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let actual = invisible_exit.evaluate(&sources()).unwrap();
    let expected = original.evaluate(&sources()).unwrap();
    close(actual.pressure_pa, expected.pressure_pa);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            expected.particle_velocity_m_s[i],
        );
    }
    assert_eq!(invisible_exit.evaluations(), 0);
    assert_eq!(invisible_exit.layer().incident_state(), state);
    assert_eq!(invisible_exit.layer().layer_state(), source_fluid);
    for i in 0..2 {
        let d = invisible_exit.source_diagnostics(i).unwrap();
        assert!(d.monopole.reduced_to_single_interface && d.force.reduced_to_single_interface);
    }
}

#[test]
fn finite_layer_source_fields_retain_scalar_composition_coherence_moments_and_work() {
    let slab = layer(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let positions = [[0.0, 0.0, 0.011], [-0.01, 0.005, 0.006]];
    let receiver = [0.04, -0.012, -0.025];
    let amplitudes = sources();
    let array = PreparedInteriorSourceFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut fields = [[C64::ZERO; 4]; 2];
    let mut work = 0;
    for (index, (&position, &amplitude)) in positions.iter().zip(&amplitudes).enumerate() {
        let q = PreparedInteriorSourceFullWaveMonopole::new(
            slab,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedInteriorSourceFullWaveForce::new(
            slab,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        let d = array.source_diagnostics(index).unwrap();
        assert_eq!(d.monopole.evaluations, q.diagnostics().evaluations);
        assert_eq!(d.force.evaluations, f.diagnostics().evaluations);
        assert!(!d.monopole.reduced_to_single_interface && !d.force.reduced_to_single_interface);
        work += q.diagnostics().evaluations + f.diagnostics().evaluations;
        let q = q
            .evaluate_harmonic_line(amplitude.volume_flow_m3_s)
            .unwrap();
        let f = f
            .evaluate_harmonic_line(amplitude.force_on_fluid_n)
            .unwrap();
        fields[index][0] = q.pressure_pa + f.pressure_pa;
        for i in 0..3 {
            fields[index][i + 1] = q.particle_velocity_m_s[i] + f.particle_velocity_m_s[i];
        }
    }
    let actual = array.evaluate(&amplitudes).unwrap();
    let coherent: [C64; 4] = core::array::from_fn(|i| fields[0][i] + fields[1][i]);
    close(actual.pressure_pa, coherent[0]);
    for i in 0..3 {
        close(actual.particle_velocity_m_s[i], coherent[i + 1]);
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            (coherent[0] * coherent[i + 1].conj()).scale(0.5),
        );
    }
    check_moments(
        &array.evaluate_independent_sources(&amplitudes).unwrap(),
        &fields,
    );
    // Shared complex factors correlate physical sources, not separate internal
    // returns. All Q/F cross terms remain inside each source and factor.
    let shared_weight = C64::new(0.3, 0.4);
    let other_weight = C64::new(-0.2, 0.1);
    let first_mode = [amplitudes[0], scaled(amplitudes[1], shared_weight)];
    let second_mode = [
        CompactSourcePhasor::default(),
        scaled(amplitudes[1], other_weight),
    ];
    let moments = array
        .evaluate_coherence_modes(
            &[&first_mode, &second_mode],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 4,
            },
        )
        .unwrap();
    let mode_fields = [
        core::array::from_fn(|i| fields[0][i] + fields[1][i] * shared_weight),
        fields[1].map(|v| v * other_weight),
    ];
    check_moments(&moments, &mode_fields);
    assert_eq!(moments.scope, FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE);
    assert_eq!(array.evaluations(), work);
    assert!(work > 0 && work <= limits().maximum_total_evaluations);
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &positions,
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: work - 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn source_side_domains_nulls_empty_arrays_limits_and_cancellation_are_explicit() {
    let state = fluid(1.2, 300.0);
    let slab = layer([state; 3], 0.08);
    let source = [0.0, 0.0, 0.025];
    let receiver = [0.035, 0.01, -0.03];
    let array = PreparedInteriorSourceFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &[source, source],
        limits(),
        &mut || false,
    )
    .unwrap();
    let first = sources()[0];
    let opposite = scaled(first, C64::from_re(-1.0));
    let null = array.evaluate(&[first, opposite]).unwrap();
    assert_eq!(null.pressure_pa, C64::ZERO);
    assert_eq!(null.particle_velocity_m_s, [C64::ZERO; 3]);
    assert!(
        array
            .evaluate_independent_sources(&[first, opposite])
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
    assert!(array.evaluate(&[first]).is_err());
    assert!(
        array
            .evaluate(&[
                first,
                CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0))
            ])
            .is_err()
    );
    assert!(
        array
            .evaluate_coherence_modes(
                &[&[first, opposite]],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1
                }
            )
            .is_err()
    );
    let empty = PreparedInteriorSourceFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &[],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(empty.evaluations(), 0);
    for observer in [
        [0.0; 3],
        source,
        [0.0, 0.0, 0.08],
        [0.0, 0.0, 0.09],
        [f64::NAN, 0.0, -0.04],
    ] {
        assert!(
            PreparedInteriorSourceFullWaveCompactSourceObserver::new(
                slab,
                700.0,
                observer,
                &[],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    for source in [
        [0.0; 3],
        receiver,
        [0.0, 0.0, 0.08],
        [0.0, 0.0, 0.09],
        [f64::NAN, 0.0, 0.01],
    ] {
        assert!(
            PreparedInteriorSourceFullWaveCompactSourceObserver::new(
                slab,
                700.0,
                receiver,
                &[source],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            layer([state; 3], 0.0),
            700.0,
            receiver,
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            slab,
            f64::NAN,
            receiver,
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    for budget in [
        FullWaveCompactSourceLimits {
            max_sources: 0,
            ..limits()
        },
        FullWaveCompactSourceLimits {
            maximum_total_evaluations: 0,
            ..limits()
        },
    ] {
        assert!(
            PreparedInteriorSourceFullWaveCompactSourceObserver::new(
                slab,
                700.0,
                receiver,
                &[source],
                budget,
                &mut || false,
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &[source],
            limits(),
            &mut || {
                polls += 1;
                polls >= 4
            },
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    // Equal speeds do not make two density jumps an exact single interface.
    let density_layer = layer([state, fluid(5.0, 300.0), fluid(2.4, 300.0)], 0.08);
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            density_layer,
            700.0,
            receiver,
            &[source],
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );

    // DD side admission survives cancellation in the plane dot product even
    // for an empty array, which has no scalar owner to validate the observer.
    let diagonal = PlanarFluidLayer::new(state, state, state, [0.0; 3], [1.0; 3], 1.0).unwrap();
    let outside = [-1e16, -1.0, 1e16];
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            diagonal,
            700.0,
            outside,
            &[],
            limits(),
            &mut || false,
        )
        .is_ok()
    );
    assert!(
        PreparedInteriorSourceFullWaveCompactSourceObserver::new(
            diagonal,
            700.0,
            outside.map(|v| -v),
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
}
