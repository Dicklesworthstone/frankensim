use super::*;

fn observer(positions: &[[f64; 3]]) -> PreparedQuadrupoleObserver {
    PreparedQuadrupoleObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        127.0,
        [1.3, -0.8, 0.6],
        positions,
        positions.len(),
    )
    .unwrap()
}

fn source(scale: f64) -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::from_re(scale),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ],
    }
}

#[test]
fn g3_coherent_cancellation_retains_the_small_resolved_source() {
    let one = observer(&[[0.0; 3]]).evaluate(&[source(1.0)]).unwrap();
    let three = observer(&[[0.0; 3]; 3]);
    let result = three
        .evaluate(&[source(1e16), source(1.0), source(-1e16)])
        .unwrap();
    assert!((result.pressure_pa - one.pressure_pa).abs() < 2e-15 * one.pressure_pa.abs());
    for (actual, expected) in result
        .particle_velocity_m_s
        .into_iter()
        .zip(one.particle_velocity_m_s)
    {
        assert!((actual - expected).abs() < 2e-15 * expected.abs());
    }
    let zero = observer(&[[0.0; 3]; 2])
        .evaluate(&[source(1.0), source(-1.0)])
        .unwrap();
    assert_eq!(zero.pressure_pa, C64::ZERO);
    assert_eq!(zero.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(zero.active_intensity_w_m2, [0.0; 3]);
}

#[test]
fn g1_source_coherence_is_propagated_before_forming_quadratic_power() {
    let obs = observer(&[[0.0; 3]; 2]);
    let budget = CoherenceBudget {
        max_modes: 2,
        max_source_mode_pairs: 4,
    };
    let both = [source(1.0), source(1.0)];
    let first = [source(1.0), source(0.0)];
    let second = [source(0.0), source(1.0)];
    let coherent = obs.evaluate_coherence_modes(&[&both], budget).unwrap();
    let independent = obs
        .evaluate_coherence_modes(&[&first, &second], budget)
        .unwrap();
    assert!(
        (coherent.mean_square_pressure_pa2 / independent.mean_square_pressure_pa2 - 2.0).abs()
            < 1e-14
    );
    let direct = obs.evaluate(&both).unwrap();
    assert!(
        (coherent.mean_square_pressure_pa2 / (0.5 * direct.pressure_pa.norm_sq()) - 1.0).abs()
            < 1e-14
    );
    assert!(
        obs.evaluate_coherence_modes(
            &[&first, &second],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 3
            }
        )
        .is_err()
    );
}

#[test]
fn g4_cancelled_or_malformed_observations_do_not_change_the_prepared_transfer() {
    let obs = observer(&[[0.0; 3]]);
    let before = obs.evaluate(&[source(1.0)]).unwrap();
    assert!(matches!(
        obs.evaluate_cancellable(&[source(1.0)], || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(obs.evaluate(&[]).is_err());
    assert!(obs.evaluate(&[source(f64::NAN)]).is_err());
    let after = obs.evaluate(&[source(1.0)]).unwrap();
    assert_eq!(before.pressure_pa, after.pressure_pa);
    assert_eq!(before.particle_velocity_m_s, after.particle_velocity_m_s);
    assert!(
        PreparedQuadrupoleObserver::new_cancellable(
            obs.acoustic(),
            127.0,
            [1.0; 3],
            &[[0.0; 3]],
            1,
            || true
        )
        .is_err()
    );
    assert!(
        PreparedQuadrupoleObserver::new(obs.acoustic(), 127.0, [1.0; 3], &[[0.0; 3]], 0).is_err()
    );
    assert!(
        PreparedQuadrupoleObserver::new(obs.acoustic(), 127.0, [0.0; 3], &[[0.0; 3]], 1).is_err()
    );
}

#[test]
fn g1_original_harmonic_readout_receives_the_complete_propagated_field() {
    let obs = observer(&[[0.0; 3], [0.2, 0.1, -0.1]]);
    let sources = [source(0.7), source(-0.3)];
    let field = obs.evaluate(&sources).unwrap();
    let line = obs.evaluate_harmonic_line(&sources).unwrap();
    assert_eq!(line.frequency_hz, 127.0);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    let zero = observer(&[]).evaluate(&[]).unwrap();
    assert_eq!(zero.pressure_pa, C64::ZERO);
    assert_eq!(zero.particle_velocity_m_s, [C64::ZERO; 3]);
}
