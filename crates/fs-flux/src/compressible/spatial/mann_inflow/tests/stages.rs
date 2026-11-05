use super::*;

fn inlet(source: &HomogeneousMannField) -> MannReservoirInflow3 {
    MannReservoirInflow3::new(
        source.model(),
        grid(),
        &[BoundarySide::XMin],
        config(),
        MannReservoirObservation3::Point,
        0,
        || false,
    )
    .unwrap()
}

fn mean(time: f64) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0 + time,
        velocity: [0.2 + 0.1 * time, -0.1, 0.05],
        pressure_pa: 1.0,
    }
}

fn prepared(inlet: &MannReservoirInflow3, source: &HomogeneousMannField) -> Euler3d {
    let mut flow = flow(grid());
    inlet
        .apply_frozen(
            &mut flow,
            frozen(source),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            |_, _, _| Ok(mean(0.0)),
            || false,
        )
        .unwrap();
    flow
}

#[test]
fn stage_adapter_matches_native_time_sampled_field_without_advancing_source() {
    let source = field();
    let inlet = inlet(&source);
    let mut adapter = prepared(&inlet, &source);
    let mut native = prepared(&inlet, &source);
    let sample_before = source.sample([0.7, 0.8, 0.9]).unwrap();
    let bytes = adapter.reservoir_profile_workspace_bytes();
    let mut times = Vec::new();
    let report = inlet
        .step_frozen(
            &mut adapter,
            frozen(&source),
            Time::new(0.0001),
            StepControl::default(),
            |time, _, _, _| {
                times.push(time);
                Ok(mean(time))
            },
            || false,
        )
        .unwrap();
    let direct = native
        .step_with_reservoir_profiles(
            Time::new(0.0001),
            StepControl::default(),
            1000,
            |time, _, _, p| {
                // Independent quarter-turn and physical convection map.
                let u = source
                    .sample([1.0 + p[1], 2.0 - p[0] + 0.5 * time, 3.0 + p[2]])
                    .unwrap()
                    .u;
                let mut state = mean(time);
                for (value, fluctuation) in state.velocity.iter_mut().zip([-u[1], u[0], u[2]]) {
                    *value += fluctuation;
                }
                Ok(state)
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.flow.elapsed, direct.elapsed);
    assert_eq!(report.evaluated_faces, times.len());
    assert_eq!(
        report.mode_samples,
        times.len() * source.model().mode_count()
    );
    assert!(times.contains(&0.0));
    assert!(times.contains(&adapter.time().value()));
    for (actual, expected) in adapter.states().iter().zip(native.states()) {
        for (a, b) in actual.array().into_iter().zip(expected.array()) {
            assert!((a - b).abs() < 1e-13);
        }
    }
    assert_eq!(bytes, adapter.reservoir_profile_workspace_bytes());
    let endpoint = adapter.reservoir_face(BoundarySide::XMin, 0).unwrap();
    assert_eq!(endpoint.density, mean(adapter.time().value()).density);
    assert_eq!(source.tick(), 0);
    assert_eq!(sample_before, source.sample([0.7, 0.8, 0.9]).unwrap());
    for residual in report.flow.conservation_residual.array() {
        assert!(residual.abs() < 2e-12);
    }
}

#[test]
fn stage_observation_and_late_mean_refusal_keep_all_accepted_state() {
    let source = field();
    let inlet = inlet(&source);
    let mut flow = prepared(&inlet, &source);
    inlet
        .step_frozen(
            &mut flow,
            frozen(&source),
            Time::new(0.0001),
            StepControl::default(),
            |time, _, _, _| Ok(mean(time)),
            || false,
        )
        .unwrap();
    let states = flow.states().to_vec();
    let time = flow.time();
    let report = *flow.last_step().unwrap();
    let profile = flow.reservoir_face(BoundarySide::XMin, 0).unwrap();
    let loads = flow.boundary_impulses(BoundarySide::XMin).to_vec();
    let failed = inlet.step_frozen(
        &mut flow,
        frozen(&source),
        Time::new(0.0001),
        StepControl::default(),
        |sample_time, _, face, _| {
            if sample_time > time.value() && face == 3 {
                Err(Error::Invalid("authored endpoint outside source domain"))
            } else {
                Ok(mean(sample_time))
            }
        },
        || false,
    );
    assert!(failed.is_err());
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.reservoir_face(BoundarySide::XMin, 0), Some(profile));
    assert_eq!(flow.boundary_impulses(BoundarySide::XMin), loads);

    let evolving = FixedStepMannField::from_snapshot(source.clone(), 119, 0.001, 10_000).unwrap();
    assert!(matches!(
        inlet.step_frozen(
            &mut flow,
            MannReservoirSource3::Evolving(&evolving),
            Time::new(0.0001),
            StepControl::default(),
            |time, _, _, _| Ok(mean(time)),
            || false,
        ),
        Err(MannReservoirError3::Input(_))
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&report));
}

#[test]
fn whole_step_work_budget_requires_both_stages_before_sampling() {
    let source = field();
    let mut config = config();
    config.maximum_mode_samples = 4 * 12; // One complete plane, not two.
    let inlet = MannReservoirInflow3::new(
        source.model(),
        grid(),
        &[BoundarySide::XMin],
        config,
        MannReservoirObservation3::Point,
        0,
        || false,
    )
    .unwrap();
    let mut flow = prepared(&inlet, &source);
    let before = flow.states().to_vec();
    let mut callbacks = 0;
    let failed = inlet.step_frozen(
        &mut flow,
        frozen(&source),
        Time::new(0.0001),
        StepControl::default(),
        |time, _, _, _| {
            callbacks += 1;
            Ok(mean(time))
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(MannReservoirError3::Work {
            required: 96,
            maximum: 48
        })
    ));
    assert_eq!(callbacks, 0);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(flow.states(), before);
}

#[test]
fn mixture_stage_profile_preserves_authored_species_at_accepted_time() {
    let source = field();
    let inlet = inlet(&source);
    let gas = crate::compressible::mixture::tests::standard(0.2);
    let grid = grid();
    let count: usize = grid.shape().iter().product();
    let mut states = vec![0.0; count * 7];
    for row in states.chunks_exact_mut(7) {
        super::super::super::mixture::state_into(&gas, &[0.6, 0.4, 0.0], [0.0; 3], 1.0, row)
            .unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, states, BOUNDARIES, bytes, || false).unwrap();
    let baseline = |time: f64, densities: &mut [f64]| {
        densities.copy_from_slice(&[0.6 + time, 0.4 - time, 0.0]);
        Ok(MixtureReservoirFaceState {
            velocity: [0.2, 0.0, 0.0],
            pressure_pa: 1.0,
        })
    };
    inlet
        .apply_mixture(
            &mut flow,
            frozen(&source),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            |_, _, _, densities| baseline(0.0, densities),
            || false,
        )
        .unwrap();
    let bytes = flow.reservoir_profile_workspace_bytes();
    let report = inlet
        .step_mixture(
            &mut flow,
            frozen(&source),
            Time::new(0.0001),
            StepControl::default(),
            |time, _, _, _, densities| baseline(time, densities),
            || false,
        )
        .unwrap();
    let face = flow.reservoir_face(BoundarySide::XMin, 0).unwrap();
    assert_eq!(
        face.partial_densities,
        &[0.6 + flow.time().value(), 0.4 - flow.time().value(), 0.0]
    );
    assert_eq!(face.pressure_pa, 1.0);
    assert_eq!(bytes, flow.reservoir_profile_workspace_bytes());
    assert!(report.mode_samples <= config().maximum_mode_samples);
    assert!(report.flow.absolute_energy_residual_j.abs() < 2e-12);
    assert_eq!(source.tick(), 0);
}
