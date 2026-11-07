use super::*;

#[test]
fn g3_lagrangian_three_fits_sample_accepted_start_and_preserve_source_gauge() {
    let (mut active, _) = independent_owner(0.0);
    let (mut shifted, _) = independent_owner(1e12);
    let mut history = reference_components(active.flow());
    assert_eq!(
        active.lagrangian_dynamic_species_subgrid_config(),
        Some(config())
    );
    assert_eq!(active.lagrangian_dynamic_subgrid_config(), None);
    assert_eq!(active.lagrangian_subgrid_sample_time(), None);
    assert_eq!(fits(&active), (None, None, None));
    let bytes = active.workspace_bytes();
    for duration in [0.001, 0.0005, 0.0015] {
        let start = active.flow().time();
        let (rho, velocity, temperature) = scalar_fields(active.flow());
        let trial = history
            .prepare(
                config().momentum,
                config().heat,
                config().species,
                start.value(),
                &rho,
                &velocity,
                &temperature,
                |component, cell| active.flow().states()[7 * cell + component] / rho[cell],
                || false,
            )
            .unwrap();
        let expected = (
            Some(trial.momentum_report()),
            Some(trial.scalar_report()),
            Some(trial.components_report()),
        );
        let report = active
            .step(Time::new(duration), MolecularStepControl::default(), || {
                false
            })
            .unwrap();
        shifted
            .step(Time::new(duration), MolecularStepControl::default(), || {
                false
            })
            .unwrap();
        conserved(&report);
        assert_eq!(fits(&active), expected);
        assert_eq!(fits(&shifted), expected);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(
            active.flow().time(),
            Time::new(start.value() + report.elapsed.value())
        );
        assert_eq!(active.flow().states(), shifted.flow().states());
        assert!(
            active
                .flow()
                .states()
                .chunks_exact(7)
                .all(|row| row[2] == 0.0)
        );
        assert_eq!(active.workspace_bytes(), bytes);
        trial.commit();
    }
}

#[test]
fn g4_lagrangian_three_histories_retry_shortened_forcing_without_aging_rejections() {
    let force = BodyForce3d::UniformAcceleration([200.0, 0.0, 0.0]);
    let force_controls = ForceControl3d { maximum_mach: 10.0 };
    let make = || {
        let (mut owner, mut work) = independent_owner(0.0);
        owner
            .step_with_forcing(
                Time::new(0.0001),
                MolecularStepControl::default(),
                force,
                force_controls,
                &mut work,
                || false,
            )
            .unwrap();
        (owner, work)
    };
    let (mut active, mut work) = make();
    let start = active.flow().time();
    let prior_states = active.flow().states().to_vec();
    let prior_sample = active.lagrangian_subgrid_sample_time();
    let prior_fits = fits(&active);
    let prior_step = active.last_step().copied();
    let prior_molecular = active.molecular().last_step().copied();
    let prior_advection = active.flow().last_step().copied();
    let prior_loads = active.molecular().boundary_impulses.clone();
    let bytes = active.workspace_bytes();
    let mut polls = 0;
    let accepted = active
        .step_with_forcing(
            Time::new(1.0),
            MolecularStepControl::default(),
            force,
            force_controls,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(accepted.coupling_retries > 0);
    assert!(accepted.elapsed.value() < 1.0);
    conserved(&accepted);
    assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
    let accepted_fits = fits(&active);
    let accepted_states = active.flow().states().to_vec();
    let next = active
        .step(Time::new(0.0005), MolecularStepControl::default(), || false)
        .unwrap();
    let next_fits = fits(&active);

    // This reaches the same physical endpoint without rejected long probes.
    // Its following sample must use that actual duration for all three memories.
    let (mut direct, mut direct_work) = make();
    let short = direct
        .step_with_forcing(
            accepted.elapsed,
            MolecularStepControl::default(),
            force,
            force_controls,
            &mut direct_work,
            || false,
        )
        .unwrap();
    assert_eq!(short.coupling_retries, 0);
    assert_eq!(short.first, accepted.first);
    assert_eq!(short.advection, accepted.advection);
    assert_eq!(short.second, accepted.second);
    assert_eq!(direct.flow().states(), accepted_states);
    assert_eq!(fits(&direct), accepted_fits);
    assert_eq!(
        direct
            .step(Time::new(0.0005), MolecularStepControl::default(), || false,)
            .unwrap(),
        next
    );
    assert_eq!(fits(&direct), next_fits);
    assert_eq!(direct.flow().states(), active.flow().states());

    for cutoff in [polls / 2, polls] {
        let (mut retry, mut retry_work) = make();
        let mut calls = 0;
        assert_eq!(
            retry.step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                force_controls,
                &mut retry_work,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(retry.flow().states(), prior_states);
        assert_eq!(retry.flow().time(), start);
        assert_eq!(retry.last_step().copied(), prior_step);
        assert_eq!(retry.molecular().last_step().copied(), prior_molecular);
        assert_eq!(retry.flow().last_step().copied(), prior_advection);
        assert_eq!(retry.molecular().boundary_impulses, prior_loads);
        assert_eq!(retry.lagrangian_subgrid_sample_time(), prior_sample);
        assert_eq!(fits(&retry), prior_fits);
        assert_eq!(retry.workspace_bytes(), bytes);
        assert_eq!(
            retry
                .step_with_forcing(
                    Time::new(1.0),
                    MolecularStepControl::default(),
                    force,
                    force_controls,
                    &mut retry_work,
                    || false,
                )
                .unwrap(),
            accepted
        );
        assert_eq!(retry.flow().states(), accepted_states);
        assert_eq!(fits(&retry), accepted_fits);
        assert_eq!(
            retry
                .step(Time::new(0.0005), MolecularStepControl::default(), || false,)
                .unwrap(),
            next
        );
        assert_eq!(fits(&retry), next_fits);
        assert_eq!(retry.flow().states(), active.flow().states());
    }
}

fn stabilized_independent_owner() -> Coupled {
    // Reuse the original demanding molecular material and split setup.
    let (flow, mut molecular) = stabilized_owner().into_parts();
    molecular.sgs = None;
    let bytes = Coupled::required_stabilized_workspace_bytes(&flow, &molecular).unwrap()
        + Coupled::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow).unwrap();
    Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, config(), bytes, || false)
        .unwrap()
}

#[test]
fn g4_lagrangian_three_fits_are_held_through_rkl_and_publish_together() {
    let mut active = stabilized_independent_owner();
    let mut replay = stabilized_independent_owner();
    let mut history = reference_components(active.flow());
    for _ in 0..2 {
        let start = active.flow().time();
        let (rho, velocity, temperature) = scalar_fields(active.flow());
        let trial = history
            .prepare(
                config().momentum,
                config().heat,
                config().species,
                start.value(),
                &rho,
                &velocity,
                &temperature,
                |component, cell| active.flow().states()[7 * cell + component] / rho[cell],
                || false,
            )
            .unwrap();
        let expected = (
            Some(trial.momentum_report()),
            Some(trial.scalar_report()),
            Some(trial.components_report()),
        );
        let mut polls = 0;
        let accepted = active
            .step_stabilized(
                Time::new(0.002),
                StabilizedMolecularStepControl::default(),
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert!(accepted.maximum_stages_used > 2);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(fits(&active), expected);
        conserved(&accepted.transport);
        trial.commit();

        let prior_states = replay.flow().states().to_vec();
        let prior_sample = replay.lagrangian_subgrid_sample_time();
        let prior_fits = fits(&replay);
        let prior_step = replay.last_step().copied();
        let mut calls = 0;
        assert_eq!(
            replay.step_stabilized(
                Time::new(0.002),
                StabilizedMolecularStepControl::default(),
                || {
                    calls += 1;
                    calls == polls
                },
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(replay.flow().states(), prior_states);
        assert_eq!(replay.flow().time(), start);
        assert_eq!(replay.last_step().copied(), prior_step);
        assert_eq!(replay.lagrangian_subgrid_sample_time(), prior_sample);
        assert_eq!(fits(&replay), prior_fits);
        assert_eq!(
            replay
                .step_stabilized(
                    Time::new(0.002),
                    StabilizedMolecularStepControl::default(),
                    || false,
                )
                .unwrap(),
            accepted
        );
        assert_eq!(fits(&replay), expected);
        assert_eq!(replay.flow().states(), active.flow().states());
        assert_eq!(replay.workspace_bytes(), active.workspace_bytes());
        assert!(
            active
                .flow()
                .states()
                .chunks_exact(7)
                .all(|row| row[2] == 0.0)
        );
    }
}
