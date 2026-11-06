use super::*;
use crate::compressible::spatial::molecular::sgs::{
    LagrangianSgsConfig, LagrangianSubgridConfig, SubgridTransport,
};
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};
use fs_lbm::sgs3::dynamic::lagrangian::LagrangianDynamicSgs;

fn memory_config() -> LagrangianSubgridConfig {
    LagrangianSubgridConfig {
        model: LagrangianSgsConfig {
            memory_time: 0.01,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.8,
    }
}

fn active_flow(mu: f64, conductivity: f64) -> MolecularEuler3d<ConstantTransport> {
    let flow = cube(6, mixed);
    let source = law(flow.gas(), mu, 0.0, conductivity);
    MolecularEuler3d::with_lagrangian_dynamic_sgs(
        flow,
        source,
        memory_config(),
        208 * 6usize.pow(3),
        || false,
    )
    .unwrap()
}

fn physical_fields(flow: &Euler3d) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut density = Vec::with_capacity(flow.states().len());
    let mut velocity = Vec::with_capacity(flow.states().len());
    for state in flow.states() {
        let primitive = state.primitive(flow.gas(), || false).unwrap();
        density.push(primitive.density);
        velocity.push(primitive.velocity);
    }
    (density, velocity)
}

fn reference(flow: &Euler3d) -> LagrangianDynamicSgs {
    LagrangianDynamicSgs::new(flow.grid.shape, flow.grid.spacing, usize::MAX).unwrap()
}

#[test]
fn lagrangian_start_samples_match_the_shared_physical_history() {
    let mut flow = active_flow(0.002, 0.004);
    let mut history = reference(flow.flow());
    assert_eq!(
        flow.lagrangian_dynamic_subgrid_config(),
        Some(memory_config())
    );
    assert_eq!(flow.dynamic_subgrid_config(), None);
    assert_eq!(flow.subgrid_config(), None);
    assert_eq!(flow.lagrangian_subgrid_sample_time(), None);
    assert_eq!(flow.lagrangian_subgrid_report(), None);
    for requested in [0.002, 0.001, 0.003] {
        let start = flow.flow().time();
        let (density, velocity) = physical_fields(flow.flow());
        let trial = history
            .prepare(
                memory_config().model,
                start.value(),
                &density,
                &velocity,
                || false,
            )
            .unwrap();
        let expected = trial.report();
        let report = flow
            .step(
                Time::new(requested),
                MolecularStepControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(flow.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(flow.lagrangian_subgrid_report(), Some(expected));
        assert_eq!(
            flow.flow().time(),
            Time::new(start.value() + report.elapsed.value())
        );
        assert_ne!(
            flow.lagrangian_subgrid_sample_time(),
            Some(flow.flow().time())
        );
        trial.commit();
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 4e-13);
        }
    }
}

#[test]
fn lagrangian_faces_use_held_coefficients_with_actual_stress_and_heat() {
    let flow = cube(6, mixed);
    let mut molecular =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.002, 0.0, 0.004), || false).unwrap();
    let mut active =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.002, 0.0, 0.004), || false).unwrap();
    active.sgs = Some(
        SubgridTransport::new_lagrangian(&flow, memory_config(), 208 * 6usize.pow(3)).unwrap(),
    );
    let mut history = active
        .sgs
        .as_mut()
        .unwrap()
        .take_lagrangian_filter()
        .unwrap();
    let width = history.filter_width();
    let trial = active
        .sgs
        .as_mut()
        .unwrap()
        .prepare_lagrangian(&mut history, &flow, &mut || false)
        .unwrap();
    let fit = trial.report();
    let coefficients: Vec<_> = trial.coefficients_squared().collect();
    assert_eq!(coefficients.len(), flow.states().len());
    assert!(fit.maximum_coefficient_squared > 0.0);
    let molecular_rate = prepare(&mut molecular, &flow);
    let active_rate = prepare(&mut active, &flow);
    assert!(active_rate > molecular_rate);
    let mut active_faces = 0;
    for right in 0..flow.states().len() {
        let coordinate = flow.grid.coordinates(right);
        let mut left_coordinate = coordinate;
        left_coordinate[0] = (coordinate[0] + 5) % 6;
        let left = flow.grid.index(left_coordinate);
        let mu_t = face_viscosity_reference(
            &active,
            [coefficients[left], coefficients[right]],
            width,
            left,
            right,
            0,
        );
        if mu_t == 0.0 {
            continue;
        }
        active_faces += 1;
        let pl = active.work.primitive[left];
        let pr = active.work.primitive[right];
        let h = flow.grid.spacing[0];
        let face = flow.grid.cell_face(0, coordinate);
        let a = active.work.flux[0][face];
        let b = molecular.work.flux[0][face];
        let mean =
            std::array::from_fn::<_, 3, _>(|axis| 0.5 * (pl.velocity[axis] + pr.velocity[axis]));
        let heat = |flux: SpatialFaceFlux| {
            flux.energy
                - mean[0] * flux.normal_momentum
                - mean[1] * flux.transverse_momentum[0]
                - mean[2] * flux.transverse_momentum[1]
        };
        close(
            heat(a) - heat(b),
            -mu_t * 3.5 / memory_config().turbulent_prandtl * (pr.temperature - pl.temperature) / h,
            2e-14,
        );
        close(
            a.transverse_momentum[0] - b.transverse_momentum[0],
            -mu_t
                * ((pr.velocity[1] - pl.velocity[1]) / h
                    + 0.5 * active.work.gradient[left][0][1]
                    + 0.5 * active.work.gradient[right][0][1]),
            2e-14,
        );
    }
    assert!(active_faces > 0);
    drop(trial);
    active
        .sgs
        .as_mut()
        .unwrap()
        .restore_lagrangian_filter(history, None);
    assert_eq!(active.sgs.as_ref().unwrap().lagrangian_sample_time(), None);

    // Exercise publication through the real outer transaction as well as the
    // independent face oracle above. Total energy includes viscous stress work.
    let initial = cube(6, mixed);
    let source = law(initial.gas(), 0.002, 0.0, 0.004);
    let mut base = MolecularEuler3d::new(initial, source, || false).unwrap();
    let mut coupled = active_flow(0.002, 0.004);
    let before = coupled.flow().totals();
    base.step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    let report = coupled
        .step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(coupled.lagrangian_subgrid_report(), Some(fit));
    assert_ne!(coupled.flow().states(), base.flow().states());
    close(coupled.flow().totals().energy, before.energy, 4e-13);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 4e-13);
    }
}

#[test]
fn lagrangian_shortened_forced_steps_and_late_cancellation_replay_history() {
    let force = BodyForce3d::UniformAcceleration([20.0, 0.0, 0.0]);
    let cap = ForceControl3d { maximum_mach: 10.0 };
    let make = || {
        let mut flow = active_flow(0.002, 0.004);
        let mut work = flow.flow().forcing_workspace(40 * 6usize.pow(3)).unwrap();
        flow.step_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let initial = cube(6, mixed);
    let mut history = reference(&initial);
    let (density, velocity) = physical_fields(&initial);
    history
        .prepare(memory_config().model, 0.0, &density, &velocity, || false)
        .unwrap()
        .commit();
    let (mut probe, mut work) = make();
    let start = probe.flow().time();
    let (density, velocity) = physical_fields(probe.flow());
    let trial = history
        .prepare(
            memory_config().model,
            start.value(),
            &density,
            &velocity,
            || false,
        )
        .unwrap();
    let expected_fit = trial.report();
    let mut polls = 0;
    let accepted = probe
        .step_with_forcing(
            Time::new(1.0),
            MolecularStepControl::default(),
            force,
            cap,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    trial.commit();
    assert!(accepted.coupling_retries > 0);
    assert!(accepted.elapsed.value() < 1.0);
    assert_eq!(probe.lagrangian_subgrid_sample_time(), Some(start));
    assert_eq!(probe.lagrangian_subgrid_report(), Some(expected_fit));
    let accepted_states = probe.flow().states().to_vec();

    let next_start = probe.flow().time();
    let (density, velocity) = physical_fields(probe.flow());
    let trial = history
        .prepare(
            memory_config().model,
            next_start.value(),
            &density,
            &velocity,
            || false,
        )
        .unwrap();
    let next_fit = trial.report();
    let next = probe
        .step(Time::new(0.0005), MolecularStepControl::default(), || false)
        .unwrap();
    trial.commit();
    assert_eq!(probe.lagrangian_subgrid_sample_time(), Some(next_start));
    assert_eq!(probe.lagrangian_subgrid_report(), Some(next_fit));

    // The same admitted duration, supplied directly, must produce identical
    // physical fields and subsequent memory despite the discarded long trials.
    let (mut direct, mut direct_work) = make();
    let short = direct
        .step_with_forcing(
            accepted.elapsed,
            MolecularStepControl::default(),
            force,
            cap,
            &mut direct_work,
            || false,
        )
        .unwrap();
    assert_eq!(short.coupling_retries, 0);
    assert_eq!(short.first, accepted.first);
    assert_eq!(short.advection, accepted.advection);
    assert_eq!(short.second, accepted.second);
    assert_eq!(direct.flow().states(), accepted_states);
    assert_eq!(direct.lagrangian_subgrid_report(), Some(expected_fit));
    assert_eq!(
        direct
            .step(Time::new(0.0005), MolecularStepControl::default(), || false)
            .unwrap(),
        next
    );
    assert_eq!(direct.lagrangian_subgrid_report(), Some(next_fit));
    assert_eq!(direct.flow().states(), probe.flow().states());

    for cutoff in [polls / 2, polls] {
        let (mut interrupted, mut interrupted_work) = make();
        let old = interrupted.flow().states().to_vec();
        let time = interrupted.flow().time();
        let previous = interrupted.last_step().copied();
        let sample_time = interrupted.lagrangian_subgrid_sample_time();
        let previous_fit = interrupted.lagrangian_subgrid_report();
        let loads = BoundarySide::ALL.map(|side| interrupted.boundary_impulses(side).to_vec());
        let mut calls = 0;
        assert_eq!(
            interrupted.step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                cap,
                &mut interrupted_work,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(interrupted.flow().states(), old);
        assert_eq!(interrupted.flow().time(), time);
        assert_eq!(interrupted.last_step().copied(), previous);
        assert_eq!(interrupted.lagrangian_subgrid_sample_time(), sample_time);
        assert_eq!(interrupted.lagrangian_subgrid_report(), previous_fit);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(interrupted.boundary_impulses(side), loads[i]);
        }
        assert_eq!(
            interrupted
                .step_with_forcing(
                    Time::new(1.0),
                    MolecularStepControl::default(),
                    force,
                    cap,
                    &mut interrupted_work,
                    || false,
                )
                .unwrap(),
            accepted
        );
        assert_eq!(interrupted.flow().states(), accepted_states);
        assert_eq!(interrupted.lagrangian_subgrid_report(), Some(expected_fit));
        assert_eq!(
            interrupted
                .step(Time::new(0.0005), MolecularStepControl::default(), || false)
                .unwrap(),
            next
        );
        assert_eq!(interrupted.lagrangian_subgrid_report(), Some(next_fit));
        for (actual, expected) in interrupted
            .flow()
            .states()
            .iter()
            .zip(probe.flow().states())
        {
            assert_eq!(
                actual.array().map(f64::to_bits),
                expected.array().map(f64::to_bits)
            );
        }
    }
}

#[test]
fn lagrangian_rkl_stages_commit_one_physical_start_sample() {
    let mut flow = active_flow(2.0, 4.0);
    let mut replay = active_flow(2.0, 4.0);
    let mut history = reference(flow.flow());
    for _ in 0..2 {
        let start = flow.flow().time();
        let (density, velocity) = physical_fields(flow.flow());
        let trial = history
            .prepare(
                memory_config().model,
                start.value(),
                &density,
                &velocity,
                || false,
            )
            .unwrap();
        let expected = trial.report();
        let mut polls = 0;
        let report = flow
            .step_stabilized(
                Time::new(0.004),
                MolecularStepControl::default(),
                MolecularStabilization::default(),
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert!(report.maximum_stages_used > 2);
        assert_eq!(flow.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(flow.lagrangian_subgrid_report(), Some(expected));
        assert_eq!(
            flow.flow().time(),
            Time::new(start.value() + report.transport.elapsed.value())
        );
        trial.commit();
        for residual in report.transport.conservation_residual.array() {
            close(residual, 0.0, 4e-13);
        }

        // Abort at the final poll, after multistage work has consumed the held
        // coefficient. Neither the flow nor history may publish that trial.
        let old = replay.flow().states().to_vec();
        let old_sample = replay.lagrangian_subgrid_sample_time();
        let old_fit = replay.lagrangian_subgrid_report();
        let old_step = replay.last_step().copied();
        let mut calls = 0;
        assert_eq!(
            replay.step_stabilized(
                Time::new(0.004),
                MolecularStepControl::default(),
                MolecularStabilization::default(),
                || {
                    calls += 1;
                    calls == polls
                },
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(replay.flow().states(), old);
        assert_eq!(replay.flow().time(), start);
        assert_eq!(replay.lagrangian_subgrid_sample_time(), old_sample);
        assert_eq!(replay.lagrangian_subgrid_report(), old_fit);
        assert_eq!(replay.last_step().copied(), old_step);
        let retried = replay
            .step_stabilized(
                Time::new(0.004),
                MolecularStepControl::default(),
                MolecularStabilization::default(),
                || false,
            )
            .unwrap();
        assert_eq!(retried, report);
        assert_eq!(replay.lagrangian_subgrid_report(), Some(expected));
        assert_eq!(replay.lagrangian_subgrid_sample_time(), Some(start));
        for (actual, expected) in replay.flow().states().iter().zip(flow.flow().states()) {
            assert_eq!(
                actual.array().map(f64::to_bits),
                expected.array().map(f64::to_bits)
            );
        }
    }
}

#[test]
fn lagrangian_budget_refusals_and_uniform_molecular_limit() {
    let flow = cube(6, mixed);
    let bytes =
        MolecularEuler3d::<ConstantTransport>::required_lagrangian_dynamic_sgs_workspace_bytes(
            &flow,
        )
        .unwrap();
    assert_eq!(bytes, 208 * 6usize.pow(3));
    let source = law(flow.gas(), 0.002, 0.0, 0.004);
    assert!(matches!(
        MolecularEuler3d::with_lagrangian_dynamic_sgs(
            flow,
            source,
            memory_config(),
            bytes - 1,
            || false
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    for prandtl in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            LagrangianSubgridConfig {
                turbulent_prandtl: prandtl,
                ..memory_config()
            }
            .validate()
            .is_err()
        );
    }
    let mut forbidden = memory_config();
    forbidden.model.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    for (shape, boundaries) in [([2, 4, 4], PERIODIC), ([4; 3], WALLS)] {
        let [nx, ny, nz] = shape;
        let flow = field(
            gas(0.0),
            CartesianGrid::new(nx, ny, nz, 0.25, 0.25, 0.25).unwrap(),
            boundaries,
            mixed,
        );
        let source = law(flow.gas(), 0.002, 0.0, 0.004);
        assert!(
            MolecularEuler3d::with_lagrangian_dynamic_sgs(
                flow,
                source,
                memory_config(),
                usize::MAX,
                || false
            )
            .is_err()
        );
    }

    // A refused fit must restore the owner too, including on repeated calls.
    let flow = cube(6, mixed);
    let old = flow.states().to_vec();
    let source = law(flow.gas(), 0.002, 0.0, 0.004);
    forbidden = memory_config();
    forbidden.model.maximum_coefficient_squared = 0.0;
    let mut refused =
        MolecularEuler3d::with_lagrangian_dynamic_sgs(flow, source, forbidden, bytes, || false)
            .unwrap();
    for _ in 0..2 {
        assert!(matches!(
            refused.step(Time::new(0.001), MolecularStepControl::default(), || false),
            Err(MolecularError::DynamicSubgrid(
                DynamicSgsError::Coefficient { .. }
            ))
        ));
        assert_eq!(refused.flow().states(), old);
        assert_eq!(refused.flow().time(), Time::new(0.0));
        assert_eq!(refused.lagrangian_subgrid_sample_time(), None);
        assert_eq!(refused.lagrangian_subgrid_report(), None);
    }

    let initial = || cube(4, |_| [1.0, 0.125, -0.0625, 0.03125, 2.0]);
    let flow = initial();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
    let flow = initial();
    let source = law(flow.gas(), 0.002, 0.003, 0.004);
    let mut active = MolecularEuler3d::with_lagrangian_dynamic_sgs(
        flow,
        source,
        memory_config(),
        208 * 4usize.pow(3),
        || false,
    )
    .unwrap();
    for _ in 0..3 {
        let start = active.flow().time();
        let a = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let b = active
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(
            active
                .lagrangian_subgrid_report()
                .unwrap()
                .maximum_coefficient_squared,
            0.0
        );
        for (qa, qb) in active.flow().states().iter().zip(base.flow().states()) {
            assert_eq!(qa.array().map(f64::to_bits), qb.array().map(f64::to_bits));
        }
    }
}

#[test]
fn lagrangian_vorticity_observation_does_not_advance_or_publish_history() {
    let mut observed = active_flow(0.002, 0.004);
    let mut control = active_flow(0.002, 0.004);
    let mut work = observed.vorticity_workspace(usize::MAX).unwrap();
    // Observe both the initial state (no committed history) and an accepted
    // endpoint. Each observation must leave the next physical step unchanged.
    for _ in 0..2 {
        let states = observed.flow().states().to_vec();
        let time = observed.flow().time();
        let sample = observed.lagrangian_subgrid_sample_time();
        let fit = observed.lagrangian_subgrid_report();
        let previous = observed.last_step().copied();
        let mut polls = 0;
        let first = observed
            .vorticity(&mut work, || {
                polls += 1;
                false
            })
            .unwrap()
            .to_vec();
        assert!(
            first
                .iter()
                .any(|budget| budget.stress.iter().any(|value| *value != 0.0))
        );
        assert_eq!(observed.flow().states(), states);
        assert_eq!(observed.flow().time(), time);
        assert_eq!(observed.lagrangian_subgrid_sample_time(), sample);
        assert_eq!(observed.lagrangian_subgrid_report(), fit);
        assert_eq!(observed.last_step().copied(), previous);

        let mut calls = 0;
        assert_eq!(
            observed.vorticity(&mut work, || {
                calls += 1;
                calls == polls
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(observed.flow().states(), states);
        assert_eq!(observed.flow().time(), time);
        assert_eq!(observed.lagrangian_subgrid_sample_time(), sample);
        assert_eq!(observed.lagrangian_subgrid_report(), fit);
        assert_eq!(observed.last_step().copied(), previous);
        assert_eq!(observed.vorticity(&mut work, || false).unwrap(), first);

        let expected = control
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let actual = observed
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            observed.lagrangian_subgrid_sample_time(),
            control.lagrangian_subgrid_sample_time()
        );
        assert_eq!(
            observed.lagrangian_subgrid_report(),
            control.lagrangian_subgrid_report()
        );
        for (actual, expected) in observed.flow().states().iter().zip(control.flow().states()) {
            assert_eq!(
                actual.array().map(f64::to_bits),
                expected.array().map(f64::to_bits)
            );
        }
    }
}
