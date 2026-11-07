//! Physical-time SGS belongs to the full reacting split, never to a source stage.

use super::*;
use crate::compressible::spatial::mixture::molecular::DissipationError;
use crate::compressible::spatial::mixture::molecular::sgs::{
    DynamicScalarReport, DynamicSgsReport, LagrangianComponentsConfig, LagrangianScalarConfig,
    LagrangianSgsConfig, LagrangianSubgridConfig, MixtureLagrangianSpeciesSubgridConfig,
    MixtureLagrangianSubgridConfig, SgsUse,
};
use fs_lbm::sgs3::dynamic::lagrangian::scalar::components::LagrangianDynamicComponents;

type Operator = MolecularTransport3d<ConstantTransport>;
type Fits = (
    Option<DynamicSgsReport>,
    Option<DynamicScalarReport>,
    Option<DynamicScalarReport>,
);

fn config() -> MixtureLagrangianSpeciesSubgridConfig {
    MixtureLagrangianSpeciesSubgridConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 0.01,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 0.017,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        species: LagrangianComponentsConfig {
            memory_time: 0.023,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn fixed_config() -> MixtureLagrangianSubgridConfig {
    MixtureLagrangianSubgridConfig {
        transport: LagrangianSubgridConfig {
            model: config().momentum,
            turbulent_prandtl: 0.7,
        },
        turbulent_schmidt: 0.8,
    }
}

fn reacting(reference: f64) -> ReactingEuler3d {
    let (gas, mechanism) = mechanism(5.0, 8.0, reference);
    let grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.20).unwrap();
    coupled(
        field(gas, grid, PERIODIC, |point| {
            let x = std::f64::consts::TAU * point[0] / (9.0 * 0.09);
            let y = std::f64::consts::TAU * point[1] / (7.0 * 0.14);
            let fraction = 0.65 + 0.125 * (x + y).sin();
            (
                [fraction, 1.0 - fraction],
                [0.3 * y.sin(), 0.0, 0.0],
                1.0 + 0.1 * (x + y).sin(),
            )
        }),
        mechanism,
    )
}

fn independent(reference: f64, diffusion: f64) -> (ReactingEuler3d, Operator) {
    let flow = reacting(reference);
    let mut molecular = operator(flow.flow(), diffusion);
    let bytes = flow
        .required_lagrangian_dynamic_species_sgs_workspace_bytes(&molecular)
        .unwrap();
    flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes, || false)
        .unwrap();
    assert_eq!(budget(&flow, &molecular), bytes);
    (flow, molecular)
}

fn fits(operator: &Operator) -> Fits {
    (
        operator.lagrangian_subgrid_report(),
        operator.lagrangian_subgrid_heat_report(),
        operator.lagrangian_subgrid_species_report(),
    )
}

fn reference_history(flow: &MixtureEuler3d) -> LagrangianDynamicComponents {
    LagrangianDynamicComponents::new(
        flow.grid().shape(),
        flow.grid().spacing(),
        flow.gas().species_count(),
        usize::MAX,
    )
    .unwrap()
}

fn expected_start(history: &mut LagrangianDynamicComponents, flow: &MixtureEuler3d) -> Fits {
    let primitive: Vec<_> = (0..flow.cell_count())
        .map(|i| {
            flow.cell(i)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap()
        })
        .collect();
    let density: Vec<_> = primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = primitive.iter().map(|p| p.velocity).collect();
    let temperature: Vec<_> = primitive.iter().map(|p| p.temperature).collect();
    let trial = history
        .prepare(
            config().momentum,
            config().heat,
            config().species,
            flow.time().value(),
            &density,
            &velocity,
            &temperature,
            |component, cell| flow.states()[6 * cell + component] / density[cell],
            || false,
        )
        .unwrap();
    let expected = (
        Some(trial.momentum_report()),
        Some(trial.scalar_report()),
        Some(trial.components_report()),
    );
    trial.commit();
    expected
}

#[derive(Debug, PartialEq)]
struct Published {
    state: Vec<f64>,
    time: Time,
    reacting: Option<ReactiveStepReport3d>,
    advection: Option<StepReport>,
    molecular: Option<MolecularReport>,
    loads: Vec<[f64; 6]>,
    molecular_loads: Vec<[f64; 6]>,
    sample: Option<Time>,
    fits: Fits,
    bytes: usize,
}

fn published(flow: &ReactingEuler3d, molecular: &Operator) -> Published {
    Published {
        state: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        reacting: flow.last_step().copied(),
        advection: flow.flow().last_step().copied(),
        molecular: molecular.last_step().copied(),
        loads: loads(flow),
        molecular_loads: molecular_loads(molecular),
        sample: molecular.lagrangian_subgrid_sample_time(),
        fits: fits(molecular),
        bytes: budget(flow, molecular),
    }
}

#[test]
fn g3_reacting_three_histories_sample_physical_start_with_original_chemical_energy() {
    let (mut flow, mut molecular) = independent(0.0, 0.02);
    let (mut shifted, mut shifted_molecular) = independent(37.0, 0.02);
    let mut history = reference_history(flow.flow());
    let bytes = budget(&flow, &molecular);
    assert_eq!(
        molecular.lagrangian_dynamic_species_subgrid_config(),
        Some(config())
    );
    assert_eq!(molecular.lagrangian_dynamic_subgrid_config(), None);
    assert_eq!(fits(&molecular), (None, None, None));
    for duration in [0.001, 0.0005] {
        let start = flow.flow().time();
        let expected = expected_start(&mut history, flow.flow());
        if start.value() == 0.0 {
            assert_eq!(expected.0.unwrap().maximum_coefficient_squared, 0.0);
            assert!(expected.1.unwrap().maximum_coefficient > 0.0);
            assert!(expected.2.unwrap().maximum_coefficient > 0.0);
        }
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(duration),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .unwrap();
        let shifted_bytes = budget(&shifted, &shifted_molecular);
        let shifted_report = shifted
            .step_with_dissipation(
                &mut shifted_molecular,
                Time::new(duration),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                shifted_bytes,
                || false,
            )
            .unwrap();
        balance(&report);
        balance(&shifted_report);
        assert_eq!(fits(&molecular), expected);
        assert_eq!(molecular.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(
            flow.flow().time().value(),
            start.value() + report.elapsed.value()
        );
        assert_eq!(report.chemistry.cell_solves, 2 * flow.flow().cell_count());
        assert!(report.chemistry.species_change_kg[0] < 0.0);
        assert!(report.chemistry.referenced_energy_change_j > 0.0);
        assert_eq!(
            molecular.last_step(),
            Some(&report.dissipation.unwrap().second)
        );
        assert_eq!(report.elapsed, shifted_report.elapsed);
        for i in 0..flow.flow().cell_count() {
            let left = flow.flow().cell(i).unwrap();
            let right = shifted.flow().cell(i).unwrap();
            for (a, b) in left.partial_densities.iter().zip(right.partial_densities) {
                close(*a, *b, 2e-10);
            }
            let left = left.primitive(flow.flow().gas(), || false).unwrap();
            let right = right.primitive(shifted.flow().gas(), || false).unwrap();
            close(left.temperature, right.temperature, 2e-10);
            for (a, b) in left.velocity.into_iter().zip(right.velocity) {
                close(a, b, 2e-10);
            }
        }
        assert_eq!(budget(&flow, &molecular), bytes);
    }
}

#[test]
fn g4_reacting_histories_retry_shortened_force_and_cancel_without_aging() {
    let force = BodyForce3d::UniformAcceleration([200.0, 0.0, 0.0]);
    let force_controls = ForceControl3d { maximum_mach: 10.0 };
    let make = || {
        let (mut flow, mut molecular) = independent(0.0, 0.02);
        let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
        let bytes = budget(&flow, &molecular) + work.workspace_bytes();
        flow.step_with_dissipation_and_forcing(
            &mut molecular,
            Time::new(0.0001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            force,
            force_controls,
            &mut work,
            || false,
        )
        .unwrap();
        (flow, molecular, work, bytes)
    };
    let (mut flow, mut molecular, mut work, bytes) = make();
    let before = published(&flow, &molecular);
    let mut polls = 0;
    let accepted = flow
        .step_with_dissipation_and_forcing(
            &mut molecular,
            Time::new(1.0),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            force,
            force_controls,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&accepted);
    assert!(accepted.coupling_retries > 0);
    assert!(accepted.elapsed.value() < 1.0);
    assert_eq!(
        molecular.lagrangian_subgrid_sample_time(),
        Some(before.time)
    );
    let accepted_states = flow.flow().states().to_vec();
    let accepted_fits = fits(&molecular);
    let next = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.0005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();

    // A direct short interval and its next sample use the same physical
    // duration, although the other execution discarded longer reacting trials.
    let (mut direct, mut direct_op, mut direct_work, direct_bytes) = make();
    let short = direct
        .step_with_dissipation_and_forcing(
            &mut direct_op,
            accepted.elapsed,
            ReactiveStepControl::default(),
            DissipationControl::default(),
            direct_bytes,
            force,
            force_controls,
            &mut direct_work,
            || false,
        )
        .unwrap();
    assert_eq!(short.coupling_retries, 0);
    assert_eq!(short.transport, accepted.transport);
    assert_eq!(short.chemistry, accepted.chemistry);
    assert_eq!(short.dissipation, accepted.dissipation);
    assert_eq!(direct.flow().states(), accepted_states);
    assert_eq!(fits(&direct_op), accepted_fits);
    assert_eq!(
        direct
            .step_with_dissipation(
                &mut direct_op,
                Time::new(0.0005),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                direct_bytes,
                || false,
            )
            .unwrap(),
        next
    );
    assert_eq!(direct.flow().states(), flow.flow().states());
    assert_eq!(fits(&direct_op), fits(&molecular));

    let (mut replay, mut replay_op, mut replay_work, replay_bytes) = make();
    let mut calls = 0;
    assert_eq!(
        replay.step_with_dissipation_and_forcing(
            &mut replay_op,
            Time::new(1.0),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            replay_bytes,
            force,
            force_controls,
            &mut replay_work,
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(published(&replay, &replay_op), before);
    assert_eq!(
        replay
            .step_with_dissipation_and_forcing(
                &mut replay_op,
                Time::new(1.0),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                replay_bytes,
                force,
                force_controls,
                &mut replay_work,
                || false,
            )
            .unwrap(),
        accepted
    );
    assert_eq!(replay.flow().states(), accepted_states);
    assert_eq!(fits(&replay_op), accepted_fits);
}

#[test]
fn g4_reacting_rkl_holds_all_fits_and_publishes_only_with_both_chemical_halves() {
    let make = || {
        let (mut flow, mut molecular) = independent(0.0, 20.0);
        let bytes = flow
            .required_stabilized_workspace_bytes(&molecular)
            .unwrap();
        // Warm the original lazily allocated RKL buffers before snapshotting.
        flow.step_with_stabilized_dissipation(
            &mut molecular,
            Time::new(0.0001),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
        (flow, molecular, bytes)
    };
    let (mut flow, mut molecular, bytes) = make();
    let (mut replay, mut replay_op, replay_bytes) = make();
    let before = published(&replay, &replay_op);
    let mut polls = 0;
    let report = flow
        .step_with_stabilized_dissipation(
            &mut molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&report.coupled);
    assert!(report.maximum_stages_used > 2);
    assert_eq!(
        report.coupled.chemistry.cell_solves,
        2 * flow.flow().cell_count()
    );
    assert_eq!(
        molecular.lagrangian_subgrid_sample_time(),
        Some(before.time)
    );
    assert!(fits(&molecular).1.unwrap().maximum_coefficient > 0.0);
    assert!(fits(&molecular).2.unwrap().maximum_coefficient > 0.0);
    let mut calls = 0;
    assert_eq!(
        replay.step_with_stabilized_dissipation(
            &mut replay_op,
            Time::new(0.002),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            replay_bytes,
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(published(&replay, &replay_op), before);
    assert_eq!(
        replay
            .step_with_stabilized_dissipation(
                &mut replay_op,
                Time::new(0.002),
                ReactiveStepControl::default(),
                StabilizedDissipationControl::default(),
                replay_bytes,
                || false,
            )
            .unwrap(),
        report
    );
    assert_eq!(published(&replay, &replay_op), published(&flow, &molecular));
}

#[test]
fn g0_reacting_history_attachment_admits_total_capacity_and_refuses_wrong_owners() {
    let mut flow = reacting(0.0);
    let mut molecular = operator(flow.flow(), 0.02);
    let original_bytes = budget(&flow, &molecular);
    let bytes = flow
        .required_lagrangian_dynamic_species_sgs_workspace_bytes(&molecular)
        .unwrap();
    assert_eq!(bytes - original_bytes, 448 * flow.flow().cell_count());
    let before = published(&flow, &molecular);
    assert!(
        flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes - 1, || false,)
            .is_err()
    );
    assert_eq!(published(&flow, &molecular), before);
    assert_eq!(molecular.lagrangian_dynamic_species_subgrid_config(), None);
    let mut calls = 0;
    assert_eq!(
        flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes, || {
            calls += 1;
            calls == 2
        },),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(published(&flow, &molecular), before);
    assert_eq!(molecular.lagrangian_dynamic_species_subgrid_config(), None);
    flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes, || false)
        .unwrap();
    assert!(
        flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes, || false)
            .is_err()
    );

    // The external operator cannot be indexed with another reactor's grid.
    let (gas, mechanism) = mechanism(5.0, 8.0, 0.0);
    let mut wrong = coupled(
        field(
            gas,
            CartesianGrid::new(3, 3, 3, 0.2, 0.2, 0.2).unwrap(),
            PERIODIC,
            |_| ([0.65, 0.35], [0.0; 3], 1.0),
        ),
        mechanism,
    );
    let wrong_before = published(&wrong, &molecular);
    assert!(
        wrong
            .step_with_dissipation(
                &mut molecular,
                Time::new(0.001),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .is_err()
    );
    assert_eq!(published(&wrong, &molecular), wrong_before);

    // Supplying enough memory does not turn a raw molecular clock into an
    // accepted reacting physical interval.
    let mut raw = reacting(0.0).into_flow();
    assert_eq!(
        molecular.advance_unpublished(
            &mut raw,
            Time::new(0.001),
            DissipationControl::default(),
            bytes,
            &mut || false,
        ),
        Err(DissipationError::Invalid(
            "Lagrangian mixture SGS requires a prepared coupled physical step",
        ))
    );
    assert_eq!(molecular.lagrangian_subgrid_sample_time(), None);
    assert_eq!(fits(&molecular), (None, None, None));
    let report = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    balance(&report);
}

#[test]
fn g0_uniform_reacting_fixed_and_independent_memory_preserve_native_molecular_bits() {
    let make = || {
        let (gas, mechanism) = mechanism(40.0, 500.0, 0.0);
        coupled(
            field(
                gas,
                CartesianGrid::new(3, 3, 3, 0.3, 0.3, 0.3).unwrap(),
                PERIODIC,
                |_| ([1.0, 0.0], [0.7, -0.4, 0.2], 1.0),
            ),
            mechanism,
        )
    };
    let mut native = make();
    let mut native_op = operator(native.flow(), 0.02);
    let native_bytes = budget(&native, &native_op);
    let expected = native
        .step_with_dissipation(
            &mut native_op,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            native_bytes,
            || false,
        )
        .unwrap();
    assert!(expected.coupling_retries > 0);
    for independent in [false, true] {
        let mut flow = make();
        let mut molecular = operator(flow.flow(), 0.02);
        let old_bytes = budget(&flow, &molecular);
        let bytes = if independent {
            let bytes = flow
                .required_lagrangian_dynamic_species_sgs_workspace_bytes(&molecular)
                .unwrap();
            flow.attach_lagrangian_dynamic_species_sgs(&mut molecular, config(), bytes, || false)
                .unwrap();
            bytes
        } else {
            let bytes = flow
                .required_lagrangian_dynamic_sgs_workspace_bytes(&molecular)
                .unwrap();
            flow.attach_lagrangian_dynamic_sgs(&mut molecular, fixed_config(), bytes, || false)
                .unwrap();
            assert_eq!(
                molecular.lagrangian_dynamic_subgrid_config(),
                Some(fixed_config())
            );
            bytes
        };
        assert_eq!(
            bytes - old_bytes,
            if independent { 448 } else { 216 } * flow.flow().cell_count()
        );
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(0.2),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .unwrap();
        assert_eq!(report, expected);
        assert_eq!(flow.flow().states(), native.flow().states());
        assert_eq!(molecular.last_step(), native_op.last_step());
        assert_eq!(
            molecular.lagrangian_subgrid_sample_time(),
            Some(Time::new(0.0))
        );
        assert_eq!(fits(&molecular).0.unwrap().maximum_coefficient_squared, 0.0);
        if independent {
            assert_eq!(fits(&molecular).1.unwrap().maximum_coefficient, 0.0);
            assert_eq!(fits(&molecular).2.unwrap().maximum_coefficient, 0.0);
        } else {
            assert_eq!(fits(&molecular).1, None);
            assert_eq!(fits(&molecular).2, None);
        }
    }
}
