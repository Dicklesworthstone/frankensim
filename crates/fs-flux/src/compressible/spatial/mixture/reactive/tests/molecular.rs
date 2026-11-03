use super::*;
use crate::compressible::spatial::mixture::molecular::{
    DissipationControl, MolecularTransport3d, SpatialTransport,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn operator(flow: &MixtureEuler3d, diffusion: f64) -> MolecularTransport3d<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "manufactured-spatial-reacting-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.03),
            bulk_viscosity: DynViscosity::new(0.01),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.02),
        },
        vec![0.0, diffusion, diffusion, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    MolecularTransport3d::new(flow, source, bytes, || false).unwrap()
}

fn budget(flow: &ReactingEuler3d, molecular: &MolecularTransport3d<ConstantTransport>) -> usize {
    flow.workspace_bytes() + molecular.workspace_bytes()
}

fn loads(flow: &ReactingEuler3d) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            flow.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn molecular_loads(molecular: &MolecularTransport3d<ConstantTransport>) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            molecular.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

#[cfg(feature = "experimental-sgs")]
#[test]
fn sgs_mixture_reactive_split_uses_explicit_turbulent_mixing_and_conserves_energy() {
    use crate::compressible::spatial::mixture::molecular::sgs::{
        MixtureSubgridConfig, SgsModel, SgsUse, SubgridConfig,
    };
    let make = || {
        let (gas, source) = mechanism(0.0, 2.0, 0.0);
        coupled(
            field(
                gas,
                CartesianGrid::new(6, 2, 2, 1.0 / 6.0, 0.3, 0.4).unwrap(),
                PERIODIC,
                |x| {
                    let phase = 2.0 * std::f64::consts::PI * x[0];
                    let first = 0.5 + 0.2 * phase.sin();
                    (
                        [first, 1.0 - first],
                        [0.0, 0.2 * phase.sin(), 0.2 * phase.cos()],
                        1.0,
                    )
                },
            ),
            source,
        )
    };
    let mut baseline = make();
    let mut base_op = operator(baseline.flow(), 0.02);
    let mut flow = make();
    let source =
        SpatialTransport::new(flow.flow().gas().clone(), base_op.law().clone(), || false).unwrap();
    let required =
        MolecularTransport3d::required_workspace_bytes_with_sgs(flow.flow(), &source).unwrap();
    let mut op = MolecularTransport3d::with_sgs(
        flow.flow(),
        source,
        MixtureSubgridConfig {
            transport: SubgridConfig {
                model: SgsModel::Smagorinsky { cs: 0.25 },
                filter_width_factor: 1.0,
                turbulent_prandtl: 0.7,
                maximum_eddy_to_molecular: 100.0,
                usage: SgsUse::Exploration,
            },
            turbulent_schmidt: 0.6,
        },
        required,
        || false,
    )
    .unwrap();
    let bytes = budget(&flow, &op);
    let base_bytes = budget(&baseline, &base_op);
    let report = flow
        .step_with_dissipation(
            &mut op,
            Time::new(0.002),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    let base_report = baseline
        .step_with_dissipation(
            &mut base_op,
            Time::new(0.002),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            base_bytes,
            || false,
        )
        .unwrap();
    balance(&report);
    assert_eq!(report.elapsed, base_report.elapsed);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(budget(&flow, &op), bytes);
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    assert_eq!(op.last_step(), Some(&split.second));
    let variance = |flow: &ReactingEuler3d| {
        (0..flow.flow().cell_count())
            .map(|i| (flow.flow().cell(i).unwrap().partial_densities[0] - 0.5).powi(2))
            .sum::<f64>()
    };
    assert!(variance(&flow) < variance(&baseline));
}

#[test]
fn g3_uniform_reacting_molecular_field_retains_chemistry_and_heat_driven_cfl_restarts() {
    let mut reference = uniform();
    let expected = reference
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut flow = uniform();
    let mut molecular = operator(flow.flow(), 0.02);
    let bytes = budget(&flow, &molecular);
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
    balance(&report);
    assert!(report.coupling_retries > 0);
    assert_eq!(report.elapsed, expected.elapsed);
    assert_eq!(report.transport, expected.transport);
    assert_eq!(report.chemistry, expected.chemistry);
    assert_eq!(flow.flow().states(), reference.flow().states());
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    assert_eq!(molecular.last_step(), Some(&split.second));
    assert_eq!(budget(&flow, &molecular), bytes);
}

fn reacting_diffusion(n: usize) -> f64 {
    let (gas, source) = mechanism(0.0, 2.0, 0.0);
    let dx = 1.0 / n as f64;
    let z = std::f64::consts::PI * dx;
    let average = (z.sin() / z).powi(3);
    let diffusion = 0.02;
    let velocity = [0.4, 0.25, -0.15];
    // Original reversible A <-> B source gives d(Y_A-1/2)/dt=-4(Y_A-1/2).
    // Three equal wave numbers add 12*pi^2*D to that exact decay rate.
    let decay_rate = 4.0 + diffusion * 12.0 * std::f64::consts::PI.powi(2);
    let exact = |x: [f64; 3], time: f64| {
        0.5 + 0.2
            * average
            * (-decay_rate * time).exp()
            * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2] - 0.5 * time)).sin()
    };
    let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
    let mut flow = coupled(
        field(gas, grid, PERIODIC, |x| {
            let a = exact(x, 0.0);
            ([a, 1.0 - a], velocity, 1.0)
        }),
        source,
    );
    let mut molecular = operator(flow.flow(), diffusion);
    let bytes = budget(&flow, &molecular);
    let end = 0.03;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .unwrap();
        balance(&report);
        let split = report.dissipation.unwrap();
        assert!(split.first.maximum_parabolic_number <= 0.40000000000001);
        assert!(split.second.maximum_parabolic_number <= 0.40000000000001);
        steps += 1;
        assert!(steps < 100);
    }
    let mut error = 0.0;
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.pressure, 1.0, 2e-10);
        close(p.temperature, 1.0, 2e-10);
        for axis in 0..3 {
            close(p.velocity[axis], velocity[axis], 2e-10);
        }
        error += (q.partial_densities[0] - exact(grid.cell_center(i).unwrap(), end)).abs();
    }
    error / flow.flow().cell_count() as f64
}

#[test]
fn g1_oblique_three_dimensional_reaction_advection_diffusion_refines_to_exact_solution() {
    let coarse = reacting_diffusion(12);
    let fine = reacting_diffusion(24);
    let order = (coarse / fine).ln() / 2.0_f64.ln();
    eprintln!("3-D reacting diffusion L1: n12={coarse:e}, n24={fine:e}, order={order:e}");
    assert!(coarse < 0.01);
    assert!(
        order > 1.8,
        "coarse={coarse:e}, fine={fine:e}, order={order:e}"
    );
}

#[test]
fn g4_all_molecular_reacting_phases_and_support_loads_publish_atomically() {
    let create = || {
        let (gas, source) = mechanism(5.0, 8.0, 0.0);
        let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.25, 0.3).unwrap();
        let input = field(
            gas,
            grid,
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let a = 0.5 + 0.2 * x[0] - 0.1 * x[2];
                (
                    [a, 1.0 - a],
                    [0.2 + 0.1 * x[1], -0.15, 0.12],
                    1.0 + 0.1 * x[2],
                )
            },
        );
        // An oversized incoming buffer must remain counted as it moves between
        // the two flow fields and the molecular operator during accepted swaps.
        let mut states = Vec::with_capacity(input.states.len() + 31);
        states.extend_from_slice(&input.states);
        let input = MixtureEuler3d::new(
            input.gas().clone(),
            *input.grid(),
            states,
            input.boundaries(),
            crate::compressible::spatial::mixture::MAX_MIXTURE_WORKSPACE_BYTES,
            || false,
        )
        .unwrap();
        let flow = coupled(input, source);
        let molecular = operator(flow.flow(), 0.02);
        (flow, molecular)
    };
    let (mut expected, mut expected_molecular) = create();
    let (mut flow, mut molecular) = create();
    let bytes = budget(&flow, &molecular);
    let controls = DissipationControl {
        maximum_step_s: 0.0005,
        ..DissipationControl::default()
    };
    for (flow, molecular) in [
        (&mut expected, &mut expected_molecular),
        (&mut flow, &mut molecular),
    ] {
        flow.step_with_dissipation(
            molecular,
            Time::new(0.001),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || false,
        )
        .unwrap();
        assert_eq!(budget(flow, molecular), bytes);
    }
    let mut polls = 0;
    let reference = expected
        .step_with_dissipation(
            &mut expected_molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&reference);
    let split = reference.dissipation.unwrap();
    assert_ne!(
        reference.boundary_change.momentum,
        reference.transport.boundary_change.momentum
    );
    let mut outward = [0.0; 3];
    for load in loads(&expected) {
        assert_eq!(load[0], 0.0);
        assert_eq!(load[1], 0.0);
        assert_eq!(load[5], 0.0);
        for axis in 0..3 {
            outward[axis] += load[axis + 2];
        }
    }
    for axis in 0..3 {
        close(
            reference.after.momentum[axis] - reference.before.momentum[axis] + outward[axis],
            0.0,
            3e-12,
        );
        close(
            reference.boundary_change.momentum[axis],
            reference.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    let states = flow.flow().states().to_vec();
    let time = flow.flow().time();
    let report = *flow.last_step().unwrap();
    let wave_report = *flow.flow().last_step().unwrap();
    let molecular_report = *molecular.last_step().unwrap();
    let accepted_loads = loads(&flow);
    let accepted_molecular_loads = molecular_loads(&molecular);
    let mut count = 0;
    assert_eq!(
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    for refused in [
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.001),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl {
                maximum_step_s: 1e-6,
                maximum_steps: 1,
                ..controls
            },
            bytes,
            || false,
        ),
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.001),
            ReactiveStepControl::default(),
            controls,
            bytes - 1,
            || false,
        ),
    ] {
        assert!(refused.is_err());
    }
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.flow().last_step(), Some(&wave_report));
    assert_eq!(molecular.last_step(), Some(&molecular_report));
    assert_eq!(loads(&flow), accepted_loads);
    assert_eq!(molecular_loads(&molecular), accepted_molecular_loads);
    assert_eq!(budget(&flow, &molecular), bytes);
    let replay = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(replay, reference);
    assert_eq!(flow.flow().states(), expected.flow().states());
    assert_eq!(loads(&flow), loads(&expected));
    assert_eq!(
        molecular_loads(&molecular),
        molecular_loads(&expected_molecular)
    );
    assert_eq!(budget(&flow, &molecular), bytes);
}

#[test]
fn g1_forced_molecular_reacting_restart_replays_one_admitted_interval() {
    let mut flow = uniform();
    let mut molecular = operator(flow.flow(), 0.02);
    let mut scratch = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let bytes = budget(&flow, &molecular) + scratch.workspace_bytes();
    let acceleration = [400.0, -200.0, 100.0];
    let forces = vec![acceleration; flow.flow().cell_count()];
    let cap = ForceControl3d { maximum_mach: 0.8 };
    let controls = ReactiveStepControl::default();
    let mut polls = 0;
    let report = flow
        .step_with_dissipation_and_forcing(
            &mut molecular,
            Time::new(0.2),
            controls,
            DissipationControl::default(),
            bytes,
            BodyForce3d::Acceleration(&forces),
            cap,
            &mut scratch,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&report);
    assert!(report.coupling_retries > 0);
    assert_eq!(report.chemistry.cell_solves, 2 * flow.flow().cell_count());
    let split = report.dissipation.unwrap();
    let dt = report.elapsed.value();
    assert_eq!(split.first.elapsed.value(), 0.5 * dt);
    assert_eq!(split.second.elapsed, split.first.elapsed);
    assert_eq!(molecular.last_step(), Some(&split.second));
    let mut rate: f64 = 0.0;
    for index in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(index)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        assert!(
            p.velocity[0].hypot(p.velocity[1]).hypot(p.velocity[2])
                <= cap.maximum_mach * p.sound_speed * (1.0 + 1e-13)
        );
        rate = rate.max(
            (0..3)
                .map(|axis| {
                    (p.velocity[axis].abs() + p.sound_speed) / flow.flow().grid().spacing()[axis]
                })
                .sum(),
        );
        for axis in 0..3 {
            close(
                p.velocity[axis],
                [0.7, -0.4, 0.2][axis] + dt * acceleration[axis],
                2e-12,
            );
        }
    }
    assert!(rate * dt <= controls.transport.cfl * (1.0 + 1e-13));
    for axis in 0..3 {
        close(
            report.source_change.impulse[axis],
            report.before.mass() * dt * acceleration[axis],
            2e-12,
        );
    }
    let work = report.before.mass()
        * (0..3)
            .map(|axis| {
                dt * acceleration[axis] * ([0.7, -0.4, 0.2][axis] + 0.5 * dt * acceleration[axis])
            })
            .sum::<f64>();
    close(report.source_change.work, work, 2e-12);
    let mut replay = uniform();
    let mut replay_molecular = operator(replay.flow(), 0.02);
    let mut replay_scratch = replay.flow().forcing_workspace(usize::MAX).unwrap();
    let expected = replay
        .step_with_dissipation_and_forcing(
            &mut replay_molecular,
            report.elapsed,
            controls,
            DissipationControl::default(),
            bytes,
            BodyForce3d::Acceleration(&forces),
            cap,
            &mut replay_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(expected.coupling_retries, 0);
    assert_eq!(expected.transport, report.transport);
    assert_eq!(expected.chemistry, report.chemistry);
    assert_eq!(expected.dissipation, report.dissipation);
    assert_eq!(expected.source_change, report.source_change);
    assert_eq!(replay.flow().states(), flow.flow().states());
    assert_eq!(loads(&replay), loads(&flow));
    let mut cancelled_flow = uniform();
    let initial = cancelled_flow.flow().states().to_vec();
    let mut cancelled_molecular = operator(cancelled_flow.flow(), 0.02);
    let mut cancelled_scratch = cancelled_flow.flow().forcing_workspace(usize::MAX).unwrap();
    let initial_loads = loads(&cancelled_flow);
    let initial_molecular_loads = molecular_loads(&cancelled_molecular);
    let mut count = 0;
    assert_eq!(
        cancelled_flow.step_with_dissipation_and_forcing(
            &mut cancelled_molecular,
            Time::new(0.2),
            controls,
            DissipationControl::default(),
            bytes,
            BodyForce3d::Acceleration(&forces),
            cap,
            &mut cancelled_scratch,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(cancelled_flow.flow().states(), initial);
    assert_eq!(cancelled_flow.flow().time(), Time::new(0.0));
    assert_eq!(cancelled_flow.last_step(), None);
    assert_eq!(cancelled_molecular.last_step(), None);
    assert_eq!(loads(&cancelled_flow), initial_loads);
    assert_eq!(
        molecular_loads(&cancelled_molecular),
        initial_molecular_loads
    );
}

#[test]
fn g0_zero_force_retains_native_spatial_reacting_and_molecular_paths() {
    let create = || {
        let (gas, mechanism) = mechanism(5.0, 8.0, 0.0);
        coupled(
            field(
                gas,
                CartesianGrid::new(3, 2, 2, 0.2, 0.25, 0.3).unwrap(),
                Boundaries3d {
                    x: WALLS,
                    y: WALLS,
                    z: WALLS,
                },
                |x| {
                    let a = 0.5 + 0.2 * x[0] - 0.1 * x[2];
                    (
                        [a, 1.0 - a],
                        [0.2 + 0.1 * x[1], -0.15, 0.12],
                        1.0 + 0.1 * x[2],
                    )
                },
            ),
            mechanism,
        )
    };
    for with_molecular in [false, true] {
        let mut native = create();
        let mut flow = create();
        let mut native_molecular = operator(native.flow(), 0.02);
        let mut molecular = operator(flow.flow(), 0.02);
        let mut scratch = flow.flow().forcing_workspace(usize::MAX).unwrap();
        let bytes = budget(&flow, &molecular) + scratch.workspace_bytes();
        let zero = vec![[0.0; 3]; flow.flow().cell_count()];
        // Inactive forcing must retain the original path even when an active
        // source with this control would reject the initial velocity.
        let cap = ForceControl3d { maximum_mach: 1e-6 };
        let expected = if with_molecular {
            native
                .step_with_dissipation(
                    &mut native_molecular,
                    Time::new(0.001),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    || false,
                )
                .unwrap()
        } else {
            native
                .step(Time::new(0.001), ReactiveStepControl::default(), || false)
                .unwrap()
        };
        let report = if with_molecular {
            flow.step_with_dissipation_and_forcing(
                &mut molecular,
                Time::new(0.001),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                BodyForce3d::Acceleration(&zero),
                cap,
                &mut scratch,
                || false,
            )
            .unwrap()
        } else {
            flow.step_with_forcing(
                Time::new(0.001),
                ReactiveStepControl::default(),
                BodyForce3d::Acceleration(&zero),
                cap,
                &mut scratch,
                || false,
            )
            .unwrap()
        };
        assert_eq!(report, expected);
        assert_eq!(report.source_change, ForceChange3d::default());
        for (left, right) in flow.flow().states().iter().zip(native.flow().states()) {
            assert_eq!(left.to_bits(), right.to_bits());
        }
        assert_eq!(flow.flow().time(), native.flow().time());
        assert_eq!(loads(&flow), loads(&native));
        assert_eq!(molecular.last_step(), native_molecular.last_step());
        assert_eq!(
            molecular_loads(&molecular),
            molecular_loads(&native_molecular)
        );
        if with_molecular {
            assert!(matches!(
                flow.step_with_dissipation_and_forcing(
                    &mut molecular,
                    Time::new(0.001),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes - 1,
                    BodyForce3d::Acceleration(&zero),
                    cap,
                    &mut scratch,
                    || false,
                ),
                Err(ReactiveError::Invalid(_))
            ));
            assert_eq!(flow.flow().states(), native.flow().states());
            assert_eq!(flow.last_step(), Some(&report));
            assert_eq!(molecular.last_step(), native_molecular.last_step());
            assert_eq!(loads(&flow), loads(&native));
            assert_eq!(
                molecular_loads(&molecular),
                molecular_loads(&native_molecular)
            );
        }
    }
}
