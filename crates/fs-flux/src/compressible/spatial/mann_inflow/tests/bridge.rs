use super::*;
use crate::compressible::spatial::mixture::molecular as mixing;
use crate::compressible::spatial::molecular as viscous;
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

mod forcing;
mod reactive;

const DRAWS: usize = 100_000;

fn source() -> FixedStepMannField {
    FixedStepMannField::from_snapshot(field(), 211, 1.0 / 1024.0, 10_000).unwrap()
}
fn bridge(source: &FixedStepMannField) -> MannOuBridge {
    MannOuBridge::prepare(source, 173, 128, 10_000, 1000, || false).unwrap()
}
fn inlet(source: &FixedStepMannField) -> MannReservoirInflow3 {
    MannReservoirInflow3::new(
        source.model(),
        grid(),
        &BoundarySide::ALL,
        config(),
        MannReservoirObservation3::FaceAverage,
        10_000,
        || false,
    )
    .unwrap()
}
fn mean(time: f64) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0 + time,
        velocity: [0.2, -0.1 + time, 0.05],
        pressure_pa: 1.0,
    }
}
fn mixture_mean(
    time: f64,
    densities: &mut [f64],
) -> Result<MixtureReservoirFaceState, MixtureError> {
    densities.copy_from_slice(&[0.6 + time, 0.4 - time, 0.0]);
    Ok(MixtureReservoirFaceState {
        velocity: mean(time).velocity,
        pressure_pa: 1.0,
    })
}
fn prepared(inlet: &MannReservoirInflow3, frame: &MannOuBridgeFrame) -> Euler3d {
    let mut flow = flow(grid());
    if frame.time_s() > 0.0 {
        flow.step(Time::new(frame.time_s()), StepControl::default(), || false)
            .unwrap();
    }
    inlet
        .apply_frozen(
            &mut flow,
            MannReservoirSource3::BridgeFrame(frame),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            |_, _, _| Ok(mean(frame.time_s())),
            || false,
        )
        .unwrap();
    flow
}
fn prepared_mixture(inlet: &MannReservoirInflow3, frame: &MannOuBridgeFrame) -> MixtureEuler3d {
    let gas = crate::compressible::mixture::tests::standard(0.2);
    let grid = grid();
    let mut states = vec![0.0; grid.len() * 7];
    for row in states.chunks_exact_mut(7) {
        crate::compressible::spatial::mixture::state_into(
            &gas,
            &[0.6, 0.4, 0.0],
            [0.0; 3],
            1.0,
            row,
        )
        .unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, states, BOUNDARIES, bytes, || false).unwrap();
    if frame.time_s() > 0.0 {
        flow.step(Time::new(frame.time_s()), StepControl::default(), || false)
            .unwrap();
    }
    inlet
        .apply_mixture(
            &mut flow,
            MannReservoirSource3::BridgeFrame(frame),
            0.001,
            MannReservoirAction3::Prepare {
                maximum_profile_bytes: 100_000,
            },
            |_, _, _, densities| mixture_mean(frame.time_s(), densities),
            || false,
        )
        .unwrap();
    flow
}
fn law(source: IdealGasEquilibriumV1) -> ConstantTransport {
    let n = source.models().len();
    let mut binary = vec![0.001; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        source,
        "native-mann-bridge-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.002),
            bulk_viscosity: DynViscosity::new(0.0001),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.002),
        },
        binary,
    )
    .unwrap()
}
fn coupled(flow: Euler3d) -> viscous::MolecularEuler3d<ConstantTransport> {
    let models = flow.gas().source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    viscous::MolecularEuler3d::new(flow, law(source), || false)
        .unwrap()
        .with_reservoir_exchange()
}
fn coupled_mixture(flow: MixtureEuler3d) -> mixing::MolecularEuler3d<ConstantTransport> {
    let transport =
        mixing::SpatialTransport::new(flow.gas().clone(), law(flow.gas().source().clone()), || {
            false
        })
        .unwrap();
    let bytes = mixing::MolecularTransport3d::required_workspace_bytes(&flow, &transport).unwrap();
    let transport = mixing::MolecularTransport3d::new(&flow, transport, bytes, || false).unwrap();
    let bytes =
        mixing::MolecularEuler3d::required_stabilized_workspace_bytes(&flow, &transport).unwrap();
    mixing::MolecularEuler3d::new(flow, transport, bytes, || false).unwrap()
}

#[test]
fn evolving_frozen_faces_match_direct_original_bridge_synthesis() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let mut frame = bridge.frame(10_000, || false).unwrap();
    let mut direct_frame = bridge.frame(10_000, || false).unwrap();
    let mut adapter = prepared(&inlet, &frame);
    let mut native = prepared(&inlet, &frame);
    let duration = Time::new(bridge.dt_s() / 4.0);
    let bytes = adapter.reservoir_profile_workspace_bytes();
    let report = inlet
        .step_frozen_bridge(
            &mut adapter,
            &bridge,
            &mut frame,
            duration,
            StepControl::default(),
            DRAWS,
            |time, _, _, _| Ok(mean(time)),
            || false,
        )
        .unwrap();
    native
        .step_with_reservoir_profiles(
            duration,
            StepControl::default(),
            1000,
            |time, side, _, p| {
                bridge
                    .sample_frame(time, &mut direct_frame, DRAWS, || false)
                    .unwrap();
                // Independent quarter-turn and exact rectangular observation.
                let u = direct_frame
                    .sample_filtered(
                        [1.0 + p[1], 2.0 - p[0], 3.0 + p[2]],
                        inlet.spatial_filter(side).unwrap(),
                    )
                    .unwrap()
                    .u;
                let mut state = mean(time);
                for (v, u) in state.velocity.iter_mut().zip([-u[1], u[0], u[2]]) {
                    *v += u;
                }
                Ok(state)
            },
            || false,
        )
        .unwrap();
    assert_eq!(adapter.states(), native.states());
    assert_eq!(adapter.time(), duration);
    assert_eq!(report.frame_updates, 2);
    assert_eq!(
        report.normal_draws,
        bridge.required_normal_draws(duration.value()).unwrap()
    );
    assert_eq!(report.evaluated_faces, 2 * inlet.faces);
    assert_eq!(
        report.mode_samples,
        report.evaluated_faces * source.model().mode_count()
    );
    assert_eq!(adapter.reservoir_profile_workspace_bytes(), bytes);
    assert_eq!(source.checkpoint(), checkpoint);
    for r in report.flow.conservation_residual.array() {
        assert!(r.abs() < 3e-12);
    }
}

#[test]
fn evolving_mixture_faces_match_native_species_and_enthalpy_owner() {
    let source = source();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let mut frame = bridge.frame(10_000, || false).unwrap();
    let mut direct_frame = bridge.frame(10_000, || false).unwrap();
    let mut adapter = prepared_mixture(&inlet, &frame);
    let mut native = prepared_mixture(&inlet, &frame);
    let duration = Time::new(bridge.dt_s() / 4.0);
    let report = inlet
        .step_mixture_bridge(
            &mut adapter,
            &bridge,
            &mut frame,
            duration,
            StepControl::default(),
            DRAWS,
            |time, _, _, _, densities| mixture_mean(time, densities),
            || false,
        )
        .unwrap();
    native
        .step_with_reservoir_profiles(
            duration,
            StepControl::default(),
            1000,
            |time, side, _, p, densities| {
                bridge
                    .sample_frame(time, &mut direct_frame, DRAWS, || false)
                    .unwrap();
                let u = direct_frame
                    .sample_filtered(
                        [1.0 + p[1], 2.0 - p[0], 3.0 + p[2]],
                        inlet.spatial_filter(side).unwrap(),
                    )
                    .unwrap()
                    .u;
                let mut state = mixture_mean(time, densities)?;
                for (v, u) in state.velocity.iter_mut().zip([-u[1], u[0], u[2]]) {
                    *v += u;
                }
                Ok(state)
            },
            || false,
        )
        .unwrap();
    assert_eq!(adapter.states(), native.states());
    let face = adapter.reservoir_face(BoundarySide::XMin, 0).unwrap();
    assert_eq!(
        face.partial_densities,
        &[0.6 + duration.value(), 0.4 - duration.value(), 0.0]
    );
    assert_eq!(report.frame_updates, 2);
    assert_eq!(
        report.normal_draws,
        bridge.required_normal_draws(duration.value()).unwrap()
    );
    assert!(report.flow.absolute_energy_residual_j.abs() < 3e-12);
}

#[test]
fn molecular_ssp_and_rkl_share_the_bridge_and_roll_back_exhausted_draws() {
    for stabilized in [false, true] {
        let source = source();
        let checkpoint = source.checkpoint();
        let bridge = bridge(&source);
        let inlet = inlet(&source);
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut flow = coupled(prepared(&inlet, &frame));
        let mut reference = coupled(prepared(&inlet, &frame));
        let original = flow.flow().states().to_vec();
        let profile = flow.flow().reservoir_face(BoundarySide::XMin, 0);
        let duration = Time::new(bridge.dt_s() / 4.0);
        let step =
            |flow: &mut viscous::MolecularEuler3d<_>, frame: &mut MannOuBridgeFrame, budget| {
                if stabilized {
                    inlet
                        .step_stabilized_molecular_frozen_bridge(
                            flow,
                            &bridge,
                            frame,
                            duration,
                            viscous::MolecularStepControl::default(),
                            viscous::MolecularStabilization::default(),
                            budget,
                            |time, _, _, _| Ok(mean(time)),
                            || false,
                        )
                        .map(|r| r.flow.transport)
                } else {
                    inlet
                        .step_molecular_frozen_bridge(
                            flow,
                            &bridge,
                            frame,
                            duration,
                            viscous::MolecularStepControl::default(),
                            budget,
                            |time, _, _, _| Ok(mean(time)),
                            || false,
                        )
                        .map(|r| r.flow)
                }
            };
        assert!(matches!(
            step(&mut flow, &mut frame, 0),
            Err(MannReservoirError3::BridgeWork { .. })
        ));
        assert_eq!(flow.flow().states(), original);
        assert_eq!(flow.flow().time().value(), 0.0);
        assert_eq!(flow.flow().reservoir_face(BoundarySide::XMin, 0), profile);
        assert!(flow.last_step().is_none());
        let actual = step(&mut flow, &mut frame, DRAWS).unwrap();
        let mut independent = bridge.frame(10_000, || false).unwrap();
        let expected = step(&mut reference, &mut independent, DRAWS).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.flow().states(), reference.flow().states());
        assert_eq!(frame.time_s(), flow.flow().time().value());
        assert!(actual.first.boundary_change.energy != 0.0);
        for r in actual.conservation_residual.array() {
            assert!(r.abs() < 3e-12);
        }
        assert_eq!(source.checkpoint(), checkpoint);
    }
}

#[test]
fn molecular_species_ssp_and_rkl_use_one_replayable_turbulent_bath() {
    for stabilized in [false, true] {
        let source = source();
        let bridge = bridge(&source);
        let inlet = inlet(&source);
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut flow = coupled_mixture(prepared_mixture(&inlet, &frame));
        let duration = Time::new(bridge.dt_s() / 4.0);
        let mut times = Vec::new();
        let mut baseline = |time, side, face, _, densities: &mut [f64]| {
            if side == BoundarySide::XMin && face == 0 {
                times.push(time);
            }
            mixture_mean(time, densities)
        };
        let report = if stabilized {
            inlet
                .step_stabilized_molecular_mixture_bridge(
                    &mut flow,
                    &bridge,
                    &mut frame,
                    duration,
                    mixing::StabilizedMolecularStepControl::default(),
                    DRAWS,
                    &mut baseline,
                    || false,
                )
                .unwrap()
                .flow
                .transport
        } else {
            inlet
                .step_molecular_mixture_bridge(
                    &mut flow,
                    &bridge,
                    &mut frame,
                    duration,
                    mixing::MolecularStepControl::default(),
                    DRAWS,
                    &mut baseline,
                    || false,
                )
                .unwrap()
                .flow
        };
        assert!(times.windows(2).any(|t| t[1] < t[0])); // Split operators replay physical times.
        assert_eq!(frame.time_s(), flow.flow().time().value());
        let face = flow.flow().reservoir_face(BoundarySide::XMin, 0).unwrap();
        assert_eq!(
            face.partial_densities,
            &[0.6 + duration.value(), 0.4 - duration.value(), 0.0]
        );
        for value in report.conservation_residual.species_masses {
            assert!(value.abs() < 3e-12);
        }
        for value in report.conservation_residual.momentum {
            assert!(value.abs() < 3e-12);
        }
        assert!(report.conservation_residual.energy.abs() < 3e-12);
        assert!(report.absolute_energy_residual_j.abs() < 3e-10);
    }
}

#[test]
fn bridge_interval_caps_native_steps_without_advancing_its_source() {
    let source = source();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let mut frame = bridge.frame(10_000, || false).unwrap();
    let mut flow = prepared(&inlet, &frame);
    let report = inlet
        .step_frozen_bridge(
            &mut flow,
            &bridge,
            &mut frame,
            Time::new(2.0 * bridge.dt_s()),
            StepControl::default(),
            0,
            |time, _, _, _| Ok(mean(time)),
            || false,
        )
        .unwrap();
    assert_eq!(report.maximum_duration.value(), bridge.dt_s());
    assert_eq!(flow.time().value(), bridge.end_time_s());
    assert_eq!(report.normal_draws, 0);
    assert_eq!(source.step(), 0);
    let states = flow.states().to_vec();
    assert!(
        inlet
            .step_frozen_bridge(
                &mut flow,
                &bridge,
                &mut frame,
                Time::new(bridge.dt_s()),
                StepControl::default(),
                DRAWS,
                |time, _, _, _| Ok(mean(time)),
                || false,
            )
            .is_err()
    );
    assert_eq!(flow.states(), states);
}

#[test]
fn nonzero_epoch_split_endpoint_is_exactly_the_retained_bridge_endpoint() {
    let mut source = FixedStepMannField::from_snapshot(field(), 211, 1e-6, 10_000).unwrap();
    source.advance_to(6, 1000, || false).unwrap();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let start = bridge.start_time_s();
    let end = bridge.end_time_s();
    let duration = end - start;
    assert!((start + duration / 2.0) + duration / 2.0 > end);
    for stabilized in [false, true] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut frozen = coupled(prepared(&inlet, &frame));
        if stabilized {
            inlet
                .step_stabilized_molecular_frozen_bridge(
                    &mut frozen,
                    &bridge,
                    &mut frame,
                    Time::new(2e-6),
                    viscous::MolecularStepControl::default(),
                    viscous::MolecularStabilization::default(),
                    DRAWS,
                    |time, _, _, _| Ok(mean(time)),
                    || false,
                )
                .unwrap();
        } else {
            inlet
                .step_molecular_frozen_bridge(
                    &mut frozen,
                    &bridge,
                    &mut frame,
                    Time::new(2e-6),
                    viscous::MolecularStepControl::default(),
                    DRAWS,
                    |time, _, _, _| Ok(mean(time)),
                    || false,
                )
                .unwrap();
        }
        assert_eq!(frozen.flow().time().value().to_bits(), end.to_bits());
        assert_eq!(frame.time_s().to_bits(), end.to_bits());
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut mixture = coupled_mixture(prepared_mixture(&inlet, &frame));
        if stabilized {
            inlet
                .step_stabilized_molecular_mixture_bridge(
                    &mut mixture,
                    &bridge,
                    &mut frame,
                    Time::new(2e-6),
                    mixing::StabilizedMolecularStepControl::default(),
                    DRAWS,
                    |time, _, _, _, densities| mixture_mean(time, densities),
                    || false,
                )
                .unwrap();
        } else {
            inlet
                .step_molecular_mixture_bridge(
                    &mut mixture,
                    &bridge,
                    &mut frame,
                    Time::new(2e-6),
                    mixing::MolecularStepControl::default(),
                    DRAWS,
                    |time, _, _, _, densities| mixture_mean(time, densities),
                    || false,
                )
                .unwrap();
        }
        assert_eq!(mixture.flow().time().value().to_bits(), end.to_bits());
        assert_eq!(frame.time_s().to_bits(), end.to_bits());
    }
    assert_eq!(source.step(), 6);
}
