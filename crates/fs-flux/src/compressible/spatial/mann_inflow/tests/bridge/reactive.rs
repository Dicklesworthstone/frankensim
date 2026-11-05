use super::*;
use crate::compressible::mixture::tests::{gas as make_gas, source as species};
use crate::compressible::spatial::mixture::{self as gas, reactive as chemistry};
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};
use fs_gas::chemistry::IsochoricKineticsMethod;
use fs_qty::Temperature;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as GAS_R,
};

const FORCE: BodyForce3d<'static> = BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]);
const FORCE_CONTROL: ForceControl3d = ForceControl3d { maximum_mach: 5.0 };

fn reactive_mean(
    time: f64,
    densities: &mut [f64],
) -> Result<MixtureReservoirFaceState, MixtureError> {
    densities.copy_from_slice(&[0.8 + time, 0.2 - time]);
    Ok(MixtureReservoirFaceState {
        velocity: mean(time).velocity,
        pressure_pa: 1.0,
    })
}

fn reacting(inlet: &MannReservoirInflow3, frame: &MannOuBridgeFrame) -> chemistry::ReactingEuler3d {
    let gas = make_gas(vec![
        species("A", GAS_R, 3.5, 0.0, 20.0),
        species("B", GAS_R, 3.5, 0.0, 12.0),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("mann-reactive-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(50.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = chemistry::IdealGasKineticsV1::new(
        gas.source().clone(),
        "manufactured-mann-reactive-inlet",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let grid = grid();
    let mut states = vec![0.0; grid.len() * 6];
    for row in states.chunks_exact_mut(6) {
        gas::state_into(&gas, &[0.8, 0.2], [0.0; 3], 1.0, row).unwrap();
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
            |_, _, _, densities| reactive_mean(frame.time_s(), densities),
            || false,
        )
        .unwrap();
    let controls = chemistry::IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::Sdirk2,
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-11,
        relative_tolerance: 2e-6,
        absolute_amount_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 0.001,
        maximum_step_s: 1.0,
        maximum_steps: 10_000,
        ..chemistry::IsochoricKineticsSolve::default()
    };
    let reactor =
        chemistry::CellChemistry::new(flow.gas().clone(), mechanism, controls, || false).unwrap();
    let bytes = chemistry::ReactingEuler3d::required_workspace_bytes(&flow, &reactor).unwrap();
    chemistry::ReactingEuler3d::with_chemistry(flow, reactor, bytes, || false).unwrap()
}

fn operator(flow: &MixtureEuler3d) -> mixing::MolecularTransport3d<ConstantTransport> {
    let transport =
        mixing::SpatialTransport::new(flow.gas().clone(), law(flow.gas().source().clone()), || {
            false
        })
        .unwrap();
    let bytes = mixing::MolecularTransport3d::required_workspace_bytes(flow, &transport).unwrap();
    mixing::MolecularTransport3d::new(flow, transport, bytes, || false).unwrap()
}

fn balanced(report: &chemistry::ReactiveStepReport3d, forced: bool) {
    assert!(report.chemistry.species_change_kg[0] < 0.0);
    assert!(report.chemistry.species_change_kg[1] > 0.0);
    assert!(report.chemistry.referenced_energy_change_j > 0.0);
    assert_eq!(report.source_change.work != 0.0, forced);
    assert!(report.mass_residual_kg.abs() < 3e-11);
    assert!(report.absolute_energy_residual_j.abs() < 3e-10);
    for value in report.conservation_residual.species_masses {
        assert!(value.abs() < 3e-11);
    }
    for value in report.conservation_residual.momentum {
        assert!(value.abs() < 3e-11);
    }
    assert!(report.conservation_residual.energy.abs() < 3e-10);
}

#[test]
fn reacting_bridge_matches_direct_native_chemistry_with_and_without_force() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    for forced in [false, true] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut direct_frame = bridge.frame(10_000, || false).unwrap();
        let mut flow = reacting(&inlet, &frame);
        let mut native = reacting(&inlet, &frame);
        let mut workspace = flow.flow().forcing_workspace(100_000).unwrap();
        let mut native_workspace = native.flow().forcing_workspace(100_000).unwrap();
        let duration = Time::new(bridge.dt_s() / 4.0);
        let report = inlet
            .step_reacting_bridge(
                &mut flow,
                &bridge,
                &mut frame,
                duration,
                chemistry::ReactiveStepControl::default(),
                forced.then_some((FORCE, FORCE_CONTROL, &mut workspace)),
                DRAWS,
                |time, _, _, _, densities| reactive_mean(time, densities),
                || false,
            )
            .unwrap();
        let sample = |time, side, _, p: [f64; 3], densities: &mut [f64]| {
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
            let mut state = reactive_mean(time, densities)?;
            for (v, u) in state.velocity.iter_mut().zip([-u[1], u[0], u[2]]) {
                *v += u;
            }
            Ok(state)
        };
        let expected = if forced {
            native.step_with_forcing_and_reservoir_profiles(
                duration,
                chemistry::ReactiveStepControl::default(),
                FORCE,
                FORCE_CONTROL,
                &mut native_workspace,
                1000,
                sample,
                || false,
            )
        } else {
            native.step_with_reservoir_profiles(
                duration,
                chemistry::ReactiveStepControl::default(),
                1000,
                sample,
                || false,
            )
        }
        .unwrap();
        assert_eq!(report.flow, expected);
        assert_eq!(flow.flow().states(), native.flow().states());
        assert_eq!(flow.flow().time(), duration);
        assert_eq!(frame.time_s(), duration.value());
        assert_eq!(report.frame_updates, 2);
        balanced(&report.flow, forced);
        assert_eq!(source.checkpoint(), checkpoint);
    }
}

#[allow(clippy::too_many_arguments)]
fn molecular_step(
    inlet: &MannReservoirInflow3,
    bridge: &MannOuBridge,
    frame: &mut MannOuBridgeFrame,
    flow: &mut chemistry::ReactingEuler3d,
    op: &mut mixing::MolecularTransport3d<ConstantTransport>,
    workspace: &mut gas::ForcingWorkspace3d,
    duration: Time,
    stabilized: bool,
    forced: bool,
    draws: usize,
) -> Result<(chemistry::ReactiveStepReport3d, usize), MannReservoirError3> {
    let bytes = flow.required_stabilized_workspace_bytes(op).unwrap() + workspace.workspace_bytes();
    let forcing = forced.then_some((FORCE, FORCE_CONTROL, workspace));
    if stabilized {
        inlet
            .step_reacting_with_stabilized_dissipation_bridge(
                flow,
                op,
                bridge,
                frame,
                duration,
                chemistry::ReactiveStepControl::default(),
                mixing::StabilizedDissipationControl::default(),
                bytes,
                forcing,
                draws,
                |time, _, _, _, densities| reactive_mean(time, densities),
                || false,
            )
            .map(|r| (r.flow.coupled, r.normal_draws))
    } else {
        inlet
            .step_reacting_with_dissipation_bridge(
                flow,
                op,
                bridge,
                frame,
                duration,
                chemistry::ReactiveStepControl::default(),
                mixing::DissipationControl::default(),
                bytes,
                forcing,
                draws,
                |time, _, _, _, densities| reactive_mean(time, densities),
                || false,
            )
            .map(|r| (r.flow, r.normal_draws))
    }
}

fn molecular_loads(op: &mixing::MolecularTransport3d<ConstantTransport>) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            op.boundary_impulses(side).map(|q| {
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

#[test]
fn reacting_molecular_bridge_rolls_back_late_draw_failure_and_replays_both_owners() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let duration = Time::new(bridge.dt_s() / 4.0);
    for stabilized in [false, true] {
        for forced in [false, true] {
            let mut frame = bridge.frame(10_000, || false).unwrap();
            let mut reference_frame = bridge.frame(10_000, || false).unwrap();
            let mut flow = reacting(&inlet, &frame);
            let mut reference = reacting(&inlet, &frame);
            let mut op = operator(flow.flow());
            let mut reference_op = operator(reference.flow());
            let mut workspace = flow.flow().forcing_workspace(100_000).unwrap();
            let mut reference_workspace = reference.flow().forcing_workspace(100_000).unwrap();
            let expected = molecular_step(
                &inlet,
                &bridge,
                &mut reference_frame,
                &mut reference,
                &mut reference_op,
                &mut reference_workspace,
                duration,
                stabilized,
                forced,
                DRAWS,
            )
            .unwrap();
            assert!(expected.1 > 0);
            let states = flow.flow().states().to_vec();
            let profile = flow.flow().reservoir_face(BoundarySide::XMin, 0).unwrap();
            let densities = profile.partial_densities.to_vec();
            let velocity = profile.velocity;
            let loads = molecular_loads(&op);
            let failure = molecular_step(
                &inlet,
                &bridge,
                &mut frame,
                &mut flow,
                &mut op,
                &mut workspace,
                duration,
                stabilized,
                forced,
                expected.1 - 1,
            );
            assert!(matches!(
                failure,
                Err(MannReservoirError3::BridgeWork { .. })
            ));
            assert_eq!(flow.flow().states(), states);
            assert_eq!(flow.flow().time().value(), 0.0);
            let profile = flow.flow().reservoir_face(BoundarySide::XMin, 0).unwrap();
            assert_eq!(profile.partial_densities, densities);
            assert_eq!(profile.velocity, velocity);
            assert!(flow.last_step().is_none());
            assert!(op.last_step().is_none());
            assert_eq!(molecular_loads(&op), loads);
            let actual = molecular_step(
                &inlet,
                &bridge,
                &mut frame,
                &mut flow,
                &mut op,
                &mut workspace,
                duration,
                stabilized,
                forced,
                DRAWS,
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(flow.flow().states(), reference.flow().states());
            assert_eq!(molecular_loads(&op), molecular_loads(&reference_op));
            assert_eq!(op.last_step(), reference_op.last_step());
            assert_eq!(frame.time_s(), flow.flow().time().value());
            let split = actual.0.dissipation.unwrap();
            assert_eq!(split.first.elapsed.value(), 0.5 * duration.value());
            assert_eq!(split.second.elapsed, split.first.elapsed);
            assert!(split.first.boundary_change.energy != 0.0);
            balanced(&actual.0, forced);
            assert_eq!(source.checkpoint(), checkpoint);
        }
    }
}

#[test]
fn reacting_split_stops_at_exact_rounded_bridge_endpoint() {
    let mut source = FixedStepMannField::from_snapshot(field(), 211, 1e-6, 10_000).unwrap();
    source.advance_to(6, 1000, || false).unwrap();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let start = bridge.start_time_s();
    let end = bridge.end_time_s();
    let duration = end - start;
    assert!((start + duration / 2.0) + duration / 2.0 > end);
    for stabilized in [false, true] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut flow = reacting(&inlet, &frame);
        let mut op = operator(flow.flow());
        let mut workspace = flow.flow().forcing_workspace(100_000).unwrap();
        let report = molecular_step(
            &inlet,
            &bridge,
            &mut frame,
            &mut flow,
            &mut op,
            &mut workspace,
            Time::new(2e-6),
            stabilized,
            true,
            DRAWS,
        )
        .unwrap()
        .0;
        assert_eq!(flow.flow().time().value(), end);
        assert_eq!(frame.time_s(), end);
        assert_eq!(report.elapsed.value(), duration);
        assert_eq!(
            flow.flow()
                .reservoir_face(BoundarySide::XMin, 0)
                .unwrap()
                .partial_densities,
            &[0.8 + end, 0.2 - end]
        );
        balanced(&report, true);
        assert_eq!(source.checkpoint(), checkpoint);
    }
}
