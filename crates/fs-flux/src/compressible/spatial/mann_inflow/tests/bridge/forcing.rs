use super::*;
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};
use crate::compressible::tests::close;

fn cap() -> ForceControl3d {
    ForceControl3d { maximum_mach: 10.0 }
}

fn owned_mixture_loads<'a>(
    loads: impl Iterator<Item = crate::compressible::spatial::mixture::BoundaryImpulse<'a>>,
) -> Vec<(Vec<f64>, [f64; 3], f64)> {
    loads
        .map(|load| (load.species_masses.to_vec(), load.momentum, load.energy))
        .collect()
}

fn forced_coupled_mixture(flow: MixtureEuler3d) -> mixing::MolecularEuler3d<ConstantTransport> {
    let (flow, molecular) = coupled_mixture(flow).into_parts();
    // The complete coupled cap includes externally supplied forcing storage
    // as well as the original SSP/RKL field, face and molecular scratch.
    let bytes = mixing::MolecularEuler3d::required_stabilized_workspace_bytes(&flow, &molecular)
        .unwrap()
        .checked_add(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    mixing::MolecularEuler3d::new(flow, molecular, bytes, || false).unwrap()
}

// Independent realization observation in this fixture's authored quarter-turn.
fn fluctuation(
    bridge: &MannOuBridge,
    frame: &mut MannOuBridgeFrame,
    inlet: &MannReservoirInflow3,
    time: f64,
    side: BoundarySide,
    point: [f64; 3],
) -> [f64; 3] {
    bridge.sample_frame(time, frame, DRAWS, || false).unwrap();
    let u = frame
        .sample_filtered(
            [1.0 + point[1], 2.0 - point[0], 3.0 + point[2]],
            inlet.spatial_filter(side).unwrap(),
        )
        .unwrap()
        .u;
    [-u[1], u[0], u[2]]
}

#[test]
fn forced_frozen_and_species_bridges_match_direct_native_force_and_realization() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let duration = Time::new(bridge.dt_s() / 4.0);
    for force in [
        BodyForce3d::UniformAcceleration([0.7, -0.2, 0.1]),
        BodyForce3d::UniformForceDensity([0.9, -0.3, 0.2]),
    ] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut direct = bridge.frame(10_000, || false).unwrap();
        let mut flow = prepared(&inlet, &frame);
        let mut native = prepared(&inlet, &frame);
        let mut work = flow.forcing_workspace(100_000).unwrap();
        let mut native_work = native.forcing_workspace(100_000).unwrap();
        let report = inlet
            .step_forced_frozen_bridge(
                &mut flow,
                &bridge,
                &mut frame,
                duration,
                StepControl::default(),
                force,
                cap(),
                &mut work,
                DRAWS,
                |time, _, _, _| Ok(mean(time)),
                || false,
            )
            .unwrap();
        let expected = native
            .step_with_forcing_and_reservoir_profiles(
                duration,
                StepControl::default(),
                force,
                cap(),
                &mut native_work,
                1000,
                |time, side, _, point| {
                    let u = fluctuation(&bridge, &mut direct, &inlet, time, side, point);
                    let mut state = mean(time);
                    for (v, u) in state.velocity.iter_mut().zip(u) {
                        *v += u;
                    }
                    Ok(state)
                },
                || false,
            )
            .unwrap();
        assert_eq!(report.flow, expected);
        assert_eq!(flow.states(), native.states());
        assert_eq!(flow.time(), duration);
        assert_eq!(frame.time_s(), flow.time().value());
        assert_eq!(report.evaluated_faces, 2 * inlet.faces);
        assert_ne!(report.flow.source_change.work, 0.0);
        for residual in report.flow.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
        for side in BoundarySide::ALL {
            assert_eq!(flow.boundary_impulses(side), native.boundary_impulses(side));
        }

        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut direct = bridge.frame(10_000, || false).unwrap();
        let mut flow = prepared_mixture(&inlet, &frame);
        let mut native = prepared_mixture(&inlet, &frame);
        let mut work = flow.forcing_workspace(100_000).unwrap();
        let mut native_work = native.forcing_workspace(100_000).unwrap();
        let report = inlet
            .step_forced_mixture_bridge(
                &mut flow,
                &bridge,
                &mut frame,
                duration,
                StepControl::default(),
                force,
                cap(),
                &mut work,
                DRAWS,
                |time, _, _, _, densities| mixture_mean(time, densities),
                || false,
            )
            .unwrap();
        let expected = native
            .step_with_forcing_and_reservoir_profiles(
                duration,
                StepControl::default(),
                force,
                cap(),
                &mut native_work,
                1000,
                |time, side, _, point, densities| {
                    let u = fluctuation(&bridge, &mut direct, &inlet, time, side, point);
                    let mut state = mixture_mean(time, densities)?;
                    for (v, u) in state.velocity.iter_mut().zip(u) {
                        *v += u;
                    }
                    Ok(state)
                },
                || false,
            )
            .unwrap();
        assert_eq!(report.flow, expected);
        assert_eq!(flow.states(), native.states());
        assert_eq!(flow.time(), duration);
        assert_eq!(frame.time_s(), flow.time().value());
        assert_eq!(report.evaluated_faces, 2 * inlet.faces);
        assert_ne!(report.flow.source_change.work, 0.0);
        close(report.flow.absolute_energy_residual_j, 0.0, 3e-10);
        for side in BoundarySide::ALL {
            assert_eq!(
                owned_mixture_loads(flow.boundary_impulses(side)),
                owned_mixture_loads(native.boundary_impulses(side))
            );
        }
    }
    assert_eq!(source.checkpoint(), checkpoint);
}

#[test]
fn forced_frozen_molecular_bridges_share_native_ssp_rkl_and_roll_back_late_draws() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let duration = Time::new(bridge.dt_s() / 4.0);
    let force = BodyForce3d::UniformForceDensity([0.9, -0.3, 0.2]);
    for stabilized in [false, true] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut direct = bridge.frame(10_000, || false).unwrap();
        let mut flow = coupled(prepared(&inlet, &frame));
        let mut native = coupled(prepared(&inlet, &frame));
        let mut work = flow.flow().forcing_workspace(100_000).unwrap();
        let mut native_work = native.flow().forcing_workspace(100_000).unwrap();
        let step = |flow: &mut viscous::MolecularEuler3d<_>,
                    frame: &mut MannOuBridgeFrame,
                    work: &mut crate::compressible::spatial::ForcingWorkspace3d,
                    budget| {
            if stabilized {
                inlet
                    .step_stabilized_forced_molecular_frozen_bridge(
                        flow,
                        &bridge,
                        frame,
                        duration,
                        viscous::MolecularStepControl::default(),
                        viscous::MolecularStabilization::default(),
                        force,
                        cap(),
                        work,
                        budget,
                        |time, _, _, _| Ok(mean(time)),
                        || false,
                    )
                    .map(|r| (r.flow.transport, r.normal_draws))
            } else {
                inlet
                    .step_forced_molecular_frozen_bridge(
                        flow,
                        &bridge,
                        frame,
                        duration,
                        viscous::MolecularStepControl::default(),
                        force,
                        cap(),
                        work,
                        budget,
                        |time, _, _, _| Ok(mean(time)),
                        || false,
                    )
                    .map(|r| (r.flow, r.normal_draws))
            }
        };
        let (report, draws) = step(&mut flow, &mut frame, &mut work, DRAWS).unwrap();
        let sample = |time, side, _, point| {
            let u = fluctuation(&bridge, &mut direct, &inlet, time, side, point);
            let mut state = mean(time);
            for (v, u) in state.velocity.iter_mut().zip(u) {
                *v += u;
            }
            Ok(state)
        };
        let expected = if stabilized {
            native
                .step_stabilized_with_forcing_and_reservoir_profiles(
                    duration,
                    viscous::MolecularStepControl::default(),
                    viscous::MolecularStabilization::default(),
                    force,
                    cap(),
                    &mut native_work,
                    10000,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            native
                .step_with_forcing_and_reservoir_profiles(
                    duration,
                    viscous::MolecularStepControl::default(),
                    force,
                    cap(),
                    &mut native_work,
                    10000,
                    sample,
                    || false,
                )
                .unwrap()
        };
        assert_eq!(report, expected);
        assert_eq!(flow.flow().states(), native.flow().states());
        assert_eq!(frame.time_s(), flow.flow().time().value());
        assert_eq!(report.source_change, report.advection.source_change);
        let volume = flow.flow().grid().cell_volume() * flow.flow().states().len() as f64;
        for (axis, density_force) in [0.9, -0.3, 0.2].into_iter().enumerate() {
            close(
                report.source_change.impulse[axis],
                volume * density_force * report.elapsed.value(),
                3e-13,
            );
        }
        assert_ne!(report.source_change.work, 0.0);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
        let mut refused_frame = bridge.frame(10_000, || false).unwrap();
        let mut refused = coupled(prepared(&inlet, &refused_frame));
        let mut refused_work = refused.flow().forcing_workspace(100_000).unwrap();
        let original = refused.flow().states().to_vec();
        let profile = refused.flow().reservoir_face(BoundarySide::XMin, 0);
        let loads = BoundarySide::ALL.map(|s| refused.boundary_impulses(s).to_vec());
        assert!(draws > 0);
        assert!(matches!(
            step(
                &mut refused,
                &mut refused_frame,
                &mut refused_work,
                draws - 1
            ),
            Err(MannReservoirError3::BridgeWork { .. })
        ));
        assert_eq!(refused.flow().states(), original);
        assert_eq!(refused.flow().time().value(), 0.0);
        assert_eq!(
            refused.flow().reservoir_face(BoundarySide::XMin, 0),
            profile
        );
        assert!(refused.last_step().is_none());
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(refused.boundary_impulses(side), loads[i]);
        }
        let (replayed, _) =
            step(&mut refused, &mut refused_frame, &mut refused_work, DRAWS).unwrap();
        assert_eq!(replayed, report);
        assert_eq!(refused.flow().states(), flow.flow().states());
        assert_eq!(source.checkpoint(), checkpoint);
    }
}

#[test]
fn forced_species_molecular_bridges_keep_species_enthalpy_and_late_draw_transaction() {
    let source = source();
    let checkpoint = source.checkpoint();
    let bridge = bridge(&source);
    let inlet = inlet(&source);
    let duration = Time::new(bridge.dt_s() / 4.0);
    let force = BodyForce3d::UniformAcceleration([0.7, -0.2, 0.1]);
    for stabilized in [false, true] {
        let mut frame = bridge.frame(10_000, || false).unwrap();
        let mut direct = bridge.frame(10_000, || false).unwrap();
        let mut flow = forced_coupled_mixture(prepared_mixture(&inlet, &frame));
        let mut native = forced_coupled_mixture(prepared_mixture(&inlet, &frame));
        let mut work = flow.flow().forcing_workspace(100_000).unwrap();
        let mut native_work = native.flow().forcing_workspace(100_000).unwrap();
        let step = |flow: &mut mixing::MolecularEuler3d<_>,
                    frame: &mut MannOuBridgeFrame,
                    work: &mut crate::compressible::spatial::mixture::ForcingWorkspace3d,
                    budget| {
            if stabilized {
                inlet
                    .step_stabilized_forced_molecular_mixture_bridge(
                        flow,
                        &bridge,
                        frame,
                        duration,
                        mixing::StabilizedMolecularStepControl::default(),
                        force,
                        cap(),
                        work,
                        budget,
                        |time, _, _, _, densities| mixture_mean(time, densities),
                        || false,
                    )
                    .map(|r| (r.flow.transport, r.normal_draws))
            } else {
                inlet
                    .step_forced_molecular_mixture_bridge(
                        flow,
                        &bridge,
                        frame,
                        duration,
                        mixing::MolecularStepControl::default(),
                        force,
                        cap(),
                        work,
                        budget,
                        |time, _, _, _, densities| mixture_mean(time, densities),
                        || false,
                    )
                    .map(|r| (r.flow, r.normal_draws))
            }
        };
        let (report, draws) = step(&mut flow, &mut frame, &mut work, DRAWS).unwrap();
        let sample = |time, side, _, point, densities: &mut [f64]| {
            let u = fluctuation(&bridge, &mut direct, &inlet, time, side, point);
            let mut state = mixture_mean(time, densities)?;
            for (v, u) in state.velocity.iter_mut().zip(u) {
                *v += u;
            }
            Ok(state)
        };
        let expected = if stabilized {
            native
                .step_stabilized_with_forcing_and_reservoir_profiles(
                    duration,
                    mixing::StabilizedMolecularStepControl::default(),
                    force,
                    cap(),
                    &mut native_work,
                    10000,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            native
                .step_with_forcing_and_reservoir_profiles(
                    duration,
                    mixing::MolecularStepControl::default(),
                    force,
                    cap(),
                    &mut native_work,
                    10000,
                    sample,
                    || false,
                )
                .unwrap()
        };
        assert_eq!(report, expected);
        assert_eq!(flow.flow().states(), native.flow().states());
        assert_eq!(frame.time_s(), flow.flow().time().value());
        assert_eq!(report.source_change, report.advection.source_change);
        assert_ne!(report.source_change.work, 0.0);
        for &residual in &report.conservation_residual.species_masses {
            close(residual, 0.0, 3e-12);
        }
        for residual in report.conservation_residual.momentum {
            close(residual, 0.0, 3e-12);
        }
        close(report.conservation_residual.energy, 0.0, 3e-12);
        close(report.absolute_energy_residual_j, 0.0, 3e-10);
        let mut refused_frame = bridge.frame(10_000, || false).unwrap();
        let mut refused = forced_coupled_mixture(prepared_mixture(&inlet, &refused_frame));
        let mut refused_work = refused.flow().forcing_workspace(100_000).unwrap();
        let original = refused.flow().states().to_vec();
        let profile = refused
            .flow()
            .reservoir_face(BoundarySide::XMin, 0)
            .unwrap()
            .partial_densities
            .to_vec();
        let loads = BoundarySide::ALL.map(|s| owned_mixture_loads(refused.boundary_impulses(s)));
        assert!(draws > 0);
        assert!(matches!(
            step(
                &mut refused,
                &mut refused_frame,
                &mut refused_work,
                draws - 1
            ),
            Err(MannReservoirError3::BridgeWork { .. })
        ));
        assert_eq!(refused.flow().states(), original);
        assert_eq!(refused.flow().time().value(), 0.0);
        assert_eq!(
            refused
                .flow()
                .reservoir_face(BoundarySide::XMin, 0)
                .unwrap()
                .partial_densities,
            profile
        );
        assert!(refused.last_step().is_none());
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(
                owned_mixture_loads(refused.boundary_impulses(side)),
                loads[i]
            );
        }
        let (replayed, _) =
            step(&mut refused, &mut refused_frame, &mut refused_work, DRAWS).unwrap();
        assert_eq!(replayed, report);
        assert_eq!(refused.flow().states(), flow.flow().states());
        assert_eq!(source.checkpoint(), checkpoint);
    }
}
