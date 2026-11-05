use super::*;
use crate::compressible::spatial::mixture::{
    ReservoirFaceState,
    molecular::{
        MolecularEuler3d, MolecularStepControl, StabilizedDissipationControl,
        StabilizedMolecularStepControl,
    },
};
use fs_gas::multicomponent::dissipation::Rkl2Scheme;

#[path = "forced_reservoir_stages.rs"]
mod forced;

const SIDES: [BoundarySide; 2] = [BoundarySide::XMin, BoundarySide::XMax];
const TOPOLOGY: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};

fn input(rho: &mut [f64], species: f64, pressure_pa: f64) -> ReservoirFaceState {
    rho.copy_from_slice(&[species, 1.0 - species, 0.0]);
    ReservoirFaceState {
        velocity: [0.0; 3],
        pressure_pa,
    }
}

fn prepared(
    grid: CartesianGrid,
    gas: MixtureGas,
    species: f64,
    diffusion: f64,
    heat: f64,
) -> (MixtureEuler3d, MolecularTransport3d<ConstantTransport>) {
    let mut flow = field(gas, grid, TOPOLOGY, |_| {
        ([species, 1.0 - species, 0.0], [0.0; 3], 1.0)
    });
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&SIDES)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &SIDES,
        bytes,
        |_, _, _, rho| Ok(input(rho, species, 1.0)),
        || false,
    )
    .unwrap();
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), diffusion, 0.0, 0.0, heat),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
    let molecular = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
    (flow, molecular)
}

fn coupled() -> MolecularEuler3d<ConstantTransport> {
    let (flow, molecular) = prepared(
        CartesianGrid::new(2, 2, 1, 0.4, 0.5, 0.6).unwrap(),
        standard(0.2),
        0.6,
        0.02,
        0.03,
    );
    let bytes = MolecularEuler3d::required_stabilized_workspace_bytes(&flow, &molecular).unwrap();
    MolecularEuler3d::new(flow, molecular, bytes, || false).unwrap()
}

fn balance(report: &MolecularReport) {
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 2e-11);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 2e-11);
    }
    close(report.conservation_residual.energy, 0.0, 2e-11);
    close(report.absolute_energy_residual_j, 0.0, 2e-9);
}

fn split_balance(report: &crate::compressible::spatial::mixture::molecular::MolecularStepReport3d) {
    balance(&report.first);
    balance(&report.second);
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 2e-11);
    }
    close(report.conservation_residual.energy, 0.0, 2e-11);
    close(report.absolute_energy_residual_j, 0.0, 2e-9);
}

#[test]
fn g0_timed_molecular_constant_profiles_match_held_original_operators() {
    for stabilized in [false, true] {
        let make = || {
            prepared(
                CartesianGrid::new(2, 2, 1, 0.4, 0.5, 0.6).unwrap(),
                standard(0.2),
                0.6,
                0.02,
                0.03,
            )
        };
        let (mut held, mut old) = make();
        let (mut timed, mut new) = make();
        for flow in [&mut held, &mut timed] {
            flow.update_reservoir_profiles(|_, _, _, rho| Ok(input(rho, 0.4, 1.2)), || false)
                .unwrap();
        }
        let controls = DissipationControl {
            maximum_step_s: 0.002,
            ..DissipationControl::default()
        };
        let mut calls = Vec::new();
        let mut sample = |time, side, face, point, rho: &mut [f64]| {
            assert!(rho.iter().all(|value| value.is_nan()));
            calls.push((time, side, face, point));
            Ok(input(rho, 0.4, 1.2))
        };
        let actual = if stabilized {
            let control = StabilizedDissipationControl {
                transport: controls,
                ..StabilizedDissipationControl::default()
            };
            let reference = old
                .advance_stabilized(&mut held, Time::new(0.005), control, || false)
                .unwrap();
            let result = new
                .advance_stabilized_with_reservoir_profiles(
                    &mut timed,
                    Time::new(4.0),
                    Time::new(0.005),
                    control,
                    200,
                    &mut sample,
                    || false,
                )
                .unwrap();
            assert_eq!(result, reference);
            result.transport
        } else {
            let reference = old
                .advance(&mut held, Time::new(0.005), controls, || false)
                .unwrap();
            let result = new
                .advance_with_reservoir_profiles(
                    &mut timed,
                    Time::new(4.0),
                    Time::new(0.005),
                    controls,
                    200,
                    &mut sample,
                    || false,
                )
                .unwrap();
            assert_eq!(result, reference);
            result
        };
        assert_eq!(timed.states(), held.states());
        assert_eq!(timed.time().value(), 0.0);
        assert_eq!(timed.last_step(), held.last_step());
        assert_eq!(timed.workspace_bytes(), held.workspace_bytes());
        assert_eq!(new.workspace_bytes(), old.workspace_bytes());
        assert!(calls.iter().all(|(t, _, _, _)| *t >= 4.0 && *t <= 4.005));
        assert_eq!(calls.last().unwrap().0, 4.005);
        for (time, side, face, point) in calls {
            assert_eq!(point, timed.grid().face_center(side, face).unwrap());
            assert!(time.is_finite());
        }
        balance(&actual);
    }
}

fn smooth_error(steps: usize, heat: bool, stabilized: bool) -> f64 {
    let grid = CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let initial = if heat { 1.0 } else { 0.4 };
    let (mut flow, mut molecular) = prepared(
        grid,
        equal_mass_gas(),
        initial,
        0.2,
        if heat { 0.5 } else { 0.0 },
    );
    let dt = 1.0 / steps as f64;
    let controls = DissipationControl {
        maximum_step_s: dt,
        ..DissipationControl::default()
    };
    for step in 0..steps {
        let sample = |absolute: f64, _, _, _, rho: &mut [f64]| {
            let wave = 0.2 * (2.0 * (absolute - 3.0)).sin();
            Ok(input(
                rho,
                if heat { 1.0 } else { 0.4 + wave },
                if heat { 1.0 + wave } else { 1.0 },
            ))
        };
        let report = if stabilized {
            molecular
                .advance_stabilized_with_reservoir_profiles(
                    &mut flow,
                    Time::new(3.0 + step as f64 * dt),
                    Time::new(dt),
                    StabilizedDissipationControl {
                        transport: controls,
                        ..StabilizedDissipationControl::default()
                    },
                    32,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            molecular
                .advance_with_reservoir_profiles(
                    &mut flow,
                    Time::new(3.0 + step as f64 * dt),
                    Time::new(dt),
                    controls,
                    32,
                    sample,
                    || false,
                )
                .unwrap()
        };
        balance(&report);
        if !heat {
            // Different original species calorics and enthalpy reference gauges
            // must retain isothermal diffusion and its conservative energy flow.
            close(
                report.boundary_change.energy,
                -0.999 * report.boundary_change.species_masses[0],
                3e-11,
            );
        }
    }
    // One cell, two half-cell Dirichlet faces: z'=k(b(t)-z), k=0.8.
    // For heat rho*Cv=2.5 and kappa=0.5; for species D=0.2.
    let k: f64 = 0.8;
    let omega: f64 = 2.0;
    let exact = initial
        + 0.2 * k / (k * k + omega * omega)
            * (k * omega.sin() - omega * omega.cos() + omega * (-k).exp());
    let cell = flow.cell(0).unwrap();
    let primitive = cell.primitive(flow.gas(), || false).unwrap();
    assert_eq!(cell.partial_densities[2], 0.0);
    let actual = if heat {
        primitive.temperature
    } else {
        close(primitive.temperature, 1.0, 2e-10);
        cell.partial_densities[0]
    };
    (actual - exact).abs()
}

#[test]
fn g1_changing_species_and_heat_baths_converge_at_second_order() {
    for stabilized in [false, true] {
        for heat in [false, true] {
            let errors = [20, 40, 80].map(|steps| smooth_error(steps, heat, stabilized));
            eprintln!("timed molecular stabilized={stabilized} heat={heat}: {errors:?}");
            assert!(errors[0] / errors[1] > 3.7);
            assert!(errors[1] / errors[2] > 3.7);
            assert!(errors[2] < 2e-5);
        }
    }
}

#[test]
fn g1_rkl_profiles_follow_legendre_stage_abscissae_and_endpoint() {
    let (mut flow, mut molecular) = prepared(
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        equal_mass_gas(),
        0.4,
        0.2,
        0.0,
    );
    let mut times = Vec::new();
    let report = molecular
        .advance_stabilized_with_reservoir_profiles(
            &mut flow,
            Time::new(2.0),
            Time::new(0.5),
            StabilizedDissipationControl::default(),
            100,
            |time, side, _, _, rho| {
                if side == BoundarySide::XMin {
                    times.push(time);
                }
                Ok(input(rho, 0.4 + 0.02 * (time - 2.0), 1.0))
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    let scheme = Rkl2Scheme::new(report.maximum_stages_used).unwrap();
    assert_eq!(times.len(), scheme.stage_count() + 1);
    for (stage, time) in times.into_iter().enumerate() {
        close(time, 2.0 + 0.5 * scheme.abscissa(stage).unwrap(), 1e-14);
    }
    for side in SIDES {
        close(
            flow.reservoir_face(side, 0).unwrap().partial_densities[0],
            0.41,
            1e-15,
        );
    }
    balance(&report.transport);
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    cells: Vec<f64>,
    profiles: Vec<f64>,
    loads: Vec<f64>,
    time: f64,
    flow_report: Option<StepReport>,
    molecular_report: Option<MolecularReport>,
    report: Option<crate::compressible::spatial::mixture::molecular::MolecularStepReport3d>,
    flow_bytes: usize,
}

fn snapshot(owner: &MolecularEuler3d<ConstantTransport>) -> Snapshot {
    let flow = owner.flow();
    let mut profiles = Vec::new();
    let mut loads = Vec::new();
    for side in BoundarySide::ALL {
        for face in 0..flow.grid().face_count(side) {
            if let Some(profile) = flow.reservoir_face(side, face) {
                profiles.extend_from_slice(profile.partial_densities);
                profiles.extend_from_slice(&profile.velocity);
                profiles.push(profile.pressure_pa);
            }
        }
        for load in flow.boundary_impulses(side) {
            loads.extend_from_slice(load.species_masses);
            loads.extend_from_slice(&load.momentum);
            loads.push(load.energy);
        }
        for load in owner.molecular().boundary_impulses(side) {
            loads.extend_from_slice(load.species_masses);
            loads.extend_from_slice(&load.momentum);
            loads.push(load.energy);
        }
        for load in owner.boundary_impulses(side) {
            loads.extend_from_slice(load.species_masses);
            loads.extend_from_slice(&load.momentum);
            loads.push(load.energy);
        }
    }
    Snapshot {
        cells: flow.states().to_vec(),
        profiles,
        loads,
        time: flow.time().value(),
        flow_report: flow.last_step().copied(),
        molecular_report: owner.molecular().last_step().copied(),
        report: owner.last_step().copied(),
        flow_bytes: flow.workspace_bytes(),
    }
}

#[test]
fn g1_coupled_profiles_use_physical_strang_intervals_and_one_budget() {
    for stabilized in [false, true] {
        let mut owner = coupled();
        // Establish a nonzero accepted clock, then probe absolute stage times.
        owner
            .step(Time::new(0.002), MolecularStepControl::default(), || false)
            .unwrap();
        let start = owner.flow().time().value();
        let dt = 0.001;
        let mut times = Vec::new();
        let sample = |time, side, face, _, rho: &mut [f64]| {
            if side == BoundarySide::XMin && face == 0 {
                times.push(time);
            }
            Ok(input(rho, 0.6 + time, 1.0 + 2.0 * time))
        };
        let report = if stabilized {
            owner
                .step_stabilized_with_reservoir_profiles(
                    Time::new(dt),
                    StabilizedMolecularStepControl::default(),
                    36,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            owner
                .step_with_reservoir_profiles(
                    Time::new(dt),
                    MolecularStepControl::default(),
                    28,
                    sample,
                    || false,
                )
                .unwrap()
        };
        assert_eq!(report.elapsed.value(), dt);
        assert_eq!(report.coupling_retries, 0);
        let offsets: &[f64] = if stabilized {
            &[0.0, 0.0, 1.0 / 6.0, 0.5, 0.0, 1.0, 0.5, 2.0 / 3.0, 1.0]
        } else {
            &[0.0, 0.0, 0.5, 0.0, 1.0, 0.5, 1.0]
        };
        assert_eq!(times.len(), offsets.len());
        for (time, offset) in times.iter().zip(offsets) {
            close(*time, start + offset * dt, 1e-15);
        }
        for side in SIDES {
            for face in 0..owner.flow().grid().face_count(side) {
                let profile = owner.flow().reservoir_face(side, face).unwrap();
                close(profile.partial_densities[0], 0.6 + start + dt, 1e-15);
                close(profile.pressure_pa, 1.0 + 2.0 * (start + dt), 1e-15);
            }
        }
        split_balance(&report);
    }
}

#[test]
fn g1_split_profile_endpoint_matches_accepted_clock_bits_at_rounded_half_origin() {
    for stabilized in [false, true] {
        let mut owner = coupled();
        owner
            .step(Time::new(6e-6), MolecularStepControl::default(), || false)
            .unwrap();
        let start = owner.flow().time().value();
        let end: f64 = 7e-6;
        let dt = end - start;
        assert_eq!((start + dt).to_bits(), end.to_bits());
        assert!((start + 0.5 * dt) + 0.5 * dt > end);
        let mut last: f64 = start;
        let sample = |time: f64, _, _, _, rho: &mut [f64]| {
            if time < start || time > end {
                return Err(Error::Invalid(
                    "prescribed history queried outside its physical interval",
                ));
            }
            last = time;
            let mut profile = input(rho, 0.6, 1.0);
            profile.velocity[1] = time;
            Ok(profile)
        };
        let report = if stabilized {
            owner
                .step_stabilized_with_reservoir_profiles(
                    Time::new(dt),
                    StabilizedMolecularStepControl::default(),
                    1000,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            owner
                .step_with_reservoir_profiles(
                    Time::new(dt),
                    MolecularStepControl::default(),
                    1000,
                    sample,
                    || false,
                )
                .unwrap()
        };
        assert_eq!(report.elapsed.value().to_bits(), dt.to_bits());
        assert_eq!(owner.flow().time().value().to_bits(), end.to_bits());
        assert_eq!(last.to_bits(), end.to_bits());
        for side in SIDES {
            for face in 0..owner.flow().grid().face_count(side) {
                assert_eq!(
                    owner.flow().reservoir_face(side, face).unwrap().velocity[1].to_bits(),
                    end.to_bits()
                );
            }
        }
        split_balance(&report);
    }
}

#[test]
fn g4_late_profile_budget_provider_and_source_failures_are_terminal_and_atomic() {
    for failure in 0..3 {
        let mut owner = coupled();
        let before = snapshot(&owner);
        let mut calls = 0;
        let budget = if failure == 0 { 24 } else { 1000 };
        let error = owner
            .step_with_reservoir_profiles(
                Time::new(0.001),
                MolecularStepControl::default(),
                budget,
                |time, _, _, _, rho| {
                    calls += 1;
                    let result = input(rho, 0.6 + time, 1.0 + time);
                    if calls == 25 {
                        if failure == 1 {
                            return Err(Error::Invalid("late prescribed bath refusal"));
                        }
                        if failure == 2 {
                            rho[0] = f64::NAN;
                        }
                    }
                    Ok(result)
                },
                || false,
            )
            .unwrap_err();
        assert!(matches!(error, DissipationError::Flow(_)), "{error:?}");
        assert_eq!(calls, if failure == 0 { 24 } else { 25 });
        assert_eq!(snapshot(&owner), before);
        // Private trial profiles had advanced through the first half and Euler.
        // A following held call must restore the accepted bank before use.
        let mut reference = coupled();
        let actual = owner
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let expected = reference
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(snapshot(&owner), snapshot(&reference));
    }
}

#[test]
fn g4_outer_retries_keep_one_profile_budget_and_replay_shortened_times() {
    let mut reference = coupled();
    let mut calls = 0;
    let history = |time: f64, rho: &mut [f64]| {
        let mut profile = input(rho, 0.6, 1.0);
        profile.velocity[0] = 1e6 * time;
        profile
    };
    let report = reference
        .step_with_reservoir_profiles(
            Time::new(0.001),
            MolecularStepControl::default(),
            10000,
            |time, _, _, _, rho| {
                calls += 1;
                Ok(history(time, rho))
            },
            || false,
        )
        .unwrap();
    assert!(report.coupling_retries > 0);
    assert!(report.elapsed.value() < 0.001);
    split_balance(&report);
    let mut limited = coupled();
    let original = snapshot(&limited);
    let mut limited_calls = 0;
    let error = limited
        .step_with_reservoir_profiles(
            Time::new(0.001),
            MolecularStepControl::default(),
            calls - 1,
            |time, _, _, _, rho| {
                limited_calls += 1;
                Ok(history(time, rho))
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        DissipationError::Flow(Error::Invalid(
            "mixture reservoir profile face evaluation budget"
        ))
    ));
    assert_eq!(limited_calls, calls - 1);
    assert_eq!(snapshot(&limited), original);
}

#[test]
fn g4_final_cancellation_preserves_timed_split_cells_profiles_reports_and_loads() {
    for stabilized in [false, true] {
        let run = |owner: &mut MolecularEuler3d<ConstantTransport>,
                   cancelled: &mut dyn FnMut() -> bool| {
            let sample = |time, _, _, _, rho: &mut [f64]| Ok(input(rho, 0.6 + time, 1.0 + time));
            if stabilized {
                owner
                    .step_stabilized_with_reservoir_profiles(
                        Time::new(0.001),
                        StabilizedMolecularStepControl::default(),
                        1000,
                        sample,
                        cancelled,
                    )
                    .map(|report| report.transport)
            } else {
                owner.step_with_reservoir_profiles(
                    Time::new(0.001),
                    MolecularStepControl::default(),
                    1000,
                    sample,
                    cancelled,
                )
            }
        };
        let mut reference = coupled();
        let mut polls = 0;
        let expected = run(&mut reference, &mut || {
            polls += 1;
            false
        })
        .unwrap();
        let mut owner = coupled();
        let before = snapshot(&owner);
        let mut observed = 0;
        assert!(matches!(
            run(&mut owner, &mut || {
                observed += 1;
                observed == polls
            }),
            Err(DissipationError::Cancelled)
        ));
        assert_eq!(snapshot(&owner), before);
        let actual = run(&mut owner, &mut || false).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(snapshot(&owner), snapshot(&reference));
    }
}
