use super::*;
use crate::compressible::spatial::mixture::molecular::{
    MolecularWall3d, WallAcceleration3d, WallThermalCondition,
};

const MU: f64 = 0.75;
const CAP: ForceControl3d = ForceControl3d { maximum_mach: 5.0 };

struct Case {
    flow: ReactingEuler3d,
    molecular: MolecularTransport3d<ConstantTransport>,
    scratch: ForcingWorkspace3d,
}

impl Case {
    fn bytes(&self) -> usize {
        budget(&self.flow, &self.molecular) + self.scratch.workspace_bytes()
    }

    fn driven(
        &mut self,
        duration: f64,
        acceleration: f64,
        body: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let bytes = self.bytes();
        self.flow
            .step_with_dissipation_forcing_and_wall_acceleration(
                &mut self.molecular,
                Time::new(duration),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                BodyForce3d::UniformAcceleration(body),
                CAP,
                &mut self.scratch,
                &drive(acceleration),
                cancelled,
            )
    }
}

fn drive(acceleration: f64) -> [WallAcceleration3d; 1] {
    [WallAcceleration3d {
        side: BoundarySide::YMax,
        acceleration: [acceleration, 0.0, 0.0],
    }]
}

fn wall_speed(operator: &MolecularTransport3d<ConstantTransport>) -> f64 {
    match operator.wall(BoundarySide::YMax) {
        MolecularWall3d::NoSlip { velocity, .. } => velocity[0],
        _ => panic!("expected no-slip wall"),
    }
}

fn case(heat: f64, rate: f64, activation: f64) -> Case {
    case_with_limit(heat, rate, activation, 50.0)
}

fn case_with_limit(heat: f64, rate: f64, activation: f64, maximum_temperature: f64) -> Case {
    let (gas, mut mechanism) = mechanism(heat, rate, 37.0);
    if activation != 0.0 {
        let reaction = ElementaryReactionV1::new(
            ReactionId::new("activated-wall-driven-isomerization").unwrap(),
            vec![(SpeciesId::new("A").unwrap(), 1)],
            vec![(SpeciesId::new("B").unwrap(), 1)],
            ArrheniusRateV1::new(
                Temperature::new(1.0),
                ReactionProgressRateQuantityV1::new(rate),
                0.0,
                MolarEnergyQuantityV1::new(activation * R),
            )
            .unwrap(),
        )
        .unwrap();
        mechanism = IdealGasKineticsV1::new(
            gas.source().clone(),
            "manufactured-wall-driven-reacting-gas",
            1,
            vec![reaction],
            1e-12,
        )
        .unwrap();
    }
    let input = field(
        gas,
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([0.9, 0.1], [0.0; 3], 1.0),
    );
    let law = ConstantTransport::new(
        input.gas().source().clone(),
        "manufactured-reacting-Couette",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(MU),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
        },
        vec![0.0, 0.001, 0.001, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(input.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&input, &source).unwrap();
    let mut molecular = MolecularTransport3d::new(&input, source, bytes, || false).unwrap();
    for side in [BoundarySide::YMin, BoundarySide::YMax] {
        molecular = molecular
            .with_wall(
                side,
                MolecularWall3d::NoSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap();
    }
    let scratch = input.forcing_workspace(usize::MAX).unwrap();
    let controls = IsochoricKineticsSolve {
        maximum_temperature_k: maximum_temperature,
        ..chemistry()
    };
    let reactor = CellChemistry::new(input.gas().clone(), mechanism, controls, || false).unwrap();
    let bytes = ReactingEuler3d::required_workspace_bytes(&input, &reactor).unwrap();
    let flow = ReactingEuler3d::with_chemistry(input, reactor, bytes, || false).unwrap();
    Case {
        flow,
        molecular,
        scratch,
    }
}

fn check(report: &ReactiveStepReport3d) {
    balance(report);
    let molecular = report.dissipation.unwrap();
    assert_eq!(
        molecular.first.elapsed.value(),
        0.5 * report.elapsed.value()
    );
    assert_eq!(molecular.second.elapsed, molecular.first.elapsed);
    assert!(report.chemistry.species_change_kg[1].abs() > 0.0);
    assert_eq!(report.chemistry.cell_solves, 2);
}

#[test]
fn g1_accelerating_reacting_couette_matches_momentum_work_heat_and_species() {
    let end: f64 = 0.2;
    let acceleration: f64 = 2.0;
    let lambda = 4.0 * MU; // rho=H=A=1; both no-slip half-cell stresses participate.
    let exp = (-lambda * end).exp();
    let exact_velocity = 0.5 * acceleration * (end - (1.0 - exp) / lambda);
    let exact_work = MU
        * acceleration.powi(2)
        * (end.powi(3) / 3.0 + end * end / (2.0 * lambda)
            - (1.0 - (1.0 + lambda * end) * exp) / lambda.powi(3));
    let exact_heat = exact_work - 0.5 * exact_velocity * exact_velocity;
    let expected = [exact_velocity, exact_work, exact_heat];
    let mut errors = Vec::new();
    for dt in [0.01, 0.005, 0.0025] {
        let mut c = case(0.0, 2.0, 0.0);
        let bytes = c.bytes();
        let initial = c.flow.flow().totals().energy;
        let mut work = 0.0;
        while c.flow.flow().time().value() < end {
            let duration = (end - c.flow.flow().time().value()).min(dt);
            let report = c
                .flow
                .step_with_dissipation_and_wall_acceleration(
                    &mut c.molecular,
                    Time::new(duration),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    &drive(acceleration),
                    || false,
                )
                .unwrap();
            check(&report);
            assert_eq!(report.coupling_retries, 0);
            work += report.boundary_change.energy;
        }
        let q = c.flow.flow().cell(0).unwrap();
        let velocity = q.momentum[0];
        let heat = c.flow.flow().totals().energy - initial - 0.5 * velocity * velocity;
        close(c.flow.flow().totals().energy - initial, work, 2e-9);
        close(wall_speed(&c.molecular), acceleration * end, 2e-14);
        close(c.flow.flow().time().value(), end, 1e-14);
        // Equal calorics/references give an independent reversible A <-> B
        // solution while the real reactor and wall-driven temperature evolve.
        close(q.partial_densities[0], 0.5 + 0.4 * (-4.0 * end).exp(), 3e-6);
        close(q.partial_densities[0] + q.partial_densities[1], 1.0, 3e-12);
        assert_eq!(c.bytes(), bytes);
        assert!(work > 0.0 && heat > 0.0);
        errors.push(std::array::from_fn::<_, 3, _>(|i| {
            ([velocity, work, heat][i] - expected[i]).abs()
        }));
    }
    eprintln!("reacting accelerating Couette dt .01/.005/.0025: {errors:?}");
    for component in 0..3 {
        for pair in errors.windows(2) {
            let ratio = pair[0][component] / pair[1][component];
            assert!(ratio > 3.5 && ratio < 4.8, "{component}: {ratio}");
        }
        assert!(errors[2][component] < 8e-7);
    }
}

#[test]
fn g2_wall_work_changes_activated_chemistry_with_complete_body_and_reference_energy_ledgers() {
    let end: f64 = 0.15;
    let body = [0.25, 0.0, -0.15];
    let run = |acceleration: f64| {
        let mut c = case(4.0, 1.5, 3.0);
        let bytes = c.bytes();
        let before = c.flow.flow().totals();
        let mut boundary_work = 0.0;
        let mut body_work = 0.0;
        let mut chemical_change = 0.0;
        let mut body_impulse = [0.0; 3];
        while c.flow.flow().time().value() < end {
            let duration = (end - c.flow.flow().time().value()).min(0.01);
            let report = c.driven(duration, acceleration, body, || false).unwrap();
            check(&report);
            boundary_work += report.boundary_change.energy;
            body_work += report.source_change.work;
            chemical_change += report.chemistry.referenced_energy_change_j;
            for axis in 0..3 {
                body_impulse[axis] += report.source_change.impulse[axis];
            }
            let mut outward = [0.0; 6];
            for row in loads(&c.flow) {
                for d in 0..6 {
                    outward[d] += row[d];
                }
            }
            for axis in 0..3 {
                close(
                    report.after.momentum[axis] - report.before.momentum[axis] + outward[2 + axis]
                        - report.source_change.impulse[axis],
                    0.0,
                    3e-11,
                );
            }
            close(outward[5] + report.boundary_change.energy, 0.0, 3e-11);
        }
        let after = c.flow.flow().totals();
        close(
            after.energy - before.energy - chemical_change,
            boundary_work + body_work,
            3e-9,
        );
        close(
            after.absolute_energy(c.flow.flow().gas()).unwrap()
                - before.absolute_energy(c.flow.flow().gas()).unwrap(),
            boundary_work + body_work,
            3e-9,
        );
        for axis in 0..3 {
            close(body_impulse[axis], body[axis] * end, 3e-12);
        }
        close(wall_speed(&c.molecular), acceleration * end, 3e-13);
        assert!(chemical_change > 0.01 && body_work > 0.0);
        assert_eq!(c.bytes(), bytes);
        let q = c.flow.flow().cell(0).unwrap();
        let temperature = q
            .primitive(c.flow.flow().gas(), || false)
            .unwrap()
            .temperature;
        (q.partial_densities[1], temperature, boundary_work)
    };
    let held = run(0.0);
    let driven = run(30.0);
    eprintln!("activated exothermic held/driven (B,T,wall work): {held:?} / {driven:?}");
    assert!(driven.0 > held.0 + 1e-3);
    assert!(driven.1 > held.1 + 0.05);
    assert!(driven.2 > 0.1);
    close(held.2, 0.0, 1e-13);
}

#[test]
fn g3_heat_driven_shortening_and_final_rollback_replay_the_original_wall_history() {
    let mut c = case(40.0, 500.0, 0.0);
    let report = c.driven(0.2, 3.0, [0.0; 3], || false).unwrap();
    check(&report);
    assert!(report.coupling_retries > 0);
    close(
        wall_speed(&c.molecular),
        3.0 * report.elapsed.value(),
        2e-14,
    );
    let mut reference = case(40.0, 500.0, 0.0);
    let replay = reference
        .driven(report.elapsed.value(), 3.0, [0.0; 3], || false)
        .unwrap();
    assert_eq!(replay.coupling_retries, 0);
    assert_eq!(replay.transport, report.transport);
    assert_eq!(replay.chemistry, report.chemistry);
    assert_eq!(replay.dissipation, report.dissipation);
    assert_eq!(reference.flow.flow().states(), c.flow.flow().states());
    assert_eq!(loads(&reference.flow), loads(&c.flow));
    assert_eq!(
        molecular_loads(&reference.molecular),
        molecular_loads(&c.molecular)
    );
    assert_eq!(wall_speed(&reference.molecular), wall_speed(&c.molecular));
    let before = (
        c.flow.flow().states().to_vec(),
        c.flow.flow().time(),
        c.flow.last_step().copied(),
        c.flow.flow().last_step().copied(),
        c.molecular.last_step().copied(),
        loads(&c.flow),
        molecular_loads(&c.molecular),
        BoundarySide::ALL.map(|side| c.molecular.wall(side)),
        c.bytes(),
    );
    let mut polls = 0;
    let expected = reference
        .driven(0.001, -2.0, [0.2, 0.0, -0.1], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        c.driven(0.001, -2.0, [0.2, 0.0, -0.1], || {
            count += 1;
            count >= polls
        }),
        Err(ReactiveError::Cancelled)
    ));
    let bytes = c.bytes();
    assert!(
        c.flow
            .step_with_dissipation_forcing_and_wall_acceleration(
                &mut c.molecular,
                Time::new(0.001),
                ReactiveStepControl {
                    maximum_coupling_retries: 0,
                    ..ReactiveStepControl::default()
                },
                DissipationControl {
                    maximum_step_s: 1e-12,
                    maximum_steps: 1,
                    ..DissipationControl::default()
                },
                bytes,
                BodyForce3d::UniformAcceleration([0.2, 0.0, -0.1]),
                CAP,
                &mut c.scratch,
                &drive(-2.0),
                || false,
            )
            .is_err()
    );
    let after = (
        c.flow.flow().states().to_vec(),
        c.flow.flow().time(),
        c.flow.last_step().copied(),
        c.flow.flow().last_step().copied(),
        c.molecular.last_step().copied(),
        loads(&c.flow),
        molecular_loads(&c.molecular),
        BoundarySide::ALL.map(|side| c.molecular.wall(side)),
        c.bytes(),
    );
    assert_eq!(before, after);
    let actual = c.driven(0.001, -2.0, [0.2, 0.0, -0.1], || false).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(c.flow.flow().states(), reference.flow.flow().states());
    assert_eq!(loads(&c.flow), loads(&reference.flow));
    assert_eq!(
        molecular_loads(&c.molecular),
        molecular_loads(&reference.molecular)
    );
    assert_eq!(wall_speed(&c.molecular), wall_speed(&reference.molecular));
}

#[test]
fn g4_final_accelerating_wall_heat_cannot_escape_the_chemical_temperature_domain() {
    let mut c = case_with_limit(0.0, 2.0, 0.0, 1.04);
    let initial = c.flow.flow().states().to_vec();
    let initial_loads = loads(&c.flow);
    let initial_molecular = molecular_loads(&c.molecular);
    let bytes = c.bytes();
    // The first D half stays below 1.04 K and both C halves are admissible.
    // The second D half completes above it: the final chemical-domain guard
    // must refuse before publishing either the field or accelerated wall.
    let error = c
        .flow
        .step_with_dissipation_forcing_and_wall_acceleration(
            &mut c.molecular,
            Time::new(0.05),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl::default(),
            bytes,
            BodyForce3d::UniformAcceleration([0.0; 3]),
            CAP,
            &mut c.scratch,
            &drive(60.0),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, ReactiveError::RetryLimit { last, .. }
        if matches!(*last, ReactiveError::Invalid("final molecular state is outside chemistry interval"))));
    assert_eq!(c.flow.flow().states(), initial);
    assert_eq!(c.flow.flow().time().value(), 0.0);
    assert!(c.flow.last_step().is_none() && c.molecular.last_step().is_none());
    assert_eq!(loads(&c.flow), initial_loads);
    assert_eq!(molecular_loads(&c.molecular), initial_molecular);
    assert_eq!(wall_speed(&c.molecular), 0.0);
    let accepted = c.driven(0.005, 60.0, [0.0; 3], || false).unwrap();
    check(&accepted);
    close(wall_speed(&c.molecular), 0.3, 1e-14);
}

#[test]
fn g0_zero_wall_acceleration_recovers_both_held_reacting_paths_and_counts_force_scratch() {
    for forced in [false, true] {
        let mut held = case(4.0, 8.0, 0.0);
        let mut driven = case(4.0, 8.0, 0.0);
        let bytes = held.bytes();
        let body = BodyForce3d::UniformAcceleration([0.2, 0.0, -0.1]);
        let expected = if forced {
            held.flow
                .step_with_dissipation_and_forcing(
                    &mut held.molecular,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    body,
                    CAP,
                    &mut held.scratch,
                    || false,
                )
                .unwrap()
        } else {
            held.flow
                .step_with_dissipation(
                    &mut held.molecular,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    || false,
                )
                .unwrap()
        };
        let actual = if forced {
            driven
                .driven(0.002, -0.0, [0.2, 0.0, -0.1], || false)
                .unwrap()
        } else {
            driven
                .flow
                .step_with_dissipation_and_wall_acceleration(
                    &mut driven.molecular,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    &drive(-0.0),
                    || false,
                )
                .unwrap()
        };
        assert_eq!(actual, expected);
        for (a, b) in held
            .flow
            .flow()
            .states()
            .iter()
            .zip(driven.flow.flow().states())
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(loads(&held.flow), loads(&driven.flow));
        assert_eq!(
            molecular_loads(&held.molecular),
            molecular_loads(&driven.molecular)
        );
        assert_eq!(
            wall_speed(&held.molecular).to_bits(),
            wall_speed(&driven.molecular).to_bits()
        );
    }
    let mut c = case(4.0, 8.0, 0.0);
    let states = c.flow.flow().states().to_vec();
    let bytes = c.bytes();
    assert!(
        c.flow
            .step_with_dissipation_forcing_and_wall_acceleration(
                &mut c.molecular,
                Time::new(0.002),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes - 1,
                BodyForce3d::UniformAcceleration([0.0; 3]),
                CAP,
                &mut c.scratch,
                &drive(0.0),
                || false,
            )
            .is_err()
    );
    assert_eq!(c.flow.flow().states(), states);
    assert_eq!(wall_speed(&c.molecular), 0.0);
    assert!(c.flow.last_step().is_none() && c.molecular.last_step().is_none());
}
