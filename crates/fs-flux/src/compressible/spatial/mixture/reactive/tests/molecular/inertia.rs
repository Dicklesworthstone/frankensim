use super::*;
use crate::compressible::spatial::mixture::molecular::{
    InertialWall3d, MolecularWall3d, WallInertiaControl3d, WallThermalCondition,
};
use fs_qty::Mass;

const MASS: f64 = 0.4;
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
    fn advance(
        &mut self,
        dt: f64,
        body: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        let bytes = self.bytes();
        self.flow.step_with_dissipation_and_forcing(
            &mut self.molecular,
            Time::new(dt),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            BodyForce3d::UniformAcceleration(body),
            CAP,
            &mut self.scratch,
            cancelled,
        )
    }
}
fn case(heat: f64, rate: f64, activation: f64, velocity: [f64; 3]) -> Case {
    let (gas, mut mechanism) = mechanism(heat, rate, 37.0);
    if activation != 0.0 {
        let reaction = ElementaryReactionV1::new(
            ReactionId::new("finite-wall-activated-isomerization").unwrap(),
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
            "manufactured-finite-wall-reacting-gas",
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
        "finite-wall-reacting-transport",
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
    let bytes = MolecularTransport3d::required_workspace_bytes(&input, &source).unwrap()
        + WallInertiaControl3d::solver_workspace_bytes();
    let molecular = MolecularTransport3d::new(&input, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity,
                thermal: WallThermalCondition::Adiabatic,
            },
        )
        .unwrap()
        .with_wall_inertia(
            &[InertialWall3d {
                side: BoundarySide::YMax,
                mass: Mass::new(MASS),
            }],
            WallInertiaControl3d::default(),
        )
        .unwrap();
    let scratch = input.forcing_workspace(usize::MAX).unwrap();
    Case {
        flow: coupled(input, mechanism),
        molecular,
        scratch,
    }
}
fn wall(operator: &MolecularTransport3d<ConstantTransport>) -> [f64; 3] {
    match operator.wall(BoundarySide::YMax) {
        MolecularWall3d::NoSlip { velocity, .. } => velocity,
        _ => panic!("finite no-slip wall"),
    }
}
fn kinetic(mass: f64, velocity: [f64; 3]) -> f64 {
    0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>()
}
fn check(report: &ReactiveStepReport3d) {
    balance(report);
    let d = report.dissipation.unwrap();
    let a = *d
        .first
        .wall_inertia
        .unwrap()
        .wall(BoundarySide::YMax)
        .unwrap();
    let b = *d
        .second
        .wall_inertia
        .unwrap()
        .wall(BoundarySide::YMax)
        .unwrap();
    assert_eq!(a.velocity_after, b.velocity_before);
    for axis in [0, 2] {
        close(
            report.after.momentum[axis] - report.before.momentum[axis]
                + MASS * (b.velocity_after[axis] - a.velocity_before[axis])
                - report.source_change.impulse[axis],
            0.0,
            5e-11,
        );
    }
    close(
        report.after.energy - report.before.energy - report.chemistry.referenced_energy_change_j
            + a.kinetic_energy_change_j
            + b.kinetic_energy_change_j
            - report.source_change.work,
        0.0,
        3e-9,
    );
    close(a.thermal_exchange_j + b.thermal_exchange_j, 0.0, 5e-11);
    assert!(
        d.second.wall_inertia.unwrap().field_evaluations
            > d.first.wall_inertia.unwrap().field_evaluations
    );
}

#[test]
fn g1_reacting_finite_wall_matches_recoil_and_species_with_distinct_time_rungs() {
    let end: f64 = 0.1;
    let initial = [0.6, 0.0, -0.2];
    let decay = (-2.0 * MU * (1.0 + MASS.recip()) * end).exp();
    let expected = initial.map(|v| MASS / (1.0 + MASS) * v * (1.0 - decay));
    let expected_wall = initial.map(|v| v * (MASS + decay) / (1.0 + MASS));
    let expected_heat =
        kinetic(MASS, initial) - kinetic(MASS, expected_wall) - kinetic(1.0, expected);
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for step in [0.01, 0.005, 0.0025] {
        let mut c = case(0.0, 2.0, 0.0, initial);
        let bytes = c.bytes();
        let initial_energy = c.flow.flow().totals().energy;
        let mut count = 0;
        for index in 0..(end / step).round() as usize {
            let target = (((index + 1) as f64) * step).min(end);
            let dt = target - c.flow.flow().time().value();
            let report = c.advance(dt, [0.0; 3], || false).unwrap();
            check(&report);
            assert_eq!(report.coupling_retries, 0);
            count += 1;
        }
        let q = c.flow.flow().cell(0).unwrap();
        let p = q.primitive(c.flow.flow().gas(), || false).unwrap();
        let heat = c.flow.flow().totals().energy - initial_energy - kinetic(1.0, p.velocity);
        close(q.partial_densities[0], 0.5 + 0.4 * (-4.0 * end).exp(), 3e-6);
        close(q.partial_densities.iter().sum::<f64>(), 1.0, 2e-12);
        close(
            c.flow.flow().totals().energy + kinetic(MASS, wall(&c.molecular)),
            initial_energy + kinetic(MASS, initial),
            3e-9,
        );
        assert_eq!(c.bytes(), bytes);
        errors.push([
            (p.velocity[0] - expected[0]).abs(),
            (wall(&c.molecular)[0] - expected_wall[0]).abs(),
            (heat - expected_heat).abs(),
        ]);
        counts.push(count);
    }
    eprintln!("reacting finite-wall counts {counts:?}, errors {errors:?}");
    assert!(counts[1] >= 2 * counts[0] - 1 && counts[2] >= 2 * counts[1] - 1);
    for component in 0..3 {
        for pair in errors.windows(2) {
            let ratio = pair[0][component] / pair[1][component];
            assert!(ratio > 3.4 && ratio < 4.8, "component {component}: {ratio}");
        }
        assert!(errors[2][component] < 1e-6);
    }
}

#[test]
fn g2_finite_wall_heat_changes_activated_chemistry_without_losing_body_or_reference_energy() {
    let body = [0.25, 0.0, -0.15];
    let end: f64 = 0.15;
    let run = |speed: f64| {
        let initial = [speed, 0.0, -0.4 * speed];
        let mut c = case(4.0, 1.5, 3.0, initial);
        let before = c.flow.flow().totals();
        let mut work = 0.0;
        let mut chemical = 0.0;
        while c.flow.flow().time().value() < end {
            let dt = (end - c.flow.flow().time().value()).min(0.005);
            let report = c.advance(dt, body, || false).unwrap();
            check(&report);
            work += report.source_change.work;
            chemical += report.chemistry.referenced_energy_change_j;
        }
        let after = c.flow.flow().totals();
        let wall_change = kinetic(MASS, wall(&c.molecular)) - kinetic(MASS, initial);
        close(
            after.energy - before.energy - chemical + wall_change,
            work,
            3e-9,
        );
        close(
            after.absolute_energy(c.flow.flow().gas()).unwrap()
                - before.absolute_energy(c.flow.flow().gas()).unwrap()
                + wall_change,
            work,
            3e-9,
        );
        for axis in [0, 2] {
            close(
                after.momentum[axis] - before.momentum[axis]
                    + MASS * (wall(&c.molecular)[axis] - initial[axis]),
                body[axis] * end,
                4e-11,
            );
        }
        let q = c.flow.flow().cell(0).unwrap();
        let p = q.primitive(c.flow.flow().gas(), || false).unwrap();
        assert!(chemical > 0.0);
        (q.partial_densities[1], p.temperature, wall_change)
    };
    let resting = run(0.0);
    let recoiling = run(3.0);
    eprintln!("reacting finite-wall (B,T,carrier dK) resting={resting:?}, recoiling={recoiling:?}");
    assert!(recoiling.0 > resting.0 + 1e-4);
    assert!(recoiling.1 > resting.1 + 0.1);
    assert!(recoiling.2 < -0.1);
}

#[test]
fn g3_reacting_finite_wall_heat_shortening_and_final_cancellation_replay_atomically() {
    let initial = [0.3, 0.0, -0.1];
    let mut c = case(40.0, 500.0, 0.0, initial);
    let report = c.advance(0.2, [0.0; 3], || false).unwrap();
    check(&report);
    assert!(report.coupling_retries > 0);
    let mut reference = case(40.0, 500.0, 0.0, initial);
    let replay = reference
        .advance(report.elapsed.value(), [0.0; 3], || false)
        .unwrap();
    assert_eq!(replay.coupling_retries, 0);
    assert_eq!(replay.transport, report.transport);
    assert_eq!(replay.chemistry, report.chemistry);
    assert_eq!(c.flow.flow().states(), reference.flow.flow().states());
    assert_eq!(wall(&c.molecular), wall(&reference.molecular));
    assert_eq!(loads(&c.flow), loads(&reference.flow));
    assert_eq!(
        molecular_loads(&c.molecular),
        molecular_loads(&reference.molecular)
    );
    let d = report.dissipation.unwrap();
    let r = replay.dissipation.unwrap();
    for (actual, expected) in [(d.first, r.first), (d.second, r.second)] {
        assert_eq!(
            actual.wall_inertia.unwrap().wall(BoundarySide::YMax),
            expected.wall_inertia.unwrap().wall(BoundarySide::YMax)
        );
        assert!(
            actual.wall_inertia.unwrap().field_evaluations
                > expected.wall_inertia.unwrap().field_evaluations
        );
    }
    let before = (
        c.flow.flow().states().to_vec(),
        c.flow.flow().time(),
        c.flow.last_step().copied(),
        c.flow.flow().last_step().copied(),
        c.molecular.last_step().copied(),
        wall(&c.molecular),
        loads(&c.flow),
        molecular_loads(&c.molecular),
        c.bytes(),
    );
    let mut polls = 0;
    let expected = reference
        .advance(0.001, [0.2, 0.0, -0.1], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        c.advance(0.001, [0.2, 0.0, -0.1], || {
            count += 1;
            count >= polls
        }),
        Err(ReactiveError::Cancelled)
    ));
    let after = (
        c.flow.flow().states().to_vec(),
        c.flow.flow().time(),
        c.flow.last_step().copied(),
        c.flow.flow().last_step().copied(),
        c.molecular.last_step().copied(),
        wall(&c.molecular),
        loads(&c.flow),
        molecular_loads(&c.molecular),
        c.bytes(),
    );
    assert_eq!(before, after);
    let recovered = c.advance(0.001, [0.2, 0.0, -0.1], || false).unwrap();
    assert_eq!(expected, recovered);
    assert_eq!(c.flow.flow().states(), reference.flow.flow().states());
}
