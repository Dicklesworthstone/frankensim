//! G1/G3/G4 native loaded-solid, conduction and irreversible heat regressions.
//! All coefficients are manufactured; these are not measured-life predictions.
use std::collections::BTreeMap;

use fs_couple::solid_enthalpy::creep::{
    CreepCell, CreepConfig, CreepCouplingError, CreepSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::source::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::*;

const L: f64 = 0.01;
const V: f64 = L * L * L / 6.0;
const E: f64 = 1.2e9;
const RHO: f64 = 1000.0;
const STRESS: f64 = 2e6;
const RATE: f64 = 2e-4;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const NODES: [[f64; 3]; 4] = [[0.0; 3], [L, 0.0, 0.0], [0.0, L, 0.0], [0.0, 0.0, L]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ROLLERS: [(usize, f64); 9] = [
    (0, 0.0),
    (1, 0.0),
    (2, 0.0),
    (4, 0.0),
    (5, 0.0),
    (6, 0.0),
    (8, 0.0),
    (9, 0.0),
    (10, 0.0),
];

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn constant(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, value), (380.0, value)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            280.0,
            380.0,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured spatial Norton fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    cp: f64,
    rate: f64,
    q: f64,
}
impl Fixture {
    fn new(cp: f64, rate: f64, q: f64) -> Self {
        let mut claims = ClaimSet::new();
        for claim in [
            constant(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
            constant(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
            constant(CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
            constant(CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, rate),
            constant(CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
            constant(
                CREEP_ACTIVATION_ENERGY_PROPERTY,
                CREEP_ACTIVATION_ENERGY_DIMS,
                q,
            ),
            constant(
                CREEP_REFERENCE_TEMPERATURE_PROPERTY,
                Temperature::DIMS,
                300.0,
            ),
            constant("density", Dims([-3, 1, 0, 0, 0, 0]), RHO),
            constant("specific_heat_capacity", CP_DIMS, cp),
        ] {
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic Norton solid".into(),
                phase: "solid".into(),
                process: "manufactured".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
                density: &PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(380.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * cp,
            1e-5,
        )
        .unwrap();
        Self {
            card,
            caloric,
            cp,
            rate,
            q,
        }
    }
    fn source(&self) -> NortonCreepSource<'_> {
        NortonCreepSource::new(
            &self.card,
            NortonSourceQuery {
                lower: &point(280.0),
                upper: &point(380.0),
                temperature_axis: "T",
                selection: MaterialPropertySelection::SingleClaimOnly,
                strain_limit: 0.05,
            },
            || false,
        )
        .unwrap()
    }
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate(&gate, f)
}
fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 0x43524854,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn config() -> CreepConfig {
    CreepConfig {
        maximum_iterations: 100,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-9,
        stress_tolerance_pa: 1e-3,
        enthalpy_relaxation: 0.7,
    }
}
fn heat_config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-11,
    }
}
fn initial<'a>(
    fixture: &'a Fixture,
    temperatures: &[f64],
    volume: f64,
    links: &[(usize, usize, f64)],
    reservoir: bool,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, &t) in temperatures.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: RHO * volume,
                specific_h_j_kg: (t - 300.0) * fixture.cp,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L; 3],
                    volume_m3: volume,
                    surface_area_m2: 4.0 * L * L,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: reservoir.then_some(Reservoir {
                        area_m2: 0.1 * L * L,
                        conductance_w_k: 0.0,
                        temperature_k: 350.0,
                    }),
                },
            }),
        });
    }
    for (i, &(a, b, g)) in links.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(a as u64),
                b: CellId(b as u64),
                area_m2: 0.1 * L * L,
                conductance_w_k: g,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 16,
            max_links: 32,
            max_events: 64,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, heat_config(0.01), || Ok(()))
    .unwrap()
    .0
}
fn bind<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    cx: &Cx<'_>,
) -> CreepSolid<'a> {
    CreepSolid::try_new(
        thermal,
        nodes,
        tets,
        (0..tets.len())
            .map(|i| CreepCell {
                cell: CellId(i as u64),
                card: &fixture.card,
                source: fixture.source(),
                caloric: &fixture.caloric,
            })
            .collect(),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "{actual} != {expected}; tol={tol}"
    );
}
// Independent stress-control/heat equation, with no production return-map call.
fn reference(f: &Fixture, initial_t: f64, dt: f64) -> (f64, f64) {
    let rate = |t: f64| f.rate * (f.q / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / t)).exp();
    let mut lo = initial_t;
    let mut hi = initial_t + 5.0;
    assert!(hi - initial_t - STRESS * dt * rate(hi) / (RHO * f.cp) > 0.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if t - initial_t - STRESS * dt * rate(t) / (RHO * f.cp) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, dt * rate(t))
}

#[test]
fn g1_source_to_loaded_tetrahedron_creeps_heats_and_continues_once_per_physical_step() {
    with_cx(|cx| {
        for multiplier in [1.0, 2.0] {
            let f = Fixture::new(5.0, RATE * multiplier, 30_000.0);
            let mut endpoints = Vec::new();
            for heat_dt in [0.5, 0.025] {
                let mut thermal = initial(&f, &[300.0], V, &[], false);
                let mut solid = bind(&f, &thermal, &NODES, &TETS, cx);
                let reference_claims: Vec<_> = solid.materials()[0]
                    .material()
                    .properties()
                    .iter()
                    .map(|p| p.answer().receipt.selected)
                    .collect();
                let mut temperature = 300.0;
                let mut alpha = 0.0;
                for n in 1..=4 {
                    let (expected_t, delta) = reference(&f, temperature, 0.5);
                    alpha += delta;
                    let candidate = solid
                        .advance(
                            &thermal,
                            &(),
                            n as f64 * 0.5,
                            &ROLLERS,
                            &[(3, STRESS * V / L)],
                            IncrementSettings::default(),
                            heat_config(heat_dt),
                            config(),
                            cx,
                        )
                        .unwrap();
                    let solution = candidate.solid().solution();
                    close(solution.element_stress[0][0], STRESS, 1e-2);
                    assert!(
                        solution.element_stress[0][1..]
                            .iter()
                            .all(|s| s.abs() < 1e-2)
                    );
                    close(
                        solution.element_states[0].equivalent_creep_strain,
                        alpha,
                        2e-10,
                    );
                    close(solution.element_strain[0][0], STRESS / E + alpha, 2e-10);
                    close(
                        solution.displacement_m[1][0],
                        L * (STRESS / E + alpha),
                        2e-12,
                    );
                    close(
                        candidate.thermal().temperature_k(CellId(0)).unwrap(),
                        expected_t,
                        2e-7,
                    );
                    let report = candidate.report();
                    close(
                        report.mechanical_energy.physical_dissipation_j,
                        STRESS * V * delta,
                        1e-10,
                    );
                    close(
                        report.thermal_audit.stored_energy_change_j,
                        STRESS * V * delta,
                        1e-10,
                    );
                    close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                    if n == 1 {
                        assert!(
                            report.mechanical_energy.numerical_dissipation_j
                                > report.applied_heat_j
                        );
                    }
                    let material = candidate.solid().materials()[0].material();
                    close(material.state_coordinate("T").unwrap(), expected_t, 2e-7);
                    assert_eq!(
                        material
                            .properties()
                            .iter()
                            .map(|p| p.answer().receipt.selected)
                            .collect::<Vec<_>>(),
                        reference_claims
                    );
                    assert_eq!(solid.time_s(), (n - 1) as f64 * 0.5);
                    (thermal, solid, _, _) = candidate.into_parts();
                    temperature = expected_t;
                }
                endpoints.push((
                    solid.solution().element_states[0].equivalent_creep_strain,
                    thermal.temperature_k(CellId(0)).unwrap(),
                ));
            }
            close(endpoints[0].0, endpoints[1].0, 2e-10);
            close(endpoints[0].1, endpoints[1].1, 2e-7);
        }
    });
}

fn bar_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    (
        vec![
            [0.0, 0.0, 0.0],
            [2.0 * L, 0.0, 0.0],
            [0.0, L, 0.0],
            [2.0 * L, L, 0.0],
            [0.0, 0.0, L],
            [2.0 * L, 0.0, L],
            [0.0, L, L],
            [2.0 * L, L, L],
        ],
        vec![
            [0, 1, 3, 7],
            [0, 1, 7, 5],
            [0, 2, 7, 3],
            [0, 2, 6, 7],
            [0, 4, 5, 7],
            [0, 4, 7, 6],
        ],
    )
}

#[test]
fn g3_nonuniform_conduction_changes_native_loaded_bar_shape_and_local_creep_history() {
    let (nodes, tets) = bar_mesh();
    let supports: Vec<_> = nodes
        .iter()
        .enumerate()
        .flat_map(|(i, x)| (0..3).filter_map(move |d| (x[d] == 0.0).then_some((3 * i + d, 0.0))))
        .collect();
    let force = STRESS * L * L;
    let loads = [
        (3, force / 3.0),
        (9, force / 6.0),
        (15, force / 6.0),
        (21, force / 3.0),
    ];
    with_cx(|cx| {
        let f = Fixture::new(100.0, RATE, 30_000.0);
        let mut outcomes = Vec::new();
        for conductance in [0.0, 0.1] {
            // Adjacent tetrahedra share these faces in the original cube mesh.
            let links = [
                (0, 1, conductance),
                (0, 2, conductance),
                (2, 3, conductance),
                (3, 5, conductance),
                (5, 4, conductance),
                (4, 1, conductance),
            ];
            let thermal = initial(
                &f,
                &[300.0, 300.0, 320.0, 320.0, 340.0, 340.0],
                2.0 * V,
                &links,
                false,
            );
            let solid = bind(&f, &thermal, &nodes, &tets, cx);
            let candidate = solid
                .advance(
                    &thermal,
                    &(),
                    0.2,
                    &supports,
                    &loads,
                    IncrementSettings::default(),
                    heat_config(0.002),
                    config(),
                    cx,
                )
                .unwrap();
            let temperatures: Vec<_> = (0..6)
                .map(|i| candidate.thermal().temperature_k(CellId(i)).unwrap())
                .collect();
            let histories: Vec<_> = candidate
                .solid()
                .solution()
                .element_states
                .iter()
                .map(|s| s.equivalent_creep_strain)
                .collect();
            assert!(histories.iter().all(|x| *x > 0.0));
            assert!(
                histories.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                    - histories.iter().copied().fold(f64::INFINITY, f64::min)
                    > 1e-6
            );
            close(
                candidate.report().thermal_audit.stored_energy_change_j,
                candidate.report().mechanical_energy.physical_dissipation_j,
                1e-9,
            );
            close(
                candidate.report().combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            for i in 0..6 {
                close(
                    candidate.solid().materials()[i]
                        .material()
                        .state_coordinate("T")
                        .unwrap(),
                    temperatures[i],
                    1e-10,
                );
            }
            outcomes.push((
                temperatures,
                candidate.solid().solution().displacement_m.clone(),
                histories,
            ));
        }
        assert!(outcomes[1].0[0] > outcomes[0].0[0]);
        assert!(outcomes[1].0[4] < outcomes[0].0[4]);
        let shape_difference: f64 = outcomes[0]
            .1
            .iter()
            .flatten()
            .zip(outcomes[1].1.iter().flatten())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(shape_difference > 1e-9);
        assert_ne!(outcomes[0].2, outcomes[1].2);
    });
}

#[derive(Clone, Debug, PartialEq)]
struct HeatObserver {
    steps: usize,
    interval_energy_j: f64,
    reservoir_calls: usize,
    conductance: f64,
    interval_power: f64,
    fail_after: Option<usize>,
}
impl DepositionObserver for HeatObserver {
    type Error = DepositionError;
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        dt: f64,
        _: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        self.interval_energy_j += dt * self.interval_power;
        Ok(BTreeMap::from([(CellId(0), self.interval_power)]))
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        self.reservoir_calls += 1;
        Ok(BTreeMap::from([(CellId(0), self.conductance)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        self.steps += 1;
        if self.fail_after.is_some_and(|n| self.steps >= n) {
            Err(DepositionError::InvalidInput(
                "manufactured late observer failure",
            ))
        } else {
            Ok(())
        }
    }
}
fn observer() -> HeatObserver {
    HeatObserver {
        steps: 0,
        interval_energy_j: 0.0,
        reservoir_calls: 0,
        conductance: 0.05,
        interval_power: 0.01,
        fail_after: None,
    }
}

#[test]
fn g1_finite_boundary_heat_and_interval_observer_reach_the_actual_loaded_creep_trajectory() {
    with_cx(|cx| {
        let f = Fixture::new(100.0, RATE, 0.0);
        let thermal = initial(&f, &[300.0], V, &[], true);
        let solid = bind(&f, &thermal, &NODES, &TETS, cx);
        let observer = observer();
        let result = solid
            .advance(
                &thermal,
                &observer,
                0.2,
                &ROLLERS,
                &[(3, STRESS * V / L)],
                IncrementSettings::default(),
                heat_config(0.01),
                config(),
                cx,
            )
            .unwrap();
        let mut expected = 300.0;
        for _ in 0..20 {
            expected += 0.01
                * (STRESS * V * RATE
                    + observer.interval_power
                    + observer.conductance * (350.0 - expected))
                / (RHO * V * f.cp);
        }
        close(
            result.thermal().temperature_k(CellId(0)).unwrap(),
            expected,
            1e-7,
        );
        close(
            result.solid().solution().element_states[0].equivalent_creep_strain,
            RATE * 0.2,
            1e-10,
        );
        close(
            result.observer().interval_energy_j,
            observer.interval_power * 0.2,
            1e-14,
        );
        assert!(result.observer().reservoir_calls >= 20);
        assert_eq!(observer.steps, 0);
        close(
            result.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert!(result.report().external_heat_j > result.report().applied_heat_j);
    });
}

#[test]
fn g4_feedback_observer_source_and_cancellation_refusals_preserve_all_accepted_state() {
    let gate = CancelGate::new();
    with_gate(&gate, |cx| {
        let f = Fixture::new(5.0, RATE, 30_000.0);
        let thermal = initial(&f, &[300.0], V, &[], true);
        let solid = bind(&f, &thermal, &NODES, &TETS, cx);
        let original_cells = thermal.cells().clone();
        let original_solution = solid.solution().clone();
        let mut observer = observer();
        observer.fail_after = Some(2);
        let original_observer = observer.clone();
        let forces = [(3, STRESS * V / L)];
        assert!(matches!(
            solid.advance(
                &thermal,
                &observer,
                0.02,
                &ROLLERS,
                &forces,
                IncrementSettings::default(),
                heat_config(0.001),
                config(),
                cx
            ),
            Err(CreepCouplingError::Observer(_))
        ));
        let exhausted = CreepConfig {
            maximum_iterations: 1,
            ..config()
        };
        assert!(matches!(
            solid.advance(
                &thermal,
                &(),
                0.1,
                &ROLLERS,
                &forces,
                IncrementSettings::default(),
                heat_config(0.01),
                exhausted,
                cx
            ),
            Err(CreepCouplingError::NotConverged { .. })
        ));
        let foreign = Fixture::new(5.0, 2.0 * RATE, 30_000.0);
        let foreign_thermal = initial(&foreign, &[300.0], V, &[], true);
        assert!(matches!(
            solid.advance(
                &foreign_thermal,
                &(),
                0.02,
                &ROLLERS,
                &forces,
                IncrementSettings::default(),
                heat_config(0.001),
                config(),
                cx
            ),
            Err(CreepCouplingError::Invalid(_))
        ));
        let retried = solid
            .advance(
                &thermal,
                &(),
                0.1,
                &ROLLERS,
                &forces,
                IncrementSettings::default(),
                heat_config(0.01),
                config(),
                cx,
            )
            .unwrap();
        assert!(retried.solid().solution().element_states[0].equivalent_creep_strain > 0.0);
        gate.request();
        assert!(matches!(
            solid.advance(
                &thermal,
                &(),
                0.1,
                &ROLLERS,
                &forces,
                IncrementSettings::default(),
                heat_config(0.01),
                config(),
                cx
            ),
            Err(CreepCouplingError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(solid.solution(), &original_solution);
        assert_eq!(observer, original_observer);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
    });
}
