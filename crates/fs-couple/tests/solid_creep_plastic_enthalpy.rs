//! Source-driven yielding, permanent creep and actual finite thermal feedback.
//! Coefficients are manufactured implementation references, not measured data.

use fs_couple::solid_enthalpy::creep::{
    CreepCell, CreepConfig, CreepCouplingError, ZeroExpansionCreepSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::source::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::{IncrementSettings, TetAssemblyBudget};
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;

const L: f64 = 0.01;
const V: f64 = L * L * L / 6.0;
const E: f64 = 1.2e9;
const RHO: f64 = 1000.0;
const CP: f64 = 20.0;
const STRESS: f64 = 2e6;
const RATE: f64 = 2e-4;
const HARDENING: f64 = 1e8;
const ACTIVATION: f64 = 30_000.0;
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
fn yield_stress(t: f64) -> f64 {
    2.1e6 - 1e4 * (t - 300.0)
}
fn claim(name: &str, dims: Dims, lower: f64, upper: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, lower), (380.0, upper)],
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
            source: "manufactured combined creep/plastic thermal fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new() -> Self {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
            (POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
            (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
            (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE),
            (CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
            (
                CREEP_ACTIVATION_ENERGY_PROPERTY,
                CREEP_ACTIVATION_ENERGY_DIMS,
                ACTIVATION,
            ),
            (
                CREEP_REFERENCE_TEMPERATURE_PROPERTY,
                Temperature::DIMS,
                300.0,
            ),
            (
                ISOTROPIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                HARDENING,
            ),
            ("density", Dims([-3, 1, 0, 0, 0, 0]), RHO),
            ("specific_heat_capacity", CP_DIMS, CP),
        ] {
            claims
                .insert_claim(claim(name, dims, value, value))
                .unwrap();
        }
        claims
            .insert_claim(claim(
                YIELD_STRESS_PROPERTY,
                Pressure::DIMS,
                yield_stress(280.0),
                yield_stress(380.0),
            ))
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic combined creep/plastic solid".into(),
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
            -20.0 * CP,
            1e-5,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self) -> NortonJ2Source<'_> {
        NortonJ2Source::new(
            &self.card,
            NortonSourceQuery {
                lower: &point(280.0),
                upper: &point(380.0),
                temperature_axis: "T",
                selection: MaterialPropertySelection::SingleClaimOnly,
                strain_limit: 0.05,
            },
            HardeningPropertyBasis::PlasticStrainModulus,
            || false,
        )
        .unwrap()
    }
    fn thermal(&self) -> DepositionSimulation<'_> {
        let cell = DepositedCell {
            id: CellId(0),
            segment_id: 0,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg: RHO * V,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [L; 3],
                volume_m3: V,
                surface_area_m2: 4.0 * L * L,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                source_w: 0.0,
                reservoir: None,
            },
        };
        DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }],
            &[(MaterialId(1), self.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 4,
                max_links: 4,
                max_events: 8,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(0.0, heat_config(), || Ok(()))
        .unwrap()
        .0
    }
    fn bind<'a>(
        &'a self,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> ZeroExpansionCreepSolid<'a, NortonJ2Source<'a>> {
        ZeroExpansionCreepSolid::try_new(
            thermal,
            &NODES,
            &TETS,
            vec![CreepCell {
                cell: CellId(0),
                card: &self.card,
                source: self.source(),
                caloric: &self.caloric,
            }],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap()
    }
}
fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 0x43525048,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
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
fn heat_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-11,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a}, expected={b}, tolerance={tolerance}"
    );
}
#[derive(Clone, Debug)]
struct Heating {
    watts: f64,
    steps: usize,
    fail: bool,
}
impl DepositionObserver for Heating {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(CellId(0), self.watts)]))
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.fail {
            Err(DepositionError::InvalidInput(
                "injected thermal observer failure",
            ))
        } else {
            Ok(())
        }
    }
}

// Independent stress-controlled first law: hardening storage and BE loss are
// excluded from caloric heat. No production return map supplies this reference.
fn reference(old_t: f64, old_plastic: f64, dt: f64, external_j: f64) -> (f64, f64, f64) {
    let increments = |t: f64| {
        let dc = dt * RATE * (ACTIVATION / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / t)).exp();
        let dp = ((STRESS - yield_stress(t)) / HARDENING - old_plastic).max(0.0);
        (dc, dp)
    };
    let residual = |t: f64| {
        let (dc, dp) = increments(t);
        RHO * CP * (t - old_t) - external_j / V - STRESS * dc - yield_stress(t) * dp
    };
    let (mut lo, mut hi) = (old_t, old_t + external_j / (RHO * V * CP) + 2.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if residual(t) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    let (dc, dp) = increments(t);
    (t, dc, dp)
}

#[test]
fn finite_heating_softens_yield_while_creep_and_both_heat_channels_continue() {
    with_gate(&CancelGate::new(), |cx| {
        let fixture = Fixture::new();
        let mut thermal = fixture.thermal();
        let mut solid = fixture.bind(&thermal, cx);
        let original_claims: Vec<_> = solid.materials()[0]
            .material()
            .properties()
            .iter()
            .map(|p| p.answer().receipt.selected)
            .collect();
        let (mut temperature, mut creep, mut plastic) = (300.0, 0.0, 0.0);
        for (i, external_rise) in [0.0, 40.0, 0.0].into_iter().enumerate() {
            let dt = 0.5;
            let external_j = RHO * V * CP * external_rise;
            let observer = Heating {
                watts: external_j / dt,
                steps: 0,
                fail: false,
            };
            let (expected_t, dc, dp) = reference(temperature, plastic, dt, external_j);
            creep += dc;
            plastic += dp;
            let old_solution = solid.solution().clone();
            let candidate = solid
                .advance(
                    &thermal,
                    &observer,
                    (i + 1) as f64 * dt,
                    &ROLLERS,
                    &[(3, STRESS * V / L)],
                    IncrementSettings::default(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(solid.solution(), &old_solution);
            assert_eq!(observer.steps, 0);
            assert!(candidate.observer().steps > 1);
            let solution = candidate.solid().solution();
            let history = &solution.element_states[0];
            close(history.creep.equivalent_creep_strain, creep, 2e-10);
            close(history.plastic.alpha, plastic, 2e-10);
            close(solution.element_stress[0][0], STRESS, 1e-2);
            close(
                solution.displacement_m[1][0],
                L * (STRESS / E + creep + plastic),
                3e-12,
            );
            close(
                candidate.thermal().temperature_k(CellId(0)).unwrap(),
                expected_t,
                2e-7,
            );
            if i == 0 {
                assert_eq!(history.plastic.alpha, 0.0);
            } else {
                assert!(history.plastic.alpha > 0.003);
            }
            let report = candidate.report();
            let energy = report.element_energy[0];
            close(energy.creep_dissipation_j_m3 * V, STRESS * dc * V, 2e-10);
            close(
                energy.plastic_dissipation_j_m3 * V,
                yield_stress(expected_t) * dp * V,
                2e-10,
            );
            close(
                report.thermal_audit.stored_energy_change_j,
                external_j + V * (STRESS * dc + yield_stress(expected_t) * dp),
                3e-10,
            );
            close(report.external_heat_j, external_j, 2e-10);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            if i == 1 {
                assert!(report.mechanical_energy.numerical_dissipation_j > 1e-5);
            }
            let material = candidate.solid().materials()[0].material();
            close(material.state_coordinate("T").unwrap(), expected_t, 2e-7);
            close(
                material.property(YIELD_STRESS_PROPERTY).unwrap().value_si(),
                yield_stress(expected_t),
                1e-2,
            );
            assert_eq!(
                material
                    .properties()
                    .iter()
                    .map(|p| p.answer().receipt.selected)
                    .collect::<Vec<_>>(),
                original_claims
            );
            let (next_thermal, next_solid, _, _) = candidate.into_parts();
            thermal = next_thermal;
            solid = next_solid;
            temperature = expected_t;
        }
    });
}

#[test]
fn failed_feedback_observation_and_cancellation_do_not_publish_either_memory() {
    let fixture = Fixture::new();
    let thermal = fixture.thermal();
    let solid = with_gate(&CancelGate::new(), |cx| fixture.bind(&thermal, cx));
    let original = solid.solution().clone();
    let original_cells = thermal.cells().clone();
    let observer = Heating {
        watts: RHO * V * CP * 80.0,
        steps: 0,
        fail: false,
    };
    with_gate(&CancelGate::new(), |cx| {
        let advance = |observer: &Heating, config| {
            solid.advance(
                &thermal,
                observer,
                0.5,
                &ROLLERS,
                &[(3, STRESS * V / L)],
                IncrementSettings::default(),
                heat_config(),
                config,
                cx,
            )
        };
        assert!(matches!(
            advance(
                &observer,
                CreepConfig {
                    maximum_iterations: 1,
                    ..config()
                }
            ),
            Err(CreepCouplingError::NotConverged { .. })
        ));
        let failing = Heating {
            fail: true,
            ..observer.clone()
        };
        assert!(matches!(
            advance(&failing, config()),
            Err(CreepCouplingError::Observer(_))
        ));
        assert_eq!(failing.steps, 0);
        let replay = advance(&observer, config()).unwrap();
        assert!(replay.solid().solution().element_states[0].plastic.alpha > 0.0);
        assert!(
            replay.solid().solution().element_states[0]
                .creep
                .equivalent_creep_strain
                > 0.0
        );
    });
    let cancelled = CancelGate::new();
    cancelled.request();
    with_gate(&cancelled, |cx| {
        assert!(matches!(
            solid.advance(
                &thermal,
                &observer,
                0.5,
                &ROLLERS,
                &[(3, STRESS * V / L)],
                IncrementSettings::default(),
                heat_config(),
                config(),
                cx
            ),
            Err(CreepCouplingError::Cancelled)
        ))
    });
    assert_eq!(solid.solution(), &original);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(observer.steps, 0);
}
