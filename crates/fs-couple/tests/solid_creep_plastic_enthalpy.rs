//! Source-driven yielding, permanent creep and actual finite thermal feedback.
//! Coefficients are manufactured implementation references, not measured data.

use fs_couple::solid_enthalpy::creep::{CreepCell, CreepConfig, CreepCouplingError, CreepSolid};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::recovery::CreepRecoveryParameters;
use fs_material::creep::recovery::plastic::source::RecoveringNortonJ2Source;
use fs_material::creep::recovery::source::{
    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
};
use fs_material::creep::source::*;
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::creep::{GAS_CONSTANT_J_MOL_K, StrainHardeningCreep};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
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
        .with_quantity("T", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        )), t)
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
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static,
            )),
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
        Self::with_expansion(None)
    }
    fn with_expansion(alpha: Option<f64>) -> Self {
        Self::with_recovery(alpha, None)
    }
    fn with_recovery(
        alpha: Option<f64>,
        recovery: Option<(StrainHardeningCreep, CreepRecoveryParameters)>,
    ) -> Self {
        Self::with_recovery_softening(alpha, recovery, false)
    }
    fn with_recovery_softening(
        alpha: Option<f64>,
        recovery: Option<(StrainHardeningCreep, CreepRecoveryParameters)>,
        softening: bool,
    ) -> Self {
        let heat_capacity = if softening { 500.0 } else { CP };
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
            ("specific_heat_capacity", CP_DIMS, heat_capacity),
        ] {
            let slope = if softening {
                match name {
                    YOUNG_MODULUS_PROPERTY => -4e6,
                    POISSON_RATIO_PROPERTY => 5e-4,
                    _ => 0.0,
                }
            } else {
                0.0
            };
            claims
                .insert_claim(claim(name, dims, value - 20.0 * slope, value + 80.0 * slope))
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
        if let Some(alpha) = alpha {
            claims
                .insert_claim(claim(
                    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                    Dims([0, 0, 0, -1, 0, 0]),
                    alpha,
                    alpha,
                ))
                .unwrap();
        }
        if let Some((primary, recovery)) = recovery {
            for (name, dims, value) in [
                (
                    CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY,
                    Dims([0; 6]),
                    primary.reference_strain,
                ),
                (
                    CREEP_HARDENING_EXPONENT_PROPERTY,
                    Dims([0; 6]),
                    primary.exponent,
                ),
                (
                    CREEP_RESIDUAL_RATE_FRACTION_PROPERTY,
                    Dims([0; 6]),
                    primary.residual_rate_fraction,
                ),
                (
                    CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
                    CREEP_RATE_DIMS,
                    recovery.reference_rate_s,
                ),
                (
                    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
                    CREEP_ACTIVATION_ENERGY_DIMS,
                    recovery.activation_energy_j_mol,
                ),
                (
                    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
                    Temperature::DIMS,
                    recovery.reference_temperature_k,
                ),
            ] {
                claims.insert_claim(claim(name, dims, value, value)).unwrap();
            }
        }
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
            -20.0 * heat_capacity,
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
    ) -> CreepSolid<'a, NortonJ2Source<'a>> {
        CreepSolid::try_new(
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

    fn expansion_reference(&self) -> ResolvedMaterialStatePoint {
        resolve_material_state_point(
            &self.card,
            &point(300.0),
            &[ScalarPropertyRequirement::try_new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                Dims([0, 0, 0, -1, 0, 0]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }
    fn expanding_source(&self) -> ExpandingCreepSource<NortonJ2Source<'_>> {
        let expansion = self.expansion_reference();
        ExpandingCreepSource::new(&self.card, self.source(), &expansion, || false).unwrap()
    }
    fn bind_expanding<'a>(
        &'a self,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> CreepSolid<'a, ExpandingCreepSource<NortonJ2Source<'a>>> {
        CreepSolid::try_new(
            thermal,
            &NODES,
            &TETS,
            vec![CreepCell {
                cell: CellId(0),
                card: &self.card,
                source: self.expanding_source(),
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

const EXPANSION: f64 = 8e-5;

#[test]
fn free_heated_tetrahedron_follows_the_original_expansion_source_without_stress() {
    with_gate(&CancelGate::new(), |cx| {
        let fixture = Fixture::with_expansion(Some(EXPANSION));
        let thermal = fixture.thermal();
        let solid = fixture.bind_expanding(&thermal, cx);
        let old = solid.solution().clone();
        let selected = solid.materials()[0]
            .thermoelastic()
            .expansion
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        let external_j = RHO * V * CP * 20.0;
        let observer = Heating {
            watts: external_j / 0.5,
            steps: 0,
            fail: false,
        };
        let candidate = solid
            .advance(
                &thermal,
                &observer,
                0.5,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                heat_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(solid.solution(), &old);
        let temperature = candidate.thermal().temperature_k(CellId(0)).unwrap();
        close(temperature, 320.0, 2e-7);
        let solution = candidate.solid().solution();
        let free = EXPANSION * (temperature - 300.0);
        for i in 0..3 {
            close(solution.element_strain[0][i], free, 2e-10);
            close(solution.displacement_m[i + 1][i], L * free, 2e-12);
        }
        assert!(solution.element_stress[0].iter().all(|s| s.abs() < 1e-2));
        let state = &solution.element_states[0];
        close(state.temperature_k(), temperature, 1e-10);
        close(state.mechanical.creep.equivalent_creep_strain, 0.0, 1e-12);
        close(state.mechanical.plastic.alpha, 0.0, 1e-12);
        let material = &candidate.solid().materials()[0].thermoelastic().expansion;
        assert_eq!(material.card_identity(), fixture.card.content_hash());
        assert_eq!(
            material
                .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected,
            selected
        );
        close(material.state_coordinate("T").unwrap(), temperature, 1e-10);
        let report = candidate.report();
        close(report.mechanical_energy.physical_dissipation_j, 0.0, 1e-10);
        close(report.mechanical_energy.reversible_heat_j, 0.0, 1e-9);
        close(
            report.thermal_audit.stored_energy_change_j,
            external_j,
            1e-9,
        );
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);

        // Individually valid sources cannot bypass the combined capacity gate:
        // rho*Cp=20,000, but 9*K*alpha^2*Tmax=22,800 J/(m3 K).
        let excessive = Fixture::with_expansion(Some(1e-4));
        let excessive_thermal = excessive.thermal();
        let old_inventory = excessive_thermal.cells().clone();
        let refused = CreepSolid::try_new(
            &excessive_thermal,
            &NODES,
            &TETS,
            vec![CreepCell {
                cell: CellId(0),
                card: &excessive.card,
                source: excessive.expanding_source(),
                caloric: &excessive.caloric,
            }],
            TetAssemblyBudget::standard(),
            cx,
        );
        assert!(matches!(refused, Err(CreepCouplingError::Invalid(message))
            if message.contains("heat capacity")));
        assert_eq!(excessive_thermal.cells(), &old_inventory);
    });
}

// Axially restrained, transversely free compression: q/E+c+p=alpha(T-300).
// The scalar force balance and finite caloric first law are independent of the
// production tensor return map, finite-element solve and thermoelastic potential.
fn restrained_reference(
    old_t: f64,
    old_q: f64,
    old_c: f64,
    old_p: f64,
    dt: f64,
    external_j: f64,
) -> (f64, f64, f64, f64) {
    let mechanical = |t: f64| {
        let rate = RATE * (ACTIVATION / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / t)).exp();
        let free = EXPANSION * (t - 300.0);
        let increments = |q: f64| {
            (
                dt * rate * (q / STRESS).powi(3),
                ((q - yield_stress(t)) / HARDENING - old_p).max(0.0),
            )
        };
        let (mut lo, mut hi) = (0.0, (E * (free - old_c - old_p)).max(0.0));
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            let (dc, dp) = increments(q);
            if q / E + old_c + dc + old_p + dp > free {
                hi = q;
            } else {
                lo = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let (dc, dp) = increments(q);
        (q, dc, dp)
    };
    let residual = |t: f64| {
        let (q, dc, dp) = mechanical(t);
        RHO * CP * (t - old_t)
            - external_j / V
            - q * dc
            - yield_stress(t) * dp
            - old_t * EXPANSION * (q - old_q)
    };
    let (mut lo, mut hi) = (
        (300.0 + (old_c + old_p) / EXPANSION).max(old_t - 5.0),
        old_t + external_j / (RHO * V * CP) + 5.0,
    );
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
    let (q, dc, dp) = mechanical(t);
    (t, q, dc, dp)
}

#[test]
fn restrained_heating_yields_and_creeps_then_relaxation_returns_heat_atomically() {
    with_gate(&CancelGate::new(), |cx| {
        let fixture = Fixture::with_expansion(Some(EXPANSION));
        let mut thermal = fixture.thermal();
        let mut solid = fixture.bind_expanding(&thermal, cx);
        let mut supports = ROLLERS.to_vec();
        supports.push((3, 0.0));
        supports.sort_unstable_by_key(|(dof, _)| *dof);
        let (mut t0, mut q0, mut creep, mut plastic) = (300.0, 0.0, 0.0, 0.0);
        for (i, supplied_rise) in [40.0, 0.0].into_iter().enumerate() {
            let dt = 0.5;
            let external_j = RHO * V * CP * supplied_rise;
            let target = (i + 1) as f64 * dt;
            let observer = Heating {
                watts: external_j / dt,
                steps: 0,
                fail: false,
            };
            let (t, q, dc, dp) = restrained_reference(t0, q0, creep, plastic, dt, external_j);
            let old_solution = solid.solution().clone();
            let old_cells = thermal.cells().clone();
            if i == 1 {
                // The failure starts from nonzero creep AND plastic memory and
                // the actual accepted temperature, not a virgin fixture.
                let failing = Heating {
                    fail: true,
                    ..observer.clone()
                };
                assert!(matches!(
                    solid.advance(
                        &thermal,
                        &failing,
                        target,
                        &supports,
                        &[],
                        IncrementSettings::default(),
                        heat_config(),
                        config(),
                        cx
                    ),
                    Err(CreepCouplingError::Observer(_))
                ));
                assert_eq!(solid.solution(), &old_solution);
                assert_eq!(thermal.cells(), &old_cells);
                assert_eq!(failing.steps, 0);
            }
            let candidate = solid
                .advance(
                    &thermal,
                    &observer,
                    target,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(solid.solution(), &old_solution);
            creep += dc;
            plastic += dp;
            let solution = candidate.solid().solution();
            let state = &solution.element_states[0];
            close(
                candidate.thermal().temperature_k(CellId(0)).unwrap(),
                t,
                3e-7,
            );
            close(state.temperature_k(), t, 3e-7);
            close(solution.element_stress[0][0], -q, 2e-2);
            assert!(
                solution.element_stress[0][1..]
                    .iter()
                    .all(|s| s.abs() < 2e-2)
            );
            close(solution.displacement_m[1][0], 0.0, 1e-14);
            let lateral = EXPANSION * (t - 300.0) + 0.2 * q / E + 0.5 * (creep + plastic);
            close(solution.displacement_m[2][1], L * lateral, 3e-12);
            close(state.mechanical.creep.creep_strain[0], -creep, 3e-10);
            close(state.mechanical.creep.equivalent_creep_strain, creep, 3e-10);
            close(state.mechanical.plastic.plastic_strain[0], -plastic, 3e-10);
            close(state.mechanical.plastic.alpha, plastic, 3e-10);
            assert!(creep > 0.0 && plastic > 0.0);
            let report = candidate.report();
            let energy = report.element_energy[0];
            let physical = q * dc + yield_stress(t) * dp;
            let reversible = t0 * EXPANSION * (q - q0);
            close(energy.physical_dissipation_j_m3 * V, physical * V, 1e-9);
            close(energy.reversible_heat_j_m3 * V, reversible * V, 1e-9);
            close(
                energy.thermoelastic_storage_change_j_m3 * V,
                EXPANSION * (t0 * q0 - t * q) * V,
                1e-9,
            );
            close(
                report.thermal_audit.stored_energy_change_j,
                external_j + V * (physical + reversible),
                1e-9,
            );
            close(
                report.mechanical_energy.applied_work_j + report.mechanical_energy.support_work_j,
                0.0,
                1e-12,
            );
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            if i == 1 {
                assert!(t < t0 && q < q0);
                assert!(energy.creep_dissipation_j_m3 > 0.0);
                assert!(energy.reversible_heat_j_m3 < 0.0);
                assert!(energy.caloric_heat_j_m3() < 0.0);
                assert!(report.applied_heat_j < 0.0);
                assert!(report.thermal_audit.stored_energy_change_j < 0.0);
            }
            let (next_thermal, next_solid, _, _) = candidate.into_parts();
            thermal = next_thermal;
            solid = next_solid;
            t0 = t;
            q0 = q;
        }
    });
}

fn recovering_joint_reference(softening: bool) {
    let expansion = 1e-6;
    let cp = if softening { 500.0 } else { CP };
    let youngs_slope = if softening { -4e6 } else { 0.0 };
    let poisson_slope = if softening { 5e-4 } else { 0.0 };
    let youngs = |t: f64| E + youngs_slope * (t - 300.0);
    let poisson = |t: f64| 0.2 + poisson_slope * (t - 300.0);
    let bulk = |t: f64| youngs(t) / (3.0 * (1.0 - 2.0 * poisson(t)));
    let shear = |t: f64| youngs(t) / (2.0 * (1.0 + poisson(t)));
    // At a free-sided uniaxial endpoint, nu'(T) cancels from e:C'(T):e.
    // This is entropy at fixed total strain/history, not along stress control.
    let entropy = |t: f64, stress: f64| {
        expansion * stress - 0.5 * youngs_slope * (stress / youngs(t)).powi(2)
    };
    let primary = StrainHardeningCreep {
        reference_strain: 2e-4,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    };
    let recovery = CreepRecoveryParameters {
        reference_rate_s: 2.0,
        activation_energy_j_mol: 15_000.0,
        reference_temperature_k: 300.0,
    };
    let fixture = Fixture::with_recovery_softening(
        Some(expansion), Some((primary, recovery)), softening,
    );
    with_gate(&CancelGate::new(), |cx| {
        let mut thermal = fixture.thermal();
        let (lower, upper) = (point(280.0), point(380.0));
        // CreepSolid admits a softening source only if its conservative
        // mechanical capacity bound, taken over the source span and three
        // strain fields at this limit, leaves rho*cp positive. At 0.05 that
        // bound is about -7.0e5 J/m3/K against 5e5; at 0.02 it is -1.1e5.
        // The trajectory stays inside 0.02, which the law's domain enforces.
        let query = NortonSourceQuery {
            lower: &lower,
            upper: &upper,
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        };
        let source = if softening {
            ExpandingCreepSource::from_recovering_norton_j2(
                &fixture.card, query, HardeningPropertyBasis::PlasticStrainModulus,
                &fixture.expansion_reference(), || false,
            ).unwrap()
        } else {
            let mechanical = RecoveringNortonJ2Source::new(
                &fixture.card, query, HardeningPropertyBasis::PlasticStrainModulus, || false,
            ).unwrap();
            ExpandingCreepSource::new(
                &fixture.card, mechanical, &fixture.expansion_reference(), || false,
            ).unwrap()
        };
        let mut solid = CreepSolid::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![CreepCell {
                cell: CellId(0),
                card: &fixture.card,
                source,
                caloric: &fixture.caloric,
            }],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let claims: Vec<_> = solid.materials()[0]
            .mechanical()
            .material()
            .properties()
            .iter()
            .map(|p| p.answer().receipt.selected)
            .collect();
        assert_eq!(claims.len(), 15);
        let (mut temperature, mut previous_stress, mut time) = (300.0, 0.0, 0.0);
        let (mut creep, mut plastic, mut resistance) = (0.0, 0.0_f64, 0.0);
        for (step, (stress, dt, external_rise)) in [
            (2.3e6_f64, 0.5, 0.0),
            (2.3e6, 0.5, if softening { 15.0 } else { 0.0 }),
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 20.0),
            (2.4e6, 0.5, 0.0),
        ]
        .into_iter()
        .enumerate()
        {
            // Solve the independent constant-stress first law. Primary creep
            // has an exact quadratic; plastic consistency uses its own alpha.
            // The recovered resistance carries neither stored energy nor heat.
            let increments = |t: f64| {
                let k = recovery.reference_rate_s
                    * (recovery.activation_energy_j_mol / GAS_CONSTANT_J_MOL_K
                        * (1.0 / recovery.reference_temperature_k - 1.0 / t))
                        .exp();
                let g = 1.0 + k * dt;
                let rate = RATE
                    * (stress / STRESS).powi(3)
                    * (ACTIVATION / GAS_CONSTANT_J_MOL_K * (1.0 / 300.0 - 1.0 / t)).exp();
                let b = primary.reference_strain * g + resistance;
                let c = dt * rate * primary.reference_strain * g;
                let dc = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
                let dp = ((stress - yield_stress(t)) / HARDENING - plastic).max(0.0);
                (dc, dp, (resistance + dc) / g)
            };
            // Change stiffness at the ORIGINAL elastic tensor, including its
            // original Poisson contraction. No constitutive routine is used.
            let old_axial = previous_stress / youngs(temperature);
            let old_trace = (1.0 - 2.0 * poisson(temperature)) * old_axial;
            let old_dev_norm = 2.0 / 3.0
                * ((1.0 + poisson(temperature)) * old_axial).powi(2);
            let coefficient_storage = |t: f64| {
                0.5 * (bulk(t) - bulk(temperature)) * old_trace.powi(2)
                    + (shear(t) - shear(temperature)) * old_dev_norm
            };
            let entropy_storage = |t: f64| {
                t * entropy(t, stress) - temperature * entropy(temperature, previous_stress)
            };
            let reversible = |t: f64| {
                stress * expansion * (t - temperature)
                    - coefficient_storage(t) - entropy_storage(t)
            };
            let residual = |t: f64| {
                let (dc, dp, _) = increments(t);
                RHO * cp * (t - temperature - external_rise)
                    - stress * dc
                    - yield_stress(t) * dp
                    - reversible(t)
            };
            let (mut lo, mut hi) = (
                temperature + external_rise - 2.0,
                temperature + external_rise + 2.0,
            );
            assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let next_temperature = 0.5 * (lo + hi);
            let (dc, dp, next_resistance) = increments(next_temperature);
            let reciprocal_heat = reversible(next_temperature);
            let external_j = RHO * V * cp * external_rise;
            let observer = Heating {
                watts: external_j / dt,
                steps: 0,
                fail: false,
            };
            let old_solution = solid.solution().clone();
            let old_cells = thermal.cells().clone();
            let candidate = solid
                .advance(
                    &thermal,
                    &observer,
                    time + dt,
                    &ROLLERS,
                    &[(3, stress * V / L)],
                    IncrementSettings::default(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(solid.solution(), &old_solution);
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(observer.steps, 0);
            let solution = candidate.solid().solution();
            let state = &solution.element_states[0];
            let history = &state.mechanical;
            close(
                history.native.creep.equivalent_creep_strain,
                creep + dc,
                3e-10,
            );
            close(history.native.creep.creep_strain[0], creep + dc, 3e-10);
            close(history.native.plastic.alpha, plastic + dp, 3e-10);
            close(history.native.plastic.plastic_strain[0], plastic + dp, 3e-10);
            close(history.resistance_strain, next_resistance, 3e-10);
            close(state.temperature_k(), next_temperature, 3e-7);
            close(
                candidate.thermal().temperature_k(CellId(0)).unwrap(),
                next_temperature,
                3e-7,
            );
            close(solution.element_stress[0][0], stress, 2e-2);
            close(
                solution.displacement_m[1][0],
                L * (stress / youngs(next_temperature)
                    + creep
                    + dc
                    + plastic
                    + dp
                    + expansion * (next_temperature - 300.0)),
                4e-12,
            );
            let transverse = -poisson(next_temperature) * stress / youngs(next_temperature)
                - 0.5 * (creep + dc + plastic + dp)
                + expansion * (next_temperature - 300.0);
            close(solution.displacement_m[2][1], L * transverse, 4e-12);
            close(solution.displacement_m[3][2], L * transverse, 4e-12);
            let report = candidate.report();
            let energy = report.element_energy[0];
            let physical = stress * dc + yield_stress(next_temperature) * dp;
            close(energy.creep_dissipation_j_m3 * V, stress * dc * V, 1e-9);
            close(
                energy.plastic_dissipation_j_m3 * V,
                yield_stress(next_temperature) * dp * V,
                1e-9,
            );
            close(energy.reversible_heat_j_m3 * V, reciprocal_heat * V, 1e-9);
            close(
                energy.thermoelastic_storage_change_j_m3 * V,
                entropy_storage(next_temperature) * V,
                1e-9,
            );
            close(
                energy.final_stored_j_m3 * V,
                (0.5 * stress * stress / youngs(next_temperature)
                    + 0.5 * HARDENING * (plastic + dp).powi(2)) * V,
                1e-9,
            );
            close(
                energy.initial_stored_j_m3 * V,
                (0.5 * previous_stress * previous_stress / youngs(temperature)
                    + 0.5 * HARDENING * plastic.powi(2)) * V,
                1e-9,
            );
            close(report.applied_heat_j, V * (physical + reciprocal_heat), 1e-9);
            close(
                report.thermal_audit.stored_energy_change_j,
                external_j + V * (physical + reciprocal_heat),
                1e-9,
            );
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            if stress == 0.0 {
                assert!(next_resistance < resistance);
                close(energy.physical_dissipation_j_m3 * V, 0.0, 1e-12);
                for i in 0..6 {
                    let before = &old_solution.element_states[0].mechanical.native;
                    close(
                        history.native.creep.creep_strain[i],
                        before.creep.creep_strain[i],
                        3e-10,
                    );
                    close(
                        history.native.plastic.plastic_strain[i],
                        before.plastic.plastic_strain[i],
                        3e-10,
                    );
                }
                if step == 2 {
                    assert!(energy.numerical_dissipation_j_m3 > 0.0);
                    assert!(energy.reversible_heat_j_m3 > 0.0);
                } else {
                    close(report.applied_heat_j, 0.0, 1e-12);
                }
            } else if step == 0 || step == 4 {
                assert!(dc > 0.0 && dp > 0.0);
            }
            assert_eq!(
                candidate.solid().materials()[0]
                    .mechanical()
                    .material()
                    .properties()
                    .iter()
                    .map(|p| p.answer().receipt.selected)
                    .collect::<Vec<_>>(),
                claims
            );
            (thermal, solid, _, _) = candidate.into_parts();
            creep += dc;
            plastic += dp;
            resistance = next_resistance;
            temperature = next_temperature;
            previous_stress = stress;
            time += dt;
        }
    });
}

#[test]
fn recovering_joint_solid_keeps_yield_history_through_heated_dwell_and_reload() {
    recovering_joint_reference(false);
}

#[test]
fn recovering_joint_solid_softens_under_load_and_closes_sourced_thermal_feedback() {
    recovering_joint_reference(true);
}
