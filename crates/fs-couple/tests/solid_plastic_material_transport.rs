//! Native global equilibrium, original caloric/source transport and energy feedback.
#[path = "solid_plastic_material_transport/heat.rs"]
mod heat;

use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialNetwork};
use fs_couple::solid_enthalpy::plastic::{
    ThermoplasticCell, ThermoplasticConfig, ThermoplasticError, ThermoplasticSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::thermal::SourceQuery;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, THERMAL_CONDUCTIVITY_DIMS,
    YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent,
};
use std::collections::BTreeMap;

const L: f64 = 0.01;
const VOLUME: f64 = L * L * L / 6.0;
const MASS: f64 = 1000.0 * VOLUME;
const LINK_AREA: f64 = L * L / 2.0;
const SHEAR: f64 = 0.01;
const SHEAR_STRESS_ZX: f64 = 0.5e6;
const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [L, 0.0, 0.0],
    [0.0, L, 0.0],
    [0.0, 0.0, L],
    [0.0, 0.0, -L],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 4]];
const CELLS: [CellId; 2] = [CellId(17), CellId(3)];
const MATS: [MaterialId; 2] = [MaterialId(11), MaterialId(29)];
const CP: &str = "specific_heat_capacity";
const K: &str = "thermal_conductivity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const PRESSURE_PA: f64 = 200000.0;
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic spatial thermoplastic regression".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn pressure_surface(property: &mut PropertyClaim, scales: [f64; 2]) {
    let knots = match &property.value {
        PropertyValue::Scalar { value, .. } => vec![(280.0, *value), (320.0, *value)],
        PropertyValue::Curve { knots, .. } => knots.clone(),
        PropertyValue::Surface2D { .. } => unreachable!("fixture starts from scalar/T data"),
    };
    property.value = PropertyValue::Surface2D {
        axes: ["T".into(), "p".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: knots
            .iter()
            .flat_map(|&(t, value)| [([t, 1e5], scales[0] * value), ([t, 3e5], scales[1] * value)])
            .collect(),
        triangles: (0..knots.len() as u32 - 1)
            .flat_map(|index| {
                let left = 2 * index;
                [[left, left + 2, left + 1], [left + 2, left + 3, left + 1]]
            })
            .collect(),
        dims: property.key.quantity().dims(),
    };
    property.validity = ValidityDomain::unconstrained()
        .with("T", 280.0, 320.0)
        .with_quantity("p", pressure_quantity(), 1e5, 3e5);
    property.interpolation = InterpolationPolicy::LinearInside;
}

struct Fixture {
    cards: [MaterialCard; 2],
    caloric: [ResolvedHeatCapacityCurve; 2],
    point: QueryPoint,
    conductivity: PropertyKey,
}
impl Fixture {
    fn new(cp: [f64; 2], k: f64, conductivity_gap: bool) -> Self {
        Self::with_transport_axis(cp, k, conductivity_gap, None)
    }
    fn with_transport_axis(
        cp: [f64; 2],
        k: f64,
        conductivity_gap: bool,
        mechanical_axis: Option<&str>,
    ) -> Self {
        Self::with_conditions(cp, k, conductivity_gap, mechanical_axis, None)
    }
    fn with_pressure(cp: [f64; 2], k: f64, pressure_pa: f64) -> Self {
        Self::with_conditions(cp, k, false, None, Some(pressure_pa))
    }
    fn with_conditions(
        cp: [f64; 2],
        k: f64,
        conductivity_gap: bool,
        mechanical_axis: Option<&str>,
        pressure_pa: Option<f64>,
    ) -> Self {
        let cards = [0, 1].map(|side| {
            let mut claims = ClaimSet::new();
            for mut property in [
                scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
                scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
                scalar(BASIS.property_name(), Pressure::DIMS, 1e8),
                scalar(DENSITY_PROPERTY, Density::DIMS, 1000.0),
                scalar(CP, CP_DIMS, cp[side]),
            ] {
                if pressure_pa.is_some() {
                    let scales = if property.key.name() == CP {
                        [1.0, 2.0]
                    } else {
                        [1.0, 1.0]
                    };
                    pressure_surface(&mut property, scales);
                }
                claims.insert_claim(property).unwrap();
            }
            let mut strength = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 40e6);
            if side == 0 {
                strength.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 5e6), (300.0, 3e6), (320.0, 1e6)],
                    dims: Pressure::DIMS,
                };
                strength.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
                strength.interpolation = InterpolationPolicy::LinearInside;
            }
            if pressure_pa.is_some() {
                pressure_surface(&mut strength, [1.0, 1.0]);
            }
            claims.insert_claim(strength).unwrap();
            let mut transport = scalar(K, THERMAL_CONDUCTIVITY_DIMS, k);
            transport.validity = ValidityDomain::unconstrained().with("p", 199999.0, 200001.0);
            if let Some(axis) = mechanical_axis {
                transport.validity = transport.validity.with(axis, 0.0, 0.001);
            }
            if conductivity_gap && side == 0 {
                transport.validity = transport.validity.with("T", 280.0, 300.1);
            }
            if pressure_pa.is_some() {
                pressure_surface(&mut transport, [1.0, 10.0]);
            }
            claims.insert_claim(transport).unwrap();
            if conductivity_gap && side == 0 {
                let mut alternate = scalar(K, THERMAL_CONDUCTIVITY_DIMS, 10.0 * k);
                alternate.validity = ValidityDomain::unconstrained()
                    .with("p", 199999.0, 200001.0)
                    .with("T", 300.11, 320.0);
                claims.insert_claim(alternate).unwrap();
            }
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("synthetic tet {side}"),
                    phase: "solid".into(),
                    process: "numerical fixture".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap()
        });
        let mut point = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with_quantity("p", pressure_quantity(), pressure_pa.unwrap_or(PRESSURE_PA))
            .unwrap();
        if let Some(axis) = mechanical_axis {
            point = point.with(axis, 0.0).unwrap();
        }
        let lower = point.clone().with("T", 280.0).unwrap();
        let upper = point.clone().with("T", 320.0).unwrap();
        let cp_key = PropertyKey::new(CP, CP_DIMS);
        let density_key = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
        let caloric = [0, 1].map(|side| {
            let reference_cp = cards[side]
                .claims()
                .query_typed(&cp_key, &point, SelectionPolicy::SingleClaimOnly)
                .unwrap()
                .evidence
                .value
                .value;
            ResolvedHeatCapacityCurve::try_new(
                &cards[side],
                HeatCapacitySourceQuery {
                    heat_capacity: &cp_key,
                    density: &density_key,
                    temperature_axis: "T",
                    lower: &lower,
                    upper: &upper,
                    selection: SelectionPolicy::SingleClaimOnly,
                },
                SolidLiquidPhase::Solid,
                -20.0 * reference_cp,
                1e-4,
            )
            .unwrap()
        });
        Self {
            cards,
            caloric,
            point,
            conductivity: PropertyKey::new(K, THERMAL_CONDUCTIVITY_DIMS),
        }
    }
    fn sources(&self) -> Vec<ThermoplasticCell<'_>> {
        [0, 1]
            .map(|side| ThermoplasticCell {
                cell: CELLS[side],
                query: SourceQuery {
                    card: &self.cards[side],
                    caloric: &self.caloric[side],
                    context: self.point.clone(),
                    temperature_axis: "T".into(),
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    hardening_basis: BASIS,
                    thermal_expansion: None,
                    strain_limit: 0.04,
                },
            })
            .into()
    }
    fn thermal(&self) -> DepositionSimulation<'_> {
        thermal_with([self.caloric[0].curve(), self.caloric[1].curve()])
    }
    fn solid<'a>(
        &'a self,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> ThermoplasticSolid<'a> {
        ThermoplasticSolid::try_new(
            thermal,
            &NODES,
            &TETS,
            self.sources(),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap()
    }
}

fn surface() -> f64 {
    (3.0 + 3.0_f64.sqrt()) * L * L / 2.0
}
fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn thermal_with(curves: [&dyn ThermalLaw; 2]) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = [0, 1]
        .map(|side| {
            let center = [L / 4.0, L / 4.0, if side == 0 { L / 4.0 } else { -L / 4.0 }];
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CELLS[side],
                    segment_id: side as u64,
                    body_id: 1,
                    material_id: MATS[side],
                    mass_kg: MASS,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: center,
                        end_m: center,
                        volume_m3: VOLUME,
                        surface_area_m2: surface(),
                    },
                    boundary: boundary(),
                }),
            }
        })
        .into();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 91,
            a: CELLS[0],
            b: CELLS[1],
            area_m2: LINK_AREA,
            // Native material observer replaces this on every accepted substep.
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MATS[0], curves[0]), (MATS[1], curves[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-10,
    }
}
fn control() -> ThermoplasticConfig {
    ThermoplasticConfig {
        maximum_iterations: 80,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        enthalpy_relaxation: 1.0,
    }
}
fn supports(shear: f64) -> Vec<(usize, f64)> {
    // Only the upper apex's x displacement is free. Its actual equilibrium
    // changes with local temperature/yield data; it is not a prescribed strain.
    (0..15)
        .filter(|&d| d != 9)
        .map(|d| {
            (
                d,
                if d % 3 == 0 {
                    2.0 * shear * NODES[d / 3][1]
                } else {
                    0.0
                },
            )
        })
        .collect()
}
fn applied() -> [(usize, f64); 1] {
    [(9, SHEAR_STRESS_ZX * VOLUME / L)]
}
fn links() -> [ConductionLink; 1] {
    [ConductionLink {
        link: 91,
        distance_a_m: L / 4.0,
        distance_b_m: L / 4.0,
    }]
}
fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 0x53504a32,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}
#[derive(Clone, Default, Debug, PartialEq)]
struct Recorder {
    durations: Vec<f64>,
}
impl DepositionObserver for Recorder {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.durations.push(dt);
        Ok(())
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

// Independent reduction for this fixture: two nonzero tensor shear components,
// one enforced shear traction at the free apex, and scalar energetic hardening.
// This solves the closed-form stress relation rather than calling native J2.
fn apex_oracle(temperature_k: f64) -> (f64, f64, f64) {
    let fy = 3e6 - 1e5 * (temperature_k - 300.0);
    let g = 5e8;
    let h = 1e8;
    let sigma_zx = |gamma: f64| {
        let q_trial = 2.0 * g * 3.0_f64.sqrt() * (SHEAR * SHEAR + gamma * gamma).sqrt();
        2.0 * g * gamma * (h + 3.0 * g * fy / q_trial) / (3.0 * g + h)
    };
    let (mut low, mut high) = (0.0, 0.03);
    assert!(sigma_zx(high) > SHEAR_STRESS_ZX);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if sigma_zx(mid) < SHEAR_STRESS_ZX {
            low = mid;
        } else {
            high = mid;
        }
    }
    let gamma = 0.5 * (low + high);
    let q_trial = 2.0 * g * 3.0_f64.sqrt() * (SHEAR * SHEAR + gamma * gamma).sqrt();
    (gamma, (q_trial - fy) / (3.0 * g + h), fy)
}

// The two-cell conservative heat matrix has one zero eigenvalue and one
// decay rate G*(1/Ca+1/Cb). Its exact 20-step FE response gives K per input J.
fn heat_response(k: f64, cp: [f64; 2]) -> [f64; 2] {
    let g = LINK_AREA / (L / (4.0 * k) + L / (4.0 * k));
    let [ca, cb] = cp.map(|c| MASS * c);
    let sum = ca + cb;
    let retained = (1.0 - 0.05 * g * (1.0 / ca + 1.0 / cb)).powi(20);
    [
        1.0 / sum + cb * cb * (1.0 - retained) / (g * sum * sum),
        1.0 / sum - ca * cb * (1.0 - retained) / (g * sum * sum),
    ]
}
fn coupled_oracle(k: f64, cp: [f64; 2]) -> (f64, f64, f64, f64) {
    let response = heat_response(k, cp);
    let (mut low, mut high) = (300.0, 320.0);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        let (_, alpha, fy) = apex_oracle(mid);
        if mid - 300.0 < response[0] * VOLUME * fy * alpha {
            low = mid;
        } else {
            high = mid;
        }
    }
    let ta = 0.5 * (low + high);
    let (gamma, alpha, fy) = apex_oracle(ta);
    (
        ta,
        300.0 + response[1] * VOLUME * fy * alpha,
        2.0 * L * gamma,
        alpha,
    )
}

#[test]
fn sourced_conductivity_swap_changes_spatial_heat_exchange_and_native_free_displacement() {
    with_cx(|_, cx| {
        let mut observations = Vec::new();
        for k in [0.01, 0.1] {
            let fixture = Fixture::new([5.0, 5.0], k, false);
            let thermal = fixture.thermal();
            let solid = fixture.solid(&thermal, cx);
            let conductivity = vec![
                Some((
                    fixture.conductivity.clone(),
                    MaterialPropertySelection::SingleClaimOnly
                ));
                2
            ];
            let materials = solid.thermal_material_cells(&conductivity).unwrap();
            let links = links();
            let network = MaterialNetwork {
                cells: &materials,
                resistors: &[],
                links: &links,
            };
            let prepared = network.prepare(&thermal, || false).unwrap();
            let initial_sample = prepared.resolve(&thermal, || false).unwrap().source;
            let observer = prepared.observe(Recorder::default()).unwrap();
            let supports = supports(SHEAR);
            let force = applied();
            let step = solid
                .advance(
                    &thermal,
                    &observer,
                    1.0,
                    &supports,
                    &force,
                    IncrementSettings::default(),
                    thermal_config(),
                    control(),
                    cx,
                )
                .unwrap();
            let expected = coupled_oracle(k, [5.0, 5.0]);
            let ta = step.thermal().temperature_k(CELLS[0]).unwrap();
            let tb = step.thermal().temperature_k(CELLS[1]).unwrap();
            let solution = step.solid().solution();
            let report = step.report();
            close(ta, expected.0, 3e-7);
            close(tb, expected.1, 3e-7);
            close(solution.displacement_m[3][0], expected.2, 2e-11);
            close(solution.element_states[0].alpha, expected.3, 2e-10);
            close(solution.element_stress[0][5], SHEAR_STRESS_ZX, 1e-2);
            assert!(ta > tb && tb > 300.0);
            assert!(solution.element_states[0].alpha > 0.0);
            assert_eq!(solution.element_states[1].alpha, 0.0);
            assert_eq!(
                report.element_energy[1].energy.physical_dissipation_j_m3,
                0.0
            );
            assert!(report.numerical_loss_j > report.physical_heat_j);
            close(report.external_heat_j, 0.0, 1e-10);
            close(
                report.thermal_audit.stored_energy_change_j,
                report.physical_heat_j,
                control().energy_tolerance_j,
            );
            let local_work: f64 = report
                .element_energy
                .iter()
                .map(|e| VOLUME * e.energy.mechanical_work_j_m3)
                .sum();
            close(
                report.mechanical_work_j,
                local_work,
                control().energy_tolerance_j,
            );
            assert!(report.heat_defect_j <= control().energy_tolerance_j);
            assert!(report.combined_residual_j.abs() <= control().energy_tolerance_j);
            assert!(report.iterations > 1);
            close(step.observer().second.durations.iter().sum(), 1.0, 1e-12);
            assert_eq!(
                step.observer().second.durations.len(),
                report.thermal_audit.substeps
            );
            assert!(observer.second.durations.is_empty());
            let sample = step.observer().first.sample().unwrap();
            assert_eq!(sample.links[0].cells, CELLS);
            close(
                sample.links[0].conductance_w_k,
                LINK_AREA / (L / (2.0 * k)),
                1e-14,
            );
            for (element, material) in step.solid().materials().iter().enumerate() {
                let id = CELLS[element];
                let temperature = step.thermal().temperature_k(id).unwrap();
                let cell = sample.cells.iter().find(|c| c.cell == id).unwrap();
                let initial = initial_sample.cells.iter().find(|c| c.cell == id).unwrap();
                let transport = cell.conductivity.as_ref().unwrap();
                assert_eq!(
                    material.material().card_identity(),
                    fixture.cards[element].content_hash()
                );
                assert_eq!(
                    transport.card_identity(),
                    material.material().card_identity()
                );
                assert_eq!(
                    cell.phase.specific_enthalpy_j_kg(),
                    step.thermal().cells()[&id].specific_h_j_kg
                );
                assert_eq!(cell.phase.temperature_k(), temperature);
                assert_eq!(
                    transport.property(K).unwrap().answer().receipt.selected,
                    initial
                        .conductivity
                        .as_ref()
                        .unwrap()
                        .property(K)
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
                for bundle in [material.material(), transport] {
                    for property in bundle.properties() {
                        let point = &property.answer().receipt.query_point;
                        assert_eq!(
                            point.iter().find(|(axis, _)| axis == "T").unwrap().1,
                            temperature
                        );
                        assert_eq!(
                            point.iter().find(|(axis, _)| axis == "p").unwrap().1,
                            PRESSURE_PA
                        );
                    }
                }
            }
            assert_eq!(step.solid().time_s(), step.thermal().time_s());
            assert_eq!(solid.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            observations.push((ta, tb, solution.displacement_m[3][0]));
        }
        assert!(observations[1].0 < observations[0].0 - 0.1);
        assert!(observations[1].1 > observations[0].1 + 0.1);
        assert!(observations[1].2 < observations[0].2);
    });
}

#[test]
fn caloric_data_swap_changes_equilibrium_and_unloading_preserves_local_plastic_history() {
    with_cx(|_, cx| {
        let mut temperatures = Vec::new();
        for cp in [[5.0, 5.0], [10.0, 5.0]] {
            let fixture = Fixture::new(cp, 0.01, false);
            let thermal = fixture.thermal();
            // Keep this immutable source/context anchor alive across candidates.
            let solid = fixture.solid(&thermal, cx);
            let conductivity = vec![
                Some((
                    fixture.conductivity.clone(),
                    MaterialPropertySelection::SingleClaimOnly
                ));
                2
            ];
            let materials = solid.thermal_material_cells(&conductivity).unwrap();
            let links = links();
            let network = MaterialNetwork {
                cells: &materials,
                resistors: &[],
                links: &links,
            };
            let prepared = network.prepare(&thermal, || false).unwrap();
            let observer = prepared.observe(Recorder::default()).unwrap();
            let prescribed = supports(SHEAR);
            let forces = applied();
            let loaded = solid
                .advance(
                    &thermal,
                    &observer,
                    1.0,
                    &prescribed,
                    &forces,
                    IncrementSettings::default(),
                    thermal_config(),
                    control(),
                    cx,
                )
                .unwrap();
            let oracle = coupled_oracle(0.01, cp);
            close(
                loaded.thermal().temperature_k(CELLS[0]).unwrap(),
                oracle.0,
                3e-7,
            );
            close(
                loaded.solid().solution().displacement_m[3][0],
                oracle.2,
                2e-11,
            );
            temperatures.push(oracle.0);
            let histories = loaded.solid().solution().element_states.clone();
            let prescribed = supports(SHEAR - 0.0005);
            let unloaded = loaded
                .solid()
                .advance(
                    loaded.thermal(),
                    loaded.observer(),
                    2.0,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    control(),
                    cx,
                )
                .unwrap();
            let solution = unloaded.solid().solution();
            let report = unloaded.report();
            assert_eq!(solution.element_states, histories);
            assert!(
                report
                    .element_energy
                    .iter()
                    .all(|e| e.energy.physical_dissipation_j_m3 == 0.0)
            );
            assert!(report.numerical_loss_j > 0.0);
            assert!(report.mechanical_work_j < 0.0);
            close(
                solution.displacement_m[3][0],
                2.0 * L * histories[0].plastic_strain[5],
                2e-11,
            );
            assert!(solution.displacement_m[3][0] > 0.0);
            assert!(
                unloaded.thermal().temperature_k(CELLS[0]).unwrap()
                    < loaded.thermal().temperature_k(CELLS[0]).unwrap()
            );
            assert!(
                unloaded.thermal().temperature_k(CELLS[1]).unwrap()
                    > loaded.thermal().temperature_k(CELLS[1]).unwrap()
            );
            close(report.thermal_audit.stored_energy_change_j, 0.0, 1e-10);
            close(
                unloaded.observer().second.durations.iter().sum(),
                2.0,
                1e-12,
            );
            close(loaded.observer().second.durations.iter().sum(), 1.0, 1e-12);
            assert_eq!(loaded.solid().time_s(), 1.0);
            assert_eq!(unloaded.solid().time_s(), 2.0);
            assert!(report.combined_residual_j.abs() <= control().energy_tolerance_j);
        }
        assert!(temperatures[1] < temperatures[0] - 0.5);
    });
}

#[test]
fn heating_past_the_original_conductivity_claim_refuses_the_whole_coupled_interval() {
    with_cx(|_, cx| {
        let fixture = Fixture::new([5.0, 5.0], 0.01, true);
        let thermal = fixture.thermal();
        let solid = fixture.solid(&thermal, cx);
        let conductivity = vec![
            Some((
                fixture.conductivity.clone(),
                MaterialPropertySelection::SingleClaimOnly
            ));
            2
        ];
        let materials = solid.thermal_material_cells(&conductivity).unwrap();
        let links = links();
        let network = MaterialNetwork {
            cells: &materials,
            resistors: &[],
            links: &links,
        };
        let prepared = network.prepare(&thermal, || false).unwrap();
        let observer = prepared.observe(Recorder::default()).unwrap();
        let prescribed = supports(SHEAR);
        let forces = applied();
        let old_cells = thermal.cells().clone();
        let old_solution = solid.solution().clone();
        let old_materials: Vec<_> = solid
            .materials()
            .iter()
            .map(|m| m.material().identity())
            .collect();
        // Cp and yield cover the whole chart; only the original selected k
        // leaves support. A hot k claim exists but may not replace the pin.
        let result = solid.advance(
            &thermal,
            &observer,
            1.0,
            &prescribed,
            &forces,
            IncrementSettings::default(),
            thermal_config(),
            control(),
            cx,
        );
        assert!(matches!(
            result,
            Err(ThermoplasticError::Observer(
                fs_thermal::observers::ObserverCouplingError::First(_)
            ))
        ));
        assert!(observer.second.durations.is_empty());
        assert!(observer.first.sample().is_none());
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.solution(), &old_solution);
        assert_eq!(
            solid
                .materials()
                .iter()
                .map(|m| m.material().identity())
                .collect::<Vec<_>>(),
            old_materials
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
        for cell in CELLS {
            assert_eq!(thermal.temperature_k(cell).unwrap(), 300.0);
        }
    });
}

#[test]
fn pressure_surface_data_changes_actual_plastic_equilibrium_and_material_conduction() {
    with_cx(|_, cx| {
        let mut observations = Vec::new();
        let mut card_ids = Vec::new();
        for pressure in [1.5e5, 2.5e5] {
            let fixture = Fixture::with_pressure([5.0, 5.0], 0.01, pressure);
            let thermal = fixture.thermal();
            let solid = fixture.solid(&thermal, cx);
            card_ids.push(fixture.cards.each_ref().map(MaterialCard::content_hash));
            let conductivity = vec![
                Some((
                    fixture.conductivity.clone(),
                    MaterialPropertySelection::SingleClaimOnly
                ));
                2
            ];
            let materials = solid.thermal_material_cells(&conductivity).unwrap();
            let links = links();
            let network = MaterialNetwork {
                cells: &materials,
                resistors: &[],
                links: &links,
            };
            let prepared = network.prepare(&thermal, || false).unwrap();
            let observer = prepared.observe(()).unwrap();
            let step = solid
                .advance(
                    &thermal,
                    &observer,
                    1.0,
                    &supports(SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    thermal_config(),
                    control(),
                    cx,
                )
                .unwrap();
            let fraction = (pressure - 1e5) / 2e5;
            let cp = 5.0 * (1.0 + fraction);
            let k = 0.01 * (1.0 + 9.0 * fraction);
            let expected = coupled_oracle(k, [cp, cp]);
            let actual = (
                step.thermal().temperature_k(CELLS[0]).unwrap(),
                step.thermal().temperature_k(CELLS[1]).unwrap(),
                step.solid().solution().displacement_m[3][0],
                step.solid().solution().element_states[0].alpha,
            );
            close(actual.0, expected.0, 3e-7);
            close(actual.1, expected.1, 3e-7);
            close(actual.2, expected.2, 2e-11);
            close(actual.3, expected.3, 2e-10);
            assert!(actual.0 > actual.1 && actual.1 > 300.0);
            assert!(actual.3 > 0.0);
            assert_eq!(step.solid().solution().element_states[1].alpha, 0.0);
            let report = step.report();
            assert!(report.physical_heat_j > 0.0);
            assert!(report.numerical_loss_j > 0.0);
            close(report.external_heat_j, 0.0, 1e-10);
            close(
                report.thermal_audit.stored_energy_change_j,
                report.physical_heat_j,
                control().energy_tolerance_j,
            );
            close(
                report.mechanical_work_j,
                report.stored_mechanical_change_j
                    + report.physical_heat_j
                    + report.numerical_loss_j,
                control().energy_tolerance_j,
            );
            assert!(report.combined_residual_j.abs() <= control().energy_tolerance_j);
            let transport_sample = step.observer().first.sample().unwrap();
            for (index, material) in step.solid().materials().iter().enumerate() {
                let cell = transport_sample
                    .cells
                    .iter()
                    .find(|cell| cell.cell == CELLS[index])
                    .unwrap();
                let transport = cell.conductivity.as_ref().unwrap();
                close(transport.property(K).unwrap().value_si(), k, 1e-14);
                assert_eq!(
                    transport.card_identity(),
                    material.material().card_identity()
                );
                for bundle in [material.material(), transport] {
                    for property in bundle.properties() {
                        let receipt = &property.answer().receipt;
                        fixture.cards[index]
                            .claims()
                            .verify_receipt(receipt)
                            .unwrap();
                        assert!(receipt.query_point.contains(&("p".into(), pressure)));
                        assert!(
                            receipt
                                .query_point
                                .contains(&("T".into(), cell.phase.temperature_k()))
                        );
                        assert_eq!(receipt.axis_quantities.get("p"), Some(&pressure_quantity()));
                    }
                }
                for property in material.material().properties() {
                    assert_eq!(
                        property.answer().receipt.selected,
                        solid.materials()[index]
                            .material()
                            .property(property.requirement().name())
                            .unwrap()
                            .answer()
                            .receipt
                            .selected
                    );
                }
            }
            assert_eq!(solid.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            observations.push(actual);
        }
        assert_eq!(
            card_ids[0], card_ids[1],
            "the same material cards serve both pressure contexts"
        );
        assert!((observations[0].0 - observations[1].0).abs() > 0.01);
        assert!((observations[0].2 - observations[1].2).abs() > 1e-8);
    });
}

#[test]
fn transport_bounded_by_frozen_mechanical_coordinates_requires_an_evolving_response() {
    with_cx(|_, cx| {
        for axis in ["strain", "stress", "strain_rate", "omega"] {
            let fixture = Fixture::with_transport_axis([5.0, 5.0], 0.01, false, Some(axis));
            let thermal = fixture.thermal();
            // The six native coefficients do not depend on this coordinate.
            // Only k has a mechanical validity bound; freezing it is not valid
            // while this owner's actual strain, stress and history evolve.
            let solid = fixture.solid(&thermal, cx);
            let conductivity = vec![
                Some((
                    fixture.conductivity.clone(),
                    MaterialPropertySelection::SingleClaimOnly
                ));
                2
            ];
            assert!(matches!(
                solid.thermal_material_cells(&conductivity),
                Err(ThermoplasticError::Material(_))
            ));
            assert_eq!(solid.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert!(
                solid
                    .solution()
                    .element_states
                    .iter()
                    .all(|state| state.alpha == 0.0)
            );
        }
    });
}
