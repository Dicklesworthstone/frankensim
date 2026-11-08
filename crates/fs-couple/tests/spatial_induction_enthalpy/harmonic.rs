//! Independent AC diffusion and two-body caloric balance, with source feedback.
#[path = "harmonic/circuit.rs"]
mod circuit;
#[path = "harmonic/terminals.rs"]
mod terminals;

#[cfg(feature = "harmonic-winding-heating")]
#[path = "harmonic/coherent.rs"]
mod coherent;

use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    contacts::{ContactThermalCoefficients, heat::*, material::ElectricalContactNetwork},
    magnetic::MagneticThermalError,
    material::{ConductionLink, MaterialNetwork},
    spatial_induction::{
        SpatialInductionHeatError,
        harmonic::{SpatialHarmonicThermalCandidate, heat_assembly_candidate},
    },
};
use fs_em::spatial::vector_potential::induction::harmonic::SpatialHarmonicInductionControls;
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::*;
use fs_material::{
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::interface::ThermalInterfacePoint,
};
use fs_math::c64::C64;
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::cell::Cell;

const BODY: CellId = CellId(7);
const SINK: CellId = CellId(8);
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;
const RHO: f64 = 0.2;
const ALPHA: f64 = 0.009;
const DT: f64 = 2.0;
const OMEGA: f64 = 80.0;
const PEAK: f64 = 8.0;
const EMISSIVITY: [f64; 2] = [0.2, 0.4];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn card(variable_mu: bool, rho: f64, emissivity: f64) -> MaterialCard {
    card_with_alpha(variable_mu, rho, emissivity, ALPHA)
}

fn card_with_alpha(variable_mu: bool, rho: f64, emissivity: f64, alpha: f64) -> MaterialCard {
    let original = if variable_mu {
        source::magnetic_curve_card(rho, 1.0, BASIS)
    } else {
        source::card(rho, alpha, 1.0, BASIS)
    };
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        claim.validity = claim
            .validity
            .with_quantity("T", absolute_temperature(), 200.0, 2000.0);
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, value) in [
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature(),
                    200.0,
                    2000.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic settled AC thermal transport".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn with_case<R>(
    variable_mu: bool,
    resistance: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        &MaterialCard,
    ) -> R,
) -> R {
    with_case_with_alpha(variable_mu, resistance, ALPHA, run)
}

fn with_case_with_alpha<R>(
    variable_mu: bool,
    resistance: f64,
    alpha: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        &MaterialCard,
    ) -> R,
) -> R {
    let cards = [
        card_with_alpha(variable_mu, RHO, EMISSIVITY[0], alpha),
        card_with_alpha(false, RHO, EMISSIVITY[1], alpha),
    ];
    let curves = cards.iter().map(curve).collect::<Vec<_>>();
    let context = source::absolute_point(300.0);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: [BODY, SINK][i],
        material_id: MaterialId(9 + i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
    });
    let links = [ConductionLink {
        link: 17,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &links,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Curve {
                abscissa: "aT".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 0.9 * resistance), (2000.0, 2.7 * resistance)],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("aT", temperature, 200.0, 2000.0)
                .with_quantity("bT", temperature, 200.0, 2000.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic AC heater thermal contact".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "body".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "sink".into(),
        },
        SystemContext {
            medium: "contact".into(),
            third_body: None,
            environment: "authored".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let contact_context = QueryPoint::new()
        .with_quantity("aT", temperature, 300.0)
        .unwrap()
        .with_quantity("bT", temperature, 340.0)
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 17,
        sides: [0, 1].map(|i| ContactCell {
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: curves[i].curve(),
        }),
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &contact_context,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let surfaces = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: cells[i].cell,
        material_id: cells[i].material_id,
        curve: curves[i].curve(),
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let h = [400.0, 560.0];
    let enclosure = radiative::enclosure([0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original AC thermal surface",
            cells[i].point,
            curves[i].curve().state_at_specific_enthalpy(h[i]).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    }));
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surfaces,
        power_tolerance_w: 1e-9,
        radiosity_tolerance_w_m2: 1e-9,
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let mut events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cells[i].cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: cells[i].material_id,
                mass_kg: 1.0,
                specific_h_j_kg: h[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [1.0 + i as f64, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 17,
            a: BODY,
            b: SINK,
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(9), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(10), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 43,
                kernel_id: 97,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &heat, thermal, &cards[0])
    })
}

fn ac_controls() -> SpatialHarmonicInductionControls {
    SpatialHarmonicInductionControls {
        relative_tolerance: 1e-12,
        absolute_residual_a: 1e-12,
        restart: 32,
        maximum_cycles: 100,
        maximum_krylov_scalars: 100000,
        maximum_field_scalars: 100000,
        current_tolerance_a: 1e-9,
        flux_tolerance_wb: 1e-10,
        power_tolerance_w: 1e-9,
        maximum_cell_magnetic_peclet: 100.0,
        maximum_cell_skin_ratio: 1000.0,
    }
}
fn drive() -> [[C64; 3]; 6] {
    [[C64::new(0.8 * PEAK, -0.6 * PEAK); 3]; 6]
}
fn advance<'t, A: HeatAssembly>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    template: &SpatialInduction<'_>,
    controls: CircuitEnthalpySolve,
) -> Result<
    SpatialHarmonicThermalCandidate<'t, A::Sample>,
    CircuitEnthalpyError<SpatialInductionHeatError<A::Error>>,
> {
    heat_assembly_candidate(
        heat,
        cx,
        thermal,
        target,
        template,
        &drive(),
        OMEGA,
        &regions(),
        1e-10,
        controls,
        ac_controls(),
        assembly(),
        &[0.0; 2],
    )
}

// Independent one-mode integration gives K=4/mu, M=3*sigma/40 and F=J/2.
// This reference uses no field, contact, radiation or enthalpy solver.
fn power(temperature: f64) -> f64 {
    let sigma = 1.0 / (RHO * (1.0 + ALPHA * (temperature - 300.0)));
    let m = 3.0 * sigma / 40.0;
    0.5 * m * OMEGA * OMEGA * (0.5 * PEAK).powi(2) / (16.0 + (OMEGA * m).powi(2))
}
fn transfer(ta: f64, tb: f64, resistance: f64) -> f64 {
    let conductance = 1.0 / (1.0 + resistance * (1.0 + 0.001 * (ta - 300.0)));
    let radiation = 5.670_374_419e-8 * (tb.powi(4) - ta.powi(4))
        / (1.0 / EMISSIVITY[0] + 1.0 / EMISSIVITY[1] - 1.0);
    conductance * (tb - ta) + radiation
}
fn root(mut a: f64, mut b: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(a) < 0.0 && f(b) > 0.0);
    for _ in 0..100 {
        let c = 0.5 * (a + b);
        if f(c) > 0.0 {
            b = c;
        } else {
            a = c;
        }
    }
    0.5 * (a + b)
}
fn reference(old: [f64; 2], resistance: f64, frozen_power: Option<f64>) -> [f64; 2] {
    let sink = |ta| {
        root(200.0, 1000.0, |tb| {
            4.0 * (tb - old[1]) + DT * transfer(ta, tb, resistance)
        })
    };
    let ta = root(200.0, 1000.0, |ta| {
        4.0 * (ta - old[0])
            - DT * (frozen_power.unwrap_or_else(|| power(ta)) + transfer(ta, sink(ta), resistance))
    });
    [ta, sink(ta)]
}

#[test]
fn settled_ac_heating_refreshes_conductivity_and_complete_heat_across_steps() {
    let mut outcomes = Vec::new();
    for resistance in [0.5, 8.0] {
        outcomes.push(with_case(
            false,
            resistance,
            |cx, _, heat, mut thermal, card| {
                let prepared = heat.prepare(cx, &thermal).unwrap();
                let problem = problem(card, BASIS);
                let template = initial(&problem, BASIS, 0.7, 1.0);
                let old_a = template.potential_wb().to_vec();
                let old_sources = template.conductors().to_vec();
                let mut temperatures = [300.0, 340.0];
                for interval in 1..=2 {
                    let expected = reference(temperatures, resistance, None);
                    let frozen = reference(temperatures, resistance, Some(power(temperatures[0])));
                    assert!((expected[0] - frozen[0]).abs() > 0.01);
                    let old_cells = thermal.cells().clone();
                    let next = advance(
                        &prepared,
                        cx,
                        &thermal,
                        f64::from(interval) * DT,
                        &template,
                        thermal_controls(),
                    )
                    .unwrap();
                    for i in 0..2 {
                        near(
                            next.thermal().temperature_k([BODY, SINK][i]).unwrap(),
                            expected[i],
                            1e-8,
                        );
                    }
                    let report = next.report();
                    let field = &report.electrical;
                    let expected_power = power(expected[0]);
                    near(field.balance.mean_ohmic_power_w, expected_power, 2e-9);
                    near(field.balance.mean_supplied_power_w, expected_power, 2e-9);
                    near(field.balance.power_residual_w, 0.0, 1e-9);
                    near(field.balance.reactive_power_residual_var, 0.0, 1e-9);
                    near(field.balance.mean_mechanical_power_w, 0.0, 1e-15);
                    near(field.heat_j, DT * expected_power, 5e-9);
                    near(report.deposited_joule_j, field.heat_j, 5e-9);
                    near(report.non_electrical_heat_j, 0.0, 2e-8);
                    near(
                        next.thermal().energy_j() - thermal.energy_j(),
                        field.supplied_work_j,
                        2e-8,
                    );
                    near(report.total_energy_residual_j, 0.0, 2e-8);
                    assert!(field.field.converged_euclidean());
                    assert!(report.circuit_trials > 1);
                    assert_eq!(field.cells.len(), 6);
                    for cell in &field.cells {
                        assert_eq!(cell.cell, BODY);
                        near(cell.mean_ohmic_power_w, expected_power / 6.0, 2e-9);
                    }
                    let contact = &report.source.network.thermal_contacts[0];
                    near(
                        contact.conductance_w_k,
                        1.0 / (1.0 + resistance * (1.0 + 0.001 * (expected[0] - 300.0))),
                        1e-12,
                    );
                    let materials: Vec<_> = report
                        .source
                        .magnetic
                        .iter()
                        .map(|s| {
                            near(
                                *s.material.query_point().get("T").unwrap(),
                                expected[0],
                                1e-8,
                            );
                            assert_eq!(
                                s.material.card_identity(),
                                template.conductors()[s.tetrahedron]
                                    .as_ref()
                                    .unwrap()
                                    .material()
                                    .card_identity()
                            );
                            SpatialMagneticMaterial::Linear(s.material.clone())
                        })
                        .collect();
                    let conductors = report
                        .source
                        .magnetic
                        .iter()
                        .map(|s| {
                            Some(
                                template.conductors()[s.tetrahedron]
                                    .as_ref()
                                    .unwrap()
                                    .with_material(s.material.clone())
                                    .unwrap(),
                            )
                        })
                        .collect();
                    let fresh = template
                        .prepare_materials_trial(materials, conductors, assembly(), &mut || false)
                        .unwrap()
                        .into_harmonic(&mut || false)
                        .unwrap();
                    let response = fresh
                        .solve_harmonic(&drive(), &[[0.0; 3]; 8], OMEGA, ac_controls(), &mut || {
                            false
                        })
                        .unwrap();
                    near(
                        response.balance().mean_magnetic_energy_j,
                        field.balance.mean_magnetic_energy_j,
                        1e-12,
                    );
                    for (native, deposited) in response.cells().iter().zip(&field.cells) {
                        assert_eq!(native.cell, deposited.tetrahedron);
                        near(
                            native.mean_ohmic_power_w,
                            deposited.mean_ohmic_power_w,
                            1e-12,
                        );
                    }
                    assert_eq!(thermal.cells(), &old_cells);
                    assert_eq!(template.time_s(), 0.0);
                    assert_eq!(template.potential_wb(), old_a);
                    assert_eq!(template.conductors(), old_sources);
                    thermal = next.into_parts().0;
                    temperatures = expected;
                }
                temperatures
            },
        ));
    }
    for i in 0..2 {
        assert!((outcomes[0][i] - outcomes[1][i]).abs() > 0.01);
    }
}

#[test]
fn temperature_dependent_permeability_and_foreign_sources_refuse_without_publication() {
    // The curve equals the starting permeability. Source dependence must still
    // refuse before a coincident endpoint could hide magnetic parameter work.
    with_case(true, 0.5, |cx, _, heat, thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let template = initial(&problem, BASIS, 0.0, 1.0);
        let old_cells = thermal.cells().clone();
        assert!(matches!(
            advance(&prepared, cx, &thermal, DT, &template, thermal_controls()),
            Err(CircuitEnthalpyError::Source(
                SpatialInductionHeatError::Native(MagneticThermalError::Invalid(
                    "settled AC permeability must be independent of temperature"
                ))
            ))
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(template.time_s(), 0.0);
    });
    with_case(false, 0.5, |cx, _, heat, thermal, _| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let foreign = card(false, 0.3, EMISSIVITY[0]);
        let problem = problem(&foreign, BASIS);
        let template = initial(&problem, BASIS, 0.0, 1.0);
        assert!(advance(&prepared, cx, &thermal, DT, &template, thermal_controls()).is_err());
        assert_eq!(thermal.time_s(), 0.0);
    });
}

struct TaggedSample<S> {
    heat: S,
    resolution: usize,
}
impl<S: AsRef<HeatNetworkSample>> AsRef<HeatNetworkSample> for TaggedSample<S> {
    fn as_ref(&self) -> &HeatNetworkSample {
        self.heat.as_ref()
    }
}
#[derive(Debug)]
enum TaggedError<E> {
    Heat(E),
    Refused,
}
impl<E: From<HeatNetworkError>> From<HeatNetworkError> for TaggedError<E> {
    fn from(error: HeatNetworkError) -> Self {
        Self::Heat(error.into())
    }
}
struct Tagged<'a, A> {
    inner: &'a A,
    gate: &'a CancelGate,
    calls: Cell<usize>,
    fail_at: Cell<Option<usize>>,
    cancel_at: Cell<Option<usize>>,
}
impl<A: HeatAssembly> HeatAssembly for Tagged<'_, A> {
    type Sample = TaggedSample<A::Sample>;
    type Error = TaggedError<A::Error>;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.inner.prepared_heat_network()
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        let resolution = self.calls.get() + 1;
        self.calls.set(resolution);
        let result = self
            .inner
            .resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(TaggedError::Heat))?;
        if self.fail_at.get() == Some(resolution) {
            return Err(CircuitEnthalpyError::Source(TaggedError::Refused));
        }
        if self.cancel_at.get() == Some(resolution) {
            self.gate.request();
        }
        Ok(ContactThermalCoefficients {
            contacts: result.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: result.bulk.resistors,
                sources_w: result.bulk.sources_w,
                conductances_w_k: result.bulk.conductances_w_k,
                source: TaggedSample {
                    heat: result.bulk.source,
                    resolution,
                },
            },
        })
    }
}

#[test]
fn late_heat_refusal_budget_cancellation_and_retry_preserve_original_histories() {
    with_case(false, 0.5, |cx, gate, heat, thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let template = initial(&problem, BASIS, 0.7, 1.0);
        let old_cells = thermal.cells().clone();
        let old_a = template.potential_wb().to_vec();
        let tagged = Tagged {
            inner: &prepared,
            gate,
            calls: Cell::new(0),
            fail_at: Cell::new(None),
            cancel_at: Cell::new(None),
        };
        let baseline = advance(&tagged, cx, &thermal, DT, &template, thermal_controls()).unwrap();
        assert!(baseline.report().source.network.resolution > 2);
        assert_eq!(
            baseline.report().source.network.resolution,
            tagged.calls.get()
        );
        tagged.calls.set(0);
        tagged.fail_at.set(Some(3));
        assert!(matches!(
            advance(&tagged, cx, &thermal, DT, &template, thermal_controls()),
            Err(CircuitEnthalpyError::Source(
                SpatialInductionHeatError::Heat(TaggedError::Refused)
            ))
        ));
        tagged.calls.set(0);
        tagged.fail_at.set(None);
        let retry = advance(&tagged, cx, &thermal, DT, &template, thermal_controls()).unwrap();
        assert_eq!(retry.thermal().cells(), baseline.thermal().cells());
        assert_eq!(
            retry.report().electrical.heat_j.to_bits(),
            baseline.report().electrical.heat_j.to_bits()
        );
        let mut budget = thermal_controls();
        budget.maximum_circuit_trials = 1;
        assert!(matches!(
            advance(&prepared, cx, &thermal, DT, &template, budget),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        tagged.calls.set(0);
        tagged.cancel_at.set(Some(4));
        assert!(matches!(
            advance(&tagged, cx, &thermal, DT, &template, thermal_controls()),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(template.potential_wb(), old_a);
        assert_eq!(template.time_s(), 0.0);
    });
}
