//! G1/G3/G4: original phase conductivity reaches native diffusion and contact heat.
use fs_conduction::transient::deposition::contact::{
    ContactCell, MaterialContactLink, MaterialContactObserver, MaterialContactObserverError,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{
    THERMAL_CONDUCTIVITY_DIMS as KD, THERMAL_CONDUCTIVITY_PROPERTY as KP,
};
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionError,
    ImplicitConductionGraph, ImplicitConductionLaw, ImplicitConductionLinearControls,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::*;

const CENTER: CellId = CellId(17);
const HOT: [CellId; 3] = [CellId(3), CellId(40), CellId(91)];

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic phase-conduction oracle, not measured data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

fn card(name: &str, scale: f64, phase_dependent: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(KP, KD),
            value: if phase_dependent {
                PropertyValue::Curve {
                    abscissa: "f".into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, scale), (1.0, 2.0 * scale)],
                    dims: KD,
                }
            } else {
                PropertyValue::Scalar {
                    value: scale,
                    dims: KD,
                }
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if phase_dependent {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "declared solid/liquid chart".into(),
            process: "synthetic uniform cell".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (500.0, 800.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
}

fn material_law<'a>(
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    card: &MaterialCard,
    context: &QueryPoint,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new(KP, KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

// One finite melting body exchanges heat with one or three finite hot bodies.
// The three-leaf topology is a real branch, not a reordered direct chain.
fn initial<'a>(
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    leaves: usize,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, id) in std::iter::once(CENTER)
        .chain(HOT[..leaves].iter().copied())
        .enumerate()
    {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(u64::from(i != 0)),
                mass_kg: 1.0,
                specific_h_j_kg: if i == 0 { 25.0 } else { 200.0 },
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for (i, &hot) in HOT[..leaves].iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 500 + i as u64,
                a: CENTER,
                b: hot,
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: leaves + 1,
            max_links: leaves,
            max_events: 2 * leaves + 1,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn conduction<'a>(
    thermal: &DepositionSimulation<'_>,
    laws: &'a [PreparedPhaseConductivity<'_>; 2],
    leaves: usize,
) -> ConductionObserver<'a> {
    let geometry = (0..leaves)
        .map(|i| ConductiveLink {
            link_id: 500 + i as u64,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        })
        .collect::<Vec<_>>();
    ConductionObserver::new(
        thermal,
        &[
            (MaterialId(0), &laws[0] as &dyn ConductivityLaw),
            (MaterialId(1), &laws[1] as &dyn ConductivityLaw),
        ],
        &geometry,
        &[],
        || Ok(()),
    )
    .unwrap()
}

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..Default::default()
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.15e} != {b:.15e}");
}

// Eliminate every finite hot-body endpoint from the backward-Euler equations.
// While the middle body melts, k=s*(1+(25+Q)/100), G=4k/(k+2), and
// Q=leaves*100*dt*G/(1+dt*G). This independent quadratic has one positive root.
fn plateau_heat(scale: f64, leaves: usize, dt: f64) -> f64 {
    let a = 1.25 * scale;
    let b = 0.01 * scale;
    let c = 4.0 * dt * 100.0 * leaves as f64;
    let aa = (1.0 + 4.0 * dt) * b;
    let bb = (1.0 + 4.0 * dt) * a + 2.0 - c * b;
    let cc = c * a;
    2.0 * cc / (bb + (bb * bb + 4.0 * aa * cc).sqrt())
}

#[test]
fn original_phase_sources_drive_branched_latent_heat_and_material_only_front_changes() {
    let mut front_h = Vec::new();
    for scale in [0.5, 2.0] {
        let cards = [
            card("melting patch", scale, true),
            card("finite hot body", 2.0, false),
        ];
        let curves = [caloric(&cards[0]), caloric(&cards[1])];
        let point = context();
        let laws = [
            material_law(&curves[0], &cards[0], &point, [-100.0, 500.0]),
            material_law(&curves[1], &cards[1], &point, [-100.0, 500.0]),
        ];
        let old = initial(&curves, 3);
        let bulk = conduction(&old, &laws, 3);
        let bindings = [
            (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &laws[1] as &dyn ImplicitConductionLaw),
        ];
        assert!(ImplicitConductionChain::new(&old, &bulk, &bindings, || Ok(())).is_err());
        let graph = ImplicitConductionGraph::new(
            &old,
            &bulk,
            &bindings,
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap();
        let mut work = ImplicitConductionWorkspace::new(4).unwrap();
        let (next, report) = graph
            .advance_to(&old, 0.1, controls(0.1), &mut work, || Ok(()))
            .unwrap();
        let q = plateau_heat(scale, 3, 0.1);
        assert!(q > 0.0 && q < 75.0);
        near(next.cells()[&CENTER].specific_h_j_kg, 25.0 + q, 2e-8);
        near(next.temperature_k(CENTER).unwrap(), 400.0, 0.0);
        for hot in HOT {
            near(next.cells()[&hot].specific_h_j_kg, 200.0 - q / 3.0, 2e-8);
        }
        let state = laws[0].state_at_specific_enthalpy(25.0 + q).unwrap();
        near(state.liquid_mass_fraction(), (25.0 + q) / 100.0, 1e-14);
        near(state.temperature_enthalpy_derivative_kg_k_j(), 0.0, 0.0);
        near(
            state.conductivity_enthalpy_derivative_kg_w_m_k_j(),
            scale / 100.0,
            1e-14,
        );
        near(next.energy_j(), old.energy_j(), 2e-8);
        near(report.energy.external_heat_j, 0.0, 0.0);
        assert_eq!(report.energy.substeps, 1);
        assert!(report.linear_iterations > 0);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert!(std::ptr::eq(laws[0].caloric_curve(), &curves[0]));

        let (front, _) = graph
            .advance_to(&old, 0.18, controls(0.005), &mut work, || Ok(()))
            .unwrap();
        front_h.push(front.cells()[&CENTER].specific_h_j_kg);
        near(front.energy_j(), old.energy_j(), 2e-7);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CENTER].specific_h_j_kg, 25.0);
    }
    // Bounds are independent of the nonlinear implementation: before melt
    // completion slow G<=4/3 and fast G>=20/9. Even their BE comparison bounds
    // place Q at t=.18 below/above the 75 J still needed for this phase front.
    assert!(
        front_h[0] < 100.0 && front_h[1] > 100.0,
        "material-only front states {front_h:?}"
    );
}

#[test]
fn direct_chain_and_graph_share_the_original_caloric_object_and_nonlinear_plateau() {
    let cards = [
        card("melting patch", 1.0, true),
        card("finite hot body", 2.0, false),
    ];
    let curves = [caloric(&cards[0]), caloric(&cards[1])];
    let point = context();
    let laws = [
        material_law(&curves[0], &cards[0], &point, [-100.0, 500.0]),
        material_law(&curves[1], &cards[1], &point, [-100.0, 500.0]),
    ];
    let old = initial(&curves, 1);
    let bulk = conduction(&old, &laws, 1);
    let bindings = [
        (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
        (MaterialId(1), &laws[1] as &dyn ImplicitConductionLaw),
    ];
    let chain = ImplicitConductionChain::new(&old, &bulk, &bindings, || Ok(())).unwrap();
    let graph = ImplicitConductionGraph::new(
        &old,
        &bulk,
        &bindings,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (direct, report) = chain
        .advance_to(&old, 0.1, controls(0.1), &mut work, || Ok(()))
        .unwrap();
    let (sparse, _) = graph
        .advance_to(&old, 0.1, controls(0.1), &mut work, || Ok(()))
        .unwrap();
    let q = plateau_heat(1.0, 1, 0.1);
    near(direct.cells()[&CENTER].specific_h_j_kg, 25.0 + q, 2e-8);
    near(direct.cells()[&HOT[0]].specific_h_j_kg, 200.0 - q, 2e-8);
    near(
        direct.cells()[&CENTER].specific_h_j_kg,
        sparse.cells()[&CENTER].specific_h_j_kg,
        2e-8,
    );
    assert_eq!(report.linear_iterations, 0);
    let copied = curves[0].clone();
    let foreign = material_law(&copied, &cards[0], &point, [-100.0, 500.0]);
    let wrong = [
        (MaterialId(0), &foreign as &dyn ImplicitConductionLaw),
        bindings[1],
    ];
    assert!(
        ImplicitConductionGraph::new(
            &old,
            &bulk,
            &wrong,
            ImplicitConductionLinearControls::default(),
            || Ok(())
        )
        .is_err()
    );
}

fn contact_source(
    cards: &[MaterialCard; 2],
    scale: f64,
    maximum_fraction: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Curve {
                abscissa: "fa".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, scale), (1.0, 2.0 * scale)],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "TA",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity(
                    "TB",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity(
                    "fa",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                    maximum_fraction,
                )
                .with_quantity("fb", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "patch".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "hot body".into(),
        },
        SystemContext {
            medium: "declared conductive interface".into(),
            third_body: None,
            environment: "insulated exterior".into(),
            history: "fixed area".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn contact_point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("TA", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
        .unwrap()
        .with_quantity("TB", QuantitySpec::dimensional(Temperature::DIMS), 500.0)
        .unwrap()
        .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.25)
        .unwrap()
        .with_quantity("fb", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
}

fn contact<'a>(
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    cards: &'a [MaterialCard; 2],
    source: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
) -> MaterialContactLink<'a> {
    MaterialContactLink {
        link: 500,
        sides: [0, 1].map(|i| ContactCell {
            cell: if i == 0 { CENTER } else { HOT[0] },
            material_id: MaterialId(i as u64),
            curve: &curves[i],
        }),
        point: ThermalInterfacePoint {
            interface: source,
            materials: [&cards[0], &cards[1]],
            context: point,
            temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [Some("fa"), Some("fb")],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }
}

#[test]
fn phase_conductivity_and_phase_contact_resistance_drive_one_native_finite_body_exchange() {
    let cards = [
        card("melting patch", 1.0, true),
        card("finite hot body", 2.0, false),
    ];
    let curves = [caloric(&cards[0]), caloric(&cards[1])];
    let point = context();
    let laws = [
        material_law(&curves[0], &cards[0], &point, [-100.0, 500.0]),
        material_law(&curves[1], &cards[1], &point, [-100.0, 500.0]),
    ];
    let old = initial(&curves, 1);
    let contact_point = contact_point();
    let mut delivered = Vec::new();
    for resistance in [0.2, 0.8] {
        let source = contact_source(&cards, resistance, 1.0);
        let bulk = conduction(&old, &laws, 1);
        let coupled = MaterialContactObserver::new(
            bulk,
            &old,
            vec![contact(&curves, &cards, &source, &contact_point)],
            || Ok(()),
        )
        .unwrap();
        let (next, accepted, audit) = old
            .advance_observed(
                &coupled,
                0.1,
                StepConfig {
                    max_dt_s: 0.1,
                    energy_tolerance_j: 1e-9,
                },
                || Ok(()),
            )
            .unwrap();
        let g = 1.0 / (0.5 / 1.25 + 0.5 / 2.0 + 1.25 * resistance);
        let q = 0.1 * g * 100.0;
        near(next.cells()[&CENTER].specific_h_j_kg, 25.0 + q, 1e-10);
        near(next.cells()[&HOT[0]].specific_h_j_kg, 200.0 - q, 1e-10);
        near(next.energy_j(), old.energy_j(), 1e-10);
        near(audit.external_heat_j, 0.0, 0.0);
        let applied = &accepted.last_samples()[0];
        near(
            applied.bulk_conductance_w_k,
            1.0 / (0.5 / 1.25 + 0.25),
            1e-14,
        );
        near(applied.conductance_w_k, g, 1e-14);
        let endpoint = &accepted.endpoint_sources()[0];
        near(
            endpoint.phases[0].liquid_mass_fraction(),
            (25.0 + q) / 100.0,
            1e-14,
        );
        near(
            endpoint.interface.value_m2_k_per_w(),
            resistance * (1.0 + (25.0 + q) / 100.0),
            1e-14,
        );
        near(
            laws[0].conductivity_w_m_k(25.0 + q).unwrap(),
            1.0 + (25.0 + q) / 100.0,
            1e-14,
        );
        assert_eq!(
            applied.interface.receipt().selected,
            endpoint.interface.receipt().selected
        );
        delivered.push(q);
    }
    assert!(delivered[0] > delivered[1]);
}

#[test]
fn implicit_domain_and_explicit_contact_exits_roll_back_and_short_retry_is_clean() {
    let cards = [
        card("melting patch", 1.0, true),
        card("finite hot body", 2.0, false),
    ];
    let curves = [caloric(&cards[0]), caloric(&cards[1])];
    let point = context();
    let bounded = [
        material_law(&curves[0], &cards[0], &point, [-100.0, 40.0]),
        material_law(&curves[1], &cards[1], &point, [-100.0, 500.0]),
    ];
    let old = initial(&curves, 3);
    let bulk = conduction(&old, &bounded, 3);
    let graph = ImplicitConductionGraph::new(
        &old,
        &bulk,
        &[
            (MaterialId(0), &bounded[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &bounded[1] as &dyn ImplicitConductionLaw),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut work = ImplicitConductionWorkspace::new(4).unwrap();
    let limits = ImplicitConductionControls {
        minimum_dt_s: 1e-4,
        maximum_halvings: 8,
        maximum_steps: 64,
        ..controls(0.1)
    };
    assert!(
        graph
            .advance_to(&old, 0.1, limits, &mut work, || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    let (short, _) = graph
        .advance_to(&old, 0.01, controls(0.01), &mut work, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut cancelled = 0;
    let result = graph.advance_to(&old, 0.01, controls(0.01), &mut work, || {
        cancelled += 1;
        if cancelled == polls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    let (retried, _) = graph
        .advance_to(&old, 0.01, controls(0.01), &mut work, || Ok(()))
        .unwrap();
    assert_eq!(short.cells(), retried.cells());
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CENTER].specific_h_j_kg, 25.0);

    let old = initial(&curves, 1);
    let contact_point = contact_point();
    let source = contact_source(&cards, 0.2, 0.3);
    let coupled = MaterialContactObserver::new(
        conduction(&old, &bounded, 1),
        &old,
        vec![contact(&curves, &cards, &source, &contact_point)],
        || Ok(()),
    )
    .unwrap();
    let config = StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-9,
    };
    assert!(matches!(
        old.advance_observed(&coupled, 0.1, config, || Ok(())),
        Err(MaterialContactObserverError::Contact(_))
    ));
    let (short, accepted, _) = old
        .advance_observed(&coupled, 0.01, config, || Ok(()))
        .unwrap();
    assert!(short.cells()[&CENTER].specific_h_j_kg < 30.0);
    assert!(accepted.endpoint_sources()[0].phases[0].liquid_mass_fraction() < 0.3);
    near(
        coupled.endpoint_sources()[0].phases[0].liquid_mass_fraction(),
        0.25,
        0.0,
    );
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CENTER].specific_h_j_kg, 25.0);
}
