//! Original source resistance and latent conductivity in one native implicit solve.
use fs_conduction::transient::deposition::contact::implicit::{
    PreparedContactError, PreparedMaterialContact,
};
use fs_conduction::transient::deposition::contact::{ContactCell, MaterialContactLink};
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
use fs_material::thermal::{interface::ThermalInterfacePoint, source::ThermalMaterialPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionError,
    ImplicitConductionGraph, ImplicitConductionLaw, ImplicitConductionLinearControls,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::*;

const C: CellId = CellId(17);
const HOT: [CellId; 3] = [CellId(3), CellId(40), CellId(91)];
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured heat-transfer reference".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn card(name: &str, phase_dependent: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(KP, KD),
            value: if phase_dependent {
                PropertyValue::Curve {
                    abscissa: "f".into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, 1.0), (1.0, 2.0)],
                    dims: KD,
                }
            } else {
                PropertyValue::Scalar {
                    value: 2.0,
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
            interpolation: if phase_dependent {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "declared phase chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
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
    let mut point = QueryPoint::new();
    for (axis, value) in [("T", 400.0), ("TA", 400.0), ("TB", 500.0), ("TM", 450.0)] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Temperature::DIMS), value)
            .unwrap();
    }
    for (axis, value) in [("f", 0.0), ("fa", 0.25), ("fb", 1.0)] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}
fn bulk_laws<'a>(
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    cards: &[MaterialCard; 2],
    point: &QueryPoint,
) -> [PreparedPhaseConductivity<'a>; 2] {
    [0, 1].map(|side| {
        PreparedPhaseConductivity::new(
            &curves[side],
            ThermalMaterialPoint {
                material: &cards[side],
                context: point,
                temperature_axis: "T",
                liquid_fraction_axis: Some("f"),
            },
            &PropertyKey::new(KP, KD),
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            [-100.0, 500.0],
            || Ok(()),
        )
        .unwrap()
    })
}
fn initial(curves: &[EquilibriumEnthalpyPhaseCurve; 2], leaves: usize) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for (i, id) in std::iter::once(C)
        .chain(HOT[..leaves].iter().copied())
        .enumerate()
    {
        let mass = if i == 0 { 2.0 } else { 3.0 };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(u64::from(i != 0)),
                mass_kg: mass,
                specific_h_j_kg: if i == 0 { 25.0 } else { 200.0 },
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: mass,
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
                id: 100 + i as u64,
                a: C,
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
            max_substeps: 1000,
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
    fixed: f64,
) -> ConductionObserver<'a> {
    let links = thermal
        .links()
        .keys()
        .map(|&link_id| ConductiveLink {
            link_id,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: fixed,
        })
        .collect::<Vec<_>>();
    ConductionObserver::new(
        thermal,
        &[
            (MaterialId(0), &laws[0] as &dyn ConductivityLaw),
            (MaterialId(1), &laws[1] as &dyn ConductivityLaw),
        ],
        &links,
        &[],
        || Ok(()),
    )
    .unwrap()
}
fn source(
    cards: &[MaterialCard; 2],
    axis: Option<&str>,
    base: f64,
    slope: f64,
    minimum_tb: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let value = if let Some(axis) = axis {
        let fraction = axis == "fa" || axis == "fb";
        let (lo, hi) = if fraction { (0.0, 1.0) } else { (300.0, 800.0) };
        PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: if fraction {
                Dims::NONE
            } else {
                Temperature::DIMS
            },
            knots: vec![(lo, base), (hi, base + slope * (hi - lo))],
            dims: RD,
        }
    } else {
        PropertyValue::Scalar {
            value: base,
            dims: RD,
        }
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value,
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
                    minimum_tb,
                    800.0,
                )
                .with_quantity(
                    "TM",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                .with_quantity("fb", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: if axis.is_some() {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B".into(),
        },
        SystemContext {
            medium: "declared contact".into(),
            third_body: None,
            environment: "closed heat exchange".into(),
            history: "fixed interface".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn binding<'a>(
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    cards: &'a [MaterialCard; 2],
    source: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
    leaf: usize,
) -> MaterialContactLink<'a> {
    MaterialContactLink {
        link: 100 + leaf as u64,
        sides: [0, 1].map(|side| ContactCell {
            cell: if side == 0 { C } else { HOT[leaf] },
            material_id: MaterialId(side as u64),
            curve: &curves[side],
        }),
        point: ThermalInterfacePoint {
            interface: source,
            materials: [&cards[0], &cards[1]],
            context: point,
            temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [Some("fa"), Some("fb")],
            mean_temperature_axis: Some("TM"),
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }
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
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn phase_contact_source_changes_actual_branched_heat_with_both_latent_and_sensible_masses() {
    let mut transfers = Vec::new();
    for (nonlinear_bulk, slope) in [(false, 0.2), (false, 0.8), (true, 0.5)] {
        let cards = [card("melting", nonlinear_bulk), card("hot", false)];
        let curves = cards.each_ref().map(chart);
        let point = context();
        let laws = bulk_laws(&curves, &cards, &point);
        let old = initial(&curves, 3);
        let source = source(&cards, Some("fa"), 0.25, slope, 300.0);
        let contacts = (0..3)
            .map(|i| {
                PreparedMaterialContact::new(
                    &old,
                    binding(&curves, &cards, &source, &point, i),
                    [[-100.0, 500.0]; 2],
                    || Ok(()),
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let bound = contacts
            .iter()
            .enumerate()
            .map(|(i, c)| (100 + i as u64, c as &dyn ImplicitContactLaw))
            .collect::<Vec<_>>();
        let graph = ImplicitConductionGraph::new(
            &old,
            &conduction(&old, &laws, 0.0),
            &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_contacts(&old, &bound, || Ok(()))
        .unwrap();
        let mut work = ImplicitConductionWorkspace::new(4).unwrap();
        let dt = 0.1;
        let (next, audit) = graph
            .advance_to(&old, dt, controls(dt), &mut work, || Ok(()))
            .unwrap();
        let q = 2.0 * (next.cells()[&C].specific_h_j_kg - 25.0);
        let fraction = (25.0 + q / 2.0) / 100.0;
        assert!(fraction > 0.25 && fraction < 1.0);
        let ka = if nonlinear_bulk { 1.0 + fraction } else { 2.0 };
        let g = 1.0 / (0.5 / ka + 0.25 + 0.25 + slope * fraction);
        near(q, 3.0 * dt * g * (100.0 - q / 9.0), 2e-8);
        if !nonlinear_bulk {
            // Independent scalar elimination: B Q²+A Q=3*100*dt.
            // C_cold=latent, m_cold=2; each hot body has m Cp=3 J/K.
            let a = 0.5 + 0.25 + 0.25 * slope + dt / 3.0;
            let b = slope / 200.0;
            let exact = 60.0 / (a + (a * a + 120.0 * b).sqrt());
            near(q, exact, 2e-8);
            transfers.push(q);
        }
        for hot in HOT {
            near(next.cells()[&hot].specific_h_j_kg, 200.0 - q / 9.0, 2e-8);
        }
        near(next.energy_j(), old.energy_j(), 2e-8);
        near(audit.energy.external_heat_j, 0.0, 0.0);
        assert!(audit.entropy.conduction_generation_j_k > 0.0);
        let response = contacts[0]
            .state_at_specific_enthalpies([25.0 + q / 2.0, 200.0 - q / 9.0])
            .unwrap();
        near(response.resistance_m2_k_w, 0.25 + slope * fraction, 1e-14);
        near(
            response.resistance_enthalpy_derivatives_m2_kg_k_w_j[0],
            slope / 100.0,
            1e-14,
        );
        assert_eq!(response.resistance_enthalpy_derivatives_m2_kg_k_w_j[1], 0.0);
        assert_eq!(
            contacts[0].source_line().knots[0].1.receipt.selected,
            contacts[0].reference().interface.receipt().selected
        );
        assert_eq!(old.cells()[&C].specific_h_j_kg, 25.0);
    }
    assert!(transfers[0] > transfers[1]);
}

#[test]
fn constant_source_preserves_fixed_resistance_chain_and_graph_results_exactly() {
    let cards = [card("melting", false), card("hot", false)];
    let curves = cards.each_ref().map(chart);
    let point = context();
    let laws = bulk_laws(&curves, &cards, &point);
    let old = initial(&curves, 1);
    let source = source(&cards, None, 0.75, 0.0, 300.0);
    let contact = PreparedMaterialContact::new(
        &old,
        binding(&curves, &cards, &source, &point, 0),
        [[-100.0, 500.0]; 2],
        || Ok(()),
    )
    .unwrap();
    let bulk = conduction(&old, &laws, 0.0);
    let fixed = conduction(&old, &laws, 0.75);
    let materials = [
        (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
        (MaterialId(1), &laws[1] as &dyn ImplicitConductionLaw),
    ];
    let attached = [(100, &contact as &dyn ImplicitContactLaw)];
    let legacy = ImplicitConductionChain::new(&old, &fixed, &materials, || Ok(())).unwrap();
    let coupled = ImplicitConductionChain::new(&old, &bulk, &materials, || Ok(()))
        .unwrap()
        .with_contacts(&old, &attached, || Ok(()))
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &old,
        &bulk,
        &materials,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_contacts(&old, &attached, || Ok(()))
    .unwrap();
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (a, _) = legacy
        .advance_to(&old, 0.2, controls(0.2), &mut work, || Ok(()))
        .unwrap();
    let (b, _) = coupled
        .advance_to(&old, 0.2, controls(0.2), &mut work, || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    let (c, _) = graph
        .advance_to(&old, 0.2, controls(0.2), &mut work, || Ok(()))
        .unwrap();
    near(
        c.cells()[&C].specific_h_j_kg,
        b.cells()[&C].specific_h_j_kg,
        2e-8,
    );
    assert!(
        ImplicitConductionChain::new(&old, &fixed, &materials, || Ok(()))
            .unwrap()
            .with_contacts(&old, &attached, || Ok(()))
            .is_err()
    );
    assert!(
        ImplicitConductionChain::new(&old, &bulk, &materials, || Ok(()))
            .unwrap()
            .with_contacts(&old, &[attached[0], attached[0]], || Ok(()))
            .is_err()
    );
}

#[test]
fn independent_temperature_and_phase_axes_keep_their_complete_enthalpy_chain_rule() {
    let cards = [card("melting", false), card("hot", false)];
    let curves = cards.each_ref().map(chart);
    let point = context();
    let old = initial(&curves, 1);
    for (axis, h, expected) in [
        ("TA", [150.0, 220.0], [0.002, 0.0]),
        ("TB", [150.0, 220.0], [0.0, 0.002]),
        ("TM", [150.0, 220.0], [0.001, 0.001]),
        ("fa", [25.0, 50.0], [0.00002, 0.0]),
        ("fb", [25.0, 50.0], [0.0, 0.00002]),
    ] {
        let source = source(&cards, Some(axis), 0.25, 0.002, 300.0);
        let contact = PreparedMaterialContact::new(
            &old,
            binding(&curves, &cards, &source, &point, 0),
            [[-100.0, 500.0]; 2],
            || Ok(()),
        )
        .unwrap();
        let actual = contact.state_at_specific_enthalpies(h).unwrap();
        for side in 0..2 {
            near(
                actual.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
                expected[side],
                1e-14,
            );
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.125;
            minus[side] -= 0.125;
            let difference = (contact
                .state_at_specific_enthalpies(plus)
                .unwrap()
                .resistance_m2_k_w
                - contact
                    .state_at_specific_enthalpies(minus)
                    .unwrap()
                    .resistance_m2_k_w)
                / 0.25;
            near(difference, expected[side], 2e-14);
        }
        assert!(std::ptr::eq(contact.caloric_curves()[0], &curves[0]));
    }
}

#[test]
fn final_contact_domain_exit_and_late_cancellation_preserve_both_original_histories() {
    let cards = [card("melting", true), card("hot", false)];
    let curves = cards.each_ref().map(chart);
    let point = context();
    let old = initial(&curves, 3);
    let laws = bulk_laws(&curves, &cards, &point);
    let source = source(&cards, Some("fa"), 0.25, 0.5, 499.0);
    let contacts = (0..3)
        .map(|i| {
            PreparedMaterialContact::new(
                &old,
                binding(&curves, &cards, &source, &point, i),
                [[-100.0, 500.0]; 2],
                || Ok(()),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    assert!(matches!(
        contacts[0].state_at_specific_enthalpies([25.0, 198.0]),
        Err(PreparedContactError::OutsideSourceDomain { .. })
    ));
    let bound = contacts
        .iter()
        .enumerate()
        .map(|(i, c)| (100 + i as u64, c as &dyn ImplicitContactLaw))
        .collect::<Vec<_>>();
    let graph = ImplicitConductionGraph::new(
        &old,
        &conduction(&old, &laws, 0.0),
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_contacts(&old, &bound, || Ok(()))
    .unwrap();
    let mut work = ImplicitConductionWorkspace::new(4).unwrap();
    let limited = ImplicitConductionControls {
        minimum_dt_s: 1e-4,
        maximum_halvings: 8,
        maximum_steps: 64,
        ..controls(0.1)
    };
    assert!(
        graph
            .advance_to(&old, 0.1, limited, &mut work, || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    let (short, _) = graph
        .advance_to(&old, 0.001, controls(0.001), &mut work, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut retry = 0;
    assert!(matches!(
        graph.advance_to(&old, 0.001, controls(0.001), &mut work, || {
            retry += 1;
            if retry == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    let (retried, _) = graph
        .advance_to(&old, 0.001, controls(0.001), &mut work, || Ok(()))
        .unwrap();
    assert_eq!(short.cells(), retried.cells());
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&C].specific_h_j_kg, 25.0);
    for contact in &contacts {
        assert_eq!(contact.reference().phases[0].specific_enthalpy_j_kg(), 25.0);
    }
}
