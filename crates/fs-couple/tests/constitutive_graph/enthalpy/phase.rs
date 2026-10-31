//! G1/G3/G4: actual latent state closes the original material/electrical graph.
use super::*;
use fs_material::conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis};
use fs_qty::QuantitySpec;

pub(super) fn phase_material(
    basis: ElectricalPropertyBasis,
    slope: f64,
    surface: bool,
    reverse_axes: bool,
) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let dims = basis.dimensions();
    let value = if surface {
        let coordinates = [[300.0, 0.0], [600.0, 0.0], [300.0, 1.0], [600.0, 1.0]];
        PropertyValue::Surface2D {
            axes: if reverse_axes {
                ["f".into(), "T".into()]
            } else {
                ["T".into(), "f".into()]
            },
            axis_dims: if reverse_axes {
                [Dims::NONE, Temperature::DIMS]
            } else {
                [Temperature::DIMS, Dims::NONE]
            },
            samples: coordinates
                .into_iter()
                .map(|p| {
                    (
                        if reverse_axes { [p[1], p[0]] } else { p },
                        2.0 + slope * p[1] + 0.01 * (p[0] - 400.0),
                    )
                })
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims,
        }
    } else if slope == 0.0 {
        PropertyValue::Scalar { value: 2.0, dims }
    } else {
        PropertyValue::Curve {
            abscissa: "f".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 2.0), (1.0, 2.0 + slope)],
            dims,
        }
    };
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (PropertyKey::new(basis.property_name(), dims), value),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            PropertyValue::Scalar {
                value: 1.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
        ),
    ] {
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            600.0,
        );
        if key.name() == basis.property_name() && slope != 0.0 {
            validity = validity.with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
        }
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured original phase transport".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "confined solid-liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (300.0, 600.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap();
    (card, curve)
}

pub(super) fn phase_node(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
    maximum_fraction: f64,
) -> MaterialOhmicNode {
    let point = point(400.0)
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.2)
        .unwrap();
    let source = resolve_conductor_state_with_basis(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        basis,
    )
    .unwrap()
    .with_liquid_mass_fraction_axis("f")
    .unwrap();
    MaterialOhmicNode::try_new_with_liquid_mass_fraction(
        card,
        source,
        "T",
        [300.0, 600.0],
        [0.0, maximum_fraction],
        1.0,
        1.0,
    )
    .unwrap()
}

pub(super) fn phase_graph(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
    max: f64,
) -> ConstitutiveGraph {
    let mut graph = ConstitutiveGraph::new();
    graph
        .add_node("resistor", Box::new(phase_node(card, basis, max)))
        .unwrap();
    graph
}

pub(super) fn phase_feedback(current: Option<f64>, voltage: f64) -> Feedback {
    let mut feedback = voltage_feedback();
    feedback.temperatures[0].coordinate = ThermalCoordinate::ChangeFrom(400.0);
    feedback.equations[0].target_si = voltage;
    if let Some(current) = current {
        feedback.fixed.insert(key("resistor", "current"), current);
        feedback.unknowns.clear();
        feedback.equations.clear();
    }
    feedback
}

fn phase_input<'a>(
    graph: &'a ConstitutiveGraph,
    feedback: &'a Feedback,
    fractions: &'a [GraphLiquidMassFraction<'a>],
) -> GraphThermalInput<'a> {
    let mut input = feedback.input(graph, &[]);
    input.liquid_mass_fractions = fractions;
    if feedback.unknowns.is_empty() {
        input.initial_coordinates = &[];
    }
    input
}

fn fraction_binding(curve: &EquilibriumEnthalpyPhaseCurve) -> [GraphLiquidMassFraction<'_>; 1] {
    [GraphLiquidMassFraction {
        input: key("resistor", "liquid_mass_fraction_change"),
        cell: HEATER,
        curve,
        reference: 0.2,
    }]
}

#[test]
fn original_fraction_changes_solved_resistance_current_and_latent_heat() {
    with_cx(|cx, _| {
        for basis in [
            ElectricalPropertyBasis::Resistivity,
            ElectricalPropertyBasis::Conductivity,
        ] {
            for (slope, surface, reversed) in [
                (0.0, false, false),
                (2.0, false, false),
                (6.0, true, false),
                (6.0, true, true),
            ] {
                let (card, curve) = phase_material(basis, slope, surface, reversed);
                let old = timeline(
                    &[(MaterialId(0), &curve)],
                    vec![birth(HEATER, MaterialId(0), 20.0)],
                );
                let graph = phase_graph(&card, basis, 1.0);
                let bindings = fraction_binding(&curve);
                let current = if reversed { -2.0 } else { 2.0 };
                let fixed_current = phase_feedback(Some(current), 0.0);
                let dt = 0.5;
                let candidate = step_candidate(
                    cx,
                    &old,
                    dt,
                    phase_input(&graph, &fixed_current, &bindings),
                    thermal_controls(),
                    &[0.0],
                    no_additional_heat,
                )
                .unwrap();
                let a = 2.0 + 0.2 * slope;
                let supplied = dt * 4.0;
                // Independent scalar BE elimination on m=1, latent span=100.
                let q = match basis {
                    ElectricalPropertyBasis::Resistivity => {
                        supplied * a / (1.0 - supplied * slope / 100.0)
                    }
                    ElectricalPropertyBasis::Conductivity => {
                        2.0 * supplied / (a + (a * a + 4.0 * supplied * slope / 100.0).sqrt())
                    }
                };
                let f = (20.0 + q) / 100.0;
                let source = 2.0 + slope * f;
                let resistance = if basis == ElectricalPropertyBasis::Resistivity {
                    source
                } else {
                    1.0 / source
                };
                close(
                    candidate.thermal().cells()[&HEATER].specific_h_j_kg,
                    20.0 + q,
                );
                close(candidate.thermal().temperature_k(HEATER).unwrap(), 400.0);
                close(
                    candidate.graph().external[&key("resistor", "liquid_mass_fraction_change")],
                    f - 0.2,
                );
                close(
                    candidate.graph().graph.outputs[&key("resistor", "voltage")],
                    current * resistance,
                );
                close(candidate.report().graph_heat_j, q);
                close(candidate.report().thermal.external_heat_j, q);
                let fixed_voltage = phase_feedback(None, current * resistance);
                let dual = step_candidate(
                    cx,
                    &old,
                    dt,
                    phase_input(&graph, &fixed_voltage, &bindings),
                    thermal_controls(),
                    &[0.0],
                    no_additional_heat,
                )
                .unwrap();
                close(dual.thermal().cells()[&HEATER].specific_h_j_kg, 20.0 + q);
                close(dual.graph().external[&key("resistor", "current")], current);
                let next = step_candidate(
                    cx,
                    candidate.thermal(),
                    2.0 * dt,
                    phase_input(&graph, &fixed_current, &bindings),
                    thermal_controls(),
                    &[0.0],
                    no_additional_heat,
                )
                .unwrap();
                let q2 = match basis {
                    ElectricalPropertyBasis::Resistivity => {
                        supplied * source / (1.0 - supplied * slope / 100.0)
                    }
                    ElectricalPropertyBasis::Conductivity => {
                        2.0 * supplied
                            / (source + (source * source + 4.0 * supplied * slope / 100.0).sqrt())
                    }
                };
                close(
                    next.thermal().cells()[&HEATER].specific_h_j_kg,
                    20.0 + q + q2,
                );
                if slope == 0.0 {
                    // Constant source data gives identical old two-input physics.
                    let query = phase_node(&card, basis, 1.0)
                        .source()
                        .original()
                        .material()
                        .query();
                    let old_source = resolve_conductor_state_with_basis(
                        &card,
                        &query,
                        MaterialPropertySelection::SingleClaimOnly,
                        basis,
                    )
                    .unwrap();
                    let old_node = MaterialOhmicNode::try_new(
                        &card,
                        old_source,
                        "T",
                        [300.0, 600.0],
                        1.0,
                        1.0,
                    )
                    .unwrap();
                    let mut old_graph = ConstitutiveGraph::new();
                    old_graph.add_node("resistor", Box::new(old_node)).unwrap();
                    let mut old_input = fixed_current.input(&old_graph, &[]);
                    old_input.initial_coordinates = &[];
                    let legacy = step_candidate(
                        cx,
                        &old,
                        dt,
                        old_input,
                        thermal_controls(),
                        &[0.0],
                        no_additional_heat,
                    )
                    .unwrap();
                    assert_eq!(legacy.thermal().cells(), candidate.thermal().cells());
                }
                assert_eq!(old.cells()[&HEATER].specific_h_j_kg, 20.0);
                assert_eq!(old.time_s(), 0.0);
            }
        }
    });
}

#[test]
fn two_coordinate_surface_tangents_keep_both_basis_chain_rules() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for reverse in [false, true] {
            let (card, _) = phase_material(basis, 3.0, true, reverse);
            let law = phase_node(&card, basis, 1.0);
            let source = 2.0 + 3.0 * 0.6 + 0.01 * 50.0;
            let (r, factor) = if basis == ElectricalPropertyBasis::Resistivity {
                (source, 1.0)
            } else {
                (1.0 / source, -1.0 / (source * source))
            };
            let actual = law.tangent(&[], &[50.0, -2.0, 0.4]).unwrap();
            let expected = [
                -2.0 * 0.01 * factor,
                r,
                -6.0 * factor,
                4.0 * 0.01 * factor,
                -4.0 * r,
                12.0 * factor,
            ];
            for (a, b) in actual.iter().zip(expected) {
                close(*a, b);
            }
            check_consistent_tangent("phase", &law, &[], &[50.0, -2.0, 0.4], 1e-7).unwrap();
        }
    }
}

#[test]
fn phase_domain_foreign_caloric_binding_and_late_cancellation_leave_the_initial_state() {
    let (card, curve) = phase_material(ElectricalPropertyBasis::Resistivity, 2.0, false, false);
    let old = timeline(
        &[(MaterialId(0), &curve)],
        vec![birth(HEATER, MaterialId(0), 20.0)],
    );
    let graph = phase_graph(&card, ElectricalPropertyBasis::Resistivity, 1.0);
    let bounded_graph = phase_graph(&card, ElectricalPropertyBasis::Resistivity, 0.23);
    let feedback = phase_feedback(Some(2.0), 0.0);
    let bindings = fraction_binding(&curve);
    let calls = std::cell::Cell::new(0usize);
    let accepted = with_cx(|cx, _| {
        let foreign = curve.clone();
        let foreign_bindings = fraction_binding(&foreign);
        assert!(matches!(
            step_candidate(
                cx,
                &old,
                0.5,
                phase_input(&graph, &feedback, &foreign_bindings),
                thermal_controls(),
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::Invalid(_))
        ));
        assert!(
            step_candidate(
                cx,
                &old,
                0.5,
                phase_input(&bounded_graph, &feedback, &bindings),
                thermal_controls(),
                &[0.0],
                no_additional_heat
            )
            .is_err()
        );
        step_candidate(
            cx,
            &old,
            0.5,
            phase_input(&graph, &feedback, &bindings),
            thermal_controls(),
            &[0.0],
            |endpoint| {
                calls.set(calls.get() + 1);
                no_additional_heat(endpoint)
            },
        )
        .unwrap()
    });
    assert!(calls.get() > 2);
    with_cx(|cx, gate| {
        let seen = std::cell::Cell::new(0usize);
        assert!(
            step_candidate(
                cx,
                &old,
                0.5,
                phase_input(&graph, &feedback, &bindings),
                thermal_controls(),
                &[0.0],
                |endpoint| {
                    seen.set(seen.get() + 1);
                    if seen.get() == calls.get() - 1 {
                        gate.request();
                    }
                    no_additional_heat(endpoint)
                }
            )
            .is_err()
        );
        assert!(cx.is_cancel_requested());
    });
    let retry = with_cx(|cx, _| {
        step_candidate(
            cx,
            &old,
            0.5,
            phase_input(&graph, &feedback, &bindings),
            thermal_controls(),
            &[0.0],
            no_additional_heat,
        )
        .unwrap()
    });
    assert_eq!(accepted.thermal().cells(), retry.thermal().cells());
    assert_eq!(old.cells()[&HEATER].specific_h_j_kg, 20.0);
    assert_eq!(old.time_s(), 0.0);
}
