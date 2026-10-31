//! Actual material queries and native circuit dynamics; synthetic constitutive
//! data exercises the mathematical coupling, not an experimental material fit.
#[path = "electrical_contacts/source.rs"]
mod contact_source;
#[path = "material_circuit_conductivity/contact_response.rs"]
mod contact_response;
use fs_couple::material_circuit::conductivity::{
    ConductionAxis, ConductionCoordinate, ConductivitySolve, ConductorResponse,
};
use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{
    ELECTRICAL_CONDUCTIVITY_DIMS, ElectricalPropertyBasis, ResolvedConductor,
    ResolvedConductorState,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Dims, QuantitySpec, Temperature};

const E: Dims = Dims([1, 1, -3, 0, -1, 0]);
const J: Dims = Dims([-2, 0, 0, 0, 1, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const AXES: [ConductionAxis<'static>; 2] = [
    ConductionAxis {
        name: "gradient",
        coordinate: ConductionCoordinate::ElectricFieldMagnitude,
    },
    ConductionAxis {
        name: "current_density",
        coordinate: ConductionCoordinate::CurrentDensityMagnitude,
    },
];
const CAPACITANCE: f64 = 0.5;
const DT: f64 = 0.1;

fn point(temperature: f64, field: f64, current: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity("gradient", QuantitySpec::dimensional(E), field)
        .unwrap()
        .with_quantity("current_density", QuantitySpec::dimensional(J), current)
        .unwrap()
}
fn card(
    slope: f64,
    current_dependent: bool,
    maximum_field: f64,
    maximum_density: f64,
) -> MaterialCard {
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            290.0,
            340.0,
        )
        .with_quantity("gradient", QuantitySpec::dimensional(E), 0.0, maximum_field)
        .with_quantity(
            "current_density",
            QuantitySpec::dimensional(J),
            0.0,
            maximum_density,
        );
    let mut claims = ClaimSet::new();
    let (axis, axis_dims, upper) = if current_dependent {
        ("current_density", J, maximum_density)
    } else {
        ("gradient", E, maximum_field)
    };
    for (name, dims, abscissa, abscissa_dims, knots) in [
        (
            "electrical_conductivity",
            ELECTRICAL_CONDUCTIVITY_DIMS,
            axis,
            axis_dims,
            vec![(0.0, 1.0), (upper, 1.0 + slope * upper)],
        ),
        (
            "density",
            DENSITY,
            "T",
            Temperature::DIMS,
            vec![(290.0, 1000.0), (340.0, 1100.0)],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: abscissa.into(),
                    abscissa_dims,
                    knots,
                    dims,
                },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic nonlinear conductor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic electrical material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn conductor(card: &MaterialCard) -> ResolvedConductor {
    let requirements = [
        ElectricalPropertyBasis::Conductivity.requirement().unwrap(),
        ScalarPropertyRequirement::try_new(
            "density",
            DENSITY,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    let bundle = resolve_material_state_point(
        card,
        &point(300.0, 0.0, 0.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedConductorState::try_from_material_with_basis(
        bundle,
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap()
    .into_uniform_conductor(1.0, 1.0)
    .unwrap()
}
fn circuit(card: &MaterialCard, voltage: f64) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: CAPACITANCE,
                    },
                ),
                (1, 0, Branch::Resistor { ohms: 999.0 }),
            ],
            transformers: vec![],
        },
        vec![MaterialResistor {
            branch: 1,
            conductor: conductor(card),
        }],
        vec![],
        vec![CAPACITANCE * voltage],
        &[],
    )
    .unwrap()
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-10,
        power_tolerance_w: 1e-10,
    }
}
fn solve() -> ConductivitySolve {
    ConductivitySolve {
        maximum_iterations: 64,
        absolute_resistance_tolerance_ohm: 1e-13,
        relative_resistance_tolerance: 1e-12,
    }
}
fn bindings(card: &MaterialCard) -> [ConductorResponse<'_>; 1] {
    [ConductorResponse {
        branch: 1,
        card,
        axes: &AXES,
    }]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn constant_conductivity_is_the_native_step_with_actual_field_receipts() {
    let source = card(0.0, false, 8.0, 16.0);
    let circuit = circuit(&source, 2.0);
    let old = circuit
        .step_candidate(circuit.resistors().to_vec(), &[], DT, controls(), || false)
        .unwrap();
    let new = circuit
        .conductivity_candidate(&bindings(&source), &[], DT, controls(), solve(), || false)
        .unwrap();
    assert_eq!(new.circuit().charges_c(), old.circuit().charges_c());
    assert_eq!(new.circuit().fluxes_wb(), old.circuit().fluxes_wb());
    assert_eq!(new.step().electrical.native.x, old.step().native.x);
    assert_eq!(new.step().electrical.heat_j, old.step().heat_j);
    assert_eq!(
        new.step().electrical.native.supplied,
        old.step().native.supplied
    );
    assert_eq!(new.step().iterations, 1);
    assert_eq!(new.step().source_evaluations, 1);
    let evaluated = &new.step().evaluated[0];
    assert_eq!(
        evaluated.conductor.material().query(),
        point(
            300.0,
            old.step().resistors[0].voltage_v.abs(),
            old.step().resistors[0].current_a.abs()
        )
    );
    assert_eq!(evaluated.resistance_defect_ohm, 0.0);
    assert_eq!(evaluated.conductor.material().properties().len(), 2);
}

#[test]
fn data_only_field_law_changes_actual_discharge_and_matches_the_constitutive_quadratic() {
    let mut endpoints = Vec::new();
    for slope in [0.0, 0.2] {
        let source = card(slope, false, 8.0, 16.0);
        let circuit = circuit(&source, 2.0);
        let original = circuit.clone();
        let proposal = circuit
            .conductivity_candidate(&bindings(&source), &[], DT, controls(), solve(), || false)
            .unwrap();
        // q1=q0-dt*I, Vmid=(q0+q1)/(2C), I=(1+slope*Vmid)*Vmid.
        // Stable positive quadratic root, independent of the production solve.
        let a = DT * slope / (2.0 * CAPACITANCE);
        let b = 1.0 + DT / (2.0 * CAPACITANCE);
        let voltage = 4.0 / (b + (b * b + 8.0 * a).sqrt());
        let current = (1.0 + slope * voltage) * voltage;
        let expected_charge = CAPACITANCE * (2.0 * voltage - 2.0);
        let step = proposal.step();
        near(proposal.circuit().charges_c()[0], expected_charge, 2e-11);
        near(step.electrical.resistors[0].voltage_v, voltage, 2e-11);
        near(step.electrical.resistors[0].current_a, current, 2e-11);
        near(step.electrical.heat_j, DT * voltage * current, 2e-11);
        near(
            step.electrical.stored_after_j + step.electrical.heat_j,
            0.5 * CAPACITANCE * 4.0,
            2e-11,
        );
        let sample = &step.evaluated[0];
        assert_eq!(
            sample.conductor.material().query(),
            point(
                300.0,
                sample.point.electric_field_v_m,
                sample.point.current_density_a_m2
            )
        );
        near(
            sample.conductor.resistance_ohm(),
            1.0 / (1.0 + slope * voltage),
            2e-11,
        );
        assert!(sample.power_defect_w.abs() <= controls().power_tolerance_w);
        for (before, after) in original.resistors()[0]
            .conductor
            .material()
            .properties()
            .iter()
            .zip(sample.conductor.material().properties())
        {
            assert_eq!(before.requirement(), after.requirement());
            assert_eq!(
                before.answer().receipt.selected,
                after.answer().receipt.selected
            );
        }
        assert_eq!(
            sample.conductor.material().card_identity(),
            source.content_hash()
        );
        assert_eq!(circuit, original);
        endpoints.push(proposal.circuit().charges_c()[0]);
    }
    assert!(
        endpoints[1] < endpoints[0],
        "the changed material data must drain more actual charge"
    );
}

#[test]
fn current_density_feedback_and_magnitude_axes_preserve_polarity() {
    let slope = 0.05;
    let source = card(slope, true, 8.0, 16.0);
    let k = DT / (2.0 * CAPACITANCE);
    let a = slope * k;
    let b = 1.0 - 2.0 * slope + k;
    let current = 4.0 / (b + (b * b + 8.0 * a).sqrt());
    for sign in [-1.0, 1.0] {
        let circuit = circuit(&source, sign * 2.0);
        let proposal = circuit
            .conductivity_candidate(&bindings(&source), &[], DT, controls(), solve(), || false)
            .unwrap();
        let sample = &proposal.step().evaluated[0];
        near(sample.point.current_a, sign * current, 3e-11);
        near(sample.point.voltage_v, sign * (2.0 - k * current), 3e-11);
        near(
            sample.conductor.material().query().axes()["current_density"],
            current,
            3e-11,
        );
        assert_eq!(
            sample.point.current_a.is_sign_positive(),
            sign.is_sign_positive()
        );
        assert!(proposal.step().electrical.heat_j > 0.0);
    }
}

#[test]
fn refreshed_thermal_context_retains_full_source_bundle_and_fixed_geometry() {
    let source = card(0.2, false, 8.0, 16.0);
    let circuit = circuit(&source, 2.0);
    let material = circuit.resistors()[0]
        .conductor
        .material()
        .at_query_point(&source, &point(320.0, 0.0, 0.0))
        .unwrap();
    let seed = ResolvedConductorState::try_from_material(material)
        .unwrap()
        .into_uniform_conductor(1.0, 1.0)
        .unwrap();
    let proposal = circuit
        .conductivity_with_materials_candidate(
            vec![MaterialResistor {
                branch: 1,
                conductor: seed,
            }],
            &bindings(&source),
            &[],
            DT,
            controls(),
            solve(),
            || false,
        )
        .unwrap();
    let sample = &proposal.step().evaluated[0];
    assert_eq!(sample.conductor.material().query().axes()["T"], 320.0);
    near(
        sample
            .conductor
            .material()
            .property("density")
            .unwrap()
            .value_si(),
        1060.0,
        1e-10,
    );
    assert_eq!(sample.conductor.length_m(), 1.0);
    assert_eq!(sample.conductor.area_m2(), 1.0);
    assert_eq!(
        circuit.resistors()[0].conductor.material().query().axes()["T"],
        300.0
    );
}

#[test]
fn source_domain_nonconvergence_identity_and_late_cancellation_preserve_original() {
    let source = card(0.2, false, 8.0, 16.0);
    let circuit = circuit(&source, 2.0);
    let original = circuit.clone();
    assert!(matches!(
        circuit.conductivity_candidate(
            &bindings(&source),
            &[],
            DT,
            controls(),
            ConductivitySolve {
                maximum_iterations: 1,
                ..solve()
            },
            || false
        ),
        Err(MaterialCircuitError::ConductivityConvergence { iterations: 1, .. })
    ));
    let other = card(0.3, false, 8.0, 16.0);
    assert!(
        circuit
            .conductivity_candidate(&bindings(&other), &[], DT, controls(), solve(), || false)
            .is_err()
    );
    let incomplete_axes = [ConductorResponse {
        branch: 1,
        card: &source,
        axes: &AXES[..1],
    }];
    assert!(
        circuit
            .conductivity_candidate(&incomplete_axes, &[], DT, controls(), solve(), || false)
            .is_err()
    );
    let changed_geometry = ResolvedConductorState::try_from_material(
        circuit.resistors()[0].conductor.material().clone(),
    )
    .unwrap()
    .into_uniform_conductor(2.0, 1.0)
    .unwrap();
    assert!(
        circuit
            .conductivity_with_materials_candidate(
                vec![MaterialResistor {
                    branch: 1,
                    conductor: changed_geometry
                }],
                &bindings(&source),
                &[],
                DT,
                controls(),
                solve(),
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    circuit
        .conductivity_candidate(&bindings(&source), &[], DT, controls(), solve(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        circuit.conductivity_candidate(&bindings(&source), &[], DT, controls(), solve(), || {
            count += 1;
            count == polls
        }),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(circuit, original);
    let limited = card(0.2, false, 1.0, 16.0);
    let limited_circuit = self::circuit(&limited, 2.0);
    let saved = limited_circuit.clone();
    assert!(matches!(
        limited_circuit.conductivity_candidate(
            &bindings(&limited),
            &[],
            DT,
            controls(),
            solve(),
            || false
        ),
        Err(MaterialCircuitError::Conductor(_))
    ));
    assert_eq!(limited_circuit, saved);
}

#[test]
fn nonlinear_bulk_and_series_contacts_use_distinct_voltage_drops_and_heat() {
    let slope = 0.2;
    let source = card(slope, false, 8.0, 16.0);
    let mut charges = Vec::new();
    for resistance in [0.0, 0.5, 2.0] {
        let interface = contact_source::card(&source, &source, resistance, 0.0);
        let contact = contact_source::contact(7, 1, &interface, 1.0);
        for sign in [-1.0, 1.0] {
            let circuit = circuit(&source, sign * 2.0)
                .rebind_contacts_candidate(vec![contact.clone()], &[]).unwrap();
            let saved = circuit.clone();
            let candidate = circuit.conductivity_with_contacts_candidate(
                circuit.resistors().to_vec(), vec![contact.clone()], &bindings(&source),
                &[], DT, controls(), solve(), || false,
            ).unwrap();
            // V0 = E + (Rc + dt/(2C)) I, I = (1 + slope E) E.
            // This analytic quadratic uses BULK E, not the whole branch drop.
            let series = resistance + DT / (2.0 * CAPACITANCE);
            let a = series * slope;
            let b = 1.0 + series;
            let field = 4.0 / (b + (b * b + 8.0 * a).sqrt());
            let current = (1.0 + slope * field) * field;
            let report = candidate.step();
            let bulk = &report.electrical.resistors[0];
            let surface = &report.electrical.contacts[0];
            near(bulk.voltage_v, sign * field, 4e-11);
            near(bulk.current_a, sign * current, 4e-11);
            near(surface.voltage_v, sign * current * resistance, 4e-11);
            near(surface.heat_j, DT * current * current * resistance, 4e-11);
            near(bulk.heat_j, DT * field * current, 4e-11);
            near(report.evaluated[0].point.electric_field_v_m, sign * field, 4e-11);
            near(report.evaluated[0].conductor.material().query().axes()["gradient"], field, 4e-11);
            near(report.electrical.heat_j, bulk.heat_j + surface.heat_j, 4e-11);
            near(report.electrical.stored_after_j + report.electrical.heat_j, 1.0, 4e-11);
            near(candidate.circuit().charges_c()[0], sign * (1.0 - DT * current), 4e-11);
            assert_eq!(surface.contact, contact);
            assert_eq!(circuit, saved);
            if sign > 0.0 { charges.push(candidate.circuit().charges_c()[0]); }
            let next = candidate.circuit();
            let continued = next.conductivity_with_contacts_candidate(
                next.resistors().to_vec(), next.contacts().to_vec(), &bindings(&source),
                &[], DT, controls(), solve(), || false,
            ).unwrap();
            assert!(continued.circuit().charges_c()[0].abs() < next.charges_c()[0].abs());
        }
    }
    assert!(charges.windows(2).all(|pair| pair[1] > pair[0]));
}

#[test]
fn contact_nonlinear_closure_preserves_native_constant_limit_and_refuses_unowned_changes() {
    use fs_couple::material_circuit::contacts::SeriesElectricalContact;
    let source = card(0.0, false, 8.0, 16.0);
    let interface = contact_source::card(&source, &source, 0.5, 0.0);
    let contact = contact_source::contact(7, 1, &interface, 1.0);
    let circuit = circuit(&source, 2.0)
        .rebind_contacts_candidate(vec![contact.clone()], &[]).unwrap();
    let native = circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),
        vec![contact.clone()], &[], DT, controls(), || false).unwrap();
    let run = |contacts, cancel: &mut dyn FnMut() -> bool| circuit.conductivity_with_contacts_candidate(
        circuit.resistors().to_vec(), contacts, &bindings(&source), &[], DT,
        controls(), solve(), cancel,
    );
    let candidate = run(vec![contact.clone()], &mut || false).unwrap();
    assert_eq!(candidate.circuit().charges_c(), native.circuit().charges_c());
    assert_eq!(candidate.step().electrical.resistors, native.step().resistors);
    assert_eq!(candidate.step().electrical.contacts, native.step().contacts);
    assert_eq!(candidate.step().iterations, 1);
    let other_interface = contact_source::card(&source, &source, 0.7, 0.0);
    for contacts in [
        vec![],
        vec![contact.clone(), contact.clone()],
        vec![contact_source::contact(7, 1, &other_interface, 1.0)],
        vec![SeriesElectricalContact::new(7, 1, contact.source().clone(), 2.0).unwrap()],
    ] {
        assert!(run(contacts, &mut || false).is_err());
    }
    let mut polls = 0;
    run(vec![contact.clone()], &mut || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(run(vec![contact], &mut || { count += 1; count == polls }),
        Err(MaterialCircuitError::Cancelled)));
    assert_eq!(circuit.charges_c(), &[1.0]);
}

#[test]
fn contact_electrical_response_and_applicability_cannot_be_frozen_in_bulk_iteration() {
    use fs_couple::material_circuit::contacts::SeriesElectricalContact;
    use fs_matdb::InterfaceSystemCard;
    use fs_material::electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, resolve_electrical_interface_state,
    };
    let source = card(0.2, false, 8.0, 16.0);
    let template = contact_source::card(&source, &source, 0.5, 0.0);
    for dims in [E, J, Dims([2, 1, -3, 0, -1, 0]), Dims([0, 0, 0, 0, 1, 0])] {
        for response in [false, true] {
            let mut claim = template.claims_for(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)[0].1.clone();
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("drive", QuantitySpec::dimensional(dims), 0.0, 10.0);
            let units = claim.key.quantity().dims();
            claim.value = if response { PropertyValue::Curve {
                abscissa: "drive".into(), abscissa_dims: dims,
                knots: vec![(0.0, 0.5), (10.0, 1.0)], dims: units,
            }} else { PropertyValue::Scalar { value: 0.5, dims: units } };
            claim.interpolation = if response { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity };
            let mut claims = ClaimSet::new();
            claims.insert_claim(claim).unwrap();
            let interface = InterfaceSystemCard::assemble(template.surface_a().clone(),
                template.surface_b().clone(), template.context().clone(), claims, vec![]).unwrap();
            let point = QueryPoint::new().with_quantity("drive", QuantitySpec::dimensional(dims), 0.0).unwrap();
            let contact = SeriesElectricalContact::new(7, 1, resolve_electrical_interface_state(
                &interface, &point, MaterialPropertySelection::SingleClaimOnly).unwrap(), 1.0).unwrap();
            let circuit = circuit(&source, 2.0).rebind_contacts_candidate(vec![contact.clone()], &[]).unwrap();
            assert!(matches!(circuit.conductivity_with_contacts_candidate(
                circuit.resistors().to_vec(), vec![contact], &bindings(&source), &[], DT,
                controls(), solve(), || false,
            ), Err(MaterialCircuitError::Invalid("contact electrical operating coordinates require their own nonlinear closure"))));
        }
    }
}
