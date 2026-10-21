//! Simultaneous temperature/field data reaches the existing native circuit solve.
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
    ELECTRICAL_CONDUCTIVITY_DIMS, ElectricalPropertyBasis, ResolvedConductorState,
};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Dims, QuantitySpec, Temperature};

const FIELD: Dims = Dims([1, 1, -3, 0, -1, 0]);
const AXES: [ConductionAxis<'static>; 1] = [ConductionAxis {
    name: "field",
    coordinate: ConductionCoordinate::ElectricFieldMagnitude,
}];

fn source(typed: bool) -> MaterialCard {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        340.0,
    );
    validity = if typed {
        validity.with_quantity("field", QuantitySpec::dimensional(FIELD), 0.0, 8.0)
    } else {
        validity.with("field", 0.0, 8.0)
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("electrical_conductivity", ELECTRICAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "field".into()],
                axis_dims: [Temperature::DIMS, FIELD],
                samples: [[300.0, 0.0], [340.0, 0.0], [300.0, 8.0], [340.0, 8.0]]
                    .into_iter()
                    .map(|at| (at, 1.0 + 0.01 * (at[0] - 300.0) + 0.2 * at[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_CONDUCTIVITY_DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured sigma(T,E); not a calibrated physical material".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "solid".into(),
            process: "analytic numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn circuit(card: &MaterialCard, temperature: f64, typed: bool) -> MaterialCircuit {
    let query = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap();
    let query = if typed {
        query
            .with_quantity("field", QuantitySpec::dimensional(FIELD), 0.0)
            .unwrap()
    } else {
        query.with("field", 0.0).unwrap()
    };
    let material = resolve_material_state_point(
        card,
        &query,
        &[ElectricalPropertyBasis::Conductivity.requirement().unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material(material)
        .unwrap()
        .into_uniform_conductor(1.0, 1.0)
        .unwrap();
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 0.5 }),
                (1, 0, Branch::Resistor { ohms: 999.0 }),
            ],
            transformers: vec![],
        },
        vec![MaterialResistor {
            branch: 1,
            conductor,
        }],
        vec![],
        vec![1.0],
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

#[test]
fn two_coordinate_conductivity_changes_actual_current_and_charge() {
    let source = source(true);
    let mut charge = Vec::new();
    for temperature in [300.0, 320.0] {
        let circuit = circuit(&source, temperature, true);
        let proposal = circuit
            .conductivity_candidate(
                &[ConductorResponse {
                    branch: 1,
                    card: &source,
                    axes: &AXES,
                }],
                &[],
                0.1,
                controls(),
                solve(),
                || false,
            )
            .unwrap();
        // Exact manufactured RC midpoint law:
        // 0.02 V^2 + (1 + 0.1 sigma(T,0)) V - 2 = 0.
        let base = 1.0 + 0.01 * (temperature - 300.0);
        let b = 1.0 + 0.1 * base;
        let voltage = 4.0 / (b + (b * b + 0.16).sqrt());
        let current = (base + 0.2 * voltage) * voltage;
        let branch = &proposal.step().electrical.resistors[0];
        assert!((branch.voltage_v - voltage).abs() < 3e-11);
        assert!((branch.current_a - current).abs() < 3e-11);
        assert!((proposal.circuit().charges_c()[0] - (1.0 - 0.1 * current)).abs() < 3e-11);
        assert!((proposal.step().electrical.heat_j - 0.1 * current * voltage).abs() < 3e-11);
        let resolved = &proposal.step().evaluated[0].conductor;
        assert_eq!(resolved.material().query().axes()["T"], temperature);
        assert!((resolved.material().query().axes()["field"] - voltage).abs() < 3e-11);
        for property in resolved.material().properties() {
            source
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
            assert_eq!(property.answer().receipt.schema_version, 4);
        }
        charge.push(proposal.circuit().charges_c()[0]);
    }
    assert!(
        charge[1] < charge[0],
        "hotter source data must change the actual discharge"
    );
}

#[test]
fn legacy_surface_field_coordinate_cannot_remain_frozen() {
    let source = source(false);
    let circuit = circuit(&source, 300.0, false);
    let original = circuit.clone();
    let result = circuit.conductivity_candidate(
        &[ConductorResponse {
            branch: 1,
            card: &source,
            axes: &[],
        }],
        &[],
        0.1,
        controls(),
        solve(),
        || false,
    );
    assert!(matches!(
        result,
        Err(MaterialCircuitError::Invalid(
            "a source law's electrical coordinate cannot remain frozen"
        ))
    ));
    assert_eq!(circuit, original);
}
