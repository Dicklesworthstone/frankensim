//! Manufactured T/E and T/J sources, independent finite-capacity BE balances.
use super::{cells, config, ends, field_controls, initial, limits, mesh, near, network};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    conduction::{ConductionDrive, ConductionThermalError},
};
use fs_em::conduction::{
    ConductionCurrentTrace,
    response::{ConductionResponseControls, PlanarConductionResponse},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState,
        response::{ConductorMagnitude, ScalarConductorResponse},
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
        ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
    },
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::CellId;

const CP: [f64; 2] = [2.0, 3.0];
const ALPHA: [f64; 2] = [0.002, 0.004];

fn point(t: f64, coordinate: ConductorMagnitude) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap()
}
fn controls() -> ConductionResponseControls {
    ConductionResponseControls {
        maximum_iterations: 100,
        absolute_conductivity_tolerance_s_m: 1e-12,
        relative_conductivity_tolerance: 1e-12,
        current_density_tolerance_a_m2: 1e-10,
    }
}

fn material(
    r: f64,
    cp: f64,
    alpha: f64,
    beta: f64,
    coordinate: ConductorMagnitude,
    field_bounded_heat: bool,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    material_with_caloric(r, cp, alpha, beta, coordinate, field_bounded_heat, None)
}

fn material_with_caloric(
    r: f64,
    cp: f64,
    alpha: f64,
    beta: f64,
    coordinate: ConductorMagnitude,
    field_bounded_heat: bool,
    caloric: Option<(&str, bool)>,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let basis = ElectricalPropertyBasis::Resistivity;
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let k = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let electrical = PropertyKey::new("original_transport", basis.dimensions());
    let mut claims = ClaimSet::new();
    for (key, mut value) in [
        (
            heat.clone(),
            PropertyValue::Scalar {
                value: cp,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
        ),
        (
            density.clone(),
            PropertyValue::Scalar {
                value: 1.0,
                dims: Density::DIMS,
            },
        ),
        (
            k.clone(),
            PropertyValue::Scalar {
                value: 1.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
        ),
        (
            PropertyKey::with_quantity(
                "melting_point",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
            ),
            PropertyValue::Scalar {
                value: 1000.0,
                dims: Temperature::DIMS,
            },
        ),
        (
            PropertyKey::new("latent_heat_fusion", fs_qty::Dims([2, 0, -2, 0, 0, 0])),
            PropertyValue::Scalar {
                value: 100.0,
                dims: fs_qty::Dims([2, 0, -2, 0, 0, 0]),
            },
        ),
        (
            electrical.clone(),
            PropertyValue::Surface2D {
                axes: ["T".into(), "operating".into()],
                axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                samples: [[200.0, 0.0], [2000.0, 0.0], [200.0, 20.0], [2000.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, r + alpha * (p[0] - 300.0) + beta * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
        ),
    ] {
        if caloric.is_some_and(|(name, response)| name == key.name() && response) {
            let base = match &value {
                PropertyValue::Scalar { value, .. } => *value,
                _ => unreachable!("only authored scalar caloric fixture properties are replaced"),
            };
            value = PropertyValue::Surface2D {
                axes: ["T".into(), "operating".into()],
                axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                samples: [[200.0, 0.0], [2000.0, 0.0], [200.0, 20.0], [2000.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, base + 0.01 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: key.quantity().dims(),
            };
        }
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            200.0,
            2000.0,
        );
        if key == electrical
            || field_bounded_heat && key == k
            || caloric.is_some_and(|(name, _)| name == key.name())
        {
            validity = validity.with_quantity("operating", coordinate.quantity(), 0.0, 20.0);
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
                    source: "manufactured nonlinear electrothermal material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0, coordinate),
            upper: &point(2000.0, coordinate),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn problem(cards: &[MaterialCard; 2], coordinate: ConductorMagnitude) -> PlanarConductionResponse {
    problem_at(cards, coordinate, &point(300.0, coordinate))
}
fn problem_at(
    cards: &[MaterialCard; 2],
    coordinate: ConductorMagnitude,
    context: &QueryPoint,
) -> PlanarConductionResponse {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("original_transport", basis.dimensions());
    let sources = [0, 0, 1, 1]
        .into_iter()
        .map(|i| {
            let material = resolve_material_state_point(
                &cards[i],
                context,
                &[
                    basis.requirement_for(&key).unwrap(),
                    ScalarPropertyRequirement::try_new(
                        "density",
                        Density::DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                    ScalarPropertyRequirement::try_new(
                        "specific_heat_capacity",
                        SPECIFIC_HEAT_CAPACITY_DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                ],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let conductor =
                ResolvedConductorState::try_from_material_with_property(material, basis, &key)
                    .unwrap();
            ScalarConductorResponse::new(&cards[i], conductor, &[("operating", coordinate)])
                .unwrap()
        })
        .collect();
    PlanarConductionResponse::new(mesh(1.0), sources, limits(), &mut || false).unwrap()
}

// Independent two-body balance with G=1 W/K, V_i=1 m3, mass_i=1 kg.
// rho=r+alpha*dT+beta*E gives E=J*(r+alpha*dT)/(1-beta*J).
// At a prescribed J, the remaining TWO enthalpy equations are linear. These
// formulas are fixture references only; runtime uses neither this reduction
// nor a uniform-resistor approximation.
fn endpoint_at_current(
    j: f64,
    dt: f64,
    old: [f64; 2],
    r: [f64; 2],
    beta: [f64; 2],
    field_law: bool,
) -> ([f64; 2], [f64; 2]) {
    let d = beta.map(|b| if field_law { 1.0 - b * j } else { 1.0 });
    let base = std::array::from_fn::<_, 2, _>(|i| r[i] + if field_law { 0.0 } else { beta[i] * j });
    let a = std::array::from_fn::<_, 2, _>(|i| CP[i] + dt - dt * j * j * ALPHA[i] / d[i]);
    let rhs = std::array::from_fn::<_, 2, _>(|i| CP[i] * old[i] + dt * j * j * base[i] / d[i]);
    let det = a[0] * a[1] - dt * dt;
    let theta = [
        (rhs[0] * a[1] + dt * rhs[1]) / det,
        (rhs[1] * a[0] + dt * rhs[0]) / det,
    ];
    let e = std::array::from_fn(|i| j * (base[i] + ALPHA[i] * theta[i]) / d[i]);
    (theta, e)
}
fn voltage_reference(
    voltage: f64,
    dt: f64,
    old: [f64; 2],
    r: [f64; 2],
    beta: [f64; 2],
) -> (f64, [f64; 2], [f64; 2]) {
    let (mut low, mut high) = (0.0, 4.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        let (_, e) = endpoint_at_current(mid, dt, old, r, beta, true);
        if e[0] + e[1] < voltage {
            low = mid;
        } else {
            high = mid;
        }
    }
    let j = 0.5 * (low + high);
    let (theta, e) = endpoint_at_current(j, dt, old, r, beta, true);
    (j, theta, e)
}

#[test]
fn nonlinear_voltage_material_swap_closes_nonuniform_current_finite_heat_and_continuation() {
    let coordinate = ConductorMagnitude::ElectricField;
    let beta = [0.05, 0.08];
    let mut totals = Vec::new();
    for first in [1.0, 2.0] {
        let r = [first, 2.0];
        let (a, ca) = material(r[0], CP[0], ALPHA[0], beta[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], beta[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let thermal = initial(&curves, 999.0, None);
        let original = problem(&cards, coordinate);
        let context = point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        let voltage = ends(6.0, 13.0);
        let next = layout
            .step_response_candidate(
                &thermal,
                0.1,
                &original,
                ConductionDrive::Voltage(&voltage),
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (j, theta, e) = voltage_reference(6.0, 0.1, [0.0; 2], r, beta);
        let step = &next.report().electrical;
        for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
            near(
                next.thermal().temperature_k(cell).unwrap(),
                300.0 + theta[i],
            );
            for triangle in 2 * i..2 * i + 2 {
                near(step.field.electric_field_v_m[triangle][0], -e[i]);
                near(step.field.current_density_a_m2[triangle][0], -j);
                near(step.field.cell_power_w[triangle], 0.5 * e[i] * j);
                let source = &step.field.evaluated[triangle].conductor;
                near(source.material().query().axes()["T"], 300.0 + theta[i]);
                near(source.material().query().axes()["operating"], e[i]);
                assert_eq!(source.source_key().name(), "original_transport");
                for property in source.material().properties() {
                    let before = original.conductors()[triangle]
                        .material()
                        .property(property.requirement().name())
                        .unwrap();
                    assert_eq!(property.requirement(), before.requirement());
                    assert_eq!(
                        property.answer().receipt.selected,
                        before.answer().receipt.selected
                    );
                }
            }
        }
        near(step.heat_j, 0.1 * 6.0 * j);
        near(step.supplied_work_j, step.heat_j);
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            CP[0] * theta[0] + CP[1] * theta[1],
        );
        near(next.report().deposited_joule_j, step.heat_j);
        assert!(next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j);
        assert!(step.field.field_solves > 1);
        assert!(next.report().circuit_trials > step.field.field_solves);
        near(next.report().source.links[0].conductance_w_k, 1.0);
        let continued = layout
            .step_response_candidate(
                next.thermal(),
                0.2,
                next.problem(),
                ConductionDrive::Voltage(&voltage),
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (j2, theta2, _) = voltage_reference(6.0, 0.1, theta, r, beta);
        near(
            continued.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta2[0],
        );
        near(
            continued.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + theta2[1],
        );
        near(
            continued.report().electrical.field.current_density_a_m2[0][0],
            -j2,
        );
        assert!(j2 < j);
        assert_eq!(context.axes()["operating"], 0.0);
        totals.push(step.heat_j);
    }
    assert!(totals[1] < totals[0]);
}

#[test]
fn current_density_source_law_changes_solved_voltage_and_same_material_heating() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let r = [1.0, 2.0];
    let mut heats = Vec::new();
    for b0 in [0.05, 0.2] {
        let beta = [b0, 0.08];
        let (a, ca) = material(r[0], CP[0], ALPHA[0], beta[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], beta[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let thermal = initial(&curves, 999.0, None);
        let original = problem(&cards, coordinate);
        let context = point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let traces = [ConductionCurrentTrace {
            vertices: [2, 5],
            inward_current_density_a_m2: 2.0,
        }];
        let next = network(&bindings)
            .step_response_candidate(
                &thermal,
                0.1,
                &original,
                ConductionDrive::CurrentDensity {
                    voltage: &[(0, 13.0), (3, 13.0)],
                    currents: &traces,
                },
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (theta, e) = endpoint_at_current(2.0, 0.1, [0.0; 2], r, beta, false);
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta[0],
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + theta[1],
        );
        let field = &next.report().electrical.field;
        near(field.voltage_v[2], 13.0 + e[0] + e[1]);
        for i in 0..4 {
            near(field.current_density_a_m2[i][0], -2.0);
            near(field.cell_power_w[i], e[i / 2]);
            near(
                field.evaluated[i].conductor.material().query().axes()["operating"],
                2.0,
            );
        }
        near(
            next.report().electrical.supplied_work_j,
            0.2 * (e[0] + e[1]),
        );
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            next.report().electrical.heat_j,
        );
        heats.push(next.report().electrical.heat_j);
    }
    assert!(heats[1] > heats[0]);
}

#[test]
fn nested_work_budget_late_cancellation_and_source_refusal_keep_both_original_histories() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], 0.05, coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let thermal = initial(&curves, 999.0, None);
    let original = problem(&cards, coordinate);
    let context = point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let saved = thermal.cells().clone();
    let sources = original.conductors().to_vec();
    let voltage = ends(6.0, 0.0);
    let mut polls = 0;
    let first = layout
        .step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Voltage(&voltage),
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Voltage(&voltage),
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut short = config();
    short.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Voltage(&voltage),
            limits(),
            field_controls(),
            controls(),
            short,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Current {
                positive: &[2, 5],
                negative: &[0, 3],
                current_a: 2.0
            },
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::Source(
            ConductionThermalError::Invalid(_)
        ))
    ));
    assert!(
        layout
            .step_response_candidate(
                &thermal,
                0.1,
                &original,
                ConductionDrive::Voltage(&ends(100.0, 0.0)),
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    // Accepted thermal state and the OLD electrical endpoint cannot be paired.
    assert!(
        layout
            .step_response_candidate(
                first.thermal(),
                0.2,
                &original,
                ConductionDrive::Voltage(&voltage),
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.conductors(), sources);
    let retry = layout
        .step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Voltage(&voltage),
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal().cells(), first.thermal().cells());
    assert_eq!(retry.problem().conductors(), first.problem().conductors());
    assert_eq!(retry.report().circuit_trials, first.report().circuit_trials);
}

#[test]
fn thermal_source_applicability_cannot_hide_an_unsolved_frozen_electrical_condition() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], 0.05, coordinate, true);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let thermal = initial(&curves, 1.0, None);
    let original = problem(&cards, coordinate);
    let context = point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    assert!(matches!(
        network(&bindings).step_response_candidate(
            &thermal,
            0.1,
            &original,
            ConductionDrive::Voltage(&ends(6.0, 0.0)),
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::Source(
            ConductionThermalError::Invalid(_)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        original.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}

#[test]
fn sourced_caloric_response_and_applicability_survive_plain_curve_handoff_and_refuse() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        for property in ["specific_heat_capacity", "density"] {
            for response in [false, true] {
                let (a, ca) = material_with_caloric(
                    1.0,
                    CP[0],
                    ALPHA[0],
                    0.05,
                    coordinate,
                    false,
                    Some((property, response)),
                );
                let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
                let (cards, curves) = ([a, b], [ca, cb]);
                let thermal = initial(&curves, 1.0, None);
                let original = problem(&cards, coordinate);
                let context = point(300.0, coordinate);
                let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
                let bindings = cells(&cards, &curves, &context, &key);
                // The supplied chart is an ordinary valid curve at reference
                // E/J=0. Its exact source receipts remain present even through
                // MaterialCell's plain EquilibriumEnthalpyPhaseCurve reference.
                assert!(!bindings[0].curve.source_receipts().is_empty());
                assert_eq!(bindings[0].curve.source_receipts(), curves[0].receipts());
                let before = thermal.cells().clone();
                assert!(matches!(
                    network(&bindings).step_response_candidate(
                        &thermal,
                        0.1,
                        &original,
                        ConductionDrive::Voltage(&ends(6.0, 0.0)),
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                        || false
                    ),
                    Err(CircuitEnthalpyError::Source(
                        ConductionThermalError::Invalid(
                            "caloric source response and applicability must be independent of E/J"
                        )
                    ))
                ));
                assert_eq!(thermal.cells(), &before);
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(
                    original.conductors()[0].material().query().axes()["operating"],
                    0.0
                );
            }
        }
    }
}

#[test]
fn original_melting_receipts_prevent_frozen_field_dependent_fusion_energy() {
    use fs_couple::circuit_enthalpy::{
        conduction::ConductionThermalNetwork, material::MaterialCell,
    };
    use fs_material::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};
    use fs_material::thermal::source::ThermalMaterialPoint;
    use fs_thermal::enthalpy::*;

    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        for condition in [None, Some(false), Some(true)] {
            let (a, _) = material_with_caloric(
                1.0,
                CP[0],
                ALPHA[0],
                0.05,
                coordinate,
                false,
                condition.map(|response| ("latent_heat_fusion", response)),
            );
            let (b, _) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
            let cards = [a, b];
            let context = point(300.0, coordinate)
                .with_quantity(
                    "liquid_fraction",
                    QuantitySpec::dimensional(fs_qty::Dims::NONE),
                    0.0,
                )
                .unwrap();
            let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
            let density = PropertyKey::new("density", Density::DIMS);
            let tm = PropertyKey::with_quantity(
                "melting_point",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
            );
            let latent = PropertyKey::new("latent_heat_fusion", fs_qty::Dims([2, 0, -2, 0, 0, 0]));
            let points = [200.0, 1000.0, 2000.0].map(|t| context.clone().with("T", t).unwrap());
            let curves = std::array::from_fn::<_, 2, _>(|i| {
                ResolvedMeltingCurve::try_new(
                    &cards[i],
                    MeltingSourceQuery {
                        solid: HeatCapacitySourceQuery {
                            heat_capacity: &cp,
                            density: &density,
                            temperature_axis: "T",
                            lower: &points[0],
                            upper: &points[1],
                            selection: SelectionPolicy::SingleClaimOnly,
                        },
                        liquid: HeatCapacitySourceQuery {
                            heat_capacity: &cp,
                            density: &density,
                            temperature_axis: "T",
                            lower: &points[1],
                            upper: &points[2],
                            selection: SelectionPolicy::SingleClaimOnly,
                        },
                        melting_point: &tm,
                        latent_heat_fusion: &latent,
                        transition: &points[1],
                        transition_selection: SelectionPolicy::SingleClaimOnly,
                        phase_selector: None,
                    },
                    0.0,
                    1e-6,
                )
                .unwrap()
            });
            let events = [CellId(30), CellId(10)]
                .into_iter()
                .enumerate()
                .map(|(i, id)| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id,
                        segment_id: i as u64,
                        body_id: i as u64,
                        material_id: MaterialId(i as u64),
                        mass_kg: 1.0,
                        specific_h_j_kg: 100.0 * CP[i],
                        geometry: CellGeometry {
                            start_m: [i as f64, 0.0, 0.0],
                            end_m: [i as f64 + 1.0, 0.0, 0.0],
                            volume_m3: 1.0,
                            surface_area_m2: 6.0,
                        },
                        boundary: super::common::boundary(),
                    }),
                })
                .collect();
            let laws = [
                (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
                (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
            ];
            let thermal = DepositionSimulation::new(
                events,
                &laws,
                DepositionBudget {
                    max_cells: 2,
                    max_links: 1,
                    max_events: 4,
                    max_substeps: 100,
                },
            )
            .unwrap()
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: 0.1,
                    energy_tolerance_j: 1e-7,
                },
                || Ok(()),
            )
            .unwrap()
            .0;
            let bindings: Vec<_> = [CellId(30), CellId(10)]
                .into_iter()
                .enumerate()
                .map(|(i, cell)| MaterialCell {
                    cell,
                    material_id: MaterialId(i as u64),
                    curve: curves[i].curve(),
                    point: ThermalMaterialPoint {
                        material: &cards[i],
                        context: &context,
                        temperature_axis: "T",
                        liquid_fraction_axis: Some("liquid_fraction"),
                    },
                    conductivity: None,
                })
                .collect();
            let layout = ConductionThermalNetwork {
                cells: &bindings,
                links: &[],
                regions: &super::REGIONS,
                volume_tolerance_m3: 1e-12,
            };
            let original = problem_at(&cards, coordinate, &context);
            assert_eq!(bindings[0].curve.source_receipts(), curves[0].receipts());
            assert!(
                bindings[0]
                    .curve
                    .source_receipts()
                    .iter()
                    .any(|receipt| receipt.property == "latent_heat_fusion")
            );
            let result = layout.step_response_candidate(
                &thermal,
                0.1,
                &original,
                ConductionDrive::Voltage(&ends(6.0, 0.0)),
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            );
            if condition.is_some() {
                assert!(matches!(
                    result,
                    Err(CircuitEnthalpyError::Source(
                        ConductionThermalError::Invalid(
                            "caloric source response and applicability must be independent of E/J"
                        )
                    ))
                ));
            } else {
                let accepted = result.unwrap();
                near(
                    accepted.thermal().energy_j() - thermal.energy_j(),
                    accepted.report().electrical.heat_j,
                );
                assert_eq!(accepted.thermal().time_s(), 0.1);
            }
            assert_eq!(thermal.time_s(), 0.0);
        }
    }
}
