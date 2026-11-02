//! Independent Faraday/Ohm + sliding-play balances; manufactured data only.
use super::*;
use crate::material_circuit::hysteretic::thermal::{
    HystereticEnthalpyCandidate, lamination::HystereticLamination,
};
use fs_material::{
    conductor::resolve_conductor_state_with_basis,
    field_storage::hysteresis::lamination::{LaminatedCoreEddyCurrent, ThinLaminationGeometry},
};

const CORE: CellId = CellId(10);
const COIL: CellId = CellId(11);
const DT: f64 = 50.0;
const VOLTAGE: f64 = 0.005;
const WINDING_R: f64 = 0.001;
const BASE_THICKNESS: f64 = 0.01;
const BASE_RHO: f64 = 2.0 * BASE_THICKNESS * BASE_THICKNESS / 12.0;
const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity),
    (3, ElectricalPropertyBasis::Resistivity),
];

fn conducting_core(rho: f64, alpha: f64, upper: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = core_material(2.0, upper);
    let mut claims = ClaimSet::new();
    for (_, original_claim) in original.claims().claims_ordered() {
        let mut claim = original_claim.clone();
        if claim.key.name() == "electrical_resistivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (290.0, rho * (1.0 - 10.0 * alpha)),
                    (310.0, rho * (1.0 + 10.0 * alpha)),
                ],
                dims: claim.key.dims(),
            };
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                290.0,
                310.0,
            );
            claim.provenance.source = "manufactured insulated-lamination rho(T)".into();
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn with_lamination<R>(
    rho: f64,
    alpha: f64,
    thickness: f64,
    upper: f64,
    run: impl FnOnce(
        &DepositionSimulation<'_>,
        &ThermalHystereticCircuit<'_, '_>,
        &MaterialNetwork<'_>,
    ) -> R,
) -> R {
    let (a, ca) = conducting_core(rho, alpha, upper);
    let (b, cb) = material(WINDING_R, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let thermal = thermal_initial(&curves);
    let core = sourced_core(&cards[0]);
    let conductor = resolve_conductor_state_with_basis(
        &cards[0],
        &context,
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    let lamination = HystereticLamination {
        branch: 3,
        loss: LaminatedCoreEddyCurrent::new(
            &core,
            conductor,
            ThinLaminationGeometry {
                thickness_m: thickness,
                width_m: 1.0,
                maximum_thickness_ratio: 0.05,
            },
        )
        .unwrap(),
    };
    let [coil] = winding();
    let resistors = [coil, lamination.resistor_binding(CORE).unwrap()];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let sources = network
        .prepare_with_bases(&thermal, &BASES, || false)
        .unwrap()
        .resolve(&thermal, || false)
        .unwrap()
        .resistors;
    let mut graph = rl_graph();
    graph
        .branches
        .push((2, 0, Branch::Resistor { ohms: 999.0 }));
    // Low-frequency monotone segment, always below saturation H=7.
    let start = initial(&core, 3.0, vec![VOLTAGE], vec![]);
    let circuit = HystereticCircuit::new(graph, 2, core, sources, start).unwrap();
    let driver = ThermalHystereticCircuit::bind_with_lamination(
        circuit,
        &thermal,
        &network,
        &BASES,
        CORE,
        [290.0, upper],
        lamination,
        || false,
    )
    .unwrap();
    run(&thermal, &driver, &network)
}

struct Reference {
    core_t: f64,
    coil_t: f64,
    field: f64,
    flux_increment: f64,
    pinning: f64,
    eddy: f64,
    coil_heat: f64,
    coil_current: f64,
    resistance: f64,
}

// Eliminate the native circuit analytically on the sliding segment L=mu+chi.
// The remaining scalar is the finite core's caloric balance; no production
// material evaluator, circuit solver or thermal residual appears in this oracle.
fn reference(rho: f64, alpha: f64, thickness: f64, old_t: [f64; 2], old_h: f64) -> Reference {
    let evaluate = |core_t: f64| {
        let resistance = 12.0 * rho / thickness.powi(2) * (1.0 + alpha * (core_t - 300.0));
        let flux_increment = (VOLTAGE - WINDING_R * old_h)
            / (WINDING_R / (2.0 * 0.21) + (1.0 + WINDING_R / resistance) / DT);
        let field = old_h + flux_increment / 0.21;
        let pinning = 0.4 * (field - old_h);
        let eddy = flux_increment.powi(2) / (DT * resistance);
        let coil_current = (VOLTAGE - flux_increment / DT) / WINDING_R;
        let coil_heat = DT * WINDING_R * coil_current.powi(2);
        Reference {
            core_t,
            coil_t: old_t[1] + coil_heat,
            field,
            flux_increment,
            pinning,
            eddy,
            coil_heat,
            coil_current,
            resistance,
        }
    };
    let residual = |t: f64| {
        let r = evaluate(t);
        r.core_t - old_t[0] - r.pinning - r.eddy
    };
    let (mut lo, mut hi) = (old_t[0], 310.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let r = evaluate(0.5 * (lo + hi));
    assert!(residual(r.core_t).abs() < 1e-11);
    assert!(r.field > old_h && r.field < 7.0);
    r
}

fn check_step(
    before_t: &DepositionSimulation<'_>,
    before_c: &ThermalHystereticCircuit<'_, '_>,
    next: &HystereticEnthalpyCandidate<'_, '_, '_>,
    r: &Reference,
) {
    let native = next.circuit().circuit();
    let step = &next.report().electrical;
    near(native.current_a(), r.field, 3e-8);
    near(
        native.flux_linkage_wb() - before_c.circuit().flux_linkage_wb(),
        r.flux_increment,
        3e-9,
    );
    near(step.terminal_voltage_v, r.flux_increment / DT, 3e-10);
    near(step.core.core_heat_j, r.pinning, 3e-9);
    near(next.thermal().temperature_k(CORE).unwrap(), r.core_t, 3e-8);
    near(next.thermal().temperature_k(COIL).unwrap(), r.coil_t, 3e-8);
    let coil = step
        .circuit
        .resistors
        .iter()
        .find(|r| r.material.branch == 1)
        .unwrap();
    let eddy = step
        .circuit
        .resistors
        .iter()
        .find(|r| r.material.branch == 3)
        .unwrap();
    near(coil.current_a, r.coil_current, 3e-8);
    near(coil.heat_j, r.coil_heat, 3e-9);
    near(eddy.voltage_v, r.flux_increment / DT, 3e-10);
    near(
        eddy.current_a,
        r.flux_increment / (DT * r.resistance),
        3e-10,
    );
    near(eddy.heat_j, r.eddy, 3e-10);
    near(eddy.material.conductor.resistance_ohm(), r.resistance, 3e-8);
    // Both distinct loss receipts enter the SAME physical core exactly once.
    near(
        next.thermal().temperature_k(CORE).unwrap() - before_t.temperature_k(CORE).unwrap(),
        step.core.core_heat_j + eddy.heat_j,
        3e-8,
    );
    near(
        next.report().deposited_joule_j,
        r.pinning + r.eddy + r.coil_heat,
        3e-8,
    );
    near(
        next.thermal().energy_j() - before_t.energy_j() + step.stored_change_j,
        step.supplied_work_j,
        3e-8,
    );
    near(next.report().total_energy_residual_j, 0.0, 2e-8);
    let accepted = next
        .circuit()
        .lamination()
        .unwrap()
        .loss
        .equivalent_conductor();
    assert_eq!(accepted, &eddy.material.conductor);
    assert_eq!(
        accepted.material().query(),
        point(next.thermal().temperature_k(CORE).unwrap())
    );
    assert_eq!(
        native.core().source().unwrap().material().query(),
        accepted.material().query()
    );
    assert_eq!(
        accepted.source_property().answer().receipt.selected,
        before_c
            .lamination()
            .unwrap()
            .loss
            .equivalent_conductor()
            .source_property()
            .answer()
            .receipt
            .selected
    );
}

#[test]
fn source_resistivity_and_sheet_thickness_react_on_circuit_and_heat_the_original_core_once() {
    let variants = [
        (BASE_RHO, 0.05, BASE_THICKNESS),
        (BASE_RHO, 0.0, BASE_THICKNESS),
        (BASE_RHO, 0.05, 2.0 * BASE_THICKNESS),
        (2.0 * BASE_RHO, 0.05, BASE_THICKNESS),
    ];
    let results: Vec<_> = variants
        .into_iter()
        .map(|(rho, alpha, d)| {
            with_lamination(rho, alpha, d, 310.0, |thermal, driver, _| {
                let saved_cells = thermal.cells().clone();
                let saved_circuit = driver.circuit().clone();
                let saved_loss = driver.lamination().unwrap().clone();
                let expected = reference(rho, alpha, d, [300.0; 2], 3.0);
                let next = driver
                    .step_candidate(
                        thermal,
                        DT,
                        &[VOLTAGE],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                check_step(thermal, driver, &next, &expected);
                assert_eq!(thermal.cells(), &saved_cells);
                assert_eq!(driver.circuit(), &saved_circuit);
                assert_eq!(driver.lamination().unwrap(), &saved_loss);
                let result = (next.circuit().circuit().current_a(), expected.eddy);
                let (t, c, _) = next.into_parts();
                // A second accepted interval retains warm source pins and the
                // original non-virgin polarization; no lamination history is reset.
                let expected = reference(
                    rho,
                    alpha,
                    d,
                    [
                        t.temperature_k(CORE).unwrap(),
                        t.temperature_k(COIL).unwrap(),
                    ],
                    c.circuit().current_a(),
                );
                let saved_memory = c.circuit().magnetic_state().clone();
                let continued = c
                    .step_candidate(
                        &t,
                        2.0 * DT,
                        &[VOLTAGE],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                check_step(&t, &c, &continued, &expected);
                assert_eq!(c.circuit().magnetic_state(), &saved_memory);
                result
            })
        })
        .collect();
    assert!(results[0].0 - results[1].0 > 1e-7); // Actual rho(T) feedback changes field.
    assert!(results[2].0 < results[0].0 && results[2].1 > results[0].1);
    assert!(results[3].0 > results[0].0 && results[3].1 < results[0].1);
}

#[test]
fn only_the_physical_core_port_is_admitted_and_late_thermal_refusal_is_transactional() {
    with_lamination(
        BASE_RHO,
        0.05,
        BASE_THICKNESS,
        300.1,
        |thermal, driver, network| {
            let saved_circuit = driver.circuit().clone();
            let saved_cells = thermal.cells().clone();
            let saved_loss = driver.lamination().unwrap().clone();
            // Legacy entry points retain their separate core/winding inventory rule.
            assert!(
                ThermalHystereticCircuit::bind(
                    saved_circuit.clone(),
                    thermal,
                    network,
                    &BASES,
                    CORE,
                    [290.0, 300.1],
                    || false,
                )
                .is_err()
            );
            let core = saved_circuit.core().clone();
            let mut wrong_graph = rl_graph();
            wrong_graph
                .branches
                .push((1, 0, Branch::Resistor { ohms: 999.0 }));
            let wrong = HystereticCircuit::new(
                wrong_graph,
                2,
                core.clone(),
                saved_circuit.native_circuit().resistors().to_vec(),
                initial(&core, 3.0, vec![VOLTAGE], vec![]),
            )
            .unwrap();
            assert!(
                ThermalHystereticCircuit::bind_with_lamination(
                    wrong,
                    thermal,
                    network,
                    &BASES,
                    CORE,
                    [290.0, 300.1],
                    saved_loss.clone(),
                    || false,
                )
                .is_err()
            );
            // Independent endpoint needs Tcore=300.1702165 K, outside this fit.
            assert!(
                driver
                    .step_candidate(
                        thermal,
                        DT,
                        &[VOLTAGE],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                        || false,
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(driver.circuit(), &saved_circuit);
            assert_eq!(driver.lamination().unwrap(), &saved_loss);
        },
    );
}
