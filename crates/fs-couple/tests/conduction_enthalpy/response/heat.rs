//! G1/G3/G4: nonlinear spatial material response and complete surface heat.
//! All source numbers are manufactured; no empirical conductor claim is made.
use super::*;
use crate::{
    LINKS, REGIONS,
    heat::{interface, radiative},
};
use fs_conduction::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_thermal::enthalpy::DepositionSimulation;

#[path = "heat/circuit.rs"]
mod circuit;

const RHO: [f64; 2] = [1.0, 3.0];
const BETA: [f64; 2] = [0.03, 0.04];

#[derive(Clone, Copy, PartialEq, Eq)]
enum Support {
    Ordinary,
    ShortEmissivity,
    EmissivityField,
    ContactField,
}

fn emitting_material(
    i: usize,
    beta: f64,
    emissivity: f64,
    support: Support,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let field = ConductorMagnitude::ElectricField;
    let (base, _) = material(RHO[i], CP[i], ALPHA[i], beta, field, false);
    let mut claims = base.claims().clone();
    let upper = if support == Support::ShortEmissivity && i == 0 {
        300.1
    } else {
        2000.0
    };
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        upper,
    );
    if support == Support::EmissivityField && i == 0 {
        validity = validity.with_quantity("operating", field.quantity(), 0.0, 20.0);
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, 0.99 * emissivity),
                    (upper, emissivity * (1.0 + 0.0001 * (upper - 300.0))),
                ],
                dims: EMISSIVITY_DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured emissive nonlinear conductor".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0, field),
            upper: &point(2000.0, field),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn with_heat<R>(
    e: f64,
    resistance: f64,
    beta: [f64; 2],
    support: Support,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        PlanarConductionResponse,
    ) -> R,
) -> R {
    let (a, ca) = emitting_material(0, beta[0], e, support);
    let (b, cb) = emitting_material(1, beta[1], 0.6, support);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0, ConductorMagnitude::ElectricField);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let bulk = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let original_interface = interface(&cards, resistance);
    let interface = if support == Support::ContactField {
        let mut claims = ClaimSet::new();
        for (_, original) in original_interface.claims().claims_ordered() {
            let mut claim = original.clone();
            claim.validity = claim.validity.with_quantity(
                "operating",
                ConductorMagnitude::ElectricField.quantity(),
                0.0,
                20.0,
            );
            claims.insert_claim(claim).unwrap();
        }
        InterfaceSystemCard::assemble(
            original_interface.surface_a().clone(),
            original_interface.surface_b().clone(),
            original_interface.context().clone(),
            claims,
            vec![],
        )
        .unwrap()
    } else {
        original_interface
    };
    let interface_context = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 300.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
        .with_quantity(
            "operating",
            ConductorMagnitude::ElectricField.quantity(),
            0.0,
        )
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 7,
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &interface_context,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        sides: [0, 1].map(|i| ContactCell {
            cell: bindings[i].cell,
            material_id: bindings[i].material_id,
            curve: curves[i].curve(),
        }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let surfaces = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: bindings[i].cell,
        material_id: bindings[i].material_id,
        curve: curves[i].curve(),
        material: bindings[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let thermal = initial(&curves, 999.0, None);
    let emissivities = [0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "nonlinear conductor surface",
            bindings[i].point,
            curves[i]
                .curve()
                .state_at_specific_enthalpy(100.0 * CP[i])
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    });
    let enclosure = radiative::enclosure(emissivities);
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surfaces,
        power_tolerance_w: 1e-8,
        radiosity_tolerance_w_m2: 1e-8,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let field = problem(&cards, ConductorMagnitude::ElectricField);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal, field)
    })
}

// Independent scalar reduction of two finite unit-volume bodies. At fixed J,
// E_i=J*(rho_i+alpha_i*dT_i)/(1-beta_i*|J|). Total energy eliminates T2;
// the remaining root includes original Fourier/contact and gray radiation laws.
fn reference(j: f64, e: f64, resistance: f64, beta: [f64; 2]) -> ([f64; 2], f64) {
    let dt = 0.5;
    let denom = beta.map(|b| 1.0 - b * j.abs());
    let a = [0, 1].map(|i| CP[i] - dt * j * j * ALPHA[i] / denom[i]);
    let supply = dt * j * j * (RHO[0] / denom[0] + RHO[1] / denom[1]);
    let equation = |d1: f64| {
        let d2 = (supply - a[0] * d1) / a[1];
        let t1 = 300.0 + d1;
        let t2 = 300.0 + d2;
        let g = 1.0 / (1.0 + resistance * (1.0 + 0.001 * d1));
        let e1 = e * (1.0 + 0.0001 * d1);
        let e2 = 0.6 * (1.0 + 0.0001 * d2);
        let transfer = g * (d1 - d2)
            + 5.670_374_419e-8 * (t1.powi(4) - t2.powi(4)) / (1.0 / e1 + 1.0 / e2 - 1.0);
        let power = j * j * (RHO[0] + ALPHA[0] * d1) / denom[0];
        (d2, CP[0] * d1 - dt * (power - transfer))
    };
    let (mut lo, mut hi) = (0.0, 20.0);
    assert!(equation(lo).1 <= 0.0 && equation(hi).1 >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if equation(mid).1 < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let d = [0.5 * (lo + hi), equation(0.5 * (lo + hi)).0];
    let voltage: f64 = (0..2)
        .map(|i| j * (RHO[i] + ALPHA[i] * d[i]) / denom[i])
        .sum();
    (d.map(|d| 300.0 + d), voltage)
}

#[test]
fn sourced_surface_transfer_and_nonlinear_field_close_the_same_energy_balance() {
    for j in [-2.0, 2.0] {
        with_heat(
            0.2,
            0.2,
            BETA,
            Support::Ordinary,
            |cx, _, assembly, thermal, field| {
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let next = selected
                    .step_conduction_response_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        ConductionDrive::Current {
                            positive: &[2, 5],
                            negative: &[0, 3],
                            current_a: j,
                        },
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let (temperature, voltage) = reference(j, 0.2, 0.2, BETA);
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(next.thermal().temperature_k(cell).unwrap(), temperature[i]);
                }
                let report = next.report();
                near(report.electrical.field.voltage_v[2], voltage);
                near(report.electrical.supplied_work_j, 0.5 * j * voltage);
                near(report.deposited_joule_j, report.electrical.supplied_work_j);
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    0.5 * j * voltage,
                );
                assert!(report.circuit_trials > 1);
                for surface in report.source.radiation.as_ref().unwrap().surfaces() {
                    let material = report
                        .source
                        .material
                        .cells
                        .iter()
                        .find(|c| c.cell == surface.cell)
                        .unwrap();
                    assert_eq!(surface.phase, material.phase);
                }
                assert!(
                    (report.total_energy_residual_j + 0.5 * report.source.radiation_net_power_w())
                        .abs()
                        < config().total_energy_tolerance_j
                );
                let continued = selected
                    .step_conduction_response_candidate(
                        cx,
                        next.thermal(),
                        1.0,
                        next.problem(),
                        ConductionDrive::Current {
                            positive: &[2, 5],
                            negative: &[0, 3],
                            current_a: j,
                        },
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                assert!(continued.thermal().temperature_k(CellId(30)).unwrap() > temperature[0]);
                assert_eq!(thermal.time_s(), 0.0);
            },
        );
    }
}

#[test]
fn data_only_contact_emissivity_and_field_response_changes_reach_actual_current() {
    let mut results = Vec::new();
    for (e, resistance, beta) in [
        (0.2, 0.2, BETA),
        (0.8, 0.2, BETA),
        (0.2, 2.0, BETA),
        (0.2, 0.2, [0.08, BETA[1]]),
    ] {
        let result = with_heat(
            e,
            resistance,
            beta,
            Support::Ordinary,
            |cx, _, assembly, thermal, field| {
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let next = selected
                    .step_conduction_response_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        ConductionDrive::Voltage(&ends(8.0, 0.0)),
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let (mut lo, mut hi) = (0.0, 3.0);
                for _ in 0..60 {
                    let j = 0.5 * (lo + hi);
                    if reference(j, e, resistance, beta).1 < 8.0 {
                        lo = j;
                    } else {
                        hi = j;
                    }
                }
                let j = 0.5 * (lo + hi);
                let expected = reference(j, e, resistance, beta).0;
                let found = [
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    -next.report().electrical.field.current_density_a_m2[0][0],
                ];
                for (a, b) in found.into_iter().zip([expected[0], expected[1], j]) {
                    near(a, b);
                }
                near(next.thermal().energy_j() - thermal.energy_j(), 4.0 * j);
                found
            },
        );
        results.push(result);
    }
    for changed in &results[1..] {
        for i in 0..3 {
            assert!((changed[i] - results[0][i]).abs() > 1e-8);
        }
    }
}

#[test]
fn absent_surface_paths_retain_original_solution_and_refusals_keep_all_histories() {
    with_heat(
        0.2,
        0.2,
        BETA,
        Support::Ordinary,
        |cx, gate, assembly, thermal, field| {
            let bare = MaterialHeatNetwork {
                electrical: assembly.electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let selected = bare.prepare(cx, &thermal).unwrap();
            let full = selected
                .step_conduction_response_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            let original = network(assembly.electrical.bulk.cells)
                .step_response_candidate(
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            assert_eq!(full.thermal().cells(), original.thermal().cells());
            assert_eq!(
                full.report().electrical.field.voltage_v,
                original.report().electrical.field.voltage_v
            );
            assert_eq!(
                full.report().electrical.field.cell_power_w,
                original.report().electrical.field.cell_power_w
            );
            let mut exhausted = config();
            exhausted.maximum_circuit_trials = 1;
            assert!(matches!(
                selected.step_conduction_response_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    exhausted,
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            gate.request();
            assert!(matches!(
                selected.step_conduction_response_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.time_s(), 0.0);
            near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
            assert!(
                field
                    .conductors()
                    .iter()
                    .all(|c| c.material().query().axes()["operating"] == 0.0)
            );
        },
    );
    for support in [
        Support::ShortEmissivity,
        Support::EmissivityField,
        Support::ContactField,
    ] {
        with_heat(
            0.2,
            0.2,
            BETA,
            support,
            |cx, _, assembly, thermal, field| {
                let selected = assembly.prepare(cx, &thermal).unwrap();
                assert!(
                    selected
                        .step_conduction_response_candidate(
                            cx,
                            &thermal,
                            0.5,
                            &field,
                            ConductionDrive::Voltage(&ends(8.0, 0.0)),
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 2]
                        )
                        .is_err()
                );
                assert_eq!(thermal.time_s(), 0.0);
                near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
                assert!(
                    field
                        .conductors()
                        .iter()
                        .all(|c| c.material().query().axes()["operating"] == 0.0)
                );
            },
        );
    }
}
