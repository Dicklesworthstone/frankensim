//! Source-only changes to thermal contacts and emissivity reach the real field.
use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    ResistorCell,
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::{CellResistor, MaterialNetwork},
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::*;
use fs_material::{
    phase::{HeatCapacitySourceQuery, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::interface::ThermalInterfacePoint,
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::DepositionSimulation;

fn emissive_material(
    rho: f64,
    cp: f64,
    alpha: f64,
    e: f64,
    upper: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = material(rho, cp, alpha);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, 0.99 * e),
                    (upper, e * (1.0 + 0.0001 * (upper - 300.0))),
                ],
                dims: EMISSIVITY_DIMS,
            },
            validity: fs_evidence::ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    upper,
                )
                .with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic electrothermal surface law".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    // The caloric source is bound to the complete new card, not the old card.
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
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

fn interface(cards: &[MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
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
            validity: fs_evidence::ValidityDomain::unconstrained()
                .with("aT", 200.0, 2000.0)
                .with("bT", 200.0, 2000.0)
                .with("pressure", 100000.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic thermal interface law".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "a".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "b".into(),
        },
        SystemContext {
            medium: "fixed contact".into(),
            third_body: None,
            environment: "authored".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_heat<R>(
    e: f64,
    resistance: f64,
    emissivity_upper: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        PlanarConduction,
    ) -> R,
) -> R {
    let (a, ca) = emissive_material(1.0, 1.0, 0.002, e, emissivity_upper);
    let (b, cb) = emissive_material(3.0, 2.0, 0.0, 0.6, 2000.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
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
    let interface = interface(&cards, resistance);
    let interface_context = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 300.0)
        .unwrap()
        .with("pressure", 100000.0)
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
    let emissivities = [0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original surface",
            bindings[i].point,
            curves[i]
                .curve()
                .state_at_specific_enthalpy(100.0 * (i + 1) as f64)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    });
    // Original infinite-parallel-plate limiting geometry, one square metre.
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
    let thermal = initial(&curves, 999.0, None);
    let field = problem(
        &cards,
        [ElectricalPropertyBasis::Resistivity; 2],
        300.0,
        1.0,
    );
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

#[test]
fn current_field_heat_contact_and_radiation_satisfy_independent_two_body_balance() {
    for (e, resistance) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
        with_heat(e, resistance, 2000.0, |cx, _, assembly, thermal, field| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let next = prepared
                .step_conduction_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Current {
                        positive: &[2, 5],
                        negative: &[0, 3],
                        current_a: 2.0,
                    },
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            // Independent scalar root: total energy eliminates the second
            // temperature. C1=1, C2=2; I=2; rho1=1+.002*d1, rho2=3.
            // Q12=G(T1-T2)+sigma(T1^4-T2^4)/(1/e1+1/e2-1).
            let state = |d1: f64| {
                let d2 = (8.0 - 0.996 * d1) / 2.0;
                let t1 = 300.0 + d1;
                let t2 = 300.0 + d2;
                let conductance = 1.0 / (1.0 + resistance * (1.0 + 0.001 * d1));
                let e1 = e * (1.0 + 0.0001 * d1);
                let e2 = 0.6 * (1.0 + 0.0001 * d2);
                let radiation =
                    5.670_374_419e-8 * (t1.powi(4) - t2.powi(4)) / (1.0 / e1 + 1.0 / e2 - 1.0);
                let transfer = conductance * (d1 - d2) + radiation;
                (d2, d1 - 0.5 * (4.0 * (1.0 + 0.002 * d1) - transfer))
            };
            let (mut lo, mut hi) = (0.0, 8.0);
            assert!(state(lo).1 < 0.0 && state(hi).1 > 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if state(mid).1 < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let d1 = 0.5 * (lo + hi);
            let d2 = state(d1).0;
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + d1,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + d2,
            );
            let report = next.report();
            near(report.electrical.voltage_v[2], 2.0 * (4.0 + 0.002 * d1));
            near(report.electrical.supplied_work_j, d1 + 2.0 * d2);
            near(report.deposited_joule_j, report.electrical.supplied_work_j);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                report.electrical.supplied_work_j,
            );
            let contact = &report.source.thermal_contacts[0];
            near(
                contact.interface.value_m2_k_per_w(),
                resistance * (1.0 + 0.001 * d1),
            );
            for surface in report.source.radiation.as_ref().unwrap().surfaces() {
                let cell = report
                    .source
                    .material
                    .cells
                    .iter()
                    .find(|c| c.cell == surface.cell)
                    .unwrap();
                assert_eq!(cell.phase, surface.phase);
            }
            assert!(
                (report.total_energy_residual_j + 0.5 * report.source.radiation_net_power_w())
                    .abs()
                    < config().total_energy_tolerance_j
            );
            assert_eq!(thermal.time_s(), 0.0);
            // Retain the prepared sources and actual accepted spatial material
            // snapshots across a second physical step.
            let continuing = prepared
                .step_conduction_candidate(
                    cx,
                    next.thermal(),
                    1.0,
                    next.problem(),
                    ConductionDrive::Current {
                        positive: &[2, 5],
                        negative: &[0, 3],
                        current_a: 2.0,
                    },
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            assert!(continuing.thermal().temperature_k(CellId(30)).unwrap() > 300.0 + d1);
        });
    }
}

#[test]
fn surface_and_interface_source_swaps_change_real_voltage_driven_current() {
    let run = |e, resistance| {
        with_heat(e, resistance, 2000.0, |cx, _, assembly, thermal, field| {
            let selected = assembly.prepare(cx, &thermal).unwrap();
            let next = selected
                .step_conduction_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            let t1 = next.thermal().temperature_k(CellId(30)).unwrap();
            let t2 = next.thermal().temperature_k(CellId(10)).unwrap();
            let current = -next.report().electrical.current_density_a_m2[0][0];
            near(current, 8.0 / (4.0 + 0.002 * (t1 - 300.0)));
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                0.5 * 8.0 * current,
            );
            [t1, t2, current]
        })
    };
    let original = run(0.2, 0.2);
    for substituted in [run(0.8, 0.2), run(0.2, 2.0)] {
        for i in 0..3 {
            assert!((original[i] - substituted[i]).abs() > 1e-8);
        }
    }
}

#[test]
fn empty_heat_paths_share_the_original_field_step_and_cancellation_keeps_history() {
    with_heat(0.2, 0.2, 2000.0, |cx, gate, assembly, thermal, field| {
        let bare = MaterialHeatNetwork {
            electrical: assembly.electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let selected = bare.prepare(cx, &thermal).unwrap();
        let next = selected
            .step_conduction_candidate(
                cx,
                &thermal,
                0.5,
                &field,
                ConductionDrive::Voltage(&ends(8.0, 0.0)),
                &REGIONS,
                1e-12,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
            )
            .unwrap();
        let original = network(assembly.electrical.bulk.cells)
            .step_candidate(
                &thermal,
                0.5,
                &field,
                &ends(8.0, 0.0),
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_eq!(next.thermal().cells(), original.thermal().cells());
        assert_eq!(
            next.report().electrical.voltage_v,
            original.report().electrical.voltage_v
        );
        assert_eq!(
            next.report().electrical.cell_power_w,
            original.report().electrical.cell_power_w
        );
        gate.request();
        assert!(matches!(
            selected.step_conduction_candidate(
                cx,
                &thermal,
                0.5,
                &field,
                ConductionDrive::Voltage(&ends(8.0, 0.0)),
                &REGIONS,
                1e-12,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
    });
}

#[test]
fn endpoint_emissivity_support_and_unstepped_circuit_bindings_refuse_atomically() {
    with_heat(0.2, 0.2, 300.1, |cx, _, assembly, thermal, field| {
        let selected = assembly.prepare(cx, &thermal).unwrap();
        assert!(
            selected
                .step_conduction_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    ConductionDrive::Voltage(&ends(8.0, 0.0)),
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2]
                )
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
        // A prepared graph binding must never silently vanish just because a
        // spatial field replaces the external electrical driver.
        let resistor = [CellResistor {
            destination: ResistorCell {
                branch: 0,
                cell: CellId(30),
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        }];
        let bulk = MaterialNetwork {
            cells: assembly.electrical.bulk.cells,
            resistors: &resistor,
            links: &LINKS,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[(0, ElectricalPropertyBasis::Resistivity)],
        };
        let unstepped = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let selected = unstepped.prepare(cx, &thermal).unwrap();
        assert!(matches!(
            selected.step_conduction_candidate(
                cx,
                &thermal,
                0.5,
                &field,
                ConductionDrive::Voltage(&ends(8.0, 0.0)),
                &REGIONS,
                1e-12,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Source(
                fs_couple::circuit_enthalpy::conduction::ConductionThermalError::Heat(
                    HeatNetworkError::Invalid(_)
                )
            ))
        ));
    });
}
