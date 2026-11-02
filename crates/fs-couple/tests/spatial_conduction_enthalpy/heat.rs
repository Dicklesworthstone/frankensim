//! G1/G3/G4: sourced contact pressure and surface radiation feed back into 3D current.
use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{
    Pressure,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn pressure_quantity() -> QuantitySpec {
    // The external source coordinate is held pressure, not a stress response.
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn emissive_material(
    rho: f64,
    cp: f64,
    alpha: f64,
    emissivity: f64,
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
                    (200.0, 0.99 * emissivity),
                    (upper, emissivity * (1.0 + 0.0001 * (upper - 300.0))),
                ],
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                200.0,
                upper,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic 3D electrothermal surface".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
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

fn interface(cards: &[MaterialCard; 2]) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    // Synthetic pressure law, solely to test source propagation. This is not a
    // measured contact model; no pressure dependence is inferred by the solver.
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Curve {
                abscissa: "pressure".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![(100000.0, 2.0), (200000.0, 0.2)],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "pressure",
                pressure_quantity(),
                100000.0,
                200000.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic pressure-dependent thermal interface".into(),
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
    emissivity: f64,
    pressure: f64,
    upper: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        SpatialConduction,
        &[SpatialConductionRegion],
    ) -> R,
) -> R {
    let (a, ca) = emissive_material(1.0, 1.0, 0.002, emissivity, upper);
    let (b, cb) = emissive_material(3.0, 2.0, 0.0, 0.6, 2000.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let interface = interface(&cards);
    let interface_context = QueryPoint::new()
        .with_quantity("aT", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("bT", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("pressure", pressure_quantity(), pressure)
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
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: curves[i].curve(),
        }),
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
    let enclosure = radiative::enclosure([0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original surface",
            cells[i].point,
            curves[i]
                .curve()
                .state_at_specific_enthalpy(100.0 * (i + 1) as f64)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    }));
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surfaces,
        power_tolerance_w: 1e-8,
        radiosity_tolerance_w_m2: 1e-8,
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let thermal = initial(&curves, 999.0, None);
    let conductors = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            &context,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    let field = SpatialConduction::new(
        positions,
        tets,
        (0..12).map(|i| conductors[i / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let regions: Vec<_> = (0..12)
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: cells[tetrahedron / 6].cell,
        })
        .collect();
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
        run(&cx, &gate, &heat, thermal, field, &regions)
    })
}

#[test]
fn full_heat_assembly_closes_independent_three_dimensional_current_and_body_balances() {
    for (emissivity, pressure) in [(0.2, 100000.0), (0.8, 100000.0), (0.2, 200000.0)] {
        with_heat(
            emissivity,
            pressure,
            2000.0,
            |cx, _, heat, thermal, field, regions| {
                let prepared = heat.prepare(cx, &thermal).unwrap();
                let [positive, negative] = electrodes(&field, 2.0);
                let drive = SpatialConductionDrive::total_current(&positive, &negative, 2.0);
                let next = prepared
                    .step_spatial_conduction_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        drive,
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let trace_next = prepared
                    .step_spatial_conduction_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        &boundary(&field, 2.0, 2.0, None),
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                for cell in [CellId(30), CellId(10)] {
                    near(
                        trace_next.thermal().temperature_k(cell).unwrap(),
                        next.thermal().temperature_k(cell).unwrap(),
                    );
                }
                near(
                    trace_next.report().electrical.supplied_work_j,
                    next.report().electrical.supplied_work_j,
                );
                // Independent two-body energy elimination and scalar bisection,
                // unrelated to the field assembly and coupled Newton implementation.
                let resistance = 2.0 - 1.8 * (pressure - 100000.0) / 100000.0;
                let root = |d1: f64| {
                    let d2 = (8.0 - 0.996 * d1) / 2.0;
                    let (t1, t2) = (300.0 + d1, 300.0 + d2);
                    let (e1, e2) = (emissivity * (1.0 + 0.0001 * d1), 0.6 * (1.0 + 0.0001 * d2));
                    let transfer = (d1 - d2) / (1.0 + resistance)
                        + 5.670_374_419e-8 * (t1.powi(4) - t2.powi(4))
                            / (1.0 / e1 + 1.0 / e2 - 1.0);
                    (d2, d1 - 0.5 * (4.0 * (1.0 + 0.002 * d1) - transfer))
                };
                let (mut lo, mut hi) = (0.0, 8.0);
                assert!(root(lo).1 < 0.0 && root(hi).1 > 0.0);
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    if root(mid).1 < 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let d1 = 0.5 * (lo + hi);
                let d2 = root(d1).0;
                near(
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    300.0 + d1,
                );
                near(
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    300.0 + d2,
                );
                let report = next.report();
                near(
                    report.electrical.field.voltage_v[8],
                    2.0 * (4.0 + 0.002 * d1),
                );
                near(
                    report.source.thermal_contacts[0]
                        .interface
                        .value_m2_k_per_w(),
                    resistance,
                );
                near(report.deposited_joule_j, d1 + 2.0 * d2);
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    report.electrical.supplied_work_j,
                );
                assert!(
                    (report.total_energy_residual_j + 0.5 * report.source.radiation_net_power_w())
                        .abs()
                        < config().total_energy_tolerance_j
                );
                assert_eq!(thermal.time_s(), 0.0);
                // Continuing from the jointly accepted field/heat endpoint retains
                // original sources; the same prepared assembly is not frozen at T0.
                let next2 = prepared
                    .step_spatial_conduction_candidate(
                        cx,
                        next.thermal(),
                        1.0,
                        next.problem(),
                        drive,
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                assert!(next2.thermal().temperature_k(CellId(30)).unwrap() > 300.0 + d1);
            },
        );
    }
}

#[test]
fn pressure_and_surface_data_changes_reach_voltage_driven_current() {
    let run = |emissivity, pressure| {
        with_heat(
            emissivity,
            pressure,
            2000.0,
            |cx, _, heat, thermal, field, regions| {
                let next = heat
                    .prepare(cx, &thermal)
                    .unwrap()
                    .step_spatial_conduction_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        &boundary(&field, 2.0, 0.0, Some(8.0)),
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let t1 = next.thermal().temperature_k(CellId(30)).unwrap();
                let current = -next.report().electrical.field.current_density_a_m2[0][0];
                near(current, 8.0 / (4.0 + 0.002 * (t1 - 300.0)));
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    0.5 * 8.0 * current,
                );
                [
                    t1,
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    current,
                ]
            },
        )
    };
    let original = run(0.2, 100000.0);
    for substituted in [run(0.8, 100000.0), run(0.2, 200000.0)] {
        for i in 0..3 {
            assert!((original[i] - substituted[i]).abs() > 1e-8);
        }
    }
}

#[test]
fn refused_surface_trials_and_cancellation_keep_field_and_enthalpy_unchanged() {
    with_heat(
        0.2,
        100000.0,
        300.2,
        |cx, gate, heat, thermal, field, regions| {
            let prepared = heat.prepare(cx, &thermal).unwrap();
            let energy = thermal.energy_j();
            let sources = field.conductors().to_vec();
            let drive = boundary(&field, 2.0, 2.0, None);
            let refused = prepared.step_spatial_conduction_candidate(
                cx,
                &thermal,
                0.5,
                &field,
                &drive,
                regions,
                1e-12,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
            );
            assert!(
                refused.is_err(),
                "a hot unsupported emissivity cannot use its 300 K slice"
            );
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(thermal.energy_j(), energy);
            assert_eq!(field.conductors(), sources);
            gate.request();
            assert!(matches!(
                prepared.step_spatial_conduction_candidate(
                    cx,
                    &thermal,
                    0.5,
                    &field,
                    &drive,
                    regions,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.energy_j(), energy);
            assert_eq!(field.conductors(), sources);
        },
    );
}
