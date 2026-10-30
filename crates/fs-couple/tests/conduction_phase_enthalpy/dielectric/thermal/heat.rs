//! G1: independent two-body electrical/internal-energy and finite heat balance.
//! Source numbers and infinite-plate geometry are explicit synthetic fixtures.
use super::*;
#[path = "../../../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::{ConductionLink, MaterialNetwork},
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SelectionPolicy, SurfaceSpec, SystemContext};
use fs_material::{
    state_point::THERMAL_CONDUCTIVITY_DIMS, thermal::interface::ThermalInterfacePoint,
};

fn material(slope: Option<f64>, e: f64, upper: f64) -> MaterialCard {
    let base = source(slope.unwrap_or(0.0), 1.0);
    let mut claims = ClaimSet::new();
    for name in [
        ElectricalPropertyBasis::Resistivity.property_name(),
        ELECTRIC_PERMITTIVITY_PROPERTY,
        PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[0],
        PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1],
    ] {
        let mut claim = base.claims().claims_for(name)[0].1.clone();
        if slope.is_none() && name != ElectricalPropertyBasis::Resistivity.property_name() {
            claim.value = PropertyValue::Scalar {
                value: if name == PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1] {
                    4.0
                } else {
                    2.0
                },
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, value, derivative, high) in [
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            1.0,
            0.0,
            360.0,
        ),
        (
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            e,
            0.001 * e,
            upper,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, value - 20.0 * derivative),
                        (high, value + (high - 300.0) * derivative),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature(),
                    280.0,
                    high,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic dielectric heat assembly".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn with_heat<R>(
    slope: Option<f64>,
    e: f64,
    resistance: f64,
    upper: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        PlanarLeakyDielectric<ResolvedConductorState>,
    ) -> R,
) -> R {
    let cards = [material(slope, e, upper), material(slope, 0.6, 360.0)];
    let curves = [caloric(&cards[0], false), caloric(&cards[1], false)];
    let points = [context(300.0), context(320.0)];
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = [0, 1].map(|i| MaterialCell {
        cell: CellId(7 + i as u64),
        material_id: MaterialId(3 + i as u64),
        curve: &curves[i],
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &points[i],
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let mut events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: bindings[i].cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: bindings[i].material_id,
                mass_kg: 1.0,
                specific_h_j_kg: (20.0 + 20.0 * i as f64) * CAPACITY,
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
                    source_w: if i == 0 { 0.4 } else { 0.0 },
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 9,
            a: CellId(7),
            b: CellId(8),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(3), &curves[0] as &dyn ThermalLaw),
            (MaterialId(4), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let links = [ConductionLink {
        link: 9,
        distance_a_m: 0.25,
        distance_b_m: 0.25,
    }];
    let bulk = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &links,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Curve {
                abscissa: "aT".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.98 * resistance), (360.0, 1.06 * resistance)],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("aT", absolute_temperature(), 280.0, 360.0)
                .with_quantity("bT", absolute_temperature(), 280.0, 360.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic dielectric thermal interface".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "a".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "b".into(),
        },
        SystemContext {
            medium: "fixed interface".into(),
            third_body: None,
            environment: "numerical fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let interface_point = QueryPoint::new()
        .with_quantity("aT", absolute_temperature(), 300.0)
        .unwrap()
        .with_quantity("bT", absolute_temperature(), 320.0)
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 9,
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &interface_point,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        sides: [0, 1].map(|i| ContactCell {
            cell: bindings[i].cell,
            material_id: bindings[i].material_id,
            curve: &curves[i],
        }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let surfaces = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: bindings[i].cell,
        material_id: bindings[i].material_id,
        curve: &curves[i],
        material: bindings[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let emissivities = [0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original dielectric enclosure",
            bindings[i].point,
            curves[i]
                .state_at_specific_enthalpy((20.0 + 20.0 * i as f64) * CAPACITY)
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
    let field = if slope.is_some() {
        thermal_field(&cards[0], None, false, [280.0, 360.0])
    } else {
        let conduction = PlanarConduction::new(
            mesh(),
            vec![source_conductor(&cards[0], &points[0]); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        PlanarLeakyDielectric::new(
            conduction,
            vec![0.0, 4.0, 4.0, 0.0, 2.0],
            &[0, 3],
            0.0,
            limits(),
            electric_controls(),
            &mut || false,
        )
        .unwrap()
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 23,
                    kernel_id: 19,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
            &assembly,
            thermal,
            field,
        )
    })
}

// Independent scalar charge equation + whole-system energy eliminates the
// passive temperature. Neither source evaluators nor FE matrices define this.
fn reference(t0: f64, b0: f64, v0: f64, slope: f64, e: f64, resistance: f64) -> [f64; 3] {
    let state = |t: f64| {
        let sigma = 1.0 / (1.0 + 0.01 * (t - 300.0));
        let v =
            v0 * (epsilon(t0, 2.0, slope) - 0.25 * sigma) / (epsilon(t, 2.0, slope) + 0.25 * sigma);
        let delta_u = field_u(t, v, 2.0, slope) - field_u(t0, v0, 2.0, slope);
        let b = b0 + (0.2 - CAPACITY * (t - t0) - delta_u) / CAPACITY;
        let contact = (t - b) / (0.5 + resistance * (1.0 + 0.001 * (t - 300.0)));
        let radiation = 5.670_374_419e-8 * (t.powi(4) - b.powi(4))
            / (1.0 / (e * (1.0 + 0.001 * (t - 300.0))) + 1.0 / (0.6 * (1.0 + 0.001 * (b - 300.0)))
                - 1.0);
        (
            b,
            v,
            CAPACITY * (t - t0) + delta_u - 0.2 + 0.5 * (contact + radiation),
        )
    };
    let t = bisect(290.0, 330.0, |t| state(t).2);
    [t, state(t).0, state(t).1]
}

#[test]
fn g1_dielectric_charge_and_complete_heat_follow_independent_two_body_balance() {
    let mut outputs = Vec::new();
    for slope in [None, Some(0.004)] {
        for (e, resistance) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
            with_heat(
                slope,
                e,
                resistance,
                360.0,
                |cx, _, assembly, thermal, field| {
                    let prepared = assembly.prepare(cx, &thermal).unwrap();
                    let (mut heat, mut charge) = (thermal, field);
                    let (mut t, mut b, mut v) = (300.0, 320.0, 4.0);
                    for step in 1..=2 {
                        let expected = reference(t, b, v, slope.unwrap_or(0.0), e, resistance);
                        let target = 0.5 * f64::from(step);
                        let (next_heat, next_charge, external, source) = if slope.is_some() {
                            let next = prepared
                                .step_dielectric_thermal_candidate(
                                    cx,
                                    &heat,
                                    target,
                                    &charge,
                                    &[(0, 0.0), (3, 0.0)],
                                    &[24.0; 4],
                                    &ROUTES,
                                    1e-12,
                                    limits(),
                                    electric_controls(),
                                    config(),
                                    &[0.0; 2],
                                )
                                .unwrap();
                            let (h, f, r) = next.into_parts();
                            near(r.total_energy_residual_j, 0.0);
                            near(
                                r.coupled.source.field_internal_energy_j[1],
                                field_u(expected[0], expected[2], 2.0, 0.004),
                            );
                            assert!(r.coupled.source.reversible_heat_j.abs() > 1e-7);
                            (h, f, r.external_heat_j, r.coupled.source.material)
                        } else {
                            let next = prepared
                                .step_dielectric_candidate(
                                    cx,
                                    &heat,
                                    target,
                                    &charge,
                                    &[(0, 0.0), (3, 0.0)],
                                    &ROUTES,
                                    1e-12,
                                    limits(),
                                    electric_controls(),
                                    config(),
                                    &[0.0; 2],
                                )
                                .unwrap();
                            let (h, f, r) = next.into_parts();
                            let external =
                                r.non_electrical_heat_j - 0.5 * r.source.radiation_net_power_w();
                            (h, f, external, r.source)
                        };
                        near(next_heat.temperature_k(CellId(7)).unwrap(), expected[0]);
                        near(next_heat.temperature_k(CellId(8)).unwrap(), expected[1]);
                        near(next_charge.voltage_v()[1], expected[2]);
                        near(external, 0.2);
                        near(
                            source.thermal_contacts[0].interface.value_m2_k_per_w(),
                            resistance * (1.0 + 0.001 * (expected[0] - 300.0)),
                        );
                        for surface in source.radiation.as_ref().unwrap().surfaces() {
                            let original = source
                                .material
                                .cells
                                .iter()
                                .find(|s| s.cell == surface.cell)
                                .unwrap();
                            assert_eq!(surface.phase, original.phase);
                        }
                        near(
                            next_charge.conduction().conductors()[0]
                                .material()
                                .query()
                                .axes()["T"],
                            expected[0],
                        );
                        assert_eq!(next_charge.time_s(), next_heat.time_s());
                        (t, b, v) = (expected[0], expected[1], expected[2]);
                        (heat, charge) = (next_heat, next_charge);
                    }
                    outputs.push([t, b, v]);
                },
            );
        }
    }
    for start in [0, 3] {
        for variant in [start + 1, start + 2] {
            for axis in 0..3 {
                assert!((outputs[start][axis] - outputs[variant][axis]).abs() > 1e-6);
            }
        }
    }
}

#[test]
fn dielectric_heat_late_source_refusal_and_cancellation_preserve_both_histories() {
    for slope in [None, Some(0.004)] {
        with_heat(
            slope,
            0.2,
            0.2,
            300.1,
            |cx, gate, assembly, thermal, field| {
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let before = field.nodal_charge_c().to_vec();
                if slope.is_some() {
                    assert!(
                        prepared
                            .step_dielectric_thermal_candidate(
                                cx,
                                &thermal,
                                0.5,
                                &field,
                                &[(0, 0.0), (3, 0.0)],
                                &[24.0; 4],
                                &ROUTES,
                                1e-12,
                                limits(),
                                electric_controls(),
                                config(),
                                &[0.0; 2]
                            )
                            .is_err()
                    );
                } else {
                    assert!(
                        prepared
                            .step_dielectric_candidate(
                                cx,
                                &thermal,
                                0.5,
                                &field,
                                &[(0, 0.0), (3, 0.0)],
                                &ROUTES,
                                1e-12,
                                limits(),
                                electric_controls(),
                                config(),
                                &[0.0; 2]
                            )
                            .is_err()
                    );
                }
                assert_eq!(field.nodal_charge_c(), before.as_slice());
                assert_eq!(field.time_s(), 0.0);
                assert_eq!(thermal.time_s(), 0.0);
                near(thermal.temperature_k(CellId(7)).unwrap(), 300.0);
                gate.request();
                assert!(matches!(
                    prepared.step_dielectric_candidate(
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        &[(0, 0.0), (3, 0.0)],
                        &ROUTES,
                        1e-12,
                        limits(),
                        electric_controls(),
                        config(),
                        &[0.0; 2]
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
            },
        );
    }
}
