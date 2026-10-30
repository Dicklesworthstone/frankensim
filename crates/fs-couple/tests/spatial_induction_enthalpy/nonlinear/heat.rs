//! G1/G3/G4: complete sourced thermal transport reaches both native 3D histories.
//! Independent two-temperature balances use no field, heat or radiosity solver.
use super::super as linear;
use super::*;
#[path = "../../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    magnetic::MagneticThermalError,
    material::{ConductionLink, MaterialNetwork},
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::*;
use fs_material::{
    state_point::THERMAL_CONDUCTIVITY_DIMS, thermal::interface::ThermalInterfacePoint,
};
use fs_qty::{Dims, Temperature};

const CELL: CellId = CellId(7);
const SINK: CellId = CellId(8);
const DT: f64 = 0.05;
const A0: f64 = 1.2;
const RHO: f64 = 0.2;
const ALPHA: f64 = 0.009;
const MU_SLOPE: f64 = 0.0004;
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;

fn claim(name: &str, dims: Dims, value: PropertyValue, upper: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        interpolation: if matches!(value, PropertyValue::Curve { .. }) {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        },
        value,
        validity: ValidityDomain::unconstrained().with("T", 200.0, upper),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic spatial induction thermal assembly".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn with_heat<R>(
    nonlinear: bool,
    epsilon: f64,
    resistance: f64,
    upper: f64,
    magnetic_restriction: Option<Dims>,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        &MaterialCard,
    ) -> R,
) -> R {
    let base = if nonlinear {
        super::card(RHO, ALPHA, BASIS, true)
    } else {
        source::absolute_affine_magnetic_card(RHO, ALPHA, 1.0, MU_SLOPE, BASIS)
    };
    let cards = [
        base,
        source::absolute_affine_magnetic_card(0.5, 0.0, 1.0, 0.0, BASIS),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, base)| {
        let mut claims = base.claims().clone();
        claims
            .insert_claim(claim(
                "thermal_conductivity",
                THERMAL_CONDUCTIVITY_DIMS,
                PropertyValue::Scalar {
                    value: 1.0,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
                500.0,
            ))
            .unwrap();
        let e = if i == 0 { epsilon } else { 0.6 };
        let limit = if i == 0 { upper } else { 500.0 };
        let mut emissivity = claim(
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [200.0, limit]
                    .map(|t| (t, e * (1.0 + 0.0001 * (t - 300.0))))
                    .to_vec(),
                dims: EMISSIVITY_DIMS,
            },
            limit,
        );
        if let Some(dims) = magnetic_restriction {
            emissivity.validity = emissivity.validity.with_quantity(
                "field",
                QuantitySpec::dimensional(dims),
                0.0,
                4.0,
            );
        }
        claims.insert_claim(emissivity).unwrap();
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
    })
    .collect::<Vec<_>>();
    let curves = cards.iter().map(linear::curve).collect::<Vec<_>>();
    let mut context = if nonlinear {
        super::context(300.0)
    } else {
        source::absolute_point(300.0)
    };
    if let Some(dims) = magnetic_restriction {
        context = context
            .with_quantity("field", QuantitySpec::dimensional(dims), 0.0)
            .unwrap();
    }
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: if i == 0 { CELL } else { SINK },
        material_id: MaterialId(9 + i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let links = [ConductionLink {
        link: 17,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &links,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let mut interface_claim = claim(
        RP,
        RD,
        PropertyValue::Curve {
            abscissa: "aT".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.9 * resistance), (500.0, 1.2 * resistance)],
            dims: RD,
        },
        500.0,
    );
    interface_claim.validity = ValidityDomain::unconstrained()
        .with("aT", 200.0, 500.0)
        .with("bT", 200.0, 500.0);
    let mut interface_claims = ClaimSet::new();
    interface_claims.insert_claim(interface_claim).unwrap();
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
            medium: "fixed contact".into(),
            third_body: None,
            environment: "authored".into(),
            history: "fixed".into(),
        },
        interface_claims,
        vec![],
    )
    .unwrap();
    let interface_context = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 340.0)
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 17,
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
    let initial_h = [400.0, 560.0];
    let enclosure = radiative::enclosure([0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original spatial material surface",
            cells[i].point,
            curves[i]
                .curve()
                .state_at_specific_enthalpy(initial_h[i])
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
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let mut events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cells[i].cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: cells[i].material_id,
                mass_kg: 1.0,
                specific_h_j_kg: initial_h[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 17,
            a: CELL,
            b: SINK,
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(9), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(10), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 97,
                kernel_id: 16,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal, &cards[0])
    })
}

fn transfer(ta: f64, tb: f64, epsilon: f64, contact: f64) -> f64 {
    let g = 1.0 / (1.0 + contact * (1.0 + 0.001 * (ta - 300.0)));
    let radiated = 5.670_374_419e-8 * (tb.powi(4) - ta.powi(4))
        / (1.0 / (epsilon * (1.0 + 0.0001 * (ta - 300.0)))
            + 1.0 / (0.6 * (1.0 + 0.0001 * (tb - 300.0)))
            - 1.0);
    g * (tb - ta) + radiated
}

fn reference(
    t0: [f64; 2],
    a0: f64,
    initial_drive: f64,
    nonlinear: bool,
    epsilon: f64,
    contact: f64,
) -> ([f64; 2], f64, f64) {
    let mu = |t| 1.0 + MU_SLOPE * (t - 300.0);
    let internal = |t, a: f64| {
        if nonlinear {
            energy(a, true)
        } else {
            2.0 * a * a / mu(t) * (1.0 + t * MU_SLOPE / mu(t))
        }
    };
    let amplitude = |t| {
        let m = 3.0 * sigma(t, RHO, ALPHA, BASIS) / 40.0;
        if nonlinear {
            bracket(0.0, a0, |a| m * (a - a0) / DT + 2.0 * h(2.0 * a, true))
        } else {
            (m * a0 + 0.5 * DT * (initial_drive - 4.0 * a0 / mu(t0[0]))) / (m + 2.0 * DT / mu(t))
        }
    };
    let sink = |ta| {
        bracket(295.0, 350.0, |tb| {
            4.0 * (tb - t0[1]) + DT * transfer(ta, tb, epsilon, contact)
        })
    };
    let ta = bracket(295.0, 350.0, |ta| {
        let a = amplitude(ta);
        let d = if nonlinear {
            2.0 * h(2.0 * a, true) * (a - a0) - energy(a, true) + energy(a0, true)
        } else {
            0.0
        };
        4.0 * (ta - t0[0]) + 4.0 * (sink(ta) - t0[1]) + internal(ta, a) - internal(t0[0], a0) + d
            - if nonlinear {
                0.0
            } else {
                0.5 * initial_drive * (a - a0)
            }
    });
    let a = amplitude(ta);
    (
        [ta, sink(ta)],
        a,
        3.0 * sigma(ta, RHO, ALPHA, BASIS) / 40.0 * (a - a0).powi(2) / DT,
    )
}

fn source_checks<D>(
    report: &fs_couple::circuit_enthalpy::spatial_induction::SpatialInductionThermalReport<
        D,
        HeatNetworkSample,
    >,
    thermal: &DepositionSimulation<'_>,
    expected: [f64; 2],
    expected_heat: f64,
) {
    near(thermal.temperature_k(CELL).unwrap(), expected[0], 3e-8);
    near(thermal.temperature_k(SINK).unwrap(), expected[1], 3e-8);
    near(report.coupled.deposited_joule_j, expected_heat, 3e-8);
    near(report.external_heat_j, 0.0, 2e-8);
    near(report.total_energy_residual_j, 0.0, 2e-8);
    let source = &report.coupled.source.network;
    assert_eq!(source.thermal_contacts.len(), 1);
    for surface in source.radiation.as_ref().unwrap().surfaces() {
        let body = source
            .material
            .cells
            .iter()
            .find(|cell| cell.cell == surface.cell)
            .unwrap();
        assert_eq!(body.phase, surface.phase);
    }
}

#[test]
fn interface_emissivity_and_magnetic_entropy_share_the_original_spatial_history() {
    let mut outputs = Vec::new();
    for (epsilon, contact) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
        outputs.push(with_heat(
            false,
            epsilon,
            contact,
            500.0,
            None,
            |cx, _, assembly, mut thermal, card| {
                let problem = linear::problem(card, BASIS);
                let mut field = linear::initial(&problem, BASIS, A0, 1.0);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let laws = linear::laws(card, &problem, 3.0);
                for step in 1..=2 {
                    let a0 = field.potential_wb()[linear::diagonal(&problem)];
                    let expected = reference(
                        [
                            thermal.temperature_k(CELL).unwrap(),
                            thermal.temperature_k(SINK).unwrap(),
                        ],
                        a0,
                        if step == 1 { 4.0 * A0 } else { 0.0 },
                        false,
                        epsilon,
                        contact,
                    );
                    let next = prepared
                        .step_spatial_induction_candidate(
                            cx,
                            &thermal,
                            step as f64 * DT,
                            &field,
                            &[[0.0; 3]; 6],
                            &[[0.0; 3]; 8],
                            &linear::regions(),
                            1e-12,
                            &laws,
                            thermal_controls(),
                            linear::field_controls(),
                            linear::assembly(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    source_checks(next.report(), next.thermal(), expected.0, expected.2);
                    near(
                        next.field().potential_wb()[linear::diagonal(&problem)],
                        expected.1,
                        3e-9,
                    );
                    assert!(next.report().coupled.source.reversible_heat_j.abs() > 1e-8);
                    let (next_thermal, next_field, _) = next.into_parts();
                    thermal = next_thermal;
                    field = next_field;
                }
                [
                    thermal.temperature_k(CELL).unwrap(),
                    thermal.temperature_k(SINK).unwrap(),
                    field.potential_wb()[linear::diagonal(&problem)],
                ]
            },
        ));
    }
    for changed in &outputs[1..] {
        for i in 0..3 {
            assert!((outputs[0][i] - changed[i]).abs() > 1e-6);
        }
    }
}

#[test]
fn source_contact_and_radiation_change_saturation_history_without_numerical_heating() {
    let mut outputs = Vec::new();
    for (epsilon, contact) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
        outputs.push(with_heat(
            true,
            epsilon,
            contact,
            500.0,
            None,
            |cx, _, assembly, mut thermal, card| {
                let problem = super::problem(card, BASIS);
                let mut field = super::initial(&problem, card, BASIS, A0, true);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                for step in 1..=2 {
                    let expected = reference(
                        [
                            thermal.temperature_k(CELL).unwrap(),
                            thermal.temperature_k(SINK).unwrap(),
                        ],
                        field.potential_wb()[super::diagonal(&problem)],
                        0.0,
                        true,
                        epsilon,
                        contact,
                    );
                    let before = thermal.energy_j();
                    let next = prepared
                        .step_spatial_nonlinear_induction_candidate(
                            cx,
                            &thermal,
                            step as f64 * DT,
                            &field,
                            &[[0.0; 3]; 6],
                            &[[0.0; 3]; 8],
                            &linear::regions(),
                            1e-12,
                            thermal_controls(),
                            super::field_controls(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    source_checks(next.report(), next.thermal(), expected.0, expected.2);
                    near(
                        next.field().potential_wb()[super::diagonal(&problem)],
                        expected.1,
                        3e-9,
                    );
                    near(next.thermal().energy_j() - before, expected.2, 3e-8);
                    assert!(next.report().coupled.electrical.numerical_dissipation_j > 1e-5);
                    let (next_thermal, next_field, _) = next.into_parts();
                    thermal = next_thermal;
                    field = next_field;
                }
                [
                    thermal.temperature_k(CELL).unwrap(),
                    thermal.temperature_k(SINK).unwrap(),
                    field.potential_wb()[super::diagonal(&problem)],
                ]
            },
        ));
    }
    for changed in &outputs[1..] {
        for i in 0..3 {
            assert!((outputs[0][i] - changed[i]).abs() > 1e-6);
        }
    }
}

#[test]
fn unsupported_surface_fields_source_exit_and_cancel_publish_neither_history() {
    for restriction in [Some(FLUX_DENSITY_DIMS), Some(FIELD_STRENGTH_DIMS), None] {
        with_heat(
            true,
            0.2,
            0.2,
            if restriction.is_none() { 300.0 } else { 500.0 },
            restriction,
            |cx, gate, assembly, thermal, card| {
                let problem = super::problem(card, BASIS);
                let field = super::initial(&problem, card, BASIS, A0, true);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let before = field.potential_wb().to_vec();
                let run = || {
                    prepared.step_spatial_nonlinear_induction_candidate(
                        cx,
                        &thermal,
                        DT,
                        &field,
                        &[[0.0; 3]; 6],
                        &[[0.0; 3]; 8],
                        &linear::regions(),
                        1e-12,
                        thermal_controls(),
                        super::field_controls(),
                        &[0.0; 2],
                    )
                };
                if restriction.is_some() {
                    assert!(matches!(
                        run(),
                        Err(CircuitEnthalpyError::Source(MagneticThermalError::Heat(
                            HeatNetworkError::Invalid(_)
                        )))
                    ));
                } else {
                    assert!(run().is_err());
                }
                gate.request();
                assert!(matches!(run(), Err(CircuitEnthalpyError::Cancelled)));
                assert_eq!(field.potential_wb(), before);
                assert_eq!(field.time_s(), 0.0);
                assert_eq!(thermal.time_s(), 0.0);
                near(thermal.temperature_k(CELL).unwrap(), 300.0, 1e-12);
                near(thermal.temperature_k(SINK).unwrap(), 340.0, 1e-12);
            },
        );
    }
}
