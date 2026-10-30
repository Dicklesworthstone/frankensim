//! The same ceramic card supplies caloric, resistive and reciprocal field laws.
//! Coolant substitutions reach its temperature-dependent stiffness and coupling.
use super::*;
use fs_couple::circuit_enthalpy::piezoelectric::thermal::{
    PiezoelectricThermalBinding, PiezoelectricThermalError, heat as ceramic_heat,
};
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::{
    COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, LongitudinalPiezoelectric,
    PERMITTIVITY_33_T_PROPERTY, STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS,
    thermal::ThermalLongitudinalPiezoelectric,
};

const SLOPES: [f64; 3] = [0.0004, 0.0015, 0.006];
const CERAMIC_DENSITY: f64 = 3.0;
const CERAMIC_CP: f64 = 1.0 / 3.0;

fn ceramic_material(scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        let value = match claim.key.name() {
            SPECIFIC_HEAT_CAPACITY_PROPERTY => Some(CERAMIC_CP),
            DENSITY_PROPERTY => Some(CERAMIC_DENSITY),
            _ => None,
        };
        if let Some(value) = value {
            claim.value = PropertyValue::Scalar {
                value,
                dims: claim.key.quantity().dims(),
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, value, slope) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, 0.1, SLOPES[0]),
        (
            STRAIN_CHARGE_33_PROPERTY,
            STRAIN_CHARGE_DIMS,
            0.2,
            SLOPES[1],
        ),
        (
            PERMITTIVITY_33_T_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            2.0,
            SLOPES[2],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, value - 20.0 * scale * slope),
                        (320.0, value + 20.0 * scale * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-coolant ceramic fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &body_point(280.0),
            upper: &body_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn thermal_pair<'a>(
    ceramic: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    cp: f64,
) -> DepositionSimulation<'a> {
    let mut events = [
        (BODY, 0, CERAMIC_DENSITY, 1.0, 30.0 * CERAMIC_CP),
        (FLUID, 1, MASS, 1e-6, 10.0 * cp),
    ]
    .map(
        |(id, material, mass_kg, volume_m3, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: material,
                body_id: material,
                material_id: MaterialId(material),
                mass_kg,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, material as f64, 0.0],
                    end_m: [2.0, material as f64, 0.0],
                    volume_m3,
                    surface_area_m2: SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        },
    )
    .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: BODY,
            b: FLUID,
            area_m2: SURFACE,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), ceramic.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
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
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

// Independent 2x2 constitutive inverse and thermodynamic derivative.
// Reference volume is 1 m^3; the original consistent rod mass is 1 kg.
fn potentials(t: f64, x: f64, q: f64, scale: f64) -> [f64; 3] {
    let [sp, dp, ep] = SLOPES.map(|s| scale * s);
    let (s, d, e) = (
        0.1 + sp * (t - 300.0),
        0.2 + dp * (t - 300.0),
        2.0 + ep * (t - 300.0),
    );
    let (strain, displacement) = (x / 2.0, q / 0.5);
    let determinant = s * e - d * d;
    let stress = (e * strain - d * displacement) / determinant;
    let field = (s * displacement - d * strain) / determinant;
    let free = 0.5 * (strain * stress + displacement * field);
    let entropy = 0.5 * (sp * stress * stress + 2.0 * dp * stress * field + ep * field * field);
    [free, entropy, free + t * entropy]
}

fn midpoint(t: f64, old: [f64; 3], scale: f64) -> [f64; 5] {
    let [x, v, q] = old;
    let [sp, dp, ep] = SLOPES.map(|s| scale * s);
    let (s, d, e) = (
        0.1 + sp * (t - 300.0),
        0.2 + dp * (t - 300.0),
        2.0 + ep * (t - 300.0),
    );
    let (k, c, theta) = (0.25 / s, 0.25 * (e - d * d / s), 0.25 * d / s);
    let r = 1.0 + 0.02 * (t - 300.0);
    let denominator = 1.0 + DT * DT * k / 4.0;
    let a = (x + 0.5 * DT * v) / denominator;
    let b = DT * DT * theta / (4.0 * denominator);
    let voltage = (q - theta * a) / (c + DT / (2.0 * r) + theta * b);
    let middle = a + b * voltage;
    [
        2.0 * middle - x,
        4.0 * (middle - x) / DT - v,
        q - DT * voltage / r,
        voltage / r,
        DT * voltage * voltage / r,
    ]
}

// Conservation of caloric + ceramic INTERNAL + kinetic energy eliminates
// coolant T. The remaining original backward-Euler fluid heat equation is scalar.
fn reference(old: [f64; 5], viscosity: f64, cp: f64, scale: f64) -> [f64; 9] {
    let [t0, f0, x0, v0, q0] = old;
    let before = potentials(t0, x0, q0, scale);
    let evaluate = |t: f64| {
        let [x, v, q, current, joule] = midpoint(t, [x0, v0, q0], scale);
        let after = potentials(t, x, q, scale);
        let change = after[2] - before[2] + 0.5 * (v * v - v0 * v0);
        let fluid = f0 - ((t - t0) + change) / (MASS * cp);
        let film = 0.5 * (t + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (t - fluid);
        let frozen = potentials(t, x0, q0, scale);
        let reversible = -(frozen[0] - before[0]) - (t * after[1] - t0 * before[1]);
        (
            (t - t0) + change + DT * power,
            [t, fluid, x, v, q, current, joule, power, reversible],
        )
    };
    let (mut low, mut high) = (280.0, 320.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let t = 0.5 * (low + high);
        if evaluate(t).0 > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    let (residual, state) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    state
}

#[test]
fn finite_coolant_and_ceramic_sources_share_temperature_motion_and_internal_energy() {
    let run = |viscosity, cp, scale| {
        let (card, curve) = ceramic_material(scale);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = thermal_pair(&curve, &fluid_curve, cp);
        let fluid = liquid_law(&fluid_card);
        let source = LongitudinalPiezoelectric::resolve(
            &card,
            &context,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let ceramic = PiezoelectricThermalBinding {
            cell: BODY,
            law: ThermalLongitudinalPiezoelectric::new(&card, &source, "T", [280.0, 320.0])
                .unwrap(),
            maximum_abs_displacement_c_m2: 0.22,
            volume_tolerance_m3: 1e-12,
        };
        let cells = [
            MaterialCell {
                cell: BODY,
                material_id: MaterialId(0),
                curve: curve.curve(),
                point: ThermalMaterialPoint {
                    material: &card,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: FLUID,
                material_id: MaterialId(1),
                curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card,
                    context: &fluid_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
        let resistors = [CellResistor {
            destination: ResistorCell {
                branch: 1,
                cell: BODY,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 0.1,
            area_m2: 1e-5,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[(1, ElectricalPropertyBasis::Resistivity)],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let links = [binding()];
        let initial = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let prepared_bulk = bulk
            .prepare_with_link_conductances(
                &thermal,
                electrical.electrical_bases,
                initial.conductances_w_k(),
                || false,
            )
            .unwrap();
        let model = PiezoelectricCircuit::new(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Capacitor { farads: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            prepared_bulk
                .resolve_with_link_conductances(&thermal, initial.conductances_w_k(), || false)
                .unwrap()
                .resistors,
            vec![],
            vec![0.1],
            0,
            PiezoelectricRod::new(
                source,
                2.0,
                0.5,
                SmallSignalLimits {
                    max_abs_strain: 0.03,
                    max_abs_electric_field_v_m: 10.0,
                },
            )
            .unwrap(),
            0.02,
            0.01,
            &[],
            || false,
        )
        .unwrap();
        let before_cells = thermal.cells().clone();
        let before_model = model.clone();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 3];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 179,
                    kernel_id: 33,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 2;
            let (mut accepted_thermal, mut accepted_model) = (thermal.clone(), model.clone());
            let mut state = [310.0, 290.0, 0.02, 0.01, 0.1];
            for step in 1..=2 {
                let expected = reference(state, viscosity, cp, scale);
                let saved_cells = accepted_thermal.cells().clone();
                let saved_model = accepted_model.clone();
                let next = ceramic_heat::step_candidate(
                    &assembly,
                    &cx,
                    &accepted_thermal,
                    f64::from(step) * DT,
                    &accepted_model,
                    &[],
                    0.0,
                    &ceramic,
                    solve,
                    &[0.0, 0.0],
                )
                .unwrap();
                let report = next.report();
                let coupled = &report.coupled;
                for (actual, wanted) in [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.thermal().temperature_k(FLUID).unwrap(),
                    next.model().extension_m(),
                    next.model().velocity_m_s(),
                    next.model().charge_c(),
                ]
                .into_iter()
                .zip(expected)
                {
                    near(actual, wanted, 2e-8);
                }
                near(coupled.electrical.resistors[0].current_a, expected[5], 2e-9);
                near(coupled.deposited_joule_j, expected[6], 2e-10);
                near(coupled.source.reversible_heat_j, expected[8], 2e-10);
                near(
                    expected[0] - state[0],
                    expected[6] + expected[8] - DT * expected[7],
                    1e-9,
                );
                near(MASS * cp * (expected[1] - state[1]), DT * expected[7], 1e-9);
                near(
                    coupled.thermal.stored_energy_change_j,
                    coupled.deposited_joule_j + coupled.source.reversible_heat_j,
                    2e-9,
                );
                near(report.external_heat_j, 0.0, 1e-10);
                near(report.total_energy_residual_j, 0.0, 2e-9);
                near(report.closed_radiation_net_j, 0.0, 1e-14);
                let source = &coupled.source;
                let sample = &source.network.finite_convection.as_ref().unwrap().samples()[0];
                near(sample.heat_into_fluid_w, expected[7], 2e-7);
                near(
                    sample.property_temperature_k,
                    0.5 * (expected[0] + expected[1]),
                    2e-8,
                );
                assert_eq!(
                    sample.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(sample.material.material().query().axes()["p"], 100000.0);
                assert_eq!(
                    source.material.material().card_identity(),
                    card.content_hash()
                );
                assert_eq!(source.material.material().query().axes()["lot"], 7.0);
                assert_eq!(next.model().rod().material(), &source.material);
                near(source.temperature_k[1], expected[0], 2e-8);
                near(
                    source.state[1].internal_energy_j,
                    potentials(expected[0], expected[2], expected[4], scale)[2],
                    1e-9,
                );
                assert_eq!(accepted_thermal.cells(), &saved_cells);
                assert_eq!(accepted_model, saved_model);
                assert_eq!(next.thermal().links(), thermal.links());
                if scale == 0.0 {
                    assert_eq!(source.reversible_heat_j, 0.0);
                    assert_eq!(coupled.electrical.parameter_work_j, 0.0);
                } else {
                    assert!(source.reversible_heat_j.abs() > 1e-5);
                    assert!(coupled.electrical.parameter_work_j.abs() > 1e-6);
                }
                result = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    coupled.electrical.resistors[0].current_a,
                    next.model().extension_m(),
                ];
                state = [
                    expected[0],
                    expected[1],
                    expected[2],
                    expected[3],
                    expected[4],
                ];
                (accepted_thermal, accepted_model, _) = next.into_parts();
            }
            if viscosity == 1.0 && cp == 4000.0 && scale == 1.0 {
                let mut insufficient = solve;
                insufficient.maximum_circuit_trials = 1;
                assert!(matches!(
                    ceramic_heat::step_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        DT,
                        &model,
                        &[],
                        0.0,
                        &ceramic,
                        insufficient,
                        &[0.0, 0.0]
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
                let queries = std::cell::Cell::new(0);
                let refused = SourceConvectionHeatNetwork {
                    heat: &prepared,
                    cells: &[],
                    maximum_cells: 0,
                    resolve_fluid: |_: &SourceConvectiveCell,
                                    _: &DepositionSimulation<'_>,
                                    t: Temperature| {
                        Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                    },
                };
                let failing = refused.with_links(&links, 1, |_, _, t| {
                    queries.set(queries.get() + 1);
                    if queries.get() > 1 {
                        Err("ceramic coolant trial refused")
                    } else {
                        Ok(fluid.at_temperature(t).unwrap())
                    }
                });
                assert!(matches!(
                    ceramic_heat::step_candidate(
                        &failing,
                        &cx,
                        &thermal,
                        DT,
                        &model,
                        &[],
                        0.0,
                        &ceramic,
                        solve,
                        &[0.0, 0.0]
                    ),
                    Err(CircuitEnthalpyError::Source(
                        PiezoelectricThermalError::Heat(ConvectionHeatError::Convection(
                            SourceCellConvectionError::Source {
                                cell: FLUID,
                                source: "ceramic coolant trial refused",
                                ..
                            }
                        ))
                    ))
                ));
                assert!(queries.get() > 1);
                let retry = ceramic_heat::step_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    &model,
                    &[],
                    0.0,
                    &ceramic,
                    solve,
                    &[0.0, 0.0],
                )
                .unwrap();
                let repeated = ceramic_heat::step_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    &model,
                    &[],
                    0.0,
                    &ceramic,
                    solve,
                    &[0.0, 0.0],
                )
                .unwrap();
                assert_eq!(retry.model(), repeated.model());
                assert_eq!(retry.thermal().cells(), repeated.thermal().cells());
                gate.request();
                assert!(matches!(
                    ceramic_heat::step_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        DT,
                        &model,
                        &[],
                        0.0,
                        &ceramic,
                        solve,
                        &[0.0, 0.0]
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
            }
            assert_eq!(thermal.cells(), &before_cells);
            assert_eq!(model, before_model);
        });
        result
    };
    let baseline = run(1.0, 4000.0, 1.0);
    let viscous = run(4.0, 4000.0, 1.0);
    let larger_capacity = run(1.0, 8000.0, 1.0);
    let frozen_ceramic = run(1.0, 4000.0, 0.0);
    assert!(viscous[0] > baseline[0] && viscous[1] < baseline[1] && viscous[2] > baseline[2]);
    assert!(
        larger_capacity[0] < baseline[0]
            && larger_capacity[1] > baseline[1]
            && larger_capacity[2] < baseline[2]
    );
    assert!((baseline[2] - frozen_ceramic[2]).abs() > 1e-5);
}
