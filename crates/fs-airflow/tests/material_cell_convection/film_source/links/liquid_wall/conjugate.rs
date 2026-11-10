//! Original material wall -> actual liquid-facing surface -> reciprocal heat.
use super::*;

#[path = "conjugate/implicit.rs"]
mod implicit;
use fs_airflow::graph::thermal::medium::convection::deposition::source::{
    heat::conjugate::{SourceConjugateHeatAssembly, SourceConjugateHeatError},
    links::conjugate::{
        PreparedSourceConjugateLinks, SourceConjugateError, SourceConjugateLinkBinding,
    },
};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::{
        heat::{HeatAssembly, MaterialHeatNetwork},
        material::ElectricalContactNetwork,
    },
    material::{MaterialCell, MaterialNetwork},
};
use fs_matdb::ClaimSelection;
use fs_material::thermal::{
    conduction::{LayerConductionConfig, ReservoirLayerGeometry, source::ResolvedConductionLayer},
    source::ThermalMaterialPoint,
};
use std::sync::Arc;

const DISTANCE: f64 = 0.0002;

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}

fn wall_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}

// Synthetic source data: k = scale * [0.6 + 0.015(T-300)], not a fluid/solid fit.
fn wall_source(scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, [1000.0, 1000.0]),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [CP, CP],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            [0.3 * scale, 0.9 * scale],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, values[0]), (320.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    280.0,
                    320.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured source-integrated wall".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
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
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &wall_point(280.0),
            upper: &wall_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-8,
    )
    .unwrap();
    (card, curve)
}

fn wall_layer<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    point: &'a QueryPoint,
    span: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        span,
        || Ok(()),
    )
    .unwrap()
}

fn binding<'a, 'm>(layer: &'a ResolvedConductionLayer<'m>) -> SourceConjugateLinkBinding<'a, 'm> {
    SourceConjugateLinkBinding {
        convection: corrected_link(),
        wall_layer: layer,
        geometry: ReservoirLayerGeometry::planar(LINK_AREA, DISTANCE).unwrap(),
    }
}

fn pair_at<'a>(
    wall: &'a EquilibriumEnthalpyPhaseCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    temperatures: [f64; 2],
    reversed: bool,
) -> DepositionSimulation<'a> {
    let original = pair(wall, fluid, reversed, 0.0);
    let mut events = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = if cell.id == IDS[0] {
                CP * (temperatures[0] - 300.0)
            } else {
                FLUID_CP * (temperatures[1] - 300.0)
            };
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect::<Vec<_>>();
    events.extend(original.links().values().cloned().map(|link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    }));
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), wall as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 200,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn solid_power(wall: f64, surface: f64, scale: f64) -> f64 {
    LINK_AREA / DISTANCE
        * scale
        * (0.6 * (wall - surface) + 0.0075 * ((wall - 300.0).powi(2) - (surface - 300.0).powi(2)))
}

// Independent quadratic conductivity integral and the parent's independent
// Gnielinski section formula. This oracle never calls a production surface solver.
fn surface_oracle(wall: f64, bulk: f64, scale: f64) -> (f64, f64) {
    let (mut lo, mut hi) = (wall.min(bulk), wall.max(bulk));
    for _ in 0..90 {
        let surface = 0.5 * lo + 0.5 * hi;
        let residual = solid_power(wall, surface, scale)
            - LINK_AREA * coefficient(bulk, surface, true) * (surface - bulk);
        if residual > 0.0 {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * lo + 0.5 * hi;
    (
        surface,
        LINK_AREA * coefficient(bulk, surface, true) * (surface - bulk),
    )
}

#[test]
fn actual_surface_closes_variable_material_and_wall_pr_heat_in_both_directions() {
    let (solid, solid_caloric) = wall_source(1.0);
    let point = wall_point(300.0);
    let layer = wall_layer(&solid, solid_caloric.curve(), &point, [280.0, 320.0]);
    let (fluid_card, fluid_caloric) = fluid_source(1.0);
    let liquid = liquid_law(&fluid_card);
    for temperatures in [[310.0, 290.0], [290.0, 310.0], [300.0, 300.0]] {
        let thermal = pair_at(solid_caloric.curve(), &fluid_caloric, temperatures, false);
        let owner =
            PreparedSourceConjugateLinks::new(&thermal, &[binding(&layer)], 1, || Ok(())).unwrap();
        let result = owner
            .resolve(
                &thermal,
                config(),
                |_, _, t| liquid.at_temperature(t),
                || Ok(()),
            )
            .unwrap();
        let sample = &result.samples()[0];
        let (surface, power) = surface_oracle(temperatures[0], temperatures[1], 1.0);
        near(sample.layer.surface_temperature_k, surface, 3e-10);
        near(sample.layer.power_to_bath_w, power, 3e-9);
        near(sample.applied_heat_into_fluid_w, power, 3e-9);
        near(
            sample.layer.solid_power_to_bath_w,
            solid_power(temperatures[0], surface, 1.0),
            3e-9,
        );
        assert!(sample.layer.heat_residual_w.abs() < 3e-9);
        assert!(sample.representation_error_w.abs() <= sample.representation_allowance_w);
        assert!(Arc::ptr_eq(&sample.source_line, layer.source_line()));
        let film = &sample.layer.film;
        assert_eq!(sample.wall_temperature_k, temperatures[0]);
        assert_eq!(film.fluid_temperature_k, temperatures[1]);
        near(film.wall_temperature_k, surface, 3e-10);
        near(
            film.surface.as_ref().unwrap().thermal.temperature_k(),
            surface,
            3e-10,
        );
        assert_eq!(
            film.material.material().card_identity(),
            fluid_card.content_hash()
        );
        assert_eq!(
            film.surface
                .as_ref()
                .unwrap()
                .material
                .material()
                .card_identity(),
            fluid_card.content_hash()
        );
        near(
            film.conductance_w_k,
            LINK_AREA * coefficient(temperatures[1], surface, true),
            3e-9,
        );
        assert!(result.conductances_w_k()[&5] < film.conductance_w_k);
        if temperatures[0] == temperatures[1] {
            let solid_g = LINK_AREA * 0.6 / DISTANCE;
            let expected = 1.0 / (1.0 / solid_g + 1.0 / film.conductance_w_k);
            near(result.conductances_w_k()[&5], expected, 2e-12);
            assert_eq!(sample.layer.iterations, 0);
        } else {
            assert!((surface - temperatures[0]).abs() > 1.0);
            let bare = LINK_AREA
                * coefficient(temperatures[1], temperatures[0], true)
                * (temperatures[0] - temperatures[1]);
            assert!(power.abs() < 0.75 * bare.abs());
        }
    }
}

fn endpoint(scale: f64, dt: f64) -> [f64; 2] {
    let (mut lo, mut hi) = (294.0, 310.0);
    for _ in 0..90 {
        let wall = 0.5 * lo + 0.5 * hi;
        let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let power = surface_oracle(wall, bulk, scale).1;
        if MASS * CP * (wall - 310.0) + dt * power > 0.0 {
            hi = wall;
        } else {
            lo = wall;
        }
    }
    let wall = 0.5 * lo + 0.5 * hi;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn native_heat_assembly_closes_both_finite_stores_and_responds_to_original_wall_conductivity() {
    type Error = CircuitEnthalpyError<SourceConjugateHeatError<MaterialStatePointError>>;
    let dt = 0.2;
    let (fluid_card, fluid_caloric) = fluid_source(1.0);
    let liquid = liquid_law(&fluid_card);
    let fluid_point = context(290.0);
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let mut walls = Vec::new();
    for scale in [1.0, 4.0] {
        let (solid_card, solid_caloric) = wall_source(scale);
        let solid_point = wall_point(310.0);
        let layer = wall_layer(
            &solid_card,
            solid_caloric.curve(),
            &solid_point,
            [280.0, 320.0],
        );
        let expected = endpoint(scale, dt);
        for reversed in [false, true] {
            let thermal = pair(solid_caloric.curve(), &fluid_caloric, reversed, 0.0);
            let saved = thermal.cells().clone();
            let owner =
                PreparedSourceConjugateLinks::new(&thermal, &[binding(&layer)], 1, || Ok(()))
                    .unwrap();
            let initial = owner
                .resolve(
                    &thermal,
                    config(),
                    |_, _, t| liquid.at_temperature(t),
                    || Ok(()),
                )
                .unwrap();
            let cells = [
                MaterialCell {
                    cell: IDS[0],
                    material_id: MaterialId(0),
                    curve: solid_caloric.curve(),
                    point: ThermalMaterialPoint {
                        material: &solid_card,
                        context: &solid_point,
                        temperature_axis: "T",
                        liquid_fraction_axis: None,
                    },
                    conductivity: None,
                },
                MaterialCell {
                    cell: IDS[1],
                    material_id: MaterialId(1),
                    curve: fluid_caloric.curve(),
                    point: ThermalMaterialPoint {
                        material: &fluid_card,
                        context: &fluid_point,
                        temperature_axis: "T",
                        liquid_fraction_axis: None,
                    },
                    conductivity: None,
                },
            ];
            let bulk = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            pool.scope(|arena| {
                let cx = fs_exec::Cx::new(
                    &gate,
                    arena,
                    fs_exec::StreamKey {
                        seed: 173,
                        kernel_id: 29,
                        tile: 0,
                        iteration: 0,
                    },
                    fs_exec::Budget::INFINITE,
                    fs_exec::ExecMode::Deterministic,
                );
                let prepared = heat
                    .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                    .unwrap();
                let assembly = SourceConjugateHeatAssembly {
                    heat: &prepared,
                    links: &owner,
                    layer_config: config(),
                    resolve_fluid: |_: &SourceConvectiveLink,
                                    _: &DepositionSimulation<'_>,
                                    t: Temperature| {
                        liquid.at_temperature(t)
                    },
                };
                let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
                let increments = step
                    .cell_ids()
                    .iter()
                    .map(|&id| {
                        if id == IDS[0] {
                            MASS * CP * (expected[0] - 310.0)
                        } else {
                            MASS * FLUID_CP * (expected[1] - 290.0)
                        }
                    })
                    .collect::<Vec<_>>();
                let trial = step
                    .evaluate(
                        &increments,
                        |state| {
                            let coefficients = assembly.resolve_heat(&cx, state)?;
                            assert!(
                                coefficients
                                    .bulk
                                    .sources_w
                                    .values()
                                    .all(|&power| power == 0.0)
                            );
                            assert_eq!(coefficients.bulk.conductances_w_k.len(), 1);
                            Ok::<_, Error>(StepCoefficients {
                                sources_w: coefficients.bulk.sources_w,
                                conductances_w_k: coefficients.bulk.conductances_w_k,
                                payload: coefficients.bulk.source,
                            })
                        },
                        || Ok(()),
                    )
                    .unwrap();
                let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
                let sample = &source.conjugate.samples()[0];
                for (cell, temperature) in IDS.into_iter().zip(expected) {
                    near(next.temperature_k(cell).unwrap(), temperature, 1e-9);
                }
                near(
                    MASS * CP * (expected[0] - 310.0) + dt * sample.applied_heat_into_fluid_w,
                    0.0,
                    1e-8,
                );
                near(
                    MASS * FLUID_CP * (expected[1] - 290.0) - dt * sample.applied_heat_into_fluid_w,
                    0.0,
                    1e-8,
                );
                near(next.energy_j(), thermal.energy_j(), 1e-8);
                near(audit.external_heat_j, 0.0, 1e-12);
                near(audit.residual_j, 0.0, 1e-8);
                assert!(
                    assembly
                        .explicit_source_conductance_bounds_w_k(&source, &mut || Ok(()),)
                        .unwrap()
                        .is_empty(),
                    "finite-link G must not appear again as a source-rate bound"
                );
                assert_eq!(thermal.cells(), &saved);
            });
        }
        walls.push(expected[0]);
    }
    assert!(
        walls[1] < walls[0] - 1.0,
        "sourced wall conductivity must change actual finite-body cooling"
    );
}

#[test]
fn original_caloric_area_contact_source_support_and_cancellation_are_retained() {
    let (solid, solid_caloric) = wall_source(1.0);
    let point = wall_point(300.0);
    let layer = wall_layer(&solid, solid_caloric.curve(), &point, [280.0, 320.0]);
    let (fluid_card, fluid_caloric) = fluid_source(1.0);
    let liquid = liquid_law(&fluid_card);
    let thermal = pair(solid_caloric.curve(), &fluid_caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let wrong_area = SourceConjugateLinkBinding {
        geometry: ReservoirLayerGeometry::planar(2.0 * LINK_AREA, DISTANCE).unwrap(),
        ..binding(&layer)
    };
    assert!(PreparedSourceConjugateLinks::new(&thermal, &[wrong_area], 1, || Ok(())).is_err());
    let copied = solid_caloric.curve().clone();
    let foreign = wall_layer(&solid, &copied, &point, [280.0, 320.0]);
    assert!(matches!(
        PreparedSourceConjugateLinks::new(&thermal, &[binding(&foreign)], 1, || Ok(())),
        Err(SourceConjugateError::Invalid(
            "conjugate wall requires its original geometry and caloric object"
        ))
    ));
    let owner =
        PreparedSourceConjugateLinks::new(&thermal, &[binding(&layer)], 1, || Ok(())).unwrap();
    let changed = pair(solid_caloric.curve(), &fluid_caloric, true, 0.0);
    let mut calls = 0;
    assert!(
        owner
            .resolve(
                &changed,
                config(),
                |_, _, t| {
                    calls += 1;
                    liquid.at_temperature(t)
                },
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(
        calls, 0,
        "changed native contact must refuse before property queries"
    );
    let short = wall_layer(&solid, solid_caloric.curve(), &point, [305.0, 320.0]);
    let limited =
        PreparedSourceConjugateLinks::new(&thermal, &[binding(&short)], 1, || Ok(())).unwrap();
    assert!(matches!(
        limited.resolve(
            &thermal,
            config(),
            |_, _, t| liquid.at_temperature(t),
            || Ok(())
        ),
        Err(SourceConjugateError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    // A restricted source endpoint is allowed a small h(T)->T(h) round trip.
    // Native link heat must still use the original UNSNAPPED caloric trace.
    let shifted = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        solid.content_hash(),
        fs_material::phase::SolidLiquidPhase::Solid,
        [280.0, 320.0]
            .map(|t| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1e7 + CP * (t - 280.0),
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap();
    for (lower, inside) in [(305.23456789, false), (305.000000001, true)] {
        let restricted = wall_layer(&solid, &shifted, &point, [lower, 320.0]);
        let h = restricted.enthalpy_range_j_kg()[0];
        let native = shifted.temperature_at_specific_enthalpy(h).unwrap();
        assert_eq!(restricted.temperature_at_enthalpy(h).unwrap(), lower);
        assert_ne!(native, lower);
        assert_eq!(native > lower, inside);
        let mut events = thermal
            .cells()
            .values()
            .cloned()
            .map(|mut cell| {
                cell.specific_h_j_kg = if cell.id == IDS[0] {
                    h
                } else {
                    FLUID_CP * (lower + 1.0 - 300.0)
                };
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                }
            })
            .collect::<Vec<_>>();
        events.extend(thermal.links().values().cloned().map(|link| TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link),
        }));
        let state = DepositionSimulation::new(
            events,
            &[
                (MaterialId(0), &shifted as &dyn ThermalLaw),
                (MaterialId(1), fluid_caloric.curve() as &dyn ThermalLaw),
            ],
            DepositionBudget {
                max_cells: 2,
                max_links: 1,
                max_events: 3,
                max_substeps: 200,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let prepared =
            PreparedSourceConjugateLinks::new(&state, &[binding(&restricted)], 1, || Ok(()));
        if inside {
            let resolved = prepared
                .unwrap()
                .resolve(
                    &state,
                    config(),
                    |_, _, t| liquid.at_temperature(t),
                    || Ok(()),
                )
                .unwrap();
            let sample = &resolved.samples()[0];
            assert_eq!(sample.wall_temperature_k.to_bits(), native.to_bits());
            assert_eq!(
                sample.applied_heat_into_fluid_w,
                resolved.conductances_w_k()[&5]
                    * (state.temperature_k(IDS[0]).unwrap() - state.temperature_k(IDS[1]).unwrap())
            );
        } else {
            assert!(
                matches!(
                    prepared,
                    Err(SourceConjugateError::Thermal(
                        DepositionError::MaterialDomain
                    ))
                ),
                "an unsnapped trace outside original source support must refuse"
            );
        }
    }
    let count = std::cell::Cell::new(0usize);
    let interrupted = owner.resolve(
        &thermal,
        config(),
        |_, _, t| {
            count.set(count.get() + 1);
            liquid.at_temperature(t)
        },
        || {
            if count.get() >= 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        interrupted,
        Err(SourceConjugateError::Thermal(DepositionError::Cancelled))
            | Err(SourceConjugateError::Convection(
                SourceCellConvectionError::Convection(ConvectionSourceError::Cancelled)
            ))
    ));
    let accepted = owner
        .resolve(
            &thermal,
            config(),
            |_, _, t| liquid.at_temperature(t),
            || Ok(()),
        )
        .unwrap();
    let mut map = BTreeMap::from([(99, 3.0)]);
    accepted.apply_to(&mut map, || Ok(())).unwrap();
    assert_eq!(map[&99], 3.0);
    let saved_map = map.clone();
    assert!(accepted.apply_to(&mut map, || Ok(())).is_err());
    assert_eq!(map, saved_map);
    assert_eq!(thermal.cells(), &saved);
}
