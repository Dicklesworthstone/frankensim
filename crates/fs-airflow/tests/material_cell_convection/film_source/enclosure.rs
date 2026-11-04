//! Full gas-gap conductance on the original two-wall thermal owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::enclosure::{
    SourceEnclosedGasLayer, evaluate_source_enclosed_gas_layers,
    observer::SourceEnclosedGasLayerObserver,
};
use fs_material::MaterialError;
use fs_material::thermochemical_gas::ThermochemicalGasModel;

#[path = "../../finite_gas_convection/source.rs"]
mod gas_source;

#[path = "enclosure/inclined.rs"]
mod inclined;
#[path = "enclosure/vertical.rs"]
mod vertical;

const GAP: f64 = 0.018;
const SPANS: [f64; 2] = [0.2, 0.3];
const GAP_AREA: f64 = SPANS[0] * SPANS[1];
const PRESSURE: f64 = 100000.0;

fn binding() -> SourceEnclosedGasLayer {
    SourceEnclosedGasLayer {
        link: 5,
        lower_cell: IDS[0],
        gap: Length::new(GAP),
        plate_spans: SPANS.map(Length::new),
        gravity_m_s2: 9.81,
        inclination_degrees: 0.0,
        maximum_relative_density_departure: 0.1,
    }
}

fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}

// Original immutable caloric chart, with unequal masses so the mean gas
// property temperature changes during otherwise closed reciprocal exchange.
fn pair<'a>(
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    temperatures: [f64; 2],
    reversed: bool,
    authored_g: f64,
) -> DepositionSimulation<'a> {
    let mut events = IDS
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                body_id: i as u64,
                segment_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: MASS * (i + 1) as f64,
                specific_h_j_kg: CP * (temperatures[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 * GAP],
                    end_m: [SPANS[0], 0.0, i as f64 * GAP],
                    volume_m3: 1e-6 * (i + 1) as f64,
                    surface_area_m2: 2.0 * GAP_AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        })
        .collect::<Vec<_>>();
    let [a, b] = if reversed { [IDS[1], IDS[0]] } else { IDS };
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 5,
            a,
            b,
            area_m2: GAP_AREA,
            conductance_w_k: authored_g,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), chart as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}

// Independent dimensional algebra and source gas-layer expressions.
fn expected(model: &ThermochemicalGasModel, lower: f64, upper: f64) -> (f64, f64, f64) {
    expected_at_inclination(model, lower, upper, 0.0)
}

fn expected_at_inclination(
    model: &ThermochemicalGasModel,
    lower: f64,
    upper: f64,
    inclination_degrees: f64,
) -> (f64, f64, f64) {
    let gas = model
        .evaluate(&gas_source::point(0.5 * lower + 0.5 * upper, PRESSURE))
        .unwrap();
    let s = gas.state();
    let beta = gas
        .thermal_transport_state()
        .unwrap()
        .volumetric_expansion_per_k()
        .unwrap();
    let ra = 9.81 * beta * (lower - upper) * GAP.powi(3) * s.density.powi(2) * s.specific_heat_cp
        / (s.dynamic_viscosity * s.thermal_conductivity);
    let angle = inclination_degrees.to_radians();
    let normal_ra = ra * angle.cos();
    let nu = if inclination_degrees == 90.0 {
        vertical::expected_nusselt(ra.abs(), SPANS[1] / GAP)
    } else if normal_ra <= 1708.0 {
        1.0
    } else {
        1.0 + 1.44
            * (1.0 - 1708.0 / normal_ra)
            * (1.0 - 1708.0 * (1.8 * angle).sin().powf(1.6) / normal_ra)
            + ((normal_ra / 5830.0).cbrt() - 1.0).max(0.0)
    };
    (ra, nu, GAP_AREA * nu * s.thermal_conductivity / GAP)
}

#[test]
fn horizontal_layer_uses_physical_lower_wall_and_retains_finite_conduction_at_rest() {
    let chart = gas_source::body_curve();
    // Pure manufactured A has Pr about 0.7195; this is a source-pack test of
    // the air-policy API, not an assertion of experimental fluid identity.
    let model = gas_source::model(1.0);
    for reversed in [false, true] {
        for temperatures in [[410.0, 390.0], [390.0, 410.0], [400.0, 400.0]] {
            let thermal = pair(&chart, temperatures, reversed, 0.0);
            let saved = thermal.cells().clone();
            let mut queries = Vec::new();
            let result = evaluate_source_enclosed_gas_layers(
                &thermal,
                &[binding()],
                1,
                |_, _, t| {
                    queries.push(t.value());
                    model.evaluate(&gas_source::point(t.value(), PRESSURE))
                },
                || false,
            )
            .unwrap();
            let sample = &result.samples()[0];
            let (ra, nu, g) = expected(&model, temperatures[0], temperatures[1]);
            assert_eq!(queries, vec![400.0]);
            assert_eq!((sample.lower_cell, sample.upper_cell), (IDS[0], IDS[1]));
            near(sample.evaluation.nusselt.groups()["Ra"], ra, 1e-9);
            near(sample.evaluation.nusselt.evidence().value, nu, 1e-12);
            near(sample.conductance_w_k, g, 1e-12);
            near(
                sample.heat_into_upper_w,
                g * (temperatures[0] - temperatures[1]),
                1e-11,
            );
            assert_eq!(sample.material.source_packs(), model.source_packs());
            assert_eq!(sample.material.state().pressure, PRESSURE);
            assert_eq!(sample.material.state().temperature, 400.0);
            if temperatures[0] <= temperatures[1] {
                assert_eq!(sample.evaluation.nusselt.evidence().value, 1.0);
                near(sample.conductance_w_k, GAP_AREA * 0.026 / GAP, 1e-14);
            }
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.links()[&5].conductance_w_k, 0.0);
        }
    }
}

#[test]
fn implicit_layer_resolves_actual_endpoint_mean_and_closes_both_original_wall_energies() {
    let chart = gas_source::body_curve();
    let model = gas_source::model_transport(1.0, true);
    for (initial, inclination_degrees) in [
        ([410.0, 390.0], 0.0),
        ([390.0, 410.0], 0.0),
        ([410.0, 390.0], 45.0),
        ([410.0, 390.0], 90.0),
        ([390.0, 410.0], 90.0),
    ] {
        for reversed in [false, true] {
            let thermal = pair(&chart, initial, reversed, 0.0);
            let saved = thermal.cells().clone();
            let capacities = [MASS * CP, 2.0 * MASS * CP];
            let direction = (initial[0] - initial[1]).signum();
            let dt = 0.5;
            let mut lo = 0.0;
            let mut hi =
                (initial[0] - initial[1]).abs() / (1.0 / capacities[0] + 1.0 / capacities[1]);
            for _ in 0..80 {
                let transfer = 0.5 * lo + 0.5 * hi;
                let lower = initial[0] - direction * transfer / capacities[0];
                let upper = initial[1] + direction * transfer / capacities[1];
                if transfer
                    - dt * expected_at_inclination(&model, lower, upper, inclination_degrees).2
                        * (lower - upper).abs()
                    > 0.0
                {
                    hi = transfer;
                } else {
                    lo = transfer;
                }
            }
            let transferred = direction * (0.5 * lo + 0.5 * hi);
            let endpoint = [
                initial[0] - transferred / capacities[0],
                initial[1] + transferred / capacities[1],
            ];
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let increments = step
                .cell_ids()
                .iter()
                .map(|id| {
                    if *id == IDS[0] {
                        -transferred
                    } else {
                        transferred
                    }
                })
                .collect::<Vec<_>>();
            let trial = step
                .evaluate(
                    &increments,
                    |state| {
                        let result = evaluate_source_enclosed_gas_layers(
                            state,
                            &[SourceEnclosedGasLayer {
                                inclination_degrees,
                                ..binding()
                            }],
                            1,
                            |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
                            || false,
                        )?;
                        let mut conductances_w_k = BTreeMap::new();
                        result.apply_to(&mut conductances_w_k, || false)?;
                        Ok::<_, SourceCellConvectionError<MaterialError>>(StepCoefficients {
                            sources_w: BTreeMap::new(),
                            conductances_w_k,
                            payload: result,
                        })
                    },
                    || Ok(()),
                )
                .unwrap();
            let (next, result, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
            let sample = &result.samples()[0];
            assert_eq!(sample.inclination_degrees, inclination_degrees);
            near(
                sample.property_temperature_k,
                0.5 * endpoint[0] + 0.5 * endpoint[1],
                1e-12,
            );
            assert!((sample.property_temperature_k - 400.0).abs() > 0.01);
            assert_eq!(
                sample.material.state().temperature,
                sample.property_temperature_k
            );
            assert_eq!(sample.material.source_packs(), model.source_packs());
            near(sample.heat_into_upper_w * dt, transferred, 1e-8);
            for (id, expected_t) in IDS.into_iter().zip(endpoint) {
                near(next.temperature_k(id).unwrap(), expected_t, 1e-10);
            }
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            near(audit.residual_j, 0.0, 1e-8);
            assert_eq!(thermal.cells(), &saved);
        }
    }
}

#[test]
fn native_layer_conductance_tightens_original_cfl_and_retains_endpoint_sources() {
    let chart = gas_source::body_curve();
    let model = gas_source::model_transport(1.0, true);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    let bindings = [binding()];
    let observer = SourceEnclosedGasLayerObserver::new(
        &bindings,
        1,
        |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
            model.evaluate(&gas_source::point(t.value(), PRESSURE))
        },
    );
    let (next, accepted, audit) = thermal
        .advance_observed(&observer, 20.0, controls(20.0), || Ok(()))
        .unwrap();
    assert!(
        audit.substeps > 1,
        "the generated total gap conductance must tighten the original timestep"
    );
    near(next.energy_j(), thermal.energy_j(), 1e-8);
    near(audit.external_heat_j, 0.0, 1e-12);
    for id in IDS {
        assert!((390.0..=410.0).contains(&next.temperature_k(id).unwrap()));
    }
    let sample = &accepted.sample().unwrap().samples()[0];
    let applied = &accepted.applied_sample().unwrap().samples()[0];
    near(
        sample.property_temperature_k,
        0.5 * next.temperature_k(IDS[0]).unwrap() + 0.5 * next.temperature_k(IDS[1]).unwrap(),
        1e-12,
    );
    assert_ne!(
        sample.property_temperature_k,
        applied.property_temperature_k
    );
    assert_eq!(
        sample.material.state().temperature,
        sample.property_temperature_k
    );
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());
}

#[test]
fn enclosed_geometry_source_policy_and_competing_conductance_refuse_atomically() {
    let chart = gas_source::body_curve();
    let model = gas_source::model(1.0);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    for b in [
        SourceEnclosedGasLayer {
            gap: Length::new(0.021),
            ..binding()
        },
        SourceEnclosedGasLayer {
            plate_spans: [Length::new(0.2); 2],
            ..binding()
        },
        SourceEnclosedGasLayer {
            lower_cell: CellId(999),
            ..binding()
        },
        SourceEnclosedGasLayer {
            maximum_relative_density_departure: 0.01,
            ..binding()
        },
    ] {
        assert!(
            evaluate_source_enclosed_gas_layers(
                &thermal,
                &[b],
                1,
                |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
                || false
            )
            .is_err()
        );
    }
    let full_patch = thermal.contact_patch(5).unwrap();
    for covered_area_m2 in [
        [GAP_AREA / 2.0; 2],
        [GAP_AREA, GAP_AREA / 2.0],
        [GAP_AREA / 2.0, GAP_AREA],
    ] {
        let partial = thermal
            .replace_contact_patches(
                &[ContactPatchUpdate {
                    id: 5,
                    expected: Some(full_patch),
                    replacement: Some(ContactPatch {
                        covered_area_m2,
                        ..full_patch
                    }),
                }],
                || Ok(()),
            )
            .unwrap();
        assert!(
            evaluate_source_enclosed_gas_layers(
                &partial,
                &[binding()],
                1,
                |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
                || false
            )
            .is_err()
        );
    }
    let with_fourier = pair(&chart, [410.0, 390.0], false, 0.104);
    assert!(
        evaluate_source_enclosed_gas_layers(
            &with_fourier,
            &[binding()],
            1,
            |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
            || false
        )
        .is_err()
    );
    let foreign_pr = gas_source::model(0.75);
    assert!(
        evaluate_source_enclosed_gas_layers(
            &thermal,
            &[binding()],
            1,
            |_, _, t| foreign_pr.evaluate(&gas_source::point(t.value(), PRESSURE)),
            || false
        )
        .is_err()
    );
    assert!(
        evaluate_source_enclosed_gas_layers(
            &thermal,
            &[binding()],
            1,
            |_, _, _| model.evaluate(&gas_source::point(405.0, PRESSURE)),
            || false
        )
        .is_err()
    );
    let result = evaluate_source_enclosed_gas_layers(
        &thermal,
        &[binding()],
        1,
        |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
        || false,
    )
    .unwrap();
    let mut destination = BTreeMap::from([(5, 2.0), (99, 3.0)]);
    let saved = destination.clone();
    assert!(result.apply_to(&mut destination, || false).is_err());
    assert_eq!(destination, saved);
    destination.remove(&5);
    let saved = destination.clone();
    assert!(result.apply_to(&mut destination, || true).is_err());
    assert_eq!(destination, saved);
}

#[test]
fn late_actual_source_refusal_and_final_checkpoint_leave_both_walls_and_observer_unchanged() {
    // The original gas pack ends at 1000 K. Its starting mean is admitted;
    // unequal wall capacities move the candidate mean outside that REAL domain.
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "solid".into(),
            process: "enclosure source-refusal fixture".into(),
            revision: 0,
        },
        ClaimSet::new(),
        vec![],
    )
    .unwrap();
    let wide_chart = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (990.0, 690000.0, 0.0),
            (1005.0, 705000.0, 0.0),
            (1005.0, 765000.0, 1.0),
            (1010.0, 770000.0, 1.0),
        ]
        .map(|(t, h, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: 1000.0,
        })
        .to_vec(),
    )
    .unwrap();
    let model = gas_source::model(1.0);
    let thermal = pair(&wide_chart, [999.6, 1000.2], false, 0.0);
    let bindings = [binding()];
    let observer = SourceEnclosedGasLayerObserver::new(
        &bindings,
        1,
        |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
            model.evaluate(&gas_source::point(t.value(), PRESSURE))
        },
    );
    let saved = thermal.cells().clone();
    assert!(
        matches!(thermal.advance_observed(&observer, 20.0, controls(20.0), || Ok(())),
        Err(SourceCellConvectionError::Source { property_temperature_k, .. }) if property_temperature_k > 1000.0)
    );
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());

    let chart = gas_source::body_curve();
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    let saved = thermal.cells().clone();
    let mut polls = 0;
    thermal
        .advance_observed(&observer, 0.2, controls(0.1), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 0.2, controls(0.1), || {
            current += 1;
            if current == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());
}
