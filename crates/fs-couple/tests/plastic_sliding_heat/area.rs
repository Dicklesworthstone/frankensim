//! Source-backed hardness recruitment followed by actual wear-induced loss of
//! contact. Cooling does not heal plastic dents or replace removed material.
use super::*;

const CENTRE: usize = 12;
const NEIGHBOUR: usize = 13;
const COLD_CLEARANCE: f64 = 1e-8;
const FACE_AREA: f64 = 25.0 * AREA;
const SLIDE: f64 = 0.1;

// Independent rectangular 1/r integral in asinh form. The source cards give
// E=20 MPa and nu=.25 on both sides; this is the classical combined compliance.
fn compliance(x: f64) -> f64 {
    fn primitive(x: f64, y: f64) -> f64 {
        (if x == 0.0 {
            0.0
        } else {
            x * (y / x.abs()).asinh()
        }) + if y == 0.0 {
            0.0
        } else {
            y * (x / y.abs()).asinh()
        }
    }
    let h = 0.5 * SIDE;
    let modulus = 2e7 / (2.0 * (1.0 - 0.25 * 0.25));
    (primitive(x + h, h) - primitive(x - h, h) - primitive(x + h, -h) + primitive(x - h, -h))
        / (std::f64::consts::PI * modulus)
}

fn two_cell_fixture(cx: &Cx<'_>) -> Fixture {
    let mut data = fixture(cx);
    let original = data.source.wear_source();
    let mut gap = vec![1e-3; 25];
    gap[CENTRE] = 0.0;
    // The cold centre is capped at 1 MPa. Its nonlocal elastic opening leaves
    // the neighbour separated by COLD_CLEARANCE at the accepted approach.
    gap[NEIGHBOUR] = 5e-5 - compliance(SIDE) * 1e6 + COLD_CLEARANCE;
    let normal = SourcePlasticContact::new(
        original.contact().normal_source().clone(),
        original.hardness_sources().clone(),
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: SIDE,
            cell_depth_m: SIDE,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 128,
            complementarity_tolerance_m: 1e-12,
            boundary_clearance_cells: 1,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1e-4,
            maximum_plastic_depth_m: 1e-4,
            half_space_depth_m: [LENGTH; 2],
            maximum_patch_to_depth: [0.25; 2],
        },
        cx,
    )
    .unwrap();
    let wear = SourcePlasticArchardWear::new(
        normal,
        original.coefficient_source().clone(),
        original.coefficient_keys().clone(),
        original.hardness_sources().clone(),
        [5e-5; 2],
    )
    .unwrap();
    data.source =
        SourcePlasticSlidingWear::new(wear, data.source.friction_source().clone()).unwrap();
    data
}

fn full_face_network(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let seed = network(curves, T0);
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.mass_kg = RHO * FACE_AREA * LENGTH;
            cell.geometry.volume_m3 = FACE_AREA * LENGTH;
            cell.geometry.surface_area_m2 = 2.0 * FACE_AREA + 4.0 * (5.0 * SIDE) * LENGTH;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(seed.links()[&LINK]),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        seed.budget(),
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn area_controls() -> PlasticSlidingHeatControls {
    PlasticSlidingHeatControls {
        temperature_scale_j_per_k: RHO * FACE_AREA * LENGTH * CP,
        ..controls()
    }
}

fn both_slide() -> Vec<PlasticSlidingInput> {
    let mut values = motion(SLIDE);
    values[NEIGHBOUR].released_increment_m[0] = SLIDE;
    values
}

// With two loaded cells at the root, symmetric finite heat storage gives
// theta_body=b*theta_interface. H(T) fixes the capped centre, and the exact
// two-cell elastic equation supplies the newly recruited neighbour pressure.
fn recruited_root() -> (f64, f64, [f64; 2]) {
    let capacity = RHO * FACE_AREA * LENGTH * CP;
    let g = CONDUCTIVITY * (2.0 * AREA) / (0.5 * LENGTH);
    let b = g / (capacity + g); // dt=1 s
    let c0 = compliance(0.0);
    let c1 = compliance(SIDE);
    let p0 = 1e6 - COLD_CLEARANCE / c0;
    let slope = 1e6 * 0.01 * b * (1.0 - c1 / c0);
    let qa = AREA * SLIDE * 0.01 * slope;
    let qb = 2.0 * capacity * b - AREA * SLIDE * (0.01 * p0 - 0.4 * slope);
    let supply = AREA * SLIDE * 0.4 * p0;
    let theta = 2.0 * supply / (qb + (qb * qb + 4.0 * qa * supply).sqrt());
    let cap = 1e6 * (1.0 - 0.01 * b * theta);
    let neighbour = (c1 * (1e6 - cap) - COLD_CLEARANCE) / c0;
    assert!(neighbour > 0.0 && neighbour < cap);
    (theta, b * theta, [cap, neighbour])
}

#[test]
fn hardness_recruits_contact_and_final_wear_footprint_controls_the_next_heat_step() {
    with_cx(|cx, _| {
        let data = two_cell_fixture(cx);
        let thermal = full_face_network(&data.caloric);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let state = loaded(&data, cx);
        let cold = state.wear_state().contact_state().contact_state();
        close(cold.pressure_pa()[CENTRE], 1e6, 1e-12, 1e-6);
        assert_eq!(cold.pressure_pa()[NEIGHBOUR], 0.0);
        let contact = PlasticSlidingHeatContact::new_evolving_area(
            &thermal,
            data.source.clone(),
            state.clone(),
            [&data.cards[0], &data.cards[1]],
            [&data.caloric[0], &data.caloric[1]],
            &route,
            limits(),
            cx,
        )
        .unwrap();
        assert_eq!(
            contact.footprint(),
            PlasticSlidingHeatFootprint::EvolvingActiveCells
        );
        let capacity = RHO * FACE_AREA * LENGTH * CP;
        let (theta, bulk_rise, pressure) = recruited_root();
        let original_cells = thermal.cells().clone();
        let first = contact
            .step_candidate(
                &thermal,
                1.0,
                &both_slide(),
                PlasticSlidingHeatBoundary::HeldApproach,
                area_controls(),
                &[0.8 * capacity * bulk_rise; 2],
                T0 + 0.8 * theta,
                |endpoint, t| {
                    assert_eq!(
                        endpoint.contact_patch(LINK).unwrap().covered_area_m2,
                        [2.0 * AREA; 2]
                    );
                    transport(endpoint, t)
                },
                cx,
            )
            .unwrap();
        let report = first.report();
        close(report.heat.interface_temperature_k, T0 + theta, 0.0, 2e-6);
        assert_eq!(report.previous_contact_area_m2, AREA);
        assert_eq!(report.refreshed_contact_area_m2, 2.0 * AREA);
        assert_eq!(report.endpoint_contact_area_m2, 2.0 * AREA);
        assert_eq!(report.conductivity_w_m_k, [CONDUCTIVITY; 2]);
        assert!(report.trials > 1);
        for (local, index) in [CENTRE, NEIGHBOUR].into_iter().enumerate() {
            close(
                report.mechanics.cells()[index].pressure_pa,
                pressure[local],
                2e-6,
                1e-5,
            );
            close(
                first.contact().state().tangential_states()[index].slip_m()[0],
                SLIDE,
                0.0,
                1e-15,
            );
        }
        let pressure_sum = pressure.iter().sum::<f64>();
        let heat = AREA * SLIDE * pressure_sum * (0.4 + 0.01 * theta);
        close(report.frictional_work_j, heat, 2e-6, 1e-12);
        close(
            report.heat.thermal.stored_energy_change_j,
            heat,
            2e-6,
            3e-11,
        );
        for (side, id) in CELLS.into_iter().enumerate() {
            close(
                first.thermal().temperature_k(id).unwrap(),
                T0 + bulk_rise,
                0.0,
                2e-6,
            );
            let volume =
                [1e-6, 2e-6][side] * AREA * SLIDE * pressure_sum / [pressure[0], 2e6][side];
            close(report.removed_volume_m3[side], volume, 2e-6, 1e-23);
            close(
                report.removed_volume_fraction[side],
                volume / (FACE_AREA * LENGTH),
                2e-6,
                1e-15,
            );
            close(
                report.heat.conductances_w_k[side],
                CONDUCTIVITY * 2.0 * AREA / (0.5 * LENGTH),
                1e-12,
                1e-18,
            );
            assert_eq!(
                first.thermal().cells()[&id].mass_kg,
                original_cells[&id].mass_kg
            );
        }
        assert_eq!(first.contact().route().contact_area_m2(), 2.0 * AREA);
        assert!(FrictionHeatRoute::from_link(first.thermal(), LINK, [CONDUCTIVITY; 2]).is_ok());
        assert_eq!(contact.state(), &state);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(
            thermal.contact_patch(LINK).unwrap().covered_area_m2,
            [AREA; 2]
        );

        let (hot, hot_thermal, _) = first.into_parts();
        let saved_hot = hot.state().clone();
        let saved_normal = saved_hot.wear_state().contact_state().contact_state();
        let p = saved_normal.pressure_pa()[NEIGHBOUR];
        assert!(p > pressure[1]);
        // The low-pressure neighbour can wear through its small elastic closure
        // while the centre remains loaded. Its entire interval uses its genuine
        // beginning pressure and reverse actual slip; it is detached only at the
        // endpoint. No lost heat or fictitious zero-pressure source is invented.
        let reverse = 6.0;
        let mut sliding = motion(0.0);
        sliding[NEIGHBOUR].released_increment_m[0] = -reverse;
        let g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
        let b = g / (capacity + g);
        let measure = AREA * p * reverse;
        let theta_second = (2.0 * capacity * b * bulk_rise + 0.4 * measure)
            / (2.0 * capacity * b - 0.01 * measure);
        let bulk_second = bulk_rise + b * (theta_second - bulk_rise);
        let cap_second = 1e6 * (1.0 - 0.01 * bulk_second);
        assert!(saved_normal.pressure_pa()[CENTRE] < cap_second);
        let recession = reverse * p * (1e-6 / cap_second + 2e-6 / 2e6);
        let elastic_closure = (compliance(0.0) - compliance(SIDE).powi(2) / compliance(0.0)) * p;
        assert!(recession > 1.2 * elastic_closure);
        let second = hot
            .step_candidate(
                &hot_thermal,
                2.0,
                &sliding,
                PlasticSlidingHeatBoundary::HeldApproach,
                area_controls(),
                &[0.8 * capacity * (bulk_second - bulk_rise); 2],
                T0 + 0.8 * theta_second,
                |endpoint, t| {
                    assert_eq!(
                        endpoint.contact_patch(LINK).unwrap().covered_area_m2,
                        [AREA; 2]
                    );
                    transport(endpoint, t)
                },
                cx,
            )
            .unwrap();
        let report = second.report();
        assert_eq!(report.previous_contact_area_m2, 2.0 * AREA);
        assert_eq!(report.refreshed_contact_area_m2, 2.0 * AREA);
        assert_eq!(report.endpoint_contact_area_m2, AREA);
        close(
            report.heat.interface_temperature_k,
            T0 + theta_second,
            0.0,
            2e-6,
        );
        close(
            report.mechanics.cells()[NEIGHBOUR].pressure_pa,
            p,
            2e-6,
            1e-5,
        );
        close(
            report.frictional_work_j,
            measure * (0.4 + 0.01 * theta_second),
            2e-6,
            1e-12,
        );
        let next = second.contact().state();
        let normal = next.wear_state().contact_state().contact_state();
        assert_eq!(normal.pressure_pa()[NEIGHBOUR], 0.0);
        assert!(normal.pressure_pa()[CENTRE] > 0.0);
        assert_eq!(normal.plastic_depth_m(), saved_normal.plastic_depth_m());
        assert_eq!(
            normal.cumulative_plastic_work_j(),
            saved_normal.cumulative_plastic_work_j()
        );
        close(
            next.tangential_states()[CENTRE].slip_m()[0],
            SLIDE,
            0.0,
            1e-15,
        );
        close(
            next.tangential_states()[NEIGHBOUR].slip_m()[0],
            SLIDE - reverse,
            0.0,
            1e-14,
        );
        for (side, id) in CELLS.into_iter().enumerate() {
            close(
                second.thermal().temperature_k(id).unwrap(),
                T0 + bulk_second,
                0.0,
                2e-6,
            );
            close(report.heat.conductances_w_k[side], g, 1e-12, 1e-18);
            let increment = [1e-6, 2e-6][side] * AREA * reverse * p / [cap_second, 2e6][side];
            close(
                next.wear_state().volumes_m3()[side] - saved_hot.wear_state().volumes_m3()[side],
                increment,
                2e-6,
                1e-23,
            );
            assert_eq!(
                next.wear_state().recession_m()[side][CENTRE],
                saved_hot.wear_state().recession_m()[side][CENTRE]
            );
            assert!(
                next.wear_state().recession_m()[side][NEIGHBOUR]
                    > saved_hot.wear_state().recession_m()[side][NEIGHBOUR]
            );
            assert_eq!(
                second.thermal().cells()[&id].geometry,
                original_cells[&id].geometry
            );
        }
        assert_eq!(second.contact().route().contact_area_m2(), AREA);
        assert!(FrictionHeatRoute::from_link(second.thermal(), LINK, [CONDUCTIVITY; 2]).is_ok());
        assert_eq!(hot.state(), &saved_hot);
        assert_eq!(
            hot_thermal.contact_patch(LINK).unwrap().covered_area_m2,
            [2.0 * AREA; 2]
        );
        assert_eq!(second.contact().time_s(), 2.0);
    });
}

#[test]
fn evolving_area_is_explicit_and_refusal_retains_all_original_histories() {
    with_cx(|cx, _| {
        let data = two_cell_fixture(cx);
        let state = loaded(&data, cx);
        let small = network(&data.caloric, T0);
        let small_route = FrictionHeatRoute::from_link(&small, LINK, [CONDUCTIVITY; 2]).unwrap();
        // The present A-sized patch fits, but its entire future 25*A support
        // does not. An evolving owner must reject that incomplete registration.
        assert!(matches!(
            PlasticSlidingHeatContact::new_evolving_area(
                &small,
                data.source.clone(),
                state.clone(),
                [&data.cards[0], &data.cards[1]],
                [&data.caloric[0], &data.caloric[1]],
                &small_route,
                limits(),
                cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let thermal = full_face_network(&data.caloric);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let fixed = bind(&data, &thermal, &route, state.clone(), limits(), cx);
        let (theta, bulk, _) = recruited_root();
        let guess = [RHO * FACE_AREA * LENGTH * CP * bulk; 2];
        assert!(matches!(
            fixed.step_candidate(
                &thermal,
                1.0,
                &both_slide(),
                PlasticSlidingHeatBoundary::HeldApproach,
                area_controls(),
                &guess,
                T0 + theta,
                transport,
                cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let evolving = PlasticSlidingHeatContact::new_evolving_area(
            &thermal,
            data.source.clone(),
            state.clone(),
            [&data.cards[0], &data.cards[1]],
            [&data.caloric[0], &data.caloric[1]],
            &route,
            limits(),
            cx,
        )
        .unwrap();
        let reached_covered_source = Cell::new(false);
        let cells = thermal.cells().clone();
        assert!(matches!(
            evolving.step_candidate(
                &thermal,
                1.0,
                &both_slide(),
                PlasticSlidingHeatBoundary::HeldApproach,
                area_controls(),
                &guess,
                T0 + theta,
                |endpoint, _| {
                    assert_eq!(
                        endpoint.contact_patch(LINK).unwrap().covered_area_m2,
                        [2.0 * AREA; 2]
                    );
                    reached_covered_source.set(true);
                    Err::<PlasticSlidingHeatSource<()>, _>(DepositionError::MaterialDomain)
                },
                cx,
            ),
            Err(PlasticSlidingHeatError::Source(
                DepositionError::MaterialDomain
            ))
        ));
        assert!(reached_covered_source.get());
        assert_eq!(evolving.state(), &state);
        assert_eq!(evolving.source(), &data.source);
        assert_eq!(fixed.state(), &state);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(
            thermal.contact_patch(LINK).unwrap().covered_area_m2,
            [AREA; 2]
        );
        assert_eq!(evolving.route().contact_area_m2(), AREA);
    });
}
