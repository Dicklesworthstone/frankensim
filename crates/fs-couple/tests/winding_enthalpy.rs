#[path = "winding_enthalpy/common.rs"]
mod common;
#[path = "../../fs-em/tests/winding_ohmic/common.rs"]
mod field;
use common::*;
use fs_couple::circuit_enthalpy::CircuitEnthalpyError;
use fs_couple::circuit_enthalpy::winding::WindingThermalNetwork;
use fs_em::winding::ohmic::OhmicWindingCircuit;
use fs_matdb::PropertyKey;
use fs_material::{conductor::ElectricalPropertyBasis, state_point::THERMAL_CONDUCTIVITY_DIMS};
use fs_thermal::enthalpy::CellId;

#[test]
fn field_material_current_heat_and_native_thermal_storage_share_one_endpoint() {
    let mut outcomes = Vec::new();
    for scale in [1.0, 2.0] {
        let alpha = 0.002;
        let cards = cards(scale, alpha);
        let curves = curves(&cards);
        let thermal = initial(&curves, 0.25, None);
        let basis = ElectricalPropertyBasis::Resistivity;
        let (problem, profiles, boundary) = field::problem(&cards, basis, 0.5);
        let magnetic = field::winding(&problem, profiles, &boundary);
        let circuit = OhmicWindingCircuit::new(
            &magnetic,
            field::sources(&magnetic, basis),
            0.0,
            &[10.0, 0.0],
            &[0.2],
            &[14.0],
            &mut || false,
        )
        .unwrap();
        let context = field::point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let routes = routes(&circuit);
        let layout = network(&cells, &routes);
        let dt = 0.1;
        let next = layout
            .step_candidate(&thermal, dt, &circuit, &[14.0], config(), &[0.0; 3], || {
                false
            })
            .unwrap();
        // at_query_point pins the original claims; a changed receipt policy
        // must not masquerade as an initial physical-temperature mismatch.
        for (old, new) in circuit.conductors().iter().zip(next.circuit().conductors()) {
            let old_receipt = &old.conductor.source_property().answer().receipt;
            let new_receipt = &new.conductor.source_property().answer().receipt;
            assert_eq!(old_receipt.policy, "single-claim-only");
            assert_eq!(new_receipt.policy, "pinned-claim");
            assert_eq!(old_receipt.selected, new_receipt.selected);
        }
        let t0 = next.thermal().temperature_k(CellId(30)).unwrap();
        let t1 = next.thermal().temperature_k(CellId(10)).unwrap();
        let final_r = next.circuit().resistance_ohm();
        near(final_r[0], 0.8 * scale * (1.0 + alpha * (t0 - 300.0)));
        near(final_r[1], 1.2 * (1.0 + alpha * (t1 - 300.0)));
        assert!(t0 > 300.0 && t1 > 300.0);
        // Independent coupled RL midpoint identity using actual endpoint R(T).
        let l = magnetic.basis().inductance_h();
        let expected = field::two_by_two(
            [
                l[0] + dt * (1.0 + final_r[0]) / 2.0,
                l[1],
                l[2],
                l[3] + dt * (3.0 + final_r[1]) / 2.0,
            ],
            [
                10.0 * (l[0] - dt * (1.0 + final_r[0]) / 2.0) + dt * 14.0,
                10.0 * l[2],
            ],
        );
        let currents = next.circuit().currents_a().unwrap();
        for i in 0..2 {
            near(currents[i], expected[i]);
        }
        let report = next.report();
        let electrical = &report.electrical;
        let mut coil = [0.0; 2];
        for loss in &electrical.cells {
            coil[loss.winding] += loss.heat_j;
            let owner = loss.winding;
            assert_eq!(
                loss.source.conductor.material().card_identity(),
                cards[owner].content_hash()
            );
            let material_t = loss
                .source
                .conductor
                .material()
                .query_point()
                .iter()
                .find(|(name, _)| name == "T")
                .unwrap()
                .1;
            near(material_t, if owner == 0 { t0 } else { t1 });
        }
        let branch = |index| {
            electrical
                .graph_resistors
                .iter()
                .find(|r| r.branch == index)
                .unwrap()
                .heat_j
        };
        let g = report.source.network.links[0].conductance_w_k;
        near(g, 2.0);
        // The actual caloric card has Cp=4 J/(kg K), density=1 kg/m3.
        // .25 m3 conductor cells therefore each store exactly 1 J/K.
        near(t0 - 300.0, coil[0] + branch(1) + dt * g * (t1 - t0));
        near(t1 - 300.0, coil[1] + branch(5) + dt * g * (t0 - t1));
        near(
            0.4 * (next.thermal().temperature_k(CellId(50)).unwrap() - 300.0),
            branch(3) + dt * 0.7,
        );
        near(report.deposited_joule_j, electrical.heat_j);
        near(report.non_electrical_heat_j, dt * 0.7);
        near(
            electrical.stored_after_j - electrical.stored_before_j
                + report.thermal.stored_energy_change_j,
            electrical.native.supplied + dt * 0.7,
        );
        assert!(report.maximum_cell_residual_j <= config().cell_energy_tolerance_j);
        assert!(report.total_energy_residual_j.abs() <= config().total_energy_tolerance_j);
        let frozen = circuit
            .step_to_candidate(
                dt,
                circuit.conductors().to_vec(),
                &[14.0],
                field::controls(),
                &mut || false,
            )
            .unwrap();
        assert!(currents[0] < frozen.circuit().currents_a().unwrap()[0]);
        assert!(report.circuit_trials > 1);
        assert_eq!(next.circuit().time_s(), dt);
        assert_eq!(next.thermal().time_s(), dt);
        assert_eq!(circuit.time_s(), 0.0);
        near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
        let continued = layout
            .step_candidate(
                next.thermal(),
                0.2,
                next.circuit(),
                &[14.0],
                config(),
                &[0.0; 3],
                || false,
            )
            .unwrap();
        near(
            continued.report().electrical.stored_before_j,
            electrical.stored_after_j,
        );
        assert_eq!(continued.circuit().time_s(), 0.2);
        assert!(continued.thermal().temperature_k(CellId(30)).unwrap() > t0);
        assert_eq!(next.circuit().time_s(), dt);
        outcomes.push((currents[0], t0, coil[0]));
    }
    assert!(outcomes[1].0 < outcomes[0].0);
    assert!(outcomes[1].1 > outcomes[0].1);
    assert!(outcomes[1].2 > outcomes[0].2);
}

#[test]
fn incomplete_graph_heat_wrong_geometry_and_caloric_substitution_refuse() {
    let cards = cards(1.0, 0.002);
    let curves = curves(&cards);
    let thermal = initial(&curves, 0.25, None);
    let basis = ElectricalPropertyBasis::Resistivity;
    let (problem, profiles, boundary) = field::problem(&cards, basis, 0.5);
    let magnetic = field::winding(&problem, profiles, &boundary);
    let circuit = OhmicWindingCircuit::new(
        &magnetic,
        field::sources(&magnetic, basis),
        0.0,
        &[10.0, 0.0],
        &[0.2],
        &[14.0],
        &mut || false,
    )
    .unwrap();
    let context = field::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let routes = routes(&circuit);
    let layout = network(&cells, &routes);
    // Receipt-policy changes are allowed; an actually inconsistent accepted
    // physical temperature still refuses before any nonlinear trial.
    let warmer_sources = circuit
        .conductors()
        .iter()
        .zip(circuit.conductor_windings())
        .map(|(source, &winding)| {
            let snapshot = source
                .conductor
                .material()
                .at_query_point(&cards[winding], &field::point(301.0))
                .unwrap();
            fs_em::winding::ohmic::WindingConductor {
                triangle: source.triangle,
                conductor:
                    fs_material::conductor::ResolvedConductorState::try_from_material_with_basis(
                        snapshot, basis,
                    )
                    .unwrap(),
            }
        })
        .collect();
    let warmer = OhmicWindingCircuit::new(
        &magnetic,
        warmer_sources,
        0.0,
        &[10.0, 0.0],
        &[0.2],
        &[14.0],
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        layout.step_candidate(&thermal, 0.1, &warmer, &[14.0], config(), &[0.0; 3], || {
            false
        }),
        Err(CircuitEnthalpyError::Source(
            fs_couple::circuit_enthalpy::winding::WindingThermalError::Invalid(
                "accepted winding source snapshots must match the actual thermal endpoint"
            )
        ))
    ));
    let missing = WindingThermalNetwork {
        graph_resistors: &GRAPH[..2],
        ..layout
    };
    assert!(
        missing
            .step_candidate(
                &thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    let wrong_routes = [GRAPH[0], GRAPH[0], GRAPH[2]];
    let duplicate = WindingThermalNetwork {
        graph_resistors: &wrong_routes,
        ..network(&cells, &routes)
    };
    assert!(
        duplicate
            .step_candidate(
                &thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    let mut lost_heat = routes.clone();
    if let fs_couple::circuit_enthalpy::winding::MagneticRegion::Cell { triangle, .. } =
        lost_heat[0]
    {
        lost_heat[0] = fs_couple::circuit_enthalpy::winding::MagneticRegion::Fixed { triangle };
    }
    assert!(
        network(&cells, &lost_heat)
            .step_candidate(
                &thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    let wrong_thermal = initial(&curves, 0.2, None);
    assert!(
        network(&cells, &routes)
            .step_candidate(
                &wrong_thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    let cloned_curves = curves.clone();
    let cloned_cells = common::cells(&cards, &cloned_curves, &context, &key);
    assert!(
        network(&cloned_cells, &routes)
            .step_candidate(
                &thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    let foreign_card = common::cards(2.0, 0.002);
    let foreign_cells = common::cells(&foreign_card, &curves, &context, &key);
    assert!(
        network(&foreign_cells, &routes)
            .step_candidate(
                &thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    assert_eq!(circuit.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn native_event_budget_and_late_cancellation_preserve_both_accepted_histories() {
    let cards = cards(1.0, 0.002);
    let curves = curves(&cards);
    let thermal = initial(&curves, 0.25, None);
    let basis = ElectricalPropertyBasis::Resistivity;
    let (problem, profiles, boundary) = field::problem(&cards, basis, 0.5);
    let magnetic = field::winding(&problem, profiles, &boundary);
    let circuit = OhmicWindingCircuit::new(
        &magnetic,
        field::sources(&magnetic, basis),
        0.0,
        &[10.0, 0.0],
        &[0.2],
        &[14.0],
        &mut || false,
    )
    .unwrap();
    let context = field::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let routes = routes(&circuit);
    let layout = network(&cells, &routes);
    let event = initial(&curves, 0.25, Some(0.05));
    assert!(matches!(
        layout.step_candidate(&event, 0.1, &circuit, &[14.0], config(), &[0.0; 3], || {
            false
        }),
        Err(CircuitEnthalpyError::Thermal(_))
    ));
    let mut budget = config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_candidate(&thermal, 0.1, &circuit, &[14.0], budget, &[0.0; 3], || {
            false
        }),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut calls = 0;
    let expected = layout
        .step_candidate(
            &thermal,
            0.1,
            &circuit,
            &[14.0],
            config(),
            &[0.0; 3],
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        layout.step_candidate(
            &thermal,
            0.1,
            &circuit,
            &[14.0],
            config(),
            &[0.0; 3],
            || {
                seen += 1;
                seen == calls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let replay = layout
        .step_candidate(
            &thermal,
            0.1,
            &circuit,
            &[14.0],
            config(),
            &[0.0; 3],
            || false,
        )
        .unwrap();
    assert_eq!(expected.circuit().state(), replay.circuit().state());
    for id in [CellId(10), CellId(30), CellId(50)] {
        assert_eq!(
            expected.thermal().cells()[&id].specific_h_j_kg,
            replay.thermal().cells()[&id].specific_h_j_kg
        );
    }
    assert_eq!(
        expected.report().electrical.cells,
        replay.report().electrical.cells
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.time_s(), 0.0);
}

#[test]
fn every_passive_magnetic_region_is_fixed_or_refreshed_from_its_actual_thermal_cell() {
    use fs_couple::circuit_enthalpy::material::MaterialCell;
    use fs_couple::circuit_enthalpy::winding::{MagneticRegion, WindingThermalError};
    use fs_em::planar::PlanarMagnetostatics;
    use fs_feec::weighted::WeightedAssemblyLimits;
    use fs_material::thermal::source::ThermalMaterialPoint;
    use fs_thermal::enthalpy::MaterialId;
    for source_kind in ["constant", "temperature", "pressure"] {
        let temperature_dependent = source_kind == "temperature";
        let pressure_dependent = source_kind == "pressure";
        let cards = cards(1.0, 0.002);
        let curves = curves(&cards);
        let basis = ElectricalPropertyBasis::Resistivity;
        let core = if temperature_dependent {
            field::magnetic_curve_card(4e-5, 2e-4, basis)
        } else if pressure_dependent {
            field::pressure_magnetic_card(4e-5, 2e-4, basis)
        } else {
            field::card(4e-5, 0.001, 2e-4, basis)
        };
        let core_context = if pressure_dependent {
            field::pressure_point(300.0, 1e5)
        } else {
            field::point(300.0)
        };
        let core_curves = if pressure_dependent {
            use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
            use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
            use fs_qty::Density;
            [core.clone(), core.clone()].each_ref().map(|card| {
                ResolvedHeatCapacityCurve::try_new(
                    card,
                    HeatCapacitySourceQuery {
                        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
                        density: &PropertyKey::new("density", Density::DIMS),
                        temperature_axis: "T",
                        lower: &field::pressure_point(200.0, 1e5),
                        upper: &field::pressure_point(2000.0, 1e5),
                        selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
                    },
                    SolidLiquidPhase::Solid,
                    0.0,
                    1e-6,
                ).unwrap()
            })
        } else {
            common::curves(&[core.clone(), core.clone()])
        };
        let (original, mut profiles, boundary) = field::problem(&cards, basis, 0.5);
        let core_triangles = [0, 1, 14, 15]; // equal outgoing/return area, so current still closes.
        let mut materials = original.materials().to_vec();
        for &triangle in &core_triangles {
            profiles[0][triangle] = 0.0;
            materials[triangle] = field::bundle_at(&core, &core_context, basis);
        }
        let problem = PlanarMagnetostatics::new(
            original.mesh().clone(),
            materials,
            WeightedAssemblyLimits {
                max_cells: 1000,
                max_dofs: 1000,
                max_triplets: 10000,
            },
            &mut || false,
        )
        .unwrap();
        let magnetic = field::winding(&problem, profiles, &boundary);
        let sources = field::sources(&magnetic, basis)
            .into_iter()
            .filter(|s| {
                magnetic
                    .basis()
                    .profiles_per_m2()
                    .iter()
                    .any(|p| p[s.triangle] != 0.0)
            })
            .collect();
        let circuit = OhmicWindingCircuit::new(
            &magnetic,
            sources,
            0.0,
            &[10.0, 0.0],
            &[0.2],
            &[14.0],
            &mut || false,
        )
        .unwrap();
        let context = field::point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let fixed_thermal = initial(&curves, 0.234375, None);
        let mut regions = routes(&circuit);
        assert!(
            network(&cells, &regions)
                .step_candidate(
                    &fixed_thermal,
                    0.1,
                    &circuit,
                    &[14.0],
                    config(),
                    &[0.0; 3],
                    || false
                )
                .is_err()
        );
        regions.extend(core_triangles.map(|triangle| MagneticRegion::Fixed { triangle }));
        let fixed = network(&cells, &regions)
            .step_candidate(
                &fixed_thermal,
                0.1,
                &circuit,
                &[14.0],
                config(),
                &[0.0; 3],
                || false,
            )
            .unwrap();
        for receipt in &fixed.report().source.magnetic {
            if core_triangles.contains(&receipt.triangle) {
                assert_eq!(receipt.cell, None);
                assert_eq!(receipt.material, problem.materials()[receipt.triangle]);
            }
        }
        let mut heated_regions = routes(&circuit);
        heated_regions.extend(core_triangles.map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CellId(40),
        }));
        let mut heated_cells = cells;
        heated_cells.push(MaterialCell {
            cell: CellId(40),
            material_id: MaterialId(3),
            curve: core_curves[0].curve(),
            point: ThermalMaterialPoint {
                material: &core,
                context: &core_context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        });
        let heated_thermal = initial_with_core(&curves, 0.234375, None, Some(&core_curves[0]));
        let result = network(&heated_cells, &heated_regions).step_candidate(
            &heated_thermal,
            0.1,
            &circuit,
            &[14.0],
            config(),
            &[0.0; 4],
            || false,
        );
        if temperature_dependent {
            assert!(matches!(
                result,
                Err(CircuitEnthalpyError::Source(WindingThermalError::Invalid(
                    "heated magnetic region changes prepared storage; new extraction and magnetic parameter work required"
                )))
            ));
        } else {
            let next = result.unwrap();
            near(next.thermal().temperature_k(CellId(40)).unwrap(), 301.6);
            near(next.report().non_electrical_heat_j, 0.17);
            for receipt in &next.report().source.magnetic {
                if core_triangles.contains(&receipt.triangle) {
                    assert_eq!(receipt.cell, Some(CellId(40)));
                    let t = receipt
                        .material
                        .query_point()
                        .iter()
                        .find(|(name, _)| name == "T")
                        .unwrap()
                        .1;
                    near(t, 301.6);
                    assert_eq!(receipt.material.card_identity(), core.content_hash());
                    if pressure_dependent {
                        let permeability = receipt.material.property("magnetic_permeability").unwrap();
                        assert_eq!(permeability.source_axes()[0].0, "pressure");
                        near(permeability.value_si(), 2e-4);
                        let pressure = receipt.material.query_point().iter()
                            .find(|(name, _)| name == "pressure").unwrap().1;
                        near(pressure, 1e5);
                    }
                }
            }
        }
        assert_eq!(heated_thermal.time_s(), 0.0);
        assert_eq!(circuit.time_s(), 0.0);
    }
}
