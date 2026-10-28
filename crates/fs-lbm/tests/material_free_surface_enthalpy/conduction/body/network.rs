//! G1/G4: finite wall cells spread heat on the original liquid timeline.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceBodyLink;

fn cell<'a>(
    id: u64,
    card: &'a MaterialCard,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    temperature_k: f64,
) -> FreeSurfaceThermalBody<'a> {
    let end_x = -0.001 * (id as f64 - 42.0);
    FreeSurfaceThermalBody::try_new(
        id,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        chart,
        CellGeometry {
            start_m: [end_x - 0.001, 0.0035, 0.0035],
            end_m: [end_x, 0.0035, 0.0035],
            volume_m3: 1e-9,
            surface_area_m2: 6e-6,
        },
        (temperature_k - 300.0) * 500.0,
        &PropertyKey::with_quantity(
            "thermal_conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn link() -> FreeSurfaceBodyLink {
    FreeSurfaceBodyLink {
        other_body_id: 43,
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

#[test]
fn g1_wall_spreading_matches_independent_three_store_heat_balance() {
    let liquid = conducting_card(1.0);
    let front = wall_card(0.5, 500.0);
    let rear = wall_card(1.5, 500.0);
    let front_chart = wall_chart(&front);
    let rear_chart = wall_chart(&rear);
    let front_context = point(400.0);
    let rear_context = point(480.0);
    with_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [
            cell(42, &front, &front_chart, &front_context, 400.0)
                .with_conduction_links(vec![link()]),
            cell(43, &rear, &rear_chart, &rear_context, 480.0),
        ];
        let before = thermal.total_enthalpy_j().unwrap()
            + bodies.iter().map(|body| body.enthalpy_j()).sum::<f64>();
        let report = thermal
            .step_with_body_conduction(&[], &[], controls(), &mut bodies, &[contact()], gate)
            .unwrap();

        // Three independent finite heat capacities, with both rates evaluated
        // at the SAME old state: liquid 0.002 J/K, each solid 0.001 J/K.
        // The liquid initially matches the front wall; only conduction through
        // that wall can deliver heat from the remote rear cell into the fluid.
        let mut t = [400.0, 400.0, 480.0];
        let mut last_solid_g = 0.0;
        for _ in 0..4 {
            let liquid_g = 1e-6
                / (0.0005 / (t[0] - 200.0) + 0.0005 / (0.5 * (t[1] - 200.0)));
            last_solid_g = 1e-6
                / (0.0005 / (0.5 * (t[1] - 200.0)) + 0.0005 / (1.5 * (t[2] - 200.0)));
            let into_liquid = 0.00025 * liquid_g * (t[1] - t[0]);
            let into_front = 0.00025 * last_solid_g * (t[2] - t[1]);
            t[0] += into_liquid / 0.002;
            t[1] += (into_front - into_liquid) / 0.001;
            t[2] -= into_front / 0.001;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(report.flow.external_heat_j, 0.0);
        near(temperature(thermal.material_at(WALL_LIQUID).unwrap()), t[0], 2e-13);
        assert!(t[0] > 400.0);
        for (body, expected) in bodies.iter().zip(&t[1..]) {
            near(body.specific_enthalpy_j_kg(), (expected - 300.0) * 500.0, 2e-13);
            near(body.mass_kg(), 2e-6, 2e-15);
        }
        near(report.flow.body_heat_j, 0.002 * (t[0] - 400.0), 2e-12);
        let body_report = report.bodies.unwrap();
        assert_eq!(body_report.solid_last_applied.len(), 1);
        let applied = &body_report.solid_last_applied[0];
        assert_eq!(applied.body_id, 42);
        assert_eq!(applied.link, link());
        assert_eq!(applied.materials[0].card_identity(), front.content_hash());
        assert_eq!(applied.materials[1].card_identity(), rear.content_hash());
        near(applied.conductance_w_k, last_solid_g, 2e-13);
        for (material, body) in body_report.solid_thermal_endpoint[0]
            .materials.iter().zip(&bodies)
        {
            assert_eq!(material, body.material());
        }
        assert!(body_report.combined_balance_defect_j.abs() < 1e-9);
        near(
            thermal.total_enthalpy_j().unwrap()
                + bodies.iter().map(|body| body.enthalpy_j()).sum::<f64>(),
            before,
            2e-12,
        );
    });
}

#[test]
fn g1_stiff_solid_links_use_native_stability_and_remain_internal_heat() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(20.0, 500.0);
    let chart = wall_chart(&solid);
    let context = point(400.0);
    with_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [
            cell(42, &solid, &chart, &context, 350.0)
                .with_conduction_links(vec![link()]),
            cell(43, &solid, &chart, &context, 450.0),
        ];
        let liquid_before = thermal.total_enthalpy_j().unwrap();
        let solid_before: f64 = bodies.iter().map(|body| body.enthalpy_j()).sum();
        let report = thermal
            .step_with_body_conduction(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut bodies,
                &[],
                gate,
            )
            .unwrap();

        // Independent two-capacity explicit reference, using each body's
        // heat-capacity/conductance bound; no liquid-wall contact exists.
        let mut t = [350.0, 450.0];
        let mut elapsed = 0.0_f64;
        let mut steps = 0;
        while elapsed < 0.001 {
            let g = 1e-6
                / (0.0005 / (20.0 * (t[0] - 200.0))
                    + 0.0005 / (20.0 * (t[1] - 200.0)));
            let dt = (0.001 - elapsed).min(0.8 * 0.001 / g);
            let q = dt * g * (t[1] - t[0]);
            t[0] += q / 0.001;
            t[1] -= q / 0.001;
            elapsed += dt;
            steps += 1;
            assert!(steps < 100);
        }
        assert_eq!(report.thermal.substeps, steps);
        assert!(steps > 1);
        for (body, expected) in bodies.iter().zip(t) {
            near(body.specific_enthalpy_j_kg(), (expected - 300.0) * 500.0, 2e-12);
        }
        assert_eq!(report.flow.external_heat_j, 0.0);
        assert!(report.flow.body_heat_j.abs() < 1e-12);
        assert!(report.bodies.unwrap().heat_to_liquid_j.abs() < 1e-12);
        near(thermal.total_enthalpy_j().unwrap(), liquid_before, 2e-13);
        near(
            bodies.iter().map(|body| body.enthalpy_j()).sum(),
            solid_before,
            2e-13,
        );
    });
}

#[test]
fn g4_body_graph_and_late_source_refusals_preserve_every_history() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(1.0, 500.0);
    let chart = wall_chart(&solid);
    let context = point(400.0);
    with_liquid(&liquid, true, |thermal, gate| {
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        for case in 0..5 {
            let mut first_link = link();
            if case == 0 {
                first_link.other_body_id = 99;
            } else if case == 1 {
                first_link.other_body_id = 42;
            } else if case == 3 {
                // This fits the solid by itself, but overcommits its surface
                // once the original 1e-6 m2 wetted liquid patch is included.
                first_link.area_m2 = 6e-6;
            }
            let mut bodies = [
                cell(42, &solid, &chart, &context, 400.0)
                    .with_conduction_links(vec![first_link]),
                cell(43, &solid, &chart, &context, 450.0),
            ];
            if case == 2 {
                bodies[1] = bodies[1].clone().with_conduction_links(vec![
                    FreeSurfaceBodyLink { other_body_id: 42, ..link() },
                ]);
            }
            let original: Vec<_> = bodies.iter()
                .map(|body| (body.specific_enthalpy_j_kg(), body.material().clone()))
                .collect();
            let mut control = controls();
            if case == 4 {
                control.maximum_faces = 1;
            }
            let error = thermal.step_with_body_conduction(
                &[], &[], control, &mut bodies, &[contact()], gate,
            ).unwrap_err();
            if case == 3 {
                assert!(matches!(error, FreeSurfaceEnthalpyError::Thermal(_)));
            } else {
                assert!(matches!(error, FreeSurfaceEnthalpyError::Invalid(_)));
            }
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            for (body, (h, source)) in bodies.iter().zip(&original) {
                assert_eq!(body.specific_enthalpy_j_kg(), *h);
                assert_eq!(body.material(), source);
            }
        }

        // Only solid/solid conduction warms this initially supported cell.
        // Its caloric chart remains valid, but its pinned k source exits.
        let limited = wall_card(1.0, 350.0);
        let limited_chart = wall_chart(&limited);
        let limited_context = point(350.0);
        let mut bodies = [
            cell(42, &limited, &limited_chart, &limited_context, 350.0)
                .with_conduction_links(vec![link()]),
            cell(43, &solid, &chart, &context, 450.0),
        ];
        let original: Vec<_> = bodies.iter()
            .map(|body| (body.specific_enthalpy_j_kg(), body.material().clone()))
            .collect();
        assert!(matches!(
            thermal.step_with_body_conduction(&[], &[], controls(), &mut bodies, &[], gate),
            Err(FreeSurfaceEnthalpyError::Material(_))
        ));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step_with_body_conduction(
            &[], &[], controls(), &mut bodies, &[], &cancelled,
        ).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        for (body, (h, source)) in bodies.iter().zip(&original) {
            assert_eq!(body.specific_enthalpy_j_kg(), *h);
            assert_eq!(body.material(), source);
        }
    });
}
