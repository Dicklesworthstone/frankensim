//! Original wall emissivity reaches finite gas, contact heat and native motion.
use super::*;
use fs_conduction::enthalpy::CellGeometry;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyRadiationFace, FreeSurfaceThermalBody,
};
use fs_material::thermal::source::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint,
};

fn emissive_wall(scale: f64, upper_k: f64) -> MaterialCard {
    let base = body::wall_card(0.5, 500.0);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.2 * scale), (500.0, 0.8 * scale)],
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, upper_k)
                .with("pressure", 100000.0, 100000.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic original solid emissivity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn radiant_face(surroundings_k: f64) -> FreeSurfaceBodyRadiationFace {
    let film = solid_face();
    FreeSurfaceBodyRadiationFace {
        wall: film.wall,
        gas: film.gas,
        area_m2: film.area_m2,
        surroundings_k,
    }
}

#[test]
fn g1_g3_wall_emissivity_closes_radiation_gas_contact_and_native_feedback() {
    let run = |scale: f64| {
        let liquid = conducting_card(1.0);
        let solid = emissive_wall(scale, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = parcel(350.0, 1e-7, 100.0);
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
                .with_gas_faces(vec![solid_face()])
                .with_radiation_faces(
                    vec![radiant_face(1500.0)],
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap()];
            let original = bodies[0]
                .material()
                .property(SURFACE_EMISSIVITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected;
            let before_liquid = thermal.total_enthalpy_j().unwrap();
            let before_solid = bodies[0].enthalpy_j();
            let before_gas = gas.state().specific_internal_energy_j_kg();
            let report = thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[body::contact()],
                    &[],
                    &mut gas,
                    &[exposed],
                    gas_controls(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;

            // Independent four-store recurrence. All rates use the same old
            // temperatures. The gas uses the integral of cv/R = 2 + 0.001T;
            // radiation uses the direct fourth-power difference, independently
            // of the production secant and existing thermal solver.
            let (mut twet, mut tl, mut ts, mut tg) = (400.0_f64, 400.0_f64, 450.0_f64, 350.0_f64);
            let (mut qgl_total, mut qgs_total, mut qsl_total, mut qr_total) = (0.0, 0.0, 0.0, 0.0);
            let r = 8.314_462_618_153_24 / 0.028;
            let mg = 100000.0 * 1e-7 / (r * 350.0);
            for _ in 0..4 {
                let kg = 2.6 * (0.5 + (100000.0 * tg / 350.0 - 1000.0) / 249000.0);
                let ks = 0.5 * (ts - 200.0);
                let gsl = 1e-6 / (0.0005 / (twet - 200.0) + 0.0005 / ks);
                let ggl = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / kg);
                let ggs = 1e-6 / (0.0005 / ks + 0.0005 / kg);
                let qsl = 0.00025 * gsl * (ts - twet);
                let qgl = 0.00025 * ggl * (tg - tl);
                let qgs = 0.00025 * ggs * (tg - ts);
                let eps = scale * (0.2 + 0.003 * (ts - 300.0));
                let qr = 0.00025 * eps * 5.670374419e-8 * 1e-6 * (1500.0_f64.powi(4) - ts.powi(4));
                twet += qsl / 0.002;
                tl += qgl / 0.002;
                ts += (qgs + qr - qsl) / 0.001;
                tg = (((2.0 + 0.001 * tg).powi(2) - 0.002 * (qgl + qgs) / (mg * r)).sqrt() - 2.0)
                    / 0.001;
                qgl_total += qgl;
                qgs_total += qgs;
                qsl_total += qsl;
                qr_total += qr;
            }
            let radiant = report.body_radiation.as_ref().unwrap();
            let endpoint = &radiant.thermal_endpoint[0];
            let walls = report.bodies.as_ref().unwrap();
            assert_eq!(report.thermal.substeps, 4);
            assert!(
                report.radiation.is_none(),
                "wall radiation is not direct liquid supply"
            );
            near(endpoint.solid_temperature_k, ts, 2e-11);
            near(
                endpoint.emissivity,
                scale * (0.2 + 0.003 * (ts - 300.0)),
                2e-11,
            );
            assert_eq!(
                endpoint.emissivity_source.answer().receipt.selected,
                original
            );
            assert_eq!(
                bodies[0]
                    .material()
                    .property(SURFACE_EMISSIVITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                original
            );
            assert_eq!(bodies[0].material().query().axes()["pressure"], 100000.0);
            near(radiant.heat_to_solids_j, qr_total, 2e-11);
            near(walls.heat_from_radiation_j, qr_total, 2e-11);
            near(walls.heat_from_gas_j, qgs_total, 2e-11);
            near(walls.heat_to_liquid_j, qsl_total, 2e-11);
            near(report.flow.body_heat_j, qsl_total, 2e-11);
            near(report.flow.external_heat_j, qgl_total, 2e-11);
            near(
                report.thermal.external_heat_j,
                qgl_total + qgs_total + qr_total,
                2e-11,
            );
            near(
                bodies[0].enthalpy_j() - before_solid,
                qgs_total + qr_total - qsl_total,
                2e-11,
            );
            near(gas.state().state().temperature, tg, 2e-10);
            near(gas.state().state().pressure, 100000.0 * tg / 350.0, 2e-10);
            let gas_change =
                gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - before_gas);
            assert!((gas_change + qgl_total + qgs_total).abs() < 1e-10);
            assert!(
                (thermal.total_enthalpy_j().unwrap() - before_liquid + bodies[0].enthalpy_j()
                    - before_solid
                    + gas_change
                    - qr_total)
                    .abs()
                    < 1e-9
            );
            assert!(
                report
                    .finite_gas
                    .unwrap()
                    .combined_thermal_balance_defect_j
                    .abs()
                    < 1e-9
            );
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .material_at(body::contact().liquid)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                qr_total,
            )
        })
    };
    let reflecting = run(0.0);
    let dim = run(0.5);
    let bright = run(1.0);
    assert_eq!(reflecting.2, 0.0);
    assert!(bright.2 > dim.2 && bright.1 < dim.1);
    assert_ne!(
        bright.0, dim.0,
        "wall emissivity data must reach actual VOF populations"
    );
}

#[test]
fn g1_g4_radiation_capacity_controls_substeps_and_failed_budget_keeps_history() {
    let liquid = conducting_card(1.0);
    let solid = emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(400.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [FreeSurfaceThermalBody::try_new(
            42,
            ThermalMaterialPoint {
                material: &solid,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            &chart,
            CellGeometry {
                start_m: [-1e-9, 0.0035, 0.0035],
                end_m: [0.0, 0.0035, 0.0035],
                volume_m3: 1e-15,
                surface_area_m2: 2e-6,
            },
            50000.0,
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
        .with_radiation_faces(
            vec![radiant_face(480.0)],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()];
        let before = bodies[0].enthalpy_j();
        let source = bodies[0].material().clone();
        let bits = thermal.flow().lattice().grid().state_bits();
        let limited = FreeSurfaceConductionControls {
            maximum_substeps: 1,
            ..controls()
        };
        assert!(
            thermal
                .step_with_body_conduction(&[], &[], limited, &mut bodies, &[], gate)
                .is_err()
        );
        assert_eq!(bodies[0].enthalpy_j(), before);
        assert_eq!(bodies[0].material(), &source);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.accepted_steps(), 0);
        let report = thermal
            .step_with_body_conduction(&[], &[], controls(), &mut bodies, &[], gate)
            .unwrap();
        let exchange = report.body_radiation.unwrap();
        assert!(
            report.thermal.substeps > 4,
            "solid radiation must enter the original combined CFL"
        );
        let temperature = exchange.thermal_endpoint[0].solid_temperature_k;
        assert!((400.0..=480.0).contains(&temperature));
        near(
            bodies[0].enthalpy_j() - before,
            exchange.heat_to_solids_j,
            2e-12,
        );
        near(
            report.thermal.external_heat_j,
            exchange.heat_to_solids_j,
            2e-12,
        );
        assert!(report.flow.body_heat_j.abs() < 1e-18);
        assert_eq!(report.flow.external_heat_j, 0.0);
    });
}

#[test]
fn g4_missing_source_geometry_and_late_emissivity_exit_are_atomic() {
    let missing = body::wall_card(0.5, 500.0);
    let missing_chart = body::wall_chart(&missing);
    let context = point(450.0);
    assert!(
        body::wall(&missing, &missing_chart, &context, 450.0)
            .with_radiation_faces(
                vec![radiant_face(1500.0)],
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
    );
    let liquid = conducting_card(1.0);
    let solid = emissive_wall(1.0, 450.00001);
    let chart = body::wall_chart(&solid);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = parcel(450.0, 1e-7, 100.0);
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
            .with_gas_faces(vec![solid_face()])
            .with_radiation_faces(
                vec![radiant_face(1500.0)],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()];
        let wall_before = bodies[0].enthalpy_j();
        let source_before = bodies[0].material().clone();
        let gas_before = gas.state().clone();
        let liquid_before = thermal.total_enthalpy_j().unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        assert!(
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[],
                    &[],
                    &mut gas,
                    &[],
                    gas_controls(),
                    false,
                    gate,
                )
                .is_err(),
            "the original emissivity domain must refuse after the first heated candidate"
        );
        assert_eq!(bodies[0].enthalpy_j(), wall_before);
        assert_eq!(bodies[0].material(), &source_before);
        assert_eq!(gas.state(), &gas_before);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), liquid_before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.accepted_steps(), 0);
    });
    let solid = emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let patches: Vec<_> = (1..=6)
            .map(|z| FreeSurfaceBodyRadiationFace {
                wall: [-1, 0, z],
                gas: [0, 0, z as usize],
                ..radiant_face(480.0)
            })
            .collect();
        for faces in [vec![radiant_face(480.0), radiant_face(480.0)], patches] {
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
                .with_radiation_faces(faces, MaterialPropertySelection::SingleClaimOnly)
                .unwrap()];
            let h = bodies[0].specific_enthalpy_j_kg();
            assert!(
                thermal
                    .step_with_body_conduction(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        gate
                    )
                    .is_err()
            );
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), h);
            assert_eq!(thermal.accepted_steps(), 0);
        }
    });
}
