//! G1/G3 gradients through the real mixed-convection import/solve pipeline.
//! Finite differences perturb physical controls and re-solve every state law.

use super::*;

fn request(spec: &mut ProjectSpec, name: &str) {
    spec.outputs.as_mut().unwrap().push(OutputRequest {
        name: name.to_string(), kind: "report".to_string(), region: None,
    });
}

fn coefficient(receipt: &J, target: &str) -> f64 {
    let rows = receipt.get("nominal_adjoint").unwrap()
        .get("parameters").unwrap().as_array().unwrap();
    let selected: Vec<_> = rows.iter()
        .filter(|row| row.str_field("target") == Some(target)).collect();
    assert_eq!(selected.len(), 1, "{target}: {rows:?}");
    n(selected[0], "derivative")
}

fn node(field: &J, selected: usize) -> f64 {
    field.get("temperature").unwrap().as_array().unwrap()[selected].as_f64().unwrap()
}

fn perturb(spec: &mut ProjectSpec, target: &str, delta: f64) {
    if target == "power" {
        spec.power.as_mut().unwrap()[0].watts.value += delta;
        return;
    }
    let cooling = spec.cooling.as_mut().unwrap();
    if target == "fan-speed-ratio" {
        cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio += delta;
        return;
    }
    let setup = cooling.conduction.as_mut().unwrap();
    if target == "radiation-reservoir-temperature" {
        setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta;
        return;
    }
    for boundary in &mut setup.boundaries {
        match (&mut boundary.condition, target) {
            (ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. },
                "natural-convection-ambient") => {
                ambient_temperature.value += delta;
                return;
            }
            (ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. },
                "air-inlet-temperature") => {
                inlet_temperature.value += delta;
                return;
            }
            (ThermalBoundaryCondition::FixedTemperature { temperature },
                "fixed-temperature") => {
                temperature.value += delta;
                return;
            }
            _ => {}
        }
    }
    panic!("unknown fixture control {target}");
}

fn check_control(
    spec: &ProjectSpec, bytes: &[u8], cards: &CardPackSet,
    nominal: &J, selected: usize, target: &str, step: f64,
) {
    let values: Vec<_> = [-1.0, 1.0].into_iter().map(|sign| {
        let mut varied = spec.clone();
        perturb(&mut varied, target, sign * step);
        node(&run(&varied, bytes, cards).1, selected)
    }).collect();
    // Keep the SAME selected node even if a nearby maximum changes location.
    let expected = (values[1] - values[0]) / (2.0 * step);
    let actual = coefficient(nominal, target);
    assert!(expected.is_finite() && actual.is_finite());
    assert!((actual - expected).abs() <= 5e-4 * expected.abs().max(0.005),
        "{target}: complete mixed adjoint {actual:e}, physical re-solves {expected:e}");
}

#[test]
fn g1_mixed_adjoint_retains_the_primal_and_differentiates_both_heat_paths() {
    let bytes = tetra_stl();
    // The second case exercises a COOLED natural wall, smooth k(T), and a
    // separate hotter radiation reservoir in the same nonsymmetric Jacobian.
    for (ambient, radiating) in [(293.15, false), (300.0, true)] {
        let mut spec = project(&bytes, ambient);
        spec.power.as_mut().unwrap()[0].duty = 0.37;
        let base_cards = if radiating {
            CardPackSet::admit(vec![raw_pack(
                CardPackKind::Material, "mixed-adjoint-smooth-k",
                material_pack_bytes_with_property(
                    "AA6061", "mixed-adjoint-smooth-k", 200.0, 450.0,
                    fs_matdb::PropertyValue::Curve {
                        abscissa: "T".to_string(),
                        abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
                        knots: vec![(200.0, 8.0), (450.0, 80.0)],
                        dims: CONDUCTIVITY_DIMS,
                    },
                    fs_matdb::InterpolationPolicy::LinearInside,
                ),
            )]).unwrap()
        } else { fixture_cards() };
        rebind_to(&mut spec, &base_cards);
        let cards = if radiating {
            let (cards, card) = radiation_product::emissivity_cards(&base_cards, 0.85);
            radiation_product::declare(&mut spec, &card, "natural-wall", 310.0);
            let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .radiation.as_mut().unwrap();
            radiation.temperature_tolerance.value = 1e-10;
            radiation.heat_tolerance.value = 1e-10;
            cards
        } else { base_cards };
        spec.solver.as_mut().unwrap().tolerance_rel = 1e-11;
        let (baseline, baseline_field) = run(&spec, &bytes, &cards);
        let mut requested = spec.clone();
        request(&mut requested, "temperature-max-adjoint");
        let (nominal, field) = run(&requested, &bytes, &cards);
        assert_eq!(baseline_field.get("temperature"), field.get("temperature"),
            "the adjoint must use the unchanged accepted primal");
        for key in ["energy", "natural_convection", "conjugate", "radiation"] {
            assert_eq!(baseline.get(key), nominal.get(key), "{key}");
        }
        let adjoint = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(adjoint.str_field("authority"), Some("Estimated"));
        assert_eq!(adjoint.str_field("mode"), Some(if radiating {
            "natural-radiation-full-air-feedback"
        } else { "natural-convection-full-air-feedback" }));
        assert!(n(adjoint, "true_relative_residual") < 1e-10);
        let natural = &nominal.get("natural_convection").unwrap()
            .get("laws").unwrap().as_array().unwrap()[0];
        assert_eq!(n(natural, "heat_rate_w") < 0.0, radiating);
        let selected = n(adjoint, "selected_vertex") as usize;
        for (target, step) in [
            ("power", 0.02),
            ("natural-convection-ambient", 0.005),
            ("air-inlet-temperature", 0.005),
            ("fan-speed-ratio", 0.002),
        ] {
            check_control(&spec, &bytes, &cards, &nominal, selected, target, step);
        }
        if radiating {
            check_control(&spec, &bytes, &cards, &nominal, selected,
                "radiation-reservoir-temperature", 0.005);
        }
    }
}

#[test]
fn g3_mixed_boundary_adjoint_includes_prescribed_nodes_in_the_full_lift() {
    let bytes = tetra_stl();
    let mut spec = project(&bytes, 300.0);
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-11;
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".to_string(), name: "cold-anchor".to_string(),
        display: "Manufactured prescribed-temperature face".to_string(), expect_id: None,
    });
    spec.assignments.as_mut().unwrap().push(GeometryAssignment {
        artifact: "enclosure".to_string(), target: "cold-anchor".to_string(),
        length_unit: "m".to_string(),
        selector: MeshSelector::HalfSpace {
            normal: [1.0, 0.0, 0.0], offset: 0.0,
            side: HalfSpaceSide::AtMost, tolerance: 0.0,
        },
        allow_overlap: true,
    });
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries
        .push(ThermalBoundary {
            target: "cold-anchor".to_string(),
            condition: ThermalBoundaryCondition::FixedTemperature {
                temperature: QtyAny::new(280.0, fs_project::spec::dims::TEMPERATURE),
            },
        });
    let cards = fixture_cards();
    let (baseline, before) = run(&spec, &bytes, &cards);
    let mut requested = spec.clone();
    request(&mut requested, "temperature-max-boundary-adjoint");
    let (nominal, after) = run(&requested, &bytes, &cards);
    assert_eq!(before.get("temperature"), after.get("temperature"));
    assert_eq!(baseline.get("energy"), nominal.get("energy"));
    let adjoint = nominal.get("nominal_adjoint").unwrap();
    let selected = n(adjoint, "selected_vertex") as usize;
    assert!(node(&after, selected) > 280.0,
        "a free-node goal must exercise the eliminated-temperature lift");
    let derivative = coefficient(&nominal, "fixed-temperature");
    assert!(derivative > 0.0 && derivative < 1.0);
    check_control(&spec, &bytes, &cards, &nominal, selected, "fixed-temperature", 0.005);
}
