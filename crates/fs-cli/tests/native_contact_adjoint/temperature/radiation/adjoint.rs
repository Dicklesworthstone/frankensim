//! Complete live-contact/ambient-radiation derivatives through native solves.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

fn source_material(
    name: &str,
    property: &str,
    dims: fs_qty::Dims,
    value: PropertyValue,
    interpolation: InterpolationPolicy,
) -> (Vec<u8>, String, String) {
    let bytes = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: name.into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
        },
        pack(name, property, dims, value, interpolation),
    ).unwrap().to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Material, source: name.into(), bytes: bytes.clone(), expect: None,
    }]).unwrap();
    (bytes, admitted.materials()[0].card().to_hex(), admitted.materials()[0].identity().to_string())
}

fn contact_source(scale: f64, surface: bool) -> (Vec<u8>, String, String) {
    if !surface {
        return temperature_pack_from_knots(vec![(200.0, 0.2 * scale), (450.0, 0.45 * scale)]);
    }
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [([200.0, 0.8e6], 0.18), ([450.0, 0.8e6], 0.43),
            ([200.0, 1.2e6], 0.22), ([450.0, 1.2e6], 0.47)]
            .into_iter().map(|(point, value)| (point, scale * value)).collect(),
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

/// Re-author source cards and import the ordinary project for every physical
/// perturbation. Only the hot material changes under its regional control.
fn solve(
    f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool,
    contact_scale: f64, conductivity_scale: f64, emissivity: f64,
) -> (J, J) {
    let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
    let kd = fs_conduction::CONDUCTIVITY_DIMS;
    let (bytes, card, state) = source_material("ambient-contact-hot", "thermal-conductivity", kd,
        PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: dims::TEMPERATURE,
            knots: vec![(200.0, 2.0 * conductivity_scale), (450.0, 17.0 * conductivity_scale)], dims: kd },
        InterpolationPolicy::LinearInside);
    let hot = f.dir.join(format!("ambient-hot-{ordinal}.fsmcdpk"));
    std::fs::write(&hot, bytes).unwrap();
    let binding = changed.project.materials.as_mut().unwrap().iter_mut()
        .find(|binding| binding.region == "hot").unwrap();
    binding.card = card;
    binding.state = state;
    binding.claim = None;
    let (bytes, card, _) = source_material("ambient-contact-gray", fs_conduction::SURFACE_EMISSIVITY_PROPERTY,
        fs_qty::Dims::NONE, PropertyValue::Scalar { value: emissivity, dims: fs_qty::Dims::NONE },
        InterpolationPolicy::ConstantWithinValidity);
    let gray = f.dir.join(format!("ambient-gray-{ordinal}.fsmcdpk"));
    std::fs::write(&gray, bytes).unwrap();
    let patch = &mut changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap().surfaces[0];
    patch.card = card;
    patch.claim = None;
    let (source, contact, ledger) = prepare_fixed_pressure_surface(&changed, ordinal, output,
        contact_source(contact_scale, surface));
    let result = fs_cli::run(vec!["--json".into(), "solve".into(),
        source.display().to_string(), ledger.display().to_string(),
        "--materials".into(), f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(), hot.display().to_string(),
        "--materials".into(), gray.display().to_string(),
        "--interfaces".into(), contact.display().to_string()]);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let run = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn value(field: &J, vertex: usize) -> f64 {
    field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap()
}

fn derivative(report: &J, target: &str, entity: &str) -> f64 {
    let rows: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target") == Some(target)
            && row.str_field("entity") == Some(entity)).collect();
    assert_eq!(rows.len(), 1, "one original physical control: {target}/{entity}");
    number(rows[0], "derivative")
}

#[test]
fn ambient_live_contact_controls_match_physical_resolves_in_vacuum_and_mixed_cooling() {
    // The vacuum case exercises curve/Newton contact on a radiation-only patch.
    // The source-surface/Picard case couples radiation, natural walls and air.
    for mixed in [false, true] {
        let mut f = super::fixture(!mixed);
        if mixed {
            f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[1].condition =
                B::NaturalConvection {
                    characteristic_length: QtyAny::new(1.0, dims::LENGTH),
                    ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
                    correlation: "convection.churchill-chu-vertical-plate".into(),
                };
            // Fix gas properties/hydraulics while independently changing the
            // physical inlet and natural reservoir within the same envelope.
            let envelope = f.project.envelope.as_mut().unwrap();
            envelope.ambient_lo.value = 296.0;
            envelope.ambient_hi.value = 298.0;
        }
        let (plain, plain_field) = solve(&f, 220, None, mixed, 1.0, 1.0, 0.73);
        let (nominal, field) = solve(&f, 221, Some(OUTPUT), mixed, 1.0, 1.0, 0.73);
        assert_eq!(plain_field.get("temperature"), field.get("temperature"));
        for key in ["energy", "interfaces", "temperature_contact", "radiation", "conjugate", "natural_convection"] {
            assert_eq!(
                forward_section(&plain, key),
                forward_section(&nominal, key),
                "adjoint changed primal {key}"
            );
        }
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("output"), Some(OUTPUT));
        assert_eq!(report.str_field("authority"), Some("Estimated"));
        assert!(number(report, "true_relative_residual") < 1e-10);
        let selected = number(report, "selected_vertex") as usize;
        assert_eq!(number(report, "value_k"), value(&field, selected));
        assert_eq!(nominal.get("temperature_contact").unwrap().str_field("algorithm"),
            Some(if mixed { "bounded-contact-picard" } else { "coupled-newton" }));
        let radiation = nominal.get("radiation").unwrap();
        let convective = number(radiation, "convective_out_w");
        let radiative = number(radiation, "nonlinear_radiative_out_w");
        assert!(radiative.abs() > 0.01);
        assert!((convective + radiative - 4.44).abs() < 1e-7);
        if mixed {
            let natural: f64 = nominal.get("natural_convection").unwrap().get("laws").unwrap()
                .as_array().unwrap().iter().map(|row| number(row, "heat_rate_w")).sum();
            let air = nominal.get("conjugate").unwrap();
            assert!((number(air, "solid_total_w") + natural - convective).abs() < 1e-7);
            assert!((number(air, "air_total_w") - number(air, "solid_total_w")).abs()
                <= number(air, "balance_tolerance_w"));
        } else {
            assert_eq!(convective, 0.0);
            assert_eq!(nominal.get("conjugate"), Some(&J::Null));
        }
        let mut controls = vec![(TARGET, "cold-hot-joint", 0.002),
            ("power", "hot", 0.01), ("conductivity-multiplier", "hot", 0.002),
            ("radiation-emissivity", "cold-radiator", 0.0005),
            ("radiation-reservoir-temperature", "cold-radiator", 0.01)];
        if mixed {
            controls.extend([("air-inlet-temperature", "cold", 0.002),
                ("fan-speed-ratio", "fixture-bank", 0.002),
                ("natural-convection-ambient", "hot", 0.002)]);
        }
        for (case, (target, entity, step)) in controls.into_iter().enumerate() {
            let actual = derivative(report, target, entity);
            assert!(actual.abs() > 1e-7, "{target} must exercise the complete dual");
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let delta = sign * step;
                match target {
                    "power" => changed.project.power.as_mut().unwrap()[0].watts.value += delta,
                    "fan-speed-ratio" => changed.project.cooling.as_mut().unwrap().fan_system.as_mut()
                        .unwrap().banks[0].speed_ratio += delta,
                    "radiation-reservoir-temperature" => changed.project.cooling.as_mut().unwrap()
                        .conduction.as_mut().unwrap().radiation.as_mut().unwrap().surfaces[0]
                        .reservoir_temperature.value += delta,
                    "air-inlet-temperature" => {
                        let B::AirflowConvection { inlet_temperature, .. } = &mut changed.project.cooling
                            .as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                            else { panic!("original transported-air inlet"); };
                        inlet_temperature.value += delta;
                    }
                    "natural-convection-ambient" => {
                        let B::NaturalConvection { ambient_temperature, .. } = &mut changed.project.cooling
                            .as_mut().unwrap().conduction.as_mut().unwrap().boundaries[1].condition
                            else { panic!("original natural reservoir"); };
                        ambient_temperature.value += delta;
                    }
                    _ => {}
                }
                let (_, field) = solve(&changed, 222 + 2 * case + side, None, mixed,
                    1.0 + if target == TARGET { delta } else { 0.0 },
                    1.0 + if target == "conductivity-multiplier" { delta } else { 0.0 },
                    0.73 + if target == "radiation-emissivity" { delta } else { 0.0 });
                // Keep the nominal selected node, including when a perturbed
                // resolve changes which other node is the hottest.
                values.push(value(&field, selected));
            }
            let expected = (values[1] - values[0]) / (2.0 * step);
            assert!((actual - expected).abs() < 5e-5 * expected.abs().max(1e-3),
                "mixed={mixed} {target}/{entity}: complete adjoint {actual:e} vs resolve {expected:e}");
        }
        let (ordinary, ordinary_field) = solve(&f, 240, Some("temperature-max-adjoint"), mixed, 1.0, 1.0, 0.73);
        assert_eq!(ordinary_field.get("temperature"), field.get("temperature"));
        let ordinary = ordinary.get("nominal_adjoint").unwrap();
        let common: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") != Some(TARGET)).collect();
        assert_eq!(common, ordinary.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
    }
}

#[test]
fn live_contact_ambient_admission_keeps_steady_pointwise_derivatives_refused() {
    let mut f = super::fixture(true);
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap().spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    let (source, contact, ledger) = prepare(&f, 242, Some(OUTPUT));
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let result = fs_cli::run(vec!["--json".into(), "solve".into(), source.display().to_string(),
        ledger.display().to_string(), "--materials".into(), f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(), root.join("data/reference-project/gray-surface.fsmcdpk").display().to_string(),
        "--interfaces".into(), contact.display().to_string()]);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    let diagnostic = format!("{}{}", result.stdout, result.stderr);
    assert!(diagnostic.contains("cli-solve-conduction-pointwise-radiation"), "{diagnostic}");
    assert!(!diagnostic.contains("run_receipt"));
}
