//! Complete capacity/contact/radiation reports against fresh physical histories.
use super::*;

fn ambient_fixture(peak: bool) -> Fixture {
    let mut f = fixture(true, 0.5);
    f.project.power.as_mut().unwrap()[0].duty = 1.0;
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let reference = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
        root.join("data/reference-project/cooling-radiation.fsim")).unwrap()).unwrap().decoded.spec;
    let mut radiation = reference.cooling.unwrap().conduction.unwrap().radiation.unwrap();
    radiation.surfaces[0].name = "cold-radiator".into();
    radiation.surfaces[0].target = "cold".into();
    radiation.surfaces[0].reservoir_temperature.value = 290.0;
    radiation.temperature_tolerance.value = 1e-11;
    radiation.heat_tolerance.value = 1e-10;
    radiation.max_iterations = 256;
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.radiation = Some(radiation);
    setup.boundaries[0].condition = if peak {
        ThermalBoundaryCondition::Convection {
            coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
            reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
        }
    } else {
        ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        }
    };
    if !peak {
        // A distinct prescribed trace exercises the complete fixed-node lift
        // while the objective and radiator remain on the free cold body.
        setup.boundaries[1].condition = ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(310.0, dims::TEMPERATURE),
        };
        f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
        f.project.power.as_mut().unwrap()[0].region = "cold".into();
        f.project.power.as_mut().unwrap()[0].watts.value = 24.0;
    }
    let policy = time(&mut f);
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(300.0, dims::TEMPERATURE));
    if peak {
        policy.power_schedules[0].steps[0].watts.value = 80.0;
        policy.power_schedules[0].steps[1].watts.value = 0.0;
    } else {
        policy.power_schedules.clear();
    }
    f
}

fn solve_ambient(f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool,
    emissivity: f64) -> (J, J)
{
    let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
    let name = "capacity-contact-gray";
    let bytes = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: name.into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, pack(name, fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_qty::Dims::NONE,
        PropertyValue::Scalar { value: emissivity, dims: fs_qty::Dims::NONE },
        InterpolationPolicy::ConstantWithinValidity)).unwrap().to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Material, source: name.into(), bytes: bytes.clone(), expect: None,
    }]).unwrap();
    let patch = &mut changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap().surfaces[0];
    patch.card = admitted.materials()[0].card().to_hex();
    patch.claim = None;
    let gray = f.dir.join(format!("capacity-gray-{ordinal}.fsmcdpk"));
    std::fs::write(&gray, bytes).unwrap();
    let (source, contact, ledger) = prepare_with_pack(&changed, ordinal, output, contact_pack(surface));
    let result = fs_cli::run(vec!["--json".into(), "solve".into(),
        source.display().to_string(), ledger.display().to_string(),
        "--materials".into(), f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(), gray.display().to_string(),
        "--interfaces".into(), contact.display().to_string()]);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let completed = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn radiation_row(derivative: &J) -> &J {
    let radiation = derivative.get("radiation").unwrap();
    assert_eq!(radiation.str_field("method"), Some("complete-area-mean-gray-capacity-tangent"));
    assert_eq!(radiation.str_field("residual_target"),
        Some("original-accepted-forward-physical-tolerance"));
    assert_eq!(n(radiation, "control_steps"), n(derivative, "state_index"));
    assert_eq!(n(radiation, "patch_count"), 1.0);
    let rows = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row.str_field("name"), Some("cold-radiator"));
    assert_eq!(row.str_field("target"), Some("cold"));
    assert_eq!(n(row, "emissivity"), 0.73);
    assert_eq!(n(row, "ambient_temperature_k"), 290.0);
    row
}

fn source_scale(derivative: &J) -> f64 {
    derivative.get("intervals").unwrap().as_array().unwrap().iter()
        .map(|row| n(row, "dtemperature_dsource_multiplier_k")).sum()
}

fn future_is_zero(derivative: &J) {
    let selected = n(derivative, "state_index") as usize;
    let intervals = derivative.get("intervals").unwrap().as_array().unwrap();
    assert!(selected < intervals.len(), "the selected peak must precede the final interval");
    for interval in &intervals[selected..] {
        assert_eq!(n(interval, "dtemperature_dsource_multiplier_k"), 0.0);
        for row in interval.get("regional_power").unwrap().as_array().unwrap() {
            assert_eq!(n(row, "dtemperature_ddelivered_watts_k_per_w"), 0.0);
        }
    }
}

#[test]
fn ambient_live_contact_final_and_earlier_peak_controls_match_full_physical_histories() {
    // Final: curve/Newton, vacuum radiator and an independent fixed bath.
    // Peak: source-surface/Picard, original Robin and radiation on one trace.
    for peak in [false, true] {
        let f = ambient_fixture(peak);
        let surface = peak;
        let output = if peak { PEAK } else { FINAL };
        let (plain, plain_field) = solve_ambient(&f, 170, None, surface, 0.73);
        let (receipt, field) = solve_ambient(&f, 171, Some(output), surface, 0.73);
        unchanged(&plain, &plain_field, &receipt, &field);
        assert_eq!(plain.get("transient").unwrap().get("radiation"),
            receipt.get("transient").unwrap().get("radiation"));
        let derivative = report(&receipt);
        let radiation = radiation_row(derivative);
        assert_eq!(n(derivative, "value_k"), qoi(&plain, peak));
        let mut controls = vec![("source", source_scale(derivative), 0.002),
            ("capacity", n(derivative, "dtemperature_dcapacity_multiplier_k"), 0.002),
            ("initial", n(derivative, "dtemperature_duniform_initial_k_per_k"), 0.002),
            ("emissivity", n(radiation, "dtemperature_demissivity_k"), 0.0005),
            ("radiation-bath", n(radiation, "dtemperature_dreservoir_temperature_k_per_k"), 0.002)];
        if peak {
            assert!(n(derivative, "time_s") > 0.0 && n(derivative, "time_s") < 2.0);
            assert!(qoi(&plain, true) > qoi(&plain, false) + 1e-3);
            future_is_zero(derivative);
            let robin = derivative.get("robin").unwrap().as_array().unwrap();
            assert_eq!(robin.len(), 1);
            assert_eq!(robin[0].str_field("region"), Some("cold"));
            controls.push(("log-htc", n(&robin[0], "dtemperature_dlog_htc_k"), 0.002));
        } else {
            assert!(derivative.get("robin").unwrap().as_array().unwrap().is_empty(),
                "radiation-only feedback must not invent a convection control");
            if let Some(flux) = derivative.get("prescribed_heat_flux") {
                assert!(flux.get("boundaries").unwrap().as_array().unwrap().iter()
                    .all(|row| row.str_field("target") != Some("cold")),
                    "a radiation-only zero-flux trace has no admitted nonzero flux control");
            }
            let baths = derivative.get("prescribed_temperature").unwrap()
                .get("boundaries").unwrap().as_array().unwrap();
            assert_eq!(baths.len(), 1);
            assert_eq!(baths[0].str_field("target"), Some("hot"));
            controls.push(("fixed-bath", n(&baths[0],
                "dtemperature_dprescribed_temperature_k_per_k"), 0.002));
        }
        let selected = n(derivative, "active_vertex") as usize;
        for (coordinate, (control, actual, step)) in controls.into_iter().enumerate() {
            assert!(actual.abs() > 1e-7, "active physical {control} control, peak={peak}");
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let delta = sign * step;
                let mut emissivity = 0.73;
                match control {
                    "source" => {
                        for row in changed.project.power.as_mut().unwrap() { row.watts.value *= 1.0 + delta; }
                        for schedule in &mut time(&mut changed).power_schedules {
                            for row in &mut schedule.steps { row.watts.value *= 1.0 + delta; }
                        }
                    }
                    "capacity" => for row in &mut time(&mut changed).capacities {
                        row.volumetric_heat_capacity.value *= 1.0 + delta;
                    },
                    "initial" => time(&mut changed).initial_temperature =
                        TransientInitialTemperature::Uniform(QtyAny::new(300.0 + delta, dims::TEMPERATURE)),
                    "emissivity" => emissivity += delta,
                    "radiation-bath" => changed.project.cooling.as_mut().unwrap().conduction.as_mut()
                        .unwrap().radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta,
                    "log-htc" => {
                        let ThermalBoundaryCondition::Convection { coefficient, .. } = &mut changed.project
                            .cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                            else { panic!("original Robin law"); };
                        coefficient.value *= delta.exp();
                    }
                    "fixed-bath" => {
                        let ThermalBoundaryCondition::FixedTemperature { temperature } = &mut changed.project
                            .cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[1].condition
                            else { panic!("original prescribed trace"); };
                        temperature.value += delta;
                    }
                    _ => unreachable!(),
                }
                let (changed_receipt, changed_field) = solve_ambient(&changed,
                    172 + 2 * coordinate + side, None, surface, emissivity);
                if peak {
                    let sample = changed_receipt.get("transient").unwrap().get("sampled_peak").unwrap();
                    assert_eq!(n(sample, "time_s"), n(derivative, "time_s"), "physical peak branch changed");
                    assert_eq!(n(sample, "state"), n(derivative, "state_index"));
                    values.push(qoi(&changed_receipt, true));
                } else {
                    // Preserve the nominal selected node if another node wins
                    // the maximum in one independently perturbed final field.
                    values.push(changed_field.get("temperature").unwrap().as_array().unwrap()[selected]
                        .as_f64().unwrap());
                }
            }
            agrees(actual, (values[1] - values[0]) / (2.0 * step));
        }
    }
}

#[test]
fn ambient_live_contact_initial_peak_has_zero_future_radiation_and_storage_controls() {
    let mut f = ambient_fixture(true);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let policy = time(&mut f);
    policy.power_schedules.clear();
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(320.0, dims::TEMPERATURE));
    let (plain, plain_field) = solve_ambient(&f, 190, None, false, 0.73);
    let (receipt, field) = solve_ambient(&f, 191, Some(PEAK), false, 0.73);
    unchanged(&plain, &plain_field, &receipt, &field);
    let derivative = report(&receipt);
    assert_eq!(n(derivative, "state_index"), 0.0);
    assert_eq!(n(derivative, "value_k"), 320.0);
    assert_eq!(n(derivative, "adjoint_krylov_iterations"), 0.0);
    assert_eq!(n(derivative, "dtemperature_duniform_initial_k_per_k"), 1.0);
    assert_eq!(n(derivative, "dtemperature_dcapacity_multiplier_k"), 0.0);
    assert_eq!(source_scale(derivative), 0.0);
    let radiation = radiation_row(derivative);
    assert_eq!(n(radiation, "dtemperature_demissivity_k"), 0.0);
    assert_eq!(n(radiation, "dtemperature_dreservoir_temperature_k_per_k"), 0.0);
    future_is_zero(derivative);
}
