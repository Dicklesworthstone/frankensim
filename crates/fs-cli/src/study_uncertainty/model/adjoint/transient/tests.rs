use super::*;
use fs_project::{TransientEnthalpy, TransientPowerStep, TransientRegionPower, spec::dims};
use fs_qty::QtyAny;

fn base() -> ProjectSpec {
    fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-transient.fsim"
    )))
    .unwrap()
    .decoded
    .spec
}

fn parameter(target: Target, entity: &str) -> UniformParameter {
    UniformParameter {
        name: format!("{target:?}"),
        entity: entity.into(),
        target,
        low: 1.0,
        high: 2.0,
    }
}

fn request(name: &str) -> OutputRequest {
    OutputRequest {
        name: name.into(),
        kind: "report".into(),
        region: None,
    }
}

fn receipt(enthalpy: bool, controls: &str) -> J {
    let schema = if enthalpy {
        "fs-cli-native-enthalpy-adjoint-v1"
    } else {
        "fs-cli-native-transient-adjoint-v1"
    };
    let time_schema = if enthalpy {
        "fs-cli-transient-enthalpy-v1"
    } else {
        "fs-cli-transient-conduction-v1"
    };
    J::parse(&format!(r#"{{"transient":{{"schema":"{time_schema}","status":"completed","qoi_time":"final","final_time_s":2,"fine_steps":8,"temporal_error":{{"qoi":"temperature-max","region":"solid","fine_final_k":301}},"adjoint":{{"schema":"{schema}","output":"temperature-final-adjoint","value_k":301,"time_s":2,"state_index":8,"max_reverse_steps":8,"maximum_dual_relative_residual":0.000000001{controls}}}}}}}"#)).unwrap()
}

const POWER: &str =
    r#", "regional_power":[{"region":"solid","dtemperature_duniform_added_watts_k_per_w":2.0}]"#;

#[test]
fn g0_calibration_requests_final_history_without_changing_physics_or_existing_intent() {
    for existing in [None, Some(FINAL), Some(FINAL_MATERIAL)] {
        let mut spec = base();
        if let Some(name) = existing {
            spec.outputs.as_mut().unwrap().push(request(name));
        }
        let mut expected = spec.clone();
        if existing.is_none() {
            expected.outputs.as_mut().unwrap().push(request(FINAL));
        }
        configure(&mut spec, &[parameter(Target::Power, "solid")]).unwrap();
        assert_eq!(spec, expected);
        configure(&mut spec, &[parameter(Target::Power, "solid")]).unwrap();
        assert_eq!(spec, expected, "recovery selects identical report bytes");
    }
    for names in [
        vec![PEAK],
        vec![PEAK_MATERIAL],
        vec![super::super::OUTPUT],
        vec!["temperature-volume-mean-adjoint"],
        vec![FINAL, FINAL],
        vec![FINAL, FINAL_MATERIAL],
    ] {
        let mut spec = base();
        spec.outputs
            .as_mut()
            .unwrap()
            .extend(names.into_iter().map(request));
        let before = spec.clone();
        assert!(configure(&mut spec, &[]).is_err());
        assert_eq!(
            spec, before,
            "a different objective is never silently replaced"
        );
    }
}

#[test]
fn g0_trajectory_coefficients_convert_delivered_power_and_log_htc_at_sampled_point() {
    let mut spec = base();
    spec.power.as_mut().unwrap()[0].duty = 0.37;
    let ThermalBoundaryCondition::Convection { coefficient, .. } = &mut spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition
    else {
        panic!("reference convection");
    };
    coefficient.value = 12.0;
    let controls = format!(
        r#"{POWER},"robin":[{{"region":"solid","dtemperature_dreference_k_per_k":0.6,"dtemperature_dlog_htc_k":-3}}]"#
    );
    let parameters = [
        parameter(Target::Power, "solid"),
        parameter(Target::ConvectionCoefficient, "solid"),
        parameter(Target::ConvectionTemperature, "solid"),
    ];
    assert_eq!(
        coefficients(
            &receipt(false, &controls),
            &spec,
            &parameters,
            301.0,
            "solid"
        )
        .unwrap(),
        [0.74, -0.25, 0.6]
    );
    let time = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap();
    time.power_schedules.push(TransientRegionPower {
        region: "solid".into(),
        source: "test workload".into(),
        steps: vec![TransientPowerStep {
            until: QtyAny::new(2.0, dims::TIME),
            watts: QtyAny::new(5.0, dims::POWER),
        }],
    });
    assert_eq!(
        coefficients(&receipt(false, ""), &spec, &parameters[..1], 301.0, "solid").unwrap(),
        [0.0],
        "scheduled watts replace the static coordinate entirely"
    );
}

#[test]
fn g0_history_controls_use_physical_reservoir_branch_bank_and_surface_identities() {
    let controls = r#", "convection":{"natural_reservoirs":[{"target":"natural-wall","dtemperature_dambient_temperature_k_per_k":0.3}],"air_inlets":[{"branch":"air-path","dtemperature_dinlet_temperature_k_per_k":0.4}],"fan_banks":[{"bank":"fans","status":"available","dtemperature_dspeed_ratio_k":-2}]},"radiation":{"surfaces":[{"name":"radiator","target":"different-geometric-target","dtemperature_dreservoir_temperature_k_per_k":0.2}]}"#;
    let parameters = [
        parameter(Target::NaturalConvectionAmbient, "natural-wall"),
        parameter(Target::AirInletTemperature, "air-path"),
        parameter(Target::FanSpeedRatio, "fans"),
        parameter(Target::RadiationReservoirTemperature, "radiator"),
    ];
    assert_eq!(
        coefficients(
            &receipt(false, controls),
            &base(),
            &parameters,
            301.0,
            "solid"
        )
        .unwrap(),
        [0.3, 0.4, -2.0, 0.2]
    );
    for bad in [
        controls.replace("available", "unsupported"),
        controls.replace("radiator", "another-surface"),
        controls.replace("air-path", "wrong-branch"),
        controls.replace("natural-wall", "wrong-target"),
    ] {
        assert!(coefficients(&receipt(false, &bad), &base(), &parameters, 301.0, "solid").is_err());
    }
    let duplicate = r#", "regional_power":[{"region":"solid","dtemperature_duniform_added_watts_k_per_w":2},{"region":"solid","dtemperature_duniform_added_watts_k_per_w":3}]"#;
    assert!(
        coefficients(
            &receipt(false, duplicate),
            &base(),
            &[parameter(Target::Power, "solid")],
            301.0,
            "solid"
        )
        .is_err()
    );
}

#[test]
fn g0_calibration_cannot_consume_another_time_grid_objective_or_storage_law() {
    let spec = base();
    let parameters = [parameter(Target::Power, "solid")];
    let good = receipt(false, POWER);
    // Keep mutations at the wire boundary used by retained calibration recovery.
    let source = format!(
        r#"{{"transient":{{"schema":"fs-cli-transient-conduction-v1","status":"completed","qoi_time":"final","final_time_s":2,"fine_steps":8,"temporal_error":{{"qoi":"temperature-max","region":"solid","fine_final_k":301}},"adjoint":{{"schema":"fs-cli-native-transient-adjoint-v1","output":"temperature-final-adjoint","value_k":301,"time_s":2,"state_index":8,"max_reverse_steps":8,"maximum_dual_relative_residual":0.000000001{POWER}}}}}}}"#
    );
    for bad in [
        source.replace(FINAL, PEAK),
        source.replace("\"time_s\":2", "\"time_s\":1"),
        source.replace("\"state_index\":8", "\"state_index\":7"),
        source.replace("\"value_k\":301", "\"value_k\":302"),
        source.replace("\"region\":\"solid\"", "\"region\":\"other\""),
        source.replace("0.000000001", "-1"),
        source.replace("\"fine_steps\":8", "\"fine_steps\":6"),
    ] {
        assert!(
            coefficients(&J::parse(&bad).unwrap(), &spec, &parameters, 301.0, "solid").is_err()
        );
    }
    assert!(coefficients(&good, &spec, &parameters, 300.0, "solid").is_err());
    assert!(coefficients(&receipt(true, POWER), &spec, &parameters, 301.0, "solid").is_err());
}

#[test]
fn g0_enthalpy_uses_its_own_producer_and_refuses_unprovided_coordinates() {
    let mut spec = base();
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
        .enthalpy = Some(TransientEnthalpy {
        materials: Vec::new(),
        max_newton_iterations: 20,
        max_backtracks: 12,
    });
    assert_eq!(
        coefficients(
            &receipt(true, POWER),
            &spec,
            &[parameter(Target::Power, "solid")],
            301.0,
            "solid"
        )
        .unwrap(),
        [2.0]
    );
    for target in [
        Target::HeatFlux,
        Target::ContactPressure,
        Target::ContactThickness,
        Target::ContactGap,
        Target::ContactTorque,
    ] {
        let p = parameter(target, "solid");
        let before = spec.clone();
        assert!(configure(&mut spec, std::slice::from_ref(&p)).is_err());
        assert_eq!(spec, before);
        let constant = UniformParameter { high: p.low, ..p };
        assert_eq!(
            coefficients(&receipt(true, ""), &spec, &[constant], 301.0, "solid").unwrap(),
            [0.0]
        );
    }
}

#[test]
fn g0_reference_transient_study_binds_existing_native_inputs() {
    let study = fs_project::uncertainty::UncertaintyStudy::parse(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-transient-uncertainty.fsim"
    )))
    .unwrap();
    let mut project = base();
    study.clone().bind(&project).unwrap();
    configure(&mut project, study.parameters()).unwrap();
    assert_eq!(project.outputs.unwrap().last().unwrap().name, FINAL);
}

#[test]
fn g0_bath_coefficients_bind_temperature_ownership_and_entire_history_for_both_storage_laws() {
    for enthalpy in [false, true] {
    let mut spec = base();
    let parameter = parameter(Target::FixedTemperature, "solid");
    assert!(configure(&mut spec, std::slice::from_ref(&parameter)).is_err(),
        "a convection reservoir is not a prescribed solid temperature");
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = ThermalBoundaryCondition::FixedTemperature {
        temperature: QtyAny::new(310.0, dims::TEMPERATURE),
    };
    if enthalpy {
        setup.transient.as_mut().unwrap().enthalpy = Some(TransientEnthalpy {
            materials: Vec::new(), max_newton_iterations: 20, max_backtracks: 12,
        });
    }
    configure(&mut spec, std::slice::from_ref(&parameter)).unwrap();
    let controls = r#", "prescribed_temperature":{"control_steps":8,"parameter_scope":"independent-uniform-bath-temperatures-held-throughout-selected-trajectory","boundaries":[{"target":"solid","temperature_k":310,"status":"available","mapped_vertices":3,"dtemperature_dprescribed_temperature_k_per_k":0.625}]}"#;
    assert_eq!(coefficients(&receipt(enthalpy, controls), &spec,
        std::slice::from_ref(&parameter), 301.0, "solid").unwrap(), [0.625]);
    for changed in [
        controls.replace("\"control_steps\":8", "\"control_steps\":0"),
        controls.replace("held-throughout-selected-trajectory", "single-endpoint"),
        controls.replace("\"target\":\"solid\"", "\"target\":\"other\""),
        controls.replace("\"temperature_k\":310", "\"temperature_k\":309"),
        controls.replace("\"status\":\"available\"", "\"status\":\"unsupported\",\"reason\":\"shared prescribed vertices\""),
        controls.replace("dtemperature_dprescribed_temperature_k_per_k", "unrelated_same_unit_control"),
        controls.replace("0.625", "null"),
    ] {
        assert!(coefficients(&receipt(enthalpy, &changed), &spec,
            std::slice::from_ref(&parameter), 301.0, "solid").is_err());
    }
    assert!(coefficients(&receipt(enthalpy, ""), &spec,
        std::slice::from_ref(&parameter), 301.0, "solid").is_err());
    let duplicate = controls.replace(
        "\"boundaries\":[{",
        "\"boundaries\":[{\"target\":\"solid\",\"temperature_k\":310,\"status\":\"available\",\"dtemperature_dprescribed_temperature_k_per_k\":0.25},{",
    );
    assert!(coefficients(&receipt(enthalpy, &duplicate), &spec,
        std::slice::from_ref(&parameter), 301.0, "solid").is_err());
    let constant = UniformParameter { high: parameter.low, ..parameter };
    assert_eq!(coefficients(&receipt(enthalpy, ""), &spec, &[constant], 301.0, "solid").unwrap(), [0.0]);
    }
}

#[test]
fn g0_outward_flux_coefficients_bind_original_flux_and_complete_history_for_both_storage_laws() {
    for enthalpy in [false, true] {
        let mut spec = base();
        let parameter = parameter(Target::HeatFlux, "solid");
        assert!(configure(&mut spec, std::slice::from_ref(&parameter)).is_err(),
            "a Robin boundary cannot supply a prescribed Neumann control");
        let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries[0].condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(-2.5, dims::HEAT_FLUX),
        };
        if enthalpy {
            setup.transient.as_mut().unwrap().enthalpy = Some(TransientEnthalpy {
                materials: Vec::new(), max_newton_iterations: 20, max_backtracks: 12,
            });
        }
        configure(&mut spec, std::slice::from_ref(&parameter)).unwrap();
        let row = r#"{"target":"solid","flux_w_m2":-2.5,"dtemperature_doutward_flux_k_per_w_m2":-0.125}"#;
        let controls = format!(r#", "prescribed_heat_flux":{{"control_steps":8,"parameter_scope":"independent-uniform-outward-heat-fluxes-held-throughout-selected-trajectory","boundaries":[{row}]}}"#);
        assert_eq!(coefficients(&receipt(enthalpy, &controls), &spec,
            std::slice::from_ref(&parameter), 301.0, "solid").unwrap(), [-0.125],
            "absolute outward-flux derivative retains its sign and timestep convention");
        for changed in [
            controls.replace("\"control_steps\":8", "\"control_steps\":7"),
            controls.replace("held-throughout-selected-trajectory", "single-endpoint"),
            controls.replace("\"target\":\"solid\"", "\"target\":\"other\""),
            controls.replace("\"flux_w_m2\":-2.5", "\"flux_w_m2\":2.5"),
            controls.replace("dtemperature_doutward_flux_k_per_w_m2", "unrelated_same_unit_control"),
            controls.replace("-0.125", "null"),
            controls.replace(row, &format!("{row},{row}")),
        ] {
            assert!(coefficients(&receipt(enthalpy, &changed), &spec,
                std::slice::from_ref(&parameter), 301.0, "solid").is_err());
        }
        assert!(coefficients(&receipt(enthalpy, ""), &spec,
            std::slice::from_ref(&parameter), 301.0, "solid").is_err());
        let constant = UniformParameter { high: parameter.low, ..parameter };
        assert_eq!(coefficients(&receipt(enthalpy, ""), &spec, &[constant], 301.0, "solid")
            .unwrap(), [0.0]);
    }
}

#[test]
fn g0_enthalpy_convection_coefficients_bind_absolute_coordinates_and_complete_history() {
    let mut spec = base();
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(12.0, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(310.0, dims::TEMPERATURE),
    };
    setup.transient.as_mut().unwrap().enthalpy = Some(TransientEnthalpy {
        materials: Vec::new(), max_newton_iterations: 20, max_backtracks: 12,
    });
    let parameters = [
        parameter(Target::ConvectionCoefficient, "solid"),
        parameter(Target::ConvectionTemperature, "solid"),
    ];
    let original_physics = spec.cooling.clone();
    configure(&mut spec, &parameters).unwrap();
    assert_eq!(spec.cooling, original_physics);
    let row = r#"{"target":"solid","coefficient_w_m2_k":12,"reference_temperature_k":310,"dtemperature_dreference_temperature_k_per_k":0.625,"dtemperature_dcoefficient_k_per_w_m2_k":-0.25,"dtemperature_dlog_coefficient_k":-3}"#;
    let controls = format!(r#","prescribed_convection":{{"control_steps":8,"parameter_scope":"independent-original-reference-and-positive-coefficient-held-throughout-selected-trajectory","boundaries":[{row}]}}"#);
    assert_eq!(coefficients(&receipt(true, &controls), &spec, &parameters, 301.0, "solid").unwrap(),
        [-0.25, 0.625], "absolute coefficient partial is not divided by h a second time");
    for changed in [
        controls.replace("\"control_steps\":8", "\"control_steps\":7"),
        controls.replace("held-throughout-selected-trajectory", "single-endpoint"),
        controls.replace("\"target\":\"solid\"", "\"target\":\"other\""),
        controls.replace("\"coefficient_w_m2_k\":12", "\"coefficient_w_m2_k\":13"),
        controls.replace("\"reference_temperature_k\":310", "\"reference_temperature_k\":309"),
        controls.replace("dtemperature_dcoefficient_k_per_w_m2_k", "unrelated_coefficient"),
        controls.replace("dtemperature_dreference_temperature_k_per_k", "unrelated_temperature"),
        controls.replace("-0.25", "null"),
        controls.replace("0.625", "null"),
        controls.replace(row, &format!("{row},{row}")),
    ] {
        assert!(coefficients(&receipt(true, &changed), &spec, &parameters, 301.0, "solid").is_err());
    }
    assert!(coefficients(&receipt(true, ""), &spec, &parameters, 301.0, "solid").is_err());
    let constants = parameters.map(|p| UniformParameter { high: p.low, ..p });
    assert_eq!(coefficients(&receipt(true, ""), &spec, &constants, 301.0, "solid").unwrap(), [0.0; 2]);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(310.0, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    assert!(coefficients(&receipt(true, &controls), &spec,
        &[parameter(Target::ConvectionTemperature, "solid")], 301.0, "solid").is_err(),
        "a dynamic natural reservoir is not an independently prescribed reference");
}

