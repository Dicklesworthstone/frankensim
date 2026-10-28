use super::*;
use fs_project::{TransientEnthalpy, TransientPowerStep, TransientRegionPower, spec::dims};
use fs_qty::QtyAny;

fn base() -> ProjectSpec {
    fs_project::parse_sexpr(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-transient.fsim"
    )))
    .unwrap()
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
        Target::ConvectionCoefficient,
        Target::ConvectionTemperature,
        Target::HeatFlux,
        Target::FixedTemperature,
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
