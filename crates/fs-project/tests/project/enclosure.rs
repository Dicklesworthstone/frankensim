//! G0: the native declaration preserves distinct radiation physics and identity.
use super::*;
use fs_project::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};

fn enclosure_project() -> ProjectSpec {
    let mut spec = radiation_project();
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.radiation = None;
    setup.boundaries[1].condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(7.0, fs_project::spec::dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: kelvin(300.0),
    };
    setup.enclosure_radiation = Some(ConductionEnclosure {
        spatial_model: None,
        ambient_temperature: None,
        surfaces: [("sink-finish", "sink-base"), ("cpu-finish", "cpu")]
            .into_iter()
            .map(|(name, target)| EnclosureRadiatingSurface {
                name: name.into(),
                target: target.into(),
                card: "ab".repeat(32),
                claim: Some("cd".repeat(32)),
                query_temperature: kelvin(320.0),
            })
            .collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "authored equal-area parallel-plate limiting model".into(),
        },
        closure_tolerance: 1e-12,
        reciprocity_relative_tolerance: 1e-12,
        max_iterations: 200,
        temperature_tolerance: kelvin(1e-6),
        heat_tolerance: watts(1e-6),
        relaxation: 0.5,
    });
    spec
}

fn declaration(spec: &mut ProjectSpec) -> &mut ConductionEnclosure {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .enclosure_radiation
        .as_mut()
        .unwrap()
}

#[test]
fn g0_enclosure_radiation_spatial_model_roundtrip_identity_and_refusal() {
    check_radiation_spatial_models(enclosure_project());
}

#[test]
fn g0_enclosure_roundtrips_without_inventing_external_reservoirs() {
    let mut spec = enclosure_project();
    for evidence in [
        EnclosureViewFactorEvidence::Analytic {
            geometry: "parallel-plates-limit".into(),
        },
        EnclosureViewFactorEvidence::ExternalQmc {
            seed: u64::MAX,
            samples: 9_007_199_254_740_995,
            generator: "external-vf/3".into(),
        },
    ] {
        declaration(&mut spec).evidence = evidence;
        assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        let text = print_sexpr(&spec).unwrap();
        assert!(text.contains(":enclosure-radiation (enclosure-radiation"));
        assert!(!text.contains("reservoir-temperature"));
        let decoded = parse_sexpr(&text).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        assert_eq!(decoded.spec, spec);
        assert_eq!(parse_json(&print_json(&spec).unwrap()).unwrap().spec, spec);
        assert_eq!(decoded.canonical, text);
    }
}

#[test]
fn g0_pointwise_enclosure_admits_vacuum_and_authored_heat_flux_without_fake_convection() {
    for outward_flux in [0.0, -3.0] {
        let mut spec = enclosure_project();
        let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.enclosure_radiation.as_mut().unwrap().spatial_model =
            Some(fs_project::RadiationSpatialModel::Pointwise);
        setup.transient = enthalpy_project().cooling.unwrap().conduction.unwrap().transient;
        for boundary in &mut setup.boundaries {
            if ["sink-base", "cpu"].contains(&boundary.target.as_str()) {
                boundary.condition = ThermalBoundaryCondition::HeatFlux {
                    outward_flux: QtyAny::new(outward_flux, fs_project::spec::dims::HEAT_FLUX),
                };
            }
        }
        assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        for model in [None, Some(fs_project::RadiationSpatialModel::PatchMean)] {
            declaration(&mut spec).spatial_model = model;
            assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        }
    }
}

#[test]
fn g0_enclosure_matrix_shape_source_and_physical_coordinates_are_required() {
    let mut spec = enclosure_project();
    let bad = declaration(&mut spec);
    bad.view_factors = vec![vec![0.0, 0.4]];
    bad.closure_tolerance = -1.0;
    bad.reciprocity_relative_tolerance = 1.0;
    bad.relaxation = 0.0;
    bad.max_iterations = 0;
    bad.temperature_tolerance = watts(1e-6);
    bad.heat_tolerance = watts(0.0);
    bad.surfaces[1].name = bad.surfaces[0].name.clone();
    bad.surfaces[1].target = bad.surfaces[0].target.clone();
    bad.surfaces[0].card = "unknown-card".into();
    bad.surfaces[0].query_temperature = kelvin(0.0);
    bad.evidence = EnclosureViewFactorEvidence::ExternalQmc {
        seed: 0,
        samples: 0,
        generator: String::new(),
    };
    let findings = spec.validate();
    for code in [
        "project-enclosure-radiation-matrix",
        "project-enclosure-radiation-closure",
        "project-enclosure-radiation-factor-tolerance",
        "project-enclosure-radiation-relaxation",
        "project-enclosure-radiation-iterations",
        "project-enclosure-radiation-dims",
        "project-enclosure-radiation-tolerance",
        "project-enclosure-radiation-identity",
        "project-enclosure-radiation-card",
        "project-enclosure-radiation-temperature",
        "project-enclosure-radiation-evidence",
    ] {
        assert!(
            findings.iter().any(|row| row.code == code),
            "missing {code}: {findings:?}"
        );
    }
    let mut spec = enclosure_project();
    declaration(&mut spec).view_factors[0][0] = f64::NAN;
    assert!(
        spec.validate()
            .iter()
            .any(|row| row.code == "project-enclosure-radiation-factor")
    );
}

#[test]
fn g0_open_exchange_roundtrips_its_physical_deficits_and_explicit_black_bath() {
    for count in [1, 2] {
        let mut spec = enclosure_project();
        let enclosure = declaration(&mut spec);
        enclosure.ambient_temperature = Some(kelvin(3.0));
        enclosure.surfaces.truncate(count);
        enclosure.view_factors = if count == 1 { vec![vec![0.0]] }
            else { vec![vec![0.0, 0.3], vec![0.3, 0.0]] };
        assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        let text = print_sexpr(&spec).unwrap();
        assert!(text.contains(":ambient-temperature"));
        assert_eq!(parse_sexpr(&text).unwrap().spec, spec);
        assert_eq!(parse_json(&print_json(&spec).unwrap()).unwrap().spec, spec);
    }
}

#[test]
fn g0_open_exchange_refuses_overfull_rows_absent_escape_and_invalid_bath() {
    let mut spec = enclosure_project();
    declaration(&mut spec).ambient_temperature = Some(kelvin(300.0));
    assert!(spec.validate().iter().any(|row| row.code == "project-enclosure-radiation-no-escape"));
    let enclosure = declaration(&mut spec);
    enclosure.view_factors = vec![vec![0.0, 1.0 + 1e-10], vec![0.4, 0.0]];
    enclosure.closure_tolerance = 1e-2;
    assert!(spec.validate().iter().any(|row| row.code == "project-enclosure-radiation-closure"),
        "a numerical closure tolerance cannot authorize negative escape");
    let enclosure = declaration(&mut spec);
    enclosure.view_factors = vec![vec![0.0, 0.4], vec![0.4, 0.0]];
    enclosure.ambient_temperature = Some(watts(0.0));
    let findings = spec.validate();
    for code in ["project-enclosure-radiation-dims", "project-enclosure-radiation-temperature"] {
        assert!(findings.iter().any(|row| row.code == code), "{findings:?}");
    }
}

#[test]
fn g0_closed_exchange_is_not_an_absolute_temperature_anchor() {
    let mut spec = enclosure_project();
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    for boundary in &mut setup.boundaries {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, fs_project::spec::dims::HEAT_FLUX),
        };
    }
    let findings = spec.validate();
    assert!(
        findings
            .iter()
            .any(|row| row.code == "project-conduction-boundary-unanchored")
    );
    assert!(
        !findings
            .iter()
            .any(|row| row.code == "project-enclosure-radiation-boundary"),
        "zero flux is an admissible enclosure law; it simply supplies no steady anchor"
    );
}

#[test]
fn g0_enclosure_refuses_ambient_overlap() {
    let mut spec = enclosure_project();
    let ambient = radiation_project()
        .cooling
        .unwrap()
        .conduction
        .unwrap()
        .radiation;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.radiation = ambient;
    let findings = spec.validate();
    assert!(
        findings
            .iter()
            .any(|row| row.code == "project-enclosure-radiation-overlap")
    );
}

#[test]
fn g0_v15_migration_preserves_historical_radiation_and_refuses_new_payload() {
    let historical = radiation_project();
    let canonical = print_sexpr(&historical).unwrap();
    let old = canonical
        .replacen(
            &format!("(fsim-project :version {FSIM_VERSION}"),
            "(fsim-project :version 15",
            1,
        )
        .replacen(
            &format!("(versions :schema {FSIM_VERSION}"),
            "(versions :schema 15",
            1,
        );
    let migrated = parse_sexpr_migrating(&old).unwrap();
    assert_eq!(migrated.decoded.spec, historical);
    assert!(
        migrated
            .decoded
            .spec
            .cooling
            .as_ref()
            .unwrap()
            .conduction
            .as_ref()
            .unwrap()
            .enclosure_radiation
            .is_none()
    );
    let receipt = migrated.migration.unwrap();
    assert_eq!(receipt.rule.label(), "conduction-enclosure-radiation-v16");
    assert!(receipt.verifies(old.as_bytes(), canonical.as_bytes()));
    let enclosure = print_sexpr(&enclosure_project())
        .unwrap()
        .replacen(
            &format!("(fsim-project :version {FSIM_VERSION}"),
            "(fsim-project :version 15",
            1,
        )
        .replacen(
            &format!("(versions :schema {FSIM_VERSION}"),
            "(versions :schema 15",
            1,
        );
    assert_eq!(
        parse_sexpr_migrating(&enclosure).unwrap_err().code,
        "fsim-migration-payload"
    );
}

#[test]
fn g0_enclosure_wire_does_not_drop_reservoir_or_matrix_intent() {
    let canonical = print_sexpr(&enclosure_project()).unwrap();
    let unknown = canonical.replacen(
        ":enclosure-radiation (enclosure-radiation",
        ":enclosure-radiation (enclosure-radiation :reservoir-temperature 300K",
        1,
    );
    assert!(
        parse_sexpr_lenient(&unknown)
            .unwrap()
            .findings()
            .iter()
            .any(|row| row.code == "project-unknown-field")
    );
    let missing = canonical.replacen("(row 0.0 1.0)", "(row 1.0)", 1);
    assert_ne!(missing, canonical);
    assert!(
        parse_sexpr_lenient(&missing)
            .unwrap()
            .findings()
            .iter()
            .any(|row| row.code == "project-enclosure-radiation-shape")
    );
}

#[test]
fn g0_shipped_enclosure_example_is_an_admissible_native_project() {
    let project = parse_sexpr(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-enclosure.fsim"
    )))
    .unwrap();
    assert!(project.findings().is_empty(), "{:?}", project.findings());
    let enclosure = project
        .spec
        .cooling
        .unwrap()
        .conduction
        .unwrap()
        .enclosure_radiation
        .unwrap();
    assert_eq!(enclosure.surfaces.len(), 2);
    assert_eq!(enclosure.view_factors, vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
}

#[test]
fn g0_enclosure_flux_laws_roundtrip_with_an_external_anchor_or_transient_storage() {
    use fs_project::spec::dims;

    for value in [0.0, -12.0, 2.0] {
        let mut spec = enclosure_project();
        // The other radiating patch retains convection as an independent anchor.
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::HeatFlux {
                outward_flux: QtyAny::new(value, dims::HEAT_FLUX),
            };
        assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        let text = print_sexpr(&spec).unwrap();
        assert_eq!(parse_sexpr(&text).unwrap().spec, spec);
        assert_eq!(parse_json(&print_json(&spec).unwrap()).unwrap().spec, spec);
    }

    for storage in [transient_project(), enthalpy_project()] {
        let mut spec = enclosure_project();
        let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.transient = storage.cooling.unwrap().conduction.unwrap().transient;
        // Storage retains the absolute initial state without a fluid reservoir.
        for (index, boundary) in setup.boundaries.iter_mut().enumerate() {
            boundary.condition = ThermalBoundaryCondition::HeatFlux {
                outward_flux: QtyAny::new(if index == 0 { -12.0 } else { 2.0 }, dims::HEAT_FLUX),
            };
        }
        assert!(spec.validate().is_empty(), "{:?}", spec.validate());
        let text = print_sexpr(&spec).unwrap();
        assert_eq!(parse_sexpr(&text).unwrap().spec, spec);
        assert_eq!(parse_json(&print_json(&spec).unwrap()).unwrap().spec, spec);
    }
}

#[test]
fn g0_enclosure_flux_admission_keeps_quantity_and_boundary_ownership_checks() {
    use fs_project::spec::dims;

    for (value, dimension, code) in [
        (f64::NAN, dims::HEAT_FLUX, "project-conduction-flux-finite"),
        (f64::INFINITY, dims::HEAT_FLUX, "project-conduction-flux-finite"),
        (0.0, dims::POWER, "project-conduction-boundary-dims"),
    ] {
        let mut spec = enclosure_project();
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::HeatFlux { outward_flux: QtyAny::new(value, dimension) };
        let findings = spec.validate();
        assert!(findings.iter().any(|row| row.code == code), "{findings:?}");
    }

    let mut spec = enclosure_project();
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::FixedTemperature { temperature: kelvin(300.0) };
    assert!(spec.validate().iter().any(|row|
        row.code == "project-enclosure-radiation-boundary"),
        "the bath must be a separately declared face, not a replaced radiating law");
}
