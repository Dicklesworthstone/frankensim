//! Actual import, native FEM field and retained continuum mean enclosure.
#[path = "../src/json_read.rs"]
mod json_read;
use fs_cli::{CardPackKind, CardPackSet, RawCardPack};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use json_read::JsonValue as J;
use std::path::{Path, PathBuf};

const OUTPUT: &str = "temperature-volume-mean-bound";

fn command(args: &[&str]) -> J {
    let result = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(
        result.exit_code,
        fs_cli::exit::SUCCESS,
        "{args:?}\n{}\n{}",
        result.stdout,
        result.stderr
    );
    J::parse(&result.stdout).unwrap()
}
fn artifact(ledger: &Path, hash: &str) -> J {
    let hash = fs_blake3::ContentHash::from_hex(hash).unwrap();
    let bytes = fs_ledger::Ledger::open(ledger.to_str().unwrap())
        .unwrap()
        .get_artifact(&hash)
        .unwrap()
        .unwrap();
    J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}
fn constant_pack(bounded: bool) -> (Vec<u8>, String, String) {
    let source = fs_blake3::hash_bytes(b"declared-constant-mean-bound-fixture");
    let provenance = || Provenance {
        source: "explicit constant numerical fixture; not physical validation".into(),
        license: "CC-BY-4.0; redistribution permitted".into(),
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "numerical-fixture".into(),
            method: "declared constant PDE".into(),
            artifact: source,
            caveats: "synthetic; not experimental evidence".into(),
            provenance: provenance(),
        })
        .unwrap();
    let domain = fs_evidence::ValidityDomain::unconstrained();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal-conductivity", fs_conduction::CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 10.0,
                dims: fs_conduction::CONDUCTIVITY_DIMS,
            },
            validity: if bounded {
                domain.with("T", 200.0, 450.0)
            } else {
                domain
            },
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let pack = NormalizedPack::new(
        "global-constant-fixture",
        "declared-v1",
        source,
        "CC-BY-4.0; redistribution permitted",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let bytes = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "constant-solid".into(),
            phase: "solid".into(),
            process: "numerical".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
    .to_bytes();
    let set = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Material,
        source: "declared-constant-fixture".into(),
        bytes: bytes.clone(),
        expect: None,
    }])
    .unwrap();
    (
        bytes,
        set.materials()[0].card().to_hex(),
        set.materials()[0].identity().to_string(),
    )
}

struct Fixture {
    dir: PathBuf,
    project: fs_project::ProjectSpec,
    meshes: Vec<PathBuf>,
    interface_pack: Option<PathBuf>,
}
impl Fixture {
    fn new(bounded: bool) -> Self {
        let mut ordinal = 0;
        let dir = loop {
            let path = std::env::temp_dir().join(format!(
                "fs-native-mean-bound-{}-{ordinal}",
                std::process::id()
            ));
            match std::fs::create_dir(&path) {
                Ok(()) => break path,
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => ordinal += 1,
                Err(e) => panic!("test directory: {e}"),
            }
        };
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let mut project = fs_project::parse_sexpr_migrating(
            &std::fs::read_to_string(root.join("data/reference-project/cooling-reference.fsim"))
                .unwrap(),
        )
        .unwrap()
        .decoded
        .spec;
        let (bytes, card, state) = constant_pack(bounded);
        std::fs::write(dir.join("constant.fsmcdpk"), bytes).unwrap();
        for binding in project.materials.as_mut().unwrap() {
            binding.card = card.clone();
            binding.state = state.clone();
            binding.claim = None;
            binding.conductivity_tolerance = None;
        }
        project.solver.as_mut().unwrap().tolerance_rel = 1e-11;
        Self {
            dir,
            project,
            meshes: vec![root.join("data/reference-project/plate.stl")],
            interface_pack: None,
        }
    }
    fn import(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (PathBuf, PathBuf) {
        let source = self.dir.join(format!("mean-{ordinal}.fsim"));
        let ledger = self.dir.join("runs.db");
        std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
        let mut args = vec!["--json", "import", source.to_str().unwrap()];
        args.extend(self.meshes.iter().map(|path| path.to_str().unwrap()));
        args.extend([
            ledger.to_str().unwrap(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]);
        command(&args);
        (source, ledger)
    }
    fn solve(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (J, J, J) {
        let (source, ledger) = self.import(project, ordinal);
        let material = self.dir.join("constant.fsmcdpk");
        let mut args = vec![
            "--json",
            "solve",
            source.to_str().unwrap(),
            ledger.to_str().unwrap(),
            "--materials",
            material.to_str().unwrap(),
        ];
        if let Some(interface) = &self.interface_pack {
            args.extend(["--interfaces", interface.to_str().unwrap()]);
        }
        let result = command(&args);
        let run = artifact(&ledger, result.str_field("run_receipt").unwrap());
        let stage = run
            .get("stages")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s.str_field("stage") == Some("conduction"))
            .unwrap();
        let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
        let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
        (result, receipt, field)
    }
    fn refused(&self, project: &fs_project::ProjectSpec, ordinal: usize, reason: &str) {
        let (source, ledger) = self.import(project, ordinal);
        let mut args = vec![
            "--json".into(),
            "solve".into(),
            source.to_string_lossy().into_owned(),
            ledger.to_string_lossy().into_owned(),
            "--materials".into(),
            self.dir
                .join("constant.fsmcdpk")
                .to_string_lossy()
                .into_owned(),
        ];
        if let Some(interface) = &self.interface_pack {
            args.extend([
                "--interfaces".into(),
                interface.to_string_lossy().into_owned(),
            ]);
        }
        let result = fs_cli::run(args);
        assert_ne!(
            result.exit_code,
            fs_cli::exit::SUCCESS,
            "a refused bound cannot be a success"
        );
        assert!(
            format!("{}{}", result.stdout, result.stderr).contains(reason),
            "{}{}",
            result.stdout,
            result.stderr
        );
    }
}
fn request(project: &mut fs_project::ProjectSpec) {
    project
        .outputs
        .as_mut()
        .unwrap()
        .push(fs_project::spec::OutputRequest {
            name: OUTPUT.into(),
            kind: "report".into(),
            region: None,
        });
}

#[cfg(feature = "thermal-verification")]
fn contact_pack(resistance: f64, bounded: bool) -> (Vec<u8>, String) {
    use fs_matdb::{NormalizedInterfacePack, SurfaceSpec, SystemContext};
    let source = fs_blake3::hash_bytes(b"declared-constant-contact-mean-fixture");
    let provenance = || Provenance {
        source: "explicit constant contact PDE; not experimental data".into(),
        license: "CC-BY-4.0; redistribution permitted",
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "numerical-contact-fixture".into(),
            method: "declared finite-resistance transmission law".into(),
            artifact: source,
            caveats: "synthetic; not experimental evidence".into(),
            provenance: provenance(),
        })
        .unwrap();
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let domain = fs_evidence::ValidityDomain::unconstrained();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                dims,
            ),
            value: PropertyValue::Scalar {
                value: resistance,
                dims,
            },
            validity: if bounded {
                domain.with("T", 200.0, 450.0)
            } else {
                domain
            },
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let pack = NormalizedPack::new(
        "constant-contact-mean",
        "declared-v1",
        source,
        "CC-BY-4.0; redistribution permitted",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let state = |chemistry: &str| MaterialStateId {
        chemistry: chemistry.into(),
        phase: "solid".into(),
        process: "numerical".into(),
        revision: 0,
    };
    let bytes = NormalizedInterfacePack::new(
        SurfaceSpec {
            material: state("cold-body"),
            texture_frame: "normal-plus-x".into(),
        },
        SurfaceSpec {
            material: state("hot-body"),
            texture_frame: "normal-minus-x".into(),
        },
        SystemContext {
            medium: "dry-contact".into(),
            third_body: None,
            environment: "fixture-air".into(),
            history: "unaged".into(),
        },
        pack,
    )
    .unwrap()
    .to_bytes();
    let cards = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Interface,
        source: "declared-contact-mean-fixture".into(),
        bytes: bytes.clone(),
        expect: None,
    }])
    .unwrap();
    (bytes, cards.interfaces()[0].card().to_hex())
}

#[cfg(feature = "thermal-verification")]
fn contact_fixture(resistance: f64, bounded: bool) -> Fixture {
    use fs_project::{
        EntityDecl, GeometryAssignment, HalfSpaceSide, MeshSelector, ThermalBoundary,
        ThermalBoundaryCondition as B, spec::dims,
    };
    use fs_qty::QtyAny;
    let mut fixture = Fixture::new(false);
    let material = fixture.project.materials.as_ref().unwrap()[0].clone();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut project = fs_project::parse_sexpr_migrating(
        &std::fs::read_to_string(root.join("examples/contact-pair/contact-pair.fsim")).unwrap(),
    )
    .unwrap()
    .decoded
    .spec;
    for binding in project.materials.as_mut().unwrap() {
        binding.card = material.card.clone();
        binding.state = material.state.clone();
        binding.claim = None;
        binding.conductivity_tolerance = None;
    }
    let (bytes, card) = contact_pack(resistance, bounded);
    let path = fixture.dir.join("contact.fsintpk");
    std::fs::write(&path, bytes).unwrap();
    let binding = &mut project.interface_cards.as_mut().unwrap()[0];
    binding.card = card;
    binding.claim = None;
    for assignment in project.assignments.as_mut().unwrap() {
        assignment.allow_overlap = true;
    }
    // Unit cubes on [0,1] and [1,2], k=10, inlet heat flux q=20:
    // T_cold=300+2x, T_hot=300+2x+20R. Contact traces differ by 20R.
    // The HOT component mean is 303+20R, not the whole-domain average.
    for (name, artifact, offset, side) in [
        ("cold-anchor", "cold-body", 0.0, HalfSpaceSide::AtMost),
        ("hot-source", "hot-body", 2.0, HalfSpaceSide::AtLeast),
    ] {
        project
            .assembly
            .as_mut()
            .unwrap()
            .push(EntityDecl::Surface {
                parent: "enclosure".into(),
                name: name.into(),
                display: name.into(),
                expect_id: None,
            });
        project
            .assignments
            .as_mut()
            .unwrap()
            .push(GeometryAssignment {
                artifact: artifact.into(),
                target: name.into(),
                length_unit: "m".into(),
                allow_overlap: true,
                selector: MeshSelector::HalfSpace {
                    normal: [1.0, 0.0, 0.0],
                    offset,
                    side,
                    tolerance: 0.0,
                },
            });
    }
    project.power.as_mut().unwrap()[0].watts.value = 0.0;
    project.solver.as_mut().unwrap().tolerance_rel = 1e-11;
    let setup = project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup.adiabatic_remainder = true;
    setup.boundaries = vec![
        ThermalBoundary {
            target: "cold-anchor".into(),
            condition: B::FixedTemperature {
                temperature: QtyAny::new(300.0, dims::TEMPERATURE),
            },
        },
        ThermalBoundary {
            target: "hot-source".into(),
            condition: B::HeatFlux {
                outward_flux: QtyAny::new(-20.0, dims::HEAT_FLUX),
            },
        },
    ];
    fixture.project = project;
    fixture.meshes = vec![
        root.join("examples/contact-pair/cold-body.stl"),
        root.join("examples/contact-pair/hot-body.stl"),
    ];
    fixture.interface_pack = Some(path);
    fixture
}

#[cfg(not(feature = "thermal-verification"))]
#[test]
fn disabled_feature_refuses_the_named_mean_request() {
    let f = Fixture::new(false);
    let mut p = f.project.clone();
    request(&mut p);
    f.refused(&p, 0, "cli-solve-mean-bound-feature");
}

#[cfg(feature = "thermal-verification")]
fn affine_project(project: &mut fs_project::ProjectSpec) {
    use fs_project::{
        EntityDecl, GeometryAssignment, HalfSpaceSide, MeshSelector, ThermalBoundary,
        ThermalBoundaryCondition as B, spec::dims,
    };
    use fs_qty::QtyAny;
    // Exact unit right tetrahedron: T = 300 + 2 x. The x=0 face is
    // prescribed, the sloping face injects 20/sqrt(3) W/m², y/z faces insulate.
    // Its volume mean is 300.5 K; the maximum is 302 K.
    for row in project.assignments.as_mut().unwrap() {
        row.allow_overlap = true;
    }
    for (name, normal, offset, side) in [
        ("cold", [1.0, 0.0, 0.0], 0.0, HalfSpaceSide::AtMost),
        ("heated", [1.0, 1.0, 1.0], 1.0, HalfSpaceSide::AtLeast),
    ] {
        project
            .assembly
            .as_mut()
            .unwrap()
            .push(EntityDecl::Surface {
                parent: "enclosure".into(),
                name: name.into(),
                display: name.into(),
                expect_id: None,
            });
        project
            .assignments
            .as_mut()
            .unwrap()
            .push(GeometryAssignment {
                artifact: "enclosure".into(),
                target: name.into(),
                length_unit: "m".into(),
                allow_overlap: true,
                selector: MeshSelector::HalfSpace {
                    normal,
                    offset,
                    side,
                    tolerance: 0.0,
                },
            });
    }
    project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let setup = project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup.adiabatic_remainder = true;
    setup.boundaries = vec![
        ThermalBoundary {
            target: "cold".into(),
            condition: B::FixedTemperature {
                temperature: QtyAny::new(300.0, dims::TEMPERATURE),
            },
        },
        ThermalBoundary {
            target: "heated".into(),
            condition: B::HeatFlux {
                outward_flux: QtyAny::new(-20.0 / 3.0_f64.sqrt(), dims::HEAT_FLUX),
            },
        },
    ];
}

#[cfg(feature = "thermal-verification")]
#[test]
fn g1_native_affine_mean_bound_preserves_field_and_reaches_report_and_package() {
    let f = Fixture::new(false);
    let mut p = f.project.clone();
    affine_project(&mut p);
    let (_, baseline, original) = f.solve(&p, 0);
    request(&mut p);
    let (result, receipt, field) = f.solve(&p, 1);
    assert_eq!(original.get("temperature"), field.get("temperature"));
    assert_eq!(baseline.get("energy"), receipt.get("energy"));
    let bound = receipt
        .get("volume_mean_bound")
        .expect("real numerical producer retained");
    assert_eq!(
        bound.str_field("functional"),
        Some("region-volume-mean-temperature")
    );
    assert_eq!(bound.str_field("authority"), Some("Verified"));
    let interval = bound.get("enclosure_k").unwrap();
    assert!(
        interval.f64_field("lower").unwrap() <= 300.5
            && interval.f64_field("upper").unwrap() >= 300.5,
        "analytic affine mean must be enclosed: {interval:?}"
    );
    assert!((bound.f64_field("value_k").unwrap() - 300.5).abs() < 1e-7);
    assert!(bound.f64_field("error_upper_k").unwrap() < 1e-4);
    assert!(
        (receipt
            .get("temperature")
            .unwrap()
            .f64_field("max")
            .unwrap()
            - 302.0)
            .abs()
            < 1e-7,
        "volume mean is not the maximum"
    );
    let row = &bound.get("materials").unwrap().as_array().unwrap()[0];
    assert!(row.str_field("receipt_bytes_hex").unwrap().len() > 100);
    let (_, replay, _) = f.solve(&p, 2);
    assert_eq!(bound, replay.get("volume_mean_bound").unwrap());
    let ledger = f.dir.join("runs.db");
    let run = result.str_field("run").unwrap();
    let exported = command(&["--json", "report", run, ledger.to_str().unwrap()]);
    let json = std::fs::read_to_string(exported.str_field("report_json").unwrap()).unwrap();
    let exported_report = J::parse(&json).unwrap();
    let mean = exported_report
        .get("qois")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|qoi| qoi.str_field("name") == Some(OUTPUT))
        .expect("separate mean QoI");
    assert!((mean.f64_field("value").unwrap() - 300.5).abs() < 1e-7);
    assert_eq!(mean.str_field("color"), Some("Verified"));
    assert_eq!(
        mean.get("total_budget"),
        Some(&J::Null),
        "physical uncertainty is not certified"
    );
    let html = std::fs::read_to_string(exported.str_field("report_html").unwrap()).unwrap();
    assert!(html.contains(OUTPUT));
    let exported = command(&["--json", "package", run, ledger.to_str().unwrap()]);
    let text = std::fs::read_to_string(exported.str_field("package").unwrap()).unwrap();
    assert!(text.contains(OUTPUT));
    let package = fs_package::EvidencePackage::from_json(&text).unwrap();
    assert!(fs_checker::check(&package).passed());
}

#[cfg(feature = "thermal-verification")]
#[test]
fn g1_native_contact_mean_preserves_traces_and_encloses_the_resistance_dependent_component() {
    for resistance in [0.05, 0.2] {
        let fixture = contact_fixture(resistance, false);
        let mut project = fixture.project.clone();
        let (_, baseline, original) = fixture.solve(&project, 0);
        request(&mut project);
        let (run, receipt, field) = fixture.solve(&project, 1);
        assert_eq!(original.get("temperature"), field.get("temperature"));
        for key in ["energy", "interfaces"] {
            assert_eq!(
                baseline.get(key),
                receipt.get(key),
                "original {key} changed"
            );
        }
        let bound = receipt.get("volume_mean_bound").unwrap();
        let expected = 303.0 + 20.0 * resistance;
        assert_eq!(bound.str_field("region"), Some("hot"));
        assert_eq!(bound.str_field("authority"), Some("Verified"));
        let interval = bound.get("enclosure_k").unwrap();
        assert!(
            interval.f64_field("lower").unwrap() <= expected
                && interval.f64_field("upper").unwrap() >= expected,
            "R={resistance}: {interval:?} must enclose analytic hot mean {expected}"
        );
        assert!((bound.f64_field("value_k").unwrap() - expected).abs() < 1e-7);
        assert!(bound.f64_field("error_upper_k").unwrap() < 1e-3);
        let contacts = bound.get("contacts").unwrap().as_array().unwrap();
        assert_eq!(contacts.len(), 1);
        let contact = &contacts[0];
        assert_eq!(contact.str_field("interface"), Some("cold-hot-joint"));
        assert_eq!(
            contact.str_field("card"),
            Some(contact_pack(resistance, false).1.as_str())
        );
        assert_eq!(contact.f64_field("resistance_m2_k_per_w"), Some(resistance));
        assert!(
            !contact
                .get("face_pairs")
                .unwrap()
                .as_array()
                .unwrap()
                .is_empty()
        );
        assert!(contact.str_field("receipt_bytes_hex").unwrap().len() > 100);
        assert!(
            bound
                .str_field("scope")
                .unwrap()
                .contains("resistance-weighted jumps")
        );
        let (_, replay, _) = fixture.solve(&project, 2);
        assert_eq!(bound, replay.get("volume_mean_bound").unwrap());

        let ledger = fixture.dir.join("runs.db");
        let exported = command(&[
            "--json",
            "report",
            run.str_field("run").unwrap(),
            ledger.to_str().unwrap(),
        ]);
        let report =
            J::parse(&std::fs::read_to_string(exported.str_field("report_json").unwrap()).unwrap())
                .unwrap();
        let mean = report
            .get("qois")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .find(|qoi| qoi.str_field("name") == Some(OUTPUT))
            .unwrap();
        assert!((mean.f64_field("value").unwrap() - expected).abs() < 1e-7);
        assert_eq!(mean.str_field("color"), Some("Verified"));
        assert_eq!(
            mean.get("total_budget"),
            Some(&J::Null),
            "contact uncertainty remains unknown"
        );
        let exported = command(&[
            "--json",
            "package",
            run.str_field("run").unwrap(),
            ledger.to_str().unwrap(),
        ]);
        let package = fs_package::EvidencePackage::from_json(
            &std::fs::read_to_string(exported.str_field("package").unwrap()).unwrap(),
        )
        .unwrap();
        assert!(fs_checker::check(&package).passed());
    }
}

#[cfg(feature = "thermal-verification")]
#[test]
fn g0_finite_contact_validity_does_not_become_a_continuum_temperature_proof() {
    let fixture = contact_fixture(0.05, true);
    let mut project = fixture.project.clone();
    fixture.solve(&project, 0);
    request(&mut project);
    // This original fixture has nonzero Neumann heating. A finite contact
    // card still cannot use the layered comparison for that unsupported model.
    fixture.refused(&project, 1, "original contact/material domains");
}

#[cfg(feature = "thermal-verification")]
#[test]
fn g0_bounded_scalar_and_feedback_models_cannot_publish_a_frozen_mean_bound() {
    let bounded = Fixture::new(true);
    let mut p = bounded.project.clone();
    request(&mut p);
    bounded.refused(&p, 10, "finite spans");
    let f = Fixture::new(false);
    let mut p = f.project.clone();
    request(&mut p);
    p.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = fs_project::ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: fs_qty::QtyAny::new(1.0, fs_project::spec::dims::LENGTH),
        ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    for (ordinal, fidelity) in [(11, "auto"), (12, "ladder"), (13, "adaptive")] {
        p.solver.as_mut().unwrap().fidelity = fidelity.into();
        f.refused(&p, ordinal, "steady fixed thermal boundary laws");
    }
    let mut p = f.project.clone();
    request(&mut p);
    p.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = fs_project::ThermalBoundaryCondition::AirflowConvection {
        branch: "air".into(),
        order: 0,
        inlet_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
        hydraulic_diameter: fs_qty::QtyAny::new(0.02, fs_project::spec::dims::LENGTH),
        flow_area: fs_qty::QtyAny::new(0.004, fs_project::spec::dims::AREA),
        channel_length: fs_qty::QtyAny::new(0.3, fs_project::spec::dims::LENGTH),
        correlation: "convection.gnielinski".into(),
    };
    f.refused(&p, 14, "steady fixed thermal boundary laws");
    let mut p = f.project.clone();
    request(&mut p);
    p.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient = Some(fs_project::spec::ConductionTransient {
        initial_temperature: fs_project::TransientInitialTemperature::Uniform(
            fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE)),
        horizon: fs_qty::QtyAny::new(1.0, fs_project::spec::dims::TIME),
        max_step: fs_qty::QtyAny::new(0.5, fs_project::spec::dims::TIME),
        max_steps: 12,
        energy_tolerance: fs_qty::QtyAny::new(1e-7, fs_qty::Dims([2, 1, -2, 0, 0, 0])),
        power_schedules: Vec::new(),
        enthalpy: None,
        capacities: vec![fs_project::spec::TransientRegionCapacity {
            region: "air".into(),
            volumetric_heat_capacity: fs_qty::QtyAny::new(
                2e6,
                fs_conduction::transient::VOLUMETRIC_HEAT_CAPACITY_DIMS,
            ),
            source: "declared fixture capacity".into(),
        }],
    });
    f.refused(&p, 15, "steady fixed thermal boundary laws");
}

#[cfg(feature = "thermal-verification")]
mod finite_contacts {
    use super::*;

    fn heated_stack(r: f64, finite_material: bool, finite_contact: bool) -> Fixture {
        let mut fixture = contact_fixture(r, finite_contact);
        let (bytes, card, state) = constant_pack(finite_material);
        std::fs::write(fixture.dir.join("constant.fsmcdpk"), bytes).unwrap();
        for binding in fixture.project.materials.as_mut().unwrap() {
            binding.card = card.clone(); binding.state = state.clone();
            binding.claim = None; binding.conductivity_tolerance = None;
            binding.temp_lo.value = 290.0; binding.temp_hi.value = 330.0;
        }
        // Only the original hot region is heated. The source-free cold region
        // still conducts through its own volume and the unchanged contact.
        let power = &mut fixture.project.power.as_mut().unwrap()[0];
        assert_eq!(power.region, "hot");
        power.watts.value = 20.0;
        let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries.truncate(1);
        setup.boundaries[0].condition = fs_project::ThermalBoundaryCondition::Convection {
            coefficient: fs_qty::QtyAny::new(20.0, fs_qty::Dims([0,1,-3,-1,0,0])),
            reference_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
        };
        fixture
    }

    #[test]
    fn finite_contact_and_material_sources_reach_native_bound_and_replay() {
        for (finite_material, finite_contact) in [(false,true),(true,false),(true,true)] {
            let fixture = heated_stack(0.2, finite_material, finite_contact);
            let mut p = fixture.project.clone();
            let (_, baseline, original) = fixture.solve(&p, 0);
            request(&mut p);
            let (_, receipt, field) = fixture.solve(&p, 1);
            assert_eq!(original.get("temperature"), field.get("temperature"));
            for key in ["energy", "interfaces"] { assert_eq!(baseline.get(key), receipt.get(key)); }
            let bound = receipt.get("volume_mean_bound").unwrap();
            let proof = bound.get("material_validity").unwrap();
            assert_eq!(proof.str_field("method"), Some("continuum-layered-contact-comparison"));
            let range = proof.get("continuum_temperature_range_k").unwrap();
            // Bound all temperatures by replacing each original P1 source with
            // the positive envelope. This is conservative, not the actual peak.
            assert!(range.f64_field("lower").unwrap() <= 300.0);
            let upper = range.f64_field("upper").unwrap();
            assert!(upper >= 310.0 && upper < 310.000001);
            assert_eq!(proof.get("retained_material_domain_k") == Some(&J::Null), !finite_material);
            assert_eq!(proof.get("retained_contact_domain_k") == Some(&J::Null), !finite_contact);
            let interval = bound.get("enclosure_k").unwrap();
            let exact_hot_mean = 303.0 + 20.0*0.2 + 2.0/3.0;
            assert!(interval.f64_field("lower").unwrap() <= exact_hot_mean);
            assert!(interval.f64_field("upper").unwrap() >= exact_hot_mean);
            let contact = &bound.get("contacts").unwrap().as_array().unwrap()[0];
            assert_eq!(contact.str_field("card"), Some(contact_pack(0.2, finite_contact).1.as_str()));
            assert_eq!(contact.f64_field("resistance_m2_k_per_w"), Some(0.2));
            assert_eq!(contact.get("temperature_domain_k") == Some(&J::Null), !finite_contact);
            assert!(contact.str_field("receipt_bytes_hex").unwrap().len() > 100);
            let (_, replay, _) = fixture.solve(&p, 2);
            assert_eq!(bound, replay.get("volume_mean_bound").unwrap());
        }
    }

    #[test]
    fn finite_contact_numerical_bound_flows_into_existing_report_and_package() {
        let fixture = heated_stack(0.05, true, true);
        let mut p = fixture.project.clone(); request(&mut p);
        let (run, receipt, _) = fixture.solve(&p, 0);
        let bound = receipt.get("volume_mean_bound").unwrap();
        let ledger = fixture.dir.join("runs.db"); let id = run.str_field("run").unwrap();
        let exported = command(&["--json", "report", id, ledger.to_str().unwrap()]);
        let report = J::parse(&std::fs::read_to_string(exported.str_field("report_json").unwrap()).unwrap()).unwrap();
        let mean = report.get("qois").unwrap().as_array().unwrap().iter()
            .find(|q| q.str_field("name") == Some(OUTPUT)).unwrap();
        assert_eq!(mean.f64_field("value"), bound.f64_field("value_k"));
        assert_eq!(mean.str_field("color"), Some("Verified"));
        assert_eq!(mean.get("total_budget"), Some(&J::Null), "physical uncertainty stays unknown");
        let exported = command(&["--json", "package", id, ledger.to_str().unwrap()]);
        let package = fs_package::EvidencePackage::from_json(&std::fs::read_to_string(
            exported.str_field("package").unwrap()).unwrap()).unwrap();
        assert!(fs_checker::check(&package).passed());
    }

    #[test]
    fn conservative_contact_range_is_not_replaced_by_an_admissible_nodal_field() {
        let fixture = heated_stack(7.25, false, true);
        let mut p = fixture.project.clone();
        for binding in p.materials.as_mut().unwrap() { binding.temp_hi.value = 500.0; }
        // Actual exact hot maximum = 449 K, whereas this source-envelope
        // comparison reaches 451 K, outside the contact card's 450 K domain.
        // Failure of this particular proof must not be called physical failure.
        fixture.solve(&p, 0);
        request(&mut p);
        fixture.refused(&p, 1, "original contact/material domains");
    }
}
