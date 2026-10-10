//! Actual import, native FEM field and retained continuum mean enclosure.
#[path = "../src/json_read.rs"]
mod json_read;
use fs_cli::cards::{CardPackKind, CardPackSet, RawCardPack};
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
        Self { dir, project }
    }
    fn import(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (PathBuf, PathBuf) {
        let source = self.dir.join(format!("mean-{ordinal}.fsim"));
        let ledger = self.dir.join("runs.db");
        std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
        let mesh =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project/plate.stl");
        command(&[
            "--json",
            "import",
            source.to_str().unwrap(),
            mesh.to_str().unwrap(),
            ledger.to_str().unwrap(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]);
        (source, ledger)
    }
    fn solve(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (J, J, J) {
        let (source, ledger) = self.import(project, ordinal);
        let result = command(&[
            "--json",
            "solve",
            source.to_str().unwrap(),
            ledger.to_str().unwrap(),
            "--materials",
            self.dir.join("constant.fsmcdpk").to_str().unwrap(),
        ]);
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
        let result = fs_cli::run(vec![
            "--json".into(),
            "solve".into(),
            source.to_string_lossy().into_owned(),
            ledger.to_string_lossy().into_owned(),
            "--materials".into(),
            self.dir
                .join("constant.fsmcdpk")
                .to_string_lossy()
                .into_owned(),
        ]);
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
        initial_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
        horizon: fs_qty::QtyAny::new(1.0, fs_project::spec::dims::TIME),
        max_step: fs_qty::QtyAny::new(0.5, fs_project::spec::dims::TIME),
        max_steps: 12,
        energy_tolerance: fs_qty::QtyAny::new(1e-7, fs_qty::Dims([2, 1, -2, 0, 0, 0])),
        power_schedules: Vec::new(),
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
