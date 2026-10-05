//! Native .fsim pipeline with a synthetic, explicitly nonlinear material card.
//! Control differences use the selected vertex, not a relocated maximum.
#[path = "../src/json_read.rs"]
mod json_read;
use json_read::JsonValue;
use std::path::{Path, PathBuf};
use fs_cli::cards::{CardPackKind, CardPackSet, RawCardPack};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};

#[path = "native_nonlinear_adjoint/material_controls.rs"]
mod material_controls;

fn scratch() -> PathBuf {
    for ordinal in 0..10000 {
        let path = std::env::temp_dir().join(format!("fs-native-nonlinear-adjoint-{}-{ordinal}", std::process::id()));
        match std::fs::create_dir(&path) {
            Ok(()) => return path,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {},
            Err(e) => panic!("test directory: {e}"),
        }
    }
    panic!("test directory namespace exhausted");
}
fn command(args: &[&str]) -> JsonValue {
    let output = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{args:?}\n{}\n{}", output.stdout, output.stderr);
    JsonValue::parse(&output.stdout).unwrap()
}
fn artifact(ledger: &Path, identity: &str) -> JsonValue {
    let hash = fs_blake3::ContentHash::from_hex(identity).unwrap();
    let bytes = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap()
        .get_artifact(&hash).unwrap().unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

fn nonlinear_card() -> (Vec<u8>, String, String) {
    let source = fs_blake3::hash_bytes(b"synthetic k(T) native-adjoint fixture: 250K=1W/mK,450K=11W/mK");
    let provenance = || Provenance { source: "synthetic nonlinear conductivity; not experimental data".into(),
        license: "CC-BY-4.0; redistribution permitted".into(), artifact: Some(source) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "native-adjoint synthetic solid".into(), method: "declared affine k(T)".into(),
        artifact: source, caveats: "numerical verification only; not a validated material".into(),
        provenance: provenance(),
    }).unwrap();
    let dims = fs_conduction::CONDUCTIVITY_DIMS;
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("thermal-conductivity", dims),
        value: PropertyValue::Curve { abscissa: "T".into(),
            abscissa_dims: fs_project::spec::dims::TEMPERATURE,
            knots: vec![(250.0,1.0),(450.0,11.0)], dims },
        validity: fs_evidence::ValidityDomain::unconstrained().with("T",250.0,450.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![observation], provenance: provenance(),
    }).unwrap();
    let pack = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "nonlinear-adjoint".into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, NormalizedPack::new("nonlinear-adjoint", "synthetic-v1", source,
        "CC-BY-4.0; redistribution permitted", claims, Vec::new(), Vec::new()).unwrap()).unwrap();
    let bytes = pack.to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack { kind: CardPackKind::Material,
        source: "synthetic-nonlinear".into(), bytes: bytes.clone(), expect: None }]).unwrap();
    let card = &admitted.materials()[0];
    (bytes, card.card().to_hex(), card.identity().to_string())
}

struct Fixture { dir: PathBuf, project: fs_project::ProjectSpec }
impl Fixture {
    fn new() -> Self {
        let dir = scratch();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        let mut project = fs_project::parse_sexpr_migrating(
            &std::fs::read_to_string(reference.join("cooling-reference.fsim")).unwrap()).unwrap().decoded.spec;
        std::fs::copy(reference.join("plate.stl"), dir.join("plate.stl")).unwrap();
        let (bytes, card, state) = nonlinear_card();
        std::fs::write(dir.join("nonlinear.fsmcdpk"), bytes).unwrap();
        let material = &mut project.materials.as_mut().unwrap()[0];
        material.card = card; material.state = state; material.claim = None;
        material.temp_lo.value = 250.0; material.temp_hi.value = 450.0;
        material.source = "synthetic nonlinear adjoint regression".into();
        material.conductivity_tolerance = None;
        let power = &mut project.power.as_mut().unwrap()[0];
        power.watts.value = 1000.0; power.duty = 0.37;
        let envelope = project.envelope.as_mut().unwrap();
        envelope.ambient_lo.value = 290.0; envelope.ambient_hi.value = 310.0;
        let fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. }
            = &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
            else { panic!("reference convection law"); };
        reference_temperature.value = 300.0;
        project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
        Self { dir, project }
    }
    fn solve(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (JsonValue, JsonValue) {
        let path = self.dir.join(format!("project-{ordinal}.fsim"));
        std::fs::write(&path, fs_project::print_sexpr(project).unwrap()).unwrap();
        let ledger = self.dir.join("results.db");
        command(&["--json","import",path.to_str().unwrap(),self.dir.join("plate.stl").to_str().unwrap(),
            ledger.to_str().unwrap(),"--unit","m","--max-hole-edges","0"]);
        let run = command(&["--json","solve",path.to_str().unwrap(),ledger.to_str().unwrap(),
            "--materials",self.dir.join("nonlinear.fsmcdpk").to_str().unwrap()]);
        let run_receipt = artifact(&ledger, run.str_field("run_receipt").unwrap());
        let stage = run_receipt.get("stages").unwrap().as_array().unwrap().iter()
            .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
        let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
        let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
        (receipt, field)
    }
}

#[test]
fn nonlinear_material_native_adjoint_matches_physical_controls_on_the_retained_vertex() {
    let fixture = Fixture::new();
    let (baseline, base_field) = fixture.solve(&fixture.project, 0);
    let mut requested = fixture.project.clone();
    requested.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
        name: "temperature-max-adjoint".into(), kind: "report".into(), region: None,
    });
    let (receipt, field) = fixture.solve(&requested, 1);
    assert_eq!(baseline.get("energy"), receipt.get("energy"));
    assert_eq!(base_field.get("temperature"), field.get("temperature"));
    let adjoint = receipt.get("nominal_adjoint").unwrap();
    assert_eq!(adjoint.str_field("mode"), Some("nonlinear-solid-full-material-feedback"));
    assert_eq!(adjoint.str_field("authority"), Some("Estimated"));
    assert!(adjoint.f64_field("true_relative_residual").unwrap() < 1e-12);
    let vertex = adjoint.f64_field("selected_vertex").unwrap() as usize;
    let rows = adjoint.get("parameters").unwrap().as_array().unwrap();
    for (case, target, delta) in [(0,"power",1.0),(1,"convection-coefficient",0.001),(2,"convection-temperature",0.01)] {
        let gradient = rows.iter().find(|row| row.str_field("target") == Some(target)).unwrap()
            .f64_field("derivative").unwrap();
        let mut values = Vec::new();
        for (side, sign) in [-1.0,1.0].into_iter().enumerate() {
            let mut spec = fixture.project.clone();
            if case == 0 { spec.power.as_mut().unwrap()[0].watts.value += sign*delta; }
            else {
                let fs_project::ThermalBoundaryCondition::Convection { coefficient, reference_temperature }
                    = &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                    else { panic!("reference convection law"); };
                if case == 1 { coefficient.value += sign*delta; }
                else { reference_temperature.value += sign*delta; }
            }
            let (_, field) = fixture.solve(&spec, 2+2*case+side);
            values.push(field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap());
        }
        let expected = (values[1]-values[0])/(2.0*delta);
        assert!((gradient-expected).abs() < 1e-4*expected.abs().max(1e-4),
            "{target}: full material adjoint {gradient:e}, native resolve {expected:e}");
    }
}
