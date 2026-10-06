//! Actual imported contact pair, card selection and complete native thermal solves.
#[path = "../src/json_read.rs"]
mod json_read;
use json_read::JsonValue as J;
use std::path::{Path, PathBuf};
use fs_cli::cards::{CardPackKind, CardPackSet, RawCardPack};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    NormalizedMaterialCardPack, NormalizedPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel};

#[path = "native_contact_adjoint/prescribed.rs"]
mod prescribed;

const OUTPUT: &str = "temperature-max-contact-adjoint";
const TARGET: &str = "contact-resistance-multiplier";

fn scratch() -> PathBuf {
    for ordinal in 0..10000 {
        let path = std::env::temp_dir().join(format!("fs-native-contact-adjoint-{}-{ordinal}", std::process::id()));
        match std::fs::create_dir(&path) {
            Ok(()) => return path,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {},
            Err(e) => panic!("test directory: {e}"),
        }
    }
    panic!("test directory namespace exhausted");
}
fn command(args: &[&str]) -> J {
    let output = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{args:?}\n{}\n{}", output.stdout, output.stderr);
    J::parse(&output.stdout).unwrap()
}
fn artifact(ledger: &Path, identity: &str) -> J {
    let hash = fs_blake3::ContentHash::from_hex(identity).unwrap();
    let bytes = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap().get_artifact(&hash).unwrap().unwrap();
    J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}
fn pack(name: &str, property: &str, dims: fs_qty::Dims,
    value: PropertyValue, interpolation: InterpolationPolicy) -> NormalizedPack {
    let source = fs_blake3::hash_bytes(name.as_bytes());
    let provenance = || Provenance { source: format!("synthetic {name}; not experimental data"),
        license: "CC-BY-4.0; redistribution permitted".into(), artifact: Some(source) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: name.into(), method: "explicit numerical contact-control fixture".into(),
        artifact: source, caveats: "numerical verification only".into(), provenance: provenance(),
    }).unwrap();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(property, dims), value,
        validity: fs_evidence::ValidityDomain::unconstrained().with("T",200.0,450.0),
        uncertainty: UncertaintyModel::Unstated, interpolation,
        observations: vec![observation], provenance: provenance(),
    }).unwrap();
    NormalizedPack::new(name,"synthetic-contact-control-v1",source,
        "CC-BY-4.0; redistribution permitted",claims,Vec::new(),Vec::new()).unwrap()
}
fn contact_pack(resistance: f64) -> (Vec<u8>, String) {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let claims = pack(&format!("contact-resistance-{resistance}"),
        fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims,
        PropertyValue::Scalar { value: resistance, dims }, InterpolationPolicy::ConstantWithinValidity);
    let state = |chemistry: &str| MaterialStateId {
        chemistry: chemistry.into(), phase:"solid".into(), process:"fixture".into(), revision:0,
    };
    let bytes = NormalizedInterfacePack::new(
        SurfaceSpec { material:state("cold-body"),texture_frame:"normal-plus-x".into() },
        SurfaceSpec { material:state("hot-body"),texture_frame:"normal-minus-x".into() },
        SystemContext {medium:"dry-contact".into(),third_body:None,
            environment:"fixture-air".into(),history:"unaged".into()}, claims).unwrap().to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack {kind:CardPackKind::Interface,
        source:"synthetic-contact-control".into(),bytes:bytes.clone(),expect:None}]).unwrap();
    (bytes,admitted.interfaces()[0].card().to_hex())
}
fn material_pack(nonlinear: bool) -> (Vec<u8>, String, String) {
    let dims = fs_conduction::CONDUCTIVITY_DIMS;
    let name = if nonlinear {"nonlinear-contact-solid"} else {"constant-contact-solid"};
    let value = if nonlinear { PropertyValue::Curve {abscissa:"T".into(),
        abscissa_dims:fs_project::spec::dims::TEMPERATURE,knots:vec![(200.0,2.0),(450.0,17.0)],dims} }
        else {PropertyValue::Scalar {value:10.0,dims}};
    let claims = pack(name,"thermal-conductivity",dims,value,if nonlinear {InterpolationPolicy::LinearInside}
        else {InterpolationPolicy::ConstantWithinValidity});
    let bytes = NormalizedMaterialCardPack::new(MaterialStateId {chemistry:name.into(),
        phase:"solid".into(),process:"synthetic".into(),revision:0},claims).unwrap().to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack {kind:CardPackKind::Material,
        source:"synthetic-contact-solid".into(),bytes:bytes.clone(),expect:None}]).unwrap();
    (bytes,admitted.materials()[0].card().to_hex(),admitted.materials()[0].identity().to_string())
}

struct Fixture { dir:PathBuf, project:fs_project::ProjectSpec }
impl Fixture {
    fn new(nonlinear: bool, coupled: bool) -> Self {
        let dir = scratch();
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let mut project = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
            root.join("examples/contact-pair/contact-pair.fsim")).unwrap()).unwrap().decoded.spec;
        let (bytes,card,state) = material_pack(nonlinear);
        std::fs::write(dir.join("solid.fsmcdpk"),bytes).unwrap();
        for binding in project.materials.as_mut().unwrap() {
            binding.card = card.clone(); binding.state = state.clone(); binding.claim = None;
            binding.conductivity_tolerance = None;
        }
        project.power.as_mut().unwrap()[0].watts.value = 12.0;
        project.power.as_mut().unwrap()[0].duty = 0.37;
        project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
        if coupled {
            use fs_project::{ThermalBoundaryCondition as B,spec::dims};
            use fs_qty::QtyAny;
            let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            setup.boundaries[0].condition = B::AirflowConvection {
                branch:"cold".into(),order:0,inlet_temperature:QtyAny::new(297.0,dims::TEMPERATURE),
                hydraulic_diameter:QtyAny::new(0.02,dims::LENGTH),flow_area:QtyAny::new(0.004,dims::AREA),
                channel_length:QtyAny::new(0.3,dims::LENGTH),correlation:"convection.gnielinski".into(),
            };
            let reference = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
                root.join("data/reference-project/cooling-radiation.fsim")).unwrap()).unwrap().decoded.spec;
            let mut radiation = reference.cooling.unwrap().conduction.unwrap().radiation.unwrap();
            radiation.surfaces[0].name = "cold-radiator".into();
            radiation.surfaces[0].target = "cold".into();
            radiation.surfaces[0].reservoir_temperature.value = 300.0;
            radiation.temperature_tolerance.value = 1e-11;
            radiation.heat_tolerance.value = 1e-10;
            radiation.max_iterations = 256;
            setup.radiation = Some(radiation);
        }
        Self {dir,project}
    }
    fn solve(&self, resistance:f64, output:Option<&str>, ordinal:usize) -> (J,J,String) {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let mut project = self.project.clone();
        let (bytes,card) = contact_pack(resistance);
        let binding = &mut project.interface_cards.as_mut().unwrap()[0];
        binding.card = card; binding.claim = None;
        if let Some(name) = output {
            project.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
                name:name.into(),kind:"report".into(),region:None,
            });
        }
        let source = self.dir.join(format!("project-{ordinal}.fsim"));
        let contact = self.dir.join(format!("contact-{ordinal}.fsintpk"));
        std::fs::write(&source,fs_project::print_sexpr(&project).unwrap()).unwrap();
        std::fs::write(&contact,bytes).unwrap();
        let ledger = self.dir.join("results.db");
        command(&["--json","import",source.to_str().unwrap(),
            root.join("examples/contact-pair/cold-body.stl").to_str().unwrap(),
            root.join("examples/contact-pair/hot-body.stl").to_str().unwrap(),
            ledger.to_str().unwrap(),"--unit","m","--max-hole-edges","0"]);
        let result = command(&["--json","solve",source.to_str().unwrap(),ledger.to_str().unwrap(),
            "--materials",self.dir.join("solid.fsmcdpk").to_str().unwrap(),
            "--materials",root.join("data/reference-project/gray-surface.fsmcdpk").to_str().unwrap(),
            "--interfaces",contact.to_str().unwrap()]);
        let run = result.str_field("run").unwrap().to_string();
        let receipt = artifact(&ledger,result.str_field("run_receipt").unwrap());
        let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
            .find(|s| s.str_field("stage") == Some("conduction")).unwrap();
        let receipt = artifact(&ledger,stage.str_field("receipt").unwrap());
        let field = artifact(&ledger,receipt.str_field("solution_artifact").unwrap());
        (receipt,field,run)
    }
}

#[test]
fn native_contact_controls_preserve_the_accepted_field_and_match_card_resolves() {
    for nonlinear in [false,true] { for coupled in [false,true] {
        let f = Fixture::new(nonlinear,coupled);
        let (plain,plain_field,plain_run) = f.solve(0.13,None,0);
        let (nominal,field,run) = f.solve(0.13,Some(OUTPUT),1);
        assert_ne!(plain_run,run,"explicit derivative request has a different project/run identity");
        assert_eq!(plain_field.get("temperature"),field.get("temperature"));
        for key in ["energy","conjugate","radiation","interfaces"] {assert_eq!(plain.get(key),nominal.get(key),"{key}");}
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("output"),Some(OUTPUT));
        assert_eq!(report.str_field("authority"),Some("Estimated"));
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        let parameters = report.get("parameters").unwrap().as_array().unwrap();
        let contacts:Vec<_> = parameters.iter().filter(|r|r.str_field("target")==Some(TARGET)).collect();
        assert_eq!(contacts.len(),1);
        let contact = contacts[0];
        assert_eq!(contact.str_field("entity"),Some("cold-hot-joint"));
        assert_eq!(contact.str_field("parameter_unit"),Some("1"));
        assert_eq!(contact.f64_field("reference_value"),Some(1.0));
        assert_eq!(contact.str_field("interface_card"),Some(contact_pack(0.13).1.as_str()));
        let actual = contact.f64_field("derivative").unwrap();
        let delta = 2e-4_f64;
        let temperature = |r:f64,n| f.solve(r,None,n).1.get("temperature").unwrap()
            .as_array().unwrap()[selected].as_f64().unwrap();
        let expected = (temperature(0.13*delta.exp(),2)-temperature(0.13*(-delta).exp(),3))/(2.0*delta);
        assert!(actual > 0.0 && expected > 0.0,"larger resistance impedes this hot body's only heat exit");
        assert!((actual-expected).abs()<5e-4*expected.abs().max(0.01),
            "nonlinear={nonlinear} coupled={coupled}: {actual:e} != {expected:e}");
        let (legacy,legacy_field,legacy_run) = f.solve(0.13,Some("temperature-max-adjoint"),4);
        assert_ne!(legacy_run,run);
        assert_eq!(legacy_field.get("temperature"),field.get("temperature"));
        let legacy = legacy.get("nominal_adjoint").unwrap();
        assert_eq!(legacy.str_field("output"),Some("temperature-max-adjoint"));
        let retained:Vec<_> = parameters.iter().filter(|r|r.str_field("target")!=Some(TARGET)).collect();
        assert_eq!(retained,legacy.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        for key in ["selected_vertex","value_k","true_relative_residual","dual_iterations","mode","unsupported"] {
            assert_eq!(report.get(key),legacy.get(key),"{key}");
        }
    } }
}
