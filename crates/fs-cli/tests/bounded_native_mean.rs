//! Real import/solve/report/package path for finite-validity material cards.
#![cfg(feature = "thermal-verification")]
#[path = "../src/json_read.rs"] mod json_read;
use json_read::JsonValue as J;
use fs_cli::cards::{CardPackKind, CardPackSet, RawCardPack};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use std::path::{Path, PathBuf};

fn command(args: &[&str]) -> J {
    let result = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS,
        "{args:?}\n{}\n{}", result.stdout, result.stderr);
    J::parse(&result.stdout).unwrap()
}
fn artifact(ledger: &Path, hash: &str) -> J {
    let bytes = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap()
        .get_artifact(&fs_blake3::ContentHash::from_hex(hash).unwrap()).unwrap().unwrap();
    J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}
fn pack(finite: bool) -> (Vec<u8>, String, String) {
    let source = fs_blake3::hash_bytes(b"bounded-native-mean-constant-fixture");
    let provenance = || Provenance { source: "synthetic constant PDE, not measured material".into(),
        license: "CC-BY-4.0; redistribution permitted".into(), artifact: Some(source) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "numerical-fixture".into(), method: "declared constant PDE".into(),
        artifact: source, caveats: "not experimental evidence".into(), provenance: provenance(),
    }).unwrap();
    let domain = fs_evidence::ValidityDomain::unconstrained();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("thermal-conductivity", fs_conduction::CONDUCTIVITY_DIMS),
        value: PropertyValue::Scalar { value:10.0, dims:fs_conduction::CONDUCTIVITY_DIMS },
        validity: if finite { domain.with("T",200.0,450.0) } else { domain },
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation], provenance: provenance(),
    }).unwrap();
    let pack = NormalizedPack::new("finite-constant-fixture", "declared-v1", source,
        "CC-BY-4.0; redistribution permitted", claims, Vec::new(), Vec::new()).unwrap();
    let bytes = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry:"constant-solid".into(), phase:"solid".into(),
        process:"numerical".into(), revision:0,
    }, pack).unwrap().to_bytes();
    let set = CardPackSet::admit(vec![RawCardPack { kind:CardPackKind::Material,
        source:"constant-fixture".into(), bytes:bytes.clone(), expect:None }]).unwrap();
    (bytes, set.materials()[0].card().to_hex(), set.materials()[0].identity().to_string())
}
struct Fixture { dir: PathBuf, project: fs_project::ProjectSpec, mesh: PathBuf }
impl Fixture {
    fn new(finite: bool, cube: bool) -> Self {
        let mut index=0;
        let dir=loop {
            let p=std::env::temp_dir().join(format!("fs-native-finite-mean-{}-{index}",std::process::id()));
            match std::fs::create_dir(&p) {
                Ok(())=>break p,
                Err(e) if e.kind()==std::io::ErrorKind::AlreadyExists=>index+=1,
                Err(e)=>panic!("fixture directory: {e}"),
            }
        };
        let root=Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let mut project=fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
            root.join("data/reference-project/cooling-reference.fsim")).unwrap()).unwrap().decoded.spec;
        let (bytes,card,state)=pack(finite);
        std::fs::write(dir.join("constant.fsmcdpk"),bytes).unwrap();
        for b in project.materials.as_mut().unwrap() {
            b.card=card.clone(); b.state=state.clone(); b.claim=None; b.conductivity_tolerance=None;
        }
        project.solver.as_mut().unwrap().tolerance_rel=1e-11;
        project.power.as_mut().unwrap()[0].watts.value=20.0;
        project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition=
            fs_project::ThermalBoundaryCondition::Convection {
                coefficient:fs_qty::QtyAny::new(20.0,fs_project::spec::dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature:fs_qty::QtyAny::new(300.0,fs_project::spec::dims::TEMPERATURE),
            };
        let mesh=root.join(if cube { "examples/contact-pair/cold-body.stl" }
            else { "data/reference-project/plate.stl" });
        Self {dir,project,mesh}
    }
    fn run(&self, p:&fs_project::ProjectSpec, ordinal:usize) -> fs_cli::CommandOutput {
        let source=self.dir.join(format!("case-{ordinal}.fsim")); let ledger=self.dir.join("runs.db");
        let material=self.dir.join("constant.fsmcdpk");
        std::fs::write(&source,fs_project::print_sexpr(p).unwrap()).unwrap();
        command(&["--json","import",source.to_str().unwrap(),self.mesh.to_str().unwrap(),
            ledger.to_str().unwrap(),"--unit","m","--max-hole-edges","0"]);
        fs_cli::run(["--json","solve",source.to_str().unwrap(),ledger.to_str().unwrap(),
            "--materials",material.to_str().unwrap()].iter().map(|s|s.to_string()).collect::<Vec<_>>())
    }
    fn solve(&self,p:&fs_project::ProjectSpec,ordinal:usize)->(J,J,J) {
        let output=self.run(p,ordinal);
        assert_eq!(output.exit_code,fs_cli::exit::SUCCESS,"{}{}",output.stdout,output.stderr);
        let result=J::parse(&output.stdout).unwrap(); let ledger=self.dir.join("runs.db");
        let run=artifact(&ledger,result.str_field("run_receipt").unwrap());
        let stage=run.get("stages").unwrap().as_array().unwrap().iter()
            .find(|s|s.str_field("stage")==Some("conduction")).unwrap();
        let receipt=artifact(&ledger,stage.str_field("receipt").unwrap());
        let field=artifact(&ledger,receipt.str_field("solution_artifact").unwrap());
        (result,receipt,field)
    }
}
fn request(p:&mut fs_project::ProjectSpec) {
    p.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
        name:"temperature-volume-mean-bound".into(),kind:"report".into(),region:None,
    });
}

#[test]
fn finite_card_bound_preserves_the_physical_field_and_reaches_retained_reports() {
    let fixture=Fixture::new(true,true); let mut p=fixture.project.clone();
    for binding in p.materials.as_mut().unwrap() {
        binding.temp_lo.value=290.0; binding.temp_hi.value=301.0;
    }
    // Unit cube: the old one-ended bound is 302 K; using both cooled ends
    // gives 300.75 K. Only the latter fits the retained 301 K domain.
    let (_,baseline,field)=fixture.solve(&p,0); request(&mut p);
    let (run,receipt,bounded_field)=fixture.solve(&p,1);
    assert_eq!(field.get("temperature"),bounded_field.get("temperature"));
    assert_eq!(baseline.get("energy"),receipt.get("energy"));
    let bound=receipt.get("volume_mean_bound").unwrap();
    let validity=bound.get("material_validity").unwrap();
    assert_eq!(validity.str_field("method"),Some("continuum-slab-comparison"));
    let range=validity.get("continuum_temperature_range_k").unwrap();
    let span=validity.get("retained_material_domain_k").unwrap();
    assert!(range.f64_field("lower").unwrap()>=span.f64_field("lower").unwrap());
    assert!(range.f64_field("upper").unwrap()<=span.f64_field("upper").unwrap());
    assert!(range.f64_field("upper").unwrap()>300.0);
    assert!(range.f64_field("upper").unwrap()<300.7500001);
    assert_eq!(validity.get("two_sided"),Some(&J::Bool(true)));
    assert_eq!(bound.str_field("authority"),Some("Verified"));
    let material=&bound.get("materials").unwrap().as_array().unwrap()[0];
    assert_eq!(material.str_field("card"),Some(p.materials.as_ref().unwrap()[0].card.as_str()));
    assert!(material.str_field("receipt_bytes_hex").unwrap().len()>100);
    assert_ne!(material.get("temperature_domain_k"),Some(&J::Null));
    let (_,replay,_)=fixture.solve(&p,2);
    assert_eq!(bound,replay.get("volume_mean_bound").unwrap());
    let ledger=fixture.dir.join("runs.db"); let id=run.str_field("run").unwrap();
    let exported=command(&["--json","report",id,ledger.to_str().unwrap()]);
    let report=J::parse(&std::fs::read_to_string(exported.str_field("report_json").unwrap()).unwrap()).unwrap();
    let mean=report.get("qois").unwrap().as_array().unwrap().iter()
        .find(|q|q.str_field("name")==Some("temperature-volume-mean-bound")).unwrap();
    assert_eq!(mean.str_field("color"),Some("Verified"));
    assert_eq!(mean.get("total_budget"),Some(&J::Null),"unknown physical uncertainty is not certified");
    let exported=command(&["--json","package",id,ledger.to_str().unwrap()]);
    let package=fs_package::EvidencePackage::from_json(&std::fs::read_to_string(
        exported.str_field("package").unwrap()).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
}

#[test]
fn finite_and_global_claims_use_the_same_nominal_equations() {
    let finite=Fixture::new(true,true); let global=Fixture::new(false,true);
    let mut a=finite.project.clone(); let mut b=global.project.clone(); request(&mut a); request(&mut b);
    let (_,ar,af)=finite.solve(&a,0); let (_,br,bf)=global.solve(&b,0);
    assert_eq!(af.get("temperature"),bf.get("temperature"));
    let a=ar.get("volume_mean_bound").unwrap(); let b=br.get("volume_mean_bound").unwrap();
    assert_eq!(a.get("enclosure_k"),b.get("enclosure_k"));
    assert_eq!(b.get("material_validity").unwrap().str_field("method"),Some("global-constant-source"));
}

#[test]
fn unsupported_geometry_and_insufficient_barriers_remain_refusals() {
    for cube in [false,true] {
        let f=Fixture::new(true,cube); let mut p=f.project.clone();
        if cube { for b in p.materials.as_mut().unwrap() { b.temp_lo.value=290.0; b.temp_hi.value=300.7; } }
        f.solve(&p,0); // the actual nominal field is admissible before the bound request
        request(&mut p);
        let result=f.run(&p,1);
        assert_ne!(result.exit_code,fs_cli::exit::SUCCESS);
        assert!(format!("{}{}",result.stdout,result.stderr).contains("finite spans"));
    }
}
