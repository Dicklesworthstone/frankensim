//! Real native geometry/card/thermal solves; no mocked calibration producer.
use super::*;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_ledger::Ledger;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, SurfaceSpec, SystemContext, UncertaintyModel};

const SOURCE: &str = r#"(fsim-uncertainty-study :version 1 :project "contact.fsim"
 :samples 4 :seed 29 :wall-time 120s :method monte-carlo :correlation independent
 :qoi "temperature-max" :geometry (
  (mesh :role "cold-body" :path "cold-body.stl" :unit "m" :max-hole-edges 0)
  (mesh :role "hot-body" :path "hot-body.stl" :unit "m" :max-hole-edges 0))
 :materials ("solid.fsmcdpk") :interfaces ("pressure.fsintpk")
 :parameters (
  (uniform :name "clamp" :target contact-pressure :entity "cold-hot-joint" :low 700000Pa :high 1100000Pa)
  (uniform :name "power" :target power :entity "hot" :low 4W :high 6W)))"#;

fn pressure_pack() -> Vec<u8> {
    let source = fs_blake3::hash_bytes(b"synthetic manufactured-state calibration test");
    let provenance = || Provenance { source:"synthetic pressure curve; not experimental data".into(),
        license:"CC0-1.0".into(),artifact:Some(source) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen:"numerical fixture".into(),method:"piecewise affine verification law".into(),
        artifact:source,caveats:"not measured contact physics".into(),provenance:provenance(),
    }).unwrap();
    let dims=fs_project::CONTACT_RESISTANCE_DIMS;
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(fs_project::CONTACT_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Curve {abscissa:"normal_pressure".into(),abscissa_dims:fs_project::spec::dims::PRESSURE,
            knots:vec![(200000.0,0.20),(1200000.0,0.10),(2000000.0,0.06)],dims},
        validity:fs_evidence::ValidityDomain::unconstrained().with("T",200.0,450.0)
            .with("normal_pressure",200000.0,2000000.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![observation],provenance:provenance(),
    }).unwrap();
    let pack=NormalizedPack::new("pressure-calibration","synthetic-v1",source,"CC0-1.0",
        claims,Vec::new(),Vec::new()).unwrap();
    let surface=|name:&str|SurfaceSpec {material:MaterialStateId {chemistry:name.into(),phase:"solid".into(),
        process:"fixture".into(),revision:0},texture_frame:"fixture".into()};
    NormalizedInterfacePack::new(surface("cold-body"),surface("hot-body"),SystemContext {
        medium:"dry-contact".into(),third_body:None,environment:"air".into(),history:"unaged".into()},pack)
        .unwrap().to_bytes()
}
struct Fixture { root:PathBuf, inputs:PathBuf }
impl Fixture {
    fn new() -> Self {
        static NEXT:AtomicU64=AtomicU64::new(0);
        let root=loop {
            let p=std::env::temp_dir().join(format!("fs-joint-uq-{}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
            match std::fs::create_dir(&p) {Ok(())=>break p,
                Err(e) if e.kind()==std::io::ErrorKind::AlreadyExists=>continue,Err(e)=>panic!("{e}")}
        };
        let inputs=root.join("inputs");std::fs::create_dir(&inputs).unwrap();
        let repo=Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        for name in ["cold-body.stl","hot-body.stl"] {
            std::fs::copy(repo.join("examples/contact-pair").join(name),inputs.join(name)).unwrap();
        }
        std::fs::copy(repo.join("data/reference-project/aa6061.fsmcdpk"),inputs.join("solid.fsmcdpk")).unwrap();
        let bytes=pressure_pack();
        let cards=crate::CardPackSet::admit(vec![crate::RawCardPack {kind:crate::CardPackKind::Interface,
            source:"synthetic-pressure".into(),bytes:bytes.clone(),expect:None}]).unwrap();
        std::fs::write(inputs.join("pressure.fsintpk"),bytes).unwrap();
        let mut project=fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
            repo.join("examples/contact-pair/contact-pair.fsim")).unwrap()).unwrap().decoded.spec;
        let binding=&mut project.interface_cards.as_mut().unwrap()[0];
        binding.card=cards.interfaces()[0].card().to_hex();binding.claim=None;
        if let fs_project::InterfaceState::DryContact {pressure,..}=&mut binding.state {
            pressure.value=1500000.0; // Deliberately NOT the probability mean; different source slope.
        }
        std::fs::write(inputs.join("contact.fsim"),fs_project::print_sexpr(&project).unwrap()).unwrap();
        Self {root,inputs}
    }
    fn run(&self, source:&str, name:&str, expected:u8) -> (J,PathBuf) {
        let path=self.inputs.join(format!("{name}.fsim"));std::fs::write(&path,source).unwrap();
        let ledger=self.root.join(format!("{name}.db"));
        let output=crate::study::study_path(&path,&ledger,None,crate::OutputMode::Json);
        assert_eq!(output.exit_code,expected,"{}\n{}",output.stdout,output.stderr);
        if crate::SOLVE_DRIVER_VERSION < 44 {
            assert!(output.stderr.contains("cli-uncertainty-contact-driver"),"{}",output.stderr);
        }
        (J::parse(&output.stdout).unwrap(),ledger)
    }
    fn report(output:&J,path:&Path) -> J {
        let id=ContentHash::from_hex(output.get("receipt").unwrap().str_field("report_json").unwrap()).unwrap();
        let bytes=Ledger::open(path.to_str().unwrap()).unwrap().get_artifact(&id).unwrap().unwrap();
        J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
    }
}

#[test]
fn joint_studies_preserve_raw_observations_and_use_the_sampled_card_chain_rule() {
    let f=Fixture::new();
    // The dependency integration is atomic with driver 44. Before it lands,
    // exercising the public command MUST refuse without creating a ledger.
    if crate::SOLVE_DRIVER_VERSION < 44 {
        let (_,ledger)=f.run(SOURCE,"not-yet-state-aware",crate::exit::REFUSED);
        assert!(!ledger.exists());
        return;
    }
    for qmc in [false,true] {
        let raw=if qmc {SOURCE.replace("monte-carlo","quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
            .replace("independent","(gaussian-copula :latent-correlation ((1 0.5) (0.5 1)))")}
            else {SOURCE.to_string()};
        let controlled=raw.replace(":version 1",":version 3 :mean-control (nominal-adjoint :max-solves 1)");
        let tag=if qmc {"qmc"} else {"mc"};
        let (a,ap)=f.run(&raw,&format!("raw-{tag}"),crate::exit::SUCCESS);
        let (b,bp)=f.run(&controlled,&format!("control-{tag}"),crate::exit::SUCCESS);
        let a=Fixture::report(&a,&ap);let b=Fixture::report(&b,&bp);
        assert_eq!(a.get("observations"),b.get("observations"));
        let control=b.get("mean_control").unwrap();
        assert_eq!(control.f64_field("probe_solves_completed"),Some(1.0));
        assert_eq!(control.get("parameter_units").unwrap().as_array().unwrap()[0].as_str(),Some("Pa"));
        let gradient=control.get("gradient").unwrap().as_array().unwrap()[0].as_f64().unwrap();
        assert!(gradient<0.0,"increasing pressure improves this card's contact");
        let model=Model::load(&f.inputs.join(format!("control-{tag}.fsim"))).unwrap();
        let ledger=Ledger::open(bp.to_str().unwrap()).unwrap();let gate=CancelGate::new_clock_free();
        let lo=model.sample(&ledger,&gate,&[899990.0,5.0],120.0).unwrap().unwrap();
        let hi=model.sample(&ledger,&gate,&[900010.0,5.0],120.0).unwrap().unwrap();
        let fd=(hi.value_k-lo.value_k)/20.0;
        assert!((gradient-fd).abs()<5e-4*fd.abs().max(1e-9),"{gradient:e} != {fd:e}");
    }
}

#[test]
fn joint_calibration_resume_is_source_free_and_nonsmooth_sources_refuse_before_sampling() {
    let f=Fixture::new();
    if crate::SOLVE_DRIVER_VERSION < 44 { return; } // Covered by the explicit public gate test above.
    let source=SOURCE.replace(":version 1",":version 3 :mean-control (nominal-adjoint :max-solves 1)");
    let path=f.inputs.join("study.fsim");std::fs::write(&path,&source).unwrap();
    let ledger=f.root.join("resume.db");
    let prefix=crate::study::study_path(&path,&ledger,Some("1"),crate::OutputMode::Json);
    assert_eq!(prefix.exit_code,crate::exit::BUDGET,"{}",prefix.stderr);
    let prefix=J::parse(&prefix.stdout).unwrap();
    let report=Fixture::report(&prefix,&ledger);
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    // An unequal-slope source knot is a valid primal state but cannot supply
    // a nominal two-sided derivative or silently trigger secant fallback.
    let kink=source.replace("700000Pa","1100000Pa").replace("1100000Pa)","1300000Pa)");
    let (k,kp)=f.run(&kink,"kink",crate::exit::REFUSED);
    let failed=Fixture::report(&k,&kp);
    assert!(failed.get("observations").unwrap().as_array().unwrap().is_empty());
    assert!(failed.str_field("failure").unwrap().contains("unequal-slope"));
    std::fs::rename(&f.inputs,f.root.join("moved-inputs")).unwrap();
    let complete=crate::study::resume_path(prefix.str_field("run").unwrap(),&ledger,None,crate::OutputMode::Json);
    assert_eq!(complete.exit_code,crate::exit::SUCCESS,"{}",complete.stderr);
    let complete=J::parse(&complete.stdout).unwrap();
    let resumed=Fixture::report(&complete,&ledger);
    assert_eq!(resumed.get("mean_control").unwrap().get("gradient"),report.get("mean_control").unwrap().get("gradient"));
    assert_eq!(resumed.get("observations").unwrap().as_array().unwrap().len(),4);
}

#[path = "mixed.rs"]
mod mixed;

#[path = "bands.rs"]
mod bands;
