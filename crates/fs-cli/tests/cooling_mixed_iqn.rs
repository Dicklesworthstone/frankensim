//! Actual-binary regressions for the production mixed-network IQN adapter.
//! The stiff steady case has only twelve primal and adjoint sweeps. Reference
//! temperatures and derivatives come from the independent slab/NTU equations.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::fs;
use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::sync::atomic::{AtomicUsize,Ordering};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static STIFF: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/stiff-mixed-contact.json"))));
// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static PULSE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/nonlinear-contact-pulse.json"))));
static NEXT:AtomicUsize=AtomicUsize::new(0);
fn scratch()->PathBuf {
    let nanos=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("frankensim-mixed-iqn-{}-{nanos}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    fs::create_dir(&path).unwrap(); path
}
fn run(dir:&Path,name:&str,text:&str)->Output {
    let path=dir.join(name); fs::write(&path,text).unwrap();
    Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json","cooling-network"])
        .arg(path).output().unwrap()
}
fn document(output:&Output)->J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn value(root:&J,path:&[&str])->f64 {root.path(path).and_then(J::as_f64).unwrap()}
fn close(a:f64,b:f64,tolerance:f64) {assert!((a-b).abs()<tolerance,"{a:.15e} != {b:.15e}");}
fn exact(h:f64,r:f64)->f64 {
    let c0=1.2*0.003*1007.0_f64; let ct=1.2*0.004*1007.0_f64;
    let f0=c0 * -(-50.0/c0).exp_m1(); let f1=ct * -(-0.01*h/ct).exp_m1();
    // Two 25-mm solids contribute 0.5 K/W, contact r/A, plus films.
    let total=0.5+r/0.01+1.0/f0+1.0/f1-1.0/ct;
    330.0-10.0/(total*f0)
}

#[test]
fn real_command_closes_stiff_mixing_and_keeps_total_contact_and_htc_gradients() {
    let dir=scratch(); let output=run(&dir,"stiff.json",&STIFF); let result=document(&output);
    assert_eq!(result.path(&["coupling_solver","method"]).and_then(J::as_str),Some("iqn-ils"));
    assert_eq!(result.path(&["coupling_solver","adjoint_method"]).and_then(J::as_str),Some("iqn-ils"));
    assert!(value(&result,&["coupling_iterations"])<=12.0);
    close(value(&result,&["objective","value_k"]),exact(800.0,0.01),3e-6);
    close(value(&result,&["robin_out_w"]),0.0,1e-7);
    let contacts=result.path(&["contact_sensitivities","rows"]).unwrap().as_array().unwrap();
    let delta=1e-4_f64;
    let expected_contact=(exact(800.0,0.01*delta.exp())-exact(800.0,0.01*(-delta).exp()))/(2.0*delta);
    close(contacts[0].f64_field("dobjective_dlog_resistance_k").unwrap(),expected_contact,3e-6);
    let walls=result.get("walls").unwrap().as_array().unwrap();
    let last=walls.iter().find(|row| row.str_field("region")==Some("last-face")).unwrap();
    let expected_h=(exact(800.0*delta.exp(),0.01)-exact(800.0*(-delta).exp(),0.01))/(2.0*delta);
    close(last.f64_field("dobjective_dlog_htc").unwrap(),expected_h,3e-6);
    // Same executable/request/runtime: fresh deterministic secants reproduce output.
    let replay=run(&dir,"replay.json",&STIFF); document(&replay); assert_eq!(output.stdout,replay.stdout);
}

#[test]
fn nonlinear_contact_pulse_retains_its_physical_peak_and_energy_under_acceleration() {
    let dir=scratch(); let result=document(&run(&dir,"pulse.json",&PULSE));
    assert_eq!(result.path(&["coupling_solver","method"]).and_then(J::as_str),Some("iqn-ils"));
    assert_eq!(result.path(&["coupling_solver","adjoint_method"]),Some(&J::Null));
    assert_eq!(result.get("contact_sensitivities"),Some(&J::Null));
    // Independent dense P1/contact/air calculation already used by the
    // nonlinear-transient regression; no change of timestep or physics here.
    close(value(&result,&["transient","sampled_peak_objective_k"]),306.1651033635,5e-5);
    close(value(&result,&["objective","value_k"]),301.41523251,5e-5);
    assert!(value(&result,&["transient","nonlinear","solid_solves"])>0.0);
}

#[test]
fn neither_primal_nor_adjoint_budget_exhaustion_can_publish_a_cooling_result() {
    let dir=scratch();
    for (index,(old,new)) in [
        ("\"coupling_iterations\":12","\"coupling_iterations\":1"),
        ("\"derivative_iterations\":12","\"derivative_iterations\":1"),
    ].into_iter().enumerate() {
        assert!(STIFF.contains(old));
        let output=run(&dir,&format!("short-{index}.json"),&STIFF.replace(old,new));
        assert!(!output.status.success()); assert!(output.stdout.is_empty());
        assert!(!output.stderr.is_empty());
    }
}
