#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/radiative-contact-pulse.json"));
const SIGMA: f64 = 5.670_374_419e-8;
fn num(v: f64) -> J { J::Number { value: v, raw: v.to_string() } }
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object") };
    if let Some((_, slot)) = fields.iter_mut().find(|(name, _)| name == key) { *slot = value; }
    else { fields.push((key.into(), value)); }
}
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object") };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object") }; fields.retain(|(name, _)| name != key);
}
fn text(root: &J) -> String {
    match root {
        J::Null => "null".into(), J::Bool(b) => b.to_string(), J::Number { raw, .. } => raw.clone(),
        J::Str(s) => format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")),
        J::Array(a) => format!("[{}]", a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(o) => format!("{{{}}}", o.iter().map(|(k,v)| format!("{}:{}", text(&J::Str(k.clone())), text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(root: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(root).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn success(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run(root: &J) -> J { success(&output(root)) }
fn n(root: &J, key: &str) -> f64 { root.f64_field(key).unwrap() }
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} versus {b}"); }
fn short() -> J {
    let mut root = J::parse(BASE).unwrap();
    let schedule = member(&mut root, "transient");
    put(schedule, "intervals", J::parse(r#"[{"duration_s":4,"power_scale":1,"fan_speed_ratio":1},{"duration_s":6,"power_scale":0,"fan_speed_ratio":1.5}]"#).unwrap());
    root
}
fn energy(result: &J, key: &str) {
    let r = result.get(key).unwrap();
    let residual = n(r,"stored_energy_change_j")-n(r,"input_energy_j")
        +n(r,"air_energy_gain_j")+n(r,"radiative_energy_loss_j");
    near(residual,n(r,"energy_residual_j"),1e-9);
    assert!(residual.abs() < 1e-4);
}

/// One tetrahedron, all faces on the same air/radiation patch. Choose the
/// nodal P1 source so T_new=310 K is EXACT in the discrete weak form. The
/// consistent mass inverse is (20/V)*(I-11^T/5); no production matrix is used.
fn manufactured(ambient: f64) -> (J, f64) {
    let opposite = [3.0_f64.sqrt()/2.0, 0.5, 0.5, 0.5];
    let area: f64 = opposite.iter().sum();
    let ntu = 2.0*area;
    let q_air_density = 10.0*(1.0-(-ntu).exp())/area;
    let q_rad_density = 0.8*SIGMA*(310.0_f64.powi(4)-ambient.powi(4));
    let surface_load: Vec<f64> = opposite.iter().map(|a| (area-a)/3.0*(q_air_density+q_rad_density)).collect();
    let sum: f64 = surface_load.iter().sum();
    let source = surface_load.iter().map(|v|num(1000.0+120.0*(v-sum/5.0))).collect();
    let mut root = J::parse(r#"{
      "schema":"frankensim.cooling-network.v1","units":"SI","seed":"17",
      "budgets":{"graph_sweeps":100,"coupling_iterations":100,"linear_iterations":10000,"derivative_iterations":100,"wall_seconds":120},
      "tolerances":{"flow_m3_s":1e-12,"heat_w":1e-7,"temperature_k":1e-9,"linear_relative":1e-12,"relaxation":1},
      "air":{"density_kg_m3":1,"specific_heat_j_kg_k":1},
      "hydraulics":{"node_count":2,"boundaries":[{"node":0,"pressure_pa":1,"temperature_k":300},{"node":1,"pressure_pa":0}],
        "branches":[{"name":"air","from":0,"to":1,"resistance_pa_s2_m6":1,"source":"analytic","regions":["wall"]}]},
      "solid":{"vertices_m":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],"tetrahedra":[[0,1,2,3]],"conductivity_w_m_k":10,"nodal_source_w_m3":[],
        "adiabatic_remainder":false,"surfaces":[{"name":"wall","faces":[[1,2,3],[0,2,3],[0,1,3],[0,1,2]],"htc_w_m2_k":2}]},
      "objective":{"max_solid_temperature":true,"gradient":false},
      "transient":{"initial_temperature_k":300,"volumetric_heat_capacity_j_m3_k":1000,"max_step_s":10,"max_steps":8,
        "intervals":[{"duration_s":10,"power_scale":1}]},
      "radiation":{"max_iterations":128,"temperature_tolerance_k":1e-10,"relaxation":0.5,
        "surfaces":[{"surface":"wall","emissivity":0.8,"ambient_temperature_k":280,"source":"manufactured endpoint"}]}
    }"#).unwrap();
    put(member(&mut root,"solid"),"nodal_source_w_m3",J::Array(source));
    let J::Array(rows)=member(member(&mut root,"radiation"),"surfaces") else {panic!()};
    put(&mut rows[0],"ambient_temperature_k",num(ambient));
    (root,area*q_rad_density)
}

#[test]
fn manufactured_endpoint_uses_new_temperature_radiation_and_immutable_old_storage() {
    for ambient in [280.0,350.0] {
        let (input,expected_radiation)=manufactured(ambient);
        let result=run(&input);energy(&result,"transient");
        for value in result.get("solid_temperatures_k").unwrap().as_array().unwrap() {
            near(value.as_f64().unwrap(),310.0,2e-6);
        }
        let time=result.get("transient").unwrap();
        near(n(time,"stored_energy_change_j"),1000.0/6.0*10.0,1e-4);
        near(n(time,"radiative_energy_loss_j"),10.0*expected_radiation,1e-5);
        near(n(result.get("radiation").unwrap(),"radiative_out_w"),expected_radiation,1e-5);
        assert_eq!(n(time,"steps"),1.0);
        let mut frozen=input.clone();put(member(&mut frozen,"radiation"),"max_iterations",num(1.0));
        let refused=output(&frozen);assert_eq!(refused.status.code(),Some(6));assert!(refused.stdout.is_empty());
    }
}

#[test]
fn nonlinear_contact_trajectory_keeps_radiation_out_of_the_air_ledger() {
    let input=short();let result=run(&input);energy(&result,"transient");
    let t=result.get("transient").unwrap();
    let rows=t.get("history").unwrap().as_array().unwrap();
    let integrated: f64=rows.iter().skip(1).map(|r|n(r,"dt_s")*n(r,"radiative_heat_w")).sum();
    near(integrated,n(t,"radiative_energy_loss_j"),1e-10);
    assert!(integrated>1.0);
    let outer:f64=rows.iter().skip(1).map(|r|n(r,"coupling_iterations")).sum();
    assert!(n(t,"forward_solid_solves")>outer);
    assert_eq!(n(t,"forward_solid_solves"),n(t.get("nonlinear").unwrap(),"solid_solves"));
    let mut plain=input.clone();remove(&mut plain,"radiation");let control=run(&plain);
    assert!(n(t,"sampled_peak_objective_k")<n(control.get("transient").unwrap(),"sampled_peak_objective_k"));
    assert!(control.get("radiation").is_none());
    assert!(control.get("transient").unwrap().get("radiative_energy_loss_j").is_none());
}

#[test]
fn repeated_radiation_equals_one_explicit_unrolled_trajectory() {
    let mut repeated=short();
    put(member(&mut repeated,"transient"),"repeat",J::parse(r#"{"cycles":2,"max_total_steps":100}"#).unwrap());
    let mut unrolled=repeated.clone();remove(member(&mut unrolled,"transient"),"repeat");
    let J::Array(rows)=member(member(&mut unrolled,"transient"),"intervals") else {panic!()};
    rows.extend(rows.clone());
    let a=run(&repeated);let b=run(&unrolled);energy(&a,"repeated_cycles");energy(&b,"transient");
    assert_eq!(a.get("solid_temperatures_k"),b.get("solid_temperatures_k"));
    let cycles=a.get("repeated_cycles").unwrap();
    near(n(cycles,"sampled_peak_objective_k"),n(b.get("transient").unwrap(),"sampled_peak_objective_k"),1e-10);
    let sum:f64=cycles.get("cycles").unwrap().as_array().unwrap().iter().map(|r|n(r,"radiative_energy_loss_j")).sum();
    near(sum,n(cycles,"radiative_energy_loss_j"),1e-10);
    near(sum,n(b.get("transient").unwrap(),"radiative_energy_loss_j"),1e-8);
    assert!(sum>n(a.get("transient").unwrap(),"radiative_energy_loss_j"));
}

#[test]
fn adaptive_coarse_trials_do_not_enter_radiation_or_physical_history() {
    let (mut adaptive,_)=manufactured(280.0);
    put(member(&mut adaptive,"transient"),"adaptive",J::parse(r#"{"absolute_tolerance_k":1000,"relative_tolerance":0,"minimum_trial_step_s":0.01,"max_trials":100}"#).unwrap());
    let mut fine=adaptive.clone();remove(member(&mut fine,"transient"),"adaptive");
    put(member(&mut fine,"transient"),"max_step_s",num(5.0));
    let a=run(&adaptive);let b=run(&fine);energy(&a,"transient");
    assert_eq!(a.get("solid_temperatures_k"),b.get("solid_temperatures_k"));
    let ta=a.get("transient").unwrap();let tb=b.get("transient").unwrap();
    near(n(ta,"radiative_energy_loss_j"),n(tb,"radiative_energy_loss_j"),1e-8);
    assert_eq!(n(ta,"steps"),2.0);
    assert!(n(ta,"forward_solid_solves")>n(tb,"forward_solid_solves"));
    put(member(member(&mut adaptive,"transient"),"adaptive"),"absolute_tolerance_k",num(0.002));
    put(member(&mut adaptive,"transient"),"max_steps",num(1000.0));
    let rejected=run(&adaptive);energy(&rejected,"transient");
    assert!(n(rejected.path(&["transient","adaptive"]).unwrap(),"rejected_trials")>0.0);
}

#[test]
fn hotter_surroundings_can_add_energy_with_no_workload() {
    let mut input=short();
    let J::Array(patches)=member(member(&mut input,"radiation"),"surfaces") else {panic!()};
    for patch in patches {put(patch,"ambient_temperature_k",num(350.0));}
    let J::Array(rows)=member(member(&mut input,"transient"),"intervals") else {panic!()};
    for row in rows {put(row,"power_scale",num(0.0));}
    let result=run(&input);energy(&result,"transient");
    let t=result.get("transient").unwrap();
    assert!(n(t,"radiative_energy_loss_j")<0.0);assert!(n(t,"stored_energy_change_j")>0.0);
    assert!(n(t,"sampled_peak_objective_k")>300.0);near(n(t,"input_energy_j"),0.0,1e-12);
}

#[test]
fn radiative_transient_replay_and_patch_order_are_identical() {
    let mut input=short();let full=output(&input);success(&full);
    let replay=output(&input);success(&replay);assert_eq!(full.stdout,replay.stdout);
    let J::Array(rows)=member(member(&mut input,"radiation"),"surfaces") else {panic!()};rows.reverse();
    let reordered=output(&input);success(&reordered);assert_eq!(full.stdout,reordered.stdout);
}

#[test]
fn unsupported_time_derivatives_missing_material_policy_and_budgets_do_not_publish() {
    let input=short();
    let mut derivative=input.clone();
    put(member(&mut derivative,"transient"),"adjoint",J::parse(r#"{"qoi":"sampled-peak","max_checkpoint_bytes":1048576}"#).unwrap());
    put(member(&mut derivative,"transient"),"adaptive",J::parse(r#"{"absolute_tolerance_k":0.01,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":100}"#).unwrap());
    let result=output(&derivative);assert!(!result.status.success());assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("fixed timesteps"));
    let mut missing=input.clone();remove(member(&mut missing,"transient"),"nonlinear");
    let result=output(&missing);assert!(!result.status.success());assert!(result.stdout.is_empty());
    let mut cancelled=input.clone();put(member(&mut cancelled,"budgets"),"wall_seconds",num(1e-12));
    let result=output(&cancelled);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
}

static NEXT:AtomicU64=AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path=std::env::temp_dir().join(format!("fs-transient-radiation-uq-{}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap();Self(path)
    }
    fn uq(&self,plan:&str,extras:&[&str])->Output {
        std::fs::write(self.0.join("uq.json"),plan).unwrap();
        Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json","cooling-network-uq","base.json","uq.json"]).args(extras).output().unwrap()
    }
}
impl Drop for Scratch {fn drop(&mut self){let _=std::fs::remove_dir_all(&self.0);}}
fn uq_plan(vary:bool)->String {
    let (lo,hi)=if vary {(0.5,0.95)}else{(0.85,0.85)};
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":600,"qoi":{{"kind":"transient-sampled-peak"}},"temperature_limit_k":310,"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"radiation-emissivity","surface":"first-face"}},"distribution":{{"kind":"uniform","lo":{lo},"hi":{hi}}}}]}}"#)
}

#[test]
fn repeated_peak_uq_uses_real_radiative_trajectories_and_replays_exactly() {
    let dir=Scratch::new();let mut input=short();
    put(member(&mut input,"transient"),"repeat",J::parse(r#"{"cycles":2,"max_total_steps":100}"#).unwrap());
    std::fs::write(dir.0.join("base.json"),text(&input)).unwrap();
    let direct=run(&input);let zero=success(&dir.uq(&uq_plan(false),&[]));
    near(n(&zero,"mean_k"),n(direct.get("repeated_cycles").unwrap(),"sampled_peak_objective_k"),1e-10);
    assert_eq!(n(&zero,"std_dev_k"),0.0);
    let full=dir.uq(&uq_plan(true),&["--checkpoint","full.uqcp"]);let parsed=success(&full);
    assert!(n(&parsed,"std_dev_k")>1e-7);
    let part=dir.uq(&uq_plan(true),&["--checkpoint","part.uqcp","--max-new-samples","2"]);
    assert_eq!(part.status.code(),Some(6));
    let saved=std::fs::read(dir.0.join("part.uqcp")).unwrap();
    let done=dir.uq(&uq_plan(true),&["--resume","part.uqcp","--checkpoint","done.uqcp"]);success(&done);
    assert_eq!(full.stdout,done.stdout);
    assert_eq!(std::fs::read(dir.0.join("full.uqcp")).unwrap(),std::fs::read(dir.0.join("done.uqcp")).unwrap());
    assert_eq!(saved,std::fs::read(dir.0.join("part.uqcp")).unwrap());
}

#[test]
fn interrupted_radiative_trajectory_does_not_become_a_partial_peak_observation() {
    let dir=Scratch::new();std::fs::write(dir.0.join("base.json"),text(&short())).unwrap();
    let plan=uq_plan(true).replace("\"wall_seconds\":600","\"wall_seconds\":0.000000000001");
    let interrupted=dir.uq(&plan,&["--checkpoint","empty.uqcp"]);
    assert_eq!(interrupted.status.code(),Some(6));
    let progress=J::parse(std::str::from_utf8(&interrupted.stdout).unwrap()).unwrap();
    assert_eq!(n(&progress,"samples_evaluated"),0.0);
    let full=dir.uq(&uq_plan(true),&[]);success(&full);
    let done=dir.uq(&uq_plan(true),&["--resume","empty.uqcp"]);success(&done);
    assert_eq!(full.stdout,done.stdout);
}

#[path = "cooling_transient_radiation/adjoint.rs"]
mod adjoint;
