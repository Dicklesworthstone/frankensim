//! Independent fan controls through real import, flow, thermal and study paths.
use super::*;
use fs_airflow::FanArrangement;
use fs_project::fansystem::FanSystemTopology;

fn multibank(parallel: bool) -> Fixture {
    let mut f=fixture();
    let system=f.project.cooling.as_mut().unwrap().fan_system.as_mut().unwrap();
    system.banks[0].count=2;
    system.banks[0].arrangement=FanArrangement::Series;
    let mut second=system.banks[0].clone();
    second.bank_id="second-bank".into();
    second.count=2; second.arrangement=FanArrangement::Parallel;
    second.speed_ratio=1.05; second.rated_point=None;
    for point in &mut second.curve.points { point.static_pressure.value*=1.15; }
    second.curve.source="synthetic heterogeneous fan regression; not manufacturer data".into();
    second.curve.source_id="native-multibank-regression-v1".into();
    system.banks.push(second);
    let names=vec!["fixture-bank".into(),"second-bank".into()];
    system.topology=if parallel { FanSystemTopology::Parallel(names) } else { FanSystemTopology::Series(names) };
    std::fs::write(f.sources.join("cooling-reference.fsim"),fs_project::print_sexpr(&f.project).unwrap()).unwrap();
    f
}

#[test]
fn native_independent_bank_adjoints_match_selected_vertex_physical_resolves() {
    for parallel in [false,true] {
        let mut f=multibank(parallel);
        f.project.power.as_mut().unwrap()[0].watts.value=300.0;
        let baseline=solve(&f,&f.project,0);
        let mut project=f.project.clone(); request(&mut project);
        let receipt=solve(&f,&project,1);
        let field=|receipt:&JsonValue| json_artifact(&f.dir.join("adjoint.db"),
            receipt.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone();
        assert_eq!(field(&baseline),field(&receipt));
        for key in ["energy","conjugate"] { assert_eq!(baseline.get(key),receipt.get(key)); }
        let goal=receipt.get("nominal_adjoint").unwrap();
        let vertex=goal.f64_field("selected_vertex").unwrap() as usize;
        let rows:Vec<_>=goal.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target")==Some("fan-speed-ratio")).collect();
        assert_eq!(rows.len(),2);
        assert!(!goal.get("unsupported").unwrap().as_array().unwrap().iter()
            .any(|row| row.str_field("target")==Some("fan-speed-ratio")));
        for (i,entity) in ["fixture-bank","second-bank"].into_iter().enumerate() {
            let matches:Vec<_>=rows.iter().filter(|row| row.str_field("entity")==Some(entity)).collect();
            assert_eq!(matches.len(),1);
            let row=matches[0];
            assert_eq!(row.f64_field("ordinal"),Some(i as f64));
            assert_eq!(row.str_field("parameter_unit"),Some("1"));
            let actual=row.f64_field("derivative").unwrap();
            let step=0.0005;
            let mut values=Vec::new();
            for (side,sign) in [-1.0,1.0].into_iter().enumerate() {
                let mut shifted=f.project.clone();
                shifted.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[i].speed_ratio+=sign*step;
                let receipt=solve(&f,&shifted,2+2*i+side);
                values.push(field(&receipt).as_array().unwrap()[vertex].as_f64().unwrap());
            }
            let expected=(values[1]-values[0])/(2.0*step);
            assert!(actual<0.0,"an independently accelerated bank cools this fixture");
            assert!((actual-expected).abs()<2e-3*expected.abs().max(1e-4),
                "parallel={parallel} {entity}: adjoint {actual:e} vs physical difference {expected:e}");
        }
    }
}

fn two_bank_study(controlled: bool, copula: bool, qmc: bool) -> String {
    study_text(copula,qmc,controlled).replace(
        "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
        "(uniform :name \"first\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.8 :high 0.9)")
        .replace(":target convection-temperature :entity \"air\" :low 294K :high 300K",
            ":target fan-speed-ratio :entity \"second-bank\" :low 1.0 :high 1.1")
}

#[test]
fn two_independent_bank_variables_share_one_calibration_without_changing_samples() {
    for (copula,qmc) in [(false,false),(true,false),(false,true),(true,true)] {
        let raw=multibank(true);
        std::fs::write(raw.sources.join("study.fsim"),two_bank_study(false,copula,qmc)).unwrap();
        let result=raw.study(None,fs_cli::exit::SUCCESS);
        let (_,baseline)=raw.report(&result);
        let f=multibank(true);
        std::fs::write(f.sources.join("study.fsim"),two_bank_study(true,copula,qmc)).unwrap();
        let result=f.study(None,fs_cli::exit::SUCCESS);
        let (_,report)=f.report(&result);
        for key in ["observations","statistics","qmc"] { assert_eq!(report.get(key),baseline.get(key),"{key}"); }
        let control=report.get("mean_control").unwrap();
        assert_eq!(control.f64_field("probe_solves_completed"),Some(1.0));
        let values=control.get("gradient").unwrap().as_array().unwrap();
        assert_eq!(values.len(),2);
        for value in values { assert!(value.as_f64().unwrap()<0.0); }
        assert_ne!(values[0],values[1],"heterogeneous banks retain independent controls");
    }
}
