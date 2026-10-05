//! Real native fan/network/thermal re-solves and pre-sampling calibration.
use super::*;

fn speed(project: &mut fs_project::ProjectSpec, value: f64) {
    project.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio=value;
}

fn fixture() -> Fixture {
    let mut f=Fixture::new();
    use fs_project::{ThermalBoundaryCondition, spec::dims};
    use fs_qty::QtyAny;
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::AirflowConvection {
            branch:"air".into(), order:0, inlet_temperature:QtyAny::new(297.0,dims::TEMPERATURE),
            hydraulic_diameter:QtyAny::new(0.02,dims::LENGTH), flow_area:QtyAny::new(0.004,dims::AREA),
            channel_length:QtyAny::new(0.3,dims::LENGTH), correlation:"convection.gnielinski".into(),
        };
    f.project.solver.as_mut().unwrap().tolerance_rel=1e-10;
    f.project.power.as_mut().unwrap()[0].duty=0.37;
    speed(&mut f.project,0.85);
    std::fs::write(f.sources.join("cooling-reference.fsim"),fs_project::print_sexpr(&f.project).unwrap()).unwrap();
    f
}

#[test]
fn native_fan_adjoint_is_an_absolute_ratio_derivative_at_nonunit_speed() {
    let mut f=fixture();
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
    let parameters:Vec<_>=goal.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target")==Some("fan-speed-ratio")).collect();
    assert_eq!(parameters.len(),1);
    assert_eq!(parameters[0].str_field("entity"),Some("fixture-bank"));
    assert_eq!(parameters[0].str_field("parameter_unit"),Some("1"));
    let actual=parameters[0].f64_field("derivative").unwrap();
    let step=0.001;
    let mut values=Vec::new();
    for (side,sign) in [-1.0,1.0].into_iter().enumerate() {
        let mut shifted=f.project.clone(); speed(&mut shifted,0.85+sign*step);
        let receipt=solve(&f,&shifted,2+side);
        values.push(field(&receipt).as_array().unwrap()[vertex].as_f64().unwrap());
    }
    let expected=(values[1]-values[0])/(2.0*step);
    assert!(actual<0.0,"higher fan speed cools this fixed-power fixture");
    assert!((actual-expected).abs()<2e-3*expected.abs().max(1e-4),"{actual:e} vs {expected:e}");
    assert!((0.85*actual-expected).abs()>0.1*expected.abs(),"must not report dJ/dln(speed) as dJ/dspeed");
}

fn fan_study(controlled: bool, copula: bool, qmc: bool) -> String {
    study_text(copula,qmc,controlled).replace(
        "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
        "(uniform :name \"speed\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.8 :high 1.2)")
        .replace("convection-temperature","air-inlet-temperature")
}

#[test]
fn native_fan_mean_control_uses_one_calibration_without_changing_probability_children() {
    for (copula,qmc) in [(false,false),(true,false),(false,true),(true,true)] {
        let raw=fixture();
        std::fs::write(raw.sources.join("study.fsim"),fan_study(false,copula,qmc)).unwrap();
        let result=raw.study(None,fs_cli::exit::SUCCESS);
        let (_,baseline)=raw.report(&result);
        let f=fixture();
        std::fs::write(f.sources.join("study.fsim"),fan_study(true,copula,qmc)).unwrap();
        let result=f.study(None,fs_cli::exit::SUCCESS);
        let (_,report)=f.report(&result);
        for key in ["observations","statistics","qmc"] {
            assert_eq!(report.get(key),baseline.get(key),"{key}");
        }
        let control=report.get("mean_control").unwrap();
        assert_eq!(control.str_field("method"),Some("nominal-adjoint"));
        assert_eq!(control.f64_field("probe_solves_completed"),Some(1.0));
        let gradient=control.get("gradient").unwrap().as_array().unwrap();
        assert_eq!(gradient.len(),2);
        assert!(gradient[0].as_f64().unwrap()<0.0,"fan coefficient is not replaced with zero");
        assert!((gradient[1].as_f64().unwrap()-1.0).abs()<1e-6);
    }
}
