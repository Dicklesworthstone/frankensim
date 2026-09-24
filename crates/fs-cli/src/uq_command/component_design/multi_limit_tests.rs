use super::*;
use thermal_constraints::{Assessment, Constraint, Envelope};
use constraint_execution::Work;

fn complete(p: &Plan, v: &[f64], work: &mut Work) -> Result<Evaluation> {
    let primary = 300.0 + 0.5*v[0] + v[1];
    let memory = 290.0 + 2.0*v[0] + 0.25*v[1];
    let constraints = Envelope::new(vec![
        Assessment { constraint: Constraint { name:"primary".into(),limit:p.limit,objective:J::Null },
            peak:primary,time:20.0,slopes:vec![Some(0.5),Some(1.0)] },
        Assessment { constraint: Constraint { name:"memory".into(),limit:296.0,objective:J::Null },
            peak:memory,time:18.0,slopes:vec![Some(2.0),Some(0.25)] },
    ])?;
    for _ in 0..2 { work.begin()?;work.record(p.planned_steps,7)?; }
    Ok(Evaluation {document:J::Null,peak:primary,peak_time:20.0,
        steps:p.planned_steps*2,solves:14,slopes:constraints.slopes(),constraints:Some(constraints)})
}

#[test]
fn allocation_tracks_the_intersection_when_a_colder_constraint_takes_over() {
    let p=plan();
    let result=search::allocate_with_work(&p,deadline(),2,|v,w|complete(&p,v,w)).unwrap();
    assert!(result.reason.is_none());assert_eq!(result.completed,2);
    assert!(result.values[0]>2.9999 && result.values[0]<=3.0);
    assert!(result.values[1]<0.001);
    assert!(result.passing.peak<p.limit-2.0);
    assert!(result.passing.margin(&p).unwrap()<=0.0);
    assert_eq!(result.passing.active_constraint(),"memory");
    assert_eq!(result.work.completed,2*result.history.len());
    assert_eq!(result.steps,2*p.planned_steps*result.history.len());
    assert_eq!(result.history[0].str_field("active_constraint"),Some("primary"));
    assert!(result.history.iter().any(|row|row.str_field("active_constraint")==Some("memory")));
    assert!(result.newton_trials>0);
}

#[test]
fn whole_constraint_set_is_admitted_before_launch_and_partial_set_is_never_selected() {
    let mut p=plan();p.max_total_steps=p.planned_steps;
    let mut calls=0;
    assert!(search::allocate_with_work(&p,deadline(),2,|v,w|{calls+=1;complete(&p,v,w)}).is_err());
    assert_eq!(calls,0);
    p.max_total_steps=100000;
    let mut calls=0;
    let partial=search::allocate_with_work(&p,deadline(),2,|v,w| {
        calls+=1;
        if calls==2 {
            w.begin()?;w.record(p.planned_steps,7)?;
            w.begin()?;
            return Err(budget("second temperature objective interrupted"));
        }
        complete(&p,v,w)
    }).unwrap();
    assert!(partial.reason.is_some());assert_eq!(partial.values,vec![0.0,0.0]);
    assert_eq!(partial.attempted,2);assert_eq!(partial.history.len(),1);
    assert_eq!(partial.work.attempted,4);assert_eq!(partial.work.completed,3);
    assert_eq!(partial.steps,3*p.planned_steps);assert_eq!(partial.solves,21);
    let mut calls=0;
    assert!(search::allocate_with_work(&p,deadline(),2,|v,w| {
        calls+=1;if calls==2 { return Err(model_failure("secondary physics failed")); }
        complete(&p,v,w)
    }).is_err());
}

#[test]
fn incomplete_success_and_unusable_active_derivatives_do_not_bypass_any_limit() {
    let p=plan();
    assert!(search::allocate_with_work(&p,deadline(),2,|v,w| {
        w.begin()?;let e=fake(&p,v);w.record(e.steps,e.solves)?;Ok(e)
    }).is_err());
    let result=search::allocate_with_work(&p,deadline(),2,|v,w| {
        let mut e=complete(&p,v,w)?;e.slopes.fill(None);Ok(e)
    }).unwrap();
    assert_eq!(result.newton_trials,0);assert!(result.reason.is_none());
    assert!(result.passing.margin(&p).unwrap()<=0.0);
    assert!((result.values[0]-3.0).abs()<0.0001);
}

#[test]
fn optional_constraint_parsing_does_not_weaken_the_original_allocation_admission() {
    let base=J::parse(&BASE).unwrap();let mut spec=J::parse(SPEC).unwrap();
    input::put(&mut spec,"thermal_constraints",J::parse(r#"[{"name":"memory","component":"memory","temperature_limit_k":301}]"#).unwrap()).unwrap();
    let (p,limits)=parse_plan(&base,&spec).unwrap();assert_eq!(limits.len(),1);
    assert_eq!(p.planned_steps,14);assert_eq!(limits[0].objective.get("gradient"),Some(&J::Bool(false)));
    let request=p.request(&[1.0,2.0]).unwrap();let observed=limits[0].request(&request).unwrap();
    for key in ["solid","hydraulics","air","radiation","budgets"] {assert_eq!(request.get(key),observed.get(key));}
    assert_eq!(request.path(&["transient","intervals"]),observed.path(&["transient","intervals"]));
    input::put(&mut spec,"unknown_policy",J::Bool(true)).unwrap();
    assert!(parse_plan(&base,&spec).is_err());
}
