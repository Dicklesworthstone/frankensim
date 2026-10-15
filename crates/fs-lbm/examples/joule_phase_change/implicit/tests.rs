use super::*;

fn mode()->IntegrationMode {IntegrationMode::Implicit {difference_step_j:1e-10}}
fn run(table:&Table,contact:Option<&ContactTable>,dt:f64)->Vec<Row> {
    let mut rows=Vec::new();
    simulate_with_mode(table,contact,10.0,60.0,dt,mode(),|r|rows.push(r)).unwrap();
    rows
}

#[test]
fn cli_suffix_keeps_legacy_arguments_and_requires_explicit_positive_joule_scale() {
    let strings=|s:&[&str]|s.iter().map(|s|s.to_string()).collect::<Vec<_>>();
    let mut legacy=strings(&["example","bulk.txt","10","60","0.25","contact.txt"]);
    let original=legacy.clone();
    assert_eq!(IntegrationMode::parse_suffix(&mut legacy).unwrap(),IntegrationMode::Explicit);
    assert_eq!(legacy,original);
    let mut selected=original.clone();selected.extend(strings(&["--implicit","1e-10"]));
    assert_eq!(IntegrationMode::parse_suffix(&mut selected).unwrap(),mode());
    assert_eq!(selected,original);
    for arguments in [
        vec!["example","bulk.txt","--implicit"],
        vec!["example","bulk.txt","--implicit","0"],
        vec!["example","bulk.txt","--implicit","NaN"],
        vec!["example","bulk.txt","--implicit","-1"],
        vec!["example","bulk.txt","--implicit","1e-10","contact.txt"],
        vec!["example","bulk.txt","--implicit","1e-10","--implicit","1e-10"],
    ] {assert!(IntegrationMode::parse_suffix(&mut strings(&arguments)).is_err());}
}

#[test]
fn data_only_latent_change_reaches_implicit_melting_and_retains_power_off_history() {
    let a=Table::parse(include_str!("../latent_20k.txt")).unwrap();
    let b=Table::parse(include_str!("../latent_40k.txt")).unwrap();
    let a=run(&a,None,0.25);let b=run(&b,None,0.25);
    assert_eq!(a.len(),60);assert_eq!(b.len(),60);
    assert_eq!(a[29].time,30.0);assert_eq!(a[59].time,60.0);
    assert_eq!(a[59].temperature,310.0);assert_eq!(b[59].temperature,310.0);
    assert!(a[59].liquid_fraction>b[59].liquid_fraction+0.1);
    assert!(a[29].current<b[29].current);
    for rows in [&a,&b] {
        assert!(rows[59].stored>0.0 && rows[59].enthalpy>0.0);
        assert!(rows[30..].iter().all(|r|r.current==0.0 && r.electrode==rows[29].electrode));
        assert!(rows.iter().all(|r|r.maximum_cell_residual_j<=1e-8 && r.residual.abs()<1e-6));
        assert!(rows[59].newton_iterations>0 && rows[59].solves>rows[59].newton_iterations);
    }
}

#[test]
fn source_card_contact_changes_implicit_current_and_phase_without_changing_bulk() {
    let table=Table::parse(include_str!("../transport_phase_dependent.txt")).unwrap();
    let low=ContactTable::parse(include_str!("../contact_low.txt"),&table).unwrap();
    let high=ContactTable::parse(include_str!("../contact_high.txt"),&table).unwrap();
    let a=run(&table,Some(&low),0.25);let b=run(&table,Some(&high),0.25);
    assert_eq!(a[59].temperature,310.0);assert_eq!(b[59].temperature,310.0);
    assert!(a[59].liquid_fraction>b[59].liquid_fraction+0.2);
    assert!(a[29].current>b[29].current);
    for rows in [&a,&b] {
        assert!(rows[20].contact_resistance.unwrap()>rows[10].contact_resistance.unwrap());
        assert!(rows[30..].iter().all(|r|r.current==0.0 && r.electrode==rows[29].electrode));
        assert!(rows.iter().all(|r|r.residual.abs()<1e-6));
    }
    // Refine the actual implicit path toward an independently refined explicit
    // trajectory of the same physical owners; neither method is the exact model.
    let mut explicit=Vec::new();
    simulate_with_contact(&table,Some(&low),10.0,60.0,0.00390625,|r|explicit.push(r)).unwrap();
    let finer=run(&table,Some(&low),0.125);
    let coarse_error=(a[59].enthalpy-explicit[59].enthalpy).abs();
    let fine_error=(finer[59].enthalpy-explicit[59].enthalpy).abs();
    assert!(fine_error<0.7*coarse_error,"{fine_error} !< 0.7*{coarse_error}");
}

#[test]
fn exact_unpowered_cold_state_and_explicit_mode_keep_their_original_meaning() {
    let table=Table::parse(include_str!("../latent_20k.txt")).unwrap();
    let mut rows=Vec::new();
    let cold=simulate_with_mode(&table,None,0.0,60.0,0.25,mode(),|r|rows.push(r)).unwrap();
    assert_eq!(cold.temperature,300.0);assert_eq!(cold.enthalpy,0.0);
    assert_eq!(cold.stored,0.0);assert_eq!(cold.electrode,0.0);
    assert_eq!(cold.newton_iterations,0);
    assert!(rows.iter().all(|r|r.maximum_cell_residual_j==0.0));
    let legacy=simulate(&table,10.0,60.0,0.25,|_|{}).unwrap();
    let selected=simulate_with_mode(&table,None,10.0,60.0,0.25,IntegrationMode::Explicit,|_|{}).unwrap();
    assert_eq!(legacy.enthalpy,selected.enthalpy);
    assert_eq!(legacy.electrode,selected.electrode);
    assert_eq!(legacy.solves,selected.solves);
}

#[test]
fn adaptive_suffix_requires_explicit_error_and_retry_floor_without_changing_data_arguments() {
    let strings=|s:&[&str]|s.iter().map(|s|s.to_string()).collect::<Vec<_>>();
    let original=strings(&["example","bulk.txt","10","60","0.5","contact.txt"]);
    let mut args=original.clone();
    args.extend(strings(&["--adaptive","1e-10","0.03","1e-7"]));
    let selected=IntegrationMode::parse_suffix(&mut args).unwrap();
    assert_eq!(selected,IntegrationMode::Adaptive {
        difference_step_j:1e-10,cell_error_budget_j:0.03,minimum_step_s:1e-7,
    });
    assert!(selected.is_adaptive() && selected.is_implicit());assert_eq!(args,original);
    for bad in [
        vec!["--adaptive","1e-10","0.03"],
        vec!["--adaptive","1e-10","0","1e-7"],
        vec!["--adaptive","1e-10","0.03","NaN"],
        vec!["--adaptive","0","0.03","1e-7"],
        vec!["--adaptive","1e-10","0.03","1e-7","--implicit","1e-10"],
    ] {
        let mut args=original.clone();args.extend(strings(&bad));let unchanged=args.clone();
        assert!(IntegrationMode::parse_suffix(&mut args).is_err());assert_eq!(args,unchanged);
    }
}

#[test]
fn adaptive_material_files_change_phase_and_keep_one_whole_run_error_budget() {
    let t=Table::parse(include_str!("../transport_phase_dependent.txt")).unwrap();
    let mode=IntegrationMode::Adaptive {
        difference_step_j:1e-10,cell_error_budget_j:0.03,minimum_step_s:1e-7,
    };
    let mut outputs=Vec::new();
    for source in [include_str!("../contact_low.txt"),include_str!("../contact_high.txt")] {
        let contact=ContactTable::parse(source,&t).unwrap();let mut rows=Vec::new();
        simulate_with_mode(&t,Some(&contact),10.0,60.0,0.5,mode,|r|rows.push(r)).unwrap();
        assert_eq!(rows.len(),60);assert_eq!(rows[29].time,30.0);assert_eq!(rows[59].time,60.0);
        assert_eq!(rows[59].temperature,310.0);assert!(rows[59].stored>0.0);
        assert!(rows.iter().any(|r|r.rejected_intervals>0));
        assert!(rows.iter().all(|r|r.residual.abs()<1e-6 && r.maximum_cell_residual_j<=1e-8));
        assert!(rows.iter().map(|r|r.temporal_indicator_j).sum::<f64>()<=0.03*(1.0+1e-10));
        assert!(rows[30..].iter().all(|r|r.current==0.0 && r.electrode==rows[29].electrode));
        assert!(rows[59].newton_iterations>0 && rows[59].solves>rows[59].newton_iterations);
        outputs.push(rows);
    }
    // Independently solved three-volume equations, not named-metal measurements.
    assert!((outputs[0][59].liquid_fraction-0.60984188).abs()<5e-4);
    assert!((outputs[1][59].liquid_fraction-0.32909529).abs()<5e-4);
    assert!(outputs[0][29].current>outputs[1][29].current);
}

#[test]
fn adaptive_cold_equilibrium_is_exact_and_invalid_floor_emits_no_observations() {
    let t=Table::parse(include_str!("../latent_20k.txt")).unwrap();
    let mode=IntegrationMode::Adaptive {
        difference_step_j:1e-10,cell_error_budget_j:0.03,minimum_step_s:1e-7,
    };
    let last=simulate_with_mode(&t,None,0.0,1.0,0.25,mode,|r| {
        assert_eq!(r.temporal_indicator_j,0.0);assert_eq!(r.rejected_intervals,0);
    }).unwrap();
    assert_eq!(last.stored,0.0);assert_eq!(last.enthalpy,0.0);assert_eq!(last.newton_iterations,0);
    let bad=IntegrationMode::Adaptive {
        difference_step_j:1e-10,cell_error_budget_j:0.03,minimum_step_s:1.0,
    };
    let mut emitted=0;
    assert!(simulate_with_mode(&t,None,10.0,1.0,0.25,bad,|_|emitted+=1).is_err());
    assert_eq!(emitted,0);
}
