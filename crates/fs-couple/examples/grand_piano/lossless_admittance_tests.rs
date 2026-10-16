use super::*;

#[test]
fn lossless_selection_is_explicit_admittance_only_and_preserves_value_boundaries() {
    let args=|s:&str|s.split_whitespace().map(str::to_owned).collect::<Vec<_>>();
    assert!(admittance_options(&[]).unwrap().1);
    let (options,damping)=admittance_options(&args("--modes 12 --lossless-structure --substeps 8")).unwrap();
    assert!(!damping);assert_eq!(options.modes,12);assert_eq!(options.substeps,8);
    for input in ["--lossless-structure --lossless-structure","--lossless-structure 1",
        "--modes --lossless-structure 12","--substeps --lossless-structure",
        "--lossless-structure --note 69","--lossless-structure --string-stretching axial.fsps"] {
        assert!(admittance_options(&args(input)).is_err(),"{input}");
    }
    for command in ["response","render","render-loaded"] {
        let mut command=args(&format!("{command} missing.fsb missing.csv missing.obj missing.fspe unused.wav"));
        if command[0]!="response" {command.push("0.1".into());}
        command.push("--lossless-structure".into());
        let error=run(&command).unwrap_err();assert!(error.contains("--lossless-structure"),"{error}");
    }
}

#[test]
fn conservative_structure_sweeps_the_actual_partial_with_bem_radiation_still_present() {
    let (board,mut courses,obj,mut spec)=tests::small_source_inputs();
    courses[0].unison=1;courses[0].detune_cents=0.0;courses[0].duplex_length_m=0.0;
    // Lengthen the physical test string so its actual partial stays inside
    // this small fixture's resolved BEM wavelength band; do not shift the solve.
    courses[0].length_m*=2.0;
    let hz=courses[0].partial_hz(1,courses[0].tension_n);
    // The center sample is the physical fixed-interface partial, including
    // floating-point grid roundoff. A supplied force cannot shift it away.
    spec.band_hz=(0.98*hz,1.02*hz);spec.frequencies=17;
    let (options,loss)=admittance_options(&["--lossless-structure".into(),"--modes".into(),"4".into()]).unwrap();
    let conservative=admittance_controlled(&board,&courses,&obj,&spec,69,&options,loss).unwrap();
    let damped=admittance_controlled(&board,&courses,&obj,&spec,69,&options,true).unwrap();
    assert!(conservative.contains("explicitly disabled by --lossless-structure"));
    let mut radiates=false;let mut rows=0;
    for row in conservative.lines().filter(|l|!l.starts_with('#')).skip(1) {
        let f=row.split(',').collect::<Vec<_>>();assert_eq!(f.len(),15);
        let values=f[3..].iter().map(|s|s.parse::<f64>().unwrap()).collect::<Vec<_>>();
        assert!(values.iter().all(|v|v.is_finite()));
        assert_eq!(values[5],0.0);assert_eq!(values[6],0.0); // wood, string watts
        assert!((values[4]-values[7]).abs()<1e-10+1e-7*values[4].abs()); // input = air loss
        radiates|=values[7]>0.0;rows+=1;
    }
    assert_eq!(rows,17*(1+spec.receivers.len()));assert!(radiates);
    assert!(damped.lines().filter(|l|!l.starts_with('#')).skip(1).any(|line| {
        let f=line.split(',').collect::<Vec<_>>();f[8].parse::<f64>().unwrap()>0.0 || f[9].parse::<f64>().unwrap()>0.0
    }));
}
