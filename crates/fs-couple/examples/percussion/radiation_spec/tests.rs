use super::*;

const VALID: &str = "frankensim-radiation-preparation-v1\nband_hz,40,6000\ntraining_intervals,64\nmax_order,24\nsubdivisions,1\nmax_panels,4096\nmax_dense_work,20000000000000\n";

#[test]
fn the_declared_lattice_partitions_training_selection_and_audit_without_overlap() {
    let s=Spec::parse(VALID).unwrap();let w=s.frequencies().unwrap();
    assert_eq!(w.len(),257);assert_eq!(w[0],40.0*core::f64::consts::TAU);
    assert_eq!(*w.last().unwrap(),6000.0*core::f64::consts::TAU);
    assert_eq!(w.iter().step_by(4).count(),65);
    assert_eq!(w.iter().skip(2).step_by(4).count(),64);
    assert_eq!(w.iter().skip(1).step_by(2).count(),128);
    for command in ["splash-mic","splash-wav","drum-modal-mic","snare-off-mic"] {s.admit_command(command).unwrap();}
    for command in ["splash","snare","export-shell-mesh","unknown"] {assert!(s.admit_command(command).is_err());}
}

#[test]
fn complete_input_and_work_bounds_refuse_without_hidden_defaults() {
    for (a,b) in [("band_hz,40,6000","band_hz,40,24000"),
        ("band_hz,40,6000","band_hz,NaN,6000"),("band_hz,40,6000","band_hz,0,6000"),
        ("max_order,24","max_order,25"),("training_intervals,64","training_intervals,32"),
        ("subdivisions,1","subdivisions,5"),("max_panels,4096","max_panels,8192"),
        ("max_dense_work,20000000000000","max_dense_work,0")] {
        assert!(Spec::parse(&VALID.replace(a,b)).is_err());
    }
    for line in VALID.lines().skip(1) {
        assert!(Spec::parse(&VALID.replace(&format!("{line}\n"),"")).is_err());
        assert!(Spec::parse(&format!("{VALID}{line}\n")).is_err());
    }
    assert!(Spec::parse(&format!("{VALID}gain,2\n")).is_err());
    assert!(Spec::parse(&" ".repeat(MAX_BYTES+1)).is_err());
    let s=Spec::parse(VALID).unwrap();
    let (n,work)=s.work(100,2,2).unwrap();assert_eq!(n,400);
    assert_eq!(work,257*(400_u64.pow(3)+4*400_u64.pow(2)+4*400));
    assert!(s.work(1025,2,2).is_err());assert!(s.work(usize::MAX,2,2).is_err());
    assert!(s.work(4,MAX_INPUTS+1,2).is_err());assert!(s.work(4,1,3).is_err());
    assert!(Spec {max_dense_work:work-1,..s}.work(100,2,2).is_err());
    assert_eq!(Spec {max_dense_work:work,..s}.work(100,2,2).unwrap().1,work);
}

#[test]
fn option_reads_one_complete_file_and_leaves_arguments_unchanged_on_refusal() {
    let path=std::env::temp_dir().join(format!("frankensim-radiation-spec-{}.fra",std::process::id()));
    // Never overwrite a coincident file, including after an interrupted test run.
    let mut file=std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    use std::io::Write;file.write_all(VALID.as_bytes()).unwrap();drop(file);
    let name=path.to_str().unwrap();
    let mut args=vec!["splash-mic".into(),"64".into(),"--radiation-spec".into(),name.into(),"--analytic-newton".into()];
    assert_eq!(option(&mut args).unwrap(),Some(Spec::parse(VALID).unwrap()));
    assert_eq!(args,["splash-mic","64","--analytic-newton"]);
    for mut args in [vec!["--radiation-spec".into()],
        vec!["--radiation-spec".into(),name.into(),"--radiation-spec".into(),name.into()]] {
        let old=args.clone();assert!(option(&mut args).is_err());assert_eq!(args,old);
    }
    std::fs::remove_file(&path).unwrap();
    let mut args=vec!["--radiation-spec".into(),name.into()];let old=args.clone();
    assert!(option(&mut args).is_err());assert_eq!(args,old);
}
