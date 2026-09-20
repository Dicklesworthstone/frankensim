use fs_couple::render::schedule::force::file::{ModalPerformance, MODAL_CONTACT_PERFORMANCE_SCHEMA};
use fs_couple::render::schedule::force::file::design::{
    DesignFileError, EquilibriumDesignFile, EQUILIBRIUM_DESIGN_HASH_DOMAIN,
};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, DesignError, DesignWork,
};
use fs_blake3::hash_domain;
use fs_exec::CancelGate;

const MODEL: &str = include_str!("../examples/equilibrium-design.model");
const FIT: &str = include_str!("../examples/equilibrium-design.fit");

// Physical-coordinate quadratic equilibrium, independent of the modal solver.
fn physical(k: f64, contact: f64, gap: f64, load: f64) -> [f64; 2] {
    let c = 1.0 / k + 1.0 / 10000.0;
    let excess = load / k - gap;
    let penetration = 2.0 * excess / (1.0 + (1.0 + 4.0 * contact * c * excess).sqrt());
    let reaction = contact * penetration * penetration;
    [(load - reaction) / k, reaction / 10000.0]
}
fn value(x: &[f64; 3]) -> f64 {
    [0.8, 1.6, 2.5].into_iter().map(|load| {
        let expected = physical(600.0, 1.8e8, 0.0002, load);
        let got = physical(400.0 + 400.0*x[0], 1e8 + 1e8*x[1], 0.0001 + 0.0002*x[2], load);
        got.iter().zip(expected).map(|(a,b)| 0.5*((a-b)/0.0002).powi(2)).sum::<f64>()
    }).sum()
}

#[test]
fn authored_model_loads_targets_and_scaled_fields_reach_the_actual_equilibrium_adjoint() {
    let gate = CancelGate::new();
    let loaded = EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), FIT.as_bytes(), &gate).unwrap();
    assert_eq!(loaded.model_info().schema, MODAL_CONTACT_PERFORMANCE_SCHEMA);
    assert_eq!(loaded.design_hash(), hash_domain(EQUILIBRIUM_DESIGN_HASH_DOMAIN, FIT.as_bytes()));
    let mut work = DesignControl::new(8, 24);
    for x in [[0.0,0.0,0.0], [0.15,0.1,0.2], [0.5,0.8,0.5]] {
        let result = loaded.problem().evaluate(&x, &mut work, &gate).unwrap();
        assert!((result.value-value(&x)).abs() < 1e-9);
        assert_eq!(result.cases.len(), 3);
        for j in 0..3 {
            let mut plus = x; let mut minus = x;
            plus[j] += 1e-5; minus[j] -= 1e-5;
            let fd = (value(&plus)-value(&minus))/2e-5;
            assert!((result.gradient[j]-fd).abs() < 1e-6*fd.abs().max(1e-3));
        }
    }
    assert_eq!(work.work(), DesignWork { evaluations:3, case_solves:9 });
}

#[test]
fn a_different_spring_only_model_and_shared_case_loads_are_not_a_hardcoded_rig() {
    let source = "frankensim-modal-performance-v2\nsample_rate_hz 24000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 1 100 1000 1000\ncompile_limits 0 4\nvoices 1\nvoice retain-state 2 2\nmode 20 0.1 0 0 0 0\nmode 30 0.1 0 0 0 0\nport 0 2 -1\nport 0 -0.5 1.5\ncoupling_limits 0 10000 0.9 1000 1000 1000 1e-10 1e-10 1e-8\nconnections 0\nevents 0\n";
    let fit = "frankensim-equilibrium-design-v1\npreload_limits 0 1 10000\nsensitivity_limits 0 10000 10000 0\ndesign_limits 2 1 2 3\ncases 2\ncase first 2 1\nload 0 0 1\nload 0 1 -2\ntarget 0 1 0.001 0.01 1\ncase second 2 1\nload 0 0 -0.5\nload 0 1 -2\ntarget 0 0 0.002 0.01 1\nvariables 1\nvariable common-N -2 0.5 -4 0 2\nbind actuator-force 0 1\nbind actuator-force 1 1\n";
    let loaded = EquilibriumDesignFile::from_bytes(source.as_bytes(), fit.as_bytes(), &CancelGate::new()).unwrap();
    let result = loaded.problem().evaluate(&[0.4], &mut DesignControl::new(1,2), &CancelGate::new()).unwrap();
    let forces = [[1.0,-1.8],[-0.5,-1.8]];
    let columns = [[2.0,-1.0],[-0.5,1.5]];
    let mut objective = 0.0; let mut gradient = 0.0;
    for i in 0..2 {
        let q = [(2.0*forces[i][0]-0.5*forces[i][1])/400.0,
                 (-forces[i][0]+1.5*forces[i][1])/900.0];
        let sensor = columns[1-i];
        let observed = sensor[0]*q[0]+sensor[1]*q[1];
        let error = (observed-[0.001,0.002][i])/0.01;
        objective += 0.5*error*error;
        gradient += error/0.01 * (sensor[0]*(-0.5)/400.0+sensor[1]*1.5/900.0)*0.5;
        assert!((result.cases[i].observations_m[0]-observed).abs() < 1e-12);
    }
    assert!((result.value-objective).abs() < 1e-10);
    assert!((result.gradient[0]-gradient).abs() < 1e-10);
}

#[test]
fn source_vibration_preload_forces_events_and_friction_are_not_silently_discarded() {
    for source in [
        MODEL.replace("mass 0.04 0 0", "mass 0.04 0 1"),
        MODEL.replace("port 0 5", "port 1 5"),
        MODEL.replace("voice free-mass", "voice free-mass-preload"),
        MODEL.replace("events 0", "events 1\nforce 0 0 0 1"),
        MODEL.replace("frankensim-modal-performance-v3", "frankensim-modal-performance-v1"),
    ] {
        assert!(matches!(EquilibriumDesignFile::from_bytes(source.as_bytes(), FIT.as_bytes(), &CancelGate::new()),
            Err(DesignFileError::Model(_))));
    }
    let friction_source = MODEL.replace("frankensim-modal-performance-v3", "frankensim-modal-performance-v5")
        .replace("contact_limits", "multi_contact_limits 4 128 16384\ncontacts 1\ncontact_limits")
        .replace("events 0", "frictions 1\nfriction none\nevents 0");
    assert!(matches!(EquilibriumDesignFile::from_bytes(friction_source.as_bytes(), FIT.as_bytes(), &CancelGate::new()),
        Err(DesignFileError::Model(_))));
    // Audio keeps its own admission: a moving input remains legal there.
    let moving = MODEL.replace("mass 0.04 0 0", "mass 0.04 0 1");
    assert!(ModalPerformance::from_bytes(moving.as_bytes(), 1).is_ok());
}

#[test]
fn complete_design_records_valid_references_and_original_joint_limits_are_required() {
    for text in [
        FIT.replace("cases 3", "cases 65"), FIT.replace("target 1 0", "target 1 9"),
        FIT.replace("load 0 0", "load 7 0"), FIT.replace("bind contact-gap 0", "bind damping 0"),
        FIT.replace("bind contact-gap 0", "bind contact-gap 0 extra"),
        FIT.replace("bind contact-gap 0", "bind contact-stiffness 0"),
        FIT.replace("1e-9", "NaN"), format!("{FIT}ignored\n"),
    ] { assert!(EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), text.as_bytes(), &CancelGate::new()).is_err()); }
    for (end, _) in FIT.match_indices('\n') {
        if end+1 < FIT.len() { assert!(EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), &FIT.as_bytes()[..end+1], &CancelGate::new()).is_err()); }
    }
    let multi = MODEL.replace("frankensim-modal-performance-v3", "frankensim-modal-performance-v4")
        .replace("contact_limits", "multi_contact_limits 4 64 16384\ncontacts 1\ncontact_limits");
    assert!(EquilibriumDesignFile::from_bytes(multi.as_bytes(), FIT.as_bytes(), &CancelGate::new()).is_err());
    let tighter = FIT.replace("preload_limits 4 128", "preload_limits 4 64");
    assert!(EquilibriumDesignFile::from_bytes(multi.as_bytes(), tighter.as_bytes(), &CancelGate::new()).is_ok());
}

#[test]
fn repeated_port_references_cannot_amplify_a_small_file_into_unbounded_shape_copies() {
    let modes = 256;
    let mut source = format!("frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 1 100 1000 1000\ncompile_limits 0 256\nvoices 1\nvoice retain-state {modes} 1\n");
    for i in 0..modes { source.push_str(&format!("mode {} 0.1 0 0 0 0\n",100+i)); }
    source.push_str(&format!("port 0 {}\n", vec!["1";modes].join(" ")));
    source.push_str("coupling_limits 0 10000 0.9 1000 1000 1000 1e-10 1e-10 1e-8\nconnections 0\nevents 0\n");
    let mut fit = String::from("frankensim-equilibrium-design-v1\npreload_limits 0 1 10000\nsensitivity_limits 0 10000 10000 0\ndesign_limits 1 1 1 1024\ncases 1\ncase many 1 257\nload 0 0 0\n");
    fit.push_str(&"target 0 0 0 1 1\n".repeat(257));
    fit.push_str("variables 1\nvariable force-N 0 1 -1 1 1\nbind actuator-force 0 0\n");
    let error = match EquilibriumDesignFile::from_bytes(source.as_bytes(), fit.as_bytes(), &CancelGate::new()) {
        Ok(_) => panic!("shape copy budget must refuse"), Err(e) => e,
    };
    assert!(error.to_string().contains("shape coefficients exceed"));
}

#[test]
fn cancellation_and_source_identity_survive_without_partial_evaluation() {
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), FIT.as_bytes(), &gate),
        Err(DesignFileError::Design(DesignError::Cancelled))));
    let a = EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), FIT.as_bytes(), &CancelGate::new()).unwrap();
    let b = EquilibriumDesignFile::from_bytes(MODEL.as_bytes(), FIT.as_bytes(), &CancelGate::new()).unwrap();
    assert_eq!(a.model_info(), b.model_info()); assert_eq!(a.design_hash(), b.design_hash());
    let mut work = DesignControl::new(1,3);
    let result = a.problem().evaluate(&[0.0;3], &mut work, &CancelGate::new()).unwrap();
    let again = b.problem().evaluate(&[0.0;3], &mut DesignControl::new(1,3), &CancelGate::new()).unwrap();
    assert_eq!(result, again);
}
