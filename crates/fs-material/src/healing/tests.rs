use super::*;
use super::strength::{StrengthKnot, StrengthMap};

fn table() -> MobilityTable {
    MobilityTable::new("synthetic-mobility-v1".into(), vec![
        MobilityKnot { temperature_k: 300.0, rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 400.0, rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 500.0, rate_per_s: 1.0 },
        MobilityKnot { temperature_k: 600.0, rate_per_s: 3.0 },
    ]).unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-13, "{a:e} != {b:e}"); }

#[test]
fn isothermal_dose_is_rate_times_actual_contact_duration() {
    let law = table();
    let virgin = ContactHistory::new(&law).unwrap();
    let hot = virgin.advance(500.0,500.0,0.4).unwrap();
    close(hot.dose(),0.4);
    close(hot.contact_time_s(),0.4);
    assert_eq!(virgin.dose(),0.0);
    assert_eq!(virgin.contact_time_s(),0.0);
    assert_eq!(hot.law_identity(),"synthetic-mobility-v1");
}

#[test]
fn knot_crossing_cooling_and_reheating_keep_the_old_dose() {
    let law = table();
    // Over 3 s, spend one second per 100 K interval. Areas: 0 + 1/2 + 2.
    close(law.dose_increment(300.0,600.0,3.0).unwrap(),2.5);
    close(law.dose_increment(600.0,300.0,3.0).unwrap(),2.5);
    let warm = ContactHistory::new(&law).unwrap().advance(450.0,550.0,0.5).unwrap();
    let cold = warm.advance(350.0,350.0,9.0).unwrap();
    assert_eq!(cold.dose().to_bits(),warm.dose().to_bits());
    let reheated = cold.advance(450.0,550.0,0.5).unwrap();
    close(reheated.dose(),2.0*warm.dose());
    close(reheated.contact_time_s(),10.0);
}

#[test]
fn linear_path_partitioning_preserves_exposure() {
    let law = table();
    for &(start,end) in &[(301.0,599.0),(599.0,301.0),(425.0,475.0),(500.0,500.0)] {
        let whole = law.dose_increment(start,end,2.0).unwrap();
        for count in [2,7,31] {
            let mut dose = 0.0;
            for i in 0..count {
                let a = start+(end-start)*f64::from(i)/f64::from(count);
                let b = start+(end-start)*f64::from(i+1)/f64::from(count);
                dose += law.dose_increment(a,b,2.0/f64::from(count)).unwrap();
            }
            close(dose,whole);
        }
    }
}

#[test]
fn bad_or_out_of_domain_paths_never_change_accepted_history() {
    let law = table();
    let state = ContactHistory::new(&law).unwrap().advance(500.0,500.0,0.2).unwrap();
    for (a,b,dt) in [(299.0,450.0,1.0),(500.0,601.0,1.0),
        (500.0,500.0,-1.0),(f64::NAN,500.0,1.0),(500.0,500.0,f64::INFINITY)] {
        assert!(state.advance(a,b,dt).is_err());
        assert_eq!(state.dose(),0.2);
    }
    assert!(MobilityTable::new("".into(),law.knots().to_vec()).is_err());
    let mut bad = law.knots().to_vec();
    bad[1].rate_per_s = -1.0;
    assert!(MobilityTable::new("bad".into(),bad).is_err());
}

fn strength_map() -> StrengthMap {
    StrengthMap::new("synthetic-strength-v1".into(),"synthetic-mobility-v1".into(),vec![
        StrengthKnot { dose: 0.0, peak_traction_pa: 0.0, fracture_energy_j_m2: 0.0 },
        StrengthKnot { dose: 1.0, peak_traction_pa: 2e5, fracture_energy_j_m2: 100.0 },
    ],1e9,2e9,0.01).unwrap()
}

#[test]
fn strength_requires_an_explicit_compatible_map_and_nonzero_exposure() {
    let law = table();
    let map = strength_map();
    let virgin = ContactHistory::new(&law).unwrap();
    assert!(map.freeze(&virgin).unwrap().law.is_none());
    let exposed = virgin.advance(500.0,500.0,0.5).unwrap();
    let frozen = map.freeze(&exposed).unwrap();
    let cohesive = frozen.law.unwrap();
    assert_eq!(cohesive.strength_pa(),1e5);
    assert_eq!(cohesive.fracture_energy_j_m2(),50.0);
    let step = cohesive.bonded_point().trial([cohesive.failure_separation_m(),0.0,0.0],1e-8).unwrap();
    assert!(step.point.failed());
    close(step.energy.dissipated,50.0);
    assert!(map.freeze(&virgin.advance(500.0,500.0,2.0).unwrap()).is_err());
    let foreign = MobilityTable::new("another-material".into(),law.knots().to_vec()).unwrap();
    assert!(map.freeze(&ContactHistory::new(&foreign).unwrap()).is_err());
}

#[test]
fn invalid_custom_increment_is_refused_and_restart_is_exact() {
    #[derive(Debug)]
    struct Broken;
    impl HealingLaw for Broken {
        fn identity(&self) -> &str { "broken" }
        fn dose_increment(&self,_a: f64,_b: f64,_dt: f64) -> Result<f64,HealingError> { Ok(f64::NAN) }
    }
    assert!(ContactHistory::new(&Broken).unwrap().advance(500.0,500.0,1.0).is_err());
    let law = table();
    let mut a = ContactHistory::new(&law).unwrap();
    for _ in 0..7 { a = a.advance(450.0,550.0,0.1).unwrap(); }
    let mut b = a.clone();
    for _ in 0..9 {
        a = a.advance(550.0,450.0,0.1).unwrap();
        b = b.advance(550.0,450.0,0.1).unwrap();
    }
    assert_eq!(a.dose().to_bits(),b.dose().to_bits());
    assert_eq!(a.contact_time_s().to_bits(),b.contact_time_s().to_bits());
}
