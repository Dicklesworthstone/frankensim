//! Real film reduction consumed by the existing nonlinear time owner.
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactError, ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::membrane::MembranePotential;
use fs_exec::CancelGate;
use fs_plate::ModePair;
use fs_plate::shell::head::{TensionedDisk, TensionedDiskSpec};
use fs_plate::shell::head::nonlinear::MembraneReductionBudget;
use fs_plate::shell::profile::ProfileBudget;

fn film() -> (TensionedDisk, Vec<ModePair>, Vec<usize>) {
    let d = TensionedDisk::new(TensionedDiskSpec { radius_m: 0.17,
        thickness_m: 0.000254, young_pa: 4e9, poisson: 0.38, density_kg_m3: 1390.0,
        tension_n_m: 3000.0, radial_intervals: 2, azimuths: 8 },
        ProfileBudget { max_nodes: 100, max_triangles: 100, max_feature_evaluations: 0 }).unwrap();
    let n = d.model.free;
    let mut k = vec![0.0; n*n]; let mut m = k.clone();
    for i in 0..n { for j in 0..n { k[i*n+j] = d.model.k.get(i,j); m[i*n+j] = d.model.m.get(i,j); }}
    let mut modes = fs_modal::eigh_gen_dense(&k, &m, n).unwrap();
    // Explicit single-mode experiment; no full-band adequacy claim.
    modes.truncate(1);
    let rim = (0..d.mesh.nodes.len()).filter(|&i| d.model.dof_map[3*i].is_none()).collect();
    (d, modes, rim)
}
fn reduced(limit: f64) -> MembranePotential {
    let (d, modes, rim) = film();
    MembranePotential::from_pencil(&d.mesh, &d.section, &d.model, &modes, &rim,
        MembraneReductionBudget { max_modes: 4, max_nodes: 100, max_facet_pairs: 1000,
            max_solve_entries: 100000, relative_tolerance: 1e-6 }, limit).unwrap()
}
fn body(potential: MembranePotential, q: f64, v: f64) -> ImpactBody {
    ImpactBody { potential: BodyPotential::Membrane(potential),
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }],
        damping_per_s: vec![0.0] }
}
fn config() -> ImpactConfig {
    ImpactConfig { dt_s: 1.0/192000.0, max_steps: 20000, maximum_energy_j: 100.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 1e6 }
}
// Independent continuous Duffing period. k and beta come from the ACTUAL
// reduced film; this quadrature does not use the numerical time integrator.
fn quarter_period(k: f64, beta: f64, amplitude: f64) -> f64 {
    let count = 1024;
    let h = std::f64::consts::FRAC_PI_2 / f64::from(count);
    let mut sum = 0.0;
    for i in 0..=count {
        let theta = f64::from(i)*h;
        let weight = if i == 0 || i == count { 1.0 } else if i%2 == 0 { 2.0 } else { 4.0 };
        sum += weight / (k + 0.5*beta*amplitude*amplitude*(1.0 + theta.cos().powi(2))).sqrt();
    }
    sum*h/3.0
}

#[test]
fn real_film_strong_and_weak_release_follow_the_geometry_derived_nonlinear_period() {
    let material = reduced(0.25);
    let k = material.reduction().omegas()[0].powi(2);
    let beta = 4.0*material.reduction().stretching_energy(&[1.0]);
    let slope_per_q = material.reduction().maximum_slope(&[1.0]);
    let mut frequencies = Vec::new();
    for slope in [0.01, 0.10] {
        let amplitude = slope / slope_per_q;
        let expected = quarter_period(k, beta, amplitude);
        let mut s = ImpactSystem::new(vec![body(material.clone(), amplitude, 0.0)],
            vec![], vec![], vec![], config()).unwrap();
        let initial = s.stored_energy_j();
        let observation = s.membrane_observation(0).unwrap();
        assert!((observation.maximum_slope-slope).abs() < 1e-13);
        assert!(observation.stretching_energy_j > 0.0);
        let gate = CancelGate::new_clock_free();
        let mut old = amplitude; let mut crossing = None;
        for _ in 0..10000 {
            let f = s.step(&[0.0], &gate).unwrap();
            let q = s.state()[0];
            if old > 0.0 && q <= 0.0 {
                crossing = Some(f.time_s - config().dt_s + config().dt_s * old/(old-q));
                break;
            }
            old = q;
        }
        let actual = crossing.expect("positive release must reach its first zero");
        assert!((actual-expected).abs() < 0.005*expected, "{actual} != {expected}");
        assert!((s.stored_energy_j()-initial).abs() < 1e-7*initial + 1e-10);
        frequencies.push(0.25/actual);
    }
    assert!(frequencies[1] > 1.01*frequencies[0], "strong-hit frequency must change through strain, not an audio control");
}

#[test]
fn sealed_air_reciprocally_loads_a_nonlinear_head_and_rest_is_silent() {
    let material = reduced(0.25);
    let make = |connected| ImpactSystem::new(vec![body(material.clone(), 3e-5, 0.0), body(material.clone(), 0.0, 0.0)],
        vec![], vec![], if connected { vec![VolumeSpring { bulk_modulus_pa: 1.4e5,
            volume_m3: 0.01, areas: vec![0.1,-0.1] }] } else { vec![] }, config()).unwrap();
    let mut a = make(true); let mut off = make(false);
    let gate = CancelGate::new_clock_free();
    for _ in 0..128 {
        a.step(&[0.0;2], &gate).unwrap(); off.step(&[0.0;2], &gate).unwrap();
        assert_eq!(off.state()[2], 0.0); assert_eq!(off.state()[3], 0.0);
    }
    assert!(a.state()[2].abs() > 1e-9);
    assert!(a.membrane_observation(1).unwrap().stretching_energy_j > 0.0);
    assert_ne!(a.state()[0].to_bits(), off.state()[0].to_bits());
}

#[test]
fn failed_slope_admission_and_cancellation_do_not_publish_motion() {
    let strict = reduced(1e-7);
    let excessive = 1e-6/strict.reduction().maximum_slope(&[1.0]);
    assert!(ImpactSystem::new(vec![body(strict.clone(), excessive, 0.0)], vec![], vec![], vec![], config()).is_err());
    let mut s = ImpactSystem::new(vec![body(strict, 0.0, 1.0)], vec![], vec![], vec![], config()).unwrap();
    let state = s.state().to_vec(); let energy = s.stored_energy_j();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(s.step(&[0.0], &cancelled), Err(ImpactError::Cancelled)));
    assert!(s.step(&[0.0], &CancelGate::new_clock_free()).is_err());
    assert_eq!(s.state(), state); assert_eq!(s.samples(), 0); assert_eq!(s.stored_energy_j(), energy);
    assert_eq!(s.membrane_observation(0).unwrap().maximum_slope, 0.0);
}

#[test]
fn membrane_admission_keeps_the_original_mode_and_rejects_unphysical_limits() {
    let (disk,modes,rim) = film();
    let budget = MembraneReductionBudget { max_modes: 4,max_nodes: 100,max_facet_pairs: 1000,
        max_solve_entries: 100000,relative_tolerance: 1e-6 };
    for limit in [0.0,1.0,f64::NAN,f64::INFINITY] {
        assert!(MembranePotential::from_pencil(&disk.mesh,&disk.section,&disk.model,&modes,&rim,budget,limit).is_err());
    }
    let p = MembranePotential::from_pencil(&disk.mesh,&disk.section,&disk.model,&modes,&rim,budget,0.2).unwrap();
    assert_eq!(p.reduction().omegas()[0].to_bits(),modes[0].lambda.sqrt().to_bits());
    assert!(p.observe(&[]).is_err()); assert!(p.observe(&[f64::NAN]).is_err());
    assert_eq!(p.observe(&[0.0]).unwrap().stretching_energy_j,0.0);
}
