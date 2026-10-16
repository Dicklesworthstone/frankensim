use super::*;
use crate::bridge::tests::{birth, config as beam_config, material, temperatures};
use fs_solver::orthant::GramOrthantConfig;

fn config() -> Config {
    Config { max_bodies: 8, energy_tolerance_j: 1e-9,
        contact: supports::Config { contact: GramOrthantConfig {
            max_variables: 8, max_columns: 256, max_setup_products: 16_384,
            max_sweeps: 4096, gradient_tolerance: 1e-13 },
            gap_tolerance_m: 1e-11, complementarity_tolerance_j: 1e-13 } }
}
fn body<'s, 'm>(id: u64, beam: &'s Bridge<'m>, t: &'s BTreeMap<(u64,u64),f64>,
    extra: &'s [f64]) -> Body<'s, 'm> {
    Body { id, beam, temperatures_k: t, gravity_m_s2: 9.81,
        additional_downward_n_m: extra, config: beam_config() }
}
fn pair(id: u64, upper: u64, lower: u64, gap: f64) -> Contact {
    Contact { id, upper: Station { body: upper, reference_x_m: 0.01 },
        lower: Station { body: lower, reference_x_m: 0.01 }, reference_gap_m: gap }
}
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a:e} != {b:e}"); }
fn w(b: &Bridge<'_>) -> f64 { b.deformed_position_m(0.01, 0.0).unwrap()[1] }
fn solve<'m>(bodies: &[Body<'_, 'm>], contacts: &[Contact]) -> Solution<'m> {
    advance_controlled(bodies,contacts,1.0,config(),&mut || false).unwrap()
}

// G1: two independently flexible simply-supported beams, not a rigid-stop limit.
#[test]
fn two_body_reaction_and_both_deflections_refine_to_continuum_solution() {
    let upper_law = material(1e7, vec![]);
    let lower_law = material(1e8, vec![]);
    let q = 1200.0*1e-6*9.81;
    let l = 0.02_f64;
    let inertia = 1e-6*0.001_f64.powi(2)/12.0;
    let cu = l.powi(3)/(48.0*1e7*inertia);
    let cl = l.powi(3)/(48.0*1e8*inertia);
    let wu = -5.0*q*l.powi(4)/(384.0*1e7*inertia);
    let wl = -5.0*q*l.powi(4)/(384.0*1e8*inertia);
    let gap = 1e-6;
    let force = -(gap+wu-wl)/(cu+cl);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let upper = birth(n, Support::SimplySupported, &upper_law);
        let lower = birth(n, Support::SimplySupported, &lower_law);
        let tu = temperatures(&upper, 300.0); let tl = temperatures(&lower, 300.0);
        let zeros = vec![0.0; n as usize];
        let solved = advance_controlled(&[body(1,&upper,&tu,&zeros),body(2,&lower,&tl,&zeros)],
            &[pair(9,1,2,gap)], 1.0, config(), &mut || false).unwrap();
        let r = solved.report.reactions[0];
        errors.push((r.force_on_upper_n-force).abs());
        near(w(&solved.bodies[&1]), wu+cu*force, 2e-7/(f64::from(n)/8.0).powi(2));
        near(w(&solved.bodies[&2]), wl-cl*force, 2e-7/(f64::from(n)/8.0).powi(2));
        near(r.gap_m, 0.0, 1e-12);
        near(solved.report.beams[&1].vertical_reactions_n.iter().sum::<f64>()+r.force_on_upper_n, q*l, 1e-12);
        near(solved.report.beams[&2].vertical_reactions_n.iter().sum::<f64>()-r.force_on_upper_n, q*l, 1e-12);
        near(solved.report.contact_work_j, -r.force_on_upper_n*gap, 1e-17);
        near(solved.report.energy_residual_j, 0.0, 1e-17);
    }
    assert!(errors[0]/errors[1] > 3.9 && errors[1]/errors[2] > 3.9, "{errors:?}");
}

// G0/G3: the middle body appears with opposite signs in the two contact rows.
#[test]
fn three_road_stack_shares_load_without_losing_action_reaction_or_ordering() {
    let law = material(1e8, vec![]);
    let a = birth(12, Support::SimplySupported, &law);
    let b = a.clone(); let c = a.clone();
    let t = temperatures(&a,300.0);
    let qa = [0.06;12]; let qb = [0.02;12]; let qc = [0.0;12];
    let inputs = [body(1,&a,&t,&qa),body(2,&b,&t,&qb),body(3,&c,&t,&qc)];
    let pairs = [pair(10,1,2,0.0),pair(20,2,3,0.0)];
    let solved = solve(&inputs,&pairs);
    let reversed = solve(&[inputs[2],inputs[0],inputs[1]],&[pairs[1],pairs[0]]);
    assert_eq!(solved.bodies,reversed.bodies);
    assert_eq!(solved.report.reactions,reversed.report.reactions);
    assert!(solved.report.reactions.iter().all(|r| r.force_on_upper_n>0.0));
    let free = solve(&inputs,&[]);
    let mean = free.bodies.values().map(w).sum::<f64>()/3.0;
    for beam in solved.bodies.values() { near(w(beam),mean,3e-13); }
    let total_reaction: f64 = solved.report.beams.values().flat_map(|r| r.vertical_reactions_n).sum();
    let external = 3.0*a.mass_kg()*9.81+(0.06+0.02)*0.02;
    near(total_reaction,external,1e-12);
    near(solved.report.energy_j.mechanical_work,
        solved.report.external_work_j+solved.report.contact_work_j,1e-17);
}

#[test]
fn different_meshes_and_contact_removal_keep_independent_material_states() {
    let upper_law = material(1e7,vec![]); let lower_law = material(1e8,vec![]);
    let u = birth(7,Support::CantileverLeft,&upper_law);
    let l = birth(13,Support::SimplySupported,&lower_law);
    let tu = temperatures(&u,330.0); let tl = temperatures(&l,390.0);
    let inputs = [body(1,&u,&tu,&[0.0;7]),body(2,&l,&tl,&[0.0;13])];
    let free = advance_controlled(&inputs,&[],1.0,config(),&mut || false).unwrap();
    for input in inputs {
        let direct = input.beam.advance_controlled(1.0,input.temperatures_k,9.81,
            input.additional_downward_n_m,beam_config(),&mut || false).unwrap();
        assert_eq!(free.bodies[&input.id],direct.0);
    }
    let held = advance_controlled(&inputs,&[pair(1,1,2,0.0)],1.0,config(),&mut || false).unwrap();
    assert!(held.report.reactions[0].force_on_upper_n>0.0);
    assert!(w(&held.bodies[&2]) < w(&free.bodies[&2]));
    let after = advance_controlled(&[body(1,&held.bodies[&1],&tu,&[0.0;7]),
        body(2,&held.bodies[&2],&tl,&[0.0;13])],&[],2.0,config(),&mut || false).unwrap();
    assert!(after.report.reactions.is_empty());
    // Pure elasticity returns to the free shape; the actual temperature is retained.
    near(w(&after.bodies[&1]),w(&free.bodies[&1]),1e-12);
    near(w(&after.bodies[&2]),w(&free.bodies[&2]),1e-12);
}

// G1/G5: creep changes which neighbour is carrying load; no virgin-history replay.
#[test]
fn creep_engagement_nonisothermal_replay_and_lift_off_preserve_memory() {
    let soft = material(1e7,vec![(9e7,0.2)]); let stiff = material(1e8,vec![]);
    let mut u = birth(8,Support::SimplySupported,&soft);
    let mut l = birth(8,Support::SimplySupported,&stiff);
    let pairs = [pair(1,1,2,5e-6)];
    for i in 1..=10 {
        let tu = temperatures(&u,300.0+f64::from(i)); let tl = temperatures(&l,310.0);
        let inputs = [body(1,&u,&tu,&[0.0;8]),body(2,&l,&tl,&[0.0;8])];
        let solved = advance_controlled(&inputs,&pairs,f64::from(i)/10.0,config(),&mut || false).unwrap();
        if i == 1 { assert_eq!(solved.report.reactions[0].force_on_upper_n,0.0); }
        if i == 10 { assert!(solved.report.reactions[0].force_on_upper_n>0.0); }
        let replay = advance_controlled(&inputs,&pairs,f64::from(i)/10.0,config(),&mut || false).unwrap();
        assert_eq!(solved.bodies,replay.bodies);
        assert_eq!(solved.report.reactions,replay.report.reactions);
        u = solved.bodies[&1].clone(); l = solved.bodies[&2].clone();
    }
    assert!(u.cells()[0].section().layers()[&7].material_points()[0].branch_strains()[0] != 0.0);
    let tu = temperatures(&u,310.0); let tl = temperatures(&l,310.0);
    let upper = Body { gravity_m_s2: 0.0, ..body(1,&u,&tu,&[-0.1;8]) };
    let lower = Body { gravity_m_s2: 0.0, ..body(2,&l,&tl,&[0.0;8]) };
    let lifted = advance_controlled(&[upper,lower],&pairs,1.1,config(),&mut || false).unwrap();
    assert_eq!(lifted.report.reactions[0].force_on_upper_n,0.0);
    assert!(lifted.report.reactions[0].gap_m>0.0);
}

// G4: a final body's refusal cannot publish the already staged first body.
#[test]
fn late_failure_cancellation_and_shared_budget_exhaustion_publish_nothing() {
    let law = material(1e8,vec![(2e8,0.5)]);
    let a = birth(8,Support::SimplySupported,&law); let b = a.clone();
    let t = temperatures(&a,350.0);
    let inputs = [body(1,&a,&t,&[0.01;8]),body(2,&b,&t,&[0.0;8])];
    let pairs = [pair(1,1,2,0.0)];
    let before = (a.clone(),b.clone());
    let mut failed = inputs;
    failed[1].config.maximum_deflection_m = 1e-30;
    assert!(matches!(advance_controlled(&failed,&pairs,1.0,config(),&mut || false),Err(Error::Beam { body:2,.. })));
    let small = Config { contact: supports::Config { contact: GramOrthantConfig {
        max_columns: 8, ..config().contact.contact }, ..config().contact }, ..config() };
    assert!(advance_controlled(&inputs,&pairs,1.0,small,&mut || false).is_err());
    let mut polls = 0;
    advance_controlled(&inputs,&pairs,1.0,config(),&mut || { polls+=1; false }).unwrap();
    for limit in [0,polls/2,polls-1] {
        let mut count = 0;
        assert!(advance_controlled(&inputs,&pairs,1.0,config(),&mut || { count+=1; count>limit }).is_err());
    }
    assert_eq!((a,b),before);
}

#[test]
fn stale_times_unknown_pairs_and_missing_temperatures_refuse() {
    let law = material(1e8,vec![]);
    let a = birth(8,Support::SimplySupported,&law); let b = a.clone();
    let t = temperatures(&a,300.0);
    let inputs = [body(1,&a,&t,&[0.0;8]),body(2,&b,&t,&[0.0;8])];
    for pairs in [vec![pair(1,1,1,0.0)],vec![pair(1,1,9,0.0)],
        vec![pair(1,1,2,0.0),pair(1,1,2,0.0)],vec![pair(1,1,2,f64::NAN)]] {
        assert!(advance_controlled(&inputs,&pairs,1.0,config(),&mut || false).is_err());
    }
    assert!(advance_controlled(&[inputs[0],inputs[0]],&[],1.0,config(),&mut || false).is_err());
    let advanced = b.advance_controlled(0.5,&t,9.81,&[0.0;8],beam_config(),&mut || false).unwrap().0;
    assert!(advance_controlled(&[inputs[0],body(2,&advanced,&t,&[0.0;8])],&[],1.0,config(),&mut || false).is_err());
    let missing = BTreeMap::new();
    assert!(advance_controlled(&[inputs[0],body(2,&b,&missing,&[0.0;8])],&[],1.0,config(),&mut || false).is_err());
}
