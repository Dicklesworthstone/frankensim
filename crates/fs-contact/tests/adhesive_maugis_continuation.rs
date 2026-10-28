//! Independent work, elastic energy and noncontact geometry, with original sources.
#[path = "source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::*;
use fs_contact::normal_patch::maugis::continuation::*;
use fs_tribo::InputAuthority;

fn model(lambda: f64) -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let r = 0.01;
    let e = normal.elastic().reduced_modulus_pa();
    let a0 = (PI*1e-7*r*r/(4.0*e/3.0)).cbrt();
    let d0 = a0*a0/r;
    SourceMaugisSphereContact::new(normal, r, 2.0*d0/(PI*lambda), MaugisEnvelope {
        maximum_indentation_m: 20.0*d0, half_space_depth_m: 0.1, layer_thickness_m: 0.1,
        maximum_interaction_radius_to_radius: 0.1, maximum_interaction_radius_to_depth: 0.1,
        maximum_interaction_radius_to_layer: 0.1, maximum_nominal_strain: 0.1,
        maximum_repulsive_pressure_to_yield: 0.5, maximum_adhesive_stress_to_yield: 0.9,
        maximum_interaction_range_to_radius: 0.01,
    }).unwrap()
}
fn scales(law: &SourceMaugisSphereContact) -> (f64, f64, f64) {
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let r = law.effective_radius_m();
    let a0 = (PI*w*r*r/(4.0*e/3.0)).cbrt();
    (a0, a0*a0/r, PI*w*r)
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*scale, "{a:e} != {b:e}, scale {scale:e}");
}
fn integral(f: impl Fn(f64)->f64, a: f64, b: f64) -> f64 {
    let n = 4096;
    let h = (b-a)/f64::from(n);
    let mut sum = f(a)+f(b);
    for i in 1..n { sum += (if i%2==0 {2.0} else {4.0}) * f(a+f64::from(i)*h); }
    sum*h/3.0
}

#[test]
fn absolute_contact_energy_has_the_original_force_gradient_and_work() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        for d in [0.0, 0.5, 3.0, 10.0] {
            let h = 1e-5*d0;
            let center = law.evaluate_contact_energy(d*d0).unwrap();
            let left = law.evaluate_contact_energy(d*d0-h).unwrap();
            let right = law.evaluate_contact_energy(d*d0+h).unwrap();
            close((right.potential_j()-left.potential_j())/(2.0*h), center.normal_force_n(), p0, 1e-6);
            assert_eq!(center.contact_response(), Some(&law.evaluate_attached(d*d0).unwrap()));
            assert_eq!(center.authority(), InputAuthority::SyntheticFixture);
            assert_eq!(center.model_identity(), law.identity());
            close(center.elastic_energy_j()+center.surface_energy_j(), center.potential_j(), p0*d0, 1e-11);
        }
        let work = law.reversible_work_between(0.0, 3.0*d0, MaugisWorkBudget {
            absolute_tolerance_j: 1e-10*p0*d0, relative_tolerance: 1e-10, maximum_evaluations: 10000,
        }).unwrap();
        close(law.evaluate_contact_energy(3.0*d0).unwrap().potential_j()
            - law.evaluate_contact_energy(0.0).unwrap().potential_j(), work.work_j, p0*d0, 2e-8);
    }
}

#[test]
fn elastic_storage_matches_independent_abel_line_load_quadrature() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let energy = law.evaluate_contact_energy(2.0*d0).unwrap();
        let response = energy.contact_response().unwrap();
        let a = response.contact_radius_m;
        let c = response.cohesive_radius_m;
        let e = law.normal_state().elastic().reduced_modulus_pa();
        let sigma = law.normal_state().adhesion_energy_j_per_m2()/law.interaction_range_m();
        let r = law.effective_radius_m();
        let contact = integral(|x| e*(response.indentation_m-x*x/r).powi(2), 0.0, a);
        let annulus = integral(|x| (4.0*sigma*sigma/e)*(c-x)*(c+x), a, c);
        close(energy.elastic_energy_j(), contact+annulus, p0*d0, 1e-8);
        assert!(energy.surface_energy_j() < 0.0);
    }
}

#[test]
fn noncontact_disk_obeys_deformed_gap_force_tangent_and_energy() {
    for lambda in [0.02, 0.2, 0.8] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let limits = law.noncontact_limits().unwrap().unwrap();
        let d = 0.25*limits.minimum_indentation_m+0.75*limits.first_contact_indentation_m;
        let result = law.evaluate_noncontact(d).unwrap();
        let c = result.cohesive_radius_m();
        let e = law.normal_state().elastic().reduced_modulus_pa();
        let sigma = law.normal_state().adhesion_energy_j_per_m2()/law.interaction_range_m();
        let r = law.effective_radius_m();
        close(d, c*c/(2.0*r)-4.0*sigma*c/(PI*e)-law.interaction_range_m(), d0, 1e-12);
        close(result.central_gap_m(), -d-2.0*sigma*c/e, d0, 1e-12);
        close(result.normal_force_n(), -PI*sigma*c*c, p0, 1e-12);
        close(result.elastic_energy_j(), 8.0*sigma*sigma*c*c*c/(3.0*e), p0*d0, 1e-12);
        let h = 1e-5*(limits.first_contact_indentation_m-limits.minimum_indentation_m);
        let lo = law.evaluate_noncontact(d-h).unwrap();
        let hi = law.evaluate_noncontact(d+h).unwrap();
        close((hi.potential_j()-lo.potential_j())/(2.0*h), result.normal_force_n(), p0, 1e-7);
        close((hi.normal_force_n()-lo.normal_force_n())/(2.0*h), result.tangent_n_per_m().unwrap(), p0/d0, 1e-7);
        assert!(result.normal_force_n() < 0.0 && result.central_gap_m() > 0.0);
        assert_eq!(result.friction_contact_state().unwrap().area_m2(), 0.0);
        assert!(result.contact_response().is_none());
        let port = result.point_port([0.0, 0.0, 1.0], [0.01, 0.0, 0.0]).unwrap();
        assert_eq!(port.action_force_n[2], result.normal_force_n());
        assert_eq!(port.reaction_force_n[2], -result.normal_force_n());
        assert_eq!(port.irreversible_work_j, 0.0);
    }
}

#[test]
fn zero_area_junction_has_no_energy_or_force_jump_and_no_shear() {
    for lambda in [0.02, 0.2, 0.8] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let limits = law.noncontact_limits().unwrap().unwrap();
        let disk = law.evaluate_noncontact(limits.first_contact_indentation_m).unwrap();
        let contact = law.evaluate_contact_energy(limits.first_contact_indentation_m).unwrap();
        close(contact.potential_j(), disk.potential_j(), p0*d0, 1e-11);
        close(contact.elastic_energy_j(), disk.elastic_energy_j(), p0*d0, 1e-11);
        close(contact.normal_force_n(), disk.normal_force_n(), p0, 1e-12);
        assert_eq!(disk.central_gap_m(), 0.0);
        assert_eq!(contact.friction_contact_state().unwrap().area_m2(), 0.0);
        assert!(disk.point_port([0.0, 0.0, 1.0], [0.0; 3]).is_err());
        let fold = law.evaluate_noncontact(limits.minimum_indentation_m).unwrap();
        assert!(fold.tangent_n_per_m().is_none());
        assert!(fold.potential_j() > 0.0);
    }
}

#[test]
fn branches_refuse_wrong_domains_without_inventing_contact_or_adhesion() {
    let law = model(0.2);
    let (_, d0, _) = scales(&law);
    let limits = law.noncontact_limits().unwrap().unwrap();
    for d in [f64::NAN, limits.minimum_indentation_m-d0, limits.first_contact_indentation_m+d0] {
        assert!(law.evaluate_noncontact(d).is_err());
    }
    assert!(model(5.0).noncontact_limits().unwrap().is_none());
    assert!(model(5.0).evaluate_noncontact(-d0).is_err());
    assert!(law.evaluate_detached(0.0).is_err());
    let detached = law.evaluate_detached(-2.0*law.interaction_range_m()).unwrap();
    assert_eq!(detached.potential_j(), 0.0);
    assert_eq!(detached.normal_force_n(), 0.0);
    assert_eq!(detached.phase(), MaugisPhase::Detached);
    assert_eq!(detached.friction_contact_state().unwrap().area_m2(), 0.0);
}

fn continuation_budget(law: &SourceMaugisSphereContact) -> MaugisContinuationBudget {
    let (_, d0, p0) = scales(law);
    MaugisContinuationBudget {
        work: MaugisWorkBudget { absolute_tolerance_j: 1e-11*p0*d0,
            relative_tolerance: 1e-10, maximum_evaluations: 10000 },
        energy_absolute_tolerance_j: 2e-9*p0*d0, energy_relative_tolerance: 1e-10,
    }
}

#[test]
fn history_localizes_both_snaps_and_cycle_work_matches_original_force_integration() {
    for lambda in [0.02, 0.2, 0.6, 0.8, 1.0, 5.0, 1000.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let controls = continuation_budget(&law);
        let owner = MaugisDisplacementContinuation::new(law.clone()).unwrap();
        let far = owner.separation().indentation_m()-law.interaction_range_m()-d0;
        let initial = owner.initialize(MaugisPhase::Detached, far).unwrap();
        let up = owner.advance(&initial, 3.0*d0, controls).unwrap();
        assert_eq!(up.snap.unwrap().kind(), MaugisSnapKind::AttractionOnset);
        assert_eq!(up.snap.unwrap().indentation_m(), -law.interaction_range_m());
        assert_eq!(up.next.equilibrium().phase(), MaugisPhase::Contact);
        assert!(up.next.equilibrium().friction_contact_state().unwrap().area_m2() > 0.0);
        let down = owner.advance(&up.next, far, controls).unwrap();
        assert_eq!(down.snap.unwrap().kind(), MaugisSnapKind::Separation);
        assert_eq!(down.next, initial);
        let release = up.released_energy_j+down.released_energy_j;
        assert!(up.released_energy_j > 0.0 && down.released_energy_j > 0.0);
        close(up.work.work_j+down.work.work_j, release, p0*d0, 4e-9);
        // Independent original-force quadrature. Squared displacement resolves
        // the square-root branch near a fold, without using either energy API.
        let on = owner.onset().indentation_m();
        let off = owner.separation().indentation_m();
        let mut independent = 0.0;
        let mut contact_from = off;
        if let Some(limits) = law.noncontact_limits().unwrap() {
            let end = on.min(limits.first_contact_indentation_m);
            let width = end-off;
            independent += integral(|s| {
                let d = if s == 0.0 { off } else if s == 1.0 { end } else { off+width*s*s };
                -law.evaluate_noncontact(d).unwrap().normal_force_n()*2.0*width*s
            }, 0.0, 1.0);
            contact_from = end;
        }
        if contact_from < on {
            let width = on-contact_from;
            independent += integral(|s| {
                let d = if s == 0.0 { contact_from } else if s == 1.0 { on } else { contact_from+width*s*s };
                -law.evaluate_attached(d).unwrap().normal_force_n*2.0*width*s
            }, 0.0, 1.0);
        }
        close(release, independent, release, 2e-5);
        let again = owner.advance(&down.next, far-d0, controls).unwrap();
        assert!(again.snap.is_none());
        assert_eq!(again.released_energy_j, 0.0);
        assert_eq!(again.work.work_j, 0.0);
    }
}

#[test]
fn low_tabor_disk_snaps_match_analytic_release_and_touch_is_continuous() {
    for lambda in [0.02, 0.2, 0.6] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let controls = continuation_budget(&law);
        let owner = MaugisDisplacementContinuation::new(law.clone()).unwrap();
        let fold = 8.0*lambda/(3.0*PI);
        let factor = PI*lambda*fold.powi(4)*p0*d0;
        close(owner.onset().released_energy_j(), (2.0/3.0)*factor, factor, 1e-9);
        close(owner.separation().released_energy_j(), factor/24.0, factor, 1e-10);
        let initial = owner.initialize(MaugisPhase::Detached, -2.0*law.interaction_range_m()).unwrap();
        let on = owner.advance(&initial, owner.onset().indentation_m(), controls).unwrap();
        assert_eq!(on.next.equilibrium().phase(), MaugisPhase::Noncontact);
        assert_eq!(on.work.work_j, 0.0);
        let stationary = owner.advance(&on.next, on.next.indentation_m(), controls).unwrap();
        assert!(stationary.snap.is_none());
        assert_eq!(stationary.released_energy_j, 0.0);
        let touch = law.noncontact_limits().unwrap().unwrap().first_contact_indentation_m;
        let at_touch = owner.advance(&on.next, touch, controls).unwrap();
        assert!(at_touch.snap.is_none());
        assert_eq!(at_touch.next.equilibrium().phase(), MaugisPhase::Contact);
        assert_eq!(at_touch.next.equilibrium().friction_contact_state().unwrap().area_m2(), 0.0);
        let above = owner.advance(&at_touch.next, 0.0, controls).unwrap();
        let below = owner.advance(&above.next, on.next.indentation_m(), controls).unwrap();
        assert!(above.snap.is_none() && below.snap.is_none());
        close(at_touch.work.work_j+above.work.work_j+below.work.work_j, 0.0, p0*d0, 4e-9);
    }
}

#[test]
fn continuation_preserves_metastability_reversal_and_segment_partitioning() {
    for lambda in [0.2, 0.8, 1.0, 5.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let controls = continuation_budget(&law);
        let owner = MaugisDisplacementContinuation::new(law.clone()).unwrap();
        let far = owner.separation().indentation_m()-law.interaction_range_m()-d0;
        let midpoint = 0.5*owner.onset().indentation_m()+0.5*owner.separation().indentation_m();
        let initial = owner.initialize(MaugisPhase::Detached, far).unwrap();
        let unloaded = owner.advance(&initial, midpoint, controls).unwrap();
        assert_eq!(unloaded.next.equilibrium().phase(), MaugisPhase::Detached);
        assert!(unloaded.snap.is_none());
        let all = owner.advance(&initial, 2.0*d0, controls).unwrap();
        let loaded = owner.advance(&all.next, midpoint, controls).unwrap();
        assert_ne!(loaded.next.equilibrium().phase(), MaugisPhase::Detached);
        assert!(loaded.snap.is_none());
        let fraction = 0.1*(owner.onset().indentation_m()-owner.separation().indentation_m());
        let reverse = owner.advance(&loaded.next, midpoint+fraction, controls).unwrap();
        assert!(reverse.snap.is_none());
        let at_on = owner.advance(&initial, owner.onset().indentation_m(), controls).unwrap();
        let rest = owner.advance(&at_on.next, 2.0*d0, controls).unwrap();
        assert!(rest.snap.is_none());
        assert_eq!(rest.next, all.next);
        close(at_on.work.work_j+rest.work.work_j, all.work.work_j, p0*d0, 1e-9);
        assert_eq!(at_on.released_energy_j, all.released_energy_j);
    }
}

#[test]
fn history_refusal_does_not_advance_and_retry_is_identical() {
    let law = model(1.0);
    let (_, d0, _) = scales(&law);
    let controls = continuation_budget(&law);
    let owner = MaugisDisplacementContinuation::new(law.clone()).unwrap();
    let old = owner.initialize(MaugisPhase::Contact, 0.0).unwrap();
    let saved = old;
    assert!(owner.advance(&old, f64::NAN, controls).is_err());
    assert!(owner.advance(&old, 2.0*law.envelope().maximum_indentation_m, controls).is_err());
    let mut small = controls;
    small.work.maximum_evaluations = 1;
    assert!(owner.advance(&old, d0, small).is_err());
    let other = MaugisDisplacementContinuation::new(model(2.0)).unwrap();
    assert!(matches!(other.advance(&old, d0, controls), Err(MaugisContinuationError::HistoryMismatch)));
    assert_eq!(old, saved);
    assert_eq!(owner.advance(&old, d0, controls).unwrap(), owner.advance(&saved, d0, controls).unwrap());
}
