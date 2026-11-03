//! Original finite-zone adhesion drives the actual native mechanical owners.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::continuation::{
    MaugisContinuationBudget, MaugisContinuationError, MaugisPhase,
};
use fs_contact::normal_patch::maugis::mechanical::MaugisAttachedPotential;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, MaugisError, MaugisWorkBudget, SourceMaugisSphereContact,
};
use fs_contact::normal_patch::{NormalLawSource, NormalLawSourceError};
use fs_couple::unilateral_contact::binding::{
    BoundNormalContact, BoundNormalContactError, BoundNormalStorage,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage};

fn model(lambda: f64) -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let radius: f64 = 0.01;
    let a0 =
        (PI * 1e-7 * radius * radius / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    SourceMaugisSphereContact::new(
        normal,
        radius,
        2.0 * d0 / (PI * lambda),
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}
fn scales(law: &SourceMaugisSphereContact) -> (f64, f64) {
    let r = law.effective_radius_m();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let a0 =
        (PI * w * r * r / (4.0 * law.normal_state().elastic().reduced_modulus_pa() / 3.0)).cbrt();
    (a0 * a0 / r, PI * w * r)
}
fn budget(law: &SourceMaugisSphereContact) -> MaugisContinuationBudget {
    let (d0, f0) = scales(law);
    MaugisContinuationBudget {
        work: MaugisWorkBudget {
            absolute_tolerance_j: 1e-12 * f0 * d0,
            relative_tolerance: 1e-11,
            maximum_evaluations: 2049,
        },
        energy_absolute_tolerance_j: 1e-10 * f0 * d0,
        energy_relative_tolerance: 1e-10,
    }
}
fn potential(law: &SourceMaugisSphereContact) -> MaugisAttachedPotential {
    MaugisAttachedPotential::new(law.clone(), budget(law)).unwrap()
}
fn storage(source: MaugisAttachedPotential) -> BoundNormalStorage {
    // Explicit unequal support masses [1 kg, 3 kg], not inferred tip masses.
    let mut q = vec![0.0; 16];
    q[5] = 1.0;
    q[15] = 1.0 / 3.0;
    BoundNormalStorage::new(
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
        2,
        vec![BoundNormalContact::new(source.into(), vec![1.0, -1.0], 0.0).unwrap()],
    )
    .unwrap()
}
fn system(source: MaugisAttachedPotential) -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    PortHamiltonian::new(4, 0, j, vec![0.0; 16], vec![], Box::new(storage(source))).unwrap()
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= scale * tolerance,
        "{actual:e} != {expected:e}; scaled error {:e}",
        (actual - expected).abs() / scale
    );
}
fn kinetic(x: &[f64]) -> f64 {
    0.5 * x[1] * x[1] + x[3] * x[3] / 6.0
}

#[test]
fn original_branch_and_path_work_drive_reciprocal_native_efforts() {
    for lambda in [0.02, 0.2, 5.0, 1000.0] {
        let law = model(lambda);
        let (d0, f0) = scales(&law);
        let owner = potential(&law);
        let source = NormalLawSource::from(owner.clone());
        assert_eq!(source.identity(), law.identity());
        assert!(source.supports_signed_approach());
        let from = 0.5 * law.minimum_attached_indentation_m();
        let to = 2.0 * d0;
        let expected = law.evaluate_contact_energy(from).unwrap();
        assert!(expected.normal_force_n() < 0.0);
        assert_eq!(
            source.energy_and_force(from).unwrap(),
            (expected.potential_j(), expected.normal_force_n())
        );
        let path = owner.path(from, to).unwrap();
        // Independent original-coordinate integral: this does not use the new
        // symmetric physical-indentation quadrature or subtract its energies.
        let work = law
            .reversible_work_between(from, to, budget(&law).work)
            .unwrap();
        close(path.work.work_j, work.work_j, f0 * d0, 2e-10);
        close(
            owner.path(to, from).unwrap().mean_force_n,
            path.mean_force_n,
            f0,
            1e-14,
        );
        close(
            owner.path(to, from).unwrap().work.work_j,
            -path.work.work_j,
            f0 * d0,
            1e-14,
        );
        let attached = storage(owner);
        let p0 = (0.75 * f0 * d0).sqrt();
        let a = [from, 0.2 * p0, 0.0, -0.3 * p0];
        let b = [to, 0.3 * p0, 0.0, -0.1 * p0];
        let effort = fs_phs::discrete_gradient(&attached, &a, &b);
        assert_eq!(effort[0], -effort[2]);
        close(effort[0], path.mean_force_n, f0, 1e-13);
        close(effort[1], 0.5 * (a[1] + b[1]), p0, 1e-14);
        close(effort[3], (a[3] + b[3]) / 6.0, p0, 1e-14);
        let actual_work: f64 = effort
            .iter()
            .zip(a.iter().zip(b))
            .map(|(e, (a, b))| e * (b - a))
            .sum();
        let energy = law.evaluate_contact_energy(to).unwrap().potential_j()
            - expected.potential_j()
            + kinetic(&b)
            - kinetic(&a);
        close(actual_work, energy, f0 * d0, 2e-10);
        let mut loss = [f64::NAN; 4];
        attached
            .dissipative_flow_into(&a, &effort, 0.001, &mut loss)
            .unwrap();
        assert_eq!(loss, [0.0; 4]);
    }
}

#[test]
fn original_adhesion_accelerates_and_reverses_unequal_bodies_without_losing_momentum() {
    let law = model(5.0);
    let (d0, f0) = scales(&law);
    let system = system(potential(&law));
    let time = (0.75 * d0 / f0).sqrt();
    let p0 = (0.75 * d0 * f0).sqrt();
    let mut x = vec![0.1875 * d0, 0.0, -0.0625 * d0, 0.0];
    let initial = law
        .evaluate_contact_energy(x[0] - x[2])
        .unwrap()
        .potential_j();
    let mut closing = false;
    let mut reversed = false;
    let mut repulsive = false;
    for index in 0..180 {
        let step = fs_phs::step(&system, &x, &[], 0.04 * time).unwrap();
        if index == 0 {
            assert!(step.x[1] > 0.0 && step.x[3] < 0.0);
        }
        x = step.x;
        let response = law.evaluate_contact_energy(x[0] - x[2]).unwrap();
        closing |= x[1] - x[3] / 3.0 > 0.1 * p0;
        reversed |= closing && x[1] - x[3] / 3.0 < -0.1 * p0;
        repulsive |= response.normal_force_n() > 0.0;
        close(x[1] + x[3], 0.0, p0, 2e-9);
        close(x[0] + 3.0 * x[2], 0.0, d0, 2e-8);
        close(response.potential_j() + kinetic(&x), initial, f0 * d0, 2e-8);
        assert_eq!(step.dissipated, 0.0);
        assert_eq!(step.supplied, 0.0);
    }
    assert!(closing && reversed && repulsive);
}

#[test]
fn prepared_motion_matches_reference_and_branch_refusals_preserve_accepted_state() {
    let law = model(0.2);
    let (d0, f0) = scales(&law);
    let owner = potential(&law);
    let system = system(owner.clone());
    let mut prepared = StepWorkspace::new(&system).unwrap();
    let dt = 0.02 * (0.75 * d0 / f0).sqrt();
    let p0 = (0.75 * d0 * f0).sqrt();
    let mut a = vec![0.1875 * d0, 0.0, -0.0625 * d0, 0.0];
    let mut b = a.clone();
    for _ in 0..24 {
        a = fs_phs::step(&system, &a, &[], dt).unwrap().x;
        let mut next = [f64::NAN; 4];
        prepared
            .step_into(&system, &b, &[], dt, &mut next, &mut [])
            .unwrap();
        b = next.to_vec();
        for i in 0..4 {
            close(a[i], b[i], if i % 2 == 0 { d0 } else { p0 }, 2e-8);
        }
    }
    let endpoint = law.minimum_attached_indentation_m();
    let source = NormalLawSource::from(owner.clone());
    assert!(source.energy_and_force(endpoint).is_ok());
    assert!(matches!(
        source.potential(endpoint),
        Err(NormalLawSourceError::Maugis(
            MaugisContinuationError::Contact(MaugisError::EndpointTangent)
        ))
    ));
    let invalid = [endpoint.next_down(), 0.0, 0.0, 0.0];
    let rejected = storage(owner);
    assert!(matches!(
        rejected.try_hamiltonian(&invalid),
        Err(BoundNormalContactError::Source(
            NormalLawSourceError::Maugis(MaugisContinuationError::Contact(
                MaugisError::NoContactEquilibrium { .. }
            ))
        ))
    ));
    assert!(fs_phs::step(&system, &invalid, &[], dt).is_err());
    let saved = b.clone();
    let mut output = [17.0; 4];
    assert!(
        prepared
            .step_into(&system, &invalid, &[], dt, &mut output, &mut [])
            .is_err()
    );
    assert_eq!(output, [17.0; 4]);
    assert_eq!(b, saved);
    let mut polls = 0;
    assert!(matches!(
        prepared.step_into_controlled(&system, &b, &[], dt, &mut output, &mut [], || {
            polls += 1;
            polls >= 3
        }),
        Err(fs_phs::PreparedStepError::Cancelled)
    ));
    assert_eq!(output, [17.0; 4]);
    assert_eq!(b, saved);
    let expected = fs_phs::step(&system, &b, &[], dt).unwrap();
    prepared
        .step_into(&system, &b, &[], dt, &mut output, &mut [])
        .unwrap();
    for i in 0..4 {
        close(
            output[i],
            expected.x[i],
            if i % 2 == 0 { d0 } else { p0 },
            2e-8,
        );
    }
}

#[test]
fn tiny_paths_remain_resolved_and_missing_work_cannot_fall_back_silently() {
    let law = model(5.0);
    let (d0, f0) = scales(&law);
    let owner = potential(&law);
    for exponent in [12, 24, 40, 50] {
        let a = d0;
        let b = a * (1.0 + 2.0_f64.powi(-exponent));
        let path = owner.path(a, b).unwrap();
        let force = law.evaluate_attached(a.midpoint(b)).unwrap().normal_force_n;
        close(path.mean_force_n, force, f0, 2e-8);
        assert!(path.work.work_j.is_finite());
    }
    let mut short = budget(&law);
    short.work.maximum_evaluations = 5;
    short.work.absolute_tolerance_j = 1e-30 * f0 * d0;
    short.work.relative_tolerance = 0.0;
    let owner = MaugisAttachedPotential::new(law.clone(), short).unwrap();
    assert!(matches!(
        owner.path(d0, 10.0 * d0),
        Err(MaugisContinuationError::Contact(
            MaugisError::WorkBudgetExceeded
        ))
    ));
    let rejected = storage(owner);
    let a = [d0, 0.0, 0.0, 0.0];
    let b = [10.0 * d0, 0.0, 0.0, 0.0];
    let mid = [5.5 * d0, 0.0, 0.0, 0.0];
    assert!(
        fs_phs::discrete_gradient(&rejected, &a, &b)
            .iter()
            .all(|v| v.is_nan())
    );
    assert!(
        rejected
            .midpoint_energy_remainder(&a, &b, &mid)
            .unwrap()
            .is_nan()
    );
    struct Opaque;
    impl Storage for Opaque {
        fn hamiltonian(&self, _: &[f64]) -> f64 {
            0.0
        }
        fn gradient(&self, _: &[f64], out: &mut [f64]) {
            out.fill(0.0);
        }
    }
    let site = BoundNormalContact::new(potential(&law).into(), vec![1.0], 0.0).unwrap();
    assert!(BoundNormalStorage::new(Box::new(Opaque), 1, vec![site]).is_err());
}

#[test]
fn continuous_noncontact_work_matches_original_disk_and_contact_integrals() {
    for lambda in [0.02, 0.2, 0.8] {
        let law = model(lambda);
        let (d0, f0) = scales(&law);
        let owner = potential(&law).with_noncontact_continuation().unwrap();
        let limits = owner.noncontact_continuation().unwrap();
        assert_eq!(owner.identity(), law.identity());
        assert!(potential(&law).noncontact_continuation().is_none());
        let lo = limits.minimum_indentation_m;
        let touch = limits.first_contact_indentation_m;
        let from = lo + 0.3 * (touch - lo);
        let disk_to = lo + 0.8 * (touch - lo);
        let disk = owner.path(from, disk_to).unwrap();
        assert_eq!(disk.work.evaluations, 0);
        assert_eq!(disk.work.estimated_absolute_error_j, 0.0);

        // Independent quadrature in PHYSICAL disk radius: original uniform
        // traction F=-pi*sigma*c^2 and delta'(c)=c/R-4*sigma/(pi*E*).
        let r = law.effective_radius_m();
        let e = law.normal_state().elastic().reduced_modulus_pa();
        let sigma = law.normal_state().adhesion_energy_j_per_m2() / law.interaction_range_m();
        let force_work = |a: f64, b: f64| {
            let c0 = law.evaluate_noncontact(a).unwrap().cohesive_radius_m();
            let c1 = law.evaluate_noncontact(b).unwrap().cohesive_radius_m();
            let f = |c: f64| -PI * sigma * c * c * (c / r - 4.0 * sigma / (PI * e));
            (c1 - c0) * (f(c0) + 4.0 * f(c0.midpoint(c1)) + f(c1)) / 6.0
        };
        close(disk.work.work_j, force_work(from, disk_to), f0 * d0, 2e-11);
        let to = 0.5 * d0;
        let path = owner.path(from, to).unwrap();
        let original_contact = law
            .reversible_work_between(touch, to, budget(&law).work)
            .unwrap();
        close(
            path.work.work_j,
            force_work(from, touch) + original_contact.work_j,
            f0 * d0,
            2e-10,
        );
        let reverse = owner.path(to, from).unwrap();
        close(reverse.mean_force_n, path.mean_force_n, f0, 2e-14);
        close(reverse.work.work_j, -path.work.work_j, f0 * d0, 2e-14);
        close(
            owner.path(from, touch).unwrap().work.work_j
                + owner.path(touch, to).unwrap().work.work_j,
            path.work.work_j,
            f0 * d0,
            2e-10,
        );
        let gap = owner.evaluate(from).unwrap();
        assert_eq!(gap.phase(), MaugisPhase::Noncontact);
        assert!(gap.normal_force_n() < 0.0 && gap.central_gap_m() > 0.0);
        assert_eq!(gap.friction_contact_state().unwrap().area_m2(), 0.0);
        assert_eq!(owner.evaluate(touch).unwrap().phase(), MaugisPhase::Contact);
        assert_eq!(
            owner
                .evaluate(touch)
                .unwrap()
                .friction_contact_state()
                .unwrap()
                .area_m2(),
            0.0
        );
        assert!(potential(&law).evaluate(from).is_err());
        // Contacting values remain the exact original branch under either policy.
        assert_eq!(
            owner.evaluate(to).unwrap(),
            potential(&law).evaluate(to).unwrap()
        );
    }
}

#[test]
fn native_inertia_makes_and_loses_solid_contact_without_impulse_or_shear_area_in_gap() {
    let law = model(0.2);
    let (d0, f0) = scales(&law);
    let owner = potential(&law).with_noncontact_continuation().unwrap();
    let limits = owner.noncontact_continuation().unwrap();
    let start = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let system = system(owner.clone());
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let mut x = [0.75 * start, 0.0, -0.25 * start, 0.0];
    let initial = owner.evaluate(x[0] - x[2]).unwrap().potential_j();
    let time = (0.75 * d0 / f0).sqrt();
    let p0 = (0.75 * f0 * d0).sqrt();
    let mut made_contact = false;
    let mut lost_contact = false;
    let mut gap_samples = 0;
    for _ in 0..400 {
        let mut next = [f64::NAN; 4];
        let report = workspace
            .step_into(&system, &x, &[], 0.03 * time, &mut next, &mut [])
            .unwrap();
        x = next;
        let response = owner.evaluate(x[0] - x[2]).unwrap();
        if response.phase() == MaugisPhase::Noncontact {
            gap_samples += 1;
            lost_contact |= made_contact;
            assert!(response.central_gap_m() > 0.0);
            assert!(response.normal_force_n() < 0.0);
            assert_eq!(response.friction_contact_state().unwrap().area_m2(), 0.0);
        } else {
            made_contact = true;
            assert!(response.friction_contact_state().unwrap().area_m2() > 0.0);
        }
        close(x[1] + x[3], 0.0, p0, 2e-9);
        close(x[0] + 3.0 * x[2], 0.0, d0, 2e-8);
        close(response.potential_j() + kinetic(&x), initial, f0 * d0, 3e-8);
        assert_eq!(report.dissipated, 0.0);
        assert_eq!(report.supplied, 0.0);
    }
    assert!(gap_samples > 2 && made_contact && lost_contact);
}

#[test]
fn noncontact_tiny_paths_and_fold_refusal_preserve_native_state() {
    let law = model(0.2);
    let (d0, f0) = scales(&law);
    let owner = potential(&law).with_noncontact_continuation().unwrap();
    let limits = owner.noncontact_continuation().unwrap();
    let fold = limits.minimum_indentation_m;
    let middle = fold.midpoint(limits.first_contact_indentation_m);
    for from in [fold, middle] {
        let to = from.next_up();
        let path = owner.path(from, to).unwrap();
        close(
            path.mean_force_n,
            owner.evaluate(from.midpoint(to)).unwrap().normal_force_n(),
            f0,
            2e-8,
        );
        assert!(path.work.work_j < 0.0);
        assert_eq!(path.work.evaluations, 0);
    }
    assert!(matches!(
        owner.evaluate(fold.next_down()),
        Err(MaugisContinuationError::NoNoncontactEquilibrium { .. })
    ));
    assert!(owner.path(middle, fold.next_down()).is_err());
    assert!(owner.path(f64::NAN, middle).is_err());
    assert!(
        potential(&model(5.0))
            .with_noncontact_continuation()
            .is_err()
    );
    let system = system(owner);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let dt = 0.01 * (0.75 * d0 / f0).sqrt();
    let invalid = [fold.next_down(), 0.0, 0.0, 0.0];
    let mut output = [17.0; 4];
    assert!(
        workspace
            .step_into(&system, &invalid, &[], dt, &mut output, &mut [])
            .is_err()
    );
    assert_eq!(output, [17.0; 4]);
    let accepted = [middle, 0.0, 0.0, 0.0];
    let reference = fs_phs::step(&system, &accepted, &[], dt).unwrap();
    workspace
        .step_into(&system, &accepted, &[], dt, &mut output, &mut [])
        .unwrap();
    for i in 0..4 {
        close(
            output[i],
            reference.x[i],
            if i % 2 == 0 {
                d0
            } else {
                (0.75 * f0 * d0).sqrt()
            },
            2e-8,
        );
    }
}
