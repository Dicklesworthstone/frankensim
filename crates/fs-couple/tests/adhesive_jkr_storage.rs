//! Original high-Tabor adhesion drives native body motion and preserves work.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, JkrError, SourceJkrSphereContact};
use fs_contact::normal_patch::{NormalLawSource, NormalLawSourceError};
use fs_couple::unilateral_contact::binding::{
    BoundNormalContact, BoundNormalContactError, BoundNormalStorage,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage};

fn model() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.05).normal_state().clone();
    SourceJkrSphereContact::new(
        normal,
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn scales(law: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let force = PI * law.normal_state().adhesion_energy_j_per_m2() * law.effective_radius_m() / 6.0;
    let displacement = law.indentation_scale_m();
    (displacement, force, (0.75 * displacement * force).sqrt())
}
fn storage(source: SourceJkrSphereContact) -> BoundNormalStorage {
    // Caller-owned unequal masses 1 kg and 3 kg, with physical relative motion.
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
fn system(source: SourceJkrSphereContact) -> PortHamiltonian {
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
fn original_jkr_source_and_path_project_only_onto_physical_contact_coordinates() {
    let law = model();
    let (d0, f0, p0) = scales(&law);
    let source = NormalLawSource::from(law.clone());
    assert_eq!(source.identity(), law.identity());
    assert!(source.supports_signed_approach());
    let from = 0.5 * law.minimum_attached_indentation_m();
    let to = 5.0 * d0;
    let response = law.evaluate_attached(from).unwrap();
    assert!(response.normal_force_n < 0.0 && response.contact_area_m2 > 0.0);
    assert_eq!(
        source.energy_and_force(from).unwrap(),
        (response.reversible_energy_j, response.normal_force_n)
    );
    let attached = storage(law.clone());
    let a = [from, 0.2 * p0, 0.0, -0.3 * p0];
    let b = [to, 0.3 * p0, 0.0, -0.1 * p0];
    let effort = fs_phs::discrete_gradient(&attached, &a, &b);
    assert_eq!(effort[0], -effort[2]);
    close(
        effort[0],
        law.attached_path_average(from, to).unwrap(),
        f0,
        2e-14,
    );
    close(effort[1], 0.5 * (a[1] + b[1]), p0, 2e-14);
    close(effort[3], (a[3] + b[3]) / 6.0, p0, 2e-14);
    let work: f64 = effort
        .iter()
        .zip(a.iter().zip(b))
        .map(|(e, (a, b))| e * (b - a))
        .sum();
    let energy = law.evaluate_attached(to).unwrap().reversible_energy_j
        - response.reversible_energy_j
        + kinetic(&b)
        - kinetic(&a);
    close(work, energy, f0 * d0, 2e-12);
    let mut loss = [f64::NAN; 4];
    attached
        .dissipative_flow_into(&a, &effort, 0.001, &mut loss)
        .unwrap();
    assert_eq!(loss, [0.0; 4]);

    // Representable microscopic motion retains a finite effort and changes
    // no unrelated momentum. A global Gonzalez correction need not do this.
    let a = [d0, 0.2 * p0, 0.0, -0.3 * p0];
    let b = [d0.next_up(), a[1], 0.0, a[3]];
    let effort = fs_phs::discrete_gradient(&attached, &a, &b);
    close(
        effort[0],
        law.evaluate_attached(d0).unwrap().normal_force_n,
        f0,
        2e-13,
    );
    assert_eq!(effort[0], -effort[2]);
    assert_eq!(effort[1], a[1]);
    close(effort[3], a[3] / 3.0, p0, 1e-15);
}

#[test]
fn tensile_contact_accelerates_and_reverses_unequal_bodies_with_conserved_energy() {
    let law = model();
    let (d0, f0, p0) = scales(&law);
    let system = system(law.clone());
    let from = -1.5 * d0;
    let mut x = vec![0.75 * from, 0.0, -0.25 * from, 0.0];
    let energy = law.evaluate_attached(from).unwrap().reversible_energy_j;
    let time = (0.75 * d0 / f0).sqrt();
    let mut closing = false;
    let mut reversed = false;
    let mut compression = false;
    for index in 0..180 {
        let step = fs_phs::step(&system, &x, &[], 0.04 * time).unwrap();
        x = step.x;
        if index == 0 {
            assert!(x[1] > 0.0 && x[3] < 0.0);
        }
        let response = law.evaluate_attached(x[0] - x[2]).unwrap();
        closing |= x[1] - x[3] / 3.0 > 0.1 * p0;
        reversed |= closing && x[1] - x[3] / 3.0 < -0.1 * p0;
        compression |= response.normal_force_n > 0.0;
        close(x[1] + x[3], 0.0, p0, 2e-9);
        close(x[0] + 3.0 * x[2], 0.0, d0, 2e-8);
        close(
            response.reversible_energy_j + kinetic(&x),
            energy,
            f0 * d0,
            2e-8,
        );
        assert_eq!(step.dissipated, 0.0);
        assert_eq!(step.supplied, 0.0);
    }
    assert!(closing && reversed && compression);
}

#[test]
fn prepared_motion_matches_native_steps_and_branch_refusal_is_atomic() {
    let law = model();
    let (d0, f0, p0) = scales(&law);
    let system = system(law.clone());
    let mut prepared = StepWorkspace::new(&system).unwrap();
    let dt = 0.02 * (0.75 * d0 / f0).sqrt();
    let mut reference = vec![0.0; 4];
    let mut accepted = reference.clone();
    for _ in 0..24 {
        reference = fs_phs::step(&system, &reference, &[], dt).unwrap().x;
        let mut next = [f64::NAN; 4];
        prepared
            .step_into(&system, &accepted, &[], dt, &mut next, &mut [])
            .unwrap();
        accepted = next.to_vec();
        for i in 0..4 {
            close(
                accepted[i],
                reference[i],
                if i % 2 == 0 { d0 } else { p0 },
                2e-8,
            );
        }
    }
    let endpoint = law.minimum_attached_indentation_m();
    let source = NormalLawSource::from(law.clone());
    assert!(source.energy_and_force(endpoint).is_ok());
    assert!(matches!(
        source.potential(endpoint),
        Err(NormalLawSourceError::Jkr(JkrError::SingularTangent))
    ));
    let invalid = [endpoint.next_down(), 0.0, 0.0, 0.0];
    assert!(matches!(
        storage(law.clone()).try_hamiltonian(&invalid),
        Err(BoundNormalContactError::Source(NormalLawSourceError::Jkr(
            JkrError::NoAttachedEquilibrium { .. }
        )))
    ));
    assert!(fs_phs::step(&system, &invalid, &[], dt).is_err());
    let saved = accepted.clone();
    let mut output = [17.0; 4];
    assert!(
        prepared
            .step_into(&system, &invalid, &[], dt, &mut output, &mut [])
            .is_err()
    );
    assert_eq!(output, [17.0; 4]);
    let mut polls = 0;
    assert!(matches!(
        prepared.step_into_controlled(&system, &accepted, &[], dt, &mut output, &mut [], || {
            polls += 1;
            polls >= 3
        }),
        Err(fs_phs::PreparedStepError::Cancelled)
    ));
    assert_eq!(output, [17.0; 4]);
    assert_eq!(accepted, saved);
    let expected = fs_phs::step(&system, &accepted, &[], dt).unwrap();
    prepared
        .step_into(&system, &accepted, &[], dt, &mut output, &mut [])
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
fn incomplete_inner_energy_cannot_silently_use_a_global_contact_correction() {
    struct OnlyGradient;
    impl Storage for OnlyGradient {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * x.iter().map(|v| v * v).sum::<f64>()
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out.copy_from_slice(x);
        }
    }
    let contact = BoundNormalContact::new(model().into(), vec![1.0, -1.0], 0.0).unwrap();
    assert!(BoundNormalStorage::new(Box::new(OnlyGradient), 2, vec![contact]).is_err());
}

#[test]
fn near_fold_attached_motion_uses_the_original_path_jacobian() {
    let law = model();
    let (d0, f0, p0) = scales(&law);
    let system = system(law.clone());
    // Both the initial state and its inward motion are physically attached.
    // A central 1e-6*d0 numerical probe would cross the fold at this state.
    let approach = law.minimum_attached_indentation_m() + 1e-8 * d0;
    let initial = [0.75 * approach, 0.0, -0.25 * approach, 0.0];
    let energy = law
        .evaluate_attached(initial[0] - initial[2])
        .unwrap()
        .reversible_energy_j;
    let dt = 1e-3 * (0.75 * d0 / f0).sqrt();
    let native = fs_phs::step(&system, &initial, &[], dt).unwrap();
    assert!(native.x[0] - native.x[2] > initial[0] - initial[2]);
    assert!(native.x[1] > 0.0 && native.x[3] < 0.0);
    close(native.x[1] + native.x[3], 0.0, p0, 2e-9);
    close(native.x[0] + 3.0 * native.x[2], 0.0, d0, 2e-8);
    close(
        law.evaluate_attached(native.x[0] - native.x[2])
            .unwrap()
            .reversible_energy_j
            + kinetic(&native.x),
        energy,
        f0 * d0,
        2e-8,
    );
    let mut prepared = StepWorkspace::new(&system).unwrap();
    let mut output = [f64::NAN; 4];
    prepared
        .step_into(&system, &initial, &[], dt, &mut output, &mut [])
        .unwrap();
    for i in 0..4 {
        close(
            output[i],
            native.x[i],
            if i % 2 == 0 { d0 } else { p0 },
            2e-8,
        );
    }
}

#[test]
fn complete_jkr_jacobian_matches_directional_work_and_refuses_the_singular_fold() {
    let law = model();
    let (d0, f0, p0) = scales(&law);
    let attached = storage(law.clone());
    let rule = attached.discrete_gradient_rule().unwrap();
    let a = [-d0, 0.3 * p0, 0.0, -0.2 * p0];
    let b = [3.0 * d0, 0.1 * p0, -d0, 0.2 * p0];
    let direction = [0.8 * d0, 0.3 * p0, -0.2 * d0, -0.4 * p0];
    let mut action = [f64::NAN; 4];
    assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
    assert_eq!(action[0], -action[2]);
    let epsilon = 1e-5;
    let plus = core::array::from_fn::<_, 4, _>(|i| b[i] + epsilon * direction[i]);
    let minus = core::array::from_fn::<_, 4, _>(|i| b[i] - epsilon * direction[i]);
    let positive = fs_phs::discrete_gradient(&attached, &a, &plus);
    let negative = fs_phs::discrete_gradient(&attached, &a, &minus);
    for i in 0..4 {
        close(
            action[i],
            (positive[i] - negative[i]) / (2.0 * epsilon),
            if i % 2 == 0 { f0 } else { p0 },
            3e-8,
        );
    }
    let fold = [law.minimum_attached_indentation_m(), 0.0, 0.0, 0.0];
    assert!(rule.end_jacobian_action(&fold, &fold, &direction, &mut action));
    assert!(action.iter().all(|v| v.is_nan()));
    assert!(rule.end_jacobian_action(&a, &b, &direction[..3], &mut action));
    assert!(action.iter().all(|v| v.is_nan()));
}

#[test]
fn original_inner_rule_cannot_be_replaced_by_its_quadratic_view_for_newton() {
    struct NativeWithoutDerivative;
    impl Storage for NativeWithoutDerivative {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * x.iter().map(|v| v * v).sum::<f64>()
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out.copy_from_slice(x);
        }
        fn discrete_gradient_rule(&self) -> Option<&dyn fs_phs::ConstitutiveDiscreteGradient> {
            Some(self)
        }
        fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
            Some(&[
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ])
        }
    }
    impl fs_phs::ConstitutiveDiscreteGradient for NativeWithoutDerivative {
        fn evaluate(&self, _: &[f64], _: &[f64], midpoint: &[f64], out: &mut [f64]) {
            out.copy_from_slice(midpoint);
        }
    }
    let contact = BoundNormalContact::new(model().into(), vec![1.0, -1.0], 0.0).unwrap();
    let bound =
        BoundNormalStorage::new(Box::new(NativeWithoutDerivative), 2, vec![contact]).unwrap();
    let rule = bound.discrete_gradient_rule().unwrap();
    assert!(!rule.end_jacobian_action(&[0.0; 4], &[0.0; 4], &[1.0; 4], &mut [0.0; 4]));
}
