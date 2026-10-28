//! Native DMT storage and lossless dynamics, using the real source-card fixture.
#[path = "../../../../fs-contact/tests/source_point/common.rs"]
mod common;
use super::*;
use fs_contact::normal_patch::dmt::{DmtEnvelope, DmtError, SourceDmtSphereContact};
use fs_math::dd::Dd;
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn model() -> SourceDmtSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    SourceDmtSphereContact::new(
        common::bind(&a, &b, None, 0.01).normal_state().clone(),
        1e-6,
        1e-9,
        DmtEnvelope {
            maximum_indentation_m: 1e-9,
            half_space_depth_m: 1e-3,
            layer_thickness_m: 1e-3,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}
fn masses(count: usize) -> QuadraticStorage {
    // Explicit 1 kg effective support masses, not an inferred microscopic mass.
    let n = 2 * count;
    let mut q = vec![0.0; n * n];
    for i in 0..count {
        q[(2 * i + 1) * n + 2 * i + 1] = 1.0;
    }
    QuadraticStorage::new(q, n).unwrap()
}
fn storage(law: &SourceDmtSphereContact, row: Vec<f64>, gap: f64) -> BoundNormalStorage {
    let n = row.len();
    BoundNormalStorage::new(
        Box::new(masses(n)),
        n,
        vec![BoundNormalContact::new(law.clone().into(), row, gap).unwrap()],
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scaled error {:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn signed_gap_and_contact_supply_equal_opposite_force_without_fictitious_damping() {
    let law = model();
    let storage = storage(&law, vec![1.0, -1.0], 0.0);
    for d in [-2e-9, -1e-9, -6e-10, 0.0, 2e-11, 2e-10] {
        let x = [d, 3e-9, 0.0, -2e-9];
        let response = law.evaluate(d).unwrap();
        close(
            storage.try_hamiltonian(&x).unwrap(),
            response.reversible_energy_j + 0.5 * (x[1] * x[1] + x[3] * x[3]),
            1e-15,
            1e-14,
        );
        let mut gradient = [0.0; 4];
        storage.try_gradient(&x, &mut gradient).unwrap();
        assert_eq!(gradient[0], response.normal_force_n);
        assert_eq!(gradient[2], -gradient[0]);
        assert_eq!(gradient[1], x[1]);
        assert_eq!(gradient[3], x[3]);
        let mut loss = [1.0; 4];
        storage
            .dissipative_flow_into(&x, &gradient, 1e-3, &mut loss)
            .unwrap();
        assert_eq!(loss, [0.0; 4]);
    }
    assert!(law.evaluate(-6e-10).unwrap().normal_force_n < 0.0);
    assert!(storage.exact_quadratic_matrix().is_none());
}

#[test]
fn corners_keep_continuous_energy_and_force_but_refuse_a_smooth_tangent() {
    let law = model();
    let source = NormalLawSource::from(law.clone());
    assert_eq!(source.identity(), law.identity());
    assert!(source.supports_signed_approach());
    for d in [-law.interaction_range_m(), 0.0] {
        let (u, f) = source.energy_and_force(d).unwrap();
        let response = law.evaluate(d).unwrap();
        assert_eq!(
            (u, f),
            (response.reversible_energy_j, response.normal_force_n)
        );
        assert!(matches!(
            source.potential(d),
            Err(NormalLawSourceError::Dmt(DmtError::CornerTangent))
        ));
    }
    let storage = storage(&law, vec![1.0], 0.0);
    assert!(matches!(
        storage.try_hamiltonian(&[law.envelope().maximum_indentation_m.next_up(), 0.0]),
        Err(BoundNormalContactError::Source(NormalLawSourceError::Dmt(
            DmtError::OutsideApplicability { .. }
        )))
    ));
    assert!(storage.try_hamiltonian(&[f64::NAN, 0.0]).is_err());
    assert!(storage.hamiltonian(&[f64::NAN, 0.0]).is_nan());
    assert!(storage.try_gradient(&[0.0, 0.0], &mut [0.0]).is_err());
}

fn contact_energy(law: &SourceDmtSphereContact, x: f64) -> Dd {
    let delta = Dd::from_f64(x);
    let coefficient = (4.0 / 3.0)
        * law.normal_state().elastic().reduced_modulus_pa()
        * law.effective_radius_m().sqrt();
    let elastic = Dd::from_f64(2.0) * Dd::from_f64(coefficient) / Dd::from_f64(5.0)
        * delta
        * delta
        * delta.sqrt();
    let f0 = Dd::from_f64(-law.force_controlled_pull_off_n());
    elastic - f0 * (delta + Dd::from_f64(0.5 * law.interaction_range_m()))
}
fn contact_force(law: &SourceDmtSphereContact, x: f64) -> Dd {
    let coefficient = (4.0 / 3.0)
        * law.normal_state().elastic().reduced_modulus_pa()
        * law.effective_radius_m().sqrt();
    let delta = Dd::from_f64(x);
    Dd::from_f64(coefficient) * delta * delta.sqrt()
        + Dd::from_f64(law.force_controlled_pull_off_n())
}

#[test]
fn short_step_work_and_cubic_remainder_agree_with_double_double_reference() {
    let law = model();
    let a = 2e-10;
    for exponent in [8, 14, 20] {
        let b = a * (1.0 + 2.0_f64.powi(-exponent));
        let mid = a.midpoint(b);
        let work = contact_energy(&law, b) - contact_energy(&law, a);
        let remainder = work - contact_force(&law, mid) * (Dd::from_f64(b) - Dd::from_f64(a));
        close(
            law.reversible_work_between(a, b).unwrap(),
            work.to_f64(),
            work.to_f64().abs(),
            2e-14,
        );
        close(
            law.midpoint_energy_remainder(a, b, mid).unwrap(),
            remainder.to_f64(),
            remainder.to_f64().abs(),
            2e-6,
        );
        assert_eq!(
            law.reversible_work_between(b, a).unwrap(),
            -law.reversible_work_between(a, b).unwrap()
        );
    }
    for (a, b) in [(-1.2e-9, -0.8e-9), (-2e-10, 2e-10), (-2e-9, 1e-9)] {
        let expected = law.evaluate(b).unwrap().reversible_energy_j
            - law.evaluate(a).unwrap().reversible_energy_j;
        close(
            law.reversible_work_between(a, b).unwrap(),
            expected,
            1e-15,
            2e-14,
        );
        let mid = a.midpoint(b);
        let rem = expected - law.evaluate(mid).unwrap().normal_force_n * (b - a);
        close(
            law.midpoint_energy_remainder(a, b, mid).unwrap(),
            rem,
            1e-15,
            3e-14,
        );
    }
}

#[test]
fn native_discrete_gradient_preserves_work_with_a_two_body_affine_contact_map() {
    let law = model();
    let storage = storage(&law, vec![2.0, -0.75], 1e-10);
    for (a, b) in [
        ([0.0, 0.0, 2e-10, 0.0], [3e-10, 2e-9, 0.0, -2e-9]),
        ([2e-10, 1e-9, 0.0, -1e-9], [2.000001e-10, 1e-9, 0.0, -1e-9]),
    ] {
        let grad = fs_phs::discrete_gradient(&storage, &a, &b);
        let mut work = Dd::ZERO;
        for i in 0..4 {
            work = work + Dd::from_f64(grad[i]) * (Dd::from_f64(b[i]) - Dd::from_f64(a[i]));
        }
        let da = 2.0 * a[0] - 0.75 * a[2] - 1e-10;
        let db = 2.0 * b[0] - 0.75 * b[2] - 1e-10;
        let kinetic = 0.5 * (b[1] * b[1] + b[3] * b[3] - a[1] * a[1] - a[3] * a[3]);
        let expected = law.reversible_work_between(da, db).unwrap() + kinetic;
        close(work.to_f64(), expected, expected.abs(), 2e-9);
    }
    struct Opaque;
    impl Storage for Opaque {
        fn hamiltonian(&self, _: &[f64]) -> f64 {
            0.0
        }
        fn gradient(&self, _: &[f64], out: &mut [f64]) {
            out.fill(0.0);
        }
    }
    let opaque = BoundNormalStorage::new(
        Box::new(Opaque),
        1,
        vec![BoundNormalContact::new(law.into(), vec![1.0], 0.0).unwrap()],
    )
    .unwrap();
    assert!(
        opaque
            .midpoint_energy_remainder(&[1e-10, 0.0], &[2e-10, 0.0], &[1.5e-10, 0.0])
            .is_none()
    );
}

#[test]
fn native_lossless_trajectory_accelerates_across_the_gap_then_makes_contact() {
    let law = model();
    let gap = 6e-10;
    let storage = storage(&law, vec![1.0, -1.0], gap);
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let system = PortHamiltonian::new(4, 0, j, vec![0.0; 16], vec![], Box::new(storage)).unwrap();
    let mut state = vec![0.0; 4];
    let initial = law.evaluate(-gap).unwrap().reversible_energy_j;
    let mut contact = false;
    for index in 0..600 {
        let step = fs_phs::step(&system, &state, &[], 0.001).unwrap();
        if index == 0 {
            assert!(step.x[0] > 0.0 && step.x[1] > 0.0);
            assert!(step.x[2] < 0.0 && step.x[3] < 0.0);
            assert!(step.x[0] - step.x[2] - gap < 0.0);
        }
        state = step.x;
        let approach = state[0] - state[2] - gap;
        contact |= approach > 0.0;
        close(state[0] + state[2], 0.0, 1e-9, 1e-10);
        close(state[1] + state[3], 0.0, 1e-8, 1e-10);
        let energy = law.evaluate(approach).unwrap().reversible_energy_j
            + 0.5 * (state[1] * state[1] + state[3] * state[3]);
        close(energy, initial, initial.abs(), 2e-8);
        assert_eq!(step.dissipated, 0.0);
        assert_eq!(step.supplied, 0.0);
    }
    assert!(
        contact,
        "the signed attractive tail must drive actual contact"
    );
}

#[test]
fn prepared_and_reference_steps_agree_under_a_coordinate_translation() {
    let law = model();
    let make = || {
        let mut j = vec![0.0; 16];
        j[1] = 1.0;
        j[4] = -1.0;
        j[11] = 1.0;
        j[14] = -1.0;
        PortHamiltonian::new(
            4,
            0,
            j,
            vec![0.0; 16],
            vec![],
            Box::new(storage(&law, vec![1.0, -1.0], 6e-10)),
        )
        .unwrap()
    };
    let system = make();
    let translated = make();
    let mut prepared = fs_phs::StepWorkspace::new(&system).unwrap();
    let mut a = vec![0.0; 4];
    let shift = 1e-9;
    let mut b = vec![shift, 0.0, shift, 0.0];
    let mut c = vec![0.0; 4];
    for _ in 0..240 {
        a = fs_phs::step(&system, &a, &[], 0.001).unwrap().x;
        b = fs_phs::step(&translated, &b, &[], 0.001).unwrap().x;
        let mut next = vec![f64::NAN; 4];
        prepared
            .step_into(&system, &c, &[], 0.001, &mut next, &mut [])
            .unwrap();
        c = next;
        for i in 0..4 {
            close(c[i], a[i], 1e-8, 2e-10);
            close(
                b[i] - if i % 2 == 0 { shift } else { 0.0 },
                a[i],
                1e-8,
                2e-10,
            );
        }
    }
}
