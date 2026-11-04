//! G1/G3: mixed contact work supplies one complete native Newton action.
#[path = "../../fs-contact/tests/normal_work_derivative/fixtures.rs"]
mod fixtures;
use fixtures::{close, curve, dmt, jkr, point};
use fs_contact::normal_patch::adhesive::SourceJkrSphereContact;
use fs_couple::unilateral_contact::binding::{BoundNormalContact, BoundNormalStorage};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage};

fn scales(law: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let d = law.indentation_scale_m();
    let f = core::f64::consts::PI
        * law.effective_radius_m()
        * law.normal_state().adhesion_energy_j_per_m2()
        / 6.0;
    (d, f, (0.75 * d * f).sqrt())
}

fn storage(law: SourceJkrSphereContact) -> BoundNormalStorage {
    let (d, f, _) = scales(&law);
    let mut q = vec![0.0; 36];
    for (i, mass) in [1.0, 3.0, 2.0].iter().enumerate() {
        q[(2 * i + 1) * 6 + 2 * i + 1] = 1.0 / mass;
    }
    BoundNormalStorage::new(
        Box::new(QuadraticStorage::new(q, 6).unwrap()),
        3,
        vec![
            BoundNormalContact::new(law.into(), vec![1.0, -1.0, 0.0], 0.0).unwrap(),
            BoundNormalContact::new(point(None).into(), vec![0.0, 1.0, -1.0], -0.4 * d).unwrap(),
            BoundNormalContact::new(curve(4.0 * d, f).into(), vec![0.0, 1.0, -1.0], -0.3 * d)
                .unwrap(),
            BoundNormalContact::new(dmt().into(), vec![1.0, 0.0, -1.0], 0.0).unwrap(),
        ],
    )
    .unwrap()
}

fn system(law: SourceJkrSphereContact) -> PortHamiltonian {
    let mut j = vec![0.0; 36];
    for i in 0..3 {
        j[2 * i * 6 + 2 * i + 1] = 1.0;
        j[(2 * i + 1) * 6 + 2 * i] = -1.0;
    }
    PortHamiltonian::new(6, 0, j, vec![0.0; 36], vec![], Box::new(storage(law))).unwrap()
}

#[test]
fn mixed_source_action_matches_complete_directional_work_and_reciprocal_rows() {
    let law = jkr();
    let (d, f, p) = scales(&law);
    let bound = storage(law);
    let a = [0.0, 0.2 * p, -0.2 * d, -0.3 * p, -0.4 * d, 0.1 * p];
    let b = [1.2 * d, 0.4 * p, 0.25 * d, -0.1 * p, -0.3 * d, 0.2 * p];
    let direction = [0.3 * d, 0.2 * p, -0.2 * d, 0.4 * p, 0.1 * d, -0.3 * p];
    let mut action = [f64::NAN; 6];
    assert!(bound.discrete_gradient_rule().unwrap().end_jacobian_action(
        &a,
        &b,
        &direction,
        &mut action
    ));
    close(action[0] + action[2] + action[4], 0.0, f, 3e-14);
    for (i, mass) in [1.0, 3.0, 2.0].iter().enumerate() {
        close(
            action[2 * i + 1],
            0.5 * direction[2 * i + 1] / mass,
            p,
            3e-15,
        );
    }
    let h = 1e-5;
    let plus = core::array::from_fn::<_, 6, _>(|i| b[i] + h * direction[i]);
    let minus = core::array::from_fn::<_, 6, _>(|i| b[i] - h * direction[i]);
    let ep = fs_phs::discrete_gradient(&bound, &a, &plus);
    let em = fs_phs::discrete_gradient(&bound, &a, &minus);
    for i in 0..6 {
        close(
            action[i],
            (ep[i] - em[i]) / (2.0 * h),
            if i % 2 == 0 { f } else { p },
            2e-8,
        );
    }
}

#[test]
fn mixed_native_and_prepared_motion_can_start_just_inside_the_jkr_fold() {
    let law = jkr();
    let (d, f, p) = scales(&law);
    let approach = law.minimum_attached_indentation_m() + 1e-8 * d;
    let system = system(law);
    let initial = [approach, 0.0, 0.0, 0.0, 0.0, 0.0];
    let dt = 1e-3 * (0.75 * d / f).sqrt();
    let energy = system.hamiltonian(&initial);
    // A global central probe crosses the JKR fold here even though the
    // original forces drive an admitted inward step of the assembled bodies.
    let native = fs_phs::step(&system, &initial, &[], dt).unwrap();
    assert!(native.x[0] - native.x[2] > approach);
    close(native.x[1] + native.x[3] + native.x[5], 0.0, p, 2e-9);
    close(
        native.x[0] + 3.0 * native.x[2] + 2.0 * native.x[4],
        approach,
        d,
        2e-8,
    );
    close(system.hamiltonian(&native.x), energy, f * d, 3e-8);
    assert_eq!(native.dissipated, 0.0);
    assert_eq!(native.supplied, 0.0);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let mut output = [f64::NAN; 6];
    workspace
        .step_into(&system, &initial, &[], dt, &mut output, &mut [])
        .unwrap();
    for i in 0..6 {
        close(output[i], native.x[i], if i % 2 == 0 { d } else { p }, 2e-8);
    }
}
