//! Relaxed initial conditions against independent static elasticity and energy.
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl, TensorStressControlBudget,
};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch, TensorPronyError};
use fs_material::{Tangent6, Voigt};

const LAMBDA: f64 = 2.0e9;
const MU: f64 = 3.0e9;
const WEIGHTS: Voigt = [
    1.0, 1.0, 1.0,
    core::f64::consts::SQRT_2, core::f64::consts::SQRT_2, core::f64::consts::SQRT_2,
];

fn law(kelvin: bool) -> TensorProny {
    let elastic: Tangent6 = core::array::from_fn(|i| core::array::from_fn(|j| {
        if i < 3 && j < 3 { LAMBDA + if i == j { 2.0 * MU } else { 0.0 } }
        else if i == j { 2.0 * MU } else { 0.0 }
    }));
    let branch = TensorPronyBranch::from_gram_rows(
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 3.0e4 } else { 0.0 })),
        0.3,
    ).unwrap();
    let original = TensorProny::new_mandel(elastic, vec![branch], 0.03).unwrap();
    if kelvin {
        original.with_kelvin_viscosity_mandel(
            core::array::from_fn(|i| core::array::from_fn(|j| {
                if i == j { (i + 1) as f64 * 1.0e7 } else { 0.0 }
            })),
        ).unwrap()
    } else { original }
}

fn axes(plane: MaterialPlane) -> ([usize; 3], [usize; 3]) {
    match plane {
        MaterialPlane::Axes12 => ([0, 1, 3], [2, 4, 5]),
        MaterialPlane::Axes21 => ([1, 0, 3], [2, 4, 5]),
        MaterialPlane::Axes13 => ([0, 2, 5], [1, 3, 4]),
        MaterialPlane::Axes23 => ([1, 2, 4], [0, 3, 5]),
    }
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2.0e-12 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}");
}

#[test]
fn relaxed_plane_state_matches_static_elasticity_for_every_plane_and_viscosity() {
    let input = [0.001, -0.0004, 0.0006];
    for kelvin in [false, true] {
        for plane in [MaterialPlane::Axes12, MaterialPlane::Axes21,
            MaterialPlane::Axes13, MaterialPlane::Axes23]
        {
            let original = law(kelvin);
            let virgin = original.history();
            let before = virgin.clone();
            let mut initial = virgin.initial_relaxed_plane_stress(plane, input, || false).unwrap();
            let (inside, outside) = axes(plane);
            let mut expected_strain = [0.0; 6];
            expected_strain[inside[0]] = input[0];
            expected_strain[inside[1]] = input[1];
            expected_strain[inside[2]] = 0.5 * input[2];
            expected_strain[outside[0]] = -LAMBDA / (LAMBDA + 2.0 * MU) * (input[0] + input[1]);
            let trace: f64 = expected_strain[..3].iter().sum();
            let expected_stress: Voigt = core::array::from_fn(|i| {
                2.0 * MU * expected_strain[i] + if i < 3 { LAMBDA * trace } else { 0.0 }
            });
            let stored: f64 = (0..6).map(|i| {
                0.5 * WEIGHTS[i] * WEIGHTS[i] * expected_stress[i] * expected_strain[i]
            }).sum();
            for i in 0..6 {
                near(initial.strain()[i], expected_strain[i], 0.001);
                near(initial.stress()[i], expected_stress[i], 1.0e7);
            }
            near(initial.stored_energy(), stored, stored);
            assert_eq!(initial.branch_stored_energy(), 0.0);
            assert_eq!(initial.work(), 0.0);
            assert_eq!(initial.dissipated(), 0.0);
            assert_eq!(virgin, before);
            // A held relaxed state supplies no incremental work or artificial
            // Maxwell/Kelvin heating even though its initial storage is nonzero.
            let held = initial.advance(initial.strain(), 0.7, || false).unwrap();
            near(held.work, 0.0, stored);
            near(held.dissipated, 0.0, stored);
            near(held.stored_energy, stored, stored);
            near(held.balance_residual, 0.0, stored);
        }
    }
}

#[test]
fn source_bound_initial_state_continues_original_mixed_memories() {
    let original = law(false);
    let virgin = original.history();
    let mask = [false, false, true, false, true, true];
    let control = PreparedTensorMixedControl::new(
        &virgin, mask, TensorStressControlBudget::default(), || false,
    ).unwrap();
    let mut initial = virgin.initial_relaxed_plane_stress(
        MaterialPlane::Axes12, [0.001, -0.0004, 0.0006], || false,
    ).unwrap();
    let mut prescribed = initial.strain();
    for i in 0..6 {
        if mask[i] { prescribed[i] = 0.0; } else { prescribed[i] *= 1.2; }
    }
    let report = control.advance(&mut initial, prescribed, 0.2, || false).unwrap();
    assert!(report.dissipated > 0.0);
    assert!(initial.branch_stored_energy() > 0.0);
    near(report.balance_residual, 0.0, initial.stored_energy() + report.work.abs());
    // Equal coefficients do not substitute for the template's original Arc.
    let mut foreign = original.history().initial_relaxed_plane_stress(
        MaterialPlane::Axes12, [0.001, -0.0004, 0.0006], || false,
    ).unwrap();
    let before = foreign.clone();
    assert!(control.advance(&mut foreign, prescribed, 0.2, || false).is_err());
    assert_eq!(foreign, before);
}

#[test]
fn full_initial_state_preserves_tensor_shear_energy_without_initial_heat() {
    let original = law(true);
    let virgin = original.history();
    let strain = [0.001, -0.0003, 0.0002, 0.0004, -0.0005, 0.0006];
    let initial = virgin.initial_relaxed_strain(strain, || false).unwrap();
    let e: Voigt = core::array::from_fn(|i| WEIGHTS[i] * strain[i]);
    let stress: Voigt = core::array::from_fn(|i| {
        (0..6).map(|j| original.equilibrium_mandel()[i][j] * e[j]).sum()
    });
    let energy = 0.5 * (0..6).map(|i| stress[i] * e[i]).sum::<f64>();
    for i in 0..6 {
        near(initial.strain()[i], strain[i], 0.001);
        near(initial.stress()[i], stress[i] / WEIGHTS[i], 1.0e7);
    }
    near(initial.stored_energy(), energy, energy);
    assert_eq!(initial.branch_stored_energy(), 0.0);
    assert_eq!(initial.work(), 0.0);
    assert_eq!(initial.dissipated(), 0.0);
}

#[test]
fn initialization_refuses_state_reset_dc_extrapolation_full_strain_and_cancellation() {
    let original = law(true);
    let virgin = original.history();
    let before = virgin.clone();
    for input in [[f64::NAN, 0.0, 0.0], [0.03, 0.0, 0.0]] {
        assert!(virgin.initial_relaxed_plane_stress(MaterialPlane::Axes12, input, || false).is_err());
    }
    assert!(virgin.initial_relaxed_strain([f64::INFINITY; 6], || false).is_err());
    let mut checks = 0;
    assert_eq!(
        virgin.initial_relaxed_plane_stress(MaterialPlane::Axes12, [0.001, 0.0, 0.0], || {
            checks += 1; checks >= 4
        }).unwrap_err(),
        TensorPronyError::Cancelled,
    );
    assert_eq!(virgin, before);
    let mut evolved = virgin.clone();
    evolved.advance([0.0001; 6], 0.1, || false).unwrap();
    let before = evolved.clone();
    assert!(evolved.initial_relaxed_strain([0.0; 6], || false).is_err());
    assert!(evolved.initial_relaxed_plane_stress(MaterialPlane::Axes12, [0.0; 3], || false).is_err());
    assert_eq!(evolved, before);
    let ac = original.with_band((1.0, 20.0)).unwrap().history();
    assert!(ac.initial_relaxed_plane_stress(MaterialPlane::Axes12, [0.0; 3], || false).is_ok());
    assert!(ac.initial_relaxed_strain([0.0; 6], || false).is_ok());
    assert!(ac.initial_relaxed_plane_stress(MaterialPlane::Axes12, [0.001, 0.0, 0.0], || false).is_err());
    assert!(ac.initial_relaxed_strain([0.0001; 6], || false).is_err());
}
