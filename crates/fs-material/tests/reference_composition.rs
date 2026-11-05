//! G0/G3: ordered natural maps, work-conjugate pullbacks and frame covariance.
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

type Matrix = [f64; 9];

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    core::array::from_fn(|ij| {
        let (i, j) = (ij / 3, ij % 3);
        a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j]
    })
}

fn transpose(a: &Matrix) -> Matrix {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6.0,
            c01: 2.0,
            kappa: 80.0,
        },
        5.0,
    )
    .unwrap()
}

fn maps() -> (ReferenceDeformation, ReferenceDeformation) {
    let manufactured =
        ReferenceDeformation::try_new([1.1, 0.23, 0.04, 0.0, 0.93, -0.08, 0.0, 0.0, 1.05]).unwrap();
    let thermal = ReferenceDeformation::from_principal_stretches(
        [1.03, 0.97, 1.02],
        [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    (manufactured, thermal)
}

#[test]
fn ordered_composition_preserves_relaxed_state_volume_and_identity() {
    let (a, b) = maps();
    let ab = a.followed_by(b).unwrap();
    let ba = b.followed_by(a).unwrap();
    assert_ne!(ab.deformation_gradient(), ba.deformation_gradient());
    for (actual, expected) in ab
        .deformation_gradient()
        .iter()
        .zip(product(b.deformation_gradient(), a.deformation_gradient()))
    {
        close(*actual, expected);
    }
    close(ab.volume_ratio(), a.volume_ratio() * b.volume_ratio());
    assert_eq!(
        ab.elastic_deformation(ab.deformation_gradient()),
        *ReferenceDeformation::identity().deformation_gradient()
    );
    assert_eq!(a.followed_by(ReferenceDeformation::identity()).unwrap(), a);
    assert_eq!(ReferenceDeformation::identity().followed_by(b).unwrap(), b);
    let material = law();
    let relaxed = ReferencedHyperelastic::new(&material, ab);
    close(relaxed.energy_density(ab.deformation_gradient()), 0.0);
    for p in relaxed.first_piola(ab.deformation_gradient()).unwrap() {
        close(p, 0.0);
    }
}

#[test]
fn composed_energy_piola_and_tangent_match_successive_reference_pullbacks() {
    let material = law();
    let (a, b) = maps();
    let f = [1.19, 0.07, -0.03, 0.02, 0.91, 0.04, 0.01, -0.02, 1.07];
    // The inner map B and outer map A give F A^-1 B^-1 = F (B A)^-1.
    let inner = ReferencedHyperelastic::new(&material, b);
    let nested = ReferencedHyperelastic::new(&inner, a);
    let composed = ReferencedHyperelastic::new(&material, a.followed_by(b).unwrap());
    close(composed.energy_density(&f), nested.energy_density(&f));
    let pc = composed.first_piola(&f).unwrap();
    let pn = nested.first_piola(&f).unwrap();
    let kc = composed.consistent_tangent(&f).unwrap();
    let kn = nested.consistent_tangent(&f).unwrap();
    for i in 0..9 {
        close(pc[i], pn[i]);
        for j in 0..9 {
            close(kc[i][j], kn[i][j]);
        }
    }
}

#[test]
fn rotating_both_declared_maps_preserves_the_composed_physical_response() {
    let material = law();
    let (a, b) = maps();
    let q = [0.36, -0.48, 0.8, 0.8, 0.6, 0.0, -0.48, 0.64, 0.6];
    let qt = transpose(&q);
    let rotate = |m: &Matrix| product(&product(&q, m), &qt);
    let ar = ReferenceDeformation::try_new(rotate(a.deformation_gradient())).unwrap();
    let br = ReferenceDeformation::try_new(rotate(b.deformation_gradient())).unwrap();
    let original = a.followed_by(b).unwrap();
    let rotated = ar.followed_by(br).unwrap();
    for (actual, expected) in rotated
        .deformation_gradient()
        .iter()
        .zip(rotate(original.deformation_gradient()))
    {
        close(*actual, expected);
    }
    let f = [1.19, 0.07, -0.03, 0.02, 0.91, 0.04, 0.01, -0.02, 1.07];
    let x = ReferencedHyperelastic::new(&material, original);
    let y = ReferencedHyperelastic::new(&material, rotated);
    close(x.energy_density(&f), y.energy_density(&rotate(&f)));
    for (actual, expected) in y
        .first_piola(&rotate(&f))
        .unwrap()
        .iter()
        .zip(rotate(&x.first_piola(&f).unwrap()))
    {
        close(*actual, expected);
    }
}

#[test]
fn unrepresentable_composed_volume_is_refused_without_changing_inputs() {
    let large =
        ReferenceDeformation::try_new([1e100, 0.0, 0.0, 0.0, 1e100, 0.0, 0.0, 0.0, 1e100]).unwrap();
    let before = large;
    assert!(large.followed_by(large).is_err());
    assert_eq!(large, before);
}
