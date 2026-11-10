//! Lame's thick cylinder against its textbook closed forms, the printed continuity amounts at the
//! membrane gate, the routing in `base_patch`, and the geometry refusals. The fixture is a synthetic
//! card carrying MIL-HDBK-698A's brass density and Young modulus with the derived Poisson ratio
//! E/(2G) - 1 = 1/3; it is not an instrument.
use super::super::super::{source::WallMaterialLossSource, AxialConstraint, LocalMaterialWall, WallMaterialLoss};
use super::thick_cylinder;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};

const BAND: [f64; 2] = [20.0, 4000.0];
const RHO: f64 = 8525.4;
const YOUNG: f64 = 110316e6;
const NU: f64 = 1.0 / 3.0;
const BOTH: [AxialConstraint; 2] = [AxialConstraint::StressFree, AxialConstraint::StrainFixed];

fn brass(nu: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, RHO),
        ("young_modulus", Pressure::DIMS, YOUNG),
        ("poisson_ratio", Dims::NONE, nu),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic brass wall card; not an instrument material".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-brass-wall".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let elastic = WallMaterialLossSource::read(&mut "elastic".split_ascii_whitespace()).unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &elastic.requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn wall(thickness: f64, axial: AxialConstraint) -> LocalMaterialWall {
    LocalMaterialWall::new(brass(NU), thickness, axial, WallMaterialLoss::Elastic, BAND).unwrap()
}

/// The membrane law's `(stiffness_per_area / E, surface_density)` as `base_patch` writes it.
fn membrane(a: f64, h: f64, axial: AxialConstraint) -> (f64, f64) {
    let axial_factor = match axial {
        AxialConstraint::StressFree => 1.0,
        AxialConstraint::StrainFixed => 1.0 / (1.0 - NU * NU),
    };
    ((h / a) / a * axial_factor, RHO * h)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:e} != {expected:e} (relative tolerance {tolerance:e})"
    );
}

/// Timoshenko and Goodier, inner-face displacement per unit pressure on a free outer face.
fn textbook_inner_displacement(a: f64, b: f64, axial: AxialConstraint) -> f64 {
    match axial {
        AxialConstraint::StressFree => a / YOUNG * ((a * a + b * b) / (b * b - a * a) + NU),
        AxialConstraint::StrainFixed => {
            (1.0 + NU) * a / YOUNG * ((1.0 - 2.0 * NU) * a * a + b * b) / (b * b - a * a)
        }
    }
}

#[test]
fn lame_stiffness_matches_the_textbook_plane_stress_and_plane_strain_forms_at_two_ratios() {
    let a = 0.01;
    for ratio in [0.1, 0.25] {
        let h = ratio * a;
        for axial in BOTH {
            let (stiffness_over_young, _) = thick_cylinder(a, h, RHO, NU, axial);
            let expected = 1.0 / textbook_inner_displacement(a, a + h, axial);
            close(stiffness_over_young * YOUNG, expected, 1e-12);
        }
    }
}

#[test]
fn lame_field_meets_both_boundary_conditions_and_its_kinetic_energy_matches_quadrature() {
    let (a, h) = (0.01, 0.0025);
    let b = a + h;
    for axial in BOTH {
        // Independent route: the 3-D constitutive law with the axial condition imposed, then the
        // two boundary conditions, then Simpson quadrature of the kinetic-energy integral.
        let (alpha, beta) = match axial {
            AxialConstraint::StressFree => {
                let c = YOUNG / (1.0 - NU * NU);
                (c * (1.0 + NU), c * (1.0 - NU))
            }
            AxialConstraint::StrainFixed => {
                let k = YOUNG / ((1.0 + NU) * (1.0 - 2.0 * NU));
                (k, k * (1.0 - 2.0 * NU))
            }
        };
        let big_b = a * a * b * b / (beta * (b * b - a * a));
        let big_a = beta * big_b / (alpha * b * b);
        let sigma_r = |r: f64| alpha * big_a - beta * big_b / (r * r);
        close(sigma_r(a), -1.0, 1e-12);
        assert!(sigma_r(b).abs() <= 1e-9, "free outer face: {}", sigma_r(b));
        let inner = big_a * a + big_b / a;
        let n = 4000;
        let step = h / n as f64;
        let mut integral = 0.0;
        for i in 0..=n {
            let r = a + step * i as f64;
            let weight = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
            integral += weight * (big_a * r + big_b / r).powi(2) * r;
        }
        integral *= step / 3.0;
        let (stiffness_over_young, mass) = thick_cylinder(a, h, RHO, NU, axial);
        close(stiffness_over_young * YOUNG, 1.0 / inner, 1e-12);
        close(mass, RHO * integral / (a * inner * inner), 1e-9);
    }
}

#[test]
fn membrane_and_lame_stiffness_differ_by_4_06_percent_and_mass_by_0_68_percent_stress_free_at_the_gate() {
    let (a, h) = (1.0, 0.05);
    let (lame_k, lame_m) = thick_cylinder(a, h, RHO, NU, AxialConstraint::StressFree);
    let (membrane_k, membrane_m) = membrane(a, h, AxialConstraint::StressFree);
    assert!((lame_k / membrane_k - 0.959438).abs() < 2e-6, "{}", lame_k / membrane_k);
    assert!((lame_m / membrane_m - 1.006828).abs() < 2e-6, "{}", lame_m / membrane_m);
}

#[test]
fn membrane_and_lame_stiffness_differ_by_4_82_percent_and_mass_by_0_12_percent_strain_fixed_at_the_gate() {
    let (a, h) = (1.0, 0.05);
    let (lame_k, lame_m) = thick_cylinder(a, h, RHO, NU, AxialConstraint::StrainFixed);
    let (membrane_k, membrane_m) = membrane(a, h, AxialConstraint::StrainFixed);
    assert!((lame_k / membrane_k - 0.951828).abs() < 2e-6, "{}", lame_k / membrane_k);
    assert!((lame_m / membrane_m - 0.998825).abs() < 2e-6, "{}", lame_m / membrane_m);
}

#[test]
fn lame_tends_to_the_membrane_law_as_thickness_over_radius_goes_to_zero() {
    for axial in BOTH {
        let mut previous = f64::INFINITY;
        for ratio in [1e-3, 1e-4, 1e-6] {
            let (lame_k, lame_m) = thick_cylinder(1.0, ratio, RHO, NU, axial);
            let (membrane_k, membrane_m) = membrane(1.0, ratio, axial);
            let gap = (lame_k / membrane_k - 1.0).abs().max((lame_m / membrane_m - 1.0).abs());
            assert!(gap < previous && gap <= 2.0 * ratio, "ratio {ratio}: gap {gap:e}");
            previous = gap;
        }
    }
}

#[test]
fn base_patch_runs_the_membrane_expressions_unchanged_under_the_gate_and_lame_over_it() {
    let h = 0.0007;
    let w = wall(h, AxialConstraint::StressFree);
    // 0.0007 / 0.015 = 0.0467: under the gate, the original two expressions, bit for bit.
    let r = 0.015;
    let pin = w.base_patch(r, 1.0).unwrap().wall;
    assert_eq!(pin.stiffness_per_area, YOUNG * ((h / r) / r * 1.0));
    assert_eq!(pin.surface_density, RHO * h);
    assert_eq!(pin.resistance, 0.0);
    // 0.0007 / 0.00625 = 0.112: over the gate, the thick law (the saxophone's mouthpiece bore).
    let r = 0.00625;
    let pin = w.base_patch(r, 1.0).unwrap().wall;
    let (stiffness_over_young, mass) = thick_cylinder(r, h, RHO, NU, AxialConstraint::StressFree);
    assert_eq!(pin.stiffness_per_area, YOUNG * stiffness_over_young);
    assert_eq!(pin.surface_density, mass);
    assert_eq!(pin.resistance, 0.0);
    assert!(pin.stiffness_per_area < YOUNG * ((h / r) / r), "Lame is softer than the membrane here");
}

#[test]
fn non_positive_thickness_radius_or_length_refuse() {
    for thickness in [0.0, -1e-3, f64::NAN] {
        assert!(
            LocalMaterialWall::new(brass(NU), thickness, AxialConstraint::StressFree, WallMaterialLoss::Elastic, BAND)
                .is_err(),
            "thickness {thickness}"
        );
    }
    let w = wall(0.0007, AxialConstraint::StressFree);
    for (radius, length) in [(0.0, 1.0), (-0.006, 1.0), (0.006, 0.0), (f64::NAN, 1.0), (0.006, f64::INFINITY)] {
        assert!(w.base_patch(radius, length).is_err(), "radius {radius} length {length}");
    }
}
