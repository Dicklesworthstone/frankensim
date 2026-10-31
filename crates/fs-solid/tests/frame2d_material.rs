//! G1/G3: source data reaches the actual static, gravity, modal and time owners.
//! Constitutive data are synthetic and do not qualify a measured material.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Frequency, Pressure, Temperature};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh};

fn material(
    density: f64,
    young: f64,
    temperature: f64,
    frequency: Option<bool>,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", density, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        let value = if name == "young_modulus" && frequency == Some(true) {
            validity = validity.with("omega", 0.0, 100.0);
            PropertyValue::Curve {
                abscissa: "omega".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, value), (100.0, 0.5 * value)],
                dims,
            }
        } else {
            if name == "young_modulus" && frequency == Some(false) {
                validity = validity.with("omega", 0.0, 100.0);
            }
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value), (400.0, 0.5 * value)],
                dims,
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic frame substitution".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic frame solid".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut point = QueryPoint::new().with("T", temperature).unwrap();
    if frequency.is_some() {
        point = point.with("omega", 50.0).unwrap();
    }
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn section(density: f64, young: f64) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.02,
            material: material(density, young, 300.0, None),
        }],
    )
    .unwrap()
}

fn cantilever(sections: &[ConcentricElasticSection]) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let mut previous = frame.add_node(0.0, 0.0);
    frame.fix(previous, [true; 3]);
    for (i, source) in sections.iter().enumerate() {
        let next = frame.add_node((i + 1) as f64 / sections.len() as f64, 0.0);
        frame
            .add_material_element(previous, next, 5, source.clone())
            .unwrap();
        previous = next;
    }
    frame
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual}, expected={expected}"
    );
}

#[test]
fn data_only_stiffness_and_density_swaps_change_compliance_weight_and_modes() {
    let a = std::f64::consts::PI * 0.02_f64.powi(2);
    let i = std::f64::consts::PI * 0.02_f64.powi(4) / 4.0;
    let mut results = Vec::new();
    for (density, young) in [(1000.0, 2e9), (4000.0, 2e9), (1000.0, 8e9)] {
        let source = section(density, young);
        let mut frame = cantilever(&vec![source.clone(); 4]);
        let modes = frame.natural_frequencies().unwrap();
        let gravity: Vec<_> = (0..frame.ndof())
            .map(|d| if d % 3 == 1 { -9.81 } else { 0.0 })
            .collect();
        let weight = frame.inertial_forces(&gravity).unwrap();
        near(
            weight.iter().skip(1).step_by(3).sum(),
            -density * a * 9.81,
            1e-13,
        );
        let mut force = vec![0.0; frame.ndof()];
        force[12] = 0.02;
        force[13] = 0.01;
        frame.static_load(&force, 1).unwrap();
        near(frame.displacements()[12], 0.02 / (young * a), 1e-10);
        near(frame.displacements()[13], 0.01 / (3.0 * young * i), 1e-10);
        assert_eq!(frame.material_section(0), Some(&source));
        assert_eq!(frame.material_section(3), Some(&source));
        results.push((frame.displacements()[13], modes));
    }
    near(results[1].0, results[0].0, 1e-12);
    near(results[2].0, results[0].0 / 4.0, 1e-12);
    for mode in 0..results[0].1.len() {
        near(results[1].1[mode], 0.5 * results[0].1[mode], 1e-8);
        near(results[2].1[mode], 2.0 * results[0].1[mode], 1e-8);
    }
}

#[test]
fn layered_and_spatial_assignments_keep_their_individual_material_response() {
    let inner = 0.01_f64;
    let outer = 0.02_f64;
    let core = material(1000.0, 2e9, 300.0, None);
    // The second material is queried at its own physical state, not averaged T.
    let skin = material(6000.0, 8e9, 400.0, None);
    let composite = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: inner,
                material: core.clone(),
            },
            ConcentricMaterialLayer {
                outer_radius_m: outer,
                material: skin.clone(),
            },
        ],
    )
    .unwrap();
    let soft = section(500.0, 1e9);
    let mut frame = cantilever(&[composite.clone(), soft]);
    let pi = std::f64::consts::PI;
    let core_area = pi * inner.powi(2);
    let skin_area = pi * (outer.powi(2) - inner.powi(2));
    let ea = 2e9 * core_area + 4e9 * skin_area;
    let ei = pi / 4.0 * (2e9 * inner.powi(4) + 4e9 * (outer.powi(4) - inner.powi(4)));
    let soft_ea = 1e9 * pi * outer.powi(2);
    let soft_ei = 1e9 * pi / 4.0 * outer.powi(4);
    let mut force = vec![0.0; frame.ndof()];
    force[6] = 0.02;
    force[7] = 0.01;
    frame.static_load(&force, 1).unwrap();
    near(
        frame.displacements()[6],
        0.02 * 0.5 * (ea.recip() + soft_ea.recip()),
        1e-10,
    );
    near(
        frame.displacements()[7],
        0.01 / 3.0 * (0.875 / ei + 0.125 / soft_ei),
        1e-10,
    );
    let source = frame.material_section(0).unwrap();
    assert_eq!(source.layers()[0].material, core);
    assert_eq!(source.layers()[1].material, skin);
    let inertia = frame
        .inertial_forces(&[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        .unwrap();
    near(
        inertia.iter().step_by(3).sum(),
        0.5 * (1000.0 * core_area + 3000.0 * skin_area + 500.0 * pi * outer.powi(2)),
        1e-13,
    );
}

#[test]
fn material_inertia_and_stiffness_drive_the_same_newmark_trajectory() {
    for (density, young) in [(1000.0, 2e9), (4000.0, 2e9), (1000.0, 8e9)] {
        let mut frame = cantilever(&[section(density, young)]);
        // One free bending coordinate: guided tip, clamped base. The exact
        // FE oscillator has k=12EI/L^3 and
        // m=156*rho*A*L/420 + 6*rho*I/(5*L), including section spin.
        frame.fix(1, [true, false, true]);
        let pi = std::f64::consts::PI;
        let k = 12.0 * young * pi * 0.02_f64.powi(4) / 4.0;
        let m = 156.0 / 420.0 * density * pi * 0.02_f64.powi(2)
            + 6.0 / 5.0 * density * pi * 0.02_f64.powi(4) / 4.0;
        let omega = (k / m).sqrt();
        let dt = 0.002 / omega;
        let mut force = vec![0.0; frame.ndof()];
        force[4] = 0.01;
        let history = frame
            .newmark_with_load(
                &vec![0.0; 501],
                dt,
                (0.0, 0.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
                &force,
            )
            .unwrap();
        for step in [100, 250, 500] {
            let exact = 0.01 / k * (1.0 - (omega * history.t[step]).cos());
            near(history.u[step][4], exact, 1e-6);
        }
        assert!(history.energy_balance_error() < 1e-8);
    }
}

#[test]
fn g1_layered_material_includes_each_density_weighted_section_spin() {
    let (inner, outer, length) = (0.01_f64, 0.02_f64, 0.3);
    let source = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: inner,
                material: material(1000.0, 2e9, 300.0, None),
            },
            ConcentricMaterialLayer {
                outer_radius_m: outer,
                material: material(6000.0, 8e9, 400.0, None),
            },
        ],
    )
    .unwrap();
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(length, 0.0);
    frame.add_material_element(0, 1, 5, source.clone()).unwrap();
    let velocity = [0.0, 0.0, 1.0, 0.0, length, 1.0];
    let momentum = frame.inertial_forces(&velocity).unwrap();
    let pi = std::f64::consts::PI;
    // The outer layer's density at 400 K is 3000, not the 300 K value 6000.
    let line_mass = pi * (1000.0 * inner.powi(2) + 3000.0 * (outer.powi(2) - inner.powi(2)));
    let line_spin = pi / 4.0 * (1000.0 * inner.powi(4) + 3000.0 * (outer.powi(4) - inner.powi(4)));
    let actual: f64 = velocity.iter().zip(momentum).map(|(v, p)| v * p).sum();
    near(
        actual,
        line_mass * length.powi(3) / 3.0 + line_spin * length,
        1e-12,
    );
    assert_eq!(frame.material_section(0), Some(&source));
}

#[test]
fn unsupported_sources_and_bad_geometry_do_not_publish_partial_members_or_mass() {
    let source = section(1000.0, 2e9);
    let mut frame = cantilever(&[source.clone()]);
    let acceleration = [1.0, -1.0, 0.0, 1.0, -1.0, 0.0];
    let before = frame.inertial_forces(&acceleration).unwrap();
    for frequency in [Some(false), Some(true)] {
        let unsupported = ConcentricElasticSection::try_new(
            0.0,
            vec![ConcentricMaterialLayer {
                outer_radius_m: 0.02,
                material: material(1000.0, 2e9, 300.0, frequency),
            }],
        )
        .unwrap();
        assert!(frame.add_material_element(0, 1, 5, unsupported).is_err());
    }
    for (a, b, np) in [(0, 1, 2), (0, 0, 5), (0, 2, 5), (usize::MAX, 1, 5)] {
        assert!(
            frame
                .add_material_element(a, b, np, source.clone())
                .is_err()
        );
    }
    assert_eq!(frame.inertial_forces(&acceleration).unwrap(), before);
    assert_eq!(frame.material_section(0), Some(&source));
    assert_eq!(frame.material_section(1), None);
    assert_eq!(frame.displacements(), &[0.0; 6]);

    let invalid = frame.add_node(f64::INFINITY, 0.0);
    assert!(frame.add_material_element(0, invalid, 5, source).is_err());
    assert_eq!(frame.material_section(1), None);
}
