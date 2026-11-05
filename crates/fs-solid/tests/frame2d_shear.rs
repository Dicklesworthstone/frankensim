//! G1/G3: original source-bound transverse shear across native frame owners.
//! Synthetic coefficients exercise mechanics; they do not qualify materials.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseControls;
use fs_qty::{Density, Frequency, Pressure};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh};

fn material(density: f64, young: f64, shear: f64, restricted: bool) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", density, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
        ("shear_modulus", shear, Pressure::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained();
        let value = if restricted && name == "shear_modulus" {
            validity = validity.with("omega", 0.0, 100.0);
            PropertyValue::Curve {
                abscissa: "omega".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, value), (100.0, value * 0.5)],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if restricted && name == "shear_modulus" {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic frame shear".into(),
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
    let point = QueryPoint::new().with("omega", 50.0).unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn section(shear: f64) -> ConcentricShearSection {
    ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.08,
            material: material(1000.0, 2.0e9, shear, false),
        }],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.85)
    .unwrap()
}

fn line_spin(section: &ConcentricShearSection) -> f64 {
    section
        .elastic()
        .layer_integrals()
        .iter()
        .map(|layer| layer.density_kg_m3() * layer.second_moment_m4())
        .sum()
}

fn frame(length: f64, direction: (f64, f64), source: &ConcentricShearSection) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(length * direction.0, length * direction.1);
    frame.fix(0, [true; 3]);
    frame
        .add_material_shear_element(0, 1, 5, source.clone())
        .unwrap();
    frame
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0e-14),
        "actual={actual}, expected={expected}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 22,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

#[test]
fn cantilever_shear_compliance_work_and_reversal_keep_original_layer_sources() {
    let source = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: 0.04,
                material: material(1000.0, 2e9, 0.8e9, false),
            },
            ConcentricMaterialLayer {
                outer_radius_m: 0.08,
                material: material(2000.0, 4e9, 1.6e9, false),
            },
        ],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.85)
    .unwrap();
    let ei = source.elastic().bending_stiffness_n_m2();
    let ga = source.shear_stiffness_n();
    for length in [0.12_f64, 0.6] {
        for (c, s) in [(1.0, 0.0), (0.6, 0.8)] {
            let mut frame = frame(length, (c, s), &source);
            for force in [1.0, -0.7, 0.2] {
                let load = [0.0, 0.0, 0.0, -s * force, c * force, 0.0];
                frame.static_load(&load, 1).unwrap();
                let u = frame.displacements();
                let transverse = -s * u[3] + c * u[4];
                near(
                    transverse,
                    force * (length.powi(3) / (3.0 * ei) + length / ga),
                    2e-10,
                );
                near(u[5], force * length * length / (2.0 * ei), 2e-10);
                let shear = frame.element_shear_response(0).unwrap();
                near(shear.force_n, force, 2e-10);
                near(shear.engineering_strain, force / ga, 2e-10);
                near(
                    shear.stored_energy_j,
                    force * force * length / (2.0 * ga),
                    4e-10,
                );
                // External work equals independent bending plus shear storage.
                near(
                    0.5 * force * transverse,
                    force * force * length.powi(3) / (6.0 * ei) + shear.stored_energy_j,
                    3e-10,
                );
            }
            assert_eq!(frame.material_shear_section(0), Some(&source));
            assert_eq!(frame.material_section(0), Some(source.elastic()));
        }
    }
}

#[test]
fn pure_bending_and_axial_load_do_not_invent_shear_energy() {
    let source = section(0.8e9);
    let length: f64 = 0.3;
    let mut frame = frame(length, (1.0, 0.0), &source);
    let (axial, moment) = (0.4, 0.2);
    frame
        .static_load(&[0.0, 0.0, 0.0, axial, 0.0, moment], 1)
        .unwrap();
    near(
        frame.displacements()[3],
        axial * length / source.elastic().axial_stiffness_n(),
        2e-10,
    );
    near(
        frame.displacements()[4],
        moment * length * length / (2.0 * source.elastic().bending_stiffness_n_m2()),
        2e-10,
    );
    near(
        frame.displacements()[5],
        moment * length / source.elastic().bending_stiffness_n_m2(),
        2e-10,
    );
    let shear = frame.element_shear_response(0).unwrap();
    assert!(shear.force_n.abs() < 1e-12);
    assert!(shear.stored_energy_j < 1e-28);
}

#[test]
fn shear_corrected_inertia_retains_rigid_mass_and_section_spin_in_any_orientation() {
    for modulus in [0.02e9, 0.8e9, 1000e9] {
        let source = section(modulus);
        let length: f64 = 0.7;
        let (mass, spin) = (source.elastic().line_density_kg_m(), line_spin(&source));
        for (c, s) in [(1.0, 0.0), (0.6, 0.8)] {
            let frame = frame(length, (c, s), &source);
            for (vx, vy, omega) in [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (-0.7, 0.3, 1.4)] {
                let v = [
                    vx,
                    vy,
                    omega,
                    vx - omega * length * s,
                    vy + omega * length * c,
                    omega,
                ];
                let momentum = frame.inertial_forces(&v).unwrap();
                let actual: f64 = v.iter().zip(momentum).map(|(v, p)| v * p).sum();
                let expected = mass * length * (vx * vx + vy * vy)
                    + mass * omega * length * length * (-vx * s + vy * c)
                    + omega * omega * (mass * length.powi(3) / 3.0 + spin * length);
                near(actual, expected, 3e-14);
            }
        }
    }
}

#[test]
fn one_native_shear_pencil_drives_modes_harmonic_response_and_newmark() {
    let source = section(0.8e9);
    let length: f64 = 0.3;
    let mut frame = frame(length, (1.0, 0.0), &source);
    frame.fix(1, [true, false, true]);
    let ei = source.elastic().bending_stiffness_n_m2();
    let beta = 1.0 / (1.0 + 12.0 * ei / (source.shear_stiffness_n() * length * length));
    let k = 12.0 * ei * beta / length.powi(3);
    // Integrals of H=(1-beta)*xi+beta*(3*xi²-2*xi³) and the independent
    // section rotation 6*beta*xi*(1-xi)/L, with both endpoint angles fixed.
    let m = source.elastic().line_density_kg_m()
        * length
        * (1.0 / 3.0 + beta / 30.0 + beta * beta / 210.0)
        + 6.0 * line_spin(&source) * beta * beta / (5.0 * length);
    let omega = (k / m).sqrt();
    near(frame.natural_frequencies().unwrap()[0], omega, 3e-12);
    with_cx(|cx| {
        for frequency in [0.0, 0.6 * omega, 1.3 * omega] {
            let response = frame
                .prepare_harmonic_response(
                    frequency,
                    Rayleigh { a0: 0.0, a1: 0.0 },
                    HarmonicResponseControls {
                        maximum_coordinates: 6,
                        maximum_payload_bytes: 1 << 18,
                        maximum_relative_residual: 1e-12,
                    },
                    cx,
                )
                .unwrap()
                .solve(&[C64::ONE], cx)
                .unwrap();
            near(
                response.displacement[0].re,
                1.0 / (k - frequency * frequency * m),
                3e-12,
            );
            assert!(response.displacement[0].im.abs() < 1e-24);
        }
    });
    let history = frame
        .newmark_with_load(
            &vec![0.0; 501],
            0.002 / omega,
            (0.0, 0.0),
            Rayleigh { a0: 0.0, a1: 0.0 },
            &[0.0, 0.0, 0.0, 0.0, 0.01, 0.0],
        )
        .unwrap();
    for step in [100, 250, 500] {
        near(
            history.u[step][4],
            0.01 / k * (1.0 - (omega * history.t[step]).cos()),
            1e-6,
        );
    }
    assert!(history.energy_balance_error() < 1e-8);
    let shear = frame.element_shear_response(0).unwrap();
    near(shear.force_n, k * history.u[500][4], 3e-10);
}

#[test]
fn simply_supported_spectrum_converges_to_the_timoshenko_continuum() {
    let source = section(0.8e9);
    let length: f64 = 0.8;
    let (ei, ga) = (
        source.elastic().bending_stiffness_n_m2(),
        source.shear_stiffness_n(),
    );
    let (m, j) = (source.elastic().line_density_kg_m(), line_spin(&source));
    // w=A sin(kx), theta=B cos(kx). The two-field continuum pencil is
    // [[GA k², -GA k],[-GA k, EI k²+GA]] versus diag(m,J).
    let wave = core::f64::consts::PI / length;
    let trace = ga * wave * wave / m + (ei * wave * wave + ga) / j;
    let determinant = ga * ei * wave.powi(4) / (m * j);
    let exact =
        (2.0 * determinant / (trace + (trace * trace - 4.0 * determinant).sqrt())).sqrt();
    let mut errors = Vec::new();
    for count in [2, 4, 8] {
        let mut frame =
            Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
        for node in 0..=count {
            frame.add_node(length * node as f64 / count as f64, 0.0);
            frame.fix(node, [true, node == 0 || node == count, false]);
            if node > 0 {
                frame
                    .add_material_shear_element(node - 1, node, 3, source.clone())
                    .unwrap();
            }
        }
        errors.push((frame.natural_frequencies().unwrap()[0] / exact - 1.0).abs());
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn rigid_shear_limit_recovers_euler_bernoulli_without_changing_existing_members() {
    let source = section(1.0e21);
    let mut shear = frame(0.3, (1.0, 0.0), &source);
    let mut original =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    original.add_node(0.0, 0.0);
    original.add_node(0.3, 0.0);
    original.fix(0, [true; 3]);
    original
        .add_material_element(0, 1, 5, source.elastic().clone())
        .unwrap();
    let load = [0.0, 0.0, 0.0, 0.2, 0.4, -0.03];
    original.static_load(&load, 1).unwrap();
    shear.static_load(&load, 1).unwrap();
    for (a, b) in original.displacements().iter().zip(shear.displacements()) {
        near(*a, *b, 2e-11);
    }
    for (a, b) in original
        .natural_frequencies()
        .unwrap()
        .iter()
        .zip(shear.natural_frequencies().unwrap())
    {
        near(*a, b, 2e-11);
    }
    assert!(original.material_shear_section(0).is_none());
    assert!(original.element_shear_response(0).is_err());
    assert!(original.reference_strain_traces(0).is_ok());
    assert!(shear.reference_strain_traces(0).is_err());
}

#[test]
fn refused_shear_sources_and_members_leave_native_state_and_inventory_unchanged() {
    let source = section(0.8e9);
    let mut frame = frame(0.3, (1.0, 0.0), &source);
    let acceleration = [0.2, 1.0, 0.4, -0.3, 1.2, 0.6];
    let before = frame.inertial_forces(&acceleration).unwrap();
    let restricted = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.08,
            material: material(1000.0, 2e9, 0.8e9, true),
        }],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.85)
    .unwrap();
    assert!(
        frame
            .add_material_shear_element(0, 1, 5, restricted)
            .is_err()
    );
    for (i, j, np) in [(0, 1, 2), (0, 0, 5), (0, 2, 5), (usize::MAX, 1, 5)] {
        assert!(
            frame
                .add_material_shear_element(i, j, np, source.clone())
                .is_err()
        );
    }
    assert_eq!(frame.inertial_forces(&acceleration).unwrap(), before);
    assert_eq!(frame.material_shear_section(0), Some(&source));
    assert_eq!(frame.material_shear_section(1), None);
    assert!(frame.element_shear_response(1).is_err());
    assert_eq!(frame.displacements(), &[0.0; 6]);
}
