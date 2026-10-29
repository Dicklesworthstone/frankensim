//! Source coating -> original regional shell -> stiffness, inertia and modal loss.
//!
//! The support is an authored rectangular coupon with the existing violin wood
//! material pack. It is not a measured violin or an acoustic-quality ranking.
//! The complete TSV compiler path is checked in xtask/tests/varnish_coating.rs.
#![deny(unsafe_code)]

use fs_couple::shear_shell::{
    CONSTANT_LOSS_REFERENCE_OMEGA, CONSTANT_LOSS_REFERENCE_OMEGA_DIMS, CoatingSide,
    ResolvedShearShell, ShellCoating, ShellRegionMaterial, TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES,
    compile_material_shear_shell, constant_loss::ConstantLossSpectrum,
};
use fs_couple::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    coating::ResolvedCoating,
    plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES,
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_modal::eigh_gen_dense;
use fs_plate::{
    PlateRegion,
    shell::{
        ShellMesh,
        mitc::{
            ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
            assemble_mitc3p, viscosity::project_viscosity,
        },
    },
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use std::f64::consts::TAU;

fn source_card(name: &str, claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "declared material source".into(),
            process: "frozen source realization with authored test geometry".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn wood() -> ResolvedMaterialStatePoint {
    let bytes = include_bytes!(
        "../examples/material_strings/packs/violin-spruce-norway-viala2020-shell.fsmatpk"
    );
    let decoded = NormalizedPack::from_bytes(bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), bytes).unwrap();
    let card = source_card(pack.pack_id(), pack.claims().clone());
    let mut keys = vec![
        ("density", Density::DIMS),
        ("young_modulus_1", Pressure::DIMS),
        ("young_modulus_2", Pressure::DIMS),
        ("poisson_ratio_12", Dims::NONE),
        ("shear_modulus_12", Pressure::DIMS),
        ("shear_modulus_13", Pressure::DIMS),
        ("shear_modulus_23", Pressure::DIMS),
        (
            CONSTANT_LOSS_REFERENCE_OMEGA,
            CONSTANT_LOSS_REFERENCE_OMEGA_DIMS,
        ),
    ];
    keys.extend(PLANE_STRESS_VISCOSITY_PROPERTIES.map(|key| (key, DynViscosity::DIMS)));
    keys.extend(TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES.map(|key| (key, DynViscosity::DIMS)));
    let requirements: Vec<_> = keys
        .into_iter()
        .map(|(key, dims)| {
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    let point = QueryPoint::new()
        .with("T", 293.15)
        .unwrap()
        .with("omega", TAU * 1000.)
        .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn varnish() -> ResolvedCoating {
    // Numerical fixture uses the published O1 coefficients without inventing
    // observations or uncertainty. Separate source-pack integration tests
    // verify their attributed TSV transport and actual conditioning schema.
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 1097.),
        ("varnish_e0", Pressure::DIMS, 20e6),
        ("varnish_einf", Pressure::DIMS, 1.37e9),
        ("varnish_alpha", Dims::NONE, 0.391),
        ("varnish_beta", Dims::NONE, 0.360),
        ("varnish_tau", Dims([0, 0, 1, 0, 0, 0]), 0.127),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 293.15, 293.15)
                    .with("omega", TAU * 20., TAU * 20000.),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source:
                        "numerical fixture: Laemmlein 2021 O1 fit, doi:10.1016/j.rinma.2020.100137"
                            .into(),
                    license: "Attributed numerical facts; no source article redistribution".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
    }
    let card = source_card("modern copal-oil O1 numerical fixture", claims);
    let point = QueryPoint::new()
        .with("T", 293.15)
        .unwrap()
        .with("omega", TAU * 1000.)
        .unwrap();
    let state = resolve_material_state_point(
        &card,
        &point,
        &ResolvedCoating::requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    // Isotropy, real constant nu and the strain limit are explicit assumptions.
    ResolvedCoating::from_material(state, 0.3, 0.0003).unwrap()
}

pub(crate) fn coupon(gate: &CancelGate) -> (ResolvedShearShell, NodalDirectors, Vec<ShellConstraint>) {
    let n = 3;
    let mut nodes = Vec::new();
    let mut tris = Vec::new();
    let mut constraints = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let node = nodes.len();
            nodes.push([0.3 * i as f64 / n as f64, 0.2 * j as f64 / n as f64, 0.]);
            if i == 0 || j == 0 || i == n || j == n {
                constraints.push(ShellConstraint {
                    node,
                    mask: [true; 6],
                    frame: ConstraintFrame::Global,
                });
            }
        }
    }
    for j in 0..n {
        for i in 0..n {
            let a = j * (n + 1) + i;
            tris.extend([[a, a + 1, a + n + 2], [a, a + n + 2, a + n + 1]]);
        }
    }
    let mesh = ShellMesh::new(nodes, tris).unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0., 0., 1.]; mesh.nodes.len()], 0.9).unwrap();
    let budget = MitcBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
    };
    let region = ShellRegionMaterial {
        region: PlateRegionMaterial {
            region: PlateRegion {
                name: "wood".into(),
                triangle_indices: (0..mesh.tris.len()).collect(),
            },
            material: wood(),
            model: PlateMaterialModel::Orthotropic12,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.003),
            material_angle_rad: 0.17,
        },
        shear_correction: 5. / 6.,
        axis1_global: [1., 0., 0.],
    };
    (
        compile_material_shear_shell(mesh, vec![region], budget, gate).unwrap(),
        directors,
        constraints,
    )
}

pub(crate) fn budget(shell: &ResolvedShearShell) -> MitcBudget {
    MitcBudget {
        max_nodes: shell.mesh().nodes.len(),
        max_triangles: shell.mesh().tris.len(),
    }
}

pub(crate) fn layer(side: CoatingSide, thickness_m: f64) -> ShellCoating {
    ShellCoating {
        region: "wood".into(),
        material: varnish(),
        side,
        thickness_m,
        reference_omega_rad_s: TAU * 1000.,
        max_storage_relative_drift: 0.9,
    }
}

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

fn modes(
    shell: &ResolvedShearShell,
    directors: &NodalDirectors,
    constraints: &[ShellConstraint],
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let model = assemble_mitc3p(
        shell.mesh(),
        directors,
        shell.sections(),
        constraints,
        ShellMass::Consistent,
        budget(shell),
        1e-3,
        || false,
    )
    .unwrap();
    let pencil = model.pencil();
    let eigensystem =
        eigh_gen_dense(&pencil.k.to_dense(), &pencil.m.to_dense(), pencil.free).unwrap();
    let selected: Vec<_> = eigensystem
        .into_iter()
        .filter(|p| p.lambda > 0.)
        .take(3)
        .collect();
    assert_eq!(selected.len(), 3);
    for p in &selected {
        assert!(p.residual <= 1e-7 * p.lambda);
    }
    let omega = selected.iter().map(|p| p.lambda.sqrt()).collect();
    let physical = selected
        .iter()
        .map(|p| model.expand(&p.phi).unwrap())
        .collect();
    (omega, physical)
}

fn memory_options() -> fs_couple::shear_shell::coating::relaxation::CoatingRelaxationOptions {
    fs_couple::shear_shell::coating::relaxation::CoatingRelaxationOptions {
        terms: 16, storage_relative_tolerance: 0.001, loss_relative_tolerance: 0.001,
        max_memory_coordinates: 128,
    }
}

#[test]
fn g1_hereditary_coating_uses_equilibrium_once_and_preserves_physical_mass() {
    let gate = CancelGate::new_clock_free();
    let (bare, directors, supports) = coupon(&gate);
    let sliced = bare.with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate).unwrap();
    let shell = sliced.with_coating_relaxation(memory_options(), &gate).unwrap();
    near(shell.mass_kg(), sliced.mass_kg(), 1e-14);
    assert_eq!(shell.coatings()[0].faces_m(), sliced.coatings()[0].faces_m());
    assert_eq!(shell.coatings()[0].reference_storage_pa(), 20e6);
    assert!(sliced.coatings()[0].reference_storage_pa() > 30.0 * 20e6);
    let (omega, physical) = modes(&shell, &directors, &supports);
    let c = shell.project_material_viscosity(&directors, &physical, &omega,
        budget(&shell), usize::MAX, &gate).unwrap();
    let mut wood_matrix = project_viscosity(shell.mesh(), &directors, shell.sections(),
        shell.substrate_viscosities(), &physical, budget(&shell), usize::MAX, || false).unwrap();
    let anchor = shell.regions()[0].constant_loss_reference_omega_rad_s().unwrap();
    let spectrum = ConstantLossSpectrum::fit();
    for i in 0..omega.len() { for j in 0..omega.len() {
        let ratio = (anchor / omega[i] * spectrum.shape(omega[i])
            * anchor / omega[j] * spectrum.shape(omega[j])).sqrt();
        wood_matrix[i * omega.len() + j] *= ratio;
    }}
    for (&a, &b) in c.iter().zip(&wood_matrix) { near(a, b, 2e-13); }
    let rebound = shell.with_region_material("wood", &wood(), &gate).unwrap();
    assert_eq!(rebound.coating_relaxation_options(), Some(memory_options()));
    assert_eq!(rebound.coatings()[0].reference_storage_pa(), 20e6);
    assert_eq!(rebound.coatings()[0].faces_m(), shell.coatings()[0].faces_m());
}

#[test]
fn g1_full_coating_memory_reconstructs_source_storage_and_loss_without_per_film_state_growth() {
    let gate = CancelGate::new_clock_free();
    let (bare, directors, supports) = coupon(&gate);
    let shell = bare.with_coatings(vec![layer(CoatingSide::Positive, 20e-6),
        layer(CoatingSide::Negative, 15e-6)], &gate).unwrap()
        .with_coating_relaxation(memory_options(), &gate).unwrap();
    let (omega, mut physical) = modes(&shell, &directors, &supports);
    // A rotated mass-orthonormal retained basis makes the cross terms live.
    let (a, b) = (physical[0].clone(), physical[1].clone());
    for i in 0..a.len() {
        physical[0][i] = (a[i] + b[i]) / 2.0_f64.sqrt();
        physical[1][i] = (b[i] - a[i]) / 2.0_f64.sqrt();
    }
    let image = shell.project_coating_relaxation(&directors, &physical, &omega,
        budget(&shell), usize::MAX, &gate).unwrap().unwrap();
    assert_eq!(image.films().len(), 2);
    assert_eq!(image.fit_certificates().len(), 1);
    let fit = &image.fit_certificates()[0];
    assert!(image.relaxation_times_s().len() <= omega.len() * fit.model.terms.len());
    let n = omega.len();
    let mut h = vec![0.0; n*n];
    for film in image.films() {
        for (total, value) in h.iter_mut().zip(film.unit_stiffness_per_pa_s2()) { *total += value; }
    }
    assert!(h[1].abs() > 1e-6 * h[0].abs());
    for hz in [100.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0] {
        let w = TAU * hz;
        let (storage, loss) = shell.coatings()[0].input().material.law().modulus(w).unwrap();
        let mut kr: Vec<_> = h.iter().map(|v| 20e6 * v).collect();
        let mut ki = vec![0.0; n*n];
        for (row, tau) in image.coupling_rows_per_s().iter().zip(image.relaxation_times_s()) {
            let x = w*tau;
            let re = x*x / (1.0+x*x); let im = x / (1.0+x*x);
            for i in 0..n { for j in 0..n {
                kr[i*n+j] += row[i]*row[j]*re;
                ki[i*n+j] += row[i]*row[j]*im;
            }}
        }
        for k in 0..n*n {
            near(kr[k], storage*h[k], 0.0011);
            near(ki[k], loss*h[k], 0.0011);
        }
        for vector in [[1.0, -2.0, 0.3], [-3.0, 0.4, 0.1]] {
            let loss: f64 = (0..n).flat_map(|i| (0..n).map(move |j| (i,j)))
                .map(|(i,j)| vector[i] * ki[i*n+j] * vector[j]).sum();
            assert!(loss >= 0.0);
        }
    }
}

#[test]
fn g0_coating_memory_refuses_fit_work_and_state_overruns_without_changing_the_source() {
    let gate = CancelGate::new_clock_free();
    let (bare, directors, supports) = coupon(&gate);
    let sliced = bare.with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate).unwrap();
    let mut options = memory_options(); options.terms = 0;
    assert!(sliced.with_coating_relaxation(options, &gate).is_err());
    options = memory_options(); options.max_memory_coordinates = 1;
    let shell = sliced.with_coating_relaxation(options, &gate).unwrap();
    let (omega, physical) = modes(&shell, &directors, &supports);
    assert!(shell.project_coating_relaxation(&directors, &physical, &omega,
        budget(&shell), usize::MAX, &gate).is_err());
    assert!(shell.project_coating_relaxation(&directors, &physical, &omega,
        budget(&shell), 1, &gate).is_err());
    assert!(sliced.coating_relaxation_options().is_none());
    assert!(sliced.coatings()[0].reference_storage_pa() > 20e6);
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(sliced.with_coating_relaxation(memory_options(), &cancelled).is_err());
}

#[test]
fn g0_coating_preserves_the_wood_reference_and_keeps_signed_b_and_first_inertia() {
    let gate = CancelGate::new_clock_free();
    let (bare, _, _) = coupon(&gate);
    let positive = bare
        .with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate)
        .unwrap();
    let negative = bare
        .with_coatings(vec![layer(CoatingSide::Negative, 20e-6)], &gate)
        .unwrap();
    let repeat = positive
        .with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate)
        .unwrap();
    assert_eq!(positive.mesh(), bare.mesh());
    assert_eq!(
        positive.regions()[0].thickness_m(),
        bare.regions()[0].thickness_m()
    );
    assert_eq!(positive.material_at_element(0), bare.material_at_element(0));
    near(
        positive.mass_kg() - bare.mass_kg(),
        0.3 * 0.2 * 1097. * 20e-6,
        1e-12,
    );
    near(negative.mass_kg(), positive.mass_kg(), 1e-14);
    near(repeat.mass_kg(), positive.mass_kg(), 1e-14);
    assert!(bare.coatings().is_empty());
    assert_eq!(positive.coatings().len(), 1);
    let rebound = positive
        .with_region_material("wood", &wood(), &gate)
        .unwrap();
    assert_eq!(rebound.coatings().len(), positive.coatings().len());
    assert_eq!(
        rebound.coatings()[0].input().material.source(),
        positive.coatings()[0].input().material.source()
    );
    assert_eq!(
        rebound.coatings()[0].faces_m(),
        positive.coatings()[0].faces_m()
    );
    near(rebound.mass_kg(), positive.mass_kg(), 1e-14);
    for (after, before) in rebound.sections().iter().zip(positive.sections()) {
        assert_eq!(after.coupling_pa_m2(), before.coupling_pa_m2());
        assert_eq!(
            after.mass_first_moment_kg_m(),
            before.mass_first_moment_kg_m()
        );
    }
    for (plus, minus) in positive.sections().iter().zip(negative.sections()) {
        assert!(plus.mass_first_moment_kg_m() > 0.);
        near(
            plus.mass_first_moment_kg_m(),
            -minus.mass_first_moment_kg_m(),
            1e-12,
        );
        near(plus.rotary_inertia_kg(), minus.rotary_inertia_kg(), 1e-12);
        for (a, b) in plus.coupling_pa_m2().iter().zip(minus.coupling_pa_m2()) {
            if a.abs().max(b.abs()) > 1e-15 {
                near(*a, -*b, 1e-11);
            }
        }
    }
    let restored = positive.with_coatings(vec![], &gate).unwrap();
    near(restored.mass_kg(), bare.mass_kg(), 1e-14);
    for (after, before) in restored.sections().iter().zip(bare.sections()) {
        assert_eq!(after.principal_matrices(), before.principal_matrices());
        assert_eq!(after.coupling_pa_m2(), before.coupling_pa_m2());
        assert_eq!(
            after.mass_first_moment_kg_m(),
            before.mass_first_moment_kg_m()
        );
    }
}

#[test]
fn g1_real_wood_and_varnish_change_the_original_shell_eigenproblem_and_damping() {
    let gate = CancelGate::new_clock_free();
    let (bare, directors, constraints) = coupon(&gate);
    let coated = bare
        .with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate)
        .unwrap();
    let (bare_omega, _) = modes(&bare, &directors, &constraints);
    let (omega, physical) = modes(&coated, &directors, &constraints);
    assert!((omega[0] / bare_omega[0] - 1.).abs() > 1e-6);
    coated.validate_coating_frequencies(&omega, &gate).unwrap();
    let c = coated
        .project_material_viscosity(
            &directors,
            &physical,
            &omega,
            budget(&coated),
            64_000_000,
            &gate,
        )
        .unwrap();
    let wood_c = project_viscosity(
        coated.mesh(),
        &directors,
        coated.sections(),
        coated.substrate_viscosities(),
        &physical,
        budget(&coated),
        64_000_000,
        || false,
    )
    .unwrap();
    let film_c = project_viscosity(
        coated.mesh(),
        &directors,
        coated.sections(),
        coated.coating_viscosities(0).unwrap(),
        &physical,
        budget(&coated),
        64_000_000,
        || false,
    )
    .unwrap();
    let law = varnish();
    let reference = TAU * 1000.;
    let eta_ref = law.law().modulus(reference).unwrap().1 / reference;
    let factors: Vec<_> = omega
        .iter()
        .map(|&w| (law.law().modulus(w).unwrap().1 / w / eta_ref).sqrt())
        .collect();
    // The wood's original anchored law must not rescale the film spectrum.
    let wood_reference = coated.regions()[0]
        .constant_loss_reference_omega_rad_s()
        .unwrap();
    let spectrum = ConstantLossSpectrum::fit();
    let wood_factors: Vec<_> = omega
        .iter()
        .map(|&w| (wood_reference * spectrum.shape(w) / w).sqrt())
        .collect();
    let n = omega.len();
    for i in 0..n {
        for j in 0..n {
            let expected = wood_factors[i] * wood_c[n * i + j] * wood_factors[j]
                + factors[i] * film_c[n * i + j] * factors[j];
            near(c[n * i + j], expected, 2e-11);
            near(c[n * i + j], c[n * j + i], 2e-11);
        }
    }
    for i in 0..n {
        assert!(c[n * i + i] > wood_factors[i].powi(2) * wood_c[n * i + i]);
    }
    // The complete operator is passive; no off-diagonal replacement is needed.
    let (eigenvalues, _) = fs_la::eigen::jacobi_eigh(&c, n);
    assert!(
        eigenvalues
            .iter()
            .all(|&v| v >= -1e-12 * c.iter().map(|v| v.abs()).fold(0., f64::max))
    );
}

#[test]
fn g3_failed_stack_replacement_and_unsupported_spectrum_leave_original_wood_untouched() {
    let gate = CancelGate::new_clock_free();
    let (bare, directors, constraints) = coupon(&gate);
    let mut input = layer(CoatingSide::Positive, 20e-6);
    input.max_storage_relative_drift = 0.;
    let strict = bare.with_coatings(vec![input.clone()], &gate).unwrap();
    strict
        .validate_coating_frequencies(&[TAU * 1000.], &gate)
        .unwrap();
    assert!(
        strict
            .validate_coating_frequencies(&[TAU * 1100.], &gate)
            .is_err()
    );
    let (omega, physical) = modes(&bare, &directors, &constraints);
    let coated = bare
        .with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &gate)
        .unwrap();
    let mut outside = omega.clone();
    outside[0] = TAU * 10.;
    assert!(
        coated
            .project_material_viscosity(
                &directors,
                &physical,
                &omega,
                budget(&coated),
                32_000_000,
                &gate
            )
            .is_err()
    );
    assert!(
        coated
            .project_material_viscosity(
                &directors,
                &physical,
                &outside,
                budget(&coated),
                64_000_000,
                &gate
            )
            .is_err()
    );
    input.region = "absent".into();
    assert!(coated.with_coatings(vec![input], &gate).is_err());
    assert!(
        coated
            .with_coatings(vec![layer(CoatingSide::Positive, 0.)], &gate)
            .is_err()
    );
    assert!(bare.coatings().is_empty());
    assert_eq!(coated.coatings().len(), 1);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        bare.with_coatings(vec![layer(CoatingSide::Positive, 20e-6)], &cancelled)
            .is_err()
    );
}

#[test]
fn g1_curved_directors_have_a_continuous_vanishing_film_limit() {
    let gate = CancelGate::new_clock_free();
    let (flat, _, constraints) = coupon(&gate);
    let mut normals = Vec::new();
    let nodes = flat
        .mesh()
        .nodes
        .iter()
        .map(|&[x, y, _]| {
            let (u, v, rise) = (x / 0.3, y / 0.2, 0.012);
            let z = 16. * rise * u * (1. - u) * v * (1. - v);
            let nx = -16. * rise * (1. - 2. * u) * v * (1. - v) / 0.3;
            let ny = -16. * rise * u * (1. - u) * (1. - 2. * v) / 0.2;
            let norm = (1. + nx * nx + ny * ny).sqrt();
            normals.push([nx / norm, ny / norm, 1. / norm]);
            [x, y, z]
        })
        .collect();
    let mesh = ShellMesh::new(nodes, flat.mesh().tris.clone()).unwrap();
    let directors = NodalDirectors::new(&mesh, normals, 0.8).unwrap();
    let bare = compile_material_shear_shell(
        mesh,
        vec![flat.regions()[0].input().clone()],
        budget(&flat),
        &gate,
    )
    .unwrap();
    let coated = bare
        .with_coatings(vec![layer(CoatingSide::Positive, 1e-12)], &gate)
        .unwrap();
    let assemble = |shell: &ResolvedShearShell| {
        assemble_mitc3p(
            shell.mesh(),
            &directors,
            shell.sections(),
            &constraints,
            ShellMass::Consistent,
            budget(shell),
            1e-3,
            || false,
        )
        .unwrap()
    };
    let before = assemble(&bare);
    let after = assemble(&coated);
    // A tiny film must not switch the substrate's interpolation for theta x d.
    // Both sides must use the same varying-director material velocity field.
    for (a, b) in [
        (before.pencil().m.to_dense(), after.pencil().m.to_dense()),
        (before.pencil().k.to_dense(), after.pencil().k.to_dense()),
    ] {
        let scale = a.iter().map(|v| v.abs()).fold(0., f64::max);
        let difference = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0., f64::max);
        assert!(
            difference < 1e-8 * scale,
            "vanishing film changes operator by {}",
            difference / scale
        );
    }
}
