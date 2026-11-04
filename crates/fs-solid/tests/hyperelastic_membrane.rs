//! Native hyperelastic source cards consumed by the existing 3-D membrane owner.
//! Synthetic verification sources; no fitted or empirical material claim.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::sheet::{
    AREAL_DENSITY, ResolvedSheet, SheetCapabilities, THICKNESS, resolve_sheet,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Dims;
use fs_solid::membrane::{
    MembraneBudget, MembraneDomain, MembraneError, MembranePatch, MembraneSolveSettings,
};
use fs_solid::mesh2::Mesh2;
use fs_solver::nonlinear::NewtonKrylovConfig;

const MU: f64 = 3.0e5;
const LAMBDA: f64 = 8.0e5;
const THICKNESS_M: f64 = 0.0013;
const MASS_KG_M2: f64 = 0.24;

fn source(fibres: bool, modulus_scale: f64, mass: f64) -> (MaterialCard, ResolvedSheet) {
    let provenance = Provenance {
        source: "synthetic hyperelastic membrane integration fixture; not measured data".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let mut parameters = vec![("stretch_limit", 5.0, Dims::NONE)];
    if fibres {
        parameters.extend([
            ("c10", 0.7e5 * modulus_scale, pa),
            ("k1", 1.8e5 * modulus_scale, pa),
            ("k2", 3.2, Dims::NONE),
            ("dispersion", 0.08, Dims::NONE),
            ("fibre_invariants", 1.0, Dims::NONE),
            ("fibre_count", 1.0, Dims::NONE),
            ("fibre_0_x", 0.8, Dims::NONE),
            ("fibre_0_y", 0.0, Dims::NONE),
            ("fibre_0_z", 0.6, Dims::NONE),
            ("volumetric_kind", 0.0, Dims::NONE),
            ("kappa", 1.2e6 * modulus_scale, pa),
        ]);
    } else {
        parameters.extend([
            ("mu", MU * modulus_scale, pa),
            ("lambda", LAMBDA * modulus_scale, pa),
        ]);
    }
    let model = ConstitutiveModelCard {
        law: LawId(
            if fibres {
                "hgo-dispersed"
            } else {
                "neo-hookean-log"
            }
            .into(),
        ),
        law_version: 1,
        parameters: parameters
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![],
        provenance: provenance.clone(),
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (AREAL_DENSITY, mass, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, THICKNESS_M, Dims([1, 0, 0, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hyperelastic membrane".into(),
            phase: "solid".into(),
            process: "native consumer verification".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap();
    let sheet = resolve_sheet(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_hyperelastic_membrane(
        &card,
        0,
        PlaneStressSettings {
            relative_tolerance: 2e-13,
            ..PlaneStressSettings::default()
        },
    )
    .unwrap();
    (card, sheet)
}

fn triangle() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    }
}

fn diamond() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [-0.1, 0.0], [0.0, -0.1]],
        elems: vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 4], vec![0, 4, 1]],
        patches: vec![],
    }
}

fn patch(mesh: Mesh2, sheet: ResolvedSheet, grain: f64) -> MembranePatch {
    let n = mesh.elems.len();
    MembranePatch::try_new(
        mesh,
        vec![sheet; n],
        vec![grain; n],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 1.0,
        },
        MembraneBudget {
            max_nodes: 16,
            max_triangles: 16,
        },
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn deformation(p: &MembranePatch) -> Vec<[f64; 3]> {
    p.mesh()
        .nodes
        .iter()
        .map(|x| {
            [
                0.18 * x[0] + 0.12 * x[1],
                0.06 * x[0] + 0.07 * x[1],
                0.1 * x[0] - 0.04 * x[1],
            ]
        })
        .collect()
}

#[test]
fn native_membrane_force_and_geometric_tangent_are_exact_at_finite_deformation() {
    for fibres in [false, true] {
        let p = patch(triangle(), source(fibres, 1.0, MASS_KG_M2).1, 0.34);
        let u = deformation(&p);
        let base = p.evaluate(&u, || false).unwrap();
        let h = 1e-7;
        for dof in 0..9 {
            let mut lo = u.clone();
            let mut hi = u.clone();
            lo[dof / 3][dof % 3] -= h;
            hi[dof / 3][dof % 3] += h;
            let a = p.evaluate(&lo, || false).unwrap();
            let b = p.evaluate(&hi, || false).unwrap();
            close(
                (b.strain_energy_j - a.strain_energy_j) / (2.0 * h),
                base.internal_gradient_n[dof / 3][dof % 3],
                2e-7,
                2e-6,
            );
            let mut v = vec![[0.0; 3]; 3];
            v[dof / 3][dof % 3] = 1.0;
            let tangent = p.tangent_apply(&u, &v, || false).unwrap();
            for node in 0..3 {
                for axis in 0..3 {
                    close(
                        (b.internal_gradient_n[node][axis] - a.internal_gradient_n[node][axis])
                            / (2.0 * h),
                        tangent[node][axis],
                        3e-5,
                        3e-6,
                    );
                }
            }
        }
        assert!(matches!(
            p.evaluate(&u, || true),
            Err(MembraneError::Cancelled)
        ));
    }
}

#[test]
fn anisotropic_sheet_is_objective_in_three_dimensions_and_keeps_source_mass() {
    let (card, original) = source(true, 1.0, MASS_KG_M2);
    let p = patch(triangle(), original.clone(), 0.34);
    let u = deformation(&p);
    let base = p.evaluate(&u, || false).unwrap();
    let (s, c) = 0.73_f64.sin_cos();
    let rotate = |v: [f64; 3]| [c * v[0] - s * v[2], v[1], s * v[0] + c * v[2]];
    let moved: Vec<_> = p
        .mesh()
        .nodes
        .iter()
        .zip(&u)
        .map(|(x, u)| {
            let y = rotate([x[0] + u[0], x[1] + u[1], u[2]]);
            [y[0] - x[0] + 0.2, y[1] - x[1] - 0.1, y[2] + 0.07]
        })
        .collect();
    let changed = p.evaluate(&moved, || false).unwrap();
    close(changed.strain_energy_j, base.strain_energy_j, 1e-12, 3e-13);
    for (a, b) in base
        .internal_gradient_n
        .iter()
        .zip(&changed.internal_gradient_n)
    {
        for (expected, actual) in rotate(*a).into_iter().zip(*b) {
            close(actual, expected, 2e-11, 3e-12);
        }
    }
    let direction = vec![[0.4, -0.1, 0.7], [-0.3, 0.5, -0.1], [0.2, 0.1, -0.4]];
    let rotated_direction: Vec<_> = direction.iter().copied().map(rotate).collect();
    let tangent = p.tangent_apply(&u, &direction, || false).unwrap();
    let rotated_tangent = p
        .tangent_apply(&moved, &rotated_direction, || false)
        .unwrap();
    for (a, b) in tangent.iter().zip(&rotated_tangent) {
        for (expected, actual) in rotate(*a).into_iter().zip(*b) {
            close(actual, expected, 2e-9, 3e-12);
        }
    }
    let rigid = p
        .mesh()
        .nodes
        .iter()
        .map(|x| {
            let y = rotate([x[0], x[1], 0.0]);
            [y[0] - x[0] + 0.2, y[1] - x[1] - 0.1, y[2] + 0.07]
        })
        .collect::<Vec<_>>();
    close(
        p.evaluate(&rigid, || false).unwrap().strain_energy_j,
        0.0,
        1e-14,
        0.0,
    );

    let heavy = patch(triangle(), source(true, 1.0, 3.0 * MASS_KG_M2).1, 0.34);
    assert_eq!(
        base.strain_energy_j,
        heavy.evaluate(&u, || false).unwrap().strain_energy_j
    );
    close(p.reference_mass_kg(), 0.005 * MASS_KG_M2, 1e-17, 1e-14);
    close(
        heavy.reference_mass_kg(),
        3.0 * p.reference_mass_kg(),
        1e-17,
        1e-14,
    );
    let body = heavy.body_load_n([0.0, 0.0, -9.81]).unwrap();
    close(
        body.iter().map(|f| f[2]).sum(),
        -9.81 * heavy.reference_mass_kg(),
        1e-15,
        1e-14,
    );
    assert_eq!(&p.sheets()[0], &original);
    assert_eq!(
        p.sheets()[0].material().card_identity(),
        card.content_hash()
    );
    assert!(p.sheets()[0].bending().is_err());
    let different_grain = patch(triangle(), original, 1.1);
    let changed_energy = different_grain
        .evaluate(&u, || false)
        .unwrap()
        .strain_energy_j;
    assert!((changed_energy - base.strain_energy_j).abs() > 0.01 * base.strain_energy_j);
}

fn fixed_rim(stretch: f64) -> Vec<(usize, f64)> {
    let mut fixed = vec![(0, 0.0), (1, 0.0)];
    for (node, x) in diamond().nodes.into_iter().enumerate().skip(1) {
        fixed.extend([
            (3 * node, (stretch - 1.0) * x[0]),
            (3 * node + 1, (stretch - 1.0) * x[1]),
            (3 * node + 2, 0.0),
        ]);
    }
    fixed
}

fn settings() -> MembraneSolveSettings {
    MembraneSolveSettings {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-9,
            relative_tolerance: 2e-11,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
        max_iterations: 30,
        max_evaluations: 1000,
        max_krylov_restart: 8,
    }
}

// Independent closed-form C and scalar transverse stationarity for the four
// congruent pyramid facets. This never calls either membrane or material code.
fn pyramid_oracle(stretch: f64, height: f64) -> (f64, f64) {
    let side = 0.1;
    let q2 = (height / side).powi(2);
    let area_ratio = stretch * (stretch * stretch + 2.0 * q2).sqrt();
    let (mut lo, mut hi): (f64, f64) = (1e-8, 10.0);
    for _ in 0..100 {
        let z = 0.5 * (lo + hi);
        if MU * (z * z - 1.0) + LAMBDA * (area_ratio * z).ln() > 0.0 {
            hi = z;
        } else {
            lo = z;
        }
    }
    let z = 0.5 * (lo + hi);
    let log_j = (area_ratio * z).ln();
    let force = 4.0 * MU * THICKNESS_M * height * (1.0 - z * z / (stretch * stretch + 2.0 * q2));
    let energy = 2.0
        * side
        * side
        * THICKNESS_M
        * (0.5 * MU * (2.0 * stretch * stretch + 2.0 * q2 + z * z - 3.0) - MU * log_j
            + 0.5 * LAMBDA * log_j * log_j);
    (force, energy)
}

#[test]
fn finite_strain_dead_load_and_pressure_equilibria_match_independent_pyramid() {
    let p = patch(diamond(), source(false, 1.0, MASS_KG_M2).1, 0.0);
    let stretch = 1.08;
    let height = 0.02;
    let (force, energy) = pyramid_oracle(stretch, height);
    let mut loads = vec![[0.0; 3]; 5];
    loads[0][2] = -force;
    let initial = vec![[0.0; 3]; 5];
    let fixed = fixed_rim(stretch);
    let solution = p
        .solve(&loads, &fixed, &initial, settings(), || false)
        .unwrap();
    assert!(solution.newton.converged);
    close(solution.displacement_m[0][2], -height, 2e-10, 2e-10);
    close(solution.strain_energy_j, energy, 2e-10, 2e-10);
    close(solution.residual_n[0][2], 0.0, 2e-9, 0.0);
    close(
        solution.residual_n.iter().skip(1).map(|r| r[2]).sum(),
        force,
        2e-9,
        2e-10,
    );
    for &(dof, value) in &fixed {
        assert_eq!(solution.displacement_m[dof / 3][dof % 3], value);
    }

    // The original current-area pressure owner assigns one third of the
    // total oriented pressure force to the center. Rim pressure stays in its
    // reactions; the exact reference-area hyperelastic energy is unchanged.
    let pressure = -3.0 * force / (0.02 * stretch * stretch);
    let pressurized = p
        .solve_with_pressure(&initial, pressure, &fixed, &initial, settings(), || false)
        .unwrap();
    assert!(pressurized.newton.converged);
    close(pressurized.displacement_m[0][2], -height, 2e-10, 2e-10);
    close(pressurized.strain_energy_j, energy, 2e-10, 2e-10);
    close(pressurized.pressure_work_j, force * height, 2e-10, 2e-10);
    close(
        pressurized.residual_n.iter().skip(1).map(|r| r[2]).sum(),
        3.0 * force,
        2e-9,
        2e-10,
    );
    close(p.reference_mass_kg(), 0.02 * MASS_KG_M2, 1e-16, 1e-14);

    let stiffer = patch(diamond(), source(false, 2.0, MASS_KG_M2).1, 0.0);
    let changed = stiffer
        .solve(&loads, &fixed, &initial, settings(), || false)
        .unwrap();
    assert!(changed.displacement_m[0][2].abs() < height);
    assert_eq!(stiffer.reference_mass_kg(), p.reference_mass_kg());
    assert!(matches!(
        p.solve(&loads, &fixed, &initial, settings(), || true),
        Err(MembraneError::Cancelled)
    ));
}
