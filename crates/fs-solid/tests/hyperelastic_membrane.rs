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
use fs_material::state_point::{MaterialPropertySelection, resolve_mean_thermal_expansion_strain};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solid::membrane::vibration::VibrationSettings;
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
    source_with_constraint(fibres, modulus_scale, mass, false)
}

fn source_with_constraint(
    fibres: bool,
    modulus_scale: f64,
    mass: f64,
    incompressible: bool,
) -> (MaterialCard, ResolvedSheet) {
    let provenance = Provenance {
        source: "synthetic hyperelastic membrane integration fixture; not measured data".into(),
        license: "test-only".into(),
        artifact: None,
    };
    // "T" is typed in kelvin. A thermal reference binds the sheet's mass, thickness and model
    // at a typed temperature context, and a typed query matches only typed validity.
    let kelvin = QuantitySpec::dimensional(Temperature::DIMS);
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
        validity: ValidityDomain::unconstrained().with_quantity("T", kelvin, 290.0, 310.0),
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
                validity: ValidityDomain::unconstrained().with_quantity("T", kelvin, 290.0, 310.0),
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
        &QueryPoint::new().with_quantity("T", kelvin, 300.0).unwrap(),
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let controls = PlaneStressSettings {
        relative_tolerance: 2e-13,
        ..PlaneStressSettings::default()
    };
    let sheet = if incompressible {
        sheet.with_incompressible_hyperelastic_membrane(&card, 0, controls)
    } else {
        sheet.with_hyperelastic_membrane(&card, 0, controls)
    }
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
    for (fibres, incompressible) in [(false, false), (true, false), (false, true), (true, true)] {
        let p = patch(
            triangle(),
            source_with_constraint(fibres, 1.0, MASS_KG_M2, incompressible).1,
            0.34,
        );
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

// Independent incompressible NH potential of the same four affine facets.
// det(F)=1 sets the normal thickness stretch to 1 / deformed-area ratio;
// differentiating this scalar energy gives the center force and stiffness.
// No production material, membrane, AD or matrix routine enters this oracle.
fn incompressible_pyramid_oracle(stretch: f64, height: f64) -> (f64, f64, f64) {
    let side: f64 = 0.1;
    let q2 = (height / side).powi(2);
    let s2 = stretch * stretch;
    let b = s2 + 2.0 * q2;
    let z2 = 1.0 / (s2 * b);
    let force = 4.0 * MU * THICKNESS_M * height * (1.0 - z2 / b);
    let energy = side * side * THICKNESS_M * MU * (2.0 * s2 + 2.0 * q2 + z2 - 3.0);
    let stiffness = 4.0 * MU * THICKNESS_M * (1.0 - z2 / b + 8.0 * q2 * z2 / (b * b));
    (force, energy, stiffness)
}

#[test]
fn incompressible_source_drives_native_pressure_equilibrium_and_loaded_vibration() {
    let (card, sheet) = source_with_constraint(false, 1.0, MASS_KG_M2, true);
    let p = patch(diamond(), sheet, 0.0);
    let stretch = 1.08;
    let height = 0.02;
    let (force, energy, stiffness) = incompressible_pyramid_oracle(stretch, height);
    let zero = vec![[0.0; 3]; 5];
    let fixed = fixed_rim(stretch);
    let pressure = -3.0 * force / (0.02 * stretch * stretch);
    let spectrum = p
        .solve_vibration_with_pressure(
            &zero,
            pressure,
            &fixed,
            &zero,
            settings(),
            VibrationSettings {
                max_free_dofs: 1,
                relative_residual_tolerance: 1e-10,
            },
            || false,
        )
        .unwrap();
    let equilibrium = spectrum.equilibrium();
    assert!(equilibrium.newton.converged);
    close(equilibrium.displacement_m[0][2], -height, 2e-10, 2e-10);
    close(equilibrium.strain_energy_j, energy, 2e-10, 2e-10);
    close(equilibrium.pressure_work_j, force * height, 2e-10, 2e-10);
    close(equilibrium.residual_n[0][2], 0.0, 2e-9, 0.0);
    close(
        equilibrium.residual_n.iter().skip(1).map(|r| r[2]).sum(),
        3.0 * force,
        2e-9,
        2e-10,
    );
    for &(dof, value) in &fixed {
        assert_eq!(equilibrium.displacement_m[dof / 3][dof % 3], value);
    }
    assert_eq!(spectrum.free_dofs(), &[2]);
    assert_eq!(spectrum.modes().len(), 1);
    // The center's vertical pressure force is independent of its height for
    // this fixed planar rim. Its follower stiffness is exactly zero here.
    let mass = 0.02 * MASS_KG_M2 / 3.0;
    close(spectrum.modes()[0].lambda, stiffness / mass, 1e-6, 2e-8);
    close(p.nodal_mass_kg()[0], mass, 1e-17, 2e-14);
    close(p.reference_mass_kg(), 0.02 * MASS_KG_M2, 1e-17, 2e-14);
    assert_eq!(
        p.sheets()[0].material().card_identity(),
        card.content_hash()
    );

    let mut dead_load = zero.clone();
    dead_load[0][2] = -force;
    let dead = p
        .solve(&dead_load, &fixed, &zero, settings(), || false)
        .unwrap();
    close(dead.displacement_m[0][2], -height, 2e-10, 2e-10);
    close(dead.strain_energy_j, energy, 2e-10, 2e-10);
    close(
        dead.residual_n.iter().skip(1).map(|r| r[2]).sum(),
        force,
        2e-9,
        2e-10,
    );
    assert!(matches!(
        p.solve_with_pressure(&zero, pressure, &fixed, &zero, settings(), || true),
        Err(MembraneError::Cancelled)
    ));
}

#[test]
fn natural_reference_preserves_anisotropic_objectivity_and_assembled_tangent() {
    let sheet = source_with_constraint(true, 1.0, MASS_KG_M2, true)
        .1
        .with_hyperelastic_reference_stretches([1.04, 0.98, 1.02])
        .unwrap();
    let p = patch(triangle(), sheet, 0.34);
    let u = deformation(&p);
    let direction = vec![[0.4, -0.1, 0.7], [-0.3, 0.5, -0.1], [0.2, 0.1, -0.4]];
    let tangent = p.tangent_apply(&u, &direction, || false).unwrap();
    let h = 1e-7;
    let offset = |sign: f64| {
        u.iter()
            .zip(&direction)
            .map(|(u, v)| core::array::from_fn(|i| u[i] + sign * h * v[i]))
            .collect::<Vec<_>>()
    };
    let low = p.evaluate(&offset(-1.0), || false).unwrap();
    let high = p.evaluate(&offset(1.0), || false).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            close(
                (high.internal_gradient_n[i][j] - low.internal_gradient_n[i][j]) / (2.0 * h),
                tangent[i][j],
                3e-5,
                3e-6,
            );
        }
    }
    let (s, c) = 0.73_f64.sin_cos();
    let rotate = |v: [f64; 3]| [c * v[0] - s * v[2], v[1], s * v[0] + c * v[2]];
    let moved = p
        .mesh()
        .nodes
        .iter()
        .zip(&u)
        .map(|(x, u)| {
            let y = rotate([x[0] + u[0], x[1] + u[1], u[2]]);
            [y[0] - x[0] + 0.2, y[1] - x[1] - 0.1, y[2] + 0.07]
        })
        .collect::<Vec<_>>();
    let base = p.evaluate(&u, || false).unwrap();
    let changed = p.evaluate(&moved, || false).unwrap();
    close(changed.strain_energy_j, base.strain_energy_j, 2e-12, 5e-13);
    let rotated_direction = direction.iter().copied().map(rotate).collect::<Vec<_>>();
    let rotated_tangent = p
        .tangent_apply(&moved, &rotated_direction, || false)
        .unwrap();
    for i in 0..3 {
        for (a, b) in rotate(base.internal_gradient_n[i])
            .into_iter()
            .zip(changed.internal_gradient_n[i])
        {
            close(a, b, 3e-11, 5e-12);
        }
        for (a, b) in rotate(tangent[i]).into_iter().zip(rotated_tangent[i]) {
            close(a, b, 3e-9, 5e-12);
        }
    }
}

#[test]
fn mean_thermal_reference_changes_native_pressure_equilibrium_and_modes_at_fixed_mass() {
    let (base, _) = source_with_constraint(false, 1.0, MASS_KG_M2, true);
    let mut claims = base.claims().clone();
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let mean_key = PropertyKey::new(
        "mean_engineering_length_expansion",
        Dims([0, 0, 0, -1, 0, 0]),
    );
    let alpha = 0.003;
    claims
        .insert_claim(PropertyClaim {
            key: mean_key.clone(),
            value: PropertyValue::Scalar {
                value: alpha,
                dims: mean_key.dims(),
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            validity: ValidityDomain::unconstrained()
                .with_quantity("T_ref", temperature, 300.0, 310.0)
                .with_quantity("T_hot", temperature, 300.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: base.models()[0].provenance.clone(),
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
    let original = resolve_sheet(
        &card,
        &QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap(),
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_incompressible_hyperelastic_membrane(
        &card,
        0,
        PlaneStressSettings {
            relative_tolerance: 2e-13,
            ..PlaneStressSettings::default()
        },
    )
    .unwrap();
    let mean = resolve_mean_thermal_expansion_strain(
        &card,
        &mean_key,
        &QueryPoint::new()
            .with_quantity("T_ref", temperature, 300.0)
            .unwrap()
            .with_quantity("T_hot", temperature, 310.0)
            .unwrap(),
        "T_ref",
        "T_hot",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let selected = original
        .clone()
        .with_hyperelastic_mean_thermal_expansion(
            &card,
            [mean.clone(), mean.clone(), mean.clone()],
            "T_ref",
            "T_hot",
        )
        .unwrap();
    let p = patch(diamond(), selected.clone(), 0.0);
    let cold = patch(diamond(), original, 0.0);
    let g = 1.0 + alpha * 10.0;
    let stretch = 1.08;
    let height = 0.02;
    // W_ref(h)=g^3 W_native(h/g). Differentiate with respect to physical
    // center height; the original pressure geometry and reference mass stay.
    let (f, w, k) = incompressible_pyramid_oracle(stretch / g, height / g);
    let (force, energy, stiffness) = (g * g * f, g * g * g * w, g * k);
    let pressure = -3.0 * force / (0.02 * stretch * stretch);
    let zero = vec![[0.0; 3]; 5];
    let fixed = fixed_rim(stretch);
    let spectrum = p
        .solve_vibration_with_pressure(
            &zero,
            pressure,
            &fixed,
            &zero,
            settings(),
            VibrationSettings {
                max_free_dofs: 1,
                relative_residual_tolerance: 1e-10,
            },
            || false,
        )
        .unwrap();
    let equilibrium = spectrum.equilibrium();
    close(equilibrium.displacement_m[0][2], -height, 2e-10, 2e-10);
    close(equilibrium.strain_energy_j, energy, 3e-10, 3e-10);
    close(equilibrium.pressure_work_j, force * height, 3e-10, 3e-10);
    close(equilibrium.residual_n[0][2], 0.0, 2e-9, 0.0);
    let mass = 0.02 * MASS_KG_M2 / 3.0;
    close(spectrum.modes()[0].lambda, stiffness / mass, 1e-6, 2e-8);
    assert_eq!(p.nodal_mass_kg(), cold.nodal_mass_kg());
    assert_eq!(p.reference_mass_kg(), cold.reference_mass_kg());
    assert_eq!(&p.sheets()[0], &selected);
    let fs_material::sheet::MembraneLaw::Hyperelastic(bound) =
        p.sheets()[0].membrane_law().unwrap()
    else {
        panic!("wrong law")
    };
    assert_eq!(bound.mean_thermal_expansion().unwrap()[0], mean);
    assert_eq!(bound.material().card_identity(), card.content_hash());
    assert_eq!(bound.model().context().axes()["T"], 310.0);
    assert_eq!(bound.material().query().axes()["T"], 300.0);
    let unchanged_reference = cold
        .solve_with_pressure(&zero, pressure, &fixed, &zero, settings(), || false)
        .unwrap();
    assert!(
        unchanged_reference.displacement_m[0][2].abs() < 0.95 * height,
        "finite thermal expansion must change the actual prestressed equilibrium"
    );
}
