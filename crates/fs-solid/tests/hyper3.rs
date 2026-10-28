//! Three-dimensional finite-strain hyperelasticity on P1 tetrahedra.
//!
//! G0: the assembled tangent is the exact derivative of the internal force
//! (central-difference directional check); at the reference state it equals
//! the independently written small-strain `linear3` stiffness.
//! G1: the affine patch test is reproduced exactly at finite strain
//! (interior nodes and stored energy).
//! G2: compressible Neo-Hookean uniaxial stretch to lambda = 1.6 against the
//! closed-form lateral stretch and nominal stress.
//! G3: a rigid rotation of the boundary stores no energy and needs no force.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::linear3::{TetLinearElasticProblem, TetMaterialField};
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget,
    TetElasticMaterial,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Conforming box mesh: every cube split into six tetrahedra around its
/// main diagonal (Kuhn), positively oriented.
fn box_mesh(n: [usize; 3], size: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let id = |i: usize, j: usize, k: usize| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let mut nodes = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                nodes.push([
                    size[0] * i as f64 / n[0] as f64,
                    size[1] * j as f64 / n[1] as f64,
                    size[2] * k as f64 / n[2] as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                for perm in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let mut c = [i, j, k];
                    let mut tet = [id(c[0], c[1], c[2]); 4];
                    for (step, axis) in perm.iter().enumerate() {
                        c[*axis] += 1;
                        tet[step + 1] = id(c[0], c[1], c[2]);
                    }
                    let x = tet.map(|v| nodes[v]);
                    let e = |a: usize| [0, 1, 2].map(|d| x[a][d] - x[0][d]);
                    let (a, b, cc) = (e(1), e(2), e(3));
                    let det = (a[1] * b[2] - a[2] * b[1]) * cc[0]
                        + (a[2] * b[0] - a[0] * b[2]) * cc[1]
                        + (a[0] * b[1] - a[1] * b[0]) * cc[2];
                    if det < 0.0 {
                        tet.swap(2, 3);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (nodes, tets)
}

fn neo_hookean(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 5.0).unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a Hyperelastic,
    prescribed: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
    load_steps: usize,
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: prescribed,
        nodal_forces_n: forces,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings {
            load_steps,
            ..HyperTetSettings::default()
        },
    }
}

fn on_boundary(x: [f64; 3], size: [f64; 3]) -> bool {
    (0..3).any(|d| x[d].abs() < 1e-12 || (x[d] - size[d]).abs() < 1e-12)
}

#[test]
fn affine_patch_test_is_exact_at_finite_strain() {
    let size = [1.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([3, 3, 3], size);
    let material = neo_hookean(1.0, 2.0);
    let f0 = [1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.02, 0.0, 1.1];
    let affine = |x: [f64; 3], i: usize| {
        (0..3)
            .map(|j| (f0[3 * i + j] - f64::from(u8::from(i == j))) * x[j])
            .sum::<f64>()
    };
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if on_boundary(x, size) {
            for i in 0..3 {
                prescribed.push((3 * node + i, affine(x, i)));
            }
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 3);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let mut worst = 0.0f64;
    for (node, &x) in nodes.iter().enumerate() {
        for i in 0..3 {
            worst = worst.max((solution.displacement_m[node][i] - affine(x, i)).abs());
        }
    }
    let exact_energy = material.energy(&f0);
    println!(
        "patch: worst interior error {worst:e}, energy {} vs {exact_energy}",
        solution.strain_energy_j
    );
    assert!(worst < 1e-10, "affine field not reproduced: {worst}");
    assert!((solution.strain_energy_j - exact_energy).abs() < 1e-12 * exact_energy);
    // Global equilibrium: the support forces balance.
    for i in 0..3 {
        let total: f64 = solution
            .reactions_n
            .iter()
            .filter(|(dof, _)| dof % 3 == i)
            .map(|(_, r)| r)
            .sum();
        assert!(
            total.abs() < 1e-10,
            "component {i} reaction resultant {total}"
        );
    }
}

/// Lateral stretch `t` of compressible Neo-Hookean uniaxial stress:
/// `mu (t^2 - 1) + lambda ln(stretch t^2) = 0`.
fn lateral_stretch(mu: f64, lambda: f64, stretch: f64) -> f64 {
    let mut t = 1.0f64;
    for _ in 0..60 {
        let g = mu * (t * t - 1.0) + lambda * (stretch * t * t).ln();
        let dg = 2.0 * mu * t + 2.0 * lambda / t;
        t -= g / dg;
    }
    t
}

#[test]
fn neo_hookean_uniaxial_stretch_matches_closed_form() {
    let size = [2.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([4, 2, 2], size);
    let (mu, lambda) = (1.0, 1.5);
    let material = neo_hookean(mu, lambda);
    let stretch = 1.6;
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if x[0].abs() < 1e-12 {
            prescribed.push((3 * node, 0.0));
        }
        if (x[0] - size[0]).abs() < 1e-12 {
            prescribed.push((3 * node, (stretch - 1.0) * size[0]));
        }
        if x[1].abs() < 1e-12 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2].abs() < 1e-12 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 4);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let t = lateral_stretch(mu, lambda, stretch);
    let j = stretch * t * t;
    let nominal = mu * (stretch - 1.0 / stretch) + lambda * j.ln() / stretch;
    let axial: f64 = prescribed
        .iter()
        .zip(&solution.reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
        .map(|(_, (_, r))| r)
        .sum();
    let corner = nodes
        .iter()
        .position(|x| x.iter().zip(&size).all(|(a, b)| (a - b).abs() < 1e-12))
        .unwrap();
    println!(
        "uniaxial lambda={stretch}: lateral {:.12} vs {t:.12}, force {axial:.12} vs {nominal:.12}, steps {:?}",
        1.0 + solution.displacement_m[corner][1],
        solution
            .steps
            .iter()
            .map(|s| s.residual_history.len())
            .collect::<Vec<_>>()
    );
    assert!((1.0 + solution.displacement_m[corner][1] - t).abs() < 1e-9);
    assert!((1.0 + solution.displacement_m[corner][2] - t).abs() < 1e-9);
    assert!(
        (axial - nominal).abs() < 1e-9 * nominal,
        "{axial} vs {nominal}"
    );
    // Newton's local quadratic rate: few iterations per load step.
    for step in &solution.steps {
        assert!(step.residual_history.len() <= 8, "{step:?}");
        assert_eq!(step.indefinite_fallbacks, 0);
    }
}

#[test]
fn ogden_card_solves_uniaxial_bar() {
    use fs_material::hyper::{OgdenTerm, Volumetric};
    // Authored constitutive fixture in Pa, retaining the published exponent
    // signs. Moderate bulk/shear ratio avoids making a locking claim for P1.
    let terms = [
        OgdenTerm {
            mu_pa: 0.63,
            alpha: 1.3,
        },
        OgdenTerm {
            mu_pa: 0.0012,
            alpha: 5.,
        },
        OgdenTerm {
            mu_pa: -0.01,
            alpha: -2.,
        },
    ];
    let kappa = 20. * 0.4225;
    let material = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: terms.to_vec(),
            volumetric: Volumetric::Quadratic { kappa },
        },
        5.,
    )
    .unwrap();
    let stretch = 1.4;
    // Independent traction-free compressible material-point root. These
    // scalar principal stresses do not call the production AD/eigen path.
    let stresses = |j: f64| {
        let s = [stretch, (j / stretch).sqrt(), (j / stretch).sqrt()];
        let tau = s.map(|l| {
            terms
                .iter()
                .map(|p| {
                    p.mu_pa
                        * fs_math::det::exp(
                            p.alpha
                                * fs_math::det::ln(
                                    l * fs_math::det::exp(-fs_math::det::ln(j) / 3.),
                                ),
                        )
                })
                .sum::<f64>()
        });
        let mean = tau.iter().sum::<f64>() / 3.;
        core::array::from_fn::<_, 3, _>(|i| (tau[i] - mean + kappa * j * (j - 1.)) / s[i])
    };
    let (mut lo, mut hi) = (1., 1.2);
    assert!(stresses(lo)[1] < 0. && stresses(hi)[1] > 0.);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if stresses(mid)[1] > 0. {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let j = 0.5 * (lo + hi);
    let nominal = stresses(j)[0];
    let lateral = (j / stretch).sqrt();
    // The exact solution is affine, so P1 represents it exactly at every
    // rung. The measured band covers nonlinear/linear solver error only.
    for n in [[1, 1, 1], [2, 1, 1], [4, 2, 2]] {
        let size = [2., 1., 1.];
        let (nodes, tets) = box_mesh(n, size);
        let mut prescribed = Vec::new();
        for (node, &x) in nodes.iter().enumerate() {
            if x[0].abs() < 1e-12 {
                prescribed.push((3 * node, 0.));
            }
            if (x[0] - size[0]).abs() < 1e-12 {
                prescribed.push((3 * node, (stretch - 1.) * size[0]));
            }
            if x[1].abs() < 1e-12 {
                prescribed.push((3 * node + 1, 0.));
            }
            if x[2].abs() < 1e-12 {
                prescribed.push((3 * node + 2, 0.));
            }
        }
        let solution =
            with_cx(|cx| problem(&nodes, &tets, &material, &prescribed, &[], 4).solve(cx)).unwrap();
        let axial: f64 = solution
            .reactions_n
            .iter()
            .filter(|(dof, _)| *dof % 3 == 0 && (nodes[*dof / 3][0] - size[0]).abs() < 1e-12)
            .map(|(_, force)| *force)
            .sum();
        let error = (axial - nominal).abs() / nominal.abs();
        assert!(error < 1e-8, "{n:?}: reaction {axial} vs {nominal}");
        for (x, u) in nodes.iter().zip(&solution.displacement_m) {
            assert!((u[0] - (stretch - 1.) * x[0]).abs() < 1e-8);
            for d in [1, 2] {
                assert!((u[d] - (lateral - 1.) * x[d]).abs() < 1e-8);
            }
        }
        println!(
            "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"ogden_card_solves_uniaxial_bar\",\"tier\":\"G1\",\"mesh\":{n:?},\"observed_relative_reaction_error\":{error},\"tolerance\":1e-8,\"locking_claim\":false,\"verdict\":\"pass\"}}"
        );
    }
}

#[test]
fn tangent_is_the_derivative_of_the_internal_force() {
    let (nodes, tets) = box_mesh([2, 2, 2], [1.0, 1.0, 1.0]);
    for material in [
        neo_hookean(1.0, 3.0),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 0.4,
                c01: 0.1,
                kappa: 5.0,
            },
            5.0,
        )
        .unwrap(),
    ] {
        let p = problem(&nodes, &tets, &material, &[], &[], 1);
        // A smooth, finite deformation and a pseudo-random direction.
        let u: Vec<f64> = nodes
            .iter()
            .flat_map(|x| {
                [
                    0.2 * x[1] * x[2],
                    0.1 * x[0] * x[0] - 0.05 * x[2],
                    0.15 * x[0] * x[1],
                ]
            })
            .collect();
        let d: Vec<f64> = (0..u.len())
            .map(|i| ((i * 7919 % 101) as f64 / 101.0) - 0.5)
            .collect();
        let (k, plus, minus) = with_cx(|cx| {
            let eps = 1e-6;
            let up: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a + eps * b).collect();
            let um: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a - eps * b).collect();
            (
                p.tangent_matrix(&u, cx).unwrap(),
                p.internal_force(&up, cx).unwrap(),
                p.internal_force(&um, cx).unwrap(),
            )
        });
        let mut kd = vec![0.0; u.len()];
        k.spmv(&d, &mut kd);
        let fd: Vec<f64> = plus
            .iter()
            .zip(&minus)
            .map(|(a, b)| (a - b) / 2e-6)
            .collect();
        let err = kd
            .iter()
            .zip(&fd)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        let scale = kd.iter().map(|a| a * a).sum::<f64>().sqrt();
        println!(
            "{:?}: |K d - FD| / |K d| = {:e}",
            material.model,
            err / scale
        );
        assert!(err / scale < 1e-6);
    }
}

#[test]
fn reference_tangent_equals_the_small_strain_stiffness() {
    let (nodes, tets) = box_mesh([2, 1, 2], [1.0, 0.5, 1.0]);
    let (mu, lambda) = (0.8, 1.7);
    let material = neo_hookean(mu, lambda);
    let young = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
    let poisson = lambda / (2.0 * (lambda + mu));
    let linear = TetElasticMaterial::try_new(1.0, young, poisson, [7; 32]).unwrap();
    let (k_hyper, k_linear) = with_cx(|cx| {
        let hyper = problem(&nodes, &tets, &material, &[], &[], 1);
        let reference = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&linear),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        };
        (
            hyper
                .tangent_matrix(&vec![0.0; 3 * nodes.len()], cx)
                .unwrap(),
            reference.assemble(cx).unwrap().stiffness,
        )
    });
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for r in 0..k_linear.nrows() {
        for c in 0..k_linear.ncols() {
            worst = worst.max((k_hyper.get(r, c) - k_linear.get(r, c)).abs());
            scale = scale.max(k_linear.get(r, c).abs());
        }
    }
    println!("reference tangent vs linear3: max |dK| = {worst:e} (scale {scale:e})");
    assert!(worst < 1e-12 * scale);
}

#[test]
fn rigid_rotation_stores_no_energy_and_needs_no_force() {
    let size = [1.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([2, 2, 2], size);
    let material = neo_hookean(1.0, 1.0);
    let (c, s) = (0.6f64, 0.8f64); // 53 degrees about z
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if on_boundary(x, size) {
            for i in 0..3 {
                let rx: f64 = (0..3).map(|j| rotation[3 * i + j] * x[j]).sum();
                prescribed.push((3 * node + i, rx - x[i]));
            }
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 4);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let largest = solution
        .reactions_n
        .iter()
        .map(|(_, r)| r.abs())
        .fold(0.0, f64::max);
    println!(
        "rotation: energy {:e}, largest reaction {largest:e}",
        solution.strain_energy_j
    );
    assert!(solution.strain_energy_j.abs() < 1e-12);
    assert!(largest < 1e-10);
    assert!((solution.min_det_f - 1.0).abs() < 1e-12);
}

#[test]
fn refusals_are_structured() {
    let (nodes, tets) = box_mesh([1, 1, 1], [1.0, 1.0, 1.0]);
    let material = neo_hookean(1.0, 1.0);
    // Inverted element.
    let mut flipped = tets.clone();
    flipped[0].swap(0, 1);
    let err =
        with_cx(|cx| problem(&nodes, &flipped, &material, &[], &[], 1).solve(cx)).unwrap_err();
    assert!(
        matches!(err, HyperTetError::DegenerateElement { element: 0, .. }),
        "{err}"
    );
    // Repeated prescribed DOF.
    let twice = [(0, 0.0), (0, 0.1)];
    let err =
        with_cx(|cx| problem(&nodes, &tets, &material, &twice, &[], 1).solve(cx)).unwrap_err();
    assert!(matches!(err, HyperTetError::InvalidInput { .. }), "{err}");
    // A tripped context publishes nothing.
    let gate = CancelGate::new();
    gate.request();
    let err = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0)];
        problem(&nodes, &tets, &material, &fixed, &[(21, 0.1)], 1).solve(&cx)
    });
    assert_eq!(err.unwrap_err(), HyperTetError::Cancelled);
}

fn layered_source_card(scale: f64) -> fs_matdb::ConstitutiveModelCard {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use fs_qty::Dims;
    let pressure = Dims([-1, 1, -2, 0, 0, 0]);
    let artifact = fs_blake3::hash_bytes(
        format!("synthetic Ogden alpha=2 layer fixture; modulus scale={scale}").as_bytes(),
    );
    ConstitutiveModelCard {
        law: LawId("ogden-decoupled".into()),
        law_version: 1,
        parameters: [
            ("stretch_limit", 5.0, Dims::NONE),
            ("volumetric_kind", 0.0, Dims::NONE),
            ("term_count", 1.0, Dims::NONE),
            ("mu_0", 3.0 * scale, pressure),
            ("alpha_0", 2.0, Dims::NONE),
            ("kappa", 8.0 * scale, pressure),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: fs_evidence::ValidityDomain::unconstrained(),
        sources: vec![artifact],
        provenance: Provenance {
            source: "manufactured two-layer finite-strain material; not measured data".into(),
            license: "test-only".into(),
            artifact: Some(artifact),
        },
    }
}

fn layered_constraints(nodes: &[[f64; 3]], extension: f64) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        // The transverse rollers impose F22=F33=1. The internal x=1 face
        // remains free axially, so its displacement must come from equilibrium.
        prescribed.extend([(3 * node + 1, 0.0), (3 * node + 2, 0.0)]);
        if x[0] == 0.0 {
            prescribed.push((3 * node, 0.0));
        } else if x[0] == 2.0 {
            prescribed.push((3 * node, extension));
        }
    }
    prescribed
}

#[test]
fn heterogeneous_source_cards_drive_the_independent_layer_equilibrium() {
    use fs_material::hyper::cards::ResolvedHyperelastic;
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let extension = 0.4;
    let prescribed = layered_constraints(&nodes, extension);
    // For the alpha=2 Ogden source and F=diag(lambda,1,1), differentiating
    // W=mu/2*(lambda^(4/3)+2 lambda^(-2/3)-3)+kappa/2*(lambda-1)^2 gives
    // P11=2 mu/3*(lambda^(1/3)-lambda^(-5/3))+kappa*(lambda-1).
    // Neither reference calls the constitutive implementation or its AD path.
    let stress = |stretch: f64, scale: f64| {
        scale * (2.0 * (stretch.powf(1.0 / 3.0) - stretch.powf(-5.0 / 3.0)) + 8.0 * (stretch - 1.0))
    };
    let energy = |stretch: f64, scale: f64| {
        scale
            * (1.5 * (stretch.powf(4.0 / 3.0) + 2.0 * stretch.powf(-2.0 / 3.0) - 3.0)
                + 4.0 * (stretch - 1.0).powi(2))
    };
    let mut prior_interface = 0.0;
    let mut prior_force = 0.0;
    for right_scale in [2.0, 5.0] {
        let cards = [layered_source_card(1.0), layered_source_card(right_scale)];
        let sources: Vec<_> = tets
            .iter()
            .map(|tet| {
                let x = tet.iter().map(|&n| nodes[n][0]).sum::<f64>() / 4.0;
                ResolvedHyperelastic::try_from_model_card(&cards[usize::from(x > 1.0)]).unwrap()
            })
            .collect();
        let p = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::PerElement(&sources),
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 4,
                ..HyperTetSettings::default()
            },
        };
        let solved = with_cx(|cx| p.solve(cx)).unwrap();
        // Two unit-length, unit-area layers in series: lambda1+lambda2=2+du
        // and P1(lambda1)=P2(lambda2). Both stress laws are strictly increasing.
        let (mut lo, mut hi) = (0.0, extension);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if stress(1.0 + mid, 1.0) > stress(1.0 + extension - mid, right_scale) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let interface = 0.5 * (lo + hi);
        for (node, x) in nodes.iter().enumerate() {
            if x[0] == 1.0 {
                assert!((solved.displacement_m[node][0] - interface).abs() < 1e-9);
            }
        }
        let force: f64 = solved
            .reactions_n
            .iter()
            .filter(|(dof, _)| dof % 3 == 0 && nodes[dof / 3][0] == 2.0)
            .map(|(_, force)| force)
            .sum();
        assert!((force - stress(1.0 + interface, 1.0)).abs() < 1e-9);
        let exact_energy =
            energy(1.0 + interface, 1.0) + energy(1.0 + extension - interface, right_scale);
        assert!((solved.strain_energy_j - exact_energy).abs() < 1e-9);
        assert!(interface > prior_interface && force > prior_force);
        prior_interface = interface;
        prior_force = force;
        for (tet, source) in tets.iter().zip(&sources) {
            let x = tet.iter().map(|&n| nodes[n][0]).sum::<f64>() / 4.0;
            assert_eq!(source.model_card(), &cards[usize::from(x > 1.0)]);
            assert_eq!(
                source.model_card().content_hash(),
                cards[usize::from(x > 1.0)].content_hash()
            );
        }
    }
}

#[test]
fn source_fields_preserve_uniform_results_and_each_elements_domain() {
    use fs_material::hyper::cards::ResolvedHyperelastic;
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let prescribed = layered_constraints(&nodes, 0.4);
    let card = layered_source_card(1.0);
    let source = ResolvedHyperelastic::try_from_model_card(&card).unwrap();
    let native = problem(&nodes, &tets, source.law(), &prescribed, &[], 4);
    let mut sourced = HyperTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        material: HyperelasticLawField::Uniform(&source),
        prescribed_m: &prescribed,
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: native.budget,
        settings: native.settings,
    };
    let field = vec![source.clone(); tets.len()];
    with_cx(|cx| {
        let baseline = native.solve(cx).unwrap();
        assert_eq!(baseline, sourced.solve(cx).unwrap());
        sourced.material = HyperelasticLawField::PerElement(&field);
        assert_eq!(baseline, sourced.solve(cx).unwrap());
        sourced.material = HyperelasticLawField::PerElement(&field[..field.len() - 1]);
        assert_eq!(
            sourced.solve(cx).unwrap_err(),
            HyperTetError::LawCountMismatch {
                expected: tets.len(),
                got: tets.len() - 1,
            }
        );
    });
    let mut limited = card.clone();
    limited.parameters.get_mut("stretch_limit").unwrap().value = 1.05;
    let mut limited_field = field.clone();
    limited_field[5] = ResolvedHyperelastic::try_from_model_card(&limited).unwrap();
    sourced.material = HyperelasticLawField::PerElement(&limited_field);
    let u: Vec<_> = nodes.iter().flat_map(|x| [0.2 * x[0], 0.0, 0.0]).collect();
    with_cx(|cx| {
        assert!(matches!(
            sourced.internal_force(&u, cx),
            Err(HyperTetError::MaterialRefused { element: 5, .. })
        ));
        assert!(matches!(
            sourced.tangent_matrix(&u, cx),
            Err(HyperTetError::MaterialRefused { element: 5, .. })
        ));
    });
    let mut wrong_units = card;
    wrong_units.parameters.get_mut("mu_0").unwrap().dims = fs_qty::Dims::NONE;
    assert!(ResolvedHyperelastic::try_from_model_card(&wrong_units).is_err());
}

#[test]
fn mixed_isotropic_source_cards_drive_actual_layer_displacement_and_force() {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use fs_material::hyper::cards::ResolvedHyperelastic;
    use fs_qty::Dims;
    let make_card = |name: &str, moduli: &[(&str, f64)]| {
        let source = fs_blake3::hash_bytes(
            format!("synthetic two-layer {name} fixture: {moduli:?}").as_bytes(),
        );
        let mut parameters: std::collections::BTreeMap<String, LawParameter> = moduli
            .iter()
            .map(|&(name, value)| {
                (
                    name.into(),
                    LawParameter {
                        value,
                        dims: Dims([-1, 1, -2, 0, 0, 0]),
                    },
                )
            })
            .collect();
        parameters.insert(
            "stretch_limit".into(),
            LawParameter {
                value: 5.0,
                dims: Dims::NONE,
            },
        );
        ConstitutiveModelCard {
            law: LawId(name.into()),
            law_version: 1,
            parameters,
            state_schema_version: 1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: fs_evidence::ValidityDomain::unconstrained(),
            sources: vec![source],
            provenance: Provenance {
                source: "manufactured isotropic two-layer fixture; not measured data".into(),
                license: "test-only".into(),
                artifact: Some(source),
            },
        }
    };
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let extension = 0.4;
    let prescribed = layered_constraints(&nodes, extension);
    // Independent one-dimensional derivatives at F=diag(s,1,1), with
    // transverse rollers. The left layer is coupled-log Neo-Hookean;
    // the right is decoupled Mooney–Rivlin with quadratic bulk energy.
    let left_stress = |s: f64| 3.0 * (s - 1.0 / s) + 8.0 * s.ln() / s;
    let left_energy = |s: f64| 1.5 * (s * s - 1.0) - 3.0 * s.ln() + 4.0 * s.ln().powi(2);
    let right_stress = |s: f64, scale: f64| {
        scale
            * ((4.0 / 3.0)
                * (1.1 * (s.powf(1.0 / 3.0) - s.powf(-5.0 / 3.0))
                    + 0.4 * (s.powf(-1.0 / 3.0) - s.powf(-7.0 / 3.0)))
                + 8.0 * (s - 1.0))
    };
    let right_energy = |s: f64, scale: f64| {
        scale
            * (1.1 * (s.powf(4.0 / 3.0) + 2.0 * s.powf(-2.0 / 3.0) - 3.0)
                + 0.4 * (2.0 * s.powf(2.0 / 3.0) + s.powf(-4.0 / 3.0) - 3.0)
                + 4.0 * (s - 1.0).powi(2))
    };
    let (mut previous_displacement, mut previous_force) = (0.0, 0.0);
    for right_scale in [2.0, 5.0] {
        let cards = [
            make_card("neo-hookean-log", &[("mu", 3.0), ("lambda", 8.0)]),
            make_card(
                "mooney-rivlin-decoupled",
                &[
                    ("c10", 1.1 * right_scale),
                    ("c01", 0.4 * right_scale),
                    ("kappa", 8.0 * right_scale),
                ],
            ),
        ];
        let sources: Vec<_> = tets
            .iter()
            .map(|tet| {
                let x = tet.iter().map(|&node| nodes[node][0]).sum::<f64>() / 4.0;
                ResolvedHyperelastic::try_from_model_card(&cards[usize::from(x > 1.0)]).unwrap()
            })
            .collect();
        let problem = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::PerElement(&sources),
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 4,
                ..HyperTetSettings::default()
            },
        };
        let solved = with_cx(|cx| problem.solve(cx)).unwrap();
        // Unit-length/unit-area layers have equal nominal traction and
        // s_left+s_right=2+extension. Solve that independent scalar equation.
        let (mut low, mut high) = (0.0, extension);
        for _ in 0..100 {
            let mid = 0.5 * (low + high);
            if left_stress(1.0 + mid) > right_stress(1.0 + extension - mid, right_scale) {
                high = mid;
            } else {
                low = mid;
            }
        }
        let interface = 0.5 * (low + high);
        for (node, x) in nodes.iter().enumerate().filter(|(_, x)| x[0] == 1.0) {
            assert!(
                (solved.displacement_m[node][0] - interface).abs() < 1e-9,
                "interface at {x:?}: {:?}, expected {interface}",
                solved.displacement_m[node]
            );
        }
        let force: f64 = solved
            .reactions_n
            .iter()
            .filter(|(dof, _)| dof % 3 == 0 && nodes[dof / 3][0] == 2.0)
            .map(|(_, force)| force)
            .sum();
        assert!((force - left_stress(1.0 + interface)).abs() < 1e-9);
        let expected_energy =
            left_energy(1.0 + interface) + right_energy(1.0 + extension - interface, right_scale);
        assert!((solved.strain_energy_j - expected_energy).abs() < 1e-9);
        assert!(interface > previous_displacement && force > previous_force);
        previous_displacement = interface;
        previous_force = force;
        for (tet, source) in tets.iter().zip(&sources) {
            let x = tet.iter().map(|&node| nodes[node][0]).sum::<f64>() / 4.0;
            assert_eq!(source.model_card(), &cards[usize::from(x > 1.0)]);
        }
    }
}

#[test]
fn source_context_requires_supported_fixed_conditions_inside_the_exact_domain() {
    use fs_matdb::QueryPoint;
    use fs_material::hyper::cards::ResolvedHyperelastic;
    use fs_qty::QuantitySpec;
    use fs_qty::semantic::{QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm};
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut card = layered_source_card(1.0);
    card.validity = card
        .validity
        .with_quantity("temperature", absolute, 290.0, 310.0);
    let inside = QueryPoint::new()
        .with_quantity("temperature", absolute, 300.0)
        .unwrap();
    let source = ResolvedHyperelastic::try_from_model_card_at(&card, &inside).unwrap();
    assert_eq!(source.context(), &inside);
    assert_eq!(source.model_card(), &card);
    assert!(ResolvedHyperelastic::try_from_model_card(&card).is_err());
    assert!(
        ResolvedHyperelastic::try_from_model_card_at(
            &card,
            &inside.clone().with("temperature", 311.0).unwrap()
        )
        .is_err()
    );
    let difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    let wrong_kind = QueryPoint::new()
        .with_quantity("temperature", difference, 300.0)
        .unwrap();
    assert!(ResolvedHyperelastic::try_from_model_card_at(&card, &wrong_kind).is_err());
    // A caller-supplied strain is not the strain evolving inside hyper3.
    let strain = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Strain {
            basis: StrainBasis::Tensor,
            component: StrainComponent::Normal,
        },
        ValueForm::Static,
    ));
    card.validity = card
        .validity
        .with_quantity("axial_strain", strain, 0.0, 0.1);
    let frozen = inside.with_quantity("axial_strain", strain, 0.05).unwrap();
    assert!(
        card.validity
            .contains_typed(frozen.axes(), frozen.axis_quantities())
    );
    assert!(ResolvedHyperelastic::try_from_model_card_at(&card, &frozen).is_err());
}

/// G3: translation cannot change finite-strain energy, force, or tangent.
/// The dyadic motion remains exactly representable at the large offset.
#[test]
fn skew_tet_preserves_force_energy_and_tangent_under_large_translation() {
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.75, 0.25, 0.5],
        [0.5, 1.5, 0.25],
        [0.25, 0.75, 1.25],
    ];
    let tets = [[0, 1, 2, 3]];
    let material = neo_hookean(2.0, 5.0);
    let p = problem(&nodes, &tets, &material, &[], &[], 1);
    let translation = [
        1_099_511_627_776.0,
        -549_755_813_888.0,
        274_877_906_944.0,
    ];

    for scale in [0.0, 1.0] {
        let u: Vec<f64> = nodes
            .iter()
            .flat_map(|&[x, y, z]| {
                [
                    scale * (0.25 * x + 0.125 * y),
                    scale * (-0.25 * y + 0.125 * z),
                    scale * (0.125 * x + 0.5 * z),
                ]
            })
            .collect();
        let translated: Vec<f64> = u
            .iter()
            .enumerate()
            .map(|(i, value)| value + translation[i % 3])
            .collect();

        with_cx(|cx| {
            assert_eq!(
                p.internal_force(&u, cx).unwrap(),
                p.internal_force(&translated, cx).unwrap()
            );
            assert_eq!(
                p.stored_energy_j(&u, cx).unwrap().to_bits(),
                p.stored_energy_j(&translated, cx).unwrap().to_bits()
            );
            let tangent = p.tangent_matrix(&u, cx).unwrap();
            let shifted_tangent = p.tangent_matrix(&translated, cx).unwrap();
            for i in 0..u.len() {
                for j in 0..u.len() {
                    assert_eq!(
                        tangent.get(i, j).to_bits(),
                        shifted_tangent.get(i, j).to_bits()
                    );
                }
            }
        });
    }
}
