//! Koiter post-buckling battery (bead tfz.15, feature
//! `koiter-asymptotics` — [F], off by default).
//!
//! stab-006: the Euler column classifies SYMMETRIC STABLE (a ≈ 0,
//! b > 0), and the SAMPLED-CONTINUATION FALLBACK ORACLE agrees: an
//! imperfect column shows no limit point below the critical load —
//! the imperfection-tolerant signature.

use fs_solid::continuation::{ArcSettings, PathEvent, PathResidual, PathState, advance};
use fs_solid::koiter::{Bifurcation, koiter_coefficients};
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::{
    HyperProblem, Mesh2, NewtonSettings, Patch, buckling_loads, expand_mode, reduced_pencil,
};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

#[test]
fn stab_006_koiter_column_symmetric_stable() {
    let l = 8.0;
    let t = 0.4;
    let p = 1e-4;
    let mesh = Mesh2::quads(l, t, 30, 3);
    let linear = LinearProblem {
        mesh: &mesh,
        youngs: 1.0,
        poisson: 0.0,
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-p, 0.0])],
    };
    let (k, kg, dof_map, _) = reduced_pencil(&linear).expect("pencil builds");
    let pencil = buckling_loads(&k, &kg, &dof_map, 1, 400).expect("pencil solves");
    let lambda_cr = pencil.loads[0];
    let mode = expand_mode(&pencil.modes[0], &pencil.dof_map);
    let card = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 0.5,
            lambda: 0.0,
        },
        4.0,
    )
    .expect("card");
    let problem = HyperProblem {
        mesh: &mesh,
        material: &card,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-p, 0.0])],
        settings: NewtonSettings {
            load_steps: 4,
            ..NewtonSettings::default()
        },
    };
    // Critical state on the fundamental branch (axial compression at
    // λ_cr, which for the column is essentially the linear state).
    let (u_cr_nodal, _) = {
        let scaled = HyperProblem {
            mesh: &mesh,
            material: &card,
            dirichlet: problem.dirichlet.clone(),
            traction: vec![(Patch::Right, &|_, _| [-p, 0.0])],
            settings: NewtonSettings {
                load_steps: 4,
                ..NewtonSettings::default()
            },
        };
        // Solve at unit load then scale: the fundamental branch is
        // linear to high accuracy well below wrinkling strains.
        let (u1, rep) = scaled.solve().expect("fundamental state");
        (u1, rep)
    };
    let mut u_cr = vec![0.0f64; 2 * mesh.node_count()];
    for (node, v) in u_cr_nodal.iter().enumerate() {
        u_cr[2 * node] = lambda_cr * v[0];
        u_cr[2 * node + 1] = lambda_cr * v[1];
    }
    let coeff = koiter_coefficients(&problem, &u_cr, lambda_cr, &mode, 5e-3, 0.05)
        .expect("coefficients evaluate");
    // Fallback oracle: an imperfect column traced through λ_cr shows
    // no limit point below it (imperfection-tolerant).
    // Imperfect GEOMETRY (a bowed centerline — state perturbations
    // would just relax back to the perfect fundamental branch): the
    // clamped-free first-mode shape at 5% of the thickness.
    let mut bowed = mesh.clone();
    for node in &mut bowed.nodes {
        let x = node[0];
        node[1] += 0.05 * t * (1.0 - (std::f64::consts::PI * x / (2.0 * l)).cos());
    }
    let imperfect = HyperProblem {
        mesh: &bowed,
        material: &card,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-p, 0.0])],
        settings: NewtonSettings::default(),
    };
    let settings = ArcSettings {
        ds: 0.4,
        ds_max: 1.0,
        ..ArcSettings::default()
    };
    let mut path = PathState::start(PathResidual::ndof(&imperfect), &settings);
    advance(&imperfect, &mut path, &settings, 50).expect("imperfect path advances");
    let early_limit = path.events.iter().any(|e| match e {
        PathEvent::LimitPoint { lambda, .. } => *lambda < 0.98 * lambda_cr,
        PathEvent::BranchPoint { .. } => false,
    });
    let reached = path.lambda > 0.9 * lambda_cr;
    let pass = coeff.class == Bifurcation::SymmetricStable
        && coeff.b > 0.0
        && !early_limit
        && reached;
    verdict(
        "stab-006",
        pass,
        &format!(
            "\"detail\":\"Koiter FD expansion + sampled-continuation oracle\",\
             \"a\":{:.3e},\"b\":{:.3e},\"class\":\"{:?}\",\
             \"lambda_cr\":{lambda_cr:.4},\"path_lambda_end\":{:.4},\
             \"early_limit_point\":{early_limit}",
            coeff.a, coeff.b, coeff.class, path.lambda
        ),
    );
}
