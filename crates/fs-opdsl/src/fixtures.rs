//! The acceptance-criteria operators, defined IN the IR: Poisson
//! (pure FEEC derivation), convection–diffusion (FEEC diffusion +
//! hand advection through the escape hatch — nonsymmetric, so the
//! transpose gates do real work), and linear elasticity (constitutive
//! integration: the element stiffness comes from fs-material's
//! `IsotropicElastic` tangent, entering as a hand atom with the same
//! obligations). Plus a nonlinear reaction–diffusion built on the
//! Poisson chain with the cubic pointwise law.

use crate::atoms::{Atom, Transpose};
use crate::expr::{Expr, OperatorDef, Space};
use crate::law::CubicReaction;
use fs_feec::ElementGeometry;
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::{Coo, Csr};

/// Poisson: R(u) = d0ᵀ·M1·d0·u — the pure FEEC derivation (every
/// atom `derived`).
///
/// # Panics
/// Never for a valid complex (the checked constructors cannot fail on
/// this chain).
#[must_use]
pub fn poisson(complex: &TetComplex, geo: &ElementGeometry) -> (OperatorDef, Expr) {
    let dims = Dims::NONE;
    let space = Space {
        degree: 0,
        n: complex.vertex_count,
        dims,
    };
    let mut def = OperatorDef::new(space);
    let d0 = def.add_atom(Atom::d(complex, 0, dims));
    let m1 = def.add_atom(Atom::mass(complex, geo, 1, dims));
    let d0t = def.add_atom(Atom::d_transposed(complex, 0, dims));
    let u = def.field();
    let expr = def
        .apply(d0, u)
        .and_then(|e| def.apply(m1, e))
        .and_then(|e| def.apply(d0t, e))
        .expect("poisson chain is well-typed");
    (def, expr)
}

/// Hand-assembled advection matrix A[i][j] = Σ_T (c·∇λ_j)·|V|/4 for
/// constant velocity c (P1 Galerkin, ∫λ_i = V/4). Deterministic COO
/// accumulation.
#[must_use]
pub fn advection_matrix(complex: &TetComplex, geo: &ElementGeometry, velocity: [f64; 3]) -> Csr {
    let n = complex.vertex_count;
    let mut coo = Coo::new(n, n);
    for (m, tet) in complex.tets.iter().enumerate() {
        let quarter_vol = geo.vol_signed[m].abs() / 4.0;
        for i in 0..4 {
            for j in 0..4 {
                let g = geo.grads[m][j];
                let cdotg =
                    velocity[0].mul_add(g[0], velocity[1].mul_add(g[1], velocity[2] * g[2]));
                coo.push(tet[i] as usize, tet[j] as usize, quarter_vol * cdotg);
            }
        }
    }
    coo.assemble()
}

/// Convection–diffusion: R(u) = ν·(d0ᵀ M1 d0)u + A_c·u. The advection
/// atom is HAND (escape hatch): its transpose is derived at
/// registration and it passes the same adjoint gates as everything
/// else — nonsymmetric, so those gates are not vacuous here.
///
/// # Panics
/// Never for a valid complex.
#[must_use]
pub fn convection_diffusion(
    complex: &TetComplex,
    geo: &ElementGeometry,
    nu: f64,
    velocity: [f64; 3],
) -> (OperatorDef, Expr) {
    let dims = Dims::NONE;
    let space = Space {
        degree: 0,
        n: complex.vertex_count,
        dims,
    };
    let mut def = OperatorDef::new(space);
    let d0 = def.add_atom(Atom::d(complex, 0, dims));
    let m1 = def.add_atom(Atom::mass(complex, geo, 1, dims));
    let d0t = def.add_atom(Atom::d_transposed(complex, 0, dims));
    let adv = def.add_atom(Atom::external(
        "advection",
        advection_matrix(complex, geo, velocity),
        Transpose::Derived,
        space,
        space,
    ));
    let diffusion = def
        .apply(d0, def.field())
        .and_then(|e| def.apply(m1, e))
        .and_then(|e| def.apply(d0t, e))
        .expect("diffusion chain is well-typed");
    let expr = def
        .add(
            def.scale(nu, diffusion),
            def.apply(adv, def.field()).expect("typed"),
        )
        .expect("spaces agree");
    (def, expr)
}

/// Reaction–diffusion: R(u) = (d0ᵀ M1 d0)u + M0·N(u) with the cubic
/// pointwise law N(u) = α·u³ — the nonlinear chain-rule fixture (mass-
/// weighted so the reaction enters the residual consistently).
///
/// # Panics
/// Never for a valid complex.
#[must_use]
pub fn reaction_diffusion(
    complex: &TetComplex,
    geo: &ElementGeometry,
    alpha: f64,
) -> (OperatorDef, Expr) {
    let dims = Dims::NONE;
    let space = Space {
        degree: 0,
        n: complex.vertex_count,
        dims,
    };
    let mut def = OperatorDef::new(space);
    let d0 = def.add_atom(Atom::d(complex, 0, dims));
    let m1 = def.add_atom(Atom::mass(complex, geo, 1, dims));
    let d0t = def.add_atom(Atom::d_transposed(complex, 0, dims));
    let m0 = def.add_atom(Atom::mass(complex, geo, 0, dims));
    let law = def.add_law(Box::new(CubicReaction { alpha }));
    let diffusion = def
        .apply(d0, def.field())
        .and_then(|e| def.apply(m1, e))
        .and_then(|e| def.apply(d0t, e))
        .expect("diffusion chain is well-typed");
    let reaction = def
        .pointwise(law, def.field())
        .and_then(|e| def.apply(m0, e))
        .expect("reaction chain is well-typed");
    let expr = def.add(diffusion, reaction).expect("spaces agree");
    (def, expr)
}

/// Hand-assembled linear-elasticity stiffness: K = Σ_T |V|·B_aᵀ·C·B_b
/// with C from fs-material's `IsotropicElastic` tangent (evaluated at
/// zero strain — the linear regime) and B the P1 strain-displacement
/// blocks from the barycentric gradients. Vector dofs are numbered
/// 3·vertex + component; Voigt order (xx, yy, zz, xy, yz, zx) with
/// engineering shear strains, matching fs-material's convention.
#[must_use]
pub fn elasticity_stiffness(
    complex: &TetComplex,
    geo: &ElementGeometry,
    law: &IsotropicElastic,
) -> Csr {
    let c = law.tangent(&[0.0; 6], &());
    let n = 3 * complex.vertex_count;
    let mut coo = Coo::new(n, n);
    for (m, tet) in complex.tets.iter().enumerate() {
        let vol = geo.vol_signed[m].abs();
        // Bᵀ_a stored column-major (3 rows of 6): row i is the strain
        // sensitivity of displacement component i at node a.
        let bt = |a: usize| -> [[f64; 6]; 3] {
            let g = geo.grads[m][a];
            [
                [g[0], 0.0, 0.0, g[1], 0.0, g[2]],
                [0.0, g[1], 0.0, g[0], g[2], 0.0],
                [0.0, 0.0, g[2], 0.0, g[1], g[0]],
            ]
        };
        for a in 0..4 {
            let bta = bt(a);
            for bb_idx in 0..4 {
                let btb = bt(bb_idx);
                // K_ab = vol · B_aᵀ C B_b (3×3).
                for (i, bai) in bta.iter().enumerate() {
                    for (j, bbj) in btb.iter().enumerate() {
                        let mut acc = 0.0f64;
                        for (p, baip) in bai.iter().enumerate() {
                            for (q, bbjq) in bbj.iter().enumerate() {
                                acc = (baip * c[p][q]).mul_add(*bbjq, acc);
                            }
                        }
                        coo.push(
                            3 * tet[a] as usize + i,
                            3 * tet[bb_idx] as usize + j,
                            vol * acc,
                        );
                    }
                }
            }
        }
    }
    coo.assemble()
}

/// Raw-vector dof space marker (not a cochain degree).
pub const VECTOR_DOFS: u8 = 255;

/// Linear elasticity: R(u) = K·u with K constitutive-assembled (hand
/// atom, symmetric by declaration — and VERIFIED symmetric by the
/// adjoint gate, which is the point of making hand atoms pass the
/// same checks).
///
/// # Panics
/// If the material parameters are invalid.
#[must_use]
pub fn elasticity(
    complex: &TetComplex,
    geo: &ElementGeometry,
    youngs: f64,
    poisson_ratio: f64,
) -> (OperatorDef, Expr) {
    let law = IsotropicElastic::new(youngs, poisson_ratio, 1.0).expect("valid material");
    let space = Space {
        degree: VECTOR_DOFS,
        n: 3 * complex.vertex_count,
        dims: Dims::NONE,
    };
    let mut def = OperatorDef::new(space);
    let k = def.add_atom(Atom::external(
        "elastic-stiffness",
        elasticity_stiffness(complex, geo, &law),
        Transpose::Symmetric,
        space,
        space,
    ));
    let expr = def.apply(k, def.field()).expect("typed");
    (def, expr)
}
