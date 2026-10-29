//! IRP M1.2-W7: classical-lamination-theory section coefficients against
//! closed-form layer integrals. Every fixture is a named synthetic layup; no
//! material data is read, so these are pure CLT identities.
//!
//! The viscoelastic passivity gate (`g1_complex_abd_passivity`), the fs-plate
//! B-coupling/mass-coupling integration, and the data-gated Schelleng-1968
//! varnish trend are separate units (see the bead); this file certifies the
//! elastic section core.

use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::{
    LaminaLayer, LaminateError, LaminateSection, LayerLaw, MAX_LAMINA_LAYERS, ReferenceSurface,
};

/// Reduced plane-stress stiffness of an orthotropic lamina in the 1–2 plane,
/// the same formula `LaminateSection` uses (`OrthotropicPlaneElastic::stiffness`).
fn reduced_q(e1: f64, e2: f64, nu12: f64, g12: f64) -> [f64; 9] {
    let denom = 1.0 - nu12 * (nu12 * e2 / e1);
    let cross = nu12 * e2 / denom;
    [
        e1 / denom,
        cross,
        0.0,
        cross,
        e2 / denom,
        0.0,
        0.0,
        0.0,
        g12,
    ]
}

/// An orthotropic spruce-like lamina (L along axis 1); moduli in Pa.
fn spruce() -> OrthotropicElastic {
    // E_L, E_R, E_T; ν_LR, ν_LT, ν_RT; G_LR, G_RT, G_TL (Pa).
    OrthotropicElastic::new(
        [11.0e9, 0.9e9, 0.5e9],
        [0.37, 0.42, 0.47],
        [0.7e9, 0.65e9, 0.04e9],
        1e-3,
    )
    .unwrap()
}

/// A stiffer maple-like lamina.
fn maple() -> OrthotropicElastic {
    OrthotropicElastic::new(
        [12.0e9, 2.2e9, 1.3e9],
        [0.43, 0.49, 0.77],
        [2.0e9, 1.9e9, 0.5e9],
        1e-3,
    )
    .unwrap()
}

fn layer(
    e: OrthotropicElastic,
    thickness_m: f64,
    orientation_rad: f64,
    density: f64,
) -> LaminaLayer {
    LaminaLayer {
        thickness_m,
        orientation_rad,
        density_kg_m3: density,
        material: None,
        law: LayerLaw::Elastic(e),
    }
}

fn max_abs_9(m: &[f64; 9]) -> f64 {
    m.iter().fold(0.0_f64, |a, v| a.max(v.abs()))
}

#[test]
fn g0_single_layer_reduces_to_homogeneous_bitwise() {
    // One spruce layer about the midsurface: A = Q·h, B = 0, D = Q·h³/12,
    // integrated by the same layer formulas as the homogeneous plate section.
    let h = 0.003;
    let rho = 420.0;
    let q = reduced_q(11.0e9, 0.9e9, 0.37, 0.7e9);
    let section = LaminateSection::try_new(
        vec![layer(spruce(), h, 0.0, rho)],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    // z faces about the midsurface, formed exactly as the section forms them.
    let (z0, z1) = (0.0 - 0.5 * h, h - 0.5 * h);
    let m0 = z1 - z0;
    let m2 = (z1 * z1 * z1 - z0 * z0 * z0) / 3.0;
    for k in 0..9 {
        assert_eq!(
            section.a_pa_m()[k].to_bits(),
            (q[k] * m0).to_bits(),
            "A[{k}]"
        );
        assert_eq!(
            section.d_pa_m3()[k].to_bits(),
            (q[k] * m2).to_bits(),
            "D[{k}]"
        );
        assert_eq!(section.b_pa_m2()[k], 0.0, "B[{k}]");
    }
    // The homogeneous plate identity A = 12 D / h² (division, so relative).
    for k in 0..9 {
        let a = section.a_pa_m()[k];
        if a != 0.0 {
            assert!((12.0 * section.d_pa_m3()[k] / (h * h) / a - 1.0).abs() <= 1e-13);
        }
    }
    // Inertia: I0 = ρh, I1 = 0, I2 = ρh³/12.
    assert_eq!(section.i0_kg_m2().to_bits(), (rho * m0).to_bits());
    assert_eq!(section.i1_kg_m(), 0.0);
    assert_eq!(section.i2_kg().to_bits(), (rho * m2).to_bits());
}

#[test]
fn g0_symmetric_layup_has_zero_b() {
    // [maple(90°) / spruce(0°) / maple(90°)] about the midsurface is symmetric,
    // so its coupling B vanishes to rounding.
    let layers = vec![
        layer(maple(), 0.4e-3, core::f64::consts::FRAC_PI_2, 640.0),
        layer(spruce(), 3.0e-3, 0.0, 420.0),
        layer(maple(), 0.4e-3, core::f64::consts::FRAC_PI_2, 640.0),
    ];
    let section = LaminateSection::try_new(layers, ReferenceSurface::Midsurface).unwrap();
    let scale = max_abs_9(section.a_pa_m()) * section.total_thickness_m();
    assert!(
        max_abs_9(section.b_pa_m2()) <= 1e-16 * scale,
        "|B| {:e} > 1e-16·|A|·h {:e}",
        max_abs_9(section.b_pa_m2()),
        1e-16 * scale
    );
    assert!(section.i1_kg_m().abs() <= 1e-16 * section.i0_kg_m2() * section.total_thickness_m());
}

#[test]
fn g1_bimaterial_abd_closed_form() {
    // A spruce/maple bimaterial strip, 0° layers, about the midsurface: A, B, D
    // and the neutral-surface offset match the closed-form layer integrals.
    let (h0, h1) = (2.0e-3, 1.0e-3);
    let total = h0 + h1;
    let q0 = reduced_q(11.0e9, 0.9e9, 0.37, 0.7e9);
    let q1 = reduced_q(12.0e9, 2.2e9, 0.43, 2.0e9);
    let section = LaminateSection::try_new(
        vec![
            layer(spruce(), h0, 0.0, 420.0),
            layer(maple(), h1, 0.0, 640.0),
        ],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    // z about the midsurface: bottom layer [-t/2, -t/2+h0], top [.., t/2].
    let zb = [-0.5 * total, -0.5 * total + h0];
    let zt = [-0.5 * total + h0, 0.5 * total];
    let qs = [q0, q1];
    let (mut a, mut b, mut d) = ([0.0; 9], [0.0; 9], [0.0; 9]);
    for l in 0..2 {
        let (z0, z1) = (zb[l], zt[l]);
        let m0 = z1 - z0;
        let m1 = 0.5 * (z1 - z0) * (z1 + z0);
        let m2 = (z1 * z1 * z1 - z0 * z0 * z0) / 3.0;
        for k in 0..9 {
            a[k] += qs[l][k] * m0;
            b[k] += qs[l][k] * m1;
            d[k] += qs[l][k] * m2;
        }
    }
    let rel9 = |got: &[f64; 9], want: &[f64; 9]| {
        (0..9).fold(0.0_f64, |acc, k| {
            let denom = want[k].abs().max(max_abs_9(want));
            acc.max((got[k] - want[k]).abs() / denom)
        })
    };
    assert!(
        rel9(section.a_pa_m(), &a) <= 1e-13,
        "A rel {:e}",
        rel9(section.a_pa_m(), &a)
    );
    assert!(
        rel9(section.b_pa_m2(), &b) <= 1e-13,
        "B rel {:e}",
        rel9(section.b_pa_m2(), &b)
    );
    assert!(
        rel9(section.d_pa_m3(), &d) <= 1e-13,
        "D rel {:e}",
        rel9(section.d_pa_m3(), &d)
    );
    // Membrane neutral-surface offset e = B11/A11 (about the midsurface).
    let want_offset = b[0] / a[0];
    let got_offset = section.b_pa_m2()[0] / section.a_pa_m()[0];
    assert!((got_offset - want_offset).abs() <= 1e-13 * want_offset.abs().max(total));
    // B must be nonzero for this unsymmetric strip (mutation guard).
    assert!(max_abs_9(section.b_pa_m2()) > 0.0);
}

#[test]
fn g0_transverse_shear_from_layers() {
    // One homogeneous orthotropic layer: the energy-consistent Whitney factor
    // is exactly 5/6 in each in-plane direction, and H = (5/6)·G·h.
    let h = 0.003;
    let e = spruce();
    let section = LaminateSection::try_new(
        vec![layer(e.clone(), h, 0.0, 420.0)],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    let kappa = section.shear_factor();
    assert!((kappa[0] - 5.0 / 6.0).abs() <= 1e-14, "κ_x = {}", kappa[0]);
    assert!((kappa[1] - 5.0 / 6.0).abs() <= 1e-14, "κ_y = {}", kappa[1]);
    // H_xz = κ·G13·h, H_yz = κ·G23·h (3-axis normal: G13 = g[2], G23 = g[1]).
    let hxz = section.shear_pa_m()[0];
    let hyz = section.shear_pa_m()[3];
    assert!((hxz / (5.0 / 6.0 * e.g[2] * h) - 1.0).abs() <= 1e-13);
    assert!((hyz / (5.0 / 6.0 * e.g[1] * h) - 1.0).abs() <= 1e-13);
    assert_eq!(section.shear_pa_m()[1], 0.0);
    assert_eq!(section.shear_pa_m()[2], 0.0);
}

#[test]
fn g0_layer_cap_128_and_129() {
    // Exactly MAX_LAMINA_LAYERS is accepted; one more refuses by type.
    let make = |n: usize| {
        (0..n)
            .map(|_| layer(spruce(), 1.0e-4, 0.0, 420.0))
            .collect::<Vec<_>>()
    };
    assert!(
        LaminateSection::try_new(make(MAX_LAMINA_LAYERS), ReferenceSurface::Midsurface).is_ok()
    );
    assert!(matches!(
        LaminateSection::try_new(make(MAX_LAMINA_LAYERS + 1), ReferenceSurface::Midsurface),
        Err(LaminateError::LayerCount { count }) if count == MAX_LAMINA_LAYERS + 1
    ));
    assert!(matches!(
        LaminateSection::try_new(Vec::new(), ReferenceSurface::Midsurface),
        Err(LaminateError::LayerCount { count: 0 })
    ));
}

#[test]
fn g0_zero_or_negative_thickness_refuses() {
    for bad in [0.0, -1.0e-3, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                LaminateSection::try_new(
                    vec![layer(spruce(), bad, 0.0, 420.0)],
                    ReferenceSurface::Midsurface
                ),
                Err(LaminateError::Thickness { layer: 0 })
            ),
            "thickness {bad} must refuse"
        );
    }
    // A non-finite reference offset refuses.
    assert!(matches!(
        LaminateSection::try_new(
            vec![layer(spruce(), 1.0e-3, 0.0, 420.0)],
            ReferenceSurface::Offset(f64::NAN)
        ),
        Err(LaminateError::Reference)
    ));
    // Non-finite density or orientation refuses.
    assert!(
        LaminateSection::try_new(
            vec![layer(spruce(), 1.0e-3, f64::NAN, 420.0)],
            ReferenceSurface::Midsurface
        )
        .is_err()
    );
}

#[test]
fn g0_reference_surface_shift_rule() {
    // Shifting the reference surface up by e sends B to B − e·A and leaves A
    // unchanged (CLT first-moment transform). Compare midsurface vs an offset.
    let layers = || {
        vec![
            layer(spruce(), 2.0e-3, 0.0, 420.0),
            layer(maple(), 1.0e-3, core::f64::consts::FRAC_PI_4, 640.0),
        ]
    };
    let mid = LaminateSection::try_new(layers(), ReferenceSurface::Midsurface).unwrap();
    let total = mid.total_thickness_m();
    // Offset measured from the bottom; the midsurface sits at total/2, so an
    // offset of total/2 + e places the reference e above the midsurface.
    let e = 0.4e-3;
    let shifted =
        LaminateSection::try_new(layers(), ReferenceSurface::Offset(0.5 * total + e)).unwrap();
    // A is reference-independent.
    for k in 0..9 {
        assert!((shifted.a_pa_m()[k] - mid.a_pa_m()[k]).abs() <= 1e-10 * max_abs_9(mid.a_pa_m()));
    }
    // B' = B − e·A.
    let scale = max_abs_9(mid.b_pa_m2()).max(e * max_abs_9(mid.a_pa_m()));
    for k in 0..9 {
        let want = mid.b_pa_m2()[k] - e * mid.a_pa_m()[k];
        assert!(
            (shifted.b_pa_m2()[k] - want).abs() <= 1e-13 * scale,
            "B'[{k}] = {} want {want}",
            shifted.b_pa_m2()[k]
        );
    }
}

#[test]
fn g3_layer_reorder_symmetry() {
    // A stack that is a palindrome in (material, thickness, orientation) about
    // the midsurface is geometrically symmetric: its coupling B and inertia
    // coupling I1 vanish, and reversing the layer order reproduces every
    // coefficient bitwise (the builder is order-faithful and the symmetry is
    // exact).
    let forward = vec![
        layer(maple(), 0.4e-3, 0.3, 640.0),
        layer(spruce(), 3.0e-3, 0.0, 420.0),
        layer(maple(), 0.4e-3, 0.3, 640.0),
    ];
    let mut reversed = forward.clone();
    reversed.reverse();
    let a = LaminateSection::try_new(forward, ReferenceSurface::Midsurface).unwrap();
    let b = LaminateSection::try_new(reversed, ReferenceSurface::Midsurface).unwrap();
    for k in 0..9 {
        assert_eq!(a.a_pa_m()[k].to_bits(), b.a_pa_m()[k].to_bits(), "A[{k}]");
        assert_eq!(a.b_pa_m2()[k].to_bits(), b.b_pa_m2()[k].to_bits(), "B[{k}]");
        assert_eq!(a.d_pa_m3()[k].to_bits(), b.d_pa_m3()[k].to_bits(), "D[{k}]");
    }
    assert_eq!(a.i0_kg_m2().to_bits(), b.i0_kg_m2().to_bits());
    assert_eq!(a.i1_kg_m().to_bits(), b.i1_kg_m().to_bits());
    assert_eq!(a.i2_kg().to_bits(), b.i2_kg().to_bits());
    // The symmetry makes the coupling vanish.
    let scale = max_abs_9(a.a_pa_m()) * a.total_thickness_m();
    assert!(
        max_abs_9(a.b_pa_m2()) <= 1e-16 * scale,
        "|B| {:e}",
        max_abs_9(a.b_pa_m2())
    );
    assert!(a.i1_kg_m().abs() <= 1e-16 * a.i0_kg_m2() * a.total_thickness_m());
}

#[test]
fn g0_complex_abd_elastic_layers_are_real_and_match() {
    // Elastic layers carry no loss: complex_abd's real part equals the real
    // ABD and its imaginary part is exactly zero. (The viscoelastic passivity
    // gate with Prony layers is a separate unit.)
    let section = LaminateSection::try_new(
        vec![
            layer(spruce(), 2.0e-3, 0.0, 420.0),
            layer(maple(), 1.0e-3, core::f64::consts::FRAC_PI_6, 640.0),
        ],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    let (a, b, d) = section
        .complex_abd(2.0 * core::f64::consts::PI * 440.0)
        .unwrap();
    for k in 0..9 {
        assert_eq!(a[k].re.to_bits(), section.a_pa_m()[k].to_bits());
        assert_eq!(b[k].re.to_bits(), section.b_pa_m2()[k].to_bits());
        assert_eq!(d[k].re.to_bits(), section.d_pa_m3()[k].to_bits());
        assert_eq!(a[k].im, 0.0);
        assert_eq!(b[k].im, 0.0);
        assert_eq!(d[k].im, 0.0);
    }
}

// --- viscoelastic passivity ------------------------------------------------

use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};

/// Crate Mandel order [11,22,33,12,13,23]; a diagonal positive Mandel tensor is
/// a valid (Poisson-free) orthotropic SPD law, with shear diagonals = 2G.
fn diag6(d: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}

/// A passive orthotropic tensor-Prony law: SPD equilibrium plus one PSD Maxwell
/// branch, admitted over the audio band.
fn prony(scale: f64) -> TensorProny {
    let c_eq = diag6([
        9.0e9 * scale,
        2.0e9 * scale,
        1.0e9 * scale,
        1.5e9 * scale,
        1.3e9 * scale,
        0.9e9 * scale,
    ]);
    let branch = PronyBranch::new(
        diag6([1.0e9, 3.0e8, 2.0e8, 2.5e8, 2.2e8, 1.6e8].map(|v| v * scale)),
        1.0e-4,
    )
    .expect("PSD branch");
    TensorProny::new(
        c_eq,
        None,
        vec![branch],
        (
            2.0 * core::f64::consts::PI * 20.0,
            2.0 * core::f64::consts::PI * 20_000.0,
        ),
    )
    .expect("passive orthotropic Prony")
}

fn prony_layer(
    p: TensorProny,
    thickness_m: f64,
    orientation_rad: f64,
    density: f64,
) -> LaminaLayer {
    LaminaLayer {
        thickness_m,
        orientation_rad,
        density_kg_m3: density,
        material: None,
        law: LayerLaw::Prony(p),
    }
}

#[test]
fn g1_complex_abd_passivity() {
    use fs_math::det;
    // An unsymmetric three-layer viscoelastic stack. The assembled complex
    // laminate operator [[A,B],[B,D]] must have a PSD imaginary part (passive
    // dissipation) at every frequency in the source band: each layer gives
    // Im Q_k ⊗ [[m0,m1],[m1,m2]], a Kronecker of a PSD reduced-loss matrix and a
    // PSD second-moment matrix, summed over layers. B's own block need not be
    // sign-definite; the 6×6 operator is the physical passivity statement.
    let section = LaminateSection::try_new(
        vec![
            prony_layer(prony(1.0), 2.0e-3, 0.0, 420.0),
            prony_layer(prony(1.7), 1.0e-3, core::f64::consts::PI / 5.0, 640.0),
            prony_layer(prony(0.6), 0.5e-3, -core::f64::consts::FRAC_PI_3, 500.0),
        ],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    let (lo, hi) = (
        2.0 * core::f64::consts::PI * 20.0,
        2.0 * core::f64::consts::PI * 20_000.0,
    );
    let mut s = fs_rand::StreamKey {
        seed: 0x4c_414d_5f56_31,
        kernel: 1,
        tile: 0,
    }
    .stream();
    let (mut worst_eval, mut max_loss) = (f64::INFINITY, 0.0_f64);
    for _ in 0..1000 {
        // Log-uniform frequency strictly inside the band.
        let omega = det::exp(det::ln(lo) + s.next_f64() * (det::ln(hi) - det::ln(lo)));
        let (a, b, d) = section.complex_abd(omega).unwrap();
        // Assemble the symmetric 6×6 imaginary operator, engineering-shear
        // throughout (the per-layer plane_stress_complex check uses the same
        // convention). Top-left Im A, off-diagonal Im B, bottom-right Im D.
        let mut im = [0.0f64; 36];
        let mut scale = 0.0f64;
        for r in 0..3 {
            for c in 0..3 {
                let (ia, ib, id) = (a[3 * r + c].im, b[3 * r + c].im, d[3 * r + c].im);
                im[6 * r + c] = ia;
                im[6 * r + (c + 3)] = ib;
                im[6 * (r + 3) + c] = ib;
                im[6 * (r + 3) + (c + 3)] = id;
                scale = scale.max(ia.abs()).max(id.abs());
                max_loss = max_loss.max(ia.abs()).max(id.abs());
            }
        }
        let (evals, _) = fs_la::eigen::jacobi_eigh(&im, 6);
        let min_eval = evals.iter().copied().fold(f64::INFINITY, f64::min);
        // PSD to a relative floor: the per-layer loss is PSD by construction and
        // the moment matrices are PSD, so the only negativity is rounding.
        assert!(
            min_eval >= -1e-9 * scale.max(1.0),
            "Im[[A,B],[B,D]] not PSD at ω = {omega}: min eval {min_eval:e}, scale {scale:e}"
        );
        worst_eval = worst_eval.min(min_eval / scale.max(1.0));
    }
    // Liveness: the stack actually dissipates (nonzero loss) in band.
    assert!(max_loss > 0.0, "viscoelastic stack shows no loss");
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W7\",\"test\":\"g1_complex_abd_passivity\",\"tier\":\"G1\",\"inputs\":{{\"layers\":3,\"band_rad_s\":[{lo:e},{hi:e}],\"frequencies\":1000}},\"provenance\":\"synthetic orthotropic Prony (not measured)\",\"observed\":{{\"worst_relative_min_eigenvalue\":{worst_eval:e},\"max_loss_pa_m\":{max_loss:e}}},\"tolerance\":{{\"psd_relative_floor\":1e-9}},\"verdict\":\"pass\"}}"
    );
}
