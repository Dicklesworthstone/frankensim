//! A1 unsteady LTI assembly (bead wf-root-guzez.5.8.1, E4.3b2-i): the
//! FOM state-space the shared-basis rational-Krylov reduction (5.8.2)
//! projects. At each FROZEN operating point (E4.3b1 grid) the A1 lane
//! is exactly linear: per-station Wagner–Jones lag pairs in reduced
//! time, spanwise/inter-surface coupling through the IMAGE-AWARE
//! prescribed-wake operator, thin-airfoil circulation closure.
//!
//!   states  x  = 2 lag states per station (Wagner–Jones a/b pairs)
//!   inputs  u  = [δ_canard, α_gust]
//!   outputs y  = [wing lift, canard lift, canard hinge moment]
//!
//! Wiring (the algebraic loop eliminated in closed form):
//!   α_eff = (1−ā)·v + a₁x₁ + a₂x₂         (per station, ā = a₁+a₂)
//!   Γ     = diag(π c V)·α_eff             (thin airfoil)
//!   w/V   = M·α_eff,  M = W·diag(π c)     (prescribed-wake operator)
//!   v     = α_geo + w/V                    (lag input feedback)
//!   ⇒ (I − (1−ā)M)·α_eff = (1−ā)·α_geo + a₁x₁ + a₂x₂
//!   ẋᵢⱼ  = (2V/cᵢ)·bⱼ·(vᵢ − xᵢⱼ)
//!
//! The LTI is CIRCULATORY-ONLY (declared): the apparent-mass share is
//! a static feed-through handled outside the reduction, so the ROM
//! never has to approximate a derivative feed-through.

use crate::images::CertifiedGround;
use crate::prescribedwake::{WakeOperatingPoint, assemble_operator};
use crate::{Refusal, refuse};
use fs_airfoil::indicial::WAGNER_JONES;
use fs_blake3::hash_domain;

/// One A1 lifting-line layout (FRD: +x forward, +z down; the wake
/// convects toward −x; the canard sits AHEAD at +x, mounted low).
#[derive(Clone, Debug, PartialEq)]
pub struct A1Layout {
    /// Shed points (quarter-chord line), FRD.
    pub shed_points: Vec<[f64; 3]>,
    /// Collocation probes (three-quarter-chord line), FRD.
    pub probes: Vec<[f64; 3]>,
    /// Station chords [m].
    pub chords_m: Vec<f64>,
    /// Station widths [m].
    pub widths_m: Vec<f64>,
    /// Canard flags.
    pub is_canard: Vec<bool>,
    /// Reference span [m] (h/b scaling).
    pub span_b_m: f64,
}

/// The registered v1 layout (flyer-reference lineage: wing 12.29 ×
/// 1.981 m as 8 stations; canard 3.66 × 0.61 m as 4 stations, qc
/// 2.23 m ahead, 0.7 m below the wing plane).
#[must_use]
pub fn wright_a1_layout_v1() -> A1Layout {
    let mut shed = Vec::new();
    let mut probes = Vec::new();
    let mut chords = Vec::new();
    let mut widths = Vec::new();
    let mut canard = Vec::new();
    let mut surface = |span: f64, chord: f64, x: f64, z: f64, n: usize, is_c: bool| {
        let dy = span / n as f64;
        for i in 0..n {
            let y = -span / 2.0 + (i as f64 + 0.5) * dy;
            shed.push([x, y, z]);
            probes.push([x - 0.5 * chord, y, z]);
            chords.push(chord);
            widths.push(dy);
            canard.push(is_c);
        }
    };
    surface(12.29, 1.981, 0.0, 0.0, 8, false);
    surface(3.66, 0.61, 2.23, 0.7, 4, true);
    A1Layout {
        shed_points: shed,
        probes,
        chords_m: chords,
        widths_m: widths,
        is_canard: canard,
        span_b_m: 12.29,
    }
}

/// The assembled FOM at one frozen point.
#[derive(Clone, Debug, PartialEq)]
pub struct A1Lti {
    /// Stations.
    pub n_stations: usize,
    /// State order (2 per station).
    pub order: usize,
    /// A, row-major order×order.
    pub a: Vec<f64>,
    /// B, row-major order×2.
    pub b: Vec<f64>,
    /// C, row-major 3×order.
    pub c: Vec<f64>,
    /// D, row-major 3×2 (direct circulatory feed-through).
    pub d: [f64; 6],
    /// The operating point.
    pub point: WakeOperatingPoint,
    /// Content digest over (A, B, C, D) bits.
    pub digest: String,
}

/// Gauss solve AX = B for X (B holds multiple columns), in place.
pub(crate) fn solve_multi(
    a: &mut [f64],
    b: &mut [f64],
    n: usize,
    ncols: usize,
) -> Result<(), Refusal> {
    for col in 0..n {
        let mut piv = col;
        for row in (col + 1)..n {
            if a[row * n + col].abs() > a[piv * n + col].abs() {
                piv = row;
            }
        }
        if a[piv * n + col].abs() < 1e-12 {
            return Err(refuse(
                "a1-loop-singular",
                format!("algebraic loop pivot {col} below floor"),
                "the (I − (1−ā)M) loop matrix is singular at this point",
            ));
        }
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            for k in 0..ncols {
                b.swap(col * ncols + k, piv * ncols + k);
            }
        }
        let d = a[col * n + col];
        for row in (col + 1)..n {
            let f = a[row * n + col] / d;
            if f != 0.0 {
                for k in col..n {
                    a[row * n + k] -= f * a[col * n + k];
                }
                for k in 0..ncols {
                    b[row * ncols + k] -= f * b[col * ncols + k];
                }
            }
        }
    }
    for col in (0..n).rev() {
        for k in 0..ncols {
            let mut s = b[col * ncols + k];
            for j in (col + 1)..n {
                s -= a[col * n + j] * b[j * ncols + k];
            }
            b[col * ncols + k] = s / a[col * n + col];
        }
    }
    Ok(())
}

/// Assemble the A1 FOM at one frozen operating point.
///
/// # Errors
/// `a1-speed-invalid` (v outside [5, 40] — cap AND cap+1);
/// `a1-loop-singular`; prescribed-wake refusals pass through.
pub fn assemble_a1_lti(
    layout: &A1Layout,
    point: &WakeOperatingPoint,
    ground: &CertifiedGround,
    v_mps: f64,
    wake_rows: usize,
) -> Result<A1Lti, Refusal> {
    if !(v_mps.is_finite() && (5.0..=40.0).contains(&v_mps)) {
        return Err(refuse(
            "a1-speed-invalid",
            format!("v {v_mps} outside [5, 40]"),
            "the A1 grid is registered for the flight-speed class",
        ));
    }
    let n = layout.shed_points.len();
    let row_dx = v_mps * (1.0 / 120.0);
    let op = assemble_operator(
        &layout.shed_points,
        &layout.probes,
        ground,
        point,
        wake_rows,
        row_dx,
        layout.span_b_m,
    )?;
    // M = W · diag(π c) (dimensionless downwash-per-α map). The
    // operator's w_normal is +z-down FRD; positive circulation lifts,
    // and downwash REDUCES effective α — the operator rows already
    // carry the FRD sign, consumed as w/V added to α_geo.
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = op.w_normal[i * n + j] * core::f64::consts::PI * layout.chords_m[j];
        }
    }
    let a_w = WAGNER_JONES.a;
    let b_w = WAGNER_JONES.b;
    let abar = a_w[0] + a_w[1];
    // Loop matrix L = I − (1−ā)·M; α_eff = L⁻¹[(1−ā)α_geo + a₁x₁ + a₂x₂].
    // Solve L · X = RHS for the two operators we need:
    //   Ge = L⁻¹·(1−ā)   (n×n, then applied to α_geo)
    //   Gx = L⁻¹          (n×n, applied to the lag-state sums)
    let mut loop_m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            loop_m[i * n + j] = f64::from(u8::from(i == j)) - (1.0 - abar) * m[i * n + j];
        }
    }
    // RHS = identity → X = L⁻¹.
    let mut linv = vec![0.0f64; n * n];
    for (i, row) in linv.chunks_mut(n).enumerate() {
        row[i] = 1.0;
    }
    solve_multi(&mut loop_m.clone(), &mut linv, n, n)?;
    // α_geo per input: station i gets u₀ (δc) iff canard, u₁ (gust) always.
    let geo = |i: usize, input: usize| -> f64 {
        match input {
            0 => f64::from(u8::from(layout.is_canard[i])),
            _ => 1.0,
        }
    };
    // α_eff = linv · [(1−ā)·α_geo + a₁x₁ + a₂x₂]
    // v      = α_geo + M·α_eff
    // ẋᵢⱼ   = λᵢⱼ(vᵢ − xᵢⱼ),  λᵢⱼ = (2V/cᵢ)bⱼ
    // State packing: x[2i] = station i pole 1, x[2i+1] = pole 2.
    let order = 2 * n;
    let mut a_mat = vec![0.0f64; order * order];
    let mut b_mat = vec![0.0f64; order * 2];
    // ∂α_eff/∂x[2k+j] = linv[:,k]·a_w[j]; ∂α_eff/∂u = linv·(1−ā)·geo.
    for i in 0..n {
        for jpole in 0..2 {
            let lam = (2.0 * v_mps / layout.chords_m[i]) * b_w[jpole];
            let row = 2 * i + jpole;
            // v_i = α_geo,i + Σ_p M[i,p] α_eff,p
            // ∂v_i/∂x[2k+q] = Σ_p M[i,p] linv[p,k] a_w[q]
            for k in 0..n {
                let mut mlik = 0.0;
                for p in 0..n {
                    mlik += m[i * n + p] * linv[p * n + k];
                }
                for q in 0..2 {
                    a_mat[row * order + 2 * k + q] += lam * mlik * a_w[q];
                }
            }
            // The −x self term.
            a_mat[row * order + row] -= lam;
            // ∂v_i/∂u = geo(i,·) + Σ_p M[i,p]·(1−ā)·Σ_r linv[p,r] geo(r,·)
            for input in 0..2 {
                let mut s = geo(i, input);
                for p in 0..n {
                    let mut lr = 0.0;
                    for r in 0..n {
                        lr += linv[p * n + r] * geo(r, input);
                    }
                    s += m[i * n + p] * (1.0 - abar) * lr;
                }
                b_mat[row * 2 + input] = lam * s;
            }
        }
    }
    // Outputs: y = T·α_eff with T rows = [wing lift, canard lift, hinge].
    // L_i = ρV·Γ_i·width = ρ V² π c_i width_i α_eff,i (ρ = 1.294 reference,
    // declared — the reduction is scale-invariant, ρ only scales C/D).
    let rho = 1.294;
    let mut t = vec![0.0f64; 3 * n];
    for i in 0..n {
        let li =
            rho * v_mps * v_mps * core::f64::consts::PI * layout.chords_m[i] * layout.widths_m[i];
        if layout.is_canard[i] {
            t[n + i] = li;
            t[2 * n + i] = -li * 0.15 * layout.chords_m[i]; // 40 %-chord hinge arm
        } else {
            t[i] = li;
        }
    }
    // C = T·linv·[a₁, a₂ per station]; D = T·linv·(1−ā)·geo.
    let mut c_mat = vec![0.0f64; 3 * order];
    let mut d_mat = [0.0f64; 6];
    for out in 0..3 {
        for k in 0..n {
            let mut tl = 0.0;
            for p in 0..n {
                tl += t[out * n + p] * linv[p * n + k];
            }
            for q in 0..2 {
                c_mat[out * order + 2 * k + q] = tl * a_w[q];
            }
        }
        for input in 0..2 {
            let mut s = 0.0;
            for p in 0..n {
                let mut lr = 0.0;
                for r in 0..n {
                    lr += linv[p * n + r] * geo(r, input);
                }
                s += t[out * n + p] * (1.0 - abar) * lr;
            }
            d_mat[out * 2 + input] = s;
        }
    }
    let mut bytes = Vec::new();
    for v in a_mat.iter().chain(b_mat.iter()).chain(c_mat.iter()) {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    for v in d_mat {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest = hash_domain("org.frankensim.fs-wing.a1-lti.v1", &bytes).to_hex();
    Ok(A1Lti {
        n_stations: n,
        order,
        a: a_mat,
        b: b_mat,
        c: c_mat,
        d: d_mat,
        point: *point,
        digest,
    })
}

impl A1Lti {
    /// Explicit-Euler simulation (battery oracle driver — dense FOM,
    /// deterministic; dt must resolve the fastest lag).
    ///
    /// # Errors
    /// `a1-sim-invalid` (dt/steps caps at cap AND cap+1).
    pub fn simulate(
        &self,
        u: &dyn Fn(usize) -> [f64; 2],
        dt_s: f64,
        n_steps: usize,
    ) -> Result<Vec<[f64; 3]>, Refusal> {
        if !(dt_s.is_finite() && dt_s > 0.0 && dt_s <= 0.01 && (1..=200_000).contains(&n_steps)) {
            return Err(refuse(
                "a1-sim-invalid",
                format!("dt {dt_s}, steps {n_steps}"),
                "dt (0, 0.01]; steps [1, 200000]",
            ));
        }
        let n = self.order;
        let mut x = vec![0.0f64; n];
        let mut out = Vec::with_capacity(n_steps);
        for k in 0..n_steps {
            let uk = u(k);
            let mut y = [0.0f64; 3];
            for (o, yo) in y.iter_mut().enumerate() {
                let mut s = self.d[o * 2] * uk[0] + self.d[o * 2 + 1] * uk[1];
                for j in 0..n {
                    s += self.c[o * n + j] * x[j];
                }
                *yo = s;
            }
            out.push(y);
            let mut dx = vec![0.0f64; n];
            for i in 0..n {
                let mut s = self.b[i * 2] * uk[0] + self.b[i * 2 + 1] * uk[1];
                for j in 0..n {
                    s += self.a[i * n + j] * x[j];
                }
                dx[i] = s;
            }
            for i in 0..n {
                x[i] += dt_s * dx[i];
            }
        }
        Ok(out)
    }

    /// Steady output for a held input (the DC oracle's independent
    /// recompute: α_eff = (I − M)⁻¹ α_geo directly — no A/B/C/D).
    ///
    /// # Errors
    /// `a1-loop-singular`.
    pub fn dc_direct(
        layout: &A1Layout,
        point: &WakeOperatingPoint,
        ground: &CertifiedGround,
        v_mps: f64,
        wake_rows: usize,
        u: [f64; 2],
    ) -> Result<[f64; 3], Refusal> {
        let n = layout.shed_points.len();
        let row_dx = v_mps * (1.0 / 120.0);
        let op = assemble_operator(
            &layout.shed_points,
            &layout.probes,
            ground,
            point,
            wake_rows,
            row_dx,
            layout.span_b_m,
        )?;
        let mut im = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                im[i * n + j] = f64::from(u8::from(i == j))
                    - op.w_normal[i * n + j] * core::f64::consts::PI * layout.chords_m[j];
            }
        }
        let mut alpha = vec![0.0f64; n];
        for (i, a) in alpha.iter_mut().enumerate() {
            *a = f64::from(u8::from(layout.is_canard[i])) * u[0] + u[1];
        }
        solve_multi(&mut im, &mut alpha, n, 1)?;
        let rho = 1.294;
        let mut y = [0.0f64; 3];
        for i in 0..n {
            let li = rho
                * v_mps
                * v_mps
                * core::f64::consts::PI
                * layout.chords_m[i]
                * layout.widths_m[i]
                * alpha[i];
            if layout.is_canard[i] {
                y[1] += li;
                y[2] += -li * 0.15 * layout.chords_m[i];
            } else {
                y[0] += li;
            }
        }
        Ok(y)
    }
}
