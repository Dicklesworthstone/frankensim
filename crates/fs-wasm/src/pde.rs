//! pde.rs — Tier-2 PDE / FEM demos, all driven by the REAL kernels:
//!
//! * [`topopt_frames`]   — SIMP topology optimization (fs-sparse assembly +
//!   Jacobi-preconditioned CG on `spmv` each iteration).
//! * [`wave2d_frames`]   — a spectral 2D wave solve (fs-fft `Fft`, Laplacian in
//!   Fourier space, leapfrog in time).
//! * [`gray_scott_frames`] — Gray–Scott reaction–diffusion (fs-sparse periodic
//!   Laplacian applied by `spmv`).
//! * [`fluid_frames`]    — a 2D stable-fluids smoke sim whose pressure
//!   projection is a real fs-sparse CG Poisson solve.
//!
//! Every input is clamped and every loop is capped: nothing here can trap.

use crate::{cg, dot, laplacian_5pt};
use fs_fft::{Fft, C64};
use fs_math::det;
use fs_rand::StreamKey;
use fs_sparse::{Coo, Csr};

/* ----------------------------------------------------------------------- */
/*  Shared helpers                                                          */
/* ----------------------------------------------------------------------- */

/// Jacobi (diagonal) preconditioned conjugate gradients on a `Csr` operator.
/// `inv_diag[i] = 1/A[i,i]`. Matrix-free via `spmv`; the iteration count is
/// hard-capped so a badly conditioned system can never spin forever.
pub(crate) fn pcg(a: &Csr, b: &[f64], inv_diag: &[f64], maxit: usize, tol: f64) -> Vec<f64> {
    let m = b.len();
    let mut x = vec![0.0f64; m];
    let mut r = b.to_vec();
    let mut z: Vec<f64> = r.iter().zip(inv_diag).map(|(&ri, &d)| ri * d).collect();
    let mut p = z.clone();
    let mut ap = vec![0.0f64; m];
    let mut rz = dot(&r, &z);
    let bnorm = dot(b, b).sqrt().max(1e-30);
    for _ in 0..maxit {
        a.spmv(&p, &mut ap);
        let denom = dot(&p, &ap);
        if denom.abs() < 1e-300 {
            break;
        }
        let alpha = rz / denom;
        for i in 0..m {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        if dot(&r, &r).sqrt() / bnorm < tol {
            break;
        }
        for i in 0..m {
            z[i] = inv_diag[i] * r[i];
        }
        let rz_new = dot(&r, &z);
        if rz.abs() < 1e-300 {
            break;
        }
        let beta = rz_new / rz;
        for i in 0..m {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    x
}

/// The periodic 5-point Laplacian `Δ` (diagonal −4, neighbours +1, wrap-around)
/// as a real fs-sparse CSR — the diffusion operator applied by `spmv`.
pub(crate) fn laplacian_5pt_periodic(n: usize) -> Csr {
    let m = n * n;
    let mut coo = Coo::new(m, m);
    let idx = |i: usize, j: usize| i * n + j;
    for i in 0..n {
        for j in 0..n {
            let k = idx(i, j);
            coo.push(k, k, -4.0);
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            let jp = (j + 1) % n;
            let jm = (j + n - 1) % n;
            coo.push(k, idx(ip, j), 1.0);
            coo.push(k, idx(im, j), 1.0);
            coo.push(k, idx(i, jp), 1.0);
            coo.push(k, idx(i, jm), 1.0);
        }
    }
    coo.assemble()
}

/* ----------------------------------------------------------------------- */
/*  1 · topopt_frames — SIMP topology optimization (the marquee)            */
/* ----------------------------------------------------------------------- */

/// Bilinear-quad plane-stress element stiffness (E = 1, ν = 0.3), node order
/// (top-left, top-right, bottom-right, bottom-left) — the classic 8×8 matrix.
fn element_stiffness() -> [[f64; 8]; 8] {
    let nu = 0.3f64;
    let k = [
        0.5 - nu / 6.0,
        0.125 + nu / 8.0,
        -0.25 - nu / 12.0,
        -0.125 + 3.0 * nu / 8.0,
        -0.25 + nu / 12.0,
        -0.125 - nu / 8.0,
        nu / 6.0,
        0.125 - 3.0 * nu / 8.0,
    ];
    // Symmetric index pattern (0-based indices into `k`) — top99 KE layout.
    const P: [[usize; 8]; 8] = [
        [0, 1, 2, 3, 4, 5, 6, 7],
        [1, 0, 7, 6, 5, 4, 3, 2],
        [2, 7, 0, 5, 6, 3, 4, 1],
        [3, 6, 5, 0, 7, 2, 1, 4],
        [4, 5, 6, 7, 0, 1, 2, 3],
        [5, 4, 3, 2, 1, 0, 7, 6],
        [6, 3, 4, 1, 2, 7, 0, 5],
        [7, 2, 1, 4, 3, 6, 5, 0],
    ];
    let coef = 1.0 / (1.0 - nu * nu);
    let mut ke = [[0.0f64; 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            ke[i][j] = coef * k[P[i][j]];
        }
    }
    ke
}

/// 2D SIMP topology optimization of a cantilever beam (left edge clamped, a unit
/// downward point load at the middle of the right edge). Each iteration:
/// assembles the real linear-elasticity stiffness with fs-sparse, solves
/// `K u = f` with a Jacobi-preconditioned CG on `spmv`, computes compliance
/// sensitivities, applies a density (sensitivity) filter, and runs an
/// Optimality-Criteria update — then snapshots the density field.
///
/// Output layout: `iters` density frames concatenated. Each frame is
/// `nx*ny` values in `[0,1]`; element `(ex,ey)` sits at index `ey*nx + ex`
/// (row-major, x fastest, rows `ey` running top → bottom). A grey blob resolves
/// into a truss; total length is `iters * nx * ny`.
///
/// `nx,ny` clamped to `[2,80]`, `iters` to `[1,80]`, `volfrac` to `[0.05,0.95]`.
pub fn topopt_frames(nx_in: usize, ny_in: usize, iters_in: usize, volfrac_in: f64) -> Vec<f64> {
    let nelx = nx_in.clamp(2, 80);
    let nely = ny_in.clamp(2, 80);
    let iters = iters_in.clamp(1, 80);
    let volfrac = volfrac_in.clamp(0.05, 0.95);
    let penal = 3.0f64;
    let rmin = 1.5f64;
    let e_min = 1.0e-9f64;
    let move_lim = 0.2f64;

    let nel = nelx * nely;
    let nnode_y = nely + 1;
    let nd = 2 * (nelx + 1) * (nely + 1);
    let ke = element_stiffness();

    let node_id = |col: usize, row: usize| col * nnode_y + row;
    let edof_of = |ex: usize, ey: usize| -> [usize; 8] {
        let tl = node_id(ex, ey);
        let tr = node_id(ex + 1, ey);
        let br = node_id(ex + 1, ey + 1);
        let bl = node_id(ex, ey + 1);
        [
            2 * tl,
            2 * tl + 1,
            2 * tr,
            2 * tr + 1,
            2 * br,
            2 * br + 1,
            2 * bl,
            2 * bl + 1,
        ]
    };

    // Clamp the entire left edge (col 0).
    let mut is_fixed = vec![false; nd];
    for row in 0..=nely {
        let id = node_id(0, row);
        is_fixed[2 * id] = true;
        is_fixed[2 * id + 1] = true;
    }
    // Free DOF <-> reduced index map.
    let mut reduced_of = vec![usize::MAX; nd];
    let mut free_dofs: Vec<usize> = Vec::new();
    for d in 0..nd {
        if !is_fixed[d] {
            reduced_of[d] = free_dofs.len();
            free_dofs.push(d);
        }
    }
    let nf = free_dofs.len();

    // Unit downward load at the middle of the right edge.
    let load_node = node_id(nelx, nely / 2);
    let load_dof = 2 * load_node + 1;
    let mut f_red = vec![0.0f64; nf];
    if load_dof < nd && !is_fixed[load_dof] {
        f_red[reduced_of[load_dof]] = -1.0;
    }

    // Density (sensitivity) filter weights.
    let rmin_i = rmin.ceil() as isize;
    let mut filt_neigh: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nel];
    let mut filt_hs = vec![0.0f64; nel];
    for ey in 0..nely {
        for ex in 0..nelx {
            let e = ey * nelx + ex;
            let mut hs = 0.0;
            for dy in -rmin_i..=rmin_i {
                for dx in -rmin_i..=rmin_i {
                    let jx = ex as isize + dx;
                    let jy = ey as isize + dy;
                    if jx < 0 || jy < 0 || jx >= nelx as isize || jy >= nely as isize {
                        continue;
                    }
                    let dist = ((dx * dx + dy * dy) as f64).sqrt();
                    let w = rmin - dist;
                    if w > 0.0 {
                        filt_neigh[e].push(((jy as usize) * nelx + (jx as usize), w));
                        hs += w;
                    }
                }
            }
            filt_hs[e] = hs.max(1e-30);
        }
    }

    let mut x = vec![volfrac; nel];
    let mut u = vec![0.0f64; nd];
    let mut out: Vec<f64> = Vec::with_capacity(iters * nel);

    for _iter in 0..iters {
        // ---- Assemble the SIMP stiffness on the free DOFs (fs-sparse). ----
        let mut coo = Coo::new(nf, nf);
        for ey in 0..nely {
            for ex in 0..nelx {
                let e = ey * nelx + ex;
                let estiff = e_min + x[e].powf(penal) * (1.0 - e_min);
                let edof = edof_of(ex, ey);
                for a in 0..8 {
                    let ga = edof[a];
                    if is_fixed[ga] {
                        continue;
                    }
                    let ra = reduced_of[ga];
                    for b in 0..8 {
                        let gb = edof[b];
                        if is_fixed[gb] {
                            continue;
                        }
                        coo.push(ra, reduced_of[gb], estiff * ke[a][b]);
                    }
                }
            }
        }
        let a = coo.assemble();
        let mut inv_diag = vec![1.0f64; nf];
        for i in 0..nf {
            let d = a.get(i, i);
            inv_diag[i] = if d.abs() > 1e-30 { 1.0 / d } else { 1.0 };
        }
        let u_red = pcg(&a, &f_red, &inv_diag, 2000, 1e-7);
        for slot in u.iter_mut() {
            *slot = 0.0;
        }
        for (ri, &gd) in free_dofs.iter().enumerate() {
            u[gd] = u_red[ri];
        }

        // ---- Compliance sensitivities. ----
        let mut dc = vec![0.0f64; nel];
        for ey in 0..nely {
            for ex in 0..nelx {
                let e = ey * nelx + ex;
                let edof = edof_of(ex, ey);
                let mut ue = [0.0f64; 8];
                for a in 0..8 {
                    ue[a] = u[edof[a]];
                }
                let mut ce = 0.0;
                for a in 0..8 {
                    let mut s = 0.0;
                    for b in 0..8 {
                        s += ke[a][b] * ue[b];
                    }
                    ce += ue[a] * s;
                }
                dc[e] = -penal * x[e].powf(penal - 1.0) * (1.0 - e_min) * ce;
            }
        }

        // ---- Sensitivity filter (mesh-independence). ----
        let mut dcn = vec![0.0f64; nel];
        for e in 0..nel {
            let mut acc = 0.0;
            for &(j, w) in &filt_neigh[e] {
                acc += w * x[j] * dc[j];
            }
            dcn[e] = acc / (filt_hs[e] * x[e].max(1e-3));
        }

        // ---- Optimality-Criteria update (bisection on the multiplier). ----
        let mut l1 = 0.0f64;
        let mut l2 = 1.0e9f64;
        let mut xnew = x.clone();
        for _ in 0..80 {
            if l2 - l1 <= 1e-9 * (l1 + l2 + 1e-30) {
                break;
            }
            let lmid = 0.5 * (l1 + l2);
            let mut vol = 0.0;
            for e in 0..nel {
                let be = (-dcn[e] / lmid).max(0.0).sqrt();
                let xe = (x[e] * be)
                    .min(x[e] + move_lim)
                    .min(1.0)
                    .max(x[e] - move_lim)
                    .max(0.001);
                xnew[e] = xe;
                vol += xe;
            }
            if vol > volfrac * nel as f64 {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
        }
        x = xnew;

        out.extend_from_slice(&x);
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  2 · wave2d_frames — spectral 2D wave (fs-fft)                           */
/* ----------------------------------------------------------------------- */

/// Apply a 1D fs-fft transform along every row then every column of an `n×n`
/// row-major complex grid (a separable 2D FFT). `inverse` scales by `1/n²`.
fn fft2d(data: &mut [C64], n: usize, plan: &Fft, inverse: bool) {
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    let mut line = vec![C64::new(0.0, 0.0); n];
    for row in 0..n {
        let base = row * n;
        let slice = &mut data[base..base + n];
        if inverse {
            plan.inverse(slice, &mut scratch);
        } else {
            plan.forward(slice, &mut scratch);
        }
    }
    for col in 0..n {
        for row in 0..n {
            line[row] = data[row * n + col];
        }
        if inverse {
            plan.inverse(&mut line, &mut scratch);
        } else {
            plan.forward(&mut line, &mut scratch);
        }
        for row in 0..n {
            data[row * n + col] = line[row];
        }
    }
}

/// `Δu` on the `n×n` periodic grid, evaluated spectrally: FFT → multiply by
/// `−(kx²+ky²)` → inverse FFT, taking the real part.
fn spectral_laplacian(u: &[f64], n: usize, plan: &Fft, kv: &[f64], cbuf: &mut [C64], out: &mut [f64]) {
    let m = n * n;
    for i in 0..m {
        cbuf[i] = C64::new(u[i], 0.0);
    }
    fft2d(cbuf, n, plan, false);
    for row in 0..n {
        let ky2 = kv[row] * kv[row];
        for col in 0..n {
            let mult = -(kv[col] * kv[col] + ky2);
            let idx = row * n + col;
            cbuf[idx] = C64::new(cbuf[idx].re * mult, cbuf[idx].im * mult);
        }
    }
    fft2d(cbuf, n, plan, true);
    for i in 0..m {
        out[i] = cbuf[i].re;
    }
}

/// A spectral 2D wave solve `u_tt = c²Δu` on a 2π-periodic square: the
/// Laplacian is applied in Fourier space with the real fs-fft transform, and
/// time is advanced by an explicit leapfrog from a centred Gaussian pulse
/// (rest initial velocity). The height field is snapshotted every frame.
///
/// Output layout: `frames` snapshots concatenated, each `n*n` real heights,
/// row-major `y*n + x` (`y` the row). Total length `frames * n * n`.
///
/// `n` is rounded to a power of two in `[8,128]`; `frames` clamped to
/// `[1,180]`; `steps_per_frame` to `[1,8]`.
pub fn wave2d_frames(n_in: usize, frames_in: usize, steps_per_frame_in: usize) -> Vec<f64> {
    let mut n = n_in.clamp(8, 128).next_power_of_two();
    if n > 128 {
        n = 128;
    }
    let frames = frames_in.clamp(1, 180);
    let spf = steps_per_frame_in.clamp(1, 8);
    let m = n * n;
    let plan = Fft::new(n);

    // Integer wavenumbers on the 2π-periodic domain.
    let kv: Vec<f64> = (0..n)
        .map(|j| if j <= n / 2 { j as f64 } else { j as f64 - n as f64 })
        .collect();
    let c = 1.0f64;
    let kmax = (n as f64) / 2.0;
    let dt = 1.4 / (c * kmax.max(1.0)); // leapfrog CFL: dt·c·kmax < 2

    let mut u_prev = vec![0.0f64; m];
    let mut u_cur = vec![0.0f64; m];
    let two_pi = 2.0 * std::f64::consts::PI;
    let sigma = 0.55f64;
    for row in 0..n {
        for col in 0..n {
            let x = two_pi * col as f64 / n as f64 - std::f64::consts::PI;
            let y = two_pi * row as f64 / n as f64 - std::f64::consts::PI;
            let v = det::exp(-(x * x + y * y) / (2.0 * sigma * sigma));
            u_prev[row * n + col] = v;
            u_cur[row * n + col] = v;
        }
    }

    let mut out = Vec::with_capacity(frames * m);
    let mut cbuf = vec![C64::new(0.0, 0.0); m];
    let mut lap = vec![0.0f64; m];
    let dt2c2 = dt * dt * c * c;
    for _f in 0..frames {
        out.extend_from_slice(&u_cur);
        for _ in 0..spf {
            spectral_laplacian(&u_cur, n, &plan, &kv, &mut cbuf, &mut lap);
            for i in 0..m {
                let u_next = 2.0 * u_cur[i] - u_prev[i] + dt2c2 * lap[i];
                u_prev[i] = u_cur[i];
                u_cur[i] = u_next;
            }
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  3 · gray_scott_frames — reaction–diffusion (fs-sparse Laplacian)        */
/* ----------------------------------------------------------------------- */

/// Gray–Scott reaction–diffusion. Diffusion is applied with the real fs-sparse
/// periodic Laplacian via `spmv`; the `u·v²` reaction and feed/kill terms are
/// forward-Euler stepped; the `v` field (the Turing pattern) is snapshotted
/// every frame.
///
/// Output layout: `frames` snapshots concatenated, each `n*n` values of the `v`
/// concentration in `[0,~1]`, row-major `y*n + x`. Total `frames * n * n`.
/// Each frame advances the sim by 14 internal steps.
///
/// `n` clamped to `[16,160]`, `frames` to `[1,180]`, `feed`/`kill` to `[0,0.1]`.
pub fn gray_scott_frames(n_in: usize, frames_in: usize, feed_in: f64, kill_in: f64) -> Vec<f64> {
    let n = n_in.clamp(16, 160);
    let frames = frames_in.clamp(1, 180);
    let feed = feed_in.clamp(0.0, 0.1);
    let kill = kill_in.clamp(0.0, 0.1);
    let m = n * n;
    let lap = laplacian_5pt_periodic(n);
    let du = 0.16f64;
    let dv = 0.08f64;
    let dt = 1.0f64;
    let steps_per_frame = 14usize;

    let mut u = vec![1.0f64; m];
    let mut v = vec![0.0f64; m];
    let mut nz = StreamKey {
        seed: 0x6A57,
        kernel: 0xF5_10,
        tile: 0,
    }
    .stream();
    let c = n / 2;
    let r = (n / 10).max(3) as isize;
    for i in 0..n {
        for j in 0..n {
            let di = i as isize - c as isize;
            let dj = j as isize - c as isize;
            if di.abs() <= r && dj.abs() <= r {
                let jitter = 0.02 * nz.next_f64();
                u[i * n + j] = 0.5 + jitter;
                v[i * n + j] = 0.25 + jitter;
            }
        }
    }

    let mut lu = vec![0.0f64; m];
    let mut lv = vec![0.0f64; m];
    let mut out = Vec::with_capacity(frames * m);
    for _f in 0..frames {
        out.extend_from_slice(&v);
        for _ in 0..steps_per_frame {
            lap.spmv(&u, &mut lu);
            lap.spmv(&v, &mut lv);
            for i in 0..m {
                let uvv = u[i] * v[i] * v[i];
                let un = u[i] + dt * (du * lu[i] - uvv + feed * (1.0 - u[i]));
                let vn = v[i] + dt * (dv * lv[i] + uvv - (feed + kill) * v[i]);
                u[i] = un.clamp(0.0, 1.5);
                v[i] = vn.clamp(0.0, 1.5);
            }
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  4 · fluid_frames — 2D stable fluids (fs-sparse CG pressure solve)       */
/* ----------------------------------------------------------------------- */

/// Semi-Lagrangian bilinear advection of `field` by `(vx,vy)` (cells/step).
fn advect(field: &[f64], vx: &[f64], vy: &[f64], n: usize, dt: f64, out: &mut [f64]) {
    let nf = n as f64;
    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            let mut x = i as f64 - dt * vx[k];
            let mut y = j as f64 - dt * vy[k];
            x = x.clamp(0.0, nf - 1.001);
            y = y.clamp(0.0, nf - 1.001);
            let i0 = x.floor() as usize;
            let j0 = y.floor() as usize;
            let i1 = (i0 + 1).min(n - 1);
            let j1 = (j0 + 1).min(n - 1);
            let sx = x - i0 as f64;
            let sy = y - j0 as f64;
            let a0 = field[j0 * n + i0] * (1.0 - sx) + field[j0 * n + i1] * sx;
            let a1 = field[j1 * n + i0] * (1.0 - sx) + field[j1 * n + i1] * sx;
            out[k] = a0 * (1.0 - sy) + a1 * sy;
        }
    }
}

/// Inject smoke + upward velocity from two wobbling bottom emitters.
fn add_sources(dens: &mut [f64], vy: &mut [f64], n: usize, frame: usize) {
    let t = frame as f64 * 0.15;
    let jrow = ((n as f64 * 0.08) as isize + 1).max(1);
    let rad = (n as f64 * 0.05).max(1.5) as isize;
    for (ci, &cx) in [0.35f64, 0.65].iter().enumerate() {
        let phase = if ci == 0 { t } else { -t };
        let wob = 0.06 * det::sin(phase);
        let cix = ((cx + wob) * n as f64) as isize;
        for dj in -2..=2isize {
            for di in -rad..=rad {
                let ii = cix + di;
                let jj = jrow + dj;
                if ii < 0 || jj < 0 || ii >= n as isize || jj >= n as isize {
                    continue;
                }
                let k = (jj as usize) * n + (ii as usize);
                dens[k] = (dens[k] + 0.6).min(1.0);
                vy[k] += 0.4;
            }
        }
    }
}

/// A 2D stable-fluids smoke simulation: two buoyant smoke sources, semi-
/// Lagrangian advection, and a divergence-free **pressure projection** whose
/// Poisson equation is solved with the real fs-sparse CG (`laplacian_5pt` +
/// `spmv`). The density field is snapshotted every frame.
///
/// Output layout: `frames` snapshots concatenated, each `n*n` density values in
/// `[0,1]`, row-major `j*n + i` with `j=0` at the BOTTOM (smoke rises toward
/// larger `j`). Total length `frames * n * n`.
///
/// `n` clamped to `[16,96]`, `frames` to `[1,200]`.
pub fn fluid_frames(n_in: usize, frames_in: usize) -> Vec<f64> {
    let n = n_in.clamp(16, 96);
    let frames = frames_in.clamp(1, 200);
    let m = n * n;
    // Pressure Poisson operator A = −Δ (SPD, Dirichlet) — real fs-sparse.
    let a = laplacian_5pt(n);

    let dt = 1.0f64;
    let buoyancy = 0.12f64;
    let vmax = n as f64 / 4.0;

    let mut vx = vec![0.0f64; m];
    let mut vy = vec![0.0f64; m];
    let mut dens = vec![0.0f64; m];
    let mut div = vec![0.0f64; m];
    let mut out = Vec::with_capacity(frames * m);

    for frame in 0..frames {
        add_sources(&mut dens, &mut vy, n, frame);
        for k in 0..m {
            vy[k] += dt * buoyancy * dens[k];
        }
        // Advect velocity.
        let vx0 = vx.clone();
        let vy0 = vy.clone();
        advect(&vx0, &vx0, &vy0, n, dt, &mut vx);
        advect(&vy0, &vx0, &vy0, n, dt, &mut vy);
        // Divergence.
        for j in 0..n {
            for i in 0..n {
                let ip = (i + 1).min(n - 1);
                let im = i.saturating_sub(1);
                let jp = (j + 1).min(n - 1);
                let jm = j.saturating_sub(1);
                div[j * n + i] = 0.5 * ((vx[j * n + ip] - vx[j * n + im]) + (vy[jp * n + i] - vy[jm * n + i]));
            }
        }
        // Pressure: solve A p = −div (real fs-sparse CG).
        let rhs: Vec<f64> = div.iter().map(|d| -d).collect();
        let p = cg(&a, &rhs, 400, 1e-5);
        // Subtract the pressure gradient → divergence-free velocity.
        for j in 0..n {
            for i in 0..n {
                let ip = (i + 1).min(n - 1);
                let im = i.saturating_sub(1);
                let jp = (j + 1).min(n - 1);
                let jm = j.saturating_sub(1);
                vx[j * n + i] -= 0.5 * (p[j * n + ip] - p[j * n + im]);
                vy[j * n + i] -= 0.5 * (p[jp * n + i] - p[jm * n + i]);
            }
        }
        // Advect density, then dissipate + clamp for a bounded, safe field.
        let dens0 = dens.clone();
        advect(&dens0, &vx, &vy, n, dt, &mut dens);
        for k in 0..m {
            vx[k] = (vx[k] * 0.999).clamp(-vmax, vmax);
            vy[k] = (vy[k] * 0.999).clamp(-vmax, vmax);
            dens[k] = (dens[k] * 0.994).clamp(0.0, 1.0);
        }
        out.extend_from_slice(&dens);
    }
    out
}
