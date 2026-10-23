//! Reflection at the open end of a thin unflanged pipe, on the dense open P1
//! screen owner ([`super::open_screen`]): Levine and Schwinger's own geometry,
//! a zero-thickness tube.
//!
//! The pipe is a "cup": a zero-thickness cylinder of radius `a` from `z = 0`
//! to the mouth at `z = L`, closed at `z = 0` by a zero-thickness hemisphere.
//! It is one oriented open two-manifold whose only free edge is the mouth rim.
//! The cylinder's rings sit at `z_k = L sin(pi k / 2K)`, which grades the
//! distance to the rim quadratically, as the disc test data does at its rim.
//!
//! Source: two axial monopoles, at `z_s` and `z_s + delta`, with strengths 1
//! and `-exp(-i k delta)`. Their backward plane waves cancel, so the forward
//! plane-wave amplitude never vanishes. A single monopole's would vanish
//! wherever the cap's reflection cancels its direct wave. The screen is
//! rigid: its drive velocity is minus the incident normal velocity, taken
//! with the smooth surface's normal at each vertex.
//!
//! Nothing radiates outside except through the mouth. So everything behind
//! the stations (cap, sources and their discretization) changes only the
//! incident amplitude `A`, never `R = B/A`.
//!
//! Receivers: the total pressure (incident plus the screen's double layer) on
//! a ring of declared radius, averaged over four azimuths.
//! - At the node ring `r = (2.4048/3.8317) a` the first radial mode
//!   `J0(3.8317 r/a)` vanishes.
//! - On the axis that mode is barely evanescent near ka 3.8 (it decays as
//!   `exp(-0.49 z/a)`), so it would bias any plane-wave fit there.
//! - The next mode decays as `exp(-5.9 z/a)`.
//!
//! Fit (`exp(-i omega t)`): `p(z) = A e^{ikz} + B e^{-ikz}` by least squares
//! over the stations. This gives `R = (B/A) e^{-2ikL}` at the rim plane, and
//! `l/a = arg(-R) / (2ka)` with `arg` in `[-pi/2, 3pi/2)`.
//!
//! The second estimate is the zero crossing of the linear interior profile:
//! `Re(-alpha/beta) - L` for `p ~ alpha + beta z`. It is the low-frequency
//! definition of the end correction, and its bias grows as `(k span)^2`.
//!
//! Every row carries the fit's residual and condition number, and `k` times
//! the station span. At low ka the two-wave fit is ill-conditioned by
//! construction.
use super::open_screen::{OpenOptions, OpenPanels, PreparedScreen};
use super::*;
use crate::helmholtz::bc::{Admittance, SurfaceBc};

/// The node of `J0(3.8317 r / a)`, in units of the radius: `j_{0,1} / j'_{0,1}`.
pub const NODE_RING: f64 = 2.404_825_557_695_773 / 3.831_705_970_207_512;

/// A declared cup, source pair and station window. Lengths in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PipeEndSpec {
    /// Tube radius `a`.
    pub radius: f64,
    /// Cap equator (`z = 0`) to the mouth.
    pub length: f64,
    /// Vertices per cylinder ring (even, at least 8).
    pub circumferential: usize,
    /// Cylinder ring intervals `K`.
    pub axial: usize,
    /// `z` of the first monopole, on the axis inside the cylinder.
    pub source: f64,
    /// Spacing `delta` of the second monopole above the first.
    pub source_spacing: f64,
    /// Station window `[z0, z1]`, above the sources and below the mouth.
    pub window: (f64, f64),
    /// Stations, equally spaced over the window inclusive (at least 3).
    pub stations: usize,
    /// Receiver ring radius; `NODE_RING * radius` cancels the first radial mode.
    pub ring_radius: f64,
}

/// One frequency of one cup.
#[derive(Debug, Clone, PartialEq)]
pub struct PipeEndRow {
    /// `k a`.
    pub ka: f64,
    /// Pressure reflection at the rim plane.
    pub reflection: C64,
    /// `l / a` from the two-wave fit.
    pub end_correction: f64,
    /// `l / a` from the zero crossing of the linear interior profile.
    pub zero_crossing_end_correction: f64,
    /// Relative residual of the two-wave fit.
    pub fit_residual: f64,
    /// Condition number of the two-wave design matrix.
    pub fit_condition: f64,
    /// `k` times the station span.
    pub k_span: f64,
    /// P1 unknowns of the screen.
    pub unknowns: usize,
    /// The declared edge level of the mesh.
    pub edge_level: u32,
    /// Backward residual of the screen's LU solve.
    pub backward_residual: f64,
    /// Station `z` [m].
    pub stations: Vec<f64>,
    /// Ring-averaged total pressure at each station [Pa per unit source].
    pub pressures: Vec<C64>,
}

/// A least-squares standing wave `A e^{ikz} + B e^{-ikz}` referred to a mouth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StandingWaveFit {
    /// `R = (B/A) e^{-2ikL}`.
    pub reflection: C64,
    /// `arg(-R) / (2 |k| a)`, with `arg` in `[-pi/2, 3pi/2)`.
    pub end_correction: f64,
    /// `|p - fit| / |p|`.
    pub residual: f64,
    /// Condition number of the design matrix `[e^{ikz}, e^{-ikz}]`.
    pub condition: f64,
}

fn cis(phase: f64) -> C64 {
    C64::new(det::cos(phase), det::sin(phase))
}

/// Fit `p(z_j) = A e^{ikz_j} + B e^{-ikz_j}` and refer `R = B/A` to `mouth`.
/// A negative `k` fits under the conjugate time convention.
///
/// # Errors
/// Fewer than three stations, mismatched lengths, a singular design or a
/// non-finite input.
pub fn standing_wave_fit(
    k: f64,
    radius: f64,
    mouth: f64,
    z: &[f64],
    p: &[C64],
) -> Result<StandingWaveFit, HelmholtzError> {
    if z.len() < 3 || z.len() != p.len() || !k.is_finite() || k == 0.0 || !(radius > 0.0) {
        return Err(bad(
            "standing-wave fit needs three stations and a nonzero k",
        ));
    }
    let n = z.len() as f64;
    // G = M^H M = [[n, s], [conj(s), n]] with s = sum e^{-2ikz}.
    let mut s = C64::ZERO;
    let (mut up, mut vp) = (C64::ZERO, C64::ZERO);
    for (&zj, &pj) in z.iter().zip(p) {
        let u = cis(k * zj);
        s = s + cis(-2.0 * k * zj);
        up = up + u.conj() * pj;
        vp = vp + u * pj;
    }
    let det_g = n * n - s.norm_sq();
    if !(det_g > 0.0) {
        return Err(bad("singular standing-wave design"));
    }
    let a = (up.scale(n) - s * vp).scale(1.0 / det_g);
    let b = (vp.scale(n) - s.conj() * up).scale(1.0 / det_g);
    let (mut defect, mut size) = (0.0, 0.0);
    for (&zj, &pj) in z.iter().zip(p) {
        let u = cis(k * zj);
        defect += (pj - a * u - b * u.conj()).norm_sq();
        size += pj.norm_sq();
    }
    let reflection = b / a * cis(-2.0 * k * mouth);
    let minus = -reflection;
    let mut arg = det::atan2(minus.im, minus.re);
    if arg < -0.5 * PI {
        arg += TAU;
    }
    // |k|: a conjugate-convention fit (negative k) reports what an analyst
    // would read from its R = 1/R_true, the end correction's negative.
    let result = StandingWaveFit {
        reflection,
        end_correction: arg / (2.0 * k.abs() * radius),
        residual: det::sqrt(defect / size),
        condition: det::sqrt((n + s.abs()) / (n - s.abs())),
    };
    if !finite(result.reflection) || !result.end_correction.is_finite() {
        return Err(bad("non-finite standing-wave fit"));
    }
    Ok(result)
}

/// `l / a` from the zero crossing `Re(-alpha/beta) - mouth` of the least-squares
/// line `p ~ alpha + beta z`.
///
/// # Errors
/// Fewer than two distinct stations, mismatched lengths or a flat profile.
pub fn zero_crossing_end_correction(
    radius: f64,
    mouth: f64,
    z: &[f64],
    p: &[C64],
) -> Result<f64, HelmholtzError> {
    if z.len() < 2 || z.len() != p.len() || !(radius > 0.0) {
        return Err(bad("zero crossing needs two stations"));
    }
    let n = z.len() as f64;
    let zm = z.iter().sum::<f64>() / n;
    let pm = p.iter().fold(C64::ZERO, |s, &v| s + v).scale(1.0 / n);
    let (mut szz, mut szp) = (0.0, C64::ZERO);
    for (&zj, &pj) in z.iter().zip(p) {
        szz += (zj - zm) * (zj - zm);
        szp = szp + pj.scale(zj - zm);
    }
    if !(szz > 0.0) || szp.abs() == 0.0 {
        return Err(bad("flat interior profile"));
    }
    let beta = szp.scale(1.0 / szz);
    let crossing = C64::from_re(zm) - pm / beta;
    Ok((crossing.re - mouth) / radius)
}

fn validate(spec: &PipeEndSpec) -> Result<(), HelmholtzError> {
    let a = spec.radius;
    let lengths = [
        a,
        spec.length,
        spec.source,
        spec.source_spacing,
        spec.window.0,
        spec.window.1,
        spec.ring_radius,
    ];
    if lengths.iter().any(|v| !v.is_finite())
        || !(a > 0.0)
        || spec.circumferential < 8
        || spec.circumferential % 2 != 0
        || spec.axial < 2
        || spec.stations < 3
        || !(spec.source > 0.0)
        || !(spec.source_spacing > 0.0)
        || !(spec.source + spec.source_spacing < spec.window.0)
        || !(spec.window.0 < spec.window.1)
        || !(spec.window.1 < spec.length)
        || !(spec.ring_radius >= 0.0 && spec.ring_radius < a)
    {
        return Err(bad(
            "pipe end needs a > 0, an even circumferential count >= 8, axial >= 2, three stations, \
             0 < source < source + spacing < window < length, and a ring inside the tube",
        ));
    }
    Ok(())
}

/// The cup's mesh: cylinder rings at `z_k = L sin(pi k / 2K)`, then a
/// hemisphere whose rings halve their vertex count (while even and at least
/// 8) as their circumference shrinks, closed by a pole fan. Normals point out
/// of the cup.
///
/// # Errors
/// An invalid spec, or a mesh `OpenPanels` refuses.
pub fn pipe_cup(spec: &PipeEndSpec) -> Result<OpenPanels, HelmholtzError> {
    validate(spec)?;
    let (a, n, rings) = (spec.radius, spec.circumferential, spec.axial);
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles = Vec::new();
    for k in 0..=rings {
        let z = spec.length * det::sin(PI * k as f64 / (2.0 * rings as f64));
        for i in 0..n {
            let theta = TAU * i as f64 / n as f64;
            vertices.push([a * det::cos(theta), a * det::sin(theta), z]);
        }
    }
    let at = |k: usize, i: usize| k * n + i % n;
    for k in 0..rings {
        for i in 0..n {
            triangles.push([at(k, i), at(k, i + 1), at(k + 1, i + 1)]);
            triangles.push([at(k, i), at(k + 1, i + 1), at(k + 1, i)]);
        }
    }
    // Hemisphere bands: ring `upper` (toward the equator) over ring `lower`.
    let bands = (n / 4).max(2);
    let mut upper: Vec<usize> = (0..n).map(|i| at(0, i)).collect();
    for j in 1..bands {
        let psi = 0.5 * PI * j as f64 / bands as f64;
        let previous = upper.len();
        let count = if previous % 2 == 0
            && previous >= 8
            && n as f64 * det::cos(psi) <= 0.5 * previous as f64
        {
            previous / 2
        } else {
            previous
        };
        let start = vertices.len();
        for i in 0..count {
            let theta = TAU * i as f64 / count as f64;
            let rho = a * det::cos(psi);
            vertices.push([
                rho * det::cos(theta),
                rho * det::sin(theta),
                -a * det::sin(psi),
            ]);
        }
        let lower: Vec<usize> = (start..start + count).collect();
        let (nu, nl) = (upper.len(), lower.len());
        if nu == nl {
            for i in 0..nl {
                triangles.push([lower[i], lower[(i + 1) % nl], upper[(i + 1) % nu]]);
                triangles.push([lower[i], upper[(i + 1) % nu], upper[i]]);
            }
        } else {
            for i in 0..nl {
                triangles.push([lower[i], upper[2 * i + 1], upper[2 * i]]);
                triangles.push([lower[i], lower[(i + 1) % nl], upper[2 * i + 1]]);
                triangles.push([
                    lower[(i + 1) % nl],
                    upper[(2 * i + 2) % nu],
                    upper[2 * i + 1],
                ]);
            }
        }
        upper = lower;
    }
    let pole = vertices.len();
    vertices.push([0.0, 0.0, -a]);
    let nu = upper.len();
    for i in 0..nu {
        triangles.push([pole, upper[(i + 1) % nu], upper[i]]);
    }
    OpenPanels::new(vertices, triangles)
}

/// The outward normal of the smooth cup at a mesh vertex.
fn cup_normal(x: [f64; 3]) -> [f64; 3] {
    let r = if x[2] >= 0.0 { [x[0], x[1], 0.0] } else { x };
    let length = norm(r);
    r.map(|v| v / length)
}

/// The source pair's pressure and gradient at `x`.
fn incident(k: f64, spec: &PipeEndSpec, x: [f64; 3]) -> (C64, [C64; 3]) {
    let sources = [
        (spec.source, C64::ONE),
        (
            spec.source + spec.source_spacing,
            -cis(-k * spec.source_spacing),
        ),
    ];
    let mut p = C64::ZERO;
    let mut g = [C64::ZERO; 3];
    for (zs, strength) in sources {
        let d = [x[0], x[1], x[2] - zs];
        let r = norm(d);
        let green = strength * cis(k * r).scale(1.0 / (4.0 * PI * r));
        p = p + green;
        // grad(e^{ikr}/(4 pi r)) = G (ik - 1/r) r_hat.
        let radial = green * C64::new(-1.0 / r, k);
        for (gc, dc) in g.iter_mut().zip(d) {
            *gc = *gc + radial.scale(dc / r);
        }
    }
    (p, g)
}

/// Solve the cup at every `ka` and fit each standing wave.
///
/// # Errors
/// An invalid spec or `ka`, a refusal of the open-screen owner (dense cap,
/// coarse mesh, budgets, solve residual), a receiver refusal, or a degenerate fit.
pub fn pipe_end_reflection(
    spec: &PipeEndSpec,
    ka: &[f64],
    options: OpenOptions,
    medium: Medium,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<PipeEndRow>, HelmholtzError> {
    validate(spec)?;
    if ka.is_empty() || ka.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(bad("pipe end needs finite positive ka"));
    }
    let a = spec.radius;
    let highest = ka.iter().fold(0.0_f64, |m, v| m.max(*v)) / a;
    let surface = pipe_cup(spec)?;
    let prepared = PreparedScreen::new(&surface, highest, options, &mut cancelled)?;
    let stations: Vec<f64> = (0..spec.stations)
        .map(|j| {
            spec.window.0 + (spec.window.1 - spec.window.0) * j as f64 / (spec.stations - 1) as f64
        })
        .collect();
    let azimuths = 4;
    let points: Vec<[f64; 3]> = stations
        .iter()
        .flat_map(|&z| {
            (0..azimuths).map(move |m| {
                let theta = TAU * m as f64 / azimuths as f64;
                [
                    spec.ring_radius * det::cos(theta),
                    spec.ring_radius * det::sin(theta),
                    z,
                ]
            })
        })
        .collect();
    let receivers = crate::near_field::Options {
        relative_tolerance: 1e-9,
        maximum_depth: 14,
        maximum_kernel_evaluations: 20_000_000,
    };
    let mut rows = Vec::with_capacity(ka.len());
    for &kaj in ka {
        let k = kaj / a;
        let omega = k * medium.sound_speed;
        let drive: Vec<C64> = surface
            .vertices()
            .iter()
            .map(|&x| {
                let (_, g) = incident(k, spec, x);
                let n = cup_normal(x);
                let normal = g[0].scale(n[0]) + g[1].scale(n[1]) + g[2].scale(n[2]);
                // v = grad p / (i omega rho); a rigid screen cancels it.
                -(normal / C64::new(0.0, omega * medium.density))
            })
            .collect();
        let solution = prepared.solve(
            k,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Rigid,
            },
            &mut cancelled,
        )?;
        // One receiver request per station: each stays inside the
        // near-field owner's per-request evaluation ceiling at every rung.
        let mut scattered = Vec::with_capacity(points.len());
        for station in points.chunks(azimuths) {
            scattered.extend(solution.pressure_at(station, 0.05 * a, receivers, &mut cancelled)?);
        }
        let pressures: Vec<C64> = stations
            .iter()
            .enumerate()
            .map(|(j, _)| {
                (0..azimuths)
                    .fold(C64::ZERO, |sum, m| {
                        let index = j * azimuths + m;
                        sum + incident(k, spec, points[index]).0 + scattered[index]
                    })
                    .scale(1.0 / azimuths as f64)
            })
            .collect();
        let fit = standing_wave_fit(k, a, spec.length, &stations, &pressures)?;
        rows.push(PipeEndRow {
            ka: kaj,
            reflection: fit.reflection,
            end_correction: fit.end_correction,
            zero_crossing_end_correction: zero_crossing_end_correction(
                a,
                spec.length,
                &stations,
                &pressures,
            )?,
            fit_residual: fit.residual,
            fit_condition: fit.condition,
            k_span: k * (spec.window.1 - spec.window.0),
            unknowns: prepared.unknowns(),
            edge_level: options.edge_level,
            backward_residual: solution.backward_residual,
            stations: stations.clone(),
            pressures,
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standing_wave_fit_recovers_a_planted_reflection_and_flips_under_conjugation() {
        let (a, mouth, k) = (1.0, 4.6, 2.0);
        let planted = C64::new(-0.8, -0.3);
        let z: Vec<f64> = (0..12).map(|j| 2.0 + 0.1 * j as f64).collect();
        let amplitude = C64::new(0.7, 0.2);
        let p: Vec<C64> = z
            .iter()
            .map(|&zj| amplitude * (cis(k * zj) + planted * cis(2.0 * k * mouth) * cis(-k * zj)))
            .collect();
        let fit = standing_wave_fit(k, a, mouth, &z, &p).unwrap();
        assert!(
            (fit.reflection - planted).abs() < 1e-13,
            "{:?}",
            fit.reflection
        );
        assert!(fit.residual < 1e-13);
        let minus = -planted;
        let l = det::atan2(minus.im, minus.re) / (2.0 * k * a);
        assert!((fit.end_correction - l).abs() < 1e-13);
        // The conjugate convention swaps the waves: R -> 1/R, l -> -l.
        let flipped = standing_wave_fit(-k, a, mouth, &z, &p).unwrap();
        assert!((flipped.reflection * planted - C64::ONE).abs() < 1e-12);
        assert!((flipped.end_correction + l).abs() < 1e-12);
    }

    #[test]
    fn zero_crossing_of_a_line() {
        let z = [1.0, 1.5, 2.0, 2.5];
        let p: Vec<C64> = z
            .iter()
            .map(|&zj| C64::new(0.0, 2.0).scale(zj - 5.25))
            .collect();
        let l = zero_crossing_end_correction(0.5, 5.0, &z, &p).unwrap();
        assert!((l - 0.5).abs() < 1e-14, "{l}");
    }

    fn spec() -> PipeEndSpec {
        PipeEndSpec {
            radius: 1.0,
            length: 4.6,
            circumferential: 16,
            axial: 12,
            source: 0.3,
            source_spacing: 0.4,
            window: (2.0, 3.1),
            stations: 8,
            ring_radius: NODE_RING,
        }
    }

    #[test]
    fn the_cup_is_one_open_manifold_with_outward_normals() {
        let s = spec();
        let cup = pipe_cup(&s).unwrap();
        let rim = (0..cup.vertices().len())
            .filter(|&i| cup.is_free_edge_vertex(i))
            .count();
        assert_eq!(rim, s.circumferential, "only the mouth ring is free");
        for t in cup.triangles() {
            let [p, q, r] = t.map(|i| cup.vertices()[i]);
            let normal = cross(sub(q, p), sub(r, p));
            let centre = [
                (p[0] + q[0] + r[0]) / 3.0,
                (p[1] + q[1] + r[1]) / 3.0,
                (p[2] + q[2] + r[2]) / 3.0,
            ];
            assert!(
                dot(normal, cup_normal(centre)) > 0.0,
                "inward facet at {centre:?}"
            );
        }
        let mut bad_spec = s;
        bad_spec.circumferential = 15;
        assert!(pipe_cup(&bad_spec).is_err());
        bad_spec = s;
        bad_spec.window = (0.5, 3.1);
        assert!(pipe_cup(&bad_spec).is_err(), "window over the sources");
    }

    #[test]
    fn the_source_pair_sends_no_backward_plane_wave() {
        // In free space the pair's far field along -z cancels: its pressure
        // falls as 1/r^2 there, against 1/r along +z.
        let s = spec();
        let k = 1.7;
        let far = |z: f64| incident(k, &s, [0.0, 0.0, z]).0.abs() * (z - s.source).abs();
        assert!(far(-1e4) < 1e-3 * far(1e4), "{} {}", far(-1e4), far(1e4));
    }
}
