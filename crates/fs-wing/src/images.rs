//! FlatPlaneVortexImageExact (bead wf-root-guzez.5.10, E4.4a). Plan
//! §5.2/§4: ground effect enters through an EXACT vortex-image system
//! across the certified flat plane — every horseshoe gains a mirror
//! (corners reflected in z, a/b SWAPPED so the image circulation opposes,
//! trailing legs along the MIRRORED stream), which enforces zero normal
//! flow through the plane analytically; the battery measures the residual
//! at the plane at machine precision (V-06a exactness clause).
//!
//! CERTIFICATE GATE: the image model is only admissible over a plane
//! carrying FlatnessCertificate numbers inside the declared bands (the
//! same bands fs-flyer::prelaunch issues under; cross-linked constants).
//! An uncertified plane is a typed refusal — the image model may never
//! silently pretend a dune is flat.

use crate::{MAX_PANELS, Panel, Refusal, SolveReport, SurfaceId};

/// Certificate bands (MUST match fs-flyer::prelaunch's issuance bands;
/// the battery cross-checks the values).
pub const CERT_MAX_SLOPE: f64 = 0.005;
/// RMS band [m] (same source).
pub const CERT_MAX_RMS_M: f64 = 1.2;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// The certified ground plane the image system reflects across.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CertifiedGround {
    /// Plane height z (frd/NED: positive DOWN; aircraft z < z_plane).
    pub z_m: f64,
    /// The FlatnessCertificate's fitted slope magnitude.
    pub certificate_slope: f64,
    /// The certificate's RMS residual [m].
    pub certificate_rms_m: f64,
}

impl CertifiedGround {
    /// Admit the plane: certificate numbers inside the bands.
    ///
    /// # Errors
    /// `ground-uncertified` naming the violated band (slope or RMS,
    /// tested at band AND the next float).
    pub fn admit(&self) -> Result<(), Refusal> {
        if !self.z_m.is_finite()
            || !self.certificate_slope.is_finite()
            || !self.certificate_rms_m.is_finite()
        {
            return Err(refuse(
                "ground-uncertified",
                "non-finite".into(),
                "real certificate",
            ));
        }
        if self.certificate_slope > CERT_MAX_SLOPE {
            return Err(refuse(
                "ground-uncertified",
                format!("slope {} beyond {CERT_MAX_SLOPE}", self.certificate_slope),
                "the flat-plane image model may not be used over uncertified terrain",
            ));
        }
        if self.certificate_rms_m > CERT_MAX_RMS_M {
            return Err(refuse(
                "ground-uncertified",
                format!("RMS {} m beyond {CERT_MAX_RMS_M}", self.certificate_rms_m),
                "dune texture too strong for the certified plane",
            ));
        }
        Ok(())
    }
}

fn mirror_point(p: [f64; 3], zg: f64) -> [f64; 3] {
    [p[0], p[1], 2.0 * zg - p[2]]
}

/// Combined induced velocity of horseshoe j AND its exact image at `p`.
pub(crate) fn horseshoe_with_image(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    stream: [f64; 3],
    zg: f64,
) -> [f64; 3] {
    let v = crate::horseshoe_velocity_pub(p, a, b, stream);
    // Image: mirrored corners, a/b swapped (opposite circulation),
    // trailing legs along the MIRRORED stream.
    let ms = [stream[0], stream[1], -stream[2]];
    let vi = crate::horseshoe_velocity_pub(p, mirror_point(b, zg), mirror_point(a, zg), ms);
    [v[0] + vi[0], v[1] + vi[1], v[2] + vi[2]]
}

/// The ground-effect Weissinger solve: identical to the free-air solve
/// except every influence entry carries the exact image term. The
/// free-air function is untouched (its pinned goldens stand).
///
/// # Errors
/// Free-air refusal set plus `ground-uncertified` and
/// `aircraft-below-ground`.
pub fn solve_weissinger_ground(
    panels: &[Panel],
    freestream_mps: [f64; 3],
    rho_kg_m3: f64,
    ground: &CertifiedGround,
) -> Result<SolveReport, Refusal> {
    ground.admit()?;
    let n = panels.len();
    if n == 0 || n > MAX_PANELS {
        return Err(refuse("panel-count-invalid", format!("{n}"), "cap"));
    }
    for p in panels {
        // frd z-down: the aircraft must be ABOVE the plane (z < z_ground).
        if p.ctrl[2] >= ground.z_m || p.a[2] >= ground.z_m || p.b[2] >= ground.z_m {
            return Err(refuse(
                "aircraft-below-ground",
                format!("{:?} panel at/below the certified plane", p.surface),
                "contact and crash regimes belong to fs-flyer, not the image model",
            ));
        }
    }
    let vmag = (freestream_mps.iter().map(|v| v * v).sum::<f64>()).sqrt();
    if !freestream_mps.iter().all(|v| v.is_finite()) || vmag < 1.0e-6 || rho_kg_m3 <= 0.0 {
        return Err(refuse(
            "freestream-invalid",
            format!("{freestream_mps:?}"),
            "finite",
        ));
    }
    let stream = [
        freestream_mps[0] / vmag,
        freestream_mps[1] / vmag,
        freestream_mps[2] / vmag,
    ];
    let dotn = |v: [f64; 3], w: [f64; 3]| v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    let mut a = vec![0.0f64; n * n];
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            let v =
                horseshoe_with_image(panels[i].ctrl, panels[j].a, panels[j].b, stream, ground.z_m);
            a[i * n + j] = dotn(v, panels[i].normal);
        }
        rhs[i] = -dotn(freestream_mps, panels[i].normal);
    }
    // Deterministic LU (same fixed tie rule as the free-air path).
    let mut lu = a.clone();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut piv = k;
        let mut best = lu[perm[k] * n + k].abs();
        for r in (k + 1)..n {
            let m = lu[perm[r] * n + k].abs();
            if m > best {
                best = m;
                piv = r;
            }
        }
        if best == 0.0 {
            return Err(refuse(
                "influence-singular",
                format!("column {k}"),
                "degenerate geometry",
            ));
        }
        perm.swap(k, piv);
        let pk = perm[k];
        for r in (k + 1)..n {
            let pr = perm[r];
            let f = lu[pr * n + k] / lu[pk * n + k];
            lu[pr * n + k] = f;
            for c in (k + 1)..n {
                lu[pr * n + c] -= f * lu[pk * n + c];
            }
        }
    }
    let mut gamma = vec![0.0f64; n];
    {
        let mut y = vec![0.0f64; n];
        for r in 0..n {
            let mut s = rhs[perm[r]];
            for c in 0..r {
                s -= lu[perm[r] * n + c] * y[c];
            }
            y[r] = s;
        }
        for r in (0..n).rev() {
            let mut s = y[r];
            for c in (r + 1)..n {
                s -= lu[perm[r] * n + c] * gamma[c];
            }
            gamma[r] = s / lu[perm[r] * n + r];
        }
    }
    // Forces (identical Kutta–Joukowsky reduction to the free-air path).
    let mut per: Vec<(SurfaceId, f64)> = Vec::new();
    let mut total = 0.0;
    for (j, p) in panels.iter().enumerate() {
        let seg = [p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]];
        let f = [
            freestream_mps[1] * seg[2] - freestream_mps[2] * seg[1],
            freestream_mps[2] * seg[0] - freestream_mps[0] * seg[2],
            freestream_mps[0] * seg[1] - freestream_mps[1] * seg[0],
        ];
        let lift = -rho_kg_m3 * gamma[j] * f[2];
        total += lift;
        match per.iter_mut().find(|(s, _)| *s == p.surface) {
            Some((_, acc)) => *acc += lift,
            None => per.push((p.surface, lift)),
        }
    }
    Ok(SolveReport {
        gamma,
        condition_est: 0.0,
        surface_lift_n: per,
        total_lift_n: total,
    })
}

/// Total induced velocity of the SOLVED system (real + images) at any
/// point — the battery's plane-residual probe.
#[must_use]
pub fn induced_velocity_with_images(
    p: [f64; 3],
    panels: &[Panel],
    gamma: &[f64],
    freestream_mps: [f64; 3],
    zg: f64,
) -> [f64; 3] {
    let vmag = (freestream_mps.iter().map(|v| v * v).sum::<f64>()).sqrt();
    let stream = [
        freestream_mps[0] / vmag,
        freestream_mps[1] / vmag,
        freestream_mps[2] / vmag,
    ];
    let mut out = [0.0f64; 3];
    for (j, panel) in panels.iter().enumerate() {
        let v = horseshoe_with_image(p, panel.a, panel.b, stream, zg);
        for k in 0..3 {
            out[k] += gamma[j] * v[k];
        }
    }
    out
}
