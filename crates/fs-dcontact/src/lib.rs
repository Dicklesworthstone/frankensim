//! # fs-dcontact — distributed unilateral contact
//!
//! Distributed collisions along lines and profiles — string-fretboard
//! rattle, sitar/tanpura bridge (jawari), snare wires, reed lays — as
//! ONE-SIDED power-law penalty potentials over collocation points:
//!
//! `Phi_contact = sum_i w_i * K/(alpha+1) * [p_i]_+^(alpha+1)`,
//! `p_i = (Phi q)_i - c_i`  (penetration of the displacement field
//! past the obstacle profile),
//!
//! added to any fs-phs `Storage` ([`ContactStorage`]). Because the
//! contact enters as a POTENTIAL with its exact analytic gradient, the
//! Gonzalez discrete-gradient stepper makes collisions energy-exact by
//! construction (the Bilbao-Chatziioannou energy-consistent collision
//! doctrine, realized through the pHS machinery instead of a bespoke
//! scheme) — no LCP, provably stable at rattle stiffnesses. The
//! implicit contact Newton, its iteration budget, and the achieved
//! residual are fs-phs's disclosed caps (`StepRecord::newton_iters`,
//! `solver_residual`) — one frozen solver law, not a second one.
//!
//! Honest scope (stated): NORMAL contact only — tangential friction
//! (bowing) is its own future bead; contact-internal viscous loss
//! (Hunt-Crossley-style damping inside the collision) is a recorded
//! follow-up, so restitution of the bare potential is exactly 1 and
//! losses enter through the modal damping `R`; contact-law parameters
//! `(K, alpha)` are caller-supplied with a PROVENANCE string logged —
//! a matdb lookup is deferred until packs carrying contact-law
//! parameters exist (no fake wiring).

use fs_math::det;
use fs_phs::Storage;

/// A distributed obstacle: collocation points with per-point gap,
/// quadrature weight, and the shared power-law contact law. Fields
/// are PRIVATE so admission is mandatory (review finding: pub fields
/// made every check advisory); mutation batteries use the documented
/// [`Obstacle::from_raw_parts`] trust escape.
#[derive(Debug, Clone)]
pub struct Obstacle {
    collocation: Vec<f64>,
    n_points: usize,
    gaps: Vec<f64>,
    weights: Vec<f64>,
    stiffness: f64,
    alpha: f64,
    provenance: String,
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum DContactError {
    /// Shape/length mismatch.
    Shape {
        /// What disagreed.
        what: &'static str,
    },
    /// Non-physical contact-law parameter.
    Parameter {
        /// Which one.
        what: &'static str,
    },
}

impl core::fmt::Display for DContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DContactError::Shape { what } => write!(f, "shape mismatch: {what}"),
            DContactError::Parameter { what } => write!(f, "bad parameter: {what}"),
        }
    }
}

impl std::error::Error for DContactError {}

impl Obstacle {
    /// Admit an obstacle: consistent lengths, finite entries,
    /// `K >= 0`, `alpha >= 1`. The load-bearing property is that the
    /// contact POTENTIAL stays C^1 for alpha >= 1 (what the Gonzalez
    /// scheme needs); the force itself is C^1 only for alpha > 1, and
    /// its Hessian is unbounded for 1 < alpha < 2 (stiff-graze
    /// guidance: alpha >= 2).
    ///
    /// # Errors
    /// Typed [`DContactError`].
    #[allow(clippy::too_many_arguments)] // one admission door for one flat data record
    pub fn new(
        collocation: Vec<f64>,
        n_points: usize,
        n_modes: usize,
        gaps: Vec<f64>,
        weights: Vec<f64>,
        stiffness: f64,
        alpha: f64,
        provenance: String,
    ) -> Result<Self, DContactError> {
        if collocation.len() != n_points * n_modes
            || gaps.len() != n_points
            || weights.len() != n_points
        {
            return Err(DContactError::Shape {
                what: "collocation/gaps/weights vs (n_points, n_modes)",
            });
        }
        if collocation
            .iter()
            .chain(&gaps)
            .chain(&weights)
            .any(|v| !v.is_finite())
        {
            return Err(DContactError::Parameter {
                what: "non-finite obstacle entry",
            });
        }
        if weights.iter().any(|&w| w < 0.0) {
            return Err(DContactError::Parameter {
                what: "negative quadrature weight",
            });
        }
        if !stiffness.is_finite() || stiffness < 0.0 {
            return Err(DContactError::Parameter {
                what: "contact stiffness must be non-negative and finite",
            });
        }
        if !alpha.is_finite() || alpha < 1.0 {
            return Err(DContactError::Parameter {
                what: "contact exponent must be at least 1",
            });
        }
        Ok(Obstacle {
            collocation,
            n_points,
            gaps,
            weights,
            stiffness,
            alpha,
            provenance,
        })
    }

    /// Bypass admission — for mutation batteries and trusted callers
    /// ONLY (the fs-phs `from_raw_parts` pattern); nothing downstream
    /// re-checks.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn from_raw_parts(
        collocation: Vec<f64>,
        n_points: usize,
        gaps: Vec<f64>,
        weights: Vec<f64>,
        stiffness: f64,
        alpha: f64,
        provenance: String,
    ) -> Self {
        Obstacle {
            collocation,
            n_points,
            gaps,
            weights,
            stiffness,
            alpha,
            provenance,
        }
    }

    /// Row-major collocation matrix `Phi[i][k]` (n_points x n_modes).
    #[must_use]
    pub fn collocation(&self) -> &[f64] {
        &self.collocation
    }

    /// Collocation point count.
    #[must_use]
    pub fn n_points(&self) -> usize {
        self.n_points
    }

    /// Per-point gaps.
    #[must_use]
    pub fn gaps(&self) -> &[f64] {
        &self.gaps
    }

    /// Per-point quadrature weights.
    #[must_use]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Contact stiffness `K`.
    #[must_use]
    pub fn stiffness(&self) -> f64 {
        self.stiffness
    }

    /// Contact exponent `alpha`.
    #[must_use]
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Contact-law provenance (logged, never invented).
    #[must_use]
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Penetration `p_i = (Phi q)_i - c_i`, clamped one-sided, for the
    /// interleaved `[q, p]` state.
    fn penetrations(&self, n_modes: usize, x: &[f64]) -> Vec<f64> {
        (0..self.n_points)
            .map(|i| {
                let mut disp = 0.0;
                for k in 0..n_modes {
                    disp += self.collocation[i * n_modes + k] * x[2 * k];
                }
                (disp - self.gaps[i]).max(0.0)
            })
            .collect()
    }
}

/// Contact diagnostics over one inspection call.
#[derive(Debug, Clone)]
pub struct ContactProbe {
    /// Points currently in contact.
    pub active_points: usize,
    /// Maximum penetration depth over the obstacle.
    pub max_penetration: f64,
    /// Stored contact energy.
    pub contact_energy: f64,
}

/// An fs-phs `Storage` wrapping an inner storage with distributed
/// contact potentials. State layout is the inner storage's
/// (`[q_0, p_0, ...]` interleaved, `n_modes` pairs).
pub struct ContactStorage {
    inner: Box<dyn Storage>,
    n_modes: usize,
    obstacles: Vec<Obstacle>,
}

impl ContactStorage {
    /// Wrap an inner storage; obstacles must match `n_modes`.
    ///
    /// # Errors
    /// [`DContactError::Shape`] on collocation-width mismatch.
    pub fn new(
        inner: Box<dyn Storage>,
        n_modes: usize,
        obstacles: Vec<Obstacle>,
    ) -> Result<Self, DContactError> {
        for ob in &obstacles {
            if ob.collocation.len() != ob.n_points * n_modes {
                return Err(DContactError::Shape {
                    what: "obstacle collocation width vs n_modes",
                });
            }
        }
        Ok(ContactStorage {
            inner,
            n_modes,
            obstacles,
        })
    }

    /// The wrapped inner storage.
    #[must_use]
    pub fn inner_storage(&self) -> &dyn Storage {
        self.inner.as_ref()
    }

    /// The obstacles.
    #[must_use]
    pub fn obstacles(&self) -> &[Obstacle] {
        &self.obstacles
    }

    /// Modal count.
    #[must_use]
    pub fn n_modes(&self) -> usize {
        self.n_modes
    }

    /// Contact diagnostics at a state (max penetration feeds the
    /// authored penetration ceiling — a FLAG for stiffness
    /// inadequacy, deliberately not a refusal).
    #[must_use]
    pub fn probe(&self, x: &[f64]) -> ContactProbe {
        let mut active = 0usize;
        let mut max_p = 0.0f64;
        let mut energy = 0.0f64;
        for ob in &self.obstacles {
            for (i, &p) in ob.penetrations(self.n_modes, x).iter().enumerate() {
                if p > 0.0 {
                    active += 1;
                    max_p = max_p.max(p);
                    energy += ob.weights[i] * ob.stiffness / (ob.alpha + 1.0)
                        * det::pow(p, ob.alpha + 1.0);
                }
            }
        }
        ContactProbe {
            active_points: active,
            max_penetration: max_p,
            contact_energy: energy,
        }
    }
}

impl Storage for ContactStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let mut h = self.inner.hamiltonian(x);
        for ob in &self.obstacles {
            for (i, &p) in ob.penetrations(self.n_modes, x).iter().enumerate() {
                if p > 0.0 {
                    h += ob.weights[i] * ob.stiffness / (ob.alpha + 1.0)
                        * det::pow(p, ob.alpha + 1.0);
                }
            }
        }
        h
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.inner.gradient(x, out);
        for ob in &self.obstacles {
            for (i, &p) in ob.penetrations(self.n_modes, x).iter().enumerate() {
                if p > 0.0 {
                    let f = ob.weights[i] * ob.stiffness * det::pow(p, ob.alpha);
                    for k in 0..self.n_modes {
                        out[2 * k] += f * ob.collocation[i * self.n_modes + k];
                    }
                }
            }
        }
    }
}

/// Gap profile of a polyline obstacle sampled under given abscissas:
/// linear interpolation of `(x, height)` vertices, refusing samples
/// outside the polyline span. The GAP convention is the caller's
/// (heights are returned as sampled; subtract from the rest position
/// to make gaps).
///
/// # Errors
/// [`DContactError`] on unsorted/short polylines or out-of-span
/// samples.
pub fn polyline_heights(
    vertices: &[(f64, f64)],
    samples: &[f64],
) -> Result<Vec<f64>, DContactError> {
    if vertices.len() < 2 {
        return Err(DContactError::Shape {
            what: "polyline needs at least 2 vertices",
        });
    }
    // NaN slips through ordering comparisons (review finding:
    // executed Ok(NaN)); refuse non-finite coordinates up front.
    if vertices
        .iter()
        .any(|&(x, y)| !x.is_finite() || !y.is_finite())
    {
        return Err(DContactError::Parameter {
            what: "non-finite polyline vertex",
        });
    }
    for pair in vertices.windows(2) {
        if pair[1].0 <= pair[0].0 {
            return Err(DContactError::Parameter {
                what: "polyline x must strictly ascend",
            });
        }
    }
    let (x0, xn) = (vertices[0].0, vertices[vertices.len() - 1].0);
    let mut out = Vec::with_capacity(samples.len());
    for &s in samples {
        if s.is_nan() || s < x0 || s > xn {
            return Err(DContactError::Parameter {
                what: "sample outside polyline span",
            });
        }
        // Find the segment (linear scan: obstacle profiles are small).
        let mut h = vertices[vertices.len() - 1].1;
        for pair in vertices.windows(2) {
            if s <= pair[1].0 {
                let t = (s - pair[0].0) / (pair[1].0 - pair[0].0);
                h = pair[0].1 + t * (pair[1].1 - pair[0].1);
                break;
            }
        }
        out.push(h);
    }
    Ok(out)
}

/// Mass-normalized sine-mode collocation matrix for a fixed-fixed
/// string: `Phi[i][k] = sqrt(2/(mu L)) sin((k+1) pi x_i / L)` — the
/// standard bridge from modal string states to physical displacement
/// at the obstacle points.
/// # Errors
/// [`DContactError::Parameter`] on non-positive length/density (the
/// silent-NaN alternative violates the refusal doctrine — review
/// finding).
pub fn string_collocation(
    length: f64,
    lin_density: f64,
    points: &[f64],
    n_modes: usize,
) -> Result<Vec<f64>, DContactError> {
    if length.is_nan() || length <= 0.0 || lin_density.is_nan() || lin_density <= 0.0 {
        return Err(DContactError::Parameter {
            what: "string length and density must be positive",
        });
    }
    let norm = det::sqrt(2.0 / (lin_density * length));
    let pi = core::f64::consts::PI;
    let mut m = vec![0.0; points.len() * n_modes];
    for (i, &x) in points.iter().enumerate() {
        for k in 0..n_modes {
            m[i * n_modes + k] = norm * det::sin((k + 1) as f64 * pi * x / length);
        }
    }
    Ok(m)
}
