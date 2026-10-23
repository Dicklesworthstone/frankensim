//! Complete solid-mode radiation loading, with exp(-i omega t).
//!
//! The existing Helmholtz batch and power-consistent Rayleigh owners produce
//! every matrix entry. The same physical trace appears in velocity and work:
//! Z_rs = integral phi_r p_s dS. No mode, cross term or negative eigenvalue is
//! discarded. A nonreciprocal or active sampled result refuses; this adapter
//! does not repair collocation data or claim a positive-real time realization.
use crate::helmholtz::rayleigh_modal::{self, AreaSample, PreparedSurface};
use crate::helmholtz::{
    self, Formulation, HelmholtzError, IterativeError, IterativeOptions, IterativeRadiation,
    Medium, PanelRule, RadiationSolution,
};
use crate::panel3d::SpherePanels;
use crate::radiation_policy::GeometryPolicy;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use std::ops::Range;

/// Physical boundary and complete outward normal-velocity traces. A structure
/// must supply a work-consistent transfer; this module does not interpolate it.
#[derive(Clone, Copy)]
pub enum SurfaceModeTraces<'a> {
    /// Closed, outward, retained triangle surface; rows are modes, columns panels.
    ClosedBody {
        surface: &'a SpherePanels,
        shapes: &'a [Vec<f64>],
    },
    /// One radiating half-space; each cell retains every mode and its SI area.
    PlanarBaffled { samples: &'a [AreaSample] },
}

/// Explicit numerical owner; selecting one never falls back to the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadingRoute {
    /// Resonance-safe combined formulation. Only BurtonMiller is admitted;
    /// neither the plain interior-resonant nor mutation arm is a loading route.
    ClosedBodyBem(Formulation),
    /// Per frequency row, [`GeometryPolicy`] selects plain CBIE below 80% of
    /// the smallest closed component's interior Dirichlet box bound, and
    /// Burton-Miller above it. Below ka of about 0.5, Burton-Miller's
    /// alpha = i/k converges slowly for non-uniform modes. On the unit
    /// icosphere at ka = 0.05 and 5120 panels, its dipole resistance is
    /// 11.6 times the closed form, against plain CBIE's 0.2% (ts1,
    /// cdf20cf10).
    ClosedBodyGeometryPolicy,
    /// Original complete Rayleigh matrix with outgoing-power Gram resistance.
    PlanarBaffledRayleigh,
}

/// Source or loading-admission refusal. No partial bake is returned. A refused
/// sampled matrix travels as a diagnostic only; it is never a bake row.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalLoadingError {
    /// Original numerical owner's refusal.
    Helmholtz(HelmholtzError),
    /// Invalid route, mode, frequency, provenance or output budget.
    Invalid(&'static str),
    /// Dense geometry must be handled by the named larger-solver work item.
    SolverRoute { panels: usize, route: &'static str },
    /// The routed iterative owner refused a row or did not converge a field.
    Iterative(String),
    /// Original data fail reciprocity, measured relative to the Frobenius norm.
    NonReciprocal {
        omega: f64,
        relative_defect: f64,
        /// The refused row-major matrix, unchanged.
        impedance: Vec<C64>,
    },
    /// Symmetric resistance has a resolved negative eigenvalue.
    NonPassive {
        omega: f64,
        minimum_eigenvalue: f64,
        norm: f64,
        /// The refused row-major matrix, unchanged.
        impedance: Vec<C64>,
    },
}
impl From<HelmholtzError> for ModalLoadingError {
    fn from(e: HelmholtzError) -> Self {
        Self::Helmholtz(e)
    }
}
impl core::fmt::Display for ModalLoadingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Helmholtz(e) => core::fmt::Display::fmt(e, f),
            Self::Invalid(s) => write!(f, "modal radiation loading: {s}"),
            Self::SolverRoute { panels, route } => write!(f, "{panels} panels require {route}"),
            Self::Iterative(e) => write!(f, "routed modal loading: {e}"),
            Self::NonReciprocal {
                omega,
                relative_defect,
                ..
            } => write!(
                f,
                "nonreciprocal loading at {omega} rad/s: {relative_defect}"
            ),
            Self::NonPassive {
                omega,
                minimum_eigenvalue,
                norm,
                ..
            } => write!(
                f,
                "active loading at {omega} rad/s: min(Re Z)={minimum_eigenvalue}, norm={norm}"
            ),
        }
    }
}
impl std::error::Error for ModalLoadingError {}

/// Allocation-free cap admission. A 257-mode BEM request is two complete source
/// batches; it never becomes a 256-mode physical model. The output budget counts
/// complex matrix entries over the ENTIRE sweep, not bytes or wall time.
pub fn admit_loading_work(
    sites: usize,
    modes: usize,
    frequencies: usize,
    route: LoadingRoute,
    maximum_matrix_entries: usize,
) -> Result<usize, ModalLoadingError> {
    admit(
        sites,
        modes,
        frequencies,
        route,
        maximum_matrix_entries,
        false,
    )
}

/// `routed` lifts only the dense cap, and only for closed bodies.
fn admit(
    sites: usize,
    modes: usize,
    frequencies: usize,
    route: LoadingRoute,
    maximum_matrix_entries: usize,
    routed: bool,
) -> Result<usize, ModalLoadingError> {
    if sites == 0 || modes == 0 || frequencies == 0 {
        return Err(ModalLoadingError::Invalid(
            "nonempty surface, basis and frequency grid required",
        ));
    }
    if sites > helmholtz::MAX_DENSE_PANELS
        && !(routed && route != LoadingRoute::PlanarBaffledRayleigh)
    {
        return Err(ModalLoadingError::SolverRoute {
            panels: sites,
            route: "M2.5-W5c",
        });
    }
    let entries = modes
        .checked_mul(modes)
        .and_then(|n| n.checked_mul(frequencies))
        .ok_or(ModalLoadingError::Invalid(
            "modal sweep output shape overflow",
        ))?;
    if entries > maximum_matrix_entries {
        return Err(ModalLoadingError::Invalid(
            "modal sweep exceeds declared matrix-entry budget",
        ));
    }
    fs_la::eigen::admit_jacobi_eigh(modes).map_err(|_| {
        ModalLoadingError::Invalid("modal resistance eigenproblem exceeds work envelope")
    })?;
    match route {
        LoadingRoute::PlanarBaffledRayleigh => {
            if modes > rayleigh_modal::MAX_MODAL_RAYLEIGH_FIELDS {
                return Err(ModalLoadingError::Invalid(
                    "Rayleigh loading requires at most 128 modes",
                ));
            }
            let (n, m) = (sites as u128, modes as u128);
            rayleigh_modal::admit_modal_work(n * m * (n + m))?;
            Ok(1)
        }
        LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller)
        | LoadingRoute::ClosedBodyGeometryPolicy => {
            Ok(modes.div_ceil(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH))
        }
        LoadingRoute::ClosedBodyBem(_) => Err(ModalLoadingError::Invalid(
            "closed-body loading requires the resonance-safe Burton-Miller owner",
        )),
    }
}

/// Deterministic, exhaustive source-column partition.
pub fn source_batches(modes: usize) -> impl Iterator<Item = Range<usize>> {
    (0..modes)
        .step_by(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH)
        .map(move |start| {
            start
                ..start
                    .saturating_add(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH)
                    .min(modes)
        })
}

/// Original closed-body solve diagnostics over every source field in one row.
/// These are resolution and conditioning warnings, not error enclosures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalLoadingBemDiagnostics {
    /// Formulation actually selected at this row's frequency.
    pub formulation: Formulation,
    /// Smallest original panels-per-wavelength diagnostic over all source batches.
    pub minimum_panels_per_wavelength: f64,
    /// Largest original probe-based condition lower bound over all source batches.
    /// The iterative owner supplies the universal lower bound 1, without probes.
    pub maximum_condition_lower_bound: f64,
}

/// Sampled work impedance [generalized force / generalized velocity]. For a
/// length coordinate Z is kg/s, resistance is Re Z, and added mass is -Im Z/w.
/// This is a frequency sample, not a fitted causal realization.
#[derive(Debug, Clone)]
pub struct ModalLoadingRow {
    omega: f64,
    z_modal: Vec<C64>,
    minimum_resistance_eigenvalue: f64,
    eigen_residual: f64,
    reciprocity_relative_defect: f64,
    bem_diagnostics: Option<ModalLoadingBemDiagnostics>,
}
impl ModalLoadingRow {
    pub const fn omega(&self) -> f64 {
        self.omega
    }
    pub fn impedance(&self) -> &[C64] {
        &self.z_modal
    }
    pub const fn minimum_resistance_eigenvalue(&self) -> f64 {
        self.minimum_resistance_eigenvalue
    }
    /// Frobenius residual of the complete resistance eigendecomposition, in Z units.
    pub const fn eigen_residual(&self) -> f64 {
        self.eigen_residual
    }
    pub const fn reciprocity_relative_defect(&self) -> f64 {
        self.reciprocity_relative_defect
    }
    /// Diagnostics retained from the same closed-body solves that supplied Z.
    /// Rayleigh rows have no BEM solve and return `None`.
    pub const fn bem_diagnostics(&self) -> Option<ModalLoadingBemDiagnostics> {
        self.bem_diagnostics
    }
}

/// Immutable complete sweep. Declaration identity binds all source geometry,
/// traces, medium, frequencies, route, provenance and output budget. Consumers
/// can pass a row verbatim to VibroacousticModel or the W3b realization owner.
#[derive(Debug, Clone)]
pub struct ModalLoadingBake {
    modes: usize,
    declaration_hash: ContentHash,
    rows: Vec<ModalLoadingRow>,
}
impl ModalLoadingBake {
    pub const fn modes(&self) -> usize {
        self.modes
    }
    pub const fn declaration_hash(&self) -> ContentHash {
        self.declaration_hash
    }
    pub fn rows(&self) -> &[ModalLoadingRow] {
        &self.rows
    }
    /// Only exact declared frequencies are admitted; interpolation and band
    /// extension must be assessed by the consuming realization owner.
    pub fn at_omega(&self, omega: f64) -> Result<&ModalLoadingRow, ModalLoadingError> {
        self.rows
            .binary_search_by(|r| r.omega.total_cmp(&omega))
            .map(|i| &self.rows[i])
            .map_err(|_| ModalLoadingError::Invalid("frequency is not a declared loading sample"))
    }
}

fn reserved<T>(n: usize) -> Result<Vec<T>, ModalLoadingError> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| HelmholtzError::AllocationFailed {
            what: "complete modal loading",
        })?;
    Ok(v)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ModalLoadingError> {
    if cancelled() {
        Err(HelmholtzError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn assess(omega: f64, z: Vec<C64>, n: usize) -> Result<ModalLoadingRow, ModalLoadingError> {
    if z.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
        return Err(ModalLoadingError::Invalid("nonfinite projected impedance"));
    }
    let norm = z.iter().fold(0.0_f64, |s, v| s.hypot(v.abs()));
    if !norm.is_finite() {
        return Err(ModalLoadingError::Invalid("impedance norm overflow"));
    }
    let mut defect = 0.0_f64;
    let mut resistance = reserved(n * n)?;
    for i in 0..n {
        for j in 0..n {
            defect = defect.hypot((z[i * n + j] - z[j * n + i]).abs());
            // Symmetric part is used only for the Hermitian-work diagnostic;
            // the returned physical matrix is unchanged, including every sign.
            resistance.push(0.5 * z[i * n + j].re + 0.5 * z[j * n + i].re);
        }
    }
    if defect > 1e-12 * norm {
        return Err(ModalLoadingError::NonReciprocal {
            omega,
            relative_defect: defect / norm,
            impedance: z,
        });
    }
    if norm == 0.0 {
        return Ok(ModalLoadingRow {
            omega,
            z_modal: z,
            minimum_resistance_eigenvalue: 0.0,
            eigen_residual: 0.0,
            reciprocity_relative_defect: 0.0,
            bem_diagnostics: None,
        });
    }
    for value in &mut resistance {
        *value /= norm;
    }
    let (values, vectors) = fs_la::eigen::jacobi_eigh(&resistance, n);
    let minimum = values[0];
    let mut residual = 0.0_f64;
    let mut orthogonality = 0.0_f64;
    for column in 0..n {
        for i in 0..n {
            let av = (0..n)
                .map(|j| resistance[i * n + j] * vectors[j * n + column])
                .sum::<f64>();
            residual = residual.hypot(av - values[column] * vectors[i * n + column]);
            let dot = (0..n)
                .map(|j| vectors[j * n + i] * vectors[j * n + column])
                .sum::<f64>();
            orthogonality = orthogonality.hypot(dot - if i == column { 1.0 } else { 0.0 });
        }
    }
    if values.iter().any(|v| !v.is_finite())
        || !residual.is_finite()
        || residual > 1e-13
        || !orthogonality.is_finite()
        || orthogonality > 1e-10
    {
        return Err(ModalLoadingError::Invalid(
            "resistance eigendecomposition is unresolved",
        ));
    }
    let eigen_error = residual / (1.0 - orthogonality).sqrt();
    if minimum - eigen_error - 0.5 * defect / norm < -1e-12 {
        return Err(ModalLoadingError::NonPassive {
            omega,
            minimum_eigenvalue: minimum * norm,
            norm,
            impedance: z,
        });
    }
    Ok(ModalLoadingRow {
        omega,
        z_modal: z,
        minimum_resistance_eigenvalue: minimum * norm,
        eigen_residual: residual * norm,
        reciprocity_relative_defect: defect / norm,
        bem_diagnostics: None,
    })
}

/// Bake all modes at every declared positive increasing angular frequency.
/// Cancellation is polled before geometry preparation, at frequency and source
/// batch boundaries, and before publication. The original dense LU has no
/// mid-factorization cancellation; no tighter latency guarantee is made.
pub fn modal_radiation_loading(
    traces: SurfaceModeTraces<'_>,
    route: LoadingRoute,
    omegas: &[f64],
    medium: Medium,
    provenance: ContentHash,
    maximum_matrix_entries: usize,
    cancelled: impl FnMut() -> bool,
) -> Result<ModalLoadingBake, ModalLoadingError> {
    load(
        traces,
        route,
        omegas,
        medium,
        provenance,
        maximum_matrix_entries,
        None,
        cancelled,
    )
}

/// [`modal_radiation_loading`] on the M2.5-W5 route past the dense cap.
/// - At or below [`helmholtz::MAX_DENSE_PANELS`] it is the dense bake under
///   `iterative`'s declared panel rule, and the rest of `iterative` is unused.
///   Under the centroid rule (or no options) that is the dense bake bitwise,
///   declaration identity included. Any other rule is
///   [`helmholtz::solve_radiation_batch_with_rule`]'s bake, and the rule is
///   bound into the declaration.
/// - Above it, a closed body's rows go to [`IterativeRadiation`] with
///   `iterative`, each in its own formulation: plain CBIE on the dense
///   owner's triangle quadrature, Burton-Miller on centroid kernels. Every
///   field must converge to the declared GMRES tolerance, or the bake
///   refuses. The declaration identity then binds the numerical options,
///   panel rule included (not workers, caches or the rule's pair budget,
///   which change no bits).
///
/// Cancellation is also polled between fields; one GMRES solve is not
/// interrupted by `cancelled`.
pub fn modal_radiation_loading_routed(
    traces: SurfaceModeTraces<'_>,
    route: LoadingRoute,
    omegas: &[f64],
    medium: Medium,
    provenance: ContentHash,
    maximum_matrix_entries: usize,
    iterative: Option<&IterativeOptions>,
    cancelled: impl FnMut() -> bool,
) -> Result<ModalLoadingBake, ModalLoadingError> {
    load(
        traces,
        route,
        omegas,
        medium,
        provenance,
        maximum_matrix_entries,
        iterative,
        cancelled,
    )
}

/// Declaration tag of a panel rule; its pair budget changes no bits.
const fn rule_tag(rule: PanelRule) -> u64 {
    match rule {
        PanelRule::Centroid => 0,
        PanelRule::TriangleNear { .. } => 1,
        PanelRule::Triangle { .. } => 2,
    }
}

fn iterative_error(e: IterativeError) -> ModalLoadingError {
    match e {
        IterativeError::Helmholtz(e) => ModalLoadingError::Helmholtz(e),
        other => ModalLoadingError::Iterative(other.to_string()),
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one admission audit
fn load(
    traces: SurfaceModeTraces<'_>,
    route: LoadingRoute,
    omegas: &[f64],
    medium: Medium,
    provenance: ContentHash,
    maximum_matrix_entries: usize,
    iterative: Option<&IterativeOptions>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ModalLoadingBake, ModalLoadingError> {
    poll(&mut cancelled)?;
    if !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
        || omegas.iter().any(|w| !w.is_finite() || *w <= 0.0)
        || omegas.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(ModalLoadingError::Invalid(
            "positive finite medium and strictly increasing positive frequencies required",
        ));
    }
    let (sites, modes) = match traces {
        SurfaceModeTraces::ClosedBody { surface, shapes } => (surface.areas().len(), shapes.len()),
        SurfaceModeTraces::PlanarBaffled { samples } => (
            samples.len(),
            samples.first().map_or(0, |s| s.mode_shape.len()),
        ),
    };
    admit(
        sites,
        modes,
        omegas.len(),
        route,
        maximum_matrix_entries,
        iterative.is_some(),
    )?;
    // Below the cap only the declared panel rule reaches the dense bake.
    let dense_rule = iterative
        .filter(|_| sites <= helmholtz::MAX_DENSE_PANELS)
        .map_or(PanelRule::Centroid, |o| o.panel_rule);
    let iterative = iterative.filter(|_| sites > helmholtz::MAX_DENSE_PANELS);
    let mut hash = DomainHasher::new("fs-bem/complete-modal-loading/v1");
    hash.update(provenance.as_bytes());
    for v in [
        sites as u64,
        modes as u64,
        omegas.len() as u64,
        maximum_matrix_entries as u64,
    ] {
        hash.update(&v.to_le_bytes());
    }
    for v in [medium.density, medium.sound_speed]
        .into_iter()
        .chain(omegas.iter().copied())
    {
        hash.update(&v.to_bits().to_le_bytes());
    }
    let planar = match (traces, route) {
        (
            SurfaceModeTraces::ClosedBody { surface, shapes },
            LoadingRoute::ClosedBodyBem(_) | LoadingRoute::ClosedBodyGeometryPolicy,
        ) => {
            if shapes
                .iter()
                .any(|s| s.len() != sites || s.iter().any(|v| !v.is_finite()))
            {
                return Err(ModalLoadingError::Invalid(
                    "complete finite panel traces required",
                ));
            }
            hash.update(&[match route {
                LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller) => 1,
                LoadingRoute::ClosedBodyGeometryPolicy => 4,
                _ => 2,
            }]);
            // The centroid dense bake keeps its original identity.
            if dense_rule != PanelRule::Centroid {
                hash.update(&[6]);
                hash.update(&rule_tag(dense_rule).to_le_bytes());
            }
            if let Some(o) = iterative {
                hash.update(&[5]);
                let precond = match o.precond {
                    helmholtz::Preconditioner::None => 0,
                    helmholtz::Preconditioner::NearFieldBlockJacobi { max_block_panels } => {
                        max_block_panels as u64
                    }
                };
                for v in [
                    o.tol.to_bits(),
                    o.restart as u64,
                    o.max_iters as u64,
                    o.fmm_order as u64,
                    o.leaf_capacity as u64,
                    precond,
                    o.max_translation_work as u64,
                    u64::from(o.directional.is_some()),
                    // The rule changes the numbers; its pair budget does not.
                    rule_tag(o.panel_rule),
                ] {
                    hash.update(&v.to_le_bytes());
                }
                if let Some(d) = &o.directional {
                    for v in [
                        d.tolerance.to_bits(),
                        d.order as u64,
                        d.leaf_capacity as u64,
                        d.low_freq_threshold.to_bits(),
                        d.cone_resolution as u64,
                        d.max_cones_per_level as u64,
                        d.compression_tolerance.to_bits(),
                        d.max_rank as u64,
                    ] {
                        hash.update(&v.to_le_bytes());
                    }
                }
            }
            // Reuse the original exact closure/orientation admission. The
            // exterior probe is geometry only, never an acoustic observation.
            let triangles = surface.triangles().ok_or(ModalLoadingError::Invalid(
                "closed-body loading requires retained triangles",
            ))?;
            if sites > helmholtz::MAX_DENSE_PANELS {
                // Routed past the dense cap, which near_field's receiver
                // admission shares: the same closed outward-component owner.
                helmholtz::bc::validate_closed(surface)?;
            } else {
                let mut low = [f64::INFINITY; 3];
                let mut high = [f64::NEG_INFINITY; 3];
                for p in triangles.iter().flatten() {
                    for i in 0..3 {
                        low[i] = low[i].min(p[i]);
                        high[i] = high[i].max(p[i]);
                    }
                }
                let length = (0..3).map(|i| high[i] - low[i]).fold(0.0_f64, f64::max);
                let outside = [high[0] + 2.0 * length, high[1], high[2]];
                crate::near_field::Geometry::new(
                    surface,
                    &[outside],
                    64.0 * f64::EPSILON * length,
                )?;
            }
            for v in triangles
                .iter()
                .flatten()
                .flatten()
                .chain(shapes.iter().flatten())
            {
                hash.update(&v.to_bits().to_le_bytes());
            }
            None
        }
        (SurfaceModeTraces::PlanarBaffled { samples }, LoadingRoute::PlanarBaffledRayleigh) => {
            hash.update(&[3]);
            for s in samples {
                for v in s.position_m.iter().chain([&s.area_m2]).chain(&s.mode_shape) {
                    hash.update(&v.to_bits().to_le_bytes());
                }
            }
            Some(PreparedSurface::new(samples)?)
        }
        _ => {
            return Err(ModalLoadingError::Invalid(
                "surface and declared loading route disagree",
            ));
        }
    };
    let policy = match (traces, route) {
        (SurfaceModeTraces::ClosedBody { surface, .. }, LoadingRoute::ClosedBodyGeometryPolicy) => {
            Some(GeometryPolicy::new(surface)?)
        }
        _ => None,
    };
    let mut rows = reserved(omegas.len())?;
    for &omega in omegas {
        poll(&mut cancelled)?;
        let mut bem_diagnostics = None;
        let z = match traces {
            SurfaceModeTraces::PlanarBaffled { .. } => planar
                .as_ref()
                .ok_or(ModalLoadingError::Invalid("missing planar preparation"))?
                .power_consistent_impedance(omega, medium)?,
            SurfaceModeTraces::ClosedBody { surface, shapes } => {
                let formulation = match (route, &policy) {
                    (LoadingRoute::ClosedBodyBem(formulation), _) => formulation,
                    (LoadingRoute::ClosedBodyGeometryPolicy, Some(policy)) => {
                        policy.formulation(omega / medium.sound_speed)?
                    }
                    _ => return Err(ModalLoadingError::Invalid("closed-body route changed")),
                };
                let mut z = reserved(modes * modes)?;
                z.resize(modes * modes, C64::ZERO);
                let mut diagnostics = ModalLoadingBemDiagnostics {
                    formulation,
                    minimum_panels_per_wavelength: f64::INFINITY,
                    maximum_condition_lower_bound: 0.0,
                };
                let mut project = |column: usize, solution: &RadiationSolution| {
                    diagnostics.minimum_panels_per_wavelength = diagnostics
                        .minimum_panels_per_wavelength
                        .min(solution.panels_per_wavelength);
                    diagnostics.maximum_condition_lower_bound = diagnostics
                        .maximum_condition_lower_bound
                        .max(solution.condition_lower_bound);
                    for (row, shape) in shapes.iter().enumerate() {
                        let mut force = C64::ZERO;
                        for ((&phi, &p), &area) in
                            shape.iter().zip(&solution.pressure).zip(surface.areas())
                        {
                            force = force + p.scale(phi * area);
                        }
                        z[row * modes + column] = force;
                    }
                };
                let k = omega / medium.sound_speed;
                if let Some(options) = iterative {
                    let gate = CancelGate::new_clock_free();
                    let operator = IterativeRadiation::prepare(
                        surface,
                        k,
                        medium,
                        formulation,
                        options,
                        &gate,
                    )
                    .map_err(iterative_error)?;
                    for (column, shape) in shapes.iter().enumerate() {
                        poll(&mut cancelled)?;
                        let field: Vec<C64> = shape.iter().map(|&v| C64::from_re(v)).collect();
                        let mut state = operator.start(&field, &gate).map_err(iterative_error)?;
                        operator
                            .advance(&mut state, options.max_iters.saturating_add(1), &gate)
                            .map_err(iterative_error)?;
                        let (solution, _, _) = operator.finish(&state).map_err(iterative_error)?;
                        project(column, &solution);
                    }
                } else {
                    for range in source_batches(modes) {
                        poll(&mut cancelled)?;
                        let fields: Vec<Vec<C64>> = shapes[range.clone()]
                            .iter()
                            .map(|r| r.iter().map(|v| C64::from_re(*v)).collect())
                            .collect();
                        let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
                        // The centroid rule is solve_radiation_batch, bit for bit.
                        let solved = helmholtz::solve_radiation_batch_with_rule(
                            surface,
                            k,
                            medium,
                            &refs,
                            formulation,
                            dense_rule,
                        )?;
                        for (column, solution) in range.zip(solved) {
                            project(column, &solution);
                        }
                    }
                }
                bem_diagnostics = Some(diagnostics);
                z
            }
        };
        let mut row = assess(omega, z, modes)?;
        row.bem_diagnostics = bem_diagnostics;
        rows.push(row);
    }
    poll(&mut cancelled)?;
    Ok(ModalLoadingBake {
        modes,
        declaration_hash: hash.finalize(),
        rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn active_and_nonreciprocal_data_are_refused_without_repair() {
        assert!(matches!(
            assess(1.0, vec![C64::from_re(-1.0)], 1),
            Err(ModalLoadingError::NonPassive { .. })
        ));
        assert!(matches!(
            assess(1.0, vec![C64::ONE, C64::ONE, C64::ZERO, C64::ONE], 2),
            Err(ModalLoadingError::NonReciprocal { .. })
        ));
        let z = vec![
            C64::from_re(1.0),
            C64::from_re(-1.0),
            C64::from_re(-1.0),
            C64::from_re(1.0),
        ];
        assert_eq!(assess(1.0, z.clone(), 2).unwrap().impedance(), z);
    }
}
