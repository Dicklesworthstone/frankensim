//! Free decay on the original hereditary spatial equations, with exp(s t).
//!
//! Maxwell strain q obeys (1+s*tau)q=s*tau*B*u. Keep these ORIGINAL six
//! Mandel coordinates at every source point/branch instead of dividing by
//! 1+s*tau: the resulting quadratic pencil is holomorphic even at material
//! relaxation rates. Mechanical rows are (K_eq+s*C_kv+s²*M)u + Σ V BᵀC_j q.
//! No material averaging, branch merging, modal truncation or damping fit.
//! The only eigenvalue iterations are the shared fs-la Newton/Beyn owners.

use super::{ViscoElementError, ViscoElements, bilinear, checkpoint};
use fs_exec::Cx;
use fs_la::nonlinear_eig::beyn::{BeynError, BeynOptions, BeynReport, Ellipse, beyn_with_cancel};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, NonlinearPole, PoleEvaluation, PoleOptions,
    newton_pole_with_cancel,
};
use fs_math::c64::C64;

mod load;
use load::PreparedPortLoad;
pub use load::{ViscoPolePortBudget, ViscoPolePortLoad, ViscoPolePortResponse};
mod dynamic;
use dynamic::PreparedDynamicLoad;
pub use dynamic::{ViscoPoleDynamicBudget, ViscoPoleDynamicLoad, ViscoPoleDynamicResponse};

mod response;
pub use response::{
    ViscoJointHarmonicBudget, ViscoJointHarmonicError, ViscoJointHarmonicFactor,
    ViscoJointHarmonicLoad, ViscoJointHarmonicResponse,
};

/// Source-owned preparation and callback resources, separate from fs-la's
/// numerical iteration budgets. Allocator overhead/process RSS are excluded.
#[derive(Clone, Copy, Debug)]
pub struct ViscoPoleBudget {
    /// Every free displacement, original branch strain, and independent load state.
    pub maximum_coordinates: usize,
    /// Coefficient caches, fixed scales, map, and one live callback evaluation.
    pub maximum_payload_bytes: usize,
    /// Conservative preparation work units; checked before dense allocation.
    pub maximum_assembly_work: u64,
    /// Callback evaluations per query, including final physical acceptance.
    pub maximum_evaluations: usize,
    /// Positive fixed balancing scale [1/s], independent of contour nodes.
    pub reference_rate_per_s: f64,
    /// Blockwise coefficient-scaled residual ceiling in EACH physical family.
    pub physical_residual_tolerance: f64,
}

/// Location of six consecutive original Mandel Maxwell-strain coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ViscoMemorySite {
    /// Original element index.
    pub element: usize,
    /// Original material quadrature-point index in that element.
    pub point: usize,
    /// Original material branch index; equal-time branches remain separate.
    pub branch: usize,
}

/// A finite pole of the complete retained hereditary equations.
#[derive(Clone, Debug)]
pub struct ViscoPole {
    /// Laplace rate [1/s] under exp(s*t); negative real part means decay.
    /// Pure real roots describe relaxation, not an oscillatory resonance.
    pub rate_per_s: C64,
    /// Shape in original free displacement coordinates, followed by the six
    /// Mandel strains for each memory_sites() entry, then every independent
    /// dynamic-load state in its original order. Arbitrary common modal
    /// normalization; this is NOT a finite-amplitude admitted material history.
    /// A memory-only or unobservable pole is not evidence of audible resonance.
    pub physical_vector: Vec<C64>,
    /// Mechanical force-row residual, scaled by the original coefficient terms.
    pub mechanical_residual: f64,
    /// Original Maxwell-strain equation residual; zero with no memory rows.
    pub memory_residual: f64,
    /// Original independent load-state equations; zero without a dynamic load.
    pub dynamic_load_residual: f64,
    /// Original attached port motion, storage and loss at the same arbitrary
    /// modal amplitude. Absent for an unloaded preparation.
    pub port_load: Option<ViscoPolePortResponse>,
    /// Independent load coordinates, physical work, storage and dissipation.
    pub dynamic_load: Option<ViscoPoleDynamicResponse>,
    /// Shared solver's scaled final Newton correction.
    pub relative_step: f64,
    /// Number of shared Newton evaluations.
    pub iterations: usize,
}

/// Numerical contour discovery, not an interval completeness certificate.
#[derive(Clone, Debug)]
pub struct ViscoPoleSpectrum {
    /// Original-coordinate poles passing both physical equation families.
    pub poles: Vec<ViscoPole>,
    /// Seed rates before local refinement, for a quadrature convergence ladder.
    pub unrefined_rates_per_s: Vec<C64>,
    /// Fine and nested coarse argument-principle integrals.
    pub argument_integrals: [C64; 2],
    /// Numerically resolved count and moment rank; the shared owner requires equality.
    pub argument_count: usize,
    /// Moment rank returned by the shared owner.
    pub rank: usize,
    /// Fine-grid node count.
    pub nodes: usize,
    /// Relative change of the two nested moment estimates.
    pub moment_relative_change: f64,
    /// Smallest numerical reciprocal conditioning margin along the contour.
    pub minimum_reciprocal_condition: f64,
}

/// Original source, resources, numerical solver, or physical-row refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ViscoPoleError {
    /// Invalid or unrepresentable data; no silent rescaling of the physical model.
    Invalid(&'static str),
    /// Original source/cancellation admission.
    Source(ViscoElementError),
    /// Source-owned resource cap.
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    /// |s| lies outside the intersection of original material frequency bands.
    SourceBand { rate: C64, band_rad_s: (f64, f64) },
    /// A supplied boundary coefficient is not numerically positive semidefinite.
    LoadNotPassive {
        coefficient: &'static str,
        minimum_scaled_eigenvalue: f64,
    },
    /// Displacement-only initialization cannot divide through this branch rate.
    SeedAtRelaxationRate { site: usize },
    /// An original mechanical, memory or independent load equation family failed.
    PhysicalResidual {
        mechanical: f64,
        memory: f64,
        dynamic_load: f64,
        tolerance: f64,
    },
    /// Shared local simple-pole solver refusal.
    Pole(NonlinearEigError),
    /// Shared contour resolution/count/conditioning/refinement refusal.
    Contour(BeynError),
    /// Cooperative cancellation; no partial spectrum is published.
    Cancelled,
}
impl core::fmt::Display for ViscoPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-VISCO-POLES: {self:?}")
    }
}
impl std::error::Error for ViscoPoleError {}
impl From<ViscoElementError> for ViscoPoleError {
    fn from(error: ViscoElementError) -> Self {
        Self::Source(error)
    }
}
fn invalid(what: &'static str) -> ViscoPoleError {
    ViscoPoleError::Invalid(what)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), ViscoPoleError> {
    if cancelled() {
        Err(ViscoPoleError::Cancelled)
    } else {
        Ok(())
    }
}
fn cap(quantity: &'static str, requested: u128, maximum: u128) -> Result<(), ViscoPoleError> {
    if requested > maximum {
        Err(ViscoPoleError::Budget {
            quantity,
            requested,
            cap: maximum,
        })
    } else {
        Ok(())
    }
}

/// Prepared full spatial hereditary pencil. The three real coefficient matrices
/// and fixed balancing are computed once; each subsequent evaluation is O(N²).
/// Dense numerical refinement/discovery remains governed by the shared owners.
pub struct PreparedViscoPoles<'a> {
    source: &'a ViscoElements,
    sites: Vec<ViscoMemorySite>,
    coefficients: [Vec<f64>; 3],
    row_scales: Vec<f64>,
    column_scales: Vec<f64>,
    dimension: usize,
    band: (f64, f64),
    budget: ViscoPoleBudget,
    port_load: Option<PreparedPortLoad>,
    dynamic_load: Option<PreparedDynamicLoad>,
}

impl ViscoElements {
    /// Prepare free decay without freezing K'(omega) or fitting modal damping.
    /// Full source quadrature, consistent mass, Kelvin viscosity, and every
    /// original six-component Maxwell strain survive. Fixed supports are
    /// homogeneous; this linearized image does not choose excitation amplitude.
    ///
    /// The source's real-frequency bands are retained as a conservative |s|
    /// scale gate. Analytic continuation is a property of the declared Prony
    /// model, not an experimental validation at complex frequency. Pure-memory
    /// poles, multiplicities and rank-deficient branch directions are retained.
    pub fn prepare_resonance_poles(
        &self,
        budget: ViscoPoleBudget,
        cx: &Cx<'_>,
    ) -> Result<PreparedViscoPoles<'_>, ViscoPoleError> {
        self.prepare_poles(budget, None, None, cx)
    }

    /// Include a physical passive boundary in the ORIGINAL hereditary pencil.
    /// For port displacement d=B^T*u the opposing effort is
    /// K_port*d + C_port*d' + M_port*d''. Thus the complete mechanical
    /// rows receive B*(K_port+s*C_port+s²*M_port)*B^T before balancing.
    ///
    /// This is the causal spring/damper/inertance law K_port/s+C_port+s*M_port,
    /// not a continuation of samples of Z(omega). Each full port coefficient
    /// must be symmetric PSD; mutual entries and singular passive directions
    /// remain present. All source Maxwell coordinates survive. Port motion is
    /// kinematically tied to the body; an attachment with independent moving
    /// coordinates requires its own coupled state and is not represented here.
    pub fn prepare_resonance_poles_with_port_load(
        &self,
        load: ViscoPolePortLoad<'_>,
        budget: ViscoPoleBudget,
        load_budget: ViscoPolePortBudget,
        cx: &Cx<'_>,
    ) -> Result<PreparedViscoPoles<'_>, ViscoPoleError> {
        self.prepare_poles(budget, Some((load, load_budget)), None, cx)
    }

    /// Couple every state of an admitted exact quadratic pHS impedance load.
    /// With original load law z'=(J-R)Qz+Gv, v=B^T*u', the body receives
    /// opposing force -B*G^T*Qz. The joint polynomial pencil retains z even at
    /// a bare load pole; no standalone impedance or state resolvent is inverted.
    /// The load supplies its physical units and declared continuation band.
    pub fn prepare_resonance_poles_with_dynamic_load(
        &self,
        load: ViscoPoleDynamicLoad<'_>,
        budget: ViscoPoleBudget,
        load_budget: ViscoPoleDynamicBudget,
        cx: &Cx<'_>,
    ) -> Result<PreparedViscoPoles<'_>, ViscoPoleError> {
        self.prepare_poles(budget, None, Some((load, load_budget)), cx)
    }

    fn prepare_poles(
        &self,
        budget: ViscoPoleBudget,
        load: Option<(ViscoPolePortLoad<'_>, ViscoPolePortBudget)>,
        dynamic: Option<(ViscoPoleDynamicLoad<'_>, ViscoPoleDynamicBudget)>,
        cx: &Cx<'_>,
    ) -> Result<PreparedViscoPoles<'_>, ViscoPoleError> {
        checkpoint(cx)?;
        let m = self.equilibrium.free_dofs.len();
        let mut blocks = 0usize;
        let mut work = 0u128;
        let mut band = (0.0_f64, f64::MAX);
        for element in &self.elements {
            checkpoint(cx)?;
            blocks = blocks
                .checked_add(
                    element
                        .points
                        .len()
                        .checked_mul(element.law.branches().len())
                        .ok_or_else(|| invalid("material state extent overflow"))?,
                )
                .ok_or_else(|| invalid("material state extent overflow"))?;
            let d = element.dofs.len() as u128;
            work += element.points.len() as u128
                * (80 * d * d + element.law.branches().len() as u128 * (80 * d + 36));
            let b = element.law.band_rad_s();
            band.0 = band.0.max(b.0);
            band.1 = band.1.min(b.1);
        }
        let body_n = blocks
            .checked_mul(6)
            .and_then(|q| q.checked_add(m))
            .ok_or_else(|| invalid("complete hereditary dimension overflow"))?;
        let n = body_n
            .checked_add(
                dynamic
                    .as_ref()
                    .map_or(0, |(load, _)| load.system.state_dim()),
            )
            .ok_or_else(|| invalid("joint hereditary/load dimension overflow"))?;
        if m == 0
            || budget.maximum_coordinates > MAX_NONLINEAR_POLE_DIMENSION
            || budget.maximum_evaluations == 0
            || band.0 > band.1
            || !budget.reference_rate_per_s.is_finite()
            || budget.reference_rate_per_s <= 0.0
            || !(budget.reference_rate_per_s * budget.reference_rate_per_s).is_finite()
            || !budget.physical_residual_tolerance.is_finite()
            || budget.physical_residual_tolerance <= 0.0
            || budget.physical_residual_tolerance >= 1.0
        {
            return Err(invalid(
                "nonempty source and finite positive pole controls required",
            ));
        }
        cap(
            "complete coordinates",
            n as u128,
            budget.maximum_coordinates as u128,
        )?;
        let nn = n as u128 * n as u128;
        // Three f64 caches, one T/T' pair, scales/maps/sites and bounded scratch.
        cap(
            "source numerical payload bytes",
            80 * nn + 64 * n as u128 + 8 * self.full_dofs as u128,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "preparation work",
            work + 24 * nn,
            budget.maximum_assembly_work as u128,
        )?;
        let port_load = load
            .map(|(load, budget)| PreparedPortLoad::new(m, load, budget, cx))
            .transpose()?;
        let dynamic_load = dynamic
            .map(|(load, budget)| PreparedDynamicLoad::new(m, load, budget, cx))
            .transpose()?;
        if let Some(load) = &port_load {
            band.0 = band.0.max(load.band_rad_s.0);
            band.1 = band.1.min(load.band_rad_s.1);
            if band.0 > band.1 {
                return Err(invalid("material and boundary-load bands do not intersect"));
            }
        }
        if let Some(load) = &dynamic_load {
            band.0 = band.0.max(load.band_rad_s.0);
            band.1 = band.1.min(load.band_rad_s.1);
            if band.0 > band.1 {
                return Err(invalid("material and dynamic-load bands do not intersect"));
            }
        }
        let mut reduced = vec![usize::MAX; self.full_dofs];
        for (i, &dof) in self.equilibrium.free_dofs.iter().enumerate() {
            reduced[dof] = i;
        }
        let mut a = [vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]];
        for (power, matrix) in [
            (0, &self.equilibrium.stiffness),
            (2, &self.equilibrium.mass),
        ] {
            for i in 0..m {
                checkpoint(cx)?;
                let (columns, entries) = matrix.row(i);
                for (&j, &value) in columns.iter().zip(entries) {
                    a[power][i * n + j] = value;
                }
            }
        }
        let mut sites = Vec::with_capacity(blocks);
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                if let Some(viscosity) = element.law.kv_pa_s() {
                    for (i, &dof_i) in element.dofs.iter().enumerate() {
                        let row = reduced[dof_i];
                        if row == usize::MAX {
                            continue;
                        }
                        for (j, &dof_j) in element.dofs.iter().enumerate().skip(i) {
                            let col = reduced[dof_j];
                            if col == usize::MAX {
                                continue;
                            }
                            let value =
                                point.volume_m3 * bilinear(&point.b[i], viscosity, &point.b[j]);
                            a[1][row * n + col] += value;
                            if row != col {
                                a[1][col * n + row] += value;
                            }
                        }
                    }
                }
                for (branch, law) in element.law.branches().iter().enumerate() {
                    checkpoint(cx)?;
                    let at = m + 6 * sites.len();
                    sites.push(ViscoMemorySite {
                        element: e,
                        point: q,
                        branch,
                    });
                    for component in 0..6 {
                        a[0][(at + component) * n + at + component] = 1.0;
                        a[1][(at + component) * n + at + component] = law.tau_s;
                        for (j, &dof) in element.dofs.iter().enumerate() {
                            let row = reduced[dof];
                            if row == usize::MAX {
                                continue;
                            }
                            let coupling = (0..6).fold(0.0, |sum, k| {
                                point.b[j][k].mul_add(law.c_pa[k][component], sum)
                            });
                            a[0][row * n + at + component] = point.volume_m3 * coupling;
                            a[1][(at + component) * n + row] = -law.tau_s * point.b[j][component];
                        }
                    }
                }
            }
        }
        if let Some(load) = &port_load {
            load.add_to_pencil(&mut a, n, cx)?;
        }
        if let Some(load) = &dynamic_load {
            load.add_to_pencil(&mut a, n, body_n, cx)?;
        }
        if a.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid(
                "original hereditary coefficient is not representable",
            ));
        }
        let mut result = PreparedViscoPoles {
            source: self,
            sites,
            coefficients: a,
            row_scales: vec![1.0; n],
            column_scales: vec![1.0; n],
            dimension: n,
            band,
            budget,
            port_load,
            dynamic_load,
        };
        let rate = budget.reference_rate_per_s;
        for i in 0..n {
            checkpoint(cx)?;
            let scale = (0..n)
                .map(|j| result.envelope(i * n + j, rate))
                .sum::<f64>();
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("invalid reference row scale"));
            }
            result.row_scales[i] = scale;
        }
        for j in 0..n {
            let col = (0..n)
                .map(|i| result.envelope(i * n + j, rate) / result.row_scales[i])
                .fold(0.0_f64, f64::max);
            let scale = 1.0 / col;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("invalid reference column scale"));
            }
            result.column_scales[j] = scale;
        }
        for i in 0..n {
            let scale = (0..n)
                .map(|j| result.envelope(i * n + j, rate) * result.column_scales[j])
                .sum::<f64>();
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("invalid balanced reference row"));
            }
            result.row_scales[i] = scale;
        }
        checkpoint(cx)?;
        Ok(result)
    }
}

impl PreparedViscoPoles<'_> {
    /// Immutable original source; no fitted or averaged material replaces it.
    pub fn source(&self) -> &ViscoElements {
        self.source
    }
    /// Total free displacement, original memory and independent load coordinates.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Original full nodal DOF indices for the leading displacement block.
    pub fn free_dofs(&self) -> &[usize] {
        &self.source.equilibrium.free_dofs
    }
    /// Element/point/branch order of the subsequent six-component blocks.
    pub fn memory_sites(&self) -> &[ViscoMemorySite] {
        &self.sites
    }
    /// Original independent load-state suffix; no load state is condensed out.
    pub fn dynamic_load_range(&self) -> Option<core::ops::Range<usize>> {
        self.dynamic_load
            .as_ref()
            .map(|_| self.body_dimension()..self.dimension)
    }

    fn body_dimension(&self) -> usize {
        self.free_dofs().len() + 6 * self.sites.len()
    }

    fn envelope(&self, at: usize, magnitude: f64) -> f64 {
        self.coefficients[0][at].abs()
            + magnitude * self.coefficients[1][at].abs()
            + magnitude * magnitude * self.coefficients[2][at].abs()
    }
    fn argument(
        &self,
        s: C64,
        count: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), ViscoPoleError> {
        poll(cancelled)?;
        if !finite(s) || !s.abs().is_finite() {
            return Err(invalid("finite representable pole argument required"));
        }
        if s.abs() < self.band.0 || s.abs() > self.band.1 {
            return Err(ViscoPoleError::SourceBand {
                rate: s,
                band_rad_s: self.band,
            });
        }
        cap(
            "operator evaluations",
            *count as u128 + 1,
            self.budget.maximum_evaluations as u128,
        )?;
        *count += 1;
        Ok(())
    }
    fn evaluate(
        &self,
        s: C64,
        count: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, ViscoPoleError> {
        self.argument(s, count, cancelled)?;
        let n = self.dimension;
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        let square = s * s;
        let mut residual_scale = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut row = 0.0;
            for j in 0..n {
                let k = i * n + j;
                let scale = self.column_scales[j] / self.row_scales[i];
                matrix[k] = (C64::from_re(self.coefficients[0][k])
                    + s.scale(self.coefficients[1][k])
                    + square.scale(self.coefficients[2][k]))
                .scale(scale);
                derivative[k] = (C64::from_re(self.coefficients[1][k])
                    + s.scale(2.0 * self.coefficients[2][k]))
                .scale(scale);
                row += self.envelope(k, s.abs()) * scale;
            }
            residual_scale = residual_scale.max(row);
        }
        if !residual_scale.is_finite()
            || residual_scale <= 0.0
            || matrix.iter().chain(&derivative).any(|&z| !finite(z))
        {
            return Err(invalid("hereditary pencil evaluation is not representable"));
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        })
    }
    fn observe(
        &self,
        pole: NonlinearPole,
        count: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ViscoPole, ViscoPoleError> {
        let s = pole.value;
        self.argument(s, count, cancelled)?;
        let n = self.dimension;
        let vector: Vec<_> = pole
            .vector
            .iter()
            .zip(&self.column_scales)
            .map(|(&v, &c)| v.scale(c))
            .collect();
        if vector.iter().any(|&z| !finite(z)) {
            return Err(invalid("original pole vector is not representable"));
        }
        let mut numerators = [0.0_f64; 3];
        let mut denominators = [0.0_f64; 3];
        for i in 0..n {
            poll(cancelled)?;
            let mut sum = C64::ZERO;
            let mut correction = C64::ZERO;
            let mut envelope = 0.0;
            for (j, &v) in vector.iter().enumerate() {
                let k = i * n + j;
                let value = C64::from_re(self.coefficients[0][k])
                    + s.scale(self.coefficients[1][k])
                    + (s * s).scale(self.coefficients[2][k]);
                let term = value * v - correction;
                let next = sum + term;
                correction = (next - sum) - term;
                sum = next;
                envelope += self.envelope(k, s.abs()) * v.abs();
            }
            if !finite(sum) || !envelope.is_finite() {
                return Err(invalid("original pole residual overflow"));
            }
            let group = if i < self.free_dofs().len() {
                0
            } else if i < self.body_dimension() {
                1
            } else {
                2
            };
            numerators[group] = numerators[group].max(sum.abs());
            denominators[group] = denominators[group].max(envelope);
        }
        let residual: [f64; 3] = core::array::from_fn(|i| {
            if denominators[i] == 0.0 {
                if numerators[i] == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                numerators[i] / denominators[i]
            }
        });
        let tolerance = self.budget.physical_residual_tolerance;
        if residual.iter().any(|&r| !r.is_finite() || r > tolerance) {
            return Err(ViscoPoleError::PhysicalResidual {
                mechanical: residual[0],
                memory: residual[1],
                dynamic_load: residual[2],
                tolerance,
            });
        }
        Ok(ViscoPole {
            rate_per_s: s,
            port_load: self
                .port_load
                .as_ref()
                .map(|load| load.observe(s, &vector[..self.free_dofs().len()], cancelled))
                .transpose()?,
            dynamic_load: self
                .dynamic_load
                .as_ref()
                .map(|load| load.observe(self, s, &vector, cancelled))
                .transpose()?,
            physical_vector: vector,
            mechanical_residual: residual[0],
            memory_residual: residual[1],
            dynamic_load_residual: residual[2],
            relative_step: pole.relative_step,
            iterations: pole.iterations,
        })
    }

    /// Refine a local simple pole from a physical displacement seed. The seed's
    /// Maxwell strains are initialized by the exact branch equation; a seed
    /// exactly at an uncoupled relaxation rate refuses. Subsequent iterations
    /// retain ALL strain coordinates and never divide by a material denominator.
    /// Independent load states start at zero; a load-dominated mode can instead
    /// use resonance_pole_from_state with an authored full-state seed.
    /// Modes are numerical original-equation roots, not damping enclosures or
    /// a guarantee that a pole is mechanically observable. Real roots relax.
    pub fn resonance_pole(
        &self,
        rate_seed: C64,
        displacement_seed: &[C64],
        options: &PoleOptions,
        cx: &Cx<'_>,
    ) -> Result<ViscoPole, ViscoPoleError> {
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        let m = self.free_dofs().len();
        let n = self.dimension;
        if !finite(rate_seed)
            || displacement_seed.len() != m
            || displacement_seed.iter().any(|&z| !finite(z))
            || displacement_seed.iter().all(|z| z.abs() == 0.0)
        {
            return Err(invalid(
                "finite nonzero original displacement seed required",
            ));
        }
        let mut seed = vec![C64::ZERO; n];
        seed[..m].copy_from_slice(displacement_seed);
        for i in m..self.body_dimension() {
            poll(&mut cancelled)?;
            let denominator = C64::ONE + rate_seed.scale(self.coefficients[1][i * n + i]);
            if denominator.abs() == 0.0 {
                return Err(ViscoPoleError::SeedAtRelaxationRate { site: (i - m) / 6 });
            }
            let sum = (0..m).fold(C64::ZERO, |v, j| {
                v + displacement_seed[j].scale(self.coefficients[1][i * n + j])
            });
            seed[i] = -(rate_seed * sum) / denominator;
        }
        self.resonance_pole_from_state(rate_seed, &seed, options, cx)
    }

    /// Refine a simple pole from the complete ORIGINAL physical state, including
    /// each Maxwell strain and independent load coordinate. This also admits
    /// mechanically dark load/memory modes without inventing a displacement.
    /// No material or bare-load denominator is used to initialize this seed.
    pub fn resonance_pole_from_state(
        &self,
        rate_seed: C64,
        physical_seed: &[C64],
        options: &PoleOptions,
        cx: &Cx<'_>,
    ) -> Result<ViscoPole, ViscoPoleError> {
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        let n = self.dimension;
        if !finite(rate_seed)
            || physical_seed.len() != n
            || physical_seed.iter().any(|&z| !finite(z))
            || physical_seed.iter().all(|z| z.abs() == 0.0)
        {
            return Err(invalid("finite nonzero complete physical seed required"));
        }
        let mut seed = physical_seed.to_vec();
        for (v, &scale) in seed.iter_mut().zip(&self.column_scales) {
            *v = C64::new(v.re / scale, v.im / scale);
        }
        let mut count = 0;
        let mut origin = None;
        let result = newton_pole_with_cancel(
            n,
            rate_seed,
            &seed,
            options,
            |s, cancel| {
                self.evaluate(s, &mut count, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: s }
                })
            },
            &mut cancelled,
        );
        if let Some(error) = origin {
            return Err(error);
        }
        self.observe(
            result.map_err(ViscoPoleError::Pole)?,
            &mut count,
            &mut cancelled,
        )
    }

    /// Discover simple poles inside an ellipse, using the original polynomial
    /// pencil through every material relaxation rate. The ENTIRE ellipse must
    /// fit the conservative source |s| band. Multiple/defective modes may refuse.
    /// Two-moment Beyn can resolve at most its probe rank: a quadratic pencil
    /// can have more poles than rows, so a wide contour may need subdivision.
    /// Numerical quadrature/count agreement is not a spectral completeness proof.
    pub fn resonance_spectrum(
        &self,
        contour: Ellipse,
        options: &BeynOptions,
        cx: &Cx<'_>,
    ) -> Result<ViscoPoleSpectrum, ViscoPoleError> {
        let radius = contour.real_radius.max(contour.imaginary_radius);
        if !finite(contour.center)
            || !contour.real_radius.is_finite()
            || contour.real_radius <= 0.0
            || !contour.imaginary_radius.is_finite()
            || contour.imaginary_radius <= 0.0
            || (contour.center.abs() - radius).max(0.0) < self.band.0
            || contour.center.abs() + radius > self.band.1
        {
            return Err(invalid(
                "complete finite positive contour must fit the source band",
            ));
        }
        let mut count = 0;
        let mut origin = None;
        let mut cancelled = || cx.checkpoint().is_err();
        let result = beyn_with_cancel(
            self.dimension,
            contour,
            options,
            |s, cancel| {
                self.evaluate(s, &mut count, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: s }
                })
            },
            &mut cancelled,
        );
        if let Some(error) = origin {
            return Err(error);
        }
        let BeynReport {
            poles,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
        } = result.map_err(ViscoPoleError::Contour)?;
        let mut physical = Vec::with_capacity(poles.len());
        for pole in poles {
            physical.push(self.observe(pole, &mut count, &mut cancelled)?);
        }
        Ok(ViscoPoleSpectrum {
            poles: physical,
            unrefined_rates_per_s: unrefined_values,
            argument_integrals: [argument_integral, coarse_argument_integral],
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
        })
    }
}
