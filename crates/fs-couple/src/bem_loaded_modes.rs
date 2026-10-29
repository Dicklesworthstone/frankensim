//! Exterior-fluid-loaded structural resonances in the original boundary equations.
//!
//! For exp(-i omega t), unit-modal-mass displacement x and outward pressure p:
//! ```text
//! [ K - omega^2 I - i omega C,       W^T area ] [x] = 0
//! [ -rho omega^2 B(omega/c) W,       A(omega/c)] [p]
//! ```
//! The existing BEM owner supplies A p = B dp/dn. Keeping pressure unknowns
//! avoids treating the meromorphic quotient A^-1 B as a holomorphic pencil.
//! The shared fs-la Newton/Beyn owners perform every eigenvalue iteration.
//! Fixed diagonal balancing preserves holomorphy; acceptance checks all of the
//! ORIGINAL mechanical and acoustic rows independently of that balancing.
//!
//! These are roots of the selected discrete boundary representation. Neither
//! Burton-Miller nor a small residual excludes fictitious complex BIE roots:
//! Zheng et al. (2015), doi:10.1016/j.enganabound.2015.04.014. Panel, source-mode
//! and formulation convergence remain physical validation axes. No fitted
//! radiation resistance, imposed modal Q, or second eigensolver is introduced.
//!
//! Optional Maxwell internal strains follow the existing fs-phs storage law:
//! H_j = (c_j.x-z_j)^2/2, tau_j dz_j/dt = c_j.x-z_j. The complete
//! state order is [displacement, boundary pressure, elastic memory e=c.x-z].
//! Keeping e rather than dividing by 1-i*omega*tau retains material relaxation
//! poles, avoids cancellation of a large instantaneous stiffness against z,
//! and keeps the pencil analytic at each uncoupled material relaxation rate.

use fs_bem::helmholtz::laplace::complex::{ComplexError, ComplexOptions, assemble_with_cancel};
use fs_bem::helmholtz::{Formulation, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynOptions, BeynReport, Ellipse, beyn_with_cancel};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, NonlinearPole, PoleEvaluation, PoleOptions,
    newton_pole_with_cancel,
};
use fs_math::{c64::C64, det};

mod shell;
mod memory;

/// Explicit physical accuracy, discretization and callback resources.
#[derive(Debug, Clone, Copy)]
pub struct BemModesOptions {
    pub formulation: Formulation,
    pub boundary: ComplexOptions,
    /// Longest-edge bisections retaining each original panel's constant velocity.
    /// This refines BEM integration, not the source's independent spatial basis.
    pub panel_refinements: usize,
    /// Positive fixed balancing reference [rad/s], independent of contour nodes.
    pub reference_angular_frequency_rad_s: f64,
    /// Maximum componentwise backward error in EACH original physical row.
    pub physical_residual_tolerance: f64,
    /// Includes every structural coordinate, boundary pressure and material memory.
    pub max_dimension: usize,
    /// Conservative callback/model payload, excluding separately budgeted fs-la
    /// iteration, source geometry, allocator overhead and process RSS.
    pub max_dense_bytes: usize,
    /// BEM assemblies in one request, including seed and physical acceptance.
    /// Every assembly additionally obeys the boundary owner's work limit.
    pub max_operator_evaluations: usize,
    /// Per-assembly deterministic projection-work admission units:
    /// 8*panels^2*modes + 16*dimension^2, not a machine FLOP count.
    pub max_projection_work: u64,
}

#[derive(Debug)]
pub enum BemModesError {
    Invalid {
        what: &'static str,
    },
    Source {
        what: String,
    },
    Boundary(ComplexError),
    Pole(NonlinearEigError),
    Contour(BeynError),
    SeedSolve(LuComplexError),
    /// Only the local initializer eliminates memory. Contour discovery does not.
    SeedMemorySingular { branch: usize },
    EvaluationBudget {
        maximum: usize,
    },
    PhysicalResidual {
        mechanical: f64,
        boundary: f64,
        memory: f64,
    },
    /// A passive structural source does not justify publishing a growing mode.
    GrowingMode {
        angular_frequency_rad_s: C64,
    },
    Cancelled,
}
impl core::fmt::Display for BemModesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-BEM-MODES: {self:?}")
    }
}
impl std::error::Error for BemModesError {}

/// Original-coordinate eigenstate; arbitrary amplitude is chosen so maximum
/// displacement magnitude is one, or the maximum remaining coordinate is one
/// for a zero-motion mode. All fields refer to this SAME arbitrary amplitude.
#[derive(Debug, Clone)]
pub struct BemModePole {
    pub angular_frequency_rad_s: C64,
    pub modal_displacement: Vec<C64>,
    pub surface_pressure_pa: Vec<C64>,
    /// e_j=c_j.x-z_j in H_j=|e_j|^2/2 [sqrt(J)], in original branch order.
    pub elastic_memory_sqrt_j: Vec<C64>,
    pub mechanical_residual: f64,
    pub boundary_residual: f64,
    pub memory_residual: f64,
    pub relative_step: f64,
    pub iterations: usize,
}
impl BemModePole {
    pub fn decay_rate_per_s(&self) -> f64 {
        -self.angular_frequency_rad_s.im
    }
    /// Q=|Re omega|/(-2 Im omega). No finite Q is invented for undamped or
    /// nonoscillatory roots.
    pub fn quality_factor(&self) -> Option<f64> {
        let w = self.angular_frequency_rad_s;
        (w.im < 0.0 && w.re != 0.0).then(|| w.re.abs() / (-2.0 * w.im))
    }
}

#[derive(Debug, Clone)]
pub struct BemModeSpectrum {
    pub poles: Vec<BemModePole>,
    pub unrefined_values: Vec<C64>,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub argument_count: usize,
    pub rank: usize,
    pub nodes: usize,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
    pub operator_evaluations: usize,
}

/// Immutable source-derived mass-normalized mechanics and reciprocal skin.
pub struct BemLoadedModes {
    surface: SpherePanels,
    weights: Vec<Vec<f64>>,
    stiffness: Vec<f64>,
    viscosity: Vec<f64>,
    memory: memory::RetainedMemory,
    medium: Medium,
    settings: BemModesOptions,
    row_scales: Vec<f64>,
    column_scales: Vec<f64>,
}
struct PhysicalPencil {
    value: PoleEvaluation,
    // Coefficient magnitudes BEFORE cancellation of K-omega^2 I or BW.
    envelope: Vec<f64>,
}
fn bad(what: &'static str) -> BemModesError {
    BemModesError::Invalid { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), BemModesError> {
    if cancelled() {
        Err(BemModesError::Cancelled)
    } else {
        Ok(())
    }
}
fn divide(z: C64, s: f64) -> C64 {
    C64::new(z.re / s, z.im / s)
}

impl BemLoadedModes {
    // Only source adapters supply these mechanics; unrelated shape/frequency
    // lists are deliberately not exposed as a public physical constructor.
    fn build(
        surface: SpherePanels,
        weights: Vec<Vec<f64>>,
        stiffness: Vec<f64>,
        viscosity: Vec<f64>,
        medium: Medium,
        settings: BemModesOptions,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        poll(cancelled)?;
        let m = weights.len();
        let p = surface.centroids().len();
        Self::admit_size(m, p, &settings)?;
        if stiffness.len() != m * m
            || viscosity.len() != m * m
            || weights
                .iter()
                .any(|w| w.len() != p || w.iter().any(|x| !x.is_finite()))
            || stiffness.iter().chain(&viscosity).any(|x| !x.is_finite())
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(bad("source mechanics, skin or medium is invalid"));
        }
        let mut model = Self {
            surface,
            weights,
            stiffness,
            viscosity,
            memory: memory::RetainedMemory::default(),
            medium,
            settings,
            row_scales: Vec::new(),
            column_scales: Vec::new(),
        };
        model.rebalance(cancelled)?;
        Ok(model)
    }

    // Both the original model and a cold material attachment use one fixed
    // normalization. Source coordinates are unchanged by this solver scaling.
    fn rebalance(&mut self, cancelled: &mut dyn FnMut() -> bool) -> Result<(), BemModesError> {
        let n = self.dimension();
        self.row_scales = vec![1.0; n];
        self.column_scales = vec![1.0; n];
        // Fixed scaling preserves the analytic matrix function. Absolute-value
        // scaling recomputed at each contour node would violate that premise.
        let raw = self.assemble(
            C64::from_re(self.settings.reference_angular_frequency_rad_s),
            &mut 0,
            cancelled,
        )?;
        for i in 0..n {
            let row: f64 = raw.envelope[i * n..(i + 1) * n].iter().sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("zero or unrepresentable reference row"));
            }
            self.row_scales[i] = row;
        }
        for j in 0..n {
            let col = (0..n)
                .map(|i| raw.envelope[i * n + j] / self.row_scales[i])
                .fold(0.0, f64::max);
            let scale = 1.0 / col;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(bad("zero or unrepresentable reference coordinate"));
            }
            self.column_scales[j] = scale;
        }
        for i in 0..n {
            let row: f64 = (0..n)
                .map(|j| raw.envelope[i * n + j] * self.column_scales[j])
                .sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("unrepresentable balanced reference row"));
            }
            self.row_scales[i] = row;
        }
        poll(cancelled)?;
        Ok(())
    }

    fn admit_size(m: usize, p: usize, settings: &BemModesOptions) -> Result<(), BemModesError> {
        Self::admit_dimensions(m, p, 0, settings)
    }

    fn admit_dimensions(m: usize, p: usize, memories: usize, settings: &BemModesOptions) -> Result<(), BemModesError> {
        let n = m
            .checked_add(p)
            .and_then(|d| d.checked_add(memories))
            .ok_or_else(|| bad("coupled dimension overflow"))?;
        let bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(224))
            .and_then(|v| n.checked_mul(512).and_then(|s| v.checked_add(s)))
            .ok_or_else(|| bad("coupled payload overflow"))?;
        let work = (p as u64)
            .checked_mul(p as u64)
            .and_then(|v| v.checked_mul(m as u64))
            .and_then(|v| v.checked_mul(8))
            .and_then(|v| {
                (n as u64)
                    .checked_mul(n as u64)
                    .and_then(|s| s.checked_mul(16))
                    .and_then(|s| v.checked_add(s))
            })
            .ok_or_else(|| bad("coupled projection work overflow"))?;
        if m == 0
            || p == 0
            || n > MAX_NONLINEAR_POLE_DIMENSION
            || n > settings.max_dimension
            || bytes > settings.max_dense_bytes
            || work > settings.max_projection_work
            || settings.max_operator_evaluations == 0
            || !settings.reference_angular_frequency_rad_s.is_finite()
            || settings.reference_angular_frequency_rad_s <= 0.0
            || !settings.physical_residual_tolerance.is_finite()
            || settings.physical_residual_tolerance <= 0.0
            || settings.physical_residual_tolerance >= 1.0
        {
            return Err(bad(
                "invalid coupled limits or insufficient dimension, payload or work budget",
            ));
        }
        Ok(())
    }

    pub fn mode_count(&self) -> usize {
        self.weights.len()
    }
    pub fn panel_count(&self) -> usize {
        self.surface.centroids().len()
    }
    pub fn dimension(&self) -> usize {
        self.mode_count() + self.panel_count() + self.memory_count()
    }
    pub fn memory_count(&self) -> usize {
        self.memory.times.len()
    }
    pub fn surface(&self) -> &SpherePanels {
        &self.surface
    }
    pub fn normal_velocity_weights(&self) -> &[Vec<f64>] {
        &self.weights
    }
    pub fn stiffness_per_s2(&self) -> &[f64] {
        &self.stiffness
    }
    pub fn viscosity_per_s(&self) -> &[f64] {
        &self.viscosity
    }
    pub fn medium(&self) -> Medium {
        self.medium
    }

    fn assemble(
        &self,
        omega: C64,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PhysicalPencil, BemModesError> {
        poll(cancelled)?;
        if !finite(omega) || omega.abs() == 0.0 {
            return Err(bad("nonzero finite angular frequency required"));
        }
        if *evaluations >= self.settings.max_operator_evaluations {
            return Err(BemModesError::EvaluationBudget {
                maximum: self.settings.max_operator_evaluations,
            });
        }
        *evaluations += 1;
        let op = assemble_with_cancel(
            &self.surface,
            divide(omega, self.medium.sound_speed),
            self.settings.formulation,
            &self.settings.boundary,
            &mut *cancelled,
        )
        .map_err(BemModesError::Boundary)?;
        let (m, p, n) = (self.mode_count(), self.panel_count(), self.dimension());
        if op.dimension != p
            || [op.a.len(), op.b.len(), op.da_dk.len(), op.db_dk.len()]
                .iter()
                .any(|&s| s != p * p)
        {
            return Err(bad("BEM operator cardinality mismatch"));
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        let mut envelope = vec![0.0; n * n];
        let minus_i_omega = C64::new(omega.im, -omega.re);
        let w2 = omega * omega;
        let aw = omega.abs();
        for i in 0..m {
            poll(cancelled)?;
            for j in 0..m {
                let k = self.stiffness[i * m + j];
                let c = self.viscosity[i * m + j];
                matrix[i * n + j] = C64::from_re(k) + minus_i_omega.scale(c);
                derivative[i * n + j] = C64::new(0.0, -c);
                envelope[i * n + j] = k.abs() + aw * c.abs();
                if i == j {
                    matrix[i * n + j] = matrix[i * n + j] - w2;
                    derivative[i * n + j] = derivative[i * n + j] - omega.scale(2.0);
                    envelope[i * n + j] += aw * aw;
                }
            }
            for j in 0..p {
                let force = self.weights[i][j] * self.surface.areas()[j];
                matrix[i * n + m + j] = C64::from_re(force);
                envelope[i * n + m + j] = force.abs();
            }
        }
        for i in 0..p {
            poll(cancelled)?;
            for j in 0..m {
                let mut bw = C64::ZERO;
                let mut dbw = C64::ZERO;
                let mut magnitude = 0.0;
                for l in 0..p {
                    let weight = self.weights[j][l];
                    bw = bw + op.b[i * p + l].scale(weight);
                    dbw = dbw + op.db_dk[i * p + l].scale(weight / self.medium.sound_speed);
                    magnitude += op.b[i * p + l].abs() * weight.abs();
                }
                matrix[(m + i) * n + j] = (w2 * bw).scale(-self.medium.density);
                derivative[(m + i) * n + j] =
                    (omega.scale(2.0) * bw + w2 * dbw).scale(-self.medium.density);
                envelope[(m + i) * n + j] = self.medium.density * aw * aw * magnitude;
            }
            for j in 0..p {
                matrix[(m + i) * n + m + j] = op.a[i * p + j];
                derivative[(m + i) * n + m + j] =
                    divide(op.da_dk[i * p + j], self.medium.sound_speed);
                envelope[(m + i) * n + m + j] = op.a[i * p + j].abs();
            }
        }
        self.memory.add_to_pencil(omega, m, p, &mut matrix, &mut derivative, &mut envelope, cancelled)?;
        if matrix.iter().chain(&derivative).any(|&z| !finite(z))
            || envelope.iter().any(|x| !x.is_finite())
        {
            return Err(bad("coupled operator or derivative overflow"));
        }
        let scale = envelope
            .chunks(n)
            .map(|r| r.iter().sum::<f64>())
            .fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(bad("coupled coefficient scale is invalid"));
        }
        Ok(PhysicalPencil {
            value: PoleEvaluation {
                matrix,
                derivative,
                residual_scale: scale,
            },
            envelope,
        })
    }

    /// Original SI-coordinate matrix and analytic derivative for independent
    /// residual checks, perturbation studies and external forcing.
    pub fn original_pencil(
        &self,
        omega: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoleEvaluation, BemModesError> {
        Ok(self.assemble(omega, &mut 0, &mut cancelled)?.value)
    }

    fn scaled(
        &self,
        omega: C64,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, BemModesError> {
        let mut raw = self.assemble(omega, evaluations, cancelled)?;
        let n = self.dimension();
        let mut scale = 0.0_f64;
        for i in 0..n {
            let mut row = 0.0;
            for j in 0..n {
                let k = i * n + j;
                let factor = self.column_scales[j] / self.row_scales[i];
                raw.value.matrix[k] = raw.value.matrix[k].scale(factor);
                raw.value.derivative[k] = raw.value.derivative[k].scale(factor);
                row += raw.envelope[k] * factor;
            }
            scale = scale.max(row);
        }
        if !scale.is_finite()
            || scale <= 0.0
            || raw
                .value
                .matrix
                .iter()
                .chain(&raw.value.derivative)
                .any(|&z| !finite(z))
        {
            return Err(bad("balanced coupled operator is not representable"));
        }
        raw.value.residual_scale = scale;
        Ok(raw.value)
    }

    fn observe(
        &self,
        pole: NonlinearPole,
        evaluations: &mut usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<BemModePole, BemModesError> {
        let (m, n) = (self.mode_count(), self.dimension());
        let mut vector: Vec<_> = pole
            .vector
            .iter()
            .zip(&self.column_scales)
            .map(|(&v, &s)| v.scale(s))
            .collect();
        let mut amplitude = vector[..m].iter().map(|v| v.abs()).fold(0.0, f64::max);
        if amplitude == 0.0 {
            amplitude = vector.iter().map(|v| v.abs()).fold(0.0, f64::max);
        }
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(bad("original eigenstate scale is not representable"));
        }
        for v in &mut vector {
            *v = divide(*v, amplitude);
        }
        if vector.iter().any(|&v| !finite(v)) {
            return Err(bad("original eigenstate overflow"));
        }
        let raw = self.assemble(pole.value, evaluations, cancelled)?;
        let mut residuals = [0.0_f64; 3];
        for i in 0..n {
            let mut sum = C64::ZERO;
            let mut correction = C64::ZERO;
            let mut denominator = 0.0;
            for j in 0..n {
                let term = raw.value.matrix[i * n + j] * vector[j] - correction;
                let next = sum + term;
                correction = (next - sum) - term;
                sum = next;
                denominator += raw.envelope[i * n + j] * vector[j].abs();
            }
            if !finite(sum) || !denominator.is_finite() {
                return Err(bad("physical residual contraction overflow"));
            }
            let relative = if denominator == 0.0 {
                if sum.abs() == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                sum.abs() / denominator
            };
            if !relative.is_finite() {
                return Err(bad("physical residual is not representable"));
            }
            let block = if i < m { 0 } else if i < m + self.panel_count() { 1 } else { 2 };
            residuals[block] = residuals[block].max(relative);
        }
        if residuals
            .iter()
            .any(|&r| r > self.settings.physical_residual_tolerance)
        {
            return Err(BemModesError::PhysicalResidual {
                mechanical: residuals[0],
                boundary: residuals[1],
                memory: residuals[2],
            });
        }
        if pole.value.im > 0.0 {
            return Err(BemModesError::GrowingMode {
                angular_frequency_rad_s: pole.value,
            });
        }
        poll(cancelled)?;
        let memories = vector.split_off(m + self.panel_count());
        let pressure = vector.split_off(m);
        Ok(BemModePole {
            angular_frequency_rad_s: pole.value,
            modal_displacement: vector,
            surface_pressure_pa: pressure,
            elastic_memory_sqrt_j: memories,
            mechanical_residual: residuals[0],
            boundary_residual: residuals[1],
            memory_residual: residuals[2],
            relative_step: pole.relative_step,
            iterations: pole.iterations,
        })
    }

    /// Local simple-pole refinement. Only the initializer solves pressure from
    /// the displacement seed; Newton retains both fields in the shared bordered
    /// solver. A singular boundary block at the seed refuses. Contour discovery
    /// does not need this initializer or an invertible isolated acoustic block.
    pub fn resonance_pole(
        &self,
        omega0: C64,
        modal_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BemModePole, BemModesError> {
        poll(&mut cancelled)?;
        let (m, p, n) = (self.mode_count(), self.panel_count(), self.dimension());
        if modal_seed.len() != m
            || modal_seed.iter().any(|&z| !finite(z))
            || modal_seed.iter().all(|z| z.abs() == 0.0)
        {
            return Err(bad("finite nonzero source-mode seed required"));
        }
        let mut evaluations = 0;
        let raw = self.assemble(omega0, &mut evaluations, &mut cancelled)?;
        let mut a = vec![C64::ZERO; p * p];
        let mut rhs = vec![C64::ZERO; p];
        for i in 0..p {
            for j in 0..p {
                a[i * p + j] = raw.value.matrix[(m + i) * n + m + j];
            }
            for j in 0..m {
                rhs[i] = rhs[i] - raw.value.matrix[(m + i) * n + j] * modal_seed[j];
            }
        }
        let lu = lu_complex_with_cancel(&a, p, &mut cancelled).map_err(BemModesError::SeedSolve)?;
        lu.solve_with_cancel(&mut rhs, &mut cancelled)
            .map_err(BemModesError::SeedSolve)?;
        let mut seed = modal_seed.to_vec();
        seed.extend(rhs);
        seed.extend(self.memory.seed(omega0, modal_seed, &mut cancelled)?);
        drop(raw);
        drop(lu);
        drop(a);
        for (v, &scale) in seed.iter_mut().zip(&self.column_scales) {
            *v = divide(*v, scale);
        }
        let mut origin = None;
        let result = newton_pole_with_cancel(
            n,
            omega0,
            &seed,
            options,
            |w, cancel| {
                self.scaled(w, &mut evaluations, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: w }
                })
            },
            &mut cancelled,
        );
        if let Some(error) = origin {
            return Err(error);
        }
        self.observe(
            result.map_err(BemModesError::Pole)?,
            &mut evaluations,
            &mut cancelled,
        )
    }

    /// Discover roots of the full retained boundary/structure pencil. The ellipse
    /// must exclude zero, where Burton-Miller has i/k. Numerical count/rank and
    /// quadrature agreement is neither an interval completeness proof nor a
    /// classification of physical versus fictitious BIE roots. No boundary or
    /// zero-opening internal coordinates are silently discarded.
    pub fn resonance_spectrum(
        &self,
        region: Ellipse,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BemModeSpectrum, BemModesError> {
        poll(&mut cancelled)?;
        let x = region.center.re / region.real_radius;
        let y = region.center.im / region.imaginary_radius;
        if !finite(region.center)
            || !region.real_radius.is_finite()
            || region.real_radius <= 0.0
            || !region.imaginary_radius.is_finite()
            || region.imaginary_radius <= 0.0
            || x * x + y * y <= 1.0
        {
            return Err(bad("finite positive contour must exclude zero"));
        }
        let mut evaluations = 0;
        let mut origin = None;
        let result = beyn_with_cancel(
            self.dimension(),
            region,
            options,
            |w, cancel| {
                self.scaled(w, &mut evaluations, cancel).map_err(|error| {
                    origin = Some(error);
                    NonlinearEigError::FunctionRefused { at: w }
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
        } = result.map_err(BemModesError::Contour)?;
        let mut observed = Vec::with_capacity(poles.len());
        for pole in poles {
            observed.push(self.observe(pole, &mut evaluations, &mut cancelled)?);
        }
        Ok(BemModeSpectrum {
            poles: observed,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
            operator_evaluations: evaluations,
        })
    }
}

#[cfg(test)]
mod tests;
