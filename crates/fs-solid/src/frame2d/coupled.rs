//! Continuous native frame resonance with a complete physical boundary load.
//!
//! A load takes conjugate flow `u_l = B^T v` and returns opposing effort
//! `y_l = G_l^T e_l`. The original equations are
//! `M v_dot + C v + K q + B y_l = f`, `q_dot = v`, and
//! `E_l x_dot = (J_l-R_l)e_l + G_l B^T v`.
//! Their port exchange cancels exactly in the continuous power balance.
//! Every load state and algebraic multiplier remains, including internal modes
//! with zero port output. No isolated impedance, massless constraint or material
//! memory is inverted to obtain a frequency-dependent fitted frame damping.

use super::linear::{FrameLinearError, FrameLinearPencil};
use super::{Frame2d, Rayleigh};
use crate::SolidError;
use fs_exec::Cx;
use fs_la::factor::{FactorError, cholesky};
use fs_math::c64::C64;
use fs_phs::harmonic::{
    BeynOptions, HarmonicBudget, HarmonicError, HarmonicFactor, HarmonicModel, HarmonicPole,
    HarmonicPoleError, HarmonicResponse, HarmonicSpectrum, HarmonicSpectrumError,
    MAX_HARMONIC_DIMENSION, PoleOptions, ResonanceRegion,
};
use fs_phs::{DescriptorPortHamiltonian, PhsError, PortHamiltonian, QuadraticStorage, Storage};

/// Existing physical load owner. Ordinary and constrained laws use one assembly.
#[derive(Clone, Copy)]
pub enum FramePortDynamics<'a> {
    /// All source coordinates are differential states.
    Ode(&'a PortHamiltonian),
    /// Preserve the source's differential prefix and trailing effort multipliers.
    Descriptor(&'a DescriptorPortHamiltonian),
}

impl FramePortDynamics<'_> {
    fn parts(&self) -> (usize, usize, usize, &[f64], &[f64], &[f64], &dyn Storage) {
        match self {
            Self::Ode(s) => {
                let (j, r, g) = s.structure();
                (
                    s.state_dim(),
                    s.state_dim(),
                    s.port_dim(),
                    j,
                    r,
                    g,
                    s.storage(),
                )
            }
            Self::Descriptor(s) => {
                let (j, r, g) = s.structure();
                (
                    s.state_dim(),
                    s.differential_dim(),
                    s.port_dim(),
                    j,
                    r,
                    g,
                    s.storage(),
                )
            }
        }
    }
}

/// Physical motion/effort coupling, in the source load's original port order.
#[derive(Clone, Copy)]
pub struct FrameDynamicPortLoad<'a> {
    /// Original row-major `frame.ndof() * load.port_dim()` force map, including
    /// fixed rows. Its transpose maps nodal velocities to conjugate load input.
    /// Areas and lever arms belong here once, with their physical signs/units.
    pub force_map: &'a [f64],
    /// Exact declared linear storage is required. Constitutive authority and
    /// complex-frequency applicability remain with the source owner.
    pub dynamics: FramePortDynamics<'a>,
}

/// Bounds for native assembly, physical coupling, and one AC factor/response.
/// Pole/contour queries have their existing separate numerical options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameCoupledControls {
    /// Original frame count including supported coordinates, before assembly.
    pub maximum_coordinates: usize,
    /// Physical load ports; no mutual entry is discarded to meet this limit.
    pub maximum_ports: usize,
    /// Complete `2 * free_frame_coordinates + load.state_dim()` count.
    pub maximum_states: usize,
    /// Conservative simultaneous numerical payload; excludes source metadata,
    /// caller inputs and allocator overhead. No state truncation on exhaustion.
    pub maximum_payload_bytes: usize,
    /// Conservative arithmetic envelope including bounded native admission.
    /// This is an operation budget, not a measured elapsed-time guarantee.
    pub maximum_work: u64,
    /// Original mass-inverse and forced-response componentwise residual gate.
    /// Homogeneous poles use original row coefficient norms times the complete
    /// state/effort norm, so zero modal components may contain roundoff.
    pub maximum_relative_residual: f64,
}

/// Physical, numerical, resource or cancellation refusal; no source is mutated.
#[derive(Debug)]
pub enum FrameCoupledError {
    Frame(SolidError),
    Mass(FactorError),
    Physical(PhsError),
    Harmonic(HarmonicError),
    Pole(HarmonicPoleError),
    Spectrum(HarmonicSpectrumError),
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    Residual {
        relative: f64,
        tolerance: f64,
    },
    Cancelled,
}

impl core::fmt::Display for FrameCoupledError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FRAME-COUPLED: {self:?}")
    }
}
impl std::error::Error for FrameCoupledError {}
impl From<FrameLinearError> for FrameCoupledError {
    fn from(e: FrameLinearError) -> Self {
        match e {
            FrameLinearError::Frame(e) => Self::Frame(e),
            FrameLinearError::Invalid(s) => Self::Invalid(s),
            FrameLinearError::Budget {
                quantity,
                requested,
                cap,
            } => Self::Budget {
                quantity,
                requested: requested as u128,
                cap: cap as u128,
            },
            FrameLinearError::Cancelled => Self::Cancelled,
        }
    }
}
impl From<HarmonicError> for FrameCoupledError {
    fn from(e: HarmonicError) -> Self {
        match e {
            HarmonicError::Cancelled => Self::Cancelled,
            e => Self::Harmonic(e),
        }
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameCoupledError> {
    cx.checkpoint().map_err(|_| FrameCoupledError::Cancelled)
}
fn limit(quantity: &'static str, requested: u128, cap: u128) -> Result<(), FrameCoupledError> {
    if requested > cap {
        Err(FrameCoupledError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn relative(residual: C64, scale: f64) -> Result<f64, FrameCoupledError> {
    if !finite(residual) || !scale.is_finite() {
        return Err(FrameCoupledError::Invalid(
            "original physical residual is not representable",
        ));
    }
    Ok(if scale == 0.0 {
        if residual == C64::ZERO {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        residual.abs() / scale
    })
}

/// Immutable original-coordinate model, prepared once for AC and pole queries.
///
/// State order is `[q, p=M*v, original load state]` with q/p in `free_dofs()`
/// order. The load's trailing multipliers remain algebraic. In every returned
/// `HarmonicPole` or `HarmonicResponse`, `effort[nf..2*nf]` is physical velocity;
/// `port_output` / `output` is the same velocity in external force-port order.
/// No normalization of a modal vector is a claim of unit physical energy.
pub struct PreparedFrameResonance {
    free_dofs: Vec<usize>,
    load_dimension: usize,
    load_differential_dimension: usize,
    ports: usize,
    model: HarmonicModel,
    mass: Vec<f64>,
    damping: Vec<f64>,
    stiffness: Vec<f64>,
    force_map: Vec<f64>,
    load_g: Vec<f64>,
    load_coupling_row_norms: Vec<f64>,
    material_band: Option<(f64, f64)>,
    tolerance: f64,
}

/// Reuses the shared complex factor for repeated external forces at one frequency.
pub struct FrameCoupledHarmonicFactor<'a> {
    owner: &'a PreparedFrameResonance,
    factor: HarmonicFactor<'a>,
    omega: f64,
}

impl Frame2d {
    /// Couple the full native linear pencil to a source's complete linear pHS.
    ///
    /// Only the same elastic, fixed-geometry frames as native harmonic response
    /// are admitted. Every free coordinate needs original positive inertia;
    /// stiffness may be semidefinite and rigid motion is retained. The native
    /// full material C and separately authored Rayleigh model each enter once.
    ///
    /// The canonical momentum chart uses the native full M inverse. It is
    /// prepared through diagonally equilibrated native Cholesky solves, then
    /// independently checked against every ORIGINAL M row after symmetric
    /// rounding. Every later response/pole is checked again in original q/v
    /// equations, including p=M*v. This does not diagonalize mass or damping.
    ///
    /// A load has flow input and opposing effort output in the map's conjugate
    /// units. A pure resistance or inertance can use a descriptor load; its
    /// algebraic rows are retained. No isolated-load invertibility is assumed.
    /// The source's exact quadratic declaration is required, never inferred by
    /// probing or substituting a fitted oscillator. State/history is not copied
    /// into a purported steady state; this API prepares a continuous linear law.
    ///
    /// Construction polls between bounded native mass/material/admission calls,
    /// per assembly member, inverse column and coupling row. Shared descriptor
    /// admission is indivisible within its explicit cubic-work envelope. This
    /// numerical path is not a stability, spectral-completeness or material-band
    /// certificate. No load state, rigid coordinate or multiplier is discarded.
    pub fn prepare_coupled_resonance(
        &self,
        rayleigh: Rayleigh,
        load: FrameDynamicPortLoad<'_>,
        controls: FrameCoupledControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedFrameResonance, FrameCoupledError> {
        poll(cx)?;
        let (nl, ndl, ports, jl, rl, gl, storage) = load.dynamics.parts();
        let n = self.ndof();
        let nf = self.fixed.iter().filter(|&&fixed| !fixed).count();
        let full = (nf as u128) * 2 + nl as u128;
        limit(
            "full frame coordinates",
            n as u128,
            controls.maximum_coordinates as u128,
        )?;
        limit("load ports", ports as u128, controls.maximum_ports as u128)?;
        limit(
            "coupled states",
            full,
            controls.maximum_states.min(MAX_HARMONIC_DIMENSION) as u128,
        )?;
        if nf == 0
            || ports == 0
            || ndl > nl
            || load.force_map.len() as u128 != n as u128 * ports as u128
            || jl.len() as u128 != (nl as u128).pow(2)
            || rl.len() != jl.len()
            || gl.len() as u128 != nl as u128 * ports as u128
            || load.force_map.iter().any(|x| !x.is_finite())
            || !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(FrameCoupledError::Invalid(
                "nonempty frame/ports, finite original map and residual tolerance in (0,1) required",
            ));
        }
        let native_bytes = (n as u128)
            .pow(2)
            .saturating_mul(512)
            .saturating_add(512 * n as u128)
            .saturating_add(2048 * self.elements.len() as u128);
        let coupled_bytes =
            512 * full.pow(2) + 512 * full + 32 * (n as u128 + nl as u128) * ports as u128;
        limit(
            "coupled numerical payload bytes",
            native_bytes.saturating_add(coupled_bytes),
            controls.maximum_payload_bytes as u128,
        )?;
        // Two at-most-60-sweep cyclic Jacobi admissions fit within 2048*d^3;
        // reserve native assembly, mass inversion and original-row observations.
        let work = 2048 * full.pow(3)
            + 512 * (n as u128).pow(3)
            + 64 * full.pow(2) * ports as u128
            + 4096 * self.elements.len() as u128;
        limit(
            "coupled arithmetic work",
            work,
            controls.maximum_work as u128,
        )?;
        let FrameLinearPencil {
            free_dofs,
            mass,
            stiffness,
            damping,
            material_omega_band_rad_s,
        } = self.prepare_linear_pencil(
            rayleigh,
            controls.maximum_coordinates,
            native_bytes as usize,
            cx,
        )?;
        let inverse = mass_inverse(&mass, nf, controls.maximum_relative_residual, cx)?;
        let part = storage.quadratic_part().ok_or(FrameCoupledError::Invalid(
            "load must explicitly declare its complete quadratic storage",
        ))?;
        if !part.exact
            || part.matrix.len() != nl * nl
            || part
                .matrix
                .iter()
                .chain(jl)
                .chain(rl)
                .chain(gl)
                .any(|x| !x.is_finite())
        {
            return Err(FrameCoupledError::Invalid(
                "load quadratic law and structure must be complete and finite",
            ));
        }
        for i in 0..nl {
            poll(cx)?;
            for j in 0..nl {
                if (i >= ndl || j >= ndl) && part.matrix[i * nl + j] != 0.0 {
                    return Err(FrameCoupledError::Invalid(
                        "load algebraic multipliers cannot carry stored energy",
                    ));
                }
            }
        }
        let d = full as usize;
        let dd = 2 * nf + ndl;
        let mut q = vec![0.0; dd * dd];
        let mut j = vec![0.0; d * d];
        let mut r = vec![0.0; d * d];
        let mut g = vec![0.0; d * nf];
        let mut force_map = Vec::with_capacity(nf * ports);
        let mut load_coupling_row_norms = vec![0.0; nf];
        for (i, &native_row) in free_dofs.iter().enumerate() {
            poll(cx)?;
            j[i * d + nf + i] = 1.0;
            j[(nf + i) * d + i] = -1.0;
            g[(nf + i) * nf + i] = 1.0;
            force_map
                .extend_from_slice(&load.force_map[native_row * ports..(native_row + 1) * ports]);
            for k in 0..nf {
                q[i * dd + k] = stiffness[i * nf + k];
                q[(nf + i) * dd + nf + k] = inverse[i * nf + k];
                r[(nf + i) * d + nf + k] = damping[i * nf + k];
            }
            for k in 0..nl {
                let coupling: f64 = (0..ports)
                    .map(|p| force_map[i * ports + p] * gl[k * ports + p])
                    .sum();
                j[(nf + i) * d + 2 * nf + k] = -coupling;
                j[(2 * nf + k) * d + nf + i] = coupling;
                load_coupling_row_norms[i] += coupling.abs();
            }
        }
        for i in 0..nl {
            poll(cx)?;
            for k in 0..nl {
                j[(2 * nf + i) * d + 2 * nf + k] = jl[i * nl + k];
                r[(2 * nf + i) * d + 2 * nf + k] = rl[i * nl + k];
                if i < ndl && k < ndl {
                    q[(2 * nf + i) * dd + 2 * nf + k] = part.matrix[i * nl + k];
                }
            }
        }
        if q.iter()
            .chain(&j)
            .chain(&r)
            .chain(&load_coupling_row_norms)
            .any(|x| !x.is_finite())
        {
            return Err(FrameCoupledError::Invalid(
                "coupled physical operator is not representable",
            ));
        }
        poll(cx)?;
        let storage = QuadraticStorage::new(q, dd).map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let system = DescriptorPortHamiltonian::new(d, dd, nf, j, r, g, Box::new(storage))
            .map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let model = HarmonicModel::from_descriptor(
            &system,
            HarmonicBudget {
                max_dimension: controls.maximum_states,
                max_dense_bytes: coupled_bytes as usize,
                max_dense_work: controls.maximum_work,
                relative_residual_tolerance: controls.maximum_relative_residual,
            },
            || cx.checkpoint().is_err(),
        )?;
        poll(cx)?;
        Ok(PreparedFrameResonance {
            free_dofs,
            load_dimension: nl,
            load_differential_dimension: ndl,
            ports,
            model,
            mass,
            damping,
            stiffness,
            force_map,
            load_g: gl.to_vec(),
            load_coupling_row_norms,
            material_band: material_omega_band_rad_s,
            tolerance: controls.maximum_relative_residual,
        })
    }
}

fn mass_inverse(
    m: &[f64],
    n: usize,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, FrameCoupledError> {
    let scales: Vec<_> = (0..n).map(|i| m[i * n + i].sqrt()).collect();
    if scales.iter().any(|&x| !x.is_finite() || x <= 0.0) {
        return Err(FrameCoupledError::Invalid(
            "positive representable native mass scales required",
        ));
    }
    let mut equilibrated = vec![0.0; n * n];
    for i in 0..n {
        poll(cx)?;
        for k in 0..n {
            equilibrated[i * n + k] = m[i * n + k] / scales[i] / scales[k];
        }
    }
    if equilibrated.iter().any(|x| !x.is_finite()) {
        return Err(FrameCoupledError::Invalid(
            "native mass equilibration is not representable",
        ));
    }
    let factor = cholesky(&equilibrated, n).map_err(FrameCoupledError::Mass)?;
    let mut inverse = vec![0.0; n * n];
    let mut column = vec![0.0; n];
    for k in 0..n {
        poll(cx)?;
        column.fill(0.0);
        column[k] = 1.0;
        factor.solve(&mut column);
        for i in 0..n {
            inverse[i * n + k] = column[i] / scales[i] / scales[k];
        }
    }
    for i in 0..n {
        poll(cx)?;
        for k in 0..i {
            let a = inverse[i * n + k];
            let b = inverse[k * n + i];
            let mean = if a == b { a } else { 0.5 * a + 0.5 * b };
            inverse[i * n + k] = mean;
            inverse[k * n + i] = mean;
        }
    }
    if inverse.iter().any(|x| !x.is_finite()) {
        return Err(FrameCoupledError::Invalid(
            "native momentum storage is not representable",
        ));
    }
    let mut error = 0.0_f64;
    for i in 0..n {
        poll(cx)?;
        for k in 0..n {
            let unit = if i == k { 1.0 } else { 0.0 };
            let mut residual = -unit;
            let mut scale = unit;
            for l in 0..n {
                residual = m[i * n + l].mul_add(inverse[l * n + k], residual);
                scale += m[i * n + l].abs() * inverse[l * n + k].abs();
            }
            error = error.max(relative(C64::from_re(residual), scale)?);
        }
    }
    if error > tolerance {
        return Err(FrameCoupledError::Residual {
            relative: error,
            tolerance,
        });
    }
    Ok(inverse)
}

impl PreparedFrameResonance {
    /// Original node-major [ux, uy, theta] map. Supported coordinates stay zero.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Full state count, including original load multipliers.
    pub fn dimension(&self) -> usize {
        2 * self.free_dofs.len() + self.load_dimension
    }
    /// The original load's complete state occupies this unchanged contiguous range.
    pub fn load_state_range(&self) -> core::ops::Range<usize> {
        2 * self.free_dofs.len()..self.dimension()
    }
    /// Retained algebraic efforts; no fictitious storage was added to these rows.
    pub fn load_algebraic_range(&self) -> core::ops::Range<usize> {
        2 * self.free_dofs.len() + self.load_differential_dimension..self.dimension()
    }
    /// Original material-viscosity real-frequency support. It is provenance for
    /// complex poles, not certification of analytic continuation outside it.
    pub fn material_omega_band_rad_s(&self) -> Option<(f64, f64)> {
        self.material_band
    }

    /// Refine one isolated pole of the complete coupled pencil under exp(s t).
    /// The seed uses this model's original state order, including all load states.
    /// Repeated/defective modes and unresolved precision can refuse. This is a
    /// local query, not a complete spectrum or transfer-observability claim.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        cx: &Cx<'_>,
    ) -> Result<HarmonicPole, FrameCoupledError> {
        poll(cx)?;
        let mut pole = self
            .model
            .resonance_pole(laplace_seed_per_s, state_seed, options, || {
                cx.checkpoint().is_err()
            })
            .map_err(FrameCoupledError::Pole)?;
        pole.relative_residual = pole.relative_residual.max(self.frame_residual(
            pole.laplace_pole_per_s,
            &pole.state,
            &pole.effort,
            None,
            cx,
        )?);
        poll(cx)?;
        Ok(pole)
    }

    /// Discover poles in a declared Laplace-plane ellipse with the existing
    /// shared contour owner. All states enter; only the requested region is
    /// searched. Rank/quadrature/count agreement is numerical evidence, not an
    /// interval completeness proof. Infinite descriptor roots are not returned.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        cx: &Cx<'_>,
    ) -> Result<HarmonicSpectrum, FrameCoupledError> {
        poll(cx)?;
        let mut spectrum = self
            .model
            .resonance_spectrum(region, options, || cx.checkpoint().is_err())
            .map_err(FrameCoupledError::Spectrum)?;
        for pole in &mut spectrum.poles {
            pole.relative_residual = pole.relative_residual.max(self.frame_residual(
                pole.laplace_pole_per_s,
                &pole.state,
                &pole.effort,
                None,
                cx,
            )?);
        }
        poll(cx)?;
        Ok(spectrum)
    }

    /// Reuse the same original coupled operator for AC forcing, with the native
    /// material's actual frequency gate. Load source-domain admission remains
    /// the load owner's responsibility. Phasors use exp(-i omega t).
    pub fn factor(
        &self,
        omega_rad_s: f64,
        cx: &Cx<'_>,
    ) -> Result<FrameCoupledHarmonicFactor<'_>, FrameCoupledError> {
        poll(cx)?;
        if let Some((lo, hi)) = self.material_band {
            if !omega_rad_s.is_finite() || omega_rad_s < lo || omega_rad_s > hi {
                return Err(FrameCoupledError::Invalid(
                    "forcing frequency is outside native material viscosity support",
                ));
            }
        }
        let factor = self
            .model
            .factor(omega_rad_s, || cx.checkpoint().is_err())?;
        Ok(FrameCoupledHarmonicFactor {
            owner: self,
            factor,
            omega: omega_rad_s,
        })
    }

    // Independent ORIGINAL q/v equations, not the assembled canonical generator.
    // Preserve the two original port maps in this observation as well as M/C/K.
    // For homogeneous modes use row coefficient envelopes in the augmented
    // variables [original state, physical v, original load efforts]. Their
    // full infinity norm handles dark modes without exact-zero requirements.
    // Forced AC retains the original componentwise force-relative denominator.
    fn frame_residual(
        &self,
        s: C64,
        state: &[C64],
        effort: &[C64],
        force: Option<&[C64]>,
        cx: &Cx<'_>,
    ) -> Result<f64, FrameCoupledError> {
        let n = self.free_dofs.len();
        let q = &state[..n];
        let v = &effort[n..2 * n];
        let homogeneous_norm = force.is_none().then(|| {
            state
                .iter()
                .chain(&effort[n..])
                .map(|z| z.abs())
                .fold(0.0_f64, f64::max)
        });
        let mut load_effort = vec![C64::ZERO; self.ports];
        for (p, value) in load_effort.iter_mut().enumerate() {
            poll(cx)?;
            for k in 0..self.load_dimension {
                *value = *value + effort[2 * n + k].scale(self.load_g[k * self.ports + p]);
            }
        }
        let mut maximum = 0.0_f64;
        for i in 0..n {
            poll(cx)?;
            let qdot = s * q[i];
            let kinematic_scale =
                homogeneous_norm.map_or(qdot.abs() + v[i].abs(), |norm| (s.abs() + 1.0) * norm);
            maximum = maximum.max(relative(qdot - v[i], kinematic_scale)?);
            let mut momentum = state[n + i];
            let mut momentum_scale = momentum.abs();
            let mut momentum_coefficients = 1.0;
            let applied = force.map_or(C64::ZERO, |f| f[i]);
            let mut row = -applied;
            let mut scale = applied.abs();
            let mut coefficients = self.load_coupling_row_norms[i];
            for k in 0..n {
                let mv = v[k].scale(self.mass[i * n + k]);
                momentum = momentum - mv;
                momentum_scale += mv.abs();
                momentum_coefficients += self.mass[i * n + k].abs();
                coefficients += s.abs() * self.mass[i * n + k].abs()
                    + self.damping[i * n + k].abs()
                    + self.stiffness[i * n + k].abs();
                for term in [
                    s * mv,
                    v[k].scale(self.damping[i * n + k]),
                    q[k].scale(self.stiffness[i * n + k]),
                ] {
                    row = row + term;
                    scale += term.abs();
                }
            }
            for (p, &value) in load_effort.iter().enumerate() {
                let term = value.scale(self.force_map[i * self.ports + p]);
                row = row + term;
                scale += term.abs();
            }
            if let Some(norm) = homogeneous_norm {
                scale = coefficients * norm;
                momentum_scale = momentum_coefficients * norm;
            }
            maximum = maximum
                .max(relative(momentum, momentum_scale)?)
                .max(relative(row, scale)?);
        }
        if maximum > self.tolerance {
            return Err(FrameCoupledError::Residual {
                relative: maximum,
                tolerance: self.tolerance,
            });
        }
        Ok(maximum)
    }
}

impl FrameCoupledHarmonicFactor<'_> {
    /// Peak nodal forces/moments in `free_dofs()` order. Shared physical storage
    /// and dissipation include every frame and load state exactly once. Port
    /// exchange is internal; net loss is not automatically acoustic radiation
    /// or deposited material heat. The factor changes no accepted trajectory.
    pub fn response(
        &self,
        forces: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, FrameCoupledError> {
        poll(cx)?;
        let mut response = self.factor.response(forces, || cx.checkpoint().is_err())?;
        response.relative_residual = response.relative_residual.max(self.owner.frame_residual(
            C64::new(0.0, -self.omega),
            &response.state,
            &response.effort,
            Some(forces),
            cx,
        )?);
        poll(cx)?;
        Ok(response)
    }
}
