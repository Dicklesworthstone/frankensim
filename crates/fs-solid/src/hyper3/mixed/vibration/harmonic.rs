//! Full forced response about one accepted mixed finite-strain equilibrium.
//!
//! Pressure has no inertia. Finite compliance uses the original Schur pencil;
//! the exact constraint uses its COMPLETE nullspace, without modal truncation.
//! Its pressure is recovered from the full forced momentum equation, including
//! force components blocked by incompressibility. The same original K and M
//! supply all displacement rows and support reactions.
//!
//! The response is a lossless, frozen-equilibrium linear perturbation under
//! exp(-i omega t). It is not a finite-amplitude constitutive admission, startup
//! trajectory, evolving thermal/history solve or damping model. Prescribed
//! support increments lift BOTH momentum and the original weak pressure
//! equations. Singular lossless poles and incompatible volume motion refuse.

use fs_exec::Cx;
use fs_la::factor::{Lu, Qr, lu, qr};
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::MixedFiniteStrainPotential;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponse, HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_sparse::{Csr, LdltFactor};

use super::super::{MixedUpError, checkpoint};
use super::MixedUpVibration;
use super::incompressible::IncompressibleMixedVibration;

#[derive(Debug)]
pub enum MixedHarmonicError {
    Mixed(MixedUpError),
    Harmonic(HarmonicResponseError),
    Residual {
        equation: &'static str,
        relative: f64,
        limit: f64,
    },
}
impl From<MixedUpError> for MixedHarmonicError {
    fn from(error: MixedUpError) -> Self {
        Self::Mixed(error)
    }
}
impl From<HarmonicResponseError> for MixedHarmonicError {
    fn from(error: HarmonicResponseError) -> Self {
        Self::Harmonic(error)
    }
}
impl core::fmt::Display for MixedHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mixed solid harmonic response: {self:?}")
    }
}
impl core::error::Error for MixedHarmonicError {}
type Result<T> = core::result::Result<T, MixedHarmonicError>;

fn invalid(message: &'static str) -> MixedHarmonicError {
    HarmonicResponseError::Invalid(message).into()
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

/// Peak perturbations about the retained static equilibrium, in original units.
#[derive(Clone, Debug)]
pub struct MixedHarmonicResponse {
    /// Original node-major xyz ordering, including prescribed support increments.
    pub displacement_m: Vec<C64>,
    pub velocity_m_s: Vec<C64>,
    /// Original equilibrium.pressure_nodes order; the original gauge has zero
    /// increment. This includes pressure excited without any displacement.
    pub pressure_pa: Vec<C64>,
    /// Force exerted by each prescribed support ON the solid:
    /// (Kuu-omega² M) du + Kup dp - authored_force, in original DOF order.
    /// An authored force at a support is subtracted exactly once.
    pub support_reaction_n: Vec<(usize, C64)>,
    /// Original free momentum and ALL weak pressure equations, respectively.
    /// Each block uses max(abs(residual)) / max(row sum of abs(terms)). The
    /// force and volume blocks have separate physical scales, without a unit
    /// floor; rows whose physical action vanishes do not amplify roundoff.
    pub momentum_relative_residual: f64,
    pub pressure_relative_residual: f64,
    /// Re(F^H v)/2 from the original full force and velocity [W]. This is
    /// signed numerical work; no dissipation or material heat is inferred.
    pub average_input_power_w: f64,
    /// Re(R_support^H v_support)/2 [W], signed work supplied BY the supports
    /// ON the solid. Applied support-coordinate forces are counted above.
    pub average_support_power_w: f64,
    /// Applied plus support power [W]. This lossless model has no heat term;
    /// signed numerical imbalance is retained, never reassigned to damping.
    pub power_balance_residual_w: f64,
    /// Native solve in free-displacement coordinates for finite compliance,
    /// or complete Z coordinates for the exact constraint. None means Z has
    /// zero columns: only the prescribed affine constraint lift can move.
    /// With moving supports this describes the homogeneous correction, whose
    /// force and power are not the complete original physical work account.
    pub projected_response: Option<HarmonicResponse>,
}

#[derive(Clone, Copy)]
struct PhysicalPencil<'a> {
    free: &'a [usize],
    tangent: &'a Csr,
    mass: &'a Csr,
    pressure_tolerance: f64,
}
#[derive(Clone, Copy)]
enum Reduction<'a> {
    Finite {
        pressure_response: &'a [f64],
        pressure_factor: &'a LdltFactor,
    },
    Exact {
        basis: &'a [f64],
        pressure_indices: &'a [usize],
    },
}

struct PressureRecovery {
    qr: Qr,
    // Native LU of QR's small upper factor, used for R^T y=b. No normal
    // equations, pressure penalty or separately invented triangular solver.
    upper: Lu,
    scale: f64,
}

/// Reusable native factor and pressure recovery at one physical frequency.
/// The immutable original solid supplies every retained coefficient.
pub struct PreparedMixedHarmonic<'a> {
    physical: PhysicalPencil<'a>,
    reduction: Reduction<'a>,
    factor: Option<PreparedHarmonicResponse>,
    pressure_factor: Option<PressureRecovery>,
    omega: f64,
    coordinates: usize,
    momentum_tolerance: f64,
}

impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedUpVibration<'_, L> {
    /// Prepare the complete original finite-compliance forced response.
    ///
    /// Controls count original FREE displacement DOFs before any reduction.
    /// maximum_payload_bytes covers this adapter plus its native harmonic
    /// preparation and one response; borrowed source operators are excluded.
    /// Forces passed to solve use ALL original displacement DOFs, never pressure
    /// coordinates or modal forces. No damping or additional prestress is added.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedMixedHarmonic<'_>> {
        prepare(
            PhysicalPencil {
                free: &self.free_dofs,
                tangent: &self.mixed_tangent,
                mass: self.reference_mass.operator.matrix(),
                pressure_tolerance: self.settings.pressure_residual_tolerance,
            },
            Reduction::Finite {
                pressure_response: &self.pressure_response_pa_per_m,
                pressure_factor: &self.pressure_factor,
            },
            &self.stiffness_n_m,
            &self.mass_kg,
            omega_rad_s,
            controls,
            cx,
        )
    }
}
impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> IncompressibleMixedVibration<'_, L> {
    /// Prepare all original exact-constraint coordinates, including a zero-
    /// dimensional admissible space. The original pressure gauge and full weak
    /// constraints remain in force. Controls and force ordering match the
    /// finite-compliance method; no pressure mass or bulk penalty is introduced.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedMixedHarmonic<'_>> {
        prepare(
            PhysicalPencil {
                free: self.free_dofs(),
                tangent: self.mixed_tangent(),
                mass: self.reference_mass_matrix(),
                pressure_tolerance: self.settings().pressure_residual_tolerance,
            },
            Reduction::Exact {
                basis: self.basis(),
                pressure_indices: self.independent_pressure_indices(),
            },
            self.projected_stiffness_n_m(),
            self.projected_mass_kg(),
            omega_rad_s,
            controls,
            cx,
        )
    }
}

fn prepare<'a>(
    physical: PhysicalPencil<'a>,
    reduction: Reduction<'a>,
    stiffness: &[f64],
    mass: &[f64],
    omega: f64,
    controls: HarmonicResponseControls,
    cx: &Cx<'_>,
) -> Result<PreparedMixedHarmonic<'a>> {
    checkpoint(cx)?;
    if !omega.is_finite()
        || omega < 0.0
        || !(omega * omega).is_finite()
        || !controls.maximum_relative_residual.is_finite()
        || controls.maximum_relative_residual <= 0.0
        || controls.maximum_relative_residual >= 1.0
    {
        return Err(invalid(
            "finite representable frequency and residual budget in (0,1) required",
        ));
    }
    let n = physical.free.len();
    let nu = physical.mass.nrows();
    let np = physical.tangent.nrows() - nu;
    let m = match reduction {
        Reduction::Finite { .. } => 0,
        Reduction::Exact {
            pressure_indices, ..
        } => pressure_indices.len(),
    };
    let r = n - m;
    if n > controls.maximum_coordinates {
        return Err(HarmonicResponseError::Budget {
            quantity: "original free displacement coordinates",
            requested: n,
            cap: controls.maximum_coordinates,
        }
        .into());
    }
    // Native dense factor + solve: 144*r²+128*r bytes. Extra payload covers
    // zero damping input, QR input/factor/scratch, native R factor/solve,
    // affine lifting, full physical outputs and all reaction pairs. QR/LU are bounded by
    // the original <=128-free-displacement owner and polled around its call.
    let extra = r
        .checked_mul(r)
        .and_then(|v| v.checked_mul(8))
        .and_then(|v| {
            n.checked_mul(m)
                .and_then(|q| q.checked_mul(24))
                .and_then(|q| v.checked_add(q))
        })
        .and_then(|v| {
            m.checked_mul(m)
                .and_then(|q| q.checked_mul(48))
                .and_then(|q| v.checked_add(q))
        })
        .and_then(|v| {
            nu.checked_add(np)
                .and_then(|q| q.checked_add(n))
                .and_then(|q| q.checked_mul(256))
                .and_then(|q| v.checked_add(q))
        })
        .and_then(|v| v.checked_add(4096))
        .ok_or_else(|| invalid("mixed harmonic payload extent overflow"))?;
    let native = r
        .checked_mul(r)
        .and_then(|v| v.checked_mul(144))
        .and_then(|v| r.checked_mul(128).and_then(|q| v.checked_add(q)))
        .ok_or_else(|| invalid("native harmonic payload extent overflow"))?;
    let bytes = extra
        .checked_add(native)
        .ok_or_else(|| invalid("combined harmonic payload extent overflow"))?;
    if bytes > controls.maximum_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "combined mixed harmonic payload bytes",
            requested: bytes,
            cap: controls.maximum_payload_bytes,
        }
        .into());
    }
    let pressure_factor = if let Reduction::Exact {
        pressure_indices, ..
    } = reduction
    {
        let mut transpose = Vec::with_capacity(n * m);
        for &dof in physical.free {
            checkpoint(cx)?;
            for &p in pressure_indices {
                transpose.push(physical.tangent.get(dof, nu + p));
            }
        }
        let scale = transpose.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale == 0.0 {
            return Err(invalid("unrepresentable original pressure constraint"));
        }
        for value in &mut transpose {
            *value /= scale;
        }
        let factor = qr(&transpose, n, m);
        checkpoint(cx)?;
        for i in 0..m {
            if !factor.r(i, i).is_finite() || factor.r(i, i) == 0.0 {
                return Err(invalid("original pressure recovery lost rank"));
            }
        }
        let mut upper = Vec::with_capacity(m * m);
        for i in 0..m {
            checkpoint(cx)?;
            for j in 0..m {
                upper.push(factor.r(i, j));
            }
        }
        let upper = lu(&upper, m).map_err(|_| invalid("original pressure lift lost rank"))?;
        checkpoint(cx)?;
        Some(PressureRecovery {
            qr: factor,
            upper,
            scale,
        })
    } else {
        None
    };
    let factor = if r == 0 {
        None
    } else {
        Some(PreparedHarmonicResponse::try_new(
            mass,
            &vec![0.0; r * r],
            stiffness,
            r,
            omega,
            HarmonicResponseControls {
                maximum_payload_bytes: controls.maximum_payload_bytes - extra,
                ..controls
            },
            cx,
        )?)
    };
    checkpoint(cx)?;
    Ok(PreparedMixedHarmonic {
        physical,
        reduction,
        factor,
        pressure_factor,
        omega,
        coordinates: r,
        momentum_tolerance: controls.maximum_relative_residual,
    })
}

impl PreparedMixedHarmonic<'_> {
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Solve a peak force phasor in ALL original node-major xyz coordinates.
    /// Fixed support entries influence reactions but have no admissible motion.
    /// The original equilibrium, materials, quadrature, supports and factors
    /// remain immutable; no partial response is published after a refusal.
    pub fn solve(&self, force_n: &[C64], cx: &Cx<'_>) -> Result<MixedHarmonicResponse> {
        self.solve_with_support_motion(force_n, &[], cx)
    }

    /// Force and base excitation about the SAME retained equilibrium.
    ///
    /// `support_motion_m` contains strictly increasing original displacement
    /// DOFs which are already prescribed by the source problem. Omitted
    /// supports have zero increment. Entries are peak displacement phasors
    /// [m], not absolute static positions, velocities or accelerations.
    ///
    /// Finite bulk response retains H^-1 B_support g. Exact incompressibility
    /// first lifts B_free u0 = -B_support g, then solves in the COMPLETE
    /// homogeneous nullspace. All pressure equations, including a removed
    /// gauge row, must remain compatible. Full K and M act on this lift; a
    /// static support load or projection of force alone is insufficient.
    pub fn solve_with_support_motion(
        &self,
        force_n: &[C64],
        support_motion_m: &[(usize, C64)],
        cx: &Cx<'_>,
    ) -> Result<MixedHarmonicResponse> {
        checkpoint(cx)?;
        let PhysicalPencil {
            free,
            tangent,
            mass,
            pressure_tolerance,
        } = self.physical;
        let nu = mass.nrows();
        let np = tangent.nrows() - nu;
        let n = free.len();
        let r = self.coordinates;
        if force_n.len() != nu {
            return Err(invalid(
                "one force phasor per ORIGINAL displacement DOF required",
            ));
        }
        for chunk in force_n.chunks(128) {
            checkpoint(cx)?;
            if chunk.iter().any(|&v| !finite(v)) {
                return Err(invalid("nonfinite force phasor"));
            }
        }
        let mut previous = None;
        for &(dof, value) in support_motion_m {
            checkpoint(cx)?;
            if dof >= nu
                || free.binary_search(&dof).is_ok()
                || previous.is_some_and(|p| p >= dof)
                || !finite(value)
            {
                return Err(invalid(
                    "finite support phasors on strictly increasing original prescribed DOFs required",
                ));
            }
            previous = Some(dof);
        }
        let mut displacement_m = vec![C64::ZERO; nu];
        for &(dof, value) in support_motion_m {
            displacement_m[dof] = value;
        }
        let mut pressure_pa = vec![C64::ZERO; np];
        if !support_motion_m.is_empty() {
            match self.reduction {
                Reduction::Finite {
                    pressure_factor, ..
                } => {
                    // B_support g is a pressure-equation RHS [m^3]. The
                    // original positive H=-Kpp maps it to pressure [Pa].
                    let mut real = vec![0.0; np];
                    let mut imaginary = vec![0.0; np];
                    for p in 0..np {
                        checkpoint(cx)?;
                        let mut rhs = C64::ZERO;
                        let (columns, values) = tangent.row(nu + p);
                        for (&dof, &value) in columns.iter().zip(values) {
                            if dof < nu {
                                rhs = rhs + displacement_m[dof].scale(value);
                            }
                        }
                        if !finite(rhs) {
                            return Err(invalid("support pressure right hand side overflow"));
                        }
                        real[p] = rhs.re;
                        imaginary[p] = rhs.im;
                    }
                    pressure_factor.solve_in_place(&mut real);
                    checkpoint(cx)?;
                    pressure_factor.solve_in_place(&mut imaginary);
                    checkpoint(cx)?;
                    for p in 0..np {
                        pressure_pa[p] = C64::new(real[p], imaginary[p]);
                    }
                }
                Reduction::Exact {
                    pressure_indices, ..
                } => {
                    let recovery = self
                        .pressure_factor
                        .as_ref()
                        .ok_or_else(|| invalid("missing original pressure factor"))?;
                    let m = pressure_indices.len();
                    let mut real = vec![0.0; n];
                    let mut imaginary = vec![0.0; n];
                    for (i, &p) in pressure_indices.iter().enumerate() {
                        checkpoint(cx)?;
                        let mut rhs = C64::ZERO;
                        let (columns, values) = tangent.row(nu + p);
                        for (&dof, &value) in columns.iter().zip(values) {
                            if dof < nu {
                                rhs = rhs - displacement_m[dof].scale(value);
                            }
                        }
                        real[i] = rhs.re / recovery.scale;
                        imaginary[i] = rhs.im / recovery.scale;
                        if !real[i].is_finite() || !imaginary[i].is_finite() {
                            return Err(invalid("support constraint right hand side overflow"));
                        }
                    }
                    // B_free^T/scale = Q1 R. Thus the minimum Euclidean
                    // norm particular displacement is Q1 R^-T rhs/scale.
                    recovery.upper.solve_transpose(&mut real[..m]);
                    checkpoint(cx)?;
                    recovery.upper.solve_transpose(&mut imaginary[..m]);
                    recovery.qr.apply_q(&mut real);
                    checkpoint(cx)?;
                    recovery.qr.apply_q(&mut imaginary);
                    checkpoint(cx)?;
                    for (i, &dof) in free.iter().enumerate() {
                        displacement_m[dof] = C64::new(real[i], imaginary[i]);
                    }
                }
            }
        }
        let omega2 = self.omega * self.omega;
        let mut projected_force = vec![C64::ZERO; r];
        for (i, &dof) in free.iter().enumerate() {
            checkpoint(cx)?;
            let mut rhs = force_n[dof];
            if !support_motion_m.is_empty() {
                let (columns, values) = tangent.row(dof);
                for (&col, &value) in columns.iter().zip(values) {
                    let motion = if col < nu {
                        displacement_m[col]
                    } else {
                        pressure_pa[col - nu]
                    };
                    rhs = rhs - motion.scale(value);
                }
                let (columns, values) = mass.row(dof);
                for (&col, &value) in columns.iter().zip(values) {
                    rhs = rhs + displacement_m[col].scale(omega2 * value);
                }
                if !finite(rhs) {
                    return Err(invalid("lifted harmonic right hand side overflow"));
                }
            }
            match self.reduction {
                Reduction::Finite { .. } => projected_force[i] = rhs,
                Reduction::Exact { basis, .. } => {
                    for j in 0..r {
                        projected_force[j] = projected_force[j] + rhs.scale(basis[i * r + j]);
                    }
                }
            }
        }
        let projected_response = self
            .factor
            .as_ref()
            .map(|factor| factor.solve(&projected_force, cx))
            .transpose()?;
        if let Some(response) = &projected_response {
            for (i, &dof) in free.iter().enumerate() {
                checkpoint(cx)?;
                let correction = match self.reduction {
                    Reduction::Finite { .. } => response.displacement[i],
                    Reduction::Exact { basis, .. } => (0..r).fold(C64::ZERO, |sum, j| {
                        sum + response.displacement[j].scale(basis[i * r + j])
                    }),
                };
                displacement_m[dof] = displacement_m[dof] + correction;
            }
        }
        match self.reduction {
            Reduction::Finite {
                pressure_response, ..
            } => {
                for p in 0..np {
                    checkpoint(cx)?;
                    pressure_pa[p] =
                        free.iter()
                            .enumerate()
                            .fold(pressure_pa[p], |sum, (j, &dof)| {
                                sum + displacement_m[dof].scale(pressure_response[p * n + j])
                            });
                }
            }
            Reduction::Exact {
                pressure_indices, ..
            } => {
                let recovery = self
                    .pressure_factor
                    .as_ref()
                    .ok_or_else(|| invalid("missing original pressure factor"))?;
                let mut real = Vec::with_capacity(n);
                let mut imaginary = Vec::with_capacity(n);
                for &dof in free {
                    checkpoint(cx)?;
                    let mut rhs = force_n[dof];
                    let (columns, values) = tangent.row(dof);
                    for (&col, &value) in columns.iter().zip(values) {
                        if col < nu {
                            rhs = rhs - displacement_m[col].scale(value);
                        }
                    }
                    let (columns, values) = mass.row(dof);
                    for (&col, &value) in columns.iter().zip(values) {
                        rhs = rhs + displacement_m[col].scale(omega2 * value);
                    }
                    if !finite(rhs) {
                        return Err(invalid("pressure momentum right hand side overflow"));
                    }
                    real.push(rhs.re);
                    imaginary.push(rhs.im);
                }
                let real = recovery.qr.solve_ls(&real);
                let imaginary = recovery.qr.solve_ls(&imaginary);
                checkpoint(cx)?;
                for ((&p, real), imaginary) in pressure_indices.iter().zip(real).zip(imaginary) {
                    pressure_pa[p] = C64::new(real / recovery.scale, imaginary / recovery.scale);
                }
            }
        }
        for values in [&displacement_m, &pressure_pa] {
            for chunk in values.chunks(128) {
                checkpoint(cx)?;
                if chunk.iter().any(|&v| !finite(v)) {
                    return Err(invalid("physical mixed harmonic response overflow"));
                }
            }
        }
        let (mut momentum_residual_norm, mut momentum_scale) = (0.0_f64, 0.0_f64);
        let (mut pressure_residual_norm, mut pressure_scale) = (0.0_f64, 0.0_f64);
        let mut support_reaction_n = Vec::with_capacity(nu - n);
        for row in 0..nu + np {
            checkpoint(cx)?;
            let mut residual = if row < nu { -force_n[row] } else { C64::ZERO };
            let mut magnitude = if row < nu { force_n[row].abs() } else { 0.0 };
            let (columns, values) = tangent.row(row);
            for (&col, &coefficient) in columns.iter().zip(values) {
                let value = if col < nu {
                    displacement_m[col]
                } else {
                    pressure_pa[col - nu]
                };
                let term = value.scale(coefficient);
                residual = residual + term;
                magnitude += term.abs();
            }
            if row < nu {
                let (columns, values) = mass.row(row);
                for (&col, &coefficient) in columns.iter().zip(values) {
                    let term = displacement_m[col].scale(omega2 * coefficient);
                    residual = residual - term;
                    magnitude += term.abs();
                }
            }
            if !finite(residual) || !magnitude.is_finite() {
                return Err(invalid("original mixed equation or reaction overflow"));
            }
            if row >= nu {
                pressure_residual_norm = pressure_residual_norm.max(residual.abs());
                pressure_scale = pressure_scale.max(magnitude);
            } else if free.binary_search(&row).is_ok() {
                momentum_residual_norm = momentum_residual_norm.max(residual.abs());
                momentum_scale = momentum_scale.max(magnitude);
            } else {
                support_reaction_n.push((row, residual));
            }
        }
        // Match the original mixed vibration owner: separately normalized
        // physical equation blocks, not ratios on nearly zero individual rows.
        // A pressure-only traction can have exactly zero force/displacement on
        // some rows; their floating-point noise is measured against the same
        // original force scale as the other free rows, never an absolute floor.
        let momentum_relative_residual = if momentum_scale == 0.0 {
            0.0
        } else {
            momentum_residual_norm / momentum_scale
        };
        let pressure_relative_residual = if pressure_scale == 0.0 {
            0.0
        } else {
            pressure_residual_norm / pressure_scale
        };
        for (equation, relative, limit) in [
            (
                "free momentum",
                momentum_relative_residual,
                self.momentum_tolerance,
            ),
            (
                "all original pressure equations",
                pressure_relative_residual,
                pressure_tolerance,
            ),
        ] {
            if !relative.is_finite() || relative > limit {
                return Err(MixedHarmonicError::Residual {
                    equation,
                    relative,
                    limit,
                });
            }
        }
        let mut velocity_m_s = Vec::with_capacity(nu);
        let mut average_input_power_w = 0.0;
        for (i, (&u, &force)) in displacement_m.iter().zip(force_n).enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            let velocity = C64::new(self.omega * u.im, -self.omega * u.re);
            average_input_power_w += 0.5 * (force.conj() * velocity).re;
            if !finite(velocity) || !average_input_power_w.is_finite() {
                return Err(invalid("mixed harmonic velocity or power overflow"));
            }
            velocity_m_s.push(velocity);
        }
        let mut average_support_power_w = 0.0;
        for &(dof, reaction) in &support_reaction_n {
            checkpoint(cx)?;
            average_support_power_w += 0.5 * (reaction.conj() * velocity_m_s[dof]).re;
        }
        let power_balance_residual_w = average_input_power_w + average_support_power_w;
        if !average_support_power_w.is_finite() || !power_balance_residual_w.is_finite() {
            return Err(invalid("mixed harmonic support power overflow"));
        }
        checkpoint(cx)?;
        Ok(MixedHarmonicResponse {
            displacement_m,
            velocity_m_s,
            pressure_pa,
            support_reaction_n,
            momentum_relative_residual,
            pressure_relative_residual,
            average_input_power_w,
            average_support_power_w,
            power_balance_residual_w,
            projected_response,
        })
    }
}
