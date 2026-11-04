//! Resolved pressure waves in stationary, spatially varying fluids.
//!
//! The physical pressure equation is
//! `div((1/rho) grad p) + omega^2/(rho c^2) p = 0`, with `exp(+i omega t)`.
//! Density belongs INSIDE the divergence: varying only the local wavenumber
//! loses the correct normal-velocity continuity and impedance reflection.
//! The inhomogeneous pressure formulation is also documented in the COMSOL
//! Acoustics Module User's Guide, "Pressure Acoustics, Frequency Domain
//! Equations": <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.164.html>.
//!
//! This operator borrows the original [`HelmholtzFe`] solely for geometry,
//! hierarchical pressure basis and boundary quadrature. Its homogeneous
//! medium and geometric matrices remain unchanged. One coherent density/sound
//! speed pair is evaluated at each volume quadrature point, then both physical
//! matrices are assembled together. The material caller owns source identity,
//! equilibrium, and primitive-field interpolation; no derived property is
//! interpolated here. Jumps must conform to tetrahedron faces. Smooth fields
//! require independent mesh, polynomial-order and quadrature convergence.
//!
//! [`HeterogeneousHelmholtzFe::modes`] uses the original sparse modal owner.
//! [`HeterogeneousHelmholtzFe::prepare`] instead factors the FULL pressure space
//! with the original complex LU owner; no modal reduction is involved. With
//! inward imposed velocity `v` and outward specific admittance `Y`, the system
//! is `(K_rho - omega^2 M_beta + i omega C_Y) p = i omega b_v`. No reference
//! density, speed, or port-area normalization is introduced.
//!
//! Scope: linear, isotropic, lossless bulk fluids at a prescribed frozen
//! equilibrium without mean flow; straight-sided tetrahedra. Passive local
//! boundary admittances may dissipate energy. This does not evolve the source
//! temperature/composition, add thermoviscous or molecular bulk losses, or
//! turn a homogeneous wall model into a heterogeneous one.

use super::{AcousticFeError, HelmholtzFe, MAX_MODES, Medium};
use crate::highorder::simplex::duffy_quadrature;
use fs_la::eigen_complex::LuComplexError;
use fs_math::{c64::C64, det};
use fs_modal::{SliceOptions, slice_window};
use fs_sparse::{Coo, Csr};

mod harmonic;
pub use harmonic::{
    BoundaryAdmittance, HarmonicLimits, HeterogeneousPressureField, MAX_HARMONIC_DOFS,
    PreparedHeterogeneous, VelocityTrace,
};

/// Broadband evolution in the original pressure-potential coordinates.
pub mod time;

/// Largest admitted number of one-dimensional Gauss points per Duffy axis.
pub const MAX_QUADRATURE_ORDER: usize = 16;

/// Refusals of modal, harmonic and transient pressure solves.
#[derive(Debug, Clone, PartialEq)]
pub enum HeterogeneousError {
    /// The original geometry, trace or modal owner refused.
    Acoustic(AcousticFeError),
    /// The original complex LU owner refused.
    Factor(LuComplexError),
    /// The original passive mechanical wall owner refused its physical data.
    Wall(fs_phs::PhsError),
    /// The caller cancelled; no partial field/factor is returned.
    Cancelled,
    /// A declared preparation limit would be exceeded.
    Limit {
        /// Which resource.
        what: &'static str,
        /// Required amount.
        requested: usize,
        /// Admitted amount.
        maximum: usize,
    },
    /// The original, unscaled equations did not meet the requested tolerance.
    Residual {
        /// Componentwise physical backward error.
        observed: f64,
        /// Requested tolerance.
        limit: f64,
    },
    /// A transient candidate did not balance physical storage, source work
    /// and independently integrated boundary dissipation.
    EnergyBalance {
        /// Signed energy-change-plus-loss-minus-input residual [J].
        observed_j: f64,
        /// Admitted absolute defect, including evaluated roundoff [J].
        limit_j: f64,
    },
}

impl core::fmt::Display for HeterogeneousError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "heterogeneous pressure acoustics refused: {self:?}")
    }
}

impl std::error::Error for HeterogeneousError {}

impl From<AcousticFeError> for HeterogeneousError {
    fn from(value: AcousticFeError) -> Self {
        Self::Acoustic(value)
    }
}

impl From<LuComplexError> for HeterogeneousError {
    fn from(value: LuComplexError) -> Self {
        match value {
            LuComplexError::Cancelled => Self::Cancelled,
            _ => Self::Factor(value),
        }
    }
}

fn invalid(what: &'static str) -> AcousticFeError {
    AcousticFeError::InvalidInput { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HeterogeneousError> {
    if cancelled() {
        Err(HeterogeneousError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

/// Physical eigenpairs of the heterogeneous rigid-wall cavity.
///
/// These are deliberately distinct from the homogeneous owner's `k^2` modes.
#[derive(Debug, Clone)]
pub struct HeterogeneousModes {
    /// Eigenvalues `omega^2` [rad^2/s^2], ascending.
    pub omega_squared_rad2_s2: Vec<f64>,
    /// Eigenfrequencies `sqrt(max(omega^2, 0))/(2 pi)` [Hz].
    pub frequency_hz: Vec<f64>,
    /// Shapes normalized by the physical compressibility mass, `phi^T M_beta phi = 1`.
    pub shapes: Vec<Vec<f64>>,
    /// Original modal-owner residuals in the physical mass norm.
    pub residual: Vec<f64>,
    /// Inertia-certified eigenvalue count strictly below the low endpoint.
    pub below_low: usize,
    /// Inertia-certified eigenvalue count strictly below the high endpoint.
    pub below_high: usize,
}

/// Coefficient-weighted pressure operator on an existing, unchanged FE space.
pub struct HeterogeneousHelmholtzFe<'f, 'c> {
    geometry: &'f HelmholtzFe<'c>,
    stiffness: Csr,
    mass: Csr,
    quadrature_order: usize,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Evaluate original material properties once at every physical quadrature
    /// point and assemble `K_rho` and `M_beta` in one pass.
    ///
    /// The callback receives the original tetrahedron index and SI position.
    /// Its error is propagated unchanged. `quadrature_order` must lie between
    /// `geometry.order()+2` and [`MAX_QUADRATURE_ORDER`]; the lower bound avoids
    /// underintegrating the existing pressure basis even for constant media.
    /// Geometry's homogeneous medium is not sampled or used for these matrices.
    /// Cancellation is polled at element and quadrature boundaries, and is an
    /// `AcousticFeError::InvalidInput` in this source-generic constructor.
    pub fn new<E: From<AcousticFeError>>(
        geometry: &'f HelmholtzFe<'c>,
        quadrature_order: usize,
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<Medium, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, E> {
        let checkpoint = |cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() {
                Err(invalid("heterogeneous pressure assembly cancelled"))
            } else {
                Ok(())
            }
        };
        checkpoint(cancelled)?;
        if quadrature_order < geometry.order + 2 || quadrature_order > MAX_QUADRATURE_ORDER {
            return Err(invalid("heterogeneous quadrature order outside admitted range").into());
        }
        if geometry.ndof() == 0 || geometry.tetrahedra().is_empty() {
            return Err(invalid("heterogeneous pressure requires a nonempty volume").into());
        }
        let metrics =
            crate::whitney::try_element_geometry(geometry.space.complex, &geometry.positions)
                .map_err(|_| invalid("unrepresentable heterogeneous volume geometry"))?;
        let quad = duffy_quadrature(quadrature_order);
        let mut stiffness = Coo::new(geometry.ndof(), geometry.ndof());
        let mut mass = Coo::new(geometry.ndof(), geometry.ndof());
        for (element, vertices) in geometry.tetrahedra().iter().enumerate() {
            checkpoint(cancelled)?;
            let dofs = geometry.space.element_dofs(element);
            let count = dofs.len();
            let jacobian = 6.0 * metrics.vol_signed[element].abs();
            let mut local_k = vec![0.0; count * count];
            let mut local_m = vec![0.0; count * count];
            let mut values = vec![0.0; count];
            let mut gradients = vec![[0.0; 3]; count];
            for &(barycentric, weight) in &quad {
                checkpoint(cancelled)?;
                let mut point = [0.0; 3];
                for (&vertex, &bary) in vertices.iter().zip(&barycentric) {
                    for (axis, coordinate) in point.iter_mut().enumerate() {
                        *coordinate += bary * geometry.positions[vertex as usize][axis];
                    }
                }
                let medium = material_at(element, point)?;
                if !medium.density_kg_m3.is_finite() || medium.density_kg_m3 <= 0.0 {
                    return Err(AcousticFeError::InvalidMedium { what: "density" }.into());
                }
                if !medium.sound_speed_m_s.is_finite() || medium.sound_speed_m_s <= 0.0 {
                    return Err(AcousticFeError::InvalidMedium {
                        what: "sound speed",
                    }
                    .into());
                }
                let inverse_density = 1.0 / medium.density_kg_m3;
                let compressibility =
                    inverse_density / medium.sound_speed_m_s / medium.sound_speed_m_s;
                let wk = weight * jacobian * inverse_density;
                let wm = weight * jacobian * compressibility;
                if !wk.is_finite() || wk <= 0.0 || !wm.is_finite() || wm <= 0.0 {
                    return Err(
                        invalid("unrepresentable heterogeneous acoustic coefficients").into(),
                    );
                }
                for (i, (_, shape)) in dofs.iter().enumerate() {
                    values[i] = shape.eval(barycentric, geometry.order);
                    gradients[i] = [0.0; 3];
                    for (a, derivative) in shape
                        .d_lambda(barycentric, geometry.order)
                        .iter()
                        .enumerate()
                    {
                        for axis in 0..3 {
                            gradients[i][axis] += derivative * metrics.grads[element][a][axis];
                        }
                    }
                }
                // Assemble one triangle and mirror it, retaining exact discrete
                // reciprocity even after floating-point quadrature accumulation.
                for i in 0..count {
                    for j in i..count {
                        let gradient_dot: f64 = gradients[i]
                            .iter()
                            .zip(&gradients[j])
                            .map(|(a, b)| a * b)
                            .sum();
                        local_k[i * count + j] += wk * gradient_dot;
                        local_m[i * count + j] += wm * values[i] * values[j];
                    }
                }
            }
            for (i, &(row, _)) in dofs.iter().enumerate() {
                for (j, &(col, _)) in dofs.iter().enumerate().skip(i) {
                    let k = local_k[i * count + j];
                    let m = local_m[i * count + j];
                    if !k.is_finite() || !m.is_finite() {
                        return Err(invalid("heterogeneous matrix assembly overflow").into());
                    }
                    stiffness.push(row, col, k);
                    mass.push(row, col, m);
                    if i != j {
                        stiffness.push(col, row, k);
                        mass.push(col, row, m);
                    }
                }
            }
        }
        checkpoint(cancelled)?;
        let stiffness = stiffness.assemble();
        let mass = mass.assemble();
        for row in 0..geometry.ndof() {
            checkpoint(cancelled)?;
            if stiffness
                .row(row)
                .1
                .iter()
                .chain(mass.row(row).1)
                .any(|v| !v.is_finite())
            {
                return Err(invalid("heterogeneous global matrix assembly overflow").into());
            }
        }
        Ok(Self {
            geometry,
            stiffness,
            mass,
            quadrature_order,
        })
    }

    /// Original geometry/space owner; its uniform medium is NOT this medium.
    #[must_use]
    pub const fn geometry(&self) -> &'f HelmholtzFe<'c> {
        self.geometry
    }

    /// Original pressure degrees of freedom.
    #[must_use]
    pub fn ndof(&self) -> usize {
        self.geometry.ndof()
    }

    /// Physical stiffness `integral((1/rho) grad phi_i . grad phi_j)` [m^4/kg].
    #[must_use]
    pub const fn stiffness(&self) -> &Csr {
        &self.stiffness
    }

    /// Physical compressibility mass `integral(phi_i phi_j/(rho c^2))` [m^3/Pa].
    #[must_use]
    pub const fn mass(&self) -> &Csr {
        &self.mass
    }

    /// Declared Gauss-point count on each Duffy axis.
    #[must_use]
    pub const fn quadrature_order(&self) -> usize {
        self.quadrature_order
    }

    /// Rigid-wall modes in the frequency window. A negative lower endpoint
    /// includes the constant mode. The original sparse modal owner supplies
    /// inertia counts and residuals in the physical matrices.
    ///
    /// Cancellation brackets `slice_window`; that existing single-window
    /// owner does not expose an internal cancellation callback.
    pub fn modes(
        &self,
        window_hz: (f64, f64),
        options: &SliceOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousModes, HeterogeneousError> {
        poll(cancelled)?;
        let to_lambda = |frequency: f64| {
            let omega = core::f64::consts::TAU * frequency;
            omega * omega.abs()
        };
        let report = slice_window(
            &self.stiffness,
            &self.mass,
            (to_lambda(window_hz.0), to_lambda(window_hz.1)),
            options,
        )
        .map_err(AcousticFeError::from)?;
        poll(cancelled)?;
        if report.modes.len() > MAX_MODES {
            return Err(AcousticFeError::TooManyModes {
                found: report.modes.len(),
                cap: MAX_MODES,
            }
            .into());
        }
        Ok(HeterogeneousModes {
            omega_squared_rad2_s2: report.modes.iter().map(|m| m.lambda).collect(),
            frequency_hz: report
                .modes
                .iter()
                .map(|m| det::sqrt(m.lambda.max(0.0)) / core::f64::consts::TAU)
                .collect(),
            residual: report.modes.iter().map(|m| m.residual).collect(),
            shapes: report.modes.into_iter().map(|m| m.phi).collect(),
            below_low: report.below_low,
            below_high: report.below_high,
        })
    }

    /// Evaluate original hierarchical pressure coefficients inside one original
    /// tetrahedron. `barycentric` follows that tetrahedron's vertex order.
    pub fn pressure_at(
        &self,
        pressure_pa: &[C64],
        element: usize,
        barycentric: [f64; 4],
    ) -> Result<C64, HeterogeneousError> {
        if pressure_pa.len() != self.ndof()
            || element >= self.geometry.tetrahedra().len()
            || barycentric
                .iter()
                .any(|b| !b.is_finite() || *b < 0.0 || *b > 1.0)
            || (barycentric.iter().sum::<f64>() - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(invalid(
                "pressure sample requires original coefficients and barycentric point",
            )
            .into());
        }
        let mut result = C64::ZERO;
        for (dof, shape) in self.geometry.space.element_dofs(element) {
            result = result + pressure_pa[dof].scale(shape.eval(barycentric, self.geometry.order));
        }
        if !finite(result) {
            return Err(invalid("nonfinite pressure sample").into());
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
