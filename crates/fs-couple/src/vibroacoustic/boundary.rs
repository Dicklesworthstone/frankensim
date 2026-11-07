//! Local pressure/normal-velocity boundary laws in the existing modal solve.
//!
//! With exp(-i omega t), outward v = Y p contributes -i omega B to the
//! unnormalized pressure weak form, B_qs = integral(psi_q Y psi_s dA).
//! Reference: the reactive boundary and acoustic weak form in
//! <https://docs.software.vt.edu/abaqusv2026/English/SIMACAETHERefMap/simathe-c-acouststruct.htm>.
//! That source uses the opposite time convention; its imaginary terms are
//! conjugated here. No new propagation or constitutive model is introduced.

use super::{C64, VibroError, VibroacousticModel};
use fs_exec::CancelGate;

/// A locally reacting cavity boundary evaluated at ONE angular frequency.
///
/// The caller supplies specific admittance Y [m/(Pa s)] in exp(-i omega t),
/// positive-area quadrature and the original pressure eigenfunctions sampled
/// there. These samples may cover several disjoint patches with different
/// materials. Signed shapes and all off-diagonal modal terms are retained.
/// Zero admittance is rigid; Re(Y) >= 0 is passive at this frequency.
///
/// This is a projection of an existing boundary law, not a material fit or a
/// broadband passivity/causality certificate. In particular a normalized
/// admittance must first be divided by its original rho*c; exp(+i omega t)
/// data must be conjugated. The producer owns material validity, spatial
/// resolution and local-reaction applicability. A nonlocal exterior needs
/// its full radiation operator, not independent local samples.
#[derive(Debug, Clone)]
pub struct CavityBoundaryAdmittance {
    pub(super) omega: f64,
    pub(super) matrix: Vec<C64>,
    modes: usize,
    // sqrt(A Re Y) psi_q, one row per quadrature point. Recovering work from
    // these physical channels preserves nonnegative loss at modal nullspaces.
    conductance_channels: Vec<Vec<f64>>,
}

fn invalid(what: &'static str) -> VibroError {
    VibroError::BadParameter { what }
}

fn poll(gate: &CancelGate) -> Result<(), VibroError> {
    if gate.is_requested() {
        Err(VibroError::Cancelled)
    } else {
        Ok(())
    }
}

impl CavityBoundaryAdmittance {
    /// Project actual surface samples onto the ORIGINAL cavity pressure basis.
    /// `pressure_modes[q][i]` is psi_q at boundary point i, before division by
    /// its volume norm. Areas [m^2] and specific admittances share that order.
    /// No absolute value, diagonal approximation or positive loss floor is used.
    ///
    /// # Errors
    /// Empty/mismatched inputs, nonfinite samples, negative areas/conductance,
    /// unrepresentable coefficients or a cancellation request. Nothing is
    /// attached to a model until every sample has been admitted and projected.
    pub fn from_samples(
        omega: f64,
        pressure_modes: &[Vec<f64>],
        areas_m2: &[f64],
        specific_admittance_m_pa_s: &[C64],
        gate: &CancelGate,
    ) -> Result<Self, VibroError> {
        poll(gate)?;
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid(
                "cavity boundary angular frequency must be finite and positive",
            ));
        }
        let modes = pressure_modes.len();
        let points = areas_m2.len();
        if modes == 0
            || points == 0
            || specific_admittance_m_pa_s.len() != points
            || pressure_modes.iter().any(|shape| shape.len() != points)
        {
            return Err(VibroError::ShapeMismatch {
                what: "boundary requires nonempty original pressure modes, areas and admittances on one quadrature",
            });
        }
        if areas_m2.iter().any(|a| !a.is_finite() || *a < 0.0)
            || pressure_modes.iter().flatten().any(|p| !p.is_finite())
            || specific_admittance_m_pa_s
                .iter()
                .any(|y| !y.re.is_finite() || !y.im.is_finite() || y.re < 0.0)
        {
            return Err(invalid(
                "boundary samples must be finite with nonnegative areas and conductance",
            ));
        }
        let cells = modes
            .checked_mul(modes)
            .ok_or_else(|| invalid("cavity boundary matrix size is unrepresentable"))?;
        let mut matrix = Vec::new();
        matrix
            .try_reserve_exact(cells)
            .map_err(|_| invalid("cavity boundary matrix allocation failed"))?;
        matrix.resize(cells, C64::ZERO);
        let mut conductance_channels = Vec::new();
        conductance_channels
            .try_reserve_exact(points)
            .map_err(|_| invalid("cavity boundary channel allocation failed"))?;
        for (i, (&area, &y)) in areas_m2.iter().zip(specific_admittance_m_pa_s).enumerate() {
            poll(gate)?;
            let root_area = area.sqrt();
            let root_y = y.re.sqrt();
            let channel: Vec<f64> = pressure_modes
                .iter()
                .map(|shape| (shape[i] * root_area) * root_y)
                .collect();
            if channel.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "cavity boundary conductance projection is unrepresentable",
                ));
            }
            for q in 0..modes {
                poll(gate)?;
                for s in q..modes {
                    let term = C64::new(
                        channel[q] * channel[s],
                        (pressure_modes[q][i] * root_area)
                            * (pressure_modes[s][i] * root_area)
                            * y.im,
                    );
                    let value = matrix[q * modes + s] + term;
                    if !value.re.is_finite() || !value.im.is_finite() {
                        return Err(invalid(
                            "cavity boundary modal projection is unrepresentable",
                        ));
                    }
                    matrix[q * modes + s] = value;
                    matrix[s * modes + q] = value;
                }
            }
            conductance_channels.push(channel);
        }
        poll(gate)?;
        Ok(Self {
            omega,
            matrix,
            modes,
            conductance_channels,
        })
    }

    /// Original angular frequency [rad/s]; reusing this sample at another
    /// frequency refuses instead of freezing a dispersive material silently.
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }

    /// Full original-basis B_qs = integral(psi_q Y psi_s dA), row-major.
    /// This is volume-flow/pressure admittance, not specific admittance.
    pub fn modal_admittance(&self) -> &[C64] {
        &self.matrix
    }

    pub(super) fn outward_power(&self, pressure: &[C64]) -> f64 {
        self.conductance_channels
            .iter()
            .map(|channel| {
                let p = pressure
                    .iter()
                    .zip(channel)
                    .fold(C64::ZERO, |p, (&a, &weight)| p + a.scale(weight));
                0.5 * p.norm_sq()
            })
            .sum()
    }
}

impl VibroacousticModel {
    /// Attach the original locally reacting boundary to every subsequent FRF.
    ///
    /// The sample must use this exact cavity basis and normalization. Its
    /// pressure reaction is solved simultaneously with structure, molecular
    /// bulk relaxation and exterior structural radiation. A truncated solve
    /// uses the same principal boundary block and original physical channels.
    /// Boundary power is included once in `PowerBreakdown::cavity` and is also
    /// reported separately as `PowerBreakdown::cavity_boundary`.
    ///
    /// Supply all disjoint patches in one sample. This rejects a second
    /// boundary and unspecified scalar cavity loss, avoiding silent replacement
    /// or double counting. Molecular relaxation remains an independent bulk
    /// mechanism. A boundary must not also be a prescribed copy of the same
    /// wall motion/source. This remains a modal approximation: basis and
    /// quadrature convergence are required independently of power closure.
    ///
    /// # Errors
    /// A different mode count, active scalar cavity loss, or an existing
    /// boundary. Subsequent solves also reject a different angular frequency.
    pub fn with_cavity_boundary_admittance(
        mut self,
        boundary: CavityBoundaryAdmittance,
    ) -> Result<Self, VibroError> {
        if boundary.modes != self.omegas_a.len() {
            return Err(VibroError::ShapeMismatch {
                what: "cavity boundary must use the original pressure mode count and basis",
            });
        }
        if self.eta_a != 0.0 || self.cavity_boundary.is_some() {
            return Err(invalid(
                "combine physical boundary patches once; scalar cavity loss is exclusive",
            ));
        }
        self.cavity_boundary = Some(boundary);
        Ok(self)
    }
}

#[cfg(test)]
mod tests;
