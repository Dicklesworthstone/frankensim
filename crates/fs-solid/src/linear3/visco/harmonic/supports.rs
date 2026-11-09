//! Source-material harmonic response to prescribed support motion.
//!
//! The original support solver receives M, K'(omega) and K''(omega)/omega,
//! including all free/support cross terms. The same original material points
//! then admit the complete displacement cycle before a response is published.

use super::{ViscoElements, ViscoHarmonicError, checkpoint};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedSupportedHarmonicResponse,
};

/// Original-coordinate material response and separately evaluated actuator work.
#[derive(Clone, Debug)]
pub struct ViscoSupportedHarmonicResponse {
    /// Peak exp(-i omega t) displacement [m] in `represented_dofs()` order.
    /// Prescribed support coordinates are retained in this complete vector.
    pub displacement_m: Vec<C64>,
    /// Peak velocity [m/s], including the moving supports.
    pub velocity_m_s: Vec<C64>,
    /// Force BY each support ON the structure [N], in `support_dofs()` order.
    /// Applied forces on support coordinates are subtracted exactly once.
    pub support_reaction_n: Vec<C64>,
    /// Componentwise backward error in the original free equations.
    pub relative_residual: f64,
    /// Signed cycle-mean work rate of every applied nodal force [W].
    pub average_input_power_w: f64,
    /// Signed cycle-mean work rate supplied by the moving supports [W].
    pub average_support_power_w: f64,
    /// Original Kelvin and Maxwell heat on ALL represented coordinates [W].
    /// Cross-support constitutive loss is retained without modal refitting.
    pub average_material_heat_w: f64,
    /// Applied power plus support power minus independently computed heat [W].
    pub power_balance_residual_w: f64,
}

/// One material-admitted frequency factor, reusable across force/support cases.
/// Original supports are selected by their global labels, never reduced rows.
pub struct ViscoSupportedHarmonicOperator<'a> {
    source: &'a ViscoElements,
    free_dofs: Vec<usize>,
    support_dofs: Vec<usize>,
    response: PreparedSupportedHarmonicResponse,
}

impl ViscoElements {
    /// Prepare base/actuator excitation with the actual dispersive material law.
    ///
    /// Assemble the source with every moving support coordinate retained, then
    /// supply its strictly increasing global DOF labels here. Homogeneous fixed
    /// coordinates already removed from the source remain at rest; a request to
    /// move an absent coordinate refuses. Empty or all-represented support sets
    /// are allowed. All original material frequency bands must admit omega > 0.
    ///
    /// The existing native solver applies A_ff u_f = f_f - A_fc u_c with
    /// A = K'(omega) - omega^2 M - i K''(omega). The original consistent mass,
    /// every tensor-Prony branch and Kelvin viscosity contribute their complete
    /// cross terms. No equilibrium-stiffness substitution, penalty or additional
    /// damping is introduced. Only the free principal block is factored.
    ///
    /// Dense numerical payload is bounded before assembly by
    /// 344*n*n + 544*n bytes, covering source conversion, coordinate maps,
    /// the shared solver and one solve. Native sparse/constitutive assembly
    /// retains the source's existing limits. This is a bounded dense particular
    /// response at fixed material conditions, not startup or thermal feedback.
    pub fn prepare_supported_harmonic_response(
        &self,
        support_dofs: &[usize],
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<ViscoSupportedHarmonicOperator<'_>, ViscoHarmonicError> {
        checkpoint(cx)?;
        let invalid = |message| ViscoHarmonicError::from(HarmonicResponseError::Invalid(message));
        let represented = &self.equilibrium.free_dofs;
        let n = represented.len();
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || n == 0
            || represented.windows(2).any(|pair| pair[0] >= pair[1])
            || !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid(
                "positive harmonic frequency, ordered represented DOFs and residual budget in (0,1) required",
            ));
        }
        if support_dofs.len() > n || support_dofs.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(invalid(
                "material support labels must be strictly increasing and represented",
            ));
        }
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "supported material coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            }
            .into());
        }
        let source_bytes = n
            .checked_mul(n)
            .and_then(|square| square.checked_mul(24))
            .and_then(|bytes| n.checked_mul(32)?.checked_add(bytes))
            .ok_or_else(|| invalid("supported material source payload overflow"))?;
        let total = n
            .checked_mul(n)
            .and_then(|square| square.checked_mul(344))
            .and_then(|bytes| n.checked_mul(544)?.checked_add(bytes))
            .ok_or_else(|| invalid("supported material payload overflow"))?;
        if total > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "supported material numerical payload bytes",
                requested: total,
                cap: controls.maximum_payload_bytes,
            }
            .into());
        }
        let mut support_rows = Vec::with_capacity(support_dofs.len());
        for &label in support_dofs {
            checkpoint(cx)?;
            support_rows.push(represented.binary_search(&label).map_err(|_| {
                invalid("moving material support was already eliminated or is absent")
            })?);
        }
        let [mass, damping, storage] = self.harmonic_pencil(omega_rad_s, cx)?;
        let response = PreparedSupportedHarmonicResponse::try_new(
            &mass,
            &damping,
            &storage,
            n,
            &support_rows,
            omega_rad_s,
            HarmonicResponseControls {
                maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
                ..controls
            },
            cx,
        )?;
        let free_dofs = response
            .free_dofs()
            .iter()
            .map(|&row| represented[row])
            .collect();
        checkpoint(cx)?;
        Ok(ViscoSupportedHarmonicOperator {
            source: self,
            free_dofs,
            support_dofs: support_dofs.to_vec(),
            response,
        })
    }
}

impl<'a> ViscoSupportedHarmonicOperator<'a> {
    /// Exact original geometry, immutable material laws and reference inertia.
    pub const fn source(&self) -> &'a ViscoElements {
        self.source
    }

    /// Original global labels of ALL displacement/velocity and applied-force rows.
    pub fn represented_dofs(&self) -> &[usize] {
        &self.source.equilibrium.free_dofs
    }

    /// Original global labels of coordinates solved by force balance.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Original global labels of prescribed coordinates and their reactions.
    pub fn support_dofs(&self) -> &[usize] {
        &self.support_dofs
    }

    pub fn omega_rad_s(&self) -> f64 {
        self.response.angular_frequency_rad_s()
    }

    /// Solve full applied-force [N] and support-displacement [m] peak phasors.
    ///
    /// Every original constitutive point and every compatible B-bar strain
    /// independently admits the whole cycle, including the support motion.
    /// Successful displacement vectors can enter this same source's single-tone
    /// or coherent-spectrum material observers directly. Material heat is only
    /// observed here; no thermal state or hereditary history is changed.
    /// Failure or cancellation publishes no partial response and leaves the
    /// cached factor usable for another admitted case.
    pub fn solve(
        &self,
        applied_force_n: &[C64],
        support_displacement_m: &[C64],
        cx: &Cx<'_>,
    ) -> Result<ViscoSupportedHarmonicResponse, ViscoHarmonicError> {
        let response = self
            .response
            .solve(applied_force_n, support_displacement_m, cx)?;
        self.source
            .validate_harmonic_displacement(&response.displacement, cx)?;
        checkpoint(cx)?;
        Ok(ViscoSupportedHarmonicResponse {
            displacement_m: response.displacement,
            velocity_m_s: response.velocity,
            support_reaction_n: response.support_reaction,
            relative_residual: response.relative_residual,
            average_input_power_w: response.average_input_power_w,
            average_support_power_w: response.average_support_power_w,
            average_material_heat_w: response.average_viscous_power_w,
            power_balance_residual_w: response.power_balance_residual_w,
        })
    }
}
