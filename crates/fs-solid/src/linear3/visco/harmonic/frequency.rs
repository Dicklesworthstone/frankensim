//! Actual material dispersion through the original cached structural factor.

use super::{
    ViscoElementError, ViscoElements, ViscoHarmonicError, ViscoHarmonicOperator,
    ViscoHarmonicResponse, checkpoint, dense, material_response,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicParameterControls, HarmonicPortParameterDirection, HarmonicPortParameterResponse,
    HarmonicResponseError,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;
use fs_sparse::Csr;

/// Original material-admitted response and its derivative per rad/s.
#[derive(Clone, Debug)]
pub struct ViscoHarmonicFrequencyResponse {
    /// Original motion and separate material/boundary work, admitted over the
    /// complete strain cycle at every original material quadrature point.
    pub response: ViscoHarmonicResponse,
    /// du/domega [m s], in the same original free-DOF order.
    pub displacement_derivative: Vec<C64>,
    /// dv/domega [m], including -i*u as well as -i*omega*du/domega.
    pub velocity_derivative: Vec<C64>,
    /// Original port flow, effort and signed absorbed-work slopes per rad/s.
    /// These retain authored attachment and impedance frequency derivatives.
    pub port_load: Option<HarmonicPortParameterResponse>,
    /// d Re(f^H v)/2 per rad/s [W s], including an authored force slope.
    pub average_input_power_derivative: f64,
    /// Derivative of the original constitutive heat rate per rad/s [W s].
    /// Includes the complete Prony dispersion and Kelvin viscosity once.
    pub average_material_heat_derivative: f64,
    /// Input minus independently differentiated material heat and boundary work.
    pub power_balance_residual_derivative: f64,
    /// Maximum tangent-solve and published-coordinate equation backward error;
    /// not a forward-error or neighboring-frequency applicability bound.
    pub relative_derivative_residual: f64,
}

/// Retain one physical factor and its assembled material-frequency direction.
/// No spectral assembly or factorization is repeated for another excitation.
pub struct ViscoHarmonicFrequencyOperator<'a> {
    base: ViscoHarmonicOperator<'a>,
    mass_derivative: Vec<f64>,
    damping_derivative: Vec<f64>,
    stiffness_derivative: Vec<f64>,
}

impl<'a> ViscoHarmonicOperator<'a> {
    /// Add the source solid's actual frequency dispersion to this cached factor.
    ///
    /// At this same admitted frequency, the original element quadrature assembles
    /// dK'/domega and d(K''/omega)/domega from each tensor-Prony law. The latter
    /// is evaluated directly by the material owner, without cancelling two
    /// nearly equal loss/omega terms. Reference mass and supports are fixed.
    ///
    /// The additional retained dense source payload is exactly 24*n*n bytes:
    /// zero mass slope, effective-viscosity slope and storage-stiffness slope.
    /// This gate precedes assembly/allocation. Original sparse mesh/COO and
    /// constitutive temporaries retain their existing ownership and limits.
    /// Each later solve separately admits its shared harmonic-tangent scratch.
    pub fn with_frequency_derivative(
        self,
        maximum_derivative_payload_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicFrequencyOperator<'a>, ViscoHarmonicError> {
        checkpoint(cx)?;
        let n = self.source.equilibrium.free_dofs.len();
        let entries = n.checked_mul(n).ok_or(HarmonicResponseError::Invalid(
            "material frequency derivative extent overflow",
        ))?;
        let bytes = entries
            .checked_mul(24)
            .ok_or(HarmonicResponseError::Invalid(
                "material frequency derivative payload overflow",
            ))?;
        if bytes > maximum_derivative_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "material frequency derivative source bytes",
                requested: bytes,
                cap: maximum_derivative_payload_bytes,
            }
            .into());
        }
        let [storage, viscosity] = self
            .source
            .spectral_frequency_derivative(self.omega_rad_s, cx)?;
        let stiffness_derivative = dense(&storage, n, cx)?;
        let damping_derivative = dense(&viscosity, n, cx)?;
        let mass_derivative = vec![0.0; entries];
        checkpoint(cx)?;
        Ok(ViscoHarmonicFrequencyOperator {
            base: self,
            mass_derivative,
            damping_derivative,
            stiffness_derivative,
        })
    }
}

impl<'a> ViscoHarmonicFrequencyOperator<'a> {
    /// The same original factor, still usable for ordinary fixed-frequency solves.
    pub fn base(&self) -> &ViscoHarmonicOperator<'a> {
        &self.base
    }

    /// Original node-major DOF indices corresponding to every response row.
    pub fn free_dofs(&self) -> &[usize] {
        self.base.free_dofs()
    }

    /// The physical angular frequency at which this factor/direction is prepared.
    pub fn omega_rad_s(&self) -> f64 {
        self.base.omega_rad_s()
    }

    /// Solve and differentiate the complete dispersive structural response.
    ///
    /// `force_derivative_n_s` is the authored df/domega; None fixes the force
    /// phasor. `port_direction` supplies complete dB/domega and dZ/domega in
    /// the original conjugate units; None explicitly fixes B and Z while their
    /// velocity-dependent load still receives its kinematic frequency term.
    /// No impedance dispersion or changing geometry is guessed from one sample.
    ///
    /// The shared total-parameter owner receives dM=0, the assembled material
    /// dC_eff/domega and dK'/domega, and d(omega)/d(omega)=1. Motion, intrinsic
    /// heat, source work and boundary work retain their independent derivatives.
    /// Before publication, the original full-cycle strain admission is applied
    /// to the primal displacement; a derivative is not itself a strain cycle.
    /// A new frequency still requires preparation and fresh source admission.
    pub fn solve(
        &self,
        force_n: &[C64],
        force_derivative_n_s: Option<&[C64]>,
        port_direction: Option<HarmonicPortParameterDirection<'_>>,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicFrequencyResponse, ViscoHarmonicError> {
        let result = self.base.response.solve_with_total_parameter_derivative(
            force_n,
            QuadraticParameterDirection {
                mass: &self.mass_derivative,
                damping: &self.damping_derivative,
                stiffness: &self.stiffness_derivative,
            },
            force_derivative_n_s,
            port_direction,
            1.0,
            controls,
            cx,
        )?;
        self.base
            .source
            .validate_harmonic_displacement(&result.response.displacement, cx)?;
        checkpoint(cx)?;
        Ok(ViscoHarmonicFrequencyResponse {
            response: material_response(result.response),
            displacement_derivative: result.displacement_derivative,
            velocity_derivative: result.velocity_derivative,
            port_load: result.port_load,
            average_input_power_derivative: result.average_input_power_derivative,
            average_material_heat_derivative: result.average_viscous_power_derivative,
            power_balance_residual_derivative: result.power_balance_residual_derivative,
            relative_derivative_residual: result.relative_derivative_residual,
        })
    }
}

impl ViscoElements {
    fn spectral_frequency_derivative(
        &self,
        omega: f64,
        cx: &Cx<'_>,
    ) -> Result<[Csr; 2], ViscoElementError> {
        // Every constitutive band is checked before any element integration.
        let tensors = self
            .elements
            .iter()
            .map(|element| {
                checkpoint(cx)?;
                let derivative = element.law.complex_modulus_frequency_derivative(omega)?;
                Ok([
                    derivative.storage_pa_s,
                    derivative.effective_viscosity_pa_s2,
                ])
            })
            .collect::<Result<Vec<_>, ViscoElementError>>()?;
        self.assemble_spectral_tensors(tensors, false, cx)
    }
}
