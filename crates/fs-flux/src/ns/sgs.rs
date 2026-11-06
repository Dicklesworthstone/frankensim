//! Implicit Smagorinsky stress on the native planar BDM1/P0 flow.
//!
//! This is an Estimated, dissipative planar screen. A two-dimensional velocity
//! field cannot resolve a three-dimensional turbulent cascade; this module
//! supplies neither DNS/LES validation, wall modelling nor transition prediction.
//!
//! The closure is nu_t = (C_s Delta)^2 gamma, gamma = sqrt(2 D:D), and augments
//! the original symmetric viscous stress by 2 nu_t D. See NIST TN 1944,
//! https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.1944.pdf, section 2.
//! The original nonlinear owner evaluates strain, convection and the accepted
//! residual at the same candidate. Filter widths and molecular viscosity are
//! explicitly authored, frozen cell fields; no mesh length or wall damping is
//! guessed. Pressure remains kinematic, with no separately evolved SGS energy.
//!
//! Cell-interior strain comes directly from the affine BDM basis. Tangential
//! jumps retain the existing harmonic SIP treatment. Reported transfer is only
//! the volumetric SGS contribution, not SIP/upwind numerical dissipation or a
//! whole-step energy balance. Multiplying it by an authored density is required
//! before interpreting it as heat. Deterministic traversal and bounded nonlinear
//! work inherit the native solver; no new cancellation context is introduced.

use super::boundary::{
    NavierStokesBoundary, NavierStokesControls, NavierStokesStep, StokesBoundaryError,
};
use super::{FluxParams, FluxSystem};
use std::cell::RefCell;

/// Admission, arithmetic or native nonlinear-solve refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarSgsError {
    /// The coefficient, field sizes or velocity vector is inadmissible.
    InvalidInput(&'static str),
    /// A cell's authored geometry or transport input is inadmissible.
    InvalidCell {
        cell: usize,
        reason: &'static str,
    },
    /// The local strain, eddy viscosity or effective viscosity is unrepresentable.
    UnrepresentableTransport { cell: usize },
    /// A local or integrated modeled transfer is unrepresentable.
    UnrepresentableDissipation { cell: usize },
    /// The original boundary, history, linear or nonlinear solve refused.
    Flow(StokesBoundaryError),
}

impl core::fmt::Display for PlanarSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar Smagorinsky: {self:?}")
    }
}

impl std::error::Error for PlanarSgsError {}

impl From<StokesBoundaryError> for PlanarSgsError {
    fn from(error: StokesBoundaryError) -> Self {
        Self::Flow(error)
    }
}

/// Cell-interior diagnostics at one explicitly supplied velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarSgsCell {
    /// Native symmetric strain invariant sqrt(2 D:D) [1/s].
    pub strain_rate_per_s: f64,
    /// Modeled kinematic eddy viscosity [m²/s].
    pub eddy_viscosity_m2_s: f64,
    /// Molecular plus modeled kinematic viscosity used by momentum [m²/s].
    pub effective_viscosity_m2_s: f64,
    /// Resolved-to-modeled transfer nu_t * gamma², per unit mass [m²/s³].
    /// This is an instantaneous nonnegative rate, not accumulated heat.
    pub dissipation_m2_s3: f64,
}

/// Native-order cell diagnostics, with the actual triangle-area integral.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarSgsField {
    /// One record per original triangle, with no resampling.
    pub cells: Vec<PlanarSgsCell>,
    /// Integral of nu_t * gamma² over area [m⁴/s³].
    /// Multiply by constant density [kg/m³] for power per span [W/m].
    /// This excludes molecular, SIP and upwind contributions.
    pub integrated_dissipation_m4_s3: f64,
}

/// An accepted implicit momentum step and its own endpoint SGS diagnostics.
pub struct PlanarSgsStep {
    /// Original accepted field, pressure datums, boundary policies and residuals.
    pub step: NavierStokesStep,
    /// Evaluated from the returned velocity, never from the previous iterate.
    pub sgs: PlanarSgsField,
}

/// Prepared, mesh-bound planar Smagorinsky stress with frozen molecular inputs.
///
/// The filter widths and viscosities are copied at construction. A zero
/// coefficient gives the exact molecular-only constitutive law. The closure is
/// nonnegative and isotropic; there is no dynamic fit or energy backscatter.
pub struct PlanarSmagorinsky<'s, 'm> {
    system: &'s FluxSystem<'m>,
    mixing_lengths_m: Vec<f64>,
    molecular_viscosity_m2_s: Vec<f64>,
}

impl<'s, 'm> PlanarSmagorinsky<'s, 'm> {
    /// Bind C_s, explicit filter widths [m] and molecular viscosity [m²/s] to
    /// the original triangle order. C_s must be finite and nonnegative; both
    /// cell fields must be finite and strictly positive. No values are clipped.
    ///
    /// # Errors
    /// Refuses malformed sizes, nonpositive inputs and unrepresentable C_s Delta.
    pub fn new(
        system: &'s FluxSystem<'m>,
        coefficient: f64,
        filter_widths_m: &[f64],
        molecular_viscosity_m2_s: &[f64],
    ) -> Result<Self, PlanarSgsError> {
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err(PlanarSgsError::InvalidInput(
                "finite nonnegative Smagorinsky coefficient required",
            ));
        }
        let count = system.mesh.tris.len();
        if count == 0
            || filter_widths_m.len() != count
            || molecular_viscosity_m2_s.len() != count
        {
            return Err(PlanarSgsError::InvalidInput(
                "one filter width and molecular viscosity per nonempty triangle mesh required",
            ));
        }
        let mut mixing_lengths_m = Vec::with_capacity(count);
        for cell in 0..count {
            let width = filter_widths_m[cell];
            let nu = molecular_viscosity_m2_s[cell];
            let length = coefficient * width;
            let reason = if !width.is_finite() || width <= 0.0 {
                Some("positive finite filter width required")
            } else if !nu.is_finite() || nu <= 0.0 {
                Some("positive finite molecular viscosity required")
            } else if !system.mesh.areas[cell].is_finite() || system.mesh.areas[cell] <= 0.0 {
                Some("positive finite triangle area required")
            } else if !length.is_finite() || (coefficient > 0.0 && length == 0.0) {
                Some("Smagorinsky mixing length is unrepresentable")
            } else {
                None
            };
            if let Some(reason) = reason {
                return Err(PlanarSgsError::InvalidCell { cell, reason });
            }
            mixing_lengths_m.push(length);
        }
        Ok(Self {
            system,
            mixing_lengths_m,
            molecular_viscosity_m2_s: molecular_viscosity_m2_s.to_vec(),
        })
    }

    // Constant-work, allocation-free cell response. This multiplication order
    // avoids squaring a length before its compensating strain scale is applied.
    fn viscosities(&self, cell: usize, rate: f64) -> Result<(f64, f64), PlanarSgsError> {
        let length = self.mixing_lengths_m[cell];
        let eddy = length * (length * rate);
        let effective = self.molecular_viscosity_m2_s[cell] + eddy;
        if !rate.is_finite()
            || rate < 0.0
            || !eddy.is_finite()
            || !effective.is_finite()
            || (length > 0.0 && rate > 0.0 && eddy == 0.0)
        {
            return Err(PlanarSgsError::UnrepresentableTransport { cell });
        }
        Ok((eddy, effective))
    }

    /// Evaluate the modeled transport and transfer on original velocity DOFs.
    /// Supply exactly n_u coefficients; pressure entries are not velocity.
    ///
    /// # Errors
    /// Refuses invalid velocity or unrepresentable local/integrated diagnostics.
    pub fn field(&self, velocity: &[f64]) -> Result<PlanarSgsField, PlanarSgsError> {
        if velocity.len() != self.system.n_u || velocity.iter().any(|value| !value.is_finite()) {
            return Err(PlanarSgsError::InvalidInput(
                "finite original velocity coefficients required",
            ));
        }
        let rates = self.system.shear_rates(velocity);
        let mut cells = Vec::with_capacity(rates.len());
        let mut integrated = 0.0;
        for (cell, rate) in rates.into_iter().enumerate() {
            let (eddy, effective) = self.viscosities(cell, rate)?;
            let dissipation = (eddy * rate) * rate;
            let contribution = self.system.mesh.areas[cell] * dissipation;
            integrated += contribution;
            if !dissipation.is_finite()
                || !contribution.is_finite()
                || !integrated.is_finite()
                || (eddy > 0.0 && rate > 0.0 && dissipation == 0.0)
                || (dissipation > 0.0 && contribution == 0.0)
            {
                return Err(PlanarSgsError::UnrepresentableDissipation { cell });
            }
            cells.push(PlanarSgsCell {
                strain_rate_per_s: rate,
                eddy_viscosity_m2_s: eddy,
                effective_viscosity_m2_s: effective,
                dissipation_m2_s3: dissipation,
            });
        }
        Ok(PlanarSgsField {
            cells,
            integrated_dissipation_m4_s3: integrated,
        })
    }

    /// Advance BDF1 with implicit SGS stress and convection at the same endpoint.
    ///
    /// Uses the original mixed velocity/total-traction/backflow policies, mass
    /// matrix, harmonic SIP stress and independently reassembled residual.
    /// params.nu is unused: the prepared cell molecular field is authoritative.
    /// Force and boundary callbacks must be pure and describe the new time level.
    ///
    /// # Errors
    /// Preserves native solve refusals and returns no partial flow on SGS failure.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<PlanarSgsStep, PlanarSgsError> {
        self.bdf2_step(
            params,
            force,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
        )
    }

    /// Advance equal-step BDF2, or exactly BDF1 when older_velocity is None.
    ///
    /// With two histories the derivative is (3 u_new - 4 u_n + u_old)/(2 dt);
    /// both history intervals must equal timestep_s. The original owner uses
    /// actual u_n to initialize the closure, never the artificial mass-history
    /// combination. Endpoint diagnostics are instantaneous rates, not BDF2's
    /// two-history energy or a claim of whole-scheme energy monotonicity.
    ///
    /// # Errors
    /// Refuses invalid history/control/boundary data, unrepresentable transport
    /// and unconverged native solves. Neither supplied history is modified.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<PlanarSgsStep, PlanarSgsError> {
        let source_error = RefCell::new(None);
        let result = self.system.bdf2_navier_stokes_mixed_rheology_step_checked(
            params,
            &|cell, _, rate| match self.viscosities(cell, rate) {
                Ok((_, effective)) => effective,
                Err(error) => {
                    if source_error.borrow().is_none() {
                        *source_error.borrow_mut() = Some(error);
                    }
                    // A finite fallback permits the native assembly to drain;
                    // its source check refuses before a solve or publication.
                    self.molecular_viscosity_m2_s[cell]
                }
            },
            force,
            boundaries,
            previous_velocity,
            older_velocity,
            timestep_s,
            controls,
            &|| source_error.borrow().is_none(),
        );
        if let Some(error) = source_error.into_inner() {
            return Err(error);
        }
        let step = result?;
        let sgs = self.field(&step.solution.field.x[..self.system.n_u])?;
        Ok(PlanarSgsStep { step, sgs })
    }
}
