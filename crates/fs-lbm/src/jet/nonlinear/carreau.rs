//! Original thermal Carreau-Yasuda viscosity in the native slender-jet step.
//!
//! In incompressible uniaxial extension, D=diag(-e/2,-e/2,e), e=du/dz.
//! Thus gamma=sqrt(2 D:D)=sqrt(3)|e|, tau_zz-tau_rr=3 mu(gamma)e,
//! and tau:D=3 mu(gamma)e². Axial force, dissipation and the viscous step
//! bound therefore consume one original material sample at each RK stage.
//!
//! This is the generalized-Newtonian extension of the original slender
//! equations, not a separately calibrated extensional viscosity or an elastic
//! normal stress. See the generalized-Newtonian discussion in Eggers &
//! Villermaux (2008), Physics of liquid jets, equation (320):
//! <https://people.maths.bris.ac.uk/~majge/rpp8_3_036601.pdf>.
//! The existing radius/slope/resolution limits remain essential near a neck.

use super::*;
use fs_material::fluid_transport::rheology::{ThermalCarreauYasuda, ThermalRheologyState};

/// Fixed physical conditions and the original generalized-Newtonian law.
/// The temperature shift and calibrated rate support are never refitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlenderJetParameters {
    pub density_kg_m3: f64,
    pub surface_tension_n_m: f64,
    pub temperature_k: f64,
    pub rheology: ThermalCarreauYasuda,
}

impl CarreauSlenderJetParameters {
    pub(super) fn is_newtonian(&self) -> bool {
        let [mu0, mu_inf, time, index, _] = self.rheology.curve().parameters();
        time == 0.0 || index == 1.0 || mu0 == mu_inf
    }

    pub(super) fn sample(
        &self,
        axial_rate_per_s: f64,
        cell: Option<usize>,
    ) -> Result<ThermalRheologyState, NonlinearJetError> {
        let rate = finite(3.0_f64.sqrt() * axial_rate_per_s.abs())?;
        self.rheology
            .evaluate(self.temperature_k, rate)
            .map_err(|source| NonlinearJetError::Rheology { cell, source })
    }
}

impl PeriodicSlenderJet {
    /// Construct the same periodic jet with source-resolved Carreau viscosity.
    ///
    /// Every stress and loss evaluates the original law at the local
    /// sqrt(3)*|du/dz|, in inverse seconds. The source must support zero rate
    /// and every actual RK stage at this fixed temperature. Unsupported
    /// algebraic Newton line-search candidates may backtrack within the work
    /// budget; source refusal of an actual stage rolls the entire advance back.
    /// Neither path clips a rate or substitutes another viscosity. Density and surface tension are
    /// authored constants; this owner adds no caloric evolution, gravity,
    /// viscoelastic stress, gas drag or continuation through pinch-off.
    ///
    /// The native conservative flux supports explicit SSPRK3 and the original
    /// third-order IMEX viscosity policy. Nonlinear implicit stages retain the
    /// original material stress and differential tangent, with bounded Newton
    /// work selected by [`Self::set_implicit_viscosity_controls`]. Scratch is
    /// prepared on policy selection and reused by successful steps.
    pub fn new_carreau(
        parameters: CarreauSlenderJetParameters,
        length_m: f64,
        radii_m: &[f64],
        velocities_m_s: &[f64],
        limits: SlenderJetLimits,
        gate: &CancelGate,
    ) -> Result<Self, NonlinearJetError> {
        poll(gate)?;
        if !parameters.density_kg_m3.is_finite()
            || parameters.density_kg_m3 <= 0.0
            || !parameters.surface_tension_n_m.is_finite()
            || parameters.surface_tension_n_m < 0.0
        {
            return Err(NonlinearJetError::InvalidInput(
                "Carreau liquid coefficients",
            ));
        }
        let reference = parameters.sample(0.0, None)?;
        let native = SlenderJetParameters {
            density_kg_m3: parameters.density_kg_m3,
            surface_tension_n_m: parameters.surface_tension_n_m,
            kinematic_viscosity_m2_s: positive(
                reference.dynamic_viscosity_pa_s() / parameters.density_kg_m3,
            )?,
        };
        Self::new_constitutive(
            native,
            Some(parameters),
            length_m,
            radii_m,
            velocities_m_s,
            limits,
            gate,
        )
    }

    /// Original law and fixed SI conditions, or None for the Newtonian path.
    #[must_use]
    pub const fn carreau_parameters(&self) -> Option<&CarreauSlenderJetParameters> {
        self.carreau.as_ref()
    }

    /// Source response at an accepted cell's actual axial velocity gradient.
    /// Returns None for a Newtonian jet; an out-of-range cell always refuses.
    /// A successful observation allocates nothing and changes neither state nor time.
    pub fn carreau_sample_at(
        &self,
        cell: usize,
    ) -> Result<Option<ThermalRheologyState>, NonlinearJetError> {
        let n = self.accepted.a.len();
        if cell >= n {
            return Err(NonlinearJetError::InvalidInput("jet rheology cell"));
        }
        let Some(parameters) = self.carreau.as_ref() else {
            return Ok(None);
        };
        let left = previous(cell, n);
        let velocity =
            |i| self.accepted.q[i] / average(self.accepted.a[i], self.accepted.a[next(i, n)]);
        let strain = finite((velocity(cell) - velocity(left)) / self.spacing_m)?;
        parameters.sample(strain, Some(cell)).map(Some)
    }
}

#[cfg(test)]
mod tests;
