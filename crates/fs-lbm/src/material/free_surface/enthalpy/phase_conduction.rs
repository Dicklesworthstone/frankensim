//! The retained phase law on the original, strictly liquid VOF caloric state.

use super::*;

/// Source conductivity and the explicitly declared phase correction actually used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceConductivity {
    /// Original pinned material source at the actual temperature/rate [W/(m K)].
    pub source_w_m_k: f64,
    /// Declared dimensionless liquid-mass-fraction multiplier; one if absent.
    pub phase_multiplier: f64,
    /// Product used by the original conductance and thermal CFL [W/(m K)].
    pub effective_w_m_k: f64,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Retain a declared phase-conductivity law on this original liquid owner.
    ///
    /// Every admitted cell remains exactly liquid according to its original h
    /// chart, so the same shared law is evaluated at liquid mass fraction one.
    /// This is not a new partial-melt flow or a measured mixture rule. Original
    /// material snapshots keep their raw source conductivity and receipts.
    /// Bulk, finite-body, gas-film and solved-surface conduction all use the
    /// effective coefficient, including their original thermal CFL checks.
    ///
    /// This consumed builder replaces any prior law. All source-state and
    /// representability checks precede attachment; later endpoint refusal
    /// remains part of the original all-or-nothing native step.
    pub fn with_phase_conductivity(
        mut self,
        law: LiquidMassFractionConductivity,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let density = self.flow.scaling.material().density_kg_m3();
        for (location, material) in &self.states {
            poll(gate)?;
            let h = *self.h.get(location).ok_or(Error::Invalid(
                "phase conductivity requires the original enthalpy history",
            ))?;
            buoyancy::liquid_temperature(self.caloric, h, density)?;
            resolve(material, Some(law))?;
        }
        poll(gate)?;
        self.phase_conduction = Some(law);
        Ok(self)
    }

    /// Original declared law, including its retained solid endpoint.
    #[must_use]
    pub const fn phase_conductivity_law(&self) -> Option<LiquidMassFractionConductivity> {
        self.phase_conduction
    }

    /// Source and effective conductivity at an accepted native wet cell.
    pub fn conductivity_at(
        &self,
        coordinates: [usize; 3],
    ) -> Result<Option<FreeSurfaceConductivity>, Error> {
        self.material_at(coordinates)
            .map(|state| resolve(state, self.phase_conduction))
            .transpose()
    }
}

/// Callers obtain `material` only after the original pure-liquid caloric
/// admission (including the rheological source sampler). No independent phase
/// coordinate, caloric model or altered source receipt is introduced here.
pub(super) fn resolve(
    material: &ResolvedLiquidState,
    law: Option<LiquidMassFractionConductivity>,
) -> Result<FreeSurfaceConductivity, Error> {
    let source_w_m_k = material.thermal_conductivity_w_m_k();
    let phase_multiplier = match law {
        Some(law) => law.multiplier_at(1.0).map_err(|_| {
            Error::Invalid("liquid phase-conductivity multiplier is not representable")
        })?,
        None => 1.0,
    };
    let effective_w_m_k = source_w_m_k * phase_multiplier;
    if !effective_w_m_k.is_finite() || effective_w_m_k <= 0.0 {
        return Err(Error::Invalid(
            "effective liquid conductivity must be finite and positive",
        ));
    }
    Ok(FreeSurfaceConductivity {
        source_w_m_k,
        phase_multiplier,
        effective_w_m_k,
    })
}
