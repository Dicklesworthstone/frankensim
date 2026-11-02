//! Temperature changes of an existing common- or mixed-mode cohesive interface.
//!
//! The declared constitutive ordering is coefficient change at fixed jump and
//! damage, followed by irreversible relaxation onto the new envelope at that
//! same jump. For `psi = (1-d) g_active.K.g_active/2 + Kc <g_n>_-^2/2`, the first
//! stage supplies signed parameter work `psi(new K, old d)-psi(old K, old d)`.
//! The second dissipates the released storage `psi(new K, old d)-psi(new K, d)`.
//! Neither stage does mechanical displacement work or rewrites past fracture
//! dissipation. Cooling or increased toughness cannot heal a damaged bond.
//!
//! This is an explicit staggered constitutive model, not a temperature-time
//! interpolation, thermal fatigue/cure law, heat partition or entropy model.
//! The caller owns the thermal source and resolves temperature/load histories.
//! A fixed-temperature law and its existing force/work integration stay intact.
//!
//! Mixed-mode transfer freezes the PHYSICAL active-separation direction `e`
//! at first initiation. At new coefficients its quadratic strength surface
//! gives `r0 = sqrt(Kn en^2+Kt |et|^2) / sqrt((Kn en/Tn)^2+(Kt |et|/Tt)^2)`.
//! The current weighting is `beta(T) = Kt |et|^2/(Kn en^2+Kt |et|^2)` and
//! `Gc(T) = GIc(T)+(GIIc(T)-GIc(T))*beta(T)^eta(T)`, `rf = 2 Gc/r0`.
//! Original onset beta remains a separate historical observation. Changing Kn
//! relative to Kt changes energy weighting without rotating the physical mode.
//! An uninitiated point first crossing the changed strength surface adopts the
//! held active gap's direction. No unseen temperature path is reconstructed.
//!
//! Temperature-dependent cohesive stiffness and fracture data are standard
//! inputs (Abaqus, "Contact Cohesive Behavior", 2025). The history-preserving
//! ordering and energy accounts above are the explicit model adopted here;
//! source data alone do not establish their calibration for a real adhesive.
//! https://docs.software.vt.edu/abaqusv2025/English/SIMACAEITNRefMap/simaitn-c-cohesivebehavior.htm

use super::mixed_mode::{coordinate, initiate, strength_coordinate};
use super::{BilinearCohesive, CohesiveError, CohesivePoint, active, finite, norm, positive};

/// Accepted fixed-jump temperature transfer and its surface-energy accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveTemperatureStep<'a> {
    /// Original separation and irreversible history in the supplied new law.
    pub point: CohesivePoint<'a>,
    /// Signed work supplied by changing coefficients, J/m2; not heat by itself.
    pub parameter_work_j_m2: f64,
    /// Nonnegative release from additional damage at fixed separation, J/m2.
    pub damage_dissipated_j_m2: f64,
    /// Complete storage change across both constitutive stages, J/m2.
    pub stored_change_j_m2: f64,
    /// Parameter work minus storage change minus new damage loss, J/m2.
    pub defect_j_m2: f64,
}

impl CohesivePoint<'_> {
    /// Transfer the SAME point to an explicitly supplied temperature law, keeping
    /// its common- or frozen-physical-mode envelope. The source adapter checks
    /// original card/claims and the temperature coordinate; this numerical
    /// operation does not confer source
    /// authority. The immutable input remains retryable on every refusal.
    ///
    /// The new envelope may increase damage at the held jump, but may never
    /// decrease it. The inverse damage map carries the retained unloading
    /// stiffness into the new onset/failure coordinates. Accumulated fracture
    /// dissipation is anchored separately, so changing Gc cannot rescale history.
    /// No time passes and dashpots do no work during this fixed-jump transfer.
    pub fn transfer_temperature<'b>(
        &self,
        law: &'b BilinearCohesive,
        energy_tolerance_j_m2: f64,
    ) -> Result<CohesiveTemperatureStep<'b>, CohesiveError> {
        positive(energy_tolerance_j_m2)?;
        if self.law.mixed_mode.is_some() != law.mixed_mode.is_some() {
            return Err(CohesiveError::Invalid(
                "temperature transfer cannot change the cohesive envelope model",
            ));
        }
        if norm(self.jump)? > law.maximum_separation {
            return Err(CohesiveError::SeparationLimit);
        }
        if self.law == law {
            return Ok(CohesiveTemperatureStep {
                point: CohesivePoint { law, ..*self },
                parameter_work_j_m2: 0.0,
                damage_dissipated_j_m2: 0.0,
                stored_change_j_m2: 0.0,
                defect_j_m2: 0.0,
            });
        }
        let mut point = CohesivePoint { law, ..*self };
        let before_damage = if let Some(p) = law.mixed_mode {
            let intact = self.mixed_integrity();
            // Viscosity-only/compression-only updates preserve the exact
            // existing envelope and damage coordinate without roundoff drift.
            if self.law.mixed_mode != law.mixed_mode {
                point.mixed_history = match self.mixed_history {
                    Some(history) => Some(history.at_parameters(p)?),
                    None if strength_coordinate(p, self.jump)? >= 1.0 => {
                        Some(initiate(p, [0.0; 3], self.jump)?.0)
                    }
                    None => None,
                };
                point.maximum = if let Some(history) = point.mixed_history {
                    equivalent_maximum(
                        history.onset_sqrt_energy,
                        history.failure_sqrt_energy,
                        intact,
                        |value| history.intact(value),
                    )?
                    .max(coordinate(p, self.jump)?)
                    .min(history.failure_sqrt_energy)
                } else {
                    0.0
                };
            }
            point.fracture_reference = point
                .mixed_history
                .map_or(0.0, |history| history.fracture(point.maximum));
            let g = active(self.jump);
            let recoverable = 0.5 * (p.normal_stiffness_pa_m * g[0]) * g[0]
                + 0.5 * (p.tangential_stiffness_pa_m * g[1]) * g[1]
                + 0.5 * (p.tangential_stiffness_pa_m * g[2]) * g[2];
            intact * recoverable
        } else {
            let intact = self.law.intact(self.maximum);
            let retained =
                equivalent_maximum(law.onset, law.failure, intact, |value| law.intact(value))?;
            point.maximum = retained.max(norm(active(self.jump))?).min(law.failure);
            point.fracture_reference = law.fracture(point.maximum);
            let stiffness = law.stiffness * intact;
            active(self.jump)
                .iter()
                .map(|&x| (0.5 * stiffness * x) * x)
                .sum::<f64>()
        };
        point.dissipated_reference = self.dissipated_j_m2();
        let compression = self.jump[0].min(0.0);
        let before_damage =
            finite(before_damage + (0.5 * law.compression_stiffness * compression) * compression)?;
        let old_storage = self.stored_energy_j_m2();
        let next_storage = point.stored_energy_j_m2();
        let released = finite(before_damage - next_storage)?;
        if released < 0.0 || point.damage() < self.damage() {
            return Err(CohesiveError::Invalid(
                "temperature transfer cannot heal damage or create negative fracture loss",
            ));
        }
        point.dissipated_reference = finite(point.dissipated_reference + released)?;
        let parameter_work = finite(before_damage - old_storage)?;
        let stored_change = finite(next_storage - old_storage)?;
        let defect = finite(parameter_work - stored_change - released)?;
        if defect.abs() > energy_tolerance_j_m2 {
            return Err(CohesiveError::EnergyDefect(defect));
        }
        Ok(CohesiveTemperatureStep {
            point,
            parameter_work_j_m2: parameter_work,
            damage_dissipated_j_m2: released,
            stored_change_j_m2: stored_change,
            defect_j_m2: defect,
        })
    }
}

fn equivalent_maximum(
    onset: f64,
    failure: f64,
    intact: f64,
    envelope_integrity: impl Fn(f64) -> f64,
) -> Result<f64, CohesiveError> {
    if intact == 1.0 {
        return Ok(0.0);
    }
    if intact == 0.0 {
        return Ok(failure);
    }
    // Invert intact = onset/r * (failure-r)/(failure-onset), avoiding
    // onset*failure overflow. Select the non-healing side of the inverse.
    let ratio = onset / failure;
    let mut maximum = positive(onset / (ratio + intact * (1.0 - ratio)))?
        .max(onset)
        .min(failure);
    for _ in 0..4 {
        if envelope_integrity(maximum) <= intact {
            return Ok(maximum);
        }
        maximum = maximum.next_up().min(failure);
    }
    Err(CohesiveError::Invalid(
        "nonrepresentable retained thermal damage",
    ))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod mixed_tests;
