//! Original finite-solution scalar balance, shared by one- and two-bath fields.

use super::*;
use crate::moisture::reservoir::water_sum;

/// One finite bath's original-to-final solvent account.
///
/// `water_to_solid_kg` contains only native face transfers assigned to THIS
/// bath. With several baths it is not the specimen's total storage change.
#[derive(Clone, Debug)]
pub struct SolutionMoistureExchange {
    pub water_activity: [f64; 2],
    /// Thermodynamic diagnostic; it is not an applied mechanical traction.
    pub osmotic_pressure_pa: [f64; 2],
    pub solvent_mass_kg: [f64; 2],
    pub water_to_solid_kg: f64,
    /// Actual represented solvent-mass change plus its assigned face transfers.
    pub water_balance_residual_kg: f64,
    /// Interior scalar trials; initial and bracket-endpoint trials are additional.
    pub iterations: usize,
}

struct Trial<T> {
    reservoir: FiniteMoistureSolution,
    value: T,
    water_to_solid_kg: f64,
    residual_kg: f64,
}

impl FiniteMoistureSolution {
    /// Close this original bath against an actual native signed water exchange.
    /// All calls restart its original solute/phase stocks. The callback owns
    /// the field proposal and supplies only this bath's native face transfers.
    pub(super) fn advance_exchange<T>(
        &self,
        dt_s: f64,
        sealed: bool,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
        mut propose: impl FnMut(f64) -> Result<(T, f64), SolutionMoistureError>,
    ) -> Result<(Self, T, SolutionMoistureExchange), SolutionMoistureError> {
        poll(cx)?;
        let initial_volume = self.solution().solvent_volume_m3();
        let [mut lower, mut upper] = controls.solvent_volume_bracket_m3;
        if !positive(lower)
            || !positive(upper)
            || lower >= upper
            || initial_volume < lower
            || initial_volume > upper
            || !positive(controls.mass_tolerance_kg)
            || controls.maximum_iterations == 0
            || controls.maximum_iterations > 256
        {
            return Err(SolutionMoistureError::Invalid(
                "positive source-supported volume bracket containing the start, finite mass tolerance and bounded iterations required",
            ));
        }
        let mut evaluate = |volume: f64| -> Result<Trial<T>, SolutionMoistureError> {
            poll(cx)?;
            let reservoir = self.at_endpoint_volume(volume, dt_s, || cx.checkpoint().is_err())?;
            let (value, water_to_solid_kg) = propose(reservoir.solution().water_activity()?)?;
            let residual_kg =
                (volume - initial_volume) * self.solvent_density_kg_m3 + water_to_solid_kg;
            if !residual_kg.is_finite() || !water_to_solid_kg.is_finite() {
                return Err(SolutionMoistureError::Invalid(
                    "finite solution water balance is unrepresentable",
                ));
            }
            Ok(Trial {
                reservoir,
                value,
                water_to_solid_kg,
                residual_kg,
            })
        };

        let mut iterations = 0;
        let start = evaluate(initial_volume)?;
        let accepted = if start.residual_kg.abs() <= controls.mass_tolerance_kg {
            start
        } else if sealed {
            return Err(SolutionMoistureError::MassBalance {
                residual_kg: start.residual_kg,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        } else if self.kinetic_solubility().is_none()
            && self
                .solution()
                .solute_moles()
                .iter()
                .all(|&amount| amount == 0.0)
        {
            // Original pure-solvent shortcut: activity is one at every positive
            // volume, so the field proposal is already correct. Kinetic solid
            // stocks cannot use this shortcut because dissolution changes a_w.
            let volume = initial_volume - start.water_to_solid_kg / self.solvent_density_kg_m3;
            if !positive(volume) || volume < lower || volume > upper {
                return Err(SolutionMoistureError::Unbracketed {
                    residual_kg: [
                        (lower - initial_volume) * self.solvent_density_kg_m3
                            + start.water_to_solid_kg,
                        (upper - initial_volume) * self.solvent_density_kg_m3
                            + start.water_to_solid_kg,
                    ],
                });
            }
            Trial {
                reservoir: self.at_endpoint_volume(volume, dt_s, || cx.checkpoint().is_err())?,
                residual_kg: (volume - initial_volume) * self.solvent_density_kg_m3
                    + start.water_to_solid_kg,
                ..start
            }
        } else {
            drop(start);
            let low = evaluate(lower)?;
            let high = evaluate(upper)?;
            if low.residual_kg.abs() <= controls.mass_tolerance_kg {
                low
            } else if high.residual_kg.abs() <= controls.mass_tolerance_kg {
                high
            } else {
                if low.residual_kg.is_sign_positive() == high.residual_kg.is_sign_positive() {
                    return Err(SolutionMoistureError::Unbracketed {
                        residual_kg: [low.residual_kg, high.residual_kg],
                    });
                }
                let mut low_residual = low.residual_kg;
                let mut high_residual = high.residual_kg;
                // Retain bracket signs, not two complete extra field proposals.
                drop(low);
                drop(high);
                loop {
                    poll(cx)?;
                    let scale = low_residual.abs().max(high_residual.abs());
                    let low_weight = low_residual.abs() / scale;
                    let high_weight = high_residual.abs() / scale;
                    let fraction = (low_weight / (low_weight + high_weight)).clamp(0.05, 0.95);
                    let volume = lower + fraction * (upper - lower);
                    if volume <= lower || volume >= upper {
                        return Err(SolutionMoistureError::NotConverged {
                            residual_kg: low_residual.abs().min(high_residual.abs()),
                            iterations,
                        });
                    }
                    let trial = evaluate(volume)?;
                    iterations += 1;
                    if trial.residual_kg.abs() <= controls.mass_tolerance_kg {
                        break trial;
                    }
                    if iterations == controls.maximum_iterations {
                        return Err(SolutionMoistureError::NotConverged {
                            residual_kg: trial.residual_kg,
                            iterations,
                        });
                    }
                    if trial.residual_kg.is_sign_positive() == low_residual.is_sign_positive() {
                        lower = volume;
                        low_residual = trial.residual_kg;
                    } else {
                        upper = volume;
                        high_residual = trial.residual_kg;
                    }
                }
            }
        };

        let reservoir = accepted.reservoir;
        self.admit_exchanged_components(&reservoir, cx)?;
        let solvent_mass_kg = [self.solvent_mass_kg(), reservoir.solvent_mass_kg()];
        let residual = solvent_mass_kg[1] - solvent_mass_kg[0] + accepted.water_to_solid_kg;
        for residual_kg in [accepted.residual_kg, residual] {
            if !residual_kg.is_finite() || residual_kg.abs() > controls.mass_tolerance_kg {
                return Err(SolutionMoistureError::MassBalance {
                    residual_kg,
                    tolerance_kg: controls.mass_tolerance_kg,
                });
            }
        }
        let step = SolutionMoistureExchange {
            water_activity: [
                self.solution().water_activity()?,
                reservoir.solution().water_activity()?,
            ],
            osmotic_pressure_pa: [
                self.solution().osmotic_pressure_pa()?,
                reservoir.solution().osmotic_pressure_pa()?,
            ],
            solvent_mass_kg,
            water_to_solid_kg: accepted.water_to_solid_kg,
            water_balance_residual_kg: residual,
            iterations,
        };
        poll(cx)?;
        Ok((reservoir, accepted.value, step))
    }

    // These original single-bath checks also apply independently to each bath
    // of a pair. Finite phases never become a different bath's solute stock.
    pub(in crate::moisture) fn admit_exchanged_components(
        &self,
        reservoir: &Self,
        cx: &Cx<'_>,
    ) -> Result<(), SolutionMoistureError> {
        if let Some(initial) = self.solubility() {
            let endpoint = reservoir
                .solubility()
                .ok_or(SolutionMoistureError::Invalid(
                    "finite solid inventories were not retained by the moisture candidate",
                ))?;
            if initial.components().len() != endpoint.components().len() {
                return Err(SolutionMoistureError::Invalid(
                    "finite solubility component layout changed during moisture exchange",
                ));
            }
            for ((original, component), partition) in initial
                .components()
                .iter()
                .zip(endpoint.components())
                .zip(endpoint.phase_partitions())
            {
                let total = original.total().value();
                let residual =
                    water_sum([partition.dissolved.value(), partition.solid.value(), -total]);
                if component != original
                    || !residual.is_finite()
                    || residual.abs() > 8.0 * f64::EPSILON * total
                {
                    return Err(SolutionMoistureError::Invalid(
                        "dissolved plus solid solute failed its original component balance",
                    ));
                }
            }
        }
        if let Some(initial) = self.kinetic_solubility() {
            let endpoint = reservoir
                .kinetic_solubility()
                .ok_or(SolutionMoistureError::Invalid(
                    "kinetic phase inventories were not retained by the moisture candidate",
                ))?;
            if initial.models() != endpoint.models()
                || initial.last_steps().len() != endpoint.last_steps().len()
            {
                return Err(SolutionMoistureError::Invalid(
                    "kinetic source laws or component layout changed during moisture exchange",
                ));
            }
            for (before, after) in initial.last_steps().iter().zip(endpoint.last_steps()) {
                poll(cx)?;
                let total = before.state().total_solute().value();
                let phase = after.state();
                let residual =
                    water_sum([phase.dissolved().value(), phase.solid().value(), -total]);
                if phase.total_solute() != before.state().total_solute()
                    || !residual.is_finite()
                    || residual.abs() > 8.0 * f64::EPSILON * total
                {
                    return Err(SolutionMoistureError::Invalid(
                        "kinetic dissolved plus solid solute failed its original component balance",
                    ));
                }
            }
        }
        if let Some(initial) = self.electrolyte_solubility() {
            let endpoint =
                reservoir
                    .electrolyte_solubility()
                    .ok_or(SolutionMoistureError::Invalid(
                        "source-bound ionic phases were not retained by the moisture candidate",
                    ))?;
            if initial.component_total() != endpoint.component_total()
                || initial.family_totals() != endpoint.family_totals()
                || initial.spectators() != endpoint.spectators()
                || initial.species() != endpoint.species()
                || initial.solvent() != endpoint.solvent()
            {
                return Err(SolutionMoistureError::Invalid(
                    "ionic analytical inventories or species layout changed during solvent exchange",
                ));
            }
        }
        Ok(())
    }
}
