//! Pure-pair-calibrated neutral phase exchange on the original virial potential.
//!
//! `mu_i - mu_s,i = RT[ln(c_i/csat_i) + (B c)_i - B_ii csat_i]`.
//! The offset retains each original binary saturation exactly. Cross-species
//! effects come only from the explicitly supplied symmetric interaction model;
//! neither a measured mixture fit nor concentrated-solution validity is inferred.
//! Pure crystals do not mix and solute partial volumes remain negligible.

use super::{KineticSolubilitySolution, MAX_OSMOTIC_SOLUTES, Result, invalid};
use crate::osmosis::{
    AqueousSolution, Sum,
    nonideal::{VirialSolutionModel, log_ratio},
};
use fs_math::det::{exp, ln, sqrt};
use fs_qty::{Amount, AmountConcentration, Temperature, Time, Volume};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::{SolubilityErrorV1, kinetics::SolubilityKineticsV1},
};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Calibration {
    pure_saturation: Box<[f64]>,
}

impl Calibration {
    fn new(models: &[SolubilityKineticsV1], thermodynamics: &VirialSolutionModel) -> Result<Self> {
        if models.len() != thermodynamics.species_count() {
            return Err(invalid(
                "virial and original kinetic species layouts must match",
            ));
        }
        let mut pure_saturation = Vec::with_capacity(models.len());
        for (model, &cap) in models
            .iter()
            .zip(thermodynamics.maximum_concentrations_mol_m3())
        {
            let pure = model
                .law()
                .saturation_concentration(Temperature::new(thermodynamics.temperature_k()))?
                .value();
            if pure > cap {
                return Err(invalid(
                    "original pure-pair saturation is outside the virial concentration support",
                ));
            }
            pure_saturation.push(pure);
        }
        Ok(Self {
            pure_saturation: pure_saturation.into_boxed_slice(),
        })
    }

    pub(super) fn heap_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + core::mem::size_of_val(self.pure_saturation.as_ref())
            + 2 * core::mem::size_of::<usize>()
            + core::mem::align_of::<Self>()
    }
}

impl KineticSolubilitySolution {
    /// Bind an explicitly qualified virial model to the SAME actual phases and
    /// original kinetic laws. Temperature, solvent molar volume and species
    /// order must match. Every original pure-pair saturation must lie inside the
    /// model's concentration box; no extrapolated binary calibration is used.
    ///
    /// The pure-solid reference is calibrated by `B_ii*csat_i`, preserving the
    /// original pure-pair saturation. Other dissolved species shift the
    /// conditional saturation through their original cross interactions. This
    /// is a neutral, convex single-solution model with separate pure crystals;
    /// it does not model ions, association, mixed crystals or a second liquid.
    ///
    /// Binding prepares shared immutable calibration and refreshes zero-time
    /// phase targets without moving either stock. Existing bindings cannot be
    /// replaced by a different thermodynamic source. A fixed-temperature
    /// virial fit supplies neither excess calorics nor a temperature derivative.
    pub fn with_nonideal_thermodynamics(
        mut self,
        thermodynamics: Arc<VirialSolutionModel>,
    ) -> Result<Self> {
        if self.solution.temperature_k().to_bits() != thermodynamics.temperature_k().to_bits()
            || self.solution.solvent_molar_volume_m3_mol().to_bits()
                != thermodynamics.solvent_molar_volume_m3_mol().to_bits()
            || self
                .solution
                .nonideal_model()
                .is_some_and(|original| original != thermodynamics.as_ref())
        {
            return Err(invalid(
                "nonideal phase binding must retain the original temperature, solvent convention and thermodynamic source",
            ));
        }
        let calibration = Calibration::new(&self.models, &thermodynamics)?;
        let solution = AqueousSolution::nonideal(
            thermodynamics,
            self.solution.solvent_volume_m3(),
            self.solution.solute_moles().to_vec(),
        )?;
        self.solution = solution;
        self.calibration = Some(Arc::new(calibration));
        self.advance(Time::new(0.0), || false)?;
        Ok(self)
    }

    /// Saturation of one neutral component while every OTHER supplied
    /// concentration is held fixed [mol/m³]. The original convex box guarantees
    /// a strictly increasing chemical potential along this coordinate. A root
    /// outside that box, or below the smallest positive value, refuses instead
    /// of flooring a trace or extrapolating the material source.
    pub fn conditional_saturation(
        &self,
        component: usize,
        concentrations_mol_m3: &[f64],
    ) -> Result<AmountConcentration> {
        if component >= self.models.len() {
            return Err(invalid(
                "conditional saturation component is outside the original layout",
            ));
        }
        self.solution
            .osmotic_pressure_at(1.0, concentrations_mol_m3)?;
        self.conditional_saturation_unchecked(component, concentrations_mol_m3)
    }

    pub(super) fn phase_concentrations(
        &self,
        volume: f64,
        dissolved: &[f64],
    ) -> Result<[f64; MAX_OSMOTIC_SOLUTES]> {
        let mut result = [0.0; MAX_OSMOTIC_SOLUTES];
        if self.calibration.is_some() {
            for (out, &amount) in result.iter_mut().zip(dissolved) {
                *out = finite(amount / volume)?;
            }
        }
        Ok(result)
    }

    pub(super) fn conditional_saturation_unchecked(
        &self,
        component: usize,
        concentrations: &[f64],
    ) -> Result<AmountConcentration> {
        let Some(calibration) = &self.calibration else {
            return Ok(self.models[component]
                .law()
                .saturation_concentration(Temperature::new(self.solution.temperature_k()))?);
        };
        let thermodynamics = self
            .solution
            .nonideal_model()
            .ok_or_else(|| invalid("missing original virial source for phase calibration"))?;
        let n = self.models.len();
        let row = &thermodynamics.interaction_matrix_m3_mol()[component * n..(component + 1) * n];
        let mut other = Sum::default();
        for (j, (&bij, &concentration)) in row.iter().zip(concentrations).enumerate() {
            if j != component {
                other.add(finite(bij * concentration)?);
            }
        }
        Ok(AmountConcentration::new(conditional_root(
            calibration.pure_saturation[component],
            row[component],
            other.finish()?,
            thermodynamics.maximum_concentrations_mol_m3()[component],
        )?))
    }

    /// Aqueous-minus-original-solid chemical potentials [J/mol]. Exactly absent
    /// species return `None`, the ideal negative-infinite limit. The complete
    /// output stays unchanged on refusal, and the common concentration reference
    /// of the virial owner cancels before arithmetic.
    pub fn solid_relative_chemical_potentials_into(
        &self,
        concentrations_mol_m3: &[f64],
        output: &mut [Option<f64>],
    ) -> Result<()> {
        if output.len() != self.models.len() {
            return Err(invalid("solid-relative chemical potential output layout"));
        }
        self.solution
            .osmotic_pressure_at(1.0, concentrations_mol_m3)?;
        let mut candidate = [None; MAX_OSMOTIC_SOLUTES];
        let temperature = Temperature::new(self.solution.temperature_k());
        for (i, (&concentration, model)) in concentrations_mol_m3
            .iter()
            .zip(self.models.iter())
            .enumerate()
        {
            if concentration == 0.0 {
                continue;
            }
            let mut potential = Sum::default();
            potential.add(
                model
                    .law()
                    .solid_relative_chemical_potential_j_mol(
                        temperature,
                        AmountConcentration::new(concentration),
                    )?
                    .expect("positive concentration has a finite ideal chemical potential"),
            );
            if let Some(calibration) = &self.calibration {
                let source = self.solution.nonideal_model().ok_or_else(|| {
                    invalid("missing original virial source for phase calibration")
                })?;
                let row = &source.interaction_matrix_m3_mol()
                    [i * self.models.len()..(i + 1) * self.models.len()];
                for (j, (&bij, &other)) in row.iter().zip(concentrations_mol_m3).enumerate() {
                    let coordinate = if i == j {
                        other - calibration.pure_saturation[i]
                    } else {
                        other
                    };
                    potential.add(energy_product(R * temperature.value(), bij, coordinate)?);
                }
            }
            candidate[i] = Some(potential.finish()?);
        }
        output.copy_from_slice(&candidate[..output.len()]);
        Ok(())
    }

    /// Isothermal solid-relative chemical work [J] between actual inventories,
    /// allowing independent solvent-volume changes. May be used for any cell
    /// sharing these original phase laws and aqueous thermodynamics. Complete
    /// dissolved-plus-solid component conservation must be checked separately
    /// before cancelling the common pure-solid datum across a closed system.
    ///
    /// Uses each original pair's stable ideal work and the SAME virial
    /// potential's bilinear excess difference, minus `RT B_ii csat_i dn_i`.
    /// This calibration is linear in component moles and adds no pressure term
    /// at fixed amounts. It is chemical work, not a heat-of-solution balance.
    pub fn solid_relative_inventory_free_energy_change_j(
        &self,
        before_volume_m3: f64,
        before_amounts_mol: &[f64],
        after_volume_m3: f64,
        after_amounts_mol: &[f64],
    ) -> Result<f64> {
        self.solution
            .osmotic_pressure_at(before_volume_m3, before_amounts_mol)?;
        self.solution
            .osmotic_pressure_at(after_volume_m3, after_amounts_mol)?;
        let temperature = Temperature::new(self.solution.temperature_k());
        let mut work = Sum::default();
        for ((model, &before), &after) in self
            .models
            .iter()
            .zip(before_amounts_mol)
            .zip(after_amounts_mol)
        {
            work.add(model.law().solid_relative_inventory_free_energy_change_j(
                temperature,
                Volume::new(before_volume_m3),
                Amount::new(before),
                Volume::new(after_volume_m3),
                Amount::new(after),
            )?);
        }
        let Some(calibration) = &self.calibration else {
            return Ok(work.finish()?);
        };
        let thermodynamics = self
            .solution
            .nonideal_model()
            .ok_or_else(|| invalid("missing original virial source for phase calibration"))?;
        let n = self.models.len();
        let before = self.phase_concentrations(before_volume_m3, before_amounts_mol)?;
        let after = self.phase_concentrations(after_volume_m3, after_amounts_mol)?;
        let relative_volume_change = (after_volume_m3 - before_volume_m3) / after_volume_m3;
        let nearby_volumes = (-1.0..=0.5).contains(&relative_volume_change);
        let mut change_c = [0.0; MAX_OSMOTIC_SOLUTES];
        for i in 0..n {
            change_c[i] = finite(
                if nearby_volumes && before_amounts_mol[i] > 0.0 && after_amounts_mol[i] > 0.0 {
                    (after_amounts_mol[i] - before_amounts_mol[i]) / after_volume_m3
                        - before[i] * relative_volume_change
                } else {
                    after[i] - before[i]
                },
            )?;
        }
        for i in 0..n {
            let delta_n = after_amounts_mol[i] - before_amounts_mol[i];
            let row = &thermodynamics.interaction_matrix_m3_mol()[i * n..(i + 1) * n];
            let mut mean_interaction = Sum::default();
            let mut changed_interaction = Sum::default();
            for (j, &bij) in row.iter().enumerate() {
                let mean = before[j] + 0.5 * (after[j] - before[j]);
                if before_volume_m3 == after_volume_m3 {
                    mean_interaction.add(finite(
                        bij * if i == j {
                            mean - calibration.pure_saturation[i]
                        } else {
                            mean
                        },
                    )?);
                } else {
                    mean_interaction.add(finite(bij * mean)?);
                    changed_interaction.add(finite(bij * change_c[j])?);
                }
            }
            let rt = R * temperature.value();
            if before_volume_m3 == after_volume_m3 {
                work.add(energy_product(rt, delta_n, mean_interaction.finish()?)?);
            } else {
                let reference = finite(row[i] * calibration.pure_saturation[i])?;
                work.add(energy_product(
                    rt,
                    delta_n,
                    0.5 * mean_interaction.finish()? - reference,
                )?);
                let mean_n = before_amounts_mol[i] + 0.5 * delta_n;
                work.add(energy_product(
                    rt,
                    mean_n,
                    0.5 * changed_interaction.finish()?,
                )?);
            }
        }
        Ok(work.finish()?)
    }
}

// A monotone scalar equilibrium on the ORIGINAL convex concentration box.
// Alternating safeguarded Newton and geometric/arithmetic bisection bounds
// work even for trace concentrations spanning the complete f64 exponent range.
fn conditional_root(pure: f64, diagonal: f64, other: f64, cap: f64) -> Result<f64> {
    if other == 0.0 {
        return Ok(pure);
    }
    let force = |c: f64| -> Result<f64> {
        let mut sum = Sum::default();
        sum.add(log_ratio(c, pure));
        sum.add(finite(diagonal * (c - pure))?);
        sum.add(other);
        Ok(sum.finish()?)
    };
    if diagonal == 0.0 {
        let direct = pure * exp(-other);
        let root = if direct.is_normal() {
            direct
        } else {
            exp(ln(pure) - other)
        };
        if !root.is_finite() || root <= 0.0 || root > cap {
            return Err(invalid(
                "conditional saturation is outside the positive representable virial support",
            ));
        }
        return Ok(root);
    }
    let minimum = f64::from_bits(1);
    let (mut low, mut high, mut f_low, mut f_high) = if other > 0.0 {
        (minimum, pure, force(minimum)?, other)
    } else {
        (pure, cap, other, force(cap)?)
    };
    if f_low > 0.0 || f_high < 0.0 {
        return Err(invalid(
            "conditional saturation is outside the positive representable virial support",
        ));
    }
    if f_low == 0.0 {
        return Ok(low);
    }
    if f_high == 0.0 {
        return Ok(high);
    }
    let estimate = pure * exp(-other / (1.0 + diagonal * pure));
    let mut candidate = if estimate.is_finite() && estimate > low && estimate < high {
        estimate
    } else {
        midpoint(low, high)
    };
    for iteration in 0..128 {
        let value = force(candidate)?;
        if value == 0.0 {
            return Ok(candidate);
        }
        if value < 0.0 {
            low = candidate;
            f_low = value;
        } else {
            high = candidate;
            f_high = value;
        }
        if high.to_bits() <= low.to_bits() + 1 || high - low <= 16.0 * f64::EPSILON * high {
            return Ok(if f_low.abs() <= f_high.abs() {
                low
            } else {
                high
            });
        }
        let derivative = 1.0 + diagonal * candidate;
        let newton = candidate - candidate * (value / derivative);
        candidate = if iteration % 2 == 0 && newton.is_finite() && newton > low && newton < high {
            newton
        } else {
            midpoint(low, high)
        };
        if candidate <= low || candidate >= high {
            return Ok(if f_low.abs() <= f_high.abs() {
                low
            } else {
                high
            });
        }
    }
    Err(invalid(
        "conditional saturation exhausted its bounded scalar solve",
    ))
}

fn midpoint(low: f64, high: f64) -> f64 {
    if high - low <= low {
        low + 0.5 * (high - low)
    } else {
        sqrt(low) * sqrt(high)
    }
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolubilityErrorV1::Unrepresentable("calibrated nonideal solubility quantity").into())
    }
}

fn energy_product(rt: f64, amount: f64, factor: f64) -> Result<f64> {
    let ordinary = (amount * factor) * rt;
    if ordinary.is_finite() && (ordinary != 0.0 || amount == 0.0 || factor == 0.0) {
        return Ok(ordinary);
    }
    for alternative in [(rt * amount) * factor, (rt * factor) * amount] {
        if alternative.is_finite() && alternative != 0.0 {
            return Ok(alternative);
        }
    }
    finite(ordinary)
}

#[cfg(test)]
mod tests;
