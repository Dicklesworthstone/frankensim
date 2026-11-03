//! Original gas calorics and transmission on the shared drifting-Maxwellian law.
//!
//! The outward normal points from the receiving vessel into the source gas.
//! Source bulk velocity is relative to that aperture. Incoming particles carry
//! the shared half-space conditional VECTOR velocity, including tangential
//! drift, and original `h_i(T)` plus its shared crossing-energy correction.
//! No extra thermalization heat or continuum pressure thrust is added.
//!
//! Gaussian moments and positive tails remain owned by
//! `fs_thermal::sublimation::molecular::drift::incident_maxwellian`. This adapter
//! only admits original material sources, applies the original transmission and
//! exposed area once, and retains the original rest-incidence beam separately.
//! A moving consumer must resample the relative drift as its velocity changes.

use super::momentum::EffusiveMomentumInterface;
use super::{
    ClosedGasParcel, EffusiveBeam, GasEnergyError, R, SpeciesId, invalid, poll,
    positive_product_ratio, signed_product, source_error, valid_provenance,
};
use fs_math::dd::Dd;
use fs_thermal::sublimation::molecular::drift::incident_maxwellian;

/// Reciprocal finite molecular streams through one shared moving aperture.
pub mod exchange;

/// Additional authored support for drift through one original aperture law.
///
/// The declaration must cover locally undisturbed Maxwellian incidence and
/// collisionless, angle/energy-independent, nondeflecting transmission over
/// every original gas mixture and T/p state AND the declared relative-speed
/// ball. Finite thickness, grazing incidence, source exposure and collisions
/// between counterstreaming molecules must be covered by that qualification;
/// the original rest-frame thinness/mean-free-path gates do not infer it.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveDriftInterface {
    momentum: EffusiveMomentumInterface,
    maximum_relative_speed_m_s: f64,
    provenance: String,
}
impl EffusiveDriftInterface {
    pub fn try_new(
        momentum: EffusiveMomentumInterface,
        maximum_relative_speed_m_s: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = provenance.into();
        if !maximum_relative_speed_m_s.is_finite()
            || maximum_relative_speed_m_s < 0.0
            || !valid_provenance(&provenance)
        {
            return Err(invalid(
                "drifting effusion needs bounded speed support and explicit qualification provenance",
            ));
        }
        Ok(Self {
            momentum,
            maximum_relative_speed_m_s,
            provenance,
        })
    }

    pub const fn momentum_interface(&self) -> &EffusiveMomentumInterface {
        &self.momentum
    }
    pub const fn maximum_relative_speed_m_s(&self) -> f64 {
        self.maximum_relative_speed_m_s
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Admit an actual relative velocity, including an enclosing solver's
    /// accepted endpoint. Zero speed support admits only exactly zero drift.
    pub fn admit_relative_velocity(&self, velocity: [f64; 3]) -> Result<(), GasEnergyError> {
        let limit = self.maximum_relative_speed_m_s;
        if velocity.iter().any(|v| !v.is_finite() || v.abs() > limit) {
            return Err(invalid(
                "relative gas drift is outside authored effusive support",
            ));
        }
        if limit == 0.0 {
            return Ok(());
        }
        // Normalize before squaring: large finite Cartesian components must
        // not overflow a speed test, nor may subnormal limits admit diagonals.
        let squared = velocity.into_iter().fold(Dd::ZERO, |sum, value| {
            let scaled = Dd::from_f64(value / limit);
            sum + scaled * scaled
        });
        if Dd::ONE.lt(squared) {
            return Err(invalid(
                "relative gas drift is outside authored effusive support",
            ));
        }
        Ok(())
    }
}

/// One original species' incident beam and its distinct drift-dependent moments.
#[derive(Clone, Debug, PartialEq)]
pub struct DriftingEffusiveSpeciesFlux {
    pub species: SpeciesId,
    /// Original source vacuum beam, interpreted as zero-drift incoming incidence.
    pub baseline_incoming: EffusiveBeam,
    /// Actual gross incoming beam, with original h and molecular correction apart.
    pub incoming: EffusiveBeam,
    /// Signed change retained even when adding it to the baseline would round
    /// it away. Do not reconstruct a far-tail gross flux by baseline + change.
    pub incoming_mass_rate_change_from_zero_drift_kg_s: f64,
    /// Conditional full velocity relative to the aperture. None for inactive
    /// species/transmissions; the normal component of an active beam is inward.
    pub mean_relative_velocity_m_s: Option<[f64; 3]>,
}
impl DriftingEffusiveSpeciesFlux {
    /// Momentum from the mass ACTUALLY admitted by the finite inventory owner.
    /// This is an incoming vector impulse, not an independently rounded rate*dt.
    pub fn impulse_for_mass_kg(&self, mass_kg: f64) -> Result<[f64; 3], GasEnergyError> {
        if !mass_kg.is_finite() || mass_kg < 0.0 {
            return Err(invalid(
                "drifting effusive impulse needs nonnegative accepted mass",
            ));
        }
        if mass_kg == 0.0 {
            return Ok([0.0; 3]);
        }
        let velocity = self.mean_relative_velocity_m_s.ok_or(invalid(
            "inactive drifting beam cannot carry an accepted impulse",
        ))?;
        Ok([
            signed_product(mass_kg, velocity[0])?,
            signed_product(mass_kg, velocity[1])?,
            signed_product(mass_kg, velocity[2])?,
        ])
    }
}

/// Source-bearing instantaneous aperture incidence, before any finite update.
#[derive(Clone, Debug, PartialEq)]
pub struct DriftingEffusiveFlux {
    pub interface: EffusiveDriftInterface,
    /// Canonical original source species, empty for zero geometric area.
    pub species: Vec<DriftingEffusiveSpeciesFlux>,
    pub relative_bulk_velocity_m_s: [f64; 3],
    pub outward_unit_normal: [f64; 3],
}

fn admit_normal(normal: [f64; 3]) -> Result<(), GasEnergyError> {
    let tolerance = 64.0 * f64::EPSILON;
    if normal
        .iter()
        .any(|v| !v.is_finite() || v.abs() > 1.0 + tolerance)
    {
        return Err(invalid(
            "drifting effusion needs a finite unit outward normal",
        ));
    }
    let squared = normal.into_iter().fold(Dd::ZERO, |sum, value| {
        let value = Dd::from_f64(value);
        sum + value * value
    });
    if (squared.to_f64() - 1.0).abs() > tolerance {
        return Err(invalid(
            "drifting effusion needs a finite unit outward normal",
        ));
    }
    Ok(())
}

impl ClosedGasParcel {
    /// Sample a source gas through an independently drift-qualified aperture.
    ///
    /// Supply source velocity MINUS aperture velocity in one common inertial
    /// basis. Zero area is exactly inactive and does not inspect source, drift
    /// or normal qualification. An active zero-transmission path still admits
    /// its original T/p/catalog and drift support, but carries no beam moments.
    ///
    /// The original vacuum sampler owns source h, zero-drift mass and support.
    /// The shared kinetic owner supplies full conditional velocity and the
    /// correction to that SAME h. Internal modes and formation energy are kept
    /// once, with no constant-Cv approximation or energy-derived directed speed.
    #[allow(clippy::too_many_arguments)]
    pub fn drifting_effusive_flux(
        &self,
        interface: &EffusiveDriftInterface,
        area_m2: f64,
        relative_bulk_velocity_m_s: [f64; 3],
        outward_unit_normal: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DriftingEffusiveFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        if !area_m2.is_finite() || area_m2 < 0.0 {
            return Err(invalid(
                "drifting effusive area must be finite and nonnegative",
            ));
        }
        let mut result = DriftingEffusiveFlux {
            interface: interface.clone(),
            species: Vec::new(),
            relative_bulk_velocity_m_s,
            outward_unit_normal,
        };
        if area_m2 == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        interface.admit_relative_velocity(relative_bulk_velocity_m_s)?;
        admit_normal(outward_unit_normal)?;
        let original = interface.momentum_interface().interface();
        let baseline = self.effusive_flux_to_vacuum(original, area_m2, &mut cancelled)?;
        let state = self.state().state();
        let mixture = self.model().thermodynamics();
        result.species.reserve(baseline.species.len());
        for ((row, source), &fraction) in baseline
            .species
            .into_iter()
            .zip(mixture.components())
            .zip(mixture.mole_composition().fractions())
        {
            poll(&mut cancelled)?;
            let mut incoming = DriftingEffusiveSpeciesFlux {
                species: row.species,
                baseline_incoming: row.outward,
                incoming: row.outward,
                incoming_mass_rate_change_from_zero_drift_kg_s: 0.0,
                mean_relative_velocity_m_s: None,
            };
            if row.outward.mass_rate_kg_s == 0.0 {
                result.species.push(incoming);
                continue;
            }
            let index = original
                .transmissions()
                .binary_search_by(|t| t.species().cmp(source.species()))
                .map_err(|_| invalid("drifting aperture lost its original species transmission"))?;
            let tau = original.transmissions()[index].probability();
            let molar_mass = source.molar_mass().value();
            let specific_r = positive_product_ratio(&[R], &[molar_mass])?;
            let density = positive_product_ratio(
                &[fraction, state.pressure, molar_mass],
                &[R, state.temperature],
            )?;
            let incidence = incident_maxwellian(
                density,
                specific_r,
                state.temperature,
                relative_bulk_velocity_m_s,
                outward_unit_normal,
            )
            .map_err(source_error)?;
            let moments = incidence.moments.ok_or(invalid(
                "active drifting source lost its conditional molecular moments",
            ))?;
            incoming.mean_relative_velocity_m_s = Some(moments.mean_relative_velocity_m_s);
            let change = incidence.incoming_mass_flux_change_from_zero_drift_kg_m2_s;
            if change != 0.0 {
                incoming.incoming_mass_rate_change_from_zero_drift_kg_s =
                    positive_product_ratio(&[area_m2, tau, change.abs()], &[])?.copysign(change);
                // Direct gross incidence preserves rare counterstreaming tails.
                // Recovering this from baseline+change could erase a real beam.
                incoming.incoming.mass_rate_kg_s = positive_product_ratio(
                    &[area_m2, tau, incidence.incoming_mass_flux_kg_m2_s],
                    &[],
                )?;
                incoming.incoming.molar_rate_mol_s =
                    positive_product_ratio(&[incoming.incoming.mass_rate_kg_s], &[molar_mass])?;
            }
            if relative_bulk_velocity_m_s != [0.0; 3] {
                let correction = if change == 0.0 {
                    // Keep the original rest thermal coordinate exactly. A
                    // purely tangential drift only adds its kinetic energy;
                    // resampling the thermal part could round it differently.
                    let baseline = row.outward.additional_specific_energy_j_kg.ok_or(invalid(
                        "active drifting source lost its original energy correction",
                    ))?;
                    let scale = relative_bulk_velocity_m_s
                        .into_iter()
                        .map(f64::abs)
                        .fold(fs_math::det::sqrt(-baseline), f64::max);
                    let scale = Dd::from_f64(scale);
                    let energy = relative_bulk_velocity_m_s.into_iter().fold(
                        Dd::from_f64(baseline) / scale / scale,
                        |sum, component| {
                            let component = Dd::from_f64(component) / scale;
                            sum + Dd::from_f64(0.5) * component * component
                        },
                    );
                    let correction = (energy * scale * scale).to_f64();
                    if !correction.is_finite() || (energy != Dd::ZERO && correction == 0.0) {
                        return Err(invalid(
                            "tangential effusive energy correction is unrepresentable",
                        ));
                    }
                    correction
                } else {
                    moments.crossing_energy_minus_enthalpy_j_kg
                };
                let h = row
                    .outward
                    .source_enthalpy_j_kg
                    .ok_or(invalid("active drifting source lost its original enthalpy"))?;
                let crossing = h + correction;
                if !crossing.is_finite() {
                    return Err(invalid(
                        "drifting molecular crossing energy is unrepresentable",
                    ));
                }
                let mass = incoming.incoming.mass_rate_kg_s;
                incoming.incoming.additional_specific_energy_j_kg = Some(correction);
                incoming.incoming.crossing_energy_j_kg = Some(crossing);
                incoming.incoming.energy_rate_w = signed_product(mass, crossing)?;
                incoming.incoming.additional_energy_rate_w = signed_product(mass, correction)?;
            }
            result.species.push(incoming);
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

#[cfg(test)]
#[path = "drift/tests.rs"]
mod tests;
