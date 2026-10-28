//! Native data adapter for the original molecular-relaxation owner.
//!
//! `iso9613-terrestrial-air` selects the existing, composition-limited law.
//! `declared N` is followed by N species records `ID RESERVOIRS`, each with
//! that many reservoir/rate pairs. Reservoirs are `vibration MODES THETA_K G …`
//! or `rotation DOF`; rates are `power P_TAU_REF T_REF EXPONENT T_MIN T_MAX`
//! or `partners N ID LOG_P_TAU SLOPE T_MIN T_MAX …`.
//!
//! All quantities are the original law's SI parameters. Declared coefficients
//! are authored constitutive data, not a measurement or a source certificate.
//! No missing species, collision partner or zero-reservoir declaration is
//! inferred. The owning law checks physical admissibility and temperature
//! domains. Counts are bounded before allocation; refusal leaves input intact.

use super::{
    InternalReservoir, MAX_COLLISION_PARTNERS, MAX_MODES_PER_PROCESS, MAX_RELAXATION_PROCESSES,
    MAX_RELAXING_SPECIES, PartnerRate, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use crate::MaterialError;
use crate::gas_species::input::GasRecord;
use fs_thermochem::SpeciesId;
use std::str::SplitAsciiWhitespace;

fn invalid(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

struct Fields<'a>(SplitAsciiWhitespace<'a>);
impl<'a> Fields<'a> {
    fn word(&mut self) -> Result<&'a str, MaterialError> {
        self.0
            .next()
            .ok_or_else(|| invalid("missing molecular-relaxation field"))
    }

    fn count(&mut self, maximum: usize, zero: bool) -> Result<usize, MaterialError> {
        let value = self
            .word()?
            .parse::<usize>()
            .map_err(|_| invalid("invalid molecular-relaxation count"))?;
        if value > maximum || (!zero && value == 0) {
            return Err(invalid(
                "molecular-relaxation count exceeds its explicit finite allowance",
            ));
        }
        Ok(value)
    }

    fn scalar(&mut self) -> Result<f64, MaterialError> {
        self.word()?
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| invalid("molecular-relaxation coefficients must be finite SI values"))
    }

    fn species(&mut self) -> Result<SpeciesId, MaterialError> {
        SpeciesId::new(self.word()?).map_err(|error| invalid(error.to_string()))
    }

    fn process(&mut self) -> Result<RelaxationProcess, MaterialError> {
        let reservoir = match self.word()? {
            "vibration" => {
                let count = self.count(MAX_MODES_PER_PROCESS, false)?;
                let mut modes = Vec::with_capacity(count);
                for _ in 0..count {
                    modes.push(VibrationalMode {
                        theta_k: self.scalar()?,
                        degeneracy: self.count(usize::from(u8::MAX), false)? as u8,
                    });
                }
                InternalReservoir::Vibrational(modes)
            }
            "rotation" => InternalReservoir::Rotational {
                degrees_of_freedom: self.count(3, false)? as u8,
            },
            _ => {
                return Err(invalid(
                    "molecular reservoir requires vibration or rotation",
                ));
            }
        };
        let rate = match self.word()? {
            "power" => RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: self.scalar()?,
                t_ref_k: self.scalar()?,
                exponent: self.scalar()?,
                t_min_k: self.scalar()?,
                t_max_k: self.scalar()?,
            },
            "partners" => {
                let count = self.count(MAX_COLLISION_PARTNERS, false)?;
                let mut partners = Vec::with_capacity(count);
                for _ in 0..count {
                    partners.push(PartnerRate {
                        partner: self.species()?,
                        ln_p_tau_pa_s: self.scalar()?,
                        landau_teller_slope_k13: self.scalar()?,
                        t_min_k: self.scalar()?,
                        t_max_k: self.scalar()?,
                    });
                }
                RelaxationRateLaw::CollisionPartners(partners)
            }
            _ => return Err(invalid("molecular rate requires power or partners")),
        };
        Ok(RelaxationProcess { reservoir, rate })
    }
}

/// Read one complete constitutive declaration, leaving any caller-owned suffix.
/// This parses data only; evaluate it against the actual material gas separately.
pub fn read_relaxation(
    fields: &mut SplitAsciiWhitespace<'_>,
) -> Result<RelaxationAbsorption, MaterialError> {
    let mut row = Fields(fields.clone());
    let law = match row.word()? {
        "iso9613-terrestrial-air" => RelaxationAbsorption::iso9613_terrestrial_air()?,
        "declared" => {
            let count = row.count(MAX_RELAXING_SPECIES, false)?;
            let mut species = Vec::with_capacity(count);
            let mut remaining = MAX_RELAXATION_PROCESSES;
            for _ in 0..count {
                let id = row.species()?;
                let processes = row.count(remaining, true)?;
                remaining -= processes;
                let mut reservoirs = Vec::with_capacity(processes);
                for _ in 0..processes {
                    reservoirs.push(row.process()?);
                }
                species.push(SpeciesRelaxation {
                    species: id,
                    processes: reservoirs,
                });
            }
            RelaxationAbsorption::new(species)?
        }
        _ => {
            return Err(invalid(
                "molecular law requires iso9613-terrestrial-air or declared",
            ));
        }
    };
    *fields = row.0;
    Ok(law)
}

/// Preserve the numerical gas and its explicitly declared species authority.
/// A legacy pseudocomponent record cannot identify molecular reservoirs.
pub fn medium_from_record(record: &GasRecord) -> Result<RelaxingMedium, MaterialError> {
    match record {
        GasRecord::Legacy(_) => Err(invalid(
            "molecular relaxation requires explicit mixture or thermochemical species; legacy air has no species record",
        )),
        GasRecord::Mixture(gas, composition) => {
            let components = composition
                .iter()
                .map(|(species, fraction)| {
                    SpeciesId::new(species.symbol())
                        .map(|id| (id, *fraction))
                        .map_err(|error| invalid(error.to_string()))
                })
                .collect::<Result<Vec<_>, _>>()?;
            RelaxingMedium::new(*gas, components)
        }
        GasRecord::Thermochemical(gas) => RelaxingMedium::from_thermochemical(gas),
    }
}
