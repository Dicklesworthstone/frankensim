//! The SAME molecular streams for two finite gases at one shared aperture.
//!
//! Each original drifting Maxwellian is sampled once in the common aperture
//! frame. Positive L-to-R and R-to-L masses retain separate conditional vector
//! impulses and complete crossing energies. The two caloric coordinates are
//! built from one signed source-relative transfer, so formation references
//! cancel before either finite gas is reconstructed. This sampler neither
//! spends inventories nor supplies a particle-flight or wall-force model.

use super::{
    ClosedGasParcel, DriftingEffusiveFlux, DriftingEffusiveSpeciesFlux, EffusiveDriftInterface,
    GasEnergyError, SpeciesId, invalid, poll, positive_product_ratio, signed_product,
};
use fs_math::dd::Dd;

/// One pair of original gross streams, integrated over the declared interval.
#[derive(Clone, Debug, PartialEq)]
pub struct DriftingEffusiveSpeciesExchange {
    pub species: SpeciesId,
    /// Actual gross masses [L-to-R, R-to-L], never canceled before withdrawal.
    pub mass_kg: [f64; 2],
    /// Signed transfer, preserving small original drift-rate increments even
    /// when the displayed gross masses round to the same floating-point value.
    pub net_mass_left_to_right_kg: f64,
    /// Each gross beam's vector impulse in the COMMON aperture frame.
    /// The reverse beam normally has a negative normal component.
    pub momentum_relative_kg_m_s: [[f64; 3]; 2],
    /// Complete gross energies in that same frame and original source datum.
    /// These are reports; do not subtract them to reconstruct small heat flow.
    pub energy_j: [f64; 2],
    /// Signed energy INTO [left, right], beyond each signed species mass
    /// change times its original h_i(T). The composing moving-gas owner must
    /// subtract its actual bulk kinetic change before its ONE caloric solve.
    pub relative_energy_into_gas_j: [f64; 2],
}

/// Both source samples and their common conservative transfer coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct DriftingEffusiveExchange {
    pub duration_s: f64,
    pub left_to_right: DriftingEffusiveFlux,
    pub right_to_left: DriftingEffusiveFlux,
    /// Canonical union of original active source species, including selected
    /// zero-transmission species; empty only for a geometrically closed port.
    pub species: Vec<DriftingEffusiveSpeciesExchange>,
}

impl DriftingEffusiveExchange {
    /// Retain original source, drift and reciprocal transmission support on
    /// both accepted endpoints without sampling or spending any stream again.
    /// Sampled and endpoint admission do not enclose a continuous trajectory.
    pub fn admit_endpoints(
        &self,
        gases: [&ClosedGasParcel; 2],
        relative_bulk_velocity_m_s: [[f64; 3]; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        poll(&mut cancelled)?;
        if self.species.is_empty() {
            return Ok(());
        }
        let drift = &self.left_to_right.interface;
        let original = drift.momentum_interface().interface();
        for (gas, velocity) in gases.into_iter().zip(relative_bulk_velocity_m_s) {
            poll(&mut cancelled)?;
            drift.admit_relative_velocity(velocity)?;
            original.admit_parcel(gas, &mut cancelled)?;
        }
        for row in &self.species {
            poll(&mut cancelled)?;
            let index = original
                .transmissions()
                .binary_search_by(|t| t.species().cmp(&row.species))
                .map_err(|_| invalid("reciprocal endpoint lost its original transmission"))?;
            for gas in gases {
                let state = gas.state().state();
                original.transmissions()[index].admit(state.temperature, state.pressure)?;
            }
        }
        poll(&mut cancelled)
    }
}

fn represented(value: Dd) -> Result<f64, GasEnergyError> {
    let result = value.to_f64();
    if !result.is_finite() || (value != Dd::ZERO && result == 0.0) {
        Err(invalid("paired molecular coordinate is unrepresentable"))
    } else {
        Ok(result)
    }
}

fn integrate(rate: f64, duration: f64) -> Result<f64, GasEnergyError> {
    if rate == 0.0 {
        Ok(0.0)
    } else {
        positive_product_ratio(&[rate, duration], &[])
    }
}

fn find<'a>(
    flux: &'a DriftingEffusiveFlux,
    species: &SpeciesId,
) -> Option<&'a DriftingEffusiveSpeciesFlux> {
    flux.species
        .binary_search_by(|row| row.species.cmp(species))
        .ok()
        .map(|index| &flux.species[index])
}

fn signed_mass(
    beams: [Option<&DriftingEffusiveSpeciesFlux>; 2],
    mass: [f64; 2],
    duration: f64,
) -> Result<f64, GasEnergyError> {
    // Near rest, retain each signed original perturbation separately. In a
    // far tail, use the positive direct gross arrival instead: baseline plus
    // its nearly opposite change would erase a real molecular beam.
    if beams.iter().all(|beam| {
        beam.is_none_or(|b| {
            b.incoming_mass_rate_change_from_zero_drift_kg_s.abs()
                <= 0.5 * b.baseline_incoming.mass_rate_kg_s
        })
    }) {
        let mut net = Dd::ZERO;
        for (side, beam) in beams.into_iter().enumerate() {
            if let Some(beam) = beam {
                let base = integrate(beam.baseline_incoming.mass_rate_kg_s, duration)?;
                let change = signed_product(
                    duration,
                    beam.incoming_mass_rate_change_from_zero_drift_kg_s,
                )?;
                let value = Dd::from_f64(base) + Dd::from_f64(change);
                net = if side == 0 { net + value } else { net - value };
            }
        }
        represented(net)
    } else {
        represented(Dd::from_f64(mass[0]) - Dd::from_f64(mass[1]))
    }
}

/// Sample reciprocal molecular exchange at one explicitly shared interface.
///
/// `normal_left_to_right` is oriented from left gas to right gas. The velocities
/// are each gas's bulk velocity MINUS the same aperture-frame velocity. The
/// original incidence sampler points its normal out of the receiving side, so
/// L-to-R uses -normal and R-to-L uses +normal. Both streams use the same original
/// transmission law and retain their distinct species, temperature and drift.
///
/// The finite owner must admit both gross original inventories and publish both
/// endpoints atomically. A common inertial frame's velocity V converts the net
/// impulse to the laboratory frame as J_relative + V * net_mass. Apply these
/// same particles with opposite signs; independent donor/receiver resampling
/// would create a different transfer. Fixed sampled rates remain first order.
/// Zero area is inactive without inspecting gas, normal or drift support.
#[allow(clippy::too_many_arguments)]
pub fn sample(
    left: &ClosedGasParcel,
    right: &ClosedGasParcel,
    interface: &EffusiveDriftInterface,
    area_m2: f64,
    relative_bulk_velocity_m_s: [[f64; 3]; 2],
    normal_left_to_right: [f64; 3],
    duration_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DriftingEffusiveExchange, GasEnergyError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(invalid(
            "paired molecular exchange needs positive finite duration",
        ));
    }
    let left_to_right = left.drifting_effusive_flux(
        interface,
        area_m2,
        relative_bulk_velocity_m_s[0],
        normal_left_to_right.map(|value| -value),
        &mut cancelled,
    )?;
    let right_to_left = right.drifting_effusive_flux(
        interface,
        area_m2,
        relative_bulk_velocity_m_s[1],
        normal_left_to_right,
        &mut cancelled,
    )?;
    let mut result = DriftingEffusiveExchange {
        duration_s,
        left_to_right,
        right_to_left,
        species: Vec::new(),
    };
    if area_m2 == 0.0 {
        poll(&mut cancelled)?;
        return Ok(result);
    }
    let original = interface.momentum_interface().interface();
    let temperatures = [left, right].map(|gas| gas.state().state().temperature);
    result.species.reserve(original.transmissions().len());
    for transmission in original.transmissions() {
        poll(&mut cancelled)?;
        let id = transmission.species();
        let beams = [
            find(&result.left_to_right, id),
            find(&result.right_to_left, id),
        ];
        if beams.iter().all(Option::is_none) {
            continue;
        }
        // A selected species qualifies its reciprocal passage at BOTH source
        // states even if absent from one side or assigned zero transmission.
        for gas in [left, right] {
            let state = gas.state().state();
            transmission.admit(state.temperature, state.pressure)?;
        }
        let mut mass = [0.0; 2];
        let mut impulse = [[0.0; 3]; 2];
        let mut energy = [0.0; 2];
        let mut correction = [None; 2];
        for (side, beam) in beams.into_iter().enumerate() {
            if let Some(beam) = beam {
                mass[side] = integrate(beam.incoming.mass_rate_kg_s, duration_s)?;
                if mass[side] > 0.0 {
                    impulse[side] = beam.impulse_for_mass_kg(mass[side])?;
                    let crossing = beam.incoming.crossing_energy_j_kg.ok_or(invalid(
                        "active reciprocal beam lost original crossing energy",
                    ))?;
                    energy[side] = signed_product(mass[side], crossing)?;
                    correction[side] = Some(beam.incoming.additional_specific_energy_j_kg.ok_or(
                        invalid("active reciprocal beam lost original enthalpy crosswalk"),
                    )?);
                }
            }
        }
        let net = signed_mass(beams, mass, duration_s)?;
        let mut relative_energy = [0.0; 2];
        if mass != [0.0; 2] {
            let delta_h = original
                .gas_catalog()
                .species_specific_enthalpy_change_j_kg(id, temperatures[0], temperatures[1])?;
            // One common left-reference energy coordinate:
            // C = net*a_L + m_RL*(a_L-a_R-Delta h).
            // If the left beam is absent, choosing a_L=a_R gives the SAME
            // identity because net=-m_RL. No fictitious left beam is created.
            let a_right = correction[1].unwrap_or(0.0);
            let a_left = correction[0].unwrap_or(a_right);
            let common = Dd::from_f64(net) * Dd::from_f64(a_left)
                + Dd::from_f64(mass[1])
                    * (Dd::from_f64(a_left) - Dd::from_f64(a_right) - Dd::from_f64(delta_h));
            relative_energy[0] = represented(Dd::ZERO - common)?;
            relative_energy[1] = represented(common - Dd::from_f64(net) * Dd::from_f64(delta_h))?;
        }
        result.species.push(DriftingEffusiveSpeciesExchange {
            species: id.clone(),
            mass_kg: mass,
            net_mass_left_to_right_kg: net,
            momentum_relative_kg_m_s: impulse,
            energy_j: energy,
            relative_energy_into_gas_j: relative_energy,
        });
    }
    poll(&mut cancelled)?;
    Ok(result)
}
