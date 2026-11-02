//! Explicit local complete mixing of incoming liquid with an occupied cell.

use super::*;

/// Positive liquid inflow with its own carried momentum and caloric state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMassInflow {
    /// Positive mass, original liquid enthalpy datum and separately supplied
    /// heat. The source owner must debit the actual represented mass.
    pub exchange: FreeSurfaceMassExchange,
    /// Incoming physical velocity in the original Cartesian grid frame [m/s].
    pub velocity_m_s: [f64; 3],
}

/// Mass, momentum and enthalpy-plus-kinetic balance of complete local mixing.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceMassInflowReport {
    /// Original native inventory account. Its heat_j contains only the caller's
    /// external heat; the conversion from kinetic energy is separate below.
    pub exchange: FreeSurfaceMassExchangeReport,
    /// Kinetic energy carried by the represented incoming mass [J].
    pub carried_kinetic_energy_j: f64,
    /// Actual endpoint change of resolved liquid kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Nonnegative loss of relative kinetic energy under the explicitly
    /// selected perfectly inelastic mixing law, deposited once as enthalpy [J].
    pub mixing_heat_j: f64,
    /// Sum of absolute local enthalpy-plus-kinetic balance defects [J].
    pub absolute_total_energy_defect_j: f64,
}

#[derive(Default)]
pub(super) struct MixingBalance {
    pub(super) carried_kinetic_energy_j: f64,
    pub(super) kinetic_energy_change_j: f64,
    pub(super) mixing_heat_j: f64,
    pub(super) absolute_total_energy_defect_j: f64,
}

impl MixingBalance {
    pub(super) fn admitted(&self, tolerance_j: f64) -> bool {
        [
            self.carried_kinetic_energy_j,
            self.kinetic_energy_change_j,
            self.mixing_heat_j,
            self.absolute_total_energy_defect_j,
        ]
        .into_iter()
        .all(f64::is_finite)
            && self.mixing_heat_j >= 0.0
            && self.absolute_total_energy_defect_j <= tolerance_j
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Mix positive liquid inflow into existing interfaces without changing
    /// pressure density or erasing physical nonequilibrium stress.
    ///
    /// This selects a local perfectly inelastic mixing model: represented mass
    /// carries its supplied velocity, and the loss of relative kinetic energy
    /// heats the SAME original liquid caloric owner. No numerical remap loss is
    /// relabeled as heat. Both incoming and accepted speeds must fit the
    /// original scenario's lattice-speed bound. The independent mass, momentum,
    /// enthalpy and total-energy checks use actual represented transfers.
    ///
    /// This is an instantaneous open-system update at the accepted clock, not
    /// resolved jet-impact dynamics, droplet nucleation, pressure work, arbitrary
    /// evaporation recoil or a global VOF momentum-conservation claim. Native
    /// interface capacity, phase/density source admission and coupled
    /// surface/contact-history restrictions are unchanged. Failure publishes
    /// neither populations, fill, enthalpy nor refreshed material state.
    pub fn mix_interface_inflow(
        &mut self,
        inflows: &[FreeSurfaceMassInflow],
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceMassInflowReport, Error> {
        poll(gate)?;
        if inflows.len() > self.controls.max_cells {
            return Err(Error::Invalid("inflow exceeds the original cell budget"));
        }
        let scale = self.flow.scaling.velocity_si(1.0)?;
        let mut velocities = BTreeMap::new();
        let mut exchanges = Vec::with_capacity(inflows.len());
        for inflow in inflows {
            poll(gate)?;
            if !inflow.exchange.mass_kg.is_finite() || inflow.exchange.mass_kg <= 0.0 {
                return Err(Error::Invalid(
                    "complete mixing requires positive incoming mass",
                ));
            }
            // This boundary supplies liquid, so the incoming state itself
            // must fit the original chart and fixed storage-density owner.
            buoyancy::liquid_temperature(
                self.caloric,
                inflow.exchange.specific_enthalpy_j_kg,
                self.flow.scaling.material().density_kg_m3(),
            )?;
            admit_velocity(
                inflow.velocity_m_s.map(|u| u / scale),
                self.flow.scaling.flow().u_lattice,
            )?;
            let location = self
                .location_at(inflow.exchange.coordinates)
                .ok_or(Error::Invalid("inflow needs an existing liquid interface"))?;
            if velocities.insert(location, inflow.velocity_m_s).is_some() {
                return Err(Error::Invalid("inflow interfaces must be unique"));
            }
            exchanges.push(inflow.exchange);
        }
        let (exchange, balance) =
            self.exchange_interface_mass_model(&exchanges, Some(&velocities), controls, gate)?;
        Ok(FreeSurfaceMassInflowReport {
            exchange,
            carried_kinetic_energy_j: balance.carried_kinetic_energy_j,
            kinetic_energy_change_j: balance.kinetic_energy_change_j,
            mixing_heat_j: balance.mixing_heat_j,
            absolute_total_energy_defect_j: balance.absolute_total_energy_defect_j,
        })
    }
}

pub(super) fn admit_velocity(native: [f64; 3], maximum: f64) -> Result<(), Error> {
    let speed = native[0].hypot(native[1]).hypot(native[2]);
    if !speed.is_finite() || speed > maximum {
        return Err(Error::Invalid(
            "inflow exceeds the original scenario velocity bound",
        ));
    }
    Ok(())
}

/// Relative arithmetic preserves boosts and avoids subtracting large kinetic
/// energies. The native mass projection, not requested dm, supplies the masses.
pub(super) fn mixing_state(
    old_mass: f64,
    dm: f64,
    new_mass: f64,
    velocity: [f64; 3],
    incoming: [f64; 3],
) -> Option<([f64; 3], f64)> {
    if [old_mass, dm, new_mass]
        .iter()
        .any(|m| !m.is_finite() || *m <= 0.0)
        || dm > new_mass
        || velocity.iter().chain(&incoming).any(|u| !u.is_finite())
    {
        return None;
    }
    let difference: [f64; 3] = core::array::from_fn(|d| incoming[d] - velocity[d]);
    let fraction = dm / new_mass;
    let next = core::array::from_fn(|d| velocity[d] + fraction * difference[d]);
    let heat = 0.5 * dm * (old_mass / new_mass) * difference.iter().map(|v| v * v).sum::<f64>();
    (heat.is_finite() && next.iter().all(|v| v.is_finite())).then_some((next, heat))
}

pub(super) fn kinetic_change(dm: f64, new_mass: f64, velocity: [f64; 3], delta: [f64; 3]) -> f64 {
    0.5 * dm * velocity.iter().map(|v| v * v).sum::<f64>()
        + new_mass
            * (0..3)
                .map(|d| (velocity[d] + 0.5 * delta[d]) * delta[d])
                .sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_complete_mixing_closes_momentum_and_energy_under_density_scaling_and_boosts() {
        for scale in [0.4, 1.0, 2.5] {
            let (old, dm) = (0.4 * scale, 0.1 * scale);
            let new = old + dm;
            let u = [0.02, -0.01, 0.03];
            let v = [-0.04, 0.05, -0.02];
            let (baseline, heat) = mixing_state(old, dm, new, u, v).unwrap();
            for boost in [[0.0; 3], [0.03, -0.02, 0.01]] {
                let u = core::array::from_fn(|d| u[d] + boost[d]);
                let v = core::array::from_fn(|d| v[d] + boost[d]);
                let (next, boosted_heat) = mixing_state(old, dm, new, u, v).unwrap();
                for d in 0..3 {
                    assert!((new * next[d] - old * u[d] - dm * v[d]).abs() < 2e-17 * scale);
                    assert!((next[d] - baseline[d] - boost[d]).abs() < 2e-17);
                }
                assert!((boosted_heat - heat).abs() < 2e-18 * scale);
                let delta = core::array::from_fn(|d| next[d] - u[d]);
                let carried = 0.5 * dm * v.iter().map(|v| v * v).sum::<f64>();
                assert!((kinetic_change(dm, new, u, delta) + heat - carried).abs() < 2e-18 * scale);
                assert!(heat > 0.0);
            }
            assert_eq!(mixing_state(old, dm, new, u, u).unwrap(), (u, 0.0));
        }
        assert!(mixing_state(1.0, -0.1, 0.9, [0.0; 3], [0.0; 3]).is_none());
        assert!(mixing_state(1.0, 0.1, 1.1, [0.0; 3], [f64::NAN; 3]).is_none());
    }
}
