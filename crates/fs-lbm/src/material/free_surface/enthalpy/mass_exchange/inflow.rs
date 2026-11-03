//! Explicit local complete mixing of incoming liquid with an occupied cell.

use super::*;

mod mixture;
pub(super) use mixture::IncomingMixture;

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
    /// Signed kinetic energy carried by represented streams [J], including
    /// incoming velocity dispersion and any comoving withdrawals.
    pub carried_kinetic_energy_j: f64,
    /// Actual endpoint change of resolved liquid kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Nonnegative loss of relative kinetic energy under the explicitly
    /// selected perfectly inelastic mixing law, deposited once as enthalpy [J].
    pub mixing_heat_j: f64,
    /// Sum of absolute local enthalpy-plus-kinetic balance defects [J].
    pub absolute_total_energy_defect_j: f64,
}

/// Accepted native inflow held for joint publication with another owner.
/// Dropping this value leaves the original liquid unchanged. Its exclusive
/// borrow prevents any intervening change to the original state or clock.
#[must_use = "a prepared inflow changes the original liquid only when published"]
pub struct PreparedFreeSurfaceMassInflow<'a, 's, 'm> {
    owner: &'a mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    candidate: MassExchangeCandidate,
    report: FreeSurfaceMassInflowReport,
}

impl PreparedFreeSurfaceMassInflow<'_, '_, '_> {
    /// Original native balance against the exclusively borrowed starting state.
    #[must_use]
    pub fn report(&self) -> &FreeSurfaceMassInflowReport {
        &self.report
    }

    /// Publish the already accepted proposal without allocation or failure.
    /// A multi-owner caller must finish its joint checks and cancellation poll
    /// before publishing any participant; publication does not poll again.
    pub fn publish(self) -> FreeSurfaceMassInflowReport {
        self.candidate.publish(self.owner);
        self.report
    }
}

#[derive(Default)]
pub(super) struct KineticBalance {
    pub(super) carried_kinetic_energy_j: f64,
    pub(super) kinetic_energy_change_j: f64,
    pub(super) mixing_heat_j: f64,
    pub(super) mechanical_work_j: f64,
    pub(super) absolute_total_energy_defect_j: f64,
}

impl KineticBalance {
    pub(super) fn admitted(&self, tolerance_j: f64) -> bool {
        [
            self.carried_kinetic_energy_j,
            self.kinetic_energy_change_j,
            self.mixing_heat_j,
            self.mechanical_work_j,
            self.absolute_total_energy_defect_j,
        ]
        .into_iter()
        .all(f64::is_finite)
            && self.mixing_heat_j >= 0.0
            && self.mechanical_work_j >= 0.0
            && self.absolute_total_energy_defect_j <= tolerance_j
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
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
    /// Multiple streams may enter the same interface. They are combined once
    /// using mass-weighted velocity and enthalpy, while their relative kinetic
    /// energy remains explicit and becomes mixing heat exactly once. Every
    /// stream must independently fit this liquid's original caloric source and
    /// velocity domain; different material/composition laws are not inferred.
    /// The report has one sample per interface. Co-located streams share that
    /// sample's represented/requested mass ratio, which the source owners use
    /// to debit their individual mass, momentum and carried energy. Reduction
    /// follows authored order; permutations agree to roundoff, not bit for bit.
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
        self.prepare_interface_inflow(inflows, controls, gate)
            .map(|prepared| prepared.publish())
    }

    /// Prepare the same native mixing transaction for joint acceptance.
    /// Exactly one lattice/enthalpy proposal is built; the complete owner and
    /// its unchanged source/history bindings are not cloned. All native source,
    /// momentum, heat, capacity and cancellation checks run before returning.
    pub fn prepare_interface_inflow<'a>(
        &'a mut self,
        inflows: &[FreeSurfaceMassInflow],
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<PreparedFreeSurfaceMassInflow<'a, 's, 'm>, Error> {
        self.prepare_interface_inflow_with_exchange(inflows, &[], controls, gate)
    }

    /// Mix positive incoming streams together with comoving withdrawals and
    /// zero-mass heat patches in ONE original native transaction. Nonpositive
    /// exchanges retain their actual old liquid velocity; no powered ejection
    /// or artificial negative mixing heat is inferred. Coordinates must be
    /// disjoint between the two batches. All source owners publish together.
    pub fn mix_interface_inflow_with_exchange(
        &mut self,
        inflows: &[FreeSurfaceMassInflow],
        comoving: &[FreeSurfaceMassExchange],
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceMassInflowReport, Error> {
        self.prepare_interface_inflow_with_exchange(inflows, comoving, controls, gate)
            .map(PreparedFreeSurfaceMassInflow::publish)
    }

    /// Prepare the mixed signed batch without changing any accepted owner.
    /// The original mass, caloric, population and force-refresh path runs once.
    pub fn prepare_interface_inflow_with_exchange<'a>(
        &'a mut self,
        inflows: &[FreeSurfaceMassInflow],
        comoving: &[FreeSurfaceMassExchange],
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<PreparedFreeSurfaceMassInflow<'a, 's, 'm>, Error> {
        poll(gate)?;
        if inflows
            .len()
            .checked_add(comoving.len())
            .is_none_or(|n| n > self.controls.max_cells)
        {
            return Err(Error::Invalid("inflow exceeds the original cell budget"));
        }
        let scale = self.flow.scaling.velocity_si(1.0)?;
        let mut mixtures = BTreeMap::<_, IncomingMixture>::new();
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
            let stream = IncomingMixture::new(
                inflow.exchange.mass_kg,
                inflow.velocity_m_s,
                inflow.exchange.specific_enthalpy_j_kg,
                inflow.exchange.heat_j,
            )
            .ok_or(Error::Invalid("incoming stream is unrepresentable"))?;
            if let Some(mixture) = mixtures.get_mut(&location) {
                mixture.combine(stream).ok_or(Error::Invalid(
                    "combined incoming streams are unrepresentable",
                ))?;
            } else {
                mixtures.insert(location, stream);
            }
        }
        let mut exchanges: Vec<_> = mixtures
            .iter()
            .map(|(&location, mixture)| FreeSurfaceMassExchange {
                coordinates: FreeSurface3::transport_coordinates(location)
                    .expect("admitted native interface"),
                mass_kg: mixture.mass_kg,
                specific_enthalpy_j_kg: mixture.specific_enthalpy_j_kg,
                heat_j: mixture.heat_j,
            })
            .collect();
        for exchange in comoving {
            poll(gate)?;
            if !exchange.mass_kg.is_finite() || exchange.mass_kg > 0.0 {
                return Err(Error::Invalid(
                    "comoving batch requires nonpositive transfers",
                ));
            }
            let location = self
                .location_at(exchange.coordinates)
                .ok_or(Error::Invalid(
                    "comoving batch needs an existing liquid interface",
                ))?;
            if mixtures.contains_key(&location) {
                return Err(Error::Invalid(
                    "inflow and comoving batches must use distinct cells",
                ));
            }
            exchanges.push(*exchange);
        }
        let (candidate, exchange, balance) = self.prepare_interface_mass_model(
            &exchanges,
            Some(TransferMomentum::Mixing(&mixtures)),
            controls,
            gate,
        )?;
        let report = FreeSurfaceMassInflowReport {
            exchange,
            carried_kinetic_energy_j: balance.carried_kinetic_energy_j,
            kinetic_energy_change_j: balance.kinetic_energy_change_j,
            mixing_heat_j: balance.mixing_heat_j,
            absolute_total_energy_defect_j: balance.absolute_total_energy_defect_j,
        };
        Ok(PreparedFreeSurfaceMassInflow {
            owner: self,
            candidate,
            report,
        })
    }
}

pub(super) fn admit_velocity(native: [f64; 3], maximum: f64) -> Result<(), Error> {
    let speed = native[0].hypot(native[1]).hypot(native[2]);
    if !speed.is_finite() || speed > maximum {
        return Err(Error::Invalid(
            "interface stream exceeds the original scenario velocity bound",
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
    fn signed_batch_keeps_withdrawals_and_thermal_only_cells_comoving() {
        let input = MassLocation3::Cell(0, 0);
        let output = MassLocation3::Cell(0, 1);
        let thermal = MassLocation3::Cell(0, 2);
        let mixtures = BTreeMap::from([(
            input,
            IncomingMixture::new(0.2, [0.1, -0.2, 0.0], 100.0, 0.0).unwrap(),
        )]);
        let model = TransferMomentum::Mixing(&mixtures);
        let velocity = [0.02, 0.03, -0.01];
        assert_eq!(
            model.state(&output, 1.0, -0.2, 0.8, velocity),
            Some((velocity, 0.0, 0.0))
        );
        assert_eq!(
            model.state(&thermal, 1.0, 0.0, 1.0, velocity),
            Some((velocity, 0.0, 0.0))
        );
        assert!(model.state(&output, 1.0, 0.2, 1.2, velocity).is_none());
        let mixed = model.state(&input, 1.0, 0.2, 1.2, velocity).unwrap();
        assert!(mixed.1 > 0.0);
        for d in 0..3 {
            assert!((1.2 * mixed.0[d] - velocity[d] - 0.2 * [0.1, -0.2, 0.0][d]).abs() < 1e-17);
        }
        let outgoing = kinetic_change(-0.2, 0.8, velocity, [0.0; 3]);
        assert_eq!(outgoing, -0.1 * velocity.iter().map(|u| u * u).sum::<f64>());
        assert_eq!(model.velocity_at(&output), None);
        assert_eq!(model.kinetic_variance_at(&thermal), 0.0);
    }

    #[test]
    fn g1_shared_cell_state_closes_all_stream_and_resident_momenta_and_kinetics() {
        let location = MassLocation3::Cell(0, 0);
        for scale in [1e-10, 0.4, 1.0, 1e10] {
            for boost in [[0.0; 3], [0.4, -0.3, 0.2]] {
                let sources = [
                    (0.2 * scale, [0.08, -0.04, 0.02]),
                    (0.3 * scale, [-0.06, 0.03, -0.01]),
                    (0.1 * scale, [0.02, 0.05, -0.04]),
                ]
                .map(|(m, v)| {
                    IncomingMixture::new(
                        m,
                        core::array::from_fn(|d| v[d] + boost[d]),
                        200_000.0,
                        0.0,
                    )
                    .unwrap()
                });
                let mut mixed = sources[0];
                mixed.combine(sources[1]).unwrap();
                mixed.combine(sources[2]).unwrap();
                let mixtures = BTreeMap::from([(location, mixed)]);
                let model = TransferMomentum::Mixing(&mixtures);
                let old = 0.4 * scale;
                let velocity = core::array::from_fn(|d| [0.02, -0.01, 0.03][d] + boost[d]);
                for fraction in [0.37, 1.0, 1.0 + 2.0 * f64::EPSILON] {
                    let dm = mixed.mass_kg * fraction;
                    let new = old + dm;
                    let (next, heat, work) =
                        model.state(&location, old, dm, new, velocity).unwrap();
                    for d in 0..3 {
                        let supplied: f64 = sources
                            .iter()
                            .map(|s| fraction * s.mass_kg * s.velocity_m_s[d])
                            .sum();
                        assert!(
                            (new * next[d] - old * velocity[d] - supplied).abs() < 2e-15 * scale
                        );
                    }
                    let supplied: f64 = sources
                        .iter()
                        .map(|s| {
                            0.5 * fraction
                                * s.mass_kg
                                * s.velocity_m_s.iter().map(|u| u * u).sum::<f64>()
                        })
                        .sum();
                    let delta = core::array::from_fn(|d| next[d] - velocity[d]);
                    let changed = kinetic_change(dm, new, velocity, delta);
                    assert!((changed + heat - supplied).abs() < 2e-15 * scale);
                    assert!(heat > dm * mixed.specific_kinetic_variance_j_kg);
                    assert_eq!(work, 0.0);
                }
            }
        }
    }

    #[test]
    fn g3_single_stream_and_ejection_keep_original_state_arithmetic() {
        let location = MassLocation3::Cell(0, 0);
        let velocity = [0.02, -0.01, 0.03];
        let carried = [-0.04, 0.05, -0.02];
        let mixtures = BTreeMap::from([(
            location,
            IncomingMixture::new(0.1, carried, 200_000.0, 0.0).unwrap(),
        )]);
        let original = mixing_state(0.4, 0.1, 0.5, velocity, carried).unwrap();
        assert_eq!(
            TransferMomentum::Mixing(&mixtures).state(&location, 0.4, 0.1, 0.5, velocity),
            Some((original.0, original.1, 0.0))
        );
        let outgoing = TransferMomentum::Ejection {
            mixtures: &mixtures,
            maximum_work_j: 1.0,
        };
        let original = outflow::ejection_state(0.4, -0.1, 0.3, velocity, carried).unwrap();
        assert_eq!(outgoing.kinetic_variance_at(&location), 0.0);
        assert_eq!(
            outgoing.state(&location, 0.4, -0.1, 0.3, velocity),
            Some((original.0, 0.0, original.1))
        );
    }

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
