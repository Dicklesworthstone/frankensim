//! Prescribed-velocity liquid ejection with an explicit mechanical work debit.

use super::*;

/// Negative liquid mass exchange through an ideal powered outlet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMassOutflow {
    /// Strictly negative mass, carried liquid enthalpy and separately supplied
    /// heat. The external owner credits the actual represented outgoing mass.
    pub exchange: FreeSurfaceMassExchange,
    /// Outgoing physical velocity in the original Cartesian grid frame [m/s].
    pub velocity_m_s: [f64; 3],
}

/// Actual recoil and work of the explicitly selected ideal ejection law.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceMassOutflowReport {
    pub exchange: FreeSurfaceMassExchangeReport,
    /// Signed kinetic energy carried into the liquid [J]; negative on exit.
    pub carried_kinetic_energy_j: f64,
    /// Actual endpoint change of resolved liquid kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Nonnegative mechanical work supplied by the external actuator [J].
    /// The caller must debit this amount from its work source. It is never heat.
    pub mechanical_work_j: f64,
    /// Absolute local enthalpy-plus-kinetic defects after supplied work [J].
    pub absolute_total_energy_defect_j: f64,
}

/// Accepted native ejection held for joint publication with another owner.
/// Dropping this value leaves the original liquid unchanged. Its exclusive
/// borrow prevents any intervening change to the original state or clock.
#[must_use = "a prepared outflow changes the original liquid only when published"]
pub struct PreparedFreeSurfaceMassOutflow<'a, 's, 'm> {
    owner: &'a mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    candidate: MassExchangeCandidate,
    report: FreeSurfaceMassOutflowReport,
}

impl PreparedFreeSurfaceMassOutflow<'_, '_, '_> {
    /// Original native balance against the exclusively borrowed starting state.
    #[must_use]
    pub fn report(&self) -> &FreeSurfaceMassOutflowReport {
        &self.report
    }

    /// Publish the already accepted proposal without allocation or failure.
    /// A multi-owner caller must finish its joint checks and cancellation poll
    /// before publishing any participant; publication does not poll again.
    pub fn publish(self) -> FreeSurfaceMassOutflowReport {
        self.candidate.publish(self.owner);
        self.report
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Remove liquid at a prescribed velocity and apply its reciprocal recoil.
    ///
    /// For represented removal r, old mass m and remaining mass M, the accepted
    /// velocity is u' = u - (r/M)(v_out-u). An ideal actuator supplies
    /// W = 0.5*r*(m/M)*|v_out-u|², with zero net impulse on liquid plus stream.
    /// `maximum_mechanical_work_j` explicitly authorizes a finite nonnegative
    /// work debit across the entire batch; insufficient work refuses the whole
    /// update. The actual debit is returned for the external source owner.
    /// No part of W becomes thermal energy, and no negative mixing heat exists.
    /// The supplied stream enthalpy and separate heat retain their original
    /// caloric meaning. Both outlet and remaining speeds must fit the original
    /// material scenario. Native pressure density and physical nonequilibrium
    /// stress are retained through the original Guo force-frame owner.
    ///
    /// This is an instantaneous, local ideal actuator law. It does not infer
    /// evaporation recoil, nozzle pressure work, or resolved jet dynamics.
    /// Positive remaining interface mass, original liquid phase/density and
    /// source domains, and existing coupled surface/contact restrictions still
    /// apply. Every refusal preserves the complete accepted liquid and clock.
    pub fn eject_interface_outflow(
        &mut self,
        outflows: &[FreeSurfaceMassOutflow],
        maximum_mechanical_work_j: f64,
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceMassOutflowReport, Error> {
        self.prepare_interface_outflow(outflows, maximum_mechanical_work_j, controls, gate)
            .map(|prepared| prepared.publish())
    }

    /// Prepare the same native ejection transaction for joint acceptance.
    /// Exactly one lattice/enthalpy proposal is built; the complete owner and
    /// its unchanged source/history bindings are not cloned. The original
    /// recoil, work, source and cancellation checks run before returning.
    pub fn prepare_interface_outflow<'a>(
        &'a mut self,
        outflows: &[FreeSurfaceMassOutflow],
        maximum_mechanical_work_j: f64,
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<PreparedFreeSurfaceMassOutflow<'a, 's, 'm>, Error> {
        poll(gate)?;
        if outflows.len() > self.controls.max_cells
            || !maximum_mechanical_work_j.is_finite()
            || maximum_mechanical_work_j < 0.0
        {
            return Err(Error::Invalid(
                "invalid powered-outflow work or cell budget",
            ));
        }
        let scale = self.flow.scaling.velocity_si(1.0)?;
        let mut velocities = BTreeMap::new();
        let mut exchanges = Vec::with_capacity(outflows.len());
        for outflow in outflows {
            poll(gate)?;
            if !outflow.exchange.mass_kg.is_finite() || outflow.exchange.mass_kg >= 0.0 {
                return Err(Error::Invalid(
                    "powered ejection requires negative liquid mass",
                ));
            }
            buoyancy::liquid_temperature(
                self.caloric,
                outflow.exchange.specific_enthalpy_j_kg,
                self.flow.scaling.material().density_kg_m3(),
            )?;
            inflow::admit_velocity(
                outflow.velocity_m_s.map(|u| u / scale),
                self.flow.scaling.flow().u_lattice,
            )?;
            let location = self
                .location_at(outflow.exchange.coordinates)
                .ok_or(Error::Invalid("outflow needs an existing liquid interface"))?;
            if velocities.insert(location, outflow.velocity_m_s).is_some() {
                return Err(Error::Invalid("outflow interfaces must be unique"));
            }
            exchanges.push(outflow.exchange);
        }
        let (candidate, exchange, balance) = self.prepare_interface_mass_model(
            &exchanges,
            Some(TransferMomentum::Ejection {
                velocities: &velocities,
                maximum_work_j: maximum_mechanical_work_j,
            }),
            controls,
            gate,
        )?;
        let report = FreeSurfaceMassOutflowReport {
            exchange,
            carried_kinetic_energy_j: balance.carried_kinetic_energy_j,
            kinetic_energy_change_j: balance.kinetic_energy_change_j,
            mechanical_work_j: balance.mechanical_work_j,
            absolute_total_energy_defect_j: balance.absolute_total_energy_defect_j,
        };
        Ok(PreparedFreeSurfaceMassOutflow {
            owner: self,
            candidate,
            report,
        })
    }
}

/// Actual represented signed mass supplies the recoil, never a requested rate.
/// Relative arithmetic retains boosts and avoids subtracting large energies.
pub(super) fn ejection_state(
    old_mass: f64,
    dm: f64,
    new_mass: f64,
    velocity: [f64; 3],
    outgoing: [f64; 3],
) -> Option<([f64; 3], f64)> {
    if !old_mass.is_finite()
        || !dm.is_finite()
        || !new_mass.is_finite()
        || old_mass <= 0.0
        || new_mass <= 0.0
        || dm >= 0.0
        || -dm >= old_mass
        || velocity.iter().chain(&outgoing).any(|u| !u.is_finite())
    {
        return None;
    }
    let relative: [f64; 3] = core::array::from_fn(|d| outgoing[d] - velocity[d]);
    let fraction = dm / new_mass;
    let delta: [f64; 3] = relative.map(|v| fraction * v);
    let next = core::array::from_fn(|d| velocity[d] + delta[d]);
    // v_out-u' = relative-delta. Each product is nonnegative. This equals
    // r*m/M*|relative|²/2 without cancellation or an unnecessary large m/M.
    let work = -0.5
        * dm
        * (0..3)
            .map(|d| relative[d] * (relative[d] - delta[d]))
            .sum::<f64>();
    (work.is_finite() && work >= 0.0 && next.iter().all(|v| v.is_finite())).then_some((next, work))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::freesurface3::{ContactModel3, FreeSurfaceInitialCell3};

    #[test]
    fn g1_powered_ejection_closes_mass_momentum_and_work_under_density_and_boosts() {
        for density in [0.4, 1.0, 2.5] {
            let m = 0.4 * density;
            let dm = -0.1 * density;
            let remaining = m + dm;
            let u = [0.02, -0.01, 0.03];
            let outgoing = [-0.04, 0.05, -0.02];
            let (baseline, work) = ejection_state(m, dm, remaining, u, outgoing).unwrap();
            for boost in [[0.0; 3], [0.03, -0.02, 0.01]] {
                let u = core::array::from_fn(|d| u[d] + boost[d]);
                let v = core::array::from_fn(|d| outgoing[d] + boost[d]);
                let (next, boosted_work) = ejection_state(m, dm, remaining, u, v).unwrap();
                for d in 0..3 {
                    assert!((remaining * next[d] - m * u[d] - dm * v[d]).abs() < 3e-17 * density);
                    assert!((next[d] - baseline[d] - boost[d]).abs() < 3e-17);
                }
                assert!((boosted_work - work).abs() < 3e-18 * density);
                let delta = core::array::from_fn(|d| next[d] - u[d]);
                let carried = 0.5 * dm * v.iter().map(|v| v * v).sum::<f64>();
                let kinetic = inflow::kinetic_change(dm, remaining, u, delta);
                assert!((kinetic - carried - work).abs() < 3e-18 * density);
                assert!(work > 0.0);
            }
            assert_eq!(ejection_state(m, dm, remaining, u, u).unwrap(), (u, 0.0));
        }
        assert!(ejection_state(1.0, 0.1, 1.1, [0.0; 3], [0.0; 3]).is_none());
        assert!(ejection_state(1.0, -1.0, 0.0, [0.0; 3], [0.0; 3]).is_none());
        assert!(ejection_state(1.0, -0.1, 0.9, [0.0; 3], [f64::NAN; 3]).is_none());
    }

    #[test]
    fn g1_native_ejection_retains_density_and_force_frame_with_actual_recoil_and_work() {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 71,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            for rho in [0.4, 0.8, 1.3, 2.5] {
                let mut flow = FreeSurface3::from_cells(
                    [8; 3],
                    0.8,
                    [0.0003, -0.0002, 0.0001],
                    0.0,
                    ContactModel3::Neutral,
                    &[],
                    &[FreeSurfaceInitialCell3 {
                        coordinates: [3; 3],
                        fill_fraction: 0.4,
                        density: rho,
                        velocity: [0.02, -0.01, 0.03],
                    }],
                    &cx,
                )
                .unwrap();
                let location = FreeSurface3::transport_location([3; 3]).unwrap();
                let before = flow.sample_at(3, 3, 3).unwrap();
                let (m, density) = flow.transport_interface_mass(location).unwrap();
                let remaining = m - 0.1 * density;
                let dm = remaining - m;
                let outgoing = [-0.04, 0.05, -0.02];
                let (next, work) =
                    ejection_state(m, dm, remaining, before.velocity, outgoing).unwrap();
                flow.transport_set_interface_mass(location, remaining, &gate)
                    .unwrap();
                flow.transport_set_interface_velocity(location, next, &gate)
                    .unwrap();
                let after = flow.sample_at(3, 3, 3).unwrap();
                assert_eq!(before.density.to_bits(), after.density.to_bits());
                assert_eq!(flow.ledger_mass(), remaining);
                assert!((after.fill_fraction - 0.3).abs() < 2e-14);
                for d in 0..3 {
                    assert!(
                        (remaining * after.velocity[d] - m * before.velocity[d] - dm * outgoing[d])
                            .abs()
                            < 2e-14 * rho
                    );
                }
                let delta = core::array::from_fn(|d| after.velocity[d] - before.velocity[d]);
                let kinetic = inflow::kinetic_change(dm, remaining, before.velocity, delta);
                let carried = 0.5 * dm * outgoing.iter().map(|u| u * u).sum::<f64>();
                assert!((kinetic - carried - work).abs() < 2e-15 * rho);
                let almost_empty = remaining * 1e-5;
                let excessive = ejection_state(
                    remaining,
                    almost_empty - remaining,
                    almost_empty,
                    after.velocity,
                    outgoing,
                )
                .unwrap()
                .0;
                assert!(
                    flow.transport_set_interface_velocity(location, excessive, &gate)
                        .is_err()
                );
                assert_eq!(flow.sample_at(3, 3, 3).unwrap(), after);
            }
        });
    }
}
