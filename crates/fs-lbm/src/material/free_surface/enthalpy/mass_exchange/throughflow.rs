//! Simultaneous inlet mixing and powered outlet recoil on one native state.

use super::*;

mod balance;
pub(super) use balance::ThroughflowCell;

/// Gross represented stream masses for one original native interface.
///
/// Each direction's co-located streams share its represented/requested ratio.
/// The external source debits incoming mass and credits outgoing mass separately,
/// including equal-throughput transfers whose native net mass change is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMassThroughflowSample {
    pub coordinates: [usize; 3],
    pub requested_incoming_mass_kg: f64,
    pub requested_outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub outgoing_mass_kg: f64,
}

/// One accepted simultaneous open-system mass, momentum and energy transaction.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceMassThroughflowReport {
    /// Net original-liquid inventory, source and actual native endpoint account.
    pub exchange: FreeSurfaceMassExchangeReport,
    /// Actual net incoming minus outgoing kinetic energy, including dispersion [J].
    pub carried_kinetic_energy_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Nonnegative inlet relative-motion loss deposited in the SAME caloric owner [J].
    pub mixing_heat_j: f64,
    /// Nonnegative outlet work that the external actuator must debit [J].
    pub mechanical_work_j: f64,
    pub absolute_total_energy_defect_j: f64,
    /// Stable native-coordinate order, one entry per occupied interface.
    pub samples: Vec<FreeSurfaceMassThroughflowSample>,
}

/// Accepted throughflow held for joint publication with its source owners.
#[must_use = "a prepared throughflow changes the original liquid only when published"]
pub struct PreparedFreeSurfaceMassThroughflow<'a, 's, 'm> {
    owner: &'a mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    candidate: MassExchangeCandidate,
    report: FreeSurfaceMassThroughflowReport,
}

impl PreparedFreeSurfaceMassThroughflow<'_, '_, '_> {
    #[must_use]
    pub fn report(&self) -> &FreeSurfaceMassThroughflowReport {
        &self.report
    }

    /// Publish without allocation or failure. All participants and the final
    /// cancellation poll must be accepted before any owner's publication.
    pub fn publish(self) -> FreeSurfaceMassThroughflowReport {
        self.candidate.publish(self.owner);
        self.report
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Evolve simultaneous liquid streams without an inlet-first/outlet-first split.
    ///
    /// Streams hold their authored velocity, enthalpy and relative mass-flow
    /// proportions over one integrated local transfer. For incoming mass a,
    /// outgoing b and progress s in [0,1], m=m0+s*(a-b) and
    /// p=m0*u0+s*(a*v_in-b*v_out). The existing infinitesimal complete-mixing
    /// loss a*|v_in-u(s)|²/2 and powered-outlet work b*|v_out-u(s)|²/2 are
    /// integrated analytically along that SAME path. Stream velocity dispersion
    /// contributes to each direction independently. Thus zero net mass can
    /// still transfer momentum, caloric energy, mixing heat and actuator work.
    ///
    /// Native rounding changes only the larger directional mass batch; the
    /// smaller batch is retained. Each returned directional mass ratio owns
    /// its external stream debits/credits. Gross rounding and net closure both
    /// consume the original mass-defect budget. No division by the net flux or
    /// artificial effective stream velocity/enthalpy is used. Single-direction
    /// cells retain the original mixing/ejection arithmetic.
    ///
    /// This is a held-stream local mixing/actuation law, not a nozzle, resolved
    /// jet, pressure-work or finite donor dynamics model. The clock does not
    /// advance. One original lattice/enthalpy candidate and one material/force
    /// refresh retain the original capacity, phase, stress and source domains.
    /// The final liquid and all external owners must be accepted together.
    ///
    /// # Errors
    /// Refuses invalid streams, exhausted original budgets, unrepresentable
    /// state, unavailable actuator work or cancellation without publishing.
    pub fn exchange_interface_throughflow(
        &mut self,
        inflows: &[FreeSurfaceMassInflow],
        outflows: &[FreeSurfaceMassOutflow],
        maximum_mechanical_work_j: f64,
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceMassThroughflowReport, Error> {
        self.prepare_interface_throughflow(
            inflows,
            outflows,
            maximum_mechanical_work_j,
            controls,
            gate,
        )
        .map(PreparedFreeSurfaceMassThroughflow::publish)
    }

    /// Prepare the complete transaction under an exclusive original-state borrow.
    /// Dropping the proposal preserves the original liquid and accepted clock.
    ///
    /// # Errors
    /// Uses the same admissions and atomic refusal as `exchange_interface_throughflow`.
    #[allow(clippy::too_many_lines)]
    pub fn prepare_interface_throughflow<'a>(
        &'a mut self,
        inflows: &[FreeSurfaceMassInflow],
        outflows: &[FreeSurfaceMassOutflow],
        maximum_mechanical_work_j: f64,
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<PreparedFreeSurfaceMassThroughflow<'a, 's, 'm>, Error> {
        poll(gate)?;
        if inflows
            .len()
            .checked_add(outflows.len())
            .is_none_or(|n| n > self.controls.max_cells)
            || !maximum_mechanical_work_j.is_finite()
            || maximum_mechanical_work_j < 0.0
        {
            return Err(Error::Invalid("invalid throughflow work or stream budget"));
        }
        let scale = self.flow.scaling.velocity_si(1.0)?;
        let mut cells = BTreeMap::<_, ThroughflowCell>::new();
        let streams = inflows
            .iter()
            .map(|s| (0, s.exchange, s.velocity_m_s))
            .chain(outflows.iter().map(|s| (1, s.exchange, s.velocity_m_s)));
        for (side, exchange, velocity) in streams {
            poll(gate)?;
            let mass = if side == 0 {
                exchange.mass_kg
            } else {
                -exchange.mass_kg
            };
            if !mass.is_finite() || mass <= 0.0 {
                return Err(Error::Invalid(
                    "throughflow requires correctly signed nonzero streams",
                ));
            }
            buoyancy::liquid_temperature(
                self.caloric,
                exchange.specific_enthalpy_j_kg,
                self.flow.scaling.material().density_kg_m3(),
            )?;
            inflow::admit_velocity(
                velocity.map(|u| u / scale),
                self.flow.scaling.flow().u_lattice,
            )?;
            let location = self
                .location_at(exchange.coordinates)
                .ok_or(Error::Invalid(
                    "throughflow needs an existing liquid interface",
                ))?;
            let stream = inflow::IncomingMixture::new(
                mass,
                velocity,
                exchange.specific_enthalpy_j_kg,
                exchange.heat_j,
            )
            .ok_or(Error::Invalid("throughflow stream is unrepresentable"))?;
            cells
                .entry(location)
                .or_default()
                .add(side, stream)
                .ok_or(Error::Invalid(
                    "combined throughflow streams are unrepresentable",
                ))?;
        }
        let mut exchanges = Vec::with_capacity(cells.len());
        let mut samples = Vec::with_capacity(cells.len());
        for (&location, cell) in &mut cells {
            poll(gate)?;
            let coordinates =
                FreeSurface3::transport_coordinates(location).expect("admitted native interface");
            let requested = cell.requested();
            let net = requested[0] - requested[1];
            let projection = self.project_interface_mass_exchange(coordinates, net)?;
            cell.represent(projection.mass_change_kg)
                .ok_or(Error::Invalid(
                    "represented throughflow masses are unrepresentable",
                ))?;
            exchanges.push(FreeSurfaceMassExchange {
                coordinates,
                mass_kg: net,
                // The shared owner uses the separate directional enthalpies;
                // this datum is only the existing liquid's native reference.
                specific_enthalpy_j_kg: self.h[&location],
                heat_j: cell.heat(),
            });
            samples.push(FreeSurfaceMassThroughflowSample {
                coordinates,
                requested_incoming_mass_kg: requested[0],
                requested_outgoing_mass_kg: requested[1],
                incoming_mass_kg: cell.mass_kg[0],
                outgoing_mass_kg: cell.mass_kg[1],
            });
        }
        let (candidate, exchange, balance) = self.prepare_interface_mass_model(
            &exchanges,
            Some(TransferMomentum::Throughflow {
                cells: &cells,
                maximum_work_j: maximum_mechanical_work_j,
            }),
            controls,
            gate,
        )?;
        Ok(PreparedFreeSurfaceMassThroughflow {
            owner: self,
            candidate,
            report: FreeSurfaceMassThroughflowReport {
                exchange,
                carried_kinetic_energy_j: balance.carried_kinetic_energy_j,
                kinetic_energy_change_j: balance.kinetic_energy_change_j,
                mixing_heat_j: balance.mixing_heat_j,
                mechanical_work_j: balance.mechanical_work_j,
                absolute_total_energy_defect_j: balance.absolute_total_energy_defect_j,
                samples,
            },
        })
    }
}
