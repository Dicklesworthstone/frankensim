//! Electrode-free conducting components use a gauge, not a fictitious ground.

use super::*;

impl ElectricSetup {
    /// Induced circulating currents, magnetic reaction and Joule heat in one or
    /// more electrically isolated bodies. `electrodes` must be empty; the solid
    /// preset table remains ignored. Both solids and fluids use the same law.
    ///
    /// The lowest-index conductor in each face-connected component fixes the
    /// arbitrary potential offset to zero. It retains its complete resistance,
    /// motional voltage, contact heat and magnetic force. No physical electrode,
    /// leakage current or external electrical power is introduced. Insulators
    /// retain NaN potential; a single disconnected conductor has zero current.
    ///
    /// All face currents are re-accumulated at EVERY cell, including the gauge,
    /// and must close within `current_tolerance_a` [A]. Power still obeys the
    /// separate tolerance in `motion`. Neither limit is an error estimator.
    /// Frozen-state, imposed-field, low-induction Ohmic physics only; no charging,
    /// capacitance, time-dependent field induction or mechanical time integration.
    /// Inputs remain unchanged on all failures; external query receipts are owned
    /// by the caller. Missing material data must be an error, never insulation.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_isolated_cell_motion(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        motion: ElectricMotionInput<'_>,
        current_tolerance_a: f64,
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<MovingConductorSolution, ChtError> {
        self.solve_isolated_cell_motion_topology(
            domain, temperature, contacts, &[], motion, current_tolerance_a, resolve, gate,
        )
    }

    /// Eddy/braking currents in electrically isolated, segmented conductors.
    ///
    /// Explicit open faces interrupt current loops without removing either
    /// bulk cell, changing its resistivity or changing mechanical/thermal
    /// contact. Numerical gauges follow the ACTUAL conducting graph after
    /// openings: each separate component receives exactly one deterministic
    /// potential reference and all of its cells still receive current, Joule
    /// heat and force. An insulating slit can therefore extinguish a loop;
    /// separated laminations can retain independent circulating currents.
    ///
    /// The original current and transpose-force assembly consume the same
    /// face map as gauge selection. Closed faces and finite resistive contacts
    /// retain their original laws. All electrode-free, continuity, power,
    /// material, cancellation and low-induction restrictions of
    /// [`Self::solve_isolated_cell_motion`] apply. No infinite-resistance
    /// surrogate, capacitive gap current or topological history is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_isolated_cell_motion_topology(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        open_faces: &[[usize; 2]],
        motion: ElectricMotionInput<'_>,
        current_tolerance_a: f64,
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<MovingConductorSolution, ChtError> {
        poll(gate)?;
        if !self.electrodes.is_empty() || !current_tolerance_a.is_finite() || current_tolerance_a < 0.0 {
            return Err(invalid("isolated motion requires no electrodes and a finite nonnegative current tolerance"));
        }
        let mut transfer = MotionCoupling::new(domain, motion, gate)?;
        let (nodes, rho) = super::super::contact::resolve_admitted_cells(
            domain, temperature, vec![Node::Off; domain.cell_count()], resolve, gate,
        )?;
        let field = super::super::contact::ContactField::new(domain, &nodes, contacts, gate)?
            .with_open_faces(domain, open_faces, gate)?;
        transfer.isolation = Some(Isolation::new(domain, &nodes, &field, current_tolerance_a, gate)?);
        let electric = solve_electric_field_with_contact_field(
            domain, self, nodes, &rho, field, Some(&mut transfer), gate,
        )?;
        transfer.finish(domain, electric, gate)
    }
}

pub(super) struct Isolation {
    pub(super) potential_gauges: Vec<bool>,
    current_balance_a: Vec<f64>,
    tolerance_a: f64,
}

impl Isolation {
    fn new(domain: &VoxelDomain, nodes: &[Node], contacts: &super::super::contact::ContactField,
        tolerance_a: f64, gate: &CancelGate) -> Result<Self, ChtError> {
        let n = domain.cell_count();
        let mut visited = vec![false; n];
        let mut potential_gauges = vec![false; n];
        let mut stack = Vec::new();
        let mut visits = 0;
        for c in 0..n {
            poll_every(gate, c)?;
            if visited[c] || nodes[c] != Node::Conductor { continue; }
            potential_gauges[c] = true;
            visited[c] = true;
            stack.push(c);
            while let Some(cell) = stack.pop() {
                poll_every(gate, visits)?;
                visits += 1;
                for face in 0..6 {
                    if let Some(next) = domain.neighbor(cell, face)
                        && !visited[next] && nodes[next] == Node::Conductor
                        && !contacts.is_open(cell, next)
                    {
                        visited[next] = true;
                        stack.push(next);
                    }
                }
            }
        }
        Ok(Self { potential_gauges, current_balance_a: vec![0.0; n], tolerance_a })
    }

    pub(super) fn record(&mut self, c: usize, n: usize, current: f64) -> Result<(), ChtError> {
        self.current_balance_a[c] += current;
        self.current_balance_a[n] -= current;
        if !self.current_balance_a[c].is_finite() || !self.current_balance_a[n].is_finite() {
            return Err(invalid("isolated current balance is unrepresentable"));
        }
        Ok(())
    }

    pub(super) fn check(&self, gate: &CancelGate) -> Result<(), ChtError> {
        for (cell, &residual) in self.current_balance_a.iter().enumerate() {
            poll_every(gate, cell)?;
            if !residual.is_finite() || residual.abs() > self.tolerance_a {
                return Err(invalid(format!("isolated cell {cell} has current residual {residual} A exceeding {} A", self.tolerance_a)));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gauge_cell_current_balance_is_checked_not_just_unknown_rows() {
        let domain=VoxelDomain::new(2,2,1,1.0).unwrap();let gate=CancelGate::new();
        let contacts=super::super::super::contact::ContactField::new(&domain,&[Node::Conductor;4],&[],&gate).unwrap();
        let mut isolation=Isolation::new(&domain,&[Node::Conductor;4],&contacts,1e-12,&gate).unwrap();
        assert_eq!(isolation.potential_gauges,vec![true,false,false,false]);
        isolation.current_balance_a[0]=1.0;
        assert!(isolation.check(&gate).is_err());
    }
}
