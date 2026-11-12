//! Window ledger: what crossed every port, connection and boundary.
//!
//! The driver fills it from the receipts and the wire values of each accepted
//! step. Components write no accounting code. Sums are compensated so that a
//! small term is not lost beside a large one over many steps.

/// A running sum with Neumaier compensation.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Sum {
    sum: f64,
    compensation: f64,
}

impl Sum {
    /// Add one term.
    pub fn add(&mut self, term: f64) {
        let total = self.sum + term;
        if self.sum.abs() >= term.abs() {
            self.compensation += (self.sum - total) + term;
        } else {
            self.compensation += (term - total) + self.sum;
        }
        self.sum = total;
    }

    /// The compensated total.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.sum + self.compensation
    }
}

/// Independently measured mass accounts of one stream component.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MassAccount {
    /// Net mass entering the component through stream ports [kg].
    pub port_mass_in_kg: Sum,
    /// Change in independently measured stored mass [kg].
    pub stored_change_kg: Sum,
    /// Sum of accepted port-transfer minus stored-change defects [kg].
    pub balance_defect_kg: Sum,
    /// Accepted steps for which the component supplied both measurements.
    pub audited_steps: u64,
}

/// Accounts of one component over the window.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct NodeAccount {
    /// Component name.
    pub name: String,
    /// Whether this component has a material stream port.
    pub has_streams: bool,
    /// Measured stream mass, or unavailable. A partial account has fewer
    /// audited steps than the ledger and never proves window closure.
    pub mass: Option<MassAccount>,
    /// Net energy that entered through its ports [J].
    pub port_energy_in_j: Sum,
    /// Change of stored energy, from receipts [J].
    pub stored_change_j: Sum,
    /// Physical dissipation [J].
    pub dissipated_j: Sum,
    /// Numerical loss [J]; never added to dissipation.
    pub numerical_loss_j: Sum,
    /// Entropy generated inside the component [J/K].
    pub entropy_generated_j_per_k: Sum,
    /// Sum over steps of the first-law defect the audit accepted [J].
    pub balance_defect_j: Sum,
}

/// Accounts of one connection over the window.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConnectionAccount {
    /// Net stream mass entering the connection's members [kg].
    pub mass_defect_kg: Sum,
    /// Sum of the energies that entered the members, plus heat sent to the
    /// connection's sink [J]. Zero for an exact connection; it is the
    /// interface residual of an iterative solve otherwise. Heat routes that
    /// share a receiver are summed in the row of the first of them, and the
    /// rows of the others stay zero.
    pub energy_defect_j: Sum,
    /// Entropy generated in the connection [J/K].
    pub entropy_generated_j_per_k: Sum,
}

/// Accounts of one boundary over the window.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BoundaryAccount {
    /// The driven port, as `node.port`.
    pub port: String,
    /// Mass entering the system through this stream boundary [kg].
    pub mass_in_kg: Sum,
    /// Energy that entered the system through it [J].
    pub energy_in_j: Sum,
    /// Entropy that entered the system through it [J/K].
    pub entropy_in_j_per_k: Sum,
}

/// Accounts of one named heat sink outside the system.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SinkAccount {
    /// Sink name.
    pub name: String,
    /// Heat delivered to it [J].
    pub heat_j: Sum,
}

/// Energy, entropy and declared stream-mass accounts of accepted steps.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Ledger {
    /// One account per component, in index order.
    pub nodes: Vec<NodeAccount>,
    /// One account per connection, in declaration order.
    pub connections: Vec<ConnectionAccount>,
    /// One account per boundary, in declaration order.
    pub boundaries: Vec<BoundaryAccount>,
    /// One account per declared heat sink.
    pub sinks: Vec<SinkAccount>,
    /// Accepted steps in the window.
    pub steps: u64,
    /// Simulated time covered by the window [s].
    pub duration_s: Sum,
}

impl Ledger {
    /// Stream-mass residual [kg]: measured storage change minus boundary
    /// inflow and connection defect, or minus the accepted component defects.
    ///
    /// Returns None before any accepted step, without stream components, or
    /// when any stream component omitted an inventory on any accepted step.
    /// Missing measurements are never interpreted as zero. Species,
    /// momentum, charge and mass of non-stream components are outside this
    /// account.
    #[must_use]
    pub fn mass_residual_kg(&self) -> Option<f64> {
        if self.steps == 0 || !self.nodes.iter().any(|node| node.has_streams) {
            return None;
        }
        let mut total = Sum::default();
        for node in &self.nodes {
            if node.has_streams {
                let account = node.mass.as_ref()?;
                if account.audited_steps != self.steps {
                    return None;
                }
                total.add(account.stored_change_kg.value());
            }
        }
        for boundary in &self.boundaries {
            total.add(-boundary.mass_in_kg.value());
        }
        for connection in &self.connections {
            total.add(-connection.mass_defect_kg.value());
        }
        Some(total.value())
    }

    /// Energy that entered the system through its boundaries [J].
    #[must_use]
    pub fn boundary_energy_in_j(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.boundaries {
            total.add(account.energy_in_j.value());
        }
        total.value()
    }

    /// Change of the energy stored in all components [J].
    #[must_use]
    pub fn stored_change_j(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.nodes {
            total.add(account.stored_change_j.value());
        }
        total.value()
    }

    /// Heat delivered to sinks outside the system [J].
    #[must_use]
    pub fn sink_heat_j(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.sinks {
            total.add(account.heat_j.value());
        }
        total.value()
    }

    /// Numerical loss of all components [J].
    #[must_use]
    pub fn numerical_loss_j(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.nodes {
            total.add(account.numerical_loss_j.value());
        }
        total.value()
    }

    /// Interface residual energy of all connections [J].
    #[must_use]
    pub fn connection_defect_j(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.connections {
            total.add(account.energy_defect_j.value());
        }
        total.value()
    }

    /// Entropy generated in components and connections [J/K].
    #[must_use]
    pub fn entropy_generated_j_per_k(&self) -> f64 {
        let mut total = Sum::default();
        for account in &self.nodes {
            total.add(account.entropy_generated_j_per_k.value());
        }
        for account in &self.connections {
            total.add(account.entropy_generated_j_per_k.value());
        }
        total.value()
    }

    /// First-law residual of the window [J]:
    /// `stored change + numerical loss + sink heat - boundary energy in
    ///  - connection defect`.
    ///
    /// It equals minus the sum of the per-step balance defects the component
    /// audits accepted (each defect is port energy in, minus stored change,
    /// minus numerical loss), so it is roundoff for exact components and it
    /// shows a leak that a component hid behind a wide tolerance.
    #[must_use]
    pub fn energy_residual_j(&self) -> f64 {
        let mut total = Sum::default();
        total.add(self.stored_change_j());
        total.add(self.numerical_loss_j());
        total.add(self.sink_heat_j());
        total.add(-self.boundary_energy_in_j());
        total.add(-self.connection_defect_j());
        total.value()
    }
}

#[cfg(test)]
mod tests {
    use super::Sum;

    #[test]
    fn a_small_term_survives_beside_large_ones() {
        let mut sum = Sum::default();
        sum.add(1.0e16);
        sum.add(1.0);
        sum.add(-1.0e16);
        assert_eq!(sum.value().to_bits(), 1.0_f64.to_bits());
    }
}
