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

    /// The running sum and its compensation, as bits: the whole sum, for
    /// the codec of a state that holds one.
    #[must_use]
    pub fn to_words(&self) -> [u64; 2] {
        [self.sum.to_bits(), self.compensation.to_bits()]
    }

    /// The sum whose [`Self::to_words`] are `words`, to the bit.
    #[must_use]
    pub fn from_words(words: [u64; 2]) -> Self {
        Self {
            sum: f64::from_bits(words[0]),
            compensation: f64::from_bits(words[1]),
        }
    }

    fn push_words(&self, words: &mut Vec<u64>) {
        words.extend(self.to_words());
    }

    fn read(words: &mut Words<'_>) -> Option<Self> {
        Some(Self::from_words([words.u64()?, words.u64()?]))
    }
}

/// Plain 64-bit words read from the front, as a snapshot holds them.
pub(crate) struct Words<'a>(pub(crate) &'a [u64]);

impl Words<'_> {
    pub(crate) fn u64(&mut self) -> Option<u64> {
        let (&first, rest) = self.0.split_first()?;
        self.0 = rest;
        Some(first)
    }

    /// A count written as a word, refused above the words that are left,
    /// since every counted item takes at least one.
    pub(crate) fn count(&mut self) -> Option<usize> {
        let count = usize::try_from(self.u64()?).ok()?;
        (count <= self.0.len()).then_some(count)
    }

    fn flag(&mut self) -> Option<bool> {
        match self.u64()? {
            0 => Some(false),
            1 => Some(true),
            _ => None,
        }
    }

    fn sums(&mut self) -> Option<Vec<Sum>> {
        let count = self.count()?;
        (0..count).map(|_| Sum::read(self)).collect()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

fn push_sums(sums: &[Sum], words: &mut Vec<u64>) {
    words.push(sums.len() as u64);
    for sum in sums {
        sum.push_words(words);
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

/// Independently measured account for one carried stream coordinate.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CarriedAccount {
    /// Net transported amount entering the owner, in the coordinate's unit.
    pub port_in: Sum,
    /// Change in independently measured inventory.
    pub stored_change: Sum,
    /// Accepted transport-minus-inventory defects.
    pub balance_defect: Sum,
    /// Steps with both endpoint measurements.
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
    /// One account per carried coordinate; unavailable measurements remain
    /// `None`. Ordinals are transport coordinates, not inferred species names.
    pub carried: Vec<Option<CarriedAccount>>,
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
    /// Uncancelled transport per carried coordinate, in its declared unit.
    pub carried_defect: Vec<Sum>,
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
    /// Signed transported amounts per carried coordinate.
    pub carried_in: Vec<Sum>,
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
    /// Every number of the ledger as plain words, which is what a snapshot
    /// holds of it. Names come from the graph and are not written.
    pub(crate) fn push_words(&self, words: &mut Vec<u64>) {
        words.push(self.nodes.len() as u64);
        for node in &self.nodes {
            match &node.mass {
                None => words.push(0),
                Some(mass) => {
                    words.push(1);
                    for sum in [
                        &mass.port_mass_in_kg,
                        &mass.stored_change_kg,
                        &mass.balance_defect_kg,
                    ] {
                        sum.push_words(words);
                    }
                    words.push(mass.audited_steps);
                }
            }
            words.push(node.carried.len() as u64);
            for account in &node.carried {
                match account {
                    None => words.push(0),
                    Some(account) => {
                        words.push(1);
                        for sum in [
                            &account.port_in,
                            &account.stored_change,
                            &account.balance_defect,
                        ] {
                            sum.push_words(words);
                        }
                        words.push(account.audited_steps);
                    }
                }
            }
            for sum in [
                &node.port_energy_in_j,
                &node.stored_change_j,
                &node.dissipated_j,
                &node.numerical_loss_j,
                &node.entropy_generated_j_per_k,
                &node.balance_defect_j,
            ] {
                sum.push_words(words);
            }
        }
        words.push(self.connections.len() as u64);
        for connection in &self.connections {
            push_sums(&connection.carried_defect, words);
            for sum in [
                &connection.mass_defect_kg,
                &connection.energy_defect_j,
                &connection.entropy_generated_j_per_k,
            ] {
                sum.push_words(words);
            }
        }
        words.push(self.boundaries.len() as u64);
        for boundary in &self.boundaries {
            push_sums(&boundary.carried_in, words);
            for sum in [
                &boundary.mass_in_kg,
                &boundary.energy_in_j,
                &boundary.entropy_in_j_per_k,
            ] {
                sum.push_words(words);
            }
        }
        words.push(self.sinks.len() as u64);
        for sink in &self.sinks {
            sink.heat_j.push_words(words);
        }
        words.push(self.steps);
        self.duration_s.push_words(words);
    }

    /// The ledger that [`Self::push_words`] wrote as exactly `words`, named
    /// as `self`, which must have as many accounts of each kind. `None`
    /// refuses the words.
    pub(crate) fn read(&self, words: &[u64]) -> Option<Self> {
        let words = &mut Words(words);
        if words.count()? != self.nodes.len() {
            return None;
        }
        let mut nodes = Vec::with_capacity(self.nodes.len());
        for shape in &self.nodes {
            let mass = if words.flag()? {
                Some(MassAccount {
                    port_mass_in_kg: Sum::read(words)?,
                    stored_change_kg: Sum::read(words)?,
                    balance_defect_kg: Sum::read(words)?,
                    audited_steps: words.u64()?,
                })
            } else {
                None
            };
            let count = words.count()?;
            let carried = (0..count)
                .map(|_| {
                    if words.flag()? {
                        Some(Some(CarriedAccount {
                            port_in: Sum::read(words)?,
                            stored_change: Sum::read(words)?,
                            balance_defect: Sum::read(words)?,
                            audited_steps: words.u64()?,
                        }))
                    } else {
                        Some(None)
                    }
                })
                .collect::<Option<Vec<_>>>()?;
            nodes.push(NodeAccount {
                name: shape.name.clone(),
                has_streams: shape.has_streams,
                mass,
                carried,
                port_energy_in_j: Sum::read(words)?,
                stored_change_j: Sum::read(words)?,
                dissipated_j: Sum::read(words)?,
                numerical_loss_j: Sum::read(words)?,
                entropy_generated_j_per_k: Sum::read(words)?,
                balance_defect_j: Sum::read(words)?,
            });
        }
        if words.count()? != self.connections.len() {
            return None;
        }
        let mut connections = Vec::with_capacity(self.connections.len());
        for _ in &self.connections {
            connections.push(ConnectionAccount {
                carried_defect: words.sums()?,
                mass_defect_kg: Sum::read(words)?,
                energy_defect_j: Sum::read(words)?,
                entropy_generated_j_per_k: Sum::read(words)?,
            });
        }
        if words.count()? != self.boundaries.len() {
            return None;
        }
        let mut boundaries = Vec::with_capacity(self.boundaries.len());
        for shape in &self.boundaries {
            boundaries.push(BoundaryAccount {
                carried_in: words.sums()?,
                port: shape.port.clone(),
                mass_in_kg: Sum::read(words)?,
                energy_in_j: Sum::read(words)?,
                entropy_in_j_per_k: Sum::read(words)?,
            });
        }
        if words.count()? != self.sinks.len() {
            return None;
        }
        let mut sinks = Vec::with_capacity(self.sinks.len());
        for shape in &self.sinks {
            sinks.push(SinkAccount {
                name: shape.name.clone(),
                heat_j: Sum::read(words)?,
            });
        }
        let ledger = Self {
            nodes,
            connections,
            boundaries,
            sinks,
            steps: words.u64()?,
            duration_s: Sum::read(words)?,
        };
        words.is_empty().then_some(ledger)
    }

    /// Closure of one carried coordinate over the accepted window.
    ///
    /// Returns `None` before a step, for an absent coordinate, or if any owner
    /// carrying it omitted a measurement on any step. Owners must agree on
    /// the meaning and units of the ordinal; no chemistry is inferred here.
    #[must_use]
    pub fn carried_residual(&self, carried: usize) -> Option<f64> {
        if self.steps == 0 || !self.nodes.iter().any(|node| carried < node.carried.len()) {
            return None;
        }
        let mut total = Sum::default();
        for node in &self.nodes {
            if let Some(account) = node.carried.get(carried) {
                let account = account.as_ref()?;
                if account.audited_steps != self.steps {
                    return None;
                }
                total.add(account.stored_change.value());
            }
        }
        for boundary in &self.boundaries {
            if let Some(input) = boundary.carried_in.get(carried) {
                total.add(-input.value());
            }
        }
        for connection in &self.connections {
            if let Some(defect) = connection.carried_defect.get(carried) {
                total.add(-defect.value());
            }
        }
        Some(total.value())
    }

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
