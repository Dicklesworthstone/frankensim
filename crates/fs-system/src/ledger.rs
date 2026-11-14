//! Window ledger: what crossed every port, connection and boundary.
//!
//! The driver fills it from the receipts and the wire values of each accepted
//! step. Components write no accounting code. Every account is kept in
//! double-double ([`Account`]), so that a small term is not lost beside a
//! large one over many steps, and the totals that combine accounts are taken
//! in double-double too and rounded once: a residual far below the rounding
//! of one large account survives. Beside each signed balance defect an
//! account sums its absolute value per step, so that defects of opposite sign
//! cannot cancel unseen. The clock and the window's duration are times, not
//! accounts, and stay [`Sum`]s.

use fs_math::dd::Dd;

/// A running sum with Neumaier compensation: the system's clock, source
/// clocks and the window's duration.
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

    /// Distance to an absolute time, retaining the clock's low bits.
    pub(crate) fn interval_to(self, time_s: f64) -> f64 {
        let mut from_event = self;
        from_event.add(-time_s);
        -from_event.value()
    }

    /// Preserve a strict cap without leaving only clock rounding before
    /// an event. Used by both source clocks and scheduled boundaries.
    pub(crate) fn limit_before_event(self, proposed_s: f64, event_s: f64) -> f64 {
        let mut at_event = self;
        at_event.add(event_s);
        let mut after_cap = self;
        after_cap.add(proposed_s);
        let leftover_s = at_event.value() - after_cap.value();
        let roundoff_s = 8.0 * f64::EPSILON * self.value().abs().max(at_event.value().abs());
        if event_s > proposed_s && leftover_s > 0.0 && leftover_s <= roundoff_s {
            proposed_s.min(0.5 * event_s)
        } else {
            proposed_s
        }
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

/// One ledger account: the running total of one quantity over accepted
/// steps, kept in double-double (fs-math's [`Dd`], normalized `hi + lo`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Account(Dd);

impl Account {
    /// Add one term (double-double addition, relative error at most 2^-104).
    pub fn add(&mut self, term: f64) {
        self.0 = self.0 + Dd::from_f64(term);
    }

    /// The total, rounded once to the nearest `f64`.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.0.to_f64()
    }

    /// The total in double-double, for a sum across accounts that is rounded
    /// once at its end.
    #[must_use]
    pub fn dd(&self) -> Dd {
        self.0
    }

    fn push_words(&self, words: &mut Vec<u64>) {
        words.extend([self.0.hi.to_bits(), self.0.lo.to_bits()]);
    }

    /// The account written as `hi` and `lo` bits; a pair that is not
    /// normalized is no account the ledger holds and is refused.
    fn read(words: &mut Words<'_>) -> Option<Self> {
        let (hi, lo) = (f64::from_bits(words.u64()?), f64::from_bits(words.u64()?));
        let normal = Dd::from_pair(hi, lo);
        (normal.hi.to_bits() == hi.to_bits() && normal.lo.to_bits() == lo.to_bits())
            .then_some(Self(normal))
    }
}

/// The sum of `accounts` in double-double.
fn total<'a>(accounts: impl IntoIterator<Item = &'a Account>) -> Dd {
    accounts
        .into_iter()
        .fold(Dd::ZERO, |sum, account| sum + account.0)
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

    fn accounts(&mut self) -> Option<Vec<Account>> {
        let count = self.count()?;
        (0..count).map(|_| Account::read(self)).collect()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

fn push_accounts(accounts: &[Account], words: &mut Vec<u64>) {
    words.push(accounts.len() as u64);
    for account in accounts {
        account.push_words(words);
    }
}

/// Independently measured mass accounts of one stream component.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MassAccount {
    /// Net mass entering the component through stream ports [kg].
    pub port_mass_in_kg: Account,
    /// Change in independently measured stored mass [kg].
    pub stored_change_kg: Account,
    /// Sum of accepted port-transfer minus stored-change defects [kg].
    pub balance_defect_kg: Account,
    /// Sum over steps of the absolute value of each step's defect \[kg\].
    pub absolute_defect_kg: Account,
    /// Accepted steps for which the component supplied both measurements.
    pub audited_steps: u64,
}

/// Independently measured account for one carried stream coordinate.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CarriedAccount {
    /// Net transported amount entering the owner, in the coordinate's unit.
    pub port_in: Account,
    /// Change in independently measured inventory.
    pub stored_change: Account,
    /// Accepted transport-minus-inventory defects.
    pub balance_defect: Account,
    /// Sum over steps of the absolute value of each step's defect.
    pub absolute_defect: Account,
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
    pub port_energy_in_j: Account,
    /// Change of stored energy, from receipts [J].
    pub stored_change_j: Account,
    /// Physical dissipation [J].
    pub dissipated_j: Account,
    /// Numerical loss [J]; never added to dissipation.
    pub numerical_loss_j: Account,
    /// Entropy generated inside the component [J/K].
    pub entropy_generated_j_per_k: Account,
    /// Sum over steps of the first-law defect the audit accepted [J].
    pub balance_defect_j: Account,
    /// Sum over steps of the absolute value of each step's accepted
    /// first-law defect \[J\]: defects of opposite sign do not cancel in it.
    pub absolute_defect_j: Account,
}

/// Accounts of one connection over the window.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConnectionAccount {
    /// Uncancelled transport per carried coordinate, in its declared unit.
    pub carried_defect: Vec<Account>,
    /// Net stream mass entering the connection's members [kg].
    pub mass_defect_kg: Account,
    /// Sum of the energies that entered the members, plus heat sent to the
    /// connection's sink [J]. Zero for an exact connection; it is the
    /// interface residual of an iterative solve otherwise. Heat routes that
    /// share a receiver are summed in the row of the first of them, and the
    /// rows of the others stay zero.
    pub energy_defect_j: Account,
    /// Entropy generated in the connection [J/K].
    pub entropy_generated_j_per_k: Account,
    /// Sum over steps of the absolute value of each step's energy defect \[J\].
    pub absolute_energy_defect_j: Account,
    /// Sum over steps of the absolute value of each step's mass defect \[kg\].
    pub absolute_mass_defect_kg: Account,
}

/// Accounts of one boundary over the window.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BoundaryAccount {
    /// Signed transported amounts per carried coordinate.
    pub carried_in: Vec<Account>,
    /// The driven port, as `node.port`.
    pub port: String,
    /// Mass entering the system through this stream boundary [kg].
    pub mass_in_kg: Account,
    /// Energy that entered the system through it [J].
    pub energy_in_j: Account,
    /// Entropy that entered the system through it [J/K].
    pub entropy_in_j_per_k: Account,
}

/// Accounts of one named heat sink outside the system.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SinkAccount {
    /// Sink name.
    pub name: String,
    /// Heat delivered to it [J].
    pub heat_j: Account,
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
                    for account in [
                        &mass.port_mass_in_kg,
                        &mass.stored_change_kg,
                        &mass.balance_defect_kg,
                        &mass.absolute_defect_kg,
                    ] {
                        account.push_words(words);
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
                        for account in [
                            &account.port_in,
                            &account.stored_change,
                            &account.balance_defect,
                            &account.absolute_defect,
                        ] {
                            account.push_words(words);
                        }
                        words.push(account.audited_steps);
                    }
                }
            }
            for account in [
                &node.port_energy_in_j,
                &node.stored_change_j,
                &node.dissipated_j,
                &node.numerical_loss_j,
                &node.entropy_generated_j_per_k,
                &node.balance_defect_j,
                &node.absolute_defect_j,
            ] {
                account.push_words(words);
            }
        }
        words.push(self.connections.len() as u64);
        for connection in &self.connections {
            push_accounts(&connection.carried_defect, words);
            for account in [
                &connection.mass_defect_kg,
                &connection.energy_defect_j,
                &connection.entropy_generated_j_per_k,
                &connection.absolute_energy_defect_j,
                &connection.absolute_mass_defect_kg,
            ] {
                account.push_words(words);
            }
        }
        words.push(self.boundaries.len() as u64);
        for boundary in &self.boundaries {
            push_accounts(&boundary.carried_in, words);
            for account in [
                &boundary.mass_in_kg,
                &boundary.energy_in_j,
                &boundary.entropy_in_j_per_k,
            ] {
                account.push_words(words);
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
                    port_mass_in_kg: Account::read(words)?,
                    stored_change_kg: Account::read(words)?,
                    balance_defect_kg: Account::read(words)?,
                    absolute_defect_kg: Account::read(words)?,
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
                            port_in: Account::read(words)?,
                            stored_change: Account::read(words)?,
                            balance_defect: Account::read(words)?,
                            absolute_defect: Account::read(words)?,
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
                port_energy_in_j: Account::read(words)?,
                stored_change_j: Account::read(words)?,
                dissipated_j: Account::read(words)?,
                numerical_loss_j: Account::read(words)?,
                entropy_generated_j_per_k: Account::read(words)?,
                balance_defect_j: Account::read(words)?,
                absolute_defect_j: Account::read(words)?,
            });
        }
        if words.count()? != self.connections.len() {
            return None;
        }
        let mut connections = Vec::with_capacity(self.connections.len());
        for _ in &self.connections {
            connections.push(ConnectionAccount {
                carried_defect: words.accounts()?,
                mass_defect_kg: Account::read(words)?,
                energy_defect_j: Account::read(words)?,
                entropy_generated_j_per_k: Account::read(words)?,
                absolute_energy_defect_j: Account::read(words)?,
                absolute_mass_defect_kg: Account::read(words)?,
            });
        }
        if words.count()? != self.boundaries.len() {
            return None;
        }
        let mut boundaries = Vec::with_capacity(self.boundaries.len());
        for shape in &self.boundaries {
            boundaries.push(BoundaryAccount {
                carried_in: words.accounts()?,
                port: shape.port.clone(),
                mass_in_kg: Account::read(words)?,
                energy_in_j: Account::read(words)?,
                entropy_in_j_per_k: Account::read(words)?,
            });
        }
        if words.count()? != self.sinks.len() {
            return None;
        }
        let mut sinks = Vec::with_capacity(self.sinks.len());
        for shape in &self.sinks {
            sinks.push(SinkAccount {
                name: shape.name.clone(),
                heat_j: Account::read(words)?,
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
        let mut stored = Vec::new();
        for node in &self.nodes {
            if let Some(account) = node.carried.get(carried) {
                let account = account.as_ref()?;
                if account.audited_steps != self.steps {
                    return None;
                }
                stored.push(&account.stored_change);
            }
        }
        let entered = total(
            self.boundaries
                .iter()
                .filter_map(|b| b.carried_in.get(carried)),
        );
        let defect = total(
            self.connections
                .iter()
                .filter_map(|c| c.carried_defect.get(carried)),
        );
        Some((total(stored) - entered - defect).to_f64())
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
        let mut stored = Vec::new();
        for node in &self.nodes {
            if node.has_streams {
                let account = node.mass.as_ref()?;
                if account.audited_steps != self.steps {
                    return None;
                }
                stored.push(&account.stored_change_kg);
            }
        }
        let entered = total(self.boundaries.iter().map(|b| &b.mass_in_kg));
        let defect = total(self.connections.iter().map(|c| &c.mass_defect_kg));
        Some((total(stored) - entered - defect).to_f64())
    }

    fn boundary_energy_in(&self) -> Dd {
        total(self.boundaries.iter().map(|account| &account.energy_in_j))
    }

    fn stored_change(&self) -> Dd {
        total(self.nodes.iter().map(|account| &account.stored_change_j))
    }

    fn sink_heat(&self) -> Dd {
        total(self.sinks.iter().map(|account| &account.heat_j))
    }

    fn numerical_loss(&self) -> Dd {
        total(self.nodes.iter().map(|account| &account.numerical_loss_j))
    }

    fn connection_defect(&self) -> Dd {
        total(
            self.connections
                .iter()
                .map(|account| &account.energy_defect_j),
        )
    }

    /// Energy that entered the system through its boundaries [J].
    #[must_use]
    pub fn boundary_energy_in_j(&self) -> f64 {
        self.boundary_energy_in().to_f64()
    }

    /// Change of the energy stored in all components [J].
    #[must_use]
    pub fn stored_change_j(&self) -> f64 {
        self.stored_change().to_f64()
    }

    /// Heat delivered to sinks outside the system [J].
    #[must_use]
    pub fn sink_heat_j(&self) -> f64 {
        self.sink_heat().to_f64()
    }

    /// Numerical loss of all components [J].
    #[must_use]
    pub fn numerical_loss_j(&self) -> f64 {
        self.numerical_loss().to_f64()
    }

    /// Interface residual energy of all connections [J].
    #[must_use]
    pub fn connection_defect_j(&self) -> f64 {
        self.connection_defect().to_f64()
    }

    /// Entropy generated in components and connections [J/K].
    #[must_use]
    pub fn entropy_generated_j_per_k(&self) -> f64 {
        let nodes = total(
            self.nodes
                .iter()
                .map(|account| &account.entropy_generated_j_per_k),
        );
        let connections = total(
            self.connections
                .iter()
                .map(|account| &account.entropy_generated_j_per_k),
        );
        (nodes + connections).to_f64()
    }

    /// First-law residual of the window [J]:
    /// `stored change + numerical loss + sink heat - boundary energy in
    ///  - connection defect`.
    ///
    /// It equals minus the sum of the per-step balance defects the component
    /// audits accepted (each defect is port energy in, minus stored change,
    /// minus numerical loss), so it is roundoff for exact components and it
    /// shows a leak that a component hid behind a wide tolerance. It is
    /// taken in double-double from the accounts and rounded once, so a
    /// residual far below the rounding of a large account is kept.
    #[must_use]
    pub fn energy_residual_j(&self) -> f64 {
        (self.stored_change() + self.numerical_loss() + self.sink_heat()
            - self.boundary_energy_in()
            - self.connection_defect())
        .to_f64()
    }

    /// The largest absolute energy defect of any component or connection
    /// \[J\]: the sum over steps of the absolute value of each step's accepted
    /// defect, so defects of opposite sign that cancel in the signed accounts
    /// show here. Zero before a step. It is a report; the kernel refuses no
    /// step on it.
    #[must_use]
    pub fn largest_absolute_energy_defect_j(&self) -> f64 {
        self.nodes
            .iter()
            .map(|account| account.absolute_defect_j.value())
            .chain(
                self.connections
                    .iter()
                    .map(|account| account.absolute_energy_defect_j.value()),
            )
            .fold(0.0, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::{Account, Sum, Words};

    #[test]
    fn a_small_term_survives_beside_large_ones() {
        let mut sum = Sum::default();
        sum.add(1.0e16);
        sum.add(1.0);
        sum.add(-1.0e16);
        assert_eq!(sum.value().to_bits(), 1.0_f64.to_bits());
        let mut account = Account::default();
        account.add(1.0e16);
        account.add(1.0);
        account.add(-1.0e16);
        assert_eq!(account.value().to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn an_account_reads_back_its_words_and_refuses_a_pair_it_cannot_hold() {
        let mut account = Account::default();
        account.add(1.0e6);
        account.add(3.0e-11);
        assert!(account.dd().lo != 0.0, "the trailing word is in use");
        let mut words = Vec::new();
        account.push_words(&mut words);
        assert_eq!(Account::read(&mut Words(&words)), Some(account));
        // hi + lo would round to another hi: not normalized.
        let unnormal = [1.0_f64.to_bits(), 1.0_f64.to_bits()];
        assert_eq!(Account::read(&mut Words(&unnormal)), None);
        assert_eq!(Account::read(&mut Words(&words[..1])), None, "one word");
    }
}
