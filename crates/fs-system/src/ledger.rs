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

use fs_evidence::Color;
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

/// Independently measured entropy account of one component that declares
/// its stored entropy (`Component::stored_entropy_j_per_k`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EntropyAccount {
    /// Net entropy entering the component through its ports [J/K].
    pub port_entropy_in_j_per_k: Account,
    /// Entropy the receipts of the audited steps say was generated [J/K].
    pub entropy_generated_j_per_k: Account,
    /// Change in independently measured stored entropy [J/K].
    pub stored_change_j_per_k: Account,
    /// Sum of the accepted defects, port entropy in plus generation minus
    /// stored change [J/K].
    pub balance_defect_j_per_k: Account,
    /// Sum over steps of the absolute value of each step's defect \[J/K\].
    pub absolute_defect_j_per_k: Account,
    /// Accepted steps for which the component supplied both measurements.
    pub audited_steps: u64,
}

impl EntropyAccount {
    /// A present account from the words `Ledger::push_words` wrote for it.
    fn read(words: &mut Words<'_>) -> Option<Self> {
        Some(Self {
            port_entropy_in_j_per_k: Account::read(words)?,
            entropy_generated_j_per_k: Account::read(words)?,
            stored_change_j_per_k: Account::read(words)?,
            balance_defect_j_per_k: Account::read(words)?,
            absolute_defect_j_per_k: Account::read(words)?,
            audited_steps: words.u64()?,
        })
    }
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
    /// Measured stored entropy, or unavailable: `None` for a component that
    /// does not declare it. A partial account has fewer audited steps than
    /// the ledger.
    pub entropy: Option<EntropyAccount>,
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
            match &node.entropy {
                None => words.push(0),
                Some(entropy) => {
                    words.push(1);
                    for account in [
                        &entropy.port_entropy_in_j_per_k,
                        &entropy.entropy_generated_j_per_k,
                        &entropy.stored_change_j_per_k,
                        &entropy.balance_defect_j_per_k,
                        &entropy.absolute_defect_j_per_k,
                    ] {
                        account.push_words(words);
                    }
                    words.push(entropy.audited_steps);
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
            let entropy = if words.flag()? {
                Some(EntropyAccount::read(words)?)
            } else {
                None
            };
            nodes.push(NodeAccount {
                name: shape.name.clone(),
                has_streams: shape.has_streams,
                mass,
                carried,
                entropy,
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

    /// Physical dissipation summed across all components [J].
    ///
    /// Numerical loss is tracked separately in [`Self::numerical_loss_j`] and
    /// is never added to physical dissipation.
    #[must_use]
    pub fn total_dissipated_j(&self) -> f64 {
        total(self.nodes.iter().map(|account| &account.dissipated_j)).to_f64()
    }

    /// Number of carried stream coordinates present across all components.
    #[must_use]
    pub fn carried_count(&self) -> usize {
        self.nodes
            .iter()
            .map(|node| node.carried.len())
            .max()
            .unwrap_or(0)
    }

    /// Attribute the largest per-component energy balance defect in the window
    /// to its owning component.
    ///
    /// Ranks components by their unsigned `absolute_defect_j` (breaking ties by
    /// `|balance_defect_j|` and then declaration order) so opposite-sign step
    /// defects cannot hide an owning component. Returns `None` when no
    /// component has a nonzero energy balance defect.
    #[must_use]
    pub fn worst_energy_leak_node(&self) -> Option<EnergyLeakAttribution<'_>> {
        let mut best: Option<EnergyLeakAttribution<'_>> = None;
        for node in &self.nodes {
            let signed = node.balance_defect_j.value();
            let abs_sum = node.absolute_defect_j.value();
            let key = abs_sum.max(signed.abs());
            if !(key > 0.0) && signed.is_finite() && abs_sum.is_finite() {
                continue;
            }
            let replace = match &best {
                None => true,
                Some(prev) => {
                    let prev_key = prev.absolute_defect_j.max(prev.balance_defect_j.abs());
                    !key.is_finite()
                        || key > prev_key
                        || (key.total_cmp(&prev_key).is_eq()
                            && signed.abs() > prev.balance_defect_j.abs())
                }
            };
            if replace {
                best = Some(EnergyLeakAttribution {
                    node: node.name.as_str(),
                    balance_defect_j: signed,
                    absolute_defect_j: abs_sum,
                    numerical_loss_j: node.numerical_loss_j.value(),
                });
            }
        }
        best
    }

    /// Evaluate the window balances and attach honest `fs-evidence` colors.
    ///
    /// Closure to IEEE-754 roundoff with zero numerical loss and complete
    /// inventory coverage yields [`Color::Verified`]; any numerical loss,
    /// non-roundoff balance defect, missing inventory measurement, or empty
    /// window stays [`Color::Estimated`] and is never laundered into a
    /// certificate.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn outcome(&self) -> WindowLedgerOutcome {
        const ROUNDOFF_ULPS: f64 = 128.0;

        let leak = self.worst_energy_leak_node();
        let worst_energy_leak_node = leak.as_ref().map(|entry| entry.node.to_owned());
        let worst_energy_leak_j = leak.as_ref().map_or(0.0, |entry| {
            entry.absolute_defect_j.max(entry.balance_defect_j.abs())
        });

        if self.steps == 0 {
            let empty = Color::Estimated {
                estimator: "window-ledger:empty".to_owned(),
                dispersion: f64::INFINITY,
            };
            let mass = self
                .nodes
                .iter()
                .any(|node| node.has_streams)
                .then(|| empty.clone());
            let carried = (0..self.carried_count())
                .map(|index| CarriedInventoryOutcome {
                    index,
                    residual: None,
                    color: empty.clone(),
                })
                .collect();
            return WindowLedgerOutcome {
                overall: empty.clone(),
                energy: empty.clone(),
                entropy: empty,
                mass,
                carried,
                worst_energy_leak_node,
                worst_energy_leak_j,
            };
        }

        let steps_f = self.steps as f64;
        let sqrt_steps = steps_f.sqrt();

        // 1. First law (energy).
        let mut energy_scale = 0.0_f64;
        for node in &self.nodes {
            energy_scale += node.port_energy_in_j.value().abs()
                + node.stored_change_j.value().abs()
                + node.dissipated_j.value().abs()
                + node.numerical_loss_j.value().abs();
        }
        for boundary in &self.boundaries {
            energy_scale += boundary.energy_in_j.value().abs();
        }
        for sink in &self.sinks {
            energy_scale += sink.heat_j.value().abs();
        }
        let energy_scale = if energy_scale.is_finite() {
            energy_scale.max(1.0)
        } else {
            f64::INFINITY
        };
        let energy_residual = self.energy_residual_j();
        let max_abs_energy_defect = self.largest_absolute_energy_defect_j();
        let total_num_loss = self.numerical_loss_j().abs();
        let max_node_num_loss = self
            .nodes
            .iter()
            .map(|node| node.numerical_loss_j.value().abs())
            .fold(0.0_f64, f64::max);
        let min_node_dissipated = self
            .nodes
            .iter()
            .map(|node| node.dissipated_j.value())
            .fold(
                0.0_f64,
                |acc, v| if v.is_nan() || v < acc { v } else { acc },
            );
        let signed_energy_bound = ROUNDOFF_ULPS * f64::EPSILON * energy_scale * sqrt_steps;
        let abs_energy_bound = ROUNDOFF_ULPS * f64::EPSILON * energy_scale * steps_f;

        let energy = if !min_node_dissipated.is_finite() || min_node_dissipated < 0.0 {
            Color::Estimated {
                estimator: "window-ledger:negative-dissipation".to_owned(),
                dispersion: f64::INFINITY,
            }
        } else if !total_num_loss.is_finite()
            || !max_node_num_loss.is_finite()
            || total_num_loss > 0.0
            || max_node_num_loss > 0.0
        {
            let spread = if total_num_loss.is_finite()
                && max_node_num_loss.is_finite()
                && energy_scale.is_finite()
            {
                (total_num_loss.max(max_node_num_loss) / energy_scale).max(0.0)
            } else {
                f64::INFINITY
            };
            Color::Estimated {
                estimator: "window-ledger:numerical-loss".to_owned(),
                dispersion: spread,
            }
        } else if !energy_residual.is_finite()
            || !max_abs_energy_defect.is_finite()
            || !signed_energy_bound.is_finite()
            || energy_residual.abs() > signed_energy_bound
            || max_abs_energy_defect > abs_energy_bound
        {
            let spread = if energy_residual.is_finite()
                && max_abs_energy_defect.is_finite()
                && energy_scale.is_finite()
            {
                (energy_residual.abs().max(max_abs_energy_defect) / energy_scale).max(0.0)
            } else {
                f64::INFINITY
            };
            Color::Estimated {
                estimator: "window-ledger:energy-defect".to_owned(),
                dispersion: spread,
            }
        } else {
            let half = signed_energy_bound.max(energy_residual.abs());
            Color::Verified {
                lo: energy_residual - half,
                hi: energy_residual + half,
            }
        };

        // 2. Second law (non-negative entropy production per component and total).
        let total_entropy = self.entropy_generated_j_per_k();
        let mut entropy_scale = total_entropy.abs();
        let mut min_node_entropy = 0.0_f64;
        for node in &self.nodes {
            let v = node.entropy_generated_j_per_k.value();
            entropy_scale += v.abs();
            if v.is_nan() || v < min_node_entropy {
                min_node_entropy = v;
            }
        }
        let mut min_conn_entropy = 0.0_f64;
        for connection in &self.connections {
            let v = connection.entropy_generated_j_per_k.value();
            entropy_scale += v.abs();
            if v.is_nan() || v < min_conn_entropy {
                min_conn_entropy = v;
            }
        }
        let entropy_scale = if entropy_scale.is_finite() {
            entropy_scale.max(1.0)
        } else {
            f64::INFINITY
        };
        let entropy_bound = ROUNDOFF_ULPS * f64::EPSILON * entropy_scale * sqrt_steps;
        let entropy = if !total_entropy.is_finite()
            || !min_node_entropy.is_finite()
            || !min_conn_entropy.is_finite()
            || !entropy_bound.is_finite()
            || min_node_entropy < 0.0
            || min_conn_entropy < -entropy_bound
            || total_entropy < -entropy_bound
        {
            Color::Estimated {
                estimator: "window-ledger:entropy-defect".to_owned(),
                dispersion: f64::INFINITY,
            }
        } else {
            let lo = (total_entropy - entropy_bound).max(0.0);
            let hi = (total_entropy + entropy_bound).max(lo);
            Color::Verified { lo, hi }
        };

        // 3. Stream mass balance (when stream ports exist).
        let mass = if self.nodes.iter().any(|node| node.has_streams) {
            let color = match self.mass_residual_kg() {
                None => Color::Estimated {
                    estimator: "window-ledger:mass-unmeasured".to_owned(),
                    dispersion: f64::INFINITY,
                },
                Some(residual) => {
                    let mut mass_scale = 0.0_f64;
                    let mut max_abs_mass_defect = 0.0_f64;
                    for node in &self.nodes {
                        if let Some(m) = &node.mass {
                            mass_scale +=
                                m.port_mass_in_kg.value().abs() + m.stored_change_kg.value().abs();
                            max_abs_mass_defect =
                                max_abs_mass_defect.max(m.absolute_defect_kg.value());
                        }
                    }
                    for boundary in &self.boundaries {
                        mass_scale += boundary.mass_in_kg.value().abs();
                    }
                    for connection in &self.connections {
                        max_abs_mass_defect =
                            max_abs_mass_defect.max(connection.absolute_mass_defect_kg.value());
                    }
                    let mass_scale = if mass_scale.is_finite() {
                        mass_scale.max(1.0)
                    } else {
                        f64::INFINITY
                    };
                    let signed_bound = ROUNDOFF_ULPS * f64::EPSILON * mass_scale * sqrt_steps;
                    let abs_bound = ROUNDOFF_ULPS * f64::EPSILON * mass_scale * steps_f;
                    if residual.is_finite()
                        && max_abs_mass_defect.is_finite()
                        && signed_bound.is_finite()
                        && residual.abs() <= signed_bound
                        && max_abs_mass_defect <= abs_bound
                    {
                        let half = signed_bound.max(residual.abs());
                        Color::Verified {
                            lo: residual - half,
                            hi: residual + half,
                        }
                    } else {
                        let spread = if residual.is_finite()
                            && max_abs_mass_defect.is_finite()
                            && mass_scale.is_finite()
                        {
                            (residual.abs().max(max_abs_mass_defect) / mass_scale).max(0.0)
                        } else {
                            f64::INFINITY
                        };
                        Color::Estimated {
                            estimator: "window-ledger:mass-defect".to_owned(),
                            dispersion: spread,
                        }
                    }
                }
            };
            Some(color)
        } else {
            None
        };

        // 4. Carried conserved inventories (species, elements, momentum, charge).
        let carried_len = self.carried_count();
        let mut carried = Vec::with_capacity(carried_len);
        for index in 0..carried_len {
            let residual = self.carried_residual(index);
            let color = match residual {
                None => Color::Estimated {
                    estimator: format!("window-ledger:carried-{index}-unmeasured"),
                    dispersion: f64::INFINITY,
                },
                Some(res) => {
                    let mut scale = 0.0_f64;
                    let mut max_abs = 0.0_f64;
                    for node in &self.nodes {
                        if let Some(Some(acc)) = node.carried.get(index) {
                            scale += acc.port_in.value().abs() + acc.stored_change.value().abs();
                            max_abs = max_abs.max(acc.absolute_defect.value());
                        }
                    }
                    for boundary in &self.boundaries {
                        if let Some(acc) = boundary.carried_in.get(index) {
                            scale += acc.value().abs();
                        }
                    }
                    for connection in &self.connections {
                        if let Some(acc) = connection.carried_defect.get(index) {
                            max_abs = max_abs.max(acc.value().abs());
                        }
                    }
                    let scale = if scale.is_finite() {
                        scale.max(1.0)
                    } else {
                        f64::INFINITY
                    };
                    let signed_bound = ROUNDOFF_ULPS * f64::EPSILON * scale * sqrt_steps;
                    let abs_bound = ROUNDOFF_ULPS * f64::EPSILON * scale * steps_f;
                    if res.is_finite()
                        && max_abs.is_finite()
                        && signed_bound.is_finite()
                        && res.abs() <= signed_bound
                        && max_abs <= abs_bound
                    {
                        let half = signed_bound.max(res.abs());
                        Color::Verified {
                            lo: res - half,
                            hi: res + half,
                        }
                    } else {
                        let spread = if res.is_finite() && max_abs.is_finite() && scale.is_finite()
                        {
                            (res.abs().max(max_abs) / scale).max(0.0)
                        } else {
                            f64::INFINITY
                        };
                        Color::Estimated {
                            estimator: format!("window-ledger:carried-{index}-defect"),
                            dispersion: spread,
                        }
                    }
                }
            };
            carried.push(CarriedInventoryOutcome {
                index,
                residual,
                color,
            });
        }

        // 5. Conservative overall window color (no laundering).
        let mut worst_est: Option<(String, f64)> = None;
        for color in core::iter::once(&energy)
            .chain(core::iter::once(&entropy))
            .chain(mass.as_ref())
            .chain(carried.iter().map(|item| &item.color))
        {
            if let Color::Estimated {
                estimator,
                dispersion,
            } = color
            {
                let replace = match &worst_est {
                    None => true,
                    Some((_, prev)) => {
                        dispersion.is_infinite() || (prev.is_finite() && *dispersion > *prev)
                    }
                };
                if replace {
                    worst_est = Some((estimator.clone(), *dispersion));
                }
            }
        }
        let overall = match worst_est {
            None => energy.clone(),
            Some((estimator, dispersion)) => Color::Estimated {
                estimator,
                dispersion,
            },
        };

        WindowLedgerOutcome {
            overall,
            energy,
            entropy,
            mass,
            carried,
            worst_energy_leak_node,
            worst_energy_leak_j,
        }
    }

    /// Build a deterministic export package for `fs-ledger` and run reports.
    ///
    /// Contains the evidence-colored [`WindowLedgerOutcome`], flat finite metric
    /// rows ready for `fs-ledger`'s `metrics` table, canonical JSON artifact
    /// bytes, and a Markdown summary table for run reports.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn export(&self) -> LedgerExport {
        use core::fmt::Write as _;

        let outcome = self.outcome();
        let mut metrics = Vec::with_capacity(12 + self.nodes.len() * 5);
        let mut push_metric = |name: String, value: f64| {
            if value.is_finite() {
                metrics.push((name, if value == 0.0 { 0.0 } else { value }));
            }
        };
        push_metric("window.steps".to_owned(), self.steps as f64);
        push_metric("window.duration_s".to_owned(), self.duration_s.value());
        push_metric(
            "window.energy.stored_change_j".to_owned(),
            self.stored_change_j(),
        );
        push_metric(
            "window.energy.boundary_in_j".to_owned(),
            self.boundary_energy_in_j(),
        );
        push_metric("window.energy.sink_heat_j".to_owned(), self.sink_heat_j());
        push_metric(
            "window.energy.dissipated_j".to_owned(),
            self.total_dissipated_j(),
        );
        push_metric(
            "window.energy.numerical_loss_j".to_owned(),
            self.numerical_loss_j(),
        );
        push_metric(
            "window.energy.connection_defect_j".to_owned(),
            self.connection_defect_j(),
        );
        push_metric(
            "window.energy.residual_j".to_owned(),
            self.energy_residual_j(),
        );
        push_metric(
            "window.energy.largest_abs_defect_j".to_owned(),
            self.largest_absolute_energy_defect_j(),
        );
        push_metric(
            "window.entropy.generated_j_per_k".to_owned(),
            self.entropy_generated_j_per_k(),
        );
        if let Some(mass_res) = self.mass_residual_kg() {
            push_metric("window.mass.residual_kg".to_owned(), mass_res);
        }
        for item in &outcome.carried {
            if let Some(res) = item.residual {
                push_metric(format!("window.carried.{}.residual", item.index), res);
            }
        }
        for node in &self.nodes {
            push_metric(
                format!("node.{}.port_energy_in_j", node.name),
                node.port_energy_in_j.value(),
            );
            push_metric(
                format!("node.{}.stored_change_j", node.name),
                node.stored_change_j.value(),
            );
            push_metric(
                format!("node.{}.dissipated_j", node.name),
                node.dissipated_j.value(),
            );
            push_metric(
                format!("node.{}.numerical_loss_j", node.name),
                node.numerical_loss_j.value(),
            );
            push_metric(
                format!("node.{}.balance_defect_j", node.name),
                node.balance_defect_j.value(),
            );
            push_metric(
                format!("node.{}.absolute_defect_j", node.name),
                node.absolute_defect_j.value(),
            );
        }

        // Canonical JSON artifact.
        let mut json = String::with_capacity(512);
        json.push_str("{\"schema\":\"frankensim.window_ledger.v1\",\"steps\":");
        let _ = write!(json, "{}", self.steps);
        json.push_str(",\"duration_s\":");
        push_json_f64(&mut json, self.duration_s.value());
        json.push_str(",\"overall_color\":");
        push_color_json(&mut json, &outcome.overall);
        json.push_str(",\"energy_color\":");
        push_color_json(&mut json, &outcome.energy);
        json.push_str(",\"entropy_color\":");
        push_color_json(&mut json, &outcome.entropy);
        json.push_str(",\"mass_color\":");
        match &outcome.mass {
            Some(c) => push_color_json(&mut json, c),
            None => json.push_str("null"),
        }
        json.push_str(",\"carried_colors\":[");
        for (i, item) in outcome.carried.iter().enumerate() {
            if i > 0 {
                json.push(',');
            }
            let _ = write!(json, "{{\"index\":{},\"residual\":", item.index);
            match item.residual {
                Some(v) => push_json_f64(&mut json, v),
                None => json.push_str("null"),
            }
            json.push_str(",\"color\":");
            push_color_json(&mut json, &item.color);
            json.push('}');
        }
        json.push_str("],\"totals\":{\"stored_change_j\":");
        push_json_f64(&mut json, self.stored_change_j());
        json.push_str(",\"boundary_energy_in_j\":");
        push_json_f64(&mut json, self.boundary_energy_in_j());
        json.push_str(",\"sink_heat_j\":");
        push_json_f64(&mut json, self.sink_heat_j());
        json.push_str(",\"dissipated_j\":");
        push_json_f64(&mut json, self.total_dissipated_j());
        json.push_str(",\"numerical_loss_j\":");
        push_json_f64(&mut json, self.numerical_loss_j());
        json.push_str(",\"connection_defect_j\":");
        push_json_f64(&mut json, self.connection_defect_j());
        json.push_str(",\"energy_residual_j\":");
        push_json_f64(&mut json, self.energy_residual_j());
        json.push_str(",\"largest_absolute_energy_defect_j\":");
        push_json_f64(&mut json, self.largest_absolute_energy_defect_j());
        json.push_str(",\"entropy_generated_j_per_k\":");
        push_json_f64(&mut json, self.entropy_generated_j_per_k());
        json.push_str(",\"mass_residual_kg\":");
        match self.mass_residual_kg() {
            Some(v) => push_json_f64(&mut json, v),
            None => json.push_str("null"),
        }
        json.push_str("},\"worst_energy_leak\":");
        match self.worst_energy_leak_node() {
            Some(leak) => {
                json.push_str("{\"node\":");
                push_json_str(&mut json, leak.node);
                json.push_str(",\"balance_defect_j\":");
                push_json_f64(&mut json, leak.balance_defect_j);
                json.push_str(",\"absolute_defect_j\":");
                push_json_f64(&mut json, leak.absolute_defect_j);
                json.push_str(",\"numerical_loss_j\":");
                push_json_f64(&mut json, leak.numerical_loss_j);
                json.push('}');
            }
            None => json.push_str("null"),
        }
        json.push_str(",\"nodes\":[");
        for (i, node) in self.nodes.iter().enumerate() {
            if i > 0 {
                json.push(',');
            }
            json.push_str("{\"name\":");
            push_json_str(&mut json, &node.name);
            json.push_str(",\"port_energy_in_j\":");
            push_json_f64(&mut json, node.port_energy_in_j.value());
            json.push_str(",\"stored_change_j\":");
            push_json_f64(&mut json, node.stored_change_j.value());
            json.push_str(",\"dissipated_j\":");
            push_json_f64(&mut json, node.dissipated_j.value());
            json.push_str(",\"numerical_loss_j\":");
            push_json_f64(&mut json, node.numerical_loss_j.value());
            json.push_str(",\"entropy_generated_j_per_k\":");
            push_json_f64(&mut json, node.entropy_generated_j_per_k.value());
            json.push_str(",\"balance_defect_j\":");
            push_json_f64(&mut json, node.balance_defect_j.value());
            json.push_str(",\"absolute_defect_j\":");
            push_json_f64(&mut json, node.absolute_defect_j.value());
            json.push('}');
        }
        json.push_str("]}");

        // Human-readable Markdown summary for run reports.
        let mut md = String::with_capacity(512);
        let _ = writeln!(
            md,
            "## Window Ledger ({} steps, {:.6e} s, overall: `{}`)",
            self.steps,
            self.duration_s.value(),
            outcome.overall.name()
        );
        let _ = writeln!(
            md,
            "- First law (`{}`): stored {:.6e} J, boundary in {:.6e} J, sink heat {:.6e} J, dissipated {:.6e} J, numerical loss {:.6e} J, connection defect {:.3e} J, residual {:.3e} J (max abs defect {:.3e} J)",
            outcome.energy.name(),
            self.stored_change_j(),
            self.boundary_energy_in_j(),
            self.sink_heat_j(),
            self.total_dissipated_j(),
            self.numerical_loss_j(),
            self.connection_defect_j(),
            self.energy_residual_j(),
            self.largest_absolute_energy_defect_j()
        );
        let _ = writeln!(
            md,
            "- Second law (`{}`): entropy generated {:.6e} J/K",
            outcome.entropy.name(),
            self.entropy_generated_j_per_k()
        );
        if let Some(ref mass_color) = outcome.mass {
            match self.mass_residual_kg() {
                Some(r) => {
                    let _ = writeln!(
                        md,
                        "- Stream mass (`{}`): residual {:.3e} kg",
                        mass_color.name(),
                        r
                    );
                }
                None => {
                    let _ = writeln!(
                        md,
                        "- Stream mass (`{}`): unmeasured / incomplete coverage",
                        mass_color.name()
                    );
                }
            }
        }
        for item in &outcome.carried {
            match item.residual {
                Some(r) => {
                    let _ = writeln!(
                        md,
                        "- Carried[{}] (`{}`): residual {:.3e}",
                        item.index,
                        item.color.name(),
                        r
                    );
                }
                None => {
                    let _ = writeln!(
                        md,
                        "- Carried[{}] (`{}`): unmeasured / incomplete coverage",
                        item.index,
                        item.color.name()
                    );
                }
            }
        }
        if let Some(leak) = self.worst_energy_leak_node() {
            let _ = writeln!(
                md,
                "- Worst component energy defect: `{}` (signed {:.3e} J, abs {:.3e} J, numerical loss {:.3e} J)",
                leak.node, leak.balance_defect_j, leak.absolute_defect_j, leak.numerical_loss_j
            );
        }
        md.push_str("\n| Component | Port In [J] | Stored [J] | Dissipated [J] | Num. Loss [J] | Defect [J] | Abs Defect [J] |\n");
        md.push_str("| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n");
        for node in &self.nodes {
            let _ = writeln!(
                md,
                "| {} | {:.6e} | {:.6e} | {:.6e} | {:.6e} | {:.3e} | {:.3e} |",
                node.name,
                node.port_energy_in_j.value(),
                node.stored_change_j.value(),
                node.dissipated_j.value(),
                node.numerical_loss_j.value(),
                node.balance_defect_j.value(),
                node.absolute_defect_j.value()
            );
        }

        LedgerExport {
            outcome,
            metrics,
            canonical_json: json,
            report_markdown: md,
        }
    }
}

/// Owning-component attribution for the worst energy balance defect in a window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyLeakAttribution<'a> {
    /// Declared component name.
    pub node: &'a str,
    /// Cumulative signed first-law balance defect [J].
    pub balance_defect_j: f64,
    /// Cumulative unsigned absolute step balance defect [J].
    pub absolute_defect_j: f64,
    /// Cumulative numerical loss declared by this component [J].
    pub numerical_loss_j: f64,
}

/// Evidence-colored closure outcome for one carried stream inventory coordinate.
#[derive(Debug, Clone, PartialEq)]
pub struct CarriedInventoryOutcome {
    /// Carried coordinate ordinal.
    pub index: usize,
    /// Window balance residual, or `None` if measurement coverage was incomplete.
    pub residual: Option<f64>,
    /// Evidence color (`Verified` on roundoff closure, `Estimated` otherwise).
    pub color: Color,
}

/// Evidence-colored window-ledger closure outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowLedgerOutcome {
    /// Conservative overall window color (`Verified` iff every applicable balance
    /// closes to roundoff with zero numerical loss and complete measurements).
    pub overall: Color,
    /// First-law energy balance color.
    pub energy: Color,
    /// Second-law entropy production color.
    pub entropy: Color,
    /// Stream-mass balance color (`Some` when any component has stream ports).
    pub mass: Option<Color>,
    /// Per-coordinate carried inventory closure colors (species, elements,
    /// momentum, charge).
    pub carried: Vec<CarriedInventoryOutcome>,
    /// Name of the owning component with the largest energy balance defect, if any.
    pub worst_energy_leak_node: Option<String>,
    /// Magnitude of that worst component energy defect [J] (`0.0` when none).
    pub worst_energy_leak_j: f64,
}

/// Deterministic export of a window [`Ledger`] for `fs-ledger` and run reports.
#[derive(Debug, Clone, PartialEq)]
pub struct LedgerExport {
    /// Evidence-colored closure summary.
    pub outcome: WindowLedgerOutcome,
    /// Flat finite `(metric_name, value)` rows ready for `fs-ledger::Ledger::record_metric`.
    pub metrics: Vec<(String, f64)>,
    /// Canonical JSON document (`schema: "frankensim.window_ledger.v1"`).
    pub canonical_json: String,
    /// Human-readable Markdown ledger table for run reports.
    pub report_markdown: String,
}

fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                use core::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn push_json_f64(out: &mut String, v: f64) {
    use core::fmt::Write as _;
    if v.is_finite() {
        let normalized = if v == 0.0 { 0.0 } else { v };
        let _ = write!(out, "{normalized:.17e}");
    } else {
        out.push_str("null");
    }
}

fn push_color_json(out: &mut String, color: &Color) {
    out.push_str("{\"color\":");
    push_json_str(out, color.name());
    out.push_str(",\"payload\":");
    out.push_str(&color.payload_json());
    out.push('}');
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
