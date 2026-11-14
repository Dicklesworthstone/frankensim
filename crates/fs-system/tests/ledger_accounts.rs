//! The ledger's accounts in double-double, and the absolute defect accounts
//! beside the signed ones (bead rc-1109 .5.4, item 5).
//!
//! A member without ports has no energy in, so the first-law defect the
//! audit accepts from it in a step is exactly minus the stored change its
//! receipt states. `Leak` states a planned change per step, which makes
//! every defect, and every sum of them, known exactly in the test.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::dd::Dd;
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_system::elements::{GasContents, GasVolume, Restriction};
use fs_system::ledger::{Account, Ledger};
use fs_system::{NodeId, Phase, System, SystemBuilder, SystemError};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let gate = CancelGate::new_clock_free();
        let key = StreamKey {
            seed: 29,
            kernel_id: 41,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// The tolerance `Leak` claims for its balance: wide enough for every
/// planned change, so the audit accepts each one as a defect.
const LEAK_TOLERANCE_J: f64 = 1.0e7;

#[derive(Debug)]
struct Refused;

/// A member without ports whose receipt states, step by step, the stored
/// changes of `plan`, and which refuses in `admit` the step that would make
/// its count `refuse_at`. Its state is its stored energy and its count.
struct Leak {
    plan: Vec<f64>,
    refuse_at: Option<u64>,
}

impl Leak {
    fn new(plan: Vec<f64>) -> Self {
        Self {
            plan,
            refuse_at: None,
        }
    }
}

impl Component for Leak {
    type State = (f64, u64);
    type Error = Refused;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &(f64, u64)) -> Result<f64, Refused> {
        Ok(state.0)
    }

    fn propose(
        &self,
        old: &(f64, u64),
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut (f64, u64),
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Refused> {
        let change = self.plan[usize::try_from(old.1).unwrap() % self.plan.len()];
        *new = (old.0 + change, old.1 + 1);
        Ok(Proposal::Step(Receipt::lossless(change, LEAK_TOLERANCE_J)))
    }

    fn admit(
        &self,
        _old: &(f64, u64),
        new: &(f64, u64),
        _receipt: &Receipt,
    ) -> Result<(), Refused> {
        if Some(new.1) == self.refuse_at {
            Err(Refused)
        } else {
            Ok(())
        }
    }
}

fn leaks(members: Vec<(&str, Leak)>) -> (System, Vec<NodeId>) {
    let mut builder = SystemBuilder::new();
    let nodes = members
        .into_iter()
        .map(|(name, leak)| builder.add(name, leak, (0.0, 0)).unwrap())
        .collect();
    (builder.seal().unwrap(), nodes)
}

fn steps(system: &mut System, count: usize) {
    with_cx(|cx| {
        for index in 0..count {
            system
                .step(1.0e-3, cx)
                .unwrap_or_else(|error| panic!("step {index}: {error:?}"));
        }
    });
}

/// A double-double as an exact integer count of `2^-scale`. Both words must
/// be multiples of that unit, which the plans below guarantee.
fn units(value: Dd, scale: i32) -> i128 {
    let count = |word: f64| {
        let scaled = word * 2.0_f64.powi(scale);
        assert_eq!(
            scaled.fract().abs().to_bits(),
            0,
            "{word:e} is not a multiple of 2^-{scale}"
        );
        scaled as i128
    };
    count(value.hi) + count(value.lo)
}

#[test]
fn defects_of_opposite_sign_cancel_in_the_signed_account_and_not_in_the_absolute_one() {
    let (mut system, _) = leaks(vec![("leak", Leak::new(vec![0.25, -0.25]))]);
    assert_eq!(
        system.ledger().largest_absolute_energy_defect_j().to_bits(),
        0.0_f64.to_bits()
    );
    steps(&mut system, 2);
    let ledger = system.ledger();
    let account = &ledger.nodes[0];
    // The defects were -0.25 J and +0.25 J.
    assert_eq!(
        account.balance_defect_j.value().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(ledger.energy_residual_j().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        account.absolute_defect_j.value().to_bits(),
        0.5_f64.to_bits()
    );
    assert_eq!(
        account.absolute_defect_j.dd().lo.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        ledger.largest_absolute_energy_defect_j().to_bits(),
        0.5_f64.to_bits()
    );
}

#[test]
fn a_refused_step_adds_nothing_to_any_account() {
    let leak = Leak {
        plan: vec![0.25, -0.125, 0.5],
        refuse_at: Some(3),
    };
    let (mut system, nodes) = leaks(vec![("leak", leak), ("other", Leak::new(vec![1.0]))]);
    steps(&mut system, 2);
    let before = format!("{:?}", system.ledger());
    let error = with_cx(|cx| system.step(1.0e-3, cx)).expect_err("the third step is refused");
    assert!(
        matches!(
            &error,
            SystemError::Component { node, phase: Phase::Admit, .. } if node == "leak"
        ),
        "{error:?}"
    );
    assert_eq!(format!("{:?}", system.ledger()), before);
    assert_eq!(system.state::<(f64, u64)>(nodes[0]).unwrap(), &(0.125, 2));
    let account = &system.ledger().nodes[0];
    assert_eq!(
        account.absolute_defect_j.value().to_bits(),
        0.375_f64.to_bits()
    );
    assert_eq!(system.ledger().steps, 2);
}

/// splitmix64: the same numbers on every build.
struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// A signed whole number below 2^30 in magnitude.
    fn signed(&mut self) -> i64 {
        let value = (self.next() >> 34) as i64;
        if self.next() & 1 == 0 { value } else { -value }
    }
}

#[test]
fn the_absolute_accounts_are_the_exact_sum_of_each_steps_absolute_defect() {
    // 200 planned changes, alternately whole multiples of 2^-10 J and of
    // 2^-70 J, each below 2^30 units: the sums span 97 bits, more than an
    // f64 holds and less than a double-double does, so the accounts must
    // equal the exact sums, counted below in units of 2^-70 J.
    let mut draw = Draw(0x5eed_1e57);
    let mut plan = Vec::new();
    let (mut signed, mut absolute) = (0_i128, 0_i128);
    for index in 0..200 {
        let whole = draw.signed();
        let shift = if index % 2 == 0 { 60 } else { 0 };
        plan.push(whole as f64 * 2.0_f64.powi(shift - 70));
        // The defect is minus the change.
        signed -= i128::from(whole) << shift;
        absolute += i128::from(whole.abs()) << shift;
    }
    let (mut system, _) = leaks(vec![("leak", Leak::new(plan))]);
    steps(&mut system, 200);
    let account = &system.ledger().nodes[0];
    assert_eq!(units(account.balance_defect_j.dd(), 70), signed);
    assert_eq!(units(account.absolute_defect_j.dd(), 70), absolute);
    // The residual is minus the defects, rounded once.
    let exact = |count: i128| count as f64 * 2.0_f64.powi(-70);
    assert_eq!(
        system.ledger().energy_residual_j().to_bits(),
        exact(-signed).to_bits()
    );
    eprintln!(
        "200 defects: signed {:e} J, absolute {:e} J, both to the last of 97 bits",
        exact(signed),
        exact(absolute)
    );
}

#[test]
fn a_residual_below_the_rounding_of_its_accounts_survives_in_double_double() {
    // Two members store about 1e6 J each, one gaining and one losing, and
    // the first also 2^-40 J: the window's residual is 2^-40 J, about 9.1e-13
    // J, far below eps * 1e6 J (about 2.2e-10 J).
    let (mut system, _) = leaks(vec![
        ("store", Leak::new(vec![1.0e6, 2.0_f64.powi(-40)])),
        ("drain", Leak::new(vec![-1.0e6, 0.0])),
    ]);
    steps(&mut system, 2);
    let ledger: &Ledger = system.ledger();
    // Exactly: (1e6 * 2^40 + 1 - 1e6 * 2^40) units of 2^-40 J.
    let whole = 1_000_000_i128 << 40;
    let exact = (whole + 1 - whole) as f64 * 2.0_f64.powi(-40);
    assert_eq!(ledger.energy_residual_j().to_bits(), exact.to_bits());
    assert_eq!(ledger.stored_change_j().to_bits(), exact.to_bits());
    // Rounding each account to f64 first, as the ledger did before, loses
    // it: 1e6 + 2^-40 rounds to 1e6.
    let rounded_first: f64 = ledger
        .nodes
        .iter()
        .map(|account| account.stored_change_j.value())
        .sum();
    assert_eq!(rounded_first.to_bits(), 0.0_f64.to_bits());
    eprintln!(
        "residual of two 1e6 J accounts: {:e} J in double-double, {:e} J rounding each \
         account first; exact {exact:e} J",
        ledger.energy_residual_j(),
        rounded_first
    );
}

const AIR_J_KG_K: f64 = 287.0;
const RATIO: f64 = 1.4;

/// 800 W into `a` for the first second, then 4 kW into `b` for the next.
fn alternate_heat(t: f64) -> [f64; 2] {
    if t < 1.0 {
        [800.0, 0.0]
    } else if t < 2.0 {
        [0.0, 4_000.0]
    } else {
        [0.0, 0.0]
    }
}

/// The flow that turns round of tests/elements.rs: two vessels of air
/// carrying two quantities, joined by a restriction.
fn vessels() -> System {
    let air = || GasVolume::new(AIR_J_KG_K, RATIO, 2).unwrap();
    let start: [GasContents; 2] = [
        air().contents(2.0e5, 300.0, 0.01, &[0.3, 0.0]).unwrap(),
        air().contents(2.0e5, 300.0, 0.03, &[0.0, 0.6]).unwrap(),
    ];
    let mut builder = SystemBuilder::new();
    let pipe = builder
        .add_with_tangent("pipe", Restriction::for_gas(2.0e-8, &air()).unwrap(), ())
        .unwrap();
    for (side, (name, contents)) in ["a", "b"].into_iter().zip(start).enumerate() {
        let node = builder.add_with_tangent(name, air(), contents).unwrap();
        builder
            .bond(
                builder.port(node, "stream").unwrap(),
                builder.port(pipe, name).unwrap(),
            )
            .unwrap();
        let port = builder.port(node, "heat").unwrap();
        builder
            .boundary(
                port,
                Box::new(move |t0, dt, values| values[0] = alternate_heat(t0 + 0.5 * dt)[side]),
            )
            .unwrap();
    }
    builder.seal().unwrap()
}

/// Every pair of a signed defect account and its absolute account in the
/// ledger: per component energy, mass and each carried coordinate, and per
/// connection energy and mass.
fn pairs(ledger: &Ledger) -> Vec<(String, Account, Account)> {
    let mut out = Vec::new();
    for node in &ledger.nodes {
        out.push((
            format!("{} energy", node.name),
            node.balance_defect_j,
            node.absolute_defect_j,
        ));
        if let Some(mass) = &node.mass {
            out.push((
                format!("{} mass", node.name),
                mass.balance_defect_kg,
                mass.absolute_defect_kg,
            ));
        }
        for (carried, account) in node.carried.iter().enumerate() {
            if let Some(account) = account {
                out.push((
                    format!("{} carried {carried}", node.name),
                    account.balance_defect,
                    account.absolute_defect,
                ));
            }
        }
    }
    for (index, connection) in ledger.connections.iter().enumerate() {
        out.push((
            format!("connection {index} energy"),
            connection.energy_defect_j,
            connection.absolute_energy_defect_j,
        ));
        out.push((
            format!("connection {index} mass"),
            connection.mass_defect_kg,
            connection.absolute_mass_defect_kg,
        ));
    }
    out
}

/// One pair of accounts as the judge follows it from step to step.
struct Row {
    name: String,
    signed: Account,
    absolute: Account,
    judge: Dd,
    allowance: f64,
    moved: u32,
}

#[test]
fn every_absolute_account_of_a_turning_flow_sums_the_size_of_each_steps_defect() {
    // The judge recovers each step's defect as the change of the signed
    // account, in double-double, and sums its size. A double-double addition
    // errs by at most 2^-104 of its result, so the change differs from the
    // step's defect by at most 2^-104 (|old| + |new| + |change|), and the two
    // sums of sizes by at most 2^-104 of each running total per step: twice
    // that is the allowance. Mass and carried accounts appear at the first
    // step that audits them; before it they count as zero.
    let unit = 2.0_f64.powi(-104);
    let mut system = vessels();
    let mut rows: Vec<Row> = Vec::new();
    with_cx(|cx| {
        for index in 0..600 {
            system
                .step(3.0 / 600.0, cx)
                .unwrap_or_else(|error| panic!("step {index}: {error:?}"));
            for (name, new, absolute) in pairs(system.ledger()) {
                if !rows.iter().any(|row| row.name == name) {
                    rows.push(Row {
                        name: name.clone(),
                        signed: Account::default(),
                        absolute: Account::default(),
                        judge: Dd::ZERO,
                        allowance: 0.0,
                        moved: 0,
                    });
                }
                let row = rows
                    .iter_mut()
                    .find(|row| row.name == name)
                    .expect("the row is there");
                let old = row.signed;
                let change = new.dd() - old.dd();
                if change.hi != 0.0 {
                    row.moved += 1;
                }
                row.judge = row.judge + change.abs();
                row.allowance += 2.0
                    * unit
                    * (old.value().abs()
                        + new.value().abs()
                        + change.hi.abs()
                        + absolute.value().abs()
                        + row.judge.hi.abs());
                row.signed = new;
                row.absolute = absolute;
            }
        }
    });
    let mut worst = 0.0_f64;
    for row in &rows {
        let off = (row.absolute.dd() - row.judge).hi.abs();
        assert!(
            off <= row.allowance,
            "{}: absolute {:e}, judge {:e}, off {off:e} beyond {:e}",
            row.name,
            row.absolute.value(),
            row.judge.hi,
            row.allowance
        );
        assert!(
            row.absolute.value() >= row.signed.value().abs(),
            "{}",
            row.name
        );
        if row.allowance > 0.0 {
            worst = worst.max(off / row.allowance);
        }
        eprintln!(
            "{}: signed {:e}, absolute {:e}, {} of 600 steps moved it",
            row.name,
            row.signed.value(),
            row.absolute.value(),
            row.moved
        );
    }
    // Rows that moved in some step, connections among them, so the check is
    // not of zeros.
    let moved = |prefix: &str| {
        rows.iter()
            .filter(|row| row.name.starts_with(prefix) && row.moved > 0)
            .count()
    };
    assert!(moved("connection") > 0, "no connection row moved");
    assert!(moved("") > moved("connection"), "no component row moved");
    eprintln!("worst: {worst:.3} of the allowance");
}
