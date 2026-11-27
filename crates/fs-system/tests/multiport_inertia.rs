//! fs-system's multi-port inertia, `elements::MultiportInertia`: one row and one port against
//! `Mass` to the bit.
//! Run with `cargo test -p fs-system --test multiport_inertia`.

// Exact float comparisons here assert bitwise agreement with `Mass`.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::{Component, Proposal, Receipt};
use fs_system::elements::{Mass, Motion, MultiportInertia};

const KIND: PortKind = PortKind::MechanicalForceVelocity;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 12,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Deterministic numbers in `[lo, hi)` (splitmix64), so a sweep is the same on every run.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        lo + (hi - lo) * ((z >> 11) as f64 / (1_u64 << 53) as f64)
    }
}

/// One step proposed directly: the new state, the outputs and the receipt.
fn step<C: Component>(
    element: &C,
    old: &C::State,
    dt_s: f64,
    inputs: &[f64],
) -> (C::State, Vec<f64>, Receipt)
where
    C::State: Clone,
{
    let (mut new, mut outputs) = (old.clone(), vec![0.0; inputs.len()]);
    let proposal = with_cx(|cx| element.propose(old, dt_s, inputs, &mut new, &mut outputs, cx));
    let Ok(Proposal::Step(receipt)) = proposal else {
        panic!("no step at {inputs:?}");
    };
    (new, outputs, receipt)
}

#[test]
fn one_row_and_one_port_is_the_mass_step_to_the_bit() {
    let mut numbers = Numbers(1_234);
    for _ in 0..400 {
        let (m, v0) = (numbers.next(0.01, 100.0), numbers.next(-5.0, 5.0));
        let (dt, force) = (numbers.next(1.0e-6, 0.1), numbers.next(-50.0, 50.0));
        let inertia = MultiportInertia::new(KIND, vec![m], vec![1.0], 1).expect("one row");
        let mass = Mass::new(m).expect("mass");
        let start = Motion {
            position_m: 0.0,
            velocity_m_s: v0,
        };
        let (motion, by_mass, receipt_mass) = step(&mass, &start, dt, &[force]);
        let (velocity, by_inertia, receipt) = step(&inertia, &vec![v0], dt, &[force]);
        assert_eq!(velocity[0].to_bits(), motion.velocity_m_s.to_bits());
        assert_eq!(by_inertia[0].to_bits(), by_mass[0].to_bits());
        let bits = |r: Receipt| [r.stored_change_j.to_bits(), r.energy_tolerance_j.to_bits()];
        assert_eq!(bits(receipt), bits(receipt_mass), "receipt");
        assert_eq!(receipt, receipt_mass);
        let energy = inertia.stored_energy_j(&velocity).expect("energy");
        let by_mass = mass.stored_energy_j(&motion).expect("energy");
        assert_eq!(energy.to_bits(), by_mass.to_bits(), "stored energy");
    }
}
