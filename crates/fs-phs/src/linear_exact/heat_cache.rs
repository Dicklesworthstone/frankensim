//! Instance-owned reuse of the positive-heat factor of the exact ZOH step.
//!
//! `PositiveHeat::new` and `from_factor` read only the system's J and R and its resistance witness, the split, the
//! clock and the resistance channel being integrated. They never read the state or the held input. A coupled owner
//! that steps one fixed system at one fixed clock therefore rebuilds the same factor on every step: nine dense
//! exponentials and a QR each time. A [`HeatCache`] keeps the factors its holder has built. The key is a full copy of
//! everything the construction can read, and a lookup succeeds only when every input matches bit for bit
//! (`f64::to_bits`, so `-0.0` and `0.0` differ and a NaN equals itself). There is no hash: a lookup either finds an
//! entry whose every input is bit-equal or it rebuilds. The cache belongs to the instance that holds the fixed
//! system and is never global, so renders that share a thread or run in parallel cannot see each other's entries.
//!
//! The same discipline holds the held-input transition of the exact step: the d x d matrix that carries the state over
//! one step reads the system's J, R and G, the whole split, the clock and the force vector, and never the state. While
//! a struck bar rings the force is the constant rest force, so the transition is the same matrix on all but the few
//! hundred substeps of the contact, and every resistance channel of a step shares it. Its key is the same full,
//! bit-for-bit copy plus the force vector's bits, compared first; a changed force is a miss and rebuilds.
//!
//! The uncached entry points of [`super`] are unchanged and build both on every call. A holder opts in through the
//! `*_cached` functions, which return exactly the bits of their uncached twins.

use super::{LinearSplit, LinearSplitBudget, QuadraticPart, positive_heat::PositiveHeat};
use crate::{PhsError, PortHamiltonian, ResistanceWitness, StepRecord};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

/// Distinct keys one cache keeps. A coupled block needs four (its R, the component and connection channels and one
/// dashpot row). The oldest entry is dropped first; a dropped key is simply rebuilt.
const CAPACITY: usize = 8;

fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

fn same_scalar(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

// Exhaustive destructuring: a field added to `LinearSplit`, `QuadraticPart` or `LinearSplitBudget` stops this from
// compiling until it is compared here, so the key cannot silently miss an input.
fn same_split(a: &LinearSplit, b: &LinearSplit) -> bool {
    let LinearSplit {
        n, m, j, r, g, quadratic, values, basis, frequencies, frequency_basis, budget, pole_rate_bound,
    } = a;
    let LinearSplit {
        n: n2,
        m: m2,
        j: j2,
        r: r2,
        g: g2,
        quadratic: quadratic2,
        values: values2,
        basis: basis2,
        frequencies: frequencies2,
        frequency_basis: frequency_basis2,
        budget: budget2,
        pole_rate_bound: pole_rate_bound2,
    } = b;
    let QuadraticPart { matrix, exact } = quadratic;
    let QuadraticPart { matrix: matrix2, exact: exact2 } = quadratic2;
    let LinearSplitBudget { max_modes, max_omega_dt } = budget;
    let LinearSplitBudget { max_modes: max_modes2, max_omega_dt: max_omega_dt2 } = budget2;
    n == n2
        && m == m2
        && same(j, j2)
        && same(r, r2)
        && same(g, g2)
        && exact == exact2
        && same(matrix, matrix2)
        && same(values, values2)
        && same(basis, basis2)
        && same(frequencies, frequencies2)
        && same(frequency_basis, frequency_basis2)
        && max_modes == max_modes2
        && same_scalar(*max_omega_dt, *max_omega_dt2)
        && same_scalar(*pole_rate_bound, *pole_rate_bound2)
}

/// The numerical content of a [`PortHamiltonian`] that the heat construction can read. Its storage law is a trait
/// object that the construction never touches: the storage's quadratic part is already inside the split, and the
/// stepping entry points check it against the split before any factor is built.
struct SystemKey {
    n: usize,
    m: usize,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    witness: Option<(Vec<f64>, Vec<f64>)>,
}

impl SystemKey {
    fn of(sys: &PortHamiltonian) -> Self {
        let PortHamiltonian { n, m, j, r, g, storage: _, resistance_witness } = sys;
        Self {
            n: *n,
            m: *m,
            j: j.clone(),
            r: r.clone(),
            g: g.clone(),
            witness: resistance_witness
                .as_ref()
                .map(|w| (w.matrix.clone(), w.root.clone())),
        }
    }

    fn matches(&self, sys: &PortHamiltonian) -> bool {
        let PortHamiltonian { n, m, j, r, g, storage: _, resistance_witness } = sys;
        self.n == *n
            && self.m == *m
            && same(&self.j, j)
            && same(&self.r, r)
            && same(&self.g, g)
            && match (&self.witness, resistance_witness) {
                (None, None) => true,
                (Some((matrix, root)), Some(w)) => same(matrix, &w.matrix) && same(root, &w.root),
                _ => false,
            }
    }
}

/// Every input of one `PositiveHeat` construction.
struct HeatKey {
    dt_bits: u64,
    system: SystemKey,
    split: LinearSplit,
    resistance: Vec<f64>,
    /// The explicit factor of a witnessed channel (`from_factor`), absent for `new`.
    factor_root: Option<Vec<f64>>,
}

impl HeatKey {
    fn of(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        resistance: &[f64],
        witness: Option<&ResistanceWitness>,
    ) -> Self {
        Self {
            dt_bits: dt.to_bits(),
            system: SystemKey::of(sys),
            split: split.clone(),
            resistance: resistance.to_vec(),
            factor_root: witness.map(|w| w.root.clone()),
        }
    }

    fn matches(
        &self,
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        resistance: &[f64],
        witness: Option<&ResistanceWitness>,
    ) -> bool {
        self.dt_bits == dt.to_bits()
            && same(&self.resistance, resistance)
            && match (&self.factor_root, witness) {
                (None, None) => true,
                (Some(root), Some(w)) => same(root, &w.root),
                _ => false,
            }
            && self.system.matches(sys)
            && same_split(&self.split, split)
    }
}

/// Every input of one held-input transition (`held_input_transition`): the system's J, R and G, the whole split, the
/// clock and the force vector. The resistance channel being integrated and the resistance witness are NOT inputs:
/// the transition is shared by every channel's call of a step, and only the force vector changes it from step to step.
struct TransitionKey {
    dt_bits: u64,
    force_bits: Vec<u64>,
    system: SystemKey,
    split: LinearSplit,
}

impl TransitionKey {
    fn of(sys: &PortHamiltonian, split: &LinearSplit, dt: f64, force: &[f64]) -> Self {
        Self {
            dt_bits: dt.to_bits(),
            force_bits: force.iter().map(|v| v.to_bits()).collect(),
            system: SystemKey::of(sys),
            split: split.clone(),
        }
    }

    // The force is compared first: it is the one input that changes between steps, so a changed force misses on the
    // shortest comparison and an unchanged one pays for the full comparison it needs.
    fn matches(&self, sys: &PortHamiltonian, split: &LinearSplit, dt: f64, force: &[f64]) -> bool {
        self.force_bits.len() == force.len()
            && self.force_bits.iter().zip(force).all(|(bits, v)| *bits == v.to_bits())
            && self.dt_bits == dt.to_bits()
            && self.system.matches(sys)
            && same_split(&self.split, split)
    }
}

/// Positive-heat factors and held-input transitions built by one holder of a fixed system, reused only on a
/// bit-for-bit equal key. See the module documentation. A clone starts empty: the entries belong to the instance,
/// not to its copies.
pub struct HeatCache {
    entries: Mutex<Vec<(HeatKey, Arc<PositiveHeat>)>>,
    transitions: Mutex<Vec<(TransitionKey, Arc<Vec<f64>>)>>,
    builds: AtomicU64,
    reuses: AtomicU64,
    transition_builds: AtomicU64,
    transition_reuses: AtomicU64,
    enabled: bool,
}

impl HeatCache {
    /// An empty cache that reuses factors on an exactly equal key.
    #[must_use]
    pub fn new() -> Self {
        Self {
            entries: Mutex::new(Vec::new()),
            transitions: Mutex::new(Vec::new()),
            builds: AtomicU64::new(0),
            reuses: AtomicU64::new(0),
            transition_builds: AtomicU64::new(0),
            transition_reuses: AtomicU64::new(0),
            enabled: true,
        }
    }

    /// A cache that never stores or reuses anything: every call builds, exactly as the uncached functions do. It
    /// counts its builds, so a test can hold it beside an enabled twin and compare their outputs bit for bit.
    #[must_use]
    pub fn disabled() -> Self {
        Self { enabled: false, ..Self::new() }
    }

    /// Factors built through this cache (a miss, or every call of a disabled cache).
    #[must_use]
    pub fn builds(&self) -> u64 {
        self.builds.load(Ordering::Relaxed)
    }

    /// Factors returned from an entry whose key matched bit for bit.
    #[must_use]
    pub fn reuses(&self) -> u64 {
        self.reuses.load(Ordering::Relaxed)
    }

    /// Held-input transitions built through this cache (a miss, or every call of a disabled cache).
    #[must_use]
    pub fn transition_builds(&self) -> u64 {
        self.transition_builds.load(Ordering::Relaxed)
    }

    /// Held-input transitions returned from an entry whose key matched bit for bit.
    #[must_use]
    pub fn transition_reuses(&self) -> u64 {
        self.transition_reuses.load(Ordering::Relaxed)
    }

    /// Held-input transitions currently held.
    #[must_use]
    pub fn transition_len(&self) -> usize {
        self.transitions.lock().unwrap_or_else(PoisonError::into_inner).len()
    }

    /// Positive-heat factors currently held.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.lock().unwrap_or_else(PoisonError::into_inner).len()
    }

    /// True when no entry is held.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Default for HeatCache {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for HeatCache {
    fn clone(&self) -> Self {
        if self.enabled { Self::new() } else { Self::disabled() }
    }
}

impl core::fmt::Debug for HeatCache {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HeatCache")
            .field("enabled", &self.enabled)
            .field("entries", &self.len())
            .field("builds", &self.builds())
            .field("reuses", &self.reuses())
            .field("transitions", &self.transition_len())
            .field("transition_builds", &self.transition_builds())
            .field("transition_reuses", &self.transition_reuses())
            .finish()
    }
}

/// The factor for one construction: from `cache` on a bit-equal key, otherwise built exactly as the uncached path
/// builds it (and stored when the cache is enabled). A failed construction stores nothing.
pub(super) fn obtain(
    cache: Option<&HeatCache>,
    sys: &PortHamiltonian,
    split: &LinearSplit,
    dt: f64,
    resistance: &[f64],
    witness: Option<&ResistanceWitness>,
) -> Result<Arc<PositiveHeat>, PhsError> {
    let build = || -> Result<PositiveHeat, PhsError> {
        match witness {
            Some(w) => PositiveHeat::from_factor(sys, split, dt, resistance, w.root.clone()),
            None => PositiveHeat::new(sys, split, dt, resistance),
        }
    };
    let Some(cache) = cache else {
        return build().map(Arc::new);
    };
    if cache.enabled {
        let entries = cache.entries.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(entry) = entries
            .iter()
            .find(|entry| entry.0.matches(sys, split, dt, resistance, witness))
        {
            cache.reuses.fetch_add(1, Ordering::Relaxed);
            return Ok(Arc::clone(&entry.1));
        }
    }
    let heat = Arc::new(build()?);
    cache.builds.fetch_add(1, Ordering::Relaxed);
    if cache.enabled {
        let mut entries = cache.entries.lock().unwrap_or_else(PoisonError::into_inner);
        if entries.len() >= CAPACITY {
            entries.remove(0);
        }
        entries.push((
            HeatKey::of(sys, split, dt, resistance, witness),
            Arc::clone(&heat),
        ));
    }
    Ok(heat)
}

/// The held-input transition for one step: from `cache` on a bit-equal key (system, split, clock and force vector),
/// otherwise built by `build` exactly as the uncached path builds it (and stored when the cache is enabled). A failed
/// build stores nothing. `build` runs outside the lock.
pub(super) fn transition(
    cache: Option<&HeatCache>,
    sys: &PortHamiltonian,
    split: &LinearSplit,
    dt: f64,
    force: &[f64],
    build: impl FnOnce() -> Result<Vec<f64>, PhsError>,
) -> Result<Arc<Vec<f64>>, PhsError> {
    let Some(cache) = cache else {
        return build().map(Arc::new);
    };
    if cache.enabled {
        let entries = cache.transitions.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(entry) = entries
            .iter()
            .find(|entry| entry.0.matches(sys, split, dt, force))
        {
            cache.transition_reuses.fetch_add(1, Ordering::Relaxed);
            return Ok(Arc::clone(&entry.1));
        }
    }
    let transition = Arc::new(build()?);
    cache.transition_builds.fetch_add(1, Ordering::Relaxed);
    if cache.enabled {
        let mut entries = cache.transitions.lock().unwrap_or_else(PoisonError::into_inner);
        if entries.len() >= CAPACITY {
            entries.remove(0);
        }
        entries.push((
            TransitionKey::of(sys, split, dt, force),
            Arc::clone(&transition),
        ));
    }
    Ok(transition)
}

/// [`super::step_render`] that takes its positive-heat factor from `cache`. The returned record is bit-identical
/// to the uncached call's.
///
/// # Errors
/// Those of [`super::step_render`].
pub fn step_render_cached(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
    cache: &HeatCache,
) -> Result<StepRecord, PhsError> {
    super::step_render_with(sys, split, x0, u, dt, Some(cache))
}

/// [`super::integrated_resistance_loss`] that takes its positive-heat factor from `cache`. Bit-identical to the
/// uncached call.
///
/// # Errors
/// Those of [`super::integrated_resistance_loss`].
pub fn integrated_resistance_loss_cached(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
    resistance: &[f64],
    cache: &HeatCache,
) -> Result<f64, PhsError> {
    super::integrated_resistance_loss_inner(sys, split, x0, u, dt, resistance, None, Some(cache))
}

/// [`super::integrated_resistance_loss_with_witness`] that takes its positive-heat factor from `cache`.
/// Bit-identical to the uncached call.
///
/// # Errors
/// Those of [`super::integrated_resistance_loss_with_witness`].
pub fn integrated_resistance_loss_with_witness_cached(
    sys: &PortHamiltonian,
    split: &LinearSplit,
    x0: &[f64],
    u: &[f64],
    dt: f64,
    witness: &ResistanceWitness,
    cache: &HeatCache,
) -> Result<f64, PhsError> {
    super::integrated_resistance_loss_inner(
        sys,
        split,
        x0,
        u,
        dt,
        &witness.matrix,
        Some(witness),
        Some(cache),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_exact::{
        integrated_resistance_loss, integrated_resistance_loss_with_witness, split_quadratic,
        step_render,
    };
    use crate::mass_spring_damper;

    const DT: f64 = 0.05;
    // Exactly representable square roots keep the witnesses' factors exact: 0.25^2 = 0.0625 and 0.5^2 = 0.25.
    const DAMPING: f64 = 0.25;
    const CHANNEL: [f64; 4] = [0.0, 0.0, 0.0, 0.0625];

    type Fixture = (PortHamiltonian, LinearSplit, ResistanceWitness);

    fn split_of(sys: &PortHamiltonian) -> LinearSplit {
        split_quadratic(
            sys,
            LinearSplitBudget { max_modes: 2, max_omega_dt: core::f64::consts::PI },
        )
        .unwrap()
    }

    fn channel_witness() -> ResistanceWitness {
        ResistanceWitness::try_from_source_rows(2, CHANNEL.to_vec(), vec![0.0, 0.25]).unwrap()
    }

    fn fresh() -> Fixture {
        let sys_witness =
            ResistanceWitness::try_from_source_rows(2, vec![0.0, 0.0, 0.0, DAMPING], vec![0.0, 0.5])
                .unwrap();
        let sys = mass_spring_damper(1.3, 5.1, DAMPING)
            .unwrap()
            .with_resistance_witness(sys_witness)
            .unwrap();
        let split = split_of(&sys);
        (sys, split, channel_witness())
    }

    fn input(k: usize) -> [f64; 1] {
        [((k * 37 % 101) as f64) / 101.0 - 0.5]
    }

    fn bits(a: &[f64]) -> Vec<u64> {
        a.iter().map(|v| v.to_bits()).collect()
    }

    fn assert_same_record(a: &StepRecord, b: &StepRecord, k: usize) {
        assert_eq!(bits(&a.x), bits(&b.x), "x at step {k}");
        assert_eq!(bits(&a.y), bits(&b.y), "y at step {k}");
        assert_eq!(a.delta_h.to_bits(), b.delta_h.to_bits(), "delta_h at step {k}");
        assert_eq!(a.dissipated.to_bits(), b.dissipated.to_bits(), "dissipated at step {k}");
        assert_eq!(a.supplied.to_bits(), b.supplied.to_bits(), "supplied at step {k}");
        assert_eq!(a.newton_iters, b.newton_iters, "newton_iters at step {k}");
        assert_eq!(
            a.solver_residual.to_bits(),
            b.solver_residual.to_bits(),
            "solver_residual at step {k}"
        );
        assert!(a.solver_certificate.is_none() && b.solver_certificate.is_none());
    }

    #[test]
    fn the_cache_can_be_shared_across_threads() {
        fn send_sync<T: Send + Sync>() {}
        send_sync::<HeatCache>();
    }

    // G2: a long replay of the stepping and both loss integrals, with the state and the held force changing every
    // step, memoised against unmemoised. Every field of every record, and every loss, is compared with to_bits.
    #[test]
    fn cached_steps_and_loss_integrals_are_bit_identical_to_the_plain_functions() {
        let (sys, split, witness) = fresh();
        let cache = HeatCache::new();
        let steps = 120;
        let mut x = vec![0.3, -0.1];
        for k in 0..steps {
            let u = input(k);
            let plain = step_render(&sys, &split, &x, &u, DT).unwrap();
            let cached = step_render_cached(&sys, &split, &x, &u, DT, &cache).unwrap();
            assert_same_record(&plain, &cached, k);
            let a = integrated_resistance_loss(&sys, &split, &x, &u, DT, &CHANNEL).unwrap();
            let b = integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &CHANNEL, &cache)
                .unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "channel loss at step {k}");
            let c = integrated_resistance_loss_with_witness(&sys, &split, &x, &u, DT, &witness)
                .unwrap();
            let d = integrated_resistance_loss_with_witness_cached(
                &sys, &split, &x, &u, DT, &witness, &cache,
            )
            .unwrap();
            assert_eq!(c.to_bits(), d.to_bits(), "witnessed loss at step {k}");
            x = plain.x;
        }
        // Three distinct keys (the system's R, the plain channel, the witnessed channel), three calls a step.
        assert_eq!(cache.builds(), 3);
        assert_eq!(cache.reuses(), 3 * steps as u64 - 3);
        assert_eq!(cache.len(), 3);
        // The force changes every step here: one transition per step, shared by the step's three calls.
        assert_eq!(cache.transition_builds(), steps as u64);
        assert_eq!(cache.transition_reuses(), 2 * steps as u64);
    }

    // A channel value for which one ulp on its diagonal really moves the plain output. A value like 0.0625 = 2^-4 does
    // not: the square root that builds the factor absorbs the ulp, so a negative planted there could never fail.
    fn moving_channel(sys: &PortHamiltonian, split: &LinearSplit) -> ([f64; 4], [f64; 4]) {
        let (x, u) = ([0.3, -0.1], input(3));
        for k in 1..=64 {
            let value = 0.0625 + 0.0137 * f64::from(k);
            let base = [0.0, 0.0, 0.0, value];
            let mut bumped = base;
            bumped[3] = f64::from_bits(value.to_bits() + 1);
            let a = integrated_resistance_loss(sys, split, &x, &u, DT, &base).unwrap();
            let b = integrated_resistance_loss(sys, split, &x, &u, DT, &bumped).unwrap();
            if a.to_bits() != b.to_bits() {
                return (base, bumped);
            }
        }
        panic!("no channel value in the scan is moved by one ulp");
    }

    // The planted negative: the resistance changes by one ulp between steps. The cache must rebuild, and the output
    // must follow the plain path bit for bit, including through the change and back.
    #[test]
    fn planted_one_ulp_resistance_change_rebuilds_and_follows_the_plain_path() {
        let (sys, split, _) = fresh();
        let cache = HeatCache::new();
        let (channel, bumped) = moving_channel(&sys, &split);
        assert_ne!(bumped[3].to_bits(), channel[3].to_bits());
        let mut moved = 0;
        let mut x = vec![0.3, -0.1];
        let steps = 40;
        for k in 0..steps {
            let u = input(k);
            // Two steps on the channel, two on the bumped channel, and so on.
            let (r, other) = if (k / 2) % 2 == 0 { (&channel, &bumped) } else { (&bumped, &channel) };
            let plain = integrated_resistance_loss(&sys, &split, &x, &u, DT, r).unwrap();
            let moved_plain = integrated_resistance_loss(&sys, &split, &x, &u, DT, other).unwrap();
            let cached =
                integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, r, &cache).unwrap();
            assert_eq!(plain.to_bits(), cached.to_bits(), "loss at step {k}");
            if plain.to_bits() != moved_plain.to_bits() {
                moved += 1;
            }
            x = step_render(&sys, &split, &x, &u, DT).unwrap().x;
        }
        // The negative can fail: the one-ulp change really does move the plain output on most steps.
        assert!(moved > 0, "a one-ulp resistance change never moved the plain output");
        // Exactly two keys: the first bump rebuilds, and the return to the first channel reuses its entry.
        assert_eq!(cache.builds(), 2);
        assert_eq!(cache.reuses(), steps as u64 - 2);
    }

    // One planted change per input the construction can read.
    #[derive(Clone, Copy, Debug)]
    enum Plant {
        Dt,
        Resistance,
        FactorRoot,
        SysJ,
        SysJNegativeZero,
        SysR,
        SysG,
        SysWitnessMatrix,
        SysWitnessRoot,
        SplitJ,
        SplitR,
        SplitG,
        SplitQuadraticMatrix,
        SplitValues,
        SplitBasis,
        SplitFrequencies,
        SplitFrequencyBasis,
        SplitMaxOmegaDt,
        SplitPoleRateBound,
    }

    fn ulp(v: &mut f64) {
        *v = f64::from_bits(v.to_bits() + 1);
    }

    fn plant(which: Plant, f: &mut Fixture, dt: &mut f64, resistance: &mut [f64]) {
        match which {
            Plant::Dt => ulp(dt),
            Plant::Resistance => ulp(&mut resistance[3]),
            Plant::FactorRoot => ulp(&mut f.2.root[1]),
            Plant::SysJ => ulp(&mut f.0.j[1]),
            Plant::SysJNegativeZero => {
                assert_eq!(f.0.j[0].to_bits(), 0);
                f.0.j[0] = -0.0;
            }
            Plant::SysR => ulp(&mut f.0.r[3]),
            Plant::SysG => ulp(&mut f.0.g[1]),
            Plant::SysWitnessMatrix => ulp(&mut f.0.resistance_witness.as_mut().unwrap().matrix[3]),
            Plant::SysWitnessRoot => ulp(&mut f.0.resistance_witness.as_mut().unwrap().root[1]),
            Plant::SplitJ => ulp(&mut f.1.j[1]),
            Plant::SplitR => ulp(&mut f.1.r[3]),
            Plant::SplitG => ulp(&mut f.1.g[1]),
            Plant::SplitQuadraticMatrix => ulp(&mut f.1.quadratic.matrix[0]),
            Plant::SplitValues => ulp(&mut f.1.values[0]),
            Plant::SplitBasis => ulp(&mut f.1.basis[0]),
            Plant::SplitFrequencies => ulp(&mut f.1.frequencies[0]),
            Plant::SplitFrequencyBasis => ulp(&mut f.1.frequency_basis[0]),
            Plant::SplitMaxOmegaDt => ulp(&mut f.1.budget.max_omega_dt),
            Plant::SplitPoleRateBound => ulp(&mut f.1.pole_rate_bound),
        }
    }

    // Every input the construction can read is part of the key: change each by one ulp (or a zero's sign) and the
    // cache must rebuild; the original key must still be served afterwards.
    #[test]
    fn every_input_the_construction_reads_is_part_of_the_key() {
        let plants = [
            Plant::Dt,
            Plant::Resistance,
            Plant::FactorRoot,
            Plant::SysJ,
            Plant::SysJNegativeZero,
            Plant::SysR,
            Plant::SysG,
            Plant::SysWitnessMatrix,
            Plant::SysWitnessRoot,
            Plant::SplitJ,
            Plant::SplitR,
            Plant::SplitG,
            Plant::SplitQuadraticMatrix,
            Plant::SplitValues,
            Plant::SplitBasis,
            Plant::SplitFrequencies,
            Plant::SplitFrequencyBasis,
            Plant::SplitMaxOmegaDt,
            Plant::SplitPoleRateBound,
        ];
        for which in plants {
            let cache = HeatCache::new();
            let (sys, split, witness) = fresh();
            let resistance = witness.matrix.clone();
            obtain(Some(&cache), &sys, &split, DT, &resistance, Some(&witness)).unwrap();
            assert_eq!((cache.builds(), cache.reuses()), (1, 0), "{which:?}: first call");
            let mut planted = fresh();
            let mut dt = DT;
            let mut planted_resistance = resistance.clone();
            plant(which, &mut planted, &mut dt, &mut planted_resistance);
            obtain(
                Some(&cache),
                &planted.0,
                &planted.1,
                dt,
                &planted_resistance,
                Some(&planted.2),
            )
            .unwrap_or_else(|e| panic!("{which:?}: construction refused the planted input: {e:?}"));
            assert_eq!(
                (cache.builds(), cache.reuses()),
                (2, 0),
                "{which:?}: a changed input must rebuild, not reuse"
            );
            obtain(Some(&cache), &sys, &split, DT, &resistance, Some(&witness)).unwrap();
            assert_eq!(
                (cache.builds(), cache.reuses()),
                (2, 1),
                "{which:?}: the original key is still served"
            );
        }
    }

    #[test]
    fn the_oldest_entry_is_dropped_at_capacity_and_every_output_still_matches() {
        let (sys, split, _) = fresh();
        let cache = HeatCache::new();
        let distinct = CAPACITY + 3;
        let x = [0.3, -0.1];
        let u = input(5);
        for i in 0..distinct {
            let channel = [0.0, 0.0, 0.0, 0.01 * (i + 1) as f64];
            let plain = integrated_resistance_loss(&sys, &split, &x, &u, DT, &channel).unwrap();
            let cached =
                integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &channel, &cache)
                    .unwrap();
            assert_eq!(plain.to_bits(), cached.to_bits(), "channel {i}");
        }
        assert_eq!(cache.len(), CAPACITY);
        assert_eq!(cache.builds(), distinct as u64);
        // The newest key is held; the oldest was dropped and is rebuilt, still bit-identical.
        let newest = [0.0, 0.0, 0.0, 0.01 * distinct as f64];
        integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &newest, &cache).unwrap();
        assert_eq!(cache.reuses(), 1);
        let oldest = [0.0, 0.0, 0.0, 0.01];
        let plain = integrated_resistance_loss(&sys, &split, &x, &u, DT, &oldest).unwrap();
        let cached =
            integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &oldest, &cache).unwrap();
        assert_eq!(plain.to_bits(), cached.to_bits());
        assert_eq!(cache.builds(), distinct as u64 + 1);
    }

    #[test]
    fn a_disabled_cache_never_reuses_and_a_clone_starts_empty() {
        let (sys, split, _) = fresh();
        let x = [0.3, -0.1];
        let u = input(7);
        let disabled = HeatCache::disabled();
        for _ in 0..5 {
            let plain = step_render(&sys, &split, &x, &u, DT).unwrap();
            let cached = step_render_cached(&sys, &split, &x, &u, DT, &disabled).unwrap();
            assert_same_record(&plain, &cached, 0);
        }
        assert_eq!((disabled.builds(), disabled.reuses(), disabled.len()), (5, 0, 0));
        assert_eq!(
            (disabled.transition_builds(), disabled.transition_reuses(), disabled.transition_len()),
            (5, 0, 0)
        );
        let live = HeatCache::new();
        step_render_cached(&sys, &split, &x, &u, DT, &live).unwrap();
        assert_eq!(live.len(), 1);
        let copy = live.clone();
        assert_eq!((copy.len(), copy.builds(), copy.reuses()), (0, 0, 0));
        assert!(copy.is_empty());
        assert_eq!(disabled.clone().builds(), 0);
    }

    #[test]
    fn a_refused_input_returns_the_plain_error_and_stores_nothing() {
        let (sys, split, _) = fresh();
        let cache = HeatCache::new();
        let x = [0.3, -0.1];
        let u = input(1);
        let asymmetric = [0.0, 0.5, 0.0, 0.0625];
        let plain = integrated_resistance_loss(&sys, &split, &x, &u, DT, &asymmetric).unwrap_err();
        let cached =
            integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &asymmetric, &cache)
                .unwrap_err();
        assert_eq!(format!("{plain:?}"), format!("{cached:?}"));
        assert_eq!((cache.builds(), cache.reuses(), cache.len()), (0, 0, 0));
        assert_eq!((cache.transition_builds(), cache.transition_len()), (0, 0));
    }

    fn records_differ(a: &StepRecord, b: &StepRecord) -> bool {
        bits(&a.x) != bits(&b.x)
            || bits(&a.y) != bits(&b.y)
            || a.delta_h.to_bits() != b.delta_h.to_bits()
            || a.dissipated.to_bits() != b.dissipated.to_bits()
            || a.supplied.to_bits() != b.supplied.to_bits()
    }

    // G2: the replay a struck bar produces. The force changes on each of 256 steps (the contact: every step is a miss,
    // and the 8-entry cache evicts), then is held for 344 steps (one miss, then hits). The cached path must equal the
    // plain path field by field and loss by loss under to_bits on every step, miss or hit.
    #[test]
    fn the_transition_replay_is_bit_identical_through_force_changes_and_holds() {
        let (sys, split, witness) = fresh();
        let cache = HeatCache::new();
        let (changing, held) = (256_usize, 344_usize);
        let hold = [0.123_456_789];
        let mut x = vec![0.3, -0.1];
        for k in 0..changing + held {
            let u = if k < changing { input(k) } else { hold };
            let plain = step_render(&sys, &split, &x, &u, DT).unwrap();
            let cached = step_render_cached(&sys, &split, &x, &u, DT, &cache).unwrap();
            assert_same_record(&plain, &cached, k);
            let a = integrated_resistance_loss(&sys, &split, &x, &u, DT, &CHANNEL).unwrap();
            let b = integrated_resistance_loss_cached(&sys, &split, &x, &u, DT, &CHANNEL, &cache)
                .unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "channel loss at step {k}");
            let c = integrated_resistance_loss_with_witness(&sys, &split, &x, &u, DT, &witness)
                .unwrap();
            let d = integrated_resistance_loss_with_witness_cached(
                &sys, &split, &x, &u, DT, &witness, &cache,
            )
            .unwrap();
            assert_eq!(c.to_bits(), d.to_bits(), "witnessed loss at step {k}");
            x = plain.x;
        }
        // One build per changed force and one for the held force; the step's other two calls (three channels share
        // one transition) and every later held step reuse it.
        let calls = 3 * (changing + held) as u64;
        assert_eq!(cache.transition_builds(), changing as u64 + 1);
        assert_eq!(cache.transition_reuses(), calls - changing as u64 - 1);
        assert_eq!(cache.transition_len(), CAPACITY);
        // The positive-heat factors are untouched by the force: three keys for the whole run.
        assert_eq!(cache.builds(), 3);
    }

    // The planted negative for a stale force key: a force that differs by one ulp must change the plain output, so a
    // cache that ignored the force would be caught; the cached path alternates between the two forces and must follow
    // the plain path bit for bit, rebuilding exactly twice and reusing from then on.
    #[test]
    fn a_one_ulp_force_change_moves_the_output_and_the_cache_follows_it() {
        let (sys, split, _) = fresh();
        let cache = HeatCache::new();
        let x0 = [0.3, -0.1];
        let (force, bumped) = (|| {
            for k in 1..=64 {
                let value = 0.37 + 0.0137 * f64::from(k);
                let bumped = f64::from_bits(value.to_bits() + 1);
                let a = step_render(&sys, &split, &x0, &[value], DT).unwrap();
                let b = step_render(&sys, &split, &x0, &[bumped], DT).unwrap();
                if records_differ(&a, &b) {
                    return ([value], [bumped]);
                }
            }
            panic!("no force in the scan is moved by one ulp");
        })();
        assert_ne!(force[0].to_bits(), bumped[0].to_bits());
        let mut moved = 0;
        let mut x = x0.to_vec();
        let steps = 40;
        for k in 0..steps {
            // Two steps on the force, two on the bumped force, and so on.
            let (u, other) = if (k / 2) % 2 == 0 { (&force, &bumped) } else { (&bumped, &force) };
            let plain = step_render(&sys, &split, &x, u, DT).unwrap();
            let stale = step_render(&sys, &split, &x, other, DT).unwrap();
            let cached = step_render_cached(&sys, &split, &x, u, DT, &cache).unwrap();
            assert_same_record(&plain, &cached, k);
            if records_differ(&plain, &stale) {
                moved += 1;
            }
            x = plain.x;
        }
        // The negative can fail: serving the other force's transition would have moved the output on most steps.
        assert!(moved > 0, "a one-ulp force change never moved the plain output");
        assert_eq!(cache.transition_builds(), 2);
        assert_eq!(cache.transition_reuses(), steps as u64 - 2);
    }

    // Every input of the transition is part of its key, the force vector first: one ulp, a zero's sign or a different
    // length is a miss. The resistance channel and the factor root are not inputs (the transition is shared by every
    // channel of a step), so planting them must reuse.
    #[test]
    fn every_input_of_the_transition_is_part_of_its_key() {
        let dummy = || Ok(vec![0.0; 4]);
        let one = [0.25_f64];
        let ulp_up = [f64::from_bits(one[0].to_bits() + 1)];
        let ulp_down = [f64::from_bits(one[0].to_bits() - 1)];
        let forces: [(&str, &[f64], &[f64]); 5] = [
            ("one ulp up", &one, &ulp_up),
            ("one ulp down", &one, &ulp_down),
            ("sign of the force", &one, &[-0.25]),
            ("sign of a zero force", &[0.0], &[-0.0]),
            ("length of the force", &one, &[0.25, 0.0]),
        ];
        for (name, base, planted) in forces {
            let cache = HeatCache::new();
            let (sys, split, _) = fresh();
            transition(Some(&cache), &sys, &split, DT, base, dummy).unwrap();
            assert_eq!((cache.transition_builds(), cache.transition_reuses()), (1, 0), "{name}");
            transition(Some(&cache), &sys, &split, DT, planted, dummy).unwrap();
            assert_eq!(
                (cache.transition_builds(), cache.transition_reuses()),
                (2, 0),
                "{name}: a changed force must miss"
            );
            transition(Some(&cache), &sys, &split, DT, base, dummy).unwrap();
            assert_eq!(
                (cache.transition_builds(), cache.transition_reuses()),
                (2, 1),
                "{name}: the original force is still served"
            );
        }
        let plants = [
            Plant::Dt,
            Plant::Resistance,
            Plant::FactorRoot,
            Plant::SysJ,
            Plant::SysJNegativeZero,
            Plant::SysR,
            Plant::SysG,
            Plant::SysWitnessMatrix,
            Plant::SysWitnessRoot,
            Plant::SplitJ,
            Plant::SplitR,
            Plant::SplitG,
            Plant::SplitQuadraticMatrix,
            Plant::SplitValues,
            Plant::SplitBasis,
            Plant::SplitFrequencies,
            Plant::SplitFrequencyBasis,
            Plant::SplitMaxOmegaDt,
            Plant::SplitPoleRateBound,
        ];
        for which in plants {
            let cache = HeatCache::new();
            let (sys, split, _) = fresh();
            transition(Some(&cache), &sys, &split, DT, &one, dummy).unwrap();
            let mut planted = fresh();
            let mut dt = DT;
            let mut resistance = vec![0.0; 4];
            plant(which, &mut planted, &mut dt, &mut resistance);
            transition(Some(&cache), &planted.0, &planted.1, dt, &one, dummy).unwrap();
            let shared = matches!(which, Plant::Resistance | Plant::FactorRoot);
            assert_eq!(
                (cache.transition_builds(), cache.transition_reuses()),
                if shared { (1, 1) } else { (2, 0) },
                "{which:?}: {}",
                if shared { "not an input of the transition: must reuse" } else { "an input: must miss" }
            );
        }
    }
}
