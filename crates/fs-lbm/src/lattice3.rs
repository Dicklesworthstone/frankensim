//! Checked lattice-generic dense reference grid with fixed64-cell tiles.
//! D3Q27 and the unmodified D3Q19 cell owner share this transport/scheduling
//! seam. Existing specialized D3Q19 sparse, boundary and SIMD paths are not
//! redirected, so their traversal/arithmetic contracts remain untouched.
//! Three population buffers buy failure-atomic publication after cancellation.

use core::{fmt, marker::PhantomData, ops::ControlFlow};
use fs_exec::{CancelGate, Cancelled, Cx, KernelRunner, Reduce, TileKernel, TilePlan};
use std::sync::Mutex;

/// A fixed population lattice and its caller-selected collision law.
pub trait Lattice3<const Q: usize>: Send + Sync + 'static {
    type Collision: Copy + Send + Sync;
    type Error: Clone + fmt::Debug + fmt::Display + Send + Sync;
    const E: [[i32; 3]; Q];
    const W: [f64; Q];
    const OPP: [usize; Q];
    const CS2: f64;
    fn validate(collision: Self::Collision) -> Result<(), Self::Error>;
    fn collide(
        f: [f64; Q],
        collision: Self::Collision,
        force: [f64; 3],
    ) -> Result<[f64; Q], Self::Error>;
}

/// Compatibility adapter calls the existing D3Q19 kernel, with no new
/// default collision model and no change to the historical grids/goldens.
pub struct D3q19;
impl Lattice3<19> for D3q19 {
    type Collision = crate::d3q19::CollisionModel3;
    type Error = crate::d3q19::CollisionError3;
    const E: [[i32; 3]; 19] = d19_directions();
    const W: [f64; 19] = crate::d3q19::W3;
    const OPP: [usize; 19] = crate::d3q19::OPP3;
    const CS2: f64 = 1.0 / 3.0;
    fn validate(c: Self::Collision) -> Result<(), Self::Error> {
        c.validate()
    }
    fn collide(f: [f64; 19], c: Self::Collision, g: [f64; 3]) -> Result<[f64; 19], Self::Error> {
        crate::d3q19::collide_cell3(f, c, g)
    }
}
const fn d19_directions() -> [[i32; 3]; 19] {
    let mut out = [[0; 3]; 19];
    let mut i = 0;
    while i < 19 {
        let e = crate::d3q19::E3[i];
        out[i] = [e.0, e.1, e.2];
        i += 1;
    }
    out
}

/// Caller-owned resource and weak-compressibility envelope.
#[derive(Clone, Copy, Debug)]
pub struct GridBudget3 {
    pub maximum_cells: usize,
    pub maximum_mach: f64,
}

/// Grid admission/step refusal. A failed step leaves published populations
/// and completed_steps unchanged; scratch may contain an unpublished prefix.
#[derive(Clone, Debug)]
pub enum GridError3<E> {
    Input(&'static str),
    Cap {
        requested: usize,
        maximum: usize,
    },
    Allocation,
    Collision {
        cell: usize,
        source: E,
    },
    NonFiniteCell {
        cell: usize,
    },
    Mach {
        cell: usize,
        observed: f64,
        maximum: f64,
    },
    Cancelled,
    Pool(String),
    Poisoned,
}
impl<E: fmt::Debug> fmt::Display for GridError3<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-LBM-GRID3: {self:?}")
    }
}
impl<E: fmt::Debug> std::error::Error for GridError3<E> {}

/// Immutable topology, cell-major populations and transactional scratch.
pub struct Grid3<L: Lattice3<Q>, const Q: usize> {
    dims: [usize; 3],
    periodic: [bool; 3],
    solid: Vec<bool>,
    collision: L::Collision,
    force: [f64; 3],
    budget: GridBudget3,
    f: Vec<[f64; Q]>,
    post: Vec<[f64; Q]>,
    next: Vec<[f64; Q]>,
    completed_steps: usize,
    marker: PhantomData<L>,
}

fn macro_cell<L: Lattice3<Q>, const Q: usize>(
    f: &[f64; Q],
    force: [f64; 3],
    cell: usize,
    maximum: f64,
) -> Result<(f64, [f64; 3]), GridError3<L::Error>> {
    let rho: f64 = f.iter().sum();
    if !rho.is_finite() || rho <= 0.0 || f.iter().any(|v| !v.is_finite()) {
        return Err(GridError3::NonFiniteCell { cell });
    }
    let u: [f64; 3] = core::array::from_fn(|a| {
        (f.iter()
            .enumerate()
            .map(|(q, p)| f64::from(L::E[q][a]) * p)
            .sum::<f64>()
            + 0.5 * force[a])
            / rho
    });
    let mach = (u.iter().map(|v| v * v).sum::<f64>() / L::CS2).sqrt();
    if !mach.is_finite() || mach > maximum {
        return Err(GridError3::Mach {
            cell,
            observed: mach,
            maximum,
        });
    }
    Ok((rho, u))
}

impl<L: Lattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Every physical choice is explicit. Nonperiodic external links and
    /// supplied solid cells use stationary halfway bounce-back. No open,
    /// moving, interpolated or guessed wall is introduced by this constructor.
    pub fn try_new(
        dims: [usize; 3],
        collision: L::Collision,
        force: [f64; 3],
        periodic: [bool; 3],
        budget: GridBudget3,
        solid: impl Fn([usize; 3]) -> bool,
        initial: impl Fn([usize; 3]) -> [f64; Q],
    ) -> Result<Self, GridError3<L::Error>> {
        L::validate(collision).map_err(|source| GridError3::Collision { cell: 0, source })?;
        if Q == 0
            || dims.contains(&0)
            || dims.iter().any(|&d| d > isize::MAX as usize)
            || force.iter().any(|f| !f.is_finite())
            || !budget.maximum_mach.is_finite()
            || budget.maximum_mach <= 0.0
            || budget.maximum_mach > crate::MACH_LIMIT
        {
            return Err(GridError3::Input("dimensions, force or Mach envelope"));
        }
        if !(L::CS2.is_finite() && L::CS2 > 0.0)
            || (0..Q).any(|q| {
                L::OPP[q] >= Q
                    || L::E[q].iter().any(|c| !(-1..=1).contains(c))
                    || !L::W[q].is_finite()
                    || L::W[q] <= 0.0
            })
        {
            return Err(GridError3::Input("lattice tables"));
        }
        for q in 0..Q {
            let o = L::OPP[q];
            if L::OPP[o] != q || L::W[q] != L::W[o] || (0..3).any(|a| L::E[o][a] != -L::E[q][a]) {
                return Err(GridError3::Input("lattice opposite involution"));
            }
        }
        let count = dims[0]
            .checked_mul(dims[1])
            .and_then(|v| v.checked_mul(dims[2]))
            .ok_or(GridError3::Input("cell-count overflow"))?;
        if count > budget.maximum_cells {
            return Err(GridError3::Cap {
                requested: count,
                maximum: budget.maximum_cells,
            });
        }
        count
            .checked_mul(Q)
            .and_then(|n| n.checked_mul(24))
            .ok_or(GridError3::Input("three-buffer bytes overflow"))?;
        let mut f = Vec::new();
        f.try_reserve_exact(count)
            .map_err(|_| GridError3::Allocation)?;
        let mut mask = Vec::new();
        mask.try_reserve_exact(count)
            .map_err(|_| GridError3::Allocation)?;
        for cell in 0..count {
            let p = coordinates(cell, dims);
            let populations = initial(p);
            macro_cell::<L, Q>(&populations, force, cell, budget.maximum_mach)?;
            f.push(populations);
            mask.push(solid(p));
        }
        let mut post = Vec::new();
        post.try_reserve_exact(count)
            .map_err(|_| GridError3::Allocation)?;
        post.resize(count, [0.0; Q]);
        let mut next = Vec::new();
        next.try_reserve_exact(count)
            .map_err(|_| GridError3::Allocation)?;
        next.resize(count, [0.0; Q]);
        Ok(Self {
            dims,
            periodic,
            solid: mask,
            collision,
            force,
            budget,
            f,
            post,
            next,
            completed_steps: 0,
            marker: PhantomData,
        })
    }
    #[must_use]
    pub fn dimensions(&self) -> [usize; 3] {
        self.dims
    }
    #[must_use]
    pub fn populations(&self) -> &[[f64; Q]] {
        &self.f
    }
    #[must_use]
    pub fn completed_steps(&self) -> usize {
        self.completed_steps
    }
    pub fn macroscopic(&self, cell: usize) -> Result<(f64, [f64; 3]), GridError3<L::Error>> {
        let f = self.f.get(cell).ok_or(GridError3::Input("cell index"))?;
        macro_cell::<L, Q>(f, self.force, cell, self.budget.maximum_mach)
    }
    /// Serial reference, same64-cell partition and same arithmetic as pooled.
    pub fn step_serial(&mut self, cx: &Cx<'_>) -> Result<(), GridError3<L::Error>> {
        self.step_serial_impl(cx, None)
    }
    /// Advance with one caller-owned collision law per cell, in x-fastest order.
    /// The full field (including solid entries) is admitted before any work.
    /// It applies only to this step; a later `step_serial` uses the constructor's
    /// uniform law. Material selection, units and coefficient epochs belong to
    /// the caller. Spatially varying rates do not supply an energy equation.
    /// Errors/cancellation leave populations and completed_steps unchanged.
    pub fn step_serial_with_collisions(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
    ) -> Result<(), GridError3<L::Error>> {
        self.admit_collisions(collisions, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.step_serial_impl(cx, Some(collisions))
    }
    fn step_serial_impl(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
    ) -> Result<(), GridError3<L::Error>> {
        let data = StepData::<L, Q> {
            dims: self.dims,
            periodic: self.periodic,
            solid: &self.solid,
            collision: self.collision,
            collisions,
            force: self.force,
            budget: self.budget,
            marker: PhantomData,
        };
        for (tile, chunk) in self.post.chunks_mut(64).enumerate() {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
            data.process(&self.f, chunk, tile * 64, false)?;
        }
        for (tile, chunk) in self.next.chunks_mut(64).enumerate() {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
            data.process(&self.post, chunk, tile * 64, true)?;
        }
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        std::mem::swap(&mut self.f, &mut self.next);
        self.completed_steps += 1;
        Ok(())
    }
    /// Real executor seam. The fixed partition never depends on worker count;
    /// each mutable chunk has one logical owner. Neither pass publishes early.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.step_pooled_impl(pool, gate, None)
    }
    /// Executor-backed counterpart of `Self::step_serial_with_collisions`.
    /// The borrowed field is immutable across both passes and all workers;
    /// partitioning, collision arithmetic and publication use the same owner.
    pub fn step_pooled_with_collisions<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
    ) -> Result<(), GridError3<L::Error>> {
        self.admit_collisions(collisions, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.step_pooled_impl(pool, gate, Some(collisions))
    }
    fn admit_collisions(
        &self,
        collisions: &[L::Collision],
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(), GridError3<L::Error>> {
        checkpoint()?;
        if collisions.len() != self.f.len() {
            return Err(GridError3::Input("collision field length"));
        }
        for (tile, chunk) in collisions.chunks(64).enumerate() {
            checkpoint()?;
            for (offset, &collision) in chunk.iter().enumerate() {
                L::validate(collision).map_err(|source| GridError3::Collision {
                    cell: tile * 64 + offset,
                    source,
                })?;
            }
        }
        checkpoint()
    }
    fn step_pooled_impl<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
    ) -> Result<(), GridError3<L::Error>> {
        let data = StepData::<L, Q> {
            dims: self.dims,
            periodic: self.periodic,
            solid: &self.solid,
            collision: self.collision,
            collisions,
            force: self.force,
            budget: self.budget,
            marker: PhantomData,
        };
        run_pass(&data, &self.f, &mut self.post, false, pool, gate)?;
        run_pass(&data, &self.post, &mut self.next, true, pool, gate)?;
        if gate.is_requested() {
            return Err(GridError3::Cancelled);
        }
        std::mem::swap(&mut self.f, &mut self.next);
        self.completed_steps += 1;
        Ok(())
    }
}

fn coordinates(cell: usize, d: [usize; 3]) -> [usize; 3] {
    [cell % d[0], (cell / d[0]) % d[1], cell / (d[0] * d[1])]
}

struct StepData<'a, L: Lattice3<Q>, const Q: usize> {
    dims: [usize; 3],
    periodic: [bool; 3],
    solid: &'a [bool],
    collision: L::Collision,
    collisions: Option<&'a [L::Collision]>,
    force: [f64; 3],
    budget: GridBudget3,
    marker: PhantomData<L>,
}
impl<L: Lattice3<Q>, const Q: usize> StepData<'_, L, Q> {
    fn process(
        &self,
        src: &[[f64; Q]],
        dst: &mut [[f64; Q]],
        base: usize,
        stream: bool,
    ) -> Result<(), GridError3<L::Error>> {
        for (offset, output) in dst.iter_mut().enumerate() {
            let cell = base + offset;
            if self.solid[cell] {
                *output = src[cell];
                continue;
            }
            if !stream {
                let collision = self.collisions.map_or(self.collision, |field| field[cell]);
                *output = L::collide(src[cell], collision, self.force)
                    .map_err(|source| GridError3::Collision { cell, source })?;
                macro_cell::<L, Q>(output, self.force, cell, self.budget.maximum_mach)?;
            } else {
                let p = coordinates(cell, self.dims);
                for q in 0..Q {
                    let mut from = [0; 3];
                    let mut wall = false;
                    for a in 0..3 {
                        let v = p[a] as isize - L::E[q][a] as isize;
                        if v < 0 || v >= self.dims[a] as isize {
                            if !self.periodic[a] {
                                wall = true;
                                break;
                            }
                            from[a] = v.rem_euclid(self.dims[a] as isize) as usize;
                        } else {
                            from[a] = v as usize;
                        }
                    }
                    let source = (from[2] * self.dims[1] + from[1]) * self.dims[0] + from[0];
                    output[q] = if wall || self.solid[source] {
                        src[cell][L::OPP[q]]
                    } else {
                        src[source][q]
                    };
                }
                macro_cell::<L, Q>(output, self.force, cell, self.budget.maximum_mach)?;
            }
        }
        Ok(())
    }
}

struct FirstError<E>(Option<GridError3<E>>);
impl<E: Send> Reduce for FirstError<E> {
    fn identity() -> Self {
        Self(None)
    }
    fn merge(self, other: Self) -> Self {
        if self.0.is_some() { self } else { other }
    }
}
struct StepKernel<'a, 'b, L: Lattice3<Q>, const Q: usize> {
    data: &'a StepData<'b, L, Q>,
    src: &'a [[f64; Q]],
    chunks: Vec<Mutex<&'a mut [[f64; Q]]>>,
    stream: bool,
}
impl<L: Lattice3<Q>, const Q: usize> TileKernel for StepKernel<'_, '_, L, Q> {
    type Out = FirstError<L::Error>;
    fn tiles(&self) -> TilePlan {
        TilePlan::new(
            if self.stream {
                "fs-lbm/grid3-stream"
            } else {
                "fs-lbm/grid3-collide"
            },
            self.chunks.len() as u64,
        )
    }
    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let outcome = match self.chunks[tile as usize].lock() {
            Ok(mut chunk) => {
                self.data
                    .process(self.src, &mut chunk, tile as usize * 64, self.stream)
            }
            Err(_) => Err(GridError3::Poisoned),
        };
        ControlFlow::Continue(FirstError(outcome.err()))
    }
}
fn run_pass<L: Lattice3<Q>, const Q: usize, P: KernelRunner>(
    data: &StepData<'_, L, Q>,
    src: &[[f64; Q]],
    dst: &mut [[f64; Q]],
    stream: bool,
    pool: &P,
    gate: &CancelGate,
) -> Result<(), GridError3<L::Error>> {
    let kernel = StepKernel {
        data,
        src,
        chunks: dst.chunks_mut(64).map(Mutex::new).collect(),
        stream,
    };
    match pool.run_with_gate(&kernel, gate).0 {
        Ok(FirstError(None)) => Ok(()),
        Ok(FirstError(Some(e))) => Err(e),
        Err(_) if gate.is_requested() => Err(GridError3::Cancelled),
        Err(e) => Err(GridError3::Pool(format!("{e:?}"))),
    }
}

/// Physical-to-lattice acoustic scaling; no relaxation or speed is defaulted.
#[derive(Clone, Copy, Debug)]
pub struct AcousticScalingPlan {
    pub dx_m: f64,
    pub dt_s: f64,
    pub nu_m2_s: f64,
    pub sound_speed_m_s: f64,
    pub lattice: crate::ScalingPlan,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ScalingError3 {
    Input(&'static str),
    Mach { observed: f64, maximum: f64 },
    Unrepresentable,
}
impl fmt::Display for ScalingError3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-LBM-ACOUSTIC-SCALING: {self:?}")
    }
}
impl std::error::Error for ScalingError3 {}

pub fn acoustic_scaling(
    nu_m2_s: f64,
    c_m_s: f64,
    dx_m: f64,
    maximum_speed_m_s: f64,
    maximum_mach: f64,
) -> Result<AcousticScalingPlan, ScalingError3> {
    for (v, field) in [
        (nu_m2_s, "nu_m2_s"),
        (c_m_s, "sound_speed_m_s"),
        (dx_m, "dx_m"),
        (maximum_mach, "maximum_mach"),
    ] {
        if !v.is_finite() || v <= 0.0 {
            return Err(ScalingError3::Input(field));
        }
    }
    if !maximum_speed_m_s.is_finite() || maximum_speed_m_s < 0.0 || maximum_mach > crate::MACH_LIMIT
    {
        return Err(ScalingError3::Input("speed/Mach envelope"));
    }
    let mach = maximum_speed_m_s / c_m_s;
    if mach > maximum_mach {
        return Err(ScalingError3::Mach {
            observed: mach,
            maximum: maximum_mach,
        });
    }
    let dt_s = crate::CS2.sqrt() * dx_m / c_m_s;
    let viscosity = nu_m2_s * dt_s / (dx_m * dx_m);
    let tau_margin = viscosity / crate::CS2;
    let tau = 0.5 + tau_margin;
    let u_lattice = maximum_speed_m_s * dt_s / dx_m;
    if !dt_s.is_finite() || dt_s <= 0.0 || !tau.is_finite() || tau <= 0.5 || !u_lattice.is_finite()
    {
        return Err(ScalingError3::Unrepresentable);
    }
    Ok(AcousticScalingPlan {
        dx_m,
        dt_s,
        nu_m2_s,
        sound_speed_m_s: c_m_s,
        lattice: crate::ScalingPlan {
            tau,
            viscosity,
            u_lattice,
            mach,
            tau_margin,
            stable: mach < crate::MACH_LIMIT,
        },
    })
}

impl AcousticScalingPlan {
    /// Derive nu=mu/rho and c from the caller's material gas state.
    pub fn from_gas(
        gas: &fs_material::gas::GasState,
        dx_m: f64,
        maximum_speed_m_s: f64,
        maximum_mach: f64,
    ) -> Result<Self, ScalingError3> {
        if !gas.density.is_finite() || gas.density <= 0.0 {
            return Err(ScalingError3::Input("gas density"));
        }
        acoustic_scaling(
            gas.dynamic_viscosity / gas.density,
            gas.sound_speed,
            dx_m,
            maximum_speed_m_s,
            maximum_mach,
        )
    }
}

#[cfg(test)]
mod collision_field_tests {
    use super::*;
    use crate::d3q19::CollisionModel3;
    use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
    use fs_exec::{RunError, RunReport, TilePool};
    use std::sync::atomic::{AtomicUsize, Ordering};

    const DIMS: [usize; 3] = [9, 5, 3];
    const CELLS: usize = 135; // two complete tiles plus a ragged tile
    const SEED: u64 = 0xC011_1D3;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        arenas.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey { seed: SEED, kernel_id: 3, tile: 0, iteration: 0 },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            f(&cx)
        })
    }
    fn make<L: Lattice3<Q>, const Q: usize>(collision: L::Collision) -> Grid3<L, Q> {
        Grid3::try_new(
            DIMS,
            collision,
            [1e-7, 0.0, 0.0],
            [true, false, true],
            GridBudget3 { maximum_cells: CELLS, maximum_mach: 0.2 },
            |p| p == [4, 2, 1],
            |p| {
                let u = 0.001 * (p[1] + 1) as f64;
                let mut f = core::array::from_fn(|q| {
                    let eu = f64::from(L::E[q][0]) * u;
                    L::W[q] * (1.0 + eu / L::CS2
                        + 0.5 * eu * eu / (L::CS2 * L::CS2) - 0.5 * u * u / L::CS2)
                });
                // Non-equilibrium stress makes the very first collision sensitive
                // to the rate, without changing density or momentum.
                let stress = 1e-5 * (p[0] + 1) as f64;
                f[0] -= 2.0 * stress;
                f[1] += stress;
                f[2] += stress;
                f
            },
        ).unwrap()
    }
    fn bits<L: Lattice3<Q>, const Q: usize>(g: &Grid3<L, Q>) -> Vec<u64> {
        g.populations().iter().flatten().map(|v| v.to_bits()).collect()
    }
    fn bgk27(tau: f64) -> CollisionConfig27 {
        CollisionConfig27 {
            model: CollisionModel27::Bgk { tau, equilibrium: Equilibrium27::Quadratic },
            maximum_mach: 0.2,
        }
    }
    fn constant_parity<L: Lattice3<Q>, const Q: usize>(collision: L::Collision) {
        let mut uniform = make::<L, Q>(collision);
        let mut field = make::<L, Q>(collision);
        let laws = vec![collision; CELLS];
        with_cx(|cx| {
            for _ in 0..5 {
                uniform.step_serial(cx).unwrap();
                field.step_serial_with_collisions(cx, &laws).unwrap();
                assert_eq!(bits(&uniform), bits(&field));
            }
        });
        assert_eq!(field.completed_steps(), 5);
    }
    #[test]
    fn constant_fields_preserve_d19_and_all_d27_collision_paths() {
        constant_parity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
        constant_parity::<D3q27, 27>(bgk27(0.8));
        let rates = Rates27 { omega: [1.25, 1.0, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5] };
        for model in [CollisionModel27::CentralMoment { rates }, CollisionModel27::Cumulant { rates }] {
            constant_parity::<D3q27, 27>(CollisionConfig27 { model, maximum_mach: 0.2 });
        }
    }
    #[test]
    fn heterogeneous_field_matches_independent_push_stencil_and_does_not_latch() {
        let mut grid = make::<D3q27, 27>(bgk27(0.8));
        let laws: Vec<_> = (0..CELLS).map(|i| bgk27(if i % 3 == 0 { 0.7 } else { 1.2 })).collect();
        let mut expected = grid.f.clone();
        // Independent PUSH traversal, not StepData's pull/index implementation.
        for i in 0..CELLS {
            if grid.solid[i] { continue; }
            let post = D3q27::collide(grid.f[i], laws[i], grid.force).unwrap();
            let p = [i % 9, (i / 9) % 5, i / 45];
            for q in 0..27 {
                let e = D3q27::E[q];
                let y = p[1] as i32 + e[1];
                if !(0..5).contains(&y) {
                    expected[i][D3q27::OPP[q]] = post[q];
                    continue;
                }
                let x = (p[0] as i32 + e[0]).rem_euclid(9) as usize;
                let z = (p[2] as i32 + e[2]).rem_euclid(3) as usize;
                let target = (z * 5 + y as usize) * 9 + x;
                if grid.solid[target] {
                    expected[i][D3q27::OPP[q]] = post[q];
                } else {
                    expected[target][q] = post[q];
                }
            }
        }
        let before_mass: f64 = grid.f.iter().flatten().sum();
        with_cx(|cx| grid.step_serial_with_collisions(cx, &laws).unwrap());
        let expected_bits: Vec<_> = expected.iter().flatten().map(|v| v.to_bits()).collect();
        assert_eq!(bits(&grid), expected_bits);
        let after_mass: f64 = grid.f.iter().flatten().sum();
        assert!((after_mass - before_mass).abs() < 1e-10);
        let mut uniform = make::<D3q27, 27>(bgk27(0.8));
        with_cx(|cx| uniform.step_serial(cx).unwrap());
        assert_ne!(bits(&grid), bits(&uniform));
        uniform.f.clone_from(&grid.f);
        with_cx(|cx| {
            grid.step_serial(cx).unwrap();
            uniform.step_serial_with_collisions(cx, &vec![bgk27(0.8); CELLS]).unwrap();
        });
        assert_eq!(bits(&grid), bits(&uniform));
    }
    #[test]
    fn malformed_fields_refuse_before_publication_including_solid_entries() {
        let mut grid = make::<D3q27, 27>(bgk27(0.8));
        let before = bits(&grid);
        with_cx(|cx| {
            for n in [0, CELLS - 1, CELLS + 1] {
                assert!(matches!(
                    grid.step_serial_with_collisions(cx, &vec![bgk27(0.8); n]),
                    Err(GridError3::Input("collision field length"))
                ));
            }
            for bad in [CELLS - 1, 67] { // late fluid entry, then the solid cell
                let mut laws = vec![bgk27(0.8); CELLS];
                laws[bad] = bgk27(f64::NAN);
                assert!(matches!(grid.step_serial_with_collisions(cx, &laws),
                    Err(GridError3::Collision { cell, .. }) if cell == bad));
            }
        });
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);
    }
    #[test]
    fn heterogeneous_pooled_steps_match_serial_for_ragged_tiles() {
        let laws: Vec<_> = (0..CELLS).map(|i| bgk27(0.7 + 0.05 * (i % 7) as f64)).collect();
        let mut serial = make::<D3q27, 27>(bgk27(0.8));
        with_cx(|cx| {
            for _ in 0..5 { serial.step_serial_with_collisions(cx, &laws).unwrap(); }
        });
        for workers in [1, 2, 4] {
            let pool = TilePool::for_host(workers, SEED);
            let gate = CancelGate::new_clock_free();
            let mut parallel = make::<D3q27, 27>(bgk27(0.8));
            for _ in 0..5 {
                parallel.step_pooled_with_collisions(&pool, &gate, &laws).unwrap();
            }
            assert_eq!(bits(&serial), bits(&parallel), "workers={workers}");
        }
    }
    #[test]
    fn field_cancellation_between_passes_rolls_back_and_retries() {
        struct CancelAfterCollision { pool: TilePool, calls: AtomicUsize }
        impl KernelRunner for CancelAfterCollision {
            fn workers(&self) -> usize { self.pool.workers() }
            fn run_with_gate<K: TileKernel>(&self, kernel: &K, gate: &CancelGate)
                -> (Result<K::Out, RunError>, RunReport) {
                let result = self.pool.run_with_gate(kernel, gate);
                if self.calls.fetch_add(1, Ordering::AcqRel) == 0 { gate.request(); }
                result
            }
        }
        let runner = CancelAfterCollision {
            pool: TilePool::for_host(2, SEED), calls: AtomicUsize::new(0),
        };
        let mut grid = make::<D3q27, 27>(bgk27(0.8));
        let before = bits(&grid);
        let laws: Vec<_> = (0..CELLS).map(|i| bgk27(0.7 + 0.1 * (i % 4) as f64)).collect();
        let cancelled = CancelGate::new_clock_free();
        assert!(matches!(grid.step_pooled_with_collisions(&runner, &cancelled, &laws),
            Err(GridError3::Cancelled)));
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);
        let fresh = CancelGate::new_clock_free();
        grid.step_pooled_with_collisions(&runner.pool, &fresh, &laws).unwrap();
        let mut reference = make::<D3q27, 27>(bgk27(0.8));
        with_cx(|cx| reference.step_serial_with_collisions(cx, &laws).unwrap());
        assert_eq!(bits(&grid), bits(&reference));
        assert_eq!(grid.completed_steps(), 1);
    }
}
