//! Fixed-site compact wall diffusion on the original physical neighbor graph.
//!
//! In an admitted half-plane, image distance squared is real distance squared
//! plus 4 h_i h_j. Every supported image pair is therefore already a real
//! neighbor; no reflected graph or radius approximation is needed. Self-image
//! loss remains present even when the physical graph has no edges.

use super::{PreparedPse, admit_count, count, payload};
use crate::diffusion::wall::{CachedLoss, pair_rates};
use crate::diffusion::{PairKernel, exchange_fraction, validate_inputs};
use crate::planar_wall::PlanarWall2;
use crate::{
    AdmittedBudget, Cx, VortexParticle, VpmBudget, VpmError, VpmRun, WorkMeter, admit_cap,
    invalid_scalar, reserved_vec, validate_particle,
};

struct PairMap {
    exchange: f64,
    losses: [CachedLoss; 2],
}

impl PairMap {
    fn loss(&self, state: &mut [VortexParticle], [i, j]: [usize; 2]) -> Result<(), VpmError> {
        state[i].circulation = self.losses[0].apply(state[i].circulation, i)?;
        state[j].circulation = self.losses[1].apply(state[j].circulation, j)?;
        Ok(())
    }
}

impl PreparedPse {
    /// Diffuse fixed physical sites beside a stationary planar stress-free wall.
    ///
    /// This is the exact compact law and exchange/loss composition of
    /// [crate::planar_wall::diffuse_with_cx] on the existing prepared graph.
    /// Image interactions are a subset of its real neighbors; no image cloud,
    /// cutoff approximation or second search is introduced. Self-image losses
    /// are included even for an edgeless graph. Gaussian laws remain excluded
    /// by preparation rather than having their tails truncated.
    ///
    /// Each run validates strict positive wall clearance, recomputes N self
    /// losses and E pair maps from the current law, and caches their relaxation
    /// fractions and decay exponentials once. Rebinding viscosity or changing
    /// the wall between runs cannot reuse stale coefficients. The canonical
    /// forward/reverse order includes interleaved self losses, matching the
    /// dense wall path's particle bits for identical admitted inputs.
    ///
    /// An active S-step run charges N + (N+E) + S(2N+4E+1) work units and
    /// S(2N+4E) pair/map visits. Owned payload is N particles, N pairs of f64
    /// loss data, and E quintuplets of f64 pair data (40N+40E bytes on a target
    /// with eight-byte f64). The borrowed graph is excluded. No step allocates
    /// or evaluates a pair kernel. Zero steps, dt or viscosity skip coefficient
    /// allocation and all pair work, charging N+S units and one state vector.
    ///
    /// Bulk circulation can diffuse through the wall. Positivity, enstrophy
    /// contraction, finite-cloud and temporal-error limits are those of the
    /// original wall law; this is not a no-slip or moving-wall boundary model.
    ///
    /// # Errors
    /// Refuses periodic graphs, original prepared-run failures, invalid wall
    /// clearance or reflection, and nonrepresentable wall maps. Inputs and
    /// graph remain unchanged; no partially advanced cloud is published.
    #[allow(clippy::too_many_lines)]
    pub fn diffuse_at_wall_with_cx(
        &self,
        cx: &Cx<'_>,
        particles: &[VortexParticle],
        dt: f64,
        steps: usize,
        wall: PlanarWall2,
        budget: VpmBudget,
    ) -> Result<VpmRun, VpmError> {
        if self.periodic_domain().is_some() {
            return Err(invalid_scalar(
                "planar wall on periodic diffusion graph",
                dt,
                None,
            ));
        }
        let n = self.positions.len();
        if particles.len() != n {
            return Err(VpmError::LengthMismatch {
                field: "prepared particles",
                expected: n,
                actual: particles.len(),
            });
        }
        validate_inputs(n, &self.areas, dt)?;
        admit_cap("particles", n, budget.max_particles)?;
        admit_cap("steps", steps, budget.max_steps)?;
        let active = steps != 0 && dt > 0.0 && self.law.kinematic_viscosity_m2_s() > 0.0;
        let self_count = if active { n } else { 0 };
        let pair_count = if active { self.edges.len() } else { 0 };
        let h = count(self_count, "prepared wall self images")?;
        let e = count(pair_count, "prepared wall edges")?;
        let s = count(steps, "prepared wall steps")?;
        let pair_evaluations = e
            .checked_mul(4)
            .and_then(|v| h.checked_mul(2).and_then(|h| v.checked_add(h)))
            .and_then(|v| v.checked_mul(s))
            .ok_or(VpmError::PlanOverflow {
                resource: "prepared wall pair evaluations",
            })?;
        admit_count(
            "pair evaluations",
            pair_evaluations,
            budget.max_pair_evaluations,
        )?;
        let work_units = count(n, "prepared wall particles")?
            .checked_add(h)
            .and_then(|v| v.checked_add(e))
            .and_then(|v| v.checked_add(pair_evaluations))
            .and_then(|v| v.checked_add(s))
            .ok_or(VpmError::PlanOverflow {
                resource: "prepared wall total work",
            })?;
        let peak_live_bytes = payload(
            n,
            core::mem::size_of::<VortexParticle>(),
            "prepared wall state bytes",
        )?
        .checked_add(payload(
            self_count,
            core::mem::size_of::<CachedLoss>(),
            "prepared wall self bytes",
        )?)
        .and_then(|v| {
            pair_count
                .checked_mul(core::mem::size_of::<PairMap>())
                .and_then(|p| v.checked_add(p))
        })
        .ok_or(VpmError::PlanOverflow {
            resource: "prepared wall live bytes",
        })?;
        admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
        let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
        accountant.checkpoint("fs-vpm.wall-prepared-admission", cx)?;
        let mut state = reserved_vec("wall-prepared-state", n)?;
        let mut self_losses = reserved_vec("wall-prepared-self-losses", self_count)?;
        let mut maps = reserved_vec("wall-prepared-pair-maps", pair_count)?;
        let mut meter = WorkMeter::default();
        for (i, particle) in particles.iter().copied().enumerate() {
            meter.unit(&mut accountant, cx, "fs-vpm.wall-prepared-input", || {
                validate_particle(particle, i)?;
                for axis in 0..2 {
                    if particle.pos[axis].to_bits() != self.positions[i][axis].to_bits() {
                        let field = if axis == 0 {
                            "prepared position.x"
                        } else {
                            "prepared position.y"
                        };
                        return Err(invalid_scalar(field, particle.pos[axis], Some(i)));
                    }
                }
                wall.admit_position(particle.pos, i)?;
                state.push(particle);
                Ok(())
            })?;
        }
        meter.boundary(&mut accountant, cx, "fs-vpm.wall-prepared-input-complete")?;
        if active {
            let kernel = PairKernel::new(self.law, dt, 0.5);
            for (i, &particle) in state.iter().enumerate() {
                meter.unit(&mut accountant, cx, "fs-vpm.wall-prepared-self-map", || {
                    let rate = pair_rates(&kernel, wall, particle, particle, [self.areas[i]; 2]);
                    self_losses.push(CachedLoss::new(rate.loss[0], i)?);
                    Ok(())
                })?;
            }
            for &[i, j] in &self.edges {
                meter.unit(&mut accountant, cx, "fs-vpm.wall-prepared-pair-map", || {
                    let rate = pair_rates(
                        &kernel,
                        wall,
                        state[i],
                        state[j],
                        [self.areas[i], self.areas[j]],
                    );
                    if !rate.exchange.is_finite() || !(0.0..=1.0).contains(&rate.exchange) {
                        return Err(VpmError::NonFiniteResult {
                            phase: "fs-vpm.wall-prepared-pair-map",
                            particle: i,
                        });
                    }
                    maps.push(PairMap {
                        exchange: rate.exchange,
                        losses: [
                            CachedLoss::new(rate.loss[0], i)?,
                            CachedLoss::new(rate.loss[1], j)?,
                        ],
                    });
                    Ok(())
                })?;
            }
            meter.boundary(&mut accountant, cx, "fs-vpm.wall-prepared-maps-complete")?;
        }
        let mut completed_pairs = 0_u64;
        for _ in 0..steps {
            if active {
                sweep(
                    &mut state,
                    &self.areas,
                    &self.edges,
                    &self_losses,
                    &maps,
                    &mut completed_pairs,
                    &mut meter,
                    &mut accountant,
                    cx,
                )?;
            }
            meter.unit(&mut accountant, cx, "fs-vpm.wall-prepared-step", || Ok(()))?;
            meter.boundary(&mut accountant, cx, "fs-vpm.wall-prepared-step-complete")?;
        }
        meter.flush(&mut accountant, "fs-vpm.wall-prepared-finalize")?;
        if meter.completed != work_units
            || completed_pairs != pair_evaluations
            || accountant.consumption().cost_charged != work_units
        {
            return Err(VpmError::AccountingMismatch {
                planned: work_units,
                completed: meter.completed,
                charged: accountant.consumption().cost_charged,
            });
        }
        accountant.checkpoint("fs-vpm.wall-prepared-publication", cx)?;
        Ok(VpmRun {
            particles: state,
            steps_completed: steps,
            pair_evaluations,
            work_units,
            peak_live_bytes,
            budget: accountant.consumption(),
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn sweep(
    state: &mut [VortexParticle],
    areas: &[f64],
    edges: &[[usize; 2]],
    self_losses: &[CachedLoss],
    maps: &[PairMap],
    visits: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    // Keep self maps at exactly their dense palindrome positions. Grouping all
    // diagonal losses at either end would change the time-splitting operator.
    let mut at = 0;
    for (i, &loss) in self_losses.iter().enumerate() {
        visit(visits, meter, accountant, cx, || {
            state[i].circulation = loss.apply(state[i].circulation, i)?;
            Ok(())
        })?;
        while at < edges.len() && edges[at][0] == i {
            let [i, j] = edges[at];
            let map = &maps[at];
            visit(visits, meter, accountant, cx, || {
                exchange_fraction(state, areas, map.exchange, i, j)
            })?;
            visit(visits, meter, accountant, cx, || map.loss(state, [i, j]))?;
            at += 1;
        }
    }
    for (i, &loss) in self_losses.iter().enumerate().rev() {
        while at > 0 && edges[at - 1][0] == i {
            at -= 1;
            let [i, j] = edges[at];
            let map = &maps[at];
            visit(visits, meter, accountant, cx, || map.loss(state, [i, j]))?;
            visit(visits, meter, accountant, cx, || {
                exchange_fraction(state, areas, map.exchange, i, j)
            })?;
        }
        visit(visits, meter, accountant, cx, || {
            state[i].circulation = loss.apply(state[i].circulation, i)?;
            Ok(())
        })?;
    }
    meter.boundary(accountant, cx, "fs-vpm.wall-prepared-sweep-complete")
}

fn visit(
    visits: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
    apply: impl FnOnce() -> Result<(), VpmError>,
) -> Result<(), VpmError> {
    meter.unit(accountant, cx, "fs-vpm.wall-prepared-pair", || {
        *visits = visits.checked_add(1).ok_or(VpmError::PlanOverflow {
            resource: "executed prepared wall pairs",
        })?;
        apply()
    })
}
