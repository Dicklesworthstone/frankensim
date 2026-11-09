//! Conservative bound-preserving flux correction of the ORIGINAL DG operator.
//!
//! M dc/dt + A c = B g is the prepared scalar equation. Symmetric graph
//! diffusion d_ij=max(A_ij,A_ji,0) and lumped P1 mass make its Euler predictor
//! monotone under the admitted CFL. Signed weak boundary coefficients receive
//! the same treatment. Limited antisymmetric corrections restore consistent
//! mass and high-order transport where the available scalar inventory permits.
//! SSPRK2 composes two bounded Euler stages; subcycling respects diffusion as
//! well as advection. No endpoint clipping or inventory repair is performed.
//!
//! Algebraic FCT and the consistent-mass correction follow Kuzmin (2009),
//! doi:10.1016/j.jcp.2008.11.019, equations (30)--(34) at theta=0; the paired
//! inventory limiter is Zalesak (1979), doi:10.1016/0021-9991(79)90051-2.
//! Here the graph is the native discontinuous P1/SIP operator, not a replacement
//! geometric two-point stencil. Bounds hold up to explicitly admitted floating
//! point roundoff. Limiting supplies numerical diffusion, not physical SGS heat.

use super::{PreparedScalarTransport, ScalarTransportError, finite, poll, sum};
use fs_math::dd::Dd;
use std::collections::BTreeSet;

#[cfg(test)]
mod tests;

/// Finite authored bounds in the scalar's original units. P1 nodal bounds also
/// bound the complete affine field inside each original triangle.
#[derive(Clone, Copy, Debug)]
pub struct ScalarBounds {
    pub lower: f64,
    pub upper: f64,
}
impl ScalarBounds {
    fn roundoff(self) -> f64 {
        128.0
            * f64::EPSILON
            * self
                .lower
                .abs()
                .max(self.upper.abs())
                .max(f64::MIN_POSITIVE)
    }
    fn contains(self, value: f64) -> bool {
        value.is_finite()
            && value >= self.lower - self.roundoff()
            && value <= self.upper + self.roundoff()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BoundedScalarControls {
    /// Allowed outgoing fraction in each mass-lumped Euler predictor, (0,0.5].
    pub maximum_courant: f64,
    /// Maximum dyadic substeps for the COMPLETE requested interval.
    pub maximum_substeps: usize,
}
impl Default for BoundedScalarControls {
    fn default() -> Self {
        Self {
            maximum_courant: 0.45,
            maximum_substeps: 4096,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BoundedScalarError {
    Native(ScalarTransportError),
    Invalid(&'static str),
    /// A prescribed heat flux needs a separately admitted source-bound policy.
    PrescribedFlux {
        edge: usize,
    },
    BoundaryBounds {
        edge: usize,
        value: f64,
    },
    StateBounds {
        dof: usize,
        value: f64,
    },
    /// The original discrete incompressible constant-mode identity failed.
    ConstantMode {
        dof: usize,
        defect: f64,
        allowance: f64,
    },
    SubstepLimit {
        required: usize,
        maximum: usize,
    },
}
impl core::fmt::Display for BoundedScalarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bounded DG scalar: {self:?}")
    }
}
impl std::error::Error for BoundedScalarError {}
impl From<ScalarTransportError> for BoundedScalarError {
    fn from(error: ScalarTransportError) -> Self {
        Self::Native(error)
    }
}
type Result<T> = core::result::Result<T, BoundedScalarError>;

/// Time-integrated outward scalar amount [scalar*m²] per out-of-plane depth.
#[derive(Clone, Copy, Debug)]
pub struct BoundedScalarBoundaryTransfer {
    pub edge: usize,
    pub advection: f64,
    pub diffusion: f64,
    /// Extra boundary exchange from graph diffusion NOT removed by limiting.
    /// This is numerical transport, separate from authored molecular/SGS mixing.
    pub limiter_diffusion: f64,
}

pub struct BoundedScalarStep {
    pub values: Vec<f64>,
    pub timestep_s: f64,
    pub substeps: usize,
    /// Pair/boundary corrections reduced by a limiter across both RK stages.
    pub limited_fluxes: usize,
    pub maximum_courant: f64,
    pub minimum: f64,
    pub maximum: f64,
    pub integral_before: f64,
    pub integral_after: f64,
    pub boundaries: Vec<BoundedScalarBoundaryTransfer>,
    /// Direct integrated storage change plus all outward boundary transfers.
    /// Roundoff is retained and never fed back into the scalar state.
    pub balance_residual: f64,
}

struct Pair {
    i: usize,
    j: usize,
    diffusion: f64,
    mass: f64,
    into_i: f64,
    into_j: f64,
}
struct Reservoir {
    dof: usize,
    face: usize,
    value: f64,
    diffusion: f64,
    incoming: f64,
}

/// Sparse graph prepared once for frozen velocity, diffusivity and boundary
/// data. It borrows the original scalar operator without changing its BDF path.
pub struct PreparedBoundedScalarTransport<'a, 'm> {
    native: &'a PreparedScalarTransport<'m>,
    bounds: ScalarBounds,
    mass: Vec<f64>,
    diagonal: Vec<f64>,
    pairs: Vec<Pair>,
    reservoirs: Vec<Reservoir>,
    maximum_rate_per_s: f64,
}

impl<'m> PreparedScalarTransport<'m> {
    /// Admit source-free bounded mixing with value boundaries or homogeneous
    /// diffusive-flux boundaries. Incoming reservoirs retain the original
    /// affine trace, including its endpoint extrema and both reversing streams.
    /// Nonzero imposed flux is refused, not silently treated as a reservoir.
    ///
    /// Every original weak boundary coefficient is retained, including signed
    /// SIP and affine-trace contributions. Positive graph corrections preserve
    /// the constant mode; admission checks the original row identity before
    /// claiming an invariant interval. No diagonal or velocity repair occurs.
    /// Setup uses O(nnz) storage and O(nnz log(nnz)) canonical pair ordering.
    pub fn prepare_bounded_transport(
        &self,
        bounds: ScalarBounds,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedBoundedScalarTransport<'_, 'm>> {
        poll(&mut cancelled)?;
        if !bounds.lower.is_finite()
            || !bounds.upper.is_finite()
            || bounds.upper < bounds.lower
            || !(bounds.upper - bounds.lower).is_finite()
        {
            return Err(BoundedScalarError::Invalid(
                "finite ordered scalar bounds required",
            ));
        }
        let n = self.operator.nrows();
        let mut mass = Vec::with_capacity(n);
        let mut diagonal = vec![0.0; n];
        let mut row_sum = vec![0.0; n];
        let mut scale = vec![0.0; n];
        let mut edges = BTreeSet::new();
        for i in 0..n {
            poll(&mut cancelled)?;
            mass.push(self.mesh.areas[i / 3] / 3.0);
            let (columns, values) = self.operator.row(i);
            row_sum[i] = sum(values.iter().copied());
            scale[i] = sum(values.iter().map(|x| x.abs()));
            for (&j, &coefficient) in columns.iter().zip(values) {
                if i == j {
                    diagonal[i] = coefficient;
                } else {
                    edges.insert((i.min(j), i.max(j)));
                }
            }
            // The consistent P1 mass couples every pair in the SAME triangle,
            // even if its spatial coefficient happens to vanish exactly.
            for j in 3 * (i / 3)..3 * (i / 3) + 3 {
                if i != j {
                    edges.insert((i.min(j), i.max(j)));
                }
            }
        }
        let coefficient = |i, j| {
            let (columns, values) = self.operator.row(i);
            columns.binary_search(&j).map_or(0.0, |index| values[index])
        };
        let mut pairs = Vec::with_capacity(edges.len());
        for (i, j) in edges {
            poll(&mut cancelled)?;
            let a = coefficient(i, j);
            let b = coefficient(j, i);
            let diffusion = a.max(b).max(0.0);
            diagonal[i] += diffusion;
            diagonal[j] += diffusion;
            pairs.push(Pair {
                i,
                j,
                diffusion,
                mass: if i / 3 == j / 3 { mass[i] / 4.0 } else { 0.0 },
                into_i: diffusion - a,
                into_j: diffusion - b,
            });
        }
        let mut reservoirs = Vec::new();
        for (face_index, face) in self.boundary_fluxes.iter().enumerate() {
            poll(&mut cancelled)?;
            if face.prescribed_flux {
                return Err(BoundedScalarError::PrescribedFlux { edge: face.edge });
            }
            if face
                .reservoir_coefficients
                .iter()
                .flatten()
                .any(|&b| b != 0.0)
            {
                // Original Gauss fractions at the physical endpoints. Checking
                // only the interior samples could admit an invalid inflow trace.
                for fraction in [0.5 - 0.5 * 3.0_f64.sqrt(), 0.5 + 0.5 * 3.0_f64.sqrt()] {
                    let value = (1.0 - fraction) * face.reservoir_values[0]
                        + fraction * face.reservoir_values[1];
                    if !bounds.contains(value) {
                        return Err(BoundedScalarError::BoundaryBounds {
                            edge: face.edge,
                            value,
                        });
                    }
                }
            }
            for local in 0..3 {
                for q in 0..2 {
                    let b = finite(
                        face.reservoir_coefficients[local][q],
                        "scalar reservoir coefficient",
                    )?;
                    if b == 0.0 {
                        continue;
                    }
                    let i = 3 * face.cell + local;
                    row_sum[i] -= b;
                    scale[i] += b.abs();
                    let diffusion = (-b).max(0.0);
                    diagonal[i] += diffusion;
                    reservoirs.push(Reservoir {
                        dof: i,
                        face: face_index,
                        value: face.reservoir_values[q],
                        diffusion,
                        incoming: b.max(0.0),
                    });
                }
            }
        }
        let mut maximum_rate_per_s = 0.0_f64;
        for i in 0..n {
            poll(&mut cancelled)?;
            let allowance = 512.0 * f64::EPSILON * scale[i];
            if !row_sum[i].is_finite() || !allowance.is_finite() || row_sum[i].abs() > allowance {
                return Err(BoundedScalarError::ConstantMode {
                    dof: i,
                    defect: row_sum[i],
                    allowance,
                });
            }
            if !diagonal[i].is_finite()
                || diagonal[i] < 0.0
                || !mass[i].is_finite()
                || mass[i] <= 0.0
            {
                return Err(BoundedScalarError::Invalid(
                    "finite nonnegative monotone diagonal and positive P1 storage required",
                ));
            }
            maximum_rate_per_s =
                maximum_rate_per_s.max(finite(diagonal[i] / mass[i], "bounded scalar rate")?);
        }
        poll(&mut cancelled)?;
        Ok(PreparedBoundedScalarTransport {
            native: self,
            bounds,
            mass,
            diagonal,
            pairs,
            reservoirs,
            maximum_rate_per_s,
        })
    }
}

struct Work {
    high_increment: Vec<f64>,
    low_sum: Vec<Dd>,
    low: Vec<f64>,
    positive: Vec<f64>,
    negative: Vec<f64>,
    r_plus: Vec<f64>,
    r_minus: Vec<f64>,
    correction: Vec<Dd>,
    pair_flux: Vec<f64>,
    reservoir_flux: Vec<f64>,
    boundary_extra: Vec<Dd>,
}
impl PreparedBoundedScalarTransport<'_, '_> {
    /// Advance the complete physical interval with bounded dyadic SSPRK2
    /// substeps. Both Euler stages and the returned affine fields retain the
    /// authored bounds to rounding accuracy. The consistent P1 mass and native
    /// high-order flux are recovered wherever limiting is inactive.
    ///
    /// All pair corrections conserve scalar inventory. Boundary corrections
    /// appear exactly once in the returned exchange account. This is passive
    /// source-free transport; it is NOT an implicit BDF step or a thermal source
    /// policy. Time-dependent boundary/velocity data require a newly prepared
    /// endpoint. Input values are read-only, and failure publishes no substeps.
    /// Work is O(substeps*nnz), with reusable O(nnz+DOFs) scratch per invocation;
    /// there are no stage allocations or dense/global nonlinear solves.
    pub fn ssprk2_step(
        &self,
        previous: &[f64],
        timestep_s: f64,
        controls: BoundedScalarControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoundedScalarStep> {
        poll(&mut cancelled)?;
        self.admit_values(previous)?;
        if !timestep_s.is_finite()
            || timestep_s <= 0.0
            || !controls.maximum_courant.is_finite()
            || controls.maximum_courant <= 0.0
            || controls.maximum_courant > 0.5
            || controls.maximum_substeps == 0
        {
            return Err(BoundedScalarError::Invalid(
                "positive interval, Courant in (0,0.5] and nonzero work budget required",
            ));
        }
        let mut count = 1_usize;
        while timestep_s / count as f64 * self.maximum_rate_per_s > controls.maximum_courant {
            count = count.checked_mul(2).ok_or(BoundedScalarError::Invalid(
                "bounded substep count overflow",
            ))?;
            if count > controls.maximum_substeps {
                return Err(BoundedScalarError::SubstepLimit {
                    required: count,
                    maximum: controls.maximum_substeps,
                });
            }
        }
        let dt = timestep_s / count as f64;
        if dt <= 0.0 || !dt.is_finite() {
            return Err(BoundedScalarError::Invalid(
                "representable positive substep required",
            ));
        }
        let n = self.mass.len();
        let zero = Dd::from_f64(0.0);
        let mut work = Work {
            high_increment: vec![0.0; n],
            low_sum: vec![zero; n],
            low: vec![0.0; n],
            positive: vec![0.0; n],
            negative: vec![0.0; n],
            r_plus: vec![0.0; n],
            r_minus: vec![0.0; n],
            correction: vec![zero; n],
            pair_flux: vec![0.0; self.pairs.len()],
            reservoir_flux: vec![0.0; self.reservoirs.len()],
            boundary_extra: vec![zero; self.native.boundary_fluxes.len()],
        };
        let mut current = previous.to_vec();
        let mut stage = vec![0.0; n];
        let mut second = vec![0.0; n];
        let mut boundary_integrals = vec![[zero; 3]; self.native.boundary_fluxes.len()];
        let mut limited_fluxes = 0;
        for _ in 0..count {
            limited_fluxes += self.euler(&current, dt, &mut stage, &mut work, &mut cancelled)?;
            self.accumulate_boundaries(
                &current,
                0.5 * dt,
                &work.boundary_extra,
                &mut boundary_integrals,
            )?;
            limited_fluxes += self.euler(&stage, dt, &mut second, &mut work, &mut cancelled)?;
            self.accumulate_boundaries(
                &stage,
                0.5 * dt,
                &work.boundary_extra,
                &mut boundary_integrals,
            )?;
            for i in 0..n {
                current[i] = 0.5 * current[i] + 0.5 * second[i];
            }
            self.admit_values(&current)?;
            poll(&mut cancelled)?;
        }
        let boundaries: Vec<_> = self
            .native
            .boundary_fluxes
            .iter()
            .zip(boundary_integrals)
            .map(|(face, integrals)| BoundedScalarBoundaryTransfer {
                edge: face.edge,
                advection: integrals[0].to_f64(),
                diffusion: integrals[1].to_f64(),
                limiter_diffusion: integrals[2].to_f64(),
            })
            .collect();
        let integral_before = finite(
            sum(previous.iter().zip(&self.mass).map(|(c, m)| c * m)),
            "bounded initial inventory",
        )?;
        let integral_after = finite(
            sum(current.iter().zip(&self.mass).map(|(c, m)| c * m)),
            "bounded final inventory",
        )?;
        let storage = sum(current
            .iter()
            .zip(previous)
            .zip(&self.mass)
            .map(|((new, old), m)| (new - old) * m));
        let balance_residual = finite(
            sum(std::iter::once(storage).chain(
                boundaries
                    .iter()
                    .flat_map(|face| [face.advection, face.diffusion, face.limiter_diffusion]),
            )),
            "bounded scalar balance",
        )?;
        let minimum = current.iter().copied().fold(f64::INFINITY, f64::min);
        let maximum = current.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        poll(&mut cancelled)?;
        Ok(BoundedScalarStep {
            values: current,
            timestep_s,
            substeps: count,
            limited_fluxes,
            maximum_courant: dt * self.maximum_rate_per_s,
            minimum,
            maximum,
            integral_before,
            integral_after,
            boundaries,
            balance_residual,
        })
    }

    fn admit_values(&self, values: &[f64]) -> Result<()> {
        if values.len() != self.mass.len() {
            return Err(BoundedScalarError::Invalid(
                "three P1 values per original triangle required",
            ));
        }
        if let Some(dof) = values.iter().position(|&v| !self.bounds.contains(v)) {
            return Err(BoundedScalarError::StateBounds {
                dof,
                value: values[dof],
            });
        }
        Ok(())
    }

    fn euler(
        &self,
        old: &[f64],
        dt: f64,
        output: &mut [f64],
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize> {
        poll(cancelled)?;
        self.native.operator.spmv(old, &mut work.high_increment);
        // Exact inverse of the ORIGINAL 3x3 P1 triangle mass A/12*(I+11').
        // Keep this consistent-mass increment in the antidiffusive fluxes.
        for cell in 0..self.native.mesh.tris.len() {
            poll(cancelled)?;
            let r: [f64; 3] = core::array::from_fn(|i| {
                self.native.boundary_load[3 * cell + i] - work.high_increment[3 * cell + i]
            });
            let average = 0.25 * sum(r.into_iter());
            for (i, &residual) in r.iter().enumerate() {
                work.high_increment[3 * cell + i] = finite(
                    (dt / self.mass[3 * cell]) * 4.0 * (residual - average),
                    "consistent scalar increment",
                )?;
            }
        }
        let zero = Dd::from_f64(0.0);
        work.positive.fill(0.0);
        work.negative.fill(0.0);
        work.correction.fill(zero);
        work.boundary_extra.fill(zero);
        for i in 0..old.len() {
            work.low_sum[i] = Dd::from_f64((1.0 - dt * self.diagonal[i] / self.mass[i]) * old[i]);
        }
        for (index, pair) in self.pairs.iter().enumerate() {
            poll(cancelled)?;
            let (i, j) = (pair.i, pair.j);
            work.low_sum[i] =
                work.low_sum[i] + Dd::from_f64((dt * pair.into_i / self.mass[i]) * old[j]);
            work.low_sum[j] =
                work.low_sum[j] + Dd::from_f64((dt * pair.into_j / self.mass[j]) * old[i]);
            let flux = finite(
                pair.mass * (work.high_increment[i] - work.high_increment[j])
                    + (dt * pair.diffusion) * (old[i] - old[j]),
                "scalar pair correction",
            )?;
            work.pair_flux[index] = flux;
            work.positive[i] += flux.max(0.0);
            work.negative[i] += (-flux).max(0.0);
            work.positive[j] += (-flux).max(0.0);
            work.negative[j] += flux.max(0.0);
        }
        for (index, reservoir) in self.reservoirs.iter().enumerate() {
            poll(cancelled)?;
            let i = reservoir.dof;
            work.low_sum[i] = work.low_sum[i]
                + Dd::from_f64((dt * reservoir.incoming / self.mass[i]) * reservoir.value);
            let flux = finite(
                (dt * reservoir.diffusion) * (old[i] - reservoir.value),
                "scalar boundary correction",
            )?;
            work.reservoir_flux[index] = flux;
            work.positive[i] += flux.max(0.0);
            work.negative[i] += (-flux).max(0.0);
        }
        for i in 0..old.len() {
            work.low[i] = finite(work.low_sum[i].to_f64(), "monotone scalar predictor")?;
        }
        self.admit_values(&work.low)?;
        for i in 0..old.len() {
            // Clip only admissible FLUX FRACTIONS, never physical scalar states.
            let fraction = |budget: f64, proposed: f64| -> Result<f64> {
                finite(proposed, "scalar correction budget")?;
                Ok(if proposed == 0.0 {
                    1.0
                } else {
                    (budget.max(0.0) / proposed).min(1.0)
                })
            };
            work.r_plus[i] = fraction(
                self.mass[i] * (self.bounds.upper - work.low[i]),
                work.positive[i],
            )?;
            work.r_minus[i] = fraction(
                self.mass[i] * (work.low[i] - self.bounds.lower),
                work.negative[i],
            )?;
        }
        let mut limited = 0;
        for (pair, &flux) in self.pairs.iter().zip(&work.pair_flux) {
            poll(cancelled)?;
            let factor = if flux >= 0.0 {
                work.r_plus[pair.i].min(work.r_minus[pair.j])
            } else {
                work.r_minus[pair.i].min(work.r_plus[pair.j])
            };
            limited += usize::from(flux != 0.0 && factor < 1.0);
            let corrected = Dd::from_f64(factor * flux);
            work.correction[pair.i] = work.correction[pair.i] + corrected;
            work.correction[pair.j] = work.correction[pair.j] - corrected;
        }
        for (reservoir, &flux) in self.reservoirs.iter().zip(&work.reservoir_flux) {
            poll(cancelled)?;
            let factor = if flux >= 0.0 {
                work.r_plus[reservoir.dof]
            } else {
                work.r_minus[reservoir.dof]
            };
            limited += usize::from(flux != 0.0 && factor < 1.0);
            work.correction[reservoir.dof] =
                work.correction[reservoir.dof] + Dd::from_f64(factor * flux);
            // Actual additional outward RATE. SSPRK2 integrates this alongside
            // the original native physical flux at the SAME Euler-stage state.
            work.boundary_extra[reservoir.face] = work.boundary_extra[reservoir.face]
                + Dd::from_f64(
                    (1.0 - factor) * reservoir.diffusion * (old[reservoir.dof] - reservoir.value),
                );
        }
        for i in 0..old.len() {
            output[i] = work.low[i] + work.correction[i].to_f64() / self.mass[i];
        }
        self.admit_values(output)?;
        poll(cancelled)?;
        Ok(limited)
    }

    fn accumulate_boundaries(
        &self,
        state: &[f64],
        weight_s: f64,
        extra: &[Dd],
        total: &mut [[Dd; 3]],
    ) -> Result<()> {
        for (index, face) in self.native.boundary_fluxes.iter().enumerate() {
            let local = &state[3 * face.cell..3 * face.cell + 3];
            let advection =
                sum(face.advection.iter().zip(local).map(|(a, c)| a * c)) + face.advective_load;
            let diffusion =
                sum(face.diffusion.iter().zip(local).map(|(a, c)| a * c)) + face.diffusive_load;
            for (slot, rate) in [advection, diffusion, extra[index].to_f64()]
                .into_iter()
                .enumerate()
            {
                total[index][slot] = total[index][slot]
                    + Dd::from_f64(finite(weight_s * rate, "bounded boundary integral")?);
            }
        }
        Ok(())
    }
}
