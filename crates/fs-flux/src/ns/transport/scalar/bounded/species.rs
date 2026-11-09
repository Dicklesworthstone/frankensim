//! Joint passive species transport on the original scalar correction graph.
//!
//! A common diffusivity and velocity make every species use the SAME linear
//! spatial operator. When all nodal and reservoir compositions sum to one,
//! the low-order predictor preserves that sum and the raw correction vectors
//! sum to zero. One common factor per interior/reservoir vector preserves this
//! identity while satisfying every species' available positive/negative budget.
//! Independently limiting, clipping or renormalizing species does not have this
//! property. All species here evolve independently; none is reconstructed as a
//! remainder. Identities are admitted to the explicit floating-point allowance.

use super::{
    BoundedScalarControls, BoundedScalarError, BoundedScalarStep, Dd,
    PreparedBoundedScalarTransport, PreparedScalarTransport, ScalarBounds, ScalarTransportError,
    Work as ScalarWork, sum,
};

/// Absolute roundoff allowance for a unit-sum composition at every admitted
/// node, original boundary trace point and completed Euler/RK stage.
pub const SPECIES_SUM_ROUNDOFF: f64 = 1024.0 * f64::EPSILON;

#[derive(Clone, Debug, PartialEq)]
pub enum BoundedSpeciesError {
    Scalar {
        species: usize,
        error: BoundedScalarError,
    },
    Invalid(&'static str),
    /// A different spatial operator cannot share a conservative zero-sum
    /// correction vector. Species-specific molecular diffusion is not admitted.
    OperatorMismatch {
        species: usize,
    },
    StateSum {
        dof: usize,
        sum: f64,
    },
    BoundarySum {
        edge: usize,
        sum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for BoundedSpeciesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bounded DG species: {self:?}")
    }
}
impl std::error::Error for BoundedSpeciesError {}
type Result<T> = core::result::Result<T, BoundedSpeciesError>;

fn scalar_error(species: usize, error: BoundedScalarError) -> BoundedSpeciesError {
    if error == BoundedScalarError::Native(ScalarTransportError::Cancelled) {
        BoundedSpeciesError::Cancelled
    } else {
        BoundedSpeciesError::Scalar { species, error }
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(BoundedSpeciesError::Cancelled)
    } else {
        Ok(())
    }
}

/// Species retain the input order. Each independent inventory includes its
/// actual original advective/SIP and unremoved numerical boundary exchanges.
pub struct BoundedSpeciesStep {
    pub species: Vec<BoundedScalarStep>,
    /// Largest absolute unit-sum defect over input, both Euler stages and all
    /// accepted RK endpoints; it is reported without modifying any fraction.
    pub maximum_sum_residual: f64,
}

/// One shared correction graph for frozen, constant-density passive mixing.
/// Original native operators retain species-specific boundary compositions.
/// Molecular plus SGS diffusivity must be identical for every species: this is
/// an equal-diffusivity model, not Maxwell--Stefan, reacting or variable-density
/// transport. Multiply fraction-area inventories by common rho and depth for
/// species masses. No thermodynamic or species identity is inferred here.
pub struct PreparedBoundedSpeciesTransport<'a, 'm> {
    graph: PreparedBoundedScalarTransport<'a, 'm>,
    native: Vec<&'a PreparedScalarTransport<'m>>,
}

impl<'a, 'm> PreparedBoundedSpeciesTransport<'a, 'm> {
    /// Bind at least two original scalar operators. Their mesh, diffusivity,
    /// spatial matrix and complete weak boundary coefficients must agree.
    /// Value/inflow traces may differ by species but must describe a simplex at
    /// the original Gauss samples AND physical affine endpoint extrema.
    /// Nonzero prescribed fluxes and volumetric sources are outside this lane.
    /// Setup stores one O(nnz) graph plus O(species) borrowed native operators.
    pub fn new(
        species: &[&'a PreparedScalarTransport<'m>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        checkpoint(&mut cancelled)?;
        if species.len() < 2 {
            return Err(BoundedSpeciesError::Invalid(
                "at least two transported species required",
            ));
        }
        let first = species[0];
        let bounds = ScalarBounds {
            lower: 0.0,
            upper: 1.0,
        };
        let graph = first
            .prepare_bounded_transport(bounds, &mut cancelled)
            .map_err(|e| scalar_error(0, e))?;
        for (s, native) in species.iter().enumerate().skip(1) {
            checkpoint(&mut cancelled)?;
            if !core::ptr::eq(first.mesh, native.mesh)
                || first.diffusivity_m2_s != native.diffusivity_m2_s
                || first.operator.nrows() != native.operator.nrows()
                || first.boundary_fluxes.len() != native.boundary_fluxes.len()
            {
                return Err(BoundedSpeciesError::OperatorMismatch { species: s });
            }
            for row in 0..first.operator.nrows() {
                checkpoint(&mut cancelled)?;
                if first.operator.row(row) != native.operator.row(row) {
                    return Err(BoundedSpeciesError::OperatorMismatch { species: s });
                }
            }
            for (a, b) in first.boundary_fluxes.iter().zip(&native.boundary_fluxes) {
                checkpoint(&mut cancelled)?;
                if b.prescribed_flux {
                    return Err(scalar_error(
                        s,
                        BoundedScalarError::PrescribedFlux { edge: b.edge },
                    ));
                }
                if a.edge != b.edge
                    || a.cell != b.cell
                    || a.advection != b.advection
                    || a.diffusion != b.diffusion
                    || a.reservoir_coefficients != b.reservoir_coefficients
                {
                    return Err(BoundedSpeciesError::OperatorMismatch { species: s });
                }
            }
        }
        for (face_index, face) in first.boundary_fluxes.iter().enumerate() {
            if !face
                .reservoir_coefficients
                .iter()
                .flatten()
                .any(|&b| b != 0.0)
            {
                continue;
            }
            for fraction in [
                0.0,
                1.0,
                0.5 - 0.5 * 3.0_f64.sqrt(),
                0.5 + 0.5 * 3.0_f64.sqrt(),
            ] {
                checkpoint(&mut cancelled)?;
                let mut total = Dd::from_f64(0.0);
                for (s, native) in species.iter().enumerate() {
                    let values = native.boundary_fluxes[face_index].reservoir_values;
                    let value = (1.0 - fraction) * values[0] + fraction * values[1];
                    if !bounds.contains(value) {
                        return Err(scalar_error(
                            s,
                            BoundedScalarError::BoundaryBounds {
                                edge: face.edge,
                                value,
                            },
                        ));
                    }
                    total = total + Dd::from_f64(value);
                }
                let total = total.to_f64();
                if !total.is_finite() || (total - 1.0).abs() > SPECIES_SUM_ROUNDOFF {
                    return Err(BoundedSpeciesError::BoundarySum {
                        edge: face.edge,
                        sum: total,
                    });
                }
            }
        }
        checkpoint(&mut cancelled)?;
        Ok(Self {
            graph,
            native: species.to_vec(),
        })
    }

    /// Joint SSPRK2 with the original bounded scalar CFL, dyadic subcycling and
    /// accepted boundary accounting. `previous[s][i]` is species s at P1 DOF i.
    /// All fractions stay in [0,1] and sum to one to declared roundoff; each
    /// species conserves its own inventory including accepted boundary fluxes.
    /// Work is O(species*substeps*nnz), with O(species*(nnz+DOFs)) scratch reused
    /// across stages. No stage allocation, scalar normalization or partial
    /// publication occurs. Changing velocity/diffusivity/boundaries needs a new
    /// prepared owner. Shared limiting may reduce order at composition fronts.
    pub fn ssprk2_step(
        &self,
        previous: &[&[f64]],
        timestep_s: f64,
        controls: BoundedScalarControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoundedSpeciesStep> {
        let mut maximum_sum_residual = self.admit_composition(previous, &mut cancelled)?;
        let (count, dt) = self
            .graph
            .substeps(timestep_s, controls)
            .map_err(|e| scalar_error(0, e))?;
        let n = self.graph.mass.len();
        let zero = Dd::from_f64(0.0);
        let mut current: Vec<_> = previous.iter().map(|v| v.to_vec()).collect();
        let mut stage = vec![vec![0.0; n]; self.native.len()];
        let mut second = stage.clone();
        let mut work = SpeciesWork {
            scalar: self.native.iter().map(|_| self.graph.work()).collect(),
            pair_factors: vec![1.0; self.graph.pairs.len()],
            reservoir_factors: vec![1.0; self.graph.reservoirs.len()],
            boundaries: vec![
                vec![[zero; 3]; self.native[0].boundary_fluxes.len()];
                self.native.len()
            ],
            limited: vec![0; self.native.len()],
        };
        for _ in 0..count {
            maximum_sum_residual = maximum_sum_residual.max(self.euler(
                &current,
                dt,
                &mut stage,
                &mut work,
                &mut cancelled,
            )?);
            maximum_sum_residual = maximum_sum_residual.max(self.euler(
                &stage,
                dt,
                &mut second,
                &mut work,
                &mut cancelled,
            )?);
            for (old, new) in current.iter_mut().zip(&second) {
                checkpoint(&mut cancelled)?;
                for (old, &new) in old.iter_mut().zip(new) {
                    *old = 0.5 * *old + 0.5 * new;
                }
            }
            maximum_sum_residual =
                maximum_sum_residual.max(self.admit_composition(&current, &mut cancelled)?);
        }
        let mut species = Vec::with_capacity(self.native.len());
        for (s, (current, boundaries)) in current.into_iter().zip(work.boundaries).enumerate() {
            species.push(
                self.graph
                    .finish_step(
                        self.native[s],
                        previous[s],
                        current,
                        timestep_s,
                        count,
                        dt,
                        work.limited[s],
                        boundaries,
                        &mut cancelled,
                    )
                    .map_err(|e| scalar_error(s, e))?,
            );
        }
        checkpoint(&mut cancelled)?;
        Ok(BoundedSpeciesStep {
            species,
            maximum_sum_residual,
        })
    }

    fn admit_composition(
        &self,
        values: &[impl AsRef<[f64]>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        checkpoint(cancelled)?;
        if values.len() != self.native.len() {
            return Err(BoundedSpeciesError::Invalid(
                "one P1 field for each prepared species required",
            ));
        }
        for (s, field) in values.iter().enumerate() {
            self.graph
                .admit_values(field.as_ref())
                .map_err(|e| scalar_error(s, e))?;
        }
        let mut maximum = 0.0_f64;
        for dof in 0..self.graph.mass.len() {
            checkpoint(cancelled)?;
            let total = sum(values.iter().map(|v| v.as_ref()[dof]));
            let defect = (total - 1.0).abs();
            if !total.is_finite() || defect > SPECIES_SUM_ROUNDOFF {
                return Err(BoundedSpeciesError::StateSum { dof, sum: total });
            }
            maximum = maximum.max(defect);
        }
        Ok(maximum)
    }

    fn euler(
        &self,
        old: &[Vec<f64>],
        dt: f64,
        output: &mut [Vec<f64>],
        work: &mut SpeciesWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        work.pair_factors.fill(1.0);
        work.reservoir_factors.fill(1.0);
        for (s, scalar) in work.scalar.iter_mut().enumerate() {
            self.graph
                .prepare_euler(self.native[s], &old[s], dt, scalar, cancelled)
                .map_err(|e| scalar_error(s, e))?;
            for (index, factor) in work.pair_factors.iter_mut().enumerate() {
                checkpoint(cancelled)?;
                if scalar.pair_flux[index] != 0.0 {
                    *factor = factor.min(self.graph.pair_factor(index, scalar));
                }
            }
            for (index, factor) in work.reservoir_factors.iter_mut().enumerate() {
                checkpoint(cancelled)?;
                if scalar.reservoir_flux[index] != 0.0 {
                    *factor = factor.min(self.graph.reservoir_factor(index, scalar));
                }
            }
        }
        for (s, scalar) in work.scalar.iter_mut().enumerate() {
            work.limited[s] += self
                .graph
                .apply_euler(
                    self.native[s],
                    &old[s],
                    &mut output[s],
                    scalar,
                    Some((&work.pair_factors, &work.reservoir_factors)),
                    cancelled,
                )
                .map_err(|e| scalar_error(s, e))?;
            self.graph
                .accumulate_boundaries(
                    self.native[s],
                    &old[s],
                    0.5 * dt,
                    &scalar.boundary_extra,
                    &mut work.boundaries[s],
                )
                .map_err(|e| scalar_error(s, e))?;
        }
        self.admit_composition(output, cancelled)
    }
}

struct SpeciesWork {
    scalar: Vec<ScalarWork>,
    pair_factors: Vec<f64>,
    reservoir_factors: Vec<f64>,
    boundaries: Vec<Vec<[Dd; 3]>>,
    limited: Vec<usize>,
}
