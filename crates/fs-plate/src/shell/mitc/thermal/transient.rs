//! Retained interior director inertia for thermally driven shell dynamics.
//!
//! The original static element writes its twenty coordinates as `y = T q`,
//! with eighteen nodal coordinates and an elastic bubble map. Here the two
//! extra coordinates remain independent: `y = T q + B r`. This is an
//! invertible change of coordinates on the ORIGINAL enriched element, not a
//! moving prescribed thermal lift. The INITIAL reference blocks are
//! `K = diag(K_condensed, K_bb)` and `M = [[M_condensed,L],[L^T,M_bb]]`.
//! A material refresh retains the generally nonzero stiffness cross blocks in
//! this original basis; it never diagonalizes them by changing the trajectory.
//! Thermal load is `[T^T f_th, f_th,b]`. Thus a heated internal rotation has
//! its own momentum equation and cannot inject unaccounted actuator work.
//!
//! The native nodal support frames, reference mass, director interpolation,
//! quadrature and numerical drilling stiffness remain unchanged. Source
//! stiffness may change in this ORIGINAL fixed coordinate basis. Physical
//! reference mass is assembled directly from the original director velocity
//! field; numerical drilling inertia in a regularized modal pencil is excluded.
//! This lane does not claim finite rotation, prestressed tangents, damping,
//! moving geometry, or temperature-dependent reference inertia.
//! Director-axis drilling has no physical inertia. Its algebraic constraint
//! must hold for initial coordinates and rates; unsupported director-axis
//! efforts are refused instead of supplying fictitious gauge mass.

use super::{PreparedThermalShell, ThermalShellResponse};
use crate::PlateError;
use crate::shell::mitc::element::PreparedStrainField;
use crate::shell::mitc::static_solver::StaticFactor;
use crate::shell::mitc::strain::RealElementStrain;
use crate::shell::mitc::{MitcShellModel, ShellMass, dd_add, two_prod};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;
use std::sync::Arc;

mod material;
mod step;

/// An immutable full enriched pencil sharing the original nodal shell model.
#[derive(Clone, Debug)]
pub struct EnrichedMitcShell {
    // Permanent coordinate/mass identity carried by every retained state.
    model: Arc<MitcShellModel>,
    material_model: Arc<MitcShellModel>,
    reference_fields: Arc<Vec<PreparedStrainField>>,
    reference_stiffness: Arc<Csr>,
    gauge_stiffness: Arc<Csr>,
    stiffness: Arc<Csr>,
    mass: Arc<Csr>,
    gauges: Arc<Vec<AlgebraicGauge>>,
}

#[derive(Debug)]
struct AlgebraicGauge {
    direction: Vec<(usize, f64)>,
    // (column, z^T K coefficient, absolute pre-cancellation row bound).
    constraint: Vec<(usize, f64, f64)>,
}

/// Accepted or proposed full mechanical state. Values cannot be changed
/// independently of their original operator, physical support map or energy.
#[derive(Clone, Debug)]
pub struct MitcTransientState {
    model: Arc<MitcShellModel>,
    time_s: f64,
    coordinates: Vec<f64>,
    velocity: Vec<f64>,
    nodal_displacement_6n: Vec<f64>,
    nodal_velocity_6n: Vec<f64>,
    bubble_coordinates: Vec<[f64; 2]>,
    bubble_velocity: Vec<[f64; 2]>,
    kinetic_energy_j: f64,
}

impl MitcTransientState {
    /// Physical endpoint time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Free nodal support coordinates followed by two residual rotations per facet.
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }
    /// Time derivatives in exactly the same generalized coordinate ordering.
    pub fn velocity(&self) -> &[f64] {
        &self.velocity
    }
    /// World displacement and physical nodal rotation, six entries per node.
    pub fn nodal_displacement_6n(&self) -> &[f64] {
        &self.nodal_displacement_6n
    }
    /// World translation and nodal rotation rates, six entries per node.
    pub fn nodal_velocity_6n(&self) -> &[f64] {
        &self.nodal_velocity_6n
    }
    /// Independent residual rotation r, added to the original elastic map Tq.
    pub fn bubble_coordinates(&self) -> &[[f64; 2]] {
        &self.bubble_coordinates
    }
    /// Rates of the independent residual rotations, not total director rates.
    pub fn bubble_velocity(&self) -> &[[f64; 2]] {
        &self.bubble_velocity
    }
    /// Complete consistent kinetic energy, including nodal/bubble cross inertia.
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }
}

/// Physical response and conjugate generalized forces at the actual retained
/// director coordinates. The old static thermal bubble is never substituted.
pub struct EnrichedThermalResponse {
    /// Facet resultants and storage. Nodal forces are conjugate to q in Tq+Br.
    pub response: ThermalShellResponse,
    /// Exact native strain polynomials, including the independent rotations.
    pub strain_fields: Vec<RealElementStrain>,
    /// Material forces conjugate to each original facet's two residual rotations.
    pub bubble_material_force: Vec<[f64; 2]>,
    /// Complete material force in the enriched reduced coordinate ordering.
    pub material_generalized_force: Vec<f64>,
    /// Original numerical drilling force, separated from physical material work.
    pub numerical_gauge_force: Vec<f64>,
    /// Numerical drilling storage; this is never a thermal source.
    pub drilling_gauge_energy_j: f64,
    /// Conservative whole-facet bound on actual director motion, including
    /// independent interior rotations. Nodal rotations alone are insufficient.
    pub maximum_director_motion_rad: f64,
}

/// A fixed-step midpoint factor retaining shared immutable operators. It may
/// be cached across accepted intervals with the same step size and shell.
pub struct PreparedMitcMidpoint {
    shell: EnrichedMitcShell,
    dt_s: f64,
    effective: Csr,
    factor: StaticFactor,
}

/// The linear mechanical first-law identity, including explicit thermal loads.
/// This is not the coupled caloric/internal-energy account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MitcTransientReport {
    /// Complete consistent kinetic-energy difference [J].
    pub kinetic_energy_change_j: f64,
    /// Endpoint-average TOTAL generalized load times actual displacement [J].
    pub load_work_j: f64,
    /// Change in one-half q^T K q, before thermal eigenstrain subtraction [J].
    pub quadratic_energy_change_j: f64,
    /// Polarized work from changing the physical stiffness at fixed reference
    /// coordinates: q0^T (K1-K0) q1 / 2 [J]. It is not a caloric heat source.
    pub stiffness_parameter_work_j: f64,
    /// Kinetic plus quadratic storage change, minus stiffness-parameter work
    /// and total load work [J].
    pub energy_residual_j: f64,
}

/// An immutable proposal; callers publish it only after their coupled gates pass.
pub struct MitcTransientStep {
    /// Proposed mechanical endpoint.
    pub state: MitcTransientState,
    /// Mechanical midpoint balance on this interval.
    pub report: MitcTransientReport,
}

impl EnrichedMitcShell {
    /// Restore the original two interior director DOFs per facet. The original
    /// nodal stiffness and drilling gauge are retained. Physical consistent
    /// mass is assembled from the native translation/director velocity field,
    /// excluding any numerical drilling inertia from the modal pencil.
    /// Consistent native mass is required; lumping explicitly removes interior
    /// inertia and therefore cannot supply these two momentum equations.
    pub fn try_new(
        model: impl Into<Arc<MitcShellModel>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        check(&mut cancelled)?;
        let model = model.into();
        if model.mass_model != ShellMass::Consistent {
            return Err(invalid(
                "retained thermal director dynamics requires native consistent mass",
            ));
        }
        let nodal = model.pencil.free;
        let dimension = model
            .mesh
            .tris
            .len()
            .checked_mul(2)
            .and_then(|n| nodal.checked_add(n))
            .ok_or_else(|| invalid("enriched shell coordinate extent overflow"))?;
        let mut k = Coo::new(dimension, dimension);
        let mut m = Coo::new(dimension, dimension);
        let mut gauge = Coo::new(nodal, nodal);
        let mut reference_fields = Vec::with_capacity(model.mesh.tris.len());
        for row in 0..nodal {
            row_check(row, &mut cancelled)?;
            let (cols, vals) = model.pencil.k.row(row);
            for (&col, &val) in cols.iter().zip(vals) {
                k.push(row, col, val);
                gauge.push(row, col, val);
            }
        }
        for (e, tri) in model.mesh.tris.iter().enumerate() {
            check(&mut cancelled)?;
            let field = PreparedStrainField::new(
                &model.mesh,
                tri,
                model.mesh.facet(e)?,
                tri.map(|n| model.directors.vectors[n]),
                &model.sections[e],
                &mut cancelled,
            )?;
            let inertia = field.physical_inertia(&model.sections[e], &mut cancelled)?;
            material::scatter(&mut m, &model, e, &inertia, 1.0, &mut cancelled)?;
            let physical = field.condensed_stiffness(&model.sections[e], &mut cancelled)?;
            material::scatter(&mut gauge, &model, e, &physical, -1.0, &mut cancelled)?;
            let kb = field.bubble_stiffness();
            let first = nodal + 2 * e;
            for i in 0..2 {
                for j in 0..2 {
                    k.push(first + i, first + j, kb[i][j]);
                }
            }
            reference_fields.push(field);
        }
        check(&mut cancelled)?;
        let stiffness = Arc::new(k.assemble());
        let mass = Arc::new(m.assemble());
        let gauge_stiffness = Arc::new(gauge.assemble());
        let gauges = Arc::new(algebraic_gauges(&model, &stiffness, &mut cancelled)?);
        check(&mut cancelled)?;
        Ok(Self {
            material_model: Arc::clone(&model),
            model,
            reference_fields: Arc::new(reference_fields),
            reference_stiffness: Arc::clone(&stiffness),
            gauge_stiffness,
            stiffness,
            mass,
            gauges,
        })
    }

    /// Current immutable material model, with the original geometry/supports.
    /// Thermal loading must be prepared on this model; retained states keep
    /// their original reference coordinate owner even when this model changes.
    pub fn model(&self) -> &MitcShellModel {
        &self.material_model
    }
    /// Share the current material model without copying its mesh or sparse pencil.
    pub fn shared_model(&self) -> Arc<MitcShellModel> {
        Arc::clone(&self.material_model)
    }
    /// Enriched stiffness, including the original nodal drilling gauge.
    pub fn stiffness(&self) -> &Csr {
        &self.stiffness
    }
    /// Complete native consistent inertia in the enriched coordinate basis.
    pub fn mass(&self) -> &Csr {
        &self.mass
    }
    /// Number of free nodal coordinates plus two retained rotations per facet.
    pub fn dimension(&self) -> usize {
        self.mass.nrows()
    }

    /// Admit an explicit complete state. Nodal values must already obey the
    /// original supports; projecting a forbidden displacement or velocity away
    /// would destroy its physical initial condition and is refused.
    #[allow(clippy::too_many_arguments)]
    pub fn state(
        &self,
        time_s: f64,
        nodal_displacement_6n: &[f64],
        nodal_velocity_6n: &[f64],
        bubble_coordinates: &[[f64; 2]],
        bubble_velocity: &[[f64; 2]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MitcTransientState, PlateError> {
        check(&mut cancelled)?;
        let mut q = admit_nodal(&self.model, nodal_displacement_6n)?;
        let mut v = admit_nodal(&self.model, nodal_velocity_6n)?;
        let count = self.model.mesh.tris.len();
        if bubble_coordinates.len() != count || bubble_velocity.len() != count {
            return Err(invalid(
                "two independent director coordinates and rates per original facet required",
            ));
        }
        q.extend(bubble_coordinates.iter().flatten().copied());
        v.extend(bubble_velocity.iter().flatten().copied());
        self.state_reduced(time_s, q, v, &mut cancelled)
    }

    /// Seed the thermal static bubble correction at a supplied nodal state.
    /// All velocities start at zero. The caller owns the nodal equilibrium
    /// solve; this constructor neither solves it nor overwrites its result.
    pub fn equilibrium_state(
        &self,
        time_s: f64,
        thermal: &PreparedThermalShell<'_>,
        nodal_displacement_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MitcTransientState, PlateError> {
        self.admit_thermal(thermal)?;
        check(&mut cancelled)?;
        admit_nodal(&self.model, nodal_displacement_6n)?;
        let mut bubbles = Vec::with_capacity(thermal.facets.len());
        for (e, (facet, tri)) in thermal.facets.iter().zip(&self.model.mesh.tris).enumerate() {
            check(&mut cancelled)?;
            let nodal = core::array::from_fn(|i| nodal_displacement_6n[6 * tri[i / 6] + i % 6]);
            bubbles.push(facet.equilibrium_bubble_in_basis(&self.reference_fields[e], &nodal)?);
        }
        self.state(
            time_s,
            nodal_displacement_6n,
            &vec![0.0; nodal_displacement_6n.len()],
            &bubbles,
            &vec![[0.0; 2]; bubbles.len()],
            &mut cancelled,
        )
    }

    /// Thermal load plus original physical nodal loads, in work-conjugate
    /// support/residual-bubble coordinates. Loads at constrained nodes do no work.
    pub fn generalized_load(
        &self,
        thermal: &PreparedThermalShell<'_>,
        external_load_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        check(&mut cancelled)?;
        self.admit_thermal(thermal)?;
        self.model.project_load(external_load_6n)?;
        let mut world = external_load_6n.to_vec();
        for (e, (facet, tri)) in thermal.facets.iter().zip(&self.model.mesh.tris).enumerate() {
            check(&mut cancelled)?;
            let load = facet.load_in_basis(&self.reference_fields[e])?;
            for i in 0..18 {
                let slot = 6 * tri[i / 6] + i % 6;
                world[slot] = finite(world[slot] + load[i])?;
            }
        }
        let mut out = self.model.project_load(&world)?;
        for facet in &thermal.facets {
            check(&mut cancelled)?;
            out.extend(facet.bubble_load());
        }
        self.admit_gauge_load(&out)?;
        check(&mut cancelled)?;
        Ok(out)
    }

    /// Recover source thermal stress at the actual independent director motion.
    pub fn recover(
        &self,
        thermal: &PreparedThermalShell<'_>,
        state: &MitcTransientState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnrichedThermalResponse, PlateError> {
        check(&mut cancelled)?;
        self.admit_thermal(thermal)?;
        self.admit_state(state)?;
        let mut facets = Vec::with_capacity(thermal.facets.len());
        let mut strain_fields = Vec::with_capacity(thermal.facets.len());
        let mut bubble_material_force = Vec::with_capacity(thermal.facets.len());
        let mut nodal_force = vec![0.0; state.nodal_displacement_6n.len()];
        let mut storage = (0.0, 0.0);
        let mut maximum_director_motion_rad = 0.0_f64;
        for (e, (prepared, tri)) in thermal.facets.iter().zip(&self.model.mesh.tris).enumerate() {
            check(&mut cancelled)?;
            let nodal =
                core::array::from_fn(|i| state.nodal_displacement_6n[6 * tri[i / 6] + i % 6]);
            let bubble = state.bubble_coordinates[e];
            let basis = &self.reference_fields[e];
            maximum_director_motion_rad = maximum_director_motion_rad.max(
                prepared.director_motion_bound_in_basis(basis, &nodal, bubble, &mut cancelled)?,
            );
            let (facet, forces) =
                prepared.recover_in_basis(basis, &nodal, bubble, &mut cancelled)?;
            strain_fields.push(prepared.strain_field_in_basis(
                basis,
                &nodal,
                bubble,
                &mut cancelled,
            )?);
            storage = dd_add(storage, (facet.strain_energy_j, 0.0));
            for i in 0..18 {
                nodal_force[6 * tri[i / 6] + i % 6] += facet.material_nodal_force[i];
            }
            facets.push(facet);
            bubble_material_force.push(forces);
        }
        let mut material_generalized_force = self.model.project_load(&nodal_force)?;
        material_generalized_force.extend(bubble_material_force.iter().flatten().copied());
        // The numerical gauge is the ORIGINAL fixed nodal operator. In
        // particular, no subtraction of temperature-dependent material forces
        // can manufacture an internal-coordinate gauge force or heat source.
        let mut numerical_gauge_force = action(
            &self.gauge_stiffness,
            &state.coordinates[..self.model.pencil.free],
            &mut cancelled,
        )?;
        numerical_gauge_force.resize(self.dimension(), 0.0);
        let drilling_gauge_energy_j = 0.5 * product(&state.coordinates, &numerical_gauge_force)?;
        let strain_energy_j = finite(storage.0 + storage.1)?;
        let response = ThermalShellResponse {
            displacement_6n: state.nodal_displacement_6n.clone(),
            facets,
            strain_energy_j,
            material_nodal_force_6n: nodal_force,
        };
        check(&mut cancelled)?;
        Ok(EnrichedThermalResponse {
            response,
            strain_fields,
            bubble_material_force,
            material_generalized_force,
            numerical_gauge_force,
            drilling_gauge_energy_j,
            maximum_director_motion_rad,
        })
    }

    /// Factor M + (dt/2)^2 K once, sharing the static owner's equilibrated
    /// sparse LDLT, compensated refinement and true-residual gate. Factorization
    /// and individual sparse kernels remain indivisible cancellation regions.
    pub fn prepare_midpoint(
        &self,
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMitcMidpoint, PlateError> {
        check(&mut cancelled)?;
        let h = 0.5 * dt_s;
        let h2 = h * h;
        if !dt_s.is_finite() || dt_s <= 0.0 || !h2.is_finite() || h2 <= 0.0 {
            return Err(invalid("positive representable midpoint interval required"));
        }
        let n = self.dimension();
        let mut coo = Coo::new(n, n);
        for row in 0..n {
            row_check(row, &mut cancelled)?;
            for (matrix, scale) in [(&self.mass, 1.0), (&self.stiffness, h2)] {
                let (cols, vals) = matrix.row(row);
                for (&col, &val) in cols.iter().zip(vals) {
                    coo.push(row, col, finite(scale * val)?);
                }
            }
        }
        check(&mut cancelled)?;
        let effective = coo.assemble();
        let factor = StaticFactor::try_new(&effective, &mut cancelled)?;
        check(&mut cancelled)?;
        Ok(PreparedMitcMidpoint {
            shell: self.clone(),
            dt_s,
            effective,
            factor,
        })
    }

    fn admit_thermal(&self, thermal: &PreparedThermalShell<'_>) -> Result<(), PlateError> {
        if !std::ptr::eq(self.material_model.as_ref(), thermal.model) {
            return Err(invalid(
                "thermal loading must belong to the exact original enriched shell",
            ));
        }
        Ok(())
    }
    fn admit_state(&self, state: &MitcTransientState) -> Result<(), PlateError> {
        if !Arc::ptr_eq(&self.model, &state.model) || state.coordinates.len() != self.dimension() {
            return Err(invalid(
                "mechanical state must retain the exact original enriched shell",
            ));
        }
        Ok(())
    }
    fn admit_gauge_load(&self, load: &[f64]) -> Result<(), PlateError> {
        for gauge in self.gauges.iter() {
            let mut value = 0.0;
            let mut scale = 0.0;
            for &(row, direction) in &gauge.direction {
                let term = direction * load[row];
                value += term;
                // Include tangential effort in the scale: support-frame and
                // world-load projection can leave a roundoff normal component.
                scale += load[row].abs();
            }
            if !value.is_finite() || !scale.is_finite() || value.abs() > 1e-8 * scale {
                return Err(invalid(
                    "director-axis effort is unsupported by the massless drilling gauge",
                ));
            }
        }
        Ok(())
    }
    fn admit_gauge_state(
        &self,
        q: &[f64],
        v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PlateError> {
        for gauge in self.gauges.iter() {
            check(cancelled)?;
            for state in [q, v] {
                let (mut value, mut scale, mut roundoff) = (0.0, 0.0, 0.0);
                for &(column, coefficient, absolute) in &gauge.constraint {
                    let term = coefficient * state[column];
                    value += term;
                    scale += term.abs();
                    roundoff += absolute * state[column].abs();
                }
                let tolerance = 1e-8 * scale + 256.0 * f64::EPSILON * roundoff;
                if !value.is_finite() || !tolerance.is_finite() || value.abs() > tolerance {
                    return Err(invalid(
                        "shell coordinates and rates must satisfy the massless drilling constraint",
                    ));
                }
            }
        }
        Ok(())
    }
    fn state_reduced(
        &self,
        time_s: f64,
        q: Vec<f64>,
        v: Vec<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MitcTransientState, PlateError> {
        check(cancelled)?;
        if !time_s.is_finite()
            || q.len() != self.dimension()
            || v.len() != self.dimension()
            || q.iter().chain(&v).any(|x| !x.is_finite())
        {
            return Err(invalid(
                "finite complete enriched mechanical state required",
            ));
        }
        self.admit_gauge_state(&q, &v, cancelled)?;
        let nodal = self.model.pencil.free;
        let nodal_displacement_6n = self.model.expand(&q[..nodal])?;
        let nodal_velocity_6n = self.model.expand(&v[..nodal])?;
        let bubble_coordinates = q[nodal..].chunks_exact(2).map(|r| [r[0], r[1]]).collect();
        let bubble_velocity = v[nodal..].chunks_exact(2).map(|r| [r[0], r[1]]).collect();
        let mv = action(&self.mass, &v, cancelled)?;
        let kinetic_energy_j = 0.5 * product(&v, &mv)?;
        if kinetic_energy_j < 0.0 {
            return Err(invalid("negative enriched kinetic energy"));
        }
        check(cancelled)?;
        Ok(MitcTransientState {
            model: Arc::clone(&self.model),
            time_s,
            coordinates: q,
            velocity: v,
            nodal_displacement_6n,
            nodal_velocity_6n,
            bubble_coordinates,
            bubble_velocity,
            kinetic_energy_j,
        })
    }
}

impl PreparedMitcMidpoint {
    /// Fixed positive step size used by this immutable effective factor.
    pub const fn dt_s(&self) -> f64 {
        self.dt_s
    }

    /// Solve a full enriched implicit-midpoint interval. Endpoint loads include
    /// thermal and external forces in the exact generalized ordering. This
    /// method never mutates the initial state or a caller's accepted timeline.
    pub fn advance(
        &self,
        initial: &MitcTransientState,
        load_before: &[f64],
        load_after: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MitcTransientStep, PlateError> {
        self.advance_from(&self.shell, initial, load_before, load_after, cancelled)
    }
}

fn algebraic_gauges(
    model: &MitcShellModel,
    stiffness: &Csr,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<AlgebraicGauge>, PlateError> {
    let mut gauges = Vec::new();
    for (node, director) in model.directors.vectors.iter().enumerate() {
        check(cancelled)?;
        let mut direction = Vec::with_capacity(3);
        let (mut free_norm, mut fixed_projection) = (0.0, 0.0_f64);
        for axis in 0..3 {
            let value: f64 = (0..3)
                .map(|j| model.frames[node][axis][j] * director[j])
                .sum();
            if let Some(row) = model.pencil.dof_map[6 * node + 3 + axis] {
                direction.push((row, value));
                free_norm += value * value;
            } else {
                fixed_projection = fixed_projection.max(value.abs());
            }
        }
        // Only a complete free director-axis rotation is a mass-null mode.
        // Partially constrained rotation subspaces retain physical inertia.
        if fixed_projection > 4e-12 || (free_norm - 1.0).abs() > 8e-12 {
            continue;
        }
        let mut row = BTreeMap::<usize, (f64, f64)>::new();
        for &(index, value) in &direction {
            let (cols, vals) = stiffness.row(index);
            for (&column, &coefficient) in cols.iter().zip(vals) {
                let term = value * coefficient;
                let entry = row.entry(column).or_default();
                entry.0 += term;
                entry.1 += term.abs();
            }
        }
        let constraint = row
            .into_iter()
            .map(|(column, (coefficient, absolute))| (column, coefficient, absolute))
            .collect();
        gauges.push(AlgebraicGauge {
            direction,
            constraint,
        });
    }
    Ok(gauges)
}

fn admit_nodal(model: &MitcShellModel, world: &[f64]) -> Result<Vec<f64>, PlateError> {
    let reduced = model.project_load(world)?;
    let expanded = model.expand(&reduced)?;
    for (a, b) in world.chunks_exact(3).zip(expanded.chunks_exact(3)) {
        let scale = a.iter().chain(b).fold(0.0_f64, |m, x| m.max(x.abs()));
        if a.iter().zip(b).any(|(a, b)| (a - b).abs() > 4e-12 * scale) {
            return Err(invalid(
                "initial shell motion violates the original homogeneous supports",
            ));
        }
    }
    Ok(reduced)
}

fn action(
    matrix: &Csr,
    x: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, PlateError> {
    let mut out = vec![0.0; matrix.nrows()];
    for (row, value) in out.iter_mut().enumerate() {
        row_check(row, cancelled)?;
        let (cols, vals) = matrix.row(row);
        let mut sum = (0.0, 0.0);
        for (&col, &a) in cols.iter().zip(vals) {
            sum = dd_add(sum, two_prod(a, x[col]));
        }
        *value = finite(sum.0 + sum.1)?;
    }
    Ok(out)
}
fn product(a: &[f64], b: &[f64]) -> Result<f64, PlateError> {
    let mut sum = (0.0, 0.0);
    for (&a, &b) in a.iter().zip(b) {
        sum = dd_add(sum, two_prod(a, b));
    }
    finite(sum.0 + sum.1)
}
fn finite(value: f64) -> Result<f64, PlateError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable enriched shell state or energy"))
    }
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}
fn row_check(row: usize, cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if row % 256 == 0 {
        check(cancelled)
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}
