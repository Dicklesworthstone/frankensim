//! Prepared tensor-memory mechanics on the original solid and force map.
//!
//! Every quadrature contribution remains attached to its original element.
//! Equal relaxation times share exact displacement-history coordinates in the
//! existing material runtime; this adapter neither changes the constitutive
//! equations nor removes a displacement coordinate from the supplied basis.

use crate::tensor_material::InitialMemory;
use crate::tensor_material::prepared::{
    PreparedTensorBody, PreparedTensorControls, PreparedTensorError, PreparedTensorHistory,
    PreparedTensorStep,
};
use fs_exec::CancelGate;
use fs_solid::linear3::visco::{ProjectedViscoSolid, ViscoElementError};

/// Source or physical-runtime refusal before accepted nodal motion changes.
#[derive(Debug)]
pub enum PreparedSolidError {
    Source(ViscoElementError),
    Dynamics(PreparedTensorError),
    Invalid(&'static str),
    /// Original nodal work and retained-coordinate work must agree [J].
    Work {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for PreparedSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prepared tensor solid refused: {self:?}")
    }
}
impl std::error::Error for PreparedSolidError {}
impl From<ViscoElementError> for PreparedSolidError {
    fn from(error: ViscoElementError) -> Self {
        Self::Source(error)
    }
}
impl From<PreparedTensorError> for PreparedSolidError {
    fn from(error: PreparedTensorError) -> Self {
        Self::Dynamics(error)
    }
}
type Result<T> = core::result::Result<T, PreparedSolidError>;

fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(PreparedSolidError::Cancelled)
    } else {
        Ok(())
    }
}

/// Fixed original mesh projection with a reusable exact tensor-memory flow.
///
/// A full identity basis retains every represented nodal coordinate. Other
/// bases keep the original projection owner's explicit approximation scope.
/// This owner fixes geometry, material coefficients, support map and clock;
/// changing them requires a separately admitted history-preserving transfer.
pub struct PreparedTensorSolidBody<'s> {
    projection: ProjectedViscoSolid<'s>,
    body: PreparedTensorBody,
    work_allowance_j: f64,
}

impl<'s> PreparedTensorSolidBody<'s> {
    /// Group independently integrated positive material loss by original
    /// element, avoiding one dense observation cache per quadrature point.
    /// The source projection's own storage/work admission remains in force;
    /// `controls` bounds the new prepared material runtime and its caches.
    pub fn new(
        projection: ProjectedViscoSolid<'s>,
        dt_s: f64,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        if projection.coordinates() > controls.maximum_coordinates
            || projection.material_points().len() > controls.maximum_material_points
            || projection.element_count() > controls.maximum_heat_groups
        {
            return Err(PreparedSolidError::Invalid(
                "original solid exceeds prepared coordinate, point or element capacity",
            ));
        }
        let body = PreparedTensorBody::from_mass_matrix(
            projection.coordinates(),
            projection.mass_matrix_kg().to_vec(),
            projection.material_points().to_vec(),
            projection.point_elements().to_vec(),
            dt_s,
            controls,
            gate,
        )?;
        if body.heat_group_count() != projection.element_count() {
            return Err(PreparedSolidError::Invalid(
                "each original element must retain its material heat channel",
            ));
        }
        poll(gate)?;
        Ok(Self {
            projection,
            body,
            work_allowance_j: controls.maximum_work_error_j,
        })
    }

    pub fn projection(&self) -> &ProjectedViscoSolid<'s> {
        &self.projection
    }
    pub fn material_body(&self) -> &PreparedTensorBody {
        &self.body
    }

    /// Initial motion and compatible constitutive prehistory are explicit.
    pub fn history<'a>(
        &'a self,
        displacement_m: &[f64],
        velocity_mps: &[f64],
        memory: InitialMemory,
        gate: &CancelGate,
    ) -> Result<PreparedTensorSolidHistory<'a, 's>> {
        let inner = self
            .body
            .history(displacement_m, velocity_mps, memory, gate)?;
        Ok(PreparedTensorSolidHistory { solid: self, inner })
    }
}

/// Accepted material step and its original nodal power-conjugate observation.
#[derive(Clone, Debug)]
pub struct PreparedTensorSolidStep {
    pub material: PreparedTensorStep,
    pub mean_nodal_velocity_mps: Vec<[f64; 3]>,
    /// Nodal held-force work minus retained-coordinate work [J].
    pub nodal_work_error_j: f64,
}
impl PreparedTensorSolidStep {
    /// Physical Kelvin plus Maxwell heat in original element order [J].
    pub fn element_heat_j(&self) -> &[f64] {
        &self.material.heat_by_group_j
    }
}

/// Accepted state borrows one immutable source projection and prepared cache.
#[derive(Clone)]
pub struct PreparedTensorSolidHistory<'a, 's> {
    solid: &'a PreparedTensorSolidBody<'s>,
    inner: PreparedTensorHistory<'a>,
}
impl PreparedTensorSolidHistory<'_, '_> {
    pub fn material_history(&self) -> &PreparedTensorHistory<'_> {
        &self.inner
    }
    pub fn nodal_displacement_m(&self, gate: &CancelGate) -> Result<Vec<[f64; 3]>> {
        Ok(self
            .solid
            .projection
            .expand_nodal(self.inner.displacement_m(), gate)?)
    }
    pub fn nodal_velocity_mps(&self, gate: &CancelGate) -> Result<Vec<[f64; 3]>> {
        Ok(self
            .solid
            .projection
            .expand_nodal(self.inner.velocity_mps(), gate)?)
    }

    /// Advance under a complete held nodal-force field. Only O(state) accepted
    /// history is forked for the candidate; no material or propagation cache is
    /// cloned. Expansion and original nodal-work checks precede publication.
    pub fn advance(
        &mut self,
        nodal_force_n: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<PreparedTensorSolidStep> {
        poll(gate)?;
        let force = self
            .solid
            .projection
            .project_nodal_force(nodal_force_n, gate)?;
        let mut candidate = self.inner.clone();
        let material = candidate.advance(&force, gate)?;
        let velocity = self
            .solid
            .projection
            .expand_nodal(&material.mean_velocity_mps, gate)?;
        let dt = self.solid.body.duration_s();
        let (mut work, mut correction) = (0.0, 0.0);
        for (force, velocity) in nodal_force_n.iter().zip(&velocity) {
            poll(gate)?;
            for axis in 0..3 {
                let term = force[axis] * (dt * velocity[axis]);
                let increment = term - correction;
                let next = work + increment;
                correction = (next - work) - increment;
                work = next;
            }
        }
        let defect_j = work - material.energy.supplied_work_j;
        if !defect_j.is_finite() || defect_j.abs() > self.solid.work_allowance_j {
            return Err(PreparedSolidError::Work {
                defect_j,
                allowance_j: self.solid.work_allowance_j,
            });
        }
        poll(gate)?;
        self.inner = candidate;
        Ok(PreparedTensorSolidStep {
            material,
            mean_nodal_velocity_mps: velocity,
            nodal_work_error_j: defect_j,
        })
    }
}
