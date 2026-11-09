//! Frozen native initial stress with an inseparable displacement reference.
use fs_cutfem::elastic3::{ElasticityError3, PreparedPrestress3};

use super::{CutElasticityRef3, CutSecondOrderProblem3, checkpoint, invalid, poll};

pub(super) struct HeldPrestress3 {
    pub(super) prepared: PreparedPrestress3,
    /// Original independent coordinates, with exactly zero strong support rows.
    pub(super) reference: Vec<f64>,
}

impl CutSecondOrderProblem3<'_> {
    /// Freeze signed initial-stress stiffness at an admitted held-eigenstrain
    /// reference. The reference is a full native displacement in meters;
    /// adaptive values are independent master coordinates. Strong support
    /// entries must be exactly zero. Use [`Self::expand_field`] to expand an
    /// already reduced reference.
    ///
    /// This builder requires [`Self::new_with_eigenstrain`]. It prepares G from
    /// that exact eigenstrain, material, scales, geometry and reference; callers
    /// cannot accidentally pair a geometric cache with another displacement.
    /// The original preparation owner admits actual/projected elastic strains,
    /// total small-deformation limits and selected embedded support points.
    /// `maximum_prestress_points` bounds bulk PLUS selected support evaluations;
    /// `maximum_geometric_entries` bounds scalar pairs before xyz replication.
    /// Two additional native-sized vectors hold the reference and action scratch.
    ///
    /// The equation becomes `M a + K q - b_free + G(q-q_0) = f_external`, in
    /// absolute displacement coordinates. G is applied to the displacement
    /// difference directly, so the reference gains no spurious geometric force.
    /// Tangents and state/force adjoints use K+G with q_0 and G held fixed.
    /// The reference's equilibrium under the intended external loads remains
    /// caller-owned; preparation establishes material admission, not equilibrium.
    ///
    /// The original consistent mass, occupancy and time integrator are unchanged.
    /// Compression can make G or K+G indefinite and is never clipped. This is
    /// the native frozen initial-stress approximation, not a finite-deformation
    /// tangent, follower-load correction, evolving prestress or stability proof.
    /// Failed preparation or cancellation publishes no successor adapter.
    pub fn with_reference_prestress(
        mut self,
        native_reference_displacement: &[f64],
        maximum_prestress_points: usize,
        maximum_geometric_entries: usize,
    ) -> Result<Self, ElasticityError3> {
        self.admit_vector(native_reference_displacement, self.mass.n())?;
        let held = self
            .eigenstrain
            .as_ref()
            .ok_or(invalid("CutFEM prestress dynamics requires held eigenstrain"))?;
        for (dof, value) in native_reference_displacement.iter().enumerate() {
            if dof % 192 == 0 {
                poll(self.gate)?;
            }
            if self.source.fixed()[dof / 3] && *value != 0.0 {
                return Err(invalid("CutFEM prestress reference moves a strong support"));
            }
        }
        let prepared = match self.source {
            CutElasticityRef3::Cartesian(op) => {
                let prepared = op.prepare_eigenstrain_prestress(
                    held.prepared,
                    native_reference_displacement,
                    held.maximum_kinematic_strain,
                    maximum_prestress_points,
                    maximum_geometric_entries,
                    || checkpoint(self.gate),
                )?;
                op.validate_prestress(&prepared)?;
                prepared
            }
            CutElasticityRef3::Adaptive(op) => {
                let prepared = op.prepare_eigenstrain_prestress(
                    held.prepared,
                    native_reference_displacement,
                    held.maximum_kinematic_strain,
                    maximum_prestress_points,
                    maximum_geometric_entries,
                    || checkpoint(self.gate),
                )?;
                op.validate_prestress(&prepared)?;
                prepared
            }
        };
        let allocation = || invalid("CutFEM prestress dynamics allocation failed");
        let mut reference = Vec::new();
        reference
            .try_reserve_exact(native_reference_displacement.len())
            .map_err(|_| allocation())?;
        reference.extend_from_slice(native_reference_displacement);
        let mut geometric = Vec::new();
        geometric
            .try_reserve_exact(reference.len())
            .map_err(|_| allocation())?;
        geometric.resize(reference.len(), 0.0);
        poll(self.gate)?;
        self.scratch.get_mut().geometric = geometric;
        self.prestress = Some(HeldPrestress3 {
            prepared,
            reference,
        });
        Ok(self)
    }

    /// Borrow the exact signed geometric operator used by this trajectory.
    /// This permits native modal queries to reuse the same reference preparation.
    #[must_use]
    pub fn reference_prestress(&self) -> Option<&PreparedPrestress3> {
        self.prestress.as_ref().map(|held| &held.prepared)
    }
}
