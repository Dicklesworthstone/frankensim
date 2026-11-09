//! Complete tensor material on the native corotated spatial connector port.
//!
//! Explicit tensor-shear strain columns and positive quadrature volume map the
//! original equilibrium, Maxwell and Kelvin tensors to work-conjugate SI
//! wrenches. Every original factor row and relaxation time survives. The
//! equilibrium channel is lossless; Maxwell/Kelvin loss alone becomes heat.
//!
//! Coordinates are integrals of the native first-body-frame relative twist.
//! This is a small-strain corotated rate law, not a finite Log(R) spring or a
//! rest geometry inferred from the current pose. All six source coordinates
//! must continue independently of the energy amplitudes: a semidefinite map
//! cannot recover strain history from a memory bank.
//!
//! No plane-stress approximation, density-derived stiffness, fitted damping,
//! thermal expansion or temperature-dependent stored tensor is introduced.

use fs_material::{
    Voigt,
    visco_tensor::{
        ProjectedTensorProny, TensorProjectionBudget, TensorProny, TensorPronyError,
        source::{ResolvedTensorMaterial, TensorSourceError},
    },
};
use fs_math::dd::Dd;
use fs_mbd::{
    Vec3,
    articulated::motion::spatial_loads::bushing::coupled::{
        BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix,
    },
    spatial_loads::{
        internal_impulse::viscous::network::{
            ViscousNetworkLink,
            relaxation::{
                RelaxationNetworkError, RelaxationNetworkLink, RelaxationNetworkReceipt,
                SpatialRelaxationBranch,
            },
        },
        internal_motion::viscous::network::relaxation::RelaxationNetworkMotionStep,
    },
};

/// Original geometry and six work-conjugate strain coordinates.
///
/// Columns are tensor-shear Voigt (11,22,33,12,23,31), not engineering shear.
/// Columns 0..3 multiply translation [m], columns 3..6 rotation [rad]. Their
/// strain components and the material tensors are expressed in the first
/// body's material frame. Volume includes the physical quadrature weight.
/// Several contributions may share bodies, retaining separate original heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorSpatialGeometry {
    pub bodies: [usize; 2],
    pub arms_body_m: [Vec3; 2],
    pub strain_columns: [Voigt; 6],
    pub volume_m3: f64,
    pub coordinate_scales: BushingCoordinateScales,
}

/// Inadmissible material, geometry, continuation or original native source.
#[derive(Debug)]
pub enum TensorSpatialError {
    Invalid(&'static str),
    Material(TensorPronyError),
    Source(TensorSourceError),
    Matrix(BushingMatrixError),
    Native(RelaxationNetworkError),
    Cancelled,
}
impl core::fmt::Display for TensorSpatialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tensor spatial material refused: {self:?}")
    }
}
impl std::error::Error for TensorSpatialError {}
impl From<TensorPronyError> for TensorSpatialError {
    fn from(error: TensorPronyError) -> Self {
        match error {
            TensorPronyError::Cancelled => Self::Cancelled,
            error => Self::Material(error),
        }
    }
}
impl From<TensorSourceError> for TensorSpatialError {
    fn from(error: TensorSourceError) -> Self {
        Self::Source(error)
    }
}
impl From<BushingMatrixError> for TensorSpatialError {
    fn from(error: BushingMatrixError) -> Self {
        Self::Matrix(error)
    }
}
impl From<RelaxationNetworkError> for TensorSpatialError {
    fn from(error: RelaxationNetworkError) -> Self {
        Self::Native(error)
    }
}
type Result<T> = core::result::Result<T, TensorSpatialError>;

/// One complete material quadrature contribution on the original native graph.
///
/// Native channel zero is the original lossless equilibrium factor; channel
/// j+1 is original Maxwell branch j. Kelvin viscosity stays memoryless.
/// The immutable source, declared rotation and complete projection are retained.
/// This binding does not add material mass to either original rigid body.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorSpatialRelaxation<'a> {
    geometry: TensorSpatialGeometry,
    projection: ProjectedTensorProny,
    native: RelaxationNetworkLink,
    source: Option<&'a ResolvedTensorMaterial>,
    material_to_first_body: Option<[[f64; 3]; 3]>,
}
impl<'a> TensorSpatialRelaxation<'a> {
    /// Bind an original tensor law already expressed in the first body frame.
    /// The existing material projection owns Mandel/tensor-shear conversion,
    /// volume scaling and all original cross terms. Preparation is bounded by
    /// the caller's existing tensor projection budget and publishes atomically.
    pub fn from_material(
        material: &TensorProny,
        geometry: TensorSpatialGeometry,
        budget: TensorProjectionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if geometry.bodies[0] == geometry.bodies[1]
            || geometry
                .arms_body_m
                .iter()
                .any(|v| ![v.x, v.y, v.z].iter().all(|x| x.is_finite()))
        {
            return Err(TensorSpatialError::Invalid(
                "distinct bodies and finite original attachment arms required",
            ));
        }
        let projection = material.project_strain(
            &geometry.strain_columns,
            geometry.volume_m3,
            budget,
            &mut cancelled,
        )?;
        let rows = |rows: &[Vec<f64>]| -> [[f64; 6]; 6] {
            core::array::from_fn(|i| core::array::from_fn(|j| rows[i][j]))
        };
        let damping = PassiveBushingMatrix::from_gram_rows(
            projection
                .viscosity_gram_rows()
                .map(rows)
                .unwrap_or([[0.0; 6]; 6]),
            geometry.coordinate_scales,
        )?;
        let mut branches = Vec::new();
        branches
            .try_reserve_exact(projection.branches().len() + 1)
            .map_err(|_| TensorSpatialError::Invalid("material memory allocation refused"))?;
        branches.push(SpatialRelaxationBranch::from_equilibrium_gram_rows(rows(
            projection.equilibrium_gram_rows(),
        ))?);
        for branch in projection.branches() {
            poll(&mut cancelled)?;
            branches.push(SpatialRelaxationBranch::from_gram_rows(
                rows(branch.rows()),
                branch.relaxation_time_s(),
            )?);
        }
        poll(&mut cancelled)?;
        let native = RelaxationNetworkLink {
            viscous: ViscousNetworkLink {
                bodies: geometry.bodies,
                arms_body_m: geometry.arms_body_m,
                damping,
            },
            branches,
        };
        Ok(Self {
            geometry,
            projection,
            native,
            source: None,
            material_to_first_body: None,
        })
    }

    /// Bind a calibrated same-card tensor source and an explicit proper
    /// material-to-first-body rotation. The original source selection, physical
    /// query, validity, density provenance and every tensor/time remain reachable.
    /// Source density is provenance only; rigid mass stays with its native owner.
    pub fn from_source(
        source: &'a ResolvedTensorMaterial,
        material_to_first_body: [[f64; 3]; 3],
        geometry: TensorSpatialGeometry,
        budget: TensorProjectionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let law = source.rotated_law(material_to_first_body)?;
        let mut next = Self::from_material(&law, geometry, budget, &mut cancelled)?;
        next.source = Some(source);
        next.material_to_first_body = Some(material_to_first_body);
        poll(&mut cancelled)?;
        Ok(next)
    }

    pub const fn geometry(&self) -> &TensorSpatialGeometry {
        &self.geometry
    }
    pub const fn projection(&self) -> &ProjectedTensorProny {
        &self.projection
    }
    pub const fn native_link(&self) -> &RelaxationNetworkLink {
        &self.native
    }
    pub const fn source(&self) -> Option<&'a ResolvedTensorMaterial> {
        self.source
    }
    pub const fn material_to_first_body(&self) -> Option<[[f64; 3]; 3]> {
        self.material_to_first_body
    }

    /// Admit the explicitly declared application spectrum against every
    /// original source restriction. A timestep never implies a frequency band.
    /// Transient applications remain responsible for a justified source band.
    pub fn admit_frequency_band(&self, band_rad_s: [f64; 2]) -> Result<()> {
        let [low, high] = band_rad_s;
        let original = self.projection.material().band_rad_s();
        if !low.is_finite()
            || !high.is_finite()
            || low < 0.0
            || high < low
            || low < original.0
            || high > original.1
        {
            return Err(TensorSpatialError::Invalid(
                "declared spectrum exceeds original material frequency domain",
            ));
        }
        Ok(())
    }

    /// Original work-conjugate tensor-shear strain and source-domain admission.
    /// Coordinates stay explicit even if the projected stiffness is singular.
    pub fn strain(&self, coordinates: [f64; 6]) -> Result<Voigt> {
        Ok(self.projection.strain(&coordinates)?)
    }

    /// Initialize equilibrium at a specified total strain, with every Maxwell
    /// branch relaxed at that same state. This is an explicit preparation of
    /// history, not a reset to apply during a continuation call. A nonzero
    /// relaxed strain requires original DC support, as in the material history
    /// owner's relaxed initializer; zero strain remains valid for an AC-only law.
    pub fn relaxed_memory(&self, coordinates: [f64; 6]) -> Result<Vec<[f64; 6]>> {
        let strain = self.strain(coordinates)?;
        if strain.iter().any(|x| *x != 0.0) {
            self.projection.material().modulus_mandel(0.0)?;
        }
        let mut memory = vec![[0.0; 6]; self.native.branches.len()];
        memory[0] = self.native.branches[0].elastic_memory(coordinates)?;
        Ok(memory)
    }

    /// Initialize a distinct physically specified elastic displacement for
    /// every original Maxwell spring. Channel zero always uses total strain.
    /// Subsequent native endpoint memories must be retained in full; never
    /// reconstruct Maxwell state from the endpoint pose or total strain.
    pub fn memory_from_elastic_coordinates(
        &self,
        coordinates: [f64; 6],
        maxwell_elastic_coordinates: &[[f64; 6]],
    ) -> Result<Vec<[f64; 6]>> {
        self.strain(coordinates)?;
        if maxwell_elastic_coordinates.len() != self.projection.branches().len() {
            return Err(TensorSpatialError::Invalid(
                "one elastic displacement per original Maxwell branch required",
            ));
        }
        let mut memory = Vec::with_capacity(self.native.branches.len());
        memory.push(self.native.branches[0].elastic_memory(coordinates)?);
        for (branch, &elastic) in self.native.branches[1..]
            .iter()
            .zip(maxwell_elastic_coordinates)
        {
            self.strain(elastic)?;
            memory.push(branch.elastic_memory(elastic)?);
        }
        Ok(memory)
    }

    /// Continue all six source coordinates using this connector's accepted
    /// exact native source integral. The caller must select this binding's
    /// original connector index and retain the receipt's complete memory.
    ///
    /// Both supplied endpoints satisfy the source strain domain. This is not a
    /// certificate for every interior time; resolve large oscillations by
    /// subdivision. No factor inverse or pose-derived strain is used.
    pub fn coordinates_after_source(
        &self,
        before: [f64; 6],
        connector: usize,
        receipt: &RelaxationNetworkReceipt,
    ) -> Result<[f64; 6]> {
        self.strain(before)?;
        let twist =
            receipt
                .average_relative_twist
                .get(connector)
                .ok_or(TensorSpatialError::Invalid(
                    "original connector index is absent from source receipt",
                ))?;
        let memory = receipt
            .memory_after
            .get(connector)
            .ok_or(TensorSpatialError::Invalid(
                "source receipt has no complete material memory",
            ))?;
        let dt = receipt.source_duration_s;
        if !dt.is_finite()
            || dt < 0.0
            || twist.iter().any(|x| !x.is_finite())
            || memory.len() != self.native.branches.len()
            || memory.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(TensorSpatialError::Invalid(
                "source receipt does not retain the complete original material coordinates",
            ));
        }
        let next = core::array::from_fn(|i| {
            (Dd::from_f64(before[i]) + Dd::from_f64(dt) * Dd::from_f64(twist[i])).to_f64()
        });
        self.strain(next)?;
        Ok(next)
    }

    /// Continue from both original fixed-pose half-source integrals. The
    /// intervening native drift carries corotated coordinates unchanged.
    /// Every half endpoint is admitted before returning the complete candidate;
    /// a refusal leaves the caller's body, coordinate and memory banks untouched.
    pub fn coordinates_after_motion(
        &self,
        before: [f64; 6],
        connector: usize,
        receipt: &RelaxationNetworkMotionStep,
    ) -> Result<[f64; 6]> {
        let middle = self.coordinates_after_source(before, connector, &receipt.half_sources[0])?;
        self.coordinates_after_source(middle, connector, &receipt.half_sources[1])
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(TensorSpatialError::Cancelled)
    } else {
        Ok(())
    }
}
