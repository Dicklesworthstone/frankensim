//! Source-bearing Ohmic current on the original tetrahedral Whitney complex.
//!
//! K=d0^T M_sigma d0, E=-grad(V), J=sigma E. Geometry, P1 face moments,
//! component gauges and checked CG come from the existing spatial owner;
//! conductivity, current and irreversible E.J power keep electrical units.
//! Cellwise fields are not smoothed across material interfaces. This is a
//! stationary, fixed-geometry, symmetric positive DC law, without induction,
//! Hall response, dielectric storage or field-dependent sigma. Matching material
//! contacts use the same operator with exact triangular voltage-jump traces.

mod contact;
pub use contact::{SpatialConductionContact, SpatialConductionContactResponse};
mod terminals;
pub use terminals::SpatialResistance;
mod matrix;
pub use matrix::{SpatialConductanceControls, SpatialConductanceMatrix};
pub mod motion;

use super::{
    ScalarBoundary, SpatialGeometry, SpatialMagneticError, dot, prepare_geometry,
    prepare_scalar_boundary_traces, solve_scalar_system,
};
use crate::conduction::ConductionControls;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::conductor::{
    ELECTRICAL_CONDUCTIVITY_DIMS, ResolvedConductorState,
    spatial::{SpatialConductivity, SpatialConductorMaterial},
};
use fs_rep_mesh::TetComplex;
use fs_solver::SolveReport;
use fs_sparse::Csr;
use std::{collections::BTreeSet, sync::Arc};

/// P1 voltage trace on one original EXTERIOR face. Values follow the mesh's
/// sorted face vertices; adjacent traces must agree at shared vertices.
#[derive(Clone, Copy, Debug)]
pub struct ConductionVoltageFace {
    pub face: usize,
    pub voltage_v: [f64; 3],
}

/// Uniform current density INTO the body on an original exterior face.
/// This is a distributed Neumann trace, not an equipotential electrode.
#[derive(Clone, Copy, Debug)]
pub struct ConductionCurrentFace {
    pub face: usize,
    pub inward_current_density_a_m2: f64,
}

/// Unlisted faces are insulating. A current-only connected component needs
/// zero net imposed current and receives one deterministic voltage gauge.
#[derive(Clone, Debug, Default)]
pub struct SpatialConductionBoundary {
    pub voltage: Vec<ConductionVoltageFace>,
    pub current: Vec<ConductionCurrentFace>,
}

#[derive(Debug)]
pub enum SpatialConductionError {
    /// Original shared tetrahedral geometry, gauge or numerical refusal.
    Native(SpatialMagneticError),
    Assembly(WeightedError),
    Invalid(&'static str),
    Cancelled,
    Balance {
        current_residual_a: f64,
        maximum_free_current_residual_a: f64,
        power_residual_w: f64,
        hodge_power_residual_w: f64,
    },
}
impl From<SpatialMagneticError> for SpatialConductionError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SpatialConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral Ohmic conduction: {self:?}")
    }
}
impl std::error::Error for SpatialConductionError {}
type Result<T> = core::result::Result<T, SpatialConductionError>;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialConductionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(SpatialConductionError::Invalid(
            "unrepresentable electrical arithmetic",
        ))
    }
}

/// One complete immutable source bundle per original tetrahedron. Scalar and
/// oriented principal-source laws share this same native Hodge and solve.
/// Material refresh reuses the admitted topology/geometry through an Arc.
#[derive(Clone)]
pub struct SpatialConduction<C = ResolvedConductorState> {
    geometry: Arc<SpatialGeometry>,
    conductors: Vec<C>,
    conductivity: Vec<SpatialConductivity>,
    stiffness: Csr,
    contacts: Vec<SpatialConductionContact>,
    contact_liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
    contact_faces: BTreeSet<usize>,
    contact_cells: Vec<[usize; 2]>,
    contact_conductance_s: Vec<f64>,
    /// Electrical connectivity includes conducting contacts between disjoint
    /// volume complexes; the geometry itself retains its original components.
    components: Vec<usize>,
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        conductors: Vec<C>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if tetrahedra.len() != conductors.len() {
            return Err(SpatialConductionError::Invalid(
                "one conductor per original tetrahedron required",
            ));
        }
        let geometry = Arc::new(prepare_geometry(
            positions_m,
            tetrahedra,
            limits,
            cancelled,
        )?);
        Self::assemble(geometry, conductors, limits, cancelled)
    }

    /// Prepare a new endpoint operator on the SAME original volume geometry.
    /// The caller owns physical source refresh and its temporal/energy account.
    pub fn with_conductors(
        &self,
        conductors: Vec<C>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::assemble(Arc::clone(&self.geometry), conductors, limits, cancelled)?.bind_contacts(
            self.contacts.clone(),
            self.contact_liquid_mass_fraction_axes.clone(),
            limits,
            cancelled,
        )
    }

    fn assemble(
        geometry: Arc<SpatialGeometry>,
        conductors: Vec<C>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if conductors.len() != geometry.mesh.tets.len()
            || geometry.mesh.vertex_count > limits.max_dofs
        {
            return Err(SpatialConductionError::Invalid(
                "complete bounded conductor geometry required",
            ));
        }
        let mut conductivity = Vec::with_capacity(conductors.len());
        let mut weights = Vec::with_capacity(conductors.len());
        for source in &conductors {
            check(cancelled)?;
            if !source.has_admitted_ohmic_sources() {
                return Err(SpatialConductionError::Invalid(
                    "original electrical sources need fixed external or explicitly owned phase conditions",
                ));
            }
            let coefficient = source.conductivity();
            weights.push(match coefficient {
                SpatialConductivity::Scalar(s) => CellWeight::Scalar(s),
                SpatialConductivity::Tensor(s) => CellWeight::Tensor(s),
            });
            conductivity.push(coefficient);
        }
        let mass = weighted_mass_matrix(
            &geometry.mesh,
            &geometry.geometry,
            1,
            &weights,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|e| {
            if e == WeightedError::Cancelled {
                SpatialConductionError::Cancelled
            } else {
                SpatialConductionError::Assembly(e)
            }
        })?;
        check(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&geometry.mesh.d0()), mass.matrix());
        for row in 0..stiffness.nrows() {
            check(cancelled)?;
            for &value in stiffness.row(row).1 {
                finite(value)?;
            }
        }
        let components = geometry.components.clone();
        Ok(Self {
            geometry,
            conductors,
            conductivity,
            stiffness,
            contacts: Vec::new(),
            contact_liquid_mass_fraction_axes: Vec::new(),
            contact_faces: BTreeSet::new(),
            contact_cells: Vec::new(),
            contact_conductance_s: Vec::new(),
            components,
        })
    }

    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn conductors(&self) -> &[C] {
        &self.conductors
    }
    pub fn conductivity_coefficients(&self) -> &[SpatialConductivity] {
        &self.conductivity
    }
    pub fn cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.geometry.geometry.vol_signed.get(cell).map(|v| v.abs())
    }
    pub fn exterior_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.geometry
            .faces
            .iter()
            .enumerate()
            .filter_map(|(id, face)| {
                (face.cells.len() == 1 && !self.contact_faces.contains(&id)).then_some(id)
            })
    }
    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.geometry
            .faces
            .get(face)
            .filter(|f| f.cells.len() == 1 && !self.contact_faces.contains(&face))
            .map(|f| f.outward_normal)
    }
    pub fn face_area_m2(&self, face: usize) -> Option<f64> {
        self.geometry.faces.get(face).map(|f| f.area_m2)
    }

    /// Solve -div(sigma grad(V))=0 with exact P1 face-current moments. Current
    /// continuity, recovered E.J versus Hodge power, and independently supplied
    /// terminal work are accepted separately. A component gauge supplies no
    /// current. These residual gates are not a spatial error certificate.
    pub fn solve(
        &self,
        boundary: &SpatialConductionBoundary,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        let traces = self.prepare_boundary(boundary, controls, cancelled)?;
        // The load is positive INTO the conductor. The shared trace helper
        // only integrates signed data; unlike magnetic B.n no sign reversal
        // belongs between the electrical load and K V.
        let (relative, report) = solve_scalar_system(
            &self.stiffness,
            &traces.flux_load_wb,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_field(
            relative,
            traces,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )
    }

    fn prepare_boundary(
        &self,
        boundary: &SpatialConductionBoundary,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ScalarBoundary> {
        check(cancelled)?;
        if !controls.current_tolerance_a.is_finite()
            || controls.current_tolerance_a <= 0.0
            || !controls.power_tolerance_w.is_finite()
            || controls.power_tolerance_w <= 0.0
            || !controls.linear.tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.linear.tolerance)
            || controls.linear.tolerance == 0.0
            || boundary.voltage.len() > self.geometry.faces.len()
            || boundary.current.len() > self.geometry.faces.len()
        {
            return Err(SpatialConductionError::Invalid(
                "positive finite electrical controls and bounded face traces required",
            ));
        }
        for face in boundary
            .voltage
            .iter()
            .map(|b| b.face)
            .chain(boundary.current.iter().map(|b| b.face))
        {
            check(cancelled)?;
            if self.contact_faces.contains(&face) {
                return Err(SpatialConductionError::Invalid(
                    "an internal material contact cannot also be an external electrical trace",
                ));
            }
        }
        let g = &self.geometry;
        Ok(prepare_scalar_boundary_traces(
            &g.mesh,
            &g.faces,
            &self.components,
            boundary.voltage.iter().map(|t| (t.face, t.voltage_v)),
            boundary
                .current
                .iter()
                .map(|t| (t.face, t.inward_current_density_a_m2)),
            controls.current_tolerance_a,
            cancelled,
        )?)
    }

    fn recover_field(
        &self,
        relative: Vec<f64>,
        traces: ScalarBoundary,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        self.recover_field_with_motion(
            relative,
            traces,
            None,
            current_tolerance_a,
            power_tolerance_w,
            report,
            cancelled,
        )
        .map(|(field, _)| field)
    }

    #[allow(clippy::too_many_arguments)]
    fn recover_field_with_motion(
        &self,
        relative: Vec<f64>,
        traces: ScalarBoundary,
        motion: Option<&motion::PreparedMotion<'_>>,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(SpatialConductionField, Option<motion::MotionPower>)> {
        let ScalarBoundary {
            physical_prescribed,
            flux_load_wb: loads,
            gauges,
            references,
            ..
        } = traces;
        let g = &self.geometry;
        let n = g.mesh.vertex_count;
        let mut applied = vec![0.0; n];
        let mut nodal_boundary_current_a = vec![None; n];
        let mut current_residual_a = 0.0;
        let mut maximum_free_current_residual_a = 0.0_f64;
        let mut supplied_power_w = 0.0;
        let mut hodge_power_w = 0.0;
        for v in 0..n {
            check(cancelled)?;
            let (columns, values) = self.stiffness.row(v);
            applied[v] = finite(
                columns
                    .iter()
                    .zip(values)
                    .map(|(&j, &value)| value * relative[j])
                    .sum(),
            )?;
            if let Some(motion) = motion {
                applied[v] = finite(applied[v] - motion.load_a[v])?;
            }
            hodge_power_w = finite(relative[v].mul_add(applied[v], hodge_power_w))?;
            let current = if physical_prescribed[v] {
                Some(applied[v])
            } else {
                maximum_free_current_residual_a =
                    maximum_free_current_residual_a.max(finite(applied[v] - loads[v])?.abs());
                (loads[v] != 0.0).then_some(loads[v])
            };
            if let Some(current) = current {
                current_residual_a = finite(current_residual_a + current)?;
                supplied_power_w = finite(relative[v].mul_add(current, supplied_power_w))?;
            }
            nodal_boundary_current_a[v] = current;
        }
        let mut electric_field_v_m = Vec::with_capacity(self.conductors.len());
        let mut current_density_a_m2 = Vec::with_capacity(self.conductors.len());
        let mut cell_power_w = Vec::with_capacity(self.conductors.len());
        let mut transfer = motion.map(|_| motion::MotionPower::new(self.conductors.len()));
        let mut power_w = 0.0;
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let mut e = [0.0; 3];
            for local in 1..4 {
                let difference = finite(relative[tet[local] as usize] - relative[tet[0] as usize])?;
                for (j, out) in e.iter_mut().enumerate() {
                    *out = finite((-difference).mul_add(g.geometry.grads[cell][local][j], *out))?;
                }
            }
            let mut conductor_e = e;
            if let Some(motion) = motion
                && motion.field_v_m[cell] != [0.0; 3]
            {
                for (axis, value) in conductor_e.iter_mut().enumerate() {
                    *value = finite(*value + motion.field_v_m[cell][axis])?;
                }
            }
            let j = self.conductivity[cell].current_density_a_m2(conductor_e);
            for value in j {
                finite(value)?;
            }
            let volume = g.geometry.vol_signed[cell].abs();
            let power = finite(dot(conductor_e, j) * volume)?;
            if power < 0.0 || (conductor_e.iter().any(|&x| x != 0.0) && power == 0.0) {
                return Err(SpatialConductionError::Invalid(
                    "positive conductivity lost nonzero irreversible heat",
                ));
            }
            power_w = finite(power_w + power)?;
            electric_field_v_m.push(e);
            current_density_a_m2.push(j);
            cell_power_w.push(power);
            if let (Some(motion), Some(transfer)) = (motion, transfer.as_mut()) {
                transfer.record(motion.prescribed, cell, e, conductor_e, j, power, volume)?;
            }
        }
        let contacts = self.contact_response(&relative, cancelled)?;
        for contact in &contacts {
            check(cancelled)?;
            power_w = finite(power_w + contact.power_w)?;
        }
        let (accounted_power_w, electrical_power_w, work_pairing_residual_w) = match &transfer {
            Some(transfer) => (
                finite(power_w + transfer.mechanical_power_w)?,
                transfer.electrical_power_w,
                transfer.work_pairing_residual_w,
            ),
            None => (power_w, power_w, 0.0),
        };
        let power_residual_w = finite(accounted_power_w - supplied_power_w)?;
        let hodge_power_residual_w = finite(electrical_power_w - hodge_power_w)?;
        if current_residual_a.abs() > current_tolerance_a
            || maximum_free_current_residual_a > current_tolerance_a
            || power_residual_w.abs() > power_tolerance_w
            || hodge_power_residual_w.abs() > power_tolerance_w
            || work_pairing_residual_w > power_tolerance_w
        {
            return Err(SpatialConductionError::Balance {
                current_residual_a,
                maximum_free_current_residual_a,
                power_residual_w: if work_pairing_residual_w > power_tolerance_w {
                    work_pairing_residual_w
                } else {
                    power_residual_w
                },
                hodge_power_residual_w,
            });
        }
        let voltage_v = relative
            .iter()
            .enumerate()
            .map(|(v, &value)| finite(value + references[self.components[v]]))
            .collect::<Result<Vec<_>>>()?;
        check(cancelled)?;
        let field = SpatialConductionField {
            voltage_v,
            electric_field_v_m,
            current_density_a_m2,
            cell_power_w,
            contacts,
            power_w,
            supplied_power_w,
            nodal_boundary_current_a,
            nodal_applied_current_a: loads,
            current_residual_a,
            maximum_free_current_residual_a,
            power_residual_w,
            hodge_power_residual_w,
            gauges,
            linear: report,
        };
        Ok((field, transfer))
    }
}

/// Actual native P1 solution. Boundary current is positive INTO the body;
/// prescribed-current corners are included in each voltage reaction once.
/// Current-only components retain their explicitly listed algebraic gauges.
#[derive(Clone, Debug)]
pub struct SpatialConductionField {
    pub voltage_v: Vec<f64>,
    pub electric_field_v_m: Vec<[f64; 3]>,
    pub current_density_a_m2: Vec<[f64; 3]>,
    /// Bulk Joule heat only; surface contact heat is retained separately below.
    pub cell_power_w: Vec<f64>,
    /// Exact contact trace work in the same canonical order as `contacts()`.
    pub contacts: Vec<SpatialConductionContactResponse>,
    /// Total irreversible power, including bulk and contact contributions once.
    pub power_w: f64,
    pub supplied_power_w: f64,
    pub nodal_boundary_current_a: Vec<Option<f64>>,
    pub nodal_applied_current_a: Vec<f64>,
    pub current_residual_a: f64,
    pub maximum_free_current_residual_a: f64,
    pub power_residual_w: f64,
    pub hodge_power_residual_w: f64,
    pub gauges: Vec<usize>,
    pub linear: Option<SolveReport>,
}
