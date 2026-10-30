//! Finite conductive solids attached to the resolved gas-temperature field.
//!
//! For `theta=T'/T0`, continuity of temperature and normal heat flux adds
//! `T0 * Y_w(omega) * integral(theta_test * theta dS)` to the original heat
//! equation. Its cycle-averaged work is
//! `Re(Y_w)/(2*T0) * integral(|T'|^2 dS)`. The gas has no additional thermal
//! boundary-layer term: its complete volume heat equation already owns that
//! diffusion. Each solid uses the original passive normal-profile admittance.
//! This is a local normal solid model, not a lateral or multilayer solid mesh.

use super::{ThermoviscousError, ThermoviscousFe, canonical, finite, poll};
use crate::acoustics::helmholtz::Medium;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::ConductiveThermalWall;
use crate::acoustics::helmholtz::tri_quadrature;
use crate::highorder::simplex::LocalFn;
use fs_math::c64::C64;
use std::collections::BTreeSet;

/// One homogeneous solid patch, perfectly contacting stationary no-slip gas
/// faces. Different patches may have different material, thickness and backing.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveWallRegion {
    /// Unique exterior faces from the original no-slip set, disjoint from
    /// isothermal faces and every other conductive region.
    pub faces: Vec<[u32; 3]>,
    /// Original solid normal conduction owner; coefficients are in SI.
    pub wall: ConductiveThermalWall,
    /// Solid normal-profile admission. Only the solid profile is tested; the
    /// gas may have overlapping thermal and viscous layers. Half-gap is not a
    /// solid thickness and does not restrict the resolved gas here.
    pub geometry: BoundaryLayerGeometry,
    /// Mean solid temperature [K], exactly the gas equilibrium temperature.
    pub temperature_k: f64,
}

#[derive(Clone)]
struct FaceTrace {
    dofs: [usize; 6],
    area: f64,
}

pub(super) struct WallEntry {
    pub row: usize,
    pub col: usize,
    pub region: usize,
    pub mass: f64,
}

#[derive(Default)]
pub(super) struct PreparedWallTraces {
    regions: Vec<ConductiveWallRegion>,
    faces: Vec<Vec<FaceTrace>>,
    pub entries: Vec<WallEntry>,
    row_ptr: Vec<usize>,
}

impl ThermoviscousFe<'_> {
    /// Attach or replace finite conductive wall regions. This retains every
    /// original velocity, pressure and temperature coordinate and every port.
    /// Exact isothermal constraints cannot be replaced after construction.
    ///
    /// Face membership, overlap and common mean temperature are admitted here.
    /// Solid transport/depth and normal-profile geometry are admitted at each
    /// frequency by the original `ConductiveThermalWall::normal_admittance`.
    /// Sparse trace work is bounded by the already admitted mesh: at most 36
    /// entries per face before deterministic accumulation, with no additional
    /// dense matrix. Empty `regions` restores the original insulated remainder.
    pub fn with_conductive_walls(
        mut self,
        regions: &[ConductiveWallRegion],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        if regions.len() > self.wall_faces.len() {
            return Err(ThermoviscousError::InvalidInput("conductive region count"));
        }
        let mut seen = BTreeSet::new();
        let mut traces = PreparedWallTraces::default();
        for (region_index, region) in regions.iter().enumerate() {
            poll(cancelled)?;
            if region.faces.is_empty()
                || region.faces.len() > self.wall_faces.len()
                || region.temperature_k.to_bits() != self.medium.temperature_k.to_bits()
            {
                return Err(ThermoviscousError::InvalidInput(
                    "conductive faces or common gas/solid temperature",
                ));
            }
            let mut faces = Vec::with_capacity(region.faces.len());
            let entry_start = traces.entries.len();
            for &face in &region.faces {
                poll(cancelled)?;
                let key = canonical(face);
                if !self.wall_faces.contains(&key)
                    || self.isothermal_faces.contains(&key)
                    || !seen.insert(key)
                {
                    return Err(ThermoviscousError::InvalidInput(
                        "conductive/isothermal/no-slip face partition",
                    ));
                }
                let (_, _, area) = self.fe.boundary_geometry(key)?;
                let mut dofs = [0; 6];
                for a in 0..3 {
                    dofs[a] = key[a] as usize;
                }
                for (offset, (a, b)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                    let edge = self
                        .space
                        .complex
                        .edges
                        .binary_search(&[key[a], key[b]])
                        .expect("original canonical boundary edge");
                    dofs[3 + offset] = self.space.edge_off + edge;
                }
                // Reuse the original face quadrature and hierarchical P2
                // assembly; constrain only the exact isothermal coordinates.
                let mass = self.fe.boundary_mass(&[key])?;
                for &original_row in &dofs {
                    let Some(row) = self.temperature_map[original_row] else {
                        continue;
                    };
                    let (cols, values) = mass.row(original_row);
                    for (&original_col, &value) in cols.iter().zip(values) {
                        if let Some(col) = self.temperature_map[original_col] {
                            if !value.is_finite() {
                                return Err(ThermoviscousError::NonFinite);
                            }
                            if value != 0.0 {
                                traces.entries.push(WallEntry {
                                    row: self.temperature_start + row,
                                    col: self.temperature_start + col,
                                    region: region_index,
                                    mass: value,
                                });
                            }
                        }
                    }
                }
                faces.push(FaceTrace { dofs, area });
            }
            if traces.entries.len() == entry_start {
                return Err(ThermoviscousError::InvalidInput(
                    "unresolved conductive temperature trace",
                ));
            }
            traces.faces.push(faces);
            traces.regions.push(region.clone());
        }
        // Keep region identity while accumulating shared nodes/edges. Stable
        // sorting preserves the declared face summation order for replay.
        traces.entries.sort_by_key(|e| (e.row, e.col, e.region));
        let mut accumulated: Vec<WallEntry> = Vec::with_capacity(traces.entries.len());
        for entry in traces.entries {
            if let Some(last) = accumulated.last_mut() {
                if (last.row, last.col, last.region) == (entry.row, entry.col, entry.region) {
                    last.mass += entry.mass;
                    if !last.mass.is_finite() {
                        return Err(ThermoviscousError::NonFinite);
                    }
                    continue;
                }
            }
            accumulated.push(entry);
        }
        traces.entries = accumulated;
        traces.row_ptr = vec![0; self.n + 1];
        for entry in &traces.entries {
            traces.row_ptr[entry.row + 1] += 1;
        }
        for row in 0..self.n {
            traces.row_ptr[row + 1] += traces.row_ptr[row];
        }
        poll(cancelled)?;
        self.conductive = traces;
        Ok(self)
    }

    /// Original finite solid descriptions in declared region order.
    #[must_use]
    pub fn conductive_walls(&self) -> &[ConductiveWallRegion] {
        &self.conductive.regions
    }
}

impl PreparedWallTraces {
    pub fn admittances(
        &self,
        medium: Medium,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousError> {
        let mut result = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            poll(cancelled)?;
            let [re, im] = region
                .wall
                .normal_admittance(medium, region.geometry, frequency_hz)?;
            let y = C64::new(re, im);
            if !finite(y) || re < 0.0 {
                return Err(ThermoviscousError::NonFinite);
            }
            result.push(y);
        }
        Ok(result)
    }

    pub fn row(&self, row: usize) -> &[WallEntry] {
        if self.entries.is_empty() {
            &[]
        } else {
            &self.entries[self.row_ptr[row]..self.row_ptr[row + 1]]
        }
    }

    // Independent pointwise integration, rather than reusing the assembled
    // surface quadratic form. Vertex and edge shapes are the original P2
    // LocalFn descriptors, including their hierarchical normalization.
    pub fn dissipated_w(
        &self,
        temperature: &[C64],
        mean_temperature: f64,
        admittances: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ThermoviscousError> {
        let quad = tri_quadrature(3);
        let shapes = [
            LocalFn::Vertex(0),
            LocalFn::Vertex(1),
            LocalFn::Vertex(2),
            LocalFn::Edge(0, 1, 0),
            LocalFn::Edge(0, 2, 0),
            LocalFn::Edge(1, 2, 0),
        ];
        let mut power = 0.0;
        for (faces, admittance) in self.faces.iter().zip(admittances) {
            for face in faces {
                poll(cancelled)?;
                for &(bary, weight) in &quad {
                    let lam = [bary[0], bary[1], bary[2], 0.0];
                    let temp = face
                        .dofs
                        .iter()
                        .zip(&shapes)
                        .fold(C64::ZERO, |sum, (&dof, shape)| {
                            sum + temperature[dof].scale(shape.eval(lam, 2))
                        });
                    power += admittance.re / mean_temperature * face.area * weight * temp.norm_sq();
                }
            }
        }
        if !power.is_finite() {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(power)
    }
}
