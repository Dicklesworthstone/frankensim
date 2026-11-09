//! Finite-volume inlet data from the existing constrained turbulent frame.
//!
//! The two strictly increasing coordinate partitions cover [-1,1] on the
//! original patch, so their tensor product covers that patch without gaps or
//! overlap. Its constant area Jacobian makes the existing analytic Fourier
//! box filter the exact face average, including for skew parallelograms.
//! No face-center approximation or mesh-dependent flux repair is applied.
//!
//! For the existing anelastic field, rho*u'=rho_ref*v. Each face therefore
//! receives q'=rho_ref*<v> and the Favre velocity u'_F=q'/<rho>. Its positive
//! mean density is integrated analytically using the original bounded
//! exponential kernel. The sum of area*q' dot normal recovers the frame's
//! zero whole-patch turbulent mass flow, to representable roundoff.
//!
//! Callers supply their mean inflow and thermal state. In particular, a
//! spatially uniform mean U is added directly to u'_F, giving prescribed
//! normal mass flow sum(area*<rho>)*U dot normal. These mass-consistent face
//! states do not carry the unresolved momentum or kinetic-energy flux; use
//! the original moment observations if those are required by a solver.

use super::{MannInflowFrame, MannInflowPatch, storage};
use crate::fetch_mass_consistent::ExponentialDensityProfile;
use crate::mann_field::wall::filter::scaled_edge_gain;
use crate::mann_field::{MannBoxFilter, dot};
use crate::{MAX_MODES, Refusal, refuse};

/// Authored physical density used to convert the constrained perturbation into
/// mass-consistent reservoir data. Anelastic height is the Mann frame's +z.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MannInflowDensity {
    /// Spatially uniform positive density [kg/m³]; velocity is unchanged.
    ConstantKgM3(f64),
    /// The existing finite-domain exponential-density anelastic construction.
    Exponential(ExponentialDensityProfile),
}

/// One fixed physical face. Its normal is the parent patch's ordered normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannInflowFace {
    /// Physical center in the original Mann frame [m].
    pub center_m: [f64; 3],
    /// Half edges in the same order as the two coordinate partitions [m].
    pub half_edges_m: [[f64; 3]; 2],
    /// Physical face area [m²].
    pub area_m2: f64,
}

/// Exact area-integrated mass data, normalized by the physical face area.
/// All velocities and fluxes are turbulent perturbations; add the caller's
/// mean separately. For uniform density, Favre and uniform velocity agree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannInflowFaceSample {
    /// Uniform area average of physical density [kg/m³].
    pub mean_density_kg_m3: f64,
    /// Density-weighted velocity <rho*u'>/<rho> [m/s].
    pub favre_velocity_m_s: [f64; 3],
    /// Uniform area average <rho*u'> [kg/(m² s)].
    pub mean_mass_flux_kg_m2_s: [f64; 3],
    /// Signed outward flux under the patch's ordered normal [kg/(m² s)].
    pub normal_mass_flux_kg_m2_s: f64,
}

const EMPTY: MannInflowFaceSample = MannInflowFaceSample {
    mean_density_kg_m3: 0.0,
    favre_velocity_m_s: [0.0; 3],
    mean_mass_flux_kg_m2_s: [0.0; 3],
    normal_mass_flux_kg_m2_s: 0.0,
};

#[derive(Clone, Debug, PartialEq)]
struct PreparedFace {
    geometry: MannInflowFace,
    filter: MannBoxFilter,
    mean_density_kg_m3: f64,
    velocity_scale: f64,
}

/// Reusable conservative inlet plane. The first coordinate varies fastest:
/// flat_index([i,j]) = i + dimensions()[0]*j. Choosing its first half edge
/// along the lowest tangential Cartesian axis matches Cartesian cell order.
/// Geometry and density are fixed at preparation. Refresh reads one accepted
/// MannInflowFrame and owns neither a stochastic process nor a second clock.
#[derive(Clone, Debug)]
pub struct MannInflowPlane {
    patch: MannInflowPatch,
    dimensions: [usize; 2],
    reference_density_kg_m3: f64,
    faces: Vec<PreparedFace>,
    samples: Vec<MannInflowFaceSample>,
    trial: Vec<MannInflowFaceSample>,
    time_s: f64,
    turbulent_mass_flow_kg_s: f64,
}

impl PartialEq for MannInflowPlane {
    fn eq(&self, other: &Self) -> bool {
        self.patch == other.patch
            && self.dimensions == other.dimensions
            && self.reference_density_kg_m3 == other.reference_density_kg_m3
            && self.faces == other.faces
            && self.samples == other.samples
            && self.time_s == other.time_s
            && self.turbulent_mass_flow_kg_s == other.turbulent_mass_flow_kg_s
    }
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use a matching constrained Mann frame, complete increasing coordinate partitions, positive representable density, and sufficient inlet storage",
    )
}

fn numerical() -> Refusal {
    error(
        "mann-inlet-plane-range",
        "face geometry, density, velocity or mass flow is not representable",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-inlet-plane-cancelled",
            "finite-volume inlet sampling cancelled",
        ))
    } else {
        Ok(())
    }
}

fn face_count(dimensions: [usize; 2]) -> Result<usize, Refusal> {
    if dimensions.contains(&0) {
        return Err(numerical());
    }
    dimensions[0]
        .checked_mul(dimensions[1])
        .ok_or_else(numerical)
}

fn partition(edges: &[f64]) -> Result<usize, Refusal> {
    if edges.len() < 2
        || edges.first() != Some(&-1.0)
        || edges.last() != Some(&1.0)
        || !edges.iter().all(|x| x.is_finite())
        || !edges.windows(2).all(|pair| pair[0] < pair[1])
    {
        return Err(error(
            "mann-inlet-plane-partition",
            "each coordinate partition must increase strictly from -1 to 1",
        ));
    }
    Ok(edges.len() - 1)
}

impl MannInflowDensity {
    fn reference(self) -> Result<f64, Refusal> {
        let value = match self {
            Self::ConstantKgM3(value) => value,
            Self::Exponential(profile) => profile.reference_density_kg_m3(),
        };
        if !value.is_finite() || value <= 0.0 {
            return Err(numerical());
        }
        Ok(value)
    }

    fn face_density(self, face: MannInflowFace) -> Result<(f64, f64), Refusal> {
        match self {
            Self::ConstantKgM3(value) => Ok((value, 1.0)),
            Self::Exponential(profile) => {
                let radius = face.half_edges_m.iter().map(|e| e[2].abs()).sum::<f64>();
                let low = profile.density_kg_m3(face.center_m[2] - radius)?;
                profile.density_kg_m3(face.center_m[2] + radius)?;
                let kappa = profile.inverse_scale_height_per_m();
                let mut mean = profile.density_kg_m3(face.center_m[2])?;
                if kappa != 0.0 && radius != 0.0 {
                    // exp(-kappa*radius)*<exp(-kappa*delta_z)> is a product
                    // of bounded factors. The low endpoint supplies the
                    // opposite normalization, avoiding a growing exponential.
                    let mut factor = 1.0;
                    for edge in face.half_edges_m {
                        let exponent = kappa * edge[2];
                        if !exponent.is_finite() {
                            return Err(numerical());
                        }
                        factor *= scaled_edge_gain(exponent, 0.0)[0];
                    }
                    mean = low * factor;
                }
                // Preserve the original velocity arithmetic when kappa=0.
                let scale = if kappa == 0.0 {
                    1.0
                } else {
                    profile.reference_density_kg_m3() / mean
                };
                if ![mean, scale].into_iter().all(|x| x.is_finite() && x > 0.0) {
                    return Err(numerical());
                }
                Ok((mean, scale))
            }
        }
    }
}

impl MannInflowPlane {
    /// Required heap capacity for face metadata, all modal filter gains and
    /// both accepted/trial sample arrays. Fixed owner metadata and stack are
    /// excluded. The same budget is checked against actual allocated capacity.
    pub fn required_workspace_bytes(
        mode_count: usize,
        dimensions: [usize; 2],
    ) -> Result<usize, Refusal> {
        if !(1..=MAX_MODES).contains(&mode_count) {
            return Err(numerical());
        }
        face_count(dimensions)?
            .checked_mul(
                core::mem::size_of::<PreparedFace>()
                    + 2 * core::mem::size_of::<MannInflowFaceSample>()
                    + mode_count * core::mem::size_of::<f64>(),
            )
            .ok_or_else(numerical)
    }

    /// Prepare exact face averages and sample the supplied frame immediately.
    /// Coordinate edges are dimensionless patch coordinates, including both
    /// -1 and +1 endpoints. The complete physical patch must lie within an
    /// exponential profile's declared height interval. The supplied byte
    /// budget bounds retained storage and thus faces*modes work per refresh.
    /// Preparation polls cancellation before allocation and between modes/faces.
    pub fn prepare(
        frame: &MannInflowFrame,
        coordinate_edges: [&[f64]; 2],
        density: MannInflowDensity,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(&mut cancelled)?;
        let dimensions = [
            partition(coordinate_edges[0])?,
            partition(coordinate_edges[1])?,
        ];
        let count = face_count(dimensions)?;
        let n = frame.model.mode_count();
        if Self::required_workspace_bytes(n, dimensions)? > maximum_workspace_bytes {
            return Err(error(
                "mann-inlet-plane-budget",
                "inlet storage exceeds its byte budget",
            ));
        }
        let reference_density_kg_m3 = density.reference()?;
        let patch = frame.patch();
        density.face_density(MannInflowFace {
            center_m: patch.center_m(),
            half_edges_m: patch.half_edges_m(),
            area_m2: patch.area_m2(),
        })?;
        let mut result = Self {
            patch,
            dimensions,
            reference_density_kg_m3,
            faces: storage(count)?,
            samples: storage(count)?,
            trial: storage(count)?,
            time_s: frame.time_s(),
            turbulent_mass_flow_kg_s: 0.0,
        };
        let filter_bytes = n * core::mem::size_of::<f64>();
        if result
            .workspace_bytes()
            .checked_add(count * filter_bytes)
            .is_none_or(|bytes| bytes > maximum_workspace_bytes)
        {
            return Err(error(
                "mann-inlet-plane-budget",
                "allocated inlet capacity exceeds its byte budget",
            ));
        }
        for j in 0..dimensions[1] {
            for i in 0..dimensions[0] {
                poll(&mut cancelled)?;
                let indices = [i, j];
                let mid: [f64; 2] = core::array::from_fn(|d| {
                    0.5 * (coordinate_edges[d][indices[d]] + coordinate_edges[d][indices[d] + 1])
                });
                let half: [f64; 2] = core::array::from_fn(|d| {
                    0.5 * (coordinate_edges[d][indices[d] + 1] - coordinate_edges[d][indices[d]])
                });
                let edges = patch.half_edges_m();
                let center = core::array::from_fn(|d| {
                    patch.center_m()[d] + mid[0] * edges[0][d] + mid[1] * edges[1][d]
                });
                let half_edges = core::array::from_fn(|d| edges[d].map(|x| x * half[d]));
                // Reuse the physical geometry admission, including tiny-edge
                // degeneracy. Parent-area fractions share one exact Jacobian.
                MannInflowPatch::new(center, half_edges)?;
                let area_m2 = (patch.area_m2() * half[0]) * half[1];
                if !area_m2.is_finite() || area_m2 <= 0.0 {
                    return Err(numerical());
                }
                let geometry = MannInflowFace {
                    center_m: center,
                    half_edges_m: half_edges,
                    area_m2,
                };
                let (mean_density_kg_m3, velocity_scale) = density.face_density(geometry)?;
                let filter = MannBoxFilter::prepare(
                    &frame.model,
                    [half_edges[0], half_edges[1], [0.0; 3]],
                    filter_bytes,
                    &mut cancelled,
                )?;
                result.faces.push(PreparedFace {
                    geometry,
                    filter,
                    mean_density_kg_m3,
                    velocity_scale,
                });
                result.samples.push(EMPTY);
                result.trial.push(EMPTY);
            }
        }
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-inlet-plane-budget",
                "allocated inlet capacity exceeds its byte budget",
            ));
        }
        result.refresh(frame, &mut cancelled)?;
        Ok(result)
    }

    /// Atomically replace all accepted face data from one frame. Success
    /// allocates nothing and costs O(faces*modes), using the original analytic
    /// sampler. Cancellation is polled at each face (at most MAX_MODES modes)
    /// and immediately before publication. Refusal/cancellation retains every
    /// prior sample, total mass flow and time; the source frame is read-only.
    pub fn refresh(
        &mut self,
        frame: &MannInflowFrame,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        poll(&mut cancelled)?;
        if self.patch != frame.patch() {
            return Err(error(
                "mann-inlet-plane-source-mismatch",
                "inlet frame uses a different physical patch",
            ));
        }
        self.faces[0].filter.gains(&frame.model)?;
        let mut total = 0.0;
        let mut correction = 0.0;
        for (face, trial) in self.faces.iter().zip(&mut self.trial) {
            poll(&mut cancelled)?;
            let velocity = frame
                .sample_filtered(face.geometry.center_m, &face.filter)?
                .u;
            let favre_velocity_m_s = velocity.map(|u| u * face.velocity_scale);
            let mean_mass_flux_kg_m2_s = velocity.map(|u| u * self.reference_density_kg_m3);
            let normal_mass_flux_kg_m2_s = dot(mean_mass_flux_kg_m2_s, self.patch.normal());
            let flow = normal_mass_flux_kg_m2_s * face.geometry.area_m2;
            if !favre_velocity_m_s
                .iter()
                .chain(&mean_mass_flux_kg_m2_s)
                .chain([normal_mass_flux_kg_m2_s, flow].iter())
                .all(|x| x.is_finite())
            {
                return Err(numerical());
            }
            *trial = MannInflowFaceSample {
                mean_density_kg_m3: face.mean_density_kg_m3,
                favre_velocity_m_s,
                mean_mass_flux_kg_m2_s,
                normal_mass_flux_kg_m2_s,
            };
            let increment = flow - correction;
            let next = total + increment;
            correction = (next - total) - increment;
            total = next;
            if !total.is_finite() || !correction.is_finite() {
                return Err(numerical());
            }
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.samples, &mut self.trial);
        self.time_s = frame.time_s();
        self.turbulent_mass_flow_kg_s = total;
        Ok(())
    }

    pub const fn dimensions(&self) -> [usize; 2] {
        self.dimensions
    }
    pub const fn patch(&self) -> MannInflowPatch {
        self.patch
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Signed integral over all faces [kg/s], without residual correction.
    pub const fn turbulent_mass_flow_kg_s(&self) -> f64 {
        self.turbulent_mass_flow_kg_s
    }

    pub fn flat_index(&self, indices: [usize; 2]) -> Option<usize> {
        if indices[0] < self.dimensions[0] && indices[1] < self.dimensions[1] {
            Some(indices[0] + self.dimensions[0] * indices[1])
        } else {
            None
        }
    }

    /// Physical geometry in the same flattened order as samples().
    pub fn face(&self, flat_index: usize) -> Option<MannInflowFace> {
        self.faces.get(flat_index).map(|face| face.geometry)
    }

    pub fn samples(&self) -> &[MannInflowFaceSample] {
        &self.samples
    }

    /// Actual retained heap bytes, including nested filter gain capacities.
    pub fn workspace_bytes(&self) -> usize {
        self.faces.capacity() * core::mem::size_of::<PreparedFace>()
            + (self.samples.capacity() + self.trial.capacity())
                * core::mem::size_of::<MannInflowFaceSample>()
            + self
                .faces
                .iter()
                .map(|face| face.filter.workspace_bytes())
                .sum::<usize>()
    }
}

#[cfg(test)]
mod tests;
