//! Uniform hydrostatic fluid loading on the existing material rod owner.
//!
//! Integrating `-p n` over a fully submerged closed displaced volume with
//! `grad(p) = rho_fluid g` gives `-rho_fluid V g`. The same linear centerline
//! shape functions as reference weight integrate each segment's displaced
//! volume; an attached straight strand puts half its volume at each endpoint.
//! Its rotating offsets supply moments through the SAME potential and tangent.
//! Fluid density never replaces solid mass or creates a negative inertial mass.

use fs_exec::CancelGate;
use fs_ga::{So3, Vec3};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::section::ConcentricElasticSection;

use super::super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use super::GravityLoad;
use crate::rod::Rod;

/// Explicit displaced-volume choice for all rod sections and attached strands.
/// Geometry and material are independent: a hollow wall's density cannot say
/// whether fluid enters its cavity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceDisplacementVolume {
    /// Fluid reaches central voids; only the actual annular material displaces
    /// fluid. A solid section gives the same result as `SealedEnvelope`.
    FloodedVoids,
    /// Closed voids exclude ambient fluid, so use the outer circular envelope.
    /// Any contents of the void need their own material mass; none is invented.
    SealedEnvelope,
}

/// Accepted shared fluid snapshot and its separately accounted hydrostatic load.
/// Source receipts stay on the original `FluidTransport` input. This numerical
/// snapshot does not turn authored coefficients into sourced measurements.
#[derive(Clone, Debug, PartialEq)]
pub struct SubmergedFluidLoad {
    transport: FluidTransportState,
    volume_model: ReferenceDisplacementVolume,
    segment_volumes_m3: Vec<f64>,
    link_volumes_m3: Vec<f64>,
    nodal_volumes_m3: Vec<f64>,
    nodal_forces_n: Vec<[f64; 3]>,
    link_endpoint_forces_n: Vec<[f64; 3]>,
    total_volume_m3: f64,
}

impl SubmergedFluidLoad {
    pub(super) fn try_new(
        rod: &SpatialMaterialRod,
        transport: FluidTransportState,
        volume_model: ReferenceDisplacementVolume,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        let acceleration = rod.gravity_m_s2().unwrap_or([0.0; 3]);
        let inventory = rod.reference_mass();
        let mut segment_volumes_m3 = Vec::with_capacity(inventory.sections().len());
        let mut nodal_volumes_m3 = vec![0.0; rod.positions_m().len()];
        let mut total_volume_m3 = 0.0;
        for (index, section) in inventory.sections().iter().enumerate() {
            checkpoint(gate)?;
            let volume =
                displaced_area(section.elastic(), volume_model) * inventory.segment_length_m();
            add_volume(
                volume,
                [index, index + 1],
                &mut nodal_volumes_m3,
                &mut total_volume_m3,
            )?;
            segment_volumes_m3.push(volume);
        }
        let mut link_volumes_m3 = Vec::with_capacity(rod.tension_links.len());
        let mut link_endpoint_forces_n = Vec::with_capacity(rod.tension_links.len());
        for link in &rod.tension_links {
            checkpoint(gate)?;
            let volume = displaced_area(link.section(), volume_model)
                * link.natural_length_m()
                * f64::from(link.parallel_strands());
            add_volume(
                volume,
                (*link.attachments()).map(|a| a.node),
                &mut nodal_volumes_m3,
                &mut total_volume_m3,
            )?;
            link_volumes_m3.push(volume);
            link_endpoint_forces_n.push(force(0.5 * volume, transport, acceleration)?);
        }
        let nodal_forces_n = nodal_volumes_m3
            .iter()
            .map(|&volume| force(volume, transport, acceleration))
            .collect::<Result<Vec<_>, _>>()?;
        checkpoint(gate)?;
        Ok(Self {
            transport,
            volume_model,
            segment_volumes_m3,
            link_volumes_m3,
            nodal_volumes_m3,
            nodal_forces_n,
            link_endpoint_forces_n,
            total_volume_m3,
        })
    }

    /// Original coherent shared transport state; hydrostatics consumes its density.
    pub const fn transport_state(&self) -> FluidTransportState {
        self.transport
    }
    /// Explicit flooded-void or sealed-envelope geometry policy.
    pub const fn volume_model(&self) -> ReferenceDisplacementVolume {
        self.volume_model
    }
    /// Reference volumes [m³], in original rod-segment order.
    pub fn segment_volumes_m3(&self) -> &[f64] {
        &self.segment_volumes_m3
    }
    /// Complete strand-bundle volumes [m³], in original attachment order.
    pub fn link_volumes_m3(&self) -> &[f64] {
        &self.link_volumes_m3
    }
    /// Centerline-integrated volume weights [m³], including the clamped root.
    pub fn nodal_volumes_m3(&self) -> &[f64] {
        &self.nodal_volumes_m3
    }
    /// Total displaced reference volume of the complete specimen [m³].
    pub const fn total_volume_m3(&self) -> f64 {
        self.total_volume_m3
    }
    /// Buoyant forces ON the nodes in world axes [N], excluding solid weight.
    pub fn nodal_forces_n(&self) -> &[[f64; 3]] {
        &self.nodal_forces_n
    }
    pub(super) fn link_endpoint_forces_n(&self) -> &[[f64; 3]] {
        &self.link_endpoint_forces_n
    }
}

impl SpatialMaterialRod {
    /// Bind one fully submerged, uniform-density hydrostatic environment through
    /// the SAME material transport interface used by flow and heat consumers.
    /// Absent gravity remains zero; replacing gravity or attached links rebuilds
    /// this load. Rebinding a fluid replaces its old buoyancy instead of adding
    /// another force, while `without_submerged_fluid` restores dry loading.
    ///
    /// All sections/strands are assumed fully immersed and their displaced
    /// volumes nonoverlapping. The explicit reference-volume approximation keeps
    /// volume fixed under extension, bending and heating. It does not infer
    /// radial expansion, fluid entering a sealed void, a free surface, partial
    /// immersion, capillarity, drag, added mass, pressure-dependent solid strain
    /// or a compressible-fluid hydrostatic profile. This is static construction,
    /// not a work-conserving live change of fluid. The native motion adapter
    /// sees the same conservative load but does not acquire fluid dynamics.
    ///
    /// # Errors
    /// Coherent-fluid, volume, force, potential or cancellation failures return
    /// no changed specimen. The original fluid retains its source authority.
    pub fn with_submerged_fluid(
        &self,
        fluid: &impl FluidTransport,
        volume_model: ReferenceDisplacementVolume,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        let transport = fluid.transport_state().map_err(|_| {
            SpatialRodError::Invalid(
                "submerged fluid requires a coherent admitted material transport state",
            )
        })?;
        let load = SubmergedFluidLoad::try_new(self, transport, volume_model, gate)?;
        let mut candidate = self.clone();
        candidate
            .gravity
            .get_or_insert_with(|| GravityLoad {
                acceleration_m_s2: [0.0; 3],
                force_n: vec![[0.0; 3]; self.rod.positions.len()],
                buoyancy: None,
            })
            .buoyancy = Some(load);
        candidate.equilibrium_energy_for(&candidate.rod)?;
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Remove only hydrostatic fluid loading; keep mass, gravity and the pose.
    pub fn without_submerged_fluid(&self, gate: &CancelGate) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        let mut candidate = self.clone();
        if let Some(gravity) = &mut candidate.gravity {
            gravity.buoyancy = None;
        }
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Current accepted fluid snapshot and loads; `None` means no ambient fluid.
    pub fn submerged_fluid_load(&self) -> Option<&SubmergedFluidLoad> {
        self.gravity
            .as_ref()
            .and_then(|gravity| gravity.buoyancy.as_ref())
    }

    /// Reference-relative hydrostatic potential [J], separately from weight and
    /// stored elasticity. Its negative gradient supplies every buoyant load.
    pub fn buoyancy_potential_j(&self) -> Result<f64, SpatialRodError> {
        self.buoyancy_potential_for(&self.rod)
    }

    pub(in super::super) fn buoyancy_potential_for(
        &self,
        rod: &Rod,
    ) -> Result<f64, SpatialRodError> {
        let Some(load) = self.submerged_fluid_load() else {
            return Ok(0.0);
        };
        let mut potential = 0.0;
        for (node, (position, force)) in rod.positions.iter().zip(&load.nodal_forces_n).enumerate()
        {
            let reference = self
                .reference_geometry
                .as_ref()
                .map_or([node as f64 * rod.l0, 0.0, 0.0], |geometry| {
                    geometry.positions_m()[node]
                });
            for axis in 0..3 {
                potential -= force[axis] * (position[axis] - reference[axis]);
            }
        }
        // Endpoint translation already appears above; only rotating offsets remain.
        for (link, force) in self.tension_links.iter().zip(&load.link_endpoint_forces_n) {
            for attachment in link.attachments() {
                let node = attachment.node;
                let reference = self
                    .reference_geometry
                    .as_ref()
                    .map_or(So3::identity(), |geometry| geometry.frames()[node]);
                let offset = vector(attachment.offset_body_m);
                let actual = rod.frames[node]
                    .rotate(offset)
                    .map_err(crate::rod::group_error)?;
                let original = reference.rotate(offset).map_err(crate::rod::group_error)?;
                potential -= vector(*force).dot(Vec3::new(
                    actual.x - original.x,
                    actual.y - original.y,
                    actual.z - original.z,
                ));
            }
        }
        if !potential.is_finite() {
            return Err(SpatialRodError::Invalid(
                "hydrostatic potential is unrepresentable",
            ));
        }
        Ok(potential)
    }

    /// Attached-strand buoyant moments ON the nodes, in their body frames [N m].
    /// This excludes weight and duplicates no centerline force.
    pub fn buoyancy_nodal_moments_body_nm(&self) -> Result<Option<Vec<[f64; 3]>>, SpatialRodError> {
        let Some(load) = self.submerged_fluid_load() else {
            return Ok(None);
        };
        let mut moments = vec![[0.0; 3]; self.rod.positions.len()];
        for (link, force) in self.tension_links.iter().zip(&load.link_endpoint_forces_n) {
            for attachment in link.attachments() {
                let local = self.rod.frames[attachment.node]
                    .inverse()
                    .rotate(vector(*force))
                    .map_err(crate::rod::group_error)?;
                let moment = vector(attachment.offset_body_m).cross(local);
                for (sum, value) in moments[attachment.node]
                    .iter_mut()
                    .zip([moment.x, moment.y, moment.z])
                {
                    *sum += value;
                }
            }
        }
        if moments.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid(
                "hydrostatic attachment moment is unrepresentable",
            ));
        }
        Ok(Some(moments))
    }
}

fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

fn displaced_area(section: &ConcentricElasticSection, model: ReferenceDisplacementVolume) -> f64 {
    match model {
        ReferenceDisplacementVolume::FloodedVoids => section.area_m2(),
        ReferenceDisplacementVolume::SealedEnvelope => {
            std::f64::consts::PI * section.outer_radius_m() * section.outer_radius_m()
        }
    }
}

fn add_volume(
    volume: f64,
    nodes: [usize; 2],
    nodal: &mut [f64],
    total: &mut f64,
) -> Result<(), SpatialRodError> {
    if !volume.is_finite() || 0.5 * volume <= 0.0 {
        return Err(SpatialRodError::Invalid(
            "positive representable displaced reference volume required",
        ));
    }
    *total += volume;
    for node in nodes {
        nodal[node] += 0.5 * volume;
    }
    if !total.is_finite() || nodes.iter().any(|&node| !nodal[node].is_finite()) {
        return Err(SpatialRodError::Invalid(
            "displaced-volume sum is unrepresentable",
        ));
    }
    Ok(())
}

fn force(
    volume: f64,
    fluid: FluidTransportState,
    gravity: [f64; 3],
) -> Result<[f64; 3], SpatialRodError> {
    let displaced_mass = fluid.density_kg_m3() * volume;
    let value = gravity.map(|g| -displaced_mass * g);
    if !displaced_mass.is_finite()
        || displaced_mass <= 0.0
        || value
            .iter()
            .zip(gravity)
            .any(|(&f, g)| !f.is_finite() || (g != 0.0 && f == 0.0))
    {
        return Err(SpatialRodError::Invalid(
            "hydrostatic force is unrepresentable",
        ));
    }
    Ok(value)
}
