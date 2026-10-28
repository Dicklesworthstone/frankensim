//! Material tensile members attached to the same deforming Cosserat rod.
//!
//! Each link represents identical straight, pin-ended strands with a common
//! natural length. Its energy is EA/(2 L0) max(L-L0,0)^2. Attachment offsets
//! rotate with the original nodal frames, so the original energy gradient also
//! supplies their moments. This is a quasistatic taut-link reduction, not a
//! sagging cable, a dynamic hair ribbon, or a friction/thermal-history model.

use fs_exec::CancelGate;
use fs_ga::{So3, Vec3};
use fs_material::section::ConcentricElasticSection;
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use fs_math::det;

use super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use crate::rod::{MAX_MATERIAL_ROD_SEGMENTS, Rod};

/// Material-frame point attached rigidly to an existing node, including root.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensionLinkAttachment {
    /// Node in 0..=segment_count.
    pub node: usize,
    /// Offset from that node's centerline in its rotating body frame [m].
    pub offset_body_m: [f64; 3],
}

/// Immutable source, geometry and strain domain of one tensile bundle.
/// Density and Young's modulus belong to the explicitly fixed material state.
/// Heating the rod does not silently heat this separate attached member.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialTensionLink {
    section: ConcentricElasticSection,
    attachments: [TensionLinkAttachment; 2],
    natural_length_m: f64,
    parallel_strands: u32,
    maximum_tensile_strain: f64,
    axial_stiffness_n: f64,
    reference_mass_kg: f64,
}

/// Physical observations at the actual rod configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct TensionLinkState {
    /// Current endpoint coordinates [m], in attachment order.
    pub endpoints_m: [[f64; 3]; 2],
    /// Current endpoint distance [m].
    pub length_m: f64,
    /// max((length-natural_length)/natural_length, 0).
    pub tensile_strain: f64,
    /// Nonnegative bundle force magnitude [N].
    pub tension_n: f64,
    /// Stored energy of the complete bundle [J].
    pub stored_energy_j: f64,
    /// Equal-and-opposite forces exerted ON the rod [N], in the world frame.
    pub force_n: [[f64; 3]; 2],
    /// Corresponding offset moments ON the rod [N m], in each node's body frame.
    pub moment_body_nm: [[f64; 3]; 2],
}

fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn rotate(frame: So3, a: [f64; 3]) -> Result<[f64; 3], SpatialRodError> {
    Ok(array(
        frame.rotate(vector(a)).map_err(crate::rod::group_error)?,
    ))
}

impl MaterialTensionLink {
    /// Bind real layered section data to a bounded, tension-only linear law.
    /// `parallel_strands` multiplies section EA and reference mass exactly once.
    /// No installed-length density rescaling or material-name preset is used.
    ///
    /// # Errors
    /// Refuses coincident node IDs, invalid geometry/count/strain limits,
    /// nonstatic selected claims, or unrepresentable stiffness/reference mass.
    pub fn try_new(
        section: ConcentricElasticSection,
        attachments: [TensionLinkAttachment; 2],
        natural_length_m: f64,
        parallel_strands: u32,
        maximum_tensile_strain: f64,
    ) -> Result<Self, SpatialRodError> {
        if attachments[0].node == attachments[1].node
            || attachments
                .iter()
                .flat_map(|a| a.offset_body_m)
                .any(|v| !v.is_finite())
            || !natural_length_m.is_finite()
            || natural_length_m <= 0.0
            || parallel_strands == 0
            || !maximum_tensile_strain.is_finite()
            || maximum_tensile_strain <= 0.0
            || maximum_tensile_strain >= 1.0
        {
            return Err(SpatialRodError::Invalid(
                "tension link needs distinct nodes, finite offsets, positive natural length/count and a strain limit in (0,1)",
            ));
        }
        // The same resolved-property admission as the static rod: constants or
        // fixed-temperature curves, with zero included in any frequency band.
        for layer in section.layers() {
            for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY] {
                let property = layer
                    .material
                    .property(name)
                    .ok_or(SpatialRodError::Invalid(
                        "tension link lacks a required material property",
                    ))?;
                if !property.is_constant_at_fixed_temperature()
                    || property
                        .angular_frequency_band_rad_s()?
                        .is_some_and(|(lo, hi)| lo > 0.0 || hi < 0.0)
                {
                    return Err(SpatialRodError::Invalid(
                        "tension link properties must support a frozen static elastic law",
                    ));
                }
            }
        }
        let axial_stiffness_n = section.axial_stiffness_n() * f64::from(parallel_strands);
        let reference_mass_kg =
            section.line_density_kg_m() * natural_length_m * f64::from(parallel_strands);
        if !axial_stiffness_n.is_finite()
            || axial_stiffness_n <= 0.0
            || !reference_mass_kg.is_finite()
            || 0.5 * reference_mass_kg <= 0.0
        {
            return Err(SpatialRodError::Invalid(
                "tension link stiffness or reference mass is unrepresentable",
            ));
        }
        Ok(Self {
            section,
            attachments,
            natural_length_m,
            parallel_strands,
            maximum_tensile_strain,
            axial_stiffness_n,
            reference_mass_kg,
        })
    }

    /// Original concentric source and all selected claim receipts.
    pub const fn section(&self) -> &ConcentricElasticSection {
        &self.section
    }
    /// Ordered attachment geometry.
    pub const fn attachments(&self) -> &[TensionLinkAttachment; 2] {
        &self.attachments
    }
    /// Authored unloaded strand length [m].
    pub const fn natural_length_m(&self) -> f64 {
        self.natural_length_m
    }
    /// Number of identical strands in parallel.
    pub const fn parallel_strands(&self) -> u32 {
        self.parallel_strands
    }
    /// Explicit linear tensile-strain domain.
    pub const fn maximum_tensile_strain(&self) -> f64 {
        self.maximum_tensile_strain
    }
    /// Reference bundle mass [kg], independent of current span and rod temperature.
    pub const fn reference_mass_kg(&self) -> f64 {
        self.reference_mass_kg
    }

    pub(super) fn state_for(&self, rod: &Rod) -> Result<TensionLinkState, SpatialRodError> {
        let mut endpoints_m = [[0.0; 3]; 2];
        for (i, attachment) in self.attachments.iter().enumerate() {
            let position = rod
                .positions
                .get(attachment.node)
                .ok_or(SpatialRodError::Invalid(
                    "tension link attachment is outside the rod",
                ))?;
            let offset = rotate(rod.frames[attachment.node], attachment.offset_body_m)?;
            endpoints_m[i] = std::array::from_fn(|c| position[c] + offset[c]);
        }
        let delta: [f64; 3] = std::array::from_fn(|c| endpoints_m[1][c] - endpoints_m[0][c]);
        let length_m = det::hypot(det::hypot(delta[0], delta[1]), delta[2]);
        if !length_m.is_finite() || length_m <= 0.0 {
            return Err(SpatialRodError::Invalid(
                "tension link current endpoint separation must be finite and positive",
            ));
        }
        let tensile_strain = ((length_m - self.natural_length_m) / self.natural_length_m).max(0.0);
        if !tensile_strain.is_finite() || tensile_strain > self.maximum_tensile_strain {
            return Err(SpatialRodError::Invalid(
                "tension link exceeds its explicit linear tensile-strain domain",
            ));
        }
        let tension_n = self.axial_stiffness_n * tensile_strain;
        let stored_energy_j = 0.5 * tension_n * (self.natural_length_m * tensile_strain);
        let a: [f64; 3] = delta.map(|d| tension_n * (d / length_m));
        let force_n = [a, a.map(|v| -v)];
        let mut moment_body_nm = [[0.0; 3]; 2];
        for (i, attachment) in self.attachments.iter().enumerate() {
            let local = rotate(rod.frames[attachment.node].inverse(), force_n[i])?;
            moment_body_nm[i] = array(vector(attachment.offset_body_m).cross(vector(local)));
        }
        if !tension_n.is_finite()
            || !stored_energy_j.is_finite()
            || force_n
                .iter()
                .flatten()
                .chain(moment_body_nm.iter().flatten())
                .any(|v| !v.is_finite())
        {
            return Err(SpatialRodError::Invalid(
                "tension link force, moment or stored energy is unrepresentable",
            ));
        }
        Ok(TensionLinkState {
            endpoints_m,
            length_m,
            tensile_strain,
            tension_n,
            stored_energy_j,
            force_n,
            moment_body_nm,
        })
    }
}

impl SpatialMaterialRod {
    /// Replace this static specimen's attached tensile members (at most 64).
    /// Every trial in the original rod equilibrium sees their actual endpoint
    /// geometry. The same source inventory also supplies explicit gravity.
    /// This is specimen construction, not an energy-conserving live replacement.
    ///
    /// # Errors
    /// Count/index, source-domain, arithmetic or cancellation failures leave the
    /// original specimen unchanged. Existing support and load fields are retained.
    pub fn with_tension_links(
        &self,
        links: Vec<MaterialTensionLink>,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if links.len() > MAX_MATERIAL_ROD_SEGMENTS {
            return Err(SpatialRodError::Invalid(
                "static material rod admits at most 64 tension links",
            ));
        }
        for link in &links {
            checkpoint(gate)?;
            link.state_for(&self.rod)?;
        }
        let mut candidate = self.clone();
        candidate.tension_links = links;
        if !candidate.total_reference_mass_kg().is_finite() {
            return Err(SpatialRodError::Invalid(
                "rod and tension-link reference mass sum is unrepresentable",
            ));
        }
        candidate.energy_j()?;
        if let Some(acceleration) = self.gravity_m_s2() {
            candidate = candidate.with_gravity(acceleration, gate)?;
        }
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Original immutable links in authored order.
    pub fn tension_links(&self) -> &[MaterialTensionLink] {
        &self.tension_links
    }
    /// Tension, storage and reciprocal endpoint forces at the accepted pose.
    pub fn tension_link_states(&self) -> Result<Vec<TensionLinkState>, SpatialRodError> {
        self.tension_links
            .iter()
            .map(|link| link.state_for(&self.rod))
            .collect()
    }
    /// Rod plus attached-member reference inventory [kg]. `reference_mass()`
    /// retains the original rod-only receipt; links retain their own sources.
    pub fn total_reference_mass_kg(&self) -> f64 {
        self.reference_mass.total_mass_kg()
            + self
                .tension_links
                .iter()
                .map(MaterialTensionLink::reference_mass_kg)
                .sum::<f64>()
    }

    pub(super) fn link_gravity_offset_potential(
        &self,
        rod: &Rod,
        acceleration_m_s2: [f64; 3],
    ) -> Result<f64, SpatialRodError> {
        let mut potential = 0.0;
        for link in &self.tension_links {
            let force = acceleration_m_s2.map(|g| 0.5 * link.reference_mass_kg * g);
            for attachment in &link.attachments {
                let reference = self
                    .reference_geometry
                    .as_ref()
                    .map_or(So3::identity(), |g| g.frames()[attachment.node]);
                let actual = rotate(rod.frames[attachment.node], attachment.offset_body_m)?;
                let original = rotate(reference, attachment.offset_body_m)?;
                for c in 0..3 {
                    potential -= force[c] * (actual[c] - original[c]);
                }
            }
        }
        Ok(potential)
    }

    /// Body-frame moments [N m] from attached-member weight at rotating offsets.
    /// The force getter includes the same link masses; no weight is added twice.
    pub fn gravity_nodal_moments_body_nm(&self) -> Result<Option<Vec<[f64; 3]>>, SpatialRodError> {
        let Some(acceleration) = self.gravity_m_s2() else {
            return Ok(None);
        };
        let mut moments = vec![[0.0; 3]; self.rod.positions.len()];
        for link in &self.tension_links {
            let force = acceleration.map(|g| 0.5 * link.reference_mass_kg * g);
            for attachment in &link.attachments {
                let local = rotate(self.rod.frames[attachment.node].inverse(), force)?;
                let moment = array(vector(attachment.offset_body_m).cross(vector(local)));
                for (sum, value) in moments[attachment.node].iter_mut().zip(moment) {
                    *sum += value;
                }
            }
        }
        if moments.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid(
                "tension link gravity moment is unrepresentable",
            ));
        }
        Ok(Some(moments))
    }
}
