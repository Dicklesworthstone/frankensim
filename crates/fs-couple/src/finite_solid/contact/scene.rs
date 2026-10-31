//! Prepared automatic geometry for one simultaneous rigid/solid contact system.
//!
//! Original surface and rigid-sphere owners supply the actual paired witnesses.
//! Complete explicit sphere-pair laws are bound once, independently of each
//! sphere's solid-interface law. Closed outward AABBs prune distant rigid pairs;
//! the retained rows then join the solid rows in the same mixed response.

use super::mixed::{FiniteSolidContactCounterpart, FiniteSolidMixedContact};
use super::surface::{
    FiniteSolidContactSurface, FiniteSolidSphere, FiniteSolidSurfaceControls,
    FiniteSolidSurfaceError, FiniteSolidSurfaceQuery,
};
use crate::finite_solid::FiniteSolidBody;
use fs_contact::broad_phase::{
    SweepConfig, SweepError, SweepItem, sweep_candidates, swept_sphere_aabb,
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_mbd::frictional_impact::spheres::{SphereContactError, sphere_pair_witness};
use fs_mbd::{RigidBodyState, Vec3};

/// Uses the same explicitly resolved ordered-pair law as the native rigid owner.
pub use fs_mbd::frictional_impact::spheres::SpherePairLaw;

/// Original feature identities, independent of input-array order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteSolidSphereSceneKey {
    Solid {
        sphere: u64,
        /// Original surface contact ID, before the disjoint scene namespace.
        contact_id: u64,
        face: usize,
    },
    Spheres {
        first: u64,
        second: u64,
    },
}

/// Complete admitted current geometry; no state has been advanced or repaired.
#[derive(Clone, Debug)]
pub struct FiniteSolidSphereSceneQuery {
    /// Sorted by scene contact ID. Surface IDs are even; rigid-pair IDs are odd.
    pub contacts: Vec<FiniteSolidMixedContact>,
    pub keys: Vec<FiniteSolidSphereSceneKey>,
    pub maximum_contact_gap_m: f64,
    pub surface: FiniteSolidSurfaceQuery,
    pub sphere_pair_candidates: usize,
    pub sphere_pair_comparisons: usize,
}

#[derive(Debug)]
pub enum FiniteSolidSphereSceneError {
    Invalid(&'static str),
    Surface(FiniteSolidSurfaceError),
    Sphere(SphereContactError),
    BroadPhase(SweepError),
    Penetration { first: u64, second: u64, gap_m: f64 },
    ContactBudget { maximum_contacts: usize },
    ContactIdOverflow,
    Cancelled,
}
impl From<FiniteSolidSurfaceError> for FiniteSolidSphereSceneError {
    fn from(e: FiniteSolidSurfaceError) -> Self {
        Self::Surface(e)
    }
}
impl From<SphereContactError> for FiniteSolidSphereSceneError {
    fn from(e: SphereContactError) -> Self {
        Self::Sphere(e)
    }
}
impl From<SweepError> for FiniteSolidSphereSceneError {
    fn from(e: SweepError) -> Self {
        Self::BroadPhase(e)
    }
}
impl core::fmt::Display for FiniteSolidSphereSceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mixed finite-solid sphere scene: {self:?}")
    }
}
impl std::error::Error for FiniteSolidSphereSceneError {}

/// Fixed sphere topology, canonical pair laws and bounded broad-phase work.
///
/// Each native rigid body has exactly one COM-centred spherical boundary.
/// Its original full inertia is used by mechanics; no sphere-inertia shortcut
/// is substituted. Rebuilding this object changes the contact-ID namespace and
/// requires the caller to rebind any per-contact constitutive history.
#[derive(Clone, Debug)]
pub struct FiniteSolidSphereScene {
    spheres: Vec<FiniteSolidSphere>,
    pair_laws: Vec<SpherePairLaw>,
    maximum_pairs: usize,
}
impl FiniteSolidSphereScene {
    /// Bind exactly one explicit law for every unordered sphere pair.
    /// Sphere IDs need not be contiguous. Their individual surface laws are
    /// never averaged or reused as a sphere/sphere material interface.
    pub fn new(
        spheres: &[FiniteSolidSphere],
        pair_laws: &[SpherePairLaw],
        maximum_pairs: usize,
    ) -> Result<Self, FiniteSolidSphereSceneError> {
        let count = spheres
            .len()
            .checked_mul(spheres.len().saturating_sub(1))
            .and_then(|n| n.checked_div(2))
            .ok_or(FiniteSolidSphereSceneError::Invalid("pair count overflow"))?;
        if spheres.len() > 4096
            || maximum_pairs > 65_536
            || count > maximum_pairs
            || count != pair_laws.len()
        {
            return Err(FiniteSolidSphereSceneError::Invalid(
                "bounded complete sphere-pair laws required",
            ));
        }
        let mut spheres = spheres.to_vec();
        spheres.sort_unstable_by_key(|s| s.id);
        let mut occupied = vec![false; spheres.len()];
        for (i, s) in spheres.iter().enumerate() {
            if s.rigid_body >= occupied.len()
                || occupied[s.rigid_body]
                || s.center_body_m != Vec3::ZERO
                || !s.radius_m.is_finite()
                || s.radius_m <= 0.0
                || (i > 0 && spheres[i - 1].id == s.id)
            {
                return Err(FiniteSolidSphereSceneError::Invalid(
                    "one distinct positive COM-centred sphere per native body required",
                ));
            }
            occupied[s.rigid_body] = true;
        }
        let mut pair_laws = pair_laws.to_vec();
        pair_laws.sort_unstable_by_key(|l| (l.first, l.second));
        let mut index = 0;
        for (i, a) in spheres.iter().enumerate() {
            for b in &spheres[i + 1..] {
                let law = pair_laws[index];
                if law.first != a.id
                    || law.second != b.id
                    || !law.friction_coefficient.is_finite()
                    || law.friction_coefficient < 0.0
                    || !law.restitution.is_finite()
                    || !(0.0..=1.0).contains(&law.restitution)
                    || !(a.radius_m + b.radius_m).is_finite()
                {
                    return Err(FiniteSolidSphereSceneError::Invalid(
                        "one valid explicit law per canonical sphere pair required",
                    ));
                }
                index += 1;
            }
        }
        Ok(Self {
            spheres,
            pair_laws,
            maximum_pairs,
        })
    }

    pub fn spheres(&self) -> &[FiniteSolidSphere] {
        &self.spheres
    }

    /// Detect every admitted current solid and rigid-pair contact together.
    /// Broad-phase exclusions are conservative for the declared activation gap.
    /// Both analytic sphere separation and actual reconstructed material gaps
    /// must be nonnegative; no penetration tolerance is used.
    pub fn contacts(
        &self,
        surface: &FiniteSolidContactSurface,
        solid: &FiniteSolidBody<'_>,
        bodies: &[RigidBodyState],
        controls: FiniteSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteSolidSphereSceneQuery, FiniteSolidSphereSceneError> {
        poll(cx)?;
        if bodies.len() != self.spheres.len() {
            return Err(FiniteSolidSphereSceneError::Invalid(
                "scene/native body count mismatch",
            ));
        }
        let surface_query =
            surface.sphere_separations(solid, bodies, &self.spheres, controls, cx)?;
        surface_query.admitted_contacts()?;
        let mut rows = Vec::with_capacity(surface_query.contacts.len());
        let mut maximum_gap = 0.0_f64;
        for witness in &surface_query.contacts {
            let mut contact = FiniteSolidMixedContact::from(witness.contact);
            contact.id = contact
                .id
                .checked_mul(2)
                .ok_or(FiniteSolidSphereSceneError::ContactIdOverflow)?;
            rows.push((
                contact,
                FiniteSolidSphereSceneKey::Solid {
                    sphere: witness.sphere_id,
                    contact_id: witness.contact.id,
                    face: witness.face,
                },
            ));
            maximum_gap = maximum_gap.max(witness.gap_m);
        }
        let mut items = Vec::with_capacity(self.spheres.len());
        for s in &self.spheres {
            poll(cx)?;
            let p = bodies[s.rigid_body].pose().position_world();
            let xyz = [p.x, p.y, p.z];
            items.push(SweepItem {
                id: s.id,
                bounds: swept_sphere_aabb(xyz, xyz, s.radius_m, controls.activation_tolerance_m)?,
            });
        }
        let broad = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: self.maximum_pairs,
                max_comparisons: self.maximum_pairs,
            },
            || cx.checkpoint().is_err(),
        )?;
        for pair in &broad.pairs {
            poll(cx)?;
            let a = self.spheres[pair.index_a];
            let b = self.spheres[pair.index_b];
            let witness = sphere_pair_witness(
                a.rigid_body,
                bodies[a.rigid_body],
                a.radius_m,
                b.rigid_body,
                bodies[b.rigid_body],
                b.radius_m,
            )?;
            let minimum_gap = witness.centre_gap_m.min(witness.material_gap_m);
            if minimum_gap < 0.0 {
                return Err(FiniteSolidSphereSceneError::Penetration {
                    first: a.id,
                    second: b.id,
                    gap_m: minimum_gap,
                });
            }
            if witness.centre_gap_m > controls.activation_tolerance_m
                || witness.material_gap_m > controls.activation_tolerance_m
            {
                continue;
            }
            let law_index = self
                .pair_laws
                .binary_search_by_key(&(a.id, b.id), |l| (l.first, l.second))
                .map_err(|_| FiniteSolidSphereSceneError::Invalid("prepared pair law missing"))?;
            let law = self.pair_laws[law_index];
            let id = u64::try_from(law_index)
                .ok()
                .and_then(|i| i.checked_mul(2))
                .and_then(|i| i.checked_add(1))
                .ok_or(FiniteSolidSphereSceneError::ContactIdOverflow)?;
            rows.push((
                FiniteSolidMixedContact {
                    id,
                    point: witness.first,
                    counterpart: FiniteSolidContactCounterpart::Rigid(witness.second),
                    normal_world: witness.normal_world,
                    tangent_world: witness.tangent_world,
                    activation_tolerance_m: controls.activation_tolerance_m,
                    friction_coefficient: law.friction_coefficient,
                    normal_law: NormalVelocityLaw::NewtonRestitution {
                        coefficient: law.restitution,
                    },
                },
                FiniteSolidSphereSceneKey::Spheres {
                    first: a.id,
                    second: b.id,
                },
            ));
            maximum_gap = maximum_gap.max(witness.material_gap_m);
            if rows.len() > controls.maximum_contacts {
                return Err(FiniteSolidSphereSceneError::ContactBudget {
                    maximum_contacts: controls.maximum_contacts,
                });
            }
        }
        rows.sort_unstable_by_key(|(c, _)| c.id);
        let (contacts, keys) = rows.into_iter().unzip();
        Ok(FiniteSolidSphereSceneQuery {
            contacts,
            keys,
            maximum_contact_gap_m: maximum_gap,
            surface: surface_query,
            sphere_pair_candidates: broad.pairs.len(),
            sphere_pair_comparisons: broad.checked_pairs,
        })
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), FiniteSolidSphereSceneError> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidSphereSceneError::Cancelled)
}
