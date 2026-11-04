//! Explicit native lifting surfaces, normal-section loads and grouped induced drag.
//!
//! Geometry is body-fixed. A control gain changes the section polar's incidence;
//! it does not rotate panels or move a hinge. The supplied strip chord is the
//! physical chord normal to the bound span, including for swept surfaces.

use super::{AerodynamicFields, GearedPropulsion, Refusal, reject, vec, wing_error};
use fs_ga::{Vec3, Wrench};
use fs_wing::nonlinear::{SectionFrame, StripSpec};
use fs_wing::{MAX_PANELS, Panel};
use std::borrow::Cow;

/// Per-strip section and control data; no horizontal-wing defaults are inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionProperties {
    /// Camber ratio used by the existing attached/separated lift law and Cm0.
    pub camber_ratio: f64,
    /// Attached profile Cd; separated pressure drag uses the same accepted alpha.
    pub attached_cd: f64,
    /// Add this gain times `PropulsionInputs::canard_rad` to the base strip twist.
    /// Zero fixes the incidence. This is a polar control, not panel articulation.
    pub canard_gain: f64,
    /// Index of the explicitly declared classical induced-drag group.
    pub induced_drag_group: usize,
}

/// Classical engineering closure for one declared set of lifting strips.
///
/// The group's complete circulatory force is projected transverse to CG wind:
/// `D = |F_perp|² / (q_CG * pi * e * span²)`. Reference area cancels from this
/// dimensional form. Group related biplane surfaces together when one effective
/// span/efficiency is intended; give a vertical surface its own group when its
/// load must not cancel a different surface's drag. This closure excludes
/// cross-group induced-drag interference and is not a Trefftz-plane solution.
/// Drag acts opposite CG wind at the CG and adds no extra moment. The coupled
/// circulation still includes all surfaces' mutual induced velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InducedDragGroup {
    pub span_m: f64,
    pub oswald_e: f64,
}

/// Immutable borrowed geometry with admitted section frames and physical widths.
/// Each panel belongs to exactly one strip and each strip to exactly one group.
pub struct LiftingSurfaces<'a> {
    pub(super) panels: &'a [Panel],
    strips: &'a [StripSpec],
    pub(super) properties: &'a [SectionProperties],
    pub(super) groups: &'a [InducedDragGroup],
    pub(super) frames: Vec<SectionFrame>,
    widths_m: Vec<f64>,
}

impl<'a> LiftingSurfaces<'a> {
    /// Prepare geometry once; all native wake owners borrow this same model.
    ///
    /// # Errors
    /// `lifting-surfaces-invalid` for incomplete/duplicate panel or group
    /// ownership, invalid section parameters, or inconsistent chordwise rows;
    /// frame admission refusals preserve their native diagnostic.
    pub fn new(
        panels: &'a [Panel],
        strips: &'a [StripSpec],
        properties: &'a [SectionProperties],
        groups: &'a [InducedDragGroup],
    ) -> Result<Self, Refusal> {
        let invalid = || {
            reject(
                "lifting-surfaces-invalid",
                "provide complete compatible geometry, sections and positive drag references",
            )
        };
        if panels.is_empty()
            || panels.len() > MAX_PANELS
            || strips.is_empty()
            || properties.len() != strips.len()
            || groups.is_empty()
            || groups.len() > strips.len()
            || groups.iter().any(|g| {
                !g.span_m.is_finite()
                    || g.span_m <= 0.0
                    || !g.oswald_e.is_finite()
                    || g.oswald_e <= 0.0
                    || !(g.span_m * g.span_m * g.oswald_e).is_finite()
                    || g.span_m * g.span_m * g.oswald_e <= 0.0
            })
        {
            return Err(invalid());
        }
        let mut owners = vec![false; panels.len()];
        let mut group_used = vec![false; groups.len()];
        let mut frames = Vec::with_capacity(strips.len());
        let mut widths_m = Vec::with_capacity(strips.len());
        for (strip, properties) in strips.iter().zip(properties) {
            if strip.panel_indices.is_empty()
                || !strip.chord_m.is_finite()
                || strip.chord_m <= 0.0
                || !strip.twist_rad.is_finite()
                || !properties.camber_ratio.is_finite()
                || properties.camber_ratio.abs() > 0.15
                || !properties.attached_cd.is_finite()
                || properties.attached_cd < 0.0
                || !properties.canard_gain.is_finite()
                || properties.induced_drag_group >= groups.len()
            {
                return Err(invalid());
            }
            let front = panels.get(strip.panel_indices[0]).ok_or_else(invalid)?;
            let frame = SectionFrame::from_panel(front).map_err(wing_error)?;
            let width = (vec(front.b) - vec(front.a)).norm();
            let chord_line = (vec(front.a) + vec(front.b)).scale(0.5) - vec(front.ctrl);
            let chord_line = chord_line.scale(1.0 / chord_line.norm());
            if !(strip.chord_m * width).is_finite() || strip.chord_m * width <= 0.0 {
                return Err(invalid());
            }
            for &index in &strip.panel_indices {
                let panel = panels.get(index).ok_or_else(invalid)?;
                frame.admit_panel(panel).map_err(wing_error)?;
                let row_width = (vec(panel.b) - vec(panel.a)).norm();
                // Register the chordwise rows along their physical panel chord.
                // Its span component is retained for swept streamwise layouts.
                let offset = vec(panel.a) - vec(front.a);
                let mismatch = offset - chord_line.scale(offset.dot(chord_line));
                if owners[index]
                    || panel.surface != front.surface
                    || !panel.width_m.is_finite()
                    || panel.width_m <= 0.0
                    || (row_width - width).abs() > 1e-9 * width
                    || mismatch.norm() > 1e-9 * (width + offset.norm())
                {
                    return Err(invalid());
                }
                owners[index] = true;
            }
            group_used[properties.induced_drag_group] = true;
            frames.push(frame);
            widths_m.push(width);
        }
        if owners.contains(&false) || group_used.contains(&false) {
            return Err(invalid());
        }
        Ok(Self {
            panels,
            strips,
            properties,
            groups,
            frames,
            widths_m,
        })
    }

    pub(super) fn controlled_strips(&self, control: f64) -> Result<Cow<'_, [StripSpec]>, Refusal> {
        if control == 0.0
            || self
                .properties
                .iter()
                .all(|section| section.canard_gain == 0.0)
        {
            return Ok(Cow::Borrowed(self.strips));
        }
        let mut strips = self.strips.to_vec();
        for (strip, properties) in strips.iter_mut().zip(self.properties) {
            strip.twist_rad += properties.canard_gain * control;
            if !strip.twist_rad.is_finite() {
                return Err(reject(
                    "lifting-surfaces-invalid",
                    "section incidence overflow",
                ));
            }
        }
        Ok(Cow::Owned(strips))
    }

    pub(super) fn additional_wrench(
        &self,
        rho: f64,
        cg_air: Vec3,
        origin: [f64; 3],
        fields: &AerodynamicFields<'_>,
        group_forces: &[Vec3],
    ) -> Result<Wrench, Refusal> {
        let local = fields.strip_velocity_mps.ok_or_else(|| {
            reject(
                "lifting-surfaces-invalid",
                "accepted section velocities are missing",
            )
        })?;
        if local.len() != self.strips.len()
            || fields.strip_alpha_eff_rad.len() != self.strips.len()
            || !origin.iter().all(|v| v.is_finite())
        {
            return Err(reject(
                "lifting-surfaces-invalid",
                "one local velocity and effective incidence per strip and a finite wrench origin are required",
            ));
        }
        let mut force = Vec3::default();
        let mut moment = Vec3::default();
        for (s, strip) in self.strips.iter().enumerate() {
            let frame = self.frames[s];
            let flow = frame.flow(local[s]).map_err(wing_error)?;
            let q = 0.5 * rho * flow.speed_mps * flow.speed_mps;
            if !q.is_finite() || !fields.strip_alpha_eff_rad[s].is_finite() {
                return Err(reject(
                    "lifting-surfaces-invalid",
                    "nonfinite strip dynamic pressure or effective incidence",
                ));
            }
            let section = self.properties[s];
            let cd = crate::sections::section_profile_drag(
                fields.strip_alpha_eff_rad[s],
                section.attached_cd,
            );
            let area = strip.chord_m * self.widths_m[s];
            let profile_force = vec(flow.velocity_mps).scale(-q * area * cd / flow.speed_mps);
            let panel = self.panels[strip.panel_indices[0]];
            let arm = (vec(panel.a) + vec(panel.b)).scale(0.5) - vec(origin);
            let cm0 = -core::f64::consts::PI * section.camber_ratio;
            force = force + profile_force;
            moment = moment
                + arm.cross(profile_force)
                + vec(frame.span_direction()).scale(q * area * strip.chord_m * cm0);
        }
        // The induced-drag closure has a finite, strictly positive CG
        // dynamic-pressure domain. Local strip flow may be admissible while
        // the CG reference is not; never turn a zero/invalid reference speed
        // into an unbounded or NaN load.
        let speed2 = cg_air.dot(cg_air);
        let q = 0.5 * rho * speed2;
        if !rho.is_finite()
            || rho <= 0.0
            || !cg_air.to_array().iter().all(|v| v.is_finite())
            || !speed2.is_finite()
            || speed2 <= 1.0
            || !q.is_finite()
            || group_forces.len() != self.groups.len()
        {
            return Err(reject(
                "lifting-surfaces-invalid",
                "finite positive density, CG airspeed above 1 m/s and one force per induced-drag group are required",
            ));
        }
        let direction = cg_air.scale(1.0 / speed2.sqrt());
        for (group, &load) in self.groups.iter().zip(group_forces) {
            if !load.to_array().iter().all(|v| v.is_finite()) {
                return Err(reject("lifting-surfaces-invalid", "nonfinite grouped circulatory load"));
            }
            let transverse = load - direction.scale(load.dot(direction));
            let drag = transverse.dot(transverse)
                / (q * core::f64::consts::PI * group.oswald_e * group.span_m * group.span_m);
            if !drag.is_finite() || drag < 0.0 {
                return Err(reject("lifting-surfaces-invalid", "unrepresentable induced drag"));
            }
            force = force - direction.scale(drag);
        }
        if !force.to_array().iter().all(|v| v.is_finite())
            || !moment.to_array().iter().all(|v| v.is_finite())
        {
            return Err(reject("lifting-surfaces-invalid", "unrepresentable total aerodynamic wrench"));
        }
        Ok(Wrench::new(moment, force))
    }
}

impl<'a, const N: usize> GearedPropulsion<'a, N> {
    /// Select an admitted custom layout for every steady, swirling, dynamic,
    /// harmonic and explicitly continued native aerodynamic stage.
    /// The original design still supplies rotors, disk locations, CG and the
    /// parasite ledger; its built-in lifting layout is replaced in this owner.
    /// CG speed must exceed 1 m/s, each section must have positive local
    /// chordwise flow, and each disk must meet its native projected-flow domain.
    /// No global body-x incidence or velocity sign is imposed on this layout.
    ///
    /// # Errors
    /// `lifting-surfaces-invalid` unless physical three-dimensional disk wash
    /// was selected with `new_spatial`.
    pub fn with_lifting_surfaces(
        mut self,
        surfaces: &'a LiftingSurfaces<'a>,
    ) -> Result<Self, Refusal> {
        if self.spatial_axes.is_none() {
            return Err(reject(
                "lifting-surfaces-invalid",
                "custom surfaces require spatial disk wash",
            ));
        }
        self.surfaces = Some(surfaces);
        Ok(self)
    }
}
