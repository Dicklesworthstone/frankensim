//! Circulatory loads from the original articulated map and nonlinear wing owner.
//!
//! Each strip is rigid in one native link. Actual link translation, rotation
//! and inherited joint motion determine its material-point air velocity.
//! Geometry, incident air and the original Kutta--Joukowsky wrench assembly
//! use one world Cartesian frame. No second motion or aerodynamic solver is
//! introduced, and fluid inertia remains in the separate added-mass owner.

use super::MAX_ARTICULATED_LINKS;
use crate::Refusal;
use fs_ga::{Vec3, Wrench};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_wing::nonlinear::{
    InfluenceOperator, NonlinearReport, SectionClosure, SectionFrame, StripSpec,
    solve_nonlinear_with_panel_incident,
};
use fs_wing::wrench::CirculatoryWrenchAccumulator;
use fs_wing::{MAX_PANELS, Panel, PanelIncidentFlow};

/// Native panel/strip geometry, with one original articulated link per strip.
///
/// Panel coordinates are local to their owning link. Every panel has exactly
/// one strip owner; all rows within that strip move with the same link. Chord
/// and twist retain their native normal-section meaning. Moving joints supply
/// geometry and velocity through the admitted map, not through a polar gain.
pub struct LinkedLiftingSurfaces<'a> {
    panels: &'a [Panel],
    strips: &'a [StripSpec],
    strip_links: &'a [usize],
    panel_strips: Vec<usize>,
    frames: Vec<SectionFrame>,
}

/// One immutable evaluation at the supplied native articulated stage.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedCirculation {
    /// Original nonlinear circulation, effective incidences and branch receipt.
    /// Force components use world axes; total_lift_n retains the native -z
    /// projection and is not an aircraft wind-axis lift measurement.
    pub nonlinear: NonlinearReport,
    /// One WORLD wrench per native link, torque about that link's world origin.
    /// Move these directly into ArticulatedLoads::world_wrenches, adding any
    /// other external wrenches about the same origins exactly once. The native
    /// generalized owner projects the resulting joint efforts by virtual work.
    pub world_wrenches: Vec<Wrench>,
    /// Actual front bound-row midpoint sampled for each native strip [m].
    pub strip_points_world_m: Vec<[f64; 3]>,
    /// Material-point velocity minus prescribed fluid velocity, in world axes
    /// [m/s]. This representative velocity drives the strip's section closure;
    /// the actual row velocities below drive its boundary loading and forces.
    pub strip_relative_velocity_world_mps: Vec<[f64; 3]>,
    /// Actual material-relative velocity at every transformed control point,
    /// in native panel order [m/s]. These determine the chord-row loading.
    pub panel_control_relative_velocity_world_mps: Vec<[f64; 3]>,
    /// Actual material-relative velocity at every transformed bound midpoint,
    /// in native panel order [m/s]. These determine panel forces and moments.
    pub panel_bound_relative_velocity_world_mps: Vec<[f64; 3]>,
}

impl<'a> LinkedLiftingSurfaces<'a> {
    /// Bind link-local geometry and prepare its original native section frames.
    /// Link existence is checked against the supplied stage map on evaluation.
    ///
    /// # Errors
    /// Invalid/incomplete strip ownership, nonfinite chord/twist, native panel
    /// or section-frame refusal, or the existing panel/link count caps.
    pub fn new(
        panels: &'a [Panel],
        strips: &'a [StripSpec],
        strip_links: &'a [usize],
    ) -> Result<Self, Refusal> {
        if panels.is_empty()
            || panels.len() > MAX_PANELS
            || strips.is_empty()
            || strips.len() > panels.len()
            || strips.len() != strip_links.len()
            || strip_links
                .iter()
                .any(|&link| link >= MAX_ARTICULATED_LINKS)
        {
            return Err(invalid(
                "bounded panels and one native link per strip required",
            ));
        }
        let mut panel_strips = vec![usize::MAX; panels.len()];
        let mut frames = Vec::with_capacity(strips.len());
        for (s, strip) in strips.iter().enumerate() {
            if strip.panel_indices.is_empty()
                || !strip.chord_m.is_finite()
                || strip.chord_m <= 0.0
                || !strip.twist_rad.is_finite()
            {
                return Err(invalid(
                    "nonempty strips with finite positive chord and finite twist required",
                ));
            }
            let first = panels
                .get(strip.panel_indices[0])
                .ok_or_else(|| invalid("strip refers to a missing panel"))?;
            let frame = SectionFrame::from_panel(first).map_err(wing_error)?;
            for &j in &strip.panel_indices {
                let panel = panels
                    .get(j)
                    .ok_or_else(|| invalid("strip refers to a missing panel"))?;
                if panel_strips[j] != usize::MAX {
                    return Err(invalid("a panel cannot belong to multiple moving strips"));
                }
                frame.admit_panel(panel).map_err(wing_error)?;
                panel_strips[j] = s;
            }
            frames.push(frame);
        }
        if panel_strips.contains(&usize::MAX) {
            return Err(invalid("every panel needs one moving strip owner"));
        }
        Ok(Self {
            panels,
            strips,
            strip_links,
            panel_strips,
            frames,
        })
    }

    /// Evaluate actual articulated circulation without changing state or history.
    ///
    /// All inputs describe the same physical stage. The caller's complete world
    /// fluid field is sampled at each control point and bound midpoint. In link
    /// axes the relative velocity is v_link - R^-1 U + omega_link cross r;
    /// inherited serial-joint motion is already in the original link twist.
    /// Subtracting the translational flow before adding the rotational term
    /// avoids erasing a small local rotational slip under a common large drift.
    ///
    /// wake_reference_velocity_world_mps is a declared common relative-air
    /// vector. It sets the original straight trailing-leg geometry; actual
    /// control velocities determine linear row ratios. Native
    /// positive chordwise-flow, conditioning and nonlinear-iteration admissions
    /// remain authoritative. Geometry is rebuilt at this stage, so no stale
    /// operator can survive a changed link pose. An optional accepted warm
    /// circulation uses the original fixed panel ownership.
    ///
    /// The section polar still uses its representative front-bound velocity;
    /// resolved row loading does not supply a new unsteady/pitch-rate polar.
    /// This quasi-steady circulatory model excludes profile/induced-drag
    /// closures, prescribed wake roll-up, added mass and unsteady pressure.
    /// It adds no hinge effort separately from the physical link wrenches.
    /// Cancellation is checked during staging and before returning; the bounded
    /// native influence factorization/nonlinear solve keeps its original API.
    /// Callback side effects are not rolled back. Use this directly inside the
    /// existing articulated midpoint load callback, alongside its fluid blocks.
    ///
    /// # Errors
    /// Missing links, incompatible geometry/air, native solve refusal,
    /// nonfinite transformed values or cancellation. No partial receipt returns.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        kinematics: &GeneralizedKinematics,
        wake_reference_velocity_world_mps: [f64; 3],
        density_kg_m3: f64,
        closure: SectionClosure<'_>,
        mut wind_world: impl FnMut([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm_start: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedCirculation, Refusal> {
        checkpoint(&mut cancelled)?;
        let links = kinematics.link_count();
        if links > MAX_ARTICULATED_LINKS || self.strip_links.iter().any(|&link| link >= links) {
            return Err(invalid(
                "surface refers to a missing or unsupported native link",
            ));
        }
        let mut panels = Vec::with_capacity(self.panels.len());
        let mut control_incident = Vec::with_capacity(self.panels.len());
        let mut bound_incident = Vec::with_capacity(self.panels.len());
        for (j, panel) in self.panels.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            let link = self.strip_links[self.panel_strips[j]];
            let pose = kinematics.world_from_link(link);
            let rotation = pose.rotation();
            let velocity = kinematics.body_twist(link);
            let point = |value| {
                pose.transform_point(vector(value))
                    .map(array)
                    .map_err(geometry)
            };
            panels.push(Panel {
                a: point(panel.a)?,
                b: point(panel.b)?,
                ctrl: point(panel.ctrl)?,
                normal: array(
                    pose.rotation()
                        .rotate(vector(panel.normal))
                        .map_err(geometry)?,
                ),
                ..*panel
            });
            let mut relative_at = |local: Vec3| -> Result<[f64; 3], Refusal> {
                let position = pose.transform_point(local).map_err(geometry)?;
                let wind = vector(wind_world(array(position))?);
                let slip = velocity.linear - rotation.inverse().rotate(wind).map_err(geometry)?
                    + velocity.angular.cross(local);
                let relative = rotation.rotate(slip).map_err(geometry)?;
                if !array(position)
                    .into_iter()
                    .chain(array(relative))
                    .all(f64::is_finite)
                {
                    return Err(invalid(
                        "unrepresentable material point or relative air velocity",
                    ));
                }
                Ok(array(relative))
            };
            control_incident.push(relative_at(vector(panel.ctrl))?);
            bound_incident.push(relative_at((vector(panel.a) + vector(panel.b)).scale(0.5))?);
        }
        let mut frames = Vec::with_capacity(self.strips.len());
        let mut points = Vec::with_capacity(self.strips.len());
        let mut incident = Vec::with_capacity(self.strips.len());
        for (s, strip) in self.strips.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            let link = self.strip_links[s];
            let pose = kinematics.world_from_link(link);
            let rotation = pose.rotation();
            let panel = self.panels[strip.panel_indices[0]];
            let local_point = (vector(panel.a) + vector(panel.b)).scale(0.5);
            let point = pose.transform_point(local_point).map_err(geometry)?;
            let local_frame = self.frames[s];
            frames.push(
                SectionFrame::new(
                    array(
                        rotation
                            .rotate(vector(local_frame.chord_direction()))
                            .map_err(geometry)?,
                    ),
                    array(
                        rotation
                            .rotate(vector(local_frame.span_direction()))
                            .map_err(geometry)?,
                    ),
                )
                .map_err(wing_error)?,
            );
            points.push(array(point));
            incident.push(bound_incident[strip.panel_indices[0]]);
        }
        checkpoint(&mut cancelled)?;
        let panel_incident = PanelIncidentFlow {
            control_mps: &control_incident,
            bound_mps: &bound_incident,
        };
        let operator = InfluenceOperator::build_with_panel_incident(
            &panels,
            wake_reference_velocity_world_mps,
            density_kg_m3,
            panel_incident,
        )
        .map_err(wing_error)?;
        let nonlinear = solve_nonlinear_with_panel_incident(
            &operator,
            &panels,
            self.strips,
            wake_reference_velocity_world_mps,
            density_kg_m3,
            closure,
            &frames,
            &incident,
            panel_incident,
            warm_start,
        )
        .map_err(wing_error)?;
        checkpoint(&mut cancelled)?;
        let mut owners = (0..links)
            .map(|link| {
                CirculatoryWrenchAccumulator::new(
                    density_kg_m3,
                    array(kinematics.world_from_link(link).translation()),
                )
                .map_err(wing_error)
            })
            .collect::<Result<Vec<_>, _>>()?;
        // Preserve native strip/panel accumulation order. Each physical force
        // has its own transformed bound midpoint and owning link-origin arm.
        for (s, strip) in self.strips.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            let owner = &mut owners[self.strip_links[s]];
            for &j in &strip.panel_indices {
                owner
                    .add_panel(&panels[j], nonlinear.gamma[j], bound_incident[j])
                    .map_err(wing_error)?;
            }
        }
        checkpoint(&mut cancelled)?;
        Ok(ArticulatedCirculation {
            nonlinear,
            world_wrenches: owners.into_iter().map(|owner| owner.total()).collect(),
            strip_points_world_m: points,
            strip_relative_velocity_world_mps: incident,
            panel_control_relative_velocity_world_mps: control_incident,
            panel_bound_relative_velocity_world_mps: bound_incident,
        })
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn geometry(error: fs_ga::GaError) -> Refusal {
    invalid(&error.to_string())
}
fn wing_error(error: fs_wing::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}
fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "articulated-circulation-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply matching link-local surfaces, native stage motion and physical incident air"
                .into(),
        ],
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(Refusal {
            code: "articulated-circulation-cancelled",
            ..invalid("circulatory load evaluation cancelled")
        })
    } else {
        Ok(())
    }
}
