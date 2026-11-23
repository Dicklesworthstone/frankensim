//! A field stored in a uniform gap of variable width between a fixed
//! boundary and a guided one: the dielectric gap of a parallel-plate
//! capacitor, or the air gap of a magnetic circuit whose return path stores
//! no field.
//!
//! With field coordinate `x`, the charge on the plates [C] or the flux
//! linkage of the winding [Wb], and gap `g` [m], the stored energy is
//! `H = x^2 g / (2 s)`. The coupling scale `s` is `epsilon A` [F m] for a
//! capacitor and `mu N^2 A` [H m] for an inductor. Its gradient is the
//! terminal voltage (capacitor) or current (inductor) `x g / s`, and the
//! force `x^2 / (2 s)` with which the field resists the gap's opening: the
//! attraction between the plates or poles.
//!
//! Port `field` receives the rate of `x` and returns its conjugate: a
//! capacitor's port is `FlowIn` (current into the positive plate, voltage
//! back), an inductor's `EffortIn` (voltage across the winding, current
//! back). Port `gap` receives the rate at which the gap opens [m/s] and
//! returns the force [N] the field holds against that opening, positive as
//! the plates or poles attract. Each returned effort is the discrete
//! gradient of `H` between the step's two end points, the exact mean of the
//! gradient along the straight segment between them:
//!
//! ```text
//! e_x = x_m g_m / s + (dx / 12) (dg / s)
//! e_g = (x_m^2 + dx^2 / 12) / (2 s)
//! ```
//!
//! with `x_m`, `g_m` the midpoints and `dx`, `dg` the increments. So
//! `e_x dx + e_g dg` is `H(end) - H(start)` for every step, and joined to a
//! circuit and a mass the device is integrated by the same discrete-gradient
//! rule as one monolithic system.
//!
//! The gap is admitted on an open interval; a state or a step end outside
//! it is refused, never clamped. The coupling scale is a constant: a field
//! coefficient that changed with temperature would need the field's
//! entropy, which this element does not hold.

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::PortKind;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::{ElementError, inputs_finite, port, positive, roundoff};

/// The state of a [`GapField`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GapFieldState {
    /// Charge on the positive plate [C], or flux linkage of the winding [Wb].
    pub coordinate: f64,
    /// Gap width [m].
    pub gap_m: f64,
}

/// A field in a gap of variable width: see the [module documentation](self).
#[derive(Debug, Clone)]
pub struct GapField {
    /// `epsilon A` [F m] or `mu N^2 A` [H m].
    scale: f64,
    /// The open interval of admitted gaps [m].
    limits_m: [f64; 2],
    inductor: bool,
    ports: [PortSpec; 2],
}

impl GapField {
    /// The gap of a parallel-plate capacitor with permittivity times plate
    /// area `epsilon_area_f_m` [F m], its gap admitted strictly between the
    /// two `gap_limits_m` [m].
    ///
    /// # Errors
    /// A coupling scale that is not positive and finite or whose reciprocal
    /// is not, or limits that are not finite with `0 <= lower < upper`.
    pub fn capacitor(epsilon_area_f_m: f64, gap_limits_m: [f64; 2]) -> Result<Self, ElementError> {
        Self::build(epsilon_area_f_m, gap_limits_m, false)
    }

    /// The air gap of a winding of `N` turns with permeability times pole
    /// area times `N^2`, `mu_turns_area_h_m` [H m], its gap admitted strictly
    /// between the two `gap_limits_m` [m]. The coordinate is the flux
    /// linkage, not the magnetic flux.
    ///
    /// # Errors
    /// As for [`Self::capacitor`].
    pub fn inductor(mu_turns_area_h_m: f64, gap_limits_m: [f64; 2]) -> Result<Self, ElementError> {
        Self::build(mu_turns_area_h_m, gap_limits_m, true)
    }

    fn build(scale: f64, limits_m: [f64; 2], inductor: bool) -> Result<Self, ElementError> {
        let scale = positive("coupling_scale", scale)?;
        positive("1 / coupling_scale", 1.0 / scale)?;
        let [lower, upper] = limits_m;
        if !(lower.is_finite() && lower >= 0.0) {
            return Err(ElementError::Parameter {
                name: "lower gap_limits_m",
                value: lower,
            });
        }
        if !(upper.is_finite() && upper > lower) {
            return Err(ElementError::Parameter {
                name: "upper gap_limits_m",
                value: upper,
            });
        }
        let field = if inductor {
            PortRole::EffortIn
        } else {
            PortRole::FlowIn
        };
        Ok(Self {
            scale,
            limits_m,
            inductor,
            ports: [
                port("field", Some(PortKind::ElectricalVoltageCurrent), field),
                port(
                    "gap",
                    Some(PortKind::MechanicalForceVelocity),
                    PortRole::FlowIn,
                ),
            ],
        })
    }

    /// The coupling scale `s` [F m or H m].
    #[must_use]
    pub fn coupling_scale(&self) -> f64 {
        self.scale
    }

    /// The open interval of admitted gaps [m].
    #[must_use]
    pub fn gap_limits_m(&self) -> [f64; 2] {
        self.limits_m
    }

    /// Whether the field is a winding's (an inductor) rather than a
    /// capacitor's.
    #[must_use]
    pub fn is_inductor(&self) -> bool {
        self.inductor
    }

    /// The state, refused unless finite with its gap inside the limits.
    fn admitted(&self, state: &GapFieldState) -> Result<(), ElementError> {
        let [lower, upper] = self.limits_m;
        if !state.coordinate.is_finite() {
            return Err(ElementError::Law {
                reason: format!("the field coordinate {} is not finite", state.coordinate),
            });
        }
        if !(state.gap_m > lower && state.gap_m < upper) {
            return Err(ElementError::Law {
                reason: format!(
                    "the gap {} m is outside its admitted range ({lower}, {upper}) m",
                    state.gap_m
                ),
            });
        }
        Ok(())
    }

    /// The discrete gradient between `a` and `b`: field effort, then force.
    fn efforts(&self, a: &GapFieldState, b: &GapFieldState) -> [f64; 2] {
        let s = self.scale;
        let xm = f64::midpoint(a.coordinate, b.coordinate);
        let gm = f64::midpoint(a.gap_m, b.gap_m);
        let dx = b.coordinate - a.coordinate;
        let dg = b.gap_m - a.gap_m;
        [
            xm * (gm / s) + (dx / 12.0) * (dg / s),
            0.5 * (xm * (xm / s) + (dx / 12.0) * (dx / s)),
        ]
    }

    fn end(old: &GapFieldState, dt_s: f64, inputs: &[f64]) -> GapFieldState {
        GapFieldState {
            coordinate: old.coordinate + dt_s * inputs[0],
            gap_m: old.gap_m + dt_s * inputs[1],
        }
    }
}

impl Component for GapField {
    type State = GapFieldState;
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::GapField",
        state_family: "gap-field",
        schema_grammar: "coordinate f64; gap_m f64; finite, the gap inside the open admitted \
            range, with finite stored energy",
        codec_grammar: "two u64 words: the bits of coordinate then gap_m",
        codec_version: 1,
    });

    fn encode_state(&self, state: &GapFieldState, words: &mut Vec<u64>) -> bool {
        if !self.stored_energy_j(state).is_ok_and(f64::is_finite) {
            return false;
        }
        words.extend([state.coordinate.to_bits(), state.gap_m.to_bits()]);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<GapFieldState> {
        let &[coordinate, gap] = words else {
            return None;
        };
        let state = GapFieldState {
            coordinate: f64::from_bits(coordinate),
            gap_m: f64::from_bits(gap),
        };
        self.stored_energy_j(&state)
            .is_ok_and(f64::is_finite)
            .then_some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            u64::from(self.inductor),
            self.scale.to_bits(),
            self.limits_m[0].to_bits(),
            self.limits_m[1].to_bits(),
        ]);
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &GapFieldState) -> Result<f64, ElementError> {
        self.admitted(state)?;
        let x = state.coordinate;
        Ok(0.5 * (x * (x / self.scale)) * state.gap_m)
    }

    fn propose(
        &self,
        old: &GapFieldState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut GapFieldState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let end = Self::end(old, dt_s, inputs);
        let before_j = self.stored_energy_j(old)?;
        let after_j = self.stored_energy_j(&end)?;
        let [field, force] = self.efforts(old, &end);
        let field_work_j = field * (end.coordinate - old.coordinate);
        let gap_work_j = force * (end.gap_m - old.gap_m);
        let stored_change_j = field_work_j + gap_work_j;
        if !(before_j.is_finite() && after_j.is_finite() && stored_change_j.is_finite()) {
            return Err(ElementError::Law {
                reason: "the field's energy or efforts are not finite".to_owned(),
            });
        }
        *new = end;
        outputs[0] = field;
        outputs[1] = force;
        // The kernel books the two port works apart; one may nearly cancel
        // the other, so the rounding is set by each, not by their sum.
        let scale = before_j + after_j + field_work_j.abs() + gap_work_j.abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            roundoff(scale),
        )))
    }
}

impl PortTangent for GapField {
    fn jvp(
        &self,
        old: &GapFieldState,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let end = Self::end(old, dt_s, inputs);
        let s = self.scale;
        let xm = f64::midpoint(old.coordinate, end.coordinate);
        let gm = f64::midpoint(old.gap_m, end.gap_m);
        let dx = end.coordinate - old.coordinate;
        let dg = end.gap_m - old.gap_m;
        // The derivatives of the two efforts with respect to the end point:
        // d e_x / d x = (g_m / 2 + dg / 12) / s, d e_x / d g = d e_g / d x
        // = (x_m / 2 + dx / 12) / s, d e_g / d g = 0. The end moves by
        // `dt` times the inputs.
        let by_coordinate = (0.5 * gm + dg / 12.0) / s;
        let mixed = (0.5 * xm + dx / 12.0) / s;
        d_outputs[0] = dt_s * (by_coordinate * d_inputs[0] + mixed * d_inputs[1]);
        d_outputs[1] = dt_s * (mixed * d_inputs[0]);
        Ok(())
    }

    fn vjp(
        &self,
        old: &GapFieldState,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // The Jacobian is symmetric.
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}
