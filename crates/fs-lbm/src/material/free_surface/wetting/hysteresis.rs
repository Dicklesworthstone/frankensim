//! Independently sourced static thresholds around the original Young angle.
//!
//! The equilibrium angle/work remains a reversible source. Advancing/receding
//! claims supply separate rate-independent resistance; they never replace that
//! potential or fabricate its hysteresis width. Native first-cell geometry and
//! accepted contact motion select the branch. This explicit angle closure does
//! not constrain the VOF contact position exactly or close a surface-energy
//! budget. Per-length resistance/loss is diagnostic, never caloric heat.

use super::*;

/// Equilibrium adhesion and separately selected static pinning thresholds.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumHysteresisSample {
    pub equilibrium_angle_radians: f64,
    pub receding_angle_radians: f64,
    pub advancing_angle_radians: f64,
    /// The same local liquid/gas tension used by native capillary pressure.
    pub surface_tension_n_m: f64,
    /// A nonzero tangential fill gradient resolves a first-cell contact line.
    pub contact_line_resolved: bool,
    /// Speed consumed by the rate model [m/s]. It is zero inside the numerical
    /// deadband, including geometric depinning, while the outer sample retains
    /// the signed observed motion. A zero-width bracket recovers smooth wetting.
    pub model_speed_m_s: f64,
    /// Signed static resistance sigma*(cos(theta_Y)-cos(theta_branch)) [N/m].
    /// A stuck branch carries the bounded reaction at its geometric angle.
    pub pinning_force_per_length_n_m: f64,
    /// Nonnegative resistance times modeled speed [W/m]; not an energy deposit.
    pub pinning_dissipation_per_length_w_m: f64,
    /// Original independently selected advancing/receding claims and receipts.
    /// The outer sample's material remains the original equilibrium claim.
    pub thresholds: ResolvedInterfaceStatePoint,
}

#[derive(Clone)]
pub(super) struct EquilibriumHysteresisBinding {
    source: ResolvedInterfaceStatePoint,
    pinning_speed_m_s: f64,
    reference_surface_tension_n_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct HysteresisBranch {
    regime: ContactLineRegime2,
    angle_radians: f64,
    model_speed_m_s: f64,
}

// The observed direction governs resolved motion. At rest/deadband, the actual
// geometric angle supplies the static reaction or selects depinning. In
// particular a tiny opposite-signed speed must not reverse the depinning rate.
fn select_hysteresis(
    speed: f64,
    deadband: f64,
    geometric: f64,
    contact_line_resolved: bool,
    young: f64,
    receding: f64,
    advancing: f64,
) -> HysteresisBranch {
    if !contact_line_resolved {
        return HysteresisBranch {
            regime: ContactLineRegime2::Equilibrium,
            angle_radians: young,
            model_speed_m_s: 0.0,
        };
    }
    let (regime, angle_radians, model_speed_m_s) = if receding == advancing {
        // No static resistance remains. Preserve the original smooth-wall
        // response even when the caller declared a finite numerical deadband.
        (
            if speed > 0.0 {
                ContactLineRegime2::Advancing
            } else if speed < 0.0 {
                ContactLineRegime2::Receding
            } else {
                ContactLineRegime2::Equilibrium
            },
            young,
            speed,
        )
    } else if speed > deadband {
        (ContactLineRegime2::Advancing, advancing, speed)
    } else if speed < -deadband {
        (ContactLineRegime2::Receding, receding, speed)
    } else if geometric > advancing {
        (ContactLineRegime2::Advancing, advancing, 0.0)
    } else if geometric < receding {
        (ContactLineRegime2::Receding, receding, 0.0)
    } else {
        (ContactLineRegime2::Pinned, geometric, 0.0)
    };
    HysteresisBranch {
        regime,
        angle_radians,
        model_speed_m_s,
    }
}

// The sine product avoids subtracting nearly equal cosines for narrow sourced
// hysteresis. Positive resistance opposes advancing; negative opposes receding.
fn pinning_force(young: f64, branch: f64, sigma: f64) -> f64 {
    (sigma * fs_math::det::sin(0.5 * (branch - young)))
        * (2.0 * fs_math::det::sin(0.5 * (young + branch)))
}

impl WettingPhysics<'_> {
    pub(in crate::material::free_surface) fn with_fixed_equilibrium_hysteresis(
        self,
        selection: MaterialPropertySelection,
        pinning_speed_m_s: f64,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        self.with_equilibrium_hysteresis_at_context(
            selection, pinning_speed_m_s, None, scaling, gate,
        )
    }

    pub(super) fn with_equilibrium_hysteresis(
        self,
        selection: MaterialPropertySelection,
        pinning_speed_m_s: f64,
        temperature_axis: &str,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        self.with_equilibrium_hysteresis_at_context(
            selection, pinning_speed_m_s, Some(temperature_axis), scaling, gate,
        )
    }

    fn with_equilibrium_hysteresis_at_context(
        mut self,
        selection: MaterialPropertySelection,
        pinning_speed_m_s: f64,
        temperature_axis: Option<&str>,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        let equilibrium = self.equilibrium.as_ref().ok_or(Error::Invalid(
            "independent pinning thresholds require an original equilibrium angle or work",
        ))?;
        if self.equilibrium_hysteresis.is_some() {
            return Err(Error::Invalid(
                "equilibrium wetting already owns pinning thresholds",
            ));
        }
        if !pinning_speed_m_s.is_finite() || pinning_speed_m_s < 0.0 {
            return Err(Error::Invalid(
                "pinning speed deadband must be finite and nonnegative",
            ));
        }
        let requirements = [
            RECEDING_CONTACT_ANGLE_PROPERTY,
            ADVANCING_CONTACT_ANGLE_PROPERTY,
        ]
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: 0.0,
                    upper: core::f64::consts::PI,
                },
            )
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
        let source = resolve_interface_state_point(
            self.card,
            &self.source.query(),
            &requirements,
            selection,
        )?;
        if source.properties().iter().any(|property| {
            !property.has_only_fixed_external_conditions(temperature_axis, false)
        }) {
            return Err(Error::Invalid(
                "pinning thresholds require static claims with owned temperature and fixed external conditions",
            ));
        }
        let (receding, advancing) = bracket(&source)?;
        let young = self
            .equilibrium_sample(equilibrium, 0.0, None, gate)?
            .applied_angle_radians;
        if !(receding <= young && young <= advancing) {
            return Err(Error::Invalid(
                "pinning thresholds must bracket the original equilibrium angle",
            ));
        }
        let sigma = scaling
            .material()
            .surface_tension_n_m()
            .filter(|value| value.is_finite() && *value > 0.0)
            .ok_or(Error::Invalid(
                "pinning requires the same positive tension as native pressure",
            ))?;
        for angle in [receding, advancing] {
            if let Some(binding) = &self.hydrodynamic {
                binding.law.admit_angle(angle)?;
            }
            if let Some(binding) = &self.molecular_kinetic {
                binding.law.admit_angle(angle)?;
            }
        }
        self.equilibrium_hysteresis = Some(EquilibriumHysteresisBinding {
            source,
            pinning_speed_m_s,
            reference_surface_tension_n_m: sigma,
        });
        checkpoint(gate)?;
        Ok(self)
    }
}

impl EquilibriumHysteresisBinding {
    pub(super) fn angles(&self) -> Result<(f64, f64), Error> {
        bracket(&self.source)
    }

    pub(super) fn sample(
        &self,
        card: &InterfaceSystemCard,
        young: f64,
        geometric: f64,
        contact_line_resolved: bool,
        speed: f64,
        temperature: Option<(&str, WettingThermalPoint)>,
        gate: &CancelGate,
    ) -> Result<(ContactLineRegime2, f64, EquilibriumHysteresisSample), Error> {
        checkpoint(gate)?;
        if !speed.is_finite()
            || !geometric.is_finite()
            || !(0.0..=core::f64::consts::PI).contains(&geometric)
        {
            return Err(Error::Invalid(
                "nonfinite or invalid hysteresis contact geometry/motion",
            ));
        }
        let thresholds = if let Some((axis, state)) = temperature
            && self.source.state_coordinate(axis) != Some(state.temperature_k)
        {
            self.source
                .at_query_point(card, &self.source.query().with(axis, state.temperature_k)?)?
        } else {
            self.source.clone()
        };
        let (receding, advancing) = bracket(&thresholds)?;
        if !(receding <= young && young <= advancing) {
            return Err(Error::Invalid(
                "thermal pinning thresholds no longer bracket equilibrium adhesion",
            ));
        }
        let sigma = temperature
            .map_or(Some(self.reference_surface_tension_n_m), |(_, state)| {
                state.surface_tension_n_m
            })
            .filter(|value| value.is_finite() && *value > 0.0)
            .ok_or(Error::Invalid(
                "pinning requires the same positive tension as native pressure",
            ))?;
        let branch = select_hysteresis(
            speed,
            self.pinning_speed_m_s,
            geometric,
            contact_line_resolved,
            young,
            receding,
            advancing,
        );
        let force = pinning_force(young, branch.angle_radians, sigma);
        let loss = force * branch.model_speed_m_s;
        if !force.is_finite() || !loss.is_finite() || loss < 0.0 {
            return Err(Error::Invalid(
                "pinning resistance or dissipation is unrepresentable",
            ));
        }
        Ok((
            branch.regime,
            branch.angle_radians,
            EquilibriumHysteresisSample {
                equilibrium_angle_radians: young,
                receding_angle_radians: receding,
                advancing_angle_radians: advancing,
                surface_tension_n_m: sigma,
                contact_line_resolved,
                model_speed_m_s: branch.model_speed_m_s,
                pinning_force_per_length_n_m: force,
                pinning_dissipation_per_length_w_m: loss.max(0.0),
                thresholds,
            },
        ))
    }
}

impl WettingBinding<'_> {
    pub(super) fn with_fixed_equilibrium_hysteresis(
        mut self,
        selection: MaterialPropertySelection,
        pinning_speed_m_s: f64,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        for physics in self.physics_mut() {
            *physics = physics.clone().with_fixed_equilibrium_hysteresis(
                selection.clone(), pinning_speed_m_s, scaling, gate,
            )?;
        }
        let mut samples = BTreeMap::new();
        for (&position, current) in &self.geometry {
            checkpoint(gate)?;
            let speed = self.samples[&position].signed_speed_m_s;
            samples.insert(position, self.physics_at(position).sample(
                current.angle_radians, current.tangent_gradient > 0.0, speed, None, gate,
            )?);
        }
        self.samples = samples;
        Ok(self)
    }

    pub(in crate::material::free_surface) fn with_equilibrium_hysteresis(
        mut self,
        selection: MaterialPropertySelection,
        pinning_speed_m_s: f64,
        temperature_axis: &str,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        for physics in self.physics_mut() {
            *physics = physics.clone().with_equilibrium_hysteresis(
                selection.clone(),
                pinning_speed_m_s,
                temperature_axis,
                scaling,
                gate,
            )?;
        }
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_sticking_and_rest_depinning_follow_the_measured_geometric_angle() {
        for speed in [-0.01, 0.0, 0.01] {
            let stuck = select_hysteresis(speed, 0.01, 1.5, true, 1.4, 1.2, 1.7);
            assert_eq!(
                stuck,
                HysteresisBranch {
                    regime: ContactLineRegime2::Pinned,
                    angle_radians: 1.5,
                    model_speed_m_s: 0.0,
                }
            );
            let advancing = select_hysteresis(speed, 0.01, 1.8, true, 1.4, 1.2, 1.7);
            assert_eq!(advancing.regime, ContactLineRegime2::Advancing);
            assert_eq!(advancing.angle_radians, 1.7);
            assert_eq!(advancing.model_speed_m_s, 0.0);
            let receding = select_hysteresis(speed, 0.01, 1.1, true, 1.4, 1.2, 1.7);
            assert_eq!(receding.regime, ContactLineRegime2::Receding);
            assert_eq!(receding.angle_radians, 1.2);
            assert_eq!(receding.model_speed_m_s, 0.0);
        }
        for angle in [1.2, 1.7] {
            assert_eq!(
                select_hysteresis(0.0, 0.01, angle, true, 1.4, 1.2, 1.7).regime,
                ContactLineRegime2::Pinned
            );
        }
    }

    #[test]
    fn g1_observed_motion_and_reversal_select_positive_static_loss() {
        for speed in [-0.02, 0.02] {
            let branch = select_hysteresis(speed, 0.01, 1.5, true, 1.4, 1.2, 1.7);
            assert_eq!(branch.model_speed_m_s, speed);
            assert_eq!(branch.angle_radians, if speed > 0.0 { 1.7 } else { 1.2 });
            let force = pinning_force(1.4, branch.angle_radians, 0.004);
            let expected = 0.004 * (1.4_f64.cos() - branch.angle_radians.cos());
            assert!((force - expected).abs() < 1e-17);
            assert!(force * speed > 0.0);
        }
        let narrow = pinning_force(1.4, 1.4 + 1e-12, 0.004);
        assert!(narrow > 0.0);
        assert!((narrow / (0.004 * 1.4_f64.sin() * ((1.4 + 1e-12) - 1.4)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn g3_no_contact_and_collapsed_thresholds_do_not_invent_pinning() {
        let dry = select_hysteresis(0.0, 0.01, 2.0, false, 1.4, 1.2, 1.7);
        assert_eq!(
            dry,
            HysteresisBranch {
                regime: ContactLineRegime2::Equilibrium,
                angle_radians: 1.4,
                model_speed_m_s: 0.0,
            }
        );
        for speed in [-0.001, 0.0, 0.001] {
            let smooth = select_hysteresis(speed, 0.01, 2.0, true, 1.4, 1.4, 1.4);
            assert_eq!(smooth.angle_radians, 1.4);
            assert_eq!(smooth.model_speed_m_s, speed);
            assert_eq!(pinning_force(1.4, smooth.angle_radians, 0.004), 0.0);
        }
    }
}
