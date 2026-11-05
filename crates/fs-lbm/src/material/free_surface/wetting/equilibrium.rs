//! Local thermal equilibrium uses the original wall and source transaction.
//!
//! The native thermal owner retains and checks the actual liquid/gas tension.
//! It is constant or admitted by a compatible thermal-capillarity owner.
//! The same local source temperature supplies work and tension; reversible
//! work is not reinterpreted as a measured angle,
//! contact-line friction, hysteresis or heat.

use super::*;
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, YoungDupreLaw, resolve_young_dupre_wetting,
};
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, SURFACE_TENSION_PROPERTY, resolve_liquid_wetting_state,
};

#[derive(Clone, Copy)]
pub(super) struct EquilibriumWetting {
    origin: EquilibriumAngleOrigin,
    surface_tension_n_m: Option<f64>,
}

impl<'m> WettingBinding<'m> {
    pub(in crate::material::free_surface) fn try_equilibrium(
        flow: &mut MaterialFreeSurface3<'_>,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        origin: EquilibriumAngleOrigin,
        selection: MaterialPropertySelection,
        temperature_axis: &str,
        evolving_surface_tension: bool,
        maximum_wall_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        if maximum_wall_cells == 0 || maximum_wall_cells > 262_144 {
            return Err(Error::Invalid("invalid equilibrium wall-cell budget"));
        }
        let liquid = flow.scaling.material();
        let wetting = match origin {
            EquilibriumAngleOrigin::SourcedAngle => {
                resolve_liquid_wetting_state(card, &liquid.material().query(), selection)?
            }
            EquilibriumAngleOrigin::YoungDupre => {
                resolve_young_dupre_wetting(card, liquid, selection)?
            }
        };
        wetting.require_liquid_state(liquid)?;
        // A scalar conditioned on age, load, speed, or another temperature
        // still needs an owner for that coordinate. Initialization must not
        // silently turn those applicability restrictions into constants.
        if wetting.resolved().properties().iter().any(|property| {
            !property.has_only_fixed_external_conditions(Some(temperature_axis), false)
        }) {
            return Err(Error::Invalid(
                "equilibrium wetting requires its owned temperature and fixed external conditions",
            ));
        }
        if origin == EquilibriumAngleOrigin::YoungDupre {
            let tension =
                liquid
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .ok_or(Error::Invalid(
                        "equilibrium adhesion requires original liquid/gas tension",
                    ))?;
            // The thermal owner has already pinned and admitted the SAME
            // source's temperature axis. Only its uniform or sidewall graph
            // mode can evolve it with the admitted contact-angle geometry.
            if !evolving_surface_tension
                && (!tension.source_axes().is_empty()
                    || tension.answer().receipt.decision
                        != fs_matdb::EvaluationDecision::ConstantWithinValidity)
            {
                return Err(Error::Invalid(
                    "thermal adhesion needs constant tension or its compatible thermal owner",
                ));
            }
        }
        let physics = WettingPhysics {
            card,
            source: wetting.resolved().clone(),
            controls: FreeSurfaceWettingControls {
                pinning_speed_m_s: 0.0,
                maximum_wall_cells,
            },
            hydrodynamic: None,
            equilibrium: Some(EquilibriumWetting {
                origin,
                surface_tension_n_m: liquid.surface_tension_n_m(),
            }),
        };
        Self::with_physics(flow, wall, physics, gate)
    }
}

impl WettingPhysics<'_> {
    pub(super) fn equilibrium_sample(
        &self,
        equilibrium: &EquilibriumWetting,
        geometric_angle_radians: f64,
        temperature: Option<(&str, WettingThermalPoint)>,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceWettingSample, Error> {
        checkpoint(gate)?;
        // The complete point and selected original claim are retained. No
        // fictitious speed axis is appended to an equilibrium source receipt.
        // Reference-temperature cells and dry ghosts reuse the admitted source
        // verbatim without rerunning claim resolution for each wall cell.
        let material = if let Some((axis, value)) = temperature
            && self.source.state_coordinate(axis) != Some(value.temperature_k)
        {
            let point = self.source.query().with(axis, value.temperature_k)?;
            self.source.at_query_point(self.card, &point)?
        } else {
            self.source.clone()
        };
        let name = match equilibrium.origin {
            EquilibriumAngleOrigin::SourcedAngle => EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
            EquilibriumAngleOrigin::YoungDupre => EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
        };
        let value = material
            .property(name)
            .ok_or(Error::Invalid(
                "equilibrium wetting lost its original source",
            ))?
            .value_si();
        let applied_angle_radians = match equilibrium.origin {
            EquilibriumAngleOrigin::SourcedAngle => value,
            EquilibriumAngleOrigin::YoungDupre => {
                let sigma = temperature
                    .map_or(equilibrium.surface_tension_n_m, |(_, state)| {
                        state.surface_tension_n_m
                    })
                    .ok_or(Error::Invalid(
                        "equilibrium adhesion requires liquid/gas tension",
                    ))?;
                YoungDupreLaw::new(value, sigma)
                    .map_err(|_| {
                        Error::Invalid(
                            "thermal adhesion leaves finite-angle partial-wetting equilibrium",
                        )
                    })?
                    .contact_angle_radians()
            }
        };
        if !(applied_angle_radians > 0.0 && applied_angle_radians < core::f64::consts::PI) {
            return Err(Error::Invalid(
                "thermal equilibrium angle must remain strictly in (0, pi)",
            ));
        }
        Ok(FreeSurfaceWettingSample {
            signed_speed_m_s: 0.0,
            regime: ContactLineRegime2::Equilibrium,
            geometric_angle_radians,
            applied_angle_radians,
            equilibrium: Some(equilibrium.origin),
            hydrodynamic: None,
            outer_gas: None,
            material,
        })
    }
}
