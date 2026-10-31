//! Source-backed, isothermal Newtonian journal shear in articulated dynamics.
//!
//! The exact concentric-cylinder laminar resistance is
//! `c = 4 pi mu L ri² ro² / (ro²-ri²)`, with joint torque `-c qdot`.
//! Equal and opposite cylinder torques give loss `c qdot²`; the joint's native
//! work-conjugate coordinate accounts for both moving links. See Fitzpatrick,
//! *Taylor-Couette Flow*, equations 10.38–10.40:
//! <https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node137.html>.
//!
//! A film has negligible mass/thermal storage and equilibrates to its explicit
//! receiving wall cell. No pressure-supported bearing load, eccentricity,
//! cavitation, turbulence, end effects or finite-time fluid spin-up is modeled.
//! Geometry may be prescribed or follow the same source-resolved free thermal
//! expansion as the two cylindrical material regions. Its slenderness and
//! quasisteady validity remain caller modeling obligations.
//! The creeping-flow limit is checked at every native mechanical RK stage.
//! Viscosity and geometry are frozen only within one explicitly capped thermal
//! substep, then refreshed from the same accepted material configuration.
//! The explicit rheological attachment instead freezes its whole original
//! T/shear source line; actual finite-gap shear and torque evolve at every stage.

use super::{ArticulatedModel, CellId, DepositionSimulation, Error, finite, positive};
use crate::thermal::{ExpandedArticulated, ExpandedRegion};
use fs_ga::Vec3;
use fs_matdb::ClaimSelection;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::state_point::MaterialStatePointError;
use fs_mbd::articulated::JointKind;
use fs_qty::Temperature;
use fs_thermal::enthalpy::DepositionError;
use fs_tribo::journal::JournalShearResponse;

/// Original-source generalized-Newtonian shear on this same thermal timeline.
pub mod rheology;
use rheology::{BoundRheology, JournalRheologyState, RheologicalJournalBearing, RheologyBinding};

/// Shared native concentric-cylinder dimensions, in SI units.
pub use fs_tribo::journal::JournalGeometry;

/// Explicit physical-domain and numerical-accuracy obligations for one film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalLimits {
    /// `rho |qdot| ri (ro-ri) / mu`; positive and at most one.
    pub maximum_gap_reynolds: f64,
    /// Maximum fractional endpoint viscosity change within ONE thermal step.
    /// Positive and at most one. Refine motion `maximum_step_s` on refusal.
    pub maximum_relative_viscosity_change: f64,
    /// Absolute tolerance for this film's additional RK loss coordinate [J].
    pub absolute_loss_tolerance_j: f64,
}

/// One physical film; multiple distinct films may act on the same joint.
#[derive(Clone, Debug)]
pub struct JournalBearing<'a> {
    /// Child-link index, whose relative joint MUST be revolute.
    pub link: usize,
    /// Active thermal cell receiving all loss and setting film temperature.
    pub heat_cell: CellId,
    pub geometry: JournalGeometry,
    pub limits: JournalLimits,
    /// Immutable card, fixed nonthermal context and originally selected claims.
    pub liquid: LiquidTemperatureLaw<'a>,
}

/// Bind reference cylinder dimensions to the original two material regions.
///
/// The inner cylinder belongs to the journal's child link and the outer cylinder
/// to its direct parent. Both are coaxial with the revolute joint and have the
/// reference `JournalGeometry.length_m` centered at that joint. Their expansion
/// fixed points must coincide with the joint origin in their own body frames.
/// This declared cylindrical geometry is not inferred from volume moments.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalJournalBinding {
    pub inner_region: String,
    pub outer_region: String,
    /// Maximum fractional change of the geometric resistance factor in one
    /// thermal substep, positive and at most one. Viscosity has its own bound.
    pub maximum_relative_geometry_factor_change: f64,
}

/// Existing film law with optional source-owned cylindrical expansion enabled.
#[derive(Clone, Debug)]
pub struct ExpandingJournalBearing<'a> {
    /// Geometry is at the original material reference state, not the current
    /// temperature. Repeated updates always use these same reference dimensions.
    pub bearing: JournalBearing<'a>,
    pub expansion: ThermalJournalBinding,
}

/// Material and physical observations at the actual accepted joint endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct JournalBearingState {
    pub link: usize,
    pub heat_cell: CellId,
    /// Actual admitted cylindrical radii and overlapping axial length [m].
    pub geometry: JournalGeometry,
    /// Original liquid at the finite wall temperature. In rheological mode
    /// this is the zero-shear reference; actual two-wall sources are below.
    pub liquid: ResolvedLiquidState,
    /// Torque/angular velocity at this endpoint; the zero-rate value is the
    /// source's zero-shear Newtonian tangent, not a fitted constant damper.
    pub resistance_nm_s: f64,
    pub resisting_torque_nm: f64,
    pub dissipated_power_w: f64,
    /// Uses the minimum viscosity over the complete source shear interval in
    /// rheological mode, giving a conservative creeping-flow bound.
    pub gap_reynolds: f64,
    pub rheology: Option<JournalRheologyState>,
}

#[derive(Clone)]
pub(super) struct BoundBearing<'a> {
    pub definition: JournalBearing<'a>,
    pub dof: usize,
    pub state: JournalBearingState,
    geometry_factor_m3: f64,
    expansion: Option<ThermalJournalBinding>,
    rheology: Option<BoundRheology<'a>>,
}

impl<'a> BoundBearing<'a> {
    pub fn new(
        definition: JournalBearing<'a>,
        model: &ArticulatedModel,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
    ) -> Result<Self, Error> {
        let geometry = definition.geometry;
        Self::with_geometry(definition, model, thermal, velocity, geometry, None, None,
            &mut || Ok(()))
    }

    pub fn new_expanding(
        definition: ExpandingJournalBearing<'a>,
        configuration: &ExpandedArticulated<'_>,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
    ) -> Result<Self, Error> {
        // Admit the reference as well: an initially heated configuration cannot
        // turn an overlapping or nonphysical authored journal into a valid one.
        geometry_factor(definition.bearing.geometry)?;
        let geometry = definition.expansion.resolve(
            definition.bearing.link,
            definition.bearing.geometry,
            configuration,
        )?;
        Self::with_geometry(
            definition.bearing,
            configuration.model(),
            thermal,
            velocity,
            geometry,
            Some(definition.expansion),
            None,
            &mut || Ok(()),
        )
    }

    pub fn new_rheological(
        definition: RheologicalJournalBearing<'a>,
        configuration: &ExpandedArticulated<'_>,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        geometry_factor(definition.bearing.geometry)?;
        let geometry = match &definition.expansion {
            Some(binding) => binding.resolve(definition.bearing.link,
                definition.bearing.geometry, configuration)?,
            None => definition.bearing.geometry,
        };
        let rheology = RheologyBinding { rate_axis: definition.shear_rate_axis,
            maximum_rate: definition.maximum_shear_rate_per_s, solve: definition.solve };
        Self::with_geometry(definition.bearing, configuration.model(), thermal,
            velocity, geometry, definition.expansion, Some(rheology), checkpoint)
    }

    #[allow(clippy::too_many_arguments)]
    fn with_geometry(
        definition: JournalBearing<'a>,
        model: &ArticulatedModel,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
        geometry: JournalGeometry,
        expansion: Option<ThermalJournalBinding>,
        rheology: Option<RheologyBinding<'a>>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        let limits = definition.limits;
        if !positive(limits.maximum_gap_reynolds)
            || limits.maximum_gap_reynolds > 1.0
            || !positive(limits.maximum_relative_viscosity_change)
            || limits.maximum_relative_viscosity_change > 1.0
            || !positive(limits.absolute_loss_tolerance_j)
        {
            return Err(Error::Invalid(
                "journal geometry and creeping-flow/loss bounds",
            ));
        }
        let link = model
            .links()
            .get(definition.link)
            .ok_or(Error::Invalid("journal child link does not exist"))?;
        if link.joint().kind() != JointKind::Revolute {
            return Err(Error::Invalid("journal shear requires a revolute joint"));
        }
        let dof = model
            .dof_index(definition.link)
            .ok_or(Error::Invalid("journal joint has no motion coordinate"))?;
        let factor = geometry_factor(geometry)?;
        let reference = definition.liquid.reference_state();
        if let Some(binding) = &rheology {
            binding.admit(&definition.liquid)?;
        } else if reference
            .material()
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(Error::Invalid(
                "Newtonian film cannot freeze strain, rate or field dependent source claims",
            ));
        }
        let temperature = thermal.temperature_k(definition.heat_cell)?;
        let liquid = definition
            .liquid
            .at_temperature(Temperature::new(temperature))?;
        support(&definition.liquid, reference, &liquid)?;
        let rheology = rheology.map(|binding| BoundRheology::compile(binding,
            &definition.liquid, &liquid, geometry, checkpoint)).transpose()?;
        let state = JournalBearingState {
            link: definition.link,
            heat_cell: definition.heat_cell,
            geometry,
            resistance_nm_s: factor * liquid.dynamic_viscosity_pa_s(),
            liquid,
            resisting_torque_nm: 0.0,
            dissipated_power_w: 0.0,
            gap_reynolds: 0.0,
            rheology: None,
        };
        if !positive(state.resistance_nm_s) {
            return Err(Error::Invalid(
                "journal material resistance is unrepresentable",
            ));
        }
        let mut result = Self {
            definition,
            dof,
            state,
            geometry_factor_m3: factor,
            expansion,
            rheology,
        };
        result.update_rate(velocity[dof])?;
        Ok(result)
    }

    /// This exact resistance supplies both native torque and heat derivative.
    pub fn response(&self, rate: f64) -> Result<(f64, f64, f64), Error> {
        let (torque, power, reynolds, _) = self.response_values(rate)?;
        Ok((torque, power, reynolds))
    }

    fn response_values(&self, rate: f64)
        -> Result<(f64, f64, f64, Option<JournalShearResponse>), Error> {
        let g = self.state.geometry;
        let kinematic_viscosity = self.rheology.as_ref().map_or(
            self.state.liquid.kinematic_viscosity_m2_s(),
            |rheology| rheology.native.minimum_viscosity_pa_s()
                / self.state.liquid.density_kg_m3());
        let reynolds = finite(
            rate.abs() * g.inner_radius_m * (g.outer_radius_m - g.inner_radius_m)
                / kinematic_viscosity,
        )?;
        if reynolds > self.definition.limits.maximum_gap_reynolds {
            return Err(Error::Budget("journal creeping-flow Reynolds limit"));
        }
        if let Some(rheology) = &self.rheology {
            let response = rheology.native.response(rate)?;
            return Ok((response.resisting_torque_nm, response.dissipated_power_w,
                reynolds, Some(response)));
        }
        let torque = finite(-self.state.resistance_nm_s * rate)?;
        let power = finite(-torque * rate)?;
        if rate != 0.0 && (torque == 0.0 || power == 0.0) {
            return Err(Error::Invalid("nonzero journal work underflow"));
        }
        Ok((torque, power, reynolds, None))
    }

    fn update_rate(&mut self, rate: f64) -> Result<(), Error> {
        let (torque, power, reynolds, response) = self.response_values(rate)?;
        if let (Some(rheology), Some(response)) = (&self.rheology, response) {
            self.state.rheology = Some(rheology.sample(&self.definition.liquid,
                &self.state.liquid, response)?);
            self.state.resistance_nm_s = if rate == 0.0 {
                self.geometry_factor_m3 * self.state.liquid.dynamic_viscosity_pa_s()
            } else { finite((torque / rate).abs())? };
        }
        self.state.resisting_torque_nm = torque;
        self.state.dissipated_power_w = power;
        self.state.gap_reynolds = reynolds;
        Ok(())
    }

    pub fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
        configuration: &ExpandedArticulated<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let temperature = thermal.temperature_k(self.definition.heat_cell)?;
        let liquid = self
            .definition
            .liquid
            .at_temperature(Temperature::new(temperature))?;
        support(&self.definition.liquid, &self.state.liquid, &liquid)?;
        let old = self.state.liquid.dynamic_viscosity_pa_s();
        let change = finite((liquid.dynamic_viscosity_pa_s() - old).abs() / old)?;
        if change > self.definition.limits.maximum_relative_viscosity_change {
            return Err(Error::Budget(
                "journal viscosity change requires a smaller thermal substep",
            ));
        }
        let geometry = match &self.expansion {
            Some(binding) => binding.resolve(
                self.definition.link,
                self.definition.geometry,
                configuration,
            )?,
            None => self.definition.geometry,
        };
        let factor = geometry_factor(geometry)?;
        if let Some(binding) = &self.expansion {
            let change =
                finite((factor - self.geometry_factor_m3).abs() / self.geometry_factor_m3)?;
            if change > binding.maximum_relative_geometry_factor_change {
                return Err(Error::Budget(
                    "journal clearance change requires a smaller thermal substep",
                ));
            }
        }
        if let Some(rheology) = &self.rheology {
            let rates = self.state.rheology.as_ref().expect("admitted rheological endpoint")
                .shear_rate_per_s;
            self.rheology = Some(rheology.at_temperature(&self.definition.liquid,
                &self.state.liquid, &liquid, geometry, rates,
                self.definition.limits.maximum_relative_viscosity_change, checkpoint)?);
        }
        let resistance = factor * liquid.dynamic_viscosity_pa_s();
        if !positive(resistance) {
            return Err(Error::Invalid(
                "journal material resistance is unrepresentable",
            ));
        }
        self.state.liquid = liquid;
        self.state.geometry = geometry;
        self.state.resistance_nm_s = resistance;
        self.geometry_factor_m3 = factor;
        self.update_rate(velocity[self.dof])
    }
}

fn geometry_factor(g: JournalGeometry) -> Result<f64, Error> {
    Ok(g.resistance_factor_m3()?)
}

impl ThermalJournalBinding {
    fn resolve(
        &self,
        link: usize,
        reference: JournalGeometry,
        configuration: &ExpandedArticulated<'_>,
    ) -> Result<JournalGeometry, Error> {
        if !positive(self.maximum_relative_geometry_factor_change)
            || self.maximum_relative_geometry_factor_change > 1.0
        {
            return Err(Error::Invalid(
                "journal geometry-change bound must lie in (0,1]",
            ));
        }
        let child = configuration
            .model()
            .links()
            .get(link)
            .ok_or(Error::Invalid("journal child link does not exist"))?;
        let parent = child.parent().ok_or(Error::Invalid(
            "expanding journal requires a material parent link",
        ))?;
        let inner = region(configuration, link, &self.inner_region)?;
        let outer = region(configuration, parent, &self.outer_region)?;
        if inner.fixed_point_m() != Vec3::default()
            || outer.fixed_point_m() != child.parent_from_child_zero().translation()
        {
            return Err(Error::Invalid(
                "expanding cylinders must share the revolute-joint expansion fixed point",
            ));
        }
        // The same integrated source strains have already configured link
        // geometry and inertia. Do not query another alpha or compound stretch.
        let inner_stretch = 1.0 + inner.expansion().free_linear_strain();
        let outer_stretch = 1.0 + outer.expansion().free_linear_strain();
        let geometry = JournalGeometry {
            inner_radius_m: reference.inner_radius_m * inner_stretch,
            outer_radius_m: reference.outer_radius_m * outer_stretch,
            length_m: reference.length_m * inner_stretch.min(outer_stretch),
        };
        geometry_factor(geometry)?;
        Ok(geometry)
    }
}

fn region<'s>(
    configuration: &'s ExpandedArticulated<'_>,
    link: usize,
    key: &str,
) -> Result<&'s ExpandedRegion, Error> {
    let body = configuration
        .bodies()
        .get(link)
        .ok_or(Error::Invalid("journal material link does not exist"))?;
    body.regions()
        .iter()
        .find(|region| body.reference().regions()[region.reference_index()].key() == key)
        .ok_or(Error::Invalid(
            "journal material region does not exist on its cylinder link",
        ))
}

// Ask the original source owner for continuous support of EVERY retained role;
// endpoint queries cannot jump over a hole or change the selected formulation.
fn support(
    law: &LiquidTemperatureLaw<'_>,
    from: &ResolvedLiquidState,
    to: &ResolvedLiquidState,
) -> Result<(), Error> {
    let from = from.material().query();
    let to = to.material().query();
    let axis = law.temperature_axis();
    let (low, high) = if from.axes()[axis] <= to.axes()[axis] {
        (&from, &to)
    } else {
        (&to, &from)
    };
    for property in law.reference_state().material().properties() {
        let selected = property.answer().receipt.selected;
        let claim = law
            .card()
            .claims()
            .claim(selected)
            .ok_or(Error::Invalid("journal liquid selected claim missing"))?;
        let path = law
            .card()
            .claims()
            .query_line_typed(
                &claim.key,
                axis,
                low,
                high,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.requirement().name().to_owned(),
                source: Box::new(source),
            })?;
        if path
            .knots
            .iter()
            .any(|(_, answer)| !positive(answer.evidence.value.value))
        {
            return Err(Error::Invalid(
                "journal liquid path contains an inadmissible coefficient",
            ));
        }
    }
    Ok(())
}
