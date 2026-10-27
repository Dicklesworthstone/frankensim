//! Source-bound dry joint brakes on the existing articulated/enthalpy timeline.
//!
//! A held normal load acts through one prescribed parent-child sliding
//! coordinate: v = r*qdot for a revolute joint, v = qdot for a prismatic joint.
//! Native fs-tribo regularized traction supplies BOTH the generalized resisting
//! effort and its power. Integrated work is partitioned to the actual child and
//! parent thermal regions. Their distinct caloric temperatures refresh the same
//! ordered interface claims before subsequent motion.
//!
//! This is a regularized sliding model, not a static sticking/contact solver.
//! Preload, contact area, lever arm and heat fractions are explicit inputs;
//! normal compliance, wear, contact detection, directional texture response,
//! flash temperature and thermal contact conductance need their own owners.
use super::{Error, finite, nonnegative, positive};
use crate::thermal::heating::{HeatingError, ThermalRegion};
use fs_contact::interface_binding::{BoundDryInterfaceSystem, bind_dry_interface_system_card};
use fs_matdb::{ClaimSelection, QueryPoint};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::articulated::{ArticulatedModel, JointKind};
use fs_qty::{Dims, Pressure, Temperature};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::{DissipationStep, FrictionLaw, HeatPartition, InputAuthority, WorkLedger};

pub const STATIC_FRICTION_PROPERTY: &str = "static-friction-coefficient";
pub const KINETIC_FRICTION_PROPERTY: &str = "kinetic-friction-coefficient";

/// Prescribed slip Jacobian; the generalized effort is its transpose times force.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JointSlipGeometry {
    Revolute { effective_radius_m: f64 },
    Prismatic,
}

/// Independent integration and coefficient-freezing controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DryBrakeLimits {
    /// Native linear transition to zero traction below this slip speed [m/s].
    /// This is an explicit model regularization, not an inferred material law.
    pub regularization_m_s: f64,
    /// Maximum absolute change in either friction coefficient in one substep.
    /// Refusal asks the caller to refine the existing thermal split interval.
    pub maximum_coefficient_change: f64,
    pub absolute_loss_tolerance_j: f64,
}

/// One ordered child/parent dry interface with an explicitly held preload.
#[derive(Clone, Debug)]
pub struct DryJointBrake<'a> {
    pub link: usize,
    /// A is a region of the child link; B is a region of its direct parent.
    pub cells: [CellId; 2],
    pub geometry: JointSlipGeometry,
    pub normal_force_n: f64,
    pub nominal_area_m2: f64,
    /// Existing pressure coordinate replaced by normal_force / nominal_area.
    pub pressure_axis: &'a str,
    pub point: ThermalInterfacePoint<'a>,
    pub selection: MaterialPropertySelection,
    pub authority: InputAuthority,
    /// Explicit fractions into these TWO finite bodies; other must be zero.
    /// This does not infer a universal partition from bulk material names.
    pub partition: HeatPartition,
    pub limits: DryBrakeLimits,
}

/// Accepted physical state and the original source answers, in input order.
#[derive(Clone, Debug, PartialEq)]
pub struct DryBrakeState {
    pub link: usize,
    pub cells: [CellId; 2],
    pub interface: BoundDryInterfaceSystem,
    pub resolved: ResolvedInterfaceStatePoint,
    pub phases: [EquilibriumPhaseState; 2],
    pub law: FrictionLaw,
    pub pressure_pa: f64,
    pub slip_velocity_m_s: f64,
    /// N m for a revolute joint, N for a prismatic joint.
    pub resisting_effort: f64,
    pub dissipated_power_w: f64,
    pub partition: HeatPartition,
    /// Cumulative accepted loss and its two named thermal destinations.
    pub work: WorkLedger,
}

#[derive(Clone)]
pub(super) struct BoundBrake<'a> {
    pub definition: DryJointBrake<'a>,
    pub dof: usize,
    pub state: DryBrakeState,
    regions: [ThermalRegion<'a>; 2],
    slip_jacobian: f64,
    response_axis: Option<String>,
}

impl<'a> BoundBrake<'a> {
    pub fn new(
        definition: DryJointBrake<'a>,
        model: &ArticulatedModel,
        regions: &[Vec<ThermalRegion<'a>>],
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        let limits = definition.limits;
        if ![
            definition.normal_force_n,
            definition.nominal_area_m2,
            limits.regularization_m_s,
            limits.maximum_coefficient_change,
            limits.absolute_loss_tolerance_j,
        ]
        .into_iter()
        .all(positive)
            || definition.cells[0] == definition.cells[1]
            || definition.partition.other() != 0.0
        {
            return Err(Error::Invalid(
                "dry brake needs positive preload/area/bounds and two distinct receiving bodies",
            ));
        }
        let link = model
            .links()
            .get(definition.link)
            .ok_or(Error::Invalid("dry brake child link does not exist"))?;
        let parent = link
            .parent()
            .ok_or(Error::Invalid("dry brake requires a parent-child joint"))?;
        let slip_jacobian = match (link.joint().kind(), definition.geometry) {
            (JointKind::Revolute, JointSlipGeometry::Revolute { effective_radius_m })
                if positive(effective_radius_m) =>
            {
                effective_radius_m
            }
            (JointKind::Prismatic, JointSlipGeometry::Prismatic) => 1.0,
            _ => {
                return Err(Error::Invalid(
                    "dry brake slip geometry must match the moving joint",
                ));
            }
        };
        let dof = model
            .dof_index(definition.link)
            .ok_or(Error::Invalid("dry brake joint has no motion coordinate"))?;
        let mut selected = Vec::with_capacity(2);
        for (side, index) in [definition.link, parent].into_iter().enumerate() {
            let region = regions[index]
                .iter()
                .find(|r| r.cell_id == definition.cells[side])
                .ok_or(Error::Invalid(
                    "dry interface A/B cells must belong to this child/parent",
                ))?;
            if region.card.content_hash() != definition.point.materials[side].content_hash() {
                return Err(Error::Invalid(
                    "dry interface must retain each actual region material card",
                ));
            }
            selected.push(region.clone());
        }
        let regions = [selected[0].clone(), selected[1].clone()];
        let axis = definition.pressure_axis;
        if axis.trim().is_empty()
            || definition.point.temperature_axes.contains(&axis)
            || definition.point.liquid_fraction_axes.contains(&Some(axis))
            || definition.point.mean_temperature_axis == Some(axis)
            || !definition.point.context.axes().contains_key(axis)
            || definition
                .point
                .context
                .axis_quantities()
                .get(axis)
                .is_none_or(|q| q.dims() != Pressure::DIMS)
        {
            return Err(Error::Invalid(
                "preload pressure needs a distinct existing typed pressure coordinate",
            ));
        }
        let pressure = finite(definition.normal_force_n / definition.nominal_area_m2)?;
        if !positive(pressure) {
            return Err(Error::Invalid("dry nominal pressure underflow"));
        }
        let interface =
            bind_dry_interface_system_card(definition.point.interface, definition.authority)?;
        let phases = phases(&regions, thermal)?;
        let point = query(&definition, phases, pressure)?;
        let requirements = [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY]
            .into_iter()
            .map(|name| {
                ScalarPropertyRequirement::try_new(
                    name,
                    Dims::NONE,
                    ScalarAdmissibility::NonNegative,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let resolved = resolve_interface_state_point(
            definition.point.interface,
            &point,
            &requirements,
            definition.selection.clone(),
        )?;
        let response_axis = admit_source(&definition, &resolved)?;
        let law = law(&resolved)?;
        let state = DryBrakeState {
            link: definition.link,
            cells: definition.cells,
            interface,
            resolved,
            phases,
            law,
            pressure_pa: pressure,
            slip_velocity_m_s: 0.0,
            resisting_effort: 0.0,
            dissipated_power_w: 0.0,
            partition: definition.partition,
            work: WorkLedger::default(),
        };
        let mut result = Self {
            definition,
            dof,
            state,
            regions,
            slip_jacobian,
            response_axis,
        };
        result.support(&result.state.resolved, checkpoint)?;
        result.update_rate(velocity[dof])?;
        checkpoint()?;
        Ok(result)
    }

    /// Identical traction supplies the ABA effort and the RK heat derivative.
    pub fn response(&self, rate: f64) -> Result<(f64, f64), Error> {
        let slip = finite(self.slip_jacobian * rate)?;
        let force = self.state.law.regularized_traction_1d(
            slip,
            self.definition.normal_force_n,
            self.definition.limits.regularization_m_s,
        )?;
        let effort = finite(self.slip_jacobian * force)?;
        let power = finite(-effort * rate)?;
        if !nonnegative(power) {
            return Err(Error::Invalid(
                "dry brake traction must dissipate mechanical work",
            ));
        }
        if rate != 0.0
            && self.state.law.zero_speed_kinetic_coefficient()? > 0.0
            && (slip == 0.0 || effort == 0.0 || power == 0.0)
        {
            return Err(Error::Invalid("nonzero dry brake response underflow"));
        }
        Ok((effort, power))
    }

    pub fn partition_loss(&self, loss_j: f64, dt_s: f64) -> Result<DissipationStep, Error> {
        if !nonnegative(loss_j) {
            return Err(Error::Budget(
                "dry brake integrated loss must be nonnegative; refine RK tolerances",
            ));
        }
        Ok(DissipationStep::from_power(
            finite(loss_j / dt_s)?,
            dt_s,
            self.definition.partition,
        )?)
    }

    pub fn accept(
        &mut self,
        loss: DissipationStep,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        checkpoint()?;
        let phases = phases(&self.regions, thermal)?;
        let point = query(&self.definition, phases, self.state.pressure_pa)?;
        let resolved = self
            .state
            .resolved
            .at_query_point(self.definition.point.interface, &point)?;
        self.support(&resolved, checkpoint)?;
        let next_law = law(&resolved)?;
        for name in [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY] {
            let old = self
                .state
                .resolved
                .property(name)
                .expect("resolved requirement")
                .value_si();
            let next = resolved
                .property(name)
                .expect("resolved requirement")
                .value_si();
            if (next - old).abs() > self.definition.limits.maximum_coefficient_change {
                return Err(Error::Budget(
                    "dry friction coefficient changed too far; refine thermal substeps",
                ));
            }
        }
        self.state.work.record(loss)?;
        self.state.resolved = resolved;
        self.state.phases = phases;
        self.state.law = next_law;
        self.update_rate(velocity[self.dof])?;
        checkpoint()?;
        Ok(())
    }

    fn update_rate(&mut self, rate: f64) -> Result<(), Error> {
        let (effort, power) = self.response(rate)?;
        self.state.slip_velocity_m_s = finite(self.slip_jacobian * rate)?;
        self.state.resisting_effort = effort;
        self.state.dissipated_power_w = power;
        Ok(())
    }

    // Both coefficients may be scalar or depend on ONE shared evolving
    // temperature coordinate (A, B or explicitly declared mean). A source
    // surface may additionally depend on held pressure/composition. This lets
    // the native exact source-line owner prove continuous support, including
    // crossed table knots and holes; no sampled 2-temperature diagonal is
    // passed off as a continuous-domain check.
    fn support(
        &self,
        next: &ResolvedInterfaceStatePoint,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let Some(axis) = &self.response_axis else {
            return Ok(());
        };
        let old = self.state.resolved.query();
        let a = old.axes()[axis];
        let b = next.query().axes()[axis];
        let lower = old.clone().with(axis, a.min(b)).map_err(|e| {
            Error::ThermalMaterial(
                fs_material::thermal::source::ThermalMaterialError::Coordinate(e.to_string()),
            )
        })?;
        let upper = old.clone().with(axis, a.max(b)).map_err(|e| {
            Error::ThermalMaterial(
                fs_material::thermal::source::ThermalMaterialError::Coordinate(e.to_string()),
            )
        })?;
        let claims = self.definition.point.interface.claims();
        let mut coordinates = Vec::new();
        for property in self.state.resolved.properties() {
            checkpoint()?;
            let selected = property.answer().receipt.selected;
            let claim = claims
                .claim(selected)
                .expect("retained original interface claim");
            let line = claims
                .query_line_typed(
                    &claim.key,
                    axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: claim.key.name().to_owned(),
                    source,
                })?;
            coordinates.extend(line.knots.into_iter().map(|(coordinate, _)| coordinate));
        }
        coordinates.sort_by(f64::total_cmp);
        coordinates.dedup();
        for temperature in coordinates {
            checkpoint()?;
            let point = old.clone().with(axis, temperature).map_err(|e| {
                Error::ThermalMaterial(
                    fs_material::thermal::source::ThermalMaterialError::Coordinate(e.to_string()),
                )
            })?;
            // Native Coulomb admission must hold at the union of BOTH source
            // breakpoint sets. Endpoint-only checks can miss mu_k > mu_s.
            law(&self
                .state
                .resolved
                .at_query_point(self.definition.point.interface, &point)?)?;
        }
        Ok(())
    }
}

fn phases(
    regions: &[ThermalRegion<'_>; 2],
    thermal: &DepositionSimulation<'_>,
) -> Result<[EquilibriumPhaseState; 2], Error> {
    let mut states = Vec::with_capacity(2);
    for region in regions {
        let cell = thermal
            .cells()
            .get(&region.cell_id)
            .ok_or(Error::Invalid("dry brake thermal cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, region.caloric.curve()) {
            return Err(Error::Invalid(
                "dry brake must retain each original caloric object",
            ));
        }
        let phase = region
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(HeatingError::from)?;
        if !phase.is_fully_solid() {
            return Err(Error::Invalid(
                "dry brake cannot use a melted contact surface",
            ));
        }
        states.push(phase);
    }
    Ok([states[0], states[1]])
}

fn query(
    definition: &DryJointBrake<'_>,
    phases: [EquilibriumPhaseState; 2],
    pressure: f64,
) -> Result<QueryPoint, Error> {
    Ok(definition
        .point
        .at(phases)?
        .with(definition.pressure_axis, pressure)
        .map_err(|e| {
            fs_material::thermal::source::ThermalMaterialError::Coordinate(e.to_string())
        })?)
}

fn law(resolved: &ResolvedInterfaceStatePoint) -> Result<FrictionLaw, Error> {
    let law = FrictionLaw::Coulomb {
        static_mu: resolved
            .property(STATIC_FRICTION_PROPERTY)
            .expect("resolved requirement")
            .value_si(),
        kinetic_mu: resolved
            .property(KINETIC_FRICTION_PROPERTY)
            .expect("resolved requirement")
            .value_si(),
    };
    law.validate()?;
    Ok(law)
}

fn admit_source(
    definition: &DryJointBrake<'_>,
    resolved: &ResolvedInterfaceStatePoint,
) -> Result<Option<String>, Error> {
    let mut response_axis: Option<String> = None;
    for property in resolved.properties() {
        if property.source_axes().is_empty()
            && property.answer().receipt.decision
                != fs_matdb::EvaluationDecision::ConstantWithinValidity
        {
            return Err(Error::Invalid(
                "dry friction scalar must explicitly hold throughout its validity domain",
            ));
        }
        if property
            .angular_frequency_band_rad_s()?
            .is_some_and(|(minimum, _)| minimum > 0.0)
        {
            return Err(Error::Invalid(
                "dry friction source must admit static response",
            ));
        }
        let claim = definition
            .point
            .interface
            .claims()
            .claim(property.answer().receipt.selected)
            .expect("resolved original interface claim");
        let allowed = |axis: &str, dims: Dims| {
            definition.point.is_response_axis(axis, dims)
                || axis == definition.pressure_axis && dims == Pressure::DIMS
        };
        for (axis, dims) in property.source_axes() {
            if !allowed(axis, *dims) {
                return Err(Error::Invalid(
                    "dry brake coefficients must be independent of unowned rate/field/history coordinates",
                ));
            }
            if definition.point.is_thermal_axis(axis, *dims) {
                if *dims != Temperature::DIMS {
                    return Err(Error::Invalid(
                        "dry brake source cannot evolve a non-solid phase fraction",
                    ));
                }
                if response_axis.as_deref().is_some_and(|old| old != axis) {
                    return Err(Error::Invalid(
                        "dry brake coefficients support one shared evolving temperature coordinate",
                    ));
                }
                response_axis = Some(axis.clone());
            }
        }
        // Applicability-only restrictions need the same ownership as response
        // axes: a scalar claim with a frozen speed/strain bound is not admitted.
        for axis in claim.validity.bounds().keys() {
            let dims = if definition.point.temperature_axes.contains(&axis.as_str())
                || definition.point.mean_temperature_axis == Some(axis.as_str())
            {
                Some(Temperature::DIMS)
            } else if definition
                .point
                .liquid_fraction_axes
                .contains(&Some(axis.as_str()))
            {
                Some(Dims::NONE)
            } else {
                definition
                    .point
                    .context
                    .axis_quantities()
                    .get(axis)
                    .map(|q| q.dims())
            };
            if dims.is_none_or(|dims| !allowed(axis, dims)) {
                return Err(Error::Invalid(
                    "dry friction applicability contains an unowned physical coordinate",
                ));
            }
        }
    }
    Ok(response_axis)
}
