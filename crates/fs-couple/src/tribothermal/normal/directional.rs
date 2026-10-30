//! Evolving normal contact and sourced directional sliding at one thermal endpoint.
//!
//! Prescribed normal indentation follows the original Hertz/Hunt--Crossley
//! segment. Its endpoint load and actual circular/elliptic area supply friction
//! and pressure. Tangent displacement is prescribed in one fixed material frame.
//! Friction work heats the original finite bodies inside the SAME enthalpy solve,
//! immediately changing both source friction and reciprocal normal stiffness.
//! Normal reversible heat remains separate from physical frictional dissipation.
//!
//! Tangential load/pressure and thermal composition are endpoint approximations:
//! refine time when load, temperature or source coefficients change. This owner
//! does not locate contact/separation events or solve body inertia, normal/tangent
//! velocity constraints, rotation, moving footprints or experimental calibration.

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    EllipticFrictionMode, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxisKind, SourceThermalCoulombError,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialError;
use fs_math::det;
use fs_qty::Temperature;
use fs_tribo::{
    DissipationStep, TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput,
    TriboError, WorkLedger,
};
use std::sync::Arc;

/// Explicit interface coordinates, fixed tangent frame and heat destinations.
#[derive(Clone, Copy, Debug)]
pub struct NormalDirectionalInterface<'a> {
    pub point: ThermalInterfacePoint<'a>,
    pub tangent_frame_id: &'a str,
    pub partition: HeatPartition,
    /// Required when a nonzero fraction leaves the two finite caloric bodies.
    pub outside_destination: Option<&'a str>,
}

#[derive(Debug)]
pub enum DirectionalEndpointError<E = Infallible> {
    Invalid(&'static str),
    Friction(SourceThermalCoulombError),
    Material(ThermalMaterialError),
    Phase(PhaseStateError),
    Tribology(TriboError),
    Thermal(DepositionError),
    Source(E),
}
pub type NormalDirectionalError<E = Infallible> = NormalThermalError<DirectionalEndpointError<E>>;

/// Original native friction and caller transport at the accepted shared endpoint.
pub struct NormalDirectionalEndpoint<S> {
    pub friction: TangentialStep<SourceThermalEllipticState>,
    pub partition: DissipationStep,
    pub contact_area_m2: f64,
    pub pressure_pa: f64,
    pub tangential_work_j: f64,
    pub external_heat_j: f64,
    pub source: S,
}

pub struct NormalDirectionalReport<S> {
    /// Its INNER external heat includes friction heat. The OUTER external heat
    /// below excludes that internal exchange and closes the combined account.
    pub normal: NormalThermalReport<NormalDirectionalEndpoint<S>>,
    pub tangential_work_j: f64,
    pub total_mechanical_work_j: f64,
    pub external_heat_j: f64,
    /// Delta(H + U_normal) + exported friction work - W_normal - W_tangent - Q_other.
    pub total_energy_residual_j: f64,
}

#[derive(Clone)]
pub struct NormalDirectionalContact<'a, M: NormalModulus = AffineNormalModulus> {
    normal: CaloricNormalContact<'a, M>,
    data: NormalDirectionalInterface<'a>,
    law: Arc<SourceThermalEllipticCoulomb>,
    state: SourceThermalEllipticState,
    distance_m: f64,
    work: WorkLedger,
}

pub struct NormalDirectionalCandidate<'a, 't, S, M: NormalModulus = AffineNormalModulus> {
    thermal: DepositionSimulation<'t>,
    contact: NormalDirectionalContact<'a, M>,
    report: NormalDirectionalReport<S>,
}
impl<'a, 't, S, M: NormalModulus> NormalDirectionalCandidate<'a, 't, S, M> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &NormalDirectionalContact<'a, M> {
        &self.contact
    }
    pub const fn report(&self) -> &NormalDirectionalReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        NormalDirectionalContact<'a, M>,
        NormalDirectionalReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

fn value<E>(x: f64) -> Result<f64, DirectionalEndpointError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DirectionalEndpointError::Invalid(
            "nonfinite normal/sliding coordinate or work",
        ))
    }
}

impl<'a, M: NormalModulus> NormalDirectionalContact<'a, M> {
    /// Attach a fresh directional history to an existing original normal owner.
    /// The exact ordered interface and original finite caloric objects must agree.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        normal: CaloricNormalContact<'a, M>,
        data: NormalDirectionalInterface<'a>,
        law: SourceThermalEllipticCoulomb,
        initial_speed_mps: f64,
        tangential: TangentialSolveConfig,
    ) -> Result<Self, NormalDirectionalError> {
        if thermal.time_s() != normal.time_s
            || data.tangent_frame_id.trim().is_empty()
            || data.partition.other() > 0.0
                && data.outside_destination.is_none_or(|s| s.trim().is_empty())
        {
            return Err(NormalThermalError::Invalid(
                "joint normal/sliding contact requires synchronized owners and named work destinations",
            ));
        }
        if temperatures::<DirectionalEndpointError>(thermal, normal.curves, normal.cells)?
            != normal.temperature_k
        {
            return Err(NormalThermalError::Invalid(
                "normal and caloric temperatures must retain the same accepted state",
            ));
        }
        for source in law.principal_sources() {
            if source.interface().card_identity() != data.point.interface.content_hash()
                || source.interface().interface() != &normal.model.reference_request().interface
                || !data
                    .point
                    .is_thermal_axis(&source.axes().temperature, Temperature::DIMS)
                || source.query_axes().iter().any(|a| {
                    a.kind == FrictionQueryAxisKind::Temperature
                        && !data.point.is_thermal_axis(&a.name, Temperature::DIMS)
                })
            {
                return Err(NormalThermalError::Invalid(
                    "normal and directional source must retain the same ordered interface and actual thermal owners",
                ));
            }
        }
        let mut contact = Self {
            normal,
            data,
            state: law.initial_state(),
            law: Arc::new(law),
            distance_m: 0.0,
            work: WorkLedger::default(),
        };
        contact.state = contact
            .evaluate(
                thermal,
                &contact.normal.endpoint,
                [0.0; 2],
                1.0,
                initial_speed_mps,
                0.0,
                tangential,
            )
            .map_err(NormalThermalError::Source)?
            .0
            .state;
        Ok(contact)
    }
    pub const fn normal(&self) -> &CaloricNormalContact<'a, M> {
        &self.normal
    }
    pub const fn binding(&self) -> &NormalDirectionalInterface<'a> {
        &self.data
    }
    pub const fn state(&self) -> &SourceThermalEllipticState {
        &self.state
    }
    pub fn source_law(&self) -> &SourceThermalEllipticCoulomb {
        &self.law
    }
    pub const fn accumulated_distance_m(&self) -> f64 {
        self.distance_m
    }
    pub const fn work(&self) -> &WorkLedger {
        &self.work
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        patch: &PointNormalPatchReceipt,
        slip: [f64; 2],
        dt: f64,
        speed: f64,
        distance: f64,
        solve: TangentialSolveConfig,
    ) -> Result<(TangentialStep<SourceThermalEllipticState>, f64, f64), DirectionalEndpointError<E>>
    {
        if !speed.is_finite() || speed < 0.0 || !distance.is_finite() || distance < self.distance_m
        {
            return Err(DirectionalEndpointError::Invalid(
                "invalid physical sliding speed or accumulated path",
            ));
        }
        let area = value(
            core::f64::consts::PI
                * patch
                    .elliptic_patch_axes
                    .map_or(patch.patch_radius_m * patch.patch_radius_m, |axes| {
                        axes.semi_major_axis_m * axes.semi_minor_axis_m
                    }),
        )?;
        if area < 0.0 || patch.normal_force_n > 0.0 && area <= 0.0 {
            return Err(DirectionalEndpointError::Invalid(
                "positive normal load requires its actual finite patch area",
            ));
        }
        let pressure = if patch.normal_force_n == 0.0 {
            0.0
        } else {
            value(patch.normal_force_n / area)?
        };
        if patch.normal_force_n > 0.0 && pressure == 0.0 {
            return Err(DirectionalEndpointError::Invalid(
                "normal pressure underflow",
            ));
        }
        let phase = |i: usize| {
            let cell = thermal.cells().get(&self.normal.cells[i]).ok_or(
                DirectionalEndpointError::Invalid("missing original caloric cell"),
            )?;
            self.normal.curves[i]
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(DirectionalEndpointError::Phase)
        };
        let mut point = self
            .data
            .point
            .at([phase(0)?, phase(1)?])
            .map_err(DirectionalEndpointError::Material)?;
        let source = &self.law.principal_sources()[0];
        if let Some(axis) = &source.axes().nominal_pressure {
            point = point.with(axis, pressure).map_err(|e| {
                DirectionalEndpointError::Material(ThermalMaterialError::Coordinate(e.to_string()))
            })?;
        }
        for axis in source.query_axes() {
            let coordinate = match axis.kind {
                FrictionQueryAxisKind::Temperature => continue,
                FrictionQueryAxisKind::NominalPressure => pressure,
                FrictionQueryAxisKind::SlidingSpeed => speed,
                FrictionQueryAxisKind::AccumulatedDistance => distance,
            };
            point = point.with(&axis.name, coordinate).map_err(|e| {
                DirectionalEndpointError::Material(ThermalMaterialError::Coordinate(e.to_string()))
            })?;
        }
        let temperature = point
            .axes()
            .get(&source.axes().temperature)
            .copied()
            .ok_or(DirectionalEndpointError::Invalid(
                "missing original friction temperature coordinate",
            ))?;
        let input = TangentialStepInput {
            normal_n: patch.normal_force_n,
            slip_increment_m: slip,
            dt_s: dt,
            compliance_m_per_n: 0.0,
            contact_temperature_k: Some(temperature),
            nominal_pressure_pa: source.axes().nominal_pressure.as_ref().map(|_| pressure),
            solve,
        };
        let friction = self
            .law
            .advance_at_query(
                &self.state,
                &input,
                &point,
                EllipticFrictionMode::StaticOrSliding,
            )
            .map_err(DirectionalEndpointError::Friction)?;
        Ok((friction, area, pressure))
    }

    /// Prescribed indentation and actual tangent travel share one accepted clock.
    /// `resolve` supplies only OTHER heat/conductances; neither contact channel
    /// may be supplied a second time. Every source/refusal preserves all owners.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        indentation_m: f64,
        relative_increment_m: [f64; 2],
        controls: NormalThermalSolve,
        tangential: TangentialSolveConfig,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalDirectionalCandidate<'a, 't, S, M>, NormalDirectionalError<E>> {
        if cancelled() {
            return Err(NormalThermalError::Cancelled);
        }
        let dt = finite(target_s - self.normal.time_s)?;
        let travel = finite(det::hypot(relative_increment_m[0], relative_increment_m[1]))?;
        if dt <= 0.0 {
            return Err(NormalThermalError::Invalid(
                "positive joint normal/sliding duration required",
            ));
        }
        let speed = finite(travel / dt)?;
        let distance = finite(self.distance_m + travel)?;
        if travel > 0.0 && (speed == 0.0 || distance == self.distance_m) {
            return Err(NormalThermalError::Invalid(
                "joint sliding history increment is unrepresentable",
            ));
        }
        let cancel = RefCell::new(cancelled);
        let candidate = self.normal.step_candidate_with_patch(
            thermal,
            target_s,
            indentation_m,
            controls,
            guess_j,
            |state, patch, _| {
                let (friction, area, pressure) = self.evaluate(
                    state,
                    patch,
                    relative_increment_m,
                    dt,
                    speed,
                    distance,
                    tangential,
                )?;
                let work = value(
                    value(friction.traction_n[0] * friction.slip_increment_m[0])?
                        + value(friction.traction_n[1] * friction.slip_increment_m[1])?,
                )?;
                if value(friction.stored_change_j + friction.dissipated_j - work)?.abs()
                    > controls.total_energy_tolerance_j
                {
                    return Err(DirectionalEndpointError::Invalid(
                        "physical directional traction/slip work does not close",
                    ));
                }
                let partition = friction
                    .dissipation_step(self.data.partition)
                    .map_err(DirectionalEndpointError::Tribology)?;
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = resolve(state).map_err(DirectionalEndpointError::Source)?;
                let external_heat = value(
                    dt * state
                        .external_step_power_w(&sources_w, &conductances_w_k, || {
                            if (*cancel.borrow_mut())() {
                                Err(DepositionError::Cancelled)
                            } else {
                                Ok(())
                            }
                        })
                        .map_err(DirectionalEndpointError::Thermal)?,
                )?;
                let mut deposited = 0.0;
                for (id, heat) in self
                    .normal
                    .cells
                    .into_iter()
                    .zip([partition.surface_a_heat_j(), partition.surface_b_heat_j()])
                {
                    let power = value(heat / dt)?;
                    if heat != 0.0 && power == 0.0 {
                        return Err(DirectionalEndpointError::Invalid(
                            "frictional heat underflows power conversion",
                        ));
                    }
                    deposited = value(deposited + value(dt * power)?)?;
                    let total = sources_w.entry(id).or_insert(0.0);
                    *total = value(*total + power)?;
                }
                if value(deposited + partition.other_work_j() - partition.total_work_j())?.abs()
                    > controls.total_energy_tolerance_j
                {
                    return Err(DirectionalEndpointError::Invalid(
                        "directional heat destinations do not close",
                    ));
                }
                Ok(StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: NormalDirectionalEndpoint {
                        friction,
                        partition,
                        contact_area_m2: area,
                        pressure_pa: pressure,
                        tangential_work_j: work,
                        external_heat_j: external_heat,
                        source: payload,
                    },
                })
            },
            || (*cancel.borrow_mut())(),
        )?;
        let (next, normal, report) = candidate.into_parts();
        let tangent = &report.source;
        let total_mechanical_work_j =
            finite(report.work.work_into_contact_j + tangent.tangential_work_j)?;
        let residual = finite(
            report.thermal.stored_energy_change_j
                + report.work.internal_energy_change_j
                + tangent.friction.stored_change_j
                + tangent.partition.other_work_j()
                - total_mechanical_work_j
                - tangent.external_heat_j,
        )?;
        if residual.abs() > controls.total_energy_tolerance_j {
            return Err(NormalThermalError::Invalid(
                "joint normal/sliding/finite-heat work balance failed",
            ));
        }
        let mut contact = self.clone();
        contact.normal = normal;
        contact.state = tangent.friction.state.clone();
        contact.distance_m = distance;
        contact
            .work
            .record(tangent.partition)
            .map_err(|e| NormalThermalError::Source(DirectionalEndpointError::Tribology(e)))?;
        let report = NormalDirectionalReport {
            tangential_work_j: tangent.tangential_work_j,
            total_mechanical_work_j,
            external_heat_j: tangent.external_heat_j,
            total_energy_residual_j: residual,
            normal: report,
        };
        if (*cancel.borrow_mut())() {
            return Err(NormalThermalError::Cancelled);
        }
        Ok(NormalDirectionalCandidate {
            thermal: next,
            contact,
            report,
        })
    }
}
