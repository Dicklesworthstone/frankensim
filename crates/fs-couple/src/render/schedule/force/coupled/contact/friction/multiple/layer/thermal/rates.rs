//! Source-selected relaxation times on the original finite-body/contact clock.
//!
//! Accepted rates are held over the mechanical interval. The same staged
//! trajectory supplies physical material heat to the original enthalpy solve;
//! endpoint temperatures select rates for the NEXT interval. This first-order
//! split needs time refinement. It is not fully implicit tau(T) within a step.
//! Maxwell storage and its internal strain do not depend on tau, so changing
//! only this mobility supplies no parameter work or reversible caloric heat.
use super::*;
use crate::render::schedule::force::coupled::coefficients::{
    heat::MaterialHeatIntegral,
    rates::{
        MaterialRelaxationTime, MaterialRelaxationTimeReport, PreparedMaterialRelaxationTimes,
    },
};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Explicit native-memory to source-time association for ONE original finite
/// material body. `cell` names an existing MaterialRelaxationHeatDestination;
/// `times` covers its complete original linear memory range in native order.
/// The physical producer supplies these selectors, never a guessed tau name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaterialRelaxationTemperature {
    pub cell: CellId,
    pub times: Vec<(usize, String)>,
}

pub(super) struct RateContinuation<'a> {
    bindings: &'a [MaterialRelaxationTemperature],
    publication: RefCell<Option<RatePublication>>,
}
struct RatePublication {
    report: MaterialRelaxationTimeReport,
    sources: Vec<(usize, ResolvedMaterialStatePoint)>,
    channels: Vec<Option<MaterialHeatIntegral>>,
}
pub(super) struct RateEndpoint {
    native: PreparedMaterialRelaxationTimes,
    maps: Option<(Vec<f64>, Vec<f64>)>,
    sources: Vec<(usize, ResolvedMaterialStatePoint)>,
    channels: Vec<Option<MaterialHeatIntegral>>,
}

impl RateContinuation<'_> {
    fn new<'a, E>(
        system: &MultiFrictionModalSystem,
        routes: &[MaterialRelaxationHeatDestination<'_>],
        bindings: &'a [MaterialRelaxationTemperature],
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<RateContinuation<'a>, LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        if bindings.is_empty()
            || bindings.len() > routes.len()
            || system.obstacles.is_some()
            || system.endpoint_controls.is_some()
            || system.core.has_nonlinear_potentials()
        {
            return Err(bad(
                "material rate feedback requires explicit finite-body routes and fixed linear contact mechanics",
            ));
        }
        let mut cells = BTreeSet::new();
        let mut branches = BTreeSet::new();
        for binding in bindings {
            poll(Some(gate)).map_err(ModalFrictionError::from)?;
            let route = routes
                .iter()
                .find(|route| route.cell == binding.cell)
                .ok_or(bad(
                    "material rate source lacks its original finite heat destination",
                ))?;
            if !cells.insert(binding.cell)
                || !route.nonlinear_channels.is_empty()
                || route
                    .caloric
                    .knots()
                    .iter()
                    .any(|k| k.liquid_mass_fraction != 0.0)
            {
                return Err(bad(
                    "rate feedback needs distinct fully solid bodies and original linear material histories",
                ));
            }
            let channel = route.channel.as_ref().ok_or(bad(
                "rate feedback requires the original linear heat integral",
            ))?;
            let [range] = channel.groups() else {
                return Err(bad(
                    "rate feedback must retain one complete physical material group",
                ));
            };
            if range.is_empty()
                || binding.times.len() != range.len()
                || binding.times.iter().enumerate().any(|(i, (branch, name))| {
                    *branch != range.start + i
                        || name.trim().is_empty()
                        || !branches.insert(*branch)
                })
            {
                return Err(bad(
                    "every original material memory needs one explicit ordered time selector",
                ));
            }
            if route
                .source
                .state_coordinate(route.temperature_axis)
                .map(f64::to_bits)
                != Some(thermal.temperature_k(route.cell)?.to_bits())
            {
                return Err(bad(
                    "accepted material source temperature must match its original finite body",
                ));
            }
        }
        for route in routes {
            if let Some(channel) = &route.channel {
                if !channel
                    .is_current(&system.core.network, gate)
                    .map_err(ModalFrictionError::from)?
                {
                    return Err(bad(
                        "material rate feedback cannot replace an externally changed source operator",
                    ));
                }
            }
        }
        let owner = RateContinuation {
            bindings,
            publication: RefCell::new(None),
        };
        let (updates, _) = owner.propose(routes, thermal, gate)?;
        for update in updates {
            if system
                .material_relaxation_branches()
                .get(update.branch)
                .is_none_or(|branch| {
                    branch.relaxation_time_s.to_bits() != update.relaxation_time_s.to_bits()
                })
            {
                return Err(bad(
                    "accepted native relaxation rate differs from its pinned material source",
                ));
            }
        }
        Ok(owner)
    }

    pub(super) fn admit_route<E>(
        &self,
        route: &MaterialRelaxationHeatDestination<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), LayerThermalError<E>> {
        let properties: Vec<_> = self
            .bindings
            .iter()
            .find(|binding| binding.cell == route.cell)
            .into_iter()
            .flat_map(|binding| binding.times.iter().map(|(_, name)| name.as_str()))
            .collect();
        route.admit_with_rates(thermal, &properties)
    }

    fn propose<E>(
        &self,
        routes: &[MaterialRelaxationHeatDestination<'_>],
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<
        (
            Vec<MaterialRelaxationTime>,
            Vec<(usize, ResolvedMaterialStatePoint)>,
        ),
        LayerThermalError<E>,
    > {
        let bad = LayerThermalError::Invalid;
        let mut updates = Vec::new();
        let mut sources = Vec::with_capacity(self.bindings.len());
        for binding in self.bindings {
            poll(Some(gate)).map_err(ModalFrictionError::from)?;
            let (index, route) = routes
                .iter()
                .enumerate()
                .find(|(_, route)| route.cell == binding.cell)
                .ok_or(bad("material rate source lost its original finite body"))?;
            self.admit_route(route, thermal)?;
            let query = route
                .source
                .query()
                .with(route.temperature_axis, thermal.temperature_k(route.cell)?)
                .map_err(|_| {
                    bad("material rate temperature must retain its original source coordinates")
                })?;
            let source = route.source.at_query_point(route.material, &query)?;
            for (branch, name) in &binding.times {
                let property = source.property(name).ok_or(bad(
                    "material rate selector is absent from its original source",
                ))?;
                let time = property.value_si();
                if property.requirement().quantity()
                    != QuantitySpec::dimensional(fs_qty::Time::DIMS)
                    || !time.is_finite()
                    || time <= 0.0
                    || !(1.0 / time).is_finite()
                {
                    return Err(bad(
                        "source relaxation time and rate must be positive representable values",
                    ));
                }
                updates.push(MaterialRelaxationTime {
                    branch: *branch,
                    relaxation_time_s: time,
                });
            }
            sources.push((index, source));
        }
        updates.sort_unstable_by_key(|update| update.branch);
        Ok((updates, sources))
    }

    pub(super) fn prepare_endpoint<E>(
        &self,
        system: &mut MultiFrictionModalSystem,
        routes: &[MaterialRelaxationHeatDestination<'_>],
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<RateEndpoint, LayerThermalError<E>> {
        let (updates, sources) = self.propose(routes, thermal, gate)?;
        let mut native = system
            .core
            .network
            .prepare_staged_material_relaxation_times(&updates, gate)
            .map_err(ModalFrictionError::from)?;
        let (maps, channels) = if native.changed() {
            // Use the COMPLETE retained columns: a condensed contact may have
            // string, bridge/body and fiber participants. Reconstructing only
            // the public local left/right geometry would lose bridge work.
            let columns: Vec<_> = system
                .points
                .iter()
                .map(|point| point.column.as_slice())
                .collect();
            system.core.network.preview_material_relaxation_times(&mut native, |network| {
                let maps = port_maps(network, &columns, gate,
                    "updated material rates require positive contact compliance and endpoint mobility")?;
                let channels = routes.iter().map(|route| route.channel.as_ref().map(|channel|
                    channel.prepare_after_material_rates(network, &updates, gate)
                        .map_err(ModalFrictionError::from)).transpose()).collect::<Result<Vec<_>, ModalFrictionError>>()?;
                Ok::<_, ModalFrictionError>((Some(maps), channels))
            })?
        } else {
            (None, (0..routes.len()).map(|_| None).collect())
        };
        Ok(RateEndpoint {
            native,
            maps,
            sources,
            channels,
        })
    }

    pub(super) fn publish(&self, system: &mut MultiFrictionModalSystem, endpoint: RateEndpoint) {
        let report = system
            .core
            .network
            .install_material_relaxation_times(endpoint.native);
        if let Some((compliance, mobility)) = endpoint.maps {
            system.core.compliance = mobility[0];
            system.core.displacement_compliance = compliance[0];
            system.compliance = compliance;
            system.endpoint = mobility;
        }
        *self.publication.borrow_mut() = Some(RatePublication {
            report,
            sources: endpoint.sources,
            channels: endpoint.channels,
        });
    }
}

impl MultiFrictionModalSystem {
    /// Read-only source/body/rate admission for the same thermal continuation.
    /// Reuse this before opening output or binding a performance. It validates
    /// the original operator and pinned material roles without staging motion,
    /// refreshing coefficients, depositing heat or advancing either clock.
    /// Contact-film and dry-heat bindings keep their own admission owners.
    pub fn admit_material_rate_thermal_feedback(
        &self,
        thermal: &DepositionSimulation<'_>,
        material_heat: &[MaterialRelaxationHeatDestination<'_>],
        rates: &[MaterialRelaxationTemperature],
        gate: &CancelGate,
    ) -> Result<(), LayerThermalError<Infallible>> {
        RateContinuation::new(self, material_heat, rates, thermal, gate).map(|_| ())
    }

    /// Original finite-body/contact heat with explicit source tau(T) feedback.
    /// Accepted rates drive this complete interval; endpoint rates drive the
    /// next. All propagation, contact response and source heat-integral caches
    /// prepare before mechanics, heat, source snapshots or rates publish.
    /// Stored moduli, density, geometry and material internal strains remain
    /// unchanged. No endpoint-power approximation supplies material heat.
    /// Refine this first-order partition in time. Moving/normal contacts and
    /// nonlinear material storage require their separate update compositions.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rate_thermal_feedback<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        dry: &[DryContactHeat<'_>],
        material_heat: &mut [MaterialRelaxationHeatDestination<'_>],
        rates: &[MaterialRelaxationTemperature],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
        resolve: impl Fn(
            &DepositionSimulation<'t>,
            &[FrictionPortFrame],
        ) -> Result<StepCoefficients<S>, E>,
    ) -> Result<(LayerThermalFrame<S>, MaterialRelaxationTimeReport), LayerThermalError<E>> {
        let continuation = RateContinuation::new(self, material_heat, rates, thermal, gate)?;
        let frame = self.step_with_contact_thermal_inner(
            thermal,
            bindings,
            series_heat,
            dry,
            material_heat,
            Some(&continuation),
            controls,
            energy_guess_j,
            external,
            drives,
            None,
            None,
            false,
            gate,
            resolve,
        )?;
        // Native acceptance and the final cancellation check already occurred.
        // These source/cache moves cannot fail or split the accepted timeline.
        let accepted = continuation
            .publication
            .into_inner()
            .expect("accepted material rate transaction");
        for (index, source) in accepted.sources {
            material_heat[index].source = source;
        }
        for (route, channel) in material_heat.iter_mut().zip(accepted.channels) {
            if let Some(channel) = channel {
                route.channel = Some(channel);
            }
        }
        Ok((frame, accepted.report))
    }
}
