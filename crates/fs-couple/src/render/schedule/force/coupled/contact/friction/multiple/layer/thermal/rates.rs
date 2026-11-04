//! Source-selected material mobilities on the original finite-body/contact clock.
//!
//! Accepted rates are held over the mechanical interval. The same staged
//! trajectory supplies physical material heat to the original enthalpy solve;
//! endpoint temperatures select rates for the NEXT interval. This first-order
//! split needs time refinement. It is not fully implicit tau(T) or eta(T).
//! Maxwell storage does not depend on tau; viscosity is purely dissipative.
//! Updating these mobilities leaves every mechanical and material state intact
//! and supplies no parameter work or reversible caloric heat.
use super::*;
use crate::render::schedule::force::coupled::coefficients::{
    heat::{MaterialHeatIntegral, viscosity::MaterialViscousHeatSource},
    rates::{
        MaterialRelaxationTime, MaterialRelaxationTimeReport, PreparedMaterialRelaxationTimes,
        mobility::ModalDampingRate,
    },
};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Explicit native-memory to source-time association for ONE original finite
/// material body. `cell` names an existing MaterialRelaxationHeatDestination;
/// `times` covers all of its original linear memory ranges in native order.
/// The physical producer supplies these selectors, never a guessed tau name.
/// This binding also selects the channel's original private viscosity/property
/// rows; `times` may be empty for a purely Kelvin-Voigt material body.
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
    normals: Option<obstacles::geometry::PreparedMaterialMobility>,
    sources: Vec<(usize, ResolvedMaterialStatePoint)>,
    channels: Vec<Option<MaterialHeatIntegral>>,
}
struct MobilityProposal {
    times: Vec<MaterialRelaxationTime>,
    damping: Vec<ModalDampingRate>,
    sources: Vec<(usize, ResolvedMaterialStatePoint)>,
    viscosity: Vec<Vec<MaterialViscousHeatSource>>,
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
            || system.endpoint_controls.is_some()
            || system.core.has_nonlinear_potentials()
        {
            return Err(bad(
                "material rate feedback requires explicit finite-body routes and original linear material storage",
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
            let ranges = channel.groups();
            // Native bending/twist/axial blocks need not be adjacent. The
            // original integral already admits ordered disjoint ranges; the
            // source selectors must cover their exact union, including gaps.
            let has_viscosity = channel.viscosity_sources().iter()
                .any(|source| source.property_names().next().is_some());
            if (ranges.is_empty() && !has_viscosity) || ranges.iter().any(|range| range.is_empty())
                || binding.times.len() != ranges.iter().map(|range| range.len()).sum::<usize>()
                || binding.times.iter().zip(ranges.iter().flat_map(|range| range.clone()))
                    .any(|((branch, name), expected)| {
                    *branch != expected
                        || name.trim().is_empty()
                        || !branches.insert(*branch)
                })
            {
                return Err(bad(
                    "material mobility needs original viscosity rows and/or every original memory's ordered time selector",
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
        let proposal = owner.propose(system, routes, thermal, gate)?;
        for update in proposal.times {
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
        if !proposal.damping.is_empty() || routes.iter().zip(&proposal.viscosity)
            .any(|(route, viscosity)| route.channel.as_ref()
                .is_some_and(|channel| channel.viscosity_sources() != viscosity.as_slice())) {
            return Err(bad("accepted native viscosity differs from its pinned material source"));
        }
        Ok(owner)
    }

    pub(super) fn admit_route<E>(
        &self,
        route: &MaterialRelaxationHeatDestination<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), LayerThermalError<E>> {
        let binding = self.bindings.iter().find(|binding| binding.cell == route.cell);
        let properties: Vec<_> = binding
            .into_iter()
            .flat_map(|binding| binding.times.iter().map(|(_, name)| name.as_str()))
            .collect();
        let viscosity: Vec<_> = binding.into_iter()
            .flat_map(|_| route.channel.iter())
            .flat_map(|channel| channel.viscosity_sources())
            .flat_map(|source| source.property_names()).collect();
        route.admit_with_rates(thermal, &properties, &viscosity)
    }

    fn propose<E>(
        &self,
        system: &MultiFrictionModalSystem,
        routes: &[MaterialRelaxationHeatDestination<'_>],
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<MobilityProposal, LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let mut updates = Vec::new();
        let mut sources = Vec::with_capacity(self.bindings.len());
        let mut viscosity: Vec<Vec<_>> = routes.iter().map(|route| route.channel.as_ref()
            .map_or_else(Vec::new, |channel| channel.viscosity_sources().to_vec())).collect();
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
            viscosity[index] = viscosity[index].iter().map(|row| row.at_material_source(&source))
                .collect::<Result<_, _>>().map_err(ModalFrictionError::from)?;
            sources.push((index, source));
        }
        updates.sort_unstable_by_key(|update| update.branch);
        let damping = system.core.network.material_viscosity_updates(
            routes.iter().zip(&viscosity).flat_map(|(route, proposed)| route.channel.iter()
                .flat_map(move |channel| channel.viscosity_sources().iter().zip(proposed))), gate)
            .map_err(ModalFrictionError::from)?;
        Ok(MobilityProposal { times: updates, damping, sources, viscosity })
    }

    pub(super) fn prepare_endpoint<E>(
        &self,
        system: &mut MultiFrictionModalSystem,
        routes: &[MaterialRelaxationHeatDestination<'_>],
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<RateEndpoint, LayerThermalError<E>> {
        let proposal = self.propose(system, routes, thermal, gate)?;
        let mut native = system
            .core
            .network
            .prepare_staged_material_mobility(&proposal.times, &proposal.damping, gate)
            .map_err(ModalFrictionError::from)?;
        let viscosity_changed = routes.iter().zip(&proposal.viscosity).any(|(route, proposed)|
            route.channel.as_ref().is_some_and(|channel| channel.viscosity_sources() != proposed.as_slice()));
        let (maps, normals, channels) = if native.changed() || viscosity_changed {
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
                let normals = system.obstacles.as_ref().map(|set|
                    set.prepare_material_mobility(network, &system.points, gate)).transpose()?;
                let channels = routes.iter().zip(&proposal.viscosity).map(|(route, viscosity)| route.channel.as_ref().map(|channel|
                    channel.prepare_after_material_mobility(network, &proposal.times, &proposal.damping, viscosity, gate)
                        .map_err(ModalFrictionError::from)).transpose()).collect::<Result<Vec<_>, ModalFrictionError>>()?;
                Ok::<_, ModalFrictionError>((Some(maps), normals, channels))
            })?
        } else {
            (None, None, (0..routes.len()).map(|_| None).collect())
        };
        Ok(RateEndpoint {
            native,
            maps,
            normals,
            sources: proposal.sources,
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
        if let Some(normals) = endpoint.normals {
            system.obstacles.as_mut().expect("prepared original normal contacts")
                .install_material_mobility(normals);
        }
        *self.publication.borrow_mut() = Some(RatePublication {
            report,
            sources: endpoint.sources,
            channels: endpoint.channels,
        });
    }

    fn complete<S>(self, material_heat: &mut [MaterialRelaxationHeatDestination<'_>],
        frame: LayerThermalFrame<S>) -> (LayerThermalFrame<S>, MaterialRelaxationTimeReport) {
        // Native acceptance and the final cancellation check already occurred.
        // These source/cache moves cannot fail or split the accepted timeline.
        let accepted = self.publication.into_inner().expect("accepted material rate transaction");
        for (index, source) in accepted.sources { material_heat[index].source = source; }
        for (route, channel) in material_heat.iter_mut().zip(accepted.channels) {
            if let Some(channel) = channel { route.channel = Some(channel); }
        }
        (frame, accepted.report)
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

    /// Original finite-body/contact heat with explicit source tau(T)/eta(T) feedback.
    /// Accepted rates drive this complete interval; endpoint rates drive the
    /// next. All propagation, contact response and source heat-integral caches
    /// prepare before mechanics, heat, source snapshots or rates publish.
    /// Stored moduli, density, geometry and material internal strains remain
    /// unchanged. No endpoint-power approximation supplies material heat.
    /// Refine this first-order partition in time. Existing solved normal
    /// contacts retain their laws, stored energy and shared traction solve.
    /// Moving tangential maps use `step_with_material_rate_thermal_geometry`;
    /// changing nonlinear material storage needs a separate work composition.
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
        self.step_with_material_rate_thermal_feedback_and_obstacles(thermal, bindings,
            series_heat, dry, material_heat, rates, controls, energy_guess_j, external, drives,
            None, None, gate, resolve)
    }

    /// The same source-mobility transaction with original normal endpoint maps
    /// and support translations. Normal force, storage, loss and actuator work
    /// remain in the existing joint contact solve. A contact-loaded normal may
    /// migrate only with its tangent through the complete geometry entry point.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rate_thermal_feedback_and_obstacles<'t, S, E>(
        &mut self, thermal: &mut DepositionSimulation<'t>, bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination], dry: &[DryContactHeat<'_>],
        material_heat: &mut [MaterialRelaxationHeatDestination<'_>],
        rates: &[MaterialRelaxationTemperature], controls: LayerThermalSolve,
        energy_guess_j: &[f64], external: &[f64], drives: &[FrictionDrive],
        normal_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        normal_offsets_m: Option<&[f64]>, gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
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
            normal_geometry,
            normal_offsets_m,
            false,
            gate,
            resolve,
        )?;
        Ok(continuation.complete(material_heat, frame))
    }

    /// Source-temperature body mobilities with the original held tangential
    /// map transaction. Its actual converged contact effort drives the SAME
    /// source-heat trajectory. Endpoint eta/tau refresh rebuilds responses from
    /// these current complete maps before any state or clock is published.
    ///
    /// The geometry owner retains source-domain checks at the actual start/end
    /// maps and restores the accepted geometry on every later thermal/source/
    /// cache refusal. Motion is a power-conjugate velocity transformer, not an
    /// extra displacement or Bdot*q excitation. Normal maps and offsets hold
    /// their accepted values; the complete companion accepts declared normal
    /// motion. Retained film patches and source/body admission are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rate_thermal_geometry<'t, S, E>(
        &mut self, thermal: &mut DepositionSimulation<'t>, bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination], dry: &[DryContactHeat<'_>],
        material_heat: &mut [MaterialRelaxationHeatDestination<'_>],
        rates: &[MaterialRelaxationTemperature], controls: LayerThermalSolve,
        energy_guess_j: &[f64], external: &[f64], drives: &[FrictionDrive],
        geometry: &[surfaces::geometry::FrictionGeometry],
        endpoints: Option<&[surfaces::geometry::FrictionEndpointGeometry]>,
        config: surfaces::geometry::FrictionGeometryConfig, gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
    ) -> Result<(LayerThermalFrame<S>, MaterialRelaxationTimeReport), LayerThermalError<E>> {
        self.step_with_material_rate_thermal_geometry_and_obstacles(thermal, bindings,
            series_heat, dry, material_heat, rates, controls, energy_guess_j, external, drives,
            geometry, endpoints, None, None, config, gate, resolve)
    }

    /// Complete moving normal/tangent patches and source-temperature body heat
    /// on one accepted clock. Endpoint eta/tau refresh prepares BOTH mixed
    /// response directions and normal compliance at the accepted endpoint
    /// geometry; no lagged normal force or previous-midpoint cache is reused.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rate_thermal_geometry_and_obstacles<'t, S, E>(
        &mut self, thermal: &mut DepositionSimulation<'t>, bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination], dry: &[DryContactHeat<'_>],
        material_heat: &mut [MaterialRelaxationHeatDestination<'_>],
        rates: &[MaterialRelaxationTemperature], controls: LayerThermalSolve,
        energy_guess_j: &[f64], external: &[f64], drives: &[FrictionDrive],
        geometry: &[surfaces::geometry::FrictionGeometry],
        endpoints: Option<&[surfaces::geometry::FrictionEndpointGeometry]>,
        normal_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        normal_offsets_m: Option<&[f64]>, config: surfaces::geometry::FrictionGeometryConfig,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
    ) -> Result<(LayerThermalFrame<S>, MaterialRelaxationTimeReport), LayerThermalError<E>> {
        // Admit the accepted source/operator before staging endpoint controls.
        // Tangential maps do not change the native linear storage generator;
        // their work-conjugate contact forces enter its actual held input.
        let continuation = RateContinuation::new(self, material_heat, rates, thermal, gate)?;
        let frame = self.with_friction_geometry_and_endpoints(geometry, endpoints, normal_geometry,
            config, gate, |system| system.step_with_contact_thermal_inner(thermal, bindings,
                series_heat, dry, material_heat, Some(&continuation), controls, energy_guess_j,
                external, drives, normal_geometry, normal_offsets_m, normal_geometry.is_some(), gate, resolve))?;
        Ok(continuation.complete(material_heat, frame))
    }
}
