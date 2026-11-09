//! Finite-volume dry-contact heating in material coordinates.
//!
//! The interface is a zero-capacity, isothermal plane. Its temperature is fixed
//! by the actual friction power and both bodies' k*A/d half-cell fluxes. Signed
//! cell heat therefore sums to friction loss; bodies can exchange heat at zero
//! slip. This is a resolved finite-cell approximation, not a semi-infinite
//! effusivity split or a flash-temperature formula. A translating footprint
//! visits the SAME original cells; accepted heat is never reset or advected to
//! an unrelated cell. Refine cells, footprint and time together.
//!
//! Original Robin boundaries cool only the uncovered material surface. A loaded
//! footprint shields its actual overlap; lifting restores that area to the
//! original ambient bath. The native boundary owner retains h, bath temperature
//! and external heat accounting. Reservoir faces require a separate geometric
//! correspondence and remain unsupported by this strip adapter.

use super::*;
use fs_material::state_point::{ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// One original solid rectangular finite cell beneath a material surface strip.
/// `half_depth_m` is its centre-to-interface distance, not a fitted resistance.
/// The declared strip width times interval length times full depth must equal
/// the original thermal cell's volume. Geometry is in the body's own frame.
#[derive(Clone, Copy, Debug)]
pub struct DryThermalCell<'a> {
    pub cell: CellId,
    pub material: ThermalMaterialPoint<'a>,
    pub caloric: &'a EquilibriumEnthalpyPhaseCurve,
    /// Original pinned conductivity AND density claims. Every trial refreshes
    /// these same claims; a neighboring claim cannot take over outside support.
    pub transport: &'a ResolvedMaterialStatePoint,
    pub conductivity_property: &'a str,
    pub density_property: &'a str,
    pub interval_m: [f64; 2],
    pub width_m: f64,
    pub half_depth_m: f64,
}

/// Ordered material-coordinate footprints on the two original participants.
/// Their areas must both equal the dry interface's original nominal area.
/// They are held over one sample, like the work-conjugate mechanical map.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DryContactFootprint {
    pub intervals_m: [[f64; 2]; 2],
    pub widths_m: [f64; 2],
}

/// A dry port and two finite receiving surfaces. Each strip is sorted,
/// contiguous and nonoverlapping. Multiple ports may heat the same original
/// cell, but one route cannot alias a cell on both sides. An original thermal
/// link between these sides would duplicate contact conduction and refuses.
/// The producer maps ordered source-interface participants to these bodies.
#[derive(Clone, Debug)]
pub struct DryContactHeat<'a> {
    pub port: usize,
    pub surfaces: [Vec<DryThermalCell<'a>>; 2],
    pub footprint: DryContactFootprint,
}

/// Accepted physical heat of a single dry contact. A released contact has no
/// interface temperature or heat exchange. Negative cell heat is cooling.
#[derive(Clone, Debug, PartialEq)]
pub struct DryContactHeatFrame {
    pub port: usize,
    pub interface_temperature_k: Option<f64>,
    pub friction_heat_j: f64,
    /// Sum of actual body fluxes minus the same staged friction work [J].
    /// Independently admitted before either mechanics or heat is published.
    pub interface_energy_residual_j: f64,
    pub cell_heat_j: Vec<(CellId, f64)>,
    pub footprint: DryContactFootprint,
}

fn positive<E>(v: f64) -> Result<f64, LayerThermalError<E>> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(LayerThermalError::Invalid(
            "dry thermal geometry and conductance must be positive and finite",
        ))
    }
}

impl DryThermalCell<'_> {
    pub(super) fn conductance<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        area: f64,
    ) -> Result<(f64, f64), LayerThermalError<E>> {
        self.conductance_with_ambient(thermal, area, false)
    }

    fn conductance_with_ambient<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        area: f64,
        allow_ambient: bool,
    ) -> Result<(f64, f64), LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(bad("dry heat cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, self.caloric)
            || self.transport.card_identity() != self.material.material.content_hash()
            || !allow_ambient && cell.boundary.htc_w_m2_k != 0.0
            || cell.boundary.reservoir.is_some()
        {
            return Err(bad(
                "contact heat requires original caloric/card identity, no unresolved reservoir face, and explicit ambient coverage",
            ));
        }
        let phase = self
            .caloric
            .state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if !phase.is_fully_solid() {
            return Err(bad("moving dry thermal cells must remain solid"));
        }
        let point = self
            .material
            .at(phase)
            .map_err(LayerThermalError::ThermalMaterial)?;
        let resolved = self
            .transport
            .at_query_point(self.material.material, &point)?;
        let k = resolved
            .property(self.conductivity_property)
            .ok_or(bad("pinned scalar thermal conductivity is absent"))?;
        let rho = resolved
            .property(self.density_property)
            .ok_or(bad("pinned reference density is absent"))?;
        let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        if k.requirement().quantity() != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && k.requirement().quantity() != semantic(QuantityKind::ThermalConductivity)
            || rho.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
        {
            return Err(bad(
                "dry heat needs scalar static conductivity and density quantities",
            ));
        }
        self.material
            .admit_static(k)
            .map_err(LayerThermalError::ThermalMaterial)?;
        self.material
            .admit_static(rho)
            .map_err(LayerThermalError::ThermalMaterial)?;
        if !rho.source_axes().is_empty() {
            return Err(bad("moving dry heat requires fixed reference density"));
        }
        let density = positive(rho.value_si())?;
        if self
            .caloric
            .knots()
            .iter()
            .any(|k| !same_size(k.bulk_density_kg_m3, density))
            || !same_size(cell.mass_kg, density * cell.geometry.volume_m3)
        {
            return Err(bad(
                "dry heat cell and caloric chart must retain the source reference mass",
            ));
        }
        let length = positive(self.interval_m[1] - self.interval_m[0])?;
        if self.interval_m.iter().any(|x| !x.is_finite())
            || !same_size(
                cell.geometry.volume_m3,
                length * positive(self.width_m)? * 2.0 * positive(self.half_depth_m)?,
            )
        {
            return Err(bad(
                "dry thermal strip must describe the actual finite-cell volume",
            ));
        }
        Ok((
            positive(positive(k.value_si())? * area / self.half_depth_m)?,
            phase.temperature_k(),
        ))
    }
}

impl DryContactHeat<'_> {
    // Coordinate translation can place a geometrically coincident edge a few
    // ulps into its neighbor. Snap only representational roundoff to ORIGINAL
    // cell boundaries, so a revisited cell does not acquire a fictitious sliver.
    fn resolved_footprint(&self) -> DryContactFootprint {
        let mut footprint = self.footprint;
        for (side, cells) in self.surfaces.iter().enumerate() {
            for edge in &mut footprint.intervals_m[side] {
                if !edge.is_finite() {
                    continue;
                }
                for boundary in cells.iter().flat_map(|cell| cell.interval_m) {
                    let scale = edge.abs().max(boundary.abs());
                    if boundary.is_finite()
                        && (*edge - boundary).abs() <= 8.0 * f64::EPSILON * scale
                    {
                        *edge = boundary;
                        break;
                    }
                }
            }
        }
        footprint
    }
    /// Validate the original interface, finite caloric cells, source masses,
    /// material-coordinate coverage and nonduplicated conduction before output
    /// or stepping. The thermal solve repeats these same checks on each step.
    pub fn admit<E>(
        &self,
        system: &MultiFrictionModalSystem,
        thermal: &DepositionSimulation<'_>,
        maximum_cells: usize,
    ) -> Result<(), LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let binding = system
            .interface(self.port)
            .ok_or(bad("dry heating needs the original dry-interface source"))?;
        if system
            .shear_layers
            .get(self.port)
            .and_then(Option::as_ref)
            .is_some()
        {
            return Err(bad(
                "a dry heat route cannot replace retained shear-film history",
            ));
        }
        let footprint = self.resolved_footprint();
        let mut all = BTreeSet::new();
        let mut side_cells = [BTreeSet::new(), BTreeSet::new()];
        for (side, cells) in self.surfaces.iter().enumerate() {
            if cells.is_empty() || cells.len() > maximum_cells {
                return Err(bad("dry contact requires bounded finite surface strips"));
            }
            let interval = footprint.intervals_m[side];
            let area =
                positive(interval[1] - interval[0])? * positive(self.footprint.widths_m[side])?;
            let coordinate_roundoff = 64.0
                * f64::EPSILON
                * self.footprint.widths_m[side]
                * interval[0]
                    .abs()
                    .max(interval[1].abs())
                    .max(interval[1] - interval[0]);
            if interval.iter().any(|x| !x.is_finite())
                || (area - binding.nominal_area_m2()).abs() > coordinate_roundoff
            {
                return Err(bad(
                    "both dry thermal footprints must equal the source nominal contact area",
                ));
            }
            let mut previous = None;
            for cell in cells {
                if !all.insert(cell.cell)
                    || previous.is_some_and(|end| end != cell.interval_m[0])
                    || cell.width_m < self.footprint.widths_m[side]
                {
                    return Err(bad(
                        "dry thermal strips must have unique cells, contiguous intervals and sufficient width",
                    ));
                }
                side_cells[side].insert(cell.cell);
                // Admit every cell, including currently uncovered history, and
                // pin its source before any mechanical candidate is published.
                cell.conductance_with_ambient::<E>(thermal, area, true)?;
                previous = Some(cell.interval_m[1]);
            }
            if interval[0] < cells[0].interval_m[0]
                || interval[1] > cells.last().expect("nonempty").interval_m[1]
            {
                return Err(bad("dry footprint leaves the resolved material surface"));
            }
        }
        if thermal.links().values().any(|link| {
            side_cells[0].contains(&link.a) && side_cells[1].contains(&link.b)
                || side_cells[0].contains(&link.b) && side_cells[1].contains(&link.a)
        }) {
            return Err(bad(
                "an existing link duplicates the dry interface's half-cell conduction",
            ));
        }
        Ok(())
    }

    pub(super) fn evaluate<E>(
        &self,
        endpoint: &DepositionSimulation<'_>,
        normal_force_n: f64,
        friction_heat_j: f64,
        dt: f64,
        tolerance_j: f64,
    ) -> Result<DryContactHeatFrame, LayerThermalError<E>> {
        let footprint = self.resolved_footprint();
        let mut frame = DryContactHeatFrame {
            port: self.port,
            interface_temperature_k: None,
            friction_heat_j,
            interface_energy_residual_j: 0.0,
            cell_heat_j: Vec::new(),
            footprint,
        };
        if !friction_heat_j.is_finite() || friction_heat_j < 0.0 {
            return Err(LayerThermalError::Invalid(
                "dry physical friction heat must be nonnegative",
            ));
        }
        if normal_force_n == 0.0 {
            if friction_heat_j != 0.0 {
                return Err(LayerThermalError::Invalid(
                    "separated dry contact cannot generate friction heat",
                ));
            }
            return Ok(frame);
        }
        let conductors = self.conductors(endpoint)?;
        let (mut sum_g, mut sum_gt) = (0.0, 0.0);
        for &(_, g, temperature) in &conductors {
            sum_g = number(sum_g + g)?;
            sum_gt = number(sum_gt + g * temperature)?;
        }
        let temperature = positive((friction_heat_j / dt + sum_gt) / positive(sum_g)?)?;
        let mut sum = 0.0;
        for (cell, g, body_temperature) in conductors {
            let heat = number(dt * g * (temperature - body_temperature))?;
            sum = number(sum + heat)?;
            frame.cell_heat_j.push((cell, heat));
        }
        frame.interface_energy_residual_j = number(sum - friction_heat_j)?;
        if frame.interface_energy_residual_j.abs() > tolerance_j {
            return Err(LayerThermalError::Invalid(
                "dry interface fluxes do not close the actual friction loss",
            ));
        }
        frame.interface_temperature_k = Some(temperature);
        Ok(frame)
    }

    fn conductors<E>(&self, endpoint: &DepositionSimulation<'_>)
        -> Result<Vec<(CellId, f64, f64)>, LayerThermalError<E>>
    {
        let footprint = self.resolved_footprint();
        let mut conductors = Vec::new();
        for (side, cells) in self.surfaces.iter().enumerate() {
            let interval = footprint.intervals_m[side];
            for cell in cells {
                let overlap = cell.interval_m[1].min(interval[1]) - cell.interval_m[0].max(interval[0]);
                if overlap > 0.0 {
                    let (g, temperature) = cell.conductance_with_ambient(
                        endpoint, overlap * footprint.widths_m[side], true)?;
                    conductors.push((cell.cell, g, temperature));
                }
            }
        }
        Ok(conductors)
    }

    /// Numerical J/K scaling from the original flux derivative. This creates
    /// no interface heat capacity or additional physical energy store.
    pub(super) fn temperature_scale<E>(&self, thermal: &DepositionSimulation<'_>, dt: f64)
        -> Result<f64, LayerThermalError<E>>
    {
        let mut conductance = 0.0;
        for (_, g, _) in self.conductors(thermal)? { conductance = number(conductance + g)?; }
        positive(dt * positive(conductance)?)
    }

    /// Fluxes at a JOINT nonlinear temperature coordinate. Unlike the fixed-law
    /// analytic elimination above, its interface residual may be nonzero during
    /// Newton trials; the caller must admit it independently before publishing.
    pub(super) fn evaluate_at_temperature<E>(
        &self, endpoint: &DepositionSimulation<'_>, normal_force_n: f64,
        friction_heat_j: f64, dt: f64, temperature: f64,
    ) -> Result<DryContactHeatFrame, LayerThermalError<E>> {
        if normal_force_n == 0.0 {
            return self.evaluate(endpoint, normal_force_n, friction_heat_j, dt, 0.0);
        }
        positive(temperature)?;
        if !friction_heat_j.is_finite() || friction_heat_j < 0.0 {
            return Err(LayerThermalError::Invalid("dry physical friction heat must be nonnegative"));
        }
        let mut frame = DryContactHeatFrame {
            port: self.port, interface_temperature_k: Some(temperature), friction_heat_j,
            interface_energy_residual_j: 0.0, cell_heat_j: Vec::new(),
            footprint: self.resolved_footprint(),
        };
        let mut sum = 0.0;
        for (cell, g, body_temperature) in self.conductors(endpoint)? {
            let heat = number(dt * g * (temperature - body_temperature))?;
            sum = number(sum + heat)?;
            frame.cell_heat_j.push((cell, heat));
        }
        frame.interface_energy_residual_j = number(sum - friction_heat_j)?;
        Ok(frame)
    }
}

struct AmbientFootprint {
    port: usize,
    interval_m: [f64; 2],
    area_m2: f64,
}

struct AmbientSurface {
    interval_m: [f64; 2],
    width_m: f64,
    exposed_area_m2: f64,
    htc_w_m2_k: f64,
    footprints: Vec<AmbientFootprint>,
}

/// Prepared geometry only: temperatures and heat remain with the original
/// thermal owner. Native topology/area is traversed once per physical step.
pub(super) struct DryAmbientBoundaries {
    surfaces: BTreeMap<CellId, AmbientSurface>,
}

impl DryAmbientBoundaries {
    pub(super) fn prepare<E>(
        routes: &[DryContactHeat<'_>],
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let mut surfaces = BTreeMap::<CellId, AmbientSurface>::new();
        for route in routes {
            checkpoint()?;
            let footprint = route.resolved_footprint();
            for (side, cells) in route.surfaces.iter().enumerate() {
                for cell in cells {
                    checkpoint()?;
                    let htc = thermal.cells()[&cell.cell].boundary.htc_w_m2_k;
                    let interval = [
                        cell.interval_m[0].max(footprint.intervals_m[side][0]),
                        cell.interval_m[1].min(footprint.intervals_m[side][1]),
                    ];
                    if htc == 0.0 || interval[1] <= interval[0] {
                        continue;
                    }
                    let surface = surfaces.entry(cell.cell).or_insert_with(|| AmbientSurface {
                        interval_m: cell.interval_m,
                        width_m: cell.width_m,
                        exposed_area_m2: 0.0,
                        htc_w_m2_k: htc,
                        footprints: Vec::new(),
                    });
                    if surface.interval_m != cell.interval_m || surface.width_m != cell.width_m {
                        return Err(bad(
                            "ambient contact coverage needs one original material face per cell",
                        ));
                    }
                    surface.footprints.push(AmbientFootprint {
                        port: route.port,
                        interval_m: interval,
                        area_m2: positive((interval[1] - interval[0]) * footprint.widths_m[side])?,
                    });
                }
            }
        }
        let ids: Vec<_> = surfaces.keys().copied().collect();
        let areas = thermal.exposed_areas_m2(&ids, &mut checkpoint)?;
        for (&id, surface) in &mut surfaces {
            checkpoint()?;
            surface.exposed_area_m2 = areas[&id];
            let face = positive((surface.interval_m[1] - surface.interval_m[0]) * surface.width_m)?;
            if face > surface.exposed_area_m2 && !same_size(face, surface.exposed_area_m2) {
                return Err(bad(
                    "original thermal contacts leave insufficient exposed area for the dry surface",
                ));
            }
            surface.footprints.sort_by(|a, b| {
                a.interval_m[0]
                    .total_cmp(&b.interval_m[0])
                    .then(a.interval_m[1].total_cmp(&b.interval_m[1]))
                    .then(a.port.cmp(&b.port))
            });
        }
        Ok(Self { surfaces })
    }

    pub(super) fn at_ports<E>(
        &self,
        ports: &[FrictionPortFrame],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BoundaryConductances, LayerThermalError<E>> {
        let bad = LayerThermalError::Invalid;
        let mut boundaries = BoundaryConductances::default();
        for (&id, surface) in &self.surfaces {
            checkpoint()?;
            let (mut covered, mut end) = (0.0, None);
            for footprint in &surface.footprints {
                checkpoint()?;
                let normal = ports
                    .get(footprint.port)
                    .ok_or(bad("ambient contact port is absent"))?
                    .normal_force_n;
                if !normal.is_finite() || normal < 0.0 {
                    return Err(bad(
                        "ambient contact coverage needs finite nonnegative load",
                    ));
                }
                if normal == 0.0 {
                    continue;
                }
                // Width offsets are not in this strip representation. Two
                // loaded overlapping intervals cannot claim independent faces.
                if end.is_some_and(|previous| footprint.interval_m[0] < previous) {
                    return Err(bad(
                        "overlapping loaded dry footprints share one ambient surface",
                    ));
                }
                end = Some(footprint.interval_m[1]);
                covered = number(covered + footprint.area_m2)?;
            }
            let uncovered = if covered <= surface.exposed_area_m2 {
                surface.exposed_area_m2 - covered
            } else if same_size(covered, surface.exposed_area_m2) {
                // Only representation roundoff at a completely covered face.
                0.0
            } else {
                return Err(bad(
                    "loaded dry footprints exceed the exposed material surface",
                ));
            };
            boundaries
                .ambient_w_k
                .insert(id, number(surface.htc_w_m2_k * uncovered)?);
        }
        Ok(boundaries)
    }
}

impl MultiFrictionModalSystem {
    /// Original dry friction, material-coordinate half-cell conduction, finite
    /// heat capacity and optional series-drive heating on ONE accepted clock.
    /// A released contact switches off both friction heating and conduction.
    /// Original Robin boundaries then cool the fully exposed material surface;
    /// loaded footprints shield only their actual covered area from that bath.
    /// Existing films remain covered only by their fixed-film route; this dry
    /// continuation never transports or resets material shear history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_contact_thermal_network<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        dry: &[DryContactHeat<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        normal_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        normal_offsets_m: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.step_with_contact_thermal_inner(
            thermal,
            &[],
            series_heat,
            dry,
            &[],
            None,
            controls,
            energy_guess_j,
            external,
            drives,
            normal_geometry,
            normal_offsets_m,
            false,
            gate,
            |endpoint, _ports| resolve(endpoint),
        )
    }

    /// The same dry/thermal solve with moving work-conjugate contact maps.
    /// Footprints are in original material coordinates, independently of modal
    /// maps. The geometry producer must supply their physical co-location.
    /// The existing endpoint/normal checks and all rollback rules are retained.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_contact_thermal_geometry<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        dry: &[DryContactHeat<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        geometry: &[surfaces::geometry::FrictionGeometry],
        endpoints: Option<&[surfaces::geometry::FrictionEndpointGeometry]>,
        normal_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        normal_offsets_m: Option<&[f64]>,
        config: surfaces::geometry::FrictionGeometryConfig,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.with_friction_geometry_and_endpoints(
            geometry,
            endpoints,
            normal_geometry,
            config,
            gate,
            |system| {
                system.step_with_contact_thermal_inner(
                    thermal,
                    &[],
                    series_heat,
                    dry,
                    &[],
                    None,
                    controls,
                    energy_guess_j,
                    external,
                    drives,
                    normal_geometry,
                    normal_offsets_m,
                    normal_geometry.is_some(),
                    gate,
                    |endpoint, _ports| resolve(endpoint),
                )
            },
        )
    }

    pub(super) fn validate_dry_temperature<E>(
        &mut self,
        port: usize,
        temperature: f64,
        drive: FrictionDrive,
    ) -> Result<interface::DryFrictionBinding, LayerThermalError<E>> {
        let binding = self
            .interface(port)
            .expect("admitted dry source")
            .at_temperature(temperature)?
            .at_average_motion(self.staged_ports[port].normal_force_n,
                self.staged_ports[port].slip_distance_m / self.sample_period_s())?;
        let original = if port == 0 {
            self.core.interface.replace(binding.clone())
        } else {
            self.interfaces[port - 1].replace(binding.clone())
        };
        let result = (|| {
            let drive = FrictionDrive {
                normal_force_n: self.staged_ports[port].normal_force_n,
                ..drive
            };
            self.check_interface_endpoint(port, &self.core.network.models, drive, false)?;
            interface::check_average(
                self.interface(port),
                drive,
                self.staged_ports[port].slip_distance_m,
                self.sample_period_s(),
            )?;
            self.check_interface_endpoint(port, &self.core.network.candidates, drive, true)
        })();
        if port == 0 {
            self.core.interface = original;
        } else {
            self.interfaces[port - 1] = original;
        }
        result?;
        Ok(binding)
    }
}
