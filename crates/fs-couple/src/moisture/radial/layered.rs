//! Coatings and composite cylinders/spheres with source-owned partitioning.
//!
//! Every layer is an original radial field with its own dry skeleton, nonlinear
//! diffusivity and monotone sorption source. Adjacent materials share activity
//! at ideal contact, not dry-basis moisture. Finite interface conductance allows
//! an activity jump and carries one equal/opposite water flux across the actual
//! curved area. The shared layered solver owns nonlinear storage and transport;
//! geometry supplies exact integral(dr / area) paths and unequal shell masses.
//!
//! This is fixed-geometry, isothermal bound-water diffusion with radial symmetry.
//! Axial/angular transport, moving interfaces, sorption heat, hysteresis and
//! liquid capillarity require additional admitted physical laws.

use super::*;
use crate::moisture::layered::{
    MoistureActivityTransfer, MoistureLayerInterfaceReport, MoistureLayeredFace,
    MoistureLayeredReport, admit_sorption,
    interface::{self, ActivityPath},
    kernel::{self, LayeredSource, LayeredStorage, LayeredSystem},
};
use fs_material::state_point::MonotoneMoistureSorption;

/// Original finite vapor inventories coupled to complete coated radial bodies.
pub mod atmosphere;
mod sdirk;

/// One native radial material field and its explicitly associated isotherm.
#[derive(Clone, Debug)]
pub struct MoistureRadialLayer {
    radial: MoistureRadial,
    sorption: MonotoneMoistureSorption,
}

impl MoistureRadialLayer {
    /// Both sources must come from the same original material and held external
    /// coordinates. Source support must contain every accepted shell content.
    pub fn new(
        radial: MoistureRadial,
        sorption: MonotoneMoistureSorption,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        admit_sorption(
            &radial.card,
            &radial.transport,
            radial.geometry.maximum_moisture,
            &radial.moisture,
            &sorption,
            cx,
        )?;
        Ok(Self { radial, sorption })
    }

    /// The original native field, available to its material/solid consumers.
    pub const fn radial(&self) -> &MoistureRadial {
        &self.radial
    }

    pub const fn sorption(&self) -> &MonotoneMoistureSorption {
        &self.sorption
    }

    fn dry_density(&self) -> f64 {
        self.radial
            .reference
            .property("density")
            .expect("admitted density")
            .value_si()
            / (1.0 + self.radial.reference.query().axes()["moisture_content"])
    }

    fn source(&self) -> LayeredSource<'_> {
        LayeredSource {
            card: &self.radial.card,
            transport: &self.radial.transport,
            property: &self.radial.property,
            sorption: &self.sorption,
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum RadialLayerPath {
    /// A single material integral over the complete center-to-center path.
    WithinLayer { metric: f64 },
    /// Two original half paths plus an explicitly declared interface law.
    BetweenLayers { metrics: [f64; 2], area: f64 },
}

/// Contiguous original radial layers ordered from the inner to outer surface.
/// A solid center has zero area and always remains sealed. Cylinder lengths,
/// temperatures and accepted clocks agree exactly; gaps/overlaps are refused.
#[derive(Clone, Debug)]
pub struct MoistureLayeredRadial {
    layers: Vec<MoistureRadialLayer>,
    layer_for_cell: Vec<usize>,
    activity: Vec<f64>,
    activity_bounds: [f64; 2],
    paths: Vec<RadialLayerPath>,
    boundary_metrics: [f64; 2],
    boundary_areas: [f64; 2],
    exterior_transfer: [MoistureActivityTransfer; 2],
    interface_transfer: Vec<MoistureActivityTransfer>,
    identity: ContentHash,
}

impl MoistureLayeredRadial {
    /// Join admitted native shells without replacing their sources or storage.
    /// Ideal contact is the default. A finite coating/interface resistance must
    /// be supplied explicitly with [`Self::with_transfer`].
    pub fn new(
        layers: Vec<MoistureRadialLayer>,
        maximum_cells: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let Some(first) = layers.first() else {
            return Err(MoistureError::Invalid(
                "a nonempty radial material stack is required",
            ));
        };
        if maximum_cells == 0 || maximum_cells > MAX_MOISTURE_CELLS {
            return Err(MoistureError::Invalid(
                "bounded layered radial cell budget required",
            ));
        }
        let symmetry = first.radial.geometry.symmetry;
        let temperature = first.sorption.temperature_kelvin();
        let time = first.radial.time_s;
        let mut bounds: [f64; 2] = [0.0, 1.0];
        let mut activity = Vec::new();
        let mut layer_for_cell = Vec::new();
        let mut paths = Vec::new();
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.moisture-layered-radial.v1");
        for (k, layer) in layers.iter().enumerate() {
            poll(cx)?;
            let radial = &layer.radial;
            let n = radial.moisture.len();
            let density = layer.dry_density();
            if n > maximum_cells.saturating_sub(activity.len())
                || radial.geometry.symmetry != symmetry
                || layer.sorption.temperature_kelvin().to_bits() != temperature.to_bits()
                || radial.time_s.to_bits() != time.to_bits()
            {
                return Err(MoistureError::Invalid(
                    "radial layers need common symmetry, temperature and time within the total cell budget",
                ));
            }
            if k > 0 {
                let prior = &layers[k - 1];
                let face = radial.geometry.faces_m[0];
                if face.to_bits() != prior.radial.geometry.faces_m.last().unwrap().to_bits() {
                    return Err(MoistureError::Invalid(
                        "radial material interfaces must be contiguous without gaps or overlaps",
                    ));
                }
                let metrics = [
                    prior.dry_density()
                        / symmetry.resistance(*prior.radial.centers_m.last().unwrap(), face),
                    density / symmetry.resistance(face, radial.centers_m[0]),
                ];
                let area = symmetry.area(face);
                if metrics.iter().any(|&m| !positive(m)) || !positive(area) {
                    return Err(MoistureError::Invalid(
                        "radial interface geometry is unrepresentable",
                    ));
                }
                paths.push(RadialLayerPath::BetweenLayers { metrics, area });
            }
            let supported = layer.sorption.activity_bounds();
            bounds[0] = bounds[0].max(supported[0]);
            bounds[1] = bounds[1].min(supported[1]);
            hash.update(radial.identity.as_bytes());
            hash.update(layer.sorption.identity().as_bytes());
            hash.update(&(n as u64).to_le_bytes());
            for (i, &m) in radial.moisture.iter().enumerate() {
                poll(cx)?;
                activity.push(layer.sorption.activity(m)?);
                layer_for_cell.push(k);
                if i > 0 {
                    let metric =
                        density / symmetry.resistance(radial.centers_m[i - 1], radial.centers_m[i]);
                    if !positive(metric) {
                        return Err(MoistureError::Invalid(
                            "radial within-material path is unrepresentable",
                        ));
                    }
                    paths.push(RadialLayerPath::WithinLayer { metric });
                }
            }
        }
        if bounds[0] >= bounds[1] || activity.iter().any(|&a| a < bounds[0] || a > bounds[1]) {
            return Err(MoistureError::Invalid(
                "radial layers need common source activity support containing their initial state",
            ));
        }
        let last = layers.last().expect("nonempty admitted layers");
        let inner = first.radial.geometry.faces_m[0];
        let outer = *last.radial.geometry.faces_m.last().unwrap();
        let boundary_areas = [symmetry.area(inner), symmetry.area(outer)];
        let boundary_metrics = [
            if inner == 0.0 {
                0.0
            } else {
                first.dry_density() / symmetry.resistance(inner, first.radial.centers_m[0])
            },
            last.dry_density() / symmetry.resistance(*last.radial.centers_m.last().unwrap(), outer),
        ];
        if !positive(boundary_metrics[1])
            || !positive(boundary_areas[1])
            || (inner > 0.0 && (!positive(boundary_metrics[0]) || !positive(boundary_areas[0])))
        {
            return Err(MoistureError::Invalid(
                "radial exterior geometry is unrepresentable",
            ));
        }
        let interface_transfer = vec![MoistureActivityTransfer::IdealContact; layers.len() - 1];
        let result = Self {
            layers,
            layer_for_cell,
            activity,
            activity_bounds: bounds,
            paths,
            boundary_metrics,
            boundary_areas,
            exterior_transfer: [MoistureActivityTransfer::IdealContact; 2],
            interface_transfer,
            identity: hash.finalize(),
        };
        if !positive(result.dry_mass_kg()) || !positive(result.wet_mass_kg()) {
            return Err(MoistureError::Invalid(
                "layered radial inventory is unrepresentable",
            ));
        }
        Ok(result)
    }

    /// Explicit area-specific activity conductance at each physical surface and
    /// material join. Each law uses that face's own curved area; interface water
    /// has no separate storage and is debited/credited exactly once.
    pub fn with_transfer(
        mut self,
        exterior: [MoistureActivityTransfer; 2],
        interfaces: Vec<MoistureActivityTransfer>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        if interfaces.len() != self.layers.len() - 1 {
            return Err(MoistureError::Invalid(
                "one activity-transfer law is required for each radial material join",
            ));
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-couple.moisture-layered-radial-transfer.v1");
        hash.update(self.identity.as_bytes());
        let mut record = |transfer: MoistureActivityTransfer, area| -> Result<(), MoistureError> {
            match transfer.conductance(area)? {
                None => hash.update(&[0]),
                Some(g) => {
                    hash.update(&[1]);
                    hash.update(&g.to_bits().to_le_bytes());
                }
            }
            Ok(())
        };
        for (side, &transfer) in exterior.iter().enumerate() {
            poll(cx)?;
            // A center has no physical exchange area, regardless of a source.
            if self.boundary_areas[side] == 0.0
                && matches!(transfer, MoistureActivityTransfer::ConductanceKgM2S(_))
            {
                return Err(MoistureError::Invalid(
                    "finite activity transfer requires a nonzero radial surface area",
                ));
            }
            record(transfer, self.boundary_areas[side])?;
        }
        for (k, &transfer) in interfaces.iter().enumerate() {
            poll(cx)?;
            let radius = *self.layers[k].radial.geometry.faces_m.last().unwrap();
            record(
                transfer,
                self.layers[k].radial.geometry.symmetry.area(radius),
            )?;
        }
        self.exterior_transfer = exterior;
        self.interface_transfer = interfaces;
        self.identity = hash.finalize();
        Ok(self)
    }

    pub fn layers(&self) -> &[MoistureRadialLayer] {
        &self.layers
    }
    pub fn activity(&self) -> &[f64] {
        &self.activity
    }
    pub const fn activity_bounds(&self) -> [f64; 2] {
        self.activity_bounds
    }
    pub const fn exterior_transfer(&self) -> &[MoistureActivityTransfer; 2] {
        &self.exterior_transfer
    }
    pub fn interface_transfer(&self) -> &[MoistureActivityTransfer] {
        &self.interface_transfer
    }
    pub const fn definition_identity(&self) -> ContentHash {
        self.identity
    }
    pub fn time_s(&self) -> f64 {
        self.layers[0].radial.time_s
    }
    pub fn dry_mass_kg(&self) -> f64 {
        self.layers.iter().map(|l| l.radial.dry_mass_kg()).sum()
    }
    pub fn water_mass_kg(&self) -> f64 {
        self.layers.iter().map(|l| l.radial.water_mass_kg()).sum()
    }
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }

    /// One fully implicit conservative proposal using the same nonlinear
    /// activity/storage/interface solver as heterogeneous slabs. Rejection
    /// publishes no partial layer, face transfer or advanced clock.
    pub fn advance(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredReport), MoistureError> {
        poll(cx)?;
        if self.boundary_areas[0] == 0.0 && !matches!(faces[0], MoistureLayeredFace::Sealed) {
            return Err(MoistureError::Invalid(
                "a solid radial center has zero area and must remain sealed",
            ));
        }
        let reservoir = std::array::from_fn(|side| match faces[side] {
            _ if self.exterior_transfer[side] == MoistureActivityTransfer::Sealed => None,
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        });
        let original: Vec<_> = self
            .layers
            .iter()
            .flat_map(|l| l.radial.moisture.iter().copied())
            .collect();
        let dry_masses: Vec<_> = self
            .layers
            .iter()
            .flat_map(|l| l.radial.dry_masses_kg.iter().copied())
            .collect();
        let sources: Vec<_> = self
            .layers
            .iter()
            .map(MoistureRadialLayer::source)
            .collect();
        let candidate = kernel::advance(
            LayeredStorage {
                sources: &sources,
                layer_for_cell: &self.layer_for_cell,
                activity: &self.activity,
                original: &original,
                dry_masses: &dry_masses,
                activity_bounds: self.activity_bounds,
                time_s: self.time_s(),
                dry_mass_kg: self.dry_mass_kg(),
                water_mass_kg: self.water_mass_kg(),
            },
            dt,
            reservoir,
            step,
            nonlinear,
            |profiles, activity, reference, system, cx| {
                self.assemble_paths(dt, reservoir, profiles, activity, reference, system, cx)
            },
            cx,
        )?;
        let mut next = self.clone();
        next.activity = candidate.activity;
        let mut offset = 0;
        for layer in &mut next.layers {
            poll(cx)?;
            let n = layer.radial.moisture.len();
            layer
                .radial
                .moisture
                .copy_from_slice(&candidate.moisture[offset..offset + n]);
            layer.radial.time_s = candidate.report.interval_s[1];
            offset += n;
        }
        poll(cx)?;
        Ok((next, candidate.report))
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_paths(
        &self,
        dt: f64,
        reservoir: [Option<f64>; 2],
        profiles: &[KirchhoffMoistureDiffusivity],
        activity: &[f64],
        reference: f64,
        system: &mut LayeredSystem,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        let path = |layer: usize, metric: f64| ActivityPath {
            sorption: &self.layers[layer].sorption,
            diffusion: &profiles[layer],
            metric,
            pullback: 1.0,
        };
        for (i, face) in self.paths.iter().enumerate() {
            poll(cx)?;
            let left = self.layer_for_cell[i];
            let right = self.layer_for_cell[i + 1];
            let g = match *face {
                RadialLayerPath::WithinLayer { metric } => {
                    path(left, metric).conductance(activity[i], activity[i + 1])?
                }
                RadialLayerPath::BetweenLayers { metrics, area } => {
                    let (g, surfaces) = interface::transfer_paths(
                        [path(left, metrics[0]), path(right, metrics[1])],
                        [activity[i], activity[i + 1]],
                        self.interface_transfer[left].conductance(area)?,
                        cx,
                    )?;
                    system.interfaces.push(MoistureLayerInterfaceReport {
                        layers: [left, right],
                        relative_activity: surfaces[0],
                        right_relative_activity: surfaces[1],
                        moisture: [
                            self.layers[left].sorption.moisture(surfaces[0])?,
                            self.layers[right].sorption.moisture(surfaces[1])?,
                        ],
                        water_left_to_right_kg: dt * g * (activity[i] - activity[i + 1]),
                    });
                    g
                }
            };
            let coefficient = kirchhoff::scaled_diffusivity(dt, g)?;
            system.operator.off[i] = coefficient;
            system.operator.diagonal[i] += coefficient;
            system.operator.diagonal[i + 1] += coefficient;
        }
        for (side, cell) in [0, activity.len() - 1].into_iter().enumerate() {
            poll(cx)?;
            if let Some(a) = reservoir[side] {
                let layer = self.layer_for_cell[cell];
                let (g, surface) = interface::surface_path(
                    path(layer, self.boundary_metrics[side]),
                    [activity[cell], a],
                    self.exterior_transfer[side].conductance(self.boundary_areas[side])?,
                    cx,
                )?;
                let coefficient = kirchhoff::scaled_diffusivity(dt, g)?;
                system.boundary[side] = coefficient;
                system.surfaces[side] = Some(surface);
                system.operator.diagonal[cell] += coefficient;
                system.operator.storage_boundary[cell] += coefficient;
                system.rhs[cell] += coefficient * (a - reference);
            }
        }
        Ok(())
    }
}
