//! Isothermal multilayer transport with source-owned equilibrium partitioning.
//!
//! Each layer retains its original slab, dry skeleton, diffusivity and strictly
//! increasing sorption source M(a). At fixed temperature, a is the source's
//! dilute-vapor relative humidity/activity. Unlike materials share interface a,
//! not interface M. Their two half-cell Kirchhoff fluxes balance locally.
//!
//! Backward Euler evolves the actual water inventory. Positive source storage
//! secants and symmetric face conductances retain the existing tridiagonal
//! solver. Explicit surface and interface activity conductances may add transfer
//! resistance without an interface store. Ordinary steps hold reference geometry;
//! the native solid consumer can explicitly pull a spatial law back through
//! compatible axial swelling. Hysteresis, sorption heat, liquid flow, general
//! moving-mesh transport and an entropy certificate remain outside this reduction.

use super::*;
use fs_material::state_point::MonotoneMoistureSorption;

/// Original finite vapor inventory coupled to one or more heterogeneous stacks.
#[path = "layered/atmosphere.rs"]
pub mod atmosphere;
#[path = "layered/interface.rs"]
pub(super) mod interface;
#[path = "layered/kernel.rs"]
pub(super) mod kernel;
#[path = "layered/sdirk.rs"]
mod sdirk;
pub use sdirk::MoistureLayeredSdirk2Report;

use kernel::{LayeredSource, LayeredStorage, LayeredSystem};

/// Bind equilibrium to the original current transport context, including a
/// supported prescribed-temperature event, while retaining reference density.
pub(in crate::moisture) fn admit_sorption(
    card: &MaterialCard,
    transport: &ResolvedMaterialStatePoint,
    maximum_moisture: f64,
    moisture: &[f64],
    sorption: &MonotoneMoistureSorption,
    cx: &Cx<'_>,
) -> Result<(), MoistureError> {
    poll(cx)?;
    if sorption.card_identity() != card.content_hash()
        || sorption.moisture_bounds()[1] > maximum_moisture
    {
        return Err(MoistureError::Invalid(
            "layer sorption must use the same original material and moisture domain",
        ));
    }
    let point = sorption.source_point();
    let coordinates_match = point
        .axes()
        .iter()
        .filter(|(name, _)| name.as_str() != "relative_humidity")
        .map(|(name, value)| (name.as_str(), value.to_bits()))
        .eq(transport
            .query_point()
            .iter()
            .filter(|(name, _)| name.as_str() != "moisture_content")
            .map(|(name, value)| (name.as_str(), value.to_bits())));
    let quantities_match = point
        .axis_quantities()
        .iter()
        .filter(|(name, _)| name.as_str() != "relative_humidity")
        .eq(transport
            .axis_quantities()
            .iter()
            .filter(|(name, _)| name.as_str() != "moisture_content"));
    if !coordinates_match || !quantities_match {
        return Err(MoistureError::Invalid(
            "layer sorption and transport must hold identical external coordinates and quantities",
        ));
    }
    for &m in moisture {
        poll(cx)?;
        sorption.activity(m)?;
    }
    Ok(())
}

/// One original material/specimen owner and its compatible equilibrium source.
#[derive(Clone, Debug)]
pub struct MoistureLayer {
    slab: MoistureSlab,
    sorption: MonotoneMoistureSorption,
}

impl MoistureLayer {
    /// Associate an isotherm with this exact material and held source context.
    /// Its complete supported moisture range must fit the slab applicability.
    pub fn new(
        slab: MoistureSlab,
        sorption: MonotoneMoistureSorption,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        admit_sorption(
            &slab.card,
            &slab.transport,
            slab.geometry.maximum_moisture,
            &slab.moisture,
            &sorption,
            cx,
        )?;
        Ok(Self { slab, sorption })
    }

    /// Accepted native slab, directly consumable by its original solid/plate owners.
    pub const fn slab(&self) -> &MoistureSlab {
        &self.slab
    }
    pub const fn sorption(&self) -> &MonotoneMoistureSorption {
        &self.sorption
    }
    fn spacing(&self) -> f64 {
        self.slab.geometry.length_m / self.slab.moisture.len() as f64
    }
}

/// The two exterior ends, ordered left/right. Prescribed RH is an infinite
/// reservoir at the layers' common fixed temperature; it is never an interior M.
#[derive(Clone, Copy, Debug)]
pub enum MoistureLayeredFace {
    Sealed,
    PrescribedRelativeHumidity(f64),
}

/// Transfer law per physical area at a surface or between two materials.
/// The finite law is j_water = h * (a_left - a_right), in kg/(m² s).
/// Its coefficient is explicit; a gas-side velocity or dry-basis moisture
/// coefficient cannot be substituted without the corresponding conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MoistureActivityTransfer {
    IdealContact,
    Sealed,
    ConductanceKgM2S(f64),
}

impl MoistureActivityTransfer {
    pub(in crate::moisture) fn conductance(self, area: f64) -> Result<Option<f64>, MoistureError> {
        match self {
            Self::IdealContact => Ok(None),
            Self::Sealed => Ok(Some(0.0)),
            Self::ConductanceKgM2S(h) if positive(h) => {
                Ok(Some(kirchhoff::scaled_diffusivity(area, h)?))
            }
            Self::ConductanceKgM2S(_) => Err(MoistureError::Invalid(
                "activity transfer needs a finite positive area-specific conductance",
            )),
        }
    }
}

/// Ordered material layers with a common area, temperature and accepted time.
#[derive(Clone, Debug)]
pub struct MoistureLayeredSlab {
    layers: Vec<MoistureLayer>,
    layer_for_cell: Vec<usize>,
    activity: Vec<f64>,
    activity_bounds: [f64; 2],
    exterior_transfer: [MoistureActivityTransfer; 2],
    interface_transfer: Vec<MoistureActivityTransfer>,
    identity: ContentHash,
}

/// Pullback of an explicitly spatial transport law onto the original dry
/// skeleton. Storage is unchanged. Only the solid consumer may construct this
/// metric, after checking that its deformation preserves the axial reduction.
pub(super) struct MoistureLayeredMetric {
    pub(super) inverse_axial_stretch_squared: Vec<f64>,
    pub(super) face_area_ratio: f64,
}

type MetricResolver<'a> =
    dyn FnMut(&MoistureLayeredSlab) -> Result<MoistureLayeredMetric, MoistureError> + 'a;

/// Two material-specific surfaces with equal/opposite water transfer.
#[derive(Clone, Debug)]
pub struct MoistureLayerInterfaceReport {
    pub layers: [usize; 2],
    /// Left material's activity; also the common activity at ideal contact.
    pub relative_activity: f64,
    /// Right material's activity. A finite resistance permits a physical jump.
    pub right_relative_activity: f64,
    pub moisture: [f64; 2],
    /// Positive from the left layer into the right; debited/credited once.
    pub water_left_to_right_kg: f64,
}

/// Physical inventory and original source/solver evidence for one whole stack.
#[derive(Clone, Debug)]
pub struct MoistureLayeredReport {
    pub interval_s: [f64; 2],
    pub boundary_water_kg: [f64; 2],
    /// Actual solid surface activities, which may differ from reservoir RH.
    pub surface_activities: [Option<f64>; 2],
    pub interfaces: Vec<MoistureLayerInterfaceReport>,
    pub stored_water_kg: [f64; 2],
    pub water_balance_residual_kg: f64,
    pub diffusivity_sources: Vec<ContentHash>,
    pub sorption_sources: Vec<ContentHash>,
    pub transport_end: Vec<ResolvedMaterialStatePoint>,
    pub linear: SolveReport,
    pub nonlinear: MoistureNonlinearReport,
}

impl MoistureLayeredSlab {
    /// Join actual slabs end to end. Layer order is physical; no material
    /// averaging or continuity of dry-basis moisture is inferred at a join.
    pub fn new(
        layers: Vec<MoistureLayer>,
        maximum_cells: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let Some(first) = layers.first() else {
            return Err(MoistureError::Invalid(
                "a nonempty moisture stack is required",
            ));
        };
        if maximum_cells == 0 || maximum_cells > MAX_MOISTURE_CELLS {
            return Err(MoistureError::Invalid(
                "bounded layered cell budget required",
            ));
        }
        let area = first.slab.geometry.cross_section_m2;
        let temperature = first.sorption.temperature_kelvin();
        let time = first.slab.time_s;
        let mut bounds: [f64; 2] = [0.0, 1.0];
        let mut activity = Vec::new();
        let mut layer_for_cell = Vec::new();
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.moisture-layered.v1");
        for (k, layer) in layers.iter().enumerate() {
            poll(cx)?;
            let n = layer.slab.moisture.len();
            if n > maximum_cells.saturating_sub(activity.len())
                || layer.slab.geometry.cross_section_m2.to_bits() != area.to_bits()
                || layer.sorption.temperature_kelvin().to_bits() != temperature.to_bits()
                || layer.slab.time_s.to_bits() != time.to_bits()
            {
                return Err(MoistureError::Invalid(
                    "layers need the same area, temperature and time within the total cell budget",
                ));
            }
            let supported = layer.sorption.activity_bounds();
            bounds[0] = bounds[0].max(supported[0]);
            bounds[1] = bounds[1].min(supported[1]);
            hash.update(layer.slab.definition_identity().as_bytes());
            hash.update(layer.sorption.identity().as_bytes());
            hash.update(&(n as u64).to_le_bytes());
            for &m in &layer.slab.moisture {
                poll(cx)?;
                activity.push(layer.sorption.activity(m)?);
                layer_for_cell.push(k);
            }
        }
        if bounds[0] >= bounds[1] || activity.iter().any(|&a| a < bounds[0] || a > bounds[1]) {
            return Err(MoistureError::Invalid(
                "all layers need a common supported activity interval containing their initial state",
            ));
        }
        let interface_transfer = vec![MoistureActivityTransfer::IdealContact; layers.len() - 1];
        let stack = Self {
            layers,
            layer_for_cell,
            activity,
            activity_bounds: bounds,
            exterior_transfer: [MoistureActivityTransfer::IdealContact; 2],
            interface_transfer,
            identity: hash.finalize(),
        };
        if !positive(stack.dry_mass_kg()) || !positive(stack.wet_mass_kg()) {
            return Err(MoistureError::Invalid(
                "layered physical inventory is unrepresentable",
            ));
        }
        Ok(stack)
    }

    /// Pin explicit transfer laws to the original geometry. The existing
    /// constructor selects ideal contact; all exposed finite-atmosphere ends
    /// automatically use these same exterior laws. There is no interface water
    /// store and changing a law does not alter the accepted physical inventory.
    pub fn with_transfer(
        mut self,
        exterior: [MoistureActivityTransfer; 2],
        interfaces: Vec<MoistureActivityTransfer>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        if interfaces.len() != self.layers.len() - 1 {
            return Err(MoistureError::Invalid(
                "one activity-transfer law is required for each material join",
            ));
        }
        let area = self.layers[0].slab.geometry.cross_section_m2;
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.moisture-layered-transfer.v1");
        hash.update(self.identity.as_bytes());
        for transfer in exterior.iter().chain(&interfaces) {
            poll(cx)?;
            let g = transfer.conductance(area)?;
            match g {
                None => hash.update(&[0]),
                Some(g) => {
                    hash.update(&[1]);
                    hash.update(&g.to_bits().to_le_bytes());
                }
            }
        }
        self.exterior_transfer = exterior;
        self.interface_transfer = interfaces;
        self.identity = hash.finalize();
        Ok(self)
    }

    pub const fn exterior_transfer(&self) -> &[MoistureActivityTransfer; 2] {
        &self.exterior_transfer
    }
    pub fn interface_transfer(&self) -> &[MoistureActivityTransfer] {
        &self.interface_transfer
    }

    pub fn layers(&self) -> &[MoistureLayer] {
        &self.layers
    }
    pub fn activity(&self) -> &[f64] {
        &self.activity
    }
    pub const fn activity_bounds(&self) -> [f64; 2] {
        self.activity_bounds
    }
    pub const fn definition_identity(&self) -> ContentHash {
        self.identity
    }
    pub fn time_s(&self) -> f64 {
        self.layers[0].slab.time_s
    }
    pub fn dry_mass_kg(&self) -> f64 {
        self.layers.iter().map(|l| l.slab.dry_mass_kg()).sum()
    }
    pub fn water_mass_kg(&self) -> f64 {
        self.layers.iter().map(|l| l.slab.water_mass_kg()).sum()
    }
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }

    /// Conservative fully implicit proposal. The original nonlinear and linear
    /// work budgets never expand; source or solver refusal publishes no layer.
    pub fn advance(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        controls: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredReport), MoistureError> {
        self.advance_with_metric(dt, faces, step, controls, None, cx)
    }

    /// The original nonlinear endpoint also owns geometry feedback. A supplied
    /// resolver must evaluate the actual proposed state: its metric enters the
    /// endpoint residual and physical face accounts before acceptance.
    pub(super) fn advance_with_metric(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        controls: MoistureNonlinearConfig,
        mut resolve_metric: Option<&mut MetricResolver<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredReport), MoistureError> {
        poll(cx)?;
        let end = self.time_s() + dt;
        let reservoir = std::array::from_fn(|side| match faces[side] {
            _ if self.exterior_transfer[side] == MoistureActivityTransfer::Sealed => None,
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        });
        let original: Vec<f64> = self
            .layers
            .iter()
            .flat_map(|l| l.slab.moisture.iter().copied())
            .collect();
        let dry_masses: Vec<_> = self
            .layer_for_cell
            .iter()
            .map(|&layer| self.layers[layer].slab.dry_mass_per_cell_kg)
            .collect();
        let sources: Vec<_> = self
            .layers
            .iter()
            .map(|layer| LayeredSource {
                card: &layer.slab.card,
                transport: &layer.slab.transport,
                property: &layer.slab.diffusivity_property,
                sorption: &layer.sorption,
            })
            .collect();
        let state = LayeredStorage {
            sources: &sources,
            layer_for_cell: &self.layer_for_cell,
            activity: &self.activity,
            original: &original,
            dry_masses: &dry_masses,
            activity_bounds: self.activity_bounds,
            time_s: self.time_s(),
            dry_mass_kg: self.dry_mass_kg(),
            water_mass_kg: self.water_mass_kg(),
        };
        let mut first_metric = true;
        let mut proposed = None;
        let candidate = kernel::advance(
            state,
            dt,
            reservoir,
            step,
            controls,
            |profiles, activity, reference, system, cx| {
                let metric = if let Some(resolve) = &mut resolve_metric {
                    if first_metric {
                        first_metric = false;
                        Some(resolve(self)?)
                    } else {
                        let next = self.with_activity(activity, end, cx)?;
                        let metric = resolve(&next)?;
                        proposed = Some(next);
                        Some(metric)
                    }
                } else {
                    None
                };
                self.assemble_paths(
                    dt,
                    reservoir,
                    profiles,
                    activity,
                    reference,
                    metric.as_ref(),
                    system,
                    cx,
                )
            },
            cx,
        )?;
        let next = match proposed {
            Some(state) => state,
            None => self.with_profile(&candidate.activity, &candidate.moisture, end, cx)?,
        };
        poll(cx)?;
        Ok((next, candidate.report))
    }

    fn with_activity(
        &self,
        activity: &[f64],
        time_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let mut moisture = Vec::with_capacity(activity.len());
        for (&a, &layer) in activity.iter().zip(&self.layer_for_cell) {
            poll(cx)?;
            moisture.push(self.layers[layer].sorption.moisture(a)?);
        }
        self.with_profile(activity, &moisture, time_s, cx)
    }

    fn with_profile(
        &self,
        activity: &[f64],
        moisture: &[f64],
        time_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let mut next = self.clone();
        next.activity.copy_from_slice(activity);
        let mut offset = 0;
        for layer in &mut next.layers {
            poll(cx)?;
            let n = layer.slab.moisture.len();
            layer
                .slab
                .moisture
                .copy_from_slice(&moisture[offset..offset + n]);
            layer.slab.time_s = time_s;
            offset += n;
        }
        Ok(next)
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_paths(
        &self,
        dt: f64,
        reservoir: [Option<f64>; 2],
        profiles: &[KirchhoffMoistureDiffusivity],
        activity: &[f64],
        reference: f64,
        metric: Option<&MoistureLayeredMetric>,
        system: &mut LayeredSystem,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        let n = activity.len();
        if metric.is_some_and(|m| {
            m.inverse_axial_stretch_squared.len() != n
                || !positive(m.face_area_ratio)
                || m.inverse_axial_stretch_squared
                    .iter()
                    .any(|&x| !positive(x))
        }) {
            return Err(MoistureError::Invalid("invalid layered deformation metric"));
        }
        for i in 0..n - 1 {
            poll(cx)?;
            let left = self.layer_for_cell[i];
            let right = self.layer_for_cell[i + 1];
            let pullback = metric.map_or([1.0; 2], |m| {
                [
                    m.inverse_axial_stretch_squared[i],
                    m.inverse_axial_stretch_squared[i + 1],
                ]
            });
            let conductance = if left == right && pullback[0] == pullback[1] {
                // One homogeneous full-cell path; no nonlinear interface root.
                kirchhoff::scaled_diffusivity(
                    0.5,
                    interface::half_conductance(
                        &self.layers[left],
                        &profiles[left],
                        activity[i],
                        activity[i + 1],
                        pullback[0],
                    )?,
                )?
            } else {
                let (g, a) = interface::transfer_join(
                    [&self.layers[left], &self.layers[right]],
                    [&profiles[left], &profiles[right]],
                    [activity[i], activity[i + 1]],
                    if left == right {
                        None
                    } else {
                        self.interface_transfer[left].conductance(
                            self.layers[left].slab.geometry.cross_section_m2
                                * metric.map_or(1.0, |m| m.face_area_ratio),
                        )?
                    },
                    pullback,
                    cx,
                )?;
                if left != right {
                    system.interfaces.push(MoistureLayerInterfaceReport {
                        layers: [left, right],
                        relative_activity: a[0],
                        right_relative_activity: a[1],
                        moisture: [
                            self.layers[left].sorption.moisture(a[0])?,
                            self.layers[right].sorption.moisture(a[1])?,
                        ],
                        water_left_to_right_kg: dt * g * (activity[i] - activity[i + 1]),
                    });
                }
                g
            };
            let coefficient = kirchhoff::scaled_diffusivity(dt, conductance)?;
            system.operator.off[i] = coefficient;
            system.operator.diagonal[i] += coefficient;
            system.operator.diagonal[i + 1] += coefficient;
        }
        for (side, cell) in [0, n - 1].into_iter().enumerate() {
            poll(cx)?;
            if let Some(a) = reservoir[side] {
                let k = self.layer_for_cell[cell];
                let (g, surface) = interface::surface(
                    &self.layers[k],
                    &profiles[k],
                    [activity[cell], a],
                    self.exterior_transfer[side].conductance(
                        self.layers[k].slab.geometry.cross_section_m2
                            * metric.map_or(1.0, |m| m.face_area_ratio),
                    )?,
                    metric.map_or(1.0, |m| m.inverse_axial_stretch_squared[cell]),
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
