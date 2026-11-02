//! Isothermal multilayer transport with source-owned equilibrium partitioning.
//!
//! Each layer retains its original slab, dry skeleton, diffusivity and strictly
//! increasing sorption source M(a). At fixed temperature, a is the source's
//! dilute-vapor relative humidity/activity. Unlike materials share interface a,
//! not interface M. Their two half-cell Kirchhoff fluxes balance locally.
//!
//! Backward Euler evolves the actual water inventory. Positive source storage
//! secants and symmetric face conductances retain the existing tridiagonal
//! solver. This ideal-contact reduction has no interface store or resistance,
//! hysteresis, sorption heat, liquid flow, moving geometry or entropy certificate.

use super::*;
use fs_la::tridiagonal::TridiagonalWorkspace;
use fs_material::state_point::MonotoneMoistureSorption;

#[path = "layered/interface.rs"]
mod interface;

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
        poll(cx)?;
        if sorption.card_identity() != slab.card.content_hash()
            || sorption.moisture_bounds()[1] > slab.geometry.maximum_moisture
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
            .eq(slab
                .reference
                .query_point()
                .iter()
                .filter(|(name, _)| name.as_str() != "moisture_content")
                .map(|(name, value)| (name.as_str(), value.to_bits())));
        let quantities_match = point
            .axis_quantities()
            .iter()
            .filter(|(name, _)| name.as_str() != "relative_humidity")
            .eq(slab
                .reference
                .axis_quantities()
                .iter()
                .filter(|(name, _)| name.as_str() != "moisture_content"));
        if !coordinates_match || !quantities_match {
            return Err(MoistureError::Invalid(
                "layer sorption and transport must hold identical external coordinates and quantities",
            ));
        }
        for &m in &slab.moisture {
            poll(cx)?;
            sorption.activity(m)?;
        }
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

/// Ordered material layers with a common area, temperature and accepted time.
#[derive(Clone, Debug)]
pub struct MoistureLayeredSlab {
    layers: Vec<MoistureLayer>,
    layer_for_cell: Vec<usize>,
    activity: Vec<f64>,
    activity_bounds: [f64; 2],
    identity: ContentHash,
}

/// Ideal-contact interface: one activity, two material-specific water contents.
#[derive(Clone, Debug)]
pub struct MoistureLayerInterfaceReport {
    pub layers: [usize; 2],
    pub relative_activity: f64,
    pub moisture: [f64; 2],
    /// Positive from the left layer into the right; debited/credited once.
    pub water_left_to_right_kg: f64,
}

/// Physical inventory and original source/solver evidence for one whole stack.
#[derive(Clone, Debug)]
pub struct MoistureLayeredReport {
    pub interval_s: [f64; 2],
    pub boundary_water_kg: [f64; 2],
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
        let stack = Self {
            layers,
            layer_for_cell,
            activity,
            activity_bounds: bounds,
            identity: hash.finalize(),
        };
        if !positive(stack.dry_mass_kg()) || !positive(stack.wet_mass_kg()) {
            return Err(MoistureError::Invalid(
                "layered physical inventory is unrepresentable",
            ));
        }
        Ok(stack)
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
        poll(cx)?;
        controls.validate()?;
        let end = self.time_s() + dt;
        if !positive(dt)
            || !end.is_finite()
            || end <= self.time_s()
            || !positive(step.water_balance_tolerance_kg)
            || !positive(step.maximum_content_increment)
        {
            return Err(MoistureError::Invalid(
                "representable forward layered interval and positive limits required",
            ));
        }
        let reservoir = faces.map(|face| match face {
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        });
        if reservoir
            .into_iter()
            .flatten()
            .any(|a| !a.is_finite() || a < self.activity_bounds[0] || a > self.activity_bounds[1])
        {
            return Err(MoistureError::Invalid(
                "boundary activity lies outside the layers' common source support",
            ));
        }
        let mut span = [f64::INFINITY, f64::NEG_INFINITY];
        for (i, a) in self
            .activity
            .iter()
            .copied()
            .chain(reservoir.into_iter().flatten())
            .enumerate()
        {
            if i % 128 == 0 {
                poll(cx)?;
            }
            span[0] = span[0].min(a);
            span[1] = span[1].max(a);
        }
        let mut profiles = Vec::with_capacity(self.layers.len());
        for layer in &self.layers {
            poll(cx)?;
            profiles.push(KirchhoffMoistureDiffusivity::from_binding(
                &layer.slab.card,
                &layer.slab.transport,
                &layer.slab.diffusivity_property,
                [
                    layer.sorption.moisture(span[0])?,
                    layer.sorption.moisture(span[1])?,
                ],
            )?);
        }
        // Solve deviations from a common activity, so even arbitrarily stiff
        // exact equilibrium has a zero RHS instead of cancelling large fluxes.
        let reference = span[0];
        let original: Vec<f64> = self
            .layers
            .iter()
            .flat_map(|l| l.slab.moisture.iter().copied())
            .collect();
        let mut iterate = self.activity.clone();
        let mut system =
            self.system(dt, reservoir, &profiles, &original, &iterate, reference, cx)?;
        let mut workspace =
            TridiagonalWorkspace::new(self.activity.len()).map_err(MoistureError::Tridiagonal)?;
        let mut solves = 0usize;
        for iteration in 1..=controls.maximum_iterations {
            poll(cx)?;
            let solution =
                system
                    .operator
                    .solve_tridiagonal(&system.rhs, step.linear, &mut workspace, cx)?;
            solves = solves
                .checked_add(solution.report.iters)
                .ok_or(MoistureError::Invalid("layered solve work count overflow"))?;
            if !solution.report.converged_euclidean() {
                return Err(MoistureError::NotConverged(solution.report));
            }
            let candidate: Vec<f64> = solution.x.iter().map(|delta| reference + delta).collect();
            let mut moisture = Vec::with_capacity(candidate.len());
            let mut transport = Vec::with_capacity(candidate.len());
            let mut water_change = 0.0;
            let roundoff = 128.0 * f64::EPSILON;
            for (i, &a) in candidate.iter().enumerate() {
                poll(cx)?;
                if !a.is_finite() || a < span[0] - roundoff || a > span[1] + roundoff {
                    return Err(MoistureError::Invalid(
                        "layered candidate violates the activity maximum principle",
                    ));
                }
                let layer = &self.layers[self.layer_for_cell[i]];
                let m = layer.sorption.moisture(a)?;
                let source = layer
                    .slab
                    .transport
                    .at_query_point(&layer.slab.card, &layer.slab.point_at(m)?)?;
                water_change += layer.slab.dry_mass_per_cell_kg
                    * (layer.sorption.moisture_change(self.activity[i], a)?
                        + (layer.sorption.moisture(self.activity[i])? - original[i]));
                moisture.push(m);
                transport.push(source);
            }
            let endpoint = self.system(
                dt, reservoir, &profiles, &original, &candidate, reference, cx,
            )?;
            let (_, relative) = endpoint.operator.defect(&solution.x, &endpoint.rhs, cx)?;
            let boundary = std::array::from_fn(|side| {
                let cell = if side == 0 { 0 } else { candidate.len() - 1 };
                reservoir[side].map_or(0.0, |a| endpoint.boundary[side] * (a - candidate[cell]))
            });
            let water_after = moisture
                .iter()
                .zip(&self.layer_for_cell)
                .map(|(m, &k)| m * self.layers[k].slab.dry_mass_per_cell_kg)
                .sum::<f64>();
            let received = boundary[0] + boundary[1];
            let balance = water_change - received;
            let stored_balance = water_after - self.water_mass_kg() - received;
            let nonlinear = MoistureNonlinearReport {
                iterations: iteration,
                tridiagonal_solves: solves,
                relative_residual: relative,
                water_balance_residual_kg: balance,
            };
            if relative < controls.relative_tolerance
                && balance.is_finite()
                && balance.abs() <= step.water_balance_tolerance_kg
                && stored_balance.is_finite()
                && stored_balance.abs() <= step.water_balance_tolerance_kg
            {
                for (cell, (&m, &old)) in moisture.iter().zip(&original).enumerate() {
                    let increment = (m - old).abs();
                    if increment > step.maximum_content_increment {
                        return Err(MoistureError::ContentIncrement {
                            cell,
                            increment,
                            maximum: step.maximum_content_increment,
                        });
                    }
                }
                if !positive(self.dry_mass_kg() + water_after) {
                    return Err(MoistureError::Invalid(
                        "layered endpoint wet inventory is unrepresentable",
                    ));
                }
                let report = MoistureLayeredReport {
                    interval_s: [self.time_s(), end],
                    boundary_water_kg: boundary,
                    interfaces: endpoint.interfaces,
                    stored_water_kg: [self.water_mass_kg(), water_after],
                    water_balance_residual_kg: balance,
                    diffusivity_sources: profiles
                        .iter()
                        .map(KirchhoffMoistureDiffusivity::identity)
                        .collect(),
                    sorption_sources: self.layers.iter().map(|l| l.sorption.identity()).collect(),
                    transport_end: transport,
                    linear: solution.report,
                    nonlinear,
                };
                let mut next = self.clone();
                next.activity = candidate;
                let mut offset = 0;
                for layer in &mut next.layers {
                    poll(cx)?;
                    let n = layer.slab.moisture.len();
                    layer
                        .slab
                        .moisture
                        .copy_from_slice(&moisture[offset..offset + n]);
                    layer.slab.time_s = end;
                    offset += n;
                }
                poll(cx)?;
                return Ok((next, report));
            }
            if iteration == controls.maximum_iterations {
                return Err(MoistureError::NonlinearNotConverged(nonlinear));
            }
            if controls.relaxation == 1.0 {
                iterate = candidate;
                system = endpoint;
            } else {
                for (a, candidate) in iterate.iter_mut().zip(candidate) {
                    *a = controls.relaxation.mul_add(candidate - *a, *a);
                }
                system =
                    self.system(dt, reservoir, &profiles, &original, &iterate, reference, cx)?;
            }
        }
        unreachable!("validated positive layered iteration budget")
    }

    #[allow(clippy::too_many_arguments)]
    fn system(
        &self,
        dt: f64,
        reservoir: [Option<f64>; 2],
        profiles: &[KirchhoffMoistureDiffusivity],
        original: &[f64],
        activity: &[f64],
        reference: f64,
        cx: &Cx<'_>,
    ) -> Result<LayeredSystem, MoistureError> {
        let n = activity.len();
        let mut system = LayeredSystem {
            operator: DiffusionStep {
                diagonal: Vec::with_capacity(n),
                off: vec![0.0; n - 1],
            },
            rhs: Vec::with_capacity(n),
            boundary: [0.0; 2],
            interfaces: Vec::new(),
        };
        for i in 0..n {
            poll(cx)?;
            let layer = &self.layers[self.layer_for_cell[i]];
            let mass = layer.slab.dry_mass_per_cell_kg;
            let capacity = mass
                * layer
                    .sorption
                    .mean_capacity(self.activity[i], activity[i])?;
            if !positive(capacity) {
                return Err(MoistureError::Invalid(
                    "layered water capacity is unrepresentable",
                ));
            }
            let offset = original[i] - layer.sorption.moisture(self.activity[i])?;
            system.operator.diagonal.push(capacity);
            system
                .rhs
                .push(capacity.mul_add(self.activity[i] - reference, mass * offset));
        }
        for i in 0..n - 1 {
            poll(cx)?;
            let left = self.layer_for_cell[i];
            let right = self.layer_for_cell[i + 1];
            let conductance = if left == right {
                // One homogeneous full-cell path; no nonlinear interface root.
                kirchhoff::scaled_diffusivity(
                    0.5,
                    interface::half_conductance(
                        &self.layers[left],
                        &profiles[left],
                        activity[i],
                        activity[i + 1],
                    )?,
                )?
            } else {
                let (g, a) = interface::join(
                    [&self.layers[left], &self.layers[right]],
                    [&profiles[left], &profiles[right]],
                    [activity[i], activity[i + 1]],
                    cx,
                )?;
                system.interfaces.push(MoistureLayerInterfaceReport {
                    layers: [left, right],
                    relative_activity: a,
                    moisture: [
                        self.layers[left].sorption.moisture(a)?,
                        self.layers[right].sorption.moisture(a)?,
                    ],
                    water_left_to_right_kg: dt * g * (activity[i] - activity[i + 1]),
                });
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
                let g =
                    interface::half_conductance(&self.layers[k], &profiles[k], activity[cell], a)?;
                let coefficient = kirchhoff::scaled_diffusivity(dt, g)?;
                system.boundary[side] = coefficient;
                system.operator.diagonal[cell] += coefficient;
                system.rhs[cell] += coefficient * (a - reference);
            }
        }
        if system.operator.diagonal.iter().any(|&x| !positive(x))
            || system.rhs.iter().any(|x| !x.is_finite())
            || system
                .interfaces
                .iter()
                .any(|face| !face.water_left_to_right_kg.is_finite())
        {
            return Err(MoistureError::Invalid("layered diffusion system overflow"));
        }
        Ok(system)
    }
}

struct LayeredSystem {
    operator: DiffusionStep,
    rhs: Vec<f64>,
    boundary: [f64; 2],
    interfaces: Vec<MoistureLayerInterfaceReport>,
}
