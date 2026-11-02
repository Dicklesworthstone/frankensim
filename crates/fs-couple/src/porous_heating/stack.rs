//! One native thermal advance for every material in a coherent porous stack.

use super::*;
use fs_duct::porous::PorousStackExcitation;
use std::ops::Range;

/// A complete heterogeneous acoustic specimen bound to its native heat cells.
///
/// Every constituent receives the coherent drive reconstructed by the original
/// scattering owner. Its local analytic viscous and thermal work is integrated
/// over the authored cell volumes. Only external incident/outgoing power enters
/// the stack energy account; waves exchanged between layers are internal.
/// Thermal links across material interfaces retain their original conductances.
///
/// This has the same fixed-drive, rigid-frame, thermally lumped constituent
/// assumptions as [`PorousHeating`]. It does not infer fluid/skeleton heat
/// partition, temperature-dependent acoustics, or interface contact resistance.
pub struct PorousStackHeating<'e, 's, 'l, 'm> {
    excitation: &'e PorousStackExcitation<'s, 'l, 'm>,
    heat: PreparedHeat,
    layer_cells: Box<[Range<usize>]>,
}

impl<'e, 's, 'l, 'm> PorousStackHeating<'e, 's, 'l, 'm> {
    /// Bind one complete local-depth partition per physical acoustic layer.
    ///
    /// `partitions[i]` covers exactly `excitation.stack().layers()[i]`, using
    /// depths measured from that layer's left face. Thermal cells may use
    /// different material/caloric laws and masses in different layers. They
    /// must be distinct, adjoining uniform prisms on one straight specimen
    /// axis, with the original full area and one common body identity.
    ///
    /// The cell budget covers the complete specimen, not each layer separately.
    /// Preparation integrates each cell once and merges cached native sources;
    /// a thermal substep never resolves acoustics or samples an interior field.
    pub fn try_new(
        excitation: &'e PorousStackExcitation<'s, 'l, 'm>,
        thermal: &DepositionSimulation<'_>,
        partitions: &[&[PorousHeatCell]],
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
        let layers = excitation.layer_excitations();
        if partitions.len() != layers.len() || partitions.is_empty() {
            return Err(PorousHeatingError::Invalid(
                "porous stack heat needs one complete partition per acoustic layer",
            ));
        }
        let count = partitions.iter().try_fold(0usize, |sum, cells| {
            sum.checked_add(cells.len())
                .ok_or(DepositionError::Budget("porous heating cells"))
        })?;
        if maximum_cells == 0
            || maximum_cells > MAX_POROUS_HEATING_CELLS
            || count == 0
            || count > maximum_cells
        {
            return Err(DepositionError::Budget("porous heating cells").into());
        }
        let mut cells: Vec<PorousCellHeating> = Vec::with_capacity(count);
        let mut sources_w = BTreeMap::new();
        let mut layer_cells = Vec::with_capacity(layers.len());
        let (mut viscous, mut relaxation, mut absorbed) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        let mut axis = [0.0; 3];
        for (layer, partition) in layers.iter().zip(partitions) {
            checkpoint()?;
            let local =
                PorousHeating::try_new(layer, thermal, partition, maximum_cells, &mut checkpoint)?;
            let first = &local.heat.cells[0];
            let delta =
                std::array::from_fn(|i| first.geometry.end_m[i] - first.geometry.start_m[i]);
            if let Some(previous) = cells.last() {
                if previous.geometry.end_m != first.geometry.start_m
                    || previous.body_id != first.body_id
                    || norm(std::array::from_fn(|i| {
                        delta[i] - axis[i] * first.end_depth_m
                    })) > 256.0 * f64::EPSILON * first.end_depth_m
                {
                    return Err(PorousHeatingError::Invalid(
                        "porous stack heat layers must adjoin along one body normal",
                    ));
                }
            } else {
                let length = norm(delta);
                axis = delta.map(|value| value / length);
            }
            for (cell, power) in local.heat.sources_w {
                checkpoint()?;
                if sources_w.insert(cell, power).is_some() {
                    return Err(PorousHeatingError::Invalid(
                        "porous stack heat cells cannot belong to more than one layer",
                    ));
                }
            }
            viscous = viscous + Dd::from_f64(local.heat.dissipation.viscous_power_w);
            relaxation = relaxation + Dd::from_f64(local.heat.dissipation.thermal_power_w);
            absorbed = absorbed + Dd::from_f64(local.heat.dissipation.absorbed_power_w);
            let start = cells.len();
            cells.extend(local.heat.cells);
            layer_cells.push(start..cells.len());
        }
        let dissipation = PorousLayerDissipation {
            viscous_power_w: finite(viscous.to_f64())?,
            thermal_power_w: finite(relaxation.to_f64())?,
            absorbed_power_w: finite(absorbed.to_f64())?,
        };
        let whole = excitation.dissipation();
        close_power(
            whole.absorbed_power_w,
            dissipation.absorbed_power_w,
            whole.absorbed_power_w.max(dissipation.absorbed_power_w),
            count + layers.len(),
        )?;
        let boundary = excitation.source_power();
        let power_residual_w = close_power(
            boundary.incident_power_w - boundary.outgoing_power_w - boundary.load_absorbed_power_w,
            dissipation.absorbed_power_w,
            boundary
                .incident_power_w
                .max(boundary.outgoing_power_w)
                .max(dissipation.absorbed_power_w),
            count + layers.len(),
        )?;
        checkpoint()?;
        Ok(Self {
            excitation,
            heat: PreparedHeat {
                cells,
                sources_w,
                dissipation,
                power_residual_w,
                boundary_power_w: [boundary.incident_power_w, boundary.outgoing_power_w],
                load_absorbed_power_w: boundary.load_absorbed_power_w,
            },
            layer_cells: layer_cells.into_boxed_slice(),
        })
    }

    /// Original prepared whole-stack field; no isolated-slab approximation.
    pub const fn excitation(&self) -> &'e PorousStackExcitation<'s, 'l, 'm> {
        self.excitation
    }

    /// All cells in physical layer/depth order. Each depth is LOCAL to its layer.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }

    /// The original cell partition of an authored acoustic layer, if present.
    pub fn layer_cells(&self, layer: usize) -> Option<&[PorousCellHeating]> {
        self.layer_cells
            .get(layer)
            .map(|range| &self.heat.cells[range.clone()])
    }

    /// Positive analytic work summed across every physical constituent.
    pub const fn dissipation(&self) -> PorousLayerDissipation {
        self.heat.dissipation
    }

    /// External incident minus outgoing minus deposited acoustic power [W].
    pub const fn power_residual_w(&self) -> f64 {
        self.heat.power_residual_w
    }

    /// Advance all receiving materials and their original interfacial thermal
    /// links together. Internal acoustic exchange never becomes an extra source;
    /// native thermal sources, latent storage, events and rollback are unchanged.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousHeatingStep<'t>, PorousHeatingError> {
        self.heat.advance_to(initial, target_s, config, checkpoint)
    }
}
