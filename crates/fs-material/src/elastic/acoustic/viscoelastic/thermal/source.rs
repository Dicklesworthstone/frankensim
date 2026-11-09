//! Preparation shared by the two original acoustic boundary owners.
//!
//! Native drives retain every source and compute positive interval work.
//! Thermal partitioning changes observation volumes, never acoustic interfaces.

use super::*;
use crate::elastic::acoustic::viscoelastic::layers::{
    DrivenSupportedViscoelasticLayers, DrivenViscoelasticLayers, SupportedViscoelasticLayers,
    SupportedViscoelasticPower, ViscoelasticLayer, ViscoelasticLayerAbsorption,
};

#[derive(Clone, Copy)]
pub(super) enum NativeStack<'s, 'a> {
    TwoFluid(&'s ViscoelasticLayers<'a>),
    Supported(&'s SupportedViscoelasticLayers<'a>),
}
impl<'s, 'a> NativeStack<'s, 'a> {
    fn allows_empty(self) -> bool {
        matches!(self, Self::Supported(_))
    }
    fn layer_count(self) -> usize {
        match self {
            Self::TwoFluid(stack) => stack.layer_count(),
            Self::Supported(stack) => stack.layer_count(),
        }
    }
    fn layer(self, index: usize) -> Option<ViscoelasticLayer<'a>> {
        match self {
            Self::TwoFluid(stack) => stack.layer(index),
            Self::Supported(stack) => stack.layer(index),
        }
    }
    fn drive(
        self,
        input: [C64; 2],
        area_m2: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<NativeDrive<'s, 'a>, ViscoelasticSlabError> {
        match self {
            Self::TwoFluid(stack) => stack.drive(input, area_m2, cancelled).map(NativeDrive::TwoFluid),
            Self::Supported(stack) => stack
                .drive(input[0], area_m2, cancelled)
                .map(NativeDrive::Supported),
        }
    }
}

enum NativeDrive<'s, 'a> {
    TwoFluid(DrivenViscoelasticLayers<'s, 'a>),
    Supported(DrivenSupportedViscoelasticLayers<'s, 'a>),
}
impl NativeDrive<'_, '_> {
    fn layer_absorption(
        &self,
        index: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<ViscoelasticLayerAbsorption, ViscoelasticSlabError> {
        match self {
            Self::TwoFluid(drive) => drive.layer_absorption(index, cancelled),
            Self::Supported(drive) => drive.layer_absorption(index, cancelled),
        }
    }
    fn absorption_between(
        &self,
        index: usize,
        depths_m: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<ViscoelasticLayerAbsorption, ViscoelasticSlabError> {
        match self {
            Self::TwoFluid(drive) => drive.absorption_between(index, depths_m, cancelled),
            Self::Supported(drive) => drive.absorption_between(index, depths_m, cancelled),
        }
    }
    fn power(&self) -> NativePower {
        match self {
            Self::TwoFluid(drive) => NativePower::TwoFluid(drive.power()),
            Self::Supported(drive) => NativePower::Supported(drive.power()),
        }
    }
}

#[derive(Clone, Copy)]
pub(super) enum NativePower {
    TwoFluid(ViscoelasticLayersPower),
    Supported(SupportedViscoelasticPower),
}
impl NativePower {
    pub(super) fn components(self) -> HeatingPower {
        match self {
            Self::TwoFluid(power) => HeatingPower {
                incident: power.incident_power_w,
                outgoing: power.outgoing_power_w,
                heat: [power.bulk_absorbed_power_w, power.shear_absorbed_power_w],
                backing: 0.0,
            },
            Self::Supported(power) => HeatingPower {
                incident: [power.incident_power_w, 0.0],
                outgoing: [power.reflected_power_w, 0.0],
                heat: [power.bulk_absorbed_power_w, power.shear_absorbed_power_w],
                backing: power.backing_absorbed_power_w,
            },
        }
    }
}

pub(super) struct Prepared {
    pub(super) bins: Arc<[BoundBin]>,
    pub(super) power: NativePower,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    stack: NativeStack<'_, '_>,
    incident_peak_pressure_pa: [C64; 2],
    area_m2: f64,
    bins: &[AcousticThermalBin],
    thermal: &DepositionSimulation<'_>,
    budget: AcousticHeatingBudget,
    checkpoint: &mut impl FnMut() -> core::result::Result<(), DepositionError>,
) -> Result<Prepared> {
    checkpoint()?;
    if budget.maximum_bins > MAX_ACOUSTIC_THERMAL_BINS
        || budget.maximum_power_segments > MAX_ACOUSTIC_POWER_SEGMENTS
        || bins.len() > budget.maximum_bins
        || !nonnegative(budget.energy_absolute_tolerance_j)
    {
        return Err(AcousticHeatingError::Budget);
    }
    if (bins.is_empty() && !stack.allows_empty()) || !positive(area_m2) {
        return Err(AcousticHeatingError::Invalid(
            "acoustic heating needs positive area and thermal bins",
        ));
    }
    let mut ordered = bins.to_vec();
    ordered.sort_by(|a, b| {
        a.layer
            .cmp(&b.layer)
            .then(a.depths_m[0].total_cmp(&b.depths_m[0]))
    });
    let mut identities = BTreeSet::new();
    for bin in &ordered {
        checkpoint()?;
        let [lo, hi] = bin.temperature_band_k;
        if bin.layer >= stack.layer_count()
            || !positive(lo)
            || !positive(hi)
            || hi < lo
            || !identities.insert(bin.cell)
        {
            return Err(AcousticHeatingError::Invalid(
                "acoustic bin index, unique cell or temperature support is invalid",
            ));
        }
    }
    let mut interruption = None;
    let drive = stack.drive(incident_peak_pressure_pa, area_m2, &mut || {
        cancelled(checkpoint, &mut interruption)
    });
    let drive = acoustic_result(drive, interruption.take())?;
    let mut bound = Vec::with_capacity(ordered.len());
    let mut cursor = 0;
    for layer_index in 0..stack.layer_count() {
        checkpoint()?;
        let layer = stack.layer(layer_index).expect("original layer");
        let mut depth = 0.0;
        let mut heat = [Dd::ZERO; 2];
        let first = cursor;
        while cursor < ordered.len() && ordered[cursor].layer == layer_index {
            checkpoint()?;
            let bin = ordered[cursor];
            let [start, end] = bin.depths_m;
            if !nonnegative(start)
                || !positive(end)
                || start != depth
                || end <= start
                || end > layer.thickness_m
            {
                return Err(AcousticHeatingError::Invalid(
                    "acoustic depth bins must partition each layer without gaps or overlap",
                ));
            }
            let cell = thermal
                .cells()
                .get(&bin.cell)
                .ok_or(AcousticHeatingError::Invalid(
                    "acoustic thermal bin is inactive",
                ))?;
            let volume = product(area_m2, end - start)?;
            let mass = product(layer.source.density_kg_m3, volume)?;
            if !roundoff_equal(cell.geometry.volume_m3, volume)
                || !roundoff_equal(cell.mass_kg, mass)
            {
                return Err(AcousticHeatingError::Invalid(
                    "acoustic bin volume or mass disagrees with the original solid",
                ));
            }
            let local = drive.absorption_between(layer_index, [start, end], &mut || {
                cancelled(checkpoint, &mut interruption)
            });
            let local = acoustic_result(local, interruption.take())?;
            let q = [local.bulk_absorbed_power_w, local.shear_absorbed_power_w];
            let heat_w = checked((Dd::from_f64(q[0]) + Dd::from_f64(q[1])).to_f64())?;
            for channel in 0..2 {
                heat[channel] = heat[channel] + Dd::from_f64(q[channel]);
            }
            let binding = BoundBin {
                power: AcousticBinPower {
                    binding: bin,
                    source_identity: local.source_identity,
                    bulk_absorbed_power_w: q[0],
                    shear_absorbed_power_w: q[1],
                },
                cell: BoundCell::new(bin.cell, bin.temperature_band_k, thermal)?,
                heat_w,
            };
            bound.push(binding);
            depth = end;
            cursor += 1;
        }
        if depth != layer.thickness_m {
            return Err(AcousticHeatingError::Invalid(
                "every positive acoustic layer needs a complete thermal depth partition",
            ));
        }
        let original = drive.layer_absorption(layer_index, &mut || {
            cancelled(checkpoint, &mut interruption)
        });
        let original = acoustic_result(original, interruption.take())?;
        for (sum, total) in heat.into_iter().zip([
            original.bulk_absorbed_power_w,
            original.shear_absorbed_power_w,
        ]) {
            let discrepancy = (sum - Dd::from_f64(total)).to_f64().abs();
            let scale = sum.to_f64().max(total);
            let allowance = 4096.0 * f64::EPSILON * (cursor - first + 1) as f64 * scale;
            if !nonnegative(discrepancy) || discrepancy > allowance {
                return Err(AcousticHeatingError::Invalid(
                    "positive depth-bin heat does not close against the original layer",
                ));
            }
        }
    }
    if cursor != ordered.len() {
        return Err(AcousticHeatingError::Invalid("unbound acoustic depth bin"));
    }
    checkpoint()?;
    Ok(Prepared {
        bins: bound.into(),
        power: drive.power(),
    })
}
