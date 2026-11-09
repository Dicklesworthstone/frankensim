//! Source-qualified intrinsic gate storage, channel response and driver losses.
//!
//! The caller explicitly decomposes intrinsic drain-source Q(V), gate-source C
//! and gate-drain C. A gate-clamped Coss measurement cannot be reused as intrinsic
//! Qds and then receive another Miller capacitance. Storage laws stay fixed;
//! only dissipative G(Vgs,T) and driver resistances refresh with temperature.
//! Separate on/off driver resistors share one explicitly authored heat receiver.

use super::capacitors::{resolve_capacitances, resolve_constant_capacitances};
use super::charge::resolve_charge_curves;
use super::*;
use fs_matdb::{ClaimSelection, PropertyLineSlice};
use fs_phs::circuit::converter::commutation::gate::{GateConductanceCurve, GateDriver};
use fs_qty::Conductance;

/// Positive channel conductance in siemens along the owned gate-source voltage.
/// The original context must already declare this axis in volts. Thermal
/// coordinates may change the complete line; unowned drain-bias or frequency
/// response cannot be treated as a fixed gate-only law.
#[derive(Clone, Debug)]
pub struct GateConductanceBinding<'a> {
    pub cell: CellId,
    pub key: &'a PropertyKey,
    pub gate_voltage_axis: &'a str,
    pub gate_voltage_range_v: [f64; 2],
    pub selection: MaterialPropertySelection,
}

/// Complete original source, current thermal slice and native positive law.
#[derive(Clone, Debug, PartialEq)]
pub struct GateConductanceMaterial {
    cell: CellId,
    key: PropertyKey,
    gate_voltage_axis: String,
    gate_voltage_range_v: [f64; 2],
    source: ResolvedMaterialStatePoint,
    path: PropertyLineSlice,
    curve: GateConductanceCurve,
}
impl GateConductanceMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }
    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }
    pub fn gate_voltage_axis(&self) -> &str {
        &self.gate_voltage_axis
    }
    pub const fn gate_voltage_range_v(&self) -> [f64; 2] {
        self.gate_voltage_range_v
    }
    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }
    pub const fn source_path(&self) -> &PropertyLineSlice {
        &self.path
    }
    pub const fn curve(&self) -> &GateConductanceCurve {
        &self.curve
    }
}

/// Explicit intrinsic decomposition and separate finite gate drivers.
/// All four-device arrays use input high/low, then output high/low order.
/// Each driver resistor pair is [on, off]; both name the same authored receiver.
/// Driver command voltages and current limits are external controls, not a
/// claim that a finite gate supply reservoir has been modeled.
#[derive(Clone, Copy)]
pub struct GateMaterialBindings<'a> {
    pub intrinsic_drain_charge: &'a [DeviceChargeBinding<'a>; 4],
    pub gate_source_capacitance: &'a [DeviceCapacitanceBinding<'a>; 4],
    pub gate_drain_capacitance: &'a [DeviceCapacitanceBinding<'a>; 4],
    pub channel_conductance: &'a [GateConductanceBinding<'a>; 4],
    pub driver_resistors: &'a [ConverterThermalChannel<'a>; 8],
    pub drivers: [GateDriver; 4],
}

/// Accepted complete storage, conduction and driver sources.
/// Gate storage and the intrinsic charge laws are immutable along the history;
/// the retained G lines and driver resistances describe its actual temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct GateMaterial {
    intrinsic_drain_charge: [DeviceChargeMaterial; 4],
    gate_source_capacitance: [DeviceCapacitanceMaterial; 4],
    gate_drain_capacitance: [DeviceCapacitanceMaterial; 4],
    channel_conductance: [GateConductanceMaterial; 4],
    driver_resistors: [ConverterChannelMaterial; 8],
    authored_drivers: [GateDriver; 4],
    drivers: [GateDriver; 4],
}
impl GateMaterial {
    pub const fn intrinsic_drain_charge(&self) -> &[DeviceChargeMaterial; 4] {
        &self.intrinsic_drain_charge
    }
    pub const fn gate_source_capacitance(&self) -> &[DeviceCapacitanceMaterial; 4] {
        &self.gate_source_capacitance
    }
    pub const fn gate_drain_capacitance(&self) -> &[DeviceCapacitanceMaterial; 4] {
        &self.gate_drain_capacitance
    }
    pub const fn channel_conductance(&self) -> &[GateConductanceMaterial; 4] {
        &self.channel_conductance
    }
    pub const fn driver_resistors(&self) -> &[ConverterChannelMaterial; 8] {
        &self.driver_resistors
    }
    pub const fn drivers(&self) -> &[GateDriver; 4] {
        &self.drivers
    }
    pub fn driver_cells(&self) -> [CellId; 4] {
        core::array::from_fn(|i| self.driver_resistors[2 * i].cell())
    }
    pub(super) fn same_electrical_coefficients(&self, other: &Self) -> bool {
        self.drivers
            .iter()
            .zip(&other.drivers)
            .all(|(a, b)| driver_bits(*a) == driver_bits(*b))
            && self
                .channel_conductance
                .iter()
                .zip(&other.channel_conductance)
                .all(|(a, b)| a.curve == b.curve)
    }
}

fn driver_bits(driver: GateDriver) -> [u64; 5] {
    [
        driver.on_voltage_v,
        driver.off_voltage_v,
        driver.on_resistance_ohm,
        driver.off_resistance_ohm,
        driver.maximum_current_a,
    ]
    .map(f64::to_bits)
}

fn resolve_conductance(
    cells: &[MaterialCell<'_>],
    bindings: &[GateConductanceBinding<'_>; 4],
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&[GateConductanceMaterial; 4]>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<[GateConductanceMaterial; 4], Error> {
    let voltage = QuantitySpec::dimensional(Voltage::DIMS);
    let mut resolved = Vec::with_capacity(4);
    for (index, binding) in bindings.iter().enumerate() {
        check(cancel)?;
        let [lower_v, upper_v] = binding.gate_voltage_range_v;
        if binding.key.quantity() != QuantitySpec::dimensional(Conductance::DIMS)
            || !lower_v.is_finite()
            || !upper_v.is_finite()
            || lower_v >= upper_v
        {
            return Err(invalid(
                "gate conductance needs a siemens property and an ordered finite voltage domain",
            ));
        }
        let body = cells
            .iter()
            .find(|cell| cell.cell == binding.cell)
            .ok_or_else(|| invalid("gate conductance needs its original semiconductor package"))?;
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or_else(|| invalid("gate conductance package is inactive"))?;
        if cell.material_id != body.material_id {
            return Err(invalid(
                "gate conductance changed its original caloric material",
            ));
        }
        if binding.gate_voltage_axis == body.point.temperature_axis
            || Some(binding.gate_voltage_axis) == body.point.liquid_fraction_axis
            || body
                .point
                .context
                .axis_quantities()
                .get(binding.gate_voltage_axis)
                != Some(&voltage)
            || !body
                .point
                .context
                .axes()
                .contains_key(binding.gate_voltage_axis)
        {
            return Err(invalid(
                "gate conductance requires a distinct existing typed gate-source voltage coordinate",
            ));
        }
        let phase = sample
            .cells
            .iter()
            .find(|row| row.cell == binding.cell)
            .ok_or_else(|| invalid("gate conductance lacks its admitted caloric phase"))?
            .phase;
        let query = body
            .point
            .at(phase)
            .map_err(point)?
            .with(binding.gate_voltage_axis, lower_v)
            .map_err(|_| invalid("gate lower voltage is unrepresentable"))?;
        let material = if let Some(originals) = retained {
            let old = &originals[index];
            if old.cell != binding.cell
                || &old.key != binding.key
                || old.gate_voltage_axis != binding.gate_voltage_axis
                || old.gate_voltage_range_v.map(f64::to_bits)
                    != binding.gate_voltage_range_v.map(f64::to_bits)
            {
                return Err(invalid(
                    "gate conductance must retain its original package, key, axis and voltage support",
                ));
            }
            for (axis, value) in old.source.query_point() {
                if axis != body.point.temperature_axis
                    && Some(axis.as_str()) != body.point.liquid_fraction_axis
                    && axis != binding.gate_voltage_axis
                    && query.axes().get(axis).map(|v| v.to_bits()) != Some(value.to_bits())
                {
                    return Err(invalid(
                        "gate conductance changed its held original source context",
                    ));
                }
            }
            old.source
                .at_query_point(body.point.material, &query)
                .map_err(source)?
        } else {
            let requirement = ScalarPropertyRequirement::try_with_key(
                binding.key,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(source)?;
            resolve_material_state_point(
                body.point.material,
                &query,
                &[requirement],
                binding.selection.clone(),
            )
            .map_err(source)?
        };
        let property = material
            .property_by_key(binding.key)
            .ok_or_else(|| invalid("gate conductance lost its original source property"))?;
        body.point
            .admit_static_with_axes(property, &[(binding.gate_voltage_axis, voltage)])
            .map_err(point)?;
        let mut owned = vec![(binding.gate_voltage_axis, voltage)];
        if let Some(quantity) = query.axis_quantities().get(body.point.temperature_axis) {
            owned.push((body.point.temperature_axis, *quantity));
        }
        if let Some(axis) = body.point.liquid_fraction_axis {
            owned.push((axis, QuantitySpec::dimensional(Dims::NONE)));
        }
        if !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned) {
            return Err(invalid(
                "gate-only conductance cannot freeze an unowned drain-bias, frequency or other response coordinate",
            ));
        }
        let upper = query
            .clone()
            .with(binding.gate_voltage_axis, upper_v)
            .map_err(|_| invalid("gate upper voltage is unrepresentable"))?;
        let path = body
            .point
            .material
            .claims()
            .query_line_typed(
                binding.key,
                binding.gate_voltage_axis,
                &query,
                &upper,
                ClaimSelection::Pinned(property.answer().receipt.selected),
            )
            .map_err(|error| {
                source(fs_material::state_point::MaterialStatePointError::Support {
                    property: binding.key.name().to_owned(),
                    source: Box::new(error),
                })
            })?;
        let mut knots = Vec::with_capacity(path.knots.len());
        for (gate_v, answer) in &path.knots {
            check(cancel)?;
            if answer.receipt.selected != property.answer().receipt.selected {
                return Err(invalid(
                    "gate conductance line changed its original source claim",
                ));
            }
            knots.push((*gate_v, answer.evidence.value.value));
        }
        let curve = if let Some(old) = retained.map(|sources| &sources[index])
            && old.curve.knots() == knots.as_slice()
        {
            old.curve.clone()
        } else {
            GateConductanceCurve::try_new(knots).map_err(model)?
        };
        resolved.push(GateConductanceMaterial {
            cell: binding.cell,
            key: binding.key.clone(),
            gate_voltage_axis: binding.gate_voltage_axis.to_owned(),
            gate_voltage_range_v: binding.gate_voltage_range_v,
            source: material,
            path,
            curve,
        });
    }
    check(cancel)?;
    resolved
        .try_into()
        .map_err(|_| invalid("four complete gate conductance sources required"))
}

pub(super) fn resolve_gate_material(
    cells: &[MaterialCell<'_>],
    devices: &[DeviceConductionBinding<'_>; 4],
    bindings: &GateMaterialBindings<'_>,
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&GateMaterial>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<GateMaterial, Error> {
    for (i, device) in devices.iter().enumerate() {
        check(cancel)?;
        if bindings.intrinsic_drain_charge[i].cell != device.cell
            || bindings.gate_source_capacitance[i].cell != device.cell
            || bindings.gate_drain_capacitance[i].cell != device.cell
            || bindings.channel_conductance[i].cell != device.cell
        {
            return Err(invalid(
                "intrinsic drain charge, gate storage, channel law and diode must share their original semiconductor package",
            ));
        }
        if bindings.driver_resistors[2 * i].cell != bindings.driver_resistors[2 * i + 1].cell {
            return Err(invalid(
                "one gate driver's on/off resistors need the same explicitly authored heat receiver",
            ));
        }
        bindings.drivers[i].validate().map_err(model)?;
        if retained.is_some_and(|old| {
            driver_bits(old.authored_drivers[i]) != driver_bits(bindings.drivers[i])
        }) {
            return Err(invalid(
                "gate driver history must retain its authored command voltages, resistance references and current limit",
            ));
        }
    }
    let intrinsic_drain_charge = resolve_charge_curves(
        cells,
        bindings.intrinsic_drain_charge,
        thermal,
        sample,
        retained.map(|old| &old.intrinsic_drain_charge),
        cancel,
    )?;
    let gate_source_capacitance = resolve_capacitances(
        cells,
        bindings.gate_source_capacitance,
        thermal,
        sample,
        retained.map(|old| &old.gate_source_capacitance),
        cancel,
    )?;
    let gate_drain_capacitance = resolve_constant_capacitances(
        cells,
        bindings.gate_drain_capacitance,
        thermal,
        sample,
        retained.map(|old| &old.gate_drain_capacitance),
        true,
        cancel,
    )?;
    let channel_conductance = resolve_conductance(
        cells,
        bindings.channel_conductance,
        thermal,
        sample,
        retained.map(|old| &old.channel_conductance),
        cancel,
    )?;
    let fixed = core::array::from_fn(|i| {
        if i % 2 == 0 {
            bindings.drivers[i / 2].on_resistance_ohm
        } else {
            bindings.drivers[i / 2].off_resistance_ohm
        }
    });
    let driver_resistors = resolve_channels(
        cells,
        bindings.driver_resistors,
        thermal,
        sample,
        fixed,
        retained.map(|old| &old.driver_resistors),
        false,
        cancel,
    )?;
    let mut drivers = bindings.drivers;
    for i in 0..4 {
        check(cancel)?;
        drivers[i].on_resistance_ohm = driver_resistors[2 * i].resistance_ohm();
        drivers[i].off_resistance_ohm = driver_resistors[2 * i + 1].resistance_ohm();
        drivers[i].validate().map_err(model)?;
    }
    Ok(GateMaterial {
        intrinsic_drain_charge,
        gate_source_capacitance,
        gate_drain_capacitance,
        channel_conductance,
        driver_resistors,
        authored_drivers: bindings.drivers,
        drivers,
    })
}
