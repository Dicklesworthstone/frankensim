//! Original material coefficients for four finite-resistance switching devices.
//!
//! Every device keeps its authored property keys and selected claims. Only the
//! receiving body's temperature and declared phase coordinate evolve here;
//! current, voltage, gate and frequency dependent laws need their own closure.

use super::*;

/// Explicit source roles for one channel and its antiparallel body diode.
/// All properties are static SI scalars at the actual thermal endpoint.
#[derive(Clone, Copy, Debug)]
pub struct DeviceConductionKeys<'a> {
    pub on_resistance: &'a PropertyKey,
    pub off_resistance: &'a PropertyKey,
    pub diode_forward_voltage: &'a PropertyKey,
    pub diode_differential_resistance: &'a PropertyKey,
    pub maximum_reverse_voltage: &'a PropertyKey,
    pub maximum_channel_current: &'a PropertyKey,
    pub maximum_diode_current: &'a PropertyKey,
}

impl<'a> DeviceConductionKeys<'a> {
    fn ordered(self) -> [&'a PropertyKey; 7] {
        [
            self.on_resistance,
            self.off_resistance,
            self.diode_forward_voltage,
            self.diode_differential_resistance,
            self.maximum_reverse_voltage,
            self.maximum_channel_current,
            self.maximum_diode_current,
        ]
    }
}

/// One finite device's sole material heat receiver and explicit source selection.
/// Several devices may share a body; each device's accepted heat is routed once.
#[derive(Clone, Debug)]
pub struct DeviceConductionBinding<'a> {
    pub cell: CellId,
    pub keys: DeviceConductionKeys<'a>,
    pub selection: MaterialPropertySelection,
}

/// Complete admitted device source and its original material/thermal binding.
#[derive(Clone, Debug)]
pub struct DeviceConductionMaterial {
    cell: CellId,
    keys: [PropertyKey; 7],
    source: ResolvedMaterialStatePoint,
    parameters: SwitchDevice,
}

/// Equal = the same material state (card, material, axes, query point, requirement, source
/// axes, selected claim, value bits; `same_material`) behind the same cell, keys and device
/// parameters. The receipt's selection-policy tag and hash are left out on purpose, because a
/// refresh at the same state re-resolves under `PinnedByProperty`: a source bound under
/// `SingleClaimOnly` and the same claim refreshed at the same state differ in exactly those two.
impl PartialEq for DeviceConductionMaterial {
    fn eq(&self, other: &Self) -> bool {
        self.cell == other.cell
            && self.keys == other.keys
            && self.parameters == other.parameters
            && same_material(&self.source, &other.source)
    }
}

impl DeviceConductionMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }

    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }

    pub const fn parameters(&self) -> SwitchDevice {
        self.parameters
    }
}

/// Resolve every device against a freshly admitted original material network.
/// The caller retains the same caloric/card bindings when preparing `sample`.
pub(super) fn resolve_devices(
    cells: &[MaterialCell<'_>],
    bindings: &[DeviceConductionBinding<'_>; 4],
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&[DeviceConductionMaterial; 4]>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<[DeviceConductionMaterial; 4], Error> {
    let dimensions = [
        Resistance::DIMS,
        Resistance::DIMS,
        Voltage::DIMS,
        Resistance::DIMS,
        Voltage::DIMS,
        Current::DIMS,
        Current::DIMS,
    ];
    let mut devices = Vec::with_capacity(bindings.len());
    for (index, binding) in bindings.iter().enumerate() {
        check(cancel)?;
        let body = cells
            .iter()
            .find(|body| body.cell == binding.cell)
            .ok_or_else(|| invalid("switching device needs its original material heat receiver"))?;
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or_else(|| invalid("switching device heat receiver is inactive"))?;
        if cell.material_id != body.material_id {
            return Err(invalid(
                "switching device receiver changed its original caloric material",
            ));
        }
        let phase = sample
            .cells
            .iter()
            .find(|entry| entry.cell == binding.cell)
            .ok_or_else(|| invalid("switching device receiver lacks its admitted caloric phase"))?
            .phase;
        let query = body.point.at(phase).map_err(point)?;
        let keys = binding.keys.ordered();
        for (key, dims) in keys.iter().zip(dimensions) {
            if key.quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid(
                    "device coefficients require their exact static resistance, voltage or current quantity",
                ));
            }
        }
        let material = if let Some(originals) = retained {
            let original = &originals[index];
            if original.cell != binding.cell
                || original
                    .keys
                    .iter()
                    .zip(&keys)
                    .any(|(old, key)| old != *key)
            {
                return Err(invalid(
                    "switching device cannot change its receiver or original source keys",
                ));
            }
            for (axis, value) in original.source.query_point() {
                if axis != body.point.temperature_axis
                    && Some(axis.as_str()) != body.point.liquid_fraction_axis
                    && query.axes().get(axis).map(|v| v.to_bits()) != Some(value.to_bits())
                {
                    return Err(invalid(
                        "switching device cannot change its held original source context",
                    ));
                }
            }
            original
                .source
                .at_query_point(body.point.material, &query)
                .map_err(source)?
        } else {
            let mut requirements = Vec::with_capacity(keys.len());
            for (role, key) in keys.iter().enumerate() {
                if keys[..role].contains(key) {
                    continue;
                }
                // Explicit reuse of a source key across roles is allowed. Any
                // strictly positive role retains the stronger common requirement.
                let positive = keys
                    .iter()
                    .enumerate()
                    .any(|(other, candidate)| other != 2 && candidate == key);
                let admissibility = if positive {
                    ScalarAdmissibility::StrictlyPositive
                } else {
                    ScalarAdmissibility::NonNegative
                };
                requirements.push(
                    ScalarPropertyRequirement::try_with_key(key, admissibility).map_err(source)?,
                );
            }
            resolve_material_state_point(
                body.point.material,
                &query,
                &requirements,
                binding.selection.clone(),
            )
            .map_err(source)?
        };
        let mut values = [0.0; 7];
        for (role, key) in keys.iter().enumerate() {
            check(cancel)?;
            let property = material
                .property_by_key(key)
                .ok_or_else(|| invalid("switching device lost an original coefficient property"))?;
            body.point.admit_static(property).map_err(point)?;
            let static_source = match body.point.liquid_fraction_axis {
                Some(axis) => property.has_only_owned_and_fixed_external_conditions(
                    axis,
                    QuantitySpec::dimensional(Dims::NONE),
                ),
                None => property.is_static_at_fixed_external_conditions(),
            };
            if !static_source {
                return Err(invalid(
                    "device coefficients cannot freeze unowned current, voltage, gate or frequency conditions",
                ));
            }
            let value = property.value_si();
            if !value.is_finite() || if role == 2 { value < 0.0 } else { value <= 0.0 } {
                return Err(invalid(
                    "device resistances and operating bounds must be positive; diode threshold must be nonnegative",
                ));
            }
            values[role] = value;
        }
        let parameters = SwitchDevice {
            on_resistance_ohm: values[0],
            off_resistance_ohm: values[1],
            body_diode: BodyDiode {
                forward_voltage_v: values[2],
                differential_resistance_ohm: values[3],
            },
            bounds: DeviceOperatingBounds {
                maximum_reverse_voltage_v: values[4],
                maximum_channel_current_a: values[5],
                maximum_diode_current_a: values[6],
            },
        };
        parameters.validate().map_err(model)?;
        devices.push(DeviceConductionMaterial {
            cell: binding.cell,
            keys: keys.map(|key| (*key).clone()),
            source: material,
            parameters,
        });
    }
    check(cancel)?;
    devices
        .try_into()
        .map_err(|_| invalid("four original switching device material bindings required"))
}
