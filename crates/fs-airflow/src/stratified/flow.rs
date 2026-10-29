use super::{
    Acceleration, Density, Layer, LayerMassRates, Length, MAX_VENT_SLABS, MassFlowRate, Pressure,
    SlabFlow, Sum, TwoLayerColumn, VentError, VentFlow, WallOpening, det, finite, poll, positive,
    validate_gravity,
};

/// Integrate a vertical opening between two source-declared hydrostatic columns.
///
/// Positive slab rates travel left to right. Each transfer retains both layer
/// labels and uses its actual donor's density. Work is bounded by six slabs;
/// successful evaluation performs no heap allocation. Cancellation is checked
/// before evaluation, at each interval/slab, and immediately before publication.
pub fn wall_vent_flow(
    opening: &WallOpening,
    left: TwoLayerColumn,
    right: TwoLayerColumn,
    gravity: Acceleration,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VentFlow, VentError> {
    poll(&mut cancelled)?;
    validate_gravity(gravity)?;
    let mut result = VentFlow {
        slabs: [SlabFlow::EMPTY; MAX_VENT_SLABS],
        count: 0,
        rates: LayerMassRates {
            left: [MassFlowRate::new(0.0); 2],
            right: [MassFlowRate::new(0.0); 2],
        },
    };
    if opening.width.value() == 0.0 || opening.bottom == opening.top {
        poll(&mut cancelled)?;
        return Ok(result);
    }
    // Constant positive layer densities make pressure monotone between these
    // endpoints, so endpoint admission covers the complete opening.
    for column in [left, right] {
        column.pressure_at(opening.bottom, gravity)?;
        column.pressure_at(opening.top, gravity)?;
    }
    let bottom = opening.bottom.value();
    let top = opening.top.value();
    let mut cuts = [bottom, top, top, top];
    let mut count = 2;
    for column in [left, right] {
        let interface = column.interface_height.value();
        if interface > bottom && interface < top {
            cuts[count] = interface;
            count += 1;
        }
    }
    cuts[..count].sort_unstable_by(f64::total_cmp);
    let mut unique = 1;
    for i in 1..count {
        if cuts[i] != cuts[unique - 1] {
            cuts[unique] = cuts[i];
            unique += 1;
        }
    }
    // Subtract reference pressures before adding hydrostatic heads. Computing
    // two large absolute pressures first would discard small buoyant drives.
    let mut difference = left.reference_pressure.value() - right.reference_pressure.value();
    if gravity.value() != 0.0 {
        difference = finite(
            difference
                - gravity.value()
                    * (left.density_height(opening.bottom)?
                        - right.density_height(opening.bottom)?),
            "opening pressure difference",
        )?;
    }
    let inverse_sqrt_resistance =
        1.0 / det::sqrt(opening.discharge.reference_loss.resistance.value());
    let mut rates = [Sum::default(); 4];
    for i in 0..unique - 1 {
        poll(&mut cancelled)?;
        let b = cuts[i];
        let t = cuts[i + 1];
        let height = t - b;
        let l = layer_below_top(left, t);
        let r = layer_below_top(right, t);
        let ld = left.density(l);
        let rd = right.density(r);
        let slope = finite(
            gravity.value() * (rd.value() - ld.value()),
            "pressure slope",
        )?;
        let end_difference = finite(difference + slope * height, "slab pressure difference")?;
        let crossing = (difference < 0.0 && end_difference > 0.0)
            || (difference > 0.0 && end_difference < 0.0);
        if crossing {
            // Ratios avoid overflow in |dp_bottom|+|dp_top|. Keep local extents
            // independently of rounded global neutral-plane coordinates, so a
            // very small counterflow is never snapped to zero at an interface.
            let scale = difference.abs().max(end_difference.abs());
            let a = difference.abs() / scale;
            let z = end_difference.abs() / scale;
            let lower_height = height * (a / (a + z));
            let upper_height = height * (z / (a + z));
            if !positive(lower_height) || !positive(upper_height) {
                return Err(VentError::NonFinite(
                    "nonzero neutral-plane slab extent underflowed",
                ));
            }
            let neutral = if lower_height <= upper_height {
                b + lower_height
            } else {
                t - upper_height
            };
            append(
                &mut result,
                &mut rates,
                opening.width.value(),
                inverse_sqrt_resistance,
                b,
                neutral,
                lower_height,
                difference,
                0.0,
                l,
                r,
                ld,
                rd,
                &mut cancelled,
            )?;
            append(
                &mut result,
                &mut rates,
                opening.width.value(),
                inverse_sqrt_resistance,
                neutral,
                t,
                upper_height,
                0.0,
                end_difference,
                l,
                r,
                ld,
                rd,
                &mut cancelled,
            )?;
        } else {
            append(
                &mut result,
                &mut rates,
                opening.width.value(),
                inverse_sqrt_resistance,
                b,
                t,
                height,
                difference,
                end_difference,
                l,
                r,
                ld,
                rd,
                &mut cancelled,
            )?;
        }
        difference = end_difference;
    }
    result.rates = LayerMassRates {
        left: [
            MassFlowRate::new(rates[0].value),
            MassFlowRate::new(rates[1].value),
        ],
        right: [
            MassFlowRate::new(rates[2].value),
            MassFlowRate::new(rates[3].value),
        ],
    };
    poll(&mut cancelled)?;
    Ok(result)
}

fn layer_below_top(column: TwoLayerColumn, top: f64) -> Layer {
    // Each interval was already split at this interface. This comparison also
    // works when the midpoint would round onto a boundary at large elevations.
    if top <= column.interface_height.value() {
        Layer::Lower
    } else {
        Layer::Upper
    }
}

#[allow(clippy::too_many_arguments)]
fn append(
    result: &mut VentFlow,
    rates: &mut [Sum; 4],
    width: f64,
    inverse_sqrt_resistance: f64,
    bottom: f64,
    top: f64,
    height: f64,
    pressure_bottom: f64,
    pressure_top: f64,
    left_layer: Layer,
    right_layer: Layer,
    left_density: Density,
    right_density: Density,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), VentError> {
    poll(cancelled)?;
    if pressure_bottom == 0.0 && pressure_top == 0.0 {
        return Ok(());
    }
    let sign = if pressure_bottom == 0.0 {
        pressure_top
    } else {
        pressure_bottom
    };
    let upstream_density = if sign > 0.0 {
        left_density
    } else {
        right_density
    };
    let large = det::sqrt(pressure_bottom.abs().max(pressure_top.abs()));
    let small = det::sqrt(pressure_bottom.abs().min(pressure_top.abs()));
    let ratio = small / large;
    // Rationalized exact integral, avoiding cancellation as both endpoint
    // pressures approach equality. It also supplies the constant-dp limit:
    // mean(sqrt(|dp|)) = 2/3 * (a²+ab+b²)/(a+b), a,b=sqrt(|dp|).
    let mean_root = large * ((2.0 / 3.0) * (1.0 + ratio + ratio * ratio) / (1.0 + ratio));
    let magnitude = finite(
        (width * height)
            * (inverse_sqrt_resistance * mean_root)
            * det::sqrt(upstream_density.value()),
        "integrated quadratic mass flow",
    )?;
    if magnitude == 0.0 {
        return Err(VentError::NonFinite("nonzero slab mass flow underflowed"));
    }
    let signed = magnitude.copysign(sign);
    rates[left_layer.index()].add(-signed)?;
    rates[2 + right_layer.index()].add(signed)?;
    if result.count == MAX_VENT_SLABS {
        return Err(VentError::Invalid("two-layer slab bound exceeded"));
    }
    result.slabs[result.count] = SlabFlow {
        bottom: Length::new(bottom),
        top: Length::new(top),
        height: Length::new(height),
        left_layer,
        right_layer,
        mass_rate: MassFlowRate::new(signed),
        upstream_density,
        pressure_difference_bottom: Pressure::new(pressure_bottom),
        pressure_difference_top: Pressure::new(pressure_top),
    };
    result.count += 1;
    Ok(())
}
