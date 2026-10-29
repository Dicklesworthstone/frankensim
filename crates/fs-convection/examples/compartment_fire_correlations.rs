//! Compare actual MQH predictions with the six retained Steckler layer means.
//!
//! From the repository root, explicitly declare an effective wall coefficient:
//! `cargo run -p fs-convection --example compartment_fire_correlations -- \
//! data/reference/compartment_fire/steckler_1982_cases.csv 10`
//! Here 10 W/(m² K) is an engineering illustration, not a measured or fitted h.
//! The compact CSV has a fixed unquoted schema; provenance is retained beside it.

use fs_convection::HeatTransferCoefficient;
use fs_convection::fire::validation::{FirePredictionPairV1, compare_fire_predictions_v1};
use fs_convection::fire::{FireCorrelationIdV1, FireOpeningV1, mqh_hot_layer_temperature};
use fs_qty::{Area, Length, Power, Temperature};
use std::error::Error;

const HEADER: &str = "case_id,room_width_m,room_depth_m,room_height_m,door_width_m,door_height_m,prescribed_heat_release_rate_w,cfast_ambient_temperature_k,reduced_upper_temperature_k";
const USAGE: &str = "usage: compartment_fire_correlations CASES.csv EFFECTIVE_WALL_H_W_M2_K\n\
    Both arguments are required; h must be a finite positive value.\n\
    Example: compartment_fire_correlations data/reference/compartment_fire/steckler_1982_cases.csv 10\n\
    That h=10 is an explicit engineering illustration, not measured or calibrated.";

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    if arguments.len() != 2 {
        return Err(USAGE.into());
    }
    let h = positive(&arguments[1], "effective wall h [W/(m2 K)]")?;
    let input = std::fs::read_to_string(&arguments[0])?;
    let mut rows = input.lines().enumerate().filter(|(_, line)| {
        let line = line.trim();
        !line.is_empty() && !line.starts_with('#')
    });
    let (_, header) = rows.next().ok_or("CSV contains no header")?;
    if header.trim() != HEADER {
        return Err(format!("expected fixed CSV header: {HEADER}").into());
    }
    println!(
        "# MQH nominal hot-layer correlation; {}",
        FireCorrelationIdV1::MqhHotLayer.card().source_url
    );
    println!("# Cases: {}", arguments[0]);
    println!("# Explicit caller-selected effective wall h = {h:.9} W/(m2 K)");
    println!(
        "# This run makes no claim that h is measured, calibrated, or a validated wall model."
    );
    println!(
        "# CFAST wall properties are an authored approximation; they do not determine h here."
    );
    println!(
        "# Observations are reduced upper-layer means; the reference baseline is each CFAST input ambient."
    );
    println!(
        "# These are steady-state comparisons, not transient measurements or raw profile maxima."
    );
    println!(
        "case_id,ambient_k,door_width_m,predicted_upper_rise_k,observed_upper_rise_k,model_over_observed"
    );
    let mut pairs = Vec::new();
    for (index, row) in rows {
        let line = index + 1;
        let columns: Vec<_> = row.split(',').map(str::trim).collect();
        if columns.len() != 9 || columns[0].is_empty() || row.contains('"') {
            return Err(format!("line {line}: expected nine unquoted CSV fields").into());
        }
        let mut values = [0.0; 8];
        for (index, value) in values.iter_mut().enumerate() {
            *value = positive(
                columns[index + 1],
                &format!("line {line}, column {}", index + 2),
            )?;
        }
        let [
            width,
            depth,
            height,
            door_width,
            door_height,
            q,
            ambient,
            observed_upper,
        ] = values;
        if door_width > width || door_height > height {
            return Err(format!("line {line}: door dimensions exceed the front wall").into());
        }
        let opening_area = door_width * door_height;
        let enclosing_area = 2.0 * (width * depth + width * height + depth * height) - opening_area;
        let opening = FireOpeningV1::new(Area::new(opening_area), Length::new(door_height))?;
        let prediction = mqh_hot_layer_temperature(
            Power::new(q),
            Temperature::new(ambient),
            opening,
            Area::new(enclosing_area),
            HeatTransferCoefficient::new(h),
        )?
        .value();
        let observed_rise = observed_upper - ambient;
        if !observed_rise.is_finite() || observed_rise <= 0.0 {
            return Err(format!(
                "line {line}: observed upper-layer rise must be finite and positive"
            )
            .into());
        }
        let predicted_rise = prediction.temperature_rise.value();
        println!(
            "{},{ambient:.2},{door_width:.2},{predicted_rise:.9},{observed_rise:.9},{:.9}",
            columns[0],
            predicted_rise / observed_rise
        );
        pairs.push(FirePredictionPairV1 {
            model: predicted_rise,
            observed: observed_rise,
        });
    }
    // The retained reduced-QoI uncertainty is unknown. Instrument error limits
    // are not substituted for a combined experimental standard uncertainty.
    let comparison = compare_fire_predictions_v1(&pairs, None)?;
    println!("# pair_count={}", comparison.pair_count);
    println!("# mean_log_ratio={:.9}", comparison.mean_log_ratio);
    println!(
        "# geometric_model_observation_ratio={:.9}",
        comparison.geometric_ratio
    );
    println!(
        "# combined_log_scatter={:.9}",
        comparison.combined_log_scatter
    );
    println!(
        "# separation={:?}; model_relative_sigma={:?}; fds_bias_factor={:?}",
        comparison.separation, comparison.model_relative_sigma, comparison.fds_bias_factor
    );
    println!(
        "# Reduced-QoI experimental uncertainty is unknown; no separated model uncertainty or adjusted FDS bias is claimed."
    );
    println!(
        "# These are descriptive comparisons, without a fitted h, confidence interval, or validation pass."
    );
    Ok(())
}

fn positive(text: &str, field: &str) -> Result<f64, Box<dyn Error>> {
    let value: f64 = text.parse().map_err(|error| format!("{field}: {error}"))?;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(format!("{field} must be finite and positive; found {value}").into())
    }
}
