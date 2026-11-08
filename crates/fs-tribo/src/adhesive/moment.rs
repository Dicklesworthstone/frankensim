//! Private exact increment of an elastic/perfectly plastic scalar moment.
//!
//! Rolling and torsion supply their own physical stiffness, calibration and
//! history identity. This arithmetic owner only integrates the common bounded
//! spring/slider path; it cannot select either physical model.

use fs_math::dd::Dd;

pub(crate) struct MomentIncrement {
    pub fraction: Dd,
    pub angle: Dd,
    pub endpoint_moment: f64,
    pub mean_moment: f64,
    pub storage_before: f64,
    pub storage_after: f64,
    pub storage_change: f64,
    pub loss: f64,
    pub work: f64,
    pub plastic_fraction: Dd,
}

pub(crate) fn advance(
    old_fraction: Dd,
    old_angle: Dd,
    angle_increment: f64,
    critical_angle: f64,
    critical_moment: f64,
    yield_energy: f64,
) -> Result<MomentIncrement, &'static str> {
    let angle = old_angle + Dd::from_f64(angle_increment);
    finite(angle.hi, "accumulated angle")?;
    finite(angle.lo, "accumulated angle residual")?;
    if angle_increment != 0.0 && angle == old_angle {
        return Err("accumulated angle increment");
    }
    let increment = Dd::from_f64(angle_increment) / Dd::from_f64(critical_angle);
    finite(increment.hi, "normalized angular increment")?;
    finite(increment.lo, "normalized angular increment residual")?;
    if angle_increment != 0.0 && increment == Dd::ZERO {
        return Err("angular increment");
    }
    let trial = old_fraction + increment;
    if increment != Dd::ZERO && trial == old_fraction {
        return Err("elastic history increment");
    }
    let mut fraction = trial;
    if Dd::ONE.lt(trial) {
        fraction = Dd::ONE;
    }
    if trial.lt(-Dd::ONE) {
        fraction = -Dd::ONE;
    }
    let plastic = trial - fraction;
    let change = fraction - old_fraction;
    let energy = Dd::from_f64(yield_energy);
    let storage_change = energy * change * (fraction + old_fraction);
    let loss = energy * Dd::from_f64(2.0) * plastic.abs();
    let work = storage_change + loss;
    // Use the analytic midpoint on the elastic branch, including a zero step.
    // The yielding numerator retains small plastic work after elastic reversal.
    let mean_fraction = if plastic == Dd::ZERO {
        (old_fraction + fraction) * Dd::from_f64(0.5)
    } else {
        (change * (fraction + old_fraction) * Dd::from_f64(0.5) + plastic.abs()) / increment
    };
    let step = MomentIncrement {
        fraction,
        angle,
        endpoint_moment: -critical_moment * fraction.to_f64(),
        mean_moment: -critical_moment * mean_fraction.to_f64(),
        storage_before: (energy * old_fraction * old_fraction).to_f64(),
        storage_after: (energy * fraction * fraction).to_f64(),
        storage_change: storage_change.to_f64(),
        loss: loss.to_f64(),
        work: work.to_f64(),
        plastic_fraction: plastic,
    };
    for (field, value) in [
        ("endpoint moment", step.endpoint_moment),
        ("mean moment", step.mean_moment),
        ("storage before", step.storage_before),
        ("storage after", step.storage_after),
        ("storage change", step.storage_change),
        ("irreversible loss", step.loss),
        ("applied work", step.work),
    ] {
        finite(value, field)?;
    }
    Ok(step)
}

fn finite(value: f64, field: &'static str) -> Result<(), &'static str> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(field)
    }
}
