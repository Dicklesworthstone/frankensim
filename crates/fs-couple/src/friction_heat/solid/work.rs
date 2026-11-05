//! Explicit plane-strain work conversion and declared heat generation.

pub(crate) fn partition(
    work_j_per_m: f64,
    depth_m: f64,
    heat_fraction: f64,
) -> Result<(f64, f64), &'static str> {
    if !work_j_per_m.is_finite()
        || work_j_per_m < 0.0
        || !depth_m.is_finite()
        || depth_m <= 0.0
        || !heat_fraction.is_finite()
        || !(0.0..=1.0).contains(&heat_fraction)
    {
        return Err("finite nonnegative native work, positive depth and a heat fraction required");
    }
    let work = work_j_per_m * depth_m;
    let heat = work * heat_fraction;
    let other = work - heat;
    if !work.is_finite()
        || !heat.is_finite()
        || !other.is_finite()
        || (work_j_per_m > 0.0 && work == 0.0)
        || (work > 0.0 && heat_fraction > 0.0 && heat == 0.0)
        || (work > 0.0 && heat_fraction < 1.0 && other == 0.0)
    {
        return Err("plane-strain friction work or its declared partition is unrepresentable");
    }
    Ok((heat, other))
}

#[cfg(test)]
mod tests {
    use super::partition;

    #[test]
    fn explicit_depth_and_heat_conversion_close_the_original_work() {
        let (heat, other) = partition(8.0, 0.125, 0.75).unwrap();
        assert_eq!((heat, other), (0.75, 0.25));
        assert_eq!(heat + other, 8.0 * 0.125);
        assert_eq!(
            partition(8.0, 0.25, 0.75).unwrap(),
            (2.0 * heat, 2.0 * other)
        );
        assert_eq!(partition(8.0, 0.125, 0.0).unwrap(), (0.0, 1.0));
        assert_eq!(partition(8.0, 0.125, 1.0).unwrap(), (1.0, 0.0));
        assert_eq!(partition(0.0, 0.125, 0.75).unwrap(), (0.0, 0.0));
    }

    #[test]
    fn mechanical_energy_creation_and_unrepresentable_conversion_refuse() {
        for input in [
            (-1.0, 1.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.1),
            (f64::MAX, 2.0, 0.5),
            (f64::MIN_POSITIVE, f64::MIN_POSITIVE, 1.0),
            (f64::NAN, 1.0, 1.0),
        ] {
            assert!(partition(input.0, input.1, input.2).is_err());
        }
    }
}
