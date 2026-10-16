//! Physical input isolation and joint-corner accounting for native budgets.
//! A declared radiative reservoir is not an alias of the fluid ambient.
use super::{ProjectSpec, PropagatedTerm, SolveRefusal, conduction_error};

mod joint;
pub(super) use joint::{MAX_CORNERS, account_joint_design};

/// Replace only the fluid temperatures covered by the ambient envelope.
/// Reservoirs, prescribed solid temperatures, material query temperatures,
/// cards, geometry and all other inputs retain their original declarations.
/// The caller has already admitted the finite, positive envelope endpoint.
pub(super) fn with_fluid_temperature(spec: &ProjectSpec, temperature_k: f64) -> ProjectSpec {
    let mut perturbed = spec.clone();
    if let Some(setup) = perturbed.cooling.as_mut().and_then(|c| c.conduction.as_mut()) {
        for boundary in &mut setup.boundaries {
            use fs_project::ThermalBoundaryCondition as B;
            match &mut boundary.condition {
                B::Convection { reference_temperature, .. } => reference_temperature.value = temperature_k,
                B::AirflowConvection { inlet_temperature, .. } => inlet_temperature.value = temperature_k,
                B::NaturalConvection { ambient_temperature, .. } => ambient_temperature.value = temperature_k,
                _ => {}
            }
        }
    }
    perturbed
}

/// Charge a successfully evaluated joint corner, or retain an explicit gap.
/// `joint` contains only the recoverable physical refusal: cancellation and
/// work-envelope errors must propagate from the caller BEFORE this function.
/// Unknown interaction is returned as NaN for the existing JSON-null encoding,
/// never as a measured zero. The individually measured model term is unchanged.
pub(super) fn account_joint_corner(
    boundary: PropagatedTerm,
    model_half_width_k: f64,
    nominal_k: f64,
    joint: Result<(String, f64), String>,
) -> Result<(PropagatedTerm, f64), SolveRefusal> {
    let bad = |what| conduction_error("cli-solve-conduction-propagation", what,
        "retain finite nonnegative input widths and the actual joint-corner result");
    let PropagatedTerm::Measured { half_width_k, method, detail, vertices } = boundary else {
        return Err(bad("a joint correction requires measured individual boundary vertices"));
    };
    if ![half_width_k, model_half_width_k].iter().all(|v| v.is_finite() && *v >= 0.0)
        || !nominal_k.is_finite()
    {
        return Err(bad("invalid individual widths or nominal value in joint-corner accounting"));
    }
    let (label, value) = match joint {
        Ok(row) => row,
        Err(reason) => {
            let isolated = vertices.iter().map(|(label, value)| format!("{label} = {value} K"))
                .collect::<Vec<_>>().join("; ");
            return Ok((PropagatedTerm::Unmeasured { reason: format!(
                "joint boundary/model-form corner could not be evaluated: {reason}; individual boundary half-width {half_width_k} K and model-form half-width {model_half_width_k} K do not establish the combined response; isolated boundary vertices: {isolated}"
            ) }, f64::NAN));
        }
    };
    let deviation = (value - nominal_k).abs();
    let separate = half_width_k + model_half_width_k;
    if !value.is_finite() || !deviation.is_finite() || !separate.is_finite() {
        return Err(bad("nonfinite joint-corner value or width arithmetic"));
    }
    let excess = (deviation - separate).max(0.0);
    let total = half_width_k + excess;
    if !total.is_finite() { return Err(bad("joint-corrected boundary width overflow")); }
    Ok((PropagatedTerm::Measured {
        half_width_k: total, method,
        detail: format!("{detail}; {label} = {value} K; includes {excess} K joint-corner interaction excess"),
        vertices,
    }, excess))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_project::ThermalBoundaryCondition as B;

    fn project() -> ProjectSpec {
        fs_project::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-radiation.fsim"))).unwrap().decoded.spec
    }
    fn measured() -> PropagatedTerm {
        PropagatedTerm::Measured { half_width_k: 2.0, method: "interval-vertex-resolve",
            detail: "individual boundary solves".into(),
            vertices: vec![("cold".into(), 298.0), ("hot".into(), 302.0)] }
    }

    #[test]
    fn ambient_vertices_preserve_independent_reservoirs_and_other_physical_inputs() {
        for natural in [false, true] {
            let mut base = project();
            let setup = base.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            let radiation = setup.radiation.as_mut().unwrap();
            radiation.surfaces[0].reservoir_temperature.value = 350.0;
            let mut second = radiation.surfaces[0].clone();
            second.name = "independent-cold-reservoir".into();
            second.reservoir_temperature.value = 270.0;
            radiation.surfaces.push(second);
            if natural {
                setup.boundaries[0].condition = B::NaturalConvection {
                    characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
                    ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
                    correlation: "convection.churchill-chu-vertical-plate".into(),
                };
            }
            let mut fixed = setup.boundaries[0].clone();
            fixed.target = "fixed-solid-support".into();
            fixed.condition = B::FixedTemperature {
                temperature: fs_qty::QtyAny::new(280.0, fs_project::spec::dims::TEMPERATURE),
            };
            setup.boundaries.push(fixed);
            for value in [290.0, 300.0, 310.0] {
                let mut expected = base.clone();
                let condition = &mut expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition;
                match condition {
                    B::Convection { reference_temperature, .. } => reference_temperature.value = value,
                    B::NaturalConvection { ambient_temperature, .. } => ambient_temperature.value = value,
                    _ => unreachable!(),
                }
                // Exercise the actual driver entry, not only an unconnected helper.
                assert_eq!(super::super::with_inlet_temperature(&base, value), expected);
            }
        }
    }

    #[test]
    fn joint_refusal_is_unknown_not_a_zero_interaction_or_measured_boundary() {
        let (term, excess) = account_joint_corner(measured(), 3.0, 300.0,
            Err("weakened cooling leaves the material temperature span".into())).unwrap();
        assert!(excess.is_nan());
        assert!(term.half_width().is_none());
        let PropagatedTerm::Unmeasured { reason } = term else { panic!("joint refusal must remain NO-DATA") };
        assert!(reason.contains("material temperature span"));
        assert!(reason.contains("hot = 302 K"));
    }

    #[test]
    fn successful_joint_corner_charges_only_the_excess_and_preserves_isolated_vertices() {
        for (value, expected_excess) in [(304.0, 0.0), (305.0, 0.0), (309.0, 4.0), (291.0, 4.0)] {
            let (term, excess) = account_joint_corner(measured(), 3.0, 300.0,
                Ok(("joint worst corner".into(), value))).unwrap();
            assert_eq!(excess, expected_excess);
            let PropagatedTerm::Measured { half_width_k, vertices, detail, .. } = term else { panic!("measured") };
            assert_eq!(half_width_k, 2.0 + expected_excess);
            assert_eq!(vertices, vec![("cold".into(), 298.0), ("hot".into(), 302.0)]);
            assert!(detail.contains(&format!("joint worst corner = {value} K")));
        }
    }

    #[test]
    fn joint_arithmetic_never_turns_nonfinite_evidence_into_zero() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(account_joint_corner(measured(), 3.0, 300.0,
                Ok(("joint".into(), value))).is_err());
            assert!(account_joint_corner(measured(), value, 300.0,
                Ok(("joint".into(), 300.0))).is_err());
        }
        assert!(account_joint_corner(measured(), -1.0, 300.0, Ok(("joint".into(), 300.0))).is_err());
        assert!(account_joint_corner(measured(), 3.0, -f64::MAX, Ok(("joint".into(), f64::MAX))).is_err());
    }
}
