//! Frozen-source scalar inversion of original Reynolds-dependent resistances.

use super::*;

pub(super) fn invert<F>(
    sections: &[PressureSection<F>],
    snapshots: &[FluidTransportState],
    boussinesq: Option<BoussinesqDensity>,
    pressure: f64,
    budget: SeriesPressureBudget,
    gate: &CancelGate,
) -> Result<f64> {
    // Re = |m| D/(mu A), including a declared Boussinesq reference density.
    // Intersect actual model support before evaluating any pressure candidate.
    let (mut lower, mut upper) = (0.0_f64, f64::MAX);
    for (section, snapshot) in sections.iter().zip(snapshots) {
        poll(gate)?;
        let law = section.resistance;
        let mass_per_re = finite(snapshot.dynamic_viscosity_pa_s()
            * law.flow_area.value() / law.reynolds_length.value())?;
        if mass_per_re <= 0.0 {
            return Err(invalid("unrepresentable Reynolds-to-mass conversion"));
        }
        lower = lower.max(finite(mass_per_re * law.reynolds_range.0)?);
        upper = upper.min(finite(mass_per_re * law.reynolds_range.1)?);
    }
    if lower >= upper {
        return Err(invalid("series resistance Reynolds intervals have no common flow interval"));
    }
    let evaluate = |mass| -> Result<(f64, f64)> {
        let (mut drop, mut tangent) = (Sum::default(), Sum::default());
        for (index, (section, snapshot)) in sections.iter().zip(snapshots).enumerate() {
            poll(gate)?;
            let rho = boussinesq.map_or(snapshot.density_kg_m3(), BoussinesqDensity::reference_kg_m3);
            let value = section.resistance.evaluate_snapshot(*snapshot, MassFlowRate::new(mass), rho)
                .map_err(|source| SeriesPressureError::Section { index, source })?;
            drop.add(value.pressure_drop_pa)?;
            tangent.add(value.dpressure_dmass_flow)?;
        }
        Ok((drop.total()?, tangent.total()?))
    };
    // A rounded Re -> mass -> speed -> Re conversion may land one ULP outside
    // a closed bound. Move only a refused endpoint inward, with a fixed cap;
    // never widen model support or modify a successfully admitted endpoint.
    let endpoint = |mut mass: f64, inward: bool| -> Result<(f64, f64)> {
        for _ in 0..8 {
            match evaluate(mass) {
                Ok((drop, _)) => return Ok((mass, drop)),
                Err(SeriesPressureError::Section { source: ResistanceError::Reynolds { value, admitted }, .. })
                    if (inward && value < admitted.0) || (!inward && value > admitted.1) => {
                        mass = if inward { mass.next_up() } else { mass.next_down() };
                    }
                Err(error) => return Err(error),
            }
        }
        Err(invalid("Reynolds support endpoint is not representable"))
    };
    let (admitted_lower, low_drop) = endpoint(lower, true)?;
    let (admitted_upper, high_drop) = endpoint(upper, false)?;
    lower = admitted_lower;
    upper = admitted_upper;
    if lower >= upper || low_drop >= high_drop {
        return Err(invalid("series resistance has no representable monotone flow interval"));
    }
    let target = pressure.abs();
    if target < low_drop || target > high_drop {
        return Err(invalid("driving pressure has no root within the original Reynolds support"));
    }
    let tolerance = |actual: f64| -> Result<f64> {
        Ok(budget.absolute_pressure_tolerance_pa.max(finite(
            budget.relative_pressure_tolerance * target.max(actual))?))
    };
    if (low_drop - target).abs() <= tolerance(low_drop)? {
        return Ok(lower.copysign(pressure));
    }
    if (high_drop - target).abs() <= tolerance(high_drop)? {
        return Ok(upper.copysign(pressure));
    }
    // Start with the pressure secant; use the actual implicit-friction tangent
    // at every trial. A rejected Newton step bisects the admitted mass interval.
    let mut mass = lower + (upper - lower) * ((target - low_drop) / (high_drop - low_drop));
    let (mut residual, mut allowed) = (high_drop - target, tolerance(high_drop)?);
    for _ in 0..96 {
        let (actual, tangent) = evaluate(mass)?;
        residual = finite(actual - target)?;
        allowed = tolerance(actual)?;
        if residual.abs() <= allowed { return Ok(mass.copysign(pressure)); }
        if residual < 0.0 { lower = mass; } else { upper = mass; }
        let candidate = mass - residual / tangent;
        let next = if candidate > lower && candidate < upper { candidate }
            else { lower + 0.5 * (upper - lower) };
        if next == mass || next <= lower || next >= upper { break; }
        mass = next;
    }
    Err(SeriesPressureError::Closure { residual_pa: residual * pressure.signum(), tolerance_pa: allowed })
}
