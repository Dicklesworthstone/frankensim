//! Admit the actual P1 source strain through the original T/rate source tables.
//!
//! Each affine velocity field has one constant symmetric strain tensor, while
//! its P1 enthalpy spans a convex interval. This is source admission only;
//! receiving collision and heat use the original native constitutive owner.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn admit_source(
    source: TetrahedralLiquidEnthalpyState<'_>,
    density: f64,
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    rate_axis: &str,
    limits: TetEnthalpyTransferLimits,
    queries: &mut usize,
    buoyancy: Option<BoussinesqConfig<3>>,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<(f64, f64, Vec<f64>), TetEnthalpyTransferError> {
    if source.motion.positions_m.len() < 4
        || source.motion.positions_m.len() != source.motion.velocities_m_s.len()
    {
        return Err(
            invalid("P1 rheology needs the actual matching geometry and velocity fields").into(),
        );
    }
    let mut minimum_h = f64::INFINITY;
    let mut maximum_h = f64::NEG_INFINITY;
    let mut rates = Vec::with_capacity(source.motion.tetrahedra.len());
    for indices in source.motion.tetrahedra {
        checkpoint(cx)?;
        if gate.is_requested() {
            return Err(TetLiquidTransferError::Cancelled.into());
        }
        if indices
            .iter()
            .any(|&i| i as usize >= source.motion.positions_m.len())
        {
            return Err(invalid("tetrahedron index lies outside the P1 rheology field").into());
        }
        let p = indices.map(|i| source.motion.positions_m[i as usize]);
        let v = indices.map(|i| source.motion.velocities_m_s[i as usize]);
        if p.iter()
            .flatten()
            .chain(v.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("P1 rheology requires finite original positions and velocities").into(),
            );
        }
        let gradient = super::super::field_gradient(p, v)?;
        // gamma_dot = sqrt(2 D:D), D = (grad(v) + grad(v)^T)/2.
        // A rigid translation/rotation contributes no constitutive rate.
        let mut rate = 0.0_f64;
        for a in 0..3 {
            for b in 0..3 {
                let strain = 0.5 * gradient[a][b] + 0.5 * gradient[b][a];
                rate = rate.hypot(core::f64::consts::SQRT_2 * strain);
            }
        }
        if !rate.is_finite() {
            return Err(invalid("unrepresentable physical P1 shear rate").into());
        }
        let h = indices.map(|i| source.specific_enthalpy_j_kg[i as usize]);
        let lower = h.into_iter().fold(f64::INFINITY, f64::min);
        let upper = h.into_iter().fold(f64::NEG_INFINITY, f64::max);
        let coordinate = Some((rate_axis, rate));
        // Shared vertices can have distinct material rates on adjacent tets,
        // so they must be re-admitted at each actual constitutive coordinate.
        for specific in h {
            checkpoint(cx)?;
            admit_source_sample(specific, density, law, caloric, limits, queries, coordinate, buoyancy)?;
        }
        for knot in caloric.knots() {
            checkpoint(cx)?;
            if knot.specific_enthalpy_j_kg >= lower && knot.specific_enthalpy_j_kg <= upper {
                admit_source_sample(
                    knot.specific_enthalpy_j_kg,
                    density,
                    law,
                    caloric,
                    limits,
                    queries,
                    coordinate,
                    buoyancy,
                )?;
            }
        }
        admit_source_lines(
            lower, upper, density, law, caloric, limits, queries, coordinate, buoyancy, cx,
        )?;
        minimum_h = minimum_h.min(lower);
        maximum_h = maximum_h.max(upper);
        rates.push(rate);
    }
    Ok((minimum_h, maximum_h, rates))
}
