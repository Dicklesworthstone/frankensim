//! Bind the accepted Churchill-Chu Robin rows to the complete native adjoint.
//! The primal already gates the card and converges h(T_mean, T_ambient).
//! Reconstruct its Nusselt value from the retained h; do not rerun a primal,
//! change the field, freeze buoyancy, or import a second correlation formula.

use super::{bad, finite, lower, poll, ConductionProblem, Cx, ProjectSpec,
    RungSolved, SolveRefusal};
use super::super::conjugate::AIR_THERMAL_CONDUCTIVITY_W_M_K;
use fs_conduction::{ScalarField, ThermalBc};
use fs_conduction::adjoint::{RobinGradient, RobinResponse};
use fs_project::ThermalBoundaryCondition;
use std::collections::BTreeMap;

pub(super) const CARD: &str = "convection.churchill-chu-vertical-plate";

/// At fixed pressure, geometry and frozen transport properties:
/// Ra is proportional to |Tw-Ta|/Tfilm^3 and Nu=(0.825+x)^2, x~Ra^(1/6).
/// Therefore d(log h)/d(log Ra)=(1-0.825/sqrt(Nu))/3. Include BOTH film
/// temperature derivatives: dTfilm/dTw=dTfilm/dTa=1/2. The ambient partial
/// here holds Tw fixed; its effect through Tw is already in the Jacobian.
/// Away from equilibrium, d(log|Tw-Ta|)/dTw = 1/(Tw-Ta) on BOTH branches.
/// Taking an absolute value in that denominator would reverse cooled-wall
/// feedback. The correlation's non-smooth equilibrium is not extrapolated.
pub(super) fn partials(length: f64, wall: f64, ambient: f64, h: f64)
    -> Result<[f64; 2], SolveRefusal>
{
    if ![length, wall, ambient, h].iter().all(|v| v.is_finite() && *v > 0.0)
        || wall == ambient { return Err(bad("natural adjoint needs a finite non-equilibrium wall and positive retained law values")); }
    let delta = finite(wall-ambient)?;
    let film = finite(ambient+0.5*delta)?;
    let nu = finite(h*length/AIR_THERMAL_CONDUCTIVITY_W_M_K)?;
    let exponent = finite((1.0-0.825/fs_math::det::sqrt(nu))/3.0)?;
    if !(exponent > 0.0 && exponent < 1.0/3.0) {
        return Err(bad("retained natural coefficient is outside the Churchill-Chu differential branch"));
    }
    Ok([finite(h*exponent*(1.0/delta-1.5/film))?,
        finite(h*exponent*(-1.0/delta-1.5/film))?])
}

pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64],
) -> Result<(RobinGradient, BTreeMap<String, f64>), SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing natural adjoint operator"))?;
    if !data.air_paths.is_empty() || data.radiating_boundary.is_some() {
        return Err(bad("natural adjoints do not freeze combined airflow or radiation feedback"));
    }
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing natural conduction setup"))?;
    let mut names = Vec::new();
    let mut slopes = Vec::new();
    let mut ambient_log_slopes = Vec::new();
    for declared in &setup.boundaries {
        let ThermalBoundaryCondition::NaturalConvection { characteristic_length,
            ambient_temperature, correlation } = &declared.condition else { continue; };
        if correlation != CARD { return Err(bad("natural adjoint requires the differentiated Churchill-Chu card")); }
        let region = data.boundary.region_names().iter().position(|name| name == &declared.target)
            .ok_or_else(|| bad("natural adjoint target has no retained boundary"))?;
        let ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(ambient) }
            = &data.boundary.conditions()[region] else { return Err(bad("natural adjoint needs its actual uniform Robin row")); };
        if *ambient != ambient_temperature.value { return Err(bad("natural reference differs from its declared ambient")); }
        let mut area = 0.0;
        let mut integral = 0.0;
        for (slot, face) in solved.mesh.boundary().iter().enumerate() {
            if slot % 512 == 0 { poll(cx)?; }
            if data.boundary.region_for(slot) != Some(region) { continue; }
            area = finite(area + face.area)?;
            for &v in &face.vertices {
                integral = finite(integral + (face.area/3.0)*solved.solution.temperature[v as usize])?;
            }
        }
        if area <= 0.0 { return Err(bad("natural adjoint target has no positive area")); }
        let [wall_slope, ambient_slope] = partials(characteristic_length.value,
            finite(integral/area)?, *ambient, *h)?;
        names.push(declared.target.as_str());
        slopes.push(wall_slope);
        ambient_log_slopes.push(finite(ambient_slope / h)?);
    }
    if names.is_empty() { return Err(bad("natural feedback adjoint has no declared law")); }
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory/32).unwrap_or(usize::MAX);
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let gradient = RobinResponse::pullback_mean_htc(cx, problem, data.interfaces.as_ref(),
        data.linear, &solved.solution.temperature, &names, &slopes, weights, entries).map_err(lower)?;
    let mut ambient = BTreeMap::new();
    for (i, name) in names.into_iter().enumerate() {
        ambient.insert(name.to_string(), finite(gradient.references[i]
            + gradient.log_htc[i]*ambient_log_slopes[i])?);
    }
    poll(cx)?;
    Ok((gradient, ambient))
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::natural;

    #[test]
    fn analytic_partials_match_the_actual_card_and_film_density_law() {
        let spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/heatsink-fan/heatsink-natural.fsim"
        ))).unwrap().decoded.spec;
        let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
        let law = natural::natural_laws(setup).unwrap().remove(0);
        let length = setup.boundaries.iter().find_map(|b| match &b.condition {
            ThermalBoundaryCondition::NaturalConvection { characteristic_length, .. } => Some(characteristic_length.value),
            _ => None,
        }).unwrap();
        for ambient in [290.0, 320.0] { for delta in [-80.0, -25.0, -1.0, 1.0, 25.0, 80.0] { for pressure in [80000.0, 101325.0] {
            let wall = ambient + delta;
            let evaluate = |tw: f64, ta: f64| {
                let mut input = law.clone(); input.ambient_k = ta;
                natural::coefficient(&input, tw-ta, pressure).unwrap().htc_w_m2_k
            };
            let [dw, da] = partials(length, wall, ambient, evaluate(wall, ambient)).unwrap();
            let step = 1e-3;
            let expected_w = (evaluate(wall+step, ambient)-evaluate(wall-step, ambient))/(2.0*step);
            let expected_a = (evaluate(wall, ambient+step)-evaluate(wall, ambient-step))/(2.0*step);
            for (actual, expected) in [(dw, expected_w), (da, expected_a)] {
                assert!((actual-expected).abs() < 2e-6*expected.abs().max(1e-3), "{actual:e} vs {expected:e}");
            }
            // Treating film density as fixed incorrectly cancels these terms.
            assert!(dw+da < 0.0);
        } } }
    }

    #[test]
    fn a_natural_derivative_never_extrapolates_equilibrium_or_nonphysical_state() {
        for values in [[0.06, 0.0, 300.0, 5.0], [0.06, 300.0, 300.0, 5.0],
            [0.06, 320.0, 300.0, f64::NAN], [0.0, 320.0, 300.0, 5.0],
            [0.06, 320.0, 300.0, 0.01]] {
            assert!(partials(values[0],values[1],values[2],values[3]).is_err());
        }
    }
}
