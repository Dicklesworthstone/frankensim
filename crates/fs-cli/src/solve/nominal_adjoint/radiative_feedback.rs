//! Adjoint of the actual consistent-trace ambient-radiation fixed point.
//!
//! The completed in-memory radiation receipt retains the exact selected
//! emissivity. Use it, not an inferred bulk-material value. Re-evaluate the
//! physical laws at the FINAL field and require their complete primal residual
//! to converge. Differentiating the last frozen outer iterate is insufficient.
//! No additional primal, field replacement, or material-card re-selection.

use std::collections::BTreeMap;
use fs_conduction::{ScalarField, ThermalBc, STEFAN_BOLTZMANN_W_M2_K4};
use fs_conduction::adjoint::{RobinGradient, RobinResponse};
use fs_project::ThermalBoundaryCondition as B;
use crate::json_read::JsonValue;
use super::{bad, finite, lower, natural_feedback, poll, row, zeros,
    ConductionProblem, Cx, ProjectSpec, RungSolved, SolveRefusal};
use super::super::natural;

pub(super) struct Pullback {
    pub gradient: RobinGradient,
    pub natural_ambient: BTreeMap<String, f64>,
    pub radiation_rows: Vec<String>,
}

// h_rad = epsilon sigma (Tw + Tr)(Tw^2 + Tr^2), exactly the producer's
// secant model. Its derivatives stay defined at equal temperatures and when
// the reservoir heats the wall. This is NOT the pointwise 4 epsilon sigma T^3
// boundary Jacobian; the nonuniform trace and the mean feedback both matter.
fn secant(epsilon: f64, wall: f64, reservoir: f64) -> Result<[f64; 3], SolveRefusal> {
    if !(epsilon.is_finite() && epsilon > 0.0 && epsilon <= 1.0
        && wall.is_finite() && wall > 0.0 && reservoir.is_finite() && reservoir > 0.0) {
        return Err(bad("radiation differential requires admitted emissivity and positive finite temperatures"));
    }
    let c = epsilon * STEFAN_BOLTZMANN_W_M2_K4;
    Ok([finite(c*(wall+reservoir)*wall.mul_add(wall, reservoir*reservoir))?,
        finite(c*(3.0*wall*wall+2.0*wall*reservoir+reservoir*reservoir))?,
        finite(c*(wall*wall+2.0*wall*reservoir+3.0*reservoir*reservoir))?])
}

// Compose one physical component's partials at FIXED wall temperature with
// the independent effective-reference/log(h) partials returned by the owner.
fn chain(reference_bar: f64, log_h_bar: f64, total_h: f64, reference: f64,
    component_h: f64, component_t: f64, dh: f64, dt: f64) -> Result<f64, SolveRefusal>
{
    finite(reference_bar*(component_h/total_h)*dt
        + (log_h_bar+reference_bar*(component_t-reference))*(dh/total_h))
}

fn number(value: &JsonValue, name: &str) -> Result<f64, SolveRefusal> {
    finite(value.f64_field(name).ok_or_else(|| bad(format!("retained radiation is missing {name}")))?)
}

pub(super) fn pullback(cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64])
    -> Result<Pullback, SolveRefusal>
{
    poll(cx)?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing radiative adjoint operator"))?;
    if !data.air_paths.is_empty() { return Err(bad("the complete radiation/airflow adjoint is not supplied")); }
    let combined = data.radiating_boundary.as_ref().ok_or_else(|| bad("missing retained combined radiation boundary"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing radiative conduction setup"))?;
    let declaration = setup.radiation.as_ref().ok_or_else(|| bad("missing radiation declaration"))?;
    if declaration.surfaces.is_empty() || declaration.surfaces.len() > 64 {
        return Err(bad("radiative adjoint needs one to 64 declared surfaces"));
    }
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let text = solved.radiation_fragment.as_deref().ok_or_else(|| bad("missing completed radiation receipt"))?;
    if text.len() as u64 > memory/4 { return Err(bad("retained radiation exceeds the adjoint read budget")); }
    let receipt = JsonValue::parse(text).map_err(|_| bad("invalid retained radiation receipt"))?;
    if receipt.str_field("model") != Some("area-mean-gray-surface-to-black-reservoir") {
        return Err(bad("unsupported retained radiation model"));
    }
    let recorded = receipt.get("surfaces").and_then(JsonValue::as_array)
        .ok_or_else(|| bad("missing retained radiation surfaces"))?;
    if recorded.len() != declaration.surfaces.len() { return Err(bad("radiation surface count changed")); }
    let mut radiation = BTreeMap::new();
    for (ordinal, surface) in declaration.surfaces.iter().enumerate() {
        poll(cx)?;
        let matches: Vec<_> = recorded.iter().filter(|r| r.str_field("target") == Some(surface.target.as_str())).collect();
        if matches.len() != 1 { return Err(bad("missing or duplicate retained radiating target")); }
        let record = matches[0];
        if record.str_field("name") != Some(surface.name.as_str())
            || record.str_field("card") != Some(surface.card.as_str())
            || number(record, "query_temperature_k")? != surface.query_temperature.value
            || number(record, "ambient_temperature_k")? != surface.reservoir_temperature.value
            || radiation.insert(surface.target.as_str(), (ordinal, surface, record)).is_some() {
            return Err(bad("radiation material, reservoir or surface identity changed"));
        }
    }
    let nr = data.boundary.region_names().len();
    if combined.region_names() != data.boundary.region_names() {
        return Err(bad("combined radiation boundary changed region ownership"));
    }
    let mut areas = zeros(nr)?;
    let mut integrals = zeros(nr)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let region = data.boundary.region_for(slot);
        if combined.region_for(slot) != region { return Err(bad("combined radiation boundary changed face ownership")); }
        let Some(region) = region else { continue; };
        areas[region] = finite(areas[region]+face.area)?;
        for &v in &face.vertices {
            integrals[region] = finite(integrals[region]+face.area/3.0*solved.solution.temperature[v as usize])?;
        }
    }
    let natural_laws = natural::natural_laws(setup)?;
    let pressure = spec.envelope.as_ref().ok_or_else(|| bad("missing pressure envelope"))?.pressure.value;
    let mut names = Vec::new();
    let mut points = Vec::new();
    let mut controls = Vec::new();
    let mut radiation_count = 0;
    for boundary in &setup.boundaries {
        poll(cx)?;
        let natural = matches!(boundary.condition, B::NaturalConvection { .. });
        let radiating = radiation.get(boundary.target.as_str());
        if !natural && radiating.is_none() { continue; }
        let region = data.boundary.region_names().iter().position(|r| r == &boundary.target)
            .ok_or_else(|| bad("radiating or natural target has no native trace"))?;
        let ThermalBc::Robin { htc: ScalarField::Uniform(retained_h), t_ref: ScalarField::Uniform(ambient) }
            = &data.boundary.conditions()[region] else { return Err(bad("radiation needs uniform original Robin laws")); };
        if areas[region] <= 0.0 { return Err(bad("radiating or natural target has zero area")); }
        let wall = finite(integrals[region]/areas[region])?;
        let (hc, hc_wall, hc_ambient) = if let B::NaturalConvection { characteristic_length,
            ambient_temperature, correlation } = &boundary.condition {
            if correlation != natural_feedback::CARD || *ambient != ambient_temperature.value {
                return Err(bad("natural/radiation adjoint has a different or unsupported native law"));
            }
            let law = natural_laws.iter().find(|law| law.target == boundary.target)
                .ok_or_else(|| bad("missing natural convection card"))?;
            let h = natural::coefficient(law, wall-ambient, pressure)?.htc_w_m2_k;
            let [dw, da] = natural_feedback::partials(characteristic_length.value, wall, *ambient, h)?;
            (h, dw, da)
        } else { (*retained_h, 0.0, 0.0) };
        let (hr, hr_wall, hr_reservoir, reservoir, epsilon) = if let Some((_, surface, record)) = radiating {
            radiation_count += 1;
            if (number(record, "mean_temperature_k")?-wall).abs() > 1e-9*wall.abs().max(1.0) {
                return Err(bad("radiation receipt is not bound to the retained wall field"));
            }
            let epsilon = number(record, "emissivity")?;
            let reservoir = surface.reservoir_temperature.value;
            let [h, dw, dr] = secant(epsilon, wall, reservoir)?;
            (h, dw, dr, reservoir, epsilon)
        } else { (0.0, 0.0, 0.0, *ambient, 1.0) };
        let h = finite(hc+hr)?;
        let reference = finite(hc.mul_add(*ambient, hr*reservoir)/h)?;
        let dh = finite(hc_wall+hr_wall)?;
        let dr = finite((hc_wall*(ambient-reference)+hr_wall*(reservoir-reference))/h)?;
        names.push(boundary.target.as_str());
        points.push([h, reference, dh, dr]);
        controls.push((natural, hc, hc_ambient, *ambient, hr, hr_reservoir, reservoir, epsilon, radiating.copied()));
    }
    if radiation_count != declaration.surfaces.len() { return Err(bad("a declared radiation surface was not differentiated")); }
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    // Rebinding preserves every nonselected law and prescribed vertex. It
    // checks the FULL constitutive residual, rather than the accepted secant
    // iterate's residual. Loose outer convergence refuses with no fallback.
    let gradient = RobinResponse::pullback_mean_robin_at(cx, problem, data.interfaces.as_ref(),
        data.linear, &solved.solution.temperature, &names, &points, weights,
        usize::try_from(memory/32).unwrap_or(usize::MAX)).map_err(|error| {
            let mut refusal = lower(error);
            refusal.fix = "tighten the declared radiation temperature/heat tolerances and thermal solve tolerance; the actual nonlinear residual and complete adjoint must both converge".into();
            refusal
        })?;
    let mut natural_ambient = BTreeMap::new();
    let mut radiation_rows = Vec::new();
    for (i, (natural, hc, hc_ambient, ambient, hr, hr_reservoir, reservoir, epsilon, radiation)) in controls.into_iter().enumerate() {
        poll(cx)?;
        let [h, reference, _, _] = points[i];
        let r = gradient.references[i];
        let l = gradient.log_htc[i];
        if natural {
            natural_ambient.insert(names[i].to_string(), chain(r, l, h, reference, hc, ambient, hc_ambient, 1.0)?);
        }
        if let Some((ordinal, surface, _)) = radiation {
            radiation_rows.push(row("radiation-reservoir-temperature", &surface.name, ordinal, "K",
                chain(r, l, h, reference, hr, reservoir, hr_reservoir, 1.0)?)?);
            radiation_rows.push(row("radiation-emissivity", &surface.name, ordinal, "1",
                chain(r, l, h, reference, hr, reservoir, finite(hr/epsilon)?, 0.0)?)?);
        }
    }
    Ok(Pullback { gradient, natural_ambient, radiation_rows })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn radiative_secant_partials_include_equal_and_hot_reservoirs() {
        for wall in [280.0, 300.0, 400.0] { for reservoir in [280.0, 300.0, 450.0] {
            let [h, dw, dr] = secant(0.73, wall, reservoir).unwrap();
            let e = 1e-3;
            let fdw = (secant(0.73,wall+e,reservoir).unwrap()[0]-secant(0.73,wall-e,reservoir).unwrap()[0])/(2.0*e);
            let fdr = (secant(0.73,wall,reservoir+e).unwrap()[0]-secant(0.73,wall,reservoir-e).unwrap()[0])/(2.0*e);
            assert!((dw-fdw).abs() < 1e-8*dw.abs());
            assert!((dr-fdr).abs() < 1e-8*dr.abs());
            // Derivative of the integrated heat law agrees only AFTER the
            // derivative of the secant and of (Tw-Tr) are BOTH included.
            let scale = 4.0*0.73*STEFAN_BOLTZMANN_W_M2_K4;
            assert!((h+dw*(wall-reservoir)-scale*wall.powi(3)).abs() < 1e-11);
            assert!((-h+dr*(wall-reservoir)+scale*reservoir.powi(3)).abs() < 1e-11);
        } }
    }

    #[test]
    fn weighted_reference_control_chain_matches_direct_component_changes() {
        let (hc, ta, tr, wall, eps) = (7.0, 295.0, 360.0, 330.0, 0.8);
        let [hr, _, dhr] = secant(eps, wall, tr).unwrap();
        let h = hc+hr;
        let reference = (hc*ta+hr*tr)/h;
        let (rb, lb) = (0.6, -1.7);
        let objective = |h: f64, r: f64| rb*r+lb*h.ln();
        let e = 1e-3;
        let eval = |a: f64| {
            let rad = secant(eps, wall, a).unwrap()[0];
            objective(hc+rad,(hc*ta+rad*a)/(hc+rad))
        };
        let actual = chain(rb,lb,h,reference,hr,tr,dhr,1.0).unwrap();
        assert!((actual-(eval(tr+e)-eval(tr-e))/(2.0*e)).abs() < 1e-8);
        assert!((actual-rb*hr/h).abs() > 1e-3, "freezing the secant loses a material control term");
        for values in [[0.0,300.0,300.0],[1.1,300.0,300.0],[0.8,0.0,300.0],[0.8,f64::NAN,300.0]] {
            assert!(secant(values[0],values[1],values[2]).is_err());
        }
    }
}
