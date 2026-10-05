//! Keep the exact card-backed patches from the accepted producer, and lower
//! the complete natural/convection/radiation derivative on its unchanged field.
use super::{bad, finite, lower, natural_feedback, poll, row, zeros, ConductionProblem,
    Cx, ProjectSpec, RungSolved, SolveRefusal};
use fs_conduction::{ScalarField, ThermalBc};
use fs_conduction::adjoint::RobinGradient;
use fs_conduction::radiation::pullback_ambient_radiation;
use fs_project::ThermalBoundaryCondition as B;
use std::collections::BTreeMap;
use super::super::natural;

pub(super) struct Pullback {
    pub gradient: RobinGradient,
    pub natural_ambient: BTreeMap<String, f64>,
    pub radiation_rows: Vec<String>,
}


pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64],
) -> Result<Pullback, SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing radiative adjoint operator"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing radiative conduction setup"))?;
    let declaration = setup.radiation.as_ref().ok_or_else(|| bad("missing physical radiation declaration"))?;
    let combined = data.radiating_boundary.as_ref()
        .ok_or_else(|| bad("missing retained combined radiation boundary"))?;
    if !data.air_paths.is_empty() {
        return Err(bad("radiation adjoints do not freeze airflow feedback"));
    }
    if declaration.surfaces.is_empty() || declaration.surfaces.len() > 64 {
        return Err(bad("radiative adjoint needs one to 64 declared surfaces"));
    }
    if declaration.surfaces.len() != data.radiation_patches.len() {
        return Err(bad("radiation adjoint did not retain every declared emissivity patch"));
    }
    if combined.region_names() != data.boundary.region_names() {
        return Err(bad("combined radiation boundary changed region ownership"));
    }
    let mut areas = zeros(data.boundary.region_names().len())?;
    let mut integrals = zeros(areas.len())?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let region = data.boundary.region_for(slot);
        if combined.region_for(slot) != region {
            return Err(bad("combined radiation boundary changed face ownership"));
        }
        let Some(region) = region else { continue; };
        areas[region] = finite(areas[region] + face.area)?;
        for &v in &face.vertices {
            integrals[region] = finite(integrals[region]
                + (face.area/3.0)*solved.solution.temperature[v as usize])?;
        }
    }
    let natural_laws = natural::natural_laws(setup)?;
    let pressure = spec.envelope.as_ref().ok_or_else(|| bad("missing pressure envelope"))?.pressure.value;
    let mut replacements = Vec::new();
    let mut names = Vec::new();
    let mut slopes = Vec::new();
    let mut ambient_log_slopes = Vec::new();
    for declared in &setup.boundaries {
        poll(cx)?;
        let is_natural = matches!(declared.condition, B::NaturalConvection { .. });
        let is_radiating = data.radiation_patches.iter().any(|p| p.region() == declared.target);
        // Ordinary nonradiating controls are already contracted from the same
        // corrected nodal dual. Do not charge them for unused feedback factors.
        if !is_natural && !is_radiating { continue; }
        match &declared.condition {
            B::Convection { .. } | B::NaturalConvection { .. } => {},
            B::AirflowConvection { .. } => return Err(bad("radiative adjoints do not freeze air feedback")),
            _ => continue,
        }
        let region = data.boundary.region_names().iter().position(|name| name == &declared.target)
            .ok_or_else(|| bad("radiative adjoint control has no original boundary trace"))?;
        let ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(reference) }
            = &data.boundary.conditions()[region]
        else { return Err(bad("radiative adjoints require original uniform convective rows")); };
        let (wall_slope, ambient_log_slope) = match &declared.condition {
            B::NaturalConvection { characteristic_length, ambient_temperature, correlation } => {
                if correlation != natural_feedback::CARD || *reference != ambient_temperature.value {
                    return Err(bad("natural/radiative adjoint must retain the declared Churchill-Chu law and ambient"));
                }
                if areas[region] <= 0.0 { return Err(bad("natural/radiative adjoint target has no area")); }
                let mean = finite(integrals[region]/areas[region])?;
                let law = natural_laws.iter().find(|law| law.target == declared.target)
                    .ok_or_else(|| bad("missing natural convection card"))?;
                let h = natural::coefficient(law, mean-reference, pressure)?.htc_w_m2_k;
                let [wall, ambient] = natural_feedback::partials(characteristic_length.value,
                    mean, *reference, h)?;
                replacements.push((region, h, *reference));
                (wall, Some(finite(ambient/h)?))
            }
            _ => (0.0, None),
        };
        names.push(declared.target.as_str());
        slopes.push(wall_slope);
        ambient_log_slopes.push(ambient_log_slope);
    }
    // Check the immutable constitutive source, not a coefficient inferred from
    // a rounded JSON receipt or from subtracting two nearly equal Robin rows.
    for (declared, patch) in declaration.surfaces.iter().zip(&data.radiation_patches) {
        poll(cx)?;
        if declared.target != patch.region()
            || declared.reservoir_temperature.value != patch.ambient_temperature_k()
            || declared.card != patch.emissivity().card_identity().to_hex()
            || declared.query_temperature.value != patch.emissivity().temperature_k()
        {
            return Err(bad("retained radiation law differs from its declared physical source"));
        }
    }
    // Re-evaluate natural h at the actual retained wall, then let the core
    // owner re-evaluate radiation and check the COMPLETE physical residual.
    let boundary = data.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory/32).unwrap_or(usize::MAX);
    let full = pullback_ambient_radiation(cx, problem, data.interfaces.as_ref(), data.linear,
        &solved.solution.temperature, &names, &slopes, &data.radiation_patches, weights, entries)
        .map_err(|error| {
            let mut refusal = lower(error);
            refusal.fix = "tighten the declared radiation temperature/heat tolerances and thermal solve tolerance; the actual nonlinear residual and complete adjoint must both converge".into();
            refusal
        })?;
    let mut ambient = BTreeMap::new();
    for (i, name) in names.into_iter().enumerate() {
        if let Some(slope) = ambient_log_slopes[i] {
            ambient.insert(name.to_string(), finite(full.convection.references[i]
                + full.convection.log_htc[i]*slope)?);
        }
    }
    let mut rows = Vec::with_capacity(2*declaration.surfaces.len());
    for (i, surface) in declaration.surfaces.iter().enumerate() {
        poll(cx)?;
        rows.push(row("radiation-reservoir-temperature", &surface.name, i, "K",
            full.reservoir_temperatures[i])?);
        rows.push(row("radiation-emissivity", &surface.name, i, "1", full.emissivities[i])?);
    }
    poll(cx)?;
    Ok(Pullback { gradient: full.convection, natural_ambient: ambient, radiation_rows: rows })
}

#[cfg(test)]
use fs_conduction::STEFAN_BOLTZMANN_W_M2_K4;

// h_rad = epsilon sigma (Tw + Tr)(Tw^2 + Tr^2), exactly the producer's
// secant model. Its derivatives stay defined at equal temperatures and when
// the reservoir heats the wall. This is NOT the pointwise 4 epsilon sigma T^3
// boundary Jacobian; the nonuniform trace and the mean feedback both matter.
#[cfg(test)]
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
#[cfg(test)]
fn chain(reference_bar: f64, log_h_bar: f64, total_h: f64, reference: f64,
    component_h: f64, component_t: f64, dh: f64, dt: f64) -> Result<f64, SolveRefusal>
{
    finite(reference_bar*(component_h/total_h)*dt
        + (log_h_bar+reference_bar*(component_t-reference))*(dh/total_h))
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
