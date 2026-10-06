//! Re-evaluate existing physical boundary laws on the unchanged native field.
//! No substitute field, finite-difference coefficient, or report-JSON physics.
use super::{BTreeSet, Cx, ProjectSpec, RungSolved, ScalarField, SolveRefusal,
    ThermalBc, ThermalBoundary, bad, finite, lower, poll, zeros};
use super::super::natural_feedback;
use super::super::super::natural;
use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};

pub(super) struct Law {
    pub boundary: ThermalBoundary,
    pub names: Vec<String>,
    pub h_slopes: Vec<f64>,
    pub reference_slopes: Vec<f64>,
    pub feedback: Vec<f64>,
}

pub(super) fn bind(cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved,
    max_entries: usize) -> Result<Law, SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("prescribed lift lost its operator"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("prescribed lift has no conduction declarations"))?;
    let mut regions = Vec::new();
    let mut points = Vec::new();
    for (region, bc) in data.boundary.conditions().iter().enumerate() {
        poll(cx)?;
        match bc {
            ThermalBc::Robin {htc:ScalarField::Uniform(h),t_ref:ScalarField::Uniform(r)} => {
                regions.push(region); points.push([*h,*r,0.0,0.0]);
            }
            ThermalBc::Robin {..} => return Err(bad("native lift requires uniform Robin laws")),
            _ => {},
        }
    }
    let m = regions.len();
    if m > 64 || m*m > max_entries || solved.mesh.vertex_count() > max_entries {
        return Err(bad("prescribed boundary feedback exceeds its trace/entry budget"));
    }
    let names: Vec<String> = regions.iter().map(|&i| data.boundary.region_names()[i].clone()).collect();
    let mut areas = zeros(m)?;
    let mut means = zeros(m)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let region = data.boundary.region_for(slot);
        if data.radiating_boundary.as_ref().is_some_and(|b|
            b.region_names() != data.boundary.region_names() || b.region_for(slot) != region) {
            return Err(bad("radiative lift changed the original boundary ownership"));
        }
        let Some(i) = region.and_then(|r| regions.iter().position(|&j| j == r)) else { continue; };
        areas[i] = finite(areas[i]+face.area)?;
        for &v in &face.vertices { means[i] = finite(means[i]
            + (face.area/3.0)*solved.solution.temperature[v as usize])?; }
    }
    for i in 0..m {
        if areas[i] <= 0.0 { return Err(bad("prescribed lift has an empty Robin trace")); }
        means[i] = finite(means[i]/areas[i])?;
    }
    let natural_laws = natural::natural_laws(setup)?;
    for law in &natural_laws {
        poll(cx)?;
        let i = names.iter().position(|name| name == &law.target)
            .ok_or_else(|| bad("natural lift has no original Robin trace"))?;
        let declared = setup.boundaries.iter().find(|b| b.target == law.target)
            .ok_or_else(|| bad("natural lift lost its declared law"))?;
        let fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length,ambient_temperature,correlation} = &declared.condition
            else { return Err(bad("natural lift law differs from the declaration")); };
        if correlation != natural_feedback::CARD || ambient_temperature.value != points[i][1] {
            return Err(bad("natural lift requires the original Churchill-Chu law and ambient"));
        }
        let pressure = spec.envelope.as_ref().ok_or_else(|| bad("natural lift has no pressure envelope"))?.pressure.value;
        let h = natural::coefficient(law,means[i]-points[i][1],pressure)?.htc_w_m2_k;
        let [slope,_] = natural_feedback::partials(characteristic_length.value,means[i],points[i][1],h)?;
        points[i][0] = h; points[i][2] = slope;
    }
    let mut feedback = zeros(m*m)?;
    if !data.air_paths.is_empty() {
        if !natural_laws.is_empty() { return Err(bad("combined natural/airflow lift is not supplied")); }
        let law = affine_reference_law(cx,&data.air_paths,64,max_entries).map_err(|e| match e {
            CoupledGoalError::Interrupted => super::super::conduction_error("cli-solve-cancelled",
                "prescribed air lift interrupted","resume the accepted pipeline prefix"),
            e => bad(format!("prescribed air-reference law refused: {e}")),
        })?;
        let indices: Vec<_> = law.regions().iter().map(|name| names.iter().position(|n| n == name)
            .ok_or_else(|| bad("prescribed air lift lost a physical trace"))).collect::<Result<_,_>>()?;
        let mut start = 0;
        for path in &data.air_paths {
            poll(cx)?;
            let end = start+path.segments().len();
            let selected = indices.get(start..end).ok_or_else(|| bad("prescribed air lift lost segment ordering"))?;
            let walls: Vec<_> = selected.iter().map(|&i| means[i]).collect();
            let march = path.march(&walls).map_err(|e| bad(format!("prescribed air lift march refused: {e}")))?;
            for ((segment,state),&i) in path.segments().iter().zip(&march.segments).zip(selected) {
                if segment.region() != names[i] || segment.htc_w_per_m2_k() != points[i][0]
                    || (areas[i]-segment.area_m2()).abs() > 128.0*f64::EPSILON*areas[i].max(segment.area_m2()) {
                    return Err(bad("prescribed air lift differs from the retained coefficient/area/region"));
                }
                points[i][1] = state.reference_temperature_k;
            }
            start = end;
        }
        if start != indices.len() { return Err(bad("prescribed air lift omitted segments")); }
        for (i,&a) in indices.iter().enumerate() { for (j,&b) in indices.iter().enumerate() {
            feedback[a*m+b] = law.wall_matrix()[i*indices.len()+j];
        } }
    }
    let surfaces = setup.radiation.as_ref().map_or(&[][..],|r| r.surfaces.as_slice());
    if surfaces.len() != data.radiation_patches.len() || surfaces.len() > 64 {
        return Err(bad("prescribed lift lost original radiation patches"));
    }
    let mut seen = BTreeSet::new();
    for (surface,patch) in surfaces.iter().zip(&data.radiation_patches) {
        poll(cx)?;
        if surface.target != patch.region() || !seen.insert(patch.region())
            || surface.reservoir_temperature.value != patch.ambient_temperature_k()
            || surface.card != patch.emissivity().card_identity().to_hex()
            || surface.query_temperature.value != patch.emissivity().temperature_k() {
            return Err(bad("prescribed radiative lift changed its original card or physical law"));
        }
        let i = names.iter().position(|name| name == patch.region())
            .ok_or_else(|| bad("prescribed radiative lift has no original convection trace"))?;
        let [hc,tc,hcs,_] = points[i];
        let hr = patch.secant_coefficient_w_m2_k(means[i]).map_err(lower)?;
        let [hrs,_] = patch.secant_partials_w_m2_k2(means[i]).map_err(lower)?;
        let tr = patch.ambient_temperature_k();
        let h = finite(hc+hr)?;
        if h <= 0.0 { return Err(bad("prescribed lift has a nonpositive combined coefficient")); }
        let reference = finite((hc/h)*tc+(hr/h)*tr)?;
        points[i] = [h,reference,finite(hcs+hrs)?,finite((hcs/h)*(tc-reference)+(hrs/h)*(tr-reference))?];
        for entry in &mut feedback[i*m..(i+1)*m] { *entry = finite((hc/h)*(*entry))?; }
    }
    let replacements: Vec<_> = regions.iter().zip(&points).map(|(&r,p)| (r,p[0],p[1])).collect();
    let boundary = data.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    poll(cx)?;
    Ok(Law {boundary,names,h_slopes:points.iter().map(|p| p[2]).collect(),
        reference_slopes:points.iter().map(|p| p[3]).collect(),feedback})
}
