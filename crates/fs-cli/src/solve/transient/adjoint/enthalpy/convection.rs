//! Physical convection at an unchanged accepted enthalpy endpoint.
//!
//! The natural card supplies h and its wall derivative. The actual air march
//! supplies reference temperatures; its analytic affine map supplies every
//! upstream wall derivative. Radiation remains with its independent owner.

use std::collections::BTreeSet;

use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};
use fs_conduction::transient::enthalpy::EnthalpyEndpointEvaluation;
use fs_conduction::{ScalarField, ThermalBc, ThermalBoundary};

use super::{
    ConductionProblem, Cx, Instant, SolveRefusal, airflow, bad, buoyancy, finite, lower, poll,
    zeros,
};

mod controls;
pub(super) use controls::Controls;
pub(in crate::solve::transient) use controls::storage_bytes as control_storage_bytes;

pub(super) struct Bound {
    pub(super) boundary: ThermalBoundary,
    pub(super) names: Vec<String>,
    pub(super) references: Vec<Vec<f64>>,
    pub(super) wall_slopes: Vec<f64>,
    ambient_slopes: Vec<f64>,
    air_regions: BTreeSet<String>,
    air_heat_w: Option<f64>,
}

pub(super) fn ports(
    natural: &[buoyancy::NaturalLaw],
    air: Option<&airflow::ConjugatePath>,
) -> Result<usize, SolveRefusal> {
    natural
        .len()
        .checked_add(air.map_or(0, |path| path.segments.len()))
        .filter(|&count| count <= 64)
        .ok_or_else(|| bad("enthalpy convection adjoints admit at most 64 distinct physical ports"))
}

/// Two rank-one terms per port at most, each with two nodal factors.
pub(super) fn feedback_entries(n: usize, ports: usize) -> Result<usize, SolveRefusal> {
    if n == 0 || ports > 64 {
        return Err(bad(
            "enthalpy convection feedback needs a nonempty mesh and at most 64 ports",
        ));
    }
    n.checked_mul(4)
        .and_then(|v| v.checked_mul(ports))
        .filter(|&v| v.checked_mul(8).is_some_and(|b| isize::try_from(b).is_ok()))
        .ok_or_else(|| bad("enthalpy convection feedback allowance overflows"))
}

/// One active boundary reconstruction/physical evaluation and complete
/// low-rank tangent, with bounded port matrices and trace metadata. No nodal
/// dense matrix or per-time-step convection matrices are retained.
pub(in crate::solve::transient) fn storage_bytes(
    n: usize,
    elements: usize,
    ports: usize,
) -> Result<usize, SolveRefusal> {
    if ports == 0 {
        return Ok(0);
    }
    feedback_entries(n, ports)?
        .checked_mul(24)
        .and_then(|v| {
            ports
                .checked_mul(ports)
                .and_then(|p| p.checked_mul(48))
                .and_then(|b| v.checked_add(b))
        })
        .and_then(|v| n.checked_mul(128).and_then(|b| v.checked_add(b)))
        .and_then(|v| elements.checked_mul(2048).and_then(|b| v.checked_add(b)))
        .and_then(|v| ports.checked_mul(8192).and_then(|b| v.checked_add(b)))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("enthalpy convection adjoint workspace allowance overflows"))
}

fn air_error(error: CoupledGoalError) -> SolveRefusal {
    match error {
        CoupledGoalError::Interrupted => crate::solve::conduction_error(
            "cli-solve-cancelled",
            "enthalpy air-reference differentiation was cancelled",
            "resume the accepted pipeline prefix",
        ),
        other => bad(format!(
            "enthalpy air-reference derivative refused: {other}"
        )),
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    temperature: &[f64],
    natural: &[buoyancy::NaturalLaw],
    air: Option<&airflow::ConjugatePath>,
    pressure_pa: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<Option<Bound>, SolveRefusal> {
    poll(cx, deadline)?;
    let count = ports(natural, air)?;
    if count == 0 {
        return Ok(None);
    }
    if temperature.len() != problem.mesh.vertex_count() {
        return Err(bad(
            "enthalpy convection derivative endpoint has the wrong shape",
        ));
    }
    let paths = air.map(|path| path.air_paths()).unwrap_or_default();
    let mut names = Vec::new();
    names
        .try_reserve_exact(count)
        .map_err(|_| bad("enthalpy convection port allocation refused"))?;
    names.extend(natural.iter().map(|law| law.target.clone()));
    names.extend(
        paths
            .iter()
            .flat_map(|path| path.segments())
            .map(|segment| segment.region().to_string()),
    );
    if names.len() != count || names.iter().collect::<BTreeSet<_>>().len() != count {
        return Err(bad(
            "natural and forced convection need distinct, complete port ownership",
        ));
    }
    let regions = names
        .iter()
        .map(|name| {
            let region = problem
                .boundary
                .region_names()
                .iter()
                .position(|found| found == name)
                .ok_or_else(|| bad("physical convection has no original boundary region"))?;
            match &problem.boundary.conditions()[region] {
                ThermalBc::Robin {
                    htc: ScalarField::Uniform(h),
                    t_ref: ScalarField::Uniform(t),
                } if h.is_finite() && *h > 0.0 && t.is_finite() && *t > 0.0 => Ok(region),
                _ => Err(bad(
                    "enthalpy feedback requires original positive uniform convection rows",
                )),
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut areas = zeros(count)?;
    let mut means = zeros(count)?;
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 {
            poll(cx, deadline)?;
        }
        let Some(index) = problem
            .boundary
            .region_for(slot)
            .and_then(|region| regions.iter().position(|&found| found == region))
        else {
            continue;
        };
        // Reproduce the forward Robin report's P1 mean in its original
        // face/vertex order. A differently rounded trace can change the
        // physical natural coefficient or air march during revalidation.
        let mut wall_sum = 0.0;
        for &vertex in &face.vertices {
            wall_sum = finite(wall_sum + finite(temperature[vertex as usize])?)?;
        }
        areas[index] = finite(areas[index] + face.area)?;
        means[index] = finite(face.area.mul_add(wall_sum / 3.0, means[index]))?;
    }
    for (area, mean) in areas.iter().zip(&mut means) {
        if *area <= 0.0 {
            return Err(bad(
                "enthalpy convection derivative has an empty physical trace",
            ));
        }
        *mean = finite(*mean * (1.0 / *area))?;
    }
    let mut replacements = Vec::new();
    replacements
        .try_reserve_exact(count)
        .map_err(|_| bad("enthalpy convection boundary allocation refused"))?;
    let mut wall_slopes = zeros(count)?;
    let mut ambient_slopes = zeros(count)?;
    for (index, law) in natural.iter().enumerate() {
        poll(cx, deadline)?;
        let (coefficient, partials) =
            buoyancy::differentiated_coefficient(law, means[index], pressure_pa)?;
        wall_slopes[index] = partials[0];
        ambient_slopes[index] = partials[1];
        replacements.push((regions[index], coefficient.htc_w_m2_k, law.ambient_k));
    }
    let mut references = Vec::new();
    references
        .try_reserve_exact(count)
        .map_err(|_| bad("enthalpy convection feedback matrix allocation refused"))?;
    for _ in 0..count {
        references.push(zeros(count)?);
    }
    let mut air_regions = BTreeSet::new();
    let air_heat_w = if paths.is_empty() {
        None
    } else {
        let law = affine_reference_law(cx, &paths, 64, count * count).map_err(air_error)?;
        let offset = natural.len();
        let air_count = count - offset;
        if law.regions() != &names[offset..] || law.wall_matrix().len() != air_count * air_count {
            return Err(bad(
                "enthalpy air-reference derivative lost its branch-major port order",
            ));
        }
        for i in 0..air_count {
            poll(cx, deadline)?;
            for j in 0..air_count {
                references[offset + i][offset + j] = law.wall_matrix()[i * air_count + j];
            }
        }
        let mut start = offset;
        let mut heat_w = 0.0;
        for path in &paths {
            poll(cx, deadline)?;
            let end = start + path.segments().len();
            // Use the physical march's arithmetic for the accepted residual,
            // not a separately rounded matrix-vector reconstruction.
            let march = path
                .march(&means[start..end])
                .map_err(|error| bad(format!("accepted enthalpy air march refused: {error}")))?;
            for (local, (segment, state)) in path.segments().iter().zip(&march.segments).enumerate()
            {
                let index = start + local;
                let ThermalBc::Robin {
                    htc: ScalarField::Uniform(h),
                    ..
                } = &problem.boundary.conditions()[regions[index]]
                else {
                    return Err(bad(
                        "enthalpy air derivative lost its original Robin coefficient",
                    ));
                };
                if *h != segment.htc_w_per_m2_k()
                    || (areas[index] - segment.area_m2()).abs()
                        > 128.0 * f64::EPSILON * areas[index].max(segment.area_m2())
                {
                    return Err(bad(
                        "enthalpy air and solid must share the same coefficient and wetted area",
                    ));
                }
                replacements.push((regions[index], *h, state.reference_temperature_k));
                air_regions.insert(names[index].clone());
            }
            heat_w = finite(heat_w + march.total_heat_rate_w)?;
            start = end;
        }
        Some(heat_w)
    };
    let boundary = problem
        .boundary
        .with_uniform_robin_replacements(&replacements)
        .map_err(lower)?;
    poll(cx, deadline)?;
    Ok(Some(Bound {
        boundary,
        names,
        references,
        wall_slopes,
        ambient_slopes,
        air_regions,
        air_heat_w,
    }))
}

impl Bound {
    /// Independently substitute actual air enthalpy gain for the original air
    /// convection rows. Natural reservoirs and radiation stay outside it.
    pub(super) fn check_air_energy(
        &self,
        physical: &EnthalpyEndpointEvaluation,
        dt_s: f64,
        tolerance_j: f64,
    ) -> Result<(), SolveRefusal> {
        let Some(air_heat_w) = self.air_heat_w else {
            return Ok(());
        };
        let mut off_path_w = 0.0;
        let mut matched = BTreeSet::new();
        for flux in &physical.convective_robin_fluxes {
            if self.air_regions.contains(&flux.region) {
                if !matched.insert(flux.region.as_str()) {
                    return Err(bad("enthalpy adjoint physical air audit repeats a port"));
                }
            } else {
                off_path_w = finite(off_path_w + flux.heat_rate_w)?;
            }
        }
        if matched.len() != self.air_regions.len() {
            return Err(bad("enthalpy adjoint physical air audit lost a port"));
        }
        let net_w = finite(
            physical.source_w
                + physical.dirichlet_in_w
                - physical.neumann_out_w
                - off_path_w
                - air_heat_w
                - physical.radiation_out_w,
        )?;
        let residual_j = finite(physical.stored_energy_change_j - dt_s * net_w)?;
        if residual_j.abs() > tolerance_j {
            return Err(bad(format!(
                "accepted enthalpy adjoint actual-air energy residual {residual_j:e} J exceeds {tolerance_j:e} J",
            )));
        }
        Ok(())
    }
}
