//! Event-aligned native heat inputs on the retained volume mesh.
//!
//! Scheduled watts are absolute delivered regional power. Each interval uses
//! the load at its left endpoint, so a switch supplies the following interval
//! and never retroactively changes the energy of the preceding one.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_conduction::{ConductionMesh, PowerMap, ScalarField};
use fs_exec::Cx;
use fs_project::{ConductionTransient, ProjectSpec, TransientRegionPower};

use super::{SolveRefusal, bad, json_string, lower, number, poll};

pub(super) struct TimeGrids {
    pub(super) coarse: Vec<f64>,
    pub(super) fine: Vec<f64>,
}

impl TimeGrids {
    pub(super) fn new(policy: &ConductionTransient) -> Result<Self, SolveRefusal> {
        let coarse = policy
            .coarse_step_ends_s()
            .map_err(|error| bad(error.what))?;
        let mut fine = Vec::with_capacity(coarse.len() * 2);
        if policy.power_schedules.is_empty() {
            // Preserve the original uniform-grid arithmetic for static loads.
            let count = coarse.len() * 2;
            for ordinal in 1..=count {
                fine.push(if ordinal == count {
                    policy.horizon.value
                } else {
                    policy.horizon.value * (ordinal as f64 / count as f64)
                });
            }
        } else {
            let mut start = 0.0;
            for &end in &coarse {
                let middle = start + 0.5 * (end - start);
                if !(middle > start && middle < end) {
                    return Err(bad(
                        "a workload interval cannot represent both positive half steps",
                    ));
                }
                fine.extend([middle, end]);
                start = end;
            }
        }
        let mut start = 0.0;
        for &end in &fine {
            if !(end.is_finite() && end > start) {
                return Err(bad("the fine time grid contains an unrepresentable step"));
            }
            start = end;
        }
        Ok(Self { coarse, fine })
    }

    pub(super) fn largest_step(ends: &[f64]) -> f64 {
        let mut start = 0.0;
        let mut largest = 0.0_f64;
        for &end in ends {
            largest = largest.max(end - start);
            start = end;
        }
        largest
    }
}

/// A bound regional power history. Spatial projection uses the same P1
/// source carrier and conservative regional projector as static power.
pub(super) struct Workload<'a> {
    static_watts: BTreeMap<u32, f64>,
    volumes: BTreeMap<u32, f64>,
    schedules: Vec<(u32, &'a TransientRegionPower)>,
    pub(super) receipt: String,
}

impl<'a> Workload<'a> {
    pub(super) fn bind(
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        policy: &'a ConductionTransient,
        mesh: &ConductionMesh,
        labels: &[u32],
        region_ids: &BTreeMap<String, u32>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Option<Self>, SolveRefusal> {
        if policy.power_schedules.is_empty() {
            return Ok(None);
        }
        if labels.len() != mesh.element_count() {
            return Err(bad(
                "workload region labels do not match the retained volume mesh",
            ));
        }
        let mut volumes = BTreeMap::<u32, f64>::new();
        for (element, &id) in labels.iter().enumerate() {
            if element % 256 == 0 {
                poll(cx, deadline)?;
            }
            *volumes.entry(id).or_default() += mesh.element_volume(element);
        }
        if volumes
            .values()
            .any(|volume| !volume.is_finite() || *volume <= 0.0)
        {
            return Err(bad(
                "workload projection needs positive finite regional mesh volumes",
            ));
        }
        let surfaces = super::super::surface_heat(spec)?;
        let mut static_watts: BTreeMap<_, _> = region_ids.values().map(|&id| (id, 0.0)).collect();
        for row in spec.power.iter().flatten() {
            if surfaces.names.contains(&row.region) {
                // Those watts already enter the immutable boundary partition.
                continue;
            }
            let id = region_ids.get(&row.region).ok_or_else(|| {
                bad(format!(
                    "static power names unknown workload region `{}`",
                    row.region
                ))
            })?;
            let value = row.watts.value * row.duty;
            if !(value.is_finite() && value >= 0.0) {
                return Err(bad(
                    "static regional delivered power is not finite and nonnegative",
                ));
            }
            let total = static_watts
                .get_mut(id)
                .ok_or_else(|| bad("unbound static power region"))?;
            *total += value;
            if !total.is_finite() {
                return Err(bad("summed static regional power overflows"));
            }
        }
        let mut schedules = Vec::with_capacity(policy.power_schedules.len());
        let mut rows = Vec::with_capacity(policy.power_schedules.len());
        let mut integrated = 0.0;
        for row in &policy.power_schedules {
            poll(cx, deadline)?;
            let id = *region_ids.get(&row.region).ok_or_else(|| {
                bad(format!(
                    "scheduled power must name a seeded volume region: `{}`",
                    row.region
                ))
            })?;
            if surfaces.names.contains(&row.region) || !volumes.contains_key(&id) {
                return Err(bad(
                    "scheduled power requires a retained volume region, not a surface",
                ));
            }
            let mut previous = 0.0;
            let mut steps = Vec::with_capacity(row.steps.len());
            for step in &row.steps {
                integrated += (step.until.value - previous) * step.watts.value;
                previous = step.until.value;
                steps.push(format!(
                    "{{\"until_s\":{},\"watts\":{}}}",
                    number(step.until.value)?,
                    number(step.watts.value)?
                ));
            }
            rows.push(format!(
                "{{\"region\":{},\"source\":{},\"steps\":[{}]}}",
                json_string(&row.region),
                json_string(&row.source),
                steps.join(",")
            ));
            schedules.push((id, row));
        }
        let receipt = format!(
            "{{\"mode\":\"piecewise-constant-regional-power\",\"watts_basis\":\"absolute-delivered-replaces-static-region-power\",\"interval_policy\":\"left-endpoint-until-next-switch\",\"integrated_scheduled_input_j\":{},\"schedules\":[{}],\"scope\":\"Only named volume-region power is replaced. Other regional power and declared surface heat remain static; no schedule is inferred from a duty factor.\"}}",
            number(integrated)?,
            rows.join(",")
        );
        Ok(Some(Self {
            static_watts,
            volumes,
            schedules,
            receipt,
        }))
    }

    pub(super) fn source(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        labels: &[u32],
        start_s: f64,
        deadline: Option<(Instant, f64)>,
    ) -> Result<ScalarField, SolveRefusal> {
        poll(cx, deadline)?;
        let mut watts = self.static_watts.clone();
        for (id, schedule) in &self.schedules {
            let index = schedule
                .steps
                .partition_point(|step| step.until.value <= start_s);
            let step = schedule.steps.get(index).ok_or_else(|| {
                bad("a workload time step extends beyond its declared regional power history")
            })?;
            watts.insert(*id, step.watts.value);
        }
        let mut densities = BTreeMap::new();
        for (id, power) in watts {
            let volume = self
                .volumes
                .get(&id)
                .ok_or_else(|| bad("a declared power region has no retained mesh volume"))?;
            let density = power / volume;
            if !(density.is_finite() && density >= 0.0) {
                return Err(bad(
                    "workload volumetric source is not finite and nonnegative",
                ));
            }
            densities.insert(id, density);
        }
        let source =
            PowerMap::regional_volumetric_source(mesh, labels, &densities).map_err(lower)?;
        poll(cx, deadline)?;
        Ok(source)
    }
}
