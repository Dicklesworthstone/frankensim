//! Request-owned finite-geometry sampling and explicit conservation projection.
//!
//! Binding is synchronous in the CLI. A successful record is retained for the
//! exact mesh, traces and emissivity declarations; failed or cancelled work
//! leaves the request retryable. No process-global sampling cache is used.

use std::sync::OnceLock;

use fs_conduction::radiation::view_factors::{
    DiffuseViewFactorConfig, DiffuseViewFactorEstimate, ProjectedViewFactors,
    ViewFactorProjection, estimate_diffuse_view_factors,
};

use super::*;

#[derive(Debug)]
pub(super) struct Generation {
    sampling: DiffuseViewFactorConfig,
    projection: ViewFactorProjection,
    record: OnceLock<Record>,
}

#[derive(Debug)]
struct Record {
    surfaces: Vec<RadiationSurface>,
    matrix: ViewFactorMatrix,
    evidence_json: String,
}

impl Generation {
    pub(super) fn parse(value: &J) -> Result<Self> {
        object(value, &[
            "seed", "samples_per_replicate", "replicates", "max_triangles",
            "max_ray_triangle_tests", "projection",
        ], "enclosure.generate_view_factors")?;
        let samples_per_replicate = u32::try_from(integer_raw(
            get(value, "samples_per_replicate")?, "samples_per_replicate",
        )?).map_err(|_| bad("samples_per_replicate exceeds u32"))?;
        if !samples_per_replicate.is_power_of_two() {
            return Err(bad("samples_per_replicate must be a positive power of two"));
        }
        let replicates = count(get(value, "replicates")?, "replicates", 32)? as u32;
        if replicates < 2 {
            return Err(bad("replicates must lie in 2..=32 for empirical standard errors"));
        }
        let sampling = DiffuseViewFactorConfig {
            seed: decimal_u64(get(value, "seed")?, "geometry seed")?,
            samples_per_replicate,
            replicates,
            max_triangles: count(get(value, "max_triangles")?, "max_triangles", 1_000_000)?,
            max_ray_triangle_tests: decimal_u64(
                get(value, "max_ray_triangle_tests")?, "max_ray_triangle_tests",
            )?,
        };
        if sampling.max_ray_triangle_tests == 0 {
            return Err(bad("max_ray_triangle_tests must be positive"));
        }
        let policy = object(get(value, "projection")?, &[
            "max_iterations", "row_sum_tolerance", "max_abs_adjustment",
        ], "enclosure.generate_view_factors.projection")?;
        let row_sum_tolerance = tolerance(get(policy, "row_sum_tolerance")?)?;
        if row_sum_tolerance == 0.0 {
            return Err(bad("generated row_sum_tolerance must be positive"));
        }
        let max_abs_adjustment = number(get(policy, "max_abs_adjustment")?, "max_abs_adjustment")?;
        if !(0.0..1.0).contains(&max_abs_adjustment) {
            return Err(bad("max_abs_adjustment must lie in [0,1)"));
        }
        Ok(Self {
            sampling,
            projection: ViewFactorProjection {
                max_iterations: count(get(policy, "max_iterations")?, "projection.max_iterations", 100_000)?,
                row_sum_tolerance,
                max_abs_adjustment,
            },
            record: OnceLock::new(),
        })
    }

    pub(super) fn bind(
        &self, mesh: &ConductionMesh, surfaces: &[RadiationSurface], cx: &Cx<'_>,
    ) -> Result<&ViewFactorMatrix> {
        poll(cx)?;
        if self.record.get().is_none() {
            self.admit_work(mesh, surfaces.len())?;
            let estimate = estimate_diffuse_view_factors(cx, mesh, surfaces, self.sampling);
            poll(cx)?;
            let estimate = estimate.map_err(producer)?;
            let projected = estimate.project_closed(cx, self.projection);
            poll(cx)?;
            let projected = projected.map_err(producer)?;
            let evidence_json = self.render(&estimate, &projected)?;
            let record = Record {
                surfaces: surfaces.to_vec(),
                matrix: projected.into_matrix(),
                evidence_json,
            };
            poll(cx)?;
            // Retain only complete successful work. Read back and check the
            // retained record even if a future caller races this initialization.
            let _ = self.record.set(record);
        }
        let record = self.record.get().ok_or_else(|| producer("generated view factors were not retained"))?;
        // RadiationSurface equality includes the complete physical mesh hash,
        // selected face slots, areas, names and card-backed emissivity.
        if record.surfaces.as_slice() != surfaces {
            return Err(bad("generated view factors are bound to a different mesh, patch partition or emissivity; parse a new request"));
        }
        poll(cx)?;
        Ok(&record.matrix)
    }

    pub(super) fn factor_data(&self) -> Result<(&ViewFactorMatrix, &str)> {
        let record = self.record.get().ok_or_else(|| bad("generated enclosure has not been bound"))?;
        Ok((&record.matrix, record.evidence_json.as_str()))
    }

    fn admit_work(&self, mesh: &ConductionMesh, patches: usize) -> Result<()> {
        let triangles = mesh.boundary().len();
        if triangles > self.sampling.max_triangles {
            return Err(work_budget(format!(
                "{triangles} opaque boundary triangles exceed max_triangles={}",
                self.sampling.max_triangles,
            )));
        }
        let work = (patches as u64)
            .checked_mul(u64::from(self.sampling.replicates))
            .and_then(|n| n.checked_mul(u64::from(self.sampling.samples_per_replicate)))
            .and_then(|n| n.checked_mul(triangles as u64))
            .ok_or_else(|| work_budget("geometry ray/triangle work count overflow"))?;
        if work > self.sampling.max_ray_triangle_tests {
            return Err(work_budget(format!(
                "geometry sampling requires at most {work} ray/triangle tests; declared budget is {}",
                self.sampling.max_ray_triangle_tests,
            )));
        }
        Ok(())
    }

    fn render(
        &self, estimate: &DiffuseViewFactorEstimate, projected: &ProjectedViewFactors,
    ) -> Result<String> {
        let names = estimate.names().iter().map(|name| quote(name)).collect::<Vec<_>>().join(",");
        let factors = estimate.factors().iter().map(|row| numbers(row)).collect::<Result<Vec<_>>>()?.join(",");
        let errors = estimate.standard_errors().iter().map(|row| numbers(row)).collect::<Result<Vec<_>>>()?.join(",");
        let hits = estimate.hit_counts().iter().map(|row| counts(row)).collect::<Vec<_>>().join(",");
        Ok(format!(
            "{{\"kind\":\"geometry-qmc\",\"method\":\"owen-scrambled-sobol-cosine-hemisphere\",\"seed\":{},\"samples_per_replicate\":{},\"replicates\":{},\"max_triangles\":{},\"max_ray_triangle_tests\":{},\"samples_per_surface\":{},\"ray_triangle_tests\":{},\"mesh_identity\":{},\"estimate_identity\":{},\"surface_names\":[{names}],\"areas_m2\":{},\"raw_view_factors\":[{factors}],\"empirical_standard_errors\":[{errors}],\"hit_counts\":[{hits}],\"escaped_counts\":{},\"unassigned_counts\":{},\"backface_counts\":{},\"count_encoding\":\"decimal-u64-string\",\"projection\":{{\"max_iterations\":{},\"iterations\":{},\"row_sum_tolerance\":{},\"max_abs_adjustment_budget\":{},\"max_abs_adjustment\":{}}},\"scope\":\"first-hit visibility on the exact declared tetrahedral boundary, with every boundary triangle opaque; independent scramble standard errors are empirical, not calibrated confidence bounds or geometric error certificates; zero sampled losses do not prove a closed enclosure; projection is an explicitly budgeted numerical conservation correction, not sampling-error calibration; no physical validation or geometry derivative\"}}",
            quote(&self.sampling.seed.to_string()), self.sampling.samples_per_replicate,
            self.sampling.replicates, self.sampling.max_triangles,
            quote(&self.sampling.max_ray_triangle_tests.to_string()),
            quote(&estimate.samples_per_surface().to_string()),
            quote(&estimate.ray_triangle_tests().to_string()),
            quote(&estimate.mesh_identity().to_hex()), quote(&estimate.identity().to_hex()),
            numbers(estimate.areas_m2())?, counts(estimate.escaped_counts()),
            counts(estimate.unassigned_counts()), counts(estimate.backface_counts()),
            self.projection.max_iterations, projected.iterations(),
            num(self.projection.row_sum_tolerance)?, num(self.projection.max_abs_adjustment)?,
            num(projected.max_abs_adjustment())?,
        ))
    }
}

fn decimal_u64(value: &J, field: &str) -> Result<u64> {
    let text = string(value, field)?;
    if !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(bad(format!("{field} must be a decimal u64 string")));
    }
    text.parse().map_err(|_| bad(format!("{field} must be a decimal u64 string")))
}

fn counts(values: &[u64]) -> String {
    format!("[{}]", values.iter().map(|value| quote(&value.to_string())).collect::<Vec<_>>().join(","))
}

fn work_budget(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-radiation-budget", message: message.into() }
}
