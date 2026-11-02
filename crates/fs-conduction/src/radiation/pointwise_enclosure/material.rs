//! Nodal source emissivity in the original diffuse-gray patch equations.
//!
//! Irradiation is uniform per admitted view-factor patch. With P1 epsilon and
//! T, J=<epsilon B>+(1-<epsilon>)G and G=FJ. Numerical area averages are not
//! source claims. Every state rebuilds the SAME original radiosity matrix.

use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;

use super::*;
use crate::radiation::pointwise::material::{
    TRIANGLE_DEGREE_SIX,
    source::{DerivativeMode, MaterialEmissivitySource},
};

pub mod adjoint;

/// Original caloric/material source for one complete named enclosure trace.
#[derive(Debug)]
pub struct MaterialEnclosureRadiationPatch<'a> {
    pub region: &'a str,
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub material: ThermalMaterialPoint<'a>,
    pub selection: MaterialPropertySelection,
}

/// Source-dependent gray reflection on an existing closed or open trace owner.
/// Original geometry, view factors, bath and watt gates remain authoritative.
/// Spectral data require their own explicit path. Source adjoints retain the
/// original ordinates and the current, source-dependent reflection matrix.
#[derive(Debug)]
pub struct MaterialPointwiseEnclosureRadiation<'a> {
    base: &'a PointwiseEnclosureRadiation<'a>,
    sources: Vec<MaterialEmissivitySource<'a>>,
    curves: Vec<Option<&'a EquilibriumEnthalpyPhaseCurve>>,
    view_factors: Arc<ViewFactorMatrix>,
}

#[derive(Debug, Clone)]
struct MaterialFaceTangent {
    vertices: [usize; 3],
    matrix: [[f64; 3]; 3],
}

/// Complete h tangent: local emission/material response minus the current
/// reflected irradiation response. No dense nodal matrix or frozen reflection.
#[derive(Debug, Clone)]
pub struct MaterialPointwiseEnclosureAssembly {
    pub outward_load_w: Vec<f64>,
    /// Positive frozen-emissivity emission diagonal in W/K, for preconditioning.
    /// The complete h tangent can be signed and nonsymmetric.
    pub emission_temperature_diagonal_w_k: Vec<f64>,
    pub report: PointwiseEnclosureReport,
    geometry: Arc<Geometry>,
    view_factors: Arc<ViewFactorMatrix>,
    factors: DenseRadiosityFactorization,
    local: Vec<MaterialFaceTangent>,
    /// Integrals of N_i epsilon_h, m², in each original patch's vertex order.
    absorption_weights: Vec<Vec<f64>>,
    /// <epsilon dB/dh + (B-G) depsilon/dh>, in patch vertex order.
    patch_slopes: Vec<Vec<f64>>,
    // Retain centered irradiation so source controls never subtract two large
    // absolute blackbody radiances near thermal equilibrium.
    reference_temperature_k: f64,
    irradiation_corrections_w_m2: Vec<f64>,
    material_parameters: Option<adjoint::MaterialEnclosureParameterData>,
}

impl MaterialPointwiseEnclosureAssembly {
    pub fn enthalpy_tangent_apply(&self, cx: &Cx<'_>, direction: &[f64], output: &mut [f64])
        -> Result<(), ConductionError>
    {
        self.check_fields(direction, output)?; output.fill(0.0);
        self.add_scaled_tangent(cx, direction, 1.0, |_| true, output)
    }

    pub fn enthalpy_tangent_apply_transpose(&self, cx: &Cx<'_>, direction: &[f64], output: &mut [f64])
        -> Result<(), ConductionError>
    {
        self.check_fields(direction, output)?; output.fill(0.0);
        self.add_scaled_tangent_transpose(cx, direction, 1.0, |_| true, output)
    }

    pub(crate) fn add_scaled_tangent(&self, cx: &Cx<'_>, direction: &[f64], scale: f64,
        include_row: impl Fn(usize) -> bool, output: &mut [f64]) -> Result<(), ConductionError>
    {
        poll(cx, 0)?; self.check_fields(direction, output)?;
        for (slot, face) in self.local.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            for (i, &vertex) in face.vertices.iter().enumerate() {
                if include_row(vertex) {
                    let applied = dot(&face.matrix[i], &face.vertices.map(|v| direction[v]))?;
                    output[vertex] = finite(scale.mul_add(applied, output[vertex]))?;
                }
            }
        }
        let mut rhs = vec![0.0; self.view_factors.len()];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            for (slot, (&vertex, &slope)) in geometry.vertices.iter().zip(&self.patch_slopes[patch]).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                rhs[patch] = finite(slope.mul_add(direction[vertex], rhs[patch]))?;
            }
        }
        // K dJ = d<epsilon B> - G d<epsilon>. The last open-bath RHS is zero.
        let radiosity = self.factors.solve(cx, &rhs)?;
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let irradiation = dot(&self.view_factors.factors()[patch], &radiosity)?;
            for (slot, (&vertex, &weight)) in geometry.vertices.iter().zip(&self.absorption_weights[patch]).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                if include_row(vertex) {
                    output[vertex] = finite((-scale * weight).mul_add(irradiation, output[vertex]))?;
                }
            }
        }
        poll(cx, self.local.len())
    }

    pub(crate) fn add_scaled_tangent_transpose(&self, cx: &Cx<'_>, direction: &[f64], scale: f64,
        include_row: impl Fn(usize) -> bool, output: &mut [f64]) -> Result<(), ConductionError>
    {
        poll(cx, 0)?; self.check_fields(direction, output)?;
        let value = |vertex| if include_row(vertex) { direction[vertex] } else { 0.0 };
        for (slot, face) in self.local.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            for (j, &vertex) in face.vertices.iter().enumerate() {
                let mut applied = 0.0;
                for (i, &row) in face.vertices.iter().enumerate() {
                    applied = finite(face.matrix[i][j].mul_add(value(row), applied))?;
                }
                output[vertex] = finite(scale.mul_add(applied, output[vertex]))?;
            }
        }
        let mut radiosity_bar = vec![0.0; self.view_factors.len()];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let mut irradiation_bar = 0.0;
            for (slot, (&vertex, &weight)) in geometry.vertices.iter().zip(&self.absorption_weights[patch]).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                irradiation_bar = finite(weight.mul_add(value(vertex), irradiation_bar))?;
            }
            for (column, bar) in radiosity_bar.iter_mut().enumerate() {
                *bar = finite(self.view_factors.factors()[patch][column].mul_add(irradiation_bar, *bar))?;
            }
        }
        let source_bar = self.factors.solve_transpose(cx, &radiosity_bar)?;
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            for (slot, (&vertex, &slope)) in geometry.vertices.iter().zip(&self.patch_slopes[patch]).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                output[vertex] = finite((-scale * slope).mul_add(source_bar[patch], output[vertex]))?;
            }
        }
        poll(cx, self.local.len())
    }

    pub(crate) fn admit_adjoint(&self, cx: &Cx<'_>) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        if self.material_parameters.is_none() {
            return Err(invalid("material enclosure adjoint requires explicit original-source admission"));
        }
        Ok(())
    }

    pub(crate) fn parameter_pullback(&self, cx: &Cx<'_>, _weights: &[f64])
        -> Result<PointwiseRadiationParameterGradient, ConductionError>
    {
        poll(cx, 0)?;
        Err(invalid("material enclosure controls require the original-source parameter pullback"))
    }

    fn check_fields(&self, direction: &[f64], output: &[f64]) -> Result<(), ConductionError> {
        if direction.len() != self.geometry.vertices || output.len() != self.geometry.vertices {
            return Err(invalid("material enclosure h action requires complete nodal fields"));
        }
        if self.patch_slopes.len() != self.geometry.patches.len() {
            return Err(invalid("material enclosure load-only evaluation has no tangent"));
        }
        Ok(())
    }
}

impl<'a> MaterialPointwiseEnclosureRadiation<'a> {
    pub fn new(cx: &Cx<'_>, base: &'a PointwiseEnclosureRadiation<'a>,
        patches: &[MaterialEnclosureRadiationPatch<'a>]) -> Result<Self, ConductionError>
    {
        poll(cx, 0)?;
        if base.system.band().is_some() || !base.additional_bands.is_empty()
            || patches.len() != base.geometry.patches.len() {
            return Err(invalid("material enclosure requires one original gray source per existing patch"));
        }
        let mut sources = Vec::with_capacity(patches.len());
        let mut curves: Vec<Option<&'a EquilibriumEnthalpyPhaseCurve>> = vec![None; base.mesh.vertex_count()];
        for (index, patch) in patches.iter().enumerate() {
            poll(cx, index)?;
            if patch.region != base.system.surfaces()[index].name() {
                return Err(invalid("material enclosure bindings must retain the original named surface order"));
            }
            let source = MaterialEmissivitySource::prepare(patch.region, patch.curve,
                patch.material, patch.selection.clone())?;
            for (slot, &vertex) in base.geometry.patches[index].vertices.iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                if curves[vertex].is_some_and(|curve| !std::ptr::eq(curve, patch.curve)) {
                    return Err(invalid("shared material enclosure nodes require one original caloric chart"));
                }
                curves[vertex] = Some(patch.curve);
            }
            sources.push(source);
        }
        Ok(Self { base, sources, curves, view_factors: Arc::new(base.system.view_factors().clone()) })
    }

    pub(crate) fn matches(&self, mesh: &ConductionMesh, boundary: &ThermalBoundary) -> bool {
        self.base.matches(mesh, boundary)
    }

    /// Closed internal exchange has no external black-bath heat account.
    pub const fn is_closed(&self) -> bool { self.base.is_closed() }

    pub(crate) fn matches_curves<'c>(&self, cx: &Cx<'_>,
        mut curve_at: impl FnMut(usize) -> &'c EquilibriumEnthalpyPhaseCurve) -> Result<bool, ConductionError>
    {
        for (vertex, curve) in self.curves.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            if curve.is_some_and(|curve| !std::ptr::eq(curve, curve_at(vertex))) { return Ok(false); }
        }
        Ok(true)
    }

    pub fn assemble(&self, cx: &Cx<'_>, specific_h: &[f64])
        -> Result<MaterialPointwiseEnclosureAssembly, ConductionError>
    { self.integrate::<true>(cx, specific_h, None) }

    pub(crate) fn assemble_at_temperature(&self, cx: &Cx<'_>, specific_h: &[f64], temperature: &[f64])
        -> Result<MaterialPointwiseEnclosureAssembly, ConductionError>
    { self.integrate::<true>(cx, specific_h, Some(temperature)) }

    pub(crate) fn load_and_report(&self, cx: &Cx<'_>, specific_h: &[f64], temperature: &[f64])
        -> Result<(Vec<f64>, PointwiseEnclosureReport), ConductionError>
    {
        let assembly = self.integrate::<false>(cx, specific_h, Some(temperature))?;
        Ok((assembly.outward_load_w, assembly.report))
    }

    fn integrate<const TANGENT: bool>(&self, cx: &Cx<'_>, h: &[f64], supplied_t: Option<&[f64]>)
        -> Result<MaterialPointwiseEnclosureAssembly, ConductionError>
    {
        poll(cx, 0)?;
        let n = self.base.geometry.vertices; let p = self.sources.len();
        if h.len() != n || supplied_t.is_some_and(|t| t.len() != n) {
            return Err(invalid("material enclosure requires complete original h and temperature fields"));
        }
        let mut samples = Vec::with_capacity(p);
        let mut means = vec![0.0; p]; let mut minimum = vec![f64::INFINITY; p];
        let mut maximum = vec![f64::NEG_INFINITY; p];
        let mode = if TANGENT { DerivativeMode::Primal } else { DerivativeMode::None };
        for (patch, geometry) in self.base.geometry.patches.iter().enumerate() {
            let mut nodes = Vec::with_capacity(geometry.vertices.len());
            for (slot, &vertex) in geometry.vertices.iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                let node = self.sources[patch].resolve(vertex, h[vertex], supplied_t.map(|t| t[vertex]), mode)?;
                let t = node.state.phase.temperature_k();
                minimum[patch] = minimum[patch].min(t); maximum[patch] = maximum[patch].max(t);
                nodes.push(node);
            }
            // Center at the largest admitted value. Uniform black traces stay
            // exactly one without rounding the physical average above one.
            let anchor = nodes.iter().map(|node| node.state.emissivity.value()).fold(0.0_f64, f64::max);
            let mut delta = 0.0;
            for (slot, (node, &weight)) in nodes.iter().zip(geometry.weights.iter()).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                delta = finite(weight.mul_add(node.state.emissivity.value() - anchor, delta))?;
            }
            means[patch] = finite(anchor + delta / self.base.system.surfaces()[patch].area_m2())?;
            if !(means[patch] > 0.0 && means[patch] <= 1.0) {
                return Err(invalid("integrated material emissivity is outside its physical range"));
            }
            samples.push(nodes);
        }
        let surface_reference = minimum.iter().copied().fold(f64::INFINITY, f64::min);
        let reference = self.base.ambient_temperature_k.map_or(surface_reference, |t| t.min(surface_reference));
        let background = finite(SIGMA * reference * reference * reference * reference)?;
        let mut emission = vec![0.0; p]; let mut blackbody = vec![0.0; p]; let mut load = vec![0.0; n];
        let mut absorption: Vec<_> = self.base.geometry.patches.iter().map(|g| vec![0.0; g.vertices.len()]).collect();
        let mut diagonal = if TANGENT { vec![0.0; n] } else { Vec::new() };
        let mut local = Vec::with_capacity(if TANGENT { self.base.geometry.faces.len() } else { 0 });
        for (slot, face) in self.base.geometry.faces.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            let nodes = face.local.map(|i| &samples[face.patch][i]);
            let t = nodes.map(|node| node.state.phase.temperature_k());
            let epsilon = nodes.map(|node| node.state.emissivity.value());
            let epsilon_sum: f64 = epsilon.iter().sum();
            for i in 0..3 {
                let value = &mut absorption[face.patch][face.local[i]];
                *value = finite(*value + face.area / 12.0 * (epsilon_sum + epsilon[i]))?;
            }
            let mut matrix = [[0.0; 3]; 3];
            for (bary, weight) in TRIANGLE_DEGREE_SIX {
                let tq = finite(t[0] + (t[1] - t[0]).mul_add(bary[1], (t[2] - t[0]) * bary[2]))?;
                let eq = finite(epsilon[0] + (epsilon[1] - epsilon[0]).mul_add(bary[1], (epsilon[2] - epsilon[0]) * bary[2]))?;
                let area_weight = finite(face.area * weight)?;
                let delta = finite(balance::emission_difference(SIGMA, tq, reference))?;
                let power = finite(area_weight * eq * delta)?;
                emission[face.patch] = finite(emission[face.patch] + power)?;
                blackbody[face.patch] = finite(blackbody[face.patch] + area_weight * delta)?;
                let emitted_slope = finite(area_weight * eq * 4.0 * SIGMA * tq * tq * tq)?;
                for i in 0..3 {
                    load[face.vertices[i]] = finite(bary[i].mul_add(power, load[face.vertices[i]]))?;
                    if TANGENT {
                        diagonal[face.vertices[i]] = finite(diagonal[face.vertices[i]] + bary[i] * bary[i] * emitted_slope)?;
                        for j in 0..3 {
                            matrix[i][j] = finite(matrix[i][j] + bary[i] * bary[j]
                                * (emitted_slope * nodes[j].dt_dh + area_weight * delta * nodes[j].depsilon_dh))?;
                        }
                    }
                }
            }
            if TANGENT { local.push(MaterialFaceTangent { vertices: face.vertices, matrix }); }
        }
        let views = self.view_factors.factors(); let equations = views.len();
        let mut matrix = vec![vec![0.0; equations]; equations]; let mut rhs = vec![0.0; equations];
        let mut row_sums = vec![0.0; equations];
        for row in 0..equations {
            poll(cx, row)?;
            let epsilon = if row < p { means[row] } else { 1.0 };
            for column in 0..equations {
                matrix[row][column] = if row == column { 1.0 } else { 0.0 } - (1.0 - epsilon) * views[row][column];
            }
            row_sums[row] = finite(views[row].iter().sum())?;
            if row < p {
                let area = self.base.system.surfaces()[row].area_m2();
                emission[row] = finite(emission[row] / area)?; blackbody[row] = finite(blackbody[row] / area)?;
                rhs[row] = finite(emission[row] + (1.0 - epsilon) * background * (row_sums[row] - 1.0))?;
            }
        }
        if let Some(bath) = self.base.ambient_temperature_k { rhs[p] = finite(balance::emission_difference(SIGMA, bath, reference))?; }
        let factors = DenseRadiosityFactorization::new(cx, &matrix)?;
        let radiosity = factors.solve(cx, &rhs)?;
        let mut corrections = vec![0.0; p]; let mut reports = Vec::with_capacity(p);
        let mut total = 0.0; let mut absolute = 0.0; let mut residual_max = 0.0_f64; let mut max_area = 0.0_f64;
        for (patch, nodes) in samples.iter().enumerate() {
            poll(cx, patch)?;
            let surface = &self.base.system.surfaces()[patch]; let geometry = &self.base.geometry.patches[patch];
            let irradiation = dot(&views[patch], &radiosity)?;
            corrections[patch] = finite(background.mul_add(row_sums[patch] - 1.0, irradiation))?;
            residual_max = residual_max.max(finite(radiosity[patch] - (1.0 - means[patch]) * irradiation - rhs[patch])?.abs());
            max_area = max_area.max(surface.area_m2());
            let mut mean_delta = 0.0;
            for (slot, &vertex) in geometry.vertices.iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                load[vertex] = finite(load[vertex] - absorption[patch][slot] * corrections[patch])?;
                mean_delta = finite(geometry.weights[slot].mul_add(nodes[slot].state.phase.temperature_k() - reference, mean_delta))?;
            }
            let heat = finite(surface.area_m2() * (emission[patch] - means[patch] * corrections[patch]))?;
            total = finite(total + heat)?; absolute = finite(absolute + heat.abs())?;
            let j = finite(background + radiosity[patch])?;
            let g = finite(background.mul_add(row_sums[patch], irradiation))?;
            if j < 0.0 || g < 0.0 { return Err(invalid("material enclosure produced negative radiosity or irradiation")); }
            reports.push(PointwiseEnclosurePatchReport {
                name: surface.name().to_owned(),
                emissivity: PointwisePatchEmissivity::Gray(self.sources[patch].snapshot().clone()),
                area_m2: surface.area_m2(), mean_temperature_k: finite(reference + mean_delta / surface.area_m2())?,
                minimum_temperature_k: minimum[patch], maximum_temperature_k: maximum[patch],
                mean_blackbody_emission_w_m2: finite(background + blackbody[patch])?,
                radiosity_w_m2: j, irradiation_w_m2: g, outward_heat_w: heat,
                material_states: Vec::new(), material_mean_emissivity: Some(means[patch]),
            });
        }
        let ambient = if let Some(temperature_k) = self.base.ambient_temperature_k {
            let irradiation = dot(&views[p], &radiosity)?; let area = self.view_factors.areas_m2()[p];
            let heat = finite(area * balance::outward_flux(1.0, rhs[p], irradiation, background, row_sums[p]))?;
            residual_max = residual_max.max(finite(radiosity[p] - rhs[p])?.abs()); max_area = max_area.max(area);
            let j = finite(background + radiosity[p])?; let g = finite(background.mul_add(row_sums[p], irradiation))?;
            if j < 0.0 || g < 0.0 { return Err(invalid("material enclosure black bath produced negative radiosity or irradiation")); }
            Some(PointwiseBlackBathReport { temperature_k, effective_area_m2: area,
                radiosity_w_m2: j, irradiation_w_m2: g, outward_heat_w: heat })
        } else { None };
        let bath_heat = ambient.as_ref().map_or(0.0, |bath| bath.outward_heat_w);
        let closure = finite(total + bath_heat)?;
        let tolerance = finite(self.base.tolerance.relative_heat.mul_add(absolute + bath_heat.abs(), self.base.tolerance.absolute_heat_w))?;
        if closure.abs() > tolerance || finite(residual_max * max_area)? > tolerance {
            return Err(invalid("material enclosure original radiosity equations or total heat exceed the watt gate"));
        }
        let mut slopes: Vec<_> = if TANGENT {
            self.base.geometry.patches.iter().map(|g| vec![0.0; g.vertices.len()]).collect()
        } else { Vec::new() };
        if TANGENT {
            for (slot, (face, tangent)) in self.base.geometry.faces.iter().zip(&mut local).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                let area = self.base.system.surfaces()[face.patch].area_m2();
                for j in 0..3 {
                    let eps_h = samples[face.patch][face.local[j]].depsilon_dh;
                    for i in 0..3 {
                        // Exact degree-two mass moment of Ni Nj times nodal epsilon'_j.
                        tangent.matrix[i][j] = finite(tangent.matrix[i][j] - corrections[face.patch]
                            * face.area / 12.0 * (if i == j { 2.0 } else { 1.0 }) * eps_h)?;
                        slopes[face.patch][face.local[j]] = finite(slopes[face.patch][face.local[j]] + tangent.matrix[i][j] / area)?;
                    }
                }
            }
        }
        for (report, nodes) in reports.iter_mut().zip(samples) {
            report.material_states = nodes.into_iter().map(|node| node.state).collect();
        }
        poll(cx, self.base.geometry.faces.len())?;
        Ok(MaterialPointwiseEnclosureAssembly { outward_load_w: load,
            emission_temperature_diagonal_w_k: diagonal,
            report: PointwiseEnclosureReport { patches: reports, outward_heat_w: total,
                absolute_patch_heat_w: absolute, linear_residual_max_w_m2: residual_max,
                view_factor_identity: self.view_factors.identity(), ambient,
                system_energy_closure_w: closure, spectral: None },
            geometry: Arc::clone(&self.base.geometry), view_factors: Arc::clone(&self.view_factors),
            factors, local, absorption_weights: absorption, patch_slopes: slopes,
            reference_temperature_k: reference, irradiation_corrections_w_m2: corrections,
            material_parameters: None })
    }
}
