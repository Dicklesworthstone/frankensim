//! A geometry/material-derived moving slit, reduced from the existing DKT plate.
//!
//! For retained shape phi, s is the length-weighted displacement of the declared
//! slit edges. In opening coordinates x = s q, M = phi'Mphi / s² and the uniform
//! pressure area is A = integral(phi dA) / s. Both pressure force and swept flow
//! use this SAME A: F = -A deltaP, U_swept = -A x_dot. No face length, mass,
//! stiffness, frequency, or closing pressure is prescribed independently.
//!
//! This is one linear plate mode with the existing P1 displacement-trace
//! quadrature, a uniform rest gap and one generalized lay contact. It is not a
//! distributed reed/lay closure, fluid-loaded modal reduction or exact instrument
//! reconstruction. The explicit slit-uniformity and linear-slope allowances are
//! approximation limits, not experimental validation. Damping remains an
//! independent physical input; an elastic card does not imply a loss law.

/// Spatial lay clearances and partially closed slit geometry.
pub mod closure;

use crate::thin_plate::ResolvedPlateChart;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::{BernoulliAperture, dynamic::DynamicApertureSpec};
use fs_exec::CancelGate;
use fs_plate::{AssemblyOptions, PlateChart, SliceOptions, modes};

/// Physical surface motion relative to the supplied rest chart, at one node.
/// No unit-amplitude animation or rescaling is applied after the dynamics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateApertureMotion {
    /// [transverse displacement in m, two DKT rotations in rad].
    pub displacement_rotation: [f64; 3],
    /// [transverse velocity in m/s, two angular velocities in rad/s].
    pub velocity_rotation_rate: [f64; 3],
}

/// Geometry, approximation limits and bounded offline mode selection.
#[derive(Debug, Clone)]
pub struct PlateApertureOptions {
    /// Supports act on the chart's explicit support-node set; no rim is inferred.
    pub assembly: AssemblyOptions,
    /// Search interval in squared rad/s; this is not a prescribed resonance.
    pub eigenvalue_window: (f64, f64),
    /// Ascending zero-based index WITHIN the searched window.
    pub mode_index: usize,
    /// Existing sparse modal owner's iteration, residual and factor controls.
    pub eigensolver: SliceOptions,
    /// Input-node allowance, checked before assembly.
    pub max_nodes: usize,
    /// Input-triangle allowance, checked before assembly.
    pub max_triangles: usize,
    /// Unique boundary edges forming the moving slit; their total length is width.
    pub slit_edges: Vec<[usize; 2]>,
    /// Positive mean rest opening [m].
    pub rest_opening_m: f64,
    /// Explicit nonnegative viscous damping ratio; never inferred from elasticity.
    pub damping_ratio: f64,
    /// Maximum |normalized slit-node motion - 1|. One DOF approximates the slit
    /// by its mean gap; this declares how nonuniform that retained shape may be.
    pub max_slit_mode_variation: f64,
    /// Positive bound on the largest nodal rotation/P1 displacement slope.
    pub max_slope: f64,
}

#[derive(Debug, Clone)]
enum PlateSource {
    Authored(PlateChart),
    Material(ResolvedPlateChart),
}
impl PlateSource {
    fn chart(&self) -> &PlateChart {
        match self { Self::Authored(c) => c, Self::Material(c) => c.chart() }
    }
}

/// Immutable reduction and original plate data. Modal sign/scaling cancels from
/// every physical coefficient and from the reconstructed opening-normalized shape.
#[derive(Debug, Clone)]
pub struct PlateApertureReduction {
    source: PlateSource,
    options: PlateApertureOptions,
    mass_kg: f64,
    stiffness_n_m: f64,
    area_m2: f64,
    width_m: f64,
    closing_pressure_pa: f64,
    lambda_interval: (f64, f64),
    in_window_modes: usize,
    shape: Vec<[f64; 3]>,
    slope_per_m: f64,
    slit_variation: f64,
}

fn map_plate(error: fs_plate::PlateError) -> AcousticRealizeError { AcousticRealizeError::Plate(error) }

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn checkpoint(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(invalid("plate aperture reduction cancelled")) } else { Ok(()) }
}

impl PlateApertureReduction {
    /// Reduce authored geometry/sections. A successful eigensolve certifies only
    /// its discrete eigenvalue statement, not omitted modes or physical fidelity.
    /// Cancellation is polled around the existing non-interruptible assembly and
    /// eigensolve, and during projection; a cancelled result is never published.
    ///
    /// # Errors
    /// Invalid physical inputs, slit geometry, exhausted work, unresolved modes,
    /// cancellation, non-closing pressure projection or approximation refusal.
    pub fn from_chart(chart: PlateChart, options: PlateApertureOptions, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        Self::reduce(PlateSource::Authored(chart), options, gate)
    }

    /// Same reduction with the exact immutable regional material receipts retained.
    /// Material state is the supplied resolved solid state, not ambient gas state.
    ///
    /// # Errors
    /// Same admission and solve refusals as [`Self::from_chart`].
    pub fn from_material_chart(chart: ResolvedPlateChart, options: PlateApertureOptions, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        Self::reduce(PlateSource::Material(chart), options, gate)
    }

    fn reduce(source: PlateSource, options: PlateApertureOptions, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        checkpoint(gate)?;
        let chart = source.chart();
        let n = chart.mesh.node_count();
        let (lo, hi) = options.eigenvalue_window;
        if n == 0 || n > options.max_nodes || chart.mesh.tris.len() > options.max_triangles
            || ![lo, hi, options.rest_opening_m, options.damping_ratio,
                options.max_slit_mode_variation, options.max_slope].iter().all(|x| x.is_finite())
            || lo < 0.0 || hi <= lo || options.rest_opening_m <= 0.0
            || options.damping_ratio < 0.0 || options.max_slit_mode_variation < 0.0
            || options.max_slope <= 0.0 || options.slit_edges.is_empty()
            || options.slit_edges.len() > n
        { return Err(invalid("plate aperture needs bounded geometry, physical gap/loss and explicit modal/slope limits")); }
        // Assembly owns mesh/section/support validation; do not inspect unchecked
        // triangle indices first. No duplicate plate element or eigensolver here.
        let model = chart.assemble(&[], &options.assembly).map_err(map_plate)?;
        checkpoint(gate)?;
        let boundary: std::collections::BTreeSet<_> = chart.mesh.boundary_edges().into_iter()
            .map(|(a, b)| (a.min(b), a.max(b))).collect();
        let mut seen = std::collections::BTreeSet::new();
        let mut weights = vec![0.0; n];
        let mut width = 0.0;
        for &[a, b] in &options.slit_edges {
            checkpoint(gate)?;
            let edge = (a.min(b), a.max(b));
            if a >= n || b >= n || a == b || !boundary.contains(&edge) || !seen.insert(edge) {
                return Err(invalid("slit edges must be unique existing plate boundary edges"));
            }
            let (x, y) = (chart.mesh.nodes[a], chart.mesh.nodes[b]);
            let length = (x.0 - y.0).hypot(x.1 - y.1);
            if !length.is_finite() || length <= 0.0 {
                return Err(invalid("slit edge length is not positive finite metres"));
            }
            width += length;
            weights[a] += 0.5 * length;
            weights[b] += 0.5 * length;
        }
        if !width.is_finite() || width <= 0.0 { return Err(invalid("slit width overflowed")); }
        for weight in &mut weights { *weight /= width; }
        let report = modes(&model, options.eigenvalue_window, &options.eigensolver).map_err(map_plate)?;
        checkpoint(gate)?;
        let pair = report.modes.get(options.mode_index)
            .ok_or_else(|| invalid("selected plate aperture mode is absent from the resolved window"))?;
        if !pair.lambda.is_finite() || pair.lambda <= 0.0 {
            return Err(invalid("plate aperture needs a positive restoring mode"));
        }
        if pair.phi.len() != model.free || pair.phi.iter().any(|x| !x.is_finite()) {
            return Err(invalid("plate aperture mode has invalid coordinates"));
        }
        let nodal_w = |node: usize| model.dof_map[3*node].map_or(0.0, |i| pair.phi[i]);
        let s: f64 = weights.iter().enumerate().map(|(i, w)| w * nodal_w(i)).sum();
        let mut m_phi = vec![0.0; model.free];
        model.m.spmv(&pair.phi, &mut m_phi);
        let modal_mass: f64 = pair.phi.iter().zip(m_phi).map(|(p, m)| p*m).sum();
        // Uniform face pressure acts on the P1 displacement trace of this DKT
        // mode. This is a mechanical virtual-work port, not a radiation filter.
        let mut modal_area = 0.0;
        for &[i, j, k] in &chart.mesh.tris {
            checkpoint(gate)?;
            let [a, b, c] = [i, j, k].map(|n| chart.mesh.nodes[n]);
            let da = 0.5 * ((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
            modal_area += da * (nodal_w(i)+nodal_w(j)+nodal_w(k))/3.0;
        }
        let mass = modal_mass / s / s;
        let stiffness = mass * pair.lambda;
        let area = modal_area / s;
        let closing = (stiffness / area) * options.rest_opening_m;
        if ![mass, stiffness, area, closing].iter().all(|x| x.is_finite() && *x > 0.0) {
            return Err(invalid("selected plate mode must move the slit and close under uniform positive pressure"));
        }
        let mut shape = Vec::with_capacity(n);
        let mut slope_per_m = 0.0_f64;
        let mut variation = 0.0_f64;
        for node in 0..n {
            checkpoint(gate)?;
            let mut value = [0.0; 3];
            for (c, component) in value.iter_mut().enumerate() {
                if let Some(r) = model.dof_map[3 * node + c] { *component = pair.phi[r] / s; }
            }
            slope_per_m = slope_per_m.max(value[1].hypot(value[2]));
            if weights[node] > 0.0 { variation = variation.max((value[0] - 1.0).abs()); }
            if !value.iter().all(|x| x.is_finite()) { return Err(invalid("normalized plate aperture shape overflowed")); }
            shape.push(value);
        }
        // Guard both DKT nodal rotations and the P1 trace used by pressure
        // quadrature. Checking rotations alone can miss a steep coarse triangle.
        for tri in &chart.mesh.tris {
            checkpoint(gate)?;
            let [a, b, c] = tri.map(|i| chart.mesh.nodes[i]);
            let [u, v, w] = tri.map(|i| shape[i][0]);
            let det = (b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0);
            let dx = ((v-u)*(c.1-a.1)-(w-u)*(b.1-a.1))/det;
            let dy = ((b.0-a.0)*(w-u)-(c.0-a.0)*(v-u))/det;
            let slope = dx.hypot(dy);
            if !slope.is_finite() { return Err(invalid("plate aperture trace slope overflowed")); }
            slope_per_m = slope_per_m.max(slope);
        }
        if !slope_per_m.is_finite() || slope_per_m <= 0.0 || variation > options.max_slit_mode_variation {
            return Err(invalid("plate aperture shape exceeds its slit-uniformity allowance or has no finite bending slope"));
        }
        checkpoint(gate)?;
        Ok(Self { source, options, mass_kg: mass, stiffness_n_m: stiffness, area_m2: area,
            width_m: width, closing_pressure_pa: closing, lambda_interval: pair.interval,
            in_window_modes: report.modes.len(), shape, slope_per_m, slit_variation: variation })
    }

    /// Actual numerical chart, including the section field and support nodes.
    #[must_use]
    pub fn chart(&self) -> &PlateChart { self.source.chart() }
    /// Original material-bound specimen, absent for authored numeric sections.
    #[must_use]
    pub fn material_chart(&self) -> Option<&ResolvedPlateChart> {
        match &self.source { PlateSource::Material(c) => Some(c), PlateSource::Authored(_) => None }
    }
    /// Complete caller choices retained separately from source-material authority.
    #[must_use]
    pub const fn options(&self) -> &PlateApertureOptions { &self.options }
    /// Effective mass in the mean-slit-opening coordinate [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Restoring stiffness in that same coordinate [N/m].
    #[must_use]
    pub const fn stiffness_n_m(&self) -> f64 { self.stiffness_n_m }
    /// Uniform-pressure force area, also used for swept volume [m²].
    #[must_use]
    pub const fn pressure_area_m2(&self) -> f64 { self.area_m2 }
    /// Actual sum of the supplied slit edge lengths [m].
    #[must_use]
    pub const fn width_m(&self) -> f64 { self.width_m }
    /// Static closing pressure of this ONE retained mode [Pa].
    #[must_use]
    pub const fn closing_pressure_pa(&self) -> f64 { self.closing_pressure_pa }
    /// Discrete squared-frequency interval; not an omitted-mode error bound.
    #[must_use]
    pub const fn eigenvalue_interval(&self) -> (f64, f64) { self.lambda_interval }
    /// Number resolved in the search window; only the selected mode is retained.
    #[must_use]
    pub const fn in_window_modes(&self) -> usize { self.in_window_modes }
    /// [w, rotation_x, rotation_y] per metre of mean slit opening change.
    #[must_use]
    pub fn shape_per_opening(&self) -> &[[f64; 3]] { &self.shape }
    /// Maximum slit-node departure from a uniform one-metre opening displacement.
    #[must_use]
    pub const fn slit_mode_variation(&self) -> f64 { self.slit_variation }
    /// Largest nodal/trace slope at this mean opening; zero at the rest geometry.
    #[must_use]
    pub fn max_slope_at(&self, opening_m: f64) -> f64 {
        (opening_m - self.options.rest_opening_m).abs() * self.slope_per_m
    }
    /// Check the declared linear-geometry domain without moving any state.
    ///
    /// # Errors
    /// Nonfinite opening or a nodal/trace slope above the declared limit.
    pub fn validate_opening(&self, opening_m: f64) -> Result<(), AcousticRealizeError> {
        let slope = self.max_slope_at(opening_m);
        if !opening_m.is_finite() || !slope.is_finite() || slope > self.options.max_slope {
            return Err(invalid("plate aperture exceeds its declared linear-slope domain"));
        }
        Ok(())
    }
    /// Reconstruct physical nodal motion from an accepted opening and velocity.
    /// This observes the mechanical state; it does not step or renormalize it.
    ///
    /// # Errors
    /// Absent node, invalid state, slope-domain refusal or nonfinite projection.
    pub fn nodal_motion(&self, node: usize, opening_m: f64, opening_velocity_m_s: f64)
        -> Result<PlateApertureMotion, AcousticRealizeError>
    {
        self.validate_opening(opening_m)?;
        let shape = self.shape.get(node).ok_or_else(|| invalid("plate aperture node is outside the source mesh"))?;
        let displacement = opening_m - self.options.rest_opening_m;
        let motion = PlateApertureMotion {
            displacement_rotation: shape.map(|x| x * displacement),
            velocity_rotation_rate: shape.map(|x| x * opening_velocity_m_s),
        };
        if !opening_velocity_m_s.is_finite() || !motion.displacement_rotation.iter()
            .chain(&motion.velocity_rotation_rate).all(|x| x.is_finite())
        { return Err(invalid("plate aperture physical nodal motion overflowed")); }
        Ok(motion)
    }

    /// Supply only the independent gas/load/clock budget to existing mechanics.
    /// No geometry-derived coefficient can be independently retuned here.
    /// The dynamic owner performs the full fluid/clock admission at construction.
    #[must_use]
    pub fn dynamic_spec(&self, density_kg_m3: f64, impedance_pa_s_m3: f64,
        time_step_s: f64, max_steps: u64) -> DynamicApertureSpec
    {
        DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: self.options.rest_opening_m,
                width_m: self.width_m, closing_pressure_pa: self.closing_pressure_pa },
            mass_kg: self.mass_kg, stiffness_n_m: self.stiffness_n_m,
            damping_ratio: self.options.damping_ratio, density_kg_m3,
            impedance_pa_s_m3, time_step_s, max_steps,
        }
    }
}
