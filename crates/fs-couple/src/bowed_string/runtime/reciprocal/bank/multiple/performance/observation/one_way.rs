//! One-way finite-distance observation of a material shell through its own
//! declared closed acoustic skin.
//!
//! This is an output owner only. No air reaction enters the mechanics, so the
//! accepted network may be the exact or the prepared-linear owner. The skin is
//! a separate closed outward surface of the same body. Each skin vertex moves
//! as a declared convex combination of points on the original parts, and is
//! checked against their positions; a panel's normal velocity per mode is the
//! mean of its vertices', the rule the structural facet trace uses. Aperture
//! caps carry the lumped cavity's outward neck volume flow in declared shares
//! and no structural motion.
//!
//! The transfer is baked on one ruled Burton--Miller operator at every
//! frequency. The geometry-routed owner switches to it from plain CBIE at the
//! skin's box bound, and that switch is a step in every transfer the DFT ladder
//! cannot converge across.
//!
//! Because nothing is fed back, the report states for every retained mode the
//! radiation loss and added mass this owner leaves out, against the mode's own
//! wood loss, bandwidth and spacing. Those numbers come from the same skin and
//! rule on the geometry-routed owner at each mode's own frequency (plain CBIE
//! below the bound, where Burton--Miller overstates low-ka resistance); they
//! are estimates on that skin, not a convergence certificate.
//!
//! The transfer and the receipt depend only on the transfer declaration and the
//! body-mode frequencies, so a [`OneWaySetupStore`] may keep them between
//! renders of one body; see `one_way_store.rs`.

use super::super::MultiBowPerformance;
use super::{BoundReport, BoundShellObserver, ShellObserverError, invalid, poll};
use fs_bem::helmholtz::{self, Medium, PanelRule};
use fs_bem::observation_transfers::{
    BemObservation, MAX_OBSERVATION_RUNGS, ObservationTransferBake, ProducerUncertainty,
    ReceiverSpec, SourceId, TransferDeclaration, TransferError,
};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_plate::shell::ShellMesh;
use fs_vfit::observation::MAX_OBSERVATION_ENTRIES;

#[path = "one_way_store.rs"]
mod setup_store;
pub use setup_store::{MAX_ONE_WAY_SETUP_BYTES, OneWaySetupStore};

/// Complete dense boundary work this owner admits, over every bake and bound
/// frequency: panels^2 * (panels + sources + 1) per frequency.
pub const MAX_ONE_WAY_BOUNDARY_TERMS: u64 = 1 << 50;

/// One point on an original part: barycentric coordinates in one of its
/// triangles. A skin vertex is a positive-weight combination of such points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SkinAnchor {
    pub weight: f64,
    pub component: usize,
    pub triangle: usize,
    pub barycentric: [f64; 3],
}

/// One cap over a sound hole. Its panels carry `flow_share` of the cavity's
/// outward neck volume flow, uniformly, and no structural motion.
#[derive(Clone, Debug, PartialEq)]
pub struct SkinAperture {
    pub id: String,
    pub triangles: Vec<usize>,
    pub flow_share: f64,
}

/// A closed, outward-oriented acoustic surface of the shell in the body frame.
/// Every vertex touched by a structural (non-aperture) triangle needs anchors.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellAcousticSkin {
    pub vertices: Vec<[f64; 3]>,
    pub anchors: Vec<Vec<SkinAnchor>>,
    pub triangles: Vec<[usize; 3]>,
    pub apertures: Vec<SkinAperture>,
}

/// Content identity of a declared skin: vertices, anchors, triangles and
/// apertures, bit for bit. A two-way load records it so that the observer can
/// confirm it reads the same geometry.
pub fn skin_identity(skin: &ShellAcousticSkin) -> ContentHash {
    let mut h = fs_blake3::DomainHasher::new("fs-couple:shell-acoustic-skin:v1");
    h.update(&(skin.vertices.len() as u64).to_le_bytes());
    for (vertex, anchors) in skin.vertices.iter().zip(&skin.anchors) {
        for x in vertex {
            h.update(&x.to_bits().to_le_bytes());
        }
        h.update(&(anchors.len() as u64).to_le_bytes());
        for a in anchors {
            h.update(&a.weight.to_bits().to_le_bytes());
            h.update(&(a.component as u64).to_le_bytes());
            h.update(&(a.triangle as u64).to_le_bytes());
            for x in a.barycentric {
                h.update(&x.to_bits().to_le_bytes());
            }
        }
    }
    h.update(&(skin.triangles.len() as u64).to_le_bytes());
    for t in &skin.triangles {
        for v in t {
            h.update(&(*v as u64).to_le_bytes());
        }
    }
    h.update(&(skin.apertures.len() as u64).to_le_bytes());
    for a in &skin.apertures {
        h.update(&(a.id.len() as u64).to_le_bytes());
        h.update(a.id.as_bytes());
        h.update(&(a.triangles.len() as u64).to_le_bytes());
        for t in &a.triangles {
            h.update(&(*t as u64).to_le_bytes());
        }
        h.update(&a.flow_share.to_bits().to_le_bytes());
    }
    h.finalize()
}

/// When a mode may stay one-way: the radiation damping it leaves out is at most
/// `maximum_loss_fraction` of its wood damping, and its added-mass frequency
/// shift is at most `maximum_shift_fraction` of both its wood half-power
/// bandwidth and its spacing to the nearest retained mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OneWayRule {
    pub maximum_loss_fraction: f64,
    pub maximum_shift_fraction: f64,
}

/// Complete receiver, skin quadrature, realization and work declaration.
#[derive(Clone, Debug, PartialEq)]
pub struct OneWayShellObserverOptions {
    pub receiver: ReceiverSpec,
    pub minimum_clearance_m: f64,
    pub pass_band_upper_hz: f64,
    pub taper_stop_hz: f64,
    pub filter_latency_samples: usize,
    pub fft_lengths: Vec<usize>,
    pub tail_energy_tolerance: f64,
    pub dtft_relative_tolerance: f64,
    pub maximum_frequency_entries: usize,
    /// Retained FIR multiply-add terms per accepted sample, all four outputs.
    pub maximum_replay_terms: u64,
    /// Boundary quadrature of the transfer and the receipt; it must integrate triangles.
    pub panel_rule: PanelRule,
    pub maximum_skin_panels: usize,
    /// Largest admitted distance from a skin vertex to its anchored point [m].
    pub maximum_anchor_distance_m: f64,
    /// Complete dense boundary work over every bake and bound frequency.
    pub maximum_boundary_terms: u64,
    /// Threads for dense assembly and LU; no bit depends on the count.
    pub workers: usize,
    pub rule: OneWayRule,
}

/// What this owner leaves out for one retained mode, at the mode's frequency.
#[derive(Clone, Debug, PartialEq)]
pub struct OneWayModeBound {
    pub frequency_hz: f64,
    /// The mode's own damping ratio plus its diagonal wood viscosity.
    pub wood_zeta: f64,
    /// Re Z_nn / (2 omega_n): the radiation damping ratio not fed back.
    pub radiation_zeta: f64,
    /// -M_nn / 2 with M_nn = -Im Z_nn / omega_n: the relative frequency change.
    pub added_mass_shift: f64,
    /// 2 wood_zeta f_n [Hz].
    pub half_power_bandwidth_hz: f64,
    /// Distance to the nearest other retained mode [Hz]; infinite for one mode.
    pub nearest_spacing_hz: f64,
    /// Largest omega_n |Z_nm| / |omega_n^2 - omega_m^2| over the other modes:
    /// the first-order air mixing of mode m into mode n.
    pub largest_mixing: f64,
    pub largest_mixing_mode: Option<usize>,
    pub one_way: bool,
}

/// Retained skin, realization and per-mode one-way receipt. Estimates only.
#[derive(Clone, Debug)]
pub struct OneWayShellObserverReport {
    pub options: OneWayShellObserverOptions,
    pub gas: GasState,
    pub body_modes: usize,
    pub apertures: usize,
    pub surface_panels: usize,
    pub largest_anchor_distance_m: f64,
    /// Below this wavenumber the receipt uses plain CBIE [rad/m].
    pub plain_cbie_limit_rad_per_m: f64,
    pub boundary_terms: u64,
    pub replay_terms: u64,
    pub modes: Vec<OneWayModeBound>,
    /// The store key of this setup: the transfer declaration and every body-mode frequency.
    pub setup_key: ContentHash,
    /// The transfer and the receipt were loaded from a store, not computed.
    pub setup_loaded: bool,
    /// The bound two-way load on this same skin, if any: its loading identity.
    /// Each receipt row then states what a one-way run of that mode would leave
    /// out; the load feeds it back.
    pub two_way_load: Option<ContentHash>,
    /// The ascending prefix of body modes that load feeds back; zero without it.
    pub two_way_modes: usize,
    transfer: ObservationTransferBake,
}
impl OneWayShellObserverReport {
    pub const fn transfer(&self) -> &ObservationTransferBake {
        &self.transfer
    }
    /// Indices of the modes outside the declared one-way rule.
    pub fn fed_back_modes(&self) -> Vec<usize> {
        self.modes
            .iter()
            .enumerate()
            .filter(|(_, m)| !m.one_way)
            .map(|(i, _)| i)
            .collect()
    }
}

impl OneWayShellObserverOptions {
    /// Numerical declaration admission, before any skin, shell or pack work.
    pub fn validate(&self, rate: u32) -> Result<(), ShellObserverError> {
        let finite_positive = |x: f64| x.is_finite() && x > 0.0;
        let fraction = |x: f64| x.is_finite() && x >= 0.0;
        if rate == 0
            || self.receiver.id.is_empty()
            || self.receiver.id.len() > 1024
            || self.receiver.position_m.iter().any(|x| !x.is_finite())
            || self.fft_lengths.len() > MAX_OBSERVATION_RUNGS
            || !finite_positive(self.minimum_clearance_m)
            || !finite_positive(self.maximum_anchor_distance_m)
            || !(4..=helmholtz::MAX_DENSE_PANELS).contains(&self.maximum_skin_panels)
            || !(1..=MAX_ONE_WAY_BOUNDARY_TERMS).contains(&self.maximum_boundary_terms)
            || self.maximum_replay_terms == 0
            || self.maximum_replay_terms > MAX_OBSERVATION_ENTRIES as u64
            || !(1..=256).contains(&self.workers)
            || !fraction(self.rule.maximum_loss_fraction)
            || !fraction(self.rule.maximum_shift_fraction)
            || matches!(
                self.panel_rule,
                PanelRule::Centroid
                    | PanelRule::TriangleNear { max_pairs: 0 }
                    | PanelRule::Triangle { max_pairs: 0 }
            )
        {
            return Err(invalid(
                "invalid one-way receiver, skin, rule or work declaration",
            ));
        }
        self.declaration(
            rate,
            vec![SourceId::StructuralMode("numerical admission only".into())],
            vec![0.0],
            fs_blake3::ContentHash([1; 32]),
        )
        .admit_work()?;
        let finest = self
            .fft_lengths
            .last()
            .copied()
            .ok_or_else(|| invalid("missing observer FFT ladder"))?;
        if self.pass_band_upper_hz < f64::from(rate) / (finest as f64) {
            return Err(invalid(
                "observer pass band must contain a resolved positive FFT bin",
            ));
        }
        Ok(())
    }

    fn declaration(
        &self,
        rate: u32,
        sources: Vec<SourceId>,
        geometric_delays_s: Vec<f64>,
        producer_hash: fs_blake3::ContentHash,
    ) -> TransferDeclaration {
        TransferDeclaration {
            producer_hash,
            sources,
            receivers: vec![self.receiver.clone()],
            fs_hz: f64::from(rate),
            geometric_delays_s,
            pass_band_upper_hz: self.pass_band_upper_hz,
            taper_stop_hz: self.taper_stop_hz,
            filter_latency_samples: self.filter_latency_samples,
            fft_lengths: self.fft_lengths.clone(),
            tail_energy_tolerance: self.tail_energy_tolerance,
            dtft_relative_tolerance: self.dtft_relative_tolerance,
            producer_uncertainty: ProducerUncertainty::Unresolved {
                source_hash: producer_hash,
            },
            maximum_frequency_entries: self.maximum_frequency_entries,
        }
    }
}

/// The skin's panels and its per-mode and per-aperture normal traces.
pub(in crate::bowed_string::runtime::reciprocal::bank) struct SkinTraces {
    pub surface: SpherePanels,
    /// Mode-major structural traces [m/s per unit modal velocity].
    pub modes: Vec<Vec<f64>>,
    /// Unit-flow aperture traces, 1/A on the cap and zero elsewhere.
    pub apertures: Vec<Vec<f64>>,
    pub largest_anchor_distance_m: f64,
}

/// Bind the skin to the original parts and compute every normal trace.
/// `parts` are the original midsurfaces with their mode-major six-DOF rows.
pub(in crate::bowed_string::runtime::reciprocal::bank) fn skin_traces(
    parts: &[(&ShellMesh, &[Vec<f64>])],
    modes: usize,
    skin: &ShellAcousticSkin,
    maximum_panels: usize,
    maximum_anchor_distance_m: f64,
    gate: &CancelGate,
) -> Result<SkinTraces, ShellObserverError> {
    let vertices = skin.vertices.len();
    if modes == 0
        || skin.triangles.len() < 4
        || skin.triangles.len() > maximum_panels
        || skin.anchors.len() != vertices
        || skin.vertices.iter().flatten().any(|x| !x.is_finite())
        || skin.triangles.iter().flatten().any(|&v| v >= vertices)
    {
        return Err(invalid(
            "acoustic skin needs finite vertices, anchors per vertex and admitted panels",
        ));
    }
    for (mesh, rows) in parts {
        let nodes = mesh.nodes.len();
        if rows.len() != modes
            || rows
                .iter()
                .any(|r| r.len() != 6 * nodes || r.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid(
                "acoustic skin must observe the complete original shell basis",
            ));
        }
    }
    // Aperture membership: every listed triangle exactly once, shares positive and summing to one.
    let mut aperture_of = vec![None; skin.triangles.len()];
    let mut share_sum = 0.0;
    for (index, aperture) in skin.apertures.iter().enumerate() {
        if aperture.id.is_empty()
            || aperture.triangles.is_empty()
            || !aperture.flow_share.is_finite()
            || aperture.flow_share <= 0.0
        {
            return Err(invalid(
                "acoustic skin aperture needs an id, triangles and a positive flow share",
            ));
        }
        for &t in &aperture.triangles {
            match aperture_of.get_mut(t) {
                Some(slot @ None) => *slot = Some(index),
                _ => {
                    return Err(invalid(
                        "acoustic skin aperture triangle is missing or repeated",
                    ));
                }
            }
        }
        share_sum += aperture.flow_share;
    }
    if !skin.apertures.is_empty() && (share_sum - 1.0).abs() > 1e-12 {
        return Err(invalid(
            "acoustic skin aperture flow shares must sum to one",
        ));
    }
    let mut structural_vertex = vec![false; vertices];
    for (t, tri) in skin.triangles.iter().enumerate() {
        if aperture_of[t].is_none() {
            for &v in tri {
                structural_vertex[v] = true;
            }
        }
    }
    // Each vertex's modal translation: the anchored combination of part points.
    let mut largest = 0.0_f64;
    let mut translation = vec![[0.0_f64; 3]; vertices * modes];
    for (v, anchors) in skin.anchors.iter().enumerate() {
        poll(gate)?;
        if anchors.is_empty() {
            if structural_vertex[v] {
                return Err(invalid(
                    "a structural skin vertex has no anchor on the original parts",
                ));
            }
            continue;
        }
        let weights: f64 = anchors.iter().map(|a| a.weight).sum();
        if (weights - 1.0).abs() > 1e-12 {
            return Err(invalid("skin anchor weights must sum to one"));
        }
        let mut position = [0.0_f64; 3];
        for anchor in anchors {
            let (mesh, rows) = parts
                .get(anchor.component)
                .ok_or_else(|| invalid("skin anchor names a missing part"))?;
            let nodes = *mesh
                .tris
                .get(anchor.triangle)
                .ok_or_else(|| invalid("skin anchor names a missing triangle"))?;
            let b = anchor.barycentric;
            if !anchor.weight.is_finite()
                || anchor.weight <= 0.0
                || b.iter().any(|x| !x.is_finite() || *x < -1e-12)
                || (b[0] + b[1] + b[2] - 1.0).abs() > 1e-12
            {
                return Err(invalid(
                    "skin anchor needs a positive weight and barycentric coordinates",
                ));
            }
            for (corner, &node) in nodes.iter().enumerate() {
                let scale = anchor.weight * b[corner];
                for axis in 0..3 {
                    position[axis] += scale * mesh.nodes[node][axis];
                }
                for (mode, row) in rows.iter().enumerate() {
                    let entry = &mut translation[v * modes + mode];
                    for axis in 0..3 {
                        entry[axis] += scale * row[6 * node + axis];
                    }
                }
            }
        }
        let distance = (0..3)
            .map(|i| (position[i] - skin.vertices[v][i]).powi(2))
            .sum::<f64>()
            .sqrt();
        if !distance.is_finite() || distance > maximum_anchor_distance_m {
            return Err(invalid(
                "a skin vertex lies farther from its anchored point than declared",
            ));
        }
        largest = largest.max(distance);
    }
    let triangles: Vec<_> = skin
        .triangles
        .iter()
        .map(|t| t.map(|v| skin.vertices[v]))
        .collect();
    let surface = SpherePanels::from_triangles(triangles)
        .map_err(|_| invalid("acoustic skin triangles are degenerate"))?;
    let mut traces = vec![vec![0.0; surface.areas().len()]; modes];
    for (t, (tri, normal)) in skin.triangles.iter().zip(surface.normals()).enumerate() {
        if aperture_of[t].is_some() {
            continue;
        }
        for (mode, trace) in traces.iter_mut().enumerate() {
            let mut value = 0.0;
            for axis in 0..3 {
                let mean = tri
                    .iter()
                    .map(|&v| translation[v * modes + mode][axis])
                    .sum::<f64>()
                    / 3.0;
                value += normal[axis] * mean;
            }
            if !value.is_finite() {
                return Err(invalid("nonfinite acoustic skin normal trace"));
            }
            trace[t] = value;
        }
    }
    let mut apertures = Vec::with_capacity(skin.apertures.len());
    for aperture in &skin.apertures {
        let area: f64 = aperture.triangles.iter().map(|&t| surface.areas()[t]).sum();
        if !area.is_finite() || area <= 0.0 {
            return Err(invalid("acoustic skin aperture has no area"));
        }
        let mut trace = vec![0.0; surface.areas().len()];
        for &t in &aperture.triangles {
            trace[t] = 1.0 / area;
        }
        apertures.push(trace);
    }
    Ok(SkinTraces {
        surface,
        modes: traces,
        apertures,
        largest_anchor_distance_m: largest,
    })
}

/// The one-way receipt from each mode's modal radiation impedance row
/// `z[n][m] = sum_j p_m(j) trace_n(j) A_j` at omega_n (exp(-i omega t)).
pub(super) fn mode_bounds(
    omega: &[f64],
    wood_zeta: &[f64],
    z: &[Vec<C64>],
    rule: OneWayRule,
) -> Result<Vec<OneWayModeBound>, ShellObserverError> {
    let n = omega.len();
    if n == 0
        || wood_zeta.len() != n
        || z.len() != n
        || z.iter().any(|row| row.len() != n)
        || omega.iter().any(|w| !w.is_finite() || *w <= 0.0)
        || wood_zeta.iter().any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(invalid(
            "one-way receipt needs positive modes, wood loss and a square impedance",
        ));
    }
    let tau = core::f64::consts::TAU;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let w = omega[i];
        let zii = z[i][i];
        let radiation_zeta = zii.re / (2.0 * w);
        let added_mass = -zii.im / w;
        let added_mass_shift = -0.5 * added_mass;
        let f = w / tau;
        let half_power_bandwidth_hz = 2.0 * wood_zeta[i] * f;
        let mut nearest_spacing_hz = f64::INFINITY;
        let mut largest_mixing = 0.0_f64;
        let mut largest_mixing_mode = None;
        for m in (0..n).filter(|&m| m != i) {
            nearest_spacing_hz = nearest_spacing_hz.min((omega[m] - w).abs() / tau);
            let detuning = (w * w - omega[m] * omega[m]).abs();
            let mixing = if detuning == 0.0 {
                f64::INFINITY
            } else {
                w * z[i][m].abs() / detuning
            };
            if mixing > largest_mixing || largest_mixing_mode.is_none() {
                largest_mixing = mixing;
                largest_mixing_mode = Some(m);
            }
        }
        if !radiation_zeta.is_finite()
            || !added_mass_shift.is_finite()
            || !largest_mixing.is_finite()
        {
            return Err(invalid("nonfinite one-way radiation receipt"));
        }
        let shift_hz = (added_mass_shift * f).abs();
        let one_way = radiation_zeta.abs() <= rule.maximum_loss_fraction * wood_zeta[i]
            && shift_hz
                <= rule.maximum_shift_fraction * half_power_bandwidth_hz.min(nearest_spacing_hz);
        out.push(OneWayModeBound {
            frequency_hz: f,
            wood_zeta: wood_zeta[i],
            radiation_zeta,
            added_mass_shift,
            half_power_bandwidth_hz,
            nearest_spacing_hz,
            largest_mixing,
            largest_mixing_mode,
            one_way,
        });
    }
    Ok(out)
}

/// Row n of the modal radiation impedance at omega_n, on the routed owner's own
/// formulation and rule: every mode's surface pressure projected on mode n's
/// trace, `sum_j p_m(j) trace_n(j) A_j`. Also returns the plain-CBIE limit.
pub(super) fn impedance_rows(
    surface: &SpherePanels,
    traces: &[Vec<f64>],
    omega: &[f64],
    medium: Medium,
    rule: PanelRule,
    workers: usize,
    gate: &CancelGate,
) -> Result<(Vec<Vec<C64>>, f64), ShellObserverError> {
    if traces.len() != omega.len() || traces.iter().any(|t| t.len() != surface.areas().len()) {
        return Err(invalid(
            "impedance rows need one complete trace per mode frequency",
        ));
    }
    let policy = GeometryPolicy::new(surface).map_err(TransferError::Helmholtz)?;
    let complex: Vec<Vec<C64>> = traces
        .iter()
        .map(|t| t.iter().copied().map(C64::from_re).collect())
        .collect();
    let areas = surface.areas();
    let mut rows = Vec::with_capacity(omega.len());
    for (n, &w) in omega.iter().enumerate() {
        poll(gate)?;
        let k = w / medium.sound_speed;
        let formulation = policy.formulation(k).map_err(TransferError::Helmholtz)?;
        let mut row = Vec::with_capacity(omega.len());
        for batch in complex.chunks(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH) {
            let refs: Vec<&[C64]> = batch.iter().map(Vec::as_slice).collect();
            let solved = helmholtz::with_dense_workers(workers, || {
                helmholtz::solve_radiation_batch_with_rule(
                    surface,
                    k,
                    medium,
                    &refs,
                    formulation,
                    rule,
                )
            })
            .map_err(TransferError::Helmholtz)?;
            row.extend(solved.iter().map(|s| {
                s.pressure
                    .iter()
                    .zip(&traces[n])
                    .zip(areas)
                    .fold(C64::ZERO, |sum, ((&p, &t), &a)| sum + p.scale(t * a))
            }));
        }
        rows.push(row);
    }
    Ok((rows, policy.plain_cbie_limit()))
}

/// Dense work per frequency: assembly, one LU and every right-hand side.
fn frequency_terms(panels: usize, sources: usize) -> Option<u64> {
    let p = panels as u64;
    p.checked_mul(p)?
        .checked_mul(p.checked_add(sources as u64)?.checked_add(1)?)
}

impl MultiBowPerformance {
    /// Replace compact pressure with one finite-distance receiver of the material
    /// shell through `skin`, without feeding any air reaction back. The accepted
    /// mechanics may be prepared-linear. A performance with a computed air load
    /// keeps that load's own observer; a spatial cavity keeps its own mouths.
    /// With a `store`, the transfer and the receipt come from the entry kept
    /// under this setup's key, or are computed and then kept there.
    ///
    /// # Errors
    /// Declaration, skin binding, boundary-work, realization-gate and receipt
    /// refusals; a bank that has accepted any sample; a damaged kept entry and
    /// the store's own failures.
    pub fn with_one_way_shell_observer(
        mut self,
        options: OneWayShellObserverOptions,
        skin: &ShellAcousticSkin,
        store: Option<&dyn OneWaySetupStore>,
        gate: &CancelGate,
    ) -> Result<Self, ShellObserverError> {
        poll(gate)?;
        if self.poisoned
            || self.shell_observer.is_some()
            || self.base_sample != 0
            || self.bank.system.samples_rendered() != 0
        {
            return Err(invalid(
                "shell observation requires one unstarted physical history",
            ));
        }
        let source = self.bank.source.bank();
        let rate = source.options().sample_rate_hz;
        options.validate(rate)?;
        // The skin maps skin velocity to the receiver whether or not the air
        // reacts, so a two-way load on this same skin keeps this observer.
        let (two_way_load, two_way_modes) = match source.body_radiation() {
            None => (None, 0),
            Some(binding) => {
                let loading = binding
                    .material_shell_bem()
                    .filter(|r| r.skin_identity == Some(skin_identity(skin)))
                    .ok_or_else(|| {
                        invalid("a computed air load is observed on its own skin or by its own loaded shell observer")
                    })?;
                (Some(loading.loading_identity), loading.body_modes)
            }
        };
        if source.spatial_cavity().is_some() {
            return Err(invalid(
                "one-way skin observation has no spatial cavity mouth binding",
            ));
        }
        let mut parts = Vec::new();
        let (damping, mode_count) = if let Some(image) = source.material_shell_image() {
            parts.push((image.source().mesh(), image.physical_modes()));
            (image.damping(), image.physical_modes().len())
        } else if let Some(image) = source.material_shell_assembly_image() {
            for (part, rows) in image.source().components.iter().zip(image.physical_modes()) {
                parts.push((part.specimen.mesh(), rows.as_slice()));
            }
            (
                image.damping(),
                image.physical_modes().first().map_or(0, Vec::len),
            )
        } else {
            return Err(invalid(
                "one-way skin observation requires an original material shell basis",
            ));
        };
        let body_component = source.strings().len();
        let body = self
            .bank
            .system
            .components()
            .get(body_component)
            .ok_or_else(|| invalid("missing original shell body component"))?;
        if body.modes().len() != mode_count {
            return Err(invalid(
                "shell observation must retain every original body mode in order",
            ));
        }
        let cavity = source.material_cavity();
        if cavity
            .is_some_and(|c| c.cavity.loss() == crate::cavity_modal::CavityNeckLoss::ExteriorMouth)
            && two_way_load.is_none()
        {
            return Err(invalid(
                "an exterior-owned mouth has no outer mass without its two-way skin load",
            ));
        }
        if cavity.is_some() != !skin.apertures.is_empty() {
            return Err(invalid(
                "a lumped cavity radiates through declared skin apertures, and only then",
            ));
        }
        let gas = *source.gas();
        gas.validate_coherence()
            .map_err(|_| invalid("incoherent original shell observation gas"))?;
        let traces = skin_traces(
            &parts,
            mode_count,
            skin,
            options.maximum_skin_panels,
            options.maximum_anchor_distance_m,
            gate,
        )?;
        let surface_panels = traces.surface.areas().len();
        let finest = options
            .fft_lengths
            .last()
            .copied()
            .ok_or_else(|| invalid("missing observer FFT ladder"))?;
        let source_count = mode_count + skin.apertures.len();
        let replay_terms = (source_count as u64)
            .checked_mul(4)
            .and_then(|n| n.checked_mul((finest / 2) as u64))
            .ok_or_else(|| invalid("complete observation replay extent overflow"))?;
        if replay_terms > options.maximum_replay_terms {
            return Err(invalid(
                "complete shell observer exceeds its per-sample replay allowance",
            ));
        }
        // Every positive bin below the stop band on the finest grid, plus one solve per mode.
        let frequency_samples =
            (options.taper_stop_hz * finest as f64 / f64::from(rate)).ceil() as usize;
        let boundary_terms = frequency_terms(surface_panels, source_count)
            .and_then(|t| t.checked_mul((frequency_samples + mode_count) as u64))
            .ok_or_else(|| invalid("one-way boundary work overflow"))?;
        if boundary_terms > options.maximum_boundary_terms {
            return Err(invalid("one-way skin exceeds its declared boundary work"));
        }
        let mut sources: Vec<_> = (0..mode_count)
            .map(|m| SourceId::StructuralMode(format!("one-way-shell-mode:{m}")))
            .collect();
        let mut coordinates: Vec<_> = (0..mode_count).map(|m| (body_component, m, 1.0)).collect();
        let mut fields = traces.modes.clone();
        if let Some(cavity) = cavity {
            let shape = cavity.cavity.neck_volume_shape_m2_per_sqrt_kg();
            let neck = self
                .bank
                .system
                .components()
                .get(body_component + 1)
                .ok_or_else(|| invalid("one-way observer lost the cavity neck component"))?;
            if neck.modes().len() != 1 || !shape.is_finite() || shape <= 0.0 {
                return Err(invalid(
                    "one-way observer needs the one physical neck coordinate",
                ));
            }
            for (aperture, trace) in skin.apertures.iter().zip(&traces.apertures) {
                sources.push(SourceId::ApertureVolumeFlow(format!(
                    "one-way-shell-aperture:{}",
                    aperture.id
                )));
                coordinates.push((body_component + 1, 0, aperture.flow_share * shape));
                fields.push(trace.clone());
            }
        }
        let medium = Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        let owner = BemObservation::new_ruled_burton_miller(
            &traces.surface,
            medium,
            &sources,
            &fields,
            core::slice::from_ref(&options.receiver),
            options.minimum_clearance_m,
            options.panel_rule,
        )?;
        let producer = owner.operator_hash();
        let declaration = options.declaration(rate, sources, owner.geometric_delays_s(), producer);
        declaration.admit_work()?;
        let omega: Vec<f64> = body
            .modes()
            .iter()
            .map(|m| m.angular_frequency_rad_s)
            .collect();
        let setup_key = setup_store::setup_key(&declaration, &omega);
        let kept = match store {
            Some(store) => store.load(&setup_key).map_err(ShellObserverError::Store)?,
            None => None,
        };
        let setup_loaded = kept.is_some();
        let (transfer, impedance, plain_cbie_limit_rad_per_m) = if let Some(entry) = kept {
            setup_store::KeptSetup::decode(&entry, &setup_key, mode_count)?.into_parts(&declaration)
        } else {
            let transfer = helmholtz::with_dense_workers(options.workers, || {
                owner.bake(&declaration, || gate.is_requested())
            })?;
            poll(gate)?;
            // The receipt: each mode's impedance row at its own frequency, on the same skin and rule.
            let (impedance, limit) = impedance_rows(
                &traces.surface,
                &traces.modes,
                &omega,
                medium,
                options.panel_rule,
                options.workers,
                gate,
            )?;
            (transfer, impedance, limit)
        };
        let observer = transfer.observer()?;
        poll(gate)?;
        let matrix = damping.matrix_per_s();
        let wood_zeta: Vec<f64> = body
            .modes()
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let viscous = if damping.mode_count() == mode_count {
                    matrix[i * mode_count + i]
                } else {
                    0.0
                };
                m.damping_ratio + viscous / (2.0 * m.angular_frequency_rad_s)
            })
            .collect();
        let modes = mode_bounds(&omega, &wood_zeta, &impedance, options.rule)?;
        poll(gate)?;
        if let (Some(store), false) = (store, setup_loaded) {
            let entry = setup_store::encode(
                &setup_key,
                &transfer,
                &impedance,
                plain_cbie_limit_rad_per_m,
            );
            store
                .keep(&setup_key, &entry)
                .map_err(ShellObserverError::Store)?;
        }
        let report = OneWayShellObserverReport {
            gas,
            body_modes: mode_count,
            apertures: skin.apertures.len(),
            surface_panels,
            largest_anchor_distance_m: traces.largest_anchor_distance_m,
            plain_cbie_limit_rad_per_m,
            boundary_terms,
            replay_terms,
            modes,
            setup_key,
            setup_loaded,
            two_way_load,
            two_way_modes,
            transfer,
            options,
        };
        self.shell_observer = Some(BoundShellObserver {
            report: BoundReport::OneWay(Box::new(report)),
            observer: super::BoundReceiver::Shell(observer),
            coordinates,
            velocity: vec![0.0; source_count],
            maximum_abs_pressure_pa: source.options().network.maximum_abs_pressure_pa,
            samples: 0,
        });
        Ok(self)
    }

    /// The bound one-way skin receiver and its per-mode receipt.
    pub fn one_way_shell_observer(&self) -> Option<&OneWayShellObserverReport> {
        match self.shell_observer.as_ref().map(|o| &o.report) {
            Some(BoundReport::OneWay(report)) => Some(report.as_ref()),
            _ => None,
        }
    }
}

#[cfg(test)]
#[path = "one_way_tests.rs"]
mod tests;
