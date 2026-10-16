//! Geometry -> FE heads (fs-plate) -> certified eigenmodes (fs-modal).
//!
//! This is the minimal `drum-modal` construction of the offline reference
//! (crates/fs-couple/examples/percussion/main.rs `drum_with_barrel` with
//! `prepared=true`, and drum_spec.rs `Spec::prepare`), lowered through the
//! SAME library owners with the SAME physical cards and sign conventions:
//! both head coordinates are positive downward, and the signed modal areas
//! integrate cavity COMPRESSION (+ for the batter head, - for the resonant).

use crate::params::{DrumParams, HeadParams, MAX_HEAD_MODES, Refusal};
use fs_couple::render::plate::impact::striker::{RadiusStation, StrikerProperties};
use fs_plate::shell::head::{TensionedDisk, TensionedDiskSpec};
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::{ModePair, SliceOptions};

/// Launch gap between the stick tip and the membrane surface [m] (reference: 0.2 mm).
pub const LAUNCH_GAP_M: f64 = 0.0002;
/// Hertz exponent of the reference's estimated tip.
pub const HERTZ_ALPHA: f64 = 1.5;

/// One prepared head: the film pencil, its retained modes and derived rows.
pub struct Head {
    /// The actual DKT + membrane-prestress pencil and mesh.
    pub film: TensionedDisk,
    /// Mass-normalized eigenpairs inside the declared window, ascending.
    pub modes: Vec<ModePair>,
    /// Natural angular frequencies [rad/s] = sqrt(lambda).
    pub omegas: Vec<f64>,
    /// Viscous coefficients 2 zeta omega [1/s].
    pub damping_per_s: Vec<f64>,
    /// Unsigned modal surface integrals [m^2/sqrt(kg)] (`TensionedDisk::modal_area`).
    pub areas: Vec<f64>,
    /// Transverse displacement per unit modal coordinate at each mesh node,
    /// node-major `[node * modes + mode]` [1/sqrt(kg)]. Rim nodes are zero.
    pub node_shapes: Vec<f64>,
    /// Largest certified eigen-residual bound in the window [rad^2/s^2].
    pub max_residual: f64,
}
impl Head {
    /// Retained mode count.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }
    /// In-vacuo natural frequencies [Hz].
    #[must_use]
    pub fn frequencies_hz(&self) -> Vec<f64> {
        self.omegas.iter().map(|w| w / core::f64::consts::TAU).collect()
    }
}

/// The estimated stick and its contact card (identical to the reference).
#[derive(Debug, Clone, Copy)]
pub struct StickCard {
    /// Whole stick mass [kg].
    pub full_mass_kg: f64,
    /// Rigid-rotation effective mass at the tip [kg].
    pub effective_mass_kg: f64,
    /// Hertz stiffness K [N/m^1.5].
    pub hertz_stiffness: f64,
}

/// A complete built instrument. Immutable until the next rebuild.
pub struct Instrument {
    /// The admitted declaration this instrument was built from.
    pub params: DrumParams,
    /// Batter then resonant head.
    pub heads: [Head; 2],
    /// Shared mesh node coordinates [m] (both heads use the same geometric mesh).
    pub nodes: Vec<[f64; 2]>,
    /// Shared mesh triangles (counter-clockwise).
    pub tris: Vec<[u32; 3]>,
    /// Lumped P1 nodal quadrature area [m^2] (one third of adjacent facets).
    pub nodal_area_m2: Vec<f64>,
    /// Stick and tip.
    pub stick: StickCard,
    /// Small-signal natural frequencies [Hz] of both heads joined by the
    /// sealed-air spring (dense symmetric eigen-solve of diag(w^2)+(K/V) a a^T).
    pub coupled_frequencies_hz: Vec<f64>,
}

fn mesh_budget() -> ProfileBudget {
    ProfileBudget { max_nodes: 10_000, max_triangles: 20_000, max_feature_evaluations: 100_000 }
}

/// The reference's estimated Z5A-shaped stick (main.rs `stick_with_speed`).
/// Only total length / maximum shaft diameter are published data; stations and
/// density are editable estimates, not a manufacturer drawing.
fn stick_card() -> Result<StickCard, Refusal> {
    let profile = [
        (0.0, 0.0068), (0.015, 0.007_112), (0.280, 0.007_112), (0.365, 0.0035),
        (0.380, 0.0025), (0.393, 0.0045), (0.403, 0.0035), (0.4064, 0.0),
    ]
    .map(|(position_m, radius_m)| RadiusStation { position_m, radius_m });
    let mass = StrikerProperties::from_profile(&profile, 800.0, 0.12, 0.400)
        .map_err(|e| Refusal::new("stick_refused", e.to_string()))?;
    Ok(StickCard {
        full_mass_kg: mass.mass_kg,
        effective_mass_kg: mass.contact_effective_mass_kg,
        // Hertz research approximation K = 4 E_eff sqrt(R_tip) / 3 with
        // E_eff = 0.8 GPa and R = 3 mm: estimates, not measured hickory.
        hertz_stiffness: 4.0 / 3.0 * 0.8e9 * 0.003_f64.sqrt(),
    })
}

fn build_head(params: &DrumParams, head: &HeadParams, index: usize) -> Result<Head, Refusal> {
    let name = if index == 0 { "batter" } else { "resonant" };
    // The section owner remains the authority on E, nu, h, rho.
    fs_plate::PlateSection::isotropic(params.young_pa, params.poisson, head.thickness_m, params.density_kg_m3)
        .map_err(|e| Refusal::new("geometry_refused", format!("{name} head section: {e}")))?;
    let film = TensionedDisk::new_with_tension_variation(
        TensionedDiskSpec {
            radius_m: params.radius_m,
            thickness_m: head.thickness_m,
            young_pa: params.young_pa,
            poisson: params.poisson,
            density_kg_m3: params.density_kg_m3,
            tension_n_m: head.tension_n_m,
            radial_intervals: params.radial_intervals,
            azimuths: params.azimuths,
        },
        head.tension_variation,
        mesh_budget(),
    )
    .map_err(|e| Refusal::new("geometry_refused", format!("{name} head: {e}")))?;
    let tau = core::f64::consts::TAU;
    let window = ((tau * params.band_low_hz).powi(2), (tau * params.band_high_hz).powi(2));
    let modes = fs_modal::slice_window(&film.model.k, &film.model.m, window, &SliceOptions::default())
        .map_err(|e| Refusal::new("eigen_refused", format!("{name} head eigen-solve: {e}")))?
        .modes;
    if modes.is_empty() {
        return Err(Refusal::new(
            "empty_band",
            format!("{name} head has no eigenmode inside {}..{} Hz; no fallback modes are inserted",
                params.band_low_hz, params.band_high_hz),
        ));
    }
    let omegas: Vec<f64> = modes.iter().map(|m| m.lambda.sqrt()).collect();
    if omegas.iter().any(|w| !w.is_finite() || *w <= 0.0) {
        return Err(Refusal::new("eigen_refused", format!("{name} head returned a non-positive eigenvalue")));
    }
    let damping_per_s = omegas.iter().map(|w| 2.0 * head.damping_ratio * w).collect();
    let mut areas = Vec::with_capacity(modes.len());
    for mode in &modes {
        areas.push(film.modal_area(&mode.phi)
            .map_err(|e| Refusal::new("geometry_refused", format!("{name} modal area: {e}")))?);
    }
    let count = modes.len();
    let mut node_shapes = vec![0.0; film.mesh.nodes.len() * count];
    for node in 0..film.mesh.nodes.len() {
        if let Some(dof) = film.model.dof_map[3 * node] {
            for (i, mode) in modes.iter().enumerate() {
                node_shapes[node * count + i] = mode.phi[dof];
            }
        }
    }
    let max_residual = modes.iter().map(|m| m.residual).fold(0.0, f64::max);
    Ok(Head { film, modes, omegas, damping_per_s, areas, node_shapes, max_residual })
}

impl Instrument {
    /// Build both heads from geometry and solve their eigen-windows.
    ///
    /// # Errors
    /// Parameter admission, FE assembly, eigen-solve, empty window or mode budget.
    pub fn build(params: DrumParams) -> Result<Self, Refusal> {
        params.validate()?;
        // Refuse a hopeless window before the eigen-solve spends seconds on
        // modes it could never be allowed to keep (a refused 80-mode solve
        // blocked the browser worker for about a second). Over nine windows
        // (9 to 61 retained modes) the eigen-solve count was 0.93 to 1.11
        // times this estimate, so a 12 % margin keeps every window the
        // budget admits; near-budget windows still go to the eigen-solve.
        let estimate = params.estimated_head_modes(&params.batter) + params.estimated_head_modes(&params.resonant);
        if estimate > 1.12 * MAX_HEAD_MODES as f64 {
            return Err(Refusal::new(
                "mode_budget",
                format!("the declared window holds roughly {estimate:.0} head modes (tensioned-plate Weyl estimate); the live engine admits at most {MAX_HEAD_MODES}. Narrow the window."),
            ));
        }
        let batter = build_head(&params, &params.batter, 0)?;
        let resonant = build_head(&params, &params.resonant, 1)?;
        let total = batter.mode_count() + resonant.mode_count();
        if total > MAX_HEAD_MODES {
            return Err(Refusal::new(
                "mode_budget",
                format!("the two heads retain {total} modes; the live engine admits at most {MAX_HEAD_MODES}. Narrow the declared window, never the returned mode set"),
            ));
        }
        let nodes: Vec<[f64; 2]> = batter.film.mesh.nodes.iter().map(|&(x, y)| [x, y]).collect();
        let tris: Vec<[u32; 3]> = batter.film.mesh.tris.iter().map(|t| t.map(|i| i as u32)).collect();
        if resonant.film.mesh.nodes.len() != nodes.len() || resonant.film.mesh.tris.len() != tris.len() {
            return Err(Refusal::new("geometry_refused", "the two heads did not produce the same geometric mesh"));
        }
        let mut nodal_area_m2 = vec![0.0; nodes.len()];
        for tri in &batter.film.mesh.tris {
            let [a, b, c] = tri.map(|i| batter.film.mesh.nodes[i]);
            let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1));
            if !area.is_finite() || area <= 0.0 {
                return Err(Refusal::new("geometry_refused", "head mesh has a non-positive facet"));
            }
            for &node in tri {
                nodal_area_m2[node] += area / 3.0;
            }
        }
        let stick = stick_card()?;
        let heads = [batter, resonant];
        let coupled_frequencies_hz = coupled_frequencies(&params, &heads);
        Ok(Self { params, heads, nodes, tris, nodal_area_m2, stick, coupled_frequencies_hz })
    }

    /// Adopt parameters that describe the SAME heads, mesh and window (see
    /// `DrumParams::same_structure`): depth, gas, clock and observation only.
    /// The eigen-solve is reused; the coupled diagnostic is recomputed.
    pub fn retarget(&mut self, params: DrumParams) {
        debug_assert!(params.same_structure(&self.params));
        self.params = params;
        self.coupled_frequencies_hz = coupled_frequencies(&self.params, &self.heads);
    }

    /// Batter + resonant retained coordinates.
    #[must_use]
    pub fn head_mode_count(&self) -> usize {
        self.heads[0].mode_count() + self.heads[1].mode_count()
    }

    /// Signed compression areas over the head coordinates only (batter then
    /// resonant): + batter, - resonant, both coordinates positive downward.
    #[must_use]
    pub fn compression_areas(&self) -> Vec<f64> {
        let mut out = self.heads[0].areas.clone();
        out.extend(self.heads[1].areas.iter().map(|a| -a));
        out
    }

    /// Batter-head modal participation at an explicit XY station [1/sqrt(kg)].
    /// Barycentric interpolation on the actual mesh: no snapping, no extrapolation.
    ///
    /// # Errors
    /// `strike_outside_head` when the point is off the polygonal head or has
    /// no retained participation (e.g. exactly on the fixed rim).
    pub fn batter_shapes_at(&self, position_m: [f64; 2]) -> Result<Vec<f64>, Refusal> {
        if position_m.iter().any(|v| !v.is_finite()) {
            return Err(Refusal::new("strike_outside_head", "strike position must be finite metres"));
        }
        let head = &self.heads[0];
        let row = fs_couple::render::plate::impact::linear::wire::film_shapes(&head.film, &head.modes, &[position_m])
            .map_err(|_| Refusal::new(
                "strike_outside_head",
                format!("strike position ({}, {}) m lies outside the meshed batter head (clear radius {} m, {}-gon)",
                    position_m[0], position_m[1], self.params.radius_m, self.params.azimuths),
            ))?
            .remove(0);
        if row.iter().all(|b| *b == 0.0) {
            return Err(Refusal::new(
                "strike_outside_head",
                "strike station has no retained head participation; rim strikes are not modelled",
            ));
        }
        Ok(row)
    }
}

/// Diagnostic only: what the two heads plus the sealed-air spring ring at for
/// vanishing amplitude and no stick. The time stepper never reads these.
fn coupled_frequencies(params: &DrumParams, heads: &[Head; 2]) -> Vec<f64> {
    let n = heads[0].mode_count() + heads[1].mode_count();
    let mut areas = heads[0].areas.clone();
    areas.extend(heads[1].areas.iter().map(|a| -a));
    let omegas: Vec<f64> = heads[0].omegas.iter().chain(&heads[1].omegas).copied().collect();
    let k = params.bulk_modulus_pa() / params.volume_m3();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            matrix[i * n + j] = k * areas[i] * areas[j] + if i == j { omegas[i] * omegas[i] } else { 0.0 };
        }
    }
    let (values, _) = fs_la::eigen::jacobi_eigh(&matrix, n);
    let mut out: Vec<f64> = values.iter().map(|l| l.max(0.0).sqrt() / core::f64::consts::TAU).collect();
    out.sort_by(f64::total_cmp);
    out
}
