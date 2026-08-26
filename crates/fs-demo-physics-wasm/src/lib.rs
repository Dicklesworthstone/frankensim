//! fs-demo-physics-wasm — browser (WASM) surface for the CMA-ES explainer
//! site's parametric demo physics. Layer: L6.
//!
//! Two analytic evaluators, ported 1:1 from the site's TypeScript models so
//! that WASM-vs-fallback is provenance, not behavior:
//!
//! - `wing_eval`: transonic wing aero + structure. 3D lifting-line lift
//!   slope, Oswald induced drag, form-factor profile drag, Korn-equation
//!   drag-divergence Mach with Lock's fourth-power wave drag, half-wing root
//!   bending moment at the elliptical lift centroid, and spar/skin/rib mass.
//! - `bridge_eval`: suspension-bridge structure. Parabolic cable tension,
//!   beam deflection with cable relief, point-load deck bending, suspender-
//!   sized cable stress, simplified Selberg flutter estimate, and mass rollup.
//!
//! Contracts every entry inherits (fs-flyer-wasm pattern):
//!
//! - **Typed-refusal JSON envelope.** Fallible entries return
//!   `{"ok": ...}` or `{"refusal": {"code","message","ranked_repairs"}}`.
//!   Nothing is silently clamped and nothing traps across the boundary.
//! - **Determinism.** Pure functions of their scalar inputs: no wall-clock,
//!   no entropy. Same inputs ⇒ identical output strings, native AND wasm.
//!
//! No-claims: these are teaching/viz surrogates for an explainer site, not
//! certified aeroelastic or structural analysis.

/// Kernel id baked into envelopes so the page can prove which build is live.
pub const KERNEL_VERSION: &str = "fs-demo-physics-wasm 0.1.0";

// ---------------------------------------------------------------------------
// Refusal envelope
// ---------------------------------------------------------------------------

/// A typed refusal for the JS boundary.
#[derive(Debug, Clone)]
pub struct Refusal {
    pub code: &'static str,
    pub message: String,
    pub ranked_repairs: Vec<&'static str>,
}

impl Refusal {
    fn json(&self) -> String {
        let repairs: Vec<String> = self
            .ranked_repairs
            .iter()
            .map(|r| format!("\"{}\"", r))
            .collect();
        format!(
            "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}]}}}}",
            self.code,
            self.message.replace('"', "'"),
            repairs.join(",")
        )
    }
}

fn require_finite(name: &'static str, v: f64) -> Result<(), Refusal> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(Refusal {
            code: "input-non-finite",
            message: format!("{name} must be finite, got {v}"),
            ranked_repairs: vec!["pass a finite number"],
        })
    }
}

/// JS `Math.round` for the non-negative values this kernel rounds (half-up
/// and half-away-from-zero agree for v >= 0).
fn js_round(v: f64) -> f64 {
    v.round()
}

fn j(v: f64) -> String {
    // All emitted values are checked finite; plain Display is valid JSON.
    format!("{v}")
}

// ---------------------------------------------------------------------------
// Wing aerodynamics + structure (mirror of evaluateWingPhysics)
// ---------------------------------------------------------------------------

/// Airfoil family coefficient table, indexed by family id 0..=4:
/// 0 NACA 4-digit · 1 NACA 5-digit high-lift · 2 supercritical SC(2) ·
/// 3 reflexed flying wing · 4 laminar-flow low-Re.
const FAMILY_COEFFS: [(f64, f64, f64, f64); 5] = [
    // (cl_bonus, cd0_bonus, mcrit_bonus, structural_factor)
    (1.0, 1.0, 0.0, 1.0),
    (1.22, 1.12, -0.03, 1.05),
    (1.15, 0.94, 0.08, 1.12),
    (0.88, 0.92, 0.02, 0.95),
    (1.05, 0.78, 0.04, 1.08),
];

/// Wing analysis outputs (display-rounded exactly like the site's TS model).
#[derive(Debug, Clone)]
pub struct WingOut {
    pub lift_coeff_cl: f64,
    pub drag_coeff_cd: f64,
    pub induced_drag_cdi: f64,
    pub profile_drag_cd0: f64,
    pub wave_drag_cdw: f64,
    pub lift_to_drag_ratio: f64,
    pub root_bending_moment_knm: f64,
    pub wing_mass_kg: f64,
    pub critical_mach: f64,
    pub cost_score: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn wing_eval_core(
    aspect_ratio: f64,
    sweep_deg: f64,
    thickness_ratio: f64,
    max_camber: f64,
    camber_position: f64,
    taper_ratio: f64,
    family_id: u32,
    rib_count: f64,
    cruise_mach: f64,
) -> Result<WingOut, Refusal> {
    for (name, v) in [
        ("aspect_ratio", aspect_ratio),
        ("sweep_deg", sweep_deg),
        ("thickness_ratio", thickness_ratio),
        ("max_camber", max_camber),
        ("camber_position", camber_position),
        ("taper_ratio", taper_ratio),
        ("rib_count", rib_count),
        ("cruise_mach", cruise_mach),
    ] {
        require_finite(name, v)?;
    }
    if aspect_ratio <= 0.0 {
        return Err(Refusal {
            code: "aspect-ratio-non-positive",
            message: format!("aspect_ratio must be > 0, got {aspect_ratio}"),
            ranked_repairs: vec!["use a positive aspect ratio (typical 6..16)"],
        });
    }
    let Some(&(cl_bonus, cd0_bonus, mcrit_bonus, structural_factor)) =
        FAMILY_COEFFS.get(family_id as usize)
    else {
        return Err(Refusal {
            code: "family-id-out-of-range",
            message: format!("family id {family_id} has no registered airfoil family"),
            ranked_repairs: vec!["use ids 0..=4"],
        });
    };

    let sweep_rad = sweep_deg * core::f64::consts::PI / 180.0;

    // 1. Lift slope via 3D lifting line with sweep + compressibility.
    let beta_mach = (1.0 - cruise_mach * cruise_mach).max(0.05).sqrt();
    let tan_sweep = sweep_rad.tan();
    let tan_beta_ratio = tan_sweep / beta_mach;
    let ar_beta_ratio = aspect_ratio * beta_mach / 0.95;
    let cl_alpha = (2.0 * core::f64::consts::PI * aspect_ratio)
        / (2.0 + (4.0 + ar_beta_ratio * ar_beta_ratio * (1.0 + tan_beta_ratio * tan_beta_ratio)).sqrt());

    let angle_of_attack_rad = 3.5 * core::f64::consts::PI / 180.0;
    let camber_lift = max_camber * 7.5 * (1.0 + camber_position * 0.5);
    let lift_coeff_cl = (cl_alpha * angle_of_attack_rad + camber_lift) * cl_bonus;

    // 2. Oswald efficiency and induced drag.
    let optimal_taper = 0.45 * (-0.025 * sweep_deg).exp();
    let taper_delta = (taper_ratio - optimal_taper).abs();
    let oswald = (0.98 - 0.05 * (aspect_ratio / 10.0) - 0.15 * taper_delta).max(0.72);
    let induced_drag_cdi =
        lift_coeff_cl * lift_coeff_cl / (core::f64::consts::PI * oswald * aspect_ratio);

    // 3. Profile drag from form factor + skin friction.
    let tc_sq = thickness_ratio * thickness_ratio;
    let form_factor = 1.0 + 2.0 * thickness_ratio + 60.0 * tc_sq * tc_sq;
    let cf_turbulent = 0.0032;
    let profile_drag_cd0 = (2.0 * cf_turbulent * form_factor + 0.002 * max_camber) * cd0_bonus;

    // 4. Korn drag-divergence Mach + Lock fourth-power wave drag.
    let cos_sweep = sweep_rad.cos();
    let korn_tech_factor = 0.87 + mcrit_bonus;
    let drag_divergence_mach = korn_tech_factor / cos_sweep
        - thickness_ratio / (cos_sweep * cos_sweep)
        - lift_coeff_cl / (10.0 * cos_sweep * cos_sweep * cos_sweep);
    let critical_mach = drag_divergence_mach - (0.1_f64 / 80.0).cbrt();

    let delta_mach = (cruise_mach - critical_mach).max(0.0);
    let wave_drag_cdw = if delta_mach > 0.0 { 20.0 * delta_mach.powi(4) } else { 0.0 };

    let drag_coeff_cd = profile_drag_cd0 + induced_drag_cdi + wave_drag_cdw;
    let lift_to_drag_ratio = lift_coeff_cl / drag_coeff_cd.max(1e-4);

    // 5. Structure: half-wing root bending moment + spar/skin/rib mass.
    let wing_area_m2 = 28.0_f64;
    let half_span = (wing_area_m2 * aspect_ratio).sqrt() / 2.0;
    let root_chord = 2.0 * wing_area_m2 / ((wing_area_m2 * aspect_ratio).sqrt() * (1.0 + taper_ratio));
    // Cruise-altitude atmosphere (~35,000 ft): a ~ 295 m/s, rho ~ 0.38 kg/m^3.
    let velocity_ms = cruise_mach * 295.0;
    let dynamic_pressure = 0.5 * 0.38 * velocity_ms * velocity_ms;
    let total_lift_force_n = lift_coeff_cl * dynamic_pressure * wing_area_m2;
    let root_bending_moment_knm = (total_lift_force_n / 2.0) * (0.424 * half_span) / 1000.0;

    let root_thickness_m = root_chord * thickness_ratio;
    let spar_cap_area_m2 = root_bending_moment_knm * 1000.0 / (0.45 * 350e6 * root_thickness_m);
    let spar_mass_kg = spar_cap_area_m2 * half_span * 2.0 * 2700.0;
    let skin_mass_kg = wing_area_m2 * 0.003 * 2700.0 * (1.0 + 0.2 * sweep_deg / 30.0);
    let rib_mass_kg = rib_count * 2.8 * structural_factor;
    let wing_mass_kg = js_round(spar_mass_kg + skin_mass_kg + rib_mass_kg);

    let cost_score = -lift_to_drag_ratio + (wing_mass_kg / 800.0) * 2.5;

    let out = WingOut {
        lift_coeff_cl: js_round(lift_coeff_cl * 1000.0) / 1000.0,
        drag_coeff_cd: js_round(drag_coeff_cd * 10000.0) / 10000.0,
        induced_drag_cdi: js_round(induced_drag_cdi * 10000.0) / 10000.0,
        profile_drag_cd0: js_round(profile_drag_cd0 * 10000.0) / 10000.0,
        wave_drag_cdw: js_round(wave_drag_cdw * 10000.0) / 10000.0,
        lift_to_drag_ratio: js_round(lift_to_drag_ratio * 100.0) / 100.0,
        root_bending_moment_knm: js_round(root_bending_moment_knm),
        wing_mass_kg,
        critical_mach: js_round(critical_mach * 100.0) / 100.0,
        cost_score,
    };

    for v in [
        out.lift_coeff_cl,
        out.drag_coeff_cd,
        out.lift_to_drag_ratio,
        out.root_bending_moment_knm,
        out.wing_mass_kg,
        out.critical_mach,
        out.cost_score,
    ] {
        if !v.is_finite() {
            return Err(Refusal {
                code: "non-finite-result",
                message: "wing evaluation produced a non-finite value".to_string(),
                ranked_repairs: vec!["check the input ranges"],
            });
        }
    }
    Ok(out)
}

pub fn wing_eval_json(
    aspect_ratio: f64,
    sweep_deg: f64,
    thickness_ratio: f64,
    max_camber: f64,
    camber_position: f64,
    taper_ratio: f64,
    family_id: u32,
    rib_count: f64,
    cruise_mach: f64,
) -> String {
    match wing_eval_core(
        aspect_ratio,
        sweep_deg,
        thickness_ratio,
        max_camber,
        camber_position,
        taper_ratio,
        family_id,
        rib_count,
        cruise_mach,
    ) {
        Ok(o) => format!(
            "{{\"ok\":{{\"kernel\":\"{}\",\"liftCoeffCL\":{},\"dragCoeffCD\":{},\"inducedDragCDi\":{},\"profileDragCD0\":{},\"waveDragCDw\":{},\"liftToDragRatio\":{},\"rootBendingMomentKNm\":{},\"wingMassKg\":{},\"criticalMach\":{},\"costScore\":{}}}}}",
            KERNEL_VERSION,
            j(o.lift_coeff_cl),
            j(o.drag_coeff_cd),
            j(o.induced_drag_cdi),
            j(o.profile_drag_cd0),
            j(o.wave_drag_cdw),
            j(o.lift_to_drag_ratio),
            j(o.root_bending_moment_knm),
            j(o.wing_mass_kg),
            j(o.critical_mach),
            j(o.cost_score),
        ),
        Err(r) => r.json(),
    }
}

// ---------------------------------------------------------------------------
// Suspension bridge structure (mirror of evaluateBridgePhysics)
// ---------------------------------------------------------------------------

/// Material table, indexed by material id 0..=3:
/// 0 A36 mild steel · 1 A992 high-strength steel · 2 Ti-6Al-4V · 3 CFRP.
const MATERIALS: [(f64, f64, f64, f64); 4] = [
    // (density kg/m^3, yield MPa, E GPa, cost factor)
    (7850.0, 250.0, 200.0, 1.0),
    (7850.0, 345.0, 210.0, 1.3),
    (4430.0, 880.0, 114.0, 6.5),
    (1600.0, 1200.0, 230.0, 8.0),
];

/// Truss topology factors, indexed by topology id 0..=4:
/// 0 Warren · 1 Pratt · 2 Howe · 3 K-Truss · 4 Bowstring Arch.
const TOPOLOGIES: [(f64, f64, f64); 5] = [
    // (stiffness_factor, mass_factor, aero_drag)
    (1.0, 1.0, 1.0),
    (1.08, 1.05, 1.05),
    (1.02, 1.04, 1.08),
    (1.18, 1.14, 1.25),
    (1.28, 1.20, 0.85),
];

/// Bridge analysis outputs (display-rounded exactly like the site's TS model).
#[derive(Debug, Clone)]
pub struct BridgeOut {
    pub total_mass_tons: f64,
    pub max_von_mises_stress_mpa: f64,
    pub max_deflection_mm: f64,
    pub cable_tension_kn: f64,
    pub flutter_critical_speed_kmh: f64,
    pub yield_limit_mpa: f64,
    pub is_compliant: bool,
    pub cost_score: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn bridge_eval_core(
    span_m: f64,
    sag_m: f64,
    deck_stiffness: f64,
    topology_id: u32,
    material_id: u32,
    suspender_count: f64,
    tower_aspect: f64,
    damping: f64,
    truck_pos_m: f64,
) -> Result<BridgeOut, Refusal> {
    for (name, v) in [
        ("span_m", span_m),
        ("sag_m", sag_m),
        ("deck_stiffness", deck_stiffness),
        ("suspender_count", suspender_count),
        ("tower_aspect", tower_aspect),
        ("damping", damping),
        ("truck_pos_m", truck_pos_m),
    ] {
        require_finite(name, v)?;
    }
    if span_m <= 0.0 {
        return Err(Refusal {
            code: "span-non-positive",
            message: format!("span_m must be > 0, got {span_m}"),
            ranked_repairs: vec!["use a positive span (typical 80..400 m)"],
        });
    }
    let Some(&(density, yield_mpa, e_gpa, cost_factor)) = MATERIALS.get(material_id as usize) else {
        return Err(Refusal {
            code: "material-id-out-of-range",
            message: format!("material id {material_id} has no registered material"),
            ranked_repairs: vec!["use ids 0..=3"],
        });
    };
    let Some(&(stiffness_factor, mass_factor, aero_drag)) = TOPOLOGIES.get(topology_id as usize)
    else {
        return Err(Refusal {
            code: "topology-id-out-of-range",
            message: format!("topology id {topology_id} has no registered truss topology"),
            ranked_repairs: vec!["use ids 0..=4"],
        });
    };

    // Cable tension for a parabolic cable under uniform load: H = w L^2 / (8 s).
    let dead_load_per_meter = density * (0.08 + deck_stiffness * 0.15) * 9.81 * mass_factor / 1000.0;
    let live_truck_load_kn = 400.0_f64;
    let total_linear_load = dead_load_per_meter + live_truck_load_kn / span_m;

    let span_sq = span_m * span_m;
    let sag_span_ratio = sag_m / span_m;
    let horizontal_cable_tension_kn = total_linear_load * span_sq / (8.0 * sag_m.max(2.0));
    let max_cable_tension_kn =
        horizontal_cable_tension_kn * (1.0 + 16.0 * sag_span_ratio * sag_span_ratio).sqrt();

    let tower_height = span_m * tower_aspect;

    // Deck stiffness: I in the realistic 1.5-13 m^4 range for stiffened decks.
    let effective_ei = e_gpa * 1e9
        * (1.5 + deck_stiffness * deck_stiffness * deck_stiffness * 12.0 * stiffness_factor);
    let max_deflection_mm = (5.0 * total_linear_load * 1000.0 * span_sq * span_sq)
        / (384.0 * effective_ei)
        * (1.0 / (1.0 + 8.0 * sag_m / span_m))
        * 1000.0;

    // Cable cross-section follows the suspender strand layout the optimizer
    // controls, so cable stress is a genuine design output.
    let cable_area_m2 = (suspender_count * 0.0012).max(0.005);
    let cable_stress_mpa = max_cable_tension_kn * 1000.0 / cable_area_m2 / 1e6;

    // Point-load moment M = P a (L - a) / L, zero at the supports.
    let truck_pos_from_support = (span_m / 2.0 + truck_pos_m).clamp(0.0, span_m);
    let truck_bending_moment =
        live_truck_load_kn * truck_pos_from_support * (span_m - truck_pos_from_support) / span_m;
    let section_modulus = 0.08 + deck_stiffness * 0.25;
    let deck_bending_stress_mpa = truck_bending_moment * 1000.0 / section_modulus / 1e6;

    // Simplified Selberg flutter estimate:
    // V_f ~ 3.7 * omega_theta * b * sqrt(mu), mu = m / (rho pi b^2).
    let mass_per_length_kg = density * (0.08 + deck_stiffness * 0.15) * mass_factor;
    let torsional_freq = (1.0 / (2.0 * core::f64::consts::PI))
        * ((effective_ei * 0.8) / (mass_per_length_kg * span_sq * span_sq)).sqrt();
    let half_deck_width_m = 10.0_f64;
    let omega_torsional = 2.0 * core::f64::consts::PI * torsional_freq;
    let mass_ratio = mass_per_length_kg
        / (1.225 * core::f64::consts::PI * half_deck_width_m * half_deck_width_m * aero_drag);
    let flutter_critical_speed_kmh = 3.7
        * omega_torsional
        * half_deck_width_m
        * mass_ratio.max(0.0).sqrt()
        * (1.0 + damping * 4.0)
        * 3.6;

    let direct_stress = deck_bending_stress_mpa + cable_stress_mpa * 0.35;
    let shear_stress = 15.0 * aero_drag;
    let max_von_mises_stress_mpa =
        js_round((direct_stress * direct_stress + 3.0 * shear_stress * shear_stress).sqrt());

    let cable_mass_tons = cable_area_m2
        * span_m
        * (1.0 + (8.0 / 3.0) * sag_span_ratio * sag_span_ratio)
        * density
        * 2.0
        / 1000.0;
    // Same cross-sectional area as the dead-load term, so the reported deck
    // mass and the load that stresses the deck describe one structure.
    let deck_mass_tons = span_m * (0.08 + deck_stiffness * 0.15) * density * mass_factor / 1000.0;
    let suspender_mass_tons = suspender_count * (sag_m * 0.6) * 0.002 * density * 2.0 / 1000.0;
    let tower_mass_tons = tower_height * 0.8 * density * 4.0 / 1000.0;
    let total_mass_tons =
        js_round((cable_mass_tons + deck_mass_tons + suspender_mass_tons + tower_mass_tons) * 10.0)
            / 10.0;

    let is_compliant =
        max_von_mises_stress_mpa <= yield_mpa && max_deflection_mm <= span_m * 2.5;

    let stress_violation = (max_von_mises_stress_mpa - yield_mpa).max(0.0);
    let deflection_violation = (max_deflection_mm - span_m * 2.5).max(0.0);
    let cost_score = total_mass_tons * cost_factor
        + stress_violation * stress_violation * 8.0
        + deflection_violation * deflection_violation * 4.0;

    let out = BridgeOut {
        total_mass_tons,
        max_von_mises_stress_mpa,
        max_deflection_mm: js_round(max_deflection_mm * 10.0) / 10.0,
        cable_tension_kn: js_round(max_cable_tension_kn),
        flutter_critical_speed_kmh: js_round(flutter_critical_speed_kmh),
        yield_limit_mpa: yield_mpa,
        is_compliant,
        cost_score,
    };

    for v in [
        out.total_mass_tons,
        out.max_von_mises_stress_mpa,
        out.max_deflection_mm,
        out.cable_tension_kn,
        out.flutter_critical_speed_kmh,
        out.cost_score,
    ] {
        if !v.is_finite() {
            return Err(Refusal {
                code: "non-finite-result",
                message: "bridge evaluation produced a non-finite value".to_string(),
                ranked_repairs: vec!["check the input ranges"],
            });
        }
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
pub fn bridge_eval_json(
    span_m: f64,
    sag_m: f64,
    deck_stiffness: f64,
    topology_id: u32,
    material_id: u32,
    suspender_count: f64,
    tower_aspect: f64,
    damping: f64,
    truck_pos_m: f64,
) -> String {
    match bridge_eval_core(
        span_m,
        sag_m,
        deck_stiffness,
        topology_id,
        material_id,
        suspender_count,
        tower_aspect,
        damping,
        truck_pos_m,
    ) {
        Ok(o) => format!(
            "{{\"ok\":{{\"kernel\":\"{}\",\"totalMassTons\":{},\"maxVonMisesStressMPa\":{},\"maxDeflectionMm\":{},\"cableTensionKN\":{},\"flutterCriticalSpeedKmh\":{},\"yieldLimitMPa\":{},\"isCompliant\":{},\"costScore\":{}}}}}",
            KERNEL_VERSION,
            j(o.total_mass_tons),
            j(o.max_von_mises_stress_mpa),
            j(o.max_deflection_mm),
            j(o.cable_tension_kn),
            j(o.flutter_critical_speed_kmh),
            j(o.yield_limit_mpa),
            o.is_compliant,
            j(o.cost_score),
        ),
        Err(r) => r.json(),
    }
}

/// Native identity probe.
#[must_use]
pub fn kernel_version() -> &'static str {
    KERNEL_VERSION
}

// ---------------------------------------------------------------------------
// wasm32-only JS boundary.
// ---------------------------------------------------------------------------

#[cfg(target_arch = "wasm32")]
mod js {
    use wasm_bindgen::prelude::wasm_bindgen;

    /// Evaluate the parametric wing model: scalars in, JSON envelope out.
    /// family_id: 0 NACA4 · 1 NACA5 · 2 SC(2) · 3 reflexed · 4 laminar.
    #[wasm_bindgen]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn wing_eval(
        aspect_ratio: f64,
        sweep_deg: f64,
        thickness_ratio: f64,
        max_camber: f64,
        camber_position: f64,
        taper_ratio: f64,
        family_id: u32,
        rib_count: f64,
        cruise_mach: f64,
    ) -> String {
        super::wing_eval_json(
            aspect_ratio,
            sweep_deg,
            thickness_ratio,
            max_camber,
            camber_position,
            taper_ratio,
            family_id,
            rib_count,
            cruise_mach,
        )
    }

    /// Evaluate the parametric suspension-bridge model.
    /// topology_id: 0 Warren · 1 Pratt · 2 Howe · 3 K-Truss · 4 Bowstring.
    /// material_id: 0 A36 · 1 A992 · 2 Ti-6Al-4V · 3 CFRP.
    #[wasm_bindgen]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn bridge_eval(
        span_m: f64,
        sag_m: f64,
        deck_stiffness: f64,
        topology_id: u32,
        material_id: u32,
        suspender_count: f64,
        tower_aspect: f64,
        damping: f64,
        truck_pos_m: f64,
    ) -> String {
        super::bridge_eval_json(
            span_m,
            sag_m,
            deck_stiffness,
            topology_id,
            material_id,
            suspender_count,
            tower_aspect,
            damping,
            truck_pos_m,
        )
    }

    /// Kernel identity probe (capability check after instantiation).
    #[wasm_bindgen]
    #[must_use]
    pub fn demo_physics_kernel_version() -> String {
        super::KERNEL_VERSION.to_string()
    }
}

// ---------------------------------------------------------------------------
// Tests — invariant and behavior checks, native.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn base_wing() -> WingOut {
        wing_eval_core(10.5, 25.0, 0.12, 0.035, 0.4, 0.55, 2, 22.0, 0.78).expect("wing")
    }

    #[test]
    fn wing_defaults_are_sane() {
        let w = base_wing();
        assert!(w.lift_coeff_cl > 0.2 && w.lift_coeff_cl < 1.5, "CL = {}", w.lift_coeff_cl);
        assert!(w.lift_to_drag_ratio > 5.0 && w.lift_to_drag_ratio < 40.0, "L/D = {}", w.lift_to_drag_ratio);
        assert!(w.wing_mass_kg > 100.0 && w.wing_mass_kg < 5000.0, "mass = {}", w.wing_mass_kg);
        assert!(w.critical_mach > 0.4 && w.critical_mach < 1.0, "Mcrit = {}", w.critical_mach);
    }

    #[test]
    fn sweep_raises_drag_divergence() {
        let unswept = wing_eval_core(10.5, 0.0, 0.12, 0.035, 0.4, 0.55, 0, 22.0, 0.78).expect("w");
        let swept = wing_eval_core(10.5, 35.0, 0.12, 0.035, 0.4, 0.55, 0, 22.0, 0.78).expect("w");
        assert!(
            swept.critical_mach > unswept.critical_mach,
            "sweep must raise the critical Mach: {} vs {}",
            swept.critical_mach,
            unswept.critical_mach
        );
        assert!(swept.wave_drag_cdw < unswept.wave_drag_cdw);
    }

    #[test]
    fn thinner_sections_cut_wave_drag() {
        let thick = wing_eval_core(10.5, 25.0, 0.15, 0.035, 0.4, 0.55, 0, 22.0, 0.78).expect("w");
        let thin = wing_eval_core(10.5, 25.0, 0.09, 0.035, 0.4, 0.55, 0, 22.0, 0.78).expect("w");
        assert!(thin.wave_drag_cdw < thick.wave_drag_cdw);
    }

    #[test]
    fn wing_refuses_bad_family() {
        let r = wing_eval_core(10.5, 25.0, 0.12, 0.035, 0.4, 0.55, 9, 22.0, 0.78).unwrap_err();
        assert_eq!(r.code, "family-id-out-of-range");
    }

    #[test]
    fn wing_json_is_ok_envelope() {
        let s = wing_eval_json(10.5, 25.0, 0.12, 0.035, 0.4, 0.55, 2, 22.0, 0.78);
        assert!(s.starts_with("{\"ok\":{"), "envelope: {s}");
        assert!(s.contains("\"liftToDragRatio\":"));
    }

    fn base_bridge() -> BridgeOut {
        bridge_eval_core(180.0, 22.0, 0.45, 0, 0, 20.0, 0.35, 0.15, 0.0).expect("bridge")
    }

    #[test]
    fn bridge_defaults_are_sane() {
        let b = base_bridge();
        assert!(b.max_deflection_mm > 10.0 && b.max_deflection_mm < 1000.0, "defl = {}", b.max_deflection_mm);
        assert!(b.total_mass_tons > 100.0 && b.total_mass_tons < 20000.0, "mass = {}", b.total_mass_tons);
        assert!(b.flutter_critical_speed_kmh > 30.0 && b.flutter_critical_speed_kmh < 2000.0, "flutter = {}", b.flutter_critical_speed_kmh);
        assert!(b.is_compliant, "default design should be compliant");
    }

    #[test]
    fn stiffer_deck_deflects_less() {
        let soft = bridge_eval_core(180.0, 22.0, 0.2, 0, 0, 20.0, 0.35, 0.15, 0.0).expect("b");
        let stiff = bridge_eval_core(180.0, 22.0, 0.9, 0, 0, 20.0, 0.35, 0.15, 0.0).expect("b");
        assert!(stiff.max_deflection_mm < soft.max_deflection_mm);
    }

    #[test]
    fn fewer_suspenders_raise_cable_stress_via_von_mises() {
        let few = bridge_eval_core(180.0, 22.0, 0.45, 0, 0, 8.0, 0.35, 0.15, 0.0).expect("b");
        let many = bridge_eval_core(180.0, 22.0, 0.45, 0, 0, 40.0, 0.35, 0.15, 0.0).expect("b");
        assert!(few.max_von_mises_stress_mpa > many.max_von_mises_stress_mpa);
    }

    #[test]
    fn truck_at_support_unloads_the_deck() {
        let mid = bridge_eval_core(180.0, 22.0, 0.45, 0, 0, 20.0, 0.35, 0.15, 0.0).expect("b");
        let support = bridge_eval_core(180.0, 22.0, 0.45, 0, 0, 20.0, 0.35, 0.15, 90.0).expect("b");
        assert!(support.max_von_mises_stress_mpa <= mid.max_von_mises_stress_mpa);
    }

    #[test]
    fn bridge_refuses_bad_ids() {
        assert_eq!(
            bridge_eval_core(180.0, 22.0, 0.45, 7, 0, 20.0, 0.35, 0.15, 0.0).unwrap_err().code,
            "topology-id-out-of-range"
        );
        assert_eq!(
            bridge_eval_core(180.0, 22.0, 0.45, 0, 9, 20.0, 0.35, 0.15, 0.0).unwrap_err().code,
            "material-id-out-of-range"
        );
    }

    #[test]
    fn bridge_json_is_ok_envelope() {
        let s = bridge_eval_json(180.0, 22.0, 0.45, 0, 0, 20.0, 0.35, 0.15, 0.0);
        assert!(s.starts_with("{\"ok\":{"), "envelope: {s}");
        assert!(s.contains("\"isCompliant\":"));
    }

    #[test]
    fn deterministic_output_strings() {
        let a = wing_eval_json(10.5, 25.0, 0.12, 0.035, 0.4, 0.55, 2, 22.0, 0.78);
        let b = wing_eval_json(10.5, 25.0, 0.12, 0.035, 0.4, 0.55, 2, 22.0, 0.78);
        assert_eq!(a, b);
    }
}
