//! 3D CAD geometry generated directly from FrankenSim L2 MORPH primitives:
//!
//! - [`fs_rep_frep::AxisymmetricChart`] (`MeridianSegment::Line` and exact circular
//!   `MeridianSegment::Arc` fillets) for the forged CrMoV rotor, labyrinth seal rings,
//!   inner guide-blade carrier, horizontally-split HP/IP outer casing, slate-blue
//!   axial exhaust diffuser hood, cantilevered rear bearing cone, and front bearing housing.
//! - [`AxisymmetricChart::mass_properties`] and [`AxisymmetricChart::surface_area`]
//!   compute exact Green's-theorem rotor mass `m`, center of mass `z_com`, polar inertia
//!   `J_p`, transverse inertia `J_d`, and wetted surface area directly from the retained
//!   line/arc meridian — feeding `RotorShaftComponent` and `fs_mbd::rotordynamics` with
//!   zero hand-typed inertia values.
//! - [`fs_rep_frep::tessellate_axisymmetric`] produces verified triangle meshes with
//!   evaluated Hausdorff sagitta bounds, Euler characteristic, watertightness, and outward
//!   orientation receipts.
//! - [`fs_rep_nurbs::NurbsSurface`] constructs 7 twisted, cambered, tapered 3D airfoil
//!   families (Curtis impulse, HP, IP1, IP2, LP penultimate, giant 56°-twist LP last-stage
//!   free-standing blade, and stator guide vane), evaluated via [`NurbsSurface::partials`]
//!   for exact analytic normals `normalize(S_u × S_v)`.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use fs_exec::Cx;
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricMassProperties, AxisymmetricMesh,
    AxisymmetricTessellationConfig, AxisymmetricTessellationReceipt, MeridianPoint,
    MeridianSegment, TessellationPurpose, tessellate_axisymmetric,
};
use fs_rep_nurbs::{KnotVector, NurbsSurface};

use crate::materials::library;

/// Total axial length of the model turbine turbine assembly [m] (front coupling to exhaust face).
pub const TURBINE_LENGTH_M: f64 = 5.45;
/// Center of front radial/thrust journal bearing along rotor axis z [m].
pub const FRONT_BEARING_Z_M: f64 = 0.50;
/// Center of rear exhaust radial journal bearing along rotor axis z [m].
pub const REAR_BEARING_Z_M: f64 = 4.44;
/// Rotor bearing span [m].
pub const BEARING_SPAN_M: f64 = REAR_BEARING_Z_M - FRONT_BEARING_Z_M;

/// Compact GPU-ready triangle mesh buffer exported to WebGL/Three.js.
#[derive(Debug, Clone)]
pub struct ExportMesh {
    pub name: &'static str,
    /// Interleaved `[x, y, z, ...]` vertex positions in meters (revolved around `z`).
    pub positions: Vec<f32>,
    /// Interleaved `[nx, ny, nz, ...]` unit outward normals.
    pub normals: Vec<f32>,
    /// Interleaved `[u, v, ...]` parametric machining / airfoil surface UVs for procedural PBR maps.
    pub uvs: Vec<f32>,
    /// Triangle index triples.
    pub indices: Vec<u32>,
    /// First index in `indices` belonging to cutaway cross-section end-caps (`indices.len()` if none).
    /// Lets WebGL render machined bright cross-section steel distinct from cast/enamel outer walls!
    pub cap_index_start: usize,
    /// Optional FrankenSim axisymmetric tessellation receipt (for revolved CAD solids).
    pub receipt: Option<AxisymmetricTessellationReceipt>,
}

impl ExportMesh {
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.positions.len() / 3
    }
}

/// Partitioned cutaway shell of a stationary axisymmetric housing, built directly from
/// a validated `AxisymmetricChart` with solid capped cutaway walls at `θ = 0, π/2, π`.
#[derive(Debug, Clone)]
pub struct CutawayShellMeshes {
    /// Bottom 180° half (`θ ∈ [-π, 0]`): always visible below the horizontal split joint.
    pub lower_half: ExportMesh,
    /// Upper-back 90° quadrant (`θ ∈ [π/2, π]`): visible in Half Cutaway and Closed Exterior.
    pub upper_back: ExportMesh,
    /// Upper-front 90° quadrant (`θ ∈ [0, π/2]`): visible only in Closed Exterior.
    pub upper_front: ExportMesh,
    /// Exact full-360° tessellation receipt from `fs_rep_frep::tessellate_axisymmetric`.
    pub full_receipt: AxisymmetricTessellationReceipt,
    /// Exact solid mass properties from `AxisymmetricChart::mass_properties`.
    pub mass_props: AxisymmetricMassProperties,
    /// Exact boundary area from `AxisymmetricChart::surface_area`.
    pub surface_area_m2: f64,
}

/// Specification of one rotor or stator blade row instantiated from a NURBS airfoil family.
#[derive(Debug, Clone, Copy)]
pub struct BladeRowSpec {
    pub stage_index: usize,
    pub family_id: usize,
    pub is_rotor: bool,
    pub z_center_m: f64,
    pub r_hub_m: f64,
    pub r_tip_m: f64,
    pub blade_count: u32,
    pub axial_chord_m: f64,
    pub has_shroud: bool,
    /// Mean pitch-line radius `(r_hub + r_tip) / 2` [m].
    pub r_mean_m: f64,
    /// Annulus flow area `π (r_tip² - r_hub²)` [m²].
    pub annulus_area_m2: f64,
}

/// Exact geometry-derived mass, inertia, and surface-area properties of the complete
/// rotating assembly (forged shaft + all 17 NURBS blade rows + shrouds).
#[derive(Debug, Clone, Copy)]
pub struct RotorGeometryProps {
    /// Forged CrMoV rotor shaft mass from `AxisymmetricChart::mass_properties` [kg].
    pub shaft_mass_kg: f64,
    /// Added mass of all 17 rotor blade rows and shrouds [kg].
    pub blades_mass_kg: f64,
    /// Total rotating assembly mass [kg].
    pub total_mass_kg: f64,
    /// Center of mass along the z axis [m].
    pub com_z_m: f64,
    /// Polar moment of inertia `J_p = I_z` about spin axis [kg·m²].
    pub polar_inertia_kg_m2: f64,
    /// Transverse centroidal moment of inertia `J_d = I_x = I_y` [kg·m²].
    pub transverse_inertia_kg_m2: f64,
    /// Exact wetted/exposed surface area of the rotor forging [m²].
    pub shaft_surface_area_m2: f64,
    /// Front bearing station relative to rotor COM `z_front - z_com` [m].
    pub front_bearing_rel_z_m: f64,
    /// Rear bearing station relative to rotor COM `z_rear - z_com` [m].
    pub rear_bearing_rel_z_m: f64,
    /// Exact axisymmetric receipt for the forged rotor mesh.
    pub rotor_vertex_count: usize,
    pub rotor_triangle_count: usize,
    pub rotor_hausdorff_bound_m: f64,
    pub rotor_euler_characteristic: i32,
    pub rotor_watertight: bool,
}

/// Complete FrankenSim-generated geometry catalog for the steam turbine model.
#[derive(Debug, Clone)]
pub struct TurbineGeometryCatalog {
    pub rotor_props: RotorGeometryProps,
    pub rotor_mesh: ExportMesh,
    pub labyrinth_seals: ExportMesh,
    pub inner_carrier: CutawayShellMeshes,
    pub outer_hp_ip_casing: CutawayShellMeshes,
    pub exhaust_diffuser: CutawayShellMeshes,
    pub rear_bearing_cone: CutawayShellMeshes,
    pub front_bearing_housing: CutawayShellMeshes,
    /// 7 distinct 3D NURBS airfoil families:
    /// 0: Curtis impulse, 1: HP reaction, 2: IP1, 3: IP2, 4: LP mid, 5: LP giant last stage, 6: Stator vane.
    pub nurbs_airfoils: Vec<ExportMesh>,
    /// All 17 rotor rows + 17 stator guide-vane rows (34 rows) along the exact flow path.
    pub blade_rows: Vec<BladeRowSpec>,
    /// Total FrankenSim-generated triangles across all base solids + NURBS airfoil templates.
    pub total_base_triangles: usize,
    /// Total triangles when all instanced blade rows are accounted for in the 3D scene.
    pub total_scene_triangles: usize,
}

// ---------------------------------------------------------------------------
// Axisymmetric Meridian Construction Helpers
// ---------------------------------------------------------------------------

#[inline]
const fn mp(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}

#[inline]
const fn seg_line(a: MeridianPoint, b: MeridianPoint) -> MeridianSegment {
    MeridianSegment::Line { start: a, end: b }
}

/// Construct a circular fillet arc cleanly on a circle centered at `(cr, cz)` with
/// radius `rad` sweeping from `ang0` to `ang1` (in radians in the `(rho, z)` plane).
fn seg_arc(cr: f64, cz: f64, rad: f64, ang0: f64, ang1: f64, clockwise: bool) -> MeridianSegment {
    let start = mp(cr + rad * ang0.cos(), cz + rad * ang0.sin());
    let end = mp(cr + rad * ang1.cos(), cz + rad * ang1.sin());
    MeridianSegment::Arc {
        start,
        end,
        center: mp(cr, cz),
        clockwise,
    }
}

/// Convert a closed polygon of `(rho, z)` vertices in CCW order into `MeridianSegment::Line`s,
/// skipping degenerate zero-length steps.
fn polyline_to_segments(pts: &[(f64, f64)]) -> Vec<MeridianSegment> {
    let mut segs = Vec::with_capacity(pts.len());
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        if (a.0 - b.0).hypot(a.1 - b.1) > 1.0e-9 {
            segs.push(seg_line(mp(a.0, a.1), mp(b.0, b.1)));
        }
    }
    segs
}

/// Build the exact forged CrMoV rotor shaft as a validated
/// [`AxisymmetricChart`] with true circular `MeridianSegment::Arc` fillets at high-stress
/// transitions (thrust collar, twin balance pistons, Curtis wheel, and giant LP last-stage disc).
pub fn build_rotor_chart() -> AxisymmetricChart {
    let mut segs: Vec<MeridianSegment> = Vec::with_capacity(96);
    let mut cur = mp(0.0, 0.0);

    let push_line = |segs: &mut Vec<MeridianSegment>, cur: &mut MeridianPoint, next: MeridianPoint| {
        if (cur.radius - next.radius).hypot(cur.axial - next.axial) > 1.0e-9 {
            segs.push(seg_line(*cur, next));
            *cur = next;
        }
    };

    let push_arc = |segs: &mut Vec<MeridianSegment>,
                    cur: &mut MeridianPoint,
                    cr: f64,
                    cz: f64,
                    rad: f64,
                    ang0: f64,
                    ang1: f64,
                    clockwise: bool| {
        let mut arc = seg_arc(cr, cz, rad, ang0, ang1, clockwise);
        let (start, end) = match &mut arc {
            MeridianSegment::Arc { start, end, .. } => {
                if (cur.radius - start.radius).hypot(cur.axial - start.axial) <= 1.0e-7 {
                    *start = *cur;
                }
                (*start, *end)
            }
            MeridianSegment::Line { start, end } => (*start, *end),
        };
        if (cur.radius - start.radius).hypot(cur.axial - start.axial) > 1.0e-7 {
            segs.push(seg_line(*cur, start));
        }
        segs.push(arc);
        *cur = end;
    };

    // 1. Front coupling spigot & bolted flange (z = 0.00..0.12)
    push_line(&mut segs, &mut cur, mp(0.132, 0.00));
    push_line(&mut segs, &mut cur, mp(0.132, 0.035));
    push_line(&mut segs, &mut cur, mp(0.205, 0.035));
    push_line(&mut segs, &mut cur, mp(0.205, 0.115));
    // Smooth concave transition fillet down to front shaft stub (r = 0.118)
    push_line(&mut segs, &mut cur, mp(0.136, 0.115));
    push_arc(&mut segs, &mut cur, 0.136, 0.133, 0.018, -FRAC_PI_2, -PI, true);

    // 2. Front shaft stub & Kingsbury thrust collar (z = 0.133..0.38)
    push_line(&mut segs, &mut cur, mp(0.118, 0.255));
    // Concave root fillet into thrust collar
    push_arc(&mut segs, &mut cur, 0.136, 0.255, 0.018, PI, FRAC_PI_2, true);
    push_line(&mut segs, &mut cur, mp(0.214, 0.273));
    push_line(&mut segs, &mut cur, mp(0.214, 0.345));
    push_line(&mut segs, &mut cur, mp(0.146, 0.345));
    push_arc(&mut segs, &mut cur, 0.146, 0.363, 0.018, -FRAC_PI_2, -PI, true);

    // 3. Front radial journal bearing land (z = 0.363..0.635, r = 0.128)
    push_line(&mut segs, &mut cur, mp(0.128, 0.635));

    // 4. Front labyrinth gland seal section (z = 0.635..0.895, r = 0.156 with seal land)
    push_line(&mut segs, &mut cur, mp(0.156, 0.650));
    push_line(&mut segs, &mut cur, mp(0.156, 0.880));

    // 5. Twin stepped balance pistons (z = 0.880..1.195 — matches CAD cutaway!)
    // Balance piston #1 (r = 0.258)
    push_arc(&mut segs, &mut cur, 0.176, 0.880, 0.020, PI, FRAC_PI_2, true);
    push_line(&mut segs, &mut cur, mp(0.258, 0.900));
    push_line(&mut segs, &mut cur, mp(0.258, 1.038));
    // Balance piston #2 (r = 0.304)
    push_line(&mut segs, &mut cur, mp(0.304, 1.050));
    push_line(&mut segs, &mut cur, mp(0.304, 1.185));

    // 6. Curtis impulse control wheel disc (z = 1.185..1.395, r_rim = 0.340)
    push_line(&mut segs, &mut cur, mp(0.340, 1.210));
    push_line(&mut segs, &mut cur, mp(0.340, 1.372));
    // Step down into HP reaction drum with circular fillet
    push_line(&mut segs, &mut cur, mp(0.295, 1.378));
    push_arc(&mut segs, &mut cur, 0.295, 1.398, 0.020, -FRAC_PI_2, -PI, true);

    // 7. HP 5-stage reaction drum (z = 1.398..2.025, r = 0.275)
    push_line(&mut segs, &mut cur, mp(0.275, 2.025));

    // 8. Controlled Extraction 1 scroll recess (z = 2.025..2.255, r = 0.244)
    push_line(&mut segs, &mut cur, mp(0.244, 2.045));
    push_line(&mut segs, &mut cur, mp(0.244, 2.235));

    // 9. IP1 4-stage reaction drum (z = 2.255..2.755, gentle conical taper 0.286 -> 0.296)
    push_line(&mut segs, &mut cur, mp(0.286, 2.255));
    push_line(&mut segs, &mut cur, mp(0.296, 2.755));

    // 10. Uncontrolled bleed recess (z = 2.755..2.935, r = 0.258)
    push_line(&mut segs, &mut cur, mp(0.258, 2.775));
    push_line(&mut segs, &mut cur, mp(0.258, 2.915));

    // 11. IP2 + LP flared reaction drum (z = 2.935..3.580, flared 0.298 -> 0.362)
    push_line(&mut segs, &mut cur, mp(0.298, 2.935));
    push_line(&mut segs, &mut cur, mp(0.326, 3.380));
    push_line(&mut segs, &mut cur, mp(0.362, 3.580));

    // 12. Giant LP last-stage sculpted integral forged disc (z = 3.580..3.960, rim r = 0.410)
    push_line(&mut segs, &mut cur, mp(0.375, 3.635));
    push_line(&mut segs, &mut cur, mp(0.410, 3.690));
    push_line(&mut segs, &mut cur, mp(0.410, 3.855));
    // Sculpted rear disc web sloping sharply back down to rear gland with large stress-relief fillet
    push_line(&mut segs, &mut cur, mp(0.355, 3.910));
    push_line(&mut segs, &mut cur, mp(0.215, 3.965));
    push_arc(&mut segs, &mut cur, 0.215, 4.000, 0.035, -FRAC_PI_2, -PI, true);

    // 13. Rear gland seal land (z = 4.000..4.210, r = 0.175)
    push_line(&mut segs, &mut cur, mp(0.180, 4.005));
    push_line(&mut segs, &mut cur, mp(0.175, 4.015));
    push_line(&mut segs, &mut cur, mp(0.175, 4.210));

    // 14. Rear journal bearing (z = 4.210..4.680, r = 0.132) & shaft end stub (4.68..4.78)
    push_line(&mut segs, &mut cur, mp(0.132, 4.225));
    push_line(&mut segs, &mut cur, mp(0.132, 4.680));
    push_line(&mut segs, &mut cur, mp(0.112, 4.695));
    push_line(&mut segs, &mut cur, mp(0.112, 4.780));

    // Close back to revolution axis and down to origin (CCW loop!)
    push_line(&mut segs, &mut cur, mp(0.0, 4.780));
    push_line(&mut segs, &mut cur, mp(0.0, 0.0));

    AxisymmetricChart::try_new(segs).expect("valid forged rotor AxisymmetricChart")
}

/// Build the annular multi-fin Inconel labyrinth seal carrier ring chart (front gland + balance pistons).
pub fn build_balance_labyrinth_chart() -> AxisymmetricChart {
    // Stepped serrated seal bush wrapping the twin balance pistons (z = 0.905..1.180)
    let mut pts = Vec::new();
    // Bottom radial face at z = 0.905
    pts.push((0.264, 0.905));
    pts.push((0.355, 0.905));
    // Outer cylindrical back face
    pts.push((0.355, 1.180));
    // Top radial face
    pts.push((0.310, 1.180));
    // Inner serrated labyrinth teeth returning downward (z decreasing -> keeps material on left = CCW)
    let n_teeth_2 = 6;
    for i in 0..n_teeth_2 {
        let z_hi = 1.174 - i as f64 * 0.019;
        let z_tip = z_hi - 0.005;
        let z_lo = z_hi - 0.012;
        pts.push((0.324, z_hi));
        pts.push((0.308, z_tip));
        pts.push((0.324, z_lo));
    }
    // Step down to balance piston #1 inner radius
    pts.push((0.278, 1.046));
    let n_teeth_1 = 6;
    for i in 0..n_teeth_1 {
        let z_hi = 1.032 - i as f64 * 0.020;
        let z_tip = z_hi - 0.005;
        let z_lo = z_hi - 0.012;
        pts.push((0.278, z_hi));
        pts.push((0.263, z_tip));
        pts.push((0.278, z_lo));
    }
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid labyrinth seal AxisymmetricChart")
}

/// Build the cast CrMoV inner guide-blade carrier (annular stepped stator shell holding
/// HP/IP/LP diaphragms, matching `image-1405.jpg` and the dark inner cutaway shell of `image-0031.jpg`).
pub fn build_inner_carrier_chart() -> AxisymmetricChart {
    // CCW annular profile: outer contour going z: 1.18 -> 3.63, inner contour returning z: 3.63 -> 1.18
    let pts = vec![
        // Front wall z = 1.18
        (0.360, 1.180),
        (0.470, 1.180),
        // Outer stepped carrier wall
        (0.470, 1.440),
        (0.460, 1.460),
        (0.460, 2.010),
        // Extraction scroll recess pocket on outside
        (0.485, 2.040),
        (0.485, 2.240),
        (0.545, 2.265),
        (0.575, 2.750),
        (0.610, 2.940),
        (0.765, 3.380),
        (0.840, 3.630),
        // Rear carrier lip z = 3.63
        (0.735, 3.630),
        // Inner flowpath shroud ceiling returning from LP back to HP (z decreasing!)
        (0.625, 3.400),
        (0.560, 3.340),
        (0.470, 2.960),
        (0.452, 2.740),
        (0.398, 2.270),
        (0.398, 2.020),
        (0.384, 1.980),
        (0.342, 1.430),
        (0.396, 1.380),
        (0.396, 1.210),
    ];
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid inner guide-blade carrier AxisymmetricChart")
}

/// Build the thick horizontally-split cast CrMoV outer HP/IP casing (`z = 0.62..3.68 m`),
/// matching the light-grey cast steel pressure barrel in `image-0011.jpg` and `image-0031.jpg`.
pub fn build_outer_hp_ip_casing_chart() -> AxisymmetricChart {
    let pts = vec![
        // Front gland end wall (z = 0.62)
        (0.172, 0.620),
        (0.340, 0.620),
        (0.380, 0.720),
        // HP nozzle-box bulge (z = 0.85..1.48)
        (0.555, 0.860),
        (0.590, 0.940),
        (0.590, 1.450),
        (0.560, 1.520),
        // HP barrel
        (0.560, 1.980),
        // Extraction scroll external enlargement (z = 2.00..2.38)
        (0.665, 2.050),
        (0.665, 2.340),
        (0.620, 2.400),
        // Conical IP-to-LP flare (matches the iconic tapered light-grey waist in image-0011!)
        (0.715, 2.860),
        (0.890, 3.350),
        (0.995, 3.580),
        // Heavy vertical bolted flange mating the blue LP exhaust hood (z = 3.58..3.68)
        (1.085, 3.580),
        (1.085, 3.680),
        (0.925, 3.680),
        // Inner pressure-containment bore returning z: 3.68 -> 0.62
        (0.855, 3.580),
        (0.780, 3.360),
        (0.625, 2.900),
        (0.565, 2.360),
        (0.565, 2.020),
        (0.490, 1.960),
        (0.490, 1.460),
        (0.495, 1.160),
        (0.375, 1.150),
        (0.375, 0.895),
        (0.210, 0.885),
        (0.210, 0.670),
        (0.172, 0.670),
    ];
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid outer HP/IP casing AxisymmetricChart")
}

/// Build the iconic slate-blue fabricated axial exhaust diffuser hood (`z = 3.66..5.45 m`),
/// with flared diffuser outer wall, bolted front mating flange, circumferential stiffener rings,
/// and circular condenser-mating mouth flange (matching `image-0011.jpg` and `image-0031.jpg`).
pub fn build_exhaust_diffuser_chart() -> AxisymmetricChart {
    let pts = vec![
        // Front vertical mating flange bolted to HP/IP casing (z = 3.66..3.74)
        (0.925, 3.660),
        (1.105, 3.660),
        (1.105, 3.745),
        (1.045, 3.745),
        // Outer fabricated shell with 3 stepped circumferential stiffener flange rings
        (1.090, 4.080),
        (1.160, 4.080),
        (1.160, 4.135),
        (1.098, 4.135),
        (1.155, 4.540),
        (1.235, 4.540),
        (1.235, 4.600),
        (1.165, 4.600),
        (1.225, 4.980),
        (1.295, 4.980),
        (1.295, 5.040),
        (1.232, 5.040),
        // Flared exhaust mouth rim & outer ring flange (z = 5.38..5.45)
        (1.270, 5.380),
        (1.315, 5.380),
        (1.315, 5.450),
        (1.230, 5.450),
        // Smooth flared inner exhaust diffuser steam-path wall returning z: 5.45 -> 3.66
        (1.190, 5.050),
        (1.115, 4.550),
        (1.020, 4.020),
        (0.960, 3.880),
        (0.955, 3.720),
        (0.925, 3.720),
    ];
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid axial exhaust diffuser AxisymmetricChart")
}

/// Build the cantilevered central rear bearing housing tunnel & aerodynamic exhaust cone
/// inside the axial exhaust diffuser (`z = 3.98..5.14 m`, clearly visible in `image-0031.jpg`).
pub fn build_rear_bearing_cone_chart() -> AxisymmetricChart {
    let pts = vec![
        // Front seal bulkhead facing LP last-stage disc (z = 3.98)
        (0.184, 3.980),
        (0.395, 3.980),
        // Outer diffuser hub sleeve tapering smoothly into exhaust cone
        (0.415, 4.160),
        (0.415, 4.780),
        (0.390, 4.920),
        // Tapered end bulkhead
        (0.335, 5.060),
        (0.240, 5.140),
        (0.240, 5.090),
        // Hollow interior oil/bearing access chamber returning z: 5.09 -> 3.98
        (0.320, 4.980),
        (0.330, 4.720),
        (0.245, 4.720),
        (0.245, 4.220),
        (0.205, 4.200),
        (0.205, 4.020),
        (0.184, 4.020),
    ];
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid rear bearing cone AxisymmetricChart")
}

/// Build the front combined radial + Kingsbury thrust bearing pedestal shell (`z = 0.14..0.64 m`).
pub fn build_front_bearing_housing_chart() -> AxisymmetricChart {
    let pts = vec![
        (0.132, 0.145),
        (0.285, 0.145),
        (0.320, 0.195),
        (0.320, 0.595),
        (0.285, 0.635),
        (0.142, 0.635),
        (0.142, 0.590),
        (0.245, 0.590),
        (0.245, 0.375),
        (0.260, 0.360),
        (0.260, 0.255),
        (0.185, 0.245),
        (0.185, 0.190),
        (0.132, 0.190),
    ];
    AxisymmetricChart::try_new(polyline_to_segments(&pts))
        .expect("valid front bearing housing AxisymmetricChart")
}

// ---------------------------------------------------------------------------
// Axisymmetric Mesh Conversion & Solid Cutaway Slicing
// ---------------------------------------------------------------------------

fn axisymmetric_to_export(name: &'static str, mesh: AxisymmetricMesh) -> ExportMesh {
    let mut positions = Vec::with_capacity(mesh.positions.len() * 3);
    let mut normals = Vec::with_capacity(mesh.normals.len() * 3);
    let mut uvs = Vec::with_capacity(mesh.positions.len() * 2);
    for (p, n) in mesh.positions.iter().zip(mesh.normals.iter()) {
        positions.push(p.x as f32);
        positions.push(p.y as f32);
        positions.push(p.z as f32);
        normals.push(n.x as f32);
        normals.push(n.y as f32);
        normals.push(n.z as f32);
        // Cylindrical lathe-turning coordinates: u along axial z, v around circumference
        let th = p.y.atan2(p.x);
        let r = p.x.hypot(p.y);
        uvs.push((p.z * 2.4) as f32);
        uvs.push(((th / TAU + 0.5) * (r * 6.0).max(1.0)) as f32);
    }
    let mut indices = Vec::with_capacity(mesh.triangles.len() * 3);
    for tri in &mesh.triangles {
        indices.push(tri[0]);
        indices.push(tri[1]);
        indices.push(tri[2]);
    }
    let cap_index_start = indices.len();
    ExportMesh {
        name,
        positions,
        normals,
        uvs,
        indices,
        cap_index_start,
        receipt: Some(mesh.receipt),
    }
}

/// Extract the dense ordered 2D CCW meridian polyline `(rho, z)` from a validated
/// `AxisymmetricChart`, sampling any `MeridianSegment::Arc` finely for solid cutaway caps.
fn sample_meridian_polygon(chart: &AxisymmetricChart) -> Vec<(f64, f64)> {
    let mut pts = Vec::new();
    for seg in chart.segments() {
        match *seg {
            MeridianSegment::Line { start, .. } => {
                pts.push((start.radius, start.axial));
            }
            MeridianSegment::Arc {
                start,
                end,
                center,
                clockwise,
            } => {
                let r = (start.radius - center.radius).hypot(start.axial - center.axial);
                let a0 = (start.axial - center.axial).atan2(start.radius - center.radius);
                let a1 = (end.axial - center.axial).atan2(end.radius - center.radius);
                let mut sweep = a1 - a0;
                if clockwise {
                    while sweep >= 0.0 {
                        sweep -= TAU;
                    }
                } else {
                    while sweep <= 0.0 {
                        sweep += TAU;
                    }
                }
                let steps = 8;
                for k in 0..steps {
                    let t = k as f64 / steps as f64;
                    let a = a0 + sweep * t;
                    pts.push((center.radius + r * a.cos(), center.axial + r * a.sin()));
                }
            }
        }
    }
    pts
}

/// Robust 2D ear-clipping triangulation for a simple CCW meridian polygon in `(rho, z)`.
/// Used to generate crisp solid cross-section end-caps on cutaway turbine casings.
fn triangulate_meridian_ccw(poly: &[(f64, f64)]) -> Vec<[usize; 3]> {
    let n = poly.len();
    if n < 3 {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..n).collect();
    let mut tris = Vec::with_capacity(n.saturating_sub(2));

    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| -> f64 {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let point_in_tri = |p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)| -> bool {
        let c1 = cross(a, b, p);
        let c2 = cross(b, c, p);
        let c3 = cross(c, a, p);
        c1 >= -1.0e-11 && c2 >= -1.0e-11 && c3 >= -1.0e-11
    };

    let mut guard = 0usize;
    let max_iter = n * n + 8;
    while idx.len() > 3 && guard < max_iter {
        guard += 1;
        let m = idx.len();
        let mut clipped = false;
        for i in 0..m {
            let i_prev = idx[(i + m - 1) % m];
            let i_curr = idx[i];
            let i_next = idx[(i + 1) % m];
            let (a, b, c) = (poly[i_prev], poly[i_curr], poly[i_next]);
            if cross(a, b, c) <= 1.0e-12 {
                continue;
            }
            let mut has_inside = false;
            for &k in &idx {
                if k == i_prev || k == i_curr || k == i_next {
                    continue;
                }
                if point_in_tri(poly[k], a, b, c) {
                    has_inside = true;
                    break;
                }
            }
            if !has_inside {
                tris.push([i_prev, i_curr, i_next]);
                idx.remove(i);
                clipped = true;
                break;
            }
        }
        if !clipped {
            // Fallback fan for nearly-collinear sliver
            tris.push([idx[0], idx[1], idx[2]]);
            idx.remove(1);
        }
    }
    if idx.len() == 3 {
        tris.push([idx[0], idx[1], idx[2]]);
    }
    tris
}

/// Build an exact revolved angular sector `[theta_start, theta_end]` of an admitted
/// `AxisymmetricChart`, including solid planar cross-section end-caps at `theta_start`
/// and `theta_end` so cutaway steel casings look genuinely thick and solid!
fn build_revolved_sector(
    name: &'static str,
    chart: &AxisymmetricChart,
    theta_start: f64,
    theta_end: f64,
    n_sectors: usize,
) -> ExportMesh {
    let poly = sample_meridian_polygon(chart);
    let m = poly.len();
    let mut positions: Vec<f32> = Vec::new();
    let mut normals: Vec<f32> = Vec::new();
    let mut uvs: Vec<f32> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    // Compute per-edge outward meridian normals (for CCW in (rho, z), right-hand normal is (dz, -dr))
    let mut arc_s = 0.0_f64;
    for seg_idx in 0..m {
        let p0 = poly[seg_idx];
        let p1 = poly[(seg_idx + 1) % m];
        let dr = p1.0 - p0.0;
        let dz = p1.1 - p0.1;
        let len = dr.hypot(dz);
        if len < 1.0e-10 {
            continue;
        }
        let nr = dz / len;
        let nz = -dr / len;
        let s0 = arc_s;
        let s1 = arc_s + len;
        arc_s = s1;

        let base = (positions.len() / 3) as u32;
        for k in 0..=n_sectors {
            let frac = k as f64 / n_sectors as f64;
            let th = theta_start + (theta_end - theta_start) * frac;
            let (ct, st) = (th.cos(), th.sin());
            let nx = (nr * ct) as f32;
            let ny = (nr * st) as f32;
            let nz_f = nz as f32;
            let v_arc = (th * 1.4) as f32;

            // Vertex at p0
            positions.extend_from_slice(&[(p0.0 * ct) as f32, (p0.0 * st) as f32, p0.1 as f32]);
            normals.extend_from_slice(&[nx, ny, nz_f]);
            uvs.extend_from_slice(&[(s0 * 1.6) as f32, v_arc]);
            // Vertex at p1
            positions.extend_from_slice(&[(p1.0 * ct) as f32, (p1.0 * st) as f32, p1.1 as f32]);
            normals.extend_from_slice(&[nx, ny, nz_f]);
            uvs.extend_from_slice(&[(s1 * 1.6) as f32, v_arc]);
        }

        for k in 0..n_sectors as u32 {
            let a = base + 2 * k;
            let b = a + 2 * k + 1 - 2 * k; // = a + 1
            let c = base + 2 * (k + 1);
            let d = base + 2 * (k + 1) + 1;
            // Winding produces outward normal matching (nr*cos(th), nr*sin(th), nz)
            indices.extend_from_slice(&[a, c, b, b, c, d]);
        }
    }

    // Record exact split index so WebGL renders solid cutaway end-caps in bright milled-steel!
    let cap_index_start = indices.len();

    // Solid machined cross-section end-caps at theta_start and theta_end!
    let cap_tris = triangulate_meridian_ccw(&poly);
    let add_cap = |theta: f64,
                   flip: bool,
                   positions: &mut Vec<f32>,
                   normals: &mut Vec<f32>,
                   uvs: &mut Vec<f32>,
                   indices: &mut Vec<u32>| {
        let (ct, st) = (theta.cos(), theta.sin());
        // Tangent +theta direction is (-sin(th), cos(th), 0).
        // End cap at theta_end faces +theta; start cap at theta_start faces -theta.
        let sign = if flip { -1.0 } else { 1.0 };
        let nx = (-st * sign) as f32;
        let ny = (ct * sign) as f32;
        let base = (positions.len() / 3) as u32;
        for &(r, z) in &poly {
            positions.extend_from_slice(&[(r * ct) as f32, (r * st) as f32, z as f32]);
            normals.extend_from_slice(&[nx, ny, 0.0]);
            uvs.extend_from_slice(&[(z * 2.5) as f32, (r * 2.5) as f32]);
        }
        for &[i0, i1, i2] in &cap_tris {
            if flip {
                // In (r, z) plane, e_r × e_z = -e_theta, so CCW in (r, z) naturally faces -e_theta!
                indices.extend_from_slice(&[
                    base + i0 as u32,
                    base + i1 as u32,
                    base + i2 as u32,
                ]);
            } else {
                // Reverse winding so normal faces +e_theta
                indices.extend_from_slice(&[
                    base + i0 as u32,
                    base + i2 as u32,
                    base + i1 as u32,
                ]);
            }
        }
    };
    add_cap(theta_start, true, &mut positions, &mut normals, &mut uvs, &mut indices);
    add_cap(theta_end, false, &mut positions, &mut normals, &mut uvs, &mut indices);

    ExportMesh {
        name,
        positions,
        normals,
        uvs,
        indices,
        cap_index_start,
        receipt: None,
    }
}

/// Tessellate a stationary housing with `fs_rep_frep::tessellate_axisymmetric` to obtain its
/// verified receipt, exact mass properties, and surface area, and build its 3 solid cutaway
/// render shells (`lower_half`, `upper_back`, `upper_front`).
fn build_cutaway_housing(
    base_name: &'static str,
    lower_name: &'static str,
    back_name: &'static str,
    front_name: &'static str,
    chart: &AxisymmetricChart,
    density_kg_m3: f64,
    hausdorff_m: f64,
    sectors_per_quadrant: usize,
    cx: &Cx<'_>,
) -> CutawayShellMeshes {
    let mut cfg = AxisymmetricTessellationConfig::new(hausdorff_m, TessellationPurpose::Rendering)
        .expect("valid tessellation config");
    cfg.min_azimuthal_sectors = 48;
    cfg.max_azimuthal_sectors = 96;
    let full_mesh = tessellate_axisymmetric(chart, cfg, cx)
        .unwrap_or_else(|e| panic!("tessellate_axisymmetric failed for {base_name}: {e}"));
    let mass_props = chart
        .mass_properties(density_kg_m3, cx)
        .unwrap_or_else(|e| panic!("mass_properties failed for {base_name}: {e}"));
    let surface_area_m2 = chart
        .surface_area(cx)
        .unwrap_or_else(|e| panic!("surface_area failed for {base_name}: {e}"))
        .area;

    // Coordinate mapping note:
    // When Three.js maps axisymmetric (x_a, y_a, z_a) -> world (-z_a, y_a, x_a):
    //   y_a = r sin(theta) = world UP (+Y)
    //   x_a = r cos(theta) = world FRONT (+Z, toward camera!)
    // Therefore:
    //   theta in [-pi, 0]     -> y_a <= 0 -> LOWER HALF (always visible!)
    //   theta in [pi/2, pi]   -> y_a >= 0, x_a <= 0 -> UPPER BACK quadrant (visible in Half Cutaway!)
    //   theta in [0, pi/2]    -> y_a >= 0, x_a >= 0 -> UPPER FRONT quadrant (cut away in Half Cutaway!)
    let lower_half = build_revolved_sector(lower_name, chart, -PI, 0.0, sectors_per_quadrant * 2);
    let upper_back = build_revolved_sector(back_name, chart, FRAC_PI_2, PI, sectors_per_quadrant);
    let upper_front = build_revolved_sector(front_name, chart, 0.0, FRAC_PI_2, sectors_per_quadrant);

    CutawayShellMeshes {
        lower_half,
        upper_back,
        upper_front,
        full_receipt: full_mesh.receipt,
        mass_props,
        surface_area_m2,
    }
}

// ---------------------------------------------------------------------------
// 3D Twisted, Cambered Airfoil Blades & Stator Vanes via fs-rep-nurbs
// ---------------------------------------------------------------------------

/// Parametric specification of one 3D NURBS turbomachinery blade/vane family.
///
/// The resulting [`ExportMesh`] is generated in local normalized coordinates:
/// - radial axis `+Y_local` spanning `[0.0, 1.0]` (scaled to `blade_height` and offset by `r_hub` per row),
/// - axial chord along `+Z_local` centered at `0.0` (scaled to `axial_chord_m`),
/// - tangential camber/stagger along `+X_local`.
#[derive(Debug, Clone, Copy)]
struct AirfoilFamilyParams {
    name: &'static str,
    /// Stagger angle at hub root [deg] measured from axial (+Z) toward tangential (+X).
    stagger_root_deg: f64,
    /// Stagger angle at tip [deg] (free-vortex twist increases toward tip).
    stagger_tip_deg: f64,
    /// Maximum camber-to-chord ratio at root (e.g. 0.22 for Curtis impulse, 0.14 for reaction).
    camber_root: f64,
    /// Maximum camber-to-chord ratio at tip.
    camber_tip: f64,
    /// Max thickness-to-chord ratio at root.
    thickness_root: f64,
    /// Max thickness-to-chord ratio at tip.
    thickness_tip: f64,
    /// Tip-to-root chord taper ratio (e.g. 0.50 for giant LP last-stage blade).
    chord_taper: f64,
    /// Tangential lean at tip as fraction of chord.
    lean_tip: f64,
    /// Whether camber bends in opposite direction (stator guide vane).
    camber_sign: f64,
    /// Number of chordwise samples around the closed NURBS contour.
    nu: usize,
    /// Number of spanwise samples along the blade height.
    nv: usize,
}

/// Build a closed clamped uniform knot vector of degree `p` for `n_ctrl` control points.
fn clamped_uniform_knots(n_ctrl: usize, p: usize) -> KnotVector<f64> {
    let n_knots = n_ctrl + p + 1;
    let mut knots = Vec::with_capacity(n_knots);
    let n_internal = n_ctrl - p - 1;
    for _ in 0..=p {
        knots.push(0.0);
    }
    for i in 1..=n_internal {
        knots.push(i as f64 / (n_internal + 1) as f64);
    }
    for _ in 0..=p {
        knots.push(1.0);
    }
    KnotVector::new(knots, p).expect("valid clamped uniform NURBS knot vector")
}

/// Construct one 3D twisted, cambered, tapered turbomachinery airfoil using
/// [`fs_rep_nurbs::NurbsSurface<f64>`], and tessellate it via exact [`NurbsSurface::partials`]
/// (`S_u × S_v` analytic outward normals) with closed root and tip caps.
fn build_nurbs_airfoil_mesh(p: AirfoilFamilyParams) -> ExportMesh {
    // 11 control points around the closed airfoil contour (u in [0, 1], degree 3):
    // TE -> suction side (s = 1.0 -> 0.0) -> LE -> pressure side (s = 0.0 -> 1.0) -> TE
    // Parametric normalized chord fractions along upper (suction) and lower (pressure) surfaces:
    const U_CTRL: [(f64, f64); 11] = [
        (1.00, 0.00),  // 0: TE
        (0.80, 0.60),  // 1: suction aft
        (0.52, 0.96),  // 2: suction mid
        (0.24, 1.00),  // 3: suction fore
        (0.04, 0.62),  // 4: LE upper shoulder
        (0.00, 0.00),  // 5: LE nose
        (0.06, -0.45), // 6: LE lower shoulder
        (0.28, -0.70), // 7: pressure fore
        (0.56, -0.68), // 8: pressure mid
        (0.82, -0.40), // 9: pressure aft
        (1.00, 0.00),  // 10: TE closed
    ];

    // 5 spanwise control stations (v in [0, 1], degree 3)
    const V_STATIONS: [f64; 5] = [0.0, 0.25, 0.50, 0.75, 1.0];

    let mut points: Vec<Vec<[f64; 3]>> = Vec::with_capacity(U_CTRL.len());
    let mut weights: Vec<Vec<f64>> = Vec::with_capacity(U_CTRL.len());

    for &(s_chord, side_t) in &U_CTRL {
        let mut row_pts = Vec::with_capacity(V_STATIONS.len());
        let mut row_w = Vec::with_capacity(V_STATIONS.len());
        // Slightly higher rational weight on the leading edge nose for a crisp rounded LE
        let w_u = if (s_chord - 0.0).abs() < 1.0e-6 { 1.35 } else { 1.0 };

        for &v in &V_STATIONS {
            // Smooth free-vortex twist progression from hub to tip
            let twist_frac = v.powf(0.85);
            let stagger = (p.stagger_root_deg + (p.stagger_tip_deg - p.stagger_root_deg) * twist_frac)
                .to_radians();
            let chord = 1.0 + (p.chord_taper - 1.0) * v;
            let camber = p.camber_root + (p.camber_tip - p.camber_root) * v;
            let thick = p.thickness_root + (p.thickness_tip - p.thickness_root) * v;

            // Parabolic camber line + airfoil thickness distribution in blade-section (xi, eta) coords
            let xi = (s_chord - 0.5) * chord;
            let camber_y = p.camber_sign * 4.0 * camber * s_chord * (1.0 - s_chord) * chord;
            let half_t = 0.5 * thick * chord * side_t;
            let eta = camber_y + half_t;

            // Rotate by stagger angle into (tangential X, radial Y, axial Z)
            let (cs, sn) = (stagger.cos(), stagger.sin());
            let x_tan = xi * sn + eta * cs + p.lean_tip * v * v;
            let z_ax = xi * cs - eta * sn;
            let y_rad = v; // normalized [0, 1] radial span

            row_pts.push([x_tan, y_rad, z_ax]);
            row_w.push(w_u);
        }
        points.push(row_pts);
        weights.push(row_w);
    }

    let ku = clamped_uniform_knots(U_CTRL.len(), 3);
    let kv = clamped_uniform_knots(V_STATIONS.len(), 3);
    let surface = NurbsSurface::new(ku, kv, &points, &weights)
        .unwrap_or_else(|e| panic!("valid NURBS airfoil surface for {}: {e}", p.name));

    let nu = p.nu;
    let nv = p.nv;
    let mut positions: Vec<f32> = Vec::with_capacity((nu + 1) * (nv + 1) * 3 + 2 * (nu + 2) * 3);
    let mut normals: Vec<f32> = Vec::with_capacity(positions.capacity());
    let mut uvs: Vec<f32> = Vec::with_capacity((nu + 1) * (nv + 1) * 2 + 2 * (nu + 2) * 2);
    let mut indices: Vec<u32> = Vec::with_capacity(nu * nv * 6 + 2 * nu * 3);

    // Evaluate surface points and EXACT analytic normals via NurbsSurface::partials(u, v)
    let mut grid_pts: Vec<[f64; 3]> = Vec::with_capacity((nu + 1) * (nv + 1));
    for iv in 0..=nv {
        let v = iv as f64 / nv as f64;
        for iu in 0..=nu {
            let u = iu as f64 / nu as f64;
            let (pos, su, sv) = surface
                .partials(u, v)
                .expect("NURBS airfoil partials inside [0,1]^2");
            // Cross product S_v × S_u points outward on the airfoil skin
            let mut n = [
                sv[1] * su[2] - sv[2] * su[1],
                sv[2] * su[0] - sv[0] * su[2],
                sv[0] * su[1] - sv[1] * su[0],
            ];
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if len > 1.0e-12 {
                n[0] /= len;
                n[1] /= len;
                n[2] /= len;
            } else {
                // Clean tangent fallback at exact trailing-edge cusp
                n = [0.0, 0.0, 1.0];
            }
            grid_pts.push(pos);
            positions.extend_from_slice(&[pos[0] as f32, pos[1] as f32, pos[2] as f32]);
            normals.extend_from_slice(&[n[0] as f32, n[1] as f32, n[2] as f32]);
            uvs.extend_from_slice(&[u as f32, v as f32]);
        }
    }

    let stride = (nu + 1) as u32;
    for iv in 0..nv as u32 {
        for iu in 0..nu as u32 {
            let a = iv * stride + iu;
            let b = a + 1;
            let c = (iv + 1) * stride + iu;
            let d = c + 1;
            indices.extend_from_slice(&[a, c, b, b, c, d]);
        }
    }

    // Close root cap (v = 0, normal -Y) and tip cap (v = 1, normal +Y)
    let add_end_cap = |iv: usize,
                       ny: f32,
                       flip: bool,
                       grid_pts: &[[f64; 3]],
                       positions: &mut Vec<f32>,
                       normals: &mut Vec<f32>,
                       uvs: &mut Vec<f32>,
                       indices: &mut Vec<u32>| {
        let row_start = iv * (nu + 1);
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for iu in 0..nu {
            let pt = grid_pts[row_start + iu];
            cx += pt[0];
            cy += pt[1];
            cz += pt[2];
        }
        let inv = 1.0 / nu as f64;
        let center_idx = (positions.len() / 3) as u32;
        positions.extend_from_slice(&[(cx * inv) as f32, (cy * inv) as f32, (cz * inv) as f32]);
        normals.extend_from_slice(&[0.0, ny, 0.0]);
        uvs.extend_from_slice(&[0.5, 0.5]);
        let ring_start = (positions.len() / 3) as u32;
        for iu in 0..=nu {
            let pt = grid_pts[row_start + iu];
            positions.extend_from_slice(&[pt[0] as f32, pt[1] as f32, pt[2] as f32]);
            normals.extend_from_slice(&[0.0, ny, 0.0]);
            uvs.extend_from_slice(&[(iu as f32) / (nu as f32), 0.5]);
        }
        for iu in 0..nu as u32 {
            if flip {
                indices.extend_from_slice(&[center_idx, ring_start + iu, ring_start + iu + 1]);
            } else {
                indices.extend_from_slice(&[center_idx, ring_start + iu + 1, ring_start + iu]);
            }
        }
    };
    add_end_cap(0, -1.0, false, &grid_pts, &mut positions, &mut normals, &mut uvs, &mut indices);
    add_end_cap(nv, 1.0, true, &grid_pts, &mut positions, &mut normals, &mut uvs, &mut indices);
    let cap_index_start = indices.len();

    ExportMesh {
        name: p.name,
        positions,
        normals,
        uvs,
        indices,
        cap_index_start,
        receipt: None,
    }
}

/// Build the full 17-stage rotor + 16-stage stator flow-path schedule after the
/// condensing/extraction cutaway (`image-0031.jpg` and `image-0009.jpg`).
pub fn build_blade_rows() -> Vec<BladeRowSpec> {
    let mut rows = Vec::with_capacity(33);
    let mut push_row = |stage_index: usize,
                        family_id: usize,
                        is_rotor: bool,
                        z_center_m: f64,
                        r_hub_m: f64,
                        r_tip_m: f64,
                        blade_count: u32,
                        axial_chord_m: f64,
                        has_shroud: bool| {
        rows.push(BladeRowSpec {
            stage_index,
            family_id,
            is_rotor,
            z_center_m,
            r_hub_m,
            r_tip_m,
            blade_count,
            axial_chord_m,
            has_shroud,
            r_mean_m: 0.5 * (r_hub_m + r_tip_m),
            annulus_area_m2: PI * (r_tip_m * r_tip_m - r_hub_m * r_hub_m),
        });
    };

    // Stage 0: Curtis impulse control wheel (nozzle box upstream at z = 1.24, rotor at z = 1.30)
    push_row(0, 6, false, 1.235, 0.342, 0.390, 48, 0.034, true);
    push_row(0, 0, true, 1.298, 0.340, 0.392, 56, 0.044, true);

    // Stages 1..=5: HP 5-stage reaction drum (z = 1.45..1.98)
    for k in 0..5 {
        let f = k as f64 / 4.0;
        let z_stator = 1.455 + k as f64 * 0.112;
        let z_rotor = z_stator + 0.052;
        let r_hub = 0.275;
        let r_tip = 0.336 + 0.042 * f;
        push_row(1 + k, 6, false, z_stator, r_hub + 0.003, r_tip, 60, 0.032, true);
        push_row(1 + k, 1, true, z_rotor, r_hub, r_tip - 0.002, 64, 0.036, true);
    }

    // Stages 6..=9: IP1 4-stage reaction drum (after controlled extraction 1 scroll, z = 2.29..2.72)
    for k in 0..4 {
        let f = k as f64 / 3.0;
        let z_stator = 2.285 + k as f64 * 0.112;
        let z_rotor = z_stator + 0.052;
        let r_hub = 0.286 + 0.009 * f;
        let r_tip = 0.392 + 0.052 * f;
        push_row(6 + k, 6, false, z_stator, r_hub + 0.003, r_tip, 56, 0.036, true);
        push_row(6 + k, 2, true, z_rotor, r_hub, r_tip - 0.002, 60, 0.040, true);
    }

    // Stages 10..=13: IP2 4-stage reaction drum (after bleed scroll, z = 2.96..3.36)
    for k in 0..4 {
        let f = k as f64 / 3.0;
        let z_stator = 2.955 + k as f64 * 0.106;
        let z_rotor = z_stator + 0.050;
        let r_hub = 0.299 + 0.024 * f;
        let r_tip = 0.465 + 0.088 * f;
        push_row(10 + k, 6, false, z_stator, r_hub + 0.003, r_tip, 52, 0.040, true);
        push_row(10 + k, 3, true, z_rotor, r_hub, r_tip - 0.002, 56, 0.044, true);
    }

    // Stages 14..=15: LP flared penultimate stages (z = 3.41..3.58)
    push_row(14, 6, false, 3.395, 0.332, 0.610, 48, 0.045, true);
    push_row(14, 4, true, 3.450, 0.334, 0.618, 50, 0.050, true);
    push_row(15, 6, false, 3.515, 0.352, 0.715, 46, 0.050, true);
    push_row(15, 4, true, 3.575, 0.356, 0.725, 48, 0.056, false);

    // Stage 16: Giant LP last-stage free-standing twisted blades (iconic focal point in image-0031!)
    push_row(16, 6, false, 3.655, 0.382, 0.905, 42, 0.068, true);
    push_row(16, 5, true, 3.775, 0.408, 0.935, 44, 0.118, false);

    rows
}

/// Build the complete FrankenSim geometry catalog once, returning
/// GPU-ready meshes, verified axisymmetric receipts, and exact rotor inertia properties.
pub fn build_turbine_geometry(cx: &Cx<'_>) -> TurbineGeometryCatalog {
    let lib = library();

    // 1. Forged CrMoV rotor shaft (exact Green's-theorem mass/inertia + tessellation)
    let rotor_chart = build_rotor_chart();
    let shaft_mass_props = rotor_chart
        .mass_properties(lib.rotor_steel.rho_kg_m3, cx)
        .expect("exact rotor mass properties");
    let shaft_surface_area_m2 = rotor_chart
        .surface_area(cx)
        .expect("exact rotor surface area")
        .area;
    let mut rotor_cfg =
        AxisymmetricTessellationConfig::new(0.0012, TessellationPurpose::Rendering)
            .expect("rotor tessellation config");
    rotor_cfg.min_azimuthal_sectors = 56;
    rotor_cfg.max_azimuthal_sectors = 72;
    let rotor_raw = tessellate_axisymmetric(&rotor_chart, rotor_cfg, cx)
        .expect("tessellate forged rotor");
    let rotor_receipt = rotor_raw.receipt.clone();
    let rotor_mesh = axisymmetric_to_export("rotor_forging", rotor_raw);

    // 2. Annular Inconel balance-piston labyrinth seals
    let lab_chart = build_balance_labyrinth_chart();
    let mut lab_cfg =
        AxisymmetricTessellationConfig::new(0.0018, TessellationPurpose::Rendering)
            .expect("labyrinth tessellation config");
    lab_cfg.min_azimuthal_sectors = 48;
    lab_cfg.max_azimuthal_sectors = 64;
    let lab_raw = tessellate_axisymmetric(&lab_chart, lab_cfg, cx)
        .expect("tessellate labyrinth seals");
    let labyrinth_seals = axisymmetric_to_export("labyrinth_seals", lab_raw);

    // 3. Stationary housings as validated AxisymmetricCharts + 3 solid cutaway sub-shells each
    let inner_carrier = build_cutaway_housing(
        "inner_blade_carrier",
        "inner_carrier_lower",
        "inner_carrier_upper_back",
        "inner_carrier_upper_front",
        &build_inner_carrier_chart(),
        lib.hp_casing_steel.rho_kg_m3,
        0.0020,
        16,
        cx,
    );
    let outer_hp_ip_casing = build_cutaway_housing(
        "outer_hp_ip_casing",
        "outer_casing_lower",
        "outer_casing_upper_back",
        "outer_casing_upper_front",
        &build_outer_hp_ip_casing_chart(),
        lib.hp_casing_steel.rho_kg_m3,
        0.0025,
        18,
        cx,
    );
    let exhaust_diffuser = build_cutaway_housing(
        "exhaust_diffuser_hood",
        "exhaust_diffuser_lower",
        "exhaust_diffuser_upper_back",
        "exhaust_diffuser_upper_front",
        &build_exhaust_diffuser_chart(),
        lib.lp_exhaust_steel.rho_kg_m3,
        0.0030,
        20,
        cx,
    );
    let rear_bearing_cone = build_cutaway_housing(
        "rear_bearing_cone",
        "rear_cone_lower",
        "rear_cone_upper_back",
        "rear_cone_upper_front",
        &build_rear_bearing_cone_chart(),
        lib.lp_exhaust_steel.rho_kg_m3,
        0.0020,
        14,
        cx,
    );
    let front_bearing_housing = build_cutaway_housing(
        "front_bearing_housing",
        "front_bearing_lower",
        "front_bearing_upper_back",
        "front_bearing_upper_front",
        &build_front_bearing_housing_chart(),
        lib.pedestal_cast_iron.rho_kg_m3,
        0.0020,
        14,
        cx,
    );

    // 4. 7 twisted, cambered 3D airfoil families generated via fs_rep_nurbs::NurbsSurface
    let airfoil_specs = [
        AirfoilFamilyParams {
            name: "nurbs_blade_curtis",
            stagger_root_deg: 22.0,
            stagger_tip_deg: 28.0,
            camber_root: 0.24,
            camber_tip: 0.22,
            thickness_root: 0.24,
            thickness_tip: 0.20,
            chord_taper: 0.94,
            lean_tip: 0.0,
            camber_sign: 1.0,
            nu: 14,
            nv: 6,
        },
        AirfoilFamilyParams {
            name: "nurbs_blade_hp",
            stagger_root_deg: 32.0,
            stagger_tip_deg: 40.0,
            camber_root: 0.16,
            camber_tip: 0.14,
            thickness_root: 0.20,
            thickness_tip: 0.16,
            chord_taper: 0.92,
            lean_tip: 0.02,
            camber_sign: 1.0,
            nu: 12,
            nv: 5,
        },
        AirfoilFamilyParams {
            name: "nurbs_blade_ip1",
            stagger_root_deg: 30.0,
            stagger_tip_deg: 48.0,
            camber_root: 0.15,
            camber_tip: 0.12,
            thickness_root: 0.18,
            thickness_tip: 0.13,
            chord_taper: 0.84,
            lean_tip: 0.04,
            camber_sign: 1.0,
            nu: 12,
            nv: 6,
        },
        AirfoilFamilyParams {
            name: "nurbs_blade_ip2",
            stagger_root_deg: 28.0,
            stagger_tip_deg: 56.0,
            camber_root: 0.14,
            camber_tip: 0.10,
            thickness_root: 0.16,
            thickness_tip: 0.10,
            chord_taper: 0.74,
            lean_tip: 0.06,
            camber_sign: 1.0,
            nu: 14,
            nv: 7,
        },
        AirfoilFamilyParams {
            name: "nurbs_blade_lp_mid",
            stagger_root_deg: 25.0,
            stagger_tip_deg: 64.0,
            camber_root: 0.13,
            camber_tip: 0.07,
            thickness_root: 0.15,
            thickness_tip: 0.08,
            chord_taper: 0.64,
            lean_tip: 0.08,
            camber_sign: 1.0,
            nu: 14,
            nv: 8,
        },
        AirfoilFamilyParams {
            name: "nurbs_blade_lp_last",
            stagger_root_deg: 22.0,
            stagger_tip_deg: 78.0, // Dramatic 56° free-vortex twist of a large last stage!
            camber_root: 0.14,
            camber_tip: 0.05,
            thickness_root: 0.15,
            thickness_tip: 0.055,
            chord_taper: 0.54,
            lean_tip: 0.14,
            camber_sign: 1.0,
            nu: 18,
            nv: 12,
        },
        AirfoilFamilyParams {
            name: "nurbs_vane_stator",
            stagger_root_deg: -36.0,
            stagger_tip_deg: -46.0,
            camber_root: 0.16,
            camber_tip: 0.13,
            thickness_root: 0.18,
            thickness_tip: 0.14,
            chord_taper: 0.90,
            lean_tip: 0.0,
            camber_sign: -1.0,
            nu: 12,
            nv: 5,
        },
    ];
    let nurbs_airfoils: Vec<ExportMesh> = airfoil_specs
        .into_iter()
        .map(build_nurbs_airfoil_mesh)
        .collect();

    // 5. Compute exact combined mass + inertia of shaft + all 17 rotor NURBS blade rows
    let blade_rows = build_blade_rows();
    let mut blades_mass_kg = 0.0;
    let mut blades_mz = 0.0;
    let mut blades_jp = 0.0;
    let mut blades_jd_origin = 0.0;

    for row in &blade_rows {
        if !row.is_rotor {
            continue;
        }
        let rho_b = if row.stage_index >= 14 {
            lib.lp_blade_steel.rho_kg_m3
        } else {
            lib.hp_blade_steel.rho_kg_m3
        };
        let height = row.r_tip_m - row.r_hub_m;
        // Mean airfoil section thickness factor (~0.10 * chord^2) × taper factor (~0.75)
        let mean_section_area = 0.095 * row.axial_chord_m * row.axial_chord_m;
        let blade_mass = rho_b * mean_section_area * height;
        let shroud_mass = if row.has_shroud {
            rho_b * TAU * row.r_tip_m * row.axial_chord_m * 0.0055
        } else {
            0.0
        };
        let row_mass = blade_mass * f64::from(row.blade_count) + shroud_mass;
        // Thin-annulus radial inertia of the blade ring: I_z = m * (r_hub^2 + r_hub*r_tip + r_tip^2)/3
        let r2_mean = (row.r_hub_m * row.r_hub_m
            + row.r_hub_m * row.r_tip_m
            + row.r_tip_m * row.r_tip_m)
            / 3.0;
        blades_mass_kg += row_mass;
        blades_mz += row_mass * row.z_center_m;
        blades_jp += row_mass * r2_mean;
        // Parallel-axis contribution about z = 0: I_x(0) = m * (r2_mean / 2 + z^2)
        blades_jd_origin += row_mass * (0.5 * r2_mean + row.z_center_m * row.z_center_m);
    }

    let total_mass_kg = shaft_mass_props.mass + blades_mass_kg;
    let com_z_m =
        (shaft_mass_props.mass * shaft_mass_props.center_of_mass.z + blades_mz) / total_mass_kg;
    let polar_inertia_kg_m2 = shaft_mass_props.principal_inertia.axial + blades_jp;
    // Shift combined origin transverse inertia to combined COM via parallel axis theorem
    let combined_jd_origin = shaft_mass_props.origin_inertia.transverse + blades_jd_origin;
    let transverse_inertia_kg_m2 = (combined_jd_origin - total_mass_kg * com_z_m * com_z_m).max(1.0);

    let rotor_props = RotorGeometryProps {
        shaft_mass_kg: shaft_mass_props.mass,
        blades_mass_kg,
        total_mass_kg,
        com_z_m,
        polar_inertia_kg_m2,
        transverse_inertia_kg_m2,
        shaft_surface_area_m2,
        front_bearing_rel_z_m: FRONT_BEARING_Z_M - com_z_m,
        rear_bearing_rel_z_m: REAR_BEARING_Z_M - com_z_m,
        rotor_vertex_count: rotor_receipt.vertex_count,
        rotor_triangle_count: rotor_receipt.triangle_count,
        rotor_hausdorff_bound_m: rotor_receipt.total_hausdorff_bound,
        rotor_euler_characteristic: rotor_receipt.euler_characteristic,
        rotor_watertight: rotor_receipt.is_watertight && rotor_receipt.is_outward_oriented,
    };

    let shell_tris = |s: &CutawayShellMeshes| {
        s.lower_half.triangle_count()
            + s.upper_back.triangle_count()
            + s.upper_front.triangle_count()
    };
    let total_base_triangles = rotor_mesh.triangle_count()
        + labyrinth_seals.triangle_count()
        + shell_tris(&inner_carrier)
        + shell_tris(&outer_hp_ip_casing)
        + shell_tris(&exhaust_diffuser)
        + shell_tris(&rear_bearing_cone)
        + shell_tris(&front_bearing_housing)
        + nurbs_airfoils.iter().map(ExportMesh::triangle_count).sum::<usize>();

    let instanced_blade_tris: usize = blade_rows
        .iter()
        .map(|row| nurbs_airfoils[row.family_id].triangle_count() * row.blade_count as usize)
        .sum();
    let total_scene_triangles = total_base_triangles + instanced_blade_tris;

    TurbineGeometryCatalog {
        rotor_props,
        rotor_mesh,
        labyrinth_seals,
        inner_carrier,
        outer_hp_ip_casing,
        exhaust_diffuser,
        rear_bearing_cone,
        front_bearing_housing,
        nurbs_airfoils,
        blade_rows,
        total_base_triangles,
        total_scene_triangles,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::with_clock_free_cx;
    use fs_geom::Point3;
    use fs_rep_mesh::Soup;
    use fs_topo::manifold_certificate;

    #[test]
    fn rotor_and_housings_are_watertight_axisymmetric_solids_with_physical_inertia() {
        let cat = with_clock_free_cx(build_turbine_geometry);

        // Forged CrMoV rotor + 17 blade rows: realistic industrial 50 MW rotor
        // (~6,000..12,000 kg total mass, J_p ~ 400..1,200 kg·m², J_d ~ 6,000..25,000 kg·m²)
        assert!(cat.rotor_props.rotor_watertight);
        assert_eq!(cat.rotor_props.rotor_euler_characteristic, 2);
        assert!(cat.rotor_props.rotor_hausdorff_bound_m < 0.002);
        assert!(
            cat.rotor_props.shaft_mass_kg > 4_000.0 && cat.rotor_props.shaft_mass_kg < 14_000.0,
            "unexpected shaft mass: {} kg",
            cat.rotor_props.shaft_mass_kg
        );
        assert!(
            cat.rotor_props.blades_mass_kg > 200.0 && cat.rotor_props.blades_mass_kg < 2_500.0,
            "unexpected blades mass: {} kg",
            cat.rotor_props.blades_mass_kg
        );
        assert!(
            cat.rotor_props.polar_inertia_kg_m2 > 200.0
                && cat.rotor_props.polar_inertia_kg_m2 < 2_500.0,
            "unexpected J_p: {} kg·m²",
            cat.rotor_props.polar_inertia_kg_m2
        );
        assert!(
            cat.rotor_props.transverse_inertia_kg_m2 > cat.rotor_props.polar_inertia_kg_m2,
            "slender steam turbine rotor must satisfy J_d > J_p"
        );
        assert!(
            cat.rotor_props.com_z_m > FRONT_BEARING_Z_M && cat.rotor_props.com_z_m < REAR_BEARING_Z_M,
            "COM must sit between front and rear bearings"
        );

        // Verify every stationary annular housing has bounded Hausdorff sagitta, positive
        // Green's-theorem mass & surface area, and non-degenerate capped cutaway sub-shells.
        for (name, shell) in [
            ("inner_carrier", &cat.inner_carrier),
            ("outer_hp_ip_casing", &cat.outer_hp_ip_casing),
            ("exhaust_diffuser", &cat.exhaust_diffuser),
            ("rear_bearing_cone", &cat.rear_bearing_cone),
            ("front_bearing_housing", &cat.front_bearing_housing),
        ] {
            assert!(
                shell.full_receipt.total_hausdorff_bound < 0.005
                    && shell.full_receipt.triangle_count > 400,
                "{name} tessellation receipt unexpected: {:?}",
                shell.full_receipt
            );
            assert!(shell.mass_props.mass > 100.0, "{name} mass too small");
            assert!(shell.surface_area_m2 > 0.5, "{name} surface area too small");
            assert!(shell.lower_half.triangle_count() > 200);
            assert!(shell.upper_back.triangle_count() > 100);
            assert!(shell.upper_front.triangle_count() > 100);
        }

        // Verify all 7 NURBS airfoil families have finite unit normals & non-degenerate triangles
        assert_eq!(cat.nurbs_airfoils.len(), 7);
        for airfoil in &cat.nurbs_airfoils {
            assert!(airfoil.triangle_count() > 100, "{} too coarse", airfoil.name);
            for chunk in airfoil.normals.chunks_exact(3) {
                let norm = (chunk[0] * chunk[0] + chunk[1] * chunk[1] + chunk[2] * chunk[2]).sqrt();
                assert!(
                    (norm - 1.0).abs() < 1.0e-3,
                    "{} non-unit normal {norm}",
                    airfoil.name
                );
            }
        }

        // Run fs-topo combinatorial manifold certificate on the tessellated forged rotor!
        let positions: Vec<Point3> = cat
            .rotor_mesh
            .positions
            .chunks_exact(3)
            .map(|c| Point3::new(f64::from(c[0]), f64::from(c[1]), f64::from(c[2])))
            .collect();
        let triangles: Vec<[u32; 3]> = cat
            .rotor_mesh
            .indices
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        let soup = Soup {
            positions,
            triangles,
        };
        let cert = manifold_certificate(&soup, Some(Point3::new(0.0, 0.0, 2.0)));
        assert!(
            cert.defects.is_empty(),
            "rotor mesh manifold_certificate defects: {:?}",
            cert.defects
        );
    }
}
