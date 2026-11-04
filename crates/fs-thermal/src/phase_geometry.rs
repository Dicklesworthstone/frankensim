//! Mass-owned, constant-density phase fronts and exact thermostatted retreat.
//!
//! A planar layer retreats through one face, a sphere radially, and a cylinder
//! radially through its curved wall with inert end caps. Curved coatings retain
//! an explicitly inert core and evolve only the deposited shell. These reductions assume
//! uniform flux and shape preservation. They are not a mesh-motion, deformation,
//! surface-energy, Kelvin-curvature, nucleation or contact-topology model.
//!
//! The original sublimation owner supplies the signed local flux. At fixed
//! surface/bath temperatures and prescribed vapor pressure, dx/dt=-j/rho is
//! constant. Integrating x, rather than clipping A*j*dt, resolves changing area
//! and the exact extinction time. The thermostat supplies latent heat explicitly;
//! this is NOT complete evaporation of an unheated finite-energy solid.

use crate::sublimation::{
    SublimationError, SublimationLaw, SurfaceKinetics, VaporReservoir, surface_flux,
};
use core::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Shape {
    Layer { footprint_m2: f64 },
    Sphere,
    Cylinder { length_m: f64 },
    SphericalCoating { core_radius_m: f64 },
    CylindricalCoating { core_radius_m: f64, length_m: f64 },
}

/// Fixed morphology and solid density. Remaining mass owns all moving dimensions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetreatGeometry {
    shape: Shape,
    density_kg_m3: f64,
}

/// Actual residual geometry; an extinct phase has exactly zero active area.
/// A layer's underlying substrate is NOT created or assigned a new boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrontGeometry {
    pub mass_kg: f64,
    pub volume_m3: f64,
    pub area_m2: f64,
    /// Layer/coating thickness, or radius of a core-free sphere/cylinder [m].
    pub extent_m: f64,
    /// dA/dm on the positive-mass branch; undefined at extinction.
    pub area_mass_derivative_m2_kg: Option<f64>,
}

fn invalid(message: &'static str) -> SublimationError { SublimationError::Invalid(message) }
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }

impl RetreatGeometry {
    /// One exposed face, fixed footprint, no lateral loss or substrate inference.
    pub fn layer(density_kg_m3: f64, footprint_m2: f64) -> Result<Self, SublimationError> {
        if !positive(footprint_m2) { return Err(invalid("positive phase-layer footprint required")); }
        Self::new(Shape::Layer { footprint_m2 }, density_kg_m3)
    }
    /// Fully exposed, uniformly retreating/growing sphere.
    pub fn sphere(density_kg_m3: f64) -> Result<Self, SublimationError> {
        Self::new(Shape::Sphere, density_kg_m3)
    }
    /// Uniform radial retreat; end caps have no phase or reflected-heat exchange.
    pub fn cylinder(density_kg_m3: f64, length_m: f64) -> Result<Self, SublimationError> {
        if !positive(length_m) { return Err(invalid("positive phase-cylinder length required")); }
        Self::new(Shape::Cylinder { length_m }, density_kg_m3)
    }
    /// Uniform coating on a fixed spherical substrate. Mass and volume describe
    /// only the deposited phase; exposed area is at the outer shell radius.
    /// The substrate's mass, energy and displacement remain caller-owned.
    pub fn spherical_coating(density_kg_m3: f64, core_radius_m: f64) -> Result<Self, SublimationError> {
        if !positive(core_radius_m) { return Err(invalid("positive spherical substrate radius required")); }
        Self::new(Shape::SphericalCoating { core_radius_m }, density_kg_m3)
    }
    /// Uniform radial coating of a fixed cylindrical substrate, with inert end
    /// caps. Extent is coating thickness, not the substrate-plus-coating radius.
    /// A positive seed is required by the existing phase-evolution owners;
    /// heterogeneous nucleation and heat storage in the substrate are not inferred.
    pub fn cylindrical_coating(density_kg_m3: f64, core_radius_m: f64, length_m: f64) -> Result<Self, SublimationError> {
        if !positive(core_radius_m) || !positive(length_m) {
            return Err(invalid("positive cylindrical substrate radius and length required"));
        }
        Self::new(Shape::CylindricalCoating { core_radius_m, length_m }, density_kg_m3)
    }
    fn new(shape: Shape, density_kg_m3: f64) -> Result<Self, SublimationError> {
        if !positive(density_kg_m3) { return Err(invalid("positive constant solid density required")); }
        Ok(Self { shape, density_kg_m3 })
    }
    #[must_use]
    pub const fn density_kg_m3(self) -> f64 { self.density_kg_m3 }

    /// Physical displaced volume [m3], including an inert coating substrate.
    /// Use this for gas-volume/diluteness admission. `FrontGeometry::volume_m3`
    /// instead counts only the phase mass that can deposit or sublime.
    pub fn occupied_volume_m3(self, mass_kg: f64) -> Result<f64, SublimationError> {
        if !mass_kg.is_finite() || mass_kg < 0.0 {
            return Err(invalid("nonnegative finite phase mass required"));
        }
        // Volume admission needs no radius/thickness inversion; moving owners
        // already derived the exposed geometry once for this same endpoint.
        let deposit = mass_kg / self.density_kg_m3;
        if !deposit.is_finite() || mass_kg > 0.0 && deposit == 0.0 {
            return Err(invalid("phase volume is unrepresentable"));
        }
        let core = match self.shape {
            Shape::SphericalCoating { core_radius_m } => {
                ((4.0 * PI / 3.0 * core_radius_m) * core_radius_m) * core_radius_m
            }
            Shape::CylindricalCoating { core_radius_m, length_m } => {
                ((PI * core_radius_m) * core_radius_m) * length_m
            }
            _ => 0.0,
        };
        let occupied = core + deposit;
        if !core.is_finite() || !occupied.is_finite() {
            return Err(invalid("occupied phase/substrate volume is unrepresentable"));
        }
        Ok(occupied)
    }

    /// Check the original phase law, rather than silently adding thermal expansion.
    pub fn admit_density(self, density: f64) -> Result<(), SublimationError> {
        if !positive(density) || (density / self.density_kg_m3 - 1.0).abs() > 64.0 * f64::EPSILON {
            Err(invalid("moving front requires the declared constant phase density"))
        } else { Ok(()) }
    }

    /// Derive geometry from the conserved mass; no independent area state drifts.
    pub fn at_mass(self, mass_kg: f64) -> Result<FrontGeometry, SublimationError> {
        if !mass_kg.is_finite() || mass_kg < 0.0 { return Err(invalid("nonnegative finite phase mass required")); }
        if mass_kg == 0.0 {
            return Ok(FrontGeometry { mass_kg: 0.0, volume_m3: 0.0, area_m2: 0.0,
                extent_m: 0.0, area_mass_derivative_m2_kg: None });
        }
        let volume_m3 = mass_kg / self.density_kg_m3;
        let (extent_m, area_m2, derivative) = match self.shape {
            Shape::Layer { footprint_m2 } => (volume_m3 / footprint_m2, footprint_m2, 0.0),
            Shape::Sphere => {
                let r = (volume_m3 / (4.0 * PI / 3.0)).cbrt();
                (r, (4.0 * PI * r) * r, (2.0 / self.density_kg_m3) / r)
            }
            Shape::Cylinder { length_m } => {
                let r = (volume_m3 / length_m / PI).sqrt();
                (r, (2.0 * PI * r) * length_m, (1.0 / self.density_kg_m3) / r)
            }
            Shape::SphericalCoating { core_radius_m } => {
                // Scale before cubing, then rationalize r - r_core. In
                // particular, a film thinner than one ULP of the core radius
                // still has a representable thickness and conserved mass.
                let bulk_radius = (volume_m3 / (4.0 * PI / 3.0)).cbrt();
                let scale = core_radius_m.max(bulk_radius);
                let a = core_radius_m / scale;
                let b = bulk_radius / scale;
                let r = scale * (a * a * a + b * b * b).cbrt();
                let ratio = core_radius_m / r;
                let thickness = (volume_m3 / r / r)
                    / ((4.0 * PI / 3.0) * (1.0 + ratio + ratio * ratio));
                (thickness, (4.0 * PI * r) * r, (2.0 / self.density_kg_m3) / r)
            }
            Shape::CylindricalCoating { core_radius_m, length_m } => {
                let bulk_radius = (volume_m3 / length_m / PI).sqrt();
                let scale = core_radius_m.max(bulk_radius);
                let a = core_radius_m / scale;
                let b = bulk_radius / scale;
                let r = scale * (a * a + b * b).sqrt();
                let thickness = (volume_m3 / length_m / PI / r)
                    / (1.0 + core_radius_m / r);
                (thickness, (2.0 * PI * r) * length_m, (1.0 / self.density_kg_m3) / r)
            }
        };
        if ![volume_m3, extent_m, area_m2].into_iter().all(positive)
            || !derivative.is_finite() || derivative < 0.0 {
            return Err(invalid("phase geometry or tangent is unrepresentable"));
        }
        Ok(FrontGeometry { mass_kg, volume_m3, area_m2, extent_m,
            area_mass_derivative_m2_kg: Some(derivative) })
    }

    /// Construct the same geometry from an explicitly supplied radius/thickness.
    pub fn at_extent(self, extent_m: f64) -> Result<FrontGeometry, SublimationError> {
        if !extent_m.is_finite() || extent_m < 0.0 { return Err(invalid("nonnegative finite phase extent required")); }
        if extent_m == 0.0 { return self.at_mass(0.0); }
        let volume = match self.shape {
            Shape::Layer { footprint_m2 } => footprint_m2 * extent_m,
            Shape::Sphere => ((4.0 * PI / 3.0 * extent_m) * extent_m) * extent_m,
            Shape::Cylinder { length_m } => ((PI * extent_m) * extent_m) * length_m,
            Shape::SphericalCoating { core_radius_m } => {
                let r = core_radius_m + extent_m;
                // Factored difference of cubes; never subtract near-equal
                // total/core volumes to recover a thin deposited shell.
                ((4.0 * PI / 3.0) * extent_m)
                    * (r * r + r * core_radius_m + core_radius_m * core_radius_m)
            }
            Shape::CylindricalCoating { core_radius_m, length_m } => {
                ((PI * extent_m) * (2.0 * core_radius_m + extent_m)) * length_m
            }
        };
        let mass = volume * self.density_kg_m3;
        if !positive(volume) || !positive(mass) { return Err(invalid("phase extent has unrepresentable mass")); }
        let mut result = self.at_mass(mass)?;
        // Retain the supplied coordinate, not a needless cube-root roundtrip.
        result.extent_m = extent_m;
        Ok(result)
    }

    /// Positive normal retreat, negative growth. Shell identities preserve tiny
    /// changes without subtracting nearly equal square/cube roots.
    pub fn recession_m(self, initial_mass_kg: f64, final_mass_kg: f64) -> Result<f64, SublimationError> {
        let a = self.at_mass(initial_mass_kg)?;
        let b = self.at_mass(final_mass_kg)?;
        if initial_mass_kg == final_mass_kg { return Ok(0.0); }
        let delta_volume = (initial_mass_kg - final_mass_kg) / self.density_kg_m3;
        let x = a.extent_m;
        let y = b.extent_m;
        let shell = match self.shape {
            Shape::Layer { footprint_m2 } => footprint_m2,
            Shape::Cylinder { length_m } => PI * length_m * (x + y),
            Shape::Sphere => (4.0 * PI / 3.0) * (x * x + x * y + y * y),
            Shape::SphericalCoating { core_radius_m } => {
                let x = core_radius_m + x;
                let y = core_radius_m + y;
                (4.0 * PI / 3.0) * (x * x + x * y + y * y)
            }
            Shape::CylindricalCoating { core_radius_m, length_m } => {
                PI * length_m * (2.0 * core_radius_m + x + y)
            }
        };
        let result = delta_volume / shell;
        if !positive(shell) || !result.is_finite() || result == 0.0 {
            Err(invalid("phase-front displacement is unrepresentable"))
        } else { Ok(result) }
    }
}

/// One exact geometric interval under explicitly thermostatted conditions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermostattedFrontStep {
    pub initial: FrontGeometry,
    pub final_geometry: FrontGeometry,
    pub transferred_mass_kg: f64,
    /// Positive normal retreat, negative frost growth.
    pub recession_m: f64,
    /// Effective integral of exposed area over active time, m2 s.
    pub area_time_m2_s: f64,
    /// Relative event time, only when extinction occurs in this call.
    pub extinction_time_s: Option<f64>,
    pub active_duration_s: f64,
    pub unused_duration_s: f64,
    /// Supplied by the thermostat; negative during deposition.
    pub thermostat_heat_j: f64,
    /// Signed vapor enthalpy in the supplied solid datum, includes latent heat.
    pub vapor_energy_j: f64,
    pub solid_energy_change_j: f64,
    pub energy_residual_j: f64,
}

/// Uniform fixed-temperature phase front in a prescribed infinite vapor bath.
///
/// Calls the existing Hertz--Knudsen/optional dilute-film owner, then integrates
/// exact layer/cylinder/sphere kinematics. All states and areas use the same mass.
/// Extinction stops this interval; no phantom heat, negative mass, residual
/// cutoff, automatic renucleation or post-extinction substrate law is introduced.
/// Finite reservoirs must use an evolving transport owner, not this frozen bath.
/// The vapor-energy convention is the original thermostatted enthalpy convention,
/// NOT the separate nonisothermal molecular stream energy. Neglects pressure work
/// and surface free energy, as does that original thermostatted surface model.
#[allow(clippy::too_many_arguments)]
pub fn advance_thermostatted_front<L: SublimationLaw + ?Sized>(
    law: &L,
    geometry: RetreatGeometry,
    initial_mass_kg: f64,
    surface_temperature_k: f64,
    solid_specific_enthalpy_j_kg: f64,
    kinetics: SurfaceKinetics,
    reservoir: VaporReservoir,
    duration_s: f64,
) -> Result<ThermostattedFrontStep, SublimationError> {
    if !positive(initial_mass_kg) || !duration_s.is_finite() || duration_s < 0.0
        || !solid_specific_enthalpy_j_kg.is_finite()
        || !matches!(reservoir, VaporReservoir::Prescribed { .. }) {
        return Err(invalid("positive seed mass, prescribed bath and finite nonnegative interval required"));
    }
    let initial = geometry.at_mass(initial_mass_kg)?;
    let flux = surface_flux(law, surface_temperature_k, kinetics, reservoir)?;
    let properties = law.properties(surface_temperature_k)?;
    geometry.admit_density(properties.solid_density_kg_m3)?;
    let j = flux.mass_flux_kg_m2_s;
    let speed = j / geometry.density_kg_m3;
    if !speed.is_finite() || (j != 0.0 && speed == 0.0) { return Err(invalid("phase-front speed is unrepresentable")); }
    let time_to_zero = if speed > 0.0 { initial.extent_m / speed } else { f64::INFINITY };
    let extinct = speed > 0.0 && duration_s >= time_to_zero;
    let active_duration_s = if extinct { time_to_zero } else { duration_s };
    if extinct && !positive(active_duration_s) { return Err(invalid("extinction time is unrepresentable")); }
    let final_geometry = if extinct {
        geometry.at_mass(0.0)?
    } else if speed == 0.0 || duration_s == 0.0 {
        initial
    } else {
        let displacement = speed * duration_s;
        let extent = initial.extent_m - displacement;
        if !positive(extent) || extent == initial.extent_m { return Err(invalid("phase-front increment is unrepresentable")); }
        geometry.at_extent(extent)?
    };
    let dm = initial.mass_kg - final_geometry.mass_kg;
    if speed != 0.0 && duration_s != 0.0 && (dm == 0.0 || dm.signum() != speed.signum()) {
        return Err(invalid("moving-front mass increment is unrepresentable"));
    }
    let recession_m = geometry.recession_m(initial.mass_kg, final_geometry.mass_kg)?;
    let area_time_m2_s = if j == 0.0 { initial.area_m2 * duration_s } else { dm / j };
    let thermostat_heat_j = dm * properties.latent_heat_j_kg;
    let solid_energy_change_j = (final_geometry.mass_kg - initial.mass_kg) * solid_specific_enthalpy_j_kg;
    let vapor_energy_j = dm * (solid_specific_enthalpy_j_kg + properties.latent_heat_j_kg);
    let energy_residual_j = solid_energy_change_j + vapor_energy_j - thermostat_heat_j;
    if ![recession_m, area_time_m2_s, thermostat_heat_j, solid_energy_change_j,
        vapor_energy_j, energy_residual_j].into_iter().all(f64::is_finite)
        || area_time_m2_s < 0.0 {
        return Err(invalid("moving-front energy or exposure is unrepresentable"));
    }
    Ok(ThermostattedFrontStep { initial, final_geometry, transferred_mass_kg: dm, recession_m,
        area_time_m2_s, extinction_time_s: extinct.then_some(active_duration_s), active_duration_s,
        unused_duration_s: duration_s - active_duration_s, thermostat_heat_j, vapor_energy_j,
        solid_energy_change_j, energy_residual_j })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sublimation::SublimationProperties;
    struct Phase;
    impl SublimationLaw for Phase {
        fn temperature_range_k(&self) -> [f64; 2] { [200.0, 400.0] }
        fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
            Ok(SublimationProperties { saturation_pressure_pa: 100.0,
                saturation_pressure_derivative_pa_k: 0.0, latent_heat_j_kg: 2e6,
                latent_heat_derivative_j_kg_k: 0.0, solid_density_kg_m3: 1000.0,
                molar_mass_kg_mol: 0.018 })
        }
    }
    fn shapes() -> [RetreatGeometry; 3] {
        [RetreatGeometry::layer(1000.0, 0.02).unwrap(), RetreatGeometry::cylinder(1000.0, 0.1).unwrap(),
            RetreatGeometry::sphere(1000.0).unwrap()]
    }
    fn near(a: f64, b: f64, scale: f64) { assert!((a-b).abs() <= scale, "{a:e} != {b:e}"); }
    fn run(g: RetreatGeometry, m: f64, p: f64, dt: f64) -> ThermostattedFrontStep {
        advance_thermostatted_front(&Phase, g, m, 300.0, 123456.0,
            SurfaceKinetics { accommodation_coefficient: 0.5, mass_transfer_coefficient_m_s: None },
            VaporReservoir::Prescribed { partial_pressure_pa: p, temperature_k: 300.0 }, dt).unwrap()
    }
    #[test]
    fn exact_geometry_tangents_and_small_shell_displacements() {
        for (i, g) in shapes().into_iter().enumerate() {
            let r = 0.003;
            let s = g.at_extent(r).unwrap();
            let (v, area) = match i {
                0 => (0.02*r, 0.02), 1 => (PI*r*r*0.1, 2.0*PI*r*0.1),
                _ => (4.0*PI*r*r*r/3.0, 4.0*PI*r*r),
            };
            near(s.mass_kg, 1000.0*v, 1e-15);
            near(s.area_m2, area, 1e-15);
            let dm = s.mass_kg * 1e-5;
            let fd = (g.at_mass(s.mass_kg+dm).unwrap().area_m2 - g.at_mass(s.mass_kg-dm).unwrap().area_m2)/(2.0*dm);
            near(fd, s.area_mass_derivative_m2_kg.unwrap(), 1e-8);
            let next = s.mass_kg.next_down();
            assert!(g.recession_m(s.mass_kg, next).unwrap() > 0.0);
            assert_eq!(g.recession_m(next, s.mass_kg).unwrap(), -g.recession_m(s.mass_kg, next).unwrap());
            assert_eq!(g.at_mass(0.0).unwrap().area_m2, 0.0);
        }
    }
    #[test]
    fn radius_retreat_changes_mass_and_exposure_for_all_morphologies() {
        let j = 0.5 * 100.0 * (0.018 / (2.0*PI*crate::sublimation::MOLAR_GAS_CONSTANT*300.0)).sqrt();
        for (i, g) in shapes().into_iter().enumerate() {
            let start = g.at_extent(0.003).unwrap();
            let dt = 0.2 * start.extent_m * 1000.0 / j;
            let step = run(g, start.mass_kg, 0.0, dt);
            near(step.final_geometry.extent_m, 0.8*start.extent_m, 1e-16);
            let power = i as i32 + 1;
            near(step.final_geometry.mass_kg/start.mass_kg, 0.8_f64.powi(power), 1e-14);
            near(step.area_time_m2_s*j, step.transferred_mass_kg, 1e-15);
            near(step.energy_residual_j, 0.0, 1e-9);
            assert_eq!(step.extinction_time_s, None);
        }
    }
    #[test]
    fn extinction_is_localized_and_never_spends_post_event_heat() {
        for g in shapes() {
            let start = g.at_extent(0.001).unwrap();
            let early = run(g, start.mass_kg, 0.0, 0.01);
            let speed = early.recession_m / 0.01;
            let lifetime = start.extent_m / speed;
            let end = run(g, start.mass_kg, 0.0, 2.0*lifetime);
            assert_eq!(end.final_geometry.mass_kg, 0.0);
            assert_eq!(end.final_geometry.area_m2, 0.0);
            assert_eq!(end.transferred_mass_kg, start.mass_kg);
            near(end.extinction_time_s.unwrap(), lifetime, 1e-7*lifetime);
            near(end.thermostat_heat_j, 2e6*start.mass_kg, 1e-10);
            near(end.solid_energy_change_j+end.vapor_energy_j, end.thermostat_heat_j, 1e-9);
            assert!(end.unused_duration_s > 0.0);
        }
    }
    #[test]
    fn growth_and_subdivision_share_the_same_mass_owned_surface() {
        for g in shapes() {
            let m = g.at_extent(0.003).unwrap().mass_kg;
            let whole = run(g, m, 300.0, 0.5);
            let a = run(g, m, 300.0, 0.25);
            let b = run(g, a.final_geometry.mass_kg, 300.0, 0.25);
            assert!(whole.recession_m < 0.0 && whole.thermostat_heat_j < 0.0);
            near(whole.final_geometry.mass_kg, b.final_geometry.mass_kg, 1e-14*m);
            near(whole.thermostat_heat_j, a.thermostat_heat_j+b.thermostat_heat_j, 1e-8);
        }
    }
    fn coatings() -> [RetreatGeometry; 2] {
        [RetreatGeometry::spherical_coating(1000.0, 0.01).unwrap(),
            RetreatGeometry::cylindrical_coating(1000.0, 0.01, 0.1).unwrap()]
    }
    /// G0/G1: shell mass, exposed area and dA/dm use the same physical radius.
    #[test]
    fn curved_coatings_preserve_substrate_and_exact_shell_geometry() {
        for (i, g) in coatings().into_iter().enumerate() {
            let thickness = 0.002;
            let outer: f64 = 0.01 + thickness;
            let (volume, area) = if i == 0 {
                (4.0 * PI / 3.0 * (outer.powi(3) - 0.01_f64.powi(3)),
                    4.0 * PI * outer * outer)
            } else {
                (PI * 0.1 * (outer * outer - 0.01 * 0.01), 2.0 * PI * outer * 0.1)
            };
            let state = g.at_extent(thickness).unwrap();
            near(state.mass_kg, 1000.0 * volume, 2e-15);
            near(state.area_m2, area, 2e-17);
            let core = if i == 0 { 4.0 * PI / 3.0 * 0.01_f64.powi(3) }
                else { PI * 0.01 * 0.01 * 0.1 };
            near(g.occupied_volume_m3(state.mass_kg).unwrap(), core + volume, 2e-20);
            near(g.occupied_volume_m3(0.0).unwrap(), core, 2e-20);
            near(g.at_mass(state.mass_kg).unwrap().extent_m, thickness, 1e-17);
            let dm = state.mass_kg * 1e-5;
            let derivative = (g.at_mass(state.mass_kg + dm).unwrap().area_m2
                - g.at_mass(state.mass_kg - dm).unwrap().area_m2) / (2.0 * dm);
            near(derivative, state.area_mass_derivative_m2_kg.unwrap(), 2e-10);
            let grown = g.at_extent(0.003).unwrap();
            near(g.recession_m(state.mass_kg, grown.mass_kg).unwrap(), -0.001, 1e-17);
            assert!(grown.area_m2 > state.area_m2);
        }
    }
    /// G0/G3: thin films survive even when outer and core radii round equal.
    #[test]
    fn coating_mass_inversion_does_not_subtract_nearly_equal_radii() {
        for g in coatings() {
            for thickness in [1e-4, 1e-12, 1e-20] {
                let state = g.at_extent(thickness).unwrap();
                assert!(state.mass_kg > 0.0);
                let recovered = g.at_mass(state.mass_kg).unwrap();
                near(recovered.extent_m / thickness, 1.0, 2e-14);
                near(g.recession_m(state.mass_kg, 0.0).unwrap() / thickness, 1.0, 2e-14);
                assert!(g.recession_m(state.mass_kg, state.mass_kg.next_up()).unwrap() < 0.0);
            }
        }
    }
    /// G1/G3: integrate the moving outer surface, not a frozen footprint or core mass.
    #[test]
    fn deposition_grows_curved_coatings_and_releases_only_transferred_latent_heat() {
        let kinetics = SurfaceKinetics { accommodation_coefficient: 0.5,
            mass_transfer_coefficient_m_s: None };
        let reservoir = VaporReservoir::Prescribed { partial_pressure_pa: 300.0, temperature_k: 300.0 };
        let j = surface_flux(&Phase, 300.0, kinetics, reservoir).unwrap().mass_flux_kg_m2_s;
        for g in coatings() {
            let initial = g.at_extent(0.001).unwrap();
            let whole = run(g, initial.mass_kg, 300.0, 0.5);
            let expected = g.at_extent(0.001 - j / 1000.0 * 0.5).unwrap();
            near(whole.final_geometry.mass_kg, expected.mass_kg, 2e-15);
            assert!(whole.final_geometry.area_m2 > initial.area_m2);
            assert!(whole.transferred_mass_kg < 0.0);
            near(whole.thermostat_heat_j, 2e6 * whole.transferred_mass_kg, 1e-12);
            near(whole.energy_residual_j, 0.0, 1e-10);
            near(j * whole.area_time_m2_s, whole.transferred_mass_kg, 1e-16);
            let first = run(g, initial.mass_kg, 300.0, 0.25);
            let second = run(g, first.final_geometry.mass_kg, 300.0, 0.25);
            near(whole.final_geometry.mass_kg, second.final_geometry.mass_kg, 2e-15);
            near(whole.recession_m, first.recession_m + second.recession_m, 2e-17);
        }
    }
    /// G0/G1: reaching the substrate consumes the coating, never the inert core.
    #[test]
    fn coating_extinction_stops_at_the_substrate() {
        for g in coatings() {
            let initial = g.at_extent(1e-5).unwrap();
            let result = run(g, initial.mass_kg, 0.0, 100.0);
            assert_eq!(result.final_geometry.mass_kg, 0.0);
            assert_eq!(result.final_geometry.area_m2, 0.0);
            assert_eq!(result.transferred_mass_kg, initial.mass_kg);
            near(result.recession_m, initial.extent_m, 1e-18);
            assert!(result.extinction_time_s.unwrap() < 100.0);
            assert!(result.unused_duration_s > 0.0);
            near(result.thermostat_heat_j, 2e6 * initial.mass_kg, 1e-12);
        }
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(RetreatGeometry::spherical_coating(1000.0, radius).is_err());
            assert!(RetreatGeometry::cylindrical_coating(1000.0, radius, 0.1).is_err());
        }
        assert!(RetreatGeometry::cylindrical_coating(1000.0, 0.01, 0.0).is_err());
    }
    #[test]
    fn equilibrium_zero_time_and_invalid_inputs_are_explicit() {
        let g = shapes()[2];
        let m = g.at_extent(0.001).unwrap().mass_kg;
        let eq = run(g, m, 100.0, 2.0);
        assert_eq!(eq.initial, eq.final_geometry);
        assert_eq!(eq.thermostat_heat_j, 0.0);
        assert_eq!(run(g, m, 0.0, 0.0).transferred_mass_kg, 0.0);
        assert!(g.at_mass(-1.0).is_err());
        assert!(RetreatGeometry::layer(1000.0, 0.0).is_err());
        assert!(RetreatGeometry::sphere(f64::NAN).is_err());
        assert!(g.admit_density(900.0).is_err());
        let k = SurfaceKinetics { accommodation_coefficient: 1.0, mass_transfer_coefficient_m_s: None };
        for (mass, reservoir) in [(0.0, VaporReservoir::Prescribed { partial_pressure_pa: 0.0, temperature_k: 300.0 }),
            (m, VaporReservoir::Finite { vapor_mass_kg: 0.0, volume_m3: 1.0, temperature_k: 300.0 })] {
            assert!(advance_thermostatted_front(&Phase, g, mass, 300.0, 0.0, k, reservoir, 1.0).is_err());
        }
    }
}
