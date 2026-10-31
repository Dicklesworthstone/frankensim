//! Hydrostatic pressure in a stable stack of source-bound fluid layers.
//!
//! At downward coordinate d, density is rho_top + sum(delta_rho_i H(d-d_i)).
//! Thus displaced mass and its first moment are the same positive combination
//! of whole-envelope and plane-clipped moments. This avoids subtracting nearly
//! equal cut volumes for a thin layer. Pressure potential is
//! |g| [rho_top integral d dV + sum(delta_rho_i integral (d-d_i)+ dV)].
//! Its negative pose derivative gives the same force and torque, including
//! passage through interfaces. No pressure jump or surface tension is inferred.

use super::{
    DisplacedVolume, DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticStep,
    HydrostaticStepControls, HydrostaticStepError, finite, finite_vec, gravity_direction, poll,
};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_ga::{Se3, Vec3};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{Gravity, RigidBodyState};

/// Maximum horizontal interfaces in one frozen hydrostatic reservoir.
pub const MAX_HYDROSTATIC_LAYERS: usize = 32;

/// A fluid below one horizontal interface. Interfaces are supplied in strictly
/// increasing downward depth from the common world datum plane [m].
#[derive(Clone, Copy)]
pub struct HydrostaticLayer<'a> {
    /// Interface depth measured parallel to gravity, not a world z coordinate.
    pub depth_m: f64,
    /// Original material owner of the fluid below this interface.
    pub source: &'a dyn FluidTransport,
}

struct ResolvedLayer<'a> {
    input: HydrostaticLayer<'a>,
    transport: FluidTransportState,
    density_jump: f64,
}

/// Frozen, stably stratified hydrostatic reservoir with exact clipped moments.
/// An optional upper fluid fills the space above the first interface, allowing
/// atmospheric buoyancy as well as gas-over-liquid and liquid-over-liquid
/// immersion. With no upper fluid, the first interface is a free surface.
///
/// Each layer has constant density resolved once from its original material;
/// temperature, pressure and species changes require rebuilding this medium.
/// Density must not decrease downward. Unstable stratification, compressible
/// density profiles, entrainment, waves, drag and added mass are not modeled.
pub struct LayeredHydrostaticMedium<'a> {
    above: Option<(&'a dyn FluidTransport, FluidTransportState)>,
    layers: Vec<ResolvedLayer<'a>>,
    gravity: Gravity,
    gravity_magnitude: f64,
    down_world: Vec3,
    datum_world: Vec3,
}

impl<'a> LayeredHydrostaticMedium<'a> {
    /// Resolve every supplied owner once and admit a stable stack. An empty
    /// interface list is allowed only with an upper fluid (uniform immersion).
    /// The datum sets the pressure-potential zero; uniform exterior pressure
    /// cancels on a closed displacement envelope and need not be supplied.
    pub fn new(
        above: Option<&'a dyn FluidTransport>,
        layers: &[HydrostaticLayer<'a>],
        gravity: Gravity,
        datum_world: Vec3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        finite_vec(datum_world)?;
        if layers.len() > MAX_HYDROSTATIC_LAYERS || (above.is_none() && layers.is_empty()) {
            return Err(HydrostaticError::Invalid(
                "hydrostatic layer budget exceeded or empty reservoir",
            ));
        }
        let (down_world, gravity_magnitude) = gravity_direction(gravity)?;
        let above = above
            .map(|source| source.transport_state().map(|state| (source, state)))
            .transpose()?;
        let mut density = above.map_or(0.0, |(_, state)| state.density_kg_m3());
        let mut last_depth = f64::NEG_INFINITY;
        let mut resolved = Vec::with_capacity(layers.len());
        for &input in layers {
            poll(&mut cancelled)?;
            if !input.depth_m.is_finite() || input.depth_m <= last_depth {
                return Err(HydrostaticError::Invalid(
                    "hydrostatic interface depths must be finite and strictly increasing",
                ));
            }
            let transport = input.source.transport_state()?;
            if transport.density_kg_m3() < density {
                return Err(HydrostaticError::Invalid(
                    "hydrostatic density must not decrease downward",
                ));
            }
            let density_jump = transport.density_kg_m3() - density;
            density = transport.density_kg_m3();
            last_depth = input.depth_m;
            resolved.push(ResolvedLayer {
                input,
                transport,
                density_jump,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            above,
            layers: resolved,
            gravity,
            gravity_magnitude,
            down_world,
            datum_world,
        })
    }

    /// Body weight is supplied separately and exactly once through this gravity.
    pub const fn gravity(&self) -> Gravity {
        self.gravity
    }
    /// Original upper-fluid source and the snapshot actually used in mechanics.
    pub fn above(&self) -> Option<(&'a dyn FluidTransport, FluidTransportState)> {
        self.above
    }
    /// Interfaces, original sources and consumed snapshots in supplied order.
    pub fn layers(
        &self,
    ) -> impl ExactSizeIterator<Item = (HydrostaticLayer<'a>, FluidTransportState)> + '_ {
        self.layers
            .iter()
            .map(|layer| (layer.input, layer.transport))
    }

    /// Buoyancy force and moment from the same analytic geometry as the
    /// single-fluid path. The geometric centroid in `displaced()` describes all
    /// immersed volume; `center_of_buoyancy_world()` describes displaced MASS.
    /// These centers generally differ in a layered reservoir.
    ///
    /// O(cells × nonzero density jumps); uniform ambient terms use cached whole
    /// envelope moments in O(1), with no allocations during load evaluation.
    pub fn load(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, &mut cancelled)
    }

    fn load_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        poll(cancelled)?;
        let pose: Se3 = body.world_from_body(state)?;
        let mut mass = 0.0;
        let mut center = None;
        let mut pressure_potential = 0.0;
        let mut displaced: Option<DisplacedVolume> = None;
        if let Some((_, transport)) = self.above {
            let whole = mesh.whole(pose, self.down_world, self.datum_world)?;
            add_moments(
                whole,
                transport.density_kg_m3(),
                &mut mass,
                &mut center,
                &mut pressure_potential,
            )?;
            displaced = Some(whole);
        }
        for layer in &self.layers {
            poll(cancelled)?;
            if layer.density_jump == 0.0 {
                continue;
            }
            let cut = mesh.submerged_below(
                pose,
                self.down_world,
                self.datum_world,
                layer.input.depth_m,
                cancelled,
            )?;
            if displaced.is_none() {
                displaced = Some(cut);
            }
            add_moments(
                cut,
                layer.density_jump,
                &mut mass,
                &mut center,
                &mut pressure_potential,
            )?;
        }
        let displaced =
            displaced.ok_or(HydrostaticError::Invalid("empty hydrostatic moment sum"))?;
        let load = HydrostaticLoad::from_mass_moments(
            displaced,
            mass,
            center,
            finite(pressure_potential * self.gravity_magnitude)?,
            self.gravity,
            body,
            state,
        )?;
        poll(cancelled)?;
        Ok(load)
    }

    /// Advance the native rigid-body dynamics with re-evaluated layered
    /// buoyancy, using the same second-order endpoint-load step and absolute energy
    /// budget as `HydrostaticMedium::step`. No fluid damping is invented.
    pub fn step(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticStep, HydrostaticStepError> {
        super::step::step_with_load(
            self.gravity,
            body.frame(),
            state,
            controls,
            &mut cancelled,
            |state, cancel| self.load_frame(mesh, body.frame(), state, cancel),
        )
    }
}

impl super::HydrostaticReservoir for LayeredHydrostaticMedium<'_> {
    fn gravity(&self) -> Gravity {
        LayeredHydrostaticMedium::gravity(self)
    }

    fn load_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        LayeredHydrostaticMedium::load_frame(self, mesh, body, state, &mut cancelled)
    }
}

fn add_moments(
    displaced: DisplacedVolume,
    density: f64,
    mass: &mut f64,
    center: &mut Option<Vec3>,
    potential: &mut f64,
) -> Result<(), HydrostaticError> {
    if displaced.volume_m3() == 0.0 {
        return Ok(());
    }
    let contribution = finite(density * displaced.volume_m3())?;
    if contribution <= 0.0 {
        return Err(HydrostaticError::Invalid(
            "displaced fluid mass is unrepresentable",
        ));
    }
    let next_mass = finite(*mass + contribution)?;
    let point = displaced
        .centroid_body_m()
        .ok_or(HydrostaticError::Invalid("missing displaced centroid"))?;
    *center = Some(match *center {
        None => point,
        Some(old) => finite_vec(old + (point - old).scale(contribution / next_mass))?,
    });
    *mass = next_mass;
    *potential = finite(*potential + density * displaced.depth_integral_m4)?;
    Ok(())
}
