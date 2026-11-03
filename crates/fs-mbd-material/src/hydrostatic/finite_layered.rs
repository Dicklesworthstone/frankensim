//! Conserved immiscible layers through the original finite-tank volume solve.
//!
//! Bottom-to-top cumulative inventories C_i satisfy A*h_i - D(h_i) = C_i.
//! With rho_(N+1)=0, positive density jumps rho_i-rho_(i+1) combine the
//! displaced moments and finite-bath potentials at these moving interfaces.
//! This is the gravitational potential of the actual conserved layer stack.

use super::{
    DisplacedVolume, DisplacementMesh, FiniteHydrostaticTank, FiniteTankControls, HydrostaticError,
    HydrostaticLoad, HydrostaticReservoir, HydrostaticStep, HydrostaticStepControls,
    HydrostaticStepError, HydrostaticTankGeometry, MAX_HYDROSTATIC_LAYERS, finite, finite_vec,
    poll,
};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_ga::{Se3, Vec3};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{Gravity, RigidBodyState};

mod bodies;
mod density_change;
mod miscibility;
pub use bodies::FiniteTankBodiesInterface;
pub use density_change::{FiniteLayerDensityChange, FiniteLayeredTankDensityChange};

/// One conserved immiscible liquid, supplied in bottom-to-top order.
#[derive(Clone, Copy)]
pub struct FiniteHydrostaticLayer<'a> {
    /// Original owner of the homogeneous single-phase liquid state.
    pub source: &'a dyn FluidTransport,
    /// Positive liquid volume, excluding the displacement envelope [m3].
    pub liquid_volume_m3: f64,
}

struct ResolvedFiniteLayer<'a> {
    input: FiniteHydrostaticLayer<'a>,
    // This private cumulative-volume solve retains the actual source snapshot.
    // Its fictitious cumulative mass is never reported as a layer's inventory.
    cumulative: FiniteHydrostaticTank<'a, dyn FluidTransport + 'a>,
    // Closed material inventory, never reconstructed after a density update.
    liquid_mass_kg: f64,
    density_jump: f64,
}

/// One solved interface, in the same bottom-to-top order as the liquid inputs.
/// The last interface is the free surface; the upper gas density is neglected.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankInterface {
    /// Height above the fixed tank floor [m].
    pub surface_height_m: f64,
    /// World-space point on this gravity-normal interface [m].
    pub surface_point_world_m: Vec3,
    /// Entire envelope below this interface, not just the adjacent layer.
    pub displaced_below: DisplacedVolume,
    /// Signed A*h - V_sub - cumulative liquid volume residual [m3].
    /// An individual layer's residual is the difference of its two adjacent
    /// cumulative residuals, bounded by twice the supplied volume tolerance.
    pub cumulative_volume_defect_m3: f64,
    /// Number of original bounded geometry evaluations for this interface.
    pub volume_evaluations: usize,
}

/// Finite, stably ordered immiscible liquids with conserved individual masses.
///
/// Each source is resolved once. Density is constant within each liquid and
/// must not increase upward. The same rectangular, fixed, gravity-aligned tank
/// and displacement-envelope admission as `FiniteHydrostaticTank` apply.
/// Every interface follows its own cumulative liquid inventory as the body
/// moves. Sources, volumes and masses remain unchanged during native motion.
///
/// This is quasistatic hydrostatics with negligible upper-gas density. No
/// mixing, entrainment, diffusion, capillarity, phase change, fluid inertia,
/// sloshing or moving-tank reaction is inferred. Source validity and physical
/// immiscibility remain caller obligations; the transport trait cannot prove
/// composition or interface stability. It is not a compressible density law.
pub struct FiniteLayeredHydrostaticTank<'a> {
    layers: Vec<ResolvedFiniteLayer<'a>>,
    liquid_mass_kg: f64,
}

impl<'a> FiniteLayeredHydrostaticTank<'a> {
    /// Freeze 1..=32 source states and their positive liquid inventories.
    /// The original volume controls apply separately to each cumulative solve.
    /// Unrepresentable sums and unstable density ordering refuse before return.
    pub fn new(
        layers: &[FiniteHydrostaticLayer<'a>],
        gravity: Gravity,
        geometry: HydrostaticTankGeometry,
        controls: FiniteTankControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        if layers.is_empty() || layers.len() > MAX_HYDROSTATIC_LAYERS {
            return Err(HydrostaticError::Invalid(
                "finite layered tank requires a bounded nonempty liquid stack",
            ));
        }
        let mut resolved = Vec::with_capacity(layers.len());
        let mut volume = 0.0;
        let mut mass = 0.0;
        let mut previous_density = f64::INFINITY;
        for &input in layers {
            poll(&mut cancelled)?;
            if !input.liquid_volume_m3.is_finite() || input.liquid_volume_m3 <= 0.0 {
                return Err(HydrostaticError::Invalid(
                    "finite liquid layers require positive finite inventories",
                ));
            }
            let cumulative_volume = finite(volume + input.liquid_volume_m3)?;
            if cumulative_volume <= volume {
                return Err(HydrostaticError::Invalid(
                    "finite liquid layer volume is lost in the cumulative inventory",
                ));
            }
            let cumulative = FiniteHydrostaticTank::new(
                input.source,
                gravity,
                geometry,
                cumulative_volume,
                controls,
            )?;
            let density = cumulative.transport().density_kg_m3();
            if density > previous_density {
                return Err(HydrostaticError::Invalid(
                    "finite liquid density must not increase upward",
                ));
            }
            let layer_mass = finite(density * input.liquid_volume_m3)?;
            if layer_mass <= 0.0 {
                return Err(HydrostaticError::Invalid(
                    "finite liquid layer mass is unrepresentable",
                ));
            }
            mass = finite(mass + layer_mass)?;
            volume = cumulative_volume;
            previous_density = density;
            resolved.push(ResolvedFiniteLayer {
                input,
                cumulative,
                liquid_mass_kg: layer_mass,
                density_jump: 0.0,
            });
        }
        set_density_jumps(&mut resolved);
        poll(&mut cancelled)?;
        Ok(Self {
            layers: resolved,
            liquid_mass_kg: mass,
        })
    }

    /// Current source bindings and checked snapshots, in bottom-to-top order.
    /// Use `layer_masses_kg` for exact conserved inventories: after a source
    /// update, multiplying the new density and volume can round differently.
    pub fn layers(
        &self,
    ) -> impl ExactSizeIterator<Item = (FiniteHydrostaticLayer<'a>, FluidTransportState)> + '_ {
        self.layers
            .iter()
            .map(|layer| (layer.input, layer.cumulative.transport()))
    }

    /// Exact closed masses, in the same bottom-to-top order as `layers` [kg].
    pub fn layer_masses_kg(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.layers.iter().map(|layer| layer.liquid_mass_kg)
    }

    /// Sum of the individual conserved source-derived liquid masses [kg].
    pub const fn liquid_mass_kg(&self) -> f64 {
        self.liquid_mass_kg
    }

    /// Total conserved liquid volume, excluding the solid envelope [m3].
    pub fn liquid_volume_m3(&self) -> f64 {
        self.layers
            .last()
            .expect("nonempty liquid stack")
            .cumulative
            .liquid_volume_m3()
    }

    /// The same fixed tank used by all interface solves.
    pub fn geometry(&self) -> HydrostaticTankGeometry {
        self.layers[0].cumulative.geometry()
    }

    /// Body gravity supplied separately and exactly once to native mechanics.
    pub fn gravity(&self) -> Gravity {
        self.layers[0].cumulative.gravity()
    }

    /// Resolve all actual interface positions at the original body-frame pose.
    ///
    /// This explicit geometry query allocates at most 32 interface records.
    /// Load and step evaluation remain allocation-free and skip zero-density
    /// jumps. Interfaces indistinguishable at the configured solve tolerance
    /// refuse rather than publishing a crossed or zero-thickness layer.
    pub fn interfaces_at_pose(
        &self,
        mesh: &DisplacementMesh,
        world_from_body: Se3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FiniteTankInterface>, HydrostaticError> {
        self.interfaces_and_potential_at_pose(mesh, world_from_body, &mut cancelled)
            .map(|(interfaces, _)| interfaces)
    }

    fn interfaces_and_potential_at_pose(
        &self,
        mesh: &DisplacementMesh,
        world_from_body: Se3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<FiniteTankInterface>, f64), HydrostaticError> {
        poll(cancelled)?;
        let mut interfaces = Vec::with_capacity(self.layers.len());
        let mut above_height = f64::INFINITY;
        let mut potential = 0.0;
        // Same positive contributions and reduction order as native loading.
        for layer in self.layers.iter().rev() {
            let fill = layer
                .cumulative
                .fill_at_pose(mesh, world_from_body, &mut *cancelled)?;
            if fill.surface_height_m >= above_height {
                return Err(HydrostaticError::Invalid(
                    "finite liquid interfaces are unresolved at the volume tolerance",
                ));
            }
            above_height = fill.surface_height_m;
            let weight = layer.density_jump / layer.cumulative.transport().density_kg_m3();
            potential = finite(potential + weight * fill.liquid_gravitational_energy_change_j)?;
            interfaces.push(FiniteTankInterface {
                surface_height_m: fill.surface_height_m,
                surface_point_world_m: fill.surface_point_world_m,
                displaced_below: fill.displaced,
                cumulative_volume_defect_m3: fill.liquid_volume_defect_m3,
                volume_evaluations: fill.volume_evaluations,
            });
        }
        interfaces.reverse();
        poll(cancelled)?;
        Ok((interfaces, potential))
    }

    /// Density-weighted buoyancy and torque, excluding body weight.
    /// The returned potential is the complete liquid stack's gravitational
    /// energy change from the same inventories with no displacement envelope.
    pub fn load(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, cancelled)
    }

    /// Continue the original endpoint-load Verlet mechanics, recomputing
    /// interface heights at both endpoints with the same conserved liquids.
    pub fn step(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: HydrostaticStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticStep, HydrostaticStepError> {
        self.step_frame(mesh, body.frame(), state, controls, cancelled)
    }
}

fn set_density_jumps(layers: &mut [ResolvedFiniteLayer<'_>]) {
    for index in 0..layers.len() {
        let above = layers
            .get(index + 1)
            .map_or(0.0, |layer| layer.cumulative.transport().density_kg_m3());
        layers[index].density_jump = layers[index].cumulative.transport().density_kg_m3() - above;
    }
}

impl HydrostaticReservoir for FiniteLayeredHydrostaticTank<'_> {
    fn gravity(&self) -> Gravity {
        FiniteLayeredHydrostaticTank::gravity(self)
    }

    fn load_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        poll(&mut cancelled)?;
        let top = self.layers.last().expect("nonempty liquid stack");
        if self.layers[..self.layers.len() - 1]
            .iter()
            .all(|layer| layer.density_jump == 0.0)
        {
            // A uniform stack is exactly the original total-inventory tank.
            return top.cumulative.load_frame(mesh, body, state, cancelled);
        }
        let pose = body.world_from_body(state)?;
        let mut displaced = None;
        let mut mass = 0.0;
        let mut center = None;
        let mut potential = 0.0;
        let mut above_height = f64::INFINITY;
        for layer in self.layers.iter().rev() {
            poll(&mut cancelled)?;
            if layer.density_jump == 0.0 {
                continue;
            }
            let fill = layer.cumulative.fill_at_pose(mesh, pose, &mut cancelled)?;
            if fill.surface_height_m >= above_height {
                return Err(HydrostaticError::Invalid(
                    "finite liquid interfaces are unresolved at the volume tolerance",
                ));
            }
            above_height = fill.surface_height_m;
            if displaced.is_none() {
                displaced = Some(fill.displaced);
            }
            let volume = fill.displaced.volume_m3();
            if volume > 0.0 {
                let contribution = finite(layer.density_jump * volume)?;
                if contribution <= 0.0 {
                    return Err(HydrostaticError::Invalid(
                        "displaced finite liquid mass is unrepresentable",
                    ));
                }
                let next_mass = finite(mass + contribution)?;
                let point = fill
                    .displaced
                    .centroid_body_m()
                    .ok_or(HydrostaticError::Invalid(
                        "missing displaced finite liquid centroid",
                    ))?;
                center = Some(match center {
                    None => point,
                    Some(old) => finite_vec(old + (point - old).scale(contribution / next_mass))?,
                });
                mass = next_mass;
            }
            // Each original finite-bath potential already contains its moving
            // waterline correction. Only the positive density jump is applied;
            // no nearly equal per-layer cut moments or energies are subtracted.
            let weight = layer.density_jump / layer.cumulative.transport().density_kg_m3();
            potential = finite(potential + weight * fill.liquid_gravitational_energy_change_j)?;
        }
        let load = HydrostaticLoad::from_mass_moments(
            displaced.expect("positive top density"),
            mass,
            center,
            potential,
            self.gravity(),
            body,
            state,
        )?;
        poll(&mut cancelled)?;
        Ok(load)
    }
}
