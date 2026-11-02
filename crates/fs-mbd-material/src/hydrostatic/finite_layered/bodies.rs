//! Shared immiscible inventories and density-weighted native body motion.
use super::FiniteLayeredHydrostaticTank;
use crate::hydrostatic::finite_tank::step_bodies_with;
use crate::hydrostatic::{
    DisplacedVolume, FiniteTankBodiesLoad, FiniteTankBodiesStep, FiniteTankBody,
    HydrostaticBuoyancy, HydrostaticError, HydrostaticStepControls, HydrostaticStepError, finite,
    finite_vec, poll,
};
use fs_ga::Vec3;
use fs_mbd::RigidBodyState;

/// One common interface and each body's displacement below that interface.
/// Records are bottom-to-top; the last interface is the free surface.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteTankBodiesInterface {
    /// Interface height above the fixed tank floor [m].
    pub surface_height_m: f64,
    /// Corresponding gravity-normal interface point in world metres.
    pub surface_point_world_m: Vec3,
    /// Original-frame cumulative cuts, in the supplied body order.
    /// These are not the separate volumes of the adjacent liquid layer.
    pub displaced_below: Vec<DisplacedVolume>,
    /// A*h - sum(displaced volume) - cumulative liquid volume [m3].
    /// The individual layer defect is the difference of adjacent residuals.
    pub cumulative_volume_defect_m3: f64,
    /// Actual bounded geometry evaluations for this cumulative inventory.
    pub volume_evaluations: usize,
}

impl FiniteLayeredHydrostaticTank<'_> {
    /// Resolve every conserved interface shared by all original body poses.
    ///
    /// Geometry admission and body-to-tank transforms occur once; all roots
    /// reuse the same bounded cuts buffer. The returned geometry is at most
    /// 32 interfaces by 64 bodies, within the original total-cell budget.
    /// Numeric admission does not establish envelope nonoverlap or source
    /// immiscibility; the original finite-tank caller obligations still apply.
    pub fn interfaces_at_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FiniteTankBodiesInterface>, HydrostaticError> {
        self.interfaces_and_potential_at_bodies(bodies, states, &mut cancelled)
            .map(|(interfaces, _)| interfaces)
    }

    /// Reuse each accepted interface cut for its complete stack potential.
    /// Source changes need both quantities but no individual wrench assembly.
    pub(super) fn interfaces_and_potential_at_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<FiniteTankBodiesInterface>, f64), HydrostaticError> {
        let top = self.layers.last().expect("nonempty liquid stack");
        let mut prepared = top.cumulative.prepare_bodies(bodies, states, cancelled)?;
        let mut interfaces = Vec::with_capacity(self.layers.len());
        let mut above_height = f64::INFINITY;
        let mut potential = 0.0;
        for layer in self.layers.iter().rev() {
            let surface = layer
                .cumulative
                .fill_prepared_bodies(&mut prepared, cancelled)?;
            check_interface(surface.surface_height_m, above_height)?;
            above_height = surface.surface_height_m;
            let weight = layer.density_jump / layer.cumulative.transport().density_kg_m3();
            potential = finite(potential + weight * surface.liquid_gravitational_energy_change_j)?;
            interfaces.push(FiniteTankBodiesInterface {
                surface_height_m: surface.surface_height_m,
                surface_point_world_m: surface.surface_point_world_m,
                displaced_below: prepared.cuts.clone(),
                cumulative_volume_defect_m3: surface.liquid_volume_defect_m3,
                volume_evaluations: surface.volume_evaluations,
            });
        }
        interfaces.reverse();
        poll(cancelled)?;
        Ok((interfaces, potential))
    }

    /// Density-weighted loads on all bodies in the SAME finite liquid stack.
    ///
    /// Each cumulative inventory solves A*h_i - sum(D_body(h_i)) = C_i.
    /// Positive density jumps combine original-frame displaced mass moments
    /// and the native finite-bath potentials. The square of TOTAL displacement
    /// retains cross-body coupling at each interface; fluid energy appears
    /// once in the returned surface account. Subtracting nearly equal thin
    /// layer moments or independent per-body tank energies is unnecessary.
    ///
    /// Source snapshots and individual liquid inventories stay fixed. Body
    /// transforms/bounds are prepared once and cuts are reused across roots;
    /// wrenches are built only after the complete mass-moment reduction.
    /// Uniform stacks delegate exactly to the original total-inventory tank.
    /// Intersecting bodies, collisions, waves, mixing and fluid inertia are
    /// outside this quasistatic owner. Failure returns no partial load set.
    pub fn load_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteTankBodiesLoad, HydrostaticError> {
        poll(&mut cancelled)?;
        let top = self.layers.last().expect("nonempty liquid stack");
        if self.layers[..self.layers.len() - 1]
            .iter()
            .all(|layer| layer.density_jump == 0.0)
        {
            return top.cumulative.load_bodies(bodies, states, cancelled);
        }
        let mut prepared = top
            .cumulative
            .prepare_bodies(bodies, states, &mut cancelled)?;
        let mut moments = vec![MassMoments::dry(); bodies.len()];
        let mut top_surface = None;
        let mut potential = 0.0;
        let mut above_height = f64::INFINITY;
        for layer in self.layers.iter().rev() {
            poll(&mut cancelled)?;
            if layer.density_jump == 0.0 {
                continue;
            }
            let surface = layer
                .cumulative
                .fill_prepared_bodies(&mut prepared, &mut cancelled)?;
            check_interface(surface.surface_height_m, above_height)?;
            above_height = surface.surface_height_m;
            let first = top_surface.is_none();
            if first {
                top_surface = Some(surface);
            }
            for (moments, &cut) in moments.iter_mut().zip(&prepared.cuts) {
                poll(&mut cancelled)?;
                if first {
                    moments.displaced = cut;
                }
                moments.add(cut, layer.density_jump)?;
            }
            let weight = layer.density_jump / layer.cumulative.transport().density_kg_m3();
            potential = finite(potential + weight * surface.liquid_gravitational_energy_change_j)?;
        }
        let mut surface = top_surface.expect("positive top density");
        surface.liquid_gravitational_energy_change_j = potential;
        let mut loads = Vec::with_capacity(bodies.len());
        for ((body, &state), moments) in bodies.iter().zip(states).zip(moments) {
            poll(&mut cancelled)?;
            loads.push(HydrostaticBuoyancy::from_mass_moments(
                moments.displaced,
                moments.mass,
                moments.center,
                self.gravity(),
                body.frame,
                state,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(FiniteTankBodiesLoad {
            surface,
            bodies: loads,
        })
    }

    /// Advance all original rigid bodies using the same evolving interfaces.
    /// The uniform-tank batch owner performs the native staged Verlet motion
    /// and one total-energy gate, with this stack's complete potential at both
    /// endpoints. All states are proposed together; no liquid/body history is
    /// mutated on overflow, cancellation or numerical-energy refusal.
    pub fn step_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteTankBodiesStep, HydrostaticStepError> {
        step_bodies_with(
            self.gravity(),
            bodies,
            states,
            controls,
            &mut cancelled,
            |states, cancel| self.load_bodies(bodies, states, cancel),
        )
    }
}

#[derive(Clone, Copy)]
struct MassMoments {
    displaced: DisplacedVolume,
    mass: f64,
    center: Option<Vec3>,
}
impl MassMoments {
    fn dry() -> Self {
        Self {
            displaced: DisplacedVolume::dry(),
            mass: 0.0,
            center: None,
        }
    }
    fn add(&mut self, cut: DisplacedVolume, density_jump: f64) -> Result<(), HydrostaticError> {
        if cut.volume_m3() == 0.0 {
            return Ok(());
        }
        let contribution = finite(density_jump * cut.volume_m3())?;
        if contribution <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "displaced layered-tank mass is unrepresentable",
            ));
        }
        let mass = finite(self.mass + contribution)?;
        let point = cut.centroid_body_m().ok_or(HydrostaticError::Invalid(
            "missing displaced layered-tank centroid",
        ))?;
        self.center = Some(match self.center {
            None => point,
            Some(old) => finite_vec(old + (point - old).scale(contribution / mass))?,
        });
        self.mass = mass;
        Ok(())
    }
}

fn check_interface(height: f64, above_height: f64) -> Result<(), HydrostaticError> {
    if height >= above_height {
        Err(HydrostaticError::Invalid(
            "finite liquid interfaces are unresolved at the volume tolerance",
        ))
    } else {
        Ok(())
    }
}
