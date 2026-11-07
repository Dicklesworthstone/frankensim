//! Simultaneous free translation and rotation of finite porous rigid support.
//!
//! Six coupled impulses drive the original rigid Cayley integrator and the
//! original nonlinear Darcy--Forchheimer fluid step. There is no translational
//! bearing or rotational guide. The actual candidate midpoint COM defines both
//! the moving drag support and the origin of every measured native exchange.
//! The skeleton velocity is `v_mid + (R0 * omega_body_mid) cross (x - COM_mid)`;
//! these are the native rigid work-conjugate velocities, including for arbitrary
//! full inertia and gyroscopic motion.
//!
//! Homogeneous sphere/ellipsoid support and principal resistance axes follow the
//! original rigid body at its geodesic midpoint attitude. Full dry inertia and
//! gyroscopic motion remain independent of the geometric and resistance axes.
//! Sampling is cell-centred inside stationary box walls. Fluid remains represented
//! throughout the support: no excluded volume, displaced mass, fictitious
//! interior inertia or porosity is inferred.
//! Native force-frame preparation preserves the incoming physical fluid momentum
//! when the support moves. A trial holds its midpoint support and coefficients
//! during the native endpoint force quadrature. Spatial/time refinement is still
//! necessary, especially when cells enter or leave the support; a discontinuous
//! mask can refuse the bounded nonlinear solve. This is not a sharp impermeable
//! moving wall, subcell overlap integration, or a deforming material mask.

use super::rotor::{PorousRotorControls, PorousRotorResistance};
use super::*;
use fs_lbm::d3q19::BoundaryLink3;
use fs_mbd::StepReceipt;

#[path = "free_spherical/support.rs"]
mod support;
pub use support::{PorousBodySupport, PrincipalPorousResistance};
use support::{PreparedSupport, midpoint_orientation};

#[path = "free_spherical/balance.rs"]
mod balance;
#[path = "free_spherical/solve.rs"]
mod solve;

/// Original fluid and finite dry inertia with all six rigid freedoms retained.
pub struct FreePorousBody {
    fluid: BoundaryGrid3,
    body: RigidBodyState,
    properties: MassProperties,
    units: PorousLatticeUnits3,
    radius_m: f64,
    support: PorousBodySupport,
    links: Vec<BoundaryLink3>,
    cells: usize,
    time_s: f64,
    accepted_steps: u64,
}

/// Accepted native fluid/body states and their independent physical accounts.
pub struct FreePorousStep {
    pub motion: FreePorousBody,
    /// Original lattice receipt, with all moments about the actual midpoint COM.
    pub exchange: ResistanceMomentumExchange3,
    pub exchange_origin_m: Vec3,
    pub mechanics: StepReceipt,
    pub nonlinear: NewtonReport,
    pub native_trials: usize,
    pub trial_cells: usize,
    pub midpoint_velocity_m_s: Vec3,
    /// Native Cayley work velocity `R0 * I^-1 (L0 + L1)/2` in world coordinates.
    pub angular_velocity_world_rad_s: Vec3,
    pub balance_origin_m: Vec3,
    pub fluid_linear_change_kg_m_s: Vec3,
    pub fluid_angular_change_kg_m2_s: Vec3,
    pub body_linear_change_kg_m_s: Vec3,
    pub body_spin_change_kg_m2_s: Vec3,
    /// Actual endpoint orbital plus spin change about `balance_origin_m`.
    pub body_angular_change_kg_m2_s: Vec3,
    pub body_external_impulse_kg_m_s: Vec3,
    /// COM force lever arm plus the authored independent world couple.
    pub body_external_angular_impulse_kg_m2_s: Vec3,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    /// Sum of participant residual norms: opposing defects cannot cancel.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub fluid_mass_change_kg: f64,
    pub fluid_kinetic_change_j: f64,
    pub body_kinetic_change_j: f64,
    pub body_external_work_j: f64,
    pub fluid_external_work_j: f64,
    pub fluid_drag_work_j: f64,
    pub relative_dissipation_j: f64,
    /// Native fluid transport/pressure/collision remainder, not heat or a
    /// certificate for complete fluid-energy or trajectory accuracy.
    pub other_fluid_kinetic_change_j: f64,
    pub mechanical_work_residual_j: f64,
}

struct Trial {
    mechanics: StepReceipt,
    centre: Vec3,
    velocity: Vec3,
    omega_body: Vec3,
    omega_world: Vec3,
    linear_change: Vec3,
    spin_change: Vec3,
    exchange: ResistanceMomentumExchange3,
    residual: Vec<f64>,
}

struct Step<'a> {
    motion: &'a FreePorousBody,
    resistance: PrincipalPorousResistance,
    forcing: Vec<[f64; 3]>,
    controls: PorousRotorControls,
    units: Units,
    external_linear: Vec3,
    external_couple: Vec3,
    impulse_scales: [f64; 2],
    time_s: f64,
    accepted_steps: u64,
}

impl FreePorousBody {
    /// Bind a spherical porous support to original populations and dry inertia.
    /// This retains the existing sphere entry point and native admission.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fluid: BoundaryGrid3,
        body: RigidBodyState,
        properties: MassProperties,
        units: PorousLatticeUnits3,
        radius_m: f64,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PorousMomentumError> {
        Self::new_with_support(
            fluid,
            body,
            properties,
            units,
            PorousBodySupport::Sphere { radius_m },
            maximum_cells,
            gate,
        )
    }

    /// Bind sphere or body-aligned ellipsoid support to the SAME native owners.
    ///
    /// The geometry is centred at the original dry body's COM. Its dimensions,
    /// mass and full inertia are held fixed; no density, interior-fluid mass,
    /// excluded volume or inertia tensor is inferred from the mask. The complete
    /// enclosing sphere must fit inside six stationary walls, so all intervening
    /// rotations and the straight COM segment stay inside the native domain.
    /// Voxel solids and open/periodic/moving box walls remain unsupported.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_support(
        fluid: BoundaryGrid3,
        body: RigidBodyState,
        properties: MassProperties,
        units: PorousLatticeUnits3,
        support: PorousBodySupport,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PorousMomentumError> {
        poll(gate)?;
        units.admit()?;
        body.kinetic_energy(properties)?;
        let radius_m = support.admit()?;
        contain_sphere(
            body.pose().position_world(),
            radius_m,
            units,
            fluid.dimensions(),
        )?;
        let (cells, links) = closed_domain(&fluid, maximum_cells, gate)?;
        Ok(Self {
            fluid,
            body,
            properties,
            units,
            radius_m,
            support,
            links,
            cells,
            time_s: 0.0,
            accepted_steps: 0,
        })
    }

    pub const fn fluid(&self) -> &BoundaryGrid3 {
        &self.fluid
    }
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    pub const fn units(&self) -> PorousLatticeUnits3 {
        self.units
    }
    /// Radius of the enclosing sphere; exactly the authored radius for a sphere.
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Original fixed body-frame support, carried by the accepted rigid pose.
    pub const fn support(&self) -> PorousBodySupport {
        self.support
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Actual endpoint angular velocity, resolved in world coordinates.
    pub fn angular_velocity_world_rad_s(&self) -> Result<Vec3, PorousMomentumError> {
        angular_velocity(self.body, self.properties)
    }

    /// Advance one native lattice tick with simultaneous force and torque recoil.
    /// `resistance` supplies homogeneous isotropic native lattice coefficients
    /// inside the candidate midpoint support. The independent `force(cell)` is
    /// sampled once per box cell and held through all trials. The body receives
    /// held `external_force_n` at its moving COM and `external_couple_nm` in world
    /// coordinates. Their impulses are distinct from the measured fluid drag.
    ///
    /// Every six-coordinate Newton trial executes the original rigid and fluid
    /// owners from the same accepted states. The final replay supplies the actual
    /// returned populations, force field, body momentum and attitude. One scratch
    /// grid is reused; native trial/cell ceilings include probes and final replay.
    ///
    /// `controls.flow.maximum_speed_m_s` bounds fluid speed and the conservative
    /// enclosing-sphere bound `|v| + radius * |omega|` at both rigid endpoints and the
    /// work midpoint. `maximum_rotation_rad` bounds the native Cayley increment.
    /// Complete start/end enclosing spheres must fit in the closed box, which also
    /// contains their straight midpoint translation and any intervening rotation.
    /// No accepted state changes on refusal.
    pub fn step_candidate(
        &self,
        resistance: PorousRotorResistance,
        force: impl FnMut([usize; 3]) -> [f64; 3],
        external_force_n: Vec3,
        external_couple_nm: Vec3,
        controls: PorousRotorControls,
        gate: &CancelGate,
    ) -> Result<FreePorousStep, PorousMomentumError> {
        self.step_candidate_with_principal_resistance(
            PrincipalPorousResistance::isotropic(resistance.darcy, resistance.forchheimer),
            force,
            external_force_n,
            external_couple_nm,
            controls,
            gate,
        )
    }

    /// Advance with anisotropic Darcy resistance attached to the body axes.
    ///
    /// The three nonnegative principal resistances and isotropic Forchheimer
    /// coefficient use the original lattice force-density convention. Zero
    /// principal values admit free directions; infinite values are not a no-slip
    /// constraint. Every native trial rotates both the actual ellipsoid support
    /// and the resistance axes through the geodesic midpoint of its two rigid
    /// attitudes. The work-conjugate velocity is still the original Cayley
    /// midpoint velocity, not an independently differenced or prescribed spin.
    ///
    /// Rotation/translation-dependent mask changes are resolved by the SAME
    /// bounded six-impulse solve, including all probes and its final replay.
    /// Discontinuous cell-centre support can require a smaller tick; no frozen
    /// mask, post-hoc momentum repair, or extra fluid solve is substituted.
    pub fn step_candidate_with_principal_resistance(
        &self,
        resistance: PrincipalPorousResistance,
        mut force: impl FnMut([usize; 3]) -> [f64; 3],
        external_force_n: Vec3,
        external_couple_nm: Vec3,
        controls: PorousRotorControls,
        gate: &CancelGate,
    ) -> Result<FreePorousStep, PorousMomentumError> {
        poll(gate)?;
        let units = self.units.admit()?;
        let c = controls.flow;
        c.admit(units, self.cells)?;
        resistance.admit()?;
        if !external_force_n.is_finite()
            || !external_couple_nm.is_finite()
            || !controls.maximum_rotation_rad.is_finite()
            || controls.maximum_rotation_rad <= 0.0
        {
            return Err(invalid(
                "finite body loads and positive rotation bound required",
            ));
        }
        admit_speed(
            self.body.center_of_mass_velocity_world(self.properties)?,
            self.angular_velocity_world_rad_s()?,
            self.radius_m,
            c.maximum_speed_m_s,
        )?;
        let dt = self.units.timestep_s;
        let external_linear = external_force_n.scale(dt);
        let external_couple = external_couple_nm.scale(dt);
        let time_s = finite(self.time_s + dt)?;
        let accepted_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or(invalid("step counter overflow"))?;
        if time_s <= self.time_s || !external_linear.is_finite() || !external_couple.is_finite() {
            return Err(invalid(
                "free porous clock or external impulse is not representable",
            ));
        }
        let largest_inertia = a(self.properties.principal_inertia_body())
            .into_iter()
            .fold(0.0_f64, f64::max);
        let impulse_scales = [
            2.0 * self.properties.mass() * c.velocity_scale_m_s,
            2.0 * largest_inertia * (c.velocity_scale_m_s / self.radius_m),
        ];
        if impulse_scales.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid("free porous impulse scales are not representable"));
        }
        let _ = fluid_change(&self.fluid, &self.fluid, self.units, c, gate)?;
        let [nx, ny, nz] = self.fluid.dimensions();
        let mut forcing = Vec::with_capacity(self.cells);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x % 256 == 0 {
                        poll(gate)?;
                    }
                    let value = force([x, y, z]);
                    if value.iter().any(|x| !x.is_finite()) {
                        return Err(invalid("independent fluid force must be finite"));
                    }
                    forcing.push(value);
                }
            }
        }
        Step {
            motion: self,
            resistance,
            forcing,
            controls,
            units,
            external_linear,
            external_couple,
            impulse_scales,
            time_s,
            accepted_steps,
        }
        .solve(gate)
    }
}

fn angular_velocity(
    body: RigidBodyState,
    properties: MassProperties,
) -> Result<Vec3, PorousMomentumError> {
    Ok(body.pose().orientation().rotate_body_to_world(
        properties.angular_velocity_body_checked(body.angular_momentum_body())?,
    ))
}

fn contain_sphere(
    centre: Vec3,
    radius: f64,
    units: PorousLatticeUnits3,
    dimensions: [usize; 3],
) -> Result<(), PorousMomentumError> {
    for i in 0..3 {
        let distance = dd(a(centre)[i]) - dd(a(units.grid_origin_m)[i]);
        let extent = dd(dimensions[i] as f64) * dd(units.cell_size_m);
        let low = finite((distance - dd(radius)).to_f64())?;
        let high = finite((extent - distance - dd(radius)).to_f64())?;
        if low < 0.0 || high < 0.0 {
            return Err(invalid(
                "complete moving porous sphere must stay inside the box",
            ));
        }
    }
    Ok(())
}

fn admit_speed(
    velocity: Vec3,
    omega: Vec3,
    radius: f64,
    maximum: f64,
) -> Result<(), PorousMomentumError> {
    let bound = norm(velocity) + radius * norm(omega);
    if !bound.is_finite() || bound > maximum {
        Err(invalid(
            "free porous translation plus rim speed exceeds the declared domain",
        ))
    } else {
        Ok(())
    }
}

fn average(x: Vec3, y: Vec3) -> Vec3 {
    v(std::array::from_fn(|i| {
        (dd(0.5) * dd(a(x)[i]) + dd(0.5) * dd(a(y)[i])).to_f64()
    }))
}

fn difference(x: Vec3, y: Vec3) -> Vec3 {
    v(std::array::from_fn(|i| {
        (dd(a(x)[i]) - dd(a(y)[i])).to_f64()
    }))
}

fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

/// Existing sphere construction and stepping names remain source-compatible.
pub type FreeSphericalPorousBody = FreePorousBody;
/// Existing spherical report name for the original common rigid/fluid receipt.
pub type FreeSphericalPorousStep = FreePorousStep;
