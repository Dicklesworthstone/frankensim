//! Real two-way wing/BEMT loads driving the ideal-gear six-DOF integrator.
//!
//! This consumes the current `OpenLoopDesign` aerodynamic geometry and closure,
//! with caller-supplied PHYSICAL mass/rotor inertias. It does not replace or
//! recalibrate a historical SimLoop preset. The present actuator-disk field is
//! axial along body +x by default; `new_spatial` selects physical disk-plane
//! overlap and vector wash for arbitrarily oriented, body-fixed shafts.
//! Motor torque, wind, gravity and canard setting are evaluated at each stage.
//! No stopped/reversed rotor or reversed axial inflow is silently floored.

use super::super::{Loads, TICK_HZ};
use super::geared::{self, GearedAssembly, GearedLoads, GearedState, GearedStep};
use super::{array, finite, vec};
use crate::Refusal;
use crate::aircraft::{OpenLoopDesign, wind_axis_lift};
use crate::propcoupling::spatial::{self, DiskAxes};
use crate::propcoupling::{
    CoupledStep, IncidentInflow, PropDisk, coupled_prop_airframe_step_with_incident,
};
use fs_ga::Vec3;
use fs_math::det;
use fs_wing::Panel;
use fs_wing::nonlinear::{InfluenceOperator, SectionClosure, StripSpec};

#[path = "propulsion/dynamic.rs"]
pub mod dynamic;

#[path = "propulsion/incident.rs"]
mod incident;
use incident::StageIncidentAir;

#[path = "propulsion/swirl.rs"]
pub mod swirl;

/// Blade handedness relative to a mounted rotor's declared positive spin axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PropellerHandedness {
    /// Positive absolute spin about the rotor axis uses the positive-pitch BEMT.
    Positive,
    /// Negative absolute spin about the rotor axis uses the mirrored blades.
    Negative,
}

impl PropellerHandedness {
    fn sign(self) -> f64 {
        match self {
            Self::Positive => 1.0,
            Self::Negative => -1.0,
        }
    }
}

/// Map one of the two aerodynamic disks to a physical geared rotor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropellerMount {
    /// Index into the gear assembly. Other geared rotors may represent an
    /// engine or flywheel and receive no propeller aerodynamic torque.
    pub rotor_index: usize,
    /// Blade handedness, independent of signed mechanical gear ratio.
    pub handedness: PropellerHandedness,
    /// Unit direction of positive BEMT thrust. `new` admits body `[1, 0, 0]`;
    /// `new_spatial` admits either direction along the physical rotor axis.
    /// No thrust direction is inferred from gear ratio or handedness.
    pub thrust_direction_body: [f64; 3],
}

/// Environmental and control values evaluated at one physical stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropulsionInputs {
    /// Positive air density [kg/m³].
    pub density_kg_m3: f64,
    /// Prescribed canard setting [rad], in the current open-loop domain.
    pub canard_rad: f64,
    /// Internal generalized motor torque, conjugate to common relative phase.
    pub motor_nm: f64,
    /// World-NED gravitational acceleration [m/s²], included exactly once.
    pub gravity_world_mps2: [f64; 3],
}

/// A real, converged propeller/airframe load evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct PropulsionEvaluation<const N: usize> {
    /// Carrier wrench and distinct aerodynamic shaft ports for the gear owner.
    pub loads: GearedLoads<N>,
    /// Aerodynamic body force before world gravity is added [N].
    pub aerodynamic_force_body_n: [f64; 3],
    /// Signed absolute physical spin rates for the L/R propeller rotors [rad/s].
    pub propeller_absolute_speed_rad_s: [f64; 2],
    /// Actual two-way wing/BEMT solution, including L/R thrust and shaft torque.
    pub coupled: CoupledStep,
}

/// Accepted motion and both real aerodynamic stage receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct PropulsionStep<const N: usize> {
    /// Accepted gear/body step, including physical work and momentum residuals.
    pub motion: GearedStep<N>,
    /// Start-of-step loads used only for the midpoint predictor.
    pub initial: PropulsionEvaluation<N>,
    /// Loads evaluated at the predicted midpoint and used by the accepted step.
    pub midpoint: PropulsionEvaluation<N>,
}

/// Two aerodynamic disks attached to an admitted, possibly larger gear train.
pub struct GearedPropulsion<'a, const N: usize> {
    design: &'a OpenLoopDesign,
    assembly: &'a GearedAssembly<N>,
    mounts: [PropellerMount; 2],
    spatial_axes: Option<DiskAxes>,
}

// Both wake models use exactly the same material-air sampling and body/shaft
// wrench assembly. Their distinct public results keep equilibrium residuals
// separate from physical wake acceleration.
struct AerodynamicStage<'a> {
    op: &'a InfluenceOperator,
    panels: &'a [Panel],
    strips: &'a [StripSpec],
    closure: SectionClosure<'a>,
    disks: &'a [PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    incident: IncidentInflow<'a>,
    carrier_angular_velocity_body_rad_s: [f64; 3],
    air: &'a StageIncidentAir<'a>,
    disk_wind_world_mps: [[f64; 3]; 2],
}

struct AerodynamicFields<'a> {
    gamma: &'a [f64],
    strip_alpha_eff_rad: &'a [f64],
    strip_velocity_mps: Option<&'a [[f64; 3]]>,
    thrust_n: [f64; 2],
    torque_nm: [f64; 2],
    oblique: Option<&'a [fs_airscrew::oblique::ObliqueLoads; 2]>,
}

trait AerodynamicSolution {
    fn load_fields(&self) -> AerodynamicFields<'_>;
}

impl<P> AerodynamicSolution for CoupledStep<P> {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        AerodynamicFields {
            gamma: &self.gamma,
            strip_alpha_eff_rad: &self.strip_alpha_eff_rad,
            strip_velocity_mps: self.strip_velocity_mps.as_deref(),
            thrust_n: self.thrust_n,
            torque_nm: self.torque_nm,
            oblique: None,
        }
    }
}

struct StageEvaluation<const N: usize, S> {
    loads: GearedLoads<N>,
    aerodynamic_force_body_n: [f64; 3],
    propeller_absolute_speed_rad_s: [f64; 2],
    aerodynamics: S,
}

fn reject(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use finite inputs inside the declared axial propeller/airframe domain".into(),
        ],
    }
}

fn wing_error(error: fs_wing::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

impl<'a, const N: usize> GearedPropulsion<'a, N> {
    /// Bind existing aerodynamic geometry to explicit geared mass properties.
    /// The two mounts must name distinct axial rotors aligned with body x.
    ///
    /// # Errors
    /// `propulsion-assembly-invalid` for inconsistent mass, mounts or geometry;
    /// existing rotor/drag-ledger admission refusals pass through.
    pub fn new(
        design: &'a OpenLoopDesign,
        assembly: &'a GearedAssembly<N>,
        mounts: [PropellerMount; 2],
    ) -> Result<Self, Refusal> {
        Self::build(design, assembly, mounts, None)
    }

    /// Bind canted, fixed propellers using actual three-dimensional disk
    /// overlap and wash vectors. Each thrust direction must be collinear
    /// with its mechanical spin axis; its sign is independent of handedness.
    /// The native BEMT still resolves only axial disk-average inflow.
    ///
    /// # Errors
    /// As [`Self::new`], plus `spatial-disk-invalid` for nonunit directions.
    pub fn new_spatial(
        design: &'a OpenLoopDesign,
        assembly: &'a GearedAssembly<N>,
        mounts: [PropellerMount; 2],
    ) -> Result<Self, Refusal> {
        let axes = DiskAxes::new(mounts.map(|mount| mount.thrust_direction_body))?;
        Self::build(design, assembly, mounts, Some(axes))
    }

    fn build(
        design: &'a OpenLoopDesign,
        assembly: &'a GearedAssembly<N>,
        mounts: [PropellerMount; 2],
        spatial_axes: Option<DiskAxes>,
    ) -> Result<Self, Refusal> {
        if design.gross_mass_kg != assembly.mass_kg()
            || !finite(vec(design.cg_m))
            || !design.drag_z_m.is_finite()
            || !design.camber_ratio.is_finite()
            || design.camber_ratio.abs() > 0.15
            || !design.cd0_wing.is_finite()
            || design.cd0_wing < 0.0
            || !design.oswald_e.is_finite()
            || design.oswald_e <= 0.0
            || mounts[0].rotor_index == mounts[1].rotor_index
        {
            return Err(reject(
                "propulsion-assembly-invalid",
                "mass/geometry invalid or duplicate propeller rotor",
            ));
        }
        for (mount, disk) in mounts.iter().zip(&design.disks) {
            let rotor = assembly.rotors().get(mount.rotor_index).ok_or_else(|| {
                reject(
                    "propulsion-assembly-invalid",
                    "propeller rotor index is outside the gear assembly",
                )
            })?;
            let aligned = if spatial_axes.is_some() {
                vec(mount.thrust_direction_body)
                    .cross(vec(rotor.axis_body))
                    .norm()
                    <= 1e-12
            } else {
                mount.thrust_direction_body == [1.0, 0.0, 0.0]
                    && rotor.axis_body[1] == 0.0
                    && rotor.axis_body[2] == 0.0
                    && rotor.axis_body[0].abs() == 1.0
            };
            if !aligned || !finite(vec(disk.center_m)) {
                return Err(reject(
                    "propulsion-assembly-invalid",
                    "thrust must align with a finite physical shaft; legacy wash requires body +x thrust",
                ));
            }
        }
        let spec = design.coupling;
        if ![spec.omega0, spec.clamp.0, spec.clamp.1, spec.tol]
            .iter()
            .all(|x| x.is_finite())
            || spec.omega0 <= 0.0
            || spec.omega0 > 1.0
            || spec.clamp.0 <= 0.0
            || spec.clamp.0 > spec.clamp.1
            || spec.clamp.1 > 1.0
            || spec.tol <= 0.0
        {
            return Err(reject(
                "propulsion-assembly-invalid",
                "invalid declared propeller coupling solver parameters",
            ));
        }
        design.rotor.admit().map_err(|error| Refusal {
            code: error.code,
            message: error.message,
            ranked_repairs: error.ranked_repairs,
        })?;
        design.ledger.admit()?;
        Ok(Self {
            design,
            assembly,
            mounts,
            spatial_axes,
        })
    }

    /// Evaluate the current distributed two-way wing/BEMT model at one state.
    /// `wind_world` returns material air velocity [m/s] at a WORLD position.
    /// Local incident velocity is `R^-1(V - wind) + omega × arm`; every
    /// strip and both disks receive all three rigid angular-rate components.
    /// Accepted local prop wash enters circulation, panel force, profile drag
    /// and section moment. Every force moment is taken about the assembly COM.
    ///
    /// The retained aerodynamic tier is quasi-steady: two-dimensional strip
    /// closures, axial disk inflow, straight wake, classical induced drag,
    /// and the existing parasite ledger. It does not supply a rudder model,
    /// disk harmonics/crossflow, added mass or a historical validity claim.
    ///
    /// # Errors
    /// Existing wing/BEMT/coupling refusals; `propulsion-state-invalid`;
    /// `propeller-spin-outside-domain` for stopped/opposite handedness spin.
    pub fn evaluate(
        &self,
        state: &GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<[f64; 2]>,
    ) -> Result<PropulsionEvaluation<N>, Refusal> {
        if warm.is_some_and(|w| !w.iter().all(|x| x.is_finite())) {
            return Err(reject(
                "propulsion-state-invalid",
                "nonfinite warm induction",
            ));
        }
        let result =
            self.evaluate_with(state, inputs, wind_world, |stage| match self.spatial_axes {
                Some(axes) => spatial::coupled_step(
                    &axes,
                    stage.op,
                    stage.panels,
                    stage.strips,
                    stage.closure,
                    &self.design.rotor,
                    stage.disks,
                    stage.freestream,
                    stage.rho,
                    &self.design.coupling,
                    warm,
                    stage.incident,
                ),
                None => coupled_prop_airframe_step_with_incident(
                    stage.op,
                    stage.panels,
                    stage.strips,
                    stage.closure,
                    &self.design.rotor,
                    stage.disks,
                    stage.freestream,
                    stage.rho,
                    &self.design.coupling,
                    warm,
                    stage.incident,
                ),
            })?;
        Ok(PropulsionEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            coupled: result.aerodynamics,
        })
    }

    fn evaluate_with<S: AerodynamicSolution>(
        &self,
        state: &GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        solve: impl FnOnce(AerodynamicStage<'_>) -> Result<S, Refusal>,
    ) -> Result<StageEvaluation<N, S>, Refusal> {
        let rotation = self.assembly.admit_state(state, 1.0 / TICK_HZ)?;
        if !inputs.density_kg_m3.is_finite()
            || inputs.density_kg_m3 <= 0.0
            || !inputs.canard_rad.is_finite()
            || inputs.canard_rad.abs() >= 0.6
            || !inputs.motor_nm.is_finite()
            || !finite(vec(inputs.gravity_world_mps2))
        {
            return Err(reject(
                "propulsion-state-invalid",
                "nonfinite or out-of-domain stage input",
            ));
        }
        let design = self.design;
        let speed = self.assembly.speed(state)?;
        let omega = vec(state.rigid.omega_body);
        let mut disks = design.disks;
        let mut propeller_absolute_speed_rad_s = [0.0; 2];
        for k in 0..2 {
            let mount = self.mounts[k];
            let rotor = self.assembly.rotors()[mount.rotor_index];
            let absolute = rotor.ratio * speed + vec(rotor.axis_body).dot(omega);
            let canonical = mount.handedness.sign() * absolute;
            if !canonical.is_finite() || canonical <= 0.0 {
                return Err(reject(
                    "propeller-spin-outside-domain",
                    "stopped or reversed blade rotation is outside the positive-speed BEMT branch",
                ));
            }
            propeller_absolute_speed_rad_s[k] = absolute;
            disks[k] = PropDisk {
                omega_rad_s: canonical,
                ..disks[k]
            };
        }
        let air = StageIncidentAir::new(
            design.cg_m,
            state.rigid.pos_m,
            state.rigid.vel_mps,
            state.rigid.omega_body,
            rotation,
            wind_world,
        );
        let freestream = air.sample(design.cg_m)?.relative_velocity_mps;
        let cg_air = vec(freestream);
        let speed2 = cg_air.dot(cg_air);
        let airspeed = speed2.sqrt();
        let alpha = det::atan2(freestream[2], freestream[0]);
        if !airspeed.is_finite() || airspeed <= 1.0 || freestream[0] <= 0.0 || alpha.abs() >= 0.6 {
            return Err(reject(
                "propulsion-state-invalid",
                "CG flow is outside the current forward-flight aerodynamic domain",
            ));
        }
        let rho = inputs.density_kg_m3;
        let (panels, strips) = design.layout(inputs.canard_rad)?;
        let op = InfluenceOperator::build(&panels, freestream, rho).map_err(wing_error)?;
        let camber = design.camber_ratio;
        // Same existing OpenLoopDesign attached/separated section closure.
        let closure = |_strip: usize, a: f64| crate::sections::section_lift(a, camber);
        let mut strip_velocity = Vec::with_capacity(strips.len());
        for strip in &strips {
            let panel = &panels[strip.panel_indices[0]];
            strip_velocity.push(
                air.sample(array((vec(panel.a) + vec(panel.b)).scale(0.5)))?
                    .relative_velocity_mps,
            );
        }
        let disk_air = [
            air.sample(disks[0].center_m)?,
            air.sample(disks[1].center_m)?,
        ];
        let incident = IncidentInflow {
            strip_velocity_mps: &strip_velocity,
            disk_velocity_mps: disk_air.each_ref().map(|air| air.relative_velocity_mps),
        };
        let aerodynamics = solve(AerodynamicStage {
            op: &op,
            panels: &panels,
            strips: &strips,
            closure: &closure,
            disks: &disks,
            freestream,
            rho,
            incident,
            carrier_angular_velocity_body_rad_s: state.rigid.omega_body,
            air: &air,
            disk_wind_world_mps: disk_air.map(|air| air.wind_world_mps),
        })?;
        let fields = aerodynamics.load_fields();
        let local = fields.strip_velocity_mps.ok_or_else(|| {
            reject(
                "propulsion-state-invalid",
                "distributed coupled solve omitted its accepted local velocities",
            )
        })?;
        let mut circulatory = fs_wing::wrench::CirculatoryWrenchAccumulator::new(rho, design.cg_m)
            .map_err(wing_error)?;
        for (s, strip) in strips.iter().enumerate() {
            for &j in &strip.panel_indices {
                circulatory
                    .add_panel(&panels[j], fields.gamma[j], local[s])
                    .map_err(wing_error)?;
            }
        }
        let circulatory_wrench = circulatory.total();
        let mut force = circulatory_wrench.force;
        let mut moment = circulatory_wrench.torque;
        // The induced-drag screen uses transverse lift, not its body-z
        // projection; the latter introduces a spurious cos(alpha)^2 factor.
        let lift = wind_axis_lift(force, alpha);
        let cm0 = -core::f64::consts::PI * camber;
        for (s, strip) in strips.iter().enumerate() {
            let panel = &panels[strip.panel_indices[0]];
            let v = local[s];
            let section_speed2 = v[0] * v[0] + v[2] * v[2];
            let section_speed = section_speed2.sqrt();
            let q = 0.5 * rho * section_speed2;
            let width = (panel.b[1] - panel.a[1]).abs();
            moment.y += q * strip.chord_m * strip.chord_m * width * cm0;
            let cd = crate::sections::section_profile_drag(
                fields.strip_alpha_eff_rad[s],
                design.cd0_wing,
            );
            let drag = q * strip.chord_m * width * cd;
            let profile_force = Vec3::new(
                -drag * v[0] / section_speed,
                0.0,
                -drag * v[2] / section_speed,
            );
            let arm = (vec(panel.a) + vec(panel.b)).scale(0.5) - vec(design.cg_m);
            force = force + profile_force;
            moment = moment + arm.cross(profile_force);
        }
        let area = 2.0 * design.wing_span_m * design.wing_chord_m;
        let aspect = design.wing_span_m * design.wing_span_m / (0.5 * area);
        let dynamic_pressure = 0.5 * rho * speed2;
        let cl = lift / (dynamic_pressure * area);
        let induced_drag =
            dynamic_pressure * area * cl * cl / (core::f64::consts::PI * aspect * design.oswald_e);
        let air_direction = cg_air.scale(1.0 / airspeed);
        force = force - air_direction.scale(induced_drag);
        let parasite_drag = design.ledger.evaluate(rho, airspeed)?.total_parasite_n;
        let parasite_force = air_direction.scale(-parasite_drag);
        force = force + parasite_force;
        moment =
            moment + Vec3::new(0.0, 0.0, design.drag_z_m - design.cg_m[2]).cross(parasite_force);
        let mut aerodynamic_nm = [0.0; N];
        for k in 0..2 {
            let mount = self.mounts[k];
            let propeller_force = match fields.oblique {
                Some(loads) => vec(loads[k].force_n),
                None => vec(mount.thrust_direction_body).scale(fields.thrust_n[k]),
            };
            force = force + propeller_force;
            moment = moment + (vec(disks[k].center_m) - vec(design.cg_m)).cross(propeller_force);
            if let Some(loads) = fields.oblique {
                // Transverse hub moments act on the carrier. The axial moment
                // belongs exclusively to the physical shaft port below.
                moment = moment + vec(loads[k].hub_moment_nm);
            }
            // BEMT reports required drive torque in its positive-speed branch.
            // Physical aerodynamic torque on the rotor has the opposite sign.
            aerodynamic_nm[mount.rotor_index] = -mount.handedness.sign() * fields.torque_nm[k];
        }
        let world_force = rotation.rotate(force).map_err(|_| {
            reject(
                "propulsion-state-invalid",
                "world aerodynamic force overflow",
            )
        })? + vec(inputs.gravity_world_mps2).scale(self.assembly.mass_kg());
        if !finite(force)
            || !finite(world_force)
            || !finite(moment)
            || !aerodynamic_nm.iter().all(|x| x.is_finite())
        {
            return Err(reject(
                "propulsion-state-invalid",
                "aerodynamic wrench overflow",
            ));
        }
        Ok(StageEvaluation {
            loads: GearedLoads {
                carrier: Loads {
                    force_n: array(world_force),
                    moment_nm: array(moment),
                },
                motor_nm: inputs.motor_nm,
                aerodynamic_nm,
            },
            aerodynamic_force_body_n: array(force),
            propeller_absolute_speed_rad_s,
            aerodynamics,
        })
    }

    /// Advance the actual coupled aerodynamic model using start/midpoint
    /// evaluations. Inputs receive time and full gear/body state; wind receives
    /// time and WORLD position. The second solve uses the first accepted wash
    /// as its warm start, never a held endpoint force. Returned midpoint wash
    /// is suitable as the next call's warm predictor.
    ///
    /// # Errors
    /// All model/integrator/closure refusals pass through. No partial motion is
    /// published; caller-owned callback side effects are outside that guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &GearedState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        mut inputs: C,
        wind_world: W,
    ) -> Result<PropulsionStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut first = None;
        let mut second = None;
        let mut stage_warm = warm;
        let motion = geared::step(
            self.assembly,
            state,
            t_s,
            dt_s,
            angular_origin_world,
            |time, stage| {
                let evaluation = self.evaluate(
                    stage,
                    inputs(time, stage)?,
                    &|point| wind_world(time, point),
                    stage_warm,
                )?;
                stage_warm = Some(evaluation.coupled.w_slip);
                let loads = evaluation.loads;
                if first.is_none() {
                    first = Some(evaluation);
                } else {
                    second = Some(evaluation);
                }
                Ok(loads)
            },
        )?;
        Ok(PropulsionStep {
            motion,
            initial: first
                .ok_or_else(|| reject("propulsion-state-invalid", "missing initial stage"))?,
            midpoint: second
                .ok_or_else(|| reject("propulsion-state-invalid", "missing midpoint stage"))?,
        })
    }
}

#[cfg(test)]
mod lift_projection_tests {
    use super::{Vec3, det, wind_axis_lift};

    #[test]
    fn signed_lift_projection_is_independent_of_pitch_and_streamwise_force() {
        for alpha in [-0.5_f64, 0.0, 0.5] {
            let (s, c) = (det::sin(alpha), det::cos(alpha));
            for lift in [-250.0, 0.0, 250.0] {
                for drag in [0.0, 37.0] {
                    // Construct known wind-frame loads in body coordinates.
                    // The side force is orthogonal to this lift axis too.
                    let force = Vec3::new(lift * s - drag * c, 19.0, -lift * c - drag * s);
                    assert!((wind_axis_lift(force, alpha) - lift).abs() < 1e-12);
                }
            }
        }
    }
}
