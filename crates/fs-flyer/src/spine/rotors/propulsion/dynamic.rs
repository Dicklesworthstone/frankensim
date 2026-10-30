//! Real shaft/body dynamics coupled to two evolving uniform axial wakes.
//!
//! Each disk's induced velocity is physical memory. A stage solves the wing
//! under that wash, projects its induced field into each actual disk axis,
//! and evaluates native blade forces at the stored induction. It does not
//! re-equilibrate the wake or delay a steady BEMT thrust. The shared carrier
//! wrench owner applies thrust arms and external aerodynamic shaft torques
//! exactly once, including counterrotation and canted fixed mounts.
//!
//! This tier uses `fs_airscrew::dynamic_inflow`'s constant-density apparent
//! air mass. Its impulse/work balance is an AXIAL reduced-state balance,
//! not the world momentum of a resolved wake. `new_oblique` additionally
//! resolves revolution-averaged section crossflow and transverse hub loads;
//! the wake still uses the declared uniform AXIAL transport approximation.
//! Skew-induced wake distributions, wake rotation during large maneuvers,
//! blade-phase harmonics and unsteady airfoil memory remain absent.

use super::{
    AerodynamicFields, AerodynamicSolution, AerodynamicStage, DiskAxes, GearedLoads,
    GearedPropulsion, GearedState, GearedStep, PropulsionInputs, Refusal, geared, reject,
    wing_error,
};
use fs_airscrew::dynamic_inflow::{InflowStep, RotorEvaluation, UniformAxialInflow};
use fs_airscrew::oblique::{AzimuthQuadrature, DiskFrame, ObliqueLoads, SpinDirection};
use fs_wing::induced_velocity_free;
use fs_wing::nonlinear::solve_nonlinear_with_incident;

/// All physical memory of the dynamic propulsion step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicState {
    /// Full six-DOF carrier plus canonical shared shaft coordinate/momentum.
    pub motion: GearedState,
    /// Actual L/R uniform disk-induced velocities [m/s], not solver guesses.
    pub induced_mps: [f64; 2],
}

/// Instantaneous wing and propeller loads at an evolving wake state.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicAerodynamics {
    /// Physical induction used by both wash and blade sections [m/s].
    pub induced_mps: [f64; 2],
    /// Projected incident plus wing-induced disk velocity, excluding self-wash.
    pub disk_axial_mps: [f64; 2],
    /// L/R native blade thrust, shaft torque and wake acceleration.
    pub propellers: [RotorEvaluation; 2],
    /// Full section-integrated body-frame loads when oblique blades are selected.
    /// Their transverse hub moments exclude the separately owned shaft torque.
    pub oblique: Option<[ObliqueLoads; 2]>,
    /// Wing circulation at this physical wash state [m²/s].
    pub gamma: Vec<f64>,
    /// Incident plus vector prop wash used for circulation and panel forces.
    pub strip_velocity_mps: Vec<[f64; 3]>,
    /// Native wing lift [N].
    pub wing_lift_n: f64,
    /// Actual nonlinear wing residual; this is not a wake-equilibrium defect.
    pub wing_residual: f64,
}

impl AerodynamicSolution for DynamicAerodynamics {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        AerodynamicFields {
            gamma: &self.gamma,
            strip_velocity_mps: Some(&self.strip_velocity_mps),
            thrust_n: self.propellers.map(|prop| prop.thrust_n),
            torque_nm: self.propellers.map(|prop| prop.torque_nm),
            oblique: self.oblique.as_ref(),
        }
    }
}

/// One actual force evaluation; physical wake rates need not vanish.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicEvaluation<const N: usize> {
    /// Carrier wrench and signed shaft ports, with internal motor torque separate.
    pub loads: GearedLoads<N>,
    /// Aerodynamic force in body coordinates, excluding gravity [N].
    pub aerodynamic_force_body_n: [f64; 3],
    /// Actual signed rotor spin rates [rad/s].
    pub propeller_absolute_speed_rad_s: [f64; 2],
    /// Native wing/blade loads and physical wake rates.
    pub aerodynamics: DynamicAerodynamics,
}

/// Accepted body/shaft/wake state and both real aerodynamic evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicStep<const N: usize> {
    /// Complete state for the next tick, including both evolving wakes.
    pub state: DynamicState,
    /// Native geared mechanics with world linear/angular impulse accounting.
    pub motion: GearedStep<N>,
    /// Per-disk axial fluid impulse and induced-work balances.
    pub inflow: [InflowStep; 2],
    /// Start loads used to predict the joint midpoint state.
    pub initial: DynamicEvaluation<N>,
    /// Loads at the predicted body/shaft/wake midpoint, used for acceptance.
    pub midpoint: DynamicEvaluation<N>,
}

/// Explicit dynamic-wake model choice over an existing physical gear assembly.
pub struct DynamicPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p GearedPropulsion<'a, N>,
    inflow: UniformAxialInflow,
    axes: DiskAxes,
    oblique: Option<(AzimuthQuadrature, [DiskFrame; 2])>,
}

fn air_error(error: fs_airscrew::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

impl<'p, 'a, const N: usize> DynamicPropulsion<'p, 'a, N> {
    /// Bind constant air density and uniform wake inertia. Both disks use
    /// actual three-dimensional overlap and their declared thrust axes.
    /// Physical rotor/aircraft inertias remain those of the gear assembly.
    ///
    /// # Errors
    /// Native dynamic-inflow or disk-axis admission.
    pub fn new(
        propulsion: &'p GearedPropulsion<'a, N>,
        density_kg_m3: f64,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            propulsion,
            inflow: UniformAxialInflow::new(propulsion.design.rotor.radius_m, density_kg_m3)
                .map_err(air_error)?,
            axes: DiskAxes::new(propulsion.mounts.map(|mount| mount.thrust_direction_body))?,
            oblique: None,
        })
    }

    /// Select revolution-averaged oblique blade loads and physical tilt-rate
    /// kinematics. The same scalar axial wake and straight axial wing wash are
    /// retained deliberately: this is not a skewed-wake momentum closure.
    /// The finite azimuth rule uses a fixed body-frame radial reference.
    ///
    /// # Errors
    /// Native density, axis, quadrature and orthonormal-frame admission.
    pub fn new_oblique(
        propulsion: &'p GearedPropulsion<'a, N>,
        density_kg_m3: f64,
        azimuth_points: usize,
    ) -> Result<Self, Refusal> {
        let mut result = Self::new(propulsion, density_kg_m3)?;
        let rule = AzimuthQuadrature::new(azimuth_points).map_err(air_error)?;
        let frame = |k: usize| {
            let mount = propulsion.mounts[k];
            let axis = mount.thrust_direction_body;
            // The least-aligned coordinate yields a well-conditioned material
            // reference, fixed in the carrier rather than chosen from airflow.
            let reference = (0..3)
                .min_by(|&i, &j| axis[i].abs().total_cmp(&axis[j].abs()))
                .expect("three coordinates");
            let mut radial = axis.map(|x| -axis[reference] * x);
            radial[reference] += 1.0;
            let length =
                (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            radial = radial.map(|x| x / length);
            let rotor_axis = propulsion.assembly.rotors()[mount.rotor_index].axis_body;
            let alignment =
                axis[0] * rotor_axis[0] + axis[1] * rotor_axis[1] + axis[2] * rotor_axis[2];
            let spin = if mount.handedness.sign() * alignment > 0.0 {
                SpinDirection::AlongThrust
            } else {
                SpinDirection::AgainstThrust
            };
            DiskFrame::new(axis, radial, spin).map_err(air_error)
        };
        result.oblique = Some((rule, [frame(0)?, frame(1)?]));
        Ok(result)
    }

    /// The explicit per-disk apparent-air-mass model.
    #[must_use]
    pub fn inflow_model(&self) -> UniformAxialInflow {
        self.inflow
    }

    fn solve(
        &self,
        stage: AerodynamicStage<'_>,
        induced_mps: [f64; 2],
    ) -> Result<DynamicAerodynamics, Refusal> {
        let wash = self.axes.wash_map(
            stage.strips,
            stage.panels,
            stage.disks,
            self.propulsion.design.rotor.radius_m,
        )?;
        let strip_velocity_mps: Vec<_> = stage
            .incident
            .strip_velocity_mps
            .iter()
            .zip(wash)
            .map(|(&incident, factors)| self.axes.velocity(incident, factors, induced_mps))
            .collect();
        let wing = solve_nonlinear_with_incident(
            stage.op,
            stage.panels,
            stage.strips,
            stage.freestream,
            stage.rho,
            stage.closure,
            &strip_velocity_mps,
            None,
            None,
        )
        .map_err(wing_error)?;
        let mut disk_velocity_mps = [[0.0; 3]; 2];
        let disk_axial_mps = core::array::from_fn(|k| {
            let induced = induced_velocity_free(
                stage.disks[k].center_m,
                stage.panels,
                &wing.gamma,
                stage.freestream,
            );
            disk_velocity_mps[k] =
                core::array::from_fn(|j| stage.incident.disk_velocity_mps[k][j] + induced[j]);
            self.axes
                .axial(k, stage.incident.disk_velocity_mps[k], induced)
        });
        let oblique = if let Some((rule, frames)) = &self.oblique {
            let loads = |k: usize| {
                rule.evaluate(
                    &self.propulsion.design.rotor,
                    stage.rho,
                    frames[k],
                    disk_velocity_mps[k],
                    stage.carrier_angular_velocity_body_rad_s,
                    induced_mps[k],
                    stage.disks[k].omega_rad_s,
                )
                .map_err(air_error)
            };
            Some([loads(0)?, loads(1)?])
        } else {
            None
        };
        let propeller = |k: usize| {
            if let Some(loads) = &oblique {
                return Ok(RotorEvaluation {
                    thrust_n: loads[k].thrust_n,
                    torque_nm: loads[k].torque_nm,
                    inflow: self
                        .inflow
                        .rate(disk_axial_mps[k], induced_mps[k], loads[k].thrust_n)
                        .map_err(air_error)?,
                });
            }
            self.inflow
                .evaluate(
                    &self.propulsion.design.rotor,
                    disk_axial_mps[k],
                    induced_mps[k],
                    stage.disks[k].omega_rad_s,
                )
                .map_err(air_error)
        };
        Ok(DynamicAerodynamics {
            induced_mps,
            disk_axial_mps,
            propellers: [propeller(0)?, propeller(1)?],
            oblique,
            gamma: wing.gamma,
            strip_velocity_mps,
            wing_lift_n: wing.total_lift_n,
            wing_residual: wing.residual,
        })
    }

    /// Evaluate wing, blades and all shaft/body forces at the actual stored
    /// induction. Density must equal the admitted constant wake density.
    ///
    /// # Errors
    /// Native force/wake refusals, or `dynamic-propulsion-invalid` for a
    /// nonfinite wake state or changing density.
    pub fn evaluate(
        &self,
        state: &DynamicState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<DynamicEvaluation<N>, Refusal> {
        if inputs.density_kg_m3 != self.inflow.density_kg_m3()
            || !state.induced_mps.iter().all(|x| x.is_finite())
        {
            return Err(reject(
                "dynamic-propulsion-invalid",
                "wake state must be finite and density must remain fixed",
            ));
        }
        let result = self
            .propulsion
            .evaluate_with(&state.motion, inputs, wind_world, |stage| {
                self.solve(stage, state.induced_mps)
            })?;
        Ok(DynamicEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            aerodynamics: result.aerodynamics,
        })
    }

    fn advance_wakes(
        &self,
        before: [f64; 2],
        stage: &DynamicAerodynamics,
        dt_s: f64,
    ) -> Result<[InflowStep; 2], Refusal> {
        let step = |k: usize| {
            self.inflow
                .step_held(
                    before[k],
                    stage.disk_axial_mps[k],
                    stage.propellers[k].thrust_n,
                    dt_s,
                )
                .map_err(air_error)
        };
        Ok([step(0)?, step(1)?])
    }

    /// Advance the joint body, shaft and wake with two real load evaluations.
    /// The initial loads predict ALL states to half time. Midpoint blade
    /// forces drive the accepted gear map and the implicit wake-flux step.
    /// This is a second-order partitioned midpoint method for smooth admitted
    /// cases. It is not an unsteady-wing model or a total fluid/body energy
    /// theorem. Each scalar wake storage balance and the native rigid/shaft
    /// impulse balance remain separately inspectable.
    ///
    /// # Errors
    /// All native force/integrator refusals. Caller state is immutable; even
    /// a late wake refusal publishes no partial body/shaft/wake update.
    /// Callback side effects remain the caller's responsibility.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &DynamicState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
    ) -> Result<DynamicStep<N>, Refusal>
    where
        C: FnMut(f64, &DynamicState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut initial = None;
        let mut midpoint = None;
        let mut stage_induction = state.induced_mps;
        let motion = geared::step(
            self.propulsion.assembly,
            &state.motion,
            t_s,
            dt_s,
            angular_origin_world,
            |time, mechanical| {
                let stage = DynamicState {
                    motion: *mechanical,
                    induced_mps: stage_induction,
                };
                let evaluation = self.evaluate(&stage, inputs(time, &stage)?, &|point| {
                    wind_world(time, point)
                })?;
                let loads = evaluation.loads;
                if initial.is_none() {
                    stage_induction = self
                        .advance_wakes(state.induced_mps, &evaluation.aerodynamics, 0.5 * dt_s)?
                        .map(|wake| wake.induced_mps);
                    initial = Some(evaluation);
                } else {
                    midpoint = Some(evaluation);
                }
                Ok(loads)
            },
        )?;
        let initial =
            initial.ok_or_else(|| reject("dynamic-propulsion-invalid", "missing initial stage"))?;
        let midpoint = midpoint
            .ok_or_else(|| reject("dynamic-propulsion-invalid", "missing midpoint stage"))?;
        let inflow = self.advance_wakes(state.induced_mps, &midpoint.aerodynamics, dt_s)?;
        Ok(DynamicStep {
            state: DynamicState {
                motion: motion.state,
                induced_mps: inflow.map(|wake| wake.induced_mps),
            },
            motion,
            inflow,
            initial,
            midpoint,
        })
    }
}
