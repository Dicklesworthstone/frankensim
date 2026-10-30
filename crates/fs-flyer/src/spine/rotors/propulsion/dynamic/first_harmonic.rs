//! Spatially asymmetric disk induction coupled to aircraft and shaft momentum.
//!
//! The same finite modal field reaches every native blade quadrature point and
//! every wing strip inside the existing straight-wash window. Blade thrust and
//! transverse moments drive mean and first-harmonic air inertia. Axial shaft
//! torque remains owned once by the original geared mechanics. The wake uses
//! an explicitly selected zero-skew NORMAL transport approximation even when
//! the section law resolves incident crossflow. It is not a full skew wake.

use super::{
    AerodynamicFields, AerodynamicSolution, AerodynamicStage, DynamicPropulsion, GearedLoads,
    GearedPropulsion, GearedState, GearedStep, PropulsionInputs, Refusal, air_error, geared,
    reject,
};
use fs_airscrew::dynamic_inflow::first_harmonic::{
    FirstHarmonicInflow, HarmonicInduction, HarmonicInflowRate, HarmonicInflowStep,
};
use fs_airscrew::oblique::{DiskFrame, ObliqueLoads};

/// Full mechanical state plus each disk's three physical induction velocities.
/// First modes use the fixed carrier disk frames returned by `disk_frames`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicState {
    pub motion: GearedState,
    pub induction: [HarmonicInduction; 2],
}

/// The actual blade/wing field and conjugate wake loads at one force stage.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicAerodynamics {
    pub induction: [HarmonicInduction; 2],
    pub disk_axial_mps: [f64; 2],
    /// Carrier rate about each disk normal; never physical blade speed.
    pub frame_rate_rad_s: [f64; 2],
    pub propellers: [ObliqueLoads; 2],
    pub wake_force_n: [[f64; 3]; 2],
    pub wake_rate: [HarmonicInflowRate; 2],
    pub gamma: Vec<f64>,
    pub strip_velocity_mps: Vec<[f64; 3]>,
    pub wing_lift_n: f64,
    pub wing_residual: f64,
}

impl AerodynamicSolution for HarmonicAerodynamics {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        AerodynamicFields {
            gamma: &self.gamma,
            strip_velocity_mps: Some(&self.strip_velocity_mps),
            thrust_n: self.propellers.map(|prop| prop.thrust_n),
            torque_nm: self.propellers.map(|prop| prop.torque_nm),
            oblique: Some(&self.propellers),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicEvaluation<const N: usize> {
    pub loads: GearedLoads<N>,
    pub aerodynamic_force_body_n: [f64; 3],
    pub propeller_absolute_speed_rad_s: [f64; 2],
    pub aerodynamics: HarmonicAerodynamics,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicStep<const N: usize> {
    pub state: HarmonicState,
    pub motion: GearedStep<N>,
    pub inflow: [HarmonicInflowStep; 2],
    pub initial: HarmonicEvaluation<N>,
    pub midpoint: HarmonicEvaluation<N>,
    /// Implicit wake input work minus actual predicted-midpoint blade
    /// induced work. This partitioned quadrature defect is NOT dissipation.
    pub coupling_work_residual_j: [f64; 2],
}

/// Explicit normal-transport first-harmonic choice. Constant-density modal
/// masses and all quadrature trigonometry are prepared once. The original
/// physical aircraft, rotor inertias, section law and wing solver are reused.
pub struct HarmonicPropulsion<'p, 'a, const N: usize> {
    base: DynamicPropulsion<'p, 'a, N>,
    inflow: FirstHarmonicInflow,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}

impl<'p, 'a, const N: usize> HarmonicPropulsion<'p, 'a, N> {
    pub fn new(
        propulsion: &'p GearedPropulsion<'a, N>,
        density_kg_m3: f64,
        azimuth_points: usize,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            base: DynamicPropulsion::new_oblique(propulsion, density_kg_m3, azimuth_points)?,
            inflow: FirstHarmonicInflow::new(propulsion.design.rotor.radius_m, density_kg_m3)
                .map_err(air_error)?,
        })
    }

    #[must_use]
    pub fn disk_frames(&self) -> [DiskFrame; 2] {
        self.base
            .oblique
            .as_ref()
            .expect("prepared harmonic frames")
            .1
    }

    #[must_use]
    pub fn inflow_model(&self) -> FirstHarmonicInflow {
        self.inflow
    }

    fn solve(
        &self,
        stage: AerodynamicStage<'_>,
        induction: [HarmonicInduction; 2],
    ) -> Result<HarmonicAerodynamics, Refusal> {
        let wing = self
            .base
            .solve_wing(&stage, induction.map(|u| u.mean_mps), Some(&induction))?;
        let (rule, frames) = self
            .base
            .oblique
            .as_ref()
            .expect("prepared harmonic frames");
        let blade = |k: usize| {
            rule.evaluate_first_harmonic(
                &self.base.propulsion.design.rotor,
                stage.rho,
                frames[k],
                wing.disk_velocity_mps[k],
                stage.carrier_angular_velocity_body_rad_s,
                induction[k],
                stage.disks[k].omega_rad_s,
            )
            .map_err(air_error)
        };
        let propellers = [blade(0)?, blade(1)?];
        let force = |k: usize| {
            frames[k]
                .wake_force_n(&propellers[k], self.base.propulsion.design.rotor.radius_m)
                .map_err(air_error)
        };
        let wake_force_n = [force(0)?, force(1)?];
        let frame_rate_rad_s = frames.map(|frame| {
            dot(
                frame.thrust_axis(),
                stage.carrier_angular_velocity_body_rad_s,
            )
        });
        let rate = |k: usize| {
            self.inflow
                .rate(
                    induction[k],
                    wing.disk_axial_mps[k],
                    wake_force_n[k],
                    frame_rate_rad_s[k],
                )
                .map_err(air_error)
        };
        Ok(HarmonicAerodynamics {
            induction,
            disk_axial_mps: wing.disk_axial_mps,
            frame_rate_rad_s,
            propellers,
            wake_force_n,
            wake_rate: [rate(0)?, rate(1)?],
            gamma: wing.gamma,
            strip_velocity_mps: wing.strip_velocity_mps,
            wing_lift_n: wing.lift_n,
            wing_residual: wing.residual,
        })
    }

    /// All induction, section, gear and wing domains must pass. No state is
    /// clipped, re-equilibrated or silently replaced by uniform induction.
    pub fn evaluate(
        &self,
        state: &HarmonicState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<HarmonicEvaluation<N>, Refusal> {
        if inputs.density_kg_m3 != self.inflow.density_kg_m3()
            || !state
                .induction
                .iter()
                .flat_map(|u| u.modes())
                .all(f64::is_finite)
        {
            return Err(reject(
                "harmonic-propulsion-invalid",
                "finite wake modes and fixed density required",
            ));
        }
        let result =
            self.base
                .propulsion
                .evaluate_with(&state.motion, inputs, wind_world, |stage| {
                    self.solve(stage, state.induction)
                })?;
        Ok(HarmonicEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            aerodynamics: result.aerodynamics,
        })
    }

    fn advance_wakes(
        &self,
        before: [HarmonicInduction; 2],
        stage: &HarmonicAerodynamics,
        dt_s: f64,
    ) -> Result<[HarmonicInflowStep; 2], Refusal> {
        let step = |k: usize| {
            self.inflow
                .step_held(
                    before[k],
                    stage.disk_axial_mps[k],
                    stage.wake_force_n[k],
                    stage.frame_rate_rad_s[k],
                    dt_s,
                )
                .map_err(air_error)
        };
        Ok([step(0)?, step(1)?])
    }

    /// Two real force evaluations predict and advance all original mechanical
    /// and modal states together. Each wake step is implicit in transport;
    /// blade/wing coupling uses the predicted midpoint. The coupling work
    /// defect is explicit and requires timestep refinement. No full fluid/body
    /// passivity, moving far-wake geometry or exact nonlinear trajectory follows.
    /// A late refusal leaves the borrowed state unchanged; callback effects are
    /// caller-owned, as in the original propulsion step.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &HarmonicState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
    ) -> Result<HarmonicStep<N>, Refusal>
    where
        C: FnMut(f64, &HarmonicState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut initial = None;
        let mut midpoint = None;
        let mut stage_induction = state.induction;
        let motion = geared::step(
            self.base.propulsion.assembly,
            &state.motion,
            t_s,
            dt_s,
            angular_origin_world,
            |time, mechanical| {
                let stage = HarmonicState {
                    motion: *mechanical,
                    induction: stage_induction,
                };
                let evaluation = self.evaluate(&stage, inputs(time, &stage)?, &|point| {
                    wind_world(time, point)
                })?;
                let loads = evaluation.loads;
                if initial.is_none() {
                    stage_induction = self
                        .advance_wakes(state.induction, &evaluation.aerodynamics, 0.5 * dt_s)?
                        .map(|wake| wake.induction);
                    initial = Some(evaluation);
                } else {
                    midpoint = Some(evaluation);
                }
                Ok(loads)
            },
        )?;
        let initial = initial
            .ok_or_else(|| reject("harmonic-propulsion-invalid", "missing initial stage"))?;
        let midpoint = midpoint
            .ok_or_else(|| reject("harmonic-propulsion-invalid", "missing midpoint stage"))?;
        let inflow = self.advance_wakes(state.induction, &midpoint.aerodynamics, dt_s)?;
        let coupling_work_residual_j = core::array::from_fn(|k| {
            inflow[k].induced_work_j - dt_s * midpoint.aerodynamics.propellers[k].induced_work_w
        });
        if !coupling_work_residual_j.into_iter().all(f64::is_finite) {
            return Err(reject(
                "harmonic-propulsion-invalid",
                "coupling work overflow",
            ));
        }
        Ok(HarmonicStep {
            state: HarmonicState {
                motion: motion.state,
                induction: inflow.map(|wake| wake.induction),
            },
            motion,
            inflow,
            initial,
            midpoint,
            coupling_work_residual_j,
        })
    }
}
