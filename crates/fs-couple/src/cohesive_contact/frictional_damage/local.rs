//! Shared local damaged-area law: source binding and actual-motion work only.
//! No local mechanical compliance or independent-site solver lives here.
use super::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct LocalDamageFriction<'a> {
    pub(crate) cohesive: &'a ResolvedBilinearCohesive,
    pub(crate) friction: &'a ResidualCoulombSource,
    pub(crate) area: f64,
    pub(crate) config: DamageFrictionConfig,
    pub(crate) identity: ContentHash,
    pub(crate) static_mu: f64,
    pub(crate) kinetic_mu: f64,
}
impl<'a> LocalDamageFriction<'a> {
    pub(crate) fn new(
        cohesive: &'a ResolvedBilinearCohesive,
        friction: &'a ResidualCoulombSource,
        area: f64,
        config: DamageFrictionConfig,
        identity: ContentHash,
    ) -> Result<Self, DamageFrictionError> {
        if !super::super::valid_config(config.contact)
            || !(1..=256).contains(&config.max_friction_iterations)
        {
            return Err(DamageFrictionError::Invalid(
                "bounded original material and friction budgets required",
            ));
        }
        for x in [
            area,
            config.force_absolute_tolerance_n,
            config.work_absolute_tolerance_j,
        ] {
            positive(x)?;
        }
        for x in [
            config.force_relative_tolerance,
            config.work_relative_tolerance,
        ] {
            if !x.is_finite() || x <= 0.0 || x >= 1.0 {
                return Err(DamageFrictionError::Invalid(
                    "relative work/force gates must be in (0,1)",
                ));
            }
        }
        friction.require_cohesive_source(cohesive)?;
        let static_mu = friction.law().static_coefficient()?;
        let kinetic_mu = friction.law().zero_speed_kinetic_coefficient()?;
        Ok(Self {
            cohesive,
            friction,
            area,
            config,
            identity,
            static_mu,
            kinetic_mu,
        })
    }
    pub(crate) fn adopt(
        &self,
        point: CohesivePoint<'a>,
        receiving: [&MaterialStateId; 2],
        activation: ReclosureActivation,
    ) -> Result<DamageFrictionState<'a>, DamageFrictionError> {
        self.cohesive
            .require_surface_materials(receiving[0], receiving[1])?;
        self.friction.require_surface_materials(receiving)?;
        if !core::ptr::eq(point.law(), self.cohesive.law()) {
            return Err(DamageFrictionError::StateMismatch);
        }
        let occupancy = if point.jump_m()[0] >= 0.0 || point.damage() == 0.0 {
            DamageFrictionRegime::Inactive
        } else if activation == ReclosureActivation::Sliding {
            DamageFrictionRegime::Sliding
        } else {
            DamageFrictionRegime::Sticking
        };
        Ok(DamageFrictionState {
            identity: self.identity,
            point,
            occupancy,
        })
    }

    pub(crate) fn audit(
        &self,
        previous: DamageFrictionState<'a>,
        jump: [f64; 3],
        dt_s: f64,
        applied: [f64; 3],
        regime: DamageFrictionRegime,
    ) -> Result<DamageFrictionAudit<'a>, DamageFrictionError> {
        self.require_state(previous)?;
        let old = previous.point.jump_m();
        let delta = core::array::from_fn(|i| jump[i] - old[i]);
        for i in 0..3 {
            finite(jump[i])?;
            finite(applied[i])?;
            if delta[i].abs() > self.config.contact.maximum_jump_increment_m {
                return Err(DamageFrictionError::Range("actual jump increment"));
            }
            if applied[i].abs() > self.config.contact.maximum_force_n {
                return Err(DamageFrictionError::Range("actual total force"));
            }
        }
        let (material, mean, pressure) = previous
            .point
            .trial_over_time_with_damage_weighted_compression(
                jump,
                dt_s,
                self.config.contact.material_energy_tolerance_j_m2,
            )?;
        let pressure = finite(self.area * pressure)?;
        let mut expected = mean.map(|x| self.area * x);
        let slip = delta[1].hypot(delta[2]);
        let mut friction = [0.0; 2];
        let loss = match regime {
            DamageFrictionRegime::Inactive => {
                if pressure != 0.0 {
                    return Err(DamageFrictionError::Range(
                        "inactive path carries damaged compression",
                    ));
                }
                0.0
            }
            DamageFrictionRegime::Sticking => {
                if slip > self.jump_tolerance(old, jump, jump)? {
                    return Err(DamageFrictionError::Range("actual sticking slip"));
                }
                friction = [applied[1] - expected[1], applied[2] - expected[2]];
                let magnitude = norm2(friction);
                let cap = finite(self.capture_mu(previous) * pressure)?;
                if magnitude > cap + self.force_tolerance(magnitude, cap)? {
                    return Err(DamageFrictionError::Range(
                        "actual sticking reaction exceeds capture cone",
                    ));
                }
                0.0
            }
            DamageFrictionRegime::Sliding => {
                if slip == 0.0 {
                    return Err(DamageFrictionError::Range("sliding requires actual slip"));
                }
                let cap = finite(self.kinetic_mu * pressure)?;
                friction = [cap * (delta[1] / slip), cap * (delta[2] / slip)];
                finite(cap * slip)?
            }
        };
        for i in 0..2 {
            expected[i + 1] = finite(expected[i + 1] + friction[i])?;
        }
        let mut force_residual = [0.0; 3];
        let mut force_tolerance = [0.0; 3];
        for i in 0..3 {
            force_residual[i] = finite(applied[i] - expected[i])?;
            force_tolerance[i] = self.force_tolerance(applied[i], expected[i])?;
            if force_residual[i].abs() > force_tolerance[i] {
                return Err(DamageFrictionError::ForceMismatch {
                    component: i,
                    residual_n: force_residual[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let work = finite(dot3(applied, delta))?;
        let material_work = finite(self.area * material.total_work_j_m2())?;
        let defect = finite(work - material_work - loss)?;
        let tolerance = finite(
            self.config.work_absolute_tolerance_j
                + self.config.work_relative_tolerance
                    * work.abs().max(material_work.abs()).max(loss),
        )?;
        if defect.abs() > tolerance {
            return Err(DamageFrictionError::WorkMismatch {
                residual_j: defect,
                tolerance_j: tolerance,
            });
        }
        let occupancy = if jump[0] >= 0.0 || material.equilibrium.point.damage() == 0.0 {
            DamageFrictionRegime::Inactive
        } else {
            regime
        };
        Ok(DamageFrictionAudit {
            state: DamageFrictionState {
                identity: self.identity,
                point: material.equilibrium.point,
                occupancy,
            },
            material,
            regime,
            damaged_compression_force_n: pressure,
            friction_force_n: friction,
            expected_force_n: expected,
            force_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
            friction_dissipation_j: loss,
            work_residual_j: defect,
        })
    }

    pub(crate) fn require_state(
        &self,
        s: DamageFrictionState<'a>,
    ) -> Result<(), DamageFrictionError> {
        if s.identity != self.identity || !core::ptr::eq(s.point.law(), self.cohesive.law()) {
            return Err(DamageFrictionError::StateMismatch);
        }
        Ok(())
    }
    pub(crate) fn capture_mu(&self, s: DamageFrictionState<'a>) -> f64 {
        if s.occupancy == DamageFrictionRegime::Sliding {
            self.kinetic_mu
        } else {
            self.static_mu
        }
    }
    pub(crate) fn force_tolerance(&self, a: f64, b: f64) -> Result<f64, DamageFrictionError> {
        finite(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * a.abs().max(b.abs()),
        )
    }
    pub(crate) fn jump_tolerance(
        &self,
        old: [f64; 3],
        free: [f64; 3],
        jump: [f64; 3],
    ) -> Result<f64, DamageFrictionError> {
        let scale = (0..3)
            .map(|i| (free[i] - old[i]).abs().max((jump[i] - old[i]).abs()))
            .fold(0.0, f64::max);
        finite(
            self.config.contact.jump_absolute_tolerance_m
                + self.config.contact.jump_relative_tolerance * scale,
        )
    }
}
