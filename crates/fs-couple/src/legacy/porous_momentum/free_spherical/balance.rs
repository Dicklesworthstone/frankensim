//! Independent common-origin momentum and native work checks before publication.

use super::*;

impl Step<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        self,
        fluid: BoundaryGrid3,
        trial: Trial,
        nonlinear: NewtonReport,
        native_trials: usize,
        trial_cells: usize,
        gate: &CancelGate,
    ) -> Result<FreePorousStep, PorousMomentumError> {
        let m = self.motion;
        let c = self.controls.flow;
        let body = trial.mechanics.state_after;
        let position0 = m.body.pose().position_world();
        let position1 = body.pose().position_world();
        // Together with the independent endpoint orbital inventory, this checks
        // the identity Δ(r × p) = r_mid × Δp used by the shared force moment.
        let drift_residual = v(std::array::from_fn(|i| {
            (dd(a(position1)[i])
                - dd(a(position0)[i])
                - dd(m.units.timestep_s) * dd(a(trial.velocity)[i]))
            .to_f64()
        }));
        check(
            "free porous native midpoint translation",
            norm(drift_residual),
            c.position_tolerance_m,
        )?;
        check(
            "free porous native Cayley work velocity",
            norm(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(trial.omega_body)
                    .sub(trial.omega_world),
            ) * m.radius_m,
            c.midpoint_tolerance_m_s,
        )?;
        let exchange = trial.exchange;
        let impulse = v(exchange.skeleton_impulse).scale(self.units.impulse);
        let moment = v(exchange.skeleton_angular_impulse).scale(self.units.angular);
        let arm = difference(trial.centre, c.balance_origin_m);
        let shift = |moment: Vec3, impulse: Vec3| {
            let orbital = cross(a(arm).map(dd), a(impulse).map(dd));
            v(std::array::from_fn(|i| {
                (dd(a(moment)[i]) + orbital[i]).to_f64()
            }))
        };
        let drag_angular = shift(moment, impulse);
        let external_angular = shift(self.external_couple, self.external_linear);
        let orbital = |body: RigidBodyState| {
            let r = std::array::from_fn(|i| {
                dd(a(body.pose().position_world())[i]) - dd(a(c.balance_origin_m)[i])
            });
            cross(r, a(body.linear_momentum_world()).map(dd))
        };
        let orbital0 = orbital(m.body);
        let orbital1 = orbital(body);
        let body_angular_change = v(std::array::from_fn(|i| {
            (orbital1[i] - orbital0[i] + dd(a(trial.spin_change)[i])).to_f64()
        }));
        let balance = |change: Vec3, first: Vec3, second: Vec3| {
            v(std::array::from_fn(|i| {
                (dd(a(change)[i]) - dd(a(first)[i]) - dd(a(second)[i])).to_f64()
            }))
        };
        let body_linear = balance(trial.linear_change, impulse, self.external_linear);
        let body_angular = balance(body_angular_change, drag_angular, external_angular);

        let change = fluid_change(&m.fluid, &fluid, m.units, c, gate)?;
        let fluid_external = v(exchange.external_impulse).scale(self.units.impulse);
        let fluid_external_angular = shift(
            v(exchange.external_angular_impulse).scale(self.units.angular),
            fluid_external,
        );
        let wall = v(exchange.walls.fluid_population_impulse).scale(self.units.impulse);
        let wall_angular = shift(
            v(exchange.walls.fluid_population_angular_impulse).scale(self.units.angular),
            wall,
        );
        let fluid_linear = v(std::array::from_fn(|i| {
            (dd(a(change.linear)[i]) - dd(a(fluid_external)[i]) + dd(a(impulse)[i])
                - dd(a(wall)[i]))
            .to_f64()
        }));
        let fluid_angular = v(std::array::from_fn(|i| {
            (dd(a(change.angular)[i]) - dd(a(fluid_external_angular)[i]) + dd(a(drag_angular)[i])
                - dd(a(wall_angular)[i]))
            .to_f64()
        }));
        let absolute_linear = norm(body_linear) + norm(fluid_linear);
        let absolute_angular = norm(body_angular) + norm(fluid_angular);
        check(
            "free porous fluid/body linear momentum",
            absolute_linear,
            c.linear_tolerance_kg_m_s,
        )?;
        check(
            "free porous fluid/body orbital and spin momentum",
            absolute_angular,
            c.angular_tolerance_kg_m2_s,
        )?;

        // Difference of the actual endpoint quadratic stores without subtracting
        // their large absolute energies. Full native inertia enters omega_body.
        let kinetic = (0..3)
            .fold(Dd::ZERO, |sum, i| {
                sum + dd(a(trial.velocity)[i]) * dd(a(trial.linear_change)[i])
                    + dd(a(trial.omega_body)[i])
                        * (dd(a(body.angular_momentum_body())[i])
                            - dd(a(m.body.angular_momentum_body())[i]))
            })
            .to_f64();
        let fluid_work = finite(exchange.external_work * self.units.energy)?;
        let drag_work = finite(exchange.fluid_drag_work * self.units.energy)?;
        let skeleton_work = finite(exchange.skeleton_work * self.units.energy)?;
        let loss = finite(exchange.relative_dissipation * self.units.energy)?;
        if loss < 0.0 {
            return Err(invalid("native porous dissipation is negative"));
        }
        let work = |impulse: Vec3, moment: Vec3| {
            (dd(dot(trial.velocity, impulse)) + dd(dot(trial.omega_world, moment))).to_f64()
        };
        let external_work = work(self.external_linear, self.external_couple);
        let body_drag_work = work(impulse, moment);
        let constitutive = (dd(drag_work) + dd(skeleton_work) + dd(loss)).to_f64();
        let interface = (dd(body_drag_work) - dd(skeleton_work)).to_f64();
        let rigid = (dd(kinetic) - dd(body_drag_work) - dd(external_work)).to_f64();
        check(
            "noncancelling free porous constitutive/interface/rigid work",
            constitutive.abs() + interface.abs() + rigid.abs(),
            c.work_tolerance_j,
        )?;
        let mechanical = (dd(kinetic) + dd(drag_work) + dd(loss) - dd(external_work)).to_f64();
        let remainder = (dd(change.kinetic) - dd(fluid_work) - dd(drag_work)).to_f64();
        let linear_residual = body_linear.add(fluid_linear);
        let angular_residual = body_angular.add(fluid_angular);
        if [kinetic, external_work, mechanical, remainder]
            .iter()
            .any(|x| !x.is_finite())
            || !linear_residual.is_finite()
            || !angular_residual.is_finite()
        {
            return Err(invalid(
                "free porous momentum/work account is not representable",
            ));
        }
        poll(gate)?;
        Ok(FreePorousStep {
            motion: FreePorousBody {
                fluid,
                body,
                properties: m.properties,
                units: m.units,
                radius_m: m.radius_m,
                support: m.support,
                links: m.links.clone(),
                cells: m.cells,
                time_s: self.time_s,
                accepted_steps: self.accepted_steps,
            },
            exchange,
            exchange_origin_m: trial.centre,
            mechanics: trial.mechanics,
            nonlinear,
            native_trials,
            trial_cells,
            midpoint_velocity_m_s: trial.velocity,
            angular_velocity_world_rad_s: trial.omega_world,
            balance_origin_m: c.balance_origin_m,
            fluid_linear_change_kg_m_s: change.linear,
            fluid_angular_change_kg_m2_s: change.angular,
            body_linear_change_kg_m_s: trial.linear_change,
            body_spin_change_kg_m2_s: trial.spin_change,
            body_angular_change_kg_m2_s: body_angular_change,
            body_external_impulse_kg_m_s: self.external_linear,
            body_external_angular_impulse_kg_m2_s: external_angular,
            linear_residual_kg_m_s: linear_residual,
            angular_residual_kg_m2_s: angular_residual,
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            fluid_mass_change_kg: change.mass,
            fluid_kinetic_change_j: change.kinetic,
            body_kinetic_change_j: kinetic,
            body_external_work_j: external_work,
            fluid_external_work_j: fluid_work,
            fluid_drag_work_j: drag_work,
            relative_dissipation_j: loss,
            other_fluid_kinetic_change_j: remainder,
            mechanical_work_residual_j: mechanical,
        })
    }
}
