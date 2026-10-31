//! Differentiate the source-owned finite-step effort itself. A Hessian of H
//! is not generally the end-state Jacobian of that constitutive effort.

use super::*;
use crate::ConstitutiveDiscreteGradient;

impl StepWorkspace {
    pub(super) fn constitutive_jacobian_into<F>(
        &mut self,
        sys: &PortHamiltonian,
        x0: &[f64],
        dt: f64,
        rule: &dyn ConstitutiveDiscreteGradient,
        dissipation: Option<Dissipation<'_>>,
        poll: &mut F,
    ) -> Result<(), PreparedStepError>
    where
        F: FnMut() -> Result<(), PreparedStepError>,
    {
        let n = self.n;
        for i in 0..n {
            self.midpoint[i] = f64::midpoint(x0[i], self.x[i]);
        }
        self.effort.fill(0.0);
        rule.evaluate(x0, &self.x, &self.midpoint, &mut self.effort);
        norm(&self.effort)?;
        if let Some(plan) = &mut self.condensation {
            // The actual constitutive Jacobian is assembled directly as the
            // base. No unrelated global Gonzalez rank-one term is present.
            plan.left.fill(0.0);
            plan.right.fill(0.0);
        }
        self.trial.fill(0.0);
        for col in 0..n {
            poll()?;
            self.trial[col] = 1.0;
            self.plus.fill(0.0);
            let available = rule.end_jacobian_action(x0, &self.x, &self.trial, &mut self.plus);
            self.trial[col] = 0.0;
            if !available {
                return Err(dimensions(
                    "source-owned discrete gradient lacks its analytic end-state Jacobian action",
                )
                .into());
            }
            norm(&self.plus)?;
            if let Some(port) = dissipation {
                // The dissipative law sees the original midpoint state and
                // actual constitutive effort, with independent derivatives.
                self.trial[col] = 0.5;
                let result = port.directional(
                    &self.midpoint,
                    &self.effort,
                    &self.trial,
                    &self.plus,
                    &mut self.nonlinear_loss,
                );
                self.trial[col] = 0.0;
                result?;
            }
            if let Some(plan) = &mut self.condensation {
                let products = self
                    .flow
                    .apply_direction(sys, &self.plus, &mut plan.flow, poll)?;
                self.jacobian_flow_products = self.jacobian_flow_products.saturating_add(products);
                for row in 0..n {
                    let flow = plan.flow[row]
                        - if dissipation.is_some() {
                            self.nonlinear_loss[row]
                        } else {
                            0.0
                        };
                    self.jacobian[row * n + col] = (if row == col { 1.0 } else { 0.0 }) - dt * flow;
                }
            } else {
                for row in 0..n {
                    let mut flow = 0.0;
                    for &k in self.flow.row(row) {
                        flow += (sys.j[row * n + k] - sys.r[row * n + k]) * self.plus[k];
                    }
                    self.jacobian_flow_products = self
                        .jacobian_flow_products
                        .saturating_add(self.flow.row(row).len());
                    if dissipation.is_some() {
                        flow -= self.nonlinear_loss[row];
                    }
                    self.jacobian[row * n + col] = (if row == col { 1.0 } else { 0.0 }) - dt * flow;
                }
            }
        }
        norm(&self.jacobian)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Storage;

    struct IndependentStores {
        analytic: bool,
    }
    impl Storage for IndependentStores {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.25 * x[0].powi(4) + 0.5 * x[1].powi(2) + 0.5 * x[2].powi(4) + 0.5 * x[3].powi(2)
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out[0] = x[0].powi(3);
            out[1] = x[1];
            out[2] = 2.0 * x[2].powi(3);
            out[3] = x[3];
        }
        fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
            Some(self)
        }
    }
    impl ConstitutiveDiscreteGradient for IndependentStores {
        fn evaluate(&self, a: &[f64], b: &[f64], mid: &[f64], out: &mut [f64]) {
            assert_eq!((a.len(), b.len(), mid.len(), out.len()), (4, 4, 4, 4));
            for (i, factor) in [(0, 0.25), (2, 0.5)] {
                out[i] = factor * (a[i] + b[i]) * (a[i] * a[i] + b[i] * b[i]);
                out[i + 1] = mid[i + 1];
            }
        }
        fn end_jacobian_action(&self, a: &[f64], b: &[f64], d: &[f64], out: &mut [f64]) -> bool {
            assert_eq!((a.len(), b.len(), d.len(), out.len()), (4, 4, 4, 4));
            if !self.analytic {
                return false;
            }
            for (i, factor) in [(0, 0.25), (2, 0.5)] {
                out[i] = factor * (a[i] * a[i] + 2.0 * a[i] * b[i] + 3.0 * b[i] * b[i]) * d[i];
                out[i + 1] = 0.5 * d[i + 1];
            }
            true
        }
    }
    fn system(analytic: bool) -> PortHamiltonian {
        let mut j = vec![0.0; 16];
        j[1] = 1.0;
        j[4] = -1.0;
        j[11] = 1.0;
        j[14] = -1.0;
        j[7] = 0.4;
        j[13] = -0.4;
        let mut r = vec![0.0; 16];
        r[5] = 0.2;
        r[15] = 0.3;
        PortHamiltonian::new(
            4,
            1,
            j,
            r,
            vec![0.0, 1.0, 0.0, 0.5],
            Box::new(IndependentStores { analytic }),
        )
        .unwrap()
    }
    fn resist(x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        out[1] = (0.2 + x[0] * x[0]) * e[1];
        out[3] = (0.3 + x[2] * x[2]) * e[3];
        true
    }
    fn tangent(x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        for i in [0, 2] {
            let offset = if i == 0 { 0.2 } else { 0.3 };
            out[i + 1] = 2.0 * x[i] * dx[i] * e[i + 1] + (offset + x[i] * x[i]) * de[i + 1];
        }
        true
    }
    fn forbidden_hessian(_: &[f64], _: &[f64], _: &mut [f64]) -> bool {
        panic!("a source-owned end-state Jacobian must never dispatch to the Hessian")
    }

    #[test]
    fn owned_jacobian_matches_residual_with_midpoint_dissipation_and_condensation() {
        let sys = system(true);
        let a = [0.25, 0.4, -0.13, 0.1];
        let b = [0.65, 0.31, -0.37, 0.22];
        let dt = 0.03;
        let mut dense = StepWorkspace::new(&sys).unwrap();
        let mut condensed = StepWorkspace::new(&sys).unwrap();
        condensed.set_condensed_pairs(&[[2, 3]]).unwrap();
        for use_dissipation in [false, true] {
            let dissipation = use_dissipation.then_some(Dissipation {
                force: &resist,
                tangent: Some(&tangent),
            });
            for work in [&mut dense, &mut condensed] {
                work.flow.refresh(&sys).unwrap();
                work.x.copy_from_slice(&b);
                if let Some(plan) = &mut work.condensation {
                    plan.left.fill(7.0);
                    plan.right.fill(9.0);
                }
                work.analytic_jacobian_into(
                    &sys,
                    &a,
                    dt,
                    &forbidden_hessian,
                    dissipation,
                    &mut || Ok(()),
                )
                .unwrap();
            }
            let plan = condensed.condensation.as_ref().unwrap();
            assert!(plan.left.iter().chain(&plan.right).all(|&v| v == 0.0));
            for (actual, expected) in condensed.jacobian.iter().zip(&dense.jacobian) {
                assert!((actual - expected).abs() < 1e-14);
            }
            for col in 0..4 {
                let (mut hi, mut lo) = (b, b);
                let h = 1e-6;
                hi[col] += h;
                lo[col] -= h;
                let gp = crate::discrete_gradient(sys.storage.as_ref(), &a, &hi);
                let gm = crate::discrete_gradient(sys.storage.as_ref(), &a, &lo);
                let mut dp = [0.0; 4];
                let mut dm = [0.0; 4];
                if use_dissipation {
                    let mp = std::array::from_fn::<f64, 4, _>(|i| f64::midpoint(a[i], hi[i]));
                    let mm = std::array::from_fn::<f64, 4, _>(|i| f64::midpoint(a[i], lo[i]));
                    resist(&mp, &gp, &mut dp);
                    resist(&mm, &gm, &mut dm);
                }
                for row in 0..4 {
                    let mut flow = -(dp[row] - dm[row]) / (2.0 * h);
                    for k in 0..4 {
                        flow +=
                            (sys.j[row * 4 + k] - sys.r[row * 4 + k]) * (gp[k] - gm[k]) / (2.0 * h);
                    }
                    let expected = (if row == col { 1.0 } else { 0.0 }) - dt * flow;
                    assert!(
                        (dense.jacobian[row * 4 + col] - expected).abs() < 2e-9,
                        "column {col}, row {row}"
                    );
                }
            }
        }
    }

    #[test]
    fn missing_owned_action_refuses_and_shared_storage_forwards_the_exact_rule() {
        let x = [0.25, 0.4, -0.13, 0.1];
        let unavailable = system(false);
        let mut work = StepWorkspace::new(&unavailable).unwrap();
        let (mut next, mut output) = ([123.0; 4], [456.0]);
        assert!(
            work.step_into_analytic(
                &unavailable,
                &x,
                &[0.2],
                0.01,
                &mut next,
                &mut output,
                &forbidden_hessian
            )
            .is_err()
        );
        assert_eq!(next, [123.0; 4]);
        assert_eq!(output, [456.0]);

        let mut original = system(true);
        let shared = original.share_storage();
        assert!(std::ptr::eq(
            original.storage().discrete_gradient_rule().unwrap(),
            shared.storage().discrete_gradient_rule().unwrap(),
        ));
        let mut direct = StepWorkspace::new(&original).unwrap();
        let mut forwarded = StepWorkspace::new(&shared).unwrap();
        let (mut a, mut b, mut ya, mut yb) = ([0.0; 4], [0.0; 4], [0.0], [0.0]);
        let first = direct
            .step_into_analytic(
                &original,
                &x,
                &[0.2],
                0.01,
                &mut a,
                &mut ya,
                &forbidden_hessian,
            )
            .unwrap();
        let second = forwarded
            .step_into_analytic(
                &shared,
                &x,
                &[0.2],
                0.01,
                &mut b,
                &mut yb,
                &forbidden_hessian,
            )
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(ya, yb);
        assert_eq!(first.balance_residual(), second.balance_residual());
        assert!(first.balance_residual().abs() < 1e-10);
    }

    #[test]
    fn sum_rule_keeps_owned_component_work_and_exact_quadratic_sibling_action() {
        let law = crate::SumStorage {
            a: (Box::new(IndependentStores { analytic: true }), 4),
            b: (
                Box::new(crate::QuadraticStorage::new(vec![3.0, 0.5, 0.5, 2.0], 2).unwrap()),
                2,
            ),
        };
        let a = [0.25, 0.4, -0.13, 0.1, 0.7, -0.4];
        let b = [0.65, 0.31, -0.37, 0.22, 0.2, 0.6];
        let effort = crate::discrete_gradient(&law, &a, &b);
        let owned =
            crate::discrete_gradient(&IndependentStores { analytic: true }, &a[..4], &b[..4]);
        assert_eq!(effort[..4], owned);
        let midpoint = [f64::midpoint(a[4], b[4]), f64::midpoint(a[5], b[5])];
        assert!((effort[4] - 3.0 * midpoint[0] - 0.5 * midpoint[1]).abs() < 1e-14);
        assert!((effort[5] - 0.5 * midpoint[0] - 2.0 * midpoint[1]).abs() < 1e-14);
        let work: f64 = effort
            .iter()
            .zip(a)
            .zip(b)
            .map(|((&e, a), b)| e * (b - a))
            .sum();
        assert!((work - (law.hamiltonian(&b) - law.hamiltonian(&a))).abs() < 1e-14);
        let direction = [0.2, -0.3, 0.4, 0.5, -0.6, 0.7];
        let mut action = [0.0; 6];
        assert!(law.discrete_gradient_rule().unwrap().end_jacobian_action(
            &a,
            &b,
            &direction,
            &mut action
        ));
        let (mut hi, mut lo) = (b, b);
        for i in 0..6 {
            hi[i] += 1e-6 * direction[i];
            lo[i] -= 1e-6 * direction[i];
        }
        let gp = crate::discrete_gradient(&law, &a, &hi);
        let gm = crate::discrete_gradient(&law, &a, &lo);
        for i in 0..6 {
            assert!((action[i] - (gp[i] - gm[i]) / 2e-6).abs() < 2e-9);
        }

        let incomplete = crate::SumStorage {
            a: (Box::new(IndependentStores { analytic: true }), 4),
            b: (
                Box::new(crate::SeparableStorage::new(vec![crate::ScalarLaw {
                    h: |x| 0.25 * x.powi(4),
                    dh: |x| x.powi(3),
                }])),
                1,
            ),
        };
        assert!(
            !incomplete
                .discrete_gradient_rule()
                .unwrap()
                .end_jacobian_action(&a[..5], &b[..5], &direction[..5], &mut [0.0; 5])
        );
    }

    #[test]
    fn differential_rule_keeps_the_exact_prefix_and_zero_algebraic_storage() {
        let mut j = vec![0.0; 25];
        j[1] = 1.0;
        j[5] = -1.0;
        j[13] = 1.0;
        j[17] = -1.0;
        j[9] = 1.0;
        j[21] = -1.0;
        let system = crate::DescriptorPortHamiltonian::new(
            5,
            4,
            0,
            j,
            vec![0.0; 25],
            vec![],
            Box::new(IndependentStores { analytic: true }),
        )
        .unwrap();
        let a = [0.25, 0.4, -0.13, 0.1, 1e20];
        let b = [0.65, 0.31, -0.37, 0.22, -1e20];
        let expected =
            crate::discrete_gradient(&IndependentStores { analytic: true }, &a[..4], &b[..4]);
        let actual = crate::discrete_gradient(system.storage(), &a, &b);
        assert_eq!(actual[..4], expected);
        assert_eq!(actual[4], 0.0);
        let rule = system.storage().discrete_gradient_rule().unwrap();
        let mut tail = [7.0; 5];
        assert!(rule.end_jacobian_action(&a, &b, &[0.0, 0.0, 0.0, 0.0, 1.0], &mut tail));
        assert_eq!(tail, [0.0; 5]);
        let dt = 0.03;
        let jacobian = crate::descriptor_constitutive_jacobian(&system, &a, &b, dt).unwrap();
        for col in 0..5 {
            let mut direction = [0.0; 5];
            direction[col] = 1.0;
            let mut effort_derivative = [0.0; 5];
            assert!(rule.end_jacobian_action(&a, &b, &direction, &mut effort_derivative));
            if col == 4 {
                effort_derivative[4] = 0.5;
            }
            for row in 0..5 {
                let flow: f64 = (0..5)
                    .map(|k| system.j[row * 5 + k] * effort_derivative[k])
                    .sum();
                let expected = if row < 4 {
                    (if row == col { 1.0 } else { 0.0 }) - dt * flow
                } else {
                    flow
                };
                assert!((jacobian[row * 5 + col] - expected).abs() < 1e-14);
            }
        }
    }
}
