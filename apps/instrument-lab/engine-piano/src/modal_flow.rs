//! Scalar bank flows lowered from the same prepared modal maps as active strings.
//!
//! The conservative map is the shared owner's exact undamped ZOH map. Bridge
//! forces use average displacements and are closed by the kernel's common Schur
//! solve. Intrinsic velocity loss is the reference's reciprocal midpoint half-flow:
//! absolute mass/gas drag and material strain drag retain their distinct terms.
//! This adapter allocates nothing and introduces no additional retained-state model.

#[derive(Clone, Copy, Debug)]
pub(crate) struct ModalFlow {
    pub(crate) c1: f64,
    pub(crate) c2: f64,
    pub(crate) c3: f64,
    pub(crate) qq: f64,
    pub(crate) qv: f64,
    pub(crate) vq: f64,
    pub(crate) vv: f64,
    pub(crate) ha: f64,
    pub(crate) bqa: f64,
    pub(crate) bva: f64,
    pub(crate) cross: f64,
}

impl ModalFlow {
    /// [new velocity, twice modal kinetic loss, bridge drag source].
    /// Modal loss can be negative during exchange with the bridge; the common
    /// loss solve adds the reciprocal board change before booking physical heat.
    #[inline]
    pub(crate) fn loss(&self, velocity: f64, bridge_mid_velocity: f64) -> [f64; 3] {
        let next = self.c1 * velocity + self.c2 * bridge_mid_velocity;
        [next, (velocity - next) * (velocity + next), self.c3 * velocity]
    }

    /// Free endpoint after the intrinsic pre-half-flow, followed by exact ZOH.
    /// The reaction's bridge-dependent part is included in the common Schur.
    #[inline]
    pub(crate) fn predict(&self, q: f64, v: f64, bridge_mid_velocity: f64) -> ([f64; 2], [f64; 3]) {
        let [velocity, loss, drag] = self.loss(v, bridge_mid_velocity);
        let fq = self.qq * q + self.qv * velocity;
        let fv = self.vq * q + self.vv * velocity;
        ([fq, fv], [loss, drag, self.ha * (q + fq)])
    }

    /// Finish against the SAME accepted average bridge displacement as active
    /// strings. Return the reciprocal average force and post-loss source.
    #[inline]
    pub(crate) fn finish(&self, q: f64, fq: f64, fv: f64, bridge_mid_position: f64)
        -> ([f64; 2], [f64; 2]) {
        let next_q = fq + self.bqa * bridge_mid_position;
        let next_v = fv + self.bva * bridge_mid_position;
        ([next_q, next_v], [self.ha * (q + next_q), self.cross * next_v])
    }
}

#[cfg(test)]
mod tests {
    use super::ModalFlow;
    use std::f64::consts::TAU;

    fn conservative(omega: f64, beta: f64, dt: f64) -> (ModalFlow, f64, f64) {
        // Analytic oscillator coefficients, independent of the production
        // table builder, which lowers fs-modal's admitted public ZOH map.
        let (sin, cos) = (omega * dt).sin_cos();
        let bq = 2.0 * (0.5 * omega * dt).sin().powi(2) / (omega * omega);
        let bv = sin / omega;
        let a = omega * omega * beta;
        (ModalFlow {
            c1: 1.0, c2: 0.0, c3: 0.0, qq: cos, qv: bv, vq: -omega * sin, vv: cos,
            ha: 0.5 * a, bqa: bq * a, bva: bv * a, cross: 0.0,
        }, bq, bv)
    }

    fn coupled_step(flow: &ModalFlow, bq: f64, bv: f64, state: [f64; 4]) -> [f64; 4] {
        let [q, v, b, w] = state;
        let (free, [_, _, reaction]) = flow.predict(q, v, 0.0);
        let board_q = flow.qq * b + flow.qv * w;
        let board_v = flow.vq * b + flow.vv * w;
        let schur_weight = flow.ha * flow.ha * bq;
        let next_b = (board_q / bq + reaction + b * schur_weight)
            / (1.0 / bq - schur_weight);
        let (next, [force, _]) = flow.finish(q, free[0], free[1], 0.5 * (b + next_b));
        [next[0], next[1], next_b, board_v + bv * force]
    }

    #[test]
    fn coupled_normal_modes_and_resonant_beating_converge_without_energy_injection() {
        let omega = TAU * 400.0;
        let beta = 0.08_f64;
        // H = (v^2+w^2+omega^2*(q^2+b^2)-2*a*q*b)/2 is positive.
        // The physical eigenfrequencies are omega*sqrt(1-beta) and
        // omega*sqrt(1+beta); energy transfers between the two coordinates.
        let low = omega * (1.0 - beta).sqrt();
        let high = omega * (1.0 + beta).sqrt();
        let energy0 = 0.5 * omega * omega;
        let mut errors = Vec::new();
        for rate in [32_000, 64_000, 128_000] {
            let dt = 1.0 / f64::from(rate);
            let (flow, bq, bv) = conservative(omega, beta, dt);
            let mut state = [1.0, 0.0, 0.0, 0.0];
            let mut error = 0.0_f64;
            let mut transfer = 0.0_f64;
            for sample in 1..=rate / 40 {
                state = coupled_step(&flow, bq, bv, state);
                let t = f64::from(sample) * dt;
                let (sl, cl) = (low * t).sin_cos();
                let (sh, ch) = (high * t).sin_cos();
                let exact = [0.5 * (cl + ch), -0.5 * (low * sl + high * sh),
                    0.5 * (cl - ch), -0.5 * (low * sl - high * sh)];
                for j in 0..4 {
                    let scale = if j % 2 == 0 { 1.0 } else { omega };
                    error = error.max((state[j] - exact[j]).abs() / scale);
                }
                let [q, v, b, w] = state;
                let energy = 0.5 * (v * v + w * w + omega * omega * (q * q + b * b)
                    - 2.0 * omega * omega * beta * q * b);
                assert!((energy / energy0 - 1.0).abs() < 5e-10);
                transfer = transfer.max(b.abs());
            }
            assert!(transfer > 0.9, "the unexcited resonator must receive the beating energy");
            errors.push(error);
        }
        for pair in errors.windows(2) {
            assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "errors {errors:?}");
        }
        assert!(errors[2] < 2e-4, "errors {errors:?}");
    }

    #[test]
    fn exchanging_excitation_and_receiver_preserves_reciprocity() {
        let omega = TAU * 400.0;
        let (flow, bq, bv) = conservative(omega, 0.08, 1.0 / 64_000.0);
        let mut left = [1.0, 0.0, 0.0, 0.0];
        let mut right = [0.0, 0.0, 1.0, 0.0];
        for _ in 0..1_600 {
            left = coupled_step(&flow, bq, bv, left);
            right = coupled_step(&flow, bq, bv, right);
            assert!((left[2] - right[0]).abs() < 1e-10);
            assert!((left[3] - right[1]).abs() < 1e-10 * omega);
        }
    }

    #[test]
    fn source_loss_matches_the_full_physical_resistance_and_its_positive_work() {
        // R = [[mass+strain, -strain*beta],
        //      [-strain*beta, endpoint+strain*beta^2]].
        // Compare the bank's actual scalar flow and bridge condensation with
        // the independent full 2x2 midpoint solve and h*v_mid^T*R*v_mid.
        let h = 0.001;
        let beta = 0.2;
        for (mass, strain, endpoint, initial) in [
            (20.0, 40.0, 0.5, [0.7, -0.3]),
            (0.0, 6000.0, 0.5, [0.7, -0.3]), // strong material loss, beyond a small h*rate
            (40.0, 0.0, 0.0, [0.2, 1.0]), // absolute gas drag is not a relative drag
            (0.0, 40.0, 0.0, [0.2, 1.0]), // zero material strain rate
            (0.0, 0.0, 5.0, [0.0, 0.5]),  // residual chord alone
        ] {
            let half = 0.5 * h;
            let inv = 1.0 / (1.0 + half * (mass + strain));
            let cross = half * strain * beta * inv;
            let flow = ModalFlow {
                c1: 1.0 - h * (mass + strain) * inv,
                c2: -h * ((mass + strain) * cross - strain * beta),
                c3: strain * beta * inv, cross,
                qq: 1.0, qv: 0.0, vq: 0.0, vv: 1.0, ha: 0.0, bqa: 0.0, bva: 0.0,
            };
            let c4 = strain * beta * (cross - beta);
            let loss_weight = half * strain * beta * beta * (1.0 + half * mass) * inv;
            let bm = (initial[1] + cross * initial[0]) / (1.0 + half * endpoint + loss_weight);
            let [v, twice_modal_loss, drag] = flow.loss(initial[0], bm);
            let w = initial[1] - h * (endpoint * bm - drag - bm * c4);
            let heat = 0.5 * (twice_modal_loss + (initial[1] - w) * (initial[1] + w));

            let r00 = mass + strain;
            let r01 = -strain * beta;
            let r11 = endpoint + strain * beta * beta;
            let (a, b, d) = (1.0 + half * r00, half * r01, 1.0 + half * r11);
            let rhs0 = (1.0 - half * r00) * initial[0] - half * r01 * initial[1];
            let rhs1 = -half * r01 * initial[0] + (1.0 - half * r11) * initial[1];
            let determinant = a * d - b * b;
            let expected = [(d * rhs0 - b * rhs1) / determinant,
                (a * rhs1 - b * rhs0) / determinant];
            assert!((v - expected[0]).abs() < 2e-14);
            assert!((w - expected[1]).abs() < 2e-14);
            let vm = 0.5 * (initial[0] + expected[0]);
            let wm = 0.5 * (initial[1] + expected[1]);
            let physical_heat = h * (mass * vm * vm + strain * (vm - beta * wm).powi(2)
                + endpoint * wm * wm);
            assert!(physical_heat >= 0.0);
            assert!((heat - physical_heat).abs() < 2e-14);
            if mass > 0.0 { assert!(heat > 0.0); }
        }
    }
}
