//! G1/G3: held current-area pressure, hot plastic shape and original inertia.
//! The diagonal return and pressure-volume work are independent scalar oracles.

use super::*;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;

#[derive(Clone, Copy)]
struct PressurePoint {
    axial: f64,
    transverse: f64,
    temperature: f64,
    plastic: f64,
    alpha: f64,
    mu: f64,
    bulk: f64,
    hardening: f64,
    log_theta: f64,
    volume_log: f64,
    deviator: f64,
    free: f64,
    entropy: f64,
    piola: [f64; 3],
    dissipation: f64,
}

impl PressurePoint {
    fn at(axial: f64, transverse: f64, temperature: f64, old: Option<Self>) -> Self {
        let d = temperature - 300.0;
        let (mu, bulk, yield_stress, hardening) = (
            8000.0 - 6.0 * d,
            16000.0 - 8.0 * d,
            200.0 - 0.4 * d,
            800.0 - d,
        );
        let (mut plastic, mut alpha) = old.map_or((0.0, 0.0), |p| (p.plastic, p.alpha));
        let trial = 2.0 * (axial / transverse).ln() / 3.0 - plastic;
        let increment = ((3.0 * mu * trial.abs() - yield_stress - hardening * alpha)
            / (3.0 * mu + hardening))
            .max(0.0);
        plastic += trial.signum() * increment;
        alpha += increment;
        let log_theta = 3.0e-4 * d + 0.5e-7 * d * d;
        let volume_log = axial.ln() + 2.0 * transverse.ln() - 3.0 * log_theta;
        let deviator = 2.0 * (axial / transverse).ln() / 3.0 - plastic;
        let free = 1.5 * mu * deviator * deviator
            + 0.5 * bulk * volume_log * volume_log
            + 0.5 * hardening * alpha * alpha;
        let entropy = 3.0 * bulk * (3.0e-4 + 1.0e-7 * d) * volume_log
            + 9.0 * deviator * deviator
            + 4.0 * volume_log * volume_log
            + 0.5 * alpha * alpha;
        Self {
            axial,
            transverse,
            temperature,
            plastic,
            alpha,
            mu,
            bulk,
            hardening,
            log_theta,
            volume_log,
            deviator,
            free,
            entropy,
            piola: [
                (bulk * volume_log + 2.0 * mu * deviator) / axial,
                (bulk * volume_log - mu * deviator) / transverse,
                (bulk * volume_log - mu * deviator) / transverse,
            ],
            dissipation: yield_stress * increment,
        }
    }

    fn pressure(self) -> f64 {
        -self.piola[0] / (self.transverse * self.transverse)
    }

    fn solve(old: Self, temperature: f64, transverse: f64, pressure: f64) -> Self {
        let (mut lower, mut upper) = (0.5, 1.4);
        assert!(Self::at(lower, transverse, temperature, Some(old)).pressure() > pressure);
        assert!(Self::at(upper, transverse, temperature, Some(old)).pressure() < pressure);
        for _ in 0..80 {
            let axial = 0.5 * (lower + upper);
            if Self::at(axial, transverse, temperature, Some(old)).pressure() > pressure {
                lower = axial;
            } else {
                upper = axial;
            }
        }
        Self::at(0.5 * (lower + upper), transverse, temperature, Some(old))
    }

    fn caloric_heat(self, old: Self) -> f64 {
        let coefficient = (self.mu - old.mu) * 1.5 * old.deviator * old.deviator
            + 0.5 * (self.bulk - old.bulk) * old.volume_log * old.volume_log
            + 0.5 * (self.hardening - old.hardening) * old.alpha * old.alpha;
        let thermal_work = (self.piola[0] * old.axial
            + 2.0 * self.piola[1] * old.transverse)
            * (self.log_theta - old.log_theta).exp_m1();
        self.dissipation + thermal_work - coefficient
            - (self.temperature * self.entropy - old.temperature * old.entropy)
    }

    fn jacobian(self) -> f64 {
        self.axial * self.transverse * self.transverse
    }
}

fn transverse_supports(stretch: f64) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| *dof != 3)
        .map(|dof| (dof, if dof == 7 || dof == 11 { stretch - 1.0 } else { 0.0 }))
        .collect()
}

fn surface(cx: &Cx<'_>) -> ClosedTetPressureSurface {
    ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap()
}

#[test]
fn pressure_driven_hot_member_retains_finite_plastic_shape_on_release() {
    with_cx(|cx| {
        let surface = surface(cx);
        let mut solid = initial(cx);
        let mut point = PressurePoint::at(1.0, 1.0, 300.0, None);
        let load = PressurePoint::at(0.85, 1.0, 300.0, None).pressure();
        let hot_transverse = (3.0e-4_f64 * 300.0 + 0.5e-7 * 300.0 * 300.0).exp();
        for (temperature, transverse, pressure) in [
            (300.0, 1.0, load),
            (600.0, hot_transverse, load),
            (600.0, hot_transverse, 0.0),
        ] {
            let next = PressurePoint::solve(point, temperature, transverse, pressure);
            let step = solid
                .advance_with_pressure(
                    &[temperature],
                    &transverse_supports(transverse),
                    &[],
                    [0.0; 3],
                    0.1,
                    &surface,
                    pressure,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            close(step.solid.displacement_m()[1][0], next.axial - 1.0, 2e-9);
            close(step.solid.element_states()[0].alpha, next.alpha, 2e-9);
            assert!(next.alpha > point.alpha);
            for i in 0..3 {
                close(
                    step.solid.element_states()[0].inverse_plastic_cauchy_green[i],
                    (if i == 0 { -2.0 * next.plastic } else { next.plastic }).exp(),
                    2e-9,
                );
                close(
                    step.solid.element_first_piola_pa()[0][4 * i],
                    next.piola[i],
                    1e-5,
                );
            }
            let account = step.pressure.unwrap();
            let work = -pressure * V * (next.jacobian() - point.jacobian());
            let endpoint = -pressure * V
                * (transverse * transverse * (next.axial - point.axial)
                    + 2.0 * next.axial * transverse * (transverse - point.transverse));
            close(account.work_j, work, 2e-6);
            close(account.endpoint_work_j, endpoint, 2e-6);
            close(account.endpoint_work_correction_j, work - endpoint, 2e-6);
            close(
                account.external_work_defect_j,
                step.energy.endpoint_work_defect_j + account.endpoint_work_correction_j,
                1e-10,
            );
            close(step.solid.free_energy_j(), V * next.free, 2e-6);
            close(step.energy.plastic_dissipation_j, V * next.dissipation, 2e-6);
            close(step.energy.caloric_heat_j, V * next.caloric_heat(point), 2e-5);
            close(step.energy.equilibrium_work_residual_j, 0.0, 1e-6);
            close(step.energy.balance_error_j, 0.0, 1e-6);
            if temperature > 300.0 && pressure > 0.0 {
                assert!(account.endpoint_work_correction_j > 1.0);
            }
            if pressure == 0.0 {
                assert!(next.axial < hot_transverse - 0.005);
            }
            assert_eq!(step.solid.temperatures_k(), &[temperature]);
            solid = step.solid;
            point = next;
        }
        assert_eq!(solid.families(), &[family()]);
    });
}

#[test]
fn follower_pressure_keeps_reference_inertia_and_refusals_out_of_heat() {
    with_cx(|cx| {
        let surface = surface(cx);
        let density = 1200.0;
        let duration: f64 = 0.08;
        let old = PressurePoint::at(1.0, 1.0, 300.0, None);
        let next = PressurePoint::at(0.92, 1.04, 350.0, Some(old));
        let solid = initial(cx).with_inertia(&[density], &[[0.0; 3]; 4], cx).unwrap();
        let acceleration = (next.axial - 1.0) / duration.powi(2);
        // The free node's consistent P1 mass is rho*V/10; pressure acts on
        // the current transverse area, not the reference area.
        let pressure = -(next.piola[0] + density * acceleration / 10.0)
            / next.transverse.powi(2);
        let supports = transverse_supports(next.transverse);
        let old_history = solid.element_states().to_vec();
        let old_displacement = solid.displacement_m().to_vec();
        let old_velocity = solid.velocity_mps().unwrap().to_vec();
        let mut strict = controls();
        strict.maximum_numerical_kinetic_loss_j = 0.0;
        assert!(solid.advance_with_pressure(
            &[350.0], &supports, &[], [0.0; 3], duration, &surface, pressure,
            settings(), strict, cx,
        ).is_err());
        let mut foreign_nodes = NODES;
        foreign_nodes[1][0] = 1.01;
        let foreign = ClosedTetPressureSurface::new(
            &foreign_nodes, &TETS, TetAssemblyBudget::standard(), 4, cx,
        ).unwrap();
        assert!(solid.advance_with_pressure(
            &[350.0], &supports, &[], [0.0; 3], duration, &foreign, pressure,
            settings(), controls(), cx,
        ).is_err());
        assert_eq!(solid.element_states(), old_history);
        assert_eq!(solid.displacement_m(), old_displacement);
        assert_eq!(solid.velocity_mps().unwrap(), old_velocity);
        assert_eq!(solid.temperatures_k(), &[300.0]);

        let step = solid.advance_with_pressure(
            &[350.0], &supports, &[], [0.0; 3], duration, &surface, pressure,
            settings(), controls(), cx,
        ).unwrap();
        close(step.solid.displacement_m()[1][0], next.axial - 1.0, 2e-9);
        close(step.solid.element_states()[0].alpha, next.alpha, 2e-9);
        close(step.solid.reference_mass_kg().unwrap(), density * V, 1e-12);
        let axial_velocity = (next.axial - 1.0) / duration;
        let transverse_velocity = (next.transverse - 1.0) / duration;
        let kinetic = density * V / 20.0
            * (axial_velocity.powi(2) + 2.0 * transverse_velocity.powi(2));
        let inertia = step.inertia.unwrap();
        close(inertia.kinetic_change_j, kinetic, 1e-6);
        close(inertia.numerical_loss_j, kinetic, 1e-6);
        assert!(inertia.numerical_loss_j > 10.0);
        close(step.energy.caloric_heat_j, V * next.caloric_heat(old), 2e-5);
        close(step.energy.balance_error_j, 0.0, 1e-6);
        let pressure_work = step.pressure.unwrap();
        close(
            pressure_work.work_j,
            -pressure * V * (next.jacobian() - 1.0),
            2e-6,
        );
    });
}
