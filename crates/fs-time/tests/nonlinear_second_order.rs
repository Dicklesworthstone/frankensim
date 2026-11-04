//! G3: nonlinear dissipative resistance through initialization, Newton and replay adjoints.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    initialization::{
        InitialSolveConfig, second_order_acceleration, second_order_acceleration_vjp,
    },
    second_order_adjoint::{
        SecondOrderAdjointConfig, SecondOrderVjp,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
            StructuralReplayBudget, StructuralTrajectoryModel,
        },
    },
};

const H: f64 = 0.07;
const RHO: f64 = 0.4;
const WORKSPACE: usize = 10_000;
const TERMINAL: ([f64; 1], [f64; 1], [f64; 1]) = ([0.7], [-0.4], [0.2]);

// R = k*q + c*(1+q^2)*v + d*v^3. For positive c,d, the damping
// dissipates c*(1+q^2)*v^2 + d*v^4, while depending on both q and v.
#[derive(Clone, Copy)]
struct Damper([f64; 4]); // mass, stiffness, linear damping, cubic damping

impl Damper {
    fn value(self, q: f64, v: f64) -> f64 {
        let [_, k, c, d] = self.0;
        k * q + c * (1.0 + q * q) * v + d * v * v * v
    }
    fn displacement_derivative(self, q: f64, v: f64) -> f64 {
        self.0[1] + 2.0 * self.0[2] * q * v
    }
    fn velocity_derivative(self, q: f64, v: f64) -> f64 {
        self.0[2] * (1.0 + q * q) + 3.0 * self.0[3] * v * v
    }
}
impl SecondOrderProblem for Damper {
    fn dimension(&self) -> usize {
        1
    }
    fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
        out[0] = self.0[0] * x[0];
    }
    // Legacy actions describe only the linear part. The full residual and
    // every derivative path must use the resistance overrides below.
    fn damping_apply(&self, x: &[f64], out: &mut [f64]) {
        out[0] = self.0[2] * x[0];
    }
    fn internal_force(&self, q: &[f64], out: &mut [f64]) {
        out[0] = self.0[1] * q[0];
    }
    fn tangent_apply(&self, _: &[f64], x: &[f64], out: &mut [f64]) {
        out[0] = self.0[1] * x[0];
    }
    fn resistance(&self, q: &[f64], v: &[f64], out: &mut [f64]) {
        out[0] = self.value(q[0], v[0]);
    }
    fn resistance_displacement_apply(&self, q: &[f64], v: &[f64], x: &[f64], out: &mut [f64]) {
        out[0] = self.displacement_derivative(q[0], v[0]) * x[0];
    }
    fn resistance_velocity_apply(&self, q: &[f64], v: &[f64], x: &[f64], out: &mut [f64]) {
        out[0] = self.velocity_derivative(q[0], v[0]) * x[0];
    }
}
impl SecondOrderVjp for Damper {
    fn parameter_count(&self) -> usize {
        4
    }
    fn mass_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.mass_apply(seed, out);
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.damping_apply(seed, out);
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        self.tangent_apply(q, seed, out);
        Ok(())
    }
    fn resistance_displacement_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        self.resistance_displacement_apply(q, v, seed, out);
        Ok(())
    }
    fn resistance_velocity_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        self.resistance_velocity_apply(q, v, seed, out);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        let derivatives = [a[0], q[0], (1.0 + q[0] * q[0]) * v[0], v[0].powi(3)];
        for (value, derivative) in out.iter_mut().zip(derivatives) {
            *value = derivative * seed[0];
        }
        Ok(())
    }
}
impl StructuralTrajectoryModel for Damper {
    fn forcing(&self, time: f64, out: &mut [f64]) -> Result<(), String> {
        out[0] = load(time);
        Ok(())
    }
    fn forcing_vjp(&self, _: f64, _: &[f64], out: &mut [f64]) -> Result<(), String> {
        out.fill(0.0);
        Ok(())
    }
}
struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, input: &[f64], output: &mut [f64]) {
        output.copy_from_slice(input);
    }
}
fn model() -> Damper {
    Damper([1.6, 2.3, 0.65, 0.4])
}
fn method() -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        1,
        H,
        RHO,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-13,
                relative_tolerance: 1e-12,
                linear_restart: 2,
                max_linear_cycles: 3,
                forcing_maximum: 1e-8,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 12,
        },
    )
}
fn initial_config() -> InitialSolveConfig {
    InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1e-12,
    }
}
fn adjoint_config() -> SecondOrderAdjointConfig {
    SecondOrderAdjointConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1e-12,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:.15e}, expected {expected:.15e}"
    );
}
fn difference<const N: usize>(
    mut point: [f64; N],
    component: usize,
    evaluate: impl Fn([f64; N]) -> f64,
) -> f64 {
    let center = point[component];
    let h = 2e-6 * (1.0 + center.abs());
    point[component] = center + h;
    let plus = evaluate(point);
    point[component] = center - h;
    (plus - evaluate(point)) / (2.0 * h)
}
fn objective(state: &SecondOrderState) -> f64 {
    TERMINAL.0[0] * state.q[0] + TERMINAL.1[0] * state.v[0] + TERMINAL.2[0] * state.a[0]
}
fn load(time: f64) -> f64 {
    0.25 + 0.1 * time
}

#[test]
fn nonlinear_initial_acceleration_and_all_pullbacks_match_physical_equilibrium() {
    let model = model();
    let input = [0.45, -0.65, 0.3];
    let seed = 0.8;
    let got = second_order_acceleration_vjp(
        &model,
        &[input[0]],
        &[input[1]],
        &[input[2]],
        &[seed],
        &Identity,
        &Identity,
        initial_config(),
        WORKSPACE,
        &mut || false,
    )
    .unwrap();
    let acceleration = |m: Damper, x: [f64; 3]| (x[2] - m.value(x[0], x[1])) / m.0[0];
    close(got.acceleration[0], acceleration(model, input), 1e-12);
    for (i, actual) in [got.initial_q[0], got.initial_v[0], got.forcing[0]]
        .into_iter()
        .enumerate()
    {
        close(
            actual,
            difference(input, i, |x| seed * acceleration(model, x)),
            1e-7,
        );
    }
    for i in 0..4 {
        close(
            got.parameters[i],
            difference(model.0, i, |p| seed * acceleration(Damper(p), input)),
            1e-7,
        );
    }
}
fn step_objective(model: Damper, initial: [f64; 3], force: f64) -> f64 {
    let mut state = SecondOrderState::new(0.0, &[initial[0]], &[initial[1]], &[initial[2]]);
    method().step(&mut state, &model, &[force]).unwrap();
    objective(&state)
}
#[test]
fn nonlinear_newton_stage_and_step_pullback_include_both_resistance_partials() {
    let model = model();
    let initial = [0.45, -0.65, 0.2];
    let force = 0.3;
    let start = SecondOrderState::new(0.0, &[initial[0]], &[initial[1]], &[initial[2]]);
    let got = method()
        .step_vjp(
            &start,
            &model,
            &[force],
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &Identity,
            adjoint_config(),
            WORKSPACE,
            &mut || false,
        )
        .unwrap();
    // Independently check the actual Chung-Hulbert equilibrium.
    let am = (2.0 * RHO - 1.0) / (RHO + 1.0);
    let af = RHO / (RHO + 1.0);
    let q = (1.0 - af) * got.q[0] + af * initial[0];
    let v = (1.0 - af) * got.v[0] + af * initial[1];
    let a = (1.0 - am) * got.a[0] + am * initial[2];
    close(model.0[0] * a + model.value(q, v), force, 1e-10);
    for (i, actual) in [got.initial_q[0], got.initial_v[0], got.initial_a[0]]
        .into_iter()
        .enumerate()
    {
        close(
            actual,
            difference(initial, i, |x| step_objective(model, x, force)),
            5e-6,
        );
    }
    for i in 0..4 {
        close(
            got.parameters[i],
            difference(model.0, i, |p| step_objective(Damper(p), initial, force)),
            5e-6,
        );
    }
    close(
        got.forcing[0],
        difference([force], 0, |f| step_objective(model, initial, f[0])),
        5e-6,
    );
}
fn consistent_initial(model: &Damper, qv: [f64; 2]) -> SecondOrderState {
    let acceleration = second_order_acceleration(
        model,
        &[qv[0]],
        &[qv[1]],
        &[load(0.0)],
        &Identity,
        initial_config(),
        WORKSPACE,
        &mut || false,
    )
    .unwrap();
    SecondOrderState::new(0.0, &[qv[0]], &[qv[1]], &acceleration.value)
}
fn trajectory_objective(model: Damper, qv: [f64; 2], steps: usize) -> f64 {
    let mut state = consistent_initial(&model, qv);
    for _ in 0..steps {
        let force = load(method().forcing_time(state.t).unwrap());
        method().step(&mut state, &model, &[force]).unwrap();
    }
    objective(&state)
}
#[test]
fn nonlinear_checkpoint_pullback_chains_consistent_initial_acceleration() {
    let model = model();
    let qv = [0.45, -0.65];
    let steps = 7;
    let initial = consistent_initial(&model, qv);
    let mut recording = RecordedStructural::new(
        method(),
        &model,
        &initial,
        StructuralRecordingConfig {
            steps,
            adjoint: adjoint_config(),
            max_workspace_components: WORKSPACE,
        },
    )
    .unwrap();
    assert_eq!(
        recording
            .advance(steps, steps, &mut || false)
            .unwrap()
            .status,
        StructuralRecordingStatus::ReachedEnd
    );
    close(
        objective(recording.state()),
        trajectory_objective(model, qv, steps),
        1e-12,
    );
    let reverse = recording
        .pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &[0.0; 4],
            &Identity,
            StructuralReplayBudget {
                checkpoints: 4,
                forward_steps: 64,
            },
            &mut || false,
        )
        .unwrap();
    let initial_reverse = second_order_acceleration_vjp(
        &model,
        &[qv[0]],
        &[qv[1]],
        &[load(0.0)],
        &reverse.initial_a,
        &Identity,
        &Identity,
        initial_config(),
        WORKSPACE,
        &mut || false,
    )
    .unwrap();
    let qv_gradient = [
        reverse.initial_q[0] + initial_reverse.initial_q[0],
        reverse.initial_v[0] + initial_reverse.initial_v[0],
    ];
    for i in 0..2 {
        close(
            qv_gradient[i],
            difference(qv, i, |x| trajectory_objective(model, x, steps)),
            5e-6,
        );
    }
    for i in 0..4 {
        close(
            reverse.parameters[i] + initial_reverse.parameters[i],
            difference(model.0, i, |p| trajectory_objective(Damper(p), qv, steps)),
            5e-6,
        );
    }
}
