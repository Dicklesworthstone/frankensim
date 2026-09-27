//! Adaptive FEM solid + finite-capacity well-mixed enclosure air.
//!
//! Run `cargo run -p fs-airflow --example stored_air -- --air-capacity-j-k 0.02`.
//! The 10 mm cube uses the actual P1 conduction/capacity operator, not a lumped
//! solid. Its entire exterior exchanges with one isothermal air inventory.
//! A 2 W pulse ends exactly at 1 s; the declared run ends at 5 s. The outlet
//! is at the air temperature and the inlet is fixed at 300 K. Ventilation is
//! a declared heat-capacity rate, not a fan/momentum solve. No CFD, mixing-law
//! validation, spatial error bound, or native .fsim/ledger integration is claimed.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder,
    fixtures::box_grid, transient::{VolumetricHeatCapacity, backward_euler::{BackwardEuler, StepConfig}}};
use fs_couple::iqn_ils::{IqnIlsConfig, driver::{BalanceControl, CouplingControls, CouplingMethod,
    CouplingTrial, InterfaceControl, StepInterval, march::adaptive::{AdaptiveEvolution, AdaptiveSettings}}};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const RHO_CP: f64 = 2.0e6;
const INITIAL_K: f64 = 300.0;
const HTC: f64 = 100.0;

#[derive(Clone, Debug, PartialEq)]
struct State { solid_k: Vec<f64>, air_k: f64, net_input_j: f64 }

#[derive(Clone, Copy, Debug)]
struct Inputs { air_capacity_j_k: f64, ventilation_w_k: f64, tolerance_k: f64, attempts: usize }
impl Default for Inputs {
    fn default() -> Self {
        Self { air_capacity_j_k: 0.02, ventilation_w_k: 0.01, tolerance_k: 1.0e-4, attempts: 4096 }
    }
}
impl Inputs {
    fn validate(self) -> Result<Self, String> {
        if !(self.air_capacity_j_k.is_finite() && self.air_capacity_j_k > 0.0
            && self.ventilation_w_k.is_finite() && self.ventilation_w_k >= 0.0
            && self.tolerance_k.is_finite() && self.tolerance_k > 0.0
            && self.attempts > 0 && self.attempts <= 20_000)
        { return Err("positive finite air capacity/tolerance and 1..=20000 attempts required; ventilation must be finite and nonnegative".into()); }
        Ok(self)
    }
}

fn parse(args: &[String]) -> Result<Inputs, String> {
    let mut inputs = Inputs::default();
    let mut seen = std::collections::BTreeSet::new();
    if args.len() % 2 != 0 { return Err("every flag requires a value".into()); }
    for pair in args.chunks_exact(2) {
        if !seen.insert(pair[0].as_str()) { return Err(format!("duplicate {}", pair[0])); }
        match pair[0].as_str() {
            "--air-capacity-j-k" => inputs.air_capacity_j_k = pair[1].parse().map_err(|_| "invalid air capacity")?,
            "--ventilation-w-k" => inputs.ventilation_w_k = pair[1].parse().map_err(|_| "invalid ventilation")?,
            "--tolerance-k" => inputs.tolerance_k = pair[1].parse().map_err(|_| "invalid tolerance")?,
            "--attempts" => inputs.attempts = pair[1].parse().map_err(|_| "invalid attempts")?,
            other => return Err(format!("unknown option {other}")),
        }
    }
    inputs.validate()
}

struct Model<'a> {
    mesh: &'a ConductionMesh,
    solid: BackwardEuler<'a>,
    material: ConductivityModel,
    capacity_j_k: f64,
    ventilation_w_k: f64,
    conductance_w_k: f64,
}
impl<'a> Model<'a> {
    fn new(cx: &Cx<'_>, mesh: &'a ConductionMesh, inputs: Inputs) -> Result<Self, ConductionError> {
        let solid = BackwardEuler::uniform(cx, mesh, VolumetricHeatCapacity::declared(RHO_CP)?)?;
        let conductance_w_k = HTC * mesh.boundary().iter().map(|f| f.area).sum::<f64>();
        Ok(Self { mesh, solid, material: ConductivityModel::isotropic_declared(1.5)?,
            capacity_j_k: inputs.air_capacity_j_k, ventilation_w_k: inputs.ventilation_w_k,
            conductance_w_k })
    }

    fn trial(&self, cx: &Cx<'_>, old: &State, interval: StepInterval, x: &[f64])
        -> Result<CouplingTrial<State>, ConductionError>
    {
        if x.len() != 1 || !(x[0].is_finite() && x[0] > 0.0) {
            return Err(ConductionError::Config { parameter: "stored-air reference",
                what: "one finite positive kelvin value is required".into() });
        }
        let boundary = ThermalBoundaryBuilder::new(self.mesh)
            .region("enclosure-air", |_| true, ThermalBc::robin(HTC, x[0])?)?.finish()?;
        let source = ScalarField::Uniform(if interval.start_s < 1.0 { 2.0e6 } else { 0.0 });
        let problem = ConductionProblem { mesh: self.mesh, boundary: &boundary,
            material: &self.material, element_materials: None, source: &source };
        let dt = interval.duration_s();
        let response = self.solid.advance(cx, problem, None, &old.solid_k, dt, StepConfig {
            linear: LinearConfig { tolerance: 1.0e-12, max_iterations: 512, restart: 60 }, energy_tolerance_j: 1.0e-8,
        })?;
        // Solve the air endpoint implicitly using the actual FEM wall exchange.
        // Q_solid(x) + H*(x - T_air_new) is the same wall field's exchange at
        // the air endpoint; no fitted wall temperature or independent flux is used.
        let delta = dt * (response.robin_out_w
            + self.conductance_w_k * (x[0] - old.air_k)
            - self.ventilation_w_k * (old.air_k - INITIAL_K))
            / (self.capacity_j_k + dt * (self.conductance_w_k + self.ventilation_w_k));
        let air_k = old.air_k + delta;
        let air_storage_j = self.capacity_j_k * (air_k - old.air_k);
        let ventilation_out_w = self.ventilation_w_k * (air_k - INITIAL_K);
        let net_input_j = dt * (response.source_w - ventilation_out_w);
        let energy_residual_j = response.stored_energy_change_j + air_storage_j - net_input_j;
        let exchange_residual_j = dt * self.conductance_w_k * (air_k - x[0]);
        if !(air_k.is_finite() && air_k > 0.0 && energy_residual_j.is_finite()
            && exchange_residual_j.is_finite() && net_input_j.is_finite()
            && (old.net_input_j + net_input_j).is_finite())
        { return Err(ConductionError::Config { parameter: "stored-air endpoint",
            what: "temperature, exchange or storage is not representable".into() }); }
        Ok(CouplingTrial { state: State { solid_k: response.temperature, air_k,
            net_input_j: old.net_input_j + net_input_j }, image: vec![air_k],
            balance_residuals: vec![energy_residual_j, exchange_residual_j] })
    }

    fn stored_energy_j(&self, state: &State) -> f64 {
        // Independent volume integral of the actual published P1 field.
        let mut energy = self.capacity_j_k * (state.air_k - INITIAL_K);
        for (e, tet) in self.mesh.complex().tets.iter().enumerate() {
            let sum = tet.iter().map(|&v| state.solid_k[v as usize] - INITIAL_K).sum::<f64>();
            energy += RHO_CP * self.mesh.element_volume(e) * sum / 4.0;
        }
        energy
    }
}

fn evolution(mesh: &ConductionMesh, end_s: f64) -> AdaptiveEvolution<State> {
    let endpoints = if end_s > 1.0 { vec![1.0, end_s] } else { vec![end_s] };
    AdaptiveEvolution::new(State { solid_k: vec![INITIAL_K; mesh.vertex_count()],
        air_k: INITIAL_K, net_input_j: 0.0 }, vec![INITIAL_K], 0.0, endpoints,
        CouplingControls { max_evaluations: 32, relaxation: 0.5,
            method: CouplingMethod::IqnIls(IqnIlsConfig::default()),
            interfaces: vec![InterfaceControl { scale: 1.0, absolute_tolerance: 1.0e-8, relative_tolerance: 0.0 }],
            balances: ["solid-plus-air-energy-j", "wall-exchange-j"].into_iter().map(|name|
                BalanceControl { name: name.into(), absolute_tolerance: 1.0e-7 }).collect(),
        }, AdaptiveSettings { initial_step_s: 0.5, minimum_step_s: 1.0e-5,
            maximum_step_s: 1.0, method_order: 1 }).expect("fixed admitted example controls")
}

fn distance(coarse: &State, fine: &State, tolerance_k: f64) -> f64 {
    coarse.solid_k.iter().zip(&fine.solid_k).fold((coarse.air_k - fine.air_k).abs(),
        |maximum, (a, b)| maximum.max((a - b).abs())) / tolerance_k
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let gate = CancelGate::new();
        let cx = Cx::new(&gate, arena, StreamKey { seed: 7, kernel_id: 919, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}
fn mesh(n: usize) -> Result<ConductionMesh, ConductionError> {
    let (complex, positions) = box_grid([n, n, n], [0.01, 0.01, 0.01]);
    ConductionMesh::new(complex, positions)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inputs = parse(&std::env::args().skip(1).collect::<Vec<_>>())?;
    with_cx(|cx| -> Result<(), Box<dyn std::error::Error>> {
        let mesh = mesh(3)?;
        let model = Model::new(cx, &mesh, inputs)?;
        let mut run = evolution(&mesh, 5.0);
        let (mut evaluations, mut rejected) = (0usize, 0usize);
        eprintln!("FEM-solid/well-mixed-air estimate; air_capacity_j_k={} ventilation_w_k={} tolerance_k={}; 2W pulse [0,1]s; no spatial/physical-validation claim",
            inputs.air_capacity_j_k, inputs.ventilation_w_k, inputs.tolerance_k);
        println!("time_s,solid_max_k,air_k,local_error_ratio,stored_minus_net_input_j");
        for _ in 0..inputs.attempts {
            let report = run.advance(1, &mut |old, interval, x| model.trial(cx, old, interval, x),
                &mut |_, coarse, fine, _| Ok(distance(coarse, fine, inputs.tolerance_k)),
                &mut || cx.checkpoint().is_err())?;
            evaluations += report.evaluations; rejected += report.rejected;
            if let Some(row) = report.accepted.last() {
                let state = run.state();
                let maximum = state.solid_k.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                println!("{:.17},{:.17},{:.17},{:.17e},{:.17e}", run.time_s(), maximum, state.air_k,
                    row.error_ratio, model.stored_energy_j(state) - state.net_input_j);
            }
            if report.complete { break; }
        }
        eprintln!("complete={} accepted={} rejected={} producer_evaluations={} retained_time_s={}",
            run.is_complete(), run.accepted_steps(), rejected, evaluations, run.time_s());
        if !run.is_complete() { return Err("attempt budget exhausted; CSV is an accepted partial trajectory".into()); }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn solve(cx: &Cx<'_>, model: &Model<'_>, tolerance: f64, end: f64) -> AdaptiveEvolution<State> {
        let mut run = evolution(model.mesh, end);
        let report = run.advance(20_000, &mut |old, interval, x| model.trial(cx, old, interval, x),
            &mut |_, a, b, _| Ok(distance(a, b, tolerance)), &mut || false).unwrap();
        assert!(report.complete); run
    }
    #[test]
    fn closed_enclosure_balances_actual_field_energy_and_retains_spatial_gradients() {
        with_cx(|cx| {
            let mesh = mesh(2).unwrap();
            let inputs = Inputs { ventilation_w_k: 0.0, ..Inputs::default() };
            let model = Model::new(cx, &mesh, inputs).unwrap();
            let run = solve(cx, &model, 1e-4, 2.0);
            let state = run.state();
            assert!((state.net_input_j - 2.0).abs() < 1e-12);
            assert!((model.stored_energy_j(state) - 2.0).abs() < 2e-5);
            let min = state.solid_k.iter().copied().fold(f64::INFINITY, f64::min);
            let max = state.solid_k.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(max - min > 1e-6, "solid was reduced to an isothermal substitute");
            assert!(state.air_k > INITIAL_K && state.air_k < max);
        });
    }
    #[test]
    fn declared_air_capacity_changes_actual_fem_coupled_heatup() {
        with_cx(|cx| {
            let mesh = mesh(2).unwrap();
            let light = Model::new(cx, &mesh, Inputs { ventilation_w_k: 0.0, ..Inputs::default() }).unwrap();
            let heavy = Model::new(cx, &mesh, Inputs { air_capacity_j_k: 2.0, ventilation_w_k: 0.0, ..Inputs::default() }).unwrap();
            let a = solve(cx, &light, 1e-4, 0.5);
            let b = solve(cx, &heavy, 1e-4, 0.5);
            assert!(a.state().air_k > b.state().air_k + 1e-3);
            assert!((light.stored_energy_j(a.state()) - 1.0).abs() < 2e-5);
            assert!((heavy.stored_energy_j(b.state()) - 1.0).abs() < 2e-5);
        });
    }
    #[test]
    fn actual_fem_attempt_checkpoint_reproduces_the_uninterrupted_field() {
        with_cx(|cx| {
            let mesh = mesh(1).unwrap();
            let model = Model::new(cx, &mesh, Inputs::default()).unwrap();
            let full = solve(cx, &model, 1e-4, 0.5);
            let mut prefix = evolution(&mesh, 0.5);
            prefix.advance(1, &mut |old, interval, x| model.trial(cx, old, interval, x),
                &mut |_, a, b, _| Ok(distance(a, b, 1e-4)), &mut || false).unwrap();
            let mut resumed = prefix.clone();
            for _ in 0..20_000 {
                if resumed.advance(3, &mut |old, interval, x| model.trial(cx, old, interval, x),
                    &mut |_, a, b, _| Ok(distance(a, b, 1e-4)), &mut || false).unwrap().complete { break; }
            }
            assert_eq!(full, resumed);
        });
    }
    #[test]
    fn command_refuses_unknown_duplicate_and_nonphysical_inputs() {
        for args in [vec!["--unknown", "1"], vec!["--air-capacity-j-k", "0"],
            vec!["--ventilation-w-k", "-1"], vec!["--tolerance-k", "NaN"],
            vec!["--attempts", "0"], vec!["--attempts", "1", "--attempts", "2"]] {
            assert!(parse(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
    }
}
