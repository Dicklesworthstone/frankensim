//! Two independently supplied air streams exchanging heat through a real FEM
//! slab. Uses the actual fan operating point, convection card, and shared-solid
//! conjugate driver. An explicit leakage branch bypasses the thermal surfaces.
//!
//! Run from the workspace root:
//! `cargo run -p fs-airflow --example branched_cooling -- 290 330`
//!
//! Arguments are the two inlet temperatures in kelvin. The geometry, fan curve,
//! and material below are illustrative declarations, not validated product data.
//! Air transport properties and hydraulic flows stay frozen during the exchange.

use std::error::Error;

use fs_airflow::conjugate::{AirPath, AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::thermal::solve_conjugate_branches;
use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint,
    LeakageElement, LossElement, LossNetwork, LossResistance, SourceProvenance,
    ToleranceBasis, solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, InitialGuess, SolveConfig, solve};
use fs_convection::{CorrelationId, ThermalConductivity, evaluate};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Area, Density, DynViscosity, Length, Pressure, VolumetricFlowRate};

const LENGTH: f64 = 0.02;
const AREA: f64 = 0.01;
const CONDUCTIVITY: f64 = 2.0;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let inlets = match args.as_slice() {
        [] => [290.0, 330.0],
        [left, right] => [left.parse::<f64>()?, right.parse::<f64>()?],
        _ => return Err("usage: branched_cooling [left-inlet-K right-inlet-K]".into()),
    };
    if inlets.iter().any(|t| !t.is_finite() || *t <= 0.0) {
        return Err("inlet temperatures must be finite and positive kelvin".into());
    }
    let source = SourceProvenance::new("illustrative analytic fan/loss data", "branched-cooling-v1");
    let fan = FanBank::new(FanCurve::new("illustrative", vec![
        FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(120.0)),
        FanPoint::new(VolumetricFlowRate::new(0.002), Pressure::new(0.0)),
    ], source.clone(), 0.0, ToleranceBasis::Analytic,
        VolumetricFlowRate::new(1.0e-8), (0.5, 2.0))?, 1, FanArrangement::Series, 1.0)?;
    let loss = |name, resistance| LossElement::new(name, LossResistance::new(resistance),
        0.0, source.clone(), ToleranceBasis::Analytic);
    let network = EnclosureNetwork::new(LossNetwork::parallel(vec![
        LossNetwork::Element(loss("left", 4.0e8)?),
        LossNetwork::Element(loss("right", 1.0e8)?),
    ])?, LeakageElement::new(loss("unheated-bypass", 1.0e10)?));
    let operating = solve_operating_point(&fan, &network)?;
    let card = CorrelationId::ALL.iter().copied()
        .find(|id| id.name() == "convection.circular-duct-hausen-developing")
        .ok_or("the developing laminar convection card is missing")?;
    let mut paths = Vec::new();
    let mut coefficients = Vec::new();
    for (index, name) in ["left", "right"].iter().enumerate() {
        let handoff = operating.correlation_handoff(name, Area::new(0.001),
            Density::new(1.2), DynViscosity::new(1.846e-5), Length::new(0.01), 0.707)?;
        let nu = evaluate(card, handoff.correlation_inputs.with_length_ratio(10.0))?;
        if !nu.evidence().model.in_domain {
            return Err("convection card is outside its declared domain".into());
        }
        let h = nu.heat_transfer_coefficient(ThermalConductivity::new(26.3e-3),
            Length::new(0.01))?.value.value();
        coefficients.push(h);
        paths.push(AirPath::new(inlets[index], 1.2 * handoff.branch_flow.value.value(),
            1007.0, vec![AirSegment::new(name, AREA, h)?])?);
    }
    let (complex, positions) = box_grid([2, 1, 1], [LENGTH, 0.1, 0.1]);
    let mesh = ConductionMesh::new(complex, positions)?;
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY)?;
    let source = ScalarField::Uniform(0.0);
    let gate = CancelGate::new();
    let result = ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<_, Box<dyn Error>> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 11, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let mut solid_error: Option<Box<dyn Error>> = None;
        let result = solve_conjugate_branches(&cx, &paths, &ConjugateConfig::default(), |cx, refs| {
            let solid = (|| -> Result<_, Box<dyn Error>> {
                let boundary = ThermalBoundaryBuilder::new(&mesh)
                    .region("left", |face| on_box_face(face.centroid[0], 0.0),
                        ThermalBc::robin(coefficients[0], refs[0])?)?
                    .region("right", |face| on_box_face(face.centroid[0], LENGTH),
                        ThermalBc::robin(coefficients[1], refs[1])?)?
                    .adiabatic_remainder().finish()?;
                let mut config = SolveConfig::default();
                config.initial = InitialGuess::Uniform(0.5 * inlets[0] + 0.5 * inlets[1]);
                config.linear.tolerance = 1.0e-13;
                let field = solve(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
                    material: &material, element_materials: None, source: &source }, config)?;
                ["left", "right"].iter().map(|name| {
                    field.report.robin_fluxes.iter().find(|flux| flux.region == *name)
                        .map(SolidRegionState::from_robin_flux)
                        .ok_or_else(|| "a declared Robin region is absent from the solid report".into())
                }).collect::<Result<Vec<_>, Box<dyn Error>>>()
            })();
            match solid {
                Ok(states) => Ok(states),
                Err(error) => {
                    solid_error = Some(error);
                    // Stop the iteration; report the original solid error below.
                    Err(AirflowError::Cancelled { iteration: 0, references_k: refs.to_vec() })
                }
            }
        });
        if let Some(error) = solid_error { return Err(error); }
        Ok(result?)
    })?;
    let effective: Vec<f64> = paths.iter().zip(&coefficients).map(|(path, &h)| {
        let capacity = path.capacity_rate_w_per_k();
        capacity * (1.0 - (-h * AREA / capacity).exp())
    }).collect();
    let q = (inlets[1] - inlets[0])
        / (1.0 / effective[0] + LENGTH / (CONDUCTIVITY * AREA) + 1.0 / effective[1]);
    println!("nominal fan flow: {:.9e} m^3/s; shared solid solves: {}",
        operating.flow.value.value(), result.iterations);
    for (index, branch) in result.branches.iter().enumerate() {
        println!("{}: inlet {:.6} K, outlet {:.6} K, wall {:.6} K, heat into air {:+.9e} W, branch imbalance {:+.3e} W",
            branch.solid[0].region, inlets[index], branch.march.outlet_temperature_k,
            branch.solid[0].mean_wall_temperature_k, branch.march.total_heat_rate_w,
            branch.balance.interface_imbalance_w);
    }
    println!("closed-form heat into left stream: {q:+.9e} W; numerical difference: {:+.3e} W",
        result.branches[0].march.total_heat_rate_w - q);
    println!("Nominal frozen-property model, not an uncertainty certificate or experimentally validated design.");
    Ok(())
}
