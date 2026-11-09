//! Source-backed synthetic gases, actual public momentum coupling, no measured fit.
// Explicit paths: tests/port_schema_adoption.rs also loads this file through
// #[path], where a bare child would resolve beside this file, not under tests/.
#[path = "tests/real_gas.rs"]
pub(super) mod real_gas;

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn species(name: &str, mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: Provenance {
            source: "synthetic caloric gas".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn transport(name: &str) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic transport".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
pub(super) fn gas(x: f64, slope: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, slope),
                species("B", 0.004, 2.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        vec![transport("A"), transport("B")],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, 0.001).unwrap()
}
pub(super) fn controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.0007, 0.0015],
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        // The authored 2e-7 N force resolution over the 0.00005 s clock
        // (1e-10 N s on the former 0.0005 s step). An absolute impulse
        // tolerance that did not follow the clock accepted small heated
        // or mass-changed motions as exactly no motion.
        momentum_tolerance_n_s: 2e-7 * 0.00005,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
pub(super) fn drive(heat_j: f64) -> PistonInput {
    // Resolves even the light material-inertia specimens (about 380 rad/s)
    // below the unchanged 0.04-radian gas-compressibility admission.
    PistonInput {
        duration_s: 0.00005,
        heat_j,
        exterior_pressure_pa: 100000.0,
        external_force_n: 0.0,
    }
}

#[test]
fn equilibrium_is_bitwise_stationary_and_does_not_query_a_fictitious_volume() {
    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    let before = piston.clone();
    let report = piston.step(drive(0.0), controls(), || false).unwrap();
    assert_eq!(piston, before);
    assert_eq!(report.material_trials, 1);
    assert_eq!(report.total_energy_residual_j, 0.0);
    assert_eq!(report.displacement_m, 0.0);
}

#[test]
fn heat_drives_motion_and_work_is_not_counted_twice() {
    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    let initial = piston.gas().clone();
    let report = piston.step(drive(20.0), controls(), || false).unwrap();
    assert!(piston.velocity_m_s() > 0.0);
    assert!(piston.gas().volume_m3() > initial.volume_m3());
    assert!(piston.gas().state().state().temperature > initial.state().state().temperature);
    assert!(report.gas.work_by_gas_j > 0.0);
    assert!(
        (report.kinetic_energy_change_j - report.gas.work_by_gas_j
            + report.exterior_pressure_work_j)
            .abs()
            < 1e-7
    );
    assert!(report.total_energy_residual_j.abs() < 1e-7);
    let force = 0.01 * (report.midpoint.state().pressure - 100000.0);
    assert!(
        (piston.mass_kg() * piston.velocity_m_s() - drive(20.0).duration_s * force).abs() < 1e-10
    );
    assert!(
        (report.displacement_m - drive(20.0).duration_s * 0.5 * piston.velocity_m_s()).abs()
            < 1e-15
    );
    assert_eq!(piston.gas().mass_kg(), initial.mass_kg());
    assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
    assert_eq!(report.midpoint.components().len(), 2);
}

#[test]
fn changing_caloric_and_composition_data_changes_actual_piston_trajectory() {
    let run = |x, slope| {
        let mut piston = GasPiston::new(gas(x, slope), 1.0, 0.01, 0.0).unwrap();
        piston.step(drive(20.0), controls(), || false).unwrap();
        piston
    };
    let constant = run(0.75, 0.0);
    let varying = run(0.75, 0.001);
    let mixture = run(0.25, 0.001);
    assert!(varying.velocity_m_s() < constant.velocity_m_s());
    assert!(varying.gas().state().state().sound_speed < constant.gas().state().state().sound_speed);
    assert!((mixture.velocity_m_s() - varying.velocity_m_s()).abs() > 1e-5);
    assert_ne!(
        mixture.gas().state().state().dynamic_viscosity,
        varying.gas().state().state().dynamic_viscosity
    );
}

#[test]
fn coupled_trajectory_closes_total_energy_with_exterior_and_external_work() {
    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.1).unwrap();
    let initial = piston.clone();
    let mut heat = 0.0;
    let mut outside_work = 0.0;
    for i in 0..80 {
        let q = if i < 10 { 0.2 } else { 0.0 };
        let input = PistonInput {
            external_force_n: if i < 20 { -0.5 } else { 0.0 },
            ..drive(q)
        };
        let report = piston.step(input, controls(), || false).unwrap();
        heat += q;
        outside_work += report.external_force_work_j - report.exterior_pressure_work_j;
        assert!(report.momentum_residual_n_s.abs() <= controls().momentum_tolerance_n_s);
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
    let stored = piston.gas().mass_kg()
        * (piston.gas().state().specific_internal_energy_j_kg()
            - initial.gas().state().specific_internal_energy_j_kg());
    let kinetic =
        0.5 * piston.mass_kg() * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
    assert!((stored + kinetic - heat - outside_work).abs() < 80.0 * controls().energy_tolerance_j);
}

#[test]
fn acoustic_compressibility_predicts_the_small_signal_mechanical_period() {
    let mut piston = GasPiston::new(gas(0.75, 0.0), 1.0, 0.01, 0.001).unwrap();
    let medium = piston.gas().state().state();
    let omega = (medium.density * medium.sound_speed.powi(2) * piston.area_m2().powi(2)
        / (piston.mass_kg() * piston.gas().volume_m3()))
    .sqrt();
    let dt = core::f64::consts::TAU / omega / 256.0;
    for _ in 0..256 {
        piston
            .step(
                PistonInput {
                    duration_s: dt,
                    ..drive(0.0)
                },
                controls(),
                || false,
            )
            .unwrap();
    }
    assert!((piston.velocity_m_s() / 0.001 - 1.0).abs() < 1e-3);
}

#[test]
fn invalid_work_bounds_source_failures_and_late_cancellation_roll_back() {
    use fs_ports::PortedModule;

    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    // Keep an actual accepted window, including its lazily expanded transfer
    // cache. GasPiston's physical-state PartialEq deliberately omits ports,
    // so every refusal also compares the public schemas, transfers and audit.
    // This file is included by port_schema_adoption as an integration fixture:
    // the rollback check must use the public owner seam in both test contexts.
    piston.step(drive(0.0), controls(), || false).unwrap();
    assert!(piston.last_window_transfers().is_some());
    let before = piston.clone();
    let window = |owner: &GasPiston| {
        (
            owner.port_schemas().to_vec(),
            owner.port_coordinate_charts().to_vec(),
            owner.junction_schemas().to_vec(),
            owner.last_window_duration_s(),
            owner.last_window_transfers().map(<[_]>::to_vec),
            owner.last_window_transfer_error().cloned(),
            owner.audit_last_window().unwrap(),
        )
    };
    let before_window = window(&before);
    for input in [
        PistonInput {
            duration_s: 0.0,
            ..drive(1.0)
        },
        PistonInput {
            exterior_pressure_pa: -1.0,
            ..drive(1.0)
        },
        drive(1e9),
    ] {
        assert!(piston.step(input, controls(), || false).is_err());
        assert_eq!(piston, before);
        assert_eq!(window(&piston), before_window);
    }
    // The resolved 50 us / 20 J trial now converges in one correction.
    // At 200 us / 300 J the same gas requires another Newton correction,
    // without leaving its source domain or the unchanged omega*dt <= 0.04
    // clock guard. A matched full-budget trial proves the refusal reason.
    let nonlinear = PistonInput {
        duration_s: 0.0002,
        ..drive(300.0)
    };
    let medium = piston.gas().state().state();
    let omega = (medium.density * medium.sound_speed.powi(2) * piston.area_m2().powi(2)
        / (piston.mass_kg() * piston.gas().volume_m3()))
    .sqrt();
    assert!(omega * nonlinear.duration_s < 0.04);
    let mut resolved = before.clone();
    let accepted = resolved.step(nonlinear, controls(), || false).unwrap();
    assert!(accepted.iterations > 1);
    assert!(accepted.momentum_residual_n_s.abs() <= controls().momentum_tolerance_n_s);
    assert!(accepted.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert_ne!(window(&resolved), before_window);
    assert!(matches!(
        piston.step(
            nonlinear,
            PistonSolve {
                maximum_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(PistonError::NotConverged { iterations: 1, .. })
    ));
    assert_eq!(piston, before);
    assert_eq!(window(&piston), before_window);
    assert!(matches!(
        piston.step(
            drive(20.0),
            PistonSolve {
                volume_bracket_m3: [0.001 - 1e-10, 0.001 + 1e-10],
                ..controls()
            },
            || false
        ),
        Err(PistonError::OutsideBracket { .. })
    ));
    assert_eq!(piston, before);
    assert_eq!(window(&piston), before_window);
    let mut count = 0;
    let mut probe = piston.clone();
    probe
        .step(drive(20.0), controls(), || {
            count += 1;
            false
        })
        .unwrap();
    let mut actual = 0;
    assert!(matches!(
        piston.step(drive(20.0), controls(), || {
            actual += 1;
            actual == count
        }),
        Err(PistonError::Cancelled)
    ));
    assert_eq!(piston, before);
    assert_eq!(window(&piston), before_window);
}

#[cfg(test)]
#[path = "tests/adiabatic_solid_cavity_bridge"]
mod adiabatic_solid_cavity_bridge {
    mod equilibrium;
    use super::gas;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
    use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
    use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
    use fs_solid::hyper3::pressure::cavity::ClosedTetCavity;
    use fs_solid::linear3::TetAssemblyBudget;
    use std::sync::Arc;

    fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            f(&gate, &Cx::new(&gate, arena, StreamKey { seed: 101, kernel_id: 7, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic))
        })
    }
    fn entropy() -> GasEntropySolve {
        GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
            entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
    }
    fn shell() -> (Vec<[f64;3]>, Vec<[usize;4]>) {
        let nodes=vec![[0.,0.,0.],[3.,0.,0.],[0.,3.,0.],[0.,0.,3.],
            [0.5,0.5,0.5],[1.5,0.5,0.5],[0.5,1.5,0.5],[0.5,0.5,1.5]];
        let mut tets=Vec::new();
        for [a,b,c] in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]] {
            tets.push([a+4,b+4,c+4,c]); tets.push([a+4,b+4,b,c]); tets.push([a+4,a,b,c]);
        }
        (nodes,tets)
    }
    fn cavity(nodes:&[[f64;3]],tets:&[[usize;4]],cx:&Cx<'_>)->ClosedTetCavity {
        let surface=Arc::new(ClosedTetPressureSurface::new(nodes,tets,TetAssemblyBudget::standard(),
            4*tets.len(),cx).unwrap());
        let inner:Vec<_>=surface.faces().iter().enumerate()
            .filter_map(|(i,f)| f.iter().all(|&n| n>=4).then_some(i)).collect();
        ClosedTetCavity::new(surface,&inner,cx).unwrap()
    }
    fn close(a:f64,b:f64,tol:f64){ assert!((a-b).abs()<=tol*a.abs().max(b.abs()).max(1.0),"{a} != {b}"); }

    #[test]
    fn source_gas_pressure_energy_and_acoustic_tangent_follow_the_actual_cavity() {
        with_cx(|_,cx| {
            let (nodes,tets)=shell(); let cavity=cavity(&nodes,&tets,cx);
            let seed=gas(0.75,0.001);
            let gas=ClosedGasParcel::new(seed.model().clone(),seed.point().clone(),cavity.reference_volume_m3()).unwrap();
            let exterior=30_000.0;
            let u:Vec<f64>=nodes.iter().flat_map(|x| x.map(|v| -0.04*v)).collect();
            let load=cavity.adiabatic_gas_load(&u,&gas,exterior,entropy(),cx).unwrap();
            assert!(load.volume_m3 < cavity.reference_volume_m3());
            assert!(load.gas.state().state().pressure > gas.state().state().pressure);
            assert!(load.potential_change_j > 0.0);
            let face=cavity.face_pressures_pa(load.net_pressure_pa,cx).unwrap();
            let expected=cavity.surface().force_n(&u,TetPressure::PerFace(&face),cx).unwrap();
            for (&a,&b) in load.force_n.iter().zip(&expected) { close(a,b,2e-11); }
            let mut d=vec![0.0;u.len()];
            for (i,&node) in cavity.nodes().iter().enumerate() {
                d[3*node]=0.02*(i as f64-1.5); d[3*node+1]=-0.01*(i as f64+1.0); d[3*node+2]=0.015;
            }
            let eps=1e-6; let mut plus=u.clone(); let mut minus=u.clone();
            for i in 0..u.len(){ plus[i]+=eps*d[i]; minus[i]-=eps*d[i]; }
            let fp=cavity.adiabatic_gas_load(&plus,&gas,exterior,entropy(),cx).unwrap().force_n;
            let fm=cavity.adiabatic_gas_load(&minus,&gas,exterior,entropy(),cx).unwrap().force_n;
            for row in 0..u.len() {
                let applied:f64=(0..u.len()).map(|col| load.force_tangent.get(row,col)*d[col]).sum();
                close(applied,(fp[row]-fm[row])/(2.0*eps),3e-6);
            }
            let before=gas.clone();
            assert!(cavity.adiabatic_gas_load(&u,&gas,-1.0,entropy(),cx).is_err());
            assert_eq!(gas,before);
        });
    }

    /// Uniform 1-to-8 tetrahedral refinement with globally shared edge midpoints.
    /// Inner-boundary ancestry is retained independently of generated face order.
    fn refined_shell(levels: usize) -> (Vec<[f64; 3]>, Vec<[usize; 4]>, Vec<bool>) {
        let (mut nodes, mut tets) = shell();
        let mut inner: Vec<_> = (0..nodes.len()).map(|i| i >= 4).collect();
        for _ in 0..levels {
            let mut edges = std::collections::BTreeMap::new();
            let mut refined = Vec::with_capacity(8 * tets.len());
            for &[a, b, c, d] in &tets {
                let mut middle = |left: usize, right: usize| {
                    let key = (left.min(right), left.max(right));
                    *edges.entry(key).or_insert_with(|| {
                        let point = std::array::from_fn(|axis|
                            0.5 * nodes[left][axis] + 0.5 * nodes[right][axis]);
                        let index = nodes.len();
                        let is_inner = inner[left] && inner[right];
                        nodes.push(point);
                        inner.push(is_inner);
                        index
                    })
                };
                let ab = middle(a, b);
                let ac = middle(a, c);
                let ad = middle(a, d);
                let bc = middle(b, c);
                let bd = middle(b, d);
                let cd = middle(c, d);
                refined.extend_from_slice(&[
                    [a, ab, ac, ad], [b, ab, bc, bd],
                    [c, ac, bc, cd], [d, ad, bd, cd],
                    [ab, cd, ac, ad], [ab, cd, ad, bd],
                    [ab, cd, bd, bc], [ab, cd, bc, ac],
                ]);
            }
            tets = refined;
        }
        (nodes, tets, inner)
    }

    #[test]
    fn g1_refined_cavity_retains_global_acoustic_stiffness_without_dense_assembly() {
        with_cx(|_, cx| {
            let (nodes, tets, inner_nodes) = refined_shell(3);
            let surface = Arc::new(ClosedTetPressureSurface::new(
                &nodes, &tets, TetAssemblyBudget::standard(), 4 * tets.len(), cx,
            ).unwrap());
            let faces: Vec<_> = surface.faces().iter().enumerate()
                .filter_map(|(i, face)| face.iter().all(|&node| inner_nodes[node]).then_some(i))
                .collect();
            let cavity = ClosedTetCavity::new(surface, &faces, cx).unwrap();
            assert_eq!(faces.len(), 256);
            assert!(3 * cavity.nodes().len() > 256);
            let seed = gas(0.75, 0.001);
            let original = ClosedGasParcel::new(
                seed.model().clone(), seed.point().clone(), cavity.reference_volume_m3(),
            ).unwrap();
            let u: Vec<_> = nodes.iter().flat_map(|x| x.map(|value| -0.02 * value)).collect();
            let direction: Vec<_> = nodes.iter().flat_map(|x| x.map(|value| 0.01 * value)).collect();
            let load = cavity.adiabatic_gas_load(&u, &original, 30_000.0, entropy(), cx).unwrap();
            let jacobian = &load.force_tangent;
            assert_eq!(jacobian.nrows(), u.len());
            assert_eq!(jacobian.ncols(), u.len());
            assert_eq!(jacobian.acoustic_rank_one_vector().len(), u.len());
            assert!(jacobian.geometric_part().nnz() <= 81 * cavity.surface().faces().len());
            let action = jacobian.apply(&direction, cx).unwrap();
            let mut output = vec![7.0; u.len()];
            let mut scratch = vec![0.0; u.len()];
            jacobian.apply_into(&direction, &mut output, &mut scratch, cx).unwrap();
            assert_eq!(action, output);

            // Dilation gives an independent closed-form geometry derivative:
            // dV = 3 V ds/s and d(grad V) = 2 grad(V) ds/s.
            let geometry = cavity.evaluate(&u, cx).unwrap();
            let response = load.gas.pressure_response().unwrap();
            let coefficient = 0.01 / 0.98 * (
                2.0 * load.net_pressure_pa
                    + 3.0 * load.volume_m3 * response.pressure_volume_derivative_pa_m3
            );
            for (local, &node) in cavity.nodes().iter().enumerate() {
                for axis in 0..3 {
                    close(action[3 * node + axis],
                        coefficient * geometry.gradient_m2()[3 * local + axis], 3e-10);
                }
            }

            // Distant boundary coordinates still couple through one conserved
            // gas inventory, even when their geometric sparse entry is zero.
            let weights = jacobian.acoustic_rank_one_vector();
            let populated: Vec<_> = weights.iter().enumerate()
                .filter_map(|(i, w)| (*w != 0.0).then_some(i)).collect();
            let mut nonlocal = false;
            'pairs: for &row in &populated {
                for &column in &populated {
                    if jacobian.geometric_part().get(row, column) == 0.0
                        && weights[row] * weights[column] != 0.0
                    {
                        close(jacobian.get(row, column), -weights[row] * weights[column], 1e-14);
                        nonlocal = true;
                        break 'pairs;
                    }
                }
            }
            assert!(nonlocal);

            let epsilon = 1e-4;
            let plus: Vec<_> = u.iter().zip(&direction).map(|(x, d)| x + epsilon * d).collect();
            let minus: Vec<_> = u.iter().zip(&direction).map(|(x, d)| x - epsilon * d).collect();
            let positive = cavity.adiabatic_gas_load(&plus, &original, 30_000.0, entropy(), cx).unwrap();
            let negative = cavity.adiabatic_gas_load(&minus, &original, 30_000.0, entropy(), cx).unwrap();
            for i in 0..u.len() {
                close(action[i], (positive.force_n[i] - negative.force_n[i]) / (2.0 * epsilon), 3e-6);
            }
            let mechanical_power: f64 = load.force_n.iter().zip(&direction).map(|(f, d)| f * d).sum();
            close(mechanical_power,
                -(positive.potential_change_j - negative.potential_change_j) / (2.0 * epsilon), 3e-6);
            assert_eq!(original.volume_m3(), cavity.reference_volume_m3());
        });
    }

    #[test]
    fn g4_cavity_tangent_output_is_atomic_on_bad_directions_and_cancellation() {
        with_cx(|gate, cx| {
            let (nodes, tets) = shell();
            let cavity = cavity(&nodes, &tets, cx);
            let seed = gas(0.75, 0.001);
            let original = ClosedGasParcel::new(
                seed.model().clone(), seed.point().clone(), cavity.reference_volume_m3(),
            ).unwrap();
            let load = cavity.adiabatic_gas_load(&[0.0; 24], &original, 100_000.0, entropy(), cx).unwrap();
            // Zero net force does not remove acoustic restoring stiffness.
            assert_eq!(load.force_n, vec![0.0; 24]);
            assert!(load.force_tangent.acoustic_rank_one_vector().iter().any(|&x| x != 0.0));
            let mut output = vec![71.0; 24];
            let mut scratch = vec![0.0; 24];
            assert!(load.force_tangent.apply_into(&[0.0; 23], &mut output, &mut scratch, cx).is_err());
            assert_eq!(output, vec![71.0; 24]);
            let mut bad = [0.0; 24];
            bad[23] = f64::NAN;
            assert!(load.force_tangent.apply_into(&bad, &mut output, &mut scratch, cx).is_err());
            assert_eq!(output, vec![71.0; 24]);
            gate.request();
            assert!(matches!(
                load.force_tangent.apply_into(&[0.0; 24], &mut output, &mut scratch, cx),
                Err(fs_solid::hyper3::HyperTetError::Cancelled)
            ));
            assert_eq!(output, vec![71.0; 24]);
        });
    }
}

