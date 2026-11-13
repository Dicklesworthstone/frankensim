//! Original causal layer heat, native caloric cells, and transactional handoff.
//! These are authored numerical sources, not empirical material validation.
mod laminate_thermal_common;

use fs_couple::laminate_thermal::harmonic::{
    CausalLayeredShell, HarmonicLayerSource, LayeredShellHeatError,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, PropertyKey, PropertyValue};
use fs_material::laminate::ReferenceSurface;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection,
};
use fs_material::thermal::transport::stack::transient::{
    BoundLayeredThermalBody,
    implicit::{LayeredThermalSolveControls, LayeredThermalWorkspace},
};
use fs_material::visco_tensor::{
    kelvin_control::TensorKelvinStressControlBudget,
    source::{
        ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
    },
    stress_control::{TensorPlaneStressHarmonicBudget, TensorStressControlBudget},
};
use fs_math::{c64::C64, dd::Dd};
use fs_plate::shell::{
    ShellSupport,
    viscosity::{
        ShellViscosityBudget,
        harmonic::{TensorShellMaterialBudget, TensorShellMaterialCausalBudget},
    },
};
use fs_qty::{Density, Dims, Pressure, Time};
use laminate_thermal_common as common;
use std::sync::Arc;

const FREQUENCIES: [f64; 2] = [4.0, 7.0];
const AMPLITUDES: [f64; 2] = [0.003, 0.001];
const FRACTIONS: [f64; 2] = [0.4, 0.7];
const TAUS: [f64; 2] = [0.2, 0.65];
const DT: f64 = 0.125;
const ENERGY_TOLERANCE: f64 = 1e-7;

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 317,
                    kernel_id: 317,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn properties() -> TensorPronyProperties {
    TensorPronyProperties {
        elasticity: TensorElasticity::Isotropic,
        strain_limit: "strain_limit".into(),
        branches: vec![TensorBranchProperties {
            stiffness_upper: core::array::from_fn(|i| format!("branch_c{i:02}")),
            relaxation_time: "branch_tau".into(),
        }],
    }
}
fn source_card(layer: usize, varying_tau: bool) -> Arc<MaterialCard> {
    let base = common::card_with_coefficients(
        &format!("causal harmonic layer {layer}"),
        0.0,
        0.0,
        0.0,
        common::CP,
    );
    let template = base.claims_for("young_modulus")[0].1.clone();
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let replacement = match claim.key.name() {
            DENSITY_PROPERTY => Some((Density::DIMS, common::RHO)),
            "young_modulus" => Some((Pressure::DIMS, common::E)),
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY => Some((INVERSE_TEMPERATURE_DIMS, 0.0)),
            _ => None,
        };
        if let Some((dims, value)) = replacement {
            claim.value = PropertyValue::Scalar { value, dims };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let branch = &properties().branches[0];
    let mut authored = vec![
        ("strain_limit".to_owned(), Dims::NONE, 0.02),
        (branch.relaxation_time.clone(), Time::DIMS, TAUS[layer]),
    ];
    let mut entry = 0;
    for row in 0..6 {
        for column in row..6 {
            authored.push((
                branch.stiffness_upper[entry].clone(),
                Pressure::DIMS,
                if row == column {
                    FRACTIONS[layer] * common::E
                } else {
                    0.0
                },
            ));
            entry += 1;
        }
    }
    for (name, dims, value) in authored {
        let mut claim = template.clone();
        claim.key = PropertyKey::new(&name, dims);
        claim.value = PropertyValue::Scalar { value, dims };
        claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        if varying_tau && name == "branch_tau" {
            // Even equal endpoint samples retain an authored temperature axis.
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(280.0, value), (360.0, value)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    Arc::new(MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap())
}
fn sources(cards: &[Arc<MaterialCard>]) -> Vec<HarmonicLayerSource> {
    cards
        .iter()
        .enumerate()
        .map(|(layer, card)| HarmonicLayerSource {
            card: Arc::clone(card),
            material: ResolvedTensorMaterial::resolve(
                card,
                &common::at(300.0),
                properties(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
            orientation_rad: [0.31, -0.47][layer],
        })
        .collect()
}
fn owner<'b, 'a>(
    bound: &'b BoundLayeredThermalBody<'a>,
    cards: &[Arc<MaterialCard>],
    cx: &Cx<'_>,
) -> Result<CausalLayeredShell<'b, 'a>, LayeredShellHeatError> {
    CausalLayeredShell::new(
        bound,
        sources(cards),
        common::discretization().mesh,
        &[0],
        ShellSupport::Clamped,
        ReferenceSurface::Offset(0.003),
        0.0,
        ShellViscosityBudget {
            max_elements: 2,
            max_full_dofs: 24,
            max_dofs: 24,
            max_modes: 18,
            max_work: 20_000_000,
        },
        cx,
    )
}
fn observation_budget() -> TensorShellMaterialCausalBudget {
    TensorShellMaterialCausalBudget {
        field: TensorShellMaterialBudget {
            maximum_points: 36,
            maximum_payload_bytes: 1 << 27,
            maximum_work: 10_000_000_000,
        },
        maximum_frequencies: 2,
        control: TensorStressControlBudget::default(),
        temporal: TensorPlaneStressHarmonicBudget {
            maximum_payload_bytes: 4_000_000,
            maximum_scalar_work: 5_000_000,
        },
    }
}
fn controls() -> LayeredThermalSolveControls {
    LayeredThermalSolveControls {
        cell_tolerance_j: 2e-10,
        total_tolerance_j: 1e-9,
        ..LayeredThermalSolveControls::default()
    }
}
fn motion(owner: &CausalLayeredShell<'_, '_>) -> [Vec<C64>; 2] {
    AMPLITUDES.map(|amplitude| {
        let shell = owner.shell();
        let mut displacement = vec![C64::ZERO; shell.model().free];
        for (node, point) in shell.mesh().nodes.iter().enumerate() {
            if let Some(free) = shell.model().dof_map[6 * node] {
                displacement[free] = C64::new(0.0, amplitude * point[0]);
            }
        }
        displacement
    })
}

// Independent real scalar source solution q'+q/tau=epsilon', q(0)=0.
// nu=0 and isotropic original branches make this affine membrane's tensor norm
// independent of either facet frame and the declared material-layer rotation.
fn spring_strain(t: f64, tau: f64) -> f64 {
    FREQUENCIES
        .iter()
        .zip(AMPLITUDES)
        .map(|(&omega, amplitude)| {
            let x = omega * tau;
            amplitude * (x * ((omega * t).cos() - (-t / tau).exp()) + x * x * (omega * t).sin())
                / (1.0 + x * x)
        })
        .sum()
}
fn initial_layer_heat(layer: usize) -> f64 {
    let steps = 4096;
    let h = DT / steps as f64;
    let mut integral = 0.0;
    for i in 0..=steps {
        let weight = if i == 0 || i == steps {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * spring_strain(i as f64 * h, TAUS[layer]).powi(2);
    }
    common::AREA * common::H / 2.0 * FRACTIONS[layer] * common::E / TAUS[layer] * integral * h / 3.0
}
fn cell_heat(layer_heat: [f64; 2]) -> [f64; 3] {
    [layer_heat[0], 0.5 * layer_heat[1], 0.5 * layer_heat[1]]
}
// Independent finite-volume backward Euler equations from the authored cell
// thicknesses/masses; no native residual, Jacobian or thermal receipt enters.
fn thermal_reference(old: [f64; 3], heat: [f64; 3]) -> [f64; 3] {
    let capacity = [600.0, 300.0, 300.0];
    let mut a = [
        [capacity[0] / DT + 800.0, -800.0, 0.0],
        [-800.0, capacity[1] / DT + 2000.0, -1200.0],
        [0.0, -1200.0, capacity[2] / DT + 1200.0],
    ];
    let mut rhs =
        core::array::from_fn::<_, 3, _>(|i| (capacity[i] * (old[i] - 300.0) + heat[i]) / DT);
    for k in 0..3 {
        for i in k + 1..3 {
            let ratio = a[i][k] / a[k][k];
            for j in k..3 {
                a[i][j] -= ratio * a[k][j];
            }
            rhs[i] -= ratio * rhs[k];
        }
    }
    let mut answer = [0.0; 3];
    for i in (0..3).rev() {
        answer[i] = (rhs[i] - (i + 1..3).map(|j| a[i][j] * answer[j]).sum::<f64>()) / a[i][i];
    }
    answer.map(|v| v + 300.0)
}

#[test]
fn original_layer_heat_reaches_subdivided_caloric_cells_and_retained_holds() {
    with_cx(|cx, _| {
        let cards = [source_card(0, false), source_card(1, false)];
        let body = common::prepare_body(&cards, &[1, 2], 1.0);
        let laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&laws, common::ids()).unwrap();
        let thermal = common::initial(&bound, &[300.0; 3]);
        let old_cells = thermal.cells().clone();
        let shell = owner(&bound, &cards, cx).unwrap();
        assert_eq!(shell.cell_source_layers(), &[0, 1, 1]);
        let observer = shell
            .shell()
            .prepare_viscoelastic_material_causal_observer(
                observation_budget(),
                TensorKelvinStressControlBudget::default(),
                cx,
            )
            .unwrap();
        let history = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let fields = motion(&shell);
        let refs = fields.each_ref().map(|v| v.as_slice());
        let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
        let first = shell
            .heat_candidate(
                &thermal,
                DT,
                &observer,
                &history,
                &FREQUENCIES,
                &refs,
                &vec![0.0; shell.shell().model().free],
                common::conduction(),
                controls(),
                &mut workspace,
                ENERGY_TOLERANCE,
                cx,
            )
            .unwrap();
        let expected_heat = cell_heat(core::array::from_fn(initial_layer_heat));
        let expected_temperature = thermal_reference([300.0; 3], expected_heat);
        assert_eq!(first.report().material.points.len(), 36);
        for (actual, expected) in first.report().cell_heat_j.iter().zip(expected_heat) {
            common::close(*actual, expected, 2e-11 * expected);
        }
        for (actual, expected) in bound
            .cell_temperatures(first.thermal(), || Ok(()))
            .unwrap()
            .into_iter()
            .zip(expected_temperature)
        {
            common::close(actual, expected, 5e-10);
        }
        let total_heat: f64 = expected_heat.iter().sum();
        common::close(
            first.report().deposited_heat_j,
            total_heat,
            2e-11 * total_heat,
        );
        common::close(
            first.thermal().energy_j() - thermal.energy_j(),
            total_heat,
            ENERGY_TOLERANCE,
        );
        common::close(first.report().other_external_heat_j, 0.0, ENERGY_TOLERANCE);
        common::close(first.report().heat_handoff_residual_j, 0.0, 1e-12);
        common::close(
            first.report().constitutive_thermal_energy_residual_j,
            0.0,
            ENERGY_TOLERANCE,
        );
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(history.time_s_dd(), Dd::ZERO);
        assert_eq!(history.dissipated_energy_j(), 0.0);

        let (accepted, caloric, _) = first.into_parts();
        let mut hold = vec![0.0; shell.shell().model().free];
        for (full, &free) in shell.shell().model().dof_map.iter().enumerate() {
            if let Some(free) = free {
                hold[free] = accepted.full_displacement()[full];
            }
        }
        let second = shell
            .heat_candidate(
                &caloric,
                2.0 * DT,
                &observer,
                &accepted,
                &[],
                &[],
                &hold,
                common::conduction(),
                controls(),
                &mut workspace,
                ENERGY_TOLERANCE,
                cx,
            )
            .unwrap();
        let held_heat = cell_heat(core::array::from_fn(|layer| {
            common::AREA * common::H / 2.0
                * FRACTIONS[layer]
                * common::E
                * 0.5
                * spring_strain(DT, TAUS[layer]).powi(2)
                * -(-2.0 * DT / TAUS[layer]).exp_m1()
        }));
        for (actual, expected) in second.report().cell_heat_j.iter().zip(held_heat) {
            common::close(*actual, expected, 3e-11 * expected);
        }
        for (actual, expected) in bound
            .cell_temperatures(second.thermal(), || Ok(()))
            .unwrap()
            .into_iter()
            .zip(thermal_reference(expected_temperature, held_heat))
        {
            common::close(actual, expected, 5e-10);
        }
        assert!(second.report().material.material_heat_j > 0.0);
        assert_eq!(second.report().material.stress_work_j, 0.0);
        assert_eq!(second.history().stress_work_j(), accepted.stress_work_j());
        common::close(
            second.history().dissipated_energy_j(),
            accepted.dissipated_energy_j() + held_heat.iter().sum::<f64>(),
            1e-10,
        );
        common::close(
            second.report().constitutive_thermal_energy_residual_j,
            0.0,
            ENERGY_TOLERANCE,
        );
        assert_eq!(accepted.time_s_dd(), Dd::from_f64(DT));
        assert_eq!(caloric.time_s(), DT);
        assert_eq!(second.history().time_s_dd(), Dd::from_f64(2.0 * DT));
    });
}

#[test]
fn source_inventory_and_late_heat_refusals_preserve_both_accepted_histories() {
    with_cx(|cx, gate| {
        let cards = [source_card(0, false), source_card(1, false)];
        let body = common::prepare_body(&cards, &[1, 2], 1.0);
        let laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&laws, common::ids()).unwrap();
        let thermal = common::initial(&bound, &[300.0; 3]);
        let cells = thermal.cells().clone();
        let links = thermal.links().clone();
        let shell = owner(&bound, &cards, cx).unwrap();
        let observer = shell
            .shell()
            .prepare_viscoelastic_material_causal_observer(
                observation_budget(),
                TensorKelvinStressControlBudget::default(),
                cx,
            )
            .unwrap();
        let history = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let fields = motion(&shell);
        let refs = fields.each_ref().map(|v| v.as_slice());
        let bias = vec![0.0; shell.shell().model().free];
        let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
        let expected = shell
            .heat_candidate(
                &thermal,
                DT,
                &observer,
                &history,
                &FREQUENCIES,
                &refs,
                &bias,
                common::conduction(),
                controls(),
                &mut workspace,
                ENERGY_TOLERANCE,
                cx,
            )
            .unwrap();
        // The material interval succeeds before the native heat solver refuses
        // its wrong-size workspace. Neither candidate may be partially accepted.
        let mut too_small = LayeredThermalWorkspace::new(2).unwrap();
        assert!(matches!(
            shell.heat_candidate(
                &thermal,
                DT,
                &observer,
                &history,
                &FREQUENCIES,
                &refs,
                &bias,
                common::conduction(),
                controls(),
                &mut too_small,
                ENERGY_TOLERANCE,
                cx,
            ),
            Err(LayeredShellHeatError::HeatSolve(_))
        ));
        let other = shell
            .shell()
            .prepare_viscoelastic_material_causal_observer(
                observation_budget(),
                TensorKelvinStressControlBudget::default(),
                cx,
            )
            .unwrap();
        assert!(
            shell
                .heat_candidate(
                    &thermal,
                    DT,
                    &other,
                    &history,
                    &FREQUENCIES,
                    &refs,
                    &bias,
                    common::conduction(),
                    controls(),
                    &mut workspace,
                    ENERGY_TOLERANCE,
                    cx,
                )
                .is_err()
        );
        let retry = shell
            .heat_candidate(
                &thermal,
                DT,
                &observer,
                &history,
                &FREQUENCIES,
                &refs,
                &bias,
                common::conduction(),
                controls(),
                &mut workspace,
                ENERGY_TOLERANCE,
                cx,
            )
            .unwrap();
        assert_eq!(retry.report().cell_heat_j, expected.report().cell_heat_j);
        assert_eq!(retry.thermal().cells(), expected.thermal().cells());
        assert_eq!(
            retry.history().full_displacement(),
            expected.history().full_displacement()
        );

        let wrong_mass = common::prepare_body(&cards, &[1, 2], 1.1);
        let wrong_laws = wrong_mass.thermal_laws(false).unwrap();
        let wrong_bound = wrong_mass.bind(&wrong_laws, common::ids()).unwrap();
        assert!(owner(&wrong_bound, &cards, cx).is_err());
        let varying = [source_card(0, true), source_card(1, false)];
        let varying_body = common::prepare_body(&varying, &[1, 2], 1.0);
        let varying_laws = varying_body.thermal_laws(false).unwrap();
        let varying_bound = varying_body.bind(&varying_laws, common::ids()).unwrap();
        assert!(owner(&varying_bound, &varying, cx).is_err());
        gate.request();
        assert!(matches!(
            shell.heat_candidate(
                &thermal,
                DT,
                &observer,
                &history,
                &FREQUENCIES,
                &refs,
                &bias,
                common::conduction(),
                controls(),
                &mut workspace,
                ENERGY_TOLERANCE,
                cx,
            ),
            Err(LayeredShellHeatError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.links(), &links);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(history.time_s_dd(), Dd::ZERO);
        assert_eq!(history.full_displacement(), &[0.0; 24]);
        assert_eq!(history.stored_energy_j(), 0.0);
        assert_eq!(history.stress_work_j(), 0.0);
        assert_eq!(history.dissipated_energy_j(), 0.0);
    });
}
