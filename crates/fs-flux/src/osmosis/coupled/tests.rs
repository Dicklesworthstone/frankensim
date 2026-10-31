//! Independent coupled water/species dynamics and actual mechanical-work gates.
use super::*;
use crate::osmosis::{OsmoticPressureVolume, coupled::CoupledSolutionControls};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::{
    osmosis::{
        AqueousSolution, HydraulicPermeance, MembraneTransport, ReflectionCoefficient,
        SolutePermeance,
        diffusion::{
            SolutionDiffusionError, SolutionDiffusionSolve, SolutionMaxwellStefan,
            source::{
                BINARY_SOLUTION_DIFFUSIVITY_DIMS as DIMS, BinarySolutionDiffusionSource,
                MaterialSolutionMaxwellStefan, SolutionDiffusionSourceAxes,
            },
        },
        nonideal::VirialSolutionModel,
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    transport::maxwell_stefan::MaxwellStefanMolarSolveV1,
};
use std::sync::Arc;

const T: f64 = 298.15;
const VBAR: f64 = 18e-6;
const B: f64 = 0.005;
const VL: f64 = 1e-6;
const VR: f64 = 2e-6;
const K: f64 = 2e11;
const LP: f64 = 5e-12;
const AREA: f64 = 0.01;
const LENGTH: f64 = 0.001;
fn model(reference: f64) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, VBAR, reference, vec![200.0], vec![B]).unwrap())
}
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn cell(
    model: &Arc<VirialSolutionModel>,
    volume: f64,
    c: f64,
    pressure: f64,
    k: f64,
) -> OsmoticCompartment {
    OsmoticCompartment::new(
        AqueousSolution::nonideal(model.clone(), volume, vec![volume * c]).unwrap(),
        OsmoticPressureVolume::new(volume, pressure, k).unwrap(),
        volume * 0.1,
    )
    .unwrap()
}
fn membrane(left: usize, right: usize) -> OsmoticMembrane {
    OsmoticMembrane::new(
        left,
        right,
        AREA,
        MembraneTransport::new(
            HydraulicPermeance::new(LP).unwrap(),
            vec![(
                ReflectionCoefficient::new(1.0).unwrap(),
                SolutePermeance::new(0.0).unwrap(),
            )],
        )
        .unwrap(),
    )
    .unwrap()
}
fn physical() -> CoupledSolutionControls {
    CoupledSolutionControls {
        diffusion: SolutionDiffusionControls {
            constitutive: SolutionDiffusionSolve {
                friction: MaxwellStefanMolarSolveV1::default(),
                entropy_tolerance_w_m3_k: 1e-10,
            },
            maximum_species_fraction: 0.2,
            amount_balance_tolerance_mol: 1e-16,
            free_energy_tolerance_j: 1e-11,
        },
        water_balance_tolerance_m3: 1e-18,
    }
}
fn temporal(relative: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: 20.0,
        relative_tolerance: relative,
        absolute_volume_tolerance_m3: 1e-18,
        absolute_amount_tolerance_mol: 1e-14,
        maximum_steps: 10000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}
fn constant(model: &Arc<VirialSolutionModel>) -> SolutionMaxwellStefan {
    SolutionMaxwellStefan::new(
        model.clone(),
        (id("water"), MolarMass::new(0.018)),
        vec![(id("A"), MolarMass::new(0.06))],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )
    .unwrap()
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn sourced(model: &Arc<VirialSolutionModel>, range: [f64; 2]) -> MaterialSolutionMaxwellStefan {
    let key = PropertyKey::new("binary", DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "c_A".into(),
                abscissa_dims: typed(QuantityKind::AmountConcentration).dims(),
                knots: range
                    .into_iter()
                    .map(|c| (c, 1e-9 * (1.0 + c / 100.0)))
                    .collect(),
                dims: DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), 290.0, 310.0)
                .with_quantity(
                    "c_A",
                    typed(QuantityKind::AmountConcentration),
                    range[0],
                    range[1],
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured coupled solution mobility".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let material = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "assigned aqueous pair".into(),
                phase: "liquid".into(),
                process: "analytic fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    MaterialSolutionMaxwellStefan::new(
        model.clone(),
        (id("water"), MolarMass::new(0.018)),
        vec![(id("A"), MolarMass::new(0.06))],
        vec![BinarySolutionDiffusionSource {
            species: [id("A"), id("water")],
            material,
            property: key,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        SolutionDiffusionSourceAxes {
            temperature: "T".into(),
            concentrations: vec![(id("A"), "c_A".into()), (id("water"), "c_water".into())],
        },
        cell(model, VL, 100.0, 0.0, 0.0).solution(),
        || false,
    )
    .unwrap()
}
fn network(reference: f64, pressure: f64, dynamic: bool, reverse: bool) -> OsmoticNetwork {
    let model = model(reference);
    let (left, right) = if reverse { (1, 0) } else { (0, 1) };
    let link = if dynamic {
        SolutionDiffusionLink::from_material_sources(
            left,
            right,
            AREA,
            LENGTH,
            sourced(&model, [1.0, 200.0]),
        )
        .unwrap()
    } else {
        SolutionDiffusionLink::new(left, right, AREA, LENGTH, constant(&model)).unwrap()
    };
    OsmoticNetwork::with_coupled_solution_transport(
        vec!["A".into()],
        vec![
            cell(&model, VL, 100.0, pressure, K),
            cell(&model, VR, 20.0, pressure, K),
        ],
        vec![membrane(left, right)],
        vec![link],
        physical(),
        1 << 22,
    )
    .unwrap()
}
fn bits(system: &OsmoticNetwork) -> Vec<u64> {
    let mut bits = vec![system.time_s().to_bits()];
    for cell in system.compartments() {
        bits.push(cell.solution().solvent_volume_m3().to_bits());
        bits.extend(cell.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    bits
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_coupled_nonideal_water_and_material_mobility_refine_to_independent_ode() {
    // State is actual [left solvent volume, left solute amount]. Both right
    // inventories follow independent closed-system conservation.
    let rate = |state: [f64; 2]| {
        let (vl, nl) = (state[0], state[1]);
        let (vr, nr) = (VL + VR - vl, 100.0 * VL + 20.0 * VR - nl);
        let (cl, cr) = (nl / vl, nr / vr);
        let osmotic = R * T * ((cl - cr) + 0.5 * B * (cl * cl - cr * cr));
        let pressure = K * (vl - VL) - K * (vr - VR);
        let mean = 0.5 * (cl + cr);
        let molecular = 1e-9
            * (1.0 + mean / 100.0)
            * mean
            * (1.0 + mean * VBAR)
            * ((cl / cr).ln() + B * (cl - cr))
            / LENGTH;
        [-AREA * LP * (pressure - osmotic), -AREA * molecular]
    };
    let mut exact = [VL, 100.0 * VL];
    let h = 20.0 / 10000.0;
    for _ in 0..10000 {
        let a = rate(exact);
        let b = rate([exact[0] + 0.5 * h * a[0], exact[1] + 0.5 * h * a[1]]);
        let c = rate([exact[0] + 0.5 * h * b[0], exact[1] + 0.5 * h * b[1]]);
        let d = rate([exact[0] + h * c[0], exact[1] + h * c[1]]);
        for i in 0..2 {
            exact[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for tolerance in [1e-3, 6.25e-5] {
        let mut system = network(1.0, 1e5, true, false);
        let report = system.advance(20.0, temporal(tolerance), || false).unwrap();
        let left = system.compartments()[0].solution();
        assert!(left.solvent_volume_m3() > VL && left.solute_moles()[0] < 100.0 * VL);
        let error = ((left.solvent_volume_m3() - exact[0]) / VL)
            .abs()
            .max(((left.solute_moles()[0] - exact[1]) / (100.0 * VL)).abs());
        errors.push(error);
        assert!(report.accepted_steps > 1 && report.rejected_steps > 0);
        assert!(report.water_balance_residual_m3.abs() <= physical().water_balance_tolerance_m3);
        assert!(
            report.solute_balance_residual_mol[0].abs()
                <= physical().diffusion.amount_balance_tolerance_mol
        );
        assert!(report.diffusion_free_energy_change_j.is_none());
        assert!(
            report
                .coupled_potential_change
                .unwrap()
                .total_potential_change_j
                < 0.0
        );
    }
    assert!(
        errors[0] < 0.01 && errors[1] < 0.4 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn g1_reverse_osmosis_can_raise_mixing_energy_while_total_potential_decreases() {
    let model = model(1.0);
    let cells = vec![
        cell(&model, VL, 100.0, 1e6, K),
        cell(&model, VL, 100.0, 1e5, K),
    ];
    let before = cells.clone();
    let mut system = OsmoticNetwork::with_coupled_solution_transport(
        vec!["A".into()],
        cells,
        vec![membrane(0, 1)],
        vec![],
        physical(),
        1 << 20,
    )
    .unwrap();
    let report = system.advance(1.0, temporal(1e-5), || false).unwrap();
    let account = report.coupled_potential_change.unwrap();
    assert!(account.mixing_free_energy_change_j > 0.0);
    assert!(account.mechanical_potential_change_j < 0.0);
    assert!(account.total_potential_change_j < 0.0);
    let mut mixing = 0.0;
    let mut mechanical = 0.0;
    for (i, (old, new)) in before.iter().zip(system.compartments()).enumerate() {
        let (v0, v1) = (
            old.solution().solvent_volume_m3(),
            new.solution().solvent_volume_m3(),
        );
        let n = old.solution().solute_moles()[0];
        mixing += v1 * model.mixing_free_energy_density_j_m3(&[n / v1]).unwrap()
            - v0 * model.mixing_free_energy_density_j_m3(&[n / v0]).unwrap();
        let p = if i == 0 { 1e6 } else { 1e5 };
        mechanical += p * (v1 - v0) + 0.5 * K * (v1 - v0).powi(2);
        assert_eq!(old.solution().solute_moles(), new.solution().solute_moles());
    }
    near(account.mixing_free_energy_change_j, mixing, 1e-11);
    near(account.mechanical_potential_change_j, mechanical, 1e-11);
    near(account.total_potential_change_j, mixing + mechanical, 1e-11);
}

#[test]
fn g3_common_pressure_concentration_references_and_edge_orientation_cancel() {
    let mut original = network(1.0, 0.0, false, false);
    let mut shifted = network(1e200, 1e200, false, true);
    let a = original.advance(5.0, temporal(1e-4), || false).unwrap();
    let b = shifted.advance(5.0, temporal(1e-4), || false).unwrap();
    assert_eq!(bits(&original), bits(&shifted));
    assert_eq!(a.coupled_potential_change, b.coupled_potential_change);
    assert_eq!(a.accepted_steps, b.accepted_steps);
}

#[test]
fn g4_water_dilution_refreshes_pinned_sources_and_late_failure_is_atomic() {
    let model = model(1.0);
    let link = SolutionDiffusionLink::from_material_sources(
        0,
        1,
        1e-6,
        LENGTH,
        sourced(&model, [1.0, 101.0]),
    )
    .unwrap();
    let mut system = OsmoticNetwork::with_coupled_solution_transport(
        vec!["A".into()],
        vec![
            cell(&model, VL, 100.0, 1e6, K),
            cell(&model, VR, 20.0, 1e5, K),
        ],
        vec![membrane(0, 1)],
        vec![link],
        physical(),
        1 << 22,
    )
    .unwrap();
    let controls = OsmoticStepControls {
        maximum_step_s: 0.1,
        ..temporal(1e-4)
    };
    system.advance(0.1, controls, || false).unwrap();
    let checkpoint = bits(&system);
    assert!(matches!(
        system.advance(10.0, controls, || false),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
    assert_eq!(bits(&system), checkpoint);
    let mut polls = 0;
    assert!(matches!(
        system.advance(0.1, controls, || {
            polls += 1;
            polls > 5
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&system), checkpoint);
}

#[test]
fn g4_coupled_transport_refuses_leaky_membranes_duplicate_faces_and_budget_exhaustion() {
    let model = model(1.0);
    let cells = vec![
        cell(&model, VL, 100.0, 0.0, K),
        cell(&model, VR, 20.0, 0.0, K),
    ];
    let leaky = OsmoticMembrane::new(
        0,
        1,
        AREA,
        MembraneTransport::new(
            HydraulicPermeance::new(LP).unwrap(),
            vec![(
                ReflectionCoefficient::new(0.5).unwrap(),
                SolutePermeance::new(1e-9).unwrap(),
            )],
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        OsmoticNetwork::with_coupled_solution_transport(
            vec!["A".into()],
            cells.clone(),
            vec![leaky],
            vec![],
            physical(),
            1 << 20
        )
        .is_err()
    );
    assert!(
        OsmoticNetwork::with_coupled_solution_transport(
            vec!["A".into()],
            cells.clone(),
            vec![membrane(0, 1), membrane(1, 0)],
            vec![],
            physical(),
            1 << 20
        )
        .is_err()
    );
    assert!(matches!(
        OsmoticNetwork::with_coupled_solution_transport(
            vec!["A".into()],
            cells,
            vec![membrane(0, 1)],
            vec![],
            physical(),
            1
        ),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let mut system = network(1.0, 0.0, false, false);
    let before = bits(&system);
    assert!(matches!(
        system.advance(
            1.0,
            OsmoticStepControls {
                maximum_steps: 1,
                maximum_step_s: 0.01,
                ..temporal(1e-3)
            },
            || false
        ),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&system), before);
}
