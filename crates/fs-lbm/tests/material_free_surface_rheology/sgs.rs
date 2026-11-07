//! Coupled source/SGS roots on actual native populations and changing wet cells.

use super::*;
use fs_lbm::d3q19::freesurface3::FreeSurfaceKinematicsWorkspace3;
use fs_lbm::d3q19::freesurface3::sgs::FreeSurfaceSgsWorkspace3;
use fs_lbm::material::free_surface::MaterialFreeSurfaceSgsConfig3;
use fs_lbm::sgs3::{SgsModel, SgsUse, eddy_viscosity};

mod motion;

fn controls(model: SgsModel) -> MaterialFreeSurfaceSgsConfig3 {
    MaterialFreeSurfaceSgsConfig3 {
        model,
        filter_width_m: 1.7 * DX,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn cache() -> FreeSurfaceKinematicsWorkspace3 {
    FreeSurfaceKinematicsWorkspace3::new(8, usize::MAX).unwrap()
}

fn bulk_card(slope: f64) -> MaterialCard {
    let base = card(slope, 1000.0, 0.001);
    let mut claims = base.claims().clone();
    let dims = Dims([-1, 1, -1, 0, 0, 0]);
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("bulk_eta", dims),
            value: PropertyValue::Curve {
                abscissa: "rate".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, 0.2), (1000.0, 4.2)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 500.0)
                .with("pressure", 100000.0, 100000.0)
                .with_quantity(
                    "rate",
                    QuantitySpec::dimensional(Frequency::DIMS),
                    0.0,
                    1000.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic independent bulk response".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::supersede(&base, claims, vec![]).unwrap()
}

fn bulk_law(source: &MaterialCard) -> IsothermalRheologyLaw<'_> {
    let original = bundle(source);
    let mut required: Vec<_> = original
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    required.push(
        ScalarPropertyRequirement::try_new(
            "bulk_eta",
            Dims([-1, 1, -1, 0, 0, 0]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let material = resolve_material_state_point(
        source,
        &original.query(),
        &required,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    IsothermalRheologyLaw::try_new(
        source,
        material,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
    )
    .unwrap()
}

fn primed<'s, 'm>(
    law: &'s IsothermalRheologyLaw<'m>,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> MaterialFreeSurfaceRheology3<'s, 'm> {
    let cells: Vec<_> = initial(true)
        .into_iter()
        .map(|mut cell| {
            let [x, y, z] = cell.coordinates.map(|v| v as f64 - 3.5);
            cell.velocity_m_s = [0.003 * y + 0.0002 * z * z, 0.002 * z, 0.001 * x];
            cell
        })
        .collect();
    let mut flow = MaterialFreeSurfaceRheology3::from_cells_with_pressure_reference(
        law,
        DX,
        DT,
        0.1,
        domain([0.03, -0.02, 0.01]),
        &cells,
        solve(),
        100000.0,
        cx,
        gate,
    )
    .unwrap();
    flow.set_bulk_viscosity_source("bulk_eta", 1.2, gate)
        .unwrap();
    for _ in 0..2 {
        flow.step_with_gate(gate).unwrap();
    }
    flow
}

#[test]
fn nonlinear_source_and_all_sgs_closures_solve_the_same_total_stress_relation() {
    let gate = CancelGate::new_clock_free();
    let source = bulk_card(0.003);
    let law = bulk_law(&source);
    with_cx(&gate, |cx| {
        let input = primed(&law, cx, &gate);
        let raw = numerator(&input, [0.03, -0.02, 0.01], CELL);
        assert!(raw > 1e-7);
        for model in [
            SgsModel::Smagorinsky { cs: 0.6 },
            SgsModel::Wale { cw: 0.6 },
            SgsModel::Vreman { cv: 0.2 },
        ] {
            let mut kinematics = cache();
            let mut workspace = FreeSurfaceSgsWorkspace3::default();
            let (a, b) = if let SgsModel::Smagorinsky { cs } = model {
                (0.8, 0.009 + 3.0 * (cs * 1.7).powi(2) * DT)
            } else {
                let field = input.lattice().kinematics(&gate, &mut kinematics).unwrap();
                let gradient = field
                    .sample_at(CELL[0], CELL[1], CELL[2])
                    .unwrap()
                    .unwrap()
                    .kinematics
                    .gradient();
                let eddy = eddy_viscosity(model, 1.7, gradient).unwrap();
                assert!(eddy > 0.0);
                (0.8 + 3.0 * eddy, 0.009)
            };
            let expected_rate = (2.0 * raw / DT) / (a + (a * a + 4.0 * b * raw / DT).sqrt());
            let mut flow = input
                .clone()
                .with_sgs_si(
                    controls(model),
                    &gate,
                    Some(&mut kinematics),
                    &mut workspace,
                )
                .unwrap();
            let sample = flow.rheology_at(CELL).unwrap();
            near(sample.shear_rate_per_s, expected_rate, 3e-9);
            near(sample.collision_tau, 0.8 + 0.009 * expected_rate, 3e-11);
            near(
                workspace.relaxation_at(CELL).unwrap(),
                a + b * expected_rate,
                3e-11,
            );
            near(
                sample.shear_rate_per_s * DT * workspace.relaxation_at(CELL).unwrap(),
                raw,
                3e-12,
            );
            let material = flow.material_at(CELL).unwrap();
            near(
                material.property("bulk_eta").unwrap().value_si(),
                0.2 + 0.004 * expected_rate,
                3e-11,
            );
            assert_eq!(
                material
                    .property("bulk_eta")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                law.material()
                    .property("bulk_eta")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            let report = flow
                .step_with_gate(&gate, Some(&mut kinematics), &mut workspace)
                .unwrap();
            assert!(report.sgs_dissipation_w > 0.0);
            assert!(report.resolved_dissipation_w > 0.0);
            assert!(report.resolved_bulk_dissipation_w >= 0.0);
            // Accepted endpoint material and the exposed total tau belong to
            // the same fresh population epoch, including its original source.
            let endpoint = flow.rheology_at(CELL).unwrap();
            let material = flow.material_at(CELL).unwrap();
            near(
                material.property("dynamic_viscosity").unwrap().value_si(),
                0.1 + 0.003 * endpoint.shear_rate_per_s,
                3e-12,
            );
            assert!(workspace.relaxation_at(CELL).unwrap() >= endpoint.collision_tau);
        }
    });
}

#[test]
fn zero_sgs_keeps_original_source_trajectory_and_cancellation_preserves_acceptance() {
    let gate = CancelGate::new_clock_free();
    let source = bulk_card(0.003);
    let law = bulk_law(&source);
    with_cx(&gate, |cx| {
        let input = primed(&law, cx, &gate);
        for model in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
            SgsModel::Vreman { cv: 0.0 },
        ] {
            let mut ordinary = input.clone();
            let mut workspace = FreeSurfaceSgsWorkspace3::default();
            let mut zero = input
                .clone()
                .with_sgs_si(controls(model), &gate, None, &mut workspace)
                .unwrap();
            for _ in 0..2 {
                ordinary.step_with_gate(&gate).unwrap();
                let report = zero.step_with_gate(&gate, None, &mut workspace).unwrap();
                assert_eq!(report.sgs_dissipation_w, 0.0);
                assert_eq!(
                    zero.lattice().grid().state_bits(),
                    ordinary.lattice().grid().state_bits()
                );
                assert_eq!(zero.material_at(CELL), ordinary.material_at(CELL));
                assert_eq!(zero.rheology_at(CELL), ordinary.rheology_at(CELL));
            }
            let bits = zero.lattice().grid().state_bits();
            let material = zero.material_at(CELL).unwrap().clone();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(
                zero.step_with_gate(&cancelled, None, &mut workspace)
                    .is_err()
            );
            assert_eq!(zero.lattice().grid().state_bits(), bits);
            assert_eq!(zero.material_at(CELL).unwrap(), &material);
        }
    });
}

#[test]
fn advancing_surface_rebinds_every_new_wet_bulk_source_before_the_next_collision() {
    let gate = CancelGate::new_clock_free();
    let source = bulk_card(0.003);
    let law = bulk_law(&source);
    with_cx(&gate, |cx| {
        let cells: Vec<_> = initial(false)
            .into_iter()
            .map(|mut cell| {
                cell.velocity_m_s = [0.03, 0.0, 0.0];
                cell
            })
            .collect();
        let mut molecular = pressure_flow(&law, &cells, 100000.0, cx, &gate);
        molecular
            .set_bulk_viscosity_source("bulk_eta", 1.2, &gate)
            .unwrap();
        let all: Vec<_> = (0..8)
            .flat_map(|z| (0..8).flat_map(move |y| (0..8).map(move |x| [x, y, z])))
            .collect();
        let original: Vec<_> = all
            .iter()
            .copied()
            .filter(|c| molecular.lattice().sample_at(c[0], c[1], c[2]).is_some())
            .collect();
        let mass = molecular.ledger_mass_kg().unwrap();
        let mut workspace = FreeSurfaceSgsWorkspace3::default();
        let mut flow = molecular
            .with_sgs_si(
                controls(SgsModel::Smagorinsky { cs: 0.2 }),
                &gate,
                None,
                &mut workspace,
            )
            .unwrap();
        let mut newly_wet = false;
        for _ in 0..4 {
            flow.step_with_gas_pressure(&gate, 0.1, |_| 100000.0, None, &mut workspace)
                .unwrap();
            near(flow.ledger_mass_kg().unwrap(), mass, 1e-10 * mass);
            for &c in &all {
                let wet = flow.lattice().sample_at(c[0], c[1], c[2]).is_some();
                if !wet {
                    assert!(workspace.relaxation_at(c).is_none());
                    continue;
                }
                newly_wet |= !original.contains(&c);
                let sample = flow.rheology_at(c).unwrap();
                let material = flow.material_at(c).unwrap();
                near(
                    material.property("bulk_eta").unwrap().value_si(),
                    0.2 + 0.004 * sample.shear_rate_per_s,
                    3e-11,
                );
                assert!(workspace.relaxation_at(c).unwrap() >= sample.collision_tau);
            }
        }
        assert!(newly_wet, "the fixture must exercise actual VOF conversion");
        let mut native_scratch = FreeSurfaceSgsWorkspace3::default();
        flow.lattice()
            .prepare_sgs(
                &gate,
                fs_lbm::sgs3::SgsConfig {
                    model: SgsModel::Smagorinsky { cs: 0.0 },
                    filter_width_lu: 1.7,
                    maximum_eddy_to_molecular: 100.0,
                    usage: SgsUse::Exploration,
                },
                &mut native_scratch,
            )
            .unwrap();
    });
}
