//! G1/G3/G4: finite-slope surface traction uses the original material/VOF clock.
//! Synthetic constitutive data verify equations and integration, not experiments.

use super::*;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::free_surface::enthalpy::thermocapillary::HeightGraphThermocapillarity;
use fs_material::liquid::surface::{
    ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

#[path = "height_graph/local.rs"]
mod local;

const DX: f64 = 0.001;
const DT: f64 = 0.001;

fn source_card(slope: f64) -> MaterialCard {
    let base = card(0.1, 1000.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == SURFACE_TENSION_PROPERTY {
            if slope == 0.0 {
                claim.value = PropertyValue::Scalar {
                    value: 0.001,
                    dims: SURFACE_TENSION_DIMS,
                };
                claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
            } else {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, 0.001 - 100.0 * slope),
                        (500.0, 0.001 + 100.0 * slope),
                    ],
                    dims: SURFACE_TENSION_DIMS,
                };
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn source(card: &MaterialCard) -> ResolvedLiquidSurfaceState {
    let mut requirements = [
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0])),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect::<Vec<_>>();
    requirements.push(
        ScalarPropertyRequirement::try_with_quantity(
            "thermal_conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let state = resolve_material_state_point(
        card,
        &point(400.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedLiquidSurfaceState::try_from_material(
        state,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    )
    .unwrap()
}

fn graph() -> HeightGraphThermocapillarity {
    HeightGraphThermocapillarity {
        normal_axis: 2,
        outward_sign: 1,
        maximum_slope: 0.75,
        minimum_surface_fill: 0.001,
    }
}

fn graph_cells(sloped: bool) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let mut fill = 0.0;
                for a in 0..8 {
                    for b in 0..8 {
                        let u = x as f64 + (a as f64 + 0.5) / 8.0 - 4.0;
                        let v = y as f64 + (b as f64 + 0.5) / 8.0 - 4.0;
                        let h = if sloped {
                            3.6 + 0.15 * u - 0.1 * v
                        } else {
                            3.5
                        };
                        fill += (h - z as f64).clamp(0.0, 1.0) / 64.0;
                    }
                }
                if fill > 0.0 {
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: fill,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    cells
}

fn with_graph<R>(
    slope: f64,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = source_card(slope);
    let source = source(&card);
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let caloric = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_height_graph_thermocapillarity(
                domain(),
                &graph_cells(true),
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                graph(),
                None,
                None,
                |p| 200000.0 + 2000.0 * (2.0 * p[0] as f64 - p[1] as f64),
                cx,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

#[test]
fn g1_height_graph_sourced_traction_has_physical_units_and_advects_heat_conservatively() {
    for slope in [1e-6, -1e-6] {
        with_graph(slope, |thermal, gate| {
            assert_eq!(
                thermal.height_graph_thermocapillary_surface(),
                Some(graph())
            );
            assert!(thermal.thermocapillary_surface().is_none());
            let mut totals = std::collections::BTreeMap::<[usize; 2], [f64; 3]>::new();
            for cell in graph_cells(true) {
                let [x, y, z] = cell.coordinates;
                let force = thermal.flow().lattice().force_density_at(x, y, z).unwrap();
                let total = totals.entry([x, y]).or_default();
                for d in 0..3 {
                    total[d] += cell.fill_fraction * force[d];
                }
            }
            // Independently project the physical gradient in N/m^2 onto the
            // exact plane, integrate the patch, then apply the SI/lattice scale.
            let normal = [-0.15, 0.1, 1.0];
            let square = normal.iter().map(|n| n * n).sum::<f64>();
            let ambient = [2.0 * slope / DX, -slope / DX, 0.0];
            let dot = ambient.iter().zip(normal).map(|(g, n)| g * n).sum::<f64>();
            for total in totals.values() {
                for d in 0..3 {
                    let physical =
                        (ambient[d] - normal[d] * dot / square) * DX * DX * square.sqrt();
                    let expected = physical * DT * DT / (1000.0 * DX.powi(4));
                    near(total[d], expected, 2e-10);
                }
            }
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            let energy = thermal.total_enthalpy_j().unwrap();
            let report = thermal.step(gate).unwrap();
            near(report.mass_after_kg, mass, 2e-12);
            near(report.enthalpy_after_j, energy, 2e-11);
            assert_eq!(report.external_heat_j, 0.0);
            assert_eq!(report.viscous_heat_j, 0.0);
            assert!(report.balance_defect_j.abs() < 1e-9);
        });
    }
}

#[test]
fn g3_height_graph_constant_sigma_preserves_original_thermal_trajectory_bits() {
    let card = source_card(0.0);
    let source = source(&card);
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let caloric = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, DX, DT, 0.1).unwrap();
        let initial = graph_cells(false);
        let mut original = scaling
            .free_surface3_from_cells(domain(), &initial, cx)
            .unwrap()
            .with_enthalpy(
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let mut curved = scaling
            .free_surface3_from_cells_with_height_graph_thermocapillarity(
                domain(),
                &initial,
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                graph(),
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            )
            .unwrap();
        for _ in 0..3 {
            assert_eq!(
                original.flow().lattice().grid().state_bits(),
                curved.flow().lattice().grid().state_bits()
            );
            assert_eq!(
                original.total_enthalpy_j().unwrap().to_bits(),
                curved.total_enthalpy_j().unwrap().to_bits()
            );
            original.step(&gate).unwrap();
            curved.step(&gate).unwrap();
        }
    });
}

#[test]
fn g4_height_graph_source_refusal_and_cancellation_publish_no_partial_flow_or_heat() {
    with_graph(1e-6, |thermal, _| {
        let before = thermal.flow().lattice().grid().state_bits();
        let mass = thermal.flow().ledger_mass_kg().unwrap().to_bits();
        let energy = thermal.total_enthalpy_j().unwrap().to_bits();
        let gate = CancelGate::new_clock_free();
        assert!(
            thermal
                .step_with_heat(
                    &[FreeSurfaceCellHeat {
                        coordinates: [3, 3, 3],
                        energy_j: 1000.0,
                    }],
                    &gate
                )
                .is_err()
        );
        assert_eq!(before, thermal.flow().lattice().grid().state_bits());
        assert_eq!(mass, thermal.flow().ledger_mass_kg().unwrap().to_bits());
        assert_eq!(energy, thermal.total_enthalpy_j().unwrap().to_bits());
        gate.request();
        assert!(thermal.step(&gate).is_err());
        assert_eq!(before, thermal.flow().lattice().grid().state_bits());
        assert_eq!(energy, thermal.total_enthalpy_j().unwrap().to_bits());
    });
}
