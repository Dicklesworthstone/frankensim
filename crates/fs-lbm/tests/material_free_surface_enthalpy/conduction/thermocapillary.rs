//! Numerical planar thermocapillary flow; authored coefficients, not validation.

use super::*;
use fs_lbm::material::free_surface::enthalpy::thermocapillary::PlanarThermocapillarity;

const DX: f64 = 0.001;
const DT: f64 = 0.001;
const SURFACE: [usize; 3] = [3, 3, 3];

fn thermal_card(slope: f64, upper: f64) -> MaterialCard {
    let base = source_card(0.01, false, 500.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == SURFACE_TENSION_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.01 - 100.0 * slope), (500.0, 0.01 + 100.0 * slope)],
                dims: SURFACE_TENSION_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claim.validity = ValidityDomain::unconstrained()
                .with("T", 300.0, upper)
                .with("pressure", 100000.0, 100000.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn planar_cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    (0..4)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..8).map(move |x| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if z == 3 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect()
}

fn plane() -> PlanarThermocapillarity {
    PlanarThermocapillarity {
        normal_axis: 2,
        outward_sign: 1,
        surface_position_m: 3.5 * DX,
        maximum_departure_m: 0.01 * DX,
        minimum_surface_fill: 0.25,
    }
}

fn with_planar<R>(
    card: &MaterialCard,
    gradient_k: [f64; 3],
    surface: PlanarThermocapillarity,
    controls: FreeSurfaceEnthalpyControls,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_thermocapillarity(
                domain(),
                &planar_cells(),
                &law,
                &chart,
                controls,
                surface,
                None,
                None,
                |p| 200000.0 + 2000.0 * (0..3).map(|d| gradient_k[d] * p[d] as f64).sum::<f64>(),
                cx,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

#[test]
fn g1_g3_source_gradient_drives_native_flow_with_physical_scale_and_replay() {
    let run = |slope: f64| {
        let card = thermal_card(slope, 500.0);
        with_planar(
            &card,
            [3.0, 0.0, 0.0],
            plane(),
            FreeSurfaceEnthalpyControls::default(),
            |thermal, gate| {
                let force = thermal.flow().lattice().force_density_at(3, 3, 3).unwrap();
                // dσ/dx = slope*3K/dx; divide by the occupied half-cell
                // thickness, rho0, then convert acceleration by dt²/dx.
                let expected = (slope * 3.0 / DX) / (0.5 * DX) / 1000.0 * DT * DT / DX;
                near(force[0], expected, 2e-12);
                assert_eq!(force[1], 0.0);
                assert_eq!(force[2], 0.0);
                assert!(
                    thermal
                        .flow()
                        .sample_at_si(3, 3, 3)
                        .unwrap()
                        .unwrap()
                        .velocity_m_s
                        .iter()
                        .all(|u| u.abs() < 1e-15),
                    "force binding must preserve the initially stationary liquid"
                );
                let original = thermal
                    .material_at(SURFACE)
                    .unwrap()
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                let energy = thermal.total_enthalpy_j().unwrap();
                let mass = thermal.flow().ledger_mass_kg().unwrap();
                for _ in 0..3 {
                    let report = thermal.step(gate).unwrap();
                    near(report.mass_after_kg, mass, 3e-12);
                    near(report.enthalpy_after_j, energy, 3e-11);
                    assert_eq!(report.external_heat_j, 0.0);
                    assert_eq!(report.viscous_heat_j, 0.0);
                    let state = thermal.material_at(SURFACE).unwrap().material();
                    let sigma = state.property(SURFACE_TENSION_PROPERTY).unwrap();
                    assert_eq!(sigma.answer().receipt.selected, original);
                    near(
                        thermal
                            .flow()
                            .lattice()
                            .surface_tension_at(3, 3, 3)
                            .unwrap(),
                        sigma.value_si() * DT * DT / (1000.0 * DX * DX * DX),
                        2e-14,
                    );
                }
                let velocity = thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s;
                assert!(
                    velocity[0] * slope > 0.0,
                    "surface flows toward higher tension"
                );
                (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.total_enthalpy_j().unwrap().to_bits(),
                    velocity[0],
                )
            },
        )
    };
    let positive = run(0.00001);
    assert_eq!(positive, run(0.00001));
    let negative = run(-0.00001);
    assert_ne!(positive.0, negative.0);
    assert!(positive.2 > 0.0 && negative.2 < 0.0);

    let normal = thermal_card(0.00001, 500.0);
    with_planar(
        &normal,
        [0.0, 0.0, 3.0],
        plane(),
        FreeSurfaceEnthalpyControls::default(),
        |thermal, _| {
            assert_eq!(
                thermal.flow().lattice().force_density_at(3, 3, 3),
                Some([0.0; 3])
            );
        },
    );
}

#[test]
fn g1_constant_planar_limit_preserves_original_thermal_population_bits() {
    let card = source_card(0.01, false, 500.0);
    let source = liquid_source(&card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, DX, DT, 0.1).unwrap();
        let mut original = scaling
            .free_surface3_from_cells(domain(), &planar_cells(), cx)
            .unwrap()
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let mut planar = scaling
            .free_surface3_from_cells_with_thermocapillarity(
                domain(),
                &planar_cells(),
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                plane(),
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            )
            .unwrap();
        original.step(&gate).unwrap();
        planar.step(&gate).unwrap();
        assert_eq!(
            original.flow().lattice().grid().state_bits(),
            planar.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            original.total_enthalpy_j().unwrap(),
            planar.total_enthalpy_j().unwrap()
        );
    });
}

#[test]
fn g4_late_source_geometry_and_work_refusals_keep_all_accepted_histories() {
    let card = thermal_card(-0.00001, 401.0);
    with_planar(
        &card,
        [0.0; 3],
        plane(),
        FreeSurfaceEnthalpyControls::default(),
        |thermal, gate| {
            let bits = thermal.flow().lattice().grid().state_bits();
            let energy = thermal.total_enthalpy_j().unwrap();
            let source = thermal.material_at(SURFACE).unwrap().clone();
            // One half-full cell has 0.5mg and 0.001J/K heat capacity.
            // The candidate would warm it 2K beyond the pinned sigma domain.
            assert!(
                thermal
                    .step_with_heat(
                        &[FreeSurfaceCellHeat {
                            coordinates: SURFACE,
                            energy_j: 0.002,
                        }],
                        gate
                    )
                    .is_err()
            );
            assert_eq!(bits, thermal.flow().lattice().grid().state_bits());
            assert_eq!(energy, thermal.total_enthalpy_j().unwrap());
            assert_eq!(&source, thermal.material_at(SURFACE).unwrap());
            assert_eq!(
                thermal.flow().lattice().force_density_at(3, 3, 3),
                Some([0.0; 3])
            );
            let report = thermal
                .step_with_heat(
                    &[FreeSurfaceCellHeat {
                        coordinates: SURFACE,
                        energy_j: 0.0005,
                    }],
                    gate,
                )
                .unwrap();
            near(report.enthalpy_after_j - energy, 0.0005, 2e-8);
            let accepted = thermal.flow().lattice().grid().state_bits();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(thermal.step(&cancelled).is_err());
            assert_eq!(accepted, thermal.flow().lattice().grid().state_bits());
        },
    );

    let card = thermal_card(0.00001, 500.0);
    for (surface, controls) in [
        (
            PlanarThermocapillarity {
                maximum_departure_m: DX * 1e-12,
                ..plane()
            },
            FreeSurfaceEnthalpyControls::default(),
        ),
        (
            plane(),
            FreeSurfaceEnthalpyControls {
                max_mass_routes: 1,
                ..FreeSurfaceEnthalpyControls::default()
            },
        ),
    ] {
        with_planar(
            &card,
            [3.0, 0.0, 0.0],
            surface,
            controls,
            |thermal, gate| {
                let bits = thermal.flow().lattice().grid().state_bits();
                let energy = thermal.total_enthalpy_j().unwrap();
                let force = thermal.flow().lattice().force_density_at(3, 3, 3);
                let sigma = thermal.flow().lattice().surface_tension_at(3, 3, 3);
                assert!(thermal.step(gate).is_err());
                assert_eq!(bits, thermal.flow().lattice().grid().state_bits());
                assert_eq!(energy, thermal.total_enthalpy_j().unwrap());
                assert_eq!(force, thermal.flow().lattice().force_density_at(3, 3, 3));
                assert_eq!(sigma, thermal.flow().lattice().surface_tension_at(3, 3, 3));
            },
        );
    }
}
