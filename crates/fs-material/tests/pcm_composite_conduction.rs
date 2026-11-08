//! G1/G3: original-source composite heat transport, latent storage and tangents.
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
    composite::PcmComposite,
};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::conduction::{
    composite::{PcmConductionMorphology as Morphology, PreparedPcmConductivity},
    phase::{PhaseConductivityError, PreparedPhaseConductivity},
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, Reservoir, StepConfig, ThermalLink, TimedEvent,
    conduction::{
        ConductionObserver, ConductiveLink, ConductiveReservoir,
        implicit::{
            ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
            ImplicitConductionWorkspace,
        },
    },
};

fn card(axis: &str, values: &[(f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("k", KD),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: if axis == "T" {
                    Temperature::DIMS
                } else {
                    Dims::NONE
                },
                knots: values.to_vec(),
                dims: KD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    600.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured composite conductivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: axis.into(),
            phase: "declared numerical fixture".into(),
            process: "none".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn constituents(
    matrix: ContentHash,
    pcm: ContentHash,
    offsets: [f64; 2],
) -> (EquilibriumEnthalpyPhaseCurve, EquilibriumEnthalpyPhaseCurve) {
    let knot = |h, t, f, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    };
    (
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            matrix,
            SolidLiquidPhase::Solid,
            vec![
                knot(offsets[0], 300.0, 0.0, 8000.0),
                knot(offsets[0] + 100.0, 400.0, 0.0, 8000.0),
                knot(offsets[0] + 200.0, 500.0, 0.0, 8000.0),
                knot(offsets[0] + 300.0, 600.0, 0.0, 8000.0),
            ],
        )
        .unwrap(),
        EquilibriumEnthalpyPhaseCurve::try_new(
            pcm,
            vec![
                knot(offsets[1], 300.0, 0.0, 1000.0),
                knot(offsets[1] + 200.0, 400.0, 0.0, 1000.0),
                knot(offsets[1] + 400.0, 400.0, 1.0, 800.0),
                knot(offsets[1] + 650.0, 500.0, 1.0, 800.0),
            ],
        )
        .unwrap(),
    )
}
fn prepare<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [
            curve.knots()[0].specific_enthalpy_j_kg,
            curve.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        || Ok(()),
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, allowance {tolerance:e}"
    );
}
fn matrix_card() -> MaterialCard {
    card(
        "T",
        &[(300.0, 20.0), (400.0, 30.0), (500.0, 50.0), (600.0, 150.0)],
    )
}
fn pcm_card() -> MaterialCard {
    card("f", &[(0.0, 2.0), (0.5, 6.0), (1.0, 3.0)])
}

#[test]
fn source_transport_uses_volume_fractions_and_complete_latent_tangents() {
    let (mc, pc) = (matrix_card(), pcm_card());
    let (matrix, pcm) = constituents(mc.content_hash(), pc.content_hash(), [0.0; 2]);
    let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let (mk, pk) = (prepare(&mc, &matrix), prepare(&pc, &pcm));
    for morphology in [Morphology::Parallel, Morphology::Series] {
        let law = PreparedPcmConductivity::new(&caloric, &mk, &pk, morphology, || Ok(())).unwrap();
        assert!(std::ptr::eq(law.caloric_law(), &caloric));
        for h in [62.5, 140.0, 160.0, 243.75] {
            let sample = law.state_at_specific_enthalpy(h).unwrap();
            let state = sample.caloric.state;
            let phi = (0.25 / state.pcm.bulk_density_kg_m3)
                / (0.75 / state.matrix.bulk_density_kg_m3 + 0.25 / state.pcm.bulk_density_kg_m3);
            assert!((phi - 0.25).abs() > 0.4);
            let km = if state.temperature_k <= 400.0 {
                20.0 + 0.1 * (state.temperature_k - 300.0)
            } else {
                30.0 + 0.2 * (state.temperature_k - 400.0)
            };
            let f = state.pcm.liquid_mass_fraction;
            let kp = if f <= 0.5 {
                2.0 + 8.0 * f
            } else {
                9.0 - 6.0 * f
            };
            let expected = match morphology {
                Morphology::Parallel => (1.0 - phi) * km + phi * kp,
                Morphology::Series => 1.0 / ((1.0 - phi) / km + phi / kp),
            };
            near(sample.conductivity_w_m_k, expected, 2e-13);
            let delta = 1e-4;
            let plus = law.state_at_specific_enthalpy(h + delta).unwrap();
            let minus = law.state_at_specific_enthalpy(h - delta).unwrap();
            near(
                sample.conductivity_slope_kg_w_m_k_j,
                (plus.conductivity_w_m_k - minus.conductivity_w_m_k) / (2.0 * delta),
                2e-9,
            );
            near(
                sample.caloric.pcm_volume_fraction_slope_kg_j,
                (plus.caloric.state.pcm_volume_fraction - minus.caloric.state.pcm_volume_fraction)
                    / (2.0 * delta),
                1e-11,
            );
            if (125.0..175.0).contains(&h) {
                assert_eq!(sample.caloric.temperature_slope_kg_k_j, 0.0);
                assert_eq!(sample.caloric.matrix_enthalpy_derivative, 0.0);
                assert_eq!(sample.caloric.pcm_enthalpy_derivative, 4.0);
                near(sample.caloric.liquid_mass_fraction_slope_kg_j, 0.005, 1e-16);
                assert_ne!(sample.conductivity_slope_kg_w_m_k_j, 0.0);
            }
            let observable = ImplicitConductionLaw::time_observables(&law, h).unwrap();
            assert_eq!(
                observable.phase_mass_fraction,
                Some(state.liquid_mass_fraction)
            );
        }
        // Tmax=500 K is an INTERIOR kink in the wider matrix source. Only its
        // admitted left tangent belongs to this composite's upper endpoint.
        let end = law.state_at_specific_enthalpy(312.5).unwrap();
        let left = law.state_at_specific_enthalpy(312.5 - 1e-4).unwrap();
        near(
            end.conductivity_slope_kg_w_m_k_j,
            (end.conductivity_w_m_k - left.conductivity_w_m_k) / 1e-4,
            2e-9,
        );
    }
}

#[test]
fn exact_cp_composite_endpoint_keeps_the_original_left_conductivity_tangent() {
    let (mc, pc) = (matrix_card(), pcm_card());
    for capacities in [[1.0, 3.0], [3.0, 1.0]] {
        let (_, pcm) = constituents(mc.content_hash(), pc.content_hash(), [0.0; 2]);
        let matrix = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
            mc.content_hash(),
            SolidLiquidPhase::Solid,
            0.0,
            &[
                HeatCapacityKnot {
                    temperature_k: 300.0,
                    specific_heat_capacity_j_kg_k: capacities[0],
                    bulk_density_kg_m3: 8000.0,
                },
                HeatCapacityKnot {
                    temperature_k: 500.0,
                    specific_heat_capacity_j_kg_k: capacities[1],
                    bulk_density_kg_m3: 8000.0,
                },
                HeatCapacityKnot {
                    temperature_k: 600.0,
                    specific_heat_capacity_j_kg_k: capacities[1],
                    bulk_density_kg_m3: 8000.0,
                },
            ],
        )
        .unwrap();
        let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
        let (mk, pk) = (prepare(&mc, &matrix), prepare(&pc, &pcm));
        let end_h = caloric.enthalpy_range_j_kg()[1];
        let matrix_h = caloric.values_at_specific_enthalpy(end_h)
            .unwrap().matrix.specific_enthalpy_j_kg;
        // The full matrix still selects the hotter/right conductivity slope.
        // The composite ends at 500 K and must use the 400--500 K source slope.
        near(
            mk.state_at_specific_enthalpy(matrix_h).unwrap()
                .conductivity_enthalpy_derivative_kg_w_m_k_j(),
            1.0 / capacities[1],
            1e-14,
        );
        for morphology in [Morphology::Parallel, Morphology::Series] {
            let law = PreparedPcmConductivity::new(&caloric, &mk, &pk, morphology, || Ok(()))
                .unwrap();
            let end = law.state_at_specific_enthalpy(end_h).unwrap();
            // Independent mixture derivative: both phase densities are constant,
            // Cp_PCM=2.5, d(k_matrix)/dT=0.2, and d(k_PCM)/dT=0 here.
            // Using a chord in h instead gives Cp averaged over 400--500 K.
            let phi = (0.25 / 800.0) / (0.75 / 8000.0 + 0.25 / 800.0);
            let cp = 0.75 * capacities[1] + 0.25 * 2.5;
            let (k, dk) = match morphology {
                Morphology::Parallel => (
                    (1.0 - phi) * 50.0 + phi * 3.0,
                    (1.0 - phi) * 0.2 / cp,
                ),
                Morphology::Series => {
                    let k: f64 = 1.0 / ((1.0 - phi) / 50.0 + phi / 3.0);
                    (k, (1.0 - phi) * (k / 50.0).powi(2) * 0.2 / cp)
                }
            };
            near(end.conductivity_w_m_k, k, 1e-13);
            near(end.conductivity_slope_kg_w_m_k_j, dk, 1e-14);
            let delta = 0.01;
            let left = law.state_at_specific_enthalpy(end_h - delta).unwrap();
            let twice_left = law.state_at_specific_enthalpy(end_h - 2.0 * delta).unwrap();
            near(
                end.conductivity_slope_kg_w_m_k_j,
                (3.0 * end.conductivity_w_m_k - 4.0 * left.conductivity_w_m_k
                    + twice_left.conductivity_w_m_k) / (2.0 * delta),
                1e-9,
            );
        }
    }
}

#[test]
fn affine_density_sources_keep_the_original_volume_derivative() {
    let (_, pcm) = constituents(ContentHash([1; 32]), ContentHash([2; 32]), [0.0; 2]);
    let matrix = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([3; 32]),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 8000.0,
            },
            HeatCapacityKnot {
                temperature_k: 500.0,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 4000.0,
            },
        ],
        1e-6,
    )
    .unwrap();
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let h = 62.5;
    let sample = law.response_at_specific_enthalpy(h).unwrap();
    near(sample.state.matrix.bulk_density_kg_m3, 7000.0, 1e-10);
    let a: f64 = 0.75 / 7000.0;
    let b: f64 = 0.25 / 1000.0;
    let da: f64 = 0.75 * 20.0 / (7000.0 * 7000.0) * 0.8;
    near(
        sample.pcm_volume_fraction_slope_kg_j,
        -b * da / (a + b).powi(2),
        1e-15,
    );
    let delta = 1e-3;
    near(
        sample.pcm_volume_fraction_slope_kg_j,
        (law.values_at_specific_enthalpy(h + delta)
            .unwrap()
            .pcm_volume_fraction
            - law
                .values_at_specific_enthalpy(h - delta)
                .unwrap()
                .pcm_volume_fraction)
            / (2.0 * delta),
        1e-11,
    );
}

#[test]
fn original_source_identity_domains_and_inward_zero_conductivity_are_enforced() {
    let (mc, pc) = (matrix_card(), pcm_card());
    let (matrix, pcm) = constituents(mc.content_hash(), pc.content_hash(), [0.0; 2]);
    let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let (mk, pk) = (prepare(&mc, &matrix), prepare(&pc, &pcm));
    let copy = matrix.clone();
    let different = prepare(&mc, &copy);
    assert!(
        PreparedPcmConductivity::new(&caloric, &different, &pk, Morphology::Parallel, || Ok(()))
            .is_err()
    );
    assert!(
        PreparedPcmConductivity::new(&caloric, &mk, &pk, Morphology::Parallel, || Err(
            fs_thermal::enthalpy::DepositionError::MaterialDomain
        ))
        .is_err()
    );
    let law =
        PreparedPcmConductivity::new(&caloric, &mk, &pk, Morphology::Parallel, || Ok(())).unwrap();
    for h in [f64::NAN, -1.0, 313.0] {
        assert!(law.state_at_specific_enthalpy(h).is_err());
    }
    let zero = card(
        "T",
        &[(300.0, 20.0), (400.0, 10.0), (500.0, 0.0), (600.0, 0.0)],
    );
    let (matrix, pcm) = constituents(zero.content_hash(), pc.content_hash(), [0.0; 2]);
    let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let (mk, pk) = (prepare(&zero, &matrix), prepare(&pc, &pcm));
    assert!(matches!(
        PreparedPcmConductivity::new(&caloric, &mk, &pk, Morphology::Series, || Ok(())),
        Err(PhaseConductivityError::Unrepresentable(_))
    ));
    assert!(
        PreparedPcmConductivity::new(&caloric, &mk, &pk, Morphology::Parallel, || Ok(())).is_ok()
    );
}

#[test]
fn native_implicit_conduction_melts_and_freezes_composite_without_extra_heat() {
    let (mc, pc) = (matrix_card(), pcm_card());
    let (matrix, pcm) = constituents(mc.content_hash(), pc.content_hash(), [0.0; 2]);
    let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let (mk, pk) = (prepare(&mc, &matrix), prepare(&pc, &pcm));
    for morphology in [Morphology::Parallel, Morphology::Series] {
        let law = PreparedPcmConductivity::new(&caloric, &mk, &pk, morphology, || Ok(())).unwrap();
        for bath in [350.0, 450.0] {
            let volume = caloric.volume_m3(140.0, 2.0).unwrap();
            let half_distance = 0.5 * volume;
            let dt = 0.01 * half_distance;
            let original = DepositionSimulation::new(
                vec![TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id: CellId(1),
                        material_id: MaterialId(1),
                        body_id: 1,
                        segment_id: 1,
                        mass_kg: 2.0,
                        specific_h_j_kg: 140.0,
                        geometry: CellGeometry {
                            start_m: [0.0; 3],
                            end_m: [volume, 0.0, 0.0],
                            volume_m3: volume,
                            surface_area_m2: 2.0 + 4.0 * volume,
                        },
                        boundary: Boundary {
                            ambient_k: 300.0,
                            htc_w_m2_k: 0.0,
                            source_w: 0.0,
                            reservoir: Some(Reservoir {
                                area_m2: 1.0,
                                conductance_w_k: 0.0,
                                temperature_k: bath,
                            }),
                        },
                    }),
                }],
                &[(MaterialId(1), &caloric)],
                DepositionBudget {
                    max_cells: 1,
                    max_links: 0,
                    max_events: 1,
                    max_substeps: 100,
                },
            )
            .unwrap()
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-9,
                },
                || Ok(()),
            )
            .unwrap()
            .0;
            let observer = ConductionObserver::new(
                &original,
                &[(MaterialId(1), &law)],
                &[],
                &[ConductiveReservoir {
                    cell_id: CellId(1),
                    half_distance_m: half_distance,
                    surface_resistance_m2_k_w: 0.0,
                }],
                || Ok(()),
            )
            .unwrap();
            let chain = ImplicitConductionChain::new(
                &original,
                &observer,
                &[(MaterialId(1), &law)],
                || Ok(()),
            )
            .unwrap();
            let mut work = ImplicitConductionWorkspace::new(1).unwrap();
            let (next, audit) = chain
                .advance_to(
                    &original,
                    dt,
                    ImplicitConductionControls {
                        maximum_dt_s: dt,
                        cell_tolerance_j: 1e-10,
                        total_tolerance_j: 1e-9,
                        ..ImplicitConductionControls::default()
                    },
                    &mut work,
                    || Ok(()),
                )
                .unwrap();
            // Independent scalar extensive-energy root, with analytic source k
            // and phase specific volume; no production constitutive query.
            let (mut lo, mut hi) = (130.0, 150.0);
            for _ in 0..70 {
                let h = 0.5 * (lo + hi);
                let fraction = (h - 125.0) / 50.0;
                let vp = (1.0 - fraction) / 1000.0 + fraction / 800.0;
                let phi = 0.25 * vp / (0.75 / 8000.0 + 0.25 * vp);
                let kp = 2.0 + 8.0 * fraction;
                let k = match morphology {
                    Morphology::Parallel => (1.0 - phi) * 30.0 + phi * kp,
                    Morphology::Series => 1.0 / ((1.0 - phi) / 30.0 + phi / kp),
                };
                if 2.0 * (h - 140.0) - 0.01 * k * (bath - 400.0) > 0.0 {
                    hi = h;
                } else {
                    lo = h;
                }
            }
            let h = next.cells()[&CellId(1)].specific_h_j_kg;
            near(h, 0.5 * (lo + hi), 1e-9);
            assert_eq!(next.temperature_k(CellId(1)).unwrap(), 400.0);
            near(audit.energy.external_heat_j, 2.0 * (h - 140.0), 1e-9);
            near(audit.energy.residual_j, 0.0, 1e-9);
            assert_eq!(next.mass_kg(), 2.0);
            assert_eq!(original.cells()[&CellId(1)].specific_h_j_kg, 140.0);
            assert!((h - 140.0) * (bath - 400.0) > 0.0);
        }
    }
}

#[test]
fn closed_spatial_conduction_melts_pcm_and_conserves_unequal_cell_masses_and_energy() {
    let (mc, pc) = (matrix_card(), pcm_card());
    let (matrix, pcm) = constituents(mc.content_hash(), pc.content_hash(), [0.0; 2]);
    let caloric = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let (mk, pk) = (prepare(&mc, &matrix), prepare(&pc, &pcm));
    let law =
        PreparedPcmConductivity::new(&caloric, &mk, &pk, Morphology::Parallel, || Ok(())).unwrap();
    let masses = [2.0, 1.0];
    let initial_h = [120.0, 225.0];
    // Unit-area adjacent slabs retain their initial reference geometry. The
    // cold PCM is solid and the hot PCM is liquid; both matrices stay solid.
    let volumes = [
        masses[0] * (0.75 / 8000.0 + 0.25 / 1000.0),
        masses[1] * (0.75 / 8000.0 + 0.25 / 800.0),
    ];
    let starts = [0.0, volumes[0]];
    let half_distances = [0.5 * volumes[0], 0.5 * volumes[1]];
    let dt = 0.5 * (half_distances[0] + half_distances[1]);
    let birth = |index: usize| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(index as u64 + 1),
            material_id: MaterialId(1),
            body_id: 1,
            segment_id: 1,
            mass_kg: masses[index],
            specific_h_j_kg: initial_h[index],
            geometry: CellGeometry {
                start_m: [starts[index], 0.0, 0.0],
                end_m: [starts[index] + volumes[index], 0.0, 0.0],
                volume_m3: volumes[index],
                surface_area_m2: 2.0 + 4.0 * volumes[index],
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                source_w: 0.0,
                reservoir: None,
            },
        }),
    };
    let original = DepositionSimulation::new(
        vec![
            birth(0),
            birth(1),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1.0,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        &[(MaterialId(1), &caloric)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &original,
        &[(MaterialId(1), &law)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: half_distances,
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&original, &observer, &[(MaterialId(1), &law)], || Ok(()))
            .unwrap();
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, audit) = chain
        .advance_to(
            &original,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                cell_tolerance_j: 1e-10,
                total_tolerance_j: 1e-9,
                ..ImplicitConductionControls::default()
            },
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(audit.rejected_intervals, 0);

    // Independent two-cell backward-Euler oracle. Q is extensive heat into
    // the cold cell: its first 10 J reach melting, and the remainder melts PCM.
    // The hot cell can release 50 J before it reaches the same 400 K plateau.
    let residual = |q: f64| {
        let fraction = (q / masses[0] - 5.0) / 50.0;
        let pcm_volume = (1.0 - fraction) / 1000.0 + fraction / 800.0;
        let phi_cold = 0.25 * pcm_volume / (0.75 / 8000.0 + 0.25 * pcm_volume);
        let k_cold = (1.0 - phi_cold) * 30.0 + phi_cold * (2.0 + 8.0 * fraction);
        let phi_hot = (0.25 / 800.0) / (0.75 / 8000.0 + 0.25 / 800.0);
        let temperature_gap = (50.0 - q) / 1.375;
        let k_hot = (1.0 - phi_hot) * (30.0 + 0.2 * temperature_gap) + phi_hot * 3.0;
        q - dt * temperature_gap / (half_distances[0] / k_cold + half_distances[1] / k_hot)
    };
    let (mut lo, mut hi) = (10.0, 50.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..70 {
        let q = 0.5 * (lo + hi);
        if residual(q) > 0.0 {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let cold = &next.cells()[&CellId(1)];
    let hot = &next.cells()[&CellId(2)];
    near(cold.specific_h_j_kg, initial_h[0] + q / masses[0], 1e-9);
    near(hot.specific_h_j_kg, initial_h[1] - q / masses[1], 1e-9);
    let cold_state = caloric
        .values_at_specific_enthalpy(cold.specific_h_j_kg)
        .unwrap();
    let hot_state = caloric
        .values_at_specific_enthalpy(hot.specific_h_j_kg)
        .unwrap();
    assert_eq!(cold_state.temperature_k, 400.0);
    assert!(cold_state.pcm.liquid_mass_fraction > 0.0 && cold_state.pcm.liquid_mass_fraction < 1.0);
    assert_eq!(hot_state.pcm.liquid_mass_fraction, 1.0);
    assert!(hot_state.temperature_k > 400.0);
    near(
        masses[0] * cold_state.liquid_mass_fraction * 200.0,
        q - 10.0,
        1e-9,
    );
    for (index, cell) in [cold, hot].into_iter().enumerate() {
        assert_eq!(cell.mass_kg, masses[index]);
        assert_eq!(original.cells()[&cell.id].specific_h_j_kg, initial_h[index]);
    }
    assert_eq!(next.mass_kg(), 3.0);
    assert_eq!(audit.energy.external_heat_j, 0.0);
    near(
        cold.mass_kg * cold.specific_h_j_kg + hot.mass_kg * hot.specific_h_j_kg,
        masses[0] * initial_h[0] + masses[1] * initial_h[1],
        1e-9,
    );
    near(audit.energy.residual_j, 0.0, 1e-9);
}
