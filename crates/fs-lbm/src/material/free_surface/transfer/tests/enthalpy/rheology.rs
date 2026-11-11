//! Original P1 motion, source rheology and heat survive one native transfer.
//! Manufactured source tables isolate the equations, not measured materials.

use super::*;

#[cfg(feature = "enthalpy-coupling")]
mod source_bulk_feedback;
use crate::d3q19::freesurface3::FreeSurface3;
use crate::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use crate::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_qty::Frequency;

pub(super) fn source_card(rate_slope: f64, maximum_rate: f64) -> MaterialCard {
    let original = thermal_card(1000.0, 2e-5);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let viscosity = claim.key.name() == "dynamic_viscosity";
        if viscosity || claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [300.0, 0.0],
                    [500.0, 0.0],
                    [300.0, maximum_rate],
                    [500.0, maximum_rate],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        if viscosity {
                            viscosity_at(p[0], p[1], rate_slope)
                        } else {
                            0.5 + 0.001 * (p[0] - 400.0) + 0.01 * p[1]
                        },
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claim.validity = claim.validity.with_quantity(
            "rate",
            QuantitySpec::dimensional(Frequency::DIMS),
            0.0,
            maximum_rate,
        );
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn viscosity_at(temperature: f64, rate: f64, slope: f64) -> f64 {
    0.1 - 0.0004 * (temperature - 300.0) + slope * rate
}

fn shear_point(rate: f64) -> QueryPoint {
    point()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), rate)
        .unwrap()
}

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

fn bounds() -> TetEnthalpyTransferLimits {
    TetEnthalpyTransferLimits {
        max_source_queries: 512,
        ..thermal_limits()
    }
}

#[test]
fn original_p1_shear_law_drives_receiving_motion_and_molecular_heat() {
    let (p, t, _) = cube([0.0; 3]);
    let v: Vec<_> = p
        .iter()
        .map(|x| [0.01 + 2.0 * (x[1] - 0.004), 0.0, 0.0])
        .collect();
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let before = (p.clone(), t.clone(), v.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    let mut outcomes = Vec::new();
    with_cx(&gate, |cx| {
        for slope in [0.0, 0.02] {
            let card = source_card(slope, 100.0);
            // The authored reference differs from the actual P1 rate (2).
            // Receiving stress must retain that rate, not restart at zero.
            let law = LiquidTemperatureLaw::try_new(
                &card,
                &shear_point(1.5),
                "T",
                LiquidViscosityBasis::Dynamic,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                law.reference_state().material().clone(),
                LiquidViscosityBasis::Dynamic,
                "T",
                "rate",
            )
            .unwrap();
            let caloric = chart(&card, 1000.0);
            let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
            assert!(
                scaling
                    .free_surface3_from_tetrahedra_with_enthalpy(
                        domain(),
                        [0.0; 3],
                        source(&p, &t, &v, &h),
                        &law,
                        &caloric,
                        bounds(),
                        Default::default(),
                        cx,
                        &gate,
                    )
                    .is_err()
            );
            let mut receiving = domain();
            receiving.gravity_m_s2 = [0.5, -0.75, 0.25];
            let (thermal, report) = scaling
                .free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
                    receiving,
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    solve(),
                    cx,
                    &gate,
                )
                .unwrap();
            let mut thermal = thermal.with_time_origin(1.25).unwrap();
            assert_eq!(report.source_shear_rates_per_s.len(), t.len());
            for rate in &report.source_shear_rates_per_s {
                close(*rate, 2.0, 1e-13);
            }
            assert!(report.source_queries > 4 * t.len());
            assert_eq!(report.mechanical.components.len(), 1);
            assert_eq!(report.components.len(), 1);
            let mass = 1000.0 * 0.004_f64.powi(3);
            let enthalpy = mass * affine_h([0.004; 3]);
            close(report.mechanical.source.mass_kg, mass, 1e-17);
            close(report.mechanical.target.mass_kg, mass, 1e-17);
            close(report.source_enthalpy_j, enthalpy, 1e-12);
            close(report.target_enthalpy_j, enthalpy, 1e-11);
            close(report.components[0].source_enthalpy_j, enthalpy, 1e-12);
            close(report.components[0].target_enthalpy_j, enthalpy, 1e-11);
            for axis in 0..3 {
                close(
                    report.mechanical.target.momentum_kg_m_s[axis],
                    report.mechanical.source.momentum_kg_m_s[axis],
                    1e-16,
                );
                close(
                    report.mechanical.target.angular_momentum_kg_m2_s[axis],
                    report.mechanical.source.angular_momentum_kg_m2_s[axis],
                    1e-18,
                );
            }
            close(thermal.rheology_at([3, 3, 3]).unwrap().shear_rate_per_s, 2.0, 1e-8);
            let initial = thermal.material_at([3, 3, 3]).unwrap();
            let temperature = initial.material().query().axes()["T"];
            close(initial.dynamic_viscosity_pa_s(), viscosity_at(temperature, 2.0, slope), 1e-9);
            let mut heat = 0.0;
            let mut maximum_native_rate = 0.0_f64;
            for _ in 0..4 {
                let original = thermal.clone();
                let step = thermal.step_with_viscous_heating(&[], &gate).unwrap();
                let mut expected_heat = 0.0;
                for sample in &step.applied {
                    let point = sample.material.material().query();
                    let temperature = point.axes()["T"];
                    let query_rate = point.axes()["rate"];
                    close(query_rate, sample.shear_rate_per_s, 1e-7);
                    let mu = viscosity_at(temperature, query_rate, slope);
                    close(sample.material.dynamic_viscosity_pa_s(), mu, 1e-11);
                    let location = FreeSurface3::transport_location(sample.coordinates).unwrap();
                    let tau = original.rheology_at(sample.coordinates).unwrap().collision_tau;
                    let dilation = original.flow().lattice()
                        .transport_strain_rates(location, tau).unwrap()[1] / 0.0005;
                    expected_heat +=
                        mu * (sample.shear_rate_per_s.powi(2) + dilation * dilation)
                            * sample.liquid_volume_m3 * 0.0005;
                    maximum_native_rate = maximum_native_rate.max(sample.shear_rate_per_s);
                    assert_eq!(
                        sample.material.material().card_identity(),
                        card.content_hash()
                    );
                }
                close(step.flow.viscous_heat_j, expected_heat, 1e-12);
                close(step.flow.external_heat_j, 0.0, 0.0);
                close(step.flow.balance_defect_j, 0.0, 1e-9);
                heat += step.flow.viscous_heat_j;
            }
            assert!(maximum_native_rate > 1e-3);
            assert!(heat > 0.0);
            close(thermal.total_enthalpy_j().unwrap(), enthalpy + heat, 1e-8);
            close(thermal.physical_time_s(), 1.252, 1e-14);
            let native = thermal.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
            outcomes.push((native.velocity_m_s, heat));
        }
    });
    assert!(
        outcomes[0]
            .0
            .iter()
            .zip(outcomes[1].0)
            .any(|(a, b)| (a - b).abs() > 1e-10)
    );
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 1e-13);
    assert_eq!((p, t, v, h), before);
}

#[test]
fn p1_transfer_retains_independent_bulk_work_without_changing_the_shear_law() {
    let card = source_card(0.02, 100.0);
    let law = LiquidTemperatureLaw::try_new(
        &card, &shear_point(1.5), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(), LiquidViscosityBasis::Dynamic, "T", "rate",
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (points, tets, _) = cube([0.0; 3]);
    let h = vec![200000.0; points.len()];
    let volume = 0.004_f64.powi(3);
    let timestep = 0.0005;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for shear_rate in [0.0, 20.0] {
            // The full P1 field contains both spin and isotropic expansion.
            // Neither changes the original generalized-Newtonian shear query.
            let velocities: Vec<_> = points.iter().map(|point| {
                let r = sub(*point, [0.004; 3]);
                [0.01 + shear_rate * r[1] + 5.0 * r[0],
                    -0.004 + 5.0 * r[1], 0.002 + 5.0 * r[2]]
            }).collect();
            for bulk_viscosity in [0.02, 0.24] {
                let scaling = plan_free_surface_from_material(
                    &surface, 0.001, timestep, 0.1,
                ).unwrap().with_bulk_viscosity(bulk_viscosity, 1.2).unwrap();
                let mut receiving = domain();
                receiving.gravity_m_s2 = [0.5, -0.75, 0.25];
                let (mut thermal, report) = scaling
                    .free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
                        receiving, [0.0; 3], source(&points, &tets, &velocities, &h),
                        &law, &caloric, bounds(), Default::default(), solve(), cx, &gate,
                    ).unwrap();
                assert_eq!(thermal.flow().scaling().bulk_viscosity_pa_s(), Some(bulk_viscosity));
                for &actual in &report.source_shear_rates_per_s {
                    close(actual, shear_rate, 1e-11);
                }
                close(report.mechanical.source.mass_kg, 1000.0 * volume, 1e-17);
                close(report.mechanical.target.mass_kg, 1000.0 * volume, 1e-17);
                close(report.source_enthalpy_j, 200000.0 * 1000.0 * volume, 1e-11);
                close(report.target_enthalpy_j, report.source_enthalpy_j, 1e-11);
                for axis in 0..3 {
                    close(report.mechanical.target.momentum_kg_m_s[axis],
                        report.mechanical.source.momentum_kg_m_s[axis], 1e-16);
                    close(report.mechanical.target.angular_momentum_kg_m2_s[axis],
                        report.mechanical.source.angular_momentum_kg_m2_s[axis], 1e-18);
                }
                let expected_mu = viscosity_at(400.0, shear_rate, 0.02);
                for component in &report.mechanical.components {
                    for &cell in &component.target_cells {
                        close(thermal.rheology_at(cell).unwrap().shear_rate_per_s,
                            shear_rate, 1e-8);
                        close(thermal.material_at(cell).unwrap().dynamic_viscosity_pa_s(),
                            expected_mu, 1e-9);
                    }
                }
                let collision = thermal.flow().lattice().collision_model();
                let step = thermal.step_with_viscous_heating(&[], &gate).unwrap();
                for sample in &step.applied {
                    close(sample.shear_rate_per_s, shear_rate, 1e-8);
                    close(sample.divergence_magnitude_per_s, 15.0, 1e-8);
                    close(sample.shear_power_w_m3,
                        expected_mu * shear_rate * shear_rate, 1e-7);
                    close(sample.bulk_power_w_m3, bulk_viscosity * 225.0, 1e-7);
                }
                // Independent work integral over the original affine field.
                // The heat owner must retain div(u)=15/s when tau_bulk differs
                // from both reference and locally resolved molecular tau.
                let expected_heat = (expected_mu * shear_rate * shear_rate
                    + bulk_viscosity * 15.0_f64.powi(2)) * volume * timestep;
                close(step.flow.viscous_heat_j, expected_heat, expected_heat * 1e-7);
                close(step.flow.external_heat_j, 0.0, 0.0);
                close(step.flow.balance_defect_j, 0.0, 2e-11);
                close(thermal.total_enthalpy_j().unwrap(),
                    report.target_enthalpy_j + expected_heat, 2e-11);
                assert_eq!(thermal.flow().lattice().collision_model(), collision);
            }
        }
    });
}

fn sourced_bulk_card(cold_vertex_refusal: bool) -> MaterialCard {
    let original = source_card(0.02, 100.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let dims = fs_qty::DynViscosity::DIMS;
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("bulk_eta", dims),
        value: PropertyValue::Surface2D {
            axes: ["T".into(), "rate".into()],
            axis_dims: [Temperature::DIMS, Frequency::DIMS],
            samples: [[300.0, 0.0], [500.0, 0.0], [300.0, 100.0], [500.0, 100.0]]
                .into_iter().map(|p| (p, if cold_vertex_refusal {
                    // Reference T/rate and every receiving cell average are
                    // positive, but the coldest original P1 vertex is not.
                    0.001 * (p[0] - 393.0) + 0.0001 * (p[1] - 20.0)
                } else {
                    sourced_bulk_at(p[0], p[1])
                })).collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]], dims,
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 100.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured local volumetric viscosity".into(),
            license: "test-only".into(), artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn sourced_bulk_at(temperature: f64, shear: f64) -> f64 {
    0.03 + 0.001 * (temperature - 300.0) + 0.002 * shear
}

fn sourced_bulk_law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    use fs_material::state_point::{
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    };
    let reference = shear_point(1.5);
    let transport = LiquidTemperatureLaw::try_new(
        card, &reference, "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let mut requirements: Vec<_> = transport.reference_state().material().properties()
        .iter().map(|property| property.requirement().clone()).collect();
    requirements.push(ScalarPropertyRequirement::try_new(
        "bulk_eta", fs_qty::DynViscosity::DIMS, ScalarAdmissibility::StrictlyPositive,
    ).unwrap());
    let bundle = resolve_material_state_point(
        card, &reference, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    LiquidTemperatureLaw::try_from_material(card, bundle, "T", LiquidViscosityBasis::Dynamic).unwrap()
}

#[test]
fn p1_transfer_seeds_source_local_bulk_stress_and_its_first_physical_heat() {
    let card = sourced_bulk_card(false);
    let law = sourced_bulk_law(&card);
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(), LiquidViscosityBasis::Dynamic, "T", "rate",
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (points, tets, _) = cube([0.0; 3]);
    let h: Vec<_> = points.iter().copied().map(affine_h).collect();
    let shear = 20.0;
    let velocities: Vec<_> = points.iter().map(|point| {
        let r = sub(*point, [0.004; 3]);
        [0.01 + shear * r[1] + 5.0 * r[0], -0.004 + 5.0 * r[1], 0.002 + 5.0 * r[2]]
    }).collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1)
            .unwrap().with_bulk_viscosity_source("bulk_eta", 1.2).unwrap();
        let mut receiving = domain();
        receiving.gravity_m_s2 = [0.5, -0.75, 0.25];
        let (mut liquid, report) = scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
            receiving, [0.0; 3], source(&points, &tets, &velocities, &h),
            &law, &caloric, bounds(), Default::default(), solve(), cx, &gate,
        ).unwrap();
        for &rate in &report.source_shear_rates_per_s { close(rate, shear, 1e-11); }
        close(report.target_enthalpy_j, report.source_enthalpy_j, 1e-11);
        let mut expected_heat = 0.0;
        let mut minimum_bulk = f64::INFINITY;
        let mut maximum_bulk = f64::NEG_INFINITY;
        for component in &report.mechanical.components {
            for &cell in &component.target_cells {
                let state = liquid.material_at(cell).unwrap();
                let temperature = state.material().state_coordinate("T").unwrap();
                let bulk = state.material().property("bulk_eta").unwrap().value_si();
                close(bulk, sourced_bulk_at(temperature, shear), 1e-10);
                minimum_bulk = minimum_bulk.min(bulk);
                maximum_bulk = maximum_bulk.max(bulk);
                close(liquid.rheology_at(cell).unwrap().shear_rate_per_s, shear, 1e-8);
                let native = liquid.flow().sample_at_si(cell[0], cell[1], cell[2]).unwrap().unwrap();
                let volume = native.density_kg_m3 / 1000.0 * native.fill_fraction * 0.001_f64.powi(3);
                expected_heat += (viscosity_at(temperature, shear, 0.02) * shear * shear
                    + bulk * 15.0_f64.powi(2)) * volume * 0.0005;
            }
        }
        assert!(maximum_bulk - minimum_bulk > 0.002);
        let step = liquid.step_with_viscous_heating(&[], &gate).unwrap();
        for sample in &step.applied {
            let temperature = sample.material.material().state_coordinate("T").unwrap();
            close(sample.shear_rate_per_s, shear, 1e-8);
            close(sample.divergence_magnitude_per_s, 15.0, 1e-8);
            close(sample.bulk_power_w_m3, sourced_bulk_at(temperature, shear) * 225.0, 1e-7);
        }
        close(step.flow.viscous_heat_j, expected_heat, expected_heat * 1e-7);
        close(step.flow.external_heat_j, 0.0, 0.0);
        close(step.flow.balance_defect_j, 0.0, 2e-11);
        close(liquid.total_enthalpy_j().unwrap(), report.target_enthalpy_j + expected_heat, 2e-11);
    });
}

#[test]
fn local_bulk_admission_checks_original_p1_vertices_before_averaging() {
    let card = sourced_bulk_card(true);
    let law = sourced_bulk_law(&card);
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(), LiquidViscosityBasis::Dynamic, "T", "rate",
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (points, tets, _) = cube([0.0; 3]);
    let h: Vec<_> = points.iter().copied().map(affine_h).collect();
    let velocities: Vec<_> = points.iter().map(|point| [20.0 * (point[1] - 0.004), 0.0, 0.0]).collect();
    let before = (points.clone(), tets.clone(), velocities.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1)
            .unwrap().with_bulk_viscosity_source("bulk_eta", 1.2).unwrap();
        assert!(scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
            domain(), [0.0; 3], source(&points, &tets, &velocities, &h),
            &law, &caloric, bounds(), Default::default(), solve(), cx, &gate,
        ).is_err());
    });
    assert_eq!((points, tets, velocities, h), before);
}

#[test]
fn g4_known_initial_tau_keeps_a_valid_strain_inside_its_tight_source_domain() {
    use crate::material::free_surface::enthalpy::rheology::RheologyBinding;

    let card = source_card(0.02, 2.001);
    let law = LiquidTemperatureLaw::try_new(
        &card, &shear_point(1.5), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(), LiquidViscosityBasis::Dynamic, "T", "rate",
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (points, tets, _) = cube([0.0; 3]);
    let velocities: Vec<_> = points.iter()
        .map(|p| [0.01 + 2.0 * (p[1] - 0.004), 0.0, 0.0]).collect();
    let h = vec![200000.0; points.len()];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
        let mut receiving = domain();
        receiving.gravity_m_s2 = [0.5, -0.75, 0.25];
        let (thermal, report) = scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
            receiving, [0.0; 3], source(&points, &tets, &velocities, &h),
            &law, &caloric, bounds(), Default::default(), solve(), cx, &gate,
        ).unwrap();
        for component in &report.mechanical.components {
            for &cell in &component.target_cells {
                close(thermal.rheology_at(cell).unwrap().shear_rate_per_s, 2.0, 1e-8);
            }
        }
        close(report.target_enthalpy_j, report.source_enthalpy_j, 1e-11);
        // The default attachment must retain its existing semantics. Its
        // reference tau=0.635 estimates 2.047.../s from this valid stress,
        // outside the source domain; the admitted initial tau=0.65 avoids that
        // false query while the SAME source solve still checks the endpoint.
        let binding = RheologyBinding::try_new(&surface, &law, solve()).unwrap();
        assert!(thermal.flow().clone().with_bound_thermal_transport(
            &law, &caloric, Default::default(),
            |cell| thermal.specific_enthalpy_at(cell).unwrap(),
            None, Some(binding), None, &gate,
        ).is_err());
    });
}

#[test]
fn actual_p1_rate_not_reference_rate_controls_source_admission() {
    let card = source_card(0.02, 1.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &shear_point(0.5),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(),
        LiquidViscosityBasis::Dynamic,
        "T",
        "rate",
    )
    .unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, rotation) = cube([0.0; 3]);
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let shear: Vec<_> = p
        .iter()
        .map(|x| [0.01 + 2.0 * (x[1] - 0.004), 0.0, 0.0])
        .collect();
    let before = (p.clone(), rotation.clone(), shear.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
        let call = |v: &[[f64; 3]], limits, controls, rheology| {
            scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
                domain(),
                [0.0; 3],
                source(&p, &t, v, &h),
                &law,
                &caloric,
                limits,
                controls,
                rheology,
                cx,
                &gate,
            )
        };
        assert!(matches!(
            call(&shear, bounds(), Default::default(), solve()),
            Err(TetEnthalpyTransferError::Enthalpy(
                FreeSurfaceEnthalpyError::Material(_)
            ))
        ));
        let (_, report) = call(&rotation, bounds(), Default::default(), solve()).unwrap();
        assert!(
            report
                .source_shear_rates_per_s
                .iter()
                .all(|rate| *rate < 1e-12)
        );
        let mut short = bounds();
        short.max_source_queries = 2;
        assert!(call(&rotation, short, Default::default(), solve()).is_err());
        assert!(
            call(
                &rotation,
                bounds(),
                Default::default(),
                RheologySolve {
                    maximum_iterations: 0,
                    ..solve()
                }
            )
            .is_err()
        );
        // Admission of actual P1 sources has passed, then the receiving native
        // thermal owner refuses atomically without dropping its shear law.
        let mut late = FreeSurfaceEnthalpyControls::default();
        late.max_cells = 1;
        assert!(call(&rotation, bounds(), late, solve()).is_err());
        gate.request();
        assert!(matches!(
            call(&rotation, bounds(), Default::default(), solve()),
            Err(TetEnthalpyTransferError::Transfer(
                TetLiquidTransferError::Cancelled
            ))
        ));
    });
    assert_eq!((p, rotation, shear, h), before);
}
