//! Source-bound native sheet mechanics and the original caloric timeline.
//! Synthetic fixtures check reciprocal accounts; they do not calibrate forming.

use super::*;
use fs_couple::solid_enthalpy::finite_plastic::sheet::{
    FiniteThermoplasticSheetCoupling, FiniteThermoplasticSheetReport,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_solid::hyper2d::plastic::SheetAssemblyBudget;
use fs_solid::hyper2d::plastic::thermal::FiniteThermoplasticSheet;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::mesh2::Mesh2;
use fs_solid::plastic3::IncrementSettings;

const THICKNESS: f64 = 0.001;
const MU: f64 = 80e6;
const BULK: f64 = 160e6;
const YIELD: f64 = 2e6;
const HARDENING: f64 = 8e6;
const BIAXIAL_COMPLIANCE: f64 = 1.0 / (6.0 * MU) + 2.0 / (9.0 * BULK);

fn sheet_budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 20,
        maximum_elements: 8,
        maximum_free_dofs: 40,
        minimum_scaled_jacobian: 1e-4,
    }
}

fn sheet_settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}

fn panel(triangular: bool) -> Mesh2 {
    if triangular {
        Mesh2 {
            nodes: vec![[0.0, 0.0], [L, 0.0], [0.0, L]],
            elems: vec![vec![0, 1, 2]],
            patches: vec![],
        }
    } else {
        Mesh2::quads(L, L, 1, 1)
    }
}

fn cell(
    id: CellId,
    material_id: MaterialId,
    start_x: f64,
    thickness: f64,
    volume: f64,
    power: f64,
) -> DepositedCell {
    DepositedCell {
        id,
        material_id,
        segment_id: id.0,
        body_id: 7,
        mass_kg: RHO * volume,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [start_x, 0.0, 0.0],
            end_m: [start_x + L, L, thickness],
            volume_m3: volume,
            surface_area_m2: 2.0 * L * L + 4.0 * L * thickness,
        },
        boundary: boundary(power),
    }
}

fn initial_panel<'a>(
    fixture: &'a Fixture,
    volume: f64,
    powers: &[f64],
) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(cell(CELL, MaterialId(0), 0.0, THICKNESS, volume, powers[0])),
    }];
    for (index, &power) in powers.iter().enumerate().skip(1) {
        events.push(TimedEvent {
            time_s: index as f64,
            event: Event::Boundary(CELL, boundary(power)),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: powers.len() + 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn biaxial_supports(mesh: &Mesh2, stretch: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, &[x, y])| {
            [
                (2 * node, (stretch - 1.0) * x),
                (2 * node + 1, (stretch - 1.0) * y),
            ]
        })
        .collect()
}

fn uniaxial_supports(mesh: &Mesh2, stretch: f64) -> Vec<(usize, f64)> {
    let mut values = Vec::new();
    for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
        if x == 0.0 {
            values.push((2 * node, 0.0));
        }
        if y == 0.0 {
            values.push((2 * node + 1, 0.0));
        }
        if x == L {
            values.push((2 * node, (stretch - 1.0) * L));
        }
    }
    values
}

fn constant_source(identity: u32, yield_stress: f64, hardening: f64) -> Fixture {
    let mut source = model();
    for name in [
        "shear_modulus_slope_pa_k",
        "bulk_modulus_slope_pa_k",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
        "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        source.parameters.get_mut(name).unwrap().value = 0.0;
    }
    source.parameters.get_mut("yield_stress_pa").unwrap().value = yield_stress;
    source
        .parameters
        .get_mut("hardening_modulus_pa")
        .unwrap()
        .value = hardening;
    Fixture::with_model(
        CP,
        identity,
        scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
        source,
    )
}

fn biaxial_reference(axial_log: f64, yield_stress: f64, hardening: f64) -> (f64, f64) {
    // tau = diag(q,q,0), e = C*q + alpha/2, q = Y + H*alpha.
    let alpha =
        (axial_log - BIAXIAL_COMPLIANCE * yield_stress) / (0.5 + BIAXIAL_COMPLIANCE * hardening);
    assert!(alpha > 0.0);
    (alpha, yield_stress + hardening * alpha)
}

fn closed_accounts(report: &FiniteThermoplasticSheetReport) {
    assert!(report.inertia.is_none());
    assert!(report.heat_defect_j <= config().energy_tolerance_j);
    close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
    close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
    close(
        report.material_energy.dead_load_endpoint_work_j
            + report.material_energy.support_endpoint_work_j
            + report.external_heat_j,
        report.material_energy.internal_energy_change_j
            + report.thermal_audit.stored_energy_change_j
            + report.material_energy.endpoint_work_defect_j,
        2e-8,
    );
}

#[test]
fn original_source_p1_and_q1_sheets_recover_heated_forming_and_free_lateral_motion() {
    context(|cx| {
        let fixture = Fixture::new(CP, 40);
        let material = FiniteThermoplasticMaterial::try_new(fixture.query(), 0.0).unwrap();
        let family = material.family().law();
        for triangular in [false, true] {
            let mesh = panel(triangular);
            let thickness = [THICKNESS];
            let volume = L * L * THICKNESS * if triangular { 0.5 } else { 1.0 };
            let mut expected = Vec::new();
            let mut old_state = family.initial_state();
            let mut old_f = IDENTITY;
            let mut old_t = 300.0;
            for (temperature, axial_log) in [(320.0, 0.10_f64), (340.0, 0.16)] {
                let p = family.parameters();
                let delta = temperature - p.reference_temperature_k;
                let mu = p.reference.shear_modulus_pa + delta * p.slopes.shear_modulus_pa_k;
                let bulk = p.reference.bulk_modulus_pa + delta * p.slopes.bulk_modulus_pa_k;
                let yield_stress = p.reference.yield_stress_pa + delta * p.slopes.yield_stress_pa_k;
                let hardening =
                    p.reference.hardening_modulus_pa + delta * p.slopes.hardening_modulus_pa_k;
                let young = 9.0 * bulk * mu / (3.0 * bulk + mu);
                let poisson = (3.0 * bulk - 2.0 * mu) / (2.0 * (3.0 * bulk + mu));
                let log_theta = p.logarithmic_expansion_per_k * delta
                    + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta;
                // Independent coaxial Hencky/J2 uniaxial-stress endpoint.
                let alpha =
                    (axial_log - log_theta - yield_stress / young) / (1.0 + hardening / young);
                let q = yield_stress + hardening * alpha;
                let lateral = (log_theta - poisson * q / young - 0.5 * alpha).exp();
                let f = [
                    axial_log.exp(),
                    0.0,
                    0.0,
                    0.0,
                    lateral,
                    0.0,
                    0.0,
                    0.0,
                    lateral,
                ];
                let increment = family
                    .checked_increment_for_step(
                        &old_f,
                        &old_state,
                        old_t,
                        &f,
                        temperature,
                        1e-4,
                        1.0,
                    )
                    .unwrap();
                close(increment.mechanical.state.alpha, alpha, 3e-10);
                assert!(increment.mechanical.state.alpha > old_state.alpha);
                let external = RHO * volume * CP * (temperature - old_t)
                    - volume * increment.caloric_heat_j_m3();
                old_state = increment.mechanical.state.clone();
                old_f = f;
                old_t = temperature;
                expected.push((temperature, f, alpha, increment, external));
            }
            let powers: Vec<_> = expected.iter().map(|value| value.4).collect();
            let mut thermal = initial_panel(&fixture, volume, &powers);
            let native = FiniteThermoplasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&material),
                &thickness,
                &[300.0],
                PlaneStressSettings::default(),
                sheet_budget(),
                cx,
            )
            .unwrap();
            let mut sheet =
                FiniteThermoplasticSheetCoupling::try_new(&thermal, native, &[CELL], cx).unwrap();
            for (index, (temperature, f, alpha, increment, external)) in expected.iter().enumerate()
            {
                let old_history = sheet.mechanics().quadrature_states().to_vec();
                let supports = uniaxial_supports(&mesh, f[0]);
                let candidate = sheet
                    .advance(
                        &thermal,
                        &(),
                        (index + 1) as f64,
                        &supports,
                        &[],
                        sheet_settings(),
                        controls().material,
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let actual = candidate.sheet().mechanics();
                let report = candidate.report();
                close(
                    candidate.thermal().temperature_k(CELL).unwrap(),
                    *temperature,
                    2e-7,
                );
                for &point_t in actual.quadrature_temperatures_k() {
                    close(point_t, *temperature, 2e-7);
                }
                for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                    close(actual.displacement_m()[node][0], (f[0] - 1.0) * x, 3e-11);
                    close(actual.displacement_m()[node][1], (f[4] - 1.0) * y, 3e-11);
                }
                for (state, actual_f) in actual
                    .quadrature_states()
                    .iter()
                    .zip(actual.quadrature_deformation())
                {
                    close(state.alpha, *alpha, 4e-9);
                    for (value, expected_f) in actual_f.iter().zip(f) {
                        close(*value, *expected_f, 3e-9);
                    }
                }
                for recovered in actual.quadrature_thickness_m() {
                    close(recovered, THICKNESS * f[8], 4e-12);
                }
                let FinitePlasticLawField::Uniform(bound) = actual.families() else {
                    panic!("original uniform source binding changed");
                };
                assert!(std::ptr::eq(bound, &material));
                assert!(std::ptr::eq(bound.caloric_source(), &fixture.caloric));
                assert_eq!(
                    bound.query().card.content_hash(),
                    fixture.card.content_hash()
                );
                assert_eq!(bound.family().model_card(), &fixture.card.models()[0]);
                close(
                    candidate.thermal().cells()[&CELL].mass_kg,
                    RHO * volume,
                    0.0,
                );
                close(
                    report.material_energy.plastic_dissipation_j,
                    volume * increment.mechanical.plastic_dissipation_j_m3,
                    2e-9,
                );
                close(
                    report.material_energy.reversible_heat_j,
                    volume * increment.reversible_heat_j_m3,
                    2e-9,
                );
                close(
                    report.material_energy.caloric_heat_j,
                    volume * increment.caloric_heat_j_m3(),
                    2e-9,
                );
                close(
                    report.material_energy.internal_energy_change_j,
                    volume * increment.internal_energy_change_j_m3(),
                    2e-9,
                );
                close(
                    report.material_energy.endpoint_work_defect_j,
                    volume * increment.endpoint_work_defect_j_m3,
                    2e-9,
                );
                close(report.external_heat_j, *external, 2e-9);
                assert!(report.equilibrium.residual_history.len() > 1);
                assert!(report.iterations > 1);
                closed_accounts(report);
                assert_eq!(sheet.mechanics().quadrature_states(), old_history);
                close(sheet.time_s(), index as f64, 0.0);
                close(thermal.time_s(), index as f64, 0.0);
                let (next_thermal, next_sheet, (), _) = candidate.into_parts();
                thermal = next_thermal;
                sheet = next_sheet;
            }
        }
    });
}

#[test]
fn plastic_self_heating_excludes_hardening_storage_and_endpoint_work_defect() {
    context(|cx| {
        let fixture = constant_source(41, YIELD, HARDENING);
        let material = FiniteThermoplasticMaterial::try_new(fixture.query(), 0.0).unwrap();
        let mesh = panel(false);
        let thickness = [THICKNESS];
        let volume = L * L * THICKNESS;
        let thermal = initial_panel(&fixture, volume, &[0.0]);
        let native = FiniteThermoplasticSheet::try_new(
            &mesh,
            FinitePlasticLawField::Uniform(&material),
            &thickness,
            &[300.0],
            PlaneStressSettings::default(),
            sheet_budget(),
            cx,
        )
        .unwrap();
        let sheet =
            FiniteThermoplasticSheetCoupling::try_new(&thermal, native, &[CELL], cx).unwrap();
        let axial: f64 = 0.12;
        let (alpha, q) = biaxial_reference(axial, YIELD, HARDENING);
        let heat = volume * YIELD * alpha;
        let hardening_storage = volume * 0.5 * HARDENING * alpha * alpha;
        let storage = volume * BIAXIAL_COMPLIANCE * q * q + hardening_storage;
        let temperature = 300.0 + YIELD * alpha / (RHO * CP);
        let candidate = sheet
            .advance(
                &thermal,
                &(),
                1.0,
                &biaxial_supports(&mesh, axial.exp()),
                &[],
                sheet_settings(),
                controls().material,
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let actual = candidate.sheet().mechanics();
        let report = candidate.report();
        close(
            candidate.thermal().temperature_k(CELL).unwrap(),
            temperature,
            2e-7,
        );
        close(actual.temperatures_k()[0], temperature, 2e-7);
        for state in actual.quadrature_states() {
            close(state.alpha, alpha, 4e-9);
        }
        close(actual.free_energy_j(), storage, 2e-9);
        close(
            report.material_energy.internal_energy_change_j,
            storage,
            2e-9,
        );
        close(report.material_energy.plastic_dissipation_j, heat, 2e-9);
        close(report.material_energy.caloric_heat_j, heat, 2e-9);
        close(report.material_energy.reversible_heat_j, 0.0, 1e-12);
        close(report.material_energy.creep_dissipation_j, 0.0, 1e-12);
        close(
            report.material_energy.incremental_potential_j,
            storage + heat,
            2e-9,
        );
        close(report.thermal_audit.stored_energy_change_j, heat, 2e-9);
        close(report.applied_heat_j, heat, 2e-9);
        close(report.external_heat_j, 0.0, 2e-9);
        assert!(hardening_storage > 0.01);
        assert!(report.material_energy.endpoint_work_defect_j.abs() > 0.005);
        let fictitious_temperature = temperature
            + (hardening_storage + report.material_energy.endpoint_work_defect_j)
                / (RHO * volume * CP);
        assert!((fictitious_temperature - temperature).abs() > 0.2);
        closed_accounts(report);
        assert_eq!(sheet.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(
            sheet
                .mechanics()
                .quadrature_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
    });
}
