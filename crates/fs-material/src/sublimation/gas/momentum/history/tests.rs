//! Original finite gas/calorics with independent load and thermal-history oracles.
use super::*;
use crate::sublimation::gas::{
    GasSurfaceBinding,
    chamber::{
        GasChamberSolve,
        advance::{advance_chamber_adaptive, partitioned::advance_partitioned_chamber_adaptive},
    },
    tests as fixture,
};
use fs_thermal::sublimation::{
    SolidSurfaceState, SublimationError, SublimationLaw, SublimationProperties,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn controls() -> MomentumChamberAdvanceSolve {
    MomentumChamberAdvanceSolve {
        chamber: GasChamberAdvanceSolve {
            step: GasChamberSolve {
                surface: fixture::controls(),
                maximum_surfaces: 8,
            },
            absolute_mass_tolerance_kg: 1e-10,
            absolute_temperature_tolerance_k: 1e-4,
            relative_tolerance: 0.0,
            minimum_dt_s: 1e-9,
            maximum_dt_s: 0.1,
            maximum_trials: 10000,
        },
        absolute_impulse_tolerance_n_s: 1e-6,
        absolute_angular_impulse_tolerance_n_m_s: 1e-6,
        relative_impulse_tolerance: 0.0,
    }
}
fn row<'a>(
    id: u64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState,
    temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasChamberSurface<'a> {
    GasChamberSurface {
        id,
        law,
        binding,
        solid,
        temperature_from_enthalpy: temperature,
    }
}
fn frame(id: u64) -> MolecularSurfaceFrame {
    MolecularSurfaceFrame {
        surface_id: id,
        outward_normal: [0.0, 0.0, 1.0],
        centroid_m: [2.0, 0.0, 0.0],
    }
}
fn states<'a>(
    templates: &[GasChamberSurface<'a>],
    endpoint: &GasChamberAdvanceCandidate<'_>,
) -> Vec<GasChamberSurface<'a>> {
    templates
        .iter()
        .map(|template| GasChamberSurface {
            solid: endpoint
                .endpoint()
                .surfaces()
                .iter()
                .find(|s| s.id == template.id)
                .unwrap()
                .solid,
            ..*template
        })
        .collect()
}
fn maxwell(catalog: &ThermochemicalGasModel, fraction: f64) -> MaxwellReflection {
    MaxwellReflection::from_diffuse_fractions(
        catalog,
        &[
            (SpeciesId::new("A").unwrap(), fraction),
            (fixture::species(), fraction),
        ],
        [200.0, 900.0],
    )
    .unwrap()
}

#[test]
fn original_adaptive_path_and_sources_are_retained_and_only_fine_loads_are_integrated() {
    struct Counted {
        law: fixture::Phase,
        calls: std::cell::Cell<usize>,
    }
    impl SublimationLaw for Counted {
        fn temperature_range_k(&self) -> [f64; 2] {
            self.law.temperature_range_k()
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            self.calls.set(self.calls.get() + 1);
            self.law.properties(t)
        }
    }
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let law = Counted {
        law: fixture::phase(),
        calls: std::cell::Cell::new(0),
    };
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let rows = [
        row(
            20,
            &law,
            &binding,
            fixture::solid(420.0, chart),
            &temperature,
        ),
        row(
            10,
            &law,
            &binding,
            fixture::solid(380.0, chart),
            &temperature,
        ),
    ];
    let frames = [
        frame(20),
        MolecularSurfaceFrame {
            outward_normal: [0.0, 0.0, -1.0],
            ..frame(10)
        },
    ];
    let mut c = controls();
    c.chamber.maximum_dt_s = 0.001;
    c.chamber.absolute_mass_tolerance_kg = 1.0;
    c.chamber.absolute_temperature_tolerance_k = 1e6;
    c.absolute_impulse_tolerance_n_s = 1e20;
    c.absolute_angular_impulse_tolerance_n_m_s = 1e20;
    let original =
        advance_chamber_adaptive(&rows, &gas, &catalog, 0.001, c.chamber, || false).unwrap();
    let calls = law.calls.replace(0);
    let result = advance_chamber_adaptive_with_momentum(
        &rows,
        &frames,
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.history(), &original);
    assert_eq!(
        law.calls.get(),
        calls,
        "history replayed a constitutive phase source"
    );
    assert_eq!(result.history().audit().accepted_intervals, 1);
    let first = super::super::advance_chamber_with_momentum(
        &rows,
        &frames,
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.0005,
        c.chamber.step,
        || false,
    )
    .unwrap();
    let middle: Vec<_> = rows
        .iter()
        .map(|r| GasChamberSurface {
            solid: first
                .exchange()
                .surfaces()
                .iter()
                .find(|s| s.id == r.id)
                .unwrap()
                .solid,
            ..*r
        })
        .collect();
    let second = super::super::advance_chamber_with_momentum(
        &middle,
        &frames,
        [0.0; 3],
        None,
        first.exchange().gas(),
        &catalog,
        0.0005,
        c.chamber.step,
        || false,
    )
    .unwrap();
    assert_eq!(result.history().endpoint(), second.exchange());
    for (id, integrated) in &result.momentum().surfaces {
        let a = first.momentum().surfaces[id];
        let b = second.momentum().surfaces[id];
        for ((whole, a), b) in normal(integrated)
            .into_iter()
            .zip(normal(&a))
            .zip(normal(&b))
        {
            fixture::near(whole, a + b, 2e-18);
        }
    }
    let reverse = advance_chamber_adaptive_with_momentum(
        &[rows[1], rows[0]],
        &[frames[1], frames[0]],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result, reverse);
}

// A sealed interface exchanges only carrier heat. Original constant-Cp calorics
// give a separable ODE: Tg'=a*sqrt(Tg)*(Teq-Tg), Ts=Ts0-Cg/Cs*(Tg-Tg0).
// Evaluate that exact trajectory and independently integrate the Maxwell load.
fn exact_reflection_impulse(
    solid: SolidSurfaceState,
    cp: f64,
    gas: &ClosedGasParcel,
    alpha: f64,
    duration: f64,
) -> f64 {
    let (ts0, tg0) = (400.0, 300.0);
    let r = R / 0.02897;
    let cs = solid.mass_kg * cp;
    let cg = gas.mass_kg() * 2.5 * r;
    let equilibrium = (cs * ts0 + cg * tg0) / (cs + cg);
    let k = alpha * solid.area_m2 * gas.mass_kg() / gas.volume_m3()
        * (r / (2.0 * core::f64::consts::PI)).sqrt()
        * 3.0
        * r;
    let x = (tg0 / equilibrium).sqrt();
    let z0 = 0.5 * ((1.0 + x) / (1.0 - x)).ln();
    let rate = k * (1.0 + cg / cs) * equilibrium.sqrt() / (2.0 * cg);
    let force = |t: f64| {
        let v = (z0 + rate * t).tanh();
        let tg = equilibrium * v * v;
        let ts = ts0 - (cg / cs) * (tg - tg0);
        solid.area_m2 * gas.mass_kg() * r / gas.volume_m3()
            * ((1.0 - 0.5 * alpha) * tg + 0.5 * alpha * (ts * tg).sqrt())
    };
    let n = 2048;
    let dt = duration / n as f64;
    let mut sum = force(0.0) + force(duration);
    for i in 1..n {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * force(i as f64 * dt);
    }
    sum * dt / 3.0
}

#[test]
fn impulse_error_alone_resolves_changing_reflection_and_rejected_probes_add_no_load() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let solid = fixture::solid(400.0, chart);
    let rows = [row(7, &phase, &binding, solid, &temperature)];
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let scattering = BTreeMap::from([(7, maxwell(&catalog, 0.5))]);
    let mut c = controls();
    c.chamber.absolute_mass_tolerance_kg = 1.0;
    c.chamber.absolute_temperature_tolerance_k = 1e6;
    c.absolute_angular_impulse_tolerance_n_m_s = 1e6;
    c.absolute_impulse_tolerance_n_s = 1e-2;
    let loose = advance_chamber_adaptive_with_momentum(
        &rows,
        &[frame(7)],
        [0.0; 3],
        Some(&scattering),
        &gas,
        &catalog,
        0.1,
        c,
        || false,
    )
    .unwrap();
    c.absolute_impulse_tolerance_n_s = 1e-6;
    let tight = advance_chamber_adaptive_with_momentum(
        &rows,
        &[frame(7)],
        [0.0; 3],
        Some(&scattering),
        &gas,
        &catalog,
        0.1,
        c,
        || false,
    )
    .unwrap();
    let exact = exact_reflection_impulse(solid, chart.1, &gas, 0.5, 0.1);
    let impulse = |x: &MomentumChamberAdvanceCandidate<'_>| {
        x.momentum().surfaces[&7].total_normal_impulse_n_s
    };
    let loose_error = (impulse(&loose) - exact).abs();
    let tight_error = (impulse(&tight) - exact).abs();
    assert!(
        tight_error < 0.2 * loose_error,
        "{tight_error} / {loose_error}"
    );
    assert!(tight_error < 5e-4);
    assert!(
        tight.history().audit().accepted_intervals > loose.history().audit().accepted_intervals
    );
    assert!(tight.history().audit().rejected_trials > 0);
    assert!(tight.maximum_accepted_impulse_error_ratio() <= 1.0);
    assert!(tight.history().audit().maximum_accepted_error_ratio <= 1.0);
    assert_eq!(tight.gas().mass_kg(), gas.mass_kg());
    assert_eq!(
        tight.history().endpoint().surfaces()[0].solid.mass_kg,
        solid.mass_kg
    );
    assert_eq!(tight.history().audit().surfaces[0].transferred_mass_kg, 0.0);
    let load = tight.momentum().surfaces[&7];
    assert_eq!(load.outgoing_normal_impulse_n_s, 0.0);
    assert_eq!(load.sticking_normal_impulse_n_s, 0.0);
    fixture::near(
        load.angular_impulse_on_surface_n_m_s[1],
        2.0 * impulse(&tight),
        2e-14,
    );
    fixture::near(
        tight.history().audit().energy_residual_j,
        0.0,
        c.chamber.step.surface.gas.energy_tolerance_j,
    );
}

#[test]
fn equilibrium_skips_gross_inventory_retries_but_integrates_pressure_and_preserves_work_limits() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let saved = gas.clone();
    let rows = [row(
        1,
        &phase,
        &binding,
        fixture::solid(400.0, chart),
        &temperature,
    )];
    let mut c = controls();
    c.chamber.step.surface.maximum_gross_inventory_fraction = 0.001;
    let result = advance_chamber_adaptive_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.03,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.history().audit().rejected_trials, 0);
    assert_eq!(result.history().audit().accepted_intervals, 1);
    let load = result.momentum().surfaces[&1];
    fixture::near(load.total_normal_impulse_n_s, 100000.0 * 1e-4 * 0.03, 2e-15);
    fixture::near(
        load.outgoing_normal_impulse_n_s,
        0.5 * 0.1 * 1000.0 * 1e-4 * 0.03,
        2e-18,
    );
    assert_eq!(
        result.history().audit().surfaces[0].transferred_mass_kg,
        0.0
    );
    assert!(result.history().audit().kernel_evaluations <= 3 * result.history().audit().trials);
    let mut limited = c;
    limited.chamber.maximum_dt_s = 0.001;
    limited.chamber.maximum_trials = 1;
    assert!(matches!(
        advance_chamber_adaptive_with_momentum(
            &rows,
            &[frame(1)],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.03,
            limited,
            || false
        ),
        Err(GasChamberAdvanceError::TrialBudget { .. })
    ));
    assert_eq!(gas, saved);
    let mut calls = 0;
    advance_chamber_adaptive_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.03,
        c,
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut checks = 0;
    assert!(matches!(
        advance_chamber_adaptive_with_momentum(
            &rows,
            &[frame(1)],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.03,
            c,
            || {
                checks += 1;
                checks == calls
            }
        ),
        Err(GasChamberAdvanceError::Surface(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
    assert_eq!(gas, saved);
    let zero = advance_chamber_adaptive_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.0,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(zero.momentum().surfaces[&1].total_normal_impulse_n_s, 0.0);
    assert_eq!(zero.history().audit().accepted_intervals, 0);
    let mut bad = frame(1);
    bad.outward_normal[2] = 2.0;
    assert!(
        advance_chamber_adaptive_with_momentum(
            &rows,
            &[bad],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.0,
            c,
            || false
        )
        .is_err()
    );
}

#[test]
fn sharp_fusion_history_preserves_original_phase_owner_and_continues_from_its_endpoint() {
    use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (hs, cp) = fixture::chart(&catalog, &binding);
    let knot = |h, t, f| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: 917.0,
    };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::ContentHash([98; 32]),
        vec![
            knot(hs - 200.0 * cp, 200.0, 0.0),
            knot(hs, 400.0, 0.0),
            knot(hs + 100000.0, 400.0, 1.0),
        ],
    )
    .unwrap();
    let caloric = SharpMeltingCaloric::new(&curve).unwrap();
    let temperature = |h| {
        curve
            .temperature_at_specific_enthalpy(h)
            .map_err(|_| SublimationError::Invalid("original fusion chart"))
    };
    let phase = fixture::phase();
    let rows = [row(
        1,
        &phase,
        &binding,
        SolidSurfaceState {
            mass_kg: 0.01,
            specific_enthalpy_j_kg: hs + 50000.0,
            area_m2: 1e-4,
        },
        &temperature,
    )];
    let sources = BTreeMap::from([(1, caloric)]);
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let mut c = controls();
    c.chamber.maximum_dt_s = 1.0 / 256.0;
    c.chamber.absolute_mass_tolerance_kg = 1.0;
    c.chamber.absolute_temperature_tolerance_k = 1e6;
    c.absolute_impulse_tolerance_n_s = 1e20;
    c.absolute_angular_impulse_tolerance_n_m_s = 1e20;
    let duration = 1.0 / 64.0;
    let original = advance_partitioned_chamber_adaptive(
        &rows,
        &sources,
        None,
        &gas,
        &catalog,
        duration,
        c.chamber,
        || false,
    )
    .unwrap();
    let result = advance_partitioned_chamber_adaptive_with_momentum(
        &rows,
        &sources,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        duration,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.history(), original.history());
    assert_eq!(result.phases(), original.phases());
    assert!(result.phases()[&1].melted_mass_kg < 0.0);
    let first = advance_partitioned_chamber_adaptive_with_momentum(
        &rows,
        &sources,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        duration / 2.0,
        c,
        || false,
    )
    .unwrap();
    let second = advance_partitioned_chamber_adaptive_with_momentum(
        &states(&rows, first.history()),
        &sources,
        &[frame(1)],
        [0.0; 3],
        None,
        first.gas(),
        &catalog,
        duration / 2.0,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.history().endpoint(), second.history().endpoint());
    fixture::near(
        result.momentum().surfaces[&1].total_normal_impulse_n_s,
        first.momentum().surfaces[&1].total_normal_impulse_n_s
            + second.momentum().surfaces[&1].total_normal_impulse_n_s,
        2e-16,
    );
    fixture::near(
        result.phases()[&1].melted_mass_kg,
        first.phases()[&1].melted_mass_kg + second.phases()[&1].melted_mass_kg,
        1e-17,
    );
}
