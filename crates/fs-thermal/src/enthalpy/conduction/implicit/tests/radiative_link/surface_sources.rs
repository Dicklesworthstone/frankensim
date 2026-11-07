//! G1/G3/G4: source emissivity, independent surface physics and native heat.

use super::*;
use crate::enthalpy::conduction::implicit::boundary::{
    ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
};

struct Surface<'a> {
    caloric: &'a dyn ThermalLaw,
    bounds: [f64; 2],
    emissivity: f64,
    slope: f64,
    convection: [f64; 2],
    refuse: Cell<bool>,
    queries: Cell<usize>,
    last_h: Cell<f64>,
}

impl<'a> Surface<'a> {
    fn new(caloric: &'a dyn ThermalLaw, emissivity: f64, slope: f64) -> Self {
        Self {
            caloric,
            bounds: [100.0, 1000.0],
            emissivity,
            slope,
            convection: [0.0; 2],
            refuse: Cell::new(false),
            queries: Cell::new(0),
            last_h: Cell::new(0.0),
        }
    }

    fn value(&self, h: f64) -> f64 {
        self.emissivity + self.slope * (h - 400.0)
    }
}

impl ImplicitGraySurfaceLaw for Surface<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.bounds
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last_h.set(h);
        if self.refuse.get() || !(self.bounds[0]..=self.bounds[1]).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(ImplicitGraySurfaceState {
            emissivity: self.value(h),
            emissivity_enthalpy_derivative_kg_j: self.slope,
            convection_w_m2_k: self.convection[0],
            convection_enthalpy_derivative_kg_w_m2_k_j: self.convection[1],
        })
    }
}

// Retain the independent physical fourth-power/bisection oracle, replacing
// only its authored emissivity values at the ORIGINAL bulk enthalpies.
fn source_power(r: &Reference<'_>, sources: &[Surface<'_>; 2], h: [f64; 2]) -> f64 {
    let reference = Reference {
        materials: r.materials,
        masses: r.masses,
        distances: r.distances,
        coatings: r.coatings,
        eps: [sources[0].value(h[0]), sources[1].value(h[1])],
        htc: r.htc,
        gap: r.gap,
        area: r.area,
    };
    reference.power(h)
}

fn source_stage(
    r: &Reference<'_>,
    sources: &[Surface<'_>; 2],
    base: [f64; 2],
    dt: f64,
) -> [f64; 2] {
    let mut lo = (r.masses[0] * (base[0] - 1000.0)).max(r.masses[1] * (100.0 - base[1]));
    let mut hi = (r.masses[0] * (base[0] - 100.0)).min(r.masses[1] * (1000.0 - base[1]));
    for _ in 0..90 {
        let energy = 0.5 * (lo + hi);
        let h = [
            base[0] - energy / r.masses[0],
            base[1] + energy / r.masses[1],
        ];
        if energy > dt * source_power(r, sources, h) {
            hi = energy;
        } else {
            lo = energy;
        }
    }
    let energy = 0.5 * (lo + hi);
    [
        base[0] - energy / r.masses[0],
        base[1] + energy / r.masses[1],
    ]
}

#[test]
fn source_emissivities_share_one_surface_balance_and_retain_latent_tangents() {
    let materials = [material(true, 2.0, 0.0), material(false, 1.0, 0.002)];
    let laws = [&materials[0], &materials[1]];
    let sources = [
        Surface::new(laws[0], 0.6, 0.0004),
        Surface::new(laws[1], 0.7, -0.0003),
    ];
    let gap = Gap::new(laws);
    for sourced_gap in [false, true] {
        let r = Reference {
            materials: laws,
            masses: [2.0, 3.0],
            distances: [0.04, 0.08],
            coatings: [0.01, 0.02],
            eps: [0.0; 2],
            htc: if sourced_gap { 0.0 } else { 6.0 },
            gap: sourced_gap.then_some(&gap),
            area: 0.7,
        };
        let mut exchange = GrayRadiativeLink::new(
            laws.map(|law| law as &dyn ImplicitConductionLaw),
            [0.0; 2],
            r.htc,
            r.coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap()
        .with_surface_laws([&sources[0], &sources[1]])
        .unwrap();
        if sourced_gap {
            exchange = exchange.with_gap_resistance(&gap).unwrap();
        }
        for h in [[400.0, 650.0], [700.0, 250.0], [800.0, 300.0]] {
            let before = sources.each_ref().map(|s| s.queries.get());
            let result = exchange
                .evaluate(h, geometry(r.distances), patch(r.area), &mut || Ok(()))
                .unwrap();
            near(result.power_w, source_power(&r, &sources, h), 2e-8);
            near(
                result.power_w,
                result.gap_power_w + result.radiation_power_w,
                2e-8,
            );
            for side in 0..2 {
                assert_eq!(sources[side].queries.get(), before[side] + 1);
                assert_eq!(sources[side].last_h.get(), h[side]);
                let (mut plus, mut minus) = (h, h);
                plus[side] += 0.001;
                minus[side] -= 0.001;
                let expected =
                    (source_power(&r, &sources, plus) - source_power(&r, &sources, minus)) / 0.002;
                near(
                    result.link.power_enthalpy_derivatives_w_kg_j[side],
                    expected,
                    3e-7,
                );
            }
            if h[0] == 400.0 && !sourced_gap {
                // This body's T and k are fixed on its latent plateau: only
                // the admitted emissivity derivative can create this tangent.
                assert_eq!(
                    materials[0]
                        .state_at_specific_enthalpy(h[0])
                        .unwrap()
                        .temperature_slope_kg_k_j,
                    0.0
                );
                assert!(result.link.power_enthalpy_derivatives_w_kg_j[0].abs() > 1e-4);
            }
        }
    }
}

#[test]
fn zero_emissivity_onset_retains_radiation_derivative_without_heat() {
    let materials = [material(true, 2.0, 0.0), material(false, 3.0, 0.0)];
    let sources = [
        Surface::new(&materials[0], 0.0, 0.001),
        Surface::new(&materials[1], 0.8, 0.0),
    ];
    let exchange = GrayRadiativeLink::new(
        [&materials[0], &materials[1]],
        [0.0; 2],
        0.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_surface_laws([&sources[0], &sources[1]])
    .unwrap();
    let area = 0.7;
    let result = exchange
        .evaluate([400.0, 650.0], geometry([0.1; 2]), patch(area), &mut || {
            Ok(())
        })
        .unwrap();
    assert_eq!(result.power_w, 0.0);
    assert_eq!(result.link.conductance_w_k, 0.0);
    assert_eq!(result.surface_temperatures_k, Some([300.0, 650.0]));
    near(
        result.link.power_enthalpy_derivatives_w_kg_j[0],
        area * SIGMA * (300.0_f64.powi(4) - 650.0_f64.powi(4)) * 0.001,
        2e-14,
    );
    assert_eq!(result.link.power_enthalpy_derivatives_w_kg_j[1], 0.0);
    let reversed = GrayRadiativeLink::new(
        [&materials[1], &materials[0]],
        [0.0; 2],
        0.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_surface_laws([&sources[1], &sources[0]])
    .unwrap()
    .evaluate([650.0, 400.0], geometry([0.1; 2]), patch(area), &mut || {
        Ok(())
    })
    .unwrap();
    assert_eq!(reversed.link.power_enthalpy_derivatives_w_kg_j[0], 0.0);
    near(
        reversed.link.power_enthalpy_derivatives_w_kg_j[1],
        -result.link.power_enthalpy_derivatives_w_kg_j[0],
        2e-14,
    );
}

#[test]
fn constant_surface_sources_recover_fixed_gray_limits_and_reversal() {
    let materials = [material(false, 2.0, 0.0), material(false, 3.0, 0.0)];
    for eps in [[0.6, 0.9], [1.0; 2], [1e-300, 1e-100], [0.0, 0.8], [0.0; 2]] {
        let sources = [
            Surface::new(&materials[0], eps[0], 0.0),
            Surface::new(&materials[1], eps[1], 0.0),
        ];
        let fixed = GrayRadiativeLink::new(
            [&materials[0], &materials[1]],
            eps,
            0.0,
            [0.01, 0.02],
            GrayRadiativeControls::default(),
        )
        .unwrap();
        let sourced = GrayRadiativeLink::new(
            [&materials[0], &materials[1]],
            [0.0; 2],
            0.0,
            [0.01, 0.02],
            GrayRadiativeControls::default(),
        )
        .unwrap()
        .with_surface_laws([&sources[0], &sources[1]])
        .unwrap();
        let h = [800.0, 300.0];
        let face = geometry([0.04, 0.08]);
        let expected = fixed.evaluate(h, face, patch(0.7), &mut || Ok(())).unwrap();
        let actual = sourced
            .evaluate(h, face, patch(0.7), &mut || Ok(()))
            .unwrap();
        assert_eq!(actual, expected);
        let reversed = GrayRadiativeLink::new(
            [&materials[1], &materials[0]],
            [0.0; 2],
            0.0,
            [0.02, 0.01],
            GrayRadiativeControls::default(),
        )
        .unwrap()
        .with_surface_laws([&sources[1], &sources[0]])
        .unwrap()
        .evaluate(
            [h[1], h[0]],
            geometry([0.08, 0.04]),
            patch(0.7),
            &mut || Ok(()),
        )
        .unwrap();
        near(
            reversed.power_w,
            -actual.power_w,
            actual.power_w.abs() * 1e-13,
        );
    }
}

#[test]
fn source_emissivity_drives_native_be_and_sdirk_on_unequal_finite_stores() {
    let materials = [material(false, 0.8, 0.002), material(true, 0.5, 0.002)];
    let laws = [&materials[0], &materials[1]];
    let sources = [
        Surface::new(laws[0], 0.6, 0.0004),
        Surface::new(laws[1], 0.7, -0.0003),
    ];
    let gap = Gap::new(laws);
    let r = Reference {
        materials: laws,
        masses: [2.0, 3.0],
        distances: [0.04, 0.08],
        coatings: [0.01, 0.02],
        eps: [0.0; 2],
        htc: 0.0,
        gap: Some(&gap),
        area: 0.7,
    };
    let initial = [800.0, 400.0];
    let dt = 0.05;
    let (thermal, observer) = pair(laws, initial, r.masses, geometry(r.distances), r.area);
    let exchange = GrayRadiativeLink::new(
        laws.map(|law| law as &dyn ImplicitConductionLaw),
        [0.0; 2],
        0.0,
        r.coatings,
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_surface_laws([&sources[0], &sources[1]])
    .unwrap()
    .with_gap_resistance(&gap)
    .unwrap();
    let bindings: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
        [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
    let chain = ImplicitConductionChain::new(&thermal, &observer, &bindings, || Ok(()))
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &bindings,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
    .unwrap();
    let generation = |h| {
        let temperatures = r.temperatures(h);
        source_power(&r, &sources, h) * (1.0 / temperatures[1] - 1.0 / temperatures[0])
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    for sdirk in [false, true] {
        let (expected, expected_generation) = if sdirk {
            let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
            let first = source_stage(&r, &sources, initial, gamma * dt);
            let earlier = (1.0 - gamma) * dt * source_power(&r, &sources, first);
            let base = [
                initial[0] - earlier / r.masses[0],
                initial[1] + earlier / r.masses[1],
            ];
            let expected = source_stage(&r, &sources, base, gamma * dt);
            (
                expected,
                dt * ((1.0 - gamma) * generation(first) + gamma * generation(expected)),
            )
        } else {
            let expected = source_stage(&r, &sources, initial, dt);
            (expected, dt * generation(expected))
        };
        let results = if sdirk {
            [
                chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            ]
        } else {
            [
                chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            ]
        };
        for result in results {
            let (next, report) = result.unwrap();
            for side in 0..2 {
                near(
                    next.cells()[&CellId(side as u64)].specific_h_j_kg,
                    expected[side],
                    3e-9,
                );
            }
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.residual_j, 0.0, 1e-9);
            near(report.entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
            near(
                report.entropy.conduction_generation_j_k,
                expected_generation,
                2e-10,
            );
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert_eq!(next.temperature_k(CellId(1)).unwrap(), 300.0);
        }
    }
    // A late source refusal cannot publish heat or alter the accepted cells.
    sources[1].refuse.set(true);
    assert!(
        chain
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial[0]);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial[1]);
}

#[test]
fn surface_sources_refuse_wrong_identity_duplicate_heat_and_unsupported_states() {
    let materials = [material(false, 0.0, 0.0), material(false, 3.0, 0.0)];
    let laws = [&materials[0], &materials[1]];
    let fresh = |eps| {
        GrayRadiativeLink::new(
            laws.map(|law| law as &dyn ImplicitConductionLaw),
            eps,
            0.0,
            [0.0; 2],
            GrayRadiativeControls::default(),
        )
        .unwrap()
    };
    let sources = [
        Surface::new(laws[0], 0.6, 0.0),
        Surface::new(laws[1], 0.7, 0.0),
    ];
    assert!(
        fresh([0.0; 2])
            .with_surface_laws([&sources[1], &sources[0]])
            .is_err()
    );
    assert!(
        fresh([0.0, 0.7])
            .with_surface_laws([&sources[0], &sources[1]])
            .is_err()
    );
    assert!(
        fresh([0.0; 2])
            .with_surface_laws([&sources[0], &sources[1]])
            .unwrap()
            .with_surface_laws([&sources[0], &sources[1]])
            .is_err()
    );
    for convection in [[1.0, 0.0], [0.0, 0.001]] {
        let mut bad = Surface::new(laws[0], 0.5, 0.0);
        bad.convection = convection;
        let link = fresh([0.0; 2])
            .with_surface_laws([&bad, &sources[1]])
            .unwrap();
        assert_eq!(
            link.evaluate([400.0, 650.0], geometry([0.1; 2]), patch(1.0), &mut || Ok(
                ()
            )),
            Err(DepositionError::MaterialDomain)
        );
    }
    let mut limited = Surface::new(laws[0], 0.5, 0.0);
    limited.bounds = [410.0, 600.0];
    let link = fresh([0.0; 2])
        .with_surface_laws([&limited, &sources[1]])
        .unwrap();
    // Even an insulating material may not bypass the source domain.
    assert_eq!(
        link.evaluate([400.0, 650.0], geometry([0.1; 2]), patch(1.0), &mut || Ok(
            ()
        )),
        Err(DepositionError::MaterialDomain)
    );
    let link = fresh([0.0; 2])
        .with_surface_laws([&sources[0], &sources[1]])
        .unwrap();
    assert_eq!(
        link.evaluate([400.0, 650.0], geometry([0.1; 2]), patch(1.0), &mut || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    );
}

#[test]
fn simultaneous_zero_emissivity_onsets_refuse_a_nonexistent_joint_tangent() {
    let laws = [material(true, 3.0, 0.0), material(false, 3.0, 0.0)];
    for slopes in [[0.001, 0.002], [0.001, 0.0], [0.0, 0.0]] {
        let sources = [
            Surface::new(&laws[0], 0.0, slopes[0]),
            Surface::new(&laws[1], 0.0, slopes[1]),
        ];
        let link = GrayRadiativeLink::new(
            [&laws[0], &laws[1]],
            [0.0; 2],
            0.0,
            [0.0; 2],
            GrayRadiativeControls::default(),
        )
        .unwrap()
        .with_surface_laws([&sources[0], &sources[1]])
        .unwrap();
        let response = link.evaluate([400.0; 2], geometry([0.1; 2]), patch(1.0), &mut || Ok(()));
        if slopes[1] != 0.0 {
            assert!(matches!(response, Err(DepositionError::InvalidInput(_))));
        } else {
            let response = response.unwrap();
            assert_eq!(response.power_w, 0.0);
            assert_eq!(response.link.power_enthalpy_derivatives_w_kg_j, [0.0; 2]);
        }
    }
}
