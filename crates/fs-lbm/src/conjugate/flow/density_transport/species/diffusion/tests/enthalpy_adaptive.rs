//! Original-path nonlinear caloric feedback, time refinement and physical gates.
use super::super::enthalpy::adaptive::{
    AdaptiveSpeciesEnthalpyDiffusionControls, SpeciesEnthalpyDiffusionModel,
};
use super::*;
use std::cell::Cell;

struct CaloricModel {
    diffusivity: f64,
    temperature_dependent_diffusivity: bool,
    closure_defect: f64,
    minimum_temperature: f64,
    cancel_at_preparation: usize,
    preparations: Cell<usize>,
}

struct PreparedCaloric {
    temperature: Vec<f64>,
    diffusivity: Vec<f64>,
}

impl CaloricModel {
    fn new() -> Self {
        Self {
            diffusivity: 0.2,
            temperature_dependent_diffusivity: false,
            closure_defect: 0.0,
            minimum_temperature: 1.0,
            cancel_at_preparation: usize::MAX,
            preparations: Cell::new(0),
        }
    }
}

impl SpeciesEnthalpyDiffusionModel<2> for CaloricModel {
    type Prepared = PreparedCaloric;
    fn prepare(
        &self,
        domain: &VoxelDomain,
        state: &MaterialAdvectionStep,
        species: &[[f64; 2]],
        h: &[f64],
        gate: &CancelGate,
    ) -> Result<Self::Prepared, ChtError> {
        let calls = self.preparations.get() + 1;
        self.preparations.set(calls);
        if calls == self.cancel_at_preparation {
            gate.request();
        }
        poll(gate)?;
        let mut temperature = Vec::new();
        let mut diffusivity = Vec::new();
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                temperature.push(0.0);
                diffusivity.push(0.0);
                continue;
            }
            // Original ideal-mixture law h=(1+Y_A)*T^2/2, with positive
            // temperature-dependent Cp_A=2T and Cp_B=T [J/(kg K)].
            let t = (2.0 * h[c] / (1.0 + species[c][0] / state.mass_kg[c])).sqrt();
            if !t.is_finite() || t < self.minimum_temperature {
                return Err(refuse("original quadratic caloric source is out of domain"));
            }
            temperature.push(t);
            diffusivity.push(
                self.diffusivity
                    * if self.temperature_dependent_diffusivity {
                        (t / 300.0).powi(2)
                    } else {
                        1.0
                    },
            );
        }
        Ok(PreparedCaloric {
            temperature,
            diffusivity,
        })
    }
    fn temperature_k<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
        &p.temperature
    }
    fn diffusivity_m2_s<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
        &p.diffusivity
    }
    fn face_partial_enthalpy_j_kg(
        &self,
        p: &Self::Prepared,
        _: &VoxelDomain,
        cells: [usize; 2],
        _: usize,
        gate: &CancelGate,
    ) -> Result<[f64; 2], ChtError> {
        poll(gate)?;
        // Declared centered interpolation of each original species enthalpy.
        let sensible = 0.25 * (p.temperature[cells[0]].powi(2) + p.temperature[cells[1]].powi(2));
        Ok([2.0 * sensible, sensible])
    }
    fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
        cells
            .checked_mul(16)?
            .checked_add(core::mem::size_of::<PreparedCaloric>())
    }
    fn caloric_defect_j(&self, _: &Self::Prepared) -> Result<f64, ChtError> {
        Ok(self.closure_defect)
    }
}

fn adaptive_controls() -> AdaptiveSpeciesEnthalpyDiffusionControls<2> {
    AdaptiveSpeciesEnthalpyDiffusionControls {
        solver: controls(1),
        initial_substeps: 1,
        maximum_refinements: 12,
        maximum_source_evaluations: 100_000,
        absolute_temperature_tolerance_k: 1e-4,
        absolute_specific_enthalpy_tolerance_j_kg: 0.01,
        absolute_mass_fraction_tolerance: [1e-10; 2],
        relative_mass_fraction_tolerance: 0.0,
        maximum_enthalpy_defect_j: 1e-5,
    }
}

fn hot_cold() -> ([[f64; 2]; 2], [f64; 2]) {
    (
        [[0.8, 0.2], [0.2, 0.8]],
        [1.8 * 0.5 * 300.0_f64.powi(2), 1.2 * 0.5 * 600.0_f64.powi(2)],
    )
}

// Independent continuous two-cell solution. delta=(Y_A,left-Y_A,right)/2
// decays at 2D. Writing a=1+mean(Y_A), E=h_l+h_r, the energy equation gives
// h_l=E/2 + E*delta/(2a) + C*sqrt(a^2-delta^2). This integrates the authored
// face law analytically; it does not reuse the discrete midpoint method.
fn exact_hot_cold(duration: f64) -> ([f64; 2], [f64; 2]) {
    let (_, h) = hot_cold();
    let a = 1.5_f64;
    let delta0 = 0.3_f64;
    let energy = h[0] + h[1];
    let constant =
        (h[0] - 0.5 * energy - energy * delta0 / (2.0 * a)) / (a * a - delta0 * delta0).sqrt();
    let delta = delta0 * (-0.4 * duration).exp();
    let left =
        0.5 * energy + energy * delta / (2.0 * a) + constant * (a * a - delta * delta).sqrt();
    let enthalpy = [left, energy - left];
    let temperature = [
        (2.0 * enthalpy[0] / (a + delta)).sqrt(),
        (2.0 * enthalpy[1] / (a - delta)).sqrt(),
    ];
    (enthalpy, temperature)
}

#[test]
fn adaptive_caloric_isothermal_limit_matches_original_frozen_transport() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let species = [[0.8, 0.2], [0.2, 0.8]];
    let h = [1.8 * 45_000.0, 1.2 * 45_000.0];
    let gate = CancelGate::new_clock_free();
    let model = CaloricModel::new();
    let got = state
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            2.0,
            &model,
            adaptive_controls(),
            &gate,
        )
        .unwrap();
    let reference = state
        .diffuse_species_with_enthalpy(
            &domain,
            &species,
            &h,
            &[0.2; 2],
            2.0,
            |_, _| Ok([90_000.0, 45_000.0]),
            controls(1),
            1e-7,
            &gate,
        )
        .unwrap();
    assert_eq!(got.report.accepted_substeps, 2);
    assert_eq!(got.report.total_pair_updates, 12);
    assert_eq!(got.report.total_component_updates, 24);
    assert_eq!(got.report.source_preparations, 7);
    // One extra original coefficient query admits the stiff coarse step.
    assert_eq!(got.report.face_source_evaluations, 7);
    assert_eq!(got.state.temperature, got.temperature_k);
    assert_eq!(got.state.diffusivity, [0.2; 2]);
    for c in 0..2 {
        near(got.temperature_k[c], 300.0, 1e-10);
        near(
            got.specific_enthalpy_j_kg[c],
            reference.specific_enthalpy_j_kg[c],
            1e-8,
        );
        for k in 0..2 {
            near(
                got.species_mass_kg[c][k],
                reference.species.species_mass_kg[c][k],
                1e-14,
            );
        }
    }
    assert_eq!(state.mass_kg, [1.0; 2]);
}

#[test]
fn adaptive_interval_preserves_native_heat_tails_under_large_formation_cancellation() {
    struct FormationModel;
    impl SpeciesEnthalpyDiffusionModel<2> for FormationModel {
        type Prepared = PreparedCaloric;
        fn prepare(
            &self,
            domain: &VoxelDomain,
            carrier: &MaterialAdvectionStep,
            species: &[[f64; 2]],
            h: &[f64],
            gate: &CancelGate,
        ) -> Result<Self::Prepared, ChtError> {
            let formation = 2.0_f64.powi(54);
            let mut temperature = Vec::new();
            for c in 0..domain.cell_count() {
                poll_every(gate, c)?;
                let (a, ae) =
                    fs_math::eft::two_prod(species[c][0] / carrier.mass_kg[c], -formation);
                let (b, be) =
                    fs_math::eft::two_prod(species[c][1] / carrier.mass_kg[c], formation + 4.0);
                let (sum, error) = fs_math::eft::two_sum(a, b);
                let formation_h = sum + (error + ae + be);
                // Original constant-Cp mixture h=h_formation+(T-300).
                temperature.push(300.0 + (h[c] - formation_h));
            }
            Ok(PreparedCaloric {
                temperature,
                diffusivity: vec![1e30; h.len()],
            })
        }
        fn temperature_k<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
            &p.temperature
        }
        fn diffusivity_m2_s<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
            &p.diffusivity
        }
        fn face_partial_enthalpy_j_kg(
            &self,
            _: &Self::Prepared,
            _: &VoxelDomain,
            _: [usize; 2],
            _: usize,
            gate: &CancelGate,
        ) -> Result<[f64; 2], ChtError> {
            poll(gate)?;
            let formation = 2.0_f64.powi(54);
            Ok([-formation, formation + 4.0])
        }
        fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
            cells
                .checked_mul(16)?
                .checked_add(core::mem::size_of::<PreparedCaloric>())
        }
        fn caloric_defect_j(&self, _: &Self::Prepared) -> Result<f64, ChtError> {
            Ok(0.0)
        }
    }
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let carrier = carrier(&domain, &[1.0; 2]);
    let formation = 2.0_f64.powi(54);
    let mut controls = adaptive_controls();
    controls.maximum_enthalpy_defect_j = 0.0;
    let result = carrier
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &[[1.0, 0.0], [0.0, 1.0]],
            &[-formation, formation + 4.0],
            1.0,
            &FormationModel,
            controls,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // The native transfer is (F,2), while its scalar public readback rounds
    // to F. Reusing that readback would invent a four-joule interval defect.
    assert_eq!(result.specific_enthalpy_j_kg, [2.0; 2]);
    assert_eq!(result.species_mass_kg, [[0.5; 2]; 2]);
    assert_eq!(result.temperature_k, [300.0; 2]);
    assert_eq!(result.state.temperature, result.temperature_k);
    assert_eq!(result.diffusive_enthalpy_change_j, [formation, -formation]);
    assert_eq!(result.report.absolute_enthalpy_defect_j, 0.0);
    assert_eq!(result.report.maximum_changing_source_relaxation, 0.0);
}

#[test]
fn changing_source_stiff_equilibration_cannot_alias_the_thermal_error_estimator() {
    struct Cubic {
        diffusivity: f64,
    }
    struct PreparedCubic {
        temperature: Vec<f64>,
        diffusivity: Vec<f64>,
        y: Vec<f64>,
    }
    impl SpeciesEnthalpyDiffusionModel<2> for Cubic {
        type Prepared = PreparedCubic;
        fn prepare(
            &self,
            domain: &VoxelDomain,
            carrier: &MaterialAdvectionStep,
            species: &[[f64; 2]],
            h: &[f64],
            gate: &CancelGate,
        ) -> Result<Self::Prepared, ChtError> {
            let mut y = Vec::new();
            let mut temperature = Vec::new();
            for c in 0..domain.cell_count() {
                poll_every(gate, c)?;
                let fraction = species[c][0] / carrier.mass_kg[c];
                y.push(fraction);
                temperature.push(300.0 + h[c] - 4.0 * fraction.powi(3));
            }
            Ok(PreparedCubic {
                temperature,
                diffusivity: vec![self.diffusivity; h.len()],
                y,
            })
        }
        fn temperature_k<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
            &p.temperature
        }
        fn diffusivity_m2_s<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
            &p.diffusivity
        }
        fn face_partial_enthalpy_j_kg(
            &self,
            p: &Self::Prepared,
            _: &VoxelDomain,
            cells: [usize; 2],
            _: usize,
            gate: &CancelGate,
        ) -> Result<[f64; 2], ChtError> {
            poll(gate)?;
            let mut face = [0.0; 2];
            for c in cells {
                let y = p.y[c];
                let sensible = p.temperature[c] - 300.0;
                face[0] += 0.5 * (sensible + 4.0 * (3.0 * y * y - 2.0 * y.powi(3)));
                face[1] += 0.5 * (sensible - 8.0 * y.powi(3));
            }
            Ok(face)
        }
        fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
            cells
                .checked_mul(24)?
                .checked_add(core::mem::size_of::<PreparedCubic>())
        }
        fn caloric_defect_j(&self, _: &Self::Prepared) -> Result<f64, ChtError> {
            Ok(0.0)
        }
    }
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let carrier = carrier(&domain, &[1.0; 2]);
    let species = [[1.0, 0.0], [0.0, 1.0]];
    let h = [4.0, 0.0];
    let mut controls = adaptive_controls();
    controls.maximum_refinements = 3;
    // For D=1000, both one and two midpoint steps fully equilibrate species
    // while aliasing the SAME wrong h=[2.5,1.5], T=[302,301]. Their difference
    // estimate is zero despite a real half-kelvin error. Inner frozen sweeps
    // do not resolve the missing outer property history either.
    for inner in [1, 4096] {
        controls.solver.substeps = inner;
        assert!(
            carrier
                .diffuse_species_with_enthalpy_adaptive(
                    &domain,
                    &species,
                    &h,
                    1.0,
                    &Cubic {
                        diffusivity: 1000.0
                    },
                    controls,
                    &CancelGate::new_clock_free()
                )
                .is_err()
        );
    }
    controls = adaptive_controls();
    controls.absolute_specific_enthalpy_tolerance_j_kg = 1e-4;
    let d = 8.0;
    let result = carrier
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            1.0,
            &Cubic { diffusivity: d },
            controls,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let v = 0.5 * (-2.0_f64 * d).exp();
    let right_h = 2.0 - 3.0 * v - 4.0 * v.powi(3);
    let exact_temperature = 301.5 - 6.0 * v * v;
    near(result.specific_enthalpy_j_kg[0], 4.0 - right_h, 1.1e-4);
    near(result.specific_enthalpy_j_kg[1], right_h, 1.1e-4);
    for &temperature in &result.temperature_k {
        near(temperature, exact_temperature, 1.1e-4);
    }
    assert!(result.report.maximum_changing_source_relaxation <= 0.5);
    assert!(result.report.refinements > 1);
}

#[test]
fn midpoint_caloric_feedback_has_second_order_thermal_accuracy() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let (species, h) = hot_cold();
    let (exact_h, exact_t) = exact_hot_cold(2.0);
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut ctl = adaptive_controls();
        ctl.initial_substeps = count;
        ctl.maximum_refinements = 1;
        ctl.absolute_temperature_tolerance_k = 1e6;
        ctl.absolute_specific_enthalpy_tolerance_j_kg = 1e9;
        let result = state
            .diffuse_species_with_enthalpy_adaptive(
                &domain,
                &species,
                &h,
                2.0,
                &CaloricModel::new(),
                ctl,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.report.accepted_substeps, 2 * count);
        let error = (result.temperature_k[0] - exact_t[0]).abs();
        assert!(error > 1e-8);
        errors.push(error);
        near(
            result.specific_enthalpy_j_kg.iter().sum(),
            h.iter().sum(),
            1e-7,
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.7 && pair[0] / pair[1] < 4.3,
            "{errors:?}"
        );
    }
    let resolved = state
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            2.0,
            &CaloricModel::new(),
            adaptive_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(resolved.report.refinements > 1);
    assert!(resolved.report.maximum_error_ratio <= 1.0);
    for c in 0..2 {
        near(resolved.temperature_k[c], exact_t[c], 1.1e-4);
        near(resolved.specific_enthalpy_j_kg[c], exact_h[c], 0.011);
    }
    assert!(
        resolved.report.absolute_enthalpy_defect_j < adaptive_controls().maximum_enthalpy_defect_j
    );
}

#[test]
fn adaptive_caloric_closure_is_separate_and_does_not_count_rejected_histories() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let (species, h) = hot_cold();
    let mut model = CaloricModel::new();
    model.closure_defect = 1e-6;
    let mut ctl = adaptive_controls();
    ctl.maximum_enthalpy_defect_j = 3e-6;
    let result = state
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            2.0,
            &model,
            ctl,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(result.report.refinements > 1);
    assert_eq!(result.report.initial_caloric_defect_j, 1e-6);
    assert_eq!(result.report.endpoint_caloric_defect_j, 1e-6);
    assert_eq!(result.report.maximum_stage_caloric_defect_j, 1e-6);
    assert!(result.report.absolute_enthalpy_defect_j >= 2e-6);
    assert!(result.report.absolute_enthalpy_defect_j < 3e-6);
    ctl.maximum_enthalpy_defect_j = 1.5e-6;
    assert!(
        state
            .diffuse_species_with_enthalpy_adaptive(
                &domain,
                &species,
                &h,
                2.0,
                &model,
                ctl,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
}

#[test]
fn midpoint_refreshes_temperature_dependent_diffusivity_and_partial_enthalpy_together() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let (species, h) = hot_cold();
    let energy = h.iter().sum::<f64>();
    // Independent RK4 integration of the two continuous scalar equations.
    // Their half-cell resistance is G=2D_lD_r/(D_l+D_r), D=.2(T/300)^2.
    let derivative = |state: [f64; 2]| {
        let [y, enthalpy] = state;
        let sensible = [enthalpy / (1.0 + y), (energy - enthalpy) / (2.0 - y)];
        let d = sensible.map(|e| 0.2 * 2.0 * e / 90_000.0);
        let mass_rate = (2.0 * d[0] * d[1] / (d[0] + d[1])) * (1.0 - 2.0 * y);
        [mass_rate, 0.5 * (sensible[0] + sensible[1]) * mass_rate]
    };
    let reference = |steps: usize| {
        let dt = 2.0 / steps as f64;
        let mut value = [0.8, h[0]];
        for _ in 0..steps {
            let a = derivative(value);
            let b = derivative(core::array::from_fn(|i| value[i] + 0.5 * dt * a[i]));
            let c = derivative(core::array::from_fn(|i| value[i] + 0.5 * dt * b[i]));
            let d = derivative(core::array::from_fn(|i| value[i] + dt * c[i]));
            value = core::array::from_fn(|i| {
                value[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i])
            });
        }
        value
    };
    let exact = reference(2048);
    let refined_reference = reference(4096);
    near(exact[0], refined_reference[0], 1e-12);
    near(exact[1], refined_reference[1], 1e-6);
    let mut model = CaloricModel::new();
    model.temperature_dependent_diffusivity = true;
    let mut ctl = adaptive_controls();
    ctl.absolute_mass_fraction_tolerance = [1e-6; 2];
    let result = state
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            2.0,
            &model,
            ctl,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(result.report.refinements > 1);
    near(result.species_mass_kg[0][0], exact[0], 1.1e-6);
    near(result.specific_enthalpy_j_kg[0], exact[1], 0.011);
    let fixed_d_fraction = 0.5 + 0.3 * (-0.8_f64).exp();
    assert!((result.species_mass_kg[0][0] - fixed_d_fraction).abs() > 0.01);
}

#[test]
fn adaptive_caloric_resource_source_and_accuracy_refusals_publish_no_state() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let (species, h) = hot_cold();
    let saved = (species, h, state.mass_kg.clone());
    let run = |model: &CaloricModel, ctl| {
        state.diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            2.0,
            model,
            ctl,
            &CancelGate::new_clock_free(),
        )
    };
    let mut ctl = adaptive_controls();
    ctl.solver.maximum_pair_updates = 11;
    assert!(run(&CaloricModel::new(), ctl).is_err());
    ctl = adaptive_controls();
    ctl.maximum_source_evaluations = 12;
    assert!(run(&CaloricModel::new(), ctl).is_err());
    ctl = adaptive_controls();
    ctl.solver.maximum_storage_bytes = 1;
    let source = CaloricModel::new();
    assert!(run(&source, ctl).is_err());
    assert_eq!(source.preparations.get(), 0);
    ctl = adaptive_controls();
    ctl.maximum_refinements = 1;
    ctl.absolute_temperature_tolerance_k = 1e-15;
    assert!(run(&CaloricModel::new(), ctl).is_err());
    let mut source = CaloricModel::new();
    source.minimum_temperature = 299.0;
    assert!(run(&source, adaptive_controls()).is_err());
    assert_eq!(source.preparations.get(), 2); // Midpoint leaves the original chart.
    let mut source = CaloricModel::new();
    source.cancel_at_preparation = 3;
    assert!(run(&source, adaptive_controls()).is_err());
    assert_eq!(source.preparations.get(), 3);
    assert_eq!((species, h, state.mass_kg.clone()), saved);
}
