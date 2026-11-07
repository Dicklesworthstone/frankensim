//! G0/G1/G3/G4: physical continuation through original finite ductile owners.
use super::*;
use fs_material::ductile::{DuctileError, JohnsonCookInitiation, JohnsonCookParameters};
use fs_material::plastic::finite::ductile::FiniteDuctileError;
use fs_solid::hyper3::plastic::ductile::dynamics::FiniteDuctileTetMotion;
use fs_solid::hyper3::plastic::ductile::path::{
    FiniteDuctileMotionPathControls, FiniteDuctilePathControls, FiniteDuctilePathError,
    FiniteDuctileRefinementCause,
};

const VOLUME: f64 = 1.0 / 6.0;
const RHO: f64 = 12.0;

fn path_controls() -> FiniteDuctilePathControls {
    FiniteDuctilePathControls {
        initial_fraction: 1.0,
        minimum_fraction: 1.0 / 65536.0,
        maximum_fraction: 1.0,
        maximum_accepted_steps: 512,
        maximum_trials: 1024,
        maximum_damage_increment: 0.025,
        initial_equilibrium_tolerance_n: 1e-8,
        work: FinitePlasticStepControls {
            maximum_endpoint_work_defect_j: 1e-3,
            ..controls()
        },
    }
}

#[test]
fn adaptive_static_softening_matches_independent_stress_work_and_exact_substep_replay() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let width = bands[0].band_width_m();
        let material = model(width, 0.04);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let target = stretch(0.08);
        let supports = affine(&NODES, target, [0.0; 3], 0);
        // The direct endpoint fails with its original typed material limit.
        assert!(matches!(
            solid.advance(
                &supports,
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            ),
            Err(FiniteDuctileTetError::Material {
                source: FiniteDuctileError::DamageIncrement { .. },
                ..
            })
        ));
        let controls = path_controls();
        let path = solid
            .advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                controls,
                cx,
            )
            .unwrap();
        assert!(path.rejected_trials > 0 && path.records.len() > 1);
        assert_eq!(path.trials, path.records.len() + path.rejected_trials);
        assert_eq!(path.records.last().unwrap().fraction, 1.0);
        near(
            path.records.iter().map(|r| r.duration_s).sum(),
            conditions().duration_s,
        );
        assert_eq!(path.energy.inertia, None);

        let mut alpha_old = 0.0;
        let mut g_old = 1.0;
        let mut native_old = 0.0;
        let mut stored_old = 0.0;
        let mut f_old = [1.0; 3];
        let mut fraction_old = 0.0;
        let mut work_sum = 0.0;
        let mut defect_sum = 0.0;
        let mut absolute_sum = 0.0;
        let mut plastic_sum = 0.0;
        let mut damage_sum = 0.0;
        let mut replay = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        for record in &path.records {
            let fraction = record.fraction;
            let span = fraction - fraction_old;
            near(record.duration_s, conditions().duration_s * span);
            let f = [
                1.0 + fraction * (target[0] - 1.0),
                1.0 + fraction * (target[4] - 1.0),
                1.0 + fraction * (target[8] - 1.0),
            ];
            // Axisymmetric coaxial log-J2 has an independent scalar solution.
            let logarithmic_ratio = (f[0] / f[1]).ln();
            let alpha = (2.0 / 3.0 * logarithmic_ratio - 0.002).max(0.0);
            let q = (1000.0 * logarithmic_ratio).min(3.0);
            let log_j = f.iter().map(|v| v.ln()).sum::<f64>();
            let g = (-width * 3.0 * (alpha - ONSET).max(0.0) / GP).exp();
            let native = q * q / 3000.0 + 500.0 * log_j.max(0.0).powi(2);
            let stored = g * native + 500.0 * log_j.min(0.0).powi(2);
            let volume_stress = 1000.0 * log_j * if log_j < 0.0 { 1.0 } else { g };
            let p = [
                (volume_stress + 2.0 / 3.0 * g * q) / f[0],
                (volume_stress - g * q / 3.0) / f[1],
                (volume_stress - g * q / 3.0) / f[2],
            ];
            let work = (0..3).map(|i| p[i] * (f[i] - f_old[i])).sum::<f64>();
            let plastic = g * 3.0 * (alpha - alpha_old);
            let damage = (g_old - g) * native_old;
            let defect = work - (stored - stored_old) - plastic - damage;
            let energy = record.energy.material;
            near(energy.stored_after_j, VOLUME * stored);
            near(energy.plastic_dissipation_j, VOLUME * plastic);
            near(energy.damage_release_j, VOLUME * damage);
            near(energy.internal_endpoint_work_j, VOLUME * work);
            near(energy.support_endpoint_work_j, VOLUME * work);
            near(energy.endpoint_work_defect_j, VOLUME * defect);
            near(
                energy.absolute_endpoint_work_defect_j,
                VOLUME * defect.abs(),
            );
            near(record.maximum_damage, 1.0 - g);
            assert!(record.maximum_damage_increment <= controls.maximum_damage_increment);
            assert!(
                energy.absolute_endpoint_work_defect_j
                    <= controls.work.maximum_endpoint_work_defect_j * span
            );
            work_sum += energy.internal_endpoint_work_j;
            defect_sum += energy.endpoint_work_defect_j;
            absolute_sum += energy.absolute_endpoint_work_defect_j;
            plastic_sum += energy.plastic_dissipation_j;
            damage_sum += energy.damage_release_j;

            // Replay precisely the accepted physical partition through the
            // existing endpoint owner: retries must not consume history/time.
            let current_supports: Vec<_> = supports
                .iter()
                .map(|&(dof, value)| (dof, fraction * value))
                .collect();
            let endpoint = replay
                .advance(
                    &current_supports,
                    &[],
                    [0.0; 3],
                    DuctileStepConditions {
                        duration_s: record.duration_s,
                        ..conditions()
                    },
                    IncrementSettings::default(),
                    FinitePlasticStepControls {
                        maximum_endpoint_work_defect_j: controls
                            .work
                            .maximum_endpoint_work_defect_j
                            * span,
                        equilibrium_work_tolerance_j: controls.work.equilibrium_work_tolerance_j
                            * span,
                        ..controls.work
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(endpoint.energy, energy);
            replay = endpoint.solid;
            f_old = f;
            alpha_old = alpha;
            g_old = g;
            native_old = native;
            stored_old = stored;
            fraction_old = fraction;
        }
        near(path.energy.material.internal_endpoint_work_j, work_sum);
        near(path.energy.material.endpoint_work_defect_j, defect_sum);
        near(
            path.energy.material.absolute_endpoint_work_defect_j,
            absolute_sum,
        );
        near(path.energy.material.plastic_dissipation_j, plastic_sum);
        near(path.energy.material.damage_release_j, damage_sum);
        assert!(absolute_sum <= controls.work.maximum_endpoint_work_defect_j);
        assert!(path.endpoint.solid.element_states()[0].damage() > 0.75);
        assert_eq!(
            path.endpoint.solid.element_states(),
            replay.element_states()
        );
        assert_eq!(
            path.endpoint.solid.displacement_m(),
            replay.displacement_m()
        );
        near(
            path.endpoint.equilibrium.reactions_n[3],
            VOLUME * g_old * 2.0 / target[0],
        );
        let held = path
            .endpoint
            .solid
            .advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(
            held.endpoint.solid.element_states(),
            path.endpoint.solid.element_states()
        );
        assert_eq!(held.energy.material.plastic_dissipation_j, 0.0);
        assert_eq!(held.energy.material.damage_release_j, 0.0);
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.element_states()[0].damage(), 0.0);
    });
}

#[test]
fn adaptive_free_motion_refines_kinetic_loss_with_original_mass_and_exact_physical_time() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let v0 = [0.25, -0.125, 0.0625];
        let motion = FiniteDuctileTetMotion::try_new(solid, &[RHO], &[v0; 4], cx).unwrap();
        let duration = 0.125;
        let acceleration = [0.5, -0.25, 0.125];
        let norm2 = acceleration.iter().map(|v| v * v).sum::<f64>();
        let mass = RHO * VOLUME;
        let one_step_loss = 0.5 * mass * norm2 * duration * duration;
        let limits = FiniteDuctileMotionPathControls {
            path: path_controls(),
            maximum_numerical_kinetic_loss_j: one_step_loss / 16.0,
        };
        let body = acceleration.map(|v| RHO * v);
        let conditions = DuctileStepConditions {
            duration_s: duration,
            ..conditions()
        };
        let path = motion
            .advance_path(
                &[],
                &[],
                [body; 2],
                conditions,
                IncrementSettings::default(),
                limits,
                cx,
            )
            .unwrap();
        assert!(path.records.len() >= 16 && path.rejected_trials >= 4);
        let mut t = 0.0;
        let mut dt2 = 0.0;
        let mut v = v0;
        for record in &path.records {
            let dt = record.duration_s;
            let account = record.energy.inertia.unwrap();
            let next = core::array::from_fn::<_, 3, _>(|i| v[i] + acceleration[i] * dt);
            let inertial_work = mass * dt * (0..3).map(|i| acceleration[i] * next[i]).sum::<f64>();
            near(account.reference_mass_kg, mass);
            near(account.inertial_work_j, inertial_work);
            near(account.numerical_loss_j, 0.5 * mass * norm2 * dt * dt);
            assert!(
                account.numerical_loss_j <= limits.maximum_numerical_kinetic_loss_j * dt / duration
            );
            near(
                record.energy.material.dead_load_endpoint_work_j,
                inertial_work,
            );
            near(record.energy.material.internal_endpoint_work_j, 0.0);
            near(record.energy.balance_error_j, 0.0);
            t += dt;
            dt2 += dt * dt;
            v = next;
        }
        near(t, duration);
        let kinetic = path.energy.inertia.unwrap();
        near(kinetic.numerical_loss_j, 0.5 * mass * norm2 * dt2);
        assert!(kinetic.numerical_loss_j <= limits.maximum_numerical_kinetic_loss_j);
        let k0 = 0.5 * mass * v0.iter().map(|v| v * v).sum::<f64>();
        let k1 = 0.5 * mass * v.iter().map(|v| v * v).sum::<f64>();
        near(kinetic.kinetic_before_j, k0);
        near(kinetic.kinetic_after_j, k1);
        near(kinetic.kinetic_change_j, k1 - k0);
        near(kinetic.inertial_work_j, k1 - k0 + kinetic.numerical_loss_j);
        for (u, velocity) in path
            .endpoint
            .motion
            .solid()
            .displacement_m()
            .iter()
            .zip(path.endpoint.motion.velocity_mps())
        {
            for i in 0..3 {
                near(velocity[i], v0[i] + acceleration[i] * duration);
                // Exact backward-Euler acceleration error for this actual mesh
                // of physical time intervals, independently of the controller.
                near(
                    u[i],
                    v0[i] * duration + 0.5 * acceleration[i] * (duration * duration + dt2),
                );
            }
        }
        assert_eq!(
            path.endpoint.motion.solid().element_states()[0].native_state(),
            material.initial_state().native_state()
        );
        assert!(core::ptr::eq(
            motion.reference_density_kg_m3(),
            path.endpoint.motion.reference_density_kg_m3()
        ));
        let replay = motion
            .advance_path(
                &[],
                &[],
                [body; 2],
                conditions,
                IncrementSettings::default(),
                limits,
                cx,
            )
            .unwrap();
        assert_eq!(path.records, replay.records);
        assert_eq!(
            path.endpoint.motion.velocity_mps(),
            replay.endpoint.motion.velocity_mps()
        );
        assert_eq!(motion.velocity_mps(), &[v0; 4]);
        assert_eq!(motion.solid().displacement_m(), &[[0.0; 3]; 4]);
    });
}

fn rate_model(width: f64) -> FiniteDuctileJ2<JohnsonCookInitiation> {
    let base = model(width, 0.04);
    let locus = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: ONSET,
            d2: 0.0,
            d3: 0.0,
            d4: 0.05,
            d5: 0.0,
            reference_strain_rate_s_inv: 1.0,
            transition_temperature_k: 299.0,
            melting_temperature_k: 301.0,
        },
        DuctileDomain {
            triaxiality: [-20.0, 20.0],
            strain_rate_s_inv: [1e-8, 10.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    FiniteDuctileJ2::new(base.native().clone(), locus, base.settings()).unwrap()
}

#[derive(Clone, Copy)]
struct Scalar {
    u: f64,
    v: f64,
    alpha: f64,
    utilization: f64,
    g: f64,
    native_storage: f64,
}
impl Scalar {
    // Independent axisymmetric logarithmic-J2/JC reduction; no production law
    // or material trial is called by this scalar equilibrium reference.
    fn at(self, u: f64, dt: f64, width: f64) -> (Self, f64) {
        let log_j = (1.0 + u).ln();
        let alpha = self.alpha.max((1000.0 * log_j - 3.0) / 1500.0);
        let delta = alpha - self.alpha;
        let q = (1000.0 * log_j - 1500.0 * self.alpha).min(3.0);
        let mut utilization = self.utilization;
        let mut onset = alpha;
        if self.utilization >= 1.0 {
            onset = self.alpha;
        } else if delta > 0.0 {
            let fracture = ONSET * (1.0 + 0.05 * (delta / dt).ln());
            utilization += delta / fracture;
            if utilization > 1.0 {
                onset = self.alpha + (1.0 - self.utilization) * fracture;
            }
        }
        let g = self.g * (-width * 3.0 * (alpha - onset).max(0.0) / GP).exp();
        let piola = g * (1000.0 * log_j + 2.0 * q / 3.0) / (1.0 + u);
        (
            Self {
                u,
                v: (u - self.u) / dt,
                alpha,
                utilization,
                g,
                native_storage: q * q / 3000.0 + 500.0 * log_j * log_j,
            },
            piola,
        )
    }
}

#[test]
fn adaptive_free_coordinate_matches_independent_rate_sensitive_finite_momentum_reference() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let width = bands[0].band_width_m();
        let material = rate_model(width);
        let materials = [&material];
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let effective_mass = RHO * VOLUME / 10.0;
        for duration in [0.125, 0.25] {
            let solid = FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let motion =
                FiniteDuctileTetMotion::try_new(solid, &[RHO], &[[0.0; 3]; 4], cx).unwrap();
            let limits = FiniteDuctileMotionPathControls {
                path: FiniteDuctilePathControls {
                    work: FinitePlasticStepControls {
                        maximum_endpoint_work_defect_j: 2e-3,
                        ..controls()
                    },
                    ..path_controls()
                },
                maximum_numerical_kinetic_loss_j: 1e-3,
            };
            let path = motion
                .advance_path(
                    &supports,
                    &[(3, [0.0, 3.0])],
                    [[0.0; 3]; 2],
                    DuctileStepConditions {
                        duration_s: duration,
                        ..conditions()
                    },
                    IncrementSettings::default(),
                    limits,
                    cx,
                )
                .unwrap();
            assert!(path.rejected_trials > 0 && path.records.len() > 1);
            let mut old = Scalar {
                u: 0.0,
                v: 0.0,
                alpha: 0.0,
                utilization: 0.0,
                g: 1.0,
                native_storage: 0.0,
            };
            for record in &path.records {
                let dt = record.duration_s;
                let force = 3.0 * record.fraction;
                let residual = |u: f64| {
                    let (trial, p) = old.at(u, dt, width);
                    VOLUME * p + effective_mass * (trial.v - old.v) / dt - force
                };
                let mut lower = old.u;
                let mut upper =
                    (((3.0 + 1500.0 * (old.alpha + 9.0 * dt)) / 1000.0).exp() - 1.0).min(0.05);
                assert!(residual(lower) <= 1e-10 && residual(upper) > 0.0);
                for _ in 0..90 {
                    let middle = 0.5 * (lower + upper);
                    if residual(middle) > 0.0 {
                        upper = middle;
                    } else {
                        lower = middle;
                    }
                }
                let (new, p) = old.at(0.5 * (lower + upper), dt, width);
                let stored_old = old.g * old.native_storage;
                let stored = new.g * new.native_storage;
                let plastic = new.g * 3.0 * (new.alpha - old.alpha);
                let damage = (old.g - new.g) * old.native_storage;
                let work = VOLUME * p * (new.u - old.u);
                let defect = work - VOLUME * ((stored - stored_old) + plastic + damage);
                let inertia = record.energy.inertia.unwrap();
                near(record.energy.material.stored_after_j, VOLUME * stored);
                near(
                    record.energy.material.plastic_dissipation_j,
                    VOLUME * plastic,
                );
                near(record.energy.material.damage_release_j, VOLUME * damage);
                near(record.energy.material.internal_endpoint_work_j, work);
                near(record.energy.material.endpoint_work_defect_j, defect);
                near(
                    record.energy.material.dead_load_endpoint_work_j,
                    force * (new.u - old.u),
                );
                near(record.energy.material.support_endpoint_work_j, 0.0);
                near(
                    inertia.kinetic_after_j,
                    0.5 * effective_mass * new.v * new.v,
                );
                near(
                    inertia.numerical_loss_j,
                    0.5 * effective_mass * (new.v - old.v).powi(2),
                );
                near(record.maximum_damage, 1.0 - new.g);
                near(record.energy.balance_error_j, 0.0);
                old = new;
            }
            let endpoint = &path.endpoint.motion;
            near(endpoint.solid().displacement_m()[1][0], old.u);
            near(endpoint.velocity_mps()[1][0], old.v);
            near(
                endpoint.solid().element_states()[0].native_state().alpha,
                old.alpha,
            );
            near(endpoint.solid().element_states()[0].integrity(), old.g);
            assert!(endpoint.solid().element_states()[0].damage() > 0.05);
            assert!(
                path.energy.material.absolute_endpoint_work_defect_j
                    <= limits.path.work.maximum_endpoint_work_defect_j
            );
            assert!(
                path.energy.inertia.unwrap().numerical_loss_j
                    <= limits.maximum_numerical_kinetic_loss_j
            );
        }
    });
}

#[test]
fn path_domain_initial_equilibrium_and_partial_progress_refusals_leave_the_original_state() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let supports = affine(&NODES, stretch(0.08), [0.0; 3], 0);
        assert!(matches!(
            solid.advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                DuctileStepConditions {
                    temperatures_k: DuctileTemperatures::Uniform(301.0),
                    ..conditions()
                },
                IncrementSettings::default(),
                path_controls(),
                cx,
            ),
            Err(FiniteDuctilePathError::Endpoint(
                FiniteDuctileTetError::Material {
                    source: FiniteDuctileError::Initiation(DuctileError::Domain {
                        quantity: "temperature_k",
                        ..
                    }),
                    ..
                }
            ))
        ));
        let free: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        assert!(matches!(
            solid.advance_path(
                &free,
                &[(3, [0.1, 0.2])],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                path_controls(),
                cx,
            ),
            Err(FiniteDuctilePathError::InitialEquilibrium {
                maximum_residual_n: 0.1,
                ..
            })
        ));
        let partial = FiniteDuctilePathControls {
            initial_fraction: 0.01,
            maximum_accepted_steps: 1,
            work: controls(),
            ..path_controls()
        };
        assert!(matches!(
            solid.advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                partial,
                cx,
            ),
            Err(FiniteDuctilePathError::WorkLimit {
                accepted_steps: 1,
                ..
            })
        ));
        let strict = FiniteDuctilePathControls {
            initial_fraction: 0.25,
            minimum_fraction: 0.0625,
            work: FinitePlasticStepControls {
                maximum_endpoint_work_defect_j: 0.0,
                ..controls()
            },
            ..path_controls()
        };
        assert!(matches!(
            solid.advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                strict,
                cx,
            ),
            Err(FiniteDuctilePathError::MinimumInterval {
                cause: FiniteDuctileRefinementCause::Endpoint(FiniteDuctileTetError::Energy {
                    quantity: "sum of absolute element integration defects",
                    ..
                }),
                ..
            })
        ));
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.element_states(), &[material.initial_state()]);
        gate.request();
        assert!(matches!(
            solid.advance_path(
                &supports,
                &[],
                [[0.0; 3]; 2],
                conditions(),
                IncrementSettings::default(),
                path_controls(),
                cx,
            ),
            Err(FiniteDuctilePathError::Endpoint(
                FiniteDuctileTetError::Solve(fs_solid::plastic3::SmallStrainTetError::Cancelled)
            ))
        ));
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.element_states(), &[material.initial_state()]);
    });
}
