//! G0/G1/G3 batteries for moving clearance and chamber-volume receipts.
//!
//! Complete time/spatial covers carry the proof.  Analytic values are
//! independent falsifiers and G1 comparison oracles; they never replace the
//! named-chart closure requirement.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Motor;
use fs_geom::fixtures::{BoxChart, SphereChart};
use fs_geom::{Aabb, Point3};
use fs_ivl::Interval;
use fs_motion::{
    ChamberChartFamily, ChamberDefinition, ChamberVolumeErrors, ChamberVolumeFunction,
    ClearanceConfig, ClearanceDecision, ClearanceErrors, ClearanceLowerEvidence, ClearanceOracle,
    ClearanceSidedness, ClearanceWitnessEvidence, IdealWankelVolumeOracle, MotionError, ProofState,
    ScrewParams, SpacetimeChart, SphereClearanceProxy, SpherePairClearanceOracle,
    overlap_inradius_witness, screw_tube, separation_over,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC1EA_AA7C_E001,
                kernel_id: 0xC1EA,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn with_cancelled_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC1EA_AA7C_CA11,
                kernel_id: 0xCA11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn sphere_pair(
    witnesses_enabled: bool,
) -> (
    SpacetimeChart<SphereChart>,
    SpacetimeChart<SphereChart>,
    SpherePairClearanceOracle,
) {
    let domain = Interval::new(-1.0, 1.0);
    let moving_center = Point3::new(0.4, 0.0, 0.0);
    let fixed_center = Point3::new(0.9, 0.0, 0.0);
    let moving_tube = screw_tube(
        &ScrewParams {
            axis: [0.0, 0.0, 1.0],
            center: [0.0, 0.0, 0.0],
            omega: 1.0,
            axial_velocity: 0.0,
            base_pose: Motor::identity(),
        },
        domain,
        10,
        8,
    )
    .expect("rotating eccentric sphere tube");
    let fixed_tube = screw_tube(
        &ScrewParams {
            axis: [0.0, 0.0, 1.0],
            center: [0.0, 0.0, 0.0],
            omega: 0.0,
            axial_velocity: 0.0,
            base_pose: Motor::identity(),
        },
        domain,
        4,
        2,
    )
    .expect("stationary sphere tube");
    let a = SpacetimeChart::new(
        SphereChart {
            center: moving_center,
            radius: 0.1,
        },
        moving_tube,
    );
    let b = SpacetimeChart::new(
        SphereChart {
            center: fixed_center,
            radius: 0.1,
        },
        fixed_tube,
    );
    let oracle = SpherePairClearanceOracle::new(
        SphereClearanceProxy {
            center: moving_center,
            radius_m: 0.1,
            errors: ClearanceErrors {
                chart_conversion_m: 1.0e-6,
                spatial_discretization_m: 2.0e-6,
                motion_tube_m: 3.0e-6,
                optimization_m: 0.0,
            },
        },
        SphereClearanceProxy {
            center: fixed_center,
            radius_m: 0.1,
            errors: ClearanceErrors {
                chart_conversion_m: 2.0e-6,
                spatial_discretization_m: 3.0e-6,
                motion_tube_m: 4.0e-6,
                optimization_m: 0.0,
            },
        },
        5.0e-6,
        witnesses_enabled,
    )
    .expect("sphere-pair evidence policy");
    (a, b, oracle)
}

struct PanicClearanceOracle;

impl ClearanceOracle<SphereChart, SphereChart> for PanicClearanceOracle {
    fn lower_bound_over(
        &self,
        _a: &SpacetimeChart<SphereChart>,
        _b: &SpacetimeChart<SphereChart>,
        _span: Interval,
        _cx: &Cx<'_>,
    ) -> Result<ClearanceLowerEvidence, MotionError> {
        panic!("a pre-cancelled clearance query must not enter its oracle")
    }

    fn witness_at(
        &self,
        _a: &SpacetimeChart<SphereChart>,
        _b: &SpacetimeChart<SphereChart>,
        _time: f64,
        _cx: &Cx<'_>,
    ) -> Result<Option<ClearanceWitnessEvidence>, MotionError> {
        panic!("a pre-cancelled clearance query must not enter its oracle")
    }
}

#[test]
fn pre_cancelled_clearance_refuses_before_oracle_work() {
    let (cam, follower, _) = sphere_pair(false);
    with_cancelled_cx(|cx| {
        let result = separation_over(
            &cam,
            &follower,
            Interval::new(-1.0, 1.0),
            &PanicClearanceOracle,
            ClearanceConfig::default(),
            cx,
        );
        assert!(matches!(result, Err(MotionError::Cancelled)));
    });
}

#[test]
fn rotating_cam_follower_has_two_sided_clearance_and_verified_witness() {
    with_cx(|cx| {
        let (cam, follower, oracle) = sphere_pair(true);
        let config = ClearanceConfig {
            value_tolerance_m: 5.0e-2,
            time_tolerance: 1.0e-5,
            max_subdivisions: 512,
        };
        let receipt = separation_over(
            &cam,
            &follower,
            Interval::new(-1.0, 1.0),
            &oracle,
            config,
            cx,
        )
        .expect("two-sided rotating clearance");

        // Independent center-distance derivation: at t=0 the centers are
        // collinear at radii 0.4 and 0.9, so the surface gap is 0.5 - 0.2.
        let analytic_minimum_m = 0.3;
        assert_eq!(receipt.sidedness, ClearanceSidedness::TwoSided);
        assert_eq!(receipt.decision, ClearanceDecision::Enclosure);
        assert!(receipt.lower_m <= analytic_minimum_m);
        let upper = receipt.upper_m.expect("feasible upper witness");
        assert!(upper >= analytic_minimum_m);
        assert!(upper - receipt.lower_m <= config.value_tolerance_m);

        let time = receipt.witness_time.expect("witness time");
        let dx = 0.4 * fs_math::det::cos(time) - 0.9;
        let dy = 0.4 * fs_math::det::sin(time);
        let witnessed_gap = fs_math::det::sqrt(dx * dx + dy * dy) - 0.2;
        let witness_error = receipt
            .errors
            .upper
            .expect("upper error accounting")
            .total_upper()
            .expect("finite upper error");
        assert!(witnessed_gap <= upper);
        assert!(upper <= witnessed_gap + witness_error + 1.0e-10);
        assert_eq!(receipt.witness_kind, Some("sphere-center-distance"));
        println!(
            "{{\"case\":\"rotating-cam-follower\",\"range_m\":[{:.17e},{:.17e}],\"witness_time\":{:.17e},\"errors_m\":{{\"conversion\":{:.17e},\"spatial\":{:.17e},\"motion\":{:.17e},\"optimization\":{:.17e}}},\"subdivisions\":{}}}",
            receipt.lower_m,
            upper,
            time,
            receipt.errors.lower.chart_conversion_m,
            receipt.errors.lower.spatial_discretization_m,
            receipt.errors.lower.motion_tube_m,
            receipt.errors.lower.optimization_m,
            receipt.subdivisions,
        );
    });
}

#[test]
fn disabled_witness_search_is_explicitly_one_sided() {
    with_cx(|cx| {
        let (cam, follower, oracle) = sphere_pair(false);
        let receipt = separation_over(
            &cam,
            &follower,
            Interval::new(-1.0, 1.0),
            &oracle,
            ClearanceConfig {
                value_tolerance_m: 5.0e-2,
                time_tolerance: 1.0e-3,
                max_subdivisions: 32,
            },
            cx,
        )
        .expect("lower-only rotating clearance");

        assert_eq!(receipt.sidedness, ClearanceSidedness::LowerOnly);
        assert_eq!(receipt.decision, ClearanceDecision::Unknown);
        assert!(receipt.upper_m.is_none());
        assert!(receipt.witness_time.is_none());
        assert!(receipt.witness_kind.is_none());
        assert!(receipt.errors.upper.is_none());
        assert_eq!(receipt.admitted_witnesses, 0);
        assert!(receipt.lower_m <= 0.3);
        println!(
            "{{\"case\":\"one-sided-clearance\",\"lower_m\":{:.17e},\"attempts\":{},\"admitted\":{}}}",
            receipt.lower_m, receipt.witness_attempts, receipt.admitted_witnesses,
        );
    });
}

#[test]
fn overlap_inradius_witness_is_not_penetration_depth() {
    let witness = overlap_inradius_witness(
        Point3::new(0.0, 0.0, 0.0),
        Interval::new(-0.31, -0.29),
        Interval::new(-0.22, -0.20),
    )
    .expect("finite exact-SDF bands")
    .expect("common-interior witness");
    assert!(witness.inradius_lower_m <= 0.2);
    assert!(witness.inradius_lower_m > 0.19);
    assert!(
        overlap_inradius_witness(
            Point3::new(0.0, 0.0, 0.0),
            Interval::new(-0.1, 0.01),
            Interval::new(-0.2, -0.1),
        )
        .expect("finite inconclusive bands")
        .is_none()
    );
}

#[derive(Debug, Clone, Copy)]
struct SliderCrankFamily {
    width_m: f64,
    depth_m: f64,
    deck_clearance_m: f64,
    crank_m: f64,
    rod_m: f64,
}

impl SliderCrankFamily {
    fn height_at(self, angle: f64) -> Result<f64, MotionError> {
        let sine = fs_math::det::sin(angle);
        let cosine = fs_math::det::cos(angle);
        let radicand = self.rod_m * self.rod_m - self.crank_m * self.crank_m * sine * sine;
        if radicand <= 0.0 {
            return Err(MotionError::InvalidGeometry {
                what: "slider-crank rod must remain longer than transverse crank offset",
            });
        }
        Ok(
            self.deck_clearance_m + self.crank_m * (1.0 - cosine) + self.rod_m
                - fs_math::det::sqrt(radicand),
        )
    }

    fn closed_form_volume(self, angle: f64) -> Result<f64, MotionError> {
        Ok(self.width_m * self.depth_m * self.height_at(angle)?)
    }
}

impl ChamberChartFamily for SliderCrankFamily {
    type Chamber = BoxChart;

    fn chart_at(
        &self,
        _definition: &ChamberDefinition,
        angle_radians: f64,
        _cx: &Cx<'_>,
    ) -> Result<Self::Chamber, MotionError> {
        let height = self.height_at(angle_radians)?;
        Ok(BoxChart {
            aabb: Aabb::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(self.width_m, self.depth_m, height),
            ),
        })
    }
}

#[test]
fn slider_crank_volume_quadrature_encloses_closed_form_and_logs_errors() {
    with_cx(|cx| {
        let family = SliderCrankFamily {
            width_m: 0.18,
            depth_m: 0.12,
            deck_clearance_m: 0.10,
            crank_m: 0.08,
            rod_m: 0.25,
        };
        let definition = ChamberDefinition::new(
            "slider-crank cylinder",
            vec![
                "cylinder wall".to_owned(),
                "cylinder head".to_owned(),
                "piston crown".to_owned(),
            ],
            "ideal zero-clearance rings; no ports or crevice volume",
            ProofState::Proven,
        )
        .expect("closed ideal cylinder definition");
        let maximum_height = family
            .height_at(std::f64::consts::PI)
            .expect("bottom-dead-center height");
        let errors = ChamberVolumeErrors {
            chart_conversion_m3: 1.0e-7,
            motion_tube_m3: 2.0e-7,
            boundary_closure_m3: 3.0e-7,
            model_form_m3: 4.0e-7,
        };
        let volume = ChamberVolumeFunction::new(
            family,
            definition,
            Aabb::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(family.width_m, family.depth_m, maximum_height + 1.0e-4),
            ),
            1.5e-2,
            errors,
        )
        .expect("slider-crank volume function");

        for angle in [0.0, 0.7, 1.9, std::f64::consts::PI] {
            let receipt = volume.at(angle, cx).expect("certified cylinder volume");
            let closed_form = family
                .closed_form_volume(angle)
                .expect("slider-crank G1 oracle");
            assert!(receipt.spatial_quadrature_m3.contains(closed_form));
            assert!(receipt.volume_m3.contains(closed_form));
            assert!(receipt.volume_m3.lo() <= receipt.spatial_quadrature_m3.lo());
            assert!(receipt.volume_m3.hi() >= receipt.spatial_quadrature_m3.hi());
            println!(
                "{{\"case\":\"slider-crank-volume\",\"theta\":{angle:.17e},\"oracle_m3\":{closed_form:.17e},\"quadrature_m3\":[{:.17e},{:.17e}],\"inflated_m3\":[{:.17e},{:.17e}],\"band_cells\":{},\"errors_m3\":{{\"conversion\":{:.17e},\"motion\":{:.17e},\"closure\":{:.17e},\"model\":{:.17e}}}}}",
                receipt.spatial_quadrature_m3.lo(),
                receipt.spatial_quadrature_m3.hi(),
                receipt.volume_m3.lo(),
                receipt.volume_m3.hi(),
                receipt.band_cells,
                receipt.errors.chart_conversion_m3,
                receipt.errors.motion_tube_m3,
                receipt.errors.boundary_closure_m3,
                receipt.errors.model_form_m3,
            );
        }
    });
}

#[derive(Debug, Clone, Copy)]
struct EquivalentVolumePrism {
    oracle: IdealWankelVolumeOracle,
    width_m: f64,
    depth_m: f64,
}

impl ChamberChartFamily for EquivalentVolumePrism {
    type Chamber = BoxChart;

    fn chart_at(
        &self,
        _definition: &ChamberDefinition,
        angle_radians: f64,
        _cx: &Cx<'_>,
    ) -> Result<Self::Chamber, MotionError> {
        let height = self.oracle.volume_at(angle_radians)? / (self.width_m * self.depth_m);
        Ok(BoxChart {
            aabb: Aabb::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(self.width_m, self.depth_m, height),
            ),
        })
    }
}

#[test]
fn ideal_wankel_g1_oracle_is_enclosed_by_manufactured_quadrature_fixture() {
    with_cx(|cx| {
        let oracle = IdealWankelVolumeOracle {
            minimum_volume_m3: 1.4e-5,
            eccentricity_m: 0.015,
            generating_radius_m: 0.105,
            housing_parallel_transfer_m: 0.002,
            rotor_parallel_transfer_m: 0.001,
            housing_depth_m: 0.080,
            volume_phase_radians: std::f64::consts::FRAC_PI_6,
        };
        let family = EquivalentVolumePrism {
            oracle,
            width_m: 0.10,
            depth_m: 0.08,
        };
        let definition = ChamberDefinition::new(
            "manufactured equivalent-volume prism",
            vec![
                "lower prism face".to_owned(),
                "upper prism face".to_owned(),
                "four lateral prism faces".to_owned(),
            ],
            "G1 oracle adapter only; not Wankel bore or rotor-flank geometry",
            ProofState::Proven,
        )
        .expect("manufactured closed chamber");
        let amplitude = 0.5
            * fs_math::det::sqrt(3.0)
            * oracle.eccentricity_m
            * (2.0 * oracle.generating_radius_m
                + oracle.housing_parallel_transfer_m
                + oracle.rotor_parallel_transfer_m)
            * oracle.housing_depth_m;
        let maximum_volume = oracle.minimum_volume_m3 + 2.0 * amplitude;
        let function = ChamberVolumeFunction::new(
            family,
            definition,
            Aabb::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(
                    family.width_m,
                    family.depth_m,
                    maximum_volume / (family.width_m * family.depth_m) + 1.0e-4,
                ),
            ),
            1.0e-2,
            ChamberVolumeErrors::default(),
        )
        .expect("manufactured Wankel-oracle comparison function");

        for angle in [0.0, 0.8, 2.4, 4.8] {
            let expected = oracle.volume_at(angle).expect("ideal Wankel G1 formula");
            let receipt = function.at(angle, cx).expect("manufactured quadrature");
            assert!(receipt.volume_m3.contains(expected));
            println!(
                "{{\"case\":\"ideal-wankel-g1-oracle\",\"theta\":{angle:.17e},\"oracle_m3\":{expected:.17e},\"quadrature_m3\":[{:.17e},{:.17e}],\"fixture\":\"equivalent-volume-prism-no-geometry-claim\"}}",
                receipt.volume_m3.lo(),
                receipt.volume_m3.hi(),
            );
        }
    });
}

#[derive(Debug, Clone, Copy)]
struct MustNotBuildWankelChart;

impl ChamberChartFamily for MustNotBuildWankelChart {
    type Chamber = BoxChart;

    fn chart_at(
        &self,
        _definition: &ChamberDefinition,
        _angle_radians: f64,
        _cx: &Cx<'_>,
    ) -> Result<Self::Chamber, MotionError> {
        panic!("an unproven actual Wankel closure must refuse before chart construction")
    }
}

#[test]
fn actual_wankel_volume_refuses_without_bore_flank_and_seal_closure() {
    with_cx(|cx| {
        let definition = ChamberDefinition::new(
            "actual Wankel chamber",
            vec![
                "visibility-trimmed epitrochoid bore".to_owned(),
                "conjugate rotor flank".to_owned(),
                "finite apex and side seals".to_owned(),
                "side housings".to_owned(),
            ],
            "finite seals and clearances; ports closed",
            ProofState::Unknown,
        )
        .expect("named but not yet proven Wankel chamber");
        let function = ChamberVolumeFunction::new(
            MustNotBuildWankelChart,
            definition,
            Aabb::new(Point3::new(-1.0, -1.0, -0.1), Point3::new(1.0, 1.0, 0.1)),
            0.05,
            ChamberVolumeErrors::default(),
        )
        .expect("typed fail-closed volume function");
        assert!(matches!(
            function.at(0.0, cx),
            Err(MotionError::InvalidEvidence {
                what: "chamber boundary closure must be Proven before volume authority"
            })
        ));
    });
}
