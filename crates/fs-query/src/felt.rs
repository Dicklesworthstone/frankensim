//! Felt thickness-field chart ingestion: field -> chart + receipt (music
//! program bead `frankensim-music-v8-root-3ez8g.3.4`; ingest law).
//!
//! The piano hammer's contact GEOMETRY is a felt thickness field over a
//! wood core. The sampler side already exists and is exactly right
//! (fs-contact's `sample_finite_gap_from_chart` + the convergence-checked
//! dense response family bound by `response_identity`); this module is
//! the missing ingest side: a validated thickness field, a provenance
//! receipt with thickness statistics, and a [`Chart`] the certified gap
//! sampler can consume.
//!
//! AUTHORITY (binding): a felt thickness field is `Estimate` BY
//! CONSTRUCTION — profiles are authored crown models or values
//! transcribed from a manufacturer source, and the chart's implicit
//! field carries no directed-rounding proof — so the chart publishes
//! `TraceStepClaim::NoClaim` and per-sample `Estimate` certificates whose
//! band is the station uncertainty half-width. Gap sampling therefore
//! runs under `AllowEstimate` and the receipt's authority is `Estimate`,
//! propagated without promotion (the `AxisymmetricChart` precedent).
//!
//! CORPUS LAW (lineage): every field names its source and transcription
//! basis. The committed manufacturer fixture
//! (`data/felt/steinway-us5125310-strip-taper.tsv`) is TRANSCRIBED, not
//! digitized, from the text of US patent 5,125,310 (Steinway, 1992) —
//! US patent text is a public record, so retention is lawful — and
//! records exactly the stated endpoint values; the linear interpolation
//! between stations is this module's declared rule, not the patent's.
//! The felt CONSTITUTIVE card is bead 87zbd's (its licensed-coupon hunt
//! is the registered absence row `acoustic-absent-felt-coupon`);
//! felt-as-absorber is a different bead and a different claim.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};

/// One station of a felt thickness profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeltStation {
    /// Station coordinate: angle around the core [rad] for a crowned
    /// hammer profile, or arc length along a strip [m] for a
    /// manufacturer strip section (the units declaration names which).
    pub coordinate: f64,
    /// Felt thickness at the station [m].
    pub thickness_m: f64,
    /// Uncertainty half-width of the thickness [m] (digitization or
    /// transcription band; 0 for an authored analytic fixture).
    pub half_width_m: f64,
}

/// Units declaration for a felt field — an explicit, refusable input
/// (undocumented units are a refusal, never a guess).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeltCoordinateUnits {
    /// Angle around the hammer core [rad]; thickness in metres.
    RadiansAroundCore,
    /// Arc length along a felt strip [m]; thickness in metres.
    MetersAlongStrip,
}

/// A validated felt thickness field with lineage.
#[derive(Debug, Clone, PartialEq)]
pub struct FeltThicknessField {
    stations: Vec<FeltStation>,
    units: FeltCoordinateUnits,
    source_id: String,
    source_locator: String,
    basis: String,
}

/// Ingest receipt: digest, statistics, and the honest authority.
#[derive(Debug, Clone, PartialEq)]
pub struct FeltFieldReceipt {
    /// FNV-1a digest over the units tag, every station's exact bits,
    /// and the source identity (the fs-query digest idiom).
    pub digest: u64,
    /// Caller-supplied immutable source identity.
    pub source_id: String,
    /// Where in the source the values live (figure/column/claim text).
    pub source_locator: String,
    /// Transcription/digitization basis sentence.
    pub basis: String,
    /// Station count.
    pub stations: usize,
    /// Minimum thickness across stations [m].
    pub minimum_thickness_m: f64,
    /// Maximum thickness across stations [m] (the crown for a crowned
    /// profile).
    pub maximum_thickness_m: f64,
    /// Coordinate of the maximum-thickness station (crown location).
    pub crown_coordinate: f64,
    /// Largest station uncertainty half-width [m].
    pub maximum_half_width_m: f64,
    /// Always [`NumericalKind::Estimate`] — recorded, never promoted.
    pub authority: NumericalKind,
}

/// Typed refusal from felt-field ingestion.
#[derive(Debug, Clone, PartialEq)]
pub enum FeltError {
    /// A structural input problem (empty source id, too few stations).
    Invalid {
        /// What was wrong.
        what: &'static str,
    },
    /// A station value is non-finite or a thickness is non-positive —
    /// a felt of zero or negative thickness is a data error, never a
    /// bare core.
    BadStation {
        /// Station index.
        index: usize,
        /// Offending field.
        field: &'static str,
    },
    /// Station coordinates must be strictly increasing.
    UnsortedStations {
        /// First offending index.
        index: usize,
    },
    /// The chart constructor's geometry is invalid.
    BadGeometry {
        /// Offending field.
        field: &'static str,
    },
}

impl FeltThicknessField {
    /// Validate and admit a thickness field.
    ///
    /// # Errors
    /// [`FeltError`] on an empty source identity, fewer than two
    /// stations, non-finite or non-positive values, or unsorted
    /// station coordinates.
    pub fn try_new(
        stations: Vec<FeltStation>,
        units: FeltCoordinateUnits,
        source_id: &str,
        source_locator: &str,
        basis: &str,
    ) -> Result<Self, FeltError> {
        if source_id.trim().is_empty() {
            return Err(FeltError::Invalid {
                what: "source_id must be non-empty",
            });
        }
        if source_locator.trim().is_empty() || basis.trim().is_empty() {
            return Err(FeltError::Invalid {
                what: "lineage (source_locator, basis) must be recorded",
            });
        }
        if stations.len() < 2 {
            return Err(FeltError::Invalid {
                what: "a field needs at least two stations",
            });
        }
        for (index, s) in stations.iter().enumerate() {
            if !s.coordinate.is_finite() {
                return Err(FeltError::BadStation {
                    index,
                    field: "coordinate",
                });
            }
            if !(s.thickness_m.is_finite() && s.thickness_m > 0.0) {
                return Err(FeltError::BadStation {
                    index,
                    field: "thickness_m",
                });
            }
            if !(s.half_width_m.is_finite() && s.half_width_m >= 0.0) {
                return Err(FeltError::BadStation {
                    index,
                    field: "half_width_m",
                });
            }
        }
        for index in 1..stations.len() {
            if stations[index].coordinate <= stations[index - 1].coordinate {
                return Err(FeltError::UnsortedStations { index });
            }
        }
        Ok(Self {
            stations,
            units,
            source_id: source_id.to_string(),
            source_locator: source_locator.to_string(),
            basis: basis.to_string(),
        })
    }

    /// The units declaration.
    #[must_use]
    pub fn units(&self) -> FeltCoordinateUnits {
        self.units
    }

    /// The stations.
    #[must_use]
    pub fn stations(&self) -> &[FeltStation] {
        &self.stations
    }

    /// Piecewise-linear thickness at `coordinate`, CLAMPED to the end
    /// stations outside coverage (the shoulder rule — documented, not
    /// silent extrapolation). Returns `(thickness_m, half_width_m)`.
    #[must_use]
    pub fn thickness_at(&self, coordinate: f64) -> (f64, f64) {
        let first = self.stations[0];
        let last = self.stations[self.stations.len() - 1];
        if coordinate <= first.coordinate {
            return (first.thickness_m, first.half_width_m);
        }
        if coordinate >= last.coordinate {
            return (last.thickness_m, last.half_width_m);
        }
        let mut upper = 1;
        while self.stations[upper].coordinate < coordinate {
            upper += 1;
        }
        let a = self.stations[upper - 1];
        let b = self.stations[upper];
        let f = (coordinate - a.coordinate) / (b.coordinate - a.coordinate);
        (
            a.thickness_m + f * (b.thickness_m - a.thickness_m),
            a.half_width_m + f * (b.half_width_m - a.half_width_m),
        )
    }

    /// Mint the ingest receipt: digest, thickness statistics, and the
    /// non-promotable `Estimate` authority.
    #[must_use]
    pub fn receipt(&self) -> FeltFieldReceipt {
        let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
        let mut absorb = |bytes: &[u8]| {
            for byte in bytes {
                digest ^= u64::from(*byte);
                digest = digest.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        absorb(match self.units {
            FeltCoordinateUnits::RadiansAroundCore => b"rad",
            FeltCoordinateUnits::MetersAlongStrip => b"m",
        });
        absorb(self.source_id.as_bytes());
        for s in &self.stations {
            absorb(&s.coordinate.to_bits().to_le_bytes());
            absorb(&s.thickness_m.to_bits().to_le_bytes());
            absorb(&s.half_width_m.to_bits().to_le_bytes());
        }
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        let mut crown = self.stations[0].coordinate;
        let mut widest = 0.0_f64;
        for s in &self.stations {
            minimum = minimum.min(s.thickness_m);
            if s.thickness_m > maximum {
                maximum = s.thickness_m;
                crown = s.coordinate;
            }
            widest = widest.max(s.half_width_m);
        }
        FeltFieldReceipt {
            digest,
            source_id: self.source_id.clone(),
            source_locator: self.source_locator.clone(),
            basis: self.basis.clone(),
            stations: self.stations.len(),
            minimum_thickness_m: minimum,
            maximum_thickness_m: maximum,
            crown_coordinate: crown,
            maximum_half_width_m: widest,
            authority: NumericalKind::Estimate,
        }
    }
}

impl FeltFieldReceipt {
    /// Per-field JSON log line (thickness statistics per ingested
    /// field — the bead's logging clause).
    #[must_use]
    pub fn debug_line(&self) -> String {
        format!(
            "{{\"suite\":\"fs-query\",\"case\":\"felt-field-receipt\",\"source\":\"{}\",\
             \"digest\":\"{:#018x}\",\"stations\":{},\"t_min_m\":{:.6e},\"t_max_m\":{:.6e},\
             \"crown_at\":{:.6},\"hw_max_m\":{:.3e},\"authority\":\"Estimate\"}}",
            self.source_id,
            self.digest,
            self.stations,
            self.minimum_thickness_m,
            self.maximum_thickness_m,
            self.crown_coordinate,
            self.maximum_half_width_m,
        )
    }
}

/// The hammer solid as a [`Chart`]: a wood core cylinder (axis along
/// `y` through `(0, y, z_axis)`) carrying the felt so the outer surface
/// sits at radius `core_radius_m + t(theta)`, intersected with the slab
/// `|y| <= half_width_m`. `theta` is measured from the `-z` (string
/// contact) direction, so the crown faces the string and the tangent
/// contact frame at the apex is the origin with outward normal `-z` —
/// the same frame the felt-island fixtures use.
#[derive(Debug, Clone, PartialEq)]
pub struct FeltThicknessChart {
    field: FeltThicknessField,
    core_radius_m: f64,
    half_width_m: f64,
    width_rounding_radius_m: f64,
    axis_height_m: f64,
    outer_radius_max_m: f64,
}

impl FeltThicknessChart {
    /// Build the chart from a crowned profile (units must be
    /// [`FeltCoordinateUnits::RadiansAroundCore`]).
    ///
    /// The core axis is placed at `z = core_radius_m + t(0)` so the
    /// apex of the felt touches the origin.
    ///
    /// # Errors
    /// [`FeltError::BadGeometry`] on a non-positive core radius or
    /// half width; [`FeltError::Invalid`] when the field's units are
    /// not radians-around-core.
    pub fn try_new(
        field: FeltThicknessField,
        core_radius_m: f64,
        half_width_m: f64,
        width_rounding_radius_m: f64,
    ) -> Result<Self, FeltError> {
        if field.units() != FeltCoordinateUnits::RadiansAroundCore {
            return Err(FeltError::Invalid {
                what: "a crowned chart needs a radians-around-core profile",
            });
        }
        if !(core_radius_m.is_finite() && core_radius_m > 0.0) {
            return Err(FeltError::BadGeometry {
                field: "core_radius_m",
            });
        }
        if !(half_width_m.is_finite() && half_width_m > 0.0) {
            return Err(FeltError::BadGeometry {
                field: "half_width_m",
            });
        }
        if !(width_rounding_radius_m.is_finite() && width_rounding_radius_m > 0.0) {
            return Err(FeltError::BadGeometry {
                field: "width_rounding_radius_m",
            });
        }
        let (apex_thickness, _) = field.thickness_at(0.0);
        let axis_height_m = core_radius_m + apex_thickness;
        let mut outer_radius_max_m = 0.0_f64;
        for s in field.stations() {
            outer_radius_max_m = outer_radius_max_m.max(core_radius_m + s.thickness_m);
        }
        Ok(Self {
            field,
            core_radius_m,
            half_width_m,
            width_rounding_radius_m,
            axis_height_m,
            outer_radius_max_m,
        })
    }

    /// The ingest receipt of the underlying field.
    #[must_use]
    pub fn field_receipt(&self) -> FeltFieldReceipt {
        self.field.receipt()
    }
}

impl Chart for FeltThicknessChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
            return ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            };
        }
        // Radial coordinates around the core axis (along y at height
        // axis_height_m): theta = 0 points at -z (the string).
        let dx = x.x;
        let dz = x.z - self.axis_height_m;
        let rho = (dx * dx + dz * dz).sqrt();
        let theta = dx.atan2(-dz);
        let (thickness, half_width) = self.field.thickness_at(theta);
        // Width rounding: real hammers are crowned ACROSS the width as
        // well, so the outer surface sags by y^2 / (2 R_w) toward the
        // cheeks — this keeps the contact patch compact (a
        // y-invariant cylinder would ride the sampler's boundary ring
        // at any approach).
        let sag = x.y * x.y / (2.0 * self.width_rounding_radius_m);
        let radial = rho - (self.core_radius_m + thickness - sag);
        let slab = x.y.abs() - self.half_width_m;
        // Sign-exact implicit for the modeled solid; NOT a distance
        // (angular stretch), so the certificate is an Estimate band of
        // the station uncertainty — never promoted.
        let signed = radial.max(slab);
        ChartSample {
            signed_distance: signed,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::estimate(signed - half_width, signed + half_width),
        }
    }

    fn support(&self) -> Aabb {
        let r = self.outer_radius_max_m;
        Aabb::new(
            Point3::new(-r, -self.half_width_m, self.axis_height_m - r),
            Point3::new(r, self.half_width_m, self.axis_height_m + r),
        )
    }

    fn name(&self) -> &'static str {
        "music/felt-thickness-field"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn cx_scope<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 11,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn crown_field() -> FeltThicknessField {
        // Authored crowned profile t(theta) = 9 mm - 4 mm * (theta/1)^2
        // sampled densely (authored analytic fixture: half-widths 0).
        let stations: Vec<FeltStation> = (0..41)
            .map(|i| {
                let theta = -1.0 + 0.05 * f64::from(i);
                FeltStation {
                    coordinate: theta,
                    thickness_m: 9.0e-3 - 4.0e-3 * theta * theta,
                    half_width_m: 0.0,
                }
            })
            .collect();
        FeltThicknessField::try_new(
            stations,
            FeltCoordinateUnits::RadiansAroundCore,
            "music/felt-crown/authored-parabolic/v1",
            "authored fixture (this test)",
            "authored analytic crown profile; Estimate by authorship",
        )
        .expect("authored crown admits")
    }

    #[test]
    fn gf_001_refusals_fire_by_name() {
        let ok = FeltStation {
            coordinate: 0.0,
            thickness_m: 5.0e-3,
            half_width_m: 0.0,
        };
        let ok2 = FeltStation {
            coordinate: 1.0,
            thickness_m: 4.0e-3,
            half_width_m: 0.0,
        };
        let mk = |st: Vec<FeltStation>| {
            FeltThicknessField::try_new(st, FeltCoordinateUnits::RadiansAroundCore, "t", "l", "b")
        };
        assert!(matches!(mk(vec![ok]), Err(FeltError::Invalid { .. })));
        assert!(matches!(
            mk(vec![
                FeltStation {
                    thickness_m: 0.0,
                    ..ok
                },
                ok2
            ]),
            Err(FeltError::BadStation {
                field: "thickness_m",
                ..
            })
        ));
        assert!(matches!(
            mk(vec![
                FeltStation {
                    thickness_m: f64::NAN,
                    ..ok
                },
                ok2
            ]),
            Err(FeltError::BadStation { .. })
        ));
        assert!(matches!(
            mk(vec![ok2, ok]),
            Err(FeltError::UnsortedStations { index: 1 })
        ));
        assert!(matches!(
            FeltThicknessField::try_new(
                vec![ok, ok2],
                FeltCoordinateUnits::RadiansAroundCore,
                "",
                "l",
                "b"
            ),
            Err(FeltError::Invalid { .. })
        ));
        assert!(matches!(
            FeltThicknessField::try_new(
                vec![ok, ok2],
                FeltCoordinateUnits::RadiansAroundCore,
                "t",
                " ",
                "b"
            ),
            Err(FeltError::Invalid { .. })
        ));
        // A strip field cannot build a crowned chart (units refusal —
        // undocumented/mismatched units never guess).
        let strip = FeltThicknessField::try_new(
            vec![ok, ok2],
            FeltCoordinateUnits::MetersAlongStrip,
            "t",
            "l",
            "b",
        )
        .expect("strip admits as a field");
        assert!(matches!(
            FeltThicknessChart::try_new(strip, 5.0e-3, 6.0e-3, 9.0e-3),
            Err(FeltError::Invalid { .. })
        ));
        let crown = crown_field();
        assert!(matches!(
            FeltThicknessChart::try_new(crown.clone(), -1.0, 6.0e-3, 9.0e-3),
            Err(FeltError::BadGeometry {
                field: "core_radius_m"
            })
        ));
        assert!(matches!(
            FeltThicknessChart::try_new(crown, 8.0e-3, 6.0e-3, 0.0),
            Err(FeltError::BadGeometry {
                field: "width_rounding_radius_m"
            })
        ));
        println!("{{\"suite\":\"fs-query\",\"case\":\"gf-001-refusals\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn gf_002_receipt_stats_and_digest_are_deterministic() {
        let field = crown_field();
        let receipt = field.receipt();
        assert_eq!(receipt.stations, 41);
        assert!((receipt.maximum_thickness_m - 9.0e-3).abs() < 1e-12);
        assert!((receipt.minimum_thickness_m - 5.0e-3).abs() < 1e-12);
        assert!(receipt.crown_coordinate.abs() < 1e-12);
        assert_eq!(receipt.authority, NumericalKind::Estimate);
        let again = field.receipt();
        assert_eq!(receipt.digest, again.digest, "digest must be stable");
        // A one-bit thickness change must move the digest.
        let mut perturbed = field.stations().to_vec();
        perturbed[20].thickness_m += 1.0e-12;
        let other = FeltThicknessField::try_new(
            perturbed,
            FeltCoordinateUnits::RadiansAroundCore,
            "music/felt-crown/authored-parabolic/v1",
            "authored fixture (this test)",
            "authored analytic crown profile; Estimate by authorship",
        )
        .expect("perturbed admits");
        assert_ne!(receipt.digest, other.receipt().digest);
        println!("{}", receipt.debug_line());
    }

    #[test]
    fn gf_003_chart_signs_and_shoulder_clamp() {
        let chart =
            FeltThicknessChart::try_new(crown_field(), 8.0e-3, 6.0e-3, 9.0e-3).expect("chart");
        cx_scope(|cx| {
            // Apex: the origin is ON the surface; just outside below,
            // just inside above.
            let below = chart.eval(Point3::new(0.0, 0.0, -1.0e-4), cx);
            let above = chart.eval(Point3::new(0.0, 0.0, 1.0e-4), cx);
            assert!(below.signed_distance > 0.0, "below apex is outside");
            assert!(above.signed_distance < 0.0, "above apex is inside");
            // The slab bound: same (x, z) but beyond the half width is
            // outside.
            let beyond = chart.eval(Point3::new(0.0, 7.0e-3, 1.0e-3), cx);
            assert!(beyond.signed_distance > 0.0, "beyond the width is outside");
            // Shoulder clamp: far around the core (theta ~ pi, above
            // the axis) the thickness holds the end-station value, so a
            // point just above the outer radius there is outside.
            let axis_z = 8.0e-3 + 9.0e-3;
            let shoulder_thickness = 9.0e-3 - 4.0e-3; // t(+-1), clamped
            let above_shoulder = chart.eval(
                Point3::new(0.0, 0.0, axis_z + 8.0e-3 + shoulder_thickness + 1.0e-4),
                cx,
            );
            assert!(above_shoulder.signed_distance > 0.0);
            // Certificates are Estimate bands containing the nominal.
            assert_eq!(below.error.kind, NumericalKind::Estimate);
            assert!(below.error.lo <= below.signed_distance);
            assert!(below.error.hi >= below.signed_distance);
            assert_eq!(chart.trace_step_claim(), TraceStepClaim::NoClaim);
        });
        println!("{{\"suite\":\"fs-query\",\"case\":\"gf-003-chart-signs\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn gf_004_committed_steinway_strip_ingests_with_stated_values() {
        // The committed manufacturer section: transcribed endpoint
        // values from US patent 5,125,310 (Steinway, 1992) — outer felt
        // strip height tapers 1 inch -> 1/8 inch over the ~44 inch
        // strip; under felt 1/4 inch -> 3/32 inch. Transcription of
        // stated values, so the half-width is the half-ULP of the
        // stated fraction (0 here: the fractions are exact statements;
        // the INTERPOLATION rule between endpoints is ours, not the
        // patent's, and the chart never consumes this strip field).
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("root")
            .to_path_buf();
        let text =
            std::fs::read_to_string(root.join("data/felt/steinway-us5125310-strip-taper.tsv"))
                .expect("committed strip taper");
        assert!(text.starts_with("# frankensim-felt-field-v1"));
        assert!(
            text.contains("US5125310"),
            "the lineage must name the patent"
        );
        let mut outer = Vec::new();
        let mut under = Vec::new();
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols[0] == "series" || cols.len() < 4 {
                continue;
            }
            let station = FeltStation {
                coordinate: cols[1].parse().expect("coordinate"),
                thickness_m: cols[2].parse().expect("thickness"),
                half_width_m: cols[3].parse().expect("half width"),
            };
            match cols[0] {
                "outer" => outer.push(station),
                "under" => under.push(station),
                other => panic!("unknown series {other}"),
            }
        }
        let outer_field = FeltThicknessField::try_new(
            outer,
            FeltCoordinateUnits::MetersAlongStrip,
            "music/felt-strip/steinway-us5125310/outer/v1",
            "US5125310 col. describing outer felt: height Ho from 1 inch to 1/8 inch over ~44 inch strip",
            "transcribed stated values (public patent text); linear interpolation between endpoints is the ingest rule, not the patent's",
        )
        .expect("outer strip admits");
        let receipt = outer_field.receipt();
        assert_eq!(receipt.stations, 2);
        assert!(
            (receipt.maximum_thickness_m - 0.0254).abs() < 1e-12,
            "1 inch"
        );
        assert!(
            (receipt.minimum_thickness_m - 3.175e-3).abs() < 1e-12,
            "1/8 inch"
        );
        let under_field = FeltThicknessField::try_new(
            under,
            FeltCoordinateUnits::MetersAlongStrip,
            "music/felt-strip/steinway-us5125310/under/v1",
            "US5125310 col. describing under felt: height Hu from 1/4 inch to 3/32 inch",
            "transcribed stated values (public patent text); linear interpolation between endpoints is the ingest rule, not the patent's",
        )
        .expect("under strip admits");
        let under_receipt = under_field.receipt();
        assert!(
            (under_receipt.maximum_thickness_m - 6.35e-3).abs() < 1e-12,
            "1/4 inch"
        );
        assert!(
            (under_receipt.minimum_thickness_m - 2.38125e-3).abs() < 1e-12,
            "3/32 inch"
        );
        println!("{}", receipt.debug_line());
        println!("{}", under_receipt.debug_line());
    }
}
