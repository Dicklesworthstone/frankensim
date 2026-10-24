//! Pack-backed material values for the fs-couple literature tests (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.32`).
//!
//! Every material number these tests use comes from a compiled pack under
//! `tests/fixtures/matdb`. It is either a seed pack compiled from
//! `data/matdb/seed-v1/<id>`, or, where no source carries the value, a
//! labelled authored-fixture pack compiled from `tests/fixtures/matdb/<id>/`
//! (its terms and observation read "authored test fixture; not data", and its
//! card's process is `authored-fixture`). Each `.fsmatpk` is the xtask
//! compiler's exact output for its source; `xtask/tests/piano_string_packs.rs`
//! refuses a stale one. Values resolve through `fs_material::state_point`, so
//! a missing, ambiguous or out-of-domain property refuses by name and nothing
//! defaults.
#![allow(dead_code)] // each test binary uses its own subset

use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint,
    UncertaintyModel,
};
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::state_point::{
    MaterialPropertySelection, OrthotropicPlaneElastic, OrthotropicRequirementSet,
    PoissonReciprocityRule, ReferenceRatioOrthotropicStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
    resolve_orthotropic_elastic_state_point_from_reference_ratios,
};
use fs_qty::Dims;

/// The process of every authored-fixture card.
pub const AUTHORED_FIXTURE: &str = "authored-fixture";
/// The label every authored-fixture source carries in its terms and observation.
pub const FIXTURE_LABEL: &str = "authored test fixture; not data";

/// Sitka spruce, FPL-GTR-282 ch. 5 at 12 % MC.
pub const SPRUCE_SITKA: &[u8] = include_bytes!("../fixtures/matdb/spruce-sitka-fpl-gtr282.fsmatpk");
/// Sitka spruce along-grain loss, Qiu 2026 (first mode, 62.36 Hz, 12 % MC).
pub const SPRUCE_SITKA_LOSS: &[u8] =
    include_bytes!("../fixtures/matdb/spruce-sitka-loss-qiu2026.fsmatpk");
/// Brazilian rosewood, FPL-GTR-282 (E_L and specific gravity only).
pub const ROSEWOOD_BRAZILIAN: &[u8] =
    include_bytes!("../fixtures/matdb/rosewood-brazilian-fpl-gtr282.fsmatpk");
/// Indian rosewood, FPL-GTR-282 (E_L and specific gravity only).
pub const ROSEWOOD_INDIAN: &[u8] =
    include_bytes!("../fixtures/matdb/rosewood-indian-fpl-gtr282.fsmatpk");
/// Aluminium 2024-T3, NASA TN D-6448, 300 K.
pub const ALUMINUM_2024_T3: &[u8] =
    include_bytes!("../fixtures/matdb/aluminum-2024-t3-nasa-tn-d6448.fsmatpk");
/// Music-wire density, NBS Circular 447.
pub const MUSIC_WIRE_DENSITY: &[u8] =
    include_bytes!("../fixtures/matdb/music-wire-nbs-c447.fsmatpk");
/// Music-wire spring moduli, Fuchs 1968.
pub const MUSIC_WIRE_MODULI: &[u8] =
    include_bytes!("../fixtures/matdb/music-wire-spring-moduli-fuchs-1968.fsmatpk");
/// Cane reed heel, static bending, dry and water-soaked.
pub const CANE_REED_HEEL: &[u8] =
    include_bytes!("../fixtures/matdb/cane-arundo-reedheel-mdpi-ma13204566.fsmatpk");
/// Cane reed blank density range.
pub const CANE_REED_BLANK: &[u8] =
    include_bytes!("../fixtures/matdb/cane-arundo-reedblank-mdpi-ma18122759.fsmatpk");
/// Cane whole-wall DMA damping (logarithmic decrement range).
pub const CANE_DAMPING: &[u8] =
    include_bytes!("../fixtures/matdb/cane-arundo-damping-scielo-mr20170795.fsmatpk");
/// Authored fixture: the guitar back's rosewood ratios, Poisson ratio and density.
pub const ROSEWOOD_BACK_FIXTURE: &[u8] =
    include_bytes!("../fixtures/matdb/rosewood-brazilian-back-authored-fixture.fsmatpk");
/// Authored fixture: the xylophone bar's rosewood density.
pub const ROSEWOOD_BAR_FIXTURE: &[u8] =
    include_bytes!("../fixtures/matdb/rosewood-indian-bar-authored-fixture.fsmatpk");
/// Authored fixture: the reed's cane Poisson ratio.
pub const CANE_REED_FIXTURE: &[u8] =
    include_bytes!("../fixtures/matdb/cane-reed-authored-fixture.fsmatpk");

/// FPL-GTR-282 ch. 5 tables are at 12 % moisture content; the compiler
/// normalizes '%' by 0.01.
pub const FPL_MOISTURE_CONTENT: f64 = 0.12;

/// One compiled pack as an immutable material card.
pub struct MaterialPack {
    /// The pack id the bytes declare.
    pub id: String,
    /// Hex content hash of the compiled pack, for receipts and provenance.
    pub hash: String,
    pack: NormalizedMaterialCardPack,
}

impl MaterialPack {
    /// A seed pack compiled from `data/matdb/seed-v1/<id>`.
    pub fn seed(bytes: &[u8]) -> Self {
        Self::wrap(bytes, "seed-v1 compiled pack")
    }

    /// An authored fixture compiled from `tests/fixtures/matdb/<id>/`. Its
    /// terms must carry [`FIXTURE_LABEL`], so a data pack cannot pass as one.
    pub fn fixture(bytes: &[u8]) -> Self {
        let pack = Self::wrap(bytes, AUTHORED_FIXTURE);
        let plain = NormalizedPack::from_bytes(bytes).expect("compiled fixture pack");
        assert!(
            plain.redistribution_terms().contains(FIXTURE_LABEL),
            "{}: an authored fixture's terms must read '{FIXTURE_LABEL}'",
            pack.id
        );
        pack
    }

    fn wrap(bytes: &[u8], process: &str) -> Self {
        let plain = NormalizedPack::from_bytes(bytes).expect("compiled pack");
        let (id, hash) = (plain.pack_id().to_string(), plain.content_hash().to_hex());
        let pack = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: format!("pack:{id}"),
                phase: "solid".into(),
                process: process.into(),
                revision: 0,
            },
            plain,
        )
        .unwrap_or_else(|e| panic!("{id}: {e:?}"));
        Self { id, hash, pack }
    }

    /// The card every resolution reads.
    pub fn card(&self) -> &MaterialCard {
        self.pack.card()
    }

    /// `id hash` for provenance text.
    pub fn provenance(&self) -> String {
        format!("{} {}", self.id, self.hash)
    }

    /// One scalar property [SI] at `point` through the state-point resolver.
    ///
    /// # Errors
    /// A missing, ambiguous, out-of-domain or dimensionally wrong property
    /// refuses, naming the pack and the property.
    pub fn scalar(&self, property: &str, dims: Dims, point: &QueryPoint) -> Result<f64, String> {
        let requirement =
            ScalarPropertyRequirement::try_new(property, dims, ScalarAdmissibility::Finite)
                .map_err(|e| format!("pack {} property '{property}': {e:?}", self.id))?;
        let resolved = resolve_material_state_point(
            self.card(),
            point,
            &[requirement],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|e| format!("pack {} property '{property}': {e:?}", self.id))?;
        Ok(resolved
            .property(property)
            .expect("a resolved requirement")
            .value_si())
    }

    /// The single state its one claim of `property` was measured at: every
    /// validity axis at its lower edge, in its own quantity schema. A claim
    /// constant within its validity has the same value anywhere inside it.
    ///
    /// # Errors
    /// No claim, or more than one, refuses by name.
    pub fn measured_state(&self, property: &str) -> Result<QueryPoint, String> {
        let claims = self.card().claims().claims_for(property);
        let [(_, claim)] = claims.as_slice() else {
            return Err(format!(
                "pack {} has {} claims of '{property}'; need exactly one",
                self.id,
                claims.len()
            ));
        };
        let mut point = QueryPoint::new();
        for (axis, &(lo, _)) in claim.validity.bounds() {
            point = match claim.validity.axis_quantities().get(axis) {
                Some(quantity) => point.with_quantity(axis, *quantity, lo),
                None => point.with(axis, lo),
            }
            .map_err(|e| format!("pack {} '{property}' axis {axis}: {e:?}", self.id))?;
        }
        Ok(point)
    }

    /// `scalar` at the property's own measured state.
    ///
    /// # Errors
    /// As [`Self::measured_state`] and [`Self::scalar`].
    pub fn scalar_as_measured(&self, property: &str, dims: Dims) -> Result<f64, String> {
        self.scalar(property, dims, &self.measured_state(property)?)
    }

    /// The stated absolute half-width [SI] of the one claim of `property`.
    ///
    /// # Errors
    /// A missing claim or an unstated or relative uncertainty refuses by name.
    pub fn half_width(&self, property: &str) -> Result<f64, String> {
        let claims = self.card().claims().claims_for(property);
        match claims.as_slice() {
            [(_, claim)] => match claim.uncertainty {
                UncertaintyModel::HalfWidth { half_width, .. } => Ok(half_width),
                ref other => Err(format!(
                    "pack {} '{property}' states no absolute half-width: {other:?}",
                    self.id
                )),
            },
            _ => Err(format!(
                "pack {} has {} claims of '{property}'; need exactly one",
                self.id,
                claims.len()
            )),
        }
    }

    /// The wood's longitudinal-radial plane from its reference ratios at
    /// `point`: density [kg/m^3] and the plane law (E_L, E_R, nu_LR, G_LR).
    ///
    /// # Errors
    /// Any missing ratio, modulus or density refuses by name.
    pub fn plane_lr(&self, point: &QueryPoint) -> Result<(f64, OrthotropicPlaneElastic), String> {
        match resolve_orthotropic_elastic_state_point_from_reference_ratios(
            self.card(),
            point,
            MaterialPropertySelection::SingleClaimOnly,
            1e-3, // declared small-strain domain, not a measured failure strain
            OrthotropicRequirementSet::Plane(MaterialPlane::Axes12),
            PoissonReciprocityRule::MajorRatios,
        )
        .map_err(|e| format!("pack {} L-R plane: {e:?}", self.id))?
        {
            ReferenceRatioOrthotropicStatePoint::Plane(state) => {
                Ok((state.density_kg_m3(), state.law().clone()))
            }
            ReferenceRatioOrthotropicStatePoint::Full3d(_) => Err(format!(
                "pack {}: a plane request returned a 3-D law",
                self.id
            )),
        }
    }
}

/// An untyped moisture-content query point (v1 packs).
pub fn at_moisture_content(fraction: f64) -> QueryPoint {
    QueryPoint::new()
        .with("moisture_content", fraction)
        .expect("moisture point")
}
