//! Native probability samples use the ordinary import, solve and sealed-QoI
//! producers. Recovery reads the retained model, never the original paths.

use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::{CancelGate, Cx};
use fs_ledger::{EdgeRole, Ledger};
use fs_project::DecodedProject;
use fs_project::uncertainty::{BoundStudy, UncertaintyStudy};

use super::{Result, artifact, fail, quoted};
use crate::json_read::JsonValue as J;
use crate::{
    CardPackKind, CardPackSet, GeometryImportLimits, RawCardPack, RawGeometryLibrary,
    SOLVE_DRIVER_VERSION, SolveRunId, SolveRunStatus,
};

pub(super) const MANIFEST_KIND: &str = "native-uncertainty-model";
const SOURCE_KIND: &str = "native-uncertainty-source";
const SCHEMA: &str = "frankensim.cli.native-uncertainty-model.v1";
const MANIFEST_CAP: u64 = 128 * 1024;
const INPUT_CAP: u64 = 256 * 1024 * 1024;

pub(super) struct Model {
    pub(super) bound: BoundStudy,
    base: DecodedProject,
    geometry: Vec<Vec<u8>>,
    cards: CardPackSet,
    manifest: String,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Sample {
    pub(super) run: String,
    pub(super) project_hash: String,
    pub(super) qoi_receipt: ContentHash,
    pub(super) value_k: f64,
    pub(super) parameters: Vec<f64>,
}

fn invalid(message: impl Into<String>) -> super::Failure {
    fail("cli-uncertainty-model", message)
}
fn project_error(error: fs_project::ProjectError) -> super::Failure {
    invalid(format!("{}: {}", error.code, error.detail))
}
fn ledger_error(error: fs_ledger::LedgerError) -> super::Failure {
    fail("cli-uncertainty-ledger", error.to_string())
}
fn utf8(bytes: &[u8]) -> Result<&str> {
    std::str::from_utf8(bytes).map_err(|error| invalid(format!("invalid UTF-8: {error}")))
}
fn hash_field(value: &J, name: &str) -> Result<ContentHash> {
    value
        .str_field(name)
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid(format!("missing or invalid model {name}")))
}
fn array<'a>(value: &'a J, name: &str, limit: usize) -> Result<&'a [J]> {
    let rows = value
        .get(name)
        .and_then(J::as_array)
        .ok_or_else(|| invalid(format!("missing model {name}")))?;
    if rows.len() > limit {
        return Err(invalid(format!("model {name} exceeds {limit} entries")));
    }
    Ok(rows)
}
fn relative(directory: &Path, name: &str) -> Result<PathBuf> {
    if Path::new(name).is_absolute() {
        return Err(invalid(
            "study asset paths must be relative to the study file",
        ));
    }
    Ok(directory.join(name))
}
fn read(path: &Path, cap: u64) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path)
        .map_err(|error| invalid(format!("cannot open {}: {error}", path.display())))?;
    let metadata = file
        .metadata()
        .map_err(|error| invalid(error.to_string()))?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err(invalid(format!(
            "{} must be a regular file within {cap} bytes",
            path.display()
        )));
    }
    let mut bytes = Vec::new();
    file.take(cap.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|error| invalid(format!("cannot read {}: {error}", path.display())))?;
    if bytes.len() as u64 > cap {
        return Err(invalid("input exceeded its admitted byte cap"));
    }
    Ok(bytes)
}
fn input_limit(base: &DecodedProject) -> Result<u64> {
    let memory = base
        .spec
        .budgets
        .as_ref()
        .ok_or_else(|| invalid("missing project budgets"))?
        .memory_bytes;
    // Retained bytes, decoded cards and the per-sample raw-source copy coexist.
    // This is an explicit input-storage envelope, not a measured RSS promise.
    Ok((memory / 4).min(INPUT_CAP))
}
fn charge(total: &mut u64, length: usize, cap: u64) -> Result<()> {
    *total = total
        .checked_add(length as u64)
        .filter(|value| *value <= cap)
        .ok_or_else(|| {
            invalid("aggregate native-study inputs exceed the declared memory envelope")
        })?;
    Ok(())
}

impl Model {
    pub(super) fn load(path: &Path) -> Result<Self> {
        let source = read(path, fs_project::uncertainty::MAX_SOURCE_BYTES as u64)?;
        let study = UncertaintyStudy::parse(utf8(&source)?).map_err(project_error)?;
        let directory = path.parent().unwrap_or_else(|| Path::new("."));
        let project_path = relative(directory, study.project_path())?;
        let base = crate::read_project_for_solve(&project_path, crate::OutputMode::Json)
            .map_err(|output| invalid(output.stderr))?;
        // Bind before reading assets: bad targets and probability support do
        // not spend geometry/card resources or create a ledger.
        let bound = study.bind(&base.spec).map_err(project_error)?;
        let cap = input_limit(&base)?;
        let mut used = 0;
        charge(&mut used, bound.study().canonical().len(), cap)?;
        charge(&mut used, base.canonical.len(), cap)?;
        let mut geometry = Vec::new();
        for source in bound.study().geometry() {
            let bytes = read(
                &relative(directory, &source.path)?,
                (GeometryImportLimits::DEFAULT.max_source_bytes as u64).min(cap - used),
            )?;
            charge(&mut used, bytes.len(), cap)?;
            geometry.push(bytes);
        }
        let mut packs = Vec::new();
        for (kind, paths) in [
            (CardPackKind::Material, bound.study().materials()),
            (CardPackKind::Interface, bound.study().interfaces()),
        ] {
            for path in paths {
                let bytes = read(
                    &relative(directory, path)?,
                    crate::MAX_CARD_PACK_BYTES.min(cap - used),
                )?;
                charge(&mut used, bytes.len(), cap)?;
                packs.push(RawCardPack {
                    kind,
                    source: path.clone(),
                    bytes,
                    expect: None,
                });
            }
        }
        let cards = CardPackSet::admit(packs).map_err(|error| invalid(error.to_string()))?;
        Self::from_parts(bound, base, geometry, cards)
    }

    fn from_parts(
        bound: BoundStudy,
        base: DecodedProject,
        geometry: Vec<Vec<u8>>,
        cards: CardPackSet,
    ) -> Result<Self> {
        if geometry.len() != bound.study().geometry().len() {
            return Err(invalid(
                "retained geometry does not cover the study declaration",
            ));
        }
        let cap = input_limit(&base)?;
        let mut used = 0;
        charge(&mut used, bound.study().canonical().len(), cap)?;
        charge(&mut used, base.canonical.len(), cap)?;
        for bytes in &geometry {
            charge(&mut used, bytes.len(), cap)?;
        }
        for pack in cards.iter() {
            charge(&mut used, pack.bytes().len(), cap)?;
        }
        let geometry_json = bound
            .study()
            .geometry()
            .iter()
            .zip(&geometry)
            .map(|(source, bytes)| {
                format!(
                    "{{\"role\":{},\"source\":{}}}",
                    quoted(&source.role),
                    quoted(&hash_bytes(bytes).to_hex())
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let cards_json = cards
            .iter()
            .map(|pack| {
                format!(
                    "{{\"kind\":{},\"source\":{}}}",
                    quoted(pack.kind().artifact_kind()),
                    quoted(&pack.artifact().to_hex())
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let manifest = format!(
            concat!(
                "{{\"schema\":{},\"solve_driver\":{},\"constellation_lock\":{},",
                "\"study\":{},\"base\":{},\"geometry\":[{}],\"cards\":[{}],\"cards_root\":{}}}"
            ),
            quoted(SCHEMA),
            SOLVE_DRIVER_VERSION,
            quoted(&hash_bytes(include_bytes!("../../../../constellation.lock")).to_hex()),
            quoted(&hash_bytes(bound.study().canonical().as_bytes()).to_hex()),
            quoted(&hash_bytes(base.canonical.as_bytes()).to_hex()),
            geometry_json,
            cards_json,
            quoted(&cards.root().to_hex())
        );
        if manifest.len() as u64 > MANIFEST_CAP {
            return Err(invalid("native model manifest exceeds its cap"));
        }
        Ok(Self {
            bound,
            base,
            geometry,
            cards,
            manifest,
        })
    }

    pub(super) fn identity(&self) -> ContentHash {
        hash_bytes(self.manifest.as_bytes())
    }

    pub(super) fn retain(&self, ledger: &Ledger, op: i64) -> Result<ContentHash> {
        let mut linked = std::collections::BTreeSet::new();
        let mut retain = |kind: &str, bytes: &[u8]| -> Result<ContentHash> {
            let hash = ledger
                .put_artifact(kind, bytes, None)
                .map_err(ledger_error)?
                .hash;
            if linked.insert(hash) {
                ledger.link(op, &hash, EdgeRole::In).map_err(ledger_error)?;
            }
            Ok(hash)
        };
        retain(SOURCE_KIND, self.bound.study().canonical().as_bytes())?;
        retain("solve-project-source", self.base.canonical.as_bytes())?;
        for bytes in &self.geometry {
            retain("geometry-source", bytes)?;
        }
        for pack in self.cards.iter() {
            retain(pack.kind().artifact_kind(), pack.bytes())?;
        }
        retain(MANIFEST_KIND, self.manifest.as_bytes())
    }

    pub(super) fn restore(ledger: &Ledger, manifest: ContentHash) -> Result<Self> {
        let bytes = artifact(ledger, manifest, MANIFEST_KIND, MANIFEST_CAP)?;
        let value = J::parse(utf8(&bytes)?).map_err(|error| invalid(error.to_string()))?;
        let source = artifact(
            ledger,
            hash_field(&value, "study")?,
            SOURCE_KIND,
            fs_project::uncertainty::MAX_SOURCE_BYTES as u64,
        )?;
        let study = UncertaintyStudy::parse(utf8(&source)?).map_err(project_error)?;
        let base_bytes = artifact(
            ledger,
            hash_field(&value, "base")?,
            "solve-project-source",
            crate::MAX_PROJECT_BYTES,
        )?;
        let base = fs_project::parse_sexpr(utf8(&base_bytes)?).map_err(project_error)?;
        let bound = study.bind(&base.spec).map_err(project_error)?;
        let cap = input_limit(&base)?;
        let mut used = 0;
        charge(&mut used, source.len(), cap)?;
        charge(&mut used, base_bytes.len(), cap)?;
        let rows = array(&value, "geometry", 32)?;
        if rows.len() != bound.study().geometry().len() {
            return Err(invalid("retained geometry count differs"));
        }
        let mut geometry = Vec::new();
        for (row, declaration) in rows.iter().zip(bound.study().geometry()) {
            if row.str_field("role") != Some(declaration.role.as_str()) {
                return Err(invalid("retained geometry role differs"));
            }
            let bytes = artifact(
                ledger,
                hash_field(row, "source")?,
                "geometry-source",
                (GeometryImportLimits::DEFAULT.max_source_bytes as u64).min(cap - used),
            )?;
            charge(&mut used, bytes.len(), cap)?;
            geometry.push(bytes);
        }
        let mut packs = Vec::new();
        for row in array(&value, "cards", crate::MAX_CARD_PACKS)? {
            let kind = row
                .str_field("kind")
                .and_then(CardPackKind::from_artifact_kind)
                .ok_or_else(|| invalid("retained native-study card kind is invalid"))?;
            let bytes = artifact(
                ledger,
                hash_field(row, "source")?,
                kind.artifact_kind(),
                crate::MAX_CARD_PACK_BYTES.min(cap - used),
            )?;
            charge(&mut used, bytes.len(), cap)?;
            packs.push(RawCardPack {
                kind,
                source: "retained-native-study".into(),
                bytes,
                expect: None,
            });
        }
        let cards = CardPackSet::admit(packs).map_err(|error| invalid(error.to_string()))?;
        let model = Self::from_parts(bound, base, geometry, cards)?;
        if model.manifest.as_bytes() != bytes
            || model.bound.study().canonical().as_bytes() != source
            || model.base.canonical.as_bytes() != base_bytes
        {
            return Err(invalid(
                "retained model differs from the canonical inputs or current solver/constellation",
            ));
        }
        Ok(model)
    }

    fn project(&self, values: &[f64]) -> Result<DecodedProject> {
        let project = self.bound.sample_project(values).map_err(project_error)?;
        let source = fs_project::print_sexpr(&project).map_err(project_error)?;
        fs_project::parse_sexpr(&source).map_err(project_error)
    }

    pub(super) fn sample(
        &self,
        ledger: &Ledger,
        gate: &CancelGate,
        values: &[f64],
        remaining_wall_s: f64,
    ) -> Result<Option<Sample>> {
        if !remaining_wall_s.is_finite() || remaining_wall_s < 0.0 {
            return Err(invalid(
                "sample wall allowance must be finite and nonnegative",
            ));
        }
        if remaining_wall_s == 0.0 || gate.is_requested() {
            return Ok(None);
        }
        let duration = Duration::try_from_secs_f64(remaining_wall_s)
            .map_err(|_| invalid("sample wall allowance is outside the clock range"))?;
        // The watchdog only requests existing cooperative cancellation. It
        // terminates immediately after the model call; no background work escapes.
        std::thread::scope(|scope| {
            let (done, wait) = std::sync::mpsc::channel();
            scope.spawn(move || {
                if matches!(
                    wait.recv_timeout(duration),
                    Err(std::sync::mpsc::RecvTimeoutError::Timeout)
                ) {
                    gate.request();
                }
            });
            let result = self.sample_inner(ledger, gate, values);
            let _ = done.send(());
            result
        })
    }

    fn sample_inner(
        &self,
        ledger: &Ledger,
        gate: &CancelGate,
        values: &[f64],
    ) -> Result<Option<Sample>> {
        let project = self.project(values)?;
        let run = SolveRunId::derive(&project, &self.cards).to_hex();
        let existing = crate::solve::load_completed_run(ledger, &run);
        let resume = match existing {
            Ok(sealed) => return self.extract(ledger, sealed, &project, values).map(Some),
            Err(error) if error.code == "cli-solve-unknown-run" => false,
            Err(error) if error.code == "cli-report-run-incomplete" => true,
            Err(error) => return Err(invalid(format!("{}: {}", error.code, error.what))),
        };
        if gate.is_requested() {
            return Ok(None);
        }
        if !resume {
            let mut raw = RawGeometryLibrary::new();
            for (source, bytes) in self.bound.study().geometry().iter().zip(&self.geometry) {
                let declaration = project
                    .spec
                    .geometry
                    .as_ref()
                    .expect("bound geometry")
                    .iter()
                    .find(|row| row.role == source.role)
                    .expect("bound geometry role");
                raw.insert_mesh(
                    declaration,
                    source.role.clone(),
                    bytes.clone(),
                    source.unit.clone(),
                    source.max_hole_edges,
                    Vec::new(),
                );
            }
            let mut limits = GeometryImportLimits::DEFAULT;
            let memory = usize::try_from(
                project
                    .spec
                    .budgets
                    .as_ref()
                    .expect("bound budgets")
                    .memory_bytes,
            )
            .unwrap_or(usize::MAX);
            limits.max_source_bytes = limits.max_source_bytes.min(memory);
            limits.max_total_source_bytes = limits.max_total_source_bytes.min(memory);
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            let imported = pool.scope(|arena| {
                let cx = Cx::new(
                    gate,
                    arena,
                    fs_exec::StreamKey {
                        seed: project.spec.seeds.as_ref().expect("bound seed").root,
                        kernel_id: 0x66_73_63_6c_69_69_6d_70,
                        tile: 0,
                        iteration: 0,
                    },
                    fs_exec::Budget::INFINITE,
                    fs_exec::ExecMode::Deterministic,
                );
                crate::import_project_geometry(&project.spec, &raw, ledger, limits, &cx)
            });
            match imported {
                Ok(_) => {}
                Err(error) if error.code == "cli-import-cancelled" => return Ok(None),
                Err(error) => return Err(invalid(format!("{}: {}", error.code, error.what))),
            }
        }
        if gate.is_requested() {
            return Ok(None);
        }
        let started = Instant::now();
        let mut clock = || started.elapsed().as_secs_f64();
        let mut progress = Vec::new();
        let outcome = if resume {
            crate::resume_solve(ledger, gate, &mut clock, &run, &mut progress)
        } else {
            crate::run_solve(
                ledger,
                gate,
                &mut clock,
                &project,
                &self.cards,
                &mut progress,
            )
        };
        match outcome {
            Ok(outcome) if matches!(outcome.status, SolveRunStatus::Completed) => {
                let sealed = crate::solve::load_completed_run(ledger, &run)
                    .map_err(|error| invalid(format!("{}: {}", error.code, error.what)))?;
                self.extract(ledger, sealed, &project, values).map(Some)
            }
            Ok(_) => Ok(None),
            Err(error)
                if matches!(
                    error.code,
                    "cli-solve-cancelled" | "cli-solve-resume-budget"
                ) =>
            {
                Ok(None)
            }
            Err(error) => Err(invalid(format!("{}: {}", error.code, error.what))),
        }
    }

    fn extract(
        &self,
        ledger: &Ledger,
        sealed: crate::solve::CompletedRunExport,
        project: &DecodedProject,
        values: &[f64],
    ) -> Result<Sample> {
        let run = SolveRunId::derive(project, &self.cards).to_hex();
        let project_hash = project.hash().to_hex();
        if sealed.run != run || sealed.project_hash != project_hash {
            return Err(invalid(
                "sealed sample belongs to another project or card set",
            ));
        }
        self.verify_geometry(ledger, &sealed, project)?;
        let qoi_receipt = sealed
            .stages
            .iter()
            .find(|row| row.0 == "qoi")
            .and_then(|row| ContentHash::from_hex(&row.2))
            .ok_or_else(|| invalid("sealed sample lacks its QoI receipt"))?;
        let bytes = artifact(ledger, qoi_receipt, "solve-stage-receipt", 4 * 1024 * 1024)?;
        let receipt = J::parse(utf8(&bytes)?).map_err(|error| invalid(error.to_string()))?;
        let rows = array(&receipt, "qoi", 1)?;
        let row = rows
            .first()
            .ok_or_else(|| invalid("sample has no temperature observation"))?;
        let region = &project
            .spec
            .requirements
            .as_ref()
            .expect("bound requirement")[0]
            .region;
        if receipt.str_field("run") != Some(run.as_str())
            || receipt.path(&["lineage", "project"]).and_then(J::as_str)
                != Some(project_hash.as_str())
            || row.str_field("name") != Some(self.bound.study().qoi())
            || row.str_field("region") != Some(region.as_str())
            || row.str_field("unit") != Some("kelvin")
        {
            return Err(invalid(
                "sealed sample QoI name, region, unit or project differs",
            ));
        }
        let value_k = row
            .f64_field("value")
            .filter(|value| value.is_finite() && *value >= 0.0)
            .ok_or_else(|| invalid("sealed sample temperature is not finite nonnegative Kelvin"))?;
        Ok(Sample {
            run,
            project_hash,
            qoi_receipt,
            value_k,
            parameters: values.to_vec(),
        })
    }

    fn verify_geometry(
        &self,
        ledger: &Ledger,
        sealed: &crate::solve::CompletedRunExport,
        project: &DecodedProject,
    ) -> Result<()> {
        // A native run id includes the project's legacy geometry-row identity.
        // Reuse also requires the strong source hashes and actual import policy
        // retained by that run to match THIS model's immutable source assets.
        let hash = sealed
            .stages
            .iter()
            .find(|row| row.0 == "import-verify")
            .and_then(|row| ContentHash::from_hex(&row.2))
            .ok_or_else(|| invalid("sealed sample lacks geometry verification"))?;
        let bytes = artifact(ledger, hash, "solve-stage-receipt", 4 * 1024 * 1024)?;
        let receipt = J::parse(utf8(&bytes)?).map_err(|error| invalid(error.to_string()))?;
        let rows = array(&receipt, "verified", 32)?;
        if rows.len() != self.geometry.len() {
            return Err(invalid("sealed sample geometry count differs"));
        }
        let import_op = receipt
            .get("import_op")
            .and_then(J::number_raw)
            .and_then(|value| value.parse::<i64>().ok())
            .filter(|value| *value > 0)
            .ok_or_else(|| invalid("sealed sample import operation is invalid"))?;
        let operation = ledger
            .op(import_op)
            .map_err(ledger_error)?
            .ok_or_else(|| invalid("sealed sample import operation is missing"))?;
        let ir = J::parse(&operation.ir).map_err(|error| invalid(error.to_string()))?;
        let policies = array(&ir, "sources", 32)?;
        if policies.len() != self.geometry.len() {
            return Err(invalid("sealed sample import policies differ"));
        }
        for (source, bytes) in self.bound.study().geometry().iter().zip(&self.geometry) {
            let declaration = project
                .spec
                .geometry
                .as_ref()
                .expect("bound geometry")
                .iter()
                .find(|row| row.role == source.role)
                .expect("bound geometry role");
            let identity = fs_project::geometry_source_identity(declaration);
            let row = rows
                .iter()
                .find(|row| row.str_field("role") == Some(source.role.as_str()))
                .ok_or_else(|| invalid("sealed sample geometry role differs"))?;
            let policy = policies
                .iter()
                .find(|row| row.str_field("source_identity") == Some(identity.as_str()))
                .and_then(|row| row.get("policy"))
                .ok_or_else(|| invalid("sealed sample import policy is missing"))?;
            if hash_field(row, "raw_source")? != hash_bytes(bytes)
                || row.str_field("source_identity") != Some(identity.as_str())
                || policy.str_field("kind") != Some("mesh")
                || policy.str_field("length_unit") != Some(source.unit.as_str())
                || policy
                    .get("max_hole_edges")
                    .and_then(J::number_raw)
                    .and_then(|value| value.parse::<usize>().ok())
                    != Some(source.max_hole_edges)
                || policy
                    .get("named_groups")
                    .and_then(J::as_array)
                    .is_none_or(|groups| !groups.is_empty())
            {
                return Err(invalid(
                    "sealed sample used different geometry bytes or import policy",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn verify_sample(&self, ledger: &Ledger, sample: &Sample) -> Result<()> {
        let project = self.project(&sample.parameters)?;
        let sealed = crate::solve::load_completed_run(ledger, &sample.run)
            .map_err(|error| invalid(format!("{}: {}", error.code, error.what)))?;
        let actual = self.extract(ledger, sealed, &project, &sample.parameters)?;
        if actual.run != sample.run
            || actual.project_hash != sample.project_hash
            || actual.qoi_receipt != sample.qoi_receipt
            || actual.value_k.to_bits() != sample.value_k.to_bits()
        {
            return Err(invalid(
                "retained sample differs from its sealed native solve",
            ));
        }
        Ok(())
    }
}
