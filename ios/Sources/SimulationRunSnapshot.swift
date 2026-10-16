import CoreTransferable
import Foundation
import UniformTypeIdentifiers

struct SimulationRunSnapshot: Codable, Equatable, Sendable {
    static let schema = "frankensim.apple-run-snapshot.v1"
    static let boundary = "Native visualization snapshot; not a Design Ledger certificate."

    let schema: String
    let createdAtMilliseconds: Int64
    let experiment: Experiment
    let run: Run
    let payload: [Double]
    let noClaim: String
    let authorityBoundary: String

    struct Experiment: Codable, Equatable, Sendable {
        let id: UInt32
        let name: String
        let kernel: String
        let evidence: String
    }

    struct Run: Codable, Equatable, Sendable {
        let packetSchema: UInt32
        let shape: String
        let width: Int
        let height: Int
        let frames: Int
        let seed: UInt32
        let quality: Double
        let qualityTier: String
        let elapsedMilliseconds: Double
        let finiteMinimum: Double
        let finiteMaximum: Double
    }

    init(
        result: SimulationResult,
        experiment: SimulationExperiment,
        createdAt: Date = .now
    ) throws {
        guard result.experimentID == experiment.id else {
            throw SimulationRunSnapshotError.identityMismatch
        }
        let scaledDate = createdAt.timeIntervalSince1970 * 1_000
        let elapsed = Self.seconds(result.elapsed) * 1_000
        guard scaledDate.isFinite,
              let milliseconds = Int64(exactly: scaledDate.rounded()),
              elapsed.isFinite,
              elapsed >= 0 else {
            throw SimulationRunSnapshotError.invalidMetadata
        }
        schema = Self.schema
        createdAtMilliseconds = milliseconds
        self.experiment = Experiment(
            id: experiment.id,
            name: experiment.name,
            kernel: experiment.kernel,
            evidence: experiment.evidence.rawValue
        )
        run = Run(
            packetSchema: 1,
            shape: String(describing: result.shape),
            width: result.width,
            height: result.height,
            frames: result.frames,
            seed: result.seed,
            quality: result.quality,
            qualityTier: Self.qualityTier(result.quality),
            elapsedMilliseconds: elapsed,
            finiteMinimum: result.finiteMinimum,
            finiteMaximum: result.finiteMaximum
        )
        payload = result.values
        noClaim = experiment.noClaim
        authorityBoundary = Self.boundary
    }

    func encoded() throws -> Data {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        encoder.nonConformingFloatEncodingStrategy = .convertToString(
            positiveInfinity: "+infinity",
            negativeInfinity: "-infinity",
            nan: "nan"
        )
        return try encoder.encode(self)
    }

    var filename: String {
        let stem = experiment.name.lowercased().map { character in
            character.isLetter || character.isNumber ? character : "-"
        }
        .split(separator: "-")
        .filter { !$0.isEmpty }
        .joined(separator: "-")
        .prefix(48)
        return "frankensim-\(experiment.id)-\(stem)-\(String(run.seed, radix: 16)).json"
    }

    private static func seconds(_ duration: Duration) -> Double {
        let components = duration.components
        return Double(components.seconds) + Double(components.attoseconds) / 1e18
    }

    private static func qualityTier(_ quality: Double) -> String {
        if quality < 0.3 { return "quick" }
        if quality < 0.8 { return "balanced" }
        return "deep"
    }
}

struct SimulationRunSnapshotFile: Transferable, Sendable {
    let snapshot: SimulationRunSnapshot

    static var transferRepresentation: some TransferRepresentation {
        FileRepresentation(exportedContentType: .json) { file in
            let directory = FileManager.default.temporaryDirectory
                .appendingPathComponent(UUID().uuidString, isDirectory: true)
            try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
            let url = directory.appendingPathComponent(file.snapshot.filename, isDirectory: false)
            try file.snapshot.encoded().write(to: url, options: .atomic)
            return SentTransferredFile(url)
        }
    }
}

enum SimulationRunSnapshotError: LocalizedError {
    case identityMismatch
    case invalidMetadata

    var errorDescription: String? {
        switch self {
        case .identityMismatch:
            "The result does not belong to the selected experiment."
        case .invalidMetadata:
            "The run metadata cannot be represented in a portable snapshot."
        }
    }
}
