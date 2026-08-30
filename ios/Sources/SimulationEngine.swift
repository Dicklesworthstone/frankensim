import Foundation
import FrankenSimCore

enum ResultShape: UInt32, Sendable {
    case signal = 0
    case grid = 1
    case gridFrames = 2
    case xyzPath = 3
    case triangles = 4
    case campaign = 5
}

struct SimulationResult: Sendable {
    let experimentID: UInt32
    let shape: ResultShape
    let width: Int
    let height: Int
    let frames: Int
    let values: [Double]
    let finiteMinimum: Double
    let finiteMaximum: Double
    let quality: Double
    let seed: UInt32
    let elapsed: Duration

    var elapsedText: String {
        let components = elapsed.components
        let seconds = Double(components.seconds) + Double(components.attoseconds) / 1e18
        if seconds < 1 { return String(format: "%.0f ms", seconds * 1_000) }
        return String(format: "%.2f s", seconds)
    }

    var payloadSummary: String {
        switch shape {
        case .grid: "\(width) × \(height) field"
        case .gridFrames: "\(width) × \(height) · \(frames) frames"
        case .xyzPath: "\(values.count / 3) spatial samples"
        case .triangles: "\(max(0, values.count / 18)) triangles"
        case .signal: "\(values.count) values"
        case .campaign: "\(values.count) evidence values"
        }
    }
}

enum KernelError: LocalizedError {
    case incompatibleSchema(UInt32)
    case refused(Int32)
    case malformed(String)

    var errorDescription: String? {
        switch self {
        case .incompatibleSchema(let schema): "Unsupported native result schema \(schema)."
        case .refused(let code): "The native kernel refused this run (code \(code))."
        case .malformed(let reason): "The native kernel returned a malformed result: \(reason)."
        }
    }
}

actor NativeSimulationKernel {
    private static let headerCount = 6
    private static let maximumValues = 2_000_006

    func run(experiment: SimulationExperiment, quality: Double, seed: UInt32) throws -> SimulationResult {
        let clock = ContinuousClock()
        let start = clock.now
        let reported = Int(frankensim_apple_run(experiment.id, quality, seed))
        let error = frankensim_apple_last_error()
        guard error == 0 else { throw KernelError.refused(error) }
        guard reported == Int(frankensim_apple_result_len()),
              reported >= Self.headerCount,
              reported <= Self.maximumValues else {
            throw KernelError.malformed("invalid packet length")
        }

        var packet = [Double]()
        packet.reserveCapacity(reported)
        for index in 0..<reported {
            if Task.isCancelled { throw CancellationError() }
            packet.append(frankensim_apple_result_value(UInt64(index)))
        }

        guard packet[0].isFinite, packet[0].rounded() == packet[0] else {
            throw KernelError.malformed("non-integral schema")
        }
        let schema = UInt32(packet[0])
        guard schema == frankensim_apple_schema_version() else { throw KernelError.incompatibleSchema(schema) }
        guard packet[1] == Double(experiment.id),
              let shapeCode = Self.dimension(packet[2]),
              shapeCode <= Int(UInt32.max),
              let shape = ResultShape(rawValue: UInt32(shapeCode)),
              let width = Self.dimension(packet[3]),
              let height = Self.dimension(packet[4]),
              let frames = Self.dimension(packet[5]) else {
            throw KernelError.malformed("invalid metadata")
        }
        let values = Array(packet.dropFirst(Self.headerCount))
        let finite = values.filter(\.isFinite)
        guard !finite.isEmpty else { throw KernelError.malformed("no finite payload") }
        if shape == .grid || shape == .gridFrames {
            let expected = width.multipliedReportingOverflow(by: height)
            guard !expected.overflow else { throw KernelError.malformed("dimension overflow") }
            let total = expected.partialValue.multipliedReportingOverflow(by: frames)
            guard !total.overflow, values.count == total.partialValue else {
                throw KernelError.malformed("field payload does not match its dimensions")
            }
        }
        if shape == .xyzPath, !values.count.isMultiple(of: 3) {
            throw KernelError.malformed("spatial path is not made of xyz triples")
        }
        if shape == .triangles {
            guard let triangleCount = values.first.flatMap(Self.dimension),
                  values.count == 1 + triangleCount * 18 else {
                throw KernelError.malformed("triangle payload does not match its declared count")
            }
        }
        return SimulationResult(
            experimentID: experiment.id,
            shape: shape,
            width: width,
            height: height,
            frames: frames,
            values: values,
            finiteMinimum: finite.min() ?? 0,
            finiteMaximum: finite.max() ?? 1,
            quality: quality,
            seed: seed,
            elapsed: start.duration(to: clock.now)
        )
    }

    private static func dimension(_ value: Double) -> Int? {
        guard value.isFinite, value >= 0, value <= 2_000_000, value.rounded() == value else { return nil }
        return Int(value)
    }
}

@MainActor
final class SimulationStudioModel: ObservableObject {
    @Published var selection = SimulationCatalog.initial
    @Published var result: SimulationResult?
    @Published var isRunning = false
    @Published var errorMessage: String?
    @Published var quality = 0.55
    @Published var seed: UInt32 = 0x5EED_0001

    private let kernel = NativeSimulationKernel()
    private var runTask: Task<Void, Never>?
    private var activeRunID: UUID?

    func select(_ experiment: SimulationExperiment) {
        guard selection.id != experiment.id else { return }
        runTask?.cancel()
        selection = experiment
        result = nil
        errorMessage = nil
        run()
    }

    func run() {
        runTask?.cancel()
        let experiment = selection
        let runQuality = quality
        let runSeed = seed
        let runID = UUID()
        activeRunID = runID
        isRunning = true
        errorMessage = nil
        runTask = Task { [weak self] in
            guard let self else { return }
            do {
                let output = try await kernel.run(experiment: experiment, quality: runQuality, seed: runSeed)
                try Task.checkCancellation()
                guard self.activeRunID == runID, self.selection.id == experiment.id else { return }
                self.result = output
                self.activeRunID = nil
                self.isRunning = false
            } catch is CancellationError {
                if self.activeRunID == runID {
                    self.activeRunID = nil
                    self.isRunning = false
                }
            } catch {
                guard self.activeRunID == runID, self.selection.id == experiment.id else { return }
                self.errorMessage = error.localizedDescription
                self.activeRunID = nil
                self.isRunning = false
            }
        }
    }

    func randomizeSeed() {
        seed = UInt32.random(in: UInt32.min...UInt32.max)
        run()
    }

    func cancel() {
        runTask?.cancel()
        activeRunID = nil
        isRunning = false
    }
}
