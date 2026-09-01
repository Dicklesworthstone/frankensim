import AVFoundation
import Foundation
import FrankenSimCore

enum ResultShape: UInt32, Sendable {
    case signal = 0
    case grid = 1
    case gridFrames = 2
    case xyzPath = 3
    case triangles = 4
    case campaign = 5
    /// ABI dimensions are frame count, channel count, and sample rate in hertz.
    case pcm = 6
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
        case .pcm: "\(width) mono PCM frames @ \(frames) Hz"
        }
    }
}

enum PCMPlaybackError: LocalizedError {
    case notPCM
    case unsupportedFormat(channels: Int, sampleRate: Int)
    case malformedPayload(String)
    case bufferAllocation
    case engineStart(String)

    var errorDescription: String? {
        switch self {
        case .notPCM: "The selected kernel result is not an audio PCM packet."
        case .unsupportedFormat(let channels, let sampleRate):
            "Unsupported PCM format: \(channels) channel(s) at \(sampleRate) Hz."
        case .malformedPayload(let reason): "Invalid PCM payload: \(reason)."
        case .bufferAllocation: "The audio engine could not allocate a PCM buffer."
        case .engineStart(let reason): "The audio engine could not start: \(reason)"
        }
    }
}

/// Typed Swift ownership of one Rust-produced, mono PCM block. The C ABI's
/// `width`, `height`, and `frames` dimensions carry frame count, channels, and
/// sample rate respectively for `ResultShape.pcm` only.
struct PCMBlock: Sendable {
    static let sampleRate = 48_000
    static let channels = 1

    let samples: [Float]

    init(result: SimulationResult) throws {
        guard result.shape == .pcm else { throw PCMPlaybackError.notPCM }
        guard result.width > 0, result.width == result.values.count else {
            throw PCMPlaybackError.malformedPayload("frame count does not match the payload")
        }
        guard result.height == Self.channels, result.frames == Self.sampleRate else {
            throw PCMPlaybackError.unsupportedFormat(channels: result.height, sampleRate: result.frames)
        }
        guard result.values.allSatisfy({ $0.isFinite && (-1...1).contains($0) }) else {
            throw PCMPlaybackError.malformedPayload("samples must be finite normalized PCM")
        }
        samples = result.values.map(Float.init)
    }
}

/// Real AVAudioEngine consumer for the Rust PCM packet. `scheduledBuffer`
/// retains the buffer for the node's complete scheduled lifetime; it is only
/// replaced after `stop()` interrupts the preceding block.
@MainActor
final class PCMPlayback {
    private let engine = AVAudioEngine()
    private let player = AVAudioPlayerNode()
    private let format: AVAudioFormat
    private var scheduledBuffer: AVAudioPCMBuffer?

    init() throws {
        guard let format = AVAudioFormat(
            standardFormatWithSampleRate: Double(PCMBlock.sampleRate),
            channels: AVAudioChannelCount(PCMBlock.channels)
        ) else {
            throw PCMPlaybackError.unsupportedFormat(
                channels: PCMBlock.channels,
                sampleRate: PCMBlock.sampleRate
            )
        }
        self.format = format
        engine.attach(player)
        engine.connect(player, to: engine.mainMixerNode, format: format)
    }

    func play(_ block: PCMBlock) throws {
        guard let buffer = AVAudioPCMBuffer(
            pcmFormat: format,
            frameCapacity: AVAudioFrameCount(block.samples.count)
        ), let destination = buffer.floatChannelData?.pointee else {
            throw PCMPlaybackError.bufferAllocation
        }
        for (index, sample) in block.samples.enumerated() {
            destination[index] = sample
        }
        buffer.frameLength = AVAudioFrameCount(block.samples.count)

        player.stop()
        scheduledBuffer = buffer
        do {
            if !engine.isRunning {
                try engine.start()
            }
        } catch {
            scheduledBuffer = nil
            throw PCMPlaybackError.engineStart(error.localizedDescription)
        }
        player.scheduleBuffer(buffer, at: nil, options: .interrupts)
        player.play()
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
    private var pcmPlayback: PCMPlayback?
    private var runTask: Task<Void, Never>?
    private var activeRunID: UUID?

    init() {
#if DEBUG
        let environment = ProcessInfo.processInfo.environment
        if let rawID = environment["FSIM_INITIAL_EXPERIMENT"],
           let id = UInt32(rawID),
           let experiment = SimulationCatalog.all.first(where: { $0.id == id })
        {
            selection = experiment
        }
        if let rawQuality = environment["FSIM_INITIAL_QUALITY"],
           let requested = Double(rawQuality)
        {
            quality = min(0.92, max(0.12, requested))
        }
#endif
    }

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
                if output.shape == .pcm {
                    let block = try PCMBlock(result: output)
                    let playback: PCMPlayback
                    if let existing = self.pcmPlayback {
                        playback = existing
                    } else {
                        let created = try PCMPlayback()
                        self.pcmPlayback = created
                        playback = created
                    }
                    try playback.play(block)
                }
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
