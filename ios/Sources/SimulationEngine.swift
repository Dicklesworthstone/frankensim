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
    case streamScheduling
    case starved

    var errorDescription: String? {
        switch self {
        case .notPCM: "The selected kernel result is not an audio PCM packet."
        case .unsupportedFormat(let channels, let sampleRate):
            "Unsupported PCM format: \(channels) channel(s) at \(sampleRate) Hz."
        case .malformedPayload(let reason): "Invalid PCM payload: \(reason)."
        case .bufferAllocation: "The audio engine could not allocate a PCM buffer."
        case .engineStart(let reason): "The audio engine could not start: \(reason)"
        case .streamScheduling: "The PCM stream exceeded its bounded refill schedule."
        case .starved: "The Rust PCM producer did not refill the audio queue before it drained."
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

/// A finite, one-second reed demonstration. Four queued 10 ms blocks give the
/// non-real-time Rust producer 40 ms of refill slack without unbounded audio
/// buffering or a device-dependent scheduling policy.
struct PCMStreamPlan: Sendable, Equatable {
    static let reedDemo = PCMStreamPlan(totalBlocks: 100, maximumQueuedBlocks: 4)

    let totalBlocks: Int
    let maximumQueuedBlocks: Int

    init(totalBlocks: Int, maximumQueuedBlocks: Int) {
        precondition(totalBlocks > 0 && maximumQueuedBlocks > 0)
        self.totalBlocks = totalBlocks
        self.maximumQueuedBlocks = maximumQueuedBlocks
    }
}

enum PCMStreamDrain: Equatable {
    case awaitingRefill
    case complete
    case starved
    case stopped
}

/// Device-independent bookkeeping for the bounded player queue. It is the
/// deterministic test seam for initial prefill, recurring refills, completion,
/// explicit drain, and starvation detection.
struct PCMStreamSchedule: Equatable {
    let plan: PCMStreamPlan
    private(set) var scheduledBlocks = 0
    private(set) var drainedBlocks = 0
    private(set) var queuedBlocks = 0
    private(set) var stopping = false
    private(set) var playbackStarted = false

    init(plan: PCMStreamPlan) {
        self.plan = plan
    }

    mutating func reserveNextBlock() -> Bool {
        guard !stopping,
              scheduledBlocks < plan.totalBlocks,
              queuedBlocks < plan.maximumQueuedBlocks else {
            return false
        }
        scheduledBlocks += 1
        queuedBlocks += 1
        return true
    }

    mutating func requestStopAfterDrain() {
        stopping = true
    }

    var isReadyToStartPlayback: Bool {
        !playbackStarted
            && queuedBlocks == min(plan.totalBlocks, plan.maximumQueuedBlocks)
    }

    var stopsImmediately: Bool {
        !playbackStarted
    }

    mutating func markPlaybackStarted() {
        precondition(isReadyToStartPlayback)
        playbackStarted = true
    }

    mutating func didDrainBlock() -> PCMStreamDrain {
        precondition(queuedBlocks > 0)
        queuedBlocks -= 1
        drainedBlocks += 1
        if stopping {
            return queuedBlocks == 0 ? .stopped : .awaitingRefill
        }
        if drainedBlocks == plan.totalBlocks {
            return .complete
        }
        return queuedBlocks == 0 ? .starved : .awaitingRefill
    }
}

/// The current stream identifier is captured by AVAudioPlayerNode completion
/// callbacks. A callback from a replaced stream must not mutate its successor.
struct PCMStreamGeneration: Equatable {
    private(set) var currentID: UUID?

    mutating func begin() -> UUID {
        let id = UUID()
        currentID = id
        return id
    }

    func accepts(_ id: UUID) -> Bool {
        currentID == id
    }
}

/// Thread-safe bounded handoff between the dedicated Rust producer closure and
/// the main-actor AVAudioEngine scheduler. A permit represents one queued PCM
/// buffer; stopping wakes a waiting producer without publishing more audio.
final class PCMStreamControl: @unchecked Sendable {
    private let permits: DispatchSemaphore
    private let lock = NSLock()
    private var stopped = false
    private let wakeCount: Int

    init(maximumQueuedBlocks: Int) {
        permits = DispatchSemaphore(value: maximumQueuedBlocks)
        wakeCount = maximumQueuedBlocks
    }

    func acquireProducerPermit() -> Bool {
        while true {
            if permits.wait(timeout: .now() + .milliseconds(50)) == .success {
                lock.lock()
                let isStopped = stopped
                lock.unlock()
                if isStopped {
                    permits.signal()
                    return false
                }
                return true
            }
            lock.lock()
            let isStopped = stopped
            lock.unlock()
            if isStopped { return false }
        }
    }

    func releasePlaybackPermit() {
        permits.signal()
    }

    func stop() {
        lock.lock()
        let wasStopped = stopped
        stopped = true
        lock.unlock()
        if !wasStopped {
            for _ in 0..<wakeCount {
                permits.signal()
            }
        }
    }
}

/// Real AVAudioEngine consumer for recurring Rust PCM. Each scheduled buffer
/// remains retained until the player reports data playback; producer permits
/// limit the queued lead and completion callbacks detect starvation.
@MainActor
final class PCMPlayback {
    private let engine = AVAudioEngine()
    private let player = AVAudioPlayerNode()
    private let format: AVAudioFormat
    private var scheduledBuffers: [UUID: AVAudioPCMBuffer] = [:]
    private var schedule: PCMStreamSchedule?
    private var control: PCMStreamControl?
    private var generation = PCMStreamGeneration()
    private var onTerminalError: @MainActor (PCMPlaybackError) -> Void = { _ in }
    private var onFinished: @MainActor () -> Void = {}

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

    func begin(
        plan: PCMStreamPlan,
        onTerminalError: @escaping @MainActor (PCMPlaybackError) -> Void,
        onFinished: @escaping @MainActor () -> Void
    ) -> (UUID, PCMStreamControl) {
        player.stop()
        scheduledBuffers.removeAll()
        control?.stop()
        let control = PCMStreamControl(maximumQueuedBlocks: plan.maximumQueuedBlocks)
        let id = generation.begin()
        self.control = control
        schedule = PCMStreamSchedule(plan: plan)
        self.onTerminalError = onTerminalError
        self.onFinished = onFinished
        return (id, control)
    }

    func enqueue(_ block: PCMBlock, streamID incomingStreamID: UUID) throws {
        guard generation.accepts(incomingStreamID),
              var schedule,
              schedule.reserveNextBlock() else {
            throw PCMPlaybackError.streamScheduling
        }
        self.schedule = schedule
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

        guard let control else {
            throw PCMPlaybackError.streamScheduling
        }

        let token = UUID()
        scheduledBuffers[token] = buffer
        player.scheduleBuffer(
            buffer,
            at: nil,
            options: [],
            completionCallbackType: .dataPlayedBack
        ) { [weak self] _ in
            Task { @MainActor in
                self?.didDrain(token: token, streamID: incomingStreamID, control: control)
            }
        }
        try startPlaybackWhenPrefilled()
    }

    func stopAndDrain() {
        guard var schedule else { return }
        let stopsImmediately = schedule.stopsImmediately
        schedule.requestStopAfterDrain()
        self.schedule = schedule
        control?.stop()
        if stopsImmediately || scheduledBuffers.isEmpty {
            player.stop()
            engine.pause()
            scheduledBuffers.removeAll()
            onFinished()
        }
    }

    private func startPlaybackWhenPrefilled() throws {
        guard var schedule, schedule.isReadyToStartPlayback else { return }
        do {
            if !engine.isRunning {
                try engine.start()
            }
        } catch {
            throw PCMPlaybackError.engineStart(error.localizedDescription)
        }
        schedule.markPlaybackStarted()
        self.schedule = schedule
        player.play()
    }

    private func didDrain(token: UUID, streamID incomingStreamID: UUID, control: PCMStreamControl) {
        guard generation.accepts(incomingStreamID),
              var schedule,
              scheduledBuffers.removeValue(forKey: token) != nil else {
            return
        }
        control.releasePlaybackPermit()
        let outcome = schedule.didDrainBlock()
        self.schedule = schedule
        switch outcome {
        case .awaitingRefill:
            break
        case .complete, .stopped:
            player.stop()
            engine.pause()
            onFinished()
        case .starved:
            abort(.starved)
        }
    }

    func abort(_ error: PCMPlaybackError) {
        control?.stop()
        player.stop()
        engine.pause()
        scheduledBuffers.removeAll()
        onTerminalError(error)
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
        try Self.runPacket(experiment: experiment, quality: quality, seed: seed)
    }

    /// One synchronous C-ABI round trip. The bounded PCM producer invokes this
    /// repeatedly within one dedicated dispatch closure so the Rust voice's
    /// thread-local render state remains on that producer thread.
    nonisolated static func runPacket(
        experiment: SimulationExperiment,
        quality: Double,
        seed: UInt32
    ) throws -> SimulationResult {
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
    private let pcmProducerQueue = DispatchQueue(
        label: "org.frankensim.reed-pcm-producer",
        qos: .userInitiated
    )
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
        pcmPlayback?.stopAndDrain()
        let experiment = selection
        let runQuality = quality
        let runSeed = seed
        let runID = UUID()
        activeRunID = runID
        isRunning = true
        errorMessage = nil
        if experiment.id == 43 {
            do {
                try startPCMStream(
                    experiment: experiment,
                    quality: runQuality,
                    seed: runSeed,
                    runID: runID
                )
            } catch {
                guard activeRunID == runID else { return }
                errorMessage = error.localizedDescription
                activeRunID = nil
                isRunning = false
            }
            return
        }
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
        pcmPlayback?.stopAndDrain()
        activeRunID = nil
        isRunning = false
    }

    private func startPCMStream(
        experiment: SimulationExperiment,
        quality: Double,
        seed: UInt32,
        runID: UUID
    ) throws {
        let playback: PCMPlayback
        if let existing = pcmPlayback {
            playback = existing
        } else {
            let created = try PCMPlayback()
            pcmPlayback = created
            playback = created
        }
        let plan = PCMStreamPlan.reedDemo
        let (streamID, control) = playback.begin(
            plan: plan,
            onTerminalError: { [weak self] error in
                guard let self, self.activeRunID == runID else { return }
                self.errorMessage = error.localizedDescription
                self.activeRunID = nil
                self.isRunning = false
            },
            onFinished: { [weak self] in
                guard let self, self.activeRunID == runID else { return }
                self.activeRunID = nil
                self.isRunning = false
            }
        )

        pcmProducerQueue.async {
            for blockIndex in 0..<plan.totalBlocks {
                guard control.acquireProducerPermit() else { return }
                do {
                    let result = try NativeSimulationKernel.runPacket(
                        experiment: experiment,
                        quality: quality,
                        seed: seed
                    )
                    let block = try PCMBlock(result: result)
                    Task { @MainActor [weak self] in
                        guard let self, self.activeRunID == runID else {
                            control.stop()
                            return
                        }
                        do {
                            try playback.enqueue(block, streamID: streamID)
                            if blockIndex == 0 {
                                self.result = result
                            }
                        } catch let error as PCMPlaybackError {
                            playback.abort(error)
                        } catch {
                            playback.abort(.malformedPayload(error.localizedDescription))
                        }
                    }
                } catch {
                    control.stop()
                    Task { @MainActor [weak self] in
                        guard let self, self.activeRunID == runID else { return }
                        playback.abort(.malformedPayload(error.localizedDescription))
                    }
                    return
                }
            }
        }
    }
}
