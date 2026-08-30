import SwiftUI

struct SimulationStudioView: View {
    @StateObject private var model = SimulationStudioModel()
    @State private var search = ""
    @State private var showsAtlas = false
    @State private var preferredColumn = NavigationSplitViewColumn.detail
    @Environment(\.scenePhase) private var scenePhase

    var body: some View {
        NavigationSplitView(preferredCompactColumn: $preferredColumn) {
            catalogSidebar
                .navigationSplitViewColumnWidth(min: 250, ideal: 300, max: 390)
        } detail: {
            GeometryReader { proxy in
                ZStack {
                    ForgeBackground()
                    if proxy.size.width >= 1_040 {
                        wideStudio(size: proxy.size)
                    } else {
                        compactStudio(size: proxy.size)
                    }
                }
            }
            .navigationBarTitleDisplayMode(.inline)
            .toolbar { studioToolbar }
        }
        .navigationSplitViewStyle(.balanced)
        .preferredColorScheme(.dark)
        .sheet(isPresented: $showsAtlas) { EpistemicAtlasView() }
        .onAppear {
            if model.result == nil && !model.isRunning { model.run() }
        }
        .onChange(of: scenePhase) { _, phase in
            if phase != .active, model.isRunning { model.cancel() }
        }
        .onReceive(NotificationCenter.default.publisher(for: .runSimulation)) { _ in
            if !model.isRunning { model.run() }
        }
    }

    private var catalogSidebar: some View {
        ZStack {
            ForgeBackground()
            List(selection: Binding(
                get: { Optional(model.selection.id) },
                set: { id in
                    guard let id, let experiment = SimulationCatalog.all.first(where: { $0.id == id }) else { return }
                    model.select(experiment)
                    preferredColumn = .detail
                }
            )) {
                Section {
                    appIdentity
                        .listRowBackground(Color.clear)
                        .listRowInsets(EdgeInsets(top: 8, leading: 4, bottom: 12, trailing: 4))
                }
                ForEach(SimulationCatalog.grouped(query: search), id: \.0) { tier, entries in
                    Section {
                        ForEach(entries) { experiment in
                            CatalogRow(experiment: experiment, selected: experiment.id == model.selection.id)
                                .tag(experiment.id)
                                .listRowBackground(Color.clear)
                        }
                    } header: {
                        HStack {
                            Text(tier.eyebrow)
                            Spacer()
                            Text("\(entries.count)")
                        }
                        .font(.system(size: ForgeTheme.size(10), weight: .bold, design: .monospaced))
                        .foregroundStyle(ForgeTheme.secondary)
                    }
                }
            }
            .scrollContentBackground(.hidden)
            .listStyle(.sidebar)
            .searchable(text: $search, placement: .sidebar, prompt: "Kernel, method, or study")
        }
        .safeAreaInset(edge: .bottom) {
            Button { showsAtlas = true } label: {
                Label("Open the theory atlas", systemImage: "books.vertical")
                    .frame(maxWidth: .infinity)
            }
            .buttonStyle(SecondaryForgeButtonStyle(tint: ForgeTheme.emerald))
            .padding(12)
            .background(.ultraThinMaterial)
        }
    }

    private var appIdentity: some View {
        HStack(spacing: 12) {
            Image("MonsterIcon")
                .resizable()
                .scaledToFill()
                .frame(width: 58, height: 58)
                .clipShape(RoundedRectangle(cornerRadius: 14, style: .continuous))
                .overlay(RoundedRectangle(cornerRadius: 14).stroke(ForgeTheme.emerald.opacity(0.45)))
            VStack(alignment: .leading, spacing: 3) {
                FrankenSimWordmark()
                Text("SIMULATION_STUDIO")
                    .font(.system(size: ForgeTheme.size(9.5), weight: .bold, design: .monospaced))
                    .kerning(1.4)
                    .foregroundStyle(ForgeTheme.secondary)
            }
        }
    }

    @ToolbarContentBuilder private var studioToolbar: some ToolbarContent {
        ToolbarItem(placement: .principal) {
            FrankenSimWordmark()
                .frame(maxWidth: 170)
        }
        ToolbarItem(placement: .primaryAction) {
            Button { showsAtlas = true } label: { Image(systemName: "info.circle") }
                .accessibilityLabel("Theory atlas")
        }
    }

    private func compactStudio(size: CGSize) -> some View {
        ScrollView {
            VStack(spacing: 14) {
                experimentHeader
                NativeSimulationCanvas(result: model.result, accent: accent, isRunning: model.isRunning)
                    .frame(height: min(440, max(286, size.height * 0.43)))
                if size.width >= 700 {
                    HStack(alignment: .top, spacing: 14) {
                        VStack(spacing: 14) {
                            controls
                            if let error = model.errorMessage { errorCard(error) }
                            resultSummary
                        }
                        evidenceCard
                    }
                } else {
                    controls
                    if let error = model.errorMessage { errorCard(error) }
                    resultSummary
                    evidenceCard
                }
            }
            .frame(maxWidth: 820)
            .padding(.horizontal, 14)
            .padding(.vertical, 12)
            .frame(maxWidth: .infinity)
        }
        .scrollIndicators(.hidden)
    }

    private func wideStudio(size: CGSize) -> some View {
        HStack(alignment: .top, spacing: 16) {
            VStack(spacing: 14) {
                experimentHeader
                NativeSimulationCanvas(result: model.result, accent: accent, isRunning: model.isRunning)
                controls
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
            ScrollView {
                VStack(spacing: 14) {
                    if let error = model.errorMessage { errorCard(error) }
                    resultSummary
                    evidenceCard
                    kernelDetails
                }
            }
            .scrollIndicators(.hidden)
            .frame(width: min(380, size.width * 0.31))
        }
        .padding(18)
    }

    private var accent: Color { ForgeTheme.accent(model.selection.accent) }

    private var experimentHeader: some View {
        HStack(alignment: .top, spacing: 13) {
            Image(systemName: model.selection.symbol)
                .font(.system(size: ForgeTheme.size(24), weight: .semibold))
                .foregroundStyle(accent)
                .frame(width: 42, height: 42)
                .background(accent.opacity(0.10), in: RoundedRectangle(cornerRadius: 12))
            VStack(alignment: .leading, spacing: 4) {
                HStack(spacing: 8) {
                    Text(model.selection.tier.eyebrow)
                        .font(.system(size: ForgeTheme.size(9.5), weight: .bold, design: .monospaced))
                        .kerning(1.4)
                        .foregroundStyle(accent)
                    Text(model.selection.kernel)
                        .font(.system(size: ForgeTheme.size(9.5), design: .monospaced))
                        .foregroundStyle(ForgeTheme.secondary)
                        .lineLimit(1)
                }
                Text(model.selection.name)
                    .font(.system(size: ForgeTheme.size(25), weight: .bold, design: .rounded))
                    .foregroundStyle(ForgeTheme.text)
                    .lineLimit(2)
                    .minimumScaleFactor(0.78)
                Text(model.selection.subtitle)
                    .font(.system(size: ForgeTheme.size(14), weight: .medium, design: .rounded))
                    .foregroundStyle(ForgeTheme.secondary)
            }
            Spacer(minLength: 4)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }

    private var controls: some View {
        ForgePanel(accent: accent) {
            ViewThatFits(in: .horizontal) {
                HStack(spacing: 12) { controlContents }
                VStack(alignment: .leading, spacing: 12) { controlContents }
            }
        }
    }

    @ViewBuilder private var controlContents: some View {
        Picker("Compute budget", selection: $model.quality) {
            Text("Quick").tag(0.12)
            Text("Balanced").tag(0.55)
            Text("Deep").tag(0.92)
        }
        .pickerStyle(.segmented)
        .frame(maxWidth: 360)
        .disabled(model.isRunning)

        HStack(spacing: 9) {
            Button {
                model.isRunning ? model.cancel() : model.run()
            } label: {
                Label(model.isRunning ? "Stop" : "Run kernel", systemImage: model.isRunning ? "stop.fill" : "bolt.fill")
            }
            .buttonStyle(PrimaryForgeButtonStyle(tint: model.isRunning ? ForgeTheme.coral : accent))

            Button { model.randomizeSeed() } label: {
                Label("New seed", systemImage: "dice")
            }
            .buttonStyle(SecondaryForgeButtonStyle(tint: ForgeTheme.cyan))
            .disabled(model.isRunning)
        }
    }

    @ViewBuilder private var resultSummary: some View {
        if let result = model.result {
            ForgePanel(accent: accent) {
                VStack(alignment: .leading, spacing: 12) {
                    panelLabel("RUN RECEIPT")
                    HStack(spacing: 16) {
                        Metric(value: result.elapsedText, label: "wall time")
                        Divider().overlay(Color.white.opacity(0.08))
                        Metric(value: result.payloadSummary, label: "native payload")
                    }
                    Text("seed 0x\(String(result.seed, radix: 16, uppercase: true)) · schema 1 · entirely on this device")
                        .font(.system(size: ForgeTheme.size(10.5), design: .monospaced))
                        .foregroundStyle(ForgeTheme.secondary)
                }
            }
        }
    }

    private var evidenceCard: some View {
        ForgePanel(accent: evidenceColor) {
            VStack(alignment: .leading, spacing: 10) {
                HStack {
                    panelLabel("EVIDENCE BOUNDARY")
                    Spacer()
                    Label(model.selection.evidence.rawValue, systemImage: evidenceSymbol)
                        .font(.system(size: ForgeTheme.size(10.5), weight: .bold, design: .rounded))
                        .foregroundStyle(evidenceColor)
                }
                Text(model.selection.explanation)
                    .font(.system(size: ForgeTheme.size(14), design: .rounded))
                    .foregroundStyle(ForgeTheme.text)
                    .fixedSize(horizontal: false, vertical: true)
                Divider().overlay(Color.white.opacity(0.08))
                Label(model.selection.noClaim, systemImage: "hand.raised")
                    .font(.system(size: ForgeTheme.size(12.5), design: .rounded))
                    .foregroundStyle(ForgeTheme.secondary)
                    .fixedSize(horizontal: false, vertical: true)
            }
        }
    }

    private var kernelDetails: some View {
        let receiptSeed = model.result?.seed ?? model.seed
        let receiptQuality = model.result?.quality ?? model.quality
        return ForgePanel(accent: ForgeTheme.violet) {
            VStack(alignment: .leading, spacing: 10) {
                panelLabel("FIVE EXPLICITS")
                detailRow("units", "fixture-defined")
                detailRow("seed", "0x\(String(receiptSeed, radix: 16, uppercase: true))")
                detailRow("budget", receiptQuality < 0.3 ? "quick" : receiptQuality < 0.8 ? "balanced" : "deep")
                detailRow("version", "native packet v1")
                detailRow("capability", model.selection.kernel)
            }
        }
    }

    private var evidenceColor: Color {
        switch model.selection.evidence {
        case .verified: ForgeTheme.emerald
        case .estimated: ForgeTheme.amber
        case .mixed: ForgeTheme.violet
        case .numerical: ForgeTheme.cyan
        }
    }

    private var evidenceSymbol: String {
        switch model.selection.evidence {
        case .verified: "checkmark.seal.fill"
        case .estimated: "waveform.badge.magnifyingglass"
        case .mixed: "circle.lefthalf.filled"
        case .numerical: "function"
        }
    }

    private func panelLabel(_ text: String) -> some View {
        Text(text)
            .font(.system(size: ForgeTheme.size(10.5), weight: .bold, design: .monospaced))
            .kerning(1.45)
            .foregroundStyle(ForgeTheme.secondary)
    }

    private func detailRow(_ name: String, _ value: String) -> some View {
        HStack {
            Text(name.capitalized).foregroundStyle(ForgeTheme.secondary)
            Spacer()
            Text(value).foregroundStyle(ForgeTheme.text)
        }
        .font(.system(size: ForgeTheme.size(12), design: .rounded))
    }

    private func errorCard(_ message: String) -> some View {
        ForgePanel(accent: ForgeTheme.coral) {
            Label(message, systemImage: "exclamationmark.triangle.fill")
                .font(.system(size: ForgeTheme.size(13), design: .rounded))
                .foregroundStyle(ForgeTheme.text)
                .frame(maxWidth: .infinity, alignment: .leading)
        }
    }
}

private struct CatalogRow: View {
    let experiment: SimulationExperiment
    let selected: Bool

    var body: some View {
        HStack(spacing: 11) {
            Image(systemName: experiment.symbol)
                .font(.system(size: ForgeTheme.size(15), weight: .semibold))
                .foregroundStyle(ForgeTheme.accent(experiment.accent))
                .frame(width: 30, height: 30)
                .background(ForgeTheme.accent(experiment.accent).opacity(0.09), in: RoundedRectangle(cornerRadius: 9))
            VStack(alignment: .leading, spacing: 2) {
                Text(experiment.name)
                    .font(.system(size: ForgeTheme.size(13.5), weight: .semibold, design: .rounded))
                    .foregroundStyle(ForgeTheme.text)
                    .lineLimit(1)
                Text(experiment.subtitle)
                    .font(.system(size: ForgeTheme.size(10.5), design: .rounded))
                    .foregroundStyle(ForgeTheme.secondary)
                    .lineLimit(1)
            }
            Spacer(minLength: 2)
            if selected { Circle().fill(ForgeTheme.emerald).frame(width: 6, height: 6) }
        }
        .padding(.vertical, 4)
        .contentShape(Rectangle())
    }
}

private struct Metric: View {
    let value: String
    let label: String

    var body: some View {
        VStack(alignment: .leading, spacing: 3) {
            Text(value)
                .font(.system(size: ForgeTheme.size(14), weight: .bold, design: .rounded))
                .foregroundStyle(ForgeTheme.text)
                .lineLimit(2)
                .minimumScaleFactor(0.76)
            Text(label.uppercased())
                .font(.system(size: ForgeTheme.size(9), weight: .bold, design: .monospaced))
                .foregroundStyle(ForgeTheme.secondary)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
}

struct EpistemicAtlasView: View {
    @Environment(\.dismiss) private var dismiss

    var body: some View {
        NavigationStack {
            ZStack {
                ForgeBackground()
                ScrollView {
                    LazyVStack(spacing: 14) {
                        atlasCard("Seven layers", "SUBSTRATE → BEDROCK → MORPH → FLUX → ASCENT → LUMEN → HELM", "From hardware and deterministic numerics through geometry, physics, optimization, rendering, and orchestration.", .cyan)
                        atlasCard("The Three Colors", "VERIFIED · VALIDATED · ESTIMATED", "Evidence composes conservatively. A weaker required ingredient can never be laundered into a stronger claim.", .emerald)
                        atlasCard("The Five Explicits", "UNITS · SEEDS · BUDGETS · VERSIONS · CAPABILITIES", "Every operation names the context required to reproduce and interpret it.", .violet)
                        atlasCard("The Gauntlet", "G0 LAWS → G5 DETERMINISM", "Properties, manufactured solutions, benchmarks, metamorphic checks, cancellation storms, and replay audits.", .amber)
                        atlasCard("Design Ledger", "CONTENT-ADDRESSED LINEAGE", "A result is useful only when its inputs, operation, budget, machine, and evidence can be reconstructed.", .coral)
                    }
                    .padding(16)
                    .frame(maxWidth: 820)
                    .frame(maxWidth: .infinity)
                }
            }
            .navigationTitle("Theory Atlas")
            .toolbar {
                ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } }
            }
        }
        .preferredColorScheme(.dark)
    }

    private func atlasCard(_ title: String, _ formula: String, _ body: String, _ accent: AccentFamily) -> some View {
        ForgePanel(accent: ForgeTheme.accent(accent)) {
            VStack(alignment: .leading, spacing: 9) {
                Text(title)
                    .font(.system(size: ForgeTheme.size(21), weight: .bold, design: .rounded))
                    .foregroundStyle(ForgeTheme.text)
                Text(formula)
                    .font(.system(size: ForgeTheme.size(10.5), weight: .bold, design: .monospaced))
                    .kerning(1.1)
                    .foregroundStyle(ForgeTheme.accent(accent))
                Text(body)
                    .font(.system(size: ForgeTheme.size(14), design: .rounded))
                    .foregroundStyle(ForgeTheme.secondary)
                    .fixedSize(horizontal: false, vertical: true)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
        }
    }
}
