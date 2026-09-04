import SwiftUI

@main
struct FrankenSimApp: App {
    var body: some Scene {
        WindowGroup {
            SimulationStudioView()
                .background(CatalystWindowFreedom())
#if targetEnvironment(macCatalyst)
                .frame(minWidth: 760, minHeight: 560)
#endif
        }
#if targetEnvironment(macCatalyst)
        .defaultSize(width: 1420, height: 920)
        .windowResizability(.contentMinSize)
#endif
        .commands {
            SidebarCommands()
            CommandMenu("Simulation") {
                Button("Run Current Kernel") {
                    NotificationCenter.default.post(name: .runSimulation, object: nil)
                }
                .keyboardShortcut("r", modifiers: .command)
            }
            SimulationTextSizeCommands()
        }
    }
}

private struct SimulationTextSizeCommands: Commands {
    @AppStorage(ForgeTheme.textScaleStorageKey) private var textScale = ForgeTheme.defaultTextScale

    var body: some Commands {
        CommandMenu("Text Size") {
            Button("Larger Text") {
                textScale = ForgeTheme.adjustedTextScale(from: textScale, steps: 1)
            }
            .keyboardShortcut("+", modifiers: .command)
            .disabled(ForgeTheme.normalizedTextScale(textScale) >= ForgeTheme.maximumTextScale)

            Button("Smaller Text") {
                textScale = ForgeTheme.adjustedTextScale(from: textScale, steps: -1)
            }
            .keyboardShortcut("-", modifiers: .command)
            .disabled(ForgeTheme.normalizedTextScale(textScale) <= ForgeTheme.minimumTextScale)

            Button("Actual Size") {
                textScale = ForgeTheme.defaultTextScale
            }
            .keyboardShortcut("0", modifiers: .command)
            .disabled(ForgeTheme.normalizedTextScale(textScale) == ForgeTheme.defaultTextScale)
        }
    }
}

extension Notification.Name {
    static let runSimulation = Notification.Name("FrankenSim.runSimulation")
}
