Download **COFdraw-macOS-arm64.zip** for Apple silicon (M-series) or **COFdraw-macOS-x86_64.zip** for Intel. In **Apple menu → About This Mac**, look for **Chip** or **Processor**.

Unzip the download, move **COFdraw.app** to **Applications**, and open it. The chemistry libraries and catalog are included. The first opening may require per-app approval in **System Settings → Privacy & Security → Open Anyway**. The beta has an ad-hoc signature and has not been notarized by Apple. [Apple’s instructions](https://support.apple.com/en-us/102445)

Both downloads are exercised on native macOS 15 runners. Release checks open the WebKit window, check its save bridge, render a COF-300 drawing, export its crystal as CIF, and run an argon adsorption smoke test. The Intel build uses RDKit 2025.9.3, the last Python 3.11 Intel Mac wheel. Other macOS versions have not been qualified.
