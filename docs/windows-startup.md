# Windows startup and downloaded releases

Download the current [Windows release](https://github.com/zzhenglab/COFdraw/releases/latest), extract the complete ZIP into a writable folder, then open `COFdraw.exe`. Keep `COFdraw.exe.config`, `_internal`, and `data` beside the executable. Administrator access is not required. The native window uses Microsoft's WebView2 Runtime.

Version 0.6.1 fixes a startup failure in downloaded copies of v0.6.0. Windows can preserve an Internet download mark on extracted DLLs. The .NET Framework loader then rejects a bundled assembly, sometimes reporting `0x80131515` or `Failed to resolve Python.Runtime.Loader.Initialize`. Microsoft documents this behavior for locally stored assemblies downloaded from the web. [Microsoft .NET Framework guidance](https://learn.microsoft.com/en-us/dotnet/framework/configure-apps/file-schema/runtime/loadfromremotesources-element).

The release includes an adjacent `COFdraw.exe.config` that enables `loadFromRemoteSources` for this application's CLR process. The setting allows assemblies loaded by COFdraw to run with full trust; it does not change machine-wide policy or remove Windows download marks. Keep the application folder intact and use the official release. The package manifest records the configuration file's SHA-256 hash.

If the native window still fails, check that the full ZIP was extracted and that the sidecar configuration is present. Install or repair the [Microsoft WebView2 Runtime](https://developer.microsoft.com/en-us/microsoft-edge/webview2/) if it is missing. Startup details are written to `data/cofdraw.log` inside the extracted application folder, or to the selected `--data-dir`. A native startup failure returns a nonzero exit status. Running as administrator, changing global .NET configuration, or recursively unblocking a Downloads folder is unnecessary for this fix.

The Windows executable is unsigned. SmartScreen can therefore show an unrecognized-app warning on the first launch. Verify the official release and its published hash before using **More info → Run anyway** for `COFdraw.exe`. This is a decision for that application, not a request to disable Windows security. Windows may remove the EXE's download mark after this confirmation; the bundled DLLs retain theirs. If organizational policy blocks the application, follow that policy rather than disabling protection.

## Release regression check

From a Windows checkout with an interactive Explorer desktop, run:

```powershell
python scripts/smoke_cofdraw_download.py --zip dist/COFdraw-Windows.zip --output-dir test-results/download-startup-check
```

Choose a fresh output directory. The harness extracts a separate copy into a path containing spaces, verifies that the release manifest binds the executable's configuration, and adds `Zone.Identifier` with `ZoneId=3` to the EXE, configuration, and every bundled DLL. A plain PowerShell helper creates a transparent `.lnk` containing the exact EXE, arguments, and working directory. The system `explorer.exe` opens that shortcut through the ordinary desktop launch path. The harness does not use encoded PowerShell or a COM execution/de-elevation mechanism. [Microsoft Shell Links documentation](https://learn.microsoft.com/en-us/windows/win32/shell/links).

The native smoke window stays hidden, while Windows' own first-run confirmation remains visible. Review and handle that exact application's prompt within the test timeout. The harness checks the actual COFdraw process token for medium integrity (`8192`); an elevated process fails the check even when the application opens.

Success requires the native document and save bridge, process exit, local server shutdown, unchanged file hashes, and retained Internet marks on every DLL and the configuration. The JSON report records the EXE mark separately before and after launch: `download_marks_retained` is false if that mark disappears after Windows confirmation, while `dll_config_download_marks_retained` must remain true. The harness never removes or restores marks. It never counts a headless server or browser fallback as native success. Existing archives and application folders are left unchanged; generated evidence stays in the selected output directory. If the interactive Explorer desktop or medium-token launch is unavailable, the check fails rather than silently switching to an elevated launch.

Windows download marks are alternate data streams, separate from DLL contents and their SHA-256 hashes. Microsoft's [Unblock-File documentation](https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.utility/unblock-file) describes the stream and Internet zone value. This regression test deliberately leaves those marks present.
